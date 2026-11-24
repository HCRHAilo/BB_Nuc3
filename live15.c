/*
 * live15.c - live thermal view with the 2-point NUC from nuc2 applied.
 *
 * New in live15:
 *   -o can be given several times (up to 8 offset models, e.g. the zone
 *               models ch1z_z1..z4.bin from omodel_zones.py). Models are
 *               sorted by their T range; the model in use moves up when T
 *               passes (top - 1 C) and down when T passes (bottom + 1 C) of
 *               the current model (zones overlap, so this is a hysteresis,
 *               Lynred NTG002 "tables overlap, switch with hysteresis").
 *   -tm N        SFFC / column tables belong to the N-th -o model (default 1);
 *               they are applied only while that model is in use.
 *   T70 fix: the USB3-A3 board sensor (TMP103, UG004 Table 6) returns 8-bit
 *               two's complement; 0xF7 is -9 C (live15 read it as 247 C).
 * New in live15:
 *   -o2 hot.bin  second offset model for a hotter camera (e.g. the chamber
 *               ramp omodel_ramp1.bin, 40..62 C). The model is chosen by the
 *               temperature in use, with hysteresis: up to model 2 above
 *               (model-1 max - 1 C), back below (model-1 max - 3 C).
 *               SFFC and column tables were measured with model 1, so they
 *               are applied only while model 1 is in use.
 *               Status and log show which model is active.
 * New in live15:
 *   -k col.bin  fixed column pattern table (col_table.py): narrow column
 *               offsets vs temperature, measured once with the lens cap during
 *               a warm-up; in mode 5 the pattern for the current T is
 *               interpolated and subtracted. K: column table on/off.
 * New in live15: AGC "ours" scales a curve wider than 255 instead of clipping
 *               hot objects to flat white.
 * New in live11:
 *   G (or H)  AGC mode: linear -> HEQ -> plateau -> ours (cycle)
 *               linear  : 0.5 % tails, min span 64 DL (live11 default)
 *               HEQ     : out = 255 * CDF (Gonzalez & Woods)
 *               plateau : live11 "HE" (plateau HEQ mixed 50/50 with linear)
 *               ours    : 0.5 % tails, plateau 1 % per 1-DL bin, 20 % linear,
 *                         slope limit gmax grey/DL, damping 0.85
 *                         (knob ideas taken from the FLIR Boson+ datasheet,
 *                          only as a reference; values from agc_study.py)
 *   + / -     gmax (ours) 1 .. 6 grey per DL, step 0.5
 *   Window title shows only the frame rate; the full status line is printed
 *   in the terminal after every key press.
 *   Linux: lost frames (driver sequence gaps) and broken USB buffers are
 *   printed in the terminal with the clock time (diagnostics).
 *
 *
 * New in live11 (rest = live11):
 *   -s sffc.bin   SFFC table (sffc_table.py): smooth camera shading vs
 *                 temperature, measured once with the lens cap during a
 *                 warm-up. In mode 5 the shading for the current T is
 *                 interpolated and subtracted. No cap is needed in use.
 *   F             SFFC on/off
 *   Status shows "T OUTSIDE MODEL" when the T used is outside the omodel range.
 *   live_log_<time>.csv: every 10 s  time, T70, auto correction, T used, SFFC,
 *                 display range (one new file per run, never overwritten).
 *
 * Read-only: the camera is only streamed (via FFmpeg, same as cap.bat).
 * No sensor or proxy register is written.
 *
 * Windows build (static, no extra DLLs):
 *   x86_64-w64-mingw32-gcc -O2 -Wall -static -o live11.exe live11.c -lgdi32 -lm
 *
 * Usage (cmd, in the folder with ffmpeg on PATH):
 *   live11.exe nuc2_c1.bin                      live camera
 *   live11.exe nuc2_c1.bin -c "camera name"     other DirectShow camera name
 *   live11.exe nuc2_c1.bin -i scene1.raw        replay a recorded .raw file
 *
 * Keys in the window:
 *   1  RAW                 no correction
 *   2  1-POINT             offset only, from the reference (needs R first)
 *   3  2-POINT (calib)     corr = raw*G + O   (G, O from nuc2 calc)
 *   4  2-POINT (ref)       corr = raw*G + O', O' = mean(G*C) - G*C (needs R)
 *   5  SHUTTERLESS         corr = raw*G - p(T): offset from the temperature model
 *                          (nuc3 omodel), T = board temperature T70 read every 10 s
 *                          ( -o omodel.bin, -p /dev/ttyACM0 ), [ ] trim T by 0.5 C
 *   A  auto T on/off       estimate the model temperature from the image itself
 *   R  take reference C    point the camera at a uniform surface first
 *                          (blackbody, lens cap, plain wall); averages 32 frames
 *   C  stripe filter on/off   vertical (learned) + horizontal (per frame) lines
 *   G/H  AGC mode cycle: linear, HEQ, plateau, ours;  + / -  gmax (ours)
 *   T  temporal noise filter on/off (less grain, motion aware)
 *   S  save snapshot       live_<n>.pgm (what is on screen)
 *   SPACE freeze / unfreeze
 *   ESC quit
 *
 * Display: automatic contrast (0.5 % .. 99.5 % of the corrected values,
 * smoothed over frames) so the same AGC is used for every mode.
 * Bad pixels from the NUC file are replaced by the median of good neighbours.
 *
 * Ubuntu build (same keys, SDL2 window, camera read directly with V4L2):
 *   sudo apt install libsdl2-dev
 *   gcc -O2 -Wall -o live11 live11.c -lm $(sdl2-config --cflags --libs)
 *   ./live11 nuc2_c1.bin                 (default /dev/video0)
 *   ./live11 nuc2_c1.bin -d /dev/video2  other device
 *   ./live11 nuc2_c1.bin -i scene1.raw   replay a recorded .raw file
 * The nuc2 .bin file is the same on Windows and Ubuntu (copy it as is).
 *
 * Headless test (no window): gcc -DHEADLESS -O2 -o live11_test live11.c -lm
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#endif

#define W 640
#define H 480
#define NPIX (W * H)
#define FRAME_BYTES (NPIX * 2)
#define REF_FRAMES 32

/* Same headers as nuc2.c */
struct settings {
    int32_t tint, gfid, gsk, gain_reg;
    float   t70_c;
};
struct nuchdr {
    char     magic[8];
    uint32_t w, h;
    float    m1, m2, t1_c, t2_c;
    int32_t  nbad;
    struct settings s;
};

/* ---------------------------------------------------------------------- */
static float *G, *O, *Oref, *Cref, *corr;
static double *acc;
static uint8_t *bad;
static int *badlist, nbadlist;
static uint8_t view8[NPIX];             /* frame shown on screen            */
static double agc_lo = -1, agc_hi = -1; /* smoothed display range           */

static volatile int mode = 3;           /* 1..4, see header                 */
static volatile int ref_req = 0;        /* R pressed                        */
static volatile int ref_ok = 0;         /* reference available              */
static volatile int ref_count = 0;      /* frames collected so far          */
static volatile int frozen = 0;
static volatile int colfilt = 1;   /* on from start */        /* C: column stripe filter on/off   */
static volatile int agcm = 3;           /* G/H: 0 linear 1 HEQ 2 plateau 3 ours */
static volatile float gmax = 3.0f;      /* +/-: max grey per DL in "ours"     */
static volatile int show_status = 1;    /* print the status line once         */
static const char *agc_names[] = { "linear", "HEQ", "plateau", "ours" };
static volatile int tfilt = 1;     /* on from start */          /* T: temporal noise filter         */
static float *tprev;                    /* filtered previous frame          */
static int tprev_ok = 0;
static float colcorr[W];                /* learned per-column correction    */
static char status[400] = "starting...";
static char ref_msg[100] = "no ref";
#define REF_MAX_SD 60.0                 /* max spatial sd of a valid reference */

/* ---- offset vs temperature model (from "nuc3 omodel") ----
 * pattern p_i(T) = a_i + b_i*(T-T0) + c_i*(T-T0)^2, corrected = G*raw - p_i(T) */
struct omhdr {
    char     magic[8];                  /* "TVEOMD1" */
    uint32_t w, h, nmaps;
    float    t0, tmin, tmax;
};
static float *OA, *OB, *OC, *Omod;
static int om_ok = 0;
static float om_t0, om_tmin, om_tmax;
/* up to two offset models; OA/OB/OC and om_t0/tmin/tmax point to the active one */
#define OM_MAX 8
static struct { float *A, *B, *C, t0, tmin, tmax; int ord; } oms[OM_MAX];
static int om_n = 0, om_cur = 0, om_init = 0, tab_ord = 0;
static void use_model(int k)
{
    OA = oms[k].A; OB = oms[k].B; OC = oms[k].C;
    om_t0 = oms[k].t0; om_tmin = oms[k].tmin; om_tmax = oms[k].tmax;
    om_cur = k;
}
static volatile float t_board = NAN;    /* latest T70 from the serial port  */
static volatile float t_fixed = NAN;    /* -T: fixed temperature (no serial) */
static volatile float t_trim = 0;       /* [ ] keys: manual trim, deg C     */
static float t_used = NAN;              /* temperature Omod was built for   */

/* ---- SFFC table (from sffc_table.py): smooth shading S_k at temperature T_k ---- */
static int sf_n = 0;
static float *sf_T, *sf_S, *Sfc;        /* T[n], S[n][NPIX], current shading */
static float sf_used = NAN;             /* temperature Sfc was built for     */
static volatile int sffc_on = 1;        /* F: SFFC on/off                    */

/* ---- fixed column table (from col_table.py): pattern c_k[x] at temperature T_k ---- */
static int ct_n = 0;
static float *ct_T, *ct_C, Cfc[W];      /* T[n], C[n][W], current pattern    */
static float ct_used = NAN;
static volatile int col_on = 1;         /* K: column table on/off            */

static int load_coltab(const char *fn)
{
    char magic[8]; uint32_t wn[2];
    FILE *f = fopen(fn, "rb");
    if (!f || fread(magic, 8, 1, f) != 1 || memcmp(magic, "TVECOL1", 7) ||
        fread(wn, 4, 2, f) != 2 || wn[0] != W || wn[1] < 1 || wn[1] > 64) {
        fprintf(stderr, "%s: not a column table file\n", fn);
        if (f) fclose(f);
        return -1;
    }
    ct_n = (int)wn[1];
    ct_T = malloc(ct_n * sizeof *ct_T);
    ct_C = malloc((size_t)ct_n * W * sizeof *ct_C);
    if (!ct_T || !ct_C || fread(ct_T, 4, ct_n, f) != (size_t)ct_n ||
        fread(ct_C, 4, (size_t)ct_n * W, f) != (size_t)ct_n * W) {
        fprintf(stderr, "%s: file too short\n", fn);
        fclose(f);
        return -1;
    }
    fclose(f);
    printf("column table %s: %d entries, T %.1f..%.1f C\n", fn, ct_n, ct_T[0], ct_T[ct_n - 1]);
    return 0;
}

/* column pattern for temperature t: linear between the nearest entries, held at the ends */
static void build_col(float t)
{
    int k0 = 0, k1 = 0;
    float w = 0;
    if (t <= ct_T[0]) k0 = k1 = 0;
    else if (t >= ct_T[ct_n - 1]) k0 = k1 = ct_n - 1;
    else {
        while (k1 < ct_n - 1 && ct_T[k1] < t) k1++;
        k0 = k1 - 1;
        w = (t - ct_T[k0]) / (ct_T[k1] - ct_T[k0]);
    }
    for (int x = 0; x < W; x++) Cfc[x] = (1 - w) * ct_C[k0 * W + x] + w * ct_C[k1 * W + x];
    ct_used = t;
}

static int load_sffc(const char *fn)
{
    char magic[8]; uint32_t wh[3];
    FILE *f = fopen(fn, "rb");
    if (!f || fread(magic, 8, 1, f) != 1 || memcmp(magic, "TVESFF1", 7) ||
        fread(wh, 4, 3, f) != 3 || wh[0] != W || wh[1] != H || wh[2] < 1 || wh[2] > 64) {
        fprintf(stderr, "%s: not an SFFC table file\n", fn);
        if (f) fclose(f);
        return -1;
    }
    sf_n = (int)wh[2];
    sf_T = malloc(sf_n * sizeof *sf_T);
    sf_S = malloc((size_t)sf_n * NPIX * sizeof *sf_S);
    Sfc = calloc(NPIX, sizeof *Sfc);
    if (!sf_T || !sf_S || !Sfc || fread(sf_T, 4, sf_n, f) != (size_t)sf_n ||
        fread(sf_S, 4, (size_t)sf_n * NPIX, f) != (size_t)sf_n * NPIX) {
        fprintf(stderr, "%s: file too short\n", fn);
        fclose(f);
        return -1;
    }
    fclose(f);
    printf("SFFC table %s: %d entries, T %.1f..%.1f C\n", fn, sf_n, sf_T[0], sf_T[sf_n - 1]);
    return 0;
}

/* shading for temperature t: linear between the two nearest entries
   (entries are sorted by T), held at the first / last entry outside */
static void build_sffc(float t)
{
    int k0 = 0, k1 = 0;
    float w = 0;
    if (t <= sf_T[0]) k0 = k1 = 0;
    else if (t >= sf_T[sf_n - 1]) k0 = k1 = sf_n - 1;
    else {
        while (k1 < sf_n - 1 && sf_T[k1] < t) k1++;
        k0 = k1 - 1;
        w = (t - sf_T[k0]) / (sf_T[k1] - sf_T[k0]);
    }
    const float *a = sf_S + (size_t)k0 * NPIX, *b = sf_S + (size_t)k1 * NPIX;
    for (int i = 0; i < NPIX; i++) Sfc[i] = (1 - w) * a[i] + w * b[i];
    sf_used = t;
}


static int load_omodel(const char *fn)
{
    struct omhdr oh;
    FILE *f = fopen(fn, "rb");
    if (!f) { perror(fn); return -1; }
    if (om_n >= OM_MAX) { fprintf(stderr, "at most %d offset models\n", OM_MAX); fclose(f); return -1; }
    OA = malloc(NPIX * 4); OB = malloc(NPIX * 4); OC = malloc(NPIX * 4);
    if (!Omod) Omod = malloc(NPIX * 4);
    if (fread(&oh, sizeof oh, 1, f) != 1 || memcmp(oh.magic, "TVEOMD1", 7) || oh.w != W || oh.h != H ||
        fread(OA, 4, NPIX, f) != NPIX || fread(OB, 4, NPIX, f) != NPIX || fread(OC, 4, NPIX, f) != NPIX) {
        fprintf(stderr, "%s: not an omodel file\n", fn);
        fclose(f);
        return -1;
    }
    fclose(f);
    oms[om_n].A = OA; oms[om_n].B = OB; oms[om_n].C = OC;
    oms[om_n].t0 = oh.t0; oms[om_n].tmin = oh.tmin; oms[om_n].tmax = oh.tmax; oms[om_n].ord = om_n;
    om_n++;
    use_model(0);
    om_ok = 1;
    printf("offset model %s: %u maps, T %.0f..%.0f C\n", fn, oh.nmaps, oh.tmin, oh.tmax);
    return 0;
}

/* Offset for temperature t; outside the recorded range the value is held at
   the nearest end (no extrapolation of the quadratic) */
/* ---- auto T: estimate the temperature the image itself fits best ----
 * If the model temperature is off by d, the corrected image carries
 * d * D_i, with D_i = dp_i/dT = b_i + 2 c_i (T - T0) (known from omodel).
 * Only ONE number d is estimated per frame (not 300k offsets), so a scene
 * cannot be burned in as a ghost. Fine detail (high-pass) is used: the
 * pixel-to-pixel part of D_i looks like FPN and does not follow the scene;
 * strong scene edges are skipped. d is integrated slowly (feedback), so the
 * estimate follows slow drift (minutes) and ignores single frames. */
static volatile int autot = 1;          /* A: auto T on/off                 */
static volatile float t_auto = 0;       /* learned correction to T70, deg C */
#define AUTO_EDGE_DL 12.0f              /* skip high-pass values above this */
#define AUTO_GAIN    0.02f              /* integration step per frame       */
#define AUTO_LIMIT   10.0f              /* max |correction| deg C           */
static void auto_t_update(const float *v, float t)
{
    float x = t - om_t0;
    double srd = 0, sdd = 0;
    for (int y = 2; y < H - 2; y += 3)
        for (int xx = 2; xx < W - 2; xx += 2) {
            int i = y * W + xx;
            if (bad[i] || bad[i - 2] || bad[i + 2] || bad[i - 2 * W] || bad[i + 2 * W]) continue;
            /* high-pass: pixel minus the mean of 4 neighbours at distance 2 */
            float hr = v[i] - 0.25f * (v[i - 2] + v[i + 2] + v[i - 2 * W] + v[i + 2 * W]);
            if (hr > AUTO_EDGE_DL || hr < -AUTO_EDGE_DL) continue;
            float d0 = OB[i] + 2 * OC[i] * x;
            float hd = d0 - 0.25f * (OB[i - 2] + 2 * OC[i - 2] * x + OB[i + 2] + 2 * OC[i + 2] * x +
                                     OB[i - 2 * W] + 2 * OC[i - 2 * W] * x + OB[i + 2 * W] + 2 * OC[i + 2 * W] * x);
            srd += hr * hd; sdd += hd * hd;
        }
    if (sdd <= 0) return;
    float d = (float)(srd / sdd);                 /* deg C this frame says */
    if (d > 2) d = 2;
    if (d < -2) d = -2;
    t_auto += AUTO_GAIN * d;
    if (t_auto > AUTO_LIMIT) t_auto = AUTO_LIMIT;
    if (t_auto < -AUTO_LIMIT) t_auto = -AUTO_LIMIT;
}

/* every 600 frames (10 s): one line in live_log_<start time>.csv */
static void log_line(float t70, float tused)
{
    static FILE *lf;
    static long nfr;
    if (++nfr % 600) return;
    if (!lf) {
        char fn[64];
        snprintf(fn, sizeof fn, "live_log_%ld.csv", (long)time(NULL));
        lf = fopen(fn, "w");
        if (!lf) return;
        fprintf(lf, "time_s,t70_c,auto_corr_c,t_used_c,model,sffc,agc_lo_dl,agc_hi_dl\n");
        printf("log: %s\n", fn);
    }
    fprintf(lf, "%ld,%.1f,%.2f,%.2f,%d,%s,%.1f,%.1f\n", (long)time(NULL), t70, t_auto, tused, om_cur + 1,
            (sf_n && sffc_on && oms[om_cur].ord == tab_ord) ? "on" : "off", agc_lo, agc_hi);
    fflush(lf);
}

static void build_omod(float t)
{
    /* allow a few degrees beyond the recorded range (T70 stops at ~48 C
       while the sensor keeps drifting); further out the value is held */
    if (t < om_tmin - 5) t = om_tmin - 5;
    if (t > om_tmax + 5) t = om_tmax + 5;
    float x = t - om_t0;
    for (int i = 0; i < NPIX; i++) Omod[i] = -(OA[i] + OB[i] * x + OC[i] * x * x);
    t_used = t;
}

static int load_nuc(const char *fn)
{
    struct nuchdr nh;
    FILE *f = fopen(fn, "rb");
    if (!f) { perror(fn); return -1; }
    if (fread(&nh, sizeof nh, 1, f) != 1 || memcmp(nh.magic, "TVENUC2", 7) ||
        nh.w != W || nh.h != H ||
        fread(G, 4, NPIX, f) != NPIX || fread(O, 4, NPIX, f) != NPIX ||
        fread(bad, 1, NPIX, f) != NPIX) {
        fprintf(stderr, "%s: not a nuc2 file\n", fn);
        fclose(f);
        return -1;
    }
    fclose(f);
    nbadlist = 0;
    for (int i = 0; i < NPIX; i++) if (bad[i]) badlist[nbadlist++] = i;
    return 0;
}

static int cmp_f(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

/* Replace bad pixels by the median of their good 8 neighbours */
static void fix_bad(float *v)
{
    for (int k = 0; k < nbadlist; k++) {
        int i = badlist[k], x = i % W, y = i / W, n = 0;
        float nb[8];
        for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++) {
                int xx = x + dx, yy = y + dy;
                if ((dx || dy) && xx >= 0 && xx < W && yy >= 0 && yy < H && !bad[yy * W + xx])
                    nb[n++] = v[yy * W + xx];
            }
        if (n) { qsort(nb, n, sizeof nb[0], cmp_f); v[i] = nb[n / 2]; }
    }
}

#define AGC_MIN_DL 64.0               /* min grey-scale span, DL (~2.5 C) */

/* 14-bit histogram AGC: map 0.5 %..99.5 % to 0..255, range smoothed */
static void agc(const float *v, uint8_t *out)
{
    static uint32_t hist[16384];
    memset(hist, 0, sizeof hist);
    long n = 0;
    for (int i = 0; i < NPIX; i++) {
        if (bad[i] && mode != 1) continue;
        int b = (int)v[i];
        hist[b < 0 ? 0 : b > 16383 ? 16383 : b]++;
        n++;
    }
    long lo_n = n / 200, hi_n = n - n / 200, c = 0;
    int lo = 0, hi = 16383;
    for (int b = 0; b < 16384; b++) {
        c += hist[b];
        if (c <= lo_n) lo = b;
        if (c < hi_n) hi = b + 1;
    }
    if (hi <= lo) hi = lo + 1;
    if (agc_lo < 0) { agc_lo = lo; agc_hi = hi; }
    agc_lo += 0.1 * (lo - agc_lo);
    agc_hi += 0.1 * (hi - agc_hi);
    /* Minimum display span. On a plain wall the whole image covers only
       ~16 DL (0.6 C); stretching that to full black..white made the 2 DL
       noise and sub-DL residual lines look huge. Never show less than
       AGC_MIN_DL (~2.5 C at 26 DL/C) across the full grey scale; the span
       is kept centred on the scene. */
    if (agc_hi - agc_lo < AGC_MIN_DL) {
        double c = 0.5 * (agc_hi + agc_lo);
        agc_lo = c - 0.5 * AGC_MIN_DL;
        agc_hi = c + 0.5 * AGC_MIN_DL;
    }

    if (agcm == 1 || agcm == 3) {
        /* 1 HEQ : lut = 255 * CDF over all good pixels
           3 ours: tails cut 0.5 %, plateau P = 1 % of the pixels per 1-DL bin,
                   20 % linear mix, slope <= gmax grey/DL, centred on the median
                   when narrower than 255, damped 0.85 frame to frame */
        static float lutd[16384];
        static int lutd_ok = 0;
        static float L[16384];
        long cum = 0, med = 0, half = n / 2;
        int medb = 0;
        for (int k = 0; k < 16384; k++) { med += hist[k]; if (med < half) medb = k + 1; }
        if (agcm == 1) {
            for (int k = 0; k < 16384; k++) { cum += hist[k]; L[k] = 255.0f * cum / (n ? n : 1); }
            for (int i = 0; i < NPIX; i++) {
                int k = (int)v[i];
                k = k < 0 ? 0 : k > 16383 ? 16383 : k;
                out[i] = (uint8_t)L[k];
            }
            lutd_ok = 0;
            return;
        }
        long P = n / 100 > 1 ? n / 100 : 1, tot = 0;
        for (int k = lo; k <= hi; k++) tot += (long)hist[k] < P ? (long)hist[k] : P;
        float prev = 0, acc = 0;
        for (int k = 0; k < 16384; k++) {
            float heq, lin, t;
            if (k < lo) { heq = 0; lin = 0; }
            else if (k > hi) { heq = 255; lin = 255; }
            else {
                cum += (long)hist[k] < P ? (long)hist[k] : P;
                heq = 255.0f * cum / (tot ? tot : 1);
                lin = 255.0f * (k - lo) / (hi - lo);
            }
            t = 0.8f * heq + 0.2f * lin;              /* linear percent 20 % */
            float d = (k == 0) ? 0 : t - prev;         /* slope limit (max gain) */
            prev = t;
            if (d > gmax) d = gmax;
            acc += d;
            L[k] = acc;
        }
        /* Fit the slope-limited curve into 0..255:
           - wider than 255 (hot object far above the scene): scale it down so
             the hottest object is not clipped to flat white (live11 clipped it)
           - narrower: centre the median pixel at mid grey, but shift only as far
             as the whole scene stays inside 0..255 */
        float span = L[hi] - L[lo], scl = 1.0f, off;
        if (span > 255.0f) { scl = 255.0f / span; off = -L[lo] * scl; }
        else {
            off = 127.5f - L[medb];
            if (L[lo] + off < 0) off = -L[lo];
            if (L[hi] + off > 255) off = 255.0f - L[hi];
        }
        for (int k = 0; k < 16384; k++) {
            float t = L[k] * scl + off;
            t = t < 0 ? 0 : t > 255 ? 255 : t;
            lutd[k] = lutd_ok ? 0.85f * lutd[k] + 0.15f * t : t;   /* damping */
        }
        lutd_ok = 1;
        for (int i = 0; i < NPIX; i++) {
            int k = (int)v[i];
            k = k < 0 ? 0 : k > 16383 ? 16383 : k;
            out[i] = (uint8_t)lutd[k];
        }
        return;
    }

    if (agcm == 2) {
        /* Plateau histogram equalisation: grey levels are shared by how many
           pixels have each value, but no value may take more than "plateau"
           pixels (stops a big uniform background from eating all the grey
           levels). The mapping is smoothed over frames to avoid flicker. */
        static float lut[16384];
        static int lut_ok = 0;
        int a = (int)agc_lo, b = (int)agc_hi;
        if (a < 0) a = 0;
        if (b > 16383) b = 16383;
        if (b <= a) b = a + 1;
        long plateau = n / (b - a + 1) + 1, tot = 0;   /* clip at the average bin */
        for (int k = a; k <= b; k++) tot += (long)hist[k] < plateau ? (long)hist[k] : plateau;
        long cum = 0;
        for (int k = 0; k < 16384; k++) {
            float t;
            if (k < a) t = 0;
            else if (k > b) t = 255;
            else {
                cum += (long)hist[k] < plateau ? (long)hist[k] : plateau;
                t = 255.0f * cum / (tot ? tot : 1);
            }
            lut[k] = lut_ok ? lut[k] + 0.2f * (t - lut[k]) : t;
        }
        lut_ok = 1;
        /* Mild HE: half equalised, half linear. Pure HE pushed the face to
           white and stretched the flat background so much that 1-2 DL
           residual stripes became visible. */
        double sl = 255.0 / (agc_hi - agc_lo > 4 ? agc_hi - agc_lo : 4);
        for (int i = 0; i < NPIX; i++) {
            int k = (int)v[i];
            k = k < 0 ? 0 : k > 16383 ? 16383 : k;
            double lin = (v[i] - agc_lo) * sl;
            lin = lin < 0 ? 0 : lin > 255 ? 255 : lin;
            out[i] = (uint8_t)(0.5 * lut[k] + 0.5 * lin);
        }
        return;
    }

    double s = 255.0 / (agc_hi - agc_lo > 4 ? agc_hi - agc_lo : 4);
    for (int i = 0; i < NPIX; i++) {
        double x = (v[i] - agc_lo) * s;
        out[i] = x < 0 ? 0 : x > 255 ? 255 : (uint8_t)x;
    }
}

/* Temporal noise filter (recursive, motion aware).
 * Measured temporal noise is ~2 DL per pixel per frame (nuc2 noise).
 * Still pixels: out = prev + K*(in - prev), K = 0.25 -> noise x0.38.
 * Motion test: the change (in - prev) is averaged over the 3x3 neighbourhood.
 * Averaging 9 pixels cuts the noise 3x (~0.7 DL), so even slow movement of
 * smooth, low-contrast objects (a face) is seen as motion and passed through
 * without smearing, while pure noise stays below the threshold. */
#define TF_K      0.25f
#define MOTION_DL 2.0f            /* about 3 sigma of the 3x3 mean noise */
static float *tdiff, *trow;
static void temporal_filter(float *v)
{
    if (!tprev_ok) { memcpy(tprev, v, NPIX * sizeof *v); tprev_ok = 1; return; }
    for (int i = 0; i < NPIX; i++) tdiff[i] = v[i] - tprev[i];
    /* 3x3 box sum: horizontal pass into trow, vertical pass below */
    for (int y = 0; y < H; y++) {
        const float *d = tdiff + y * W;
        float *r = trow + y * W;
        for (int x = 0; x < W; x++) {
            int xl = x > 0 ? x - 1 : x, xr = x < W - 1 ? x + 1 : x;
            r[x] = d[xl] + d[x] + d[xr];
        }
    }
    for (int y = 0; y < H; y++) {
        int yu = y > 0 ? y - 1 : y, yd = y < H - 1 ? y + 1 : y;
        for (int x = 0; x < W; x++) {
            int i = y * W + x;
            float m = (trow[yu * W + x] + trow[i] + trow[yd * W + x]) / 9.0f;
            if (m > MOTION_DL || m < -MOTION_DL) tprev[i] = v[i];     /* motion: take new */
            else tprev[i] += TF_K * tdiff[i];                          /* still: smooth    */
            v[i] = tprev[i];
        }
    }
}

#define EDGE_DL 20.0f   /* bigger steps are real scene edges, not stripes */

/* Row stripe filter (horizontal lines).
 * nuc2 noise showed rows that move together every frame (row flicker,
 * ratio ~3.6) plus a small fixed row pattern (nuc2 stripes, ~0.5 DL).
 * Because rows change every frame, this one works per frame (no learning):
 * each pixel is compared with the average of the pixels above and below;
 * real horizontal edges (|d| > EDGE_DL) are skipped; the mean of the rest is
 * that row's stripe and is subtracted. 640 pixels per row make the estimate
 * precise (~0.1 DL) even from a single frame. */
#ifndef ROW_ITER
#define ROW_ITER 4
#endif
static void row_filter(float *v)
{
    static float rc[H];
    /* One pass removes only part of a stripe (its neighbours' share stays),
       so repeat a few times within the frame; it converges quickly. */
    for (int it = 0; it < ROW_ITER; it++) {
    for (int y = 0; y < H; y++) {
        int yu = y > 0 ? y - 1 : y + 1, yd = y < H - 1 ? y + 1 : y - 1;
        double sum = 0;
        int n = 0;
        for (int x = 0; x < W; x++) {
            int i = y * W + x;
            if (bad[i] || bad[yu * W + x] || bad[yd * W + x]) continue;
            float d = v[i] - 0.5f * (v[yu * W + x] + v[yd * W + x]);
            if (d > -EDGE_DL && d < EDGE_DL) { sum += d; n++; }
        }
        rc[y] = n > W / 4 ? 0.5f * (float)(sum / n) : 0.0f;   /* 0.5: stable step */
    }
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) v[y * W + x] -= rc[y];
    }
}

/* Column stripe filter (residual column offset, Lynred UG042 5.2: secondary
 * to the offset correction).
 * For every column, look at each pixel minus the average of its left and
 * right neighbours. Large differences are real scene edges and are skipped
 * (|d| > EDGE_DL). The mean of the remaining small differences is the
 * column's stripe offset. It is learned slowly over frames (feedback), so a
 * moving scene averages out and only the fixed stripes are removed. */
static void column_filter(float *v)
{
    for (int x = 0; x < W; x++) {
        float c = colcorr[x];
        for (int y = 0; y < H; y++) v[y * W + x] -= c;
    }
    for (int x = 0; x < W; x++) {
        int xl = x > 0 ? x - 1 : x + 1, xr = x < W - 1 ? x + 1 : x - 1;
        double sum = 0;
        int n = 0;
        for (int y = 0; y < H; y++) {
            int i = y * W;
            if (bad[i + x] || bad[i + xl] || bad[i + xr]) continue;
            float d = v[i + x] - 0.5f * (v[i + xl] + v[i + xr]);
            if (d > -EDGE_DL && d < EDGE_DL) { sum += d; n++; }
        }
        if (n > H / 4) colcorr[x] += 0.05f * (float)(sum / n);
    }
    /* Keep the correction NARROW (lines up to ~16 px).
       The neighbour test above only sees column-to-column steps; the slow
       (wide) part of colcorr is not pinned by it and drifted with noise and
       smooth scene curvature, slowly building wide vertical bands over
       minutes. Each frame a small part (2 %) of the smoothed colcorr is
       removed: wide drift cannot build up, narrow stripes are still learned
       (the triangular smoother passes almost nothing of a narrow line).
       Wide shading is the job of the offset model / auto T, not of this. */
    static float s1[W], s2[W];
    for (int pass = 0; pass < 2; pass++) {             /* box twice = triangle */
        const float *in = pass ? s1 : colcorr;
        float *out = pass ? s2 : s1;
        for (int x = 0; x < W; x++) {
            float s = 0; int c = 0;
            for (int k = x - 8; k <= x + 8; k++) if (k >= 0 && k < W) { s += in[k]; c++; }
            out[x] = s / c;
        }
    }
    for (int x = 0; x < W; x++) {
        colcorr[x] -= 0.02f * s2[x];
        if (colcorr[x] > 10.0f) colcorr[x] = 10.0f;      /* stripes are a few DL */
        if (colcorr[x] < -10.0f) colcorr[x] = -10.0f;
    }
}

/* Decode one YUY2 frame, collect reference if asked, correct, AGC */
static void process(const uint8_t *raw)
{
    for (int i = 0; i < NPIX; i++)
        corr[i] = (float)(((raw[2 * i] << 8) | raw[2 * i + 1]) >> 2);

    /* Reference: average REF_FRAMES raw frames of a uniform surface */
    if (ref_req) {
        if (ref_count == 0) memset(acc, 0, NPIX * sizeof *acc);
        for (int i = 0; i < NPIX; i++) acc[i] += corr[i];
        if (++ref_count == REF_FRAMES) {
            double gcm = 0; long n = 0;
            for (int i = 0; i < NPIX; i++) {
                Cref[i] = (float)(acc[i] / REF_FRAMES);
                if (!bad[i]) { gcm += G[i] * Cref[i]; n++; }
            }
            gcm /= n;
            /* Uniformity check: after the 2-point correction a uniform surface
               leaves only ~10-20 DL (offset drift). A person or object in view
               gives hundreds of DL and would be burned into the image as a ghost. */
            double km = 0;
            for (int i = 0; i < NPIX; i++) if (!bad[i]) km += G[i] * Cref[i] + O[i];
            km /= n;
            double ksd = 0;                          /* two-pass sd (no cancellation) */
            for (int i = 0; i < NPIX; i++)
                if (!bad[i]) { double k = G[i] * Cref[i] + O[i] - km; ksd += k * k; }
            ksd = sqrt(ksd / n);
            if (ksd > REF_MAX_SD) {
                snprintf(ref_msg, sizeof ref_msg, "REF REJECTED: not uniform (%.0f DL) - point at wall/lens cap, press R", ksd);
                ref_req = 0; ref_count = 0;
                goto ref_done;
            }
            snprintf(ref_msg, sizeof ref_msg, "ref OK (%.0f DL)", ksd);
            /* O' = mean(G*C) - G*C (Lynred NUC slide 18) */
            for (int i = 0; i < NPIX; i++) Oref[i] = (float)(gcm - G[i] * Cref[i]);
            ref_ok = 1; ref_req = 0; ref_count = 0;
            memset(colcorr, 0, sizeof colcorr);   /* new offset: relearn stripes */
        }
    }
ref_done:

    int m = mode;
    if ((m == 2 || m == 4) && !ref_ok) m = 3;   /* no reference yet */
    if (m == 5 && (!om_ok || (isnan(t_board) && isnan(t_fixed)))) m = 3;   /* no model / no T yet */
    if (m == 2) {
        double cm = 0; long n = 0;
        for (int i = 0; i < NPIX; i++) if (!bad[i]) { cm += Cref[i]; n++; }
        cm /= n;
        for (int i = 0; i < NPIX; i++) corr[i] -= (float)(Cref[i] - cm);
    } else if (m == 3) {
        for (int i = 0; i < NPIX; i++) corr[i] = corr[i] * G[i] + O[i];
    } else if (m == 4) {
        for (int i = 0; i < NPIX; i++) corr[i] = corr[i] * G[i] + Oref[i];
    } else if (m == 5) {
        float t = (isnan(t_fixed) ? t_board : t_fixed) + t_trim + (autot ? t_auto : 0);
        if (!isnan(t) && om_n > 1) {
            int want = om_cur;
            int first = !om_init;
            if (!om_init) {
                /* sort the models by their T range (insertion sort, n <= 8) */
                for (int i = 1; i < om_n; i++)
                    for (int j = i; j > 0 && oms[j].tmin < oms[j - 1].tmin; j--) {
                        __typeof__(oms[0]) tmp = oms[j]; oms[j] = oms[j - 1]; oms[j - 1] = tmp;
                    }
                /* first T: the model whose centre is nearest */
                float bd = 1e9f;
                for (int k = 0; k < om_n; k++) {
                    float d = fabsf(t - 0.5f * (oms[k].tmin + oms[k].tmax));
                    if (d < bd) { bd = d; want = k; }
                }
                om_init = 1;
            } else {
                /* move up near the top, down near the bottom (zones overlap -> hysteresis) */
                if (t > oms[om_cur].tmax - 1 && om_cur + 1 < om_n) want = om_cur + 1;
                if (t < oms[om_cur].tmin + 1 && om_cur > 0) want = om_cur - 1;
            }
            if (want != om_cur || first) {
                use_model(want);
                t_used = NAN;                        /* rebuild Omod for the new model */
                printf("T %.1f C: offset model %d/%d in use (T %.0f..%.0f C)%s\n", t, want + 1, om_n, om_tmin, om_tmax,
                       oms[want].ord == tab_ord ? "" : " - SFFC / column tables paused");
            }
        }
        if (!isnan(t) && (isnan(t_used) || fabsf(t - t_used) > 0.05f)) build_omod(t);
        if (sf_n && sffc_on && oms[om_cur].ord == tab_ord) {
            /* SFFC: subtract the smooth camera shading for this temperature */
            if (isnan(sf_used) || fabsf(t_used - sf_used) > 0.05f) build_sffc(t_used);
            for (int i = 0; i < NPIX; i++) corr[i] = corr[i] * G[i] + Omod[i] - Sfc[i];
        } else {
            for (int i = 0; i < NPIX; i++) corr[i] = corr[i] * G[i] + Omod[i];
        }
        if (ct_n && col_on && oms[om_cur].ord == tab_ord) {
            /* fixed column pattern for this temperature */
            if (isnan(ct_used) || fabsf(t_used - ct_used) > 0.05f) build_col(t_used);
            for (int i = 0; i < NPIX; i++) corr[i] -= Cfc[i % W];
        }
        if (autot) { fix_bad(corr); auto_t_update(corr, t_used); }
        log_line(isnan(t_fixed) ? t_board : t_fixed, t_used);
    }
    if (m != 1) fix_bad(corr);
    static int last_m = 0;
    if (m != last_m) { memset(colcorr, 0, sizeof colcorr); last_m = m; }
    if (colfilt && m != 1) { column_filter(corr); row_filter(corr); }
    if (tfilt) temporal_filter(corr); else tprev_ok = 0;
    agc(corr, view8);

    static const char *names[] = { "", "1 RAW", "2 1-POINT (ref)", "3 2-POINT (calib offset)",
                                   "4 2-POINT (ref offset)", "5 SHUTTERLESS (T model)" };
    char tbuf[160] = "";
    if (om_ok) snprintf(tbuf, sizeof tbuf, " | T70 %.0f%+.1f auto %s%+.1f C | T %.1f MODEL %d%s | SFFC %s | COL %s",
                        isnan(t_fixed) ? t_board : t_fixed, t_trim, autot ? "ON" : "off", t_auto, t_used, om_cur + 1,
                        (m == 5 && !isnan(t_used) && (t_used < om_tmin || t_used > om_tmax)) ? " T OUTSIDE MODEL" : "",
                        !sf_n ? "none" : (oms[om_cur].ord != tab_ord ? "paused" : sffc_on ? "ON" : "off"),
                        !ct_n ? "none" : (oms[om_cur].ord != tab_ord ? "paused" : col_on ? "ON" : "off"));
    snprintf(status, sizeof status, "%s%s%s | AGC %s %.0f..%.0f DL | stripes %s | tf %s%s%s%s",
             names[m], (m != mode) ? (mode == 5 ? "  [waiting for model/T]" : "  [press R first]") : "", tbuf,
             agc_names[agcm], agc_lo, agc_hi, colfilt ? "ON" : "off", tfilt ? "ON" : "off",
             ref_req ? " | taking reference..." : " | ", ref_req ? "" : ref_msg,
             frozen ? " | FROZEN" : "");
    if (agcm == 3) { size_t L0 = strlen(status); snprintf(status + L0, sizeof status - L0, " | gmax %.1f", gmax); }
    if (show_status) { printf("%s\n", status); show_status = 0; }
}

#ifndef HEADLESS
static void save_snapshot(void)
{
    static int n = 0;
    char fn[64];
    snprintf(fn, sizeof fn, "live_%d.pgm", ++n);
    FILE *f = fopen(fn, "wb");
    if (!f) return;
    fprintf(f, "P5\n%d %d\n255\n", W, H);
    fwrite(view8, 1, NPIX, f);
    fclose(f);
}
#endif

static int alloc_all(void)
{
    G = malloc(NPIX * 4); O = malloc(NPIX * 4); Oref = malloc(NPIX * 4);
    Cref = malloc(NPIX * 4); corr = malloc(NPIX * 4); tprev = malloc(NPIX * 4); tdiff = malloc(NPIX * 4); trow = malloc(NPIX * 4);
    acc = malloc(NPIX * sizeof *acc); bad = malloc(NPIX); badlist = malloc(NPIX * sizeof *badlist);
    return (tprev && tdiff && trow && G && O && Oref && Cref && corr && acc && bad && badlist) ? 0 : -1;
}

/* ====================================================================== */
#ifdef _WIN32

static HWND hwnd;
static CRITICAL_SECTION lock;
static uint8_t shown[NPIX];
static FILE *src;
static int from_file;
static volatile int quit;

static DWORD WINAPI reader(LPVOID arg)
{
    (void)arg;
    uint8_t *raw = malloc(FRAME_BYTES);
    DWORD t0 = GetTickCount();
    int frames = 0;
    double fps = 0;
    char st[200];
    while (!quit) {
        if (fread(raw, 1, FRAME_BYTES, src) != FRAME_BYTES) {
            if (from_file) { fseek(src, 0, SEEK_SET); continue; }   /* loop replay */
            snprintf(status, sizeof status, "stream ended - check camera / ffmpeg");
            printf("%s\n", status);
            PostMessage(hwnd, WM_APP, 0, 0);
            break;
        }
        if (from_file) Sleep(16);
        static int first = 1;
        if (first) { printf("first frame received - video OK\n"); first = 0; }
        frames++;
        DWORD t = GetTickCount();
        if (t - t0 >= 1000) { fps = frames * 1000.0 / (t - t0); frames = 0; t0 = t; }
        if (frozen && !ref_req) continue;
        EnterCriticalSection(&lock);
        process(raw);
        memcpy(shown, view8, NPIX);
        snprintf(st, sizeof st, "live15 | %.1f fps", fps);
        LeaveCriticalSection(&lock);
        SetWindowTextA(hwnd, st);
        InvalidateRect(hwnd, NULL, FALSE);
    }
    free(raw);
    return 0;
}

static LRESULT CALLBACK wndproc(HWND h, UINT msg, WPARAM wp, LPARAM lp)
{
    switch (msg) {
    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC dc = BeginPaint(h, &ps);
        RECT rc; GetClientRect(h, &rc);
        static struct { BITMAPINFOHEADER h; RGBQUAD pal[256]; } bi;
        if (!bi.h.biSize) {
            bi.h.biSize = sizeof bi.h; bi.h.biWidth = W; bi.h.biHeight = -H;
            bi.h.biPlanes = 1; bi.h.biBitCount = 8; bi.h.biCompression = BI_RGB;
            for (int i = 0; i < 256; i++)
                bi.pal[i].rgbRed = bi.pal[i].rgbGreen = bi.pal[i].rgbBlue = (BYTE)i;
        }
        EnterCriticalSection(&lock);
        SetStretchBltMode(dc, COLORONCOLOR);
        StretchDIBits(dc, 0, 0, rc.right, rc.bottom, 0, 0, W, H, shown,
                      (BITMAPINFO *)&bi, DIB_RGB_COLORS, SRCCOPY);
        LeaveCriticalSection(&lock);
        EndPaint(h, &ps);
        return 0;
    }
    case WM_APP:
        SetWindowTextA(h, status);
        return 0;
    case WM_KEYDOWN:
        show_status = 1;
        if (wp >= '1' && wp <= '5') mode = (int)(wp - '0');
        else if (wp == 'A') autot = !autot;
        else if (wp == 'F') sffc_on = !sffc_on;
        else if (wp == 'K') col_on = !col_on;
        else if (wp == VK_OEM_4) t_trim -= 0.5f;
        else if (wp == VK_OEM_6) t_trim += 0.5f;
        else if (wp == 'R') { ref_count = 0; ref_req = 1; }
        else if (wp == 'C') colfilt = !colfilt;
        else if (wp == 'H' || wp == 'G') agcm = (agcm + 1) % 4;
        else if (wp == VK_OEM_PLUS || wp == VK_ADD) { if (gmax < 6) gmax += 0.5f; }
        else if (wp == VK_OEM_MINUS || wp == VK_SUBTRACT) { if (gmax > 1) gmax -= 0.5f; }
        else if (wp == 'T') tfilt = !tfilt;
        else if (wp == 'S') { EnterCriticalSection(&lock); save_snapshot(); LeaveCriticalSection(&lock); }
        else if (wp == VK_SPACE) frozen = !frozen;
        else if (wp == VK_ESCAPE) DestroyWindow(h);
        return 0;
    case WM_DESTROY:
        quit = 1;
        PostQuitMessage(0);
        return 0;
    }
    return DefWindowProcA(h, msg, wp, lp);
}

/* Console program (not -mwindows): _popen() does not work in a GUI-subsystem
   program, and the console shows FFmpeg / live11 error messages. */
int main(int argc, char **argv)
{
    HINSTANCE hi = GetModuleHandle(NULL);
    const char *cam = "FullScale_USB3-A3 v3.030", *file = NULL;
    if (argc < 2) {
        fprintf(stderr, "usage: live15.exe nuc2.bin [-o omodel.bin] [-o2 hot_omodel.bin] [-s sffc.bin] [-k col.bin] [-c camera_name] [-i file.raw]\n");
        return 1;
    }
    for (int i = 2; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-c")) cam = argv[i + 1];
        else if (!strcmp(argv[i], "-i")) file = argv[i + 1];
        else if (!strcmp(argv[i], "-o")) { if (load_omodel(argv[i + 1]) < 0) return 1; mode = 5; }
        else if (!strcmp(argv[i], "-o2")) { if (load_omodel(argv[i + 1]) < 0) return 1; mode = 5; }
        else if (!strcmp(argv[i], "-tm")) tab_ord = atoi(argv[i + 1]) - 1;
        else if (!strcmp(argv[i], "-s")) { if (load_sffc(argv[i + 1]) < 0) return 1; }
        else if (!strcmp(argv[i], "-k")) { if (load_coltab(argv[i + 1]) < 0) return 1; }
        else if (!strcmp(argv[i], "-T")) t_fixed = (float)atof(argv[i + 1]);
    }
    if (alloc_all() < 0 || load_nuc(argv[1]) < 0) {
        fprintf(stderr, "cannot load the nuc2 file (wrong name or not made by nuc2 calc)\n");
        return 1;
    }
    if (file) {
        src = fopen(file, "rb");
        from_file = 1;
    } else {
        /* Same FFmpeg capture as cap.bat, raw bytes to stdout */
        char cmd[512];
        snprintf(cmd, sizeof cmd,
                 "ffmpeg -nostdin -hide_banner -loglevel error -f dshow -video_size 640x480 "
                 "-pixel_format yuyv422 -rtbufsize 200M -i video=\"%s\" "
                 "-fps_mode passthrough -f rawvideo -", cam);
        printf("running: %s\n", cmd);
        src = _popen(cmd, "rb");
    }
    if (!src) {
        fprintf(stderr, "cannot open the video source (ffmpeg on PATH? file name?)\n");
        return 1;
    }
    InitializeCriticalSection(&lock);

    WNDCLASSA wc = {0};
    wc.lpfnWndProc = wndproc;
    wc.hInstance = hi;
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.lpszClassName = "live11";
    RegisterClassA(&wc);
    RECT r = { 0, 0, W * 2, H * 2 };       /* 2x zoom, window can be resized */
    AdjustWindowRect(&r, WS_OVERLAPPEDWINDOW, FALSE);
    hwnd = CreateWindowA("live11", "live11", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, CW_USEDEFAULT,
                         r.right - r.left, r.bottom - r.top, NULL, NULL, hi, NULL);
    ShowWindow(hwnd, SW_SHOW);
    CreateThread(NULL, 0, reader, NULL, 0, NULL);

    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0) > 0) { TranslateMessage(&msg); DispatchMessage(&msg); }
    quit = 1;
    return 0;
}

#elif defined(HEADLESS)  /* ------------- Linux headless test ------------- */

int main(int argc, char **argv)
{
    if (argc < 4 || strcmp(argv[2], "-i")) { fprintf(stderr, "usage: live15 nuc2.bin -i file.raw\n"); return 1; }
    if (alloc_all() < 0 || load_nuc(argv[1]) < 0) return 1;
    FILE *f = fopen(argv[3], "rb");
    if (!f) { perror(argv[3]); return 1; }
    uint8_t *raw = malloc(FRAME_BYTES);
    /* reference from the first REF_FRAMES frames */
    ref_req = 1;
    for (int k = 0; k < REF_FRAMES; k++) {
        if (fread(raw, 1, FRAME_BYTES, f) != FRAME_BYTES) { fprintf(stderr, "file too short\n"); return 1; }
        process(raw);
    }
    fseek(f, 0, SEEK_END);
    long nfr = ftell(f) / FRAME_BYTES;
    fseek(f, (nfr > 60 ? 60L : nfr - 1) * FRAME_BYTES, SEEK_SET);
    if (fread(raw, 1, FRAME_BYTES, f) != FRAME_BYTES) { fprintf(stderr, "file too short\n"); return 1; }
    for (int m = 1; m <= 4; m++) {
        mode = m;
        agc_lo = -1;
        process(raw);
        char fn[32];
        snprintf(fn, sizeof fn, "live_mode%d.pgm", m);
        FILE *o = fopen(fn, "wb");
        fprintf(o, "P5\n%d %d\n255\n", W, H);
        fwrite(view8, 1, NPIX, o);
        fclose(o);
        printf("%s  ->  %s\n", status, fn);
    }
    return 0;
}
#else  /* ---------------- Linux: SDL2 window ---------------- */

#include <SDL.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/select.h>
#include <unistd.h>
#include <linux/videodev2.h>
#include <poll.h>
#include <termios.h>
#include <time.h>

/* Keys typed in the terminal work too (not only in the window) */
static struct termios term_old;
static int term_raw = 0;
static void term_restore(void) { if (term_raw) tcsetattr(0, TCSANOW, &term_old); }
static void term_setup(void)
{
    if (!isatty(0) || tcgetattr(0, &term_old) < 0) return;
    struct termios t = term_old;
    t.c_lflag &= ~(ICANON | ECHO);
    t.c_cc[VMIN] = 0; t.c_cc[VTIME] = 0;
    tcsetattr(0, TCSANOW, &t);
    term_raw = 1;
    atexit(term_restore);
}

/* Direct V4L2 capture (same as fps_check.c, which measured 60.03 fps with
 * no loss). Only the format is set; no camera control is written. */
#define NBUF 4
static int vfd = -1;
static void *vmem[NBUF];

static int xioctl(int fd, unsigned long req, void *arg)
{
    int r;
    do r = ioctl(fd, req, arg); while (r < 0 && errno == EINTR);
    return r;
}

static size_t vlen[NBUF];
static int vquiet = 0;                  /* no error prints while reconnecting */
#define VERR(x) do { if (!vquiet) perror(x); } while (0)

static void video_close(void)
{
    if (vfd < 0) return;
    enum v4l2_buf_type t = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    xioctl(vfd, VIDIOC_STREAMOFF, &t);
    for (int i = 0; i < NBUF; i++)
        if (vmem[i] && vmem[i] != MAP_FAILED) { munmap(vmem[i], vlen[i]); vmem[i] = NULL; }
    close(vfd);
    vfd = -1;
}

static int video_open(const char *dev)
{
    vfd = open(dev, O_RDWR);
    if (vfd < 0) { VERR(dev); return -1; }
    struct v4l2_format fmt = {0};
    fmt.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    fmt.fmt.pix.width = W;
    fmt.fmt.pix.height = H;
    fmt.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
    if (xioctl(vfd, VIDIOC_S_FMT, &fmt) < 0) { VERR("S_FMT"); video_close(); return -1; }
    struct v4l2_requestbuffers rb = {0};
    rb.count = NBUF; rb.type = V4L2_BUF_TYPE_VIDEO_CAPTURE; rb.memory = V4L2_MEMORY_MMAP;
    if (xioctl(vfd, VIDIOC_REQBUFS, &rb) < 0) { VERR("REQBUFS"); video_close(); return -1; }
    for (unsigned i = 0; i < rb.count && i < NBUF; i++) {
        struct v4l2_buffer b = {0};
        b.type = V4L2_BUF_TYPE_VIDEO_CAPTURE; b.memory = V4L2_MEMORY_MMAP; b.index = i;
        if (xioctl(vfd, VIDIOC_QUERYBUF, &b) < 0) { VERR("QUERYBUF"); video_close(); return -1; }
        vlen[i] = b.length;
        vmem[i] = mmap(NULL, b.length, PROT_READ | PROT_WRITE, MAP_SHARED, vfd, b.m.offset);
        if (vmem[i] == MAP_FAILED) { VERR("mmap"); video_close(); return -1; }
        if (xioctl(vfd, VIDIOC_QBUF, &b) < 0) { VERR("QBUF"); video_close(); return -1; }
    }
    enum v4l2_buf_type t = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    if (xioctl(vfd, VIDIOC_STREAMON, &t) < 0) { VERR("STREAMON"); video_close(); return -1; }
    return 0;
}

/* Copy the next complete frame into raw. 0 = ok, -1 = stream stopped */
static int video_grab(uint8_t *raw)
{
    for (;;) {
        fd_set fds; FD_ZERO(&fds); FD_SET(vfd, &fds);
        struct timeval tv = { 2, 0 };
        int r = select(vfd + 1, &fds, NULL, NULL, &tv);
        if (r < 0 && errno == EINTR) continue;
        if (r <= 0) return -1;                       /* 2 s without a frame */
        struct v4l2_buffer b = {0};
        b.type = V4L2_BUF_TYPE_VIDEO_CAPTURE; b.memory = V4L2_MEMORY_MMAP;
        if (xioctl(vfd, VIDIOC_DQBUF, &b) < 0) return -1;
        int ok = b.bytesused == FRAME_BYTES && !(b.flags & V4L2_BUF_FLAG_ERROR);
        /* USB diagnostics: frames lost (gap in the driver sequence number) and
           broken buffers are printed with the time, to match them with what
           is seen on the screen */
        {
            static long last_seq = -1, n_drop = 0, n_bad = 0;
            if (last_seq >= 0 && (long)b.sequence > last_seq + 1) {
                n_drop += (long)b.sequence - last_seq - 1;
                time_t tn = time(NULL); char hm[16];
                strftime(hm, sizeof hm, "%H:%M:%S", localtime(&tn));
                printf("%s  USB: %ld frame(s) lost (total lost %ld, broken %ld)\n",
                       hm, (long)b.sequence - last_seq - 1, n_drop, n_bad);
            }
            last_seq = (long)b.sequence;
            if (!ok) {
                n_bad++;
                time_t tn = time(NULL); char hm[16];
                strftime(hm, sizeof hm, "%H:%M:%S", localtime(&tn));
                printf("%s  USB: broken frame skipped (%u bytes, flags 0x%x; total broken %ld)\n",
                       hm, b.bytesused, b.flags, n_bad);
            }
        }
        if (ok) memcpy(raw, vmem[b.index], FRAME_BYTES);
        xioctl(vfd, VIDIOC_QBUF, &b);
        if (ok) return 0;                            /* skip broken frames */
    }
}

/* ---- T70 (USB3-A3 board temperature) over the serial port, read-only ----
 * Same commands as nuc2: replies ON, R70, replies OFF (FullScale UG004). */
static int sfd = -1;
static int ser_open(const char *name)
{
    sfd = open(name, O_RDWR | O_NOCTTY);
    if (sfd < 0) { perror(name); return -1; }
    struct termios t;
    tcgetattr(sfd, &t);
    cfmakeraw(&t);
    cfsetispeed(&t, B115200);
    cfsetospeed(&t, B115200);
    t.c_cflag |= CLOCAL | CREAD;
    t.c_cflag &= ~CRTSCTS;
    tcsetattr(sfd, TCSANOW, &t);
    tcflush(sfd, TCIOFLUSH);
    return 0;
}
static int uart_cmd(const char *c, char *rep, int max)
{
    tcflush(sfd, TCIFLUSH);
    if (write(sfd, c, strlen(c)) != (ssize_t)strlen(c)) return -1;
    tcdrain(sfd);
    int n = 0, first = 1;
    for (;;) {
        struct pollfd p = { .fd = sfd, .events = POLLIN };
        if (poll(&p, 1, first ? 400 : 60) <= 0) break;
        int k = read(sfd, rep + n, max - 1 - n);
        if (k <= 0) break;
        n += k; first = 0;
        if (n >= max - 1) break;
    }
    rep[n] = 0;
    return n;
}
/* Read R70 only. UART replies are switched ON once when the thread starts
   and OFF when it ends (not around every read), to keep traffic to the
   USB3-A3 minimal while it is streaming video. */
static float read_t70(void)
{
    char rep[128];
    unsigned v;
    int n = uart_cmd("R70 0000 0000", rep, sizeof rep);
    const char *p = n > 0 && !strstr(rep, "rror") ? strstr(rep, "R70 0000 ") : NULL;
    /* TMP103: 8-bit two's complement in the low byte (0xF7 = -9 C) */
    return (p && sscanf(p + 9, "%4x", &v) == 1) ? (float)(int8_t)(v & 0xFF) : NAN;
}

static SDL_mutex *mtx;
static volatile int quit;
static const char *vdev = "/dev/v4l/by-id/usb-FullScale_FullScale_USB3-A3_v3.030-video-index0";
static const char *t70_port = "/dev/ttyACM0";
#define T70_PERIOD_S 10     /* T70 changes slowly (1 C steps); read every 10 s */
static int t70_thread(void *arg)
{
    (void)arg;
    char rep[64];
    uart_cmd("W01 0009 0001", rep, sizeof rep);              /* UART replies ON (once) */
    int fails = 0;
    while (!quit) {
        float t = read_t70();
        if (!isnan(t)) { t_board = t; fails = 0; }
        else if (++fails >= 3) {
            /* serial port gone (USB disconnect): reopen when it is back */
            if (sfd >= 0) { close(sfd); sfd = -1; }
            t_board = NAN;                       /* mode 5 waits, falls back to 3 */
            while (!quit) {
                int fd = open(t70_port, O_RDWR | O_NOCTTY);
                if (fd >= 0) { close(fd); if (ser_open(t70_port) == 0) break; }
                SDL_Delay(1000);
            }
            if (!quit) uart_cmd("W01 0009 0001", rep, sizeof rep);   /* replies ON again */
            fails = 0;
            continue;
        }
        for (int k = 0; k < T70_PERIOD_S * 10 && !quit; k++) SDL_Delay(100);
    }
    uart_cmd("W01 0009 0000", rep, sizeof rep);              /* UART replies OFF */
    return 0;
}
static void save_snapshot(void);

static void handle_key(int k)
{
    show_status = 1;
    if (k >= 'A' && k <= 'Z') k += 'a' - 'A';
    if (k >= '1' && k <= '5') mode = k - '0';
    else if (k == 'a') { autot = !autot; }
    else if (k == 'f') sffc_on = !sffc_on;
    else if (k == 'k') col_on = !col_on;
    else if (k == '[') t_trim -= 0.5f;
    else if (k == ']') t_trim += 0.5f;
    else if (k == 'r') { ref_count = 0; ref_req = 1; }
    else if (k == 'c') colfilt = !colfilt;
    else if (k == 'h' || k == 'g') agcm = (agcm + 1) % 4;
    else if (k == '+' || k == '=') { if (gmax < 6) gmax += 0.5f; }
    else if (k == '-') { if (gmax > 1) gmax -= 0.5f; }
    else if (k == 't') tfilt = !tfilt;
    else if (k == 's') { SDL_LockMutex(mtx); save_snapshot(); SDL_UnlockMutex(mtx); }
    else if (k == ' ') frozen = !frozen;
    else if (k == 'q' || k == 27) quit = 1;
}
static uint8_t shown[NPIX];
static FILE *src;
static int from_file;
static volatile int new_frame;
static char title[460];

static int reader(void *arg)
{
    (void)arg;
    int lost = 0;
    Uint32 lost_t = 0;
    uint8_t *raw = malloc(FRAME_BYTES);
    Uint32 t0 = SDL_GetTicks();
    int frames = 0;
    double fps = 0;
    while (!quit) {
        if (from_file) {
            if (fread(raw, 1, FRAME_BYTES, src) != FRAME_BYTES) { fseek(src, 0, SEEK_SET); continue; }
        } else if (video_grab(raw) < 0) {
            /* Camera lost (USB disconnect / no frames). Keep the window and
               try to open the device again every second; the stream resumes
               by itself when the camera is back (same -d path). Reported
               once, not every retry. */
            if (!lost) {
                time_t tl = time(NULL);
                fprintf(stderr, "camera lost (no frames) at %s", ctime(&tl));
                lost = 1; lost_t = SDL_GetTicks();
            }
            video_close();
            SDL_LockMutex(mtx);
            snprintf(title, sizeof title, "live11 - CAMERA LOST - no frames from %s for %u s - replug if it stays",
                     vdev, (SDL_GetTicks() - lost_t) / 1000);
            new_frame = 1;
            SDL_UnlockMutex(mtx);
            vquiet = 1;
            while (!quit && video_open(vdev) < 0) SDL_Delay(1000);
            vquiet = 0;
            tprev_ok = 0;
            continue;
        }
        if (lost && !from_file) {
            time_t tl = time(NULL);
            fprintf(stderr, "camera back after %u s at %s", (SDL_GetTicks() - lost_t) / 1000, ctime(&tl));
            lost = 0;
        }
        if (from_file) SDL_Delay(16);
        frames++;
        Uint32 t = SDL_GetTicks();
        if (t - t0 >= 1000) { fps = frames * 1000.0 / (t - t0); frames = 0; t0 = t; }
        if (frozen && !ref_req) continue;
        SDL_LockMutex(mtx);
        process(raw);
        memcpy(shown, view8, NPIX);
        snprintf(title, sizeof title, "live15 | %.1f fps", fps);
        new_frame = 1;
        SDL_UnlockMutex(mtx);
    }
    free(raw);
    return 0;
}

int main(int argc, char **argv)
{
    /* default camera: the FullScale board by its USB id (a laptop webcam can
       take /dev/video0); /dev/video0 only if the by-id link is missing */
    const char *byid = "/dev/v4l/by-id/usb-FullScale_FullScale_USB3-A3_v3.030-video-index0";
    const char *dev = access(byid, F_OK) == 0 ? byid : "/dev/video0", *file = NULL, *port = "/dev/ttyACM0";
    if (argc < 2) {
        fprintf(stderr, "usage: live15 nuc2.bin [-o omodel.bin] [-o2 hot_omodel.bin] [-s sffc.bin] [-k col.bin] [-d /dev/videoN] [-p /dev/ttyACM0|none] [-T degC] [-i file.raw]\n");
        return 1;
    }
    for (int i = 2; i + 1 < argc; i += 2) {
        if (!strcmp(argv[i], "-d")) dev = vdev = argv[i + 1];
        else if (!strcmp(argv[i], "-i")) file = argv[i + 1];
        else if (!strcmp(argv[i], "-p")) port = t70_port = argv[i + 1];
        else if (!strcmp(argv[i], "-T")) t_fixed = (float)atof(argv[i + 1]);
        else if (!strcmp(argv[i], "-o")) { if (load_omodel(argv[i + 1]) < 0) return 1; mode = 5; }
        else if (!strcmp(argv[i], "-o2")) { if (load_omodel(argv[i + 1]) < 0) return 1; mode = 5; }
        else if (!strcmp(argv[i], "-tm")) tab_ord = atoi(argv[i + 1]) - 1;
        else if (!strcmp(argv[i], "-s")) { if (load_sffc(argv[i + 1]) < 0) return 1; }
        else if (!strcmp(argv[i], "-k")) { if (load_coltab(argv[i + 1]) < 0) return 1; }
    }
    if (alloc_all() < 0 || load_nuc(argv[1]) < 0) return 1;
    vdev = dev;
    if (!file) printf("camera: %s\n", dev);
    if (file) {
        src = fopen(file, "rb");
        from_file = 1;
        if (!src) { perror(file); return 1; }
    } else if (video_open(dev) < 0) {
        /* camera not there yet (just plugged / re-enumerating): wait for it */
        fprintf(stderr, "waiting for %s ... (plug the camera, Ctrl+C to stop)\n", dev);
        vquiet = 1;
        while (video_open(dev) < 0) sleep(1);
        vquiet = 0;
        fprintf(stderr, "camera found\n");
    }

    term_setup();
    printf("keys (window or this terminal): 1-5 mode ([ ] trim T, a auto T, f SFFC, k columns), r ref, c column, t temporal, h HE, s snap, q quit\n");
    if (SDL_Init(SDL_INIT_VIDEO) < 0) { fprintf(stderr, "SDL: %s\n", SDL_GetError()); return 1; }
    SDL_Window *win = SDL_CreateWindow("live11", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       W * 2, H * 2, SDL_WINDOW_RESIZABLE);
    SDL_Renderer *ren = SDL_CreateRenderer(win, -1, 0);
    SDL_Texture *tex = SDL_CreateTexture(ren, SDL_PIXELFORMAT_ARGB8888,
                                         SDL_TEXTUREACCESS_STREAMING, W, H);
    static uint32_t argb[NPIX];
    mtx = SDL_CreateMutex();
    SDL_Thread *th = SDL_CreateThread(reader, "reader", NULL);
    SDL_Thread *tth = NULL;
    if (om_ok && isnan(t_fixed) && strcmp(port, "none") && ser_open(port) == 0)
        tth = SDL_CreateThread(t70_thread, "t70", NULL);
    else if (om_ok && isnan(t_fixed))
        fprintf(stderr, "no temperature: give -p /dev/ttyACM0 or -T <degC>\n");

    while (!quit) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_QUIT) quit = 1;
            else if (e.type == SDL_KEYDOWN) {
                SDL_Keycode k = e.key.keysym.sym;   /* ASCII keys have their ASCII code */
                if (k < 128) handle_key((int)k);
            }
        }
        /* keys typed in the terminal */
        char tk;
        while (term_raw && read(0, &tk, 1) == 1) handle_key((unsigned char)tk);
        if (new_frame) {
            SDL_LockMutex(mtx);
            for (int i = 0; i < NPIX; i++)
                argb[i] = 0xFF000000u | (uint32_t)shown[i] << 16 | (uint32_t)shown[i] << 8 | shown[i];
            SDL_SetWindowTitle(win, title);
            new_frame = 0;
            SDL_UnlockMutex(mtx);
            SDL_UpdateTexture(tex, NULL, argb, W * 4);
            SDL_RenderClear(ren);
            SDL_RenderCopy(ren, tex, NULL, NULL);
            SDL_RenderPresent(ren);
        }
        SDL_Delay(5);
    }
    SDL_WaitThread(th, NULL);
    if (tth) SDL_WaitThread(tth, NULL);
    if (from_file) fclose(src);
    else { enum v4l2_buf_type t = V4L2_BUF_TYPE_VIDEO_CAPTURE; xioctl(vfd, VIDIOC_STREAMOFF, &t); close(vfd); }
    SDL_Quit();
    return 0;
}
#endif
