/*
 * nuc2.c - 2-point NUC (Non-Uniformity Correction) tool for the ATTO640D-04W
 *          (FullScale K613-04D proxy + USB3-A3 interface).
 *
 * Works on Windows (MSYS2 / MinGW gcc) and Linux (gcc).
 * It never writes any sensor / proxy setting. On the UART it only sends
 * read commands, plus "W01 0009 0001/0000" (enable / disable UART replies,
 * as described in FullScale UG004).
 *
 * ---------------------------------------------------------------------------
 * Theory (Lynred NUC training, slides 16-18)
 *   I1 = average image in front of the cold target (T1)
 *   I2 = average image in front of the hot target  (T2)
 *   M1 = mean(I1), M2 = mean(I2)                   (over good pixels)
 *   G  = (M2 - M1) / (I2 - I1)                     per pixel  -> gain map
 *   O  = M1 - G * I1                               per pixel  -> offset map
 *   corrected = raw * G + O
 * Gain depends on GFID / CINT / TINT (slide 19), so the sensor settings are
 * stored in every map and "calc" refuses to mix maps with different settings.
 *
 * ---------------------------------------------------------------------------
 * Commands
 *   map     : average frames from a raw YUY2 file -> .map file   (Windows + Linux)
 *             nuc2 map cold1.raw -o cold1.map -t 25.0 [-n 128] [-s 30] [-p COM5|/dev/ttyACM0|none]
 *   capture : grab frames directly with V4L2 -> .map file        (Linux only)
 *             nuc2 capture -o cold1.map -t 25.0 [-n 128] [-d /dev/video0] [-p /dev/ttyACM0|none]
 *   calc    : cold + hot (+ second cold) maps -> nuc2.bin (G, O, bad pixels)
 *             nuc2 calc cold1.map hot.map [cold2.map] -o nuc2.bin
 *   check   : residual FPN of a test map: raw vs 1-point vs 2-point (slide 20)
 *             nuc2 check nuc2.bin test.map [ref.map]
 *             ref.map = uniform reference (e.g. 25 C blackbody) taken just before the
 *             test: gives 1-point and 2-point with refreshed offset (slide 18)
 *   apply   : correct one frame of a recorded scene, write a 2x2 comparison image
 *             nuc2 apply nuc2.bin scene.raw [-f 30] [-r ref.map]
 *             -r = uniform reference (lens cap / plate) taken just before the scene,
 *                  used for 1-point and for the offset refresh of slide 18
 *
 *   -t  target temperature in deg C (read it with a thermometer)
 *   -n  number of frames to average (default 128)
 *   -s  number of frames to skip at the start of the raw file (default 30)
 *   -p  UART port to read the sensor settings (TINT, GFID, GSK, GAIN, board temp)
 *
 * Raw file format (for "map"): consecutive 640x480 YUY2 frames, 614400 bytes each,
 * e.g. recorded on Windows with:
 *   ffmpeg -f dshow -video_size 640x480 -pixel_format yuyv422 -rtbufsize 200M
 *          -i video="FullScale_USB3-A3 v3.030" -frames:v 160 -f rawvideo cold1.raw
 * Pixel value = ((byte0 << 8) | byte1) >> 2   (FullScale UG004 section 9.1)
 *
 * Build
 *   Windows (MSYS2 UCRT64 shell):  gcc -O2 -Wall -o nuc2.exe nuc2.c -lm
 *   Linux:                         gcc -O2 -Wall -o nuc2 nuc2.c -lm
 *
 * Output images (PGM, open with GIMP / IrfanView / eog)
 *   map / capture : <name>.map.pgm   averaged image
 *   calc          : nuc2_gain.pgm    gain map
 *                   nuc2_resp.pgm    response I2 - I1 (target must fill the view)
 *                   nuc2_bad.pgm     bad pixels in white
 *   apply         : apply_all.pgm    raw | 2-pt / 1-pt | 2-pt + offset refresh
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <poll.h>
#include <termios.h>
#include <unistd.h>
#endif

#ifdef __linux__
#include <errno.h>
#include <linux/videodev2.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <sys/select.h>
#endif

#define W 640
#define H 480
#define NPIX (W * H)
#define FRAME_BYTES (NPIX * 2)

/* Bad pixel rules */
#define DEAD_DL   1000      /* average below this -> dead pixel           */
#define HOT_DL    15000     /* average above this -> hot / saturated      */
#define RESP_LO   0.5       /* response < 0.5 x median response -> weak   */
#define RESP_HI   2.0       /* response > 2.0 x median response -> bad    */
#define G_SIGMA   6.0       /* gain further than 6 sigma from mean -> bad */

/* ======================================================================
 * File formats (same layout on Windows and Linux 64-bit)
 * ====================================================================== */
struct settings {                     /* sensor state when a map was taken */
    int32_t tint, gfid, gsk, gain_reg; /* -1 = not read                     */
    float   t70_c;                     /* USB3-A3 board temp, NAN = not read */
};

struct maphdr {
    char     magic[8];                 /* "TVEMAP1"                          */
    uint32_t w, h, nframes;
    float    target_c;                 /* target temperature given by user   */
    float    noise_dl;                 /* temporal noise per frame (DL)      */
    int64_t  unix_time;
    struct settings s;
};

struct nuchdr {
    char     magic[8];                 /* "TVENUC2"                          */
    uint32_t w, h;
    float    m1, m2, t1_c, t2_c;
    int32_t  nbad;
    struct settings s;
    /* followed by: float G[NPIX], float O[NPIX], uint8_t bad[NPIX] */
};

/* ======================================================================
 * Small helpers
 * ====================================================================== */
static int write_pgm(const char *fn, const float *v, const uint8_t *mask, double lo, double hi)
{
    FILE *f = fopen(fn, "wb");
    if (!f) { perror(fn); return -1; }
    fprintf(f, "P5\n%d %d\n255\n", W, H);
    for (int i = 0; i < NPIX; i++) {
        double x = (mask && mask[i]) ? 0 : (v[i] - lo) / (hi - lo) * 255.0;
        fputc(x < 0 ? 0 : x > 255 ? 255 : (int)x, f);
    }
    fclose(f);
    return 0;
}

static int cmp_f(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

/* Percentile p (0..100) over good pixels */
static double pct(const float *v, const uint8_t *bad, double p)
{
    float *t = malloc(NPIX * sizeof *t);
    int n = 0;
    for (int i = 0; i < NPIX; i++) if (!bad || !bad[i]) t[n++] = v[i];
    qsort(t, n, sizeof *t, cmp_f);
    double r = n ? t[(int)(p / 100.0 * (n - 1))] : 0;
    free(t);
    return r;
}

/* Mean and standard deviation over good pixels (two passes for precision) */
static void stats(const float *v, const uint8_t *bad, double *mean, double *sd)
{
    double s = 0, s2 = 0; long n = 0;
    for (int i = 0; i < NPIX; i++) {
        if (bad && bad[i]) continue;
        s += v[i]; n++;
    }
    *mean = s / n;
    for (int i = 0; i < NPIX; i++) {
        if (bad && bad[i]) continue;
        double d = v[i] - *mean;
        s2 += d * d;
    }
    *sd = sqrt(s2 / n);
}

static void print_settings(const char *tag, const struct settings *s)
{
    printf("%s TINT=%d  GFID=0x%02X  GSK=0x%03X  GAIN_REG=0x%02X  T70=%.0f C\n", tag,
           s->tint, s->gfid & 0xFF, s->gsk & 0x3FF, s->gain_reg & 0xFF, s->t70_c);
}

/* 64-bit file size / seek (raw files can be large) */
static long long file_size(FILE *f)
{
#ifdef _WIN32
    _fseeki64(f, 0, SEEK_END);
    long long s = _ftelli64(f);
    _fseeki64(f, 0, SEEK_SET);
#else
    fseeko(f, 0, SEEK_END);
    long long s = ftello(f);
    fseeko(f, 0, SEEK_SET);
#endif
    return s;
}

static void file_seek(FILE *f, long long pos)
{
#ifdef _WIN32
    _fseeki64(f, pos, SEEK_SET);
#else
    fseeko(f, (off_t)pos, SEEK_SET);
#endif
}

/* Decode one YUY2 frame (FullScale UG004 9.1): Y byte = 8 MSB, U/V byte = 6 LSB */
static void decode_yuy2(const uint8_t *raw, uint16_t *img)
{
    for (int i = 0; i < NPIX; i++)
        img[i] = (uint16_t)(((raw[2 * i] << 8) | raw[2 * i + 1]) >> 2);
}

/* ======================================================================
 * UART (read only). Windows: COMx, Linux: /dev/ttyACMx
 * ====================================================================== */
#ifdef _WIN32
static HANDLE hser = INVALID_HANDLE_VALUE;

static int ser_open(const char *name)
{
    char path[64];
    snprintf(path, sizeof path, "\\\\.\\%s", name);          /* works for COM10+ too */
    hser = CreateFileA(path, GENERIC_READ | GENERIC_WRITE, 0, NULL, OPEN_EXISTING, 0, NULL);
    if (hser == INVALID_HANDLE_VALUE) {
        fprintf(stderr, "cannot open %s (error %lu)\n", name, GetLastError());
        return -1;
    }
    DCB dcb = { .DCBlength = sizeof dcb };
    GetCommState(hser, &dcb);
    dcb.BaudRate = CBR_115200;                               /* 115200 8N1, no flow ctrl */
    dcb.ByteSize = 8;
    dcb.Parity = NOPARITY;
    dcb.StopBits = ONESTOPBIT;
    dcb.fOutxCtsFlow = FALSE;
    dcb.fOutxDsrFlow = FALSE;
    dcb.fDtrControl = DTR_CONTROL_ENABLE;
    dcb.fRtsControl = RTS_CONTROL_ENABLE;
    SetCommState(hser, &dcb);
    COMMTIMEOUTS to = {0};
    to.ReadIntervalTimeout = 60;                             /* reply finished after 60 ms gap */
    to.ReadTotalTimeoutConstant = 400;                       /* wait max 400 ms for a reply    */
    SetCommTimeouts(hser, &to);
    PurgeComm(hser, PURGE_RXCLEAR | PURGE_TXCLEAR);
    return 0;
}
static void ser_flush_rx(void) { PurgeComm(hser, PURGE_RXCLEAR); }
static int ser_write(const char *s, int n) { DWORD w = 0; WriteFile(hser, s, n, &w, NULL); return (int)w; }
static int ser_read(char *b, int max) { DWORD r = 0; if (!ReadFile(hser, b, max, &r, NULL)) return -1; return (int)r; }
static void ser_close(void) { CloseHandle(hser); hser = INVALID_HANDLE_VALUE; }
#else
static int sfd = -1;

static int ser_open(const char *name)
{
    sfd = open(name, O_RDWR | O_NOCTTY);
    if (sfd < 0) { perror(name); return -1; }
    struct termios t;
    tcgetattr(sfd, &t);
    cfmakeraw(&t);                                           /* raw mode, no echo */
    cfsetispeed(&t, B115200);
    cfsetospeed(&t, B115200);
    t.c_cflag |= CLOCAL | CREAD;
    t.c_cflag &= ~CRTSCTS;
    tcsetattr(sfd, TCSANOW, &t);
    tcflush(sfd, TCIOFLUSH);
    return 0;
}
static void ser_flush_rx(void) { tcflush(sfd, TCIFLUSH); }
static int ser_write(const char *s, int n) { int w = write(sfd, s, n); tcdrain(sfd); return w; }
static int ser_read(char *b, int max)
{
    static int first = 1;
    struct pollfd p = { .fd = sfd, .events = POLLIN };
    if (poll(&p, 1, first ? 400 : 60) <= 0) { first = 1; return 0; }
    int k = read(sfd, b, max);
    first = (k <= 0);
    return k;
}
static void ser_close(void) { close(sfd); sfd = -1; }
#endif

/* Send one 13-character command and collect the reply text */
static int uart_cmd(const char *c, char *rep, int max)
{
    ser_flush_rx();
    if (ser_write(c, (int)strlen(c)) != (int)strlen(c)) return -1;
    int n = 0;
    for (;;) {
        int k = ser_read(rep + n, max - 1 - n);
        if (k <= 0) break;
        n += k;
        if (n >= max - 1) break;
    }
    rep[n] = 0;
    for (int i = 0; i < n; i++) if (rep[i] < 0x20) rep[i] = ' ';
    return n;
}

/* Read a register; returns the 16-bit value or -1 on error / no reply */
static long uart_read(const char *c)
{
    char rep[128], prefix[10];
    if (uart_cmd(c, rep, sizeof rep) <= 0 || strstr(rep, "rror")) return -1;
    memcpy(prefix, c, 9); prefix[9] = 0;                     /* e.g. "R12 004F " */
    const char *p = strstr(rep, prefix);
    unsigned v;
    if (!p || sscanf(p + 9, "%4x", &v) != 1) return -1;
    return v;
}

/* Read TINT, GFID, GSK, GAIN_IMAGE (sensor, I2C dest 12) and board temp (R70) */
static void read_settings(const char *port, struct settings *s)
{
    s->tint = s->gfid = s->gsk = s->gain_reg = -1;
    s->t70_c = NAN;
    if (!port || !strcmp(port, "none")) return;
    if (ser_open(port) < 0) return;
    char rep[64];
    uart_cmd("W01 0009 0001", rep, sizeof rep);              /* UART replies ON  */
    long a = uart_read("R12 004F 0000"), b = uart_read("R12 0050 0000");
    if (a >= 0 && b >= 0) s->tint = (int32_t)((a << 8) | b); /* TINT[13:0]       */
    s->gfid = (int32_t)uart_read("R12 004B 0000");          /* DAC_GFID         */
    a = uart_read("R12 004C 0000"); b = uart_read("R12 004D 0000");
    if (a >= 0 && b >= 0) s->gsk = (int32_t)(((a & 3) << 8) | b); /* DAC_GSK[9:0] */
    s->gain_reg = (int32_t)uart_read("R12 0040 0000");      /* GAIN_IMAGE       */
    long t70 = uart_read("R70 0000 0000");
    if (t70 >= 0) s->t70_c = (int16_t)t70;                   /* deg C (UG004 p.60) */
    uart_cmd("W01 0009 0000", rep, sizeof rep);              /* UART replies OFF */
    ser_close();
    if (s->tint < 0) fprintf(stderr, "WARNING: no reply on %s - settings not recorded\n", port);
}

/* ======================================================================
 * Map file I/O
 * ====================================================================== */
static int save_map(const char *fn, const struct maphdr *h, const float *v)
{
    FILE *f = fopen(fn, "wb");
    if (!f) { perror(fn); return -1; }
    fwrite(h, sizeof *h, 1, f);
    fwrite(v, sizeof *v, NPIX, f);
    fclose(f);
    return 0;
}

static int load_map(const char *fn, struct maphdr *h, float *v)
{
    FILE *f = fopen(fn, "rb");
    if (!f) { perror(fn); return -1; }
    int ok = fread(h, sizeof *h, 1, f) == 1 && !memcmp(h->magic, "TVEMAP1", 7) &&
             h->w == W && h->h == H && fread(v, sizeof *v, NPIX, f) == NPIX;
    fclose(f);
    if (!ok) fprintf(stderr, "%s: not a valid map file\n", fn);
    return ok ? 0 : -1;
}

/* Common end of "map" and "capture": statistics, save .map and preview .pgm */
static int finish_map(const char *out, struct maphdr *h, const double *acc,
                      const double *acc2, int n, int dups)
{
    float *avg = malloc(NPIX * sizeof *avg);
    double tn = 0;
    for (int i = 0; i < NPIX; i++) {
        double m = acc[i] / n, var = acc2[i] / n - m * m;
        avg[i] = (float)m;
        tn += var > 0 ? sqrt(var) : 0;
    }
    h->nframes = n;
    h->noise_dl = (float)(tn / NPIX);
    h->unix_time = (int64_t)time(NULL);

    double mean, sd, q[4] = {0};
    stats(avg, NULL, &mean, &sd);
    for (int y = 0; y < H; y++)                              /* quadrant means */
        for (int x = 0; x < W; x++) q[(y >= H / 2) * 2 + (x >= W / 2)] += avg[y * W + x];
    for (int k = 0; k < 4; k++) q[k] /= NPIX / 4;

    printf("\nframes averaged %d, duplicate frames %d\n", n, dups);
    printf("mean %.1f DL   spatial FPN %.1f DL   temporal noise %.2f DL/frame\n",
           mean, sd, h->noise_dl);
    printf("quadrants  TL %.0f  TR %.0f  BL %.0f  BR %.0f\n", q[0], q[1], q[2], q[3]);
    if (mean < 2000 || mean > 14000)
        printf("WARNING: mean outside 2000..14000 - saturation or wrong target?\n");
    if (h->noise_dl < 0.5)
        printf("WARNING: temporal noise very low - frames may be duplicated / not live\n");
    if (dups > n / 10)
        printf("WARNING: many duplicate frames - check the recording\n");

    if (save_map(out, h, avg) < 0) return 1;
    char pgm[512];
    snprintf(pgm, sizeof pgm, "%s.pgm", out);
    write_pgm(pgm, avg, NULL, pct(avg, NULL, 1), pct(avg, NULL, 99));
    printf("saved %s  (+ %s)\n", out, pgm);
    free(avg);
    return 0;
}

/* Parse the common options of "map" and "capture" */
struct opts { const char *out, *dev, *port; int n, skip; float target; };

static int parse_opts(int argc, char **argv, int first, struct opts *o)
{
    for (int i = first; i < argc; i++) {
        const char *a = argv[i];
        const char *v = (i + 1 < argc) ? argv[i + 1] : NULL;
        if      (!strcmp(a, "-o") && v) { o->out = v; i++; }
        else if (!strcmp(a, "-t") && v) {
            /* Reject a non-numeric temperature (e.g. a file name shifted into
               this position by a wrong command line) */
            char *end;
            o->target = strtof(v, &end);
            if (end == v || *end) { fprintf(stderr, "-t %s: not a temperature\n", v); exit(1); }
            i++;
        }
        else if (!strcmp(a, "-n") && v) { o->n = atoi(v); i++; }
        else if (!strcmp(a, "-s") && v) { o->skip = atoi(v); i++; }
        else if (!strcmp(a, "-d") && v) { o->dev = v; i++; }
        else if (!strcmp(a, "-p") && v) { o->port = v; i++; }
        else { fprintf(stderr, "unknown option: %s\n", a); return -1; }
    }
    if (!o->out || isnan(o->target) || o->n < 16) return -1;
    return 0;
}

/* ======================================================================
 * map: average frames from a raw YUY2 file
 * ====================================================================== */
static int cmd_map(int argc, char **argv)
{
    if (argc < 3) goto usage;
    const char *raw = argv[2];
    struct opts o = { .n = 128, .skip = 30, .target = NAN, .port = "none" };
    if (parse_opts(argc, argv, 3, &o) < 0) goto usage;

    FILE *f = fopen(raw, "rb");
    if (!f) { perror(raw); return 1; }
    long long size = file_size(f);
    long long nfile = size / FRAME_BYTES;
    printf("%s: %lld bytes = %lld frames%s\n", raw, size, nfile,
           size % FRAME_BYTES ? "  (WARNING: size is not a multiple of 614400)" : "");
    if (nfile < o.skip + o.n) {
        fprintf(stderr, "ERROR: need %d frames (skip %d + average %d), file has %lld\n",
                o.skip + o.n, o.skip, o.n, nfile);
        return 1;
    }

    struct maphdr h = { .magic = "TVEMAP1", .w = W, .h = H, .target_c = o.target };
    read_settings(o.port, &h.s);
    print_settings("settings:", &h.s);

    uint8_t *buf = malloc(FRAME_BYTES), *prev = calloc(FRAME_BYTES, 1);
    uint16_t *img = malloc(NPIX * sizeof *img);
    double *acc = calloc(NPIX, sizeof *acc), *acc2 = calloc(NPIX, sizeof *acc2);
    if (!buf || !prev || !img || !acc || !acc2) { fprintf(stderr, "out of memory\n"); return 1; }

    file_seek(f, (long long)o.skip * FRAME_BYTES);           /* skip the first frames */
    int dups = 0;
    for (int k = 0; k < o.n; k++) {
        if (fread(buf, 1, FRAME_BYTES, f) != FRAME_BYTES) { fprintf(stderr, "read error\n"); return 1; }
        if (k && !memcmp(buf, prev, FRAME_BYTES)) dups++;    /* identical frame = not live */
        memcpy(prev, buf, FRAME_BYTES);
        decode_yuy2(buf, img);
        for (int i = 0; i < NPIX; i++) { acc[i] += img[i]; acc2[i] += (double)img[i] * img[i]; }
    }
    fclose(f);
    return finish_map(o.out, &h, acc, acc2, o.n, dups);

usage:
    fprintf(stderr, "usage: nuc2 map file.raw -o file.map -t <deg C> [-n 128] [-s 30] [-p COM5|none]\n");
    return 1;
}

/* ======================================================================
 * capture: V4L2 direct capture (Linux only)
 * ====================================================================== */
#ifdef __linux__
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
    if (vfd < 0) { if (!vquiet) perror(dev); return -1; }
    struct v4l2_format f = {0};
    f.type = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    f.fmt.pix.width = W; f.fmt.pix.height = H;
    f.fmt.pix.pixelformat = V4L2_PIX_FMT_YUYV;
    if (xioctl(vfd, VIDIOC_S_FMT, &f) < 0) { if (!vquiet) perror("S_FMT"); video_close(); return -1; }
    struct v4l2_requestbuffers rb = {0};
    rb.count = NBUF; rb.type = V4L2_BUF_TYPE_VIDEO_CAPTURE; rb.memory = V4L2_MEMORY_MMAP;
    if (xioctl(vfd, VIDIOC_REQBUFS, &rb) < 0) { if (!vquiet) perror("REQBUFS"); video_close(); return -1; }
    for (unsigned i = 0; i < rb.count && i < NBUF; i++) {
        struct v4l2_buffer b = {0};
        b.type = V4L2_BUF_TYPE_VIDEO_CAPTURE; b.memory = V4L2_MEMORY_MMAP; b.index = i;
        if (xioctl(vfd, VIDIOC_QUERYBUF, &b) < 0) { if (!vquiet) perror("QUERYBUF"); video_close(); return -1; }
        vlen[i] = b.length;
        vmem[i] = mmap(NULL, b.length, PROT_READ | PROT_WRITE, MAP_SHARED, vfd, b.m.offset);
        if (vmem[i] == MAP_FAILED) { if (!vquiet) perror("mmap"); video_close(); return -1; }
        if (xioctl(vfd, VIDIOC_QBUF, &b) < 0) { if (!vquiet) perror("QBUF"); video_close(); return -1; }
    }
    enum v4l2_buf_type t = V4L2_BUF_TYPE_VIDEO_CAPTURE;
    if (xioctl(vfd, VIDIOC_STREAMON, &t) < 0) { if (!vquiet) perror("STREAMON"); video_close(); return -1; }
    return 0;
}

/* One frame into raw (may be NULL). 1 = good frame, 0 = bad frame, -1 = error */
static int grab(uint8_t *raw)
{
    fd_set fds; FD_ZERO(&fds); FD_SET(vfd, &fds);
    struct timeval tv = { 2, 0 };
    if (select(vfd + 1, &fds, NULL, NULL, &tv) <= 0) {
        if (!vquiet) fprintf(stderr, "TIMEOUT: no frame in 2 s (unplug / replug the camera)\n");
        return -1;
    }
    struct v4l2_buffer b = {0};
    b.type = V4L2_BUF_TYPE_VIDEO_CAPTURE; b.memory = V4L2_MEMORY_MMAP;
    if (xioctl(vfd, VIDIOC_DQBUF, &b) < 0) { perror("DQBUF"); return -1; }
    int ok = b.bytesused == FRAME_BYTES && !(b.flags & V4L2_BUF_FLAG_ERROR);
    if (ok && raw) memcpy(raw, vmem[b.index], FRAME_BYTES);
    if (xioctl(vfd, VIDIOC_QBUF, &b) < 0) { perror("QBUF"); return -1; }
    return ok;
}

static int cmd_capture(int argc, char **argv)
{
    struct opts o = { .n = 128, .skip = 30, .target = NAN,
                      .dev = "/dev/video0", .port = "/dev/ttyACM0" };
    if (parse_opts(argc, argv, 2, &o) < 0) {
        fprintf(stderr, "usage: nuc2 capture -o file.map -t <deg C> [-n 128] [-d /dev/video0] [-p tty|none]\n");
        return 1;
    }
    struct maphdr h = { .magic = "TVEMAP1", .w = W, .h = H, .target_c = o.target };
    read_settings(o.port, &h.s);
    print_settings("settings:", &h.s);
    if (video_open(o.dev) < 0) return 1;

    uint8_t *buf = malloc(FRAME_BYTES), *prev = calloc(FRAME_BYTES, 1);
    uint16_t *img = malloc(NPIX * sizeof *img);
    double *acc = calloc(NPIX, sizeof *acc), *acc2 = calloc(NPIX, sizeof *acc2);
    if (!buf || !prev || !img || !acc || !acc2) { fprintf(stderr, "out of memory\n"); return 1; }

    for (int i = 0; i < o.skip; i++) if (grab(NULL) < 0) return 1;   /* drop old buffers */
    printf("capturing %d frames - keep the target still ...\n", o.n);
    int dups = 0;
    for (int got = 0; got < o.n; ) {
        int r = grab(buf);
        if (r < 0) return 1;
        if (!r) continue;
        if (got && !memcmp(buf, prev, FRAME_BYTES)) dups++;
        memcpy(prev, buf, FRAME_BYTES);
        decode_yuy2(buf, img);
        for (int i = 0; i < NPIX; i++) { acc[i] += img[i]; acc2[i] += (double)img[i] * img[i]; }
        got++;
    }
    return finish_map(o.out, &h, acc, acc2, o.n, dups);
}

/* For long recordings (warmup / chamber): if the camera stops or drops off
 * USB, close it, wait until the device is back, reopen, and carry on.
 * Returns like grab(): 1 good frame, 0 no frame (caller just loops). */
static int reconnects = 0;
static int grab_keep(const char *dev, uint8_t *raw)
{
    int r = grab(raw);
    if (r >= 0) return r;
    time_t t = time(NULL);
    fprintf(stderr, "\ncamera lost at %s  waiting for %s ...\n", ctime(&t), dev);
    video_close();
    vquiet = 1;
    while (video_open(dev) < 0) sleep(1);
    for (int k = 0; k < 30; k++) if (grab(NULL) < 0) break;   /* settle */
    vquiet = 0;
    t = time(NULL);
    fprintf(stderr, "camera back at %s", ctime(&t));
    reconnects++;
    return 0;
}

/* ----------------------------------------------------------------------
 * warmup: record how the offset changes while the camera warms up.
 *   nuc2 warmup -o warm [-i 30] [-m 60] [-n 64] [-d /dev/video0] [-p /dev/ttyACM0|none]
 * Start with a COLD camera (unplugged for 30+ min) looking at a uniform,
 * stable target (lens cap or blackbody). Plug in and start immediately.
 * Every -i seconds for -m minutes it averages -n frames and saves
 * warm_000.map, warm_001.map ... (board temperature T70 in each header)
 * and one line per map in warm.csv. Read-only: only R70 is read on the UART.
 * ---------------------------------------------------------------------- */
static float read_t70(const char *port)
{
    if (!port || !strcmp(port, "none") || ser_open(port) < 0) return NAN;
    char rep[64];
    uart_cmd("W01 0009 0001", rep, sizeof rep);              /* UART replies ON  */
    long t = uart_read("R70 0000 0000");
    uart_cmd("W01 0009 0000", rep, sizeof rep);              /* UART replies OFF */
    ser_close();
    return t >= 0 ? (float)(int16_t)t : NAN;
}

static int cmd_warmup(int argc, char **argv)
{
    const char *out = NULL, *dev = "/dev/video0", *port = "/dev/ttyACM0";
    int interval = 30, minutes = 60, n = 64;
    for (int i = 2; i + 1 < argc; i += 2) {
        if      (!strcmp(argv[i], "-o")) out = argv[i + 1];
        else if (!strcmp(argv[i], "-i")) interval = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-m")) minutes = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-n")) n = atoi(argv[i + 1]);
        else if (!strcmp(argv[i], "-d")) dev = argv[i + 1];
        else if (!strcmp(argv[i], "-p")) port = argv[i + 1];
        else { out = NULL; break; }
    }
    if (!out || interval < 5 || minutes < 1 || n < 16) {
        fprintf(stderr, "usage: nuc2 warmup -o warm [-i 30] [-m 60] [-n 64] [-d /dev/video0] [-p /dev/ttyACM0|none]\n");
        return 1;
    }
    struct maphdr h = { .magic = "TVEMAP1", .w = W, .h = H, .target_c = NAN };
    read_settings(port, &h.s);                               /* sensor settings once */
    print_settings("settings:", &h.s);
    if (video_open(dev) < 0) return 1;

    char fn[512];
    snprintf(fn, sizeof fn, "%s.csv", out);
    FILE *csv = fopen(fn, "w");
    if (!csv) { perror(fn); return 1; }
    fprintf(csv, "index,time_s,t70_c,mean_dl,tl,tr,bl,br,noise_dl,reconnects\n");

    uint8_t *buf = malloc(FRAME_BYTES);
    uint16_t *img = malloc(NPIX * sizeof *img);
    double *acc = malloc(NPIX * sizeof *acc), *acc2 = malloc(NPIX * sizeof *acc2);
    float *avg = malloc(NPIX * sizeof *avg);
    time_t t0 = time(NULL);
    int total = minutes * 60 / interval + 1;
    printf("recording %d maps, one every %d s - keep the uniform target in front of the lens\n",
           total, interval);
    printf("  idx   time   T70   mean      TL     TR     BL     BR\n");

    for (int k = 0; k < total; k++) {
        /* wait for the slot; keep reading frames so the stream stays fresh */
        while (time(NULL) < t0 + (time_t)k * interval)
            grab_keep(dev, NULL);
        memset(acc, 0, NPIX * sizeof *acc); memset(acc2, 0, NPIX * sizeof *acc2);
        for (int got = 0; got < n; ) {
            int r = grab_keep(dev, buf);
            if (r <= 0) continue;
            decode_yuy2(buf, img);
            for (int i = 0; i < NPIX; i++) { acc[i] += img[i]; acc2[i] += (double)img[i] * img[i]; }
            got++;
        }
        double tn = 0, mean = 0, q[4] = {0};
        for (int i = 0; i < NPIX; i++) {
            double m = acc[i] / n, var = acc2[i] / n - m * m;
            avg[i] = (float)m; mean += m;
            tn += var > 0 ? sqrt(var) : 0;
            q[((i / W) >= H / 2) * 2 + ((i % W) >= W / 2)] += m;
        }
        mean /= NPIX;
        for (int j = 0; j < 4; j++) q[j] /= NPIX / 4;
        h.nframes = n;
        h.noise_dl = (float)(tn / NPIX);
        h.unix_time = (int64_t)time(NULL);
        h.s.t70_c = read_t70(port);
        snprintf(fn, sizeof fn, "%s_%03d.map", out, k);
        if (save_map(fn, &h, avg) < 0) return 1;
        long el = (long)(time(NULL) - t0);
        fprintf(csv, "%d,%ld,%.0f,%.1f,%.0f,%.0f,%.0f,%.0f,%.2f,%d\n",
                k, el, h.s.t70_c, mean, q[0], q[1], q[2], q[3], h.noise_dl, reconnects);
        fflush(csv);
        printf("  %3d  %4ld s  %4.0f  %7.1f  %6.0f %6.0f %6.0f %6.0f\n",
               k, el, h.s.t70_c, mean, q[0], q[1], q[2], q[3]);
        fflush(stdout);
    }
    fclose(csv);
    printf("done: %s_000.map .. %s_%03d.map and %s.csv\n", out, out, total - 1, out);
    return 0;
}
#endif

/* ======================================================================
 * calc: build the 2-point NUC tables
 * ====================================================================== */
static int same_settings(const struct settings *a, const struct settings *b)
{
    return a->tint == b->tint && a->gfid == b->gfid && a->gsk == b->gsk &&
           a->gain_reg == b->gain_reg;
}

static int cmd_calc(int argc, char **argv)
{
    const char *out = "nuc2.bin", *files[3];
    int nf = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i + 1 < argc) out = argv[++i];
        else if (nf < 3) files[nf++] = argv[i];
    }
    if (nf < 2) {
        fprintf(stderr, "usage: nuc2 calc cold1.map hot.map [cold2.map] -o nuc2.bin\n");
        return 1;
    }

    struct maphdr hc, hh, hc2;
    float *I1 = malloc(NPIX * 4), *I2 = malloc(NPIX * 4), *C2 = malloc(NPIX * 4);
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4), *R = malloc(NPIX * 4);
    uint8_t *bad = calloc(NPIX, 1);
    if (!I1 || !I2 || !C2 || !G || !O || !R || !bad) { fprintf(stderr, "out of memory\n"); return 1; }
    if (load_map(files[0], &hc, I1) < 0 || load_map(files[1], &hh, I2) < 0) return 1;

    print_settings("cold :", &hc.s);
    print_settings("hot  :", &hh.s);
    if (!same_settings(&hc.s, &hh.s)) {
        fprintf(stderr, "ERROR: cold / hot sensor settings differ (TINT/GFID/GSK/GAIN). "
                        "Gain changes with them (NUC slide 19).\n");
        return 1;
    }

    /* cold - hot - cold: average both cold maps to cancel linear drift */
    if (nf == 3) {
        if (load_map(files[2], &hc2, C2) < 0) return 1;
        if (!same_settings(&hc.s, &hc2.s)) { fprintf(stderr, "ERROR: cold2 settings differ\n"); return 1; }
        double d1, s1, d2, s2;
        stats(I1, NULL, &d1, &s1);
        stats(C2, NULL, &d2, &s2);
        printf("cold1 vs cold2 mean difference: %+.1f DL (drift during the captures)\n", d2 - d1);
        /* pattern part of the drift: spatial sd of (cold2 - cold1), mean removed */
        float *dd = malloc(NPIX * 4);
        for (int i = 0; i < NPIX; i++) dd[i] = C2[i] - I1[i];
        double dm, dsd;
        stats(dd, NULL, &dm, &dsd);
        printf("cold1 vs cold2 pattern drift : %.2f DL sd (non-uniform part; smaller is better)\n", dsd);
        free(dd);
        for (int i = 0; i < NPIX; i++) I1[i] = 0.5f * (I1[i] + C2[i]);
    }

    double t1 = hc.target_c, t2 = hh.target_c;
    printf("T1 = %.1f C, T2 = %.1f C\n", t1, t2);
    if (!(t2 > t1)) { fprintf(stderr, "ERROR: hot target must be warmer than cold target\n"); return 1; }

    /* 1. response per pixel, dead / hot pixels */
    for (int i = 0; i < NPIX; i++) {
        R[i] = I2[i] - I1[i];
        bad[i] = I1[i] < DEAD_DL || I1[i] > HOT_DL || I2[i] < DEAD_DL || I2[i] > HOT_DL;
    }
    double rmed = pct(R, bad, 50);
    printf("response median %.1f DL  (= %.2f DL per deg C)\n", rmed, rmed / (t2 - t1));
    if (rmed < 20) printf("WARNING: very small response - increase T2 - T1 or check the targets\n");
    for (int i = 0; i < NPIX; i++)
        if (R[i] < RESP_LO * rmed || R[i] > RESP_HI * rmed) bad[i] = 1;

    /* 2. M1, M2 over good pixels (slide 16) */
    double m1, s1, m2, s2;
    stats(I1, bad, &m1, &s1);
    stats(I2, bad, &m2, &s2);

    /* 3. gain and offset per pixel (slide 17) */
    for (int i = 0; i < NPIX; i++) {
        if (bad[i]) { G[i] = 1; O[i] = 0; continue; }
        G[i] = (float)((m2 - m1) / R[i]);
        O[i] = (float)(m1 - G[i] * I1[i]);
    }
    double gm, gsd;
    stats(G, bad, &gm, &gsd);
    for (int i = 0; i < NPIX; i++)
        if (!bad[i] && fabs(G[i] - gm) > G_SIGMA * gsd) { bad[i] = 1; G[i] = 1; O[i] = 0; }
    stats(G, bad, &gm, &gsd);

    int nbad = 0;
    for (int i = 0; i < NPIX; i++) nbad += bad[i];

    /* response uniformity per quadrant: did the target fill the whole view? */
    double q[4] = {0}; long qn[4] = {0};
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            int i = y * W + x, k = (y >= H / 2) * 2 + (x >= W / 2);
            if (!bad[i]) { q[k] += R[i]; qn[k]++; }
        }
    double qmin = 1e9, qmax = -1e9;
    for (int k = 0; k < 4; k++) { q[k] /= qn[k]; qmin = fmin(qmin, q[k]); qmax = fmax(qmax, q[k]); }

    printf("\nM1 %.1f  M2 %.1f  (M2-M1 = %.1f DL)\n", m1, m2, m2 - m1);
    printf("gain G: mean %.4f  sd %.4f  p1 %.3f  p99 %.3f\n", gm, gsd, pct(G, bad, 1), pct(G, bad, 99));
    printf("bad pixels: %d\n", nbad);
    printf("response quadrants TL %.0f TR %.0f BL %.0f BR %.0f  (spread %.1f%%)\n",
           q[0], q[1], q[2], q[3], 100 * (qmax - qmin) / rmed);
    if ((qmax - qmin) / rmed > 0.10)
        printf("WARNING: response differs > 10%% between quadrants - target not filling the view "
               "or not uniform? Look at nuc2_resp.pgm\n");

    /* save tables */
    struct nuchdr nh = { .magic = "TVENUC2", .w = W, .h = H, .m1 = (float)m1, .m2 = (float)m2,
                         .t1_c = (float)t1, .t2_c = (float)t2, .nbad = nbad, .s = hc.s };
    FILE *f = fopen(out, "wb");
    if (!f) { perror(out); return 1; }
    fwrite(&nh, sizeof nh, 1, f);
    fwrite(G, 4, NPIX, f);
    fwrite(O, 4, NPIX, f);
    fwrite(bad, 1, NPIX, f);
    fclose(f);

    write_pgm("nuc2_gain.pgm", G, bad, gm - 4 * gsd, gm + 4 * gsd);
    write_pgm("nuc2_resp.pgm", R, NULL, pct(R, bad, 1), pct(R, bad, 99));
    float *bm = malloc(NPIX * 4);
    for (int i = 0; i < NPIX; i++) bm[i] = bad[i] ? 255.0f : 0.0f;
    write_pgm("nuc2_bad.pgm", bm, NULL, 0, 255);
    free(bm);
    printf("\nsaved %s  (+ nuc2_gain.pgm, nuc2_resp.pgm, nuc2_bad.pgm)\n", out);
    return 0;
}

/* ======================================================================
 * check: residual FPN on a third map (slide 20)
 * ====================================================================== */
static int cmd_check(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: nuc2 check nuc2.bin test.map [ref.map]\n"); return 1; }
    struct nuchdr nh;
    struct maphdr th, ch;
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4), *T = malloc(NPIX * 4);
    float *C = malloc(NPIX * 4), *tmp = malloc(NPIX * 4);
    uint8_t *bad = malloc(NPIX);

    FILE *f = fopen(argv[2], "rb");
    if (!f) { perror(argv[2]); return 1; }
    if (fread(&nh, sizeof nh, 1, f) != 1 || memcmp(nh.magic, "TVENUC2", 7) ||
        fread(G, 4, NPIX, f) != NPIX || fread(O, 4, NPIX, f) != NPIX ||
        fread(bad, 1, NPIX, f) != NPIX) { fprintf(stderr, "%s: not a nuc2 file\n", argv[2]); return 1; }
    fclose(f);
    if (load_map(argv[3], &th, T) < 0) return 1;
    if (!same_settings(&nh.s, &th.s))
        printf("WARNING: test map settings differ from the calibration\n");

    printf("calibration T1 %.1f / T2 %.1f C (board %.0f C), test target %.1f C (board %.0f C)\n",
           nh.t1_c, nh.t2_c, nh.s.t70_c, th.target_c, th.s.t70_c);

    double m, raw_sd, one_sd = NAN, two_sd, ref_sd = NAN;
    stats(T, bad, &m, &raw_sd);

    /* 2-point with the offset from calibration time (slide 17) */
    for (int i = 0; i < NPIX; i++) tmp[i] = T[i] * G[i] + O[i];
    stats(tmp, bad, &m, &two_sd);

    /* With a uniform reference map (lens cap / blackbody) */
    if (argc >= 5 && load_map(argv[4], &ch, C) == 0) {
        double cm, cs;
        stats(C, bad, &cm, &cs);
        /* 1-point (offset only, slide 8): T - (C - mean C) */
        for (int i = 0; i < NPIX; i++) tmp[i] = T[i] - (C[i] - (float)cm);
        stats(tmp, bad, &m, &one_sd);
        /* 2-point gain + offset refreshed from the reference (slide 18):
           O' = mean(G*C) - G*C. If the reference is taken just before the test,
           offset drift cancels and only the gain quality remains. */
        double gcm = 0; long n = 0;
        for (int i = 0; i < NPIX; i++) if (!bad[i]) { gcm += G[i] * C[i]; n++; }
        gcm /= n;
        for (int i = 0; i < NPIX; i++) tmp[i] = T[i] * G[i] + (float)(gcm - G[i] * C[i]);
        stats(tmp, bad, &m, &ref_sd);
        printf("reference %s: target %.1f C (board %.0f C)\n", argv[4], ch.target_c, ch.s.t70_c);
    }

    printf("\nResidual FPN (spatial sd over good pixels):\n");
    printf("  raw                       : %8.2f DL\n", raw_sd);
    if (!isnan(one_sd))
        printf("  1-point  (offset from ref): %8.2f DL\n", one_sd);
    printf("  2-point  (calib offset)   : %8.2f DL   (includes offset drift since calibration)\n", two_sd);
    if (!isnan(ref_sd))
        printf("  2-point  (offset from ref): %8.2f DL   (gain quality if ref is fresh)\n", ref_sd);
    printf("  noise floor of this map   : %.2f DL/frame / sqrt(%u) = %.2f DL\n",
           th.noise_dl, th.nframes, th.noise_dl / sqrt((double)th.nframes));
    return 0;
}

/* ======================================================================
 * apply: correct one frame of a recorded scene and write comparison images
 * ====================================================================== */
static int cmp_u16(const void *a, const void *b)
{
    float x = *(const float *)a, y = *(const float *)b;
    return (x > y) - (x < y);
}

/* Replace bad pixels by the median of their good 8-neighbours */
static void fix_bad(float *v, const uint8_t *bad)
{
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            int i = y * W + x;
            if (!bad[i]) continue;
            float nb[8]; int n = 0;
            for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++) {
                    int xx = x + dx, yy = y + dy;
                    if ((!dx && !dy) || xx < 0 || yy < 0 || xx >= W || yy >= H) continue;
                    if (!bad[yy * W + xx]) nb[n++] = v[yy * W + xx];
                }
            if (n) { qsort(nb, n, sizeof nb[0], cmp_u16); v[i] = nb[n / 2]; }
        }
}

/* Copy one image into a quadrant of the 2x2 montage, stretched 1%..99% to 0..255 */
static void put_panel(uint8_t *mont, int qx, int qy, const float *v, const uint8_t *bad)
{
    double lo = pct(v, bad, 1), hi = pct(v, bad, 99);
    if (hi <= lo) hi = lo + 1;
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            double t = (v[y * W + x] - lo) / (hi - lo) * 255.0;
            mont[(qy * H + y) * (2 * W) + qx * W + x] = t < 0 ? 0 : t > 255 ? 255 : (uint8_t)t;
        }
}

static int cmd_apply(int argc, char **argv)
{
    if (argc < 4) goto usage;
    const char *nucf = argv[2], *raw = argv[3], *reff = NULL;
    int frame = 30;
    for (int i = 4; i < argc; i++) {
        if (!strcmp(argv[i], "-f") && i + 1 < argc) frame = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-r") && i + 1 < argc) reff = argv[++i];
        else goto usage;
    }

    /* load gain / offset / bad pixel tables */
    struct nuchdr nh;
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4);
    uint8_t *bad = malloc(NPIX);
    FILE *f = fopen(nucf, "rb");
    if (!f) { perror(nucf); return 1; }
    if (fread(&nh, sizeof nh, 1, f) != 1 || memcmp(nh.magic, "TVENUC2", 7) ||
        fread(G, 4, NPIX, f) != NPIX || fread(O, 4, NPIX, f) != NPIX ||
        fread(bad, 1, NPIX, f) != NPIX) { fprintf(stderr, "%s: not a nuc2 file\n", nucf); return 1; }
    fclose(f);

    /* read one frame of the scene */
    f = fopen(raw, "rb");
    if (!f) { perror(raw); return 1; }
    long long nfile = file_size(f) / FRAME_BYTES;
    if (frame < 0 || frame >= nfile) {
        fprintf(stderr, "frame %d out of range (file has %lld frames)\n", frame, nfile);
        return 1;
    }
    uint8_t *buf = malloc(FRAME_BYTES);
    uint16_t *img = malloc(NPIX * sizeof *img);
    file_seek(f, (long long)frame * FRAME_BYTES);
    if (fread(buf, 1, FRAME_BYTES, f) != FRAME_BYTES) { fprintf(stderr, "read error\n"); return 1; }
    fclose(f);
    decode_yuy2(buf, img);

    float *rawv = malloc(NPIX * 4), *one = malloc(NPIX * 4);
    float *two = malloc(NPIX * 4), *tworef = malloc(NPIX * 4), *C = malloc(NPIX * 4);
    for (int i = 0; i < NPIX; i++) rawv[i] = img[i];

    /* panel 2: 2-point with the offset from calibration time (slide 17) */
    for (int i = 0; i < NPIX; i++) two[i] = rawv[i] * G[i] + O[i];
    fix_bad(two, bad);

    /* panels 3 and 4 need a fresh uniform reference (lens cap / plate) taken now */
    int have_ref = 0;
    if (reff) {
        struct maphdr rh;
        if (load_map(reff, &rh, C) < 0) return 1;
        have_ref = 1;
        /* 1-point (offset only): corr = raw - (C - mean C)   (slide 8) */
        double cm, cs;
        stats(C, bad, &cm, &cs);
        for (int i = 0; i < NPIX; i++) one[i] = rawv[i] - (C[i] - (float)cm);
        fix_bad(one, bad);
        /* 2-point with refreshed offset: O' = mean(G*C) - G*C   (slide 18) */
        double gcm = 0; long n = 0;
        for (int i = 0; i < NPIX; i++) if (!bad[i]) { gcm += G[i] * C[i]; n++; }
        gcm /= n;
        for (int i = 0; i < NPIX; i++) tworef[i] = rawv[i] * G[i] + (float)(gcm - G[i] * C[i]);
        fix_bad(tworef, bad);
    }

    /* 2x2 montage: raw | 2-point (calib offset) / 1-point | 2-point (refreshed offset) */
    uint8_t *mont = calloc((size_t)4 * NPIX, 1);
    put_panel(mont, 0, 0, rawv, bad);
    put_panel(mont, 1, 0, two, bad);
    if (have_ref) { put_panel(mont, 0, 1, one, bad); put_panel(mont, 1, 1, tworef, bad); }
    FILE *o = fopen("apply_all.pgm", "wb");
    if (!o) { perror("apply_all.pgm"); return 1; }
    fprintf(o, "P5\n%d %d\n255\n", 2 * W, 2 * H);
    fwrite(mont, 1, (size_t)4 * NPIX, o);
    fclose(o);

    printf("frame %d of %s corrected with %s (calibrated T1 %.1f / T2 %.1f C)\n",
           frame, raw, nucf, nh.t1_c, nh.t2_c);
    printf("saved apply_all.pgm (1280x960):\n");
    printf("  top-left    : raw (no correction)\n");
    printf("  top-right   : 2-point, offset from calibration time\n");
    if (have_ref) {
        printf("  bottom-left : 1-point, offset from %s (like thermal_live)\n", reff);
        printf("  bottom-right: 2-point, gain from calibration + offset refreshed from %s\n", reff);
    } else {
        printf("  bottom row  : empty (give -r ref.map taken just before the scene)\n");
    }
    return 0;

usage:
    fprintf(stderr, "usage: nuc2 apply nuc2.bin scene.raw [-f frame] [-r ref.map]\n");
    return 1;
}

/* ======================================================================
 * noise: temporal noise split into pixel / column / row / whole-frame parts
 *   nuc2 noise file.raw [-s 30] [-n 128]
 * Point the camera at a uniform, stable target for the recording.
 * Each frame minus the average frame leaves only the time-varying part.
 *   pixel  : sd over time per pixel, averaged over pixels
 *   column : per frame, mean of each column (minus the frame mean); sd over
 *            columns. Pure pixel noise gives pixel/sqrt(480). Much bigger
 *            means the columns flicker together (temporal column noise),
 *            which an offset table (NUC) cannot remove.
 *   row    : same for rows (expected pixel/sqrt(640))
 *   frame  : sd over time of the whole-frame mean (global flicker)
 * ====================================================================== */
static int cmd_noise(int argc, char **argv)
{
    if (argc < 3) { fprintf(stderr, "usage: nuc2 noise file.raw [-s 30] [-n 128]\n"); return 1; }
    /* out/target are not used here; set so the common parser accepts it */
    struct opts o = { .out = "-", .n = 128, .skip = 30, .target = 0, .port = "none" };
    if (parse_opts(argc, argv, 3, &o) < 0) { fprintf(stderr, "bad options (-n must be >= 16)\n"); return 1; }
    FILE *f = fopen(argv[2], "rb");
    if (!f) { perror(argv[2]); return 1; }
    long long nfile = file_size(f) / FRAME_BYTES;
    if (nfile < o.skip + o.n) { fprintf(stderr, "file has %lld frames, need %d\n", nfile, o.skip + o.n); return 1; }

    uint8_t *buf = malloc(FRAME_BYTES);
    uint16_t *img = malloc(NPIX * sizeof *img);
    double *mean = calloc(NPIX, sizeof *mean), *var = calloc(NPIX, sizeof *var);
    double *fm = malloc(o.n * sizeof *fm);
    double col[W], row[H];

    /* pass 1: average frame */
    file_seek(f, (long long)o.skip * FRAME_BYTES);
    for (int k = 0; k < o.n; k++) {
        if (fread(buf, 1, FRAME_BYTES, f) != FRAME_BYTES) { fprintf(stderr, "read error\n"); return 1; }
        decode_yuy2(buf, img);
        for (int i = 0; i < NPIX; i++) mean[i] += img[i];
    }
    for (int i = 0; i < NPIX; i++) mean[i] /= o.n;

    /* pass 2: residuals */
    double colsd = 0, rowsd = 0;
    file_seek(f, (long long)o.skip * FRAME_BYTES);
    for (int k = 0; k < o.n; k++) {
        if (fread(buf, 1, FRAME_BYTES, f) != FRAME_BYTES) { fprintf(stderr, "read error\n"); return 1; }
        decode_yuy2(buf, img);
        memset(col, 0, sizeof col); memset(row, 0, sizeof row);
        double tot = 0;
        for (int y = 0; y < H; y++)
            for (int x = 0; x < W; x++) {
                double d = img[y * W + x] - mean[y * W + x];
                var[y * W + x] += d * d;
                col[x] += d; row[y] += d; tot += d;
            }
        tot /= NPIX;
        fm[k] = tot;
        double s2 = 0;
        for (int x = 0; x < W; x++) { double c = col[x] / H - tot; s2 += c * c; }
        colsd += sqrt(s2 / W);
        s2 = 0;
        for (int y = 0; y < H; y++) { double r = row[y] / W - tot; s2 += r * r; }
        rowsd += sqrt(s2 / H);
    }
    fclose(f);
    colsd /= o.n; rowsd /= o.n;
    double pix = 0;
    for (int i = 0; i < NPIX; i++) pix += sqrt(var[i] / (o.n - 1));
    pix /= NPIX;
    double fsd = 0;
    for (int k = 0; k < o.n; k++) fsd += fm[k] * fm[k];
    fsd = sqrt(fsd / o.n);

    printf("%s: %d frames (skip %d)\n", argv[2], o.n, o.skip);
    printf("  pixel temporal noise : %6.3f DL\n", pix);
    printf("  column flicker       : %6.3f DL   (pure pixel noise would give %.3f)  ratio %.1f\n",
           colsd, pix / sqrt((double)H), colsd / (pix / sqrt((double)H)));
    printf("  row flicker          : %6.3f DL   (pure pixel noise would give %.3f)  ratio %.1f\n",
           rowsd, pix / sqrt((double)W), rowsd / (pix / sqrt((double)W)));
    printf("  whole-frame flicker  : %6.3f DL\n", fsd);
    printf("A ratio near 1 = no extra flicker. A ratio of several = lines that change\n"
           "every frame; a fixed NUC table cannot remove those, only a per-frame filter.\n");
    return 0;
}

/* ======================================================================
 * stripes: what is left after the 2-point correction on a uniform target
 *   nuc2 stripes nuc2.bin wall.map
 * Splits the residual into
 *   shading : slow brightness change over the image (corners vs centre)
 *   columns : column-to-column steps (vertical lines)
 *   rows    : row-to-row steps (horizontal lines)
 * and checks whether a column line is the same from top to bottom
 * (4 bands of 120 rows), so we know which filter can remove it.
 * Writes stripes_view.pgm: the residual, mid grey = 0, +-20 DL = white/black.
 * ====================================================================== */
static void smooth1d(const double *in, double *out, int n, int r)
{
    for (int i = 0; i < n; i++) {
        double s = 0; int c = 0;
        for (int k = i - r; k <= i + r; k++) if (k >= 0 && k < n) { s += in[k]; c++; }
        out[i] = s / c;
    }
}

static int cmd_stripes(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: nuc2 stripes nuc2.bin wall.map\n"); return 1; }
    struct nuchdr nh; struct maphdr mh;
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4), *v = malloc(NPIX * 4);
    uint8_t *bad = malloc(NPIX);
    FILE *f = fopen(argv[2], "rb");
    if (!f) { perror(argv[2]); return 1; }
    if (fread(&nh, sizeof nh, 1, f) != 1 || memcmp(nh.magic, "TVENUC2", 7) ||
        fread(G, 4, NPIX, f) != NPIX || fread(O, 4, NPIX, f) != NPIX ||
        fread(bad, 1, NPIX, f) != NPIX) { fprintf(stderr, "%s: not a nuc2 file\n", argv[2]); return 1; }
    fclose(f);
    if (load_map(argv[3], &mh, v) < 0) return 1;

    /* 2-point with the calibration offset, then remove the mean */
    double m, sd;
    for (int i = 0; i < NPIX; i++) v[i] = v[i] * G[i] + O[i];
    stats(v, bad, &m, &sd);
    for (int i = 0; i < NPIX; i++) v[i] -= (float)m;

    /* column profile per band of 120 rows and for the whole image */
    enum { NB = 4 };
    static double colb[NB][W], col[W], row[H], tmp[W > H ? W : H];
    static int ncolb[NB][W], ncol[W], nrow[H];
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            int i = y * W + x;
            if (bad[i]) continue;
            colb[y / (H / NB)][x] += v[i]; ncolb[y / (H / NB)][x]++;
            col[x] += v[i]; ncol[x]++;
            row[y] += v[i]; nrow[y]++;
        }
    for (int x = 0; x < W; x++) {
        col[x] /= ncol[x] ? ncol[x] : 1;
        for (int b = 0; b < NB; b++) colb[b][x] /= ncolb[b][x] ? ncolb[b][x] : 1;
    }
    for (int y = 0; y < H; y++) row[y] /= nrow[y] ? nrow[y] : 1;

    /* high-pass = profile minus its +-8 smoothed version (removes shading) */
    double hpc[W], hpr[H], hpb[NB][W];
    smooth1d(col, tmp, W, 8);
    for (int x = 0; x < W; x++) hpc[x] = col[x] - tmp[x];
    smooth1d(row, tmp, H, 8);
    for (int y = 0; y < H; y++) hpr[y] = row[y] - tmp[y];
    for (int b = 0; b < NB; b++) {
        smooth1d(colb[b], tmp, W, 8);
        for (int x = 0; x < W; x++) hpb[b][x] = colb[b][x] - tmp[x];
    }
    /* stats skip 8 edge columns/rows; edges reported separately */
    double sc = 0, sr = 0, sdiff = 0; int nc = 0, nr = 0, nd = 0;
    for (int x = 8; x < W - 8; x++) {
        sc += hpc[x] * hpc[x]; nc++;
        for (int b = 0; b < NB; b++) { double d = hpb[b][x] - hpc[x]; sdiff += d * d; nd++; }
    }
    for (int y = 8; y < H - 8; y++) { sr += hpr[y] * hpr[y]; nr++; }
    sc = sqrt(sc / nc); sr = sqrt(sr / nr); sdiff = sqrt(sdiff / nd);

    /* shading: mean of 64x64 blocks at corners and centre */
    double blk[5]; const int bx[5] = { 0, W - 64, 0, W - 64, W / 2 - 32 }, by[5] = { 0, 0, H - 64, H - 64, H / 2 - 32 };
    for (int k = 0; k < 5; k++) {
        double s2 = 0; int n = 0;
        for (int y = by[k]; y < by[k] + 64; y++)
            for (int x = bx[k]; x < bx[k] + 64; x++) if (!bad[y * W + x]) { s2 += v[y * W + x]; n++; }
        blk[k] = s2 / n;
    }

    printf("2-point residual on %s (target %.1f C) with %s\n", argv[3], mh.target_c, argv[2]);
    printf("  total residual sd           : %6.2f DL  (temporal noise of the map %.2f DL)\n",
           sd, mh.noise_dl / sqrt((double)mh.nframes));
    printf("  shading (64x64 block means) : TL %+.1f  TR %+.1f  BL %+.1f  BR %+.1f  centre %+.1f DL\n",
           blk[0], blk[1], blk[2], blk[3], blk[4]);
    printf("  vertical lines   (column sd): %6.2f DL\n", sc);
    printf("    same top->bottom?         : band-to-band difference %.2f DL  (%s)\n", sdiff,
           sdiff < 0.5 * sc ? "yes - one value per column removes them" :
                              "no - lines change along the column");
    printf("  horizontal lines (row sd)   : %6.2f DL\n", sr);
    printf("  left  edge columns 0..7     :");
    for (int x = 0; x < 8; x++) printf(" %+.1f", col[x]);
    printf("\n  right edge columns 632..639 :");
    for (int x = W - 8; x < W; x++) printf(" %+.1f", col[x]);
    printf("\n  strongest columns           :");
    for (int k = 0; k < 8; k++) {
        int best = -1; double bv = 0;
        for (int x = 8; x < W - 8; x++) if (fabs(hpc[x]) > bv) { bv = fabs(hpc[x]); best = x; }
        if (best < 0) break;
        printf(" x%d:%+.1f", best, hpc[best]);
        hpc[best] = 0;
    }
    printf("\n");
    write_pgm("stripes_view.pgm", v, bad, -20, 20);
    printf("saved stripes_view.pgm (grey = 0, white = +20 DL, black = -20 DL)\n");
    return 0;
}


/* ======================================================================
 * t70: show the USB3-A3 board temperature (register R70) every 2 s.
 *   nuc3 t70 -p COM5            (Windows)
 *   nuc3 t70 -p /dev/ttyACM0    (Linux)        Ctrl+C to stop
 * Read-only: UART replies ON, read R70, replies OFF (FullScale UG004).
 * Can run while the video is streaming (serial port only, no video).
 * ====================================================================== */
static int cmd_t70(int argc, char **argv)
{
    const char *port = NULL;
    for (int i = 2; i + 1 < argc; i++) if (!strcmp(argv[i], "-p")) port = argv[i + 1];
    if (!port) { fprintf(stderr, "usage: nuc3 t70 -p COM5|/dev/ttyACM0\n"); return 1; }
    if (ser_open(port) < 0) return 1;
    char rep[64];
    uart_cmd("W01 0009 0001", rep, sizeof rep);              /* UART replies ON  */
    printf("board temperature T70 every 2 s (Ctrl+C to stop)\n");
    for (;;) {
        long t = uart_read("R70 0000 0000");
        time_t now = time(NULL);
        char ts[16];
        strftime(ts, sizeof ts, "%H:%M:%S", localtime(&now));
        if (t >= 0) printf("%s  T70 = %d C\n", ts, (int)(int16_t)t);
        else        printf("%s  no reply\n", ts);
        fflush(stdout);
#ifdef _WIN32
        Sleep(2000);
#else
        sleep(2);
#endif
    }
    return 0;
}

/* ======================================================================
 * omodel: offset vs temperature model from a warm-up recording
 *   nuc2 omodel nuc2.bin warm [-o omodel.bin]
 * Reads warm_000.map, warm_001.map ... (uniform target, made by "warmup").
 * For every map k: c = G*map (gain from the 2-point NUC), pattern
 * p = c - mean(c). The pattern is what the offset must cancel.
 * Per pixel a quadratic in board temperature T (T70) is fitted:
 *     p_i(T) = a_i + b_i*(T-T0) + c_i*(T-T0)^2
 * Runtime correction (live5): corrected = G*raw - p_i(T).
 * The design matrix is the same for every pixel, so it is inverted once.
 * Evaluation per map: residual with the old fixed 2-point offset vs with the
 * model at the recorded T70, and with the best-fitting T ("self T"), which
 * shows how much is lost because T70 is coarse (1 C steps, board not FPA).
 * ====================================================================== */
struct omhdr {
    char     magic[8];                 /* "TVEOMD1"                          */
    uint32_t w, h, nmaps;
    float    t0, tmin, tmax;
    /* followed by float a[NPIX], b[NPIX], c[NPIX] */
};

static int cmd_omodel(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: nuc2 omodel nuc2.bin warm [-o omodel.bin]\n"); return 1; }
    const char *out = "omodel.bin", *pre = argv[3];
    if (argc >= 6 && !strcmp(argv[4], "-o")) out = argv[5];

    struct nuchdr nh;
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4), *v = malloc(NPIX * 4);
    uint8_t *bad = malloc(NPIX);
    FILE *f = fopen(argv[2], "rb");
    if (!f) { perror(argv[2]); return 1; }
    if (fread(&nh, sizeof nh, 1, f) != 1 || memcmp(nh.magic, "TVENUC2", 7) ||
        fread(G, 4, NPIX, f) != NPIX || fread(O, 4, NPIX, f) != NPIX ||
        fread(bad, 1, NPIX, f) != NPIX) { fprintf(stderr, "%s: not a nuc2 file\n", argv[2]); return 1; }
    fclose(f);

    /* list maps and their T70 */
    enum { MAXM = 2000 };
    static float T[MAXM];
    int nm = 0;
    char fn[512];
    struct maphdr mh;
    double tsum = 0, tmin = 1e9, tmax = -1e9;
    int nvalid = 0;
    for (nm = 0; nm < MAXM; nm++) {
        snprintf(fn, sizeof fn, "%s_%03d.map", pre, nm);
        FILE *t = fopen(fn, "rb");
        if (!t) break;
        if (fread(&mh, sizeof mh, 1, t) != 1) { fclose(t); break; }
        fclose(t);
        T[nm] = mh.s.t70_c;               /* NAN = T70 not read: skipped below */
        if (isnan(T[nm])) continue;
        tsum += T[nm]; nvalid++;
        if (T[nm] < tmin) tmin = T[nm];
        if (T[nm] > tmax) tmax = T[nm];
    }
    if (nvalid < 5 || tmax - tmin < 3) { fprintf(stderr, "need >= 5 maps over >= 3 C (found %d maps, %.0f..%.0f C)\n", nm, tmin, tmax); return 1; }
    double t0 = tsum / nvalid;
    printf("%d maps, T70 %.0f .. %.0f C, T0 = %.2f C\n", nm, tmin, tmax, t0);

    /* normal equations: same 3x3 matrix for every pixel */
    double M[3][3] = {{0}}, *S0 = calloc(NPIX, 8), *S1 = calloc(NPIX, 8), *S2 = calloc(NPIX, 8);
    for (int k = 0; k < nm; k++) {
        snprintf(fn, sizeof fn, "%s_%03d.map", pre, k);
        if (load_map(fn, &mh, v) < 0) return 1;
        if (isnan(T[k])) continue;
        double x = T[k] - t0, xp[3] = { 1, x, x * x }, m = 0; long n = 0;
        for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) M[r][c] += xp[r] * xp[c];
        for (int i = 0; i < NPIX; i++) { v[i] *= G[i]; if (!bad[i]) { m += v[i]; n++; } }
        m /= n;
        for (int i = 0; i < NPIX; i++) {
            double pv = v[i] - m;
            S0[i] += pv; S1[i] += pv * x; S2[i] += pv * x * x;
        }
    }
    /* invert M (3x3, symmetric) */
    double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
               - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
               + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
    if (fabs(det) < 1e-9) { fprintf(stderr, "temperature values too similar for a quadratic fit\n"); return 1; }
    double Mi[3][3];
    Mi[0][0] =  (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / det;
    Mi[0][1] = -(M[0][1] * M[2][2] - M[0][2] * M[2][1]) / det;
    Mi[0][2] =  (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
    Mi[1][0] = -(M[1][0] * M[2][2] - M[1][2] * M[2][0]) / det;
    Mi[1][1] =  (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det;
    Mi[1][2] = -(M[0][0] * M[1][2] - M[0][2] * M[1][0]) / det;
    Mi[2][0] =  (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / det;
    Mi[2][1] = -(M[0][0] * M[2][1] - M[0][1] * M[2][0]) / det;
    Mi[2][2] =  (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;
    float *A = malloc(NPIX * 4), *B = malloc(NPIX * 4), *C = malloc(NPIX * 4);
    for (int i = 0; i < NPIX; i++) {
        A[i] = (float)(Mi[0][0] * S0[i] + Mi[0][1] * S1[i] + Mi[0][2] * S2[i]);
        B[i] = (float)(Mi[1][0] * S0[i] + Mi[1][1] * S1[i] + Mi[1][2] * S2[i]);
        C[i] = (float)(Mi[2][0] * S0[i] + Mi[2][1] * S1[i] + Mi[2][2] * S2[i]);
    }

    struct omhdr oh = { .magic = "TVEOMD1", .w = W, .h = H, .nmaps = (uint32_t)nm,
                        .t0 = (float)t0, .tmin = (float)tmin, .tmax = (float)tmax };
    f = fopen(out, "wb");
    if (!f || fwrite(&oh, sizeof oh, 1, f) != 1 || fwrite(A, 4, NPIX, f) != NPIX ||
        fwrite(B, 4, NPIX, f) != NPIX || fwrite(C, 4, NPIX, f) != NPIX) { perror(out); return 1; }
    fclose(f);

    /* evaluation */
    snprintf(fn, sizeof fn, "%s_eval.csv", pre);
    FILE *ev = fopen(fn, "w");
    if (ev) fprintf(ev, "index,t70_c,fixed_offset_sd,model_t70_sd,self_t_c,model_self_sd\n");
    printf("\n  idx  T70  | residual sd (DL): fixed 2-pt offset   model@T70   | self T   model@selfT\n");
    double sum_fix = 0, sum_mod = 0, sum_self = 0;
    for (int k = 0; k < nm; k++) {
        if (isnan(T[k])) continue;
        snprintf(fn, sizeof fn, "%s_%03d.map", pre, k);
        if (load_map(fn, &mh, v) < 0) return 1;
        double m = 0; long n = 0;
        for (int i = 0; i < NPIX; i++) { v[i] *= G[i]; if (!bad[i]) { m += v[i]; n++; } }
        m /= n;
        for (int i = 0; i < NPIX; i++) v[i] -= (float)m;           /* pattern p */
        /* fixed offset: p + O (minus its mean) */
        double om = 0; for (int i = 0; i < NPIX; i++) if (!bad[i]) om += O[i]; om /= n;
        double s_fix = 0, s_mod = 0, x = T[k] - t0;
        for (int i = 0; i < NPIX; i++) {
            if (bad[i]) continue;
            double e1 = v[i] + O[i] - om, e2 = v[i] - (A[i] + B[i] * x + C[i] * x * x);
            s_fix += e1 * e1; s_mod += e2 * e2;
        }
        s_fix = sqrt(s_fix / n); s_mod = sqrt(s_mod / n);
        /* best T for this map (grid 0.05 C, every 5th pixel) */
        double best = 1e300, bx = x;
        for (double xs = tmin - t0 - 5; xs <= tmax - t0 + 5; xs += 0.05) {
            double e = 0, em = 0; long ne = 0;
            for (int i = 0; i < NPIX; i += 5) if (!bad[i]) { em += v[i] - (A[i] + B[i] * xs + C[i] * xs * xs); ne++; }
            em /= ne;
            for (int i = 0; i < NPIX; i += 5) if (!bad[i]) {
                double d = v[i] - (A[i] + B[i] * xs + C[i] * xs * xs) - em; e += d * d; }
            if (e < best) { best = e; bx = xs; }
        }
        double s_self = 0;
        for (int i = 0; i < NPIX; i++) if (!bad[i]) {
            double d = v[i] - (A[i] + B[i] * bx + C[i] * bx * bx); s_self += d * d; }
        s_self = sqrt(s_self / n);
        sum_fix += s_fix; sum_mod += s_mod; sum_self += s_self;
        if (ev) fprintf(ev, "%d,%.0f,%.2f,%.2f,%.2f,%.2f\n", k, T[k], s_fix, s_mod, bx + t0, s_self);
        if (k % 10 == 0 || k == nm - 1)
            printf("  %3d  %3.0f  |        %8.1f            %8.1f     |  %5.1f    %8.1f\n",
                   k, T[k], s_fix, s_mod, bx + t0, s_self);
    }
    if (ev) fclose(ev);
    printf("\naverage residual sd: fixed 2-pt offset %.1f DL, model@T70 %.1f DL, model@self T %.1f DL\n",
           sum_fix / nvalid, sum_mod / nvalid, sum_self / nvalid);
    printf("saved %s and %s_eval.csv\n", out, pre);
    return 0;
}

/* ======================================================================
 * omcheck: test the offset model on a new uniform-target map
 *   nuc3 omcheck nuc2.bin omodel.bin test.map
 * test.map must be a uniform target (lens cap) taken with -p so that T70 is
 * in its header. Prints the residual with the fixed 2-point offset, with the
 * model at the recorded T70, and the temperature that would fit best.
 * Writes omcheck_view.pgm (model residual, grey = 0, +-40 DL = white/black).
 * ====================================================================== */
/* sd over good pixels of v - model(x), after removing the mean */
static double res_sd(const float *v, const float *A, const float *B, const float *C,
                     const uint8_t *bad, long n, double x)
{
    double e = 0, em = 0;
    for (int i = 0; i < NPIX; i++)
        if (!bad[i]) em += v[i] - (A[i] + B[i] * x + C[i] * x * x);
    em /= n;
    for (int i = 0; i < NPIX; i++)
        if (!bad[i]) { double d = v[i] - (A[i] + B[i] * x + C[i] * x * x) - em; e += d * d; }
    return sqrt(e / n);
}

static int cmd_omcheck(int argc, char **argv)
{
    if (argc < 5) { fprintf(stderr, "usage: nuc3 omcheck nuc2.bin omodel.bin test.map\n"); return 1; }
    struct nuchdr nh; struct omhdr oh; struct maphdr mh;
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4), *v = malloc(NPIX * 4);
    float *A = malloc(NPIX * 4), *B = malloc(NPIX * 4), *C = malloc(NPIX * 4), *r = malloc(NPIX * 4);
    uint8_t *bad = malloc(NPIX);
    FILE *f = fopen(argv[2], "rb");
    if (!f || fread(&nh, sizeof nh, 1, f) != 1 || memcmp(nh.magic, "TVENUC2", 7) ||
        fread(G, 4, NPIX, f) != NPIX || fread(O, 4, NPIX, f) != NPIX ||
        fread(bad, 1, NPIX, f) != NPIX) { fprintf(stderr, "%s: not a nuc2 file\n", argv[2]); return 1; }
    fclose(f);
    f = fopen(argv[3], "rb");
    if (!f || fread(&oh, sizeof oh, 1, f) != 1 || memcmp(oh.magic, "TVEOMD1", 7) ||
        fread(A, 4, NPIX, f) != NPIX || fread(B, 4, NPIX, f) != NPIX ||
        fread(C, 4, NPIX, f) != NPIX) { fprintf(stderr, "%s: not an omodel file\n", argv[3]); return 1; }
    fclose(f);
    if (load_map(argv[4], &mh, v) < 0) return 1;

    double m = 0, om = 0; long n = 0;
    for (int i = 0; i < NPIX; i++) { v[i] *= G[i]; if (!bad[i]) { m += v[i]; om += O[i]; n++; } }
    m /= n; om /= n;
    for (int i = 0; i < NPIX; i++) v[i] -= (float)m;
    double s_fix = 0;
    for (int i = 0; i < NPIX; i++) if (!bad[i]) { double d = v[i] + O[i] - om; s_fix += d * d; }
    s_fix = sqrt(s_fix / n);
    double s_t70 = NAN, x70 = mh.s.t70_c - oh.t0;
    if (!isnan(mh.s.t70_c)) s_t70 = res_sd(v, A, B, C, bad, n, x70);
    double best = 1e300, bx = 0;
    for (double xs = oh.tmin - oh.t0 - 5; xs <= oh.tmax - oh.t0 + 5; xs += 0.1) {
        double e = res_sd(v, A, B, C, bad, n, xs);
        if (e < best) { best = e; bx = xs; }
    }
    printf("model: T %.0f..%.0f C   test map T70: %.0f C\n", oh.tmin, oh.tmax, mh.s.t70_c);
    printf("  residual sd, fixed 2-point offset : %7.1f DL\n", s_fix);
    printf("  residual sd, model at T70         : %7.1f DL\n", s_t70);
    printf("  best-fitting temperature          : %7.1f C  -> residual %.1f DL\n", bx + oh.t0, best);
    if (best > 0.7 * s_fix)
        printf("  NOTE: no temperature fits well - the pattern now has a shape the model has not seen\n");
    else if (!isnan(s_t70) && s_t70 > 1.5 * best)
        printf("  NOTE: the model fits at %.1f C but T70 says %.0f C - T70 is not tracking the sensor\n",
               bx + oh.t0, mh.s.t70_c);
    double em = 0;
    for (int i = 0; i < NPIX; i++) if (!bad[i]) em += v[i] - (A[i] + B[i] * bx + C[i] * bx * bx);
    em /= n;
    for (int i = 0; i < NPIX; i++) r[i] = (float)(v[i] - (A[i] + B[i] * bx + C[i] * bx * bx) - em);
    write_pgm("omcheck_view.pgm", r, bad, -40, 40);
    printf("saved omcheck_view.pgm (residual at the best T; grey = 0, +-40 DL)\n");
    return 0;
}

/* ======================================================================
 * autot: check on recorded maps how well the temperature can be found from
 * the image's fine detail alone (same method as live7 "auto T").
 *   nuc3 autot nuc2.bin omodel.bin warm
 * ====================================================================== */
static int cmd_autot(int argc, char **argv)
{
    if (argc < 5) { fprintf(stderr, "usage: nuc3 autot nuc2.bin omodel.bin warm\n"); return 1; }
    struct nuchdr nh; struct omhdr oh; struct maphdr mh;
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4), *v = malloc(NPIX * 4), *r = malloc(NPIX * 4);
    float *A = malloc(NPIX * 4), *B = malloc(NPIX * 4), *C = malloc(NPIX * 4);
    uint8_t *bad = malloc(NPIX);
    FILE *f = fopen(argv[2], "rb");
    if (!f || fread(&nh, sizeof nh, 1, f) != 1 || fread(G, 4, NPIX, f) != NPIX ||
        fread(O, 4, NPIX, f) != NPIX || fread(bad, 1, NPIX, f) != NPIX) { fprintf(stderr, "%s: bad nuc file\n", argv[2]); return 1; }
    fclose(f);
    f = fopen(argv[3], "rb");
    if (!f || fread(&oh, sizeof oh, 1, f) != 1 || memcmp(oh.magic, "TVEOMD1", 7) || fread(A, 4, NPIX, f) != NPIX ||
        fread(B, 4, NPIX, f) != NPIX || fread(C, 4, NPIX, f) != NPIX) { fprintf(stderr, "%s: bad omodel file\n", argv[3]); return 1; }
    fclose(f);
    printf("  idx  T70   T(auto, fine detail)   residual sd: model@T70   model@autoT (DL)\n");
    char fn[512];
    double s70 = 0, sau = 0; int nm = 0;
    for (int k = 0; ; k++) {
        snprintf(fn, sizeof fn, "%s_%03d.map", argv[4], k);
        FILE *t = fopen(fn, "rb");
        if (!t) break;
        fclose(t);
        if (load_map(fn, &mh, v) < 0) return 1;
        for (int i = 0; i < NPIX; i++) v[i] *= G[i];
        double sd70 = 0, sdau = 0, T = mh.s.t70_c;
        for (int it = 0; it < 6; it++) {               /* Gauss-Newton on fine detail */
            double x = T - oh.t0, srd = 0, sdd = 0;
            for (int i = 0; i < NPIX; i++) r[i] = v[i] - (A[i] + B[i] * x + C[i] * x * x);
            for (int y = 2; y < H - 2; y += 3)
                for (int xx = 2; xx < W - 2; xx += 2) {
                    int i = y * W + xx;
                    if (bad[i] || bad[i - 2] || bad[i + 2] || bad[i - 2 * W] || bad[i + 2 * W]) continue;
                    float hr = r[i] - 0.25f * (r[i - 2] + r[i + 2] + r[i - 2 * W] + r[i + 2 * W]);
                    if (fabsf(hr) > 12) continue;
                    double d0 = B[i] + 2 * C[i] * x;
                    double hd = d0 - 0.25 * (B[i - 2] + 2 * C[i - 2] * x + B[i + 2] + 2 * C[i + 2] * x +
                                             B[i - 2 * W] + 2 * C[i - 2 * W] * x + B[i + 2 * W] + 2 * C[i + 2 * W] * x);
                    srd += hr * hd; sdd += hd * hd;
                }
            if (sdd > 0) T += srd / sdd;
        }
        for (int pass = 0; pass < 2; pass++) {
            double x = (pass ? T : mh.s.t70_c) - oh.t0, m = 0, s2 = 0; long n = 0;
            for (int i = 0; i < NPIX; i++) if (!bad[i]) { m += v[i] - (A[i] + B[i] * x + C[i] * x * x); n++; }
            m /= n;
            for (int i = 0; i < NPIX; i++) if (!bad[i]) { double e = v[i] - (A[i] + B[i] * x + C[i] * x * x) - m; s2 += e * e; }
            if (pass) sdau = sqrt(s2 / n); else sd70 = sqrt(s2 / n);
        }
        s70 += sd70; sau += sdau; nm++;
        if (k % 10 == 0) printf("  %3d  %3.0f        %6.2f                   %6.1f        %6.1f\n", k, mh.s.t70_c, T, sd70, sdau);
    }
    if (!nm) { fprintf(stderr, "no maps %s_000.map\n", argv[4]); return 1; }
    printf("average residual sd: model@T70 %.1f DL, model@autoT %.1f DL (%d maps)\n", s70 / nm, sau / nm, nm);
    return 0;
}

/* ======================================================================
 * avg: average several maps of the same target into one
 *   nuc3 avg a.map b.map [c.map ...] -o out.map
 * Used for the blackbody taken at two camera positions: the camera's own
 * pattern stays on the same pixels, the blackbody surface pattern moves,
 * so averaging reduces the blackbody's non-uniformity in the result.
 * ====================================================================== */
static int cmd_avg(int argc, char **argv)
{
    const char *out = NULL;
    int nin = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-o") && i + 1 < argc) { out = argv[++i]; continue; }
        nin++;
    }
    if (!out || nin < 2) { fprintf(stderr, "usage: nuc3 avg a.map b.map [...] -o out.map\n"); return 1; }
    struct maphdr h, h0;
    float *v = malloc(NPIX * 4), *a = calloc(NPIX, 4);
    double tsum = 0, nsum = 0;
    int k = 0;
    for (int i = 2; i < argc; i++) {
        if (!strcmp(argv[i], "-o")) { i++; continue; }
        if (load_map(argv[i], &h, v) < 0) return 1;
        if (k == 0) h0 = h;
        else if (h.s.tint != h0.s.tint || h.s.gfid != h0.s.gfid || h.s.gsk != h0.s.gsk || h.s.gain_reg != h0.s.gain_reg) {
            fprintf(stderr, "%s: sensor settings differ from %s - not averaged\n", argv[i], argv[2]); return 1;
        }
        double m = 0;
        for (int j = 0; j < NPIX; j++) { a[j] += v[j]; m += v[j]; }
        printf("  %-24s mean %8.1f DL  target %.1f C  T70 %.0f C\n", argv[i], m / NPIX, h.target_c, h.s.t70_c);
        tsum += h.s.t70_c; nsum += h.noise_dl;
        k++;
    }
    for (int j = 0; j < NPIX; j++) a[j] /= k;
    h0.nframes *= k;                      /* frames behind the average */
    h0.noise_dl = (float)(nsum / k);
    h0.s.t70_c = (float)(tsum / k);
    if (save_map(out, &h0, a) < 0) return 1;
    printf("averaged %d maps -> %s\n", k, out);
    return 0;
}

/* ======================================================================
 * bands: what is left on a UNIFORM target after the full live7 correction
 *   nuc3 bands nuc2.bin target.map [omodel.bin]
 * With omodel.bin: G*raw - p(T), T found from the image (auto T, as live7).
 * Without:         G*raw + O (fixed 2-point offset).
 * The column profile (mean of each column) is split by width:
 *   narrow : lines 1..16 px wide   (the live 'c' filter removes these)
 *   medium : bands 17..120 px wide (NOT removed by 'c')
 *   wide   : slower shading across the image
 * Same for rows. Writes bands_view.pgm (+-20 DL) and bands_cols.csv.
 * Run it on a lens-cap map (room temp) and on a warm uniform map: if the
 * medium bands grow with the target temperature they come from the gain G,
 * if they stay the same they come from the offset.
 * ====================================================================== */
static void smooth_d(const double *in, double *out, int n, int r)
{
    for (int i = 0; i < n; i++) {
        double s = 0; int c = 0;
        for (int k = i - r; k <= i + r; k++) if (k >= 0 && k < n) { s += in[k]; c++; }
        out[i] = s / c;
    }
}

static int cmd_bands(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: nuc3 bands nuc2.bin target.map [omodel.bin]\n"); return 1; }
    struct nuchdr nh; struct maphdr mh; struct omhdr oh;
    float *G = malloc(NPIX * 4), *O = malloc(NPIX * 4), *v = malloc(NPIX * 4), *r = malloc(NPIX * 4);
    float *A = NULL, *B = NULL, *C = NULL;
    uint8_t *bad = malloc(NPIX);
    FILE *f = fopen(argv[2], "rb");
    if (!f || fread(&nh, sizeof nh, 1, f) != 1 || memcmp(nh.magic, "TVENUC2", 7) || fread(G, 4, NPIX, f) != NPIX ||
        fread(O, 4, NPIX, f) != NPIX || fread(bad, 1, NPIX, f) != NPIX) { fprintf(stderr, "%s: bad nuc file\n", argv[2]); return 1; }
    fclose(f);
    if (load_map(argv[3], &mh, v) < 0) return 1;
    for (int i = 0; i < NPIX; i++) v[i] *= G[i];
    double T = NAN;
    if (argc >= 5) {
        A = malloc(NPIX * 4); B = malloc(NPIX * 4); C = malloc(NPIX * 4);
        f = fopen(argv[4], "rb");
        if (!f || fread(&oh, sizeof oh, 1, f) != 1 || memcmp(oh.magic, "TVEOMD1", 7) || fread(A, 4, NPIX, f) != NPIX ||
            fread(B, 4, NPIX, f) != NPIX || fread(C, 4, NPIX, f) != NPIX) { fprintf(stderr, "%s: bad omodel file\n", argv[4]); return 1; }
        fclose(f);
        T = isnan(mh.s.t70_c) ? oh.t0 : mh.s.t70_c;
        for (int it = 0; it < 8; it++) {             /* auto T on fine detail */
            double x = T - oh.t0, srd = 0, sdd = 0;
            for (int i = 0; i < NPIX; i++) r[i] = v[i] - (A[i] + B[i] * x + C[i] * x * x);
            for (int y = 2; y < H - 2; y += 3)
                for (int xx = 2; xx < W - 2; xx += 2) {
                    int i = y * W + xx;
                    if (bad[i] || bad[i - 2] || bad[i + 2] || bad[i - 2 * W] || bad[i + 2 * W]) continue;
                    float hr = r[i] - 0.25f * (r[i - 2] + r[i + 2] + r[i - 2 * W] + r[i + 2 * W]);
                    if (fabsf(hr) > 12) continue;
                    double hd = (B[i] + 2 * C[i] * x) - 0.25 * (B[i - 2] + 2 * C[i - 2] * x + B[i + 2] + 2 * C[i + 2] * x +
                                 B[i - 2 * W] + 2 * C[i - 2 * W] * x + B[i + 2 * W] + 2 * C[i + 2 * W] * x);
                    srd += hr * hd; sdd += hd * hd;
                }
            if (sdd > 0) T += srd / sdd;
        }
        double x = T - oh.t0;
        for (int i = 0; i < NPIX; i++) v[i] -= A[i] + B[i] * x + C[i] * x * x;
    } else {
        for (int i = 0; i < NPIX; i++) v[i] += O[i];
    }
    double m, sd;
    stats(v, bad, &m, &sd);
    for (int i = 0; i < NPIX; i++) v[i] -= (float)m;

    static double col[W], row[H], s8[W > H ? W : H], s120[W > H ? W : H];
    static int ncol[W], nrow[H];
    for (int y = 0; y < H; y++)
        for (int x = 0; x < W; x++) {
            int i = y * W + x;
            if (bad[i]) continue;
            col[x] += v[i]; ncol[x]++; row[y] += v[i]; nrow[y]++;
        }
    for (int x = 0; x < W; x++) col[x] /= ncol[x] ? ncol[x] : 1;
    for (int y = 0; y < H; y++) row[y] /= nrow[y] ? nrow[y] : 1;

    double cn = 0, cm = 0, cw = 0, rn = 0, rm = 0, rw = 0; int nc = 0, nr = 0;
    double medc[W];
    smooth_d(col, s8, W, 8); smooth_d(col, s120, W, 60);
    FILE *csv = fopen("bands_cols.csv", "w");
    if (csv) fprintf(csv, "x,column_mean,narrow,medium,wide\n");
    for (int x = 0; x < W; x++) {
        double nn = col[x] - s8[x], md = s8[x] - s120[x], wd = s120[x];
        medc[x] = md;
        if (csv) fprintf(csv, "%d,%.2f,%.2f,%.2f,%.2f\n", x, col[x], nn, md, wd);
        if (x >= 60 && x < W - 60) { cn += nn * nn; cm += md * md; cw += wd * wd; nc++; }
    }
    if (csv) fclose(csv);
    smooth_d(row, s8, H, 8); smooth_d(row, s120, H, 60);
    for (int y = 60; y < H - 60; y++) {
        double nn = row[y] - s8[y], md = s8[y] - s120[y], wd = s120[y];
        rn += nn * nn; rm += md * md; rw += wd * wd; nr++;
    }
    printf("%s on %s (target %.1f C, T70 %.0f C)%s\n", argc >= 5 ? "full live7 correction (G + O(T) + auto T)" : "2-point (fixed O)",
           argv[3], mh.target_c, mh.s.t70_c, "");
    if (argc >= 5) printf("  auto T = %.2f C\n", T);
    printf("  total residual sd : %6.2f DL\n", sd);
    printf("                      narrow(1-16px)  medium(17-120px)  wide/shading\n");
    printf("  vertical (columns): %8.2f        %8.2f          %8.2f  DL\n", sqrt(cn / nc), sqrt(cm / nc), sqrt(cw / nc));
    printf("  horizontal (rows) : %8.2f        %8.2f          %8.2f  DL\n", sqrt(rn / nr), sqrt(rm / nr), sqrt(rw / nr));
    printf("  strongest medium vertical bands:");
    for (int k = 0; k < 6; k++) {
        int best = -1; double bv = 0;
        for (int x = 60; x < W - 60; x++) if (fabs(medc[x]) > bv) { bv = fabs(medc[x]); best = x; }
        if (best < 0) break;
        printf(" x%d:%+.1f", best, medc[best]);
        for (int x = best - 20; x <= best + 20; x++) if (x >= 0 && x < W) medc[x] = 0;   /* next band */
    }
    printf("\n");
    write_pgm("bands_view.pgm", v, bad, -20, 20);
    printf("saved bands_view.pgm (grey 0, white +20, black -20 DL) and bands_cols.csv\n");
    return 0;
}

/* ======================================================================
 * same: are two maps of the same target stable?   (used by lab_bb.bat)
 *   nuc3 same a.map b.map [limit_dl=5]
 * Prints both means and the difference; exit code 0 = |diff| < limit.
 * ====================================================================== */
static int cmd_same(int argc, char **argv)
{
    if (argc < 4) { fprintf(stderr, "usage: nuc3 same a.map b.map [limit_dl]\n"); return 2; }
    double lim = argc >= 5 ? atof(argv[4]) : 5.0;
    struct maphdr ha, hb;
    float *a = malloc(NPIX * 4), *b = malloc(NPIX * 4);
    if (load_map(argv[2], &ha, a) < 0 || load_map(argv[3], &hb, b) < 0) return 2;
    double ma = 0, mb = 0;
    for (int i = 0; i < NPIX; i++) { ma += a[i]; mb += b[i]; }
    ma /= NPIX; mb /= NPIX;
    double d = mb - ma;
    printf("mean %.1f -> %.1f DL, change %+.1f DL (limit %.1f): %s\n", ma, mb, d, lim,
           fabs(d) < lim ? "STABLE" : "not stable yet");
    return fabs(d) < lim ? 0 : 1;
}

int main(int argc, char **argv)
{
    if (argc >= 2 && !strcmp(argv[1], "map"))     return cmd_map(argc, argv);
#ifdef __linux__
    if (argc >= 2 && !strcmp(argv[1], "capture")) return cmd_capture(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "warmup"))  return cmd_warmup(argc, argv);
#endif
    if (argc >= 2 && !strcmp(argv[1], "calc"))    return cmd_calc(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "check"))   return cmd_check(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "apply"))   return cmd_apply(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "noise"))   return cmd_noise(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "stripes")) return cmd_stripes(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "omodel"))  return cmd_omodel(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "autot"))   return cmd_autot(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "bands"))   return cmd_bands(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "same"))    return cmd_same(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "t70"))     return cmd_t70(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "avg"))     return cmd_avg(argc, argv);
    if (argc >= 2 && !strcmp(argv[1], "omcheck")) return cmd_omcheck(argc, argv);
    fprintf(stderr,
        "nuc3 - 2-point NUC + offset model tool\n"
        "  nuc2 map file.raw -o cold1.map -t 25 [-n 128] [-s 30] [-p COM5|/dev/ttyACM0|none]\n"
#ifdef __linux__
        "  nuc2 capture -o cold1.map -t 25 [-n 128] [-d /dev/video0] [-p /dev/ttyACM0|none]\n"
        "  nuc2 warmup -o warm [-i 30] [-m 60] [-n 64] [-d /dev/video0] [-p /dev/ttyACM0|none]\n"
#endif
        "  nuc2 calc cold1.map hot.map [cold2.map] -o nuc2.bin\n"
        "  nuc2 check nuc2.bin test.map [ref.map]\n"
        "  nuc2 apply nuc2.bin scene.raw [-f frame] [-r ref.map]\n"
        "  nuc2 noise file.raw [-s 30] [-n 128]\n"
        "  nuc2 stripes nuc2.bin wall.map\n"
        "  nuc3 omodel nuc2.bin warm [-o omodel.bin]   (warm_000.map ... from warmup)\n"
        "  nuc3 omcheck nuc2.bin omodel.bin test.map   (lens-cap map taken with -p)\n"
        "  nuc3 autot nuc2.bin omodel.bin warm\n"
        "  nuc3 bands nuc2.bin target.map [omodel.bin]\n"
        "  nuc3 same a.map b.map [limit_dl]\n"
        "  nuc3 t70 -p COM5|/dev/ttyACM0      (board temperature every 2 s)\n");
    return 1;
}
