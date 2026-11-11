@echo off
REM ------------------------------------------------------------------------
REM lab_bb2.bat - blackbody captures at ONE chamber temperature, with a
REM               COVER check (does the blackbody fill the whole image?)
REM
REM Usage (cmd, in the folder with lab_bb2.bat, nuc3.exe, ffmpeg.exe,
REM        nuc2_c1.bin):
REM     lab_bb2.bat setup1       setup check only (stop after the COVER result)
REM     lab_bb2.bat c12p5        chamber 12.5 C
REM     lab_bb2.bat c40          chamber 40 C
REM     lab_bb2.bat cm10         chamber -10 C   (m = minus)
REM
REM Steps: BB 20 -> BB 35 -> COVER check -> BB 50 -> BB 20 again (drift).
REM Every step: you set the BB, type VTEMP, the script records 160 frames
REM 60 s apart until two maps differ by < 5 DL (stable), max 8 tries.
REM Never overwrites a file. Log: <session>\log.csv
REM Run the same command again after a stop: finished BB steps are skipped.
REM Read-only on the camera (frames + register reads only).
REM ------------------------------------------------------------------------
setlocal enabledelayedexpansion
set CAM=FullScale_USB3-A3 v3.030
set PORT=COM5
set FRAMES=160
set WAIT_S=60
set MAXTRY=8

set S=%~1
if "%S%"=="" (
  echo usage: lab_bb2.bat c12p5
  exit /b 1
)
for %%f in (nuc3.exe ffmpeg.exe nuc2_c1.bin) do if not exist %%f ( echo %%f is not in this folder & exit /b 1 )
if not exist "%S%" mkdir "%S%"
if not exist "%S%\log.csv" echo time,file,bb_c,vtemp_c,result>"%S%\log.csv"

REM Resume: BB steps already stable in this session (good.txt) are skipped
if exist "%S%\good.txt" for /f "tokens=1,2" %%a in (%S%\good.txt) do (set GOOD_%%a=%%b& echo   already done: %%a  %%b)

echo Session %S%  -^>  folder %CD%\%S%
echo Before starting: chamber at its temperature, VTEMP stable (5 min change below 0.1 C).
set /p CH=  Chamber set temperature (C):
echo chamber_setpoint,%CH%>>"%S%\log.csv"

call :bb 20 bb20 || goto failed
call :bb 35 bb35 || goto failed

echo.
echo ==================================================================
echo  COVER CHECK: does the blackbody fill the whole image?
nuc3.exe cover nuc2_c1.bin "%GOOD_bb20%" "%GOOD_bb35%" > "%S%\cover.txt"
set RC=%errorlevel%
type "%S%\cover.txt"
copy /y cover.pgm "%S%\cover.pgm" >nul
if not %RC%==0 (
  echo.
  echo  !! STOP: the blackbody does not fill the image. Data from this setup is
  echo     not usable. Move the BB closer / use a bigger BB, then run this
  echo     command again with a NEW session name.
  exit /b 1
)
if /i "%S:~0,5%"=="setup" (
  echo.
  echo  Setup is OK. Now start the real session, e.g.  lab_bb2.bat c12p5
  exit /b 0
)

call :bb 50 bb50 || goto failed
echo.
echo After the hot BB the lens front is warm: wait 5 minutes before the last step.
call :bb 20 bb20end || goto failed

echo.
echo ==================================================================
echo  Session %S% done. Stable captures:
type "%S%\good.txt"
echo  Copy the whole folder %S% back for analysis.
exit /b 0

:failed
echo.
echo  !! stopped - fix the problem and run the same command again
echo     (finished files are kept, nothing is overwritten)
exit /b 1

REM ---------------------------------------------------------------- :bb T TAG
:bb
set T=%1
set TAG=%2
if defined GOOD_%TAG% ( echo. & echo  BB %T% C  ^(%TAG%^) already done - skipped & exit /b 0 )
echo.
echo ==================================================================
echo  BB %T% C   (%TAG%)
echo   - set the blackbody to %T% C (camera looks at it, BB fills the image)
echo   - wait until the BB display is stable
pause
set VT=
set /p VT=  VTEMP / FPA temperature now (C):
set PREV=
set K=0
:try
set /a K+=1
if %K% GTR %MAXTRY% (
  echo   !! not stable after %MAXTRY% captures - check BB / camera / chamber
  exit /b 1
)
set N=1
:freename
if exist "%S%\%S%_%TAG%_%N%.raw" ( set /a N+=1 & goto freename )
set F=%S%\%S%_%TAG%_%N%
echo   capture %K% -^> %F%.map
ffmpeg.exe -nostdin -hide_banner -loglevel error -y -f dshow -video_size 640x480 -pixel_format yuyv422 -rtbufsize 200M -i video="%CAM%" -frames:v %FRAMES% -f rawvideo "%F%.raw"
if errorlevel 1 ( echo   !! ffmpeg failed - camera name / camera busy & exit /b 1 )
nuc3.exe map "%F%.raw" -o "%F%.map" -t %T% -p %PORT% > "%F%.txt" 2>&1
if errorlevel 1 ( echo   !! nuc3 map failed - see %F%.txt & exit /b 1 )
for /f "tokens=2" %%m in ('findstr /b /c:"mean " "%F%.txt"') do set MEAN=%%m
for /f "tokens=6" %%t in ('findstr /c:"T70=" "%F%.txt"') do set T70=%%t
echo     mean !MEAN! DL   !T70!
if "%PREV%"=="" (
  echo %time%,%F%.map,%T%,%VT%,first>>"%S%\log.csv"
  goto next
)
nuc3.exe same "%PREV%.map" "%F%.map" 5 > "%F%_same.txt"
set RC=%errorlevel%
type "%F%_same.txt"
for /f "delims=" %%r in ('type "%F%_same.txt"') do set RES=%%r
echo %time%,%F%.map,%T%,%VT%,!RES!>>"%S%\log.csv"
if %RC%==0 (
  echo   OK: BB %T% C stable -^> use %F%.map
  echo %TAG% %F%.map>>"%S%\good.txt"
  set GOOD_%TAG%=%F%.map
  exit /b 0
)
:next
set PREV=%F%
echo   waiting %WAIT_S% s ...
timeout /t %WAIT_S% /nobreak >nul
goto try
