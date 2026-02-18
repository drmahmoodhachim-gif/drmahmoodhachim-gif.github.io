@echo off
title RNA-seq Tutorial - Running Analysis
echo.
echo ========================================
echo   RNA-seq Pipeline for Students
echo   (No Python needed - R only!)
echo ========================================
echo.

set RSCRIPT=
where Rscript >nul 2>nul
if %ERRORLEVEL% EQU 0 (
    set RSCRIPT=Rscript
)

REM If not in PATH, try common installation locations
if "%RSCRIPT%"=="" (
    if exist "C:\Program Files\R\R-4.4.2\bin\Rscript.exe" set RSCRIPT="C:\Program Files\R\R-4.4.2\bin\Rscript.exe"
)
if "%RSCRIPT%"=="" (
    if exist "C:\Program Files\R\R-4.4.1\bin\Rscript.exe" set RSCRIPT="C:\Program Files\R\R-4.4.1\bin\Rscript.exe"
)
if "%RSCRIPT%"=="" (
    if exist "C:\Program Files\R\R-4.3.3\bin\Rscript.exe" set RSCRIPT="C:\Program Files\R\R-4.3.3\bin\Rscript.exe"
)
if "%RSCRIPT%"=="" (
    for /d %%D in ("C:\Program Files\R\R-*") do (
        if exist "%%D\bin\Rscript.exe" set RSCRIPT="%%D\bin\Rscript.exe"
    )
)

REM Try 32-bit Program Files
if "%RSCRIPT%"=="" (
    if exist "C:\Program Files (x86)\R\R-4.4.2\bin\Rscript.exe" set RSCRIPT="C:\Program Files (x86)\R\R-4.4.2\bin\Rscript.exe"
)
if "%RSCRIPT%"=="" (
    for /d %%D in ("C:\Program Files (x86)\R\R-*") do (
        if exist "%%D\bin\Rscript.exe" set RSCRIPT="%%D\bin\Rscript.exe"
    )
)

if "%RSCRIPT%"=="" (
    echo ERROR: R was not found on your computer.
    echo.
    echo Do one of these:
    echo.
    echo 1. Install R: https://cran.r-project.org/
    echo    - Download R for Windows
    echo    - During install, CHECK "Add R to PATH"
    echo    - Restart your computer, then try again
    echo.
    echo 2. Use RStudio instead:
    echo    - Install RStudio from https://posit.co/download/rstudio-desktop/
    echo    - Open RStudio, then File ^> Open ^> run_all_in_R.R
    echo    - Click the "Source" button
    echo.
    pause
    exit /b 1
)

echo Found R. Running analysis... (1-2 minutes)
echo.

cd /d "%~dp0"
%RSCRIPT% run_all_in_R.R

echo.
echo ========================================
echo   Analysis complete!
echo ========================================
echo.
echo Open these folders to see your results:
echo   - 05_deg_results   (volcano plot, heatmap)
echo   - 06_systems_biology
echo.
echo Double-click volcano_plot.png to view it!
echo.
pause
