@echo off
REM ============================================
REM  Setup Script for Windows
REM  Cell Annealing Tool - Symmetric Relaxation of 2D Cellular Networks
REM  Requires Python 3.14+ managed via uv
REM ============================================

setlocal enabledelayedexpansion

echo ============================================
echo  Cell Annealing Tool - Environment Setup
echo ============================================
echo.

REM ---- Step 1: Check uv, install Python 3.14.6, create .venv ----
echo [1/4] Checking uv installation...
uv --version >nul 2>&1
if %errorlevel% neq 0 (
    echo [ERROR] uv is not found. Please install uv first:
    echo         https://docs.astral.sh/uv/getting-started/installation/
    pause
    exit /b 1
)
echo [OK] uv detected:
uv --version
echo.
echo       Installing Python 3.14.6 and creating .venv...
uv python install 3.14.6
if not exist ".venv" (
    uv venv .venv --python 3.14.6
)
echo [OK] .venv is ready.

REM ---- Step 2: Install Python dependencies ----
echo.
echo [2/4] Installing Python dependencies...
uv pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo [ERROR] Failed to install Python dependencies.
    pause
    exit /b 1
)
echo [OK] Python dependencies installed.

REM ---- Step 3: Check / Setup R ----
echo.
echo [3/4] Checking R environment...

REM Try using bundled R_Dist first
set "R_CMD="
if exist "%~dp0R_Dist\bin\x64\R.exe" (
    echo [INFO] Found bundled R_Dist, 64-bit.
    set "R_CMD=%~dp0R_Dist\bin\x64\R.exe"
    set "RSCRIPT_CMD=%~dp0R_Dist\bin\x64\Rscript.exe"
    set "R_HOME=%~dp0R_Dist"
) else if exist "%~dp0R_Dist\bin\R.exe" (
    echo [INFO] Found bundled R_Dist.
    set "R_CMD=%~dp0R_Dist\bin\R.exe"
    set "RSCRIPT_CMD=%~dp0R_Dist\bin\Rscript.exe"
    set "R_HOME=%~dp0R_Dist"
) else (
    REM Fallback: check system R
    where R >nul 2>&1
    if !errorlevel! equ 0 (
        echo [INFO] Using system R.
        set "R_CMD=R"
        set "RSCRIPT_CMD=Rscript"
    ) else (
        echo.
        echo [WARNING] R is not found - neither bundled R_Dist nor system R.
        echo.
        echo  Would you like to download and install R to R_Dist automatically?
        echo   Requires about 300 MB, internet connection needed.
        echo.
        choice /C YN /M "Download and install R to R_Dist?"
        if !errorlevel! equ 1 (
            echo.
            echo [INFO] Downloading R 4.4.2 for Windows...
            echo  Downloading from CRAN...
            powershell -Command "& {[Net.ServicePointManager]::SecurityProtocol = [Net.SecurityProtocolType]::Tls12; Invoke-WebRequest -Uri 'https://cloud.r-project.org/bin/windows/base/R-4.4.2-win.exe' -OutFile '%TEMP%\R-installer.exe'}"
            if !errorlevel! neq 0 (
                echo [ERROR] Failed to download R. Please install R manually from:
                echo         https://www.r-project.org/
                pause
                exit /b 1
            )
            echo [INFO] Installing R to R_Dist directory...
            "%TEMP%\R-installer.exe" /VERYSILENT /DIR="%~dp0R_Dist" /NOICONS
            if !errorlevel! neq 0 (
                echo [ERROR] Failed to install R. Please install R manually from:
                echo         https://www.r-project.org/
                pause
                exit /b 1
            )
            del "%TEMP%\R-installer.exe"
            echo [OK] R installed to R_Dist.
            if exist "%~dp0R_Dist\bin\x64\R.exe" (
                set "R_CMD=%~dp0R_Dist\bin\x64\R.exe"
                set "RSCRIPT_CMD=%~dp0R_Dist\bin\x64\Rscript.exe"
                set "R_HOME=%~dp0R_Dist"
            ) else (
                set "R_CMD=%~dp0R_Dist\bin\R.exe"
                set "RSCRIPT_CMD=%~dp0R_Dist\bin\Rscript.exe"
                set "R_HOME=%~dp0R_Dist"
            )
        ) else (
            echo [INFO] Skipped R installation. You can manually install R from:
            echo         https://www.r-project.org/
            echo  After installing, run the following in R console:
            echo      install.packages(c("conicfit", "sp", "shotGroups"))
            echo.
            pause
            exit /b 0
        )
    )
)

echo [OK] R detected:
%R_CMD% --version

REM ---- Step 4: Install R packages ----
echo.
echo [4/4] Installing required R packages (conicfit, sp, shotGroups)...

REM Configure R_HOME environment for the process
if defined R_HOME set "R_HOME=%R_HOME%"

%RSCRIPT_CMD% -e "if (!require('conicfit')) install.packages('conicfit', repos='https://cloud.r-project.org'); if (!require('sp')) install.packages('sp', repos='https://cloud.r-project.org'); if (!require('shotGroups')) install.packages('shotGroups', repos='https://cloud.r-project.org'); cat('[OK] All R packages installed.\n')"
if %errorlevel% neq 0 (
    echo [WARNING] Some R packages may not have installed correctly.
    echo  You can manually install them later in R console:
    echo      install.packages(c("conicfit", "sp", "shotGroups"))
) else (
    echo [OK] All R packages installed successfully.
)

echo.
echo ============================================
echo  Setup complete! You can now run:
echo      uv run python energy_annealing_main.py
echo ============================================
pause
