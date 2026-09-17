@echo off
REM ============================================
REM  Build Script for Windows - PyInstaller EXE
REM  Cell Annealing Tool
REM  Uses the project .venv (Python 3.14) created by setup.bat
REM ============================================

title Build Cell Annealing EXE

echo ============================================
echo  Cell Annealing Tool - PyInstaller Build
echo ============================================
echo.

REM Switch to script directory
cd /d "%~dp0"

REM Check .venv exists
if not exist ".venv\Scripts\python.exe" (
    echo [ERROR] .venv not found. Please run setup.bat first.
    pause
    exit /b 1
)

REM Install dependencies (skip if already installed)
echo [1/3] Ensuring dependencies are installed...
uv pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo [ERROR] Failed to install dependencies.
    pause
    exit /b 1
)

REM Install PyInstaller into .venv
echo [2/3] Ensuring PyInstaller is installed...
uv pip install pyinstaller
if %errorlevel% neq 0 (
    echo [ERROR] Failed to install PyInstaller.
    pause
    exit /b 1
)

REM Build
echo [3/3] Building...
echo.
echo Building single-file EXE (no console window)...
".venv\Scripts\python.exe" -m PyInstaller --onefile --windowed --name "CellAnnealing" ^
    --add-data "annealing;annealing" ^
    --add-data "cell;cell" ^
    --add-data "utillib;utillib" ^
    --add-data "initVoronoi.py;." ^
    --hidden-import "scipy.spatial" ^
    --hidden-import "scipy.optimize" ^
    --hidden-import "openpyxl" ^
    --hidden-import "pyenvelope" ^
    energy_annealing_main.py

if %errorlevel% equ 0 (
    echo.
    echo ============================================
    echo  Build succeeded!
    echo   Output: dist\CellAnnealing.exe
    echo ============================================
) else (
    echo.
    echo  Build failed. Check the error messages above.
)

pause
