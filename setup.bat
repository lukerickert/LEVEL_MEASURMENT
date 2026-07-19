@echo off
setlocal enabledelayedexpansion
title Precision Level Straightness - Setup and Run
cd /d "%~dp0"

echo ==================================================
echo    Precision Level Straightness  -  Setup and Run
echo ==================================================
echo.

set "PY="

REM ---- 1. Look for an existing Python (ignore the Microsoft Store stub) ----
for /f "delims=" %%P in ('where python 2^>nul') do (
    echo %%P| findstr /I "WindowsApps" >nul
    if errorlevel 1 if not defined PY set "PY=%%P"
)
if not defined PY for /f "delims=" %%P in ('where py 2^>nul') do if not defined PY set "PY=%%P"

REM ---- 2. If none found, install Python 3.13 with winget ----
if not defined PY (
    where winget >nul 2>&1 || (
        echo winget is not available on this PC.
        echo Please install Python manually from https://www.python.org/downloads/
        echo and tick "Add python.exe to PATH" on the first installer screen,
        echo then run this file again.
        echo.
        pause
        exit /b 1
    )
    echo No Python found. Installing Python 3.13 via winget...
    echo Click "Yes" if Windows asks for permission.
    echo.
    winget install -e --id Python.Python.3.13 --accept-source-agreements --accept-package-agreements
    echo.
    echo Locating the newly installed Python...
    for /f "delims=" %%P in ('dir /b /s "%LOCALAPPDATA%\Programs\Python\python.exe" 2^>nul') do if not defined PY set "PY=%%P"
    if not defined PY for /f "delims=" %%P in ('dir /b /s "%ProgramFiles%\Python3*\python.exe" 2^>nul') do if not defined PY set "PY=%%P"
)

REM ---- 3. Give up gracefully if still not found ----
if not defined PY (
    echo.
    echo Could not locate Python automatically.
    echo Please CLOSE this window, open it again, and double-click setup.bat once more.
    echo The second run will detect the Python that was just installed.
    echo.
    pause
    exit /b 1
)

echo Using Python at: !PY!
"!PY!" --version
echo.

REM ---- 4. Install the two required libraries (tkinter is built in) ----
echo Installing required libraries (numpy, matplotlib)...
"!PY!" -m pip install --upgrade pip
"!PY!" -m pip install numpy matplotlib
if errorlevel 1 (
    echo.
    echo Library installation failed. Check your internet connection and try again.
    echo.
    pause
    exit /b 1
)
echo.

REM ---- 5. Download the script if it isn't already next to this bat ----
if not exist "%~dp0level_straightness.py" (
    echo Downloading level_straightness.py ...
    powershell -NoProfile -ExecutionPolicy Bypass -Command "Invoke-WebRequest -Uri 'https://raw.githubusercontent.com/lukerickert/LEVEL_MEASURMENT/main/level_straightness.py' -OutFile '%~dp0level_straightness.py'"
    if errorlevel 1 (
        echo Download failed. Check your internet connection.
        echo.
        pause
        exit /b 1
    )
)

REM ---- 6. Run it ----
echo.
echo Launching the tool...
echo.
"!PY!" "%~dp0level_straightness.py"

echo.
echo The program has closed.
pause
endlocal
