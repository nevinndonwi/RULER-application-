@echo off
TITLE Setup and Run - FINAL_BIO_PROJECT
CLS

echo ======================================================
echo  FINAL BIO PROJECT - Environment Setup
echo ======================================================

:: 1. DETERMINE PYTHON COMMAND
:: We check if 'python' is version 3. If not, we check 'py -3'.

set PYTHON_CMD=python

:: Check if 'python' exists and is version 3
python --version 2>NUL | findstr /C:"Python 3" >NUL
IF %ERRORLEVEL% NEQ 0 (
    :: 'python' is not found or is not Python 3. Let's try the 'py' launcher.
    py -3 --version 2>NUL | findstr /C:"Python 3" >NUL
    IF %ERRORLEVEL% EQU 0 (
        set PYTHON_CMD=py -3
        echo Found Python 3 via Python Launcher.
    ) ELSE (
        goto :INSTALL_PYTHON
    )
)

goto :CHECK_DEPS

:INSTALL_PYTHON
echo Python 3 is not found on this computer.
echo.
echo Attempting to install Python 3.11 automatically...
:: This installs specifically Python 3.11
winget install -e --id Python.Python.3.11

if %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] Automatic installation failed.
    echo Please download Python 3 from python.org and install it manually.
    echo IMPORTANT: Check the box "Add Python to PATH" during installation.
    PAUSE
    EXIT /B
)
echo.
echo Python 3 installed successfully!
echo Please RESTART this script to continue (close this window and run it again).
PAUSE
EXIT

:CHECK_DEPS
echo Using command: %PYTHON_CMD%
echo.

:: 2. Check for requirements.txt
IF NOT EXIST "requirements.txt" (
    echo [ERROR] requirements.txt is missing!
    echo Please make sure requirements.txt is in the same folder as this script.
    PAUSE
    EXIT /B
)

:: 3. Install Libraries
echo Installing required libraries...
:: We use the specific python command we found to run pip
%PYTHON_CMD% -m pip install -r requirements.txt

IF %ERRORLEVEL% NEQ 0 (
    echo.
    echo [ERROR] Failed to install libraries. Check your internet connection.
    PAUSE
    EXIT /B
)

echo.
echo ======================================================
echo  Setup Complete. Launching FINAL_BIO_PROJECT.py...
echo ======================================================
echo.

:: 4. Run the project
%PYTHON_CMD% FINAL_BIO_PROJECT.py

:: Keep window open
PAUSE