@echo off
REM Run in the SAME activated Python environment used by the Spyder kernel.
python -m pip install -r "%~dp0requirements.txt"
if errorlevel 1 exit /b 1
echo Installation finished. Restart the Spyder kernel before running analysis.
pause
