@echo off
cls

rem "name" and "dirout" are named according to the testcase

set name=CaseDambreakVal2D
set dirout=%name%_out

rem "executables" are renamed and called from their directory

set tools=../../EXECS_WIN
set gencase="%tools%/GenCase4_win64.exe"
set dualsphysics="%tools%/DualSPHysics_win64.exe"
set partvtk="%tools%/PartVTK4_win64.exe"
set partvtkout="%tools%/PartVTKOut4_win64.exe"
set measuretool="%tools%/MeasureTool4_win64.exe"
set isosurface="%tools%/IsoSurface4_win64.exe"

rem "dirout" is created to store results or it is removed if it already exists

if exist %dirout% del /Q %dirout%\*.*
if not exist %dirout% mkdir %dirout%

rem CODES are executed according the selected parameters of execution in this testcase

%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

%dualsphysics% %dirout%/%name% %dirout% -svres -gpu
if not "%ERRORLEVEL%" == "0" goto fail

%partvtk% -dirin %dirout% -filexml %dirout%/%name%.xml -savevtk %dirout%/PartFluid -onlytype:-all,fluid -vars:+idp,+vel,+rhop,+press,+vor
if not "%ERRORLEVEL%" == "0" goto fail

%partvtkout% -dirin %dirout% -filexml %dirout%/%name%.xml -savevtk %dirout%/PartFluidOut -SaveResume %dirout%/ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

%isosurface% -dirin %dirout% -saveslice %dirout%/Slices 
if not "%ERRORLEVEL%" == "0" goto fail

:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause

