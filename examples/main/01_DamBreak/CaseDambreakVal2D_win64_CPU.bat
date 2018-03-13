@echo off

REM "name" and "dirout" are named according to the testcase

set name=CaseDambreakVal2D
set dirout=%name%_out

REM "executables" are renamed and called from their directory

set dirbin=../../../bin/windows
set gencase="%dirbin%/GenCase4_win64.exe"
set dualsphysicscpu="%dirbin%/DualSPHysics4.2CPU_win64.exe"
set dualsphysicsgpu="%dirbin%/DualSPHysics4.2_win64.exe"
set boundaryvtk="%dirbin%/BoundaryVTK4_win64.exe"
set partvtk="%dirbin%/PartVTK4_win64.exe"
set partvtkout="%dirbin%/PartVTKOut4_win64.exe"
set measuretool="%dirbin%/MeasureTool4_win64.exe"
set computeforces="%dirbin%/ComputeForces4_win64.exe"
set isosurface="%dirbin%/IsoSurface4_win64.exe"
set measureboxes="%dirbin%/MeasureBoxes4_win64.exe"

REM "dirout" is created to store results or it is removed if it already exists

if exist %dirout% del /Q %dirout%\*.*
if not exist %dirout% mkdir %dirout%

REM CODES are executed according the selected parameters of execution in this testcase

%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

%dualsphysicscpu% -cpu %dirout%/%name% %dirout% -svres
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

