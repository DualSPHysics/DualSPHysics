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
set flowtool="%dirbin%/FlowTool4_win64.exe"

REM "dirout" is created to store results or it is removed if it already exists

if exist %dirout% rd /s /q %dirout%
mkdir %dirout%
if not "%ERRORLEVEL%" == "0" goto fail
set diroutdata=%dirout%\data
mkdir %diroutdata%

REM CODES are executed according the selected parameters of execution in this testcase

REM Executes GenCase4 to create initial files for simulation.
%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes DualSPHysics to simulate SPH method.
%dualsphysicsgpu% -gpu %dirout%/%name% %dirout% -dirdataout data -svres
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes PartVTK4 to create VTK files with particles.
set dirout2=%dirout%\particles
mkdir %dirout2%
%partvtk% -dirin %diroutdata% -filexml %dirout%/%name%.xml -savevtk %dirout2%/PartFluid -onlytype:-all,fluid -vars:+idp,+vel,+rhop,+press,+vor
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes PartVTKOut4 to create VTK files with excluded particles.
%partvtkout% -dirin %diroutdata% -filexml %dirout%/%name%.xml -savevtk %dirout2%/PartFluidOut -SaveResume %dirout%/ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes IsoSurface4 to create VTK files with slices of surface.
set dirout2=%dirout%\surface
mkdir %dirout2%
%isosurface% -dirin %diroutdata% -saveslice %dirout2%/Slices 
if not "%ERRORLEVEL%" == "0" goto fail



:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause

