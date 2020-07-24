@echo off

rem "name" and "dirout" are named according to the testcase

set name=CaseOpenFtMove
set dirout=%name%_out
set diroutdata=%dirout%\data

rem "executables" are renamed and called from their directory

set dirbin=../../../bin/windows
set gencase="%dirbin%/GenCase_win64.exe"
set dualsphysicscpu="%dirbin%/DualSPHysics5.0CPU_win64.exe"
set dualsphysicsgpu="%dirbin%/DualSPHysics5.0_win64.exe"
set boundaryvtk="%dirbin%/BoundaryVTK_win64.exe"
set partvtk="%dirbin%/PartVTK_win64.exe"
set partvtkout="%dirbin%/PartVTKOut_win64.exe"
set measuretool="%dirbin%/MeasureTool_win64.exe"
set computeforces="%dirbin%/ComputeForces_win64.exe"
set isosurface="%dirbin%/IsoSurface_win64.exe"
set flowtool="%dirbin%/FlowTool_win64.exe"
set floatinginfo="%dirbin%/FloatingInfo_win64.exe"

rem RUN FIRST PART OF SIMULATION (0-1 seconds)

rem "dirout" is created to store results or it is removed if it already exists

if exist %dirout% rd /s /q %dirout%
if not "%ERRORLEVEL%" == "0" goto fail

rem CODES are executed according the selected parameters of execution in this testcase

rem Executes GenCase to create initial files for simulation.
%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes DualSPHysics to simulate first second.
%dualsphysicsgpu% -gpu %dirout%/%name% %dirout% -dirdataout data -svres -tmax:1
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes post-processing tools...
set dirout2=%dirout%\particles
%partvtk% -dirin %diroutdata% -savevtk %dirout2%/PartFluid -onlytype:-all,+fluid
if not "%ERRORLEVEL%" == "0" goto fail

%partvtkout% -dirin %diroutdata% -savevtk %dirout2%/PartFluidOut -SaveResume %dirout2%/_ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

set dirout2=%dirout%\boundary
%boundaryvtk% -loadvtk %dirout%/%name%__Actual.vtk -motiondata %diroutdata% -savevtkdata %dirout2%/Floatings.vtk
if not "%ERRORLEVEL%" == "0" goto fail

set dirout2=%dirout%\floatinginfo
%floatinginfo% -dirin %diroutdata% -savedata %dirout2%/FloatingMotion 
if not "%ERRORLEVEL%" == "0" goto fail

set dirout2=%dirout%\elevation
%measuretool% -dirin %diroutdata% -points %name%_PointsElevation.txt -onlytype:-all,+fluid -elevation -savevtk %dirout2%/PointsElevation -savecsv %dirout2%/_Elevation
if not "%ERRORLEVEL%" == "0" goto fail


rem RESTART SIMULATION AND RUN LAST PART OF SIMULATION (1-4 seconds)

set olddiroutdata=%diroutdata%
set dirout=%name%_restart_out
set diroutdata=%dirout%\data

rem "redirout" is created to store results of restart simulation

if exist %dirout% rd /s /q %dirout%
if not "%ERRORLEVEL%" == "0" goto fail

rem CODES are executed according the selected parameters of execution in this testcase

rem Executes GenCase to create initial files for simulation.
%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes DualSPHysics to simulate the last 3 seconds.
%dualsphysicsgpu% -gpu %dirout%/%name% %dirout% -dirdataout data -svres -partbegin:100 %olddiroutdata%
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes post-processing tools for restart simulation...
set dirout2=%dirout%\particles
%partvtk% -dirin %diroutdata% -savevtk %dirout2%/PartFluid -onlytype:-all,+fluid
if not "%ERRORLEVEL%" == "0" goto fail

%partvtkout% -dirin %diroutdata% -savevtk %dirout2%/PartFluidOut -SaveResume %dirout2%/ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

set dirout2=%dirout%\boundary
%boundaryvtk% -loadvtk %dirout%/%name%__Actual.vtk -motiondata0 %olddiroutdata% -motiondata %diroutdata% -savevtkdata %dirout2%/Floatings.vtk
if not "%ERRORLEVEL%" == "0" goto fail

set dirout2=%dirout%\fluidslices
%isosurface% -dirin %diroutdata% -saveslice %dirout2%/Slices
if not "%ERRORLEVEL%" == "0" goto fail

set dirout2=%dirout%\floatinginfo
%floatinginfo% -dirin %diroutdata%  -savedata %dirout2%/FloatingMotion 
if not "%ERRORLEVEL%" == "0" goto fail

set dirout2=%dirout%\elevation
%measuretool% -dirin %diroutdata% -points %name%_PointsElevation.txt -onlytype:-all,+fluid -elevation -savevtk %dirout2%/PointsElevation -savecsv %dirout2%/_Elevation
if not "%ERRORLEVEL%" == "0" goto fail


:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause
