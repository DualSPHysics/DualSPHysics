@echo off

REM "name" and "dirout" are named according to the testcase

set name=CaseDambreak
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
%dualsphysicscpu% -cpu %dirout%/%name% %dirout% -dirdataout data -svres
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes PartVTK4 to create VTK files with particles.
set dirout2=%dirout%\particles
mkdir %dirout2%
%partvtk% -dirin %diroutdata% -savevtk %dirout2%/PartFluid -onlytype:-all,+fluid
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes PartVTKOut4 to create VTK files with excluded particles.
%partvtkout% -dirin %diroutdata% -filexml %dirout%/%name%.xml -savevtk %dirout2%/PartFluidOut -SaveResume %dirout%/ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes MeasureTool4 to create VTK files with velocity and a CSV file with velocity at each simulation time.
set dirout2=%dirout%\velocity
mkdir %dirout2%
%measuretool% -dirin %diroutdata% -points CaseDambreak_PointsVelocity.txt -onlytype:-all,+fluid -vars:-all,+vel.x,+vel.m -savevtk %dirout2%/PointsVelocity -savecsv %dirout%/PointsVelocity
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes MeasureTool4 to create VTK files with incorrect pressure and a CSV file with value at each simulation time.
set dirout2=%dirout%\pressure
mkdir %dirout2%
%measuretool% -dirin %diroutdata% -points CaseDambreak_PointsPressure_Incorrect.txt -onlytype:-all,+fluid -vars:-all,+press,+kcorr -kcusedummy:0 -kclimit:0.5 -savevtk %dirout2%/PointsPressure_Incorrect -savecsv %dirout%/PointsPressure_Incorrect
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes MeasureTool4 to create VTK files with correct pressure and a CSV file with value at each simulation time.
%measuretool% -dirin %diroutdata% -points CaseDambreak_PointsPressure_Correct.txt -onlytype:-all,+fluid -vars:-all,+press,+kcorr -kcusedummy:0 -kclimit:0.5 -savevtk %dirout2%/PointsPressure_Correct -savecsv %dirout%/PointsPressure_Correct
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes ComputeForces to create a CSV file with force at each simulation time.
%computeforces% -dirin %diroutdata% -filexml %dirout%/%name%.xml -onlymk:21 -savecsv %dirout%/WallForce 
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes IsoSurface4 to create VTK files with surface fluid and slices of surface.
set dirout2=%dirout%\surface
mkdir %dirout2%
set planesy="-slicevec:0:0.1:0:0:1:0 -slicevec:0:0.2:0:0:1:0 -slicevec:0:0.3:0:0:1:0 -slicevec:0:0.4:0:0:1:0 -slicevec:0:0.5:0:0:1:0 -slicevec:0:0.6:0:0:1:0"
set planesx="-slicevec:0.1:0:0:1:0:0 -slicevec:0.2:0:0:1:0:0 -slicevec:0.3:0:0:1:0:0 -slicevec:0.4:0:0:1:0:0 -slicevec:0.5:0:0:1:0:0 -slicevec:0.6:0:0:1:0:0 -slicevec:0.7:0:0:1:0:0 -slicevec:0.8:0:0:1:0:0 -slicevec:0.9:0:0:1:0:0 -slicevec:1.0:0:0:1:0:0"
set planesd="-slice3pt:0:0:0:1:0.7:0:1:0.7:1"
%isosurface% -dirin %diroutdata% -saveiso %dirout2%/Surface -vars:-all,vel,rhop,idp,type -saveslice %dirout2%/Slices %planesy% %planesx% %planesd%
if not "%ERRORLEVEL%" == "0" goto fail

REM Executes FlowTool4 to create VTK files with particles assigned to different zones and a CSV file with information of each zone.
set dirout2=%dirout%\flow
mkdir %dirout2%
%flowtool% -dirin %diroutdata% -fileboxes CaseDambreak_FileBoxes.txt -savecsv %dirout%/ResultBoxes.csv -savevtk %dirout2%/Boxes.vtk
if not "%ERRORLEVEL%" == "0" goto fail



:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause
