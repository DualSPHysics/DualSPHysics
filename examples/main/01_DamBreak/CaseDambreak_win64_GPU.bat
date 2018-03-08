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
set measureboxes="%dirbin%/MeasureBoxes4_win64.exe"

REM "dirout" is created to store results or it is removed if it already exists

if exist %dirout% del /Q %dirout%\*.*
if not exist %dirout% mkdir %dirout%

REM CODES are executed according the selected parameters of execution in this testcase

%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

%dualsphysicsgpu% -gpu %dirout%/%name% %dirout% -svres
if not "%ERRORLEVEL%" == "0" goto fail

%partvtk% -dirin %dirout% -savevtk %dirout%/PartFluid -onlytype:-all,+fluid
if not "%ERRORLEVEL%" == "0" goto fail

%partvtkout% -dirin %dirout% -filexml %dirout%/%name%.xml -savevtk %dirout%/PartFluidOut -SaveResume %dirout%/ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

%measuretool% -dirin %dirout% -points CaseDambreak_PointsVelocity.txt -onlytype:-all,+fluid -vars:-all,+vel.x,+vel.m -savevtk %dirout%/PointsVelocity -savecsv %dirout%/PointsVelocity
if not "%ERRORLEVEL%" == "0" goto fail

%measuretool% -dirin %dirout% -points CaseDambreak_PointsPressure_Incorrect.txt -onlytype:-all,+fluid -vars:-all,+press,+kcorr -kcusedummy:0 -kclimit:0.5 -savevtk %dirout%/PointsPressure_Incorrect -savecsv %dirout%/PointsPressure_Incorrect
if not "%ERRORLEVEL%" == "0" goto fail

%measuretool% -dirin %dirout% -points CaseDambreak_PointsPressure_Correct.txt -onlytype:-all,+fluid -vars:-all,+press,+kcorr -kcusedummy:0 -kclimit:0.5 -savevtk %dirout%/PointsPressure_Correct -savecsv %dirout%/PointsPressure_Correct
if not "%ERRORLEVEL%" == "0" goto fail

%computeforces% -dirin %dirout% -filexml %dirout%/%name%.xml -onlymk:21 -savecsv %dirout%/WallForce 
if not "%ERRORLEVEL%" == "0" goto fail

set planesy="-slicevec:0:0.1:0:0:1:0 -slicevec:0:0.2:0:0:1:0 -slicevec:0:0.3:0:0:1:0 -slicevec:0:0.4:0:0:1:0 -slicevec:0:0.5:0:0:1:0 -slicevec:0:0.6:0:0:1:0"
set planesx="-slicevec:0.1:0:0:1:0:0 -slicevec:0.2:0:0:1:0:0 -slicevec:0.3:0:0:1:0:0 -slicevec:0.4:0:0:1:0:0 -slicevec:0.5:0:0:1:0:0 -slicevec:0.6:0:0:1:0:0 -slicevec:0.7:0:0:1:0:0 -slicevec:0.8:0:0:1:0:0 -slicevec:0.9:0:0:1:0:0 -slicevec:1.0:0:0:1:0:0"
set planesd="-slice3pt:0:0:0:1:0.7:0:1:0.7:1"
%isosurface% -dirin %dirout% -saveiso %dirout%/Surface -vars:-all,vel,rhop,idp,type -saveslice %dirout%/Slices %planesy% %planesx% %planesd%
if not "%ERRORLEVEL%" == "0" goto fail

%measureboxes% -dirin %dirout% -fileboxes CaseDambreak_FileBoxes.txt -savecsv %dirout%/ResultBoxes.csv -savevtk %dirout%/Boxes.vtk
if not "%ERRORLEVEL%" == "0" goto fail



:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause
