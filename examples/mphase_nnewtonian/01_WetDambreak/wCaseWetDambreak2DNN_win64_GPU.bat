@echo off
setlocal EnableDelayedExpansion
rem Don't remove the two jump line after than the next line [set NL=^]
set NL=^


rem "name" and "dirout" are named according to the testcase

set name=CaseWetDambreak2DNN
set dirout=%name%_out
set diroutdata=%dirout%\data

rem "executables" are renamed and called from their directory

set dirbin=../../../bin/windows
set dirbindsnn=%dirbin%/DSNNewtonian
set gencase="%dirbin%/GenCase_win64.exe"
set dualsphysicscpu="%dirbindsnn%/DualSPHysics5.0_NNewtonianCPU_win64.exe"
set dualsphysicsgpu="%dirbindsnn%/DualSPHysics5.0_NNewtonian_win64.exe"
set boundaryvtk="%dirbin%/BoundaryVTK_win64.exe"
set partvtk="%dirbin%/PartVTK_win64.exe"
set partvtkout="%dirbin%/PartVTKOut_win64.exe"
set measuretool="%dirbin%/MeasureTool_win64.exe"
set computeforces="%dirbin%/ComputeForces_win64.exe"
set isosurface="%dirbin%/IsoSurface_win64.exe"
set flowtool="%dirbin%/FlowTool_win64.exe"
set floatinginfo="%dirbin%/FloatingInfo_win64.exe"

:menu
if exist %dirout% ( 
	set /p option="The folder "%dirout%" already exists. Choose an option.!NL!  [1]- Delete it and continue.!NL!  [2]- Execute post-processing.!NL!  [3]- Abort and exit.!NL!"
	if "!option!" == "1" goto run else (
		if "!option!" == "2" goto postprocessing else (
			if "!option!" == "3" goto fail else ( 
				goto menu
			)
		)
	)
)

:run
rem "dirout" to store results is removed if it already exists
if exist %dirout% rd /s /q %dirout%

rem CODES are executed according the selected parameters of execution in this testcase

%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

%dualsphysicsgpu% -gpu %dirout%/%name% %dirout% -dirdataout data -svres
if not "%ERRORLEVEL%" == "0" goto fail

:postprocessing
rem Executes PartVTK to create VTK files with particles.
set dirout2=%dirout%\particles
%partvtk% -dirin %diroutdata% -savevtk %dirout2%/PartFluid -onlytype:-all,fluid -vars:+idp,+vel,+rhop,+press,+vor,+mk
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes PartVTKOut to create VTK files with excluded particles.
%partvtkout% -dirin %diroutdata% -savevtk %dirout2%/PartFluidOut -SaveResume %dirout2%/_ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes IsoSurface to create VTK files with slices of surface.
set dirout2=%dirout%\surface
%isosurface% -dirin %diroutdata% -saveslice %dirout2%/Slices_phase1 -onlymk:1 
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes IsoSurface to create VTK files with slices of surface.
set dirout2=%dirout%\surface
%isosurface% -dirin %diroutdata% -saveslice %dirout2%/Slices_phase2 -onlymk:2 
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes IsoSurface to create VTK files with slices of surface.
set dirout2=%dirout%\surface
%isosurface% -dirin %diroutdata% -saveslice %dirout2%/Slices_phase3 -onlymk:3 
if not "%ERRORLEVEL%" == "0" goto fail

:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause

