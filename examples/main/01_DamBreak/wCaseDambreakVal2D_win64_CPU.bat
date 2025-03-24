@echo off
setlocal EnableDelayedExpansion
rem Don't remove the two jump line after than the next line [set NL=^]
set NL=^


rem "name" and "dirout" are named according to the testcase

set name=CaseDambreakVal2D
set dirout=%name%_out
set diroutdata=%dirout%\data

rem "executables" are renamed and called from their directory

set dirbin=../../../bin/windows
set gencase="%dirbin%/GenCase_win64.exe"
set dualsphysicscpu="%dirbin%/DualSPHysics5.4CPU_win64.exe"
set dualsphysicsgpu="%dirbin%/DualSPHysics5.4_win64.exe"
set boundaryvtk="%dirbin%/BoundaryVTK_win64.exe"
set partvtk="%dirbin%/PartVTK_win64.exe"
set partvtkout="%dirbin%/PartVTKOut_win64.exe"
set measuretool="%dirbin%/MeasureTool_win64.exe"
set computeforces="%dirbin%/ComputeForces_win64.exe"
set isosurface="%dirbin%/IsoSurface_win64.exe"
set flowtool="%dirbin%/FlowTool_win64.exe"
set floatinginfo="%dirbin%/FloatingInfo_win64.exe"
set tracerparts="%dirbin%/TracerParts_win64.exe"

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

rem Executes GenCase to create initial files for simulation.
%gencase% %name%_Def %dirout%/%name% -save:all
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes DualSPHysics to simulate SPH method.
%dualsphysicscpu% %dirout%/%name% %dirout%
if not "%ERRORLEVEL%" == "0" goto fail

:postprocessing
rem Executes PartVTK to create VTK files with particles.
set dirout2=%dirout%\particles
%partvtk% -dirdata %diroutdata% -savevtk %dirout2%/PartFluid -onlytype:-all,fluid -vars:+idp,+vel,+rhop,+press,+vor,+energy
if not "%ERRORLEVEL%" == "0" goto fail

%partvtk% -dirdata %diroutdata% -savevtk %dirout2%/PartBound -onlytype:-all,bound -vars:-all -last:0
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes PartVTKOut to create VTK files with excluded particles.
%partvtkout% -dirdata %diroutdata% -savevtk %dirout2%/PartFluidOut -SaveResume %dirout2%/_ResumeFluidOut
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes PartVTK to create CSV with energy values.
%partvtk% -dirdata %diroutdata% -saveenergy %dirout%/Energy -last:90
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes TracerParts to create VTK files with trajectory of some fluid particles.
set dirout2=%dirout%\tracer
%tracerparts% -dirdata %diroutdata% -savevtk %dirout2%/BorderParts -onlytype:-all,+fluid -nearpartsdist:0.02 -nearpartsdef:pt=0.1:0:0.1,pt=0.15:0:0.15,pt=0.2:0:0.2,ptels[x=1:0:1,z=0:0.05:2],ptels[x=0:0.05:1,z=2:0:2] -tailsize:200
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes MeasureTool to create VTK and CSV files with elevation at each simulation time.
set dirout2=%dirout%\measuretool
%measuretool% -dirdata %diroutdata% -pointsdef:ptels[x=0.2:0:0.2,y=0:0:0,z=0:0.02:2.1] -onlytype:-all,+fluid -elevation -savevtk %dirout2%/EtaPoints -savecsv %dirout%/MeasuredA
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes MeasureTool to create VTK and CSV files with velocity and pressure at 2 points.
%measuretool% -dirdata %diroutdata% -pointsdef:pt=0.2:0:0.2,pt=0.4:0:0.2 -onlytype:-all,+fluid -vars:-all,vel,press -savevtk %dirout2%/CheckPoints -savecsv %dirout%/MeasuredB
if not "%ERRORLEVEL%" == "0" goto fail

rem Executes IsoSurface to create VTK files with slices of surface.
set dirout2=%dirout%\surface
%isosurface% -dirdata %diroutdata% -saveslice %dirout2%/Slices 
if not "%ERRORLEVEL%" == "0" goto fail


:success
echo All done
goto end

:fail
echo Execution aborted.

:end
pause

