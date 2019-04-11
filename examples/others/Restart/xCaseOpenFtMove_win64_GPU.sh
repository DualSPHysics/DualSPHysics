#!/bin/bash

# "name" and "dirout" are named according to the testcase

name=CaseOpenFtMove
dirout=${name}_out

# "executables" are renamed and called from their directory

dirbin=../../../bin/linux
gencase="${dirbin}/GenCase4_linux64"
dualsphysicscpu="${dirbin}/DualSPHysics4.4CPU_linux64"
dualsphysicsgpu="${dirbin}/DualSPHysics4.4_linux64"
boundaryvtk="${dirbin}/BoundaryVTK4_linux64"
partvtk="${dirbin}/PartVTK4_linux64"
partvtkout="${dirbin}/PartVTKOut4_linux64"
measuretool="${dirbin}/MeasureTool4_linux64"
computeforces="${dirbin}/ComputeForces4_linux64"
isosurface="${dirbin}/IsoSurface4_linux64"
flowtool="${dirbin}/FlowTool4_linux64"
floatinginfo="${dirbin}/FloatingInfo4_linux64"

# Library path must be indicated properly

current=$(pwd)
cd $dirbin
path_so=$(pwd)
cd $current
export LD_LIBRARY_PATH=$path_so


# "dirout" is created to store results or it is cleaned if it already exists

if [ -e $dirout ]; then
  rm -r $dirout
fi
mkdir $dirout
diroutdata=${dirout}/data; mkdir $diroutdata

# CODES are executed according the selected parameters of execution in this testcase
errcode=0

# Executes GenCase4 to create initial files for simulation.
if [ $errcode -eq 0 ]; then
  $gencase ${name}_Def $dirout/$name -save:all
  errcode=$?
fi

# Executes DualSPHysics to simulate first second.
if [ $errcode -eq 0 ]; then
  $dualsphysicsgpu -gpu $dirout/$name $dirout -dirdataout data -svres -tmax:1
  errcode=$?
fi

# Executes post-processing tools...
dirout2=${dirout}/particles; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $partvtk -dirin $diroutdata -savevtk $dirout2/PartFluid -onlytype:-all,+fluid
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  $partvtkout -dirin $diroutdata -savevtk $dirout2/PartFluidOut -SaveResume $dirout2/_ResumeFluidOut
  errcode=$?
fi

dirout2=${dirout}/floatings; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $boundaryvtk -loadvtk $dirout/${name}__Actual.vtk -motiondata $diroutdata -savevtkdata $dirout2/Floatings.vtk
  errcode=$?
fi

dirout2=${dirout}/fluidslices; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $isosurface -dirin $diroutdata -saveslice $dirout2/Slices
  errcode=$?
fi

dirout2=${dirout}/floatinginfo; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $floatinginfo -dirin $diroutdata -savemotion -savedata $dirout2/FloatingMotion
  errcode=$?
fi

dirout2=${dirout}/height; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $measuretool -dirin $diroutdata -points ${name}_PointsHeights.txt -onlytype:-all,+fluid -height -savevtk $dirout2/PointsHeight -savecsv $dirout2/_Height
  errcode=$?
fi

# RESTART SIMULATION AND RUN LAST PART OF SIMULATION (1-4 seconds)

olddiroutdata=$diroutdata
dirout=${name}_restart_out
diroutdata=$dirout/data

# "redirout" is created to store results of restart simulation

if [ -e $dirout ]; then
  rm -r $dirout
fi
mkdir $dirout
diroutdata=${dirout}/data; mkdir $diroutdata

# CODES are executed according the selected parameters of execution in this testcase

# Executes GenCase4 to create initial files for simulation.
if [ $errcode -eq 0 ]; then
  $gencase ${name}_Def $dirout/$name -save:all
  errcode=$?
fi

# Executes DualSPHysics to simulate the last 3 seconds.
if [ $errcode -eq 0 ]; then
  $dualsphysicsgpu -gpu $dirout/$name $dirout -dirdataout data -svres -partbegin:100 $olddiroutdata
  errcode=$?
fi

# Executes post-processing tools for restart simulation...
dirout2=${dirout}/particles; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $partvtk -dirin $diroutdata -savevtk $dirout2/PartFluid -onlytype:-all,+fluid
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  $partvtkout -dirin $diroutdata -savevtk $dirout2/PartFluidOut -SaveResume $dirout2/ResumeFluidOut
  errcode=$?
fi

dirout2=${dirout}/floatings; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $boundaryvtk -loadvtk $dirout/${name}__Actual.vtk -motiondata0 $olddiroutdata -motiondata $diroutdata -savevtkdata $dirout2/Floatings.vtk
  errcode=$?
fi

dirout2=${dirout}/fluidslices; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $isosurface -dirin $diroutdata -saveslice $dirout2/Slices
  errcode=$?
fi

dirout2=${dirout}/floatinginfo; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $floatinginfo -dirin $diroutdata -savemotion -savedata $dirout2/FloatingMotion
  errcode=$?
fi

dirout2=${dirout}/height; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $measuretool -dirin $diroutdata -points ${name}_PointsHeights.txt -onlytype:-all,+fluid -height -savevtk $dirout2/PointsHeight -savecsv $dirout2/_Height
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
