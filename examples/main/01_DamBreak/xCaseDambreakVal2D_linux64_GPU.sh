#!/bin/bash

# "name" and "dirout" are named according to the testcase

name=CaseDambreakVal2D
dirout=${name}_out
diroutdata=${dirout}/data

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

# "dirout" to store results is removed if it already exists
if [ -e $dirout ]; then
  rm -r $dirout
fi


# CODES are executed according the selected parameters of execution in this testcase
errcode=0

# Executes GenCase4 to create initial files for simulation.
if [ $errcode -eq 0 ]; then
  $gencase ${name}_Def $dirout/$name -save:all
  errcode=$?
fi

# Executes DualSPHysics to simulate SPH method.
if [ $errcode -eq 0 ]; then
  $dualsphysicsgpu -gpu $dirout/$name $dirout -dirdataout data -svres
  errcode=$?
fi

# Executes PartVTK4 to create VTK files with particles.
dirout2=${dirout}/particles
if [ $errcode -eq 0 ]; then
  $partvtk -dirin $diroutdata -savevtk $dirout2/PartFluid -onlytype:-all,fluid -vars:+idp,+vel,+rhop,+press,+vor
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  $partvtk -dirin $diroutdata -savevtk $dirout2/PartBound -onlytype:-all,bound -vars:-all -last:0
  errcode=$?
fi

# Executes PartVTKOut4 to create VTK files with excluded particles.
if [ $errcode -eq 0 ]; then
  $partvtkout -dirin $diroutdata -savevtk $dirout2/PartFluidOut -SaveResume $dirout2/_ResumeFluidOut
  errcode=$?
fi

# Executes IsoSurface4 to create VTK files with slices of surface.
dirout2=${dirout}/surface
if [ $errcode -eq 0 ]; then
  $isosurface -dirin $diroutdata -saveslice $dirout2/Slices 
  errcode=$?
fi


if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo
