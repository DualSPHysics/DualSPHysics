#!/bin/bash

# "name" and "dirout" are named according to the testcase

name=CaseDambreakVal2D
dirout=${name}_out

# "executables" are renamed and called from their directory

dirbin=../../../bin/linux
gencase="${dirbin}/GenCase4_linux64"
dualsphysicscpu="${dirbin}/DualSPHysics4.2CPU_linux64"
dualsphysicsgpu="${dirbin}/DualSPHysics4.2_linux64"
boundaryvtk="${dirbin}/BoundaryVTK4_linux64"
partvtk="${dirbin}/PartVTK4_linux64"
partvtkout="${dirbin}/PartVTKOut4_linux64"
measuretool="${dirbin}/MeasureTool4_linux64"
computeforces="${dirbin}/ComputeForces4_linux64"
isosurface="${dirbin}/IsoSurface4_linux64"
measureboxes="${dirbin}/MeasureBoxes4_linux64"


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

# Executes DualSPHysics to simulate SPH method.
if [ $errcode -eq 0 ]; then
  $dualsphysicscpu -cpu $dirout/$name $dirout -dirdataout data -svres -tmax:0.0005 -tout:0.0001
  errcode=$?
fi

# Executes PartVTK4 to create VTK files with particles.
dirout2=${dirout}/particles; mkdir $dirout2
if [ $errcode -eq 0 ]; then
  $partvtk -dirin $diroutdata -filexml $dirout/${name}.xml -savevtk $dirout2/PartFluid -onlytype:-all,fluid -vars:+idp,+vel,+rhop,+press,+vor
  errcode=$?
fi

# Executes PartVTKOut4 to create VTK files with excluded particles.
if [ $errcode -eq 0 ]; then
  $partvtkout -dirin $diroutdata -filexml $dirout/${name}.xml -savevtk $dirout2/PartFluidOut -SaveResume $dirout/ResumeFluidOut
  errcode=$?
fi

# Executes IsoSurface4 to create VTK files with slices of surface.
dirout2=${dirout}/surface; mkdir $dirout2
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
