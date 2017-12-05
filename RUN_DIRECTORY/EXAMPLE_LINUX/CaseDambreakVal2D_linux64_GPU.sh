#!/bin/bash


# "name" and "dirout" are named according to the testcase

name=CaseDambreakVal2D
dirout=${name}_out


# "executables" are renamed and called from their directory
tools=../../EXECS_LINUX
gencase=${tools}/GenCase4_linux64
dualsphysics=${tools}/DualSPHysics_linux64
partvtk=${tools}/PartVTK4_linux64
partvtkout=${tools}/PartVTKOut4_linux64
measuretool=${tools}/MeasureTool4_linux64
isosurface=${tools}/IsoSurface4_linux64

# Library path must be indicated properly

current=$(pwd)
cd ${tools}
path_so=$(pwd)
cd $current
export LD_LIBRARY_PATH=$path_so

# "dirout" is created to store results or it is removed if it already exists

if [ -e $dirout ]; then
  rm -f -r $dirout
fi
mkdir $dirout


# CODES are executed according the selected parameters of execution in this testcase

errcode=0

if [ $errcode -eq 0 ]; then
  $gencase ${name}_Def $dirout/$name -save:all
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  $dualsphysics $dirout/$name $dirout -svres -gpu
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  $partvtk -dirin $dirout -filexml $dirout/${name}.xml -savevtk $dirout/PartFluid -onlytype:-all,fluid -vars:+idp,+vel,+rhop,+press,+vor
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  $partvtkout -dirin $dirout -filexml $dirout/${name}.xml -savevtk $dirout/PartFluidOut -SaveResume $dirout/ResumeFluidOut
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  $isosurface -dirin $dirout -saveslice $dirout/Slices 
  errcode=$?
fi

if [ $errcode -eq 0 ]; then
  echo All done
else
  echo Execution aborted
fi
read -n1 -r -p "Press any key to continue..." key
echo

