#!/bin/bash 

fail () { 
 echo Execution aborted. 
 read -n1 -r -p "Press any key to continue..." key 
 exit 1 
}

# "name" and "dirout" are named according to the testcase

export name=CaseDambreak
export dirout=${name}_out
export diroutdata=${dirout}/data

# "executables" are renamed and called from their directory

export dirbin=../../../bin/linux
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${dirbin}
export gencase="${dirbin}/GenCase_linux64"
export dualsphysicscpu="${dirbin}/DualSPHysics5.0CPU_linux64"
export dualsphysicsgpu="${dirbin}/DualSPHysics5.0_linux64"
export boundaryvtk="${dirbin}/BoundaryVTK_linux64"
export partvtk="${dirbin}/PartVTK_linux64"
export partvtkout="${dirbin}/PartVTKOut_linux64"
export measuretool="${dirbin}/MeasureTool_linux64"
export computeforces="${dirbin}/ComputeForces_linux64"
export isosurface="${dirbin}/IsoSurface_linux64"
export flowtool="${dirbin}/FlowTool_linux64"
export floatinginfo="${dirbin}/FloatingInfo_linux64"

option=-1
 if [ -e $dirout ]; then
 while [ "$option" != 1 -a "$option" != 2 -a "$option" != 3 ] 
 do 

	echo -e "The folder "${dirout}" already exists. Choose an option.
  [1]- Delete it and continue.
  [2]- Execute post-processing.
  [3]- Abort and exit.
"
 read -n 1 option 
 done 
  else 
   option=1 
fi 

if [ $option -eq 1 ]; then
# "dirout" to store results is removed if it already exists
if [ -e ${dirout} ]; then rm -r ${dirout}; fi

# CODES are executed according the selected parameters of execution in this testcase

# Executes GenCase to create initial files for simulation.
${gencase} ${name}_Def ${dirout}/${name} -save:all
if [ $? -ne 0 ] ; then fail; fi

# Executes DualSPHysics to simulate SPH method.
${dualsphysicscpu} ${dirout}/${name} ${dirout} -dirdataout data -svres
if [ $? -ne 0 ] ; then fail; fi

fi

if [ $option -eq 2 -o $option -eq 1 ]; then
# Executes PartVTK to create VTK files with particles.
export dirout2=${dirout}/particles
${partvtk} -dirin ${diroutdata} -savevtk ${dirout2}/PartFluid -onlytype:-all,+fluid
if [ $? -ne 0 ] ; then fail; fi

# Executes PartVTKOut to create VTK files with excluded particles.
${partvtkout} -dirin ${diroutdata} -savevtk ${dirout2}/PartFluidOut -SaveResume ${dirout2}/_ResumeFluidOut
if [ $? -ne 0 ] ; then fail; fi

# Executes MeasureTool to create VTK files with velocity and a CSV file with velocity at each simulation time.
export dirout2=${dirout}/measuretool
${measuretool} -dirin ${diroutdata} -points CaseDambreak_PointsVelocity.txt -onlytype:-all,+fluid -vars:-all,+vel.x,+vel.m -savevtk ${dirout2}/PointsVelocity -savecsv ${dirout2}/_PointsVelocity
if [ $? -ne 0 ] ; then fail; fi

# Executes MeasureTool to create VTK files with incorrect pressure and a CSV file with value at each simulation time.
${measuretool} -dirin ${diroutdata} -points CaseDambreak_PointsPressure_Incorrect.txt -onlytype:-all,+fluid -vars:-all,+press,+kcorr -kcusedummy:0 -kclimit:0.5 -savevtk ${dirout2}/PointsPressure_Incorrect -savecsv ${dirout2}/_PointsPressure_Incorrect
if [ $? -ne 0 ] ; then fail; fi

# Executes MeasureTool to create VTK files with correct pressure and a CSV file with value at each simulation time.
${measuretool} -dirin ${diroutdata} -points CaseDambreak_PointsPressure_Correct.txt -onlytype:-all,+fluid -vars:-all,+press,+kcorr -kcusedummy:0 -kclimit:0.5 -savevtk ${dirout2}/PointsPressure_Correct -savecsv ${dirout2}/_PointsPressure_Correct
if [ $? -ne 0 ] ; then fail; fi

# Executes ComputeForces to create a CSV file with force at each simulation time.
export dirout2=${dirout}/forces
${computeforces} -dirin ${diroutdata} -onlymk:20 -viscoart:0.1 -savecsv ${dirout2}/_ForceBuilding
if [ $? -ne 0 ] ; then fail; fi

# Executes IsoSurface to create VTK files with surface fluid and slices of surface.
export dirout2=${dirout}/surface
export planesy="-slicevec:0:0.1:0:0:1:0 -slicevec:0:0.2:0:0:1:0 -slicevec:0:0.3:0:0:1:0 -slicevec:0:0.4:0:0:1:0 -slicevec:0:0.5:0:0:1:0 -slicevec:0:0.6:0:0:1:0"
export planesx="-slicevec:0.1:0:0:1:0:0 -slicevec:0.2:0:0:1:0:0 -slicevec:0.3:0:0:1:0:0 -slicevec:0.4:0:0:1:0:0 -slicevec:0.5:0:0:1:0:0 -slicevec:0.6:0:0:1:0:0 -slicevec:0.7:0:0:1:0:0 -slicevec:0.8:0:0:1:0:0 -slicevec:0.9:0:0:1:0:0 -slicevec:1.0:0:0:1:0:0"
export planesd="-slice3pt:0:0:0:1:0.7:0:1:0.7:1"
${isosurface} -dirin ${diroutdata} -saveiso ${dirout2}/Surface -vars:-all,vel,rhop,idp,type -saveslice ${dirout2}/Slices ${planesy} ${planesx} ${planesd}
if [ $? -ne 0 ] ; then fail; fi

# Executes FlowTool to create VTK files with particles assigned to different zones and a CSV file with information of each zone.
export dirout2=${dirout}/flow
${flowtool} -dirin ${diroutdata} -fileboxes CaseDambreak_FileBoxes.txt -savecsv ${dirout2}/_ResultFlow.csv -savevtk ${dirout2}/Boxes.vtk
if [ $? -ne 0 ] ; then fail; fi

fi
if [ $option != 3 ];then
 echo All done
 else
 echo Execution aborted
fi

read -n1 -r -p "Press any key to continue..." key
