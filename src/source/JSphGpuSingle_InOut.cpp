//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpuSingle_InOut.cpp  \brief Implements InOut functions of class \ref JSphGpuSingle.

#include "JSphGpuSingle.h"
#include "JCellDivGpuSingle.h"
#include "JArraysGpu.h"
#include "JSphMk.h"
#include "JDsPartsInit.h"
#include "JSphInOut.h"
#include "JSphBoundCorr.h"
#include "JSphGpu_InOut_iker.h"
#include "JSphInOutPoints.h"
#include "JSimpleNeigs.h"
#include "FunctionsMath.h"
#include "FunctionsCuda.h"
#include "JDebugSphGpu.h"
#include "JAppInfo.h"
#include "JTimeControl.h"
#include "JDataArrays.h"
#include "JVtkLib.h"
#include <climits>

using namespace std;

//==============================================================================
/// Mark special fluid particles to ignore.
/// Marca las particulas fluidas especiales para ignorar.
//==============================================================================
void JSphGpuSingle::InOutIgnoreFluidDef(const std::vector<unsigned> &mkfluidlist){
  const unsigned nc=unsigned(mkfluidlist.size());
  for(unsigned c=0;c<nc;c++){
    const unsigned cmk=MkInfo->GetMkBlockByMkFluid(mkfluidlist[c]);
    if(cmk<MkInfo->Size()){
      const typecode rcode=MkInfo->Mkblock(cmk)->Code;
      const typecode rcode2=CODE_SetOutIgnore(rcode);
      //-Mark special fluid particles to ignore.
      cusphinout::InOutIgnoreFluidDef(Np,rcode,rcode2,Codeg);
    }
  }
}

//==============================================================================
/// Checks proximity of inout particles to other particles and excludes fluid 
/// particles near the inout particles.
///
/// Comprueba proximidad de particulas inout con otras particulas y excluye 
/// particulas fluid cerca de particulas inout.
//==============================================================================
void JSphGpuSingle::InOutCheckProximity(unsigned newnp){
  if(Np && newnp){
    //-Obtain particle data from GPU memory.
    if(Np!=ParticlesDataDown(Np,0,true,false))Run_Exceptioon("The number of particles is invalid.");
    //-Check proximity on CPU.
    InOut->InitCheckProximity(Np,newnp,Scell,AuxPos,Idp,Code);
    //-Uploads updated code of particles to the GPU.
    cudaMemcpy(Codeg,Code,sizeof(typecode)*Np,cudaMemcpyHostToDevice);
  }
}

//==============================================================================
/// Initialises inlet/outlet conditions.
/// Inicia condiciones inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutInit(double timestepini){
  InOut->Nstep=Nstep; //-For debug.
  TmgStart(Timers,TMG_SuInOut);
  Log->Print("Initialising InOut...");
  if(PartBegin)Run_Exceptioon("Simulation restart not allowed when Inlet/Outlet is used.");

  //-Configures InOut zones and prepares new inout particles to create.
  const unsigned newnp=InOut->Config(timestepini,Stable,PeriActive
    ,MapRealPosMin,MapRealPosMax,MkInfo->GetCodeNewFluid()
    ,PartsInit,GaugeSystem,NuxLib);

  //-Mark special fluid particles to ignore. | Marca las particulas fluidas especiales para ignorar.
  InOutIgnoreFluidDef(InOut->MkFluidList);

  //Log->Printf("++> newnp:%u",newnp);
  //-Resizes memory when it is necessary (always at the beginning).
  if(true || !CheckGpuParticlesSize(Np+newnp)){
    const unsigned newnp2=newnp+InOut->GetNpResizePlus0();
    TmgStop(Timers,TMG_SuInOut);
    ResizeParticlesSize(Np+newnp2,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    TmgStart(Timers,TMG_SuInOut);
  }

  //-Creates initial inlet particles.
  unsigned idnext=IdMax+1;
  //-Compute initial data on CPU memory.
  InOut->LoadInitPartsData(idnext,newnp,Idp,Code,AuxPos,Velrhop); 
  for(unsigned p=0;p<newnp;p++){
    Posxy[p]=TDouble2(AuxPos[p].x,AuxPos[p].y);
    Posz [p]=AuxPos[p].z;
  }
  //-Copy data to GPU memory.
  cudaMemcpy(Idpg    +Np,Idp    ,sizeof(unsigned)*newnp,cudaMemcpyHostToDevice);
  cudaMemcpy(Codeg   +Np,Code   ,sizeof(typecode)*newnp,cudaMemcpyHostToDevice);
  cudaMemcpy(Posxyg  +Np,Posxy  ,sizeof(double2) *newnp,cudaMemcpyHostToDevice);
  cudaMemcpy(Poszg   +Np,Posz   ,sizeof(double)  *newnp,cudaMemcpyHostToDevice);
  cudaMemcpy(Velrhopg+Np,Velrhop,sizeof(float4)  *newnp,cudaMemcpyHostToDevice);

  //-Checks position of new particles and calculates cell.
  cusphinout::UpdatePosFluid(PeriActive,newnp,Np,Posxyg,Poszg,Dcellg,Codeg);

  //-Updates new particle values for Laminar+SPS.
  if(SpsTaug)cudaMemset(SpsTaug+Np,0,sizeof(tsymatrix3f)*newnp);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",0,Np,Np+newnp,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates number of particles.
  Np+=newnp;
  TotalNp+=newnp;
  InOut->AddNewNp(newnp);
  IdMax=unsigned(TotalNp-1);

  //-Checks proximity of inout particles to other particles and excludes fluid particles near the inout particles.
  InOutCheckProximity(newnp);

  //-Shows configuration.
  InOut->VisuConfig("\nInOut configuration:"," ");
  //-Checks invalid options for symmetry. //<vs_syymmetry>
  if(Symmetry && InOut->Use_ExtrapolatedData())Run_Exceptioon("Symmetry is not allowed with inlet/outlet conditions when extrapolate option is enabled."); //<vs_syymmetry>

  //-Updates divide information.
  TmgStop(Timers,TMG_SuInOut);
  RunCellDivide(true);
  TmgStart(Timers,TMG_SuInOut);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",1,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates Velocity data of inout zones according to current timestep.
  InOut->UpdateVelData(timestepini);
  //-Updates Zsurf data of inout zones according to current timestep.
  InOut->UpdateZsurfData(timestepini,true);

  //-Updates inout particle data according inlet configuration.
  InOutUpdatePartsData(timestepini);

  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",2,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  TmgStop(Timers,TMG_SuInOut);
}

//==============================================================================
/// ComputeStep over inlet/outlet particles:
/// - Move particles in/out according its velocity.
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new in/out particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
void JSphGpuSingle::InOutComputeStep(double stepdt){
  const double newtimestep=TimeStep+stepdt;
  InOut->Nstep=Nstep; //-For debug.
  //Log->Printf("%u>--------> [InOutComputeStep_000]",Nstep);
  //DgSaveVtkParticlesGpu("BB_ComputeStepA.vtk",DgNum,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  TmgStart(Timers,TMG_SuInOut);
  //-Resizes memory when it is necessary. InOutCount is the maximum number of new inlet particles.
  if(!CheckGpuParticlesSize(Np+InOut->GetCurrentNp())){
    if(!InOut->GetNpResizePlus1())Run_Exceptioon("Allocated memory is not enough and resizing is not allowed by XML configuration (check the value inout.memoryresize.size).");
    const unsigned newnp2=InOut->GetCurrentNp()+InOut->GetNpResizePlus1();
    TmgStop(Timers,TMG_SuInOut);
    ResizeParticlesSize(Np+newnp2,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    TmgStart(Timers,TMG_SuInOut);
  }

  //-Updates Velocity data of inout zones according to current timestep.
  InOut->UpdateVelData(newtimestep);
  //-Updates Zsurf data of inout zones according to current timestep.
  InOut->UpdateZsurfData(newtimestep,false);

  //-Create and remove inout particles.
  unsigned newnp=0;
  {
    byte *zsurfok=NULL;
    //-Creates list with current inout particles and normal fluid (no periodic) in inout zones.
    int* inoutpart=ArraysGpu->ReserveInt();
    const unsigned inoutcountpre=InOut->CreateListGpu(Np-Npb,Npb,Posxyg,Poszg,Codeg,GpuParticlesSize,inoutpart);

    //-Updates code of inout particles according its position and create new inlet particles when refilling=false.
    byte *newizoneg=ArraysGpu->ReserveByte();
    newnp=InOut->ComputeStepGpu(inoutcountpre,inoutpart,IdMax+1,GpuParticlesSize
      ,Np,Posxyg,Poszg,Dcellg,Codeg,Idpg,zsurfok,Velrhopg,newizoneg,this);
    ArraysGpu->Free(newizoneg);  newizoneg=NULL;

    //-Creates new inlet particles using advanced refilling mode.
    if(InOut->Use_RefillAdvanced()){
      //-Creates new inlet particles using advanced refilling mode.
      float   *prodistg =ArraysGpu->ReserveFloat();
      double2 *proposxyg=ArraysGpu->ReserveDouble2();
      double  *proposzg =ArraysGpu->ReserveDouble();
      newnp+=InOut->ComputeStepFillingGpu(Nstep,stepdt,inoutcountpre,inoutpart
        ,IdMax+1+newnp,GpuParticlesSize,Np+newnp,Posxyg,Poszg,Dcellg,Codeg,Idpg,Velrhopg
        ,zsurfok,prodistg,proposxyg,proposzg,Timers);
      ArraysGpu->Free(prodistg);
      ArraysGpu->Free(proposxyg);
      ArraysGpu->Free(proposzg);
    }
    //-Free arrays.
    ArraysGpu->Free(inoutpart);
    ArraysGpu->Free(zsurfok);
  }

  //-Updates new particle values for Laminar+SPS.
  if(SpsTaug)cudaMemset(SpsTaug+Np,0,sizeof(tsymatrix3f)*newnp);

  //-Updates number of particles.
  if(newnp){
    Np+=newnp;
    TotalNp+=newnp;
    InOut->AddNewNp(newnp);
    IdMax=unsigned(TotalNp-1);
  }

  //-Updates divide information.
  TmgStop(Timers,TMG_SuInOut);
  RunCellDivide(true);
  TmgStart(Timers,TMG_SuInOut);

  //-Updates inout particle data according inlet configuration.
  InOutUpdatePartsData(newtimestep);

  //-Saves files per PART.
  if(TimeStep+stepdt>=TimePartNext)InOut->SavePartFiles(Part);

  TmgStop(Timers,TMG_SuInOut);
}

//==============================================================================
/// Updates inout particle data according inlet configuration.
//==============================================================================
void JSphGpuSingle::InOutUpdatePartsData(double timestepnew){
  //-Create list of current inout particles (normal and periodic).
  int* inoutpart=ArraysGpu->ReserveInt();
  const unsigned inoutcount=InOut->CreateListSimpleGpu(Np-Npb,Npb,Codeg,GpuParticlesSize,inoutpart);
  InOut->SetCurrentNp(inoutcount);

  //-Updates velocity and rhop (with analytical solution).
  if(InOut->Use_AnalyticalData()){
    float *zsurfpart=NULL;
    //-Updates velocity and rhop (with analytical solution).
    InOut->SetAnalyticalDataGpu(float(timestepnew),inoutcount,inoutpart,Posxyg,Poszg,Codeg,Idpg,zsurfpart,Velrhopg);
    //-Free array.
    ArraysGpu->Free(zsurfpart); zsurfpart=NULL;
  }

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->Use_ExtrapolatedData())InOutExtrapolateData(inoutcount,inoutpart);

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->Use_InterpolatedVel())
    InOut->InterpolateVelGpu(float(timestepnew),inoutcount,inoutpart,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhopM1g)InOut->UpdateVelrhopM1Gpu(inoutcount,inoutpart,Velrhopg,VelrhopM1g);

  //-Free array for inoutpart list.
  ArraysGpu->Free(inoutpart); inoutpart=NULL;
}

//==============================================================================
/// Calculates extrapolated data for inlet/outlet particles from fluid domain.
/// Calcula datos extrapolados en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutExtrapolateData(unsigned inoutcount,const int *inoutpart){
  const float4 *planesg =InOut->GetPlanesg();
  const byte   *cfgzoneg=InOut->GetCfgZoneg();
  const float  *widthg  =InOut->GetWidthg();
  const float3 *dirdatag=InOut->GetDirDatag();
  const float determlimit=InOut->GetDetermLimit();
  const byte doublemode=InOut->GetExtrapolateMode();
  const byte extraprhopmask=InOut->GetExtrapRhopMask();
  const byte extrapvelmask =InOut->GetExtrapVelMask();
  cusphinout::Interaction_InOutExtrap(doublemode,Simulate2D,TKernel
    ,inoutcount,inoutpart,cfgzoneg,extraprhopmask,extrapvelmask
    ,planesg,widthg,dirdatag,determlimit
    ,DivData,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
}

//==============================================================================
/// Calculates extrapolated data for boundary particles from fluid domain.
/// Calcula datos extrapolados en el contorno para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::BoundCorrectionData(){
  TmgStart(Timers,TMG_SuBoundCorr);
  const unsigned n=BoundCorr->GetCount();
  const float determlimit=BoundCorr->GetDetermLimit();
  const byte doublemode=BoundCorr->GetExtrapolateMode();
  for(unsigned c=0;c<n;c++){
    const JSphBoundCorrZone* zo=BoundCorr->GetMkZone(c);
    const typecode boundcode=zo->GetBoundCode();
    const tfloat4 plane=TPlane3fToTFloat4(zo->GetPlane());
    const tfloat3 direction=ToTFloat3(zo->GetDirection());
    cusphinout::Interaction_BoundCorr(doublemode,Simulate2D,TKernel,NpbOk
      ,boundcode,plane,direction,determlimit
      ,DivData,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  }
  TmgStop(Timers,TMG_SuBoundCorr);
}


