//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "JSphInOut.h"
#include "JSphBoundExtrap.h"
#include "JSphGpu_InOut_ker.h"
#include "JSphInOutPoints.h"
#include <climits>

using namespace std;

//==============================================================================
/// Mark special fluid particles to ignore.
/// Marca las particulas fluidas especiales para ignorar.
//==============================================================================
void JSphGpuSingle::InOutIgnoreFluidDef(const std::vector<unsigned> &mkfluidlist,typecode *code){
  const unsigned nc=unsigned(mkfluidlist.size());
  for(unsigned c=0;c<nc;c++){
    const unsigned cmk=MkInfo->GetMkBlockByMkFluid(mkfluidlist[c]);
    if(cmk<MkInfo->Size()){
      const typecode rcode=MkInfo->Mkblock(cmk)->Code;
      const typecode rcode2=CODE_SetOutIgnore(rcode);
      //-Mark special fluid particles to ignore.
      for(unsigned p=0;p<Np;p++)if(code[p]==rcode)code[p]=rcode2;
    }
  }
  //-Uploads updated code of particles to the GPU.
  if(nc)cudaMemcpy(Codeg,code,sizeof(typecode)*Np,cudaMemcpyHostToDevice);
}

//==============================================================================
/// Initialises inlet/outlet conditions.
/// Inicia condiciones inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutInit(double timestepini){
  const char met[]="InOutInit";
  //Log->Print("--------> [InOutInit_000]");
  TmgStart(Timers,TMG_SuInOut);
  Log->Print("InOut configuration:");

  //-Prepares particle data to define inout points starting from special fluid particles.
  JSphInOutPointsParticles partdata;
  if(InOut->MkFluidList.size()>0){
    if(CellOrder!=ORDER_XYZ)RunException(met,"Only order XYZ is valid for now...");
    unsigned np=ParticlesDataDown(Np,0,true,false,false);
    if(np!=Np)RunException(met,"The number of particles is invalid.");
    partdata.Config(MkInfo,Np,AuxPos,Code);
  }

  //-Configures InOut zones and prepares new inout particles to create.
  const unsigned newnp=InOut->Config(timestepini,Stable,Simulate2D,Simulate2DPosY,PeriActive,RhopZero,CteB,Gamma,Gravity,Dp,MapRealPosMin,MapRealPosMax,MkInfo->GetCodeNewFluid(),&partdata);

  //-Mark special fluid particles to ignore. | Marca las particulas fluidas especiales para ignorar.
  InOutIgnoreFluidDef(InOut->MkFluidList,Code);
  partdata.Reset();

  //-Excludes fluid particles near the inlet particles.
  Log->Print("**** PDTE: Excludes fluid particles near the inlet particles.");

  //Log->Printf("++> newnp:%u",newnp);
  //-Resizes memory when it is necessary.
  if(!CheckGpuParticlesSize(Np+newnp)){
    const unsigned newnp2=newnp+InOut->CalcResizeNp(timestepini);
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

  //-Updates new particle values for Verlet.
  if(VelrhopM1g)cudaMemset(VelrhopM1g+Np,0,sizeof(float4)*+newnp);//-VelrhopM1c is not used for inlet particles.
  //if(VelrhopM1c)memcpy(VelrhopM1c+Np,Velrhopc+Np,sizeof(tfloat4)*+newnp);
  //-Updates new particle values for Laminar+SPS.
  if(SpsTaug)cudaMemset(SpsTaug+Np,0,sizeof(tsymatrix3f)*newnp);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",0,Np,Np+newnp,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates number of particles.
  Np+=newnp;
  TotalNp+=newnp;
  InOut->AddNewNp(newnp);
  IdMax=unsigned(TotalNp-1);

  //-Shows configuration.
  InOut->VisuConfig(""," ");
  //-Updates divide information.
  RunCellDivide(true);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",1,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates velocity and rhop (no extrapolated).
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataGpu(float(timestepini),true,InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData();

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelGpu(float(timestepini),InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",2,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  TmgStop(Timers,TMG_SuInOut);
  //Log->Print("--------> [InOutInit_fin]");
}

//==============================================================================
/// ComputeStep over inlet/outlet particles:
/// - Move particles in/out according its velocity.
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new in/out particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
void JSphGpuSingle::InOutComputeStep(double stepdt){
  const char met[]="InOutComputeStep";
  //Log->Printf("%u>--------> [InOutComputeStep_000]",Nstep);
  //Log->Printf("%u]%u> ======>> BB_ComputeStepA.vtk (Np:%u)",DgNum,Nstep,Np);
  //DgSaveVtkParticlesGpu("BB_ComputeStepA.vtk",DgNum,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //DgSaveVtkParticlesGpu("_ComputeStep_XX.vtk",0,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  TmgStart(Timers,TMG_SuInOut);
  //-Resizes memory when it is necessary. InOutCount is the maximum number of new inlet particles.
  if(!CheckGpuParticlesSize(Np+InOutCount)){
    const unsigned newnp2=InOutCount+InOut->CalcResizeNp(TimeStep+stepdt);
    TmgStop(Timers,TMG_SuInOut);
    ResizeParticlesSize(Np+newnp2,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    TmgStart(Timers,TMG_SuInOut);
  }

  //-Removes interpolated Z velocity of inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateResetZVelGpu(InOutCount,InOutPartg,Codeg,Velrhopg);

  //-Updates position of in/out particles according its velocity and create new inlet particles.
  //DgSaveVtkParticlesGpu("_ComputeStep_XX.vtk",1,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  TmgStart(Timers,TMG_SuInOut_AA);
  unsigned newnp=0;
  if(InOut->GetUseRefilling()){
    float   *prodistg =ArraysGpu->ReserveFloat();
    double2 *proposxyg=ArraysGpu->ReserveDouble2();
    double  *proposzg =ArraysGpu->ReserveDouble();
    newnp=InOut->ComputeStepFillingGpu(Nstep,stepdt,InOutCount,InOutPartg,IdMax+1,GpuParticlesSize,Np,Posxyg,Poszg,Dcellg,Codeg,Idpg,Velrhopg,prodistg,proposxyg,proposzg);
    ArraysGpu->Free(prodistg);
    ArraysGpu->Free(proposxyg);
    ArraysGpu->Free(proposzg);
  }
  else newnp=InOut->ComputeStepGpu(Nstep,stepdt,InOutCount,InOutPartg,IdMax+1,GpuParticlesSize,Np,Posxyg,Poszg,Dcellg,Codeg,Idpg,Velrhopg);
  TmgStop(Timers,TMG_SuInOut_AA);
  //DgSaveVtkParticlesGpu("_ComputeStep_XX.vtk",2,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates new particle values for Verlet.
  TmgStart(Timers,TMG_SuInOut_BB);
  if(VelrhopM1g)cudaMemset(VelrhopM1g+Np,0,sizeof(tfloat4)*+newnp);//-VelrhopM1g is not used for inlet particles.
  //-Updates new particle values for Laminar+SPS.
  if(SpsTaug)cudaMemset(SpsTaug+Np,0,sizeof(tsymatrix3f)*newnp);
  //-Updates number of particles.
  if(newnp){
    Np+=newnp;
    TotalNp+=newnp;
    InOut->AddNewNp(newnp);
    IdMax=unsigned(TotalNp-1);
  }
  TmgStop(Timers,TMG_SuInOut_BB);
  //DgSaveVtkParticlesGpu("_ComputeStep_XX.vtk",3,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  //DgSaveVtkParticlesGpu("_ComputeStep_BBB.vtk",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates divide information.
  TmgStop(Timers,TMG_SuInOut);
  RunCellDivide(true);
  TmgStart(Timers,TMG_SuInOut);
  //DgSaveVtkParticlesGpu("_ComputeStep_CCC.vtk",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  //RunException(met,"Stop");

  //-Updates zsurf.
  if(InOut->GetCalculatedZsurf())InOutCalculeZsurf();
  if(InOut->GetCalculatedZsurf() || InOut->GetVariableZsurf())InOut->UpdateZsurf(TimeStep+stepdt);
  //-Creates VTK file with Zsurf.
  if(TimeStep+stepdt>=TimePartNext)InOut->SaveVtkZsurf(Part);

  //-Updates velocity and rhop (no extrapolated).
  TmgStart(Timers,TMG_SuInOut_DD);
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataGpu(float(TimeStep+stepdt),true,InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  TmgStop(Timers,TMG_SuInOut_DD);
//  DgSaveVtkParticlesGpu("_ComputeStep_DDD.vtk",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData();

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelGpu(float(TimeStep+stepdt),InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  TmgStop(Timers,TMG_SuInOut);
  //Log->Printf("%u>--------> [InOutComputeStep_fin]",Nstep);
}

//==============================================================================
/// Calculates zsurf for inlet/outlet particles from fluid domain.
/// Calcula zsurf en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutCalculeZsurf(){
  TmgStart(Timers,TMG_SuInOut_Zsurf);
  for(unsigned ci=0;ci<InOut->GetCount();ci++)if(InOut->GetCountPtzPos(ci)){
    const unsigned nptz=InOut->GetCountPtzPos(ci);
    const float3 *ptz=InOut->GetPtzPosg(ci);
    float *auxg=InOut->GetPtzAuxg(ci);
    float *auxh=InOut->GetPtzAux(ci);
    const float maxdist=(float)InOut->GetDistPtzPos(ci);
    const float zbottom=InOut->GetZbottom(ci);
    const float zsurf=cusphinout::InOutComputeZsurf(nptz,ptz,maxdist,zbottom
      ,CellMode,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
      ,Posxyg,Poszg,Codeg,auxg,auxh);
    InOut->SetInputZsurf(ci,zsurf);
  }
  TmgStop(Timers,TMG_SuInOut_Zsurf);
}

//==============================================================================
/// Calculates extrapolated data for inlet/outlet particles from fluid domain.
/// Calcula datos extrapolados en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutExtrapolateData(){
  TmgStart(Timers,TMG_SuInOutExtrap);
  const float4 *planesg =InOut->GetPlanesg();
  const byte   *cfgzoneg=InOut->GetCfgZoneg();
  const float  *widthg  =InOut->GetWidthg();
  const float3 *dirdatag=InOut->GetDirDatag();
  const float determlimit=InOut->GetDetermLimit();
  const byte extraprhopmask=InOut->GetExtrapRhopMask();
  const byte extrapvelmask =InOut->GetExtrapVelMask();
  cusphinout::Interaction_InOutExtrap_Double(Simulate2D,TKernel,CellMode
    ,InOutCount,InOutPartg,cfgzoneg,extraprhopmask,extrapvelmask
    ,planesg,widthg,dirdatag,determlimit
    ,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
    ,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  TmgStop(Timers,TMG_SuInOutExtrap);
}

//==============================================================================
/// Calculates extrapolated data for boundary particles from fluid domain.
/// Calcula datos extrapolados en el contorno para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::BoundExtrapolateData(){
  TmgStart(Timers,TMG_SuInOutBExtrap);
  const unsigned n=BoundExtrap->GetCount();
  const float determlimit=BoundExtrap->GetDetermLimit();
  for(unsigned c=0;c<n;c++){
    const JSphBoundExtrapZone* zo=BoundExtrap->GetMkZone(c);
    const typecode boundcode=zo->GetBoundCode();
    const tfloat4 plane=zo->GetPlane();
    const tfloat3 direction=ToTFloat3(zo->GetDirection());
    cusphinout::Interaction_BoundExtrap_Double(Simulate2D,TKernel,CellMode,NpbOk
      ,boundcode,plane,direction,determlimit
      ,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
      ,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  }
  TmgStop(Timers,TMG_SuInOutBExtrap);
}



