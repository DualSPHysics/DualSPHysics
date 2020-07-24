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

/// \file JSphCpuSingle_InOut.cpp \brief Implements InOut functions of class \ref JSphCpuSingle.

#include "JSphCpuSingle.h"
#include "JCellDivCpuSingle.h"
#include "JArraysCpu.h"
#include "JSphInOut.h"
#include "JSphMk.h"
#include "JDsPartsInit.h"
#include "JSphBoundCorr.h"
#include "JSphInOutPoints.h"
#include "JSimpleNeigs.h"
#include "FunctionsMath.h"
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
void JSphCpuSingle::InOutIgnoreFluidDef(const std::vector<unsigned> &mkfluidlist){
  const unsigned nc=unsigned(mkfluidlist.size());
  for(unsigned c=0;c<nc;c++){
    const unsigned cmk=MkInfo->GetMkBlockByMkFluid(mkfluidlist[c]);
    if(cmk<MkInfo->Size()){
      const typecode rcode=MkInfo->Mkblock(cmk)->Code;
      const typecode rcode2=CODE_SetOutIgnore(rcode);
      //-Mark special fluid particles to ignore.
      for(unsigned p=0;p<Np;p++)if(Codec[p]==rcode)Codec[p]=rcode2;
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
void JSphCpuSingle::InOutCheckProximity(unsigned newnp){
  if(Np && newnp){
    InOut->InitCheckProximity(Np,newnp,Scell,Posc,Idpc,Codec);
  }
}

//==============================================================================
/// Initialises inlet/outlet conditions.
/// Inicia condiciones inlet/outlet.
//==============================================================================
void JSphCpuSingle::InOutInit(double timestepini){
  TmcStart(Timers,TMC_SuInOut);
  Log->Print("Initialising InOut...");
  if(PartBegin)Run_Exceptioon("Simulation restart not allowed when Inlet/Outlet is used.");

  //-Configures InOut zones and prepares new inout particles to create.
  const unsigned newnp=InOut->Config(timestepini,Stable,Simulate2D,Simulate2DPosY,PeriActive
    ,RhopZero,CteB,Gamma,Gravity,Dp,MapRealPosMin,MapRealPosMax,MkInfo->GetCodeNewFluid()
    ,PartsInit,GaugeSystem,NuxLib);

  //-Mark special fluid particles to ignore. | Marca las particulas fluidas especiales para ignorar.
  InOutIgnoreFluidDef(InOut->MkFluidList);

  //Log->Printf("++> newnp:%u",newnp);
  //-Resizes memory when it is necessary (always at the beginning).
  if(true || !CheckCpuParticlesSize(Np+newnp)){
    const unsigned newnp2=newnp+InOut->GetNpResizePlus0();
    TmcStop(Timers,TMC_SuInOut);
    ResizeParticlesSize(Np+newnp2,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    TmcStart(Timers,TMC_SuInOut);
  }

  //-Creates initial inlet particles with pos, idp, code and velrhop=0.
  unsigned idnext=IdMax+1;
  InOut->LoadInitPartsData(idnext,newnp,Idpc+Np,Codec+Np,Posc+Np,Velrhopc+Np);

  //-Checks position of new particles and calculates cell.
  for(unsigned p=Np;p<Np+newnp;p++)UpdatePos(Posc[p],0,0,0,false,p,Posc,Dcellc,Codec);

  //-Updates new particle values for Laminar+SPS.
  if(SpsTauc)memset(SpsTauc+Np,0,sizeof(tsymatrix3f)*newnp);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",0,Np,Np+newnp,Posc,Codec,Idpc,Velrhopc);

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
  if(Symmetry && InOut->GetExtrapolatedData())Run_Exceptioon("Symmetry is not allowed with inlet/outlet conditions when extrapolate option is enabled."); //<vs_syymmetry>

  //-Updates divide information and creates inout particles list.
  TmcStop(Timers,TMC_SuInOut);
  RunCellDivide(true);
  TmcStart(Timers,TMC_SuInOut);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",1,0,Np,Posc,Codec,Idpc,Velrhopc);

  //-Create list of current inout particles (normal and periodic).
  int* inoutpart=ArraysCpu->ReserveInt();
  const unsigned inoutcount=InOut->CreateListSimpleCpu(Nstep,Np-Npb,Npb,Codec,inoutpart);
  InOut->SetCurrentNp(inoutcount);

  //-Updates velocity and rhop (no extrapolated).
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataCpu(float(timestepini),true,inoutcount,inoutpart,Posc,Codec,Idpc,Velrhopc);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData(inoutcount,inoutpart);

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelCpu(float(timestepini),inoutcount,inoutpart,Posc,Codec,Idpc,Velrhopc);

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhopM1c)InOut->UpdateVelrhopM1Cpu(inoutcount,inoutpart,Velrhopc,VelrhopM1c);

  //-Free array for inoutpart list.
  ArraysCpu->Free(inoutpart); inoutpart=NULL;

  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",2,0,Np,Posc,Codec,Idpc,Velrhopc);
  TmcStop(Timers,TMC_SuInOut);
}

//==============================================================================
/// ComputeStep over inlet/outlet particles:
/// - Move particles in/out according its velocity.
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new in/out particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
void JSphCpuSingle::InOutComputeStep(double stepdt){
  //Log->Printf("%u>--------> [InOutComputeStep_000]",Nstep);
  //DgSaveVtkParticlesCpu("_ComputeStep_XX.vtk",0,0,Np,Posc,Codec,Idpc,Velrhopc);
  TmcStart(Timers,TMC_SuInOut);
  //-Resizes memory when it is necessary. InOutCount is the maximum number of new inlet particles.
  if(!CheckCpuParticlesSize(Np+InOut->GetCurrentNp())){
    if(!InOut->GetNpResizePlus1())Run_Exceptioon("Allocated memory is not enough and resizing is not allowed by XML configuration (check the value inout.memoryresize.size).");
    const unsigned newnp2=InOut->GetCurrentNp()+InOut->GetNpResizePlus1();
    TmcStop(Timers,TMC_SuInOut);
    ResizeParticlesSize(Np+newnp2,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    TmcStart(Timers,TMC_SuInOut);
  }

  //-Create and remove inout particles.
  unsigned newnp=0;
  {
    //-Creates list with current inout particles and normal fluid (no periodic) in inout zones.
    int *inoutpart=ArraysCpu->ReserveInt();
    const unsigned inoutcountpre=InOut->CreateListCpu(Nstep,Np-Npb,Npb,Posc,Idpc,Codec,inoutpart);

    //-Updates code of inout particles according its position and create new inlet particles when refilling=false.
    //if(1)for(unsigned p=0;p<Np;p++)if(Idpc[p]==4382)Log->Printf("%d>=CS_005>> vel[%d].x:%f",Nstep,p,Velrhopc[p].x);
    byte *newizone=ArraysCpu->ReserveByte();
    newnp=InOut->ComputeStepCpu(Nstep,stepdt,inoutcountpre,inoutpart,this,IdMax+1,CpuParticlesSize,Np,Posc,Dcellc,Codec,Idpc,Velrhopc,newizone);
    ArraysCpu->Free(newizone);  newizone=NULL;

    //-Creates new inlet particles using advanced refilling mode.
    if(InOut->GetRefillAdvanced()){
      float    *prodist=ArraysCpu->ReserveFloat();
      tdouble3 *propos =ArraysCpu->ReserveDouble3();
      newnp+=InOut->ComputeStepFillingCpu(Nstep,stepdt,inoutcountpre,inoutpart
        ,this,IdMax+1+newnp,CpuParticlesSize,Np+newnp,Posc,Dcellc,Codec,Idpc,Velrhopc
        ,prodist,propos);
      ArraysCpu->Free(prodist);
      ArraysCpu->Free(propos);
    }
    ArraysCpu->Free(inoutpart);
  }

  //-Updates new particle values for Laminar+SPS.
  if(SpsTauc)memset(SpsTauc+Np,0,sizeof(tsymatrix3f)*newnp);

  //-Updates number of particles.
  if(newnp){
    Np+=newnp;
    TotalNp+=newnp;
    InOut->AddNewNp(newnp);
    IdMax=unsigned(TotalNp-1);
  }

  //-Updates divide information.
  TmcStop(Timers,TMC_SuInOut);
  RunCellDivide(true);
  TmcStart(Timers,TMC_SuInOut);

  //-Create list of current inout particles (normal and periodic).
  int* inoutpart=ArraysCpu->ReserveInt();
  const unsigned inoutcount=InOut->CreateListSimpleCpu(Nstep,Np-Npb,Npb,Codec,inoutpart);
  InOut->SetCurrentNp(inoutcount);

  //-Updates zsurf.
  if(InOut->GetCalculatedZsurf())InOutCalculeZsurf();
  if(InOut->GetCalculatedZsurf() || InOut->GetVariableZsurf())InOut->UpdateZsurf(TimeStep+stepdt);

  //-Updates velocity and rhop (no extrapolated).
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataCpu(float(TimeStep+stepdt),true
    ,inoutcount,inoutpart,Posc,Codec,Idpc,Velrhopc);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData(inoutcount,inoutpart);

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelCpu(float(TimeStep+stepdt)
    ,inoutcount,inoutpart,Posc,Codec,Idpc,Velrhopc);

  //-Removes interpolated Z velocity of inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateResetZVelCpu(inoutcount,inoutpart,Codec,Velrhopc);

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhopM1c)InOut->UpdateVelrhopM1Cpu(inoutcount,inoutpart,Velrhopc,VelrhopM1c);

  //-Free array for inoutpart list.
  ArraysCpu->Free(inoutpart); inoutpart=NULL;

  //-Saves files per PART.
  if(TimeStep+stepdt>=TimePartNext)InOut->SavePartFiles(Part);

  //if(1)for(unsigned p=0;p<Np;p++)if(Idpc[p]==4382)Log->Printf("%d>=CS_FIN>> vel[%d].x:%f",Nstep,p,Velrhopc[p].x);
  TmcStop(Timers,TMC_SuInOut);
}

//==============================================================================
/// Calculates zsurf for inlet/outlet particles from fluid domain.
/// Calcula zsurf en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphCpuSingle::InOutCalculeZsurf(){
  TmcStart(Timers,TMC_SuInOut);
  for(unsigned ci=0;ci<InOut->GetCount();ci++)if(InOut->GetCountPtzPos(ci)){
    const unsigned nptz=InOut->GetCountPtzPos(ci);
    const tfloat3 *ptz=InOut->GetPtzPos(ci);
    const float maxdist=(float)InOut->GetDistPtzPos(ci);
    const float zbottom=InOut->GetZbottom(ci);
    const float zsurf=Interaction_InOutZsurf(nptz,ptz,maxdist,zbottom,DivData,Posc,Codec);
    InOut->SetInputZsurf(ci,zsurf);
  }
  TmcStop(Timers,TMC_SuInOut);
}

//==============================================================================
/// Calculates extrapolated data for inlet/outlet particles from fluid domain.
/// Calcula datos extrapolados en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphCpuSingle::InOutExtrapolateData(unsigned inoutcount,const int *inoutpart){
  const tplane3f *planes=InOut->GetPlanes();
  const byte    *cfgzone=InOut->GetCfgZone();
  const float   *width  =InOut->GetWidth();
  const tfloat3 *dirdata=InOut->GetDirData();
  const float determlimit=InOut->GetDetermLimit();
  const byte doublemode=InOut->GetExtrapolateMode();
  Interaction_InOutExtrap(doublemode,inoutcount,inoutpart,cfgzone,planes,width,dirdata,determlimit
    ,Dcellc,Posc,Codec,Idpc,Velrhopc);
}

//==============================================================================
/// Calculates extrapolated data for boundary particles from fluid domain.
/// Calcula datos extrapolados en el contorno para las particulas inlet/outlet.
//==============================================================================
void JSphCpuSingle::BoundCorrectionData(){
  TmcStart(Timers,TMC_SuBoundCorr);
  const unsigned n=BoundCorr->GetCount();
  const float determlimit=BoundCorr->GetDetermLimit();
  const byte doublemode=BoundCorr->GetExtrapolateMode();
  for(unsigned c=0;c<n;c++){
    const JSphBoundCorrZone* zo=BoundCorr->GetMkZone(c);
    const typecode boundcode=zo->GetBoundCode();
    const tplane3f plane=zo->GetPlane();
    const tfloat3 direction=ToTFloat3(zo->GetDirection());
    Interaction_BoundCorr(doublemode,boundcode,plane,direction,determlimit
      ,Posc,Codec,Idpc,Velrhopc);
  }
  TmcStop(Timers,TMC_SuBoundCorr);
}



