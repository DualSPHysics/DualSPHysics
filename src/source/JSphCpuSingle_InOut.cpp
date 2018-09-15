//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "JSphBoundCorr.h"
#include "JSphInOutPoints.h"
#include "JSimpleNeigs.h"
#include "FunctionsMath.h"
#include "JFormatFiles2.h"
#include "JAppInfo.h"
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
  const char met[]="InOutCheckProximity";
  //-Look for nearby particles.
  const double disterror=Dp*0.8;
  JSimpleNeigs neigs(Np,Posc,Scell);
  unsigned* errpart=(unsigned*)InOutPartc; //-Use InOutPartc like auxiliary memory.
  memset(errpart,0,sizeof(int)*Np);
  const unsigned pini=Np-newnp;
  for(unsigned p=pini;p<Np;p++){//-Only inout particles.
    const unsigned n=neigs.NearbyPositions(Posc[p],p,disterror);
    const unsigned *selpos=neigs.GetSelectPos();
    for(unsigned cp=0;cp<n;cp++)errpart[selpos[cp]]=1;
  }
  //-Obtain number and type of nearby particles.
  unsigned nfluid=0,nfluidinout=0,nbound=0;
  for(unsigned p=0;p<Np;p++)if(errpart[p]){
    const typecode cod=Codec[p];
    if(CODE_IsFluid(cod)){
      if(CODE_IsFluidNotInout(cod)){ //-Normal fluid.
        errpart[p]=1;
        nfluid++;
      }
      else{ //-Inout fluid.
        errpart[p]=2;
        nfluidinout++;
      }
    }
    else{ //-Boundary.
      errpart[p]=3;
      nbound++;
    } 
  }
  //-Saves VTK file with nearby particles and check errors.
  if(nfluid+nfluidinout+nbound>0){
    const unsigned n=nfluid+nfluidinout+nbound;
    tfloat3* vpos=new tfloat3[n];
    byte* vtype=new byte[n];
    unsigned pp=0;
    for(unsigned p=0;p<Np;p++)if(errpart[p]){
      vpos[pp]=ToTFloat3(Posc[p]);
      vtype[pp]=byte(errpart[p]);
      pp++;
    }
    std::vector<JFormatFiles2::StScalarData> fields;
    fields.push_back(JFormatFiles2::DefineField("ErrorType",JFormatFiles2::UChar8,1,vtype));
    const string filevtk=AppInfo.GetDirOut()+(n>nfluid? "CfgInOut_ErrorParticles.vtk": "CfgInOut_ExcludedParticles.vtk");
    JFormatFiles2::SaveVtk(filevtk,n,vpos,fields);
    delete[] vpos;  vpos=NULL;
    delete[] vtype; vtype=NULL;
    if(n>nfluid){
      Log->AddFileInfo(filevtk,"Saves error fluid and boundary particles too close to inout particles.");
      RunException(met,"There are inout fluid or boundary particles too close to inout particles. Check VTK file CfgInOut_ErrorParticles.vtk with excluded particles.");
    }
    else{
      Log->AddFileInfo(filevtk,"Saves excluded fluid particles too close to inout particles.");
      Log->PrintfWarning("%u fluid particles were excluded since they are too close to inout particles. Check VTK file CfgInOut_ExcludedParticles.vtk",nfluid);
      //-Mark fluid particles to ignore.
      for(unsigned p=0;p<Np;p++)if(errpart[p]==1){
        Codec[p]=CODE_SetOutIgnore(Codec[p]); //-Mark fluid particles to ignore.
      }
    }
  }
}

//==============================================================================
/// Initialises inlet/outlet conditions.
/// Inicia condiciones inlet/outlet.
//==============================================================================
void JSphCpuSingle::InOutInit(double timestepini){
  const char met[]="InOutInit";
  TmcStart(Timers,TMC_SuInOut);
  Log->Print("InOut configuration:");

  //-Prepares particle data to define inout points starting from special fluid particles.
  JSphInOutPointsParticles partdata;
  if(InOut->MkFluidList.size()>0){
    partdata.Config(MkInfo,Np,Posc,Codec);
  }

  //-Configures InOut zones and prepares new inout particles to create.
  const unsigned newnp=InOut->Config(timestepini,Stable,Simulate2D,Simulate2DPosY,PeriActive,RhopZero,CteB,Gamma,Gravity,Dp,MapRealPosMin,MapRealPosMax,MkInfo->GetCodeNewFluid(),&partdata);

  //-Mark special fluid particles to ignore. | Marca las particulas fluidas especiales para ignorar.
  InOutIgnoreFluidDef(InOut->MkFluidList);
  partdata.Reset();

  //Log->Printf("++> newnp:%u",newnp);
  //-Resizes memory when it is necessary.
  if(!CheckCpuParticlesSize(Np+newnp)){
    const unsigned newnp2=newnp+InOut->CalcResizeNp(timestepini);
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
  InOut->VisuConfig(""," ");
  //-Updates divide information and creates inout particles list.
  RunCellDivide(true);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",1,0,Np,Posc,Codec,Idpc,Velrhopc);

  //-Updates velocity and rhop (no extrapolated).
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataCpu(float(timestepini),true,InOutCount,InOutPartc,Posc,Codec,Idpc,Velrhopc);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData();

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelCpu(float(timestepini),InOutCount,InOutPartc,Posc,Codec,Idpc,Velrhopc);

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhopM1c)InOut->UpdateVelrhopM1Cpu(InOutCount,InOutPartc,Velrhopc,VelrhopM1c);

  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",2,0,Np,Posc,Codec,Idpc,Velrhopc);
  TmcStop(Timers,TMC_SuInOut);
  //Log->Print("--------> [InOutInit_fin]");
}

//==============================================================================
/// ComputeStep over inlet/outlet particles:
/// - Move particles in/out according its velocity.
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new in/out particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
void JSphCpuSingle::InOutComputeStep(double stepdt){
  const char met[]="InOutComputeStep";
  //Log->Printf("%u>--------> [InOutComputeStep_000]",Nstep);
  //DgSaveVtkParticlesCpu("_ComputeStep_XX.vtk",0,0,Np,Posc,Codec,Idpc,Velrhopc);
  TmcStart(Timers,TMC_SuInOut);
  //-Resizes memory when it is necessary. InOutCount is the maximum number of new inlet particles.
  if(!CheckCpuParticlesSize(Np+InOutCount)){
    const unsigned newnp2=InOutCount+InOut->CalcResizeNp(TimeStep+stepdt);
    TmcStop(Timers,TMC_SuInOut);
    ResizeParticlesSize(Np+newnp2,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    TmcStart(Timers,TMC_SuInOut);
  }

  //-Removes interpolated Z velocity of inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateResetZVelCpu(InOutCount,InOutPartc,Codec,Velrhopc);

  //-Updates position of inout particles according its velocity and create new inlet particles.
  //DgSaveVtkParticlesCpu("_ComputeStep_XX.vtk",1,0,Np,Posc,Codec,Idpc,Velrhopc);
  unsigned newnp=0;
  if(InOut->GetUseRefilling()){
    float    *prodist=ArraysCpu->ReserveFloat();
    tdouble3 *propos =ArraysCpu->ReserveDouble3();
    newnp=InOut->ComputeStepFillingCpu(Nstep,stepdt,InOutCount,InOutPartc,this,IdMax+1,CpuParticlesSize,Np,Posc,Dcellc,Codec,Idpc,Velrhopc,prodist,propos);
    ArraysCpu->Free(prodist);
    ArraysCpu->Free(propos);
  }
  else newnp=InOut->ComputeStepCpu(Nstep,stepdt,InOutCount,InOutPartc,this,IdMax+1,CpuParticlesSize,Np,Posc,Dcellc,Codec,Idpc,Velrhopc);
  //DgSaveVtkParticlesCpu("_ComputeStep_XX.vtk",2,0,Np,Posc,Codec,Idpc,Velrhopc);

  //-Updates new particle values for Laminar+SPS.
  if(SpsTauc)memset(SpsTauc+Np,0,sizeof(tsymatrix3f)*newnp);

  //-Updates number of particles.
  if(newnp){
    Np+=newnp;
    TotalNp+=newnp;
    InOut->AddNewNp(newnp);
    IdMax=unsigned(TotalNp-1);
  }
  //DgSaveVtkParticlesCpu("_ComputeStep_XX.vtk",3,0,Np,Posc,Codec,Idpc,Velrhopc);
  //DgSaveVtkParticlesCpu("_ComputeStep_BBB.vtk",Nstep,0,Np,Posc,Codec,Idpc,Velrhopc);
  //-Updates divide information.
  TmcStop(Timers,TMC_SuInOut);
  RunCellDivide(true);
  TmcStart(Timers,TMC_SuInOut);
  //DgSaveVtkParticlesCpu("_ComputeStep_CCC.vtk",Nstep,0,Np,Posc,Codec,Idpc,Velrhopc);
  //RunException(met,"Stop");

  //-Updates zsurf.
  if(InOut->GetCalculatedZsurf())InOutCalculeZsurf();
  if(InOut->GetCalculatedZsurf() || InOut->GetVariableZsurf())InOut->UpdateZsurf(TimeStep+stepdt);
  //-Creates VTK file with Zsurf.
  if(TimeStep+stepdt>=TimePartNext)InOut->SaveVtkZsurf(Part);

  //-Updates velocity and rhop (no extrapolated).
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataCpu(float(TimeStep+stepdt),true,InOutCount,InOutPartc,Posc,Codec,Idpc,Velrhopc);
//  DgSaveVtkParticlesCpu("_ComputeStep_DDD.vtk",Nstep,0,Np,Posc,Codec,Idpc,Velrhopc);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData();

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelCpu(float(TimeStep+stepdt),InOutCount,InOutPartc,Posc,Codec,Idpc,Velrhopc);

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhopM1c)InOut->UpdateVelrhopM1Cpu(InOutCount,InOutPartc,Velrhopc,VelrhopM1c);

  TmcStop(Timers,TMC_SuInOut);
  //Log->Printf("%u>--------> [InOutComputeStep_fin]",Nstep);
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
    const float zsurf=Interaction_InOutZsurf(nptz,ptz,maxdist,zbottom
      ,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
      ,Posc,Codec);
    InOut->SetInputZsurf(ci,zsurf);
  }
  TmcStop(Timers,TMC_SuInOut);
}

//==============================================================================
/// Calculates extrapolated data for inlet/outlet particles from fluid domain.
/// Calcula datos extrapolados en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphCpuSingle::InOutExtrapolateData(){
  TmcStart(Timers,TMC_SuInOutExtrap);
  const tfloat4 *planes =InOut->GetPlanes();
  const byte    *cfgzone=InOut->GetCfgZone();
  const float   *width  =InOut->GetWidth();
  const tfloat3 *dirdata=InOut->GetDirData();
  const float determlimit=InOut->GetDetermLimit();
  Interaction_InOutExtrap(InOutCount,InOutPartc,cfgzone,planes,width,dirdata,determlimit
    ,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
    ,Dcellc,Posc,Codec,Idpc,Velrhopc);
  TmcStop(Timers,TMC_SuInOutExtrap);
}

//==============================================================================
/// Calculates extrapolated data for boundary particles from fluid domain.
/// Calcula datos extrapolados en el contorno para las particulas inlet/outlet.
//==============================================================================
void JSphCpuSingle::BoundCorrectionData(){
  TmcStart(Timers,TMC_SuInOutBExtrap);
  const unsigned n=BoundCorr->GetCount();
  const float determlimit=BoundCorr->GetDetermLimit();
  for(unsigned c=0;c<n;c++){
    const JSphBoundCorrZone* zo=BoundCorr->GetMkZone(c);
    const typecode boundcode=zo->GetBoundCode();
    const tfloat4 plane=zo->GetPlane();
    const tfloat3 direction=ToTFloat3(zo->GetDirection());
    Interaction_BoundCorr(boundcode,plane,direction,determlimit
      ,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
      ,Posc,Codec,Idpc,Velrhopc);
  }
  TmcStop(Timers,TMC_SuInOutBExtrap);
}



