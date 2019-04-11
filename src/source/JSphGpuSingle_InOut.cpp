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
#include "JSphPartsInit.h"
#include "JSphInOut.h"
#include "JSphBoundCorr.h"
#include "JSphGpu_InOut_ker.h"
#include "JSphInOutPoints.h"
#include "JSimpleNeigs.h"
#include "FunctionsMath.h"
#include "JFormatFiles2.h"
#include "JAppInfo.h"
#include "JTimeControl.h"
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
  const char met[]="InOutCheckProximity";
  //-Obtain particle data from GPU memory.
  if(Np!=ParticlesDataDown(Np,0,true,false))RunException(met,"The number of particles is invalid.");
  //-Look for nearby particles.
  const double disterror=Dp*0.8;
  JSimpleNeigs neigs(Np,AuxPos,Scell);
  unsigned* errpart=(unsigned*)AuxRhop; //-Use AuxRhop like auxiliary memory.
  memset(errpart,0,sizeof(int)*Np);
  const unsigned pini=Np-newnp;
  JTimeControl tc(5,60);
  for(unsigned p=pini;p<Np;p++){//-Only inout particles.
    //Log->Printf("==> InOutCheckProximity_000 ---------- p:%u",Np-p);
    const unsigned n=neigs.NearbyPositions(AuxPos[p],p,disterror);
    const unsigned *selpos=neigs.GetSelectPos();
    for(unsigned cp=0;cp<n;cp++)errpart[selpos[cp]]=1;
    if(tc.CheckTime())Log->Print(string("  ")+tc.GetInfoFinish(double(p-pini)/double(Np-pini)));
  }
  //-Obtain number and type of nearby particles.
  unsigned nfluid=0,nfluidinout=0,nbound=0;
  for(unsigned p=0;p<Np;p++)if(errpart[p]){
    const typecode cod=Code[p];
    if(CODE_IsNormal(cod)){
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
    else errpart[p]=0; //-Ignores non-normal particles.
  }
  //-Saves VTK file with nearby particles and check errors.
  if(nfluid+nfluidinout+nbound>0){
    const unsigned n=nfluid+nfluidinout+nbound;
    tfloat3* vpos=new tfloat3[n];
    byte* vtype=new byte[n];
    unsigned pp=0;
    for(unsigned p=0;p<Np;p++)if(errpart[p]){
      vpos[pp]=ToTFloat3(AuxPos[p]);
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
        Code[p]=CODE_SetOutIgnore(Code[p]); //-Mark fluid particles to ignore.
      }
      //-Uploads updated code of particles to the GPU.
      cudaMemcpy(Codeg,Code,sizeof(typecode)*Np,cudaMemcpyHostToDevice);
    }
  }
}

//==============================================================================
/// Creates list with particles in inlet/outlet zones.
/// Crea lista de particulas en zonas inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutCreateList(){
  TmgStart(Timers,TMG_SuInOut);
  InOutCount=InOut->CreateListGpu(Nstep,Np-Npb,Npb,Posxyg,Poszg,Codeg,GpuParticlesSize,InOutPartg);
  TmgStop(Timers,TMG_SuInOut);
}

//==============================================================================
/// Initialises inlet/outlet conditions.
/// Inicia condiciones inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutInit(double timestepini){
  const char met[]="InOutInit";
  TmgStart(Timers,TMG_SuInOut);
  Log->Print("InOut configuration:");
  if(PartBegin)RunException(met,"Simulation restart not allowed when Inlet/Outlet is used.");

  //-Configures InOut zones and prepares new inout particles to create.
  const unsigned newnp=InOut->Config(timestepini,Stable,Simulate2D,Simulate2DPosY,PeriActive,RhopZero,CteB,Gamma,Gravity,Dp,MapRealPosMin,MapRealPosMax,MkInfo->GetCodeNewFluid(),PartsInit);

  //-Mark special fluid particles to ignore. | Marca las particulas fluidas especiales para ignorar.
  InOutIgnoreFluidDef(InOut->MkFluidList);

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
  InOut->VisuConfig(""," ");
  //-Updates divide information.
  TmgStop(Timers,TMG_SuInOut);
  RunCellDivide(true);
  TmgStart(Timers,TMG_SuInOut);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",1,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates velocity and rhop (no extrapolated).
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataGpu(float(timestepini),true,InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData();

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelGpu(float(timestepini),InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhopM1g)InOut->UpdateVelrhopM1Gpu(InOutCount,InOutPartg,Velrhopg,VelrhopM1g);

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
  const char met[]="InOutComputeStep";
  //Log->Printf("%u>--------> [InOutComputeStep_000]",Nstep);
  //DgSaveVtkParticlesGpu("BB_ComputeStepA.vtk",DgNum,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
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

  //-Updates zsurf.
  if(InOut->GetCalculatedZsurf())InOutCalculeZsurf();
  if(InOut->GetCalculatedZsurf() || InOut->GetVariableZsurf())InOut->UpdateZsurf(TimeStep+stepdt);
  //-Creates VTK file with Zsurf.
  if(TimeStep+stepdt>=TimePartNext)InOut->SaveVtkZsurf(Part);

  //-Updates velocity and rhop (no extrapolated).
  if(InOut->GetNoExtrapolatedData())InOut->UpdateDataGpu(float(TimeStep+stepdt),true,InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->GetExtrapolatedData())InOutExtrapolateData();

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->GetInterpolatedVel())InOut->InterpolateVelGpu(float(TimeStep+stepdt),InOutCount,InOutPartg,Posxyg,Poszg,Codeg,Idpg,Velrhopg);

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhopM1g)InOut->UpdateVelrhopM1Gpu(InOutCount,InOutPartg,Velrhopg,VelrhopM1g);

  TmgStop(Timers,TMG_SuInOut);
}

//==============================================================================
/// Calculates zsurf for inlet/outlet particles from fluid domain.
/// Calcula zsurf en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutCalculeZsurf(){
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
}

//==============================================================================
/// Calculates extrapolated data for inlet/outlet particles from fluid domain.
/// Calcula datos extrapolados en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutExtrapolateData(){
  const float4 *planesg =InOut->GetPlanesg();
  const byte   *cfgzoneg=InOut->GetCfgZoneg();
  const float  *widthg  =InOut->GetWidthg();
  const float3 *dirdatag=InOut->GetDirDatag();
  const float determlimit=InOut->GetDetermLimit();
  const byte doublemode=InOut->GetExtrapolateMode();
  const byte extraprhopmask=InOut->GetExtrapRhopMask();
  const byte extrapvelmask =InOut->GetExtrapVelMask();
  cusphinout::Interaction_InOutExtrap(doublemode,Simulate2D,TKernel,CellMode
    ,InOutCount,InOutPartg,cfgzoneg,extraprhopmask,extrapvelmask
    ,planesg,widthg,dirdatag,determlimit
    ,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
    ,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
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
    cusphinout::Interaction_BoundCorr(doublemode,Simulate2D,TKernel,CellMode,NpbOk
      ,boundcode,plane,direction,determlimit
      ,CellDivSingle->GetNcells(),CellDivSingle->GetBeginCell(),CellDivSingle->GetCellDomainMin()
      ,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
  }
  TmgStop(Timers,TMG_SuBoundCorr);
}


