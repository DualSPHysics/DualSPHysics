//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "JSphMk.h"
#include "JDsPartsInit.h"
#include "JSphInOut.h"
#include "JSphGpu_InOut_iker.h"
#include "JSphInOutPoints.h"
#include "JSimpleNeigs.h"
#include "FunctionsMath.h"
#include "FunctionsCuda.h"
#include "JDebugSphGpu.h"
#include "JAppInfo.h"
#include "JTimeControl.h"
#include "JDataArrays.h"
#include <climits>

using namespace std;

//==============================================================================
/// Mark special fluid particles to ignore.
/// Marca las particulas fluidas especiales para ignorar.
//==============================================================================
void JSphGpuSingle::InOutIgnoreFluidDef(const std::vector<unsigned>& mkfluidlist
 ,typecode* codeg)
{
  const unsigned nc=unsigned(mkfluidlist.size());
  for(unsigned c=0;c<nc;c++){
    const unsigned cmk=MkInfo->GetMkBlockByMkFluid(mkfluidlist[c]);
    if(cmk<MkInfo->Size()){
      const typecode rcode=MkInfo->Mkblock(cmk)->Code;
      const typecode rcode2=CODE_SetOutIgnore(rcode);
      //-Mark special fluid particles to ignore.
      cusphinout::InOutIgnoreFluidDef(Np,rcode,rcode2,codeg);
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
    unsigned nfilter=0;
    if(Np!=ParticlesDataDown(Np,0,true,false,NULL,nfilter))
      Run_Exceptioon("The number of particles is invalid.");
    //-Check proximity on CPU.
    InOut->InitCheckProximity(Np,newnp,Scell,AuxPos_c->cptr(),Idp_c->cptr()
      ,Code_c->ptr());
    //-Uploads updated code of particles to the GPU.
    Code_g->CuCopyFromHost(Code_c,Np);
  }
}

//==============================================================================
/// Initialises inlet/outlet conditions.
/// Inicia condiciones inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutInit(double timestepini){
  InOut->Nstep=Nstep; //-For debug.
  Timersg->TmStart(TMG_SuInOut,false);
  Log->Print("Initialising InOut...");
  if(PartBegin)Run_Exceptioon("Simulation restart not allowed when Inlet/Outlet is used.");

  //-Configures InOut zones and prepares new inout particles to create.
  const unsigned newnp=InOut->Config(timestepini,Stable,PeriActive
    ,MapRealPosMin,MapRealPosMax,MkInfo->GetCodeNewFluid()
    ,PartsInit,GaugeSystem);

  //-Mark special fluid particles to ignore. | Marca las particulas fluidas especiales para ignorar.
  InOutIgnoreFluidDef(InOut->MkFluidList,Code_g->ptr());

  //Log->Printf("++> newnp:%u",newnp);
  //-Resizes memory when it is necessary (always at the beginning).
  if(true || !CheckGpuParticlesSize(Np+newnp)){
    const unsigned newnp2=newnp+InOut->GetNpResizePlus0();
    Timersg->TmStop(TMG_SuInOut,true);
    const unsigned ndatacpu=0,ndatagpu=Np;
    ResizeParticlesSizeData(ndatacpu,ndatagpu,Np+newnp2,Np+newnp,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    Timersg->TmStart(TMG_SuInOut,false);
  }

  //-Creates initial inlet particles.
  unsigned idnext=IdMax+1;
  //-Compute initial data on CPU memory.
  InOut->LoadInitPartsData(idnext,newnp,Idp_c->ptr(),Code_c->ptr()
    ,AuxPos_c->ptr(),Velrho_c->ptr()); 
  Pos3ToPos21(newnp,AuxPos_c->cptr(),Posxy_c->ptr(),Posz_c->ptr());

  //-Copy data to GPU memory.
  Idp_g   ->CuCopyFromHostOffset(Np,Idp_c   ,0,newnp);
  Code_g  ->CuCopyFromHostOffset(Np,Code_c  ,0,newnp);
  Posxy_g ->CuCopyFromHostOffset(Np,Posxy_c ,0,newnp);
  Posz_g  ->CuCopyFromHostOffset(Np,Posz_c  ,0,newnp);
  Velrho_g->CuCopyFromHostOffset(Np,Velrho_c,0,newnp);

  //-Checks position of new particles and calculates cell.
  cusphinout::UpdatePosFluid(PeriActive,newnp,Np,Posxy_g->ptr(),Posz_g->ptr()
    ,Dcell_g->ptr(),Code_g->ptr());

  //-Updates new particle values for Laminar+SPS, mDBC...
  if(SpsTauRho2_g)SpsTauRho2_g->CuMemsetOffset(Np,0,newnp);
  if(BoundNor_g)BoundNor_g->CuMemsetOffset(Np,0,newnp);
  if(FSType_g)FSType_g->CuMemsetOffset(Np,3,newnp); //<vs_advshift>
  #ifdef AVAILABLE_DIVCLEAN
  if(PsiClean_g)PsiClean_g->CuMemsetOffset(Np,0,newnp);
  #endif
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",0,Np,Np+newnp
    ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->cptr());

  //-Updates number of particles.
  Np+=newnp;
  TotalNp+=newnp;
  InOut->AddNewNp(newnp);
  IdMax=unsigned(TotalNp-1);

  //-Checks proximity of inout particles to other particles and excludes fluid particles near the inout particles.
  InOutCheckProximity(newnp);

  //-Shows configuration.
  InOut->VisuConfig("\nInOut configuration:"," ");

  //-Updates divide information.
  Timersg->TmStop(TMG_SuInOut,true);
  RunCellDivide(true);
  Timersg->TmStart(TMG_SuInOut,false);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",1
    ,0,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr()
    ,Idp_g->cptr(),Velrho_g->cptr());

  //-Updates Velocity data of inout zones according to current timestep.
  InOut->UpdateVelData(timestepini);
  //-Updates Zsurf data of inout zones according to current timestep.
  InOut->UpdateZsurfData(timestepini,true);

  //-Updates inout particle data according inlet configuration.
  InOutUpdatePartsData(timestepini);

  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesGpu("CfgInOut_InletIni.vtk",2
    ,0,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr()
    ,Idp_g->cptr(),Velrho_g->cptr());
  
  //-Defines NpfMinimum according to the current fluid including inlet particles.
  const unsigned npfnormal=Np-NpbPer-NpfPer-CaseNbound;
  NpfMinimum=unsigned(MinFluidStop*npfnormal);
  Log->Printf("**MinFluidStop value was updated with inlet particles to %s (%g x %s)."
    ,KINT(NpfMinimum),MinFluidStop,KINT(npfnormal));

  Timersg->TmStop(TMG_SuInOut,true);
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
  Timersg->TmStart(TMG_SuInOut,false);
  //-Resizes memory when it is necessary. InOutCount is the maximum number of new inlet particles.
  if(!CheckGpuParticlesSize(Np+InOut->GetCurrentNp())){
    if(!InOut->GetNpResizePlus1())Run_Exceptioon("Allocated memory is not enough and resizing is not allowed by XML configuration (check the value inout.memoryresize.size).");
    const unsigned newnp2=InOut->GetCurrentNp()+InOut->GetNpResizePlus1();
    Timersg->TmStop(TMG_SuInOut,true);
    const unsigned ndatacpu=0,ndatagpu=Np;
    ResizeParticlesSizeData(ndatacpu,ndatagpu,Np+newnp2,Np+InOut->GetCurrentNp(),0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    Timersg->TmStart(TMG_SuInOut,false);
  }

  //-Updates Velocity data of inout zones according to current timestep.
  InOut->UpdateVelData(newtimestep);
  //-Updates Zsurf data of inout zones according to current timestep.
  InOut->UpdateZsurfData(newtimestep,false);

  //-Create and remove inout particles.
  unsigned newnp=0;
  {
    agbyte zsurfokg("-",Arrays_Gpu,InOut->Use_ZsurfNonUniform());
    agint inoutpartg("-",Arrays_Gpu,true);
    agbyte newizoneg("-",Arrays_Gpu,true);

    //-Creates list with current inout particles and normal fluid (no periodic) in inout zones.
    const unsigned inoutcountpre=InOut->CreateListGpu(Np-Npb,Npb
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr()
      ,GpuParticlesSize,inoutpartg.ptr());

    //-Computes Zsurf-ok of current inout particles when it is necessary. //<vs_meeshdat_ini>
    if(InOut->Use_ZsurfNonUniform()){
      InOut->ComputeZsurfokPartGpu(inoutcountpre,inoutpartg.cptr()
        ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr()
        ,Idp_g->cptr(),zsurfokg.ptr());
    } //<vs_meeshdat_end>

    //-Updates code of inout particles according its position and create new inlet particles when refilling=false.
    newnp=InOut->ComputeStepGpu(inoutcountpre,inoutpartg.ptr(),IdMax+1
      ,GpuParticlesSize,Np,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr()
      ,Code_g->ptr(),Idp_g->ptr(),zsurfokg.cptr(),Velrho_g->ptr()
      ,newizoneg.ptr(),this);

    //-Creates new inlet particles using advanced refilling mode.
    if(InOut->Use_RefillAdvanced()){
      //-Computes Zsurf-ok of inout points when it is necessary.       //<vs_meeshdat>
      if(zsurfokg.cptr())InOut->ComputeZsurfokPtosGpu(zsurfokg.ptr()); //<vs_meeshdat>
      //-Creates new inlet particles using advanced refilling mode.
      if(!InOut->RefillingRate || (Nstep%InOut->RefillingRate)==0){
        agfloat   prodistg ("-",Arrays_Gpu,true);
        agdouble2 proposxyg("-",Arrays_Gpu,true);
        agdouble  proposzg ("-",Arrays_Gpu,true);
        newnp+=InOut->ComputeStepFillingGpu(Nstep,stepdt,inoutcountpre,inoutpartg.ptr()
          ,IdMax+1+newnp,GpuParticlesSize,Np+newnp,Posxy_g->ptr()
          ,Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr(),Idp_g->ptr()
          ,Velrho_g->ptr()
          ,zsurfokg.cptr(),prodistg.ptr(),proposxyg.ptr(),proposzg.ptr(),Timersg);
      }
    }
  }

  //-Updates new particle values for Laminar+SPS, mDBC...
  if(SpsTauRho2_g)SpsTauRho2_g->CuMemsetOffset(Np,0,newnp);
  if(BoundNor_g)BoundNor_g->CuMemsetOffset(Np,0,newnp);
  if(FSType_g)FSType_g->CuMemsetOffset(Np,3,newnp); //<vs_advshift>
  #ifdef AVAILABLE_DIVCLEAN
  if(PsiClean_g)PsiClean_g->CuMemsetOffset(Np,0,newnp);
  #endif

  //-Updates number of particles.
  if(newnp){
    Np+=newnp;
    TotalNp+=newnp;
    InOut->AddNewNp(newnp);
    IdMax=unsigned(TotalNp-1);
  }

  //-Updates divide information.
  Timersg->TmStop(TMG_SuInOut,true);
  RunCellDivide(true);
  Timersg->TmStart(TMG_SuInOut,false);

  //-Updates inout particle data according inlet configuration.
  InOutUpdatePartsData(newtimestep);

  //-Saves files per PART.
  if(TimeStep+stepdt>=TimePartNext)InOut->SavePartFiles(Part);

  Timersg->TmStop(TMG_SuInOut,true);
}

//==============================================================================
/// Updates inout particle data according inlet configuration.
//==============================================================================
void JSphGpuSingle::InOutUpdatePartsData(double timestepnew){
  //-Create list of current inout particles (normal and periodic).
  agint inoutpartg("-",Arrays_Gpu,true);
  const unsigned inoutcount=InOut->CreateListSimpleGpu(Np-Npb,Npb
    ,Code_g->cptr(),GpuParticlesSize,inoutpartg.ptr());
  InOut->SetCurrentNp(inoutcount);

  //-Updates velocity and rhop (with analytical solution).
  if(InOut->Use_AnalyticalData()){
    agfloat zsurfpartg("-",Arrays_Gpu,InOut->Use_ZsurfNonUniform());
    //-Computes Zsurf of current inout particles when it is necessary. //<vs_meeshdat_ini>
    if(InOut->Use_ZsurfNonUniform()){
      InOut->ComputeZsurfPartGpu(inoutcount,inoutpartg.cptr(),Posxy_g->cptr()
        ,Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),zsurfpartg.ptr());
    } //<vs_meeshdat_end>
    //-Updates velocity and rhop (with analytical solution).
    InOut->SetAnalyticalDataGpu(float(timestepnew),inoutcount,inoutpartg.cptr()
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr()
      ,Idp_g->cptr(),zsurfpartg.cptr(),Velrho_g->ptr());
    
    //-Updates velocity of inout particles when it uses an special velocity profile.
    if(InOut->Use_SpecialProfile()){
      InOut->SetSpecialVelGpu(float(timestepnew),inoutcount,inoutpartg.cptr()
        ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->ptr());
    }
  }

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->Use_ExtrapolatedData())InOutExtrapolateData(inoutcount,inoutpartg.cptr());

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->Use_InterpolatedVel()){
    InOut->InterpolateVelGpu(float(timestepnew),inoutcount,inoutpartg.cptr()
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr()
      ,Velrho_g->ptr());
  }

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhoM1_g){
    InOut->UpdateVelrhopM1Gpu(inoutcount,inoutpartg.cptr(),Velrho_g->cptr()
      ,VelrhoM1_g->ptr());
  }
}

//==============================================================================
/// Calculates extrapolated data for inlet/outlet particles from fluid domain.
/// Calcula datos extrapolados en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphGpuSingle::InOutExtrapolateData(unsigned inoutcount,const int* inoutpart){
  const float4* planesg =InOut->GetPlanesg();
  const byte*   cfgzoneg=InOut->GetCfgZoneg();
  const float*  widthg  =InOut->GetWidthg();
  const float3* dirdatag=InOut->GetDirDatag();
  const float determlimit=InOut->GetDetermLimit();
  const byte doublemode=InOut->GetExtrapolateMode();
  const byte extraprhopmask=InOut->GetExtrapRhopMask();
  const byte extrapvelmask =InOut->GetExtrapVelMask();
  cusphinout::Interaction_InOutExtrap(doublemode,Simulate2D,TKernel
    ,inoutcount,inoutpart,cfgzoneg,extraprhopmask,extrapvelmask
    ,planesg,widthg,dirdatag,determlimit
    ,DivData,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr()
    ,Velrho_g->ptr());
}


