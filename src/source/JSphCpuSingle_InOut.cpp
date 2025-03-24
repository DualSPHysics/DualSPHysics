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

/// \file JSphCpuSingle_InOut.cpp \brief Implements InOut functions of class \ref JSphCpuSingle.

#include "JSphCpuSingle.h"
#include "JCellDivCpuSingle.h"
#include "JSphInOut.h"
#include "JSphMk.h"
#include "JDsPartsInit.h"
#include "JSphInOutPoints.h"
#include "JSimpleNeigs.h"
#include "FunctionsMath.h"
#include "JAppInfo.h"
#include "JTimeControl.h"
#include "JDataArrays.h"
#include <climits>

using namespace std;

//==============================================================================
/// Mark special fluid particles to ignore.
/// Marca las particulas fluidas especiales para ignorar.
//==============================================================================
void JSphCpuSingle::InOutIgnoreFluidDef(const std::vector<unsigned>& mkfluidlist
  ,typecode* code)
{
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
    InOut->InitCheckProximity(Np,newnp,Scell,Pos_c->cptr(),Idp_c->cptr()
      ,Code_c->ptr());
  }
}

//==============================================================================
/// Initialises inlet/outlet conditions.
/// Inicia condiciones inlet/outlet.
//==============================================================================
void JSphCpuSingle::InOutInit(double timestepini){
  InOut->Nstep=Nstep; //-For debug.
  Timersc->TmStart(TMC_SuInOut);
  Log->Print("Initialising InOut...");
  if(PartBegin)Run_Exceptioon("Simulation restart not allowed when Inlet/Outlet is used.");
  
  //-Configures InOut zones and prepares new inout particles to create.
  const unsigned newnp=InOut->Config(timestepini,Stable,PeriActive
    ,MapRealPosMin,MapRealPosMax,MkInfo->GetCodeNewFluid()
    ,PartsInit,GaugeSystem);

  //-Mark special fluid particles to ignore. | Marca las particulas fluidas especiales para ignorar.
  InOutIgnoreFluidDef(InOut->MkFluidList,Code_c->ptr());

  //Log->Printf("++> newnp:%u",newnp);
  //-Resizes memory when it is necessary (always at the beginning).
  if(true || !CheckCpuParticlesSize(Np+newnp)){
    const unsigned newnp2=newnp+InOut->GetNpResizePlus0();
    Timersc->TmStop(TMC_SuInOut);
    const unsigned ndatacpu=Np;
    ResizeParticlesSizeData(ndatacpu,Np+newnp2,Np+newnp,0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    Timersc->TmStart(TMC_SuInOut);
  }

  //-Creates initial inlet particles with pos, idp, code and velrhop=0.
  unsigned idnext=IdMax+1;
  InOut->LoadInitPartsData(idnext,newnp,Idp_c->ptr()+Np,Code_c->ptr()+Np
    ,Pos_c->ptr()+Np,Velrho_c->ptr()+Np);

  //-Checks position of new particles and calculates cell.
  for(unsigned p=Np;p<Np+newnp;p++){
    const tdouble3 ps=Pos_c->cptr()[p];
    UpdatePos(ps,0,0,0,false,p,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr());
  }

  //-Updates new particle values for Laminar+SPS, mDBC...
  if(SpsTauRho2_c)SpsTauRho2_c->MemsetOffset(Np,0,newnp);
  if(BoundNor_c)BoundNor_c->MemsetOffset(Np,0,newnp);
  if(FSType_c)FSType_c->MemsetOffset(Np,3,newnp); //<vs_advshift>
  #ifdef AVAILABLE_DIVCLEAN
  if(PsiClean_c)PsiClean_c->MemsetOffset(Np,0,newnp);
  #endif
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",0
    ,Np,Np+newnp,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr());

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
  Timersc->TmStop(TMC_SuInOut);
  RunCellDivide(true);
  Timersc->TmStart(TMC_SuInOut);
  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",1
    ,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr());

  //-Updates Velocity data of inout zones according to current timestep.
  InOut->UpdateVelData(timestepini);
  //-Updates Zsurf data of inout zones according to current timestep.
  InOut->UpdateZsurfData(timestepini,true);

  //-Updates inout particle data according inlet configuration.
  InOutUpdatePartsData(timestepini);

  if(DBG_INOUT_PARTINIT)DgSaveVtkParticlesCpu("CfgInOut_InletIni.vtk",2,0,Np
    ,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr());
  
  //-Defines NpfMinimum according to the current fluid including inlet particles.
  const unsigned npfnormal=Np-NpbPer-NpfPer-CaseNbound;
  NpfMinimum=unsigned(MinFluidStop*npfnormal);
  Log->Printf("**MinFluidStop value was updated with inlet particles to %s (%g x %s)."
    ,KINT(NpfMinimum),MinFluidStop,KINT(npfnormal));

  Timersc->TmStop(TMC_SuInOut);
}

//==============================================================================
/// ComputeStep over inlet/outlet particles:
/// - Move particles in/out according its velocity.
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new in/out particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
void JSphCpuSingle::InOutComputeStep(double stepdt){
  const double newtimestep=TimeStep+stepdt;
  InOut->Nstep=Nstep; //-For debug.
  //Log->Printf("%u>--------> [InOutComputeStep_000]",Nstep);
  //DgSaveVtkParticlesCpu("_ComputeStep_XX.vtk",0,0,Np,Posc,Codec,Idpc,Velrhopc);
  Timersc->TmStart(TMC_SuInOut);
  //-Resizes memory when it is necessary. InOutCount is the maximum number of new inlet particles.
  if(!CheckCpuParticlesSize(Np+InOut->GetCurrentNp())){
    if(!InOut->GetNpResizePlus1())Run_Exceptioon("Allocated memory is not enough and resizing is not allowed by XML configuration (check the value inout.memoryresize.size).");
    const unsigned newnp2=InOut->GetCurrentNp()+InOut->GetNpResizePlus1();
    Timersc->TmStop(TMC_SuInOut);
    const unsigned ndatacpu=Np;
    ResizeParticlesSizeData(ndatacpu,Np+newnp2,Np+InOut->GetCurrentNp(),0,false);
    CellDivSingle->SetIncreaseNp(newnp2);
    Timersc->TmStart(TMC_SuInOut);
  }

  //-Updates Velocity data of inout zones according to current timestep.
  InOut->UpdateVelData(newtimestep);
  //-Updates Zsurf data of inout zones according to current timestep.
  InOut->UpdateZsurfData(newtimestep,false);

  //-Create and remove inout particles.
  unsigned newnp=0;
  {
    acbyte zsurfok("-",Arrays_Cpu,InOut->Use_ZsurfNonUniform());
    acint inoutpart("-",Arrays_Cpu,true);
    acbyte newizone("-",Arrays_Cpu,true);

    //-Creates list with current inout particles and normal fluid (no periodic) in inout zones.
    const unsigned inoutcountpre=InOut->CreateListCpu(Np-Npb,Npb
      ,Pos_c->cptr(),Idp_c->cptr(),Code_c->ptr(),inoutpart.ptr());

    //-Computes Zsurf-ok of current inout particles when it is necessary. //<vs_meeshdat_ini>
    if(InOut->Use_ZsurfNonUniform()){
      InOut->ComputeZsurfokPartCpu(inoutcountpre,inoutpart.cptr()
        ,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),zsurfok.ptr());
    } //<vs_meeshdat_end>

    //-Updates code of inout particles according its position and create new inlet particles when refilling=false.
    newnp=InOut->ComputeStepCpu(inoutcountpre,inoutpart.ptr(),this,IdMax+1
      ,CpuParticlesSize,Np,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr(),Idp_c->ptr()
      ,zsurfok.cptr(),Velrho_c->ptr(),newizone.ptr());

    //-Creates new inlet particles using advanced refilling mode.
    if(InOut->Use_RefillAdvanced()){
      //-Computes Zsurf-ok of inout points when it is necessary.     //<vs_meeshdat>
      if(zsurfok.cptr())InOut->ComputeZsurfokPtosCpu(zsurfok.ptr()); //<vs_meeshdat>
      //-Creates new inlet particles using advanced refilling mode.
      if(!InOut->RefillingRate || (Nstep%InOut->RefillingRate)==0){
        acfloat   prodist("-",Arrays_Cpu,true);
        acdouble3 propos ("-",Arrays_Cpu,true);
        newnp+=InOut->ComputeStepFillingCpu(inoutcountpre,inoutpart.ptr()
          ,this,IdMax+1+newnp,CpuParticlesSize,Np+newnp
          ,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr(),Idp_c->ptr(),Velrho_c->ptr()
          ,zsurfok.cptr(),prodist.ptr(),propos.ptr());
      }
    }
  }

  //-Updates new particle values for Laminar+SPS, mDBC...
  if(SpsTauRho2_c)SpsTauRho2_c->MemsetOffset(Np,0,newnp);
  if(BoundNor_c)BoundNor_c->MemsetOffset(Np,0,newnp);
  if(FSType_c)FSType_c->MemsetOffset(Np,3,newnp); //<vs_advshift>
  #ifdef AVAILABLE_DIVCLEAN
  if(PsiClean_c)PsiClean_c->MemsetOffset(Np,0,newnp);
  #endif

  //-Updates number of particles.
  if(newnp){
    Np+=newnp;
    TotalNp+=newnp;
    InOut->AddNewNp(newnp);
    IdMax=unsigned(TotalNp-1);
  }

  //-Updates divide information.
  Timersc->TmStop(TMC_SuInOut);
  RunCellDivide(true);
  Timersc->TmStart(TMC_SuInOut);

  //-Updates inout particle data according inlet configuration.
  InOutUpdatePartsData(newtimestep);

  //-Saves files per PART.
  if(TimeStep+stepdt>=TimePartNext)InOut->SavePartFiles(Part);

  Timersc->TmStop(TMC_SuInOut);
}

//==============================================================================
/// Updates inout particle data according inlet configuration.
//==============================================================================
void JSphCpuSingle::InOutUpdatePartsData(double timestepnew){
  //-Create list of current inout particles (normal and periodic).
  acint inoutpart("-",Arrays_Cpu,true);
  const unsigned inoutcount=InOut->CreateListSimpleCpu(Np-Npb,Npb
    ,Code_c->cptr(),inoutpart.ptr());
  InOut->SetCurrentNp(inoutcount);

  //-Updates velocity and rhop (with analytical solution).
  if(InOut->Use_AnalyticalData()){
    acfloat zsurfpart("-",Arrays_Cpu,InOut->Use_ZsurfNonUniform());
    //-Computes Zsurf of current inout particles when it is necessary. //<vs_meeshdat_ini>
    if(InOut->Use_ZsurfNonUniform()){
      InOut->ComputeZsurfPartCpu(inoutcount,inoutpart.cptr(),Pos_c->cptr()
        ,Code_c->cptr(),Idp_c->cptr(),zsurfpart.ptr());
    } //<vs_meeshdat_end>
    //-Updates velocity and rhop (with analytical solution).
    InOut->SetAnalyticalDataCpu(float(timestepnew),inoutcount,inoutpart.cptr()
      ,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),zsurfpart.cptr()
      ,Velrho_c->ptr());

    //-Updates velocity of inout particles when it uses an special velocity profile.
    if(InOut->Use_SpecialProfile()){
      InOut->SetSpecialVelCpu(float(timestepnew),inoutcount,inoutpart.cptr()
        ,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->ptr());
    }
  }

  //-Calculates extrapolated velocity and/or rhop for inlet/outlet particles from fluid domain.
  if(InOut->Use_ExtrapolatedData())InOutExtrapolateData(inoutcount,inoutpart.cptr());

  //-Calculates interpolated velocity for inlet/outlet particles.
  if(InOut->Use_InterpolatedVel()){
    InOut->InterpolateVelCpu(float(timestepnew),inoutcount,inoutpart.cptr()
      ,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->ptr());
  }

  //-Updates velocity and rhop of M1 variables starting from current velocity and rhop when Verlet is used. 
  if(VelrhoM1_c){
    InOut->UpdateVelrhopM1Cpu(inoutcount,inoutpart.cptr(),Velrho_c->cptr()
      ,VelrhoM1_c->ptr());
  }
}

//==============================================================================
/// Calculates extrapolated data for inlet/outlet particles from fluid domain.
/// Calcula datos extrapolados en el fluido para las particulas inlet/outlet.
//==============================================================================
void JSphCpuSingle::InOutExtrapolateData(unsigned inoutcount,const int* inoutpart){
  const tplane3f* planes=InOut->GetPlanes();
  const byte*     cfgzone=InOut->GetCfgZone();
  const float*    width  =InOut->GetWidth();
  const tfloat3*  dirdata=InOut->GetDirData();
  const float determlimit=InOut->GetDetermLimit();
  const byte doublemode=InOut->GetExtrapolateMode();
  Interaction_InOutExtrap(doublemode,inoutcount,inoutpart,cfgzone,planes,width
    ,dirdata,determlimit,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,Idp_c->cptr(),Velrho_c->ptr());
}



