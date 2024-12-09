//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2023 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphCpu.cpp \brief Implements the class \ref JSphCpu.

#include "JSphCpu.h"
#include "JCellSearch_inline.h"
#include "JCellDivCpu.h"
#include "FunSphKernel.h"
#include "FunSphEos.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JDsMotion.h"
#include "JDsFixedDt.h"
#include "JWaveGen.h"
#include "JMLPistons.h"
#include "JRelaxZones.h"
#include "JChronoObjects.h"
#include "JDsFtForcePoints.h"
#include "JDsDamping.h"
#include "JXml.h"
#include "JDsSaveDt.h"
#include "JDsOutputTime.h"
#include "JDsAccInput.h"
#include "JDsGaugeSystem.h"
#include "JSphInOut.h"
#include "JSphShifting.h"
#include "JSphShiftingAdv.h" //<vs_advshift>

#include <climits>

using namespace std;

//##############################################################################
//# JSphCpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphCpu::JSphCpu(bool withmpi):JSph(0,withmpi){
  ClassName="JSphCpu";
  CellDiv=NULL;
  Arrays_Cpu=NULL;
  Timersc=new JDsTimersCpu;
  InitVars();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphCpu::~JSphCpu(){
  DestructorActive=true;
  FreeCpuMemoryParticles();
  FreeCpuMemoryFixed();
  delete Arrays_Cpu; Arrays_Cpu=NULL;
  delete Timersc;    Timersc=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphCpu::InitVars(){
  RunMode="";
  OmpThreads=1;

  DivData=DivDataCpuNull();

  Np=Npb=NpbOk=0;
  NpfPer=0;
  NpbPer=0;
  NpfPerM1=0;
  NpbPerM1=0;

  Idp_c=NULL;
  Code_c=NULL;
  Dcell_c=NULL;
  Pos_c=NULL;
  Velrho_c=NULL;

  PeriParent_c=NULL;

  BoundNor_c=NULL;    //-mDBC
  MotionVel_c=NULL;   //-mDBC2  //<vs_m2dbc>
  MotionAce_c=NULL;   //-mDBC2  //<vs_m2dbc>
  BoundMode_c=NULL;   //-mDBC2  //<vs_m2dbc>
  TangenVel_c=NULL;   //-mDBC2  //<vs_m2dbc>

  VelrhoM1_c=NULL;    //-Verlet
  PosPre_c=NULL;      //-Symplectic
  VelrhoPre_c=NULL;   //-Symplectic

  Ace_c=NULL;
  Ar_c=NULL;
  Press_c=NULL;
  Delta_c=NULL;
  NoPenShift_c=NULL;  //<vs_m2dbc>

  SpsTauRho2_c=NULL;  //-Laminar+SPS.
  Sps2Strain_c=NULL;  //-Laminar+SPS.

  ShiftPosfs_c=NULL;  //-Shifting.

  ShiftVel_c=NULL;    //-ShiftingAdvanced //<vs_advshift>
  FSType_c=NULL;      //-ShiftingAdvanced //<vs_advshift>
  FSNormal_c=NULL;    //-ShiftingAdvanced //<vs_advshift>
  FSMinDist_c=NULL;   //-ShiftingAdvanced //<vs_advshift>
  FSTresh_c=NULL;     //-ShiftingAdvanced //<vs_advshift>
  LCorr_c=NULL;       //-ShiftingAdvanced //<vs_advshift>
  PressSym_c=NULL;    //-ShiftingAdvanced //<vs_advshift>
  PressAsym_c=NULL;   //-ShiftingAdvanced //<vs_advshift>

  FreeCpuMemoryParticles();

  RidpMot=NULL;
  FreeCpuMemoryFixed();
}

//==============================================================================
/// Deallocate fixed memory on CPU for moving and floating bodies.
/// Libera memoria fija en cpu para moving y floating.
//==============================================================================
void JSphCpu::FreeCpuMemoryFixed(){
  MemCpuFixed=0;
  delete[] RidpMot;  RidpMot=NULL;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphCpu::AllocCpuMemoryFixed(){
  MemCpuFixed=0;
  try{
    //-Allocates memory for moving and floating particles.
    if(CaseNmoving || CaseNfloat){
      const unsigned sizen=CaseNmoving+CaseNfloat;
      RidpMot=new unsigned[sizen];  MemCpuFixed+=(sizeof(unsigned)*sizen);
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Frees CPU memory for the particles.
/// Libera memoria en CPU para particulas.
//==============================================================================
void JSphCpu::FreeCpuMemoryParticles(){
  //-Free array objects.
  delete Idp_c;         Idp_c=NULL;
  delete Code_c;        Code_c=NULL;
  delete Dcell_c;       Dcell_c=NULL;
  delete Pos_c;         Pos_c=NULL;
  delete Velrho_c;      Velrho_c=NULL;

  delete PeriParent_c;  PeriParent_c=NULL;

  delete BoundNor_c;    BoundNor_c=NULL;    //-mDBC
  delete MotionVel_c;   MotionVel_c=NULL;   //-mDBC2  //<vs_m2dbc>
  delete MotionAce_c;   MotionAce_c=NULL;   //-mDBC2  //<vs_m2dbc>
  delete BoundMode_c;   BoundMode_c=NULL;   //-mDBC2  //<vs_m2dbc>
  delete TangenVel_c;   TangenVel_c=NULL;   //-mDBC2  //<vs_m2dbc>

  delete VelrhoM1_c;    VelrhoM1_c=NULL;    //-Verlet
  delete PosPre_c;      PosPre_c=NULL;      //-Symplectic
  delete VelrhoPre_c;   VelrhoPre_c=NULL;   //-Symplectic

  delete Ace_c;         Ace_c=NULL;
  delete Ar_c;          Ar_c=NULL;
  delete Press_c;       Press_c=NULL;
  delete Delta_c;       Delta_c=NULL;
  delete NoPenShift_c;  NoPenShift_c=NULL;  //<vs_m2dbcNP>

  delete SpsTauRho2_c;  SpsTauRho2_c=NULL;  //-Laminar+SPS.
  delete Sps2Strain_c;  Sps2Strain_c=NULL;  //-Laminar+SPS.

  delete ShiftPosfs_c;  ShiftPosfs_c=NULL;  //-Shifting.

  delete ShiftVel_c;    ShiftVel_c=NULL;    //-ShiftingAdvanced //<vs_advshift>
  delete FSType_c;      FSType_c=NULL;      //-ShiftingAdvanced //<vs_advshift>
  delete FSNormal_c;    FSNormal_c=NULL;    //-ShiftingAdvanced //<vs_advshift>
  delete FSMinDist_c;   FSMinDist_c=NULL;   //-ShiftingAdvanced //<vs_advshift>
  delete FSTresh_c;     FSTresh_c=NULL;     //-ShiftingAdvanced //<vs_advshift>
  delete LCorr_c;       LCorr_c=NULL;       //-ShiftingAdvanced //<vs_advshift>
  delete PressSym_c;    PressSym_c=NULL;    //-ShiftingAdvanced //<vs_advshift>
  delete PressAsym_c;   PressAsym_c=NULL;   //-ShiftingAdvanced //<vs_advshift>
    
  //-Free CPU memory for array objects.
  CpuParticlesSize=0;
  if(Arrays_Cpu)Arrays_Cpu->Reset();
}

//==============================================================================
/// Allocate memory on CPU for the particles. 
/// Reserva memoria en CPU para las particulas. 
//==============================================================================
void JSphCpu::AllocCpuMemoryParticles(unsigned np){
  FreeCpuMemoryParticles();
  //-Calculate number of partices to allocate memory.
  CpuParticlesSize=np+PARTICLES_OVERMEMORY_MIN;
  //-Set size of arrays.
  Arrays_Cpu->SetArraySize(CpuParticlesSize);
  //-Create arrays for basic particle data. 
  Idp_c   =new acuint    ("Idpc"   ,Arrays_Cpu,true);
  Code_c  =new actypecode("Codec"  ,Arrays_Cpu,true);
  Dcell_c =new acuint    ("Dcellc" ,Arrays_Cpu,true);
  Pos_c   =new acdouble3 ("Posc"   ,Arrays_Cpu,true);
  Velrho_c=new acfloat4  ("Velrhoc",Arrays_Cpu,true);
  //-Arrays for mDBC.
  if(UseNormals){
    BoundNor_c=new acfloat3("BoundNorc",Arrays_Cpu,true);
    if(SlipMode>=SLIP_NoSlip){ //<vs_m2dbc_ini>
      MotionVel_c=new acfloat3("MotionVelc",Arrays_Cpu,true);
      MotionAce_c=new acfloat3("MotionAcec",Arrays_Cpu,true);
      BoundMode_c=new acbyte  ("BoundModec",Arrays_Cpu,false); //-NO INITIAL MEMORY.
      TangenVel_c=new acfloat3("TangenVelc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
    } //<vs_m2dbc_end>
  }
  //-Arrays for Verlet.
  if(TStep==STEP_Verlet){
    VelrhoM1_c=new acfloat4("VelrhoM1c",Arrays_Cpu,true);
  }
  //-Arrays for Symplectic.
  if(TStep==STEP_Symplectic){
    PosPre_c   =new acdouble3("PosPrec"   ,Arrays_Cpu,false); //-NO INITIAL MEMORY.
    VelrhoPre_c=new acfloat4 ("VelrhoPrec",Arrays_Cpu,false); //-NO INITIAL MEMORY.
  }
  //-Arrays for forces computation.
  Ace_c  =new acfloat3("Acec"  ,Arrays_Cpu,false); //-NO INITIAL MEMORY.
  Ar_c   =new acfloat ("Arc"   ,Arrays_Cpu,false); //-NO INITIAL MEMORY.
  Press_c=new acfloat ("Pressc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
  Delta_c=new acfloat ("Deltac",Arrays_Cpu,false); //-NO INITIAL MEMORY.
  //-Arrays for No Penetration.
  NoPenShift_c=new acfloat4("NoPenShiftc",Arrays_Cpu,false); //-NO INITIAL MEMORY.  //<vs_m2dbcNP>
  //-Arrays for Laminar+SPS.
  if(TVisco==VISCO_LaminarSPS){
    SpsTauRho2_c=new acsymatrix3f("SpsTauRho2c",Arrays_Cpu,true);
    Sps2Strain_c=new acsymatrix3f("Sps2Strainc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
  }
  //-Arrays for Shifting.
  ShiftPosfs_c=new acfloat4("ShiftPosfsc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
  //-Arrays for Advanced shifting. //<vs_advshift_ini>
  if(ShiftingAdv!=NULL){
    ShiftVel_c  =new acfloat4   ("ShiftVelc" ,Arrays_Cpu,true);
    FSType_c    =new acuint     ("FStypec"   ,Arrays_Cpu,true);
    FSNormal_c  =new acfloat3   ("FSNormalc" ,Arrays_Cpu,false); //-NO INITIAL MEMORY.
    FSMinDist_c =new acfloat    ("FSMinDistc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
    FSTresh_c   =new acfloat    ("FSTreshc"  ,Arrays_Cpu,false); //-NO INITIAL MEMORY.
    LCorr_c     =new acmatrix3d ("LCorrc"    ,Arrays_Cpu,false); //-NO INITIAL MEMORY.
    PressSym_c  =new acfloat3   ("PressSymc" ,Arrays_Cpu,false); //-NO INITIAL MEMORY.
    PressAsym_c =new acfloat3   ("PressAsymc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
    if(PeriActive && !PeriParent_c){
      PeriParent_c=new acuint("PeriParentc",Arrays_Cpu,true);
    }
  }//<vs_advshift_end>

  ////-Shows the allocated memory.
  //PrintSizeNp(CpuParticlesSize,MemCpuParticles,0);
}

//==============================================================================
/// Resizes CPU memory for particles saving current data (ndatacpu).
//==============================================================================
void JSphCpu::ResizeCpuMemoryParticlesData(unsigned ndatacpu
  ,unsigned npnew,unsigned npmin)
{
  npnew=npnew+PARTICLES_OVERMEMORY_MIN;
  //-Shows CPU memory allocation.
  const llong memcpuparticles=Arrays_Cpu->GetAllocMemoryCpu();
  const double mbparticle=(double(memcpuparticles)/MEBIBYTE)/CpuParticlesSize; //-MB por particula.
  const string txover=(npmin>1? fun::PrintStr(" (over-allocation: %.2fX)",double(npnew)/npmin): "");
  Log->Printf("**JSphCpu: Requesting CPU memory for %s particles%s: %.1f MiB (%u times)."
    ,KINT(npnew),txover.c_str(),mbparticle*npnew,Arrays_Cpu->GetCountResizeDataCalls()+1);
  //-Resizes CPU memory allocation.
  Arrays_Cpu->SetDataArraySize(npnew,ndatacpu);
  //-Updates values.
  CpuParticlesSize=npnew;
}

//==============================================================================
/// Returns the allocated memory on the CPU.
/// Devuelve la memoria reservada en CPU.
//==============================================================================
llong JSphCpu::GetAllocMemoryCpu()const{  
  llong s=JSph::GetAllocMemoryCpu();
  //-Allocated for particle arrays on CPU.
  if(Arrays_Cpu)s+=Arrays_Cpu->GetAllocMemoryCpu();
  //-Allocated in AllocCpuMemoryFixed().
  s+=MemCpuFixed;
  //-Allocated in other objects.
  if(MLPistons)s+=MLPistons->GetAllocMemoryCpu();
  return(s);
}

//==============================================================================
/// Visualize the reserved memory.
/// Visualiza la memoria reservada.
//==============================================================================
void JSphCpu::PrintAllocMemory(llong mcpu)const{
  if(!Nstep)Log->Print("");
  Log->Printf("%s allocated memory in CPU: %s (%.2f MiB)  (%u particle arrays)"
    ,(!Nstep? "Initial": "**Updated"),KINT(mcpu),double(mcpu)/MEBIBYTE,Arrays_Cpu->GetArrayCount());
  Arrays_Cpu->GetArrayCountUpdated();
}

//==============================================================================
/// Collect data from a range of particles and return the number of particles that 
/// will be less than n and eliminate the periodic ones
/// - onlynormal: Only keep the normal ones and eliminate the periodic particles.
///
/// Recupera datos de un rango de particulas y devuelve el numero de particulas que
/// sera menor que n si se eliminaron las periodicas.
/// - onlynormal: Solo se queda con las normales, elimina las particulas periodicas.
//==============================================================================
unsigned JSphCpu::GetParticlesData(unsigned n,unsigned pini,bool onlynormal
  ,unsigned* idp,tdouble3* pos,tfloat3* vel,float* rho,typecode* code
  ,const byte* filter,unsigned& npfilterdel)
{
  unsigned num=n;
  //-Copy selected values.
  if(code)Code_c->CopyToOffset(pini,code,0,n);
  if(idp) Idp_c ->CopyToOffset(pini,idp ,0,n);
  if(pos) Pos_c ->CopyToOffset(pini,pos ,0,n);
  const tfloat4* velrhoc=Velrho_c->cptr();
  if(vel && rho){
    for(unsigned p=0;p<n;p++){
      const tfloat4 vr=velrhoc[p+pini];
      vel[p]=TFloat3(vr.x,vr.y,vr.z);
      rho[p]=vr.w;
    }
  }
  else{
    if(vel)for(unsigned p=0;p<n;p++){ 
      const tfloat4 vr=velrhoc[p+pini];
      vel[p]=TFloat3(vr.x,vr.y,vr.z); 
    }
    if(rho)for(unsigned p=0;p<n;p++)rho[p]=velrhoc[p+pini].w;
  }

  //-Eliminate non-normal particles (periodic&  others).
  const bool usefilter=(filter!=NULL);
  unsigned nfilter=0;
  if(onlynormal || usefilter){
    if(!idp || !pos || !vel || !rho)Run_Exceptioon("Pointers without data.");
    actypecode* code2_c=NULL;
    typecode* code2=code;
    if(!code2){
      code2_c=new actypecode("-",Arrays_Cpu,true);
      code2=code2_c->ptr();
      Code_c->CopyToOffset(pini,code2,0,n);
    }
    unsigned ndel=0;
    for(unsigned p=0;p<n;p++){
      const bool isnormal=(!onlynormal || CODE_IsNormal(code2[p]));
      const bool selected=(isnormal && (!usefilter || (filter[p]&1)));
      if(ndel && selected){
        const unsigned p2=p-ndel;
        idp  [p2]=idp[p];
        pos  [p2]=pos[p];
        vel  [p2]=vel[p];
        rho  [p2]=rho[p];
        code2[p2]=code2[p];
      }
      if(!selected){
        ndel++;
        if(isnormal)nfilter++;//-Normal particles removed by filters.
      }
    }
    num-=ndel;
    delete code2_c; code2_c=NULL;
  }
  npfilterdel=nfilter;
  return(num);
}

//==============================================================================
/// Load the execution configuration with OpenMP.
/// Carga la configuracion de ejecucion con OpenMP.
//==============================================================================
void JSphCpu::ConfigOmp(const JSphCfgRun* cfg){
#ifdef OMP_USE
  //-Determine number of threads for host with OpenMP. | Determina numero de threads por host con OpenMP.
  if(Cpu && cfg->OmpThreads!=1){
    OmpThreads=cfg->OmpThreads;
    if(OmpThreads<=0)OmpThreads=max(omp_get_num_procs(),1);
    if(OmpThreads>OMP_MAXTHREADS)OmpThreads=OMP_MAXTHREADS;
    omp_set_num_threads(OmpThreads);
    Log->Printf("Threads by host for parallel execution: %d",omp_get_max_threads());
  }
  else{
    OmpThreads=1;
    omp_set_num_threads(OmpThreads);
  }
#else
  OmpThreads=1;
#endif
}

//==============================================================================
/// Configures execution mode in CPU.
/// Configura modo de ejecucion en CPU.
//==============================================================================
void JSphCpu::ConfigRunMode(){
  Hardware="CPU";
  //-Defines RunMode.
  RunMode="";
  if(Stable)RunMode=RunMode+(!RunMode.empty()? " - ": "") + "Stable";
  RunMode=RunMode+(!RunMode.empty()? " - ": "") + "Pos-Double";
  if(OmpThreads==1)RunMode=RunMode+(!RunMode.empty()? " - ": "") + "Single core";
  else             RunMode=RunMode+(!RunMode.empty()? " - ": "") + fun::PrintStr("OpenMP(Threads:%d)",OmpThreads); 
  //-Shows RunMode.
  Log->Print(" ");
  Log->Print(fun::VarStr("RunMode",RunMode));
  Log->Print(" ");
}

//==============================================================================
/// Initialisation of arrays and variables for execution.
/// Inicializa vectores y variables para la ejecucion.
//==============================================================================
void JSphCpu::InitRunCpu(){
  InitRun(Np,Idp_c->cptr(),Pos_c->cptr());
  if(TStep==STEP_Verlet)VelrhoM1_c->CopyFrom(Velrho_c,Np);
  if(TVisco==VISCO_LaminarSPS)SpsTauRho2_c->Memset(0,Np);
  if(MotionVel_c)MotionVel_c->Memset(0,Np); //<vs_m2dbc>
  if(MotionAce_c)MotionAce_c->Memset(0,Np); //<vs_m2dbc>
  if(ShiftVel_c)ShiftVel_c->Memset(0,Np);   //<vs_advshift>
  if(FSType_c)FSType_c->Memset(3,Np);       //<vs_advshift>
}

//==============================================================================
/// Prepares variables for interaction (pre-loop and forces).
/// Prepara variables para interaccion (pre-loop y fuerzas).
//==============================================================================
void JSphCpu::PreInteraction_Forces(TpInterStep interstep){
  Timersc->TmStart(TMC_CfPreForces);
  //-Assign memory.
  Ar_c->Reserve();
  Ace_c->Reserve();
  Press_c->Reserve();
  if(DDTArray)Delta_c->Reserve();
  if(Shifting)ShiftPosfs_c->Reserve();
  if(TVisco==VISCO_LaminarSPS)Sps2Strain_c->Reserve();
  if(SlipMode>=SLIP_NoSlip)NoPenShift_c->Reserve(); //<vs_m2dbcNP>
    if(ShiftingAdv){ //<vs_advshift_ini>
    FSMinDist_c->Reserve();
    FSNormal_c->Reserve();
    FSTresh_c->Reserve();
    LCorr_c->Reserve();
    PressSym_c->Reserve();
    PressAsym_c->Reserve();
  } //<vs_advshift_end>

  //-Initialise arrays.
  const unsigned npf=Np-Npb;
  Ar_c->Memset(0,Np);                                             //Arc[]=0
  Ace_c->Memset(0,Np);                                            //Acec[]=(0)
  if(AC_CPTR(Delta_c))Delta_c->Memset(0,Np);                      //Deltac[]=0
  if(AC_CPTR(Sps2Strain_c))Sps2Strain_c->MemsetOffset(Npb,0,npf); //Sps2Strainc[]=(0).
  if(AC_CPTR(NoPenShift_c))NoPenShift_c->Memset(0,Np);            //NoPenShiftc[]=(0) //<vs_m2dbcNP>
  
  //-Select particles for shifting.
  if(AC_CPTR(ShiftPosfs_c))Shifting->InitCpu(npf,Npb,Pos_c->cptr()
    ,ShiftPosfs_c->ptr());

  //<vs_advshift_ini>
  if((AC_CPTR(ShiftVel_c)) && interstep==INTERSTEP_SymPredictor)
    ShiftVel_c->Memset(0,Np);
  if(AC_CPTR(FSMinDist_c))FSMinDist_c->Memset(0,Np);
  if(AC_CPTR(FSNormal_c))FSNormal_c->Memset(0,Np);
  if(AC_CPTR(FSTresh_c))FSTresh_c->Memset(0,Np);
  if(AC_CPTR(LCorr_c))LCorr_c->Memset(0,Np);
  if(AC_CPTR(PressSym_c))PressSym_c->Memset(0,Np);
  if(AC_CPTR(PressAsym_c))PressAsym_c->Memset(0,Np);
  //<vs_advshift_end>

  //-Adds variable acceleration from input configuration.
  if(AccInput)AccInput->RunCpu(TimeStep,Gravity,npf,Npb,Code_c->cptr()
    ,Pos_c->cptr(),Velrho_c->cptr(),Ace_c->ptr());

  //-Prepare press values for interaction.
  {
    float* press=Press_c->ptr(); 
    const tfloat4* velrho=Velrho_c->cptr(); 
    const int n=int(Np);
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
    #endif
    for(int p=0;p<n;p++){
      press[p]=fsph::ComputePress(velrho[p].w,CSP);
    }
  }

  //-Calculate VelMax: Floating object particles are included and do not affect use of periodic condition.
  //-Calcula VelMax: Se incluyen las particulas floatings y no afecta el uso de condiciones periodicas.
  const unsigned pini=(DtAllParticles? 0: Npb);
  VelMax=CalcVelMaxOmp(Np-pini,Velrho_c->cptr()+pini);
  ViscDtMax=0;
  Timersc->TmStop(TMC_CfPreForces);
}

//==============================================================================
/// Returns maximum velocity from an array tfloat4.
/// Devuelve la velociad maxima de un array tfloat4.
//==============================================================================
float JSphCpu::CalcVelMaxSeq(unsigned np,const tfloat4* velrho)const{
  float velmax=0;
  for(unsigned p=0;p<np;p++){
    const tfloat4 v=velrho[p];
    const float v2=v.x*v.x+v.y*v.y+v.z*v.z;
    velmax=max(velmax,v2);
  }
  return(sqrt(velmax));
}

//==============================================================================
/// Returns maximum velocity from an array tfloat4 using OpenMP.
/// Devuelve la velociad maxima de un array tfloat4 usando OpenMP.
//==============================================================================
float JSphCpu::CalcVelMaxOmp(unsigned np,const tfloat4* velrho)const{
  float velmax=0;
  #ifdef OMP_USE
    if(np>OMP_LIMIT_COMPUTELIGHT){
      const int n=int(np);
      if(n<0)Run_Exceptioon("Number of values is too big.");
      float vmax=0;
      #pragma omp parallel 
      {
        float vmax2=0;
        #pragma omp for nowait
        for(int c=0;c<n;++c){
          const tfloat4 v=velrho[c];
          const float v2=v.x*v.x+v.y*v.y+v.z*v.z;
          if(vmax2<v2)vmax2=v2;
        }
        #pragma omp critical 
        {
          if(vmax<vmax2)vmax=vmax2;
        }
      }
      //-Saves result.
      velmax=sqrt(vmax);
    }
    else if(np)velmax=CalcVelMaxSeq(np,velrho);
  #else
    if(np)velmax=CalcVelMaxSeq(np,velrho);
  #endif
  return(velmax);
}

//==============================================================================
/// Frees memory allocated from Arrays_Cpu in PreInteraction_Forces().
/// Libera memoria asignada de Arrays_Cpu en PreInteraction_Forces().
//==============================================================================
void JSphCpu::PosInteraction_Forces(){
  Ar_c->Free();
  Ace_c->Free();
  Press_c->Free();
  Delta_c->Free();
  ShiftPosfs_c->Free();
  if(Sps2Strain_c)Sps2Strain_c->Free();
  if(BoundMode_c)BoundMode_c->Free(); //-Reserved in MdbcBoundCorrection(). //<vs_m2dbc>
  if(TangenVel_c)TangenVel_c->Free(); //-Reserved in MdbcBoundCorrection(). //<vs_m2dbc>
  if(NoPenShift_c)NoPenShift_c->Free(); //<vs_m2dbcNP>
  //<vs_advshift_ini>
  if(FSNormal_c)FSNormal_c->Free();
  if(FSMinDist_c)FSMinDist_c->Free();
  if(FSTresh_c)FSTresh_c->Free();
  if(LCorr_c)LCorr_c->Free();
  if(PressSym_c)PressSym_c->Free();
  if(PressAsym_c)PressAsym_c->Free();
  //<vs_advshift_end>
}

//==============================================================================
/// Perform interaction between particles. Bound-Fluid/Float
/// Realiza interaccion entre particulas. Bound-Fluid/Float
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> void JSphCpu::InteractionForcesBound
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const tfloat4* velrho,const typecode* code,const unsigned* idp
  ,float& viscdt,float* ar)const
{
  //-Initialize viscth to calculate max viscdt with OpenMP. | Inicializa viscth para calcular visdt maximo con OpenMP.
  float viscth[OMP_MAXTHREADS*OMP_STRIDE];
  for(int th=0;th<OmpThreads;th++)viscth[th*OMP_STRIDE]=0;
  //-Starts execution using OpenMP.
  const int pfin=int(pinit+n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<pfin;p1++){
    float visc=0,arp1=0;

    //-Load data of particle p1. | Carga datos de particula p1.
    const tdouble3 posp1=pos[p1];
    const tfloat4 velrhop1=velrho[p1];

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],false,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
      //---------------------------------------------------------------------------------------------
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float drx=float(posp1.x-pos[p2].x);
        const float dry=float(posp1.y-pos[p2].y);
        const float drz=float(posp1.z-pos[p2].z);
        const float rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass of particle p2 ===== 
          float massp2=MassFluid; //-Contains particle mass of incorrect fluid. | Contiene masa de particula por defecto fluid.
          bool compute=true;      //-Deactivate when using DEM and/or bound-float. | Se desactiva cuando se usa DEM y es bound-float.
          if(USE_FLOATING){
            bool ftp2=CODE_IsFloating(code[p2]);
            if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
            compute=!(USE_FTEXTERNAL && ftp2); //-Deactivate when using DEM/Chrono and/or bound-float. | Se desactiva cuando se usa DEM/Chrono y es bound-float.
          }

          if(compute){
            //-Density derivative (Continuity equation).
            const tfloat4 velrhop2=velrho[p2];
            const float dvx=velrhop1.x-velrhop2.x, dvy=velrhop1.y-velrhop2.y, dvz=velrhop1.z-velrhop2.z;
            if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

            {//-Viscosity.
              const float dot=drx*dvx + dry*dvy + drz*dvz;
              const float dot_rr2=dot/(rr2+Eta2);
              visc=max(dot_rr2,visc);
            }
          }
        }
      }
    }
    //-Sum results together. | Almacena resultados.
    if(arp1||visc){
      ar[p1]+=arp1;
      const int th=omp_get_thread_num();
      if(visc>viscth[th*OMP_STRIDE])viscth[th*OMP_STRIDE]=visc;
    }
  }
  //-Keep max value in viscdt. | Guarda en viscdt el valor maximo.
  for(int th=0;th<OmpThreads;th++)if(viscdt<viscth[th*OMP_STRIDE])viscdt=viscth[th*OMP_STRIDE];
}

//==============================================================================
/// van Albada limiter for Green DDT implementation.
//==============================================================================
float VanAlbadaLimiter(float beta){
  const float beta2=beta*beta;
  return(beta>0? (beta2+beta)/(1+beta2): 0);
}

//==============================================================================
/// Perform interaction between particles: Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Realiza interaccion entre particulas: Fluid/Float-Fluid/Float or Fluid/Float-Bound
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity
  ,bool shift,bool mdbc2
  ,bool shiftadv,bool aleform,bool ncpress> //<vs_advshift>
  void JSphCpu::InteractionForcesFluid(unsigned n,unsigned pinit,bool boundp2
  ,float visco,StDivDataCpu divdata,const unsigned* dcell
  ,const tsymatrix3f* tau,tsymatrix3f* two_strain
  ,const tdouble3* pos,const tfloat4* velrho,const typecode* code
  ,const unsigned* idp,const float* press,const tfloat3* dengradcorr
  ,const byte* boundmode,const tfloat3* tangenvel,const tfloat3* motionvel //<vs_m2dbc>
  ,const tfloat3* boundnorm //<vs_m2dbcNP> SHABA
  ,float& viscdt,float* ar,tfloat3* ace,float* delta
  ,TpShifting shiftmode,tfloat4* shiftposfs
  ,tfloat4* nopenshift
  ,unsigned* fstype,tfloat4* shiftvel,tmatrix3d* lcorr        //<vs_advshift>
  ,float* fstresh,tfloat3* presssym,tfloat3* pressasym)const //<vs_m2dbcNP> SHABA  //<vs_advshift>
{
  //-Initialize viscth to calculate viscdt maximo con OpenMP. | Inicializa viscth para calcular visdt maximo con OpenMP.
  float viscth[OMP_MAXTHREADS*OMP_STRIDE];
  for(int th=0;th<OmpThreads;th++)viscth[th*OMP_STRIDE]=0;
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<pfin;p1++){
    float visc=0,arp1=0,deltap1=0;
    tfloat3 acep1=TFloat3(0);
    tsymatrix3f two_strainp1={0,0,0,0,0,0};

    //-Variables for Shifting.
    tfloat4 shiftposfsp1;
    if(shift)shiftposfsp1=shiftposfs[p1];
    tfloat4 nopenshiftp1 = TFloat4(0); //<vs_m2dbcNP> SHABA

    //-Obtain data of particle p1 in case of floating objects. | Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1=false;     //-Indicate if it is floating. | Indica si es floating.
    if(USE_FLOATING){
      ftp1=CODE_IsFloating(code[p1]);
      if(ftp1 && tdensity!=DDT_None)deltap1=FLT_MAX; //-DDT is not applied to floating particles.
      if(ftp1 && shift)shiftposfsp1.x=FLT_MAX;  //-For floating objects do not calculate shifting. | Para floatings no se calcula shifting.
    }

    //-Variables for Advanced Shifting. //<vs_advshift_ini>
    float fstreshp1;
    tmatrix3d     lcorrp1=TMatrix3d(0);
    tfloat4       shiftvelp1=TFloat4(0);
    tfloat3       presssymp1=TFloat3(0);
    tfloat3       pressasymp1=TFloat3(0);
    if(ncpress)   presssymp1=presssym[p1];
    if(ncpress)   pressasymp1=pressasym[p1];
    if(aleform && !ftp1)shiftvelp1=shiftvel[p1];
    if(shiftadv)  fstreshp1=fstresh[p1];
    if(ncpress)   lcorrp1=lcorr[p1];
    //<vs_advshift_end>

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];
    const tfloat3 velp1=TFloat3(velrho[p1].x,velrho[p1].y,velrho[p1].z);
    const float rhop1=velrho[p1].w;
    const float pressp1=press[p1];
    
    //-Variables for Laminar+SPS.
    tsymatrix3f taup1; //-Note that taup1 is tau_a/rho_a^2.
    if(tvisco==VISCO_LaminarSPS)taup1=tau[p1];

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
      //------------------------------------------------------------------------------------------------
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float drx=float(posp1.x-pos[p2].x);
        const float dry=float(posp1.y-pos[p2].y);
        const float drz=float(posp1.z-pos[p2].z);
        const float rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass of particle p2 ===== 
          float massp2=(boundp2? MassBound: MassFluid); //-Contiene masa de particula segun sea bound o fluid.
          bool ftp2=false;    //-Indicate if it is floating | Indica si es floating.
          bool compute=true;  //-Deactivate when using DEM and if it is of type float-float or float-bound | Se desactiva cuando se usa DEM y es float-float o float-bound.
          if(USE_FLOATING){
            ftp2=CODE_IsFloating(code[p2]);
            if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
            #ifdef DELTA_HEAVYFLOATING
              if(ftp2 && tdensity==DDT_DDT && massp2<=(MassFluid*1.2f))deltap1=FLT_MAX;
            #else
              if(ftp2 && tdensity==DDT_DDT)deltap1=FLT_MAX;
            #endif
            if(ftp2 && shift && shiftmode==SHIFT_NoBound)shiftposfsp1.x=FLT_MAX; //-With floating objects do not use shifting. | Con floatings anula shifting.
            compute=!(USE_FTEXTERNAL && ftp1 && (boundp2 || ftp2)); //-Deactivate when using DEM and if it is of type float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
          }
          //-Changing the mass of boundary particle with boundmode. //<vs_m2dbc_ini>
          if(mdbc2 && boundp2 && !ftp2){
            if(boundmode[p2]==BMODE_MDBC2OFF)massp2=0;
          } //<vs_m2dbc_end>

          const tfloat4 velrhop2=velrho[p2];

          //-Velocity derivative (Momentum equation).
          if(compute && !ncpress){
            const float prs=(pressp1+press[p2])/(rhop1*velrhop2.w) 
              +(tker==KERNEL_Cubic? fsph::GetKernelCubic_Tensil(CSP,rr2,rhop1,pressp1,velrhop2.w,press[p2]): 0);
            const float p_vpm=-prs*massp2;
            acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
          }

          if(ncpress && compute){ //<vs_advshift_ini>
            const float prs=(pressp1+press[p2])/(rhop1*velrhop2.w)
              +(tker==KERNEL_Cubic? fsph::GetKernelCubic_Tensil(CSP,rr2,rhop1,pressp1,velrhop2.w,press[p2]): 0);
            const float p_vpm=-prs*massp2;
            const float ncprs=(-pressp1+press[p2])/(rhop1*velrhop2.w);
            const float ncp_vpm=-ncprs*massp2;
            presssym[p1].x+=p_vpm*frx; presssym[p1].y+=p_vpm*fry; presssym[p1].z+=p_vpm*frz;
            pressasym[p1].x+=ncp_vpm*frx; pressasym[p1].y+=ncp_vpm*fry; pressasym[p1].z+=ncp_vpm*frz;
          } //<vs_advshift_end>

          //-Density derivative (Continuity equation).
          float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
          #ifndef MDBC2_KEEPVEL
            if(mdbc2 && boundp2 && !ftp2){ //<vs_m2dbc_ini>
              const tfloat3 movvelp2=motionvel[p2];
              dvx=velp1.x-movvelp2.x; //-mDBC2 no slip.
              dvy=velp1.y-movvelp2.y;
              dvz=velp1.z-movvelp2.z;
            } //<vs_m2dbc_end>
          #endif
          if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz)*(rhop1/velrhop2.w);

          if(aleform && compute){ //<vs_advshift_ini>
            tfloat4 shiftvelp2=TFloat4(0);
            if(!boundp2 && !ftp2 && !ftp1)shiftvelp2=shiftvel[p2];      

            float massrhop=(massp2)/velrhop2.w;
            float rhozeroover1=RhopZero/rhop1;
            float divshiftp1=shiftvelp1.x*frx+shiftvelp1.y*fry+shiftvelp1.z*frz;
            float divshiftp2=shiftvelp2.x*frx+shiftvelp2.y*fry+shiftvelp2.z*frz;
            float div_pm1=divshiftp1*massrhop*rhozeroover1;
            float div_pm2=divshiftp2*massrhop*rhozeroover1;
            float dvx=shiftvelp1.x-shiftvelp2.x, dvy=shiftvelp1.y-shiftvelp2.y, dvz=shiftvelp1.z-shiftvelp2.z;

            acep1.x+=velp1.x*div_pm1;  acep1.y+=velp1.y*div_pm1;  acep1.z+=velp1.z*div_pm1;
            acep1.x+=velrhop2.x*div_pm2;  acep1.y+=velrhop2.y*div_pm2;  acep1.z+=velrhop2.z*div_pm2;

            float dotdv=massrhop*(-dvx*frx-dvy*fry-dvz*frz);
            acep1.x-= velp1.x*dotdv;   acep1.y-= velp1.y*dotdv; acep1.z-= velp1.z*dotdv;

            float dvx1=shiftvelp1.x*rhop1+shiftvelp2.x*velrhop2.w, dvy1=shiftvelp1.y*rhop1+shiftvelp2.y*velrhop2.w, dvz1=shiftvelp1.z*rhop1+shiftvelp2.z*velrhop2.w;
            arp1+=massrhop*(dvx1*frx+dvy1*fry+dvz1*frz);
            // //-Density derivative (Continuity equation).
            arp1+=massrhop*(dvx*frx+dvy*fry+dvz*frz)*rhop1;
          } //<vs_advshift_end>

          const float cbar=(float)Cs0;
          //-Density Diffusion Term (Molteni and Colagrossi 2009).
          if(tdensity==DDT_DDT && deltap1!=FLT_MAX){
            const float rhop1over2=rhop1/velrhop2.w;
            const float visc_densi=DDTkh*cbar*(rhop1over2-1.f)/(rr2+Eta2);
            const float dot3=(drx*frx+dry*fry+drz*frz);
            const float delta=visc_densi*dot3*massp2;
            //deltap1=(boundp2? FLT_MAX: deltap1+delta);
            deltap1=(boundp2 && TBoundary==BC_DBC? FLT_MAX: deltap1+delta);
          }
          //-Density Diffusion Term (Fourtakas et al 2019).
          if((tdensity==DDT_DDT2 || (tdensity==DDT_DDT2Full && !boundp2)) && deltap1!=FLT_MAX && !ftp2){
            const float rh=1.f+DDTgz*drz;
            const float drho=RhopZero*pow(rh,1.f/Gamma)-RhopZero;    
            const float visc_densi=DDTkh*cbar*((velrhop2.w-rhop1)-drho)/(rr2+Eta2);
            const float dot3=(drx*frx+dry*fry+drz*frz);
            const float delta=visc_densi*dot3*massp2/velrhop2.w;
            deltap1=(boundp2? FLT_MAX: deltap1-delta); //-blocks it makes it boil - bloody DBC
          }

          //-Shifting correction.
          if(shift && shiftposfsp1.x!=FLT_MAX){
            const float massrho=massp2/velrhop2.w;
            const bool noshift=(boundp2 && (shiftmode==SHIFT_NoBound || (shiftmode==SHIFT_NoFixed && CODE_IsFixed(code[p2]))));
            shiftposfsp1.x=(noshift? FLT_MAX: shiftposfsp1.x+massrho*frx); //-For boundary do not use shifting. | Con boundary anula shifting.
            shiftposfsp1.y+=massrho*fry;
            shiftposfsp1.z+=massrho*frz;
            shiftposfsp1.w-=massrho*(drx*frx+dry*fry+drz*frz);
          }

          //-No Penetration Correction SHABA
          if (boundp2 && mdbc2 && !ftp2) {//<vs_m2dbcNP_ini>
              float rrmag = sqrt(rr2);
              if (rrmag < 1.25f * Dp) { //-if fluid particle is less than 1.25dp from a boundary particle
                  float norm = sqrt(boundnorm[p2].x * boundnorm[p2].x + boundnorm[p2].y * boundnorm[p2].y + boundnorm[p2].z * boundnorm[p2].z);
                  float normx = boundnorm[p2].x / norm; float normy = boundnorm[p2].y / norm; float normz = boundnorm[p2].z / norm;
                  float normdist = (normx * drx + normy * dry + normz * drz);
                  if (normdist < 0.75f * norm && norm < 1.75f * Dp) {//-if normal distance is less than 0.75 boundary normal size and only first layer of bound
                      const tfloat3 movvelp2 = motionvel[p2];
                      dvx = velp1.x - movvelp2.x;
                      dvy = velp1.y - movvelp2.y;
                      dvz = velp1.z - movvelp2.z;
                      float vfc = dvx * normx + dvy * normy + dvz * normz; //-fluid velocity normal to boundary particle
                      if (vfc < 0.f) { //-if fluid particle velocity is pointing towards boundary add correction velocity
                          nopenshiftp1.w += 1.f; //-boundary particle counter for average
                          //-delta v = sum uij dot (nj cross nj)
                          nopenshiftp1.x -= (dvx * normx * normx + dvy * normx * normy + dvz * normx * normz);
                          nopenshiftp1.y -= (dvx * normx * normy + dvy * normy * normy + dvz * normy * normz);
                          nopenshiftp1.z -= (dvx * normx * normz + dvy * normy * normz + dvz * normz * normz);
                          if (normdist < 0.25f * norm) {// if normal distanne is less than 0.25 boundary normal size double correction velocity
                              nopenshiftp1.x -= (dvx * normx * normx + dvy * normx * normy + dvz * normx * normz);
                              nopenshiftp1.y -= (dvx * normx * normy + dvy * normy * normy + dvz * normy * normz);
                              nopenshiftp1.z -= (dvx * normx * normz + dvy * normy * normz + dvz * normz * normz);
                          }
                      }
                  }
              }
          }//<vs_m2dbcNP_end>

          //-Advanced shifting. //<vs_advshift_ini>
          if(shiftadv){
            const float massrho=(massp2)/velrhop2.w;        
            const float wab=fsph::GetKernel_Wab<tker>(CSP,rr2);;
            fstreshp1-=massrho*(drx*frx+dry*fry+drz*frz);
            if(ncpress && compute){
              // pou+=wab*massrho;              
              lcorr[p1].a11+=-drx*frx*massrho; lcorr[p1].a12+=-drx*fry*massrho; lcorr[p1].a13+=-drx*frz*massrho;
              lcorr[p1].a21+=-dry*frx*massrho; lcorr[p1].a22+=-dry*fry*massrho; lcorr[p1].a23+=-dry*frz*massrho;
              lcorr[p1].a31+=-drz*frx*massrho; lcorr[p1].a32+=-drz*fry*massrho; lcorr[p1].a33+=-drz*frz*massrho;
            }
          } //<vs_advshift_end>

          //===== Viscosity ===== 
          if(compute){
            if(mdbc2 && boundp2 && !ftp2){ //<vs_m2dbc_ini>
              const tfloat3 tangentvelp2=tangenvel[p2];
              dvx=velp1.x-tangentvelp2.x;
              dvy=velp1.y-tangentvelp2.y;
              dvz=velp1.z-tangentvelp2.z;
            } //<vs_m2dbc_end>
            const float dot=drx*dvx + dry*dvy + drz*dvz;
            const float dot_rr2=dot/(rr2+Eta2);
            visc=max(dot_rr2,visc);
            if(tvisco==VISCO_Artificial){//-Artificial viscosity.
              if(dot<0){
                const float amubar=KernelH*dot_rr2;  //amubar=CTE.h*dot/(rr2+CTE.eta2);
                const float robar=(rhop1+velrhop2.w)*0.5f;
                const float pi_visc=(-visco*cbar*amubar/robar)*massp2;
                acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
              }
            }
            if(tvisco==VISCO_Laminar || tvisco==VISCO_LaminarSPS){//-Laminar and Laminar+SPS viscosity.
              const float robar2=(rhop1+velrhop2.w);
              const float temp=4.f*visco/((rr2+Eta2)*robar2);  //-Simplification of: temp=2.0f*visco/((rr2+CTE.eta2)*robar); robar=(rhop1+velrhop2.w)*0.5f;
              const float vtemp=massp2*temp*(drx*frx+dry*fry+drz*frz);  
              acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
            }
            if(tvisco==VISCO_LaminarSPS){//-SPS contribution for Laminar viscosity. 
              //-SPS turbulence model.
              //-Note that taup1 is tau_a/rho_a^2 for interaction of particle a and b.
              //-And taup1 is always zero when p1 is not a fluid particle.
              float stau_xx=taup1.xx,stau_xy=taup1.xy,stau_xz=taup1.xz; 
              float stau_yy=taup1.yy,stau_yz=taup1.yz,stau_zz=taup1.zz;
              if(!boundp2 && !ftp2){//-When p2 is a fluid particle. 
                //-Note that tau[p2] is tau_b/rho_b^2 for interaction of particle a and b.
                stau_xx+=tau[p2].xx; stau_xy+=tau[p2].xy; stau_xz+=tau[p2].xz;
                stau_yy+=tau[p2].yy; stau_yz+=tau[p2].yz; stau_zz+=tau[p2].zz;
              }
              acep1.x+=massp2*(stau_xx*frx + stau_xy*fry + stau_xz*frz);
              acep1.y+=massp2*(stau_xy*frx + stau_yy*fry + stau_yz*frz);
              acep1.z+=massp2*(stau_xz*frx + stau_yz*fry + stau_zz*frz);
              //-Velocity gradients.
              if(!ftp1){//-When p1 is a fluid particle. 
                const float volp2=-massp2/velrhop2.w;
                float dv=dvx*volp2; two_strainp1.xx+=dv*frx; two_strainp1.xy+=dv*fry; two_strainp1.xz+=dv*frz;
                      dv=dvy*volp2; two_strainp1.xy+=dv*frx; two_strainp1.yy+=dv*fry; two_strainp1.yz+=dv*frz;
                      dv=dvz*volp2; two_strainp1.xz+=dv*frx; two_strainp1.yz+=dv*fry; two_strainp1.zz+=dv*frz;
                //-To compute tau terms we assume that two_strain.xy=dudy+dvdx, two_strain.xz=dudz+dwdx, two_strain.yz=dvdz+dwdy
                //-so only 6 elements are needed instead of 3x3.
              }
            }
          }
        }
      }
    }

    //<vs_advshift_ini>
    if(shiftadv)fstresh[p1]=fstreshp1;
    if(ncpress && boundp2){
      if(fstype[p1]==0 && shiftvelp1.w>0.95){
        tmatrix3d LCorr_inv=TMatrix3d(0);
        if(Simulate2D){
          tmatrix2d Lcorr2D;
          tmatrix2d Lcorr2D_inv;
          Lcorr2D.a11=lcorr[p1].a11; Lcorr2D.a12=lcorr[p1].a13;
          Lcorr2D.a21=lcorr[p1].a31; Lcorr2D.a22=lcorr[p1].a33;
          float lcorr_det=float(Lcorr2D.a11*Lcorr2D.a22-Lcorr2D.a12*Lcorr2D.a21);
          Lcorr2D_inv.a11=Lcorr2D.a22/lcorr_det; Lcorr2D_inv.a12=-Lcorr2D.a12/lcorr_det; Lcorr2D_inv.a22=Lcorr2D.a11/lcorr_det; Lcorr2D_inv.a21=-Lcorr2D.a21/lcorr_det;
          LCorr_inv.a11=Lcorr2D_inv.a11;  LCorr_inv.a13=Lcorr2D_inv.a12;
          LCorr_inv.a31=Lcorr2D_inv.a21;  LCorr_inv.a33=Lcorr2D_inv.a22;
        }
        else{
          const float determ=(float)fmath::Determinant3x3(lcorr[p1]);
          LCorr_inv = fmath::InverseMatrix3x3(lcorr[p1], determ);
        }
        acep1.x+=float(pressasym[p1].x*LCorr_inv.a11 + pressasym[p1].y*LCorr_inv.a12 + pressasym[p1].z*LCorr_inv.a13);
        acep1.y+=float(pressasym[p1].x*LCorr_inv.a21 + pressasym[p1].y*LCorr_inv.a22 + pressasym[p1].z*LCorr_inv.a23);
        acep1.z+=float(pressasym[p1].x*LCorr_inv.a31 + pressasym[p1].y*LCorr_inv.a32 + pressasym[p1].z*LCorr_inv.a33);
      }
      else{
        acep1.x+=presssym[p1].x; acep1.y+=presssym[p1].y; acep1.z+=presssym[p1].z;
      }
    }
    //<vs_advshift_end>

    //-Sum results together. | Almacena resultados.
    if(shift||arp1||acep1.x||acep1.y||acep1.z||visc){
      if(tdensity!=DDT_None){
        if(delta)delta[p1]=(delta[p1]==FLT_MAX || deltap1==FLT_MAX? FLT_MAX: delta[p1]+deltap1);
        else if(deltap1!=FLT_MAX)arp1+=deltap1;
      }
      ar[p1]+=arp1;
      ace[p1]=ace[p1]+acep1;
      const int th=omp_get_thread_num();
      if(visc>viscth[th*OMP_STRIDE])viscth[th*OMP_STRIDE]=visc;
      if(tvisco==VISCO_LaminarSPS){
        two_strain[p1].xx+=two_strainp1.xx;
        two_strain[p1].xy+=two_strainp1.xy;
        two_strain[p1].xz+=two_strainp1.xz;
        two_strain[p1].yy+=two_strainp1.yy;
        two_strain[p1].yz+=two_strainp1.yz;
        two_strain[p1].zz+=two_strainp1.zz;
      }
      if(shift)shiftposfs[p1]=shiftposfsp1;
    }
    //-No-Penetration correction SHABA
    if (mdbc2) { //<vs_m2dbcNP_end>
        if (nopenshiftp1.w > 0.f) {//-if correction required
            //-Average correction velocity over number of boundary particles
            nopenshift[p1].x = nopenshiftp1.x / nopenshiftp1.w;
            nopenshift[p1].y = nopenshiftp1.y / nopenshiftp1.w;
            nopenshift[p1].z = nopenshiftp1.z / nopenshiftp1.w;
            nopenshift[p1].w = 10; //-correction needed? yes
        }
        else {
            nopenshift[p1].w = 0;//-correction needed? no
        }
    }//<vs_m2dbcNP_end>
  }
  //-Keep max value in viscdt. | Guarda en viscdt el valor maximo.
  for(int th=0;th<OmpThreads;th++)if(viscdt<viscth[th*OMP_STRIDE])viscdt=viscth[th*OMP_STRIDE];
}

//==============================================================================
/// Perform DEM interaction between particles Floating-Bound & Floating-Floating //(DEM)
/// Realiza interaccion DEM entre particulas Floating-Bound & Floating-Floating //(DEM)
//==============================================================================
void JSphCpu::InteractionForcesDEM(unsigned nfloat,StDivDataCpu divdata
  ,const unsigned* dcell,const unsigned* ftridp,const StDemData* demdata
  ,const tdouble3* pos,const tfloat4* velrho
  ,const typecode* code,const unsigned* idp
  ,float& viscdt,tfloat3* ace)const
{
  //-Initialise demdtth to calculate max demdt with OpenMP. | Inicializa demdtth para calcular demdt maximo con OpenMP.
  float demdtth[OMP_MAXTHREADS*OMP_STRIDE];
  for(int th=0;th<OmpThreads;th++)demdtth[th*OMP_STRIDE]=-FLT_MAX;
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int nft=int(nfloat);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int cf=0;cf<nft;cf++){
    const unsigned p1=ftridp[cf];
    if(p1!=UINT_MAX){
      float demdtp1=0;
      tfloat3 acep1=TFloat3(0);

      //-Get data of particle p1.
      const tdouble3 posp1=pos[p1];
      const typecode tavp1=CODE_GetTypeAndValue(code[p1]);
      const float masstotp1=demdata[tavp1].mass;
      const float taup1=demdata[tavp1].tau;
      const float kfricp1=demdata[tavp1].kfric;
      const float restitup1=demdata[tavp1].restitu;
      const float ftmassp1=demdata[tavp1].massp;

      //-Search for neighbours in adjacent cells (first bound and then fluid+floating).
      for(byte tpfluid=0;tpfluid<=1;tpfluid++){
        const StNgSearch ngs=nsearch::Init(dcell[p1],!tpfluid,divdata);
        for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
          const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

          //-Interaction of Floating Object particles with type Fluid or Bound. | Interaccion de Floating con varias Fluid o Bound.
          //-----------------------------------------------------------------------------------------------------------------------
          for(unsigned p2=pif.x;p2<pif.y;p2++)if(CODE_IsNotFluid(code[p2]) && tavp1!=CODE_GetTypeAndValue(code[p2])){
            const float drx=float(posp1.x-pos[p2].x);
            const float dry=float(posp1.y-pos[p2].y);
            const float drz=float(posp1.z-pos[p2].z);
            const float rr2=drx*drx+dry*dry+drz*drz;
            if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
              const float rad=sqrt(rr2);

              //-Calculate max value of demdt. | Calcula valor maximo de demdt.
              const typecode tavp2=CODE_GetTypeAndValue(code[p2]);
              const float masstotp2=demdata[tavp2].mass;
              const float taup2=demdata[tavp2].tau;
              const float kfricp2=demdata[tavp2].kfric;
              const float restitup2=demdata[tavp2].restitu;
              //const StDemData* demp2=demobjs+CODE_GetTypeAndValue(code[p2]);

              const float nu_mass=(!tpfluid? masstotp1/2: masstotp1*masstotp2/(masstotp1+masstotp2)); //-Con boundary toma la propia masa del floating 1.
              const float kn=4.f/(3.f*(taup1+taup2))*sqrt(float(Dp)/4); //-Generalized rigidity - Lemieux 2008.
              const float dvx=velrho[p1].x-velrho[p2].x, dvy=velrho[p1].y-velrho[p2].y, dvz=velrho[p1].z-velrho[p2].z; //vji
              const float nx=drx/rad, ny=dry/rad, nz=drz/rad; //normal_ji               
              const float vn=dvx*nx+dvy*ny+dvz*nz; //vji.nji
              const float demvisc=0.2f/(3.21f*(pow(nu_mass/kn,0.4f)*pow(fabs(vn),-0.2f))/40.f);
              if(demdtp1<demvisc)demdtp1=demvisc;

              const float over_lap=1.0f*float(Dp)-rad; //-(ri+rj)-|dij|
              if(over_lap>0.0f){ //-Contact.
                //-Normal.
                const float eij=(restitup1+restitup2)/2;
                const float gn=-(2.f*log(eij)*sqrt(nu_mass*kn))/(sqrt(float(PI)+log(eij)*log(eij))); //-Generalized damping - Cummins 2010.
                //const float gn=0.08f*sqrt(nu_mass*sqrt(float(Dp)/2)/((taup1+taup2)/2)); //-Generalized damping - Lemieux 2008.
                const float rep=kn*pow(over_lap,1.5f);
                const float fn=rep-gn*pow(over_lap,0.25f)*vn;
                float acef=fn/ftmassp1; //-Divides by the mass of particle to obtain the acceleration.
                acep1.x+=(acef*nx); acep1.y+=(acef*ny); acep1.z+=(acef*nz); //-Force is applied in the normal between the particles.
                //-Tangential.
                const float dvxt=dvx-vn*nx, dvyt=dvy-vn*ny, dvzt=dvz-vn*nz; //Vji_t
                const float vt=sqrt(dvxt*dvxt + dvyt*dvyt + dvzt*dvzt);
                float tx=0, ty=0, tz=0; //-Tang vel unit vector.
                if(vt!=0){ tx=dvxt/vt; ty=dvyt/vt; tz=dvzt/vt; }
                const float ft_elast=(kn*float(DemDtForce)-gn)*vt/3.5f; //-Elastic frictional string -->  ft_elast=2*(kn*fdispl-gn*vt)/7; fdispl=dtforce*vt;
                const float kfric_ij=(kfricp1+kfricp2)/2;
                float ft=kfric_ij*fn*tanh(vt*8);  //-Coulomb.
                ft=(ft<ft_elast? ft: ft_elast);   //-Not above yield criteria, visco-elastic model.
                acef=ft/ftmassp1; //-Divides by the mass of particle to obtain the acceleration.
                acep1.x+=(acef*tx); acep1.y+=(acef*ty); acep1.z+=(acef*tz);
              } 
            }
          }
        }
      }
      //-Sum results together. | Almacena resultados.
      if(acep1.x||acep1.y||acep1.z){
        ace[p1]=ace[p1]+acep1;
        const int th=omp_get_thread_num();
        if(demdtth[th*OMP_STRIDE]<demdtp1)demdtth[th*OMP_STRIDE]=demdtp1;
      }
    }
  }
  //-Update viscdt with max value of viscdt or demdt* | Actualiza viscdt con el valor maximo de viscdt y demdt*.
  float demdt=demdtth[0];
  for(int th=1;th<OmpThreads;th++)if(demdt<demdtth[th*OMP_STRIDE])demdt=demdtth[th*OMP_STRIDE];
  if(viscdt<demdt)viscdt=demdt;
}

//<vs_advshift_ini>
//==============================================================================
/// Identify probably free'surface particles.   
//==============================================================================
void JSphCpu::ComputeFsType(unsigned n,unsigned pini,bool sim2d
  ,const float* fstresh,unsigned* fstype)const
{
  const int pfin=int(pini+n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=int(pini);p<pfin;p++){
    //-Find particles that are probably on the free-surface.
    unsigned fstypep1=0;
    if(sim2d){
      if(fstresh[p]<1.7) fstypep1=2;
      if(fstresh[p]<1.1) fstypep1=3;
    }
    else{
      if(fstresh[p]<2.75) fstypep1=2;
      if(fstresh[p]<1.8 ) fstypep1=3;
    }
    fstype[p]=fstypep1;
  }
}
//<vs_advshift_end>

//==============================================================================
/// Computes sub-particle stress tensor divided by rho^2 (tau/rho^2) for SPS 
/// turbulence model.   
//==============================================================================
void JSphCpu::ComputeSpsTau(unsigned n,unsigned pini,const tfloat4* velrho
  ,const tsymatrix3f* sps2strain,tsymatrix3f* tau_rho2)const
{
  const int pfin=int(pini+n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=int(pini);p<pfin;p++){
    const tsymatrix3f two_strain=sps2strain[p];
    const float pow1=two_strain.xx*two_strain.xx
                   + two_strain.yy*two_strain.yy
                   + two_strain.zz*two_strain.zz;
    const float prr= two_strain.xy*two_strain.xy
                   + two_strain.xz*two_strain.xz
                   + two_strain.yz*two_strain.yz + pow1+pow1;
    const float visc_sps=SpsSmag*sqrt(prr);
    const float div_u=two_strain.xx+two_strain.yy+two_strain.zz;
    const float sps_k=visc_sps*div_u; //-Factor 2/3 is included in SpsSmag constant.
    const float sps_blin=SpsBlin*prr;
    const float sumsps=-(sps_k+sps_blin);
    const float twovisc_sps=(visc_sps+visc_sps);
    const float one_rho=1.0f/velrho[p].w;
    //-Computes new values of tau/rho^2.
    tau_rho2[p].xx=one_rho*(twovisc_sps*two_strain.xx +sumsps);
    tau_rho2[p].xy=one_rho*(visc_sps  * two_strain.xy);
    tau_rho2[p].xz=one_rho*(visc_sps  * two_strain.xz);
    tau_rho2[p].yy=one_rho*(twovisc_sps*two_strain.yy +sumsps);
    tau_rho2[p].yz=one_rho*(visc_sps  * two_strain.yz);
    tau_rho2[p].zz=one_rho*(twovisc_sps*two_strain.zz +sumsps);
  }
}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid (forces and DEM).
/// Interaccion Fluid-Fluid/Bound & Bound-Fluid (forces and DEM).
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity
  ,bool shift,bool mdbc2
  ,bool shiftadv,bool aleform,bool ncpress> //<vs_advshift>
  void JSphCpu::Interaction_ForcesCpuT(const stinterparmsc& t,StInterResultc& res)const
{
  float viscdt=res.viscdt;
  if(t.npf){
    //-Interaction Fluid-Fluid.
    InteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,mdbc2,shiftadv,aleform,ncpress> 
      (t.npf,t.npb,false,Visco,t.divdata,t.dcell,t.spstaurho2,t.sps2strain
      ,t.pos,t.velrho,t.code,t.idp,t.press,t.dengradcorr
      ,t.boundmode,t.tangenvel,t.motionvel,t.boundnormal //<vs_m2dbc>
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs,t.nopenshift
      ,t.fstype,t.shiftvel,t.lcorr,t.fstresh,t.presssym,t.pressasym); //<vs_advshift>
    //-Interaction Fluid-Bound.
    const float viscb=Visco*ViscoBoundFactor;
    InteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,mdbc2,shiftadv,aleform,ncpress>
      (t.npf,t.npb,true ,viscb,t.divdata,t.dcell,t.spstaurho2,t.sps2strain
      ,t.pos,t.velrho,t.code,t.idp,t.press,NULL
      ,t.boundmode,t.tangenvel,t.motionvel,t.boundnormal //<vs_m2dbc>
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs,t.nopenshift
      ,t.fstype,t.shiftvel,t.lcorr,t.fstresh,t.presssym,t.pressasym); //<vs_advshift>

    //-Interaction of DEM Floating-Bound & Floating-Floating. //(DEM)
    if(UseDEM)InteractionForcesDEM(CaseNfloat,t.divdata,t.dcell
      ,RidpMot+CaseNmoving,DemData,t.pos,t.velrho,t.code,t.idp,viscdt,t.ace);

    //-Computes tau for Laminar+SPS.
    if(tvisco==VISCO_LaminarSPS)
      ComputeSpsTau(t.npf,t.npb,t.velrho,t.sps2strain,t.spstaurho2);
    
    //-Identify probably free-surface particles for advanced shifting.  //<vs_advshift>
    if(shiftadv && InterStep==INTERSTEP_SymCorrector)                   //<vs_advshift>
      ComputeFsType(t.npf,t.npb,Simulate2D,t.fstresh,t.fstype);         //<vs_advshift>
  }
  if(t.npbok){
    //-Interaction Bound-Fluid.
    InteractionForcesBound<tker,ftmode> (t.npbok,0,t.divdata,t.dcell
      ,t.pos,t.velrho,t.code,t.idp,viscdt,t.ar);
  }
  res.viscdt=viscdt;
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool mdbc2> 
  void JSphCpu::Interaction_Forces_ct6(const stinterparmsc& t,StInterResultc& res)const
{
  if(ShiftingAdv){ //<vs_advshift_ini>
    const bool aleform=ShiftingAdv->GetAleActive();
    const bool ncpress=ShiftingAdv->GetNcPress();
    if(aleform){ const bool ale=true;
      if(ncpress)Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,shift,mdbc2,true,ale,true>(t,res);
      else       Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,shift,mdbc2,true,ale,false>(t,res);
    }
    else{        const bool ale=false;
      if(ncpress)Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,shift,mdbc2,true,ale,true>(t,res);
      else       Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,shift,mdbc2,true,ale,false>(t,res);
    }
  } //<vs_advshift_end>
  else           Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,shift,mdbc2,false,false,false>(t,res);
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity> 
  void JSphCpu::Interaction_Forces_ct5(const stinterparmsc& t,StInterResultc& res)const
{
  const bool usemdbc2=(SlipMode>=SLIP_NoSlip); //<vs_m2dbc>
  if(usemdbc2){ const bool mdbc2=true;
    if(Shifting)Interaction_Forces_ct6<tker,ftmode,tvisco,tdensity,true ,mdbc2>(t,res);
    else        Interaction_Forces_ct6<tker,ftmode,tvisco,tdensity,false,mdbc2>(t,res);
  }
  else{         const bool mdbc2=false;
    if(Shifting)Interaction_Forces_ct6<tker,ftmode,tvisco,tdensity,true ,mdbc2>(t,res);
    else        Interaction_Forces_ct6<tker,ftmode,tvisco,tdensity,false,mdbc2>(t,res);
  }
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco>
  void JSphCpu::Interaction_Forces_ct4(const stinterparmsc& t,StInterResultc& res)const
{
       if(TDensity==DDT_None)    Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_None    >(t,res);
  else if(TDensity==DDT_DDT)     Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_DDT     >(t,res);
  else if(TDensity==DDT_DDT2)    Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_DDT2    >(t,res);
  else if(TDensity==DDT_DDT2Full)Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_DDT2Full>(t,res);
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> 
  void JSphCpu::Interaction_Forces_ct3(const stinterparmsc& t,StInterResultc& res)const
{
       if(TVisco==VISCO_Artificial)Interaction_Forces_ct4<tker,ftmode,VISCO_Artificial>(t,res);
  else if(TVisco==VISCO_Laminar   )Interaction_Forces_ct4<tker,ftmode,VISCO_Laminar   >(t,res);
  else if(TVisco==VISCO_LaminarSPS)Interaction_Forces_ct4<tker,ftmode,VISCO_LaminarSPS>(t,res);
}
//==============================================================================
template<TpKernel tker>
  void JSphCpu::Interaction_Forces_ct2(const stinterparmsc& t,StInterResultc& res)const
{
       if(FtMode==FTMODE_None)Interaction_Forces_ct3<tker,FTMODE_None>(t,res);
  else if(FtMode==FTMODE_Sph )Interaction_Forces_ct3<tker,FTMODE_Sph >(t,res);
  else if(FtMode==FTMODE_Ext )Interaction_Forces_ct3<tker,FTMODE_Ext >(t,res);
}
//==============================================================================
void JSphCpu::Interaction_Forces_ct(const stinterparmsc& t,StInterResultc& res)const{
       if(TKernel==KERNEL_Wendland)  Interaction_Forces_ct2<KERNEL_Wendland>(t,res);
  else if(TKernel==KERNEL_Cubic)     Interaction_Forces_ct2<KERNEL_Cubic   >(t,res);
}

//==============================================================================
/// Update pos, dcell and code to move with indicated displacement.
/// The value of outrho indicates is it outside of the density limits.
/// Check the limits in funcion of MapRealPosMin & MapRealSize that this is valid
/// for single-cpu because DomRealPos & MapRealPos are equal. For multi-cpu it will be 
/// necessary to mark the particles that leave the domain without leaving the map.
///
/// Actualiza pos, dcell y code a partir del desplazamiento indicado.
/// El valor de outrho indica si esta fuera de los limites de densidad.
/// Comprueba los limites en funcion de MapRealPosMin y MapRealSize esto es valido
/// para single-cpu pq DomRealPos y MapRealPos son iguales. Para multi-cpu seria 
/// necesario marcar las particulas q salgan del dominio sin salir del mapa.
//==============================================================================
void JSphCpu::UpdatePos(tdouble3 rpos,double movx,double movy,double movz
  ,bool outrho,unsigned p,tdouble3* pos,unsigned* cell,typecode* code)const
{
  //-Check validity of displacement. | Comprueba validez del desplazamiento.
  bool outmov=(fabs(float(movx))>MovLimit || fabs(float(movy))>MovLimit || fabs(float(movz))>MovLimit);
  //-Applies dsiplacement. | Aplica desplazamiento.
  rpos.x+=movx; rpos.y+=movy; rpos.z+=movz;
  //-Check limits of real domain. | Comprueba limites del dominio reales.
  double dx=rpos.x-MapRealPosMin.x;
  double dy=rpos.y-MapRealPosMin.y;
  double dz=rpos.z-MapRealPosMin.z;
  bool out=(dx!=dx || dy!=dy || dz!=dz || dx<0 || dy<0 || dz<0 || dx>=MapRealSize.x || dy>=MapRealSize.y || dz>=MapRealSize.z);
  //-Adjust position according to periodic conditions and compare domain limits. | Ajusta posicion segun condiciones periodicas y vuelve a comprobar los limites del dominio.
  if(PeriActive && out){
    if(PeriX){
      if(dx<0)             { dx-=PeriXinc.x; dy-=PeriXinc.y; dz-=PeriXinc.z; }
      if(dx>=MapRealSize.x){ dx+=PeriXinc.x; dy+=PeriXinc.y; dz+=PeriXinc.z; }
    }
    if(PeriY){
      if(dy<0)             { dx-=PeriYinc.x; dy-=PeriYinc.y; dz-=PeriYinc.z; }
      if(dy>=MapRealSize.y){ dx+=PeriYinc.x; dy+=PeriYinc.y; dz+=PeriYinc.z; }
    }
    if(PeriZ){
      if(dz<0)             { dx-=PeriZinc.x; dy-=PeriZinc.y; dz-=PeriZinc.z; }
      if(dz>=MapRealSize.z){ dx+=PeriZinc.x; dy+=PeriZinc.y; dz+=PeriZinc.z; }
    }
    bool outx=!PeriX && (dx<0 || dx>=MapRealSize.x);
    bool outy=!PeriY && (dy<0 || dy>=MapRealSize.y);
    bool outz=!PeriZ && (dz<0 || dz>=MapRealSize.z);
    out=(outx||outy||outz);
    rpos=TDouble3(dx,dy,dz)+MapRealPosMin;
  }
  //-Keep current position. | Guarda posicion actualizada.
  pos[p]=rpos;
  //-Keep cell and check. | Guarda celda y check.
  if(outrho || outmov || out){//-Particle out.
    typecode rcode=code[p];
    if(out)rcode=CODE_SetOutPos(rcode);
    else if(outrho)rcode=CODE_SetOutRho(rcode);
    else rcode=CODE_SetOutMov(rcode);
    code[p]=rcode;
    cell[p]=0xFFFFFFFF;
  }
  else{//-Particle in.
    if(PeriActive){
      dx=rpos.x-DomPosMin.x;
      dy=rpos.y-DomPosMin.y;
      dz=rpos.z-DomPosMin.z;
    }
    unsigned cx=unsigned(dx/Scell),cy=unsigned(dy/Scell),cz=unsigned(dz/Scell);
    cell[p]=DCEL_Cell(DomCellCode,cx,cy,cz);
  }
}


//==============================================================================
/// Calculate new values of position, velocity & density for fluid (using Verlet).
/// Calcula nuevos valores de posicion, velocidad y densidad para el fluido (usando Verlet).
//==============================================================================
void JSphCpu::ComputeVerletVarsFluid(bool shift,const tfloat3* indirvel
  ,const tfloat4* velrho1,const tfloat4* velrho2,const byte* boundmode
  ,double dt,double dt2,const float* ar,const tfloat3* ace,const tfloat4* shiftposfs 
  ,tdouble3* pos,unsigned* dcell,typecode* code,tfloat4* velrhonew)const
{
  const double dt205=0.5*dt*dt;
  const tdouble3 gravity=ToTDouble3(Gravity);
  const int pini=int(Npb),pfin=int(Np),npf=int(Np-Npb);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=pini;p<pfin;p++){
    //-Calculate density. | Calcula densidad.
    const float rhonew=float(double(velrho2[p].w)+dt2*ar[p]);
    if(!WithFloating || CODE_IsFluid(code[p])){//-Fluid Particles.
      const tdouble3 acegr=ToTDouble3(ace[p])+gravity; //-Adds gravity.
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=double(velrho1[p].x)*dt + acegr.x*dt205;
      double dy=double(velrho1[p].y)*dt + acegr.y*dt205;
      double dz=double(velrho1[p].z)*dt + acegr.z*dt205;
      if(shift){
        dx+=double(shiftposfs[p].x);
        dy+=double(shiftposfs[p].y);
        dz+=double(shiftposfs[p].z);
      }
      bool outrho=(rhonew<RhopOutMin || rhonew>RhopOutMax);
      //-Calculate velocity & density. | Calcula velocidad y densidad.
      tfloat4 rvelrhonew=TFloat4(
        float(double(velrho2[p].x) + acegr.x*dt2),
        float(double(velrho2[p].y) + acegr.y*dt2),
        float(double(velrho2[p].z) + acegr.z*dt2),
        rhonew);
      //-Restore data of inout particles.
      if(InOut && CODE_IsFluidInout(code[p])){
        outrho=false;
        rvelrhonew=velrho2[p];
        const tfloat3 vd=indirvel[CODE_GetIzoneFluidInout(code[p])];
        if(vd.x!=FLT_MAX){
          const float v=velrho1[p].x*vd.x + velrho1[p].y*vd.y + velrho1[p].z*vd.z;
          dx=double(v*vd.x)*  dt;
          dy=double(v*vd.y)*  dt;
          dz=double(v*vd.z)*  dt;
        }
        else{
          dx=double(velrho1[p].x)*  dt;
          dy=double(velrho1[p].y)*  dt;
          dz=double(velrho1[p].z)*  dt;
        }
      }
      //-Update particle data.
      UpdatePos(pos[p],dx,dy,dz,outrho,p,pos,dcell,code);
      velrhonew[p]=rvelrhonew;
    }
    else{//-Floating Particles.
      velrhonew[p]=velrho1[p];
      if(!boundmode || boundmode[p]<BMODE_MDBC2) //<vs_m2dbc>
        velrhonew[p].w=(rhonew<RhopZero? RhopZero: rhonew); //-Avoid fluid particles being absorved by floating ones. | Evita q las floating absorvan a las fluidas.
    }
  }
}

//==============================================================================
/// Calculate new values of density and set velocity=zero for cases of  
/// (fixed+moving, no floating).
///
/// Calcula nuevos valores de densidad y pone velocidad a cero para el contorno 
/// (fixed+moving, no floating).
//==============================================================================
void JSphCpu::ComputeVelrhoBound(const tfloat4* velrhoold,const byte* boundmode
  ,const float* ar,double armul,tfloat4* velrhonew)const
{
  const bool mdbc2=(SlipMode>=SLIP_NoSlip); //<vs_m2dbc>
  const int npb=int(Npb);
  if(mdbc2){ //<vs_m2dbc_ini>
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      if(boundmode[p]>=BMODE_MDBC2)velrhonew[p]=velrhoold[p]; //-For mDBC2.
      else{//-For DBC and mDBC (SLIP_Vel0).
        const float rhonew=float(double(velrhoold[p].w)+armul*ar[p]);
        velrhonew[p]=TFloat4(0,0,0,(rhonew<RhopZero? RhopZero: rhonew));//-Avoid fluid particles being absorved by boundary ones. | Evita q las boundary absorvan a las fluidas.
      }
    }
  } //<vs_m2dbc_end>
  else{ //-For DBC and mDBC (SLIP_Vel0). //<vs_m2dbc>
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      const float rhonew=float(double(velrhoold[p].w)+armul*ar[p]);
      velrhonew[p]=TFloat4(0,0,0,(rhonew<RhopZero? RhopZero: rhonew));//-Avoid fluid particles being absorved by boundary ones. | Evita q las boundary absorvan a las fluidas.
    }
  } //<vs_m2dbc>
}

//==============================================================================
/// Update of particles according to forces and dt using Verlet.
/// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphCpu::ComputeVerlet(double dt){
  Timersc->TmStart(TMC_SuComputeStep);
  const bool shift=(Shifting!=NULL);
  const tfloat3* indirvel=(InOut? InOut->GetDirVel(): NULL);
  const byte* boundmode=AC_CPTR(BoundMode_c); //<vs_m2dbc>
  VerletStep++;
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    ComputeVerletVarsFluid(shift,indirvel,Velrho_c->cptr(),VelrhoM1_c->cptr()
      ,boundmode,dt,twodt,Ar_c->cptr(),Ace_c->cptr(),ShiftPosfs_c->cptr()
      ,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr(),VelrhoM1_c->ptr());
    ComputeVelrhoBound(VelrhoM1_c->cptr(),boundmode,Ar_c->cptr()
      ,twodt,VelrhoM1_c->ptr());
  }
  else{
    ComputeVerletVarsFluid(shift,indirvel,Velrho_c->cptr(),Velrho_c->cptr()
      ,boundmode,dt,dt,Ar_c->cptr(),Ace_c->cptr(),ShiftPosfs_c->cptr()
      ,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr(),VelrhoM1_c->ptr());
    ComputeVelrhoBound(Velrho_c->cptr(),boundmode,Ar_c->cptr()
      ,dt,VelrhoM1_c->ptr());
    VerletStep=0;
  }
  //-New values are calculated en VelrhoM1_c. | Los nuevos valores se calculan en VelrhoM1_c.
  Velrho_c->SwapPtr(VelrhoM1_c);  //-Swap Velrho_c & VelrhoM1_c. | Intercambia Velrho_c y VelrhoM1_c.
  Timersc->TmStop(TMC_SuComputeStep);
}


//==============================================================================
/// Update of particles according to forces and dt using Symplectic-Predictor.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
void JSphCpu::ComputeSymplecticPre(double dt){
  Timersc->TmStart(TMC_SuComputeStep);
  const bool mdbc2=(SlipMode>=SLIP_NoSlip); //<vs_m2dbc>
  const bool shift=false; //(ShiftingMode!=SHIFT_None); //-We strongly recommend running the shifting correction only for the corrector. If you want to re-enable shifting in the predictor, change the value here to "true".
  const double dt05=dt*.5;
  const int np=int(Np);
  const int npb=int(Npb);
  const int npf=np-npb;
  //-Assign memory to PRE variables.
  PosPre_c->Reserve();
  VelrhoPre_c->Reserve();
  //-Move current data to PRE variables for calculating the new data.
  PosPre_c->SwapPtr(Pos_c);       //- PosPre_c[]    <= Pos_c[]
  VelrhoPre_c->SwapPtr(Velrho_c); //- VelrhoPre_c[] <= Velrho_c[]

  //-Calculate new density for boundary and copy velocity.
  {
    const float*   arc=Ar_c->cptr();
    const tfloat4* velrhoprec=VelrhoPre_c->cptr();
    tfloat4*       velrhoc=Velrho_c->ptr();
    const byte*    boundmode=AC_CPTR(BoundMode_c); //<vs_m2dbc>
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      tfloat4 vr=velrhoprec[p];
      if(!mdbc2 || boundmode[p]<BMODE_MDBC2){ //-For DBC or mDBC1 //<vs_m2dbc>
        vr.w=float(double(vr.w)+dt05*arc[p]);
        vr.w=(vr.w<RhopZero? RhopZero: vr.w); //-To prevent absorption of fluid particles by boundaries. | Evita que las boundary absorvan a las fluidas.
      }
      velrhoc[p]=vr;
    }
  }

  //-Compute displacement, velocity and density for fluid.
  acdouble3 mov_c("movc",Arrays_Cpu,true);
  {
    const tfloat3*  indirvel=(InOut? InOut->GetDirVel(): NULL);
    const typecode* codec=Code_c->cptr();
    const float*    arc=Ar_c->cptr();
    const tfloat3*  acec=Ace_c->cptr();
    const tfloat4*  velrhoprec=VelrhoPre_c->cptr();
    const tfloat4*  shiftposfc=ShiftPosfs_c->cptr();
    typecode*       codec2=Code_c->ptr();
    tdouble3*       movc=mov_c.ptr();
    tfloat4*        velrhoc=Velrho_c->ptr();
    const byte*     boundmode=AC_CPTR(BoundMode_c); //<vs_m2dbc>
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=npb;p<np;p++){
      const typecode rcode=codec[p];
      //-Calculate density.
      const float rhonew=float(double(velrhoprec[p].w)+dt05*arc[p]);
      if(!WithFloating || CODE_IsFluid(rcode)){//-Fluid Particles.
        //-Calculate displacement. | Calcula desplazamiento.
        double dx=double(velrhoprec[p].x)*dt05;
        double dy=double(velrhoprec[p].y)*dt05;
        double dz=double(velrhoprec[p].z)*dt05;
        if(shift){
          dx+=double(shiftposfc[p].x);
          dy+=double(shiftposfc[p].y);
          dz+=double(shiftposfc[p].z);
        }
        bool outrho=(rhonew<RhopOutMin || rhonew>RhopOutMax);
        //-Calculate velocity & density. | Calcula velocidad y densidad.
        tfloat4 rvelrhonew=TFloat4(
          float(double(velrhoprec[p].x) + (double(acec[p].x)+Gravity.x)*  dt05),
          float(double(velrhoprec[p].y) + (double(acec[p].y)+Gravity.y)*  dt05),
          float(double(velrhoprec[p].z) + (double(acec[p].z)+Gravity.z)*  dt05),
          rhonew);
        //-Restore data of inout particles.
        if(InOut && CODE_IsFluidInout(rcode)){
          outrho=false;
          rvelrhonew=velrhoprec[p];
          const tfloat3 vd=indirvel[CODE_GetIzoneFluidInout(rcode)];
          if(vd.x!=FLT_MAX){
            const float v=rvelrhonew.x*vd.x + rvelrhonew.y*vd.y + rvelrhonew.z*vd.z;
            dx=double(v*vd.x)*  dt05;
            dy=double(v*vd.y)*  dt05;
            dz=double(v*vd.z)*  dt05;
          }
        }
        //-Update particle data.
        movc[p]=TDouble3(dx,dy,dz);
        velrhoc[p]=rvelrhonew;
        if(outrho && CODE_IsNormal(rcode))codec2[p]=CODE_SetOutRho(rcode); //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
      }
      else{//-Floating Particles.
        velrhoc[p]=velrhoprec[p];
        if(!mdbc2 || boundmode[p]<BMODE_MDBC2) //<vs_m2dbc>
          velrhoc[p].w=(rhonew<RhopZero? RhopZero: rhonew); //-Avoid fluid particles being absorbed by floating ones. | Evita q las floating absorvan a las fluidas.
      }
    }
  }

  //-Applies displacement to non-periodic fluid particles.
  {
    const tdouble3* posprec=PosPre_c->cptr();
    const tdouble3* movc=mov_c.ptr();
    typecode*       codec=Code_c->ptr();
    tdouble3*       posc=Pos_c->ptr();
    unsigned*       dcellc=Dcell_c->ptr();
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=npb;p<np;p++){
      const typecode rcode=codec[p];
      const bool outrho=CODE_IsOutRho(rcode);
      const bool normal=(!PeriActive || outrho || CODE_IsNormal(rcode));
      if(normal){//-Does not apply to periodic particles. | No se aplica a particulas periodicas
        if(CODE_IsFluid(rcode)){//-Only applied for fluid displacement. | Solo se aplica desplazamiento al fluido.
          UpdatePos(posprec[p],movc[p].x,movc[p].y,movc[p].z,outrho,p,posc,dcellc,codec);
        }
        else posc[p]=posprec[p]; //-Copy position of floating particles.
      }
    }
  }
  //-Copy previous position of the boundary particles.
  Pos_c->CopyFrom(PosPre_c,Npb);
  Timersc->TmStop(TMC_SuComputeStep);
}

//==============================================================================
/// Update particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphCpu::ComputeSymplecticCorr(double dt){
  Timersc->TmStart(TMC_SuComputeStep);
  const bool mdbc2=(SlipMode>=SLIP_NoSlip); //<vs_m2dbc>
  const bool shift=(Shifting!=NULL);
  const bool shiftadv=(ShiftingAdv!=NULL); //<vs_advshift>
  const double dt05=dt*.5;
  const int np=int(Np);
  const int npb=int(Npb);
  const int npf=np-npb;
  
  //-Calculate rho of boundary and set velocity=0. | Calcula rho de contorno y vel igual a cero.
  {
    const float*   arc=Ar_c->cptr();
    const tfloat4* velrhoprec=VelrhoPre_c->cptr();
    tfloat4*       velrhoc=Velrho_c->ptr();
    const byte*    boundmode=AC_CPTR(BoundMode_c); //<vs_m2dbc>
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      if(!mdbc2 || boundmode[p]<BMODE_MDBC2){ //<vs_m2dbc>
        const double epsilon_rdot=(-double(arc[p])/double(velrhoc[p].w))*dt;
        const float rhonew=float(double(velrhoprec[p].w)*  (2.-epsilon_rdot)/(2.+epsilon_rdot));
        velrhoc[p]=TFloat4(0,0,0,(rhonew<RhopZero? RhopZero: rhonew));//-Avoid fluid particles being absorbed by boundary ones. | Evita q las boundary absorvan a las fluidas.
      } //<vs_m2dbc>
      else velrhoc[p]=velrhoprec[p]; //-Like in GPU implementation.
    }
  }

  //-Compute displacement, velocity and density for fluid.
  acdouble3 mov_c("movc",Arrays_Cpu,true);
  {
    const tfloat3*  indirvel=(InOut? InOut->GetDirVel(): NULL);
    const typecode* codec=Code_c->cptr();
    const float*    arc=Ar_c->cptr();
    const tfloat3*  acec=Ace_c->cptr();
    const tfloat4*  velrhoprec=VelrhoPre_c->cptr();
    const tfloat4*  shiftposfc=ShiftPosfs_c->cptr();
    const tfloat4*  shiftvel=AC_CPTR(ShiftVel_c); //<vs_advshift>
    typecode*       codec2=Code_c->ptr();
    tdouble3*       movc=mov_c.ptr();
    tfloat4*        velrhoc=Velrho_c->ptr();
    const byte*     boundmode=AC_CPTR(BoundMode_c); //<vs_m2dbc>
    const tfloat4*  nopenshift=NoPenShift_c->cptr(); //<vs_m2dbcNP> SHABA 
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=npb;p<np;p++){
      const typecode rcode=codec[p];
      const double epsilon_rdot=(-double(arc[p])/double(velrhoc[p].w))*dt;
      const float rhonew=float(double(velrhoprec[p].w)*  (2.-epsilon_rdot)/(2.+epsilon_rdot));
      if(!WithFloating || CODE_IsFluid(rcode)){//-Fluid Particles.
        //-Calculate velocity & density. | Calcula velocidad y densidad.
        tfloat4 rvelrhonew=TFloat4(
          float(double(velrhoprec[p].x) + (double(acec[p].x)+Gravity.x)*  dt), 
          float(double(velrhoprec[p].y) + (double(acec[p].y)+Gravity.y)*  dt), 
          float(double(velrhoprec[p].z) + (double(acec[p].z)+Gravity.z)*  dt),
          rhonew);
        //-Calculate displacement. | Calcula desplazamiento.
        double dx=(double(velrhoprec[p].x)+double(rvelrhonew.x))*  dt05; 
        double dy=(double(velrhoprec[p].y)+double(rvelrhonew.y))*  dt05; 
        double dz=(double(velrhoprec[p].z)+double(rvelrhonew.z))*  dt05;
        //-Adding no-penetration correction velocity SHABA
        if (mdbc2) {//<vs_m2dbcNP_start>
            if (nopenshift[p].w > 5.f) { //-check if correction should be applied or not
                //-Correcting velocity
                rvelrhonew.x = velrhoprec[p].x + nopenshift[p].x;
                rvelrhonew.y = velrhoprec[p].y + nopenshift[p].y;
                rvelrhonew.z = velrhoprec[p].z + nopenshift[p].z;
                // Adding displacement
                dx = double(rvelrhonew.x) * dt;
                dy = double(rvelrhonew.y) * dt;
                dz = double(rvelrhonew.z) * dt;
            }
        }//<vs_m2dbcNP_end>
        if(shift){
          dx+=double(shiftposfc[p].x);
          dy+=double(shiftposfc[p].y);
          dz+=double(shiftposfc[p].z);
        }
        if(shiftadv){ //<vs_advshift_ini>
          dx+=double(shiftvel[p].x)*dt;
          dy+=double(shiftvel[p].y)*dt;
          dz+=double(shiftvel[p].z)*dt;
        } //<vs_advshift_end>
        bool outrho=(rhonew<RhopOutMin || rhonew>RhopOutMax);
        //-Restore data of inout particles.
        if(InOut && CODE_IsFluidInout(rcode)){
          outrho=false;
          rvelrhonew=velrhoprec[p];
          const tfloat3 vd=indirvel[CODE_GetIzoneFluidInout(rcode)];
          if(vd.x!=FLT_MAX){
            const float v=rvelrhonew.x*vd.x + rvelrhonew.y*vd.y + rvelrhonew.z*vd.z;
            dx=double(v*vd.x)*  dt;
            dy=double(v*vd.y)*  dt;
            dz=double(v*vd.z)*  dt;
          }
          else{
            dx=double(rvelrhonew.x)*  dt; 
            dy=double(rvelrhonew.y)*  dt; 
            dz=double(rvelrhonew.z)*  dt;
          }
        }
        //-Update particle data.
        movc[p]=TDouble3(dx,dy,dz);
        if(outrho && CODE_IsNormal(rcode))codec2[p]=CODE_SetOutRho(rcode); //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
        velrhoc[p]=rvelrhonew;
      }
      else{//-Floating Particles.
        velrhoc[p]=velrhoprec[p];
        if(!mdbc2 || boundmode[p]<BMODE_MDBC2) //<vs_m2dbc>
          velrhoc[p].w=(rhonew<RhopZero? RhopZero: rhonew); //-Avoid fluid particles being absorbed by floating ones. | Evita q las floating absorvan a las fluidas.
      }
    }
  }

  //-Applies displacement to non-periodic fluid particles.
  {
    const tdouble3* posprec=PosPre_c->cptr();
    const tdouble3* movc=mov_c.ptr();
    typecode*       codec=Code_c->ptr();
    tdouble3*       posc=Pos_c->ptr();
    unsigned*       dcellc=Dcell_c->ptr();
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=npb;p<np;p++){
      const typecode rcode=codec[p];
      const bool outrho=CODE_IsOutRho(rcode);
      const bool normal=(!PeriActive || outrho || CODE_IsNormal(rcode));
      if(normal){//-Does not apply to periodic particles. | No se aplica a particulas periodicas
        if(CODE_IsFluid(rcode)){//-Only applied for fluid displacement. | Solo se aplica desplazamiento al fluido.
          UpdatePos(posprec[p],movc[p].x,movc[p].y,movc[p].z,outrho,p,posc,dcellc,codec);
        }
        else posc[p]=posprec[p]; //-Copy position of floating particles.
      }
    }
  }

  //-Free memory assigned to PRE variables in ComputeSymplecticPre().
  PosPre_c->Free();
  VelrhoPre_c->Free();
  Timersc->TmStop(TMC_SuComputeStep);
}

//==============================================================================
/// Calculate variable Dt.
/// Calcula un Dt variable.
//==============================================================================
double JSphCpu::DtVariable(bool final){
  //-dt1 depends on force per unit mass.
  const double dt1=(AceMax? (sqrt(double(KernelH)/AceMax)): DBL_MAX); 
  //-dt2 combines the Courant and the viscous time-step controls.
  const double dt2=double(KernelH)/(max(Cs0,VelMax*10.)+double(KernelH)*ViscDtMax);
  //-dt new value of time step.
  double dt=CFLnumber*min(dt1,dt2);
  //-Use dt value defined by FixedDt object.
  if(FixedDt)dt=FixedDt->GetDt(TimeStep,dt);
  //-Check that the dt is valid.
  if(fun::IsNAN(dt) || fun::IsInfinity(dt))Run_Exceptioon(fun::PrintStr(
    "The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f) is NaN or infinity at nstep=%u."
    ,dt,AceMax,VelMax,ViscDtMax,Nstep));
  //-Correct dt with minimum allowed value.
  if(dt<double(DtMin)){ 
    dt=double(DtMin); DtModif++;
    if(DtModif>=DtModifWrn){
      Log->PrintfWarning("%d DTs adjusted to DtMin (t:%g, nstep:%u)",DtModif,TimeStep,Nstep);
      DtModifWrn*=10;
    }
  }
  //-Saves information about dt.
  if(final){
    if(PartDtMin>dt)PartDtMin=dt;
    if(PartDtMax<dt)PartDtMax=dt;
    //-Saves detailed information about dt in SaveDt object.
    if(SaveDt)SaveDt->AddValues(TimeStep,dt,dt1*CFLnumber,dt2*CFLnumber
      ,AceMax,ViscDtMax,VelMax);
  }
  return(dt);
}

//==============================================================================
/// Calculate final Shifting for particles' position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphCpu::RunShifting(double dt){
  Timersc->TmStart(TMC_SuShifting);
  Shifting->RunCpu(Np-Npb,Npb,dt,Velrho_c->cptr(),ShiftPosfs_c->ptr());
  Timersc->TmStop(TMC_SuShifting);
}

//==============================================================================
/// Calculate position of particles according to idp[]. When it is not met set as UINT_MAX.
/// When periactive is False assume that there are no duplicate particles (periodic ones)
/// and all are set as CODE_NORMAL.
///
/// Calcula posicion de particulas segun idp[]. Cuando no la encuentra es UINT_MAX.
/// Cuando periactive es False supone que no hay particulas duplicadas (periodicas)
/// y todas son CODE_NORMAL.
//==============================================================================
void JSphCpu::CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini
  ,unsigned idfin,const typecode* code,const unsigned* idp,unsigned* ridp)const
{
  //-Assign values UINT_MAX. | Asigna valores UINT_MAX.
  const unsigned nsel=idfin-idini;
  memset(ridp,255,sizeof(unsigned)*nsel); 
  //-Calculate position according to id. | Calcula posicion segun id.
  const int pfin=int(pini+np);
  if(periactive){//-Calculate position according to id checking that the particles are normal (i.e. not periodic). | Calcula posicion segun id comprobando que las particulas son normales (no periodicas).
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(pfin>OMP_LIMIT_COMPUTELIGHT)
    #endif
    for(int p=int(pini);p<pfin;p++){
      const unsigned id=idp[p];
      if(idini<=id && id<idfin){
        if(CODE_IsNormal(code[p]))ridp[id-idini]=p;
      }
    }
  }
  else{//-Calculate position according to id assuming that all the particles are normal (i.e. not periodic). | Calcula posicion segun id suponiendo que todas las particulas son normales (no periodicas).
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(pfin>OMP_LIMIT_COMPUTELIGHT)
    #endif
    for(int p=int(pini);p<pfin;p++){
      const unsigned id=idp[p];
      if(idini<=id && id<idfin)ridp[id-idini]=p;
    }
  }
}

//==============================================================================
/// Applies a linear movement to a group of particles.
/// Aplica un movimiento lineal a un conjunto de particulas.
//==============================================================================
void JSphCpu::MoveLinBound(unsigned np,unsigned ini,const tdouble3& mvpos
  ,const tfloat3& mvvel,const unsigned* ridpmot,tdouble3* pos,unsigned* dcell
  ,tfloat4* velrho,typecode* code)const
{
  const unsigned fin=ini+np;
  for(unsigned id=ini;id<fin;id++){
    const unsigned pid=ridpmot[id];
    if(pid!=UINT_MAX){
      UpdatePos(pos[pid],mvpos.x,mvpos.y,mvpos.z,false,pid,pos,dcell,code);
      velrho[pid].x=mvvel.x;  velrho[pid].y=mvvel.y;  velrho[pid].z=mvvel.z;
    }
  }
}

//==============================================================================
/// Applies a matrix movement to a group of particles.
/// Aplica un movimiento matricial a un conjunto de particulas.
//==============================================================================
void JSphCpu::MoveMatBound(unsigned np,unsigned ini,tmatrix4d m,double dt
  ,const unsigned* ridpmot,tdouble3* pos,unsigned* dcell,tfloat4* velrho
  ,typecode* code,tfloat3* boundnor)const
{
  const unsigned fin=ini+np;
  for(unsigned id=ini;id<fin;id++){
    const unsigned pid=ridpmot[id];
    if(pid!=UINT_MAX){
      const tdouble3 ps=pos[pid];
      tdouble3 ps2=MatrixMulPoint(m,ps);
      if(Simulate2D)ps2.y=ps.y;
      const double dx=ps2.x-ps.x, dy=ps2.y-ps.y, dz=ps2.z-ps.z;
      UpdatePos(ps,dx,dy,dz,false,pid,pos,dcell,code);
      velrho[pid].x=float(dx/dt);  velrho[pid].y=float(dy/dt);  velrho[pid].z=float(dz/dt);
      //-Computes normal.
      if(boundnor){
        const tdouble3 gs=ps+ToTDouble3(boundnor[pid]);
        const tdouble3 gs2=MatrixMulPoint(m,gs);
        boundnor[pid]=ToTFloat3(gs2-ps2);
      }
    }
  }
}

//==============================================================================
/// Calculates predefined movement of boundary particles.
/// Calcula movimiento predefinido de boundary particles.
//==============================================================================
void JSphCpu::CalcMotion(double stepdt){
  Timersc->TmStart(TMC_SuMotion);
  JSph::CalcMotion(stepdt);
  Timersc->TmStop(TMC_SuMotion);
}

//==============================================================================
/// Process movement of boundary particles.
/// Procesa movimiento de boundary particles.
//==============================================================================
void JSphCpu::RunMotion(double stepdt){
  Timersc->TmStart(TMC_SuMotion);
  BoundChanged=false;
  //-Add motion from automatic wave generation.
  if(WaveGen)CalcMotionWaveGen(stepdt);
  //-Process particles motion.
  if(DsMotion->GetActiveMotion()){
    BoundChanged=true;
    const unsigned nref=DsMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      if(m.type==MOTT_Linear){//-Linear movement.
        MoveLinBound(m.count,m.idbegin-CaseNfixed,m.linmov,ToTFloat3(m.linvel)
          ,RidpMot,Pos_c->ptr(),Dcell_c->ptr(),Velrho_c->ptr(),Code_c->ptr());
      }
      if(m.type==MOTT_Matrix){//-Matrix movement (for rotations).
        MoveMatBound(m.count,m.idbegin-CaseNfixed,m.matmov,stepdt,RidpMot
          ,Pos_c->ptr(),Dcell_c->ptr(),Velrho_c->ptr(),Code_c->ptr()
          ,AC_PTR(BoundNor_c)); 
      }      
    }
  }
  //-Management of Multi-Layer Pistons.
  if(MLPistons){
    BoundChanged=true;
    if(MLPistons->GetPiston1dCount()){//-Process motion for pistons 1D.
      MLPistons->CalculateMotion1d(TimeStep+MLPistons->GetTimeMod()+stepdt);
      MovePiston1d(CaseNmoving,0,MLPistons->GetPoszMin(),MLPistons->GetPoszCount()
        ,MLPistons->GetPistonId(),MLPistons->GetMovx(),MLPistons->GetVelx()
        ,RidpMot,Pos_c->ptr(),Dcell_c->ptr(),Velrho_c->ptr(),Code_c->ptr());
    }
    for(unsigned cp=0;cp<MLPistons->GetPiston2dCount();cp++){//-Process motion for pistons 2D.
      const JMLPistons::StMotionInfoPiston2D mot=MLPistons->CalculateMotion2d(cp
        ,TimeStep+MLPistons->GetTimeMod()+stepdt);
      MovePiston2d(mot.np,mot.idbegin-CaseNfixed,mot.posymin,mot.poszmin
        ,mot.poszcount,mot.movyz,mot.velyz,RidpMot,Pos_c->ptr(),Dcell_c->ptr()
        ,Velrho_c->ptr(),Code_c->ptr());
    }
  }
  //<vs_m2dbc_ini>
  //-Copy motion velocity and compute acceleration of moving particles.
  if(MotionVel_c){
    CopyMotionVelAce(CaseNmoving,stepdt,RidpMot,Velrho_c->cptr()
      ,MotionVel_c->ptr(),MotionAce_c->ptr());
  } //<vs_m2dbc_end>
  Timersc->TmStop(TMC_SuMotion);
}

//==============================================================================
/// Applies movement and velocity of piston 1D to a group of particles.
/// Aplica movimiento y velocidad de piston 1D a conjunto de particulas.
//==============================================================================
void JSphCpu::MovePiston1d(unsigned np,unsigned ini
  ,double poszmin,unsigned poszcount,const byte* pistonid,const double* movx
  ,const double* velx,const unsigned* ridpmot,tdouble3* pos,unsigned* dcell
  ,tfloat4* velrho,typecode* code)const
{
  const int fin=int(ini+np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(fin>OMP_LIMIT_LIGHT)
  #endif
  for(int id=int(ini);id<fin;id++){
    const unsigned pid=ridpmot[id];
    if(pid!=UINT_MAX){
      const unsigned pisid=pistonid[CODE_GetTypeValue(code[pid])];
      if(pisid<255){
        const unsigned cz=unsigned((pos[pid].z-poszmin)/Dp);
        const double rmovx=(cz<poszcount? movx[pisid*poszcount+cz]: 0);
        const float rvelx=float(cz<poszcount? velx[pisid*poszcount+cz]: 0);
        //-Updates position.
        UpdatePos(pos[pid],rmovx,0,0,false,pid,pos,dcell,code);
        //-Updates velocity.
        velrho[pid].x=rvelx;
      }
    }
  }
}
//==============================================================================
/// Applies movement and velocity of piston 2D to a group of particles.
/// Aplica movimiento y velocidad de piston 2D a conjunto de particulas.
//==============================================================================
void JSphCpu::MovePiston2d(unsigned np,unsigned ini
  ,double posymin,double poszmin,unsigned poszcount,const double* movx
  ,const double* velx,const unsigned* ridpmot,tdouble3* pos,unsigned* dcell
  ,tfloat4* velrho,typecode* code)const
{
  const int fin=int(ini+np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(fin>OMP_LIMIT_LIGHT)
  #endif
  for(int id=int(ini);id<fin;id++){
    const unsigned pid=ridpmot[id];
    if(pid!=UINT_MAX){
      const tdouble3 ps=pos[pid];
      const unsigned cy=unsigned((ps.y-posymin)/Dp);
      const unsigned cz=unsigned((ps.z-poszmin)/Dp);
      const double rmovx=(cz<poszcount? movx[cy*poszcount+cz]: 0);
      const float rvelx=float(cz<poszcount? velx[cy*poszcount+cz]: 0);
      //-Updates position.
      UpdatePos(ps,rmovx,0,0,false,pid,pos,dcell,code);
      //-Updates velocity.
      velrho[pid].x=rvelx;
    }
  }
}

//==============================================================================
/// Applies RelaxZone to selected particles.
/// Aplica RelaxZone a las particulas indicadas.
//==============================================================================
void JSphCpu::RunRelaxZone(double dt){
  Timersc->TmStart(TMC_SuMotion);
  acbyte  rzid("rzid",Arrays_Cpu,false);
  acfloat rzfactor("rzfactor",Arrays_Cpu,false);
  RelaxZones->SetFluidVel(TimeStep,dt,Np-Npb,Npb,Pos_c->cptr(),Idp_c->cptr()
    ,Velrho_c->ptr(),rzid.ptr(),rzfactor.ptr());
  Timersc->TmStop(TMC_SuMotion);
}

//==============================================================================
/// Applies Damping to selected particles.
/// Aplica Damping a las particulas indicadas.
//==============================================================================
void JSphCpu::RunDamping(double dt){
  const typecode* codeptr=(CaseNfloat || PeriActive? Code_c->cptr(): NULL);
  Damping->ComputeDampingCpu(TimeStep,dt,Np-Npb,Npb,Pos_c->cptr()
    ,codeptr,Velrho_c->ptr());
}



