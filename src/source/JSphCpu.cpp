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

#include <climits>
#include <numeric>

using namespace std;

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

  BoundNor_c=NULL;    //-mDBC
  MotionVel_c=NULL;   //-mDBC2  //<vs_m2dbc>
  MotionAce_c=NULL;   //-mDBC2  //<vs_m2dbc>
  BoundOnOff_c=NULL;  //-mDBC2  //<vs_m2dbc>

  VelrhoM1_c=NULL;    //-Verlet
  PosPre_c=NULL;      //-Symplectic
  VelrhoPre_c=NULL;   //-Symplectic

  Ace_c=NULL;
  Ar_c=NULL;
  Press_c=NULL;
  Delta_c=NULL;
  ShiftPosfs_c=NULL;  //-Shifting.

  SpsTauRho2_c=NULL;  //-Laminar+SPS.
  Sps2Strain_c=NULL;  //-Laminar+SPS.

  //<vs_flexstruc_ini>
  FlexStrucDatac=NULL;
  FlexStrucRidpc=NULL;
  Pos0c=NULL;
  NumPairsc=NULL;
  PairIdxBufferc=NULL;
  PairIdxc=NULL;
  KerCorrc=NULL;
  DefGradc=NULL;
  BoundNor0c=NULL;
  FlexStrucDtMax=0;
  //<vs_flexstruc_end>

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
  //<vs_flexstruc_ini>
  delete[] FlexStrucDatac;  FlexStrucDatac=NULL;
  delete[] FlexStrucRidpc;  FlexStrucRidpc=NULL;
  delete[] Pos0c;           Pos0c=NULL;
  delete[] NumPairsc;       NumPairsc=NULL;
  delete[] PairIdxBufferc;  PairIdxBufferc=NULL;
  delete[] PairIdxc;        PairIdxc=NULL;
  delete[] KerCorrc;        KerCorrc=NULL;
  delete[] DefGradc;        DefGradc=NULL;
  delete[] BoundNor0c;      BoundNor0c=NULL;
  //<vs_flexstruc_end>
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

  delete BoundNor_c;    BoundNor_c=NULL;    //-mDBC
  delete MotionVel_c;   MotionVel_c=NULL;   //-mDBC2  //<vs_m2dbc>
  delete MotionAce_c;   MotionAce_c=NULL;   //-mDBC2  //<vs_m2dbc>
  delete BoundOnOff_c;  BoundOnOff_c=NULL;  //-mDBC2  //<vs_m2dbc>

  delete VelrhoM1_c;    VelrhoM1_c=NULL;    //-Verlet
  delete PosPre_c;      PosPre_c=NULL;      //-Symplectic
  delete VelrhoPre_c;   VelrhoPre_c=NULL;   //-Symplectic

  delete Ace_c;         Ace_c=NULL;
  delete Ar_c;          Ar_c=NULL;
  delete Press_c;       Press_c=NULL;
  delete Delta_c;       Delta_c=NULL;
  delete ShiftPosfs_c;  ShiftPosfs_c=NULL;  //-Shifting.

  delete SpsTauRho2_c;  SpsTauRho2_c=NULL;  //-Laminar+SPS.
  delete Sps2Strain_c;  Sps2Strain_c=NULL;  //-Laminar+SPS.
  
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
    if(SlipMode!=SLIP_Vel0){ //<vs_m2dbc_ini>
      MotionVel_c =new acfloat3("MotionVelc" ,Arrays_Cpu,true);
      MotionAce_c =new acfloat3("MotionAcec" ,Arrays_Cpu,true);
      BoundOnOff_c=new acfloat( "BoundOnOffc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
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
  //-Arrays for Shifting.
  ShiftPosfs_c=new acfloat4("ShiftPosfsc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
  //-Arrays for Laminar+SPS.
  if(TVisco==VISCO_LaminarSPS){
    SpsTauRho2_c=new acsymatrix3f("SpsTauRho2c",Arrays_Cpu,true);
    Sps2Strain_c=new acsymatrix3f("Sps2Strainc",Arrays_Cpu,false); //-NO INITIAL MEMORY.
  }


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
}

//==============================================================================
/// Prepare variables for interaction functions.
/// Prepara variables para interaccion.
//==============================================================================
void JSphCpu::PreInteraction_Forces(){
  Timersc->TmStart(TMC_CfPreForces);
  //-Assign memory.
  Ar_c->Reserve();
  Ace_c->Reserve();
  Press_c->Reserve();
  if(DDTArray)Delta_c->Reserve();
  if(Shifting)ShiftPosfs_c->Reserve();
  if(TVisco==VISCO_LaminarSPS)Sps2Strain_c->Reserve();

  //-Initialise arrays.
  const unsigned npf=Np-Npb;
  Ar_c->Memset(0,Np);                                             //Arc[]=0
  Ace_c->Memset(0,Np);                                            //Acec[]=(0)
  if(AC_CPTR(Delta_c))Delta_c->Memset(0,Np);                      //Deltac[]=0
  if(AC_CPTR(Sps2Strain_c))Sps2Strain_c->MemsetOffset(Npb,0,npf); //Sps2Strainc[]=(0).

  //-Select particles for shifting.
  if(AC_CPTR(ShiftPosfs_c))Shifting->InitCpu(npf,Npb,Pos_c->cptr()
    ,ShiftPosfs_c->ptr());

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
  if(FlexStruc)FlexStrucDtMax=0;  //<vs_flexstruc>
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
  if(BoundOnOff_c)BoundOnOff_c->Free(); //-Reserved in MdbcBoundCorrection(). //<vs_m2dbc>
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
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>
    const tfloat4 velrhop1=velrho[p1];

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],false,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
      //---------------------------------------------------------------------------------------------
      bool rsym=false; //<vs_syymmetry>
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float drx=float(posp1.x-pos[p2].x);
              float dry=float(posp1.y-pos[p2].y);
        if(rsym)    dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
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
            tfloat4 velrhop2=velrho[p2];
            if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>
            const float dvx=velrhop1.x-velrhop2.x, dvy=velrhop1.y-velrhop2.y, dvz=velrhop1.z-velrhop2.z;
            if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

            {//-Viscosity.
              const float dot=drx*dvx + dry*dvy + drz*dvz;
              const float dot_rr2=dot/(rr2+Eta2);
              visc=max(dot_rr2,visc);
            }
          }
          rsym=(rsymp1 && !rsym && float(posp1.y-dry)<=KernelSize); //<vs_syymmetry>
          if(rsym)p2--;                                             //<vs_syymmetry>
        }
        else rsym=false;                                            //<vs_syymmetry>
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
/// Perform interaction between particles: Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Realiza interaccion entre particulas: Fluid/Float-Fluid/Float or Fluid/Float-Bound
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity
  ,bool shift,bool mdbc2> void JSphCpu::InteractionForcesFluid
  (unsigned n,unsigned pinit,bool boundp2,float visco
  ,StDivDataCpu divdata,const unsigned* dcell
  ,const tsymatrix3f* tau,tsymatrix3f* two_strain
  ,const tdouble3* pos,const tfloat4* velrho,const typecode* code
  ,const unsigned* idp,const float* press,const tfloat3* dengradcorr
  ,const tfloat3* boundnor,const float* boundonoff,const tfloat3* motionvel //<vs_m2dbc>
  ,float& viscdt,float* ar,tfloat3* ace,float* delta
  ,TpShifting shiftmode,tfloat4* shiftposfs)const
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

    //-Obtain data of particle p1 in case of floating objects. | Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1=false;     //-Indicate if it is floating. | Indica si es floating.
    if(USE_FLOATING){
      ftp1=CODE_IsFloating(code[p1]);
      if(ftp1 && tdensity!=DDT_None)deltap1=FLT_MAX; //-DDT is not applied to floating particles.
      if(ftp1 && shift)shiftposfsp1.x=FLT_MAX;  //-For floating objects do not calculate shifting. | Para floatings no se calcula shifting.
    }

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];
    const tfloat3 velp1=TFloat3(velrho[p1].x,velrho[p1].y,velrho[p1].z);
    const float rhop1=velrho[p1].w;
    const float pressp1=press[p1];
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>
    
    //-Variables for Laminar+SPS.
    tsymatrix3f taup1; //-Note that taup1 is tau_a/rho_a^2.
    if(tvisco==VISCO_LaminarSPS)taup1=tau[p1];

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
      //------------------------------------------------------------------------------------------------
      bool rsym=false; //<vs_syymmetry>
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float drx=float(posp1.x-pos[p2].x);
              float dry=float(posp1.y-pos[p2].y);
        if(rsym)    dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
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

          tfloat4 velrhop2=velrho[p2];
          if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

          //<vs_m2dbc_ini>
          //-Setting up boundary normals and new boundary velocities
          tfloat3 normalp2    =TFloat3(0); //-Creating a normailsed boundary normal
          tfloat3 normalvelp2 =TFloat3(0); //-Boundary particle velocity normal to boundary used in con of mass
          tfloat3 tangentvelp2=TFloat3(0); //-Boundary particle velocity tangent to boundary used in momentu
          tfloat3 movvelp2    =motionvel[p2];
          if(mdbc2 && boundp2 && !ftp2){//-Computes velocities and change massp2 for mDBC2 (not for floating particles).
            float norm=sqrt(boundnor[p2].x*boundnor[p2].x + boundnor[p2].y*boundnor[p2].y + boundnor[p2].z*boundnor[p2].z);
            if(norm>0.f){ // mDBC
              normalp2.x=boundnor[p2].x/norm; normalp2.y=boundnor[p2].y/norm; normalp2.z=boundnor[p2].z/norm;
              // calculating boundary particle velocity components
              float veldotnorm=velrhop2.x*normalp2.x + velrhop2.y*normalp2.y + velrhop2.z*normalp2.z;
              normalvelp2.x=veldotnorm*normalp2.x; normalvelp2.y=veldotnorm*normalp2.y; normalvelp2.z=veldotnorm*normalp2.z;
              tangentvelp2.x=velrhop2.x-normalvelp2.x; tangentvelp2.y=velrhop2.y-normalvelp2.y; tangentvelp2.z=velrhop2.z-normalvelp2.z;
            }
            else{ // DBC hopefulyl this doesnt break it
              normalvelp2.x=velrhop2.x;  normalvelp2.y=velrhop2.y;  normalvelp2.z=velrhop2.z;
              tangentvelp2.x=velrhop2.x; tangentvelp2.y=velrhop2.y; tangentvelp2.z=velrhop2.z;
            }
            // changing the mass of boundary particle with boundonoff
            massp2=boundonoff[p2]*massp2; 
          }
          //<vs_m2dbc_end>

          //-Velocity derivative (Momentum equation).
          if(compute){
            const float prs=(pressp1+press[p2])/(rhop1*velrhop2.w) 
              +(tker==KERNEL_Cubic? fsph::GetKernelCubic_Tensil(CSP,rr2,rhop1,pressp1,velrhop2.w,press[p2]): 0);
            const float p_vpm=-prs*massp2;
            acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
          }

          //-Density derivative (Continuity equation).
          float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
          if(mdbc2 && boundp2 && !ftp2){ //<vs_m2dbc_ini>
            dvx=velp1.x-movvelp2.x; //-mDBC2 no slip.
            dvy=velp1.y-movvelp2.y;
            dvz=velp1.z-movvelp2.z;
          } //<vs_m2dbc_end>
          
          if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz)*(rhop1/velrhop2.w);

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

          //===== Viscosity ===== 
          if(compute){
            if(mdbc2 && boundp2 && !ftp2){ //<vs_m2dbc_ini>
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
          rsym=(rsymp1 && !rsym && float(posp1.y-dry)<=KernelSize); //<vs_syymmetry>
          if(rsym)p2--;                                             //<vs_syymmetry>
        }
        else rsym=false;                                            //<vs_syymmetry>
      }
    }
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
            const float rad=sqrt(rr2);

            //-Calculate max value of demdt. | Calcula valor maximo de demdt.
            const typecode tavp2=CODE_GetTypeAndValue(code[p2]);
            const float masstotp2=demdata[tavp2].mass;
            const float taup2=demdata[tavp2].tau;
            const float kfricp2=demdata[tavp2].kfric;
            const float restitup2=demdata[tavp2].restitu;
            //const StDemData* demp2=demobjs+CODE_GetTypeAndValue(code[p2]);

            const float nu_mass=(!tpfluid? masstotp1/2: masstotp1*masstotp2/(masstotp1+masstotp2)); //-Con boundary toma la propia masa del floating 1.
            const float kn=4/(3*(taup1+taup2))*sqrt(float(Dp)/4); //-Generalized rigidity - Lemieux 2008.
            const float dvx=velrho[p1].x-velrho[p2].x, dvy=velrho[p1].y-velrho[p2].y, dvz=velrho[p1].z-velrho[p2].z; //vji
            const float nx=drx/rad, ny=dry/rad, nz=drz/rad; //normal_ji               
            const float vn=dvx*nx+dvy*ny+dvz*nz; //vji.nji
            const float demvisc=0.2f/(3.21f*(pow(nu_mass/kn,0.4f)*pow(fabs(vn),-0.2f))/40.f);
            if(demdtp1<demvisc)demdtp1=demvisc;

            const float over_lap=1.0f*float(Dp)-rad; //-(ri+rj)-|dij|
            if(over_lap>0.0f){ //-Contact.
              //-Normal.
              const float eij=(restitup1+restitup2)/2;
              const float gn=-(2.0f*log(eij)*sqrt(nu_mass*kn))/(sqrt(float(PI)+log(eij)*log(eij))); //-Generalized damping - Cummins 2010.
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
              const float ft_elast=2*(kn*float(DemDtForce)-gn)*vt/7; //-Elastic frictional string -->  ft_elast=2*(kn*fdispl-gn*vt)/7; fdispl=dtforce*vt;
              const float kfric_ij=(kfricp1+kfricp2)/2;
              float ft=kfric_ij*fn*tanh(8*vt);  //-Coulomb.
              ft=(ft<ft_elast? ft: ft_elast);   //-Not above yield criteria, visco-elastic model.
              acef=ft/ftmassp1; //-Divides by the mass of particle to obtain the acceleration.
              acep1.x+=(acef*tx); acep1.y+=(acef*ty); acep1.z+=(acef*tz);
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
  ,bool shift,bool mdbc2> void JSphCpu::Interaction_ForcesCpuT
  (const stinterparmsc& t,StInterResultc& res)const
{
  float viscdt=res.viscdt;
  if(t.npf){
    //-Interaction Fluid-Fluid.
    InteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,mdbc2> (t.npf,t.npb,false,Visco                 
      ,t.divdata,t.dcell,t.spstaurho2,t.sps2strain,t.pos,t.velrho,t.code,t.idp,t.press,t.dengradcorr
      ,t.boundnor,t.boundonoff,t.motionvel //<vs_m2dbc>
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);
    //-Interaction Fluid-Bound.
    InteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,mdbc2> (t.npf,t.npb,true ,Visco*ViscoBoundFactor
      ,t.divdata,t.dcell,t.spstaurho2,t.sps2strain,t.pos,t.velrho,t.code,t.idp,t.press,NULL
      ,t.boundnor,t.boundonoff,t.motionvel //<vs_m2dbc>
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);

    //-Interaction of DEM Floating-Bound & Floating-Floating. //(DEM)
    if(UseDEM)InteractionForcesDEM(CaseNfloat,t.divdata,t.dcell
      ,RidpMot+CaseNmoving,DemData,t.pos,t.velrho,t.code,t.idp,viscdt,t.ace);

    //-Computes tau for Laminar+SPS.
    if(tvisco==VISCO_LaminarSPS)ComputeSpsTau(t.npf,t.npb,t.velrho,t.sps2strain,t.spstaurho2);
  }
  if(t.npbok){
    //-Interaction Bound-Fluid.
    InteractionForcesBound<tker,ftmode> (t.npbok,0,t.divdata,t.dcell
      ,t.pos,t.velrho,t.code,t.idp,viscdt,t.ar);
  }
  res.viscdt=viscdt;
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity> 
  void JSphCpu::Interaction_Forces_ct5(const stinterparmsc& t,StInterResultc& res)const
{
  const bool usemdbc2=(TBoundary==BC_MDBC && SlipMode==SLIP_NoSlip); //<vs_m2dbc>
  if(usemdbc2){ const bool mdbc2=true;
    if(Shifting)Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,true ,mdbc2>(t,res);
    else        Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,false,mdbc2>(t,res);
  }
  else{         const bool mdbc2=false;
    if(Shifting)Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,true ,mdbc2>(t,res);
    else        Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,false,mdbc2>(t,res);
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
/// Perform interaction between ghost nodes of boundaries and fluid.
//==============================================================================
template<TpKernel tker,bool sim2d> void JSphCpu::InteractionMdbcCorrectionT2
  (unsigned n,StDivDataCpu divdata,float mdbcthreshold
  ,const tdouble3* pos,const typecode* code,const unsigned* idp
  ,const tfloat3* boundnor,tfloat4* velrho)
{
  const float determlimit=1e-3f;
  const int nn=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<nn;p1++)if(boundnor[p1]!=TFloat3(0)){
    float rhofinal=FLT_MAX;
    float sumwab=0;

    //-Calculates ghost node position.
    tdouble3 gposp1=pos[p1]+ToTDouble3(boundnor[p1]);
    gposp1=(PeriActive!=0? UpdatePeriodicPos(gposp1): gposp1); //-Corrected interface Position.
    //-Initializes variables for calculation.
    float rhop1=0;
    tfloat3 gradrhop1=TFloat3(0);
    tdouble3 velp1=TDouble3(0);       //-Only for velocity.
    tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
    tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(gposp1,false,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);
      //-Interaction of boundary with type Fluid/Float.
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float drx=float(gposp1.x-pos[p2].x);
        const float dry=float(gposp1.y-pos[p2].y);
        const float drz=float(gposp1.z-pos[p2].z);
        const float rr2=(drx*drx + dry*dry + drz*drz);
        if(rr2<=KernelSize2 && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
          //-Wendland kernel.
          float fac;
          const float wab=fsph::GetKernel_WabFac<tker>(CSP,rr2,fac);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass and volume of particle p2 =====
          const tfloat4 velrhop2=velrho[p2];
          const float massp2=MassFluid;
          const float volp2=massp2/velrhop2.w;

          //===== Density and its gradient =====
          rhop1+=massp2*wab;
          gradrhop1.x+=massp2*frx;
          gradrhop1.y+=massp2*fry;
          gradrhop1.z+=massp2*frz;

          //===== Kernel values multiplied by volume =====
          const float vwab=wab*volp2;
          sumwab+=vwab;
          const float vfrx=frx*volp2;
          const float vfry=fry*volp2;
          const float vfrz=frz*volp2;

          //===== Matrix A for correction =====
          if(sim2d){
            a_corr2.a11+=vwab;  a_corr2.a12+=drx*vwab;  a_corr2.a13+=drz*vwab;
            a_corr2.a21+=vfrx;  a_corr2.a22+=drx*vfrx;  a_corr2.a23+=drz*vfrx;
            a_corr2.a31+=vfrz;  a_corr2.a32+=drx*vfrz;  a_corr2.a33+=drz*vfrz;
          }
          else{
            a_corr3.a11+=vwab;  a_corr3.a12+=drx*vwab;  a_corr3.a13+=dry*vwab;  a_corr3.a14+=drz*vwab;
            a_corr3.a21+=vfrx;  a_corr3.a22+=drx*vfrx;  a_corr3.a23+=dry*vfrx;  a_corr3.a24+=drz*vfrx;
            a_corr3.a31+=vfry;  a_corr3.a32+=drx*vfry;  a_corr3.a33+=dry*vfry;  a_corr3.a34+=drz*vfry;
            a_corr3.a41+=vfrz;  a_corr3.a42+=drx*vfrz;  a_corr3.a43+=dry*vfrz;  a_corr3.a44+=drz*vfrz;
          }
        }
      }
    }

    //-Store the results.
    //--------------------
    if(sumwab>=mdbcthreshold || (mdbcthreshold>=2 && sumwab+2>=mdbcthreshold)){
      const tfloat3 dpos=(boundnor[p1]*(-1.f)); //-Boundary particle position - ghost node position.
      if(sim2d){
        const double determ=fmath::Determinant3x3(a_corr2);
        if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
          const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
          const float rhoghost=float(invacorr2.a11*rhop1 + invacorr2.a12*gradrhop1.x + invacorr2.a13*gradrhop1.z);
          const float grx=    -float(invacorr2.a21*rhop1 + invacorr2.a22*gradrhop1.x + invacorr2.a23*gradrhop1.z);
          const float grz=    -float(invacorr2.a31*rhop1 + invacorr2.a32*gradrhop1.x + invacorr2.a33*gradrhop1.z);
          rhofinal=(rhoghost + grx*dpos.x + grz*dpos.z);
        }
        else if(a_corr2.a11>0){//-Determinant is small but a11 is nonzero (0th order).
          rhofinal=float(rhop1/a_corr2.a11);
        }
      }
      else{
        const double determ=fmath::Determinant4x4(a_corr3);
        if(fabs(determ)>=determlimit){
          const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
          const float rhoghost=float(invacorr3.a11*rhop1 + invacorr3.a12*gradrhop1.x + invacorr3.a13*gradrhop1.y + invacorr3.a14*gradrhop1.z);
          const float grx=    -float(invacorr3.a21*rhop1 + invacorr3.a22*gradrhop1.x + invacorr3.a23*gradrhop1.y + invacorr3.a24*gradrhop1.z);
          const float gry=    -float(invacorr3.a31*rhop1 + invacorr3.a32*gradrhop1.x + invacorr3.a33*gradrhop1.y + invacorr3.a34*gradrhop1.z);
          const float grz=    -float(invacorr3.a41*rhop1 + invacorr3.a42*gradrhop1.x + invacorr3.a43*gradrhop1.y + invacorr3.a44*gradrhop1.z);
          rhofinal=(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
        }
        else if(a_corr3.a11>0){//-Determinant is small but a11 is nonzero (0th order).
          rhofinal=float(rhop1/a_corr3.a11);
        }
      }
      //-Store the results.
      rhofinal=(rhofinal!=FLT_MAX? rhofinal: RhopZero);
      //-SlipMode==SLIP_Vel0 (DBC vel=0)
      velrho[p1].w=rhofinal;
    }
  }
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
 template<TpKernel tker> void JSphCpu::Interaction_MdbcCorrectionT
  (const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat3* boundnor,tfloat4* velrho)
{
  //-Interaction GhostBoundaryNodes-Fluid.
  const unsigned n=(UseNormalsFt? Np: NpbOk);
  if(Simulate2D)InteractionMdbcCorrectionT2 <tker,true >(n,divdata,MdbcThreshold,pos,code,idp,boundnor,velrho);
  else          InteractionMdbcCorrectionT2 <tker,false>(n,divdata,MdbcThreshold,pos,code,idp,boundnor,velrho);
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphCpu::Interaction_MdbcCorrection(const StDivDataCpu &divdata
  ,const tdouble3* pos,const typecode* code,const unsigned* idp
  ,const tfloat3* boundnor,tfloat4* velrho)
{
  switch(TKernel){
    case KERNEL_Cubic:       Interaction_MdbcCorrectionT <KERNEL_Cubic     > (divdata,pos,code,idp,boundnor,velrho);  break;
    case KERNEL_Wendland:    Interaction_MdbcCorrectionT <KERNEL_Wendland  > (divdata,pos,code,idp,boundnor,velrho);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

//<vs_m2dbc_ini>
//------------------------------------------------------------------------------
/// Perform Pressure cloning for mDBC2.
//------------------------------------------------------------------------------
//  const float pressfinal=KerMdbc2PressClone(sim2d,rhoghost,bnormalp1,gravity,motace,dpos);
float JSphCpu::Mdbc2PressClone(bool sim2d,const float rhoghost,tfloat3 bnormalp1
  ,const tfloat3 gravity,const tfloat3 motacep1,const tfloat3 dpos)const
{
  float pressfinal=0.f;
  if(sim2d){
    const float pghost=float(Cs0*Cs0*(rhoghost-RhopZero));
    const float norm=sqrt(bnormalp1.x*bnormalp1.x + bnormalp1.z*bnormalp1.z);
    const float normx=bnormalp1.x/norm; 
    const float normz=bnormalp1.z/norm;
    const float normpos=dpos.x*normx + dpos.z*normz;
    const tfloat3 force=TFloat3(gravity.x-motacep1.x,0,gravity.z-motacep1.z);
    const float normforce=RhopZero*(force.x*normx + force.z*normz);
    pressfinal=pghost+normforce*normpos;
  }
  else{
    const float pghost=float(Cs0*Cs0*(rhoghost-RhopZero));
    const float norm=sqrt(bnormalp1.x*bnormalp1.x + bnormalp1.y*bnormalp1.y + bnormalp1.z*bnormalp1.z);
    const float normx=bnormalp1.x/norm;
    const float normy=bnormalp1.y/norm;
    const float normz=bnormalp1.z/norm;
    const float normpos=dpos.x*normx + dpos.y*normy + dpos.z*normz;
    const tfloat3 force=gravity-motacep1;
    const float normforce=RhopZero*(force.x*normx + force.y*normy + force.z*normz);
    pressfinal=pghost+normforce*normpos;
  }
  return(pressfinal);
}
//------------------------------------------------------------------------------
/// Calculates the infinity norm of a 3x3 matrix.
//------------------------------------------------------------------------------
float JSphCpu::Mdbc2InfNorm3x3(tmatrix3d mat)const{
  const float row1=float(fabs(mat.a11) + fabs(mat.a12) + fabs(mat.a13));
  const float row2=float(fabs(mat.a21) + fabs(mat.a22) + fabs(mat.a23));
  const float row3=float(fabs(mat.a31) + fabs(mat.a32) + fabs(mat.a33));
  const float infnorm=max(row1,max(row2,row3));
  return(infnorm);
}
//------------------------------------------------------------------------------
/// Calculates the infinity norm of a 4x4 matrix.
//------------------------------------------------------------------------------
float JSphCpu::Mdbc2InfNorm4x4(tmatrix4d mat)const{
  const float row1=float(fabs(mat.a11) + fabs(mat.a12) + fabs(mat.a13) + fabs(mat.a14));
  const float row2=float(fabs(mat.a21) + fabs(mat.a22) + fabs(mat.a23) + fabs(mat.a24));
  const float row3=float(fabs(mat.a31) + fabs(mat.a32) + fabs(mat.a33) + fabs(mat.a34));
  const float row4=float(fabs(mat.a41) + fabs(mat.a42) + fabs(mat.a43) + fabs(mat.a44));
  const float infnorm=max(row1,max(row2,max(row3,row4)));
  return(infnorm);
}
//==============================================================================
/// Perform interaction between ghost nodes of boundaries and fluid.
//==============================================================================
template<TpKernel tker,bool sim2d> void JSphCpu::InteractionMdbc2CorrectionT2
  (unsigned n,const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat3* boundnor,const tfloat3* motionvel
  ,const tfloat3* motionace,tfloat4* velrho,float* boundonoff)
{
  const int nn=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<nn;p1++){
    const tfloat3 bnormalp1=boundnor[p1];
    if(bnormalp1!=TFloat3(0)){
      float rhofinal=FLT_MAX;
      tfloat3 velrhofinal=TFloat3(0);
      float sumwab=0;
      float submerged=0;                    //-Not for Vel0.
      const tfloat3 motacep1=motionace[p1]; //-Not for Vel0.

      //-Calculates ghost node position.
      tdouble3 gposp1=pos[p1]+ToTDouble3(bnormalp1);
      //-Corrected interface Position.
      gposp1=(PeriActive!=0? UpdatePeriodicPos(gposp1): gposp1);

      //-Initializes variables for calculation.
      float rhop1=0;
      tfloat3 gradrhop1=TFloat3(0);
      tdouble3 velp1=TDouble3(0);       //-Only for velocity.
      tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
      tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(gposp1,false,divdata);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);
        //-Interaction of boundary with type Fluid/Float.
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(gposp1.x-pos[p2].x);
          const float dry=float(gposp1.y-pos[p2].y);
          const float drz=float(gposp1.z-pos[p2].z);
          const float rr2=(drx*drx + dry*dry + drz*drz);
          if(rr2<=KernelSize2 && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
            //-Wendland kernel.
            float fac;
            const float wab=fsph::GetKernel_WabFac<tker>(CSP,rr2,fac);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

            //===== Get mass and volume of particle p2 =====
            const tfloat4 velrhop2=velrho[p2];
            const float massp2=MassFluid;
            const float volp2=massp2/velrhop2.w;

            //===== Check if ghost node is submerged =====
            submerged-=volp2*(drx*frx + dry*fry + drz*frz);

            //===== Density and its gradient =====
            rhop1+=massp2*wab;
            gradrhop1.x+=massp2*frx;
            gradrhop1.y+=massp2*fry;
            gradrhop1.z+=massp2*frz;

            //===== Kernel values multiplied by volume =====
            const float vwab=wab*volp2;
            sumwab+=vwab;
            const float vfrx=frx*volp2;
            const float vfry=fry*volp2;
            const float vfrz=frz*volp2;

            //===== Velocity =====
            if(SlipMode!=SLIP_Vel0){
              velp1.x+=vwab*velrhop2.x;
              velp1.y+=vwab*velrhop2.y;
              velp1.z+=vwab*velrhop2.z;
            }

            //===== Matrix A for correction =====
            if(sim2d){
              a_corr2.a11+=vwab;  a_corr2.a12-=drx*vwab;  a_corr2.a13-=drz*vwab;
              a_corr2.a21+=vfrx;  a_corr2.a22-=drx*vfrx;  a_corr2.a23-=drz*vfrx;
              a_corr2.a31+=vfrz;  a_corr2.a32-=drx*vfrz;  a_corr2.a33-=drz*vfrz;
            }
            else{
              a_corr3.a11+=vwab;  a_corr3.a12-=drx*vwab;  a_corr3.a13-=dry*vwab;  a_corr3.a14-=drz*vwab;
              a_corr3.a21+=vfrx;  a_corr3.a22-=drx*vfrx;  a_corr3.a23-=dry*vfrx;  a_corr3.a24-=drz*vfrx;
              a_corr3.a31+=vfry;  a_corr3.a32-=drx*vfry;  a_corr3.a33-=dry*vfry;  a_corr3.a34-=drz*vfry;
              a_corr3.a41+=vfrz;  a_corr3.a42-=drx*vfrz;  a_corr3.a43-=dry*vfrz;  a_corr3.a44-=drz*vfrz;
            }
          }
        }
      }

      //-Store the results.
      //--------------------
      if(submerged>0.f){
        boundonoff[p1]=1.0f;
        const tfloat3 dpos=(boundnor[p1]*(-1.f)); //-Boundary particle position - ghost node position.
        if(sim2d){//-2D simulation.
          if(sumwab<0.1f){ //-If kernel sum is small use shepherd density and vel0.
            //-Trying to avoid too negative pressures.
            const float rhoghost=max(RhopZero,float(rhop1/a_corr2.a11));
            //-Clone particle proceedure.
            const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
            //-Final values to be saved.
            rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
          }
          else{//-Chech if matrix is invertible and well conditioned.
            const double determ=fmath::Determinant3x3(a_corr2);
            if(fabs(determ)>=0.001){//-Use 1e-3f (first_order).
              const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);
              //-Calculate the scaled condition number
              const float infnorma   =Mdbc2InfNorm3x3(a_corr2);
              const float infnormainv=Mdbc2InfNorm3x3(invacorr2);
              const float condinf=float(Dp*Dp)*infnorma*infnormainv;
              if(condinf<=50){//-If matrix is well conditioned use matrix inverse for density and shepherd for velocity.
                const float rhoghost=float(invacorr2.a11*rhop1 + invacorr2.a12*gradrhop1.x + invacorr2.a13*gradrhop1.z);
                //-Clone particle proceedure.
                const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
                //-Final values to be saved.
                rhofinal=RhopZero + float(pressfinal/(Cs0*Cs0));
              }
              else{ //-If ill conditioned use shepherd.
                const float rhoghost=float(rhop1/a_corr2.a11);
                //-Clone particle proceedure.
                const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
                //-Final values to be saved.
                rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
              }
            }
            else{//-If not invertible use shepherd for density.
              const float rhoghost=float(rhop1/a_corr2.a11);
              //-Clone particle proceedure.
              const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
              //-Final values to be saved.
              rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
            }
          }
          //-velocity with Shepherd Sum.
          velrhofinal.x=float(velp1.x/a_corr2.a11);
          velrhofinal.y=0.f;
          velrhofinal.z=float(velp1.z/a_corr2.a11);
        }
        else{//-3D simulation.
          //-Density with pressure cloning.
          if(sumwab<0.1f){//-If kernel sum is small use shepherd density and vel0.
            //-Trying to avoid too negative pressures.
            const float rhoghost=max(RhopZero,float(rhop1/a_corr3.a11));
            //-Clone particle proceedure.
            const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
            //-Final values to be saved.
            rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
          }
          else{
            const double determ=fmath::Determinant4x4(a_corr3);
            if(fabs(determ)>=0.001){
              const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
              //-Calculate the scaled condition number.
              const float infnorma   =Mdbc2InfNorm4x4(a_corr3);
              const float infnormainv=Mdbc2InfNorm4x4(invacorr3);
              const float condinf=float(Dp*Dp)*infnorma*infnormainv;
              if(condinf<=50){//-If matrix is well conditioned use matrix inverse for density.
                const float rhoghost=float(invacorr3.a11*rhop1 + invacorr3.a12*gradrhop1.x + invacorr3.a13*gradrhop1.y + invacorr3.a14*gradrhop1.z);
                //-Clone particle proceedure.
                const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
                //-Final values to be saved.
                rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
              }
              else{//-Matrix is not well conditioned, use shepherd for density.
                const float rhoghost=float(rhop1/a_corr3.a11);
                //-Clone particle proceedure.
                const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
                //-Final values to be saved.
                rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
              }
            }
            else{//-Matrix is not invertible use shepherd for density.
              const float rhoghost=float(rhop1/a_corr3.a11);
              //-Clone particle proceedure.
              const float pressfinal=Mdbc2PressClone(sim2d,rhoghost,bnormalp1,Gravity,motacep1,dpos);
              //-Final values to be saved.
              rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
            }
          }
          //-Velocity with Shepherd sum.
          velrhofinal.x=float(velp1.x/a_corr3.a11);
          velrhofinal.y=float(velp1.y/a_corr3.a11);
          velrhofinal.z=float(velp1.z/a_corr3.a11);
        }

        //-Store the results.
        if(SlipMode==SLIP_NoSlip){//-No-Slip: vel = 2*motion - ghost
          const tfloat3 v=motionvel[p1];
          velrho[p1]=TFloat4(v.x+v.x - velrhofinal.x,
                             v.y+v.y - velrhofinal.y,
                             v.z+v.z - velrhofinal.z, 
                             rhofinal);
        }
        if(SlipMode==SLIP_FreeSlip){//-No-Penetration and free slip.
          tfloat3 fsvelfinal; //-Final free slip boundary velocity.
          const tfloat3 v=motionvel[p1];
          const float motion=sqrt(v.x*v.x + v.y*v.y + v.z*v.z); //-To check if boundary moving.
          const float norm=sqrt(bnormalp1.x*bnormalp1.x + bnormalp1.y*bnormalp1.y + bnormalp1.z*bnormalp1.z);
          //-Creating a normailsed boundary normal.
          const tfloat3 normal=TFloat3(fabs(bnormalp1.x)/norm,fabs(bnormalp1.y)/norm,fabs(bnormalp1.z)/norm);
          //-Finding the velocity componants normal and tangential to boundary.
          const tfloat3 normvel=TFloat3(velrhofinal.x*normal.x,velrhofinal.y*normal.y,velrhofinal.z*normal.z);//-Velocity in direction of normal pointing into fluid).
          if(motion>0){ //-If moving boundary.
            const tfloat3 normmot=TFloat3(v.x*normal.x,v.y*normal.y,v.z*normal.z); //-Boundary motion in direction normal to boundary.
            fsvelfinal=TFloat3(normmot.x+normmot.x - normvel.x,
                               normmot.y+normmot.y - normvel.y,
                               normmot.z+normmot.z - normvel.z);
            //-Only velocity in normal direction for no-penetration.
            //-Fluid sees zero velocity in the tangetial direction.
          }
          else {
            const tfloat3 tangvel=velrhofinal-normvel; //-Velocity tangential to normal.
            fsvelfinal=tangvel-normvel;
            //-Tangential velocity equal to fluid velocity for free slip.
            //-Normal velocity reversed for no-penetration.
          }
          //-Save the velocity and density.
          velrho[p1]=TFloat4(fsvelfinal.x,fsvelfinal.y,fsvelfinal.z,rhofinal);
        }
      }
      else{//-If unsubmerged switch off boundary particle.
        boundonoff[p1]=0;
        const tfloat3 v=motionvel[p1];
        velrho[p1]=TFloat4(v.x,v.y,v.z,RhopZero);
      }
    }
  }
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
 template<TpKernel tker> void JSphCpu::Interaction_Mdbc2CorrectionT
  (const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat3* boundnor,const tfloat3* motionvel
  ,const tfloat3* motionace,tfloat4* velrho,float* boundonoff)
{
  //-Interaction GhostBoundaryNodes-Fluid.
  const unsigned n=NpbOk;  //-Note that floating boidies are not supported by mDBC2.
  if(Simulate2D)InteractionMdbc2CorrectionT2 <tker,true >(n,divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundonoff);
  else          InteractionMdbc2CorrectionT2 <tker,false>(n,divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundonoff);
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphCpu::Interaction_Mdbc2Correction(const StDivDataCpu &divdata
  ,const tdouble3* pos,const typecode* code,const unsigned* idp
  ,const tfloat3* boundnor,const tfloat3* motionvel,const tfloat3* motionace
  ,tfloat4* velrho,float* boundonoff)
{
  switch(TKernel){
    case KERNEL_Cubic:       Interaction_Mdbc2CorrectionT <KERNEL_Cubic     > (divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundonoff);  break;
    case KERNEL_Wendland:    Interaction_Mdbc2CorrectionT <KERNEL_Wendland  > (divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundonoff);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}
//<vs_m2dbc_end>


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
  if(Symmetry && rpos.y<0)rpos.y=-rpos.y; //<vs_syymmetry>
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
  ,const tfloat4* velrho1,const tfloat4* velrho2,double dt,double dt2
  ,const float* ar,const tfloat3* ace,const tfloat4* shiftposfs 
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
void JSphCpu::ComputeVelrhoBound(const tfloat4* velrhoold,const float* ar
  ,double armul,tfloat4* velrhonew)const
{
  const bool mdbc2=(TBoundary==BC_MDBC && SlipMode==SLIP_NoSlip); //<vs_m2dbc>
  if(mdbc2){ //<vs_m2dbc_ini>
    const int npb=int(Npb);
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      velrhonew[p]=velrhoold[p]; //-Check...
    }
  } //<vs_m2dbc_end>
  else{ //-For DBC and mDBC (SLIP_Vel0).
    const int npb=int(Npb);
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      const float rhonew=float(double(velrhoold[p].w)+armul*ar[p]);
      velrhonew[p]=TFloat4(0,0,0,(rhonew<RhopZero? RhopZero: rhonew));//-Avoid fluid particles being absorved by boundary ones. | Evita q las boundary absorvan a las fluidas.
    }
  }
}

//==============================================================================
/// Update of particles according to forces and dt using Verlet.
/// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphCpu::ComputeVerlet(double dt){
  Timersc->TmStart(TMC_SuComputeStep);
  const bool shift=(Shifting!=NULL);
  const tfloat3* indirvel=(InOut? InOut->GetDirVel(): NULL);
  VerletStep++;
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    ComputeVerletVarsFluid(shift,indirvel,Velrho_c->cptr(),VelrhoM1_c->cptr()
      ,dt,twodt,Ar_c->cptr(),Ace_c->cptr(),ShiftPosfs_c->cptr()
      ,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr(),VelrhoM1_c->ptr());
    ComputeVelrhoBound(VelrhoM1_c->cptr(),Ar_c->cptr(),twodt,VelrhoM1_c->ptr());
  }
  else{
    ComputeVerletVarsFluid(shift,indirvel,Velrho_c->cptr(),Velrho_c->cptr()
      ,dt,dt,Ar_c->cptr(),Ace_c->cptr(),ShiftPosfs_c->cptr()
      ,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr(),VelrhoM1_c->ptr());
    ComputeVelrhoBound(Velrho_c->cptr(),Ar_c->cptr(),dt,VelrhoM1_c->ptr());
    VerletStep=0;
  }

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersc->TmStart(TMC_SuFlexStruc);
    ComputeSemiImplicitEulerFlexStruc(dt,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr());
    BoundChanged=true;
    Timersc->TmStop(TMC_SuFlexStruc);
  }
  //<vs_flexstruc_end>

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
  const bool mdbc2=(TBoundary==BC_MDBC && SlipMode==SLIP_NoSlip); //<vs_m2dbc>
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
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      const tfloat4 vr=velrhoprec[p];
      const float rhonew=float(double(vr.w)+dt05*arc[p]);
      velrhoc[p]=TFloat4(vr.x,vr.y,vr.z,(rhonew<RhopZero? RhopZero: rhonew));//-Avoid fluid particles being absorbed by boundary ones. | Evita q las boundary absorvan a las fluidas.
      if(mdbc2)velrhoc[p]=vr; //<vs_m2dbc>
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

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersc->TmStart(TMC_SuFlexStruc);
    ComputeSymplecticPreFlexStruc(dt05,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr());
    BoundChanged=true;
    Timersc->TmStop(TMC_SuFlexStruc);
  }
  //<vs_flexstruc_end>

  Timersc->TmStop(TMC_SuComputeStep);
}

//==============================================================================
/// Update particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphCpu::ComputeSymplecticCorr(double dt){
  Timersc->TmStart(TMC_SuComputeStep);
  const bool mdbc2=(TBoundary==BC_MDBC && SlipMode==SLIP_NoSlip); //<vs_m2dbc>
  const bool shift=(Shifting!=NULL);
  const double dt05=dt*.5;
  const int np=int(Np);
  const int npb=int(Npb);
  const int npf=np-npb;
  
  //-Calculate rho of boundary and set velocity=0. | Calcula rho de contorno y vel igual a cero.
  if(!mdbc2){ //<vs_m2dbc>
    const float*   arc=Ar_c->cptr();
    const tfloat4* velrhoprec=VelrhoPre_c->cptr();
    tfloat4*       velrhoc=Velrho_c->ptr();
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
    #endif
    for(int p=0;p<npb;p++){
      const double epsilon_rdot=(-double(arc[p])/double(velrhoc[p].w))*dt;
      const float rhonew=float(double(velrhoprec[p].w)*  (2.-epsilon_rdot)/(2.+epsilon_rdot));
      velrhoc[p]=TFloat4(0,0,0,(rhonew<RhopZero? RhopZero: rhonew));//-Avoid fluid particles being absorbed by boundary ones. | Evita q las boundary absorvan a las fluidas.
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
        if(shift){
          dx+=double(shiftposfc[p].x);
          dy+=double(shiftposfc[p].y);
          dz+=double(shiftposfc[p].z);
        }
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

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersc->TmStart(TMC_SuFlexStruc);
    ComputeSymplecticCorrFlexStruc(dt05,dt,Pos_c->ptr(),Dcell_c->ptr(),Code_c->ptr());
    BoundChanged=true;
    Timersc->TmStop(TMC_SuFlexStruc);
  }
  //<vs_flexstruc_end>

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
  //-dt3 uses the maximum speed of sound across all structure particles.
  const double dt3=(FlexStrucDtMax? double(KernelH)/FlexStrucDtMax: DBL_MAX); //<vs_flexstruc>
  //-dt new value of time step.
  double dt=CFLnumber*min(dt1,min(dt2,dt3));
  if(FixedDt)dt=FixedDt->GetDt(TimeStep,dt);
  if(fun::IsNAN(dt) || fun::IsInfinity(dt)){
    if(FlexStruc)Run_Exceptioon(fun::PrintStr("The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f, FlexStrucDtMax=%f) is NaN or infinity at nstep=%u.",dt,AceMax,VelMax,ViscDtMax,FlexStrucDtMax,Nstep));
    else         Run_Exceptioon(fun::PrintStr("The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f) is NaN or infinity at nstep=%u.",dt,AceMax,VelMax,ViscDtMax,Nstep));
  }
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
    if(SaveDt)SaveDt->AddValues(TimeStep,dt,dt1*CFLnumber,dt2*CFLnumber,dt3*CFLnumber,AceMax,ViscDtMax,FlexStrucDtMax,VelMax);
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
/// Copy motion velocity to MotionVel[].
/// Copia velocidad de movimiento a MotionVel[].
//==============================================================================
void JSphCpu::CopyMotionVelAce(unsigned nmoving,double dt,const unsigned* ridpmot
  ,const tfloat4* velrho,tfloat3* motionvel,tfloat3* motionace)const
{
  const int n=int(nmoving);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int id=0;id<n;id++){
    const unsigned pid=ridpmot[id];
    if(pid!=UINT_MAX){
      const tfloat3 mvel0=motionvel[pid];
      const tfloat4 v=velrho[pid];
      motionace[pid]=TFloat3(float((double(v.x)-mvel0.x)/dt),
                             float((double(v.y)-mvel0.y)/dt),
                             float((double(v.z)-mvel0.z)/dt));
      motionvel[pid]=TFloat3(v.x,v.y,v.z);
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
  //-Copy motion velocity and compute acceleration of moving particles.
  if(MotionVel_c){ //<vs_m2dbc_ini>
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

//<vs_flexstruc_ini>
//==============================================================================
/// Finds the clamp particles and updates the code.
/// Encuentra las partículas de abrazadera y actualiza el código.
//==============================================================================
void JSphCpu::SetFlexStrucClampCodes(unsigned np,const tdouble3* pos
  ,const StFlexStrucData* flexstrucdata,typecode* code)const
{
  //-Starts execution using OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<n;p++){
    const int p1=p;      //-Number of particle.
    const typecode codep1=code[p1];

    //-If potentially a clamp particle.
    if((CODE_IsFixed(codep1)||CODE_IsMoving(codep1))&&!CODE_IsFlexStrucAny(codep1)){

      //-Loads particle p1 data.
      const tdouble3 posp1=pos[p1];

      //-Loop through other boundary particles.
      for(int p2=0;p2<n;p2++){
        const typecode codep2=code[p2];
        if(CODE_IsFlexStrucFlex(codep2)){
          const unsigned nc=flexstrucdata[CODE_GetIbodyFlexStruc(codep2)].nc;
          const typecode* clampcode=flexstrucdata[CODE_GetIbodyFlexStruc(codep2)].clampcode;
          if(std::find(clampcode,clampcode+nc,codep1)!=clampcode+nc){
            const tdouble3 posp2=pos[p2];
            const float drx=float(posp1.x-posp2.x);
            const float dry=float(posp1.y-posp2.y);
            const float drz=float(posp1.z-posp2.z);
            const float rr2=drx*drx+dry*dry+drz*drz;
            if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
              code[p1]=CODE_ToFlexStrucClamp(codep1,CODE_GetIbodyFlexStruc(codep2));
              break;
            }
          }
        }
      }
    }
  }
}

//==============================================================================
/// Checks if any normals are defined for the flexible structure particles.
/// Comprueba si se han definido normales para las partículas de estructura flexible.
//==============================================================================
bool JSphCpu::FlexStrucHasNormals(unsigned npb,const typecode* code,const tfloat3* boundnor)const{
  vector<unsigned> idx(npb);
  iota(idx.begin(),idx.end(),0);
  return any_of(idx.begin(),idx.end(),[code,boundnor](unsigned i){return CODE_IsFlexStrucFlex(code[i]) && (boundnor[i].x!=0 || boundnor[i].y!=0 || boundnor[i].z!=0);});
}

//==============================================================================
/// Counts the number of flexible structure particles (includes clamps).
/// Cuenta el número de partículas de estructura flexible (incluye abrazaderas).
//==============================================================================
unsigned JSphCpu::CountFlexStrucParts(unsigned npb,const typecode* code)const{
  return(unsigned( count_if(code,code+npb,[](const typecode c){return CODE_IsFlexStrucAny(c);}) ));
}

//==============================================================================
/// Calculates indices to the main arrays for the flexible structure particles.
/// Calcula los índices de las matrices principales para las partículas de estructura flexible.
//==============================================================================
void JSphCpu::CalcFlexStrucRidp(unsigned npb,const typecode* code,unsigned* flexstrucridp)const{
  vector<unsigned> idx(npb);
  iota(idx.begin(),idx.end(),0);
  copy_if(idx.begin(),idx.end(),flexstrucridp,[code](unsigned i){return CODE_IsFlexStrucAny(code[i]);});
}

//==============================================================================
/// Gathers values from a main array into the smaller flexible structure array.
/// Reúne valores de una matriz principal en la matriz de estructura flexible más pequeña.
//==============================================================================
void JSphCpu::GatherToFlexStrucArray(unsigned npfs,const unsigned* flexstrucridp,const tdouble3* fullarray,tdouble3* flexstrucarray)const{
  //-Starts execution using OpenMP.
  const int n=int(npfs);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<n;p++)flexstrucarray[p]=fullarray[flexstrucridp[p]];
}

//==============================================================================
/// Gathers values from a main array into the smaller flexible structure array.
/// Reúne valores de una matriz principal en la matriz de estructura flexible más pequeña.
//==============================================================================
void JSphCpu::GatherToFlexStrucArray(unsigned npfs,const unsigned* flexstrucridp,const tfloat3* fullarray,tfloat3* flexstrucarray)const{
  //-Starts execution using OpenMP.
  const int n=int(npfs);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<n;p++)flexstrucarray[p]=fullarray[flexstrucridp[p]];
}

//==============================================================================
/// Counts the total number of flexible structure pairs (neighbours).
/// Cuenta el número total de pares de estructuras flexibles (vecinos).
//==============================================================================
unsigned JSphCpu::CountFlexStrucPairs(unsigned np,const tdouble3* pos0,unsigned* numpairs)const{
  //-Starts execution using OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<n;p++){
    const int pfs1=p;      //-Number of particle.
    unsigned numpairsp1=0;

    //-Loads particle p1 data.
    const tdouble3 pos0p1=pos0[pfs1];

    //-Loop through other flexible structure particles.
    for(int pfs2=0;pfs2<n;pfs2++){
      const tdouble3 pos0p2=pos0[pfs2];
      const float drx0=float(pos0p1.x-pos0p2.x);
      const float dry0=float(pos0p1.y-pos0p2.y);
      const float drz0=float(pos0p1.z-pos0p2.z);
      const float rr20=drx0*drx0+dry0*dry0+drz0*drz0;
      if(rr20<=KernelSize2&&rr20>=ALMOSTZERO)numpairsp1++;
    }
    numpairs[pfs1]=numpairsp1;
  }
  return accumulate(numpairs,numpairs+n,0);
}

//==============================================================================
/// Sets the indices for each flexible structure pair.
/// Establece los índices para cada par de estructuras flexibles.
//==============================================================================
void JSphCpu::SetFlexStrucPairs(unsigned np,const tdouble3* pos0,unsigned** pairidx)const{
  //-Starts execution using OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<n;p++){
    const int pfs1=p;      //-Number of particle.
    unsigned idx=0;

    //-Loads particle p1 data.
    const tdouble3 pos0p1=pos0[pfs1];

    //-Loop through other flexible structure particles.
    for(int pfs2=0;pfs2<n;pfs2++){
      const tdouble3 pos0p2=pos0[pfs2];
      const float drx0=float(pos0p1.x-pos0p2.x);
      const float dry0=float(pos0p1.y-pos0p2.y);
      const float drz0=float(pos0p1.z-pos0p2.z);
      const float rr20=drx0*drx0+dry0*dry0+drz0*drz0;
      if(rr20<=KernelSize2&&rr20>=ALMOSTZERO)pairidx[pfs1][idx++]=pfs2;
    }
  }
}

//==============================================================================
/// Calculates the kernel correction matrix for each flexible structure particle.
/// Calcula la matriz de corrección del kernel para cada partícula de estructura flexible.
//==============================================================================
template<TpKernel tker,bool simulate2d> void JSphCpu::CalcFlexStrucKerCorr(unsigned np,const typecode* code,const StFlexStrucData* flexstrucdata
        ,const unsigned* flexstrucridp,const tdouble3* pos0,const unsigned* numpairs,const unsigned* const* pairidx
        ,tmatrix3f* kercorr)const
{
  //-Starts execution using OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<n;p++){
    const unsigned pfs1=p;      //-Number of particle.

    //-Get number of pairs for this particle.
    const unsigned numpairsp1=numpairs[pfs1];

    //-If this particle has pairs.
    if(numpairsp1){
      tmatrix3f kercorrp1={0};

      //-Obtains basic data of particle p1.
      const float vol0p1=flexstrucdata[CODE_GetIbodyFlexStruc(code[flexstrucridp[pfs1]])].vol0;
      const tdouble3 pos0p1=pos0[pfs1];

      //-Calculate kernel correction matrix.
      for(unsigned pair=0;pair<numpairsp1;pair++){
        const unsigned pfs2=pairidx[pfs1][pair];
        const tdouble3 pos0p2=pos0[pfs2];
        const float drx0=float(pos0p1.x-pos0p2.x);
        const float dry0=float(pos0p1.y-pos0p2.y);
        const float drz0=float(pos0p1.z-pos0p2.z);
        const float rr20=drx0*drx0+dry0*dry0+drz0*drz0;
        //-Computes kernel.
        const float fac0=fsph::GetKernel_Fac<tker>(CSP,rr20);
        const float frx0=fac0*drx0,fry0=fac0*dry0,frz0=fac0*drz0; //-Gradients.
        kercorrp1.a11-=vol0p1*drx0*frx0; kercorrp1.a12-=vol0p1*drx0*fry0; kercorrp1.a13-=vol0p1*drx0*frz0;
        kercorrp1.a21-=vol0p1*dry0*frx0; kercorrp1.a22-=vol0p1*dry0*fry0; kercorrp1.a23-=vol0p1*dry0*frz0;
        kercorrp1.a31-=vol0p1*drz0*frx0; kercorrp1.a32-=vol0p1*drz0*fry0; kercorrp1.a33-=vol0p1*drz0*frz0;
      }
      if(simulate2d){
        kercorrp1.a12=kercorrp1.a21=kercorrp1.a23=kercorrp1.a32=0.0;
        kercorrp1.a22=1.0;
      }
      kercorr[pfs1]=fmath::InverseMatrix3x3(kercorrp1);
    }
  }
}

//==============================================================================
/// Calculates the kernel correction matrix for each flexible structure particle.
/// Calcula la matriz de corrección del kernel para cada partícula de estructura flexible.
//==============================================================================
template<TpKernel tker,bool simulate2d> void JSphCpu::CalcFlexStrucKerCorrT()const{
  if(CaseNflexstruc){
    CalcFlexStrucKerCorr<tker,simulate2d>
      (CaseNflexstruc,Code_c->cptr(),FlexStrucDatac,FlexStrucRidpc,Pos0c,NumPairsc,PairIdxc,KerCorrc);
  }
}

//==============================================================================
/// Calculates the kernel correction matrix for each flexible structure particle.
/// Calcula la matriz de corrección del kernel para cada partícula de estructura flexible.
//==============================================================================
template<TpKernel tker> void JSphCpu::CalcFlexStrucKerCorr_ct0()const{
  if(Simulate2D)CalcFlexStrucKerCorrT<tker,true> ();
  else          CalcFlexStrucKerCorrT<tker,false>();
}

//==============================================================================
/// Calculates the kernel correction matrix for each flexible structure particle.
/// Calcula la matriz de corrección del kernel para cada partícula de estructura flexible.
//==============================================================================
void JSphCpu::CalcFlexStrucKerCorr()const{
  if(TKernel==KERNEL_Wendland)  CalcFlexStrucKerCorr_ct0<KERNEL_Wendland>();
  else if(TKernel==KERNEL_Cubic)CalcFlexStrucKerCorr_ct0<KERNEL_Cubic>   ();
}

//==============================================================================
/// Calculates the deformation gradient matrix for each flexible structure particle.
/// Calcula la matriz de gradiente de deformación para cada partícula de estructura flexible.
//==============================================================================
template<TpKernel tker,bool simulate2d> void JSphCpu::CalcFlexStrucDefGrad(unsigned np,const tdouble3* pos,const typecode* code
    ,const StFlexStrucData* flexstrucdata,const unsigned* flexstrucridp
    ,const tdouble3* pos0,const unsigned* numpairs,const unsigned* const* pairidx,const tmatrix3f* kercorr
    ,tmatrix3f* defgrad)const
{
  //-Starts execution using OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const int pfs1=p;      //-Number of particle.

    //-Get number of pairs for this particle.
    const unsigned numpairsp1=numpairs[pfs1];

    //-If this particle has pairs.
    if(numpairsp1){
      tmatrix3f defgradp1={0};

      //-Obtains basic data of particle p1.
      const unsigned p1=flexstrucridp[pfs1];
      const float vol0p1=flexstrucdata[CODE_GetIbodyFlexStruc(code[p1])].vol0;
      const tdouble3 posp1=pos[p1];
      const tdouble3 pos0p1=pos0[pfs1];
      const tmatrix3f kercorrp1=kercorr[pfs1];

      //-Calculate deformation gradient.
      for(unsigned pair=0;pair<numpairsp1;pair++){
        const unsigned pfs2=pairidx[pfs1][pair];
        const unsigned p2=flexstrucridp[pfs2];
        const tdouble3 posp2=pos[p2];
        const tdouble3 pos0p2=pos0[pfs2];
        const float drx=float(posp1.x-posp2.x);
        const float dry=float(posp1.y-posp2.y);
        const float drz=float(posp1.z-posp2.z);
        const float drx0=float(pos0p1.x-pos0p2.x);
        const float dry0=float(pos0p1.y-pos0p2.y);
        const float drz0=float(pos0p1.z-pos0p2.z);
        const float rr20=drx0*drx0+dry0*dry0+drz0*drz0;
        const float fac0=fsph::GetKernel_Fac<tker>(CSP,rr20);
        const float frx0=fac0*drx0,fry0=fac0*dry0,frz0=fac0*drz0; //-Gradients.
        defgradp1.a11-=vol0p1*drx*frx0; defgradp1.a12-=vol0p1*drx*fry0; defgradp1.a13-=vol0p1*drx*frz0;
        defgradp1.a21-=vol0p1*dry*frx0; defgradp1.a22-=vol0p1*dry*fry0; defgradp1.a23-=vol0p1*dry*frz0;
        defgradp1.a31-=vol0p1*drz*frx0; defgradp1.a32-=vol0p1*drz*fry0; defgradp1.a33-=vol0p1*drz*frz0;
      }
      if(simulate2d){
        defgradp1.a12=defgradp1.a21=defgradp1.a23=defgradp1.a32=0.0;
        defgradp1.a22=1.0;
      }
      defgrad[pfs1]=fmath::MulMatrix3x3(defgradp1,kercorrp1);
    }
  }
}

//==============================================================================
/// Calculates the deformation gradient matrix for each flexible structure particle.
/// Calcula la matriz de gradiente de deformación para cada partícula de estructura flexible.
//==============================================================================
template<TpKernel tker,bool simulate2d> void JSphCpu::CalcFlexStrucDefGradT()const{
  if(CaseNflexstruc){
    CalcFlexStrucDefGrad<tker,simulate2d>
      (CaseNflexstruc,Pos_c->cptr(),Code_c->cptr(),FlexStrucDatac,FlexStrucRidpc,Pos0c,NumPairsc,PairIdxc,KerCorrc,DefGradc);
  }
}

//==============================================================================
/// Calculates the deformation gradient matrix for each flexible structure particle.
/// Calcula la matriz de gradiente de deformación para cada partícula de estructura flexible.
//==============================================================================
template<TpKernel tker> void JSphCpu::CalcFlexStrucDefGrad_ct0()const{
  if(Simulate2D)CalcFlexStrucDefGradT<tker,true> ();
  else          CalcFlexStrucDefGradT<tker,false>();
}

//==============================================================================
/// Calculates the surface normal for each flexible structure particle.
/// Calcula la superficie normal para cada partícula de estructura flexible.
//==============================================================================
void JSphCpu::CalcFlexStrucNormals()const
{
  const int nflex=int(CaseNflexstruc);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(nflex>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<nflex;p++){
    const tfloat3 boundnor0=BoundNor0c[p];
    if(boundnor0.x!=0 || boundnor0.y!=0 || boundnor0.z!=0){
      const tmatrix3f rot=fmath::InverseMatrix3x3(fmath::TrasMatrix3x3(DefGradc[p]));
      tfloat3 boundnornew;
      boundnornew.x=rot.a11*boundnor0.x+rot.a12*boundnor0.y+rot.a13*boundnor0.z;
      boundnornew.y=rot.a21*boundnor0.x+rot.a22*boundnor0.y+rot.a23*boundnor0.z;
      boundnornew.z=rot.a31*boundnor0.x+rot.a32*boundnor0.y+rot.a33*boundnor0.z;
      const float mag0=sqrt(boundnor0.x*boundnor0.x+boundnor0.y*boundnor0.y+boundnor0.z*boundnor0.z);
      const float magnew=sqrt(boundnornew.x*boundnornew.x+boundnornew.y*boundnornew.y+boundnornew.z*boundnornew.z);
      boundnornew.x*=mag0/magnew; boundnornew.y*=mag0/magnew; boundnornew.z*=mag0/magnew;
      BoundNor_c->ptr()[FlexStrucRidpc[p]]=boundnornew;
    }
  }
}

//==============================================================================
/// Updates the geometric information for each flexible structure particle.
/// Actualiza la información geométrica de cada partícula de estructura flexible.
//==============================================================================
void JSphCpu::UpdateFlexStrucGeometry()const{
  if(TKernel==KERNEL_Wendland)  CalcFlexStrucDefGrad_ct0<KERNEL_Wendland>();
  else if(TKernel==KERNEL_Cubic)CalcFlexStrucDefGrad_ct0<KERNEL_Cubic>   ();
  if(CaseNflexstruc&&UseNormals)CalcFlexStrucNormals();
}

//==============================================================================
/// Calculates the PK1 stress matrix for each flexible structure particle.
/// Calcula la matriz de tensión PK1 para cada partícula de estructura flexible.
//==============================================================================
inline tmatrix3f JSphCpu::CalcFlexStrucPK1Stress(const tmatrix3f& defgrad,const tmatrix6f& cmat)const{
  //-Calculate Green-Lagrange strain from deformation gradient.
  tmatrix3f gl=fmath::MulMatrix3x3(fmath::TrasMatrix3x3(defgrad),defgrad);
  gl.a11-=1.0f; gl.a22-=1.0f; gl.a33-=1.0f;
  gl.a11*=0.5f; gl.a12*=0.5f; gl.a13*=0.5f;
  gl.a21*=0.5f; gl.a22*=0.5f; gl.a23*=0.5f;
  gl.a31*=0.5f; gl.a32*=0.5f; gl.a33*=0.5f;
  //-Convert to Voigt notation.
  const float gl_1=gl.a11;
  const float gl_2=gl.a22;
  const float gl_3=gl.a33;
  const float gl_4=gl.a23+gl.a32;
  const float gl_5=gl.a31+gl.a13;
  const float gl_6=gl.a12+gl.a21;
  //-Multiply with stiffness tensor to get PK2 stress in Voigt form.
  const float pk2_1=cmat.a11*gl_1+cmat.a12*gl_2+cmat.a13*gl_3+cmat.a14*gl_4+cmat.a15*gl_5+cmat.a16*gl_6;
  const float pk2_2=cmat.a21*gl_1+cmat.a22*gl_2+cmat.a23*gl_3+cmat.a24*gl_4+cmat.a25*gl_5+cmat.a26*gl_6;
  const float pk2_3=cmat.a31*gl_1+cmat.a32*gl_2+cmat.a33*gl_3+cmat.a34*gl_4+cmat.a35*gl_5+cmat.a36*gl_6;
  const float pk2_4=cmat.a41*gl_1+cmat.a42*gl_2+cmat.a43*gl_3+cmat.a44*gl_4+cmat.a45*gl_5+cmat.a46*gl_6;
  const float pk2_5=cmat.a51*gl_1+cmat.a52*gl_2+cmat.a53*gl_3+cmat.a54*gl_4+cmat.a55*gl_5+cmat.a56*gl_6;
  const float pk2_6=cmat.a61*gl_1+cmat.a62*gl_2+cmat.a63*gl_3+cmat.a64*gl_4+cmat.a65*gl_5+cmat.a66*gl_6;
  //-Convert PK2 stress back to normal notation.
  const tmatrix3f pk2={pk2_1,pk2_6,pk2_5,pk2_6,pk2_2,pk2_4,pk2_5,pk2_4,pk2_3};
  //-Convert PK2 to PK1 and return.
  return fmath::MulMatrix3x3(defgrad,pk2);
}

//==============================================================================
/// Interaction forces for the flexible structure particles.
/// Fuerzas de interacción para las partículas de estructura flexible.
//==============================================================================
template<TpKernel tker,TpVisco tvisco> void JSphCpu::InteractionForcesFlexStruc(unsigned np,float visco
    ,StDivDataCpu divdata,const unsigned* dcell
    ,const tdouble3* pos,const tfloat4* velrhop,const float* press,const typecode* code
    ,const StFlexStrucData* flexstrucdata,const unsigned* flexstrucridp
    ,const tdouble3* pos0,const unsigned* numpairs,const unsigned* const* pairidx,const tmatrix3f* kercorr,const tmatrix3f* defgrad
    ,float& flexstrucdt,tfloat3* ace)const{
  //-Initialise flexstructh to calculate max flexstrucdt.
  float flexstructh[OMP_MAXTHREADS*OMP_STRIDE];
  for(int th=0;th<OmpThreads;th++)flexstructh[th*OMP_STRIDE]=0;
  //-Starts execution using OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const int pfs1=p;      //-Number of particle.
    const unsigned p1=flexstrucridp[pfs1];

    //-Get codep1.
    const typecode codep1=code[p1];
    if(CODE_IsFlexStrucFlex(codep1)){
      tfloat3 acep1=TFloat3(0);

      //-Obtains basic data of particle p1.
      const tdouble3 posp1=pos[p1];
      const tfloat4 velrhop1=velrhop[p1];
      const float pressp1=press[p1];
      const tdouble3 pos0p1=pos0[pfs1];
      const tmatrix3f kercorrp1=kercorr[pfs1];
      const tmatrix3f defgradp1=defgrad[pfs1];

      //-Obtains flexible structure data.
      const float vol0p1=flexstrucdata[CODE_GetIbodyFlexStruc(codep1)].vol0;
      const float rho0p1=flexstrucdata[CODE_GetIbodyFlexStruc(codep1)].rho0;
      const float youngmod=flexstrucdata[CODE_GetIbodyFlexStruc(codep1)].youngmod;
      const float poissratio=flexstrucdata[CODE_GetIbodyFlexStruc(codep1)].poissratio;
      const float hgfactor=flexstrucdata[CODE_GetIbodyFlexStruc(codep1)].hgfactor;
      const tmatrix6f cmat=flexstrucdata[CODE_GetIbodyFlexStruc(codep1)].cmat;
      const tmatrix3f pk1p1=CalcFlexStrucPK1Stress(defgradp1,cmat);
      const tmatrix3f pk1kercorrp1=fmath::MulMatrix3x3(pk1p1,kercorrp1);

      //-Get current mass of flexible structure particle.
      const float mass0p1=vol0p1*rho0p1;
      const float rhop1=rho0p1/fmath::Determinant3x3(defgradp1);

      //-Calculate structural speed of sound.
      const float csp1=float( sqrt(youngmod*(1.0-poissratio)/(rhop1*(1.0+poissratio)*(1.0-2.0*poissratio))) );

      //-Obtains neighborhood search limits.
      const StNgSearch ngs=nsearch::Init(dcell[p1],false,divdata);

      //-Flexible structure-Fluid interaction.
      for(int z=ngs.zini;z<ngs.zfin;z++){
        for(int y=ngs.yini;y<ngs.yfin;y++){
          const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);
          for(unsigned p2=pif.x;p2<pif.y;p2++){
            if(CODE_IsFluid(code[p2])){
              const tdouble3 posp2=pos[p2];
              const float drx=float(posp1.x-posp2.x);
              const float dry=float(posp1.y-posp2.y);
              const float drz=float(posp1.z-posp2.z);
              const float rr2=drx*drx+dry*dry+drz*drz;
              if(rr2<=KernelSize2&&rr2>=ALMOSTZERO){
                //-Computes kernel.
                const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
                const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
                tfloat4 velrhop2=velrhop[p2];
                const float dvx=velrhop1.x-velrhop2.x,dvy=velrhop1.y-velrhop2.y,dvz=velrhop1.z-velrhop2.z;
                //-Artificial viscosity.
                if(tvisco==VISCO_Artificial){
                  const float dot=drx*dvx+dry*dvy+drz*dvz;
                  if(dot<0){
                    const float dot_rr2=dot/(rr2+Eta2);
                    const float amubar=KernelH*dot_rr2;
                    const float robar=(velrhop1.w+velrhop2.w)*0.5f;
                    const float pi_visc=(float)(-visco*Cs0*amubar/robar)*MassFluid*(MassFluid/mass0p1);
                    acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
                  }
                }
                  //-Laminar viscosity.
                else{
                  const float robar2=(velrhop1.w+velrhop2.w);
                  const float temp=4.f*visco/((rr2+Eta2)*robar2);
                  const float vtemp=MassFluid*temp*(drx*frx+dry*fry+drz*frz)*(MassFluid/mass0p1);
                  acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
                }
                //-Velocity derivative (Momentum equation).
                const float pressp2=press[p2];
                const float prs=(pressp1+pressp2)/(velrhop1.w*velrhop2.w)+(tker==KERNEL_Cubic?fsph::GetKernelCubic_Tensil(CSP,rr2,velrhop1.w,pressp1,velrhop2.w,pressp2):0);
                const float p_vpm=-prs*MassFluid*(MassFluid/mass0p1);
                acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
              }
            }
          }
        }
      }

      //-Loop through pairs and calculate forces.
      for(unsigned pair=0;pair<numpairs[pfs1];pair++){
        const unsigned pfs2=pairidx[pfs1][pair];
        const tdouble3 pos0p2=pos0[pfs2];
        const float drx0=float(pos0p1.x-pos0p2.x);
        const float dry0=float(pos0p1.y-pos0p2.y);
        const float drz0=float(pos0p1.z-pos0p2.z);
        const float rr20=drx0*drx0+dry0*dry0+drz0*drz0;
        const float fac0=fsph::GetKernel_Fac<tker>(CSP,rr20);
        const float frx0=fac0*drx0,fry0=fac0*dry0,frz0=fac0*drz0; //-Gradients.
        //-Acceleration due to structure.
        const tmatrix3f defgradp2=defgrad[pfs2];
        const tmatrix3f pk1p2=CalcFlexStrucPK1Stress(defgradp2,cmat);
        const tmatrix3f kercorrp2=kercorr[pfs2];
        const tmatrix3f pk1kercorrp2=fmath::MulMatrix3x3(pk1p2,kercorrp2);
        tmatrix3f pk1kercorrp1p2;
        pk1kercorrp1p2.a11=pk1kercorrp1.a11+pk1kercorrp2.a11; pk1kercorrp1p2.a12=pk1kercorrp1.a12+pk1kercorrp2.a12; pk1kercorrp1p2.a13=pk1kercorrp1.a13+pk1kercorrp2.a13;
        pk1kercorrp1p2.a21=pk1kercorrp1.a21+pk1kercorrp2.a21; pk1kercorrp1p2.a22=pk1kercorrp1.a22+pk1kercorrp2.a22; pk1kercorrp1p2.a23=pk1kercorrp1.a23+pk1kercorrp2.a23;
        pk1kercorrp1p2.a31=pk1kercorrp1.a31+pk1kercorrp2.a31; pk1kercorrp1p2.a32=pk1kercorrp1.a32+pk1kercorrp2.a32; pk1kercorrp1p2.a33=pk1kercorrp1.a33+pk1kercorrp2.a33;
        tfloat3 pk1kercorrdw;
        pk1kercorrdw.x=pk1kercorrp1p2.a11*frx0+pk1kercorrp1p2.a12*fry0+pk1kercorrp1p2.a13*frz0;
        pk1kercorrdw.y=pk1kercorrp1p2.a21*frx0+pk1kercorrp1p2.a22*fry0+pk1kercorrp1p2.a23*frz0;
        pk1kercorrdw.z=pk1kercorrp1p2.a31*frx0+pk1kercorrp1p2.a32*fry0+pk1kercorrp1p2.a33*frz0;
        acep1.x+=pk1kercorrdw.x*vol0p1/rho0p1; acep1.y+=pk1kercorrdw.y*vol0p1/rho0p1; acep1.z+=pk1kercorrdw.z*vol0p1/rho0p1;
        //-Hourglass correction.
        if(hgfactor){
          const unsigned p2=flexstrucridp[pfs2];
          const tdouble3 posp2=pos[p2];
          const float drx=float(posp1.x-posp2.x);
          const float dry=float(posp1.y-posp2.y);
          const float drz=float(posp1.z-posp2.z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          const float wab=fsph::GetKernel_Wab<tker>(CSP,rr2);
          const float xijbarx=defgradp1.a11*drx0+defgradp1.a12*dry0+defgradp1.a13*drz0;
          const float xijbary=defgradp1.a21*drx0+defgradp1.a22*dry0+defgradp1.a23*drz0;
          const float xijbarz=defgradp1.a31*drx0+defgradp1.a32*dry0+defgradp1.a33*drz0;
          const float xjibarx=-(defgradp2.a11*drx0+defgradp2.a12*dry0+defgradp2.a13*drz0);
          const float xjibary=-(defgradp2.a21*drx0+defgradp2.a22*dry0+defgradp2.a23*drz0);
          const float xjibarz=-(defgradp2.a31*drx0+defgradp2.a32*dry0+defgradp2.a33*drz0);
          const tfloat3 epsij=TFloat3(drx-xijbarx,dry-xijbary,drz-xijbarz);
          const tfloat3 epsji=TFloat3(-drx-xjibarx,-dry-xjibary,-drz-xjibarz);
          const float rr=sqrt(rr2);
          const float deltaij=(epsij.x*drx+epsij.y*dry+epsij.z*drz)/rr;
          const float deltaji=-(epsji.x*drx+epsji.y*dry+epsji.z*drz)/rr;
          const float mulfac=(hgfactor*vol0p1*vol0p1*wab*youngmod/(rr20*rr*mass0p1)*0.5f*(deltaij+deltaji));
          acep1.x-=mulfac*drx; acep1.y-=mulfac*dry; acep1.z-=mulfac*drz;
        }
      }

      //-Store results.
      if(acep1.x||acep1.y||acep1.z||csp1){
        ace[p1]=ace[p1]+acep1;
        const int th=omp_get_thread_num();
        if(csp1>flexstructh[th*OMP_STRIDE])flexstructh[th*OMP_STRIDE]=csp1;
      }
    }
  }
  //-Keep max value in flexstrucdt.
  for(int th=0;th<OmpThreads;th++)if(flexstrucdt<flexstructh[th*OMP_STRIDE])flexstrucdt=flexstructh[th*OMP_STRIDE];
}

//==============================================================================
/// Interaction forces for the flexible structure particles.
/// Fuerzas de interacción para las partículas de estructura flexible.
//==============================================================================
template<TpKernel tker,TpVisco tvisco> void JSphCpu::Interaction_ForcesFlexStrucT(float& flexstrucdtmax)const{
  if(CaseNflexstruc){
    InteractionForcesFlexStruc<tker,tvisco>
      (CaseNflexstruc,Visco*ViscoBoundFactor,DivData,Dcell_c->cptr(),Pos_c->cptr(),Velrho_c->cptr(),Press_c->cptr(),Code_c->cptr(),FlexStrucDatac,FlexStrucRidpc,Pos0c,NumPairsc,PairIdxc,KerCorrc,DefGradc,flexstrucdtmax,Ace_c->ptr());
  }
}

//==============================================================================
/// Interaction forces for the flexible structure particles.
/// Fuerzas de interacción para las partículas de estructura flexible.
//==============================================================================
template<TpKernel tker> void JSphCpu::Interaction_ForcesFlexStruc_ct0(float& flexstrucdtmax)const{
       if(TVisco==VISCO_Artificial)Interaction_ForcesFlexStrucT<tker,VISCO_Artificial>(flexstrucdtmax);
  else if(TVisco==VISCO_Laminar   )Interaction_ForcesFlexStrucT<tker,VISCO_Laminar   >(flexstrucdtmax);
  else if(TVisco==VISCO_LaminarSPS)Interaction_ForcesFlexStrucT<tker,VISCO_LaminarSPS>(flexstrucdtmax);
}

//==============================================================================
/// Interaction forces for the flexible structure particles.
/// Fuerzas de interacción para las partículas de estructura flexible.
//==============================================================================
void JSphCpu::Interaction_ForcesFlexStruc(float& flexstrucdtmax)const{
  if(TKernel==KERNEL_Wendland)  Interaction_ForcesFlexStruc_ct0<KERNEL_Wendland>(flexstrucdtmax);
  else if(TKernel==KERNEL_Cubic)Interaction_ForcesFlexStruc_ct0<KERNEL_Cubic>   (flexstrucdtmax);
}

//==============================================================================
/// Updates position and velocity using semi-implicit Euler scheme (used with Verlet scheme).
/// Actualiza la posición y la velocidad usando el esquema de Euler semiimplícito (usado con el esquema de Verlet).
//==============================================================================
void JSphCpu::ComputeSemiImplicitEulerFlexStruc(double dt,tdouble3* pos,unsigned* dcell,typecode* code)const
{
  const tdouble3 gravity=ToTDouble3(Gravity);
  const int nflex=int(CaseNflexstruc);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(nflex>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<nflex;p++){
    const unsigned p1=FlexStrucRidpc[p]; //-Number of particle.
    if(CODE_IsFlexStrucFlex(code[p1])){
      const tfloat4 rvelrhop=Velrho_c->cptr()[p1];
      const tfloat3 race=Ace_c->cptr()[p1];
      const double acegrx=double(race.x)+gravity.x;
      const double acegry=double(race.y)+gravity.y;
      const double acegrz=double(race.z)+gravity.z;
      const tfloat4 rvelrhopnew=TFloat4(
          float(double(rvelrhop.x) + acegrx*dt),
          float(double(rvelrhop.y) + acegry*dt),
          float(double(rvelrhop.z) + acegrz*dt),
          VelrhoM1_c->cptr()[p1].w);
      //-Calculate displacement. | Calcula desplazamiento.
      const double dx=double(rvelrhopnew.x)*dt;
      const double dy=double(rvelrhopnew.y)*dt;
      const double dz=double(rvelrhopnew.z)*dt;
      bool outrhop=(rvelrhopnew.w<RhopOutMin || rvelrhopnew.w>RhopOutMax);
      //-Update particle data.
      UpdatePos(pos[p1],dx,dy,dz,outrhop,p1,pos,dcell,code);
      VelrhoM1_c->ptr()[p1]=rvelrhopnew;
    }
  }
}

//==============================================================================
/// Updates position and velocity using symplectic predictor scheme.
/// Actualiza la posición y la velocidad utilizando un esquema predictor simpléctico.
//==============================================================================
void JSphCpu::ComputeSymplecticPreFlexStruc(double dtm,tdouble3* pos,unsigned* dcell,typecode* code)const
{
  const tdouble3 gravity=ToTDouble3(Gravity);
  const int nflex=int(CaseNflexstruc);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(nflex>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<nflex;p++){
    const unsigned p1=FlexStrucRidpc[p]; //-Number of particle.
    if(CODE_IsFlexStrucFlex(code[p1])){
      const tfloat4 rvelrhoppre=VelrhoPre_c->cptr()[p1];
      const tfloat3 race=Ace_c->cptr()[p1];
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=double(rvelrhoppre.x)*dtm;
      double dy=double(rvelrhoppre.y)*dtm;
      double dz=double(rvelrhoppre.z)*dtm;
      //-Calculate velocity & density. | Calcula velocidad y densidad.
      const tfloat4 rvelrhopnew=TFloat4(
          float(double(rvelrhoppre.x) + (double(race.x)+gravity.x) * dtm),
          float(double(rvelrhoppre.y) + (double(race.y)+gravity.y) * dtm),
          float(double(rvelrhoppre.z) + (double(race.z)+gravity.z) * dtm),
          Velrho_c->cptr()[p1].w);
      bool outrhop=(rvelrhopnew.w<RhopOutMin || rvelrhopnew.w>RhopOutMax);
      if(outrhop && CODE_IsNormal(code[p1]))code[p1]=CODE_SetOutRho(code[p1]); //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
      //-Update particle data.
      Velrho_c->ptr()[p1]=rvelrhopnew;
      outrhop=CODE_IsOutRho(code[p1]);
      const bool flexstruc=CODE_IsFlexStrucFlex(code[p1]);
      const bool normal=(outrhop || CODE_IsNormal(code[p1]));
      if(normal && flexstruc)UpdatePos(PosPre_c->cptr()[p1],dx,dy,dz,outrhop,p1,pos,dcell,code);
    }
  }
}

//==============================================================================
/// Updates position and velocity using symplectic corrector scheme.
/// Actualiza la posición y la velocidad utilizando un esquema corrector simpléctico.
//==============================================================================
void JSphCpu::ComputeSymplecticCorrFlexStruc(double dtm,double dt,tdouble3* pos,unsigned* dcell,typecode* code)const
{
  const tdouble3 gravity=ToTDouble3(Gravity);
  const int nflex=int(CaseNflexstruc);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(nflex>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<nflex;p++){
    const unsigned p1=FlexStrucRidpc[p]; //-Number of particle.
    if(CODE_IsFlexStrucFlex(code[p1])){
      const tfloat4 rvelrhoppre=VelrhoPre_c->cptr()[p1];
      const tfloat3 race=Ace_c->cptr()[p1];
      //-Calculate velocity. | Calcula velocidad.
      const tfloat4 rvelrhopnew=TFloat4(
          float(double(rvelrhoppre.x) + (double(race.x)+gravity.x) * dt),
          float(double(rvelrhoppre.y) + (double(race.y)+gravity.y) * dt),
          float(double(rvelrhoppre.z) + (double(race.z)+gravity.z) * dt),
          Velrho_c->cptr()[p1].w);
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=(double(rvelrhoppre.x)+double(rvelrhopnew.x)) * dtm;
      double dy=(double(rvelrhoppre.y)+double(rvelrhopnew.y)) * dtm;
      double dz=(double(rvelrhoppre.z)+double(rvelrhopnew.z)) * dtm;
      bool outrhop=(rvelrhopnew.w<RhopOutMin || rvelrhopnew.w>RhopOutMax);
      if(outrhop && CODE_IsNormal(code[p1]))code[p1]=CODE_SetOutRho(code[p1]); //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
      //-Update particle data.
      Velrho_c->ptr()[p1]=rvelrhopnew;
      outrhop=CODE_IsOutRho(code[p1]);
      const bool flexstruc=CODE_IsFlexStrucFlex(code[p1]);
      const bool normal=(outrhop || CODE_IsNormal(code[p1]));
      if(normal && flexstruc)UpdatePos(PosPre_c->cptr()[p1],dx,dy,dz,outrhop,p1,pos,dcell,code);
    }
  }
}
//<vs_flexstruc_end>
