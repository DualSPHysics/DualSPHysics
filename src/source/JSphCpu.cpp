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

/// \file JSphCpu.cpp \brief Implements the class \ref JSphCpu.

#include "JSphCpu.h"
#include "JCellSearch_inline.h"
#include "JCellDivCpu.h"
#include "JPartFloatBi4.h"
#include "FunSphKernel.h"
#include "FunSphEos.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JDsMotion.h"
#include "JArraysCpu.h"
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
#include "JSphBoundCorr.h"
#include "JSphInOut.h"
#include "JSphShifting.h"

#include <climits>
#ifndef WIN32
#include <unistd.h>
#endif

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JSphCpu::JSphCpu(bool withmpi):JSph(true,false,withmpi){
  ClassName="JSphCpu";
  CellDiv=NULL;
  ArraysCpu=new JArraysCpu;
  InitVars();
  TmcCreation(Timers,false);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphCpu::~JSphCpu(){
  DestructorActive=true;
  FreeCpuMemoryParticles();
  FreeCpuMemoryFixed();
  delete ArraysCpu;
  TmcDestruction(Timers);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphCpu::InitVars(){
  RunMode="";
  OmpThreads=1;

  DivData=DivDataCpuNull();

  Np=Npb=NpbOk=0;
  NpbPer=NpfPer=0;

  Idpc=NULL; Codec=NULL; Dcellc=NULL; Posc=NULL; Velrhopc=NULL;
  BoundNormalc=NULL; MotionVelc=NULL; //-mDBC
  VelrhopM1c=NULL;                //-Verlet
  PosPrec=NULL; VelrhopPrec=NULL; //-Symplectic
  SpsTauc=NULL; SpsGradvelc=NULL; //-Laminar+SPS.
  Arc=NULL; Acec=NULL; Deltac=NULL;
  ShiftPosfsc=NULL;               //-Shifting.
  Pressc=NULL;
  RidpMove=NULL; 
  FtRidp=NULL;
  FtoForces=NULL;
  FtoForcesRes=NULL;
  FreeCpuMemoryParticles();
  FreeCpuMemoryFixed();
}

//==============================================================================
/// Deallocate fixed memory on CPU for moving and floating bodies.
/// Libera memoria fija en cpu para moving y floating.
//==============================================================================
void JSphCpu::FreeCpuMemoryFixed(){
  MemCpuFixed=0;
  delete[] RidpMove;     RidpMove=NULL;
  delete[] FtRidp;       FtRidp=NULL;
  delete[] FtoForces;    FtoForces=NULL;
  delete[] FtoForcesRes; FtoForcesRes=NULL;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphCpu::AllocCpuMemoryFixed(){
  MemCpuFixed=0;
  try{
    //-Allocates memory for moving objects.
    if(CaseNmoving){
      RidpMove=new unsigned[CaseNmoving];  MemCpuFixed+=(sizeof(unsigned)*CaseNmoving);
    }
    //-Allocates memory for floating bodies.
    if(CaseNfloat){
      FtRidp      =new unsigned[CaseNfloat];     MemCpuFixed+=(sizeof(unsigned)*CaseNfloat);
      FtoForces   =new StFtoForces[FtCount];     MemCpuFixed+=(sizeof(StFtoForces)*FtCount);
      FtoForcesRes=new StFtoForcesRes[FtCount];  MemCpuFixed+=(sizeof(StFtoForcesRes)*FtCount);
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Deallocate memory in CPU for particles.
/// Libera memoria en cpu para particulas.
//==============================================================================
void JSphCpu::FreeCpuMemoryParticles(){
  CpuParticlesSize=0;
  MemCpuParticles=0;
  ArraysCpu->Reset();
}

//==============================================================================
/// Allocte memory on CPU for the particles. 
/// Reserva memoria en Cpu para las particulas. 
//==============================================================================
void JSphCpu::AllocCpuMemoryParticles(unsigned np,float over){
  FreeCpuMemoryParticles();
  //-Calculate number of partices with reserved memory | Calcula numero de particulas para las que se reserva memoria.
  const unsigned np2=(over>0? unsigned(over*np): np);
  CpuParticlesSize=np2+PARTICLES_OVERMEMORY_MIN;
  //-Define number or arrays to use. | Establece numero de arrays a usar.
  ArraysCpu->SetArraySize(CpuParticlesSize);
  #ifdef CODE_SIZE4
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_4B,2);  //-code,code2
  #else
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_2B,2);  //-code,code2
  #endif
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_4B,5);  //-idp,ar,viscdt,dcell,prrhop
  if(DDTArray)ArraysCpu->AddArrayCount(JArraysCpu::SIZE_4B,1);  //-delta
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_12B,1); //-ace
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,2); //-velrhop,poscell
  ArraysCpu->AddArrayCount(JArraysCpu::SIZE_24B,2); //-pos
  if(TStep==STEP_Verlet){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,1); //-velrhopm1
  }
  else if(TStep==STEP_Symplectic){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_24B,1); //-pospre
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,1); //-velrhoppre
  }
  if(TVisco==VISCO_LaminarSPS){     
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_24B,2); //-SpsTau,SpsGradvel
  }
  if(Shifting){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,1); //-shiftposfs
  }
  if(UseNormals){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_12B,1); //-BoundNormal
    if(SlipMode!=SLIP_Vel0)ArraysCpu->AddArrayCount(JArraysCpu::SIZE_12B,1); //-MotionVel
  }
  if(InOut){
    //ArraysCpu->AddArrayCount(JArraysCpu::SIZE_4B,1);  //-InOutPart
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_1B,2);  //-newizone,zsurfok
  }
  //-Shows the allocated memory.
  MemCpuParticles=ArraysCpu->GetAllocMemoryCpu();
  PrintSizeNp(CpuParticlesSize,MemCpuParticles,0);
}

//==============================================================================
/// Resizes space in CPU memory for particles.
//==============================================================================
void JSphCpu::ResizeCpuMemoryParticles(unsigned npnew){
  npnew=npnew+PARTICLES_OVERMEMORY_MIN;
  //-Saves current data from CPU.
  unsigned    *idp        =SaveArrayCpu(Np,Idpc);
  typecode    *code       =SaveArrayCpu(Np,Codec);
  unsigned    *dcell      =SaveArrayCpu(Np,Dcellc);
  tdouble3    *pos        =SaveArrayCpu(Np,Posc);
  tfloat4     *velrhop    =SaveArrayCpu(Np,Velrhopc);
  tfloat4     *velrhopm1  =SaveArrayCpu(Np,VelrhopM1c);
  tdouble3    *pospre     =SaveArrayCpu(Np,PosPrec);
  tfloat4     *velrhoppre =SaveArrayCpu(Np,VelrhopPrec);
  tsymatrix3f *spstau     =SaveArrayCpu(Np,SpsTauc);
  tfloat3     *boundnormal=SaveArrayCpu(Np,BoundNormalc);
  tfloat3     *motionvel  =SaveArrayCpu(Np,MotionVelc);
  //-Frees pointers.
  ArraysCpu->Free(Idpc);
  ArraysCpu->Free(Codec);
  ArraysCpu->Free(Dcellc);
  ArraysCpu->Free(Posc);
  ArraysCpu->Free(Velrhopc);
  ArraysCpu->Free(VelrhopM1c);
  ArraysCpu->Free(PosPrec);
  ArraysCpu->Free(VelrhopPrec);
  ArraysCpu->Free(SpsTauc);
  ArraysCpu->Free(BoundNormalc);
  ArraysCpu->Free(MotionVelc);
  //-Resizes CPU memory allocation.
  const double mbparticle=(double(MemCpuParticles)/(1024*1024))/CpuParticlesSize; //-MB por particula.
  Log->Printf("**JSphCpu: Requesting cpu memory for %u particles: %.1f MB.",npnew,mbparticle*npnew);
  ArraysCpu->SetArraySize(npnew);
  //-Reserve pointers.
  Idpc    =ArraysCpu->ReserveUint();
  Codec   =ArraysCpu->ReserveTypeCode();
  Dcellc  =ArraysCpu->ReserveUint();
  Posc    =ArraysCpu->ReserveDouble3();
  Velrhopc=ArraysCpu->ReserveFloat4();
  if(velrhopm1)  VelrhopM1c  =ArraysCpu->ReserveFloat4();
  if(pospre)     PosPrec     =ArraysCpu->ReserveDouble3();
  if(velrhoppre) VelrhopPrec =ArraysCpu->ReserveFloat4();
  if(spstau)     SpsTauc     =ArraysCpu->ReserveSymatrix3f();
  if(boundnormal)BoundNormalc=ArraysCpu->ReserveFloat3();
  if(motionvel)  MotionVelc  =ArraysCpu->ReserveFloat3();
  //-Restore data in CPU memory.
  RestoreArrayCpu(Np,idp,Idpc);
  RestoreArrayCpu(Np,code,Codec);
  RestoreArrayCpu(Np,dcell,Dcellc);
  RestoreArrayCpu(Np,pos,Posc);
  RestoreArrayCpu(Np,velrhop,Velrhopc);
  RestoreArrayCpu(Np,velrhopm1,VelrhopM1c);
  RestoreArrayCpu(Np,pospre,PosPrec);
  RestoreArrayCpu(Np,velrhoppre,VelrhopPrec);
  RestoreArrayCpu(Np,spstau,SpsTauc);
  RestoreArrayCpu(Np,boundnormal,BoundNormalc);
  RestoreArrayCpu(Np,motionvel,MotionVelc);
  //-Updates values.
  CpuParticlesSize=npnew;
  MemCpuParticles=ArraysCpu->GetAllocMemoryCpu();
}

//==============================================================================
/// Saves a CPU array in CPU memory. 
//==============================================================================
template<class T> T* JSphCpu::TSaveArrayCpu(unsigned np,const T *datasrc)const{
  T *data=NULL;
  if(datasrc){
    try{
      data=new T[np];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
    memcpy(data,datasrc,sizeof(T)*np);
  }
  return(data);
}

//==============================================================================
/// Restores an array (generic) from CPU memory. 
//==============================================================================
template<class T> void JSphCpu::TRestoreArrayCpu(unsigned np,T *data,T *datanew)const{
  if(data&&datanew)memcpy(datanew,data,sizeof(T)*np);
  delete[] data;
}

//==============================================================================
/// Arrays for basic particle data. 
/// Arrays para datos basicos de las particulas. 
//==============================================================================
void JSphCpu::ReserveBasicArraysCpu(){
  Idpc=ArraysCpu->ReserveUint();
  Codec=ArraysCpu->ReserveTypeCode();
  Dcellc=ArraysCpu->ReserveUint();
  Posc=ArraysCpu->ReserveDouble3();
  Velrhopc=ArraysCpu->ReserveFloat4();
  if(TStep==STEP_Verlet)VelrhopM1c=ArraysCpu->ReserveFloat4();
  if(TVisco==VISCO_LaminarSPS)SpsTauc=ArraysCpu->ReserveSymatrix3f();
  if(UseNormals){
    BoundNormalc=ArraysCpu->ReserveFloat3();
    if(SlipMode!=SLIP_Vel0)MotionVelc=ArraysCpu->ReserveFloat3();
  }
}

//==============================================================================
/// Return memory reserved on CPU.
/// Devuelve la memoria reservada en cpu.
//==============================================================================
llong JSphCpu::GetAllocMemoryCpu()const{  
  llong s=JSph::GetAllocMemoryCpu();
  //-Reserved in AllocCpuMemoryParticles().
  s+=MemCpuParticles;
  //-Reserved in AllocCpuMemoryFixed().
  s+=MemCpuFixed;
  //-Reserved in other objects.
  if(MLPistons)s+=MLPistons->GetAllocMemoryCpu();
  return(s);
}

//==============================================================================
/// Visualize the reserved memory.
/// Visualiza la memoria reservada.
//==============================================================================
void JSphCpu::PrintAllocMemory(llong mcpu)const{
  Log->Printf("Allocated memory in CPU: %lld (%.2f MB)",mcpu,double(mcpu)/(1024*1024));
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
  ,unsigned *idp,tdouble3 *pos,tfloat3 *vel,float *rhop,typecode *code)
{
  unsigned num=n;
  //-Copy selected values.
  if(code)memcpy(code,Codec+pini,sizeof(typecode)*n);
  if(idp) memcpy(idp ,Idpc +pini,sizeof(unsigned)*n);
  if(pos) memcpy(pos ,Posc +pini,sizeof(tdouble3)*n);
  if(vel && rhop){
    for(unsigned p=0;p<n;p++){
      tfloat4 vr=Velrhopc[p+pini];
      vel[p]=TFloat3(vr.x,vr.y,vr.z);
      rhop[p]=vr.w;
    }
  }
  else{
    if(vel) for(unsigned p=0;p<n;p++){ tfloat4 vr=Velrhopc[p+pini]; vel[p]=TFloat3(vr.x,vr.y,vr.z); }
    if(rhop)for(unsigned p=0;p<n;p++)rhop[p]=Velrhopc[p+pini].w;
  }
  //-Eliminate non-normal particles (periodic & others). | Elimina particulas no normales (periodicas y otras).
  if(onlynormal){
    if(!idp || !pos || !vel || !rhop)Run_Exceptioon("Pointers without data.");
    typecode *code2=code;
    if(!code2){
      code2=ArraysCpu->ReserveTypeCode();
      memcpy(code2,Codec+pini,sizeof(typecode)*n);
    }
    unsigned ndel=0;
    for(unsigned p=0;p<n;p++){
      bool normal=CODE_IsNormal(code2[p]);
      if(ndel && normal){
        const unsigned pdel=p-ndel;
        idp[pdel]  =idp[p];
        pos[pdel]  =pos[p];
        vel[pdel]  =vel[p];
        rhop[pdel] =rhop[p];
        code2[pdel]=code2[p];
      }
      if(!normal)ndel++;
    }
    num-=ndel;
    if(!code)ArraysCpu->Free(code2);
  }
  return(num);
}

//==============================================================================
/// Load the execution configuration with OpenMP.
/// Carga la configuracion de ejecucion con OpenMP.
//==============================================================================
void JSphCpu::ConfigOmp(const JSphCfgRun *cfg){
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
void JSphCpu::ConfigRunMode(const JSphCfgRun *cfg,std::string preinfo){
  #ifndef WIN32
    const int len=128; char hname[len];
    gethostname(hname,len);
    preinfo=preinfo+(!preinfo.empty()? " - ": "")+"HostName:"+hname;
  #endif
  Hardware="Cpu";
  if(OmpThreads==1)RunMode="Single core";
  else RunMode=string("OpenMP(Threads:")+fun::IntStr(OmpThreads)+")";
  if(!preinfo.empty())RunMode=preinfo+" - "+RunMode;
  if(Stable)RunMode=string("Stable - ")+RunMode;
  RunMode=string("Pos-Double - ")+RunMode;
  Log->Print(" ");
  Log->Print(fun::VarStr("RunMode",RunMode));
  Log->Print(" ");
}

//==============================================================================
/// Initialisation of arrays and variables for execution.
/// Inicializa vectores y variables para la ejecucion.
//==============================================================================
void JSphCpu::InitRunCpu(){
  InitRun(Np,Idpc,Posc);

  if(TStep==STEP_Verlet)memcpy(VelrhopM1c,Velrhopc,sizeof(tfloat4)*Np);
  if(TVisco==VISCO_LaminarSPS)memset(SpsTauc,0,sizeof(tsymatrix3f)*Np);
  if(CaseNfloat)InitFloating();
  if(MotionVelc)memset(MotionVelc,0,sizeof(tfloat3)*Np);
}

//==============================================================================
/// Prepare variables for interaction functions.
/// Prepara variables para interaccion.
//==============================================================================
void JSphCpu::PreInteractionVars_Forces(unsigned np,unsigned npb){
  //-Initialise arrays.
  const unsigned npf=np-npb;
  memset(Arc,0,sizeof(float)*np);                                    //Arc[]=0
  if(Deltac)memset(Deltac,0,sizeof(float)*np);                       //Deltac[]=0
  memset(Acec,0,sizeof(tfloat3)*np);                                 //Acec[]=(0,0,0)
  if(SpsGradvelc)memset(SpsGradvelc+npb,0,sizeof(tsymatrix3f)*npf);  //SpsGradvelc[]=(0,0,0,0,0,0).

  //-Select particles for shifting.
  if(ShiftPosfsc)Shifting->InitCpu(npf,npb,Posc,ShiftPosfsc);

  //-Adds variable acceleration from input configuration.
  if(AccInput)AccInput->RunCpu(TimeStep,Gravity,npf,npb,Codec,Posc,Velrhopc,Acec);

  //-Prepare press values for interaction.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=0;p<n;p++){
    Pressc[p]=fsph::ComputePress(Velrhopc[p].w,CSP);
  }
}

//==============================================================================
/// Prepare variables for interaction functions.
/// Prepara variables para interaccion.
//==============================================================================
void JSphCpu::PreInteraction_Forces(){
  TmcStart(Timers,TMC_CfPreForces);
  //-Assign memory.
  Arc=ArraysCpu->ReserveFloat();
  Acec=ArraysCpu->ReserveFloat3();
  if(DDTArray)Deltac=ArraysCpu->ReserveFloat();
  if(Shifting)ShiftPosfsc=ArraysCpu->ReserveFloat4();
  Pressc=ArraysCpu->ReserveFloat();
  if(TVisco==VISCO_LaminarSPS)SpsGradvelc=ArraysCpu->ReserveSymatrix3f();

  //-Initialise arrays.
  PreInteractionVars_Forces(Np,Npb);

  //-Calculate VelMax: Floating object particles are included and do not affect use of periodic condition.
  //-Calcula VelMax: Se incluyen las particulas floatings y no afecta el uso de condiciones periodicas.
  const unsigned pini=(DtAllParticles? 0: Npb);
  VelMax=CalcVelMaxOmp(Np-pini,Velrhopc+pini);
  ViscDtMax=0;
  TmcStop(Timers,TMC_CfPreForces);
}

//==============================================================================
/// Returns maximum velocity from an array tfloat4.
/// Devuelve la velociad maxima de un array tfloat4.
//==============================================================================
float JSphCpu::CalcVelMaxSeq(unsigned np,const tfloat4* velrhop)const{
  float velmax=0;
  for(unsigned p=0;p<np;p++){
    const tfloat4 v=velrhop[p];
    const float v2=v.x*v.x+v.y*v.y+v.z*v.z;
    velmax=max(velmax,v2);
  }
  return(sqrt(velmax));
}

//==============================================================================
/// Returns maximum velocity from an array tfloat4 using OpenMP.
/// Devuelve la velociad maxima de un array tfloat4 usando OpenMP.
//==============================================================================
float JSphCpu::CalcVelMaxOmp(unsigned np,const tfloat4* velrhop)const{
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
          const tfloat4 v=velrhop[c];
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
    else if(np)velmax=CalcVelMaxSeq(np,velrhop);
  #else
    if(np)velmax=CalcVelMaxSeq(np,velrhop);
  #endif
  return(velmax);
}

//==============================================================================
/// Free memory assigned to ArraysCpu.
/// Libera memoria asignada de ArraysCpu.
//==============================================================================
void JSphCpu::PosInteraction_Forces(){
  //-Free memory assigned in PreInteraction_Forces(). | Libera memoria asignada en PreInteraction_Forces().
  ArraysCpu->Free(Arc);          Arc=NULL;
  ArraysCpu->Free(Acec);         Acec=NULL;
  ArraysCpu->Free(Deltac);       Deltac=NULL;
  ArraysCpu->Free(ShiftPosfsc);  ShiftPosfsc=NULL;
  ArraysCpu->Free(Pressc);       Pressc=NULL;
  ArraysCpu->Free(SpsGradvelc);  SpsGradvelc=NULL;
}

//==============================================================================
/// Perform interaction between particles. Bound-Fluid/Float
/// Realiza interaccion entre particulas. Bound-Fluid/Float
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> void JSphCpu::InteractionForcesBound
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
  ,float &viscdt,float *ar)const
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
    const tfloat4 velrhop1=velrhop[p1];

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
            tfloat4 velrhop2=velrhop[p2];
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
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift> 
  void JSphCpu::InteractionForcesFluid(unsigned n,unsigned pinit,bool boundp2,float visco
  ,StDivDataCpu divdata,const unsigned *dcell
  ,const tsymatrix3f* tau,tsymatrix3f* gradvel
  ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
  ,const float *press,const tfloat3 *dengradcorr
  ,float &viscdt,float *ar,tfloat3 *ace,float *delta
  ,TpShifting shiftmode,tfloat4 *shiftposfs)const
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
    tsymatrix3f gradvelp1={0,0,0,0,0,0};

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
    const tfloat3 velp1=TFloat3(velrhop[p1].x,velrhop[p1].y,velrhop[p1].z);
    const float rhopp1=velrhop[p1].w;
    const float pressp1=press[p1];
    const tsymatrix3f taup1=(tvisco==VISCO_Artificial? gradvelp1: tau[p1]);
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>

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

          tfloat4 velrhop2=velrhop[p2];
          if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

          //-Velocity derivative (Momentum equation).
          if(compute){
            const float prs=(pressp1+press[p2])/(rhopp1*velrhop2.w) + (tker==KERNEL_Cubic? fsph::GetKernelCubic_Tensil(CSP,rr2,rhopp1,pressp1,velrhop2.w,press[p2]): 0);
            const float p_vpm=-prs*massp2;
            acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
          }

          //-Density derivative (Continuity equation).
          const float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
          if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz)*(rhopp1/velrhop2.w);

          const float cbar=(float)Cs0;
          //-Density Diffusion Term (Molteni and Colagrossi 2009).
          if(tdensity==DDT_DDT && deltap1!=FLT_MAX){
            const float rhop1over2=rhopp1/velrhop2.w;
            const float visc_densi=DDTkh*cbar*(rhop1over2-1.f)/(rr2+Eta2);
            const float dot3=(drx*frx+dry*fry+drz*frz);
            const float delta=visc_densi*dot3*massp2;
            //deltap1=(boundp2? FLT_MAX: deltap1+delta);
            deltap1=(boundp2 && TBoundary==BC_DBC? FLT_MAX: deltap1+delta);
          }
          //-Density Diffusion Term (Fourtakas et al 2019).
          if((tdensity==DDT_DDT2 || (tdensity==DDT_DDT2Full && !boundp2)) && deltap1!=FLT_MAX && !ftp2){
            const float rh=1.f+DDTgz*drz;
            const float drhop=RhopZero*pow(rh,1.f/Gamma)-RhopZero;    
            const float visc_densi=DDTkh*cbar*((velrhop2.w-rhopp1)-drhop)/(rr2+Eta2);
            const float dot3=(drx*frx+dry*fry+drz*frz);
            const float delta=visc_densi*dot3*massp2/velrhop2.w;
            deltap1=(boundp2? FLT_MAX: deltap1-delta); //-blocks it makes it boil - bloody DBC
          }
          
          //-Shifting correction.
          if(shift && shiftposfsp1.x!=FLT_MAX){
            const float massrhop=massp2/velrhop2.w;
            const bool noshift=(boundp2 && (shiftmode==SHIFT_NoBound || (shiftmode==SHIFT_NoFixed && CODE_IsFixed(code[p2]))));
            shiftposfsp1.x=(noshift? FLT_MAX: shiftposfsp1.x+massrhop*frx); //-For boundary do not use shifting. | Con boundary anula shifting.
            shiftposfsp1.y+=massrhop*fry;
            shiftposfsp1.z+=massrhop*frz;
            shiftposfsp1.w-=massrhop*(drx*frx+dry*fry+drz*frz);
          }

          //===== Viscosity ===== 
          if(compute){
            const float dot=drx*dvx + dry*dvy + drz*dvz;
            const float dot_rr2=dot/(rr2+Eta2);
            visc=max(dot_rr2,visc);
            if(tvisco==VISCO_Artificial){//-Artificial viscosity.
              if(dot<0){
                const float amubar=KernelH*dot_rr2;  //amubar=CTE.h*dot/(rr2+CTE.eta2);
                const float robar=(rhopp1+velrhop2.w)*0.5f;
                const float pi_visc=(-visco*cbar*amubar/robar)*massp2;
                acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
              }
            }
            else if(tvisco==VISCO_LaminarSPS){//-Laminar+SPS viscosity. 
              {//-Laminar contribution.
                const float robar2=(rhopp1+velrhop2.w);
                const float temp=4.f*visco/((rr2+Eta2)*robar2);  //-Simplification of: temp=2.0f*visco/((rr2+CTE.eta2)*robar); robar=(rhopp1+velrhop2.w)*0.5f;
                const float vtemp=massp2*temp*(drx*frx+dry*fry+drz*frz);  
                acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
              }
              //-SPS turbulence model.
              float tau_xx=taup1.xx,tau_xy=taup1.xy,tau_xz=taup1.xz; //-taup1 is always zero when p1 is not a fluid particle. | taup1 siempre es cero cuando p1 no es fluid.
              float tau_yy=taup1.yy,tau_yz=taup1.yz,tau_zz=taup1.zz;
              if(!boundp2 && !ftp2){//-When p2 is a fluid particle. 
                tau_xx+=tau[p2].xx; tau_xy+=tau[p2].xy; tau_xz+=tau[p2].xz;
                tau_yy+=tau[p2].yy; tau_yz+=tau[p2].yz; tau_zz+=tau[p2].zz;
              }
              acep1.x+=massp2*(tau_xx*frx + tau_xy*fry + tau_xz*frz);
              acep1.y+=massp2*(tau_xy*frx + tau_yy*fry + tau_yz*frz);
              acep1.z+=massp2*(tau_xz*frx + tau_yz*fry + tau_zz*frz);
              //-Velocity gradients.
              if(!ftp1){//-When p1 is a fluid particle. 
                const float volp2=-massp2/velrhop2.w;
                float dv=dvx*volp2; gradvelp1.xx+=dv*frx; gradvelp1.xy+=dv*fry; gradvelp1.xz+=dv*frz;
                      dv=dvy*volp2; gradvelp1.xy+=dv*frx; gradvelp1.yy+=dv*fry; gradvelp1.yz+=dv*frz;
                      dv=dvz*volp2; gradvelp1.xz+=dv*frx; gradvelp1.yz+=dv*fry; gradvelp1.zz+=dv*frz;
                //-To compute tau terms we assume that gradvel.xy=gradvel.dudy+gradvel.dvdx, gradvel.xz=gradvel.dudz+gradvel.dwdx, gradvel.yz=gradvel.dvdz+gradvel.dwdy
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
        gradvel[p1].xx+=gradvelp1.xx;
        gradvel[p1].xy+=gradvelp1.xy;
        gradvel[p1].xz+=gradvelp1.xz;
        gradvel[p1].yy+=gradvelp1.yy;
        gradvel[p1].yz+=gradvelp1.yz;
        gradvel[p1].zz+=gradvelp1.zz;
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
void JSphCpu::InteractionForcesDEM(unsigned nfloat,StDivDataCpu divdata,const unsigned *dcell
  ,const unsigned *ftridp,const StDemData* demdata
  ,const tdouble3 *pos,const tfloat4 *velrhop
  ,const typecode *code,const unsigned *idp
  ,float &viscdt,tfloat3 *ace)const
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
            //const StDemData *demp2=demobjs+CODE_GetTypeAndValue(code[p2]);

            const float nu_mass=(!tpfluid? masstotp1/2: masstotp1*masstotp2/(masstotp1+masstotp2)); //-Con boundary toma la propia masa del floating 1.
            const float kn=4/(3*(taup1+taup2))*sqrt(float(Dp)/4); //-Generalized rigidity - Lemieux 2008.
            const float dvx=velrhop[p1].x-velrhop[p2].x, dvy=velrhop[p1].y-velrhop[p2].y, dvz=velrhop[p1].z-velrhop[p2].z; //vji
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
/// Computes sub-particle stress tensor (Tau) for SPS turbulence model.   
//==============================================================================
void JSphCpu::ComputeSpsTau(unsigned n,unsigned pini,const tfloat4 *velrhop,const tsymatrix3f *gradvel,tsymatrix3f *tau)const{
  const int pfin=int(pini+n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=int(pini);p<pfin;p++){
    const tsymatrix3f gradvel=SpsGradvelc[p];
    const float pow1=gradvel.xx*gradvel.xx + gradvel.yy*gradvel.yy + gradvel.zz*gradvel.zz;
    const float prr=pow1+pow1 + gradvel.xy*gradvel.xy + gradvel.xz*gradvel.xz + gradvel.yz*gradvel.yz;
    const float visc_sps=SpsSmag*sqrt(prr);
    const float div_u=gradvel.xx+gradvel.yy+gradvel.zz;
    const float sps_k=(2.0f/3.0f)*visc_sps*div_u;
    const float sps_blin=SpsBlin*prr;
    const float sumsps=-(sps_k+sps_blin);
    const float twovisc_sps=(visc_sps+visc_sps);
    const float one_rho2=1.0f/velrhop[p].w;   
    tau[p].xx=one_rho2*(twovisc_sps*gradvel.xx +sumsps);
    tau[p].xy=one_rho2*(visc_sps   *gradvel.xy);
    tau[p].xz=one_rho2*(visc_sps   *gradvel.xz);
    tau[p].yy=one_rho2*(twovisc_sps*gradvel.yy +sumsps);
    tau[p].yz=one_rho2*(visc_sps   *gradvel.yz);
    tau[p].zz=one_rho2*(twovisc_sps*gradvel.zz +sumsps);
  }
}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid (forces and DEM).
/// Interaccion Fluid-Fluid/Bound & Bound-Fluid (forces and DEM).
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
  void JSphCpu::Interaction_ForcesCpuT(const stinterparmsc &t,StInterResultc &res)const
{
  float viscdt=res.viscdt;
  if(t.npf){
    //-Interaction Fluid-Fluid.
    InteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift> (t.npf,t.npb,false,Visco                 
      ,t.divdata,t.dcell,t.spstau,t.spsgradvel,t.pos,t.velrhop,t.code,t.idp,t.press,t.dengradcorr
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);
    //-Interaction Fluid-Bound.
    InteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift> (t.npf,t.npb,true ,Visco*ViscoBoundFactor
      ,t.divdata,t.dcell,t.spstau,t.spsgradvel,t.pos,t.velrhop,t.code,t.idp,t.press,NULL
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);

    //-Interaction of DEM Floating-Bound & Floating-Floating. //(DEM)
    if(UseDEM)InteractionForcesDEM(CaseNfloat,t.divdata,t.dcell
      ,FtRidp,DemData,t.pos,t.velrhop,t.code,t.idp,viscdt,t.ace);

    //-Computes tau for Laminar+SPS.
    if(tvisco==VISCO_LaminarSPS)ComputeSpsTau(t.npf,t.npb,t.velrhop,t.spsgradvel,t.spstau);
  }
  if(t.npbok){
    //-Interaction Bound-Fluid.
    InteractionForcesBound<tker,ftmode> (t.npbok,0,t.divdata,t.dcell
      ,t.pos,t.velrhop,t.code,t.idp,viscdt,t.ar);
  }
  res.viscdt=viscdt;
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity> void JSphCpu::Interaction_Forces_ct5(const stinterparmsc &t,StInterResultc &res)const{
  if(Shifting)Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,true >(t,res);
  else        Interaction_ForcesCpuT<tker,ftmode,tvisco,tdensity,false>(t,res);
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco> void JSphCpu::Interaction_Forces_ct4(const stinterparmsc &t,StInterResultc &res)const{
       if(TDensity==DDT_None)    Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_None    >(t,res);
  else if(TDensity==DDT_DDT)     Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_DDT     >(t,res);
  else if(TDensity==DDT_DDT2)    Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_DDT2    >(t,res);
  else if(TDensity==DDT_DDT2Full)Interaction_Forces_ct5<tker,ftmode,tvisco,DDT_DDT2Full>(t,res);
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> void JSphCpu::Interaction_Forces_ct3(const stinterparmsc &t,StInterResultc &res)const{
       if(TVisco==VISCO_Artificial)Interaction_Forces_ct4<tker,ftmode,VISCO_Artificial>(t,res);
  else if(TVisco==VISCO_LaminarSPS)Interaction_Forces_ct4<tker,ftmode,VISCO_LaminarSPS>(t,res);
}
//==============================================================================
template<TpKernel tker> void JSphCpu::Interaction_Forces_ct2(const stinterparmsc &t,StInterResultc &res)const{
       if(FtMode==FTMODE_None)Interaction_Forces_ct3<tker,FTMODE_None>(t,res);
  else if(FtMode==FTMODE_Sph )Interaction_Forces_ct3<tker,FTMODE_Sph >(t,res);
  else if(FtMode==FTMODE_Ext )Interaction_Forces_ct3<tker,FTMODE_Ext >(t,res);
}
//==============================================================================
void JSphCpu::Interaction_Forces_ct(const stinterparmsc &t,StInterResultc &res)const{
       if(TKernel==KERNEL_Wendland)  Interaction_Forces_ct2<KERNEL_Wendland>(t,res);
  else if(TKernel==KERNEL_Cubic)     Interaction_Forces_ct2<KERNEL_Cubic   >(t,res);
}

//==============================================================================
/// Perform interaction between ghost nodes of boundaries and fluid.
//==============================================================================
template<TpKernel tker,bool sim2d,TpSlipMode tslip> void JSphCpu::InteractionMdbcCorrectionT2
  (unsigned n,StDivDataCpu divdata,float determlimit,float mdbcthreshold
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp
  ,const tfloat3 *boundnormal,const tfloat3 *motionvel,tfloat4 *velrhop)
{
  if(tslip==SLIP_FreeSlip)Run_Exceptioon("SlipMode=\'Free slip\' is not yet implemented...");
  const int nn=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<nn;p1++)if(boundnormal[p1]!=TFloat3(0)){
    float rhopfinal=FLT_MAX;
    tfloat3 velrhopfinal=TFloat3(0);
    float sumwab=0;

    //-Calculates ghost node position.
    tdouble3 gposp1=pos[p1]+ToTDouble3(boundnormal[p1]);
    gposp1=(PeriActive!=0? UpdatePeriodicPos(gposp1): gposp1); //-Corrected interface Position.
    //-Initializes variables for calculation.
    float rhopp1=0;
    tfloat3 gradrhopp1=TFloat3(0);
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
          const tfloat4 velrhopp2=velrhop[p2];
          const float massp2=MassFluid;
          const float volp2=massp2/velrhopp2.w;

          //===== Density and its gradient =====
          rhopp1+=massp2*wab;
          gradrhopp1.x+=massp2*frx;
          gradrhopp1.y+=massp2*fry;
          gradrhopp1.z+=massp2*frz;

          //===== Kernel values multiplied by volume =====
          const float vwab=wab*volp2;
          sumwab+=vwab;
          const float vfrx=frx*volp2;
          const float vfry=fry*volp2;
          const float vfrz=frz*volp2;

          //===== Velocity =====
          if(tslip!=SLIP_Vel0){
            velp1.x+=vwab*velrhopp2.x;
            velp1.y+=vwab*velrhopp2.y;
            velp1.z+=vwab*velrhopp2.z;
          }

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
    if(sumwab>=mdbcthreshold){
      const tfloat3 dpos=(boundnormal[p1]*(-1.f)); //-Boundary particle position - ghost node position.
      if(sim2d){
        const double determ=fmath::Determinant3x3(a_corr2);
        if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
          const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
          const float rhoghost=float(invacorr2.a11*rhopp1 + invacorr2.a12*gradrhopp1.x + invacorr2.a13*gradrhopp1.z);
          const float grx=    -float(invacorr2.a21*rhopp1 + invacorr2.a22*gradrhopp1.x + invacorr2.a23*gradrhopp1.z);
          const float grz=    -float(invacorr2.a31*rhopp1 + invacorr2.a32*gradrhopp1.x + invacorr2.a33*gradrhopp1.z);
          rhopfinal=(rhoghost + grx*dpos.x + grz*dpos.z);
        }
        else if(a_corr2.a11>0){//-Determinant is small but a11 is nonzero (0th order).
          rhopfinal=float(rhopp1/a_corr2.a11);
        }
        //-Ghost node velocity (0th order).
        if(tslip!=SLIP_Vel0){
          velrhopfinal.x=float(velp1.x/a_corr2.a11);
          velrhopfinal.z=float(velp1.z/a_corr2.a11);
          velrhopfinal.y=0;
        }
      }
      else{
        const double determ=fmath::Determinant4x4(a_corr3);
        if(fabs(determ)>=determlimit){
          const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
          const float rhoghost=float(invacorr3.a11*rhopp1 + invacorr3.a12*gradrhopp1.x + invacorr3.a13*gradrhopp1.y + invacorr3.a14*gradrhopp1.z);
          const float grx=    -float(invacorr3.a21*rhopp1 + invacorr3.a22*gradrhopp1.x + invacorr3.a23*gradrhopp1.y + invacorr3.a24*gradrhopp1.z);
          const float gry=    -float(invacorr3.a31*rhopp1 + invacorr3.a32*gradrhopp1.x + invacorr3.a33*gradrhopp1.y + invacorr3.a34*gradrhopp1.z);
          const float grz=    -float(invacorr3.a41*rhopp1 + invacorr3.a42*gradrhopp1.x + invacorr3.a43*gradrhopp1.y + invacorr3.a44*gradrhopp1.z);
          rhopfinal=(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
        }
        else if(a_corr3.a11>0){//-Determinant is small but a11 is nonzero (0th order).
          rhopfinal=float(rhopp1/a_corr3.a11);
        }
        //-Ghost node velocity (0th order).
        if(tslip!=SLIP_Vel0){
          velrhopfinal.x=float(velp1.x/a_corr3.a11);
          velrhopfinal.y=float(velp1.y/a_corr3.a11);
          velrhopfinal.z=float(velp1.z/a_corr3.a11);
        }
      }
      //-Store the results.
      rhopfinal=(rhopfinal!=FLT_MAX? rhopfinal: RhopZero);
      if(tslip==SLIP_Vel0){//-DBC vel=0
        velrhop[p1].w=rhopfinal;
      }
      if(tslip==SLIP_NoSlip){//-No-Slip
        const tfloat3 v=motionvel[p1];
        velrhop[p1]=TFloat4(v.x+v.x-velrhopfinal.x,v.y+v.y-velrhopfinal.y,v.z+v.z-velrhopfinal.z,rhopfinal);
      }
      if(tslip==SLIP_FreeSlip){//-No-Penetration and free slip    SHABA

		tfloat3 FSVelFinal; // final free slip boundary velocity
		const tfloat3 v = motionvel[p1];
		float motion = sqrt(v.x*v.x + v.y*v.y + v.z*v.z); // to check if boundary moving
		float norm = sqrt(boundnormal[p1].x*boundnormal[p1].x + boundnormal[p1].y*boundnormal[p1].y + boundnormal[p1].z*boundnormal[p1].z);
		tfloat3 normal; // creating a normailsed boundary normal
		normal.x = fabs(boundnormal[p1].x )/ norm; normal.y = fabs(boundnormal[p1].y) / norm; normal.z = fabs(boundnormal[p1].z) / norm;
		
		// finding the velocity componants normal and tangential to boundary 
		tfloat3 normvel = TFloat3(velrhopfinal.x*normal.x, velrhopfinal.y*normal.y, velrhopfinal.z*normal.z); // velocity in direction of normal pointin ginto fluid)
		tfloat3 tangvel = TFloat3(velrhopfinal.x - normvel.x, velrhopfinal.y - normvel.y, velrhopfinal.z - normvel.z); // velocity tangential to normal
		
		if (motion > 0.f) { // if moving boundary
			tfloat3 normmot = TFloat3(v.x*normal.x, v.y*normal.y, v.z*normal.z); // boundary motion in direction normal to boundary 
			FSVelFinal = TFloat3(normmot.x+normmot.x-normvel.x, normmot.y + normmot.y -normvel.y, normmot.z + normmot.z -normvel.z);
			// only velocity in normal direction for no-penetration
			// fluid sees zero velocity in the tangetial direction
		}
		else {
			FSVelFinal = TFloat3(tangvel.x - normvel.x, tangvel.y - normvel.y, tangvel.z - normvel.z);
			// tangential velocity equal to fluid velocity for free slip
			// normal velocity reversed for no-penetration
		}
		
		// Save the velocity and density
		velrhop[p1]=TFloat4(FSVelFinal.x, FSVelFinal.y, FSVelFinal.z,rhopfinal); 

      }
    }
  }
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
 template<TpKernel tker> void JSphCpu::Interaction_MdbcCorrectionT(TpSlipMode slipmode
  ,const StDivDataCpu &divdata,const tdouble3 *pos,const typecode *code,const unsigned *idp
  ,const tfloat3 *boundnormal,const tfloat3 *motionvel,tfloat4 *velrhop)
{
  const float determlimit=1e-3f;
  //-Interaction GhostBoundaryNodes-Fluid.
  unsigned n=NpbOk;
  if(Simulate2D){ const bool sim2d=true;
    if(slipmode==SLIP_Vel0    )InteractionMdbcCorrectionT2 <tker,sim2d,SLIP_Vel0    > (n,divdata,determlimit,MdbcThreshold,pos,code,idp,boundnormal,motionvel,velrhop);
    if(slipmode==SLIP_NoSlip  )InteractionMdbcCorrectionT2 <tker,sim2d,SLIP_NoSlip  > (n,divdata,determlimit,MdbcThreshold,pos,code,idp,boundnormal,motionvel,velrhop);
    if(slipmode==SLIP_FreeSlip)InteractionMdbcCorrectionT2 <tker,sim2d,SLIP_FreeSlip> (n,divdata,determlimit,MdbcThreshold,pos,code,idp,boundnormal,motionvel,velrhop);
  }else{          const bool sim2d=false;
    if(slipmode==SLIP_Vel0    )InteractionMdbcCorrectionT2 <tker,sim2d,SLIP_Vel0    > (n,divdata,determlimit,MdbcThreshold,pos,code,idp,boundnormal,motionvel,velrhop);
    if(slipmode==SLIP_NoSlip  )InteractionMdbcCorrectionT2 <tker,sim2d,SLIP_NoSlip  > (n,divdata,determlimit,MdbcThreshold,pos,code,idp,boundnormal,motionvel,velrhop);
    if(slipmode==SLIP_FreeSlip)InteractionMdbcCorrectionT2 <tker,sim2d,SLIP_FreeSlip> (n,divdata,determlimit,MdbcThreshold,pos,code,idp,boundnormal,motionvel,velrhop);
  }
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphCpu::Interaction_MdbcCorrection(TpSlipMode slipmode,const StDivDataCpu &divdata
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp
  ,const tfloat3 *boundnormal,const tfloat3 *motionvel,tfloat4 *velrhop)
{
  switch(TKernel){
    case KERNEL_Cubic:       Interaction_MdbcCorrectionT <KERNEL_Cubic     > (slipmode,divdata,pos,code,idp,boundnormal,motionvel,velrhop);  break;
    case KERNEL_Wendland:    Interaction_MdbcCorrectionT <KERNEL_Wendland  > (slipmode,divdata,pos,code,idp,boundnormal,motionvel,velrhop);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}


//==============================================================================
/// Update pos, dcell and code to move with indicated displacement.
/// The value of outrhop indicates is it outside of the density limits.
/// Check the limits in funcion of MapRealPosMin & MapRealSize that this is valid
/// for single-cpu because DomRealPos & MapRealPos are equal. For multi-cpu it will be 
/// necessary to mark the particles that leave the domain without leaving the map.
///
/// Actualiza pos, dcell y code a partir del desplazamiento indicado.
/// El valor de outrhop indica si esta fuera de los limites de densidad.
/// Comprueba los limites en funcion de MapRealPosMin y MapRealSize esto es valido
/// para single-cpu pq DomRealPos y MapRealPos son iguales. Para multi-cpu seria 
/// necesario marcar las particulas q salgan del dominio sin salir del mapa.
//==============================================================================
void JSphCpu::UpdatePos(tdouble3 rpos,double movx,double movy,double movz
  ,bool outrhop,unsigned p,tdouble3 *pos,unsigned *cell,typecode *code)const
{
  //-Check validity of displacement. | Comprueba validez del desplazamiento.
  bool outmove=(fabs(float(movx))>MovLimit || fabs(float(movy))>MovLimit || fabs(float(movz))>MovLimit);
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
  if(outrhop || outmove || out){//-Particle out.
    typecode rcode=code[p];
    if(out)rcode=CODE_SetOutPos(rcode);
    else if(outrhop)rcode=CODE_SetOutRhop(rcode);
    else rcode=CODE_SetOutMove(rcode);
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
    cell[p]=PC__Cell(DomCellCode,cx,cy,cz);
  }
}


//==============================================================================
/// Calculate new values of position, velocity & density for fluid (using Verlet).
/// Calcula nuevos valores de posicion, velocidad y densidad para el fluido (usando Verlet).
//==============================================================================
void JSphCpu::ComputeVerletVarsFluid(bool shift,const tfloat3 *indirvel
  ,const tfloat4 *velrhop1,const tfloat4 *velrhop2,double dt,double dt2
  ,tdouble3 *pos,unsigned *dcell,typecode *code,tfloat4 *velrhopnew)const
{
  const double dt205=0.5*dt*dt;
  const tdouble3 gravity=ToTDouble3(Gravity);
  const int pini=int(Npb),pfin=int(Np),npf=int(Np-Npb);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=pini;p<pfin;p++){
    //-Calculate density. | Calcula densidad.
    const float rhopnew=float(double(velrhop2[p].w)+dt2*Arc[p]);
    if(!WithFloating || CODE_IsFluid(code[p])){//-Fluid Particles.
      const tdouble3 acegr=ToTDouble3(Acec[p])+gravity; //-Adds gravity.
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=double(velrhop1[p].x)*dt + acegr.x*dt205;
      double dy=double(velrhop1[p].y)*dt + acegr.y*dt205;
      double dz=double(velrhop1[p].z)*dt + acegr.z*dt205;
      if(shift){
        dx+=double(ShiftPosfsc[p].x);
        dy+=double(ShiftPosfsc[p].y);
        dz+=double(ShiftPosfsc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin || rhopnew>RhopOutMax);
      //-Calculate velocity & density. | Calcula velocidad y densidad.
      tfloat4 rvelrhopnew=TFloat4(
        float(double(velrhop2[p].x) + acegr.x*dt2),
        float(double(velrhop2[p].y) + acegr.y*dt2),
        float(double(velrhop2[p].z) + acegr.z*dt2),
        rhopnew);
      //-Restore data of inout particles.
      if(InOut && CODE_IsFluidInout(Codec[p])){
        outrhop=false;
        rvelrhopnew=velrhop2[p];
        const tfloat3 vd=indirvel[CODE_GetIzoneFluidInout(Codec[p])];
        if(vd.x!=FLT_MAX){
          const float v=velrhop1[p].x*vd.x + velrhop1[p].y*vd.y + velrhop1[p].z*vd.z;
          dx=double(v*vd.x) * dt;
          dy=double(v*vd.y) * dt;
          dz=double(v*vd.z) * dt;
        }
        else{
          dx=double(velrhop1[p].x) * dt;
          dy=double(velrhop1[p].y) * dt;
          dz=double(velrhop1[p].z) * dt;
        }
      }
      //-Update particle data.
      UpdatePos(pos[p],dx,dy,dz,outrhop,p,pos,dcell,code);
      velrhopnew[p]=rvelrhopnew;
    }
    else{//-Floating Particles.
      velrhopnew[p]=velrhop1[p];
      velrhopnew[p].w=(rhopnew<RhopZero? RhopZero: rhopnew); //-Avoid fluid particles being absorved by floating ones. | Evita q las floating absorvan a las fluidas.
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
void JSphCpu::ComputeVelrhopBound(const tfloat4* velrhopold,double armul,tfloat4* velrhopnew)const{
  const int npb=int(Npb);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<npb;p++){
    const float rhopnew=float(double(velrhopold[p].w)+armul*Arc[p]);
    velrhopnew[p]=TFloat4(0,0,0,(rhopnew<RhopZero? RhopZero: rhopnew));//-Avoid fluid particles being absorved by boundary ones. | Evita q las boundary absorvan a las fluidas.
  }
}

//==============================================================================
/// Update of particles according to forces and dt using Verlet.
/// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphCpu::ComputeVerlet(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  const bool shift=(Shifting!=NULL);
  const tfloat3 *indirvel=(InOut? InOut->GetDirVel(): NULL);
  VerletStep++;
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    ComputeVerletVarsFluid(shift,indirvel,Velrhopc,VelrhopM1c,dt,twodt,Posc,Dcellc,Codec,VelrhopM1c);
    ComputeVelrhopBound(VelrhopM1c,twodt,VelrhopM1c);
  }
  else{
    ComputeVerletVarsFluid(shift,indirvel,Velrhopc,Velrhopc,dt,dt,Posc,Dcellc,Codec,VelrhopM1c);
    ComputeVelrhopBound(Velrhopc,dt,VelrhopM1c);
    VerletStep=0;
  }
  //-New values are calculated en VelrhopM1c. | Los nuevos valores se calculan en VelrhopM1c.
  swap(Velrhopc,VelrhopM1c);     //-Swap Velrhopc & VelrhopM1c. | Intercambia Velrhopc y VelrhopM1c.
  TmcStop(Timers,TMC_SuComputeStep);
}


//==============================================================================
/// Update of particles according to forces and dt using Symplectic-Predictor.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
void JSphCpu::ComputeSymplecticPre(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  const bool shift=false; //(ShiftingMode!=SHIFT_None); //-We strongly recommend running the shifting correction only for the corrector. If you want to re-enable shifting in the predictor, change the value here to "true".
  const double dt05=dt*.5;
  const int np=int(Np);
  const int npb=int(Npb);
  const int npf=np-npb;

  //-Assign memory to variables Pre. | Asigna memoria a variables Pre.
  PosPrec=ArraysCpu->ReserveDouble3();
  VelrhopPrec=ArraysCpu->ReserveFloat4();
  //-Change data to variables Pre to calculate new data. | Cambia datos a variables Pre para calcular nuevos datos.
  swap(PosPrec,Posc);         //Put value of Pos[] in PosPre[].         | Es decir... PosPre[] <= Pos[].
  swap(VelrhopPrec,Velrhopc); //Put value of Velrhop[] in VelrhopPre[]. | Es decir... VelrhopPre[] <= Velrhop[].

  //-Calculate new density for boundary and copy velocity. | Calcula nueva densidad para el contorno y copia velocidad.
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<npb;p++){
    const tfloat4 vr=VelrhopPrec[p];
    const float rhopnew=float(double(vr.w)+dt05*Arc[p]);
    Velrhopc[p]=TFloat4(vr.x,vr.y,vr.z,(rhopnew<RhopZero? RhopZero: rhopnew));//-Avoid fluid particles being absorbed by boundary ones. | Evita q las boundary absorvan a las fluidas.
  }

  //-Compute displacement, velocity and density for fluid.
  tdouble3 *movc=ArraysCpu->ReserveDouble3();
  const tfloat3 *indirvel=(InOut? InOut->GetDirVel(): NULL);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=npb;p<np;p++){
    const typecode rcode=Codec[p];
    //-Calculate density.
    const float rhopnew=float(double(VelrhopPrec[p].w)+dt05*Arc[p]);
    if(!WithFloating || CODE_IsFluid(rcode)){//-Fluid Particles.
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=double(VelrhopPrec[p].x)*dt05;
      double dy=double(VelrhopPrec[p].y)*dt05;
      double dz=double(VelrhopPrec[p].z)*dt05;
      if(shift){
        dx+=double(ShiftPosfsc[p].x);
        dy+=double(ShiftPosfsc[p].y);
        dz+=double(ShiftPosfsc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin || rhopnew>RhopOutMax);
      //-Calculate velocity & density. | Calcula velocidad y densidad.
      tfloat4 rvelrhopnew=TFloat4(
        float(double(VelrhopPrec[p].x) + (double(Acec[p].x)+Gravity.x) * dt05),
        float(double(VelrhopPrec[p].y) + (double(Acec[p].y)+Gravity.y) * dt05),
        float(double(VelrhopPrec[p].z) + (double(Acec[p].z)+Gravity.z) * dt05),
        rhopnew);
      //-Restore data of inout particles.
      if(InOut && CODE_IsFluidInout(rcode)){
        outrhop=false;
        rvelrhopnew=VelrhopPrec[p];
        const tfloat3 vd=indirvel[CODE_GetIzoneFluidInout(rcode)];
        if(vd.x!=FLT_MAX){
          const float v=rvelrhopnew.x*vd.x + rvelrhopnew.y*vd.y + rvelrhopnew.z*vd.z;
          dx=double(v*vd.x) * dt05;
          dy=double(v*vd.y) * dt05;
          dz=double(v*vd.z) * dt05;
        }
      }
      //-Update particle data.
      movc[p]=TDouble3(dx,dy,dz);
      Velrhopc[p]=rvelrhopnew;
      if(outrhop && CODE_IsNormal(rcode))Codec[p]=CODE_SetOutRhop(rcode); //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
    }
    else{//-Floating Particles.
      Velrhopc[p]=VelrhopPrec[p];
      Velrhopc[p].w=(rhopnew<RhopZero? RhopZero: rhopnew); //-Avoid fluid particles being absorbed by floating ones. | Evita q las floating absorvan a las fluidas.
    }
  }

  //-Applies displacement to non-periodic fluid particles.
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=npb;p<np;p++){
    const typecode rcode=Codec[p];
    const bool outrhop=CODE_IsOutRhop(rcode);
    const bool normal=(!PeriActive || outrhop || CODE_IsNormal(rcode));
    if(normal){//-Does not apply to periodic particles. | No se aplica a particulas periodicas
      if(CODE_IsFluid(rcode)){//-Only applied for fluid displacement. | Solo se aplica desplazamiento al fluido.
        UpdatePos(PosPrec[p],movc[p].x,movc[p].y,movc[p].z,outrhop,p,Posc,Dcellc,Codec);
      }
      else Posc[p]=PosPrec[p]; //-Copy position of floating particles.
    }
  }
  
  //-Frees memory allocated for the displacement.
  ArraysCpu->Free(movc);   movc=NULL;

  //-Copy previous position of boundary. | Copia posicion anterior del contorno.
  memcpy(Posc,PosPrec,sizeof(tdouble3)*Npb);

  TmcStop(Timers,TMC_SuComputeStep);
}

//==============================================================================
/// Update particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphCpu::ComputeSymplecticCorr(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  const bool shift=(Shifting!=NULL);
  const double dt05=dt*.5;
  const int np=int(Np);
  const int npb=int(Npb);
  const int npf=np-npb;
  
  //-Calculate rhop of boudary and set velocity=0. | Calcula rhop de contorno y vel igual a cero.
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<npb;p++){
    const double epsilon_rdot=(-double(Arc[p])/double(Velrhopc[p].w))*dt;
    const float rhopnew=float(double(VelrhopPrec[p].w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
    Velrhopc[p]=TFloat4(0,0,0,(rhopnew<RhopZero? RhopZero: rhopnew));//-Avoid fluid particles being absorbed by boundary ones. | Evita q las boundary absorvan a las fluidas.
  }

  //-Compute displacement, velocity and density for fluid.
  tdouble3 *movc=ArraysCpu->ReserveDouble3();
  const tfloat3 *indirvel=(InOut? InOut->GetDirVel(): NULL);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=npb;p<np;p++){
    const typecode rcode=Codec[p];
    const double epsilon_rdot=(-double(Arc[p])/double(Velrhopc[p].w))*dt;
    const float rhopnew=float(double(VelrhopPrec[p].w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
    if(!WithFloating || CODE_IsFluid(rcode)){//-Fluid Particles.
      //-Calculate velocity & density. | Calcula velocidad y densidad.
      tfloat4 rvelrhopnew=TFloat4(
        float(double(VelrhopPrec[p].x) + (double(Acec[p].x)+Gravity.x) * dt), 
        float(double(VelrhopPrec[p].y) + (double(Acec[p].y)+Gravity.y) * dt), 
        float(double(VelrhopPrec[p].z) + (double(Acec[p].z)+Gravity.z) * dt),
        rhopnew);
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=(double(VelrhopPrec[p].x)+double(rvelrhopnew.x)) * dt05; 
      double dy=(double(VelrhopPrec[p].y)+double(rvelrhopnew.y)) * dt05; 
      double dz=(double(VelrhopPrec[p].z)+double(rvelrhopnew.z)) * dt05;
      if(shift){
        dx+=double(ShiftPosfsc[p].x);
        dy+=double(ShiftPosfsc[p].y);
        dz+=double(ShiftPosfsc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin || rhopnew>RhopOutMax);
      //-Restore data of inout particles.
      if(InOut && CODE_IsFluidInout(rcode)){
        outrhop=false;
        rvelrhopnew=VelrhopPrec[p];
        const tfloat3 vd=indirvel[CODE_GetIzoneFluidInout(rcode)];
        if(vd.x!=FLT_MAX){
          const float v=rvelrhopnew.x*vd.x + rvelrhopnew.y*vd.y + rvelrhopnew.z*vd.z;
          dx=double(v*vd.x) * dt;
          dy=double(v*vd.y) * dt;
          dz=double(v*vd.z) * dt;
        }
        else{
          dx=double(rvelrhopnew.x) * dt; 
          dy=double(rvelrhopnew.y) * dt; 
          dz=double(rvelrhopnew.z) * dt;
        }
      }
      //-Update particle data.
      movc[p]=TDouble3(dx,dy,dz);
      if(outrhop && CODE_IsNormal(rcode))Codec[p]=CODE_SetOutRhop(rcode); //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
      Velrhopc[p]=rvelrhopnew;
    }
    else{//-Floating Particles.
      Velrhopc[p]=VelrhopPrec[p];
      Velrhopc[p].w=(rhopnew<RhopZero? RhopZero: rhopnew); //-Avoid fluid particles being absorbed by floating ones. | Evita q las floating absorvan a las fluidas.
    }
  }

  //-Applies displacement to non-periodic fluid particles.
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=npb;p<np;p++){
    const typecode rcode=Codec[p];
    const bool outrhop=CODE_IsOutRhop(rcode);
    const bool normal=(!PeriActive || outrhop || CODE_IsNormal(rcode));
    if(normal){//-Does not apply to periodic particles. | No se aplica a particulas periodicas
      if(CODE_IsFluid(rcode)){//-Only applied for fluid displacement. | Solo se aplica desplazamiento al fluido.
        UpdatePos(PosPrec[p],movc[p].x,movc[p].y,movc[p].z,outrhop,p,Posc,Dcellc,Codec);
      }
      else Posc[p]=PosPrec[p]; //-Copy position of floating particles.
    }
  }
 
  //-Frees memory allocated for the displacement.
  ArraysCpu->Free(movc);   movc=NULL;

  //-Free memory assigned to variables Pre and ComputeSymplecticPre(). | Libera memoria asignada a variables Pre en ComputeSymplecticPre().
  ArraysCpu->Free(PosPrec);      PosPrec=NULL;
  ArraysCpu->Free(VelrhopPrec);  VelrhopPrec=NULL;
  TmcStop(Timers,TMC_SuComputeStep);
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
  double dt=double(CFLnumber)*min(dt1,dt2);
  if(FixedDt)dt=FixedDt->GetDt(TimeStep,dt);
  if(fun::IsNAN(dt) || fun::IsInfinity(dt))Run_Exceptioon(fun::PrintStr("The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f) is NaN or infinity at nstep=%u.",dt,AceMax,VelMax,ViscDtMax,Nstep));
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
    if(SaveDt)SaveDt->AddValues(TimeStep,dt,dt1*CFLnumber,dt2*CFLnumber,AceMax,ViscDtMax,VelMax);
  }
  return(dt);
}

//==============================================================================
/// Calculate final Shifting for particles' position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphCpu::RunShifting(double dt){
  TmcStart(Timers,TMC_SuShifting);
  Shifting->RunCpu(Np-Npb,Npb,dt,Velrhopc,ShiftPosfsc);
  TmcStop(Timers,TMC_SuShifting);
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
void JSphCpu::CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini,unsigned idfin,const typecode *code,const unsigned *idp,unsigned *ridp)const{
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
void JSphCpu::MoveLinBound(unsigned np,unsigned ini,const tdouble3 &mvpos,const tfloat3 &mvvel
  ,const unsigned *ridp,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code)const
{
  const unsigned fin=ini+np;
  for(unsigned id=ini;id<fin;id++){
    const unsigned pid=RidpMove[id];
    if(pid!=UINT_MAX){
      UpdatePos(pos[pid],mvpos.x,mvpos.y,mvpos.z,false,pid,pos,dcell,code);
      velrhop[pid].x=mvvel.x;  velrhop[pid].y=mvvel.y;  velrhop[pid].z=mvvel.z;
    }
  }
}

//==============================================================================
/// Applies a matrix movement to a group of particles.
/// Aplica un movimiento matricial a un conjunto de particulas.
//==============================================================================
void JSphCpu::MoveMatBound(unsigned np,unsigned ini,tmatrix4d m,double dt,const unsigned *ridpmv
  ,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code,tfloat3 *boundnormal)const
{
  const unsigned fin=ini+np;
  for(unsigned id=ini;id<fin;id++){
    const unsigned pid=RidpMove[id];
    if(pid!=UINT_MAX){
      const tdouble3 ps=pos[pid];
      tdouble3 ps2=MatrixMulPoint(m,ps);
      if(Simulate2D)ps2.y=ps.y;
      const double dx=ps2.x-ps.x, dy=ps2.y-ps.y, dz=ps2.z-ps.z;
      UpdatePos(ps,dx,dy,dz,false,pid,pos,dcell,code);
      velrhop[pid].x=float(dx/dt);  velrhop[pid].y=float(dy/dt);  velrhop[pid].z=float(dz/dt);
      //-Computes normal.
      if(boundnormal){
        const tdouble3 gs=ps+ToTDouble3(boundnormal[pid]);
        const tdouble3 gs2=MatrixMulPoint(m,gs);
        boundnormal[pid]=ToTFloat3(gs2-ps2);
      }
    }
  }
}

//==============================================================================
/// Copy motion velocity to MotionVel[].
/// Copia velocidad de movimiento a MotionVel[].
//==============================================================================
void JSphCpu::CopyMotionVel(unsigned nmoving,const unsigned *ridp
  ,const tfloat4 *velrhop,tfloat3 *motionvel)const
{
  for(unsigned id=0;id<nmoving;id++){
    const unsigned pid=RidpMove[id];
    if(pid!=UINT_MAX){
      const tfloat4 v=velrhop[pid];
      motionvel[pid]=TFloat3(v.x,v.y,v.z);
    }
  }
}

//==============================================================================
/// Calculates predefined movement of boundary particles.
/// Calcula movimiento predefinido de boundary particles.
//==============================================================================
void JSphCpu::CalcMotion(double stepdt){
  TmcStart(Timers,TMC_SuMotion);
  JSph::CalcMotion(stepdt);
  TmcStop(Timers,TMC_SuMotion);
}

//==============================================================================
/// Process movement of boundary particles.
/// Procesa movimiento de boundary particles.
//==============================================================================
void JSphCpu::RunMotion(double stepdt){
  TmcStart(Timers,TMC_SuMotion);
  tfloat3 *boundnormal=NULL;
  boundnormal=BoundNormalc;
  const bool motsim=true;
  BoundChanged=false;
  //-Add motion from automatic wave generation.
  if(WaveGen)CalcMotionWaveGen(stepdt);
  //-Process particles motion.
  if(DsMotion->GetActiveMotion()){
    CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codec,Idpc,RidpMove);
    BoundChanged=true;
    const unsigned nref=DsMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      if(m.type==MOTT_Linear){//-Linear movement.
        if(motsim)MoveLinBound   (m.count,m.idbegin-CaseNfixed,m.linmov,ToTFloat3(m.linvel),RidpMove,Posc,Dcellc,Velrhopc,Codec);
        //else    MoveLinBoundAce(m.count,m.idbegin-CaseNfixed,m.linmov,ToTFloat3(m.linvel),ToTFloat3(m.linace),RidpMove,Posc,Dcellc,Velrhopc,Acec,Codec);
      }
      if(m.type==MOTT_Matrix){//-Matrix movement (for rotations).
        if(motsim)MoveMatBound   (m.count,m.idbegin-CaseNfixed,m.matmov,stepdt,RidpMove,Posc,Dcellc,Velrhopc,Codec,boundnormal); 
        //else    MoveMatBoundAce(m.count,m.idbegin-CaseNfixed,m.matmov,m.matmov2,stepdt,RidpMove,Posc,Dcellc,Velrhopc,Acec,Codec);
      }      
      //-Applies predefined motion to BoundCorr configuration.
      if(BoundCorr && BoundCorr->GetUseMotion())BoundCorr->RunMotion(m);
    }
  }
  //-Management of Multi-Layer Pistons.
  if(MLPistons){
    if(!BoundChanged)CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codec,Idpc,RidpMove);
    BoundChanged=true;
    if(MLPistons->GetPiston1dCount()){//-Process motion for pistons 1D.
      MLPistons->CalculateMotion1d(TimeStep+MLPistons->GetTimeMod()+stepdt);
      MovePiston1d(CaseNmoving,0,MLPistons->GetPoszMin(),MLPistons->GetPoszCount()
        ,MLPistons->GetPistonId(),MLPistons->GetMovx(),MLPistons->GetVelx()
        ,RidpMove,Posc,Dcellc,Velrhopc,Codec);
    }
    for(unsigned cp=0;cp<MLPistons->GetPiston2dCount();cp++){//-Process motion for pistons 2D.
      JMLPistons::StMotionInfoPiston2D mot=MLPistons->CalculateMotion2d(cp,TimeStep+MLPistons->GetTimeMod()+stepdt);
      MovePiston2d(mot.np,mot.idbegin-CaseNfixed,mot.posymin,mot.poszmin,mot.poszcount,mot.movyz,mot.velyz
        ,RidpMove,Posc,Dcellc,Velrhopc,Codec);
    }
  }
  if(MotionVelc)CopyMotionVel(CaseNmoving,RidpMove,Velrhopc,MotionVelc);
  TmcStop(Timers,TMC_SuMotion);
}

//==============================================================================
/// Applies movement and velocity of piston 1D to a group of particles.
/// Aplica movimiento y velocidad de piston 1D a conjunto de particulas.
//==============================================================================
void JSphCpu::MovePiston1d(unsigned np,unsigned ini
  ,double poszmin,unsigned poszcount,const byte *pistonid,const double* movx,const double* velx
  ,const unsigned *ridpmv,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code)const
{
  const int fin=int(ini+np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(fin>OMP_LIMIT_LIGHT)
  #endif
  for(int id=int(ini);id<fin;id++){
    const unsigned pid=ridpmv[id];
    if(pid!=UINT_MAX){
      const unsigned pisid=pistonid[CODE_GetTypeValue(code[pid])];
      if(pisid<255){
        const unsigned cz=unsigned((pos[pid].z-poszmin)/Dp);
        const double rmovx=(cz<poszcount? movx[pisid*poszcount+cz]: 0);
        const float rvelx=float(cz<poszcount? velx[pisid*poszcount+cz]: 0);
        //-Updates position.
        UpdatePos(pos[pid],rmovx,0,0,false,pid,pos,dcell,code);
        //-Updates velocity.
        velrhop[pid].x=rvelx;
      }
    }
  }
}
//==============================================================================
/// Applies movement and velocity of piston 2D to a group of particles.
/// Aplica movimiento y velocidad de piston 2D a conjunto de particulas.
//==============================================================================
void JSphCpu::MovePiston2d(unsigned np,unsigned ini
  ,double posymin,double poszmin,unsigned poszcount,const double* movx,const double* velx
  ,const unsigned *ridpmv,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code)const
{
  const int fin=int(ini+np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(fin>OMP_LIMIT_LIGHT)
  #endif
  for(int id=int(ini);id<fin;id++){
    const unsigned pid=ridpmv[id];
    if(pid!=UINT_MAX){
      const tdouble3 ps=pos[pid];
      const unsigned cy=unsigned((ps.y-posymin)/Dp);
      const unsigned cz=unsigned((ps.z-poszmin)/Dp);
      const double rmovx=(cz<poszcount? movx[cy*poszcount+cz]: 0);
      const float rvelx=float(cz<poszcount? velx[cy*poszcount+cz]: 0);
      //-Updates position.
      UpdatePos(ps,rmovx,0,0,false,pid,pos,dcell,code);
      //-Updates velocity.
      velrhop[pid].x=rvelx;
    }
  }
}

//==============================================================================
/// Applies RelaxZone to selected particles.
/// Aplica RelaxZone a las particulas indicadas.
//==============================================================================
void JSphCpu::RunRelaxZone(double dt){
  TmcStart(Timers,TMC_SuMotion);
  RelaxZones->SetFluidVel(TimeStep,dt,Np-Npb,Npb,Posc,Idpc,Velrhopc);
  TmcStop(Timers,TMC_SuMotion);
}

//==============================================================================
/// Applies Damping to selected particles.
/// Aplica Damping a las particulas indicadas.
//==============================================================================
void JSphCpu::RunDamping(double dt,unsigned np,unsigned npb,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const{
  if(CaseNfloat || PeriActive)Damping->ComputeDamping(TimeStep,dt,np-npb,npb,pos,code,velrhop);
  else Damping->ComputeDamping(TimeStep,dt,np-npb,npb,pos,NULL,velrhop);
}

//==============================================================================
/// Adjust variables of floating body particles.
/// Ajusta variables de particulas floating body.
//==============================================================================
void JSphCpu::InitFloating(){
  if(PartBegin){
    JPartFloatBi4Load ftdata;
    ftdata.LoadFile(PartBeginDir);
    //-Check cases of constant values. | Comprueba coincidencia de datos constantes.
    for(unsigned cf=0;cf<FtCount;cf++){
      const StFloatingData &ft=FtObjs[cf];
      ftdata.CheckHeadData(cf,ft.mkbound,ft.begin,ft.count,ft.mass,ft.massp);
    }
    //-Load PART data. | Carga datos de PART.
    ftdata.LoadPart(PartBegin);
    for(unsigned cf=0;cf<FtCount;cf++){
      FtObjs[cf].center=ftdata.GetPartCenter(cf);
      FtObjs[cf].fvel  =ftdata.GetPartVelLin(cf);
      FtObjs[cf].fomega=ftdata.GetPartVelAng(cf);
      FtObjs[cf].radius=ftdata.GetHeadRadius(cf);
    }
    DemDtForce=ftdata.GetPartDemDtForce();
  }
}

//==============================================================================
/// Show active timers.
/// Muestra los temporizadores activos.
//==============================================================================
void JSphCpu::ShowTimers(bool onlyfile){
  JLog2::TpMode_Out mode=(onlyfile? JLog2::Out_File: JLog2::Out_ScrFile);
  Log->Print("[CPU Timers]",mode);
  if(!SvTimers)Log->Print("none",mode);
  else for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c))Log->Print(TimerToText(c),mode);
}

//==============================================================================
/// Return string with names and values of active timers.
/// Devuelve string con nombres y valores de los timers activos.
//==============================================================================
void JSphCpu::GetTimersInfo(std::string &hinfo,std::string &dinfo)const{
  for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c)){
    hinfo=hinfo+";"+TimerGetName(c);
    dinfo=dinfo+";"+fun::FloatStr(TimerGetValue(c)/1000.f);
  }
}


