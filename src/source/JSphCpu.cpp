//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2019 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "JCellDivCpu.h"
#include "JPartFloatBi4.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JSphMotion.h"
#include "JArraysCpu.h"
#include "JSphDtFixed.h"
#include "JWaveGen.h"
#include "JMLPistons.h"     //<vs_mlapiston>
#include "JRelaxZones.h"    //<vs_rzone>
#include "JChronoObjects.h" //<vs_chroono>
#include "JDamping.h"
#include "JXml.h"
#include "JSaveDt.h"
#include "JTimeOut.h"
#include "JSphAccInput.h"
#include "JGaugeSystem.h"
#include "JSphBoundCorr.h"  //<vs_innlet>
#include "JShifting.h"
#include <climits>

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

  Np=Npb=NpbOk=0;
  NpbPer=NpfPer=0;

  Idpc=NULL; Codec=NULL; Dcellc=NULL; Posc=NULL; Velrhopc=NULL;
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
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_24B,1); //-SpsTau,SpsGradvel
  }
  if(Shifting){
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_16B,1); //-shiftposfs
  }
  if(InOut){  //<vs_innlet_ini>
    //ArraysCpu->AddArrayCount(JArraysCpu::SIZE_4B,1);  //-InOutPart
    ArraysCpu->AddArrayCount(JArraysCpu::SIZE_1B,1);  //-newizone
  }  //<vs_innlet_end>
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
  if(MLPistons)s+=MLPistons->GetAllocMemoryCpu();  //<vs_mlapiston>
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
void JSphCpu::ConfigOmp(const JCfgRun *cfg){
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
void JSphCpu::ConfigRunMode(const JCfgRun *cfg,std::string preinfo){
  //#ifndef WIN32  //-Error compilation when gcc5 is used.
  //  const int len=128; char hname[len];
  //  gethostname(hname,len);
  //  if(!preinfo.empty())preinfo=preinfo+", ";
  //  preinfo=preinfo+"HostName:"+hname;
  //#endif
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
    const float rhop=Velrhopc[p].w,rhop_r0=rhop/RhopZero;
    Pressc[p]=CteB*(pow(rhop_r0,Gamma)-1.0f);
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
/// Returns values of kernel Wendland, gradients: frx, fry and frz.
/// Devuelve valores de kernel Wendland, gradients: frx, fry y frz.
//==============================================================================
void JSphCpu::GetKernelWendland(float rr2,float drx,float dry,float drz
  ,float &frx,float &fry,float &frz)const
{
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  //-Wendland kernel.
  const float wqq1=1.f-0.5f*qq;
  const float fac=Bwen*qq*wqq1*wqq1*wqq1/rad; //-Kernel derivative (divided by rad).
  frx=fac*drx; fry=fac*dry; frz=fac*drz;
}

//<vs_innlet_ini>
//============================================================================== 
/// Returns values of kernel Wendland, gradients: frx, fry, frz and wab.
/// Devuelve valores de kernel Wendland, gradients: frx, fry, frz y wab.
//==============================================================================
void JSphCpu::GetKernelWendland(float rr2,float drx,float dry,float drz
  ,float &frx,float &fry,float &frz,float &wab)const
{
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  //-Wendland kernel.
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  const float fac=Bwen*qq*wqq2*wqq1/rad; //-Kernel derivative (divided by rad).
  frx=fac*drx; fry=fac*dry; frz=fac*drz;
  const float wqq=2.f*qq+1.f;
  wab=Awen*wqq*wqq2*wqq2; //-Kernel.
}  //<vs_innlet_end>

//==============================================================================
/// Returns values of kernel Gaussian, gradients: frx, fry and frz.
/// Devuelve valores de kernel Gaussian, gradients: frx, fry y frz.
//==============================================================================
void JSphCpu::GetKernelGaussian(float rr2,float drx,float dry,float drz
  ,float &frx,float &fry,float &frz)const
{
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  //-Gaussian kernel.
  const float qqexp=-4.0f*qq*qq;
  //const float wab=Agau*expf(qqexp); //-Kernel.
  const float fac=Bgau*qq*expf(qqexp)/rad; //-Kernel derivative (divided by rad).
  frx=fac*drx; fry=fac*dry; frz=fac*drz;
}

//<vs_innlet_ini>
//==============================================================================
/// Returns values of kernel Gaussian, gradients: frx, fry, frz and wab.
/// Devuelve valores de kernel Gaussian, gradients: frx, fry, frz y wab.
//==============================================================================
void JSphCpu::GetKernelGaussian(float rr2,float drx,float dry,float drz
  ,float &frx,float &fry,float &frz,float &wab)const
{
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  //-Gaussian kernel.
  const float qqexp=-4.0f*qq*qq;
  const float eqqexp=expf(qqexp);
  wab=Agau*eqqexp; //-Kernel.
  const float fac=Bgau*qq*eqqexp/rad; //-Kernel derivative (divided by rad).
  frx=fac*drx; fry=fac*dry; frz=fac*drz;
}  //<vs_innlet_end>

//==============================================================================
/// Return values of kernel Cubic without tensil correction, gradients: frx, fry and frz.
/// Devuelve valores de kernel Cubic sin correccion tensil, gradients: frx, fry y frz.
//==============================================================================
void JSphCpu::GetKernelCubic(float rr2,float drx,float dry,float drz
  ,float &frx,float &fry,float &frz)const
{
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  //-Cubic Spline kernel.
  float fac;
  if(rad>H){
    float wqq1=2.0f-qq;
    float wqq2=wqq1*wqq1;
    fac=CubicCte.c2*wqq2/rad; //-Kernel derivative (divided by rad).
  }
  else{
    float wqq2=qq*qq;
    fac=(CubicCte.c1*qq+CubicCte.d1*wqq2)/rad; //-Kernel derivative (divided by rad).
  }
  //-Gradients.
  frx=fac*drx; fry=fac*dry; frz=fac*drz;
}

//<vs_innlet_ini>
//==============================================================================
/// Return values of kernel Cubic without tensil correction, gradients: frx, fry, frz and wab.
/// Devuelve valores de kernel Cubic sin correccion tensil, gradients: frx, fry,frz y wab.
//==============================================================================
void JSphCpu::GetKernelCubic(float rr2,float drx,float dry,float drz
  ,float &frx,float &fry,float &frz,float &wab)const
{
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  //-Cubic Spline kernel.
  float fac;
  if(rad>H){
    float wqq1=2.0f-qq;
    float wqq2=wqq1*wqq1;
    fac=CubicCte.c2*wqq2/rad; //-Kernel derivative (divided by rad).
    wab=CubicCte.a24*(wqq2*wqq1); //-Kernel.
  }
  else{
    float wqq2=qq*qq;
    fac=(CubicCte.c1*qq+CubicCte.d1*wqq2)/rad; //-Kernel derivative (divided by rad).
    wab=CubicCte.a2*(1.0f+(0.75f*qq-1.5f)*wqq2); //-Kernel.
  }
  //-Gradients.
  frx=fac*drx; fry=fac*dry; frz=fac*drz;
}  //<vs_innlet_end>

//==============================================================================
/// Return tensil correction for kernel Cubic.
/// Devuelve correccion tensil para kernel Cubic.
//==============================================================================
float JSphCpu::GetKernelCubicTensil(float rr2,float rhopp1,float pressp1,float rhopp2,float pressp2)const{
  const float rad=sqrt(rr2);
  const float qq=rad/H;
  //-Cubic Spline kernel.
  float wab;
  if(rad>H){
    float wqq1=2.0f-qq;
    float wqq2=wqq1*wqq1;
    wab=CubicCte.a24*(wqq2*wqq1); //-Kernel.
  }
  else{
    float wqq2=qq*qq;
    float wqq3=wqq2*qq;
    wab=CubicCte.a2*(1.0f-1.5f*wqq2+0.75f*wqq3); //-Kernel.
  }
  //-Tensile correction.
  float fab=wab*CubicCte.od_wdeltap;
  fab*=fab; fab*=fab; //fab=fab^4
  const float tensilp1=(pressp1/(rhopp1*rhopp1))*(pressp1>0? 0.01f: -0.2f);
  const float tensilp2=(pressp2/(rhopp2*rhopp2))*(pressp2>0? 0.01f: -0.2f);
  return(fab*(tensilp1+tensilp2));
}

//==============================================================================
/// Return cell limits for interaction starting from cell coordinates.
/// Devuelve limites de celdas para interaccion a partir de coordenadas de celda.
//==============================================================================
void JSphCpu::GetInteractionCells(unsigned rcell
  ,int hdiv,const tint4 &nc,const tint3 &cellzero
  ,int &cxini,int &cxfin,int &yini,int &yfin,int &zini,int &zfin)const
{
  //-Get interaction limits. | Obtiene limites de interaccion.
  const int cx=PC__Cellx(DomCellCode,rcell)-cellzero.x;
  const int cy=PC__Celly(DomCellCode,rcell)-cellzero.y;
  const int cz=PC__Cellz(DomCellCode,rcell)-cellzero.z;
  //-Code for hdiv 1 or 2 but not zero. | Codigo para hdiv 1 o 2 pero no cero.
  cxini=cx-min(cx,hdiv);
  cxfin=cx+min(nc.x-cx-1,hdiv)+1;
  yini=cy-min(cy,hdiv);
  yfin=cy+min(nc.y-cy-1,hdiv)+1;
  zini=cz-min(cz,hdiv);
  zfin=cz+min(nc.z-cz-1,hdiv)+1;
}

//============================================================================== //<vs_innlet_ini>
/// Return cell limits for interaction starting from position.
/// Devuelve limites de celdas para interaccion a partir de posicion.
//==============================================================================
void JSphCpu::GetInteractionCells(const tdouble3 &pos
  ,int hdiv,const tint4 &nc,const tint3 &cellzero
  ,int &cxini,int &cxfin,int &yini,int &yfin,int &zini,int &zfin)const
{
  //-Get cell coordinates of position pos.
  const int cx=int((pos.x-DomPosMin.x)/Scell)-cellzero.x;
  const int cy=int((pos.y-DomPosMin.y)/Scell)-cellzero.y;
  const int cz=int((pos.z-DomPosMin.z)/Scell)-cellzero.z;
  //-code for hdiv 1 or 2 but not zero. | Codigo para hdiv 1 o 2 pero no cero.
  cxini=cx-min(cx,hdiv);
  cxfin=cx+min(nc.x-cx-1,hdiv)+1;
  yini=cy-min(cy,hdiv);
  yfin=cy+min(nc.y-cy-1,hdiv)+1;
  zini=cz-min(cz,hdiv);
  zfin=cz+min(nc.z-cz-1,hdiv)+1;
} //<vs_innlet_end>

//==============================================================================
/// Perform interaction between particles. Bound-Fluid/Float
/// Realiza interaccion entre particulas. Bound-Fluid/Float
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> void JSphCpu::InteractionForcesBound
  (unsigned n,unsigned pinit,tint4 nc,int hdiv,unsigned cellinitial
  ,const unsigned *beginendcell,tint3 cellzero,const unsigned *dcell
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
    const bool rsymp1=(Symmetry && posp1.y<=Dosh); //<vs_syymmetry>
    const tfloat3 velp1=TFloat3(velrhop[p1].x,velrhop[p1].y,velrhop[p1].z);

    //-Obtain limits of interaction. | Obtiene limites de interaccion.
    int cxini,cxfin,yini,yfin,zini,zfin;
    GetInteractionCells(dcell[p1],hdiv,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin);

    //-Search for neighbours in adjacent cells. | Busqueda de vecinos en celdas adyacentes.
    for(int z=zini;z<zfin;z++){
      const int zmod=(nc.w)*z+cellinitial; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        const unsigned pini=beginendcell[cxini+ymod];
        const unsigned pfin=beginendcell[cxfin+ymod];

        //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
        //---------------------------------------------------------------------------------------------
        bool rsym=false; //<vs_syymmetry>
        for(unsigned p2=pini;p2<pfin;p2++){
          const float drx=float(posp1.x-pos[p2].x);
                float dry=float(posp1.y-pos[p2].y);
          if(rsym)    dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=Fourh2 && rr2>=ALMOSTZERO){
            //-Cubic Spline, Wendland or Gaussian kernel.
            float frx,fry,frz;
            if(tker==KERNEL_Wendland)GetKernelWendland(rr2,drx,dry,drz,frx,fry,frz);
            else if(tker==KERNEL_Gaussian)GetKernelGaussian(rr2,drx,dry,drz,frx,fry,frz);
            else if(tker==KERNEL_Cubic)GetKernelCubic(rr2,drx,dry,drz,frx,fry,frz);

            //===== Get mass of particle p2 ===== 
            float massp2=MassFluid; //-Contains particle mass of incorrect fluid. | Contiene masa de particula por defecto fluid.
            bool compute=true;      //-Deactivate when using DEM and/or bound-float. | Se desactiva cuando se usa DEM y es bound-float.
            if(USE_FLOATING){
              bool ftp2=CODE_IsFloating(code[p2]);
              if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
              compute=!(USE_FTEXTERNAL && ftp2); //-Deactivate when using DEM/Chrono and/or bound-float. | Se desactiva cuando se usa DEM/Chrono y es bound-float.
            }

            if(compute){
              //-Density derivative.
              //const float dvx=velp1.x-velrhop[p2].x, dvy=velp1.y-velrhop[p2].y, dvz=velp1.z-velrhop[p2].z;
              tfloat4 velrhop2=velrhop[p2];
              if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>
              const float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
              if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz);

              {//-Viscosity.
                const float dot=drx*dvx + dry*dvy + drz*dvz;
                const float dot_rr2=dot/(rr2+Eta2);
                visc=max(dot_rr2,visc);
              }
            }
            rsym=(rsymp1 && !rsym && float(posp1.y-dry)<=Dosh); //<vs_syymmetry>
            if(rsym)p2--;                                       //<vs_syymmetry>
          }
          else rsym=false;                                      //<vs_syymmetry>
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
/// Perform interaction between particles: Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Realiza interaccion entre particulas: Fluid/Float-Fluid/Float or Fluid/Float-Bound
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,bool lamsps,TpDensity tdensity,bool shift> 
  void JSphCpu::InteractionForcesFluid
  (unsigned n,unsigned pinit,tint4 nc,int hdiv,unsigned cellinitial,float visco
  ,const unsigned *beginendcell,tint3 cellzero,const unsigned *dcell
  ,const tsymatrix3f* tau,tsymatrix3f* gradvel
  ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
  ,const float *press 
  ,float &viscdt,float *ar,tfloat3 *ace,float *delta
  ,TpShifting shiftmode,tfloat4 *shiftposfs)const
{
  const bool boundp2=(!cellinitial); //-Interaction with type boundary (Bound). | Interaccion con Bound.
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
    const tsymatrix3f taup1=(lamsps? tau[p1]: gradvelp1);
    const bool rsymp1=(Symmetry && posp1.y<=Dosh); //<vs_syymmetry>

    //-Obtain interaction limits.
    int cxini,cxfin,yini,yfin,zini,zfin;
    GetInteractionCells(dcell[p1],hdiv,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin);

    //-Search for neighbours in adjacent cells.
    for(int z=zini;z<zfin;z++){
      const int zmod=(nc.w)*z+cellinitial; //-Sum from start of fluid or boundary cells. | Le suma donde empiezan las celdas de fluido o bound.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        const unsigned pini=beginendcell[cxini+ymod];
        const unsigned pfin=beginendcell[cxfin+ymod];

        //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
        //------------------------------------------------------------------------------------------------
        bool rsym=false; //<vs_syymmetry>
        for(unsigned p2=pini;p2<pfin;p2++){
          const float drx=float(posp1.x-pos[p2].x);
                float dry=float(posp1.y-pos[p2].y);
          if(rsym)    dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=Fourh2 && rr2>=ALMOSTZERO){
            //-Cubic Spline, Wendland or Gaussian kernel.
            float frx,fry,frz;
            if(tker==KERNEL_Wendland)GetKernelWendland(rr2,drx,dry,drz,frx,fry,frz);
            else if(tker==KERNEL_Gaussian)GetKernelGaussian(rr2,drx,dry,drz,frx,fry,frz);
            else if(tker==KERNEL_Cubic)GetKernelCubic(rr2,drx,dry,drz,frx,fry,frz);

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
            //===== Acceleration ===== 
            if(compute){
              const float prs=(pressp1+press[p2])/(rhopp1*velrhop2.w) + (tker==KERNEL_Cubic? GetKernelCubicTensil(rr2,rhopp1,pressp1,velrhop2.w,press[p2]): 0);
              const float p_vpm=-prs*massp2;
              acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
            }

            //-Density derivative.
            const float dvx=velp1.x-velrhop2.x, dvy=velp1.y-velrhop2.y, dvz=velp1.z-velrhop2.z;
            if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz);

            const float cbar=(float)Cs0;
            //-Density Diffusion Term (Molteni and Colagrossi 2009).
            if(tdensity==DDT_DDT && deltap1!=FLT_MAX){
              const float rhop1over2=rhopp1/velrhop2.w;
              const float visc_densi=DDT2h*cbar*(rhop1over2-1.f)/(rr2+Eta2);
              const float dot3=(drx*frx+dry*fry+drz*frz);
              const float delta=visc_densi*dot3*massp2;
              //deltap1=(boundp2? FLT_MAX: deltap1+delta);
              deltap1=(boundp2 && TBoundary==BC_DBC? FLT_MAX: deltap1+delta);
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
              if(!lamsps){//-Artificial viscosity.
                if(dot<0){
                  const float amubar=H*dot_rr2;  //amubar=CTE.h*dot/(rr2+CTE.eta2);
                  const float robar=(rhopp1+velrhop2.w)*0.5f;
                  const float pi_visc=(-visco*cbar*amubar/robar)*massp2;
                  acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
                }
              }
              else{//-Laminar+SPS viscosity. 
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
            rsym=(rsymp1 && !rsym && float(posp1.y-dry)<=Dosh); //<vs_syymmetry>
            if(rsym)p2--;                                       //<vs_syymmetry>
          }
          else rsym=false;                                      //<vs_syymmetry>
        }
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
      if(lamsps){
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
void JSphCpu::InteractionForcesDEM
  (unsigned nfloat,tint4 nc,int hdiv,unsigned cellfluid
  ,const unsigned *beginendcell,tint3 cellzero,const unsigned *dcell
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

      //-Get interaction limits.
      int cxini,cxfin,yini,yfin,zini,zfin;
      GetInteractionCells(dcell[p1],hdiv,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin);

      //-Search for neighbours in adjacent cells (first bound and then fluid+floating).
      for(unsigned cellinitial=0;cellinitial<=cellfluid;cellinitial+=cellfluid){
        for(int z=zini;z<zfin;z++){
          const int zmod=(nc.w)*z+cellinitial; //-Sum from start of fluid or boundary cells. | Le suma donde empiezan las celdas de fluido o bound.
          for(int y=yini;y<yfin;y++){
            int ymod=zmod+nc.x*y;
            const unsigned pini=beginendcell[cxini+ymod];
            const unsigned pfin=beginendcell[cxfin+ymod];

            //-Interaction of Floating Object particles with type Fluid or Bound. | Interaccion de Floating con varias Fluid o Bound.
            //-----------------------------------------------------------------------------------------------------------------------
            for(unsigned p2=pini;p2<pfin;p2++)if(CODE_IsNotFluid(code[p2]) && tavp1!=CODE_GetTypeAndValue(code[p2])){
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

              const float nu_mass=(!cellinitial? masstotp1/2: masstotp1*masstotp2/(masstotp1+masstotp2)); //-Con boundary toma la propia masa del floating 1.
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
template<TpKernel tker,TpFtMode ftmode,bool lamsps,TpDensity tdensity,bool shift>
  void JSphCpu::Interaction_ForcesCpuT(const stinterparmsc &t,float &viscdt)const
{
  const tint4 nc=TInt4(int(t.ncells.x),int(t.ncells.y),int(t.ncells.z),int(t.ncells.x*t.ncells.y));
  const tint3 cellzero=TInt3(t.cellmin.x,t.cellmin.y,t.cellmin.z);
  const unsigned cellfluid=nc.w*nc.z+1;
  const int hdiv=(CellMode==CELLMODE_H? 2: 1);
  
  if(t.npf){
    //-Interaction Fluid-Fluid.
    InteractionForcesFluid<tker,ftmode,lamsps,tdensity,shift> (t.npf,t.npb,nc,hdiv,cellfluid,Visco                 ,t.begincell,cellzero,t.dcell,t.spstau,t.spsgradvel,t.pos,t.velrhop,t.code,t.idp,t.press,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);
    //-Interaction Fluid-Bound.
    InteractionForcesFluid<tker,ftmode,lamsps,tdensity,shift> (t.npf,t.npb,nc,hdiv,0        ,Visco*ViscoBoundFactor,t.begincell,cellzero,t.dcell,t.spstau,t.spsgradvel,t.pos,t.velrhop,t.code,t.idp,t.press,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);

    //-Interaction of DEM Floating-Bound & Floating-Floating. //(DEM)
    if(UseDEM)InteractionForcesDEM(CaseNfloat,nc,hdiv,cellfluid,t.begincell,cellzero,t.dcell,FtRidp,DemData,t.pos,t.velrhop,t.code,t.idp,viscdt,t.ace);

    //-Computes tau for Laminar+SPS.
    if(lamsps)ComputeSpsTau(t.npf,t.npb,t.velrhop,t.spsgradvel,t.spstau);
  }
  if(t.npbok){
    //-Interaction Bound-Fluid.
    InteractionForcesBound<tker,ftmode> (t.npbok,0,nc,hdiv,cellfluid,t.begincell,cellzero,t.dcell,t.pos,t.velrhop,t.code,t.idp,viscdt,t.ar);
  }
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,bool lamsps,TpDensity tdensity> void JSphCpu::Interaction_Forces_ct5(const stinterparmsc &t,float &viscdt)const{
  if(Shifting)Interaction_ForcesCpuT<tker,ftmode,lamsps,tdensity,true >(t,viscdt);
  else        Interaction_ForcesCpuT<tker,ftmode,lamsps,tdensity,false>(t,viscdt);
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,bool lamsps> void JSphCpu::Interaction_Forces_ct4(const stinterparmsc &t,float &viscdt)const{
       if(TDensity==DDT_None)    Interaction_Forces_ct5<tker,ftmode,lamsps,DDT_None    >(t,viscdt);
  else if(TDensity==DDT_DDT)     Interaction_Forces_ct5<tker,ftmode,lamsps,DDT_DDT     >(t,viscdt);
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> void JSphCpu::Interaction_Forces_ct3(const stinterparmsc &t,float &viscdt)const{
  if(TVisco==VISCO_LaminarSPS)Interaction_Forces_ct4<tker,ftmode,true >(t,viscdt);
  else                        Interaction_Forces_ct4<tker,ftmode,false>(t,viscdt);
}
//==============================================================================
template<TpKernel tker> void JSphCpu::Interaction_Forces_ct2(const stinterparmsc &t,float &viscdt)const{
       if(FtMode==FTMODE_None)Interaction_Forces_ct3<tker,FTMODE_None>(t,viscdt);
  else if(FtMode==FTMODE_Sph )Interaction_Forces_ct3<tker,FTMODE_Sph >(t,viscdt);
  else if(FtMode==FTMODE_Ext )Interaction_Forces_ct3<tker,FTMODE_Ext >(t,viscdt);
}
//==============================================================================
void JSphCpu::Interaction_Forces_ct(const stinterparmsc &t,float &viscdt)const{
       if(TKernel==KERNEL_Wendland)Interaction_Forces_ct2<KERNEL_Wendland>(t,viscdt);
  else if(TKernel==KERNEL_Gaussian)Interaction_Forces_ct2<KERNEL_Gaussian>(t,viscdt);
  else if(TKernel==KERNEL_Cubic)   Interaction_Forces_ct2<KERNEL_Cubic   >(t,viscdt);
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
template<bool shift> void JSphCpu::ComputeVerletVarsFluid(
  const tfloat4 *velrhop1,const tfloat4 *velrhop2,double dt,double dt2
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
      //-Calculate displacement and update position. | Calcula desplazamiento y actualiza posicion.
      double dx=double(velrhop1[p].x)*dt + acegr.x*dt205;
      double dy=double(velrhop1[p].y)*dt + acegr.y*dt205;
      double dz=double(velrhop1[p].z)*dt + acegr.z*dt205;
      if(shift){
        dx+=double(ShiftPosfsc[p].x);
        dy+=double(ShiftPosfsc[p].y);
        dz+=double(ShiftPosfsc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin||rhopnew>RhopOutMax);
      UpdatePos(pos[p],dx,dy,dz,outrhop,p,pos,dcell,code);
      //-Update velocity & density. | Actualiza velocidad y densidad.
      velrhopnew[p].x=float(double(velrhop2[p].x) + acegr.x*dt2);
      velrhopnew[p].y=float(double(velrhop2[p].y) + acegr.y*dt2);
      velrhopnew[p].z=float(double(velrhop2[p].z) + acegr.z*dt2);
      velrhopnew[p].w=rhopnew;
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
  VerletStep++;
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    if(Shifting)ComputeVerletVarsFluid<true>  (Velrhopc,VelrhopM1c,dt,twodt,Posc,Dcellc,Codec,VelrhopM1c);
    else        ComputeVerletVarsFluid<false> (Velrhopc,VelrhopM1c,dt,twodt,Posc,Dcellc,Codec,VelrhopM1c);
    ComputeVelrhopBound(VelrhopM1c,twodt,VelrhopM1c);
  }
  else{
    if(Shifting)ComputeVerletVarsFluid<true>  (Velrhopc,Velrhopc,dt,dt,Posc,Dcellc,Codec,VelrhopM1c);
    else        ComputeVerletVarsFluid<false> (Velrhopc,Velrhopc,dt,dt,Posc,Dcellc,Codec,VelrhopM1c);
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
  if(Shifting)ComputeSymplecticPreT<false>(dt); //-We strongly recommend running the shifting correction only for the corrector. If you want to re-enable shifting in the predictor, change the value here to "true".
  else        ComputeSymplecticPreT<false>(dt);
}

//==============================================================================
/// Update of particles according to forces and dt using Symplectic-Predictor.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
template<bool shift> void JSphCpu::ComputeSymplecticPreT(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  //-Assign memory to variables Pre. | Asigna memoria a variables Pre.
  PosPrec=ArraysCpu->ReserveDouble3();
  VelrhopPrec=ArraysCpu->ReserveFloat4();
  //-Change data to variables Pre to calculate new data. | Cambia datos a variables Pre para calcular nuevos datos.
  swap(PosPrec,Posc);         //Put value of Pos[] in PosPre[].         | Es decir... PosPre[] <= Pos[].
  swap(VelrhopPrec,Velrhopc); //Put value of Velrhop[] in VelrhopPre[]. | Es decir... VelrhopPre[] <= Velrhop[].
  //-Calculate new values of particles. | Calcula nuevos datos de particulas.
  const double dt05=dt*.5;
  
  //-Calculate new density for boundary and copy velocity. | Calcula nueva densidad para el contorno y copia velocidad.
  const int npb=int(Npb);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<npb;p++){
    const tfloat4 vr=VelrhopPrec[p];
    const float rhopnew=float(double(vr.w)+dt05*Arc[p]);
    Velrhopc[p]=TFloat4(vr.x,vr.y,vr.z,(rhopnew<RhopZero? RhopZero: rhopnew));//-Avoid fluid particles being absorbed by boundary ones. | Evita q las boundary absorvan a las fluidas.
  }

  //-Calculate new values of fluid. | Calcula nuevos datos del fluido.
  const int np=int(Np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(np>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=npb;p<np;p++){
    //-Calculate density.
    const float rhopnew=float(double(VelrhopPrec[p].w)+dt05*Arc[p]);
    if(!WithFloating || CODE_IsFluid(Codec[p])){//-Fluid Particles.
      //-Calculate displacement & update position. | Calcula desplazamiento y actualiza posicion.
      double dx=double(VelrhopPrec[p].x)*dt05;
      double dy=double(VelrhopPrec[p].y)*dt05;
      double dz=double(VelrhopPrec[p].z)*dt05;
      if(shift){
        dx+=double(ShiftPosfsc[p].x);
        dy+=double(ShiftPosfsc[p].y);
        dz+=double(ShiftPosfsc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin||rhopnew>RhopOutMax);
      UpdatePos(PosPrec[p],dx,dy,dz,outrhop,p,Posc,Dcellc,Codec);
      //-Update velocity & density. | Actualiza velocidad y densidad.
      Velrhopc[p].x=float(double(VelrhopPrec[p].x) + (double(Acec[p].x)+Gravity.x) * dt05);
      Velrhopc[p].y=float(double(VelrhopPrec[p].y) + (double(Acec[p].y)+Gravity.y) * dt05);
      Velrhopc[p].z=float(double(VelrhopPrec[p].z) + (double(Acec[p].z)+Gravity.z) * dt05);
      Velrhopc[p].w=rhopnew;
    }
    else{//-Floating Particles.
      Velrhopc[p]=VelrhopPrec[p];
      Velrhopc[p].w=(rhopnew<RhopZero? RhopZero: rhopnew); //-Avoid fluid particles being absorbed by floating ones. | Evita q las floating absorvan a las fluidas.
      //-Copy position. | Copia posicion.
      Posc[p]=PosPrec[p];
    }
  }

  //-Copy previous position of boundary. | Copia posicion anterior del contorno.
  memcpy(Posc,PosPrec,sizeof(tdouble3)*Npb);

  TmcStop(Timers,TMC_SuComputeStep);
}

//==============================================================================
/// Update particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphCpu::ComputeSymplecticCorr(double dt){
  if(Shifting)ComputeSymplecticCorrT<true> (dt);
  else        ComputeSymplecticCorrT<false>(dt);
}

//==============================================================================
/// Update particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
template<bool shift> void JSphCpu::ComputeSymplecticCorrT(double dt){
  TmcStart(Timers,TMC_SuComputeStep);
  
  //-Calculate rhop of boudary and set velocity=0. | Calcula rhop de contorno y vel igual a cero.
  const int npb=int(Npb);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npb>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=0;p<npb;p++){
    const double epsilon_rdot=(-double(Arc[p])/double(Velrhopc[p].w))*dt;
    const float rhopnew=float(double(VelrhopPrec[p].w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
    Velrhopc[p]=TFloat4(0,0,0,(rhopnew<RhopZero? RhopZero: rhopnew));//-Avoid fluid particles being absorbed by boundary ones. | Evita q las boundary absorvan a las fluidas.
  }

  //-Calculate fluid values. | Calcula datos de fluido.
  const double dt05=dt*.5;
  const int np=int(Np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(np>OMP_LIMIT_COMPUTESTEP)
  #endif
  for(int p=npb;p<np;p++){
    const double epsilon_rdot=(-double(Arc[p])/double(Velrhopc[p].w))*dt;
    const float rhopnew=float(double(VelrhopPrec[p].w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
    if(!WithFloating || CODE_IsFluid(Codec[p])){//-Fluid Particles.
      //-Update velocity & density. | Actualiza velocidad y densidad.
      Velrhopc[p].x=float(double(VelrhopPrec[p].x) + (double(Acec[p].x)+Gravity.x) * dt); 
      Velrhopc[p].y=float(double(VelrhopPrec[p].y) + (double(Acec[p].y)+Gravity.y) * dt); 
      Velrhopc[p].z=float(double(VelrhopPrec[p].z) + (double(Acec[p].z)+Gravity.z) * dt); 
      Velrhopc[p].w=rhopnew;
      //-Calculate displacement and update position. | Calcula desplazamiento y actualiza posicion.
      double dx=(double(VelrhopPrec[p].x)+double(Velrhopc[p].x)) * dt05; 
      double dy=(double(VelrhopPrec[p].y)+double(Velrhopc[p].y)) * dt05; 
      double dz=(double(VelrhopPrec[p].z)+double(Velrhopc[p].z)) * dt05;
      if(shift){
        dx+=double(ShiftPosfsc[p].x);
        dy+=double(ShiftPosfsc[p].y);
        dz+=double(ShiftPosfsc[p].z);
      }
      bool outrhop=(rhopnew<RhopOutMin||rhopnew>RhopOutMax);
      UpdatePos(PosPrec[p],dx,dy,dz,outrhop,p,Posc,Dcellc,Codec);
    }
    else{//-Floating Particles.
      Velrhopc[p]=VelrhopPrec[p];
      Velrhopc[p].w=(rhopnew<RhopZero? RhopZero: rhopnew); //-Avoid fluid particles being absorbed by floating ones. | Evita q las floating absorvan a las fluidas.
      //-Copy position. | Copia posicion.
      Posc[p]=PosPrec[p];
    }
  }

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
  const double dt1=(AceMax? (sqrt(double(H)/AceMax)): DBL_MAX); 
  //-dt2 combines the Courant and the viscous time-step controls.
  const double dt2=double(H)/(max(Cs0,VelMax*10.)+double(H)*ViscDtMax);
  //-dt new value of time step.
  double dt=double(CFLnumber)*min(dt1,dt2);
  if(DtFixed)dt=DtFixed->GetDt(TimeStep,dt);
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
  const bool motsim=true;
  BoundChanged=false;
  //-Add motion from automatic wave generation.
  if(WaveGen)CalcMotionWaveGen(stepdt);
  //-Process particles motion.
  if(SphMotion->GetActiveMotion()){
    CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codec,Idpc,RidpMove);
    BoundChanged=true;
    const unsigned nref=SphMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=SphMotion->GetMotionData(ref);
      if(m.type==MOTT_Linear){//-Linear movement.
        if(motsim)MoveLinBound   (m.count,m.idbegin-CaseNfixed,m.linmov,ToTFloat3(m.linvel),RidpMove,Posc,Dcellc,Velrhopc,Codec);
        //else    MoveLinBoundAce(m.count,m.idbegin-CaseNfixed,m.linmov,ToTFloat3(m.linvel),ToTFloat3(m.linace),RidpMove,Posc,Dcellc,Velrhopc,Acec,Codec);
      }
      if(m.type==MOTT_Matrix){//-Matrix movement (for rotations).
        if(motsim)MoveMatBound   (m.count,m.idbegin-CaseNfixed,m.matmov,stepdt,RidpMove,Posc,Dcellc,Velrhopc,Codec,boundnormal); 
        //else    MoveMatBoundAce(m.count,m.idbegin-CaseNfixed,m.matmov,m.matmov2,stepdt,RidpMove,Posc,Dcellc,Velrhopc,Acec,Codec);
      }      
      //-Applies predefined motion to BoundCorr configuration.           //<vs_innlet> 
      if(BoundCorr && BoundCorr->GetUseMotion())BoundCorr->RunMotion(m); //<vs_innlet> 
    }
  }
  //-Management of Multi-Layer Pistons.  //<vs_mlapiston_ini>
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
  }  //<vs_mlapiston_end>
  TmcStop(Timers,TMC_SuMotion);
}

//<vs_mlapiston_ini>
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
//<vs_mlapiston_end>

//<vs_rzone_ini>
//==============================================================================
/// Applies RelaxZone to selected particles.
/// Aplica RelaxZone a las particulas indicadas.
//==============================================================================
void JSphCpu::RunRelaxZone(double dt){
  TmcStart(Timers,TMC_SuMotion);
  RelaxZones->SetFluidVel(TimeStep,dt,Np-Npb,Npb,Posc,Idpc,Velrhopc);
  TmcStop(Timers,TMC_SuMotion);
}
//<vs_rzone_end>

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
    for(unsigned cf=0;cf<FtCount;cf++)ftdata.CheckHeadData(cf,FtObjs[cf].mkbound,FtObjs[cf].begin,FtObjs[cf].count,FtObjs[cf].mass);
    //-Load PART data. | Carga datos de PART.
    ftdata.LoadPart(PartBegin);
    for(unsigned cf=0;cf<FtCount;cf++){
      FtObjs[cf].center=ftdata.GetPartCenter(cf);
      FtObjs[cf].fvel=ftdata.GetPartFvel(cf);
      FtObjs[cf].fomega=ftdata.GetPartFomega(cf);
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


