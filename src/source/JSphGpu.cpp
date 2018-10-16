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

/// \file JSphGpu.cpp \brief Implements the class \ref JSphGpu.

#include "JSphGpu.h"
#include "JSphGpu_ker.h"
#include "JBlockSizeAuto.h"
#include "JCellDivGpu.h"
#include "JPartFloatBi4.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include "JSphMotion.h"
#include "JArraysGpu.h"
#include "JSphDtFixed.h"
#include "JSaveDt.h"
#include "JTimeOut.h"
#include "JWaveGen.h"
#include "JDamping.h"
#include "JSphAccInput.h"
#include "JXml.h"
#include "JFormatFiles2.h"
#include "JGaugeSystem.h"
#include <climits>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JSphGpu::JSphGpu(bool withmpi):JSph(false,withmpi){
  ClassName="JSphGpu";
  Idp=NULL; Code=NULL; Dcell=NULL; Posxy=NULL; Posz=NULL; Velrhop=NULL;
  AuxPos=NULL; AuxVel=NULL; AuxRhop=NULL;
  CellDiv=NULL;
  FtoAuxDouble6=NULL; FtoAuxFloat9=NULL; //-Calculates forces on floating bodies.
  ArraysGpu=new JArraysGpu;
  InitVars();
  TmgCreation(Timers,false);
  BsAuto=NULL;
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphGpu::~JSphGpu(){
  DestructorActive=true;
  FreeCpuMemoryParticles();
  FreeGpuMemoryParticles();
  FreeGpuMemoryFixed();
  delete ArraysGpu;
  TmgDestruction(Timers);
  cudaDeviceReset();
  delete BsAuto; BsAuto=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphGpu::InitVars(){
  RunMode="";
  memset(&BlockSizes,0,sizeof(StBlockSizes));
  BlockSizesStr="";
  BlockSizeMode=BSIZEMODE_Fixed;

  Np=Npb=NpbOk=0;
  NpbPer=NpfPer=0;

  FreeCpuMemoryParticles();
  FreeCpuMemoryFixed();
  Idpg=NULL; Codeg=NULL; Dcellg=NULL; Posxyg=NULL; Poszg=NULL; Velrhopg=NULL;
  VelrhopM1g=NULL;                                 //-Verlet
  PosxyPreg=NULL; PoszPreg=NULL; VelrhopPreg=NULL; //-Symplectic
  PsPospressg=NULL;                                //-Interaction Pos-Single.
  SpsTaug=NULL; SpsGradvelg=NULL;                  //-Laminar+SPS. 
  ViscDtg=NULL; 
  Arg=NULL; Aceg=NULL; Deltag=NULL;
  ShiftPosg=NULL; ShiftDetectg=NULL; //-Shifting.
  RidpMoveg=NULL;
  FtRidpg=NULL;   FtoMasspg=NULL;                      //-Floatings.
  FtoDatag=NULL;  FtoForcesSumg=NULL;  FtoForcesg=NULL;  FtoForcesResg=NULL;  FtoCenterResg=NULL; //-Calculates forces on floating bodies.
  FtoCenterg=NULL; FtoAnglesg=NULL; FtoVelg=NULL; FtoOmegag=NULL;//-Management of floating bodies.
  FtoInertiaini8g=NULL; FtoInertiaini1g=NULL;//-Management of floating bodies.
  DemDatag=NULL; //(DEM)
  FreeGpuMemoryParticles();
  FreeGpuMemoryFixed();
}

//==============================================================================
/// Throws exception for an error in the CUDA code.
/// Lanza excepcion por un error Cuda.
//==============================================================================
void JSphGpu::RunExceptionCuda(const std::string &method,const std::string &msg,cudaError_t error){
  std::string tx=fun::PrintStr("%s (CUDA error: %s).\n",msg.c_str(),cudaGetErrorString(error)); 
  Log->Print(GetExceptionText(method,tx));
  RunException(method,msg);
}

//==============================================================================
/// Checks error and throws exception.
/// Comprueba error y lanza excepcion si lo hubiera.
//==============================================================================
void JSphGpu::CheckCudaError(const std::string &method,const std::string &msg){
  cudaError_t err=cudaGetLastError();
  if(err!=cudaSuccess)RunExceptionCuda(method,msg,err);
}

//==============================================================================
/// Frees fixed memory on CPU for moving and floating bodies.
/// Libera memoria fija en CPU para moving y floating.
//==============================================================================
void JSphGpu::FreeCpuMemoryFixed(){
  MemCpuFixed=0;
  delete[] FtoAuxDouble6; FtoAuxDouble6=NULL;
  delete[] FtoAuxFloat9;  FtoAuxFloat9=NULL;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphGpu::AllocCpuMemoryFixed(){
  MemCpuFixed=0;
  if(sizeof(tfloat3)*2!=sizeof(StFtoForces))RunException("AllocCpuMemoryFixed","Error: FtoForcesg and FtoForcesSumg does not match float3*2.");
  if(sizeof(float)*9!=sizeof(tmatrix3f))RunException("AllocCpuMemoryFixed","Error: FtoInertiainig does not match float*9.");
  try{
    //-Allocates memory for floating bodies.
    if(CaseNfloat){
      FtoAuxDouble6=new tdouble3[FtCount*2];  MemCpuFixed+=(sizeof(tdouble3)*FtCount*2);
      FtoAuxFloat9 =new tfloat3 [FtCount*3];  MemCpuFixed+=(sizeof(tfloat3) *FtCount*3);
    }
  }
  catch(const std::bad_alloc){
    RunException("AllocCpuMemoryFixed","Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Frees fixed memory on the GPU for moving and floating bodies.
/// Libera memoria fija en GPU para moving y floating.
//==============================================================================
void JSphGpu::FreeGpuMemoryFixed(){
  MemGpuFixed=0;
  if(RidpMoveg)      cudaFree(RidpMoveg);       RidpMoveg=NULL;
  if(FtRidpg)        cudaFree(FtRidpg);         FtRidpg=NULL;
  if(FtoMasspg)      cudaFree(FtoMasspg);       FtoMasspg=NULL;
  if(FtoDatag)       cudaFree(FtoDatag);        FtoDatag=NULL;
  if(FtoForcesSumg)  cudaFree(FtoForcesSumg);   FtoForcesSumg=NULL;
  if(FtoForcesg)     cudaFree(FtoForcesg);      FtoForcesg=NULL;
  if(FtoForcesResg)  cudaFree(FtoForcesResg);   FtoForcesResg=NULL;
  if(FtoCenterg)     cudaFree(FtoCenterg);      FtoCenterg=NULL;
  if(FtoAnglesg)     cudaFree(FtoAnglesg);      FtoAnglesg=NULL;
  if(FtoVelg)        cudaFree(FtoVelg);         FtoVelg=NULL;
  if(FtoOmegag)      cudaFree(FtoOmegag);       FtoOmegag=NULL;
  if(FtoInertiaini8g)cudaFree(FtoInertiaini8g); FtoInertiaini8g=NULL;
  if(FtoInertiaini1g)cudaFree(FtoInertiaini1g); FtoInertiaini1g=NULL;
  if(DemDatag)       cudaFree(DemDatag);        DemDatag=NULL;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphGpu::AllocGpuMemoryFixed(){
  MemGpuFixed=0;
  //-Allocates memory for moving objects.
  if(CaseNmoving){
    size_t m=sizeof(unsigned)*CaseNmoving;
    cudaMalloc((void**)&RidpMoveg,m);   MemGpuFixed+=m;
  }
  //-Allocates memory for floating bodies.
  if(CaseNfloat){
    size_t m=0;
    m=sizeof(unsigned)*CaseNfloat;  cudaMalloc((void**)&FtRidpg        ,m);  MemGpuFixed+=m;
    m=sizeof(float)   *FtCount;     cudaMalloc((void**)&FtoMasspg      ,m);  MemGpuFixed+=m;
    m=sizeof(float4)  *FtCount;     cudaMalloc((void**)&FtoDatag       ,m);  MemGpuFixed+=m;
    m=sizeof(float3)  *FtCount*2;   cudaMalloc((void**)&FtoForcesSumg  ,m);  MemGpuFixed+=m;
    m=sizeof(float3)  *FtCount*2;   cudaMalloc((void**)&FtoForcesg     ,m);  MemGpuFixed+=m;
    m=sizeof(float3)  *FtCount*2;   cudaMalloc((void**)&FtoForcesResg  ,m);  MemGpuFixed+=m;
    m=sizeof(double3) *FtCount;     cudaMalloc((void**)&FtoCenterResg  ,m);  MemGpuFixed+=m;
    m=sizeof(double3) *FtCount;     cudaMalloc((void**)&FtoCenterg     ,m);  MemGpuFixed+=m;
    m=sizeof(float3)  *FtCount;     cudaMalloc((void**)&FtoAnglesg     ,m);  MemGpuFixed+=m;
    m=sizeof(float3)  *FtCount;     cudaMalloc((void**)&FtoVelg        ,m);  MemGpuFixed+=m;
    m=sizeof(float3)  *FtCount;     cudaMalloc((void**)&FtoOmegag      ,m);  MemGpuFixed+=m;
    m=sizeof(float4)  *FtCount*2;   cudaMalloc((void**)&FtoInertiaini8g,m);  MemGpuFixed+=m;
    m=sizeof(float)   *FtCount;     cudaMalloc((void**)&FtoInertiaini1g,m);  MemGpuFixed+=m;
  }
  if(UseDEM){ //(DEM)
    size_t m=sizeof(float4)*DemDataSize;
    cudaMalloc((void**)&DemDatag,m);     MemGpuFixed+=m;
  }
}

//==============================================================================
/// Releases memory for the main particle data.
/// Libera memoria para datos principales de particulas.
//==============================================================================
void JSphGpu::FreeCpuMemoryParticles(){
  CpuParticlesSize=0;
  MemCpuParticles=0;
  delete[] Idp;        Idp=NULL;
  delete[] Code;       Code=NULL;
  delete[] Dcell;      Dcell=NULL;
  delete[] Posxy;      Posxy=NULL;
  delete[] Posz;       Posz=NULL;
  delete[] Velrhop;    Velrhop=NULL;
  delete[] AuxPos;     AuxPos=NULL;
  delete[] AuxVel;     AuxVel=NULL;
  delete[] AuxRhop;    AuxRhop=NULL;
}

//==============================================================================
/// Allocates memory for the main particle data.
/// Reserva memoria para datos principales de particulas.
//==============================================================================
void JSphGpu::AllocCpuMemoryParticles(unsigned np){
  const char* met="AllocCpuMemoryParticles";
  FreeCpuMemoryParticles();
  if(np>0)np=np+PARTICLES_OVERMEMORY_MIN;
  CpuParticlesSize=np;
  if(np>0){
    try{
      Idp=new unsigned[np];      MemCpuParticles+=sizeof(unsigned)*np;
      Code=new typecode[np];     MemCpuParticles+=sizeof(typecode)*np;
      Dcell=new unsigned[np];    MemCpuParticles+=sizeof(unsigned)*np;
      Posxy=new tdouble2[np];    MemCpuParticles+=sizeof(tdouble2)*np;
      Posz=new double[np];       MemCpuParticles+=sizeof(double)*np;
      Velrhop=new tfloat4[np];   MemCpuParticles+=sizeof(tfloat4)*np;
      AuxPos=new tdouble3[np];   MemCpuParticles+=sizeof(tdouble3)*np; 
      AuxVel=new tfloat3[np];    MemCpuParticles+=sizeof(tfloat3)*np;
      AuxRhop=new float[np];     MemCpuParticles+=sizeof(float)*np;
    }
    catch(const std::bad_alloc){
      RunException(met,fun::PrintStr("Could not allocate the requested memory (np=%u).",np));
    }
  }
}

//==============================================================================
/// Frees GPU memory for the particles.
/// Libera memoria en Gpu para particulas.
//==============================================================================
void JSphGpu::FreeGpuMemoryParticles(){
  GpuParticlesSize=0;
  MemGpuParticles=0;
  ArraysGpu->Reset();
}

//==============================================================================
/// Allocates GPU memory for the particles.
/// Reserva memoria en Gpu para las particulas. 
//==============================================================================
void JSphGpu::AllocGpuMemoryParticles(unsigned np,float over){
  const char* met="AllocGpuMemoryParticles";
  FreeGpuMemoryParticles();
  //-Computes number of particles for which memory will be allocated.
  //-Calcula numero de particulas para las que se reserva memoria.
  const unsigned np2=(over>0? unsigned(over*np): np);
  GpuParticlesSize=np2+PARTICLES_OVERMEMORY_MIN;
  //-Define number or arrays to use. | Establece numero de arrays a usar.
  ArraysGpu->SetArraySize(GpuParticlesSize);
  #ifdef CODE_SIZE4
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,2);  //-code,code2
  #else
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_2B,2);  //-code,code2
  #endif
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,4);  //-idp,ar,viscdt,dcell
  if(TDeltaSph==DELTA_DynamicExt)ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,1);  //-delta
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_12B,1); //-ace
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_16B,4); //-velrhop,posxy
  ArraysGpu->AddArrayCount(JArraysGpu::SIZE_8B,2);  //-posz
  if(TStep==STEP_Verlet){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_16B,1); //-velrhopm1
  }
  else if(TStep==STEP_Symplectic){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_8B,1);  //-poszpre
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_16B,2); //-posxypre,velrhoppre
  }
  if(TVisco==VISCO_LaminarSPS){     
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_24B,2); //-SpsTau,SpsGradvel
  }
  if(CaseNfloat){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,4);  //-FtMasspg
  }
  if(TShifting!=SHIFT_None){
    ArraysGpu->AddArrayCount(JArraysGpu::SIZE_12B,1); //-shiftpos
    if(ShiftTFS)ArraysGpu->AddArrayCount(JArraysGpu::SIZE_4B,1); //-shiftdetectc
  }
  //-Shows the allocated memory.
  MemGpuParticles=ArraysGpu->GetAllocMemoryGpu();
  PrintSizeNp(GpuParticlesSize,MemGpuParticles);
  CheckCudaError(met,"Failed GPU memory allocation.");
}

//==============================================================================
/// Resizes space in GPU memory for particles.
//==============================================================================
void JSphGpu::ResizeGpuMemoryParticles(unsigned npnew){
  npnew=npnew+PARTICLES_OVERMEMORY_MIN;
  //-Saves current data from GPU.
  unsigned    *idp       =SaveArrayGpu(Np,Idpg);
  typecode    *code      =SaveArrayGpu(Np,Codeg);
  unsigned    *dcell     =SaveArrayGpu(Np,Dcellg);
  double2     *posxy     =SaveArrayGpu(Np,Posxyg);
  double      *posz      =SaveArrayGpu(Np,Poszg);
  float4      *velrhop   =SaveArrayGpu(Np,Velrhopg);
  float4      *velrhopm1 =SaveArrayGpu(Np,VelrhopM1g);
  double2     *posxypre  =SaveArrayGpu(Np,PosxyPreg);
  double      *poszpre   =SaveArrayGpu(Np,PoszPreg);
  float4      *velrhoppre=SaveArrayGpu(Np,VelrhopPreg);
  tsymatrix3f *spstau    =SaveArrayGpu(Np,SpsTaug);
  //-Frees pointers.
  ArraysGpu->Free(Idpg);
  ArraysGpu->Free(Codeg);
  ArraysGpu->Free(Dcellg);
  ArraysGpu->Free(Posxyg);
  ArraysGpu->Free(Poszg);
  ArraysGpu->Free(Velrhopg);
  ArraysGpu->Free(VelrhopM1g);
  ArraysGpu->Free(PosxyPreg);
  ArraysGpu->Free(PoszPreg);
  ArraysGpu->Free(VelrhopPreg);
  ArraysGpu->Free(SpsTaug);
  //-Resizes GPU memory allocation.
  const double mbparticle=(double(MemGpuParticles)/(1024*1024))/GpuParticlesSize; //-MB por particula.
  Log->Printf("**JSphGpu: Requesting gpu memory for %u particles: %.1f MB.",npnew,mbparticle*npnew);
  ArraysGpu->SetArraySize(npnew);
  //-Reserve pointers.
  Idpg    =ArraysGpu->ReserveUint();
  Codeg   =ArraysGpu->ReserveTypeCode();
  Dcellg  =ArraysGpu->ReserveUint();
  Posxyg  =ArraysGpu->ReserveDouble2();
  Poszg   =ArraysGpu->ReserveDouble();
  Velrhopg=ArraysGpu->ReserveFloat4();
  if(velrhopm1) VelrhopM1g =ArraysGpu->ReserveFloat4();
  if(posxypre)  PosxyPreg  =ArraysGpu->ReserveDouble2();
  if(poszpre)   PoszPreg   =ArraysGpu->ReserveDouble();
  if(velrhoppre)VelrhopPreg=ArraysGpu->ReserveFloat4();
  if(spstau)    SpsTaug    =ArraysGpu->ReserveSymatrix3f();
  //-Restore data in GPU memory.
  RestoreArrayGpu(Np,idp,Idpg);
  RestoreArrayGpu(Np,code,Codeg);
  RestoreArrayGpu(Np,dcell,Dcellg);
  RestoreArrayGpu(Np,posxy,Posxyg);
  RestoreArrayGpu(Np,posz,Poszg);
  RestoreArrayGpu(Np,velrhop,Velrhopg);
  RestoreArrayGpu(Np,velrhopm1,VelrhopM1g);
  RestoreArrayGpu(Np,posxypre,PosxyPreg);
  RestoreArrayGpu(Np,poszpre,PoszPreg);
  RestoreArrayGpu(Np,velrhoppre,VelrhopPreg);
  RestoreArrayGpu(Np,spstau,SpsTaug);
  //-Updates values.
  GpuParticlesSize=npnew;
  MemGpuParticles=ArraysGpu->GetAllocMemoryGpu();
}

//==============================================================================
/// Saves a GPU array in CPU memory. 
//==============================================================================
template<class T> T* JSphGpu::TSaveArrayGpu(unsigned np,const T *datasrc)const{
  T *data=NULL;
  if(datasrc){
    try{
      data=new T[np];
    }
    catch(const std::bad_alloc){
      RunException("TSaveArrayGpu","Could not allocate the requested memory.");
    }
    cudaMemcpy(data,datasrc,sizeof(T)*np,cudaMemcpyDeviceToHost);
  }
  return(data);
}

//==============================================================================
/// Restores a GPU array (generic) from CPU memory. 
//==============================================================================
template<class T> void JSphGpu::TRestoreArrayGpu(unsigned np,T *data,T *datanew)const{
  if(data&&datanew)cudaMemcpy(datanew,data,sizeof(T)*np,cudaMemcpyHostToDevice);
  delete[] data;
}

//==============================================================================
/// Arrays for basic particle data.
/// Arrays para datos basicos de las particulas. 
//==============================================================================
void JSphGpu::ReserveBasicArraysGpu(){
  Idpg=ArraysGpu->ReserveUint();
  Codeg=ArraysGpu->ReserveTypeCode();
  Dcellg=ArraysGpu->ReserveUint();
  Posxyg=ArraysGpu->ReserveDouble2();
  Poszg=ArraysGpu->ReserveDouble();
  Velrhopg=ArraysGpu->ReserveFloat4();
  if(TStep==STEP_Verlet)VelrhopM1g=ArraysGpu->ReserveFloat4();
  if(TVisco==VISCO_LaminarSPS)SpsTaug=ArraysGpu->ReserveSymatrix3f();
}

//==============================================================================
/// Returns the allocated memory on the CPU.
/// Devuelve la memoria reservada en CPU.
//==============================================================================
llong JSphGpu::GetAllocMemoryCpu()const{  
  llong s=JSph::GetAllocMemoryCpu();
  //-Allocated in AllocCpuMemoryFixed().
  s+=MemCpuFixed;
  //-Allocated in AllocMemoryParticles().
  s+=MemCpuParticles;
  //-Allocated in other objects.
  return(s);
}

//==============================================================================
/// Returns the allocated memory in the GPU.
/// Devuelve la memoria reservada en GPU.
//==============================================================================
llong JSphGpu::GetAllocMemoryGpu()const{  
  llong s=0;
  //-Allocated in AllocGpuMemoryParticles().
  s+=MemGpuParticles;
  //-Allocated in AllocGpuMemoryFixed().
  s+=MemGpuFixed;
  //-Allocated in ther objects.
  return(s);
}

//==============================================================================
/// Displays the allocated memory.
/// Visualiza la memoria reservada.
//==============================================================================
void JSphGpu::PrintAllocMemory(llong mcpu,llong mgpu)const{
  Log->Printf("Allocated memory in CPU: %lld (%.2f MB)",mcpu,double(mcpu)/(1024*1024));
  Log->Printf("Allocated memory in GPU: %lld (%.2f MB)",mgpu,double(mgpu)/(1024*1024));
}

//==============================================================================
/// Uploads constants to the GPU.
/// Copia constantes a la GPU.
//==============================================================================
void JSphGpu::ConstantDataUp(){
  StCteInteraction ctes;
  memset(&ctes,0,sizeof(StCteInteraction));
  ctes.nbound=CaseNbound;
  ctes.massb=MassBound; ctes.massf=MassFluid;
  ctes.fourh2=Fourh2; ctes.h=H;
  if(TKernel==KERNEL_Wendland){
    ctes.awen=Awen; ctes.bwen=Bwen;
  }
  else if(TKernel==KERNEL_Gaussian){
    ctes.agau=Agau; ctes.bgau=Bgau;
  }
  else if(TKernel==KERNEL_Cubic){
    ctes.cubic_a1=CubicCte.a1; ctes.cubic_a2=CubicCte.a2; ctes.cubic_aa=CubicCte.aa; ctes.cubic_a24=CubicCte.a24;
    ctes.cubic_c1=CubicCte.c1; ctes.cubic_c2=CubicCte.c2; ctes.cubic_d1=CubicCte.d1; ctes.cubic_odwdeltap=CubicCte.od_wdeltap;
  }
  ctes.cs0=float(Cs0); ctes.eta2=Eta2;
  ctes.delta2h=Delta2H;
  ctes.scell=Scell; ctes.dosh=Dosh; ctes.dp=float(Dp);
  ctes.cteb=CteB; ctes.gamma=Gamma;
  ctes.rhopzero=RhopZero;
  ctes.ovrhopzero=1.f/RhopZero;
  ctes.movlimit=MovLimit;
  ctes.maprealposminx=MapRealPosMin.x; ctes.maprealposminy=MapRealPosMin.y; ctes.maprealposminz=MapRealPosMin.z;
  ctes.maprealsizex=MapRealSize.x; ctes.maprealsizey=MapRealSize.y; ctes.maprealsizez=MapRealSize.z;
  ctes.periactive=PeriActive;
  ctes.xperincx=PeriXinc.x; ctes.xperincy=PeriXinc.y; ctes.xperincz=PeriXinc.z;
  ctes.yperincx=PeriYinc.x; ctes.yperincy=PeriYinc.y; ctes.yperincz=PeriYinc.z;
  ctes.zperincx=PeriZinc.x; ctes.zperincy=PeriZinc.y; ctes.zperincz=PeriZinc.z;
  ctes.cellcode=DomCellCode;
  ctes.domposminx=DomPosMin.x; ctes.domposminy=DomPosMin.y; ctes.domposminz=DomPosMin.z;
  cusph::CteInteractionUp(&ctes);
  CheckCudaError("ConstantDataUp","Failed copying constants to GPU.");
}

//==============================================================================
/// Uploads particle data to the GPU.
/// Sube datos de particulas a la GPU.
//==============================================================================
void JSphGpu::ParticlesDataUp(unsigned n){
  cudaMemcpy(Idpg    ,Idp    ,sizeof(unsigned)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Codeg   ,Code   ,sizeof(typecode)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Dcellg  ,Dcell  ,sizeof(unsigned)*n,cudaMemcpyHostToDevice);
  cudaMemcpy(Posxyg  ,Posxy  ,sizeof(double2)*n ,cudaMemcpyHostToDevice);
  cudaMemcpy(Poszg   ,Posz   ,sizeof(double)*n  ,cudaMemcpyHostToDevice);
  cudaMemcpy(Velrhopg,Velrhop,sizeof(float4)*n  ,cudaMemcpyHostToDevice);
  CheckCudaError("ParticlesDataUp","Failed copying data to GPU.");
}

//==============================================================================
/// Recovers particle data from the GPU and returns the particle number that
/// are less than n if the paeriodic particles are removed.
/// - code: Recovers data of Codeg.
/// - onlynormal: Onlly retains the normal particles, removes the periodic ones.
///
/// Recupera datos de particulas de la GPU y devuelve el numero de particulas que
/// sera menor que n si se eliminaron las periodicas.
/// - code: Recupera datos de Codeg.
/// - onlynormal: Solo se queda con las normales, elimina las particulas periodicas.
//==============================================================================
unsigned JSphGpu::ParticlesDataDown(unsigned n,unsigned pini,bool code,bool onlynormal){
  unsigned num=n;
  cudaMemcpy(Idp,Idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy(Posxy,Posxyg+pini,sizeof(double2)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy(Posz,Poszg+pini,sizeof(double)*n,cudaMemcpyDeviceToHost);
  cudaMemcpy(Velrhop,Velrhopg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  if(code || onlynormal)cudaMemcpy(Code,Codeg+pini,sizeof(typecode)*n,cudaMemcpyDeviceToHost);
  CheckCudaError("ParticlesDataDown","Failed copying data from GPU.");
  //-Eliminates abnormal particles (periodic and others). | Elimina particulas no normales (periodicas y otras).
  if(onlynormal){
    unsigned ndel=0;
    for(unsigned p=0;p<n;p++){
      const bool normal=CODE_IsNormal(Code[p]);
      if(ndel && normal){
        Idp[p-ndel]    =Idp[p];
        Posxy[p-ndel]  =Posxy[p];
        Posz[p-ndel]   =Posz[p];
        Velrhop[p-ndel]=Velrhop[p];
        Code[p-ndel]   =Code[p];
      }
      if(!normal)ndel++;
    }
    num-=ndel;
  }
  //-Converts data to a simple format. | Convierte datos a formato simple.
  for(unsigned p=0;p<n;p++){
    AuxPos[p]=TDouble3(Posxy[p].x,Posxy[p].y,Posz[p]);
    AuxVel[p]=TFloat3(Velrhop[p].x,Velrhop[p].y,Velrhop[p].z);
    AuxRhop[p]=Velrhop[p].w;
  }
  return(num);
}

//==============================================================================
/// Initialises CUDA device.
/// Inicializa dispositivo CUDA.
//==============================================================================
void JSphGpu::SelecDevice(int gpuid){
  const char* met="SelecDevice";
  Log->Print("[Select CUDA Device]");
  //-Get and show GPU information.
  vector<string> gpuinfo;
  vector<fcuda::StGpuInfo> gpuprops;
  const int devcount=fcuda::GetCudaDevicesInfo(&gpuinfo,&gpuprops);
  Log->Print(gpuinfo);
  Log->Print(" ");
  //-GPU selection.
  GpuSelect=-1;
  if(devcount){
    if(gpuid>=0)cudaSetDevice(gpuid);
    else{
      unsigned *ptr=NULL;
      cudaMalloc((void**)&ptr,sizeof(unsigned)*100);
      cudaFree(ptr);
    }
    cudaDeviceProp devp;
    int dev;
    cudaGetDevice(&dev);
    cudaGetDeviceProperties(&devp,dev);
    GpuSelect=dev;
    GpuName=devp.name;
    GpuGlobalMem=devp.totalGlobalMem;
    GpuSharedMem=int(devp.sharedMemPerBlock);
    GpuCompute=devp.major*10+devp.minor;
    //-Displays information on the selected hardware.
    //-Muestra informacion del hardware seleccionado.
    Log->Print("[GPU Hardware]");
    if(gpuid<0)Hardware=fun::PrintStr("Gpu_%d?=\"%s\"",GpuSelect,GpuName.c_str());
    else Hardware=fun::PrintStr("Gpu_%d=\"%s\"",GpuSelect,GpuName.c_str());
    if(gpuid<0)Log->Printf("Device default: %d  \"%s\"",GpuSelect,GpuName.c_str());
    else Log->Printf("Device selected: %d  \"%s\"",GpuSelect,GpuName.c_str());
    Log->Printf("Compute capability: %.1f",float(GpuCompute)/10);
    Log->Printf("Memory global: %d MB",int(GpuGlobalMem/(1024*1024)));
    Log->Printf("Memory shared: %u Bytes",GpuSharedMem);
  }
  else RunException(met, "There are no available CUDA devices.");
}

//==============================================================================
/// Configures BlockSizeMode to compute optimum size of CUDA blocks.
//==============================================================================
void JSphGpu::ConfigBlockSizes(bool usezone,bool useperi){
  const char met[]="ConfigBlockSizes";
  Log->Print(" ");
  BlockSizesStr="";
  if(CellMode==CELLMODE_2H || CellMode==CELLMODE_H){
    const bool lamsps=(TVisco==VISCO_LaminarSPS);
    const bool shift=(TShifting!=SHIFT_None);
    BlockSizes.forcesbound=BlockSizes.forcesfluid=BlockSizes.forcesdem=BSIZE_FIXED;
    //-Collects kernel information.
    StKerInfo kerinfo;
    memset(&kerinfo,0,sizeof(StKerInfo));
    #ifndef DISABLE_BSMODES
      cusph::Interaction_Forces(Psingle,TKernel,FtMode,lamsps,TDeltaSph,CellMode,0,0,0,0,100,50,20,TUint3(0),NULL,TUint3(0),NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,TShifting,NULL,NULL,Simulate2D,&kerinfo,NULL);
      if(UseDEM)cusph::Interaction_ForcesDem(Psingle,CellMode,BlockSizes.forcesdem,CaseNfloat,TUint3(0),NULL,TUint3(0),NULL,NULL,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,&kerinfo);
    #endif
    //Log->Printf("====> bound -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesbound_rg,kerinfo.forcesbound_bs,kerinfo.forcesbound_bsmax);
    //Log->Printf("====> fluid -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesfluid_rg,kerinfo.forcesfluid_bs,kerinfo.forcesfluid_bsmax);
    //Log->Printf("====> dem   -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesdem_rg,kerinfo.forcesdem_bs,kerinfo.forcesdem_bsmax);
    //-Defines blocsize according BlockSizeMode.
    if(BlockSizeMode==BSIZEMODE_Occupancy){
      if(!kerinfo.forcesbound_bs || !kerinfo.forcesfluid_bs){
        Log->Printf("**BlockSize calculation mode %s is invalid.",GetNameBlockSizeMode(BlockSizeMode));
        BlockSizeMode=BSIZEMODE_Fixed;
      }
      else{
        if(kerinfo.forcesbound_bs)BlockSizes.forcesbound=kerinfo.forcesbound_bs;
        if(kerinfo.forcesfluid_bs)BlockSizes.forcesfluid=kerinfo.forcesfluid_bs;
        if(kerinfo.forcesdem_bs)BlockSizes.forcesdem=kerinfo.forcesdem_bs;
      }
    }
    if(BlockSizeMode==BSIZEMODE_Empirical){
      BsAuto=new JBlockSizeAuto(Log,500);
      BsAuto->AddKernel("KerInteractionForcesFluid",64,31,32,BSIZE_FIXED);  //15:512 31:1024
      BsAuto->AddKernel("KerInteractionForcesBound",64,31,32,BSIZE_FIXED);
      BsAuto->AddKernel("KerInteractionForcesDem",64,31,32,BSIZE_FIXED);
      if(kerinfo.forcesdem_bs)BlockSizes.forcesdem=kerinfo.forcesdem_bs;
    }
    Log->Printf("BlockSize calculation mode: %s.",GetNameBlockSizeMode(BlockSizeMode));
    string txrb=(kerinfo.forcesbound_rg? fun::PrintStr("(%d regs)",kerinfo.forcesbound_rg): string("(? regs)"));
    string txrf=(kerinfo.forcesbound_rg? fun::PrintStr("(%d regs)",kerinfo.forcesfluid_rg): string("(? regs)"));
    string txrd=(kerinfo.forcesdem_rg  ? fun::PrintStr("(%d regs)",kerinfo.forcesdem_rg  ): string("(? regs)"));
    string txb=string("BsForcesBound=")+(BlockSizeMode==BSIZEMODE_Empirical? string("Dynamic"): fun::IntStr(BlockSizes.forcesbound))+" "+txrb;
    string txf=string("BsForcesFluid=")+(BlockSizeMode==BSIZEMODE_Empirical? string("Dynamic"): fun::IntStr(BlockSizes.forcesfluid))+" "+txrf;
    string txd=string("BsForcesDem="  )+fun::IntStr(BlockSizes.forcesdem)+" "+txrd;
    Log->Print(string("  ")+txb);
    Log->Print(string("  ")+txf);
    if(UseDEM)Log->Print(string("  ")+txd);
    if(!BlockSizesStr.empty())BlockSizesStr=BlockSizesStr+" - ";
    BlockSizesStr=BlockSizesStr+txb+" - "+txf;
    if(UseDEM)BlockSizesStr=BlockSizesStr+" - "+txd;
  }
  else RunException(met,"CellMode unrecognised.");
  Log->Print(" ");
}

//==============================================================================
/// Configures execution mode in the GPU.
/// Configura modo de ejecucion en GPU.
//==============================================================================
void JSphGpu::ConfigRunMode(std::string preinfo){
  //#ifndef WIN32  //-Error compilation when gcc5 is used.
  //  const int len=128; char hname[len];
  //  gethostname(hname,len);
  //  if(!preinfo.empty())preinfo=preinfo+", ";
  //  preinfo=preinfo+"HostName:"+hname;
  //#endif
  RunMode=preinfo+RunMode;
  if(Stable)RunMode=string("Stable - ")+RunMode;
  if(Psingle)RunMode=string("Pos-Single - ")+RunMode;
  else RunMode=string("Pos-Double - ")+RunMode;
  Log->Print(" ");
  Log->Print(fun::VarStr("RunMode",RunMode));
  Log->Print(" ");
  if(!RunMode.empty())RunMode=RunMode+" - "+BlockSizesStr;
}

//==============================================================================
/// Adjusts variables of particles of floating bodies.
//==============================================================================
void JSphGpu::InitFloating(){
  if(PartBegin){
    JPartFloatBi4Load ftdata;
    ftdata.LoadFile(PartBeginDir);
    //-Checks if the constant data match.
    //-Comprueba coincidencia de datos constantes.
    for(unsigned cf=0;cf<FtCount;cf++)ftdata.CheckHeadData(cf,FtObjs[cf].mkbound,FtObjs[cf].begin,FtObjs[cf].count,FtObjs[cf].mass);
    //-Loads PART data.
    ftdata.LoadPart(PartBegin);
    for(unsigned cf=0;cf<FtCount;cf++){
      FtObjs[cf].center=ftdata.GetPartCenter(cf);
      FtObjs[cf].fvel=ftdata.GetPartFvel(cf);
      FtObjs[cf].fomega=ftdata.GetPartFomega(cf);
      FtObjs[cf].radius=ftdata.GetHeadRadius(cf);
    }
    DemDtForce=ftdata.GetPartDemDtForce();
  }
  //-Copies massp values to GPU.
  {
    float *massp=new float[FtCount];
    for(unsigned cf=0;cf<FtCount;cf++)massp[cf]=FtObjs[cf].massp;
    cudaMemcpy(FtoMasspg,massp,sizeof(float)*FtCount,cudaMemcpyHostToDevice);
    delete[] massp;
  }
  //-Copies floating values to GPU.
  {
    typedef struct{
      unsigned pini;
      unsigned np;
      float radius;
      float mass;
    }stdata;

    stdata   *data  =new stdata  [FtCount];
    tdouble3 *center=new tdouble3[FtCount];
    tfloat3  *angles=new tfloat3 [FtCount];
    tfloat3  *vel   =new tfloat3 [FtCount];
    tfloat3  *omega =new tfloat3 [FtCount];
    tfloat4  *inert8=new tfloat4 [FtCount*2];
    float    *inert1=new float   [FtCount];
    for(unsigned cf=0;cf<FtCount;cf++){
      data[cf].pini=FtObjs[cf].begin-CaseNpb;
      data[cf].np=FtObjs[cf].count;
      data[cf].radius=FtObjs[cf].radius;
      data[cf].mass=FtObjs[cf].mass;
      center[cf]=FtObjs[cf].center;
      angles[cf]=FtObjs[cf].angles;
      vel   [cf]=FtObjs[cf].fvel;
      omega [cf]=FtObjs[cf].fomega;
      const tmatrix3f v=FtObjs[cf].inertiaini;
      inert8[cf*2]  =TFloat4(v.a11,v.a12,v.a13,v.a21);
      inert8[cf*2+1]=TFloat4(v.a22,v.a23,v.a31,v.a32);
      inert1[cf]    =v.a33;
    }
    cudaMemcpy(FtoDatag       ,data  ,sizeof(float4)   *FtCount  ,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoCenterg     ,center,sizeof(double3)  *FtCount  ,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoAnglesg     ,angles,sizeof(float3)   *FtCount  ,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoVelg        ,vel   ,sizeof(float3)   *FtCount  ,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoOmegag      ,omega ,sizeof(float3)   *FtCount  ,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoInertiaini8g,inert8,sizeof(float4)   *FtCount*2,cudaMemcpyHostToDevice);
    cudaMemcpy(FtoInertiaini1g,inert1,sizeof(float)    *FtCount  ,cudaMemcpyHostToDevice);
    delete[] data;
    delete[] center;
    delete[] angles;
    delete[] vel;
    delete[] omega;
    delete[] inert8;
    delete[] inert1;
  }
  //-Copies data object for DEM to GPU.
  if(UseDEM){ //(DEM)
    float4 *data=new float4[DemDataSize];
    for(unsigned c=0;c<DemDataSize;c++){
      data[c].x=DemData[c].mass;
      data[c].y=DemData[c].tau;
      data[c].z=DemData[c].kfric;
      data[c].w=DemData[c].restitu;
    }
    cudaMemcpy(DemDatag,data,sizeof(float4)*DemDataSize,cudaMemcpyHostToDevice);
    delete[] data;
  }
}

//==============================================================================
/// Initialises arrays and variables for the execution.
/// Inicializa vectores y variables para la ejecucion.
//==============================================================================
void JSphGpu::InitRunGpu(){
  ParticlesDataDown(Np,0,false,false);
  InitRun(Np,Idp,AuxPos);

  if(TStep==STEP_Verlet)cudaMemcpy(VelrhopM1g,Velrhopg,sizeof(float4)*Np,cudaMemcpyDeviceToDevice);
  if(TVisco==VISCO_LaminarSPS)cudaMemset(SpsTaug,0,sizeof(tsymatrix3f)*Np);
  if(CaseNfloat)InitFloating();
  CheckCudaError("InitRunGpu","Failed initializing variables for execution.");
}

//==============================================================================
/// Adds variable acceleration from input files.
//==============================================================================
void JSphGpu::AddAccInput(){
  for(unsigned c=0;c<AccInput->GetCount();c++){
    unsigned mkfluid;
    tdouble3 acclin,accang,centre,velang,vellin;
    bool setgravity;
    AccInput->GetAccValues(c,TimeStep,mkfluid,acclin,accang,centre,velang,vellin,setgravity);
    const typecode codesel=typecode(mkfluid);
    cusph::AddAccInput(Np-Npb,Npb,codesel,acclin,accang,centre,velang,vellin,setgravity,Gravity,Codeg,Posxyg,Poszg,Velrhopg,Aceg);
  }
}

//==============================================================================
/// Prepares variables for interaction.
/// Prepara variables para interaccion.
//==============================================================================
void JSphGpu::PreInteractionVars_Forces(unsigned np,unsigned npb){
  //-Initialises arrays.
  const unsigned npf=np-npb;
  cudaMemset(ViscDtg,0,sizeof(float)*np);                                //ViscDtg[]=0
  cudaMemset(Arg,0,sizeof(float)*np);                                    //Arg[]=0
  if(Deltag)cudaMemset(Deltag,0,sizeof(float)*np);                       //Deltag[]=0
  if(ShiftPosg)cudaMemset(ShiftPosg,0,sizeof(tfloat3)*np);               //ShiftPosg[]=0
  if(ShiftDetectg)cudaMemset(ShiftDetectg,0,sizeof(float)*np);           //ShiftDetectg[]=0
  cudaMemset(Aceg,0,sizeof(tfloat3)*npb);                                //Aceg[]=(0,0,0) para bound //Aceg[]=(0,0,0) for the boundary
  cusph::InitArray(npf,Aceg+npb,Gravity);                                //Aceg[]=Gravity para fluid //Aceg[]=Gravity for the fluid
  if(SpsGradvelg)cudaMemset(SpsGradvelg+npb,0,sizeof(tsymatrix3f)*npf);  //SpsGradvelg[]=(0,0,0,0,0,0).

  //-Apply the extra forces to the correct particle sets.
  if(AccInput)AddAccInput();
}

//==============================================================================
/// Prepares variables for interaction.
/// Prepara variables para interaccion.
//==============================================================================
void JSphGpu::PreInteraction_Forces(){
  TmgStart(Timers,TMG_CfPreForces);
  //-Allocates memory.
  ViscDtg=ArraysGpu->ReserveFloat();
  Arg=ArraysGpu->ReserveFloat();
  Aceg=ArraysGpu->ReserveFloat3();
  if(TDeltaSph==DELTA_DynamicExt)Deltag=ArraysGpu->ReserveFloat();
  if(TShifting!=SHIFT_None){
    ShiftPosg=ArraysGpu->ReserveFloat3();
    if(ShiftTFS)ShiftDetectg=ArraysGpu->ReserveFloat();
  }   
  if(TVisco==VISCO_LaminarSPS)SpsGradvelg=ArraysGpu->ReserveSymatrix3f();

  //-Prepares data for interation Pos-Single.
  if(Psingle){
    PsPospressg=ArraysGpu->ReserveFloat4();
    cusph::PreInteractionSingle(Np,Posxyg,Poszg,Velrhopg,PsPospressg,CteB,Gamma);
  }
  //-Initialises arrays.
  PreInteractionVars_Forces(Np,Npb);

  //-Computes VelMax: Includes the particles from floating bodies and does not affect the periodic conditions.
  //-Calcula VelMax: Se incluyen las particulas floatings y no afecta el uso de condiciones periodicas.
  const unsigned pini=(DtAllParticles? 0: Npb);
  cusph::ComputeVelMod(Np-pini,Velrhopg+pini,ViscDtg);
  float velmax=cusph::ReduMaxFloat(Np-pini,0,ViscDtg,CellDiv->GetAuxMem(cusph::ReduMaxFloatSize(Np-pini)));
  VelMax=sqrt(velmax);
  cudaMemset(ViscDtg,0,sizeof(float)*Np);           //ViscDtg[]=0
  ViscDtMax=0;
  CheckCudaError("PreInteraction_Forces","Failed calculating VelMax.");
  TmgStop(Timers,TMG_CfPreForces);
}

//==============================================================================
/// Frees memory allocated for ArraysGpu.
/// Libera memoria asignada de ArraysGpu.
//==============================================================================
void JSphGpu::PosInteraction_Forces(){
  //-Frees memory allocated in PreInteraction_Forces().
  ArraysGpu->Free(Arg);          Arg=NULL;
  ArraysGpu->Free(Aceg);         Aceg=NULL;
  ArraysGpu->Free(ViscDtg);      ViscDtg=NULL;
  ArraysGpu->Free(Deltag);       Deltag=NULL;
  ArraysGpu->Free(ShiftPosg);    ShiftPosg=NULL;
  ArraysGpu->Free(ShiftDetectg); ShiftDetectg=NULL;
  ArraysGpu->Free(PsPospressg);  PsPospressg=NULL;
  ArraysGpu->Free(SpsGradvelg);  SpsGradvelg=NULL;
}

//==============================================================================
/// Updates particles according to forces and dt using Verlet.
/// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphGpu::ComputeVerlet(double dt){  //pdtedom
  TmgStart(Timers,TMG_SuComputeStep);
  const bool floatings=(WithFloating);
  const bool shift=TShifting!=SHIFT_None;
  VerletStep++;
  //-Allocates memory to compute the displacement.
  //-Asigna memoria para calcular el desplazamiento.
  double2 *movxyg=ArraysGpu->ReserveDouble2();
  double *movzg=ArraysGpu->ReserveDouble();
  //-Computes displacement, velocity and density.
  //-Calcula desplazamiento, velocidad y densidad.
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    cusph::ComputeStepVerlet(floatings,shift,Np,Npb,Velrhopg,VelrhopM1g,Arg,Aceg,ShiftPosg,dt,twodt,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,VelrhopM1g);
  }
  else{
    cusph::ComputeStepVerlet(floatings,shift,Np,Npb,Velrhopg,Velrhopg,Arg,Aceg,ShiftPosg,dt,dt,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,VelrhopM1g);
    VerletStep=0;
  }
  //-The new values are calculated in VelRhopM1g.
  //-Los nuevos valores se calculan en VelrhopM1g.
  swap(Velrhopg,VelrhopM1g);   //-Exchanges Velrhopg and VelrhopM1g. | Intercambia Velrhopg y VelrhopM1g.
  //-Applies displacement to non-periodic fluid particles.
  //-Aplica desplazamiento a las particulas fluid no periodicas.
  cusph::ComputeStepPos(PeriActive,floatings,Np,Npb,movxyg,movzg,Posxyg,Poszg,Dcellg,Codeg);
  //-Frees memory allocated for the diplacement.
  ArraysGpu->Free(movxyg);   movxyg=NULL;
  ArraysGpu->Free(movzg);    movzg=NULL;
  TmgStop(Timers,TMG_SuComputeStep);
}

//==============================================================================
/// Updates particles according to forces and dt using Symplectic-Predictor.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
void JSphGpu::ComputeSymplecticPre(double dt){
  TmgStart(Timers,TMG_SuComputeStep);
  const bool floatings=(WithFloating);
  const bool shift=(TShifting!=SHIFT_None);
  //-Allocates memory to PRE variables.
  PosxyPreg=ArraysGpu->ReserveDouble2();
  PoszPreg=ArraysGpu->ReserveDouble();
  VelrhopPreg=ArraysGpu->ReserveFloat4();
  //-Changes data of PRE variables for calculating the new data.
  //-Cambia datos a variables PRE para calcular nuevos datos.
  swap(PosxyPreg,Posxyg);      //-PosxyPre[] <= Posxy[]
  swap(PoszPreg,Poszg);        //-PoszPre[] <= Posz[]
  swap(VelrhopPreg,Velrhopg);  //-VelrhopPre[] <= Velrhop[]
  //-Allocate memory to compute the diplacement.
  double2 *movxyg=ArraysGpu->ReserveDouble2();
  double *movzg=ArraysGpu->ReserveDouble();
  //-Compute displacement, velocity and density.
  const double dt05=dt*.5;
  cusph::ComputeStepSymplecticPre(floatings,shift,Np,Npb,VelrhopPreg,Arg,Aceg,ShiftPosg,dt05,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,Velrhopg);
  //-Applies displacement to non-periodic fluid particles.
  //-Aplica desplazamiento a las particulas fluid no periodicas.
  cusph::ComputeStepPos2(PeriActive,floatings,Np,Npb,PosxyPreg,PoszPreg,movxyg,movzg,Posxyg,Poszg,Dcellg,Codeg);
  //-Frees memory allocated for the displacement
  ArraysGpu->Free(movxyg);   movxyg=NULL;
  ArraysGpu->Free(movzg);    movzg=NULL;
  //-Copies previous position of the boundaries.
  //-Copia posicion anterior del contorno.
  cudaMemcpy(Posxyg,PosxyPreg,sizeof(double2)*Npb,cudaMemcpyDeviceToDevice);
  cudaMemcpy(Poszg,PoszPreg,sizeof(double)*Npb,cudaMemcpyDeviceToDevice);
  TmgStop(Timers,TMG_SuComputeStep);
}

//==============================================================================
/// Updates particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphGpu::ComputeSymplecticCorr(double dt){
  TmgStart(Timers,TMG_SuComputeStep);
  const bool floatings=(WithFloating);
  const bool shift=(TShifting!=SHIFT_None);
  //-Allocates memory to calculate the displacement.
  double2 *movxyg=ArraysGpu->ReserveDouble2();
  double *movzg=ArraysGpu->ReserveDouble();
  //-Computes displacement, velocity and density.
  const double dt05=dt*.5;
  cusph::ComputeStepSymplecticCor(floatings,shift,Np,Npb,VelrhopPreg,Arg,Aceg,ShiftPosg,dt05,dt,RhopOutMin,RhopOutMax,Codeg,movxyg,movzg,Velrhopg);
  //-Applies displacement to non-periodic fluid particles.
  //-Aplica desplazamiento a las particulas fluid no periodicas.
  cusph::ComputeStepPos2(PeriActive,floatings,Np,Npb,PosxyPreg,PoszPreg,movxyg,movzg,Posxyg,Poszg,Dcellg,Codeg);
  //-Frees memory allocated for diplacement.
  ArraysGpu->Free(movxyg);   movxyg=NULL;
  ArraysGpu->Free(movzg);    movzg=NULL;
  //-Frees memory allocated for the predictor variables in ComputeSymplecticPre().
  ArraysGpu->Free(PosxyPreg);    PosxyPreg=NULL;
  ArraysGpu->Free(PoszPreg);     PoszPreg=NULL;
  ArraysGpu->Free(VelrhopPreg);  VelrhopPreg=NULL;
  TmgStop(Timers,TMG_SuComputeStep);
}

//==============================================================================
/// Computes the variable Dt.
/// Calcula un Dt variable.
//==============================================================================
double JSphGpu::DtVariable(bool final){
  //-dt1 depends on force per unit mass.
  const double acemax=sqrt(double(AceMax));
  const double dt1=(AceMax? (sqrt(double(H)/AceMax)): DBL_MAX); 
  //-dt2 combines the Courant and the viscous time-step controls.
  const double dt2=double(H)/(max(Cs0,VelMax*10.)+double(H)*ViscDtMax);
  //-dt new value of time step.
  double dt=double(CFLnumber)*min(dt1,dt2);
  if(DtFixed)dt=DtFixed->GetDt(float(TimeStep),float(dt));
  if(dt<double(DtMin)){ 
    dt=double(DtMin); DtModif++;
    if(DtModif>=DtModifWrn){
      Log->PrintfWarning("%d DTs adjusted to DtMin (t:%g, nstep:%u)",DtModif,TimeStep,Nstep);
      DtModifWrn*=10;
    }
  }
  if(SaveDt && final)SaveDt->AddValues(TimeStep,dt,dt1*CFLnumber,dt2*CFLnumber,AceMax,ViscDtMax,VelMax);
  return(dt);
}

//==============================================================================
/// Computes final shifting distance for the particle position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphGpu::RunShifting(double dt){
  TmgStart(Timers,TMG_SuShifting);
  const double coeftfs=(Simulate2D? 2.0: 3.0)-ShiftTFS;
  cusph::RunShifting(Np,Npb,dt,ShiftCoef,ShiftTFS,coeftfs,Velrhopg,ShiftDetectg,ShiftPosg);
  TmgStop(Timers,TMG_SuShifting);
}

//==============================================================================
/// Calculates predefined movement of boundary particles.
/// Calcula movimiento predefinido de boundary particles.
//==============================================================================
void JSphGpu::CalcMotion(double stepdt){
  TmgStart(Timers,TMG_SuMotion);
  const bool motsim=true;
  const JSphMotion::TpMotionMode mode=(motsim? JSphMotion::MOMT_Simple: JSphMotion::MOMT_Ace2dt);
  SphMotion->ProcesTime(mode,TimeStep,stepdt);
}

//==============================================================================
/// Processes boundary particle movement.
/// Procesa movimiento de boundary particles.
//==============================================================================
void JSphGpu::RunMotion(double stepdt){
  const char met[]="RunMotion";
  TmgStart(Timers,TMG_SuMotion);
  const bool motsim=true;
  BoundChanged=false;
  if(SphMotion->GetActiveMotion()){
    cusph::CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codeg,Idpg,RidpMoveg);
    BoundChanged=true;
    bool typesimple;
    tdouble3 simplemov,simplevel,simpleace;
    tmatrix4d matmov,matmov2;
    unsigned nparts,idbegin;
    const unsigned nref=SphMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++)if(SphMotion->ProcesTimeGetData(ref,typesimple,simplemov,simplevel,simpleace,matmov,matmov2,nparts,idbegin)){
      const unsigned pini=idbegin-CaseNfixed;
      if(typesimple){//-Simple movement. | Movimiento simple.
        if(Simulate2D)simplemov.y=simplevel.y=simpleace.y=0;
        if(motsim)cusph::MoveLinBound(PeriActive,nparts,pini,simplemov,ToTFloat3(simplevel),RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
        //else    cusph::MoveLinBoundAce(PeriActive,nparts,pini,simplemov,ToTFloat3(simplevel),ToTFloat3(simpleace),RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
      }
      else{//-Movement using a matrix. | Movimiento con matriz.
        if(motsim)cusph::MoveMatBound   (PeriActive,Simulate2D,nparts,pini,matmov,stepdt,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
        //else    cusph::MoveMatBoundAce(PeriActive,Simulate2D,nparts,pini,matmov,matmov2,stepdt,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
      }
    }
  }
  //-Process other modes of motion. | Procesa otros modos de motion.
  if(WaveGen){
    if(!BoundChanged)cusph::CalcRidp(PeriActive!=0,Npb,0,CaseNfixed,CaseNfixed+CaseNmoving,Codeg,Idpg,RidpMoveg);
    BoundChanged=true;
    //-Control of wave generation (WaveGen). | Gestion de WaveGen.
    for(unsigned c=0;c<WaveGen->GetCount();c++){
      bool typesimple;
      tdouble3 simplemov,simplevel,simpleace;
      tmatrix4d matmov,matmov2;
      unsigned nparts,idbegin;
      //-Get movement data.
      const bool svdata=(TimeStep+stepdt>=TimePartNext);
      if(motsim)typesimple=WaveGen->GetMotion   (svdata,c,TimeStep,stepdt,simplemov,simplevel,matmov,nparts,idbegin);
      else      typesimple=WaveGen->GetMotionAce(svdata,c,TimeStep,stepdt,simplemov,simplevel,simpleace,matmov,matmov2,nparts,idbegin);
      //-Applies movement to paddle particles.
      const unsigned np=nparts,pini=idbegin-CaseNfixed;
      if(typesimple){//-Simple movement. | Movimiento simple.
        if(Simulate2D)simplemov.y=simplevel.y=simpleace.y=0;
        if(motsim)cusph::MoveLinBound(PeriActive,np,pini,simplemov,ToTFloat3(simplevel),RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
        //else    cusph::MoveLinBoundAce(PeriActive,np,pini,simplemov,ToTFloat3(simplevel),ToTFloat3(simpleace),RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
      }
      else{
        if(motsim)cusph::MoveMatBound   (PeriActive,Simulate2D,np,pini,matmov,stepdt,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
        //else    cusph::MoveMatBoundAce(PeriActive,Simulate2D,np,pini,matmov,matmov2,stepdt,RidpMoveg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);
      }
    }
  }
  TmgStop(Timers,TMG_SuMotion);
}

//==============================================================================
/// Applies Damping to indicated particles.
/// Aplica Damping a las particulas indicadas.
//==============================================================================
void JSphGpu::RunDamping(double dt,unsigned np,unsigned npb,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop){
  for(unsigned c=0;c<Damping->GetCount();c++){
    const JDamping::StDamping* da=Damping->GetDampingZone(c);
    const tdouble4 plane=da->plane;
    const float dist=da->dist;
    const float over=da->overlimit;
    const tfloat3 factorxyz=da->factorxyz;
    const float redumax=da->redumax;
    if(!da->usedomain){
      if(CaseNfloat || PeriActive)cusph::ComputeDamping(dt,plane,dist,over,factorxyz,redumax,np-npb,npb,posxy,posz,code,velrhop);
      else cusph::ComputeDamping(dt,plane,dist,over,factorxyz,redumax,np-npb,npb,posxy,posz,NULL,velrhop);
    }
    else{
      const double zmin=da->domzmin;
      const double zmax=da->domzmax;
      const tdouble4 pla0=da->dompla0;
      const tdouble4 pla1=da->dompla1;
      const tdouble4 pla2=da->dompla2;
      const tdouble4 pla3=da->dompla3;
      if(CaseNfloat || PeriActive)cusph::ComputeDampingPla(dt,plane,dist,over,factorxyz,redumax,zmin,zmax,pla0,pla1,pla2,pla3,np-npb,npb,posxy,posz,code,velrhop);
      else cusph::ComputeDampingPla(dt,plane,dist,over,factorxyz,redumax,zmin,zmax,pla0,pla1,pla2,pla3,np-npb,npb,posxy,posz,NULL,velrhop);
    }
  }
}

//==============================================================================
/// Displays the active timers.
/// Muestra los temporizadores activos.
//==============================================================================
void JSphGpu::ShowTimers(bool onlyfile){
  JLog2::TpMode_Out mode=(onlyfile? JLog2::Out_File: JLog2::Out_ScrFile);
  Log->Print("[GPU Timers]",mode);
  if(!SvTimers)Log->Print("none",mode);
  else for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c))Log->Print(TimerToText(c),mode);
}

//==============================================================================
/// Returns string with numbers and values of the active timers. 
/// Devuelve string con nombres y valores de los timers activos.
//==============================================================================
void JSphGpu::GetTimersInfo(std::string &hinfo,std::string &dinfo)const{
  for(unsigned c=0;c<TimerGetCount();c++)if(TimerIsActive(c)){
    hinfo=hinfo+";"+TimerGetName(c);
    dinfo=dinfo+";"+fun::FloatStr(TimerGetValue(c)/1000.f);
  }
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin
 ,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg
 ,const float4 *velrhopg)const
{
  const unsigned np=pfin-pini;
  //-Copy data to CPU memory.
  tdouble3 *posh=fcuda::ToHostPosd3(pini,np,posxyg,poszg);
  #ifdef CODE_SIZE4
    typecode *codeh=fcuda::ToHostUint(pini,np,codeg);
  #else
    typecode *codeh=fcuda::ToHostWord(pini,np,codeg);
  #endif
  unsigned *idph=fcuda::ToHostUint(pini,np,idpg);
  tfloat4  *velrhoph=fcuda::ToHostFloat4(pini,np,velrhopg);
  //-Creates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,np,posh,codeh,idph,velrhoph);
  //-Deallocates memory.
  delete[] posh;  
  delete[] codeh;
  delete[] idph;
  delete[] velrhoph;
}

//==============================================================================
/// Saves VTK file with particle data (debug).
/// Graba fichero vtk con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,unsigned cellcode
  ,const double2 *posxyg,const double *poszg,const unsigned *idpg,const unsigned *dcelg
  ,const typecode *codeg,const float4 *velrhopg,const float4 *velrhopm1g,const float3 *aceg)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  //-Loads position.
  tfloat3 *pos=new tfloat3[n];
  {
    tdouble2 *pxy=new tdouble2[n];
    double *pz=new double[n];
    cudaMemcpy(pxy,posxyg+pini,sizeof(double2)*n,cudaMemcpyDeviceToHost);
    cudaMemcpy(pz,poszg+pini,sizeof(double)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)pos[p]=TFloat3(float(pxy[p].x),float(pxy[p].y),float(pz[p]));
    delete[] pxy;
    delete[] pz;
  }
  //-Loads idp.
  unsigned *idp=NULL;
  if(idpg){
    idp=new unsigned[n];
    cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  }
  //-Loads dcel.
  tuint3 *dcel=NULL;
  if(dcelg){
    dcel=new tuint3[n];
    unsigned *aux=new unsigned[n];
    cudaMemcpy(aux,dcelg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)dcel[p]=TUint3(unsigned(PC__Cellx(cellcode,aux[p])),unsigned(PC__Celly(cellcode,aux[p])),unsigned(PC__Cellz(cellcode,aux[p])));
    delete[] aux;
  }
  //-Loads vel and rhop.
  tfloat3 *vel=NULL;
  float *rhop=NULL;
  if(velrhopg){
    vel=new tfloat3[n];
    rhop=new float[n];
    tfloat4 *aux=new tfloat4[n];
    cudaMemcpy(aux,velrhopg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ vel[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rhop[p]=aux[p].w; }
    delete[] aux;
  }
  //-Loads velm1 and rhopm1.
  tfloat3 *velm1=NULL;
  float *rhopm1=NULL;
  if(velrhopm1g){
    velm1=new tfloat3[n];
    rhopm1=new float[n];
    tfloat4 *aux=new tfloat4[n];
    cudaMemcpy(aux,velrhopm1g+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ velm1[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rhopm1[p]=aux[p].w; }
    delete[] aux;
  }
  //-Loads ace.
  tfloat3 *ace=NULL;
  if(aceg){
    ace=new tfloat3[n];
    cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  }
  //-Loads type.
  byte *type=NULL;
  if(codeg){
    type=new byte[n];
    typecode *aux=new typecode[n];
    cudaMemcpy(aux,codeg+pini,sizeof(typecode)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ 
      const typecode cod=aux[p];
      byte tp=99;
      if(CODE_IsFixed(cod))tp=0;
      else if(CODE_IsMoving(cod))tp=1;
      else if(CODE_IsFloating(cod))tp=2;
      else if(CODE_IsFluid(cod))tp=3;
      if(CODE_IsNormal(cod))tp+=0;
      else if(CODE_IsPeriodic(cod))tp+=10;
      else if(CODE_IsOutIgnore(cod))tp+=20;
      else tp+=30;
      type[p]=tp;
    }
    delete[] aux;
  }

  //-Deifines fields.
  JFormatFiles2::StScalarData fields[10];
  unsigned nfields=0;
  if(idp){    fields[nfields]=JFormatFiles2::DefineField("Idp"   ,JFormatFiles2::UInt32  ,1,idp);    nfields++; }
  if(dcel){   fields[nfields]=JFormatFiles2::DefineField("Dcel"  ,JFormatFiles2::UInt32  ,3,dcel);   nfields++; }
  if(vel){    fields[nfields]=JFormatFiles2::DefineField("Vel"   ,JFormatFiles2::Float32 ,3,vel);    nfields++; }
  if(rhop){   fields[nfields]=JFormatFiles2::DefineField("Rhop"  ,JFormatFiles2::Float32 ,1,rhop);   nfields++; }
  if(velm1){  fields[nfields]=JFormatFiles2::DefineField("Velm1" ,JFormatFiles2::Float32 ,3,velm1);  nfields++; }
  if(rhopm1){ fields[nfields]=JFormatFiles2::DefineField("Rhopm1",JFormatFiles2::Float32 ,1,rhopm1); nfields++; }
  if(ace){    fields[nfields]=JFormatFiles2::DefineField("Ace"   ,JFormatFiles2::Float32 ,3,ace);    nfields++; }
  if(type){   fields[nfields]=JFormatFiles2::DefineField("Typex" ,JFormatFiles2::UChar8  ,1,type);   nfields++; }

  //-Generates file.
  JFormatFiles2::SaveVtk(fun::FileNameSec(filename,numfile),n,pos,nfields,fields);

  //-Frees memory.
  delete[] pos;
  delete[] idp;
  delete[] dcel;
  delete[] vel;
  delete[] rhop;
  delete[] velm1;
  delete[] rhopm1;
  delete[] ace;
  delete[] type;
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,bool idp,bool vel,bool rhop,bool code)
{
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  if(Log->GetMpiRank()>=0)filename=string("p")+fun::IntStr(Log->GetMpiRank())+"_"+filename;
  filename=DirDataOut+filename;
  //-Allocates memory.
  const unsigned n=pfin-pini;
  ParticlesDataDown(n,pini,code,false);
  tfloat3 *pos=new tfloat3[n];
  for(unsigned p=0;p<n;p++)pos[p]=ToTFloat3(AuxPos[p]);
  byte *type=new byte[n];
  for(unsigned p=0;p<n;p++){
    const typecode cod=Code[p];
    byte tp=99;
    if(CODE_IsFixed(cod))tp=0;
    else if(CODE_IsMoving(cod))tp=1;
    else if(CODE_IsFloating(cod))tp=2;
    else if(CODE_IsFluid(cod))tp=3;
    if(CODE_IsNormal(cod))tp+=0;
    else if(CODE_IsPeriodic(cod))tp+=10;
    else if(CODE_IsOutIgnore(cod))tp+=20;
    else tp+=30;
    type[p]=tp;
  }
  JFormatFiles2::StScalarData fields[8];
  unsigned nfields=0;
  if(idp){  fields[nfields]=JFormatFiles2::DefineField("Idp"  ,JFormatFiles2::UInt32  ,1,Idp);     nfields++; }
  if(type){ fields[nfields]=JFormatFiles2::DefineField("Typex",JFormatFiles2::UChar8  ,1,type);    nfields++; }
  if(vel){  fields[nfields]=JFormatFiles2::DefineField("Vel"  ,JFormatFiles2::Float32 ,3,AuxVel);  nfields++; }
  if(rhop){ fields[nfields]=JFormatFiles2::DefineField("Rhop" ,JFormatFiles2::Float32 ,1,AuxRhop); nfields++; }
#ifdef CODE_SIZE4
  if(code){ fields[nfields]=JFormatFiles2::DefineField("Code" ,JFormatFiles2::UInt32  ,1,Code);    nfields++; }
#else
  if(code){ fields[nfields]=JFormatFiles2::DefineField("Code" ,JFormatFiles2::UShort16,1,Code);    nfields++; }
#endif
  //-Generates file.
  JFormatFiles2::SaveVtk(filename,n,pos,nfields,fields);
  //-Frees memory. 
  delete[] pos;
  delete[] type;
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,const float3 *posg,const byte *checkg,const unsigned *idpg,const float3 *velg,const float *rhopg){
  //-Allocates memory.
  const unsigned n=pfin-pini;
  tfloat3 *pos=new tfloat3[n];
  byte *check=NULL;
  unsigned *idp=NULL;
  tfloat3 *vel=NULL;
  float *rhop=NULL;
  if(checkg)check=new byte[n];
  if(idpg)idp=new unsigned[n];
  if(velg)vel=new tfloat3[n];
  if(rhopg)rhop=new float[n];
  //-Copies data from GPU.
  cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(checkg)cudaMemcpy(check,checkg+pini,sizeof(byte)*n,cudaMemcpyDeviceToHost);
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhopg)cudaMemcpy(rhop,rhopg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  //-Generates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,n,pos,check,idp,vel,rhop);
  //-Frees memory.
  delete[] pos;
  delete[] check;
  delete[] idp;
  delete[] vel;
  delete[] rhop;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const float3 *posg,const unsigned *idpg,const float3 *velg,const float *rhopg,const float *arg,const float3 *aceg,const float3 *vcorrg){
  const char met[]="DgSaveCsvParticlesGpu"; 
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned *idp=NULL;  if(idpg)idp=new unsigned[n];
  tfloat3 *pos=NULL;   if(posg)pos=new tfloat3[n];
  tfloat3 *vel=NULL;   if(velg)vel=new tfloat3[n];
  float *rhop=NULL;    if(rhopg)rhop=new float[n];
  float *ar=NULL;      if(arg)ar=new float[n];
  tfloat3 *ace=NULL;   if(aceg)ace=new tfloat3[n];
  tfloat3 *vcorr=NULL; if(vcorrg)vcorr=new tfloat3[n];
  //-Copies data from GPU.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhopg)cudaMemcpy(rhop,rhopg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(arg)cudaMemcpy(ar,arg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(aceg)cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(vcorrg)cudaMemcpy(vcorr,vcorrg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  CheckCudaError(met,"Failed copying data from GPU.");
  //-Generates CSV file.
  DgSaveCsvParticlesCpu(filename,numfile,0,n,head,pos,idp,vel,rhop,ar,ace,vcorr);
  //-Frees memory.
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rhop;
  delete[] ar;
  delete[] ace;
  delete[] vcorr;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu2(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const float3 *posg,const unsigned *idpg,const float3 *velg,const float *rhopg,const float4 *pospresg,const float4 *velrhopg){
  const char met[]="DgSaveCsvParticlesGpu2";
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned *idp=NULL;  if(idpg)idp=new unsigned[n];
  tfloat3 *pos=NULL;   if(posg)pos=new tfloat3[n];
  tfloat3 *vel=NULL;   if(velg)vel=new tfloat3[n];
  float *rhop=NULL;    if(rhopg)rhop=new float[n];
  tfloat4 *pospres=NULL;  if(pospresg)pospres=new tfloat4[n];
  tfloat4 *velrhop=NULL;  if(velrhopg)velrhop=new tfloat4[n];
  //-Copies data from GPU.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhopg)cudaMemcpy(rhop,rhopg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(pospresg)cudaMemcpy(pospres,pospresg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  if(velrhopg)cudaMemcpy(velrhop,velrhopg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  CheckCudaError(met,"Failed copying data from GPU.");
  //-Generates CSV file.
  DgSaveCsvParticles2(filename,numfile,0,n,head,pos,idp,vel,rhop,pospres,velrhop);
  //-Frees memory.
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rhop;
  delete[] pospres;
  delete[] velrhop;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticles2(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const tfloat3 *pos,const unsigned *idp,const tfloat3 *vel,const float *rhop,const tfloat4 *pospres,const tfloat4 *velrhop){
  const char met[]="DgSaveCsvParticlesCpu";
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Generates CSV file.
  ofstream pf;
  pf.open(filename.c_str());
  if(pf){
    if(!head.empty())pf << head << endl;
    pf << "Num";
    if(idp)pf << ";Idp";
    if(pos)pf << ";PosX;PosY;PosZ";
    if(vel)pf << ";VelX;VelY;VelZ";
    if(rhop)pf << ";Rhop";
    if(pospres)pf << ";Px;Py;Pz;Pp";
    if(velrhop)pf << ";Vx;Vy;Vz;Vr";
    pf << endl;
    const char fmt1[]="%f"; //="%24.16f";
    const char fmt3[]="%f;%f;%f"; //="%24.16f;%24.16f;%24.16f";
    for(unsigned p=pini;p<pfin;p++){
      pf << fun::UintStr(p-pini);
      if(idp)pf << ";" << fun::UintStr(idp[p]);
      if(pos)pf << ";" << fun::Float3Str(pos[p],fmt3);
      if(vel)pf << ";" << fun::Float3Str(vel[p],fmt3);
      if(rhop)pf << ";" << fun::FloatStr(rhop[p],fmt1);
      if(pospres)pf << ";" <<  fun::FloatStr(pospres[p].x,fmt1) << ";" << fun::FloatStr(pospres[p].y,fmt1) << ";" << fun::FloatStr(pospres[p].z,fmt1) << ";" << fun::FloatStr(pospres[p].w,fmt1);
      if(velrhop)pf << ";" <<  fun::FloatStr(velrhop[p].x,fmt1) << ";" << fun::FloatStr(velrhop[p].y,fmt1) << ";" << fun::FloatStr(velrhop[p].z,fmt1) << ";" << fun::FloatStr(velrhop[p].w,fmt1);
      pf << endl;
    }
    if(pf.fail())RunException(met,"Failed writing to file.",filename);
    pf.close();
  }
  else RunException(met,"File could not be opened.",filename);
}


