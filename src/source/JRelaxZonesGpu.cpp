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

/// \file JRelaxZonesGpu.cpp \brief Implements the class \ref JRelaxZonesGpu.

#include "JRelaxZonesGpu.h"
#ifdef _WITHGPU
  #include "JRelaxZone_ker.h"
  #include <cuda_runtime_api.h>
#endif


using namespace std;

//##############################################################################
//# JRelaxZoneRegularGpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JRelaxZoneRegularGpu::JRelaxZoneRegularGpu(){
  ClassName="JRelaxZoneRegularGpu";
}

//==============================================================================
/// Sets velocity of fluid to generate regular waves.
//==============================================================================
void JRelaxZoneRegularGpu::SetFluidVel(unsigned n,unsigned pini,bool order2,bool subdrift
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double kl,double sinhkld
  ,double wpf,double cta,double depth,double framp
  ,double ct2,double sinhkld4
  ,double ctd,double ctd2,unsigned fluidbeginidp
  ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop)
{
  #ifdef _WITHGPU
    curelaxzone::SetFluidVel(n,pini,order2,subdrift,centerx,widthhalf,coeffx,coeffz,falpha,fbeta
      ,fsub,fdiv,timewave,swl,kl,sinhkld,wpf,cta,depth,framp,ct2,sinhkld4,ctd,ctd2,fluidbeginidp
      ,(const double2*)posxy,posz,idp,(float4*)velrhop);
  #endif
}


//##############################################################################
//# JRelaxZoneSpectrumGpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JRelaxZoneSpectrumGpu::JRelaxZoneSpectrumGpu(){
  ClassName="JRelaxZoneSpectrumGpu";
  WavesOnGpu=false; 
  WaveKlg=NULL; WaveAmpg=NULL; WaveFangg=NULL; WavePhaseg=NULL;
  MemGpuFixed=0; 
}

//==============================================================================
/// Destructor.
//==============================================================================
JRelaxZoneSpectrumGpu::~JRelaxZoneSpectrumGpu(){
  DestructorActive=true; 
  FreeMemoryGpu(); 
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JRelaxZoneSpectrumGpu::FreeMemoryGpu(){
  MemGpuFixed=0;
  WavesOnGpu=false;
  #ifdef _WITHGPU
    if(WaveKlg)   cudaFree(WaveKlg);    WaveKlg=NULL;
    if(WaveAmpg)  cudaFree(WaveAmpg);   WaveAmpg=NULL;
    if(WaveFangg) cudaFree(WaveFangg);  WaveFangg=NULL;
    if(WavePhaseg)cudaFree(WavePhaseg); WavePhaseg=NULL;
  #endif 
}

//==============================================================================
/// Allocates GPU memory.
//==============================================================================
void JRelaxZoneSpectrumGpu::AllocMemoryGpu(unsigned wavecount){
  FreeMemoryGpu();
  #ifdef _WITHGPU
    const size_t m=sizeof(double)*wavecount;
    cudaMalloc((void**)&WaveKlg,m);
    cudaMalloc((void**)&WaveAmpg,m);
    cudaMalloc((void**)&WaveFangg,m);
    cudaMalloc((void**)&WavePhaseg,m);
    MemGpuFixed=m*4;
  #endif 
}

//==============================================================================
// Copy wave data on GPU memory.
//==============================================================================
void JRelaxZoneSpectrumGpu::PrepareWaveDataGpu(unsigned wavecount
  ,const double *kl,const double *amp,const double *fang,const double *phase)
{
  AllocMemoryGpu(wavecount);
  #ifdef _WITHGPU
    cudaMemcpy(WaveKlg   ,kl   ,sizeof(double)*wavecount,cudaMemcpyHostToDevice);
    cudaMemcpy(WaveAmpg  ,amp  ,sizeof(double)*wavecount,cudaMemcpyHostToDevice);
    cudaMemcpy(WaveFangg ,fang ,sizeof(double)*wavecount,cudaMemcpyHostToDevice);
    cudaMemcpy(WavePhaseg,phase,sizeof(double)*wavecount,cudaMemcpyHostToDevice);
    WavesOnGpu=true;
  #endif 
}

//==============================================================================
/// Sets velocity of fluid to generate irregular waves.
//==============================================================================
void JRelaxZoneSpectrumGpu::SetFluidVelSpectrumSub(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double depth,double framp,unsigned wavecount
  ,unsigned fluidbeginidp
  ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop
  ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2)
{
  #ifdef _WITHGPU
    curelaxzone::SetFluidVelSpectrumSub(n,pini,centerx,widthhalf,coeffx,coeffz
      ,falpha,fbeta,fsub,fdiv
      ,timewave,swl,depth,framp
      ,wavecount,WaveKlg,WaveAmpg,WaveFangg,WavePhaseg,fluidbeginidp
      ,(const double2*)posxy,posz,idp,(float4*)velrhop
      ,subdrift,fun,ctd,ctd2,ctd_2,ctd2_2);
  #endif
}


//##############################################################################
//# JRelaxZonesExternalGpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JRelaxZonesExternalGpu::JRelaxZonesExternalGpu(){
  ClassName="JRelaxZonesExternalGpu";
  GpuReady=false; 
  VelXg=NULL; VelZg=NULL;
  MemGpuFixed=0; 
}

//==============================================================================
/// Destructor.
//==============================================================================
JRelaxZonesExternalGpu::~JRelaxZonesExternalGpu(){
  DestructorActive=true; 
  FreeMemoryGpu();
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JRelaxZonesExternalGpu::FreeMemoryGpu(){
  MemGpuFixed=0;
  GpuReady=false;
  #ifdef _WITHGPU
    if(VelXg)cudaFree(VelXg); VelXg=NULL;
    if(VelZg)cudaFree(VelZg); VelZg=NULL;
#endif 
}

//==============================================================================
/// Allocates GPU memory for velocities.
//==============================================================================
void JRelaxZonesExternalGpu::AllocMemoryGpu(unsigned size,bool loadvelz){
  FreeMemoryGpu();
  #ifdef _WITHGPU
    const size_t m=sizeof(double)*size;
    cudaMalloc((void**)&VelXg,m);
    if(loadvelz)cudaMalloc((void**)&VelZg,m);
    GpuReady=true;
    MemGpuFixed=m*2;
  #endif 
}

//==============================================================================
// Copy wave data on GPU memory.
//==============================================================================
void JRelaxZonesExternalGpu::PrepareDataGpu(unsigned size,bool loadvelz,const double *velx,const double *velz){
  if(!GpuReady)AllocMemoryGpu(size,loadvelz);
  #ifdef _WITHGPU
    cudaMemcpy(VelXg,velx,sizeof(double)*size,cudaMemcpyHostToDevice);
    if(loadvelz)cudaMemcpy(VelZg,velz,sizeof(double)*size,cudaMemcpyHostToDevice);
  #endif 
}

//==============================================================================
/// Sets velocity of fluid to generate irregular waves.
//==============================================================================
void JRelaxZonesExternalGpu::SetFluidVelExternal(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double pxmin,double pymin,double pzmin
  ,double dpx,double dpy,double dpz
  ,unsigned npx1,unsigned npy1,unsigned npz1
  ,unsigned fluidbeginidp
  ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop
  ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2,double bottom)
{
  #ifdef _WITHGPU
    curelaxzone::SetFluidVelExternal(n,pini,centerx,widthhalf,coeffx,coeffz
      ,falpha,fbeta,fsub,fdiv
      ,pxmin,pymin,pzmin
      ,dpx,dpy,dpz
      ,npx1,npy1,npz1
      ,VelXg,VelZg
      ,fluidbeginidp
      ,(const double2*)posxy,posz,idp,(float4*)velrhop
      ,subdrift,fun,ctd,ctd2,ctd_2,ctd2_2,bottom);
  #endif
}


//##############################################################################
//# JRelaxZoneUniformGpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JRelaxZoneUniformGpu::JRelaxZoneUniformGpu(){
  ClassName="JRelaxZoneUniformGpu";
}

//==============================================================================
/// Sets velocity of fluid to generate regular waves.
//==============================================================================
void JRelaxZoneUniformGpu::SetFluidVelUniform(unsigned n,unsigned pini
  ,const tfloat3 &vt,const tfloat4 &cenpla
  ,const tfloat4 &dompla1,const tfloat4 &dompla2,const tfloat4 &dompla3
  ,const float domsize1,const float domsize2,const float domsize3,float widthhalf
  ,float coeff,double falpha,double fbeta,double fsub,double fdiv,unsigned fluidbeginidp
  ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop)
{
  #ifdef _WITHGPU
    curelaxzone::SetFluidVelUniform(n,pini
      ,vt,cenpla,dompla1,dompla2,dompla3
      ,domsize1,domsize2,domsize3,widthhalf
      ,coeff,falpha,fbeta,fsub,fdiv,fluidbeginidp
      ,(const double2*)posxy,posz,idp,(float4*)velrhop);
  #endif
}














