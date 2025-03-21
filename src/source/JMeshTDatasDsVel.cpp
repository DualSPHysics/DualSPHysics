//HEAD_DSCODES
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

/// \file JMeshTDatasDsVel.cpp \brief Implements the classes JMeshTDatasDsVel.

#include "JMeshTDatasDsVel.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JMeshTDatasLoad.h"
#include "JMeshData.h"
#include "JDataArrays.h"

#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JSphGpu_InOut_iker.h"
  #include "JRelaxZone_ker.h"
#endif

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <algorithm>

using namespace std;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshTDatasDsVel
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMeshTDatasDsVel::JMeshTDatasDsVel(std::string appname,bool cpu)
  :JMeshTDatas(appname,false),Cpu(cpu)
{
  ClassName="JMeshTDatasDsVel";
#ifdef _WITHGPU
  Vel1Data1g=NULL;  Vel1Data2g=NULL;  Vel1DataTg=NULL;
  Vel3Data1g=NULL;  Vel3Data2g=NULL;  Vel3DataTg=NULL;
#endif
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JMeshTDatasDsVel::~JMeshTDatasDsVel(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JMeshTDatasDsVel::Reset(){
  JMeshTDatas::Reset();
  SetPos=TDouble3(0);
  InitialTime=0;
  LoopTbeginRq=DBL_MAX;
  LoopTmaxRq=DBL_MAX;
  VelMagnitude=false;
  VelDir=TFloat3(0);
  VelReverse=false;
  SetVelMul=TDouble3(0);
  SetVelAdd=TDouble3(0);
#ifdef _WITHGPU
  if(!Cpu)ResetGpu();
#endif
}

//==============================================================================
/// Configures for using.
//==============================================================================
void JMeshTDatasDsVel::ConfigVel(std::string file,tdouble3 setpos
  ,double initialtime,double looptmax,double looptbegin
  ,bool magnitude,tfloat3 veldir,bool reverse
  ,tdouble3 velmul,tdouble3 veladd)
{
  Reset();
  LoadFile(file,(magnitude? "VelDir": "Vel")
    ,DBL_MAX,DBL_MAX,looptmax,looptbegin
    ,TDouble3(DBL_MAX),TDouble3(DBL_MAX),0,setpos);
  if(reverse)ReverseData("all");
  InitialTime=initialtime;
  LoopTbeginRq=looptbegin;
  LoopTmaxRq=looptmax;
  SetPos=setpos;
  IntpConfigFreePos();
#ifdef _WITHGPU
  if(!Cpu)AllocMemoryGpu();
#endif
  VelMagnitude=magnitude;
  VelDir=fgeo::VecUnitary(veldir);
  if(VelMagnitude && VelDir==TFloat3(0))Run_Exceptioon("Direction vector for velocity is invalid.");
  VelReverse=reverse;
  SetVelMul=velmul;
  if(SetVelMul!=TDouble3(1)){
    if(SetVelMul.x==SetVelMul.y && SetVelMul.x==SetVelMul.z)SetMulData("vel",SetVelMul.x);
    else{
      if(SetVelMul.x!=1)SetMulData("vel.x",SetVelMul.x);
      if(SetVelMul.y!=1)SetMulData("vel.y",SetVelMul.y);
      if(SetVelMul.z!=1)SetMulData("vel.z",SetVelMul.z);
    }
  }
  SetVelAdd=veladd;
  if(SetVelAdd!=TDouble3(0)){
    if(SetVelAdd.x==SetVelAdd.y && SetVelAdd.x==SetVelAdd.z)SetAddData("vel",SetVelAdd.x);
    else{
      if(SetVelAdd.x!=0)SetAddData("vel.x",SetVelAdd.x);
      if(SetVelAdd.y!=0)SetAddData("vel.y",SetVelAdd.y);
      if(SetVelAdd.z!=0)SetAddData("vel.z",SetVelAdd.z);
    }
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Frees memory on GPU.
//==============================================================================
void JMeshTDatasDsVel::ResetGpu(){
  GpuCtime1=GpuCtime2=UINT_MAX;
  FreeMemoryGpu();
}

//==============================================================================
/// Frees memory on GPU.
//==============================================================================
void JMeshTDatasDsVel::FreeMemoryGpu(){
   if(Vel1Data1g)cudaFree(Vel1Data1g);  Vel1Data1g=NULL;
   if(Vel1Data2g)cudaFree(Vel1Data2g);  Vel1Data2g=NULL;
   if(Vel1DataTg)cudaFree(Vel1DataTg);  Vel1DataTg=NULL;
   if(Vel3Data1g)cudaFree(Vel3Data1g);  Vel3Data1g=NULL;
   if(Vel3Data2g)cudaFree(Vel3Data2g);  Vel3Data2g=NULL;
   if(Vel3DataTg)cudaFree(Vel3DataTg);  Vel3DataTg=NULL;
}
//==============================================================================
/// Allocates memory on GPU.
//==============================================================================
void JMeshTDatasDsVel::AllocMemoryGpu(){
  if(FrType==TypeFloat){
    fcuda::Malloc(&Vel1Data1g,FrSize);
    fcuda::Malloc(&Vel1Data2g,FrSize);
    fcuda::Malloc(&Vel1DataTg,FrSize);
  }
  else if(FrType==TypeFloat3){
    fcuda::Malloc(&Vel3Data1g,FrSize);
    fcuda::Malloc(&Vel3Data2g,FrSize);
    fcuda::Malloc(&Vel3DataTg,FrSize);
  }
  else Run_Exceptioon("Type of data is not supported.");
}

//==============================================================================
/// Computes interpolation on requeseted time when it is necessary.
//==============================================================================
void JMeshTDatasDsVel::IntpComputeTimeGpu(double timestep){
  FindTime(timestep);
  if(FrCtime!=Position || FrTimeFactor!=TimeFactor){
    const float tfactor=float(TimeFactor);
    const bool flt1=(FrType==TypeFloat);
    //printf("\nIntpComputeTimeGpu> t:%f  tf:%f  positions:%u-%u  cts:%u-%u\n",timestep,tfactor,Position,PositionNext,GpuCtime1,GpuCtime2);
    const size_t size=(flt1? sizeof(float): sizeof(float3))*FrSize;
    //-Updates data on GPU memory.
    if(GpuCtime1!=Position || GpuCtime2!=PositionNext){
      //-Swap data when it is possible.
      if(GpuCtime2==Position || GpuCtime1==PositionNext){
        swap(GpuCtime1,GpuCtime2);
        swap(Vel1Data1g,Vel1Data2g);
        swap(Vel3Data1g,Vel3Data2g);
      }
      if(GpuCtime1!=Position){
        const void* data=Datas[Position]->GetArrays()->GetArrayCte(0).ptr;
        if(flt1)cudaMemcpy(Vel1Data1g,data,size,cudaMemcpyHostToDevice);
        else    cudaMemcpy(Vel3Data1g,data,size,cudaMemcpyHostToDevice);
        GpuCtime1=Position;
      }
      if(GpuCtime2!=PositionNext){
        if(GpuCtime2==GpuCtime1){
          if(flt1)cudaMemcpy(Vel1Data2g,Vel1Data1g,size,cudaMemcpyDeviceToDevice);
          else    cudaMemcpy(Vel3Data2g,Vel3Data1g,size,cudaMemcpyDeviceToDevice);
        }
        else{
          const void* data=Datas[PositionNext]->GetArrays()->GetArrayCte(0).ptr;
          if(flt1)cudaMemcpy(Vel1Data2g,data,size,cudaMemcpyHostToDevice);
          else    cudaMemcpy(Vel3Data2g,data,size,cudaMemcpyHostToDevice);
        }
        GpuCtime2=PositionNext;
      }
    }
    //-Interpolate data on timestep.
    if(Position==PositionNext || TimeFactor==0){
      if(flt1)cudaMemcpy(Vel1DataTg,Vel1Data1g,size,cudaMemcpyDeviceToDevice);
      else    cudaMemcpy(Vel3DataTg,Vel3Data1g,size,cudaMemcpyDeviceToDevice);
    }
    else if(TimeFactor>=1){
      if(flt1)cudaMemcpy(Vel1DataTg,Vel1Data2g,size,cudaMemcpyDeviceToDevice);
      else    cudaMemcpy(Vel3DataTg,Vel3Data2g,size,cudaMemcpyDeviceToDevice);
    }
    else{
      if(flt1)cusphinout::InOutInterpolateDataTime(FrSize,tfactor,Vel1Data1g,Vel1Data2g,Vel1DataTg);
      else    cusphinout::InOutInterpolateDataTime(FrSize,tfactor,Vel3Data1g,Vel3Data2g,Vel3DataTg);
    }
    FrCtime=Position;
    FrTimeFactor=TimeFactor;
  }
  FrTimeStep=timestep;
}
#endif

//==============================================================================
/// Interpolate velocity in time and position of selected partiles in a list.
//==============================================================================
void JMeshTDatasDsVel::InterpolateInOutVelCpu(double t,unsigned izone,unsigned np,const int* plist
  ,const tdouble3* pos,const typecode* code,const unsigned* idp,tfloat4* velrhop,float velcorr1)
{
  if(InConfigured!=IntpVar || (FrType!=TypeFloat && FrType!=TypeFloat3))
    Run_Exceptioon("Free interplation configuration is missing or data type does not match.");
  const bool flt1=(FrType==TypeFloat);
  //-Select previous and following times with data.
  if(FrTimeStep!=t){
    if(flt1)IntpComputeTime_f(t);
    else    IntpComputeTime_f3(t);
  }
  //-Interpolate in positions.
  const float*   tdat1=(float*)  (FrType==TypeFloat?  FrData: NULL);
  const tfloat3* tdat3=(tfloat3*)(FrType==TypeFloat3? FrData: NULL);
  const tfloat3 velcorr3=VelDir*velcorr1;
  const int n=int(np);
  switch(FrMode){
    case 0:{
      const tfloat3 vel=(flt1? VelDir*(tdat1[0]-velcorr1): tdat3[0]-velcorr3);
      #ifdef OMP_USE
        #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
      #endif
      for(int cp=0;cp<n;cp++){
        const unsigned p=plist[cp];
        if(izone==CODE_GetIzoneFluidInout(code[p]))velrhop[p]=TFloat4(vel.x,vel.y,vel.z,velrhop[p].w);
      }
    }break;
    case 1:{
      if(flt1){
        #ifdef OMP_USE
          #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
        #endif
        for(int cp=0;cp<n;cp++){
          const unsigned p=plist[cp];
          if(izone==CODE_GetIzoneFluidInout(code[p])){
            const float vel=IntpComputeVarFr1_f(pos[p])-velcorr1;
            velrhop[p]=TFloat4(VelDir.x*vel,VelDir.y*vel,VelDir.z*vel,velrhop[p].w);
          }
        }
      }
      else{
        #ifdef OMP_USE
          #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
        #endif
        for(int cp=0;cp<n;cp++){
          const unsigned p=plist[cp];
          if(izone==CODE_GetIzoneFluidInout(code[p])){
            const tfloat3 vel=IntpComputeVarFr1_f3(pos[p])-velcorr3;
            velrhop[p]=TFloat4(vel.x,vel.y,vel.z,velrhop[p].w);
          }
        }
      }
    }break;
    case 2:{
      if(flt1){
        #ifdef OMP_USE
          #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
        #endif
        for(int cp=0;cp<n;cp++){
          const unsigned p=plist[cp];
          if(izone==CODE_GetIzoneFluidInout(code[p])){
            const float vel=IntpComputeVarFr2_f(pos[p])-velcorr1;  
            velrhop[p]=TFloat4(VelDir.x*vel,VelDir.y*vel,VelDir.z*vel,velrhop[p].w);
          }
        }
      }
      else{
        #ifdef OMP_USE
          #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
        #endif
        for(int cp=0;cp<n;cp++){
          const unsigned p=plist[cp];
          if(izone==CODE_GetIzoneFluidInout(code[p])){
            const tfloat3 vel=IntpComputeVarFr2_f3(pos[p])-velcorr3;  
            velrhop[p]=TFloat4(vel.x,vel.y,vel.z,velrhop[p].w);
          }
        }
      }
    }break;
    case 3:{
      if(flt1){
        #ifdef OMP_USE
          #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
        #endif
        for(int cp=0;cp<n;cp++){
          const unsigned p=plist[cp];
          if(izone==CODE_GetIzoneFluidInout(code[p])){
            const float vel=IntpComputeVarFr3_f(pos[p])-velcorr1;  
            velrhop[p]=TFloat4(VelDir.x*vel,VelDir.y*vel,VelDir.z*vel,velrhop[p].w);
          }
        }
      }
      else{
        #ifdef OMP_USE
          #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
        #endif
        for(int cp=0;cp<n;cp++){
          const unsigned p=plist[cp];
          if(izone==CODE_GetIzoneFluidInout(code[p])){
            const tfloat3 vel=IntpComputeVarFr3_f3(pos[p])-velcorr3;  
            velrhop[p]=TFloat4(vel.x,vel.y,vel.z,velrhop[p].w);
          }
        }
      }
    }break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
}


//==============================================================================
/// Interpolate velocity in time and position of selected partiles.
//==============================================================================
void JMeshTDatasDsVel::InterpolateRzVelCpu(double t,byte zoneid,unsigned np
  ,unsigned pini,const byte* rzid,const float* rzfactor,const tdouble3* pos
  ,tfloat4* velrhop)
{
  if(InConfigured!=IntpVar || (FrType!=TypeFloat && FrType!=TypeFloat3))Run_Exceptioon("Free interplation configuration is missing or data type does not match.");
  const bool flt1=(FrType==TypeFloat);
  //-Select previous and following times with data.
  if(FrTimeStep!=t){
    if(flt1)IntpComputeTime_f(t);
    else    IntpComputeTime_f3(t);
  }
  //-Interpolate in positions.
  const float*   tdat1=(float*)  (FrType==TypeFloat?  FrData: NULL);
  const tfloat3* tdat3=(tfloat3*)(FrType==TypeFloat3? FrData: NULL);
  const int p0=int(pini);
  const int pfin=p0+int(np);
  switch(FrMode){
    case 0:{
      const tfloat3 vt=(flt1? VelDir*(tdat1[0]): tdat3[0]);
      #ifdef OMP_USE
        #pragma omp parallel for schedule (guided)
      #endif
      for(int p=p0;p<pfin;p++){
        if(rzid[p]==zoneid){
          const float f=rzfactor[p],ff=(1.f-f);
          const tfloat4 v=velrhop[p];
          velrhop[p]=TFloat4(f*vt.x+ff*v.x, f*vt.y+ff*v.y, f*vt.z+ff*v.z, v.w);
        }
      }
    }break;
    case 1:{
      if(flt1){
        #ifdef OMP_USE
          #pragma omp parallel for schedule (guided)
        #endif
        for(int p=p0;p<pfin;p++){
          if(rzid[p]==zoneid){
            const float f=rzfactor[p],ff=(1.f-f);
            const tfloat4 v=velrhop[p];
            const tfloat3 vt=VelDir*IntpComputeVarFr1_f(pos[p]);
            velrhop[p]=TFloat4(f*vt.x+ff*v.x, f*vt.y+ff*v.y, f*vt.z+ff*v.z, v.w);
          }
        }
      }
      else{
        #ifdef OMP_USE
          #pragma omp parallel for schedule (guided)
        #endif
        for(int p=p0;p<pfin;p++){
          if(rzid[p]==zoneid){
            const float f=rzfactor[p],ff=(1.f-f);
            const tfloat4 v=velrhop[p];
            const tfloat3 vt=IntpComputeVarFr1_f3(pos[p]);
            velrhop[p]=TFloat4(f*vt.x+ff*v.x, f*vt.y+ff*v.y, f*vt.z+ff*v.z, v.w);
          }
        }
      }
    }break;
    case 2:{
      if(flt1){
        #ifdef OMP_USE
          #pragma omp parallel for schedule (guided)
        #endif
        for(int p=p0;p<pfin;p++){
          if(rzid[p]==zoneid){
            const float f=rzfactor[p],ff=(1.f-f);
            const tfloat4 v=velrhop[p];
            const tfloat3 vt=VelDir*IntpComputeVarFr2_f(pos[p]);
            velrhop[p]=TFloat4(f*vt.x+ff*v.x, f*vt.y+ff*v.y, f*vt.z+ff*v.z, v.w);
          }
        }
      }
      else{
        #ifdef OMP_USE
          #pragma omp parallel for schedule (guided)
        #endif
        for(int p=p0;p<pfin;p++){
          if(rzid[p]==zoneid){
            const float f=rzfactor[p],ff=(1.f-f);
            const tfloat4 v=velrhop[p];
            const tfloat3 vt=IntpComputeVarFr2_f3(pos[p]);
            velrhop[p]=TFloat4(f*vt.x+ff*v.x, f*vt.y+ff*v.y, f*vt.z+ff*v.z, v.w);
          }
        }
      }
    }break;
    case 3:{
      if(flt1){
        #ifdef OMP_USE
          #pragma omp parallel for schedule (guided)
        #endif
        for(int p=p0;p<pfin;p++){
          if(rzid[p]==zoneid){
            const float f=rzfactor[p],ff=(1.f-f);
            const tfloat4 v=velrhop[p];
            const tfloat3 vt=VelDir*IntpComputeVarFr3_f(pos[p]);
            velrhop[p]=TFloat4(f*vt.x+ff*v.x, f*vt.y+ff*v.y, f*vt.z+ff*v.z, v.w);
          }
        }
      }
      else{
        #ifdef OMP_USE
          #pragma omp parallel for schedule (guided)
        #endif
        for(int p=p0;p<pfin;p++){
          if(rzid[p]==zoneid){
            const float f=rzfactor[p],ff=(1.f-f);
            const tfloat4 v=velrhop[p];
            const tfloat3 vt=IntpComputeVarFr3_f3(pos[p]);
            velrhop[p]=TFloat4(f*vt.x+ff*v.x, f*vt.y+ff*v.y, f*vt.z+ff*v.z, v.w);
          }
        }
      }
    }break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Interpolate velocity in time and position of selected partiles in a list.
//==============================================================================
void JMeshTDatasDsVel::InterpolateInOutVelGpu(double t,unsigned izone
  ,unsigned np,const int* plistg,const double2* posxyg,const double* poszg
  ,const typecode* codeg,const unsigned* idpg,float4* velrhopg,float velcorr)
{
  t+=InitialTime;
  //printf("t:%f  InitialTime:%f\n",t,InitialTime);
  if(InConfigured!=IntpVar || (FrType!=TypeFloat && FrType!=TypeFloat3))Run_Exceptioon("Free interplation configuration is missing or data type does not match.");
  const bool flt1=(FrType==TypeFloat);
  //-Select previous and following times with data.
  if(FrTimeStep!=t)IntpComputeTimeGpu(t);
  //-Interpolate in positions.
  const byte iz=byte(izone);
  const float*  data1=(flt1? Vel1DataTg: NULL);
  const float3* data3=(flt1? NULL: Vel3DataTg);
  switch(FrMode){
    case 0:{
      cusphinout::InOutIntpVelFr0(iz,velcorr,VelDir,data1,data3,np,plistg,codeg,velrhopg);
    }break;
    case 1:{
      cusphinout::InOutIntpVelFr1(iz,velcorr,VelDir,data1,data3
        ,FrPla1,FrCmax1,np,plistg,codeg,posxyg,poszg,velrhopg);
    }break;
    case 2:{
      cusphinout::InOutIntpVelFr2(iz,velcorr,VelDir,data1,data3,FrNum1
        ,FrPla1,FrCmax1,FrPla2,FrCmax2,np,plistg,codeg,posxyg,poszg,velrhopg);
    }break;
    case 3:{
      cusphinout::InOutIntpVelFr3(iz,velcorr,VelDir,data1,data3,Npt1,Npt12
        ,FrPla1,FrCmax1,FrPla2,FrCmax2,FrPla3,FrCmax3
        ,np,plistg,codeg,posxyg,poszg,velrhopg);
    }break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
}
#endif

//==============================================================================
/// Prepares data for two times, set time0 and time1 and return pointers (cpu or 
/// gpu) for the modification of velocity data.
//==============================================================================
StRnVelData JMeshTDatasDsVel::RnGetVelPtr(double time0,double time1){ 
  if(InConfigured!=IntpVar || (FrType!=TypeFloat && FrType!=TypeFloat3))
    Run_Exceptioon("Operation is invalid because only Velocity magnintude mode is supported for now.");
  if(SizeTimes<2)Run_Exceptioon("The mesh-data must be defined with at least 2 data instants.");
  const bool vf3=(FrType==TypeFloat3);
  SizeTimes=2;
  FrTimeStep=DBL_MAX;
  FrCtime=Position=UINT_MAX;
  Times[0]=time0;
  Times[1]=time1;
  float* ptr0=NULL;
  float* ptr1=NULL;
  if(Cpu){
    ptr0=(float*)Datas[0]->GetArrays()->GetArrayCte(0).ptr;
    ptr1=(float*)Datas[1]->GetArrays()->GetArrayCte(0).ptr;
  }
  else{
   #ifdef _WITHGPU
    if(vf3){
      ptr0=(float*)Vel3Data1g;
      ptr1=(float*)Vel3Data2g;
    }
    else{
      ptr0=Vel1Data1g;
      ptr1=Vel1Data2g;
    }
   #endif
  }
  StRnVelData ret{FrSize,ptr0,ptr1,!Cpu,vf3};
  return(ret);
}


}
