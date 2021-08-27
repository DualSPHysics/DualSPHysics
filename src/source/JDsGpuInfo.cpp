//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2021 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsGpuInfo.cpp \brief Implements the class \ref JDsGpuInfo.

#include "JDsGpuInfo.h"
#include "JLog2.h"
#include "FunctionsCuda.h"
#include "Functions.h"

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JDsGpuInfo::JDsGpuInfo(){
  ClassName="JDsGpuInfo";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsGpuInfo::~JDsGpuInfo(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsGpuInfo::Reset(){
  Ngpus=-1;
  GpuId=-1;
  AutoSelect=false;
  Name="";
  GlobalMem=0;
  SharedMem=0;
  Compute=0;
  Hardware="";
}

//==============================================================================
/// Shows information on available GPUs.
//==============================================================================
int JDsGpuInfo::ShowGpusInfo(JLog2 *log){
  //log->Print("[Available CUDA devices]");
  vector<string> gpuinfo;
  vector<fcuda::StGpuInfo> gpuprops;
  const int ngpus=fcuda::GetCudaDevicesInfo(&gpuinfo,&gpuprops);
  log->Print(gpuinfo);
  log->Print(" ");
  return(ngpus);
}

//==============================================================================
/// Returns number of available GPUs.
//==============================================================================
int JDsGpuInfo::GetNgpus(){
  if(Ngpus<0){
    int deviceCount=0;
    cudaGetDeviceCount(&deviceCount);
    Ngpus=deviceCount;
  }
  return(Ngpus);
}

//==============================================================================
/// Selects GPU device (-1 for automatic selection).
//==============================================================================
int JDsGpuInfo::SelectGpu(int gpuid){
  if(GetNgpus()){
    //-GPU selection.
    AutoSelect=(gpuid<0);
    if(AutoSelect){
      unsigned *ptr=NULL;
      cudaMalloc((void**)&ptr,sizeof(unsigned)*100);
      cudaFree(ptr);
    }
    else cudaSetDevice(gpuid);
    //-Get information on GPU selection.
    cudaDeviceProp devp;
    int dev;
    cudaGetDevice(&dev);
    cudaGetDeviceProperties(&devp,dev);
    GpuId=dev;
    if(!AutoSelect && gpuid!=GpuId)Run_Exceptioon("Requested GPU is not available.");
    Name=devp.name;
    GlobalMem=devp.totalGlobalMem;
    SharedMem=int(devp.sharedMemPerBlock);
    Compute=devp.major*10+devp.minor;
    Hardware=fun::PrintStr("Gpu_%d%s\"%s\"",GpuId,(AutoSelect? "?=": "="),Name.c_str());
  }
  else Run_Exceptioon("There are no available CUDA devices.");
  return(GpuId);
}

//==============================================================================
/// Shows main information on selected GPU.
//==============================================================================
void JDsGpuInfo::ShowSelectGpusInfo(JLog2 *log){
  log->Printf("Device %s: %d \"%s\"",(AutoSelect? "default": "selected"),GpuId,Name.c_str());
  log->Printf("Compute capability: %.1f",float(Compute)/10);
  log->Printf("Memory global: %d MB",int(GlobalMem/(1024*1024)));
  log->Printf("Memory shared: %u Bytes",SharedMem);
  log->Print(" ");
}

