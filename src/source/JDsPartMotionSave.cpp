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

/// \file JDsPartMotionSave.cpp \brief Implements the class JDsPartMotionSave.

#include "JDsPartMotionSave.h"
#include "JPartMotRefBi4Save.h"
#include "JComputeMotionRef.h"
#include "Functions.h"
#include "JAppInfo.h"
#include "JLog2.h"

#ifdef _WITHGPU
#include "FunctionsCuda.h"
#include "JSphGpu_ker.h"
#endif

#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JDsPartMotionSave
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsPartMotionSave::JDsPartMotionSave(bool cpu,std::string appname,std::string dir
  ,unsigned casenfixed,unsigned casenmoving,unsigned casenfloat
  ,word mkboundfirst,unsigned mkcount,double timeout,double timeout2)
  :Log(AppInfo.LogPtr()),Cpu(cpu),AppName(appname),Dir(dir)
  ,CaseNfixed(casenfixed),CaseNmoving(casenmoving),CaseNfloat(casenfloat)
  ,MkBoundFirst(mkboundfirst),MkCount(mkcount),PsCount(mkcount*3) 
  ,TimeOut(timeout),TimeOut2(timeout2)
{
  ClassName="JDsPartMotionSave";
  MkMotionData=new StMkMotionData[MkCount];

  IdpRef=NULL;  PosRef=NULL;
  #ifdef _WITHGPU
    IdpRefg=NULL; PosRefg=NULL;
  #endif

  MotData1=NULL;
  MotData2=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsPartMotionSave::~JDsPartMotionSave(){
  DestructorActive=true;
  Reset();
  delete[] MkMotionData; MkMotionData=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsPartMotionSave::Reset(){
  MkMovingCount=MkFloatCount=0;
  IdpBegin=0;
  IdpCount=0;
  LastStep=UINT_MAX;
  LastTimestep=-DBL_MAX;
  //-Auxiliary memory.
  delete[] IdpRef; IdpRef=NULL;
  delete[] PosRef; PosRef=NULL;
  //-Free GPU memory.
  #ifdef _WITHGPU
    ResetDataGpu();
  #endif
  //-Free JPartMotRefBi4Save objects.
  delete MotData1;  MotData1=NULL;
  delete MotData2;  MotData2=NULL;

  UseUnits=false;

}

//==============================================================================
/// Add configuration of moving or floating mk block.
//==============================================================================
void JDsPartMotionSave::ConfigAddMk(bool floating,word mktype,unsigned begin
  ,unsigned mknp)
{
  unsigned cmk=MkMovingCount+MkFloatCount;
  if(cmk>=MkCount)Run_Exceptioon("Number of mk blocks invalid.");
  if(!floating && MkFloatCount)Run_Exceptioon("Adding moving after floating blocks.");
  StMkMotionData& v=MkMotionData[cmk];
  v.mkbound=mktype;
  v.begin=begin;
  v.np=mknp;
  v.nid=0;
  if(floating)MkFloatCount++;
  else        MkMovingCount++;
}

//==============================================================================
/// Configures reference position data of mk blocks.
//==============================================================================
void JDsPartMotionSave::ConfigMotionRefs(unsigned np,const tdouble3* pos
  ,const unsigned* idp)
{
  //-Checks configuration.
  if(MkMovingCount+MkFloatCount!=MkCount || !MkCount)Run_Exceptioon("Number of defined mk blocks is invalid.");
  if(MkMotionData[0].begin!=CaseNfixed)Run_Exceptioon("Begin of first MK block is invalid.");
  if(MkFloatCount && MkMotionData[MkMovingCount].begin!=CaseNfixed+CaseNmoving)Run_Exceptioon("Begin of first floating block is invalid.");
  {
    unsigned nmv=0,nft=0;
    for(unsigned cmk=0;cmk<MkCount;cmk++){
      if(cmk<MkMovingCount)nmv+=MkMotionData[cmk].np;
      else                 nft+=MkMotionData[cmk].np;
    }
    if(nmv!=CaseNmoving || nft!=CaseNfloat)Run_Exceptioon("Number of moving o floating particles is invalid.");
  }
  //-Confgures range of selected idp values.
  IdpBegin=CaseNfixed;
  IdpCount=CaseNmoving+CaseNfloat;
  //-Find reference points to compute motion.
  JComputeMotionRef computemotref;
  computemotref.AddMkBlocks(MkBoundFirst,MkCount,MkMotionData);
  computemotref.ComputeRefPoints(IdpBegin,0,IdpCount,np,idp,pos,NULL,NULL);
  computemotref.GetMotionRef(MkCount,MkMotionData);
  computemotref.Reset();
  //-Allocates auxiliary memory.
  IdpRef=new unsigned[PsCount];
  PosRef=new tdouble3[PsCount];
  //-Loads initial IdpRef[] and PosRef[].
  for(unsigned cmk=0;cmk<MkCount;cmk++){
    const StMkMotionData& v=MkMotionData[cmk];
    //Log->Printf("[%d] mkbound:%d  nid:%u",cmk,v.mkbound,v.nid);
    for(unsigned ci=0;ci<3;ci++){
      //IdpRef[cmk*3+ci]=(ci<v.nid? v.id[ci]-IdpBegin: UINT_MAX);
      IdpRef[cmk*3+ci]=(ci<v.nid? v.id[ci]: UINT_MAX);
      PosRef[cmk*3+ci]=(ci<v.nid? v.ps[ci]: TDouble3(0));
      //Log->Printf("IdpRef[%u]=%u",cmk*3+ci,IdpRef[cmk*3+ci]);
    }
  }
  //-Creates JPartMotRefBi4Save objects for output data.
  MotData1=new JPartMotRefBi4Save(AppName,Dir);
  MotData1->Config(true,TimeOut,MkBoundFirst,MkMovingCount,MkFloatCount,MkMotionData);
  if(TimeOut2>=0){
    MotData2=new JPartMotRefBi4Save(AppName,Dir);
    MotData2->Config(false,TimeOut2,MkBoundFirst,MkMovingCount,MkFloatCount,MkMotionData);
  }
  //-Configure data for GPU execution.
  #ifdef _WITHGPU
  if(!Cpu)ConfigDataGpu();
  #endif
}

//==============================================================================
/// Initial files recording with info from Data.
//==============================================================================
void JDsPartMotionSave::SaveInitial(){
  if(MotData1)MotData1->SaveInitial();
  if(MotData2)MotData2->SaveInitial();
}

//==============================================================================
/// Returns file name BI4.
//==============================================================================
std::string JDsPartMotionSave::GetFileName(bool mainfile)const{
  string fname;
  if(mainfile  && MotData1)fname=JPartMotRefBi4Save::GetFileNameDef(true);
  if(!mainfile && MotData2)fname=JPartMotRefBi4Save::GetFileNameDef(false);
  return(fname);
}

//==============================================================================
/// Load current reference position data from particle data on CPU.
//==============================================================================
void JDsPartMotionSave::LoadPosRefCpu(double timestep,unsigned step,unsigned np
  ,const tdouble3* pos,const unsigned* ridpmot)
{
  //Log->Printf("===> CaseNfixed:%u  MkMotionData[0].begin:%u",CaseNfixed,MkMotionData[0].begin);
  //Log->Printf("===-> LoadPosRefCpu: %u",step);
  if(LastStep!=step || LastTimestep!=timestep){
    if(!pos || !ridpmot)Run_Exceptioon("No data pointer valid.");
    for(unsigned cp=0;cp<PsCount;cp++){
      const unsigned iref=IdpRef[cp]-CaseNfixed;
      //Log->Printf("  IdpRef[%u]:%u  iref:%u  IdpCount:%u",cp,IdpRef[cp],iref,IdpCount);
      if(iref>=IdpCount)Run_Exceptioon("Error: IdpRef value is invalid.");
      const unsigned p=ridpmot[iref];
      if(p>=np)Run_Exceptioon("Error: Index of particle is invalid.");
      PosRef[cp]=pos[p];
    }
    LastStep=step;
    LastTimestep=timestep;
  }
}

//==============================================================================
/// Saves MAIN current reference position data in BI4 file.
//==============================================================================
void JDsPartMotionSave::SaveDataMainCpu(int cpart,double timestep,unsigned step
  ,unsigned np,const tdouble3* pos,const unsigned* ridpmot)
{
  //-Obtains current reference positions.
  LoadPosRefCpu(timestep,step,np,pos,ridpmot);
  //-Saves data in BI4 file.
  MotData1->SavePart(cpart,LastTimestep,LastStep,PsCount,PosRef);
}

//==============================================================================
/// Saves EXTRA current reference position data in BI4 file.
//==============================================================================
void JDsPartMotionSave::AddDataExtraCpu(int cpart,double timestep,unsigned step
  ,unsigned np,const tdouble3* pos,const unsigned* ridpmot)
{
  if(step==0 && pos==NULL && ridpmot==NULL){
    MotData2->SavePart(cpart,timestep,step,PsCount,PosRef);
  }
  else{
    //-Obtains current reference positions (or use initial data for step=0).
    LoadPosRefCpu(timestep,step,np,pos,ridpmot);
    //-Saves data in BI4 file.
    MotData2->AddDataPart(cpart,LastTimestep,LastStep,PsCount,PosRef);
  }
}

//==============================================================================
/// Saves stored EXTRA reference position data in BI4 file.
//==============================================================================
void JDsPartMotionSave::SaveDataExtra(){
  if(MotData2)MotData2->SaveStoredData();
}

#ifdef _WITHGPU
//==============================================================================
/// Free GPU memory.
//==============================================================================
void JDsPartMotionSave::ResetDataGpu(){
  if(IdpRefg)cudaFree(IdpRefg); IdpRefg=NULL;
  if(PosRefg)cudaFree(PosRefg); PosRefg=NULL;
}
//==============================================================================
/// Configure data for GPU execution.
//==============================================================================
void JDsPartMotionSave::ConfigDataGpu(){
  if(IdpRef==NULL)Run_Exceptioon("Error: Object not ready for GPU configuration.");
  if(IdpRefg!=NULL || PosRefg!=NULL)Run_Exceptioon("Error: GPU data is already configured.");
  //-Allocates GPU memory.
  fcuda::Malloc(&IdpRefg,PsCount);
  fcuda::Malloc(&PosRefg,PsCount);
  fcuda::Check_CudaErroorFun("Memory allocation.");
  //-Copy data to GPU memory.
  cudaMemcpy(IdpRefg,IdpRef,sizeof(unsigned)*PsCount,cudaMemcpyHostToDevice);
  cudaMemset(PosRefg,0,sizeof(tdouble3)*PsCount);
  fcuda::Check_CudaErroorFun("Copy to GPU memory.");
}

//==============================================================================
/// Load current reference position data from particle data on GPU.
//==============================================================================
void JDsPartMotionSave::LoadPosRefGpu(double timestep,unsigned step,unsigned np
  ,const double2* posxy,const double* posz,const unsigned* ridpmot)
{
  if(LastStep!=step || LastTimestep!=timestep){
    if(IdpRefg==NULL)ConfigDataGpu();
    cusph::LoadPosRef(PsCount,CaseNfixed,np,posxy,posz,ridpmot,IdpRefg,PosRefg);
    cudaMemcpy(PosRef,PosRefg,sizeof(tdouble3)*PsCount,cudaMemcpyDeviceToHost);
    LastStep=step;
    LastTimestep=timestep;
  }
}

//==============================================================================
/// Saves MAIN current reference position data from GPU in BI4 file.
//==============================================================================
void JDsPartMotionSave::SaveDataMainGpu(int cpart,double timestep,unsigned step
  ,unsigned np,const double2* posxy,const double* posz,const unsigned* ridpmot)
{
  //-Obtains current reference positions.
  LoadPosRefGpu(timestep,step,np,posxy,posz,ridpmot);
  //-Saves data in BI4 file.
  MotData1->SavePart(cpart,LastTimestep,LastStep,PsCount,PosRef);
}

//==============================================================================
/// Saves EXTRA current reference position data from GPU in BI4 file.
//==============================================================================
void JDsPartMotionSave::AddDataExtraGpu(int cpart,double timestep,unsigned step
  ,unsigned np,const double2* posxy,const double* posz,const unsigned* ridpmot)
{
  //-Obtains current reference positions.
  LoadPosRefGpu(timestep,step,np,posxy,posz,ridpmot);
  //-Saves data in BI4 file.
  MotData2->AddDataPart(cpart,LastTimestep,LastStep,PsCount,PosRef);
}
#endif

