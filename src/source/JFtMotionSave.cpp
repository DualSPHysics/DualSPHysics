//HEAD_DSCODES
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

/// \file JFtMotionSave.cpp \brief Implements the class \ref JFtMotionSave.

#include "JFtMotionSave.h"
#include "Functions.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JPartFloatBi4.h"
#include "JComputeMotionRef.h"

#ifdef _WITHGPU
#include "JSphGpu_ker.h"
#endif

#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstring>
#include <cmath>

using namespace std;

//##############################################################################
//# JFtMotionSave
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JFtMotionSave::JFtMotionSave(double tout):Log(AppInfo.LogPtr()),TimeOut(tout){
  ClassName="JFtMotionSave";
  FtMks=NULL;
  IdpRef=NULL;  PosRef=NULL;
  IdpRefg=NULL; PosRefg=NULL;
  FtData=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JFtMotionSave::~JFtMotionSave(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JFtMotionSave::Reset(){
  FtCount=0;
  MkBoundFirst=0;
  CaseFtBegin=CaseNfloat=0;
  delete[] FtMks;  FtMks=NULL;
  delete FtData;  FtData=NULL;
  Num=0;
  NextTimeOutput=0;
  //-Auxiliary memory.
  delete[] IdpRef; IdpRef=NULL;
  delete[] PosRef; PosRef=NULL;
#ifdef _WITHGPU
  //-Free GPU memory.
  if(IdpRefg)cudaFree(IdpRefg); IdpRefg=NULL;
  if(PosRefg)cudaFree(PosRefg); PosRefg=NULL;
#endif
}

//==============================================================================
/// Initial configuration of floating data.
//==============================================================================
void JFtMotionSave::Config(std::string appname,std::string dirout,word mkboundfirst
  ,unsigned ftcount,const StFloatingData *ftobjs
  ,unsigned np,const tdouble3 *pos,const unsigned *idp)
{
  Reset();
  FtCount=ftcount;
  MkBoundFirst=mkboundfirst;
  if(!FtCount)Run_Exceptioon("No floating bodies available.");
  //-Loads data of FtMks[].
  FtMks=new StMkMotionData[FtCount];
  CaseNfloat=0;
  for(unsigned cf=0;cf<FtCount;cf++){
    const StFloatingData &ft=ftobjs[cf];
    StMkMotionData &v=FtMks[cf];
    v.mkbound=ft.mkbound;
    v.begin=ft.begin;
    v.np=ft.count;
    v.nid=0;
    CaseNfloat+=ft.count;
  }
  CaseFtBegin=FtMks[0].begin;
  //-Create object for output data.
  FtData=new JPartFloatBi4Save();
  FtData->Config(appname,dirout,mkboundfirst,FtCount,true);
  for(unsigned cf=0;cf<FtCount;cf++){
    const StFloatingData &ft=ftobjs[cf];
    FtData->AddHeadData(cf,ft.mkbound,ft.begin,ft.count,ft.mass,ft.massp,ft.radius);
  }
  FtData->SaveInitial("PartFloatMotion.fbi4");
  //-Allocates auxiliary memory.
  IdpRef=new unsigned[FtCount*3];
  PosRef=new tdouble3[FtCount*3];
  //-Configures reference position.
  ConfigPosRef(np,pos,idp);
  //-Saves first data.
  SaveFtData(0,0,ftobjs);
}

//==============================================================================
/// Initial configuration of floating data.
//==============================================================================
void JFtMotionSave::ConfigPosRef(unsigned np,const tdouble3 *pos,const unsigned *idp){
  //-Find reference points to compute motion.
  JComputeMotionRef computemotref;
  computemotref.AddMkBlocks(MkBoundFirst,FtCount,FtMks);
  computemotref.ComputeRefPoints(CaseFtBegin,0,CaseNfloat,np,idp,pos,NULL,NULL);
  computemotref.GetMotionRef(FtCount,FtMks);
  computemotref.Reset();
  //-Load IdpRef[].
  for(unsigned cf=0;cf<FtCount;cf++){
    const StMkMotionData &v=FtMks[cf];
    for(unsigned ci=0;ci<3;ci++){
      IdpRef[cf*3+ci]=(ci<v.nid? v.id[ci]-CaseFtBegin: UINT_MAX);
      PosRef[cf*3+ci]=(ci<v.nid? v.ps[ci]: TDouble3(0));
    }
  }
}

//==============================================================================
/// Returns next time to save PART file.
//==============================================================================
double JFtMotionSave::GetNextTime(double t)const{
  double tnext=t;
  if(TimeOut){
    unsigned ct=unsigned(t/TimeOut);
    tnext=TimeOut*ct;
    for(;tnext<=t;ct++)tnext=TimeOut*ct;
  }
  return(tnext);
}

//==============================================================================
/// Saves floating data from simulation.
//==============================================================================
void JFtMotionSave::SaveFtData(double timestep,unsigned nstep,const StFloatingData *ftobjs){
  //Log->Printf("===> SaveFtDataCpu> t:%f > %f  ns:%u",timestep,NextTimeOutput,nstep);
  //-Stores data from ftobjs.
  for(unsigned cf=0;cf<FtCount;cf++){
    const StFloatingData &v=ftobjs[cf];
    FtData->AddPartData(cf,v.center,v.fvel,v.fomega,v.facelin,v.faceang);
  }
  //-Stores position data.
  FtData->AddPartDataPosRef(FtCount,PosRef);
  //-Saves data in file.
  FtData->SavePartFloat(Num,nstep,timestep,0);
  //-Prepare next data.
  Num++;
  NextTimeOutput=GetNextTime(timestep);
}

//==============================================================================
/// Saves floating data from CPU simulation.
//==============================================================================
void JFtMotionSave::SaveFtDataCpu(double timestep,unsigned nstep,const StFloatingData *ftobjs
  ,unsigned np,const tdouble3 *pos,const unsigned *ftridp)
{
  //-Gets position data.
  const unsigned n=FtCount*3;
  for(unsigned cp=0;cp<n;cp++){
    const unsigned cid=IdpRef[cp];
    if(cid!=UINT_MAX){
      const unsigned p=ftridp[cid];
      PosRef[cp]=(p<np? pos[p]: TDouble3(DBL_MAX));
    }
  }
  //-Saves data.
  SaveFtData(timestep,nstep,ftobjs);
}

#ifdef _WITHGPU
//==============================================================================
/// Saves floating data from GPU simulation.
//==============================================================================
void JFtMotionSave::SaveFtDataGpu(double timestep,unsigned nstep,const StFloatingData *ftobjs
  ,unsigned np,const double2 *posxyg,const double *poszg,const unsigned *ftridpg)
{
  //-Allocates GPU memory when it is necessary (first time).
  if(IdpRefg==NULL){
    cudaMalloc((void**)&IdpRefg,sizeof(unsigned)*FtCount*3);
    cudaMalloc((void**)&PosRefg,sizeof(double)  *FtCount*9);
    cudaMemcpy(IdpRefg,IdpRef,sizeof(unsigned)*FtCount*3,cudaMemcpyHostToDevice);
  }
  //-Gets position data.
  cusph::FtGetPosRef(FtCount*3,IdpRefg,ftridpg,posxyg,poszg,PosRefg);
  cudaMemcpy(PosRef,PosRefg,sizeof(double)*FtCount*9,cudaMemcpyDeviceToHost);
  //-Saves data.
  SaveFtData(timestep,nstep,ftobjs);
}
#endif

