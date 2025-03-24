//HEAD_DSPH
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

/// \file JDsVresData.cpp \brief Implements the class \ref JDsVResData.

#include "JDsVresData.h"
#include "JBinaryData.h"
#include "JLog2.h"
#include "Functions.h"
#include "JRangeFilter.h"

#include <algorithm>
#include <cstring>

using namespace std;

//##############################################################################
//# JDsVResDataSave
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsVResDataSave::JDsVResDataSave(std::string appname,std::string dir,JLog2* log)
    :AppName(appname),Dir(fun::GetDirWithSlash(dir)),Log(log)
{
  ClassName="JDsVResDataSave";
  FilterParts=NULL;
  Data=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsVResDataSave::~JDsVResDataSave(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsVResDataSave::Reset(){
  SvParts=0;
  delete FilterParts; FilterParts=NULL;
  Cpart=0;
  delete Data; Data=NULL;
}

//==============================================================================
/// Defines configuration.
//==============================================================================
void JDsVResDataSave::Config(std::string svparts){
  Reset();
  if(fun::StrOnlyChars(svparts,"0123456789"))SvParts=fun::StrToInt(svparts);
  else{
    FilterParts=new JRangeFilter(svparts);
    SvParts=1;
  }
  //Log->PrintfDbg("\n[%s]",svparts.c_str());
  //for(int c=0;c<100;c++)if(CheckSave(c))Log->PrintfDbg("%02d]: %d",c,(CheckSave(c)? 1: 0));
}

//==============================================================================
/// Check PART number to save.
//==============================================================================
bool JDsVResDataSave::CheckSave(int cpart)const{
  return(cpart>0 && SvParts && cpart%SvParts==0 
    && (FilterParts==NULL || FilterParts->CheckValue(unsigned(cpart))));
}

//==============================================================================
/// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JDsVResDataSave::GetFileNamePart(int cpart){
  return(fun::PrintStr("PartVResData_%04u.bi4",cpart));
}

//==============================================================================
/// Initialise PART data to save.
//==============================================================================
void JDsVResDataSave::InitPartData(int cpart,double timestep,int nstep){
  Cpart=cpart;
  delete Data;
  Data=new JBinaryData("JPartVResBi4");
  Data->SetvText("AppName",AppName);
  Data->SetvUint("FormatVer",FormatVerDef);
  Data->SetvInt("Cpart",cpart);
  Data->SetvUint("Step",nstep);
  Data->SetvDouble("TimeStep",timestep);
}

//==============================================================================
/// Add normals data.
//==============================================================================
void JDsVResDataSave::AddArray(unsigned np,unsigned msize,const tdouble3* pos,const tfloat3* normals
  ,const tfloat3* velmot,const float* mass,const double* matarray)
{
  const unsigned nsize=np;
  //Log->Printf("AddNormals----> np:%u  npb:%u  nsize:%u",np,npb,nsize);
  tdouble3* vpos    =Data->CreateArrayDouble3("Pos",nsize,true);
  tfloat3*  vnor    =Data->CreateArrayFloat3("Normals",nsize,true);
  tfloat3*  vvm     =Data->CreateArrayFloat3("VelMot",nsize,true);
  float*    vmass   =Data->CreateArrayFloat("Mass",nsize,true);
  double*   mat     =Data->CreateArrayDouble("Mat",msize,true);

  for(unsigned p=0;p<np;p++){
    vpos[p] =pos[p];
    vnor[p] =normals[p];
    vvm[p]  =velmot[p];
    vmass[p]=mass[p];
  }
  for(unsigned p=0;p<msize;p++){
    mat[p] =matarray[p];
  }
}

//==============================================================================
/// Saves file with extra PART data.
//==============================================================================
void JDsVResDataSave::SavePartData(){
  const string file=Dir+GetFileNamePart(Cpart);
  Data->SaveFile(file,false,true);
  //-Clear data.
  Cpart=0;
  delete Data; Data=NULL;
}


//##############################################################################
//# JDsVResDataLoad
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsVResDataLoad::JDsVResDataLoad(JLog2* log):Log(log)
{
  ClassName="JDsVResDataLoad";
  Data=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsVResDataLoad::~JDsVResDataLoad(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsVResDataLoad::Reset(){
  Dir="";
  Cpart=0;
  FileData="";
  delete Data; Data=NULL;
}

//==============================================================================
/// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JDsVResDataLoad::GetFileNamePart(std::string dir,int cpart){
  return(fun::GetDirWithSlash(dir)+fun::PrintStr("PartVResData_%04u.bi4",cpart));
}

//==============================================================================
/// Check if PART file exists.
//==============================================================================
bool JDsVResDataLoad::ExistsPartData(std::string dir,int cpart){
  return(fun::FileExists(GetFileNamePart(dir,cpart)));
}

//==============================================================================
/// Loads file with extra PART data.
//==============================================================================
void JDsVResDataLoad::LoadPartData(std::string dir,int cpart){
  Reset();
  Dir=dir;
  Cpart=cpart;
  FileData=GetFileNamePart(Dir,Cpart);
  Data=new JBinaryData("JPartVResBi4");
  Data->LoadFile(FileData);
}

//==============================================================================
/// Load normals data from extra data.
//==============================================================================
void JDsVResDataLoad::LoadArray(unsigned np,unsigned msize,tfloat3* velmot
  ,float* mass,double* matarray)
{
  JBinaryDataArray* ar=Data->GetArray("Mat");
  const double* vmat=(const double* )ar->GetDataPointer();
  for(unsigned p=0;p<msize;p++)matarray[p]=vmat[p];

  ar=Data->GetArray("Mass");
  const float* vmass=(const float* )ar->GetDataPointer();
  for(unsigned p=0;p<msize;p++)mass[p]=vmass[p];

  ar=Data->GetArray("VelMot");
  const tfloat3* vvm=(const tfloat3* )ar->GetDataPointer();
  for(unsigned p=0;p<msize;p++)velmot[p]=vvm[p];
  
}

//==============================================================================
/// Return motion matrix. Only one matrix is valid.
//==============================================================================
tmatrix4d JDsVResDataLoad::GetMat(){
  tmatrix4d mat=TMatrix4d();
  JBinaryDataArray* ar=Data->GetArray("Mat");
  if(!ar || ar->GetCount()!=16)
    Run_ExceptioonFile("Matrix data is missing or invalid.",FileData);
  memcpy(&mat,ar->GetDataPointer(),sizeof(tmatrix4d));
  return(mat);
}

