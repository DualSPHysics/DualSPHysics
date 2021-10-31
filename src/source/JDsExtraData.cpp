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

/// \file JDsExtraData.cpp \brief Implements the class \ref JDsExtraData.

#include "JDsExtraData.h"
#include "JBinaryData.h"
#include "JLog2.h"
#include "Functions.h"
#include "JRangeFilter.h"
#include <algorithm>

using namespace std;

//##############################################################################
//# JDsExtraDataSave
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsExtraDataSave::JDsExtraDataSave(std::string appname,std::string dir
  ,unsigned casenbound,unsigned casenfloat,JLog2* log):AppName(appname)
  ,Dir(fun::GetDirWithSlash(dir)),CaseNbound(casenbound),CaseNfloat(casenfloat)
  ,Log(log)
{
  ClassName="JDsExtraDataSave";
  FilterParts=NULL;
  Data=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsExtraDataSave::~JDsExtraDataSave(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsExtraDataSave::Reset(){
  SvParts=0;
  delete FilterParts; FilterParts=NULL;
  Cpart=0;
  delete Data; Data=NULL;
}

//==============================================================================
/// Defines configuration.
//==============================================================================
void JDsExtraDataSave::Config(std::string svparts){
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
bool JDsExtraDataSave::CheckSave(int cpart)const{
  return(cpart>0 && SvParts && cpart%SvParts==0 && (FilterParts==NULL || FilterParts->CheckValue(unsigned(cpart))));
}

//==============================================================================
/// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JDsExtraDataSave::GetFileNamePart(int cpart){
  return(fun::PrintStr("PartExtra_%04u.bi4",cpart));
}

//==============================================================================
/// Initialise PART data to save.
//==============================================================================
void JDsExtraDataSave::InitPartData(int cpart,double timestep,int nstep){
  Cpart=cpart;
  delete Data;
  Data=new JBinaryData("JPartExtraBi4");
  Data->SetvText("AppName",AppName);
  Data->SetvUint("FormatVer",FormatVerDef);
  Data->SetvUint("CaseNbound",CaseNbound);
  Data->SetvUint("CaseNfloat",CaseNfloat);
  Data->SetvInt("Cpart",cpart);
  Data->SetvUint("Step",nstep);
  Data->SetvDouble("TimeStep",timestep);
}

//==============================================================================
/// Add normals data.
//==============================================================================
void JDsExtraDataSave::AddNormals(bool usenormalsft,unsigned np,unsigned npb
  ,const unsigned *idp,const typecode *code,const tfloat3 *boundnormal)
{
  if(!CaseNfloat)usenormalsft=false;
  Data->SetvBool("UseNormalsFt",usenormalsft);
  const unsigned nsize=(usenormalsft? CaseNbound: CaseNbound-CaseNfloat);
  //Log->Printf("AddNormals----> np:%u  npb:%u  nsize:%u",np,npb,nsize);
  tfloat3 *vnor=Data->CreateArrayFloat3("Normals",nsize,true);
  //-Fixed and moving boundary particles.
  if(code!=NULL){
    for(unsigned p=0;p<npb;p++)if(idp[p]<nsize && CODE_IsNormal(code[p]))vnor[idp[p]]=boundnormal[p];
  }
  else for(unsigned p=0;p<npb;p++)if(idp[p]<nsize)vnor[idp[p]]=boundnormal[p]; 
  //-Floating boundary particles.
  if(usenormalsft && CaseNfloat){
    if(code!=NULL){
      for(unsigned p=npb;p<np;p++)if(idp[p]<nsize && CODE_IsNormal(code[p]))vnor[idp[p]]=boundnormal[p];
    }
    else for(unsigned p=npb;p<np;p++)if(idp[p]<nsize)vnor[idp[p]]=boundnormal[p];
  }
}

//==============================================================================
/// Saves file with extra PART data.
//==============================================================================
void JDsExtraDataSave::SavePartData(){
  const string file=Dir+GetFileNamePart(Cpart);
  Data->SaveFile(file,false,true);
  //-Clear data.
  Cpart=0;
  delete Data; Data=NULL;
}


//##############################################################################
//# JDsExtraDataLoad
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsExtraDataLoad::JDsExtraDataLoad(unsigned casenbound,unsigned casenfloat
  ,JLog2* log):CaseNbound(casenbound),CaseNfloat(casenfloat),Log(log)
{
  ClassName="JDsExtraDataLoad";
  Data=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsExtraDataLoad::~JDsExtraDataLoad(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsExtraDataLoad::Reset(){
  Dir="";
  Cpart=0;
  FileData="";
  delete Data; Data=NULL;
}

//==============================================================================
/// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JDsExtraDataLoad::GetFileNamePart(std::string dir,int cpart){
  return(fun::GetDirWithSlash(dir)+fun::PrintStr("PartExtra_%04u.bi4",cpart));
}

//==============================================================================
/// Check if PART file exists.
//==============================================================================
bool JDsExtraDataLoad::ExistsPartData(std::string dir,int cpart){
  return(fun::FileExists(GetFileNamePart(dir,cpart)));
}

//==============================================================================
/// Loads file with extra PART data.
//==============================================================================
void JDsExtraDataLoad::LoadPartData(std::string dir,int cpart){
  Reset();
  Dir=dir;
  Cpart=cpart;
  FileData=GetFileNamePart(Dir,Cpart);
  Data=new JBinaryData("JPartExtraBi4");
  Data->LoadFile(FileData);
}

//==============================================================================
/// Load normals data from extra data.
//==============================================================================
bool JDsExtraDataLoad::LoadNormals(unsigned np,unsigned npb
  ,const unsigned *idp,tfloat3 *boundnormal)
{
  if(Data->GetvUint("CaseNbound")!=CaseNbound)Run_ExceptioonFile("CaseNbound value does not match.",FileData);
  if(Data->GetvUint("CaseNfloat")!=CaseNfloat)Run_ExceptioonFile("CaseNfloat value does not match.",FileData);
  const bool usenormalsft=Data->GetvBool("UseNormalsFt");
  JBinaryDataArray *ar=Data->GetArray("Normals");
  if(!ar || ar->GetType()!=JBinaryDataDef::DatFloat3)Run_ExceptioonFile("The array \'Normals\' is missing or type invalid.",FileData);
  const unsigned nsize=ar->GetCount();
  const tfloat3 *vnor=(const tfloat3 *)ar->GetDataPointer();
  //-Fixed and moving boundary particles.
  for(unsigned p=0;p<npb;p++)if(idp[p]<nsize)boundnormal[p]=vnor[idp[p]]; 
  //-Floating boundary particles.
  if(usenormalsft){
    for(unsigned p=npb;p<np;p++)if(idp[p]<nsize)boundnormal[p]=vnor[idp[p]];
  }
  return(usenormalsft);
}

