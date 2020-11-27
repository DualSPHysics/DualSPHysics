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

/// \file JPartNormalData.cpp \brief Implements the class \ref JPartNormalData.

#include "JPartNormalData.h"
#include "JBinaryData.h"
#include "Functions.h"
#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstring>

using namespace std;

//##############################################################################
//# JPartNormalData
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartNormalData::JPartNormalData(){
  ClassName="JPartNormalData";
  PartNormals=NULL;
  NormalBegin=NULL;
  Normals=NULL; NormalsDist=NULL;
  OutVecs=NULL; OutVecsDist=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartNormalData::~JPartNormalData(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartNormalData::Reset(){
  FmtVersion=FmtVersionDef;
  DirData="";
  ConfigBasic("","",false,0,0,0,0,false);
  PartNormalsName="";
  AllocNormals(0,0,false);
}

//==============================================================================
/// Devuelve true si existe el array indicado.
/// Returns true if the indicated array exists.
//==============================================================================
bool JPartNormalData::ArrayExists(JBinaryData *bd,std::string name)const{
  return(bd->GetArray(name)!=NULL);
}

//==============================================================================
/// Devuelve el puntero al array indicado.
/// Returns the pointer to the indicated array.
//==============================================================================
JBinaryDataArray* JPartNormalData::GetArray(JBinaryData *bd,std::string name)const{
  JBinaryDataArray* ar=bd->GetArray(name);
  if(!ar)Run_Exceptioon(fun::PrintStr("Array \'%s\' is not available.",name.c_str()));
  return(ar);
}

//==============================================================================
/// Devuelve el puntero al array indicado y comprueba el tipo.
/// Returns the pointer to the indicated array and checks the type.
//==============================================================================
JBinaryDataArray* JPartNormalData::GetArray(JBinaryData *bd,std::string name,JBinaryDataDef::TpData type)const{
  JBinaryDataArray* ar=bd->GetArray(name);
  if(ar->GetType()!=type)Run_Exceptioon(fun::PrintStr("Type of array \'%s\' is not %s.",name.c_str(),JBinaryDataDef::TypeToStr(type).c_str()));
  return(ar);
}

//==============================================================================
/// Configuracion de variables basicas.
/// Configuration of basic variables.
//==============================================================================
void JPartNormalData::ConfigBasic(std::string appname,std::string casename
  ,bool data2d,double data2dposy,double dp,double h,double dist,bool ftsupport)
{
  AppName=appname;
  Date=fun::GetDateTime();
  CaseName=casename;
  Data2d=data2d;
  Data2dPosY=data2dposy;
  Dp=dp;
  H=h;
  Dist=dist;
  FtSupport=ftsupport;
}

//==============================================================================
/// Allocates or frees memory for particles.
//==============================================================================
void JPartNormalData::AllocNormals(unsigned nbound,unsigned countnor,bool usepartnormals){
  Nbound=CountNormals=0;
  delete[] PartNormals; PartNormals=NULL;
  delete[] NormalBegin; NormalBegin=NULL;
  delete[] Normals;     Normals=NULL;
  delete[] NormalsDist; NormalsDist=NULL;
  delete[] OutVecs;     OutVecs=NULL;
  delete[] OutVecsDist; OutVecsDist=NULL;
  Nbound=nbound;
  CountNormals=countnor;
  try{
    if(Nbound && usepartnormals)PartNormals=new tdouble3[Nbound];
    if(Nbound && CountNormals){
      NormalBegin=new unsigned[Nbound+1];
      Normals    =new tdouble3[CountNormals];
      NormalsDist=new double  [CountNormals];
      OutVecs    =new tdouble3[CountNormals];
      OutVecsDist=new double  [CountNormals];
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Add normal data and final normals of particles.
//==============================================================================
void JPartNormalData::AddNormalData(
  std::string partnorname,unsigned nbound,const tdouble3 *partnor
  ,const unsigned *norbegin,unsigned countnor
  ,const tdouble3 *nordata,const double *nordist
  ,const tdouble3 *outdata,const double *outdist)
{
  PartNormalsName=partnorname;
  AllocNormals(nbound,countnor,partnor!=NULL);
  if(partnor)memcpy(PartNormals,partnor,sizeof(tdouble3)*Nbound);
  if(CountNormals){
    memcpy(NormalBegin,norbegin,sizeof(unsigned)*(Nbound+1));
    memcpy(Normals    ,nordata ,sizeof(tdouble3)*CountNormals);
    memcpy(NormalsDist,nordist ,sizeof(double)  *CountNormals);
    memcpy(OutVecs    ,outdata ,sizeof(tdouble3)*CountNormals);
    memcpy(OutVecsDist,outdist ,sizeof(double)  *CountNormals);
  }
}

//==============================================================================
/// Add final normals of particles.
//==============================================================================
void JPartNormalData::AddNormalData(
  std::string partnorname,unsigned nbound,const tdouble3 *partnor)
{
  AddNormalData(partnorname,nbound,partnor,NULL,0,NULL,NULL,NULL,NULL);
}

//==============================================================================
/// Returns file name and path.
//==============================================================================
std::string JPartNormalData::GetFileName(std::string casename,std::string dir){
  return(fun::GetDirWithSlash(dir)+casename+"_NormalData.nbi4");
}

//==============================================================================
/// Saves binary file Part_Head.ibi4.
//==============================================================================
void JPartNormalData::SaveFile(std::string dir){
  DirData=fun::GetDirWithSlash(dir);
  JBinaryData bdat(ClassName);
  bdat.SetvUint("FmtVersion",FmtVersion);
  //-Saves general variables.
  bdat.SetvText("AppName",AppName);
  bdat.SetvText("Date",Date);
  //bdat.SetvText("AppName","XxX"); //AppName);
  //bdat.SetvText("Date","XxX"); //Date);
  bdat.SetvText("CaseName",CaseName);
  bdat.SetvBool("Data2d",Data2d);
  bdat.SetvDouble("Data2dPosY",Data2dPosY);
  bdat.SetvDouble("Dp",Dp);
  bdat.SetvDouble("H",H);
  //-Saves normal data.
  bdat.SetvDouble("Dist",Dist);
  bdat.SetvBool("FtSupport",FtSupport);
  bdat.SetvText("PartNormalsName",PartNormalsName);
  bdat.SetvUint("Nbound",Nbound);
  bdat.SetvUint("CountNormals",CountNormals);
  if(Nbound && CountNormals){
    bdat.CreateArray("NormalBegin",JBinaryDataDef::DatUint   ,Nbound+1    ,NormalBegin,true);
    bdat.CreateArray("Normals"    ,JBinaryDataDef::DatDouble3,CountNormals,Normals    ,true);
    bdat.CreateArray("NormalsDist",JBinaryDataDef::DatDouble ,CountNormals,NormalsDist,true);
    bdat.CreateArray("OutVecs"    ,JBinaryDataDef::DatDouble3,CountNormals,OutVecs    ,true);
    bdat.CreateArray("OutVecsDist",JBinaryDataDef::DatDouble ,CountNormals,OutVecsDist,true);
  }
  if(PartNormals)bdat.CreateArray("PartNormals",JBinaryDataDef::DatDouble3,Nbound,PartNormals,true);
  bdat.SaveFile(DirData+GetFileName(CaseName),false,true);
}

//==============================================================================
/// Returns the file name with normal data.
//==============================================================================
std::string JPartNormalData::GetNormalDataFile(std::string casename){
  const string dirdata=fun::GetDirWithSlash(fun::GetDirParent(casename));
  const string casenam=fun::GetWithoutExtension(fun::GetFile(casename));
  return(GetFileName(dirdata+casenam));
}

//==============================================================================
/// Loads binary file Part_Head.ibi4.
//==============================================================================
void JPartNormalData::LoadFile(std::string casename){
  Reset();
  DirData=fun::GetDirWithSlash(fun::GetDirParent(casename));
  CaseName=fun::GetWithoutExtension(fun::GetFile(casename));
  //printf("----> dir:[%s] case:[%s]\n",DirData.c_str(),CaseName.c_str());
  JBinaryData bdat;
  string file=GetFileName(DirData+CaseName);
  //printf("----> file:[%s]\n",file.c_str());
  bdat.LoadFile(file,ClassName);
  FmtVersion=bdat.GetvUint("FmtVersion");
  //-Loads general variables.
  AppName   =bdat.GetvText("AppName");
  Date      =bdat.GetvText("Date");
  CaseName  =bdat.GetvText("CaseName");
  Data2d    =bdat.GetvBool("Data2d");
  Data2dPosY=bdat.GetvDouble("Data2dPosY");
  Dp        =bdat.GetvDouble("Dp");
  H         =bdat.GetvDouble("H");
  //-Loads normal data.
  Dist           =bdat.GetvDouble("Dist");
  FtSupport      =bdat.GetvBool("FtSupport",true,false);
  PartNormalsName=bdat.GetvText("PartNormalsName");
  Nbound         =bdat.GetvUint("Nbound");
  CountNormals   =bdat.GetvUint("CountNormals");
  bool partnormals=(bdat.GetArray("PartNormals")!=NULL);
  //printf("----> PartNormalsName:[%s]\n",PartNormalsName.c_str());
  //printf("----> Nbound:[%d]\n",Nbound);
  //printf("----> CountNormals:[%d]\n",CountNormals);
  AllocNormals(Nbound,CountNormals,partnormals);
  if(Nbound && CountNormals){
    GetArray(&bdat,"NormalBegin",JBinaryDataDef::DatUint   )->GetDataCopy(Nbound+1    ,NormalBegin);
    GetArray(&bdat,"Normals"    ,JBinaryDataDef::DatDouble3)->GetDataCopy(CountNormals,Normals);
    GetArray(&bdat,"NormalsDist",JBinaryDataDef::DatDouble )->GetDataCopy(CountNormals,NormalsDist);
    GetArray(&bdat,"OutVecs"    ,JBinaryDataDef::DatDouble3)->GetDataCopy(CountNormals,OutVecs);
    GetArray(&bdat,"OutVecsDist",JBinaryDataDef::DatDouble )->GetDataCopy(CountNormals,OutVecsDist);
  }
  if(PartNormals)GetArray(&bdat,"PartNormals",JBinaryDataDef::DatDouble3)->GetDataCopy(Nbound,PartNormals);
}


