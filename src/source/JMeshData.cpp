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

/// \file JMeshData.cpp \brief Implements the class \ref JMeshData.

#include "JMeshData.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JDataArrays.h"
#include "JSaveCsv2.h"
#include "JException.h"
#include <cfloat>
#include <climits>
#include <algorithm>
#include <cstring>

using namespace std;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshData
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMeshData::JMeshData(){
  ClassName="JMeshData";
  PtPos=NULL;
  Arrays=NULL;
  Reset(true);
}

//==============================================================================
/// Destructor.
//==============================================================================
JMeshData::~JMeshData(){
  DestructorActive=true;
  Reset(false);
}

//==============================================================================
/// Throws exception related to a file from a static method.
//==============================================================================
void JMeshData::RunExceptioonStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,const std::string& msg,const std::string& file)
{
  throw JException(srcfile,srcline,"JMeshData",method,msg,file);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JMeshData::Reset(bool createarrays){
  UseMesh=true;
  PtRef=Vdp1=Vdp2=Vdp3=TDouble3(0);
  delete[] PtPos; PtPos=NULL;
  Npt1=Npt2=Npt3=Npt=Npt12=0;
  DirDat=TFloat3(0);
  TimeStep=0;
  delete Arrays; Arrays=NULL;
  if(createarrays)Arrays=new JDataArrays();
}

//==============================================================================
/// Copies data from other object.
//==============================================================================
void JMeshData::CopyDataFrom(const JMeshData* obj){
  TimeStep=obj->TimeStep;
  Arrays->CopyDataFrom(*(obj->Arrays),true);
}

//==============================================================================
/// Adds special arrays in varlist ("name1:name2:name3).
//==============================================================================
void JMeshData::AddVarList(std::string varlist){
  const bool cleardata=true;
  std::vector<string> names;
  const unsigned n=fun::VectorSplitStr(",",varlist,names);
  for(unsigned c=0;c<n;c++){
    const string name=names[c];
    const string namelow=fun::StrLower(names[c]);
    if(Arrays->ExistsName(name) || Arrays->ExistsName(namelow))Run_Exceptioon(fun::PrintStr("The array \'%s\' already exists.",name.c_str()));
         if(namelow=="vel"   )Arrays->CreateArrayFloat3(name,Npt,cleardata);
    else if(namelow=="veldir")Arrays->CreateArrayFloat (name,Npt,cleardata);
    else if(namelow=="rhop"  )Arrays->CreateArrayFloat (name,Npt,cleardata);
    else if(namelow=="zsurf" ){
      Arrays->GetArray(Arrays->CreateArrayFloat(name,Npt12,cleardata)).tag=TAGDATA12;
    }  
    else Run_Exceptioon(fun::PrintStr("The name of array \'%s\' is not a predefined array.",name.c_str()));
  }
}

//==============================================================================
/// Reverse selected variables "name1,name2,name3" (v=-v).
//==============================================================================
void JMeshData::ReverseData(std::string varlist){
  std::vector<string> names;
  fun::VectorSplitStr(",",varlist,names);
  fun::VectorLower(names);
  const bool all=(fun::VectorFind("all",names)!=UINT_MAX);
  const unsigned na=Arrays->Count();
  for(unsigned ca=0;ca<na;ca++){
    const JDataArrays::StDataArray& ar=Arrays->GetArrayCte(ca);
    if(all || fun::VectorFind(fun::StrLower(ar.keyname),names)!=UINT_MAX){
      Arrays->ReverseArrayData(ca);
    }
  }
}

//==============================================================================
/// Set data of selected variables "name1,name1.x,name1.y,name2.z,all" (v*=v2).
//==============================================================================
void JMeshData::SetMulData(std::string varlist,double v2){
  std::vector<string> names;
  fun::VectorSplitStr(",",varlist,names);
  fun::VectorLower(names);
  const unsigned na=Arrays->Count();
  const bool all=(fun::VectorFind("all",names)!=UINT_MAX);
  if(all){
    for(unsigned ca=0;ca<na;ca++){
      const JDataArrays::StDataArray& ar=Arrays->GetArrayCte(ca);
      Arrays->SetMulArrayData(ca,' ',v2);
    }
  }
  const unsigned nv=unsigned(names.size());
  for(unsigned cv=0;cv<nv;cv++){
    string namecc=names[cv];
    const string name=fun::StrSplit(".",namecc);
    const char namec=(namecc.size()>=1? namecc[0]: ' ') ;
    for(unsigned ca=0;ca<na;ca++){
      const JDataArrays::StDataArray& ar=Arrays->GetArrayCte(ca);
      if(fun::StrLower(ar.keyname)==name){
        Arrays->SetMulArrayData(ca,namec,v2);
      }
    }
  }
}

//==============================================================================
/// Set data of selected variables "name1,name1.x,name1.y,name2.z,all" (v+=v2).
//==============================================================================
void JMeshData::SetAddData(std::string varlist,double v2){
  std::vector<string> names;
  fun::VectorSplitStr(",",varlist,names);
  fun::VectorLower(names);
  const unsigned na=Arrays->Count();
  const bool all=(fun::VectorFind("all",names)!=UINT_MAX);
  if(all){
    for(unsigned ca=0;ca<na;ca++){
      const JDataArrays::StDataArray& ar=Arrays->GetArrayCte(ca);
      Arrays->SetAddArrayData(ca,' ',v2);
    }
  }
  const unsigned nv=unsigned(names.size());
  for(unsigned cv=0;cv<nv;cv++){
    string namecc=names[cv];
    const string name=fun::StrSplit(".",namecc);
    const char namec=(namecc.size()>=1? namecc[0]: ' ') ;
    for(unsigned ca=0;ca<na;ca++){
      const JDataArrays::StDataArray& ar=Arrays->GetArrayCte(ca);
      if(fun::StrLower(ar.keyname)==name){
        Arrays->SetAddArrayData(ca,namec,v2);
      }
    }
  }
}

//==============================================================================
/// Creates array float, checks size, copy data and returns memory pointer.
//==============================================================================
float* JMeshData::CreateArrayPtrFloat(std::string fullname,int tag
  ,unsigned size,float* data)
{
  const unsigned sizeok=(tag==TAGDATA12? Npt12: Npt);
  if(size!=sizeok)Run_Exceptioon(fun::PrintStr("Size %u of array \'%s\' does not match the expected value %u.",size,fullname.c_str(),sizeok));
  const string keyname=fun::StrSplitValue(":",fullname,0);
  if(Arrays->ExistsName(keyname) || Arrays->ExistsName(fun::StrLower(keyname)))
    Run_Exceptioon(fun::PrintStr("The array \'%s\' already exists.",keyname.c_str()));
  const unsigned idx=Arrays->CreateArrayFloat(fullname,size,false);
  Arrays->GetArray(idx).tag=tag;
  float* ptr=(float*)Arrays->GetArray(idx).ptr;
  if(data!=NULL)memcpy(ptr,data,sizeof(float)*size);
  return(ptr);
}

//==============================================================================
/// Creates array tfloat3, checks size, copy data and returns memory pointer.
//==============================================================================
tfloat3* JMeshData::CreateArrayPtrFloat3(std::string fullname,int tag
  ,unsigned size,tfloat3* data)
{
  const unsigned sizeok=(tag==TAGDATA12? Npt12: Npt);
  if(size!=sizeok)Run_Exceptioon(fun::PrintStr("Size %u of array \'%s\' does not match the expected value %u.",size,fullname.c_str(),sizeok));
  const string keyname=fun::StrSplitValue(":",fullname,0);
  if(Arrays->ExistsName(keyname) || Arrays->ExistsName(fun::StrLower(keyname)))
    Run_Exceptioon(fun::PrintStr("The array \'%s\' already exists.",keyname.c_str()));
  const unsigned idx=Arrays->CreateArrayFloat3(fullname,size,false);
  Arrays->GetArray(idx).tag=tag;
  tfloat3* ptr=(tfloat3*)Arrays->GetArray(idx).ptr;
  if(data!=NULL)memcpy(ptr,data,sizeof(tfloat3)*size);
  return(ptr);
}

//==============================================================================
/// Creates array double, checks size, copy data and returns memory pointer.
//==============================================================================
double* JMeshData::CreateArrayPtrDouble(std::string fullname,int tag
  ,unsigned size,double* data)
{
  const unsigned sizeok=(tag==TAGDATA12? Npt12: Npt);
  if(size!=sizeok)Run_Exceptioon(fun::PrintStr("Size %u of array \'%s\' does not match the expected value %u.",size,fullname.c_str(),sizeok));
  const string keyname=fun::StrSplitValue(":",fullname,0);
  if(Arrays->ExistsName(keyname) || Arrays->ExistsName(fun::StrLower(keyname)))
    Run_Exceptioon(fun::PrintStr("The array \'%s\' already exists.",keyname.c_str()));
  const unsigned idx=Arrays->CreateArrayDouble(fullname,size,false);
  Arrays->GetArray(idx).tag=tag;
  double* ptr=(double*)Arrays->GetArray(idx).ptr;
  if(data!=NULL)memcpy(ptr,data,sizeof(double)*size);
  return(ptr);
}

//==============================================================================
/// Creates array tdouble3, checks size, copy data and returns memory pointer.
//==============================================================================
tdouble3* JMeshData::CreateArrayPtrDouble3(std::string fullname,int tag
  ,unsigned size,tdouble3* data)
{
  const unsigned sizeok=(tag==TAGDATA12? Npt12: Npt);
  if(size!=sizeok)Run_Exceptioon(fun::PrintStr("Size %u of array \'%s\' does not match the expected value %u.",size,fullname.c_str(),sizeok));
  const string keyname=fun::StrSplitValue(":",fullname,0);
  if(Arrays->ExistsName(keyname) || Arrays->ExistsName(fun::StrLower(keyname)))
    Run_Exceptioon(fun::PrintStr("The array \'%s\' already exists.",keyname.c_str()));
  const unsigned idx=Arrays->CreateArrayDouble3(fullname,size,false);
  Arrays->GetArray(idx).tag=tag;
  tdouble3* ptr=(tdouble3*)Arrays->GetArray(idx).ptr;
  if(data!=NULL)memcpy(ptr,data,sizeof(tdouble3)*size);
  return(ptr);
}

//==============================================================================
/// Returns varlist according to filters.
//==============================================================================
std::string JMeshData::GetVarList(bool novar12,bool var12)const{
  string varlist;
  const JDataArrays* ars=GetArrays();
  const unsigned na=ars->Count();
  for(unsigned ca=0;ca<na;ca++){
    const JDataArrays::StDataArray& ar=ars->GetArrayCte(ca);
    const bool tag12=(ar.tag==TAGDATA12);
    if((novar12 && !tag12) ||(var12 && tag12)){
      if(!varlist.empty())varlist=varlist+",";
      varlist=varlist+ar.keyname;
    }
  }
  return(varlist);
}
//==============================================================================
/// Returns true when number of points of mesh and Arrays structure mmatch.
//==============================================================================
bool JMeshData::EqualStructure(const JMeshData& mdat)const{
  return(Npt1==mdat.Npt1 && Npt2==mdat.Npt2 && Npt3==mdat.Npt3 
    && Arrays->EqualStructure(*(mdat.Arrays)));
}

//==============================================================================
/// Returns true when number of points of mesh and Arrays structure mmatch.
//==============================================================================
bool JMeshData::EqualStructure(const JMeshData* mdat)const{
  return(mdat!=NULL && EqualStructure(*mdat));
}

//==============================================================================
/// Configures mesh definition and data.
//==============================================================================
void JMeshData::ConfigMesh(const jmsh::StMeshPts& mesh,double timestep
  ,std::string varlist)
{
  Reset(true);
  UseMesh=true;
  PtRef=mesh.ptref;
  Vdp1=mesh.vdp1;
  Vdp2=mesh.vdp2;
  Vdp3=mesh.vdp3;
  Npt1=mesh.npt1;
  Npt2=mesh.npt2;
  Npt3=mesh.npt3;
  Npt=Npt1*Npt2*Npt3;
  Npt12=Npt1*Npt2;
  DirDat=mesh.dirdat;
  SetTimeStep(timestep);
  if(!varlist.empty())AddVarList(varlist);
}

//==============================================================================
/// Configures points definition and data.
//==============================================================================
//void JMeshData::ConfigPtos(unsigned npt12,unsigned npt,const tdouble3* pos
void JMeshData::ConfigPtos(unsigned npt,const tdouble3* pos
  ,tfloat3 dirdat,double timestep,std::string varlist)
{
  const unsigned npt12=npt;
  Reset(true);
  UseMesh=false;
  Npt1=npt12;
  Npt2=1;
  Npt3=unsigned(npt/npt12);
  Npt=Npt1*Npt2*Npt3;
  if(Npt!=npt)Run_Exceptioon("The division npt/npt12 is not an integer division.");
  Npt12=Npt1*Npt2;
  DirDat=dirdat;
  if(pos!=NULL){
    PtPos=new tdouble3[Npt];
    memcpy(PtPos,pos,sizeof(tdouble3)*Npt);
  }
  SetTimeStep(timestep);
  if(!varlist.empty())AddVarList(varlist);
}

//==============================================================================
/// Initialises data of arrays to zero.
//==============================================================================
void JMeshData::ClearData(){
  const unsigned nar=Arrays->Count();
  for(unsigned ca=0;ca<nar;ca++){
    JDataArrays::StDataArray& ar=Arrays->GetArray(ca);
    memset(ar.ptr,0,TypeSize(ar.type)*ar.count);
  }
}

//==============================================================================
/// Returns current mesh data definition.
//==============================================================================
StMeshPts JMeshData::GetMeshPt()const{
  StMeshPts m;
  memset(&m,0,sizeof(StMeshPts));
  m.ptref=PtRef;
  m.vdp1=Vdp1;  m.vdp2=Vdp2;  m.vdp3=Vdp3;
  m.npt1=Npt1;  m.npt2=Npt2;  m.npt3=Npt3;  m.npt=Npt;
  m.dirdat=DirDat;
  return(m);
}

//==============================================================================
/// Returns number of variables in Arrays object.
//==============================================================================
unsigned JMeshData::GetVarCount()const{
  return(Arrays->Count());
}

//==============================================================================
/// Returns name of requested variable in Arrays object.
//==============================================================================
std::string JMeshData::GetVarName(unsigned idx)const{
  if(idx>=Arrays->Count())Run_Exceptioon("Array idx is invalid.");
  return(Arrays->GetArrayCte(idx).keyname);
}

//==============================================================================
/// Returns tag of requested variable in Arrays object.
//==============================================================================
int JMeshData::GetVarTag(unsigned idx)const{
  if(idx>=Arrays->Count())Run_Exceptioon("Array idx is invalid.");
  return(Arrays->GetArrayCte(idx).tag);
}

//==============================================================================
/// Returns requested pointer data of variable.
//==============================================================================
tfloat3* JMeshData::GetVarFloat3(const std::string& name){
  const unsigned idx=Arrays->GetIdxName(name);
  return(idx==UINT_MAX? NULL: (tfloat3*)Arrays->GetArrayCte(idx).ptr);
}

//==============================================================================
/// Returns requested pointer data of variable.
//==============================================================================
float* JMeshData::GetVarFloat(const std::string& name){
  const unsigned idx=Arrays->GetIdxName(name);
  return(idx==UINT_MAX? NULL: (float*)Arrays->GetArrayCte(idx).ptr);
}

//==============================================================================
/// Stores positions X,Y,Z with single precision.
//==============================================================================
void JMeshData::GetPosf(unsigned np,tfloat3* pos)const{
  //printf("==> np:%u  Npt:%u  Npt12:%u \n",np,Npt,Npt12);
  if(np!=Npt && np!=Npt12)Run_Exceptioon("Number of points does not match.");
  const bool data12=(np!=Npt && np==Npt12);
  if(UseMesh)GetMeshPos(GetMeshPt(),(data12? Npt12: Npt),pos);
  else{
    const unsigned npt=(data12? Npt12: Npt);
    if(PtPos)for(unsigned cp=0;cp<npt;cp++)pos[cp]=ToTFloat3(PtPos[cp]);
    else memset(pos,0,sizeof(tfloat3)*npt);
  }
}

//==============================================================================
/// Stores positions X,Y,Z with double precision.
//==============================================================================
void JMeshData::GetPosd(unsigned np,tdouble3* pos)const{
  if(np!=Npt && np!=Npt12)Run_Exceptioon("Number of points does not match.");
  const bool data12=(np!=Npt && np==Npt12);
  if(UseMesh)GetMeshPos(GetMeshPt(),(data12? Npt12: Npt),pos);
  else{
    const unsigned npt=(data12? Npt12: Npt);
    if(PtPos)memcpy(pos,PtPos,sizeof(tdouble3)*npt);
    else memset(pos,0,sizeof(tdouble3)*npt);
  }
}

//==============================================================================
/// Stores positions X,Y,Z according to a mesh definition (single precision).
//==============================================================================
void JMeshData::GetMeshPos(const StMeshPts& m,unsigned np,tfloat3* pos){
  const unsigned npt12=m.npt1*m.npt2;
  //printf("==> np:%u  npt:%u  npt12:%u \n",np,m.npt,npt12);
  if(np!=m.npt && np!=npt12)Run_ExceptioonSta("Number of points does not match.");
  const bool data12=(np!=m.npt && np==npt12);
  const unsigned np3=(data12? 1: m.npt3);
  unsigned cp=0;
  for(unsigned cp3=0;cp3<np3;cp3++){
    const tdouble3 pt3=m.ptref+(m.vdp3*cp3);
    for(unsigned cp2=0;cp2<m.npt2;cp2++){
      const tdouble3 pt2=pt3+(m.vdp2*cp2);
      for(unsigned cp1=0;cp1<m.npt1;cp1++,cp++){
        pos[cp]=ToTFloat3(pt2+(m.vdp1*cp1));
      }
    }
  }
}

//==============================================================================
/// Stores positions X,Y,Z according to a mesh definition (single precision).
//==============================================================================
void JMeshData::GetMeshPos(const StMeshPts& m,unsigned np,tdouble3* pos){
  const unsigned npt12=m.npt1*m.npt2;
  //printf("==> np:%u  npt:%u  npt12:%u \n",np,m.npt,npt12);
  if(np!=m.npt && np!=npt12)Run_ExceptioonSta("Number of points does not match.");
  const bool data12=(np!=m.npt && np==npt12);
  const unsigned np3=(data12? 1: m.npt3);
  unsigned cp=0;
  for(unsigned cp3=0;cp3<np3;cp3++){
    const tdouble3 pt3=m.ptref+(m.vdp3*cp3);
    for(unsigned cp2=0;cp2<m.npt2;cp2++){
      const tdouble3 pt2=pt3+(m.vdp2*cp2);
      for(unsigned cp1=0;cp1<m.npt1;cp1++,cp++){
        pos[cp]=(pt2+(m.vdp1*cp1));
      }
    }
  }
}

//==============================================================================
/// Returns StGridPts data according to StGridBasic configuration.
//==============================================================================
StMeshPts JMeshData::MakeMeshPt(const StMeshBasic& b){
  StMeshPts m;
  memset(&m,0,sizeof(StMeshPts));
  m.ptref=b.ptref;
  m.vdp1=m.vdp2=m.vdp3=TDouble3(0);
  m.npt1=m.npt2=m.npt3=1;
  tdouble3 vec;
  double   dis,dispt;
  //-First direction.
  vec=b.vec1; dis=b.dis1; dispt=b.dispt1;
  if(vec!=TDouble3(0) && dispt>0 && dis>=dispt){
    m.vdp1=fgeo::VecUnitary(vec)*dispt;
    m.npt1=unsigned((dis+dispt*0.01)/dispt)+1;
  }
  //-Second direction.
  vec=b.vec2; dis=b.dis2; dispt=b.dispt2;
  if(vec!=TDouble3(0) && dispt>0 && dis>=dispt){
    m.vdp2=fgeo::VecUnitary(vec)*dispt;
    m.npt2=unsigned((dis+dispt*0.01)/dispt)+1;
  }
  //-Third direction.
  vec=b.vec3; dis=b.dis3; dispt=b.dispt3;
  if(vec!=TDouble3(0) && dispt>0 && dis>=dispt){
    m.vdp3=fgeo::VecUnitary(vec)*dispt;
    m.npt3=unsigned((dis+dispt*0.01)/dispt)+1;
  }
  m.npt=m.npt1*m.npt2*m.npt3;
  m.dirdat=fgeo::VecUnitary(b.dirdat);
  return(m);
}

//==============================================================================
/// Displays values of StGridPts data.
//==============================================================================
void JMeshData::PrintMeshPts(std::string text,const StMeshPts& m){
  printf("%s\n",text.c_str());
  printf(" npt.:%4u  ptref:(%s)\n",m.npt,fun::Double3gStr(m.ptref).c_str());
  printf(" npt1:%4u  vdp1.:(%s)\n",m.npt1,fun::Double3gStr(m.vdp1).c_str());
  printf(" npt2:%4u  vdp2.:(%s)\n",m.npt2,fun::Double3gStr(m.vdp2).c_str());
  printf(" npt3:%4u  vdp3.:(%s)\n",m.npt3,fun::Double3gStr(m.vdp3).c_str());
  printf(" dirdat:(%s)\n\n",fun::Float3gStr(m.dirdat).c_str());
}


}


