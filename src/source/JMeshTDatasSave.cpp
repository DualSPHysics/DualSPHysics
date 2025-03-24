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

/// \file JMeshTDatasSave.cpp \brief Implements the classes JMeshTDatasSave.

#include "JMeshTDatasSave.h"
#include "Functions.h"
#include "JMeshData.h"
#include "JSpVtkData.h"
#include "JDataArrays.h"
#include "JSpVtkShape.h"
#include "JSaveCsv2.h"
#include "JException.h"

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshTDatasSave
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMeshTDatasSave::JMeshTDatasSave(){
  ClassName="JMeshTDatasSave";
  Data=NULL;  DataTime=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JMeshTDatasSave::~JMeshTDatasSave(){
  DestructorActive=true;
  Reset();
  delete Data; Data=NULL;
}

//==============================================================================
/// Throws exception related to a file from a static method.
//==============================================================================
void JMeshTDatasSave::RunExceptioonStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,const std::string& msg,const std::string& file)
{
  throw JException(srcfile,srcline,"JMeshTDatasSave",method,msg,file);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JMeshTDatasSave::Reset(){
  ResetData();
  FormatVer=FormatVerDef;
  memset(&Mesh,0,sizeof(StMeshPts));
  Npt=Npt12=0;
  Datdefs.clear();
  InitialSaved=false;
}

//==============================================================================
/// Elimina informacion de Data.
/// Deletes information from data. 
//==============================================================================
void JMeshTDatasSave::ResetData(){
  delete Data; 
  Data=new JBinaryData("JMeshTDatas");
  DataTime=Data->CreateItem("DataTime");
  Ctime=0;
}

//==============================================================================
/// Elimina informacion del instante actual.
/// Deletes information of current time.
//==============================================================================
void JMeshTDatasSave::ResetDataTime(){
  DataTime->Clear();
}

//==============================================================================
/// Devuelve la memoria reservada.
/// Returns allocated memory
//==============================================================================
llong JMeshTDatasSave::GetAllocMemory()const{  
  return(Data->GetAllocMemory());
}

//==============================================================================
/// Devuelve nombre de instante segun su numero.
/// Returns name of datatime according to their number.
//==============================================================================
std::string JMeshTDatasSave::GetNameDataTime(unsigned ctime){
  return(fun::PrintStr("TIME_%06u",ctime));
}

//==============================================================================
/// Configuracion de datos de cabecera.
/// Configuration of header data.
//==============================================================================
void JMeshTDatasSave::Config(const std::string savefile,std::string appname
  ,const JMeshData* mdat)
{
  Reset();
  AppName=appname;
  SaveFile=savefile;
  Mesh=mdat->GetMeshPt();
  if(Mesh.npt>1 && !CheckMeshVectors(Mesh))Run_Exceptioon("Mesh data without mesh definition is invalid for binary files.");
  Npt12=Mesh.npt1*Mesh.npt2;
  Npt=Npt12*Mesh.npt3;
}

//==============================================================================
/// Grabacion inicial de fichero con info de Data.
/// Initial file recording with info from Data.
//==============================================================================
void JMeshTDatasSave::SaveInitial(){
  if(!InitialSaved){
    Data->SetvText   ("AppName",AppName);
    Data->SetvUint   ("FormatVer",FormatVer);
    Data->SetvDouble3("MeshPtRef" ,Mesh.ptref);
    Data->SetvDouble3("MeshVdp1"  ,Mesh.vdp1);
    Data->SetvDouble3("MeshVdp2"  ,Mesh.vdp2);
    Data->SetvDouble3("MeshVdp3"  ,Mesh.vdp3);
    Data->SetvUint   ("MeshNpt1"  ,Mesh.npt1);
    Data->SetvUint   ("MeshNpt2"  ,Mesh.npt2);
    Data->SetvUint   ("MeshNpt3"  ,Mesh.npt3);
    Data->SetvFloat3 ("MeshDirDat",Mesh.dirdat);
    const unsigned nd=unsigned(Datdefs.size());
    Data->SetvUint   ("Data_Count",nd);
    for(unsigned c=0;c<nd;c++){
      const StMeshData& dd=Datdefs[c];
      Data->SetvText(fun::PrintStr("DataName_%02u",c),dd.name);
      Data->SetvUint(fun::PrintStr("DataType_%02u",c),unsigned(dd.type));
      Data->SetvUint(fun::PrintStr("DataSize_%02u",c),dd.size);
      Data->SetvInt (fun::PrintStr("DataTag_%02u" ,c),dd.tag);
    }
    DataTime->SetHide(true);
    Data->SaveFile(SaveFile,true,false);
    InitialSaved=true;
  }
}

//==============================================================================
/// Checks number of points and data definitions according to previous data.
//==============================================================================
void JMeshTDatasSave::CheckDataDefinitions(const JMeshData* mdat){
  //-Checks size of data.
  if(Mesh.npt1!=mdat->GetNpt1() || Mesh.npt2!=mdat->GetNpt2() || Mesh.npt3!=mdat->GetNpt3())
    Run_ExceptioonFile("Number of points does not match.",SaveFile);
  //-Checks data definitions or creates it.
  if(Datdefs.empty()){
    const JDataArrays* ars=mdat->GetArrays();
    const unsigned na=ars->Count();
    for(unsigned ca=0;ca<na;ca++){
      const JDataArrays::StDataArray& ar=ars->GetArrayCte(ca);
      StMeshData dd={ar.keyname,ar.type,ar.count,ar.tag};
      Datdefs.push_back(dd);
    }
  }
  else{
    bool eq=(mdat->GetArrays()->Count()==unsigned(Datdefs.size()));
    const JDataArrays* ars=mdat->GetArrays();
    const unsigned na=ars->Count();
    for(unsigned ca=0;ca<na && eq;ca++){
      const JDataArrays::StDataArray& ar=ars->GetArrayCte(ca);
      const StMeshData& dd=Datdefs[ca];
      eq=(ar.keyname==dd.name && ar.type==dd.type && ar.count==dd.size && ar.tag==dd.tag);
    }
    if(!eq)Run_Exceptioon("Data definitions do not match the previous data definitions.");
  }
}

//==============================================================================
/// Incorpora datos de particulas de de nuevo part.
/// Adds data of particles to new part.
//==============================================================================
void JMeshTDatasSave::SaveDataTime(const JMeshData* mdat){
  //-Configures item DataTime.
  DataTime->Clear();
  DataTime->SetName(GetNameDataTime(Ctime));
  DataTime->SetvUint("Ctime",Ctime);
  DataTime->SetvDouble("TimeStep",mdat->GetTimeStep());
  //-Checks number of points and data definitions according to previous data.
  CheckDataDefinitions(mdat);

  //-Crea arrays con datos en DataTime object.
  const bool externalpointer=true;
  const JDataArrays* ars=mdat->GetArrays();
  const unsigned na=ars->Count();
  for(unsigned ca=0;ca<na;ca++){
    const JDataArrays::StDataArray& ar=ars->GetArrayCte(ca);
    JBinaryDataDef::TpData dtype=JBinaryDataDef::DatNull;
    switch(ar.type){
      case TypeFloat:    dtype=JBinaryDataDef::DatFloat;    break;
      case TypeDouble:   dtype=JBinaryDataDef::DatDouble;   break;
      case TypeFloat3:   dtype=JBinaryDataDef::DatFloat3;   break;
      case TypeDouble3:  dtype=JBinaryDataDef::DatDouble3;  break;
      default: Run_Exceptioon(fun::PrintStr("Type of pointer \'%s\' is invalid.",TypeToStr(ar.type)));
    }
    DataTime->CreateArray(ar.keyname,dtype,ar.count,ar.ptr,externalpointer);
  }
  //-Save file data.
  if(!InitialSaved)SaveInitial();
  DataTime->SaveFileListApp(SaveFile,"JMeshTDatas",true,true);
  DataTime->RemoveArrays();
  Ctime++;
}

//==============================================================================
/// Incorpora datos de particulas de de nuevo part.
/// Adds data of particles to new part.
//==============================================================================
void JMeshTDatasSave::SaveDataTimes(unsigned nt,JMeshData** vmdat){
  //-Configures item DataTime.
  DataTime->Clear();
  DataTime->SetName(GetNameDataTime(Ctime));
  DataTime->SetvUint("Ctime",Ctime);
  DataTime->SetvUint("Ntimes",nt);
  //DataTime->SetvDouble("TimeStep",mdat->GetTimeStep());
  double* vtimes=new double[nt];
  for(unsigned ct=0;ct<nt;ct++)vtimes[ct]=vmdat[ct]->GetTimeStep();
  DataTime->CreateArray("TimeSteps",JBinaryDataDef::DatDouble,nt,vtimes,true);
  //-Checks number of points and data definitions according to previous data.
  for(unsigned ct=0;ct<nt;ct++)CheckDataDefinitions(vmdat[ct]);
  //-Crea arrays con datos en DataTime object.
  const bool externalpointer=true;
  const JDataArrays* ars=vmdat[0]->GetArrays();
  const unsigned na=ars->Count();
  for(unsigned ca=0;ca<na;ca++){
    const JDataArrays::StDataArray& ar=ars->GetArrayCte(ca);
    JBinaryDataDef::TpData dtype=JBinaryDataDef::DatNull;
    switch(ar.type){
      case TypeFloat:    dtype=JBinaryDataDef::DatFloat;    break;
      case TypeDouble:   dtype=JBinaryDataDef::DatDouble;   break;
      case TypeFloat3:   dtype=JBinaryDataDef::DatFloat3;   break;
      case TypeDouble3:  dtype=JBinaryDataDef::DatDouble3;  break;
      default: Run_Exceptioon(fun::PrintStr("Type of pointer \'%s\' is invalid.",TypeToStr(ar.type)));
    }
    JBinaryDataArray* ar2=DataTime->CreateArray(ar.keyname,dtype);
    ar2->AllocMemory(ar.count*nt,false);
    for(unsigned ct=0;ct<nt;ct++){
      const JDataArrays* ars_ct=vmdat[ct]->GetArrays();
      const JDataArrays::StDataArray& ar_ct=ars_ct->GetArrayCte(ca);
      ar2->AddData(ar_ct.count,ar_ct.ptr,false);
    }
  }
  //-Save file data.
  if(!InitialSaved)SaveInitial();
  DataTime->SaveFileListApp(SaveFile,"JMeshTDatas",true,true);
  DataTime->RemoveArrays();
  Ctime++;
  //-Free memory.
  delete[] vtimes; vtimes=NULL;
}

//==============================================================================
/// Saves data in CSV.
//==============================================================================
void JMeshTDatasSave::SaveCsvHead(const StMeshPts& m,std::string name,std::string tunits
  ,TpTypeData type,bool data12,jcsv::JSaveCsv2& scsv)
{
  const bool usegrid=CheckMeshVectors(m);
  scsv.SetHead();
  //-Head definition.
  scsv << "Format;DataName;DataUnits;DataType;Data12;Npt1;Npt2;Npt3";
  if(usegrid)scsv << "PtRef.x [m];PtRef.y;PtRef.z";
  if(usegrid)scsv << "Vec1.x [m];Vec1.y;Vec1.z;Vec2.x;Vec2.y;Vec2.z;Vec3.x;Vec3.y;Vec3.z";
  scsv << "DirDat.x;DirDat.y;DirDat.z";
  scsv << jcsv::Endl();
  const string fmt=(usegrid? fun::PrintStr("JMeshTDatas-%u",FormatVerDef): string("SimpleCSV"));
  scsv << fmt << name << tunits << TypeToStr(type) << (data12? "true": "false") << m.npt1 << m.npt2 << m.npt3;
  if(usegrid)scsv << m.ptref << m.vdp1 << m.vdp2 << m.vdp3;
  scsv << m.dirdat;
  scsv << jcsv::Endl();
  //-Data headers.
  const int tdim=TypeDim(type);
  if(tdim!=1 && tdim!=3)Run_ExceptioonFileSta("Dimension of type of data is invalid.",scsv.GetFileName());
  scsv << "time [s]";
  if(tdim==1){
    if(data12)for(unsigned c2=0;c2<m.npt2;c2++)for(unsigned c1=0;c1<m.npt1;c1++)
      scsv << fun::PrintStr("v(%d:%d)",c1,c2);
    else for(unsigned c3=0;c3<m.npt3;c3++)for(unsigned c2=0;c2<m.npt2;c2++)for(unsigned c1=0;c1<m.npt1;c1++)
      scsv << fun::PrintStr("v(%d:%d:%d)",c1,c2,c3);
  }
  else{
    if(data12)for(unsigned c2=0;c2<m.npt2;c2++)for(unsigned c1=0;c1<m.npt1;c1++)for(int d=1;d<=tdim;d++)
      scsv << fun::PrintStr("v(%d:%d).%c",c1,c2,(d==1? 'x': (d==2? 'y': 'z')));
    else for(unsigned c3=0;c3<m.npt3;c3++)for(unsigned c2=0;c2<m.npt2;c2++)for(unsigned c1=0;c1<m.npt1;c1++)for(int d=1;d<=tdim;d++)
      scsv << fun::PrintStr("v(%d:%d:%d).%c",c1,c2,c3,(d==1? 'x': (d==2? 'y': 'z')));
  }
  scsv << jcsv::Endl();
}

//==============================================================================
/// Saves data in CSV.
//==============================================================================
void JMeshTDatasSave::SaveCsv(const JMeshData* mdat,unsigned cdata
  ,jcsv::JSaveCsv2& scsv,bool svhead)
{
  if(cdata>=mdat->GetArrays()->Count())Run_ExceptioonFileSta("Selected array does not exists.",scsv.GetFileName());
  const StMeshPts m=mdat->GetMeshPt();
  const JDataArrays::StDataArray& ar=mdat->GetArrays()->GetArrayCte(cdata);
  const bool data12=(ar.tag==JMeshData::TAGDATA12);
  //-Output format configuration.
  if(1){//-Full string format.
    scsv << jcsv::Fmt(jcsv::TpFloat1,"%15.7E") << jcsv::Fmt(jcsv::TpFloat3,"%15.7E;%15.7E;%15.7E");
    scsv << jcsv::Fmt(jcsv::TpDouble1,"%20.12E") << jcsv::Fmt(jcsv::TpDouble3,"%20.12E;%20.12E;%20.12E");
  }
  else scsv << jcsv::Fmt(jcsv::TpFloat1,"%g") << jcsv::Fmt(jcsv::TpFloat3,"%g;%g;%g");
  //-Saves head data.
  if(svhead)SaveCsvHead(m,ar.keyname,mdat->GetArrays()->GetArrayUnits(cdata)
    ,ar.type,data12,scsv);
  //-Saves current data.
  scsv.SetData();
  scsv << mdat->GetTimeStep();
  const unsigned npt1=m.npt1;
  const unsigned npt2=m.npt2;
  const unsigned npt3=(data12? 1: m.npt3);
  switch(ar.type){
    case TypeFloat:{
      const float* pdat=(float*)ar.ptr;
      for(unsigned c3=0;c3<npt3;c3++)for(unsigned c2=0;c2<npt2;c2++)for(unsigned c1=0;c1<npt1;c1++)
        scsv << *(pdat++);
    }break;
    case TypeFloat3:{
      const tfloat3* pdat=(tfloat3*)ar.ptr;
      for(unsigned c3=0;c3<npt3;c3++)for(unsigned c2=0;c2<npt2;c2++)for(unsigned c1=0;c1<npt1;c1++)
        scsv << *(pdat++);
    }break;
    case TypeDouble:{
      const double* pdat=(double*)ar.ptr;
      for(unsigned c3=0;c3<npt3;c3++)for(unsigned c2=0;c2<npt2;c2++)for(unsigned c1=0;c1<npt1;c1++)
        scsv << *(pdat++);
    }break;
    case TypeDouble3:{
      const tdouble3* pdat=(tdouble3*)ar.ptr;
      for(unsigned c3=0;c3<npt3;c3++)for(unsigned c2=0;c2<npt2;c2++)for(unsigned c1=0;c1<npt1;c1++)
        scsv << *(pdat++);
    }break;
    default: Run_ExceptioonFileSta(fun::PrintStr("Type of pointer \'%s\' is invalid.",TypeToStr(ar.type)),scsv.GetFileName());
  }
  scsv << jcsv::Endl();
}

//==============================================================================
/// Saves VTK with data.
//==============================================================================
void JMeshTDatasSave::SaveVtk(std::string file,int fnum,const JMeshData* mdat
  ,bool svdata12)
{
  if(fun::GetExtension(file).empty())file=fun::AddExtension(file,"vtk");
  if(fnum>=0)file=fun::FileNameSec(file,fnum);
  //-Loads data in arrays.
  JDataArrays arrays;
  const StMeshPts m=mdat->GetMeshPt();
  //JMeshData::PrintMeshPts("JMeshTDatasSave::SaveVtk\n",m);
  const unsigned np=(svdata12? m.npt1*m.npt2: m.npt);
  //printf("--> np:%u\n",np);
  const JDataArrays* ars=mdat->GetArrays();
  const unsigned na=ars->Count();
  for(unsigned ca=0;ca<na;ca++){
    const JDataArrays::StDataArray& ar=ars->GetArrayCte(ca);
    if(svdata12==(ar.tag==JMeshData::TAGDATA12)){
      if(ar.count!=np)Run_ExceptioonFileSta(fun::PrintStr("Size of array \'%s\' is invalid.",ar.keyname.c_str()),file);
      arrays.AddArray(ar.keyname,ar.type,ar.count,ar.ptr,false);
    }
  }
  //-Loads positions.
  if(arrays.Count()){
    tfloat3* pos=arrays.CreateArrayPtrFloat3("Pos",np,false);
    mdat->GetPosf(np,pos);
    //-Updates z of position when Zsurf data is available.
    const unsigned idx=arrays.GetIdxName("Zsurf");
    if(idx!=UINT_MAX){
      const float* zsurf=arrays.GetArrayFloat(idx);
      for(unsigned p=0;p<np;p++)pos[p].z=zsurf[p];
    }
    //-Saves VTK file.
    JSpVtkData::Save(file,arrays,"Pos");
  }
}

//==============================================================================
/// Creates VTK file with the scheme of mesh definition.
//==============================================================================
void JMeshTDatasSave::SaveVtkScheme(std::string file,StMeshPts m){
  //JMeshData::PrintMeshPts("SaveVtkSchem>\n",m);
  //if(!CheckMeshVectors(m))Run_ExceptioonFileSta("Mesh definition invalid to create scheme.",file);
  JSpVtkShape ss;
  if(CheckMeshVectors(m)){
    const int ndef=(m.npt1>1? 1: 0)+(m.npt2>1? 1: 0)+(m.npt3>1? 1: 0);
    if(!ndef)ss.AddPoint(m.ptref);
    else if(ndef==1){//-Grid-Line.
      if(m.npt1>1)ss.AddLine(m.ptref,m.ptref+(m.vdp1*(m.npt1-1)));
      if(m.npt2>1)ss.AddLine(m.ptref,m.ptref+(m.vdp2*(m.npt2-1)));
      if(m.npt3>1)ss.AddLine(m.ptref,m.ptref+(m.vdp3*(m.npt3-1)));
    }
    else if(ndef==2){//-Grid-Plane.
      if(m.npt1>1 && m.npt2>1){
        const unsigned n1=m.npt1-1,n2=m.npt2-1;
        const tdouble3 v1=m.vdp1,v2=m.vdp2;
        const tdouble3 p0=m.ptref;
        const tdouble3 p1=p0+(v1*n1);
        const tdouble3 p2=p1+(v2*n2);
        const tdouble3 p3=p0+(v2*n2);
        ss.AddQuad(p0,p1,p2,p3);
      }
      else if(m.npt1>1 && m.npt3>1){
        const unsigned n1=m.npt1-1,n2=m.npt3-1;
        const tdouble3 v1=m.vdp1,v2=m.vdp3;
        const tdouble3 p0=m.ptref;
        const tdouble3 p1=p0+(v1*n1);
        const tdouble3 p2=p1+(v2*n2);
        const tdouble3 p3=p0+(v2*n2);
        ss.AddQuad(p0,p1,p2,p3);
      }
      else{//-Grid-Box.
        const unsigned n1=m.npt2-1,n2=m.npt3-1;
        const tdouble3 v1=m.vdp2,v2=m.vdp3;
        const tdouble3 p0=m.ptref;
        const tdouble3 p1=p0+(v1*n1);
        const tdouble3 p2=p1+(v2*n2);
        const tdouble3 p3=p0+(v2*n2);
        ss.AddQuad(p0,p1,p2,p3);
      }
    }
    else if(ndef==3){
      ss.AddBoxSizeVec(m.ptref,m.vdp1*(m.npt1-1),m.vdp2*(m.npt2-1),m.vdp3*(m.npt3-1));
    }
  }
  else ss.AddPoint(m.ptref);
  ss.SaveVtk(file,"num");
}

}
