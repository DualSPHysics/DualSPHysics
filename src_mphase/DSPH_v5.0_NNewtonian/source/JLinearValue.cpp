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

/// \file JLinearValue.cpp \brief Implements the class \ref JLinearValue.

#include "JLinearValue.h"
#include "Functions.h"
#include "JReadDatafile.h"
#ifdef JLinearValue_UseJXml
  #include "JXml.h"
#endif
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JLinearValue
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JLinearValue::JLinearValue(unsigned nvalues,bool specialvalues,bool optionalvalues)
  :Nvalues(max(1u,nvalues)),SpecialValues(specialvalues),OptionalValues(optionalvalues)
{
  ClassName="JLinearValue";
  Times=NULL;
  Values=NULL;
  Reset();
  if(SpecialValues && OptionalValues)
    Run_Exceptioon("Options SpecialValues and OptionalValues at the same time are invalid.");
}

//==============================================================================
/// Constructor with input file.
//==============================================================================
JLinearValue::JLinearValue(const std::string &inputfile,unsigned nvalues
  ,bool specialvalues,bool optionalvalues)
  :Nvalues(max(1u,nvalues)),SpecialValues(specialvalues),OptionalValues(optionalvalues)
{
  ClassName="JLinearValue";
  Times=NULL;
  Values=NULL;
  Reset();
  File=inputfile;
  if(SpecialValues && OptionalValues)
    Run_ExceptioonFile("Options SpecialValues and OptionalValues at the same time are invalid.",inputfile);
}

//==============================================================================
/// Constructor for copy.
//==============================================================================
JLinearValue::JLinearValue(const JLinearValue &obj)
  :Nvalues(obj.Nvalues),SpecialValues(obj.SpecialValues),OptionalValues(obj.OptionalValues)
{
  ClassName="JLinearValue";
  Times=NULL;
  Values=NULL;
  CopyFrom(obj);
  if(SpecialValues && OptionalValues)
    Run_Exceptioon("Options SpecialValues and OptionalValues at the same time are invalid.");
}

//==============================================================================
/// Destructor.
//==============================================================================
JLinearValue::~JLinearValue(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JLinearValue::Reset(){
  SetSize(0);
  File="";
  NewInterval=false;
  TimeStep=TimePre=TimeNext=TimeFactor=0;
}

//==============================================================================
/// Copy data from other object.
//==============================================================================
void JLinearValue::CopyFrom(const JLinearValue &obj){
  if(Nvalues!=obj.Nvalues)Run_Exceptioon("The copy is invalid since Nvalues does not match.");
  if(SpecialValues!=obj.SpecialValues)Run_Exceptioon("The copy is invalid since SpecialValues does not match.");
  if(OptionalValues!=obj.OptionalValues)Run_Exceptioon("The copy is invalid since OptionalValues does not match.");
  Reset();
  File=obj.File;
  SetSize(obj.Size);
  memcpy(Times,obj.Times,sizeof(double)*Size);
  memcpy(Values,obj.Values,sizeof(double)*Nvalues*Size);
  Count        =obj.Count;
  NewInterval  =obj.NewInterval;
  TimeStep     =obj.TimeStep;
  Position     =obj.Position;
  PositionNext =obj.PositionNext;
  TimePre      =obj.TimePre;
  TimeNext     =obj.TimeNext;
  TimeFactor   =obj.TimeFactor;
}

//==============================================================================
/// Returns the allocated memory.
/// Devuelve la memoria reservada.
//==============================================================================
unsigned JLinearValue::GetAllocMemory()const{
  unsigned s=0;
  if(Times)s+=sizeof(double)*Size;
  if(Values)s+=sizeof(double)*Nvalues*Size;
  return(s);
}

//==============================================================================
/// Sets the indicated size to maintain the content.
/// Ajusta al tamano indicado manteniendo el contenido.
//==============================================================================
void JLinearValue::SetSize(unsigned size){
  if(size>=SIZEMAX)Run_Exceptioon("It has reached the maximum size allowed.");
  Size=size;
  if(Count>Size)Count=Size;
  if(Size){
    Times=fun::ResizeAlloc(Times,Count,size);
    Values=fun::ResizeAlloc(Values,Nvalues*Count,Nvalues*size);
  }
  else{
    delete[] Times;  Times=NULL;
    delete[] Values; Values=NULL;
  }
  Position=PositionNext=UINT_MAX;
}

//==============================================================================
/// Modify times and values.
/// Modifica tiempos y valores.
//==============================================================================
void JLinearValue::SetTimeValue(unsigned idx,double time,double value){
  if(idx>=Count)Run_Exceptioon("Index is invalid.");
  //printf("\n=====>> SET t[%d]:%f\n",idx,time);
  //if(idx>0)printf("=====>> t-1[%d]:%g\n",idx-1,Times[idx-1]);
  //if(true) printf("=====>> t  [%d]:%g\n",idx,Times[idx]);
  //if(idx+1<Count)printf("=====>> t+1[%d]:%g\n",idx+1,Times[idx+1]);
  if((idx>0 && Times[idx-1]>time) || (idx+1<Count && Times[idx+1]<time))
    Run_Exceptioon("New time values is lower than previous time or it is higher than next time.");
  Times[idx]=time;
  Values[Nvalues*idx]=value;
  for(unsigned cv=1;cv<Nvalues;cv++)Values[Nvalues*idx+cv]=0;
  //if(value==DBL_MAX)SpecialValues=true;
}

//==============================================================================
/// Adds values at the end of the list.
/// Anhade valores al final de la lista.
//==============================================================================
unsigned JLinearValue::AddTimeValue(double time,double value){
  if(Count==Size)SetSize(Size+SIZEINITIAL);
  const unsigned idx=Count;
  Times[idx]=time;
  Values[Nvalues*idx]=value;
  for(unsigned cv=1;cv<Nvalues;cv++)Values[Nvalues*idx+cv]=0;
  Count++;
  //if(value==DBL_MAX)SpecialValues=true;
  return(idx);
}

//==============================================================================
/// Modify value at the indicated position.
/// Modifica valor en la posicion indicada.
//==============================================================================
void JLinearValue::SetValue(unsigned idx,unsigned cvalue,double value){
  if(idx>=Count)Run_Exceptioon("Idx is not valid.");
  if(cvalue>=Nvalues)Run_Exceptioon("Number of value is not valid.");
  Values[Nvalues*idx+cvalue]=value;
  //if(value==DBL_MAX)SpecialValues=true;
}

//==============================================================================
/// Finds time before and after indicated time.
/// Busca time anterior y posterior del instante indicado.
//==============================================================================
void JLinearValue::FindTime(double timestep){
  if(!Count)Run_Exceptioon("There are not times.");
  NewInterval=false;
  if(timestep!=TimeStep || Position==UINT_MAX){
    unsigned pos=(Position==UINT_MAX? 0: Position);
    unsigned posnext=pos;
    double tpre=Times[pos];
    double tnext=tpre;
    if(Count>1){
      while(tpre>=timestep && pos>0){//-Retrocede.
        pos--;
        tpre=Times[pos];
      }
      posnext=(pos+1<Count? pos+1: pos);
      tnext=Times[posnext];
      while(tnext<timestep && posnext+1<Count){//-Avanza.
        posnext++;
        tnext=Times[posnext];
      }
      if(posnext-pos>1){
        pos=posnext-1;
        tpre=Times[pos];
      }
    }
    NewInterval=(Position!=pos || (timestep>tnext && TimeStep<=tnext));
    TimeStep=timestep;
    Position=pos;
    PositionNext=posnext;
    TimePre=tpre;
    TimeNext=tnext;
    const double tdif=TimeNext-TimePre;
    TimeFactor=(tdif? (timestep-TimePre)/tdif: 0);
  }
}

//==============================================================================
/// Returns the interpolated value for the time indicated.
/// If no values always returns 0.
/// If only one value always returns that value.
/// If the indicated t is less than the minimum returns the first value.
/// If the indicated t is greater than the maximum returns the last value.
/// 
/// Devuelve valor el valor interpolado para el instante indicado.
/// Si no hay valores siempre devuelve 0.
/// Si solo hay un valor siempre devuelve ese valor.
/// Si el t indicado es menor que el minimo devuelve el primer valor.
/// Si el t indicado es mayor que el maximo devuelve el ultimo valor.
//==============================================================================
double JLinearValue::GetValue(double timestep,unsigned cvalue){
  double ret=0;
  FindTime(timestep);
  //printf("--> t:%f  [%u - %u]  [%f - %f]\n",timestep,Position,PositionNext,TimePre,TimeNext);
  if(timestep<=TimePre)ret=Values[Nvalues*Position+cvalue];
  else if(timestep>=TimeNext)ret=Values[Nvalues*PositionNext+cvalue];
  else{
    const double vini=Values[Nvalues*Position+cvalue];
    const double vnext=Values[Nvalues*PositionNext+cvalue];
    ret=(TimeFactor*(vnext-vini)+vini);
    if(SpecialValues){
      if(vini==DBL_MAX)ret=DBL_MAX;
      else if(vnext==DBL_MAX)ret=vini;
    }
  }
  return(ret);
}

//==============================================================================
/// Returns result of GetValue() converted to float according special values.
//==============================================================================
float JLinearValue::GetValuef(double timestep,unsigned cvalue){
  const double v=GetValue(timestep,cvalue);
  return(v==DBL_MAX? FLT_MAX: float(v));
}

//==============================================================================
/// Returns the interpolated values for the time indicated.
/// If no values always returns 0.
/// If only one value always returns that value.
/// If the indicated t is less than the minimum returns the first value.
/// If the indicated t is greater than the maximum returns the last value.
/// 
/// Devuelve valor el valor interpolado para el instante indicado.
/// Si no hay valores siempre devuelve 0.
/// Si solo hay un valor siempre devuelve ese valor.
/// Si el t indicado es menor que el minimo devuelve el primer valor.
/// Si el t indicado es mayor que el maximo devuelve el ultimo valor.
//==============================================================================
tdouble3 JLinearValue::GetValue3d(double timestep){
  tdouble3 ret=TDouble3(0);
  FindTime(timestep);
  //printf("--> t:%f  [%u - %u]  [%f - %f]\n",timestep,Position,PositionNext,TimePre,TimeNext);
  if(timestep<=TimePre){
    const unsigned rpos=Nvalues*Position;
    ret=TDouble3(Values[rpos],Values[rpos+1],Values[rpos+2]);
  }
  else if(timestep>=TimeNext){
    const unsigned rpos=Nvalues*PositionNext;
    ret=TDouble3(Values[rpos],Values[rpos+1],Values[rpos+2]);
  }
  else{
    const unsigned rpos=Nvalues*Position;
    const unsigned rpos2=Nvalues*PositionNext;
    const tdouble3 vini=TDouble3(Values[rpos],Values[rpos+1],Values[rpos+2]);
    const tdouble3 vnext=TDouble3(Values[rpos2],Values[rpos2+1],Values[rpos2+2]);
    ret=((vnext-vini)*TimeFactor+vini);
    if(SpecialValues){
      if(vini.x==DBL_MAX)ret.x=DBL_MAX;
      else if(vnext.x==DBL_MAX)ret.x=vini.x;
      if(vini.y==DBL_MAX)ret.y=DBL_MAX;
      else if(vnext.y==DBL_MAX)ret.y=vini.y;
      if(vini.z==DBL_MAX)ret.z=DBL_MAX;
      else if(vnext.z==DBL_MAX)ret.z=vini.z;
    }
  }
  return(ret);
}

//==============================================================================
/// Returns the interpolated values for the time indicated.
/// If no values always returns 0.
/// If only one value always returns that value.
/// If the indicated t is less than the minimum returns the first value.
/// If the indicated t is greater than the maximum returns the last value.
/// 
/// Devuelve valor el valor interpolado para el instante indicado.
/// Si no hay valores siempre devuelve 0.
/// Si solo hay un valor siempre devuelve ese valor.
/// Si el t indicado es menor que el minimo devuelve el primer valor.
/// Si el t indicado es mayor que el maximo devuelve el ultimo valor.
//==============================================================================
void JLinearValue::GetValue3d3d(double timestep,tdouble3 &v1,tdouble3 &v2){
  v1=v2=TDouble3(0);
  FindTime(timestep);
  //printf("--> t:%f  [%u - %u]  [%f - %f]\n",timestep,Position,PositionNext,TimePre,TimeNext);
  if(timestep<=TimePre){
    const unsigned rpos=Nvalues*Position;
    v1=TDouble3(Values[rpos  ],Values[rpos+1],Values[rpos+2]);
    v2=TDouble3(Values[rpos+3],Values[rpos+4],Values[rpos+5]);
  }
  else if(timestep>=TimeNext){
    const unsigned rpos=Nvalues*PositionNext;
    v1=TDouble3(Values[rpos  ],Values[rpos+1],Values[rpos+2]);
    v2=TDouble3(Values[rpos+3],Values[rpos+4],Values[rpos+5]);
  }
  else{
    const unsigned rpos =Nvalues*Position;
    const unsigned rpos2=Nvalues*PositionNext;
    const tdouble3 v1ini =TDouble3(Values[rpos   ],Values[rpos +1],Values[rpos +2]);
    const tdouble3 v2ini =TDouble3(Values[rpos+3 ],Values[rpos +4],Values[rpos +5]);
    const tdouble3 v1next=TDouble3(Values[rpos2  ],Values[rpos2+1],Values[rpos2+2]);
    const tdouble3 v2next=TDouble3(Values[rpos2+3],Values[rpos2+4],Values[rpos2+5]);
    v1=((v1next-v1ini)*TimeFactor+v1ini);
    v2=((v2next-v2ini)*TimeFactor+v2ini);
    if(SpecialValues){
      if(v1ini.x==DBL_MAX)v1.x=DBL_MAX;
      else if(v1next.x==DBL_MAX)v1.x=v1ini.x;
      if(v1ini.y==DBL_MAX)v1.y=DBL_MAX;
      else if(v1next.y==DBL_MAX)v1.y=v1ini.y;
      if(v1ini.z==DBL_MAX)v1.z=DBL_MAX;
      else if(v1next.z==DBL_MAX)v1.z=v1ini.z;
      if(v2ini.x==DBL_MAX)v2.x=DBL_MAX;
      else if(v2next.x==DBL_MAX)v2.x=v2ini.x;
      if(v2ini.y==DBL_MAX)v2.y=DBL_MAX;
      else if(v2next.y==DBL_MAX)v2.y=v2ini.y;
      if(v2ini.z==DBL_MAX)v2.z=DBL_MAX;
      else if(v2next.z==DBL_MAX)v2.z=v2ini.z;
    }
  }
}

//==============================================================================
/// Returns result of GetValue3d() converted to tfloat3 according special values.
//==============================================================================
tfloat3 JLinearValue::GetValue3f(double timestep){
  const tdouble3 v=GetValue3d(timestep);
  return(TFloat3((v.x==DBL_MAX? FLT_MAX: float(v.x)),(v.y==DBL_MAX? FLT_MAX: float(v.y)),(v.z==DBL_MAX? FLT_MAX: float(v.z))));
}

//==============================================================================
/// Devuelve el tiempo indicado.
/// Returns the indicated time.
//==============================================================================
double JLinearValue::GetTimeByIdx(unsigned idx)const{
  if(idx>=Count)Run_Exceptioon("Idx is not valid.");
  return(Times[idx]);
}

//==============================================================================
/// Devuelve el valor indicado.
/// Returns the indicated value.
//==============================================================================
double JLinearValue::GetValueByIdx(unsigned idx,unsigned cvalue)const{
  if(idx>=Count)Run_Exceptioon("Idx is not valid.");
  if(cvalue>=Nvalues)Run_Exceptioon("Number of value is not valid.");
  return(Values[Nvalues*idx+cvalue]);
}

//==============================================================================
/// Devuelve el valor indicado.
/// Returns the indicated value.
//==============================================================================
tdouble3 JLinearValue::GetValue3ByIdx(unsigned idx,unsigned cvalue3)const{
  if(idx>=Count)Run_Exceptioon("Idx is not valid.");
  if(cvalue3*3+2>=Nvalues)Run_Exceptioon("Number of value3 is not valid.");
  const unsigned cv=Nvalues*idx+cvalue3*3;
  return(TDouble3(Values[cv],Values[cv+1],Values[cv+2]));
}


//==============================================================================
/// Reads value and checks special values when SpecialValues is true.
/// Lee un valor comprobando si es especial cuando SpecialValues es true.
//==============================================================================
double JLinearValue::ReadNextDouble(JReadDatafile &rdat,bool in_line){
  const string value=rdat.ReadNextValue(in_line);
  double v=atof(value.c_str());
  if(SpecialValues && fun::StrLower(value)=="none")v=DBL_MAX;
  return(v);
}

//==============================================================================
/// Loads values for different times.
/// Carga valores para diferentes instantes.
//==============================================================================
void JLinearValue::LoadFile(std::string file){
  Reset();
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  SetSize(rows);
  for(unsigned r=0;r<rows;r++){
    const double atime=rdat.ReadNextDouble();  //-Time.
    const double v=ReadNextDouble(rdat,true);  //-Value_1.
    const unsigned idx=AddTimeValue(atime,v);
    for(unsigned cv=1;cv<Nvalues;cv++)SetValue(idx,cv,ReadNextDouble(rdat,true));
  }
  if(Count<2)Run_ExceptioonFile("Cannot be less than two values.",file);
  File=file;
}

//==============================================================================
/// Shows data in memory.
/// Muestra los datos en memoria.
//==============================================================================
void JLinearValue::VisuData(){
  printf("Time");
  for(unsigned cv=0;cv<Nvalues;cv++)printf("  v_%u",cv);
  printf("\n");
  for(unsigned c=0;c<Count;c++){
    printf("  %f",GetTimeByIdx(c));
    for(unsigned cv=0;cv<Nvalues;cv++)printf("  %f",GetValueByIdx(c,cv));
    printf("\n");
  }
}

#ifdef JLinearValue_UseJXml
//==============================================================================
/// Reads data from XML.
/// Lee datos del XML.
//==============================================================================
void JLinearValue::ReadXmlValues(const JXml *sxml,TiXmlElement* ele,std::string name
  ,std::string subname,std::string attributes)
{
  Reset();
  if(sxml->ExistsElement(ele,name)){
    File=sxml->ReadElementStr(ele,name,"file",true);
    if(File.empty()){
      TiXmlElement* xlis=sxml->GetFirstElement(ele,name);
      if(xlis){
        vector<string> attr;
        if(fun::VectorSplitStr(":",attributes,attr)!=Nvalues+1)Run_Exceptioon("Number of values does not match.");
        SetSize(sxml->CountElements(xlis,subname));
        TiXmlElement* elet=xlis->FirstChildElement(subname.c_str()); 
        while(elet){
          double t=sxml->GetAttributeDouble(elet,attr[0]); //-Reads time.
          if(SpecialValues){
            unsigned idx=UINT_MAX;
            for(unsigned ca=0;ca<Nvalues;ca++){
              double v=0;
              if(fun::StrLower(sxml->GetAttributeStr(elet,attr[ca+1],true))=="none")v=DBL_MAX;
              else v=sxml->GetAttributeDouble(elet,attr[ca+1],true,DBL_MAX);
              if(idx==UINT_MAX)idx=AddTimeValue(t,v);
              else SetValue(idx,ca,v);
            }
          }
          else{
            const unsigned idx=AddTimeValue(t,sxml->GetAttributeDouble(elet,attr[1],OptionalValues));
            for(unsigned ca=1;ca<Nvalues;ca++)SetValue(idx,ca,sxml->GetAttributeDouble(elet,attr[ca+1],OptionalValues));
          }
          elet=elet->NextSiblingElement(subname.c_str());
        }
      }
    }
  }
}

//==============================================================================
/// Writes data on XML.
/// Escribe datos en XML.
//==============================================================================
TiXmlElement* JLinearValue::WriteXmlValues(JXml *sxml,TiXmlElement* ele,std::string name
  ,std::string subname,std::string attributes)const
{
  TiXmlElement* rele=NULL;
  if(!GetFile().empty())rele=sxml->AddElementAttrib(ele,name,"file",GetFile());
  else{
    TiXmlElement* xlis=sxml->AddElement(ele,name);
    const unsigned nv=GetCount();
    //-Checks values.
    bool *vvoid=new bool[Nvalues];
    for(unsigned ca=0;ca<Nvalues;ca++)vvoid[ca]=true;
    for(unsigned c=0;c<nv;c++){
      for(unsigned ca=0;ca<Nvalues;ca++)if(vvoid[ca] && GetValueByIdx(c,ca)!=DBL_MAX)vvoid[ca]=false;
    }
    //-Write values.
    vector<string> attr;
    if(fun::VectorSplitStr(":",attributes,attr)!=Nvalues+1)Run_Exceptioon("Number of values does not match.");
    for(unsigned c=0;c<nv;c++){
      const double t=GetTimeByIdx(c);
      TiXmlElement* elet=sxml->AddElementAttrib(xlis,subname,attr[0],t);
      for(unsigned ca=0;ca<Nvalues;ca++){
        if(!vvoid[ca]){
          const double v=GetValueByIdx(c,ca);
          if(SpecialValues && v==DBL_MAX)sxml->AddAttribute(elet,attr[ca+1],string("none"));
          else sxml->AddAttribute(elet,attr[ca+1],v);
        }
      }
    }
    rele=xlis;
  }
  return(rele);
}

#endif

//==============================================================================
/// Cambia valores despues de la configuracion inicial.
/// Set values after the initial configuration.
//==============================================================================
void JLinearValue::RnSetValues(double t0,double v0,double t1,double v1){
  //-Prepares data.
  if(GetCount()!=2){
    SetSize(2);
    memset(Times,0,sizeof(double)*Size);
    memset(Values,0,sizeof(double)*Size*Nvalues);
    Count=2;
  }
  //-Set data.
  Times[0]=t0;          Values[0      ]=v0;
  Times[1]=max(t0,t1);  Values[Nvalues]=v1;
  Position=UINT_MAX;
}


