//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
JLinearValue::JLinearValue(unsigned nvalues):Nvalues(max(1u,nvalues)){
  ClassName="JLinearValue";
  Times=NULL;
  Values=NULL;
  Reset();
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
  if(size>=SIZEMAX)RunException("SetSize","It has reached the maximum size allowed.");
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
/// Adds values at the end of the list.
/// Añade valores al final de la lista.
//==============================================================================
unsigned JLinearValue::AddTimeValue(double time,double value){
  if(Count==Size)SetSize(Size+SIZEINITIAL);
  const unsigned idx=Count;
  Times[idx]=time;
  Values[Nvalues*idx]=value;
  for(unsigned cv=1;cv<Nvalues;cv++)Values[Nvalues*idx+cv]=0;
  Count++;
  return(idx);
}

//==============================================================================
/// Modify value at the indicated position.
/// Modifica valor en la posicion indicada.
//==============================================================================
void JLinearValue::SetValue(unsigned idx,unsigned cvalue,double value){
  const char met[]="SetValue";
  if(idx>=Count)RunException(met,"idx is not valid.");
  if(cvalue>=Nvalues)RunException(met,"cvalue is not valid.");
  Values[Nvalues*idx+cvalue]=value;
}

//==============================================================================
/// Finds time before and after indicated time.
/// Busca time anterior y posterior del instante indicado.
//==============================================================================
void JLinearValue::FindTime(double timestep){
  if(!Count)RunException("FindTime","There are not times.");
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
/// Returns the interpolated value value for the time indicated.
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
  }
  return(ret);
}

//==============================================================================
/// Devuelve el tiempo indicado.
/// Returns the indicated time.
//==============================================================================
double JLinearValue::GetTimeByIdx(unsigned idx)const{
  if(idx>=Count)RunException("GetTimeByIdx","idx is not valid.");
  return(Times[idx]);
}

//==============================================================================
/// Devuelve el valor indicado.
/// Returns the indicated value.
//==============================================================================
double JLinearValue::GetValueByIdx(unsigned idx,unsigned cvalue)const{
  if(idx>=Count)RunException("GetValueByIdx","idx is not valid.");
  if(cvalue>=Nvalues)RunException("GetValueByIdx","cvalue is not valid.");
  return(Values[Nvalues*idx+cvalue]);
}

//==============================================================================
/// Loads values for different times.
/// Carga valores para diferentes instantes.
//==============================================================================
void JLinearValue::LoadFile(std::string file){
  const char met[]="LoadFile";
  Reset();
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  SetSize(rows);
  for(unsigned r=0;r<rows;r++){
    const double atime=rdat.ReadNextDouble();  //-Time.
    const double v=rdat.ReadNextDouble(true);  //-Value_1.
    const unsigned idx=AddTimeValue(atime,v);
    for(unsigned cv=1;cv<Nvalues;cv++)SetValue(idx,cv,rdat.ReadNextDouble(true));
  }
  if(Count<2)RunException(met,"Cannot be less than two values.",file);
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


