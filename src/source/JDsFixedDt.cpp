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

/// \file JDsFixedDt.cpp \brief Implements the class \ref JDsFixedDt.

#include "JDsFixedDt.h"
#include "Functions.h"
#include "JReadDatafile.h"
#include <cstring>
#include <cfloat>
#include <algorithm>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JDsFixedDt::JDsFixedDt(double fixedvalue){
  ClassName="JDsFixedDt";
  Times=NULL;
  Values=NULL;
  Reset();
  FixedValue=fixedvalue;
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsFixedDt::~JDsFixedDt(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsFixedDt::Reset(){
  FixedValue=0;
  delete[] Times;  Times=NULL;
  delete[] Values; Values=NULL;
  File="";
  Size=Count=Position=0;
  GetDtError(true);
  LastTimestepInput=LastDtInput=LastDtOutput=-1;
}

//==============================================================================
/// Resizes memory space for values.
//==============================================================================
void JDsFixedDt::Resize(unsigned size){
  Times=fun::ResizeAlloc(Times,Count,size);
  Values=fun::ResizeAlloc(Values,Count,size);
  Size=size;
}

//==============================================================================
/// Returns the allocated memory.
//==============================================================================
unsigned JDsFixedDt::GetAllocMemory()const{
  unsigned s=0;
  if(Times)s+=sizeof(double)*Size;
  if(Values)s+=sizeof(double)*Size;
  return(s);
}

//==============================================================================
/// Loads values of dt (MILISECONDS) for different instants (in SECONDS).
//==============================================================================
void JDsFixedDt::LoadFile(std::string file){
  Reset();
  JReadDatafile rdat;
  rdat.LoadFile(file,FILESIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  Resize(rows);
  for(unsigned r=0;r<rows;r++){
    Times[r]=rdat.ReadNextDouble(false);
    Values[r]=rdat.ReadNextDouble(true);
    //printf("FileData[%u]>  t:%f  ang:%f\n",r,Times[r],Values[r]);
  }
  Count=rows;
  if(Count<2)Run_ExceptioonFile("Cannot be less than two values.",file);
  File=file;
}

//==============================================================================
/// Returns the value of dt (in SECONDS) for a given instant.
//==============================================================================
double JDsFixedDt::GetDt(double timestep,double dtvar){
  if(LastTimestepInput>=0 && timestep==LastTimestepInput && dtvar==LastDtInput)return(LastDtOutput);
  LastTimestepInput=timestep;
  LastDtInput=dtvar;
  double ret=FixedValue;
  if(!ret){
    //-Busca intervalo del instante indicado.
    //-Searches indicated interval of time.
    double tini=Times[Position];
    double tnext=(Position+1<Count? Times[Position+1]: tini);
    for(;tnext<timestep&&Position+2<Count;Position++){
      tini=tnext;
      tnext=Times[Position+2];
    }
    //-Calcula dt en el instante indicado.
    //-Computes dt for the indicated instant.
    if(timestep<=tini)ret=Values[Position]/1000;
    else if(timestep>=tnext)ret=Values[Position+1]/1000;
    else{
      const double tfactor=(timestep-tini)/(tnext-tini);
      double vini=Values[Position];
      double vnext=Values[Position+1];
      ret=(tfactor*(vnext-vini)+vini)/1000;
    }
  }
  const double dterror=ret-dtvar;
  if(DtError<dterror)DtError=dterror;
  LastDtOutput=ret;
  return(ret);
}

//==============================================================================
/// Returns the maximum error regarding the dtvariable. max(FixedDt-DtVariable).
//==============================================================================
double JDsFixedDt::GetDtError(bool reset){
  double ret=DtError;
  if(reset)DtError=-DBL_MAX;
  return(ret);
}


