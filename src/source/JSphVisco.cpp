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

/// \file JSphVisco.cpp \brief Implements the class \ref JSphVisco.

#include "JSphVisco.h"
#include "Functions.h"
#include "JReadDatafile.h"
#include <cstring>
#include <float.h>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JSphVisco::JSphVisco(){
  ClassName="JSphVisco";
  Times=NULL;
  Values=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphVisco::~JSphVisco(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphVisco::Reset(){
  delete[] Times;  Times=NULL;
  delete[] Values; Values=NULL;
  File="";
  Size=Count=Position=0;
  LastTimestepInput=LastViscoOutput=-1;
}

//==============================================================================
/// Redimensiona espacio para valores.
/// Resizes allocated space for values.
//==============================================================================
void JSphVisco::Resize(unsigned size){
  Times=fun::ResizeAlloc(Times,Count,size);
  Values=fun::ResizeAlloc(Values,Count,size);
  Size=size;
}

//==============================================================================
/// Devuelve la memoria reservada.
/// Returns the allocated memory.
//==============================================================================
unsigned JSphVisco::GetAllocMemory()const{
  unsigned s=0;
  if(Times)s+=sizeof(float)*Size;
  if(Values)s+=sizeof(float)*Size;
  return(s);
}

//==============================================================================
/// Carga valores de viscosidad para diferentes instantes (en segundos).
/// Loads viscosity values for different instants (in secods).
//==============================================================================
void JSphVisco::LoadFile(std::string file){
  Reset();
  JReadDatafile rdat;
  rdat.LoadFile(file,FILESIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  Resize(rows);
  for(unsigned r=0;r<rows;r++){
    Times[r]=rdat.ReadNextFloat(false);
    Values[r]=rdat.ReadNextFloat(true);
    //printf("FileData[%u]>  t:%f  ang:%f\n",r,Times[r],Values[r]);
  }
  Count=rows;
  if(Count<2)Run_ExceptioonFile("Cannot be less than two values.",file);
  File=file;
}

//==============================================================================
/// Devuelve el valor de viscosidad para el instante indicado.
/// Returns the viscosity value for the indicated instant.
//==============================================================================
float JSphVisco::GetVisco(float timestep){
  if(LastTimestepInput>=0 && timestep==LastTimestepInput)return(LastViscoOutput);
  LastTimestepInput=timestep;
  float ret=0;
  //-Busca intervalo del instante indicado.
  //-Searches indicated interval of time.
  float tini=Times[Position];
  float tnext=(Position+1<Count? Times[Position+1]: tini);
  for(;tnext<timestep&&Position+2<Count;Position++){
    tini=tnext;
    tnext=Times[Position+2];
  }
  //-Calcula dt en el instante indicado.
  //-Computes dt for the indicated instant.
  if(timestep<=tini)ret=Values[Position];
  else if(timestep>=tnext)ret=Values[Position+1];
  else{
    const double tfactor=double(timestep-tini)/double(tnext-tini);
    float vini=Values[Position];
    float vnext=Values[Position+1];
    ret=float(tfactor*(vnext-vini)+vini);
  }
  LastViscoOutput=ret;
  return(ret);
}


