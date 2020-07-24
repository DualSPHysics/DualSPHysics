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

/// \file JDsPartsOut.cpp \brief Implements the class \ref JDsPartsOut.

#include "JDsPartsOut.h"
#include "Functions.h"
#include <algorithm>

using namespace std;

//##############################################################################
//# JDsPartsOut
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsPartsOut::JDsPartsOut(unsigned sizeunit):SizeUnit(sizeunit)
{
  ClassName="JDsPartsOut";
  Idp=NULL; Pos=NULL; Vel=NULL; Rhop=NULL;
  Motive=NULL;
  Reset();
  AllocMemory(SizeUnit,true);
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsPartsOut::~JDsPartsOut(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsPartsOut::Reset(){
  Clear();
  AllocMemory(0,true);
  MemAllocs=0;
}

//==============================================================================
/// Resizes memory space for particles data.
//==============================================================================
void JDsPartsOut::AllocMemory(unsigned size,bool reset){
  if(reset){
    MemCpuParticles=0;
    Count=0;
    delete[] Idp;    Idp=NULL;
    delete[] Pos;    Pos=NULL;
    delete[] Vel;    Vel=NULL;
    delete[] Rhop;   Rhop=NULL;
    delete[] Motive; Motive=NULL;
  }
  Size=unsigned((size+SizeUnit-1)/SizeUnit)*SizeUnit;
  Count=min(Count,Size);
  if(Size){
    MemAllocs++;
    try{
      Idp   =fun::ResizeAlloc(Idp   ,Count,Size);  MemCpuParticles+=sizeof(unsigned)*Size;
      Pos   =fun::ResizeAlloc(Pos   ,Count,Size);  MemCpuParticles+=sizeof(tdouble3)*Size;
      Vel   =fun::ResizeAlloc(Vel   ,Count,Size);  MemCpuParticles+=sizeof(tfloat3) *Size;
      Rhop  =fun::ResizeAlloc(Rhop  ,Count,Size);  MemCpuParticles+=sizeof(float)   *Size;
      Motive=fun::ResizeAlloc(Motive,Count,Size);  MemCpuParticles+=sizeof(byte)    *Size;
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Adds motive information and updates numbers.
//==============================================================================
void JDsPartsOut::AddData(unsigned np,const typecode* code){
  //-Checks reason for exclusion.
  unsigned outpos=0,outrhop=0,outmove=0;
  for(unsigned c=0;c<np;c++){
    switch(CODE_GetSpecialValue(code[c])){
      case CODE_OUTPOS:   Motive[Count+c]=1; outpos++;   break;
      case CODE_OUTRHOP:  Motive[Count+c]=2; outrhop++;  break; 
      case CODE_OUTMOVE:  Motive[Count+c]=3; outmove++;  break; 
    }
  }
  //-Updates numbers.
  Count+=np;
  OutPosCount+=outpos;
  OutRhopCount+=outrhop;
  OutMoveCount+=outmove;
  //printf("\n=====> count:%d outpos:%d\n",Count,OutPosCount);
}

//==============================================================================
/// Adds out particles data.
//==============================================================================
void JDsPartsOut::AddParticles(unsigned np,const unsigned* idp,const tdouble3* pos
  ,const tfloat3* vel,const float* rhop,const typecode* code)
{
  if(Count+np>Size)AllocMemory(Count+np+SizeUnit,false);
  memcpy(Idp +Count,idp ,sizeof(unsigned)*np);
  memcpy(Pos +Count,pos ,sizeof(tdouble3)*np);
  memcpy(Vel +Count,vel ,sizeof(tfloat3 )*np);
  memcpy(Rhop+Count,rhop,sizeof(float   )*np);
  //-Adds motive information and updates numbers.
  AddData(np,code);
}



