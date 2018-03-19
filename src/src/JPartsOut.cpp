//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JPartsOut.cpp \brief Implements the class \ref JPartsOut.

#include "JPartsOut.h"
#include "Functions.h"
#include <algorithm>

using namespace std;

//##############################################################################
//# JPartsOut
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartsOut::JPartsOut(unsigned sizeini){
  ClassName="JPartsOut";
  SizeIni=sizeini;
  Idp=NULL; Pos=NULL; Vel=NULL; Rhop=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartsOut::~JPartsOut(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartsOut::Reset(){
  Clear();
  AllocMemory(0,true);
}

//==============================================================================
/// Resizes memory space for particles data.
//==============================================================================
void JPartsOut::AllocMemory(unsigned size,bool reset){
  if(reset){
    Count=0;
    delete[] Idp;  Idp=NULL;
    delete[] Pos;  Pos=NULL;
    delete[] Vel;  Vel=NULL;
    delete[] Rhop; Rhop=NULL;
  }
  Size=(!Size && size<SizeIni? SizeIni: size);
  Count=min(Count,Size);
  if(Size){
    try{
      Idp=fun::ResizeAlloc(Idp,Count,Size);
      Pos=fun::ResizeAlloc(Pos,Count,Size);
      Vel=fun::ResizeAlloc(Vel,Count,Size);
      Rhop=fun::ResizeAlloc(Rhop,Count,Size);
    }
    catch(const std::bad_alloc){
      RunException("AllocMemory","Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Returns the allocated memory in CPU.
//==============================================================================
llong JPartsOut::GetAllocMemory()const{  
  llong s=0;
  //Reservada en AllocMemory()
  //Allocated in AllocMemory()
  if(Idp)s+=sizeof(unsigned)*Size;
  if(Pos)s+=sizeof(tdouble3)*Size;
  if(Vel)s+=sizeof(tfloat3)*Size;
  if(Rhop)s+=sizeof(float)*Size;
  return(s);
}

//==============================================================================
/// Resizes arrays for particles.
//==============================================================================
void JPartsOut::AddParticles(unsigned np,const unsigned* idp,const tdouble3* pos,const tfloat3* vel,const float* rhop,unsigned outrhop,unsigned outmove){
  if(Count+np>Size)AllocMemory(Count+np+SizeIni,false);
  memcpy(Idp+Count,idp,sizeof(unsigned)*np);
  memcpy(Pos+Count,pos,sizeof(tdouble3)*np);
  memcpy(Vel+Count,vel,sizeof(tfloat3)*np);
  memcpy(Rhop+Count,rhop,sizeof(float)*np);
  Count+=np;
  OutPosCount+=np-(outrhop+outmove);
  OutRhopCount+=outrhop;
  OutMoveCount+=outmove;
}


