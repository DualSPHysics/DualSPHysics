//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2016, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphPartsInit.cpp \brief Implements the class \ref JSphPartsInit.

#include "JSphPartsInit.h"
#include "JSphMk.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsMath.h"

#include <cfloat>
#include <algorithm>
#include <cstring>

using namespace std;

//##############################################################################
//# JSphPartsInit
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphPartsInit::JSphPartsInit(bool simulate2d,double simulate2dposy,double dp
  ,const JSphMk* mkinfo,unsigned np,const tdouble3 *pos,const typecode *code)
  :Simulate2D(simulate2d),Simulate2DPosY(simulate2dposy),Dp(dp)
{
  ClassName="JSphPartsInit";
  Pos=NULL; Code=NULL;
  Reset();
  MkInfo=mkinfo;
  LoadParticleData(np,pos,code);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphPartsInit::~JSphPartsInit(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphPartsInit::Reset(){
  MkInfo=NULL;
  FreeParticleData();
}

//==============================================================================
/// Free memory for particle data.
//==============================================================================
void JSphPartsInit::FreeParticleData(){
  Np=0;
  delete[] Pos;  Pos=NULL;
  delete[] Code; Code=NULL;
}

//==============================================================================
/// Allocates dynamic memory and loads particle data.
//==============================================================================
void JSphPartsInit::LoadParticleData(unsigned np,const tdouble3 *pos,const typecode *code){
  FreeParticleData();
  Np=np;
  try{
    Pos=new tdouble3[Np];
    Code=new typecode[Np];
    memcpy(Pos,pos,sizeof(tdouble3)*Np);
    memcpy(Code,code,sizeof(typecode)*Np);
  }
  catch(const std::bad_alloc){
    RunException("LoadParticleData","Could not allocate the requested memory.");
  }
}



