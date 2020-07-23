//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsPartsInit.cpp \brief Implements the class \ref JDsPartsInit.

#include "JDsPartsInit.h"
#include "JSphMk.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsMath.h"

#include <cfloat>
#include <algorithm>
#include <cstring>

using namespace std;

//##############################################################################
//# JDsPartsInit
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsPartsInit::JDsPartsInit(bool simulate2d,double simulate2dposy,double dp
  ,const JSphMk* mkinfo,unsigned np,const tdouble3 *pos,const typecode *code)
  :Simulate2D(simulate2d),Simulate2DPosY(simulate2dposy),Dp(dp)
{
  ClassName="JDsPartsInit";
  Pos=NULL; Code=NULL;
  Reset();
  MkInfo=mkinfo;
  LoadParticleData(np,pos,code);
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsPartsInit::~JDsPartsInit(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsPartsInit::Reset(){
  MkInfo=NULL;
  FreeParticleData();
}

//==============================================================================
/// Free memory for particle data.
//==============================================================================
void JDsPartsInit::FreeParticleData(){
  Np=0;
  delete[] Pos;  Pos=NULL;
  delete[] Code; Code=NULL;
}

//==============================================================================
/// Allocates dynamic memory and loads particle data.
//==============================================================================
void JDsPartsInit::LoadParticleData(unsigned np,const tdouble3 *pos,const typecode *code){
  FreeParticleData();
  Np=np;
  try{
    Pos=new tdouble3[Np];
    Code=new typecode[Np];
    memcpy(Pos,pos,sizeof(tdouble3)*Np);
    memcpy(Code,code,sizeof(typecode)*Np);
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}



