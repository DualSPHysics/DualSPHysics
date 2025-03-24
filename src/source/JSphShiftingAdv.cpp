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

/// \file JSphShifting.cpp \brief Implements the class \ref JSphShifting.

#include "JSphShiftingAdv.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "JXml.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JDataArrays.h"

#include <cfloat>
#include <algorithm>

using std::string;


//##############################################################################
//# JSphShiftingAdv
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphShiftingAdv::JSphShiftingAdv(bool simulate2d,double dp,float kernelh)
  :Log(AppInfo.LogPtr()),Simulate2D(simulate2d),Dp(dp),KernelH(kernelh)
{
  ClassName="JSphShifting";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphShiftingAdv::~JSphShiftingAdv(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphShiftingAdv::Reset(){
  ConfigBasic(-0.01f,false,false);
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSphShiftingAdv::VisuConfig(std::string txhead,std::string txfoot){
  if(!txhead.empty())Log->Print(txhead);
  Log->Print(fun::VarStr("Shifting","Advanced"));
  Log->Print(fun::VarStr("  ShiftAdvCoef",ShiftCoef));
  Log->Print(fun::VarStr("  ShiftAdvALE",AleActive));
  Log->Print(fun::VarStr("  ShiftAdvNCPress",NcPress));
}

//==============================================================================
/// Returns configuration information in one string line.
//==============================================================================
std::string JSphShiftingAdv::GetConfigInfo()const{
  string ret=string("Shifting(Advanced")+fun::PrintStr(",%g,%s,%s)"
    ,ShiftCoef,(AleActive? "1": "0"),(NcPress? "1": "0"));
  return(ret);
}

//==============================================================================
/// Explicit configuration from execution parameters.
//==============================================================================
void JSphShiftingAdv::ConfigBasic(float shiftcoef,bool aleactive,bool ncpress){
  ShiftCoef=shiftcoef;
  AleActive=aleactive;
  NcPress=ncpress;
}

