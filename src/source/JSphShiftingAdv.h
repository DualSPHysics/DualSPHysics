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

/// \file JSphShiftingAdv.h \brief Declares the class \ref JSphShiftingAdv.

#ifndef _JSphShiftingAdv_
#define _JSphShiftingAdv_

#include "JObject.h"
#include "DualSphDef.h"
#include <string>
#include <vector>

class JLog2;
class JXml;
class TiXmlElement;

//##############################################################################
//# JSphShiftingAdv
//##############################################################################
/// \brief Manages the Shifting advanced model.

class JSphShiftingAdv : protected JObject
{
private:
  JLog2* Log;
  const bool Simulate2D;
  const double Dp;       ///<Initial distance between particles [m].
  float KernelH;         ///<The smoothing length of SPH kernel [m].

  float ShiftCoef;       ///<Coefficient for shifting computation.

  bool  AleActive;       ///<Arbitrarian Eulerian-Lagrangian model.
  bool  NcPress;         ///<Non conservative pressure formulation.

public:
  JSphShiftingAdv(bool simulate2d,double dp,float kernelh);
  ~JSphShiftingAdv();
  void Reset();

  void ConfigBasic(float shiftcoef=-0.01f,bool aleactive=false,bool ncpress=false);
  void VisuConfig(std::string txhead="",std::string txfoot="");
  std::string GetConfigInfo()const;

  float       GetShiftCoef()const{ return(ShiftCoef); }
  bool        GetAleActive()const{ return(AleActive); }
  bool        GetNcPress  ()const{ return(NcPress); }

};

#endif
