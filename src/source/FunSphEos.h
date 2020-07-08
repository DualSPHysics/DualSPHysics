//HEAD_DSCODES
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

/// \file FunSphEos.h \brief implements inline functions for Equation of State for SPH.

#ifndef _FunSphEos_
#define _FunSphEos_

#include "TypesDef.h"
#include "DualSphDef.h"
#include <cmath>

/// Implements inline functions for Equation of State for SPH.
namespace fsph{

//##############################################################################
//# Equation of State - Monaghan, 1994
//##############################################################################
//==============================================================================
/// Returns pressure starting from density using equation of state 
/// based on [Monaghan, 1994].
//==============================================================================
inline float ComputePressMonaghan(float rhop,float rhop0,float b,float gamma){ 
  return(b*(pow(rhop/rhop0,gamma)-1.0f));
}
//==============================================================================
/// Returns pressure starting from density using equation of state 
/// based on [Monaghan, 1994].
//==============================================================================
inline float ComputePressMonaghan(float rhop,const StCteSph &csp){
  return(ComputePressMonaghan(rhop,csp.rhopzero,csp.cteb,csp.gamma));
}


//##############################################################################
//# Default Equation of State.
//##############################################################################
//------------------------------------------------------------------------------
/// Returns pressure starting from density using default equation of state.
//------------------------------------------------------------------------------
inline float ComputePress(float rhop,float rhop0,float b,float gamma,float cs0){ 
  return(ComputePressMonaghan(rhop,rhop0,b,gamma));
}

//------------------------------------------------------------------------------
/// Returns pressure starting from density using default equation of state.
//------------------------------------------------------------------------------
inline float ComputePress(float rhop,const StCteSph &csp){ 
  return(ComputePressMonaghan(rhop,csp));
}


}

#endif


