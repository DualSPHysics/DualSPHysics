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

/// \file FunSphEos_iker.h \brief implements CUDA device functions for Equation of State for SPH.

#include "TypesDef.h"
#include "DualSphDef.h"
#include <cuda_runtime_api.h>

/// Implements CUDA device functions for Equation of State for SPH.
namespace cufsph{

//##############################################################################
//# Equation of State - Monaghan, 1994
//##############################################################################
//#define CTE_AVAILABLE
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Returns pressure starting from density using equation of state 
/// based on [Monaghan, 1994].
//------------------------------------------------------------------------------
__device__ float ComputePressMonaghanCte(float rhop){ 
  return(CTE.cteb*(pow(rhop*CTE.ovrhopzero,CTE.gamma)-1.0f));
}
#endif
//------------------------------------------------------------------------------
/// Returns pressure starting from density using equation of state 
/// based on [Monaghan, 1994].
//------------------------------------------------------------------------------
__device__ float ComputePressMonaghan(float rhop,float rhop0,float b,float gamma){ 
  return(b*(pow(rhop/rhop0,gamma)-1.0f));
}


//##############################################################################
//# Default Equation of State.
//##############################################################################
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Returns pressure starting from density using default equation of state.
//------------------------------------------------------------------------------
__device__ float ComputePressCte(float rhop){ 
  return(ComputePressMonaghanCte(rhop));
}
#endif
//------------------------------------------------------------------------------
/// Returns pressure starting from density using default equation of state.
//------------------------------------------------------------------------------
__device__ float ComputePress(float rhop,float rhop0,float b,float gamma,float cs0){ 
  return(ComputePressMonaghan(rhop,rhop0,b,gamma));
}


}