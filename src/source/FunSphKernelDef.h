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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Creacion para implementacion de kernels SPH. (03-07-2020)
//:#############################################################################

/// \file FunSphKernelDef.h \brief Defines structures for SPH kenerls.

#ifndef _FunSphKernelDef_
#define _FunSphKernelDef_

/// Implements a set of basic/general functions related to SPH.
namespace fsph{

//##############################################################################
//# Cubic Spline kernel
//##############################################################################
///Structure with constants for the Cubic Spline kernel.
typedef struct {
  float a1,a2,aa,a24,c1,d1,c2;
  float od_wdeltap;        ///<Parameter for tensile instability correction.  
}StKCubicCte;

//##############################################################################
//# Wendland kernel
//##############################################################################
///Structure with constants for the Wendland kernel.
typedef struct {
  float awen;  ///<Constant to compute wab.
  float bwen;  ///<Constant to compute fac (kernel derivative).
}StKWendlandCte;

}

#endif


