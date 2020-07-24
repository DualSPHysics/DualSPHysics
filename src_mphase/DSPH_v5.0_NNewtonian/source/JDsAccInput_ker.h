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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Implementacion en fichero independiente. (22-12-2019)
//:# - Cambio de nombre de J.SphAccInput a J.DsAccInput. (28-06-2020)
//:#############################################################################

/// \file JDsAccInput_ker.h \brief Declares functions and CUDA kernels for external forces (JDsAccInput) on GPU.

#ifndef _JDsAccInput_ker_
#define _JDsAccInput_ker_

#include "DualSphDef.h"
#include <cuda_runtime_api.h>

/// Implements a set of functions and CUDA kernels for external forces (JDsAccInput) on GPU.
namespace cuaccin{

//-Kernels for external forces (JDsAccInput).
void AddAccInput(unsigned n,unsigned pini,typecode codesel1,typecode codesel2
  ,tdouble3 acclin,tdouble3 accang,tdouble3 centre,tdouble3 velang,tdouble3 vellin,bool setgravity
  ,tfloat3 gravity,const typecode *code,const double2 *posxy,const double *posz
  ,const float4 *velrhop,float3 *ace,cudaStream_t stm);

}

#endif


