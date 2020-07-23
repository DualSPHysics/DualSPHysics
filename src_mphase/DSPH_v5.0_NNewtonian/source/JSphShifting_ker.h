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
//:# - Implementacion inical. (01-12-2019)
//:# - Cambio de nombre de J.Shifting a J.SphShifting. (28-06-2020)
//:#############################################################################

/// \file JSphShifting_ker.h \brief Declares functions and CUDA kernels for Shifting correction on GPU.

#ifndef _JSphShifting_ker_
#define _JSphShifting_ker_

#include "TypesDef.h"
#include <cuda_runtime_api.h>

/// Implements a set of functions and CUDA kernels for Shifting correction on GPU.
namespace cushift{

//-Kernels for JSphShifting.
void InitGpuPosMax(bool tfirst,bool tdbl,unsigned n,unsigned pini
  ,const tdouble3& pmin1,const tdouble3& pmax1,const tdouble3& pmin2,const tdouble3& pmax2
  ,const double2* posxy,const double* posz,float4* shiftposfs,cudaStream_t stm);

void InitGpuPlanes(bool tfirst,bool tdbl,unsigned n,unsigned pini
  ,const tplane3d& plax1,const tplane3d& play1,const tplane3d& plaz1,const tdouble3& pladis1
  ,const tplane3d& plax2,const tplane3d& play2,const tplane3d& plaz2,const tdouble3& pladis2
  ,const double2* posxy,const double* posz,float4* shiftposfs,cudaStream_t stm);

void RunShifting(unsigned n,unsigned pini,double dt
  ,double coefumagn,float shifttfs,double coeftfs,float maxdist
  ,const float4 *velrhop,float4 *shiftposfs,cudaStream_t stm);

}

#endif


