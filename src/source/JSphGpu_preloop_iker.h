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

/// \file JSphGpu_preloop_iker.h \brief Declares functions and CUDA kernels for preloop.

#ifndef _JSphGpu_preloop_iker_
#define _JSphGpu_preloop_iker_

#include "DualSphDef.h"
#include "JCellDivDataGpu.h"
#include "JSphGpu_ker.h"
#include <cuda_runtime_api.h>

/// Implements a set of functions and CUDA kernels for preloop.
namespace cusph{

void ComputeFSNormals(TpKernel tkernel,bool simulate2d,unsigned bsfluid
  ,unsigned fluidini,unsigned fluidnum,StDivDataGpu& dvd,const unsigned* dcell
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const float4* velrho,const typecode* code,const float* ftomassp
  ,float4* shiftposfs,unsigned* fstype,float3* fsnormal,unsigned* listp
  ,cudaStream_t stm);
    
void ComputeUmbrellaRegion(TpKernel tkernel,bool simulate2d
  ,unsigned bsfluid,unsigned fluidini,unsigned fluidnum,StDivDataGpu& dvd
  ,const unsigned* dcell,const float4* poscell,const typecode* code
  ,const float3* fsnormal,unsigned* listp,unsigned* fstype,cudaStream_t stm);

void PreLoopInteraction(TpKernel tkernel,bool simulate2d,bool shiftadv
  ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd
  ,const unsigned* dcell,const float4* poscell,const float4* velrho
  ,const typecode* code,const float* ftomassp,float4* shiftvel,unsigned* fstype
  ,float3* fsnormal,float* fsmindist,cudaStream_t stm);

void ComputeShiftingVel(unsigned bsfluid,unsigned fluidnum,unsigned fluidini
  ,bool sim2d,float shiftcoef,bool ale,float dt,const unsigned* fstype
  ,const float3* fsnormal,const float* fsmindist,float4* shiftvel
  ,cudaStream_t stm);    

void PeriodicSaveParent(unsigned n,unsigned pini,const unsigned* listp
  ,unsigned* periparent);

void PeriPreLoopCorr(unsigned n,unsigned pinit,const unsigned* periparent
  ,unsigned* fstype,float4* shiftvel);
}

#endif

