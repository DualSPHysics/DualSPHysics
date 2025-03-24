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

/// \file JSphGpu_mdbc_iker.h \brief Declares functions and CUDA kernels for mDBC.

#ifndef _JSphGpu_mdbc_iker_
#define _JSphGpu_mdbc_iker_

#include "DualSphDef.h"
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>

/// Implements a set of functions and CUDA kernels for mDBC.
namespace cusph{

//-Kernels for the boundary correction (mDBC).
void Interaction_MdbcCorrection(TpKernel tkernel,bool simulate2d,unsigned n
  ,unsigned nbound,const StDivDataGpu& dvd,const tdouble3& mapposmin
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const typecode* code,const unsigned* idp,const float3* boundnor
  ,float4* velrho,cudaStream_t stm=NULL);

//<vs_m2dbc_ini>
//-Kernels for the boundary correction (mDBC2).
void Interaction_Mdbc2Correction(TpKernel tkernel,bool simulate2d
  ,TpSlipMode slipmode,unsigned n,unsigned nbound,const tfloat3 gravity
  ,const StDivDataGpu& dvd,const tdouble3& mapposmin,const double2* posxy
  ,const double* posz,const float4* poscell,const typecode* code
  ,const unsigned* idp,const float3* boundnor,const float3* motionvel
  ,const float3* motionace,float4* velrho,byte* boundmode,float3* tangenvel
  ,cudaStream_t stm=NULL);

void CopyMotionVelAce(unsigned nmoving,double dt,const unsigned* ridpmot
  ,const float4* velrho,float3* motionvel,float3* motionace);
//<vs_m2dbc_end>

}


#endif

