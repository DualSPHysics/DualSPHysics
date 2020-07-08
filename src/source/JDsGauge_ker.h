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

/// \file JDsGauge_ker.h \brief Declares functions and CUDA kernels for classes JGauge.

#ifndef _JDsGauge_ker_
#define _JDsGauge_ker_

#include "DualSphDef.h"
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>

/// Implements a set of functions and CUDA kernels for classes that manage gauges.
namespace cugauge{

//-Kernel for JGaugeVelocity.
void Interaction_GaugeVel(const StCteSph &CSP,const StDivDataGpu &dvd
  ,tdouble3 ptpos,const double2 *posxy,const double *posz
  ,const typecode *code,const float4 *velrhop,float3 *ptvel);

//-Kernel for JGaugeSwl.
void Interaction_GaugeSwl(const StCteSph &CSP,const StDivDataGpu &dvd
  ,tdouble3 point0,tdouble3 pointdir,unsigned pointnp,float masslimit
  ,const double2 *posxy,const double *posz,const typecode *code
  ,const float4 *velrhop,float3 *ptres);

//-Kernel for JGaugeMaxZ.
void Interaction_GaugeMaxz(tdouble3 point0,float maxdist2,const StDivDataGpu &dvd
  ,int cxini,int cxfin,int yini,int yfin,int zini,int zfin
  ,const double2 *posxy,const double *posz,const typecode *code
  ,float3 *ptres);

//-Kernel for JGaugeForce.
void Interaction_GaugeForce(const StCteSph &CSP,const StDivDataGpu &dvd
  ,unsigned n,unsigned idbegin,typecode codesel,const double2 *posxy,const double *posz
  ,const typecode *code,const unsigned *idp,const float4 *velrhop
  ,float3 *partace);


}


#endif


