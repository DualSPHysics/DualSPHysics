//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JGauge_ker.h \brief Declares functions and CUDA kernels for classes JGauge.

#ifndef _JGauge_ker_
#define _JGauge_ker_

#include "Types.h"
#include <cuda_runtime_api.h>

/// Implements a set of functions and CUDA kernels for classes that manage gauges.
namespace cugauge{

inline int4 Int4(const tint4& v){ int4 p={v.x,v.y,v.z,v.w}; return(p); }
inline float3 Float3(const tfloat3& v){ float3 p={v.x,v.y,v.z}; return(p); }
inline float3 Float3(float x,float y,float z){ float3 p={x,y,z}; return(p); }
inline tfloat3 ToTFloat3(const float3& v){ return(TFloat3(v.x,v.y,v.z)); }
inline double3 Double3(const tdouble3& v){ double3 p={v.x,v.y,v.z}; return(p); }

dim3 GetGridSize(unsigned n,unsigned blocksize);

//-Kernel for JGaugeVelocity.
void Interaction_GaugeVel(tdouble3 ptpos
  ,float awen,int hdiv,tuint3 ncells,tuint3 cellmin,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,float3 *ptvel
  ,tdouble3 domposmin,float scell,float fourh2,float h,float massf);

//-Kernel for JGaugeSwl.
void Interaction_GaugeSwl(tdouble3 point0,tdouble3 pointdir,unsigned pointnp,float masslimit
  ,float awen,int hdiv,tuint3 ncells,tuint3 cellmin,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,tdouble3 domposmin,float scell,float fourh2,float h,float massf,float3 *ptres);

//-Kernel for JGaugeMaxZ.
void Interaction_GaugeMaxz(tdouble3 point0,float maxdist2
  ,int cxini,int cxfin,int yini,int yfin,int zini,int zfin
  ,int4 nc,unsigned cellfluid,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code
  ,float3 *ptres);

//-Kernel for JGaugeForce.
void Interaction_GaugeForce(unsigned n,unsigned idbegin,typecode codesel
  ,float fourh2,float h,float bwen,float massf,float cteb,float rhopzero,float gamma
  ,int hdiv,tuint3 ncells,tuint3 cellmin,const int2 *begincell,tdouble3 domposmin,float scell
  ,const double2 *posxy,const double *posz,const typecode *code,const unsigned *idp,const float4 *velrhop
  ,float3 *partace);


}


#endif


