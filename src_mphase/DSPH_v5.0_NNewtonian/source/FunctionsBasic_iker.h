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
//:# - Implementacion inicial. (01-12-2019)
//:#############################################################################

/// \file FunctionsBasic_iker.h \brief Implements basic functions for CUDA files.

#include "TypesDef.h"
#include <cuda_runtime_api.h>

#ifndef SPHBSIZE
  #define SPHBSIZE 256  //-CUDA blocksize by default.
#endif

//==============================================================================
// Basic conversion functions.
//==============================================================================
inline int3 Int3(const tint3& v){ int3 p={v.x,v.y,v.z}; return(p); }
inline int4 Int4(const tint4& v){ int4 p={v.x,v.y,v.z,v.w}; return(p); }
inline float3 Float3(const tfloat3& v){ float3 p={v.x,v.y,v.z}; return(p); }
inline float3 Float3(float x,float y,float z){ float3 p={x,y,z}; return(p); }
inline float4 Float4(const tfloat4& v){ float4 p={v.x,v.y,v.z,v.w}; return(p); }
inline float4 Float4(const tplane3f& v){ float4 p={v.a,v.b,v.c,v.d}; return(p); }
inline double3 Double3(const tdouble3& v){ double3 p={v.x,v.y,v.z}; return(p); }
inline double4 Double4(const tdouble4& v){ double4 p={v.x,v.y,v.z,v.w}; return(p); }
inline double4 Double4(const tplane3d& v){ double4 p={v.a,v.b,v.c,v.d}; return(p); }


//==============================================================================
// Returns gridsize according parameters using X dimension.
//==============================================================================
inline dim3 GetSimpleGridSize(unsigned n,unsigned blocksize){
  const unsigned nb=unsigned(n+blocksize-1)/blocksize;//-Numero total de bloques a lanzar.
  return(dim3(nb,1,1));
}

////==============================================================================
///// Returns the dimensions of gridsize according to parameters.
///// Devuelve tamanho de gridsize segun parametros.
////==============================================================================
//inline dim3 GetGridSize_Old(unsigned n,unsigned blocksize){
//  dim3 sgrid;//=dim3(1,2,3);
//  unsigned nb=unsigned(n+blocksize-1)/blocksize; //-Total number of blocks to execute.
//  sgrid.x=(nb<=65535? nb: unsigned(sqrt(float(nb))));
//  sgrid.y=(nb<=65535? 1: unsigned((nb+sgrid.x-1)/sgrid.x));
//  sgrid.z=1;
//  return(sgrid);
//}

