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
//:# - Funciones: PointInMinMax(), PlanePoint() y PlanesDomainCheck(). (01-12-2019)
//:#############################################################################

/// \file FunctionsGeo3d_iker.h \brief Implements geometry functions for 3D on CUDA.

#include "TypesDef.h"
#include <cuda_runtime_api.h>

namespace cugeo{

//------------------------------------------------------------------------------
/// Devuelve true cuando pmin <= pt <= pmax.
/// Returns true when pmin <= pt <= pmax.
//------------------------------------------------------------------------------
__device__ bool PointInMinMax(const double3 &pt,const double3 &pmin,const double3 &pmax){
  return(pmin.x<=pt.x && pmin.y<=pt.y && pmin.z<=pt.z && pt.x<=pmax.x && pt.y<=pmax.x && pt.z<=pmax.z);
}

//------------------------------------------------------------------------------
/// Resuelve punto en el plano.
/// Solves point in the plane.
//------------------------------------------------------------------------------
__device__ double PlanePoint(const double4 &pla,const double3 &pt){ 
  return(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w);
}

//------------------------------------------------------------------------------
/// Resuelve punto en el plano.
/// Solves point in the plane.
//------------------------------------------------------------------------------
__device__ float PlanePoint(const float4 &pla,const float3 &pt){ 
  return(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w);
}

//------------------------------------------------------------------------------
/// Comprueba si el punto esta dentro del dominio definido.
/// Checks the point is inside the defined domain.
//------------------------------------------------------------------------------
__device__ bool PlanesDomainCheck(const double3 &pt,const double4 &plax
  ,const double4 &play,const double4 &plaz,const double3 &pladist)
{
  const double dx=PlanePoint(plax,pt);
  const double dy=PlanePoint(play,pt);
  const double dz=PlanePoint(plaz,pt);
  return(dx>=0 && dx<=pladist.x && dy>=0 && dy<=pladist.y && dz>=0 && dz<=pladist.z);
}

//------------------------------------------------------------------------------
/// Resuelve punto en el plano.
/// Solves point in the plane.
//------------------------------------------------------------------------------
__device__ float PlanePoint(const float4 &pla,const float &ptx,const float &pty,const float &ptz){ 
  return(pla.x*ptx+pla.y*pty+pla.z*ptz+pla.w);
}

//------------------------------------------------------------------------------
/// Devuelve la distancia entre un punto y un plano con signo.
/// Returns the distance between a point and a plane with sign.
//------------------------------------------------------------------------------
__device__ float PlaneDistSign(const float4 &pla,const float &ptx,const float &pty,const float &ptz){ 
  return(PlanePoint(pla,ptx,pty,ptz)/sqrt(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z));
}

//==============================================================================
/// Devuelve la distancia al (0,0,0).
/// Returns the distance from (0,0,0).
//==============================================================================
__device__ float PointDist(const float3 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}

//==============================================================================
/// Devuelve la distancia al (0,0,0).
/// Returns the distance from (0,0,0).
//==============================================================================
__device__ float PointDist(const float4 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}

//------------------------------------------------------------------------------
/// Devuelve la distancia^2 entre dos puntos.
/// Returns the distance^2 between two points.
//------------------------------------------------------------------------------
__device__ float PointsDist2(const float3 &p1,const float3 &p2){
  const float dx=p1.x-p2.x;
  const float dy=p1.y-p2.y;
  const float dz=p1.z-p2.z;
  return(dx*dx + dy*dy + dz*dz);
}

//------------------------------------------------------------------------------
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//------------------------------------------------------------------------------
__device__ float3 VecUnitary(const float3 &v){
  const float m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  return(m? make_float3(v.x/m,v.y/m,v.z/m): v);
}

//------------------------------------------------------------------------------
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//------------------------------------------------------------------------------
__device__ float3 VecModule(const float3 &v,float module){
  const float m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  const float m2=(m? module/m: 0);
  return(m? make_float3(v.x*m2,v.y*m2,v.z*m2): v);
}

//==============================================================================
/// Devuelve el producto escalar de 2 vectores.
/// Returns the scalar product of two vectors.
//==============================================================================
__device__ float ProductScalar(const float3 &v1,const float3 &v2){
  return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

//------------------------------------------------------------------------------
/// Devuelve vector de rebote a una normal.
/// Returns bounce vector to a normal.
//------------------------------------------------------------------------------
__device__ float3 VecBounce(const float3 &vec,const float3 &normal){
  const float pp=ProductScalar(vec,normal)/ProductScalar(normal,normal);
  const float ux=normal.x*pp;
  const float uy=normal.y*pp;
  const float uz=normal.z*pp;
  return(make_float3(vec.x-ux-ux,vec.y-uy-uy,vec.z-uz-uz));
}

//------------------------------------------------------------------------------
/// Devuelve el plano formado por un punto y un vector.
/// Returns the plane defined by a point and a vector.
//------------------------------------------------------------------------------
__device__ float4 PlanePtVec(const float3 &pt,const float3 &vec){
  const float3 v=VecUnitary(vec);//-No es necesario pero asi el modulo del vector no afecta al resultado de PointPlane().
  return(make_float4(v.x,v.y,v.z,-v.x*pt.x-v.y*pt.y-v.z*pt.z));
}

}


