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


}


