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

/// \file JSphShifting_ker.cu \brief Implements a set of functions and CUDA kernels for Shifting correction on GPU.

#include "JSphShifting_ker.h"
//#include "Functions.h"
//#include "FunctionsCuda.h"
#include <string>
#include <cstdio>
#include <cfloat>
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"

namespace cushift{
#include "FunctionsBasic_iker.h"
#include "FunctionsGeo3d_iker.h"

//##############################################################################
//# Kernels for JSphShifting.
//##############################################################################
//------------------------------------------------------------------------------
/// Select particles for shifting according min-max configuration
/// Selecciona particulas para shifting segun configuracion min-max.
//------------------------------------------------------------------------------
template<bool first,bool dbl> __global__ void KerInitGpuPosMax(unsigned n,unsigned pini
  ,const double3 pmin1,const double3 pmax1,const double3 pmin2,const double3 pmax2
  ,const double2* posxy,const double* posz,float4* shiftposfs)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    const unsigned p1=p+pini;
    const double2 pxy=posxy[p1];
    const double3 ps=make_double3(pxy.x,pxy.y,posz[p1]);
    if(cugeo::PointInMinMax(ps,pmin1,pmax1) || (dbl && cugeo::PointInMinMax(ps,pmin2,pmax2))){
      shiftposfs[p1]=make_float4(0,0,0,0);
    }
    else if(first)shiftposfs[p1]=make_float4(FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX);
  }
}

//==============================================================================
/// Select particles for shifting according min-max configuration
/// Selecciona particulas para shifting segun configuracion min-max.
//==============================================================================
void InitGpuPosMax(bool tfirst,bool tdbl,unsigned n,unsigned pini
  ,const tdouble3& pmin1,const tdouble3& pmax1,const tdouble3& pmin2,const tdouble3& pmax2
  ,const double2* posxy,const double* posz,float4* shiftposfs,cudaStream_t stm)
{
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    if(tfirst){ const bool first=true;
      if(tdbl)KerInitGpuPosMax<first,true > <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double3(pmin1),Double3(pmax1),Double3(pmin2),Double3(pmax2),posxy,posz,shiftposfs);
      else    KerInitGpuPosMax<first,false> <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double3(pmin1),Double3(pmax1),Double3(pmin2),Double3(pmax2),posxy,posz,shiftposfs);
    }else{      const bool first=false;
      if(tdbl)KerInitGpuPosMax<first,true > <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double3(pmin1),Double3(pmax1),Double3(pmin2),Double3(pmax2),posxy,posz,shiftposfs);
      else    KerInitGpuPosMax<first,false> <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double3(pmin1),Double3(pmax1),Double3(pmin2),Double3(pmax2),posxy,posz,shiftposfs);
    }
  }
}

//------------------------------------------------------------------------------
/// Select particles for shifting according planes configuration
/// Selecciona particulas para shifting segun configuracion de planos.
//------------------------------------------------------------------------------
template<bool first,bool dbl> __global__ void KerInitGpuPlanes(unsigned n,unsigned pini
  ,const double4 plax1,const double4 play1,const double4 plaz1,const double3 pladis1
  ,const double4 plax2,const double4 play2,const double4 plaz2,const double3 pladis2
  ,const double2* posxy,const double* posz,float4* shiftposfs)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    const unsigned p1=p+pini;
    const double2 pxy=posxy[p1];
    const double3 ps=make_double3(pxy.x,pxy.y,posz[p1]);
    if(cugeo::PlanesDomainCheck(ps,plax1,play1,plaz1,pladis1) || (dbl && cugeo::PlanesDomainCheck(ps,plax2,play2,plaz2,pladis2))){
      shiftposfs[p1]=make_float4(0,0,0,0);
    }
    else if(first)shiftposfs[p1]=make_float4(FLT_MAX,FLT_MAX,FLT_MAX,FLT_MAX);
  }
}

//==============================================================================
/// Select particles for shifting according planes configuration
/// Selecciona particulas para shifting segun configuracion de planos.
//==============================================================================
void InitGpuPlanes(bool tfirst,bool tdbl,unsigned n,unsigned pini
  ,const tplane3d& plax1,const tplane3d& play1,const tplane3d& plaz1,const tdouble3& pladis1
  ,const tplane3d& plax2,const tplane3d& play2,const tplane3d& plaz2,const tdouble3& pladis2
  ,const double2* posxy,const double* posz,float4* shiftposfs,cudaStream_t stm)
{
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    if(tfirst){ const bool first=true;
      if(tdbl)KerInitGpuPlanes<first,true > <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double4(plax1),Double4(play1),Double4(plaz1),Double3(pladis1),Double4(plax2),Double4(play2),Double4(plaz2),Double3(pladis2),posxy,posz,shiftposfs);
      else    KerInitGpuPlanes<first,false> <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double4(plax1),Double4(play1),Double4(plaz1),Double3(pladis1),Double4(plax2),Double4(play2),Double4(plaz2),Double3(pladis2),posxy,posz,shiftposfs);
    }else{      const bool first=false;
      if(tdbl)KerInitGpuPlanes<first,true > <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double4(plax1),Double4(play1),Double4(plaz1),Double3(pladis1),Double4(plax2),Double4(play2),Double4(plaz2),Double3(pladis2),posxy,posz,shiftposfs);
      else    KerInitGpuPlanes<first,false> <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,Double4(plax1),Double4(play1),Double4(plaz1),Double3(pladis1),Double4(plax2),Double4(play2),Double4(plaz2),Double3(pladis2),posxy,posz,shiftposfs);
    }
  }
}

//------------------------------------------------------------------------------
/// Computes final shifting for the particle position.
/// Calcula Shifting final para posicion de particulas.
//------------------------------------------------------------------------------
__global__ void KerRunShifting(unsigned n,unsigned pini,double dt
  ,double coefumagn,float shifttfs,double coeftfs,float maxdist
  ,const float4 *velrhop,float4 *shiftposfs)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    const unsigned p1=p+pini;
    float4 rs=shiftposfs[p1];
    if(rs.x!=FLT_MAX){
      const float4 rvel=velrhop[p1];
      const double vx=double(rvel.x);
      const double vy=double(rvel.y);
      const double vz=double(rvel.z);
      double umagn=coefumagn*sqrt(vx*vx+vy*vy+vz*vz);
      if(shifttfs){
        if(rs.w<shifttfs)umagn=0;
        else umagn*=(double(rs.w)-shifttfs)/coeftfs;
      }
      const float shiftdistx=float(double(rs.x)*umagn);
      const float shiftdisty=float(double(rs.y)*umagn);
      const float shiftdistz=float(double(rs.z)*umagn);
      rs.x=(shiftdistx<maxdist? shiftdistx: maxdist);
      rs.y=(shiftdisty<maxdist? shiftdisty: maxdist);
      rs.z=(shiftdistx<maxdist? shiftdistz: maxdist);
    }
    else rs=make_float4(0,0,0,0); //-Cancels shifting close to the boundaries. | Anula shifting por proximidad del contorno. 
    shiftposfs[p1]=rs;
  }
}

//==============================================================================
/// Computes final shifting for the particle position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void RunShifting(unsigned n,unsigned pini,double dt
  ,double coefumagn,float shifttfs,double coeftfs,float maxdist
  ,const float4 *velrhop,float4 *shiftposfs,cudaStream_t stm)
{
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerRunShifting <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,dt,coefumagn,shifttfs,coeftfs,maxdist,velrhop,shiftposfs);
  }
}


}


