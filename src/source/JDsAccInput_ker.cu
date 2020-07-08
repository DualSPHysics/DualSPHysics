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

/// \file JDsAccInput_ker.cu \brief Implements a set of functions and CUDA kernels for external forces (JDsAccInput) on GPU.

#include "JDsAccInput_ker.h"

namespace cuaccin{
#include "FunctionsBasic_iker.h"

//##############################################################################
//# Kernels for external forces (JDsAccInput).
//##############################################################################
//------------------------------------------------------
/// Adds variable forces to particle sets.
//------------------------------------------------------
__global__ void KerAddAccInputAng(unsigned n,unsigned pini,typecode codesel1,typecode codesel2
  ,float3 gravity,bool setgravity,double3 acclin,double3 accang,double3 centre,double3 velang,double3 vellin
  ,const typecode *code,const double2 *posxy,const double *posz,const float4 *velrhop,float3 *ace)
{
  const unsigned pp=blockIdx.x*blockDim.x + threadIdx.x;
  if(pp<n){
    const unsigned p=pp+pini;
    //Check if the current particle is part of the particle set by its Mk.
    const typecode tav=CODE_GetTypeAndValue(code[p]);
    if(codesel1<=tav && tav<=codesel2){
      const float3 accf=ace[p]; //-Gets the current particles acceleration value.
      double accx=accf.x,accy=accf.y,accz=accf.z;
      //-Adds linear acceleration.
      accx+=acclin.x;  accy+=acclin.y;  accz+=acclin.z;
      //-Subtract global gravity from the acceleration if it is set in the input file.
      if(!setgravity){
        accx-=gravity.x;  accy-=gravity.y;  accz-=gravity.z; 
      }

      //-Adds angular acceleration.
      const double2 rxy=posxy[p];
      const double dcx=rxy.x-centre.x;
      const double dcy=rxy.y-centre.y;
      const double dcz=posz[p]-centre.z;
      //-Get the current particle's velocity.
      const float4 rvel=velrhop[p];
      const double velx=rvel.x;
      const double vely=rvel.y;
      const double velz=rvel.z;

      //-Calculate angular acceleration ((Dw/Dt) x (r_i - r)) + (w x (w x (r_i - r))) + (2w x (v_i - v))
      //(Dw/Dt) x (r_i - r) (term1)
      accx+=(accang.y*dcz)-(accang.z*dcy);
      accy+=(accang.z*dcx)-(accang.x*dcz);
      accz+=(accang.x*dcy)-(accang.y*dcx);

      //-Centripetal acceleration (term2).
      //-First find w x (r_i - r)).
      const double innerx=(velang.y*dcz)-(velang.z*dcy);
      const double innery=(velang.z*dcx)-(velang.x*dcz);
      const double innerz=(velang.x*dcy)-(velang.y*dcx);
      //-Find w x inner.
      accx+=(velang.y*innerz)-(velang.z*innery);
      accy+=(velang.z*innerx)-(velang.x*innerz);
      accz+=(velang.x*innery)-(velang.y*innerx);

      //-Coriolis acceleration 2w x (v_i - v) (term3).
      accx+=((2.0*velang.y)*velz)-((2.0*velang.z)*(vely-vellin.y));
      accy+=((2.0*velang.z)*velx)-((2.0*velang.x)*(velz-vellin.z));
      accz+=((2.0*velang.x)*vely)-((2.0*velang.y)*(velx-vellin.x));

      //-Stores the new acceleration value.
      ace[p]=make_float3(float(accx),float(accy),float(accz));
    }
  }
}

//------------------------------------------------------
/// Adds variable forces to particle sets.
//------------------------------------------------------
__global__ void KerAddAccInputLin(unsigned n,unsigned pini,typecode codesel1,typecode codesel2
  ,float3 gravity,bool setgravity,double3 acclin,const typecode *code,float3 *ace)
{
  const unsigned pp=blockIdx.x*blockDim.x + threadIdx.x;
  if(pp<n){
    const unsigned p=pp+pini;
    //-Check if the current particle is part of the particle set by its Mk.
    const typecode tav=CODE_GetTypeAndValue(code[p]);
    if(codesel1<=tav && tav<=codesel2){
      const float3 accf=ace[p]; //-Gets the current particles acceleration value.
      double accx=accf.x,accy=accf.y,accz=accf.z;
      //-Adds linear acceleration.
      accx+=acclin.x;  accy+=acclin.y;  accz+=acclin.z;
      //-Subtract global gravity from the acceleration if it is set in the input file.
      if(!setgravity){
        accx-=gravity.x;  accy-=gravity.y;  accz-=gravity.z; 
      }
      //-Stores the new acceleration value.
      ace[p]=make_float3(float(accx),float(accy),float(accz));
    }
  }
}

//==================================================================================================
/// Adds external variable acceleration forces for particles according MK.
//==================================================================================================
void AddAccInput(unsigned n,unsigned pini,typecode codesel1,typecode codesel2
  ,tdouble3 acclin,tdouble3 accang,tdouble3 centre,tdouble3 velang,tdouble3 vellin,bool setgravity
  ,tfloat3 gravity,const typecode *code,const double2 *posxy,const double *posz,const float4 *velrhop,float3 *ace,cudaStream_t stm)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    const bool withaccang=(accang.x!=0 || accang.y!=0 || accang.z!=0);
    if(withaccang)KerAddAccInputAng <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,codesel1,codesel2,Float3(gravity),setgravity,Double3(acclin),Double3(accang),Double3(centre),Double3(velang),Double3(vellin),code,posxy,posz,velrhop,ace);
    else          KerAddAccInputLin <<<sgrid,SPHBSIZE,0,stm>>> (n,pini,codesel1,codesel2,Float3(gravity),setgravity,Double3(acclin),code,ace);
  }
}


}


