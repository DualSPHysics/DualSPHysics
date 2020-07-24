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

/// \file JWaveOrder2_ker.cu \brief Implements functions and CUDA kernels for second order irregular wave generation.

#include "JWaveOrder2_ker.h"
#include "JReduSum_ker.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include <string>
#include <cstdio>
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"

namespace cuwave2{
#include "FunctionsBasic_iker.h"

//##############################################################################
//# Kernels for JWaveSpectrum.
//##############################################################################
//==============================================================================
/// Returns required size of auxiliary memory according to data to be processed.
/// Devuelve tamanho necesario de memoria auxiliar segun datos a procesar.
//==============================================================================
unsigned GetSizeAux(unsigned n){
  return(curedus::GetAuxSize_ReduSumDouble(n)*2);
}

//==============================================================================
/// Reduction using sum of double values in shared memory for a warp.
/// Reduccion mediante suma de valores double en memoria shared para un warp.
//==============================================================================
template <unsigned blockSize> __device__ void KerReduSumDoubleWarp(volatile double* sdat,unsigned tid){
  if(blockSize>=64)sdat[tid]+=sdat[tid+32];
  if(blockSize>=32)sdat[tid]+=sdat[tid+16];
  if(blockSize>=16)sdat[tid]+=sdat[tid+8];
  if(blockSize>=8) sdat[tid]+=sdat[tid+4];
  if(blockSize>=4) sdat[tid]+=sdat[tid+2];
  if(blockSize>=2) sdat[tid]+=sdat[tid+1];
}

//==============================================================================
/// Calculates component of 2nd order for the piston position in irregular generation.
/// Calcula componente de 2nd orden para la posicion de piston en generacion irregular. 
//==============================================================================
template <unsigned blockSize> __global__ void KerCalcPosition(unsigned n,double time
  ,const double *dnm,const double2 *coefx,double *res)
{
  extern __shared__ double sdat[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  double e2=0;
  if(c<n){
    const double v=dnm[c]*time;
    const double2 coef=coefx[c];
    e2=coef.x*sin(v) - coef.y*cos(v);
  }
  sdat[tid]=e2;
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sdat[tid]+=sdat[tid+256];  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sdat[tid]+=sdat[tid+128];  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sdat[tid]+=sdat[tid+64];   __syncthreads(); }
  if(tid<32)KerReduSumDoubleWarp<blockSize>(sdat,tid);
  if(tid==0)res[blockIdx.x]=sdat[0];
}

//==============================================================================
/// Calculates component of 2nd order for the piston position in irregular generation.
/// Calcula componente de 2nd orden para la posicion de piston en generacion irregular. 
//==============================================================================
double CalcPosition(double time,unsigned n,const double *dnm,const double2 *coefx,double *aux)
{
  double res=0;
  if(n){
    const unsigned smemSize=sizeof(double)*WAVEBSIZE;
    const dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
    const unsigned n_blocks=sgrid.x*sgrid.y;
    KerCalcPosition<WAVEBSIZE> <<<sgrid,WAVEBSIZE,smemSize>>> (n,time,dnm,coefx,aux);
    res=curedus::ReduSumDouble(n_blocks,0,aux,aux+n_blocks);
  }
  return(res);
}

//==============================================================================
/// Calculates component of 2nd order for the elevation in irregular generation.
/// Calcula componente de 2nd orden para la elevacion en generacion irregular. 
//==============================================================================
template <unsigned blockSize> __global__ void KerCalcElevation(unsigned n,double time
  ,double x,const double4 *coefe,double *res)
{
  extern __shared__ double sdat[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  double eta2=0;
  if(c<n){
    const double4 coef=coefe[c];
    const double v=coef.x*time-x*coef.y;
    eta2=coef.z*cos(v) + coef.w*sin(v); 
  }
  sdat[tid]=eta2;
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sdat[tid]+=sdat[tid+256];  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sdat[tid]+=sdat[tid+128];  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sdat[tid]+=sdat[tid+64];   __syncthreads(); }
  if(tid<32)KerReduSumDoubleWarp<blockSize>(sdat,tid);
  if(tid==0)res[blockIdx.x]=sdat[0];
}

//==============================================================================
/// Calculates component of 2nd order for the elevation in irregular generation.
/// Calcula componente de 2nd orden para la elevacion en generacion irregular. 
//==============================================================================
double CalcElevation(double time,double x,unsigned n,const double4 *coefe,double *aux)
{
  double res=0;
  if(n){
    const unsigned smemSize=sizeof(double)*WAVEBSIZE;
    const dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
    const unsigned n_blocks=sgrid.x*sgrid.y;
    KerCalcElevation<WAVEBSIZE> <<<sgrid,WAVEBSIZE,smemSize>>> (n,time,x,coefe,aux);
    res=curedus::ReduSumDouble(n_blocks,0,aux,aux+n_blocks);
  }
  return(res);
}


/*:
////------------------------------------------------------------------------------
//// Calcula componente de 2nd orden para la posicion de piston en generacion irregular. 
////------------------------------------------------------------------------------
//__global__ void  KerCalcPosition1(unsigned n,double time
//  ,const double *dnm,const double2 *coefx,double *aux)
//{
//  unsigned c=blockIdx.x*blockDim.x + threadIdx.x; //-Number of interaction.
//  if(c<n){
//    const double v=dnm[c]*time;
//    const double2 coef=coefx[c];
//    double e2=coef.x*sin(v) - coef.y*cos(v);
//    aux[c]=e2;
//  }
//}
//
////==============================================================================
//// Calcula componente de 2nd orden para la posicion de piston en generacion irregular. 
////==============================================================================
//double CalcPosition1(double time,unsigned n,const double *dnm,const double2 *coefx,double *aux)
//{
//  double res=0;
//  if(n){
//    dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
//    KerCalcPosition1 <<<sgrid,WAVEBSIZE>>> (n,time,dnm,coefx,aux);
//    double *auxh=new double[n];
//    cudaMemcpy(auxh,aux,sizeof(double)*n,cudaMemcpyDeviceToHost);
//    for(unsigned c=0;c<n;c++)res+=auxh[c];
//    delete[] auxh;
//  }
//  return(res);
//}
//
////==============================================================================
//// Calcula componente de 2nd orden para la posicion de piston en generacion irregular. 
////==============================================================================
//double CalcPosition2(double time,unsigned n,const double *dnm,const double2 *coefx,double *aux)
//{
//  double res=0;
//  if(n){
//    dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
//    KerCalcPosition1 <<<sgrid,WAVEBSIZE>>> (n,time,dnm,coefx,aux);
//    res=curedu::ReduSumDouble(n,0,aux,aux+n);
//  }
//  return(res);
//}
:*/


}


