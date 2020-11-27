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

/// \file JCellDivGpu_ker.cu \brief Implements functions and CUDA kernels to compute operations of the Neighbour List.

#include "JCellDivGpu_ker.h"
#include "DualSphDef.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include <float.h>
#include <cmath>
#include "JLog2.h"
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"

#pragma warning(disable : 4267) //Cancels "warning C4267: conversion from 'size_t' to 'int', possible loss of data"
#pragma warning(disable : 4244) //Cancels "warning C4244: conversion from 'unsigned __int64' to 'unsigned int', possible loss of data"
#include <thrust/device_vector.h>
#include <thrust/sort.h>

namespace cudiv{
#include "FunctionsBasic_iker.h"

//------------------------------------------------------------------------------
/// Reduction of shared memory values for a warp of KerPosLimitsRedu.
/// Reduccion de valores en memoria shared para un warp de KerPosLimitsRedu.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerUintLimitsWarpRedu(volatile unsigned* sp1,volatile unsigned* sp2,const unsigned &tid){
  if(blockSize>=64){
    const unsigned tid2=tid+32;
    sp1[tid]=min(sp1[tid],sp1[tid2]);
    sp2[tid]=max(sp2[tid],sp2[tid2]);
  }
  if(blockSize>=32){
    const unsigned tid2=tid+16;
    sp1[tid]=min(sp1[tid],sp1[tid2]);
    sp2[tid]=max(sp2[tid],sp2[tid2]);
  }
  if(blockSize>=16){
    const unsigned tid2=tid+8;
    sp1[tid]=min(sp1[tid],sp1[tid2]);
    sp2[tid]=max(sp2[tid],sp2[tid2]);
  }
  if(blockSize>=8){
    const unsigned tid2=tid+4;
    sp1[tid]=min(sp1[tid],sp1[tid2]);
    sp2[tid]=max(sp2[tid],sp2[tid2]);
  }
  if(blockSize>=4){
    const unsigned tid2=tid+2;
    sp1[tid]=min(sp1[tid],sp1[tid2]);
    sp2[tid]=max(sp2[tid],sp2[tid2]);
  }
  if(blockSize>=2){
    const unsigned tid2=tid+1;
    sp1[tid]=min(sp1[tid],sp1[tid2]);
    sp2[tid]=max(sp2[tid],sp2[tid2]);
  }
}

//------------------------------------------------------------------------------
/// Reduction of shared memory values for KerPosLimits
/// Reduccion de valores en memoria shared para KerPosLimits.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerUintLimitsRedu(unsigned* sp1,unsigned* sp2,const unsigned &tid,unsigned* results){
  __syncthreads();
  if(blockSize>=512){ 
    if(tid<256){
      sp1[tid]=min(sp1[tid],sp1[tid+256]);
      sp2[tid]=max(sp2[tid],sp2[tid+256]);
    }
    __syncthreads(); 
  }
  if(blockSize>=256){ 
    if(tid<128){
      sp1[tid]=min(sp1[tid],sp1[tid+128]);
      sp2[tid]=max(sp2[tid],sp2[tid+128]);
    }
    __syncthreads(); 
  }
  if(blockSize>=128){ 
    if(tid<64){
      sp1[tid]=min(sp1[tid],sp1[tid+64]);
      sp2[tid]=max(sp2[tid],sp2[tid+64]);
    }
    __syncthreads(); 
  }
  if(tid<32)KerUintLimitsWarpRedu<blockSize>(sp1,sp2,tid);
  if(tid==0){
    const unsigned nblocks=gridDim.x;
    unsigned cr=blockIdx.x;
    results[cr]=sp1[0]; cr+=nblocks;
    results[cr]=sp2[0];
  }
}

//==============================================================================
/// Reorders the values using RadixSort thrust.
/// Ordena valores usando RadixSort de thrust.
//==============================================================================
void Sort(unsigned* keys,unsigned* values,unsigned size,bool stable){
  if(size){
    thrust::device_ptr<unsigned> dev_keysg(keys);
    thrust::device_ptr<unsigned> dev_valuesg(values);
    if(stable)thrust::stable_sort_by_key(dev_keysg,dev_keysg+size,dev_valuesg);
    else thrust::sort_by_key(dev_keysg,dev_keysg+size,dev_valuesg);
  }
}

//------------------------------------------------------------------------------
/// Reduction of values in shared memory for a warp of KerPosLimitsRedu.
/// Reduccion de valores en memoria shared para un warp de KerPosLimitsRedu.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerPosLimitsWarpRedu(volatile float* spx1,volatile float* spy1,volatile float* spz1,volatile float* spx2,volatile float* spy2,volatile float* spz2,const unsigned &tid){
  if(blockSize>=64){
    const unsigned tid2=tid+32;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=32){
    const unsigned tid2=tid+16;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=16){
    const unsigned tid2=tid+8;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=8){
    const unsigned tid2=tid+4;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=4){
    const unsigned tid2=tid+2;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=2){
    const unsigned tid2=tid+1;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
}

//------------------------------------------------------------------------------
/// Reduction of shared memory values for KerPosLimits.
/// Reduccion de valores en memoria shared para KerPosLimits.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerPosLimitsRedu(float* spx1,float* spy1,float* spz1,float* spx2,float* spy2,float* spz2,const unsigned &tid,float* results){
  __syncthreads();
  if(blockSize>=512){ 
    if(tid<256){
      spx1[tid]=min(spx1[tid],spx1[tid+256]); spy1[tid]=min(spy1[tid],spy1[tid+256]); spz1[tid]=min(spz1[tid],spz1[tid+256]);  
      spx2[tid]=max(spx2[tid],spx2[tid+256]); spy2[tid]=max(spy2[tid],spy2[tid+256]); spz2[tid]=max(spz2[tid],spz2[tid+256]);  
    }
    __syncthreads(); 
  }
  if(blockSize>=256){ 
    if(tid<128){
      spx1[tid]=min(spx1[tid],spx1[tid+128]); spy1[tid]=min(spy1[tid],spy1[tid+128]); spz1[tid]=min(spz1[tid],spz1[tid+128]);  
      spx2[tid]=max(spx2[tid],spx2[tid+128]); spy2[tid]=max(spy2[tid],spy2[tid+128]); spz2[tid]=max(spz2[tid],spz2[tid+128]);  
    }
    __syncthreads(); 
  }
  if(blockSize>=128){ 
    if(tid<64){
      spx1[tid]=min(spx1[tid],spx1[tid+64]); spy1[tid]=min(spy1[tid],spy1[tid+64]); spz1[tid]=min(spz1[tid],spz1[tid+64]);  
      spx2[tid]=max(spx2[tid],spx2[tid+64]); spy2[tid]=max(spy2[tid],spy2[tid+64]); spz2[tid]=max(spz2[tid],spz2[tid+64]);  
    }
    __syncthreads(); 
  }
  if(tid<32)KerPosLimitsWarpRedu<blockSize>(spx1,spy1,spz1,spx2,spy2,spz2,tid);
  if(tid==0){
    const unsigned nblocks=gridDim.x;
    unsigned cr=blockIdx.x;
    results[cr]=spx1[0]; cr+=nblocks;
    results[cr]=spy1[0]; cr+=nblocks;
    results[cr]=spz1[0]; cr+=nblocks;
    results[cr]=spx2[0]; cr+=nblocks;
    results[cr]=spy2[0]; cr+=nblocks;
    results[cr]=spz2[0];
  }
}


//------------------------------------------------------------------------------
/// Computes minimum and maximum position starting from the results of KerPosLimit.
/// Calcula posicion minima y maxima a partir de los resultados de KerPosLimit.
//------------------------------------------------------------------------------
template <unsigned int blockSize> __global__ void KerReduPosLimits(unsigned n,float* data,float *results)
{
  extern __shared__ float spx1[];
  float *spy1=spx1+blockDim.x;
  float *spz1=spy1+blockDim.x;
  float *spx2=spz1+blockDim.x;
  float *spy2=spx2+blockDim.x;
  float *spz2=spy2+blockDim.x;
  const unsigned tid=threadIdx.x;
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Value number.
  //-Loads values in shared memory.
  //-Carga valores en memoria shared.
  unsigned p2=p;
  spx1[tid]=(p<n? data[p2]: FLT_MAX);  p2+=n;
  spy1[tid]=(p<n? data[p2]: FLT_MAX);  p2+=n;
  spz1[tid]=(p<n? data[p2]: FLT_MAX);  p2+=n;
  spx2[tid]=(p<n? data[p2]: -FLT_MAX); p2+=n;
  spy2[tid]=(p<n? data[p2]: -FLT_MAX); p2+=n;
  spz2[tid]=(p<n? data[p2]: -FLT_MAX);
  __syncthreads();
  //-Reduction of values in shared memory.
  //-Reduce valores de memoria shared.
  KerPosLimitsRedu<blockSize>(spx1,spy1,spz1,spx2,spy2,spz2,tid,results);
}

//==============================================================================
/// Reduction of position limits starting from results[].
/// In results[] each block stores xmin,ymin,zmin,xmax,ymax,zmax
/// grouped per block.
///
/// Reduce los limites de posicion a partir de results[].
/// En results[] cada bloque graba xmin,ymin,zmin,xmax,ymax,zmax agrupando por
/// bloque.
//==============================================================================
void ReduPosLimits(unsigned nblocks,float *aux,tfloat3 &pmin,tfloat3 &pmax){
  unsigned n=nblocks;
  const unsigned smemSize=DIVBSIZE*sizeof(float)*6;
  dim3 sgrid=GetSimpleGridSize(n,DIVBSIZE);
  unsigned n_blocks=sgrid.x*sgrid.y;
  //:printf("n:%d  n_blocks:%d]\n",n,n_blocks);
  float *dat=aux;
  float *res=aux+(n_blocks*6);
  while(n>1){
    //:printf("##>ReduMaxF n:%d  n_blocks:%d]\n",n,n_blocks);
    //:printf("##>ReduMaxF>sgrid=(%d,%d,%d)\n",sgrid.x,sgrid.y,sgrid.z);
    KerReduPosLimits<DIVBSIZE><<<sgrid,DIVBSIZE,smemSize>>>(n,dat,res);
    //fcuda::Check_CudaError("#>ReduMaxF Fallo en KerReduMaxF.");
    n=n_blocks;
    sgrid=GetSimpleGridSize(n,DIVBSIZE);  
    n_blocks=sgrid.x*sgrid.y;
    float* x=dat; dat=res; res=x;
  }
  float resf[6];
  cudaMemcpy(resf,dat,sizeof(float)*6,cudaMemcpyDeviceToHost);
  //fcuda::Check_CudaError("#>ReduMaxF Fallo en cudaMemcpy.");
  pmin=TFloat3(resf[0],resf[1],resf[2]);
  pmax=TFloat3(resf[3],resf[4],resf[5]);
}

//------------------------------------------------------------------------------
/// Reduction of shared memory values for a warp of KerLimitsCellRedu.
/// Reduccion de valores en memoria shared para un warp de KerLimitsCellRedu.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerLimitsCellWarpRedu(volatile unsigned* spx1,volatile unsigned* spy1,volatile unsigned* spz1,volatile unsigned* spx2,volatile unsigned* spy2,volatile unsigned* spz2,const unsigned &tid){
  if(blockSize>=64){
    const unsigned tid2=tid+32;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=32){
    const unsigned tid2=tid+16;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=16){
    const unsigned tid2=tid+8;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=8){
    const unsigned tid2=tid+4;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=4){
    const unsigned tid2=tid+2;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
  if(blockSize>=2){
    const unsigned tid2=tid+1;
    spx1[tid]=min(spx1[tid],spx1[tid2]); spy1[tid]=min(spy1[tid],spy1[tid2]); spz1[tid]=min(spz1[tid],spz1[tid2]);
    spx2[tid]=max(spx2[tid],spx2[tid2]); spy2[tid]=max(spy2[tid],spy2[tid2]); spz2[tid]=max(spz2[tid],spz2[tid2]);
  }
}

//------------------------------------------------------------------------------
/// Reduction of shared memory values for KerLimitsCell.
/// Reduccion de valores en memoria shared para KerLimitsCell.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerLimitsCellRedu(unsigned cellcode,unsigned* spx1,unsigned* spy1,unsigned* spz1,unsigned* spx2,unsigned* spy2,unsigned* spz2,const unsigned &tid,unsigned* results){
  __syncthreads();
  if(blockSize>=512){ 
    if(tid<256){
      spx1[tid]=min(spx1[tid],spx1[tid+256]); spy1[tid]=min(spy1[tid],spy1[tid+256]); spz1[tid]=min(spz1[tid],spz1[tid+256]);  
      spx2[tid]=max(spx2[tid],spx2[tid+256]); spy2[tid]=max(spy2[tid],spy2[tid+256]); spz2[tid]=max(spz2[tid],spz2[tid+256]);  
    }
    __syncthreads(); 
  }
  if(blockSize>=256){ 
    if(tid<128){
      spx1[tid]=min(spx1[tid],spx1[tid+128]); spy1[tid]=min(spy1[tid],spy1[tid+128]); spz1[tid]=min(spz1[tid],spz1[tid+128]);  
      spx2[tid]=max(spx2[tid],spx2[tid+128]); spy2[tid]=max(spy2[tid],spy2[tid+128]); spz2[tid]=max(spz2[tid],spz2[tid+128]);  
    }
    __syncthreads(); 
  }
  if(blockSize>=128){ 
    if(tid<64){
      spx1[tid]=min(spx1[tid],spx1[tid+64]); spy1[tid]=min(spy1[tid],spy1[tid+64]); spz1[tid]=min(spz1[tid],spz1[tid+64]);  
      spx2[tid]=max(spx2[tid],spx2[tid+64]); spy2[tid]=max(spy2[tid],spy2[tid+64]); spz2[tid]=max(spz2[tid],spz2[tid+64]);  
    }
    __syncthreads(); 
  }
  if(tid<32)KerLimitsCellWarpRedu<blockSize>(spx1,spy1,spz1,spx2,spy2,spz2,tid);
  if(tid==0){
    const unsigned nblocks=gridDim.x;
    unsigned cr=blockIdx.x;
    results[cr]=PC__Cell(cellcode,spx1[0],spy1[0],spz1[0]);  cr+=nblocks;
    results[cr]=PC__Cell(cellcode,spx2[0],spy2[0],spz2[0]);
  }
}

//------------------------------------------------------------------------------
/// Computes minimum and maximum postion startibg from the results of KerPosLimit.
/// Calcula posicion minima y maxima a partir de los resultados de KerPosLimit.
//------------------------------------------------------------------------------
template <unsigned int blockSize> __global__ void KerLimitsCellReduBase(unsigned cellcode,unsigned n,unsigned* data,unsigned *results)
{
  extern __shared__ unsigned scx1[];
  unsigned *scy1=scx1+blockDim.x;
  unsigned *scz1=scy1+blockDim.x;
  unsigned *scx2=scz1+blockDim.x;
  unsigned *scy2=scx2+blockDim.x;
  unsigned *scz2=scy2+blockDim.x;
  const unsigned tid=threadIdx.x;
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Value number.
  //-Loads values in shared memory.
  //-Carga valores en memoria shared.
  unsigned p2=p;
  const unsigned celmin=(p<n? data[p2]: UINT_MAX);  p2+=n;
  const unsigned celmax=(p<n? data[p2]: 0);
  scx1[tid]=PC__Cellx(cellcode,celmin);
  scy1[tid]=PC__Celly(cellcode,celmin);
  scz1[tid]=PC__Cellz(cellcode,celmin);
  scx2[tid]=PC__Cellx(cellcode,celmax);
  scy2[tid]=PC__Celly(cellcode,celmax);
  scz2[tid]=PC__Cellz(cellcode,celmax);
  __syncthreads();
  //-Reduction of shared memory values.
  //-Reduce valores de memoria shared.
  KerLimitsCellRedu<blockSize>(cellcode,scx1,scy1,scz1,scx2,scy2,scz2,tid,results);
}

//==============================================================================
/// Reduction of cell limits starting from results[].
/// In results[] each block stores cxmin,cymin,czmin,cxmax,cymax,czmax encodes
/// the values as cells in 2 unsigned and groups them per block.
///
/// Reduce los limites de celdas a partir de results[].
/// En results[] cada bloque graba cxmin,cymin,czmin,cxmax,cymax,czmax codificando
/// los valores como celdas en 2 unsigned y agrupando por bloque.
//==============================================================================
void LimitsCellRedu(unsigned cellcode,unsigned nblocks,unsigned *aux
  ,tuint3 &celmin,tuint3 &celmax)
{
  unsigned n=nblocks;
  const unsigned smemSize=DIVBSIZE*sizeof(float)*6;
  dim3 sgrid=GetSimpleGridSize(n,DIVBSIZE);
  unsigned n_blocks=sgrid.x*sgrid.y;
  //:printf("n:%d  n_blocks:%d]\n",n,n_blocks);
  unsigned *dat=aux;
  unsigned *res=aux+(n_blocks*2); //-Minimum and maximum value.
  while(n>1){
    //:printf("##>ReduMaxF n:%d  n_blocks:%d]\n",n,n_blocks);
    //:printf("##>ReduMaxF>sgrid=(%d,%d,%d)\n",sgrid.x,sgrid.y,sgrid.z);
    KerLimitsCellReduBase<DIVBSIZE><<<sgrid,DIVBSIZE,smemSize>>>(cellcode,n,dat,res);
    //fcuda::Check_CudaError("#>ReduMaxF Fallo en KerReduMaxF.");
    n=n_blocks;
    sgrid=GetSimpleGridSize(n,DIVBSIZE);  
    n_blocks=sgrid.x*sgrid.y;
    unsigned* x=dat; dat=res; res=x;
  }
  unsigned resf[6];
  cudaMemcpy(resf,dat,sizeof(unsigned)*2,cudaMemcpyDeviceToHost);
  //fcuda::Check_CudaError("#>ReduMaxF Fallo en cudaMemcpy.");
  celmin=TUint3(PC__Cellx(cellcode,resf[0]),PC__Celly(cellcode,resf[0]),PC__Cellz(cellcode,resf[0]));
  celmax=TUint3(PC__Cellx(cellcode,resf[1]),PC__Celly(cellcode,resf[1]),PC__Cellz(cellcode,resf[1]));
}

//------------------------------------------------------------------------------
/// Computes minimun and maximum cell for valid particles.
/// The excluded particles are already marked in code[].
/// In case of having no valid particles the minimum value igreater than the maximum.
/// In results[], each block stores cxmin,cymin,czmin,cxmax,cymax,czmax encodes
/// the values as cells in 2 unsigned and groups them per block.
///
/// Calcula celda minima y maxima de las particulas validas.
/// Las particulas excluidas ya estan marcadas en code[].
/// En caso de no haber ninguna particula valida el minimo sera mayor que el maximo.
/// En results[] cada bloque graba cxmin,cymin,czmin,cxmax,cymax,czmax codificando
/// los valores como celdas en 2 unsigned y agrupando por bloque.
//------------------------------------------------------------------------------
template <unsigned int blockSize> __global__ void KerLimitsCell(unsigned n,unsigned pini
  ,unsigned cellcode,const unsigned *dcell,const typecode *code,unsigned *results)
{
  extern __shared__ unsigned scx1[];
  unsigned *scy1=scx1+blockDim.x;
  unsigned *scz1=scy1+blockDim.x;
  unsigned *scx2=scz1+blockDim.x;
  unsigned *scy2=scx2+blockDim.x;
  unsigned *scz2=scy2+blockDim.x;
  const unsigned tid=threadIdx.x;
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  //-Loads shared memory values.
  if(p<n){
    const unsigned pp=p+pini;
    const unsigned rcell=dcell[pp];
    const unsigned cx=PC__Cellx(cellcode,rcell);
    const unsigned cy=PC__Celly(cellcode,rcell);
    const unsigned cz=PC__Cellz(cellcode,rcell);
    if(CODE_GetSpecialValue(code[pp])<CODE_OUTIGNORE){ //-Particle not excluded | Particula no excluida.
      scx1[tid]=cx; scy1[tid]=cy; scz1[tid]=cz;
      scx2[tid]=cx; scy2[tid]=cy; scz2[tid]=cz;
    }
    else{
      scx1[tid]=UINT_MAX; scy1[tid]=UINT_MAX;  scz1[tid]=UINT_MAX;
      scx2[tid]=0;        scy2[tid]=0;         scz2[tid]=0;
    }
  }
  else{
    scx1[tid]=UINT_MAX; scy1[tid]=UINT_MAX;  scz1[tid]=UINT_MAX;
    scx2[tid]=0;        scy2[tid]=0;         scz2[tid]=0;
  }
  __syncthreads();
  //-Reduction of shared memory values.
  KerLimitsCellRedu<blockSize>(cellcode,scx1,scy1,scz1,scx2,scy2,scz2,tid,results);
}

//==============================================================================
/// Computes minimun and maximum cell for valid particles.
/// The excluded particles are already marked in code[].
/// In case of having no valid particles the minimum value igreater than the maximum.
/// In results[], each block stores cxmin,cymin,czmin,cxmax,cymax,czmax encodes
/// the values as cells in 2 unsigned and groups them per block.
///
/// Calcula celda minima y maxima de las particulas validas.
/// Las particulas excluidas ya estan marcadas en code[].
/// En caso de no haber ninguna particula valida el minimo sera mayor que el maximo.
/// En results[] cada bloque graba cxmin,cymin,czmin,cxmax,cymax,czmax codificando
/// los valores como celdas en 2 unsigned y agrupando por bloque.
//==============================================================================
void LimitsCell(unsigned np,unsigned pini,unsigned cellcode,const unsigned *dcell
  ,const typecode *code,unsigned *aux,tuint3 &celmin,tuint3 &celmax)
{
  if(!np){//-Execution is canceled when no particles.
    celmin=TUint3(UINT_MAX);
    celmax=TUint3(0);
    return;
  }
  //:printf("[ReduMaxF ndata:%d  BLOCKSIZE:%d]\n",ndata,BLOCKSIZE);
  const unsigned smemSize=DIVBSIZE*sizeof(unsigned)*6;
  dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
  unsigned nblocks=sgrid.x*sgrid.y;
  KerLimitsCell<DIVBSIZE><<<sgrid,DIVBSIZE,smemSize>>>(np,pini,cellcode,dcell,code,aux);
  LimitsCellRedu(cellcode,nblocks,aux,celmin,celmax);
#ifdef DG_LimitsCell  //:delbeg:
  char cad[1024];
  sprintf(cad,"LimitsPos_%s> n:%u  pini:%u",(velrhop? "Fluid": "Bound"),np,pini); log->Print(cad);
  float4 *poscellh=new float4[np];
  byte *checkh=new byte[np];
  cudaMemcpy(poscellh,poscell+pini,sizeof(float4)*np,cudaMemcpyDeviceToHost);
  cudaMemcpy(checkh,check+pini,sizeof(byte)*np,cudaMemcpyDeviceToHost);
  tuint3 pminh=TUint3(UINT_MAX);
  tuint3 pmaxh=TUint3(0);
  for(unsigned p=0;p<np;p++)if(!checkh[p]){
    unsigned cell=*((unsigned *)&(poscellh[p].w));
    unsigned px=PC__Cellx(cellcode,cell),py=PC__Celly(cellcode,cell),pz=PC__Cellz(cellcode,cell);
    //sprintf(cad,"LimitsPos> cell[%u]=%u  (%u,%u,%u)",p,cell,px,py,pz); log->Print(cad);
    if(pminh.x>px)pminh.x=px;  if(pminh.y>py)pminh.y=py;  if(pminh.z>pz)pminh.z=pz;
    if(pmaxh.x<px)pmaxh.x=px;  if(pmaxh.y<py)pmaxh.y=py;  if(pmaxh.z<pz)pmaxh.z=pz;
  }
  delete[] poscellh;
  delete[] checkh;
  if(celmin.x!=pminh.x||celmin.y!=pminh.y||celmin.z!=pminh.z||celmax.x!=pmaxh.x||celmax.y!=pmaxh.y||celmax.z!=pmaxh.z){
    sprintf(cad,"LimitsPos> GPU pmin= (%u,%u,%u)  pmax= (%u,%u,%u)",celmin.x,celmin.y,celmin.z,celmax.x,celmax.y,celmax.z); log->Print(cad);
    sprintf(cad,"LimitsPos> CPU pminh=(%u,%u,%u)  pmaxh=(%u,%u,%u)",pminh.x,pminh.y,pminh.z,pmaxh.x,pmaxh.y,pmaxh.z); log->Print(cad);
    fun::Run_ExceptionStr("Error en LimitsPos()...");
  }
#endif  //:delend:
}

//------------------------------------------------------------------------------
/// Compute first and last particle for each cell.
/// Calcula particula inicial y final de cada celda.
//------------------------------------------------------------------------------
__global__ void KerCalcBeginEndCell(unsigned n,unsigned pini,const unsigned *cellpart,int2 *begcell)
{
  extern __shared__ unsigned scell[];    // [blockDim.x+1}
  const unsigned pt=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  const unsigned p=pt+pini;
  unsigned cel;
  if(pt<n){
    cel=cellpart[p];
    scell[threadIdx.x+1]=cel;
    if(pt&&!threadIdx.x)scell[0]=cellpart[p-1];
  }
  __syncthreads();
  if(pt<n){
    if(!pt||cel!=scell[threadIdx.x]){
      begcell[cel].x=p;
      if(pt)begcell[scell[threadIdx.x]].y=p;
    }
    if(pt==n-1)begcell[cel].y=p+1;
  }
}

//==============================================================================
/// Compute first and last particle for each cell.
/// Calcula particula inicial y final de cada celda.
//==============================================================================
void CalcBeginEndCell(bool full,unsigned np,unsigned npb,unsigned sizebegcell,unsigned cellfluid,const unsigned *cellpart,int2 *begcell){
  if(full)cudaMemset(begcell,0,sizeof(int2)*sizebegcell);
  else cudaMemset(begcell+cellfluid,0,sizeof(int2)*(sizebegcell-cellfluid));
  const unsigned pini=(full? 0: npb);
  const unsigned n=np-pini;
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,DIVBSIZE);
    KerCalcBeginEndCell <<<sgrid,DIVBSIZE,sizeof(unsigned)*(DIVBSIZE+1)>>> (n,pini,cellpart,begcell);
  }
}


//------------------------------------------------------------------------------
/// Reorders particle data according to idsort[].
/// Reordena datos de particulas segun idsort[].
//------------------------------------------------------------------------------
__global__ void KerSortDataParticles(unsigned n,unsigned pini,const unsigned *sortpart
  ,const unsigned *idp,const typecode *code,const unsigned *dcell,const double2 *posxy,const double *posz,const float4 *velrhop
  ,unsigned *idp2,typecode *code2,unsigned *dcell2,double2 *posxy2,double *posz2,float4 *velrhop2)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    const unsigned oldpos=(p<pini? p: sortpart[p]);
    idp2[p]=idp[oldpos];
    code2[p]=code[oldpos];
    dcell2[p]=dcell[oldpos];
    posxy2[p]=posxy[oldpos];
    posz2[p]=posz[oldpos];
    velrhop2[p]=velrhop[oldpos];
  }
}
//------------------------------------------------------------------------------
/// Reorders particle data according to sortpart[].
/// Reordena datos de particulas segun sortpart[].
//------------------------------------------------------------------------------
__global__ void KerSortDataParticles(unsigned n,unsigned pini,const unsigned *sortpart,const float4 *a,float4 *a2)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    const unsigned oldpos=(p<pini? p: sortpart[p]);
    a2[p]=a[oldpos];
  }
}
//------------------------------------------------------------------------------
/// Reorders particle data according to sortpart[].
/// Reordena datos de particulas segun sortpart[].
//------------------------------------------------------------------------------
__global__ void KerSortDataParticles(unsigned n,unsigned pini,const unsigned *sortpart,const float *a,const float *b,float *a2,float *b2)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    const unsigned oldpos=(p<pini? p: sortpart[p]);
    a2[p]=a[oldpos];
    b2[p]=b[oldpos];
  }
}
//------------------------------------------------------------------------------
/// Reorders particle data according to sortpart[].
/// Reordena datos de particulas segun sortpart[].
//------------------------------------------------------------------------------
__global__ void KerSortDataParticles(unsigned n,unsigned pini,const unsigned *sortpart,const double2 *a,const double *b,const float4 *c,double2 *a2,double *b2,float4 *c2)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    const unsigned oldpos=(p<pini? p: sortpart[p]);
    a2[p]=a[oldpos];
    b2[p]=b[oldpos];
    c2[p]=c[oldpos];
  }
}
//------------------------------------------------------------------------------
/// Reorders particle data according to sortpart[].
/// Reordena datos de particulas segun sortpart[].
//------------------------------------------------------------------------------
__global__ void KerSortDataParticles(unsigned n,unsigned pini,const unsigned *sortpart,const tsymatrix3f *a,tsymatrix3f *a2)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    const unsigned oldpos=(p<pini? p: sortpart[p]);
    a2[p]=a[oldpos];
  }
}
//------------------------------------------------------------------------------
/// Reorders particle data according to sortpart[].
/// Reordena datos de particulas segun sortpart[].
//------------------------------------------------------------------------------
__global__ void KerSortDataParticles(unsigned n,unsigned pini,const unsigned *sortpart,const float3 *a,float3 *a2)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    const unsigned oldpos=(p<pini? p: sortpart[p]);
    a2[p]=a[oldpos];
  }
}
//------------------------------------------------------------------------------
/// Reorders particle data according to sortpart[].
/// Reordena datos de particulas segun sortpart[].
//------------------------------------------------------------------------------
__global__ void KerSortDataParticles(unsigned n,unsigned pini,const unsigned *sortpart,const float *a,float *a2)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    const unsigned oldpos=(p<pini? p: sortpart[p]);
    a2[p]=a[oldpos];
  }
}

//==============================================================================
/// Reorders particle data according to sortpart.
/// Reordena datos de particulas segun sortpart.
//==============================================================================
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart
  ,const unsigned *idp,const typecode *code,const unsigned *dcell,const double2 *posxy,const double *posz,const float4 *velrhop
  ,unsigned *idp2,typecode *code2,unsigned *dcell2,double2 *posxy2,double *posz2,float4 *velrhop2)
{
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerSortDataParticles <<<sgrid,DIVBSIZE>>>(np,pini,sortpart,idp,code,dcell,posxy,posz,velrhop,idp2,code2,dcell2,posxy2,posz2,velrhop2);
  }
}
//==============================================================================
/// Reorders particle data according to sortpart.
/// Reordena datos de particulas segun sortpart.
//==============================================================================
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float4 *a,float4 *a2){
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerSortDataParticles <<<sgrid,DIVBSIZE>>>(np,pini,sortpart,a,a2);
  }
}
//==============================================================================
/// Reorders particle data according to sortpart.
/// Reordena datos de particulas segun sortpart.
//==============================================================================
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float *a,const float *b,float *a2,float *b2){
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerSortDataParticles <<<sgrid,DIVBSIZE>>>(np,pini,sortpart,a,b,a2,b2);
  }
}
//==============================================================================
/// Reorders particle data according to sortpart.
/// Reordena datos de particulas segun sortpart.
//==============================================================================
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const double2 *a,const double *b,const float4 *c,double2 *a2,double *b2,float4 *c2){
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerSortDataParticles <<<sgrid,DIVBSIZE>>>(np,pini,sortpart,a,b,c,a2,b2,c2);
  }
}
//==============================================================================
/// Reorders particle data according to sortpart.
/// Reordena datos de particulas segun sortpart.
//==============================================================================
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const tsymatrix3f *a,tsymatrix3f *a2){
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerSortDataParticles <<<sgrid,DIVBSIZE>>>(np,pini,sortpart,a,a2);
  }
}

//==============================================================================
/// Reorders particle data according to sortpart.
/// Reordena datos de particulas segun sortpart.
//==============================================================================
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float3 *a,float3 *a2){
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerSortDataParticles <<<sgrid,DIVBSIZE>>>(np,pini,sortpart,a,a2);
  }
}

//==============================================================================
/// Reorders particle data according to sortpart.
/// Reordena datos de particulas segun sortpart.
//==============================================================================
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float *a,float *a2){
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerSortDataParticles <<<sgrid,DIVBSIZE>>>(np,pini,sortpart,a,a2);
  }
}

//------------------------------------------------------------------------------
/// Compute minimum and maximum values starting from data[].
/// Calcula valores minimo y maximo a partir de data[].
//------------------------------------------------------------------------------
template <unsigned int blockSize> __global__ void KerReduUintLimits(unsigned n,unsigned* data,unsigned *results)
{
  extern __shared__ unsigned sp1[];
  unsigned *sp2=sp1+blockDim.x;
  const unsigned tid=threadIdx.x;
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Value number.
  //-Loads variables to shared memory.
  //-Carga valores en memoria shared.
  unsigned p2=p;
  sp1[tid]=(p<n? data[p2]: UINT_MAX);  p2+=n;
  sp2[tid]=(p<n? data[p2]: 0);
  __syncthreads();
  //-Reduction of shared memory values.
  //-Reduce valores de memoria shared.
  KerUintLimitsRedu<blockSize>(sp1,sp2,tid,results);
}

//==============================================================================
/// Reduce the limits of unsigned values from results[].
/// In results[] each block stores vmin,vamx gropued per block.
///
/// Reduce los limites de valores unsigned a partir de results[].
/// En results[] cada bloque graba vmin,vmax agrupando por bloque.
//==============================================================================
void ReduUintLimits(unsigned nblocks,unsigned *aux,unsigned &vmin,unsigned &vmax){
  unsigned n=nblocks;
  const unsigned smemSize=DIVBSIZE*sizeof(unsigned)*2;
  dim3 sgrid=GetSimpleGridSize(n,DIVBSIZE);
  unsigned n_blocks=sgrid.x*sgrid.y;
  //:printf("n:%d  n_blocks:%d]\n",n,n_blocks);
  unsigned *dat=aux;
  unsigned *res=aux+(n_blocks*2);
  while(n>1){
    //:printf("##>ReduMaxF n:%d  n_blocks:%d]\n",n,n_blocks);
    //:printf("##>ReduMaxF>sgrid=(%d,%d,%d)\n",sgrid.x,sgrid.y,sgrid.z);
    KerReduUintLimits<DIVBSIZE><<<sgrid,DIVBSIZE,smemSize>>>(n,dat,res);
    //fcuda::Check_CudaError("#>ReduMaxF Fallo en KerReduMaxF.");
    n=n_blocks;
    sgrid=GetSimpleGridSize(n,DIVBSIZE);  
    n_blocks=sgrid.x*sgrid.y;
    unsigned* x=dat; dat=res; res=x;
  }
  unsigned resf[2];
  cudaMemcpy(resf,dat,sizeof(unsigned)*2,cudaMemcpyDeviceToHost);
  //fcuda::Check_CudaError("#>ReduMaxF Fallo en cudaMemcpy.");
  vmin=resf[0];
  vmax=resf[1];
}

//------------------------------------------------------------------------------
/// Reduction of shared memory values for a warp of KerReduUintSum.
/// Reduccion de valores en memoria shared para un warp de KerReduUintSum.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerUintSumWarpRedu(volatile unsigned* sp1,const unsigned &tid){
  if(blockSize>=64)sp1[tid]+=sp1[tid+32];
  if(blockSize>=32)sp1[tid]+=sp1[tid+16];
  if(blockSize>=16)sp1[tid]+=sp1[tid+ 8];
  if(blockSize>= 8)sp1[tid]+=sp1[tid+ 4];
  if(blockSize>= 4)sp1[tid]+=sp1[tid+ 2];
  if(blockSize>= 2)sp1[tid]+=sp1[tid+ 1];
}

//------------------------------------------------------------------------------
/// Reduction of shared memory values for KerReduUintSum.
/// Reduccion de valores en memoria shared para KerReduUintSum.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerUintSumRedu(unsigned* sp1,const unsigned &tid,unsigned* results){
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sp1[tid]+=sp1[tid+256]; __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sp1[tid]+=sp1[tid+128]; __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sp1[tid]+=sp1[tid+64];  __syncthreads(); }
  if(tid<32)KerUintSumWarpRedu<blockSize>(sp1,tid);
  if(tid==0)results[blockIdx.x]=sp1[0];
}

//------------------------------------------------------------------------------
/// Returns the sum of the values contained in data[].
/// Devuelve la suma de los valores contenidos en data[].
//------------------------------------------------------------------------------
template <unsigned int blockSize> __global__ void KerReduUintSum(unsigned n,unsigned* data,unsigned *results)
{
  extern __shared__ unsigned sp1[];
  const unsigned tid=threadIdx.x;
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Value number.
  //-Loads values in shared memory.
  //-Carga valores en memoria shared.
  sp1[tid]=(p<n? data[p]: 0);
  __syncthreads();
  //-Reduce shared memory values.
  //-Reduce valores de memoria shared.
  KerUintSumRedu<blockSize>(sp1,tid,results);
}

//==============================================================================
/// Returns the sum of the values contained in aux[].
/// Devuelve la suma de los valores contenidos en aux[].
//==============================================================================
unsigned ReduUintSum(unsigned nblocks,unsigned *aux){
  unsigned n=nblocks;
  const unsigned smemSize=DIVBSIZE*sizeof(unsigned);
  dim3 sgrid=GetSimpleGridSize(n,DIVBSIZE);
  unsigned n_blocks=sgrid.x*sgrid.y;
  //:printf("n:%d  n_blocks:%d]\n",n,n_blocks);
  unsigned *dat=aux;
  unsigned *res=aux+(n_blocks);
  while(n>1){
    //:printf("##>ReduMaxF n:%d  n_blocks:%d]\n",n,n_blocks);
    //:printf("##>ReduMaxF>sgrid=(%d,%d,%d)\n",sgrid.x,sgrid.y,sgrid.z);
    KerReduUintSum<DIVBSIZE><<<sgrid,DIVBSIZE,smemSize>>>(n,dat,res);
    //fcuda::Check_CudaError("#>ReduMaxF Fallo en KerReduMaxF.");
    n=n_blocks;
    sgrid=GetSimpleGridSize(n,DIVBSIZE);  
    n_blocks=sgrid.x*sgrid.y;
    unsigned* x=dat; dat=res; res=x;
  }
  unsigned resf;
  cudaMemcpy(&resf,dat,sizeof(unsigned),cudaMemcpyDeviceToHost);
  //fcuda::Check_CudaError("#>ReduMaxF Fallo en cudaMemcpy.");
  return(resf);
}

/*:
////------------------------------------------------------------------------------
//// Devuelve rango de particulas en el rango de celdas indicadas.
////------------------------------------------------------------------------------
//template <unsigned blockSize> __global__ void KerGetRangeParticlesCells(unsigned ncel,unsigned ini,const int2 *begcell,unsigned *results)
//{ //torder{ORDER_XYZ=1,ORDER_YZX=2,ORDER_XZY=3} 
//  extern __shared__ unsigned sp1[];
//  unsigned *sp2=sp1+blockDim.x;
//  const unsigned tid=threadIdx.x;
//  const unsigned cel=blockIdx.x*blockDim.x + threadIdx.x; //-Number of cell
//  //-Carga valores en memoria shared.
//  if(cel<ncel){
//    int2 rg=begcell[cel+ini];
//    sp1[tid]=(rg.x<rg.y? unsigned(rg.x): UINT_MAX);
//    sp2[tid]=(rg.x<rg.y? unsigned(rg.y): 0);
//  }
//  else{
//    sp1[tid]=UINT_MAX;
//    sp2[tid]=0;
//  }
//  __syncthreads();
//  //-Reduce valores de memoria shared.
//  KerUintLimitsRedu<blockSize>(sp1,sp2,tid,results);
//}

////==============================================================================
//// Devuelve rango de particulas en el rango de celdas indicadas.
////==============================================================================
//void GetRangeParticlesCells(unsigned celini,unsigned celfin,const int2 *begcell,unsigned *aux,unsigned &pmin,unsigned &pmax){
//  unsigned ncel=celfin-celini;
//  if(!ncel){//-Si no hay celdas cancela proceso.
//    pmin=UINT_MAX; pmax=0;
//    return;
//  }
//  const unsigned smemSize=DIVBSIZE*sizeof(unsigned)*2;
//  dim3 sgrid=GetSimpleGridSize(ncel,DIVBSIZE);
//  unsigned nblocks=sgrid.x*sgrid.y;
//  KerGetRangeParticlesCells<DIVBSIZE><<<sgrid,DIVBSIZE,smemSize>>>(ncel,celini,begcell,aux);
//  ReduUintLimits(nblocks,aux,pmin,pmax,log);
//#ifdef DG_GetRangeParticlesCells
//  char cad[1024];
//  sprintf(cad,"GetRangeParticlesCells> ncel:%u  celini:%u",ncel,celini); log->Print(cad);
//  int2 *begcellh=new int2[ncel];
//  cudaMemcpy(begcellh,begcell+celini,sizeof(int2)*ncel,cudaMemcpyDeviceToHost);
//  unsigned pminh=UINT_MAX;
//  unsigned pmaxh=0;
//  for(unsigned p=0;p<ncel;p++){
//    unsigned x=unsigned(begcellh[p].x),y=unsigned(begcellh[p].y);
//    if(x<y){
//      if(pminh>x)pminh=x;
//      if(pmaxh<y)pmaxh=y;
//    }
//  }
//  delete[] begcellh;
//  if(pmin!=pminh||pmax!=pmaxh){
//    sprintf(cad,"GetRangeParticlesCells> GPU pmin= (%u)  pmax= (%u)",pmin,pmax); log->Print(cad);
//    sprintf(cad,"GetRangeParticlesCells> CPU pminh=(%u)  pmaxh=(%u)",pminh,pmaxh); log->Print(cad);
//    fun::Run_ExceptionStr("Error en GetRangeParticlesCells()...");
//  }
//#endif
//}

////------------------------------------------------------------------------------
//// Devuelve rango de particulas en el rango de celdas indicadas.
////------------------------------------------------------------------------------
//template <unsigned blockSize> __global__ void KerGetParticlesCells(unsigned ncel,unsigned ini,const int2 *begcell,unsigned *results)
//{ //torder{ORDER_XYZ=1,ORDER_YZX=2,ORDER_XZY=3} 
//  extern __shared__ unsigned sp1[];
//  const unsigned tid=threadIdx.x;
//  const unsigned cel=blockIdx.x*blockDim.x + threadIdx.x; //-Number of cell
//  //-Carga valores en memoria shared.
//  if(cel<ncel){
//    int2 rg=begcell[cel+ini];
//    sp1[tid]=(rg.y>rg.x? unsigned(rg.y-rg.x): 0);
//  }
//  else sp1[tid]=0;
//  __syncthreads();
//  //-Reduce valores de memoria shared.
//  KerUintSumRedu<blockSize>(sp1,tid,results);
//}

////==============================================================================
//// Devuelve numero de particulas en el rango de celdas indicadas.
////==============================================================================
//unsigned GetParticlesCells(unsigned celini,unsigned celfin,const int2 *begcell,unsigned *aux){
//  unsigned ncel=celfin-celini;
//  if(!ncel)return(0);//-Si no hay celdas cancela proceso.
//  const unsigned smemSize=DIVBSIZE*sizeof(unsigned);
//  dim3 sgrid=GetSimpleGridSize(ncel,DIVBSIZE);
//  unsigned nblocks=sgrid.x*sgrid.y;
//  KerGetParticlesCells<DIVBSIZE><<<sgrid,DIVBSIZE,smemSize>>>(ncel,celini,begcell,aux);
//  unsigned sum=ReduUintSum(nblocks,aux,log);
//#ifdef DG_GetParticlesCells
//  char cad[1024];
//  //sprintf(cad,"GetParticlesCells> ncel:%u  celini:%u",ncel,celini); log->PrintDbg(cad);
//  int2 *begcellh=new int2[ncel];
//  cudaMemcpy(begcellh,begcell+celini,sizeof(int2)*ncel,cudaMemcpyDeviceToHost);
//  unsigned sumh=0;
//  for(unsigned p=0;p<ncel;p++){
//    unsigned x=unsigned(begcellh[p].x),y=unsigned(begcellh[p].y);
//    if(y>x)sumh+=(y-x);
//  }
//  delete[] begcellh;
//  if(sum!=sumh){
//    sprintf(cad,"GetParticlesCells> GPU sum= (%u)",sum); log->PrintDbg(cad);
//    sprintf(cad,"GetParticlesCells> CPU sumh=(%u)",sumh); log->PrintDbg(cad);
//    fun::Run_ExceptionStr("Error en GetParticlesCells()...");
//  }
//#endif
//  return(sum);
//}
:*/

}


