//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsPips_ker.cu \brief Implements functions and CUDA kernels for PIPS calculation on GPU.

#include "JDsPips_ker.h"
#include <cfloat>
#include <math_constants.h>

namespace cupips{
#include "FunctionsBasic_iker.h"
#include "JCellSearch_iker.h"

//##############################################################################
//# Kernels for JDsPips.
//##############################################################################

//==============================================================================
/// Reduction using sum of unsigned values in shared memory for a warp.
/// Reduccion mediante suma de valores unsigned en memoria shared para un warp.
//==============================================================================
template <unsigned blockSize> __device__ void KerReduSumUintWarp(volatile unsigned* sudat,unsigned tid){
  if(blockSize>=64)sudat[tid]+=sudat[tid+32];
  if(blockSize>=32)sudat[tid]+=sudat[tid+16];
  if(blockSize>=16)sudat[tid]+=sudat[tid+8];
  if(blockSize>=8) sudat[tid]+=sudat[tid+4];
  if(blockSize>=4) sudat[tid]+=sudat[tid+2];
  if(blockSize>=2) sudat[tid]+=sudat[tid+1];
}

//------------------------------------------------------------------------------
/// Count number of real and checked neighbours.
//------------------------------------------------------------------------------
__device__ void KerInteractionNgBox(const unsigned &pini,const unsigned &pfin
  ,float kernelsize2,float poscellsize,const float4 &pscellp1,const float4 *poscell
  ,unsigned &numr,unsigned &numc)
{
  for(int p2=pini;p2<pfin;p2++){
    const float4 pscellp2=poscell[p2];
    const float rr2=cunsearch::Distance2(pscellp1,pscellp2,poscellsize);
    if(rr2<=kernelsize2 && rr2>=ALMOSTZERO)numr++;
  }
  numc+=(pfin-pini);
}

//------------------------------------------------------------------------------
/// Count number of real and checked neighbours in particle interacion.
//------------------------------------------------------------------------------
template <unsigned blockSize> __global__ void KerInteractionNg
  (unsigned nb,unsigned pinitb,unsigned nf,unsigned pinitf
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
  ,unsigned axis,unsigned cellcode,float kernelsize2,float poscellsize,const float4 *poscell,uint4 *res)
{
  extern __shared__ unsigned snrf[];
  unsigned *snrb=snrf+blockDim.x;
  unsigned *sncf=snrb+blockDim.x;
  unsigned *sncb=sncf+blockDim.x;
  const unsigned tid=threadIdx.x;
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<nb+nf){
    const unsigned p1=(p<nb? p+pinitb: p-nb+pinitf); //-Number of particle.
    unsigned numr=0,numc=0;
    //-Obtains position of particle p1.
    const float4 pscellp1=poscell[p1]; 

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::Initsp(dcell[p1],axis,cellcode,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin)KerInteractionNgBox(pini,pfin,kernelsize2,poscellsize,pscellp1,poscell,numr,numc);
    }
    
    //-Interaction with boundaries.
    if(p>=nb){//-When p particle is fluid.
      ini3-=cellfluid; fin3-=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin)KerInteractionNgBox(pini,pfin,kernelsize2,poscellsize,pscellp1,poscell,numr,numc);
      }
    }

    //-Stores numbers in shared memory.
    if(p>=nb){ snrf[tid]=numr; sncf[tid]=numc; snrb[tid]=0;    sncb[tid]=0;    }
    else{      snrf[tid]=0;    sncf[tid]=0;    snrb[tid]=numr; sncb[tid]=numc; }

  }
  else{
    snrf[tid]=0; snrb[tid]=0; sncf[tid]=0; sncb[tid]=0;
  }
  __syncthreads();
  if(blockSize>=512){ if(tid<256){ snrf[tid]+=snrf[tid+256]; snrb[tid]+=snrb[tid+256]; sncf[tid]+=sncf[tid+256]; sncb[tid]+=sncb[tid+256]; }  __syncthreads(); }
  if(blockSize>=256){ if(tid<128){ snrf[tid]+=snrf[tid+128]; snrb[tid]+=snrb[tid+128]; sncf[tid]+=sncf[tid+128]; sncb[tid]+=sncb[tid+128]; }  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) { snrf[tid]+=snrf[tid+ 64]; snrb[tid]+=snrb[tid+ 64]; sncf[tid]+=sncf[tid+ 64]; sncb[tid]+=sncb[tid+ 64]; }  __syncthreads(); }
  if(tid<32){ KerReduSumUintWarp<blockSize>(snrf,tid);  KerReduSumUintWarp<blockSize>(snrb,tid); KerReduSumUintWarp<blockSize>(sncf,tid);  KerReduSumUintWarp<blockSize>(sncb,tid); }
  if(tid==0)res[blockIdx.x]=make_uint4(snrf[0],snrb[0],sncf[0],sncb[0]);
}

//==============================================================================
/// Count number of real and checked neighbours in particle interacion.
//==============================================================================
void InteractionNg_1st(unsigned nb,unsigned pinitb,unsigned nf,unsigned pinitf
  ,const StDivDataGpu &dvd,const unsigned *dcell,const float4 *poscell,uint4 *res
  ,cudaStream_t stm)
{
  const unsigned BSIZE=256;
  const unsigned np=nb+nf;
  if(np){
    const unsigned shmem=sizeof(unsigned)*4*BSIZE;
    dim3 sgrid=GetSimpleGridSize(np,BSIZE);
    KerInteractionNg <BSIZE> <<<sgrid,BSIZE,shmem,stm>>>(nb,pinitb,nf,pinitf
      ,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,dcell
      ,dvd.axis,dvd.domcellcode,dvd.kernelsize2,dvd.poscellsize,poscell,res);
  }
}

//==============================================================================
/// Returns number of positions to store results of InteractionNg_1st().
//==============================================================================
unsigned InteractionNgSize_1st(unsigned np){
  const unsigned BSIZE=256;
  dim3 sgrid=GetSimpleGridSize(np,BSIZE);
  return(sgrid.x);
}

//==============================================================================
/// Accumulates the sum of n values of array dat[], storing the result in 
/// the beginning of res[].(Many positions of res[] are used as blocks, 
/// storing the final result in res[0]).
///
/// Acumula la suma de n valores del vector dat[], guardando el resultado al 
/// principio de res[] (Se usan tantas posiciones del res[] como bloques, 
/// quedando el resultado final en res[0]).
//==============================================================================
template <unsigned blockSize> __global__ void KerReduSumUint4(unsigned n,unsigned ini
  ,const uint4 *dat,uint4 *res)
{
  extern __shared__ unsigned sfdatx[];
  unsigned *sfdaty=sfdatx+blockDim.x;
  unsigned *sfdatz=sfdaty+blockDim.x;
  unsigned *sfdatw=sfdatz+blockDim.x;
  const unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  uint4 value=(c<n? dat[c+ini]: make_uint4(0,0,0,0));
  sfdatx[tid]=value.x;
  sfdaty[tid]=value.y;
  sfdatz[tid]=value.z;
  sfdatw[tid]=value.w;
  __syncthreads();
  if(blockSize>=512){ if(tid<256){ sfdatx[tid]+=sfdatx[tid+256]; sfdaty[tid]+=sfdaty[tid+256]; sfdatz[tid]+=sfdatz[tid+256]; sfdatw[tid]+=sfdatw[tid+256]; }  __syncthreads(); }
  if(blockSize>=256){ if(tid<128){ sfdatx[tid]+=sfdatx[tid+128]; sfdaty[tid]+=sfdaty[tid+128]; sfdatz[tid]+=sfdatz[tid+128]; sfdatw[tid]+=sfdatw[tid+128]; }  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) { sfdatx[tid]+=sfdatx[tid+ 64]; sfdaty[tid]+=sfdaty[tid+ 64]; sfdatz[tid]+=sfdatz[tid+ 64]; sfdatw[tid]+=sfdatw[tid+ 64]; }  __syncthreads(); }
  if(tid<32){ KerReduSumUintWarp<blockSize>(sfdatx,tid);  KerReduSumUintWarp<blockSize>(sfdaty,tid);  KerReduSumUintWarp<blockSize>(sfdatz,tid);  KerReduSumUintWarp<blockSize>(sfdatw,tid); }
  if(tid==0)res[blockIdx.x]=make_uint4(sfdatx[0],sfdaty[0],sfdatz[0],sfdatw[0]);
}


//==============================================================================
/// Count number of real and checked neighbours in particle interacion.
//==============================================================================
void InteractionNg_2nd(unsigned n,const uint4 *data,uint4 *res,cudaStream_t stm)
{
  if(n){
    const unsigned BSIZE=256;
    const unsigned shmem=sizeof(unsigned)*4*BSIZE;
    dim3 sgrid=GetSimpleGridSize(n,BSIZE);
    KerReduSumUint4 <BSIZE> <<<sgrid,BSIZE,shmem,stm>>> (n,0,data,res);
  }
}

//==============================================================================
/// Returns number of positions to store results of InteractionNg_2nd().
//==============================================================================
unsigned InteractionNgSize_2nd(unsigned n){
  const unsigned BSIZE=256;
  dim3 sgrid=GetSimpleGridSize(n,BSIZE);
  return(sgrid.x);
}

//==============================================================================
/// Reduction using sum of ullong values in shared memory for a warp.
/// Reduccion mediante suma de valores ullong en memoria shared para un warp.
//==============================================================================
template <unsigned blockSize> __device__ void KerReduSumUllongWarp(volatile ullong* sudat,unsigned tid){
  if(blockSize>=64)sudat[tid]+=sudat[tid+32];
  if(blockSize>=32)sudat[tid]+=sudat[tid+16];
  if(blockSize>=16)sudat[tid]+=sudat[tid+8];
  if(blockSize>=8) sudat[tid]+=sudat[tid+4];
  if(blockSize>=4) sudat[tid]+=sudat[tid+2];
  if(blockSize>=2) sudat[tid]+=sudat[tid+1];
}

//==============================================================================
/// Accumulates the sum of n values of array dat[], storing the result in 
/// the beginning of res[].(Many positions of res[] are used as blocks, 
/// storing the final result in res[0]).
///
/// Acumula la suma de n valores del vector dat[], guardando el resultado al 
/// principio de res[] (Se usan tantas posiciones del res[] como bloques, 
/// quedando el resultado final en res[0]).
//==============================================================================
template <unsigned blockSize> __global__ void KerReduSumUintlong4(unsigned n,unsigned ini
  ,const uint4 *dat,ullong *res)
{
  extern __shared__ ullong sldatx[];
  ullong *sldaty=sldatx+blockDim.x;
  ullong *sldatz=sldaty+blockDim.x;
  ullong *sldatw=sldatz+blockDim.x;
  const unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  const uint4 value=(c<n? dat[c+ini]: make_uint4(0,0,0,0));
  sldatx[tid]=value.x;
  sldaty[tid]=value.y;
  sldatz[tid]=value.z;
  sldatw[tid]=value.w;
  __syncthreads();
  if(blockSize>=512){ if(tid<256){ sldatx[tid]+=sldatx[tid+256]; sldaty[tid]+=sldaty[tid+256]; sldatz[tid]+=sldatz[tid+256]; sldatw[tid]+=sldatw[tid+256]; }  __syncthreads(); }
  if(blockSize>=256){ if(tid<128){ sldatx[tid]+=sldatx[tid+128]; sldaty[tid]+=sldaty[tid+128]; sldatz[tid]+=sldatz[tid+128]; sldatw[tid]+=sldatw[tid+128]; }  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) { sldatx[tid]+=sldatx[tid+ 64]; sldaty[tid]+=sldaty[tid+ 64]; sldatz[tid]+=sldatz[tid+ 64]; sldatw[tid]+=sldatw[tid+ 64]; }  __syncthreads(); }
  if(tid<32){ 
    KerReduSumUllongWarp<blockSize>(sldatx,tid);  
    KerReduSumUllongWarp<blockSize>(sldaty,tid);  
    KerReduSumUllongWarp<blockSize>(sldatz,tid);  
    KerReduSumUllongWarp<blockSize>(sldatw,tid); 
  }
  if(tid==0){
    unsigned offset=blockIdx.x*4;
    res[offset  ]=sldatx[0];
    res[offset+1]=sldaty[0];
    res[offset+2]=sldatz[0];
    res[offset+3]=sldatw[0];
  }
}

//==============================================================================
/// Count number of real and checked neighbours in particle interacion.
//==============================================================================
void InteractionNg_3th(unsigned n,const uint4 *data,ullong *res,cudaStream_t stm)
{
  if(n){
    const unsigned BSIZE=256;
    const unsigned shmem=sizeof(ullong)*4*BSIZE;
    dim3 sgrid=GetSimpleGridSize(n,BSIZE);
    KerReduSumUintlong4 <BSIZE> <<<sgrid,BSIZE,shmem,stm>>> (n,0,data,res);
  }
}

//==============================================================================
/// Returns number of positions to store results of InteractionNg_3th().
//==============================================================================
unsigned InteractionNgSize_3th(unsigned n){
  const unsigned BSIZE=256;
  dim3 sgrid=GetSimpleGridSize(n,BSIZE);
  return(sgrid.x);
}



}


