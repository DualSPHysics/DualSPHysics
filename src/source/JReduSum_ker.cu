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

/// \file JReduSum_ker.cu \brief Implements functions and CUDA kernels for reduction using using the sum.

#include "JReduSum_ker.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include <cstdio>
#include <cfloat>
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"


namespace curedus{
#include "FunctionsBasic_iker.h"

//##############################################################################
//# Kernels for ReduSumDouble.
//##############################################################################

//==============================================================================
/// Reduction using sum of double values in shared memory for a warp.
/// Reduccion mediante suma de valores double en memoria shared para un warp.
//==============================================================================
template <unsigned blockSize> __device__ void KerReduSumDoubleWarp(volatile double* sddat,unsigned tid){
  if(blockSize>=64)sddat[tid]+=sddat[tid+32];
  if(blockSize>=32)sddat[tid]+=sddat[tid+16];
  if(blockSize>=16)sddat[tid]+=sddat[tid+8];
  if(blockSize>=8) sddat[tid]+=sddat[tid+4];
  if(blockSize>=4) sddat[tid]+=sddat[tid+2];
  if(blockSize>=2) sddat[tid]+=sddat[tid+1];
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
template <unsigned blockSize> __global__ void KerReduSumDouble(unsigned n,unsigned ini,const double *dat,double *res){
  extern __shared__ double sddat[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  sddat[tid]=(c<n? dat[c+ini]: 0);
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sddat[tid]+=sddat[tid+256];  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sddat[tid]+=sddat[tid+128];  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sddat[tid]+=sddat[tid+64];   __syncthreads(); }
  if(tid<32)KerReduSumDoubleWarp<blockSize>(sddat,tid);
  if(tid==0)res[blockIdx.x]=sddat[0];
}

//==============================================================================
/// Returns size of auxiliary memory for execution.
/// Devuelve el tamanho minimo de la memoria auxiliar necesaria segun ndata.
//==============================================================================
unsigned GetAuxSize_ReduSumDouble(unsigned ndata){
  const unsigned sizedata=1;
  dim3 sgrid=GetSimpleGridSize(ndata,REDUBSIZE);
  unsigned nblocks=sgrid.x*sgrid.y;  //-First iteration.
  sgrid=GetSimpleGridSize(nblocks,REDUBSIZE);
  unsigned nblocks2=sgrid.x*sgrid.y; //-Second iteration.
  return((nblocks+nblocks2)*sizedata);
}

//==============================================================================
/// Calculate sum of a vector, using resu [] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumDouble(ndata).
/// The execution can be synchronous or asynchronous and the result can be stored
/// in resu[], in pim1_sum or returned according to execution parameters.
/// Use ReduSumDouble() or ReduSumDoubleAsyn().
///
/// Calcula suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumDouble(ndata).
/// La ejecucion puede ser sincrona o asincrona y el resultado puede grabarse en
/// resu[], en pim1_sum o devolverse segun parametros de ejecucion.
/// Usar ReduSumDouble() o ReduSumDoubleAsyn().
//==============================================================================
double ReduSumDoubleBase(unsigned ndata,unsigned inidata,const double* data,double* resu,double *pim1_sum,cudaStream_t stm){
  double ret=0;
  if(ndata){
    unsigned n=ndata,ini=inidata;
    unsigned smemSize=sizeof(double)*REDUBSIZE;
    dim3 sgrid=GetSimpleGridSize(n,REDUBSIZE);
    unsigned nblocks=sgrid.x*sgrid.y;
    //:printf(">> n:%d  nblocks:%d]\n",n,nblocks);
    const double *dat=data;
    double *resu1=resu,*resu2=resu+nblocks;
    double *res=resu1;
    while(n>1){
      KerReduSumDouble<REDUBSIZE><<<sgrid,REDUBSIZE,smemSize,stm>>>(n,ini,dat,res);
      n=nblocks; ini=0;
      sgrid=GetSimpleGridSize(n,REDUBSIZE);  
      nblocks=sgrid.x*sgrid.y;
      if(n>1){
        dat=res; res=(dat==resu1? resu2: resu1); 
      }
    }
    //-Manages the result.
    const double *result=(ndata>1? res: data);
    if(!stm)cudaMemcpy(&ret,result,sizeof(double),cudaMemcpyDeviceToHost);
    else if(pim1_sum)cudaMemcpyAsync(pim1_sum,result,sizeof(double),cudaMemcpyDeviceToHost,stm);
    else if(res!=result)cudaMemcpyAsync(res,result,sizeof(double),cudaMemcpyDeviceToDevice,stm);
  }
  return(ret);
}

//==============================================================================
/// It performs the same calculations as ReduSumDouble but in CPU (for debugging).
/// Realiza los mismos calculos que ReduSumDouble pero en CPU (para debug).
//==============================================================================
double DgReduSumDouble(unsigned ndata,unsigned inidata,const double* datag){
  double res=0;
  if(ndata){
    //-Allocates CPU memory.
    double *data=new double[ndata];
    //-Gets data from GPU.
    cudaMemcpy(data,datag+inidata,sizeof(double)*ndata,cudaMemcpyDeviceToHost);
    fcuda::Check_CudaErroorFun("DgReduSumDouble");
    //-Process data.
    for(unsigned p=0;p<ndata;p++)res=res+data[p];
    //-Frees GPU memory.
    delete[] data; data=NULL;
  }
  return(res);
}

//==============================================================================
/// Returns the sum of a vector, using resu[] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumDouble(ndata).
/// Synchronous execution in stream 0.
///
/// Devuelve la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumDouble(ndata).
/// Ejecucion sincrona en el stream 0.
//==============================================================================
double ReduSumDouble(unsigned ndata,unsigned inidata,const double* data,double* resu){
  double ret=ReduSumDoubleBase(ndata,inidata,data,resu,NULL,NULL);
  #ifdef DG_curedus_ReduSumDouble
    #ifdef DG_curedus_Print
      printf("-=[DG_curedus_ReduSumDouble]=-\n");
    #endif
    double dgret=DgReduSumDouble(ndata,inidata,data);
    if(ret!=dgret)throw "ReduSumDouble: Results do not match.";
  #endif
  return(ret);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumDouble(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// the beginning of resu[].
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumDouble(ndata).
/// Ejecucion asincrona en el stream stm y el resultado final se graba al inicio
/// de resu[].
//==============================================================================
void ReduSumDoubleAsyn(unsigned ndata,unsigned inidata,const double* data,double* resu,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumDoubleAsyn().";
  if(!stm)throw "Error in parameters calling ReduSumDoubleAsyn().";
  ReduSumDoubleBase(ndata,inidata,data,resu,NULL,stm);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumDouble(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// pim1_sum (IT MUST BE PINNED MEMORY).
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumDouble(ndata).
/// Ejecucion asincrona en el stream stm y el resultado se graba en pim1_sum 
/// (QUE DEBE SER PINNED MEMORY).
//==============================================================================
void ReduSumDoubleAsyn(unsigned ndata,unsigned inidata,const double* data,double* resu,double *pim1_sum,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumDoubleAsyn().";
  if(!pim1_sum || !stm)throw "Error in parameters calling ReduSumDoubleAsyn().";
  ReduSumDoubleBase(ndata,inidata,data,resu,pim1_sum,stm);
}


//##############################################################################
//# Kernels for ReduSumFloat.
//##############################################################################

//==============================================================================
/// Reduction using sum of float values in shared memory for a warp.
/// Reduccion mediante suma de valores float en memoria shared para un warp.
//==============================================================================
template <unsigned blockSize> __device__ void KerReduSumFloatWarp(volatile float* sfdat,unsigned tid){
  if(blockSize>=64)sfdat[tid]+=sfdat[tid+32];
  if(blockSize>=32)sfdat[tid]+=sfdat[tid+16];
  if(blockSize>=16)sfdat[tid]+=sfdat[tid+8];
  if(blockSize>=8) sfdat[tid]+=sfdat[tid+4];
  if(blockSize>=4) sfdat[tid]+=sfdat[tid+2];
  if(blockSize>=2) sfdat[tid]+=sfdat[tid+1];
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
template <unsigned blockSize> __global__ void KerReduSumFloat(unsigned n,unsigned ini,const float *dat,float *res){
  extern __shared__ float sfdat[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  sfdat[tid]=(c<n? dat[c+ini]: 0);
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sfdat[tid]+=sfdat[tid+256];  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sfdat[tid]+=sfdat[tid+128];  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sfdat[tid]+=sfdat[tid+64];   __syncthreads(); }
  if(tid<32)KerReduSumFloatWarp<blockSize>(sfdat,tid);
  if(tid==0)res[blockIdx.x]=sfdat[0];
}

//==============================================================================
/// Returns size of auxiliary memory for execution.
/// Devuelve el tamanho minimo de la memoria auxiliar necesaria segun ndata.
//==============================================================================
unsigned GetAuxSize_ReduSumFloat(unsigned ndata){ return(GetAuxSize_ReduSumDouble(ndata)); }

//==============================================================================
/// Calculate sum of a vector, using resu [] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumDouble(ndata).
/// The execution can be synchronous or asynchronous and the result can be stored
/// in resu[], in pim1_sum or returned according to execution parameters.
/// Use ReduSumFloat() or ReduSumFloatAsyn().
///
/// Calcula suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// La ejecucion puede ser sincrona o asincrona y el resultado puede grabarse en
/// resu[], en pim1_sum o devolverse segun parametros de ejecucion.
/// Usar ReduSumFloat() o ReduSumFloatAsyn().
//==============================================================================
float ReduSumFloatBase(unsigned ndata,unsigned inidata,const float* data,float* resu,float *pim1_sum,cudaStream_t stm){
  float ret=0;
  if(ndata){
    unsigned n=ndata,ini=inidata;
    unsigned smemSize=sizeof(float)*REDUBSIZE;
    dim3 sgrid=GetSimpleGridSize(n,REDUBSIZE);
    unsigned nblocks=sgrid.x*sgrid.y;
    //:printf(">> n:%d  nblocks:%d]\n",n,nblocks);
    const float *dat=data;
    float *resu1=resu,*resu2=resu+nblocks;
    float *res=resu1;
    while(n>1){
      KerReduSumFloat<REDUBSIZE><<<sgrid,REDUBSIZE,smemSize,stm>>>(n,ini,dat,res);
      n=nblocks; ini=0;
      sgrid=GetSimpleGridSize(n,REDUBSIZE);  
      nblocks=sgrid.x*sgrid.y;
      if(n>1){
        dat=res; res=(dat==resu1? resu2: resu1); 
      }
    }
    //-Manages the result.
    const float *result=(ndata>1? res: data);
    if(!stm)cudaMemcpy(&ret,result,sizeof(float),cudaMemcpyDeviceToHost);
    else if(pim1_sum)cudaMemcpyAsync(pim1_sum,result,sizeof(float),cudaMemcpyDeviceToHost,stm);
    else if(res!=result)cudaMemcpyAsync(res,result,sizeof(float),cudaMemcpyDeviceToDevice,stm);
  }
  return(ret);
}

//==============================================================================
/// It performs the same calculations as ReduSumFloat but in CPU (for debugging).
/// Realiza los mismos calculos que ReduSumFloat pero en CPU (para debug).
//==============================================================================
float DgReduSumFloat(unsigned ndata,unsigned inidata,const float* datag){
  float res=0;
  if(ndata){
    //-Allocates CPU memory.
    float *data=new float[ndata];
    //-Gets data from GPU.
    cudaMemcpy(data,datag+inidata,sizeof(float)*ndata,cudaMemcpyDeviceToHost);
    fcuda::Check_CudaErroorFun("DgReduSumFloat");
    //-Process data.
    for(unsigned p=0;p<ndata;p++)res=res+data[p];
     //-Frees GPU memory.
    delete[] data; data=NULL;
  }
  return(res);
}

//==============================================================================
/// Returns the sum of a vector, using resu[] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumFloat(ndata).
/// Synchronous execution in stream 0.
///
/// Devuelve la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// Ejecucion sincrona en el stream 0.
//==============================================================================
float ReduSumFloat(unsigned ndata,unsigned inidata,const float* data,float* resu){
  float ret=ReduSumFloatBase(ndata,inidata,data,resu,NULL,NULL);
  #ifdef DG_curedus_ReduSumFloat
    #ifdef DG_curedus_Print
      printf("-=[DG_curedus_ReduSumFloat]=-\n");
    #endif
    float dgret=DgReduSumFloat(ndata,inidata,data);
    if(ret!=dgret && fabs(ret/dgret-1)>0.00001f){
      printf("-=[ %.8E == %.8E   dif:%.8f ]=-\n",ret,dgret,fabs(ret/dgret-1));
      throw "ReduSumFloat: Results do not match.";
    }
  #endif
  return(ret);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumFloat(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// the beginning of resu[].
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// Ejecucion asincrona en el stream stm y el resultado final se graba al inicio
/// de resu[].
//==============================================================================
void ReduSumFloatAsyn(unsigned ndata,unsigned inidata,const float* data,float* resu,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumFloatAsyn().";
  if(!stm)throw "Error in parameters calling ReduSumFloatAsyn().";
  ReduSumFloatBase(ndata,inidata,data,resu,NULL,stm);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumFloat(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// pim1_sum (IT MUST BE PINNED MEMORY).
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// Ejecucion asincrona en el stream stm y el resultado se graba en pim1_sum 
/// (QUE DEBE SER PINNED MEMORY).
//==============================================================================
void ReduSumFloatAsyn(unsigned ndata,unsigned inidata,const float* data,float* resu,float *pim1_sum,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumFloatAsyn().";
  if(!pim1_sum || !stm)throw "Error in parameters calling ReduSumFloatAsyn().";
  ReduSumFloatBase(ndata,inidata,data,resu,pim1_sum,stm);
}


//##############################################################################
//# Kernels for ReduSumUint.
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

//==============================================================================
/// Accumulates the sum of n values of array dat[], storing the result in 
/// the beginning of res[].(Many positions of res[] are used as blocks, 
/// storing the final result in res[0]).
///
/// Acumula la suma de n valores del vector dat[], guardando el resultado al 
/// principio de res[] (Se usan tantas posiciones del res[] como bloques, 
/// quedando el resultado final en res[0]).
//==============================================================================
template <unsigned blockSize> __global__ void KerReduSumUint(unsigned n,unsigned ini,const unsigned *dat,unsigned *res){
  extern __shared__ unsigned sudat[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  sudat[tid]=(c<n? dat[c+ini]: 0);
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sudat[tid]+=sudat[tid+256];  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sudat[tid]+=sudat[tid+128];  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sudat[tid]+=sudat[tid+64];   __syncthreads(); }
  if(tid<32)KerReduSumUintWarp<blockSize>(sudat,tid);
  if(tid==0)res[blockIdx.x]=sudat[0];
}

//==============================================================================
/// Returns size of auxiliary memory for execution.
/// Devuelve el tamanho minimo de la memoria auxiliar necesaria segun ndata.
//==============================================================================
unsigned GetAuxSize_ReduSumUint(unsigned ndata){ return(GetAuxSize_ReduSumDouble(ndata)); }

//==============================================================================
/// Calculate sum of a vector, using resu [] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumUint(ndata).
/// The execution can be synchronous or asynchronous and the result can be stored
/// in resu[], in pim1_sum or returned according to execution parameters.
/// Use ReduSumUint() or ReduSumUintAsyn().
///
/// Calcula suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumUint(ndata).
/// La ejecucion puede ser sincrona o asincrona y el resultado puede grabarse en
/// resu[], en pim1_sum o devolverse segun parametros de ejecucion.
/// Usar ReduSumUint() o ReduSumUintAsyn().
//==============================================================================
unsigned ReduSumUintBase(unsigned ndata,unsigned inidata,const unsigned* data,unsigned* resu,unsigned *pim1_sum,cudaStream_t stm){
  unsigned ret=0;
  if(ndata){
    unsigned n=ndata,ini=inidata;
    unsigned smemSize=sizeof(unsigned)*REDUBSIZE;
    dim3 sgrid=GetSimpleGridSize(n,REDUBSIZE);
    unsigned nblocks=sgrid.x*sgrid.y;
    //printf(">> n:%d  nblocks:%d]\n",n,nblocks);
    const unsigned *dat=data;
    unsigned *resu1=resu,*resu2=resu+nblocks;
    unsigned *res=resu1;
    while(n>1){
      KerReduSumUint<REDUBSIZE><<<sgrid,REDUBSIZE,smemSize,stm>>>(n,ini,dat,res);
      n=nblocks; ini=0;
      sgrid=GetSimpleGridSize(n,REDUBSIZE);  
      nblocks=sgrid.x*sgrid.y;
      if(n>1){
        dat=res; res=(dat==resu1? resu2: resu1); 
      }
    }
    //-Manages the result.
    const unsigned *result=(ndata>1? res: data);
    if(!stm)cudaMemcpy(&ret,result,sizeof(unsigned),cudaMemcpyDeviceToHost);
    else if(pim1_sum)cudaMemcpyAsync(pim1_sum,result,sizeof(unsigned),cudaMemcpyDeviceToHost,stm);
    else if(res!=result)cudaMemcpyAsync(res,result,sizeof(unsigned),cudaMemcpyDeviceToDevice,stm);
  }
  return(ret);
}

//==============================================================================
/// It performs the same calculations as ReduSumUint but in CPU (for debugging).
/// Realiza los mismos calculos que ReduSumUint pero en CPU (para debug).
//==============================================================================
unsigned DgReduSumUint(unsigned ndata,unsigned inidata,const unsigned* datag){
  unsigned res=0;
  if(ndata){
    //-Allocates CPU memory.
    unsigned *data=new unsigned[ndata];
    //-Gets data from GPU.
    cudaMemcpy(data,datag+inidata,sizeof(unsigned)*ndata,cudaMemcpyDeviceToHost);
    fcuda::Check_CudaErroorFun("DgReduSumUint");
    //-Process data.
    for(unsigned p=0;p<ndata;p++)res=res+data[p];
    //-Frees GPU memory.
    delete[] data; data=NULL;
  }
  return(res);
}

//==============================================================================
/// Returns the sum of a vector, using resu[] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumUint(ndata).
/// Synchronous execution in stream 0.
///
/// Devuelve la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumUint(ndata).
/// Ejecucion sincrona en el stream 0.
//==============================================================================
unsigned ReduSumUint(unsigned ndata,unsigned inidata,const unsigned* data,unsigned* resu){
  unsigned ret=ReduSumUintBase(ndata,inidata,data,resu,NULL,NULL);
  #ifdef DG_curedus_ReduSumUint
    #ifdef DG_curedus_Print
      printf("-=[DG_curedus_ReduSumUint]=-\n");
    #endif
    unsigned dgret=DgReduSumUint(ndata,inidata,data);
    if(ret!=dgret)throw "ReduSumUint: Results do not match.";
  #endif
  return(ret);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumUint(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// the beginning of resu[].
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumUint(ndata).
/// Ejecucion asincrona en el stream stm y el resultado final se graba al inicio
/// de resu[].
//==============================================================================
void ReduSumUintAsyn(unsigned ndata,unsigned inidata,const unsigned* data,unsigned* resu,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumUintAsyn().";
  if(!stm)throw "Error in parameters calling ReduSumUintAsyn().";
  ReduSumUintBase(ndata,inidata,data,resu,NULL,stm);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumUint(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// pim1_sum (IT MUST BE PINNED MEMORY).
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumUint(ndata).
/// Ejecucion asincrona en el stream stm y el resultado se graba en pim1_sum 
/// (QUE DEBE SER PINNED MEMORY).
//==============================================================================
void ReduSumUintAsyn(unsigned ndata,unsigned inidata,const unsigned* data,unsigned* resu,unsigned *pim1_sum,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumUintAsyn().";
  if(!pim1_sum || !stm)throw "Error in parameters calling ReduSumUintAsyn().";
  ReduSumUintBase(ndata,inidata,data,resu,pim1_sum,stm);
}


//##############################################################################
//# Kernels for ReduSumFloat3.
//##############################################################################

//==============================================================================
/// Accumulates the sum of n values of array dat[], storing the result in 
/// the beginning of res[].(Many positions of res[] are used as blocks, 
/// storing the final result in res[0]).
///
/// Acumula la suma de n valores del vector dat[], guardando el resultado al 
/// principio de res[] (Se usan tantas posiciones del res[] como bloques, 
/// quedando el resultado final en res[0]).
//==============================================================================
template <unsigned blockSize> __global__ void KerReduSumFloat3(unsigned n,unsigned ini,const float3 *dat,float3 *res){
  extern __shared__ float sfdatx[];
  float *sfdaty=sfdatx+blockDim.x;
  float *sfdatz=sfdaty+blockDim.x;
  const unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  float3 value=(c<n? dat[c+ini]: make_float3(0,0,0));
  sfdatx[tid]=value.x;
  sfdaty[tid]=value.y;
  sfdatz[tid]=value.z;
  __syncthreads();
  if(blockSize>=512){ if(tid<256){ sfdatx[tid]+=sfdatx[tid+256]; sfdaty[tid]+=sfdaty[tid+256]; sfdatz[tid]+=sfdatz[tid+256]; }  __syncthreads(); }
  if(blockSize>=256){ if(tid<128){ sfdatx[tid]+=sfdatx[tid+128]; sfdaty[tid]+=sfdaty[tid+128]; sfdatz[tid]+=sfdatz[tid+128]; }  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) { sfdatx[tid]+=sfdatx[tid+ 64]; sfdaty[tid]+=sfdaty[tid+ 64]; sfdatz[tid]+=sfdatz[tid+ 64]; }  __syncthreads(); }
  if(tid<32){ KerReduSumFloatWarp<blockSize>(sfdatx,tid);  KerReduSumFloatWarp<blockSize>(sfdaty,tid);  KerReduSumFloatWarp<blockSize>(sfdatz,tid); }
  if(tid==0)res[blockIdx.x]=make_float3(sfdatx[0],sfdaty[0],sfdatz[0]);
}

//==============================================================================
/// Returns size of auxiliary memory for execution.
/// Devuelve el tamanho minimo de la memoria auxiliar necesaria segun ndata.
//==============================================================================
unsigned GetAuxSize_ReduSumFloat3(unsigned ndata){ return(GetAuxSize_ReduSumDouble(ndata)); }

//==============================================================================
/// Calculate sum of a vector, using resu [] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumDouble(ndata).
/// The execution can be synchronous or asynchronous and the result can be stored
/// in resu[], in pim1_sum or returned according to execution parameters.
/// Use ReduSumFloat() or ReduSumFloatAsyn().
///
/// Calcula suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// La ejecucion puede ser sincrona o asincrona y el resultado puede grabarse en
/// resu[], en pim1_sum o devolverse segun parametros de ejecucion.
/// Usar ReduSumFloat() o ReduSumFloatAsyn().
//==============================================================================
float3 ReduSumFloat3Base(unsigned ndata,unsigned inidata,const float3* data,float3* resu,float3 *pim1_sum,cudaStream_t stm){
  float3 ret; ret.x=ret.y=ret.z=0;
  if(ndata){
    unsigned n=ndata,ini=inidata;
    unsigned smemSize=sizeof(float)*3*REDUBSIZE;
    dim3 sgrid=GetSimpleGridSize(n,REDUBSIZE);
    unsigned nblocks=sgrid.x*sgrid.y;
    //:printf(">> n:%d  nblocks:%d]\n",n,nblocks);
    const float3 *dat=data;
    float3 *resu1=resu,*resu2=resu+nblocks;
    float3 *res=resu1;
    while(n>1){
      KerReduSumFloat3<REDUBSIZE><<<sgrid,REDUBSIZE,smemSize,stm>>>(n,ini,dat,res);
      n=nblocks; ini=0;
      sgrid=GetSimpleGridSize(n,REDUBSIZE);  
      nblocks=sgrid.x*sgrid.y;
      if(n>1){
        dat=res; res=(dat==resu1? resu2: resu1); 
      }
    }
    //-Manages the result.
    const float3 *result=(ndata>1? res: data);
    if(!stm)cudaMemcpy(&ret,result,sizeof(float3),cudaMemcpyDeviceToHost);
    else if(pim1_sum)cudaMemcpyAsync(pim1_sum,result,sizeof(float3),cudaMemcpyDeviceToHost,stm);
    else if(res!=result)cudaMemcpyAsync(res,result,sizeof(float3),cudaMemcpyDeviceToDevice,stm);
  }
  return(ret);
}

//==============================================================================
/// It performs the same calculations as ReduSumFloat but in CPU (for debugging).
/// Realiza los mismos calculos que ReduSumFloat pero en CPU (para debug).
//==============================================================================
float3 DgReduSumFloat3(unsigned ndata,unsigned inidata,const float3* datag){
  float3 res; res.x=res.y=res.z=0;
  if(ndata){
    //-Allocates CPU memory.
    float3 *data=new float3[ndata];
    //-Gets data from GPU.
    cudaMemcpy(data,datag+inidata,sizeof(float3)*ndata,cudaMemcpyDeviceToHost);
    fcuda::Check_CudaErroorFun("DgReduSumFloat3");
    //-Process data.
    for(unsigned p=0;p<ndata;p++){
      res.x=res.x+data[p].x;
      res.y=res.y+data[p].y;
      res.z=res.z+data[p].z;
    }
     //-Frees GPU memory.
    delete[] data; data=NULL;
  }
  return(res);
}

//==============================================================================
/// Returns the sum of a vector, using resu[] as auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumFloat(ndata).
/// Synchronous execution in stream 0.
///
/// Devuelve la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// Ejecucion sincrona en el stream 0.
//==============================================================================
float3 ReduSumFloat3(unsigned ndata,unsigned inidata,const float3* data,float3* resu){
  float3 ret=ReduSumFloat3Base(ndata,inidata,data,resu,NULL,NULL);
  #ifdef DG_curedus_ReduSumFloat
    #ifdef DG_curedus_Print
      printf("-=[DG_curedus_ReduSumFloat3]=-\n");
    #endif
    float dgret=DgReduSumFloat3(ndata,inidata,data);
    if(ret!=dgret && (fabs(ret.x/dgret.x-1)>0.00001f) || fabs(ret.y/dgret.y-1)>0.00001f) || fabs(ret.z/dgret.z.-1)>0.00001f)){
      printf("-=[x %.8E == %.8E   dif:%.8f ]=-\n",ret.x,dgret.x,fabs(ret.x/dgret.x-1));
      printf("-=[y %.8E == %.8E   dif:%.8f ]=-\n",ret.y,dgret.y,fabs(ret.y/dgret.y-1));
      printf("-=[z %.8E == %.8E   dif:%.8f ]=-\n",ret.z,dgret.z,fabs(ret.z/dgret.z-1));
      throw "ReduSumFloat3: Results do not match.";
    }
  #endif
  return(ret);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumFloat(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// the beginning of resu[].
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// Ejecucion asincrona en el stream stm y el resultado final se graba al inicio
/// de resu[].
//==============================================================================
void ReduSumFloat3Asyn(unsigned ndata,unsigned inidata,const float3* data,float3* resu,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumFloat3Asyn().";
  if(!stm)throw "Error in parameters calling ReduSumFloat3Asyn().";
  ReduSumFloat3Base(ndata,inidata,data,resu,NULL,stm);
}

//==============================================================================
/// Calculates the sum of a vector, using resu[] as an auxiliary vector.
/// The size of resu[] must be >= a GetAuxSize_ReduSumFloat(ndata).
/// Asynchronous execution in the stream stm and the final result is stored at 
/// pim1_sum (IT MUST BE PINNED MEMORY).
///
/// Calcula la suma de un vector, usando resu[] como vector auxiliar.
/// El tamanho de resu[] debe ser >= a GetAuxSize_ReduSumFloat(ndata).
/// Ejecucion asincrona en el stream stm y el resultado se graba en pim1_sum 
/// (QUE DEBE SER PINNED MEMORY).
//==============================================================================
void ReduSumFloat3Asyn(unsigned ndata,unsigned inidata,const float3* data,float3* resu,float3 *pim1_sum,cudaStream_t stm){
  if(!ndata)throw "Number of values cannot be zero ReduSumFloat3Asyn().";
  if(!pim1_sum || !stm)throw "Error in parameters calling ReduSumFloat3Asyn().";
  ReduSumFloat3Base(ndata,inidata,data,resu,pim1_sum,stm);
}




}


