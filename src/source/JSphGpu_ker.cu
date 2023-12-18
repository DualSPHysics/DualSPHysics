//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2023 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpu_ker.cu \brief Implements functions and CUDA kernels for the Particle Interaction and System Update.

#include "JSphGpu_ker.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include "JLog2.h"
#include <cfloat>
#include <math_constants.h>
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"

#pragma warning(disable : 4267) //Cancels "warning C4267: conversion from 'size_t' to 'int', possible loss of data"
#pragma warning(disable : 4244) //Cancels "warning C4244: conversion from 'unsigned __int64' to 'unsigned int', possible loss of data"
#pragma warning(disable : 4503) //Cancels "warning C4503: decorated name length exceeded, name was truncated"
#include <thrust/device_vector.h>
#include <thrust/sort.h>

__constant__ StCteInteraction CTE;
#define CTE_AVAILABLE

namespace cusph{
#include "FunctionsBasic_iker.h"
#include "FunctionsMath_iker.h"
#include "FunctionsGeo3d_iker.h"
#include "FunSphKernel_iker.h"
#include "FunSphEos_iker.h"
#include "JCellSearch_iker.h"


//==============================================================================
/// Reduction using maximum of float values in shared memory for a warp.
/// Reduccion mediante maximo de valores float en memoria shared para un warp.
//==============================================================================
template <unsigned blockSize> __device__ void KerReduMaxFloatWarp(
  volatile float* sdat,unsigned tid)
{
  if(blockSize>=64)sdat[tid]=max(sdat[tid],sdat[tid+32]);
  if(blockSize>=32)sdat[tid]=max(sdat[tid],sdat[tid+16]);
  if(blockSize>=16)sdat[tid]=max(sdat[tid],sdat[tid+8]);
  if(blockSize>=8)sdat[tid]=max(sdat[tid],sdat[tid+4]);
  if(blockSize>=4)sdat[tid]=max(sdat[tid],sdat[tid+2]);
  if(blockSize>=2)sdat[tid]=max(sdat[tid],sdat[tid+1]);
}

//==============================================================================
/// Accumulates the maximum of n values of array dat[], storing the result in 
/// the beginning of res[].(Many positions of res[] are used as blocks, 
/// storing the final result in res[0]).
///
/// Acumula el maximo de n valores del vector dat[], guardando el resultado al 
/// principio de res[] (Se usan tantas posiciones del res[] como bloques, 
/// quedando el resultado final en res[0]).
//==============================================================================
template <unsigned blockSize> __global__ void KerReduMaxFloat(unsigned n
  ,unsigned ini,const float* dat,float* res)
{
  extern __shared__ float sdat[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  sdat[tid]=(c<n? dat[c+ini]: -FLT_MAX);
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sdat[tid]=max(sdat[tid],sdat[tid+256]);  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sdat[tid]=max(sdat[tid],sdat[tid+128]);  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sdat[tid]=max(sdat[tid],sdat[tid+64]);   __syncthreads(); }
  if(tid<32)KerReduMaxFloatWarp<blockSize>(sdat,tid);
  if(tid==0)res[blockIdx.x]=sdat[0];
}

//==============================================================================
/// Returns the maximum of an array, using resu[] as auxiliar array.
/// Size of resu[] must be >= a (N/SPHBSIZE+1)+(N/(SPHBSIZE*SPHBSIZE)+SPHBSIZE)
///
/// Devuelve el maximo de un vector, usando resu[] como vector auxiliar. El tamanho
/// de resu[] debe ser >= a (N/SPHBSIZE+1)+(N/(SPHBSIZE*SPHBSIZE)+SPHBSIZE)
//==============================================================================
float ReduMaxFloat(unsigned ndata,unsigned inidata,float* data,float* resu){
  float resf=0;
  if(ndata>=1){
    unsigned n=ndata,ini=inidata;
    unsigned smemSize=SPHBSIZE*sizeof(float);
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    unsigned n_blocks=sgrid.x*sgrid.y;
    float* dat=data;
    float* resu1=resu;
    float* resu2=resu+n_blocks;
    float* res=resu1;
    while(n>1){
      KerReduMaxFloat<SPHBSIZE><<<sgrid,SPHBSIZE,smemSize>>>(n,ini,dat,res);
      n=n_blocks; ini=0;
      sgrid=GetSimpleGridSize(n,SPHBSIZE);  
      n_blocks=sgrid.x*sgrid.y;
      if(n>1){
        dat=res; res=(dat==resu1? resu2: resu1); 
      }
    }
    if(ndata>1)cudaMemcpy(&resf,res,sizeof(float),cudaMemcpyDeviceToHost);
    else cudaMemcpy(&resf,data,sizeof(float),cudaMemcpyDeviceToHost);
  }
  //else{//-Using Thrust library is slower than ReduMasFloat() with ndata < 5M.
  //  thrust::device_ptr<float> dev_ptr(data);
  //  resf=thrust::reduce(dev_ptr,dev_ptr+ndata,-FLT_MAX,thrust::maximum<float>());
  //}
  return(resf);
}

//==============================================================================
/// Accumulates the sum of n values of array dat[], storing the result in 
/// the beginning of res[].(Many positions of res[] are used as blocks, 
/// storing the final result in res[0]).
///
/// Acumula la suma de n valores del vector dat[].w, guardando el resultado al 
/// principio de res[] (Se usan tantas posiciones del res[] como bloques, 
/// quedando el resultado final en res[0]).
//==============================================================================
template <unsigned blockSize> __global__ void KerReduMaxFloat_w(unsigned n
  ,unsigned ini,const float4* dat,float* res)
{
  extern __shared__ float sdat[];
  unsigned tid=threadIdx.x;
  unsigned c=blockIdx.x*blockDim.x + threadIdx.x;
  sdat[tid]=(c<n? dat[c+ini].w: -FLT_MAX);
  __syncthreads();
  if(blockSize>=512){ if(tid<256)sdat[tid]=max(sdat[tid],sdat[tid+256]);  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)sdat[tid]=max(sdat[tid],sdat[tid+128]);  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) sdat[tid]=max(sdat[tid],sdat[tid+64]);   __syncthreads(); }
  if(tid<32)KerReduMaxFloatWarp<blockSize>(sdat,tid);
  if(tid==0)res[blockIdx.x]=sdat[0];
}

//==============================================================================
/// Returns the maximum of an array, using resu[] as auxiliar array.
/// Size of resu[] must be >= a (N/SPHBSIZE+1)+(N/(SPHBSIZE*SPHBSIZE)+SPHBSIZE).
///
/// Devuelve el maximo de la componente w de un vector float4, usando resu[] como 
/// vector auxiliar. El tamanho de resu[] debe ser >= a (N/SPHBSIZE+1)+(N/(SPHBSIZE*SPHBSIZE)+SPHBSIZE).
//==============================================================================
float ReduMaxFloat_w(unsigned ndata,unsigned inidata,float4* data,float* resu){
  unsigned n=ndata,ini=inidata;
  unsigned smemSize=SPHBSIZE*sizeof(float);
  dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
  unsigned n_blocks=sgrid.x*sgrid.y;
  float* dat=NULL;
  float* resu1=resu;
  float* resu2=resu+n_blocks;
  float* res=resu1;
  while(n>1){
    if(!dat)KerReduMaxFloat_w<SPHBSIZE><<<sgrid,SPHBSIZE,smemSize>>>(n,ini,data,res);
    else KerReduMaxFloat<SPHBSIZE><<<sgrid,SPHBSIZE,smemSize>>>(n,ini,dat,res);
    n=n_blocks; ini=0;
    sgrid=GetSimpleGridSize(n,SPHBSIZE);  
    n_blocks=sgrid.x*sgrid.y;
    if(n>1){
      dat=res; res=(dat==resu1? resu2: resu1); 
    }
  }
  float resf;
  if(ndata>1)cudaMemcpy(&resf,res,sizeof(float),cudaMemcpyDeviceToHost);
  else{
    float4 resf4;
    cudaMemcpy(&resf4,data,sizeof(float4),cudaMemcpyDeviceToHost);
    resf=resf4.w;
  }
  return(resf);
}

//==============================================================================
/// Stores constants for the GPU interaction.
/// Graba constantes para la interaccion a la GPU.
//==============================================================================
void CteInteractionUp(const StCteInteraction* cte){
  cudaMemcpyToSymbol(CTE,cte,sizeof(StCteInteraction));
}

//------------------------------------------------------------------------------
/// Initialises array with the indicated value.
/// Inicializa array con el valor indicado.
//------------------------------------------------------------------------------
__global__ void KerInitArray(unsigned n,float3* v,float3 value)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n)v[p]=value;
}

//==============================================================================
/// Initialises array with the indicated value.
/// Inicializa array con el valor indicado.
//==============================================================================
void InitArray(unsigned n,float3* v,tfloat3 value){
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerInitArray <<<sgrid,SPHBSIZE>>> (n,v,Float3(value));
  }
}

//------------------------------------------------------------------------------
/// Sets v[].y to zero.
/// Pone v[].y a cero.
//------------------------------------------------------------------------------
__global__ void KerResety(unsigned n,unsigned ini,float3* v)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n)v[p+ini].y=0;
}

//==============================================================================
/// Sets v[].y to zero.
/// Pone v[].y a cero.
//==============================================================================
void Resety(unsigned n,unsigned ini,float3* v){
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerResety <<<sgrid,SPHBSIZE>>> (n,ini,v);
  }
}

//------------------------------------------------------------------------------
/// Calculates module^2 of ace.
//------------------------------------------------------------------------------
__global__ void KerComputeAceMod(unsigned n,const float3* ace,float* acemod)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const float3 r=ace[p];
    acemod[p]=r.x*r.x+r.y*r.y+r.z*r.z;
  }
}

//==============================================================================
/// Calculates module^2 of ace.
//==============================================================================
void ComputeAceMod(unsigned n,const float3* ace,float* acemod){
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeAceMod <<<sgrid,SPHBSIZE>>> (n,ace,acemod);
  }
}

//------------------------------------------------------------------------------
/// Calculates module^2 of ace, comprobando que la particula sea normal.
/// Uses zero for periodic particles.
//------------------------------------------------------------------------------
__global__ void KerComputeAceMod(unsigned n,const typecode* code
  ,const float3* ace,float* acemod)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const typecode rcod=code[p];
    const float3 r=(CODE_IsNormal(rcod) && !CODE_IsFluidInout(rcod)? ace[p]: make_float3(0,0,0));
    acemod[p]=r.x*r.x+r.y*r.y+r.z*r.z;
  }
}

//==============================================================================
/// Calculates module^2 of ace, comprobando que la particula sea normal.
/// Uses zero for periodic particles.
//==============================================================================
void ComputeAceMod(unsigned n,const typecode* code,const float3* ace
  ,float* acemod)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeAceMod <<<sgrid,SPHBSIZE>>> (n,code,ace,acemod);
  }
}


//##############################################################################
//# Other kernels...
//# Otros kernels...
//##############################################################################
//------------------------------------------------------------------------------
/// Calculates module^2 of vel.
//------------------------------------------------------------------------------
__global__ void KerComputeVelMod(unsigned n,const float4* vel,float* velmod)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const float4 r=vel[p];
    velmod[p]=r.x*r.x+r.y*r.y+r.z*r.z;
  }
}

//==============================================================================
/// Calculates module^2 of vel.
//==============================================================================
void ComputeVelMod(unsigned n,const float4* vel,float* velmod){
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeVelMod <<<sgrid,SPHBSIZE>>> (n,vel,velmod);
  }
}


//##############################################################################
//# Kernels para cambiar la posicion.
//# Kernels for changing the position.
//##############################################################################
//------------------------------------------------------------------------------
/// Updates pos, dcell and code from the indicated displacement.
/// The code may be CODE_OUTRHO because in ComputeStepVerlet / Symplectic this is evaluated
/// and is executed before ComputeStepPos.
/// Checks limits depending on maprealposmin and maprealsize, this is valid 
/// for single-GPU because maprealpos and domrealpos are equal. For multi-gpu it is
/// important to mark particles that leave the domain without leaving the map.
///
/// Actualiza pos, dcell y code a partir del desplazamiento indicado.
/// Code puede ser CODE_OUTRHO pq en ComputeStepVerlet/Symplectic se evalua esto 
/// y se ejecuta antes que ComputeStepPos.
/// Comprueba los limites en funcion de maprealposmin y maprealsize esto es valido
/// para single-gpu pq domrealpos y maprealpos son iguales. Para multi-gpu seria 
/// necesario marcar las particulas q salgan del dominio sin salir del mapa.
//------------------------------------------------------------------------------
template<bool periactive> __device__ void KerUpdatePos(
  double2 rxy,double rz,double movx,double movy,double movz
  ,bool outrho,unsigned p,double2* posxy,double* posz,unsigned* dcell
  ,typecode* code)
{
  //-Checks validity of displacement. | Comprueba validez del desplazamiento.
  const bool outmov=(fmaxf(fabsf(float(movx)),fmaxf(fabsf(float(movy)),fabsf(float(movz))))>CTE.movlimit);
  //-Applies diplacement.
  double3 rpos=make_double3(rxy.x,rxy.y,rz);
  rpos.x+=movx; rpos.y+=movy; rpos.z+=movz;
  if(rpos.y<0 && CTE.symmetry)rpos.y=-rpos.y; //<vs_syymmetry>
  //-Checks limits of real domain. | Comprueba limites del dominio reales.
  double dx=rpos.x-CTE.maprealposminx;
  double dy=rpos.y-CTE.maprealposminy;
  double dz=rpos.z-CTE.maprealposminz;
  bool out=(dx!=dx || dy!=dy || dz!=dz || dx<0 || dy<0 || dz<0 || dx>=CTE.maprealsizex || dy>=CTE.maprealsizey || dz>=CTE.maprealsizez);
  if(periactive && out){
    bool xperi=(CTE.periactive&1),yperi=(CTE.periactive&2),zperi=(CTE.periactive&4);
    if(xperi){
      if(dx<0)                { dx-=CTE.xperincx; dy-=CTE.xperincy; dz-=CTE.xperincz; }
      if(dx>=CTE.maprealsizex){ dx+=CTE.xperincx; dy+=CTE.xperincy; dz+=CTE.xperincz; }
    }
    if(yperi){
      if(dy<0)                { dx-=CTE.yperincx; dy-=CTE.yperincy; dz-=CTE.yperincz; }
      if(dy>=CTE.maprealsizey){ dx+=CTE.yperincx; dy+=CTE.yperincy; dz+=CTE.yperincz; }
    }
    if(zperi){
      if(dz<0)                { dx-=CTE.zperincx; dy-=CTE.zperincy; dz-=CTE.zperincz; }
      if(dz>=CTE.maprealsizez){ dx+=CTE.zperincx; dy+=CTE.zperincy; dz+=CTE.zperincz; }
    }
    bool outx=!xperi && (dx<0 || dx>=CTE.maprealsizex);
    bool outy=!yperi && (dy<0 || dy>=CTE.maprealsizey);
    bool outz=!zperi && (dz<0 || dz>=CTE.maprealsizez);
    out=(outx||outy||outz);
    rpos=make_double3(dx+CTE.maprealposminx,dy+CTE.maprealposminy,dz+CTE.maprealposminz);
  }
  //-Stores updated position.
  posxy[p]=make_double2(rpos.x,rpos.y);
  posz[p]=rpos.z;
  //-Stores cell and check. | Guarda celda y check.
  if(outrho || outmov || out){//-Particle out. Only brands as excluded normal particles (not periodic). | Particle out. Solo las particulas normales (no periodicas) se pueden marcar como excluidas.
    typecode rcode=code[p];
    if(out)rcode=CODE_SetOutPos(rcode);
    else if(outrho)rcode=CODE_SetOutRho(rcode);
    else rcode=CODE_SetOutMov(rcode);
    code[p]=rcode;
    dcell[p]=DCEL_CodeMapOut;
  }
  else{//-Particle in.
    if(periactive){
      dx=rpos.x-CTE.domposminx;
      dy=rpos.y-CTE.domposminy;
      dz=rpos.z-CTE.domposminz;
    }
    const unsigned cx=unsigned(dx/CTE.scell);
    const unsigned cy=unsigned(dy/CTE.scell);
    const unsigned cz=unsigned(dz/CTE.scell);
    dcell[p]=DCEL_Cell(CTE.cellcode,cx,cy,cz);
  }
}

//------------------------------------------------------------------------------
/// Returns the corrected position after applying periodic conditions.
/// Devuelve la posicion corregida tras aplicar condiciones periodicas.
//------------------------------------------------------------------------------
__device__ double3 KerUpdatePeriodicPos(double3 ps)
{
  double dx=ps.x-CTE.maprealposminx;
  double dy=ps.y-CTE.maprealposminy;
  double dz=ps.z-CTE.maprealposminz;
  const bool out=(dx!=dx || dy!=dy || dz!=dz || dx<0 || dy<0 || dz<0 || dx>=CTE.maprealsizex || dy>=CTE.maprealsizey || dz>=CTE.maprealsizez);
  //-Adjusts position according to periodic conditions and rechecks domain limits.
  //-Ajusta posicion segun condiciones periodicas y vuelve a comprobar los limites del dominio.
  if(out){
    bool xperi=(CTE.periactive&1),yperi=(CTE.periactive&2),zperi=(CTE.periactive&4);
    if(xperi){
      if(dx<0)                { dx-=CTE.xperincx; dy-=CTE.xperincy; dz-=CTE.xperincz; }
      if(dx>=CTE.maprealsizex){ dx+=CTE.xperincx; dy+=CTE.xperincy; dz+=CTE.xperincz; }
    }
    if(yperi){
      if(dy<0)                { dx-=CTE.yperincx; dy-=CTE.yperincy; dz-=CTE.yperincz; }
      if(dy>=CTE.maprealsizey){ dx+=CTE.yperincx; dy+=CTE.yperincy; dz+=CTE.yperincz; }
    }
    if(zperi){
      if(dz<0)                { dx-=CTE.zperincx; dy-=CTE.zperincy; dz-=CTE.zperincz; }
      if(dz>=CTE.maprealsizez){ dx+=CTE.zperincx; dy+=CTE.zperincy; dz+=CTE.zperincz; }
    }
    ps=make_double3(dx+CTE.maprealposminx,dy+CTE.maprealposminy,dz+CTE.maprealposminz);
  }
  return(ps);
}


//##############################################################################
//# Kernels for calculating forces (Pos-Double).
//# Kernels para calculo de fuerzas (Pos-Double).
//##############################################################################
//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles. Bound-Fluid/Float
/// Realiza la interaccion de una particula con un conjunto de ellas. Bound-Fluid/Float
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,bool symm>
  __device__ void KerInteractionForcesBoundBox
  (unsigned p1,const unsigned& pini,const unsigned& pfin
  ,const float* ftomassp
  ,const float4* poscell,const float4* velrho,const typecode* code,const unsigned* idp
  ,float massf,const float4& pscellp1,const float4& velrhop1,float& arp1,float& visc)
{
  for(int p2=pini;p2<pfin;p2++){
    const float4 pscellp2=poscell[p2];
    float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
    float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
    float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
    if(symm)dry=pscellp1.y+pscellp2.y + CTE.poscellsize*PSCEL_GetfY(pscellp2.w); //<vs_syymmetry>
    const float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

      float4 velrhop2=velrho[p2];
      if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>
      //-Obtains particle mass p2 if there are floating bodies.
      //-Obtiene masa de particula p2 en caso de existir floatings.
      float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massf si es fluid.
      bool compute=true; //-Deactivated when DEM is used and is float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
      if(USE_FLOATING){
        const typecode cod=code[p2];
        bool ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massf);
        compute=!(USE_FTEXTERNAL && ftp2); //-Deactivated when DEM or Chrono is used and is bound-float. | Se desactiva cuando se usa DEM o Chrono y es bound-float.
      }

      if(compute){
        //-Density derivative (Continuity equation).
        const float dvx=velrhop1.x-velrhop2.x, dvy=velrhop1.y-velrhop2.y, dvz=velrhop1.z-velrhop2.z;
        arp1+=(USE_FLOATING? ftmassp2: massf)*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

        {//===== Viscosity ===== 
          const float dot=drx*dvx + dry*dvy + drz*dvz;
          const float dot_rr2=dot/(rr2+CTE.eta2);
          visc=max(dot_rr2,visc); 
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Particle interaction. Bound-Fluid/Float
/// Realiza interaccion entre particulas. Bound-Fluid/Float
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,bool symm> 
  __global__ void KerInteractionForcesBound(unsigned n,unsigned pinit
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid,const unsigned* dcell
  ,const float* ftomassp
  ,const float4* poscell,const float4* velrho,const typecode* code,const unsigned* idp
  ,float* viscdt,float* ar)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of thread.
  if(p<n){
    const unsigned p1=p+pinit;      //-Number of particle.
    float visc=0,arp1=0;

    //-Loads particle p1 data.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrho[p1];
    const bool rsymp1=(symm && PSCEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>
    
    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Boundary-Fluid interaction.
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
      if(pfin){
                          KerInteractionForcesBoundBox<tker,ftmode,false> (p1,pini,pfin,ftomassp,poscell,velrho,code,idp,CTE.massf,pscellp1,velrhop1,arp1,visc);
        if(symm && rsymp1)KerInteractionForcesBoundBox<tker,ftmode,true > (p1,pini,pfin,ftomassp,poscell,velrho,code,idp,CTE.massf,pscellp1,velrhop1,arp1,visc);
      }
    }
    //-Stores results.
    if(arp1 || visc){
      ar[p1]+=arp1;
      if(visc>viscdt[p1])viscdt[p1]=visc;
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles. (Fluid/Float-Fluid/Float/Bound)
/// Realiza la interaccion de una particula con un conjunto de ellas. (Fluid/Float-Fluid/Float/Bound)
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool symm>
  __device__ void KerInteractionForcesFluidBox(bool boundp2,unsigned p1
  ,const unsigned& pini,const unsigned& pfin,float visco
  ,const float* ftomassp,const float2* tauff,const float3* dengradcorr
  ,const float4* poscell,const float4* velrho,const typecode* code,const unsigned* idp
  ,float massp2,bool ftp1
  ,const float4& pscellp1,const float4& velrhop1,float pressp1
  ,const float2& taup1_xx_xy,const float2& taup1_xz_yy,const float2& taup1_yz_zz
  ,float2& two_strainp1_xx_xy,float2& two_strainp1_xz_yy,float2& two_strainp1_yz_zz
  ,float3& acep1,float& arp1,float& visc,float& deltap1
  ,TpShifting shiftmode,float4& shiftposfsp1
  ,const float3* boundnormal, const float* boundonoff, const float3* motionvel)
{
  for(int p2=pini;p2<pfin;p2++){
    const float4 pscellp2=poscell[p2];
    float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
    float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
    float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
    if(symm)dry=pscellp1.y+pscellp2.y + CTE.poscellsize*PSCEL_GetfY(pscellp2.w); //<vs_syymmetry>
    const float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

      //-Obtains mass of particle p2 if any floating bodies exist.
      //-Obtiene masa de particula p2 en caso de existir floatings.
      bool ftp2=false;         //-Indicates if it is floating. | Indica si es floating.
      float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
      bool compute=true; //-Deactivated when DEM is used and is float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
      if(USE_FLOATING){
        const typecode cod=code[p2];
        ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);
        #ifdef DELTA_HEAVYFLOATING
          if(ftp2 && tdensity==DDT_DDT && ftmassp2<=(massp2*1.2f))deltap1=FLT_MAX;
        #else
          if(ftp2 && tdensity==DDT_DDT)deltap1=FLT_MAX;
        #endif
        if(ftp2 && shift && shiftmode==SHIFT_NoBound)shiftposfsp1.x=FLT_MAX; //-Cancels shifting with floating bodies. | Con floatings anula shifting.
        compute=!(USE_FTEXTERNAL && ftp1 && (boundp2 || ftp2)); //-Deactivated when DEM or Chrono is used and is float-float or float-bound. | Se desactiva cuando se usa DEM o Chrono y es float-float o float-bound.
      }

      float4 velrhop2=velrho[p2];
      if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

      // SHABA setting up boundary normals and new boundary velocities
      float3 normalp2 = make_float3(0, 0, 0); // creating a normailsed boundary normal
      float3 normalvelp2 = make_float3(0, 0, 0); // boundary particle velocity normal to boundary used in con of mass
      float3 tangentvelp2 = make_float3(0, 0, 0); // boundary particle velocity tangent to boundary used in momentu
      float3 movvelp2 = motionvel[p2];
      
      if (boundp2) {//changing things if fluid-bound
          float norm = sqrt(boundnormal[p2].x * boundnormal[p2].x + boundnormal[p2].y * boundnormal[p2].y + boundnormal[p2].z * boundnormal[p2].z);
          if (norm > 0.f) { // mDBC
              normalp2.x = boundnormal[p2].x / norm; normalp2.y = boundnormal[p2].y / norm; normalp2.z = boundnormal[p2].z / norm;
              // calculating boundary particle velocity components
              float veldotnorm = velrhop2.x * normalp2.x + velrhop2.y * normalp2.y + velrhop2.z * normalp2.z;
              normalvelp2.x = veldotnorm * normalp2.x; normalvelp2.y = veldotnorm * normalp2.y; normalvelp2.z = veldotnorm * normalp2.z;
              tangentvelp2.x = velrhop2.x - normalvelp2.x; tangentvelp2.y = velrhop2.y - normalvelp2.y; tangentvelp2.z = velrhop2.z - normalvelp2.z;
          }
          else { // DBC hopefulyl this doesnt break it
              normalvelp2.x = velrhop2.x; normalvelp2.y = velrhop2.y; normalvelp2.z = velrhop2.z;
              tangentvelp2.x = velrhop2.x; tangentvelp2.y = velrhop2.y; tangentvelp2.z = velrhop2.z;
          }
          // changing the mass of boundary particle with boundonoff
          massp2 = boundonoff[p2] * massp2; 
      }
      
      //-Velocity derivative (Momentum equation).
      if(compute){
        const float pressp2=cufsph::ComputePressCte(velrhop2.w);
        const float prs=(pressp1+pressp2)/(velrhop1.w*velrhop2.w)
          +(tker==KERNEL_Cubic? cufsph::GetKernelCubic_Tensil(rr2,velrhop1.w,pressp1,velrhop2.w,pressp2): 0);
        const float p_vpm=-prs*(USE_FLOATING? ftmassp2: massp2);
        acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
      }

      //-Density derivative (Continuity equation).
      float dvx=velrhop1.x-velrhop2.x, dvy=velrhop1.y-velrhop2.y, dvz=velrhop1.z-velrhop2.z;
      if (boundp2) {
          dvx = velrhop1.x - movvelp2.x; dvy = velrhop1.y - movvelp2.y; dvz = velrhop1.z - movvelp2.z; // SHABA no slip mDBC
      }
      if(compute)arp1+=(USE_FLOATING? ftmassp2: massp2)*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

      const float cbar=CTE.cs0;
      const float dot3=(tdensity!=DDT_None || shift? drx*frx+dry*fry+drz*frz: 0);
      //-Density Diffusion Term (Molteni and Colagrossi 2009).
      if(tdensity==DDT_DDT && deltap1!=FLT_MAX){
        const float rhop1over2=velrhop1.w/velrhop2.w;
        const float visc_densi=CTE.ddtkh*cbar*(rhop1over2-1.f)/(rr2+CTE.eta2);
        const float delta=visc_densi*dot3*(USE_FLOATING? ftmassp2: massp2);
        //deltap1=(boundp2? FLT_MAX: deltap1+delta);
        deltap1=(boundp2 && CTE.tboundary==BC_DBC? FLT_MAX: deltap1+delta);
      }
      //-Density Diffusion Term (Fourtakas et al 2019).
      if((tdensity==DDT_DDT2 || (tdensity==DDT_DDT2Full && !boundp2)) && deltap1!=FLT_MAX && !ftp2){
        const float rh=1.f+CTE.ddtgz*drz;
        const float drho=CTE.rhopzero*pow(rh,1.f/CTE.gamma)-CTE.rhopzero;  
        const float visc_densi=CTE.ddtkh*cbar*((velrhop2.w-velrhop1.w)-drho)/(rr2+CTE.eta2);
        const float delta=visc_densi*dot3*massp2/velrhop2.w;
        deltap1=(boundp2? FLT_MAX: deltap1-delta); //-blocks it makes it boil - bloody DBC
      }

      //-Shifting correction.
      if(shift && shiftposfsp1.x!=FLT_MAX){
        const float massrho=(USE_FLOATING? ftmassp2: massp2)/velrhop2.w;
        const bool noshift=(boundp2 && (shiftmode==SHIFT_NoBound || (shiftmode==SHIFT_NoFixed && CODE_IsFixed(code[p2]))));
        shiftposfsp1.x=(noshift? FLT_MAX: shiftposfsp1.x+massrho*frx); //-Removes shifting for the boundaries. | Con boundary anula shifting.
        shiftposfsp1.y+=massrho*fry;
        shiftposfsp1.z+=massrho*frz;
        shiftposfsp1.w-=massrho*dot3;
      }

      //===== Viscosity ===== 
      if(compute){
        if (boundp2) {
            dvx = velrhop1.x - tangentvelp2.x; dvy = velrhop1.y - tangentvelp2.y; dvz = velrhop1.z - tangentvelp2.z; // SHABA no slip mDBC
        }
        const float dot=drx*dvx + dry*dvy + drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        visc=max(dot_rr2,visc);  //ViscDt=max(dot/(rr2+Eta2),ViscDt);
        if(tvisco==VISCO_Artificial){//-Artificial viscosity.
          if(dot<0){
            const float amubar=CTE.kernelh*dot_rr2;  //amubar=CTE.kernelh*dot/(rr2+CTE.eta2);
            const float robar=(velrhop1.w+velrhop2.w)*0.5f;
            const float pi_visc=(-visco*cbar*amubar/robar)*(USE_FLOATING? ftmassp2: massp2);
            acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
          }
        }
        if(tvisco==VISCO_Laminar || tvisco==VISCO_LaminarSPS){//-Laminar and Laminar+SPS viscosity.
          const float robar2=(velrhop1.w+velrhop2.w);
          const float temp=4.f*visco/((rr2+CTE.eta2)*robar2);  //-Simplication of temp=2.0f*visco/((rr2+CTE.eta2)*robar); robar=(rhopp1+velrhop2.w)*0.5f;
          const float vtemp=(USE_FLOATING? ftmassp2: massp2)*temp*(drx*frx+dry*fry+drz*frz);  
          acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
        }
        if(tvisco==VISCO_LaminarSPS){//-SPS contribution for Laminar viscosity. 
          //-SPS turbulence model.
          //-Note that taup1 is tau_a/rho_a^2 for interaction of particle a and b.
          //-And taup1 is always zero when p1 is not a fluid particle.
          float2 stau_xx_xy=taup1_xx_xy;
          float2 stau_xz_yy=taup1_xz_yy;
          float2 stau_yz_zz=taup1_yz_zz;
          if(!boundp2 && (USE_NOFLOATING || !ftp2)){//-When p2 is fluid.
            //-Note that taup2 is tau_b/rho_b^2 for interaction of particle a and b.
            float2 taup2=tauff[p2*3];     stau_xx_xy.x+=taup2.x; stau_xx_xy.y+=taup2.y;
                   taup2=tauff[p2*3+1];   stau_xz_yy.x+=taup2.x; stau_xz_yy.y+=taup2.y;
                   taup2=tauff[p2*3+2];   stau_yz_zz.x+=taup2.x; stau_yz_zz.y+=taup2.y;
          }
          acep1.x+=(USE_FLOATING? ftmassp2: massp2)*(stau_xx_xy.x*frx + stau_xx_xy.y*fry + stau_xz_yy.x*frz);
          acep1.y+=(USE_FLOATING? ftmassp2: massp2)*(stau_xx_xy.y*frx + stau_xz_yy.y*fry + stau_yz_zz.x*frz);
          acep1.z+=(USE_FLOATING? ftmassp2: massp2)*(stau_xz_yy.x*frx + stau_yz_zz.x*fry + stau_yz_zz.y*frz);
          //-Velocity gradients.
          if(USE_NOFLOATING || !ftp1){//-When p1 is fluid.
            const float volp2=-(USE_FLOATING? ftmassp2: massp2)/velrhop2.w;
            float dv=dvx*volp2; two_strainp1_xx_xy.x+=dv*frx; two_strainp1_xx_xy.y+=dv*fry; two_strainp1_xz_yy.x+=dv*frz;
                  dv=dvy*volp2; two_strainp1_xx_xy.y+=dv*frx; two_strainp1_xz_yy.y+=dv*fry; two_strainp1_yz_zz.x+=dv*frz;
                  dv=dvz*volp2; two_strainp1_xz_yy.x+=dv*frx; two_strainp1_yz_zz.x+=dv*fry; two_strainp1_yz_zz.y+=dv*frz;
            // to compute tau terms we assume that two_strain.xy=dudy+dvdx, two_strain.xz=dudz+dwdx, two_strain.yz=dvdz+dwdy
            // so only 6 elements are needed instead of 3x3.
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction between particles. Fluid/Float-Fluid/Float or Fluid/Float-Bound.
/// Includes artificial/laminar viscosity and normal/DEM floating bodies.
///
/// Realiza interaccion entre particulas. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings normales/dem.
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool symm>
  __global__ void KerInteractionForcesFluid(unsigned n,unsigned pinit,float viscob,float viscof
  ,int scelldiv,int4 nc,int3 cellzero,const int2* begincell,unsigned cellfluid
  ,const unsigned* dcell,const float* ftomassp,const float2* tauff,float2* two_strain
  ,const float3* dengradcorr,const float4* poscell,const float4* velrho
  ,const typecode* code,const unsigned* idp,float* viscdt,float* ar,float3* ace
  ,float* delta,TpShifting shiftmode,float4* shiftposfs
  ,const float3* boundnormal,const float* boundonoff,const float3* motionvel)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pinit;      //-Number of particle.
    float visc=0,arp1=0,deltap1=0;
    float3 acep1=make_float3(0,0,0);

    //-Variables for Shifting.
    float4 shiftposfsp1;
    if(shift)shiftposfsp1=shiftposfs[p1];

    //-Obtains data of particle p1 in case there are floating bodies.
    bool ftp1;       //-Indicates if it is floating. | Indica si es floating.
    if(USE_FLOATING){
      const typecode cod=code[p1];
      ftp1=CODE_IsFloating(cod);
      if(ftp1 && tdensity!=DDT_None)deltap1=FLT_MAX; //-DDT is not applied to floating particles.
      if(ftp1 && shift)shiftposfsp1.x=FLT_MAX; //-Shifting is not calculated for floating bodies. | Para floatings no se calcula shifting.
    }

    //-Obtains basic data of particle p1.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrho[p1];
    const float pressp1=cufsph::ComputePressCte(velrhop1.w);
    const bool rsymp1=(symm && PSCEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>

    //-Variables for Laminar+SPS.
    float2 taup1_xx_xy,taup1_xz_yy,taup1_yz_zz; //-Note that taup1 is tau_a/rho_a^2.
    if(tvisco==VISCO_LaminarSPS){
      taup1_xx_xy=tauff[p1*3];
      taup1_xz_yy=tauff[p1*3+1];
      taup1_yz_zz=tauff[p1*3+2];
    }
    //-Variables for Laminar+SPS (computation).
    float2 two_strainp1_xx_xy,two_strainp1_xz_yy,two_strainp1_yz_zz;
    if(tvisco==VISCO_LaminarSPS){
      two_strainp1_xx_xy=make_float2(0,0);
      two_strainp1_xz_yy=make_float2(0,0);
      two_strainp1_yz_zz=make_float2(0,0);
    }

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
                          KerInteractionForcesFluidBox<tker,ftmode,tvisco,tdensity,shift,false> (false
                            ,p1,pini,pfin,viscof,ftomassp,tauff,dengradcorr,poscell,velrho,code,idp
                            ,CTE.massf,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz
                            ,two_strainp1_xx_xy,two_strainp1_xz_yy,two_strainp1_yz_zz,acep1,arp1,visc
                            ,deltap1,shiftmode,shiftposfsp1,boundnormal,boundonoff,motionvel);
        //<vs_syymmetry_ini>
        if(symm && rsymp1)KerInteractionForcesFluidBox<tker,ftmode,tvisco,tdensity,shift,true > (false
                            ,p1,pini,pfin,viscof,ftomassp,tauff,dengradcorr,poscell,velrho,code,idp
                            ,CTE.massf,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz
                            ,two_strainp1_xx_xy,two_strainp1_xz_yy,two_strainp1_yz_zz,acep1,arp1,visc
                            ,deltap1,shiftmode,shiftposfsp1,boundnormal,boundonoff,motionvel);
        //<vs_syymmetry_end>
      }
    }
    //-Interaction with boundaries.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
                          KerInteractionForcesFluidBox<tker,ftmode,tvisco,tdensity,shift,false> (true
                            ,p1,pini,pfin,viscob,ftomassp,tauff,NULL,poscell,velrho,code,idp,CTE.massb
                            ,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz
                            ,two_strainp1_xx_xy,two_strainp1_xz_yy,two_strainp1_yz_zz,acep1,arp1,visc
                            ,deltap1,shiftmode,shiftposfsp1,boundnormal,boundonoff,motionvel);
        //<vs_syymmetry_ini>
        if(symm && rsymp1)KerInteractionForcesFluidBox<tker,ftmode,tvisco,tdensity,shift,true > (true 
                            ,p1,pini,pfin,viscob,ftomassp,tauff,NULL,poscell,velrho,code,idp,CTE.massb
                            ,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz
                            ,two_strainp1_xx_xy,two_strainp1_xz_yy,two_strainp1_yz_zz,acep1,arp1,visc
                            ,deltap1,shiftmode,shiftposfsp1,boundnormal,boundonoff,motionvel);
        //<vs_syymmetry_end>
      }
    }

    //-Stores results.
    if(shift||arp1||acep1.x||acep1.y||acep1.z||visc){
      if(tdensity!=DDT_None){
        if(delta){
          const float rdelta=delta[p1];
          delta[p1]=(rdelta==FLT_MAX || deltap1==FLT_MAX? FLT_MAX: rdelta+deltap1);
        }
        else if(deltap1!=FLT_MAX)arp1+=deltap1;
      }
      ar[p1]+=arp1;
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      if(visc>viscdt[p1])viscdt[p1]=visc;
      if(tvisco==VISCO_LaminarSPS){
        float2 rg;
        rg=two_strain[p1*3  ];  rg=make_float2(rg.x+two_strainp1_xx_xy.x,rg.y+two_strainp1_xx_xy.y);  two_strain[p1*3  ]=rg;
        rg=two_strain[p1*3+1];  rg=make_float2(rg.x+two_strainp1_xz_yy.x,rg.y+two_strainp1_xz_yy.y);  two_strain[p1*3+1]=rg;
        rg=two_strain[p1*3+2];  rg=make_float2(rg.x+two_strainp1_yz_zz.x,rg.y+two_strainp1_yz_zz.y);  two_strain[p1*3+2]=rg;
      }
      if(shift)shiftposfs[p1]=shiftposfsp1;
    }
  }
}

#ifndef DISABLE_BSMODES
//==============================================================================
/// Collects kernel information.
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool symm> 
  void Interaction_ForcesT_KerInfo(StKerInfo* kerinfo)
{
 #if CUDART_VERSION >= 6050
  {
    typedef void (*fun_ptr)(unsigned,unsigned,float,float,int,int4,int3,const int2*,unsigned,const unsigned*,const float*,const float2*,float2*,const float3*,const float4*,const float4*,const typecode*,const unsigned*,float*,float*,float3*,float*,TpShifting,float4*
      ,const float3* boundnormal,const float* boundonoff,const float3* motionvel);
    fun_ptr ptr=&KerInteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,symm>;
    int qblocksize=0,mingridsize=0;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize,&qblocksize,(void*)ptr,0,0);
    struct cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr,(void*)ptr);
    kerinfo->forcesfluid_bs=qblocksize;
    kerinfo->forcesfluid_rg=attr.numRegs;
    kerinfo->forcesfluid_bsmax=attr.maxThreadsPerBlock;
    //printf(">> KerInteractionForcesFluid  blocksize:%u (%u)\n",qblocksize,0);
  }
  {
    typedef void (*fun_ptr)(unsigned,unsigned,int,int4,int3,const int2*,const unsigned*,const float*,const float4*,const float4*,const typecode*,const unsigned*,float*,float*);
    fun_ptr ptr=&KerInteractionForcesBound<tker,ftmode,symm>;
    int qblocksize=0,mingridsize=0;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize,&qblocksize,(void*)ptr,0,0);
    struct cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr,(void*)ptr);
    kerinfo->forcesbound_bs=qblocksize;
    kerinfo->forcesbound_rg=attr.numRegs;
    kerinfo->forcesbound_bsmax=attr.maxThreadsPerBlock;
    //printf(">> KerInteractionForcesBound  blocksize:%u (%u)\n",qblocksize,0);
  }
  fcuda::Check_CudaErroorFun("Error collecting kernel information.");
 #endif
}
#endif

//==============================================================================
/// Interaction for the force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift> 
  void Interaction_ForcesGpuT(const StInterParmsg& t)
{
  //-Collects kernel information.
#ifndef DISABLE_BSMODES
  if(t.kerinfo){
    Interaction_ForcesT_KerInfo<tker,ftmode,tvisco,tdensity,shift,false>(t.kerinfo);
    return;
  }
#endif
  const StDivDataGpu& dvd=t.divdatag;
  const int2* beginendcell=dvd.beginendcell;
  //-Interaction Fluid-Fluid & Fluid-Bound.
  if(t.fluidnum){
    //printf("[ns:%u  id:%d] halo:%d fini:%d(%d) bini:%d(%d)\n",t.nstep,t.id,t.halo,t.fluidini,t.fluidnum,t.boundini,t.boundnum);
    dim3 sgridf=GetSimpleGridSize(t.fluidnum,t.bsfluid);
    if(t.symmetry) //<vs_syymmetry_ini>
      KerInteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,true> <<<sgridf,t.bsfluid,0,t.stm>>> 
      (t.fluidnum,t.fluidini,t.viscob,t.viscof,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
      ,t.ftomassp,(const float2*)t.spstaurho2,(float2*)t.sps2strain,t.dengradcorr,t.poscell,t.velrho,t.code,t.idp
      ,t.viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs,t.boundnormal,t.boundonoff,t.motionvel);
    else //<vs_syymmetry_end>
      KerInteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,false> <<<sgridf,t.bsfluid,0,t.stm>>> 
      (t.fluidnum,t.fluidini,t.viscob,t.viscof,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
      ,t.ftomassp,(const float2*)t.spstaurho2,(float2*)t.sps2strain,t.dengradcorr,t.poscell,t.velrho,t.code,t.idp
      ,t.viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs,t.boundnormal,t.boundonoff,t.motionvel);
      //KerInteractionForcesFluid<tker,ftmode,tvisco,tdensity,shift,false> <<<sgridf,t.bsfluid,0,t.stm>>> (t.fluidnum,t.fluidini,t.scelldiv,t.nc,t.cellfluid,t.viscob,t.viscof,t.begincell,Int3(t.cellmin),t.dcell,t.ftomassp,(const float2*)t.tau,(float2*)t.gradvel,t.poscell,t.velrho,t.code,t.idp,t.viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);
  }
  //-Interaction Boundary-Fluid.
  if(t.boundnum){
    const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
    dim3 sgridb=GetSimpleGridSize(t.boundnum,t.bsbound);
    //printf("bsbound:%u\n",bsbound);
    if(t.symmetry) //<vs_syymmetry_ini>
      KerInteractionForcesBound<tker,ftmode,true > <<<sgridb,t.bsbound,0,t.stm>>> 
      (t.boundnum,t.boundini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcell+dvd.cellfluid,t.dcell
        ,t.ftomassp,t.poscell,t.velrho,t.code,t.idp,t.viscdt,t.ar);
    else //<vs_syymmetry_end>
      KerInteractionForcesBound<tker,ftmode,false> <<<sgridb,t.bsbound,0,t.stm>>> 
      (t.boundnum,t.boundini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.dcell
        ,t.ftomassp,t.poscell,t.velrho,t.code,t.idp,t.viscdt,t.ar);
  }
}

//==============================================================================
//#define FAST_COMPILATION
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco>
  void Interaction_Forces_gt2(const StInterParmsg& t)
{
#ifdef FAST_COMPILATION
  if(t.shiftmode || t.tdensity!=DDT_DDT4)throw "Shifting and extra DDT are disabled for FastCompilation...";
  Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_DDT4,false> (t);
#else
  if(t.shiftmode){               const bool shift=true;
    if(t.tdensity==DDT_None)    Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_None    ,shift> (t);
    if(t.tdensity==DDT_DDT)     Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_DDT     ,shift> (t);
    if(t.tdensity==DDT_DDT2)    Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_DDT2    ,shift> (t);
    if(t.tdensity==DDT_DDT2Full)Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_DDT2Full,shift> (t);
  }
  else{                           const bool shift=false;
    if(t.tdensity==DDT_None)    Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_None    ,shift> (t);
    if(t.tdensity==DDT_DDT)     Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_DDT     ,shift> (t);
    if(t.tdensity==DDT_DDT2)    Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_DDT2    ,shift> (t);
    if(t.tdensity==DDT_DDT2Full)Interaction_ForcesGpuT<tker,ftmode,tvisco,DDT_DDT2Full,shift> (t);
  }
#endif
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> 
  void Interaction_Forces_gt1(const StInterParmsg& t)
{
#ifdef FAST_COMPILATION
  if(t.tvisco!=VISCO_Artificial)throw "Extra viscosity options are disabled for FastCompilation...";
  Interaction_Forces_gt2<tker,ftmode,VISCO_Artificial> (t);
#else
  if(t.tvisco==VISCO_Artificial)     Interaction_Forces_gt2<tker,ftmode,VISCO_Artificial>(t);
  else if(t.tvisco==VISCO_Laminar)   Interaction_Forces_gt2<tker,ftmode,VISCO_Laminar>   (t);
  else if(t.tvisco==VISCO_LaminarSPS)Interaction_Forces_gt2<tker,ftmode,VISCO_LaminarSPS>(t);
#endif
}
//==============================================================================
template<TpKernel tker> void Interaction_Forces_gt0(const StInterParmsg& t){
#ifdef FAST_COMPILATION
  if(t.ftmode!=FTMODE_None)throw "Extra FtMode options are disabled for FastCompilation...";
  Interaction_Forces_gt1<tker,FTMODE_None> (t);
#else
  if(t.ftmode==FTMODE_None)    Interaction_Forces_gt1<tker,FTMODE_None> (t);
  else if(t.ftmode==FTMODE_Sph)Interaction_Forces_gt1<tker,FTMODE_Sph>  (t);
  else if(t.ftmode==FTMODE_Ext)Interaction_Forces_gt1<tker,FTMODE_Ext>  (t);
#endif
}
//==============================================================================
void Interaction_Forces(const StInterParmsg& t){
#ifdef FAST_COMPILATION
  if(t.tkernel!=KERNEL_Wendland)throw "Extra kernels are disabled for FastCompilation...";
  Interaction_Forces_gt0<KERNEL_Wendland> (t);
#else
  if(t.tkernel==KERNEL_Wendland)     Interaction_Forces_gt0<KERNEL_Wendland> (t);
 #ifndef DISABLE_KERNELS_EXTRA
  else if(t.tkernel==KERNEL_Cubic)   Interaction_Forces_gt0<KERNEL_Cubic   > (t);
 #endif
#endif
}

//------------------------------------------------------------------------------
/// Returns the corrected position after applying periodic conditions.
/// Devuelve la posicion corregida tras aplicar condiciones periodicas.
//------------------------------------------------------------------------------
__device__ float4 KerComputePosCell(const double3& ps,const double3& mapposmin
  ,float poscellsize)
{
  const double dx=ps.x-mapposmin.x;
  const double dy=ps.y-mapposmin.y;
  const double dz=ps.z-mapposmin.z;
  int cx=int(dx/poscellsize);
  int cy=int(dy/poscellsize);
  int cz=int(dz/poscellsize);
  cx=(cx>=0? cx: 0);
  cy=(cy>=0? cy: 0);
  cz=(cz>=0? cz: 0);
  const float px=float(dx-(double(poscellsize)*cx));
  const float py=float(dy-(double(poscellsize)*cy));
  const float pz=float(dz-(double(poscellsize)*cz));
  const float pw=__uint_as_float(PSCEL_Code(cx,cy,cz));
  return(make_float4(px,py,pz,pw));
}

//------------------------------------------------------------------------------
/// Perform interaction between ghost node of selected bondary and fluid.
//------------------------------------------------------------------------------
template<TpKernel tker,bool sim2d,TpSlipMode tslip>
  __global__ void KerInteractionMdbcCorrection_Fast
  (unsigned n,unsigned nbound,float determlimit,float mdbcthreshold
  ,double3 mapposmin,float poscellsize,const float4* poscell
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,const double2* posxy,const double* posz,const typecode* code,const unsigned* idp
  ,const float3* boundnor,const float3* motionvel,float4* velrho
  ,const float3* motionace,float* boundonoff,const tfloat3 gravity)
{
  const unsigned p1=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p1<n){
    const float3 bnormalp1=boundnor[p1];
    if(bnormalp1.x!=0 || bnormalp1.y!=0 || bnormalp1.z!=0){
      float rhofinal=FLT_MAX;
      float3 velrhofinal=make_float3(0,0,0);
      float sumwab=0.0f;
      float submerged = 0.0f; //SHABA
      float minrho=1000.0f;


      //-Calculates ghost node position.
      double3 gposp1=make_double3(posxy[p1].x+bnormalp1.x,posxy[p1].y+bnormalp1.y,posz[p1]+bnormalp1.z);
      gposp1=(CTE.periactive!=0? KerUpdatePeriodicPos(gposp1): gposp1); //-Corrected interface Position.
      const float4 gpscellp1=KerComputePosCell(gposp1,mapposmin,poscellsize);

      //-Initializes variables for calculation.
      float rhop1=0;
      float3 gradrhop1=make_float3(0,0,0);
      float3 velp1=make_float3(0,0,0);                              // -Only for velocity
      tmatrix3f a_corr2; if(sim2d) cumath::Tmatrix3fReset(a_corr2); //-Only for 2D.
      tmatrix4f a_corr3; if(!sim2d)cumath::Tmatrix4fReset(a_corr3); //-Only for 3D.
    
      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::InitCte(gposp1.x,gposp1.y,gposp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      //-Boundary-Fluid interaction.
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
        if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
          const float4 pscellp2=poscell[p2];
          float drx=gpscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(gpscellp1.w)-PSCEL_GetfX(pscellp2.w));
          float dry=gpscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(gpscellp1.w)-PSCEL_GetfY(pscellp2.w));
          float drz=gpscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(gpscellp1.w)-PSCEL_GetfZ(pscellp2.w));
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=CTE.kernelsize2 && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
            //-Computes kernel.
            float fac;
            const float wab=cufsph::GetKernel_WabFac<tker>(rr2,fac);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

            //===== Get mass and volume of particle p2 =====
            const float4 velrhop2=velrho[p2];
            float massp2=CTE.massf;
            const float volp2=massp2/velrhop2.w;
            minrho=fminf(minrho,velrhop2.w);


            //===== Check if ghost node is submerged ===== SHABA
            submerged -= volp2 * (drx * frx + dry * fry + drz * frz);

            //===== Density and its gradient =====
            rhop1+=massp2*wab;
            gradrhop1.x+=massp2*frx;
            gradrhop1.y+=massp2*fry;
            gradrhop1.z+=massp2*frz;

            //===== Kernel values multiplied by volume =====
            const float vwab=wab*volp2;
            sumwab+=vwab;
            const float vfrx=frx*volp2;
            const float vfry=fry*volp2;
            const float vfrz=frz*volp2;

            //===== Velocity =====
            if(tslip!=SLIP_Vel0) {
              velp1.x+=vwab*velrhop2.x;
              velp1.y+=vwab*velrhop2.y;
              velp1.z+=vwab*velrhop2.z;
            }

            //===== Matrix A for correction =====
            if(sim2d){
              a_corr2.a11+=vwab;  a_corr2.a12+=drx*vwab;  a_corr2.a13+=drz*vwab;
              a_corr2.a21+=vfrx;  a_corr2.a22+=drx*vfrx;  a_corr2.a23+=drz*vfrx;
              a_corr2.a31+=vfrz;  a_corr2.a32+=drx*vfrz;  a_corr2.a33+=drz*vfrz;
            }
            else{
              a_corr3.a11+=vwab;  a_corr3.a12+=drx*vwab;  a_corr3.a13+=dry*vwab;  a_corr3.a14+=drz*vwab;
              a_corr3.a21+=vfrx;  a_corr3.a22+=drx*vfrx;  a_corr3.a23+=dry*vfrx;  a_corr3.a24+=drz*vfrx;
              a_corr3.a31+=vfry;  a_corr3.a32+=drx*vfry;  a_corr3.a33+=dry*vfry;  a_corr3.a34+=drz*vfry;
              a_corr3.a41+=vfrz;  a_corr3.a42+=drx*vfrz;  a_corr3.a43+=dry*vfrz;  a_corr3.a44+=drz*vfrz;
            }
          }
        }
      }

      //-Store the results.
      //--------------------
      if (submerged > 0.f) {
          //-Store the results.
          //--------------------
          //if (sumwab >= mdbcthreshold) {
          boundonoff[p1] = 1.0f;
          const float3 dpos = make_float3(-bnormalp1.x, -bnormalp1.y, -bnormalp1.z); //-Boundary particle position - ghost node position.
          if (sim2d) {
              if (sumwab < 0.1f) { // If kernel sum is small use shepherd density and vel0
                  // trying to avoid too negative pressures
                  const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr2.a11));
                  // Clone particle proceedure
                  float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                  float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                  float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                  float normpos = dpos.x * normx + dpos.z * normz;
                  float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                  float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                  float pressfinal = pghost + normforce * normpos;
                  // final values to be saved
                  rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                  if (tslip != SLIP_Vel0) {
                      velrhofinal.x = float(velp1.x / a_corr2.a11);
                      velrhofinal.y = 0.f;
                      velrhofinal.z = float(velp1.z / a_corr2.a11);
                  }
              }
              else {// chech if matrix is invertible and well conditioned
                  const double determ = cumath::Determinant3x3dbl(a_corr2);
                  if (fabs(determ) >= 0.001) {//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
                      const tmatrix3f invacorr2 = cumath::InverseMatrix3x3dbl(a_corr2, determ);
                      // finding inifinity norm of matrix a
                      float row1 = float(fabs(a_corr2.a11) + fabs(a_corr2.a12) + fabs(a_corr2.a13));
                      float row2 = float(fabs(a_corr2.a21) + fabs(a_corr2.a22) + fabs(a_corr2.a23));
                      float row3 = float(fabs(a_corr2.a31) + fabs(a_corr2.a32) + fabs(a_corr2.a33));
                      const float infnorma = max(row1, max(row2, row3));
                      // finding infinity norm of matrix a inverse
                      float invrow1 = float(fabs(invacorr2.a11) + fabs(invacorr2.a12) + fabs(invacorr2.a13));
                      float invrow2 = float(fabs(invacorr2.a21) + fabs(invacorr2.a22) + fabs(invacorr2.a23));
                      float invrow3 = float(fabs(invacorr2.a31) + fabs(invacorr2.a32) + fabs(invacorr2.a33));
                      const float infnormainv = max(invrow1, max(invrow2, invrow3));
                      // calculate the scaled condition number
                      const float condinf = float(CTE.dp) * float(CTE.dp) * infnorma * infnormainv;
                      if (condinf <= 50.f) {// if matrix is well conditioned use matrix inverse for density and shepherd for velocity
                          // trying to avoid too negative pressure
                          const float rhoghost = max(CTE.rhopzero, float(invacorr2.a11 * rhop1 + invacorr2.a12 * gradrhop1.x + invacorr2.a13 * gradrhop1.z));
                          // clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final values to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr2.a11);
                              velrhofinal.y = 0.f;
                              velrhofinal.z = float(velp1.z / a_corr2.a11);
                          }
                      }
                      else { // if ill conditioned use shepherd
                          // trying to avoid too negative pressures
                          const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr2.a11));
                          // Clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final values to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr2.a11);
                              velrhofinal.y = 0.f;
                              velrhofinal.z = float(velp1.z / a_corr2.a11);
                          }
                      }
                  }
                  else {// if not invertible use shepherd
                      // trying to avoid too negative pressures
                      const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr2.a11));
                      // Clone particle proceedure
                      float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                      float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                      float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                      float normpos = dpos.x * normx + dpos.z * normz;
                      float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                      float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                      float pressfinal = pghost + normforce * normpos;
                      // final values to be saved
                      rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                      if (tslip != SLIP_Vel0) {
                          velrhofinal.x = float(velp1.x / a_corr2.a11);
                          velrhofinal.y = 0.f;
                          velrhofinal.z = float(velp1.z / a_corr2.a11);
                      }
                  }
              }
          }
          else {
              if (sumwab < 0.1f) { // If kernel sum is small use shepherd density and vel0
                  // trying to avoid too negative pressures
                  const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr3.a11));
                  // Clone particle proceedure
                  float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                  float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                  float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                  float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                  float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                  float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                  float pressfinal = pghost + normforce * normpos;
                  // final values to be saved
                  rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                  if (tslip != SLIP_Vel0) {
                      velrhofinal.x = float(velp1.x / a_corr3.a11);
                      velrhofinal.y = float(velp1.y / a_corr3.a11);
                      velrhofinal.z = float(velp1.z / a_corr3.a11);
                  }
              }
              else {
                  const double determ = cumath::Determinant4x4dbl(a_corr3);
                  if (fabs(determ) >= 0.001) {
                      const tmatrix4f invacorr3 = cumath::InverseMatrix4x4dbl(a_corr3, determ);
                      // finding inifinity norm of matrix a
                      float row1 = float(fabs(a_corr3.a11) + fabs(a_corr3.a12) + fabs(a_corr3.a13) + fabs(a_corr3.a14));
                      float row2 = float(fabs(a_corr3.a21) + fabs(a_corr3.a22) + fabs(a_corr3.a23) + fabs(a_corr3.a24));
                      float row3 = float(fabs(a_corr3.a31) + fabs(a_corr3.a32) + fabs(a_corr3.a33) + fabs(a_corr3.a34));
                      float row4 = float(fabs(a_corr3.a41) + fabs(a_corr3.a42) + fabs(a_corr3.a43) + fabs(a_corr3.a44));
                      const float infnorma = max(row1, max(row2, max(row3, row4)));
                      // finding infinity norm of matrix a inverse
                      float invrow1 = float(fabs(invacorr3.a11) + fabs(invacorr3.a12) + fabs(invacorr3.a13) + fabs(invacorr3.a14));
                      float invrow2 = float(fabs(invacorr3.a21) + fabs(invacorr3.a22) + fabs(invacorr3.a23) + fabs(invacorr3.a24));
                      float invrow3 = float(fabs(invacorr3.a31) + fabs(invacorr3.a32) + fabs(invacorr3.a33) + fabs(invacorr3.a34));
                      float invrow4 = float(fabs(invacorr3.a41) + fabs(invacorr3.a42) + fabs(invacorr3.a43) + fabs(invacorr3.a44));
                      const float infnormainv = max(invrow1, max(invrow2, max(invrow3, invrow4)));
                      // calculate the scaled condition number
                      const float condinf = float(CTE.dp) * float(CTE.dp) * infnorma * infnormainv;
                      if (condinf <= 50.f) { // if matrix is well conditioned use matrix inverse for density and shepherd for velocity
                          // trying to avoid too negative pressure
                          const float rhoghost = max(CTE.rhopzero, float(invacorr3.a11 * rhop1 + invacorr3.a12 * gradrhop1.x + invacorr3.a13 * gradrhop1.y + invacorr3.a14 * gradrhop1.z));
                          // clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final values to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr3.a11);
                              velrhofinal.y = float(velp1.y / a_corr3.a11);
                              velrhofinal.z = float(velp1.z / a_corr3.a11);
                          }
                      }
                      else { // matrix is not well conditioned, use shepherd for pressure and velocity
                          // trying to avoid too negative pressure
                          const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr3.a11));
                          // clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final value to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr3.a11);
                              velrhofinal.y = float(velp1.y / a_corr3.a11);
                              velrhofinal.z = float(velp1.z / a_corr3.a11);
                          }
                      }
                  }
                  else { // matrix is not invertible use shepherd for density and velocity
                      // trying to avoid too negative pressure
                      const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr3.a11));
                      // clone particle proceedure
                      float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                      float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                      float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                      float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                      float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                      float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                      float pressfinal = pghost + normforce * normpos;
                      // final values to be saved
                      rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                      if (tslip != SLIP_Vel0) {
                          velrhofinal.x = float(velp1.x / a_corr3.a11);
                          velrhofinal.y = float(velp1.y / a_corr3.a11);
                          velrhofinal.z = float(velp1.z / a_corr3.a11);
                      }
                  }
              }
          }
          //-Store the results.
          rhofinal = (rhofinal != FLT_MAX ? rhofinal : CTE.rhopzero);
          if (tslip == SLIP_Vel0) {//-DBC vel=0
              velrho[p1].w = rhofinal;
          }
          if (tslip == SLIP_NoSlip) {//-No-Slip
              const float3 v = motionvel[p1];
              velrho[p1] = make_float4(v.x + v.x - velrhofinal.x, v.y + v.y - velrhofinal.y, v.z + v.z - velrhofinal.z, rhofinal);
          }
          if (tslip == SLIP_FreeSlip) {//-No-Penetration and free slip    SHABA
              float3 FSVelFinal; // final free slip boundary velocity
              const float3 v = motionvel[p1];
              float motion = sqrt(v.x * v.x + v.y * v.y + v.z * v.z); // to check if boundary moving
              float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
              float3 normal; // creating a normailsed boundary normal
              normal.x = fabs(bnormalp1.x) / norm; normal.y = fabs(bnormalp1.y) / norm; normal.z = fabs(bnormalp1.z) / norm;

              // finding the velocity componants normal and tangential to boundary 
              float3 normvel = make_float3(velrhofinal.x * normal.x, velrhofinal.y * normal.y, velrhofinal.z * normal.z); // velocity in direction of normal pointin ginto fluid)
              float3 tangvel = make_float3(velrhofinal.x - normvel.x, velrhofinal.y - normvel.y, velrhofinal.z - normvel.z); // velocity tangential to normal

              if (motion > 0) { // if moving boundary
                  float3 normmot = make_float3(v.x * normal.x, v.y * normal.y, v.z * normal.z); // boundary motion in direction normal to boundary 
                  FSVelFinal = make_float3(normmot.x + normmot.x - normvel.x, normmot.y + normmot.y - normvel.y, normmot.z + normmot.z - normvel.z);
                  // only velocity in normal direction for no-penetration
                  // fluid sees zero velocity in the tangetial direction
              }
              else {
                  FSVelFinal = make_float3(tangvel.x - normvel.x, tangvel.y - normvel.y, tangvel.z - normvel.z);
                  // tangential velocity equal to fluid velocity for free slip
                  // normal velocity reversed for no-penetration
              }

              // Save the velocity and density
              velrho[p1] = make_float4(FSVelFinal.x, FSVelFinal.y, FSVelFinal.z, rhofinal);
          }
          //}
      }
      else { // if unsubmerged switch off boundary particle
          boundonoff[p1] = 0.0f;
          const float3 v = motionvel[p1];
          velrho[p1] = make_float4(v.x, v.y, v.z, CTE.rhopzero);
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Perform interaction between ghost node of selected bondary and fluid.
//------------------------------------------------------------------------------
template<TpKernel tker,bool sim2d,TpSlipMode tslip>
  __global__ void KerInteractionMdbcCorrection_Dbl
  (unsigned n,unsigned nbound,float determlimit,float mdbcthreshold
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,const double2* posxy,const double* posz,const typecode* code,const unsigned* idp
  ,const float3* boundnor,const float3* motionvel,float4* velrho
  ,const float3* motionace,float* boundonoff,const tfloat3 gravity)
{
  const unsigned p1=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p1<n){
    const float3 bnormalp1=boundnor[p1];
    if(bnormalp1.x!=0 || bnormalp1.y!=0 || bnormalp1.z!=0){
      float rhofinal=FLT_MAX;
      float3 velrhofinal=make_float3(0,0,0);
      float sumwab=0;
      float submerged = 0; // SHABA


      //-Calculates ghost node position.
      double3 gposp1=make_double3(posxy[p1].x+bnormalp1.x,posxy[p1].y+bnormalp1.y,posz[p1]+bnormalp1.z);
      gposp1=(CTE.periactive!=0? KerUpdatePeriodicPos(gposp1): gposp1); //-Corrected interface Position.
      //-Initializes variables for calculation.
      float rhop1=0;
      float3 gradrhop1=make_float3(0,0,0);
      float3 velp1=make_float3(0,0,0);                              // -Only for velocity
      tmatrix3d a_corr2; if(sim2d) cumath::Tmatrix3dReset(a_corr2); //-Only for 2D.
      tmatrix4d a_corr3; if(!sim2d)cumath::Tmatrix4dReset(a_corr3); //-Only for 3D.
    
      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::InitCte(gposp1.x,gposp1.y,gposp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      //-Boundary-Fluid interaction.
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
        if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
          const double2 p2xy=posxy[p2];
          const float drx=float(gposp1.x-p2xy.x);
          const float dry=float(gposp1.y-p2xy.y);
          const float drz=float(gposp1.z-posz[p2]);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=CTE.kernelsize2 && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
            //-Computes kernel.
            float fac;
            const float wab=cufsph::GetKernel_WabFac<tker>(rr2,fac);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

            //===== Get mass and volume of particle p2 =====
            const float4 velrhop2=velrho[p2];
            float massp2=CTE.massf;
            const float volp2=massp2/velrhop2.w;


            //===== Check if ghost node is submerged ===== SHABA
            submerged -= volp2 * (drx * frx + dry * fry + drz * frz);

            //===== Density and its gradient =====
            rhop1+=massp2*wab;
            gradrhop1.x+=massp2*frx;
            gradrhop1.y+=massp2*fry;
            gradrhop1.z+=massp2*frz;

            //===== Kernel values multiplied by volume =====
            const float vwab=wab*volp2;
            sumwab+=vwab;
            const float vfrx=frx*volp2;
            const float vfry=fry*volp2;
            const float vfrz=frz*volp2;

            //===== Velocity =====
            if(tslip!=SLIP_Vel0) {
              velp1.x+=vwab*velrhop2.x;
              velp1.y+=vwab*velrhop2.y;
              velp1.z+=vwab*velrhop2.z;
            }

            //===== Matrix A for correction =====
            if(sim2d){
              a_corr2.a11+=vwab;  a_corr2.a12+=drx*vwab;  a_corr2.a13+=drz*vwab;
              a_corr2.a21+=vfrx;  a_corr2.a22+=drx*vfrx;  a_corr2.a23+=drz*vfrx;
              a_corr2.a31+=vfrz;  a_corr2.a32+=drx*vfrz;  a_corr2.a33+=drz*vfrz;
            }
            else{
              a_corr3.a11+=vwab;  a_corr3.a12+=drx*vwab;  a_corr3.a13+=dry*vwab;  a_corr3.a14+=drz*vwab;
              a_corr3.a21+=vfrx;  a_corr3.a22+=drx*vfrx;  a_corr3.a23+=dry*vfrx;  a_corr3.a24+=drz*vfrx;
              a_corr3.a31+=vfry;  a_corr3.a32+=drx*vfry;  a_corr3.a33+=dry*vfry;  a_corr3.a34+=drz*vfry;
              a_corr3.a41+=vfrz;  a_corr3.a42+=drx*vfrz;  a_corr3.a43+=dry*vfrz;  a_corr3.a44+=drz*vfrz;
            }
          }
        }
      }

    //-Store the results.
      //--------------------
      if (submerged > 0.f) {
          //-Store the results.
          //--------------------
          //if (sumwab >= mdbcthreshold) {
          boundonoff[p1] = 1.0f;
          const float3 dpos = make_float3(-bnormalp1.x, -bnormalp1.y, -bnormalp1.z); //-Boundary particle position - ghost node position.
          if (sim2d) {
              if (sumwab < 0.1f) { // If kernel sum is small use shepherd density and vel0
                  // trying to avoid too negative pressures
                  const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr2.a11));
                  // Clone particle proceedure
                  float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                  float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                  float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                  float normpos = dpos.x * normx + dpos.z * normz;
                  float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                  float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                  float pressfinal = pghost + normforce * normpos;
                  // final values to be saved
                  rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                  if (tslip != SLIP_Vel0) {
                      velrhofinal.x = float(velp1.x / a_corr2.a11);
                      velrhofinal.y = 0.f;
                      velrhofinal.z = float(velp1.z / a_corr2.a11);
                  }
              }
              else {// chech if matrix is invertible and well conditioned
                  const double determ = cumath::Determinant3x3(a_corr2);
                  if (fabs(determ) >= 0.001) {//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
                      const tmatrix3d invacorr2 = cumath::InverseMatrix3x3(a_corr2, determ);
                      // finding inifinity norm of matrix a
                      float row1 = float(fabs(a_corr2.a11) + fabs(a_corr2.a12) + fabs(a_corr2.a13));
                      float row2 = float(fabs(a_corr2.a21) + fabs(a_corr2.a22) + fabs(a_corr2.a23));
                      float row3 = float(fabs(a_corr2.a31) + fabs(a_corr2.a32) + fabs(a_corr2.a33));
                      const float infnorma = max(row1, max(row2, row3));
                      // finding infinity norm of matrix a inverse
                      float invrow1 = float(fabs(invacorr2.a11) + fabs(invacorr2.a12) + fabs(invacorr2.a13));
                      float invrow2 = float(fabs(invacorr2.a21) + fabs(invacorr2.a22) + fabs(invacorr2.a23));
                      float invrow3 = float(fabs(invacorr2.a31) + fabs(invacorr2.a32) + fabs(invacorr2.a33));
                      const float infnormainv = max(invrow1, max(invrow2, invrow3));
                      // calculate the scaled condition number
                      const float condinf = float(CTE.dp) * float(CTE.dp) * infnorma * infnormainv;
                      if (condinf <= 50.f) {// if matrix is well conditioned use matrix inverse for density and shepherd for velocity
                          // trying to avoid too negative pressure
                          const float rhoghost = max(CTE.rhopzero, float(invacorr2.a11 * rhop1 + invacorr2.a12 * gradrhop1.x + invacorr2.a13 * gradrhop1.z));
                          // clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final values to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr2.a11);
                              velrhofinal.y = 0.f;
                              velrhofinal.z = float(velp1.z / a_corr2.a11);
                          }
                      }
                      else { // if ill conditioned use shepherd
                          // trying to avoid too negative pressures
                          const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr2.a11));
                          // Clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final values to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr2.a11);
                              velrhofinal.y = 0.f;
                              velrhofinal.z = float(velp1.z / a_corr2.a11);
                          }
                      }
                  }
                  else {// if not invertible use shepherd
                      // trying to avoid too negative pressures
                      const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr2.a11));
                      // Clone particle proceedure
                      float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                      float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.z * bnormalp1.z);
                      float normx = bnormalp1.x / norm; float normz = bnormalp1.z / norm;
                      float normpos = dpos.x * normx + dpos.z * normz;
                      float3 force = make_float3(gravity.x - motionace[p1].x, 0, gravity.z - motionace[p1].z);
                      float normforce = CTE.rhopzero * (force.x * normx + force.z * normz);
                      float pressfinal = pghost + normforce * normpos;
                      // final values to be saved
                      rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                      if (tslip != SLIP_Vel0) {
                          velrhofinal.x = float(velp1.x / a_corr2.a11);
                          velrhofinal.y = 0.f;
                          velrhofinal.z = float(velp1.z / a_corr2.a11);
                      }
                  }
              }
          }
          else {
              if (sumwab < 0.1f) { // If kernel sum is small use shepherd density and vel0
                  // trying to avoid too negative pressures
                  const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr3.a11));
                  // Clone particle proceedure
                  float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                  float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                  float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                  float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                  float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                  float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                  float pressfinal = pghost + normforce * normpos;
                  // final values to be saved
                  rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                  if (tslip != SLIP_Vel0) {
                      velrhofinal.x = float(velp1.x / a_corr3.a11);
                      velrhofinal.y = float(velp1.y / a_corr3.a11);
                      velrhofinal.z = float(velp1.z / a_corr3.a11);
                  }
              }
              else {
                  const double determ = cumath::Determinant4x4(a_corr3);
                  if (fabs(determ) >= 0.001) {
                      const tmatrix4d invacorr3 = cumath::InverseMatrix4x4(a_corr3, determ);
                      // finding inifinity norm of matrix a
                      float row1 = float(fabs(a_corr3.a11) + fabs(a_corr3.a12) + fabs(a_corr3.a13) + fabs(a_corr3.a14));
                      float row2 = float(fabs(a_corr3.a21) + fabs(a_corr3.a22) + fabs(a_corr3.a23) + fabs(a_corr3.a24));
                      float row3 = float(fabs(a_corr3.a31) + fabs(a_corr3.a32) + fabs(a_corr3.a33) + fabs(a_corr3.a34));
                      float row4 = float(fabs(a_corr3.a41) + fabs(a_corr3.a42) + fabs(a_corr3.a43) + fabs(a_corr3.a44));
                      const float infnorma = max(row1, max(row2, max(row3, row4)));
                      // finding infinity norm of matrix a inverse
                      float invrow1 = float(fabs(invacorr3.a11) + fabs(invacorr3.a12) + fabs(invacorr3.a13) + fabs(invacorr3.a14));
                      float invrow2 = float(fabs(invacorr3.a21) + fabs(invacorr3.a22) + fabs(invacorr3.a23) + fabs(invacorr3.a24));
                      float invrow3 = float(fabs(invacorr3.a31) + fabs(invacorr3.a32) + fabs(invacorr3.a33) + fabs(invacorr3.a34));
                      float invrow4 = float(fabs(invacorr3.a41) + fabs(invacorr3.a42) + fabs(invacorr3.a43) + fabs(invacorr3.a44));
                      const float infnormainv = max(invrow1, max(invrow2, max(invrow3, invrow4)));
                      // calculate the scaled condition number
                      const float condinf = float(CTE.dp) * float(CTE.dp) * infnorma * infnormainv;
                      if (condinf <= 50.f) { // if matrix is well conditioned use matrix inverse for density and shepherd for velocity
                          // trying to avoid too negative pressure
                          const float rhoghost = max(CTE.rhopzero, float(invacorr3.a11 * rhop1 + invacorr3.a12 * gradrhop1.x + invacorr3.a13 * gradrhop1.y + invacorr3.a14 * gradrhop1.z));
                          // clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final values to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr3.a11);
                              velrhofinal.y = float(velp1.y / a_corr3.a11);
                              velrhofinal.z = float(velp1.z / a_corr3.a11);
                          }
                      }
                      else { // matrix is not well conditioned, use shepherd for pressure and velocity
                          // trying to avoid too negative pressure
                          const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr3.a11));
                          // clone particle proceedure
                          float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                          float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                          float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                          float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                          float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                          float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                          float pressfinal = pghost + normforce * normpos;
                          // final value to be saved
                          rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                          if (tslip != SLIP_Vel0) {
                              velrhofinal.x = float(velp1.x / a_corr3.a11);
                              velrhofinal.y = float(velp1.y / a_corr3.a11);
                              velrhofinal.z = float(velp1.z / a_corr3.a11);
                          }
                      }
                  }
                  else { // matrix is not invertible use shepherd for density and velocity
                      // trying to avoid too negative pressure
                      const float rhoghost = max(CTE.rhopzero, float(rhop1 / a_corr3.a11));
                      // clone particle proceedure
                      float pghost = float(CTE.cs0 * CTE.cs0 * (rhoghost - CTE.rhopzero));
                      float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
                      float normx = bnormalp1.x / norm; float normy = bnormalp1.y / norm; float normz = bnormalp1.z / norm;
                      float normpos = dpos.x * normx + dpos.y * normy + dpos.z * normz;
                      float3 force = make_float3(gravity.x - motionace[p1].x, gravity.y - motionace[p1].y, gravity.z - motionace[p1].z);
                      float normforce = CTE.rhopzero * (force.x * normx + force.y * normy + force.z * normz);
                      float pressfinal = pghost + normforce * normpos;
                      // final values to be saved
                      rhofinal = CTE.rhopzero + float(pressfinal / (CTE.cs0 * CTE.cs0));
                      if (tslip != SLIP_Vel0) {
                          velrhofinal.x = float(velp1.x / a_corr3.a11);
                          velrhofinal.y = float(velp1.y / a_corr3.a11);
                          velrhofinal.z = float(velp1.z / a_corr3.a11);
                      }
                  }
              }
          }
          //-Store the results.
          rhofinal = (rhofinal != FLT_MAX ? rhofinal : CTE.rhopzero);
          if (tslip == SLIP_Vel0) {//-DBC vel=0
              velrho[p1].w = rhofinal;
          }
          if (tslip == SLIP_NoSlip) {//-No-Slip
              const float3 v = motionvel[p1];
              velrho[p1] = make_float4(v.x + v.x - velrhofinal.x, v.y + v.y - velrhofinal.y, v.z + v.z - velrhofinal.z, rhofinal);
          }
          if (tslip == SLIP_FreeSlip) {//-No-Penetration and free slip    SHABA
              float3 FSVelFinal; // final free slip boundary velocity
              const float3 v = motionvel[p1];
              float motion = sqrt(v.x * v.x + v.y * v.y + v.z * v.z); // to check if boundary moving
              float norm = sqrt(bnormalp1.x * bnormalp1.x + bnormalp1.y * bnormalp1.y + bnormalp1.z * bnormalp1.z);
              float3 normal; // creating a normailsed boundary normal
              normal.x = fabs(bnormalp1.x) / norm; normal.y = fabs(bnormalp1.y) / norm; normal.z = fabs(bnormalp1.z) / norm;

              // finding the velocity componants normal and tangential to boundary 
              float3 normvel = make_float3(velrhofinal.x * normal.x, velrhofinal.y * normal.y, velrhofinal.z * normal.z); // velocity in direction of normal pointin ginto fluid)
              float3 tangvel = make_float3(velrhofinal.x - normvel.x, velrhofinal.y - normvel.y, velrhofinal.z - normvel.z); // velocity tangential to normal

              if (motion > 0) { // if moving boundary
                  float3 normmot = make_float3(v.x * normal.x, v.y * normal.y, v.z * normal.z); // boundary motion in direction normal to boundary 
                  FSVelFinal = make_float3(normmot.x + normmot.x - normvel.x, normmot.y + normmot.y - normvel.y, normmot.z + normmot.z - normvel.z);
                  // only velocity in normal direction for no-penetration
                  // fluid sees zero velocity in the tangetial direction
              }
              else {
                  FSVelFinal = make_float3(tangvel.x - normvel.x, tangvel.y - normvel.y, tangvel.z - normvel.z);
                  // tangential velocity equal to fluid velocity for free slip
                  // normal velocity reversed for no-penetration
              }

              // Save the velocity and density
              velrho[p1] = make_float4(FSVelFinal.x, FSVelFinal.y, FSVelFinal.z, rhofinal);
          }
          //}
      }
      else { // if unsubmerged switch off boundary particle
          boundonoff[p1] = 0.0f;
          const float3 v = motionvel[p1];
          velrho[p1] = make_float4(v.x, v.y, v.z, CTE.rhopzero);
      }
    }
  }
}


//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
template<TpKernel tker,bool sim2d,TpSlipMode tslip> void Interaction_MdbcCorrectionT2(
  bool fastsingle,unsigned n,unsigned nbound,float mdbcthreshold
  ,const StDivDataGpu& dvd,const tdouble3& mapposmin,const double2* posxy
  ,const double* posz,const float4* poscell,const typecode* code
  ,const unsigned* idp,const float3* boundnor,const float3* motionvel
  ,float4* velrho,const float3* motionace,float* boundonoff,const tfloat3 gravity
  ,cudaStream_t stm)
{
  const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
  const float determlimit=1e-3f;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(n){
    const unsigned bsbound=128;
    dim3 sgridb=cusph::GetSimpleGridSize(n,bsbound);
    if(fastsingle){//-mDBC-Fast_v2
      KerInteractionMdbcCorrection_Fast <tker,sim2d,tslip> <<<sgridb,bsbound,0,stm>>>
        (n,nbound,determlimit,mdbcthreshold,Double3(mapposmin),dvd.poscellsize
        ,poscell,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
        ,posxy,posz,code,idp,boundnor,motionvel,velrho
        ,motionace,boundonoff,gravity);
    }
    else{//-mDBC_v0
      KerInteractionMdbcCorrection_Dbl <tker,sim2d,tslip> <<<sgridb,bsbound,0,stm>>>
        (n,nbound,determlimit,mdbcthreshold,dvd.scelldiv,dvd.nc,dvd.cellzero
        ,beginendcellfluid,posxy,posz,code,idp,boundnor,motionvel,velrho
        ,motionace,boundonoff,gravity);
    }
  }
}
//==============================================================================
template<TpKernel tker> void Interaction_MdbcCorrectionT(bool simulate2d
  ,TpSlipMode slipmode,bool fastsingle,unsigned n,unsigned nbound
  ,float mdbcthreshold,const StDivDataGpu& dvd,const tdouble3& mapposmin
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const typecode* code,const unsigned* idp,const float3* boundnor
  ,const float3* motionvel,float4* velrho
  ,const float3* motionace,float* boundonoff,const tfloat3 gravity
  ,cudaStream_t stm)
{
  switch(slipmode){
    case SLIP_Vel0:{ const TpSlipMode tslip=SLIP_Vel0;
      if(simulate2d)Interaction_MdbcCorrectionT2 <tker,true ,tslip> (fastsingle,n,nbound,mdbcthreshold,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho,motionace,boundonoff,gravity,stm);
      else          Interaction_MdbcCorrectionT2 <tker,false,tslip> (fastsingle,n,nbound,mdbcthreshold,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho,motionace,boundonoff,gravity,stm);
    }break;
#ifndef DISABLE_MDBC_EXTRAMODES
    case SLIP_NoSlip:{ const TpSlipMode tslip=SLIP_NoSlip;
      if(simulate2d)Interaction_MdbcCorrectionT2 <tker,true ,tslip> (fastsingle,n,nbound,mdbcthreshold,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho,motionace,boundonoff,gravity,stm);
      else          Interaction_MdbcCorrectionT2 <tker,false,tslip> (fastsingle,n,nbound,mdbcthreshold,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho,motionace,boundonoff,gravity,stm);
    }break;
    case SLIP_FreeSlip:{ const TpSlipMode tslip=SLIP_FreeSlip;
      if(simulate2d)Interaction_MdbcCorrectionT2 <tker,true ,tslip> (fastsingle,n,nbound,mdbcthreshold,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho,motionace,boundonoff,gravity,stm);
      else          Interaction_MdbcCorrectionT2 <tker,false,tslip> (fastsingle,n,nbound,mdbcthreshold,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho,motionace,boundonoff,gravity,stm);
    }break;
#endif
    default: throw "SlipMode unknown at Interaction_MdbcCorrectionT().";
  }
}
//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void Interaction_MdbcCorrection(TpKernel tkernel,bool simulate2d,TpSlipMode slipmode
  ,bool fastsingle,unsigned n,unsigned nbound,float mdbcthreshold
  ,const StDivDataGpu& dvd,const tdouble3& mapposmin
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const typecode* code,const unsigned* idp,const float3* boundnor
  ,const float3* motionvel,float4* velrho
  ,const float3* motionace,float* boundonoff,const tfloat3 gravity
  ,cudaStream_t stm)
{
  switch(tkernel){
    case KERNEL_Wendland:{ const TpKernel tker=KERNEL_Wendland;
      Interaction_MdbcCorrectionT <tker> (simulate2d,slipmode,fastsingle,n,nbound,mdbcthreshold
        ,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho
        ,motionace,boundonoff,gravity,stm);
    }break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:{ const TpKernel tker=KERNEL_Cubic;
      Interaction_MdbcCorrectionT <tker> (simulate2d,slipmode,fastsingle,n,nbound,mdbcthreshold
        ,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,velrho
        ,motionace,boundonoff,gravity,stm);
    }break;
#endif
    default: throw "Kernel unknown at Interaction_MdbcCorrection().";
  }
}


//##############################################################################
//# Kernels for DEM interaction.
//# Kernels para interaccion DEM.
//##############################################################################
//------------------------------------------------------------------------------
/// DEM interaction of a particle with a set of particles. (Float-Float/Bound)
/// Realiza la interaccion DEM de una particula con un conjunto de ellas. (Float-Float/Bound)
//------------------------------------------------------------------------------
__device__ void KerInteractionForcesDemBox(
  bool boundp2,const unsigned& pini,const unsigned& pfin
  ,const float4* demdata,float dtforce
  ,const float4* poscell,const float4* velrho,const typecode* code
  ,const unsigned* idp,const float4& pscellp1,const float4& velp1
  ,typecode tavp1,float masstotp1,float ftmassp1,float taup1,float kfricp1
  ,float restitup1,float3& acep1,float& demdtp1)
{
  for(int p2=pini;p2<pfin;p2++){
    const typecode codep2=code[p2];
    if(CODE_IsNotFluid(codep2) && tavp1!=CODE_GetTypeAndValue(codep2)){
      const float4 pscellp2=poscell[p2];
      const float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
      const float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
      const float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
      const float rr2=drx*drx+dry*dry+drz*drz;
      const float rad=sqrt(rr2);

      //-Computes maximum value of demdt.
      float4 demdatap2=demdata[CODE_GetTypeAndValue(codep2)];
      const float nu_mass=(boundp2? masstotp1/2: masstotp1*demdatap2.x/(masstotp1+demdatap2.x)); //-With boundary takes the actual mass of floating 1. | Con boundary toma la propia masa del floating 1.
      const float kn=4/(3*(taup1+demdatap2.y))*sqrt(CTE.dp/4); //-Generalized rigidity - Lemieux 2008.
      const float dvx=velp1.x-velrho[p2].x, dvy=velp1.y-velrho[p2].y, dvz=velp1.z-velrho[p2].z; //vji
      const float nx=drx/rad, ny=dry/rad, nz=drz/rad; //-normal_ji             
      const float vn=dvx*nx+dvy*ny+dvz*nz; //-vji.nji    
      const float demvisc=0.2f/(3.21f*(pow(nu_mass/kn,0.4f)*pow(fabs(vn),-0.2f))/40.f);
      if(demdtp1<demvisc)demdtp1=demvisc;

      const float over_lap=1.0f*CTE.dp-rad; //-(ri+rj)-|dij|
      if(over_lap>0.0f){ //-Contact.
        //-Normal.
        const float eij=(restitup1+demdatap2.w)/2;
        const float gn=-(2.0f*log(eij)*sqrt(nu_mass*kn))/(sqrt(float(PI)+log(eij)*log(eij))); //-Generalized damping - Cummins 2010.
        //const float gn=0.08f*sqrt(nu_mass*sqrt(CTE.dp/2)/((taup1+demdatap2.y)/2)); //-generalized damping - Lemieux 2008.
        const float rep=kn*pow(over_lap,1.5f);
        const float fn=rep-gn*pow(over_lap,0.25f)*vn;
        float acef=fn/ftmassp1; //-Divides by the mass of particle to obtain the acceleration.
        acep1.x+=(acef*nx); acep1.y+=(acef*ny); acep1.z+=(acef*nz); //-Force is applied in the normal between the particles.
        //-Tangencial.
        const float dvxt=dvx-vn*nx, dvyt=dvy-vn*ny, dvzt=dvz-vn*nz; //Vji_t
        const float vt=sqrt(dvxt*dvxt + dvyt*dvyt + dvzt*dvzt);
        const float tx=(vt!=0? dvxt/vt: 0), ty=(vt!=0? dvyt/vt: 0), tz=(vt!=0? dvzt/vt: 0); //-Tang vel unit vector.
        const float ft_elast=2*(kn*dtforce-gn)*vt/7; //-Elastic frictional string -->  ft_elast=2*(kn*fdispl-gn*vt)/7; fdispl=dtforce*vt;
        const float kfric_ij=(kfricp1+demdatap2.z)/2;
        float ft=kfric_ij*fn*tanh(8*vt);  //-Coulomb.
        ft=(ft<ft_elast? ft: ft_elast);   //-Not above yield criteria, visco-elastic model.
        acef=ft/ftmassp1; //-Divides by the mass of particle to obtain the acceleration.
        acep1.x+=(acef*tx); acep1.y+=(acef*ty); acep1.z+=(acef*tz);
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction between particles. Fluid/Float-Fluid/Float or Fluid/Float-Bound.
/// Includes artificial/laminar viscosity and normal/DEM floating bodies.
///
/// Realiza interaccion entre particulas. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings normales/dem.
//------------------------------------------------------------------------------
__global__ void KerInteractionForcesDem(unsigned nfloat
  ,int scelldiv,int4 nc,int3 cellzero,const int2* begincell,unsigned cellfluid
  ,const unsigned* dcell,const unsigned* ftridp,const float4* demdata
  ,const float* ftomassp,float dtforce,const float4* poscell,const float4* velrho
  ,const typecode* code,const unsigned* idp,float* viscdt,float3* ace)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<nfloat){
    const unsigned p1=ftridp[p]; //-Number of particle.
    if(p1!=UINT_MAX){
      float demdtp1=0;
      float3 acep1=make_float3(0,0,0);

      //-Obtains basic data of particle p1.
      const float4 pscellp1=poscell[p1];
      const float4 velp1=velrho[p1];
      const typecode cod=code[p1];
      const typecode tavp1=CODE_GetTypeAndValue(cod);
      const float4 rdata=demdata[tavp1];
      const float masstotp1=rdata.x;
      const float taup1=rdata.y;
      const float kfricp1=rdata.z;
      const float restitup1=rdata.w;
      const float ftmassp1=ftomassp[CODE_GetTypeValue(cod)];

      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      //-Interaction with boundaries.
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin)KerInteractionForcesDemBox (true ,pini,pfin,demdata,dtforce,poscell,velrho,code,idp,pscellp1,velp1,tavp1,masstotp1,ftmassp1,taup1,kfricp1,restitup1,acep1,demdtp1);
      }

      //-Interaction with fluids.
      ini3+=cellfluid; fin3+=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin)KerInteractionForcesDemBox (false,pini,pfin,demdata,dtforce,poscell,velrho,code,idp,pscellp1,velp1,tavp1,masstotp1,ftmassp1,taup1,kfricp1,restitup1,acep1,demdtp1);
      }

      //-Stores results.
      if(acep1.x || acep1.y || acep1.z || demdtp1){
        float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
        if(viscdt[p1]<demdtp1)viscdt[p1]=demdtp1;
      }
    }
  }
}

#ifndef DISABLE_BSMODES
//==============================================================================
/// Collects kernel information.
//==============================================================================
void Interaction_ForcesDemT_KerInfo(StKerInfo* kerinfo)
{
#if CUDART_VERSION >= 6050
  {
    typedef void (*fun_ptr)(unsigned,int,int4,int3,const int2*,unsigned,const unsigned*,const unsigned*,const float4*,const float*,float,const float4*,const float4*,const typecode*,const unsigned*,float*,float3*);
    fun_ptr ptr=&KerInteractionForcesDem;
    int qblocksize=0,mingridsize=0;
    cudaOccupancyMaxPotentialBlockSize(&mingridsize,&qblocksize,(void*)ptr,0,0);
    struct cudaFuncAttributes attr;
    cudaFuncGetAttributes(&attr,(void*)ptr);
    kerinfo->forcesdem_bs=qblocksize;
    kerinfo->forcesdem_rg=attr.numRegs;
    kerinfo->forcesdem_bsmax=attr.maxThreadsPerBlock;
    //printf(">> KerInteractionForcesDem  blocksize:%u (%u)\n",qblocksize,0);
  }
  fcuda::Check_CudaErroorFun("Error collecting kernel information.");
#endif
}
#endif

//==============================================================================
/// Interaction for the force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void Interaction_ForcesDem(unsigned bsize,unsigned nfloat
  ,const StDivDataGpu& dvd,const unsigned* dcell
  ,const unsigned* ftridp,const float4* demdata,const float* ftomassp,float dtforce
  ,const float4* poscell,const float4* velrho
  ,const typecode* code,const unsigned* idp,float* viscdt,float3* ace
  ,StKerInfo* kerinfo)
{
  const int2* beginendcell=dvd.beginendcell;
  //-Collects kernel information.
#ifndef DISABLE_BSMODES
  if(kerinfo){
    Interaction_ForcesDemT_KerInfo(kerinfo);
    return;
  }
#endif
  //-Interaction Fluid-Fluid & Fluid-Bound.
  if(nfloat){
    dim3 sgrid=GetSimpleGridSize(nfloat,bsize);
    KerInteractionForcesDem <<<sgrid,bsize>>> (nfloat
      ,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcell,dvd.cellfluid,dcell
      ,ftridp,demdata,ftomassp,dtforce,poscell,velrho,code,idp,viscdt,ace);
  }
}


//##############################################################################
//# Kernels for Laminar+SPS.
//##############################################################################
//------------------------------------------------------------------------------
/// Computes sub-particle stress tensor (Tau) for SPS turbulence model.
//------------------------------------------------------------------------------
__global__ void KerComputeSpsTau(unsigned n,unsigned pini,float smag,float blin
  ,const float4* velrho,const float2* sps2strain,float2* tau_rho2)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; 
  if(p<n){
    const unsigned p1=p+pini;
    float2 rr=sps2strain[p1*3];   const float two_strain_xx=rr.x,two_strain_xy=rr.y;
           rr=sps2strain[p1*3+1]; const float two_strain_xz=rr.x,two_strain_yy=rr.y;
           rr=sps2strain[p1*3+2]; const float two_strain_yz=rr.x,two_strain_zz=rr.y;
    const float pow1=two_strain_xx*two_strain_xx
                   + two_strain_yy*two_strain_yy
                   + two_strain_zz*two_strain_zz;
    const float prr= two_strain_xy*two_strain_xy
                   + two_strain_xz*two_strain_xz
                   + two_strain_yz*two_strain_yz + pow1+pow1;
    const float visc_sps=smag*sqrt(prr);
    const float div_u=two_strain_xx+two_strain_yy+two_strain_zz;
    const float sps_k=visc_sps*div_u; //-Factor 2/3 is included in smag constant.
    const float sps_blin=blin*prr;
    const float sumsps=-(sps_k+sps_blin);
    const float twovisc_sps=(visc_sps+visc_sps);
    float one_rho=1.0f/velrho[p1].w;
    //-Computes new values of tau/rho^2.
    const float tau_xx=one_rho*(twovisc_sps*two_strain_xx +sumsps);
    const float tau_xy=one_rho*(visc_sps   *two_strain_xy);
    tau_rho2[p1*3]=make_float2(tau_xx,tau_xy);
    const float tau_xz=one_rho*(visc_sps   *two_strain_xz);
    const float tau_yy=one_rho*(twovisc_sps*two_strain_yy +sumsps);
    tau_rho2[p1*3+1]=make_float2(tau_xz,tau_yy);
    const float tau_yz=one_rho*(visc_sps   *two_strain_yz);
    const float tau_zz=one_rho*(twovisc_sps*two_strain_zz +sumsps);
    tau_rho2[p1*3+2]=make_float2(tau_yz,tau_zz);
  }
}

//==============================================================================
/// Computes sub-particle stress tensor divided by rho^2 (tau/rho^2) for SPS 
/// turbulence model.   
//==============================================================================
void ComputeSpsTau(unsigned np,unsigned npb,float smag,float blin
  ,const float4* velrho,const tsymatrix3f* sps2strain,tsymatrix3f* tau_rho2
  ,cudaStream_t stm)
{
  const unsigned npf=np-npb;
  if(npf){
    dim3 sgridf=GetSimpleGridSize(npf,SPHBSIZE);
    KerComputeSpsTau <<<sgridf,SPHBSIZE,0,stm>>> (npf,npb,smag,blin,velrho
      ,(const float2*)sps2strain,(float2*)tau_rho2);
  }
}


//##############################################################################
//# Kernels for Delta-SPH.
//# Kernels para Delta-SPH.
//##############################################################################
//------------------------------------------------------------------------------
/// Adds value of delta[] to ar[] provided it is not FLT_MAX.
/// Anhade valor de delta[] a ar[] siempre que no sea FLT_MAX.
//------------------------------------------------------------------------------
__global__ void KerAddDelta(unsigned n,const float* delta,float* ar)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    float rdelta=delta[p];
    if(rdelta!=FLT_MAX)ar[p]+=rdelta;
  }
}

//==============================================================================
/// Adds value of delta[] to ar[] provided it is not FLT_MAX.
/// Anhade valor de delta[] a ar[] siempre que no sea FLT_MAX.
//==============================================================================
void AddDelta(unsigned n,const float* delta,float* ar,cudaStream_t stm){
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerAddDelta <<<sgrid,SPHBSIZE,0,stm>>> (n,delta,ar);
  }
}


//##############################################################################
//# Kernels para ComputeStep (position)
//# Kernels for ComputeStep (position)
//##############################################################################
//------------------------------------------------------------------------------
/// Updates particle position according to displacement.
/// Actualizacion de posicion de particulas segun desplazamiento.
//------------------------------------------------------------------------------
template<bool periactive,bool floatings> __global__ void KerComputeStepPos(
  unsigned n,unsigned pini,const double2* movxy,const double* movz
  ,double2* posxy,double* posz,unsigned* dcell,typecode* code)
{
  unsigned pt=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pt<n){
    unsigned p=pt+pini;
    const typecode rcode=code[p];
    const bool outrhop=CODE_IsOutRho(rcode);
    const bool fluid=(!floatings || CODE_IsFluid(rcode));
    const bool normal=(!periactive || outrhop || CODE_IsNormal(rcode));
    if(normal && fluid){ //-Does not apply to periodic or floating particles. | No se aplica a particulas periodicas o floating.
      const double2 rmovxy=movxy[p];
      KerUpdatePos<periactive>(posxy[p],posz[p],rmovxy.x,rmovxy.y,movz[p],outrhop,p,posxy,posz,dcell,code);
    }
    //-In case of floating maintains the original position.
    //-En caso de floating mantiene la posicion original.
  }
}

//==============================================================================
/// Updates particle position according to displacement.
/// Actualizacion de posicion de particulas segun desplazamiento.
//==============================================================================
void ComputeStepPos(byte periactive,bool floatings,unsigned np,unsigned npb
  ,const double2* movxy,const double* movz
  ,double2* posxy,double* posz,unsigned* dcell,typecode* code)
{
  const unsigned pini=npb;
  const unsigned npf=np-pini;
  if(npf){
    dim3 sgrid=GetSimpleGridSize(npf,SPHBSIZE);
    if(periactive){ const bool peri=true;
      if(floatings)KerComputeStepPos<peri,true>  <<<sgrid,SPHBSIZE>>> (npf,pini,movxy,movz,posxy,posz,dcell,code);
      else         KerComputeStepPos<peri,false> <<<sgrid,SPHBSIZE>>> (npf,pini,movxy,movz,posxy,posz,dcell,code);
    }
    else{ const bool peri=false;
      if(floatings)KerComputeStepPos<peri,true>  <<<sgrid,SPHBSIZE>>> (npf,pini,movxy,movz,posxy,posz,dcell,code);
      else         KerComputeStepPos<peri,false> <<<sgrid,SPHBSIZE>>> (npf,pini,movxy,movz,posxy,posz,dcell,code);
    }
  }
}

//------------------------------------------------------------------------------
/// Updates particle position according to displacement.
/// Actualizacion de posicion de particulas segun desplazamiento.
//------------------------------------------------------------------------------
template<bool periactive,bool floatings> __global__ void KerComputeStepPos2(
  unsigned n,unsigned pini,const double2* posxypre,const double* poszpre
  ,const double2* movxy,const double* movz,double2* posxy,double* posz
  ,unsigned* dcell,typecode* code)
{
  unsigned pt=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pt<n){
    unsigned p=pt+pini;
    const typecode rcode=code[p];
    const bool outrhop=CODE_IsOutRho(rcode);
    const bool fluid=(!floatings || CODE_IsFluid(rcode));
    const bool normal=(!periactive || outrhop || CODE_IsNormal(rcode));
    if(normal){//-Does not apply to periodic particles. | No se aplica a particulas periodicas
      if(fluid){//-Only applied for fluid displacement. | Solo se aplica desplazamiento al fluido.
        const double2 rmovxy=movxy[p];
        KerUpdatePos<periactive>(posxypre[p],poszpre[p],rmovxy.x,rmovxy.y,movz[p],outrhop,p,posxy,posz,dcell,code);
      }
      else{ //-Copy position of floating particles.
        posxy[p]=posxypre[p];
        posz[p]=poszpre[p];
      }
    }
  }
}

//==============================================================================
/// Updates particle position according to displacement.
/// Actualizacion de posicion de particulas segun desplazamiento.
//==============================================================================
void ComputeStepPos2(byte periactive,bool floatings,unsigned np,unsigned npb
  ,const double2* posxypre,const double* poszpre,const double2* movxy
  ,const double* movz,double2* posxy,double* posz,unsigned* dcell,typecode* code)
{
  const unsigned pini=npb;
  const unsigned npf=np-pini;
  if(npf){
    dim3 sgrid=GetSimpleGridSize(npf,SPHBSIZE);
    if(periactive){ const bool peri=true;
      if(floatings)KerComputeStepPos2<peri,true>  <<<sgrid,SPHBSIZE>>> (npf,pini,posxypre,poszpre,movxy,movz,posxy,posz,dcell,code);
      else         KerComputeStepPos2<peri,false> <<<sgrid,SPHBSIZE>>> (npf,pini,posxypre,poszpre,movxy,movz,posxy,posz,dcell,code);
    }
    else{ const bool peri=false;
      if(floatings)KerComputeStepPos2<peri,true>  <<<sgrid,SPHBSIZE>>> (npf,pini,posxypre,poszpre,movxy,movz,posxy,posz,dcell,code);
      else         KerComputeStepPos2<peri,false> <<<sgrid,SPHBSIZE>>> (npf,pini,posxypre,poszpre,movxy,movz,posxy,posz,dcell,code);
    }
  }
}



//##############################################################################
//# Kernels for motion.
//# Kernels para Motion
//##############################################################################
//------------------------------------------------------------------------------
/// Computes for a range of particles, their position according to idp[].
/// Calcula para un rango de particulas calcula su posicion segun idp[].
//------------------------------------------------------------------------------
__global__ void KerCalcRidp(unsigned n,unsigned ini,unsigned idini,unsigned idfin
  ,const typecode* code,const unsigned* idp,unsigned* ridp)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    p+=ini;
    unsigned id=idp[p];
    if(idini<=id && id<idfin){
      if(CODE_IsNormal(code[p]))ridp[id-idini]=p;
    }
  }
}
//------------------------------------------------------------------------------
__global__ void KerCalcRidp(unsigned n,unsigned ini,unsigned idini,unsigned idfin
  ,const unsigned* idp,unsigned* ridp)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    p+=ini;
    const unsigned id=idp[p];
    if(idini<=id && id<idfin)ridp[id-idini]=p;
  }
}

//==============================================================================
/// Calculate particle position according to idp[]. When it does not find UINT_MAX.
/// When periactive is false it means there are no duplicate particles (periodic)
/// and all are CODE_NORMAL.
///
/// Calcula posicion de particulas segun idp[]. Cuando no la encuentra es UINT_MAX.
/// Cuando periactive es False sumpone que no hay particulas duplicadas (periodicas)
/// y todas son CODE_NORMAL.
//==============================================================================
void CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini
  ,unsigned idfin,const typecode* code,const unsigned* idp,unsigned* ridp
  ,cudaStream_t stm)
{
  //-Assigns values UINT_MAX
  const unsigned nsel=idfin-idini;
  cudaMemset(ridp,255,sizeof(unsigned)*nsel); 
  //-Computes position according to id. | Calcula posicion segun id.
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
    if(periactive)KerCalcRidp <<<sgrid,SPHBSIZE,0,stm>>> (np,pini,idini,idfin,code,idp,ridp);
    else          KerCalcRidp <<<sgrid,SPHBSIZE,0,stm>>> (np,pini,idini,idfin,idp,ridp);
  }
}


//------------------------------------------------------------------------------
/// Load current reference position data from particle data.
//------------------------------------------------------------------------------
__global__ void KerLoadPosRef(unsigned pscount,unsigned casenfixed,unsigned np
  ,const double2* posxy,const double* posz,const unsigned* ridpmot
  ,const unsigned* idpref,double3* posref)
{
  unsigned cp=blockIdx.x*blockDim.x + threadIdx.x;
  if(cp<pscount){
    const unsigned iref=idpref[cp]-casenfixed;
    const unsigned p=ridpmot[iref];
    if(p<np){
      const double2 rxy=posxy[p];
      const double rz=posz[p];
      posref[cp]=make_double3(rxy.x,rxy.y,rz);
    }
    else posref[cp]=make_double3(DBL_MAX,DBL_MAX,DBL_MAX);
  }
}

//==============================================================================
/// Load current reference position data from particle data.
//==============================================================================
void LoadPosRef(unsigned pscount,unsigned casenfixed,unsigned np
  ,const double2* posxy,const double* posz,const unsigned* ridpmot
  ,const unsigned* idpref,double3* posref)
{
  //-Computes position according to id. | Calcula posicion segun id.
  if(pscount){
    const dim3 sgrid=GetSimpleGridSize(pscount,SPHBSIZE);
    KerLoadPosRef <<<sgrid,SPHBSIZE>>> (pscount,casenfixed,np,posxy,posz
      ,ridpmot,idpref,posref);
  }
}


//------------------------------------------------------------------------------
/// Applies a linear movement to a set of particles.
/// Aplica un movimiento lineal a un conjunto de particulas.
//------------------------------------------------------------------------------
template<bool periactive> __global__ void KerMoveLinBound(unsigned n,unsigned ini
  ,double3 mvpos,float3 mvvel,const unsigned* ridpmot,double2* posxy,double* posz
  ,unsigned* dcell,float4* velrho,typecode* code)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    int pid=ridpmot[p+ini];
    if(pid>=0){
      //-Computes displacement and updates position.
      KerUpdatePos<periactive>(posxy[pid],posz[pid],mvpos.x,mvpos.y,mvpos.z,false,pid,posxy,posz,dcell,code);
      //-Computes velocity.
      velrho[pid]=make_float4(mvvel.x,mvvel.y,mvvel.z,velrho[pid].w);
    }
  }
}

//==============================================================================
/// Applies a linear movement to a set of particles.
/// Aplica un movimiento lineal a un conjunto de particulas.
//==============================================================================
void MoveLinBound(byte periactive,unsigned np,unsigned ini,tdouble3 mvpos
  ,tfloat3 mvvel,const unsigned* ridpmot,double2* posxy,double* posz
  ,unsigned* dcell,float4* velrho,typecode* code)
{
  dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
  if(periactive)KerMoveLinBound<true>  <<<sgrid,SPHBSIZE>>> (np,ini,Double3(mvpos),Float3(mvvel),ridpmot,posxy,posz,dcell,velrho,code);
  else          KerMoveLinBound<false> <<<sgrid,SPHBSIZE>>> (np,ini,Double3(mvpos),Float3(mvvel),ridpmot,posxy,posz,dcell,velrho,code);
}



//------------------------------------------------------------------------------
/// Applies a matrix movement to a set of particles.
/// Aplica un movimiento matricial a un conjunto de particulas.
//------------------------------------------------------------------------------
template<bool periactive,bool simulate2d> __global__ void KerMoveMatBound(
  unsigned n,unsigned ini,tmatrix4d m,double dt,const unsigned* ridmot
  ,double2* posxy,double* posz,unsigned* dcell,float4* velrho,typecode* code
  ,float3* boundnor)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    int pid=ridmot[p+ini];
    if(pid>=0){
      double2 rxy=posxy[pid];
      double3 rpos=make_double3(rxy.x,rxy.y,posz[pid]);
      //-Computes new position.
      double3 rpos2;
      rpos2.x= rpos.x*m.a11 + rpos.y*m.a12 + rpos.z*m.a13 + m.a14;
      rpos2.y= rpos.x*m.a21 + rpos.y*m.a22 + rpos.z*m.a23 + m.a24;
      rpos2.z= rpos.x*m.a31 + rpos.y*m.a32 + rpos.z*m.a33 + m.a34;
      if(simulate2d)rpos2.y=rpos.y;
      //-Computes displacement and updates position.
      const double dx=rpos2.x-rpos.x;
      const double dy=rpos2.y-rpos.y;
      const double dz=rpos2.z-rpos.z;
      KerUpdatePos<periactive>(make_double2(rpos.x,rpos.y),rpos.z,dx,dy,dz,false,pid,posxy,posz,dcell,code);
      //-Computes velocity.
      velrho[pid]=make_float4(float(dx/dt),float(dy/dt),float(dz/dt),velrho[pid].w);
      //-Computes normal.
      if(boundnor){
        const float3 bnor=boundnor[pid];
        const double3 gs=make_double3(rpos.x+bnor.x,rpos.y+bnor.y,rpos.z+bnor.z);
        const double gs2x=gs.x*m.a11 + gs.y*m.a12 + gs.z*m.a13 + m.a14;
        const double gs2y=gs.x*m.a21 + gs.y*m.a22 + gs.z*m.a23 + m.a24;
        const double gs2z=gs.x*m.a31 + gs.y*m.a32 + gs.z*m.a33 + m.a34;
        boundnor[pid]=make_float3(gs2x-rpos2.x,gs2y-rpos2.y,gs2z-rpos2.z);
      }
    }
  }
}

//==============================================================================
/// Applies a matrix movement to a set of particles.
/// Aplica un movimiento matricial a un conjunto de particulas.
//==============================================================================
void MoveMatBound(byte periactive,bool simulate2d,unsigned np,unsigned ini
  ,tmatrix4d m,double dt,const unsigned* ridpmot,double2* posxy,double* posz
  ,unsigned* dcell,float4* velrho,typecode* code,float3* boundnor)
{
  dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
  if(periactive){ const bool peri=true;
    if(simulate2d)KerMoveMatBound<peri,true>  <<<sgrid,SPHBSIZE>>> (np,ini,m,dt,ridpmot,posxy,posz,dcell,velrho,code,boundnor);
    else          KerMoveMatBound<peri,false> <<<sgrid,SPHBSIZE>>> (np,ini,m,dt,ridpmot,posxy,posz,dcell,velrho,code,boundnor);
  }
  else{ const bool peri=false;
    if(simulate2d)KerMoveMatBound<peri,true>  <<<sgrid,SPHBSIZE>>> (np,ini,m,dt,ridpmot,posxy,posz,dcell,velrho,code,boundnor);
    else          KerMoveMatBound<peri,false> <<<sgrid,SPHBSIZE>>> (np,ini,m,dt,ridpmot,posxy,posz,dcell,velrho,code,boundnor);
  }
}




//------------------------------------------------------------------------------
/// Copy motion velocity to MotionAce[].
/// Copia velocidad de movimiento a MotionAce[].
//------------------------------------------------------------------------------
template<bool periactive> __global__ void KerCopyMotionAce(unsigned n
    , const unsigned* ridpmv, const float4* velrhop, float3* motionvel, float3* motionace, double stepdt)
{
    unsigned p = blockIdx.x * blockDim.x + threadIdx.x; //-Number of particle.
    float dt = float(stepdt);
    if (p < n) {
        int pid = ridpmv[p];
        if (pid >= 0) {
            //-Computes velocity.
            const float4 v = velrhop[pid];
            motionace[pid] = make_float3(float((v.x - motionvel[pid].x) / dt), float((v.y - motionvel[pid].y) / dt), float((v.z - motionvel[pid].z) / dt));
        }
    }
}

//==============================================================================
/// Give boundary acceleration to MotionAce[].
/// SHABA
//==============================================================================
void CopyMotionAce(unsigned nmoving, const unsigned* ridp, const float4* velrhop, float3* motionvel, float3* motionace, double stepdt)
{
    dim3 sgrid = GetSimpleGridSize(nmoving, SPHBSIZE);
    KerCopyMotionAce<true> << <sgrid, SPHBSIZE >> > (nmoving, ridp, velrhop, motionvel,motionace,stepdt);
}

//------------------------------------------------------------------------------
/// Copy motion velocity to MotionVel[].
/// Copia velocidad de movimiento a MotionVel[].
//------------------------------------------------------------------------------
template<bool periactive> __global__ void KerCopyMotionVel(unsigned n
  ,const unsigned* ridpmot,const float4* velrho,float3* motionvel)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    int pid=ridpmot[p];
    if(pid>=0){
      //-Computes velocity.
      const float4 v=velrho[pid];
      motionvel[pid]=make_float3(v.x,v.y,v.z);
    }
  }
}

//==============================================================================
/// Copy motion velocity to MotionVel[].
/// Copia velocidad de movimiento a MotionVel[].
//==============================================================================
void CopyMotionVel(unsigned nmoving,const unsigned* ridpmot
  ,const float4* velrho,float3* motionvel)
{
  dim3 sgrid=GetSimpleGridSize(nmoving,SPHBSIZE);
  KerCopyMotionVel<true>  <<<sgrid,SPHBSIZE>>> (nmoving,ridpmot,velrho,motionvel);
}


//------------------------------------------------------------------------------
/// Applies a matrix movement to a set of particles.
/// Aplica un movimiento matricial a un conjunto de particulas.
//------------------------------------------------------------------------------
__global__ void KerFtNormalsUpdate(unsigned n,unsigned fpini
  ,double a11,double a12,double a13,double a21,double a22,double a23
  ,double a31,double a32,double a33
  ,const unsigned* ridpmot,float3* boundnor)
{
  const unsigned fp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of floating particle.
  if(fp<n){
    const unsigned p=ridpmot[fp+fpini];
    if(p!=UINT_MAX){
      float3 rnor=boundnor[p];
      const double nx=rnor.x;
      const double ny=rnor.y;
      const double nz=rnor.z;
      rnor.x=float(a11*nx + a12*ny + a13*nz);
      rnor.y=float(a21*nx + a22*ny + a23*nz);
      rnor.z=float(a31*nx + a32*ny + a33*nz);
      boundnor[p]=rnor;
    }
  }
}

//==============================================================================
/// Applies a matrix movement to a set of particles.
/// Aplica un movimiento matricial a un conjunto de particulas.
//==============================================================================
void FtNormalsUpdate(unsigned np,unsigned ini,tmatrix4d m,const unsigned* ridpmot
  ,float3* boundnor)
{
  dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
  if(np)KerFtNormalsUpdate <<<sgrid,SPHBSIZE>>> (np,ini,m.a11,m.a12,m.a13
    ,m.a21,m.a22,m.a23,m.a31,m.a32,m.a33,ridpmot,boundnor);
}



//##############################################################################
//# Kernels for MLPistons motion.
//##############################################################################
//------------------------------------------------------------------------------
/// Applies movement and velocity of piston 1D to a group of particles.
/// Aplica movimiento y velocidad de piston 1D a conjunto de particulas.
//------------------------------------------------------------------------------
template<byte periactive> __global__ void KerMovePiston1d(unsigned n,unsigned idini
  ,double dp,double poszmin,unsigned poszcount,const byte* pistonid
  ,const double* movx,const double* velx,const unsigned* ridpmv,double2* posxy
  ,double* posz,unsigned* dcell,float4* velrho,typecode* code)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle
  if(p<n){
    const unsigned id=p+idini;
    int pid=ridpmv[id];
    if(pid>=0){
      const unsigned pisid=pistonid[CODE_GetTypeValue(code[pid])];
      if(pisid<255){
        const double2 rpxy=posxy[pid];
        const double rpz=posz[pid];
        const unsigned cz=unsigned((rpz-poszmin)/dp);
        const double rmovx=(cz<poszcount? movx[pisid*poszcount+cz]: 0);
        const float rvelx=float(cz<poszcount? velx[pisid*poszcount+cz]: 0);
        //-Updates position.
        KerUpdatePos<periactive>(rpxy,rpz,rmovx,0,0,false,pid,posxy,posz,dcell,code);
        //-Updates velocity.
        velrho[pid].x=rvelx;
      }
    }
  }
}

//==============================================================================
/// Applies movement and velocity of piston 1D to a group of particles.
/// Aplica movimiento y velocidad de piston 1D a conjunto de particulas.
//==============================================================================
void MovePiston1d(bool periactive,unsigned np,unsigned idini
  ,double dp,double poszmin,unsigned poszcount,const byte* pistonid
  ,const double* movx,const double* velx,const unsigned* ridpmv,double2* posxy
  ,double* posz,unsigned* dcell,float4* velrho,typecode* code)
{
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
    if(periactive)KerMovePiston1d<true>  <<<sgrid,SPHBSIZE>>> (np,idini,dp,poszmin,poszcount,pistonid,movx,velx,ridpmv,posxy,posz,dcell,velrho,code);
    else          KerMovePiston1d<false> <<<sgrid,SPHBSIZE>>> (np,idini,dp,poszmin,poszcount,pistonid,movx,velx,ridpmv,posxy,posz,dcell,velrho,code);
  }
}

//------------------------------------------------------------------------------
/// Applies movement and velocity of piston 2D to a group of particles.
/// Aplica movimiento y velocidad de piston 2D a conjunto de particulas.
//------------------------------------------------------------------------------
template<byte periactive> __global__ void KerMovePiston2d(unsigned n,unsigned idini
  ,double dp,double posymin,double poszmin,unsigned poszcount,const double* movx
  ,const double* velx,const unsigned* ridpmot,double2* posxy,double* posz
  ,unsigned* dcell,float4* velrho,typecode* code)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle
  if(p<n){
    const unsigned id=p+idini;
    int pid=ridpmot[id];
    if(pid>=0){
      const double2 rpxy=posxy[pid];
      const double rpz=posz[pid];
      const unsigned cy=unsigned((rpxy.y-posymin)/dp);
      const unsigned cz=unsigned((rpz-poszmin)/dp);
      const double rmovx=(cz<poszcount? movx[cy*poszcount+cz]: 0);
      const float rvelx=float(cz<poszcount? velx[cy*poszcount+cz]: 0);
      //-Actualiza posicion.
      KerUpdatePos<periactive>(rpxy,rpz,rmovx,0,0,false,pid,posxy,posz,dcell,code);
      //-Actualiza velocidad.
      velrho[pid].x=rvelx;
    }
  }
}

//==============================================================================
/// Applies movement and velocity of piston 2D to a group of particles.
/// Aplica movimiento y velocidad de piston 2D a conjunto de particulas.
//==============================================================================
void MovePiston2d(bool periactive,unsigned np,unsigned idini,double dp
  ,double posymin,double poszmin,unsigned poszcount,const double* movx
  ,const double* velx,const unsigned* ridpmot,double2* posxy,double* posz
  ,unsigned* dcell,float4* velrho,typecode* code)
{
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
    if(periactive)KerMovePiston2d<true>  <<<sgrid,SPHBSIZE>>> (np,idini,dp,posymin,poszmin,poszcount,movx,velx,ridpmot,posxy,posz,dcell,velrho,code);
    else          KerMovePiston2d<false> <<<sgrid,SPHBSIZE>>> (np,idini,dp,posymin,poszmin,poszcount,movx,velx,ridpmot,posxy,posz,dcell,velrho,code);
  }
}


//##############################################################################
//# Kernels for Floating bodies.
//##############################################################################
//==============================================================================
/// Computes distance between floating and centre particles according to periodic conditions.
/// Calcula distancia entre pariculas floating y centro segun condiciones periodicas.
//==============================================================================
template<bool periactive> __device__ void KerFtPeriodicDist(double px,double py
  ,double pz,double cenx,double ceny,double cenz,float radius,float& dx
  ,float& dy,float& dz)
{
  if(periactive){
    double ddx=px-cenx;
    double ddy=py-ceny;
    double ddz=pz-cenz;
    const unsigned peri=CTE.periactive;
    if(PERI_AxisX(peri) && fabs(ddx)>radius){
      if(ddx>0){ ddx+=CTE.xperincx; ddy+=CTE.xperincy; ddz+=CTE.xperincz; }
      else{      ddx-=CTE.xperincx; ddy-=CTE.xperincy; ddz-=CTE.xperincz; }
    }
    if(PERI_AxisY(peri) && fabs(ddy)>radius){
      if(ddy>0){ ddx+=CTE.yperincx; ddy+=CTE.yperincy; ddz+=CTE.yperincz; }
      else{      ddx-=CTE.yperincx; ddy-=CTE.yperincy; ddz-=CTE.yperincz; }
    }
    if(PERI_AxisZ(peri) && fabs(ddz)>radius){
      if(ddz>0){ ddx+=CTE.zperincx; ddy+=CTE.zperincy; ddz+=CTE.zperincz; }
      else{      ddx-=CTE.zperincx; ddy-=CTE.zperincy; ddz-=CTE.zperincz; }
    }
    dx=float(ddx);
    dy=float(ddy);
    dz=float(ddz);
  }
  else{
    dx=float(px-cenx);
    dy=float(py-ceny);
    dz=float(pz-cenz);
  }
}

//------------------------------------------------------------------------------
/// Calculate summation: face, fomegaace in ftoforcessum[].
/// Calcula suma de face y fomegaace a partir de particulas floating en ftoforcessum[].
//------------------------------------------------------------------------------
template<bool periactive> __global__ void KerFtPartsSumAce( //ftodatp={pini-CaseNfixed,np,radius,massp}
  const float4* ftodatp,const double3* ftocenter,const unsigned* ridpmot
  ,const double2* posxy,const double* posz,const float3* ace
  ,float3* ftoacelinang)
{
  extern __shared__ float racelinx[];
  float* raceliny=racelinx+blockDim.x;
  float* racelinz=raceliny+blockDim.x;
  float* raceangx=racelinz+blockDim.x;
  float* raceangy=raceangx+blockDim.x;
  float* raceangz=raceangy+blockDim.x;

  const unsigned tid=threadIdx.x;  //-Thread number.
  const unsigned cf=blockIdx.x;    //-Floating number.
  
  //-Loads floating data.
  const float4 rfdata=ftodatp[cf];
  const unsigned fpini=(unsigned)__float_as_int(rfdata.x);
  const unsigned fnp=(unsigned)__float_as_int(rfdata.y);
  const float fradius=rfdata.z;
  //const float fmassp=rfdata.w;
  const double3 rcenter=ftocenter[cf];

  //-Initialises shared memory to zero.
  const unsigned ntid=(fnp<blockDim.x? fnp: blockDim.x); //-Number of used threads. | Numero de threads utilizados.
  if(tid<ntid){
    racelinx[tid]=raceliny[tid]=racelinz[tid]=0;
    raceangx[tid]=raceangy[tid]=raceangz[tid]=0;
  }

  //-Computes data in shared memory. | Calcula datos en memoria shared.
  const unsigned nfor=unsigned((fnp+blockDim.x-1)/blockDim.x);
  for(unsigned cfor=0;cfor<nfor;cfor++){
    unsigned p=cfor*blockDim.x+tid;
    if(p<fnp){
      const unsigned rp=ridpmot[p+fpini];
      if(rp!=UINT_MAX){
        const float3 acep=ace[rp];
        racelinx[tid]+=acep.x; raceliny[tid]+=acep.y; racelinz[tid]+=acep.z;
        //-Computes distance from the centre. | Calcula distancia al centro.
        const double2 rposxy=posxy[rp];
        float dx,dy,dz;
        KerFtPeriodicDist<periactive>(rposxy.x,rposxy.y,posz[rp],rcenter.x,rcenter.y,rcenter.z,fradius,dx,dy,dz);
        //-Computes omegaace.
        raceangx[tid]+=(acep.z*dy - acep.y*dz);
        raceangy[tid]+=(acep.x*dz - acep.z*dx);
        raceangz[tid]+=(acep.y*dx - acep.x*dy);
      }
    }
  }

  //-Reduces data in shared memory and stores results.
  //-Reduce datos de memoria shared y guarda resultados.
  __syncthreads();
  if(!tid){
    float3 acelin=make_float3(0,0,0);
    float3 aceang=make_float3(0,0,0);
    for(unsigned c=0;c<ntid;c++){
      acelin.x+=racelinx[c];  acelin.y+=raceliny[c];  acelin.z+=racelinz[c];
      aceang.x+=raceangx[c];  aceang.y+=raceangy[c];  aceang.z+=raceangz[c];
    }
    //-Stores results in ftoacelinang[].
    const unsigned cf2=cf*2;
    ftoacelinang[cf2  ]=acelin;
    ftoacelinang[cf2+1]=aceang;
  }
}

//==============================================================================
/// Calculate summation: face, fomegaace in ftoforcessum[].
/// Calcula suma de face y fomegaace a partir de particulas floating en ftoforcessum[].
//==============================================================================
void FtPartsSumAce(bool periactive,unsigned ftcount
  ,const float4* ftodatp,const double3* ftocenter,const unsigned* ridpmot
  ,const double2* posxy,const double* posz,const float3* ace
  ,float3* ftoacelinang)
{
  if(ftcount){
    const unsigned bsize=256;
    const unsigned smem=sizeof(float)*(3+3)*bsize;
    dim3 sgrid=GetSimpleGridSize(ftcount*bsize,bsize);
    if(periactive)KerFtPartsSumAce<true>  <<<sgrid,bsize,smem>>> (ftodatp,ftocenter,ridpmot,posxy,posz,ace,ftoacelinang);
    else          KerFtPartsSumAce<false> <<<sgrid,bsize,smem>>> (ftodatp,ftocenter,ridpmot,posxy,posz,ace,ftoacelinang);
  }
}

//------------------------------------------------------------------------------
/// Updates information and particles of floating bodies.
//------------------------------------------------------------------------------
template<bool periactive> __global__ void KerFtPartsUpdate(double dt
  ,bool updatenormals,unsigned np,unsigned fpini,float fradius,tmatrix4d mat
  ,float3 fvel,float3 fomega,double3 fcenter
  ,const unsigned* ridpmot,double2* posxy,double* posz,float4* velrho
  ,unsigned* dcell,typecode* code,float3* boundnor)
{
  const unsigned fp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(fp<np){
    const int p=ridpmot[fp+fpini];
    if(p>=0){
      float4 vr=velrho[p];
      //-Computes displacement and updates position.
      const double dx=dt*double(vr.x);
      const double dy=dt*double(vr.y);
      const double dz=dt*double(vr.z);
      const double2 rxy=posxy[p];
      KerUpdatePos<periactive>(rxy,posz[p],dx,dy,dz,false,p,posxy,posz,dcell,code);
      //-Computes and updates velocity.
      const double2 rposxy=posxy[p];
      float distx,disty,distz;
      KerFtPeriodicDist<periactive>(rposxy.x,rposxy.y,posz[p],fcenter.x,fcenter.y,fcenter.z,fradius,distx,disty,distz);
      vr.x=fvel.x+(fomega.y*distz - fomega.z*disty);
      vr.y=fvel.y+(fomega.z*distx - fomega.x*distz);
      vr.z=fvel.z+(fomega.x*disty - fomega.y*distx);
      velrho[p]=vr;
      //-Updates floating normals for mDBC.
      if(updatenormals){
        const float3 norf=boundnor[p];
        const double norx=norf.x;
        const double nory=norf.y;
        const double norz=norf.z;
        const float nx=float( mat.a11*norx + mat.a12*nory + mat.a13*norz );
        const float ny=float( mat.a21*norx + mat.a22*nory + mat.a23*norz );
        const float nz=float( mat.a31*norx + mat.a32*nory + mat.a33*norz );
        boundnor[p]=make_float3(nx,ny,nz);
      }
    }
  }
}

//==============================================================================
/// Updates information and particles of floating bodies.
//==============================================================================
void FtPartsUpdate(bool periactive,double dt,bool updatenormals
  ,unsigned np,unsigned fpini,float fradius,tmatrix4d mat
  ,tfloat3 fto_vellin,tfloat3 fto_velang,tdouble3 fto_center
  ,const unsigned* ridpmot,double2* posxy,double* posz,float4* velrho
  ,unsigned* dcell,typecode* code,float3* boundnor,cudaStream_t stm)
{
  if(np){
    const unsigned bsize=128; 
    dim3 sgrid=GetSimpleGridSize(np,bsize);
    if(periactive)KerFtPartsUpdate<true>  <<<sgrid,bsize,0,stm>>> (dt,updatenormals,np,fpini,fradius,mat,Float3(fto_vellin),Float3(fto_velang),Double3(fto_center),ridpmot,posxy,posz,velrho,dcell,code,boundnor);
    else          KerFtPartsUpdate<false> <<<sgrid,bsize,0,stm>>> (dt,updatenormals,np,fpini,fradius,mat,Float3(fto_vellin),Float3(fto_velang),Double3(fto_center),ridpmot,posxy,posz,velrho,dcell,code,boundnor); 
  }
}



//##############################################################################
//# Kernels for Periodic conditions
//# Kernels para Periodic conditions
//##############################################################################
//------------------------------------------------------------------------------
/// Marks current periodics to be ignored.
/// Marca las periodicas actuales como ignorar.
//------------------------------------------------------------------------------
__global__ void KerPeriodicIgnore(unsigned n,typecode* code)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    //-Checks code of particles.
    //-Comprueba codigo de particula.
    const typecode rcode=code[p];
    if(CODE_IsPeriodic(rcode))code[p]=CODE_SetOutIgnore(rcode);
  }
}

//==============================================================================
/// Marks current periodics to be ignored.
/// Marca las periodicas actuales como ignorar.
//==============================================================================
void PeriodicIgnore(unsigned n,typecode* code){
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerPeriodicIgnore <<<sgrid,SPHBSIZE>>> (n,code);
  }
}

//------------------------------------------------------------------------------
/// Create list of new periodic particles to be duplicated and 
/// marks old periodics to be ignored.
///
/// Crea lista de nuevas particulas periodicas a duplicar y con delper activado
/// marca las periodicas viejas para ignorar.
//------------------------------------------------------------------------------
__global__ void KerPeriodicMakeList(unsigned n,unsigned pini,unsigned nmax
  ,double3 mapposmin,double3 mapposmax,double3 perinc
  ,const double2* posxy,const double* posz,const typecode* code,unsigned* listp)
{
  extern __shared__ unsigned slist[];
  if(!threadIdx.x)slist[0]=0;
  __syncthreads();
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p2=p+pini;
    //-Inteacts with normal or periodic particles.
    //-Se queda con particulas normales o periodicas.
    if(CODE_GetSpecialValue(code[p2])<=CODE_PERIODIC){
      //-Obtains particle position.
      const double2 rxy=posxy[p2];
      const double rx=rxy.x,ry=rxy.y;
      const double rz=posz[p2];
      double rx2=rx+perinc.x,ry2=ry+perinc.y,rz2=rz+perinc.z;
      if(mapposmin.x<=rx2 && mapposmin.y<=ry2 && mapposmin.z<=rz2 && rx2<mapposmax.x && ry2<mapposmax.y && rz2<mapposmax.z){
        unsigned cp=atomicAdd(slist,1);  slist[cp+1]=p2;
      }
      rx2=rx-perinc.x; ry2=ry-perinc.y; rz2=rz-perinc.z;
      if(mapposmin.x<=rx2 && mapposmin.y<=ry2 && mapposmin.z<=rz2 && rx2<mapposmax.x && ry2<mapposmax.y && rz2<mapposmax.z){
        unsigned cp=atomicAdd(slist,1);  slist[cp+1]=(p2|0x80000000);
      }
    }
  }
  __syncthreads();
  const unsigned ns=slist[0];
  __syncthreads();
  if(!threadIdx.x && ns)slist[0]=atomicAdd((listp+nmax),ns);
  __syncthreads();
  if(threadIdx.x<ns){
    unsigned cp=slist[0]+threadIdx.x;
    if(cp<nmax)listp[cp]=slist[threadIdx.x+1];
  }
  if(blockDim.x+threadIdx.x<ns){ //-There may be twice as many periodics per thread. | Puede haber el doble de periodicas que threads.
    unsigned cp=blockDim.x+slist[0]+threadIdx.x;
    if(cp<nmax)listp[cp]=slist[blockDim.x+threadIdx.x+1];
  }
}

//==============================================================================
/// Create list of new periodic particles to be duplicated.
/// With stable activated reorders perioc list.
///
/// Crea lista de nuevas particulas periodicas a duplicar.
/// Con stable activado reordena lista de periodicas.
//==============================================================================
unsigned PeriodicMakeList(unsigned n,unsigned pini,bool stable,unsigned nmax
  ,tdouble3 mapposmin,tdouble3 mapposmax,tdouble3 perinc
  ,const double2* posxy,const double* posz,const typecode* code,unsigned* listp)
{
  unsigned count=0;
  if(n){
    //-lspg size list initialized to zero.
    //-Inicializa tamanho de lista lspg a cero.
    cudaMemset(listp+nmax,0,sizeof(unsigned));
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    const unsigned smem=(SPHBSIZE*2+1)*sizeof(unsigned); //-Each particle can leave two new periodic over the counter position. | De cada particula pueden salir 2 nuevas periodicas mas la posicion del contador.
    KerPeriodicMakeList <<<sgrid,SPHBSIZE,smem>>> (n,pini,nmax,Double3(mapposmin),Double3(mapposmax),Double3(perinc),posxy,posz,code,listp);
    cudaMemcpy(&count,listp+nmax,sizeof(unsigned),cudaMemcpyDeviceToHost);
    //-Reorders list if it is valid and stable has been activated.
    //-Reordena lista si es valida y stable esta activado.
    if(stable && count && count<=nmax){
      thrust::device_ptr<unsigned> dev_list(listp);
      thrust::sort(dev_list,dev_list+count);
    }
  }
  return(count);
}

//------------------------------------------------------------------------------
/// Doubles the position of the indicated particle using a displacement.
/// Duplicate particles are considered valid and are always within
/// the domain.
/// This kernel applies to single-GPU and multi-GPU because the calculations are made
/// from domposmin.
/// It controls the cell coordinates not exceed the maximum.
///
/// Duplica la posicion de la particula indicada aplicandole un desplazamiento.
/// Las particulas duplicadas se considera que siempre son validas y estan dentro
/// del dominio.
/// Este kernel vale para single-gpu y multi-gpu porque los calculos se hacen 
/// a partir de domposmin.
/// Se controla que las coordendas de celda no sobrepasen el maximo.
//------------------------------------------------------------------------------
__device__ void KerPeriodicDuplicatePos(unsigned pnew,unsigned pcopy
  ,bool inverse,double dx,double dy,double dz,uint3 cellmax
  ,double2* posxy,double* posz,unsigned* dcell)
{
  //-Obtains position of the particle to be duplicated.
  //-Obtiene pos de particula a duplicar.
  double2 rxy=posxy[pcopy];
  double rz=posz[pcopy];
  //-Applies displacement.
  rxy.x+=(inverse? -dx: dx);
  rxy.y+=(inverse? -dy: dy);
  rz+=(inverse? -dz: dz);
  //-Computes cell coordinates within the domain.
  //-Calcula coordendas de celda dentro de dominio.
  unsigned cx=unsigned((rxy.x-CTE.domposminx)/CTE.scell);
  unsigned cy=unsigned((rxy.y-CTE.domposminy)/CTE.scell);
  unsigned cz=unsigned((rz-CTE.domposminz)/CTE.scell);
  //-Adjust cell coordinates if they exceed the maximum.
  //-Ajusta las coordendas de celda si sobrepasan el maximo.
  cx=(cx<=cellmax.x? cx: cellmax.x);
  cy=(cy<=cellmax.y? cy: cellmax.y);
  cz=(cz<=cellmax.z? cz: cellmax.z);
  //-Stores position and cell of the new particles.
  //-Graba posicion y celda de nuevas particulas.
  posxy[pnew]=rxy;
  posz[pnew]=rz;
  dcell[pnew]=DCEL_Cell(CTE.cellcode,cx,cy,cz);
}

//------------------------------------------------------------------------------
/// Creates periodic particles from a list of particles to duplicate.
/// It is assumed that all particles are valid.
/// This kernel applies to single-GPU and multi-GPU because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-gpu y multi-gpu porque usa domposmin. 
//------------------------------------------------------------------------------
__global__ void KerPeriodicDuplicateVerlet(unsigned n,unsigned pini,uint3 cellmax
  ,double3 perinc,const unsigned* listp,unsigned* idp,typecode* code,unsigned* dcell
  ,double2* posxy,double* posz,float4* velrho,tsymatrix3f* spstau,float4* velrhom1)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Adjusts cell position of the new particles.
    //-Ajusta posicion y celda de nueva particula.
    KerPeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000)
      ,perinc.x,perinc.y,perinc.z,cellmax,posxy,posz,dcell);
    //-Copies the remaining data.
    //-Copia el resto de datos.
    idp     [pnew]=idp[pcopy];
    code    [pnew]=CODE_SetPeriodic(code[pcopy]);
    velrho  [pnew]=velrho[pcopy];
    velrhom1[pnew]=velrhom1[pcopy];
    if(spstau)spstau[pnew]=spstau[pcopy];
  }
}

//==============================================================================
/// Creates periodic particles from a list of particles to duplicate.
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
//==============================================================================
void PeriodicDuplicateVerlet(unsigned n,unsigned pini,tuint3 domcells
  ,tdouble3 perinc,const unsigned* listp,unsigned* idp,typecode* code
  ,unsigned* dcell,double2* posxy,double* posz,float4* velrho
  ,tsymatrix3f* spstau,float4* velrhom1)
{
  if(n){
    uint3 cellmax=make_uint3(domcells.x-1,domcells.y-1,domcells.z-1);
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerPeriodicDuplicateVerlet <<<sgrid,SPHBSIZE>>> (n,pini,cellmax,Double3(perinc)
      ,listp,idp,code,dcell,posxy,posz,velrho,spstau,velrhom1);
  }
}

//------------------------------------------------------------------------------
/// Creates periodic particles from a list of particles to duplicate.
/// It is assumed that all particles are valid.
/// This kernel applies to single-GPU and multi-GPU because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-gpu y multi-gpu porque usa domposmin. 
//------------------------------------------------------------------------------
template<bool varspre> __global__ void KerPeriodicDuplicateSymplectic(unsigned n
  ,unsigned pini,uint3 cellmax,double3 perinc,const unsigned* listp,unsigned* idp
  ,typecode* code,unsigned* dcell,double2* posxy,double* posz,float4* velrho
  ,tsymatrix3f* spstau,double2* posxypre,double* poszpre,float4* velrhopre)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Adjusts cell position of the new particles.
    //-Ajusta posicion y celda de nueva particula.
    KerPeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000),perinc.x,perinc.y,perinc.z,cellmax,posxy,posz,dcell);
    //-Copies the remaining data.
    //-Copia el resto de datos.
    idp   [pnew]=idp[pcopy];
    code  [pnew]=CODE_SetPeriodic(code[pcopy]);
    velrho[pnew]=velrho[pcopy];
    if(varspre){
      posxypre [pnew]=posxypre[pcopy];
      poszpre  [pnew]=poszpre[pcopy];
      velrhopre[pnew]=velrhopre[pcopy];
    }
    if(spstau)spstau[pnew]=spstau[pcopy];
  }
}

//==============================================================================
/// Creates periodic particles from a list of particles to duplicate.
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
//==============================================================================
void PeriodicDuplicateSymplectic(unsigned n,unsigned pini
  ,tuint3 domcells,tdouble3 perinc,const unsigned* listp,unsigned* idp
  ,typecode* code,unsigned* dcell,double2* posxy,double* posz,float4* velrho
  ,tsymatrix3f* spstau,double2* posxypre,double* poszpre,float4* velrhopre)
{
  if(n){
    uint3 cellmax=make_uint3(domcells.x-1,domcells.y-1,domcells.z-1);
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    if(posxypre!=NULL)KerPeriodicDuplicateSymplectic<true>  <<<sgrid,SPHBSIZE>>> 
      (n,pini,cellmax,Double3(perinc),listp,idp,code,dcell,posxy,posz,velrho,spstau
        ,posxypre,poszpre,velrhopre);
    else              KerPeriodicDuplicateSymplectic<false> <<<sgrid,SPHBSIZE>>> 
      (n,pini,cellmax,Double3(perinc),listp,idp,code,dcell,posxy,posz,velrho,spstau
        ,posxypre,poszpre,velrhopre);
  }
}

//------------------------------------------------------------------------------
/// Creates periodic particles from a list of particles to duplicate.
/// It is assumed that all particles are valid.
/// This kernel applies to single-GPU and multi-GPU because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-gpu y multi-gpu porque usa domposmin. 
//------------------------------------------------------------------------------
__global__ void KerPeriodicDuplicateNormals(unsigned n,unsigned pini
  ,const unsigned* listp,float3* normals,float3* motionvel,float3* motionace)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    normals[pnew]=normals[pcopy];
    if(motionvel)motionvel[pnew]=motionvel[pcopy];
    if(motionace)motionace[pnew]=motionace[pcopy];
  }
}

//==============================================================================
/// Creates periodic particles from a list of particles to duplicate.
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
//==============================================================================
void PeriodicDuplicateNormals(unsigned n,unsigned pini,const unsigned* listp
  ,float3* normals,float3* motionvel,float3* motionace)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerPeriodicDuplicateNormals <<<sgrid,SPHBSIZE>>> (n,pini,listp,normals,motionvel,motionace);
  }
}

//##############################################################################
//# Kernels for Damping.
//##############################################################################
//------------------------------------------------------------------------------
/// Returns TRUE when code==NULL or particle is normal and fluid.
//------------------------------------------------------------------------------
__device__ bool KerIsNormalFluid(const typecode* code,unsigned p){
  if(code){//-Descarta particulas floating o periodicas.
    const typecode cod=code[p];
    return(CODE_IsNormal(cod) && CODE_IsFluid(cod));
  }
  return(true);
}
//------------------------------------------------------------------------------
/// Checks position is inside box limits.
/// Comprueba si la posicion esta dentro de los limites.
//------------------------------------------------------------------------------
__device__ bool KerPointInBox(double px,double py,double pz,const double3& p1
  ,const double3& p2)
{
  return(p1.x<=px && p1.y<=py && p1.z<=pz && px<=p2.x && py<=p2.y && pz<=p2.z);
}
//------------------------------------------------------------------------------
/// Solves point on the plane.
/// Resuelve punto en el plano.
//------------------------------------------------------------------------------
__device__ double KerPointPlane(const double4& pla,double px,double py,double pz)
{
  return(pla.x*px+pla.y*py+pla.z*pz+pla.w);
}
//------------------------------------------------------------------------------
/// Solves point on the plane.
/// Resuelve punto en el plano.
//------------------------------------------------------------------------------
__device__ double KerPointPlane(const double4& pla,const double3& pt)
{
  return(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w);
}

//------------------------------------------------------------------------------
/// Applies Damping.
/// Aplica Damping.
//------------------------------------------------------------------------------
__global__ void KerComputeDampingPlane(unsigned n,unsigned pini
  ,double dt,double4 plane,float dist,float over,float3 factorxyz,float redumax
  ,const double2* posxy,const double* posz,const typecode* code
  ,float4* velrho)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pini;
    const bool ok=KerIsNormalFluid(code,p1);//-Ignore floating and periodic particles. | Descarta particulas floating o periodicas.
    if(ok){
      const double2 rposxy=posxy[p1];
      const double rposz=posz[p1];
      double vdis=KerPointPlane(plane,rposxy.x,rposxy.y,rposz);  //fgeo::PlanePoint(plane,ps);
      if(0<vdis && vdis<=dist+over){
        const double fdis=(vdis>=dist? 1.: vdis/dist);
        const double redudt=dt*(fdis*fdis)*redumax;
        double redudtx=(1.-redudt*factorxyz.x);
        double redudty=(1.-redudt*factorxyz.y);
        double redudtz=(1.-redudt*factorxyz.z);
        redudtx=(redudtx<0? 0.: redudtx);
        redudty=(redudty<0? 0.: redudty);
        redudtz=(redudtz<0? 0.: redudtz);
        float4 rvel=velrho[p1];
        rvel.x=float(redudtx*rvel.x); 
        rvel.y=float(redudty*rvel.y); 
        rvel.z=float(redudtz*rvel.z);
        velrho[p1]=rvel;
      }
    }
  }
}
//==============================================================================
/// Applies Damping.
/// Aplica Damping.
//==============================================================================
void ComputeDampingPlane(double dt,double4 plane,float dist,float over
  ,float3 factorxyz,float redumax,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,const typecode* code,float4* velrho)
{
  if(n){
    dim3 sgridf=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeDampingPlane <<<sgridf,SPHBSIZE>>> (n,pini,dt,plane,dist,over
      ,factorxyz,redumax,posxy,posz,code,velrho);
  }
}

//------------------------------------------------------------------------------
/// Applies Damping to limited domain.
/// Aplica Damping limitado a un dominio.
//------------------------------------------------------------------------------
__global__ void KerComputeDampingPlaneDom(unsigned n,unsigned pini
  ,double dt,double4 plane,float dist,float over,float3 factorxyz,float redumax
  ,double zmin,double zmax,double4 pla0,double4 pla1,double4 pla2,double4 pla3
  ,const double2* posxy,const double* posz,const typecode* code
  ,float4* velrho)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pini;
    const bool ok=KerIsNormalFluid(code,p1);//-Ignore floating and periodic particles. | Descarta particulas floating o periodicas.
    if(ok){
      const double2 rposxy=posxy[p1];
      const double rposz=posz[p1];
      const double3 ps=make_double3(rposxy.x,rposxy.y,rposz);
      double vdis=KerPointPlane(plane,ps);  //fgeo::PlanePoint(plane,ps);
      if(0<vdis && vdis<=dist+over){
        if(ps.z>=zmin && ps.z<=zmax && KerPointPlane(pla0,ps)<=0 && KerPointPlane(pla1,ps)<=0 && KerPointPlane(pla2,ps)<=0 && KerPointPlane(pla3,ps)<=0){
          const double fdis=(vdis>=dist? 1.: vdis/dist);
          const double redudt=dt*(fdis*fdis)*redumax;
          double redudtx=(1.-redudt*factorxyz.x);
          double redudty=(1.-redudt*factorxyz.y);
          double redudtz=(1.-redudt*factorxyz.z);
          redudtx=(redudtx<0? 0.: redudtx);
          redudty=(redudty<0? 0.: redudty);
          redudtz=(redudtz<0? 0.: redudtz);
          float4 rvel=velrho[p1];
          rvel.x=float(redudtx*rvel.x); 
          rvel.y=float(redudty*rvel.y); 
          rvel.z=float(redudtz*rvel.z); 
          velrho[p1]=rvel;
        }
      }
    }
  }
}
//==============================================================================
/// Applies Damping to limited domain.
/// Aplica Damping limitado a un dominio.
//==============================================================================
void ComputeDampingPlaneDom(double dt,double4 plane,float dist,float over
  ,float3 factorxyz,float redumax
  ,double zmin,double zmax,double4 pla0,double4 pla1,double4 pla2,double4 pla3
  ,unsigned n,unsigned pini,const double2* posxy,const double* posz
  ,const typecode* code,float4* velrho)
{
  if(n){
    dim3 sgridf=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeDampingPlaneDom <<<sgridf,SPHBSIZE>>> (n,pini,dt,plane,dist,over,factorxyz
      ,redumax,zmin,zmax,pla0,pla1,pla2,pla3,posxy,posz,code,velrho);
  }
}


//------------------------------------------------------------------------------
/// Applies Damping according box configuration.
/// Aplica Damping segun cofiguracion de caja.
//------------------------------------------------------------------------------
__global__ void KerComputeDampingBox(unsigned n,unsigned pini
  ,double dt,float3 factorxyz,float redumax
  ,double3 limitmin1,double3 limitmin2,double3 limitmax1,double3 limitmax2
  ,double3 limitover1,double3 limitover2,double3 boxsize1,double3 boxsize2
  ,const double2* posxy,const double* posz,const typecode* code
  ,float4* velrho)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pini;
    const bool ok=KerIsNormalFluid(code,p1);//-Ignore floating and periodic particles. | Descarta particulas floating o periodicas.
    if(ok){
      const double2 rposxy=posxy[p1];
      const double rposz=posz[p1];
      //-Check if it is within the domain. | Comprueba si esta dentro del dominio.
      if(KerPointInBox(rposxy.x,rposxy.y,rposz,limitover1,limitover2)){//-Inside overlimit domain.
        if(!KerPointInBox(rposxy.x,rposxy.y,rposz,limitmin1,limitmin2)){//-Outside free domain.
          double fdis=1.;
          if(KerPointInBox(rposxy.x,rposxy.y,rposz,limitmax1,limitmax2)){//-Compute damping coefficient.
            fdis=0;
            if(boxsize2.z){ const double fdiss=(rposz   -limitmin2.z)/boxsize2.z; fdis=(fdis>=fdiss? fdis: fdiss); }
            if(boxsize2.y){ const double fdiss=(rposxy.y-limitmin2.y)/boxsize2.y; fdis=(fdis>=fdiss? fdis: fdiss); }
            if(boxsize2.x){ const double fdiss=(rposxy.x-limitmin2.x)/boxsize2.x; fdis=(fdis>=fdiss? fdis: fdiss); }
            if(boxsize1.z){ const double fdiss=(limitmin1.z-rposz   )/boxsize1.z; fdis=(fdis>=fdiss? fdis: fdiss); }
            if(boxsize1.y){ const double fdiss=(limitmin1.y-rposxy.y)/boxsize1.y; fdis=(fdis>=fdiss? fdis: fdiss); }
            if(boxsize1.x){ const double fdiss=(limitmin1.x-rposxy.x)/boxsize1.x; fdis=(fdis>=fdiss? fdis: fdiss); }
          }
          const double redudt=dt*(fdis*fdis)*redumax;
          double redudtx=(1.-redudt*factorxyz.x);
          double redudty=(1.-redudt*factorxyz.y);
          double redudtz=(1.-redudt*factorxyz.z);
          redudtx=(redudtx<0? 0.: redudtx);
          redudty=(redudty<0? 0.: redudty);
          redudtz=(redudtz<0? 0.: redudtz);
          float4 rvel=velrho[p1];
          rvel.x=float(redudtx*rvel.x); 
          rvel.y=float(redudty*rvel.y); 
          rvel.z=float(redudtz*rvel.z);
          //rvel.x=rvel.y=rvel.z=0;
          velrho[p1]=rvel;
        }
      }
    }
  }
}
//==============================================================================
/// Applies Damping according box configuration.
/// Aplica Damping segun cofiguracion de caja.
//==============================================================================
void ComputeDampingBox(unsigned n,unsigned pini,double dt,float3 factorxyz,float redumax
  ,double3 limitmin1,double3 limitmin2,double3 limitmax1,double3 limitmax2
  ,double3 limitover1,double3 limitover2,double3 boxsize1,double3 boxsize2
  ,const double2* posxy,const double* posz,const typecode* code,float4* velrho)
{
  if(n){
    dim3 sgridf=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeDampingBox <<<sgridf,SPHBSIZE>>> (n,pini,dt,factorxyz,redumax
      ,limitmin1,limitmin2,limitmax1,limitmax2,limitover1,limitover2,boxsize1,boxsize2
      ,posxy,posz,code,velrho);
  }
}


//------------------------------------------------------------------------------
/// Applies Damping to limited cylinder domain.
/// Aplica Damping limitado a un dominio de cilindro.
//------------------------------------------------------------------------------
__global__ void KerComputeDampingCylinder(unsigned n,unsigned pini
  ,double dt,bool isvertical,double3 point1,double3 point2,double limitmin
  ,float dist,float over,float3 factorxyz,float redumax
  ,const double2* posxy,const double* posz,const typecode* code
  ,float4* velrho)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pini;
    const bool ok=KerIsNormalFluid(code,p1);//-Ignore floating and periodic particles. | Descarta particulas floating o periodicas.
    if(ok){
      //-Check if it is within the domain. | Comprueba si esta dentro del dominio.
      const double2 rposxy=posxy[p1];
      const double rposz=posz[p1];
      const double3 ps=make_double3(rposxy.x,rposxy.y,rposz);
      const double vdis=(isvertical? 
        sqrt((ps.x-point1.x)*(ps.x-point1.x)+(ps.y-point1.y)*(ps.y-point1.y)): 
        cugeo::LinePointDist(ps,point1,point2)
        ) - limitmin;
      if(0<vdis && vdis<=dist+over){
        const double fdis=(vdis>=dist? 1.: vdis/dist);
        const double redudt=dt*(fdis*fdis)*redumax;
        double redudtx=(1.-redudt*factorxyz.x);
        double redudty=(1.-redudt*factorxyz.y);
        double redudtz=(1.-redudt*factorxyz.z);
        redudtx=(redudtx<0? 0.: redudtx);
        redudty=(redudty<0? 0.: redudty);
        redudtz=(redudtz<0? 0.: redudtz);
        float4 rvel=velrho[p1];
        rvel.x=float(redudtx*rvel.x); 
        rvel.y=float(redudty*rvel.y); 
        rvel.z=float(redudtz*rvel.z); 
        velrho[p1]=rvel;
      }
    }
  }
}
//==============================================================================
/// Applies Damping to limited cylinder domain.
/// Aplica Damping limitado a un dominio de cilindro.
//==============================================================================
void ComputeDampingCylinder(unsigned n,unsigned pini
  ,double dt,double3 point1,double3 point2,double limitmin
  ,float dist,float over,float3 factorxyz,float redumax
  ,const double2* posxy,const double* posz,const typecode* code
  ,float4* velrho)
{
  if(n){
    const bool isvertical=(point1.x==point2.x && point1.y==point2.y);
    dim3 sgridf=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeDampingCylinder <<<sgridf,SPHBSIZE>>> (n,pini,dt
      ,isvertical,point1,point2,limitmin,dist,over,factorxyz,redumax
      ,posxy,posz,code,velrho);
  }
}

 //<vs_outpaarts_ini>
//##############################################################################
//# Kernels for OutputParts.
//##############################################################################
//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
__global__ void KerComputeOutputPartsInit(unsigned n,unsigned pini
  ,byte resmask,bool selall,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    sel[p]=(selall? psel|resmask: psel&(~resmask));
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsInit(byte resmask,bool selall,unsigned n,unsigned pini
  ,byte* sel)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeOutputPartsInit <<<sgrid,SPHBSIZE>>> (n,pini,resmask,selall,sel);
  }
}

//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
__global__ void KerComputeOutputPartsGroup(unsigned n,unsigned pini
  ,byte resmask,byte resprev,bool cmband,bool inverse,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    bool r=((psel&resprev)!=0);  //-lv=0
    bool ok=((psel&resmask)!=0); //-lv+1
    if(inverse)ok=!ok;
    r=(cmband? r&&ok: r||ok);
    psel=psel&(~resprev);
    if(r)psel=psel|resprev;
    sel[p]=psel;
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsGroup(byte resmask,byte resprev,bool cmband,bool inverse
  ,unsigned n,unsigned pini,byte* sel)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeOutputPartsGroup <<<sgrid,SPHBSIZE>>> (n,pini,resmask,resprev,cmband,inverse,sel);
  }
}

//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
__global__ void KerComputeOutputPartsPos(unsigned n,unsigned pini
  ,byte resmask,bool cmband,bool inverse,double3 pmin,double3 pmax
  ,const double2* posxy,const double* posz,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==cmband){ //if((r && cmband) || (!r && !cmband)){
      const double2 rposxy=posxy[p];
      const double rposz=posz[p];
      bool ok=(pmin.x<=rposxy.x && rposxy.x <=pmax.x &&
               pmin.z<=rposz    && rposz    <=pmax.z &&
               pmin.y<=rposxy.y && rposxy.y <=pmax.y);
      if(inverse)ok=!ok;
      r=(cmband? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsPos(byte resmask,bool cmband,bool inverse
  ,double3 pmin,double3 pmax,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,byte* sel)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeOutputPartsPos <<<sgrid,SPHBSIZE>>> (n,pini,resmask,cmband,inverse
      ,pmin,pmax,posxy,posz,sel);
  }
}

//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
__global__ void KerComputeOutputPartsPlane(unsigned n,unsigned pini
  ,byte resmask,bool cmband,bool inverse,double4 plane,float maxdist
  ,const double2* posxy,const double* posz,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==cmband){ //if((r && cmband) || (!r && !cmband)){
      const double2 rposxy=posxy[p];
      const double rposz=posz[p];
      const double dist=KerPointPlane(plane,rposxy.x,rposxy.y,rposz);  //fgeo::PlanePoint(plane,ps);
      bool ok=(dist>=0 && dist<=maxdist);
      if(inverse)ok=!ok;
      r=(cmband? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsPlane(byte resmask,bool cmband,bool inverse
  ,double4 plane,float maxdist,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,byte* sel)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeOutputPartsPlane <<<sgrid,SPHBSIZE>>> (n,pini,resmask,cmband,inverse
      ,plane,maxdist,posxy,posz,sel);
  }
}

//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
__global__ void KerComputeOutputPartsSphere(unsigned n,unsigned pini
  ,byte resmask,bool cmband,bool inverse,double3 pcen,float radius2
  ,const double2* posxy,const double* posz,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==cmband){ //if((r && cmband) || (!r && !cmband)){
      const double2 rposxy=posxy[p];
      const double rposz=posz[p];
      const float dx=float(pcen.x-rposxy.x);
      const float dy=float(pcen.y-rposxy.y);
      const float dz=float(pcen.z-rposz);
      bool ok=(dx*dx+dy*dy+dz*dz <= radius2);
      if(inverse)ok=!ok;
      r=(cmband? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsSphere(byte resmask,bool cmband,bool inverse
  ,double3 pcen,float radius2,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,byte* sel)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeOutputPartsSphere <<<sgrid,SPHBSIZE>>> (n,pini,resmask,cmband,inverse
      ,pcen,radius2,posxy,posz,sel);
  }
}


//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
template<bool isvertical> __global__ void KerComputeOutputPartsCylinder(unsigned n
  ,unsigned pini,byte resmask,bool cmband,bool inverse
  ,double4 plane,float maxdist,double3 pcen1,double3 pcen2,float radius
  ,const double2* posxy,const double* posz,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==cmband){ //if((r && cmband) || (!r && !cmband)){
      const double2 rposxy=posxy[p];
      const double rposz=posz[p];
      const double dist=KerPointPlane(plane,rposxy.x,rposxy.y,rposz);  //fgeo::PlanePoint(plane,ps);
      bool ok=(dist>=0 && dist<=maxdist);
      if(ok && isvertical){
        const float dx=float(rposxy.x-pcen1.x);
        const float dy=float(rposxy.y-pcen1.y);
        ok=(dx*dx+dy*dy<=radius);
      }
      if(ok && !isvertical){
        //cugeo::LinePointDist(ps,Point1,Point2)
        const double ar=cugeo::TriangleArea(make_double3(rposxy.x,rposxy.y,rposz),pcen1,pcen2);
        const double dis=(ar*2)/maxdist;
        ok=(dis<=radius);
      }
      if(inverse)ok=!ok;
      r=(cmband? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsCylinder(byte resmask,bool cmband,bool inverse
  ,double4 plane,float maxdist,double3 pcen1,double3 pcen2,float radius
  ,unsigned n,unsigned pini,const double2* posxy,const double* posz,byte* sel)
{
  if(n){
    const bool isvertical=(pcen1.x==pcen2.x && pcen1.y==pcen2.y);
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    if(isvertical) KerComputeOutputPartsCylinder<true > <<<sgrid,SPHBSIZE>>> (n,pini
      ,resmask,cmband,inverse,plane,maxdist,pcen1,pcen2,radius,posxy,posz,sel);
    if(!isvertical)KerComputeOutputPartsCylinder<false> <<<sgrid,SPHBSIZE>>> (n,pini
      ,resmask,cmband,inverse,plane,maxdist,pcen1,pcen2,radius,posxy,posz,sel);
  }
}

//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
__global__ void KerComputeOutputPartsType(unsigned n,unsigned pini
  ,byte resmask,bool cmband,bool inverse,byte types
  ,const typecode* code,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==cmband){ //if((r && cmband) || (!r && !cmband)){
      const typecode codetp=CODE_GetType(code[p]);
      bool ok=((types&1 && codetp==CODE_TYPE_FIXED) ||
               (types&2 && codetp==CODE_TYPE_MOVING) ||
               (types&4 && codetp==CODE_TYPE_FLOATING) ||
               (types&8 && codetp==CODE_TYPE_FLUID) );
      if(inverse)ok=!ok;
      r=(cmband? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsType(byte resmask,bool cmband,bool inverse
  ,byte types,unsigned n,unsigned pini,const typecode* code,byte* sel)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeOutputPartsType <<<sgrid,SPHBSIZE>>> (n,pini,resmask,cmband,inverse
      ,types,code,sel);
  }
}

//------------------------------------------------------------------------------
/// Compute filter for particles.
//------------------------------------------------------------------------------
__global__ void KerComputeOutputPartsMk(unsigned n,unsigned pini
  ,byte resmask,bool cmband,bool inverse,typecode mkcode1,typecode mkcode2
  ,const typecode* code,byte* sel)
{
  unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==cmband){ //if((r && cmband) || (!r && !cmband)){
      const typecode codetp=CODE_GetTypeAndValue(code[p]);
      bool ok=(mkcode1<=codetp && codetp<=mkcode2);
      if(inverse)ok=!ok;
      r=(cmband? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void ComputeOutputPartsMk(byte resmask,bool cmband,bool inverse
  ,typecode mkcode1,typecode mkcode2,unsigned n,unsigned pini
  ,const typecode* code,byte* sel)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerComputeOutputPartsMk <<<sgrid,SPHBSIZE>>> (n,pini,resmask,cmband,inverse
      ,mkcode1,mkcode2,code,sel);
  }
}

//<vs_outpaarts_end>



}


//##############################################################################
//# Kernels for InOut (JSphInOut).
//# Kernels para InOut (JSphInOut).
//##############################################################################
#include "JSphGpu_InOut_iker.cu"


