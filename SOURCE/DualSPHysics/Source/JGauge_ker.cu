//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpu_ker.cu \brief Implements functions and CUDA kernels for class JGauge.

#include "JGauge_ker.h"
#include <float.h>
#include <math_constants.h>
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"
#include <cstdio>
#include <string>


namespace cugauge{

//==============================================================================
/// Checks error and ends execution.
/// Comprueba error y finaliza ejecucion.
//==============================================================================
#define CheckErrorCuda(text)  __CheckErrorCuda(text,__FILE__,__LINE__)
void __CheckErrorCuda(const char *text,const char *file,const int line){
  cudaError_t err=cudaGetLastError();
  if(cudaSuccess!=err){
    char cad[2048]; 
    sprintf(cad,"%s (CUDA error: %s -> %s:%i).\n",text,cudaGetErrorString(err),file,line); 
    throw std::string(cad);
  }
}

//==============================================================================
/// Returns size of gridsize according to parameters.
/// Devuelve tamaño de gridsize segun parametros.
//==============================================================================
dim3 GetGridSize(unsigned n,unsigned blocksize){
  dim3 sgrid;//=dim3(1,2,3);
  unsigned nb=unsigned(n+blocksize-1)/blocksize;//-Numero total de bloques a lanzar.
  sgrid.x=(nb<=65535? nb: unsigned(sqrt(float(nb))));
  sgrid.y=(nb<=65535? 1: unsigned((nb+sgrid.x-1)/sgrid.x));
  sgrid.z=1;
  return(sgrid);
}


//##############################################################################
//# Kernels for gauge interaction.
//##############################################################################

//------------------------------------------------------------------------------
/// Returns cell limits for the interaction.
/// Devuelve limites de celdas para interaccion.
//------------------------------------------------------------------------------
__device__ void KerGetInteractionCells(double px,double py,double pz
  ,int hdiv,const int4 &nc,const int3 &cellzero
  ,int &cxini,int &cxfin,int &yini,int &yfin,int &zini,int &zfin
  ,const double3 &domposmin,float scell)
{
  //-Obtains interaction limits.
  const int cx=int((px-domposmin.x)/scell)-cellzero.x;
  const int cy=int((py-domposmin.y)/scell)-cellzero.y;
  const int cz=int((pz-domposmin.z)/scell)-cellzero.z;
  //-Code for hdiv 1 or 2. The result is always within the range [0,nc-1+1].
  //-Codigo para hdiv 1 o 2. El resultado siempre esta dentro del intervalo [0,nc-1+1].
  cxini=cx-min(cx,hdiv);
  cxfin=cx+min(nc.x-cx-1,hdiv)+1;
  yini=cy-min(cy,hdiv);
  yfin=cy+min(nc.y-cy-1,hdiv)+1;
  zini=cz-min(cz,hdiv);
  zfin=cz+min(nc.z-cz-1,hdiv)+1;
}

//------------------------------------------------------------------------------
/// Performs interaction between particles. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// It includes visco artificial/laminar and floatings SPH/DEM.
///
/// Realiza interaccion entre particulas. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings SPH/DEM.
//------------------------------------------------------------------------------
__global__ void KerInteractionGaugeVel(double3 ptpos
  ,float awen,int hdiv,int4 nc,int3 cellzero,unsigned cellfluid,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,float3 *ptvel
  ,double3 domposmin,float scell,float fourh2,float h,float massf)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; //-Nº de la partícula
  if(!p){
    const double px=ptpos.x;
    const double py=ptpos.y;
    const double pz=ptpos.z;

    //-Obtains interaction limits.
    int cxini,cxfin,yini,yfin,zini,zfin;
    KerGetInteractionCells(px,py,pz,hdiv,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin,domposmin,scell);

    double sumwab=0;
    double3 sumvel=make_double3(0,0,0);
    //-Interaction with fluid particles. | Interaccion con fluidas.
    for(int z=zini;z<zfin;z++){
      int zmod=(nc.w)*z+cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        unsigned pini,pfin=0;
        for(int x=cxini;x<cxfin;x++){
          int2 cbeg=begincell[x+ymod];
          if(cbeg.y){
            if(!pfin)pini=cbeg.x;
            pfin=cbeg.y;
          }
        }
        if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
          double2 pxy2=posxy[p2];
          const float drx=float(px-pxy2.x);
          const float dry=float(py-pxy2.y);
          const float drz=float(pz-posz[p2]);
          const float rr2=(drx*drx+dry*dry+drz*drz);
          //-Interaction with real neighboring particles. 
          //-Interaccion con particulas vecinas reales.
          if(rr2<=fourh2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
            float wab;
            {//-Wendland kernel.
              const float qq=sqrt(rr2)/h;
              const float wqq=2.f*qq+1.f;
              const float wqq1=1.f-0.5f*qq;
              const float wqq2=wqq1*wqq1;
              wab=awen*wqq*wqq2*wqq2;
            }
            float4 velrhopp2=velrhop[p2];
            wab*=massf/velrhopp2.w;
            sumwab+=wab;
            sumvel.x+=wab*velrhopp2.x;
            sumvel.y+=wab*velrhopp2.y;
            sumvel.z+=wab*velrhopp2.z;
          }
        }
      }
    }
    //-Aplies kernel correction.
    //if(sumwab){
    //  sumvel.x/=sumwab;
    //  sumvel.y/=sumwab;
    //  sumvel.z/=sumwab;
    //}
    //-Stores result. | Guarda resultado.
    ptvel[0]=make_float3(float(sumvel.x),float(sumvel.y),float(sumvel.z));
  }
}
//==============================================================================
/// Calculates velocity in indicated point.
//==============================================================================
void Interaction_GaugeVel(tdouble3 ptpos
  ,float awen,int hdiv,tuint3 ncells,tuint3 cellmin,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,float3 *ptvel
  ,tdouble3 domposmin,float scell,float fourh2,float h,float massf)
{
  const int4 nc=make_int4(ncells.x,ncells.y,ncells.z,ncells.x*ncells.y);
  const unsigned cellfluid=nc.w*nc.z+1;
  const int3 cellzero=make_int3(cellmin.x,cellmin.y,cellmin.z);
  //-Interaction Fluid-Fluid & Fluid-Bound.
  if(1){
    const unsigned bsize=32;
    dim3 sgrid=GetGridSize(1,bsize);
    //:JDgKerPrint info;
    //:byte* ik=NULL; //info.GetInfoPointer(sgridf,bsfluid);
    KerInteractionGaugeVel <<<sgrid,bsize>>> (Double3(ptpos),awen,hdiv,nc,cellzero,cellfluid,begincell,posxy,posz,code,velrhop,ptvel,Double3(domposmin),scell,fourh2,h,massf);
    //:info.PrintValuesFull(true); //info.PrintValuesInfo();
  }
}

//------------------------------------------------------------------------------
/// Calculates mass value at one point by interacting with the fluid.
/// Calcula valor de masa en un punto mediante la interaccion con el fluido.
//------------------------------------------------------------------------------
__device__ float KerCalculeMass(double px,double py,double pz,float awen
  ,int hdiv,int4 nc,int3 cellzero,unsigned cellfluid,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,double3 domposmin,float scell,float fourh2,float h,float massf)
{
  //-Obtains interaction limits.
  int cxini,cxfin,yini,yfin,zini,zfin;
  KerGetInteractionCells(px,py,pz,hdiv,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin,domposmin,scell);

  double summass=0;
  //-Interaction with fluid particles. | Interaccion con fluidas.
  for(int z=zini;z<zfin;z++){
    int zmod=(nc.w)*z+cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
    for(int y=yini;y<yfin;y++){
      int ymod=zmod+nc.x*y;
      unsigned pini,pfin=0;
      for(int x=cxini;x<cxfin;x++){
        int2 cbeg=begincell[x+ymod];
        if(cbeg.y){
          if(!pfin)pini=cbeg.x;
          pfin=cbeg.y;
        }
      }
      if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
        double2 pxy2=posxy[p2];
        const float drx=float(px-pxy2.x);
        const float dry=float(py-pxy2.y);
        const float drz=float(pz-posz[p2]);
        const float rr2=(drx*drx+dry*dry+drz*drz);
        //-Interaction with real neighboring particles. 
        //-Interaccion con particulas vecinas reales
        if(rr2<=fourh2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
          float wab;
          {//-Wendland kernel.
            const float qq=sqrt(rr2)/h;
            const float wqq=2.f*qq+1.f;
            const float wqq1=1.f-0.5f*qq;
            const float wqq2=wqq1*wqq1;
            wab=awen*wqq*wqq2*wqq2;
          }
          wab*=massf/velrhop[p2].w;
          summass+=wab*massf;
        }
      }
    }
  }
  return(float(summass));
}

//------------------------------------------------------------------------------
/// Calculate position of free surface in (px,py).
/// Returns DBL_MAX if it does not find the free surface.
/// Calcula posicion de la superficie libre en (px,py).
/// Devuelve DBL_MAX si no encuentra la superficie libre.
//------------------------------------------------------------------------------
__global__ void KerInteractionGaugeZsurf(double px,double py
  ,double gaugezmin,double gaugedp,float masslimit,unsigned czmin,unsigned czmax,float awen
  ,int hdiv,int4 nc,int3 cellzero,unsigned cellfluid,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,double *zsurfres
  ,double3 domposmin,float scell,float fourh2,float h,float massf)
{
  extern __shared__ float shmass[];
  const unsigned tid=threadIdx.x;
  double zpre=gaugezmin;
  float mpre=masslimit;
  unsigned czini=czmin;
  double zsurf=DBL_MAX;
  float fin=0;
  while(czini<=czmax && fin==0){//-Calculates mass in blocks of blockDim.x measurement points. | Calcula masa en bloques de blockDim.x puntos de medida.
    //-Calcula valores de masa.
    unsigned cz=czini+tid;
    if(cz<=czmax){
      shmass[tid]=KerCalculeMass(px,py,gaugezmin+gaugedp*cz,awen,hdiv,nc,cellzero,cellfluid,begincell,posxy,posz,code,velrhop,domposmin,scell,fourh2,h,massf);
    }
    //-Calculates mass values.
    __syncthreads();
    if(!tid){
      for(unsigned th=0;th<blockDim.x && czini+th<=czmax && zsurf==DBL_MAX;th++){
        double pz=gaugezmin+gaugedp*(czini+th);
        float mass=shmass[th];
        if(mass<masslimit){
          const double fx=(masslimit-mpre)/(mass-mpre);
          zsurf=fx*(pz-zpre)+zpre;
        }
        zpre=pz;
        mpre=mass;
      }
      shmass[blockDim.x]=(zsurf==DBL_MAX? 0: 1);
    }
    __syncthreads();
    fin=shmass[blockDim.x];
    czini+=blockDim.x;
  }
  if(!tid)zsurfres[0]=zsurf;
}

//==============================================================================
/// Calculates position of free surface.
/// Calcula posicion de la superficie libre.
//==============================================================================
void Interaction_GaugeZsurf(double px,double py
  ,double gaugezmin,double gaugedp,float masslimit,unsigned czmin,unsigned czmax
  ,float awen,int hdiv,tuint3 ncells,tuint3 cellmin,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,double *zsurf
  ,tdouble3 domposmin,float scell,float fourh2,float h,float massf)
{
  const int4 nc=make_int4(ncells.x,ncells.y,ncells.z,ncells.x*ncells.y);
  const unsigned cellfluid=nc.w*nc.z+1;
  const int3 cellzero=make_int3(cellmin.x,cellmin.y,cellmin.z);
  //-Interaction Fluid-Fluid.
  if(1){
    const unsigned bsize=64;
    dim3 sgrid=GetGridSize(bsize,bsize);
    const unsigned smem=sizeof(float)*(bsize+1);
    //:JDgKerPrint info;
    //:byte* ik=NULL; //info.GetInfoPointer(sgrid,bsize);
    KerInteractionGaugeZsurf <<<sgrid,bsize,smem>>> (px,py,gaugezmin,gaugedp,masslimit,czmin,czmax,awen,hdiv,nc,cellzero,cellfluid,begincell,posxy,posz,code,velrhop,zsurf,Double3(domposmin),scell,fourh2,h,massf);
    //:info.PrintValuesFull(true); //info.PrintValuesInfo();
  }
}


}


