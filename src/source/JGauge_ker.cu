//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpu_ker.cu \brief Implements functions and CUDA kernels for classes JGauge.

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
              wab=awen*wqq*wqq2*wqq2; //-Kernel.
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
            wab=awen*wqq*wqq2*wqq2; //-Kernel.
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
/// Calculates surface water level at indicated line.
//------------------------------------------------------------------------------
__global__ void KerInteractionGaugeSwl(double p0x,double p0y,double p0z
  ,double pdirx,double pdiry,double pdirz,unsigned pointnp,float masslimit
  ,float awen,int hdiv,int4 nc,int3 cellzero,unsigned cellfluid,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,double3 domposmin,float scell,float fourh2,float h,float massf,float3 *ptres)
{
  extern __shared__ float shmass[];
  const unsigned tid=threadIdx.x;
  unsigned cpbase=0;
  float psurfx=FLT_MAX,psurfy=FLT_MAX,psurfz=FLT_MAX;
  float mpre=0;
  while(cpbase<=pointnp){
    //-Saves mass values in shared memory.
    const unsigned cp=cpbase+tid;
    if(cp<=pointnp){
      shmass[tid]=KerCalculeMass(p0x+pdirx*cp,p0y+pdiry*cp,p0z+pdirz*cp,awen,hdiv,nc,cellzero,cellfluid,begincell,posxy,posz,code,velrhop,domposmin,scell,fourh2,h,massf);
    }
    else shmass[tid]=0;
    __syncthreads();
    //-Checks mass values.
    if(!tid){
      for(unsigned c=0;c<blockDim.x;c++){
        const float mass=shmass[c];
        if(mass>masslimit)mpre=mass;
        if(mass<masslimit && mpre){
          const float fxm1=((masslimit-mpre)/(mass-mpre)-1)+float(cpbase+c);
          psurfx=p0x+pdirx*fxm1;
          psurfy=p0y+pdiry*fxm1;
          psurfz=p0z+pdirz*fxm1;
          shmass[0]=FLT_MAX;
          break;
        }
      }
    }
    __syncthreads();
    if(shmass[0]==FLT_MAX)break;
    cpbase+=blockDim.x;
  }
  //-Stores result.
  if(!tid){
    if(psurfx==FLT_MAX){
      const unsigned cp=(mpre? pointnp: 0);
      psurfx=p0x+(pdirx*cp);
      psurfy=p0y+(pdiry*cp);
      psurfz=p0z+(pdirz*cp);
    }
    ptres[0]=make_float3(psurfx,psurfy,psurfz);
  }
}
//==============================================================================
/// Calculates surface water level at indicated line.
//==============================================================================
void Interaction_GaugeSwl(tdouble3 point0,tdouble3 pointdir,unsigned pointnp,float masslimit
  ,float awen,int hdiv,tuint3 ncells,tuint3 cellmin,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop
  ,tdouble3 domposmin,float scell,float fourh2,float h,float massf,float3 *ptres)
{
  const int4 nc=make_int4(ncells.x,ncells.y,ncells.z,ncells.x*ncells.y);
  const unsigned cellfluid=nc.w*nc.z+1;
  const int3 cellzero=make_int3(cellmin.x,cellmin.y,cellmin.z);
  if(1){
    const unsigned bsize=128;
    dim3 sgrid=GetGridSize(bsize,bsize);
    const unsigned smem=sizeof(float)*(bsize+1);
    //:JDgKerPrint info;
    //:byte* ik=NULL; //info.GetInfoPointer(sgrid,bsize);
    KerInteractionGaugeSwl <<<sgrid,bsize,smem>>> (point0.x,point0.y,point0.z,pointdir.x,pointdir.y,pointdir.z,pointnp,masslimit,awen,hdiv,nc,cellzero,cellfluid,begincell,posxy,posz,code,velrhop,Double3(domposmin),scell,fourh2,h,massf,ptres);
    //:info.PrintValuesFull(true); //info.PrintValuesInfo();
  }
}


//------------------------------------------------------------------------------
/// Calculates maximum z of fluid at distance of a vertical line.
//------------------------------------------------------------------------------
__global__ void KerInteractionGaugeMaxz(double p0x,double p0y,float maxdist2
  ,int cxini,int cxfin,int yini,int yfin,int zini,int zfin
  ,int4 nc,unsigned cellfluid,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code
  ,float3 *ptres)
{
  if(threadIdx.x==0){
    unsigned pmax=UINT_MAX;
    float zmax=-FLT_MAX;
    //-Interaction with fluid particles. | Interaccion con fluidas.
    for(int z=zfin-1;z>=zini && pmax==UINT_MAX;z--){
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
          const double posz2=posz[p2];
          if(posz2>zmax){
            const double2 posxy2=posxy[p2];
            const float drx=float(p0x-posxy2.x);
            const float dry=float(p0y-posxy2.y);
            const float rr2=drx*drx+dry*dry;
            if(rr2<=maxdist2 && CODE_IsFluid(code[p2])){//-Only with fluid particles.
              zmax=float(posz2);
              pmax=p2;
            }
          }
        }
      }
    }
    //-Stores result.
    ptres[0]=make_float3(0,0,zmax);
  }
}
//==============================================================================
/// Calculates maximum z of fluid at distance of a vertical line.
//==============================================================================
void Interaction_GaugeMaxz(tdouble3 point0,float maxdist2
  ,int cxini,int cxfin,int yini,int yfin,int zini,int zfin
  ,int4 nc,unsigned cellfluid,const int2 *begincell
  ,const double2 *posxy,const double *posz,const typecode *code
  ,float3 *ptres)
{
  const unsigned bsize=128;
  dim3 sgrid=GetGridSize(1,bsize);
  KerInteractionGaugeMaxz <<<sgrid,bsize>>> (point0.x,point0.y,maxdist2,cxini,cxfin,yini,yfin,zini,zfin,nc,cellfluid,begincell,posxy,posz,code,ptres);
}


//------------------------------------------------------------------------------
/// Calculates force on selected fixed or moving particles using only fluid particles.
/// Ignores periodic boundary particles to avoid race condition problems.
//------------------------------------------------------------------------------
__global__ void KerInteractionGaugeForce(unsigned n,unsigned idbegin,typecode codesel
  ,float fourh2,float h,float bwen,float massf,float cteb,float rhopzero,float gamma
  ,int hdiv,int4 nc,int3 cellzero,unsigned cellfluid,const int2 *begincell,double3 domposmin,float scell
  ,const double2 *posxy,const double *posz,const typecode *code,const unsigned *idp,const float4 *velrhop
  ,float3 *partace)
{
  unsigned p=blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const typecode code1=code[p];
    if(CODE_GetTypeAndValue(code1)==codesel && CODE_IsNormal(code1)){
      const double2 ptposxy=posxy[p];
      const double px=ptposxy.x;
      const double py=ptposxy.y;
      const double pz=posz[p];
      const float rhop1=velrhop[p].w;
      const float press1=cteb*(pow(rhop1/rhopzero,gamma)-1.0f);
      float3 ace=make_float3(0,0,0);

      //-Obtains interaction limits.
      int cxini,cxfin,yini,yfin,zini,zfin;
      KerGetInteractionCells(px,py,pz,hdiv,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin,domposmin,scell);
   
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
            //-Interaction with real neighboring fluid particles.
            if(rr2<=fourh2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
              float frx,fry,frz;
              {//-Wendland kernel.
                const float rad=sqrt(rr2);
                const float qq=rad/h;
                const float wqq1=1.f-0.5f*qq;
                const float fac=bwen*qq*wqq1*wqq1*wqq1/rad; //-Kernel derivative (divided by rad).
                frx=fac*drx; fry=fac*dry; frz=fac*drz;
              }
              const float mass2=massf;
              const float rhop2=velrhop[p2].w;
              const float press2=cteb*(pow(rhop2/rhopzero,gamma)-1.0f);
              const float prs=(press1+press2)/(rhop1*rhop2);
              {//-Adds aceleration.
               const float p_vpm1=-prs*mass2;
                ace.x+=p_vpm1*frx;  ace.y+=p_vpm1*fry;  ace.z+=p_vpm1*frz;
              }
            }
          }
        }
      }
      //-Saves ace.
      partace[idp[p]-idbegin]=ace;
    }
  }
}

//==============================================================================
/// Calculates force on selected fixed or moving particles using only fluid particles.
/// Ignores periodic boundary particles to avoid race condition problems.
//==============================================================================
void Interaction_GaugeForce(unsigned n,unsigned idbegin,typecode codesel
  ,float fourh2,float h,float bwen,float massf,float cteb,float rhopzero,float gamma
  ,int hdiv,tuint3 ncells,tuint3 cellmin,const int2 *begincell,tdouble3 domposmin,float scell
//  ,int hdiv,int4 nc,int3 cellzero,unsigned cellfluid,const int2 *begincell,double3 domposmin,float scell
  ,const double2 *posxy,const double *posz,const typecode *code,const unsigned *idp,const float4 *velrhop
  ,float3 *partace)
{
  const int4 nc=make_int4(ncells.x,ncells.y,ncells.z,ncells.x*ncells.y);
  const unsigned cellfluid=nc.w*nc.z+1;
  const int3 cellzero=make_int3(cellmin.x,cellmin.y,cellmin.z);
  //-Interaction bound-Fluid.
  if(n){
    const unsigned bsize=128;
    dim3 sgrid=GetGridSize(n,bsize);
    //:JDgKerPrint info;
    //:byte* ik=NULL; //info.GetInfoPointer(sgridf,bsfluid);
    KerInteractionGaugeForce <<<sgrid,bsize>>> (n,idbegin,codesel,fourh2,h,bwen,massf,cteb,rhopzero,gamma
      ,hdiv,nc,cellzero,cellfluid,begincell,Double3(domposmin),scell
      ,posxy,posz,code,idp,velrhop,partace);
    //:info.PrintValuesFull(true); //info.PrintValuesInfo();
  }
}


}


