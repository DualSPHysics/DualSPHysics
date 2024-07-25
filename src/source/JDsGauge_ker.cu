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

/// \file JDsGauge_ker.cu \brief Implements functions and CUDA kernels for classes JGauge.

#include "JDsGauge_ker.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include <float.h>
#include <math_constants.h>
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"
#include <cstdio>
#include <string>

//#include "TypesDef.h"
//#include <cuda_runtime_api.h>

namespace cugauge{
#include "FunctionsBasic_iker.h"
#include "FunSphKernel_iker.h"
#include "FunSphEos_iker.h"
#include "JCellSearch_iker.h"

//##############################################################################
//# Kernels for gauge interaction.
//##############################################################################
//------------------------------------------------------------------------------
/// Performs interaction between particles. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// It includes visco artificial/laminar and floatings SPH/DEM.
///
/// Realiza interaccion entre particulas. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings SPH/DEM.
//------------------------------------------------------------------------------
template<TpKernel tker> __global__ void KerInteractionGaugeVel(float aker
  ,double3 ptpos,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,unsigned axis,unsigned cellcode,double3 domposmin,float scell,float poscellsize
  ,float kernelsize2,float kernelh,float massf
  ,const double2* posxy,const double* posz,const typecode* code,const float4* velrhop
  ,float3* ptvel)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle
  if(!p){
    const double px=ptpos.x;
    const double py=ptpos.y;
    const double pz=ptpos.z;

    double sumwab=0;
    double3 sumvel=make_double3(0,0,0);

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::Initsp(px,py,pz,axis,domposmin,scell,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    //ini3+=cellfluid; fin3+=cellfluid; //cellfluid is included in *beginendcellfluid.
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
      if(pfin)for(int p2=pini;p2<pfin;p2++){
        const double2 pxyp2=posxy[p2];
        const float drx=float(px-pxyp2.x);
        const float dry=float(py-pxyp2.y);
        const float drz=float(pz-posz[p2]);
        const float rr2=(drx*drx + dry*dry + drz*drz);
        //-Interaction with real neighboring fluid particles.
        if(rr2<=kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
          float wab=cufsph::GetKernel_Wab<tker>(rr2,kernelh,aker);
          const float4 velrhopp2=velrhop[p2];
          wab*=massf/velrhopp2.w;
          sumwab+=wab;
          sumvel.x+=wab*velrhopp2.x;
          sumvel.y+=wab*velrhopp2.y;
          sumvel.z+=wab*velrhopp2.z;

        }
      }
    }
    //-Applies kernel correction.
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
void Interaction_GaugeVel(const StCteSph& CSP,const StDivDataGpu& dvd
  ,tdouble3 ptpos,const double2* posxy,const double* posz,const typecode* code
  ,const float4* velrhop,float3* ptvel)
  //,tdouble3 domposmin,float scell,float kernelsize2,float h,float massf)
{
  //-Interaction Fluid-Fluid & Fluid-Bound.
  const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
  const unsigned bsize=32;
  dim3 sgrid=GetSimpleGridSize(1,bsize);
  //:JDgKerPrint info;
  //:byte* ik=NULL; //info.GetInfoPointer(sgridf,bsfluid);
  switch(CSP.tkernel){
    case KERNEL_Cubic:   //Kernel Cubic is not available.
    case KERNEL_Wendland:{ const float aker=CSP.kwend.awen;
      KerInteractionGaugeVel<KERNEL_Wendland> <<<sgrid,bsize>>> (aker,Double3(ptpos)
        ,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
        ,dvd.axis,dvd.domcellcode,dvd.domposmin,dvd.scell,dvd.poscellsize
        ,dvd.kernelsize2,CSP.kernelh,CSP.massfluid,posxy,posz,code,velrhop,ptvel);
    }break;
    default: throw "Kernel unknown at Interaction_GaugeVel().";
  }
  //:info.PrintValuesFull(true); //info.PrintValuesInfo();
}

//------------------------------------------------------------------------------
/// Calculates mass value at one point by interacting with the fluid.
/// Calcula valor de masa en un punto mediante la interaccion con el fluido.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float KerCalculeMass(float aker
  ,double px,double py,double pz,float kernelsize2,float kernelh,float massf
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,unsigned axis,unsigned cellcode,double3 domposmin,float scell,float poscellsize
  ,const double2* posxy,const double* posz,const typecode* code
  ,const float4* velrhop)
{
  double summass=0;

  //-Obtains neighborhood search limits.
  int ini1,fin1,ini2,fin2,ini3,fin3;
  cunsearch::Initsp(px,py,pz,axis,domposmin,scell,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

  //-Interaction with fluids.
  for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
    unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
    if(pfin)for(int p2=pini;p2<pfin;p2++){
      const double2 pxyp2=posxy[p2];
      const float drx=float(px-pxyp2.x);
      const float dry=float(py-pxyp2.y);
      const float drz=float(pz-posz[p2]);
      const float rr2=(drx*drx + dry*dry + drz*drz);
      if(rr2<=kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
        float wab=cufsph::GetKernel_Wab<tker>(rr2,kernelh,aker);
        wab*=massf/velrhop[p2].w;
        summass+=wab*massf;
      }
    }
  }
  return(float(summass));
}

//------------------------------------------------------------------------------
/// Calculates surface water level at indicated line.
//------------------------------------------------------------------------------
template<TpKernel tker> __global__ void KerInteractionGaugeSwl(float aker
  ,double p0x,double p0y,double p0z,double pdirx,double pdiry,double pdirz
  ,unsigned pointnp,float masslimit,float kernelsize2,float kernelh,float massf
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,unsigned axis,unsigned cellcode,double3 domposmin,float scell,float poscellsize
  ,const double2* posxy,const double* posz,const typecode* code
  ,const float4* velrhop,float3* ptres)
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
      shmass[tid]=KerCalculeMass<tker>(aker,p0x+pdirx*cp,p0y+pdiry*cp,p0z+pdirz*cp
        ,kernelsize2,kernelh,massf,scelldiv,nc,cellzero,beginendcellfluid
        ,axis,cellcode,domposmin,scell,poscellsize,posxy,posz,code,velrhop);
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
void Interaction_GaugeSwl(const StCteSph& CSP,const StDivDataGpu& dvd
  ,tdouble3 point0,tdouble3 pointdir,unsigned pointnp,float masslimit
  ,const double2* posxy,const double* posz,const typecode* code
  ,const float4* velrhop,float3* ptres)
{
  const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
  const unsigned bsize=128;
  dim3 sgrid=GetSimpleGridSize(bsize,bsize);
  const unsigned smem=sizeof(float)*(bsize+1);
  switch(CSP.tkernel){
    case KERNEL_Cubic:   //Kernel Cubic is not available.
    case KERNEL_Wendland:{ const float aker=CSP.kwend.awen;
      KerInteractionGaugeSwl<KERNEL_Wendland> <<<sgrid,bsize,smem>>> (aker
        ,point0.x,point0.y,point0.z,pointdir.x,pointdir.y,pointdir.z
        ,pointnp,masslimit,dvd.kernelsize2,CSP.kernelh,CSP.massfluid
        ,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
        ,dvd.axis,dvd.domcellcode,dvd.domposmin,dvd.scell,dvd.poscellsize
        ,posxy,posz,code,velrhop,ptres);
    }break;
    default: throw "Kernel unknown at Interaction_GaugeSwl().";
  }
}


//------------------------------------------------------------------------------
/// Calculates maximum z of fluid at distance of a vertical line.
//------------------------------------------------------------------------------
__global__ void KerInteractionGaugeMaxz(double p0x,double p0y,float maxdist2
  ,int cxini,int cxfin,int yini,int yfin,int zini,int zfin
  ,int4 nc,unsigned cellfluid,const int2* begincell
  ,const double2* posxy,const double* posz,const typecode* code
  ,float3* ptres)
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
void Interaction_GaugeMaxz(tdouble3 point0,float maxdist2,const StDivDataGpu& dvd
  ,int cxini,int cxfin,int yini,int yfin,int zini,int zfin
  ,const double2* posxy,const double* posz,const typecode* code
  ,float3* ptres)
{
  const unsigned bsize=128;
  dim3 sgrid=GetSimpleGridSize(1,bsize);
  KerInteractionGaugeMaxz <<<sgrid,bsize>>> (point0.x,point0.y,maxdist2
    ,cxini,cxfin,yini,yfin,zini,zfin,dvd.nc,dvd.cellfluid,dvd.beginendcell
    ,posxy,posz,code,ptres);
}



//<vs_meeshdat_ini>
//------------------------------------------------------------------------------
/// Calculates data interpolation according to fluid particle data.
//------------------------------------------------------------------------------
template<TpKernel tker,bool tvel> __global__ void KerComputeGaugeMesh(unsigned nptot
  ,unsigned npt12,unsigned npt1,double3 ptref
  ,double3 vdp1,double3 vdp2,double3 vdp3,float3 vdir
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,unsigned axis,unsigned cellcode,double3 domposmin,float scell,float poscellsize
  ,float aker,float kernelsize2,float kernelh,float massf,float kclimit,float kcdummy
  ,const double2* posxy,const double* posz,const typecode* code
  ,const float4* velrhop,float3* ptvel,float* ptvdir,float* ptrhop
  ,float* ptmass)
{
  unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of position.
  if(cp<nptot){
    const unsigned cp3=cp/npt12;
    const unsigned cp3r=cp-cp3*npt12;
    const unsigned cp2=cp3r/npt1;
    const unsigned cp1=cp3r-cp2*npt1;
    const double px=ptref.x+(vdp1.x*cp1)+(vdp2.x*cp2)+(vdp3.x*cp3);
    const double py=ptref.y+(vdp1.y*cp1)+(vdp2.y*cp2)+(vdp3.y*cp3);
    const double pz=ptref.z+(vdp1.z*cp1)+(vdp2.z*cp2)+(vdp3.z*cp3);
    float sumwab=0;
    float summass=0;
    float sumrhop=0;
    float3 sumvel=make_float3(0,0,0);

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::Initsp(px,py,pz,axis,domposmin,scell,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    //ini3+=cellfluid; fin3+=cellfluid; //cellfluid is included in *beginendcellfluid.
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
      if(pfin)for(int p2=pini;p2<pfin;p2++){
        const double2 pxyp2=posxy[p2];
        const float drx=float(px-pxyp2.x);
        const float dry=float(py-pxyp2.y);
        const float drz=float(pz-posz[p2]);
        const float rr2=(drx*drx + dry*dry + drz*drz);
        //-Interaction with real neighboring fluid particles.
        if(rr2<=kernelsize2 && CODE_IsFluid(code[p2])){
          float wab=cufsph::GetKernel_Wab<tker>(rr2,kernelh,aker);
          const float4 velrhopp2=velrhop[p2];
          wab*=massf/velrhopp2.w;
          sumwab+=wab;
          summass+=wab*massf;
          sumrhop+=wab*velrhopp2.w;
          if(tvel){
            sumvel.x+=wab*velrhopp2.x;
            sumvel.y+=wab*velrhopp2.y;
            sumvel.z+=wab*velrhopp2.z;
          }
        }
      }
    }
    //-Applies kernel correction.
    if(kclimit!=FLT_MAX){
      if(sumwab>=kclimit){
        sumvel.x/=sumwab;
        sumvel.y/=sumwab;
        sumvel.z/=sumwab;
        sumrhop/=sumwab;
      }
      else if(kcdummy!=FLT_MAX){
        sumvel=make_float3(kcdummy,kcdummy,kcdummy);
        sumrhop=kcdummy;
      }
    }
    //-Stores results.
    if(tvel){
      if(ptvel )ptvel [cp]=sumvel;
      if(ptvdir)ptvdir[cp]=(sumvel.x*vdir.x + sumvel.y*vdir.y + sumvel.z*vdir.z);
    }
    if(ptrhop)ptrhop[cp]=sumrhop;
    if(ptmass)ptmass[cp]=summass;
  }
}

//==============================================================================
/// Calculates data interpolation according to fluid particle data.
//==============================================================================
template<TpKernel tker,bool tvel> void ComputeGaugeMeshT(float aker
  ,const StCteSph& CSP,const StDivDataGpu& dvd,const jmsh::StMeshPts& mp
  ,float kclimit,float kcdummy,const double2* posxy,const double* posz
  ,const typecode* code,const float4* velrhop,float3* ptvel,float* ptvdir
  ,float* ptrhop,float* ptmass)
{
  //-Interaction with fluid particles.
  const unsigned npt12=mp.npt1*mp.npt2;
  const unsigned nptot=npt12*mp.npt3;
  if(nptot){
    const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(nptot,bsize);
    //:JDgKerPrint info;
    //:byte* ik=NULL; //info.GetInfoPointer(sgridf,bsfluid);
    KerComputeGaugeMesh<tker,tvel> <<<sgrid,bsize>>> (nptot,npt12,mp.npt1
     ,Double3(mp.ptref),Double3(mp.vdp1),Double3(mp.vdp2),Double3(mp.vdp3),Float3(mp.dirdat)
     ,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
     ,dvd.axis,dvd.domcellcode,dvd.domposmin,dvd.scell,dvd.poscellsize
     ,aker,dvd.kernelsize2,CSP.kernelh,CSP.massfluid,kclimit,kcdummy
     ,posxy,posz,code,velrhop,ptvel,ptvdir,ptrhop,ptmass);
    //:info.PrintValuesFull(true); //info.PrintValuesInfo();
  }
}

//==============================================================================
/// Calculates velocity and/or mass interpolation according to fluid particle data.
//==============================================================================
void ComputeGaugeMesh(const StCteSph& CSP,const StDivDataGpu& dvd
  ,const jmsh::StMeshPts& mp,float kclimit,float kcdummy,const double2* posxy
  ,const double* posz,const typecode* code,const float4* velrhop
  ,float3* ptvel,float* ptvdir,float* ptrhop,float* ptmass)
{
  const bool tvel=(ptvel!=NULL || ptvdir!=NULL);
  switch(CSP.tkernel){
    case KERNEL_Cubic:   //Kernel Cubic is not available.
    case KERNEL_Wendland:{ const float aker=CSP.kwend.awen;
      if(tvel)ComputeGaugeMeshT<KERNEL_Wendland,true >(aker,CSP,dvd,mp,kclimit,kcdummy,posxy,posz,code,velrhop,ptvel,ptvdir,ptrhop,ptmass);
      else    ComputeGaugeMeshT<KERNEL_Wendland,false>(aker,CSP,dvd,mp,kclimit,kcdummy,posxy,posz,code,velrhop,ptvel,ptvdir,ptrhop,ptmass);
    }break;
    default: throw "Kernel unknown at ComputeGaugeMesh().";
  }
}

//------------------------------------------------------------------------------
/// Calculates Z surface water level according mass interpolation values on fluid.
//------------------------------------------------------------------------------
__global__ void KerComputeGaugeMeshZsurf(unsigned npt12,unsigned npt1,unsigned npt3
  ,float masslimit,double ptrefz,double vdp1z,double vdp2z,double vdp3z
  ,const float* ptmass,float* ptzsurf)
{
  const unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of position.
  if(cp<npt12){
    const unsigned cp2=cp/npt1;
    const unsigned cp1=cp-cp2*npt1;
    float masspre=0;
    unsigned cpsurf=0;
    float    fsurf=0;
    unsigned cpp=cp;
    for(unsigned cp3=0;cp3<npt3 && !cpsurf;cp3++){
      const float mass=ptmass[cpp];
      if(mass>masslimit)masspre=mass;
      if(mass<masslimit && masspre){
        fsurf=(masslimit-masspre)/(mass-masspre);
        cpsurf=cp3;
      }
      cpp+=npt12;
    }
    float zsurf=float(ptrefz+(vdp1z*cp1)+(vdp2z*cp2));    //-Minimum zsurf.
    if(cpsurf==0 && masspre)zsurf+=float(vdp3z*(npt3-1)); //-Maximum zsurf.
    if(cpsurf){
      zsurf+=float((vdp3z*(cpsurf-1))+(vdp3z*fsurf));     //-Found zsurf.
    }
    ptzsurf[cp]=zsurf;
  }
}

//==============================================================================
/// Calculates Z surface water level according mass interpolation values on fluid.
//==============================================================================
void ComputeGaugeMeshZsurf(float masslimit,const jmsh::StMeshPts& mp
  ,const float* ptmass,float* ptzsurf)
{
  //-Interaction with fluid particles.
  const unsigned npt12=mp.npt1*mp.npt2;
  if(npt12){
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(npt12,bsize);
    //:JDgKerPrint info;
    //:byte* ik=NULL; //info.GetInfoPointer(sgridf,bsfluid);
    cudaMemset(ptzsurf,0,sizeof(float)*npt12);
    KerComputeGaugeMeshZsurf <<<sgrid,bsize>>> (npt12,mp.npt1,mp.npt3
      ,masslimit,mp.ptref.z,mp.vdp1.z,mp.vdp2.z,mp.vdp3.z,ptmass,ptzsurf);
    //:info.PrintValuesFull(true); //info.PrintValuesInfo();
  }
}
//<vs_meeshdat_end>


//<vs_flowdat_ini>
//------------------------------------------------------------------------------
/// Calculates data interpolation according to fluid particle data.
//------------------------------------------------------------------------------
template<TpKernel tker,bool tvel> __global__ void KerComputeGaugeFlow(unsigned nptot
  ,unsigned npt12,unsigned npt1,double3 ptref
  ,double3 vdp1,double3 vdp2,double3 vdp3,float3 vdir
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,unsigned axis,unsigned cellcode,double3 domposmin,float scell,float poscellsize
  ,float aker,float kernelsize2,float kernelh,float massf,float kclimit
  ,const double2* posxy,const double* posz,const typecode* code
  ,const float4* velrhop,float* ptvdir)
{
  unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of position.
  if(cp<nptot){
    const unsigned cp3=cp/npt12;
    const unsigned cp3r=cp-cp3*npt12;
    const unsigned cp2=cp3r/npt1;
    const unsigned cp1=cp3r-cp2*npt1;
    const double px=ptref.x+(vdp1.x*cp1)+(vdp2.x*cp2)+(vdp3.x*cp3);
    const double py=ptref.y+(vdp1.y*cp1)+(vdp2.y*cp2)+(vdp3.y*cp3);
    const double pz=ptref.z+(vdp1.z*cp1)+(vdp2.z*cp2)+(vdp3.z*cp3);
    float sumwab=0;
    float3 sumvel=make_float3(0,0,0);

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::Initsp(px,py,pz,axis,domposmin,scell,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    //ini3+=cellfluid; fin3+=cellfluid; //cellfluid is included in *beginendcellfluid.
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
      if(pfin)for(int p2=pini;p2<pfin;p2++){
        const double2 pxyp2=posxy[p2];
        const float drx=float(px-pxyp2.x);
        const float dry=float(py-pxyp2.y);
        const float drz=float(pz-posz[p2]);
        const float rr2=(drx*drx + dry*dry + drz*drz);
        //-Interaction with real neighboring fluid particles.
        if(rr2<=kernelsize2 && CODE_IsFluid(code[p2])){
          float wab=cufsph::GetKernel_Wab<tker>(rr2,kernelh,aker);
          const float4 velrhopp2=velrhop[p2];
          wab*=massf/velrhopp2.w;
          sumwab+=wab;
          if(tvel){
            sumvel.x+=wab*velrhopp2.x;
            sumvel.y+=wab*velrhopp2.y;
            sumvel.z+=wab*velrhopp2.z;
          }
        }
      }
    }
    //-Applies kernel correction.
    if(kclimit!=FLT_MAX){
      if(sumwab>=kclimit){
        sumvel.x/=sumwab;
        sumvel.y/=sumwab;
        sumvel.z/=sumwab;        
      }
      else {
        sumvel=make_float3(0,0,0);        
      }
    }
    //-Stores results.
    if(ptvdir)ptvdir[cp]=(sumvel.x*vdir.x + sumvel.y*vdir.y + sumvel.z*vdir.z);
  }
}

//==============================================================================
/// Calculates data interpolation according to fluid particle data.
//==============================================================================
template<TpKernel tker,bool tvel> void ComputeGaugeFlowT(float aker
  ,const StCteSph& CSP,const StDivDataGpu& dvd,const jmsh::StMeshPts& mp
  ,float kclimit,const double2* posxy,const double* posz
  ,const typecode* code,const float4* velrhop,float* ptvdir)
{
  //-Interaction with fluid particles.
  const unsigned npt12=mp.npt1*mp.npt2;
  const unsigned nptot=npt12*mp.npt3;
  //const unsigned nptot=mp.npt1*mp.npt2;
  if(nptot){
    const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(nptot,bsize);
    //:JDgKerPrint info;
    //:byte* ik=NULL; //info.GetInfoPointer(sgridf,bsfluid);
    KerComputeGaugeFlow<tker,tvel> <<<sgrid,bsize>>> (nptot,npt12,mp.npt1
     ,Double3(mp.ptref),Double3(mp.vdp1),Double3(mp.vdp2),Double3(mp.vdp3),Float3(mp.dirdat)
     ,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
     ,dvd.axis,dvd.domcellcode,dvd.domposmin,dvd.scell,dvd.poscellsize
     ,aker,dvd.kernelsize2,CSP.kernelh,CSP.massfluid,kclimit
     ,posxy,posz,code,velrhop,ptvdir);
    //:info.PrintValuesFull(true); //info.PrintValuesInfo();
  }
}

//==============================================================================
/// Calculates velocity according to fluid particle data.
//==============================================================================
void ComputeGaugeFlow(const StCteSph& CSP,const StDivDataGpu& dvd
  ,const jmsh::StMeshPts& mp,float kclimit,const double2* posxy
  ,const double* posz,const typecode* code,const float4* velrhop
  ,float* ptvdir)
{  
  switch(CSP.tkernel){
    case KERNEL_Cubic:   //Kernel Cubic is not available.
    case KERNEL_Wendland:{ const float aker=CSP.kwend.awen;
      ComputeGaugeFlowT<KERNEL_Wendland,true >(aker,CSP,dvd,mp,kclimit,posxy,posz,code,velrhop,ptvdir);
    }break;
    default: throw "Kernel unknown at ComputeGaugeMesh().";
  }
}
//<vs_flowdat_end>


//------------------------------------------------------------------------------
/// Calculates force on selected fixed or moving particles using only fluid particles.
/// Ignores periodic boundary particles to avoid race condition problems.
//------------------------------------------------------------------------------
template<TpKernel tker> __global__ void KerInteractionGaugeForce(float bhker
  ,unsigned n,unsigned idbegin,typecode codesel
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,unsigned axis,unsigned cellcode,double3 domposmin,float scell,float poscellsize
  ,float kernelsize2,float kernelh,float massf,float cteb,float rhopzero,float gamma,float cs0
  ,const double2* posxy,const double* posz,const typecode* code
  ,const unsigned* idp,const float4* velrhop,float3* partace)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const typecode code1=code[p];
    if(CODE_GetTypeAndValue(code1)==codesel && CODE_IsNormal(code1)){
      const double2 ptposxy=posxy[p];
      const double px=ptposxy.x;
      const double py=ptposxy.y;
      const double pz=posz[p];
      const float rhop1=velrhop[p].w;
      const float press1=cufsph::ComputePress(rhop1,rhopzero,cteb,gamma,cs0);
      float3 ace=make_float3(0,0,0);

      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::Initsp(px,py,pz,axis,domposmin,scell,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);
   
      //-Interaction with fluids.
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
        if(pfin)for(int p2=pini;p2<pfin;p2++){
          double2 pxyp2=posxy[p2];
          const float drx=float(px-pxyp2.x);
          const float dry=float(py-pxyp2.y);
          const float drz=float(pz-posz[p2]);
          const float rr2=(drx*drx + dry*dry + drz*drz);
          //-Interaction with real neighboring fluid particles.
          if(rr2<=kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
            const float fac=cufsph::GetKernel_Fac<tker>(rr2,kernelh,bhker);
            const float frx=fac*drx;
            const float fry=fac*dry;
            const float frz=fac*drz;

            //-Velocity derivative (Momentum equation).
            const float mass2=massf;
            const float rhop2=velrhop[p2].w;
            const float press2=cufsph::ComputePress(rhop2,rhopzero,cteb,gamma,cs0);
            const float prs=(press1+press2)/(rhop1*rhop2);
            {//-Adds aceleration.
              const float p_vpm1=-prs*mass2;
              ace.x+=p_vpm1*frx;  ace.y+=p_vpm1*fry;  ace.z+=p_vpm1*frz;
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
void Interaction_GaugeForce(const StCteSph& CSP,const StDivDataGpu& dvd
  ,unsigned n,unsigned idbegin,typecode codesel,const double2* posxy
  ,const double* posz,const typecode* code,const unsigned* idp
  ,const float4* velrhop,float3* partace)
{
  //const float ovrhopzero=1.f/rhopzero;
  //-Interaction bound-Fluid.
  if(n){
    const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(n,bsize);
    switch(CSP.tkernel){
      case KERNEL_Cubic:   //Kernel Cubic is not available.
      case KERNEL_Wendland:{ const float bhker=CSP.kwend.bwenh;
        KerInteractionGaugeForce<KERNEL_Wendland> <<<sgrid,bsize>>>(bhker,n,idbegin,codesel
          ,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
          ,dvd.axis,dvd.domcellcode,dvd.domposmin,dvd.scell,dvd.poscellsize
          ,dvd.kernelsize2,CSP.kernelh,CSP.massfluid,CSP.cteb,CSP.rhopzero,CSP.gamma,float(CSP.cs0)
          ,posxy,posz,code,idp,velrhop,partace);
      }break;
      default: throw "Kernel unknown at Interaction_GaugeForce().";
    }
  }
}


}


