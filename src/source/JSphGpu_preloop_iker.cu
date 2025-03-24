//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpu_shift_iker.cu \brief Implements functions and CUDA kernels for Shifting improved.

#include "JSphGpu_preloop_iker.h"

namespace cusph{

//------------------------------------------------------------------------------
/// Saves particle index to access to the parent of periodic particles.
//------------------------------------------------------------------------------
__global__ void KerPeriodicSaveParent(unsigned n,unsigned pini
  ,const unsigned* listp,unsigned* periparent)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    periparent[pnew]=pcopy;
  }
}
//==============================================================================
/// Saves particle index to access to the parent of periodic particles.
//==============================================================================
void PeriodicSaveParent(unsigned n,unsigned pini,const unsigned* listp
  ,unsigned* periparent)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerPeriodicSaveParent <<<sgrid,SPHBSIZE>>> (n,pini,listp,periparent);
  }
}

//------------------------------------------------------------------------------
/// Updates PreLoop variables in periodic particles.
//------------------------------------------------------------------------------
 __global__ void KerPeriPreLoopCorr(unsigned n,unsigned pinit
  ,const unsigned* periparent,unsigned* fstype,float4* shiftvel)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pinit;//-Number of particle.
    const unsigned pp=periparent[p1];
    if(pp!=UINT_MAX)fstype[p]=fstype[pp];
    if(pp!=UINT_MAX)shiftvel[p]=shiftvel[pp];
  }
}
//==============================================================================
/// Updates PreLoop variables in periodic particles.
//==============================================================================
void PeriPreLoopCorr(unsigned n,unsigned pinit
  ,const unsigned* periparent,unsigned* fstype,float4* shiftvel)
{
  if(n){
    const unsigned bsize=256;
    dim3 sgrid=GetSimpleGridSize(n,bsize);
    //-Calculates Kernel Gradient Correction.
    KerPeriPreLoopCorr <<<sgrid,bsize>>> (n,pinit,periparent,fstype,shiftvel);
  }
}

//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles. (Fluid/Float-Fluid/Float/Bound)
/// Realiza la interaccion de una particula con un conjunto de ellas. (Fluid/Float-Fluid/Float/Bound)
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ void KerComputeNormalsBox(bool boundp2
  ,unsigned p1,const unsigned& pini,const unsigned& pfin,const float4* poscell
  ,const float4* velrhop,const typecode* code,float massp2,const float4& pscellp1
  ,float& fs_treshold,float3& gradc,tmatrix3f& lcorr,unsigned& neigh,float& pou
  ,const float* ftomassp)
{
  const float w0=cufsph::GetKernel_Wab<tker>(CTE.dp*CTE.dp);
  for(int p2=pini;p2<pfin;p2++){
  const float4 pscellp2=poscell[p2];
    const float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
    const float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
    const float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
    const double rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
      // float4 velrhop2=velrhop[p2];
      const float rhopp2= float(velrhop[p2].w);
      //-Velocity derivative (Momentum equation).

      bool ftp2;
      float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
      const typecode cod=code[p2];
      ftp2=CODE_IsFloating(cod);
      ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);

      const float vol2=(ftp2 ? float(ftmassp2/rhopp2) : float(massp2/rhopp2));
      neigh++;

      const float dot3=drx*frx+dry*fry+drz*frz;
      gradc.x+=vol2*frx;
      gradc.y+=vol2*fry;
      gradc.z+=vol2*frz;

      fs_treshold-=vol2*dot3;
      lcorr.a11+=-drx*frx*vol2; lcorr.a12+=-drx*fry*vol2; lcorr.a13+=-drx*frz*vol2;
      lcorr.a21+=-dry*frx*vol2; lcorr.a22+=-dry*fry*vol2; lcorr.a23+=-dry*frz*vol2;
      lcorr.a31+=-drz*frx*vol2; lcorr.a32+=-drz*fry*vol2; lcorr.a33+=-drz*frz*vol2;

      const float wab=cufsph::GetKernel_Wab<tker>(rr2);
      pou+=wab*vol2;       
    }
  }
}

//==============================================================================
/// Perform interaction between particles: Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Realiza interaccion entre particulas: Fluid/Float-Fluid/Float or Fluid/Float-Bound
//==============================================================================
template<TpKernel tker> __global__ void KerComputeNormals(unsigned n,unsigned pinit
  ,int scelldiv,int4 nc,int3 cellzero,const int2* begincell,unsigned cellfluid
  ,const unsigned* dcell,const float4* poscell,const float4* velrhop
  ,const typecode* code,unsigned* fstype,float3* fsnormal,bool simulate2d
  ,float4* shiftvel,const float* ftomassp,const unsigned* listp)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=listp[p];   
    //-Obtains basic data of particle p1.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrhop[p1];
      
    float fs_treshold=0;                              //-Divergence of the position.
    float3 gradc=make_float3(0,0,0);                  //-Gradient of the concentration
    float pou=0;                                      //-Partition of unity.
    unsigned neigh=0;                                 //-Number of neighbours.
      
    tmatrix3f lcorr; cumath::Tmatrix3fReset(lcorr);             //-Correction matrix.
    tmatrix3f lcorr_inv; cumath::Tmatrix3fReset(lcorr_inv);     //-Inverse of the correction matrix.

    //-Calculate approx. number of neighbours when uniform distribution (in the future on the constant memory?)
    float Nzero=0;
    if(simulate2d)Nzero=(3.141592)*CTE.kernelsize2/(CTE.dp*CTE.dp);
    else          Nzero=(4.f/3.f)*(3.141592)*CTE.kernelsize2*CTE.kernelsize/(CTE.dp*CTE.dp*CTE.dp);

    //-Copy the value of shift to gradc. For single resolution is zero, but in Vres take in account virtual stencil.
    gradc=make_float3(shiftvel[p1].x,shiftvel[p1].y,shiftvel[p1].z); 
    
    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
        KerComputeNormalsBox<tker>(false,p1,pini,pfin,poscell,velrhop,code,CTE.massf
          ,pscellp1,fs_treshold,gradc,lcorr,neigh,pou,ftomassp);
      }
    }

    //-Interaction with boundaries.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
        KerComputeNormalsBox<tker>(false,p1,pini,pfin,poscell,velrhop,code,CTE.massb
          ,pscellp1,fs_treshold,gradc,lcorr,neigh,pou,ftomassp);
      }
    }
    //-Find particles that are probably on the free-surface.
    unsigned fstypep1=0;
    if(simulate2d){
      if(fs_treshold<1.7)fstypep1=2;
      if(fs_treshold<1.1 && Nzero/float(neigh)<0.4f)fstypep1=3;
    }
    else{
      if(fs_treshold<2.75)fstypep1=2;
      if(fs_treshold<1.8 && Nzero/float(neigh)<0.4f)fstypep1=3;
    }
    fstype[p1]=fstypep1;

    //-Add the contribution of the particle itself
    pou+=cufsph::GetKernel_Wab<tker>(0.f)*CTE.massf/velrhop1.w;

    //-Calculation of the inverse of the correction matrix (Don't think there is a better way, create function for Determinant2x2 for clarity?).
    if(simulate2d){
      tmatrix2f lcorr2d;
      tmatrix2f lcorr2d_inv;
      lcorr2d.a11=lcorr.a11; lcorr2d.a12=lcorr.a13;
      lcorr2d.a21=lcorr.a31; lcorr2d.a22=lcorr.a33;
      float lcorr_det=(lcorr2d.a11*lcorr2d.a22-lcorr2d.a12*lcorr2d.a21);
      lcorr2d_inv.a11=lcorr2d.a22/lcorr_det; lcorr2d_inv.a12=-lcorr2d.a12/lcorr_det;
      lcorr2d_inv.a22=lcorr2d.a11/lcorr_det; lcorr2d_inv.a21=-lcorr2d.a21/lcorr_det;
      lcorr_inv.a11=lcorr2d_inv.a11;  lcorr_inv.a13=lcorr2d_inv.a12;
      lcorr_inv.a31=lcorr2d_inv.a21;  lcorr_inv.a33=lcorr2d_inv.a22;
    }
    else{
      const float determ=cumath::Determinant3x3(lcorr);
      lcorr_inv=cumath::InverseMatrix3x3(lcorr,determ);
    }

    //-Correction of the gradient of concentration.
    float3 gradc1=make_float3(0,0,0);    
    gradc1.x=gradc.x*lcorr_inv.a11+gradc.y*lcorr_inv.a12+gradc.z*lcorr_inv.a13;
    gradc1.y=gradc.x*lcorr_inv.a21+gradc.y*lcorr_inv.a22+gradc.z*lcorr_inv.a23;
    gradc1.z=gradc.x*lcorr_inv.a31+gradc.y*lcorr_inv.a32+gradc.z*lcorr_inv.a33;    
    float gradc_norm=sqrt(gradc1.x*gradc1.x+gradc1.y*gradc1.y+gradc1.z*gradc1.z);
    //-Set normal.
    fsnormal[p1].x=-gradc1.x/gradc_norm;
    fsnormal[p1].y=-gradc1.y/gradc_norm;
    fsnormal[p1].z=-gradc1.z/gradc_norm;
  }
}

 
//------------------------------------------------------------------------------
/// Obtain the list of particle that are probably on the free-surface.
//------------------------------------------------------------------------------
__global__ void KerCountFreeSurface(unsigned n,unsigned pini
  ,const unsigned* fstype,unsigned* listp)
{
  extern __shared__ unsigned slist[];
  if(!threadIdx.x)slist[0]=0;
  __syncthreads();
  const unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    if(fstype[p]>1) slist[atomicAdd(slist,1)+1]=p;
  }
  __syncthreads();
  const unsigned ns=slist[0];
  __syncthreads();
  if(!threadIdx.x && ns)slist[0]=atomicAdd((listp+n),ns);
  __syncthreads();
  if(threadIdx.x<ns){
    const unsigned cp=slist[0]+threadIdx.x;
    listp[cp]=slist[threadIdx.x+1];
  }
}

//==============================================================================
/// Compute free-surface particles and their normals.
//==============================================================================
void ComputeFSNormals(TpKernel tkernel,bool simulate2d,unsigned bsfluid
  ,unsigned fluidini,unsigned fluidnum,StDivDataGpu& dvd,const unsigned* dcell
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const float4* velrho,const typecode* code,const float* ftomassp
  ,float4* shiftvel,unsigned* fstype,float3* fsnormal,unsigned* listp
  ,cudaStream_t stm)
{
  //-Creates list with free-surface particle (normal and periodic).
  unsigned count=0;
  if(fluidnum){
    cudaMemset(listp+fluidnum,0,sizeof(unsigned));
    dim3 sgridf=GetSimpleGridSize(fluidnum,bsfluid);
    const unsigned smem=(bsfluid+1)*sizeof(unsigned); 
    KerCountFreeSurface <<<sgridf,bsfluid,smem,stm>>> (fluidnum,fluidini,fstype,listp);
  }
  cudaMemcpy(&count,listp+fluidnum,sizeof(unsigned),cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();

  if(count){
    dim3 sgridf=GetSimpleGridSize(count,bsfluid);
    switch(tkernel){
      case KERNEL_Wendland:{ const TpKernel tker=KERNEL_Wendland;
        KerComputeNormals <tker> <<<sgridf,bsfluid,0,stm>>> (count,fluidini
          ,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,dcell
          ,poscell,velrho,code,fstype,fsnormal,simulate2d,shiftvel,ftomassp,listp);
      }break;
     #ifndef DISABLE_KERNELS_EXTRA
      case KERNEL_Cubic:{ const TpKernel tker=KERNEL_Cubic;
        KerComputeNormals <tker> <<<sgridf,bsfluid,0,stm>>> (count,fluidini
          ,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,dcell
          ,poscell,velrho,code,fstype,fsnormal,simulate2d,shiftvel,ftomassp,listp);
      }break;
     #endif
      default: throw "Kernel unknown at ComputeFSNormals().";
    }
  }
  cudaDeviceSynchronize();
}

//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles. (Fluid/Float-Fluid/Float/Bound)
/// Realiza la interaccion de una particula con un conjunto de ellas. (Fluid/Float-Fluid/Float/Bound)
//------------------------------------------------------------------------------
__device__ void KerScanUmbrellaRegionBox(bool boundp2,unsigned p1
  ,const unsigned& pini,const unsigned& pfin,const float4* poscell
  ,const float4& pscellp1,bool& fs_flag,const float3* fsnormal,bool simulate2d)
{
  for(int p2=pini;p2<pfin;p2++){
    const float4 pscellp2=poscell[p2];
    const float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
    const float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
    const float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
    const double rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
      const float3 posq=make_float3(fsnormal[p1].x*CTE.kernelh,fsnormal[p1].y*CTE.kernelh,fsnormal[p1].z*CTE.kernelh);
      if(rr2>2.f*CTE.kernelh*CTE.kernelh){
        const float drxq=-drx-posq.x;
        const float dryq=-dry-posq.y;
        const float drzq=-drz-posq.z;
        const float rrq=sqrt(drxq*drxq+dryq*dryq+drzq*drzq);
        if(rrq<CTE.kernelh)fs_flag=true;
      }
      else{
        if(simulate2d){
          const float drxq=-drx-posq.x;
          const float drzq=-drz-posq.z;
          const float3 normalq=make_float3(drxq*fsnormal[p1].x,0,drzq*fsnormal[p1].z);
          const float3 tangq=make_float3(-drxq*fsnormal[p1].z,0,drzq*fsnormal[p1].x);
          const float normalqnorm=sqrt(normalq.x*normalq.x+normalq.z*normalq.z);
          const float tangqnorm=sqrt(tangq.x*tangq.x+tangq.z*tangq.z);
          if(normalqnorm+tangqnorm<CTE.kernelh) fs_flag=true;
        }
        else{
          float rrr=1.f/sqrt(rr2);
          const float arccosin=acos((-drx*fsnormal[p1].x*rrr-dry*fsnormal[p1].y*rrr-drz*fsnormal[p1].z*rrr));
          if(arccosin<0.785398f)fs_flag=true;
        }
      }
    }
  }
}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid.
/// Interaccion Fluid-Fluid/Bound & Bound-Fluid.
//==============================================================================
__global__ void KerScanUmbrellaRegion(unsigned n,unsigned pinit
  ,int scelldiv,int4 nc,int3 cellzero,const int2* begincell,unsigned cellfluid
  ,const unsigned* dcell,const float4 *poscell,const typecode* code
  ,bool simulate2d,const float3* fsnormal,const unsigned* listp,unsigned* fstype)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=listp[p];   
    bool fs_flag=false;
    //-Obtains basic data of particle p1.
    const float4 pscellp1=poscell[p1];

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin)
        KerScanUmbrellaRegionBox(false,p1,pini,pfin,poscell,pscellp1,fs_flag,fsnormal,simulate2d);
    }

    //-Interaction with boundaries.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin)
        KerScanUmbrellaRegionBox(true,p1,pini,pfin,poscell,pscellp1,fs_flag,fsnormal,simulate2d);
    }
    //-If particle was present in umbrella region, change the code of the particle.
    if(fs_flag && fstype[p1]==2)fstype[p1]=0;
    //-Periodic particle are internal by default.
    if(CODE_IsPeriodic(code[p1]))fstype[p1]=0;
  }
}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
//==============================================================================
void ComputeUmbrellaRegion(TpKernel tkernel,bool simulate2d
  ,unsigned bsfluid,unsigned fluidini,unsigned fluidnum,StDivDataGpu& dvd
  ,const unsigned* dcell,const float4* poscell,const typecode* code
  ,const float3* fsnormal,unsigned* listp,unsigned* fstype,cudaStream_t stm)
{
  //-Obtain the list of particle that are probably on the free-surface (in ComputeUmbrellaRegion maybe is unnecessary).
  unsigned count=0;
  if(fluidnum){
    cudaMemset(listp+fluidnum,0,sizeof(unsigned));
    dim3 sgridf=GetSimpleGridSize(fluidnum,bsfluid);
    const unsigned smem=(bsfluid+1)*sizeof(unsigned); 
    KerCountFreeSurface <<<sgridf,bsfluid,smem,stm>>> (fluidnum,fluidini,fstype,listp);
  }
  cudaMemcpy(&count,listp+fluidnum,sizeof(unsigned),cudaMemcpyDeviceToHost);
  cudaDeviceSynchronize();
  //-Scan Umbrella region on selected free-surface particles.
  if(count){
    dim3 sgridf=GetSimpleGridSize(count,bsfluid);
    KerScanUmbrellaRegion <<<sgridf,bsfluid,0,stm>>> (count,fluidini
      ,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid
      ,dcell,poscell,code,simulate2d,fsnormal,listp,fstype);
  }
  cudaDeviceSynchronize();
}

//==============================================================================
/// Interaction of a particle with a set of particles. (Fluid/Float-Fluid/Float/Bound)
//==============================================================================
template<TpKernel tker,bool simulate2d,bool shiftadv>
__device__ void KerPreLoopInteractionBox(bool boundp2,unsigned p1
  ,const unsigned& pini,const unsigned& pfin,const float4* poscell
  ,const float4* velrhop,const typecode* code,float massp2,const float4& pscellp1
  ,const float4& velrhop1,const float* ftomassp,float4& shiftposf1,unsigned* fs
  ,float3* fsnormal,bool& nearfs,float& mindist,float& maxarccos,bool& bound_inter
  ,float3& fsnormalp1,float& pou)
{
  for(int p2=pini;p2<pfin;p2++){
    const float4 pscellp2=poscell[p2];
    const float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
    const float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
    const float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
    const double rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
      float4 velrhop2=velrhop[p2];

      bool ftp2;
      float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
      const typecode cod=code[p2];
      ftp2=CODE_IsFloating(cod);
      ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);
        
      if(shiftadv){
        const float massrho=(boundp2 ? CTE.massb/velrhop2.w : (ftmassp2)/velrhop2.w);

        //-Compute gradient of concentration and partition of unity.        
        shiftposf1.x+=massrho*frx;    
        shiftposf1.y+=massrho*fry;
        shiftposf1.z+=massrho*frz;          
        const float wab=cufsph::GetKernel_Wab<tker>(rr2);
        shiftposf1.w+=wab*massrho;

        //-Check if the particle is too close to solid or floating object
        if((boundp2 || ftp2)) bound_inter=true;

        //-Check if it close to the free-surface, calculate distance from free-surface and smoothing of free-surface normals.
        if(fs[p2]>1 && fs[p2]<3 && !boundp2 ) {
          nearfs=true;
          mindist=min(sqrt(rr2),mindist);
          pou+=wab*massrho;
          fsnormalp1.x+=fsnormal[p2].x*wab*massrho;
          fsnormalp1.y+=fsnormal[p2].y*wab*massrho;
          fsnormalp1.z+=fsnormal[p2].z*wab*massrho;
        }
        //-Check maximum curvature.
        if(fs[p1]>1 && fs[p2]>1){
          const float norm1=sqrt(fsnormal[p1].x*fsnormal[p1].x+fsnormal[p1].y*fsnormal[p1].y+fsnormal[p1].z*fsnormal[p1].z);
          const float norm2=sqrt(fsnormal[p2].x*fsnormal[p2].x+fsnormal[p2].y*fsnormal[p2].y+fsnormal[p2].z*fsnormal[p2].z);
          maxarccos=max(maxarccos,(acos((fsnormal[p1].x*fsnormal[p2].x+fsnormal[p2].y*fsnormal[p1].y+fsnormal[p2].z*fsnormal[p1].z))));
        }
      }
    }
  }
}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid for models before
/// InteractionForces
//==============================================================================
template<TpKernel tker,bool simulate2d,bool shiftadv>
__global__ void KerPreLoopInteraction(unsigned n,unsigned pinit,int scelldiv
  ,int4 nc,int3 cellzero,const int2* begincell,unsigned cellfluid
  ,const unsigned* dcell,const float4* poscell,const float4* velrhop
  ,const typecode* code,const float* ftomassp,float4* shiftvel,unsigned* fstype
  ,float3* fsnormal,float* fsmindist)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pinit;      //-Number of particle.
    //-Obtains basic data of particle p1.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrhop[p1];
    const float pressp1=cufsph::ComputePressCte(velrhop1.w);
    bool    nearfs=false;                     //-Bool for detecting near free-surface particles. <shiftImproved>
    float4  shiftposp1=make_float4(0,0,0,0);
      
    float mindist=CTE.kernelh;                //-Set Min Distance from free-surface to kernel radius. <shiftImproved>
    float maxarccos=0.0;                      //-Variable for identify high-curvature free-surface particle <shiftImproved>
    bool bound_inter=false;                   //-Variable for identify free-surface that interact with boundary <shiftImproved>
    float3 fsnormalp1=make_float3(0,0,0);     //-Normals for near free-surface particles <shiftImproved>
    unsigned fsp1=fstype[p1];                 //-Free-surface identification code: 0-internal, 1-close to free-surface, 2 free-surface, 3-isolated.
    float   pou=false;                        //-Partition of unity for normal correction.                      <ShiftingAdvanced>
    
    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
        KerPreLoopInteractionBox<tker,simulate2d,shiftadv> (false,p1,pini,pfin
          ,poscell,velrhop,code,CTE.massf,pscellp1,velrhop1,ftomassp,shiftposp1
          ,fstype,fsnormal,nearfs,mindist,maxarccos,bound_inter,fsnormalp1,pou);
      } 
    }

    //-Interaction with bound.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
        KerPreLoopInteractionBox<tker,simulate2d,shiftadv> (true,p1,pini,pfin
          ,poscell,velrhop,code,CTE.massb,pscellp1,velrhop1,ftomassp,shiftposp1
          ,fstype,fsnormal,nearfs,mindist,maxarccos,bound_inter,fsnormalp1,pou);
      }
    }

    if(shiftadv){
      shiftposp1.w+=cufsph::GetKernel_Wab<tker>(0.0)*CTE.massf/velrhop1.w;
      fsmindist[p1]=mindist;
      //-Assign correct code to near free-surface particle and correct their normals by Shepard's Correction.
      if(fsp1==0 && nearfs){
        if(pou>1e-6){
          fsnormalp1=make_float3(fsnormalp1.x,fsnormalp1.y,fsnormalp1.z);
          float norm=sqrt(fsnormalp1.x*fsnormalp1.x+fsnormalp1.y*fsnormalp1.y+fsnormalp1.z*fsnormalp1.z);
          fsnormal[p1]=make_float3(fsnormalp1.x/norm,fsnormalp1.y/norm,fsnormalp1.z/norm);
        }      
        fstype[p1]=1;
        if(bound_inter) fstype[p1]=3;
      }
      //-Check if free-surface particle interact with bound or has high-curvature.
      if(fsp1==2 && (bound_inter||maxarccos>0.52)) shiftposp1=make_float4(0,0,0,shiftposp1.w);
      //-Compute shifting when <shiftImproved> true
      shiftvel[p1]=shiftposp1;
    }
  }
}
  
//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid for models before
/// InteractionForces
//==============================================================================
template<TpKernel tker,bool simulate2d,bool shiftadv> void PreLoopInteractionT3
  (unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd
  ,const unsigned* dcell,const float4* poscell,const float4* velrho
  ,const typecode* code,const float* ftomassp,float4* shiftvel,unsigned* fstype
  ,float3* fsnormal,float* fsmindist,cudaStream_t stm)
{
  if(fluidnum){
    dim3 sgridf=GetSimpleGridSize(fluidnum,bsfluid);
    KerPreLoopInteraction <tker,simulate2d,shiftadv> <<<sgridf,bsfluid,0,stm>>> 
      (fluidnum,fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid
      ,dcell,poscell,velrho,code,ftomassp,shiftvel,fstype,fsnormal,fsmindist);
  }
}
//==============================================================================
template<TpKernel tker,bool simulate2d> void PreLoopInteractionT2(bool shiftadv
  ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd
  ,const unsigned* dcell,const float4* poscell,const float4* velrho
  ,const typecode* code,const float* ftomassp,float4* shiftvel,unsigned* fstype
  ,float3* fsnormal,float* fsmindist,cudaStream_t stm)
{
  if(shiftadv){
    PreLoopInteractionT3 <tker,simulate2d,true > (bsfluid
        ,fluidnum,fluidini,dvd,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,stm);
  }
  else{
    PreLoopInteractionT3 <tker,simulate2d,false> (bsfluid
        ,fluidnum,fluidini,dvd,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,stm);
  }
}
//==============================================================================
template<TpKernel tker> void PreLoopInteractionT(bool simulate2d,bool shiftadv
  ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd
  ,const unsigned* dcell,const float4* poscell,const float4* velrho
  ,const typecode* code,const float* ftomassp,float4* shiftvel,unsigned* fstype
  ,float3* fsnormal,float* fsmindist,cudaStream_t stm)
{
  if(simulate2d){
    PreLoopInteractionT2 <tker,true > (shiftadv,bsfluid
        ,fluidnum,fluidini,dvd,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,stm);
  }
  else{
    PreLoopInteractionT2 <tker,false> (shiftadv,bsfluid
        ,fluidnum,fluidini,dvd,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,stm);
  }
}
//==============================================================================
void PreLoopInteraction(TpKernel tkernel,bool simulate2d,bool shiftadv
  ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd
  ,const unsigned* dcell,const float4* poscell,const float4* velrho
  ,const typecode* code,const float* ftomassp,float4* shiftvel,unsigned* fstype
  ,float3* fsnormal,float* fsmindist,cudaStream_t stm)
{
  switch(tkernel){
    case KERNEL_Wendland:{ const TpKernel tker=KERNEL_Wendland;
      PreLoopInteractionT <tker> (simulate2d,shiftadv,bsfluid
        ,fluidnum,fluidini,dvd,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,stm);
    }break;
   #ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:{ const TpKernel tker=KERNEL_Cubic;
      PreLoopInteractionT <tker> (simulate2d,shiftadv,bsfluid
        ,fluidnum,fluidini,dvd,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,stm);
    }break;
   #endif
    default: throw "Kernel unknown at PreLoopInteraction().";
  }
  cudaDeviceSynchronize();
}

//------------------------------------------------------------------------------
/// Compute shifting velocity for advanced shifting model.
//------------------------------------------------------------------------------
__global__ void KerComputeShiftingVel(unsigned n,unsigned pinit,bool sim2d
  ,float shiftcoef,bool ale,float dt,const unsigned* fstype
  ,const float3* fsnormal,const float* fsmindist,float4* shiftvel)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=p+pinit;      //-Number of particle.
    //-Obtains basic data of particle p1.
    const unsigned fstypep1=fstype[p1];
    const float4   shiftp1=shiftvel[p1];
    const float    fsmindistp1=fsmindist[p1];
    const float3   fsnormalp1=fsnormal[p1];
    const float    theta=min(1.f,max(0.f,(fsmindistp1-CTE.kernelsize)/(0.5*CTE.kernelsize-CTE.kernelsize)));
    float4         shift_final=make_float4(0,0,0,0);
    const float    normshift=(fsnormalp1.x*shiftp1.x+fsnormalp1.y*shiftp1.y+fsnormalp1.z*shiftp1.z);

    if(fstypep1==0){
      shift_final=shiftp1;
    }
    else if(fstypep1==1 || fstypep1==2){
      if(ale){
        shift_final.x=shiftp1.x-theta*fsnormalp1.x*normshift;
        shift_final.y=shiftp1.y-theta*fsnormalp1.y*normshift;
        shift_final.z=shiftp1.z-theta*fsnormalp1.z*normshift;
      }
      else{
        shift_final=make_float4(0,0,0,0);
      }
    }
    else if(fstypep1==3){
      shift_final=make_float4(0,0,0,0);
    }

    if(sim2d)shift_final.y=0.f;

    const float rhovar=abs(shiftp1.x*shift_final.x+shiftp1.y*shift_final.y+shiftp1.z*shift_final.z)*min(CTE.kernelh,fsmindistp1);
    const float eps=1e-5f;
    const float umagn_1=shiftcoef*CTE.kernelh/dt;
    const float umagn_2=abs(eps/(2.f*dt*rhovar));
    const float umagn=min(umagn_1,umagn_2)*(min(CTE.kernelh,fsmindistp1)*dt);
    // const float umagn= shiftcoef*CTE.kernelh*fsmindistp1;
    const float maxdist=0.1f*CTE.dp;
    shift_final.x=(fabs(umagn*shift_final.x)<maxdist? umagn*shift_final.x: (umagn*shift_final.x>=0? maxdist: -maxdist));
    shift_final.y=(fabs(umagn*shift_final.y)<maxdist? umagn*shift_final.y: (umagn*shift_final.y>=0? maxdist: -maxdist));
    shift_final.z=(fabs(umagn*shift_final.z)<maxdist? umagn*shift_final.z: (umagn*shift_final.z>=0? maxdist: -maxdist));
    shiftvel[p1].x=(shift_final.x)/dt;
    shiftvel[p1].y=(shift_final.y)/dt;
    shiftvel[p1].z=(shift_final.z)/dt;
  }
}

//==============================================================================
/// Compute shifting velocity for advanced shifting model.
//==============================================================================
void ComputeShiftingVel(unsigned bsfluid,unsigned fluidnum,unsigned fluidini
  ,bool sim2d,float shiftcoef,bool ale,float dt,const unsigned* fstype
  ,const float3* fsnormal,const float* fsmindist,float4* shiftvel
  ,cudaStream_t stm)
{
  if(fluidnum){
    dim3 sgridf=GetSimpleGridSize(fluidnum,bsfluid);
    KerComputeShiftingVel<<<sgridf,bsfluid,0,stm>>> 
      (fluidnum,fluidini,sim2d,shiftcoef,ale,dt,fstype,fsnormal,fsmindist,shiftvel);
  }
}


}

