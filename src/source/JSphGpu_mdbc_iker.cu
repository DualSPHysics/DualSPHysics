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

/// \file JSphGpu_mdbc_iker.cu \brief Implements functions and CUDA kernels for mDBC.

#include "JSphGpu_mdbc_iker.h"

namespace cusph{

//##############################################################################
//# Kernels for mDBC and mDBC2.
//# Kernels para mDBC y mDBC2.
//##############################################################################
//------------------------------------------------------------------------------
/// Perform interaction between ghost node of selected boundary and fluid.
/// Only for SlipMode==SLIP_Vel0 (DBC vel=0)
//------------------------------------------------------------------------------
template<TpKernel tker,bool sim2d>
  __global__ void KerInteractionMdbcCorrection_Fast(unsigned n,unsigned nbound
  ,double3 mapposmin,float poscellsize,const float4* poscell
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,const double2* posxy,const double* posz,const typecode* code
  ,const unsigned* idp,const float3* boundnor,float4* velrho)
{
  const unsigned p1=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p1<n){
    const float determlimit=1e-3f;
    const float3 bnormalp1=boundnor[p1];
    if(bnormalp1.x!=0 || bnormalp1.y!=0 || bnormalp1.z!=0){
      float rhofinal=FLT_MAX;

      //-Calculates ghost node position.
      double3 gposp1=make_double3(posxy[p1].x+bnormalp1.x,posxy[p1].y+bnormalp1.y,posz[p1]+bnormalp1.z);
      gposp1=(CTE.periactive!=0? KerUpdatePeriodicPos(gposp1): gposp1); //-Corrected interface Position.
      const float4 gpscellp1=KerComputePosCell(gposp1,mapposmin,poscellsize);

      //-Initializes variables for calculation.
      float rhop1=0;
      float3 gradrhop1=make_float3(0,0,0);
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

            //===== Density and its gradient =====
            rhop1+=massp2*wab;
            gradrhop1.x+=massp2*frx;
            gradrhop1.y+=massp2*fry;
            gradrhop1.z+=massp2*frz;

            //===== Kernel values multiplied by volume =====
            const float vwab=wab*volp2;
            const float vfrx=frx*volp2;
            const float vfry=fry*volp2;
            const float vfrz=frz*volp2;

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
      {
        const float3 dpos=make_float3(-bnormalp1.x,-bnormalp1.y,-bnormalp1.z); //-Boundary particle position - ghost node position.
        if(sim2d){
          const double determ=cumath::Determinant3x3dbl(a_corr2);
          if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
            const tmatrix3f invacorr2=cumath::InverseMatrix3x3dbl(a_corr2,determ);
            //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
            const float rhoghost=float(invacorr2.a11*rhop1 + invacorr2.a12*gradrhop1.x + invacorr2.a13*gradrhop1.z);
            const float grx=    -float(invacorr2.a21*rhop1 + invacorr2.a22*gradrhop1.x + invacorr2.a23*gradrhop1.z);
            const float grz=    -float(invacorr2.a31*rhop1 + invacorr2.a32*gradrhop1.x + invacorr2.a33*gradrhop1.z);
            rhofinal=(rhoghost + grx*dpos.x + grz*dpos.z);
          }
          else if(a_corr2.a11>0){//-Determinant is small but a11 is nonzero, 0th order ANGELO.
            rhofinal=float(rhop1/a_corr2.a11);
          }
        }
        else{
          const double determ=cumath::Determinant4x4dbl(a_corr3);
          if(fabs(determ)>=determlimit){
            const tmatrix4f invacorr3=cumath::InverseMatrix4x4dbl(a_corr3,determ);
            //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
            const float rhoghost=float(invacorr3.a11*rhop1 + invacorr3.a12*gradrhop1.x + invacorr3.a13*gradrhop1.y + invacorr3.a14*gradrhop1.z);
            const float grx=    -float(invacorr3.a21*rhop1 + invacorr3.a22*gradrhop1.x + invacorr3.a23*gradrhop1.y + invacorr3.a24*gradrhop1.z);
            const float gry=    -float(invacorr3.a31*rhop1 + invacorr3.a32*gradrhop1.x + invacorr3.a33*gradrhop1.y + invacorr3.a34*gradrhop1.z);
            const float grz=    -float(invacorr3.a41*rhop1 + invacorr3.a42*gradrhop1.x + invacorr3.a43*gradrhop1.y + invacorr3.a44*gradrhop1.z);
            rhofinal=(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
          }
          else if(a_corr3.a11>0){//-Determinant is small but a11 is nonzero, 0th order ANGELO.
            rhofinal=float(rhop1/a_corr3.a11);
          }
        }
        //-Store the results.
        rhofinal=(rhofinal!=FLT_MAX? rhofinal: CTE.rhopzero);
        //-SlipMode==SLIP_Vel0 (DBC vel=0)
        velrho[p1].w=rhofinal;
      }
    }
  }
}


//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
template<TpKernel tker,bool sim2d> void Interaction_MdbcCorrectionT2(unsigned n
  ,unsigned nbound,const StDivDataGpu& dvd,const tdouble3& mapposmin
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const typecode* code,const unsigned* idp,const float3* boundnor
  ,float4* velrho,cudaStream_t stm)
{
  const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(n){
    const unsigned bsbound=128;
    dim3 sgridb=cusph::GetSimpleGridSize(n,bsbound);
    KerInteractionMdbcCorrection_Fast <tker,sim2d> <<<sgridb,bsbound,0,stm>>>
      (n,nbound,Double3(mapposmin),dvd.poscellsize
      ,poscell,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
      ,posxy,posz,code,idp,boundnor,velrho);
  }
}
//==============================================================================
template<TpKernel tker> void Interaction_MdbcCorrectionT(bool simulate2d
  ,unsigned n,unsigned nbound,const StDivDataGpu& dvd,const tdouble3& mapposmin
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const typecode* code,const unsigned* idp,const float3* boundnor
  ,float4* velrho,cudaStream_t stm)
{
  if(simulate2d){
    Interaction_MdbcCorrectionT2 <tker,true > (n,nbound,dvd
      ,mapposmin,posxy,posz,poscell,code,idp,boundnor,velrho,stm);
  }
  else{
    Interaction_MdbcCorrectionT2 <tker,false> (n,nbound,dvd
      ,mapposmin,posxy,posz,poscell,code,idp,boundnor,velrho,stm);
  }
}
//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void Interaction_MdbcCorrection(TpKernel tkernel,bool simulate2d,unsigned n
  ,unsigned nbound,const StDivDataGpu& dvd,const tdouble3& mapposmin
  ,const double2* posxy,const double* posz,const float4* poscell
  ,const typecode* code,const unsigned* idp,const float3* boundnor
  ,float4* velrho,cudaStream_t stm)
{
  switch(tkernel){
    case KERNEL_Wendland:{ const TpKernel tker=KERNEL_Wendland;
      Interaction_MdbcCorrectionT <tker> (simulate2d,n,nbound
        ,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,velrho,stm);
    }break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:{ const TpKernel tker=KERNEL_Cubic;
      Interaction_MdbcCorrectionT <tker> (simulate2d,n,nbound
        ,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,velrho,stm);
    }break;
#endif
    default: throw "Kernel unknown at Interaction_MdbcCorrection().";
  }
}

//<vs_m2dbc_ini>
//##############################################################################
//# Kernels for mDBC2 interaction.
//# Kernels para interaccion mDBC2.
//##############################################################################
//------------------------------------------------------------------------------
/// Perform Pressure cloning for mDBC2.
//------------------------------------------------------------------------------
//  const float pressfinal=KerMdbc2PressClone(sim2d,rhoghost,bnormalp1,gravity,motace,dpos);
__device__ float KerMdbc2PressClone(bool sim2d,const float rhoghost,float3 bnormalp1
  ,const float3 gravity,const float3 motacep1,const float3 dpos)
{
  float pressfinal=0.f;
  if(sim2d){
    const float pghost=float(CTE.cs0*CTE.cs0*(rhoghost-CTE.rhopzero));
    const float norm=sqrt(bnormalp1.x*bnormalp1.x + bnormalp1.z*bnormalp1.z);
    const float normx=bnormalp1.x/norm; 
    const float normz=bnormalp1.z/norm;
    const float normpos=dpos.x*normx + dpos.z*normz;
    const float3 force=make_float3(gravity.x-motacep1.x,0,gravity.z-motacep1.z);
    const float normforce=CTE.rhopzero*(force.x*normx + force.z*normz);
    pressfinal=pghost+normforce*normpos;
  }
  else{
    const float pghost=float(CTE.cs0*CTE.cs0*(rhoghost-CTE.rhopzero));
    const float norm=sqrt(bnormalp1.x*bnormalp1.x + bnormalp1.y*bnormalp1.y + bnormalp1.z*bnormalp1.z);
    const float normx=bnormalp1.x/norm;
    const float normy=bnormalp1.y/norm;
    const float normz=bnormalp1.z/norm;
    const float normpos=dpos.x*normx + dpos.y*normy + dpos.z*normz;
    const float3 force=make_float3(gravity.x-motacep1.x,gravity.y-motacep1.y,gravity.z-motacep1.z);
    const float normforce=CTE.rhopzero*(force.x*normx + force.y*normy + force.z*normz);
    pressfinal=pghost+normforce*normpos;
  }
  return(pressfinal);
}
//------------------------------------------------------------------------------
/// Calculates the infinity norm of a 3x3 matrix.
//------------------------------------------------------------------------------
__device__ float KerMdbc2InfNorm3x3(tmatrix3f mat){
  const float row1=float(fabs(mat.a11) + fabs(mat.a12) + fabs(mat.a13));
  const float row2=float(fabs(mat.a21) + fabs(mat.a22) + fabs(mat.a23));
  const float row3=float(fabs(mat.a31) + fabs(mat.a32) + fabs(mat.a33));
  const float infnorm=max(row1,max(row2,row3));
  return(infnorm);
}
//------------------------------------------------------------------------------
/// Calculates the infinity norm of a 4x4 matrix.
//------------------------------------------------------------------------------
__device__ float KerMdbc2InfNorm4x4(tmatrix4f mat){
  const float row1=float(fabs(mat.a11) + fabs(mat.a12) + fabs(mat.a13) + fabs(mat.a14));
  const float row2=float(fabs(mat.a21) + fabs(mat.a22) + fabs(mat.a23) + fabs(mat.a24));
  const float row3=float(fabs(mat.a31) + fabs(mat.a32) + fabs(mat.a33) + fabs(mat.a34));
  const float row4=float(fabs(mat.a41) + fabs(mat.a42) + fabs(mat.a43) + fabs(mat.a44));
  const float infnorm=max(row1,max(row2,max(row3,row4)));
  return(infnorm);
}
//------------------------------------------------------------------------------
/// Calculates tangent velocity
//------------------------------------------------------------------------------
__device__ float3 KerMdbc2TangenVel(const float3& boundnor,const float3& velfinal){
  const float snormal=sqrt(boundnor.x*boundnor.x + boundnor.y*boundnor.y + boundnor.z*boundnor.z);
  const float bnormalx=boundnor.x/snormal;
  const float bnormaly=boundnor.y/snormal;
  const float bnormalz=boundnor.z/snormal;
	const float veldotnorm=velfinal.x*bnormalx + velfinal.y*bnormaly + velfinal.z*bnormalz;
  const float3 tangentvel=make_float3(velfinal.x-veldotnorm*bnormalx,
                                      velfinal.y-veldotnorm*bnormaly,
                                      velfinal.z-veldotnorm*bnormalz);
  return(tangentvel);
}
//------------------------------------------------------------------------------
/// Perform interaction between ghost node of selected boundary and fluid.
//------------------------------------------------------------------------------
template<TpKernel tker,bool sim2d,TpSlipMode tslip,bool sp>
  __global__ void KerInteractionMdbc2Correction_Fast
  (unsigned n,unsigned nbound,float3 gravity
  ,double3 mapposmin,float poscellsize,const float4* poscell
  ,int scelldiv,int4 nc,int3 cellzero,const int2* beginendcellfluid
  ,const double2* posxy,const double* posz,const typecode* code
  ,const unsigned* idp,const float3* boundnor,const float3* motionvel
  ,const float3* motionace,float4* velrho,byte* boundmode,float3* tangenvel)
{
  const unsigned p1=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p1<n){
    const float3 bnormalp1=boundnor[p1];
    if(bnormalp1.x!=0 || bnormalp1.z!=0 || bnormalp1.y!=0){
      float rhofinal=FLT_MAX;
      float3 velrhofinal=make_float3(0,0,0);
      float sumwab=0;
      float submerged=0;                   //-Not for Vel0.
      const float3 motacep1=motionace[p1]; //-Not for Vel0.

      //-Calculates ghost node position.
      double3 gposp1=make_double3(posxy[p1].x+bnormalp1.x,posxy[p1].y+bnormalp1.y,posz[p1]+bnormalp1.z);
      //-Corrected interface Position.
      gposp1=(CTE.periactive!=0? KerUpdatePeriodicPos(gposp1): gposp1); 
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

            //===== Check if ghost node is submerged =====
            submerged-=volp2*(drx*frx + dry*fry + drz*frz);

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
            velp1.x+=vwab*velrhop2.x;
            velp1.y+=vwab*velrhop2.y;
            velp1.z+=vwab*velrhop2.z;

            //===== Matrix A for correction =====
            if(sim2d){
              a_corr2.a11+=vwab;  a_corr2.a12-=drx*vwab;  a_corr2.a13-=drz*vwab;
              a_corr2.a21+=vfrx;  a_corr2.a22-=drx*vfrx;  a_corr2.a23-=drz*vfrx;
              a_corr2.a31+=vfrz;  a_corr2.a32-=drx*vfrz;  a_corr2.a33-=drz*vfrz;
            }
            else{
              a_corr3.a11+=vwab;  a_corr3.a12-=drx*vwab;  a_corr3.a13-=dry*vwab;  a_corr3.a14-=drz*vwab;
              a_corr3.a21+=vfrx;  a_corr3.a22-=drx*vfrx;  a_corr3.a23-=dry*vfrx;  a_corr3.a24-=drz*vfrx;
              a_corr3.a31+=vfry;  a_corr3.a32-=drx*vfry;  a_corr3.a33-=dry*vfry;  a_corr3.a34-=drz*vfry;
              a_corr3.a41+=vfrz;  a_corr3.a42-=drx*vfrz;  a_corr3.a43-=dry*vfrz;  a_corr3.a44-=drz*vfrz;
            }
          }
        }
      }

      //-Store the results.
      //--------------------
      if(submerged>0.f){
        boundmode[p1]=BMODE_MDBC2;
        const float3 dpos=make_float3(-bnormalp1.x,-bnormalp1.y,-bnormalp1.z); //-Boundary particle position - ghost node position.
        if(sim2d){//-2D simulation.
          if(sumwab<0.1f){ //-If kernel sum is small use shepherd density and vel0.
            //-Trying to avoid too negative pressures.
            rhofinal=max(CTE.rhopzero,(rhop1/a_corr2.a11));
          }
          else{//-Chech if matrix is invertible and well conditioned.
            const double determ=(sp? cumath::Determinant3x3(a_corr2): 
                                     cumath::Determinant3x3dbl(a_corr2));
            if(fabs(determ)>=0.001){//-Use 1e-3f (first_order).
              const tmatrix3f invacorr2=(sp? cumath::InverseMatrix3x3(a_corr2,determ):
                                             cumath::InverseMatrix3x3dbl(a_corr2,determ));
              //-Calculate the scaled condition number
              const float infnorma   =KerMdbc2InfNorm3x3(a_corr2);
              const float infnormainv=KerMdbc2InfNorm3x3(invacorr2);
              const float condinf=CTE.dp*CTE.dp * infnorma * infnormainv;
              if(condinf<=50 && !CODE_IsFlexStrucFlex(code[p1]))//-If matrix is well conditioned use matrix inverse for density and shepherd for velocity.
                rhofinal=float(invacorr2.a11*rhop1 + invacorr2.a12*gradrhop1.x + invacorr2.a13*gradrhop1.z);
              else//-If ill conditioned use shepherd.
                rhofinal=float(rhop1/a_corr2.a11);
            }
            else//-If not invertible use shepherd for density.
              rhofinal=float(rhop1/a_corr2.a11);
          }
          //-Final density according to press.
          const float pressfinal=KerMdbc2PressClone(sim2d,rhofinal,bnormalp1,gravity,motacep1,dpos);
          rhofinal=CTE.rhopzero+float(pressfinal/(CTE.cs0*CTE.cs0));
          //-velocity with Shepherd Sum.
          velrhofinal.x=float(velp1.x/a_corr2.a11);
          velrhofinal.y=0.f;
          velrhofinal.z=float(velp1.z/a_corr2.a11);
        }
        else{//-3D simulation.
          //-Density with pressure cloning.
          if(sumwab<0.1f){//-If kernel sum is small use shepherd for density.
            //-Trying to avoid too negative pressures for empty kernel.
            rhofinal=max(CTE.rhopzero,float(rhop1/a_corr3.a11));
          }
          else{
            const double determ=(sp? cumath::Determinant4x4(a_corr3): 
                                     cumath::Determinant4x4dbl(a_corr3));
            if(fabs(determ)>=0.001){
              const tmatrix4f invacorr3=(sp? cumath::InverseMatrix4x4(a_corr3,determ):
                                             cumath::InverseMatrix4x4dbl(a_corr3,determ));
              //-Calculate the scaled condition number.
              const float infnorma   =KerMdbc2InfNorm4x4(a_corr3);
              const float infnormainv=KerMdbc2InfNorm4x4(invacorr3);
              const float condinf=CTE.dp*CTE.dp*infnorma*infnormainv;
              if(condinf<=50 && !CODE_IsFlexStrucFlex(code[p1]))//-If matrix is well conditioned use matrix inverse for density.
                rhofinal=float(invacorr3.a11*rhop1 + invacorr3.a12*gradrhop1.x + invacorr3.a13*gradrhop1.y + invacorr3.a14*gradrhop1.z);
              else//-Matrix is not well conditioned, use shepherd for density.
                rhofinal=float(rhop1/a_corr3.a11);
            }
            else//-Matrix is not invertible use shepherd for density.
              rhofinal=float(rhop1/a_corr3.a11);
          }
          //-Final density according to press.
          const float pressfinal=KerMdbc2PressClone(sim2d,rhofinal,bnormalp1,gravity,motacep1,dpos);
          rhofinal=CTE.rhopzero+float(pressfinal/(CTE.cs0*CTE.cs0));
          //-Velocity with Shepherd sum.
          velrhofinal.x=float(velp1.x/a_corr3.a11);
          velrhofinal.y=float(velp1.y/a_corr3.a11);
          velrhofinal.z=float(velp1.z/a_corr3.a11);
        }

        //-Store the results.
        if(tslip==SLIP_NoSlip){//-No-Slip: vel = 2*motion - ghost
          const float3 v=motionvel[p1];
          const float3 v2=make_float3(v.x+v.x - velrhofinal.x,
                                      v.y+v.y - velrhofinal.y,
                                      v.z+v.z - velrhofinal.z); 
          #ifndef MDBC2_KEEPVEL
            velrho[p1]=make_float4(v2.x,v2.y,v2.z,rhofinal);
          #else
            velrho[p1].w=rhofinal;
          #endif
          tangenvel[p1]=KerMdbc2TangenVel(bnormalp1,v2);
        }
        else if (tslip==SLIP_FreeSlip) {//-Free slip: vel = ghost vel
          // copy velocity from ghost node
          const float3 v2 = make_float3(velrhofinal.x, velrhofinal.y, velrhofinal.z);
          #ifndef MDBC2_KEEPVEL
            velrho[p1] = make_float4(v2.x, v2.y, v2.z, rhofinal);
          #else
            velrho[p1].w = rhofinal;
          #endif
            tangenvel[p1] = KerMdbc2TangenVel(bnormalp1, v2);
        }
      }
      else{//-If unsubmerged switch off boundary particle.
        boundmode[p1]=BMODE_MDBC2OFF;
        const float3 v=motionvel[p1];
        #ifndef MDBC2_KEEPVEL
          velrho[p1]=make_float4(v.x,v.y,v.z,CTE.rhopzero);
        #else
          velrho[p1].w=CTE.rhopzero;
        #endif
        tangenvel[p1]=KerMdbc2TangenVel(bnormalp1,v);
      }
    }
  }
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
template<TpKernel tker,bool sim2d,TpSlipMode tslip> void Interaction_Mdbc2CorrectionT2
  (unsigned n,unsigned nbound,const tfloat3 gravity
  ,const StDivDataGpu& dvd,const tdouble3& mapposmin,const double2* posxy
  ,const double* posz,const float4* poscell,const typecode* code
  ,const unsigned* idp,const float3* boundnor,const float3* motionvel
  ,const float3* motionace,float4* velrho,byte* boundmode,float3* tangenvel
  ,cudaStream_t stm)
{
  const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(n){
    const unsigned bsbound=128;
    dim3 sgridb=cusph::GetSimpleGridSize(n,bsbound);
    const bool usefloat=false;
    KerInteractionMdbc2Correction_Fast <tker,sim2d,tslip,usefloat> <<<sgridb,bsbound,0,stm>>>
      (n,nbound,Float3(gravity),Double3(mapposmin),dvd.poscellsize
      ,poscell,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid
      ,posxy,posz,code,idp,boundnor,motionvel,motionace
      ,velrho,boundmode,tangenvel);
  }
}
//==============================================================================
template<TpKernel tker> void Interaction_Mdbc2CorrectionT(bool simulate2d
  ,TpSlipMode slipmode,unsigned n,unsigned nbound,const tfloat3 gravity
  ,const StDivDataGpu& dvd,const tdouble3& mapposmin,const double2* posxy
  ,const double* posz,const float4* poscell,const typecode* code
  ,const unsigned* idp,const float3* boundnor,const float3* motionvel
  ,const float3* motionace,float4* velrho,byte* boundmode,float3* tangenvel
  ,cudaStream_t stm)
{
  switch(slipmode){
    case SLIP_NoSlip: { const TpSlipMode tslip = SLIP_NoSlip;
      if (simulate2d) {
          const bool sim2d = true;
          Interaction_Mdbc2CorrectionT2 <tker, sim2d, tslip>(n, nbound, gravity
              , dvd, mapposmin, posxy, posz, poscell, code, idp, boundnor, motionvel
              , motionace, velrho, boundmode, tangenvel, stm);
      }
      else {
          const bool sim2d = false;
          Interaction_Mdbc2CorrectionT2 <tker, sim2d, tslip>(n, nbound, gravity
              , dvd, mapposmin, posxy, posz, poscell, code, idp, boundnor, motionvel
              , motionace, velrho, boundmode, tangenvel, stm);
      }
    }break;
    case SLIP_FreeSlip: { const TpSlipMode tslip = SLIP_FreeSlip;
      if (simulate2d) {
          const bool sim2d = true;
          Interaction_Mdbc2CorrectionT2 <tker, sim2d, tslip>(n, nbound, gravity
              , dvd, mapposmin, posxy, posz, poscell, code, idp, boundnor, motionvel
              , motionace, velrho, boundmode, tangenvel, stm);
      }
      else {
          const bool sim2d = false;
          Interaction_Mdbc2CorrectionT2 <tker, sim2d, tslip>(n, nbound, gravity
              , dvd, mapposmin, posxy, posz, poscell, code, idp, boundnor, motionvel
              , motionace, velrho, boundmode, tangenvel, stm);
      }
    }break;
    default: throw "SlipMode unknown at Interaction_Mdbc2CorrectionT().";
  }
}
//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void Interaction_Mdbc2Correction(TpKernel tkernel,bool simulate2d
  ,TpSlipMode slipmode,unsigned n,unsigned nbound,const tfloat3 gravity
  ,const StDivDataGpu& dvd,const tdouble3& mapposmin,const double2* posxy
  ,const double* posz,const float4* poscell,const typecode* code
  ,const unsigned* idp,const float3* boundnor,const float3* motionvel
  ,const float3* motionace,float4* velrho,byte* boundmode,float3* tangenvel
  ,cudaStream_t stm)
{
  switch(tkernel){
    case KERNEL_Wendland:{ const TpKernel tker=KERNEL_Wendland;
      Interaction_Mdbc2CorrectionT <tker> (simulate2d,slipmode,n,nbound,gravity
        ,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,motionace
        ,velrho,boundmode,tangenvel,stm);
    }break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:{ const TpKernel tker=KERNEL_Cubic;
      Interaction_Mdbc2CorrectionT <tker> (simulate2d,slipmode,n,nbound,gravity
        ,dvd,mapposmin,posxy,posz,poscell,code,idp,boundnor,motionvel,motionace
        ,velrho,boundmode,tangenvel,stm);
    }break;
#endif
    default: throw "Kernel unknown at Interaction_Mdbc2Correction().";
  }
}

//------------------------------------------------------------------------------
/// Copy motion velocity and compute acceleration of moving particles.
//------------------------------------------------------------------------------
__global__ void KerCopyMotionVelAce(unsigned n,double dt,const unsigned* ridpmot
  ,const float4* velrho,float3* motionvel,float3* motionace)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const int pid=ridpmot[p];
    if(pid>=0){
      //-Computes acceleration and copy new velocity.
      const float3 mvel0=motionvel[pid];
      const float4 v=velrho[pid];
      motionace[pid]=make_float3(float((double(v.x)-mvel0.x)/dt),
                                 float((double(v.y)-mvel0.y)/dt),
                                 float((double(v.z)-mvel0.z)/dt));
      motionvel[pid]=make_float3(v.x,v.y,v.z);
    }
  }
}

//==============================================================================
/// Copy motion velocity and compute acceleration of moving particles.
//==============================================================================
void CopyMotionVelAce(unsigned nmoving,double dt,const unsigned* ridpmot
  ,const float4* velrho,float3* motionvel,float3* motionace)
{
  dim3 sgrid=GetSimpleGridSize(nmoving,SPHBSIZE);
  KerCopyMotionVelAce <<<sgrid,SPHBSIZE>>> (nmoving,dt,ridpmot,velrho
    ,motionvel,motionace);
}
//<vs_m2dbc_end>


}

