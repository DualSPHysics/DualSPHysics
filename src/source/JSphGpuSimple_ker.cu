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

/// \file JSphGpuSimple_ker.cu \brief Implements functions and CUDA kernels for the Particle Interaction and System Update.

#include "JSphGpuSimple_ker.h"
//#include "Functions.h"
//#include "FunctionsCuda.h"
//#include <math_constants.h>
//#include "JDgKerPrint.h"
//#include "JDgKerPrint_ker.h"
#include <cfloat>


namespace cusphs{
#include "FunctionsBasic_iker.h"


//##############################################################################
//# Kernels to prepare data before Interaction_Forces().
//##############################################################################
//------------------------------------------------------------------------------
/// Update PosCellg[] according to current position of particles.
/// Actualiza PosCellg[] segun la posicion de las particulas.
//------------------------------------------------------------------------------
__global__ void KerUpdatePosCell(unsigned np,double3 posmin,float poscellsize
  ,const double2* posxy,const double* posz,float4* poscell)
{
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<np){
    const double2 rxy=posxy[p];
    const double dx=rxy.x-posmin.x;
    const double dy=rxy.y-posmin.y;
    const double dz=posz[p]-posmin.z;
    const unsigned cx=unsigned(dx/poscellsize);
    const unsigned cy=unsigned(dy/poscellsize);
    const unsigned cz=unsigned(dz/poscellsize);
    const float px=float(dx-(double(poscellsize)*cx));
    const float py=float(dy-(double(poscellsize)*cy));
    const float pz=float(dz-(double(poscellsize)*cz));
    const float pw=__uint_as_float(PSCEL_Code(cx,cy,cz));
    poscell[p]=make_float4(px,py,pz,pw);
  }
}
//==============================================================================
/// Update PosCellg[] according to current position of particles.
/// Actualiza PosCellg[] segun la posicion de las particulas.
//==============================================================================
void UpdatePosCell(unsigned np,tdouble3 posmin,float poscellsize
  ,const double2* posxy,const double* posz,float4* poscell,cudaStream_t stm)
{
  const dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
  if(np)KerUpdatePosCell <<<sgrid,SPHBSIZE,0,stm>>> (np,Double3(posmin),poscellsize,posxy,posz,poscell);
}

//------------------------------------------------------------------------------
/// Initialises ace array with 0 for bound and gravity for fluid.
/// Inicializa el array ace con 0 para contorno y gravity para fluido.
//------------------------------------------------------------------------------
__global__ void KerInitAceGravity(unsigned np,unsigned npb,float3 gravity
  ,float3* ace)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<np){
    ace[p]=(p<npb? make_float3(0,0,0): gravity);
  }
}
//==============================================================================
/// Initialises ace array with 0 for bound and gravity for fluid.
/// Inicializa el array ace con 0 para contorno y gravity para fluido.
//==============================================================================
void InitAceGravity(unsigned np,unsigned npb,tfloat3 gravity,float3* ace
  ,cudaStream_t stm)
{
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
    KerInitAceGravity <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,Float3(gravity),ace);
  }
}


//##############################################################################
//# Kernels to run after Interaction_Forces().
//##############################################################################
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
void Resety(unsigned n,unsigned ini,float3* v,cudaStream_t stm){
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerResety <<<sgrid,SPHBSIZE,0,stm>>> (n,ini,v);
  }
}


//##############################################################################
//# Kernels for ComputeStep (vel & rho).
//# Kernels para ComputeStep (vel & rho).
//##############################################################################
//------------------------------------------------------------------------------
/// Computes new values for Pos, Check, Vel and Ros (using Verlet).
/// The value of Vel always set to be reset.
///
/// Calcula nuevos valores de  Pos, Check, Vel y Rhop (usando Verlet).
/// El valor de Vel para bound siempre se pone a cero.
//------------------------------------------------------------------------------
template<bool floating,bool shift,bool inout> __global__ void KerComputeStepVerlet
  (unsigned n,unsigned npb,TpMdbc2Mode mdbc2,float rhopzero,float rhopoutmin
  ,float rhopoutmax,const float4* velrho1,const float4* velrho2,const byte* boundmode
  ,const float* ar,const float3* ace,const float4* shiftposfs,const float3* indirvel
  ,const float4* nopenshift,double dt,double dt205,double dt2,float3 gravity
  ,double2* movxy,double* movz,typecode* code,float4* velrhonew)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    if(p<npb){ //-Particles: Fixed & Moving.
      if(mdbc2>=MDBC2_Std && boundmode[p]>=BMODE_MDBC2)velrhonew[p]=velrho2[p]; //-For mDBC2. //<vs_m2dbc>
      else{ //-For DBC and mDBC (SLIP_Vel0).
        float rrho=float(double(velrho2[p].w)+dt2*ar[p]);
        rrho=(rrho<rhopzero? rhopzero: rrho); //-To prevent absorption of fluid particles by boundaries. | Evita q las boundary absorvan a las fluidas.
        velrhonew[p]=make_float4(0,0,0,rrho);
      }
    }
    else{ //-Particles: Floating & Fluid.
      const typecode rcode=code[p];
      //-Updates density.
      const float4 rvelrho2=velrho2[p];
      const float rhonew=float(double(velrho2[p].w)+dt2*ar[p]);
      float4 rvel1=velrho1[p];
      if(!floating || CODE_IsFluid(rcode)){ //-Particles: Fluid.
        //-Calculate displacement. | Calcula desplazamiento.
        const float3 race=ace[p];
        const double acegrx=double(race.x)+gravity.x;
        const double acegry=double(race.y)+gravity.y;
        const double acegrz=double(race.z)+gravity.z;
        double dx=double(rvel1.x)*dt + acegrx*dt205;
        double dy=double(rvel1.y)*dt + acegry*dt205;
        double dz=double(rvel1.z)*dt + acegrz*dt205;
        if(shift){
          const float4 rshiftpos=shiftposfs[p];
          dx+=double(rshiftpos.x);
          dy+=double(rshiftpos.y);
          dz+=double(rshiftpos.z);
        }
        bool outrho=(rhonew<rhopoutmin || rhonew>rhopoutmax);
        //-Calculate velocity & density. | Calcula velocidad y densidad.
        float4 rvelrhonew=make_float4(
          float(double(rvelrho2.x) + acegrx*dt2),
          float(double(rvelrho2.y) + acegry*dt2),
          float(double(rvelrho2.z) + acegrz*dt2),
          rhonew);
        if(mdbc2==MDBC2_NoPen){//<vs_m2dbcNP_start>
          if(nopenshift[p].w>5.f){ //-check if correction should be applied or not
            if(nopenshift[p].x!=0){
                rvelrhonew.x=rvel1.x+nopenshift[p].x;//-Correcting velocity
                dx=double(rvelrhonew.x)*dt;//-Adding displacement
            }
            if(nopenshift[p].y!=0){
                rvelrhonew.y=rvel1.y+nopenshift[p].y;//-Correcting velocity
                dy=double(rvelrhonew.y)*dt;//-Adding displacement
            }
            if(nopenshift[p].z!=0){
                rvelrhonew.z=rvel1.z+nopenshift[p].z;//-Correcting velocity
                dz=double(rvelrhonew.z)*dt;//-Adding displacement
            }
          }
        }//<vs_m2dbcNP_end>
        //-Restore data of inout particles.
        if(inout && CODE_IsFluidInout(rcode)){
          outrho=false;
          rvelrhonew=rvelrho2;
          const float3 vd=indirvel[CODE_GetIzoneFluidInout(rcode)];
          if(vd.x!=FLT_MAX){
            const float v=rvel1.x*vd.x + rvel1.y*vd.y + rvel1.z*vd.z;
            dx=double(v*vd.x) * dt;
            dy=double(v*vd.y) * dt;
            dz=double(v*vd.z) * dt;
          }
          else{
            dx=double(rvel1.x) * dt;
            dy=double(rvel1.y) * dt;
            dz=double(rvel1.z) * dt;
          }
        }
        //-Update particle data.
        movxy[p]=make_double2(dx,dy);
        movz[p]=dz;
        if(outrho){ //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
          if(CODE_IsNormal(rcode))code[p]=CODE_SetOutRho(rcode);
        }
        velrhonew[p]=rvelrhonew;
      }
      else{ //-Particles: Floating.
        if(mdbc2==MDBC2_None || boundmode[p]<BMODE_MDBC2) //<vs_m2dbc>
          rvel1.w=(rhonew<rhopzero? rhopzero: rhonew); //-To prevent absorption of fluid particles by boundaries. | Evita q las floating absorvan a las fluidas.
        velrhonew[p]=rvel1;
      }
    }
  }
}
//==============================================================================
/// Updates particles according to forces and dt using Verlet. 
/// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void ComputeStepVerlet(bool floating,bool shift,bool inout,TpMdbc2Mode mdbc2
  ,unsigned np,unsigned npb,const float4* velrho1,const float4* velrho2
  ,const byte* boundmode,const float* ar,const float3* ace
  ,const float4* shiftposfs,const float3* indirvel,const float4* nopenshift,double dt,double dt2
  ,float rhopzero,float rhopoutmin,float rhopoutmax,tfloat3 gravity
  ,typecode* code,double2* movxy,double* movz,float4* velrhonew,cudaStream_t stm)
{
  double dt205=(0.5*dt*dt);
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
    if(inout){      const bool tinout=true;
      if(shift){    const bool shift=true;
        if(floating)KerComputeStepVerlet<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
        else        KerComputeStepVerlet<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
      }else{        const bool shift=false;
        if(floating)KerComputeStepVerlet<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
        else        KerComputeStepVerlet<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
      }
    }
    else{           const bool tinout=false;
      if(shift){    const bool shift=true;
        if(floating)KerComputeStepVerlet<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
        else        KerComputeStepVerlet<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
      }else{        const bool shift=false;
        if(floating)KerComputeStepVerlet<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
        else        KerComputeStepVerlet<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,rhopzero,rhopoutmin,rhopoutmax,velrho1,velrho2,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,dt,dt205,dt2,Float3(gravity),movxy,movz,code,velrhonew);
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Computes new values for Pos, Check, Vel and Ros (used with Symplectic-Predictor).
/// Calcula los nuevos valores de Pos, Vel y Rhop (usando para Symplectic-Predictor).
//------------------------------------------------------------------------------
template<bool floating,bool shift,bool inout> __global__ void KerComputeStepSymplecticPre
  (unsigned n,unsigned npb,TpMdbc2Mode mdbc2,const float4* velrhopre,const byte* boundmode
  ,const float* ar,const float3* ace,const float4* shiftposfs,const float3* indirvel
  ,double dtm,float rhopzero,float rhopoutmin,float rhopoutmax,float3 gravity
  ,typecode* code,double2* movxy,double* movz,float4* velrho
  ,float* psiclean,const float* psicleanpre,const float* psicleanrhs,bool divclean) //<vs_divclean>
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    if(p<npb){ //-Particles: Fixed & Moving.
      float4 vr=velrhopre[p];
      if(mdbc2==MDBC2_None || boundmode[p]<BMODE_MDBC2){ //-For DBC or mDBC1 //<vs_m2dbc>
        vr.w=float(double(vr.w)+dtm*ar[p]);
        vr.w=(vr.w<rhopzero? rhopzero: vr.w); //-To prevent absorption of fluid particles by boundaries. | Evita que las boundary absorvan a las fluidas.
      }
      velrho[p]=vr;
    }
    else{ //-Particles: Floating & Fluid.
      const typecode rcode=code[p];
      //-Updates density.
      const float4 rvelrhopre=velrhopre[p];
      float4 rvelrhonew=rvelrhopre;
      rvelrhonew.w=float(double(rvelrhopre.w)+dtm*ar[p]);

      #ifdef AVAILABLE_DIVCLEAN
      float rpsicleannew;
      float rpsicleanpre;
      if(divclean){
        rpsicleanpre=psicleanpre[p];        //<vs_divclean>
        rpsicleannew=rpsicleanpre;          //<vs_divclean> 
      }
      #endif
      
      if(!floating || CODE_IsFluid(rcode)){ //-Particles: Fluid.
        //-Calculate displacement. | Calcula desplazamiento.
        double dx=double(rvelrhopre.x)*dtm;
        double dy=double(rvelrhopre.y)*dtm;
        double dz=double(rvelrhopre.z)*dtm;
        if(shift){
          const float4 rshiftpos=shiftposfs[p];
          dx+=double(rshiftpos.x);
          dy+=double(rshiftpos.y);
          dz+=double(rshiftpos.z);
        }
        bool outrho=(rvelrhonew.w<rhopoutmin || rvelrhonew.w>rhopoutmax);
        //-Calculate velocity & density. | Calcula velocidad y densidad.
        const float3 race=ace[p];
        rvelrhonew.x=float(double(rvelrhopre.x) + (double(race.x)+gravity.x) * dtm);
        rvelrhonew.y=float(double(rvelrhopre.y) + (double(race.y)+gravity.y) * dtm);
        rvelrhonew.z=float(double(rvelrhopre.z) + (double(race.z)+gravity.z) * dtm);

        if(CODE_IsFluidBuffer(rcode)) rvelrhonew=rvelrhopre; //-Restore data of vres buffer particles. //<vs_vrres

        #ifdef AVAILABLE_DIVCLEAN
        if(divclean)rpsicleannew=rpsicleanpre+psicleanrhs[p]*dtm;
        #endif

        //-Restore data of inout particles.
        if(inout && CODE_IsFluidInout(rcode)){
          outrho=false;
          rvelrhonew=rvelrhopre;
          const float3 vd=indirvel[CODE_GetIzoneFluidInout(rcode)];
          if(vd.x!=FLT_MAX){
            const float v=rvelrhonew.x*vd.x + rvelrhonew.y*vd.y + rvelrhonew.z*vd.z;
            dx=double(v*vd.x) * dtm;
            dy=double(v*vd.y) * dtm;
            dz=double(v*vd.z) * dtm;
          }
        }
        //-Update particle data.
        movxy[p]=make_double2(dx,dy);
        movz[p]=dz;
        if(outrho){ //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
          if(CODE_IsNormal(rcode))code[p]=CODE_SetOutRho(rcode);
        }
      }
      else{ //-Particles: Floating.
        if(mdbc2==MDBC2_None || boundmode[p]<BMODE_MDBC2) //<vs_m2dbc>
          rvelrhonew.w=(rvelrhonew.w<rhopzero? rhopzero: rvelrhonew.w); //-To prevent absorption of fluid particles by boundaries. | Evita q las floating absorvan a las fluidas.
      }
      //-Stores new velocity and density.
      velrho[p]=rvelrhonew;
      #ifdef AVAILABLE_DIVCLEAN
      //-Update divergence cleaning scalar.
      if(divclean)psiclean[p]=rpsicleannew;
      #endif
    }
  }
}
//==============================================================================
/// Updates particles using Symplectic-Predictor.
/// Actualizacion de particulas usando Symplectic-Predictor.
//==============================================================================   
void ComputeStepSymplecticPre(bool floating,bool shift,bool inout,TpMdbc2Mode mdbc2
  ,unsigned np,unsigned npb,const float4* velrhopre,const byte* boundmode
  ,const float* ar,const float3* ace,const float4* shiftposfs,const float3* indirvel
  ,double dtm,float rhopzero,float rhopoutmin,float rhopoutmax,tfloat3 gravity
  ,typecode* code,double2* movxy,double* movz,float4* velrho
  ,float* psiclean,const float* psicleanpre,const float* psicleanrhs,bool divclean,cudaStream_t stm)  //<vs_divclean>
{
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
    if(inout){      const bool tinout=true;
      if(shift){    const bool shift=true;
        if(floating)KerComputeStepSymplecticPre<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticPre<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }else{        const bool shift=false;
        if(floating)KerComputeStepSymplecticPre<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticPre<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }
    }
    else{           const bool tinout=false;
      if(shift){    const bool shift=true;
        if(floating)KerComputeStepSymplecticPre<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticPre<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }else{        const bool shift=false;
        if(floating)KerComputeStepSymplecticPre<true ,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticPre<false,shift,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,dtm,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Computes new values for Pos, Check, Vel and Ros (using Verlet).
/// The value of Vel always set to be reset.
///
/// Calcula los nuevos valores de Pos, Vel y Rhop (usandopara Symplectic-Corrector).
/// Pone vel de contorno a cero.
//------------------------------------------------------------------------------
template<bool floating,bool shift,bool shiftadv,bool inout> __global__ void KerComputeStepSymplecticCor
  (unsigned n,unsigned npb,TpMdbc2Mode mdbc2,const float4* velrhopre,const byte* boundmode
  ,const float* ar,const float3* ace,const float4* shiftposfs,const float3* indirvel,const float4* nopenshift // SHABA
  ,const float4* shiftvel //<vs_advshift>
  ,double dtm,double dt,float rhopzero,float rhopoutmin,float rhopoutmax,float3 gravity
  ,typecode* code,double2* movxy,double* movz,float4* velrho
  ,float* psiclean,const float* psicleanpre,const float* psicleanrhs,bool divclean)  //<vs_divclean>
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    if(p<npb){ //-Particles: Fixed & Moving.
      if(mdbc2==MDBC2_None || boundmode[p]<BMODE_MDBC2){
        const double epsilon_rdot=(-double(ar[p])/double(velrho[p].w))*dt;
        float rrho=float(double(velrhopre[p].w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
        rrho=(rrho<rhopzero? rhopzero: rrho); //-To prevent absorption of fluid particles by boundaries. | Evita q las boundary absorvan a las fluidas.
        velrho[p]=make_float4(0,0,0,rrho);
      }
      else velrho[p]=velrhopre[p]; //-For mDBC2. //<vs_m2dbc>
    }
    else{ //-Particles: Floating & Fluid.
      const typecode rcode=code[p];
      //-Updates density.
      const double epsilon_rdot=(-double(ar[p])/double(velrho[p].w))*dt;
      const float4 rvelrhopre=velrhopre[p];
      float4 rvelrhonew=rvelrhopre;
      rvelrhonew.w=float(double(rvelrhopre.w) * (2.-epsilon_rdot)/(2.+epsilon_rdot));
      #ifdef AVAILABLE_DIVCLEAN
      float rpsicleannew;
      float rpsicleanpre;
      if(divclean){
        rpsicleanpre=psicleanpre[p];        //<vs_divclean>
        rpsicleannew=rpsicleanpre;          //<vs_divclean> 
      }
      #endif
      if(!floating || CODE_IsFluid(rcode)){//-Particles: Fluid.
        //-Calculate velocity. | Calcula velocidad.
        const float3 race=ace[p];
        rvelrhonew.x=float(double(rvelrhopre.x) + (double(race.x)+gravity.x) * dt);
        rvelrhonew.y=float(double(rvelrhopre.y) + (double(race.y)+gravity.y) * dt);
        rvelrhonew.z=float(double(rvelrhopre.z) + (double(race.z)+gravity.z) * dt);

        if(CODE_IsFluidBuffer(rcode)) rvelrhonew=velrho[p];   //-Restore data of vres buffer particles. //<vs_vrres
        
        #ifdef AVAILABLE_DIVCLEAN
        if(divclean)rpsicleannew=rpsicleanpre+psicleanrhs[p]*dt;
        #endif
        
        //-Calculate displacement. | Calcula desplazamiento.
        double dx=(double(rvelrhopre.x)+double(rvelrhonew.x)) * dtm;
        double dy=(double(rvelrhopre.y)+double(rvelrhonew.y)) * dtm;
        double dz=(double(rvelrhopre.z)+double(rvelrhonew.z)) * dtm;
        //-Adding no-penetration correction velocity SHABA
        if(mdbc2==MDBC2_NoPen){//<vs_m2dbcNP_start>
          if(nopenshift[p].w>5.f){ //-check if correction should be applied or not
            if(nopenshift[p].x!=0){
                rvelrhonew.x=rvelrhopre.x+nopenshift[p].x;//-Adding displacement
                dx=double(rvelrhonew.x)*dt;//-Adding displacement
            }
            if(nopenshift[p].y!=0){
                rvelrhonew.y=rvelrhopre.y+nopenshift[p].y;//-Adding displacement
                dy=double(rvelrhonew.y)*dt;//-Adding displacement
            }
            if(nopenshift[p].z!=0){
                rvelrhonew.z=rvelrhopre.z+nopenshift[p].z;//-Adding displacement
                dz=double(rvelrhonew.z)*dt;//-Adding displacement
            }
          }
        }//<vs_m2dbcNP_end>
        if(shift){
          const float4 rshiftpos=shiftposfs[p];
          dx+=double(rshiftpos.x);
          dy+=double(rshiftpos.y);
          dz+=double(rshiftpos.z);
        }
        if(shiftadv){ //<vs_advshift_ini>
          const float4 rshitvel=shiftvel[p];
          dx+=double(rshitvel.x)*dt;
          dy+=double(rshitvel.y)*dt;
          dz+=double(rshitvel.z)*dt;
        } //<vs_advshift_end>
        bool outrho=(rvelrhonew.w<rhopoutmin || rvelrhonew.w>rhopoutmax);
        //-Restore data of inout particles.
        if(inout && CODE_IsFluidInout(rcode)){
          outrho=false;
          rvelrhonew=rvelrhopre;
          const float3 vd=indirvel[CODE_GetIzoneFluidInout(rcode)];
          if(vd.x!=FLT_MAX){
            const float v=rvelrhonew.x*vd.x + rvelrhonew.y*vd.y + rvelrhonew.z*vd.z;
            dx=double(v*vd.x) * dt;
            dy=double(v*vd.y) * dt;
            dz=double(v*vd.z) * dt;
          }
          else{
            dx=double(rvelrhonew.x) * dt; 
            dy=double(rvelrhonew.y) * dt; 
            dz=double(rvelrhonew.z) * dt;
          }
        }
        //-Update particle data.
        movxy[p]=make_double2(dx,dy);
        movz[p]=dz;
        if(outrho){ //-Only brands as excluded normal particles (not periodic). | Solo marca como excluidas las normales (no periodicas).
          const typecode rcode=code[p];
          if(CODE_IsNormal(rcode))code[p]=CODE_SetOutRho(rcode);
        }
      }
      else{ //-Particles: Floating.
        if(mdbc2==MDBC2_None || boundmode[p]<BMODE_MDBC2) //<vs_m2dbc>
          rvelrhonew.w=(rvelrhonew.w<rhopzero? rhopzero: rvelrhonew.w); //-To prevent absorption of fluid particles by boundaries. | Evita q las floating absorvan a las fluidas.
      }
      //-Stores new velocity and density.
      velrho[p]=rvelrhonew;
      #ifdef AVAILABLE_DIVCLEAN
      //-Update divergence cleaning scalar.
      if(divclean)psiclean[p]=rpsicleannew;
      #endif
    }
  }
}
//==============================================================================
/// Updates particles using Symplectic-Corrector.
/// Actualizacion de particulas usando Symplectic-Corrector.
//==============================================================================   
void ComputeStepSymplecticCor(bool floating,bool shift,bool shiftadv,bool inout,TpMdbc2Mode mdbc2
  ,unsigned np,unsigned npb,const float4* velrhopre,const byte* boundmode
  ,const float* ar,const float3* ace,const float4* shiftposfs
  ,const float3* indirvel,const float4* nopenshift,const float4* shiftvel,double dtm,double dt
  ,float rhopzero,float rhopoutmin,float rhopoutmax,tfloat3 gravity
  ,typecode* code,double2* movxy,double* movz,float4* velrho
  ,float* psiclean,const float* psicleanpre,const float* psicleanrhs,bool divclean,cudaStream_t stm)  //<vs_divclean>
{
  if(np){
    dim3 sgrid=GetSimpleGridSize(np,SPHBSIZE);
    if(inout){      const bool tinout=true;
      if(shift){    const bool shift=true;  const bool shiftadv=false;
        if(floating)KerComputeStepSymplecticCor<true ,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticCor<false,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }else if(shiftadv){ const bool shift=false; const bool shiftadv=true; //<vs_advshift>
        if(floating)KerComputeStepSymplecticCor<true ,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean); //<vs_advshift>
        else        KerComputeStepSymplecticCor<false,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean); //<vs_advshift>
      }else{        const bool shift=false; const bool shiftadv=false;
        if(floating)KerComputeStepSymplecticCor<true ,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticCor<false,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }
    }
    else{           const bool tinout=false;
      if(shift){    const bool shift=true;  const bool shiftadv=false;
        if(floating)KerComputeStepSymplecticCor<true ,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticCor<false,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }else if(shiftadv){        const bool shift=false; const bool shiftadv=true; //<vs_advshift>
        if(floating)KerComputeStepSymplecticCor<true ,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean); //<vs_advshift>
        else        KerComputeStepSymplecticCor<false,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean); //<vs_advshift>
      }else{        const bool shift=false; const bool shiftadv=false;
        if(floating)KerComputeStepSymplecticCor<true ,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
        else        KerComputeStepSymplecticCor<false,shift,shiftadv,tinout> <<<sgrid,SPHBSIZE,0,stm>>> (np,npb,mdbc2,velrhopre,boundmode,ar,ace,shiftposfs,indirvel,nopenshift,shiftvel,dtm,dt,rhopzero,rhopoutmin,rhopoutmax,Float3(gravity),code,movxy,movz,velrho,psiclean,psicleanpre,psicleanrhs,divclean);
      }
    }
  }
}

//<vs_flexstruc_ini>
//==============================================================================
/// Copy motion velocity of flexible structure particles.
/// Copie la velocidad de movimiento de partículas de estructura flexible.
//==============================================================================
__global__ void KerCopyMotionVelFlexStruc(unsigned n,const typecode* code,const unsigned* flexstrucridp
    ,const float3* motionvel,float4* velrhop)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=flexstrucridp[p]; //-Number of particle.
    if(CODE_IsFlexStrucFlex(code[p1])){
      const float3 mvel=motionvel[p1];
      velrhop[p1]=make_float4(mvel.x,mvel.y,mvel.z,velrhop[p1].w);
    }
  }
}

//==============================================================================
/// Copy motion velocity of flexible structure particles.
/// Copie la velocidad de movimiento de partículas de estructura flexible.
//==============================================================================
void CopyMotionVelFlexStruc(unsigned npfs,const typecode* code,const unsigned* flexstrucridp
    ,const float3* motionvel,float4* velrhop)
{
  if(npfs){
    dim3 sgrid=GetSimpleGridSize(npfs,SPHBSIZE);
    KerCopyMotionVelFlexStruc <<<sgrid,SPHBSIZE>>> (npfs,code,flexstrucridp,motionvel,velrhop);
  }
}

//==============================================================================
/// Updates position and velocity using semi-implicit Euler scheme (used with Verlet scheme).
/// Actualiza la posición y la velocidad usando el esquema de Euler semiimplícito (usado con el esquema de Verlet).
//==============================================================================
__global__ void KerComputeStepFlexStrucSemiImplicitEuler(unsigned n,const float4* velrhop,const typecode* code,const unsigned* flexstrucridp
    ,const float3* ace,double dt,float3 gravity
    ,double2* movxy,double* movz,float4* velrhopnew)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=flexstrucridp[p]; //-Number of particle.
    if(CODE_IsFlexStrucFlex(code[p1])){
      const float4 rvelrhop=velrhop[p1];
      const float3 race=ace[p1];
      const double acegrx=double(race.x)+gravity.x;
      const double acegry=double(race.y)+gravity.y;
      const double acegrz=double(race.z)+gravity.z;
      const float4 rvelrhopnew=make_float4(
          float(double(rvelrhop.x) + acegrx*dt),
          float(double(rvelrhop.y) + acegry*dt),
          float(double(rvelrhop.z) + acegrz*dt),
          velrhopnew[p1].w);
      //-Calculate displacement. | Calcula desplazamiento.
      const double dx=double(rvelrhopnew.x)*dt;
      const double dy=double(rvelrhopnew.y)*dt;
      const double dz=double(rvelrhopnew.z)*dt;
      //-Update particle data.
      movxy[p1]=make_double2(dx,dy);
      movz[p1]=dz;
      velrhopnew[p1]=rvelrhopnew;
    }
  }
}

//==============================================================================
/// Updates position and velocity using semi-implicit Euler scheme (used with Verlet scheme).
/// Actualiza la posición y la velocidad usando el esquema de Euler semiimplícito (usado con el esquema de Verlet).
//==============================================================================
void ComputeStepFlexStrucSemiImplicitEuler(unsigned npfs,const float4* velrhop,const typecode* code,const unsigned* flexstrucridp
    ,const float3* ace,double dt,tfloat3 gravity
    ,double2* movxy,double* movz,float4* velrhopnew,cudaStream_t stm)
{
  if(npfs){
    dim3 sgrid=GetSimpleGridSize(npfs,SPHBSIZE);
    KerComputeStepFlexStrucSemiImplicitEuler <<<sgrid,SPHBSIZE,0,stm>>> (npfs,velrhop,code,flexstrucridp
        ,ace,dt,Float3(gravity)
        ,movxy,movz,velrhopnew);
  }
}

//==============================================================================
/// Updates position and velocity using symplectic predictor scheme.
/// Actualiza la posición y la velocidad utilizando un esquema predictor simpléctico.
//==============================================================================
__global__ void KerComputeStepFlexStrucSymplecticPre(unsigned n,const float4* velrhoppre,const typecode* code,const unsigned* flexstrucridp
    ,const float3* ace,double dtm,float3 gravity
    ,double2* movxy,double* movz,float4* velrhop)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=flexstrucridp[p]; //-Number of particle.
    if(CODE_IsFlexStrucFlex(code[p1])){
      const float4 rvelrhoppre=velrhoppre[p1];
      const float3 race=ace[p1];
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=double(rvelrhoppre.x)*dtm;
      double dy=double(rvelrhoppre.y)*dtm;
      double dz=double(rvelrhoppre.z)*dtm;
      //-Calculate velocity & density. | Calcula velocidad y densidad.
      const float4 rvelrhopnew=make_float4(
          float(double(rvelrhoppre.x) + (double(race.x)+gravity.x) * dtm),
          float(double(rvelrhoppre.y) + (double(race.y)+gravity.y) * dtm),
          float(double(rvelrhoppre.z) + (double(race.z)+gravity.z) * dtm),
          velrhop[p1].w);
      //-Update particle data.
      movxy[p1]=make_double2(dx,dy);
      movz[p1]=dz;
      velrhop[p1]=rvelrhopnew;
    }
  }
}

//==============================================================================
/// Updates position and velocity using symplectic predictor scheme.
/// Actualiza la posición y la velocidad utilizando un esquema predictor simpléctico.
//==============================================================================
void ComputeStepFlexStrucSymplecticPre(unsigned npfs,const float4* velrhoppre,const typecode* code,const unsigned* flexstrucridp
    ,const float3* ace,double dtm,tfloat3 gravity
    ,double2* movxy,double* movz,float4* velrhop,cudaStream_t stm)
{
  if(npfs){
    dim3 sgrid=GetSimpleGridSize(npfs,SPHBSIZE);
    KerComputeStepFlexStrucSymplecticPre <<<sgrid,SPHBSIZE,0,stm>>> (npfs,velrhoppre,code,flexstrucridp
        ,ace,dtm,Float3(gravity)
        ,movxy,movz,velrhop);
  }
}

//==============================================================================
/// Updates position and velocity using symplectic corrector scheme.
/// Actualiza la posición y la velocidad utilizando un esquema corrector simpléctico.
//==============================================================================
__global__ void KerComputeStepFlexStrucSymplecticCor(unsigned n,const float4* velrhoppre,const typecode* code,const unsigned* flexstrucridp
    ,const float3* ace,double dtm,double dt,float3 gravity
    ,double2* movxy,double* movz,float4* velrhop)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(p<n){
    const unsigned p1=flexstrucridp[p]; //-Number of particle.
    if(CODE_IsFlexStrucFlex(code[p1])){
      const float4 rvelrhoppre=velrhoppre[p1];
      const float3 race=ace[p1];
      //-Calculate velocity. | Calcula velocidad.
      const float4 rvelrhopnew=make_float4(
          float(double(rvelrhoppre.x) + (double(race.x)+gravity.x) * dt),
          float(double(rvelrhoppre.y) + (double(race.y)+gravity.y) * dt),
          float(double(rvelrhoppre.z) + (double(race.z)+gravity.z) * dt),
          velrhop[p1].w);
      //-Calculate displacement. | Calcula desplazamiento.
      double dx=(double(rvelrhoppre.x)+double(rvelrhopnew.x)) * dtm;
      double dy=(double(rvelrhoppre.y)+double(rvelrhopnew.y)) * dtm;
      double dz=(double(rvelrhoppre.z)+double(rvelrhopnew.z)) * dtm;
      //-Update particle data.
      movxy[p1]=make_double2(dx,dy);
      movz[p1]=dz;
      velrhop[p1]=rvelrhopnew;
    }
  }
}

//==============================================================================
/// Updates position and velocity using symplectic corrector scheme.
/// Actualiza la posición y la velocidad utilizando un esquema corrector simpléctico.
//==============================================================================
void ComputeStepFlexStrucSymplecticCor(unsigned npfs,const float4* velrhoppre,const typecode* code,const unsigned* flexstrucridp
    ,const float3* ace,double dtm,double dt,tfloat3 gravity
    ,double2* movxy,double* movz,float4* velrhop,cudaStream_t stm)
{
  if(npfs){
    dim3 sgrid=GetSimpleGridSize(npfs,SPHBSIZE);
    KerComputeStepFlexStrucSymplecticCor <<<sgrid,SPHBSIZE,0,stm>>> (npfs,velrhoppre,code,flexstrucridp
        ,ace,dtm,dt,Float3(gravity)
        ,movxy,movz,velrhop);
  }
}
//<vs_flexstruc_end>

}


