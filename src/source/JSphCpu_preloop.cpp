//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphCpu.cpp \brief Implements the class \ref JSphCpu.

#include "JSphCpu.h"
#include "JCellSearch_inline.h"
#include "FunSphKernel.h"
#include "FunctionsMath.h"

#include <climits>
#include <cmath>

using namespace std;

//==============================================================================
/// Creates list with free-surface particle (normal and periodic).
//==============================================================================
unsigned JSphCpu::CountFreeSurfaceParticles(unsigned npf,unsigned pini
  ,const unsigned* fstype,unsigned* listp)const
{
  unsigned count=0;
  const unsigned pfin=pini+npf;
  for(unsigned p=pini;p<pfin;p++){
    const unsigned fstypep=fstype[p];
    if(fstypep){//-It includes normal and periodic particles.
      listp[count]=p; count++;
    }
  }
  return(count);
}

//==============================================================================
/// Perform interaction between particles: Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Realiza interaccion entre particulas: Fluid/Float-Fluid/Float or Fluid/Float-Bound
//==============================================================================
template<TpKernel tker,bool sim2d> void JSphCpu::InteractionComputeFSNormals
  (unsigned np,unsigned pinit
  ,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,const unsigned* listp,unsigned* fstype,tfloat3* fsnormal)const
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const unsigned p1=listp[p];

    float fs_treshold=0;                                //-Divergence of the position.
    tfloat3 gradc=TFloat3(0);                           //-Gradient of the concentration
    unsigned neigh=0;                                   //-Number of neighbours.
      
    tmatrix3f lcorr;        lcorr=TMatrix3f(0);         //-Correction matrix.
    tmatrix3f lcorr_inv;    lcorr_inv=TMatrix3f(0);     //-Inverse of the correction matrix.

    float Nzero=0;

    if(sim2d)Nzero=float((3.141592)*KernelSize2/(Dp*Dp));
    else     Nzero=float((4.f/3.f)*(3.141592)*KernelSize2*KernelSize2/(Dp*Dp*Dp));
    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];

    bool boundp2=false;

    for(int b2=0;b2<2;b2++){
      if(b2==1)boundp2=true;

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

        //-Interaction of Fluid with type Fluid or Bound.
        //-----------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(posp1.x-pos[p2].x);
          const float dry=float(posp1.y-pos[p2].y);
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
          //-Computes kernel.
            const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
            const float rhopp2= float(velrho[p2].w);

            //===== Get mass of particle p2 ===== 
            float massp2=(boundp2? MassBound: MassFluid); //-Contiene masa de particula segun sea bound o fluid.

            // bool ftp2;
            // float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
            // const typecode cod=code[p2];
            // ftp2=CODE_IsFloating(cod);
            // ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);

            const float vol2=(float(massp2/rhopp2));
            neigh++;

            const float dot3=drx*frx+dry*fry+drz*frz;
            gradc.x+=vol2*frx;
            gradc.y+=vol2*fry;
            gradc.z+=vol2*frz;

            fs_treshold-=vol2*dot3;
            lcorr.a11+=-drx*frx*vol2; lcorr.a12+=-drx*fry*vol2; lcorr.a13+=-drx*frz*vol2;
            lcorr.a21+=-dry*frx*vol2; lcorr.a22+=-dry*fry*vol2; lcorr.a23+=-dry*frz*vol2;
            lcorr.a31+=-drz*frx*vol2; lcorr.a32+=-drz*fry*vol2; lcorr.a33+=-drz*frz*vol2;

            // const float wab=cufsph::GetKernel_Wab<KERNEL_Wendland>(rr2);
            // pou+=wab*vol2;  
          }
        }
      }
    }
    //-Find particles that are probably on the free-surface.
    unsigned fstypep1=0;
    if(sim2d){
      if(fs_treshold<1.7)fstypep1=2;
      if(fs_treshold<1.1 && Nzero/float(neigh)<0.4f)fstypep1=3;
    }
    else{
      if(fs_treshold<2.75)fstypep1=2;
      if(fs_treshold<1.8 && Nzero/float(neigh)<0.4f)fstypep1=3;
    }
    fstype[p1]=fstypep1;
    
    //-Calculation of the inverse of the correction matrix (Don't think there is a better way, create function for Determinant2x2 for clarity?).
    if(sim2d){
      tmatrix2f lcorr2d;
      tmatrix2f lcorr2d_inv;
      lcorr2d.a11=lcorr.a11; lcorr2d.a12=lcorr.a13;
      lcorr2d.a21=lcorr.a31; lcorr2d.a22=lcorr.a33;
      float lcorr_det=(lcorr2d.a11*lcorr2d.a22-lcorr2d.a12*lcorr2d.a21);
      lcorr2d_inv.a11=lcorr2d.a22/lcorr_det; lcorr2d_inv.a12=-lcorr2d.a12/lcorr_det; lcorr2d_inv.a22=lcorr2d.a11/lcorr_det; lcorr2d_inv.a21=-lcorr2d.a21/lcorr_det;
      lcorr_inv.a11=lcorr2d_inv.a11;  lcorr_inv.a13=lcorr2d_inv.a12;
      lcorr_inv.a31=lcorr2d_inv.a21;  lcorr_inv.a33=lcorr2d_inv.a22;
    }
    else{
      const float determ=fmath::Determinant3x3(lcorr);
      lcorr_inv=fmath::InverseMatrix3x3(lcorr,determ);
    }

    //-Correction of the gradient of concentration.
    tfloat3 gradc1=TFloat3(0);    
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

//==============================================================================
/// Computes free-surface particles and their normals.
//==============================================================================
void JSphCpu::CallComputeFSNormals(const StDivDataCpu& divdata
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,const tfloat4* velrho,unsigned* fstype,tfloat3* fsnormal,unsigned* listp)const
{
  const unsigned npf=Np-Npb;
  if(npf){
    //-Creates list with free-surface particle (normal and periodic).
    const unsigned count=CountFreeSurfaceParticles(npf,Npb,fstype,listp);
    //-Computes normals on selected free-surface particles.
    if(count){
      if(Simulate2D){ const bool sim2d=true;
        switch(TKernel){
          case KERNEL_Wendland:  InteractionComputeFSNormals<KERNEL_Wendland,sim2d>(count,Npb,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal);  break;
          case KERNEL_Cubic:     InteractionComputeFSNormals<KERNEL_Cubic   ,sim2d>(count,Npb,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal);  break;
          default: Run_Exceptioon("Kernel unknown.");
        }
      }
      else{ const bool sim2d=false;
        switch(TKernel){
          case KERNEL_Wendland:  InteractionComputeFSNormals<KERNEL_Wendland,sim2d>(count,Npb,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal);  break;
          case KERNEL_Cubic:     InteractionComputeFSNormals<KERNEL_Cubic   ,sim2d>(count,Npb,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal);  break;
          default: Run_Exceptioon("Kernel unknown.");
        }
      }
    }
  }
}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid.
/// Interaccion Fluid-Fluid/Bound & Bound-Fluid.
//==============================================================================
void JSphCpu::InteractionCallScanUmbrellaRegion(unsigned np,unsigned pinit
  ,StDivDataCpu divdata,const unsigned* dcell,const tdouble3* pos
  ,const typecode* code,const tfloat3* fsnormal,const unsigned* listp
  ,unsigned* fstype)const
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const unsigned p1=listp[p];
    bool fs_flag=false;
    bool boundp2=false;

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];

    for(int b2=0; b2<2;b2++){
      if(b2==1) boundp2=true;
      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

        //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
        //------------------------------------------------------------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(posp1.x-pos[p2].x);
          const float dry=float(posp1.y-pos[p2].y);
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
            //-Computes kernel.
            const tfloat3 posq=fsnormal[p1]*KernelH;
            if(rr2>2.f*KernelH*KernelH){
              const float drxq=-drx-posq.x;
              const float dryq=-dry-posq.y;
              const float drzq=-drz-posq.z;
              const float rrq=sqrt(drxq*drxq+dryq*dryq+drzq*drzq);
              if(rrq<KernelH)fs_flag=true;
            }
            else{
              if(Simulate2D){
                const float drxq=-drx-posq.x;
                const float drzq=-drz-posq.z;
                const tfloat3 normalq=TFloat3(drxq*fsnormal[p1].x,0,drzq*fsnormal[p1].z);
                const tfloat3 tangq=TFloat3(-drxq*fsnormal[p1].z,0,drzq*fsnormal[p1].x);
                const float normalqnorm=sqrt(normalq.x*normalq.x+normalq.z*normalq.z);
                const float tangqnorm=sqrt(tangq.x*tangq.x+tangq.z*tangq.z);
                if(normalqnorm+tangqnorm<KernelH)fs_flag=true;
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
void JSphCpu::CallScanUmbrellaRegion(const StDivDataCpu& divdata
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,const tfloat3* fsnormal,unsigned* listp,unsigned* fstype)const
{
  const unsigned npf=Np-Npb;
  if(npf){
    //-Obtain the list of particle that are probably on the free-surface (in ComputeUmbrellaRegion maybe is unnecessary).
    const unsigned count=CountFreeSurfaceParticles(npf,Npb,fstype,listp);
    //-Scan Umbrella region on selected free-surface particles.
    if(count)InteractionCallScanUmbrellaRegion(count,Npb,divdata,dcell,pos,code,fsnormal,listp,fstype);
  }
}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid for models before
/// InteractionForces
//==============================================================================
template<TpKernel tker,bool sim2d,bool shiftadv> void JSphCpu::PreLoopInteraction
  (unsigned n,unsigned pinit ,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat4* shiftvel,tfloat3* fsnormal,float* fsmindist)const
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<pfin;p1++){
    float fs_treshold=0.f;                    //-Divergence of the position.
    tfloat3 gradc=TFloat3(0);                 //-Gradient of the concentration
    unsigned neigh=0;                         //-Number of neighbours.
    tfloat4 shiftvelp1=shiftvel[p1];          //-Shifting Velocity vector
    tfloat3 fsnormalp1=TFloat3(0);            //-Free-surface vector
    bool    nearfs=false;                     //-Bool for detecting near free-surface particles.                <ShiftingAdvanced>
    float   pou=false;                        //-Partition of unity for normal correction.                      <ShiftingAdvanced>
      
    float mindist=KernelSize;                 //-Set Min Distance from free-surface to kernel radius.           <ShiftingAdvanced>
    float maxarccos=0.0;                      //-Variable for identify high-curvature free-surface particle     <ShiftingAdvanced>
    bool bound_inter=false;                   //-Variable for identify free-surface that interact with boundary <ShiftingAdvanced>
    unsigned fsp1=fstype[p1];         
    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];

    bool boundp2=false;

    for(int b2=0;b2<2;b2++){
      if(b2==1)boundp2=true;

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

        //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
        //------------------------------------------------------------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(posp1.x-pos[p2].x);
          const float dry=float(posp1.y-pos[p2].y);
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
            //-Computes kernel.
            const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

            const float rhopp2= float(velrho[p2].w);

            //===== Get mass of particle p2 ===== 
            float massp2=(boundp2? MassBound: MassFluid); //-Contiene masa de particula segun sea bound o fluid.
            bool ftp2;
            //float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
            const typecode cod=code[p2];
            ftp2=CODE_IsFloating(cod);
            // ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);

            if(shiftadv){
              const float massrho=massp2/rhopp2;

              //-Compute gradient of concentration and partition of unity.        
              shiftvelp1.x+=massrho*frx;    
              shiftvelp1.y+=massrho*fry;
              shiftvelp1.z+=massrho*frz;          
              const float wab=fsph::GetKernel_Wab<tker>(CSP,rr2);
              shiftvelp1.w+=wab*massrho;

              //-Check if the particle is too close to solid or floating object
              if(boundp2 || ftp2)bound_inter=true;

              //-Check if it close to the free-surface, calculate distance from free-surface and smoothing of free-surface normals.
              if(fstype[p2]==2 && !boundp2){
                nearfs=true;
                mindist=min(sqrt(rr2),mindist);
                pou+=wab*massrho;
                fsnormalp1.x+=fsnormal[p2].x*wab*massrho;
                fsnormalp1.y+=fsnormal[p2].y*wab*massrho;
                fsnormalp1.z+=fsnormal[p2].z*wab*massrho;
              }

              //-Check maximum curvature.
              if(fstype[p1]>1 && fstype[p2]>1){
                const float norm1=sqrt(fsnormal[p1].x*fsnormal[p1].x+fsnormal[p1].y*fsnormal[p1].y+fsnormal[p1].z*fsnormal[p1].z);
                const float norm2=sqrt(fsnormal[p2].x*fsnormal[p2].x+fsnormal[p2].y*fsnormal[p2].y+fsnormal[p2].z*fsnormal[p2].z);
                maxarccos=max(maxarccos,(acos((fsnormal[p1].x*fsnormal[p2].x+fsnormal[p2].y*fsnormal[p1].y+fsnormal[p2].z*fsnormal[p1].z))));
              }
            }
          }
        }
      }
    }

    if(shiftadv){
      shiftvelp1.w+=fsph::GetKernel_Wab<tker>(CSP,0)*MassFluid/velrho[p1].w;
      fsmindist[p1]=mindist;
      //-Assign correct code to near free-surface particle and correct their normals by Shepard's Correction.
      if(fsp1==0 && nearfs){
        if(pou>1e-6){
          fsnormalp1=TFloat3(fsnormalp1.x,fsnormalp1.y,fsnormalp1.z);
          float norm=sqrt(fsnormalp1.x*fsnormalp1.x+fsnormalp1.y*fsnormalp1.y+fsnormalp1.z*fsnormalp1.z);
          fsnormal[p1]=TFloat3(fsnormalp1.x/norm,fsnormalp1.y/norm,fsnormalp1.z/norm);
        }
        fstype[p1]=1;
        if(bound_inter)fstype[p1]=3;
      }
      //-Check if free-surface particle interact with bound or has high-curvature.
      if(fsp1==2 && (bound_inter||maxarccos>0.52f))shiftvelp1=TFloat4(0,0,0,shiftvelp1.w);
      //-Compute shifting when <shiftImproved> true
      shiftvel[p1]=shiftvelp1;
    }
  }
}
//==============================================================================
template<TpKernel tker,bool sim2d> void JSphCpu::PreLoopInteraction_ct1
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat4* shiftvel,tfloat3* fsnormal,float* fsmindist)const
{
  if(ShiftingAdv)PreLoopInteraction <tker,sim2d,true >
    (n,pinit,divdata,dcell,pos,code,velrho,fstype,shiftvel,fsnormal,fsmindist);  
  else           PreLoopInteraction <tker,sim2d,false>
    (n,pinit,divdata,dcell,pos,code,velrho,fstype,shiftvel,fsnormal,fsmindist);  
}
//==============================================================================
void JSphCpu::PreLoopInteraction_ct(const StDivDataCpu& divdata
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,const tfloat4* velrho,unsigned* fstype,tfloat4* shiftvel,tfloat3* fsnormal
  ,float* fsmindist)const
{
  const unsigned npf=Np-Npb;
  if(npf){
    if(Simulate2D){ const bool sim2d=true;
      switch(TKernel){
        case KERNEL_Cubic:       PreLoopInteraction_ct1 <KERNEL_Cubic     ,sim2d> (npf,Npb,divdata,dcell,pos,code,velrho,fstype,shiftvel,fsnormal,fsmindist);  break;
        case KERNEL_Wendland:    PreLoopInteraction_ct1 <KERNEL_Wendland  ,sim2d> (npf,Npb,divdata,dcell,pos,code,velrho,fstype,shiftvel,fsnormal,fsmindist);  break;
        default: Run_Exceptioon("Kernel unknown.");
      }
    }
    else{ const bool sim2d=false;
      switch(TKernel){
        case KERNEL_Cubic:       PreLoopInteraction_ct1 <KERNEL_Cubic     ,sim2d> (npf,Npb,divdata,dcell,pos,code,velrho,fstype,shiftvel,fsnormal,fsmindist);  break;
        case KERNEL_Wendland:    PreLoopInteraction_ct1 <KERNEL_Wendland  ,sim2d> (npf,Npb,divdata,dcell,pos,code,velrho,fstype,shiftvel,fsnormal,fsmindist);  break;
        default: Run_Exceptioon("Kernel unknown.");
      }
    }
  }
}

//==============================================================================
/// Compute shifting velocity for advanced shifting model.
//==============================================================================
void JSphCpu::ComputeShiftingVel(bool sim2d,float shiftcoef,bool ale,double dt
  ,const unsigned* fstype,const tfloat3* fsnormal,const float* fsmindist
  ,tfloat4* shiftvel)const
{
  const int np=int(Np);
  const int npb=int(Npb);
  const int npf=np-npb;
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int p=npb;p<np;p++){
    const unsigned fstypep1=fstype[p];
    const tfloat4  shiftp1=shiftvel[p];
    const float    fsmindistp1=fsmindist[p];
    const tfloat3  fsnormalp1=fsnormal[p];
    const float    theta=min(1.f,max(0.f,(fsmindistp1-(float)KernelSize)/(float)(0.5*KernelSize-KernelSize)));
    tfloat4        shift_final=TFloat4(0,0,0,0);
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
        shift_final=TFloat4(0,0,0,0);
      }
    }
    else if(fstypep1==3){
      shift_final=TFloat4(0,0,0,0);
    }

    if(sim2d)shift_final.y=0.0;

    const float rhovar=abs(shiftp1.x*shift_final.x+shiftp1.y*shift_final.y+shiftp1.z*shift_final.z)*min(KernelH,fsmindistp1);
    const float eps=1e-5f;
    const float umagn_1=float(shiftcoef*KernelH/dt);
    const float umagn_2=float(abs(eps/(2.f*dt*rhovar)));
    const float umagn=float(min(umagn_1,umagn_2)*min(KernelH,fsmindistp1)*dt);

    const float maxdist=float(0.1f*Dp);
    shift_final.x=(fabs(umagn*shift_final.x)<maxdist? umagn*shift_final.x: (umagn*shift_final.x>=0? maxdist: -maxdist));
    shift_final.y=(fabs(umagn*shift_final.y)<maxdist? umagn*shift_final.y: (umagn*shift_final.y>=0? maxdist: -maxdist));
    shift_final.z=(fabs(umagn*shift_final.z)<maxdist? umagn*shift_final.z: (umagn*shift_final.z>=0? maxdist: -maxdist));
    shiftvel[p].x=float((shift_final.x)/dt);
    shiftvel[p].y=float((shift_final.y)/dt);
    shiftvel[p].z=float((shift_final.z)/dt);
  }
}

