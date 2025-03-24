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

using namespace std;

//##############################################################################
//# JSphCpu for mDBC
//##############################################################################
//==============================================================================
/// Perform interaction between ghost nodes of boundaries and fluid.
//==============================================================================
template<TpKernel tker,bool sim2d> void JSphCpu::InteractionMdbcCorrectionT2
  (unsigned n,StDivDataCpu divdata,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat3* boundnor,tfloat4* velrho)
{
  const float determlimit=1e-3f;
  const int nn=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<nn;p1++)if(boundnor[p1]!=TFloat3(0)){
    float rhofinal=FLT_MAX;

    //-Calculates ghost node position.
    tdouble3 gposp1=pos[p1]+ToTDouble3(boundnor[p1]);
    gposp1=(PeriActive!=0? UpdatePeriodicPos(gposp1): gposp1); //-Corrected interface Position.
    //-Initializes variables for calculation.
    float rhop1=0;
    tfloat3 gradrhop1=TFloat3(0);
    tdouble3 velp1=TDouble3(0);       //-Only for velocity.
    tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
    tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(gposp1,false,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);
      //-Interaction of boundary with type Fluid/Float.
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float drx=float(gposp1.x-pos[p2].x);
        const float dry=float(gposp1.y-pos[p2].y);
        const float drz=float(gposp1.z-pos[p2].z);
        const float rr2=(drx*drx + dry*dry + drz*drz);
        if(rr2<=KernelSize2 && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
          //-Wendland kernel.
          float fac;
          const float wab=fsph::GetKernel_WabFac<tker>(CSP,rr2,fac);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass and volume of particle p2 =====
          const tfloat4 velrhop2=velrho[p2];
          const float massp2=MassFluid;
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
      const tfloat3 dpos=(boundnor[p1]*(-1.f)); //-Boundary particle position - ghost node position.
      if(sim2d){
        const double determ=fmath::Determinant3x3(a_corr2);
        if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
          const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
          const float rhoghost=float(invacorr2.a11*rhop1 + invacorr2.a12*gradrhop1.x + invacorr2.a13*gradrhop1.z);
          const float grx=    -float(invacorr2.a21*rhop1 + invacorr2.a22*gradrhop1.x + invacorr2.a23*gradrhop1.z);
          const float grz=    -float(invacorr2.a31*rhop1 + invacorr2.a32*gradrhop1.x + invacorr2.a33*gradrhop1.z);
          rhofinal=(rhoghost + grx*dpos.x + grz*dpos.z);
        }
        else if(a_corr2.a11>0){//-Determinant is small but a11 is nonzero (0th order).
          rhofinal=float(rhop1/a_corr2.a11);
        }
      }
      else{
        const double determ=fmath::Determinant4x4(a_corr3);
        if(fabs(determ)>=determlimit){
          const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE BOUNDARY PARTICLES.
          const float rhoghost=float(invacorr3.a11*rhop1 + invacorr3.a12*gradrhop1.x + invacorr3.a13*gradrhop1.y + invacorr3.a14*gradrhop1.z);
          const float grx=    -float(invacorr3.a21*rhop1 + invacorr3.a22*gradrhop1.x + invacorr3.a23*gradrhop1.y + invacorr3.a24*gradrhop1.z);
          const float gry=    -float(invacorr3.a31*rhop1 + invacorr3.a32*gradrhop1.x + invacorr3.a33*gradrhop1.y + invacorr3.a34*gradrhop1.z);
          const float grz=    -float(invacorr3.a41*rhop1 + invacorr3.a42*gradrhop1.x + invacorr3.a43*gradrhop1.y + invacorr3.a44*gradrhop1.z);
          rhofinal=(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
        }
        else if(a_corr3.a11>0){//-Determinant is small but a11 is nonzero (0th order).
          rhofinal=float(rhop1/a_corr3.a11);
        }
      }
      //-Store the results.
      rhofinal=(rhofinal!=FLT_MAX? rhofinal: RhopZero);
      //-SlipMode==SLIP_Vel0 (DBC vel=0)
      velrho[p1].w=rhofinal;
    }
  }
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
 template<TpKernel tker> void JSphCpu::Interaction_MdbcCorrectionT
  (const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat3* boundnor,tfloat4* velrho)
{
  //-Interaction GhostBoundaryNodes-Fluid.
  const unsigned n=(UseNormalsFt? Np: NpbOk);
  if(Simulate2D)InteractionMdbcCorrectionT2 <tker,true >(n,divdata,pos,code,idp,boundnor,velrho);
  else          InteractionMdbcCorrectionT2 <tker,false>(n,divdata,pos,code,idp,boundnor,velrho);
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphCpu::Interaction_MdbcCorrection(const StDivDataCpu &divdata
  ,const tdouble3* pos,const typecode* code,const unsigned* idp
  ,const tfloat3* boundnor,tfloat4* velrho)
{
  switch(TKernel){
    case KERNEL_Cubic:       Interaction_MdbcCorrectionT <KERNEL_Cubic     > (divdata,pos,code,idp,boundnor,velrho);  break;
    case KERNEL_Wendland:    Interaction_MdbcCorrectionT <KERNEL_Wendland  > (divdata,pos,code,idp,boundnor,velrho);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

//<vs_m2dbc_ini>
//------------------------------------------------------------------------------
/// Perform Pressure cloning for mDBC2.
//------------------------------------------------------------------------------
//  const float pressfinal=KerMdbc2PressClone(sim2d,rhoghost,bnormalp1,gravity,motace,dpos);
float JSphCpu::Mdbc2PressClone(bool sim2d,const float rhoghost,tfloat3 bnormalp1
  ,const tfloat3 gravity,const tfloat3 motacep1,const tfloat3 dpos)const
{
  float pressfinal=0.f;
  if(sim2d){
    const float pghost=float(Cs0*Cs0*(rhoghost-RhopZero));
    const float norm=sqrt(bnormalp1.x*bnormalp1.x + bnormalp1.z*bnormalp1.z);
    const float normx=bnormalp1.x/norm; 
    const float normz=bnormalp1.z/norm;
    const float normpos=dpos.x*normx + dpos.z*normz;
    const tfloat3 force=TFloat3(gravity.x-motacep1.x,0,gravity.z-motacep1.z);
    const float normforce=RhopZero*(force.x*normx + force.z*normz);
    pressfinal=pghost+normforce*normpos;
  }
  else{
    const float pghost=float(Cs0*Cs0*(rhoghost-RhopZero));
    const float norm=sqrt(bnormalp1.x*bnormalp1.x + bnormalp1.y*bnormalp1.y + bnormalp1.z*bnormalp1.z);
    const float normx=bnormalp1.x/norm;
    const float normy=bnormalp1.y/norm;
    const float normz=bnormalp1.z/norm;
    const float normpos=dpos.x*normx + dpos.y*normy + dpos.z*normz;
    const tfloat3 force=gravity-motacep1;
    const float normforce=RhopZero*(force.x*normx + force.y*normy + force.z*normz);
    pressfinal=pghost+normforce*normpos;
  }
  return(pressfinal);
}
//------------------------------------------------------------------------------
/// Calculates the infinity norm of a 3x3 matrix.
//------------------------------------------------------------------------------
float JSphCpu::Mdbc2InfNorm3x3(tmatrix3d mat)const{
  const float row1=float(fabs(mat.a11) + fabs(mat.a12) + fabs(mat.a13));
  const float row2=float(fabs(mat.a21) + fabs(mat.a22) + fabs(mat.a23));
  const float row3=float(fabs(mat.a31) + fabs(mat.a32) + fabs(mat.a33));
  const float infnorm=max(row1,max(row2,row3));
  return(infnorm);
}
//------------------------------------------------------------------------------
/// Calculates the infinity norm of a 4x4 matrix.
//------------------------------------------------------------------------------
float JSphCpu::Mdbc2InfNorm4x4(tmatrix4d mat)const{
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
tfloat3 JSphCpu::Mdbc2TangenVel(const tfloat3& boundnor,const tfloat3& velfinal)const{
  const float snormal=sqrt(boundnor.x*boundnor.x + boundnor.y*boundnor.y + boundnor.z*boundnor.z);
  const tfloat3 bnormal=boundnor/snormal;
	const float veldotnorm=velfinal.x*bnormal.x + velfinal.y*bnormal.y + velfinal.z*bnormal.z;
  const tfloat3 tangentvel=TFloat3(velfinal.x-veldotnorm*bnormal.x,
                                   velfinal.y-veldotnorm*bnormal.y,
                                   velfinal.z-veldotnorm*bnormal.z);
  return(tangentvel);
}

//==============================================================================
/// Perform interaction between ghost nodes of boundaries and fluid.
//==============================================================================
template<TpKernel tker,bool sim2d> void JSphCpu::InteractionMdbc2CorrectionT2
  (unsigned n,const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat3* boundnor,const tfloat3* motionvel
  ,const tfloat3* motionace,tfloat4* velrho,byte* boundmode,tfloat3* tangenvel)
{
  const int nn=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<nn;p1++){
    const tfloat3 bnormalp1=boundnor[p1];
    if(bnormalp1!=TFloat3(0)){
      float rhofinal=FLT_MAX;
      tfloat3 velrhofinal=TFloat3(0);
      float sumwab=0;
      float submerged=0;                    //-Not for Vel0.
      const tfloat3 motacep1=motionace[p1]; //-Not for Vel0.

      //-Calculates ghost node position.
      tdouble3 gposp1=pos[p1]+ToTDouble3(bnormalp1);
      //-Corrected interface Position.
      gposp1=(PeriActive!=0? UpdatePeriodicPos(gposp1): gposp1);

      //-Initializes variables for calculation.
      float rhop1=0;
      tfloat3 gradrhop1=TFloat3(0);
      tdouble3 velp1=TDouble3(0);       //-Only for velocity.
      tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
      tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(gposp1,false,divdata);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);
        //-Interaction of boundary with type Fluid/Float.
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(gposp1.x-pos[p2].x);
          const float dry=float(gposp1.y-pos[p2].y);
          const float drz=float(gposp1.z-pos[p2].z);
          const float rr2=(drx*drx + dry*dry + drz*drz);
          if(rr2<=KernelSize2 && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
            //-Wendland kernel.
            float fac;
            const float wab=fsph::GetKernel_WabFac<tker>(CSP,rr2,fac);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

            //===== Get mass and volume of particle p2 =====
            const tfloat4 velrhop2=velrho[p2];
            const float massp2=MassFluid;
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
        const tfloat3 dpos=(boundnor[p1]*(-1.f)); //-Boundary particle position - ghost node position.
        if(sim2d){//-2D simulation.
          if(sumwab<0.1f){ //-If kernel sum is small use shepherd density and vel0.
            //-Trying to avoid too negative pressures.
            rhofinal=max(RhopZero,float(rhop1/a_corr2.a11));
          }
          else{//-Chech if matrix is invertible and well conditioned.
            const double determ=fmath::Determinant3x3(a_corr2);
            if(fabs(determ)>=0.001){//-Use 1e-3f (first_order).
              const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);
              //-Calculate the scaled condition number
              const float infnorma   =Mdbc2InfNorm3x3(a_corr2);
              const float infnormainv=Mdbc2InfNorm3x3(invacorr2);
              const float condinf=float(Dp*Dp)*infnorma*infnormainv;
              if(condinf<=50 && !CODE_IsFlexStrucFlex(code[p1]))//-If matrix is well conditioned use matrix inverse for density and shepherd for velocity.
                rhofinal=float(invacorr2.a11*rhop1 + invacorr2.a12*gradrhop1.x + invacorr2.a13*gradrhop1.z);
              else//-If ill conditioned use shepherd.
                rhofinal=float(rhop1/a_corr2.a11);
            }
            else//-If not invertible use shepherd for density.
              rhofinal=float(rhop1/a_corr2.a11);
          }
          //-Final density according to press.
          const float pressfinal=Mdbc2PressClone(sim2d,rhofinal,bnormalp1,Gravity,motacep1,dpos);
          rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
          //-velocity with Shepherd Sum.
          velrhofinal.x=float(velp1.x/a_corr2.a11);
          velrhofinal.y=0.f;
          velrhofinal.z=float(velp1.z/a_corr2.a11);
        }
        else{//-3D simulation.
          //-Density with pressure cloning.
          if(sumwab<0.1f){//-If kernel sum is small use shepherd density and vel0.
            //-Trying to avoid too negative pressures.
            rhofinal=max(RhopZero,float(rhop1/a_corr3.a11));
          }
          else{
            const double determ=fmath::Determinant4x4(a_corr3);
            if(fabs(determ)>=0.001){
              const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
              //-Calculate the scaled condition number.
              const float infnorma   =Mdbc2InfNorm4x4(a_corr3);
              const float infnormainv=Mdbc2InfNorm4x4(invacorr3);
              const float condinf=float(Dp*Dp)*infnorma*infnormainv;
              if(condinf<=50 && !CODE_IsFlexStrucFlex(code[p1]))//-If matrix is well conditioned use matrix inverse for density.
                rhofinal=float(invacorr3.a11*rhop1 + invacorr3.a12*gradrhop1.x + invacorr3.a13*gradrhop1.y + invacorr3.a14*gradrhop1.z);
              else//-Matrix is not well conditioned, use shepherd for density.
                rhofinal=float(rhop1/a_corr3.a11);
            }
            else//-Matrix is not invertible use shepherd for density.
              rhofinal=float(rhop1/a_corr3.a11);
          }
          //-Final density according to press.
          const float pressfinal=Mdbc2PressClone(sim2d,rhofinal,bnormalp1,Gravity,motacep1,dpos);
          rhofinal=RhopZero+float(pressfinal/(Cs0*Cs0));
          //-Velocity with Shepherd sum.
          velrhofinal.x=float(velp1.x/a_corr3.a11);
          velrhofinal.y=float(velp1.y/a_corr3.a11);
          velrhofinal.z=float(velp1.z/a_corr3.a11);
        }

        //-Store the results.
        if(SlipMode==SLIP_NoSlip){//-No-Slip: vel = 2*motion - ghost
          const tfloat3 v=motionvel[p1];
          const tfloat3 v2=TFloat3(v.x+v.x - velrhofinal.x,
                                   v.y+v.y - velrhofinal.y,
                                   v.z+v.z - velrhofinal.z); 
          #ifndef MDBC2_KEEPVEL
            velrho[p1]=TFloat4(v2.x,v2.y,v2.z,rhofinal);
          #else
            velrho[p1].w=rhofinal;
          #endif
          tangenvel[p1]=Mdbc2TangenVel(bnormalp1,v2);
        }
        else if(SlipMode==SLIP_FreeSlip){//-Free slip: vel = ghost vel
            // copy velocity from ghost node
            const tfloat3 v2 = TFloat3(velrhofinal.x, velrhofinal.y, velrhofinal.z);
            #ifndef MDBC2_KEEPVEL
                velrho[p1] = TFloat4(v2.x, v2.y, v2.z, rhofinal);
            #else
                velrho[p1].w = rhofinal;
            #endif
            tangenvel[p1] = Mdbc2TangenVel(bnormalp1, v2);
        }
      }
      else{//-If unsubmerged switch off boundary particle.
        boundmode[p1]=BMODE_MDBC2OFF;
        const tfloat3 v=motionvel[p1];
        #ifndef MDBC2_KEEPVEL
          velrho[p1]=TFloat4(v.x,v.y,v.z,RhopZero);
        #else
          velrho[p1].w=RhopZero;
        #endif
        tangenvel[p1]=Mdbc2TangenVel(bnormalp1,v);
      }
    }
  }
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
 template<TpKernel tker> void JSphCpu::Interaction_Mdbc2CorrectionT
  (const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat3* boundnor,const tfloat3* motionvel
  ,const tfloat3* motionace,tfloat4* velrho,byte* boundmode,tfloat3* tangenvel)
{
  //-Interaction GhostBoundaryNodes-Fluid.
  const unsigned n=(UseNormalsFt? Np: Npb);
  if(Simulate2D)InteractionMdbc2CorrectionT2 <tker,true >(n,divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundmode,tangenvel);
  else          InteractionMdbc2CorrectionT2 <tker,false>(n,divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundmode,tangenvel);
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphCpu::Interaction_Mdbc2Correction(const StDivDataCpu &divdata
  ,const tdouble3* pos,const typecode* code,const unsigned* idp
  ,const tfloat3* boundnor,const tfloat3* motionvel,const tfloat3* motionace
  ,tfloat4* velrho,byte* boundmode,tfloat3* tangenvel)
{
  switch(TKernel){
    case KERNEL_Cubic:       Interaction_Mdbc2CorrectionT <KERNEL_Cubic     > (divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundmode,tangenvel);  break;
    case KERNEL_Wendland:    Interaction_Mdbc2CorrectionT <KERNEL_Wendland  > (divdata,pos,code,idp,boundnor,motionvel,motionace,velrho,boundmode,tangenvel);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

//==============================================================================
/// Copy motion velocity to motionvel[] and computes acceleration in motionace[].
/// Copia velocidad de movimiento a motionvel[] y aceleracion en motionace[].
//==============================================================================
void JSphCpu::CopyMotionVelAce(unsigned nmoving,double dt,const unsigned* ridpmot
  ,const tfloat4* velrho,tfloat3* motionvel,tfloat3* motionace)const
{
  const int n=int(nmoving);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int id=0;id<n;id++){
    const unsigned pid=ridpmot[id];
    if(pid!=UINT_MAX){
      const tfloat3 mvel0=motionvel[pid];
      const tfloat4 v=velrho[pid];
      motionace[pid]=TFloat3(float((double(v.x)-mvel0.x)/dt),
                             float((double(v.y)-mvel0.y)/dt),
                             float((double(v.z)-mvel0.z)/dt));
      motionvel[pid]=TFloat3(v.x,v.y,v.z);
    }
  }
}
//<vs_m2dbc_end>

