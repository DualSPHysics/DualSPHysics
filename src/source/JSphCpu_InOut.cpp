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

/// \file JSphCpu_InOut.cpp \brief Implements InOut functions of class \ref JSphCpu.

#include "JSphCpu.h"
#include "JCellSearch_inline.h"
#include "FunSphKernel.h"
#include "FunctionsMath.h"
#include "FunGeo3d.h"
#include "JSphInOut.h"
#include "JSphInOutZone.h"
#include <climits>

using namespace std;

//==============================================================================
/// Returns original position of periodic particle.
//==============================================================================
tdouble3 JSphCpu::Interaction_PosNoPeriodic(tdouble3 posp1)const{
  if(PeriX){
    if(posp1.x<MapRealPosMin.x)posp1=posp1-PeriXinc;
    if(posp1.x>MapRealPosMax.x)posp1=posp1+PeriXinc;
  }
  if(PeriY){
    if(posp1.y<MapRealPosMin.y)posp1=posp1-PeriYinc;
    if(posp1.y>MapRealPosMax.y)posp1=posp1+PeriYinc;
  }
  if(PeriZ){
    if(posp1.z<MapRealPosMin.z)posp1=posp1-PeriZinc;
    if(posp1.z>MapRealPosMax.z)posp1=posp1+PeriZinc;
  }
  return(posp1);
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<bool sim2d,TpKernel tker> void JSphCpu::InteractionInOutExtrap_Double
  (unsigned inoutcount,const int *inoutpart,const byte *cfgzone
  ,const tplane3f *planes,const float* width,const tfloat3 *dirdata,float determlimit
  ,StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,const typecode *code
  ,const unsigned *idp,tfloat4 *velrhop)
{
  //Log->Printf("%u>++> InteractionInOutGhost_Double",Nstep);
  //-Inicia ejecucion con OpenMP.
  const int n=int(inoutcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const unsigned p1=(unsigned)inoutpart[p];
    if(!CODE_IsFluidInout(code[p1]))Run_Exceptioon("InOut particle is invalid.");
    const unsigned izone=CODE_GetIzoneFluidInout(code[p1]);
    const byte cfg=cfgzone[izone];
    const bool computerhop=(JSphInOutZone::GetConfigRhopMode(cfg)==InRhop_Extrapolated);
    const bool computevel=(JSphInOutZone::GetConfigVelMode(cfg)==InVelM_Extrapolated);
    if(computerhop || computevel){
      //-Calculates ghost node position.
      tdouble3 pos_p1=pos[p1];
      if(CODE_IsPeriodic(code[p1]))pos_p1=Interaction_PosNoPeriodic(pos_p1);
      const double displane=fgeo::PlaneDist(TPlane3d(planes[izone]),pos_p1)*2;
      const tdouble3 posp1=pos_p1+TDouble3(displane*dirdata[izone].x, displane*dirdata[izone].y, displane*dirdata[izone].z); //-Ghost node position.
  //Log->Printf("%u>++> idp[%u]:%u posx(p1,node): %f %f",Nstep,p1,idp[p1],pos_p1.x,posp1.x);
  //Log->Printf("%u>++> idp[%u]:%u izone:%u rhop:%u vel:%u",Nstep,p1,idp[p1],izone,(computerhop? 1: 0),(computevel? 1: 0));

      //-Initializes variables for calculation.
      double rhopp1=0;
      tdouble3 gradrhopp1=TDouble3(0);
      tdouble3 velp1=TDouble3(0);
      tmatrix3d gradvelp1=TMatrix3d(0); //-Only for velocity.
      tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
      tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(posp1,false,dvd);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
        //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
        //---------------------------------------------------------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const double drx=double(posp1.x-pos[p2].x);
          const double dry=double(posp1.y-pos[p2].y);
          const double drz=double(posp1.z-pos[p2].z);
          const double rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=KernelSize2 && rr2>=ALMOSTZERO && CODE_IsFluidNotInout(code[p2])){//-Only with fluid particles but not inout particles.
            //-Computes kernel.
            float fac;
            const double wab=fsph::GetKernel_WabFac<tker>(CSP,float(rr2),fac);
            const double frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.

            const tfloat4 velrhopp2=velrhop[p2];
            //===== Get mass and volume of particle p2 =====
            double massp2=MassFluid;
            double volp2=massp2/velrhopp2.w;

            //===== Density and its gradient =====
            rhopp1+=massp2*wab;
            gradrhopp1.x+=massp2*frx;
            gradrhopp1.y+=massp2*fry;
            gradrhopp1.z+=massp2*frz;

            //===== Kernel values multiplied by volume =====
            const double vwab=wab*volp2;
            const double vfrx=frx*volp2;
            const double vfry=fry*volp2;
            const double vfrz=frz*volp2;

            //===== Velocity and its gradient =====
            if(computevel){
              velp1.x+=vwab*velrhopp2.x;
              velp1.y+=vwab*velrhopp2.y;
              velp1.z+=vwab*velrhopp2.z;
              gradvelp1.a11+=vfrx*velrhopp2.x;    // du/dx
              gradvelp1.a12+=vfry*velrhopp2.x;    // du/dy
              gradvelp1.a13+=vfrz*velrhopp2.x;    // du/dz
              gradvelp1.a21+=vfrx*velrhopp2.y;    // dv/dx
              gradvelp1.a22+=vfry*velrhopp2.y;    // dv/dx
              gradvelp1.a23+=vfrz*velrhopp2.y;    // dv/dx
              gradvelp1.a31+=vfrx*velrhopp2.z;    // dw/dx
              gradvelp1.a32+=vfry*velrhopp2.z;    // dw/dx
              gradvelp1.a33+=vfrz*velrhopp2.z;    // dw/dx
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
      tfloat4 velrhopfinal=velrhop[p1];
      const tdouble3 dpos=(pos_p1-posp1); //-Inlet/outlet particle position - ghost node position.
      if(sim2d){
        const double determ=fmath::Determinant3x3(a_corr2);
        if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
          const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);    //CHECKED
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          if(computerhop){
            const double rhoghost=rhopp1*invacorr2.a11 + gradrhopp1.x*invacorr2.a12 + gradrhopp1.z*invacorr2.a13;
            const double grx=-(rhopp1*invacorr2.a21 + gradrhopp1.x*invacorr2.a22 + gradrhopp1.z*invacorr2.a23);
            const double grz=-(rhopp1*invacorr2.a31 + gradrhopp1.x*invacorr2.a32 + gradrhopp1.z*invacorr2.a33);
            velrhopfinal.w=float(rhoghost + grx*dpos.x + grz*dpos.z);
          }
          //-GHOST NODE VELOCITY ARE MIRRORED BACK TO THE OUTFLOW PARTICLES.
          if(computevel){
            const double velghost_x=velp1.x*invacorr2.a11 + gradvelp1.a11*invacorr2.a12 + gradvelp1.a13*invacorr2.a13;
            const double velghost_z=velp1.z*invacorr2.a11 + gradvelp1.a31*invacorr2.a12 + gradvelp1.a33*invacorr2.a13;
            const double a11=-(velp1.x*invacorr2.a21 + gradvelp1.a11*invacorr2.a22 + gradvelp1.a13*invacorr2.a23);
            const double a13=-(velp1.z*invacorr2.a21 + gradvelp1.a31*invacorr2.a22 + gradvelp1.a33*invacorr2.a23);
            const double a31=-(velp1.x*invacorr2.a31 + gradvelp1.a11*invacorr2.a32 + gradvelp1.a13*invacorr2.a33);
            const double a33=-(velp1.z*invacorr2.a31 + gradvelp1.a31*invacorr2.a32 + gradvelp1.a33*invacorr2.a33);
            velrhopfinal.x=float(velghost_x + a11*dpos.x + a31*dpos.z);
            velrhopfinal.z=float(velghost_z + a13*dpos.x + a33*dpos.z);
            velrhopfinal.y=0;
          }
        }
        else if(a_corr2.a11>0){ // Determinant is small but a11 is nonzero, 0th order ANGELO
          if(computerhop)velrhopfinal.w=float(rhopp1/a_corr2.a11);
          if(computevel){
            velrhopfinal.x=float(velp1.x/a_corr2.a11);
            velrhopfinal.z=float(velp1.z/a_corr2.a11);
            velrhopfinal.y=0;
          }
        }
      }
      else{
        const double determ=fmath::Determinant4x4(a_corr3);
        if(fabs(determ)>=determlimit){
          const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          if(computerhop){
            const double rhoghost=rhopp1*invacorr3.a11 + gradrhopp1.x*invacorr3.a12 + gradrhopp1.y*invacorr3.a13 + gradrhopp1.z*invacorr3.a14;
            const double grx=   -(rhopp1*invacorr3.a21 + gradrhopp1.x*invacorr3.a22 + gradrhopp1.y*invacorr3.a23 + gradrhopp1.z*invacorr3.a24);
            const double gry=   -(rhopp1*invacorr3.a31 + gradrhopp1.x*invacorr3.a32 + gradrhopp1.y*invacorr3.a33 + gradrhopp1.z*invacorr3.a34);
            const double grz=   -(rhopp1*invacorr3.a41 + gradrhopp1.x*invacorr3.a42 + gradrhopp1.y*invacorr3.a43 + gradrhopp1.z*invacorr3.a44);
            velrhopfinal.w=float(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
          }
          //-GHOST NODE VELOCITY ARE MIRRORED BACK TO THE OUTFLOW PARTICLES.
          if(computevel){
            const double velghost_x=velp1.x*invacorr3.a11 + gradvelp1.a11*invacorr3.a12 + gradvelp1.a12*invacorr3.a13 + gradvelp1.a13*invacorr3.a14;
            const double velghost_y=velp1.y*invacorr3.a11 + gradvelp1.a11*invacorr3.a12 + gradvelp1.a12*invacorr3.a13 + gradvelp1.a13*invacorr3.a14;
            const double velghost_z=velp1.z*invacorr3.a11 + gradvelp1.a31*invacorr3.a12 + gradvelp1.a32*invacorr3.a13 + gradvelp1.a33*invacorr3.a14;
            const double a11=-(velp1.x*invacorr3.a21 + gradvelp1.a11*invacorr3.a22 + gradvelp1.a12*invacorr3.a23 + gradvelp1.a13*invacorr3.a24);
            const double a12=-(velp1.y*invacorr3.a21 + gradvelp1.a21*invacorr3.a22 + gradvelp1.a22*invacorr3.a23 + gradvelp1.a23*invacorr3.a24);
            const double a13=-(velp1.z*invacorr3.a21 + gradvelp1.a31*invacorr3.a22 + gradvelp1.a32*invacorr3.a23 + gradvelp1.a33*invacorr3.a24);
            const double a21=-(velp1.x*invacorr3.a31 + gradvelp1.a11*invacorr3.a32 + gradvelp1.a12*invacorr3.a33 + gradvelp1.a13*invacorr3.a34);
            const double a22=-(velp1.y*invacorr3.a31 + gradvelp1.a21*invacorr3.a32 + gradvelp1.a22*invacorr3.a33 + gradvelp1.a23*invacorr3.a34);
            const double a23=-(velp1.z*invacorr3.a31 + gradvelp1.a31*invacorr3.a32 + gradvelp1.a32*invacorr3.a33 + gradvelp1.a33*invacorr3.a34);
            const double a31=-(velp1.x*invacorr3.a41 + gradvelp1.a11*invacorr3.a42 + gradvelp1.a12*invacorr3.a43 + gradvelp1.a13*invacorr3.a44);
            const double a32=-(velp1.y*invacorr3.a41 + gradvelp1.a21*invacorr3.a42 + gradvelp1.a22*invacorr3.a43 + gradvelp1.a23*invacorr3.a44);
            const double a33=-(velp1.z*invacorr3.a41 + gradvelp1.a31*invacorr3.a42 + gradvelp1.a32*invacorr3.a43 + gradvelp1.a33*invacorr3.a44);
            velrhopfinal.x=float(velghost_x + a11*dpos.x + a21*dpos.y + a31*dpos.z);
            velrhopfinal.y=float(velghost_y + a12*dpos.x + a22*dpos.y + a32*dpos.z);
            velrhopfinal.z=float(velghost_z + a13*dpos.x + a23*dpos.y + a33*dpos.z);
          }
        }
        else if(a_corr3.a11>0){ // Determinant is small but a11 is nonzero, 0th order ANGELO
          if(computerhop)velrhopfinal.w=float(rhopp1/a_corr3.a11);
          if(computevel){
            velrhopfinal.x=float(velp1.x/a_corr3.a11);
            velrhopfinal.y=float(velp1.y/a_corr3.a11);
            velrhopfinal.z=float(velp1.z/a_corr3.a11);
          }
        }
      }
      velrhop[p1]=velrhopfinal;
    }
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<bool sim2d,TpKernel tker> void JSphCpu::InteractionInOutExtrap_Single
  (unsigned inoutcount,const int *inoutpart,const byte *cfgzone
  ,const tplane3f *planes,const float* width,const tfloat3 *dirdata,float determlimit
  ,StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,const typecode *code
  ,const unsigned *idp,tfloat4 *velrhop)
{
  //-Inicia ejecucion con OpenMP.
  const int n=int(inoutcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const unsigned p1=(unsigned)inoutpart[p];
    if(!CODE_IsFluidInout(code[p1]))Run_Exceptioon("InOut particle is invalid.");
    const unsigned izone=CODE_GetIzoneFluidInout(code[p1]);
    const byte cfg=cfgzone[izone];
    const bool computerhop=(JSphInOutZone::GetConfigRhopMode(cfg)==InRhop_Extrapolated);
    const bool computevel=(JSphInOutZone::GetConfigVelMode(cfg)==InVelM_Extrapolated);
    if(computerhop || computevel){
      //-Calculates ghost node position.
      tdouble3 pos_p1=pos[p1];
      if(CODE_IsPeriodic(code[p1]))pos_p1=Interaction_PosNoPeriodic(pos_p1);
      const double displane=fgeo::PlaneDist(TPlane3d(planes[izone]),pos_p1)*2;
      const tdouble3 posp1=pos_p1+TDouble3(displane*dirdata[izone].x, displane*dirdata[izone].y, displane*dirdata[izone].z); //-Ghost node position.

      //-Initializes variables for calculation.
      float rhopp1=0;
      tfloat3 gradrhopp1=TFloat3(0);
      tfloat3 velp1=TFloat3(0);
      tmatrix3f gradvelp1=TMatrix3f(0); //-Only for velocity.
      tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
      tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(posp1,false,dvd);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
        //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
        //---------------------------------------------------------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(posp1.x-pos[p2].x);
          const float dry=float(posp1.y-pos[p2].y);
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=KernelSize2 && rr2>=ALMOSTZERO && CODE_IsFluidNotInout(code[p2])){//-Only with fluid particles but not inout particles.
            //-Computes kernel.
            float fac;
            const float wab=fsph::GetKernel_WabFac<tker>(CSP,rr2,fac);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

            const tfloat4 velrhopp2=velrhop[p2];
            //===== Get mass and volume of particle p2 =====
            float massp2=MassFluid;
            float volp2=massp2/velrhopp2.w;

            //===== Density and its gradient =====
            rhopp1+=massp2*wab;
            gradrhopp1.x+=massp2*frx;
            gradrhopp1.y+=massp2*fry;
            gradrhopp1.z+=massp2*frz;

            //===== Kernel values multiplied by volume =====
            const float vwab=wab*volp2;
            const float vfrx=frx*volp2;
            const float vfry=fry*volp2;
            const float vfrz=frz*volp2;

            //===== Velocity and its gradient =====
            if(computevel){
              velp1.x+=vwab*velrhopp2.x;
              velp1.y+=vwab*velrhopp2.y;
              velp1.z+=vwab*velrhopp2.z;
              gradvelp1.a11+=vfrx*velrhopp2.x;    // du/dx
              gradvelp1.a12+=vfry*velrhopp2.x;    // du/dy
              gradvelp1.a13+=vfrz*velrhopp2.x;    // du/dz
              gradvelp1.a21+=vfrx*velrhopp2.y;    // dv/dx
              gradvelp1.a22+=vfry*velrhopp2.y;    // dv/dx
              gradvelp1.a23+=vfrz*velrhopp2.y;    // dv/dx
              gradvelp1.a31+=vfrx*velrhopp2.z;    // dw/dx
              gradvelp1.a32+=vfry*velrhopp2.z;    // dw/dx
              gradvelp1.a33+=vfrz*velrhopp2.z;    // dw/dx
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
      tfloat4 velrhopfinal=velrhop[p1];
      const tfloat3 dpos=ToTFloat3(pos_p1-posp1); //-Inlet/outlet particle position - ghost node position.
      if(sim2d){
        const double determ=fmath::Determinant3x3(a_corr2);
        if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
          const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);    //CHECKED
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          if(computerhop){
            const float rhoghost=float(invacorr2.a11*rhopp1 + invacorr2.a12*gradrhopp1.x + invacorr2.a13*gradrhopp1.z);
            const float grx=    -float(invacorr2.a21*rhopp1 + invacorr2.a22*gradrhopp1.x + invacorr2.a23*gradrhopp1.z);
            const float grz=    -float(invacorr2.a31*rhopp1 + invacorr2.a32*gradrhopp1.x + invacorr2.a33*gradrhopp1.z);
            velrhopfinal.w=(rhoghost + grx*dpos.x + grz*dpos.z);
          }
          //-GHOST NODE VELOCITY ARE MIRRORED BACK TO THE OUTFLOW PARTICLES.
          if(computevel){
            const float velghost_x=float(invacorr2.a11*velp1.x + invacorr2.a12*gradvelp1.a11 + invacorr2.a13*gradvelp1.a13);
            const float velghost_z=float(invacorr2.a11*velp1.z + invacorr2.a12*gradvelp1.a31 + invacorr2.a13*gradvelp1.a33);
            const float a11=-float(invacorr2.a21*velp1.x + invacorr2.a22*gradvelp1.a11 + invacorr2.a23*gradvelp1.a13);
            const float a13=-float(invacorr2.a21*velp1.z + invacorr2.a22*gradvelp1.a31 + invacorr2.a23*gradvelp1.a33);
            const float a31=-float(invacorr2.a31*velp1.x + invacorr2.a32*gradvelp1.a11 + invacorr2.a33*gradvelp1.a13);
            const float a33=-float(invacorr2.a31*velp1.z + invacorr2.a32*gradvelp1.a31 + invacorr2.a33*gradvelp1.a33);
            velrhopfinal.x=(velghost_x + a11*dpos.x + a31*dpos.z);
            velrhopfinal.z=(velghost_z + a13*dpos.x + a33*dpos.z);
            velrhopfinal.y=0;
          }
        }
        else if(a_corr2.a11>0){ // Determinant is small but a11 is nonzero, 0th order ANGELO
          if(computerhop)velrhopfinal.w=float(rhopp1/a_corr2.a11);
          if(computevel){
            velrhopfinal.x=float(velp1.x/a_corr2.a11);
            velrhopfinal.z=float(velp1.z/a_corr2.a11);
            velrhopfinal.y=0;
          }
        }
      }
      else{
        const double determ=fmath::Determinant4x4(a_corr3);
        if(fabs(determ)>=determlimit){
          const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          if(computerhop){
            const float rhoghost=float(invacorr3.a11*rhopp1 + invacorr3.a12*gradrhopp1.x + invacorr3.a13*gradrhopp1.y + invacorr3.a14*gradrhopp1.z);
            const float grx=    -float(invacorr3.a21*rhopp1 + invacorr3.a22*gradrhopp1.x + invacorr3.a23*gradrhopp1.y + invacorr3.a24*gradrhopp1.z);
            const float gry=    -float(invacorr3.a31*rhopp1 + invacorr3.a32*gradrhopp1.x + invacorr3.a33*gradrhopp1.y + invacorr3.a34*gradrhopp1.z);
            const float grz=    -float(invacorr3.a41*rhopp1 + invacorr3.a42*gradrhopp1.x + invacorr3.a43*gradrhopp1.y + invacorr3.a44*gradrhopp1.z);
            velrhopfinal.w=(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
          }
          //-GHOST NODE VELOCITY ARE MIRRORED BACK TO THE OUTFLOW PARTICLES.
          if(computevel){
            const float velghost_x=float(invacorr3.a11*velp1.x + invacorr3.a12*gradvelp1.a11 + invacorr3.a13*gradvelp1.a12 + invacorr3.a14*gradvelp1.a13);
            const float velghost_y=float(invacorr3.a11*velp1.y + invacorr3.a12*gradvelp1.a11 + invacorr3.a13*gradvelp1.a12 + invacorr3.a14*gradvelp1.a13);
            const float velghost_z=float(invacorr3.a11*velp1.z + invacorr3.a12*gradvelp1.a31 + invacorr3.a13*gradvelp1.a32 + invacorr3.a14*gradvelp1.a33);
            const float a11=      -float(invacorr3.a21*velp1.x + invacorr3.a22*gradvelp1.a11 + invacorr3.a23*gradvelp1.a12 + invacorr3.a24*gradvelp1.a13);
            const float a12=      -float(invacorr3.a21*velp1.y + invacorr3.a22*gradvelp1.a21 + invacorr3.a23*gradvelp1.a22 + invacorr3.a24*gradvelp1.a23);
            const float a13=      -float(invacorr3.a21*velp1.z + invacorr3.a22*gradvelp1.a31 + invacorr3.a23*gradvelp1.a32 + invacorr3.a24*gradvelp1.a33);
            const float a21=      -float(invacorr3.a31*velp1.x + invacorr3.a32*gradvelp1.a11 + invacorr3.a33*gradvelp1.a12 + invacorr3.a34*gradvelp1.a13);
            const float a22=      -float(invacorr3.a31*velp1.y + invacorr3.a32*gradvelp1.a21 + invacorr3.a33*gradvelp1.a22 + invacorr3.a34*gradvelp1.a23);
            const float a23=      -float(invacorr3.a31*velp1.z + invacorr3.a32*gradvelp1.a31 + invacorr3.a33*gradvelp1.a32 + invacorr3.a34*gradvelp1.a33);
            const float a31=      -float(invacorr3.a41*velp1.x + invacorr3.a42*gradvelp1.a11 + invacorr3.a43*gradvelp1.a12 + invacorr3.a44*gradvelp1.a13);
            const float a32=      -float(invacorr3.a41*velp1.y + invacorr3.a42*gradvelp1.a21 + invacorr3.a43*gradvelp1.a22 + invacorr3.a44*gradvelp1.a23);
            const float a33=      -float(invacorr3.a41*velp1.z + invacorr3.a42*gradvelp1.a31 + invacorr3.a43*gradvelp1.a32 + invacorr3.a44*gradvelp1.a33);
            velrhopfinal.x=(velghost_x + a11*dpos.x + a21*dpos.y + a31*dpos.z);
            velrhopfinal.y=(velghost_y + a12*dpos.x + a22*dpos.y + a32*dpos.z);
            velrhopfinal.z=(velghost_z + a13*dpos.x + a23*dpos.y + a33*dpos.z);
          }
        }
        else if(a_corr3.a11>0){ // Determinant is small but a11 is nonzero, 0th order ANGELO
          if(computerhop)velrhopfinal.w=float(rhopp1/a_corr3.a11);
          if(computevel){
            velrhopfinal.x=float(velp1.x/a_corr3.a11);
            velrhopfinal.y=float(velp1.y/a_corr3.a11);
            velrhopfinal.z=float(velp1.z/a_corr3.a11);
          }
        }
      }
      velrhop[p1]=velrhopfinal;
    }
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker> void JSphCpu::Interaction_InOutExtrapT(byte doublemode,unsigned inoutcount,const int *inoutpart
  ,const byte *cfgzone,const tplane3f *planes
  ,const float* width,const tfloat3 *dirdata,float determlimit
  ,const unsigned *dcell,const tdouble3 *pos,const typecode *code
  ,const unsigned *idp,tfloat4 *velrhop)
{
  const StDivDataCpu &dvd=DivData;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(doublemode==2){
    if(Simulate2D)InteractionInOutExtrap_Single<true ,tker> (inoutcount,inoutpart,cfgzone,planes,width,dirdata,determlimit,dvd,dcell,pos,code,idp,velrhop);
    else          InteractionInOutExtrap_Single<false,tker> (inoutcount,inoutpart,cfgzone,planes,width,dirdata,determlimit,dvd,dcell,pos,code,idp,velrhop);
  }
  else if(doublemode==3){
    if(Simulate2D)InteractionInOutExtrap_Double<true ,tker> (inoutcount,inoutpart,cfgzone,planes,width,dirdata,determlimit,dvd,dcell,pos,code,idp,velrhop);
    else          InteractionInOutExtrap_Double<false,tker> (inoutcount,inoutpart,cfgzone,planes,width,dirdata,determlimit,dvd,dcell,pos,code,idp,velrhop);
  }
  else Run_Exceptioon("Double mode calculation is invalid.");
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void JSphCpu::Interaction_InOutExtrap(byte doublemode,unsigned inoutcount,const int *inoutpart
  ,const byte *cfgzone,const tplane3f *planes
  ,const float* width,const tfloat3 *dirdata,float determlimit
  ,const unsigned *dcell,const tdouble3 *pos,const typecode *code
  ,const unsigned *idp,tfloat4 *velrhop)
{
  switch(TKernel){
    case KERNEL_Cubic:       Interaction_InOutExtrapT<KERNEL_Cubic>     (doublemode,inoutcount,inoutpart,cfgzone,planes,width,dirdata,determlimit,dcell,pos,code,idp,velrhop);  break;
    case KERNEL_Wendland:    Interaction_InOutExtrapT<KERNEL_Wendland>  (doublemode,inoutcount,inoutpart,cfgzone,planes,width,dirdata,determlimit,dcell,pos,code,idp,velrhop);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

//==============================================================================
/// Calculates maximum zsurf in fluid domain.
/// Calcula zsurf maximo en el fluido.
//==============================================================================
float JSphCpu::Interaction_InOutZsurf(unsigned nptz,const tfloat3 *ptzpos,float maxdist,float zbottom
  ,const StDivDataCpu &divdata,const tdouble3 *pos,const typecode *code)
{
  const float maxdist2=maxdist*maxdist;
//  Log->Printf("%u>++> InteractionInOutGhost_Double",Nstep);
  float zsurfmax=zbottom;
  for(unsigned p1=0;p1<nptz;p1++){
    const tfloat3 posp1=ptzpos[p1];

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(ToTDouble3(posp1),false,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);
      //-Interaction of boundary with type Fluid/Float.
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const tfloat3 posp2=ToTFloat3(pos[p2]);
        if(posp2.z>zsurfmax){
          const float drx=posp1.x-posp2.x;
          const float dry=posp1.y-posp2.y;
          const float drz=posp1.z-posp2.z;
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=maxdist2 && CODE_IsFluidNotInout(code[p2])){//-Only with fluid particles but not inout particles.
            zsurfmax=posp2.z;
          }
        }
      }
    }
  }
  return(zsurfmax);
}

//==============================================================================
/// Perform interaction between ghost node of selected boundary and fluid.
//==============================================================================
template<bool sim2d,TpKernel tker> void JSphCpu::InteractionBoundCorr_Double
  (unsigned npb,typecode boundcode,tplane3f plane,tfloat3 direction,float determlimit
  ,const StDivDataCpu &dvd,const tdouble3 *pos,const typecode *code
  ,const unsigned *idp,tfloat4 *velrhop)
{
  const int n=int(npb);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<n;p1++)if(CODE_GetTypeAndValue(code[p1])==boundcode){
    float rhopfinal=FLT_MAX;
    //-Calculates ghost node position.
    tdouble3 pos_p1=pos[p1];
    if(CODE_IsPeriodic(code[p1]))pos_p1=Interaction_PosNoPeriodic(pos_p1);
    const double displane=fgeo::PlaneDist(TPlane3d(plane),pos_p1)*2;
    if(displane<=KernelSize*2.f){
      const tdouble3 posp1=pos_p1+TDouble3(displane*direction.x,displane*direction.y,displane*direction.z); //-Ghost node position.
      //-Initializes variables for calculation.
      double rhopp1=0;
      tdouble3 gradrhopp1=TDouble3(0);
      tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
      tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(posp1,false,dvd);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
        //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
        //---------------------------------------------------------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const double drx=double(posp1.x-pos[p2].x);
          const double dry=double(posp1.y-pos[p2].y);
          const double drz=double(posp1.z-pos[p2].z);
          const double rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=KernelSize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
            //-Computes kernel.
            float fac;
            const double wab=fsph::GetKernel_WabFac<tker>(CSP,float(rr2),fac);
            const double frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.

            //===== Get mass and volume of particle p2 =====
            const double massp2=MassFluid;
            const double volp2=massp2/double(velrhop[p2].w);

            //===== Density and its gradient =====
            rhopp1+=massp2*wab;
            gradrhopp1.x+=massp2*frx;
            gradrhopp1.y+=massp2*fry;
            gradrhopp1.z+=massp2*frz;

            //===== Kernel values multiplied by volume =====
            const double vwab=wab*volp2;
            const double vfrx=frx*volp2;
            const double vfry=fry*volp2;
            const double vfrz=frz*volp2;

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
      const tdouble3 dpos=(pos_p1-posp1); //-Boundary particle position - ghost node position.
      if(sim2d){
        const double determ=fmath::Determinant3x3(a_corr2);
        if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
          const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          const double rhoghost=rhopp1*invacorr2.a11 + gradrhopp1.x*invacorr2.a12 + gradrhopp1.z*invacorr2.a13;
          const double grx=   -(rhopp1*invacorr2.a21 + gradrhopp1.x*invacorr2.a22 + gradrhopp1.z*invacorr2.a23);
          const double grz=   -(rhopp1*invacorr2.a31 + gradrhopp1.x*invacorr2.a32 + gradrhopp1.z*invacorr2.a33);
          rhopfinal=float(rhoghost + grx*dpos.x + grz*dpos.z);
        }
        else if(a_corr2.a11>0){//-Determinant is small but a11 is nonzero, 0th order ANGELO.
          rhopfinal=float(rhopp1/a_corr2.a11);
        }
      }
      else{
        const double determ=fmath::Determinant4x4(a_corr3);
        if(fabs(determ)>=determlimit){
          const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          const double rhoghost=rhopp1*invacorr3.a11 + gradrhopp1.x*invacorr3.a12 + gradrhopp1.y*invacorr3.a13 + gradrhopp1.z*invacorr3.a14;
          const double grx=   -(rhopp1*invacorr3.a21 + gradrhopp1.x*invacorr3.a22 + gradrhopp1.y*invacorr3.a23 + gradrhopp1.z*invacorr3.a24);
          const double gry=   -(rhopp1*invacorr3.a31 + gradrhopp1.x*invacorr3.a32 + gradrhopp1.y*invacorr3.a33 + gradrhopp1.z*invacorr3.a34);
          const double grz=   -(rhopp1*invacorr3.a41 + gradrhopp1.x*invacorr3.a42 + gradrhopp1.y*invacorr3.a43 + gradrhopp1.z*invacorr3.a44);
          rhopfinal=float(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
        }
        else if(a_corr3.a11>0){//-Determinant is small but a11 is nonzero, 0th order ANGELO.
          rhopfinal=float(rhopp1/a_corr3.a11);
        }
      }
    }
    velrhop[p1].w=(rhopfinal!=FLT_MAX? rhopfinal: RhopZero);
  }
}

//==============================================================================
/// Perform interaction between ghost node of selected boundary and fluid.
//==============================================================================
template<bool sim2d,TpKernel tker> void JSphCpu::InteractionBoundCorr_Single
  (unsigned npb,typecode boundcode,tplane3f plane,tfloat3 direction,float determlimit
  ,const StDivDataCpu &dvd,const tdouble3 *pos,const typecode *code
  ,const unsigned *idp,tfloat4 *velrhop)
{
  const int n=int(npb);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<n;p1++)if(CODE_GetTypeAndValue(code[p1])==boundcode){
    float rhopfinal=FLT_MAX;
    //-Calculates ghost node position.
    tdouble3 pos_p1=pos[p1];
    if(CODE_IsPeriodic(code[p1]))pos_p1=Interaction_PosNoPeriodic(pos_p1);
    const double displane=fgeo::PlaneDist(TPlane3d(plane),pos_p1)*2;
    if(displane<=KernelSize*2.f){
      const tdouble3 posp1=pos_p1+TDouble3(displane*direction.x,displane*direction.y,displane*direction.z); //-Ghost node position.
      //-Initializes variables for calculation.
      float rhopp1=0;
      tfloat3 gradrhopp1=TFloat3(0);
      tmatrix3d a_corr2=TMatrix3d(0);   //-Only for 2D.
      tmatrix4d a_corr3=TMatrix4d(0);   //-Only for 3D.

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(posp1,false,dvd);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
        //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
        //---------------------------------------------------------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(posp1.x-pos[p2].x);
          const float dry=float(posp1.y-pos[p2].y);
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=KernelSize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){//-Only with fluid particles (including inout).
            //-Computes kernel.
            float fac;
            const float wab=fsph::GetKernel_WabFac<tker>(CSP,rr2,fac);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

            //===== Get mass and volume of particle p2 =====
            const float massp2=MassFluid;
            const float volp2=massp2/velrhop[p2].w;

            //===== Density and its gradient =====
            rhopp1+=massp2*wab;
            gradrhopp1.x+=massp2*frx;
            gradrhopp1.y+=massp2*fry;
            gradrhopp1.z+=massp2*frz;

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
      const tfloat3 dpos=ToTFloat3(pos_p1-posp1); //-Boundary particle position - ghost node position.
      if(sim2d){
        const double determ=fmath::Determinant3x3(a_corr2);
        if(fabs(determ)>=determlimit){//-Use 1e-3f (first_order) or 1e+3f (zeroth_order).
          const tmatrix3d invacorr2=fmath::InverseMatrix3x3(a_corr2,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          const float rhoghost=float(invacorr2.a11*rhopp1 + invacorr2.a12*gradrhopp1.x + invacorr2.a13*gradrhopp1.z);
          const float grx=    -float(invacorr2.a21*rhopp1 + invacorr2.a22*gradrhopp1.x + invacorr2.a23*gradrhopp1.z);
          const float grz=    -float(invacorr2.a31*rhopp1 + invacorr2.a32*gradrhopp1.x + invacorr2.a33*gradrhopp1.z);
          rhopfinal=(rhoghost + grx*dpos.x + grz*dpos.z);
        }
        else if(a_corr2.a11>0){//-Determinant is small but a11 is nonzero, 0th order ANGELO.
          rhopfinal=float(rhopp1/a_corr2.a11);
        }
      }
      else{
        const double determ=fmath::Determinant4x4(a_corr3);
        if(fabs(determ)>=determlimit){
          const tmatrix4d invacorr3=fmath::InverseMatrix4x4(a_corr3,determ);
          //-GHOST NODE DENSITY IS MIRRORED BACK TO THE INFLOW OR OUTFLOW PARTICLES.
          const float rhoghost=float(invacorr3.a11*rhopp1 + invacorr3.a12*gradrhopp1.x + invacorr3.a13*gradrhopp1.y + invacorr3.a14*gradrhopp1.z);
          const float grx=    -float(invacorr3.a21*rhopp1 + invacorr3.a22*gradrhopp1.x + invacorr3.a23*gradrhopp1.y + invacorr3.a24*gradrhopp1.z);
          const float gry=    -float(invacorr3.a31*rhopp1 + invacorr3.a32*gradrhopp1.x + invacorr3.a33*gradrhopp1.y + invacorr3.a34*gradrhopp1.z);
          const float grz=    -float(invacorr3.a41*rhopp1 + invacorr3.a42*gradrhopp1.x + invacorr3.a43*gradrhopp1.y + invacorr3.a44*gradrhopp1.z);
          rhopfinal=(rhoghost + grx*dpos.x + gry*dpos.y + grz*dpos.z);
        }
        else if(a_corr3.a11>0){//-Determinant is small but a11 is nonzero, 0th order ANGELO.
          rhopfinal=float(rhopp1/a_corr3.a11);
        }
      }
    }
    velrhop[p1].w=(rhopfinal!=FLT_MAX? rhopfinal: RhopZero);
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker> void JSphCpu::Interaction_BoundCorrT
  (byte doublemode,typecode boundcode,tplane3f plane,tfloat3 direction,float determlimit
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop)
{
  const StDivDataCpu &dvd=DivData;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(doublemode==2){
    if(Simulate2D)InteractionBoundCorr_Single<true ,tker>(NpbOk,boundcode,plane,direction,determlimit,dvd,pos,code,idp,velrhop);
    else          InteractionBoundCorr_Single<false,tker>(NpbOk,boundcode,plane,direction,determlimit,dvd,pos,code,idp,velrhop);
  }
  else if(doublemode==3){
    if(Simulate2D)InteractionBoundCorr_Double<true ,tker>(NpbOk,boundcode,plane,direction,determlimit,dvd,pos,code,idp,velrhop);
    else          InteractionBoundCorr_Double<false,tker>(NpbOk,boundcode,plane,direction,determlimit,dvd,pos,code,idp,velrhop);
  }
  else Run_Exceptioon("Double mode calculation is invalid.");
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void JSphCpu::Interaction_BoundCorr(byte doublemode,typecode boundcode
  ,tplane3f plane,tfloat3 direction,float determlimit
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp
  ,tfloat4 *velrhop)
{
  switch(TKernel){
    case KERNEL_Cubic:       Interaction_BoundCorrT<KERNEL_Cubic>     (doublemode,boundcode,plane,direction,determlimit,pos,code,idp,velrhop);  break;
    case KERNEL_Wendland:    Interaction_BoundCorrT<KERNEL_Wendland>  (doublemode,boundcode,plane,direction,determlimit,pos,code,idp,velrhop);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

