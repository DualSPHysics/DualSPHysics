//HEAD_DSCODES
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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Creacion para implementacion de kernels SPH. (03-07-2020)
//:#############################################################################

/// \file FunSphKernel.h \brief Defines inline functions for SPH kenernels.

#ifndef _FunSphKernel_
#define _FunSphKernel_

#include "DualSphDef.h"
#include "FunSphKernelDef.h"
#include <cmath>

/// Implements a set of basic/general functions related to SPH.
namespace fsph{

//##############################################################################
//# Cubic Spline kernel
//##############################################################################
//==============================================================================
/// Returns the name of the kernel in text format.
//==============================================================================
inline const char* GetKernelCubic_Name(){ return("CubicSpline"); }
//==============================================================================
/// Returns factor value to compute KernelSize according to KernelH.
//==============================================================================
inline float GetKernelCubic_Factor(){ return(2.0f); }
//==============================================================================
/// Returns constants for kernel calculation.
//==============================================================================
inline StKCubicCte GetKernelCubic_Ctes(bool sim2d,double h){
  StKCubicCte kc={0,0,0,0,0,0,0,0};
  if(sim2d){
    const double a1=10./(PI*7.);
    const double a2=a1/(h*h);
    const double aa=a1/(h*h*h);
    const double deltap=1./1.5;
    const double wdeltap=a2*(1.-1.5*deltap*deltap+0.75*deltap*deltap*deltap);
    kc.od_wdeltap=float(1./wdeltap);
    kc.a1=float(a1);
    kc.a2=float(a2);
    kc.aa=float(aa);
    kc.a24=float(0.25*a2);
    kc.c1=float(-3.*aa);
    kc.d1=float(9.*aa/4.);
    kc.c2=float(-3.*aa/4.);
  }
  else{
    const double a1=1./PI;
    const double a2=a1/(h*h*h);
    const double aa=a1/(h*h*h*h);
    const double deltap=1./1.5;
    const double wdeltap=a2*(1.-1.5*deltap*deltap+0.75*deltap*deltap*deltap);
    kc.od_wdeltap=float(1./wdeltap);
    kc.a1=float(a1);
    kc.a2=float(a2);
    kc.aa=float(aa);
    kc.a24=float(0.25*a2);
    kc.c1=float(-3.*aa);
    kc.d1=float(9.*aa/4.);
    kc.c2=float(-3.*aa/4.);
  }
  return(kc);
}
//============================================================================== 
/// Returns wab of kernel.
//==============================================================================
inline float GetKernelCubic_Wab(const StKCubicCte &kc,float h,float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  if(rad>h){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    return(kc.a24*(wqq2*wqq1));
  }
  else{
    const float wqq2=qq*qq;
    return(kc.a2*(1.0f+(0.75f*qq-1.5f)*wqq2));
  }
}
//============================================================================== 
/// Returns fac of kernel.
//==============================================================================
inline float GetKernelCubic_Fac(const StKCubicCte &kc,float h,float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  if(rad>h){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    return(kc.c2*wqq2/rad);
  }
  else{
    const float wqq2=qq*qq;
    return((kc.c1*qq+kc.d1*wqq2)/rad);
  }
}
//============================================================================== 
/// Returns wab and fac of kernel.
//==============================================================================
inline float GetKernelCubic_WabFac(const StKCubicCte &kc,float h,float rr2,float &fac){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  if(rad>h){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    fac=(kc.c2*wqq2/rad);
    return(kc.a24*(wqq2*wqq1));
  }
  else{
    float wqq2=qq*qq;
    fac=((kc.c1*qq+kc.d1*wqq2)/rad);
    return(kc.a2*(1.0f+(0.75f*qq-1.5f)*wqq2));
  }
}
//==============================================================================
/// Return tensil correction of kernel Cubic.
//==============================================================================
inline float GetKernelCubic_Tensil(const StKCubicCte &kc,float h,float rr2
  ,float rhopp1,float pressp1,float rhopp2,float pressp2)
{
  const float wab=GetKernelCubic_Wab(kc,h,rr2);
  float fab=wab*kc.od_wdeltap;
  fab*=fab; fab*=fab; //fab=fab^4
  const float tensilp1=(pressp1/(rhopp1*rhopp1))*(pressp1>0? 0.01f: -0.2f);
  const float tensilp2=(pressp2/(rhopp2*rhopp2))*(pressp2>0? 0.01f: -0.2f);
  return(fab*(tensilp1+tensilp2));
}

//============================================================================== 
/// Returns wab of kernel.
//==============================================================================
inline float GetKernelCubic_Wab(const StCteSph &csp,float rr2){
  return(GetKernelCubic_Wab(csp.kcubic,csp.kernelh,rr2));
}
//============================================================================== 
/// Returns fac of kernel.
//==============================================================================
inline float GetKernelCubic_Fac(const StCteSph &csp,float rr2){
  return(GetKernelCubic_Fac(csp.kcubic,csp.kernelh,rr2));
}
//============================================================================== 
/// Returns wab and fac of kernel.
//==============================================================================
inline float GetKernelCubic_WabFac(const StCteSph &csp,float rr2,float &fac){
  return(GetKernelCubic_WabFac(csp.kcubic,csp.kernelh,rr2,fac));
}
//============================================================================== 
/// Return tensil correction of kernel Cubic.
//==============================================================================
inline float GetKernelCubic_Tensil(const StCteSph &csp,float rr2
  ,float rhopp1,float pressp1,float rhopp2,float pressp2)
{
  return(GetKernelCubic_Tensil(csp.kcubic,csp.kernelh,rr2,rhopp1,pressp1,rhopp2,pressp2));
}


//##############################################################################
//# Wendland kernel
//##############################################################################
//==============================================================================
/// Returns the name of the kernel in text format.
//==============================================================================
inline const char* GetKernelWendland_Name(){ return("Wendland"); }
//==============================================================================
/// Returns factor value to compute KernelSize according to KernelH.
//==============================================================================
inline float GetKernelWendland_Factor(){ return(2.0f); }
//==============================================================================
/// Returns constants for kernel calculation.
//==============================================================================
inline StKWendlandCte GetKernelWendland_Ctes(bool sim2d,double h){
  StKWendlandCte kc;
  if(sim2d){
    kc.awen=float(0.557/(h*h));
    kc.bwen=float(-2.7852/(h*h*h));
  }
  else{
    kc.awen=float(0.41778/(h*h*h));
    kc.bwen=float(-2.08891/(h*h*h*h));
  }
  return(kc);
}
//============================================================================== 
/// Returns wab of kernel.
//==============================================================================
inline float GetKernelWendland_Wab(const StKWendlandCte &kc,float h,float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  const float wqq=qq+qq+1.f;
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  return(kc.awen*wqq*wqq2*wqq2);
}
//============================================================================== 
/// Returns fac of kernel.
//==============================================================================
inline float GetKernelWendland_Fac(const StKWendlandCte &kc,float h,float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  const float wqq1=1.f-0.5f*qq;
  return(kc.bwen*qq*wqq1*wqq1*wqq1/rad);
}
//============================================================================== 
/// Returns wab and fac of kernel.
//==============================================================================
inline float GetKernelWendland_WabFac(const StKWendlandCte &kc,float h,float rr2,float &fac){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  fac=kc.bwen*qq*wqq2*wqq1/rad;
  const float wqq=qq+qq+1.f;
  return(kc.awen*wqq*wqq2*wqq2);
}

//============================================================================== 
/// Returns wab of kernel.
//==============================================================================
inline float GetKernelWendland_Wab(const StCteSph &csp,float rr2){
  return(GetKernelWendland_Wab(csp.kwend,csp.kernelh,rr2));
}
//============================================================================== 
/// Returns fac of kernel.
//==============================================================================
inline float GetKernelWendland_Fac(const StCteSph &csp,float rr2){
  return(GetKernelWendland_Fac(csp.kwend,csp.kernelh,rr2));
}
//============================================================================== 
/// Returns wab and fac of kernel.
//==============================================================================
inline float GetKernelWendland_WabFac(const StCteSph &csp,float rr2,float &fac){
  return(GetKernelWendland_WabFac(csp.kwend,csp.kernelh,rr2,fac));
}



//##############################################################################
//# Returns kernel information.
//##############################################################################
//============================================================================== 
/// Returns factor value to compute KernelSize according to KernelH.
//==============================================================================
inline float GetKernel_Factor(TpKernel tker){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Factor());
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_Factor());
  return(0);
}


//##############################################################################
//# Computes kernel values usig templates.
//##############################################################################
//============================================================================== 
/// Returns wab of kernel according to temaplate.
//==============================================================================
template<TpKernel tker> inline float GetKernel_Wab(const StCteSph &csp,float rr2){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Wab  (csp.kwend  ,csp.kernelh,rr2));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_Wab     (csp.kcubic ,csp.kernelh,rr2));
  else return(0);
}
//============================================================================== 
/// Returns fac of kernel  according to template.
//==============================================================================
template<TpKernel tker> inline float GetKernel_Fac(const StCteSph &csp,float rr2){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Fac  (csp.kwend  ,csp.kernelh,rr2));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_Fac     (csp.kcubic ,csp.kernelh,rr2));
  else return(0);
}
//============================================================================== 
/// Returns wab and fac of kernel according to template.
//==============================================================================
template<TpKernel tker> inline float GetKernel_WabFac(const StCteSph &csp,float rr2,float &fac){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_WabFac  (csp.kwend  ,csp.kernelh,rr2,fac));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_WabFac     (csp.kcubic ,csp.kernelh,rr2,fac));
  else return(0);
}


}

#endif


