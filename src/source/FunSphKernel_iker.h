//HEAD_DSCODES
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

/// \file FunSphKernel_iker.h \brief implements CUDA device functions for SPH kenernels.

#include "TypesDef.h"
#include "DualSphDef.h"
#include <cuda_runtime_api.h>

/// Implements CUDA device functions for SPH kenernels.
namespace cufsph{

//##############################################################################
//# Cubic Spline kernel
//##############################################################################
//#define CTE_AVAILABLE
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Returns the kernel (wab).
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_Wab(float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  if(rad>CTE.kernelh){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    return(CTE.cubic_a24*(wqq2*wqq1));
  }
  else{
    const float wqq2=qq*qq;
    return(CTE.cubic_a2*(1.0f+(0.75f*qq-1.5f)*wqq2));
  }
}
//------------------------------------------------------------------------------
/// Returns fac of kernel (module of the gradient of the kernel divided by the 
/// distance between particles).
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_Fac(float rr2){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  if(rad>CTE.kernelh){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    return(CTE.cubic_c2*wqq2/rad);
  }
  else{
    const float wqq2=qq*qq;
    return((CTE.cubic_c1*qq+CTE.cubic_d1*wqq2)/rad);
  }
}
//------------------------------------------------------------------------------
/// Returns wab (the kernel) and fac (module of the gradient of the kernel 
/// divided by the distance between particles).
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_WabFac(float rr2,float& fac){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  if(rad>CTE.kernelh){
    const float wqq1=2.0f-qq;
    const float wqq2=wqq1*wqq1;
    fac=(CTE.cubic_c2*wqq2/rad);
    return(CTE.cubic_a24*(wqq2*wqq1));
  }
  else{
    float wqq2=qq*qq;
    fac=((CTE.cubic_c1*qq+CTE.cubic_d1*wqq2)/rad);
    return(CTE.cubic_a2*(1.0f+(0.75f*qq-1.5f)*wqq2));
  }
}
//------------------------------------------------------------------------------
/// Return tensil correction of kernel Cubic.
//------------------------------------------------------------------------------
__device__ float GetKernelCubic_Tensil(float rr2
  ,float rhopp1,float pressp1,float rhopp2,float pressp2)
{
  const float wab=GetKernelCubic_Wab(rr2);
  float fab=wab*CTE.cubic_odwdeltap;
  fab*=fab; fab*=fab; //fab=fab^4
  const float tensilp1=(pressp1/(rhopp1*rhopp1))*(pressp1>0? 0.01f: -0.2f);
  const float tensilp2=(pressp2/(rhopp2*rhopp2))*(pressp2>0? 0.01f: -0.2f);
  return(fab*(tensilp1+tensilp2));
}
#endif

//##############################################################################
//# Wendland kernel
//##############################################################################
//------------------------------------------------------------------------------
/// Returns the kernel (wab).
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Wab(float rr2,float h,float awen){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  const float wqq=qq+qq+1.f;
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  return(awen*wqq*wqq2*wqq2);
}
//------------------------------------------------------------------------------
/// Returns fac of kernel (module of the gradient of the kernel divided by the 
/// distance between particles).
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Fac(float rr2,float h,float bwenh){
  const float rad=sqrt(rr2);
  const float qq=rad/h;
  const float wqq1=1.f-0.5f*qq;
  return(bwenh*wqq1*wqq1*wqq1);
}

#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Returns wab (the kernel) and fac (module of the gradient of the kernel 
/// divided by the distance between particles).
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_WabFac(float rr2,float& fac){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  fac=CTE.bwenh*wqq2*wqq1;
  const float wqq=qq+qq+1.f;
  return(CTE.awen*wqq*wqq2*wqq2);
}
//<vs_vrres_ini>
//------------------------------------------------------------------------------
/// Returns wab (the kernel) and fac (module of the gradient of the kernel 
/// divided by the distance between particles) and facc  (module of the second derivative of the kernel).
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_WabFacFacc(float rr2,float &fac,float &facc){
  const float rad=sqrt(rr2);
  const float qq=rad/CTE.kernelh;
  const float wqq1=1.f-0.5f*qq;
  const float wqq2=wqq1*wqq1;
  fac=CTE.bwenh*wqq2*wqq1;
  const float wqq=qq+qq+1.f;
  facc=-CTE.bwenh*wqq2*(qq+qq-1.f);
  return(CTE.awen*wqq*wqq2*wqq2);
}
//<vs_vrres_end>

//------------------------------------------------------------------------------
/// Returns the kernel (wab).
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Wab(float rr2){ 
  return(GetKernelWendland_Wab(rr2,CTE.kernelh,CTE.awen)); 
}
//------------------------------------------------------------------------------
/// Returns fac of kernel (module of the gradient of the kernel divided by the 
/// distance between particles).
//------------------------------------------------------------------------------
__device__ float GetKernelWendland_Fac(float rr2){
  return(GetKernelWendland_Fac(rr2,CTE.kernelh,CTE.bwenh));
}
#endif

//##############################################################################
//# Computes kernel values using templates.
//##############################################################################
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Returns the kernel (wab) according to temaplate.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Wab(float rr2){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Wab  (rr2));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_Wab     (rr2));
  else return(0);
}
//------------------------------------------------------------------------------
/// Returns fac of kernel (module of the gradient of the kernel divided by the 
/// distance between particles) according to template.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Fac(float rr2){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Fac  (rr2));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_Fac     (rr2));
  else return(0);
}
//------------------------------------------------------------------------------
/// Returns wab (the kernel) and fac (module of the gradient of the kernel 
/// divided by the distance between particles) according to template.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_WabFac(float rr2,float& fac){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_WabFac  (rr2,fac));
  else if(tker==KERNEL_Cubic     )return(GetKernelCubic_WabFac     (rr2,fac));
  else return(0);
}

//<vs_vrres_ini>
    //------------------------------------------------------------------------------
    /// Returns wab (the kernel) and fac (module of the gradient of the kernel 
    /// divided by the distance between particles) and facc  (module of the second derivative of the kernel).
    //------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_WabFacFacc(float rr2,float &fac,float &facc){
          if(tker==KERNEL_Wendland  )return(GetKernelWendland_WabFacFacc  (rr2,fac,facc));
      else return(0);
}
//<vs_vrres_end>

#endif
//------------------------------------------------------------------------------
/// Returns the kernel (wab) according to temaplate.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Wab(float rr2,float h,float aker){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Wab  (rr2,h,aker));
  else return(0);
}
//------------------------------------------------------------------------------
/// Returns fac of kernel (module of the gradient of the kernel divided by the 
/// distance between particles) according to template.
//------------------------------------------------------------------------------
template<TpKernel tker> __device__ float GetKernel_Fac(float rr2,float h,float bhker){
       if(tker==KERNEL_Wendland  )return(GetKernelWendland_Fac  (rr2,h,bhker));
  else return(0);
}

}

