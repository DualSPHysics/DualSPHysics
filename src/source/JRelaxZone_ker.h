//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2020 by Dr. Jose M. Dominguez
 All rights reserved.

 DualSPHysics is an international collaboration between:
 - EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 - School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
 the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the distribution.
 * Neither the name of the DualSPHysics nor the names of its contributors may be used to endorse or promote products derived 
   from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
 SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

#ifndef _JRelaxZone_ker_
#define _JRelaxZone_ker_

#include "DualSphDef.h"
#include <cuda_runtime_api.h>

#define WAVEBSIZE 256

namespace curelaxzone{

//# Kernels for JRelaxZone.
void SetFluidVelUniform(unsigned n,unsigned pini
  ,const tfloat3 &vt,const tfloat4 &cenpla
  ,const tfloat4 &dompla1,const tfloat4 &dompla2,const tfloat4 &dompla3
  ,float domsize1,float domsize2,float domsize3,float widthhalf
  ,float coeff,double falpha,double fbeta,double fsub,double fdiv,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop);

void SetFluidVel(unsigned n,unsigned pini,bool order2,bool subdrift
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double kl,double sinhkld
  ,double wpf,double cta,double depth,double framp
  ,double ct2,double sinhkld4
  ,double ctd,double ctd2,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop);

void SetFluidVelSpectrumSub(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double depth,double framp,unsigned wavecount
  ,const double *wavekl,const double *waveamp,const double *wavefang,const double *wavephase  
  ,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop
  ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2);

void SetFluidVelExternal(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double pxmin,double pymin,double pzmin
  ,double dpx,double dpy,double dpz
  ,unsigned npx1,unsigned npy1,unsigned npz1
  ,const double *velx,const double *velz
  ,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop
  ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2,double bottom);

//# Kernels for JRelaxZoneDrift.
unsigned ComputeDrift(unsigned n,unsigned pini,double xmin,double xmax
  ,const double2 *posxy,const typecode *code);


}

#endif


