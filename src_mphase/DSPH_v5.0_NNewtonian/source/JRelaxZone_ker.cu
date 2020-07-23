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

#include "JRelaxZone_ker.h"
#include "JReduSum_ker.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include <string>
#include <cstdio>
//:#include "JDgKerPrint.h"
//:#include "JDgKerPrint_ker.h"


namespace curelaxzone{
#include "FunctionsBasic_iker.h"
#include "FunctionsMath_iker.h"

//------------------------------------------------------------------------------
/// Returns TRUE when code==NULL or particle is normal and fluid.
//------------------------------------------------------------------------------
__device__ bool KerIsNormalFluid(const typecode *code,unsigned p){
  if(code){//-Descarta particulas floating o periodicas.
    const typecode cod=code[p];
    return(CODE_IsNormal(cod) && CODE_IsFluid(cod));
  }
  return(true);
}

//##############################################################################
//# Kernels for JRelaxZone uniform.
//##############################################################################

//------------------------------------------------------------------------------
/// Sets velocity of fluid to generate waves.
//------------------------------------------------------------------------------
__global__ void KerSetFluidVelUniform(unsigned n,unsigned pini
  ,float3 vt,float4 cenpla,float4 dompla1,float4 dompla2,float4 dompla3
  ,float domsize1,float domsize2,float domsize3,float widthhalf
  ,float coeff,double falpha,double fbeta,double fsub,double fdiv,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    const unsigned p1=p+pini;
    if(idp==NULL || idp[p1]>=fluidbeginidp){//-Ignore floating particles.
      const double2 pxy=posxy[p1];
      const double3 ps=make_double3(pxy.x,pxy.y,posz[p1]);
      const float pd1=float(cumath::PointPlane(dompla1,ps));
      const float pd2=float(cumath::PointPlane(dompla2,ps));
      const float pd3=float(cumath::PointPlane(dompla3,ps));
      if(pd1>0 && pd1<domsize1 && pd2>0 && pd2<domsize2 && pd3>0 && pd3<domsize3){
        const double pdis=cumath::DistPlaneSign(cenpla,ps);
        const double vdis=pdis/widthhalf;
        const double f1=(( tanh((vdis+falpha)*fbeta) - tanh((vdis-falpha)*fbeta) ) - fsub) / fdiv;
        const double f=f1*coeff;
        float4 v=velrhop[p1];
        v.x=float(f*vt.x+(1.-f)*v.x);
        v.y=float(f*vt.y+(1.-f)*v.y);
        v.z=float(f*vt.z+(1.-f)*v.z);
        velrhop[p1]=v;
      }
    }
  }
}

//==============================================================================
/// Sets velocity of fluid to generate waves.
//==============================================================================
void SetFluidVelUniform(unsigned n,unsigned pini
  ,const tfloat3 &vt,const tfloat4 &cenpla
  ,const tfloat4 &dompla1,const tfloat4 &dompla2,const tfloat4 &dompla3
  ,float domsize1,float domsize2,float domsize3,float widthhalf
  ,float coeff,double falpha,double fbeta,double fsub,double fdiv,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop)
{
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
    KerSetFluidVelUniform<<<sgrid,WAVEBSIZE>>> (n,pini,Float3(vt),Float4(cenpla)
      ,Float4(dompla1),Float4(dompla2),Float4(dompla3),domsize1,domsize2,domsize3
      ,widthhalf,coeff,falpha,fbeta,fsub,fdiv,fluidbeginidp,posxy,posz,idp,velrhop);
  }
}



//##############################################################################
//# Kernels for JRelaxZone regular.
//##############################################################################
//------------------------------------------------------------------------------
// Calcula velocidad del fluido en una posicion concreta.
// x es la distacia al punto de generacion y z la distancia al fondo.
//------------------------------------------------------------------------------
template <bool order2,bool subdrift> __device__ float KerCalcVelocityX(double x,double z
  ,double kl,double sinhkld,double wpf,double cta,double depth,double framp
  ,double ct2,double sinhkld4,double ctd,double ctd2)
{
  //const double ce=Gravity*WavePeriod/TWOPI*tanh(kl*Depth); //-1st and 2nd order wave celerity.
  //if(order2)vel.x+=3./4.*(PI*WaveHeight/WaveLength)*(PI*WaveHeight/WaveLength)*ce*cosh(2.*kl*(Depth+z))*(cos(4.*PI*f*time-2.*kl*x+2.*InitialPhase))/(sinhkld*sinhkld*sinhkld*sinhkld);
  double vx=(wpf*cosh(kl*(depth+z))*(cos(cta-kl*x))/sinhkld);
  if(order2)vx+=ct2*cosh(2.*kl*(depth+z))*(cos(2.*cta-2.*kl*x))/sinhkld4;
  if(subdrift)vx-=ctd*cosh(ctd2*(depth+z));
  return(float(vx*framp));
}

//------------------------------------------------------------------------------
// Calcula velocidad del fluido en una posicion concreta.
// x es la distacia al punto de generacion y z la distancia al fondo.
//------------------------------------------------------------------------------
template <bool order2> __device__ float KerCalcVelocityZ(double x,double z
  ,double kl,double sinhkld,double wpf,double cta,double depth,double framp
  ,double ct2,double sinhkld4)
{
  //const double ce=Gravity*WavePeriod/TWOPI*tanh(kl*Depth); //-1st and 2nd order wave celerity.
  //if(order2)vel.z-=3./4.*(PI*WaveHeight/WaveLength)*(PI*WaveHeight/WaveLength)*ce*sinh(2.*kl*(Depth+z))*(sin(4.*PI*f*time-2.*kl*x+2.*InitialPhase))/(sinhkld*sinhkld*sinhkld*sinhkld);
  double vz=(wpf*sinh(kl*(depth+z))*(sin(cta-kl*x))/sinhkld);
  if(order2)vz+=ct2*sinh(2.*kl*(depth+z))*(sin(2.*cta-2.*kl*x))/sinhkld4;
  return(-float(vz*framp));
}

//------------------------------------------------------------------------------
/// Sets velocity of fluid to generate waves.
//------------------------------------------------------------------------------
template <bool order2,bool subdrift> __global__ void KerSetFluidVel(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double kl,double sinhkld
  ,double wpf,double cta,double depth,double framp
  ,double ct2,double sinhkld4
  ,double ctd,double ctd2,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    const unsigned p1=p+pini;
    if(idp==NULL || idp[p1]>=fluidbeginidp){//-Ignore floating particles.
      const double px=posxy[p1].x-centerx;
      const double vdis=px/widthhalf;
      if(vdis>=-1 && vdis<=1){
        //const float f=((tanh((vdis+FunAlpha)*FunBeta)-tanh((vdis-FunAlpha)*FunBeta)) - (tanh((1+FunAlpha)*FunBeta)-tanh((1-FunAlpha)*FunBeta))) / ((tanh(FunAlpha*FunBeta)-tanh(-FunAlpha*FunBeta)) - (tanh((1+FunAlpha)*FunBeta) - tanh((1-FunAlpha)*FunBeta)));
        const double f=(( tanh((vdis+falpha)*fbeta) - tanh((vdis-falpha)*fbeta) ) - fsub) / fdiv;
        //const double f=cos(vdis*PIHALF);
        const double fx=f*coeffx;
        const double fz=f*coeffz;
        const double pzd=posz[p1]-swl;
        double vtx=KerCalcVelocityX<order2,subdrift>(px,pzd,kl,sinhkld,wpf,cta,depth,framp,ct2,sinhkld4,ctd,ctd2);
        double vtz=KerCalcVelocityZ<order2>(px,pzd,kl,sinhkld,wpf,cta,depth,framp,ct2,sinhkld4);
        float4 v=velrhop[p1];
        v.x=float(fx*vtx+(1.-fx)*v.x);
        v.z=float(fz*vtz+(1.-fz)*v.z);
        velrhop[p1]=v;
      }
    }
  }
}

//==============================================================================
/// Sets velocity of fluid to generate waves.
//==============================================================================
void SetFluidVel(unsigned n,unsigned pini,bool order2,bool subdrift
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double kl,double sinhkld
  ,double wpf,double cta,double depth,double framp
  ,double ct2,double sinhkld4
  ,double ctd,double ctd2,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop)
{
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
    if(!subdrift){ const bool sdrift=false;
      if(!order2)KerSetFluidVel<false,sdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,timewave,swl,kl,sinhkld,wpf,cta,depth,framp,ct2,sinhkld4,ctd,ctd2,fluidbeginidp,posxy,posz,idp,velrhop);
      else       KerSetFluidVel<true ,sdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,timewave,swl,kl,sinhkld,wpf,cta,depth,framp,ct2,sinhkld4,ctd,ctd2,fluidbeginidp,posxy,posz,idp,velrhop);
    }
    else{          const bool sdrift=true;
      if(!order2)KerSetFluidVel<false,sdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,timewave,swl,kl,sinhkld,wpf,cta,depth,framp,ct2,sinhkld4,ctd,ctd2,fluidbeginidp,posxy,posz,idp,velrhop);
      else       KerSetFluidVel<true ,sdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,timewave,swl,kl,sinhkld,wpf,cta,depth,framp,ct2,sinhkld4,ctd,ctd2,fluidbeginidp,posxy,posz,idp,velrhop);
    }
  }
}


//##############################################################################
//# Kernels for JRelaxZone irregular.
//##############################################################################
//------------------------------------------------------------------------------
// Calcula velocidad del fluido en una posicion concreta.
// x es la distacia al punto de generacion y z la distancia al fondo.
//------------------------------------------------------------------------------
__device__ double2 KerCalcVelxzSpectrum(double x,double z,double timewave
  ,double depth,unsigned wavecount
  ,const double *wavekl,const double *waveamp,const double *wavefang,const double *wavephase)
{
  double vx=0,vz=0;
  for(unsigned c=0;c<wavecount;c++){
    const double kl=wavekl[c];  //-Wave number.
    const double wf=wavefang[c];
    const double waf=waveamp[c]*wf;
    const double wft=wf*timewave-kl*x-wavephase[c];
    const double sinhkld=sinh(kl*depth);
    const double kldz=kl*(depth+z);
    vx+=waf*cosh(kldz)*(cos(wft))/sinhkld;
    vz-=waf*sinh(kldz)*(sin(wft))/sinhkld;
  }
  double2 vel=make_double2(vx,vz);
  return(vel);
}

//------------------------------------------------------------------------------
// Calcula velocidad del fluido en una posicion concreta.
// x es la distacia al punto de generacion y depth_z la distancia al fondo.
//------------------------------------------------------------------------------
template <byte subdriftmode> __device__ float KerCalcDriftX(double depth_z
  ,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2)
{
  float corr=0;
  if(subdriftmode==1)corr=float(ctd*cosh(ctd2*(depth_z)));
  if(subdriftmode==2){
    const double corr0=ctd  *cosh(ctd2  *(depth_z));
    const double corr1=ctd_2*cosh(ctd2_2*(depth_z));
    corr=float(corr1*fun+corr*(1.f-fun));
  }
  return(corr);
}

//------------------------------------------------------------------------------
/// Sets velocity of fluid to generate irregular waves.
//------------------------------------------------------------------------------
template <byte subdriftmode> __global__ void KerSetFluidVelSpectrumSub(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double depth,double framp,unsigned wavecount
  ,const double *wavekl,const double *waveamp,const double *wavefang,const double *wavephase  
  ,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop
  ,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    const unsigned p1=p+pini;
    if(idp==NULL || idp[p1]>=fluidbeginidp){//-Ignore floating particles.
      const double px=posxy[p1].x-centerx;
      const double vdis=px/widthhalf;
      if(vdis>=-1 && vdis<=1){
        //const float f=((tanh((vdis+FunAlpha)*FunBeta)-tanh((vdis-FunAlpha)*FunBeta)) - (tanh((1+FunAlpha)*FunBeta)-tanh((1-FunAlpha)*FunBeta))) / ((tanh(FunAlpha*FunBeta)-tanh(-FunAlpha*FunBeta)) - (tanh((1+FunAlpha)*FunBeta) - tanh((1-FunAlpha)*FunBeta)));
        const double f=(( tanh((vdis+falpha)*fbeta) - tanh((vdis-falpha)*fbeta) ) - fsub) / fdiv;
        //const double f=cos(vdis*PIHALF);
        const double fx=f*coeffx;
        const double fz=f*coeffz;
        const double pzd=posz[p1]-swl;
        double2 vel=KerCalcVelxzSpectrum(px,pzd,timewave,depth,wavecount,wavekl,waveamp,wavefang,wavephase);
        if(subdriftmode!=0)vel.x-=KerCalcDriftX<subdriftmode>(depth+pzd,fun,ctd,ctd2,ctd_2,ctd2_2);
        //-Applies the initial ramp to theoretical velocity.
        vel.x*=framp;
        vel.y*=framp;
        //-Upadates velocity of particles using function(alfa,beta).
        float4 v=velrhop[p1];
        v.x=float(fx*double(vel.x)+(1.-fx)*v.x);
        v.z=float(fz*double(vel.y)+(1.-fz)*v.z);
        velrhop[p1]=v;
      }
    }
  }
}

//==============================================================================
/// Sets velocity of fluid to generate irregular waves.
//==============================================================================
void SetFluidVelSpectrumSub(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double timewave,double swl,double depth,double framp,unsigned wavecount
  ,const double *wavekl,const double *waveamp,const double *wavefang,const double *wavephase  
  ,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop
  ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2)
{
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
    if(!subdrift)KerSetFluidVelSpectrumSub<0> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,timewave,swl,depth,framp,wavecount,wavekl,waveamp,wavefang,wavephase,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2);
    else{
      if(!fun)KerSetFluidVelSpectrumSub<1> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,timewave,swl,depth,framp,wavecount,wavekl,waveamp,wavefang,wavephase,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2);
      else    KerSetFluidVelSpectrumSub<2> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,timewave,swl,depth,framp,wavecount,wavekl,waveamp,wavefang,wavephase,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2);
    }
  }
}


//##############################################################################
//# Kernels for JRelaxZone using external velocity data (SWASH).
//##############################################################################
//------------------------------------------------------------------------------
// Devuelve datos de interpolacion en X, Y o Z.
//------------------------------------------------------------------------------
__device__ void KerGetInterpolationData(double pos,double posmin,double dp,unsigned np1,unsigned &cx,double &fx){
  double dx=(pos>posmin? pos-posmin: 0);
  cx=unsigned(dx/dp);
  if(cx>np1){ cx=np1; dx=0; } else dx-=dp*cx;
  fx=dx/dp;
}

//------------------------------------------------------------------------------
// Devuelve la velocidad en un punto por interpolacion de velocidades externas.
//------------------------------------------------------------------------------
template <bool usevelz> __device__ double2 KerCalcVelocityExternalXZ(double x,double y,double z
  ,double pxmin,double pymin,double pzmin
  ,double dpx,double dpy,double dpz
  ,unsigned npx1,unsigned npy1,unsigned npz1
  ,const double *velx,const double *velz)
{
  double2 vel=make_double2(0,0);
  //-Obtiene datos de interpolacion en X, Y, Z.
  unsigned cx=0,cy=0,cz=0;
  double fx=0,fy=0,fz=0;
  if(npx1)KerGetInterpolationData(x,pxmin,dpx,npx1,cx,fx);
  if(npy1)KerGetInterpolationData(y,pymin,dpy,npy1,cy,fy);
  if(npz1)KerGetInterpolationData(z,pzmin,dpz,npz1,cz,fz);
  const unsigned npx=npx1+1;
  const unsigned npz=npz1+1;
  const unsigned npxz=npx*npz;
  const unsigned p000=cy*npxz+cx*npz+cz;
  const unsigned p100=(cx<npx1? p000+npz: p000);
  const unsigned p001=(cz<npz1? p000+1:   p000);
  const unsigned p101=(cx<npx1? p001+npz: p001);
  const unsigned p010=(cy<npy1? p000+npxz: p000);
  const unsigned p110=(cy<npy1? p100+npxz: p100);
  const unsigned p011=(cy<npy1? p001+npxz: p001);
  const unsigned p111=(cy<npy1? p101+npxz: p101);
  //-Interpolacion de velocidad X.
  const double v000=fz*velx[p001]+(1.-fz)*velx[p000];
  const double v100=(npx1? fz*velx[p101]+(1.-fz)*velx[p100]: v000);
  const double v010=(npy1? fz*velx[p011]+(1.-fz)*velx[p010]: v000);
  const double v110=(npy1? fz*velx[p111]+(1.-fz)*velx[p110]: v100);
  const double v00=(npx1? fx*(v100-v000)+v000: v000);
  const double v01=(npx1? fx*(v110-v010)+v010: v010);
  vel.x=(npy1? fy*(v01-v00)+v00: v00);
  //-Interpolacion de velocidad Z.
  if(usevelz){
    const double v000=fz*velz[p001]+(1.-fz)*velz[p000];
    const double v100=(npx1? fz*velz[p101]+(1.-fz)*velz[p100]: v000);
    const double v010=(npy1? fz*velz[p011]+(1.-fz)*velz[p010]: v000);
    const double v110=(npy1? fz*velz[p111]+(1.-fz)*velz[p110]: v100);
    const double v00=(npx1? fx*(v100-v000)+v000: v000);
    const double v01=(npx1? fx*(v110-v010)+v010: v010);
    vel.y=(npy1? fy*(v01-v00)+v00: v00);
  }
  return(vel);
}

//------------------------------------------------------------------------------
/// Sets velocity of fluid to generate waves with external velocity data (SWASH).
//------------------------------------------------------------------------------
template<bool usevelz,byte subdriftmode> __global__ void KerSetFluidVelExternal(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double pxmin,double pymin,double pzmin
  ,double dpx,double dpy,double dpz
  ,unsigned npx1,unsigned npy1,unsigned npz1
  ,const double *velx,const double *velz
  ,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop
  ,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2,double bottom)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
    const unsigned p1=p+pini;
    if(idp==NULL || idp[p1]>=fluidbeginidp){//-Ignore floating particles.
      const double2 psxy=posxy[p1];
      const double px=psxy.x-centerx;
      const double vdis=px/widthhalf;
      if(vdis>=-1 && vdis<=1){
        //const float f=((tanh((vdis+FunAlpha)*FunBeta)-tanh((vdis-FunAlpha)*FunBeta)) - (tanh((1+FunAlpha)*FunBeta)-tanh((1-FunAlpha)*FunBeta))) / ((tanh(FunAlpha*FunBeta)-tanh(-FunAlpha*FunBeta)) - (tanh((1+FunAlpha)*FunBeta) - tanh((1-FunAlpha)*FunBeta)));
        const double f=(( tanh((vdis+falpha)*fbeta) - tanh((vdis-falpha)*fbeta) ) - fsub) / fdiv;
        //const double f=cos(vdis*PIHALF);
        const double fx=f*coeffx;
        const double fz=f*coeffz;
        const double psz=posz[p1];
        double2 vtxz=KerCalcVelocityExternalXZ<usevelz>(psxy.x,psxy.y,psz,pxmin,pymin,pzmin,dpx,dpy,dpz,npx1,npy1,npz1,velx,velz);
        if(subdriftmode!=0)vtxz.x-=KerCalcDriftX<subdriftmode>(psz-bottom,fun,ctd,ctd2,ctd_2,ctd2_2);
        float4 v=velrhop[p1];
        v.x=float(fx*vtxz.x+(1.-fx)*v.x);
        v.z=float(fz*vtxz.y+(1.-fz)*v.z);
        velrhop[p1]=v;
      }
    }
  }
}

//==============================================================================
/// Sets velocity of fluid to generate waves with external velocity data (SWASH).
//==============================================================================
void SetFluidVelExternal(unsigned n,unsigned pini
  ,double centerx,float widthhalf,float coeffx,float coeffz
  ,double falpha,double fbeta,double fsub,double fdiv
  ,double pxmin,double pymin,double pzmin
  ,double dpx,double dpy,double dpz
  ,unsigned npx1,unsigned npy1,unsigned npz1
  ,const double *velx,const double *velz
  ,unsigned fluidbeginidp
  ,const double2 *posxy,const double *posz,const unsigned *idp,float4 *velrhop
  ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2,double bottom)
{
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
    if(!subdrift){ const byte tdrift=0;
      if(!coeffz || !velz)KerSetFluidVelExternal<false,tdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,pxmin,pymin,pzmin,dpx,dpy,dpz,npx1,npy1,npz1,velx,velz,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2,bottom);
      else                KerSetFluidVelExternal<true ,tdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,pxmin,pymin,pzmin,dpx,dpy,dpz,npx1,npy1,npz1,velx,velz,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2,bottom);
    }
    else if(!fun){ const byte tdrift=1;
      if(!coeffz || !velz)KerSetFluidVelExternal<false,tdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,pxmin,pymin,pzmin,dpx,dpy,dpz,npx1,npy1,npz1,velx,velz,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2,bottom);
      else                KerSetFluidVelExternal<true ,tdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,pxmin,pymin,pzmin,dpx,dpy,dpz,npx1,npy1,npz1,velx,velz,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2,bottom);
    }
    else{          const byte tdrift=2;
      if(!coeffz || !velz)KerSetFluidVelExternal<false,tdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,pxmin,pymin,pzmin,dpx,dpy,dpz,npx1,npy1,npz1,velx,velz,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2,bottom);
      else                KerSetFluidVelExternal<true ,tdrift> <<<sgrid,WAVEBSIZE>>> (n,pini,centerx,widthhalf,coeffx,coeffz,falpha,fbeta,fsub,fdiv,pxmin,pymin,pzmin,dpx,dpy,dpz,npx1,npy1,npz1,velx,velz,fluidbeginidp,posxy,posz,idp,velrhop,fun,ctd,ctd2,ctd_2,ctd2_2,bottom);
    }
  }
}


//##############################################################################
//# Kernels for JRelaxZoneDrift.
//##############################################################################
//------------------------------------------------------------------------------
// Reduccion mediante suma de valores unsigned en memoria shared para un warp.
//------------------------------------------------------------------------------
template <unsigned blockSize> __device__ void KerReduSumUintWarp(volatile unsigned* sddat,unsigned tid){
  if(blockSize>=64)sddat[tid]+=sddat[tid+32];
  if(blockSize>=32)sddat[tid]+=sddat[tid+16];
  if(blockSize>=16)sddat[tid]+=sddat[tid+8];
  if(blockSize>=8) sddat[tid]+=sddat[tid+4];
  if(blockSize>=4) sddat[tid]+=sddat[tid+2];
  if(blockSize>=2) sddat[tid]+=sddat[tid+1];
}

//------------------------------------------------------------------------------
// Cuenta numero de particulas fluid entre xmin y xmax. El resultado se guarda
// en res[] por bloque.
//------------------------------------------------------------------------------
template <unsigned blockSize> __global__ void KerComputeDrift(unsigned n,unsigned pini
  ,double xmin,double xmax,const double2 *posxy,const typecode *code,unsigned *res)
{
  extern __shared__ unsigned shcount[];
  const unsigned tid=threadIdx.x;
  const unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  shcount[tid]=0;
  __syncthreads();
  if(p<n){
    const unsigned p1=p+pini;
    const double px=posxy[p1].x;
    if(xmin<=px && px<xmax){
      if(KerIsNormalFluid(code,p1))shcount[tid]=1;
    }
  }
  __syncthreads();
  if(blockSize>=512){ if(tid<256)shcount[tid]+=shcount[tid+256];  __syncthreads(); }
  if(blockSize>=256){ if(tid<128)shcount[tid]+=shcount[tid+128];  __syncthreads(); }
  if(blockSize>=128){ if(tid<64) shcount[tid]+=shcount[tid+64];   __syncthreads(); }
  if(tid<32)KerReduSumUintWarp<blockSize>(shcount,tid);
  if(tid==0)res[blockIdx.x]=shcount[0];
}

//==============================================================================
// Computes drift on CPU. Counts fluid particles between xmin and xmax.
//==============================================================================
unsigned ComputeDrift(unsigned n,unsigned pini,double xmin,double xmax
  ,const double2 *posxy,const typecode *code)
{
  unsigned ret=0;
  if(n){
    const dim3 sgrid=GetSimpleGridSize(n,WAVEBSIZE);
    const unsigned smemSize=sizeof(unsigned)*WAVEBSIZE;
    const unsigned nblocks=sgrid.x*sgrid.y;
    //-Allocates memory on GPU.
    unsigned *res=NULL;
    size_t size=sizeof(unsigned)*(nblocks+curedus::GetAuxSize_ReduSumUint(nblocks));
    cudaMalloc((void**)&res,size);
    //-Compute drift.
    KerComputeDrift<WAVEBSIZE><<<sgrid,WAVEBSIZE,smemSize>>>(n,pini,xmin,xmax,posxy,code,res);
    ret=curedus::ReduSumUint(nblocks,0,res,res+nblocks);
    //-Frees memory on GPU.
    cudaFree(res);
  }
  return(ret);
}


}


