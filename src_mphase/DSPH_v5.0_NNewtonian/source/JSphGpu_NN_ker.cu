//HEAD_DSPH
/*
<DUALSPHYSICS>  Copyright (c) 2019 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

This file is part of DualSPHysics.

DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.
*/

/// \file JSphGpu_ker.cu \brief Implements functions and CUDA kernels for the Particle Interaction and System Update.

#include "JSphGpu_NN_ker.h"
//#include <cfloat>
//#include <math_constants.h>

#define MAXNUMBERPHASE 10

__constant__ StPhaseCte PHASECTE[MAXNUMBERPHASE];
__constant__ StPhaseArray PHASEARRAY[MAXNUMBERPHASE];


namespace cusphNN {
#include "FunctionsBasic_iker.h" //orig
#include "FunctionsMath_iker.h" //orig
#include "FunSphKernel_iker.h"
#include "FunSphEos_iker.h"
#undef _JCellSearch_iker_
#include "JCellSearch_iker.h"


//==============================================================================
/// Stores constants for the GPU interaction.
/// Graba constantes para la interaccion a la GPU.
//==============================================================================
void CteInteractionUp_NN(unsigned phasecount,const StPhaseCte *phasecte,const StPhaseArray *phasearray) {
  cudaMemcpyToSymbol(PHASECTE,phasecte,sizeof(StPhaseCte)*phasecount);
  cudaMemcpyToSymbol(PHASEARRAY,phasearray,sizeof(StPhaseArray)*phasecount);
}

//------------------------------------------------------------------------------
/// Doubles the position of the indicated particle using a displacement.
/// Duplicate particles are considered valid and are always within
/// the domain.
/// This kernel applies to single-GPU and multi-GPU because the calculations are made
/// from domposmin.
/// It controls the cell coordinates not exceed the maximum.
///
/// Duplica la posicion de la particula indicada aplicandole un desplazamiento.
/// Las particulas duplicadas se considera que siempre son validas y estan dentro
/// del dominio.
/// Este kernel vale para single-gpu y multi-gpu porque los calculos se hacen 
/// a partir de domposmin.
/// Se controla que las coordendas de celda no sobrepasen el maximo.
//------------------------------------------------------------------------------
__device__ void KerPeriodicDuplicatePos(unsigned pnew,unsigned pcopy
  ,bool inverse,double dx,double dy,double dz,uint3 cellmax
  ,double2 *posxy,double *posz,unsigned *dcell)
{
  //-Obtains position of the particle to be duplicated.
  //-Obtiene pos de particula a duplicar.
  double2 rxy=posxy[pcopy];
  double rz=posz[pcopy];
  //-Applies displacement.
  rxy.x+=(inverse ? -dx : dx);
  rxy.y+=(inverse ? -dy : dy);
  rz+=(inverse ? -dz : dz);
  //-Computes cell coordinates within the domain.
  //-Calcula coordendas de celda dentro de dominio.
  unsigned cx=unsigned((rxy.x-CTE.domposminx)/CTE.scell);
  unsigned cy=unsigned((rxy.y-CTE.domposminy)/CTE.scell);
  unsigned cz=unsigned((rz-CTE.domposminz)/CTE.scell);
  //-Adjust cell coordinates if they exceed the maximum.
  //-Ajusta las coordendas de celda si sobrepasan el maximo.
  cx=(cx<=cellmax.x ? cx : cellmax.x);
  cy=(cy<=cellmax.y ? cy : cellmax.y);
  cz=(cz<=cellmax.z ? cz : cellmax.z);
  //-Stores position and cell of the new particles.
  //-Graba posicion y celda de nuevas particulas.
  posxy[pnew]=rxy;
  posz[pnew]=rz;
  dcell[pnew]=PC__Cell(CTE.cellcode,cx,cy,cz);
}
//------------------------------------------------------------------------------
/// Creates periodic particles from a list of particles to duplicate for non-Newtonian models.
/// It is assumed that all particles are valid.
/// This kernel applies to single-GPU and multi-GPU because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-gpu y multi-gpu porque usa domposmin. 
//------------------------------------------------------------------------------
template<bool varspre> __global__ void KerPeriodicDuplicateSymplectic_NN(unsigned n,unsigned pini
  ,uint3 cellmax,double3 perinc,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,float *auxnn,double2 *posxypre,double *poszpre,float4 *velrhoppre)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Adjusts cell position of the new particles.
    //-Ajusta posicion y celda de nueva particula.
    KerPeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000),perinc.x,perinc.y,perinc.z,cellmax,posxy,posz,dcell);
    //-Copies the remaining data.
    //-Copia el resto de datos.
    idp[pnew]=idp[pcopy];
    code[pnew]=CODE_SetPeriodic(code[pcopy]);
    velrhop[pnew]=velrhop[pcopy];
    if(varspre) {
      posxypre[pnew]=posxypre[pcopy];
      poszpre[pnew]=poszpre[pcopy];
      velrhoppre[pnew]=velrhoppre[pcopy];
    }
    if(auxnn)auxnn[pnew]=auxnn[pcopy];
  }
}

//==============================================================================
/// Creates periodic particles from a list of particles to duplicate for non-Newotnian formulation..
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
//==============================================================================
void PeriodicDuplicateSymplectic(unsigned n,unsigned pini
  ,tuint3 domcells,tdouble3 perinc,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,float *auxnn,double2 *posxypre,double *poszpre,float4 *velrhoppre)
{
  if(n) {
    uint3 cellmax=make_uint3(domcells.x-1,domcells.y-1,domcells.z-1);
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    if(posxypre!=NULL)KerPeriodicDuplicateSymplectic_NN<true><<<sgrid,SPHBSIZE>>>(n,pini,cellmax,Double3(perinc),listp,idp,code,dcell,posxy,posz,velrhop,auxnn,posxypre,poszpre,velrhoppre);
    else                 KerPeriodicDuplicateSymplectic_NN<false><<<sgrid,SPHBSIZE>>>(n,pini,cellmax,Double3(perinc),listp,idp,code,dcell,posxy,posz,velrhop,auxnn,posxypre,poszpre,velrhoppre);
  }
}

//------------------------------------------------------------------------------
/// Creates periodic particles from a list of particles to duplicate.
/// It is assumed that all particles are valid.
/// This kernel applies to single-GPU and multi-GPU because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-gpu y multi-gpu porque usa domposmin. 
//------------------------------------------------------------------------------
__global__ void KerPeriodicDuplicateVerlet(unsigned n,unsigned pini,uint3 cellmax,double3 perinc
  ,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,float *auxnn,float4 *velrhopm1)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    const unsigned pnew=p+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Adjusts cell position of the new particles.
    //-Ajusta posicion y celda de nueva particula.
    KerPeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000),perinc.x,perinc.y,perinc.z,cellmax,posxy,posz,dcell);
    //-Copies the remaining data.
    //-Copia el resto de datos.
    idp[pnew]=idp[pcopy];
    code[pnew]=CODE_SetPeriodic(code[pcopy]);
    velrhop[pnew]=velrhop[pcopy];
    velrhopm1[pnew]=velrhopm1[pcopy];
    if(auxnn)auxnn[pnew]=auxnn[pcopy];
  }
}

//==============================================================================
/// Creates periodic particles from a list of particles to duplicate.
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
//==============================================================================
void PeriodicDuplicateVerlet(unsigned n,unsigned pini,tuint3 domcells,tdouble3 perinc
  ,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,float *auxnn,float4 *velrhopm1)
{
  if(n) {
    uint3 cellmax=make_uint3(domcells.x-1,domcells.y-1,domcells.z-1);
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerPeriodicDuplicateVerlet<<<sgrid,SPHBSIZE>>>(n,pini,cellmax,Double3(perinc),listp,idp,code,dcell,posxy,posz,velrhop,auxnn,velrhopm1);
  }
}

//##############################################################################
//# Kernels for calculating NN Tensors.
//# Kernels para calcular tensores.
//##############################################################################
//------------------------------------------------------------------------------
/// Velocity gradients for non-Newtonian models using FDAs approach.
/// Gradientes de velocidad usando FDAs.
//------------------------------------------------------------------------------
__device__ void GetVelocityGradients_FDA(float rr2,float drx,float dry,float drz
  ,float dvx,float dvy,float dvz,tmatrix3f &dvelp1,float &div_vel)
{
  //vel gradients
  dvelp1.a11=dvx*drx/rr2; dvelp1.a12=dvx*dry/rr2; dvelp1.a13=dvx*drz/rr2; //Fan et al., 2010
  dvelp1.a21=dvy*drx/rr2; dvelp1.a22=dvy*dry/rr2; dvelp1.a23=dvy*drz/rr2;
  dvelp1.a31=dvz*drx/rr2; dvelp1.a32=dvz*dry/rr2; dvelp1.a33=dvz*drz/rr2;
  div_vel=(dvelp1.a11+dvelp1.a22+dvelp1.a33)/3.f;
}

//==============================================================================
//symetric tensors
//==============================================================================
/// Calculates the Stress Tensor (symetric)
/// Obtener tensor de velocidad de deformacion symetric.
//==============================================================================
__device__ void GetStressTensor_sym(float2 &d_p1_xx_xy,float2 &d_p1_xz_yy,float2 &d_p1_yz_zz,float visco_etap1
  ,float &I_t,float &II_t,float &J1_t,float &J2_t,float &tau_tensor_magn
  ,float2 &tau_xx_xy,float2 &tau_xz_yy,float2 &tau_yz_zz)
{
  //Stress tensor and invariant
  tau_xx_xy.x=2.f*visco_etap1*(d_p1_xx_xy.x);	tau_xx_xy.y=2.f*visco_etap1*d_p1_xx_xy.y;		tau_xz_yy.x=2.f*visco_etap1*d_p1_xz_yy.x;
  tau_xz_yy.y=2.f*visco_etap1*(d_p1_xz_yy.y);	tau_yz_zz.x=2.f*visco_etap1*d_p1_yz_zz.x;
  tau_yz_zz.y=2.f*visco_etap1*(d_p1_yz_zz.y);
  //I_t - the first invariant -
  I_t=tau_xx_xy.x+tau_xz_yy.y+tau_yz_zz.y;
  //II_t - the second invariant - expnaded form witout symetry 
  float II_t_1=tau_xx_xy.x*tau_xz_yy.y+tau_xz_yy.y*tau_yz_zz.y+tau_xx_xy.x*tau_yz_zz.y;
  float II_t_2=tau_xx_xy.y*tau_xx_xy.y+tau_yz_zz.x*tau_yz_zz.x+tau_xz_yy.x*tau_xz_yy.x;
  II_t=-II_t_1+II_t_2;
  //stress tensor magnitude
  tau_tensor_magn=sqrt(II_t);
  //if (II_t < 0.f) {
  //	printf("****tau_tensor_magn is negative**** \n");
  //}
  //Main Stress rate invariants
  J1_t=I_t; J2_t=I_t*I_t-2.f*II_t;
}

//==============================================================================
/// Calculates the Strain Rate Tensor (symetric).
/// Obtener tensor de velocidad de deformacion symetric.
//==============================================================================
__device__ void GetStrainRateTensor_tsym(float2 &dvelp1_xx_xy,float2 &dvelp1_xz_yy,float2 &dvelp1_yz_zz
  ,float &I_D,float &II_D,float &J1_D,float &J2_D,float &div_D_tensor,float &D_tensor_magn
  ,float2 &D_tensor_xx_xy,float2 &D_tensor_xz_yy,float2 &D_tensor_yz_zz)
{
  //Strain tensor and invariant	
  float div_vel=(dvelp1_xx_xy.x+dvelp1_xz_yy.y+dvelp1_yz_zz.y)/3.f;
  D_tensor_xx_xy.x=dvelp1_xx_xy.x-div_vel;		D_tensor_xx_xy.y=0.5f*(dvelp1_xx_xy.y);		D_tensor_xz_yy.x=0.5f*(dvelp1_xz_yy.x);
  D_tensor_xz_yy.y=dvelp1_xz_yy.y-div_vel;	D_tensor_yz_zz.x=0.5f*(dvelp1_yz_zz.x);
  D_tensor_yz_zz.y=dvelp1_yz_zz.y-div_vel;
  //the off-diagonal entries of velocity gradients are i.e. 0.5f*(du/dy+dvdx) with dvelp1.xy=du/dy+dvdx
  div_D_tensor=(D_tensor_xx_xy.x+D_tensor_xz_yy.y+D_tensor_yz_zz.y)/3.f;

  ////I_D - the first invariant -
  I_D=D_tensor_xx_xy.x+D_tensor_xz_yy.y+D_tensor_yz_zz.y;
  //II_D - the second invariant - expnaded form witout symetry 
  float II_D_1=D_tensor_xx_xy.x*D_tensor_xz_yy.y+D_tensor_xz_yy.y*D_tensor_yz_zz.y+D_tensor_xx_xy.x*D_tensor_yz_zz.y;
  float II_D_2=D_tensor_xx_xy.y*D_tensor_xx_xy.y+D_tensor_yz_zz.x*D_tensor_yz_zz.x+D_tensor_xz_yy.x*D_tensor_xz_yy.x;
  II_D=-II_D_1+II_D_2;
  ////deformation tensor magnitude
  D_tensor_magn=sqrt((II_D));

  //Main Strain rate invariants
  J1_D=I_D; J2_D=I_D*I_D-2.f*II_D;
}

//==============================================================================
/// Velocity gradients using SPH approach.
/// Gradientes de velocidad usando SPH.
//==============================================================================
__device__ void GetVelocityGradients_SPH_tsym(float massp2,const float4 &velrhop2,float dvx,float dvy,float dvz,float frx,float fry,float frz
  ,float2 &grap1_xx_xy,float2 &grap1_xz_yy,float2 &grap1_yz_zz)
{
  ///SPH vel gradients calculation
  const float volp2=-massp2/velrhop2.w;
  float dv=dvx*volp2;  grap1_xx_xy.x+=dv*frx; grap1_xx_xy.y+=dv*fry; grap1_xz_yy.x+=dv*frz;
  dv=dvy*volp2;  grap1_xx_xy.y+=dv*frx;	grap1_xz_yy.y+=dv*fry; grap1_yz_zz.x+=dv*frz;
  dv=dvz*volp2;  grap1_xz_yy.x+=dv*frx; grap1_yz_zz.x+=dv*fry; grap1_yz_zz.y+=dv*frz;
}

//==============================================================================
/// Calculate strain rate tensor (full matrix).
/// Obtener tensor de velocidad de deformacion (full matrix).
//==============================================================================
__device__ void GetStrainRateTensor(const tmatrix3f &dvelp1,float div_vel,float &I_D,float &II_D,float &J1_D
  ,float &J2_D,float &div_D_tensor,float &D_tensor_magn,tmatrix3f &D_tensor)
{
  //Strain tensor and invariant
  D_tensor.a11=dvelp1.a11-div_vel;          D_tensor.a12=0.5f*(dvelp1.a12+dvelp1.a21);      D_tensor.a13=0.5f*(dvelp1.a13+dvelp1.a31);
  D_tensor.a21=0.5f*(dvelp1.a21+dvelp1.a12);      D_tensor.a22=dvelp1.a22-div_vel;          D_tensor.a23=0.5f*(dvelp1.a23+dvelp1.a32);
  D_tensor.a31=0.5f*(dvelp1.a31+dvelp1.a13);      D_tensor.a32=0.5f*(dvelp1.a32+dvelp1.a23);      D_tensor.a33=dvelp1.a33-div_vel;
  div_D_tensor=(D_tensor.a11+D_tensor.a22+D_tensor.a33)/3.f;

  //I_D - the first invariant -
  I_D=D_tensor.a11+D_tensor.a22+D_tensor.a33;
  //II_D - the second invariant - expnaded form witout symetry 
  float II_D_1=D_tensor.a11*D_tensor.a22+D_tensor.a22*D_tensor.a33+D_tensor.a11*D_tensor.a33;
  float II_D_2=D_tensor.a12*D_tensor.a21+D_tensor.a23*D_tensor.a32+D_tensor.a13*D_tensor.a31;
  II_D=II_D_1-II_D_2;
  //deformation tensor magnitude
  D_tensor_magn=sqrt((II_D*II_D));

  //Main Strain rate invariants
  J1_D=I_D; J2_D=I_D*I_D-2.f*II_D;
}
//==============================================================================
/// Calculates the effective visocity.
/// Calcule la viscosidad efectiva.
//==============================================================================
__device__ void KerGetEta_Effective(const typecode ppx,float tau_yield,float D_tensor_magn,float visco
  ,float m_NN,float n_NN,float &visco_etap1)
{

  if(D_tensor_magn<=ALMOSTZERO)D_tensor_magn=ALMOSTZERO;
  float miou_yield=(PHASECTE[ppx].tau_max ? PHASECTE[ppx].tau_max/(2.0f*D_tensor_magn) : (tau_yield)/(2.0f*D_tensor_magn)); //HPB will adjust eta		

  //if tau_max exists
  bool bi_region=PHASECTE[ppx].tau_max && D_tensor_magn<=PHASECTE[ppx].tau_max/(2.f*PHASECTE[ppx].Bi_multi*visco);
  if(bi_region) { //multiplier
    miou_yield=PHASECTE[ppx].Bi_multi*visco;
  }
  //Papanastasiou
  float miouPap=miou_yield *(1.f-exp(-m_NN*D_tensor_magn));
  float visco_etap1_term1=(PHASECTE[ppx].tau_max ? miou_yield : (miouPap>m_NN*tau_yield||D_tensor_magn==ALMOSTZERO ? m_NN*tau_yield : miouPap));

  //HB
  float miouHB=visco*pow(D_tensor_magn,(n_NN-1.0f));
  float visco_etap1_term2=(bi_region ? visco : (miouPap>m_NN*tau_yield||D_tensor_magn==ALMOSTZERO ? visco : miouHB));

  visco_etap1=visco_etap1_term1+visco_etap1_term2;

  /*
  //use according to you criteria
  - Herein we limit visco_etap1 at very low shear rates
  */
}
//------------------------------------------------------------------------------
/// Calclulate stress tensor.
/// Calcular tensor de estres.
//------------------------------------------------------------------------------
__device__ void GetStressTensor(const tmatrix3f &D_tensor,float visco_etap1,float &I_t,float &II_t,float &J1_t
  ,float &J2_t,float &tau_tensor_magn,tmatrix3f &tau_tensor)
{
  //Stress tensor and invariant
  tau_tensor.a11=2.f*visco_etap1*(D_tensor.a11);	tau_tensor.a12=2.f*visco_etap1*D_tensor.a12;		tau_tensor.a13=2.f*visco_etap1*D_tensor.a13;
  tau_tensor.a21=2.f*visco_etap1*D_tensor.a21;		tau_tensor.a22=2.f*visco_etap1*(D_tensor.a22);	tau_tensor.a23=2.f*visco_etap1*D_tensor.a23;
  tau_tensor.a31=2.f*visco_etap1*D_tensor.a31;		tau_tensor.a32=2.f*visco_etap1*D_tensor.a32;		tau_tensor.a33=2.f*visco_etap1*(D_tensor.a33);

  //I_t - the first invariant -
  I_t=tau_tensor.a11+tau_tensor.a22+tau_tensor.a33;
  //II_t - the second invariant - expnaded form witout symetry 
  float II_t_1=tau_tensor.a11*tau_tensor.a22+tau_tensor.a22*tau_tensor.a33+tau_tensor.a11*tau_tensor.a33;
  float II_t_2=tau_tensor.a12*tau_tensor.a21+tau_tensor.a23*tau_tensor.a32+tau_tensor.a13*tau_tensor.a31;
  II_t=II_t_1-II_t_2;
  //stress tensor magnitude
  tau_tensor_magn=sqrt(II_t*II_t);

  //Main Strain rate invariants
  J1_t=I_t; J2_t=I_t*I_t-2.f*II_t;
}

//##############################################################################
//# Kernels for calculating forces (Pos-Double) for non-Newtonian models.
//# Kernels para calculo de fuerzas (Pos-Double) para modelos no-Newtonianos.
//##############################################################################
//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles. Bound-Fluid/Float
/// Realiza la interaccion de una particula con un conjunto de ellas. Bound-Fluid/Float
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,bool symm>
__device__ void KerInteractionForcesBoundBox_NN
(unsigned p1,const unsigned &pini,const unsigned &pfin
  ,const float *ftomassp
  ,const float4 *poscell,const float4 *velrhop,const typecode *code,const unsigned* idp
  ,float massf,const float4 &pscellp1,const float4 &velrhop1,float &arp1,float &visc)
{
  for(int p2=pini; p2<pfin; p2++) {
    const float4 pscellp2=poscell[p2];
    float drx=pscellp1.x-pscellp2.x+CTE.poscellsize*(CEL_GetX(__float_as_int(pscellp1.w))-CEL_GetX(__float_as_int(pscellp2.w)));
    float dry=pscellp1.y-pscellp2.y+CTE.poscellsize*(CEL_GetY(__float_as_int(pscellp1.w))-CEL_GetY(__float_as_int(pscellp2.w)));
    float drz=pscellp1.z-pscellp2.z+CTE.poscellsize*(CEL_GetZ(__float_as_int(pscellp1.w))-CEL_GetZ(__float_as_int(pscellp2.w)));
    if(symm)dry=pscellp1.y+pscellp2.y+CTE.poscellsize*CEL_GetY(__float_as_int(pscellp2.w)); //<vs_syymmetry>
    const float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO) {
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

      float4 velrhop2=velrhop[p2];
      if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

      //-Obtains particle mass p2 if there are floating bodies.
      //-Obtiene masa de particula p2 en caso de existir floatings.
      float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massf si es fluid.
      bool compute=true; //-Deactivated when DEM is used and is float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
      if(USE_FLOATING) {
        const typecode cod=code[p2];
        bool ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2 ? ftomassp[CODE_GetTypeValue(cod)] : massf);
        compute=!(USE_FTEXTERNAL && ftp2); //-Deactivated when DEM or Chrono is used and is bound-float. | Se desactiva cuando se usa DEM o Chrono y es bound-float.
      }

      if(compute) {
        //-Density derivative (Continuity equation).
        const float dvx=velrhop1.x-velrhop2.x,dvy=velrhop1.y-velrhop2.y,dvz=velrhop1.z-velrhop2.z;
        arp1+=(USE_FLOATING ? ftmassp2 : massf)*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

        {//===== Viscosity ===== 
          const float dot=drx*dvx+dry*dvy+drz*dvz;
          const float dot_rr2=dot/(rr2+CTE.eta2);
          visc=max(dot_rr2,visc);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Particle interaction for non-Newtonian models. Bound-Fluid/Float 
/// Realiza interaccion entre particulas para modelos no-Newtonianos. Bound-Fluid/Float
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,bool symm>
__global__ void KerInteractionForcesBound_NN(unsigned n,unsigned pinit
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *beginendcellfluid,const unsigned *dcell
  ,const float *ftomassp
  ,const float4 *poscell,const float4 *velrhop,const typecode *code,const unsigned *idp
  ,float *viscdt,float *ar)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of thread.
  if(p<n) {
    const unsigned p1=p+pinit;      //-Number of particle.
    float visc=0,arp1=0;

    //-Loads particle p1 data.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrhop[p1];
    const bool rsymp1=(symm && CEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Boundary-Fluid interaction.
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x) {
      unsigned pini,pfin=0; cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
      if(pfin) {
        KerInteractionForcesBoundBox_NN<tker,ftmode,false>(p1,pini,pfin,ftomassp,poscell,velrhop,code,idp,CTE.massf,pscellp1,velrhop1,arp1,visc);
        if(symm && rsymp1)KerInteractionForcesBoundBox_NN<tker,ftmode,true >(p1,pini,pfin,ftomassp,poscell,velrhop,code,idp,CTE.massf,pscellp1,velrhop1,arp1,visc);
      }
    }
    //-Stores results.
    if(arp1||visc) {
      ar[p1]+=arp1;
      if(visc>viscdt[p1])viscdt[p1]=visc;
    }
  }
}
//======================Start of FDA approach===================================
//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles for non-Newtonian models using the FDA approach. (Fluid/Float-Fluid/Float/Bound)
/// Realiza la interaccion de una particula con un conjunto de ellas para modelos no Newtonianos que utilizan el enfoque de la FDA. (Fluid/Float-Fluid/Float/Bound)
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool symm>
__device__ void KerInteractionForcesFluidBox_FDA(bool boundp2,unsigned p1
  ,const unsigned &pini,const unsigned &pfin,float visco,float *visco_eta
  ,const float *ftomassp,float2 *tauff
  ,const float4 *poscell,const float4 *velrhop,const typecode *code,const unsigned *idp
  ,float massp2,const typecode pp1,bool ftp1
  ,const float4 &pscellp1,const float4 &velrhop1,float pressp1
  ,float2 &taup1_xx_xy,float2 &taup1_xz_yy,float2 &taup1_yz_zz
  ,float2 &grap1_xx_xy,float2 &grap1_xz_yy,float2 &grap1_yz_zz
  ,float3 &acep1,float &arp1,float &visc,float &visceta,float &visco_etap1,float &deltap1
  ,TpShifting shiftmode,float4 &shiftposfsp1)
{
  for(int p2=pini; p2<pfin; p2++){
    const float4 pscellp2=poscell[p2];
    float drx=pscellp1.x-pscellp2.x+CTE.poscellsize*(CEL_GetX(__float_as_int(pscellp1.w))-CEL_GetX(__float_as_int(pscellp2.w)));
    float dry=pscellp1.y-pscellp2.y+CTE.poscellsize*(CEL_GetY(__float_as_int(pscellp1.w))-CEL_GetY(__float_as_int(pscellp2.w)));
    float drz=pscellp1.z-pscellp2.z+CTE.poscellsize*(CEL_GetZ(__float_as_int(pscellp1.w))-CEL_GetZ(__float_as_int(pscellp2.w)));
    if(symm)dry=pscellp1.y+pscellp2.y+CTE.poscellsize*CEL_GetY(__float_as_int(pscellp2.w)); //<vs_syymmetry>
    const float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

      //-Obtains mass of particle p2 for NN and if any floating bodies exist.
      const typecode cod=code[p2];
      const typecode pp2=(boundp2 ? pp1 : CODE_GetTypeValue(cod)); //<vs_non-Newtonian>
      float massp2=(boundp2 ? CTE.massb : PHASEARRAY[pp2].mass); //massp2 not neccesary to go in _Box function
      //Note if you masses are very different more than a ratio of 1.3 then: massp2 = (boundp2 ? PHASEARRAY[pp1].mass : PHASEARRAY[pp2].mass);

      //-Obtiene masa de particula p2 en caso de existir floatings.
      bool ftp2=false;         //-Indicates if it is floating. | Indica si es floating.
      float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
      bool compute=true; //-Deactivated when DEM is used and is float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
      if(USE_FLOATING) {
        ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2 ? ftomassp[CODE_GetTypeValue(cod)] : massp2);
#ifdef DELTA_HEAVYFLOATING
        if(ftp2 && tdensity==DDT_DDT && ftmassp2<=(massp2*1.2f))deltap1=FLT_MAX;
#else
        if(ftp2 && tdensity==DDT_DDT)deltap1=FLT_MAX;
#endif
        if(ftp2 && shift && shiftmode==SHIFT_NoBound)shiftposfsp1.x=FLT_MAX; //-Cancels shifting with floating bodies. | Con floatings anula shifting.
        compute=!(USE_FTEXTERNAL && ftp1&&(boundp2||ftp2)); //-Deactivated when DEM or Chrono is used and is float-float or float-bound. | Se desactiva cuando se usa DEM o Chrono y es float-float o float-bound.
      }

      float4 velrhop2=velrhop[p2];
      if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

      //===== Aceleration ===== 
      if(compute) {
        const float pressp2=cufsph::ComputePressCte_NN(velrhop2.w,PHASEARRAY[pp2].rho,PHASEARRAY[pp2].CteB,PHASEARRAY[pp2].Gamma,PHASEARRAY[pp2].Cs0,cod);
        const float prs=(pressp1+pressp2)/(velrhop1.w*velrhop2.w)
          +(tker==KERNEL_Cubic ? cufsph::GetKernelCubic_Tensil(rr2,velrhop1.w,pressp1,velrhop2.w,pressp2) : 0);
        const float p_vpm=-prs*(USE_FLOATING ? ftmassp2 : massp2);
        acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
      }

      //-Density derivative (Continuity equation).
      float dvx=velrhop1.x-velrhop2.x,dvy=velrhop1.y-velrhop2.y,dvz=velrhop1.z-velrhop2.z;
      if(compute)arp1+=(USE_FLOATING ? ftmassp2 : massp2)*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

      const float cbar=max(PHASEARRAY[pp1].Cs0,PHASEARRAY[pp2].Cs0);
      const float dot3=(tdensity!=DDT_None||shift ? drx*frx+dry*fry+drz*frz : 0);
      //-Density derivative (DeltaSPH Molteni).
      if(tdensity==DDT_DDT && deltap1!=FLT_MAX) {
        const float rhop1over2=velrhop1.w/velrhop2.w;
        const float visc_densi=CTE.ddtkh*cbar*(rhop1over2-1.f)/(rr2+CTE.eta2);
        const float delta=(pp1==pp2 ? visc_densi*dot3*(USE_FLOATING ? ftmassp2 : massp2) : 0); //<vs_non-Newtonian>
        //deltap1=(boundp2? FLT_MAX: deltap1+delta);
        deltap1=(boundp2 && CTE.tboundary==BC_DBC ? FLT_MAX : deltap1+delta);
      }
      //-Density Diffusion Term (Fourtakas et al 2019). //<vs_dtt2_ini>
      if((tdensity==DDT_DDT2||(tdensity==DDT_DDT2Full&&!boundp2))&&deltap1!=FLT_MAX&&!ftp2) {
        const float rh=1.f+CTE.ddtgz*drz;
        const float drhop=CTE.rhopzero*pow(rh,1.f/CTE.gamma)-CTE.rhopzero;
        const float visc_densi=CTE.ddtkh*cbar*((velrhop2.w-velrhop1.w)-drhop)/(rr2+CTE.eta2);
        const float delta=(pp1==pp2 ? visc_densi*dot3*massp2/velrhop2.w : 0); //<vs_non-Newtonian>
        deltap1=(boundp2 ? FLT_MAX : deltap1-delta);
      } //<vs_dtt2_end>		

      //-Shifting correction.
      if(shift && shiftposfsp1.x!=FLT_MAX) {
        bool heavyphase=(PHASEARRAY[pp1].mass>PHASEARRAY[pp2].mass && pp1!=pp2 ? true : false); //<vs_non-Newtonian>
        const float massrhop=(USE_FLOATING ? ftmassp2 : massp2)/velrhop2.w;
        const bool noshift=(boundp2&&(shiftmode==SHIFT_NoBound||(shiftmode==SHIFT_NoFixed && CODE_IsFixed(code[p2]))));
        shiftposfsp1.x=(noshift ? FLT_MAX : (heavyphase ? 0 : shiftposfsp1.x+massrhop*frx)); //-Removes shifting for the boundaries. | Con boundary anula shifting.
        shiftposfsp1.y+=(heavyphase ? 0 : massrhop*fry);
        shiftposfsp1.z+=(heavyphase ? 0 : massrhop*frz);
        shiftposfsp1.w-=(heavyphase ? 0 : massrhop*dot3);
      }

      //===== Viscosity ===== 
      if(compute) {
        const float dot=drx*dvx+dry*dvy+drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        visc=max(dot_rr2,visc);  //ViscDt=max(dot/(rr2+Eta2),ViscDt);
        //<vs_non-Newtonian>
        const float visco_NN=PHASECTE[pp2].visco;
        if(tvisco==VISCO_Artificial) {//-Artificial viscosity.
          if(dot<0) {
            const float amubar=CTE.kernelh*dot_rr2;  //amubar=CTE.kernelh*dot/(rr2+CTE.eta2);
            const float robar=(velrhop1.w+velrhop2.w)*0.5f;
            const float pi_visc=(-visco_NN*cbar*amubar/robar)*(USE_FLOATING ? ftmassp2 : massp2);
            acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
          }
        }
        else if(tvisco==VISCO_LaminarSPS||tvisco==VISCO_ConstEq) {
          {
            //vel gradients
            if(boundp2) { //this applies no slip on stress tensor
              dvx=2.f*velrhop1.x; dvy=2.f*velrhop1.y; dvz=2.f*velrhop1.z;  //fomraly I should use the moving BC vel as ug=2ub-uf
            }
            tmatrix3f dvelp1; float div_vel;
            GetVelocityGradients_FDA(rr2,drx,dry,drz,dvx,dvy,dvz,dvelp1,div_vel);

            //Strain rate tensor 
            tmatrix3f D_tensor; float div_D_tensor; float D_tensor_magn;
            float I_D,II_D; float J1_D,J2_D;
            GetStrainRateTensor(dvelp1,div_vel,I_D,II_D,J1_D,J2_D,div_D_tensor,D_tensor_magn,D_tensor);

            //Effective viscosity
            float m_NN=PHASECTE[pp2].m_NN; float n_NN=PHASECTE[pp2].n_NN; float tau_yield=PHASECTE[pp2].tau_yield;
            KerGetEta_Effective(pp1,tau_yield,D_tensor_magn,visco_NN,m_NN,n_NN,visco_etap1);
            visceta=max(visceta,visco_etap1);

            if(tvisco==VISCO_LaminarSPS){ //-Laminar contribution.
              //Morris Operator
              const float temp=2.f*(visco_etap1)/((rr2+CTE.eta2)*velrhop2.w);  //-Note this is the Morris operator and not Lo and Shao
              const float vtemp=(USE_FLOATING ? ftmassp2 : massp2)*temp*(drx*frx+dry*fry+drz*frz);
              acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;

            }
            else if(tvisco==VISCO_ConstEq) {
              //stress tensor tau 
              tmatrix3f tau_tensor; float tau_tensor_magn;
              float I_t,II_t; float J1_t,J2_t;
              GetStressTensor(D_tensor,visco_etap1,I_t,II_t,J1_t,J2_t,tau_tensor_magn,tau_tensor);

              //viscous forces
              float taux=(tau_tensor.a11*frx+tau_tensor.a12*fry+tau_tensor.a13*frz)/(velrhop2.w); //Morris 1997
              float tauy=(tau_tensor.a21*frx+tau_tensor.a22*fry+tau_tensor.a23*frz)/(velrhop2.w);
              float tauz=(tau_tensor.a31*frx+tau_tensor.a32*fry+tau_tensor.a33*frz)/(velrhop2.w);
              const float mtemp=(USE_FLOATING ? ftmassp2 : massp2);
              acep1.x+=taux*mtemp; acep1.y+=tauy*mtemp; acep1.z+=tauz*mtemp;
            }
          }
          //-SPS turbulence model.
          //-SPS turbulence model is disabled in v5.0 NN version
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction between particles for non-Newtonian models using the FDA approach. Fluid/Float-Fluid/Float or Fluid/Float-Bound.
/// Includes artificial/laminar/Const Eq. viscosity and normal/DEM floating bodies.
///
/// Realiza interaccion entre particulas para modelos no-Newtonianos que utilizan el enfoque de la FDA. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings normales/dem.
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool symm>
__global__ void KerInteractionForcesFluid_NN_FDA(unsigned n,unsigned pinit,float viscob,float viscof,float *visco_eta
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
  ,const float *ftomassp,float2 *tauff,float2 *gradvelff
  ,const float4 *poscell,const float4 *velrhop
  ,const typecode *code,const unsigned *idp
  ,float *viscdt,float *viscetadt,float *ar,float3 *ace,float *delta
  ,TpShifting shiftmode,float4 *shiftposfs)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    unsigned p1=p+pinit;      //-Number of particle.
    float visc=0,arp1=0,deltap1=0;
    float3 acep1=make_float3(0,0,0);

    //-Variables for Shifting.
    float4 shiftposfsp1;
    if(shift)shiftposfsp1=shiftposfs[p1];

    //-Obtains data of particle p1 in case there are floating bodies.		
    bool ftp1;       //-Indicates if it is floating. | Indica si es floating.
    const typecode cod=code[p1];
    if(USE_FLOATING) {
      ftp1=CODE_IsFloating(cod);
      if(ftp1 && tdensity!=DDT_None)deltap1=FLT_MAX; //-DDT is not applied to floating particles.
      if(ftp1 && shift)shiftposfsp1.x=FLT_MAX; //-Shifting is not calculated for floating bodies. | Para floatings no se calcula shifting.
    }

    //-Obtains basic data of particle p1.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrhop[p1];
    //<vs_non-Newtonian>
    const typecode pp1=CODE_GetTypeValue(cod);
    float visco_etap1=0;
    float visceta=0;

    //Obtain pressure		
    const float pressp1=cufsph::ComputePressCte_NN(velrhop1.w,PHASEARRAY[pp1].rho,PHASEARRAY[pp1].CteB,PHASEARRAY[pp1].Gamma,PHASEARRAY[pp1].Cs0,cod);
    const bool rsymp1=(symm && CEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>

    //-Variables for Laminar+SPS.
    float2 taup1_xx_xy,taup1_xz_yy,taup1_yz_zz;
    if(tvisco!=VISCO_Artificial) {
      taup1_xx_xy=tauff[p1*3];
      taup1_xz_yy=tauff[p1*3+1];
      taup1_yz_zz=tauff[p1*3+2];
    }
    //-Variables for Laminar+SPS (computation).
    float2 grap1_xx_xy,grap1_xz_yy,grap1_yz_zz;
    if(tvisco!=VISCO_Artificial) {
      grap1_xx_xy=make_float2(0,0);
      grap1_xz_yy=make_float2(0,0);
      grap1_yz_zz=make_float2(0,0);
    }

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
        KerInteractionForcesFluidBox_FDA<tker,ftmode,tvisco,tdensity,shift,false>(false,p1,pini,pfin,viscof,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,CTE.massf,pp1,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,visceta,visco_etap1,deltap1,shiftmode,shiftposfsp1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_FDA<tker,ftmode,tvisco,tdensity,shift,true >(false,p1,pini,pfin,viscof,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,CTE.massf,pp1,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,visceta,visco_etap1,deltap1,shiftmode,shiftposfsp1); //<vs_syymmetry>
      }
    }
    //-Interaction with boundaries.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin){
        KerInteractionForcesFluidBox_FDA<tker,ftmode,tvisco,tdensity,shift,false>(true,p1,pini,pfin,viscob,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,CTE.massf,pp1,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,visceta,visco_etap1,deltap1,shiftmode,shiftposfsp1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_FDA<tker,ftmode,tvisco,tdensity,shift,true >(true,p1,pini,pfin,viscob,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,CTE.massf,pp1,ftp1,pscellp1,velrhop1,pressp1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,visceta,visco_etap1,deltap1,shiftmode,shiftposfsp1);  //<vs_syymmetry>
      }
    }
    //-Stores results.
    if(shift||arp1||acep1.x||acep1.y||acep1.z||visc||visceta||visco_etap1) {
      if(tdensity!=DDT_None) {
        if(delta) {
          const float rdelta=delta[p1];
          delta[p1]=(rdelta==FLT_MAX||deltap1==FLT_MAX ? FLT_MAX : rdelta+deltap1);
        }
        else if(deltap1!=FLT_MAX)arp1+=deltap1;
      }
      ar[p1]+=arp1;
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      if(visc>viscdt[p1])viscdt[p1]=visc;
      if(visceta>viscetadt[p1])viscetadt[p1]=visceta;
      if(tvisco==VISCO_LaminarSPS) {
        float2 rg;
        rg=gradvelff[p1*3];		 rg=make_float2(rg.x+grap1_xx_xy.x,rg.y+grap1_xx_xy.y);  gradvelff[p1*3]=rg;
        rg=gradvelff[p1*3+1];  rg=make_float2(rg.x+grap1_xz_yy.x,rg.y+grap1_xz_yy.y);  gradvelff[p1*3+1]=rg;
        rg=gradvelff[p1*3+2];  rg=make_float2(rg.x+grap1_yz_zz.x,rg.y+grap1_yz_zz.y);  gradvelff[p1*3+2]=rg;
      }
      if(shift)shiftposfs[p1]=shiftposfsp1;
      //auxnn[p1] = visco_etap1; //to be used if an auxilary is needed for debug or otherwise.
    }
  }
}

//==============================================================================
/// Interaction for the force computation for non-Newtonian models using the FDA approach.
/// Interaccion para el calculo de fuerzas para modelos no-Newtonianos que utilizan el enfoque de la FDA.
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
void Interaction_ForcesGpuT_NN_FDA(const StInterParmsg &t)
{
  //-Collects kernel information.
#ifndef DISABLE_BSMODES
  if(t.kerinfo) {
    cusph::Interaction_ForcesT_KerInfo<tker,ftmode,true,tdensity,shift,false>(t.kerinfo);
    return;
  }
#endif
  const StDivDataGpu &dvd=t.divdatag;
  const int2* beginendcell=dvd.beginendcell;
  //-Interaction Fluid-Fluid & Fluid-Bound.
  if(t.fluidnum) {
    dim3 sgridf=GetSimpleGridSize(t.fluidnum,t.bsfluid);
    if(t.symmetry) //<vs_syymmetry_ini>
      KerInteractionForcesFluid_NN_FDA<tker,ftmode,tvisco,tdensity,shift,true ><<<sgridf,t.bsfluid,0,t.stm>>>
      (t.fluidnum,t.fluidini,t.viscob,t.viscof,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
        ,t.ftomassp,(float2*)t.tau,(float2*)t.gradvel,t.poscell,t.velrhop,t.code,t.idp
        ,t.viscdt,t.viscetadt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);
    else //<vs_syymmetry_end>
      KerInteractionForcesFluid_NN_FDA<tker,ftmode,tvisco,tdensity,shift,false><<<sgridf,t.bsfluid,0,t.stm>>>
      (t.fluidnum,t.fluidini,t.viscob,t.viscof,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
        ,t.ftomassp,(float2*)t.tau,(float2*)t.gradvel,t.poscell,t.velrhop,t.code,t.idp
        ,t.viscdt,t.viscetadt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);
  }
  //-Interaction Boundary-Fluid.
  if(t.boundnum) {
    const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
    dim3 sgridb=GetSimpleGridSize(t.boundnum,t.bsbound);
    //printf("bsbound:%u\n",bsbound);
    if(t.symmetry) //<vs_syymmetry_ini>
      KerInteractionForcesBound_NN<tker,ftmode,true ><<<sgridb,t.bsbound,0,t.stm>>>
      (t.boundnum,t.boundini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcell+dvd.cellfluid,t.dcell
        ,t.ftomassp,t.poscell,t.velrhop,t.code,t.idp,t.viscdt,t.ar);
    else //<vs_syymmetry_end>
      KerInteractionForcesBound_NN<tker,ftmode,false><<<sgridb,t.bsbound,0,t.stm>>>
      (t.boundnum,t.boundini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.dcell
        ,t.ftomassp,t.poscell,t.velrhop,t.code,t.idp,t.viscdt,t.ar);
  }
}
//======================END of FDA==============================================

//======================Start of SPH============================================
//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles for non-Newtonian models using the SPH approach with Const Eq. (Fluid/Float-Fluid/Float/Bound)
/// Realiza la interaccion de una particula con un conjunto de ellas para modelos no-Newtonianos que utilizan el enfoque de la SPH Const. Eq. (Fluid/Float-Fluid/Float/Bound)
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,bool symm>
__device__ void KerInteractionForcesFluidBox_SPH_ConsEq(bool boundp2,unsigned p1
  ,const unsigned &pini,const unsigned &pfin,float visco,float *visco_eta
  ,const float *ftomassp,float2 *tauff
  ,const float4 *poscell,const float4 *velrhop
  ,const typecode *code,const unsigned *idp
  ,const typecode pp1,bool ftp1
  ,const float4 &pscellp1,const float4 &velrhop1
  ,float2 &taup1_xx_xy,float2 &taup1_xz_yy,float2 &taup1_yz_zz
  ,float3 &acep1,float &visc,float &visco_etap1)
{
  for(int p2=pini; p2<pfin; p2++) {
    const float4 pscellp2=poscell[p2];
    float drx=pscellp1.x-pscellp2.x+CTE.poscellsize*(CEL_GetX(__float_as_int(pscellp1.w))-CEL_GetX(__float_as_int(pscellp2.w)));
    float dry=pscellp1.y-pscellp2.y+CTE.poscellsize*(CEL_GetY(__float_as_int(pscellp1.w))-CEL_GetY(__float_as_int(pscellp2.w)));
    float drz=pscellp1.z-pscellp2.z+CTE.poscellsize*(CEL_GetZ(__float_as_int(pscellp1.w))-CEL_GetZ(__float_as_int(pscellp2.w)));
    if(symm)dry=pscellp1.y+pscellp2.y+CTE.poscellsize*CEL_GetY(__float_as_int(pscellp2.w)); //<vs_syymmetry>
    const float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO) {
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

      //-Obtains mass of particle p2 for NN and if any floating bodies exist.
      const typecode cod=code[p2];
      const typecode pp2=(boundp2 ? pp1 : CODE_GetTypeValue(cod)); //<vs_non-Newtonian>
      float massp2=(boundp2 ? CTE.massb : PHASEARRAY[pp2].mass); //massp2 not neccesary to go in _Box function
      //Note if you masses are very different more than a ratio of 1.3 then: massp2 = (boundp2 ? PHASEARRAY[pp1].mass : PHASEARRAY[pp2].mass);

      //-Obtiene masa de particula p2 en caso de existir floatings.
      bool ftp2=false;         //-Indicates if it is floating. | Indica si es floating.
      float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
      bool compute=true; //-Deactivated when DEM is used and is float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
      if(USE_FLOATING) {
        const typecode cod=code[p2];
        ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2 ? ftomassp[CODE_GetTypeValue(cod)] : massp2);
        compute=!(USE_FTEXTERNAL && ftp1&&(boundp2||ftp2)); //-Deactivated when DEM or Chrono is used and is float-float or float-bound. | Se desactiva cuando se usa DEM o Chrono y es float-float o float-bound.
      }

      float4 velrhop2=velrhop[p2];
      if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

      //-velocity dvx.
      const float dvx=velrhop1.x-velrhop2.x,dvy=velrhop1.y-velrhop2.y,dvz=velrhop1.z-velrhop2.z;
      const float cbar=max(PHASEARRAY[pp2].Cs0,PHASEARRAY[pp2].Cs0);

      //===== Viscosity ===== 
      if(compute) {
        const float dot=drx*dvx+dry*dvy+drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        visc=max(dot_rr2,visc);  //ViscDt=max(dot/(rr2+Eta2),ViscDt);

        //<vs_non-Newtonian>				
        float2 tau_sum_xx_xy,tau_sum_xz_yy,tau_sum_yz_zz;
        float2 taup2_xx_xy=tauff[p2*3];
        float2 taup2_xz_yy=tauff[p2*3+1];
        float2 taup2_yz_zz=tauff[p2*3+2];
        //boundary particles only
        if(boundp2) {
          taup2_xx_xy=make_float2(taup1_xx_xy.x,taup1_xx_xy.y); // use (-) for slip and (+1) for no slip
          taup2_xz_yy=make_float2(taup1_xz_yy.x,taup1_xz_yy.y); //
          taup2_yz_zz=make_float2(taup1_yz_zz.x,taup1_yz_zz.y); //
        }

        tau_sum_xx_xy.x=taup1_xx_xy.x+taup2_xx_xy.x; tau_sum_xx_xy.y=taup1_xx_xy.y+taup2_xx_xy.y;	tau_sum_xz_yy.x=taup1_xz_yy.x+taup2_xz_yy.x;
        tau_sum_xz_yy.y=taup1_xz_yy.y+taup2_xz_yy.y;	tau_sum_yz_zz.x=taup1_yz_zz.x+taup2_yz_zz.x;
        tau_sum_yz_zz.y=taup1_yz_zz.y+taup2_yz_zz.y;

        float taux=(tau_sum_xx_xy.x*frx+tau_sum_xx_xy.y*fry+tau_sum_xz_yy.x*frz)/(velrhop2.w);
        float tauy=(tau_sum_xx_xy.y*frx+tau_sum_xz_yy.y*fry+tau_sum_yz_zz.x*frz)/(velrhop2.w);
        float tauz=(tau_sum_xz_yy.x*frx+tau_sum_yz_zz.x*fry+tau_sum_yz_zz.y*frz)/(velrhop2.w);
        //store stresses
        massp2=(USE_FLOATING ? ftmassp2 : massp2);
        acep1.x+=taux*massp2; acep1.y+=tauy*massp2; acep1.z+=tauz*massp2;
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction between particles for non-Newtonian models using the SPH approach with Const. Eq. Fluid/Float-Fluid/Float or Fluid/Float-Bound.
/// Includes Const. Eq. viscosity and normal/DEM floating bodies que utilizan el enfoque de la SPH Const. Eq..
///
/// Realiza interaccion entre particulas. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings normales/dem.
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,bool symm>
__global__ void KerInteractionForcesFluid_NN_SPH_ConsEq(unsigned n,unsigned pinit,float viscob,float viscof,float *visco_eta
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
  ,const float *ftomassp,float2 *tauff,float *auxnn,const float4 *poscell,const float4 *velrhop
  ,const typecode *code,const unsigned *idp,float3 *ace)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    unsigned p1=p+pinit;      //-Number of particle.			
    float3 acep1=make_float3(0,0,0);
    float visc=0;

    //-Obtains data of particle p1 in case there are floating bodies.
    //-Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1;       //-Indicates if it is floating. | Indica si es floating.
    const typecode cod=code[p1];
    if(USE_FLOATING) {
      const typecode cod=code[p1];
      ftp1=CODE_IsFloating(cod);
    }

    //-Obtains basic data of particle p1.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrhop[p1];
    const bool rsymp1=(symm && CEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>
    //<vs_non-Newtonian>
    const typecode pp1=CODE_GetTypeValue(cod);
    float visco_etap1=visco_eta[p1];

    //-Variables for tau.			
    float2 taup1_xx_xy=tauff[p1*3];
    float2 taup1_xz_yy=tauff[p1*3+1];
    float2 taup1_yz_zz=tauff[p1*3+2];

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x) {
      unsigned pini,pfin=0; cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin) {
        KerInteractionForcesFluidBox_SPH_ConsEq<tker,ftmode,tvisco,false>(false,p1,pini,pfin,viscof,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,acep1,visc,visco_etap1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_SPH_ConsEq<tker,ftmode,tvisco,true>(false,p1,pini,pfin,viscof,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,acep1,visc,visco_etap1); //<vs_syymmetry>
      }
    }
    //-Interaction with boundaries.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x) {
      unsigned pini,pfin=0; cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin) {
        KerInteractionForcesFluidBox_SPH_ConsEq<tker,ftmode,tvisco,false>(true,p1,pini,pfin,viscob,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,acep1,visc,visco_etap1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_SPH_ConsEq<tker,ftmode,tvisco,true>(true,p1,pini,pfin,viscob,visco_eta,ftomassp,tauff,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz,acep1,visc,visco_etap1); //<vs_syymmetry>
      }
    }

    //-Stores results.
    if(acep1.x||acep1.y||acep1.z) {
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      //auxnn[p1] = visco_etap1; // to be used if an auxilary is needed.
    }
  }
}


//==============================================================================
/// Calculates the strain rate tensor and effective viscocity for each particle for non-Newtonian models.
/// Calcula el tensor de la velocidad de deformacion y la viscosidad efectiva para cada particula para modelos no-Newtonianos.
//==============================================================================
template<TpFtMode ftmode,TpVisco tvisco,bool symm>
__global__ void KerInteractionForcesFluid_NN_SPH_Visco_Stress_tensor(unsigned n,unsigned pinit,float *visco_eta
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
  ,const float *ftomassp,float2 *tauff,float2 *d_tensorff,float *auxnn,const float4 *poscell,const float4 *velrhop
  ,const typecode *code,const unsigned *idp)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    unsigned p1=p+pinit;      //-Number of particle.
    //<vs_non-Newtonian>
    float visco_etap1=visco_eta[p1];;

    //Strain rate tensor 
    float2 dtsrp1_xx_xy=d_tensorff[p1*3];
    float2 dtsrp1_xz_yy=d_tensorff[p1*3+1];
    float2 dtsrp1_yz_zz=d_tensorff[p1*3+2];

    //Strain rate tensor 
    float2 taup1_xx_xy=make_float2(0,0);
    float2 taup1_xz_yy=make_float2(0,0);
    float2 taup1_yz_zz=make_float2(0,0);
    float I_t,II_t; float J1_t,J2_t; float tau_tensor_magn;
    GetStressTensor_sym(dtsrp1_xx_xy,dtsrp1_xz_yy,dtsrp1_yz_zz,visco_etap1,I_t,II_t,J1_t,J2_t,tau_tensor_magn,taup1_xx_xy,taup1_xz_yy,taup1_yz_zz);

    //-Stores results.
    if(tvisco!=VISCO_Artificial) {
      //save deformation tensor
      float2 rg;
      rg=tauff[p1*3];  rg=make_float2(rg.x+taup1_xx_xy.x,rg.y+taup1_xx_xy.y);  tauff[p1*3]=rg;
      rg=tauff[p1*3+1];  rg=make_float2(rg.x+taup1_xz_yy.x,rg.y+taup1_xz_yy.y);  tauff[p1*3+1]=rg;
      rg=tauff[p1*3+2];  rg=make_float2(rg.x+taup1_yz_zz.x,rg.y+taup1_yz_zz.y);  tauff[p1*3+2]=rg;
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles for non-Newtonian models using the SPH approach. (Fluid/Float-Fluid/Float/Bound)
/// Realiza la interaccion de una particula con un conjunto de ellas para modelos no-Newtonianos que utilizan el enfoque de la SPH. (Fluid/Float-Fluid/Float/Bound)
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,bool symm>
__device__ void KerInteractionForcesFluidBox_SPH_Morris(bool boundp2,unsigned p1
  ,const unsigned &pini,const unsigned &pfin,float visco,float *visco_eta
  ,const float *ftomassp
  ,const float4 *poscell,const float4 *velrhop
  ,const typecode *code,const unsigned *idp
  ,const typecode pp1,bool ftp1
  ,const float4 &pscellp1,const float4 &velrhop1
  ,float3 &acep1,float &visc,float &visco_etap1)
{
  for(int p2=pini; p2<pfin; p2++) {
    const float4 pscellp2=poscell[p2];
    float drx=pscellp1.x-pscellp2.x+CTE.poscellsize*(CEL_GetX(__float_as_int(pscellp1.w))-CEL_GetX(__float_as_int(pscellp2.w)));
    float dry=pscellp1.y-pscellp2.y+CTE.poscellsize*(CEL_GetY(__float_as_int(pscellp1.w))-CEL_GetY(__float_as_int(pscellp2.w)));
    float drz=pscellp1.z-pscellp2.z+CTE.poscellsize*(CEL_GetZ(__float_as_int(pscellp1.w))-CEL_GetZ(__float_as_int(pscellp2.w)));
    if(symm)dry=pscellp1.y+pscellp2.y+CTE.poscellsize*CEL_GetY(__float_as_int(pscellp2.w)); //<vs_syymmetry>
    const float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO) {
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

      //-Obtains mass of particle p2 for NN and if any floating bodies exist.
      const typecode cod=code[p2];
      const typecode pp2=(boundp2 ? pp1 : CODE_GetTypeValue(cod)); //<vs_non-Newtonian>
      float massp2=(boundp2 ? CTE.massb : PHASEARRAY[pp2].mass); //massp2 not neccesary to go in _Box function
      //Note if you masses are very different more than a ratio of 1.3 then: massp2 = (boundp2 ? PHASEARRAY[pp1].mass : PHASEARRAY[pp2].mass);

      bool ftp2=false;        //-Indicates if it is floating. | Indica si es floating.
      float ftmassp2;						//-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
      bool compute=true;			//-Deactivated when DEM is used and is float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
      if(USE_FLOATING) {
        const typecode cod=code[p2];
        ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2 ? ftomassp[CODE_GetTypeValue(cod)] : massp2);
        compute=!(USE_FTEXTERNAL && ftp1&&(boundp2||ftp2)); //-Deactivated when DEM or Chrono is used and is float-float or float-bound. | Se desactiva cuando se usa DEM o Chrono y es float-float o float-bound.
      }

      float4 velrhop2=velrhop[p2];
      if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

      //-velocity dvx.
      float dvx=velrhop1.x-velrhop2.x,dvy=velrhop1.y-velrhop2.y,dvz=velrhop1.z-velrhop2.z;
      if(boundp2) { //this applies no slip on stress tensor
        dvx=2.f*velrhop1.x; dvy=2.f*velrhop1.y; dvz=2.f*velrhop1.z;  //fomraly I should use the moving BC vel as ug=2ub-uf
      }
      const float cbar=max(PHASEARRAY[pp2].Cs0,PHASEARRAY[pp2].Cs0); //get max Cs0 of phases

      //===== Viscosity ===== 
      if(compute) {
        const float dot=drx*dvx+dry*dvy+drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        visc=max(dot_rr2,visc);  //ViscDt=max(dot/(rr2+Eta2),ViscDt);
        //<vs_non-Newtonian>
        const float visco_NN=PHASECTE[pp2].visco;
        if(tvisco==VISCO_Artificial) {//-Artificial viscosity.
          if(dot<0) {
            const float amubar=CTE.kernelh*dot_rr2;  //amubar=CTE.kernelh*dot/(rr2+CTE.eta2);
            const float robar=(velrhop1.w+velrhop2.w)*0.5f;
            const float pi_visc=(-visco_NN*cbar*amubar/robar)*(USE_FLOATING ? ftmassp2 : massp2);
            acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
          }
        }
        else if(tvisco!=VISCO_Artificial) {//-Laminar viscosity.
          {//-Laminar contribution.
            //vel gradients
            float visco_etap2=visco_eta[p2];
            //Morris Operator
            if(boundp2)visco_etap2=visco_etap1;
            //Morris Operator
            const float temp=(visco_etap1+visco_etap2)/((rr2+CTE.eta2)*velrhop2.w);
            const float vtemp=(USE_FLOATING ? ftmassp2 : massp2)*temp*(drx*frx+dry*fry+drz*frz);
            acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
          }
          //-SPS turbulence model.
          //-SPS turbulence model is disabled in v5.0 NN version
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction between particles for non-Newtonian models using the SPH approach. Fluid/Float-Fluid/Float or Fluid/Float-Bound.
/// Includes artificial/laminar viscosity and normal/DEM floating bodies.
///
/// Realiza interaccion entre particulas para modelos no-Newtonianos que utilizan el enfoque de la SPH. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings normales/dem.
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,bool symm>
__global__ void KerInteractionForcesFluid_NN_SPH_Morris(unsigned n,unsigned pinit,float viscob,float viscof,float *visco_eta
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
  ,const float *ftomassp,float *auxnn,const float4 *poscell,const float4 *velrhop
  ,const typecode *code,const unsigned *idp
  ,float3 *ace)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    unsigned p1=p+pinit;      //-Number of particle.			
    float3 acep1=make_float3(0,0,0);
    float visc=0;

    //-Obtains data of particle p1 in case there are floating bodies.
    //-Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1;       //-Indicates if it is floating. | Indica si es floating.		
    const typecode cod=code[p1];
    if(USE_FLOATING) {
      const typecode cod=code[p1];
      ftp1=CODE_IsFloating(cod);
    }

    //-Obtains basic data of particle p1.
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrhop[p1];
    const bool rsymp1=(symm && CEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>

    //<vs_non-Newtonian>
    const typecode pp1=CODE_GetTypeValue(cod);
    float visco_etap1=visco_eta[p1];

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x) {
      unsigned pini,pfin=0; cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin) {
        KerInteractionForcesFluidBox_SPH_Morris<tker,ftmode,tvisco,false>(false,p1,pini,pfin,viscof,visco_eta,ftomassp,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,acep1,visc,visco_etap1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_SPH_Morris<tker,ftmode,tvisco,true>(false,p1,pini,pfin,viscof,visco_eta,ftomassp,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,acep1,visc,visco_etap1);
      }
    }
    //-Interaction with boundaries.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x) {
      unsigned pini,pfin=0; cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin) {
        KerInteractionForcesFluidBox_SPH_Morris<tker,ftmode,tvisco,false>(true,p1,pini,pfin,viscob,visco_eta,ftomassp,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,acep1,visc,visco_etap1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_SPH_Morris<tker,ftmode,tvisco,true>(true,p1,pini,pfin,viscob,visco_eta,ftomassp,poscell,velrhop,code,idp,pp1,ftp1,pscellp1,velrhop1,acep1,visc,visco_etap1);
      }
    }
    //-Stores results.
    if(acep1.x||acep1.y||acep1.z) {
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      //auxnn[p1] = visco_etap1; // to be used if an auxilary is needed.
    }
  }
}

//==============================================================================
/// Calculates the strain rate tensor and effective viscocity for each particle
/// Calcula el tensor de la velocidad de deformacion y la viscosidad efectiva para cada particula.
//==============================================================================
template<TpFtMode ftmode,TpVisco tvisco,bool symm>
__global__ void KerInteractionForcesFluid_NN_SPH_Visco_eta(unsigned n,unsigned pinit,float viscob,float *visco_eta,const float4 *velrhop
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
  ,float2 *d_tensorff,float2 *gradvelff
  ,const typecode *code,const unsigned *idp
  ,float *viscetadt)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    unsigned p1=p+pinit;      //-Number of particle.
    //-Obtains basic data of particle p1.
    //const float4 pscellp1 = poscell[p1];
    //const float4 velrhop1 = velrhop[p1];

    //<vs_non-Newtonian>
    const typecode cod=code[p1];
    const typecode pp1=CODE_GetTypeValue(cod);
    float visco_etap1=0;

    //-Variables for gradients.
    float2 grap1_xx_xy,grap1_xz_yy,grap1_yz_zz;
    grap1_xx_xy=gradvelff[p1*3];
    grap1_xz_yy=gradvelff[p1*3+1];
    grap1_yz_zz=gradvelff[p1*3+2];

    //Strain rate tensor 
    float2 dtsrp1_xx_xy=make_float2(0,0);
    float2 dtsrp1_xz_yy=make_float2(0,0);
    float2 dtsrp1_yz_zz=make_float2(0,0);
    float div_D_tensor=0; float D_tensor_magn=0;
    float I_D,II_D; float J1_D,J2_D;
    GetStrainRateTensor_tsym(grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,I_D,II_D,J1_D,J2_D,div_D_tensor,D_tensor_magn,dtsrp1_xx_xy,dtsrp1_xz_yy,dtsrp1_yz_zz);

    //Effective viscosity
    float m_NN=PHASECTE[pp1].m_NN; float n_NN=PHASECTE[pp1].n_NN; float tau_yield=PHASECTE[pp1].tau_yield; float visco_NN=PHASECTE[pp1].visco;
    KerGetEta_Effective(pp1,tau_yield,D_tensor_magn,visco_NN,m_NN,n_NN,visco_etap1);

    //-Stores results.
    if(tvisco!=VISCO_Artificial) {
      //time step restriction
      if(visco_etap1>viscetadt[p1])viscetadt[p1]=visco_etap1; //no visceta necessary here
      //save deformation tensor
      float2 rg;
      rg=d_tensorff[p1*3];  rg=make_float2(rg.x+dtsrp1_xx_xy.x,rg.y+dtsrp1_xx_xy.y);  d_tensorff[p1*3]=rg;
      rg=d_tensorff[p1*3+1];  rg=make_float2(rg.x+dtsrp1_xz_yy.x,rg.y+dtsrp1_xz_yy.y);  d_tensorff[p1*3+1]=rg;
      rg=d_tensorff[p1*3+2];  rg=make_float2(rg.x+dtsrp1_yz_zz.x,rg.y+dtsrp1_yz_zz.y);  d_tensorff[p1*3+2]=rg;
      visco_eta[p1]=visco_etap1;
    }
    //auxnn[p1] = visco_etap1; // to be used if an auxilary is needed.
  }
}

//------------------------------------------------------------------------------
/// Interaction of a particle with a set of particles. (Fluid/Float-Fluid/Float/Bound)
/// Realiza la interaccion de una particula con un conjunto de ellas. (Fluid/Float-Fluid/Float/Bound)
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool symm>
__device__ void KerInteractionForcesFluidBox_NN_SPH_PressGrad(bool boundp2,unsigned p1
  ,const unsigned &pini,const unsigned &pfin
  ,const float *ftomassp
  ,const float4 *poscell
  ,const float4 *velrhop,const typecode *code,const unsigned *idp
  ,float massp2,const typecode pp1,bool ftp1
  ,const float4 &pscellp1,const float4 &velrhop1,float pressp1
  ,float2 &grap1_xx_xy,float2 &grap1_xz_yy,float2 &grap1_yz_zz
  ,float3 &acep1,float &arp1,float &visc,float &deltap1
  ,TpShifting shiftmode,float4 &shiftposfsp1)
{
  for(int p2=pini; p2<pfin; p2++) {
    const float4 pscellp2=poscell[p2];
    float drx=pscellp1.x-pscellp2.x+CTE.poscellsize*(CEL_GetX(__float_as_int(pscellp1.w))-CEL_GetX(__float_as_int(pscellp2.w)));
    float dry=pscellp1.y-pscellp2.y+CTE.poscellsize*(CEL_GetY(__float_as_int(pscellp1.w))-CEL_GetY(__float_as_int(pscellp2.w)));
    float drz=pscellp1.z-pscellp2.z+CTE.poscellsize*(CEL_GetZ(__float_as_int(pscellp1.w))-CEL_GetZ(__float_as_int(pscellp2.w)));
    if(symm)dry=pscellp1.y+pscellp2.y+CTE.poscellsize*CEL_GetY(__float_as_int(pscellp2.w)); //<vs_syymmetry>
    const float rr2=drx*drx+dry*dry+drz*drz;
    if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO) {
      //-Computes kernel.
      const float fac=cufsph::GetKernel_Fac<tker>(rr2);
      const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

      //-Obtains mass of particle p2 for NN and if any floating bodies exist.
      const typecode cod=code[p2];
      const typecode pp2=(boundp2 ? pp1 : CODE_GetTypeValue(cod)); //<vs_non-Newtonian>
      float massp2=(boundp2 ? CTE.massb : PHASEARRAY[pp2].mass); //massp2 not neccesary to go in _Box function
      //Note if you masses are very different more than a ratio of 1.3 then: massp2 = (boundp2 ? PHASEARRAY[pp1].mass : PHASEARRAY[pp2].mass);

      //-Obtiene masa de particula p2 en caso de existir floatings.
      bool ftp2=false;        //-Indicates if it is floating. | Indica si es floating.
      float ftmassp2;						//-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
      bool compute=true;			//-Deactivated when DEM is used and is float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
      if(USE_FLOATING) {
        const typecode cod=code[p2];
        ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2 ? ftomassp[CODE_GetTypeValue(cod)] : massp2);
#ifdef DELTA_HEAVYFLOATING
        if(ftp2 && tdensity==DDT_DDT && ftmassp2<=(massp2*1.2f))deltap1=FLT_MAX;
#else
        if(ftp2 && tdensity==DDT_DDT)deltap1=FLT_MAX;
#endif
        if(ftp2 && shift && shiftmode==SHIFT_NoBound)shiftposfsp1.x=FLT_MAX; //-Cancels shifting with floating bodies. | Con floatings anula shifting.
        compute=!(USE_FTEXTERNAL && ftp1&&(boundp2||ftp2)); //-Deactivated when DEM or Chrono is used and is float-float or float-bound. | Se desactiva cuando se usa DEM o Chrono y es float-float o float-bound.
      }
      float4 velrhop2=velrhop[p2];
      if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

      //===== Aceleration ===== 
      if(compute) {
        const float pressp2=cufsph::ComputePressCte_NN(velrhop2.w,PHASEARRAY[pp2].rho,PHASEARRAY[pp2].CteB,PHASEARRAY[pp2].Gamma,PHASEARRAY[pp2].Cs0,cod);
        const float prs=(pressp1+pressp2)/(velrhop1.w*velrhop2.w)+(tker==KERNEL_Cubic ? cufsph::GetKernelCubic_Tensil(rr2,velrhop1.w,pressp1,velrhop2.w,pressp2) : 0);
        const float p_vpm=-prs*(USE_FLOATING ? ftmassp2 : massp2);
        acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
      }

      //-Density derivative.
      float dvx=velrhop1.x-velrhop2.x,dvy=velrhop1.y-velrhop2.y,dvz=velrhop1.z-velrhop2.z;
      if(compute)arp1+=(USE_FLOATING ? ftmassp2 : massp2)*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

      const float cbar=max(PHASEARRAY[pp1].Cs0,PHASEARRAY[pp2].Cs0);
      const float dot3=(tdensity!=DDT_None||shift ? drx*frx+dry*fry+drz*frz : 0);
      //-Density derivative (DeltaSPH Molteni).
      if(tdensity==DDT_DDT && deltap1!=FLT_MAX) {
        const float rhop1over2=velrhop1.w/velrhop2.w;
        const float visc_densi=CTE.ddtkh*cbar*(rhop1over2-1.f)/(rr2+CTE.eta2);
        const float delta=(pp1==pp2 ? visc_densi*dot3*(USE_FLOATING ? ftmassp2 : massp2) : 0); //<vs_non-Newtonian>
        //deltap1=(boundp2? FLT_MAX: deltap1+delta);
        deltap1=(boundp2 && CTE.tboundary==BC_DBC ? FLT_MAX : deltap1+delta);
      }
      //-Density Diffusion Term (Fourtakas et al 2019). //<vs_dtt2_ini>
      if((tdensity==DDT_DDT2||(tdensity==DDT_DDT2Full&&!boundp2))&&deltap1!=FLT_MAX&&!ftp2) {
        const float rh=1.f+CTE.ddtgz*drz;
        const float drhop=CTE.rhopzero*pow(rh,1.f/CTE.gamma)-CTE.rhopzero;
        const float visc_densi=CTE.ddtkh*cbar*((velrhop2.w-velrhop1.w)-drhop)/(rr2+CTE.eta2);
        const float delta=(pp1==pp2 ? visc_densi*dot3*massp2/velrhop2.w : 0); //<vs_non-Newtonian>
        deltap1=(boundp2 ? FLT_MAX : deltap1-delta);
      } //<vs_dtt2_end>		

      //-Shifting correction.
      if(shift && shiftposfsp1.x!=FLT_MAX) {
        bool heavyphase=(PHASEARRAY[pp1].mass>PHASEARRAY[pp2].mass && pp1!=pp2 ? true : false); //<vs_non-Newtonian>
        const float massrhop=(USE_FLOATING ? ftmassp2 : massp2)/velrhop2.w;
        const bool noshift=(boundp2&&(shiftmode==SHIFT_NoBound||(shiftmode==SHIFT_NoFixed && CODE_IsFixed(code[p2]))));
        shiftposfsp1.x=(noshift ? FLT_MAX : (heavyphase ? 0 : shiftposfsp1.x+massrhop*frx)); //-Removes shifting for the boundaries. | Con boundary anula shifting.
        shiftposfsp1.y+=(heavyphase ? 0 : massrhop*fry);
        shiftposfsp1.z+=(heavyphase ? 0 : massrhop*frz);
        shiftposfsp1.w-=(heavyphase ? 0 : massrhop*dot3);
      }

      //===== Viscosity ===== 
      if(compute) {
        const float dot=drx*dvx+dry*dvy+drz*dvz;
        const float dot_rr2=dot/(rr2+CTE.eta2);
        visc=max(dot_rr2,visc);  //ViscDt=max(dot/(rr2+Eta2),ViscDt);

        if(tvisco!=VISCO_Artificial) { //&& !boundp2
          //vel gradients
          if(boundp2) {
            dvx=2.f*velrhop1.x; dvy=2.f*velrhop1.y; dvz=2.f*velrhop1.z;  //fomraly I should use the moving BC vel as ug=2ub-uf
          }
          GetVelocityGradients_SPH_tsym(massp2,velrhop2,dvx,dvy,dvz,frx,fry,frz,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz);
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
/// Interaction between particles for non-Newtonian models using the SPH approach. Fluid/Float-Fluid/Float or Fluid/Float-Bound.
/// Includes pressure calculations, velocity gradients and normal/DEM floating bodies.
///
/// Realiza interaccion entre particulas para modelos no-Newtonianos que utilizan el enfoque de la SPH. Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Incluye visco artificial/laminar y floatings normales/dem.
//------------------------------------------------------------------------------
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool symm>
__global__ void KerInteractionForcesFluid_NN_SPH_PressGrad(unsigned n,unsigned pinit
  ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
  ,const float *ftomassp,float2 *gradvelff
  ,const float4 *poscell
  ,const float4 *velrhop,const typecode *code,const unsigned *idp
  ,float *viscdt,float *ar,float3 *ace,float *delta
  ,TpShifting shiftmode,float4 *shiftposfs)
{
  const unsigned p=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(p<n) {
    unsigned p1=p+pinit;      //-Number of particle.
    float visc=0,arp1=0,deltap1=0;
    float3 acep1=make_float3(0,0,0);

    //-Variables for Shifting.
    float4 shiftposfsp1;
    if(shift)shiftposfsp1=shiftposfs[p1];

    //-Obtains data of particle p1 in case there are floating bodies.
    //-Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1;       //-Indicates if it is floating. | Indica si es floating.
    const typecode cod=code[p1];
    if(USE_FLOATING) {
      ftp1=CODE_IsFloating(cod);
      if(ftp1 && tdensity!=DDT_None)deltap1=FLT_MAX; //-DDT is not applied to floating particles.
      if(ftp1 && shift)shiftposfsp1.x=FLT_MAX; //-Shifting is not calculated for floating bodies. | Para floatings no se calcula shifting.
    }

    //-Obtains basic data of particle p1.		
    const float4 pscellp1=poscell[p1];
    const float4 velrhop1=velrhop[p1];
    //<vs_non-Newtonian>
    const typecode pp1=CODE_GetTypeValue(cod);

    //Obtain pressure
    const float pressp1=cufsph::ComputePressCte_NN(velrhop1.w,PHASEARRAY[pp1].rho,PHASEARRAY[pp1].CteB,PHASEARRAY[pp1].Gamma,PHASEARRAY[pp1].Cs0,cod);
    const bool rsymp1=(symm && CEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>

    //-Variables for vel gradients
    float2 grap1_xx_xy,grap1_xz_yy,grap1_yz_zz;
    if(tvisco!=VISCO_Artificial) {
      grap1_xx_xy=make_float2(0,0);
      grap1_xz_yy=make_float2(0,0);
      grap1_yz_zz=make_float2(0,0);
    }

    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    ini3+=cellfluid; fin3+=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x) {
      unsigned pini,pfin=0; cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin) {
        KerInteractionForcesFluidBox_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift,false>(false,p1,pini,pfin,ftomassp,poscell,velrhop,code,idp,CTE.massf,pp1,ftp1,pscellp1,velrhop1,pressp1,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,deltap1,shiftmode,shiftposfsp1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift,true >(false,p1,pini,pfin,ftomassp,poscell,velrhop,code,idp,CTE.massf,pp1,ftp1,pscellp1,velrhop1,pressp1,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,deltap1,shiftmode,shiftposfsp1); //<vs_syymmetry>
      }
    }
    //-Interaction with boundaries.
    ini3-=cellfluid; fin3-=cellfluid;
    for(int c3=ini3; c3<fin3; c3+=nc.w)for(int c2=ini2; c2<fin2; c2+=nc.x) {
      unsigned pini,pfin=0; cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
      if(pfin) {
        KerInteractionForcesFluidBox_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift,false>(true,p1,pini,pfin,ftomassp,poscell,velrhop,code,idp,CTE.massb,pp1,ftp1,pscellp1,velrhop1,pressp1,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,deltap1,shiftmode,shiftposfsp1);
        if(symm && rsymp1)	KerInteractionForcesFluidBox_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift,true >(true,p1,pini,pfin,ftomassp,poscell,velrhop,code,idp,CTE.massb,pp1,ftp1,pscellp1,velrhop1,pressp1,grap1_xx_xy,grap1_xz_yy,grap1_yz_zz,acep1,arp1,visc,deltap1,shiftmode,shiftposfsp1); //<vs_syymmetry>
      }
    }
    //-Stores results.
    if(shift||arp1||acep1.x||acep1.y||acep1.z||visc) {
      if(tdensity!=DDT_None) {
        if(delta) {
          const float rdelta=delta[p1];
          delta[p1]=(rdelta==FLT_MAX||deltap1==FLT_MAX ? FLT_MAX : rdelta+deltap1);
        }
        else if(deltap1!=FLT_MAX)arp1+=deltap1;
      }
      ar[p1]+=arp1;
      float3 r=ace[p1]; r.x+=acep1.x; r.y+=acep1.y; r.z+=acep1.z; ace[p1]=r;
      if(visc>viscdt[p1])viscdt[p1]=visc;
      if(tvisco!=VISCO_Artificial) {
        float2 rg;
        rg=gradvelff[p1*3];		 rg=make_float2(rg.x+grap1_xx_xy.x,rg.y+grap1_xx_xy.y);  gradvelff[p1*3]=rg;
        rg=gradvelff[p1*3+1];  rg=make_float2(rg.x+grap1_xz_yy.x,rg.y+grap1_xz_yy.y);  gradvelff[p1*3+1]=rg;
        rg=gradvelff[p1*3+2];  rg=make_float2(rg.x+grap1_yz_zz.x,rg.y+grap1_yz_zz.y);  gradvelff[p1*3+2]=rg;
      }
      if(shift)shiftposfs[p1]=shiftposfsp1;
      //auxnn[p1] = visco_etap1; // to be used if an auxilary is needed.
    }
  }
}

//==============================================================================
/// Interaction for the force computation using the SPH approach.
/// Interaccion para el calculo de fuerzas que utilizan el enfoque de la SPH .
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
void Interaction_ForcesGpuT_NN_SPH(const StInterParmsg &t)
{
  //-Collects kernel information.
#ifndef DISABLE_BSMODES
  if(t.kerinfo) {
    cusph::Interaction_ForcesT_KerInfo<tker,ftmode,true,tdensity,shift,false>(t.kerinfo);
    return;
  }
#endif
  const StDivDataGpu &dvd=t.divdatag;
  const int2* beginendcell=dvd.beginendcell;
  //-Interaction Fluid-Fluid & Fluid-Bound.
  if(t.fluidnum) {
    dim3 sgridf=GetSimpleGridSize(t.fluidnum,t.bsfluid);
    if(t.symmetry){ //<vs_syymmetry_ini>
      KerInteractionForcesFluid_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift,true ><<<sgridf,t.bsfluid,0,t.stm>>>
        (t.fluidnum,t.fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
          ,t.ftomassp,(float2*)t.gradvel,t.poscell,t.velrhop,t.code,t.idp
          ,t.viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);

      KerInteractionForcesFluid_NN_SPH_Visco_eta<ftmode,tvisco,true ><<<sgridf,t.bsfluid,0,t.stm>>>
        (t.fluidnum,t.fluidini,t.viscob,t.visco_eta,t.velrhop,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
          ,(float2*)t.d_tensor,(float2*)t.gradvel,t.code,t.idp
          ,t.viscetadt);
      //choice of visc formulation
      if(tvisco!=VISCO_ConstEq) KerInteractionForcesFluid_NN_SPH_Morris<tker,ftmode,tvisco,true ><<<sgridf,t.bsfluid,0,t.stm>>>
        (t.fluidnum,t.fluidini,t.viscob,t.viscof,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
          ,t.ftomassp,t.auxnn,t.poscell,t.velrhop,t.code,t.idp
          ,t.ace);
      else {
        // Build stress tensor
        KerInteractionForcesFluid_NN_SPH_Visco_Stress_tensor<ftmode,tvisco,true ><<<sgridf,t.bsfluid,0,t.stm>>>
          (t.fluidnum,t.fluidini,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
            ,t.ftomassp,(float2*)t.tau,(float2*)t.d_tensor,t.auxnn,t.poscell,t.velrhop,t.code,t.idp);
        //Get stresses
        KerInteractionForcesFluid_NN_SPH_ConsEq<tker,ftmode,tvisco,true ><<<sgridf,t.bsfluid,0,t.stm>>>
          (t.fluidnum,t.fluidini,t.viscob,t.viscof,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
            ,t.ftomassp,(float2*)t.tau,t.auxnn,t.poscell,t.velrhop,t.code,t.idp
            ,t.ace);
      }

    }
    else {//<vs_syymmetry_end>			
      KerInteractionForcesFluid_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift,false ><<<sgridf,t.bsfluid,0,t.stm>>>
        (t.fluidnum,t.fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
          ,t.ftomassp,(float2*)t.gradvel,t.poscell,t.velrhop,t.code,t.idp
          ,t.viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);

      KerInteractionForcesFluid_NN_SPH_Visco_eta<ftmode,tvisco,false ><<<sgridf,t.bsfluid,0,t.stm>>>
        (t.fluidnum,t.fluidini,t.viscob,t.visco_eta,t.velrhop,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
          ,(float2*)t.d_tensor,(float2*)t.gradvel,t.code,t.idp
          ,t.viscetadt);
      //choice of visc formulation
      if(tvisco!=VISCO_ConstEq) KerInteractionForcesFluid_NN_SPH_Morris<tker,ftmode,tvisco,false ><<<sgridf,t.bsfluid,0,t.stm>>>
        (t.fluidnum,t.fluidini,t.viscob,t.viscof,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
          ,t.ftomassp,t.auxnn,t.poscell,t.velrhop,t.code,t.idp
          ,t.ace);
      else {
        // Build stress tensor				
        KerInteractionForcesFluid_NN_SPH_Visco_Stress_tensor<ftmode,tvisco,false ><<<sgridf,t.bsfluid,0,t.stm>>>
          (t.fluidnum,t.fluidini,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
            ,t.ftomassp,(float2*)t.tau,(float2*)t.d_tensor,t.auxnn,t.poscell,t.velrhop,t.code,t.idp);
        //Get stresses
        KerInteractionForcesFluid_NN_SPH_ConsEq<tker,ftmode,tvisco,false ><<<sgridf,t.bsfluid,0,t.stm>>>
          (t.fluidnum,t.fluidini,t.viscob,t.viscof,t.visco_eta,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,t.dcell
            ,t.ftomassp,(float2*)t.tau,t.auxnn,t.poscell,t.velrhop,t.code,t.idp
            ,t.ace);
      }
    }
  }

  //-Interaction Boundary-Fluid.
  if(t.boundnum) {
    const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
    dim3 sgridb=GetSimpleGridSize(t.boundnum,t.bsbound);
    //printf("bsbound:%u\n",bsbound);
    if(t.symmetry) //<vs_syymmetry_ini>
      KerInteractionForcesBound_NN<tker,ftmode,true ><<<sgridb,t.bsbound,0,t.stm>>>
      (t.boundnum,t.boundini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcell+dvd.cellfluid,t.dcell
        ,t.ftomassp,t.poscell,t.velrhop,t.code,t.idp,t.viscdt,t.ar);
    else //<vs_syymmetry_end>
      KerInteractionForcesBound_NN<tker,ftmode,false><<<sgridb,t.bsbound,0,t.stm>>>
      (t.boundnum,t.boundini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.dcell
        ,t.ftomassp,t.poscell,t.velrhop,t.code,t.idp,t.viscdt,t.ar);
  }
}
//======================END of SPH==============================================

//======================Start of non-Newtonian Templates=======================================
//Uncomment for fast compile 
//#define FAST_COMPILATION
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift> void Interaction_ForcesGpuT_NN(const StInterParmsg &t) {
#ifdef FAST_COMPILATION
  if(t.tvelgrad!=VELGRAD_FDA)throw "Extra SPH Gradients options are disabled for FastCompilation...";
  Interaction_ForcesGpuT_NN_FDA	    < tker,ftmode,tvisco,tdensity,shift>(t);
#else	
  if(t.tvelgrad==VELGRAD_FDA) Interaction_ForcesGpuT_NN_FDA	    < tker,ftmode,tvisco,tdensity,shift>(t);
  else if(t.tvelgrad==VELGRAD_SPH)	Interaction_ForcesGpuT_NN_SPH		< tker,ftmode,tvisco,tdensity,shift>(t);
#endif
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco> void Interaction_ForcesNN_gt2(const StInterParmsg &t) {
#ifdef FAST_COMPILATION
  if(!t.shiftmode||t.tdensity!=DDT_DDT2Full)throw "Shifting and extra DDT are disabled for FastCompilation...";
  Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_DDT2Full,true>(t);
#else
  if(t.shiftmode) {
    const bool shift=true;
    if(t.tdensity==DDT_None)    Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_None,shift>(t);
    if(t.tdensity==DDT_DDT)     Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_DDT,shift>(t);
    if(t.tdensity==DDT_DDT2)    Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_DDT2,shift>(t);  //<vs_dtt2>
    if(t.tdensity==DDT_DDT2Full)Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_DDT2Full,shift>(t);  //<vs_dtt2>
  }
  else {
    const bool shift=false;
    if(t.tdensity==DDT_None)    Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_None,shift>(t);
    if(t.tdensity==DDT_DDT)     Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_DDT,shift>(t);
    if(t.tdensity==DDT_DDT2)    Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_DDT2,shift>(t);  //<vs_dtt2>
    if(t.tdensity==DDT_DDT2Full)Interaction_ForcesGpuT_NN<tker,ftmode,tvisco,DDT_DDT2Full,shift>(t);  //<vs_dtt2>
  }
#endif
}
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> void Interaction_ForcesNN_gt1(const StInterParmsg &t) {
  //GFCheck how to add fast compilation of laminar viscosity
#ifdef FAST_COMPILATION
  if(t.tvisco!=VISCO_LaminarSPS)throw "Extra viscosity options are disabled for FastCompilation...";
  Interaction_ForcesNN_gt2<tker,ftmode,VISCO_LaminarSPS>(t);
#else
  if(t.tvisco==VISCO_ConstEq)		      Interaction_ForcesNN_gt2<tker,ftmode,VISCO_ConstEq>(t);
  else if(t.tvisco==VISCO_LaminarSPS)	Interaction_ForcesNN_gt2<tker,ftmode,VISCO_LaminarSPS>(t);
  else if(t.tvisco==VISCO_Artificial)	Interaction_ForcesNN_gt2<tker,ftmode,VISCO_Artificial>(t);
#endif
}
//==============================================================================
template<TpKernel tker> void Interaction_ForcesNN_gt0(const StInterParmsg &t) {
#ifdef FAST_COMPILATION
  if(t.ftmode!=FTMODE_None)throw "Extra FtMode options are disabled for FastCompilation...";
  Interaction_ForcesNN_gt1<tker,FTMODE_None>(t);
#else
  if(t.ftmode==FTMODE_None)    Interaction_ForcesNN_gt1<tker,FTMODE_None>(t);
  else if(t.ftmode==FTMODE_Sph)Interaction_ForcesNN_gt1<tker,FTMODE_Sph>(t);
  else if(t.ftmode==FTMODE_Ext)Interaction_ForcesNN_gt1<tker,FTMODE_Ext>(t);
#endif
}
//==============================================================================
void Interaction_ForcesNN(const StInterParmsg &t) {
#ifdef FAST_COMPILATION
  if(t.tkernel!=KERNEL_Wendland)throw "Extra kernels are disabled for FastCompilation...";
  Interaction_ForcesNN_gt0<KERNEL_Wendland>(t);
#else
  if(t.tkernel==KERNEL_Wendland)     Interaction_ForcesNN_gt0<KERNEL_Wendland>(t);
#ifndef DISABLE_KERNELS_EXTRA
  else if(t.tkernel==KERNEL_Cubic)   Interaction_ForcesNN_gt0<KERNEL_Cubic   >(t);
#endif
#endif
}
//======================End of NN Templates=======================================

} //end of file