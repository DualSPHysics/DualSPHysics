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

/// \file JSphCpu.cpp \brief Implements the class \ref JSphCpu.

#include "JSphCpu.h"
#include "JCellDivCpu.h"
#include "JPartFloatBi4.h"
#include "Functions.h"
#include "JDsMotion.h"
#include "JArraysCpu.h"
#include "JDsFixedDt.h"
#include "JWaveGen.h"
#include "JMLPistons.h"     //<vs_mlapiston>
#include "JRelaxZones.h"    //<vs_rzone>
#include "JChronoObjects.h" //<vs_chroono>
#include "JDsDamping.h"
#include "JXml.h"
#include "JDsSaveDt.h"
#include "JDsOutputTime.h"
#include "JDsAccInput.h"
#include "JDsGaugeSystem.h"
#include "JSphBoundCorr.h"  //<vs_innlet>
#include <climits>

using namespace std;

//==============================================================================
/// Perform interaction between particles for NN using the SPH approach. Bound-Fluid/Float
//==============================================================================
template<TpKernel tker,TpFtMode ftmode> void JSphCpu::InteractionForcesBound_NN_SPH
(unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
  ,float &viscdt,float *ar)const
{
  //-Initialize viscth to calculate max viscdt with OpenMP. | Inicializa viscth para calcular visdt maximo con OpenMP.
  float viscth[OMP_MAXTHREADS*OMP_STRIDE];
  for(int th=0; th<OmpThreads; th++)viscth[th*OMP_STRIDE]=0;
  //-Starts execution using OpenMP.
  const int pfin=int(pinit+n);
#ifdef OMP_USE
#pragma omp parallel for schedule (guided)
#endif
  for(int p1=int(pinit); p1<pfin; p1++) {
    float visc=0,arp1=0;

    //-Load data of particle p1. | Carga datos de particula p1.
    const tdouble3 posp1=pos[p1];
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>
    const tfloat4 velrhop1=velrhop[p1];

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],false,divdata);
    for(int z=ngs.zini; z<ngs.zfin; z++)for(int y=ngs.yini; y<ngs.yfin; y++) {
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
      //---------------------------------------------------------------------------------------------
      bool rsym=false; //<vs_syymmetry>
      for(unsigned p2=pif.x; p2<pif.y; p2++) {
        const float drx=float(posp1.x-pos[p2].x);
        float dry=float(posp1.y-pos[p2].y);
        if(rsym)    dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
        const float drz=float(posp1.z-pos[p2].z);
        const float rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<=KernelSize2 && rr2>=ALMOSTZERO) {
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

                                                                   //===== Get mass of particle p2 ===== 
          float massp2=MassFluid; //-Contains particle mass of incorrect fluid. | Contiene masa de particula por defecto fluid.
          bool compute=true;      //-Deactivate when using DEM and/or bound-float. | Se desactiva cuando se usa DEM y es bound-float.
          if(USE_FLOATING) {
            bool ftp2=CODE_IsFloating(code[p2]);
            if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
            compute=!(USE_FTEXTERNAL && ftp2); //-Deactivate when using DEM/Chrono and/or bound-float. | Se desactiva cuando se usa DEM/Chrono y es bound-float.
          }

          if(compute) {
            //-Density derivative (Continuity equation).
            tfloat4 velrhop2=velrhop[p2];
            if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>
            const float dvx=velrhop1.x-velrhop2.x,dvy=velrhop1.y-velrhop2.y,dvz=velrhop1.z-velrhop2.z;
            if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz)*(velrhop1.w/velrhop2.w);

            {//-Viscosity.
              const float dot=drx*dvx+dry*dvy+drz*dvz;
              const float dot_rr2=dot/(rr2+Eta2);
              visc=max(dot_rr2,visc);
            }
          }
          rsym=(rsymp1&&!rsym && float(posp1.y-dry)<=KernelSize); //<vs_syymmetry>
          if(rsym)p2--;                                             //<vs_syymmetry>
        }
        else rsym=false;                                            //<vs_syymmetry>
      }
    }
    //-Sum results together. | Almacena resultados.
    if(arp1||visc) {
      ar[p1]+=arp1;
      const int th=omp_get_thread_num();
      if(visc>viscth[th*OMP_STRIDE])viscth[th*OMP_STRIDE]=visc;
    }
  }
  //-Keep max value in viscdt. | Guarda en viscdt el valor maximo.
  for(int th=0; th<OmpThreads; th++)if(viscdt<viscth[th*OMP_STRIDE])viscdt=viscth[th*OMP_STRIDE];
}

//==============================================================================
/// Perform calculations for NN using the SPH approach: Calculates the Stress tensor
//==============================================================================
template<TpFtMode ftmode,TpVisco tvisco> void JSphCpu::InteractionForcesFluid_NN_SPH_Visco_Stress_tensor
(unsigned n,unsigned pinit,float visco,float *visco_eta
  ,tsymatrix3f *tau,const tsymatrix3f *D_tensor,float *auxnn
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp)const
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
#ifdef OMP_USE
#pragma omp parallel for schedule (guided)
#endif
  for(int p1=int(pinit); p1<pfin; p1++) {

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];

    //stress related to p1
    const float visco_etap1=visco_eta[p1];
    const tsymatrix3f D_tensorp1=D_tensor[p1];
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>
    const typecode pp1=CODE_GetTypeValue(code[p1]);

    //setup stress related variables
    tsymatrix3f tau_tensorp1={0,0,0,0,0,0}; float tau_tensor_magn=0.f;
    float I_t,II_t=0.f; float J1_t,J2_t=0.f;
    GetStressTensor_sym(D_tensorp1,visco_etap1,I_t,II_t,J1_t,J2_t,tau_tensor_magn,tau_tensorp1);

    //debug uncomment if you want an output of effective viscosity
    //if(tau_tensor_magn){
    //  auxnn[p1] = visco_etap1; // D_tensor_magn;
    //}

    ////-Sum results together. | Almacena resultados.
    if(tvisco!=VISCO_Artificial) {
      tau[p1].xx=tau_tensorp1.xx; tau[p1].xy=tau_tensorp1.xy; tau[p1].xz=tau_tensorp1.xz;
      tau[p1].yy=tau_tensorp1.yy; tau[p1].yz=tau_tensorp1.yz;
      tau[p1].zz=tau_tensorp1.zz;
    }
  }
}

//==============================================================================
/// Perform calculations for NN using the SPH approach: Calculates the effective viscosity
//==============================================================================
template<TpFtMode ftmode,TpVisco tvisco> void JSphCpu::InteractionForcesFluid_NN_SPH_Visco_eta
(unsigned n,unsigned pinit,float visco,float *visco_eta,const tfloat4 *velrhop
  ,const tsymatrix3f* gradvel,tsymatrix3f* D_tensor,float *auxnn,float &viscetadt
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp)const
{
  //-Initialize viscth to calculate viscdt maximo con OpenMP. | Inicializa viscth para calcular visdt maximo con OpenMP.
  float viscetath[OMP_MAXTHREADS*OMP_STRIDE];
  for(int th=0; th<OmpThreads; th++)viscetath[th*OMP_STRIDE]=0;
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
#ifdef OMP_USE
#pragma omp parallel for schedule (guided)
#endif
  for(int p1=int(pinit); p1<pfin; p1++) {

    //-Obtain data of particle p1.
    const tsymatrix3f gradvelp1=gradvel[p1];
    const tdouble3 posp1=pos[p1];
    const float rhopp1=velrhop[p1].w;
    //const tsymatrix3f taup1 = (tvisco == VISCO_Artificial ? gradvelp1 : tau[p1]);
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>
    const typecode pp1=CODE_GetTypeValue(code[p1]); //<vs_non-Newtonian>

    //Strain rate tensor 
    tsymatrix3f D_tensorp1={0,0,0,0,0,0};
    float div_D_tensor=0; float D_tensor_magn=0;
    float I_D,II_D=0.f; float J1_D,J2_D=0.f;
    GetStrainRateTensor_tsym(gradvelp1,I_D,II_D,J1_D,J2_D,div_D_tensor,D_tensor_magn,D_tensorp1);

    //Effective viscosity
    float visco_etap1=0.f;
    float m_NN=PhaseCte[pp1].m_NN; float n_NN=PhaseCte[pp1].n_NN; float tau_yield=PhaseCte[pp1].tau_yield; float visco_NN=PhaseCte[pp1].visco;
    GetEta_Effective(pp1,tau_yield,D_tensor_magn,visco_NN,m_NN,n_NN,visco_etap1);

    //-Sum results together. | Almacena resultados.
    if(visco_etap1) {
      visco_eta[p1]=visco_etap1;
    }
    if(tvisco!=VISCO_Artificial) {
      D_tensor[p1].xx=D_tensorp1.xx; D_tensor[p1].xy=D_tensorp1.xy; D_tensor[p1].xz=D_tensorp1.xz;
      D_tensor[p1].yy=D_tensorp1.yy; D_tensor[p1].yz=D_tensorp1.yz;
      D_tensor[p1].zz=D_tensorp1.zz;
    }
    //debug
    //if(D_tensor_magn){
    //  auxnn[p1] = visco_etap1; // d_tensor_magn;
    //} 
    const int th=omp_get_thread_num();
    if(visco_etap1>viscetath[th*OMP_STRIDE])viscetath[th*OMP_STRIDE]=visco_etap1;
  }
  for(int th=0; th<OmpThreads; th++)if(viscetadt<viscetath[th*OMP_STRIDE])viscetadt=viscetath[th*OMP_STRIDE];
}

//==============================================================================
/// Perform interaction between particles for the SPH approcach using the Const. Eq.: Fluid/Float-Fluid/Float or Fluid/Float-Bound 
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
void JSphCpu::InteractionForcesFluid_NN_SPH_ConsEq(unsigned n,unsigned pinit,bool boundp2,float visco
  ,StDivDataCpu divdata,const unsigned *dcell
  ,float *visco_eta,const tsymatrix3f* tau,float *auxnn
  ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
  ,tfloat3 *ace)const
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
#ifdef OMP_USE
#pragma omp parallel for schedule (guided)
#endif
  for(int p1=int(pinit); p1<pfin; p1++) {
    tfloat3 acep1=TFloat3(0);
    float visc=0;
    //-Obtain data of particle p1 in case of floating objects. | Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1=false;     //-Indicate if it is floating. | Indica si es floating.
    if(USE_FLOATING) {
      ftp1=CODE_IsFloating(code[p1]);
    }

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];
    const tfloat3 velp1=TFloat3(velrhop[p1].x,velrhop[p1].y,velrhop[p1].z);
    const float rhopp1=velrhop[p1].w;
    const float visco_etap1=visco_eta[p1]; //<vs_non-Newtonian>
    const tsymatrix3f tau_tensorp1=tau[p1];
    const typecode pp1=(CODE_GetTypeValue(code[p1])); //<vs_non-Newtonian>
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
    for(int z=ngs.zini; z<ngs.zfin; z++)for(int y=ngs.yini; y<ngs.yfin; y++) {
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
      //------------------------------------------------------------------------------------------------
      bool rsym=false; //<vs_syymmetry>
      for(unsigned p2=pif.x; p2<pif.y; p2++) {
        const float drx=float(posp1.x-pos[p2].x);
        float dry=float(posp1.y-pos[p2].y);
        if(rsym)   dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
        const float drz=float(posp1.z-pos[p2].z);
        const float rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<=KernelSize2 && rr2>=ALMOSTZERO) {
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass of particle p2 ===== 
          //<vs_non-Newtonian>
          const typecode pp2=(boundp2 ? pp1 : CODE_GetTypeValue(code[p2])); //<vs_non-Newtonian>
          float massp2=(boundp2 ? MassBound : PhaseArray[pp2].mass); //-Contiene masa de particula segun sea bound o fluid.
          //Note if you masses are very different more than a ratio of 1.3 then: massp2 = (boundp2 ? PhaseArray[pp1].mass : PhaseArray[pp2].mass);

          //Floating
          bool ftp2=false;    //-Indicate if it is floating | Indica si es floating.
          bool compute=true;  //-Deactivate when using DEM and if it is of type float-float or float-bound | Se desactiva cuando se usa DEM y es float-float o float-bound.
          if(USE_FLOATING) {
            ftp2=CODE_IsFloating(code[p2]);
            if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
            compute=!(USE_FTEXTERNAL && ftp1&&(boundp2||ftp2)); //-Deactivate when using DEM and if it is of type float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
          }

          tfloat4 velrhop2=velrhop[p2];
          if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

          //-velocity dvx.
          const float dvx=velp1.x-velrhop2.x,dvy=velp1.y-velrhop2.y,dvz=velp1.z-velrhop2.z;
          const float cbar=max(PhaseArray[pp2].Cs0,PhaseArray[pp2].Cs0);

          //===== Viscosity ===== 
          if(compute) {
            const float dot=drx*dvx+dry*dvy+drz*dvz;
            const float dot_rr2=dot/(rr2+Eta2);
            visc=max(dot_rr2,visc);

            tsymatrix3f tau_tensorp2=tau[p2]; tsymatrix3f tau_sum={0,0,0,0,0,0};
            if(boundp2)tau_tensorp2=tau_tensorp1;
            tau_sum.xx=tau_tensorp1.xx+tau_tensorp2.xx;	tau_sum.xy=tau_tensorp1.xy+tau_tensorp2.xy;	tau_sum.xz=tau_tensorp1.xz+tau_tensorp2.xz;
            tau_sum.yy=tau_tensorp1.yy+tau_tensorp2.yy;	tau_sum.yz=tau_tensorp1.yz+tau_tensorp2.yz;
            tau_sum.zz=tau_tensorp1.zz+tau_tensorp2.zz;

            float taux=(tau_sum.xx*frx+tau_sum.xy*fry+tau_sum.xz*frz)/(velrhop2.w); // as per symetric tensor grad
            float tauy=(tau_sum.xy*frx+tau_sum.yy*fry+tau_sum.yz*frz)/(velrhop2.w);
            float tauz=(tau_sum.xz*frx+tau_sum.yz*fry+tau_sum.zz*frz)/(velrhop2.w);
            //store stresses
            acep1.x+=taux*massp2; acep1.y+=tauy*massp2; acep1.z+=tauz*massp2;
          }
          //-SPS turbulence model.
          //-SPS turbulence model is disabled in v5.0 NN version
          rsym=(rsymp1&&!rsym && float(posp1.y-dry)<=KernelSize); //<vs_syymmetry>
          if(rsym)p2--;																									//<vs_syymmetry>
        }
        else rsym=false;																								//<vs_syymmetry>		
      }
    }
    //-Sum results together. | Almacena resultados.
    if(acep1.x||acep1.y||acep1.z) {
      ace[p1]=ace[p1]+acep1;
    }
  }
}


//==============================================================================
/// Perform interaction between particles for the SPH approcach using the Morris operator: Fluid/Float-Fluid/Float or Fluid/Float-Bound 
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
void JSphCpu::InteractionForcesFluid_NN_SPH_Morris(unsigned n,unsigned pinit,bool boundp2,float visco
  ,StDivDataCpu divdata,const unsigned *dcell
  ,float *visco_eta,const tsymatrix3f* tau,tsymatrix3f* gradvel,float *auxnn
  ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
  ,tfloat3 *ace)const
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
#ifdef OMP_USE
#pragma omp parallel for schedule (guided)
#endif
  for(int p1=int(pinit); p1<pfin; p1++) {
    tfloat3 acep1=TFloat3(0);
    float visc=0;

    //-Obtain data of particle p1 in case of floating objects. | Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1=false;     //-Indicate if it is floating. | Indica si es floating.
    if(USE_FLOATING) {
      ftp1=CODE_IsFloating(code[p1]);
    }

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];
    const tfloat3 velp1=TFloat3(velrhop[p1].x,velrhop[p1].y,velrhop[p1].z);
    const float rhopp1=velrhop[p1].w;
    const float visco_etap1=visco_eta[p1];
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>
    const typecode pp1=CODE_GetTypeValue(code[p1]); //<vs_non-Newtonian>

                                                      //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
    for(int z=ngs.zini; z<ngs.zfin; z++)for(int y=ngs.yini; y<ngs.yfin; y++) {
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
      //------------------------------------------------------------------------------------------------
      bool rsym=false; //<vs_syymmetry>
      for(unsigned p2=pif.x; p2<pif.y; p2++) {
        const float drx=float(posp1.x-pos[p2].x);
        float dry=float(posp1.y-pos[p2].y);
        if(rsym)    dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
        const float drz=float(posp1.z-pos[p2].z);
        const float rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<=KernelSize2 && rr2>=ALMOSTZERO){
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass of particle p2 ===== 
          //<vs_non-Newtonian>
          const typecode pp2=(boundp2 ? pp1 : CODE_GetTypeValue(code[p2])); //<vs_non-Newtonian>
          float massp2=(boundp2 ? MassBound : PhaseArray[pp2].mass); //-Contiene masa de particula segun sea bound o fluid.
          //Note if you masses are very different more than a ratio of 1.3 then: massp2 = (boundp2 ? PhaseArray[pp1].mass : PhaseArray[pp2].mass);

          //Floating
          bool ftp2=false;    //-Indicate if it is floating | Indica si es floating.
          bool compute=true;  //-Deactivate when using DEM and if it is of type float-float or float-bound | Se desactiva cuando se usa DEM y es float-float o float-bound.
          if(USE_FLOATING) {
            ftp2=CODE_IsFloating(code[p2]);
            if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
            compute=!(USE_FTEXTERNAL && ftp1&&(boundp2||ftp2)); //-Deactivate when using DEM and if it is of type float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
          }

          tfloat4 velrhop2=velrhop[p2];
          if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>                                               
          //-velocity dvx.
          float dvx=velp1.x-velrhop2.x,dvy=velp1.y-velrhop2.y,dvz=velp1.z-velrhop2.z;
          if(boundp2) { //this applies no slip on tensor                                     
            dvx=2.f*velp1.x; dvy=2.f*velp1.y; dvz=2.f*velp1.z;  //fomraly I should use the moving BC vel as ug=2ub-uf
          }
          const float cbar=max(PhaseArray[pp2].Cs0,PhaseArray[pp2].Cs0);

          //===== Viscosity ===== 
          if(compute) {
            const float dot=drx*dvx+dry*dvy+drz*dvz;
            const float dot_rr2=dot/(rr2+Eta2);
            visc=max(dot_rr2,visc);
            const float visco_NN=PhaseCte[pp2].visco;
            if(tvisco==VISCO_Artificial) {//-Artificial viscosity.
              if(dot<0) {
                const float amubar=KernelH*dot_rr2;  //amubar=CTE.h*dot/(rr2+CTE.eta2);
                const float robar=(rhopp1+velrhop2.w)*0.5f;
                const float pi_visc=(-visco_NN*cbar*amubar/robar)*massp2;
                acep1.x-=pi_visc*frx; acep1.y-=pi_visc*fry; acep1.z-=pi_visc*frz;
              }
            }
            else if(tvisco==VISCO_LaminarSPS) {//-Laminar viscosity.
              {//-Laminar contribution.
                float visco_etap2=visco_eta[p2];
                //Morris Operator
                if(boundp2)visco_etap2=visco_etap1;
                const float temp=(visco_etap1+visco_etap2)/((rr2+Eta2)*velrhop2.w);
                const float vtemp=massp2*temp*(drx*frx+dry*fry+drz*frz);
                acep1.x+=vtemp*dvx; acep1.y+=vtemp*dvy; acep1.z+=vtemp*dvz;
              }
              //-SPS turbulence model.
              //-SPS turbulence model is disabled in v5.0 NN version
            }
          }
          rsym=(rsymp1&&!rsym && float(posp1.y-dry)<=KernelSize); //<vs_syymmetry>
          if(rsym)p2--;                                       //<vs_syymmetry>
        }
        else rsym=false; //<vs_syymmetry>
      }
    }
    //-Sum results together. | Almacena resultados.
    if(acep1.x||acep1.y||acep1.z) {
      ace[p1]=ace[p1]+acep1;
    }
  }
}

//==============================================================================
/// Perform interaction between particles for the SPH approcach, Pressure and vel gradients: Fluid/Float-Fluid/Float or Fluid/Float-Bound 
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
void JSphCpu::InteractionForcesFluid_NN_SPH_PressGrad(unsigned n,unsigned pinit,bool boundp2,float visco
  ,StDivDataCpu divdata,const unsigned *dcell
  ,tsymatrix3f* gradvel
  ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
  ,const float *press
  ,float &viscdt,float *ar,tfloat3 *ace,float *delta
  ,TpShifting shiftmode,tfloat4 *shiftposfs)const
{
  //-Initialize viscth to calculate viscdt maximo con OpenMP. | Inicializa viscth para calcular visdt maximo con OpenMP.
  float viscth[OMP_MAXTHREADS*OMP_STRIDE];
  for(int th=0; th<OmpThreads; th++)viscth[th*OMP_STRIDE]=0;
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int pfin=int(pinit+n);
#ifdef OMP_USE
#pragma omp parallel for schedule (guided)
#endif
  for(int p1=int(pinit); p1<pfin; p1++) {
    float visc=0,arp1=0,deltap1=0;
    tfloat3 acep1=TFloat3(0);
    tsymatrix3f gradvelp1={0,0,0,0,0,0};
    //-Variables for Shifting.
    tfloat4 shiftposfsp1;
    if(shift)shiftposfsp1=shiftposfs[p1];

    //-Obtain data of particle p1 in case of floating objects. | Obtiene datos de particula p1 en caso de existir floatings.
    bool ftp1=false;     //-Indicate if it is floating. | Indica si es floating.
    if(USE_FLOATING) {
      ftp1=CODE_IsFloating(code[p1]);
      if(ftp1 && tdensity!=DDT_None)deltap1=FLT_MAX; //-DDT is not applied to floating particles.
      if(ftp1 && shift)shiftposfsp1.x=FLT_MAX;  //-For floating objects do not calculate shifting. | Para floatings no se calcula shifting.
    }

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];
    const tfloat3 velp1=TFloat3(velrhop[p1].x,velrhop[p1].y,velrhop[p1].z);
    const float rhopp1=velrhop[p1].w;
    const float pressp1=press[p1];
    const bool rsymp1=(Symmetry && posp1.y<=KernelSize); //<vs_syymmetry>                                 
    //<vs_non-Newtonian>
    const typecode pp1=CODE_GetTypeValue(code[p1]);

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
    for(int z=ngs.zini; z<ngs.zfin; z++)for(int y=ngs.yini; y<ngs.yfin; y++) {
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

      //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
      //------------------------------------------------------------------------------------------------
      bool rsym=false; //<vs_syymmetry>
      for(unsigned p2=pif.x; p2<pif.y; p2++) {
        const float drx=float(posp1.x-pos[p2].x);
        float dry=float(posp1.y-pos[p2].y);
        if(rsym)    dry=float(posp1.y+pos[p2].y); //<vs_syymmetry>
        const float drz=float(posp1.z-pos[p2].z);
        const float rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<=KernelSize2 && rr2>=ALMOSTZERO) {
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(CSP,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients

          //===== Get mass of particle p2 ===== 
          //<vs_non-Newtonian>
          const typecode pp2=(boundp2 ? pp1 : CODE_GetTypeValue(code[p2])); //<vs_non-Newtonian>
          float massp2=(boundp2 ? MassBound : PhaseArray[pp2].mass); //-Contiene masa de particula segun sea bound o fluid.
          //Note if you masses are very different more than a ratio of 1.3 then: massp2 = (boundp2 ? PhaseArray[pp1].mass : PhaseArray[pp2].mass);

          //Floating                      
          bool ftp2=false;    //-Indicate if it is floating | Indica si es floating.
          bool compute=true;  //-Deactivate when using DEM and if it is of type float-float or float-bound | Se desactiva cuando se usa DEM y es float-float o float-bound.
          if(USE_FLOATING) {
            ftp2=CODE_IsFloating(code[p2]);
            if(ftp2)massp2=FtObjs[CODE_GetTypeValue(code[p2])].massp;
#ifdef DELTA_HEAVYFLOATING
            if(ftp2 && tdensity==DDT_DDT && massp2<=(MassFluid*1.2f))deltap1=FLT_MAX;
#else
            if(ftp2 && tdensity==DDT_DDT)deltap1=FLT_MAX;
#endif
            if(ftp2 && shift && shiftmode==SHIFT_NoBound)shiftposfsp1.x=FLT_MAX; //-With floating objects do not use shifting. | Con floatings anula shifting.
            compute=!(USE_FTEXTERNAL && ftp1&&(boundp2||ftp2)); //-Deactivate when using DEM and if it is of type float-float or float-bound. | Se desactiva cuando se usa DEM y es float-float o float-bound.
          }

          tfloat4 velrhop2=velrhop[p2];
          if(rsym)velrhop2.y=-velrhop2.y; //<vs_syymmetry>

          //===== Acceleration ===== 
          if(compute) {
            const float prs=(pressp1+press[p2])/(rhopp1*velrhop2.w)+(tker==KERNEL_Cubic ? fsph::GetKernelCubic_Tensil(CSP,rr2,rhopp1,pressp1,velrhop2.w,press[p2]) : 0);
            const float p_vpm=-prs*massp2;
            acep1.x+=p_vpm*frx; acep1.y+=p_vpm*fry; acep1.z+=p_vpm*frz;
          }

          //-Density derivative.
          const float rhop1over2=rhopp1/velrhop2.w;
          float dvx=velp1.x-velrhop2.x,dvy=velp1.y-velrhop2.y,dvz=velp1.z-velrhop2.z;
          if(compute)arp1+=massp2*(dvx*frx+dvy*fry+dvz*frz)*rhop1over2;

          const float cbar=max(PhaseArray[pp2].Cs0,PhaseArray[pp2].Cs0); //<vs_non-Newtonian>
          //-Density Diffusion Term (DeltaSPH Molteni).
          if(tdensity==DDT_DDT && deltap1!=FLT_MAX) {
            const float visc_densi=DDTkh*cbar*(rhop1over2-1.f)/(rr2+Eta2);
            const float dot3=(drx*frx+dry*fry+drz*frz);
            const float delta=(pp1==pp2 ? visc_densi*dot3*massp2 : 0); //<vs_non-Newtonian>
            //deltap1 = (boundp2 ? FLT_MAX : deltap1 + delta);
            deltap1=(boundp2 && TBoundary==BC_DBC ? FLT_MAX : deltap1+delta);
          }
          //-Density Diffusion Term (Fourtakas et al 2019).  //<vs_dtt2_ini>
          if((tdensity==DDT_DDT2||(tdensity==DDT_DDT2Full&&!boundp2))&&deltap1!=FLT_MAX&&!ftp2) {
            const float rh=1.f+DDTgz*drz;
            const float drhop=RhopZero*pow(rh,1.f/Gamma)-RhopZero;
            const float visc_densi=DDTkh*cbar*((velrhop2.w-rhopp1)-drhop)/(rr2+Eta2);
            const float dot3=(drx*frx+dry*fry+drz*frz);
            const float delta=(pp1==pp2 ? visc_densi*dot3*massp2/velrhop2.w : 0); //<vs_non-Newtonian>
            deltap1=(boundp2 ? FLT_MAX : deltap1-delta);
          }  //<vs_dtt2_end>

             //-Shifting correction.
          if(shift && shiftposfsp1.x!=FLT_MAX) {
            bool heavyphase=(PhaseArray[pp1].mass>PhaseArray[pp2].mass && pp1!=pp2 ? true : false);
            const float massrhop=massp2/velrhop2.w;
            const bool noshift=(boundp2&&(shiftmode==SHIFT_NoBound||(shiftmode==SHIFT_NoFixed && CODE_IsFixed(code[p2]))));
            shiftposfsp1.x=(noshift ? FLT_MAX : (heavyphase ? 0 : shiftposfsp1.x+massrhop*frx)); //-For boundary do not use shifting. | Con boundary anula shifting.                            
            shiftposfsp1.y+=(heavyphase ? 0 : massrhop*fry); //<vs_non-Newtonian>
            shiftposfsp1.z+=(heavyphase ? 0 : massrhop*frz); //<vs_non-Newtonian>
            shiftposfsp1.w-=(heavyphase ? 0 : massrhop*(drx*frx+dry*fry+drz*frz)); //<vs_non-Newtonian>
          }

          //===== Viscosity ===== 
          if(compute) {
            const float dot=drx*dvx+dry*dvy+drz*dvz;
            const float dot_rr2=dot/(rr2+Eta2);
            visc=max(dot_rr2,visc);
            if(tvisco!=VISCO_Artificial) { //<vs_non-Newtonian>
              {//vel gradients
                if(boundp2) {
                  dvx=2.f*velp1.x; dvy=2.f*velp1.y; dvz=2.f*velp1.z;  //fomraly I should use the moving BC vel as ug=2ub-uf
                }
                GetVelocityGradients_SPH_tsym(massp2,velrhop2,dvx,dvy,dvz,frx,fry,frz,gradvelp1);
              }
            }
          }
          rsym=(rsymp1&&!rsym && float(posp1.y-dry)<=KernelSize);   //<vs_syymmetry>
          if(rsym)p2--;																										//<vs_syymmetry>
        }
        else rsym=false;																									//<vs_syymmetry>
      }
    }
    //-Sum results together. | Almacena resultados.
    if(shift||arp1||acep1.x||acep1.y||acep1.z||visc) {
      if(tdensity!=DDT_None) {
        if(delta)delta[p1]=(delta[p1]==FLT_MAX||deltap1==FLT_MAX ? FLT_MAX : delta[p1]+deltap1);
        else if(deltap1!=FLT_MAX)arp1+=deltap1;
      }
      ar[p1]+=arp1;
      ace[p1]=ace[p1]+acep1;
      const int th=omp_get_thread_num();
      if(visc>viscth[th*OMP_STRIDE])viscth[th*OMP_STRIDE]=visc;
      if(tvisco!=VISCO_Artificial) {
        gradvel[p1].xx+=gradvelp1.xx;
        gradvel[p1].xy+=gradvelp1.xy;
        gradvel[p1].xz+=gradvelp1.xz;
        gradvel[p1].yy+=gradvelp1.yy;
        gradvel[p1].yz+=gradvelp1.yz;
        gradvel[p1].zz+=gradvelp1.zz;
      }
      if(shift)shiftposfs[p1]=shiftposfsp1;
    }
  }
  //-Keep max value in viscdt. | Guarda en viscdt el valor maximo.
  for(int th=0; th<OmpThreads; th++)if(viscdt<viscth[th*OMP_STRIDE])viscdt=viscth[th*OMP_STRIDE];
}

////==============================================================================
///// Computes sub-particle stress tensor (Tau) for SPS turbulence model.   
//// The SPS model is disabled in the v5.0 NN version
////==============================================================================
//void JSphCpu::ComputeSpsTau(unsigned n, unsigned pini, const tfloat4 *velrhop, const tsymatrix3f *gradvel, tsymatrix3f *tau)const {
//  const int pfin = int(pini + n);
//#ifdef OMP_USE
//#pragma omp parallel for schedule (static)
//#endif
//  for (int p = int(pini); p<pfin; p++) {
//      const tsymatrix3f gradvel = SpsGradvelc[p];
//      const float pow1 = gradvel.xx*gradvel.xx + gradvel.yy*gradvel.yy + gradvel.zz*gradvel.zz;
//      const float prr = pow1 + pow1 + gradvel.xy*gradvel.xy + gradvel.xz*gradvel.xz + gradvel.yz*gradvel.yz;
//      const float visc_sps = SpsSmag*sqrt(prr);
//      const float div_u = gradvel.xx + gradvel.yy + gradvel.zz;
//      const float sps_k = (2.0f / 3.0f)*visc_sps*div_u;
//      const float sps_blin = SpsBlin*prr;
//      const float sumsps = -(sps_k + sps_blin);
//      const float twovisc_sps = (visc_sps + visc_sps);
//      const float one_rho2 = 1.0f / velrhop[p].w;
//      tau[p].xx = one_rho2*(twovisc_sps*gradvel.xx + sumsps);
//      tau[p].xy = one_rho2*(visc_sps   *gradvel.xy);
//      tau[p].xz = one_rho2*(visc_sps   *gradvel.xz);
//      tau[p].yy = one_rho2*(twovisc_sps*gradvel.yy + sumsps);
//      tau[p].yz = one_rho2*(visc_sps   *gradvel.yz);
//      tau[p].zz = one_rho2*(twovisc_sps*gradvel.zz + sumsps);
//  }
//}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid (forces and DEM) for NN using the SPH approach
//==============================================================================
template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
void JSphCpu::Interaction_ForcesCpuT_NN_SPH(const stinterparmsc &t,StInterResultc &res)const
{
  float viscdt=res.viscdt;
  float viscetadt=res.viscetadt;
  if(t.npf) {
    //Pressure gradient, velocity gradients (symetric)
    //-Interaction Fluid-Fluid.
    InteractionForcesFluid_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift>(t.npf,t.npb,false,Visco
      ,t.divdata,t.dcell,t.spsgradvel,t.pos,t.velrhop,t.code,t.idp,t.press
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);
    //-Interaction Fluid-Bound.
    InteractionForcesFluid_NN_SPH_PressGrad<tker,ftmode,tvisco,tdensity,shift>(t.npf,t.npb,true,Visco*ViscoBoundFactor
      ,t.divdata,t.dcell,t.spsgradvel,t.pos,t.velrhop,t.code,t.idp,t.press
      ,viscdt,t.ar,t.ace,t.delta,t.shiftmode,t.shiftposfs);

    if(tvisco!=VISCO_Artificial) {
      //Build strain rate tensor and compute eta_visco //what will happen here for boundary particles?
      InteractionForcesFluid_NN_SPH_Visco_eta<ftmode,tvisco>(t.npf,t.npb,Visco
        ,t.visco_eta,t.velrhop,t.spsgradvel,t.d_tensor,t.auxnn,viscetadt,t.pos,t.code,t.idp);
    }

    if(tvisco!=VISCO_ConstEq) {
      //Morris
      InteractionForcesFluid_NN_SPH_Morris<tker,ftmode,tvisco,tdensity,shift>(t.npf,t.npb,false,Visco
        ,t.divdata,t.dcell,t.visco_eta,t.spstau,t.spsgradvel,t.auxnn,t.pos,t.velrhop,t.code,t.idp,t.ace);
      //-Interaction Fluid-Bound.     
      InteractionForcesFluid_NN_SPH_Morris<tker,ftmode,tvisco,tdensity,shift>(t.npf,t.npb,true,Visco*ViscoBoundFactor
        ,t.divdata,t.dcell,t.visco_eta,t.spstau,t.spsgradvel,t.auxnn,t.pos,t.velrhop,t.code,t.idp,t.ace);
    }
    else {
      //-ConsEq;
      // Build strain rate tensor and compute eta_visco //what will happen here for boundary particles?
      InteractionForcesFluid_NN_SPH_Visco_Stress_tensor<ftmode,tvisco>(t.npf,t.npb,Visco
        ,t.visco_eta,t.spstau,t.d_tensor,t.auxnn,t.pos,t.code,t.idp);
      //-Compute viscous terms and add them to the acceleration       
      InteractionForcesFluid_NN_SPH_ConsEq<tker,ftmode,tvisco,tdensity,shift>(t.npf,t.npb,false,Visco
        ,t.divdata,t.dcell,t.visco_eta,t.spstau,t.auxnn,t.pos,t.velrhop,t.code,t.idp,t.ace);
      //-Interaction Fluid-Bound.     
      InteractionForcesFluid_NN_SPH_ConsEq<tker,ftmode,tvisco,tdensity,shift>(t.npf,t.npb,true,Visco*ViscoBoundFactor
        ,t.divdata,t.dcell,t.visco_eta,t.spstau,t.auxnn,t.pos,t.velrhop,t.code,t.idp,t.ace);
    }

    //-Interaction of DEM Floating-Bound & Floating-Floating. //(DEM)
    if(UseDEM)InteractionForcesDEM(CaseNfloat,t.divdata,t.dcell
      ,FtRidp,DemData,t.pos,t.velrhop,t.code,t.idp,viscdt,t.ace);
  }
  if(t.npbok) {
    //-Interaction Bound-Fluid.
    InteractionForcesBound_NN_SPH<tker,ftmode>(t.npbok,0,t.divdata,t.dcell
      ,t.pos,t.velrhop,t.code,t.idp,viscdt,t.ar);
  }
  res.viscdt=viscdt;
  res.viscetadt=viscetadt;
}

//end_of_file