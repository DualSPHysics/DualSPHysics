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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

/// \file JSphGpu_ker.h \brief Declares functions and CUDA kernels for the Particle Interaction and System Update.

#ifndef _JSphGpu_ker_
#define _JSphGpu_ker_

#include "DualSphDef.h"
#include "JSphTimersGpu.h"
#include <cuda_runtime_api.h>

class JLog2;

#define SPHBSIZE 256

/// Structure with constants stored in the constant memory of GPU for the particle interactions.
typedef struct{
  unsigned nbound;
  float massb;              ///<Mass of a boundary particle.
  float massf;              ///<Mass of a fluid particle.
  float h;                  ///<Smoothing length (=coef*sqrt(dx*dx+dy*dy+dz*dz))
  float fourh2;             ///< \ref h * \ref h * 4 
  float awen;               ///<Cte. of Wendland kernel to compute wab.
  float bwen;               ///<Cte. of Wendland kernel to compute fac (kernel derivative).
  float cs0;                ///<Speed of sound of reference.
  float eta2;               ///<eta*eta being eta=0.1*\ref h
  float ddt2h;              ///<Constant for DDT1 & DDT2. ddt2h=DDTValue*2*H
  float ddtgz;              ///<Constant for DDT2.        ddtgz=RhopZero*Gravity.z/CteB
  float scell,dosh,dp;
  float cteb,gamma;
  float rhopzero;           ///<rhopzero=RhopZero
  float ovrhopzero;         ///<ovrhopzero=1/RhopZero
  float movlimit;
  unsigned symmetry;   //<vs_syymmetry>
  unsigned tboundary;  
  unsigned periactive;
  double xperincx,xperincy,xperincz;
  double yperincx,yperincy,yperincz;
  double zperincx,zperincy,zperincz;
  double maprealposminx,maprealposminy,maprealposminz;
  double maprealsizex,maprealsizey,maprealsizez;
  //-Values depending on the assigned domain (can change). | Valores que dependen del dominio asignado (puden cambiar).
  unsigned axis;
  unsigned cellcode;
  double domposminx,domposminy,domposminz;
  //-Ctes. of Gaussian kernel.
  float agau,bgau;
  //-Ctes. of Cubic Spline kernel.
  float cubic_a1,cubic_a2,cubic_aa,cubic_a24,cubic_c1,cubic_d1,cubic_c2,cubic_odwdeltap;
}StCteInteraction; 

/// Structure to collect kernel information.
typedef struct{
  int forcesbound_rg;
  int forcesbound_bs;
  int forcesbound_bsmax;
  int forcesfluid_rg;
  int forcesfluid_bs;
  int forcesfluid_bsmax;
  int forcesdem_rg;
  int forcesdem_bs;
  int forcesdem_bsmax;
}StKerInfo; 

///Returns the reordered value according to a given axis.
/// (x,y,z) in MGDIV_Z ---> (y,z,x) for MGDIV_X
/// (x,y,z) in MGDIV_Z ---> (x,z,y) for MGDIV_Y
inline tuint3 CodeAxisOrder(TpMgDivMode axis,tuint3 v){
  if(axis==MGDIV_X)return(TUint3(v.y,v.z,v.x));
  if(axis==MGDIV_Y)return(TUint3(v.x,v.z,v.y));
  return(v);
}

///Structure with the parameters for particle interaction on GPU.
typedef struct StrInterParmsg{
  //-Configuration options.
  bool simulate2d;
  bool symmetry; //<vs_syymmetry>
  TpKernel tkernel;
  TpFtMode ftmode;
  bool lamsps;
  TpDensity tdensity;
  TpShifting shiftmode;
  //-Execution values.
  int hdiv;  //hdiv=(cellmode==CELLMODE_H? 2: 1)
  float viscob,viscof;
  unsigned bsbound,bsfluid;
  unsigned vnp,vnpb,vnpbok;
  unsigned boundini;
  unsigned fluidini;
  unsigned boundnum;
  unsigned fluidnum;         
  unsigned id;
  TpMgDivMode axis;       ///<Axis used in current division. It is used to sort cells and particles.
  tuint3 ncells;
  int4 nc;                ///<Number of cells according axis.
  unsigned cellfluid;
  tint3 cellmin;
  //-Input data arrays.
  const int2 *begincell;
  const unsigned *dcell;
  const double2 *posxy;
  const double *posz;
  const float4 *poscell;
  const float4 *velrhop;
  const unsigned *idp;
  const typecode *code;
  const float *ftomassp;
  const tsymatrix3f *tau;
  //-Output data arrays.
  float *viscdt;
  float* ar;
  float3 *ace;
  float *delta;
  tsymatrix3f *gradvel;
  float4 *shiftposfs;
  //-Other values and objects.
  cudaStream_t stm;
  StKerInfo *kerinfo;

  ///Structure constructor.
  StrInterParmsg(
     bool simulate2d_
    ,bool symmetry_ //<vs_syymmetry>
    ,TpKernel tkernel_,TpFtMode ftmode_
    ,bool lamsps_,TpDensity tdensity_,TpShifting shiftmode_
    ,TpCellMode cellmode_
    ,float viscob_,float viscof_
    ,unsigned bsbound_,unsigned bsfluid_
    ,unsigned np_,unsigned npb_,unsigned npbok_
    ,unsigned id_,TpMgDivMode axis_
    ,tuint3 ncells_,tuint3 cellmin_
    ,const int2 *begincell_,const unsigned *dcell_
    ,const double2 *posxy_,const double *posz_,const float4 *poscell_
    ,const float4 *velrhop_,const unsigned *idp_,const typecode *code_
    ,const float *ftomassp_,const tsymatrix3f *spstau_
    ,float *viscdt_,float* ar_,float3 *ace_,float *delta_
    ,tsymatrix3f *spsgradvel_
    ,float4 *shiftposfs_
    ,cudaStream_t stm_
    ,StKerInfo *kerinfo_)
  {
    //-Configuration options.
    simulate2d=simulate2d_;
    symmetry=symmetry_; //<vs_syymmetry>
    tkernel=tkernel_; ftmode=ftmode_;
    lamsps=lamsps_; tdensity=tdensity_; shiftmode=shiftmode_;
    //-Execution values.
    hdiv=(cellmode_==CELLMODE_H? 2: 1);
    viscob=viscob_; viscof=viscof_;
    bsbound=bsbound_; bsfluid=bsfluid_;
    vnp=np_; vnpb=npb_; vnpbok=npbok_;
    boundini=0;   boundnum=vnpbok;
    fluidini=vnpb; fluidnum=vnp-vnpb;
    id=id_; axis=axis_;
    ncells=ncells_;
    const tuint3 nc3=CodeAxisOrder(axis,ncells);
    nc.x=int(nc3.x); nc.y=int(nc3.y); nc.z=int(nc3.z); nc.w=int(nc3.x*nc3.y);
    cellfluid=nc.w*nc.z+1;
    cellmin=TInt3(int(cellmin_.x),int(cellmin_.y),int(cellmin_.z));
    //-Input data arrays.
    begincell=begincell_; dcell=dcell_;
    posxy=posxy_; posz=posz_; poscell=poscell_;
    velrhop=velrhop_; idp=idp_; code=code_;
    ftomassp=ftomassp_; tau=spstau_;
    //-Output data arrays.
    viscdt=viscdt_; ar=ar_; ace=ace_; delta=delta_;
    gradvel=spsgradvel_;
    shiftposfs=shiftposfs_;
    //-Other values and objects.
    stm=stm_;
    kerinfo=kerinfo_;
  }

}StInterParmsg;


/// Implements a set of functions and CUDA kernels for the particle interaction and system update.
namespace cusph{

dim3 GetGridSize(unsigned n,unsigned blocksize);
inline unsigned ReduMaxFloatSize(unsigned ndata){ return((ndata/SPHBSIZE+1)+(ndata/(SPHBSIZE*SPHBSIZE)+SPHBSIZE)); }
float ReduMaxFloat(unsigned ndata,unsigned inidata,float* data,float* resu);
float ReduMaxFloat_w(unsigned ndata,unsigned inidata,float4* data,float* resu);

void CteInteractionUp(const StCteInteraction *cte);
void InitArray(unsigned n,float3 *v,tfloat3 value);
void Resety(unsigned n,unsigned ini,float3 *v);
void ComputeAceMod(unsigned n,const float3 *ace,float *acemod);
void ComputeAceMod(unsigned n,const typecode *code,const float3 *ace,float *acemod);

void ComputeVelMod(unsigned n,const float4 *vel,float *velmod);

//-Kernels for the force calculation.
void Interaction_Forces(const StInterParmsg &t);

//-Kernels for the calculation of the DEM forces.
void Interaction_ForcesDem(TpCellMode cellmode,unsigned bsize
  ,unsigned nfloat,tuint3 ncells,const int2 *begincell,tuint3 cellmin,const unsigned *dcell
  ,const unsigned *ftridp,const float4 *demdata,const float *ftomassp,float dtforce
  ,const float4 *poscell,const float4 *velrhop
  ,const typecode *code,const unsigned *idp,float *viscdt,float3 *ace,StKerInfo *kerinfo);

//-Kernels for calculating the Laminar+SPS viscosity.
void ComputeSpsTau(unsigned np,unsigned npb,float smag,float blin
  ,const float4 *velrhop,const tsymatrix3f *gradvelg,tsymatrix3f *tau,cudaStream_t stm=NULL);

//-Kernels for Delta-SPH.
void AddDelta(unsigned n,const float *delta,float *ar,cudaStream_t stm=NULL);

//-Kernels for ComputeStep (position).
void ComputeStepPos (byte periactive,bool floatings,unsigned np,unsigned npb
  ,const double2 *movxy,const double *movz,double2 *posxy,double *posz
  ,unsigned *dcell,typecode *code);
void ComputeStepPos2(byte periactive,bool floatings,unsigned np,unsigned npb
  ,const double2 *posxypre,const double *poszpre,const double2 *movxy,const double *movz
  ,double2 *posxy,double *posz,unsigned *dcell,typecode *code);

//-Kernels for Motion.
void CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini,unsigned idfin
  ,const typecode *code,const unsigned *idp,unsigned *ridp);
void MoveLinBound(byte periactive,unsigned np,unsigned ini,tdouble3 mvpos,tfloat3 mvvel
  ,const unsigned *ridp,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);
void MoveMatBound(byte periactive,bool simulate2d,unsigned np,unsigned ini,tmatrix4d m,double dt
  ,const unsigned *ridpmv,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code,float3 *boundnormal);

//-Kernels for MLPistons motion.  //<vs_mlapiston_ini>
void MovePiston1d(bool periactive,unsigned np,unsigned idini,double dp,double poszmin,unsigned poszcount,const byte *pistonid,const double* movx,const double* velx,const unsigned *ridpmv,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);
void MovePiston2d(bool periactive,unsigned np,unsigned idini,double dp,double posymin,double poszmin,unsigned poszcount,const double* movx,const double* velx,const unsigned *ridpmv,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);
//<vs_mlapiston_end>

//-Kernels for Floating bodies.
void FtCalcForcesSum(bool periactive,unsigned ftcount
  ,const float4 *ftodata,const double3 *ftocenter,const unsigned *ftridp
  ,const double2 *posxy,const double *posz,const float3 *ace
  ,float3 *ftoforcessum);
void FtCalcForces(unsigned ftcount,tfloat3 gravity
  ,const float *ftomass,const float3 *ftoangles
  ,const float4 *ftoinertiaini8,const float *ftoinertiaini1
  ,const float3 *ftoforcessum,float3 *ftoforces);
void FtCalcForcesRes(unsigned ftcount,bool simulate2d,double dt
  ,const float3 *ftoomega,const float3 *ftovel,const double3 *ftocenter,const float3 *ftoforces
  ,float3 *ftoforcesres,double3 *ftocenterres);
void FtApplyConstraints(unsigned ftcount,const byte *ftoconstraints
  ,float3 *ftoforces,float3 *ftoforcesres);
void FtUpdate(bool periactive,bool predictor,unsigned ftcount,double dt
  ,const float4 *ftodatp,const float3 *ftoforcesres,double3 *ftocenterres,const unsigned *ftridp
  ,double3 *ftocenter,float3 *ftoangles,float3 *ftovel,float3 *ftoomega
  ,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);

//-Kernels for periodic conditions.
void PeriodicIgnore(unsigned n,typecode *code);
unsigned PeriodicMakeList(unsigned n,unsigned pini,bool stable,unsigned nmax,tdouble3 mapposmin,tdouble3 mapposmax,tdouble3 perinc,const double2 *posxy,const double *posz,const typecode *code,unsigned *listp);
void PeriodicDuplicateVerlet(unsigned n,unsigned pini,tuint3 domcells,tdouble3 perinc
  ,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,tsymatrix3f *spstau,float4 *velrhopm1);
void PeriodicDuplicateSymplectic(unsigned n,unsigned pini
  ,tuint3 domcells,tdouble3 perinc,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,tsymatrix3f *spstau,double2 *posxypre,double *poszpre,float4 *velrhoppre);

//-Kernels for Damping.
void ComputeDamping(double dt,tdouble4 plane,float dist,float over,tfloat3 factorxyz,float redumax
  ,unsigned n,unsigned pini,const double2 *posxy,const double *posz,const typecode *code
  ,float4 *velrhop);
void ComputeDampingPla(double dt,tdouble4 plane,float dist,float over,tfloat3 factorxyz,float redumax
  ,double zmin,double zmax,tdouble4 pla0,tdouble4 pla1,tdouble4 pla2,tdouble4 pla3
  ,unsigned n,unsigned pini,const double2 *posxy,const double *posz,const typecode *code
  ,float4 *velrhop);

}


#endif


