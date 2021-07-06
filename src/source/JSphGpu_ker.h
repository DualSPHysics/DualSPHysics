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
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>

class JLog2;

#define SPHBSIZE 256

/// Structure with constants stored in the constant memory of GPU for the particle interactions.
typedef struct{
  unsigned nbound;
  float massb;              ///<Reference mass of the general boundary particle [kg].
  float massf;              ///<Reference mass of the fluid particle [kg].
  float kernelh;            ///<The smoothing length of SPH kernel [m].
  float kernelsize2;        ///<Maximum interaction distance squared (KernelSize^2).
  float poscellsize;        ///<Size of cells used for coding PosCell (it is usually KernelSize).
  float awen;               ///<Wendland kernel constant (awen) to compute wab.
  float bwen;               ///<Wendland kernel constant (bwen) to compute fac (kernel derivative).
  float cs0;                ///<Speed of sound at the reference density.
  float eta2;               ///<Constant related to H (Eta2=(h*0.1)*(h*0.1)).
  float ddtkh;              ///<Constant for DDT1 & DDT2. DDTkh=DDTValue*KernelSize
  float ddtgz;              ///<Constant for DDT2.        ddtgz=RhopZero*Gravity.z/CteB
  float scell;              ///<Cell size: KernelSize/ScellDiv (KernelSize or KernelSize/2).
  float kernelsize;         ///<Maximum interaction distance between particles (KernelK*KernelH).
  float dp;                 ///<Initial distance between particles [m].
  float cteb;               ///<Constant used in the state equation [Pa].
  float gamma;              ///<Politropic constant for water used in the state equation.
  float rhopzero;           ///<Reference density of the fluid [kg/m3].
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
  float viscob,viscof;
  unsigned bsbound,bsfluid;
  unsigned vnp,vnpb,vnpbok;
  unsigned boundini;
  unsigned fluidini;
  unsigned boundnum;
  unsigned fluidnum;         
  unsigned id;
  StDivDataGpu divdatag;
  //-Input data arrays.
  const unsigned *dcell;
  const double2 *posxy;
  const double *posz;
  const float4 *poscell;
  const float4 *velrhop;
  const unsigned *idp;
  const typecode *code;
  const float *ftomassp;
  const tsymatrix3f *tau;
  const float3 *dengradcorr;
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
    ,float viscob_,float viscof_
    ,unsigned bsbound_,unsigned bsfluid_
    ,unsigned np_,unsigned npb_,unsigned npbok_
    ,unsigned id_
    ,const StDivDataGpu &divdatag_,const unsigned *dcell_
    ,const double2 *posxy_,const double *posz_,const float4 *poscell_
    ,const float4 *velrhop_,const unsigned *idp_,const typecode *code_
    ,const float *ftomassp_,const tsymatrix3f *spstau_
    ,const float3 *dengradcorr_
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
    viscob=viscob_; viscof=viscof_;
    bsbound=bsbound_; bsfluid=bsfluid_;
    vnp=np_; vnpb=npb_; vnpbok=npbok_;
    boundini=0;   boundnum=vnpbok;
    fluidini=vnpb; fluidnum=vnp-vnpb;
    id=id_; 
    divdatag=divdatag_;
    //-Input data arrays.
    dcell=dcell_;
    posxy=posxy_; posz=posz_; poscell=poscell_;
    velrhop=velrhop_; idp=idp_; code=code_;
    ftomassp=ftomassp_; tau=spstau_;
    dengradcorr=dengradcorr_;
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

//-Kernels for the boundary correction (mDBC).
void Interaction_MdbcCorrection(TpKernel tkernel,bool simulate2d
  ,TpSlipMode slipmode,bool fastsingle,unsigned n,unsigned nbound
  ,float mdbcthreshold,const StDivDataGpu &dvd,const tdouble3 &mapposmin
  ,const double2 *posxy,const double *posz,const float4 *poscell
  ,const typecode *code,const unsigned *idp,const float3 *boundnormal
  ,const float3 *motionvel,float4 *velrhop);

//-Kernels for the calculation of the DEM forces.
void Interaction_ForcesDem(unsigned bsize,unsigned nfloat
  ,const StDivDataGpu &dvd,const unsigned *dcell
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
void CopyMotionVel(unsigned nmoving,const unsigned *ridpmv,const float4 *velrhop,float3 *motionvel);

//-Kernels for MLPistons motion.
void MovePiston1d(bool periactive,unsigned np,unsigned idini,double dp,double poszmin,unsigned poszcount,const byte *pistonid,const double* movx,const double* velx,const unsigned *ridpmv,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);
void MovePiston2d(bool periactive,unsigned np,unsigned idini,double dp,double posymin,double poszmin,unsigned poszcount,const double* movx,const double* velx,const unsigned *ridpmv,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);

//-Kernels for Floating bodies.
void FtCalcForcesSum(bool periactive,unsigned ftcount
  ,const float4 *ftodata,const double3 *ftocenter,const unsigned *ftridp
  ,const double2 *posxy,const double *posz,const float3 *ace
  ,float3 *ftoforcessum);
void FtCalcForces(unsigned ftcount,tfloat3 gravity
  ,const float *ftomass,const float3 *ftoangles
  ,const float4 *ftoinertiaini8,const float *ftoinertiaini1
  ,const float3 *ftoforcessum,float3 *ftoforces,const float3 *ftoextforces);
void FtCalcForcesRes(unsigned ftcount,bool simulate2d,double dt
  ,const float3 *ftovelace,const double3 *ftocenter,const float3 *ftoforces
  ,float3 *ftoforcesres,double3 *ftocenterres);
void FtApplyConstraints(unsigned ftcount,const byte *ftoconstraints
  ,float3 *ftoforces,float3 *ftoforcesres);
void FtUpdate(bool periactive,bool predictor,unsigned ftcount,double dt
  ,const float4 *ftodatp,const float3 *ftoforcesres,double3 *ftocenterres,const unsigned *ftridp
  ,double3 *ftocenter,float3 *ftoangles,float3 *ftovelace
  ,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);
void FtGetPosRef(unsigned np,const unsigned *idpref,const unsigned *ftridp //<vs_ftmottionsv>
  ,const double2 *posxy,const double *posz,double *posref);                 //<vs_ftmottionsv>

//-Kernels for periodic conditions.
void PeriodicIgnore(unsigned n,typecode *code);
unsigned PeriodicMakeList(unsigned n,unsigned pini,bool stable,unsigned nmax,tdouble3 mapposmin,tdouble3 mapposmax,tdouble3 perinc,const double2 *posxy,const double *posz,const typecode *code,unsigned *listp);
void PeriodicDuplicateVerlet(unsigned n,unsigned pini,tuint3 domcells,tdouble3 perinc
  ,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,tsymatrix3f *spstau,float4 *velrhopm1);
void PeriodicDuplicateSymplectic(unsigned n,unsigned pini
  ,tuint3 domcells,tdouble3 perinc,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
  ,double2 *posxy,double *posz,float4 *velrhop,tsymatrix3f *spstau,double2 *posxypre,double *poszpre,float4 *velrhoppre);
void PeriodicDuplicateNormals(unsigned n,unsigned pini,const unsigned *listp,float3 *normals,float3 *motionvel);

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


