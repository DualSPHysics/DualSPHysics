//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpu_ker.h \brief Declares functions and CUDA kernels for the Particle Interaction and System Update.

#ifndef _JSphGpu_ker_
#define _JSphGpu_ker_

#include "Types.h"
#include "JSphTimersGpu.h"
#include <cuda_runtime_api.h>

class JLog2;
class JBlockSizeAuto;

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
  float delta2h;            ///<delta2h=DeltaSph*H*2
  float scell,dosh,dp;
  float cteb,gamma;
  float rhopzero;           ///<rhopzero=RhopZero
  float ovrhopzero;         ///<ovrhopzero=1/RhopZero
  float movlimit;
  unsigned periactive;
  double xperincx,xperincy,xperincz;
  double yperincx,yperincy,yperincz;
  double zperincx,zperincy,zperincz;
  double maprealposminx,maprealposminy,maprealposminz;
  double maprealsizex,maprealsizey,maprealsizez;
  //-Values depending on the assigned domain (can change). | Valores que dependen del dominio asignado (puden cambiar).
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

/// Implements a set of functions and CUDA kernels for the particle interaction and system update.
namespace cusph{

inline float3 Float3(const tfloat3& v){ float3 p={v.x,v.y,v.z}; return(p); }
inline float3 Float3(float x,float y,float z){ float3 p={x,y,z}; return(p); }
inline tfloat3 ToTFloat3(const float3& v){ return(TFloat3(v.x,v.y,v.z)); }
inline double3 Double3(const tdouble3& v){ double3 p={v.x,v.y,v.z}; return(p); }
inline double4 Double4(const tdouble4& v){ double4 p={v.x,v.y,v.z,v.w}; return(p); }

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

//-Kernels for preparing the force calculation for Pos-Single.
//-Kernels para preparar calculo de fuerzas con Pos-Single.
void PreInteractionSingle(unsigned np,const double2 *posxy,const double *posz
  ,const float4 *velrhop,float4 *pospress,float cteb,float ctegamma);

//-Kernels for the force calculation.
void Interaction_Forces(bool psingle,TpKernel tkernel,TpFtMode ftmode,bool lamsps
  ,TpDeltaSph tdelta,TpCellMode cellmode
  ,float viscob,float viscof,unsigned bsbound,unsigned bsfluid
  ,unsigned np,unsigned npb,unsigned npbok,tuint3 ncells
  ,const int2 *begincell,tuint3 cellmin,const unsigned *dcell
  ,const double2 *posxy,const double *posz,const float4 *pospress
  ,const float4 *velrhop,const typecode *code,const unsigned *idp
  ,const float *ftomassp,const tsymatrix3f *tau,tsymatrix3f *gradvel
  ,float *viscdt,float* ar,float3 *ace,float *delta
  ,TpShifting tshifting,float3 *shiftpos,float *shiftdetect
  ,bool simulate2d,StKerInfo *kerinfo,JBlockSizeAuto *bsauto);

//-Kernels for the calculation of the DEM forces.
void Interaction_ForcesDem(bool psingle,TpCellMode cellmode,unsigned bsize
  ,unsigned nfloat,tuint3 ncells,const int2 *begincell,tuint3 cellmin,const unsigned *dcell
  ,const unsigned *ftridp,const float4 *demdata,float dtforce
  ,const double2 *posxy,const double *posz,const float4 *pospress,const float4 *velrhop
  ,const typecode *code,const unsigned *idp,float *viscdt,float3 *ace,StKerInfo *kerinfo);

//-Kernels for calculating the Laminar+SPS viscosity.
void ComputeSpsTau(unsigned np,unsigned npb,float smag,float blin
  ,const float4 *velrhop,const tsymatrix3f *gradvelg,tsymatrix3f *tau);

//-Kernels for Delta-SPH.
void AddDelta(unsigned n,const float *delta,float *ar);

//-Kernels for Shifting.
void RunShifting(unsigned np,unsigned npb,double dt
  ,double shiftcoef,float shifttfs,double coeftfs
  ,const float4 *velrhop,const float *shiftdetect,float3 *shiftpos);

//-Kernels for ComputeStep (vel & rhop).
void ComputeStepVerlet(bool floatings,bool shift,unsigned np,unsigned npb
  ,const float4 *velrhop1,const float4 *velrhop2
  ,const float *ar,const float3 *ace,const float3 *shiftpos
  ,double dt,double dt2,float rhopoutmin,float rhopoutmax
  ,typecode *code,double2 *movxy,double *movz,float4 *velrhopnew);
void ComputeStepSymplecticPre(bool floatings,bool shift,unsigned np,unsigned npb
  ,const float4 *velrhoppre,const float *ar,const float3 *ace,const float3 *shiftpos
  ,double dtm,float rhopoutmin,float rhopoutmax
  ,typecode *code,double2 *movxy,double *movz,float4 *velrhop);
void ComputeStepSymplecticCor(bool floatings,bool shift,unsigned np,unsigned npb
  ,const float4 *velrhoppre,const float *ar,const float3 *ace,const float3 *shiftpos
  ,double dtm,double dt,float rhopoutmin,float rhopoutmax
  ,typecode *code,double2 *movxy,double *movz,float4 *velrhop);

//-Kernels for ComputeStep (position).
void ComputeStepPos (byte periactive,bool floatings,unsigned np,unsigned npb,const double2 *movxy,const double *movz,double2 *posxy,double *posz,unsigned *dcell,typecode *code);
void ComputeStepPos2(byte periactive,bool floatings,unsigned np,unsigned npb,const double2 *posxypre,const double *poszpre,const double2 *movxy,const double *movz,double2 *posxy,double *posz,unsigned *dcell,typecode *code);

//-Kernels for Motion.
void CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini,unsigned idfin,const typecode *code,const unsigned *idp,unsigned *ridp);
void MoveLinBound(byte periactive,unsigned np,unsigned ini,tdouble3 mvpos,tfloat3 mvvel,const unsigned *ridp,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);
void MoveMatBound(byte periactive,bool simulate2d,unsigned np,unsigned ini,tmatrix4d m,double dt,const unsigned *ridpmv,double2 *posxy,double *posz,unsigned *dcell,float4 *velrhop,typecode *code);

//-Kernels for Floating bodies.
void FtCalcForcesSum(bool periactive,unsigned ftcount
  ,tfloat3 gravity,const float4 *ftodata,const double3 *ftocenter,const unsigned *ftridp
  ,const double2 *posxy,const double *posz,const float3 *ace
  ,float3 *ftoforcessum);
void FtCalcForces(unsigned ftcount,tfloat3 gravity,const float4 *ftodata
  ,const float3 *ftoangles,const float4 *ftoinertiaini8,const float *ftoinertiaini1
  ,const float3 *ftoforcessum,float3 *ftoforces);
void FtCalcForcesRes(unsigned ftcount,bool simulate2d,double dt
  ,const float3 *ftoomega,const float3 *ftovel,const double3 *ftocenter,const float3 *ftoforces
  ,float3 *ftoforcesres,double3 *ftocenterres);
void FtUpdate(bool periactive,bool predictor,unsigned ftcount,double dt
  ,const float4 *ftodata,const float3 *ftoforcesres,double3 *ftocenterres,const unsigned *ftridp
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

//-Kernels for external forces (JSphAccInput).
void AddAccInput(unsigned n,unsigned pini,typecode codesel
  ,tdouble3 acclin,tdouble3 accang,tdouble3 centre,tdouble3 velang,tdouble3 vellin,bool setgravity
  ,tfloat3 gravity,const typecode *code,const double2 *posxy,const double *posz,const float4 *velrhop,float3 *ace);

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


