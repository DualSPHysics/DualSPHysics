//HEAD_DSPH
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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

/// \file JSphGpu_ker.h \brief Declares functions and CUDA kernels for the Particle Interaction and System Update.

#ifndef _JSphGpu_ker_
#define _JSphGpu_ker_

#include "DualSphDef.h"
#include "JSphGpu_cte.h"
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>
#include "JSphVRes.h"       //<vs_vrres

class JLog2;

#define SPHBSIZE 256


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
  TpKernel tkernel;
  TpFtMode ftmode;
  TpVisco tvisco;
  TpDensity tdensity;
  TpShifting shiftmode;
  TpMdbc2Mode mdbc2;      //<vs_m2dbcNP>
  bool shiftadv;   //<vs_advshift>
  bool corrector;  //<vs_advshift>
  bool aleform;    //<vs_advshift>
  bool ncpress;    //<vs_advshift>
  //-Execution values.
  float viscob,viscof;
  unsigned bsbound,bsfluid;
  unsigned vnp,vnpb,vnpbok;
  unsigned boundini;
  unsigned fluidini;
  unsigned boundnum;
  unsigned fluidnum;         
  unsigned id;
  unsigned nstep;
  StDivDataGpu divdatag;
  //-Input data arrays.
  const unsigned*  dcell;
  const double2*   posxy;
  const double*    posz;
  const float4*    poscell;
  const float4*    velrho;
  const unsigned*  idp;
  const typecode*  code;
  const byte*      boundmode;    //<vs_m2dbc>
  const float3*    tangenvel;    //<vs_m2dbc>
  const float3*    motionvel;    //<vs_m2dbc>
  const float3*    boundnormal;  //<vs_m2dbcNP> 
  float4* nopenshift;            //<vs_m2dbcNP>  
  const float*     ftomassp;
  const tsymatrix3f* spstaurho2;
  const float3*    dengradcorr;
  //-Output data arrays.
  float*  viscdt;
  float*  ar;
  float3* ace;
  float*  delta;
  tsymatrix3f*  sps2strain;
  float4*       shiftposfs;
  unsigned*     fstype;          //<vs_advshift>
  const float4* shiftvel;        //<vs_advshift>
  const float*  psiclean;        //<vs_divclean> 
  float*        psicleanrhs;     //<vs_divclean>
  float*        cspsiclean;      //<vs_divclean>
  float         divcleankp;      //<vs_divclean>
  bool          divclean;        //<vs_divclean> 
  //-Other values and objects.
  cudaStream_t stm;
  StKerInfo* kerinfo;

  ///Structure constructor.
  StrInterParmsg(
     bool simulate2d
    ,TpKernel tkernel
    ,TpFtMode ftmode
    ,TpVisco tvisco
    ,TpDensity tdensity
    ,TpShifting shiftmode
    ,TpMdbc2Mode mdbc2      //<vs_m2dbcNP>
    ,bool shiftadv      //<vs_advshift>
    ,bool corrector     //<vs_advshift>
    ,bool aleform       //<vs_advshift>
    ,bool ncpress       //<vs_advshift>
    ,float viscob,float viscof
    ,unsigned bsbound,unsigned bsfluid
    ,unsigned np,unsigned npb,unsigned npbok
    ,unsigned id
    ,unsigned nstep
    ,const StDivDataGpu& divdatag
    ,const unsigned* dcell
    ,const double2* posxy,const double* posz,const float4* poscell
    ,const float4* velrho,const unsigned* idp,const typecode* code
    ,const byte*   boundmode    //<vs_m2dbc>
    ,const float3* tangenvel    //<vs_m2dbc>
    ,const float3* motionvel    //<vs_m2dbc>
    ,const float3* boundnormal,float4* nopenshift  //<vs_m2dbcNP> 
    ,const float* ftomassp
    ,const tsymatrix3f* spstaurho2
    ,const float3* dengradcorr
    ,float* viscdt
    ,float* ar
    ,float3* ace
    ,float* delta
    ,tsymatrix3f* sps2strain
    ,float4* shiftposfs
    ,unsigned* fstype           //<vs_advshift>
    ,const float4* shiftvel     //<vs_advshift>
    ,const float* psiclean,float* psicleanrhs,float* cspsiclean   //<vs_divclean>
    ,float divcleankp,bool divclean                               //<vs_divclean>
    ,cudaStream_t stm
    ,StKerInfo* kerinfo)
  {
    //-Configuration options.
    this->simulate2d=simulate2d;
    this->tkernel=tkernel; 
    this->ftmode=ftmode;
    this->tvisco=tvisco;
    this->tdensity=tdensity;
    this->shiftmode=shiftmode;
    this->mdbc2=mdbc2;         //<vs_m2dbcNP>
    this->shiftadv=shiftadv;   //<vs_advshift>
    this->corrector=corrector; //<vs_advshift>
    this->aleform=aleform;     //<vs_advshift>
    this->ncpress=ncpress;     //<vs_advshift>
    //-Execution values.
    this->viscob=viscob;   this->viscof=viscof;
    this->bsbound=bsbound; this->bsfluid=bsfluid;
    this->vnp=np; this->vnpb=npb; this->vnpbok=npbok;
    this->boundini=0;    this->boundnum=vnpbok;
    this->fluidini=vnpb; this->fluidnum=vnp-vnpb;
    this->id=id;
    this->nstep=nstep; 
    this->divdatag=divdatag;
    //-Input data arrays.
    this->dcell=dcell;
    this->posxy=posxy; this->posz=posz; this->poscell=poscell;
    this->velrho=velrho; this->idp=idp; this->code=code;
    this->boundmode=boundmode;    //<vs_m2dbc>
    this->tangenvel=tangenvel;    //<vs_m2dbc>
    this->motionvel=motionvel;    //<vs_m2dbc>
    this->boundnormal=boundnormal; this->nopenshift=nopenshift; //<vs_m2dbcNP> 
    this->ftomassp=ftomassp;
    this->spstaurho2=spstaurho2;
    this->dengradcorr=dengradcorr;
    //-Output data arrays.
    this->viscdt=viscdt;
    this->ar=ar;
    this->ace=ace;
    this->delta=delta;
    this->sps2strain=sps2strain;
    this->shiftposfs=shiftposfs;
    this->fstype=fstype;      //<vs_advshift>
    this->shiftvel=shiftvel;  //<vs_advshift>
    this->psiclean=psiclean;  this->psicleanrhs=psicleanrhs; this->cspsiclean=cspsiclean; //<vs_divclean>
    this->divcleankp=divcleankp; this->divclean=divclean;     //<vs_divclean>
    //-Other values and objects.
    this->stm=stm;
    this->kerinfo=kerinfo;
  }

}StInterParmsg;

//<vs_flexstruc_ini>
///Structure with the parameters for flexible structure interaction on GPU.
typedef struct StrInterParmsFlexStrucg{
  //-Configuration options.
  bool simulate2d;
  TpKernel tkernel;
  TpVisco tvisco;
  TpMdbc2Mode mdbc2;
  //-Execution values.
  float viscob;
  unsigned vnpfs;
  StDivDataGpu divdatag;
  //-Input data arrays.
  const unsigned* dcell;
  const float4* poscell;
  const float4* velrhop;
  const typecode* code;
  const byte* boundmode;
  const float3* tangenvel;
  const StFlexStrucData* flexstrucdata;
  const unsigned* flexstrucridp;
  const float4* poscell0;
  const unsigned* numpairs;
  const unsigned* const* pairidx;
  const tmatrix3f* kercorr;
  const float3* boundnor0;
  //-Output data arrays.
  tmatrix3f* defgrad;
  float3* boundnor;
  float* flexstrucdt;
  float3* ace;
  //-Other values and objects.
  cudaStream_t stm;

  ///Structure constructor.
  StrInterParmsFlexStrucg(
       bool simulate2d,TpKernel tkernel,TpVisco tvisco,TpMdbc2Mode mdbc2
      ,float viscob,unsigned vnpfs
      ,const StDivDataGpu& divdatag
      ,const unsigned* dcell
      ,const float4* poscell,const float4* velrhop,const typecode* code
      ,const byte* boundmode,const float3* tangenvel
      ,const StFlexStrucData* flexstrucdata
      ,const unsigned* flexstrucridp,const float4* poscell0
      ,const unsigned* numpairs,const unsigned* const* pairidx
      ,const tmatrix3f* kercorr,const float3* boundnor0
      ,tmatrix3f* defgrad,float3* boundnor,float* flexstrucdt,float3* ace
      ,cudaStream_t stm)
  {
    //-Configuration options.
    this->simulate2d=simulate2d; this->tkernel=tkernel; this->tvisco=tvisco; this->mdbc2=mdbc2;
    //-Execution values.
    this->viscob=viscob;
    this->vnpfs=vnpfs;
    this->divdatag=divdatag;
    //-Input data arrays.
    this->dcell=dcell;
    this->poscell=poscell; this->velrhop=velrhop; this->code=code;
    this->boundmode=boundmode; this->tangenvel=tangenvel;
    this->flexstrucdata=flexstrucdata;
    this->flexstrucridp=flexstrucridp; this->poscell0=poscell0;
    this->numpairs=numpairs; this->pairidx=pairidx;
    this->kercorr=kercorr; this->boundnor0=boundnor0;
    //-Output data arrays.
    this->defgrad=defgrad; this->boundnor=boundnor; this->flexstrucdt=flexstrucdt; this->ace=ace;
    //-Other values and objects.
    this->stm=stm;
  }
}StInterParmsFlexStrucg;
//<vs_flexstruc_end>


/// Implements a set of functions and CUDA kernels for the particle interaction and system update.
namespace cusph{

inline unsigned ReduMaxFloatSize(unsigned ndata){ return((ndata/SPHBSIZE+1)+(ndata/(SPHBSIZE*SPHBSIZE)+SPHBSIZE)); }
float ReduMaxFloat(unsigned ndata,unsigned inidata,float* data,float* resu);
float ReduMaxFloat_w(unsigned ndata,unsigned inidata,float4* data,float* resu);

void CteInteractionUp(const StCteInteraction* cte);
void InitArray(unsigned n,float3* v,tfloat3 value);
void Resety(unsigned n,unsigned ini,float3* v);
void ComputeAceMod(unsigned n,const float3* ace,float* acemod);
void ComputeAceMod(unsigned n,const typecode* code,const float3* ace,float* acemod);

void ComputeVelMod(unsigned n,const float4* vel,float* velmod);

//-Kernels for the force calculation.
void Interaction_Forces(const StInterParmsg& t);

//-Kernels for the calculation of the DEM forces.
void Interaction_ForcesDem(unsigned bsize,unsigned nfloat
  ,const StDivDataGpu& dvd,const unsigned* dcell
  ,const unsigned* ftridp,const float4* demdata,const float* ftomassp
  ,float dtforce,const float4* poscell,const float4* velrho
  ,const typecode* code,const unsigned* idp,float* viscdt,float3* ace
  ,StKerInfo* kerinfo,cudaStream_t stm=NULL);

//-Kernels for calculating the Laminar+SPS viscosity.
void ComputeSpsTau(unsigned np,unsigned npb,float smag,float blin
  ,const float4* velrho,const tsymatrix3f* sps2strain,tsymatrix3f* tau_rho2
  ,cudaStream_t stm=NULL);

//-Kernels for Delta-SPH.
void AddDelta(unsigned n,const float* delta,float* ar,cudaStream_t stm=NULL);

//-Kernels for ComputeStep (position).
void ComputeStepPos (byte periactive,bool floatings,unsigned np,unsigned npb
  ,const double2* movxy,const double* movz,double2* posxy,double* posz
  ,unsigned* dcell,typecode* code);
void ComputeStepPos2(byte periactive,bool floatings,unsigned np,unsigned npb
  ,const double2* posxypre,const double* poszpre,const double2* movxy,const double* movz
  ,double2* posxy,double* posz,unsigned* dcell,typecode* code);

//-Kernels for Motion.
void CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini
  ,unsigned idfin,const typecode* code,const unsigned* idp,unsigned* ridp
  ,cudaStream_t stm=NULL);
void LoadPosRef(unsigned pscount,unsigned casenfixed,unsigned np
  ,const double2* posxy,const double* posz,const unsigned* ridpmot
  ,const unsigned* idpref,double3* posref);
void MoveLinBound(byte periactive,unsigned np,unsigned ini,tdouble3 mvpos,tfloat3 mvvel
  ,const unsigned* ridpmot,double2* posxy,double* posz,unsigned* dcell,float4* velrho,typecode* code);
void MoveMatBound(byte periactive,bool simulate2d,unsigned np,unsigned ini,tmatrix4d m,double dt
  ,const unsigned* ridpmot,double2* posxy,double* posz,unsigned* dcell,float4* velrho,typecode* code,float3* boundnor);
void FtNormalsUpdate(unsigned np,unsigned ini,tmatrix4d m,const unsigned* ridpmot,float3* boundnor);

//-Kernels for MLPistons motion.
void MovePiston1d(bool periactive,unsigned np,unsigned idini,double dp,double poszmin
  ,unsigned poszcount,const byte* pistonid,const double* movx,const double* velx
  ,const unsigned* ridpmot,double2* posxy,double* posz,unsigned* dcell,float4* velrho,typecode* code);
void MovePiston2d(bool periactive,unsigned np,unsigned idini,double dp,double posymin
  ,double poszmin,unsigned poszcount,const double* movx,const double* velx,const unsigned* ridpmot
  ,double2* posxy,double* posz,unsigned* dcell,float4* velrho,typecode* code);

//-Kernels for Floating bodies NEW.
void FtPartsSumAce(bool periactive,unsigned ftcount
  ,const float4* ftodata,const double3* ftocenter,const unsigned* ridpmot
  ,const double2* posxy,const double* posz,const float3* ace
  ,float3* ftoacelinang);

void FtPartsUpdate(bool periactive,double dt,bool updatenormals
  ,unsigned np,unsigned fpini,float fradius,tmatrix4d mat
  ,tfloat3 fto_vellin,tfloat3 fto_velang,tdouble3 fto_center
  ,const unsigned* ridpmot,double2* posxy,double* posz,float4* velrho
  ,unsigned* dcell,typecode* code,float3* boundnor,float3* motionvel
  ,float3* motionace,cudaStream_t stm=NULL);


//-Kernels for periodic conditions.
void PeriodicIgnore(unsigned n,typecode* code);
unsigned PeriodicMakeList(unsigned n,unsigned pini,bool stable,unsigned nmax
  ,tdouble3 mapposmin,tdouble3 mapposmax,tdouble3 perinc,const double2* posxy
  ,const double* posz,const typecode* code,unsigned* listp);
void PeriodicDuplicateVerlet(unsigned n,unsigned pini,tuint3 domcells,tdouble3 perinc
  ,const unsigned* listp,unsigned* idp,typecode* code,unsigned* dcell
  ,double2* posxy,double* posz,float4* velrho,tsymatrix3f* spstau,float4* velrhom1);
void PeriodicDuplicateSymplectic(unsigned n,unsigned pini
  ,tuint3 domcells,tdouble3 perinc,const unsigned* listp,unsigned* idp,typecode* code
  ,unsigned* dcell,double2* posxy,double* posz,float4* velrho,tsymatrix3f* spstau
  ,double2* posxypre,double* poszpre,float4* velrhopre);
void PeriodicDuplicateNormals(unsigned n,unsigned pini,const unsigned* listp
  ,float3* normals,float3* motionvel,float3* motionace);

//-Kernels for Damping.
void ComputeDampingPlane(double dt,double4 plane,float dist,float over
  ,float3 factorxyz,float redumax,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,const typecode* code,float4* velrho);
void ComputeDampingPlaneDom(double dt,double4 plane,float dist,float over,float3 factorxyz
  ,float redumax,double zmin,double zmax,double4 pla0,double4 pla1,double4 pla2,double4 pla3
  ,unsigned n,unsigned pini,const double2* posxy,const double* posz,const typecode* code
  ,float4* velrho);
void ComputeDampingBox(unsigned n,unsigned pini,double dt,float3 factorxyz,float redumax
  ,double3 limitmin1,double3 limitmin2,double3 limitmax1,double3 limitmax2
  ,double3 limitover1,double3 limitover2,double3 boxsize1,double3 boxsize2
  ,const double2* posxy,const double* posz,const typecode* code,float4* velrho);
void ComputeDampingCylinder(unsigned n,unsigned pini
  ,double dt,double3 point1,double3 point2,double limitmin
  ,float dist,float over,float3 factorxyz,float redumax
  ,const double2* posxy,const double* posz,const typecode* code
  ,float4* velrho);

//<vs_outpaarts_ini>
//-Kernels for OutputParts.
void ComputeOutputPartsInit(byte resmask,bool selall
  ,unsigned n,unsigned pini,byte* sel);
void ComputeOutputPartsGroup(byte resmask,byte resprev
  ,bool cmband,bool inverse,unsigned n,unsigned pini,byte* sel);
void ComputeOutputPartsPos(byte resmask,bool cmband,bool inverse
  ,double3 pmin,double3 pmax,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,byte* sel);
void ComputeOutputPartsPlane(byte resmask,bool cmband,bool inverse
  ,double4 plane,float maxdist,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,byte* sel);
void ComputeOutputPartsSphere(byte resmask,bool cmband,bool inverse
  ,double3 pcen,float radius2,unsigned n,unsigned pini
  ,const double2* posxy,const double* posz,byte* sel);
void ComputeOutputPartsCylinder(byte resmask,bool cmband,bool inverse
  ,double4 plane,float maxdist,double3 pcen1,double3 pcen2,float radius
  ,unsigned n,unsigned pini,const double2* posxy,const double* posz
  ,byte* sel);
void ComputeOutputPartsType(byte resmask,bool cmband,bool inverse
  ,byte types,unsigned n,unsigned pini,const typecode* code,byte* sel);
void ComputeOutputPartsMk(byte resmask,bool cmband,bool inverse
  ,typecode mkcode1,typecode mkcode2,unsigned n,unsigned pini
  ,const typecode* code,byte* sel);

//<vs_outpaarts_end>

//<vs_flexstruc_ini>
void SetFlexStrucClampCodes(unsigned npb,const float4* poscell,const StFlexStrucData* flexstrucdata,typecode* code);
unsigned CountFlexStrucParts(unsigned npb,const typecode* code);
void CalcFlexStrucRidp(unsigned npb,const typecode* code,unsigned* flexstrucridp);
void GatherToFlexStrucArray(unsigned npfs,const unsigned* flexstrucridp,const float4* fullarray,float4* flexstrucarray);
void GatherToFlexStrucArray(unsigned npfs,const unsigned* flexstrucridp,const float3* fullarray,float3* flexstrucarray);
unsigned CountFlexStrucPairs(unsigned npfs,const float4* poscell0,unsigned* numpairs);
void SetFlexStrucPairs(unsigned npfs,const float4* poscell0,unsigned** pairidx);
void CalcFlexStrucKerCorr(const StInterParmsFlexStrucg& tfs);
void UpdateFlexStrucGeometry(const StInterParmsFlexStrucg& tfs);
void Interaction_ForcesFlexStruc(const StInterParmsFlexStrucg& tfs);
void ComputeStepPosFlexStruc(unsigned npfs,const unsigned* flexstrucridp
    ,const double2* posxypre,const double* poszpre,const double2* movxy,const double* movz
    ,double2* posxy,double* posz,unsigned* dcell,typecode* code);
bool FlexStrucStepIsValid(unsigned npb,const typecode* code);
//<vs_flexstruc_end>

}


#endif


