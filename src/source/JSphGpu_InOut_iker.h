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

/// \file JSphGpu_InOut_iker.h \brief Declares functions and CUDA kernels for InOut feature.

#ifndef _JSphGpu_InOut_iker_
#define _JSphGpu_InOut_iker_

#include "DualSphDef.h"
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>

#define INOUT_RefillAdvanced_MASK 0x01
#define INOUT_RefillSpFull_MASK 0x02
#define INOUT_RemoveInput_MASK 0x04
#define INOUT_RemoveZsurf_MASK 0x08
#define INOUT_ConvertInput_MASK 0x10

/// Implements a set of functions and CUDA kernels for InOut feature.
namespace cusphinout{

//inline tfloat3 ToTFloat3(const float3& v){ return(TFloat3(v.x,v.y,v.z)); }

//-Kernels for inlet/outlet (JSphInOut).
void InOutIgnoreFluidDef(unsigned n,typecode cod,typecode codnew,typecode *code);
void UpdatePosFluid(byte periactive,unsigned n,unsigned pini
  ,double2 *posxy,double *posz,unsigned *dcell,typecode *code);
unsigned InOutCreateListSimple(bool stable,unsigned n,unsigned pini
  ,const typecode *code,unsigned *listp);
unsigned InOutCreateList(bool stable,unsigned n,unsigned pini
  ,byte chkinputmask,byte nzone,const byte *cfgzone,const float4 *planes
  ,tfloat3 freemin,tfloat3 freemax
  ,const float2 *boxlimit,const double2 *posxy,const double *posz
  ,typecode *code,unsigned *listp);

void InOutSetAnalyticalData(unsigned n,const unsigned *listp
  ,byte izone,byte rmode,byte vmode,byte vprof,byte refillspfull
  ,float timestep,float zsurfv,tfloat4 veldata,tfloat4 veldata2,tfloat3 dirdata
  ,float coefhydro,float rhopzero,float gamma
  ,const typecode *code,const double *posz,const float *zsurfpart,float4 *velrhop);

void InoutClearInteractionVars(unsigned npf,unsigned pini,const typecode *code
    ,float3 *ace,float *ar,float *viscdt,float4 *shiftposfs);

void InOutUpdateVelrhopM1(unsigned n,const int *inoutpart
    ,const float4 *velrhop,float4 *velrhopm1);

void InOutComputeStep(unsigned n,int *inoutpart,const float4 *planes
  ,const float *width,const byte *cfgupdate,const float *zsurfv,typecode codenewpart
  ,const double2 *posxy,const double *posz,const byte *zsurfok
  ,typecode *code,byte *newizone);
unsigned InOutListCreate(bool stable,unsigned n,unsigned nmax,const byte *newizone,int *inoutpart);
void InOutCreateNewInlet(byte periactive,unsigned newn
  ,const unsigned *inoutpart,unsigned inoutcount,const byte *newizone
  ,unsigned np,unsigned idnext,typecode codenewpart,const float3 *dirdata,const float *width
  ,double2 *posxy,double *posz,unsigned *dcell,typecode *code,unsigned *idp,float4 *velrhop);

//-Kernels for inlet/outlet filling (JSphInOut).
void InOutFillMove(byte periactive,unsigned n,const unsigned *inoutpart
  ,double dt,const float4 *velrhop
  ,double2 *posxy,double *posz,unsigned *dcell,typecode *code);
void InOutFillProjection(unsigned n,const unsigned *inoutpart
  ,const byte *cfgupdate,const float4 *planes,const double2 *posxy,const double *posz
  ,const typecode *code,float *prodist,double2 *proposxy,double *proposz);
unsigned InOutFillListCreate(bool stable,unsigned npt
  ,const double2 *ptposxy,const double *ptposz,const byte *zsurfok
  ,const byte *ptzone,const byte *cfgupdate,const float *zsurf,const float *width
  ,unsigned npropt,const float *prodist,const double2 *proposxy,const double *proposz
  ,float dpmin,float dpmin2,float dp,float *ptdist,unsigned nmax,unsigned *inoutpart);
void InOutFillCreate(byte periactive,unsigned newn,const unsigned *newinoutpart
  ,const double2 *ptposxy,const double *ptposz,const byte *ptzone,const float *ptauxdist
  ,unsigned np,unsigned idnext,typecode codenewpart,const float3 *dirdata
  ,double2 *posxy,double *posz,unsigned *dcell,typecode *code,unsigned *idp,float4 *velrhop);

//-Kernels to extrapolate rhop and velocity (JSphInOut).
void Interaction_InOutExtrap(byte doublemode,bool simulate2d,TpKernel tkernel
  ,unsigned inoutcount,const int *inoutpart,const byte *cfgzone,byte computerhopmask,byte computevelmask
  ,const float4 *planes,const float* width,const float3 *dirdata,float determlimit
  ,const StDivDataGpu &dvd,const double2 *posxy,const double *posz,const typecode *code
  ,const unsigned *idp,float4 *velrhop);

//-Kernels to extrapolate rhop on boundary particles (JSphBoundCorr).
void Interaction_BoundCorr(byte doublemode,bool simulate2d,TpKernel tkernel
  ,unsigned npbok,typecode boundcode,tfloat4 plane,tfloat3 direction,float determlimit
  ,const StDivDataGpu &dvd,const double2 *posxy,const double *posz
  ,const typecode *code,const unsigned *idp,float4 *velrhop);

//-Kernels to interpolate velocity (JSphInOutGridDataTime).
void InOutInterpolateTime(unsigned npt,double time,double t0,double t1
  ,const float *velx0,const float *velx1,float *velx
  ,const float *velz0,const float *velz1,float *velz);
void InOutInterpolateZVel(unsigned izone,double posminz,double dpz,int nz1
  ,const float *velx,const float *velz,unsigned np,const int *plist
  ,const double *posz,const typecode *code,float4 *velrhop,float velcorr);
void InOutInterpolateResetZVel(unsigned izone,unsigned np,const int *plist
  ,const typecode *code,float4 *velrhop);

}


#endif


