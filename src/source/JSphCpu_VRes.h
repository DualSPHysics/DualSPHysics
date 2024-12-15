//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2023 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphCpu.h \brief Declares the class \ref JSphCpu.

#ifndef _JSphCpu_VRes_
#define _JSphCpu_VRes_


#include "DualSphDef.h"
#include "JDsTimersCpu.h"
#include "JCellDivDataCpu.h"
#include "JSph.h"
#include "JArraysCpu.h"
#include "JSphCpu.h"


typedef struct{
  unsigned np,npb,npbok,npf; // npf=np-npb
  StDivDataCpu divdata;
  const unsigned *dcell;
  const tdouble3 *pos;
  const tfloat4 *velrhop;
  const unsigned *idp;
  const typecode *code;
   StCteSph csp;
}stinterparmscb;


namespace fvres{
//-Code for VRes in JSphCpu_Buffer.cpp
  //--------------------------------------

  template<bool sim2d,TpKernel tker,unsigned order> void InteractionBufferExtrap(unsigned bufferpartcount,const int *bufferpart,
			StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,
			const typecode *code,const unsigned *idp,const tfloat4 *velrhop,const StCteSph csp,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb);

  template<TpKernel tker> void Interaction_BufferExtrapT(unsigned bufferpartcount,
	const int *bufferpart,const stinterparmscb &t,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb,unsigned order);

  void Interaction_BufferExtrap(unsigned bufferpartcount,const int *bufferpart
     ,const stinterparmscb &t,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb,unsigned order);

  template<bool sim2d,TpKernel tker,unsigned order> void InteractionBufferExtrapFlux(const unsigned n,const int pini,
		StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,
		const typecode *code,const unsigned *idp,const tfloat4 *velrhop,const StCteSph csp,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold);

  template<TpKernel tker> void Interaction_BufferExtrapFluxT(const unsigned n,const int pini
  ,const stinterparmscb &t,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold);

  void Interaction_BufferExtrapFlux(const unsigned n,const int pini
  ,const stinterparmscb &t,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold);

  inline bool InZone(const tdouble3 &ps,const tdouble3 &boxlimitmin,tdouble3 &boxlimitmax){
    return (boxlimitmin.x <= ps.x && ps.x <= boxlimitmax.x && boxlimitmin.y <= ps.y 
    && ps.y <= boxlimitmax.y && boxlimitmin.z <= ps.z && ps.z <= boxlimitmax.z);
  }
  tdouble3 MovePoint(tdouble3 oldpos,const tmatrix4d& mat);

  void BufferShiftingCpu(unsigned n,unsigned pinit,const tdouble3 *pos
    ,tfloat4 *shiftpos,const typecode *code,StrGeomVresCpu& vresdata);

  void CheckMassFlux(unsigned n,unsigned pinit,const StCteSph csp
    ,const StDivDataCpu& divdata,const unsigned* dcell,const tdouble3* pos
    ,const typecode* code,tdouble3* posb,tfloat3 *normals,float *fluxes);

  unsigned CountFreeSurfaceParticles(unsigned npf,unsigned pini
    ,const unsigned* fstype,unsigned* listp);
  
  template<TpKernel tker,bool sim2d> void InteractionComputeFSNormals
    (unsigned np,unsigned pinit,const StCteSph csp
    ,StDivDataCpu divdata,const unsigned* dcell
    ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
    ,const unsigned* listp,unsigned* fstype,tfloat3* fsnormal,StrGeomVresCpu& vresdata);
  void CallComputeFSNormals(const unsigned np,const unsigned npb
    ,const StCteSph csp,const StDivDataCpu& divdata,const unsigned* dcell
    ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
    ,unsigned* fstype,tfloat3* fsnormal,unsigned* listp,StrGeomVresCpu& vresdata);

  void InteractionCallScanUmbrellaRegion(unsigned n,unsigned pinit,const StCteSph csp
    ,StDivDataCpu divdata,const unsigned* dcell,const tdouble3* pos
    ,const typecode* code,const tfloat3* fsnormal,const unsigned* listp
    ,unsigned* fstype,StrGeomVresCpu& vresdata);
  void CallScanUmbrellaRegion(const unsigned np,const unsigned npb
    ,const StCteSph csp,const StDivDataCpu& divdata
    ,const unsigned* dcell,const tdouble3* pos,const typecode* code
    ,const tfloat3* fsnormal,unsigned* listp,unsigned* fstype,StrGeomVresCpu& vresdata);

  template<TpKernel tker,bool sim2d> void CorrectShiftBuff
    (const unsigned n,const unsigned pinit  ,const StCteSph csp
    ,const unsigned* dcell,const tdouble3* pos,const typecode* code
    ,tfloat4* shiftvel,unsigned* fstype,StrGeomVresCpu& vresdata);

  void CallCorrectShiftBuff(const unsigned np,const unsigned npb
    ,const StCteSph csp,const unsigned* dcell,const tdouble3* pos,const typecode* code
    ,tfloat4* shiftvel,unsigned* fstype,StrGeomVresCpu& vresdata);
}

  #endif