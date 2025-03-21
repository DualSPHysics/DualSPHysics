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

/// \file JSphGpu_cte.h \brief Defines constant memory structure for GPU particle interactions.

#ifndef _JSphGpu_cte_
#define _JSphGpu_cte_

#include "FunSphKernelDef.h"

/// Structure with constants stored in the constant memory of GPU for the particle interactions.
typedef struct{
  float massb;              ///<Reference mass of the general boundary particle [kg].
  float massf;              ///<Reference mass of the fluid particle [kg].
  float kernelh;            ///<The smoothing length of SPH kernel [m].
  float kernelsize2;        ///<Maximum interaction distance squared (KernelSize^2).
  float poscellsize;        ///<Size of cells used for coding PosCell (it is usually KernelSize).
  float awen;               ///<Wendland kernel constant (awen) to compute wab.
  float bwenh;              ///<Wendland kernel constant (bwenh) to compute fac (kernel derivative).
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
  //-Values on formulation options.
  unsigned tboundary;  
  //-Values for open periodic boundaries.
  unsigned periactive;
  double xperincx,xperincy,xperincz;
  double yperincx,yperincy,yperincz;
  double zperincx,zperincy,zperincz;
  //-Values for map definition.
  double maprealposminx,maprealposminy,maprealposminz;
  double maprealsizex,maprealsizey,maprealsizez;
  //-Values depending on the assigned domain (can change).
  unsigned axis;
  unsigned cellcode;
  double domposminx,domposminy,domposminz;
  //-Ctes. for Cubic Spline kernel.
  float cubic_a1,cubic_a2,cubic_aa,cubic_a24,cubic_c1,cubic_d1,cubic_c2,cubic_odwdeltap;
}StCteInteraction; 

//==============================================================================
/// Set mass particle values.
//==============================================================================
inline void SetCtegMass(StCteInteraction& cte,float massb,float massf){
  cte.massb=massb;  cte.massf=massf;
}

//==============================================================================
/// Set distance values depending on h and dp.
//==============================================================================
inline void SetCtegKsize(StCteInteraction& cte,float kernelh
  ,float kernelsize,float kernelsize2,float poscellsize
  ,float eta2,float dp,float scell,float movlimit)
{
  cte.kernelh=kernelh;
  cte.kernelsize=kernelsize;
  cte.kernelsize2=kernelsize2;
  cte.poscellsize=poscellsize;
  cte.eta2=eta2;
  cte.dp=dp;
  cte.scell=scell;
  cte.movlimit=movlimit;
}

//==============================================================================
/// Set values for density and pressure.
//==============================================================================
inline void SetCtegRho(StCteInteraction& cte,float rhopzero,float ovrhopzero
  ,float gamma,float cs0,float cteb)
{
  cte.rhopzero=rhopzero;
  cte.ovrhopzero=ovrhopzero;
  cte.gamma=gamma;
  cte.cs0=cs0;
  cte.cteb=cteb;
}

//==============================================================================
/// Set values for DDT.
//==============================================================================
inline void SetCtegDdt(StCteInteraction& cte,float ddtkh,float ddtgz){
  cte.ddtkh=ddtkh;  cte.ddtgz=ddtgz;
}

//==============================================================================
/// Set values for map definition.
//==============================================================================
inline void SetCtegMap(StCteInteraction& cte,const tdouble3& maprealposmin
  ,const tdouble3& maprealsize)
{
  cte.maprealposminx=maprealposmin.x;
  cte.maprealposminy=maprealposmin.y;
  cte.maprealposminz=maprealposmin.z;
  cte.maprealsizex=maprealsize.x;
  cte.maprealsizey=maprealsize.y;
  cte.maprealsizez=maprealsize.z;
}

//==============================================================================
/// Set values depending on the assigned domain (can change).
//==============================================================================
inline void SetCtegDomain(StCteInteraction& cte,TpMgDivMode axis
  ,unsigned cellcode,const tdouble3& dompmin)
{
  cte.axis=unsigned(axis);
  cte.cellcode=cellcode;
  cte.domposminx=dompmin.x; cte.domposminy=dompmin.y; cte.domposminz=dompmin.z;
}

//==============================================================================
/// Set values on formulation options.
//==============================================================================
inline void SetCtegOpts(StCteInteraction& cte,unsigned tboundary){
  cte.tboundary=tboundary;  
}

//==============================================================================
/// Set values for open periodic boundaries.
//==============================================================================
inline void SetCtegPeriodic(StCteInteraction& cte,byte pactive
  ,const tdouble3& pxinc,const tdouble3& pyinc,const tdouble3& pzinc)
{
  cte.periactive=pactive;
  cte.xperincx=pxinc.x; cte.xperincy=pxinc.y; cte.xperincz=pxinc.z;
  cte.yperincx=pyinc.x; cte.yperincy=pyinc.y; cte.yperincz=pyinc.z;
  cte.zperincx=pzinc.x; cte.zperincy=pzinc.y; cte.zperincz=pzinc.z;
}

//==============================================================================
/// Set ctes. for Cubic Spline kernel.
//==============================================================================
inline void SetCtegKerCubic(StCteInteraction& cte
  ,const fsph::StKCubicCte& cker)
{
  cte.cubic_a1=cker.a1; cte.cubic_a2=cker.a2;
  cte.cubic_aa=cker.aa; cte.cubic_a24=cker.a24;
  cte.cubic_c1=cker.c1; cte.cubic_c2=cker.c2;
  cte.cubic_d1=cker.d1; cte.cubic_odwdeltap=cker.od_wdeltap;
}

//==============================================================================
/// Set ctes. for Wendland kernel.
//==============================================================================
inline void SetCtegKerWendland(StCteInteraction& cte
  ,const fsph::StKWendlandCte& cker)
{
  cte.awen =cker.awen;
  cte.bwenh=cker.bwenh;
}


#endif

