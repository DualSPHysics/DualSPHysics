//HEAD_DSCODES
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

/// \file JCellSearch_iker.h \brief implements CUDA device functions for neighborhood search.

#ifndef _JCellSearch_iker_
#define _JCellSearch_iker_

#include "DualSphDef.h"

/// Implements CUDA device functions for neighborhood search.
namespace cunsearch{
//#define CTE_AVAILABLE
#ifdef CTE_AVAILABLE
//------------------------------------------------------------------------------
/// Returns cell limits for the interaction.
/// Devuelve limites de celdas para interaccion.
//------------------------------------------------------------------------------
__device__ void InitCte(unsigned rcell
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  //-Obtains interaction limits.
  const int cx=PC__Cellx(CTE.cellcode,rcell)-cellzero.x;
  const int cy=PC__Celly(CTE.cellcode,rcell)-cellzero.y;
  const int cz=PC__Cellz(CTE.cellcode,rcell)-cellzero.z;
  //-Sorts components according axis used in cell order.
  const int c1=(CTE.axis==MGDIV_X? cy: cx);
  const int c2=(CTE.axis==MGDIV_Z? cy: cz);
  const int c3=(CTE.axis==MGDIV_Z? cz: (CTE.axis==MGDIV_X? cx: cy));
  //-Code for scelldiv 1 or 2 but not zero.
  //-Codigo para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  //-Compute absolute limits.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}

//------------------------------------------------------------------------------
/// Returns cell limits for the interaction.
/// Devuelve limites de celdas para interaccion.
//------------------------------------------------------------------------------
__device__ void InitCte(const double &px,const double &py,const double &pz
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  //-Obtains interaction limits.
  const int cx=int((px-CTE.domposminx)/CTE.scell)-cellzero.x;
  const int cy=int((py-CTE.domposminy)/CTE.scell)-cellzero.y;
  const int cz=int((pz-CTE.domposminz)/CTE.scell)-cellzero.z;
  //-Sorts components according axis used in cell order.
  const int c1=(CTE.axis==MGDIV_X? cy: cx);
  const int c2=(CTE.axis==MGDIV_Z? cy: cz);
  const int c3=(CTE.axis==MGDIV_Z? cz: (CTE.axis==MGDIV_X? cx: cy));
  //-Code for scelldiv 1 or 2 but not zero.
  //-Codigo para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  //-Compute absolute limits.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}

#endif

//------------------------------------------------------------------------------
/// Returns cell limits for the interaction.
/// Devuelve limites de celdas para interaccion.
//------------------------------------------------------------------------------
__device__ void Initsp(unsigned rcell
  ,const unsigned &axis,const unsigned &cellcode
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  //-Obtains interaction limits.
  const int cx=PC__Cellx(cellcode,rcell)-cellzero.x;
  const int cy=PC__Celly(cellcode,rcell)-cellzero.y;
  const int cz=PC__Cellz(cellcode,rcell)-cellzero.z;
  //-Sorts components according axis used in cell order.
  const int c1=(axis==MGDIV_X? cy: cx);
  const int c2=(axis==MGDIV_Z? cy: cz);
  const int c3=(axis==MGDIV_Z? cz: (axis==MGDIV_X? cx: cy));
  //-Code for scelldiv 1 or 2 but not zero.
  //-Codigo para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  //-Compute absolute limits.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}

//------------------------------------------------------------------------------
/// Returns cell limits for the interaction.
/// Devuelve limites de celdas para interaccion.
//------------------------------------------------------------------------------
__device__ void Initsp(const double &px,const double &py,const double &pz
  ,const unsigned &axis,const double3 &domposmin,const float &scell
  ,const int &scelldiv,const int4 &nc,const int3 &cellzero
  ,int &ini1,int &fin1,int &ini2,int &fin2,int &ini3,int &fin3)
{
  //-Obtains interaction limits.
  const int cx=int((px-domposmin.x)/scell)-cellzero.x;
  const int cy=int((py-domposmin.y)/scell)-cellzero.y;
  const int cz=int((pz-domposmin.z)/scell)-cellzero.z;
  //-Sorts components according axis used in cell order.
  const int c1=(axis==MGDIV_X? cy: cx);
  const int c2=(axis==MGDIV_Z? cy: cz);
  const int c3=(axis==MGDIV_Z? cz: (axis==MGDIV_X? cx: cy));
  //-Code for scelldiv 1 or 2 but not zero.
  //-Codigo para scelldiv 1 o 2 pero no cero.
  ini1=c1-min(c1,scelldiv);
  fin1=c1+min(nc.x-c1-1,scelldiv)+1;
  ini2=c2-min(c2,scelldiv);
  fin2=c2+min(nc.y-c2-1,scelldiv)+1;
  ini3=c3-min(c3,scelldiv);
  fin3=c3+min(nc.z-c3-1,scelldiv)+1;
  //-Compute absolute limits.
  ini3*=nc.w; fin3*=nc.w;
  ini2*=nc.x; fin2*=nc.x;
}


//------------------------------------------------------------------------------
/// Returns range of boundary particles for neighborhood search.
/// Devuelve rango de particulas bound para busqueda de vecinos.
//------------------------------------------------------------------------------
__device__ void ParticleRange(const int &c2,const int &c3
  ,const int &ini1,const int &fin1,const int2 *begincell
  ,unsigned &pini,unsigned &pfin)
{
  const int v=c2+c3;
  for(int c1=ini1;c1<fin1;c1++){
    const int2 cbeg=begincell[c1+v];
    if(cbeg.y){
      if(!pfin)pini=cbeg.x;
      pfin=cbeg.y;
    }
  }
}

//==============================================================================
/// Returns distance squared between particles 1 and 2 (rr2).
//==============================================================================
__device__ float Distance2(const float4 &pscellp1,const float4 &pscellp2
  ,const float &poscellsize)
{
  const float drx=pscellp1.x-pscellp2.x + poscellsize*(CEL_GetfX(pscellp1.w)-CEL_GetfX(pscellp2.w));
  const float dry=pscellp1.y-pscellp2.y + poscellsize*(CEL_GetfY(pscellp1.w)-CEL_GetfY(pscellp2.w));
  const float drz=pscellp1.z-pscellp2.z + poscellsize*(CEL_GetfZ(pscellp1.w)-CEL_GetfZ(pscellp2.w));
  return(drx*drx + dry*dry + drz*drz);
}








}

#endif


