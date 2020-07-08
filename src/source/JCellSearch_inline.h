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

/// \file JCellSearch_inline.h \brief implements inline functions for neighborhood search on CPU.

#ifndef _JCellSearch_inline_
#define _JCellSearch_inline_

#include "JCellDivDataCpu.h"

/// Implements inline functions for neighborhood search on CPU.
namespace nsearch{

//==============================================================================
/// Return initial data for neighborhood search according to cell number of particle.
/// Devuelve datos iniciales para busqueda de vecinos segun numero de celda de particula.
//==============================================================================
inline StNgSearch Init(unsigned rcell,bool boundp2,const StDivDataCpu &dvd){
  //-Get cell coordinates of cell number.
  const int cx=PC__Cellx(dvd.domcellcode,rcell)-dvd.cellzero.x;
  const int cy=PC__Celly(dvd.domcellcode,rcell)-dvd.cellzero.y;
  const int cz=PC__Cellz(dvd.domcellcode,rcell)-dvd.cellzero.z;
  //-Code for scelldiv 1 or 2 but not zero. | Codigo para scelldiv 1 o 2 pero no cero.
  StNgSearch ret;
  ret.cellinit=(boundp2? 0: dvd.cellfluid);
  ret.cxini=cx-(cx<dvd.scelldiv? cx: dvd.scelldiv);
  ret.cxfin=cx+(dvd.nc.x-cx-1<dvd.scelldiv? dvd.nc.x-cx-1: dvd.scelldiv)+1;
  ret.yini=cy-(cy<dvd.scelldiv? cy: dvd.scelldiv);
  ret.yfin=cy+(dvd.nc.y-cy-1<dvd.scelldiv? dvd.nc.y-cy-1: dvd.scelldiv)+1;
  ret.zini=cz-(cz<dvd.scelldiv? cz: dvd.scelldiv);
  ret.zfin=cz+(dvd.nc.z-cz-1<dvd.scelldiv? dvd.nc.z-cz-1: dvd.scelldiv)+1;
  return(ret);
}

//==============================================================================
/// Return initial data for neighborhood search according to position.
/// Devuelve datos iniciales para busqueda de vecinos segun posicion.
//==============================================================================
inline StNgSearch Init(const tdouble3 &pos,bool boundp2,const StDivDataCpu &dvd){
  //-Get cell coordinates of position pos.
  const int cx=int((pos.x-dvd.domposmin.x)/dvd.scell)-dvd.cellzero.x;
  const int cy=int((pos.y-dvd.domposmin.y)/dvd.scell)-dvd.cellzero.y;
  const int cz=int((pos.z-dvd.domposmin.z)/dvd.scell)-dvd.cellzero.z;
  //-Code for scelldiv 1 or 2 but not zero. | Codigo para scelldiv 1 o 2 pero no cero.
  StNgSearch ret;
  ret.cellinit=(boundp2? 0: dvd.cellfluid);
  ret.cxini=cx-(cx<dvd.scelldiv? cx: dvd.scelldiv);
  ret.cxfin=cx+(dvd.nc.x-cx-1<dvd.scelldiv? dvd.nc.x-cx-1: dvd.scelldiv)+1;
  ret.yini=cy-(cy<dvd.scelldiv? cy: dvd.scelldiv);
  ret.yfin=cy+(dvd.nc.y-cy-1<dvd.scelldiv? dvd.nc.y-cy-1: dvd.scelldiv)+1;
  ret.zini=cz-(cz<dvd.scelldiv? cz: dvd.scelldiv);
  ret.zfin=cz+(dvd.nc.z-cz-1<dvd.scelldiv? dvd.nc.z-cz-1: dvd.scelldiv)+1;
  return(ret);
}

//==============================================================================
/// Returns range of particles for neighborhood search.
/// Devuelve rango de particulas para busqueda de vecinos.
//==============================================================================
inline tuint2 ParticleRange(int y,int z,const StNgSearch &ngs,const StDivDataCpu &dvd){
  const int v=dvd.nc.w*z + dvd.nc.x*y + ngs.cellinit;
  const unsigned pini=dvd.begincell[v+ngs.cxini];
  const unsigned pfin=dvd.begincell[v+ngs.cxfin];
  return(TUint2(pini,pfin));
}

//==============================================================================
/// Returns distance between particles 1 and 2 (drx,dry,drz and rr2).
//==============================================================================
inline tfloat4 Distances(const tdouble3 &pos1,const tdouble3 &pos2){
  tfloat4 dr;
  dr.x=float(pos1.x-pos2.x);
  dr.y=float(pos1.y-pos2.y);
  dr.z=float(pos1.z-pos2.z);
  dr.w=dr.x*dr.x + dr.y*dr.y + dr.z*dr.z;
  return(dr);
}

//==============================================================================
/// Returns distance squared between particles 1 and 2 (rr2).
//==============================================================================
inline float Distance2(const tdouble3 &pos1,const tdouble3 &pos2){
  const float drx=float(pos1.x-pos2.x);
  const float dry=float(pos1.y-pos2.y);
  const float drz=float(pos1.z-pos2.z);
  return(drx*drx + dry*dry + drz*drz);
}


}
#endif


