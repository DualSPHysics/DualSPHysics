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

/// \file JCellDivDataGpu.h \brief Declares structures and implements inline functions for neighborhood search.

#ifndef _JCellDivDataGpu_
#define _JCellDivDataGpu_

#include "DualSphDef.h"
#include <cuda_runtime_api.h>

///Structure with data of cell division for neighborhood search on GPU.
typedef struct{
  TpMgDivMode axis;
  tuint3 ncells;
  int scelldiv;         ///<Value to divide KernelSize (1 or 2).
  int4 nc;
  unsigned cellfluid;
  int3 cellzero;
  const int2* beginendcell; //- int2*
  float scell;
  unsigned domcellcode;
  double3 domposmin;
  float kernelsize2;   ///<Maximum interaction distance squared (KernelSize^2).
  float poscellsize;   ///<Size of cells used for coding PosCell (it is usually KernelSize).
}StDivDataGpu;

//==============================================================================
///Returns empty StDivDataGpu structure.
//==============================================================================
inline StDivDataGpu DivDataGpuNull(){
  StDivDataGpu c={MGDIV_None,TUint3(0),0,{0,0,0,0},0,{0,0,0},NULL,0,0,{0,0,0},0,0};
  return(c);
}

//==============================================================================
/// Returns structure with data for neighborhood search on Single-GPU.
//==============================================================================
inline StDivDataGpu MakeDivDataGpu(int scelldiv,const tuint3 &ncells
  ,const tuint3 &cellmin,const int2* beginendcell,float scell,unsigned domcellcode
  ,const tdouble3 &domposmin,float kernelsize2,float poscellsize)
{
  StDivDataGpu ret;
  ret.axis=MGDIV_Z;
  ret.scelldiv=scelldiv;
  ret.ncells=ncells;
  ret.nc.x=int(ncells.x);
  ret.nc.y=int(ncells.y);
  ret.nc.z=int(ncells.z);
  ret.nc.w=int(ncells.x*ncells.y);//-For Single-GPU.
  ret.cellfluid=ret.nc.w*ret.nc.z+1;
  ret.cellzero.x=int(cellmin.x);
  ret.cellzero.y=int(cellmin.y);
  ret.cellzero.z=int(cellmin.z);
  ret.beginendcell=beginendcell;
  ret.scell=scell;
  ret.domcellcode=domcellcode;
  ret.domposmin.x=domposmin.x;
  ret.domposmin.y=domposmin.y;
  ret.domposmin.z=domposmin.z;
  ret.kernelsize2=kernelsize2;
  ret.poscellsize=poscellsize;
  return(ret);
}


#endif


