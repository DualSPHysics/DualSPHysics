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

/// \file JCellDivDataCpu.h \brief Declares structures and implements inline functions for neighborhood search.

#ifndef _JCellDivDataCpu_
#define _JCellDivDataCpu_

#include "DualSphDef.h"

///Structure with data of cell division for neighborhood search on GPU.
typedef struct{
  int scelldiv; ///<Value to divide KernelSize (1 or 2).
  tint4 nc;
  unsigned cellfluid;
  tint3 cellzero;
  const unsigned* begincell;
  float scell;
  unsigned domcellcode;
  tdouble3 domposmin;
}StDivDataCpu;

//==============================================================================
///Returns empty StDivDataCpu structure.
//==============================================================================
inline StDivDataCpu DivDataCpuNull(){
  StDivDataCpu c={0,TInt4(0),0,TInt3(0),NULL,0,0,TDouble3(0)};
  return(c);
}

//==============================================================================
/// Returns structure with data for neighborhood search on Single-GPU.
//==============================================================================
inline StDivDataCpu MakeDivDataCpu(int scelldiv,const tuint3 &ncells,const tuint3 &cellmin
  ,const unsigned* begincell,float scell,unsigned domcellcode,const tdouble3 &domposmin)
{
  StDivDataCpu ret;
  ret.scelldiv=scelldiv;
  ret.nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
  ret.cellfluid=ret.nc.w*ret.nc.z+1;
  ret.cellzero=ToTInt3(cellmin);
  ret.begincell=begincell;
  ret.scell=scell;
  ret.domcellcode=domcellcode;
  ret.domposmin=domposmin;
  return(ret);
}


///Structure with data for neighborhood search.
typedef struct{
  int cellinit;
  int cxini,cxfin;
  int yini,yfin;
  int zini,zfin;
}StNgSearch;


#endif


