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

/// \file JCellDivDatacpu.h \brief Declares structures and implements inline functions for neighborhood search.

#ifndef _JCellDivDataCpu_
#define _JCellDivDataCpu_

#include "DualSphDef.h"

///Structure with data of cell division for neighborhood search on GPU.
typedef struct{
  int hdiv; //hdiv=(cellmode==CELLMODE_H? 2: 1)
  tint4 nc;
  unsigned cellfluid;
  tint3 cellzero;
  const unsigned* begincell;
  float scell;
  unsigned domcellcode;
  tdouble3 domposmin;
}StDivDataCpu;

///Structure with data for neighborhood search.
typedef struct{
  int cellinit;
  int cxini,cxfin;
  int yini,yfin;
  int zini,zfin;
}StNgSearch;


#endif


