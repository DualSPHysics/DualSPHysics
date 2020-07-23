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

/// \file JCellDivGpuSingle_ker.h \brief Declares functions and CUDA kernels to compute operations of the Neighbour List.

#ifndef _JCellDivGpuSingle_ker_
#define _JCellDivGpuSingle_ker_

#include "JCellDivGpu_ker.h"

namespace cudiv{

void PreSortFull(unsigned np,unsigned cellcode,const unsigned *dcell,const typecode *code,tuint3 cellmin,tuint3 ncells,unsigned *cellpart,unsigned *sortpart);
void PreSortFluid(unsigned npf,unsigned pini,unsigned cellcode,const unsigned *dcell,const typecode *code,tuint3 cellmin,tuint3 ncells,unsigned *cellpart,unsigned *sortpart);

}
#endif


