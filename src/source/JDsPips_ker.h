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

/// \file JDsPips_ker.h \brief Declares functions and CUDA kernels for PIPS calculation on GPU.

#ifndef _JDsPips_ker_
#define _JDsPips_ker_

#include "DualSphDef.h"
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>


/// Implements a set of functions and CUDA kernels for PIPS calculation.
namespace cupips{

unsigned InteractionNgSize_1st(unsigned n);
void InteractionNg_1st(unsigned nb,unsigned pinitb,unsigned nf,unsigned pinitf
  ,const StDivDataGpu &dvd,const unsigned *dcell,const float4 *poscell,uint4 *res
  ,cudaStream_t stm=NULL);

unsigned InteractionNgSize_2nd(unsigned n);
void InteractionNg_2nd(unsigned n,const uint4 *data,uint4 *res,cudaStream_t stm=NULL);

unsigned InteractionNgSize_3th(unsigned n);
void InteractionNg_3th(unsigned n,const uint4 *data,ullong *res,cudaStream_t stm=NULL);

}


#endif


