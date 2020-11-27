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

/// \file FunSphKernelsCfg.h \brief Declares basic/general functions for SPH kenernels.

#ifndef _FunSphKernels_
#define _FunSphKernels_

#include "FunSphKernelDef.h"
#include "DualSphDef.h"
#include <vector>


/// Implements a set of basic/general functions related to SPH.
namespace fsph{
void RunExceptioonFun(const std::string &srcfile,int srcline,const std::string &fun
  ,const std::string &msg,const std::string &file="");

float GetKernelFactor(TpKernel tkernel);
std::string GetKernelName(TpKernel tkernel);
void GetKernelConfig(const StCteSph &CSP,std::vector<std::string> &lines);

}

#endif


