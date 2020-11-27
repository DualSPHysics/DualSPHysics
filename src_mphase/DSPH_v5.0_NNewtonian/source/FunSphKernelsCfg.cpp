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

/// \file FunSphKernelsCfg.cpp \brief Implements basic/general functions for SPH kenernels.

#include "FunSphKernelsCfg.h"
#include "FunSphKernel.h"
#include "Functions.h"
//#include <limits>
//#include <cstdio>
//#include <cstdlib>
//#include <cstring>
//#include <cfloat>
//#include <cmath>
//#include <climits>
//#include <algorithm>

using namespace std;

namespace fsph{

//==============================================================================
/// Throws an exception related to a file or not.
//==============================================================================
void RunExceptioonFun(const std::string &srcfile,int srcline,const std::string &fun
  ,const std::string &msg,const std::string &file)
{ 
  fun::RunExceptioonFun(srcfile,srcline,fun,msg,file);
}

//==============================================================================
/// Returns factor value to compute KernelSize according to KernelH.
//==============================================================================
float GetKernelFactor(TpKernel tkernel){
  float kh=GetKernel_Factor(tkernel);
  if(!kh)Run_ExceptioonFun("Kernel unknown.");
  return(kh);
}

//==============================================================================
/// Returns the name of the kernel in text format.
//==============================================================================
std::string GetKernelName(TpKernel tkernel){
  string tx;
  switch(tkernel){
    case KERNEL_Cubic:      tx=GetKernelCubic_Name();       break;
    case KERNEL_Wendland:   tx=GetKernelWendland_Name();    break;
    default: Run_ExceptioonFun("Kernel unknown.");
  }
  return(tx);
}

//==============================================================================
/// Returns strings with Kernel configuration.
//==============================================================================
void GetKernelConfig(const StCteSph &CSP,std::vector<std::string> &lines){
  lines.push_back(fun::VarStr("Kernel",GetKernelName(CSP.tkernel)));
  switch(CSP.tkernel){
    case KERNEL_Cubic:{
      const StKCubicCte &kc=CSP.kcubic;
      lines.push_back(fun::VarStr("  Cubic.a1" ,kc.a1));
      lines.push_back(fun::VarStr("  Cubic.aa" ,kc.aa));
      lines.push_back(fun::VarStr("  Cubic.a24",kc.a24));
      lines.push_back(fun::VarStr("  Cubic.c1" ,kc.c1));
      lines.push_back(fun::VarStr("  Cubic.c2" ,kc.c2));
      lines.push_back(fun::VarStr("  Cubic.d1" ,kc.d1));
      lines.push_back(fun::VarStr("  Cubic.od_wdeltap",kc.od_wdeltap));
    }break;
    case KERNEL_Wendland:{
      const StKWendlandCte &kc=CSP.kwend;
      lines.push_back(fun::VarStr("  Wendland.awen" ,kc.awen));
      lines.push_back(fun::VarStr("  Wendland.bwen" ,kc.bwen));
    }break;
    default: Run_ExceptioonFun("Kernel unknown.");
  }
}



}


