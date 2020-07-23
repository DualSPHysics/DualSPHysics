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

/// \file JObjectGpu.cpp \brief Implements the class \ref JObjectGpu.

#include "JObjectGpu.h"
#include "JException.h"
#include "Functions.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

//##############################################################################
//# JObjectGpu
//##############################################################################
////==============================================================================
///// Throws exception related to a CUDA error from a static method.
////==============================================================================
//void JObjectGpu::RunExceptioonCudaStatic(const std::string &srcfile,int srcline
//  ,const std::string &method
//  ,cudaError_t cuerr,std::string msg)
//{
//  msg=msg+fun::PrintStr(" (CUDA error %d (%s)).\n",cuerr,cudaGetErrorString(cuerr));
//  throw JException(srcfile,srcline,"JObjectGpu",method,msg,"");
//}
////==============================================================================
///// Checks CUDA error and throws exception from a static method.
////==============================================================================
//void JObjectGpu::CheckCudaErroorStatic(const std::string &srcfile,int srcline
//  ,const std::string &method,std::string msg)
//{
//  cudaError_t cuerr=cudaGetLastError();
//  if(cuerr!=cudaSuccess)RunExceptioonCudaStatic(srcfile,srcline,method,cuerr,msg);
//}

//==============================================================================
/// Throws exception related to a CUDA error.
//==============================================================================
void JObjectGpu::RunExceptioonCuda(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,cudaError_t cuerr,std::string msg)const
{
  msg=msg+fun::PrintStr(" (CUDA error %d (%s)).\n",cuerr,cudaGetErrorString(cuerr));
  throw JException(srcfile,srcline,classname,method,msg,"");
}

//==============================================================================
/// Checks CUDA error and throws exception.
/// Comprueba error de CUDA y lanza excepcion si lo hubiera.
//==============================================================================
void JObjectGpu::CheckCudaErroor(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,std::string msg)const
{
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess)RunExceptioonCuda(srcfile,srcline,classname,method,cuerr,msg);
}


