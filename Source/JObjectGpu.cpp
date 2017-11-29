//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "Functions.h"
#include <cstdio>
#include <cstdlib>
#include <ctime>

//==============================================================================
/// Throws exception due to a Cuda error.
/// \param method Name of the method that throws an exception.
/// \param msg Text of the exception.
/// \param error Code of the Cuda error.
/// \throw JException 
//==============================================================================
void JObjectGpu::RunExceptionCuda(const std::string &method,const std::string &msg,cudaError_t error)const{
  RunException(method,fun::PrintStr("%s (CUDA error: %s).\n",msg.c_str(),cudaGetErrorString(error)));
}

//==============================================================================
/// Checks error if there is a Cuda error and throws exception.
/// \param method Name of the method that throws an exception.
/// \param msg Text of the exception.
//==============================================================================
void JObjectGpu::CheckCudaError(const std::string &method,const std::string &msg)const{
  cudaError_t err=cudaGetLastError();
  if(err!=cudaSuccess)RunExceptionCuda(method,msg,err);
}


