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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

/// \file JObjectGpu.h \brief Declares the class \ref JObjectGpu.

#ifndef _JObjectGpu_
#define _JObjectGpu_

#include "JObject.h"
#include "RunExceptionGpuDef.h"
#include <cuda_runtime_api.h>
#include <string>


//##############################################################################
//# JObjectGpu
//##############################################################################
/// \brief Defines objects with methods that throw exceptions for tasks on the GPU.
class JObjectGpu : protected JObject
{
protected:
  //static void RunExceptioonCudaStatic(const std::string &srcfile,int srcline
  //  ,const std::string &method
  //  ,cudaError_t cuerr,std::string msg);
  //static void CheckCudaErroorStatic(const std::string &srcfile,int srcline
  //  ,const std::string &method
  //  ,std::string msg);

  void RunExceptioonCuda(const std::string &srcfile,int srcline
    ,const std::string &classname,const std::string &method
    ,cudaError_t cuerr,std::string msg)const;
  void CheckCudaErroor(const std::string &srcfile,int srcline
    ,const std::string &classname,const std::string &method
    ,std::string msg)const;
public:  
  JObjectGpu(){ ClassName="JObjectGpu"; } ///<Constructor.
};

#endif


