//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2021 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsGpuInfo.h \brief Declares the class \ref JDsGpuInfo.

#ifndef _JDsGpuInfo_
#define _JDsGpuInfo_

#include "JObject.h"
#include "TypesDef.h"

class JLog2;

//##############################################################################
//# JDsGpuInfo
//##############################################################################
/// \brief Manages general information of selected GPU.

class JDsGpuInfo : protected JObject
{
protected:
  int Ngpus;            ///<Number of available GPUs (def=-1).
  int GpuId;            ///<GPU Selection.
  bool AutoSelect;      ///<By automatic selection.
  std::string Name;     ///<Name of the selected GPU.
  size_t GlobalMem;     ///<Size of global memory in bytes.
  unsigned SharedMem;   ///<Size of shared memory for each block in bytes.
  unsigned Compute;     ///<Compute capability: 10,11,12,20... 
  std::string Hardware; ///<Hardware description in short text.

public:
  JDsGpuInfo();
  ~JDsGpuInfo();

  void Reset();

  static int ShowGpusInfo(JLog2 *log);
  int GetNgpus();

  int SelectGpu(int gpuid);
  void ShowSelectGpusInfo(JLog2 *log);

  int         GetGpuId   ()const{ return(GpuId);    }
  std::string GetName    ()const{ return(Name);     }
  std::string GetHardware()const{ return(Hardware); }
};

#endif


