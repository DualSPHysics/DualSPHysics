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

/// \file JDsTimersGpu.h \brief Declares the class \ref JDsTimersGpu.

#ifndef _JDsTimersGpu_
#define _JDsTimersGpu_

#include "JDsTimers.h"
#include <cuda_runtime_api.h>

/// List of possible timers to define for single GPU executions.
typedef enum{
   TMG_Init=0
  ,TMG_NlLimits=1
  ,TMG_NlPreSort=2
  ,TMG_NlRadixSort=3
  ,TMG_NlCellBegin=4
  ,TMG_NlSortData=5
  ,TMG_NlOutCheck=6
  ,TMG_CfPreForces=7
  ,TMG_CfForces=8
  ,TMG_SuShifting=9
  ,TMG_SuComputeStep=10
  ,TMG_SuFloating=11
  ,TMG_SuMotion=12
  ,TMG_SuPeriodic=13
  ,TMG_SuResizeNp=14
  ,TMG_SuDownData=15
  ,TMG_SuSavePart=16
  ,TMG_SuChrono=17
  ,TMG_SuMoorings=18
  ,TMG_SuInOut=19
  ,TMG_SuGauges=20
}TpTimersGPU;

//##############################################################################
//# JDsTimersGpu
//##############################################################################
/// \brief Measures time intervals during GPU execution.

class JDsTimersGpu : public JDsTimers
{

public:
  //==============================================================================
  /// Constructor.
  //==============================================================================
  JDsTimersGpu():JDsTimers("JDsTimersGpu"){ }

  //==============================================================================
  /// Configures timers for CPU executions.
  //==============================================================================
  void Config(bool svtimers){
    Reset();
    SvTimers=svtimers;
    Add(TMG_Init              ,"VA-Init"       ,0,SvTimers);
    Add(TMG_NlLimits          ,"NL-Limits"     ,0,SvTimers);
    Add(TMG_NlPreSort         ,"NL-PreSort"    ,0,SvTimers);
    Add(TMG_NlRadixSort       ,"NL-RadixSort"  ,0,SvTimers);
    Add(TMG_NlCellBegin       ,"NL-CellBegin"  ,0,SvTimers);
    Add(TMG_NlSortData        ,"NL-SortData"   ,0,SvTimers);
    Add(TMG_NlOutCheck        ,"NL-OutCheck"   ,0,SvTimers);
    Add(TMG_CfPreForces       ,"CF-PreForces"  ,0,SvTimers);
    Add(TMG_CfForces          ,"CF-Forces"     ,0,SvTimers);
    Add(TMG_SuShifting        ,"SU-Shifting"   ,0,SvTimers);
    Add(TMG_SuComputeStep     ,"SU-ComputeStep",0,SvTimers);
    Add(TMG_SuFloating        ,"SU-Floating"   ,0,SvTimers);
    Add(TMG_SuMotion          ,"SU-Motion"     ,0,SvTimers);
    Add(TMG_SuPeriodic        ,"SU-Periodic"   ,0,SvTimers);
    Add(TMG_SuResizeNp        ,"SU-ResizeNp"   ,0,SvTimers);
    Add(TMG_SuDownData        ,"SU-DownData"   ,0,SvTimers);
    Add(TMG_SuSavePart        ,"SU-SavePart"   ,0,SvTimers);
    Add(TMG_SuChrono          ,"SU-Chrono"     ,0,SvTimers);
    Add(TMG_SuMoorings        ,"SU-Moorings"   ,0,SvTimers);
    Add(TMG_SuInOut           ,"SU-InOut"      ,0,SvTimers);
    Add(TMG_SuGauges          ,"SU-Gauges"     ,0,SvTimers);
  }
  
  //==============================================================================
  /// Add timer to list.
  //==============================================================================
  void Add(TpTimersGPU ct,std::string name,unsigned level,bool active){ 
    AddTimer(unsigned(ct),name,level,active); 
  }

  //==============================================================================
  /// Marks start of timer.
  //==============================================================================
  inline void TmStart(TpTimersGPU ct,bool synchronize){
    if(List[ct].active){
      if(synchronize)cudaDeviceSynchronize();
      List[ct].timer.Start();
    }
  }

  //==============================================================================
  /// Marks end of timer and accumulates time.
  //==============================================================================
  inline void TmStop(TpTimersGPU ct,bool synchronize){ 
    StDsTimer* t=List+unsigned(ct);
    if(t->active){
      if(synchronize)cudaDeviceSynchronize();
      t->timer.Stop();
      t->time+=t->timer.GetElapsedTimeD();
    }
  }
};

#endif


