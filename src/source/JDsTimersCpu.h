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

/// \file JDsTimersCpu.h \brief Declares the class \ref JDsTimersCpu.

#ifndef _JDsTimersCpu_
#define _JDsTimersCpu_

#include "JDsTimers.h"

/// List of possible timers to define for CPU executions.
typedef enum{
   TMC_Init=0
  ,TMC_NlLimits=1
  ,TMC_NlMakeSort=2
  ,TMC_NlSortData=3
  ,TMC_NlOutCheck=4
  ,TMC_CfPreForces=5
  ,TMC_CfForces=6
  ,TMC_SuShifting=7
  ,TMC_SuComputeStep=8
  ,TMC_SuFloating=9
  ,TMC_SuMotion=10
  ,TMC_SuPeriodic=11
  ,TMC_SuResizeNp=12
  ,TMC_SuSavePart=13
  ,TMC_SuChrono=14
  ,TMC_SuMoorings=15
  ,TMC_SuInOut=16
  ,TMC_SuGauges=17
}TpTimersCPU;

//##############################################################################
//# JDsTimersCpu
//##############################################################################
/// \brief Measures time intervals during CPU execution.

class JDsTimersCpu : public JDsTimers
{

public:
  //==============================================================================
  /// Constructor.
  //==============================================================================
  JDsTimersCpu():JDsTimers("JDsTimersCpu"){ }

  //==============================================================================
  /// Configures timers for CPU executions.
  //==============================================================================
  void Config(bool svtimers){
    Reset();
    SvTimers=svtimers;
    Add(TMC_Init         ,"VA-Init"       ,0,SvTimers);
    Add(TMC_NlLimits     ,"NL-Limits"     ,0,SvTimers);
    Add(TMC_NlMakeSort   ,"NL-MakeSort"   ,0,SvTimers);
    Add(TMC_NlSortData   ,"NL-SortData"   ,0,SvTimers);
    Add(TMC_NlOutCheck   ,"NL-OutCheck"   ,0,SvTimers);
    Add(TMC_CfPreForces  ,"CF-PreForces"  ,0,SvTimers);
    Add(TMC_CfForces     ,"CF-Forces"     ,0,SvTimers);
    Add(TMC_SuShifting   ,"SU-Shifting"   ,0,SvTimers);
    Add(TMC_SuComputeStep,"SU-ComputeStep",0,SvTimers);
    Add(TMC_SuFloating   ,"SU-Floating"   ,0,SvTimers);
    Add(TMC_SuMotion     ,"SU-Motion"     ,0,SvTimers);
    Add(TMC_SuPeriodic   ,"SU-Periodic"   ,0,SvTimers);
    Add(TMC_SuResizeNp   ,"SU-ResizeNp"   ,0,SvTimers);
    Add(TMC_SuSavePart   ,"SU-SavePart"   ,0,SvTimers);
    Add(TMC_SuChrono     ,"SU-Chrono"     ,0,SvTimers);
    Add(TMC_SuMoorings   ,"SU-Moorings"   ,0,SvTimers);
    Add(TMC_SuInOut      ,"SU-InOut"      ,0,SvTimers);
    Add(TMC_SuGauges     ,"SU-Gauges"     ,0,SvTimers);
  }
  
  //==============================================================================
  /// Add timer to list.
  //==============================================================================
  void Add(TpTimersCPU ct,std::string name,unsigned level,bool active){ 
    AddTimer(unsigned(ct),name,level,active); 
  }

  //==============================================================================
  /// Marks start of timer.
  //==============================================================================
  inline void TmStart(TpTimersCPU ct){ TimerStart(unsigned(ct)); }

  //==============================================================================
  /// Marks end of timer and accumulates time.
  //==============================================================================
  inline void TmStop(TpTimersCPU ct){ TimerStop(unsigned(ct)); }
};

#endif


