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

/// \file JSphTimersCpu.h \brief Measures time intervals during CPU execution.

#ifndef _JSphTimersCpu_
#define _JSphTimersCpu_

#ifdef DISABLE_TIMERS
  #define TmcStart(x,y) ;
  #define TmcStop(x,y) ;
#else
  #define TmcStart(x,y) _TmcStart(x,y)
  #define TmcStop(x,y) _TmcStop(x,y)
#endif

#include "JTimer.h" //"JTimerClock.h"

/// Structure with information of the timer and time value in CPU.
typedef struct{
  JTimer timer; //JTimerClock timer;
  bool active;
  double time;
}StSphTimerCpu; 

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
  ,TMC_SuBoundCorr=15
  ,TMC_SuInOut=16
}CsTypeTimerCPU;
#define TMC_COUNT 17

typedef StSphTimerCpu TimersCpu[TMC_COUNT];

//==============================================================================
/// Returns the name of the timer.
//==============================================================================
inline const char* TmcGetName(CsTypeTimerCPU ct){
  switch(ct){
    case TMC_Init:              return("VA-Init");
    case TMC_NlLimits:          return("NL-Limits");
    case TMC_NlMakeSort:        return("NL-MakeSort");
    case TMC_NlSortData:        return("NL-SortData");
    case TMC_NlOutCheck:        return("NL-OutCheck");
    case TMC_CfPreForces:       return("CF-PreForces");
    case TMC_CfForces:          return("CF-Forces");
    case TMC_SuShifting:        return("SU-Shifting");
    case TMC_SuComputeStep:     return("SU-ComputeStep");
    case TMC_SuFloating:        return("SU-Floating");
    case TMC_SuMotion:          return("SU-Motion");
    case TMC_SuPeriodic:        return("SU-Periodic");
    case TMC_SuResizeNp:        return("SU-ResizeNp");
    case TMC_SuSavePart:        return("SU-SavePart");
    case TMC_SuChrono:          return("SU-Chrono");
    case TMC_SuBoundCorr:       return("SU-BoundCorr");
    case TMC_SuInOut:           return("SU-InOut");
  }
  return("???");
}

//==============================================================================
/// Returns the number of timers.
//==============================================================================
inline unsigned TmcGetCount(){ return(TMC_COUNT); }

//==============================================================================
/// Creates timers to measure time intervals.
//==============================================================================
inline void TmcCreation(TimersCpu vtimer,bool active){
  for(unsigned c=0;c<TMC_COUNT;c++){
    StSphTimerCpu* t=vtimer+c;
    t->timer.Reset();
    t->active=active;
    t->time=0;
  }
}

//==============================================================================
/// Destroys the timers used to measure times.
//==============================================================================
inline void TmcDestruction(TimersCpu vtimer){ TmcCreation(vtimer,false); }

//==============================================================================
/// Marks start of timer.
//==============================================================================
inline void _TmcStart(TimersCpu vtimer,CsTypeTimerCPU ct){ if(vtimer[ct].active)vtimer[ct].timer.Start(); }

//==============================================================================
/// Marks end of timer and accumulates time.
//==============================================================================
inline void _TmcStop(TimersCpu vtimer,CsTypeTimerCPU ct){
  StSphTimerCpu* t=vtimer+unsigned(ct);
  if(t->active){
    t->timer.Stop();
    t->time+=t->timer.GetElapsedTimeD();
  }
}

//==============================================================================
/// Initialises the time accumulated by all timers.
//==============================================================================
inline void TmcResetValues(TimersCpu vtimer){
  for(unsigned ct=0;ct<TMC_COUNT;ct++)if(vtimer[ct].active)vtimer[ct].time=0; 
}  

//==============================================================================
/// Increases the time accumulated manually.
//==============================================================================
inline void TmcIncreaseValue(TimersCpu vtimer,CsTypeTimerCPU ct,double value){ if(vtimer[ct].active)vtimer[ct].time+=value; }  

//==============================================================================
/// Returns the time accumulated by the timer.
//==============================================================================
inline float TmcGetValue(const TimersCpu vtimer,CsTypeTimerCPU ct){ return(float(vtimer[ct].time)); }

//==============================================================================
/// Returns the time accumulated by the timer.
//==============================================================================
inline double TmcGetValueD(const TimersCpu vtimer,CsTypeTimerCPU ct){ return(vtimer[ct].time); }

//==============================================================================
/// Returns pointer pointing at time accumulated by the timer.
//==============================================================================
inline const double* TmcGetPtrValue(const TimersCpu vtimer,CsTypeTimerCPU ct){ return(&(vtimer[ct].time)); }

//==============================================================================
/// Indicates whether a timer is active or not.
//==============================================================================
inline bool TmcIsActive(const TimersCpu vtimer,CsTypeTimerCPU ct){ return(unsigned(ct)<TmcGetCount()&&vtimer[ct].active); }

//==============================================================================
/// Activates a given timer.
//==============================================================================
inline void TmcActive(TimersCpu vtimer,CsTypeTimerCPU ct,bool active){ if(unsigned(ct)<TmcGetCount())vtimer[ct].active=active; }

//==============================================================================
/// Activates o deactivates all timers.
//==============================================================================
inline void TmcActiveAll(TimersCpu vtimer,bool active){ for(unsigned c=0;c<TMC_COUNT;c++)vtimer[c].active=active; }

#endif


