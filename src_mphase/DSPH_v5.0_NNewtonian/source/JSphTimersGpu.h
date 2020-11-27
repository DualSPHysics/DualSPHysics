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

/// \file JSphTimersGpu.h \brief Measures time intervals during GPU execution.

#ifndef _JSphTimersGpu_
#define _JSphTimersGpu_

#ifdef DISABLE_TIMERS
  #define TmgStart(x,y) ;
  #define TmgStop(x,y) ;
#else
  #define TmgStart(x,y) _TmgStart(x,y)
  #define TmgStop(x,y) _TmgStop(x,y)
#endif

#include "JTimerCuda.h"

/// Structure with information of the timer and time value in CPU.
typedef struct{
  JTimerCuda timer;
  bool active;
  double time;
}StSphTimerGpu; 

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
  ,TMG_SuBoundCorr=18
  ,TMG_SuInOut=19
}CsTypeTimerGPU;
#define TMG_COUNT 20

typedef StSphTimerGpu TimersGpu[TMG_COUNT];

//==============================================================================
/// Returns the name of the timer.
//==============================================================================
inline const char* TmgGetName(CsTypeTimerGPU ct){
  switch(ct){
    case TMG_Init:              return("VA-Init");
    case TMG_NlLimits:          return("NL-Limits");
    case TMG_NlPreSort:         return("NL-PreSort");
    case TMG_NlRadixSort:       return("NL-RadixSort");
    case TMG_NlCellBegin:       return("NL-CellBegin");
    case TMG_NlSortData:        return("NL-SortData");
    case TMG_NlOutCheck:        return("NL-OutCheck");
    case TMG_CfPreForces:       return("CF-PreForces");
    case TMG_CfForces:          return("CF-Forces");
    case TMG_SuShifting:        return("SU-Shifting");
    case TMG_SuComputeStep:     return("SU-ComputeStep");
    case TMG_SuFloating:        return("SU-Floating");
    case TMG_SuMotion:          return("SU-Motion");
    case TMG_SuPeriodic:        return("SU-Periodic");
    case TMG_SuResizeNp:        return("SU-ResizeNp");
    case TMG_SuDownData:        return("SU-DownData");
    case TMG_SuSavePart:        return("SU-SavePart");
    case TMG_SuChrono:          return("SU-Chrono");
    case TMG_SuBoundCorr:       return("SU-BoundCorr");
    case TMG_SuInOut:           return("SU-InOut");
  }
  return("???");
}

//==============================================================================
/// Returns the number of timers.
//==============================================================================
inline unsigned TmgGetCount(){ return(TMG_COUNT); }

//==============================================================================
/// Creates timers to measure time intervals.
//==============================================================================
inline void TmgCreation(TimersGpu vtimer,bool active){
  for(unsigned c=0;c<TMG_COUNT;c++){
    StSphTimerGpu* t=vtimer+c;
    t->timer.Reset();
    t->active=active;
    t->time=0;
  }
}

//==============================================================================
/// Destroys the timers used to measure times.
//==============================================================================
inline void TmgDestruction(TimersGpu vtimer){ TmgCreation(vtimer,false); }

//==============================================================================
/// Marks start of timer.
//==============================================================================
inline void _TmgStart(TimersGpu vtimer,CsTypeTimerGPU ct){ if(vtimer[ct].active)vtimer[ct].timer.Start(); }

//==============================================================================
/// Marks end of timer and accumulates time.
//==============================================================================
inline void _TmgStop(TimersGpu vtimer,CsTypeTimerGPU ct){
  StSphTimerGpu* t=vtimer+unsigned(ct);
  if(t->active){
    t->timer.Stop();
    t->time+=t->timer.GetElapsedTimeD();
  }
}

//==============================================================================
/// Initialises the time accumulated by all timers.
//==============================================================================
inline void TmgResetValues(TimersGpu vtimer){
  for(unsigned ct=0;ct<TMG_COUNT;ct++)if(vtimer[ct].active)vtimer[ct].time=0; 
}  

//==============================================================================
/// Increases the time accumulated manually.
//==============================================================================
inline void TmgIncreaseValue(TimersGpu vtimer,CsTypeTimerGPU ct,double value){ if(vtimer[ct].active)vtimer[ct].time+=value; }  

//==============================================================================
/// Returns the time accumulated by the timer.
//==============================================================================
inline float TmgGetValue(const TimersGpu vtimer,CsTypeTimerGPU ct){ return(float(vtimer[ct].time)); }

//==============================================================================
/// Returns the time accumulated by the timer.
//==============================================================================
inline double TmgGetValueD(const TimersGpu vtimer,CsTypeTimerGPU ct){ return(vtimer[ct].time); }

//==============================================================================
/// Returns pointer pointing at time accumulated by the timer.
//==============================================================================
inline const double* TmgGetPtrValue(const TimersGpu vtimer,CsTypeTimerGPU ct){ return(&(vtimer[ct].time)); }

//==============================================================================
/// Indicates whether a timer is active or not.
//==============================================================================
inline bool TmgIsActive(const TimersGpu vtimer,CsTypeTimerGPU ct){ return(unsigned(ct)<TmgGetCount()&&vtimer[ct].active); }

//==============================================================================
/// Activates a given timer.
//==============================================================================
inline void TmgActive(TimersGpu vtimer,CsTypeTimerGPU ct,bool active){ if(unsigned(ct)<TmgGetCount())vtimer[ct].active=active; }

//==============================================================================
/// Activates o deactivates all timers.
//==============================================================================
inline void TmgActiveAll(TimersGpu vtimer,bool active){ for(unsigned c=0;c<TMG_COUNT;c++)vtimer[c].active=active; }

#endif


