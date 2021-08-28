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

/// \file JDsTimers.h \brief Declares the class \ref JDsTimers.

#ifndef _JDsTimers_
#define _JDsTimers_

#include "JObject.h"
#include "TypesDef.h"
#include "JTimer.h"

class JLog2;

/// Structure with information of the timer and time value in CPU.
typedef struct{
  bool active;
  JTimer timer;
  double time;
  unsigned level;
  std::string name;
}StDsTimer; 

//##############################################################################
//# JDsTimers
//##############################################################################
/// \brief Measures time intervals during tasks execution.

class JDsTimers : protected JObject
{
public:
  static const int TIMERSIZE=50;

protected:
  StDsTimer *List;
  unsigned CtMax;
  bool SvTimers;

  void ResetTimer(unsigned c);
  void AddTimer(unsigned c,std::string name,unsigned level,bool active);
  std::string TimerToText(unsigned c,unsigned maxlen)const;

  /// Marks start of timer.
  inline void TimerStart(unsigned c){ if(List[c].active)List[c].timer.Start(); }

  /// Marks end of timer and accumulates time.
  inline void TimerStop(unsigned c){
    StDsTimer* t=List+c;
    if(t->active){
      t->timer.Stop();
      t->time+=t->timer.GetElapsedTimeD();
    }
  }


public:
  JDsTimers(std::string classname);
  ~JDsTimers();

  void Reset();
  void ResetTimes();
  
  void ShowTimes(std::string title,JLog2 *log,bool onlyfile=false)const;
  void GetTimersInfo(std::string &hinfo,std::string &dinfo)const;
};

#endif


