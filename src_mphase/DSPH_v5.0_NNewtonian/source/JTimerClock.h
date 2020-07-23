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
//:# - Implementacion de una clase para medir intervalos de tiempo con la
//:#   precision (~ milisegundos) que proporciona clock(). (10-01-2011)
//:#############################################################################

/// \file JTimerClock.h \brief Declares the class \ref JTimerClock.

#ifndef _JTimerClock_
#define _JTimerClock_

#include <ctime>


//##############################################################################
//# JTimerClock
//##############################################################################
/// \brief Defines a class to measure time intervals with precision of clock().

class JTimerClock
{
private:
  bool Stopped;
  clock_t CounterIni,CounterEnd;

public:
  JTimerClock(){ Reset(); }
  void Reset(){ Stopped=false; CounterIni=0; CounterEnd=0; }
  void Start(){ Stopped=false; CounterIni=clock(); }
  void Stop(){ CounterEnd=clock(); Stopped=true; }
  //-Returns time in milliseconds.
  //-Devuelve tiempo en milisegundos.
  float GetElapsedTimeF(){ return((float(Stopped? CounterEnd-CounterIni: 0)*float(1000))/float(CLOCKS_PER_SEC)); }
  double GetElapsedTimeD(){ return((double(Stopped? CounterEnd-CounterIni: 0)*double(1000))/double(CLOCKS_PER_SEC)); }
};

#endif


