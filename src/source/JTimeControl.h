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

/// \file JTimeControl.h \brief Declares the class \ref JTimeControl.

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Implementacion de clase para controlar tiempo de ejecucion (en segundos) y
//:#   calcular tiempo estimado de finalizacion con un overhead minimo. Se calculo  
//:#   un overhead del 2.1% en bucles de iteraciones muy ligeras. (25-04-2016)
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:#############################################################################

#ifndef _JTimeControl_
#define _JTimeControl_

#include "JObject.h"
#include "TypesDef.h"
#include "JTimer.h"
#include <string>
#include <cstdlib>
#include <vector>

//##############################################################################
//# JTimeControl
//##############################################################################
/// \brief Defines class to manage information of runtime in long processes.

class JTimeControl : protected JObject
{
protected:
  bool Active;

  JTimer Timer;

  double NextTime;

  bool Periodic;        ///<For periodic intervals of TimeOut. | Para intervalos periodicos de TimeOut.
  double FirstTime;     ///<First time to evaluate (when it is less than timeout and not zero). | Primer instante a evaluar (cuando es menor que timeout y no es cero).
  double TimeOut;       ///<Periodic interval duration (in seconds). | Duracion de intervalo periodico (en segundos).

  unsigned TimeOutNum;  //<Number of intervals processed with Periodic or Times. | Numero de intervalos procesados con Periodic o Times.

  unsigned TimesSize;        ///<Number of elements of Times. | Numero de elementos de Times. 
  std::vector<double> Times; ///<List of times when it is not periodic. | Lista de tiempos cuando no es Periodic.

  double IteStart;  ///<Instant to start (IteNum is initialized). | Instante en que empieza (se inicializa IteNum).
  unsigned IteNum;  ///<Number of iterations. | Numero de iteraciones.
  unsigned NextIte; ///<Next iteration to check. | Siguiente iteracion a comprobar.

  double LastTime;


  void ConfigPeriodic(double tfirst,double tout);
  void ConfigTimes(unsigned ntimes,const double *vtimes);
  void ConfigTimes(std::string times);
  void PrepareTimes();


  double CalcNextTime();
  bool CheckRealTime();

public:
  JTimeControl(double tout);
  JTimeControl(double tfirst,double tout);
  JTimeControl(unsigned ntimes,const double *vtimes);
  JTimeControl(const std::string &times);
  //~JTimeControl();

  void Reset();

  bool CheckTime(){
    if(Active){
      IteNum++;
      if(IteNum>=NextIte)return(CheckRealTime());
    }
    return(false);
  }

  double GetLastTime()const{ return(LastTime); }
  std::string GetInfoFinish(double done);
};

#endif


