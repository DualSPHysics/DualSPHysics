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
//:# - Clases para calcular promedios simples, moviles y ponderados. (20-02-2016)
//:# - Se elimina codigo de JMeanMoving por falta de uso. (20-11-2020)
//:#############################################################################

/// \file JMeanValues.h \brief Declares the class \ref JMeanValue and class \ref JMeanMoving.

#ifndef _JMeanValues_
#define _JMeanValues_

#include "JObject.h"
#include "TypesDef.h"
#include <cfloat>

//##############################################################################
//# JMeanValue
//##############################################################################
/// \brief Calculates the average value of a sequence of values.

class JMeanValue
{
public:
  double Max;
  double Min;
  double Mean;
  ullong Values;

public:
  JMeanValue():Max(-DBL_MAX),Min(DBL_MAX),Mean(0),Values(0){ }
  void Reset(){ Max=-DBL_MAX; Min=DBL_MAX; Mean=0; Values=0; }
  void AddValue(double v){ 
    Max=(Max<v? v: Max);
    Min=(Min>v? v: Min);
    Mean=(Mean*Values+v)/(Values+1); Values++; 
  }
  double GetMax()const{ return(Max); }
  double GetMin()const{ return(Min); }
  double GetMean()const{ return(Mean); }
  ullong GetValues()const{ return(Values); }
};


#endif


