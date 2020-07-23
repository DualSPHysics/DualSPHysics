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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestiona el uso de un dt fijo (en milisegundos) a partir de los valores 
//:#   para determinados instantes en segundos, interpolando los valores 
//:#   intermedios. (10-11-2012)
//:# - Los datos float pasaron a double. (28-11-2013) 
//:# - GetNextTime() guarda entrada y salida para evitar calculos con llamadas 
//:#   consecutivas iguales. (02-09-2019)
//:# - Cambio de nombre de J.SphDtFixed a J.DsFixedDt. (28-06-2020)
//:# - Permite establecer un dt fijo en el constructor. (18-07-2020)
//:#############################################################################

/// \file JDsFixedDt.h \brief Declares the class \ref JDsFixedDt.

#ifndef _JDsFixedDt_
#define _JDsFixedDt_

#include "JObject.h"
#include "DualSphDef.h"
#include "JTimer.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

//##############################################################################
//# JDsFixedDt
//##############################################################################
/// \brief Manages the use of prefixed values of DT loaded from an input file.

class JDsFixedDt : protected JObject
{
protected:
  static const unsigned FILESIZEMAX=104857600; ///<Maximum file size (100mb).
  double FixedValue;

  std::string File;
  unsigned Size;
  unsigned Count;
  unsigned Position;
  double *Times;
  double *Values;
  double DtError; //- max(FixedDt-DtVariable)

  double LastTimestepInput;  ///<Saves the last value used with GetDt().
  double LastDtInput;        ///<Saves the last value used with GetDt().
  double LastDtOutput;       ///<Saves the last value returned by GetDt().

  void Resize(unsigned size);

public:
  JDsFixedDt(double fixedvalue=0);
  ~JDsFixedDt();
  void Reset();
  unsigned GetAllocMemory()const;
  void LoadFile(std::string file);
  double GetDt(double timestep,double dtvar);
  double GetDtError(bool reset);
  std::string GetFile()const{ return(File); };
  double GetFixedValue()const{ return(FixedValue); };
};

#endif


