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
//:#############################################################################

/// \file JSphDtFixed.h \brief Declares the class \ref JSphDtFixed.

#ifndef _JSphDtFixed_
#define _JSphDtFixed_

#include "JObject.h"
#include "DualSphDef.h"
#include "JTimer.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

//##############################################################################
//# JSphDtFixed
//##############################################################################
/// \brief Manages the use of prefixed values of DT loaded from an input file.

class JSphDtFixed : protected JObject
{
protected:
  static const unsigned FILESIZEMAX=104857600; ///<Maximum file size (100mb).

  std::string File;
  unsigned Size;
  unsigned Count;
  unsigned Position;
  double *Times;
  double *Values;
  double DtError; //- max(DtFixed-DtVariable)

  double LastTimestepInput;  ///<Saves the last value used with GetDt().
  double LastDtInput;        ///<Saves the last value used with GetDt().
  double LastDtOutput;       ///<Saves the last value returned by GetDt().

  void Resize(unsigned size);

public:
  JSphDtFixed();
  ~JSphDtFixed();
  void Reset();
  unsigned GetAllocMemory()const;
  void LoadFile(std::string file);
  double GetDt(double timestep,double dtvar);
  double GetDtError(bool reset);
  std::string GetFile()const{ return(File); };
};

#endif


