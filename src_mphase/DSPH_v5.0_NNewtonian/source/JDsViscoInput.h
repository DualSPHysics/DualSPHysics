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
//:# - Gestiona el uso de un valor de viscosidad variable a partir de los valores 
//:#   para determinados instantes en segundos, interpolando los valores 
//:#   intermedios. (12-04-2013)
//:# - Usa JReadDatafile para cargar datos de fichero. (17-12-2015)
//:# - GetVisco() guarda entrada y salida para evitar calculos con llamadas 
//:#   consecutivas iguales. (08-09-2019)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Cambio de nombre de J.SphVisco a J.DsViscoInput. (28-06-2020)
//:#############################################################################

/// \file JDsViscoInput.h \brief Declares the class \ref JDsViscoInput.

#ifndef _JDsViscoInput_
#define _JDsViscoInput_

#include "JObject.h"
#include "DualSphDef.h"
#include "JTimer.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

//##############################################################################
//# JDsViscoInput
//##############################################################################
/// \brief Manages the use of viscosity values from an input file.

class JDsViscoInput : protected JObject
{
protected:
  static const unsigned FILESIZEMAX=104857600; ///<Maximum file size (100mb).

  std::string File;
  unsigned Size;
  unsigned Count;
  unsigned Position;
  float *Times;
  float *Values;

  float LastTimestepInput;  ///<Saves the last value used with GetVisco().
  float LastViscoOutput;    ///<Saves the last value returned by GetVisco().

  void Resize(unsigned size);

public:
  JDsViscoInput();
  ~JDsViscoInput();
  void Reset();
  unsigned GetAllocMemory()const;
  void LoadFile(std::string file);
  float GetVisco(float timestep);
  std::string GetFile()const{ return(File); };
};

#endif


