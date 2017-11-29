//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Clase para Simplificar la creacion de ficheros CSV en debug. (23-03-2014)
//:# - Ahora al ejecutar Save() vacia el buffer. (06-02-2015)
//:# - Funciones AddValue() para datos triples. (26-05-2016)
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:# - New attribute CsvSepComa to configure separator in CSV files. (24-10-2017)
//:#############################################################################

/// \file JSaveCsv.h \brief Declares the class \ref JSaveCsv.

#ifndef _JSaveCsv_
#define _JSaveCsv_

#include "TypesDef.h"
#include "JObject.h"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

//##############################################################################
//# JSaveCsv
//##############################################################################
/// \brief Creates CSV files in a simple way.

class JSaveCsv : protected JObject
{
private:
  const bool CsvSepComa;   ///<Separator character in CSV files (0=semicolon, 1=coma).
  char CsvSep[2];
  std::string FileName;
  bool App;

  std::fstream *Pf;
  bool FileError;

  std::string Head;
  std::string Data;
  void Save(const std::string &tx);

public:
  JSaveCsv(std::string fname,bool app,bool csvsepcoma);
  ~JSaveCsv();
  void Reset();

  std::string GetFileName()const{ return(FileName); }

  void SaveData(bool closefile=false);

  void AddHead(const std::string &head,bool endl=true);
  void AddData(const std::string &data,bool endl=true);

  void AddValuesf(const char *format,...);
  void AddValue(const std::string &v);
  void AddValue(unsigned v);
  void AddValue(int v);
  void AddValue(float v);
  void AddValue(double v);
  void AddValue(tint3    v){ AddValue(v.x); AddValue(v.y); AddValue(v.z); }
  void AddValue(tuint3   v){ AddValue(v.x); AddValue(v.y); AddValue(v.z); }
  void AddValue(tfloat3  v){ AddValue(v.x); AddValue(v.y); AddValue(v.z); }
  void AddValue(tdouble3 v){ AddValue(v.x); AddValue(v.y); AddValue(v.z); }

  void AddEndl();
};

#endif


