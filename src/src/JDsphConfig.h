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
//:# - Carga configuracion del fichero DsphConfig.xml. (23-10-2017)
//:#############################################################################

/// \file JDsphConfig.h \brief Declares the class \ref JDsphConfig.

#ifndef _JDsphConfig_
#define _JDsphConfig_


#include "JObject.h"
#include "TypesDef.h"
#include <string>

//##############################################################################
//# JDsphConfig
//##############################################################################
/// \brief Loads configuration from DphConfig.xml

class JDsphConfig : protected JObject
{
protected:
  std::string FileCfg;

  int CsvSeparator; //-Separator character in CSV files (-1=undefined, 0=semicolon, 1=coma)

public:
  JDsphConfig();
  ~JDsphConfig();
  void Reset();
  void Init(std::string path);

  std::string GetFileCfg()const{ return(FileCfg); }

  int GetCsvSeparator()const{ return(CsvSeparator); }
};

#endif


