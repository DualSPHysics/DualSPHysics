//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSpVtkData.h \brief Declares the class \ref JSpVtkData.

#ifndef _JSpVtkData_
#define _JSpVtkData_

#include "TypesDef.h"
#include "JObject.h"
#include <string>

class JDataArrays;

//##############################################################################
//# JSpVtkData
//##############################################################################
/// \brief Creates VTK file with data arrays.

class JSpVtkData : protected JObject
{
public:
  JSpVtkData();
  ~JSpVtkData();
  void SaveVtk(std::string file,const JDataArrays& arrays,std::string posfield);
  static void Save(std::string file,const JDataArrays& arrays,std::string posfield);
};

#endif

