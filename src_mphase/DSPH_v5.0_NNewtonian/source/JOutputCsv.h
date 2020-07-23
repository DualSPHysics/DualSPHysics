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
//:# - Clase para grabar arrays de datos en distintos formatos. (05-07-2019)
//:# - Codigo actualizado para trabajar con ultima version de JDataArrays. (10-12-2019)
//:# - CsvSepComa se define en el constructor y por defecto es false. (27-12-2019)
//:# - CreatPath se define en el constructor y por defecto es true. (27-12-2019)
//:#############################################################################

/// \file JOutputCsv.h \brief Declares the class \ref JOutputCsv.

#ifndef _JOutputCsv_
#define _JOutputCsv_

#include "TypesDef.h"
#include "JObject.h"
#include "JDataArrays.h"
#include <string>
#include <cstring>
#include <string>
#include <vector>

//##############################################################################
//# JOutputCsv
//##############################################################################
/// \brief Saves data arrays in CSV files.

class JOutputCsv : protected JObject
{
public:
  //-Output configuration.
  bool CreatPath;   ///<Creates full path for output files (true by default).
  bool CsvSepComa;  ///<Separator character in CSV files (0=semicolon, 1=coma).

protected:
  std::string FileName; ///<Last file generated.
  template<typename T> void CalculateStatsArray1(unsigned ndata,T *data
    ,double &valmin,double &valmax,double &valmean)const;
  template<typename T> void CalculateStatsArray3(unsigned ndata,T *data
    ,tdouble4 &valmin,tdouble4 &valmax,tdouble4 &valmean)const;

public:
  JOutputCsv(bool csvsepcoma=false,bool createpath=true);
  ~JOutputCsv();
  void Reset();

  void SaveCsv(std::string fname,const JDataArrays &arrays,std::string head="");
  //void SaveStatsCsv(std::string fname,bool create,int part,double timestep,const JDataArrays &arrays,std::string head="");

};




#endif

