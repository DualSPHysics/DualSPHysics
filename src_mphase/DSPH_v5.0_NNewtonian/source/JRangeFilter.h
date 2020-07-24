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
//:# - Clase creada para facilitar el filtrado de valores entre una serie de 
//:#   rangos. (17-12-2010)
//:# - Cambio de nombre de JFilterValues a JRangeFilter. (27-12-2013)
//:# - Los metodos GetNextValue() y GetFirstValue devuelven UINT_MAX en caso de 
//:#   no tener otro valor para devolver. (27-12-2013)
//:# - Permite rangos con incremento configurable Ej:(1,4-10:2,99). (02-01-2017)
//:# - Error corregido en Config(). (23-01-2017)
//:# - Nueva funcion GetValues(). (30-07-2018)
//:#############################################################################

/// \file JRangeFilter.h \brief Declares the class \ref JRangeFilter.

#ifndef _JRangeFilter_
#define _JRangeFilter_

#include "TypesDef.h"
#include "JObject.h"
#include <string>
#include <vector>

//##############################################################################
//# JRangeFilter
//##############################################################################
/// \brief Facilitates filtering values within a list.

class JRangeFilter : protected JObject
{
private:
  unsigned* Ranges;          ///<Stores intervals
  unsigned Count;            ///<Number of intervals stored in \ref Ranges.
  unsigned Size;             ///<Number of intervals allocated in \ref Ranges.

  unsigned ValueMin,ValueMax;
  byte *FastValue;           ///<Array to optimise the values search.

  void ResizeRanges(unsigned size);
  void AddValue(unsigned v);
  void AddRange(unsigned v,unsigned v2);
  void AddRangeStep(unsigned v,unsigned v2,unsigned step);
  void SortRanges();
  void JoinRanges();
  bool CheckNewValue(unsigned v)const;

public:
  JRangeFilter(std::string filter="");
  ~JRangeFilter(){ DestructorActive=true; Reset(); }
  void Reset();
  void Config(std::string filter);
  bool CheckValue(unsigned v)const;
  unsigned GetFirstValue()const;
  unsigned GetNextValue(unsigned v)const;
  bool Empty()const{ return(!Count); }
  std::string ToString()const;
  void GetValues(std::vector<unsigned> &values)const;
};

#endif


