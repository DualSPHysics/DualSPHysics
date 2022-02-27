//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2021 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Implementacion de funciones basicas para uso de dcellcode. (09-09-2021)
//:#############################################################################

/// \file JDsDcell.h \brief Declares the class \ref JDsDcell.

#ifndef _JDsDcell_
#define _JDsDcell_

#include "TypesDef.h"
#include <string>

//##############################################################################
//# JDsDcell
//##############################################################################
/// \brief Defines functions to use dcell code instead of coordinates in 3D.

class JDsDcell
{
public:
  static unsigned    CalcBitsValue(unsigned v,unsigned minbits=1);
  static tuint3      CalcCellDistribution(const tuint3 &ncells,const unsigned maxbits=32);
  static unsigned    CalcCellCode(const tuint3 &ncells);
  static bool        InvalidCellCode(unsigned dcc,const tuint3 &ncells);
  static std::string DcellCodeStr(unsigned dcc);
  static void        PrintDcellCodeInfo(unsigned dcc);
};

#endif


