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

/// \file JSphVResDef.h \brief Defines basic/general definitions and functions for VRes code.

#ifndef _JSphVResDef_
#define _JSphVResDef_

//#include "TypesDef.h"


///Defines the treatment of fluid particles entering a inlet/outlet zone.
typedef enum{ 
   VrOrder_0th=0     ///<Free input mode (no changes).
  ,VrOrder_1st=1  ///<Convert fluid to inlet/outlet.
  ,VrOrder_2nd=2   ///<Remove fluid.
}TpVresOrder;

typedef enum{
   VrMethod_Liu=0
  ,VrMethod_Mls=1
}TpVresMethod;



#endif


