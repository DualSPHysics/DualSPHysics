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

/// \file JDsMotionDef.h \brief Includes definitions and specific types for DualSPHysics motion system.

#ifndef _JDsMotionDef_
#define _JDsMotionDef_

#include "TypesDef.h"

///Defines type of movement.
typedef enum{ 
  MOTT_None=0,    ///<No movement.
  MOTT_Linear=1,  ///<Linear movement.
  MOTT_Matrix=2   ///<Matrix movement (for rotations).
}TpMotionType;   

///Structure with the information for moving particles (lineal and matrix movement).
typedef struct{
  word ref;            ///<Idx of moving object.
  word mkbound;        ///<MkBound of moving particles.
  unsigned idbegin;    ///<First id of moving particles.
  unsigned count;      ///<Number of moving particles.
  TpMotionType type;   ///<Type of motion (none, linear, matrix).
  tdouble3 linmov;     ///<Linear displacement to apply to the particles position.
  tdouble3 linvel;     ///<Linear velocity for particles.
  tdouble3 linace;     ///<Linear acceleration for particles (when acceleration movement is computed).
  tmatrix4d matmov;    ///<Matrix transformation to apply to the particles position.
  tmatrix4d matmov2;   ///Matrix transformation to compute acceleration of particles (when acceleration movement is computed).
}StMotionData;


#endif

