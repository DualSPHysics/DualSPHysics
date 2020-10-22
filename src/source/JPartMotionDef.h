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

/// \file JPartMotionDef.h \brief Defines basic/general definitions and functions for particles motion code.

#ifndef _JPartMotionDef_
#define _JPartMotionDef_

#include "TypesDef.h"

/// Structure with data to compute motion from particles.
typedef struct {
  word mkbound;
  unsigned begin;
  unsigned np;
  //-Reference data to compute motion.
  unsigned nid;   //-Number of selected Idp.
  unsigned id[3]; //-Idp of selected particles.
  tdouble3 ps[3]; //-Initial position of selected particles.
  double dis[3];  //-Maximum distance to first position for periodic boundaries.
}StMkMotionData;

/// Structure with the information of floating body motions.
typedef struct{
  double surge; //-Motion in X.
  double sway;  //-Motion in Y.
  double heave; //-Motion in Z.
  double roll;  //-Rotation in X.
  double pitch; //-Rotation in Y.
  double yaw;   //-Rotation in Z.
}StBodyMotions;


#endif

