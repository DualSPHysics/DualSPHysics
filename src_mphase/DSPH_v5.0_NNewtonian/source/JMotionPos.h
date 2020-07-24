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

/// \file JMotionPos.h \brief Declares the class \ref JMotionPos.

#ifndef _JMotionPos_
#define _JMotionPos_

#include "TypesDef.h"
#include <string>
#include <cstdlib>
#include <ctime>
#include "JMatrix4.h"

//##############################################################################
//# JMotionPos
//##############################################################################
/// \brief Manages the position of objects.

class JMotionPos
{
private:
  tdouble3 PosSimple;
  JMatrix4d PosMatrix;
  bool TypeSimple;

public:

  JMotionPos();
  void Reset();
  void Move(const tdouble3 &dis);
  void Rotate(double ang,const tdouble3 &axisp1,const tdouble3 &axisp2);
  void MoveMix(const JMotionPos &modpos);
  void ToMatrix();

  tdouble3 PointMove(const tdouble3 &p) const;
  void PointsMove(tdouble3 &p1,tdouble3 &p2) const{ p1=PointMove(p1); p2=PointMove(p2); }

  bool IsSimple()const{ return(TypeSimple); }
  tdouble3 GetSimple()const{ return(PosSimple); }
  JMatrix4d GetMatrix()const{ return(PosMatrix); }
};


#endif

