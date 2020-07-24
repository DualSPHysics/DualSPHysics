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

/// \file JMotionEvent.h \brief Declares the class \ref JMotionEvent.

#ifndef _JMotionEvent_
#define _JMotionEvent_

//#include "TypesDef.h"

class JMotionObj;
class JMotionMov;

//##############################################################################
//# JMotionEvent
//##############################################################################
/// \brief Manages events for motions.

class JMotionEvent
{
public:
  JMotionObj* const Obj;
  JMotionMov* const Mov;
  const double TimeStart;
  const double TimeFinish;

  JMotionEvent(JMotionObj* obj,JMotionMov* mov,double timestart,double timefinish):Obj(obj),Mov(mov),TimeStart(timestart),TimeFinish(timefinish){}
};


#endif

