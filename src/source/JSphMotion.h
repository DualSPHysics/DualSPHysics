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
//:# - Para los ficheros de datos usa ruta absoluta si el nombre contiene alguna
//:#   barra de directorio. (11-09-2013)
//:#############################################################################

/// \file JSphMotion.h \brief Declares the class \ref JSphMotion.

#ifndef _JSphMotion_
#define _JSphMotion_

#include "TypesDef.h"
#include <string>

class JMotion;
class JXml;


//##############################################################################
//# JSphMotion
//##############################################################################
/// \brief Provides the displacement of moving objects during a time interval.  

class JSphMotion
{
public:

  ///Controls the output of information on the screen and/or log.
  typedef enum{ 
    MOMT_Simple=0,  ///<Simple mode for only forward.
    MOMT_Ace2dt=1,  ///<Calculates acceleration using one dt in the future (always from the beginning).
  }TpMotionMode;   

private:
  JMotion *Mot;

public:

  //==============================================================================
  /// Constructor.
  //==============================================================================
  JSphMotion();

  //==============================================================================
  /// Destructor.
  //==============================================================================
  ~JSphMotion();

  //==============================================================================
  /// Initialisation of variables.
  //==============================================================================
  void Reset();

  //==============================================================================
  /// Initialisation of configuration and returns number of moving objects.
  //==============================================================================
  unsigned Init(JXml *jxml,const std::string &path,const std::string &dirdata);

  //==============================================================================
  /// Returns number of moving objects.
  //==============================================================================
  unsigned GetNumObjects()const;

  //==============================================================================
  /// Processes next time interval and returns true if there are active motions.
  //==============================================================================
  bool ProcesTime(TpMotionMode mode,double timestep,double dt);

  //==============================================================================
  /// Returns data of one moving object. Returns true when the motion is active.
  //==============================================================================
  bool ProcesTimeGetData(unsigned ref,bool &typesimple,tdouble3 &simplemov,tdouble3 &simplevel,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2)const;

  ////==============================================================================
  ///// Processes next time interval and returns true if there are active motions.
  ////==============================================================================
  //bool ProcesTime(double timestep,double dt);

  ////==============================================================================
  ///// Returns the number of performed movements.
  ////==============================================================================
  //unsigned GetMovCount()const;

  ////==============================================================================
  ///// Returns data of the motion of an object.
  ////==============================================================================
  //bool GetMov(unsigned mov,unsigned &ref,tfloat3 &mvsimple,tmatrix4f &mvmatrix)const;
  //bool GetMov(unsigned mov,unsigned &ref,tdouble3 &mvsimple,tmatrix4d &mvmatrix)const;

  ////==============================================================================
  ///// Returns the number of finished movements.
  ////==============================================================================
  //unsigned GetStopCount()const;

  ////==============================================================================
  ///// Returns the reference of the stopped object.
  ////==============================================================================
  //unsigned GetStopRef(unsigned mov)const;
};

#endif


