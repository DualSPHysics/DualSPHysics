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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Uso de DestructorActive para controlar cuando esta el destructor en 
//:#   ejecucion para la generacion de excepciones. (21-03-2018)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

/// \file JObject.h \brief Declares the class \ref JObject.

#ifndef _JObject_
#define _JObject_

#include "RunExceptionDef.h"
#include <string>


//##############################################################################
//# JObject
//##############################################################################
/// \brief Defines objects with methods that throw exceptions.

class JObject
{
protected:
  std::string ClassName;   ///<Name of the class.
  bool DestructorActive;   ///<Destructor of object in execution.

  //static void RunExceptioonStatic(const std::string& srcfile,int srcline
  //  ,const std::string& method
  //  ,const std::string& msg,const std::string& file="");

  void RunExceptioon(const std::string& srcfile,int srcline
    ,const std::string& classname,const std::string& method
    ,const std::string& msg,const std::string& file="")const;

  void RunException(const std::string& method,const std::string& msg
    ,const std::string& file="")const;

public:  
  JObject():ClassName("JObject"),DestructorActive(false){} ///<Constructor of objects.
  virtual std::string ObjectId()const;

};

#endif


