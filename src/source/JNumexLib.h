//HEAD_DSPH
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
//:# - Interfaz simplificado de JNumex. (17-03-2020)
//:# - CreateVar() permite remplazar variables y en ese caso devuelve true. (20-03-2020)
//:#############################################################################

/// \file JNumexLib.h \brief Declares the class \ref JNumexLib.

#ifndef _JNumexLib_
#define _JNumexLib_

#include "TypesDef.h"
#include <string>
#include <vector>
#include "JNumexLibDef.h"   //Defines DISABLE_NUMEXLIB to compile without Numex library.

class JNumex;

//##############################################################################
//# JNumexLib
//##############################################################################
/// \brief Implements library to evalute user-defined expressions.

#ifdef DISABLE_NUMEXLIB
#include "JNumexLibUndef.h"
#else
class JNumexLib
{
private:
  JNumex* Nux; 
  bool DeleteNux;
  
public:
  /// Constructor.
  JNumexLib(JNumex* nux=NULL);
  /// Destructor.
  ~JNumexLib();
  /// Initialization of variables.
  void ResetVars();
  /// Returns true when this feature is available.
  static bool Available(){ return(true); }

  /// Returns Nux pointer.
  JNumex* GetNuxPtr()const{ return(Nux); }

  /// Creates a new double variable.
  bool CreateVar(const std::string &name,bool cte,bool replace,double value,std::string errtext="");
  /// Creates new double variables (x,y,z) from tdouble3 value.
  bool CreateVar(const std::string &name,bool cte,bool replace,const tdouble3 &value,std::string errtext="");
  /// Creates new double variables (x,y,z) from tfloat3 value.
  bool CreateVar(const std::string &name,bool cte,bool replace,const tfloat3 &value,std::string errtext=""){
    return(CreateVar(name,cte,replace,ToTDouble3(value),errtext));
  }
  /// Creates a new boolean variable.
  bool CreateVar(const std::string &name,bool cte,bool replace,bool value,std::string errtext=""){
    return(CreateVar(name,cte,replace,(value? 1.: 0.),errtext));
  };
  /// Creates a new string variable.
  bool CreateVar(const std::string &name,bool cte,bool replace,const std::string &value,std::string errtext="");

  /// Returns number of variables.
  unsigned CountVars()const;
  /// Returns name of requested variable according to index.
  std::string GetVarName(unsigned idx)const;
  /// Returns true when variable is numeric. Throws exception when it is not found.
  bool VarIsNum(unsigned idx)const;
  /// Returns value of requested variable according to index.
  double GetVarNum(unsigned idx)const;
  /// Returns value of requested variable according to index.
  std::string GetVarStr(unsigned idx)const;
  /// Returns list of selected variables to export or all variables.
  unsigned GetExportVars(std::vector<unsigned> &vars)const;
  /// Returns string with list of variables and its values in one line.
  std::string ListVarsToStr(unsigned firstvar=0)const;

  /// Evaluates numeric expression and returns double result.
  double ComputeExpr(std::string expr,const std::string &errtext="");
  /// Evaluate expression to return a boolean result.
  bool ComputeExprBool(std::string expr,const std::string &errtext="");
  /// Evaluate text expression.
  std::string ComputeExprStr(std::string expr,const std::string &errtext="");

};
#endif

#endif

