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
//:# - Gestiona lista de variables definidas por el usarios en el XML. (10-03-2020)
//:# - Cambio de nombre de J.SpaceUserVars a J.CaseUserVars. (28-06-2020)
//:#############################################################################

/// \file JCaseUserVars.h \brief Declares the class \ref JCaseUserVars.

#ifndef _JCaseUserVars_
#define _JCaseUserVars_

#include "JObject.h"
#include "TypesDef.h"
#include "JNumexLibDef.h"   //Defines DISABLE_NUMEXLIB to compile without Numex library.
#include <vector>

class JXml;
class TiXmlElement;
class JNumexLib;


//##############################################################################
//# JCaseUserVars
//##############################################################################
/// \brief Manages the user-defined variables in output XML.

class JCaseUserVars : protected JObject
{
public:
  typedef struct{
    std::string name;
    bool isnum;
    double valuenum;
    std::string valuestr;
  }StVar;
private:
  std::vector<StVar> Vars;

  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void WriteXml(JXml *sxml,TiXmlElement* lis)const;

public:
  JCaseUserVars();
  ~JCaseUserVars();
  void Reset();

  void LoadExportVars(const JNumexLib *nuxlib);
  void SaveXml(JXml *sxml,const std::string &place)const;

  void LoadFileXml(const std::string &file,const std::string &path);
  void LoadXml(const JXml *sxml,const std::string &place,bool optional);

  unsigned CountVars()const{ return(unsigned(Vars.size())); };
  StVar GetVar(unsigned idx)const;

};

#endif

