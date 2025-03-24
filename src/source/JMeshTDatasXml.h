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
//:# Descripcion:
//:# =============
//:# Clase para lectura/escitura de configuracion en XML.
//:# - Implementacion. (11-09-2020)
//:#############################################################################

/// \file JMeshTDatasXml.h \brief Declares the class \ref JMeshTDatasXml.

#ifndef _JMeshTDatasXml_
#define _JMeshTDatasXml_

#include "JObject.h"
#include "TypesDef.h"
#include "JMeshDataDef.h"
#include <string>
#include <climits>
#include <cfloat>

class JXml;
class TiXmlElement;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshTDatasXml
//##############################################################################
/// \brief Class with some useful features related to XML data.

class JMeshTDatasXml : protected JObject
{
 private:
  static void RunExceptioonStatic(const std::string& srcfile,int srcline
    ,const std::string& method
    ,const std::string& msg,const std::string& file="");

 public:

  static void ReadXmlRho(const JXml* sxml,const TiXmlElement* xele,StMeshRhoCfg& cfg);
  static TiXmlElement* WriteXmlRho(JXml* sxml,TiXmlElement* xele,const StMeshRhoCfg& cfg);

  static void ReadXmlVel(const JXml* sxml,const TiXmlElement* xele,StMeshVelCfg& cfg);
  static TiXmlElement* WriteXmlVel(JXml* sxml,TiXmlElement* xele,const StMeshVelCfg& cfg);

  static void ReadXmlVelExt(const JXml* sxml,const TiXmlElement* xmes,StMeshVelExtCfg& cfg);
};

}

#endif


