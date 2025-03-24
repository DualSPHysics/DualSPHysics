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

/// \file JMeshTDatasXml.cpp \brief Implements the classes JMeshTDatasXml.

#include "JMeshTDatasXml.h"
#include "JException.h"
#include "JXml.h"
#include "Functions.h"

using namespace std;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshTDatasXml
//##############################################################################
//==============================================================================
/// Throws exception related to a file from a static method.
//==============================================================================
void JMeshTDatasXml::RunExceptioonStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,const std::string& msg,const std::string& file)
{
  throw JException(srcfile,srcline,"JMeshTDatasXml",method,msg,file);
}

//==============================================================================
/// Reads density mesh-data configuration from XML.
//==============================================================================
void JMeshTDatasXml::ReadXmlRho(const JXml* sxml,const TiXmlElement* xele
  ,StMeshRhoCfg& cfg)
{
  cfg.Clear();
  cfg.file=sxml->ReadElementStr(xele,"meshdata","file");
  cfg.setpos.x=sxml->ReadElementDouble(xele,"meshdata","setx",true);
  cfg.setpos.y=sxml->ReadElementDouble(xele,"meshdata","sety",true);
  cfg.setpos.z=sxml->ReadElementDouble(xele,"meshdata","setz",true);
  cfg.initialtime=sxml->ReadElementDouble(xele,"meshdata","initialtime",true);
}
//==============================================================================
/// Writes density mesh-data configuration to XML.
//==============================================================================
TiXmlElement* JMeshTDatasXml::WriteXmlRho(JXml* sxml,TiXmlElement* ele
  ,const StMeshRhoCfg& cfg)
{
  TiXmlElement* ele2=sxml->AddElement(ele,"meshdata");
  sxml->AddAttribute(ele2,"file",cfg.file);
  if(cfg.setpos.x)sxml->AddAttribute(ele2,"setx",cfg.setpos.x);
  if(cfg.setpos.y)sxml->AddAttribute(ele2,"sety",cfg.setpos.y);
  if(cfg.setpos.z)sxml->AddAttribute(ele2,"setz",cfg.setpos.z);
  if(cfg.initialtime)sxml->AddAttribute(ele2,"initialtime",cfg.initialtime);
  sxml->AddAttribute(ele2,"comment","Density data from mesh-data file (CSV or binary).");
  sxml->AddAttribute(ele2,"units_comment","kg/m^3");
  return(ele);
}

//==============================================================================
/// Reads velocity mesh-data configuration from XML.
//==============================================================================
void JMeshTDatasXml::ReadXmlVel(const JXml* sxml,const TiXmlElement* xele
  ,StMeshVelCfg& cfg)
{
  cfg.Clear();
  cfg.veldir=sxml->ReadElementFloat3(xele,"direction",true);
  cfg.file=sxml->ReadElementStr(xele,"meshdata","file");
  cfg.setpos.x=sxml->ReadElementDouble(xele,"meshdata","setx",true);
  cfg.setpos.y=sxml->ReadElementDouble(xele,"meshdata","sety",true);
  cfg.setpos.z=sxml->ReadElementDouble(xele,"meshdata","setz",true);
  cfg.initialtime=sxml->ReadElementDouble(xele,"meshdata","initialtime",true);
  cfg.velmagnitude=sxml->ReadElementBool(xele,"meshdata","magnitude");
  cfg.velreverse=sxml->ReadElementBool(xele,"meshdata","reverse",true,false);
}
//==============================================================================
/// Writes velocity mesh-data configuration to XML.
//==============================================================================
TiXmlElement* JMeshTDatasXml::WriteXmlVel(JXml* sxml,TiXmlElement* ele
  ,const StMeshVelCfg& cfg)
{
  if(cfg.veldir!=TFloat3(0)){
    TiXmlElement* ele2=sxml->AddElementFloat3(ele,"direction",cfg.veldir);
    sxml->AddAttribute(ele2,"comment","Direction to apply velocity magnitude or used to calculate starting from 3-component velocity data.");
  }
  TiXmlElement* ele2=sxml->AddElement(ele,"meshdata");
  sxml->AddAttribute(ele2,"file",cfg.file);
  sxml->AddAttribute(ele2,"magnitude",cfg.velmagnitude);
  if(cfg.velreverse)sxml->AddAttribute(ele2,"reverse",cfg.velreverse);
  if(cfg.setpos.x)sxml->AddAttribute(ele2,"setx",cfg.setpos.x);
  if(cfg.setpos.y)sxml->AddAttribute(ele2,"sety",cfg.setpos.y);
  if(cfg.setpos.z)sxml->AddAttribute(ele2,"setz",cfg.setpos.z);
  if(cfg.initialtime)sxml->AddAttribute(ele2,"initialtime",cfg.initialtime);
  sxml->AddAttribute(ele2,"comment","Velocity data (magnitude or 3 components) mesh-data file (CSV or binary).");
  sxml->AddAttribute(ele2,"units_comment","m/s");
  return(ele);
}


//==============================================================================
/// Reads velocity mesh-data configuration from XML.
//==============================================================================
void JMeshTDatasXml::ReadXmlVelExt(const JXml* sxml,const TiXmlElement* xmes
  ,StMeshVelExtCfg& cfg)
{
  cfg.Clear();
  cfg.file=sxml->GetAttributeStr(xmes,"file");
  //-Old input model.
  cfg.setpos.x=sxml->GetAttributeDouble(xmes,"setx",true);
  cfg.setpos.y=sxml->GetAttributeDouble(xmes,"sety",true);
  cfg.setpos.z=sxml->GetAttributeDouble(xmes,"setz",true);
  cfg.initialtime=sxml->GetAttributeDouble(xmes,"initialtime",true);
  cfg.velmagnitude=sxml->GetAttributeBool(xmes,"magnitude",true,true);
  cfg.velreverse=sxml->GetAttributeBool(xmes,"reverse",true,false);
  cfg.velmul.x=sxml->GetAttributeDouble(xmes,"setmulx",true,1);
  cfg.velmul.y=sxml->GetAttributeDouble(xmes,"setmuly",true,1);
  cfg.velmul.z=sxml->GetAttributeDouble(xmes,"setmulz",true,1);
  cfg.veladd.x=sxml->GetAttributeDouble(xmes,"setaddx",true,0);
  cfg.veladd.y=sxml->GetAttributeDouble(xmes,"setaddy",true,0);
  cfg.veladd.z=sxml->GetAttributeDouble(xmes,"setaddz",true,0);
  //-New input model.
  sxml->CheckElementNames(xmes,true,"magnitude reverse initialtime timeloop setpos setvelmul setveladd");
  cfg.velmagnitude=sxml->ReadElementBool(xmes,"magnitude","v",true,cfg.velmagnitude);
  cfg.velreverse=sxml->ReadElementBool(xmes,"reverse","v",true,cfg.velreverse);
  cfg.initialtime=sxml->ReadElementDouble(xmes,"initialtime","v",true,cfg.initialtime);
  cfg.looptbeg=sxml->ReadElementDouble(xmes,"timeloop","tbegin",true,DBL_MAX);
  cfg.looptmax=sxml->ReadElementDouble(xmes,"timeloop","tend",true,DBL_MAX);
  if(cfg.looptbeg!=DBL_MAX || cfg.looptmax!=DBL_MAX){
    if(cfg.looptbeg==DBL_MAX || cfg.looptmax==DBL_MAX)sxml->ErrReadElement(xmes,"timeloop",false,"The value tbegin or tend was not defined.");
    if(cfg.looptbeg>=cfg.looptmax)sxml->ErrReadElement(xmes,"timeloop",false,"The value tbegin is greater than or equal to tend.");
  }
  cfg.setpos.x=sxml->ReadElementDouble(xmes,"setpos","x",true,cfg.setpos.x);
  cfg.setpos.y=sxml->ReadElementDouble(xmes,"setpos","y",true,cfg.setpos.y);
  cfg.setpos.z=sxml->ReadElementDouble(xmes,"setpos","z",true,cfg.setpos.z);
  cfg.velmul.x=sxml->ReadElementDouble(xmes,"setvelmul","x",true,cfg.velmul.x);
  cfg.velmul.y=sxml->ReadElementDouble(xmes,"setvelmul","y",true,cfg.velmul.y);
  cfg.velmul.z=sxml->ReadElementDouble(xmes,"setvelmul","z",true,cfg.velmul.z);
  cfg.veladd.x=sxml->ReadElementDouble(xmes,"setveladd","x",true,cfg.veladd.x);
  cfg.veladd.y=sxml->ReadElementDouble(xmes,"setveladd","y",true,cfg.veladd.y);
  cfg.veladd.z=sxml->ReadElementDouble(xmes,"setveladd","z",true,cfg.veladd.z);
}



}
