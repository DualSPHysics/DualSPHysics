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

/// \file JCaseUserVars.cpp \brief Implements the class \ref JCaseUserVars.

#include "JCaseUserVars.h"
#include "Functions.h"
#include "JNumexLib.h"
#include "JXml.h"
#include <climits>
#include <algorithm>  

using namespace std;

//##############################################################################
//# JCaseUserVars
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCaseUserVars::JCaseUserVars(){
  ClassName="JCaseUserVars";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JCaseUserVars::~JCaseUserVars(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JCaseUserVars::Reset(){
  Vars.clear();
}

//==============================================================================
/// Sort funtion.
//==============================================================================
bool JCaseUserVarsSort(const JCaseUserVars::StVar &a,const JCaseUserVars::StVar &b){ 
  return (a.isnum!=b.isnum? !a.isnum: a.name<b.name);
}

//==============================================================================
/// Configures MkBoundFirst and MkFluidFirst.
//==============================================================================
void JCaseUserVars::LoadExportVars(const JNumexLib *nuxlib){
  Reset();
  std::vector<unsigned> vars;
  const unsigned nv=nuxlib->GetExportVars(vars);
  for(unsigned cv=0;cv<nv;cv++){
    const unsigned idx=vars[cv];
    StVar v;
    v.name=nuxlib->GetVarName(idx);
    v.isnum=nuxlib->VarIsNum(idx);
    v.valuenum=(v.isnum? nuxlib->GetVarNum(idx): 0);
    if(!v.isnum)v.valuestr=nuxlib->GetVarStr(idx);
    Vars.push_back(v);
  }
  //-Sort variables.
  std::sort(Vars.begin(),Vars.end(),JCaseUserVarsSort);
}

//==============================================================================
/// Reads particles information in XML format.
//==============================================================================
void JCaseUserVars::ReadXml(const JXml *sxml,TiXmlElement* lis){
  sxml->CheckElementNames(lis,false,"varstr varnum");
  //-Loads user variables.
  TiXmlElement* ele=lis->FirstChildElement();
  while(ele){
    string ename=ele->Value();
    if(ename.length() && ename[0]!='_'){
      StVar v;
      if(ename=="varstr"){
        v.name=sxml->GetAttributeStr(ele,"name");
        v.isnum=false;
        v.valuenum=0;
        v.valuestr=sxml->GetAttributeStr(ele,"value");
      }
      else if(ename=="varnum"){
        v.name=sxml->GetAttributeStr(ele,"name");
        v.isnum=true;
        v.valuenum=sxml->GetAttributeDouble(ele,"value");
      }
      else sxml->ErrReadElement(ele,ename,false);
      Vars.push_back(v);
    }
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Writes information in XML format.
//==============================================================================
void JCaseUserVars::WriteXml(JXml *sxml,TiXmlElement* lis)const{
  lis->Clear();
  const unsigned nv=unsigned(Vars.size());
  for(unsigned cv=0;cv<nv;cv++){
    const StVar &v=Vars[cv];
    TiXmlElement item(v.isnum? "varnum": "varstr");
    JXml::AddAttribute(&item,"name",v.name);
    if(v.isnum)JXml::AddAttribute(&item,"value",v.valuenum);
    else       JXml::AddAttribute(&item,"value",v.valuestr);
    lis->InsertEndChild(item);
  }
}

//==============================================================================
/// Loads data in XML format from a file.
//==============================================================================
void JCaseUserVars::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path,false);
}

//==============================================================================
/// Loads information from the object XML.
//==============================================================================
void JCaseUserVars::LoadXml(const JXml *sxml,const std::string &place,bool optional){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node && !optional)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(node)ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Stores information in the object XML.
//==============================================================================
void JCaseUserVars::SaveXml(JXml *sxml,const std::string &place)const{
  if(unsigned(Vars.size())>0)WriteXml(sxml,sxml->GetNode(place,true)->ToElement());
}

//==============================================================================
/// Returns data of requested variable according to index.
//==============================================================================
JCaseUserVars::StVar JCaseUserVars::GetVar(unsigned idx)const{
  if(idx>=CountVars())Run_Exceptioon(fun::PrintStr("Index %d of variable is invalid.",idx));
  return(Vars[idx]);
}


