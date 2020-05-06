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

/// \file JSpaceProperties.cpp \brief Implements the class \ref JSpaceProperties.

#include "JSpaceProperties.h"
#include "Functions.h"
#include "JXml.h"
#include "JRangeFilter.h"
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSpacePropValue
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSpacePropValue::JSpacePropValue(std::string name):Name(name){
  ClassName="JSpacePropValue";
  Clear();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSpacePropValue::~JSpacePropValue(){
  DestructorActive=true;
  Clear();
}

//==============================================================================
/// Removes values.
//==============================================================================
void JSpacePropValue::Clear(){
  Simple=true;
  Value="";
  SubNames.clear();
  SubValues.clear();
}

//==============================================================================
/// Assigns value without name.
//==============================================================================
void JSpacePropValue::SetValue(std::string v){
  if(!Simple)Clear();
  Value=v;
}

//==============================================================================
/// Assigns value with name.
//==============================================================================
void JSpacePropValue::AddSubValue(std::string subname,std::string v){
  subname=fun::StrLower(subname);
  if(ExistsSubName(subname))Run_Exceptioon("Name-Value already exists.");
  if(Simple){ Value=""; Simple=false; }
  SubNames.push_back(subname);
  SubValues.push_back(v);
}

//==============================================================================
/// Returns position of the given name-value (-1 if it does not exist).
//==============================================================================
int JSpacePropValue::GetIndexSubName(std::string subname)const{
  subname=fun::StrLower(subname);
  unsigned c=0;
  for(;c<SubNames.size() && subname!=SubNames[c];c++);
  return(c<SubNames.size()? int(c): -1);
}

//==============================================================================
/// Retunrs value inline.
//==============================================================================
std::string JSpacePropValue::GetValue()const{
  if(!Simple)Run_Exceptioon("Value has one or more several subvalues.");
  return(Value);
}

//==============================================================================
/// Retunrs a subvalue starting from the name.
//==============================================================================
std::string JSpacePropValue::GetSubValue(std::string subname)const{
  subname=fun::StrLower(subname);
  const int idx=GetIndexSubName(subname);
  if(idx<0)Run_Exceptioon("SubName not found.");
  return(SubValues[idx]);
}

//==============================================================================
/// Returns name of the subvalue.
//==============================================================================
std::string JSpacePropValue::GetSubValueName(unsigned idx)const{
  if(idx>=GetSubValuesCount())Run_Exceptioon("Index of subvalue is invalid.");
  return(SubNames[idx]);
}

//==============================================================================
/// Returns value of the subvalue.
//==============================================================================
std::string JSpacePropValue::GetSubValue(unsigned idx)const{
  if(idx>=GetSubValuesCount())Run_Exceptioon("Index of subvalue is invalid.");
  return(SubValues[idx]);
}

//==============================================================================
/// Returns string with info of object.
//==============================================================================
std::string JSpacePropValue::ToStr()const{
  std::string tx=Name+"=[";
  if(Simple)tx=tx+Value;
  else{
    for(unsigned c=0;c<SubNames.size();c++){
      if(c)tx=tx+";";
      tx=tx+SubNames[c]+"="+SubValues[c];
    }  
  }
  return(tx+"]");
}


//##############################################################################
//# JSpacePropProperty
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSpacePropProperty::JSpacePropProperty(std::string name){
  ClassName="JSpacePropProperty";
  Clear();
  //-Stores Name without spaces.
  Name=fun::StrWithoutChar(name,' ');
}

//==============================================================================
/// Destructor.
//==============================================================================
JSpacePropProperty::~JSpacePropProperty(){
  DestructorActive=true;
  Clear();
}

//==============================================================================
/// Removes values.
//==============================================================================
void JSpacePropProperty::Clear(){
  for(unsigned c=0;c<Values.size();c++)delete Values[c];
  Values.clear();
}

//==============================================================================
/// Copies data from JSpaceProperties object.
//==============================================================================
void JSpacePropProperty::CopyFrom(const JSpacePropProperty* pro){
  Clear();
  //-Copies values.
  for(unsigned c=0;c<unsigned(pro->Values.size());c++){
    const JSpacePropValue* v=pro->Values[c];
    if(v->GetSimple())AddValue(v->GetName(),v->GetValue());
    else{
      const unsigned nv=v->GetSubValuesCount();
      for(unsigned cv=0;cv<nv;cv++)AddSubValue(v->GetName(),v->GetSubValueName(cv),v->GetSubValue(cv));
    }
  }
}

//==============================================================================
/// Returns position of the given value (-1 if it does not exist).
//==============================================================================
int JSpacePropProperty::GetIndexValue(std::string name)const{
  name=fun::StrLower(name);
  unsigned c=0;
  for(;c<Values.size() && name!=fun::StrLower(Values[c]->GetName());c++);
  return(c<Values.size()? int(c): -1);
}

//==============================================================================
/// Adds value.
//==============================================================================
void JSpacePropProperty::AddValue(std::string name,std::string v){
  if(fun::StrLower(name)=="name")Run_Exceptioon("The name of value cannot be \'name\'.");
  if(ExistsNameValue(name))Run_Exceptioon("Value already exists.");
  JSpacePropValue *va=new JSpacePropValue(name);
  Values.push_back(va);
  va->SetValue(v);
}

//==============================================================================
/// Adds value with name.
//==============================================================================
void JSpacePropProperty::AddSubValue(std::string name,std::string subname,std::string v){
  int idx=GetIndexValue(name);
  //-Creates value.
  if(idx<0){
    JSpacePropValue *va=new JSpacePropValue(name);
    Values.push_back(va);
    idx=int(Values.size())-1;
  }
  //-Adds value with name.
  Values[idx]->AddSubValue(subname,v);
}

//==============================================================================
/// Returns name of given value.
//==============================================================================
std::string JSpacePropProperty::GetValueName(unsigned idx)const{
  if(idx>=GetValuesCount())Run_Exceptioon("Index of value is invalid.");
  return(Values[idx]->GetName());
}

//==============================================================================
/// Returns type of given value.
//==============================================================================
bool JSpacePropProperty::GetValueSimple(unsigned idx)const{
  if(idx>=GetValuesCount())Run_Exceptioon("Index of value is invalid.");
  return(Values[idx]->GetSimple());
}

//==============================================================================
/// Returns value of given value.
//==============================================================================
std::string JSpacePropProperty::GetValue(unsigned idx)const{
  if(idx>=GetValuesCount())Run_Exceptioon("Index of value is invalid.");
  return(Values[idx]->GetValue());
}

//==============================================================================
/// Returns pointer to value.
//==============================================================================
const JSpacePropValue* JSpacePropProperty::GetValuePtr(unsigned idx)const{
  if(idx>=GetValuesCount())Run_Exceptioon("Index of value is invalid.");
  return(Values[idx]);
}

//==============================================================================
/// Returns number of subvalues of given value.
//==============================================================================
unsigned JSpacePropProperty::GetSubValuesCount(unsigned idx)const{
  if(idx>=GetValuesCount())Run_Exceptioon("Index of value is invalid.");
  return(Values[idx]->GetSubValuesCount());
}

//==============================================================================
/// Returns name of subvalues of given value.
//==============================================================================
std::string JSpacePropProperty::GetSubValueName(unsigned idx,unsigned subidx)const{
  if(idx>=GetValuesCount())Run_Exceptioon("Index of value is invalid.");
  if(subidx>=GetSubValuesCount(idx))Run_Exceptioon("Index of subvalue is invalid.");
  return(Values[idx]->GetSubValueName(subidx));
}

//==============================================================================
/// Returns value of subvalues of given value.
//==============================================================================
std::string JSpacePropProperty::GetSubValue(unsigned idx,unsigned subidx)const{
  if(idx>=GetValuesCount())Run_Exceptioon("Index of value is invalid.");
  if(subidx>=GetSubValuesCount(idx))Run_Exceptioon("Index of subvalue is invalid.");
  return(Values[idx]->GetSubValue(subidx));
}

//==============================================================================
/// Returns string with info of object.
//==============================================================================
std::string JSpacePropProperty::ToStr()const{
  std::string tx=Name+"={";
  for(unsigned c=0;c<Values.size();c++){
    if(c)tx=tx+", ";
    tx=tx+Values[c]->ToStr();
  }
  return(tx+"}");
}


//##############################################################################
//# JSpacePropLinks
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSpacePropLinks::JSpacePropLinks(){
  ClassName="JSpacePropLinks";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSpacePropLinks::~JSpacePropLinks(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSpacePropLinks::Reset(){
  for(unsigned c=0;c<Links.size();c++)delete Links[c];
  Links.clear();
}

//==============================================================================
/// Copies data from JSpacePropLinks object.
//==============================================================================
void JSpacePropLinks::CopyFrom(const JSpacePropLinks* links){
  Reset();
  for(unsigned c=0;c<unsigned(links->Links.size());c++){
    const JSpacePropLink* lk=links->Links[c];
    AddLink(lk->GetType(),lk->GetMks(),lk->GetProps());
  }
}

//==============================================================================
/// Adds new link.
//==============================================================================
void JSpacePropLinks::AddLink(JSpacePropLink::TpLink type,std::string mks,std::string props){
  mks=fun::StrWithoutChar(mks,' ');
  props=fun::StrWithoutChar(props,' ');
  if(!mks.empty()&&!props.empty()){
    JSpacePropLink *link=new JSpacePropLink(type,mks,props);
    Links.push_back(link);
  }
}

//==============================================================================
/// Reads the links section in XML format.
//==============================================================================
void JSpacePropLinks::ReadXml(const JXml *sxml,TiXmlElement* eprops){
  Reset();
  TiXmlElement* elinks=eprops->FirstChildElement("links");
  if(elinks){
    TiXmlElement* ele=elinks->FirstChildElement("link"); 
    while(ele){
      JSpacePropLink::TpLink type=JSpacePropLink::LINK_Mk;
      string mks;
      if(sxml->ExistsAttribute(ele,"mkbound")){
        type=JSpacePropLink::LINK_MkBound;
        mks=sxml->GetAttributeStr(ele,"mkbound");
      }
      else if(sxml->ExistsAttribute(ele,"mkfluid")){
        type=JSpacePropLink::LINK_MkFluid;
        mks=sxml->GetAttributeStr(ele,"mkfluid");
      }
      else mks=sxml->GetAttributeStr(ele,"mk");
      AddLink(type,mks,sxml->GetAttributeStr(ele,"property"));
      ele=ele->NextSiblingElement("link");
    }
  }
}

//==============================================================================
/// Writes the links section in XML format.
//==============================================================================
void JSpacePropLinks::WriteXml(JXml *sxml,TiXmlElement* eprops)const{
  if(Links.size()){
    TiXmlElement item("links");
    TiXmlElement* elinks=eprops->InsertEndChild(item)->ToElement();
    for(unsigned c=0;c<Links.size();c++){
      string typemk;
      if(Links[c]->GetType()==JSpacePropLink::LINK_Mk)typemk="mk";
      else if(Links[c]->GetType()==JSpacePropLink::LINK_MkBound)typemk="mkbound";
      else if(Links[c]->GetType()==JSpacePropLink::LINK_MkFluid)typemk="mkfluid";
      TiXmlElement* elink=JXml::AddElementAttrib(elinks,"link",typemk,Links[c]->GetMks());
      JXml::AddAttribute(elink,"property",Links[c]->GetProps());
    }
  }
}

//==============================================================================
/// Returns list of properties in order and without repeated properties.
//==============================================================================
std::string JSpacePropLinks::GetPropsSort(std::string props){
  //-Loads array with no repeated properties.
  std::vector<string> vprops;
  while(!props.empty()){
    string pro=fun::StrSplit("+",props);
    if(!pro.empty()){
      unsigned c=0;
      for(;c<vprops.size() && pro!=vprops[c];c++);
      if(c==vprops.size())vprops.push_back(pro);
    }
  }
  //-Reorders properties.
  if(vprops.size())for(unsigned c=0;c<vprops.size()-1;c++)for(unsigned c2=c+1;c2<vprops.size();c2++)if(vprops[c]>vprops[c2]){
    string aux=vprops[c];
    vprops[c]=vprops[c2];
    vprops[c2]=aux;
  }
  //-Returns list of properties in order and without repeated properties.
  props="";
  if(vprops.size())for(unsigned c=0;c<vprops.size();c++){
    if(c)props=props+"+";
    props=props+vprops[c];
  }
  return(props);
}

//==============================================================================
/// Returns the properties associated to each MK.
//==============================================================================
void JSpacePropLinks::GetPropsList(std::vector<std::string> &vprops_mk
  ,std::vector<std::string> &vprops_mkb
  ,std::vector<std::string> &vprops_mkf)const
{
  const unsigned nlinks=unsigned(Links.size());
  for(unsigned c=0;c<nlinks;c++){
    const JSpacePropLink* lk=Links[c];
    int tp=(lk->GetType()==JSpacePropLink::LINK_MkBound? 1: (lk->GetType()==JSpacePropLink::LINK_MkFluid? 2: 0));
    vector<string> &vprops=(lk->GetType()==JSpacePropLink::LINK_MkBound? vprops_mkb: (lk->GetType()==JSpacePropLink::LINK_MkFluid? vprops_mkf: vprops_mk));
    JRangeFilter rg(lk->GetMks());
    unsigned v=rg.GetFirstValue();
    while(v!=UINT_MAX){
      const word mk=word(v);
      while(mk>=unsigned(vprops.size()))vprops.push_back("");//-Add new mk entries.
      vprops[mk]=vprops[mk]+lk->GetProps()+"+";
      v=rg.GetNextValue(v);
    }
  }
  for(unsigned cmk=0;cmk<unsigned(vprops_mk .size());cmk++)vprops_mk [cmk]=GetPropsSort(vprops_mk [cmk]);
  for(unsigned cmk=0;cmk<unsigned(vprops_mkb.size());cmk++)vprops_mkb[cmk]=GetPropsSort(vprops_mkb[cmk]);
  for(unsigned cmk=0;cmk<unsigned(vprops_mkf.size());cmk++)vprops_mkf[cmk]=GetPropsSort(vprops_mkf[cmk]);
  /*:if(0){
    unsigned size=unsigned(vprops_mk.size());
    size=max(size,unsigned(vprops_mkb.size()));
    size=max(size,unsigned(vprops_mkf.size()));
    for(unsigned cmk=0;cmk<size;cmk++){
      string mk =(cmk<unsigned(vprops_mk .size())? vprops_mk [cmk]: "");
      string mkb=(cmk<unsigned(vprops_mkb.size())? vprops_mkb[cmk]: "");
      string mkf=(cmk<unsigned(vprops_mkf.size())? vprops_mkf[cmk]: "");
      printf("GetPropsList --> [%u]=[%s]-[%s]-[%s] \n",cmk,mk.c_str(),mkb.c_str(),mkf.c_str());
    }
  }:*/
}

//==============================================================================
/// Returns the properties associated to the given MK.
//==============================================================================
std::string JSpacePropLinks::GetPropsFast(word mk,word mkboundfirst,word mkfluidfirst
  ,std::vector<std::string> &vprops_mk
  ,std::vector<std::string> &vprops_mkb
  ,std::vector<std::string> &vprops_mkf)const
{
  //-Gets properties assigned to MK.
  string props;
  if(unsigned(mk)<unsigned(vprops_mk.size()))props=props+vprops_mk[mk]+"+";
  if(mk>=mkboundfirst){
    const word mk2=mk-mkboundfirst;
    if(unsigned(mk2)<unsigned(vprops_mkb.size()))props=props+vprops_mkb[mk2]+"+";
  }
  if(mk>=mkfluidfirst && mk<mkboundfirst){
    const word mk2=mk-mkfluidfirst;
    if(unsigned(mk2)<unsigned(vprops_mkf.size()))props=props+vprops_mkf[mk2]+"+";
  }
  return(GetPropsSort(props));
}

//==============================================================================
/// Returns the properties associated to the given MK.
//==============================================================================
std::string JSpacePropLinks::GetProps(word mk,word mkboundfirst,word mkfluidfirst)const{
  //-Gets properties assigned to MK.
  string props;
  for(unsigned c=0;c<unsigned(Links.size());c++){
    const JSpacePropLink* lk=Links[c];
    JRangeFilter rg(lk->GetMks());
    if(lk->GetType()==JSpacePropLink::LINK_Mk){
      if(rg.CheckValue(mk))props=props+lk->GetProps()+"+";
    }
    if(lk->GetType()==JSpacePropLink::LINK_MkBound && mk>=mkboundfirst){
      if(rg.CheckValue(mk-mkboundfirst))props=props+lk->GetProps()+"+";
    }
    if(lk->GetType()==JSpacePropLink::LINK_MkFluid && mk>=mkfluidfirst && mk<mkboundfirst){
      if(rg.CheckValue(mk-mkfluidfirst))props=props+lk->GetProps()+"+";
    }
  }
  return(GetPropsSort(props));
}

//==============================================================================
/// Returns the properties associated to the given MK.
//==============================================================================
std::string JSpacePropLinks::GetProps(word mk)const{
  //-Gets properties assigned to MK.
  string props;
  for(unsigned c=0;c<unsigned(Links.size());c++){
    const JSpacePropLink* lk=Links[c];
    JRangeFilter rg(lk->GetMks());
    if(lk->GetType()==JSpacePropLink::LINK_Mk){
      if(rg.CheckValue(mk))props=props+lk->GetProps()+"+";
    }
    else Run_Exceptioon("Type of links (mkbound or mkfluid) is invalid.");
  }
  return(GetPropsSort(props));
}

//==============================================================================
/// Returns list with all properties that will be used.
//==============================================================================
std::string JSpacePropLinks::GetAllProps()const{
  string props;
  for(unsigned c=0;c<unsigned(Links.size());c++)props=props+Links[c]->GetProps()+"+";
  return(GetPropsSort(props));
}


//##############################################################################
//# JSpaceProperties
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSpaceProperties::JSpaceProperties(){
  ClassName="JSpaceProperties";
  Links=new JSpacePropLinks();
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSpaceProperties::~JSpaceProperties(){
  DestructorActive=true;
  Reset();
  delete Links; Links=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSpaceProperties::Reset(){
  for(unsigned c=0;c<Props.size();c++)delete Props[c];
  Props.clear();
  PropsFile.clear();
  Links->Reset();
}

//==============================================================================
/// Returns index of property with given name (-1 if it does not exist).
//==============================================================================
int JSpaceProperties::GetIndexProperty(std::string name)const{
  name=fun::StrLower(name);
  unsigned c=0;
  for(;c<Props.size() && name!=fun::StrLower(Props[c]->GetName());c++);
  return(c<Props.size()? int(c): -1);
}

//==============================================================================
/// Returns given property (throws exception it does not exist).
//==============================================================================
const JSpacePropProperty* JSpaceProperties::GetProperty(std::string name)const{
  int idx=GetIndexProperty(name);
  if(idx<0)Run_Exceptioon(string("Property \'")+name+"\' not found.");
  return(Props[idx]);
}

//==============================================================================
/// Copies data from JSpaceProperties object.
//==============================================================================
void JSpaceProperties::CopyFrom(const JSpaceProperties* props){
  Reset();
  //-Copies links.
  Links->CopyFrom(props->Links);
  //-Copies properties.
  for(unsigned c=0;c<unsigned(props->Props.size());c++){
    const JSpacePropProperty* pro=props->Props[c];
    JSpacePropProperty* pro2=AddProperty(pro->GetName());
    pro2->CopyFrom(pro);
  }
}

//==============================================================================
/// Loads data in XML format from a file.
//==============================================================================
void JSpaceProperties::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path);
}

//==============================================================================
/// Stores data in XML format in a file.
//==============================================================================
void JSpaceProperties::SaveFileXml(const std::string &file,const std::string &path,bool newfile)const{
  JXml jxml;
  if(!newfile)jxml.LoadFile(file);
  SaveXml(&jxml,path);
  jxml.SaveFile(file);
}

//==============================================================================
/// Loads initial conditions from the XML object.
//==============================================================================
void JSpaceProperties::LoadXml(const JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(node)ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Stores initial conditions in the XML object.
//==============================================================================
void JSpaceProperties::SaveXml(JXml *sxml,const std::string &place)const{
  WriteXml(sxml,sxml->GetNode(place,true)->ToElement());
}

//==============================================================================
/// Reads property in XML format.
//==============================================================================
void JSpaceProperties::ReadXmlProperty(const JXml *sxml,TiXmlElement* eprop){
  //-Creates property if it does not exist.
  JSpacePropProperty *pro=AddProperty(sxml->GetAttributeStr(eprop,"name"));
  //-Loads values inline.
  TiXmlAttribute *at=eprop->FirstAttribute();
  while(at){
    string name=at->Name();
    string val=at->Value();
    if(name!="name"&&name[0]!='_'){
      pro->AddValue(name,val);
    }
    at=at->Next();
  }
  //-Loads values in items.
  TiXmlElement* ele=eprop->FirstChildElement(); 
  while(ele){
    string name=ele->Value();
    if(name[0]!='_'){
      TiXmlAttribute *at=ele->FirstAttribute();
      while(at){
        string atname=at->Name();
        string atval=at->Value();
        if(atname[0]!='_'){
          pro->AddSubValue(name,atname,atval);
        }
        at=at->Next();
      }
    }
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Reads external file with properties in XML format.
//==============================================================================
void JSpaceProperties::ReadXmlPropertyFile(const JXml *sxml,TiXmlElement* epropfile){
  std::string file=sxml->GetAttributeStr(epropfile,"file");
  std::string path=sxml->GetAttributeStr(epropfile,"path");
  //-Loads XML file.
  JXml jxml;
  jxml.LoadFile(file);
  //-Gets node with data.
  TiXmlNode* node=jxml.GetNode(path,false);
  if(!node)Run_ExceptioonFile(std::string("Cannot find the element \'")+path+"\'.",file);
  //-Loads properties of node.
  TiXmlElement *eprops=node->ToElement();
  TiXmlElement* ele=eprops->FirstChildElement(); 
  while(ele){
    std::string cmd=ele->Value();
    if(cmd=="property")ReadXmlProperty(&jxml,ele);
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Writes property in XML format.
//==============================================================================
void JSpaceProperties::WriteXmlProperty(JXml *sxml,TiXmlElement* eprops,const JSpacePropProperty* prop)const{
  TiXmlElement* eprop=JXml::AddElementAttrib(eprops,"property","name",prop->GetName());
  //-Writes values inline.
  for(unsigned c=0;c<prop->GetValuesCount();c++){
    if(prop->GetValueSimple(c))JXml::AddAttribute(eprop,prop->GetValueName(c),prop->GetValue(c));
  }
  //-Writes other values.
  for(unsigned c=0;c<prop->GetValuesCount();c++){
    if(!prop->GetValueSimple(c)){
      TiXmlElement* ele=JXml::AddElement(eprop,prop->GetValueName(c));
      for(unsigned c2=0;c2<prop->GetSubValuesCount(c);c2++){
        JXml::AddAttribute(ele,prop->GetSubValueName(c,c2),prop->GetSubValue(c,c2));
      }
    }
  }
}

//==============================================================================
/// Writes external file with properties in XML format.
//==============================================================================
void JSpaceProperties::WriteXmlPropertyFile(JXml *sxml,TiXmlElement* eprops,StPropertyFile propfile)const{
  TiXmlElement* ele=JXml::AddElement(eprops,"propertyfile");
  JXml::AddAttribute(ele,"file",propfile.file);
  JXml::AddAttribute(ele,"path",propfile.path);
}

//==============================================================================
/// Reads the list of initial conditions in XML format.
//==============================================================================
void JSpaceProperties::ReadXml(const JXml *sxml,TiXmlElement* eprops){
  //-Loads links.
  Links->ReadXml(sxml,eprops);
  //-Loads properties.
  TiXmlElement* ele=eprops->FirstChildElement(); 
  while(ele){
    std::string cmd=ele->Value();
    if(cmd.length()&&cmd[0]!='_'){
      if(cmd=="property")ReadXmlProperty(sxml,ele);
      else if(cmd=="propertyfile")ReadXmlPropertyFile(sxml,ele);
      else if(cmd!="links")sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
  //-Checks relations.
  CheckLinks();
}

//==============================================================================
/// Writes the list in XML format.
//==============================================================================
void JSpaceProperties::WriteXml(JXml *sxml,TiXmlElement* eprops)const{
  //-Stores links.
  Links->WriteXml(sxml,eprops);
  //-Stores property files.
  for(unsigned c=0;c<PropsFile.size();c++)WriteXmlPropertyFile(sxml,eprops,PropsFile[c]);
  //-Stores properties.
  for(unsigned c=0;c<Props.size();c++)WriteXmlProperty(sxml,eprops,Props[c]);
}

//==============================================================================
/// Checks that there are all properties that will be used.
//==============================================================================
void JSpaceProperties::CheckLinks(){
  string props=Links->GetAllProps();
  while(!props.empty()){
    string pro=fun::StrSplit("+",props);
    if(!pro.empty() && GetIndexProperty(pro)<0)Run_Exceptioon("Property \'"+pro+"\' not found.");
  }
}

//==============================================================================
/// Fliters information only for MK values of mkselect.
//==============================================================================
//#include "JTimer.h"
void JSpaceProperties::FilterMk(word mkboundfirst,word mkfluidfirst,std::string mkselect){
  vector<string> vprops;
  vector<string> vmks;
  //-Gets properties of each MK of mkselect.
  //JTimer timer; timer.Start();  
  if(1){//-Fast method when number of mks is very high.
    vector<string> vprops_mk;
    vector<string> vprops_mkb;
    vector<string> vprops_mkf;
    Links->GetPropsList(vprops_mk,vprops_mkb,vprops_mkf);
    JRangeFilter rg(mkselect);
    unsigned v=rg.GetFirstValue();
    while(v!=UINT_MAX){
      word mk=word(v);
      string props=Links->GetPropsFast(mk,mkboundfirst,mkfluidfirst,vprops_mk,vprops_mkb,vprops_mkf);
      ////-Checks new method results.
      //string props2=Links->GetProps(mk,mkboundfirst,mkfluidfirst);
      //if(props!=props2){
      //  printf("\nError mk: %d  [%s]!=[%s] \n",mk,props.c_str(),props2.c_str());
      //  Run_Exceptioon("Error in result of GetPropsFast().");
      //}
      if(!props.empty()){
        unsigned c=0;
        for(;c<vprops.size() && props!=vprops[c];c++);
        if(c==vprops.size()){
          vprops.push_back(props);
          vmks.push_back(fun::IntStr(mk));
        }
        else vmks[c]=vmks[c]+","+fun::IntStr(mk);
      }
      v=rg.GetNextValue(v);
    }
    //timer.Stop(); printf("\nTime of fast method: %.3f sec.\n\n",timer.GetElapsedTimeF()/1000.f); fflush(stdout);
  }
  else{//-Old method, it is too slow when number of mks is very high.
    JRangeFilter rg(mkselect);
    unsigned v=rg.GetFirstValue();
    while(v!=UINT_MAX){
      word mk=word(v);
      string props=Links->GetProps(mk,mkboundfirst,mkfluidfirst);
      if(!props.empty()){
        unsigned c=0;
        for(;c<vprops.size() && props!=vprops[c];c++);
        if(c==vprops.size()){
          vprops.push_back(props);
          vmks.push_back(fun::IntStr(mk));
        }
        else vmks[c]=vmks[c]+","+fun::IntStr(mk);
      }
      v=rg.GetNextValue(v);
    }
    //timer.Stop(); printf("\nTime of old method: %.3f sec.\n\n",timer.GetElapsedTimeF()/1000.f); fflush(stdout);
  }

  //-Reorders MK's.
  for(unsigned c=0;c<vmks.size();c++){
    JRangeFilter rg(vmks[c]);
    vmks[c]=rg.ToString();
  }
  //-Generates list of properties that will be used.
  std::vector<string> vprop;
  for(unsigned c=0;c<vprops.size();c++){
    string props=vprops[c];
    while(!props.empty()){
      string pro=fun::StrSplit("+",props);
      if(!pro.empty()){
        unsigned c=0;
        for(;c<vprop.size() && pro!=vprop[c];c++);
        if(c==vprop.size())vprop.push_back(pro);
      }
    }
  }
  //-Removes properties that will not be used.
  for(unsigned p=0;p<Props.size();p++){
    string pro=Props[p]->GetName();
    unsigned c=0;
    for(;c<vprop.size() && pro!=vprop[c];c++);
    if(c==vprop.size()){
      Props.erase(Props.begin()+p);
      p--;
    }
  }
  //-Remakes links with mkselect.
  Links->Reset();
  for(unsigned c=0;c<vmks.size();c++)Links->AddLink(JSpacePropLink::LINK_Mk,vmks[c],vprops[c]);
  //-Checks links.
  CheckLinks();
}

//==============================================================================
/// Returns properties for a given MK.
//==============================================================================
std::string JSpaceProperties::GetPropertyMk(word mk)const{
  return(Links->GetProps(mk));
}

//==============================================================================
/// Returns pointer to value with the given index.
//==============================================================================
const JSpacePropValue* JSpaceProperties::GetValuePtr(std::string props,unsigned idx)const{
  const JSpacePropValue* value=NULL;
  unsigned n=0;
  while(!props.empty() && !value){
    std::string pro=fun::StrSplit("+",props);
    if(!pro.empty()){
      const JSpacePropProperty* pr=GetProperty(pro);
      const unsigned nv=pr->GetValuesCount();
      if(nv){
        if(idx-n<nv)value=pr->GetValuePtr(idx-n);
        n+=nv;
      }
    }
  }
  return(value);
}

//==============================================================================
/// Returns value with the given name.
//==============================================================================
const JSpacePropValue* JSpaceProperties::GetValuePtr(std::string props,std::string name)const{
  const JSpacePropValue* value=NULL;
  while(!props.empty() && !value){
    std::string pro=fun::StrSplit("+",props);
    if(!pro.empty()){
      const JSpacePropProperty* pr=GetProperty(pro);
      const unsigned nv=pr->GetValuesCount();
      for(unsigned cv=0;cv<nv && !value;cv++){
        if(name==pr->GetValueName(cv))value=pr->GetValuePtr(cv);
      }
    }
  }
  return(value);
}

//==============================================================================
/// Returns value with the given index.
//==============================================================================
const JSpacePropValue* JSpaceProperties::GetValue(std::string props,unsigned idx)const{
  const JSpacePropValue* value=GetValuePtr(props,idx);
  if(!value)Run_Exceptioon("Index of value is invalid.");
  return(value);
}

//==============================================================================
/// Returns value with the name indicated.
//==============================================================================
const JSpacePropValue* JSpaceProperties::GetValue(std::string props,std::string name)const{
  const JSpacePropValue* value=GetValuePtr(props,name);
  if(!value)Run_Exceptioon(string("Value \'")+name+"\' not found.");
  return(value);
}

//==============================================================================
/// Returns if the given value exists.
//==============================================================================
bool JSpaceProperties::ExistsValue(std::string props,std::string name)const{
  const JSpacePropValue* value=GetValuePtr(props,name);
  return(value!=NULL);
}

//==============================================================================
/// Returns number of values.
//==============================================================================
unsigned JSpaceProperties::GetValuesCount(std::string props)const{
  unsigned n=0;
  while(!props.empty()){
    std::string pro=fun::StrSplit("+",props);
    if(!pro.empty())n+=GetProperty(pro)->GetValuesCount();
  }
  return(n);
}

//==============================================================================
/// Returns name of given value.
//==============================================================================
std::string JSpaceProperties::GetValueName(std::string props,unsigned idx)const{
  return(GetValue(props,idx)->GetName());
}

//==============================================================================
/// Returns value as string of given value.
//==============================================================================
std::string JSpaceProperties::GetValueStr(std::string props,unsigned idx)const{
  const JSpacePropValue* value=GetValue(props,idx);
  string val;
  if(value->GetSimple())val=value->GetValue();
  else{
    const unsigned nv=value->GetSubValuesCount();
    for(unsigned c=0;c<nv;c++){
      if(c)val=val+";";
      val=val+value->GetSubValueName(c)+"="+value->GetSubValue(c);
    }  
  }
  return(val);
}

//==============================================================================
/// Returns value as string of given value.
//==============================================================================
std::string JSpaceProperties::GetValueStr(std::string props,std::string name)const{
  const JSpacePropValue* value=GetValue(props,name);
  string val;
  if(value->GetSimple())val=value->GetValue();
  else{
    const unsigned nv=value->GetSubValuesCount();
    for(unsigned c=0;c<nv;c++){
      if(c)val=val+";";
      val=val+value->GetSubValueName(c)+"="+value->GetSubValue(c);
    }  
  }
  return(val);
}

//==============================================================================
/// Returns number of subvalues.
//==============================================================================
unsigned JSpaceProperties::GetSubValuesCount(std::string props,unsigned idx)const{
  return(GetValue(props,idx)->GetSubValuesCount());
}

//==============================================================================
/// Returns name of given subvalue.
//==============================================================================
std::string JSpaceProperties::GetSubValueName(std::string props,unsigned idx,unsigned subidx)const{
  return(GetValue(props,idx)->GetSubValueName(subidx));
}

//==============================================================================
/// Returns value as string of given subvalue.
//==============================================================================
std::string JSpaceProperties::GetSubValueStr(std::string props,unsigned idx,unsigned subidx)const{
  return(GetValue(props,idx)->GetSubValue(subidx));
}

//==============================================================================
/// Rreturns if the given subvalue exists.
//==============================================================================
bool JSpaceProperties::ExistsSubValue(std::string props,std::string name,std::string subname)const{
  const JSpacePropValue* value=GetValuePtr(props,name);
  return(value && value->ExistsSubName(subname));
}

//==============================================================================
/// Returns value as string of given subvalue..
//==============================================================================
std::string JSpaceProperties::GetSubValueStr(std::string props,std::string name,std::string subname)const{
  return(GetValue(props,name)->GetSubValue(subname));
}

//==============================================================================
/// Adds the given link for absolute MK.
//==============================================================================
void JSpaceProperties::AddLink(std::string mks,std::string props){
  Links->AddLink(JSpacePropLink::LINK_Mk,mks,props);
}

//==============================================================================
/// Adds the given link for bound MK.
//==============================================================================
void JSpaceProperties::AddLinkBound(std::string mks,std::string props){
  Links->AddLink(JSpacePropLink::LINK_MkBound,mks,props);
}

//==============================================================================
/// Adds the given link for fluid MK.
//==============================================================================
void JSpaceProperties::AddLinkFluid(std::string mks,std::string props){
  Links->AddLink(JSpacePropLink::LINK_MkFluid,mks,props);
}

//==============================================================================
/// It it does not exist, it will create a property with the given name and returns it.
//==============================================================================
JSpacePropProperty* JSpaceProperties::AddProperty(std::string name){
  name=fun::StrWithoutChar(name,' ');
  //-Checks if a property with the given name already exists.
  int idx=GetIndexProperty(name);
  //-Creates a new property with that name.
  if(idx<0){
    JSpacePropProperty *pro=new JSpacePropProperty(name);
    Props.push_back(pro);
    idx=int(Props.size())-1;
  }
  return(Props[idx]);
}

//==============================================================================
/// Creates reference to an external file.
//==============================================================================
void JSpaceProperties::AddPropertyFile(std::string file,std::string path){
  StPropertyFile data={file,path};
  PropsFile.push_back(data);
}

//==============================================================================
/// Removes list of external files.
//==============================================================================
void JSpaceProperties::ClearPropertyFile(){
  PropsFile.clear();
}


