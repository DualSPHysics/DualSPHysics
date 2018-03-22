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

/// \file JSpaceEParms.cpp \brief Implements the class \ref JSpaceEParms.

#include "JSpaceEParms.h"
#include "JXml.h"

//==============================================================================
/// Constructor.
//==============================================================================
JSpaceEParms::JSpaceEParms(){
  ClassName="JSpaceEParms";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSpaceEParms::~JSpaceEParms(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSpaceEParms::Reset(){
  List.clear();
}

//==============================================================================
/// Adds element to the list.
//==============================================================================
void JSpaceEParms::Add(const std::string &key,const std::string &value,const std::string &comment,const std::string &unitscomment){
  JSpaceEParmsItem* item=GetItemPointer(key);
  if(item){ item->value=value; item->comment=comment; item->unitscomment=unitscomment; }
  else{
    JSpaceEParmsItem ite; ite.key=key; ite.value=value; ite.comment=comment; ite.unitscomment=unitscomment;
    List.push_back(ite);
  }
}

//==============================================================================
/// Modifies the value of a pre-existing value.
//==============================================================================
void JSpaceEParms::SetValue(const std::string &key,const std::string &value){
  JSpaceEParmsItem* item=GetItemPointer(key);
  if(!item)RunException("SetValue","The parameter to modify does not exist");
  item->value=value;
}

//==============================================================================
/// Modifies the coment of a pre-existing value.
//==============================================================================
void JSpaceEParms::SetComment(const std::string &key,const std::string &comment){
  JSpaceEParmsItem* item=GetItemPointer(key);
  if(!item)RunException("SetComment","The parameter to modify does not exist");
  item->comment=comment;
}

//==============================================================================
/// Returns a given value associated to the requested key.
//==============================================================================
std::string JSpaceEParms::GetValueNum(const std::string &key,int num){
  std::string value;
  JSpaceEParmsItem* item=GetItemPointer(key);
  if(item){
    std::string ttx=item->value.c_str();
    for(int tc=0;ttx!=""&&tc<=num;tc++){
      int tpos=int(ttx.find(":"));
      std::string ttxopt=(tpos>0? ttx.substr(0,tpos): ttx);
      std::string ttxopt2;
      if(tpos>0)ttxopt2=ttx.substr(tpos+1);
      if(tc==num)value=ttxopt;
      ttx=ttxopt2;
    }   
  }
  return(value);
}

//==============================================================================
/// Returns a given value associated to the requested key.
//==============================================================================
std::string JSpaceEParms::GetValue(const std::string &key){
  std::string value;
  JSpaceEParmsItem* item=GetItemPointer(key);
  if(item)value=item->value.c_str();
  return(value);
}

//==============================================================================
/// Returns a given value associated to the requested key.
//==============================================================================
int JSpaceEParms::GetValueNumInt(const std::string &key,int num,bool optional,int valdef){
  int ret=valdef;
  std::string txval=GetValueNum(key,num);
  if(!txval.empty())ret=atoi(txval.c_str());
  else if(!optional)RunException("GetValueNumInt",std::string("The requested value \'")+key+"\' does not exist.");
  return(ret);
}
//==============================================================================
double JSpaceEParms::GetValueNumDouble(const std::string &key,int num,bool optional,double valdef){
  double ret=valdef;
  std::string txval=GetValueNum(key,num);
  if(!txval.empty())ret=atof(txval.c_str());
  else if(!optional)RunException("GetValueNumDouble",std::string("The requested value \'")+key+"\' does not exist.");
  return(ret);
}
//==============================================================================
std::string JSpaceEParms::GetValueNumStr(const std::string &key,int num,bool optional,std::string valdef){
  std::string ret=valdef;
  std::string txval=GetValueNum(key,num);
  if(!txval.empty())ret=txval;
  else if(!optional)RunException("GetValueNumStr",std::string("The requested value \'")+key+"\' does not exist.");
  return(ret);
}

//==============================================================================
/// Returns the position of a requested element.
//==============================================================================
JSpaceEParms::JSpaceEParmsItem* JSpaceEParms::GetItemPointer(const std::string &key){
  JSpaceEParmsItem *item=NULL;
  for(VecListIt it=List.begin();it<List.end()&&!item;it++)if((*it).key==key)item=&(*it);
  return(item);
}
//==============================================================================
/// Returns a parameter in text format.
//==============================================================================
std::string JSpaceEParms::ToString(unsigned pos)const{
  if(pos>=Count())RunException("ToString","The resquested parameter does not exist.");
  JSpaceEParmsItem ite=List[pos];
  std::string tx=ite.key+"="+ite.value;
  if(ite.comment!=""){
    while(tx.length()<30)tx=tx+" ";
    tx=tx+"#"+ite.comment;
  }
  if(!ite.unitscomment.empty())tx=tx+" #"+ite.unitscomment;
  return(tx);
}
//==============================================================================
/// Returns the requested information of each parameter.
//==============================================================================
JSpaceEParms::JSpaceEParmsItem JSpaceEParms::GetParm(unsigned pos)const{
  if(pos>=Count())RunException("GetParm","The requested parameter does not exist.");
  return(List[pos]);
}
//==============================================================================
/// Loads data of a file in XML format.
//==============================================================================
void JSpaceEParms::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path);
}

//==============================================================================
/// Stores data in a file in XML format.
//==============================================================================
void JSpaceEParms::SaveFileXml(const std::string &file,const std::string &path,bool newfile)const{
  JXml jxml;
  if(!newfile)jxml.LoadFile(file);
  SaveXml(&jxml,path);
  jxml.SaveFile(file);
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSpaceEParms::LoadXml(JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Stores initial conditions of XML object.
//==============================================================================
void JSpaceEParms::SaveXml(JXml *sxml,const std::string &place)const{
  WriteXml(sxml,sxml->GetNode(place,true)->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JSpaceEParms::ReadXml(JXml *sxml,TiXmlElement* lis){
  TiXmlElement* ele=lis->FirstChildElement("parameter"); 
  while(ele){
    Add(sxml->GetAttributeStr(ele,"key"),sxml->GetAttributeStr(ele,"value"),sxml->GetAttributeStr(ele,"comment",true),sxml->GetAttributeStr(ele,"units_comment",true));
    ele=ele->NextSiblingElement("parameter");
  }
}

//==============================================================================
/// Writes list in the XML node.
//==============================================================================
void JSpaceEParms::WriteXml(JXml *sxml,TiXmlElement* lis)const{
  for(unsigned c=0;c<List.size();c++){
    JSpaceEParmsItem ite=List[c];
    TiXmlElement item("parameter");
    JXml::AddAttribute(&item,"key",ite.key);
    JXml::AddAttribute(&item,"value",ite.value);
    if(!ite.comment.empty())JXml::AddAttribute(&item,"comment",ite.comment);
    if(!ite.unitscomment.empty())JXml::AddAttribute(&item,"units_comment",ite.unitscomment);
    lis->InsertEndChild(item);
  }
}


