//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2019 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "Functions.h"
#include "JXml.h"
#include <algorithm>

using namespace std;

//##############################################################################
//# JSpaceEParms
//##############################################################################
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
  Posminx=Posminy=Posminz=Posmaxx=Posmaxy=Posmaxz="default";
}

//==============================================================================
/// Adds element to the list.
//==============================================================================
void JSpaceEParms::Add(const std::string &key,const std::string &value,const std::string &comment,const std::string &unitscomment){
  JSpaceEParmsItem* item=GetItemPointer(key);
  string comment2=comment;
  //-Checks deprecated parameters.
  if(1){
    const string k=fun::StrLower(key);
    if(k=="incz" || k.substr(0,11)=="domainfixed" || k.substr(0,15)=="domainparticles" || k=="deltasph" || k=="posdouble"){
      if(int(fun::StrUpper(comment2).find("**DEPRECATED**"))<0)
        comment2=string("**DEPRECATED** ")+comment2;
    }
  }
  if(item){ item->value=value; item->comment=comment2; item->unitscomment=unitscomment; }
  else{
    JSpaceEParmsItem ite; ite.key=key; ite.value=value; ite.comment=comment2; ite.unitscomment=unitscomment;
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
/// Modifies values Posminx, Posminy and Posminz.
//==============================================================================
void JSpaceEParms::SetPosmin(std::string x,std::string y,std::string z){
  const char met[]="SetPosmin";
  const bool isposmin=true;
  JSpaceEParmsPos ps;
  if(CheckPosValue(x,isposmin,ps))RunException(met,fun::PrintStr("The value posmin.x=\"%s\" is invalid.",x.c_str()));
  Posminx=ps.textmod;
  if(CheckPosValue(y,isposmin,ps))RunException(met,fun::PrintStr("The value posmin.y=\"%s\" is invalid.",y.c_str()));
  Posminy=ps.textmod;
  if(CheckPosValue(z,isposmin,ps))RunException(met,fun::PrintStr("The value posmin.z=\"%s\" is invalid.",z.c_str()));
  Posminz=ps.textmod;
}

//==============================================================================
/// Modifies values Posmaxx, Posmaxy and Posmaxz.
//==============================================================================
void JSpaceEParms::SetPosmax(std::string x,std::string y,std::string z){
  const char met[]="SetPosmax";
  const bool isposmin=false;
  JSpaceEParmsPos ps;
  if(CheckPosValue(x,isposmin,ps))RunException(met,fun::PrintStr("The value posmax.x=\"%s\" is invalid.",x.c_str()));
  Posmaxx=ps.textmod;
  if(CheckPosValue(y,isposmin,ps))RunException(met,fun::PrintStr("The value posmax.y=\"%s\" is invalid.",y.c_str()));
  Posmaxy=ps.textmod;
  if(CheckPosValue(z,isposmin,ps))RunException(met,fun::PrintStr("The value posmax.z=\"%s\" is invalid.",z.c_str()));
  Posmaxz=ps.textmod;
}

//==============================================================================
/// Returns values of Posminx, Posminy or Posminz.
//==============================================================================
JSpaceEParms::JSpaceEParmsPos JSpaceEParms::GetPosminValue(char key)const{
  const char met[]="GetPosminValue";
  const bool isposmin=true;
  JSpaceEParmsPos ps;
  switch(key){
    case 'x':  if(CheckPosValue(Posminx,isposmin,ps))RunException(met,fun::PrintStr("The value posmin.x=\"%s\" is invalid.",Posminx.c_str()));  break;
    case 'y':  if(CheckPosValue(Posminy,isposmin,ps))RunException(met,fun::PrintStr("The value posmin.y=\"%s\" is invalid.",Posminy.c_str()));  break;
    case 'z':  if(CheckPosValue(Posminz,isposmin,ps))RunException(met,fun::PrintStr("The value posmin.z=\"%s\" is invalid.",Posminz.c_str()));  break;
    default:   RunException(met,fun::PrintStr("Key value \'%c\' is invalid.",key));
  }
  return(ps);
}

//==============================================================================
/// Returns values of Posmaxx, Posmaxy or Posmaxz.
//==============================================================================
JSpaceEParms::JSpaceEParmsPos JSpaceEParms::GetPosmaxValue(char key)const{
  const char met[]="GetPosmaxValue";
  const bool isposmin=false;
  JSpaceEParmsPos ps;
  switch(key){
    case 'x':  if(CheckPosValue(Posmaxx,isposmin,ps))RunException(met,fun::PrintStr("The value posmax.x=\"%s\" is invalid.",Posminx.c_str()));  break;
    case 'y':  if(CheckPosValue(Posmaxy,isposmin,ps))RunException(met,fun::PrintStr("The value posmax.y=\"%s\" is invalid.",Posminy.c_str()));  break;
    case 'z':  if(CheckPosValue(Posmaxz,isposmin,ps))RunException(met,fun::PrintStr("The value posmax.z=\"%s\" is invalid.",Posminz.c_str()));  break;
    default:   RunException(met,fun::PrintStr("Key value \'%c\' is invalid.",key));
  }
  return(ps);
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
/// Checks format error and returns error and results.
//==============================================================================
int JSpaceEParms::CheckPosValue(const std::string &value,bool isposmin,JSpaceEParmsPos &ps)const{
  int error=0;
  ps.mode=DC_Default;
  ps.value=0;
  ps.textmod="";
  string v=fun::StrWithoutChar(value,' ');
  v=fun::StrWithoutChar(v,'\t');
  v=fun::StrLower(v);
  if(v=="default")ps.textmod=v;
  else{
    const int posdef=max(int(v.find("default+")),int(v.find("default-")));
    const int posprc=int(v.find("%"));
    if(posdef<0 && posprc>=0)error=2;                       //-The use of %
    if(!error && posdef>0)error=1;                          //-The use of default
    if(!error && posprc>=0 && posprc!=v.length()-1)error=2; //-The use of %
    if(!error && posdef==0 && v.length()<=string("default+").length())error=1; //-The use of default
    if(!error && posprc>=0 && posprc!=v.length()-1)error=2; //-The use of %
    if(!error){
      if(posdef==0){
        const char sign=v[7];
        const bool prc=(posprc>=0);
        if((isposmin && sign=='+') || (!isposmin && sign=='-'))error=4; //-The use of +/-
        ps.mode=(prc? DC_DefPrc: DC_DefValue);
        v=v.substr(8,v.length()-(prc? 9: 8));
        ps.value=atof(v.c_str());
        ps.textmod=fun::PrintStr("default %c %g",sign,ps.value);
        if(prc)ps.textmod=ps.textmod+"%";
      }
      else{
        ps.mode=DC_Fixed;
        ps.value=atof(v.c_str());
        ps.textmod=v;
      }
      if(!fun::StrIsRealNumber(v))error=3;  //-Invalid number.
    }
  }
  //printf("--------- [%s] --> [%s] --> [%s]\n",value.c_str(),v.c_str(),value2.c_str());
  return(error);
}

//==============================================================================
/// Checks format error and returns an improved string.
//==============================================================================
std::string JSpaceEParms::ReadPosValue(JXml *sxml,TiXmlElement* ele
  ,const std::string &name,const std::string &subname)const
{
  const char met[]="ReadPosValue";
  const string value=sxml->ReadElementStr(ele,name,subname,true,"default");
  const bool isposmin=(name=="posmin");
  JSpaceEParmsPos ps;
  int error=CheckPosValue(value,isposmin,ps);
  if(error){
    const string errtext=fun::PrintStr("is invalid in %s=\"%s\"",(name+"."+subname).c_str(),value.c_str());
    const string errfile=sxml->ErrGetFileRow(ele);
    if(error==1)RunException(met,fun::PrintStr("The use of default %s",errtext.c_str()),errfile);
    else if(error==2)RunException(met,fun::PrintStr("The use of %% %s",errtext.c_str()),errfile);
    else if(error==3)RunException(met,fun::PrintStr("Number %s",errtext.c_str()),errfile);
    else if(error==4){
      if(isposmin)RunException(met,fun::PrintStr("The sign + %s to increase the domain.",errtext.c_str()),errfile);
      else        RunException(met,fun::PrintStr("The sign - %s to increase the domain.",errtext.c_str()),errfile);
    }
    else RunException(met,fun::PrintStr("Format error %s",errtext.c_str()),errfile);
  }
  return(ps.textmod);
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
  if(sxml->ExistsElement(lis,"simulationdomain")){
    TiXmlElement* elep=lis->FirstChildElement("simulationdomain"); 
    Posminx=ReadPosValue(sxml,elep,"posmin","x");
    Posminy=ReadPosValue(sxml,elep,"posmin","y");
    Posminz=ReadPosValue(sxml,elep,"posmin","z");
    Posmaxx=ReadPosValue(sxml,elep,"posmax","x");
    Posmaxy=ReadPosValue(sxml,elep,"posmax","y");
    Posmaxz=ReadPosValue(sxml,elep,"posmax","z");
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
  //-Defines simulation domain.
  TiXmlElement *dom=sxml->AddElementAttrib(lis,"simulationdomain","comment","Defines domain of simulation (default=Uses minimun and maximum position of the generated particles)");
  TiXmlElement *pmin=sxml->AddElementAttrib(dom,"posmin","x",Posminx);
  sxml->AddAttribute(pmin,"y",Posminy);
  sxml->AddAttribute(pmin,"z",Posminz);
  sxml->AddAttribute(pmin,"comment","e.g.: x=0.5, y=default-1, z=default-10%");
  TiXmlElement *pmax=sxml->AddElementAttrib(dom,"posmax","x",Posmaxx);
  sxml->AddAttribute(pmax,"y",Posmaxy);
  sxml->AddAttribute(pmax,"z",Posmaxz);
}


