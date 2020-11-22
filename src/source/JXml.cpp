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

/// \file JXml.cpp \brief Implements the class \ref JXml.

#include "JXml.h"
#include "Functions.h"

#ifndef DISABLE_NUMEXLIB
#include "JNumexLib.h"
#endif

#include <climits>
#include <ctime>
#include <sys/stat.h>

#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

//==============================================================================
/// Constructor of objects.
//==============================================================================
JXml::JXml(){
  ClassName="JXml";
  Doc=new TiXmlDocument;
  Reset();
}
//==============================================================================
/// Destructor of objects.
//==============================================================================
JXml::~JXml(){
  DestructorActive=true;
  Reset();
  delete Doc;
}
//==============================================================================
/// Reinitialises the object state, recovering its configuration by default.
//==============================================================================
void JXml::Reset(){
  Doc->Clear();
  FileReading="";
  NuxLib=NULL;
}

//==============================================================================
/// Returns the requested node and creates it if necessary.
/// \param path Path of the requested node.
/// \param createpath Allows to create the path if does not exist,
/// otherwise returns \a NULL.
//==============================================================================
TiXmlNode* JXml::GetNode(const std::string &path,bool createpath){
  std::string pathx=path;
  TiXmlNode* node=NULL;
  TiXmlNode* base=Doc;
  while(pathx.length() && base){
    int pos=int(pathx.find("."));
    std::string txword=(pos>0? pathx.substr(0,pos): pathx);
    pathx=(pos>0? pathx.substr(pos+1): "");
    if(txword.length()){
      node=base->FirstChild(txword.c_str());
      if(!node&&createpath){
        if(txword=="case")base->InsertEndChild(TiXmlDeclaration("1.0","UTF-8",""));
        base->InsertEndChild(TiXmlElement(txword.c_str()));
        node=base->FirstChild(txword.c_str());
      }
      base=node;
    }
  }
  return(node);
}

//==============================================================================
/// Returns the requested node.
/// \param path Path of the requested node.
  /// \param checkactive Checks if node is active.
/// otherwise returns \a NULL.
//==============================================================================
TiXmlNode* JXml::GetNodeSimple(const std::string &path,bool checkactive)const{
  std::string pathx=path;
  TiXmlNode* node=NULL;
  TiXmlNode* base=Doc;
  while(pathx.length() && base){
    int pos=int(pathx.find("."));
    std::string txword=(pos>0? pathx.substr(0,pos): pathx);
    pathx=(pos>0? pathx.substr(pos+1): "");
    if(txword.length()){
      node=base->FirstChild(txword.c_str());
      base=node;
    }
  }
  //Checks node is actived.
  if(node && checkactive && !GetAttributeBool(node->ToElement(),"active",true,true))node=NULL;
  return(node);
}
//==============================================================================
/// Returns the requested node and if does not exist an exception is thrown.
/// \throw The requested node does not exist.
//==============================================================================
TiXmlNode* JXml::GetNodeError(const std::string &path){
  TiXmlNode* node=GetNode(path,false);
  if(!node){
    std::string tex="Error reading xml - can not find the element \'";
    tex=tex+path+"\'";
    Run_ExceptioonFile(tex,ErrGetFileRow(Doc));
  }
  return(node);
}

//==============================================================================
/// Returns the first requested element of a node TiXmlNode.
/// \param node Xml node where the reach is performed.
/// \param name Name of filtered elements (no filter using "").
/// \param optional If it does not exist,
/// returns \a NULL instead of throwing an exception.
/// \throw JException The requested element does not exist...
//==============================================================================
TiXmlElement* JXml::GetFirstElement(const TiXmlNode* node,const std::string &name,bool optional)const{
  TiXmlElement* ret=(TiXmlElement*)node->FirstChildElement(name.c_str());
  if(!ret && !optional)ErrReadElement(node,name,true);
  return(ret);
}

//==============================================================================
/// Returns the next requested element of a node TiXmlNode.
/// \param node Xml node where the reach is performed.
/// \param name Name of filtered elements (no filter using "").
/// \param optional If it does not exist,
/// returns \a NULL instead of throwing an exception.
/// \throw JException The requested element does not exist...
//==============================================================================
TiXmlElement* JXml::GetNextElement(const TiXmlNode* node,const std::string &name,bool optional)const{
  const TiXmlElement* ret=node->NextSiblingElement(name.c_str());
  if(!ret && !optional)ErrReadElement(node,name,true);
  return((TiXmlElement*)ret);
}

//==============================================================================
/// Returns the number of elements with a requested name of a node TiXmlNode.
/// \param node Xml node where the reach is performed.
/// \param name Name of filtered elements (no filter using "").
//==============================================================================
unsigned JXml::CountElements(const TiXmlNode* node,const std::string &name)const{
  unsigned count=0;
  TiXmlElement* ele=GetFirstElement(node,name,true);
  while(ele){
    count++;
    ele=GetNextElement(ele,name,true);
  }
  return(count);
}

//==============================================================================
/// Returns true when element is valid and it is not deactivated.
/// \param lis Xml element to look for requesed element name.
/// \param name Element name to look for.
//==============================================================================
bool JXml::CheckElementActive(const TiXmlElement* lis,const std::string &name)const{
  //return(lis? ReadElementBool(lis,name,"active",true,true): false);
  bool active=false;
  TiXmlElement* ele=(lis? GetFirstElement(lis,name,true): NULL);
  if(ele)active=GetAttributeBool(ele,"active",true,true);
  return(active);
}

//==============================================================================
/// Returns true when element is not deactivated.
/// \param ele Xml element to check.
//==============================================================================
bool JXml::CheckElementActive(const TiXmlElement* ele)const{
  return(ele? GetAttributeBool(ele,"active",true,true): false);
}

//==============================================================================
/// Removes the requested node.
/// \param path Path of the requested node.
//==============================================================================
void JXml::RemoveNode(const std::string &path){
  TiXmlNode* node=GetNode(path,false);
  if(node)node->Parent()->RemoveChild(node);
}

//==============================================================================
/// Returns the filename of the current xml with row of the requested node.
//==============================================================================
std::string JXml::ErrGetFileRow(const TiXmlNode* node)const{
  return(FileReading+(node? fun::PrintStr("(row:%d)",node->Row()): string("(row:?)")));
}

//==============================================================================
/// Returns the filename of the current xml with row of the first element in the node.
//==============================================================================
std::string JXml::ErrGetFileRow(const TiXmlNode* node,const std::string &firstelement)const{
  TiXmlNode* ele=GetFirstElement(node,firstelement,true);
  return(ErrGetFileRow(ele? ele: node));
}

//==============================================================================
/// Throws an exception with the xml node and the name of the element.
/// \param node Xml node of the error.
/// \param element Name of the element where the error appears.
/// \param missing Error because it does not exist.
/// \throw JException Error in element...
//==============================================================================
void JXml::ErrReadElement(const TiXmlNode* node,const std::string &element,bool missing,std::string errortext)const{
  std::string tex="Error reading xml - ";
  if(missing)tex=tex+"Some element is missing \'"+element+"\'";
  else{
    tex=tex+"Element \'"+element+"\' is invalid.";
    if(!errortext.empty())tex=tex+" "+errortext;
  }
  Run_ExceptioonFile(tex,ErrGetFileRow(node));
}
//==============================================================================
/// Throws an exception with the xml element and the name of the attribute.
/// \param ele Xml element of the error.
/// \param atrib Name of the attribute where the error appears.
/// \param missing Error because it does not exist.
/// \throw JException Error in element...
//==============================================================================
void JXml::ErrReadAtrib(const TiXmlElement* ele,const std::string &atrib,bool missing,std::string errortext)const{
  std::string tex="Error reading xml - ";
  if(missing)tex=tex+"Attribute \'"+atrib+"\' is missing.";
  else{
    tex=tex+"Value of \'"+atrib+"\' invalid.";
    if(!errortext.empty())tex=tex+" "+errortext;
  }
  Run_ExceptioonFile(tex,ErrGetFileRow(ele));
}
//==============================================================================
/// Throws an exception with the xml element and the name of the attribute.
/// \param ele Xml element of the error.
/// \param atrib Name of the attribute where the error appears.
/// \param missing Error because it does not exist.
/// \throw JException Error in element...
//==============================================================================
void JXml::ErrUnknownAtrib(const TiXmlElement* ele,const std::string &atrib)const{
  std::string tex="Error reading xml - ";
  tex=tex+"Attribute \'"+atrib+"\' is unknown.";
  Run_ExceptioonFile(tex,ErrGetFileRow(ele));
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//==============================================================================
/// Throws an exception if there are unknown or repeated elements.
/// \param lis Xml element to check.
/// \param names List of valid names (separated by spaces and *name for repeatable names).
/// \param checkrepeated Checks if there are repeated elements.
//==============================================================================
void JXml::CheckElementNames(const TiXmlElement* lis,bool checkrepeated,std::string names)const{
  //-Create list of elements.
  vector<string> vnames;
  vector<string> vnamesr;
  const unsigned nv=fun::VectorSplitStr(" ",fun::StrTrimRepeated(names),vnames);
  for(unsigned c=0;c<nv;c++){
    if(vnames[c][0]=='*'){
      vnames[c]=vnames[c].substr(1);
      vnamesr.push_back(vnames[c]);
    }
  }
  vector<string> vused;
  const TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    std::string ename=ele->Value();
    if(!ename.empty() && ename[0]!='_'){
      if(fun::VectorFind(ename,vnames)==UINT_MAX)ErrReadElement(ele,ename,false); //-Error: unknown element.
      if(checkrepeated){
        const bool rname=(fun::VectorFind(ename,vnamesr)!=UINT_MAX); //-Repeatable element.
        if(!rname && fun::VectorFind(ename,vused)!=UINT_MAX)ErrReadElement(ele,ename,false); //-Error: repeated element.
        if(!rname)vused.push_back(ename);
      }
    }
    ele=ele->NextSiblingElement();
  }
  //std::string listenames=" ";
  //names=std::string(" ")+names+" ";
  //TiXmlElement* ele=lis->FirstChildElement(); 
  //while(ele){
  //  std::string ename=ele->Value();
  //  if(!ename.empty() && ename[0]!='_'){
  //    //printf("++> ename:[%s]\n",ename.c_str());
  //    if(int(names.find(std::string(" ")+ename+" "))<0)ErrReadElement(ele,ename,false); //-Error elemento desconocido.
  //    if(checkrepeated){
  //      if(int(listenames.find(std::string(" ")+ename+" "))>=0)ErrReadElement(ele,ename,false); //-Error elemento repetido.
  //      listenames=listenames+ename+" ";
  //    }
  //  }
  //  ele=ele->NextSiblingElement();
  //}
}

//==============================================================================
/// Checks if some or several attributes appers in the element. Returns number
/// of found attribute (1...n), 0 none found and -1 several found.
/// \param ele Xml element of the error.
/// \param names Names of the requested attributes separated by by spaces.
/// \param checkmanyatt Throw exception if several attributes exist.
/// \param checkmanyele Throw exception if several elements exist.
//==============================================================================
int JXml::CheckElementAttributes(const TiXmlElement* ele,const std::string &name,std::string attnames,bool checkmanyatt,bool checkmanyele)const{
  if(checkmanyele && CountElements(ele,name)>1)Run_ExceptioonFile(string("Element \'"+name+"\' appears several times."),ErrGetFileRow(ele));
  TiXmlElement* ele2=GetFirstElement(ele,name,true); 
  if(ele2)return(CheckAttributes(ele2,attnames,checkmanyatt));
  else return(0);
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//==============================================================================
/// Throws an exception if there are unknown attribute.
/// \param ele Xml element to check.
/// \param names List of valid names (separated by spaces).
/// \param checkrepeated Checks if there are repeated elements.
//==============================================================================
void JXml::CheckAttributeNames(TiXmlElement* ele,std::string names)const{
  //-Create list of elements.
  vector<string> vnames;
  const unsigned nv=fun::VectorSplitStr(" ",fun::StrTrimRepeated(names),vnames);
  vector<string> vused;
  TiXmlAttribute* att=ele->FirstAttribute();
  while(att){
    std::string ename=att->Name();
    if(!ename.empty() && ename[0]!='_'){
      if(fun::VectorFind(ename,vnames)==UINT_MAX)ErrUnknownAtrib(ele,ename); //-Error: unknown attribute.
    }
    att=att->Next();
  }
}

//==============================================================================
/// Throws an exception if there are unknown attribute.
/// \param lis Xml element to look for element name.
/// \param names List of valid names (separated by spaces).
/// \param checkrepeated Checks if there are repeated elements.
//==============================================================================
void JXml::CheckAttributeNames(TiXmlElement* lis,std::string elementname,std::string names)const{
  TiXmlElement* ele=GetFirstElement(lis,elementname,true);
  if(ele)CheckAttributeNames(ele,names);
}

//==============================================================================
/// Checks if the requested attribute of an element already exists.
/// \param ele Xml element of the error.
/// \param name Name of the requested attribute.
//==============================================================================
bool JXml::ExistsAttribute(const TiXmlElement* ele,const std::string &name)const{
  return(ele->Attribute(name.c_str())!=NULL);
}

//==============================================================================
/// Checks if some or several attributes appers in the element. Returns number
/// of found attribute (1...n), 0 none found and -1 several found.
/// \param ele Xml element of the error.
/// \param names Names of the requested attributes separated by spaces.
/// \param checkmanyatt Throw exception if several attributes exist.
//==============================================================================
int JXml::CheckAttributes(const TiXmlElement* ele,std::string names,bool checkmanyatt)const{
  int ret=0;
  string severaldefs;
  for(int c=1;!names.empty() && ret!=-1;c++){
    string name=fun::StrSplit(" ",names);
    if(ExistsAttribute(ele,name)){
      severaldefs=severaldefs+(severaldefs.empty()? "": ", ")+name;
      ret=(ret? -1: c);
    }
  }
  if(checkmanyatt && ret==-1)Run_ExceptioonFile(string("Several definitions (")+severaldefs+") for \'"+ele->Value()+"\'.",ErrGetFileRow(ele));
  return(ret);
}

//==============================================================================
/// Checks if some or several attributes appers in the element. Returns number
/// of found attribute (1...n), 0 none found and -1 several found.
/// \param lis List of Xml elements of the error.
/// \param elementname Name of the requested element.
/// \param names Names of the requested attributes separated by spaces.
/// \param checkmanyatt Throw exception if several attributes exist.
//==============================================================================
int JXml::CheckAttributes(const TiXmlElement* lis,std::string elementname
  ,std::string names,bool checkmanyatt)const
{
  const TiXmlElement *ele=lis->FirstChildElement(elementname.c_str());
  if(ele)return(CheckAttributes(ele,names,checkmanyatt));
  return(0);
}


//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//==============================================================================
/// Checks and returns a value of type string of an xml element without any 
/// special treatment.
/// \param ele Xml element.
/// \param name Name of the requested attribute.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
std::string JXml::GetAttributeStrSimple(const TiXmlElement* ele,const std::string &name
  ,bool optional,const std::string &valdef)const
{
  std::string ret;
  if(ele->Attribute(name.c_str()))ret=ele->Attribute(name.c_str());
  else if(optional)ret=valdef;
  else ErrReadAtrib(ele,name,true);
  return(ret);
}

//==============================================================================
/// Checks and returns a value of type string of an xml element.
/// \param ele Xml element.
/// \param name Name of the requested attribute.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
std::string JXml::GetAttributeStr(const TiXmlElement* ele,const std::string &name
  ,bool optional,const std::string &valdef)const
{
  std::string ret;
  if(ele->Attribute(name.c_str()))ret=ele->Attribute(name.c_str());
  else if(optional)ret=valdef;
  else ErrReadAtrib(ele,name,true);
  //-String should be evaluated.
  if(NuxLib!=NULL && !ret.empty() && ret[0]=='$'){
    #ifndef DISABLE_NUMEXLIB
      ret=NuxLib->ComputeExprStr(ret,std::string("\nFile: ")+ErrGetFileRow(ele));
    #endif
  }
  //-Checks first letter.
  if(!ret.empty() && ret[0]=='$')ErrReadAtrib(ele,name,false,"It is not a valid string since it starts with '$' and should be evaluated.");
  return(ret);
}

//==============================================================================
/// Checks and returns a value of type bool of an xml element, 
/// valid values in xml text are: \a '0', \a '1', \a 'false' y \a 'true'.
/// \param ele Xml element.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
bool JXml::GetAttributeBool(const TiXmlElement* ele,const std::string &name,bool optional,bool valdef)const{
  bool ret=false;
  if(ele->Attribute(name.c_str())){
    const string val=ele->Attribute(name.c_str());
    if(NuxLib!=NULL && !val.empty() && val[0]=='#'){
      #ifndef DISABLE_NUMEXLIB
        ret=int(NuxLib->ComputeExprBool(val,std::string("\nFile: ")+ErrGetFileRow(ele)));
      #endif
    }
    else if(fun::StrLower(val)=="true"  || val=="1")ret=true;
    else if(fun::StrLower(val)=="false" || val=="0")ret=false;
    else ErrReadAtrib(ele,name,false);
  }
  else if(optional)ret=valdef;
  else ErrReadAtrib(ele,name,true);
  return(ret);
}

//==============================================================================
/// Checks and returns a value of type byte of an xml element that must be (0-225).
/// \param ele Xml element.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
byte JXml::GetAttributeByte(TiXmlElement* ele,const std::string &name,bool optional,byte valdef)const{
  int ret=GetAttributeInt(ele,name,optional,int(valdef));
  if(ret<0||ret>255)ErrReadAtrib(ele,name,false);
  return(byte(ret));
}

//==============================================================================
/// Checks and returns a value of type word of an xml element that must be (0-65535).
/// \param ele Xml element.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
word JXml::GetAttributeWord(TiXmlElement* ele,const std::string &name,bool optional,word valdef)const{
  int ret=GetAttributeInt(ele,name,optional,int(valdef));
  if(ret<0||ret>65535)ErrReadAtrib(ele,name,false);
  return(word(ret));
}

//==============================================================================
/// Checks and returns a value of type unsigned of an xml element that must be positive.
/// \param ele Xml element.
/// \param name Name of the requested attribute.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
unsigned JXml::GetAttributeUnsigned(TiXmlElement* ele,const std::string &name,bool optional,unsigned valdef)const{
  int ret=GetAttributeInt(ele,name,optional,int(valdef));
  double retdbl=GetAttributeDouble(ele,name,optional,double(valdef));
  if(retdbl<0)ErrReadAtrib(ele,name,false);
  return(unsigned(ret));
}

//==============================================================================
/// Checks and returns a value of type int of an xml element. 
/// \param ele Xml element.
/// \param name Name of the requested attribute.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
int JXml::GetAttributeInt(TiXmlElement* ele,const std::string &name,bool optional,int valdef)const{
  int ret;
  const char *vchar=ele->Attribute(name.c_str());
  if(vchar==NULL){
    if(optional)ret=valdef;
    else ErrReadAtrib(ele,name,true);
  }
  else{
    if(NuxLib!=NULL && vchar[0]=='#'){
      #ifndef DISABLE_NUMEXLIB
        ret=int(NuxLib->ComputeExpr(vchar,std::string("\nFile: ")+ErrGetFileRow(ele)));
      #endif
    }
    else{
      if(!fun::StrIsIntegerNumber(vchar))ErrReadAtrib(ele,name,false,"It is not a valid integer number.");
      ret=atoi(vchar);
    }
  }
  return(ret);
}

//==============================================================================
/// Checks and returns a value of type double of an xml element. 
/// \param ele Xml element.
/// \param name Name of the requested attribute.
/// \param optional If it does not exist,
/// returns \a valdef instead of throwing an exception.
/// \param valdef Value by default if it does not exist and \a optional was activated. 
/// \throw JException The requested attribute does not exist...
//==============================================================================
double JXml::GetAttributeDouble(TiXmlElement* ele,const std::string &name
  ,bool optional,double valdef)const
{
  double ret;
  const char *vchar=ele->Attribute(name.c_str());
  if(vchar==NULL){
    if(optional)ret=valdef;
    else ErrReadAtrib(ele,name,true);
  }
  else{
    if(NuxLib!=NULL && vchar[0]=='#'){
      #ifndef DISABLE_NUMEXLIB
        ret=NuxLib->ComputeExpr(vchar,std::string("\nFile: ")+ErrGetFileRow(ele));
      #endif
    }
    else{
      if(!fun::StrIsRealNumber(vchar))ErrReadAtrib(ele,name,false,"It is not a valid real number.");
      ret=atof(vchar);
    }
  }
  return(ret);
}

//==============================================================================
/// Loads a list of xml elements of type tfloat3 in an array 
/// and returns the number of loaded elements.
/// \param node Xml node the list of values is loaded from. 
/// \param name Name of the elements of type tfloat3.
/// \param vec Array where loaded values are stored.
/// \param count Maximum number of elements to load.
/// \param readcount If activated, throws an exception 
/// if there are less elements than the requested ones.
/// \throw JException Values missing or any value is not valid...
//==============================================================================
unsigned JXml::ReadArrayFloat3(const TiXmlNode* node,const std::string &name
  ,tfloat3 *vec,unsigned count,bool readcount)const
{
  unsigned rcount=0;
  if(count){
    TiXmlElement* ele=GetFirstElement(node,name,!readcount);
    for(unsigned c=0;c<count&&ele;c++){
      vec[c]=GetAttributeFloat3(ele); rcount++;
      ele=(rcount<count? GetNextElement(ele,name,!readcount): NULL);
    }
  }
  return(rcount);
}

//==============================================================================
/// Loads a list of xml elements of type double3 in an array 
/// and returns the number of loaded elements.
/// \param node Xml node the list of values is loaded from. 
/// \param name Name of the elements of type double3.
/// \param vec Array where loaded values are stored.
/// \param count Maximum number of elements to load.
/// \param readcount If activated, throws an exception 
/// if there are less elements than the requested ones.
/// \throw JException Values missing or any value is not valid...
//==============================================================================
unsigned JXml::ReadArrayDouble3(const TiXmlNode* node,const std::string &name
  ,tdouble3 *vec,unsigned count,bool readcount)const
{
  unsigned rcount=0;
  if(count){
    TiXmlElement* ele=GetFirstElement(node,name,!readcount);
    for(unsigned c=0;c<count&&ele;c++){
      vec[c]=GetAttributeDouble3(ele); rcount++;
      ele=(rcount<count? GetNextElement(ele,name,!readcount): NULL);
    }
  }
  return(rcount);
}

//==============================================================================
/// Loads a list of xml elements of type int3 in an array 
/// and returns the number of loaded elements.
/// \param node Xml node the list of values is loaded from. 
/// \param name Name of the elements of type int3.
/// \param vec Array where loaded values are stored.
/// \param count Maximum number of elements to load.
/// \param readcount If activated, throws an exception 
/// if there are less elements than the requested ones.
/// \throw JException Values missing or any value is not valid...
//==============================================================================
unsigned JXml::ReadArrayInt3(const TiXmlNode* node,const std::string &name
  ,tint3 *vec,unsigned count,bool readcount)const
{
  unsigned rcount=0;
  if(count){
    TiXmlElement* ele=GetFirstElement(node,name,!readcount);
    for(unsigned c=0;c<count&&ele;c++){
      vec[c]=GetAttributeInt3(ele); rcount++;
      ele=(rcount<count? GetNextElement(ele,name,!readcount): NULL);
    }
  }
  return(rcount);
}

//==============================================================================
/// Loads a matrix from xml and returns the number of loaded elements.
/// \param node Xml node the list of values is loaded from. 
/// \param name Name of the elements of type double.
/// \param nrow Number of rows (the maximum allowed is 9).
/// \param ncol Number of cols (the maximum allowed is 9).
/// \param data Memory pointer where loaded values are stored.
/// \param ndata Available memory size of data pointer.
/// \param optionalvalues If some value does not exist, valdef is used.
/// \param valdef Value by default if some value does not exist and \a optional was activated. 
/// \throw JException Size of data is not enough...
/// \throw JException Element is not found...
/// \throw JException Values missing or any value is not valid...
//==============================================================================
unsigned JXml::ReadMatrixDouble(const TiXmlNode* node,const std::string &name
  ,unsigned nrows,unsigned ncols,unsigned ndata,double *data
  ,bool optionalvalues,double valdef)const
{
  TiXmlElement* ele=GetFirstElement(node,name);
  if(nrows>9 || ncols>9)ErrReadElement(ele,name,false,"Number number of columns or rows greater than 9.");
  if(nrows*ncols>ndata)ErrReadElement(ele,name,false,"Size of data is not enough to store the requested matrix.");
  unsigned nvalues=0;
  //-Look for each value.
  for(unsigned cr=1;cr<=nrows;cr++)for(unsigned cc=1;cc<=ncols;cc++){
    unsigned n=0;
    double v=valdef;
    string value=fun::PrintStr("v%d%d",cr,cc);
    if(ExistsAttribute(ele,value)){
      v=GetAttributeDouble(ele,value);
      n++;
    }
    TiXmlElement* elerow=GetFirstElement(ele,"values",true);
    while(elerow){
      if(ExistsAttribute(elerow,value)){
        v=GetAttributeDouble(elerow,value);
        n++;
      }
      elerow=GetNextElement(elerow,"values",true);
    }
    if(!n && optionalvalues)ErrReadElement(ele,name,false,string("Value \'")+value+"\' is missing.");
    if(n>1)ErrReadElement(ele,name,false,string("Value \'")+value+"\' is repeated.");
    data[(cr-1)*ncols+(cc-1)]=v;
    if(n)nvalues++;
  }
  return(nvalues);
}




//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//==============================================================================
/// Returns a value of type double in text according to the requested format.
/// \param v Real value.
/// \param fmt Format to be converted into text (used by printf()).
//==============================================================================
std::string JXml::ToStr(double v,const char* fmt){
  char cad[256];
  sprintf(cad,fmt,v);
  return(cad);
}

//==============================================================================
/// Adds attribute of type double to an xml element.
/// \param ele Xml element.
/// \param attrib Name of the attribute.
/// \param v Value of the attribute.
/// \param fmt Format to be converted into text (used by printf()).
//==============================================================================
void JXml::AddAttribute(TiXmlElement* ele,const std::string &attrib,double v,const char* fmt){
  ele->SetAttribute(attrib.c_str(),ToStr(v,fmt).c_str());
}

//==============================================================================
/// Adds attribute of type bool to an xml element.
/// \param ele Xml element.
/// \param attrib Name of the attribute.
/// \param v Value of the attribute.
//==============================================================================
void JXml::AddAttribute(TiXmlElement* ele,const std::string &attrib,bool v){
  ele->SetAttribute(attrib.c_str(),ToStr(v).c_str());
}

//==============================================================================
/// Adds attribute of type int to an xml element.
/// \param ele Xml element.
/// \param attrib Name of the attribute.
/// \param v Value of the attribute.
//==============================================================================
void JXml::AddAttribute(TiXmlElement* ele,const std::string &attrib,int v){
  ele->SetAttribute(attrib.c_str(),v);
}

//==============================================================================
/// Adds attribute of type string to an xml element.
/// \param ele Xml element.
/// \param attrib Name of the attribute.
/// \param v Value of the attribute.
//==============================================================================
void JXml::AddAttribute(TiXmlElement* ele,const std::string &attrib,const std::string &v){
  ele->SetAttribute(attrib.c_str(),v.c_str());
}

//==============================================================================
/// Creates and returns an element starting from a value of type int3.
/// \param name Name of the element.
/// \param v Value of the element.
/// \param name1 Name of the first attribute (x by default).
/// \param name2 Name of the second attribute (y by default).
/// \param name3 Name of the third attribute (z by default).
//==============================================================================
TiXmlElement JXml::MakeElementInt3(const std::string &name,const tint3 &v,const char* name1,const char* name2,const char* name3){
  TiXmlElement item(name.c_str());
  AddAttribute(&item,name1,v.x);
  AddAttribute(&item,name2,v.y);
  AddAttribute(&item,name3,v.z);
  return(item);
}

//==============================================================================
/// Creates and returns an element starting from a value of type double3.
/// \param name Name of the element.
/// \param v Value of the element.
/// \param name1 Name of the first attribute (x by default).
/// \param name2 Name of the second attribute (y by default).
/// \param name3 Name of the third attribute (z by default).
//==============================================================================
TiXmlElement JXml::MakeElementDouble3(const std::string &name,const tdouble3 &v,const char* name1,const char* name2,const char* name3){
  TiXmlElement item(name.c_str());
  AddAttribute(&item,name1,v.x);
  AddAttribute(&item,name2,v.y);
  AddAttribute(&item,name3,v.z);
  return(item);
}

//==============================================================================
/// Creates and returns an element starting from a value of type double3.
/// \param name Name of the element.
/// \param nrow Number of rows.
/// \param ncol Number of cols.
/// \param values Values of the matrix.
//==============================================================================
TiXmlElement JXml::MakeElementMatrixDouble(const std::string &name
  ,unsigned nrows,unsigned ncols,const double* values)
{
  const std::string fmt=(nrows>9 || ncols>9? "v%d_%d": "v%d%d");
  TiXmlElement item(name.c_str());
  for(unsigned cr=1;cr<=nrows;cr++){
    TiXmlElement itemrow("values");
    for(unsigned cc=1;cc<=ncols;cc++){
      string value=fun::PrintStr(fmt.c_str(),cr,cc);
      AddAttribute(&itemrow,value,values[(cr-1)*ncols+(cc-1)]);
    }
    item.InsertEndChild(itemrow);
  }
  return(item);
}

//==============================================================================
/// Creates and returns an element with attribute of type double.
/// \param name Name of the element.
/// \param attrib Name of the attribute.
/// \param v Value of the attribute.
/// \param fmt Format to be converted into text (used by printf()).
//==============================================================================
TiXmlElement JXml::MakeElementAttrib(const std::string &name,const std::string &attrib,double v,const char* fmt){
  TiXmlElement item(name.c_str());
  item.SetAttribute(attrib.c_str(),ToStr(v,fmt).c_str());
  return(item);
}

//==============================================================================
/// Creates and returns an element with attribute of type int.
/// \param name Name of the element.
/// \param attrib Name of the attribute.
/// \param v Value of the attribute.
//==============================================================================
TiXmlElement JXml::MakeElementAttrib(const std::string &name,const std::string &attrib,int v){
  TiXmlElement item(name.c_str());
  item.SetAttribute(attrib.c_str(),v);
  return(item);
}

//==============================================================================
/// Creates and returns an element with attribute of type string.
/// \param name Name of the element.
/// \param attrib Name of the attribute.
/// \param v Value of the attribute.
//==============================================================================
TiXmlElement JXml::MakeElementAttrib(const std::string &name,const std::string &attrib,const std::string &v){
  TiXmlElement item(name.c_str());
  item.SetAttribute(attrib.c_str(),v.c_str());
  return(item);
}


//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
//==============================================================================
/// Stores the xml document in a file, including in the main node 
/// atributtes with the name of the application that generates it and when.
/// \param fname Filename.
/// \param app Name of the application that generates it
/// \param date Stores date and time of creation.
/// \throw JException Problems with file access...
//==============================================================================
void JXml::SaveFile(const std::string &fname,const std::string &app,bool date){
  if((!app.empty()||date)&&Doc->FirstChild()){
    TiXmlElement* ele=Doc->FirstChildElement();
    if(!app.empty())AddAttribute(ele,"app",app);
    if(date)AddAttribute(ele,"date",fun::GetDateTime());
  }
  if(!Doc->SaveFile(fname.c_str()))Run_ExceptioonFile("Cannot save the xml document.",fname);
  CorrectFile(fname);
}
//==============================================================================
/// Creates an xml document of a file.
/// \param fname Filename.
/// \throw JException Problems with file access...
//==============================================================================
void JXml::LoadFile(const std::string &fname){
  Reset();
  FileReading=fname;
  if(!Doc->LoadFile(FileReading.c_str())){
    std::string tex="Cannot load the xml file: ";
    tex=tex+Doc->ErrorDesc();
    char cad[256];
    sprintf(cad," (row:%d col:%d)",Doc->ErrorRow(),Doc->ErrorCol());
    Run_ExceptioonFile(tex+cad,FileReading);
  }
}

//==============================================================================
/// Correct symbols in file. (&gt; => '>' and &lt; => '<').
/// \param fname Filename.
/// \throw JException Problems with file access...
//==============================================================================
void JXml::CorrectFile(const std::string &fname){
  const unsigned sizemax=1024*1024*100;
  unsigned size=0;
  char *data=NULL;
  bool modif=false;
  //-Lee datos de fichero.
  {
    ifstream pf;
    pf.open(fname.c_str(),ios::binary);
    if(pf){
      pf.seekg(0,ios::end);
      size=unsigned(pf.tellg());
      if(size<=sizemax){
        pf.seekg(0,ios::beg);
        data=new char[size];
        pf.read(data,size);
      }
      else Run_ExceptioonFile("File size is larger than 100 Mb.",fname);
      pf.close();
    }
    else Run_ExceptioonFile("Cannot load the xml file.",fname);
    //-Procesa datos.
    bool open=false;
    unsigned cp=0;
    for(unsigned c=0;c<size;c++){
      char let=data[c];
      // &lt; <  &gt; >
      if(let=='\"')open=!open;
      else if(let=='&' && open && c+3<size && data[c+2]=='t' && data[c+3]==';'){
        if(data[c+1]=='l'){
          let='<'; c+=3; modif=true; 
        }
        if(data[c+1]=='g'){
          let='>'; c+=3; modif=true; 
        }
      }
      if(modif)data[cp]=let;
      cp++;
    }
    size=cp;
  }
  //-Graba nuevo fichero.
  if(modif){
    ofstream pf2;
    pf2.open(fname.c_str(),ios::binary);
    if(pf2){
      pf2.write(data,size);
      pf2.close();
    }
    else Run_ExceptioonFile("Cannot modify the xml file.",fname);
  }
  delete[] data;
}


