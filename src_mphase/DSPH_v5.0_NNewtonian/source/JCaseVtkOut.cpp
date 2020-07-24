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

/// \file JCaseVtkOut.cpp \brief Implements the class \ref JCaseVtkOut.

#include "JCaseVtkOut.h"
#include "Functions.h"
#include "JRangeFilter.h"
#include "JXml.h"
#include <climits>

using namespace std;

//##############################################################################
//# JCaseVtkOutFile
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JCaseVtkOutFile::JCaseVtkOutFile(const std::string &file,const std::string &mks):File(file){
  Reset();
  SetMks(mks);
}

//==============================================================================
// Destructor.
//==============================================================================
JCaseVtkOutFile::~JCaseVtkOutFile(){
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JCaseVtkOutFile::Reset(){
  ListMk="";
}

//==============================================================================
// Replace list of mk values.
//==============================================================================
void JCaseVtkOutFile::SetMks(const std::string &mks){
  JRangeFilter rg(mks);
  ListMk=rg.ToString();
}

//##############################################################################
//# JCaseVtkOut
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JCaseVtkOut::JCaseVtkOut(){
  ClassName="JCaseVtkOut";
  Reset();
}

//==============================================================================
// Destructor.
//==============================================================================
JCaseVtkOut::~JCaseVtkOut(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JCaseVtkOut::Reset(){
  for(unsigned c=0;c<Files.size();c++)delete Files[c];
  Files.clear();
}

//==============================================================================
// Configures MkBoundFirst and MkFluidFirst.
//==============================================================================
void JCaseVtkOut::ConfigMkFirst(word mkboundfirst,word mkfluidfirst){
  Reset();
  MkBoundFirst=mkboundfirst;
  MkFluidFirst=mkfluidfirst;
}

//==============================================================================
// Configures MkBoundFirst and MkFluidFirst.
//==============================================================================
void JCaseVtkOut::AddFile(const std::string &fname,const std::string &mks){
  unsigned cfile=GetByFileName(fname);
  if(cfile==UINT_MAX){
    JCaseVtkOutFile* pfile=new JCaseVtkOutFile(fname,mks);
    Files.push_back(pfile);
  }
  else Files[cfile]->SetMks(mks);
}

//==============================================================================
/// Returns position in Files (UINT_MAX when it was not found).
//==============================================================================
unsigned JCaseVtkOut::GetByFileName(std::string fname)const{
  unsigned ipos=0;
  fname=fun::StrLower(fname);
  for(;ipos<Count() && fun::StrLower(Files[ipos]->File)!=fname;ipos++);
  return(ipos<Count()? ipos: UINT_MAX);
}

//==============================================================================
/// Returns the requested file by number.
//==============================================================================
std::string JCaseVtkOut::GetFile(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("Number of requested file name is invalid.");
  return(Files[idx]->File);
}

//==============================================================================
/// Returns the requested file by number.
//==============================================================================
std::string JCaseVtkOut::GetFileListMk(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("Number of requested file info is invalid.");
  return(Files[idx]->GetMks());
}

//==============================================================================
/// Returns list of files with indicated key.
//==============================================================================
unsigned JCaseVtkOut::GetFiles(std::string key,std::vector<std::string> &list)const{
  list.clear();
  for(unsigned ipos=0;ipos<Count();ipos++){
    if(int(Files[ipos]->File.find(key))>=0)list.push_back(Files[ipos]->File);
  }
  return(unsigned(Files.size()));
}

//==============================================================================
/// Returns list of files with requested mk.
//==============================================================================
unsigned JCaseVtkOut::GetFilesByMk(bool bound,word mk,std::vector<std::string> &list)const{
  list.clear();
  for(unsigned ipos=0;ipos<Count();ipos++){
    JRangeFilter rg(Files[ipos]->GetMks());
    if(rg.CheckValue(mk))list.push_back(Files[ipos]->File);
  }
  return(unsigned(Files.size()));
}

//==============================================================================
// Replace list of mk values.
//==============================================================================
std::string JCaseVtkOut::GetListMkType(bool bound,const std::string &mks)const{
  std::vector<unsigned> values;
  JRangeFilter rg(mks);
  rg.GetValues(values);
  const unsigned nv=unsigned(values.size());
  string txmks;
  for(unsigned c=0;c<nv;c++){
    const word v=word(values[c]);
    if(bound  && v>=MkBoundFirst)txmks=txmks+","+fun::UintStr(v-MkBoundFirst);//-For mkbound.
    if(!bound && v>=MkFluidFirst && v<MkBoundFirst)txmks=txmks+","+fun::UintStr(v-MkFluidFirst);//-For mkfluid.
  }
  rg.Config(txmks);
  return(rg.ToString());
}

//==============================================================================
/// Reads particles information in XML format.
//==============================================================================
void JCaseVtkOut::ReadXml(const JXml *sxml,TiXmlElement* lis){
  const word mkboundfirst=(word)sxml->GetAttributeUnsigned(lis,"mkboundfirst");
  const word mkfluidfirst=(word)sxml->GetAttributeUnsigned(lis,"mkfluidfirst");
  ConfigMkFirst(mkboundfirst,mkfluidfirst);
  //-Loads files information.
  TiXmlElement* ele=lis->FirstChildElement("vtkfile");
  while(ele){
    const string name=sxml->GetAttributeStr(ele,"name");
    const string mks=sxml->GetAttributeStr(ele,"mk");
    AddFile(name,mks);
    ele=ele->NextSiblingElement("vtkfile");
  }
}

//==============================================================================
/// Writes information in XML format.
//==============================================================================
void JCaseVtkOut::WriteXml(JXml *sxml,TiXmlElement* lis)const{
  lis->Clear();
  JXml::AddAttribute(lis,"mkboundfirst",MkBoundFirst);
  JXml::AddAttribute(lis,"mkfluidfirst",MkFluidFirst);
  for(unsigned c=0;c<Count();c++){
    TiXmlElement item("vtkfile");
    JXml::AddAttribute(&item,"name",Files[c]->File);
    const string mks=Files[c]->GetMks();
    JXml::AddAttribute(&item,"mk",mks);
    string mklist=GetListMkType(true,mks);
    if(!mklist.empty())JXml::AddAttribute(&item,"mkbound",mklist);
    mklist=GetListMkType(false,mks);
    if(!mklist.empty())JXml::AddAttribute(&item,"mkfluid",mklist);
    lis->InsertEndChild(item)->ToElement();
  }
}

//==============================================================================
/// Loads data in XML format from a file.
//==============================================================================
void JCaseVtkOut::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path,false);
}

//==============================================================================
/// Loads information from the object XML.
//==============================================================================
void JCaseVtkOut::LoadXml(const JXml *sxml,const std::string &place,bool optional){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node && !optional)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(node)ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Stores information in the object XML.
//==============================================================================
void JCaseVtkOut::SaveXml(JXml *sxml,const std::string &place)const{
  WriteXml(sxml,sxml->GetNode(place,true)->ToElement());
}


