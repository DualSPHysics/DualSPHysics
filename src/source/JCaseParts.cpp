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

/// \file JCaseParts.cpp \brief Implements the class \ref JCaseParts.

#include "JCaseParts.h"
#ifdef JCaseParts_UseProps
  #include "JCaseProperties.h"
#endif
#include "JRangeFilter.h"
#include "Functions.h"
#include "JXml.h"
#ifdef JCaseParts_UseLinearv
  #include "JLinearValue.h"
#endif
#include <algorithm>
#include <cfloat>

using std::string;

//##############################################################################
//# JCasePartBlock
//##############################################################################
//==============================================================================
/// Updates attribute Props.
//==============================================================================
#ifdef JCaseParts_UseProps
void JCasePartBlock::UpdateProperty(){
  Props=Properties->GetPropertyMk(Mk);
}
#endif

//==============================================================================
/// Returns name for XML file.
//==============================================================================
std::string JCasePartBlock::GetNameXml()const{
  std::string ret;
  for(unsigned c=BaseNameSize;c<ClassName.length();c++)ret=ret+char(tolower(ClassName[c]));
  return(ret);
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JCasePartBlock::ReadXml(const JXml *sxml,TiXmlElement* ele){
  MkType=(word)sxml->GetAttributeUnsigned(ele,(Bound? "mkbound": "mkfluid")); 
  Begin=sxml->GetAttributeUnsigned(ele,"begin");
  Count=sxml->GetAttributeUnsigned(ele,"count");
  Props=sxml->GetAttributeStr(ele,"property",true);
}

//==============================================================================
/// Writes particles information in xml format.
//==============================================================================
TiXmlElement* JCasePartBlock::WriteXml(JXml *sxml,TiXmlElement* ele)const{
  TiXmlElement item(GetNameXml().c_str());
  JXml::AddAttribute(&item,(Bound? "mkbound": "mkfluid"),MkType);
  JXml::AddAttribute(&item,"mk",Mk);
  JXml::AddAttribute(&item,"begin",Begin);
  JXml::AddAttribute(&item,"count",Count);
  if(!GetProperty().empty())JXml::AddAttribute(&item,"property",GetProperty());
  return(ele->InsertEndChild(item)->ToElement());
}

#ifdef JCaseParts_UseProps
//==============================================================================
/// Returns number of values.
//==============================================================================
unsigned JCasePartBlock::GetValuesCount()const{
  return(Properties->GetValuesCount(Props));
}

//==============================================================================
/// Returns name of value indicated.
//==============================================================================
std::string JCasePartBlock::GetValueName(unsigned idx)const{
  return(Properties->GetValueName(Props,idx));
}

//==============================================================================
/// Returns value as string of value indicated.
//==============================================================================
std::string JCasePartBlock::GetValueStr(unsigned idx)const{
  return(Properties->GetValueStr(Props,idx));
}

//==============================================================================
/// returns if exists the value indicated.
//==============================================================================
bool JCasePartBlock::ExistsValue(std::string name)const{
  return(Properties->ExistsValue(Props,name));
}

//==============================================================================
/// Returns value as string of value indicated.
//==============================================================================
std::string JCasePartBlock::GetValueStr(std::string name)const{
  return(Properties->GetValueStr(Props,name));
}

//==============================================================================
/// Returns number of subvalues.
//==============================================================================
unsigned JCasePartBlock::GetSubValuesCount(unsigned idx)const{
  return(Properties->GetSubValuesCount(Props,idx));
}

//==============================================================================
/// Returns name of subvalue indicated.
//==============================================================================
std::string JCasePartBlock::GetSubValueName(unsigned idx,unsigned subidx)const{
  return(Properties->GetSubValueName(Props,idx,subidx));
}

//==============================================================================
/// Returns value as string of subvalue indicated.
//==============================================================================
std::string JCasePartBlock::GetSubValueStr(unsigned idx,unsigned subidx)const{
  return(Properties->GetSubValueStr(Props,idx,subidx));
}

//==============================================================================
/// returns if exists the subvalue indicated.
//==============================================================================
bool JCasePartBlock::ExistsSubValue(std::string name,std::string subname)const{
  return(Properties->ExistsSubValue(Props,name,subname));
}

//==============================================================================
/// Returns value as string of subvalue indicated.
//==============================================================================
std::string JCasePartBlock::GetSubValueStr(std::string name,std::string subname,bool optional,std::string valdef)const{
  if(optional && !Properties->ExistsSubValue(Props,name,subname))return(valdef);
  else return(Properties->GetSubValueStr(Props,name,subname));
}

//==============================================================================
/// Returns value as int of subvalue indicated.
//==============================================================================
int JCasePartBlock::GetSubValueInt(std::string name,std::string subname,bool optional,int valdef)const{   
  if(optional && !Properties->ExistsSubValue(Props,name,subname))return(valdef);
  else return(atoi(GetSubValueStr(name,subname).c_str()));
}

//==============================================================================
/// Returns value as double of subvalue indicated.
//==============================================================================
double JCasePartBlock::GetSubValueDouble(std::string name,std::string subname,bool optional,double valdef)const{   
  if(optional && !Properties->ExistsSubValue(Props,name,subname))return(valdef);
  else return(atof(GetSubValueStr(name,subname).c_str()));
}
#endif

//##############################################################################
//# JCasePartBlock_Moving
//##############################################################################
//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JCasePartBlock_Moving::ReadXml(const JXml *sxml,TiXmlElement* ele){
  JCasePartBlock::ReadXml(sxml,ele);
  RefMotion=sxml->GetAttributeUnsigned(ele,"refmotion"); 
}

//==============================================================================
/// Writes particles information in xml format.
//==============================================================================
TiXmlElement* JCasePartBlock_Moving::WriteXml(JXml *sxml,TiXmlElement* ele)const{
  ele=JCasePartBlock::WriteXml(sxml,ele);
  JXml::AddAttribute(ele,"refmotion",RefMotion);
  return(ele);
}


//##############################################################################
//# JCasePartBlock_Floating
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCasePartBlock_Floating::JCasePartBlock_Floating(const JCaseProperties* properties
  ,word mktype,unsigned begin,unsigned count,double massbody,double masspart
  ,const tdouble3& center,const tmatrix3d& inertia
  ,const tint3 &translationfree,const tint3 &rotationfree
  ,const tdouble3 &linvelini,const tdouble3 &angvelini
  ,const JLinearValue *linvel  ,const JLinearValue *angvel
  ,const JLinearValue *linforce,const JLinearValue *angforce)
  :JCasePartBlock(properties,TpPartFloating,"Floating",mktype,begin,count)
  ,Massbody(massbody),Masspart(masspart),Center(center),Inertia(inertia)
  ,TranslationFree(translationfree),RotationFree(rotationfree)
  ,LinearVelini(linvelini),AngularVelini(angvelini)
{
   TranslationFree=TInt3((!TranslationFree.x? 0: 1),(!TranslationFree.y? 0: 1),(!TranslationFree.z? 0: 1));
   RotationFree   =TInt3((!RotationFree.x? 0: 1),(!RotationFree.y? 0: 1),(!RotationFree.z? 0: 1));
  #ifdef JCaseParts_UseLinearv
   LinearVel =(linvel? new JLinearValue(*linvel): NULL); 
   AngularVel=(angvel? new JLinearValue(*angvel): NULL); 
   LinearForce =(linforce? new JLinearValue(*linforce): NULL); 
   AngularForce=(angforce? new JLinearValue(*angforce): NULL); 
  #else                                                  
   LinearVel=NULL;                                        
   AngularVel=NULL;
   LinearForce=NULL;                                        
   AngularForce=NULL;
  #endif
}
//==============================================================================
/// Constructor from XML data.
//==============================================================================
JCasePartBlock_Floating::JCasePartBlock_Floating(const JCaseProperties* properties
  ,const JXml *sxml,TiXmlElement* ele):JCasePartBlock(properties,TpPartFloating,"Floating")
{
  LinearVel=NULL;   AngularVel=NULL;
  LinearForce=NULL; AngularForce=NULL;
  ReadXml(sxml,ele);
}
//==============================================================================
/// Destructor.
//==============================================================================
JCasePartBlock_Floating::~JCasePartBlock_Floating(){
  DestructorActive=true;
 #ifdef JCaseParts_UseLinearv
  delete LinearVel;  LinearVel=NULL;
  delete AngularVel; AngularVel=NULL;
  delete LinearForce;  LinearForce=NULL;
  delete AngularForce; AngularForce=NULL;
 #endif
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JCasePartBlock_Floating::ReadXml(const JXml *sxml,TiXmlElement* ele){
  JCasePartBlock::ReadXml(sxml,ele);
  Massbody=sxml->ReadElementDouble(ele,"massbody","value");
  Masspart=sxml->ReadElementDouble(ele,"masspart","value");
  Center=sxml->ReadElementDouble3(ele,"center");
  //-Reads inertia data from double3 or tmatrix3d XML element.
  TiXmlElement *item=sxml->GetFirstElement(ele,"inertia");
  if(sxml->ExistsAttribute(item,"x") && sxml->ExistsAttribute(item,"y") && sxml->ExistsAttribute(item,"z")){
    tdouble3 v3=sxml->ReadElementDouble3(ele,"inertia");
    Inertia=TMatrix3d(0);
    Inertia.a11=v3.x; Inertia.a22=v3.y; Inertia.a33=v3.z;
  }
  else Inertia=sxml->ReadElementMatrix3d(ele,"inertia");
  //-Reads motion data.
  if(sxml->ExistsElement(ele,"translation") && sxml->ExistsElement(ele,"translationDOF"))
    Run_ExceptioonFile("Only \'translation\' or \'translationDOF\' must be defined.",sxml->ErrGetFileRow(ele));
  if(sxml->ExistsElement(ele,"rotation")    && sxml->ExistsElement(ele,"rotationDOF"))
    Run_ExceptioonFile("Only \'rotation\' or \'rotationDOF\' must be defined.",sxml->ErrGetFileRow(ele));
  TranslationFree=sxml->ReadElementInt3(ele,(sxml->ExistsElement(ele,"translation")? "translation": "translationDOF"),true,TInt3(1));
  RotationFree   =sxml->ReadElementInt3(ele,(sxml->ExistsElement(ele,"rotation")   ? "rotation":    "rotationDOF")   ,true,TInt3(1));
  TranslationFree=TInt3((!TranslationFree.x? 0: 1),(!TranslationFree.y? 0: 1),(!TranslationFree.z? 0: 1));
  RotationFree   =TInt3((!RotationFree.x? 0: 1),(!RotationFree.y? 0: 1),(!RotationFree.z? 0: 1));
  LinearVelini =sxml->ReadElementDouble3(ele,(sxml->ExistsElement(ele,"linearvelini" )? "linearvelini" : "velini"  ),true);
  AngularVelini=sxml->ReadElementDouble3(ele,(sxml->ExistsElement(ele,"angularvelini")? "angularvelini": "omegaini"),true);
  //-Reads imposed velocity.
#ifdef JCaseParts_UseLinearv
  if(sxml->ExistsElement(ele,"linearvel")){
    LinearVel=new JLinearValue(3,true);
    LinearVel->ReadXmlValues(sxml,ele,"linearvel","vel","time:x:y:z");
  }
  if(sxml->ExistsElement(ele,"angularvel")){
    AngularVel=new JLinearValue(3,true);
    AngularVel->ReadXmlValues(sxml,ele,"angularvel","vel","time:x:y:z");
  } 
  if(sxml->ExistsElement(ele,"linearforce")){
    LinearForce=new JLinearValue(3,false);
    LinearForce->ReadXmlValues(sxml,ele,"linearforce","force","time:x:y:z");
  }
  if(sxml->ExistsElement(ele,"angularforce")){
    AngularForce=new JLinearValue(3,false);
    AngularForce->ReadXmlValues(sxml,ele,"angularforce","force","time:x:y:z");
  } 
#endif
}

//==============================================================================
/// Writes particles information in xml format.
//==============================================================================
TiXmlElement* JCasePartBlock_Floating::WriteXml(JXml *sxml,TiXmlElement* ele)const{
  ele=JCasePartBlock::WriteXml(sxml,ele);
  sxml->AddAttribute(sxml->AddElementAttrib(ele,"massbody","value",Massbody),"units_comment","kg");
  sxml->AddAttribute(sxml->AddElementAttrib(ele,"masspart","value",Masspart),"units_comment","kg");
  sxml->AddAttribute(sxml->AddElementDouble3(ele,"center",Center),"units_comment","metres (m)");
  //-Writes inertia data from double3 or tmatrix3d XML element.
  if(!Inertia.a12 && !Inertia.a13 && !Inertia.a21 && !Inertia.a23 && !Inertia.a31 && !Inertia.a32){
    sxml->AddAttribute(sxml->AddElementDouble3(ele,"inertia",TDouble3(Inertia.a11,Inertia.a22,Inertia.a33)),"units_comment","kg*m^2");
  }
#ifndef SPACEPARTS_USE_INERTIA_OLD
  else sxml->AddAttribute(sxml->AddElementMatrix3d(ele,"inertia",Inertia),"units_comment","kg*m^2");
#else
  else sxml->AddAttribute(sxml->AddElementDouble3(ele,"inertia",TDouble3(Inertia.a11,Inertia.a22,Inertia.a33)),"units_comment","kg*m^2");
#endif
  //-Writes motion data.
  if(TranslationFree!=TInt3(1))JXml::AddAttribute(sxml->AddElementInt3(ele,"translationDOF",TranslationFree),"comment","Use 0 to restrict Degrees Of Freedon for translation (default=(1,1,1))");
  if(RotationFree   !=TInt3(1))JXml::AddAttribute(sxml->AddElementInt3(ele,"rotationDOF"   ,RotationFree   ),"comment","Use 0 to restrict Degrees Of Freedon for rotation (default=(1,1,1))");
  if(LinearVelini !=TDouble3(0))JXml::AddAttribute(sxml->AddElementDouble3(ele,"linearvelini" ,LinearVelini),"units_comment","m/s");
  if(AngularVelini!=TDouble3(0))JXml::AddAttribute(sxml->AddElementDouble3(ele,"angularvelini",AngularVelini),"units_comment","rad/s");
#ifdef JCaseParts_UseLinearv
  //-Writes imposed velocity.
  if(LinearVel )JXml::AddAttribute(LinearVel ->WriteXmlValues(sxml,ele,"linearvel" ,"vel","time:x:y:z"),"units_comment","m/s");
  if(AngularVel)JXml::AddAttribute(AngularVel->WriteXmlValues(sxml,ele,"angularvel","vel","time:x:y:z"),"units_comment","rad/s");
  //-Writes added force.
  if(LinearForce )JXml::AddAttribute(LinearForce ->WriteXmlValues(sxml,ele,"linearforce" ,"force","time:x:y:z"),"units_comment","N/s");
  if(AngularForce)JXml::AddAttribute(AngularForce->WriteXmlValues(sxml,ele,"angularforce","force","time:x:y:z"),"units_comment","N/s");
#endif
  return(ele);
}


//##############################################################################
//# JCaseParts
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCaseParts::JCaseParts(){
  ClassName="JCaseParts";
  #ifdef JCaseParts_UseProps
    Properties=new JCaseProperties;
  #endif
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JCaseParts::~JCaseParts(){
  DestructorActive=true;
  Reset();
  #ifdef JCaseParts_UseProps
    delete Properties; Properties=NULL;
  #endif
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JCaseParts::Reset(){
  for(unsigned c=0;c<Blocks.size();c++){
    switch(Blocks[c]->Type){
      case TpPartFixed:   { JCasePartBlock_Fixed    *ptr=(JCasePartBlock_Fixed   *)Blocks[c]; delete ptr; }break;
      case TpPartMoving:  { JCasePartBlock_Moving   *ptr=(JCasePartBlock_Moving  *)Blocks[c]; delete ptr; }break;
      case TpPartFloating:{ JCasePartBlock_Floating *ptr=(JCasePartBlock_Floating*)Blocks[c]; delete ptr; }break;
      case TpPartFluid:   { JCasePartBlock_Fluid    *ptr=(JCasePartBlock_Fluid   *)Blocks[c]; delete ptr; }break;
      default: Run_Exceptioon("Type of block is unknown.");
    }
  }
  Blocks.clear();
  Begin=0;
  LastType=TpPartFixed;
  SetMkFirst(0,0);
  Posmin=Posmax=TDouble3(0);
#ifdef JCaseParts_UseProps
  Properties->Reset();
#endif
}

//==============================================================================
/// Returns the number of particles of a given type.
//==============================================================================
unsigned JCaseParts::Count(TpParticles type)const{
  unsigned n=0;
  for(unsigned c=0;c<Blocks.size();c++)if(Blocks[c]->Type==type)n+=Blocks[c]->GetCount();
  return(n);
}

//==============================================================================
/// Returns the total number of particles.
//==============================================================================
unsigned JCaseParts::Count()const{
  unsigned n=0;
  for(unsigned c=0;c<Blocks.size();c++)n+=Blocks[c]->GetCount();
  return(n);
}

//==============================================================================
/// Returns the number of particles of the requested type.
//==============================================================================
unsigned JCaseParts::CountBlocks(TpParticles type)const{
  unsigned n=0;
  for(unsigned c=0;c<Blocks.size();c++)if(Blocks[c]->Type==type)n++;
  return(n);
}

//==============================================================================
/// Returns the requested block.
//==============================================================================
const JCasePartBlock& JCaseParts::GetBlock(unsigned pos)const{
  if(pos>=CountBlocks())Run_Exceptioon("The requested particles block is missing.");
  return(*(Blocks[pos]));
}

//==============================================================================
/// Returns the block with a given MK.
//==============================================================================
JCasePartBlock* JCaseParts::GetByMkType(bool bound,word mktype)const{
  JCasePartBlock* bk=NULL;
  for(unsigned c=0;c<Blocks.size()&&!bk;c++)if(IsBound(Blocks[c]->Type)==bound && Blocks[c]->GetMkType()==mktype)bk=Blocks[c];
  return(bk);
}

//==============================================================================
/// Returns TRUE when number of particles match.
//==============================================================================
bool JCaseParts::CheckNparticles(unsigned casenfixed,unsigned casenmoving,unsigned casenfloat,unsigned casenfluid)const{
  return(casenfixed==Count(TpPartFixed) && casenmoving==Count(TpPartMoving) && casenfloat==Count(TpPartFloating) && casenfluid==Count(TpPartFluid));
}

//==============================================================================
/// Returns TRUE when number of particles match.
//==============================================================================
bool JCaseParts::CheckNparticles(ullong casenfixed,ullong casenmoving,ullong casenfloat,ullong casenfluid)const{
  return(CheckNparticles(unsigned(casenfixed),unsigned(casenmoving),unsigned(casenfloat),unsigned(casenfluid)));
}

//==============================================================================
/// Checks and adds a new block of particles.
//==============================================================================
void JCaseParts::Add(JCasePartBlock* block){
  if(GetByMkType(block->Bound,block->GetMkType()))Run_Exceptioon("Cannot add a block with a existing mk.");
  if(block->Type<LastType)Run_Exceptioon("The block type is invalid after the last type added.");
  block->ConfigMk(block->Bound? MkBoundFirst: MkFluidFirst);
  Blocks.push_back(block);
  Begin+=block->GetCount();
  LastType=block->Type;
}

//==============================================================================
/// Loads data in XML format from a file.
//==============================================================================
void JCaseParts::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path);
}

//==============================================================================
/// Stores data in XML format into a file.
//==============================================================================
void JCaseParts::SaveFileXml(const std::string &file,const std::string &path,bool newfile)const{
  JXml jxml;
  if(!newfile)jxml.LoadFile(file);
  SaveXml(&jxml,path);
  jxml.SaveFile(file);
}

//==============================================================================
/// Loads particles information from the object XML.
//==============================================================================
void JCaseParts::LoadXml(const JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Stores particles information in the object XML.
//==============================================================================
void JCaseParts::SaveXml(JXml *sxml,const std::string &place)const{
  WriteXml(sxml,sxml->GetNode(place,true)->ToElement());
}

//==============================================================================
/// Reads particles information in XML format.
//==============================================================================
void JCaseParts::ReadXml(const JXml *sxml,TiXmlElement* lis){
  unsigned np=sxml->GetAttributeUnsigned(lis,"np");
  unsigned nb=sxml->GetAttributeUnsigned(lis,"nb");
  unsigned nbf=sxml->GetAttributeUnsigned(lis,"nbf");
  word mkboundfirst=(word)sxml->GetAttributeUnsigned(lis,"mkboundfirst");
  word mkfluidfirst=(word)sxml->GetAttributeUnsigned(lis,"mkfluidfirst");
  SetMkFirst(mkboundfirst,mkfluidfirst);
  //-Loads properties information.
#ifdef JCaseParts_UseProps
  Properties->Reset();
  TiXmlElement* eprops=lis->FirstChildElement("properties"); 
  if(eprops)Properties->ReadXml(sxml,eprops);
#endif
  //-Loads particles information.
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    std::string cmd=ele->Value();
    if(cmd.length()&&cmd[0]!='_'){
      if(cmd=="fixed")          Add(new JCasePartBlock_Fixed(Properties,sxml,ele));
      else if(cmd=="moving")    Add(new JCasePartBlock_Moving(Properties,sxml,ele));
      else if(cmd=="floating")  Add(new JCasePartBlock_Floating(Properties,sxml,ele));
      else if(cmd=="fluid")     Add(new JCasePartBlock_Fluid(Properties,sxml,ele));
      else if(cmd!="properties")sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
  Begin=Count();
  if(np!=Count()||nb!=np-Count(TpPartFluid)||nbf!=Count(TpPartFixed))Run_Exceptioon("The amount of particles does not match the header.");
  //-Checks property info in blocks.
#ifdef JCaseParts_UseProps
  for(unsigned c=0;c<Blocks.size();c++){
    if(Blocks[c]->GetProperty()!=Properties->GetPropertyMk(Blocks[c]->GetMk())){
      Run_Exceptioon(std::string("Property information of mk=")+fun::IntStr(Blocks[c]->GetMk())+" does not correspond to Links configuration.");
    }
  }
#endif
}

//==============================================================================
/// Writes particles information in XML format.
//==============================================================================
void JCaseParts::WriteXml(JXml *sxml,TiXmlElement* lis)const{
  unsigned np=Count();
  lis->Clear();
  JXml::AddAttribute(lis,"np",np);
  JXml::AddAttribute(lis,"nb",np-Count(TpPartFluid));
  JXml::AddAttribute(lis,"nbf",Count(TpPartFixed));
  JXml::AddAttribute(lis,"mkboundfirst",GetMkBoundFirst());
  JXml::AddAttribute(lis,"mkfluidfirst",GetMkFluidFirst());
  WriteXmlSummary(sxml,lis);
  for(unsigned c=0;c<Blocks.size();c++)Blocks[c]->WriteXml(sxml,lis);
#ifdef JCaseParts_UseProps
  if(Properties->GetPropertyCount()){
    TiXmlElement* eprops=JXml::AddElement(lis,"properties");
    Properties->WriteXml(sxml,eprops);
  }
#endif
}

//==============================================================================
/// Writes particle summary in XML file.
//==============================================================================
void JCaseParts::WriteXmlSummary(JXml *sxml,TiXmlElement* ele)const{
  //-Loads summary data.
  StSummaryData dat=GetSummaryData();
  const unsigned ntp=4;
  const string txtp[ntp]={"fixed","moving","floating","fluid"};
  //-Writes summary data in XML.
  TiXmlElement items("_summary");
  ele=ele->InsertEndChild(items)->ToElement();
  TiXmlElement positions("positions");
  JXml::AddAttribute(&positions,"units_comment","metres (m)");
  TiXmlElement* elep=ele->InsertEndChild(positions)->ToElement();
  sxml->AddElementDouble3(elep,"posmin",Posmin);
  sxml->AddElementDouble3(elep,"posmax",Posmax);
  for(unsigned c=0;c<ntp;c++)if(dat.np[c]){
    TiXmlElement item(txtp[c].c_str());
    JXml::AddAttribute(&item,"count",dat.np[c]);
    JXml::AddAttribute(&item,"id",fun::PrintStr("%u-%u",dat.idini[c],dat.idlast[c]));
    JXml::AddAttribute(&item,"mkcount",dat.nmk[c]);
    JXml::AddAttribute(&item,"mkvalues",dat.mklist[c]);
    ele->InsertEndChild(item)->ToElement();
  }
}

//==============================================================================
/// Adjusts the value of Mk according to boundfirst and fluidfirst.
//==============================================================================
void JCaseParts::SetMkFirst(word boundfirst,word fluidfirst){
  MkBoundFirst=boundfirst; MkFluidFirst=fluidfirst;
  for(unsigned c=0;c<Blocks.size();c++)Blocks[c]->ConfigMk(Blocks[c]->Bound? MkBoundFirst: MkFluidFirst);
}

//==============================================================================
/// Changes particle number of block.
//==============================================================================
void JCaseParts::SetBlockSize(unsigned pos,unsigned np){
  if(pos>=CountBlocks())Run_Exceptioon("The requested particles block is missing.");
  JCasePartBlock* bk=Blocks[pos];
  if(bk->GetCount()!=np){
    const unsigned nsum=np-bk->GetCount();
    bk->SetCount(np);
    //-Cambia el begin del resto de bloques.
    for(unsigned c=pos+1;c<CountBlocks();c++)Blocks[c]->SetBegin(Blocks[c]->GetBegin()+nsum);
  }
}

#ifdef JCaseParts_UseProps
//==============================================================================
/// Load information of properties and filter data.
//==============================================================================
void JCaseParts::LoadProperties(const JCaseProperties *props){
  //-Load data of properties.
  Properties->CopyFrom(props);
  if(Properties->GetPropertyCount()){
    //-Get list of used mk.
    std::string mks;
    for(unsigned c=0;c<CountBlocks();c++){
      const JCasePartBlock &block=GetBlock(c);
      mks=mks+fun::IntStr(block.GetMk())+",";
    }
    //-Filter data of mks.
    Properties->FilterMk(MkBoundFirst,MkFluidFirst,mks);
  }
  //-Update property info in blocks.
  for(unsigned c=0;c<Blocks.size();c++)Blocks[c]->UpdateProperty();
}
#endif

//==============================================================================
/// Returns the list of MKs of a given type.
//==============================================================================
std::string JCaseParts::GetMkList(TpParticles type)const{
  string list;
  for(unsigned c=0;c<Blocks.size();c++)if(Blocks[c]->Type==type)list=list+fun::UintStr(Blocks[c]->GetMk())+",";
  JRangeFilter rg(list);
  return(rg.ToString());
}

//==============================================================================
/// Returns summary of particle data.
//==============================================================================
JCaseParts::StSummaryData JCaseParts::GetSummaryData()const{
  const unsigned ntp=TPPARTICLES_COUNT;
  TpParticles types[ntp]={TpPartFixed,TpPartMoving,TpPartFloating,TpPartFluid};
  StSummaryData dat;
  unsigned idbegin=0;
  for(unsigned c=0;c<ntp;c++){
    dat.np[c]=Count(types[c]);
    dat.nmk[c]=CountBlocks(types[c]);
    dat.mklist[c]=GetMkList(types[c]);
    if(dat.np[c]){
      dat.idini[c]=idbegin; idbegin+=dat.np[c]; dat.idlast[c]=idbegin-1;
    }
    else dat.idini[c]=dat.idlast[c]=0;
  }
  return(dat);
}

//==============================================================================
/// Returns vector of strings with the summary of particles.
//==============================================================================
void JCaseParts::GetParticleSummary(std::vector<std::string> &out)const{
  //-Loads summary data.
  StSummaryData dat=GetSummaryData();
  const unsigned ntp=TPPARTICLES_COUNT;
  //-Configures output format.
  string fmt="%u  id:(%u-%u)   MKs:%u (%s)";
  const bool aligned=false;
  if(aligned){
    unsigned snp=0,snm=0,sid=0,sid2=0;
    for(unsigned c=0;c<ntp;c++){
      snp=std::max(snp,dat.np[c]);
      snm=std::max(snm,dat.nmk[c]);
      sid=std::max(sid,dat.idini[c]);
      sid2=std::max(sid2,dat.idlast[c]);
    }
    snp=unsigned(fun::UintStr(snp).length());
    snm=unsigned(fun::UintStr(snm).length());
    sid=unsigned(fun::UintStr(sid).length());
    sid2=unsigned(fun::UintStr(sid2).length());
    fmt=fun::PrintStr("%%%uu  id:(%%%uu-%%%uu)   MKs:%%%uu (%%s)",snp,sid,sid2,snm);
  }
  //-Generates output.
  string txtype[ntp]={"  Fixed....: ","  Moving...: ","  Floating.: ","  Fluid....: "};
  out.push_back("Particle summary:");
  for(unsigned c=0;c<ntp;c++){
    out.push_back(txtype[c]+(dat.np[c]? fun::PrintStr(fmt.c_str(),dat.np[c],dat.idini[c],dat.idlast[c],dat.nmk[c],dat.mklist[c].c_str()): "0"));
  }
  out.push_back("");
  const unsigned nb=dat.np[0]+dat.np[1]+dat.np[2],nt=nb+dat.np[3];
  const unsigned mb=dat.nmk[0]+dat.nmk[1]+dat.nmk[2],mt=mb+dat.nmk[3];
  out.push_back(fun::PrintStr("Total particles: %u (bound=%u (fx=%u mv=%u ft=%u) fluid=%u)",nt,nb,dat.np[0],dat.np[1],dat.np[2],dat.np[3]));
  out.push_back(fun::PrintStr("Total MK blocks: %u (bound=%u (fx=%u mv=%u ft=%u) fluid=%u)",mt,mb,dat.nmk[0],dat.nmk[1],dat.nmk[2],dat.nmk[3]));
}

//==============================================================================
/// Returns string with particles information.
//==============================================================================
void JCaseParts::GetParticlesInfo(std::vector<std::string> &out)const{
  out.push_back("");
  out.push_back("Particles data:");
  out.push_back("----------------");
  const unsigned nfixed=Count(TpPartFixed);
  const unsigned nmoving=Count(TpPartMoving);
  const unsigned nfloat=Count(TpPartFloating);
  const unsigned nbound=nfixed+nmoving+nfloat;
  const unsigned nfluid=Count(TpPartFluid);
  out.push_back(fun::PrintStr("Number of particles: %u",nbound+nfluid));
  out.push_back(fun::PrintStr("- Boundary: %u",nbound));
  const unsigned nc=CountBlocks();
  out.push_back(fun::PrintStr("  - Fixed: %u",nfixed));
  for(unsigned c=0;c<nc;c++)if(GetBlock(c).Type==TpPartFixed){
    const JCasePartBlock& bk=GetBlock(c);
    const unsigned count=bk.GetCount(),begin=bk.GetBegin();
    out.push_back(fun::PrintStr("      Mk[%d]: %u (%u-%u)",bk.GetMk(),count,begin,begin+count-1));
  }
  out.push_back(fun::PrintStr("  - Moving: %u",nmoving));
  for(unsigned c=0;c<nc;c++)if(GetBlock(c).Type==TpPartMoving){
    const JCasePartBlock& bk=GetBlock(c);
    const unsigned count=bk.GetCount(),begin=bk.GetBegin();
    out.push_back(fun::PrintStr("      Mk[%d]: %u (%u-%u)",bk.GetMk(),count,begin,begin+count-1));
  }
  out.push_back(fun::PrintStr("  - Floating: %u",nfloat));
  for(unsigned c=0;c<nc;c++)if(GetBlock(c).Type==TpPartFloating){
    const JCasePartBlock& bk=GetBlock(c);
    const unsigned count=bk.GetCount(),begin=bk.GetBegin();
    out.push_back(fun::PrintStr("      Mk[%d]: %u (%u-%u)",bk.GetMk(),count,begin,begin+count-1));
  }
  out.push_back(fun::PrintStr("- Fluid: %u",nfluid));
  for(unsigned c=0;c<nc;c++)if(GetBlock(c).Type==TpPartFluid){
    const JCasePartBlock& bk=GetBlock(c);
    const unsigned count=bk.GetCount(),begin=bk.GetBegin();
    out.push_back(fun::PrintStr("      Mk[%d]: %u (%u-%u)",bk.GetMk(),count,begin,begin+count-1));
  }
}

//==============================================================================
/// Prints particles information.
//==============================================================================
void JCaseParts::VisuParticlesInfo()const{
  std::vector<string> lines;
  GetParticlesInfo(lines);
  for(unsigned c=0;c<unsigned(lines.size());c++)printf("%s\n",lines[c].c_str());
}

//==============================================================================
/// Returns true when imposed velocity was defined for some floating body.
//==============================================================================
bool JCaseParts::UseImposedFtVel()const{
  bool usevel=false;
  for(unsigned c=0;c<CountBlocks() && !usevel;c++)if(GetBlock(c).Type==TpPartFloating){
    const JCasePartBlock_Floating &fblock=(const JCasePartBlock_Floating &)GetBlock(c);
    usevel=(fblock.GetLinearVel()!=NULL || fblock.GetAngularVel()!=NULL);
  }
  return(usevel);
}

//==============================================================================
/// Returns true when added force was defined for some floating body.
//==============================================================================
bool JCaseParts::UseAddedFtForce()const{
  bool useforce=false;
  for(unsigned c=0;c<CountBlocks() && !useforce;c++)if(GetBlock(c).Type==TpPartFloating){
    const JCasePartBlock_Floating &fblock=(const JCasePartBlock_Floating &)GetBlock(c);
    useforce=(fblock.GetLinearForce()!=NULL || fblock.GetAngularForce()!=NULL);
  }
  return(useforce);
}

//##############################################################################
//# JCasePartsGetMk
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCasePartsGetMk::JCasePartsGetMk(const JCaseParts *sparts,bool splitting):Splitting(splitting){
  ClassName="JCasePartsGetMk";
  MkRange=NULL; MkValue=NULL;
  Reset();
  Config(sparts);
}

//==============================================================================
/// Destructor.
//==============================================================================
JCasePartsGetMk::~JCasePartsGetMk(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JCasePartsGetMk::Reset(){
  delete[] MkRange; MkRange=NULL;
  delete[] MkValue; MkValue=NULL;
  MkCount=0;
  MkSplitting=0;
}

//==============================================================================
/// Config object.
//==============================================================================
void JCasePartsGetMk::Config(const JCaseParts *sparts){
  MkCount=sparts->CountBlocks();
  MkRange=new unsigned[MkCount];
  MkValue=new word[MkCount];
  for(unsigned c=0;c<sparts->CountBlocks();c++){
    const JCasePartBlock &block=sparts->GetBlock(c);
    if(Splitting && block.Type==TpPartFluid)MkSplitting=block.GetMk();
    MkValue[c]=block.GetMk();
    MkRange[c]=block.GetBegin()+block.GetCount();
  }
}


