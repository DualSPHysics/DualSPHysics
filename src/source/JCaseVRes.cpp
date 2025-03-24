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

/// \file JCaseVRes.cpp \brief Implements the class \ref JCaseVRes.

#include "JCaseVRes.h"
#include "JCaseEParms.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JXml.h"
#include "JSpVtkShape.h"

#ifdef JXml_UseNux
#include "JNumexLib.h"
#endif

#include <algorithm>
#include <cmath>
#include <cfloat>

using namespace std;

//##############################################################################
//# JCaseVResExtra
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCaseVResExtra::JCaseVResExtra(){
  ClassName="JCaseVResExtra";
  Reset();
}

//==============================================================================
/// Constructor for copy.
//==============================================================================
JCaseVResExtra::JCaseVResExtra(const JCaseVResExtra& src):JCaseVResExtra()
{
  *this=src;
}

//==============================================================================
/// Overload assignment operator.
//==============================================================================
JCaseVResExtra& JCaseVResExtra::operator=(const JCaseVResExtra& src){
  if(this!=&src){
    Keys=src.Keys;
    Values=src.Values;
  }
  return(*this);
}

//==============================================================================
/// Destructor.
//==============================================================================
JCaseVResExtra::~JCaseVResExtra(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JCaseVResExtra::Reset(){
  Keys.clear();
  Values.clear();
}

//==============================================================================
/// Adds value and throws exception when it already exists.
//==============================================================================
void JCaseVResExtra::AddValue(const std::string& key,const std::string& value){
  std::vector<string>::iterator it=find(Keys.begin(),Keys.end(),key);
  if(it!=Keys.end())Run_Exceptioon("The key is already stored.");
  Keys.push_back(key);
  Values.push_back(value);
}

//==============================================================================
/// Returns key accoding to the index.
//==============================================================================
std::string JCaseVResExtra::GetKey(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("The requested index is invalid.");
  return(Keys[idx]);
}

//==============================================================================
/// Returns value accoding to the index.
//==============================================================================
std::string JCaseVResExtra::GetValue(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("The requested index is invalid.");
  return(Values[idx]);
}

//==============================================================================
/// Returns index of key or UINT_MAX.
//==============================================================================
unsigned JCaseVResExtra::FindKey(const std::string& key)const{
  unsigned c=0;
  for(;c<Count() && Keys[c]!=key;c++);
  return(c>=Count()? UINT_MAX: c);
}

//==============================================================================
/// Returns value by key as string.
//==============================================================================
std::string JCaseVResExtra::GetValueStr(const std::string& key,bool optional
  ,std::string valdef)const
{
  const unsigned idx=FindKey(key);
  if(idx==UINT_MAX && !optional)Run_Exceptioon(fun::PrintStr(
    "Requested key \'%s\' is missing.",key.c_str()));
  return(idx==UINT_MAX? valdef: Values[idx]);
}

//==============================================================================
/// Returns value by key as int.
//==============================================================================
int JCaseVResExtra::GetValueInt(const std::string& key,bool optional
  ,int valdef)const
{
  int v=valdef;
  const unsigned idx=FindKey(key);
  if(idx==UINT_MAX && !optional)Run_Exceptioon(fun::PrintStr(
    "Requested key \'%s\' is missing.",key.c_str()));
  if(idx!=UINT_MAX){
    const string tx=Values[idx];
    if(!fun::StrIsIntegerNumber(tx))Run_Exceptioon(fun::PrintStr(
      "Value of key \'%s\' is not a valid integer number.",tx.c_str()));
    v=atoi(tx.c_str());
  }
  return(v);
}

//==============================================================================
/// Returns value by key as unsigned.
//==============================================================================
unsigned JCaseVResExtra::GetValueUint(const std::string& key,bool optional
  ,unsigned valdef)const
{
  unsigned v=valdef;
  const unsigned idx=FindKey(key);
  if(idx==UINT_MAX && !optional)Run_Exceptioon(fun::PrintStr(
    "Requested key \'%s\' is missing.",key.c_str()));
  if(idx!=UINT_MAX){
    const string tx=Values[idx];
    if(!fun::StrIsIntegerNumber(tx))Run_Exceptioon(fun::PrintStr(
      "Value of key \'%s\' is not a valid integer number.",tx.c_str()));
    v=unsigned(atoi(tx.c_str()));
  }
  return(v);
}

//==============================================================================
/// Returns value by key as double.
//==============================================================================
double JCaseVResExtra::GetValueDbl(const std::string& key,bool optional
  ,double valdef)const
{
  double v=valdef;
  const unsigned idx=FindKey(key);
  if(idx==UINT_MAX && !optional)Run_Exceptioon(fun::PrintStr(
    "Requested key \'%s\' is missing.",key.c_str()));
  if(idx!=UINT_MAX){
    const string tx=Values[idx];
    if(!fun::StrIsRealNumber(tx))Run_Exceptioon(fun::PrintStr(
      "Value of key \'%s\' is not a valid real number.",tx.c_str()));
    v=atof(tx.c_str());
  }
  return(v);
}

//==============================================================================
/// Prints list of values (for debug).
//==============================================================================
void JCaseVResExtra::Print(std::string caption)const{
  if(!caption.empty())printf("%s\n",caption.c_str());
  for(unsigned c=0;c<Count();c++){
    printf("%02u. [%s] = [%s]\n",c,Keys[c].c_str(),Values[c].c_str());
  }
}


//##############################################################################
//# JCaseVResBase
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCaseVResBase::JCaseVResBase(TpMRShape shape,int id,JCaseVResBase* parent
  ,bool is2d,double posy2d,double hdp,double dp,double bsizeh,double overlaph
  ,std::string filerow):Shape(shape),Id(id),Parent(parent),Is2D(is2d)
  ,Posy2D(posy2d),Hdp(hdp),Dp(dp),BuffSizeh(bsizeh),Overlaph(overlaph)
  ,FileRow(filerow)
{
  ClassName=std::string("JCaseVRes_")+GetSubName();
  SimReset();
  TrackingDisable();
  NpReset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JCaseVResBase::~JCaseVResBase(){
  DestructorActive=true;
  SubZones.clear();
}

//==============================================================================
/// Set Extra data object.
//==============================================================================
void JCaseVResBase::SetExtraData(const JCaseVResExtra& edata){
  Extra_Data=edata;
}

//==============================================================================
// Write extra section data in XML object.
//==============================================================================
void JCaseVResBase::WriteXmlExtra(JXml* sxml,TiXmlElement* ele)const{
  if(Extra_Data.Count()){
    ele=sxml->AddElement(ele,"extra");
    for(unsigned c=0;c<Extra_Data.Count();c++){
      string k=Extra_Data.GetKey(c);
      string v=Extra_Data.GetValue(c);
      vector<std::string> vk;
      const unsigned nk=fun::VectorSplitStr(".",k,vk);
      if(nk){
        TiXmlElement* ele2=ele;
        for(unsigned ck=0;ck+1<nk;ck++)if(!vk[ck].empty()){
          if(!sxml->ExistsElement(ele2,vk[ck]))ele2=sxml->AddElement(ele2,vk[ck]);
          else ele2=sxml->GetFirstElement(ele2,vk[ck],true);
        }
        sxml->AddAttribute(ele2,vk[nk-1],v);
      }
    }
  }
}

//==============================================================================
/// Initialization of number of particles.
//==============================================================================
void JCaseVResBase::NpReset(){
  NpSet(0,0,0,0);
}

//==============================================================================
/// Set number of particles.
//==============================================================================
void JCaseVResBase::NpSet(ullong npfixed,ullong npmoving,ullong npfloating
  ,ullong npfluid)
{
  NpFixed   =npfixed;
  NpMoving  =npmoving;
  NpFloating=npfloating;
  NpFluid   =npfluid;
}

//==============================================================================
/// Initialization of simulation domain variables.
//==============================================================================
void JCaseVResBase::SimReset(){
  SimUseParent=false;
  SimPosXml=false;
  SimPosmin[0]=SimPosmin[1]=SimPosmin[2]="";
  SimPosmax[0]=SimPosmax[1]=SimPosmax[2]="";
  SimPartsPos=false;
  SimPartsPosmin=SimPartsPosmax=TDouble3(DBL_MAX);
}

//==============================================================================
/// Reads simulation domain configuration from XML node.
//==============================================================================
void JCaseVResBase::SimReadXmlDef(const JXml* sxml,const TiXmlElement* node){
  SimReset();
  TiXmlElement* ele=sxml->GetFirstElement(node,"simulationdomain",true);
  if(ele){
    sxml->CheckElementNames(ele,true,"posmin posmax");
    sxml->CheckAttributeNames(ele,"useparent comment");
    SimUseParent=sxml->GetAttributeBool(ele,"useparent",true,false);
    if(!Id && SimUseParent)Run_ExceptioonFile("Option useparent=true is invalid for main domain."
      ,sxml->ErrGetFileRow(ele));
    for(unsigned c=0;c<3;c++){
      const string attname=(!c? "x": (c==1? "y": "z"));
      SimPosmin[c]=fun::StrLower(sxml->ReadElementStrSimple(ele,"posmin",attname,true,"default"));
      SimPosmax[c]=fun::StrLower(sxml->ReadElementStrSimple(ele,"posmax",attname,true,"default"));
      if(SimPosmin[c]!="default" || SimPosmax[c]!="default")SimPosXml=true;
      printf("posmin.%s=\"%s\"  posmax.%s=\"%s\"\n",attname.c_str(),SimPosmin[c].c_str(),attname.c_str(),SimPosmax[c].c_str());
    }
    //-Checks simulation domain configuration.
    JCaseEParms eparms;
    eparms.SetPosmin(SimPosmin[0],SimPosmin[1],SimPosmin[2],sxml->ErrGetFileRow(ele));
    eparms.SetPosmax(SimPosmax[0],SimPosmax[1],SimPosmax[2],sxml->ErrGetFileRow(ele));
    //-Updates style of text configuration.
    if(SimPosXml){
      SimPosmin[0]=eparms.GetPosminValue('x').textmod;
      SimPosmin[1]=eparms.GetPosminValue('y').textmod;
      SimPosmin[2]=eparms.GetPosminValue('z').textmod;
      SimPosmax[0]=eparms.GetPosmaxValue('x').textmod;
      SimPosmax[1]=eparms.GetPosmaxValue('y').textmod;
      SimPosmax[2]=eparms.GetPosmaxValue('z').textmod;
    }
  }
}
  
//==============================================================================
// Configures minimum and maximum position of generated particles.
//==============================================================================
void JCaseVResBase::SimConfigPartsPos(tdouble3 partpmin,tdouble3 partpmax){
  SimPartsPos=true;
  SimPartsPosmin=partpmin;
  SimPartsPosmax=partpmax;
  //printf(">> SimConfigPartsPos_%d: %s\n",Id,fun::Double3gRangeStr(SimPartsPosmin,SimPartsPosmax).c_str());
}

//==============================================================================
// Compute final value for simulation domain.
//==============================================================================
double JCaseVResBase::SimComputeValue(int mode,double vdef,double vsize
  ,double vmod)const
{
  double v=vdef;
  switch(JCaseEParms::TpDomainMode(mode)){
    case JCaseEParms::DC_Fixed:     v=vmod;                     break;
    case JCaseEParms::DC_DefValue:  v=vdef+vmod;                break;
    case JCaseEParms::DC_DefPrc:    v=vdef + vsize*(vmod/100);  break;
  }
  //printf(">> SimComputeValue(mode=%d,vdef=%f,vsize=%f,vmod=%f) --> %f\n",mode,vdef,vsize,vmod,v);
  return(v);
}

//==============================================================================
// Write simulation domain configuration for execution in XML object.
//==============================================================================
void JCaseVResBase::SimWriteXmlRun(JXml* sxml,TiXmlElement* ele)const{
  if(SimUseParent || SimPosXml){
    TiXmlElement item("simulationdomain");
    //-Write basic configuration.
    sxml->AddAttribute(&item,"comment","Defines simulation domain for subdomain");
    if(SimUseParent){
      if(!Parent || !Parent->SimPartsPos)Run_Exceptioon("Particle domain limits of parent are missing.");
      tdouble3 pmin=Parent->SimPartsPosmin;
      tdouble3 pmax=Parent->SimPartsPosmax;
      const tdouble3 size=pmax-pmin;
      //printf(">> Id:%d  pmin/pmax: %s\n",Id,fun::Double3gRangeStr(SimPartsPosmin,SimPartsPosmax).c_str());
      JCaseEParms eparms;
      eparms.SetPosmin(SimPosmin[0],SimPosmin[1],SimPosmin[2],"");
      eparms.SetPosmax(SimPosmax[0],SimPosmax[1],SimPosmax[2],"");
      JCaseEParms::JCaseEParmsPos epos;
      epos=eparms.GetPosminValue('x'); pmin.x=SimComputeValue(int(epos.mode),pmin.x,size.x,-epos.value);
      epos=eparms.GetPosminValue('y'); pmin.y=SimComputeValue(int(epos.mode),pmin.y,size.y,-epos.value);
      epos=eparms.GetPosminValue('z'); pmin.z=SimComputeValue(int(epos.mode),pmin.z,size.z,-epos.value);
      epos=eparms.GetPosmaxValue('x'); pmax.x=SimComputeValue(int(epos.mode),pmax.x,size.x, epos.value);
      epos=eparms.GetPosmaxValue('y'); pmax.y=SimComputeValue(int(epos.mode),pmax.y,size.y, epos.value);
      epos=eparms.GetPosmaxValue('z'); pmax.z=SimComputeValue(int(epos.mode),pmax.z,size.z, epos.value);
      sxml->AddElementDouble3(&item,"posmin",pmin);
      sxml->AddElementDouble3(&item,"posmax",pmax);
    }
    else if(SimPosXml){
      TiXmlElement posmin("posmin");
      TiXmlElement posmax("posmax");
      for(unsigned c=0;c<3;c++){
        const string attname=(!c? "x": (c==1? "y": "z"));
        sxml->AddAttribute(&posmin,attname,SimPosmin[c]);
        sxml->AddAttribute(&posmax,attname,SimPosmax[c]);
      }
      item.InsertEndChild(posmin);
      item.InsertEndChild(posmax);
    }
    TiXmlElement* itemold=sxml->GetFirstElement(ele,"simulationdomain",true);
    if(itemold)ele->ReplaceChild(itemold,item);
    else ele->InsertEndChild(item);
  }
}

//==============================================================================
/// Defines mk for tracking (moving or floating bodies).
//==============================================================================
void JCaseVResBase::TrackingDisable(){
  TrackingMk=TrackingMkBound=-1;
}

//==============================================================================
/// Defines mk for tracking (moving or floating bodies).
//==============================================================================
void JCaseVResBase::TrackingConfig(int mkbound,int mk){
  if(mk<=0 || mkbound<0)Run_Exceptioon("Tracking configuration (mk or mkbound) is invalid.");
  TrackingMk=mk;
  TrackingMkBound=mkbound;
}

//==============================================================================
// Returns the requested zone.
//==============================================================================
const JCaseVResBase* JCaseVResBase::GetSubZone(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("The requested subdomain does not exist.");
  return((const JCaseVResBase*)SubZones[idx]);
}

//==============================================================================
// Returns the requested zone as JCaseVRes_Box object.
//==============================================================================
const JCaseVRes_Box* JCaseVResBase::GetSubZoneBox(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("The requested subdomain does not exist.");
  if(SubZones[idx]->Shape!=MRSH_Box)Run_Exceptioon("The requested subdomain is not a box domain.");
  return((const JCaseVRes_Box*)SubZones[idx]);
}

//==============================================================================
// Returns string with list of subzones.
//==============================================================================
std::string JCaseVResBase::GetSubZonesStr()const{
  string tx;
  const unsigned n=Count();
  for(unsigned c=0;c<n;c++){
    const JCaseVResBase* zo=GetSubZone(c);
    if(c)tx=tx+",";
    tx=tx+fun::IntStr(zo->Id);
  }
  return(tx);
}


//##############################################################################
//# JCaseVRes_Box
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCaseVRes_Box::JCaseVRes_Box(int id,JCaseVResBase* parent,double hdp,double dp
  ,double bsizeh,double overlaph,std::string filerow,tdouble3 configref
  ,const JBoxDef& configbox)
  :JCaseVResBase(MRSH_Box,id,parent,configbox.GetIs2D(),configbox.GetPosy2D()
  ,hdp,dp,bsizeh,overlaph,filerow)
{
  Reset();
  ConfigRef=configref;
  ConfigBox=configbox;
  if(id==0 && !ConfigBox.IsSimple())
    Run_ExceptioonFile("Position transformation is invalid for main domain.",FileRow);

  if(1){
    if(parent)printf("  Zone_Id:%d  (ParentId:%d)\n",id,parent->Id);
    else      printf("  Zone_Id:%d\n",id);
    printf("    Dp:%f\n",dp);
    string txpref=(ConfigRef.x!=DBL_MAX? fun::DoubleStr(ConfigRef.x): "undefined");
    txpref=txpref+","+(ConfigRef.y!=DBL_MAX ? fun::DoubleStr(ConfigRef.y) : "undefined");
    txpref=txpref+","+(ConfigRef.z!=DBL_MAX ? fun::DoubleStr(ConfigRef.z) : "undefined");
    printf("    ConfigRef:(%s)\n",txpref.c_str());
    printf("    ConfigMin:%s\n",fun::Double3Str(ConfigBox.GetPosMin()).c_str());
    printf("    ConfigMax:%s\n",fun::Double3Str(ConfigBox.GetPosMax()).c_str());
    if(!ConfigBox.IsSimple()){
      printf("    ConfigVx:%s\n",fun::Double3Str(ConfigBox.GetVx()).c_str());
      printf("    ConfigVy:%s\n",fun::Double3Str(ConfigBox.GetVy()).c_str());
      printf("    ConfigVz:%s\n",fun::Double3Str(ConfigBox.GetVz()).c_str());
    }
  }

  PtBox=ConfigBox;
  if(parent){
    if(Parent->Shape!=MRSH_Box)Run_Exceptioon("Only Box shape is not supported.");
    const JCaseVRes_Box* pazone=(const JCaseVRes_Box*)Parent;
    //-Fit domain according to the own ConfigRef by increasing the domain.
    PtBox=FitDomainBox(ConfigBox);
    string err1=CheckSubOutDomainBox("The subdomain limits %s are not within parent domain %s."
      ,PtBox,pazone->GetPtBox());
    if(!err1.empty())Run_ExceptioonFile(
      fun::PrintStr("Limits of the inner subdomain Id=%d are invalid. ",Id)+err1,FileRow);
    
    //-Compute buffer limits around the subdomain.
    const double h=dp*Hdp;
    BuffBox=ComputeExtendingLimitsBox(PtBox,Dp,h*BuffSizeh);
    //-Compute fixed limits around the buffer limits.
    FixedBox=ComputeExtendingLimitsBox(BuffBox,Dp,h*2.0);

    //-Check fixed limits within parent domain limits.
    string err2=CheckSubOutDomainBox("The fixed limits of subdomain %s are not within parent domain %s."
      ,FixedBox,pazone->GetPtBox());
    if(!err2.empty())Run_ExceptioonFile(
      fun::PrintStr("Fixed limits of the subdomain Id=%d are invalid. ",Id)+err2,FileRow);
    
    //-Compute parent buffer limits within the subdomain.
    ComputeParentBufferLimitsBox(ParentBuffIniBox,ParentBuffEndBox,ParentFixedBox);

    //-Adds subdomain to parent domain.
    Parent->SubZones.push_back(this);
  }
}

//==============================================================================
/// Constructor.
//==============================================================================
JCaseVRes_Box::JCaseVRes_Box(int id,JCaseVResBase* parent,bool is2d
  ,double posy2d,double dp,const JXml* sxml,TiXmlElement* item)
  :JCaseVResBase(MRSH_Box,id,parent,is2d,posy2d,0,dp,0,0
  ,sxml->ErrGetFileRow(item))
{
  Reset();
  ReadXmlRun(sxml,item);
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JCaseVRes_Box::Reset(){
  TrackingDisable();
  SubZones.clear();
  ConfigRef=TDouble3(0);
  ConfigBox.Reset();
  PtBox.Reset();
  BuffBox.Reset();
  FixedBox.Reset();
  ParentBuffIniBox.Reset();
  ParentBuffEndBox.Reset();
  ParentFixedBox.Reset();
}

//==============================================================================
// Return true when all boxes definition are simple.
//==============================================================================
bool JCaseVRes_Box::UseSimpleBoxes()const{
  bool simple=ConfigBox.IsSimple();
  simple=simple && PtBox.IsSimple();
  simple=simple && BuffBox.IsSimple();
  simple=simple && FixedBox.IsSimple();
  simple=simple && ParentBuffIniBox.IsSimple();
  simple=simple && ParentBuffEndBox.IsSimple();
  simple=simple && ParentFixedBox.IsSimple();
  for(unsigned c=0;c<Count();c++){
    simple=simple && (GetSubZone(c)->Shape==MRSH_Box && GetSubZoneBox(c)->UseSimpleBoxes());
  }
  return(simple);
}

//==============================================================================
// Checks one subdomain is within other domain and return error details.
// E.g.: tex="The subdomain %s is not within domain %s."
//==============================================================================
std::string JCaseVRes_Box::CheckSubOutDomainBox(const std::string& tex
  ,const JBoxDef& subdomain,const JBoxDef& domain)
{
  string txerr;
  const bool inside=domain.Inside(subdomain);
  if(!inside){
    txerr=fun::PrintStr(tex.c_str(),subdomain.GetDomainStr().c_str()
      ,domain.GetDomainStr().c_str());
  }
  return(txerr);
}

//==============================================================================
// Write JBoxDef data in XML object.
//==============================================================================
void JCaseVRes_Box::WriteXmlRunBox(JXml* sxml,TiXmlElement* ele
  ,JBoxDef box,double resize)const
{
  if(resize)box.Resize(resize);
  if(box.IsSimple()){
    sxml->AddElementDouble3(ele,"pointmin",box.GetPosMin());
    sxml->AddElementDouble3(ele,"pointmax",box.GetPosMax());
  }
  else{
    sxml->AddElementDouble3(ele,"pointmin",box.GetPosMin());
    sxml->AddElementDouble3(ele,"vectorx"    ,box.GetVx());
    if(!box.GetIs2D())sxml->AddElementDouble3(ele,"vectory",box.GetVy());
    sxml->AddElementDouble3(ele,"vectorz"    ,box.GetVz());
  }
}

//==============================================================================
// Read JBoxDef data in XML object.
//==============================================================================
JBoxDef JCaseVRes_Box::ReadXmlRunBox(const JXml* sxml,const TiXmlElement* ele
  ,bool is2d)const
{
  JBoxDef box;
  const tdouble3 pmin=sxml->ReadElementDouble3(ele,"pointmin");
  if(!sxml->ExistsElement(ele,"vectorx")){
    const tdouble3 pmax=sxml->ReadElementDouble3(ele,"pointmax");
    box=JBoxDef(pmin,pmax,is2d);
  }
  else{
    const tdouble3 vx=sxml->ReadElementDouble3(ele,"vectorx");
    const tdouble3 vy=(!is2d? sxml->ReadElementDouble3(ele,"vectory"): TDouble3(0));
    const tdouble3 vz=sxml->ReadElementDouble3(ele,"vectorz");
    box=JBoxDef(pmin,vx,vy,vz,is2d);
  }
  return(box);
}

//==============================================================================
// Write execution configuration in XML object.
//==============================================================================
void JCaseVRes_Box::WriteXmlRun(JXml* sxml,TiXmlElement* ele,bool halfdp)const{
  TiXmlElement item("bufferbox");
  //-Write basic configuration.
  sxml->AddAttribute(&item,"id",Id);
  if(Parent)sxml->AddAttribute(&item,"parentid",Parent->Id);
  sxml->AddAttribute(&item,"dp",Dp);
  if(Parent)sxml->AddAttribute(&item,"dpratio",Parent->Dp/Dp);
  //-Write tracking configuration.
  if(TrackingIsActive()){
    TiXmlElement item2("tracking");
    sxml->AddAttribute(&item2,"mkbound",GetTrackingMkBound());
    sxml->AddAttribute(&item2,"mk",GetTrackingMk());
    item.InsertEndChild(item2);
  }
  //-Write full domain (FixedMin/Max).
  const double dpm=(halfdp? Dp/2: 0.);
  {
    TiXmlElement item2("fulldomain");
    sxml->AddElementDouble3(&item2,"pointref",GetConfigRef());
    if(Parent)WriteXmlRunBox(sxml,&item2,FixedBox,dpm);
    else      WriteXmlRunBox(sxml,&item2,PtBox   ,0);
    item.InsertEndChild(item2);
  }
  //-Write buffer domain (BuffMin/Max).
  if(Parent){
    TiXmlElement item2("bufferdomain");
    WriteXmlRunBox(sxml,&item2,BuffBox,dpm);
    item.InsertEndChild(item2);
  }
  //-Write inner domain (PtMin/Max).
  if(Parent){
    TiXmlElement item2("innerdomain");
    WriteXmlRunBox(sxml,&item2,PtBox,dpm);
    item.InsertEndChild(item2);
  }

  //-Write parent buffer and fixed limits.
  if(Parent){
    const tdouble3 padpm3=(halfdp? TDouble3(Parent->Dp/2,(Is2D? 0: Parent->Dp/2),Parent->Dp/2): TDouble3(0));
    const double padpm=(halfdp? -Parent->Dp/2: 0);
    {
      TiXmlElement item2("parentbufferinilimits");
      WriteXmlRunBox(sxml,&item2,ParentBuffIniBox,padpm);
      item.InsertEndChild(item2);
    }
    {
      TiXmlElement item2("parentbufferendlimits");
      WriteXmlRunBox(sxml,&item2,ParentBuffEndBox,padpm);
      item.InsertEndChild(item2);
    }
    {
      TiXmlElement item2("parentfixedlimits");
      WriteXmlRunBox(sxml,&item2,ParentFixedBox,padpm);
      item.InsertEndChild(item2);
    }
  }
  //-Add new item.
  TiXmlElement* subitem=ele->InsertEndChild(item)->ToElement();
  //-Write extra data section.
  WriteXmlExtra(sxml,subitem);
  //-Write subdomains of current subdomain.
  for(unsigned c=0;c<Count();c++){
    if(GetSubZone(c)->Shape==MRSH_Box)GetSubZoneBox(c)->WriteXmlRun(sxml,subitem,halfdp);
    else Run_Exceptioon("Only Box shape is supported.");
  }
}

//==============================================================================
// Read execution configuration from XML object.
//==============================================================================
void JCaseVRes_Box::ReadXmlRun(const JXml* sxml,TiXmlElement* item){
  const string ckvalid=string("*bufferbox tracking fulldomain bufferdomain")
    +" innerdomain parentbufferinilimits parentbufferendlimits parentfixedlimits"
    +" extra";
  sxml->CheckElementNames(item,true,ckvalid);
  //-Read tracking configuration.
  if(sxml->ExistsElement(item,"tracking")){
    sxml->CheckAttributeNames(item,"tracking","mkbound mk comment");
    const int mkbound=sxml->ReadElementInt(item,"tracking","mkbound");
    const int mk     =sxml->ReadElementInt(item,"tracking","mk");
    TrackingConfig(mkbound,mk);
  }
  //-Read full domain (FixedMin/Max).
  {
    TiXmlElement* item2=sxml->GetFirstElement(item,"fulldomain");
    ConfigRef=sxml->ReadElementDouble3(item2,"pointref");
    if(!Parent)PtBox   =ReadXmlRunBox(sxml,item2,Is2D);
    else       FixedBox=ReadXmlRunBox(sxml,item2,Is2D);
  }
  //-Read buffer domain (BuffMin/Max).
  if(Parent){
    TiXmlElement* item2=sxml->GetFirstElement(item,"bufferdomain");
    BuffBox=ReadXmlRunBox(sxml,item2,Is2D);
  }
  //-Read inner domain (PtMin/Max).
  if(Parent){
    TiXmlElement* item2=sxml->GetFirstElement(item,"innerdomain");
    PtBox=ReadXmlRunBox(sxml,item2,Is2D);
  }
  //-Read parent buffer limits.
  if(Parent){
    TiXmlElement* item2=sxml->GetFirstElement(item,"parentbufferinilimits");
    ParentBuffIniBox=ReadXmlRunBox(sxml,item2,Is2D);
  }
  if(Parent){
    TiXmlElement* item2=sxml->GetFirstElement(item,"parentbufferendlimits");
    ParentBuffEndBox=ReadXmlRunBox(sxml,item2,Is2D);
  }
  if(Parent){
    TiXmlElement* item2=sxml->GetFirstElement(item,"parentfixedlimits");
    ParentFixedBox=ReadXmlRunBox(sxml,item2,Is2D);
  }
}

//==============================================================================
/// Compute domain limits for particles (PtMin,PtMax).
/// Fit domain according to the own ConfigRef by increasing the domain.
//==============================================================================
void JCaseVRes_Box::FitDomain(tdouble3 configmin,tdouble3 configmax
  ,tdouble3& ptmin,tdouble3& ptmax)const
{
  const double dmin=Dp*0.001;
  //-Fit minimum limit (rounded down).
  ptmin.x=fmath::CalcRoundPos(configmin.x,ConfigRef.x,Dp);
  ptmin.y=fmath::CalcRoundPos(configmin.y,ConfigRef.y,Dp);
  ptmin.z=fmath::CalcRoundPos(configmin.z,ConfigRef.z,Dp);
  if(ptmin.x>(configmin.x+dmin))ptmin.x-=Dp;
  if(ptmin.y>(configmin.y+dmin))ptmin.y-=Dp;
  if(ptmin.z>(configmin.z+dmin))ptmin.z-=Dp;
  //-Fit maximum limit (rounded up).
  ptmax.x=fmath::CalcRoundPos(configmax.x,ConfigRef.x,Dp);
  ptmax.y=fmath::CalcRoundPos(configmax.y,ConfigRef.y,Dp);
  ptmax.z=fmath::CalcRoundPos(configmax.z,ConfigRef.z,Dp);
  if(ptmax.x<(configmax.x-dmin))ptmax.x+=Dp;
  if(ptmax.y<(configmax.y-dmin))ptmax.y+=Dp;
  if(ptmax.z<(configmax.z-dmin))ptmax.z+=Dp;
  if(Is2D)ptmin.y=ptmax.y=Posy2D;
  //printf("FitDomain> Parent_z%d  dp:%f  pmin:%s  pmax:%s\n",Id,Dp
  //  ,fun::Double3Str(PtMin).c_str(),fun::Double3Str(PtMax).c_str());
  //printf("FitDomain> z%d  dp:%f  pmin:%s  pmax:%s\n",pzone->Id,pzone->Dp
  //  ,fun::Double3Str(ptmin).c_str(),fun::Double3Str(ptmax).c_str());
  //printf("FitDomain> Result z%d  dp:%f  pmin:%s  pmax:%s\n",Id,Dp
  //  ,fun::Double3Str(ptmin).c_str(),fun::Double3Str(ptmax).c_str());
}

//==============================================================================
/// Fit vector size according to dp units.
//==============================================================================
tdouble3 JCaseVRes_Box::FitVectorDp(const tdouble3& v,double dp)const{
  double size=fun::Length(v);
  double size2=dp*int(round(size/dp));
  return(size>0? v*(size2/size): TDouble3(0));
}

//==============================================================================
/// Compute domain limits for particles (PtBox).
/// If simple definition is used then fit domain according to the own ConfigRef 
/// by increasing the domain.
//==============================================================================
JBoxDef JCaseVRes_Box::FitDomainBox(const JBoxDef& configbox)const{
  const double dmin=Dp*0.001;
  JBoxDef ret;
  //-Fit minimum limit (rounded down).
  if(configbox.IsSimple()){
    const tdouble3 configmin=configbox.GetPosMin();
    const tdouble3 configmax=configbox.GetPosMax();
    tdouble3 ptmin=configmin;
    tdouble3 ptmax=configmax;
    FitDomain(configmin,configmax,ptmin,ptmax);
    ret=JBoxDef(ptmin,ptmax,configbox.GetIs2D());
  }
  else{
    const tdouble3 p0=configbox.GetPosMin();
    const tdouble3 vx=FitVectorDp(configbox.GetVx(),Dp);
    const tdouble3 vy=FitVectorDp(configbox.GetVy(),Dp);
    const tdouble3 vz=FitVectorDp(configbox.GetVz(),Dp);
    ret=JBoxDef(p0,vx,vy,vz,configbox.GetIs2D());
  }
  return(ret);
}

//==============================================================================
/// Compute extending limits around the box domain
/// (pmin0 - resizemin, pmax0 + resizemin).
//==============================================================================
void JCaseVRes_Box::ComputeExtendingLimits(const tdouble3& pmin0
  ,const tdouble3& pmax0,double dp,double resizemin,tdouble3& pmin,tdouble3& pmax)const
{
  const double resize=dp*int((resizemin+dp*0.999)/dp);
  pmin=pmin0-TDouble3(resize);
  pmax=pmax0+TDouble3(resize);
  if(Is2D)pmin.y=pmax.y=Posy2D;
}

//==============================================================================
/// Compute extending limits around the box domain
/// (min - resizemin, max + resizemin).
//==============================================================================
JBoxDef JCaseVRes_Box::ComputeExtendingLimitsBox(const JBoxDef& box,double dp
  ,double resizemin)const
{
  const double resize=dp*int((resizemin+dp*0.999)/dp);
  JBoxDef box2=box;
  box2.Resize(resize);
  return(box2);
}

//==============================================================================
/// Compute reduction limits around the box domain 
/// (pmin0 + resizemin, pmax0 - resizemin).
//==============================================================================
void JCaseVRes_Box::ComputeReductionLimits(const tdouble3& pmin0
  ,const tdouble3& pmax0,double dp,double resizemin,tdouble3& pmin,tdouble3& pmax)const
{
  const double resize=dp*int((resizemin+dp*0.999)/dp);
  pmin=pmin0+TDouble3(resize);
  pmax=pmax0-TDouble3(resize);
  if(Is2D)pmin.y=pmax.y=Posy2D;
}

//==============================================================================
/// Compute parent fixed and buffer limits within the subdomain.
//==============================================================================
void JCaseVRes_Box::ComputeParentBufferLimits(const JBoxDef& ptbox
  ,tdouble3& pinimin,tdouble3& pinimax
  ,tdouble3& pendmin,tdouble3& pendmax
  ,tdouble3& pfixmin,tdouble3& pfixmax)const
{
  if(Parent->Shape!=MRSH_Box)Run_Exceptioon("Only Box shape is not supported.");
  const JCaseVRes_Box* pabox=(const JCaseVRes_Box*)Parent;
  const double dp=pabox->Dp;
  const double dmin=dp*0.001;
  const tdouble3 cref=pabox->ConfigRef;

  //-Compute inner limits of parent buffer within the subdomain.
  //-Fit minimum limit (rounded up).
  const double soverlap=Overlaph*(Dp*Hdp);
  const tdouble3 soverlap3=TDouble3(soverlap,(Is2D? 0: soverlap),soverlap);
  tdouble3 pmin0=ptbox.GetPosMin() + soverlap3;
  pinimin.x=fmath::CalcRoundPos(pmin0.x,cref.x,dp);
  pinimin.y=fmath::CalcRoundPos(pmin0.y,cref.y,dp);
  pinimin.z=fmath::CalcRoundPos(pmin0.z,cref.z,dp);
  if(pinimin.x<(pmin0.x-dmin))pinimin.x+=dp;
  if(pinimin.y<(pmin0.y-dmin))pinimin.y+=dp;
  if(pinimin.z<(pmin0.z-dmin))pinimin.z+=dp;
  //-Fit maximum limit (rounded down).
  tdouble3 pmax0=ptbox.GetPosMax() - soverlap3;
  pinimax.x=fmath::CalcRoundPos(pmax0.x,cref.x,dp);
  pinimax.y=fmath::CalcRoundPos(pmax0.y,cref.y,dp);
  pinimax.z=fmath::CalcRoundPos(pmax0.z,cref.z,dp);
  if(pinimax.x>(pmax0.x+dmin))pinimax.x-=dp;
  if(pinimax.y>(pmax0.y+dmin))pinimax.y-=dp;
  if(pinimax.z>(pmax0.z+dmin))pinimax.z-=dp;
  if(Is2D)pinimin.y=pinimax.y=Posy2D;
  if(!(pinimin<=pinimax))Run_ExceptioonFile(
    "Parent inner buffer limits within the subdomain are invalid.",FileRow);

  //-Compute outer limits of parent buffer within the subdomain.
  const double h=dp*pabox->Hdp;
  ComputeReductionLimits(pinimin,pinimax,dp,h*pabox->BuffSizeh,pendmin,pendmax);
  if(!(pendmin<=pendmax))Run_ExceptioonFile(fun::PrintStr(
    "Parent outer buffer limits within the subdomain Id=%d are invalid.",Id),FileRow);

  //-Compute limits of parent fixed zone within the subdomain.
  ComputeReductionLimits(pendmin,pendmax,dp,h*2.0,pfixmin,pfixmax);
  //printf("==-> pfixmin-max:%s\n",fun::Double3gRangeStr(pfixmin,pfixmax).c_str());
  if(!(pfixmin<=pfixmax))Run_ExceptioonFile(fun::PrintStr(
    "Parent fixed limits within the subdomain Id=%d are invalid.",Id),FileRow);
}

//==============================================================================
/// Compute parent fixed and buffer limits within the subdomain.
//==============================================================================
void JCaseVRes_Box::ComputeParentBufferLimitsBox(JBoxDef& pini,JBoxDef& pend
  ,JBoxDef& pfix)const
{
  if(Parent->Shape!=MRSH_Box)Run_Exceptioon("Only Box shape is not supported.");
  if(PtBox.IsSimple()){
    tdouble3 pinimin,pinimax,pendmin,pendmax,pfixmin,pfixmax;
    ComputeParentBufferLimits(PtBox,pinimin,pinimax,pendmin,pendmax,pfixmin,pfixmax);
    pini=JBoxDef(pinimin,pinimax,Is2D);
    pend=JBoxDef(pendmin,pendmax,Is2D);
    pfix=JBoxDef(pfixmin,pfixmax,Is2D);
  }
  else{
    const JCaseVRes_Box* pabox=(const JCaseVRes_Box*)Parent;
    const double dp=pabox->Dp;

    //-Compute inner limits of parent buffer within the subdomain.
    {
      const double soverlapmin=Overlaph*(Dp*Hdp);
      const double soverlap=dp*int((soverlapmin+dp*0.999)/dp);
      pini=PtBox;
      if(pini.GetMinSize()<=soverlap*2)Run_ExceptioonFile(
        "Parent inner buffer limits within the subdomain are invalid.",FileRow);
      pini.Resize(-soverlap);
    }

    //-Compute outer limits of parent buffer within the subdomain.
    {
      const double h=dp*pabox->Hdp;
      const double resizemin=h*pabox->BuffSizeh;
      const double resize=dp*int((resizemin+dp*0.999)/dp);
      pend=pini;
      if(pend.GetMinSize()<=resize*2)Run_ExceptioonFile(fun::PrintStr(
        "Parent outer buffer limits within the subdomain Id=%d are invalid.",Id),FileRow);
      pend.Resize(-resize);
    }

    //-Compute limits of parent fixed zone within the subdomain.
    {
      const double h=dp*pabox->Hdp;
      const double resizemin=h*2;
      const double resize=dp*int((resizemin+dp*0.999)/dp);
      pfix=pend;
      if(pfix.GetMinSize()<=resize*2)Run_ExceptioonFile(fun::PrintStr(
        "Parent fixed limits within the subdomain Id=%d are invalid.",Id),FileRow);
      pfix.Resize(-resize);
    }
  }
}

//##############################################################################
//# JCaseVRes
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JCaseVRes::JCaseVRes(){
  ClassName="JCaseVRes";
  Reset();
}

//==============================================================================
// Destructor.
//==============================================================================
JCaseVRes::~JCaseVRes(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JCaseVRes::Reset(){
  Is2D=false;
  Posy2D=0;
  HSwl=false;
  DepthDismin=DBL_MAX;
  RunData=false;
  SelId=UINT_MAX;
  #ifndef _WITHMR  //-It should be corrected in VRes code.
  for(unsigned c=0;c<Zones.size();c++)delete Zones[c];
  #endif
  Zones.clear();
}

//==============================================================================
/// Loads <pointref> definition in XML according to dp.
//==============================================================================
tdouble3 JCaseVRes::LoadXmlPtref(const JXml* sxml,double dp)const{
  tdouble3 ptref=TDouble3(0);
#ifdef JXml_UseNux
  #ifndef _WITHMR  //-VRes code must be updated to enable this code.
  JNumexLib* nuxlib=sxml->GetNuxLib();
  if(!nuxlib)Run_Exceptioon("Error: JNumexLib object is not available.");
  const unsigned dpidx=nuxlib->GetVarIdx("Dp");
  if(dpidx==UINT_MAX)Run_Exceptioon("Error: Dp variable is undefined.");
  const double dp0=nuxlib->GetVarNum(dpidx);
  nuxlib->CreateVar("Dp",true,true,dp);
  const TiXmlElement* ele=sxml->GetNodeError("case.casedef.geometry.definition")->ToElement();
  ptref.x=sxml->ReadElementDouble(ele,"pointref","x",true,ptref.x);
  ptref.y=sxml->ReadElementDouble(ele,"pointref","y",true,ptref.y);
  ptref.z=sxml->ReadElementDouble(ele,"pointref","z",true,ptref.z);
  //-Restore original dp value.
  nuxlib->CreateVar("Dp",true,true,dp0);
  #endif
#else
  Run_Exceptioon("JNumexLib is not available.");
#endif
  return(ptref);
}

//==============================================================================
/// Reads transformation configuration from XML (<transform>).
//==============================================================================
void JCaseVRes::ReadXmlTransform(const JXml* sxml,const TiXmlElement* ele
  ,bool is2d,tdouble3& mov,tdouble3& rot,tdouble3& rotcen)
{
  mov=rot=rotcen=TDouble3(0);
  const TiXmlElement* elet=sxml->GetFirstElement(ele,"transform",true);
  if(elet){
    sxml->CheckElementNames(elet,true,"move rotate center");
    sxml->CheckAttributeNames(elet,"rotate","angx angy angz anglesunits");
    mov=sxml->ReadElementDouble3(elet,"move",true);
    if(is2d && mov.y)sxml->ErrReadElement(elet,"move",false
      ,"Translation in Y is not allowed for 2-D definition."); 
    rot.x=sxml->ReadElementDouble(elet,"rotate","angx",true);
    rot.y=sxml->ReadElementDouble(elet,"rotate","angy",true);
    rot.z=sxml->ReadElementDouble(elet,"rotate","angz",true);
    if(is2d && (rot.x || rot.z))sxml->ErrReadElement(elet,"rotate",false
      ,"Rotation in X or Z is not allowed for 2-D definition."); 
    if(rot!=TDouble3(0)){
      rotcen=sxml->ReadElementDouble3(elet,"center",true);
      string angunits=fun::StrLower(sxml->ReadElementStr(elet,"rotate","anglesunits"));
      if(angunits=="radians")rot=rot*TODEG;
      else if(angunits!="degrees")sxml->ErrReadElement(elet,"rotate",false
        ,"The value anglesunits must be \"degrees\" or \"radians\"."); 
    }
  }
}

//==============================================================================
/// Reads variable resolution configuration from XML node for GenCase.
//==============================================================================
void JCaseVRes::ReadXmlDef(JCaseVResBase* parent,const JXml* sxml
  ,TiXmlElement* lis,int mkboundfirst,double hdp)
{
  if(parent->Id==0)sxml->CheckElementNames(lis,true
    ,"*bufferbox buffsizeh simulationdomain extra");
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    std::string cmd=ele->Value();
    //printf("\n cmd:[%s]  find:%d  from %s\n",cmd.c_str(),int(cmd.find("buffer")),sxml->ErrGetFileRow(ele).c_str());
    if(cmd!="bufferbox")cmd="";
    if(!cmd.empty() && cmd[0]!='_' && sxml->CheckElementActive(ele)){
      //printf("\n cmd:[%s] from %s\n",cmd.c_str(),sxml->ErrGetFileRow(ele).c_str());
      //-Other configurations.
      if(cmd=="bufferbox"){
        const string ckvalid=string("*bufferbox dpratio buffsizeh overlappingh")
          +" point size endpoint tracking simulationdomain extra transform";
        sxml->CheckElementNames(ele,true,ckvalid);
        //-General values.
        const int id=Count();
        const int parentid=parent->Id;
        const double dpratio=sxml->ReadElementDouble(ele,"dpratio","v");
        const double dp=parent->Dp/dpratio;
        const double buffsizeh=sxml->ReadElementDouble(ele,"buffsizeh","v",true,2.0);
        const double overlappingh=sxml->ReadElementDouble(ele,"overlappingh","v",true,0);
        const string filerow=sxml->ErrGetFileRow(ele);
        const bool is2d=Zones[0]->Is2D;
        const double posy2d=Zones[0]->Posy2D;
        //-Loads tracking configuration.
        sxml->CheckAttributeNames(ele,"tracking","mkbound comment");
        const int trackingmkb=sxml->ReadElementInt(ele,"tracking","mkbound",true,-1);
        //-Loads pointref according to dp value.
        const tdouble3 ptref=LoadXmlPtref(sxml,dp);
        //-Loads initial limits.
        JBoxDef ptbox;
        {
          //-Loads ptmin and ptmax.
          tdouble3 ptmin=sxml->ReadElementDouble3(ele,"point");
          tdouble3 ptmax=ptmin;
          if(sxml->ExistsElement(ele,"size") && sxml->ExistsElement(ele,"endpoint"))
            Run_ExceptioonFile("Buffer zone dimensions must be defined with <size> or <endpoint>, but not both.",filerow);
          if(sxml->ExistsElement(ele,"endpoint"))ptmax=sxml->ReadElementDouble3(ele,"endpoint");
          else{
            const tdouble3 size=sxml->ReadElementDouble3(ele,"size");
            ptmax=ptmin+size;
          }
          if(is2d)ptmin.y=ptmax.y=posy2d;
          ptbox=JBoxDef(ptmin,ptmax,is2d);
          //-Loads transformation.
          tdouble3 mov,rot,rotcen;
          ReadXmlTransform(sxml,ele,is2d,mov,rot,rotcen);
          ptbox.MovRotate(mov,rot,rotcen);
        }
        //-Creates buffer zone object.
        JCaseVRes_Box* pzone=new JCaseVRes_Box(id,parent,hdp,dp,buffsizeh
          ,overlappingh,filerow,ptref,ptbox);
        if(trackingmkb>=0)pzone->TrackingConfig(trackingmkb,mkboundfirst+trackingmkb);
        pzone->SimReadXmlDef(sxml,ele);
        //-Loads extra data from XML.
        JCaseVResExtra edata;
        ReadXmlExtra(sxml,sxml->GetFirstElement(ele,"extra",true),edata);
        pzone->SetExtraData(edata);
        //-Add new zone and loads subdomains from XML.
        Zones.push_back(pzone);
        ReadXmlDef(pzone,sxml,ele,mkboundfirst,hdp);
      }
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Loads extra section from XML object for GenCase and DualSPHysics.
//==============================================================================
void JCaseVRes::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele
  ,JCaseVResExtra& edata,std::string kparent)
{
  if(ele){
    const string elename=ele->Value();
    //-Loads attributes.
    if(elename!="extra"){
      kparent=kparent+(kparent.empty()? "": ".")+elename;
      const TiXmlAttribute* att=ele->FirstAttribute();
      while(att){
        const string attname=att->Name();
        if(!attname.empty() && attname[0]!='_'){
          const string k=kparent+(kparent.empty()? "": ".")+attname;
          if(edata.FindKey(k)!=UINT_MAX)
            Run_ExceptioonFile("The key is already stored.",sxml->ErrGetFileRow(ele));
          edata.AddValue(k,att->Value());
        }
        att=att->Next();
      }
    }
    //-Loads elements.
    const TiXmlElement* ele2=ele->FirstChildElement(); 
    while(ele2){
      const string ename=ele2->Value();
      if(!ename.empty() && ename[0]!='_')ReadXmlExtra(sxml,ele2,edata,kparent);
      ele2=ele2->NextSiblingElement();
    }
  }
}

//==============================================================================
/// Loads configuration from XML object for GenCase.
//==============================================================================
void JCaseVRes::LoadXmlDef(const JXml* sxml,const std::string& place
  ,int mkboundfirst,double hdp,double dp,tdouble3 ptref,tdouble3 ptmin
  ,tdouble3 ptmax)
{
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place,true);
  //if(!node)Run_Exceptioon(string("Cannot find the element \'")+place+"\'.");
  if(node){
    Is2D=(ptmin.y==ptmax.y);
    Posy2D=(Is2D? ptmin.y: 0);
    if(Is2D)ptref.y=0;
    //-Creates domain 0 (full domain).
    TiXmlElement* ele=node->ToElement();
    const double buffsizeh=sxml->ReadElementDouble(ele,"buffsizeh","v",true,2.0);
    const double overlaph=0;
    JBoxDef boxpt(ptmin,ptmax,Is2D);
    JCaseVRes_Box* pzone=new JCaseVRes_Box(0,NULL,hdp,dp,buffsizeh
      ,overlaph,sxml->ErrGetFileRow(ele),ptref,boxpt);
    pzone->SimReadXmlDef(sxml,ele);
    //-Loads extra data from XML.
    JCaseVResExtra edata;
    ReadXmlExtra(sxml,sxml->GetFirstElement(ele,"extra",true),edata);
    pzone->SetExtraData(edata);
    //-Add new zone and loads subdomains from XML.
    Zones.push_back(pzone);
    ReadXmlDef(pzone,sxml,ele,mkboundfirst,hdp);
  }
}

//==============================================================================
// Returns the requested zone.
//==============================================================================
const JCaseVRes_Box* JCaseVRes::GetZoneBox(unsigned id)const{
  if(id>=Count())Run_Exceptioon("The requested subdomain does not exist.");
  if(Zones[id]->Shape!=MRSH_Box)Run_Exceptioon("The requested subdomain is not a box domain.");
  return((const JCaseVRes_Box*)Zones[id]);
}

//==============================================================================
/// Load vertices of 2-D quad from JBoxDef object with resize.
//==============================================================================
void JCaseVRes::GetBoxPoints2d(JBoxDef box,double resize,tdouble3* vpt){
  box.Resize(resize);
  tdouble3 points[8];
  box.GetBoxPoints(points);
  vpt[0]=points[0];
  vpt[1]=points[1];
  vpt[2]=points[5];
  vpt[3]=points[4];
}

//==============================================================================
/// Load vertice and vectors of 3-D box from JBoxDef object with resize.
//==============================================================================
void JCaseVRes::GetBoxPoints3d(JBoxDef box,double resize,tdouble3* ptvec){
  box.Resize(resize);
  ptvec[0]=box.GetPosMin();
  ptvec[1]=box.GetVx();
  ptvec[2]=box.GetVy();
  ptvec[3]=box.GetVz();
}

//==============================================================================
// Saves VTK file with domain of zones.
//==============================================================================
void JCaseVRes::SaveVtkDomains(std::string fname,bool onefile,bool halfdp)const{
  string file=fun::GetWithoutExtension(fname);
  if(Count()<2)Run_ExceptioonFile("No subdomains are available.",file+".vtk");
  const unsigned nz=Count();
  JSpVtkShape ss_1;
  if(Is2D){
    for(unsigned id=0;id<nz;id++){
      const word wid=word(id);
      JSpVtkShape ss_2;
      JSpVtkShape& ss=(onefile? ss_1: ss_2);
      const JCaseVRes_Box* pzone=GetZoneBox(id);
      const double fdpm=(halfdp && id? pzone->Dp/2: 0);
      const double dpm=(halfdp? pzone->Dp/2: 0);
      tdouble3 pt[4];
      GetBoxPoints2d(pzone->GetPtBox(),fdpm,pt);
      ss.AddQuadWire(pt[0],pt[1],pt[2],pt[3],wid);
      if(id){
        tdouble3 pt[4];
        GetBoxPoints2d(pzone->GetBuffBox(),dpm,pt);
        ss.AddQuadWire(pt[0],pt[1],pt[2],pt[3],wid);
        //GetBoxPoints2d(pzone->GetFixedBox(),dpm,pt);
        //ss.AddQuadWire(pt[0],pt[1],pt[2],pt[3],wid);
      }
      const unsigned n2=pzone->Count();
      for(unsigned c2=0;c2<n2;c2++){
        const JCaseVRes_Box* psubzone=pzone->GetSubZoneBox(c2);
        tdouble3 pt[4];
        GetBoxPoints2d(psubzone->GetParentBuffIniBox(),-dpm,pt);
        tdouble3 pb[4];
        GetBoxPoints2d(psubzone->GetParentBuffEndBox(),-dpm,pb);
        //tdouble3 pf[4];
        //GetBoxPoints2d(psubzone->GetParentFixedBox(),-dpm,pf);
        ss.AddQuadWire(pt[0],pt[1],pt[2],pt[3],wid);
        ss.AddQuadWire(pb[0],pb[1],pb[2],pb[3],wid);
        //sh.AddShapeQuadWire(pf[0],pf[1],pf[2],pf[3],id);
        for(unsigned c=0;c<4;c++){
          ss.AddLine(pt[c],pb[c],wid);
          //ss.AddLine(pb[c],pf[c],wid);
        }
      }
      if(!onefile)ss_2.SaveVtk(file+fun::PrintStr("%02d.vtk",id),"Zone");
    }
    if(onefile)ss_1.SaveVtk(file+".vtk","Zone");
  }
  else{
    for(unsigned id=0;id<nz;id++){
      const word wid=word(id);
      JSpVtkShape ss_2;
      JSpVtkShape& ss=(onefile? ss_1: ss_2);
      const JCaseVRes_Box* pzone=GetZoneBox(id);
      const double fdpm=(halfdp && id? pzone->Dp/2: 0);
      const double dpm=(halfdp? pzone->Dp/2: 0);

      tdouble3 ptvec[4];
      GetBoxPoints3d(pzone->GetPtBox(),fdpm,ptvec);
      ss.AddBoxSizeVec(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
      if(id){
        tdouble3 ptvec[4];
        GetBoxPoints3d(pzone->GetBuffBox(),dpm,ptvec);
        ss.AddBoxSizeVec(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
        //GetBoxPoints3d(pzone->GetFixedBox(),dpm,ptvec);
        //ss.AddBoxSizeVec(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
      }
      const unsigned n2=pzone->Count();
      for(unsigned c2=0;c2<n2;c2++){
        const JCaseVRes_Box* psubzone=pzone->GetSubZoneBox(c2);
        tdouble3 ptvec[4];
        GetBoxPoints3d(psubzone->GetParentBuffIniBox(),-dpm,ptvec);
        ss.AddBoxSizeVec(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
        GetBoxPoints3d(psubzone->GetParentBuffEndBox(),-dpm,ptvec);
        ss.AddBoxSizeVec(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
        //GetBoxPoints3d(psubzone->GetParentFixedBox(),-dpm,ptvec);
        //ss.AddBoxSizeVec(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
      }
      if(!onefile)ss_2.SaveVtk(file+fun::PrintStr("%02d.vtk",id),"Zone");
    }
    if(onefile)ss_1.SaveVtk(file+".vtk","Zone");
  }
}

//==============================================================================
// Saves VTK file with outer limits of zones.
//==============================================================================
void JCaseVRes::SaveVtkLimits(std::string fname,bool halfdp)const{
  string file=fun::GetWithoutExtension(fname);
  if(Count()<2)Run_ExceptioonFile("No subdomains are available.",file+".vtk");
  const unsigned nz=Count();
  JSpVtkShape ss;
  if(Is2D){
    for(unsigned id=0;id<nz;id++){
      const JCaseVRes_Box* pzone=GetZoneBox(id);
      const double fdpm=(halfdp && id? pzone->Dp/2: 0);
      const double dpm=(halfdp? pzone->Dp/2: 0);
      tdouble3 pt[4];
      GetBoxPoints2d(pzone->GetPtBox(),fdpm,pt);
      ss.AddQuadWire(pt[0],pt[1],pt[2],pt[3],word(id));
      const unsigned n2=pzone->Count();
    }
    ss.SaveVtk(file+".vtk","Zone");
  }
  else{
    for(unsigned id=0;id<nz;id++){
      const JCaseVRes_Box* pzone=GetZoneBox(id);
      const double fdpm=(halfdp && id? pzone->Dp/2: 0);
      const double dpm=(halfdp? pzone->Dp/2: 0);
      tdouble3 ptvec[4];
      GetBoxPoints3d(pzone->GetPtBox(),fdpm,ptvec);
      ss.AddBoxSizeVecWire(ptvec[0],ptvec[1],ptvec[2],ptvec[3],word(id));
    }
    ss.SaveVtk(file+".vtk","Zone");
  }
}

//==============================================================================
/// Stores configuration in XML object for DualSPHysics.
//==============================================================================
void JCaseVRes::SaveXmlRun(JXml *sxml,const std::string &place
  ,unsigned vresid,bool halfdp)const
{
  if(vresid>=Count())Run_Exceptioon("The requested subdomain does not exist.");
  TiXmlElement* mainitem=sxml->GetNode(place,true)->ToElement();
  sxml->AddAttribute(mainitem,"selid",vresid);
  if(Is2D)sxml->AddAttribute(mainitem,"posy2d",Posy2D);
  if(!Count())Run_Exceptioon("Error: No subdomains available.");
  vresid=0;
  if(Zones[vresid]->Shape==MRSH_Box)GetZoneBox(vresid)->WriteXmlRun(sxml,mainitem,halfdp);
  else Run_ExceptioonFile("Only Box shape is not supported.",sxml->ErrGetFileRow(mainitem));
}

//==============================================================================
/// Stores simulation domain configuration in XML object for DualSPHysics.
//==============================================================================
void JCaseVRes::SaveXmlSimRun(JXml* sxml,const std::string& place
  ,unsigned vresid)const
{
  if(vresid>=Count())Run_Exceptioon("The requested subdomain does not exist.");
  TiXmlElement* mainitem=sxml->GetNode(place,true)->ToElement();
  if(!Count())Run_Exceptioon("Error: No subdomains available.");
  if(Zones[vresid]->Shape==MRSH_Box)GetZoneBox(vresid)->SimWriteXmlRun(sxml,mainitem);
  else Run_ExceptioonFile("Only Box shape is not supported.",sxml->ErrGetFileRow(mainitem));
}

//==============================================================================
/// Set number of particles of zone.
//==============================================================================
void JCaseVRes::NpSet(unsigned id,ullong npfixed,ullong npmoving
  ,ullong npfloating,ullong npfluid)
{
  if(id>=Count())Run_Exceptioon("The requested subdomain does not exist.");
  Zones[id]->NpSet(npfixed,npmoving,npfloating,npfluid);
}

//==============================================================================
/// Reads variable resolution configuration from XML node for DualSPHysics.
/// **It should be improved when there are more configuration options.
//==============================================================================
void JCaseVRes::ReadXmlRun(JCaseVResBase* parent,const JXml* sxml
  ,TiXmlElement* lis)
{
  const string ckvalid=string("*bufferbox tracking fulldomain bufferdomain")
    +" innerdomain parentbufferinilimits parentbufferendlimits parentfixedlimits"
    +" extra";
  sxml->CheckElementNames(lis,true,ckvalid);
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    std::string cmd=ele->Value();
    if(cmd!="bufferbox")cmd="";
    if(cmd.length() && cmd[0]!='_' && sxml->CheckElementActive(ele)){
      if(cmd=="bufferbox"){
        //-Main attributes.
        const unsigned id=sxml->GetAttributeUnsigned(ele,"id");
        if(id!=Count())Run_ExceptioonFile("Id is invalid.",sxml->ErrGetFileRow(ele));
        const unsigned parentid=sxml->GetAttributeUnsigned(ele,"parentid",true,UINT_MAX);
        if(parentid!=(parent? parent->Id: UINT_MAX))
          Run_ExceptioonFile("ParentId is invalid.",sxml->ErrGetFileRow(ele));
        const double dp=sxml->GetAttributeDouble(ele,"dp");
        //-Creates box domain.
        JCaseVRes_Box* pzone=new JCaseVRes_Box(id,parent,Is2D,Posy2D,dp,sxml,ele);
        if(parent)parent->SubZones.push_back(pzone);
        //-Loads extra data from XML.
        JCaseVResExtra edata;
        ReadXmlExtra(sxml,sxml->GetFirstElement(ele,"extra",true),edata);
        pzone->SetExtraData(edata);
        //-Add new zone and loads subdomains from XML.
        Zones.push_back(pzone);
        ReadXmlRun(pzone,sxml,ele);
      }
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Loads execution data from XML file and returns selid value.
//==============================================================================
unsigned JCaseVRes::LoadFileXmlRun(const std::string& file
  ,const std::string& place)
{
  JXml jxml;
  jxml.LoadFile(file);
  return(LoadXmlRun(&jxml,place));
}

//==============================================================================
/// Loads execution data from XML object and returns selid value.
//==============================================================================
unsigned JCaseVRes::LoadXmlRun(const JXml* sxml,const std::string& place){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place,true);
  if(node){
    TiXmlElement* mainitem=node->ToElement();
    SelId=sxml->GetAttributeUint(mainitem,"selid");
    if(sxml->ExistsAttribute(mainitem,"posy2d")){
      Is2D=true;
      Posy2D=sxml->GetAttributeDouble(mainitem,"posy2d");
    }
    sxml->CheckElementNames(mainitem,true,"bufferbox");
    ReadXmlRun(NULL,sxml,mainitem);
    RunData=true;
  }
  return(SelId);
}

