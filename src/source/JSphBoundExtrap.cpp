//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2016, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphBoundExtrap.cpp \brief Implements the class \ref JSphBoundExtrap.

#include "JSphBoundExtrap.h"
#include "JSphCpu.h"
#include "JSphMk.h"
#include "JXml.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JMatrix4.h"
#include "JLinearValue.h"
#include "JSaveCsv2.h"
#include "JRangeFilter.h"
#include "JFormatFiles2.h"

#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSphBoundExtrapZone
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphBoundExtrapZone::JSphBoundExtrapZone(JLog2 *log,unsigned idzone,word mkbound
 ,TpDirection autodir,tdouble3 limitpos,tdouble3 direction)
 :Log(log),IdZone(idzone),MkBound(mkbound)
{
  ClassName="JSphBoundExtrapZone";
  Reset();
  AutoDir=autodir;
  LimitPos=limitpos;
  Direction=direction;
  Plane=ToTFloat4(fmath::PlanePtVec(LimitPos,Direction));
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphBoundExtrapZone::~JSphBoundExtrapZone(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphBoundExtrapZone::Reset(){
  BoundCode=0;
  AutoDir=DIR_None;
  LimitPos=Direction=TDouble3(0);
  Plane=TFloat4(0);
}

//==============================================================================
/// Configures BoundCode.
//==============================================================================
void JSphBoundExtrapZone::ConfigBoundCode(typecode boundcode){
  if(BoundCode)RunException("ConfigBoundCode",fun::PrintStr("BoundCode was already configured for mkbound=%u.",MkBound));
  BoundCode=boundcode;
}

//==============================================================================
/// Configures BoundCode.
//==============================================================================
void JSphBoundExtrapZone::ConfigAutoLimit(double halfdp,tdouble3 pmin,tdouble3 pmax){
  const tdouble3 pmed=(pmin+pmax)/2.;
  if(AutoDir!=DIR_None){
    switch(AutoDir){
      case DIR_Top:
        Direction=TDouble3(0,0,1);
        LimitPos=TDouble3(pmed.x,pmed.y,pmax.z+halfdp);
      break;
      case DIR_Bottom:
        Direction=TDouble3(0,0,-1);
        LimitPos=TDouble3(pmed.x,pmed.y,pmin.z-halfdp);
      break;
      case DIR_Left:
        Direction=TDouble3(-1,0,0);
        LimitPos=TDouble3(pmin.x-halfdp,pmed.y,pmed.z);
      break;
      case DIR_Right:
        Direction=TDouble3(1,0,0);
        LimitPos=TDouble3(pmax.x+halfdp,pmed.y,pmed.z);
      break;
      case DIR_Front:
        Direction=TDouble3(0,-1,0);
        LimitPos=TDouble3(pmed.x,pmin.y-halfdp,pmed.z);
      break;
      case DIR_Back:
        Direction=TDouble3(0,1,0);
        LimitPos=TDouble3(pmed.x,pmax.y+halfdp,pmed.z);
      break;
    }
    Plane=ToTFloat4(fmath::PlanePtVec(LimitPos,Direction));
  }
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphBoundExtrapZone::GetConfig(std::vector<std::string> &lines)const{
  if(AutoDir!=DIR_None){
    lines.push_back(fun::PrintStr("Limit position: (%g,%g,%g) (automatic)",LimitPos.x,LimitPos.y,LimitPos.z));
    lines.push_back(fun::PrintStr("Fluid Direction: (%g,%g,%g) (automatic)",Direction.x,Direction.y,Direction.z));
  }
  else{
    lines.push_back(fun::PrintStr("Limit position: (%g,%g,%g)",LimitPos.x,LimitPos.y,LimitPos.z));
    lines.push_back(fun::PrintStr("Fluid Direction: (%g,%g,%g)",Direction.x,Direction.y,Direction.z));
  }
}


//##############################################################################
//# JSphBoundExtrap
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphBoundExtrap::JSphBoundExtrap(JLog2 *log,JXml *sxml,const std::string &place,const JSphMk *mkinfo)
  :Log(log)
{
  ClassName="JSphBoundExtrap";
  Reset();
  LoadXml(sxml,place);
  UpdateMkCode(mkinfo);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphBoundExtrap::~JSphBoundExtrap(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphBoundExtrap::Reset(){
  DetermLimit=0;
  for(int c=0;c<List.size();c++)delete List[c];
  List.clear();
}

//==============================================================================
/// Returns true if mkbound value is already configured.
//==============================================================================
bool JSphBoundExtrap::ExistMk(word mkbound)const{
  bool ret=false;
  for(unsigned c=0;c<List.size() && !ret;c++)ret=(List[c]->MkBound==mkbound);
  return(ret);
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSphBoundExtrap::LoadXml(JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of configurations in the XML node.
//==============================================================================
void JSphBoundExtrap::ReadXml(const JXml *sxml,TiXmlElement* lis){
  const char met[]="ReadXml";
  //-Loads value determlimit.
  if(sxml->CountElements(lis,"determlimit")>1)sxml->ErrReadElement(lis,"determlimit",false,"Several definitions for this value.");
  DetermLimit=sxml->ReadElementFloat(lis,"determlimit","value",true,1e+3f);
  //-Loads list of inputs.
  TiXmlElement* ele=lis->FirstChildElement("mkzone"); 
  while(ele){
    const word mkbound=sxml->GetAttributeWord(ele,"mkbound");
    tdouble3 limitpos=TDouble3(0);
    tdouble3 direction=TDouble3(0);
    JSphBoundExtrapZone::TpDirection autodir=JSphBoundExtrapZone::DIR_None;
    string autodirtx=fun::StrLower(sxml->ReadElementStr(ele,"autoconfig","direction",true));
    if(autodirtx.empty()){
      limitpos=sxml->ReadElementDouble3(ele,"limitpoint");
      direction=sxml->ReadElementDouble3(ele,"direction");
    }
    else{
      if     (autodirtx=="top"   )autodir=JSphBoundExtrapZone::DIR_Top;
      else if(autodirtx=="bottom")autodir=JSphBoundExtrapZone::DIR_Bottom;
      else if(autodirtx=="left"  )autodir=JSphBoundExtrapZone::DIR_Left;
      else if(autodirtx=="right" )autodir=JSphBoundExtrapZone::DIR_Right;
      else if(autodirtx=="front" )autodir=JSphBoundExtrapZone::DIR_Front;
      else if(autodirtx=="back"  )autodir=JSphBoundExtrapZone::DIR_Back;
      if(autodir==JSphBoundExtrapZone::DIR_None)sxml->ErrReadElement(ele,"autoconfig",false,"Direction label is invalid.");
    }
    if(ExistMk(mkbound))RunException(met,"An input already exists for the same mkbound.");
    JSphBoundExtrapZone *zo=new JSphBoundExtrapZone(Log,GetCount(),mkbound,autodir,limitpos,direction);
    List.push_back(zo);
    ele=ele->NextSiblingElement("mkzone");
  }
}

//==============================================================================
/// Updates BoundCode of each configuration.
//==============================================================================
void JSphBoundExtrap::UpdateMkCode(const JSphMk *mkinfo){
  const char met[]="UpdateMkCode";
  for(unsigned c=0;c<GetCount();c++){
    const word mkbound=List[c]->MkBound;
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(List[c]->MkBound);
    if(cmk<mkinfo->Size() && (CODE_IsFixed(mkinfo->Mkblock(cmk)->Code) || CODE_IsMoving(mkinfo->Mkblock(cmk)->Code))){
      List[c]->ConfigBoundCode(mkinfo->Mkblock(cmk)->Code);
    }
    else RunException(met,fun::PrintStr("MkBound value (%u) is not a Mk fixed boundary valid.",List[c]->MkBound));
  }
}

//==============================================================================
/// Run automatic configuration of LimitPos and Direction for each configuration
/// and saves VTK file with limit configuration.
//==============================================================================
void JSphBoundExtrap::RunAutoConfig(double dp,const JSphMk *mkinfo){
  const char met[]="RunAutoConfig";
  for(unsigned c=0;c<GetCount();c++){
    const word mkbound=List[c]->MkBound;
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(List[c]->MkBound);
    if(cmk<mkinfo->Size()){
      const tdouble3 pmin=mkinfo->Mkblock(cmk)->GetPosMin();
      const tdouble3 pmax=mkinfo->Mkblock(cmk)->GetPosMax();
      List[c]->ConfigAutoLimit(dp/2.,pmin,pmax);
    }
    else RunException(met,fun::PrintStr("MkBound value (%u) is not a Mk fixed boundary valid.",List[c]->MkBound));
  }
  SaveVtkConfig(dp);
}

//==============================================================================
/// Saves VTK file with LimitPos and Direction for each configuration.
//==============================================================================
void JSphBoundExtrap::SaveVtkConfig(double dp)const{
  std::vector<JFormatFiles2::StShapeData> shapes;
  for(unsigned c=0;c<GetCount();c++){
    const JSphBoundExtrapZone* zo=List[c];
    const int mkbound=zo->MkBound;
    const tdouble3 ps=zo->GetLimitPos();
    const tdouble3 ps2=ps+(fmath::VecUnitary(zo->GetDirection())*(dp*3));
    shapes.push_back(JFormatFiles2::DefineShape_Line(ps,ps2,mkbound,0)); //-Direction line.
    tdouble3 pt1=ps-TDouble3(dp/2.);
    tdouble3 pt2=ps+TDouble3(dp/2.);
    const double dp2=dp*2;
    switch(zo->GetAutoDir()){
      case JSphBoundExtrapZone::DIR_Top:
      case JSphBoundExtrapZone::DIR_Bottom:
        pt1=pt1-TDouble3(dp2,dp2,0);
        pt2=pt2+TDouble3(dp2,dp2,0);
      break;
      case JSphBoundExtrapZone::DIR_Left:
      case JSphBoundExtrapZone::DIR_Right:
        pt1=pt1-TDouble3(0,dp2,dp2);
        pt2=pt2+TDouble3(0,dp2,dp2);
      break;
      case JSphBoundExtrapZone::DIR_Front:
      case JSphBoundExtrapZone::DIR_Back:
        pt1=pt1-TDouble3(dp2,0,dp2);
        pt2=pt2+TDouble3(dp2,0,dp2);
      break;
    }
    List[c]->MkBound;
    shapes.push_back(JFormatFiles2::DefineShape_Box(pt1,pt2-pt1,mkbound,0)); //-Limit box.
  }
  if(GetCount()){
    const string filevtk=AppInfo.GetDirOut()+"CfgBoundExtrap_Limit.vtk";
    JFormatFiles2::SaveVtkShapes(filevtk,"mkbound","",shapes);
    Log->AddFileInfo(filevtk,"Saves VTK file with BoundExtrap configurations.");
  }
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSphBoundExtrap::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  Log->Printf("DetermLimit: %g",DetermLimit);
  for(unsigned c=0;c<GetCount();c++){
    const JSphBoundExtrapZone* zo=List[c];
    Log->Printf("MkZone_%u (mkfluid:%u)",zo->IdZone,zo->MkBound);
    std::vector<std::string> lines;
    zo->GetConfig(lines);
    for(unsigned i=0;i<unsigned(lines.size());i++)Log->Print(string("  ")+lines[i]);
  }
  if(!txfoot.empty())Log->Print(txfoot);
}











