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

/// \file JSphBoundCorr.cpp \brief Implements the class \ref JSphBoundCorr.

#include "JSphBoundCorr.h"
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
//# JSphBoundCorrZone
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphBoundCorrZone::JSphBoundCorrZone(JLog2 *log,unsigned idzone,word mkbound
 ,TpDirection autodir,tdouble3 limitpos,tdouble3 direction)
 :Log(log),IdZone(idzone),MkBound(mkbound)
{
  ClassName="JSphBoundCorrZone";
  Reset();
  AutoDir=autodir;
  LimitPos=limitpos;
  Direction=direction;
  Plane=ToTFloat4(fmath::PlanePtVec(LimitPos,Direction));
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphBoundCorrZone::~JSphBoundCorrZone(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphBoundCorrZone::Reset(){
  BoundCode=0;
  AutoDir=DIR_None;
  LimitPos=Direction=TDouble3(0);
  Plane=TFloat4(0);
}

//==============================================================================
/// Configures BoundCode.
//==============================================================================
void JSphBoundCorrZone::ConfigBoundCode(typecode boundcode){
  if(BoundCode)RunException("ConfigBoundCode",fun::PrintStr("BoundCode was already configured for mkbound=%u.",MkBound));
  BoundCode=boundcode;
}

//==============================================================================
/// Configures BoundCode.
//==============================================================================
void JSphBoundCorrZone::ConfigAutoLimit(double halfdp,tdouble3 pmin,tdouble3 pmax){
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
/// Applies motion to direction data.
//==============================================================================
void JSphBoundCorrZone::RunMotion(bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix){
  if(simple)LimitPos=LimitPos+msimple;
  else{
    const tdouble3 limitpos0=LimitPos;
    LimitPos=MatrixMulPoint(mmatrix,limitpos0);
    Direction=MatrixMulPoint(mmatrix,limitpos0+Direction)-LimitPos;
  }
  Plane=ToTFloat4(fmath::PlanePtVec(LimitPos,Direction));
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphBoundCorrZone::GetConfig(std::vector<std::string> &lines)const{
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
//# JSphBoundCorr
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphBoundCorr::JSphBoundCorr(bool cpu,double dp,JLog2 *log,JXml *sxml,const std::string &place,const JSphMk *mkinfo)
  :Cpu(cpu),Dp(dp),Log(log)
{
  ClassName="JSphBoundCorr";
  Reset();
  LoadXml(sxml,place);
  UpdateMkCode(mkinfo);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphBoundCorr::~JSphBoundCorr(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphBoundCorr::Reset(){
  DetermLimit=0;
  ExtrapolateMode=0;
  for(int c=0;c<List.size();c++)delete List[c];
  List.clear();
  UseMotion=false;
}

//==============================================================================
/// Returns true if mkbound value is already configured.
//==============================================================================
bool JSphBoundCorr::ExistMk(word mkbound)const{
  bool ret=false;
  for(unsigned c=0;c<List.size() && !ret;c++)ret=(List[c]->MkBound==mkbound);
  return(ret);
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSphBoundCorr::LoadXml(JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of configurations in the XML node.
//==============================================================================
void JSphBoundCorr::ReadXml(const JXml *sxml,TiXmlElement* lis){
  const char met[]="ReadXml";
  //-Loads value determlimit.
  if(sxml->CountElements(lis,"determlimit")>1)sxml->ErrReadElement(lis,"determlimit",false,"Several definitions for this value.");
  DetermLimit=sxml->ReadElementFloat(lis,"determlimit","value",true,1e-3f);
  //-Loads ExtrapolateMode.
  ExtrapolateMode=sxml->ReadElementInt(lis,"extrapolatemode","value",true,1);
  if(ExtrapolateMode>3)ExtrapolateMode=3;
  if(ExtrapolateMode<1)ExtrapolateMode=1;
  if(ExtrapolateMode<2 && Cpu)ExtrapolateMode=2;

  //-Loads list of inputs.
  TiXmlElement* ele=lis->FirstChildElement("mkzone"); 
  while(ele){
    const word mkbound=sxml->GetAttributeWord(ele,"mkbound");
    tdouble3 limitpos=TDouble3(0);
    tdouble3 direction=TDouble3(0);
    JSphBoundCorrZone::TpDirection autodir=JSphBoundCorrZone::DIR_None;
    string autodirtx=fun::StrLower(sxml->ReadElementStr(ele,"autoconfig","direction",true));
    if(autodirtx.empty()){
      limitpos=sxml->ReadElementDouble3(ele,"limitpoint");
      direction=sxml->ReadElementDouble3(ele,"direction");
    }
    else{
      if     (autodirtx=="top"   )autodir=JSphBoundCorrZone::DIR_Top;
      else if(autodirtx=="bottom")autodir=JSphBoundCorrZone::DIR_Bottom;
      else if(autodirtx=="left"  )autodir=JSphBoundCorrZone::DIR_Left;
      else if(autodirtx=="right" )autodir=JSphBoundCorrZone::DIR_Right;
      else if(autodirtx=="front" )autodir=JSphBoundCorrZone::DIR_Front;
      else if(autodirtx=="back"  )autodir=JSphBoundCorrZone::DIR_Back;
      if(autodir==JSphBoundCorrZone::DIR_None)sxml->ErrReadElement(ele,"autoconfig",false,"Direction label is invalid.");
    }
    if(ExistMk(mkbound))RunException(met,"An input already exists for the same mkbound.");
    JSphBoundCorrZone *zo=new JSphBoundCorrZone(Log,GetCount(),mkbound,autodir,limitpos,direction);
    List.push_back(zo);
    ele=ele->NextSiblingElement("mkzone");
  }
}

//==============================================================================
/// Updates BoundCode of each configuration.
//==============================================================================
void JSphBoundCorr::UpdateMkCode(const JSphMk *mkinfo){
  const char met[]="UpdateMkCode";
  for(unsigned c=0;c<GetCount();c++){
    const word mkbound=List[c]->MkBound;
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(List[c]->MkBound);
    if(cmk<mkinfo->Size() && (CODE_IsFixed(mkinfo->Mkblock(cmk)->Code) || CODE_IsMoving(mkinfo->Mkblock(cmk)->Code))){
      List[c]->ConfigBoundCode(mkinfo->Mkblock(cmk)->Code);
      if(CODE_IsMoving(mkinfo->Mkblock(cmk)->Code))UseMotion=true;
    }
    else RunException(met,fun::PrintStr("MkBound value (%u) is not a Mk fixed boundary valid.",List[c]->MkBound));
  }
}

//==============================================================================
/// Run automatic configuration of LimitPos and Direction for each configuration
/// and saves VTK file with limit configuration.
//==============================================================================
void JSphBoundCorr::RunAutoConfig(const JSphMk *mkinfo){
  const char met[]="RunAutoConfig";
  for(unsigned c=0;c<GetCount();c++){
    const word mkbound=List[c]->MkBound;
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(List[c]->MkBound);
    if(cmk<mkinfo->Size()){
      const tdouble3 pmin=mkinfo->Mkblock(cmk)->GetPosMin();
      const tdouble3 pmax=mkinfo->Mkblock(cmk)->GetPosMax();
      List[c]->ConfigAutoLimit(Dp/2.,pmin,pmax);
    }
    else RunException(met,fun::PrintStr("MkBound value (%u) is not a Mk fixed boundary valid.",List[c]->MkBound));
  }
  SaveVtkConfig(Dp,-1);
}

//==============================================================================
/// Saves VTK file with LimitPos and Direction for each configuration.
//==============================================================================
void JSphBoundCorr::SaveVtkConfig(double dp,int part)const{
  const double sizequad=dp*16;
  const double sizedir=dp*4;
  std::vector<JFormatFiles2::StShapeData> shapes;
  for(unsigned c=0;c<GetCount();c++){
    const JSphBoundCorrZone* zo=List[c];
    const int mkbound=zo->MkBound;
    const tdouble3 ps=zo->GetLimitPos();
    const tdouble3 ve=fmath::VecUnitary(zo->GetDirection());
    const tdouble3 v1=fmath::VecOrthogonal2(ve,sizequad,true);
    const tdouble3 v2=fmath::VecUnitary(fmath::ProductVec(ve,v1))*sizequad;
    const tdouble3 p0=ps-(v1/2)-(v2/2);
    const tdouble3 p1=p0+v1;
    const tdouble3 p2=p0+v1+v2;
    const tdouble3 p3=p0+v2;
    //-Adds limit quad.
    shapes.push_back(JFormatFiles2::DefineShape_Line(p0,p1,mkbound,0));
    shapes.push_back(JFormatFiles2::DefineShape_Line(p1,p2,mkbound,0));
    shapes.push_back(JFormatFiles2::DefineShape_Line(p2,p3,mkbound,0));
    shapes.push_back(JFormatFiles2::DefineShape_Line(p3,p0,mkbound,0));
    shapes.push_back(JFormatFiles2::DefineShape_Quad(p0,p1,p2,p3,mkbound,0));
    //-Adds direction line.
    shapes.push_back(JFormatFiles2::DefineShape_Line(ps,ps+(ve*sizedir),mkbound,0)); //-Direction line.
  }
  if(GetCount()){
    string filevtk;
    if(part<0)filevtk=AppInfo.GetDirOut()+"CfgBoundCorr_Limit.vtk";
    else filevtk=AppInfo.GetDirDataOut()+fun::FileNameSec("CfgBoundCorr_Limit.vtk",part);
    JFormatFiles2::SaveVtkShapes(filevtk,"mkbound","",shapes);
    if(part<0)Log->AddFileInfo(filevtk,"Saves VTK file with BoundCorr configurations.");
  }
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSphBoundCorr::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  Log->Printf("DetermLimit.: %g %s",DetermLimit,(DetermLimit==1e-3f? "(1st order)": (DetermLimit==1e+3f? "(0th order)": " ")));
  Log->Printf("ExtrapolateMode: %s",(ExtrapolateMode==1? "FastSingle": (ExtrapolateMode==2? "Single": (ExtrapolateMode==3? "Double": "???"))));
  for(unsigned c=0;c<GetCount();c++){
    const JSphBoundCorrZone* zo=List[c];
    Log->Printf("MkZone_%u (mkbound:%u)",zo->IdZone,zo->MkBound);
    std::vector<std::string> lines;
    zo->GetConfig(lines);
    for(unsigned i=0;i<unsigned(lines.size());i++)Log->Print(string("  ")+lines[i]);
  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Applies motion to direction data.
//==============================================================================
void JSphBoundCorr::RunMotion(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix){
  for(unsigned c=0;c<GetCount();c++)if(List[c]->MkBound==mkbound)List[c]->RunMotion(simple,msimple,mmatrix);
}

//==============================================================================
/// Saves VTK file with configuration when motion is used.
//==============================================================================
void JSphBoundCorr::SaveData(int part)const{
  if(UseMotion)SaveVtkConfig(Dp,part);
}











