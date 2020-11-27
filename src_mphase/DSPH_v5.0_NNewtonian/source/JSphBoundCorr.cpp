//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "JDsPartsInit.h"
#include "JXml.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JMatrix4.h"
#include "JLinearValue.h"
#include "JSaveCsv2.h"
#include "JRangeFilter.h"
#include "JVtkLib.h"

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
JSphBoundCorrZone::JSphBoundCorrZone(unsigned idzone,word mkbound,TpDirection autodir
 ,double autodpfactor,tdouble3 limitpos,tdouble3 direction)
 :Log(AppInfo.LogPtr()),IdZone(idzone),MkBound(mkbound)
{
  ClassName="JSphBoundCorrZone";
  Reset();
  AutoDir=autodir;
  AutoDpFactor=autodpfactor;
  LimitPos=limitpos;
  Direction=direction;
  Plane=TPlane3f(fgeo::PlanePtVec(LimitPos,Direction));
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
  AutoDpFactor=0;
  LimitPos=Direction=TDouble3(0);
  Plane=TPlane3f(0);
}

//==============================================================================
/// Configures BoundCode.
//==============================================================================
void JSphBoundCorrZone::ConfigBoundCode(typecode boundcode){
  if(BoundCode)Run_Exceptioon(fun::PrintStr("BoundCode was already configured for mkbound=%u.",MkBound));
  BoundCode=boundcode;
}

//==============================================================================
/// Run automatic configuration of LimitPos and Direction.
//==============================================================================
void JSphBoundCorrZone::ConfigAuto(const JDsPartsInit *partsdata){
  if(AutoDir!=DIR_None){
    //-Calculates limits of MK particles.
    tdouble3 pmin=TDouble3(0);
    tdouble3 pmax=TDouble3(0);
    const JSphMk* mkinfo=partsdata->GetMkInfo();
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(MkBound);
    if(cmk<mkinfo->Size()){
      pmin=mkinfo->Mkblock(cmk)->GetPosMin();
      pmax=mkinfo->Mkblock(cmk)->GetPosMax();
    }
    else Run_Exceptioon(fun::PrintStr("MkBound value (%u) is not a Mk boundary valid.",MkBound));
    const tdouble3 pmed=(pmin+pmax)/2.;
    if(AutoDir==DIR_Defined){
      const typecode codesel=mkinfo->Mkblock(cmk)->Code;
      const tplane3d pla=fgeo::PlanePtVec(pmed,Direction);
      const tdouble3* pos=partsdata->GetPos();
      const typecode* code=partsdata->GetCode();
      const unsigned np=partsdata->GetNp();
      double dismax=-DBL_MAX;
      for(unsigned p=0;p<np;p++)if(code[p]==codesel){
        const double dist=fgeo::PlaneDistSign(pla,pos[p]);
        if(dist>dismax)dismax=dist;
      }
      if(dismax==-DBL_MAX)Run_Exceptioon(fun::PrintStr("It was not possible to calculate the limit position for MkBound=%u automatically.",MkBound));
      LimitPos=pmed+(Direction*(dismax+(partsdata->Dp*AutoDpFactor)));
    }
    else{
      const double halfdp=partsdata->Dp/2.;
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
    }
    Plane=TPlane3f(fgeo::PlanePtVec(LimitPos,Direction));
  }
}

//==============================================================================
/// Applies motion to direction data.
//==============================================================================
void JSphBoundCorrZone::RunMotion(const StMotionData& m){
  if(m.type==MOTT_Linear)LimitPos=LimitPos+m.linmov;
  if(m.type==MOTT_Matrix){
    const tdouble3 limitpos0=LimitPos;
    LimitPos=MatrixMulPoint(m.matmov,limitpos0);
    Direction=MatrixMulPoint(m.matmov,limitpos0+Direction)-LimitPos;
  }
  Plane=TPlane3f(fgeo::PlanePtVec(LimitPos,Direction));
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphBoundCorrZone::GetConfig(std::vector<std::string> &lines)const{
  bool autodir=false;
  bool autopos=false;
  if(AutoDir==DIR_Defined)autopos=true;
  if(AutoDir!=DIR_Defined && AutoDir!=DIR_None)autodir=autopos=true;
  lines.push_back(fun::PrintStr("Limit position.: (%g,%g,%g)%s",LimitPos.x,LimitPos.y,LimitPos.z,(autopos? " (automatic)": "")));
  lines.push_back(fun::PrintStr("Fluid direction: (%g,%g,%g)%s",Direction.x,Direction.y,Direction.z,(autodir? " (automatic)": "")));
}


//##############################################################################
//# JSphBoundCorr
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphBoundCorr::JSphBoundCorr(bool cpu,double dp,const JXml *sxml
  ,const std::string &place,const JSphMk *mkinfo)
  :Log(AppInfo.LogPtr()),Cpu(cpu),Dp(dp),SaveMotionVtk(false)
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
void JSphBoundCorr::LoadXml(const JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of configurations in the XML node.
//==============================================================================
void JSphBoundCorr::ReadXml(const JXml *sxml,TiXmlElement* lis){
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
    if(sxml->CheckElementActive(ele)){
      std::vector<unsigned> mkbounds;
      JRangeFilter rg(sxml->GetAttributeStr(ele,"mkbound"));
      rg.GetValues(mkbounds);
      //const word mkbound=sxml->GetAttributeWord(ele,"mkbound");
      const bool autoconfig=sxml->ExistsElement(ele,"autoconfig");
      const bool autolimitpoint=sxml->ExistsElement(ele,"autolimitpoint");
      if(autoconfig && autolimitpoint)sxml->ErrReadElement(lis,"determlimit",false,"Several configuration modes. \'autoconfig\' and \'autolimitpoint\' definitions are not compatible.");
      double autodpfactor=0;
      tdouble3 limitpoint=TDouble3(0);
      tdouble3 direction=TDouble3(0);
      JSphBoundCorrZone::TpDirection autodir=JSphBoundCorrZone::DIR_None;
      if(autoconfig){
        string autodirtx=fun::StrLower(sxml->ReadElementStr(ele,"autoconfig","direction",true));
        if     (autodirtx=="top"   )autodir=JSphBoundCorrZone::DIR_Top;
        else if(autodirtx=="bottom")autodir=JSphBoundCorrZone::DIR_Bottom;
        else if(autodirtx=="left"  )autodir=JSphBoundCorrZone::DIR_Left;
        else if(autodirtx=="right" )autodir=JSphBoundCorrZone::DIR_Right;
        else if(autodirtx=="front" )autodir=JSphBoundCorrZone::DIR_Front;
        else if(autodirtx=="back"  )autodir=JSphBoundCorrZone::DIR_Back;
        if(autodir==JSphBoundCorrZone::DIR_None)sxml->ErrReadElement(ele,"autoconfig",false,"Direction label is invalid.");
      }
      else if(autolimitpoint){
        direction=sxml->ReadElementDouble3(ele,"direction");
        if(direction==TDouble3(0))sxml->ErrReadElement(ele,"direction",false,"Direction vector is invalid.");
        else direction=fgeo::VecUnitary(direction);
        autodir=JSphBoundCorrZone::DIR_Defined;
        autodpfactor=sxml->ReadElementDouble(ele,"autolimitpoint","dpfactor");
      }
      else{
        direction=sxml->ReadElementDouble3(ele,"direction");
        if(direction==TDouble3(0))sxml->ErrReadElement(ele,"direction",false,"Direction vector is invalid.");
        else direction=fgeo::VecUnitary(direction);
        limitpoint=sxml->ReadElementDouble3(ele,"limitpoint");
      }
      const unsigned nmkbounds=unsigned(mkbounds.size());
      for(unsigned cmk=0;cmk<nmkbounds;cmk++){
        const word mkbound=word(mkbounds[cmk]);
        if(ExistMk(mkbound))Run_Exceptioon(fun::PrintStr("An input already exists for the same mkbound=%u.",mkbound));
        JSphBoundCorrZone *zo=new JSphBoundCorrZone(GetCount(),mkbound,autodir,autodpfactor,limitpoint,direction);
        List.push_back(zo);
      }
    }
    ele=ele->NextSiblingElement("mkzone");
  }
}

//==============================================================================
/// Updates BoundCode of each configuration.
//==============================================================================
void JSphBoundCorr::UpdateMkCode(const JSphMk *mkinfo){
  for(unsigned c=0;c<GetCount();c++){
    const word mkbound=List[c]->MkBound;
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(List[c]->MkBound);
    if(cmk<mkinfo->Size() && (CODE_IsFixed(mkinfo->Mkblock(cmk)->Code) || CODE_IsMoving(mkinfo->Mkblock(cmk)->Code))){
      List[c]->ConfigBoundCode(mkinfo->Mkblock(cmk)->Code);
      if(CODE_IsMoving(mkinfo->Mkblock(cmk)->Code))UseMotion=true;
    }
    else Run_Exceptioon(fun::PrintStr("MkBound value (%u) is not a Mk fixed boundary valid.",List[c]->MkBound));
  }
}

//==============================================================================
/// Run automatic configuration of LimitPos and Direction for each configuration
/// and saves VTK file with limit configuration.
//==============================================================================
void JSphBoundCorr::RunAutoConfig(const JDsPartsInit *partsdata){
  for(unsigned c=0;c<GetCount();c++)List[c]->ConfigAuto(partsdata);
  SaveVtkConfig(Dp,-1);
}

//==============================================================================
/// Saves VTK file with LimitPos and Direction for each configuration.
//==============================================================================
void JSphBoundCorr::SaveVtkConfig(double dp,int part)const{
  const double sizequad=dp*16;
  const double sizedir=dp*4;
  JVtkLib sh;
  for(unsigned c=0;c<GetCount();c++){
    const JSphBoundCorrZone* zo=List[c];
    const int mkbound=zo->MkBound;
    const tdouble3 ps=zo->GetLimitPos();
    const tdouble3 ve=fgeo::VecUnitary(zo->GetDirection());
    const tdouble3 v1=fgeo::VecOrthogonal2(ve,sizequad,true);
    const tdouble3 v2=fgeo::VecUnitary(fgeo::ProductVec(ve,v1))*sizequad;
    const tdouble3 p0=ps-(v1/2)-(v2/2);
    const tdouble3 p1=p0+v1;
    const tdouble3 p2=p0+v1+v2;
    const tdouble3 p3=p0+v2;
    //-Adds limit quad.
    sh.AddShapeQuad(p0,p1,p2,p3,mkbound);
    sh.AddShapeQuadWire(p0,p1,p2,p3,mkbound);
    //-Adds direction line.
    sh.AddShapeLine(ps,ps+(ve*sizedir),mkbound); //-Direction line.
  }
  if(GetCount()){
    string filevtk;
    if(part<0)filevtk=AppInfo.GetDirOut()+"CfgBoundCorr_Limit.vtk";
    else filevtk=AppInfo.GetDirDataOut()+fun::FileNameSec("CfgBoundCorr_Limit.vtk",part);
    sh.SaveShapeVtk(filevtk,"MkBound");
    if(part<0)Log->AddFileInfo(filevtk,"Saves VTK file with BoundCorr configurations.");
  }
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSphBoundCorr::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  Log->Printf("DetermLimit....: %g %s",DetermLimit,(DetermLimit==1e-3f? "(1st order)": (DetermLimit==1e+3f? "(0th order)": " ")));
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
void JSphBoundCorr::RunMotion(const StMotionData& motiondata){
  for(unsigned c=0;c<GetCount();c++)if(List[c]->MkBound==motiondata.mkbound)List[c]->RunMotion(motiondata);
}

//==============================================================================
/// Saves VTK file with configuration when motion is used.
//==============================================================================
void JSphBoundCorr::SaveData(int part)const{
  if(UseMotion && SaveMotionVtk)SaveVtkConfig(Dp,part);
}











