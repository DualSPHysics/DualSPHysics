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

/// \file JSphInOut.cpp \brief Implements the class \ref JSphInOut.

#include "JSphInOutZone.h"
#include "JSphInOutPoints.h"
#include "JSphInOutVel.h"
#include "JSphInOutZsurf.h"
#include "JSphInOutGridData.h"
#include "JXml.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JLinearValue.h"
#include "JDsGaugeSystem.h"

#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JDebugSphGpu.h"
#endif

#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;


//##############################################################################
//# JSphInOutZone
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutZone::JSphInOutZone(bool cpu,unsigned idzone,const StCteSph &csp
  ,const tdouble3 &posmin,const tdouble3 &posmax
  ,const JXml *sxml,TiXmlElement* ele,const std::string &dirdatafile
  ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem)
  :Log(AppInfo.LogPtr()),Cpu(cpu),IdZone(idzone),CSP(csp),MapRealPosMin(posmin),MapRealPosMax(posmax)
{
  ClassName="JSphInOutZone";
  Points=NULL;

  InOutVel=NULL;
  InOutZsurf=NULL;

  Reset();
  ReadXml(sxml,ele,dirdatafile,partsdata,gaugesystem);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOutZone::~JSphInOutZone(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInOutZone::Reset(){
  delete Points; Points=NULL;
  Layers=0;
  InputMode=InInput_Free;
  InputCheck=false;
  RefillingMode=InRefill_SimpleZsurf;

  for(unsigned c=0;c<PTDOMSIZE;c++)PtDom[c]=TDouble3(DBL_MAX);
  BoxLimitMin=BoxLimitMax=TFloat3(FLT_MAX);

  Direction=PtPlane=TDouble3(0);
  Plane=TPlane3f(0);
  NptInit=NpartInit=0;

  VelMode=InVelM_Fixed;
  delete InOutVel; InOutVel=NULL;

  ZsurfMode=InZsurf_Undefined;
  delete InOutZsurf; InOutZsurf=NULL;

  RhopMode=InRhop_Constant;

  ExternalVarInput=false;
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphInOutZone::GetConfig(std::vector<std::string> &lines)const{
  const bool simulate2d=CSP.simulate2d;
  lines.push_back(fun::PrintStr("Zone type: %s",InOutVel->GetInletBehaviourName().c_str()));
  Points->GetConfig(lines);
  if(simulate2d)lines.push_back(fun::PrintStr("InOut Position(x,z): (%g,%g)",PtPlane.x,PtPlane.z));
  else          lines.push_back(fun::PrintStr("InOut Position: (%g,%g,%g)",PtPlane.x,PtPlane.y,PtPlane.z));
  if(simulate2d)lines.push_back(fun::PrintStr("InOut Direction(x,z): (%g,%g)",Direction.x,Direction.z));
  else          lines.push_back(fun::PrintStr("InOut Direction: (%g,%g,%g)",Direction.x,Direction.y,Direction.z));
  lines.push_back(fun::PrintStr("InOut points: %u",NptInit));
  lines.push_back(fun::PrintStr("Initial InOut particles: %u",NpartInit));
  lines.push_back(fun::PrintStr("Layers: %u",Layers));
  lines.push_back(fun::PrintStr("Refilling mode: %s",(RefillingMode==InRefill_SimpleFull? "Simple-Full": (RefillingMode==InRefill_SimpleZsurf? "Simple-Zsurf": (RefillingMode==InRefill_Advanced? "Advanced": "???")))));
  lines.push_back(fun::PrintStr("Input treatment: %s",(InputMode==InInput_Free? "Free (no changes)": (InputMode==InInput_Convert? "Convert fluid": (InputMode==InInput_Remove? "Remove fluid": "???")))));
  lines.push_back(fun::PrintStr("Density mode: %s",(RhopMode==InRhop_Constant? "Constant": (RhopMode==InRhop_Hydrostatic? "Hydrostatic": (RhopMode==InRhop_Extrapolated? "Extrapolated": "???")))));
  InOutVel->GetConfig(lines);
  InOutZsurf->GetConfig(lines);
}

//==============================================================================
/// Checks options, shows warnings and throws exceptions.
//==============================================================================
void JSphInOutZone::CheckConfig()const{
  const bool removezsurf=InOutZsurf->GetRemoveZsurf();
  const TpInBehaviour tvelb=InOutVel->GetInletBehaviour();
  if(RefillingMode==InRefill_SimpleFull && removezsurf)Log->PrintfWarning("Inlet/outlet %d: If RemoveZsurf==true then RefillingMode==SimpleFull is equivalent to SimpleZsurf.",IdZone);
  if(RefillingMode==InRefill_SimpleFull || RefillingMode==InRefill_SimpleZsurf){
    if(tvelb==InVelB_Unknown)Log->PrintfWarning("Inlet/outlet %d: Velocity limits are unknown and RefillingMode==SimpleFull/SimpleZsurf is invalid for reverse flows.",IdZone);
    if(tvelb==InVelB_Reverse)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleFull/SimpleZsurf is invalid for reverse flows.",IdZone));
  }
  if(InputMode==InInput_Free    && tvelb==InVelB_Reverse)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Free is invalid for reverse flows (inlet & outlet).",IdZone));
  if(InputMode==InInput_Remove  && tvelb==InVelB_Reverse)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Remove is invalid for reverse flows (inlet & outlet).",IdZone));
  if(InputMode==InInput_Free    && tvelb==InVelB_Outlet )Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Free is invalid for outlet flows.",IdZone));
  if(InputMode==InInput_Remove  && tvelb==InVelB_Outlet )Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Remove is invalid for outlet flows.",IdZone));
  if(InputMode==InInput_Convert && tvelb==InVelB_Inlet  )Log->PrintfWarning("Inlet/outlet %d: InputTreatment==Convert is not recommended for inlet flows.",IdZone);
  if(InputMode==InInput_Free    && tvelb==InVelB_Unknown)Log->PrintfWarning("Inlet/outlet %d: InputTreatment==Free is not recommended for undefined flows since it is invalid for outlet flows.",IdZone);
  if(InputMode==InInput_Remove  && tvelb==InVelB_Unknown)Log->PrintfWarning("Inlet/outlet %d: InputTreatment==Remove is not recommended for undefined flows since it is invalid for outlet flows.",IdZone);
  if(RhopMode==InRhop_Hydrostatic && ZsurfMode==InZsurf_Undefined)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: Z-Surface is undefined for hydrostatic density calculation.",IdZone));
  if(VelMode==InVelM_Extrapolated)Log->PrintfWarning("Inlet/outlet %d: VelocityMode==Extrapolated is not recommended for inlet (or reverse) flows.",IdZone);
  if(ZsurfMode==InZsurf_Variable && RefillingMode==InRefill_SimpleZsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleZsurf is invalid with ZSurface==Variable.",IdZone));
  if(ZsurfMode==InZsurf_Variable && RefillingMode==InRefill_SimpleFull && removezsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleFull is invalid with ZSurface==Variable and RemoveZsurf==true.",IdZone));
  if(ZsurfMode==InZsurf_Calculated && RefillingMode==InRefill_SimpleZsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleZsurf is invalid with ZSurface==Calculated.",IdZone));
  if(ZsurfMode==InZsurf_Calculated && RefillingMode==InRefill_SimpleFull && removezsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleFull is invalid with ZSurface==Calculated and RemoveZsurf==true.",IdZone));
}

//==============================================================================
/// Loads domain and calculates BoxLimitMin/Max.
//==============================================================================
void JSphInOutZone::LoadDomain(){
  //-Obtains domain data from Points object.
  const tdouble3* ptdom=Points->GetPtDomain();
  for(unsigned c=0;c<10;c++)PtDom[c]=ptdom[c];
  //-Calculates BoxLimit.
  tdouble3 pmin=PtDom[0],pmax=PtDom[0];
  for(unsigned c=1;c<unsigned(CSP.simulate2d? 4: 8);c++){
    const tdouble3 ps=PtDom[c];
    pmin=MinValues(pmin,ps);
    pmax=MaxValues(pmax,ps);
  }
  BoxLimitMin=ToTFloat3(pmin);
  BoxLimitMax=ToTFloat3(pmax);
  //-Adds minimum border to avoid problems in float/double comparisons.
  if(CSP.simulate2d){
    BoxLimitMin.y=float(CSP.simulate2dposy-CSP.dp*.1);
    BoxLimitMax.y=float(CSP.simulate2dposy+CSP.dp*.1);
  }
  ////-Adds border of size Scell.
  //BoxLimitMin=ToTFloat3(pmin-TDouble3(Scell));
  //BoxLimitMax=ToTFloat3(pmax+TDouble3(Scell));
}

//==============================================================================
/// Reads initial configuration in the XML node.
//==============================================================================
void JSphInOutZone::ReadXml(const JXml *sxml,TiXmlElement* ele,const std::string &dirdatafile
  ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem)
{
  //-Checks old configuration.
  if(sxml->ExistsElement(ele,"userefilling"))Run_ExceptioonFile(fun::PrintStr("Inlet/outlet zone %d: <userefilling> is not supported by current inlet/outlet version.",IdZone),sxml->ErrGetFileRow(ele,"userefilling"));
  if(sxml->ExistsElement(ele,"convertfluid"))Run_ExceptioonFile(fun::PrintStr("Inlet/outlet zone %d: <convertfluid> is not supported by current inlet/outlet version.",IdZone),sxml->ErrGetFileRow(ele,"convertfluid"));
  //-Checks element names.
  sxml->CheckElementNames(ele,true,"refilling inputtreatment layers zone2d zone3d imposevelocity imposerhop imposezsurf");
  //-Loads InputMode (InputTreatment).
  {
    const unsigned imode=sxml->ReadElementUnsigned(ele,"inputtreatment","value");
    switch(imode){
      case 0:  InputMode=InInput_Free;     break;
      case 1:  InputMode=InInput_Convert;  break;
      case 2:  InputMode=InInput_Remove;   break;
      default: sxml->ErrReadElement(ele,"inputtreatment",false);
    }
  }
  //-Loads RefillinMode.
  {
    const unsigned imode=sxml->ReadElementUnsigned(ele,"refilling","value");
    switch(imode){
      case 0:  RefillingMode=InRefill_SimpleFull;   break;
      case 1:  RefillingMode=InRefill_SimpleZsurf;  break;
      case 2:  RefillingMode=InRefill_Advanced;     break;
      default: sxml->ErrReadElement(ele,"refilling",false);
    }
  }
  //-Loads layers value.
  unsigned layers=sxml->ReadElementUnsigned(ele,"layers","value");
  if(layers<1)sxml->ErrReadElement(ele,"layers",false,"Minumum number of layers is 1.");
  if(layers>250)sxml->ErrReadElement(ele,"layers",false,"Maximum number of layers is 250.");
  Layers=byte(layers);
  //-Creates inlet points.
  Points=new JSphInOutPoints(CSP.simulate2d,CSP.simulate2dposy,Layers,CSP.dp,0,MapRealPosMin,MapRealPosMax);
  if(CSP.simulate2d){
    TiXmlElement* zone2d=ele->FirstChildElement("zone2d"); 
    if(zone2d)Points->CreatePoints(sxml,zone2d,partsdata);
    else sxml->ErrReadElement(ele,"zone2d",true);
  }
  else{
    TiXmlElement* zone3d=ele->FirstChildElement("zone3d"); 
    if(zone3d)Points->CreatePoints(sxml,zone3d,partsdata);
    else sxml->ErrReadElement(ele,"zone3d",true);
  }
  //-Obtains domain data from Points object.
  LoadDomain();

  //-Obtains information inlet/outlet points.
  Direction=Points->GetDirection();
  PtPlane=PtDom[8];
  if(PtPlane==TDouble3(DBL_MAX))Run_Exceptioon("Reference point in inout plane is invalid.");
  Plane=TPlane3f(fgeo::PlanePtVec(PtPlane,Direction));
  NptInit=Points->GetCount();

  //-Velocity configuration.
  InOutVel=new JSphInOutVel(Cpu,IdZone,CSP,Direction,PtPlane,Points->GetZonePosMin(),Points->GetZonePosMax());
  VelMode=InOutVel->ReadXml(sxml,ele,dirdatafile,gaugesystem,MapRealPosMin.y);

  //-Z-surface configuration.
  InOutZsurf=new JSphInOutZsurf(Cpu,IdZone,CSP,Direction,Points->GetZonePosMin(),Points->GetZonePosMax());
  ZsurfMode=InOutZsurf->ReadXml(sxml,ele,dirdatafile,gaugesystem);

  //-Rhop configuration.
  TiXmlElement *xele=ele->FirstChildElement("imposerhop");
  if(xele){
    const unsigned mode=sxml->GetAttributeUint(xele,"mode",true);
    switch(mode){
      case 0:  RhopMode=InRhop_Constant;      break;
      case 1:  RhopMode=InRhop_Hydrostatic;   break;
      case 2:  RhopMode=InRhop_Extrapolated;  break;
      default: sxml->ErrReadAtrib(xele,"mode",false,"Value is not valid.");
    }
    //-Checks old configuration.
    if(sxml->ExistsElement(ele,"zsurf"  ))Run_ExceptioonFile(fun::PrintStr("Inlet/outlet zone %d: <zsurf> is not supported by current inlet/outlet version.",IdZone),sxml->ErrGetFileRow(ele,"zsurf"));
    if(sxml->ExistsElement(ele,"zbottom"))Run_ExceptioonFile(fun::PrintStr("Inlet/outlet zone %d: <zbottom> is not supported by current inlet/outlet version.",IdZone),sxml->ErrGetFileRow(ele,"zbottom"));
    //-Checks configuration.
    sxml->CheckElementNames(xele,true," ");
  }
  else{//-Default configuration for rhop.
    RhopMode=InRhop_Constant;
  }

  //-Defines initial valid points according to zsurf.
  if(ZsurfMode==InZsurf_Undefined || RefillingMode==InRefill_SimpleFull)Points->SetPointsInit(true);
  else InOutZsurf->SetInitialPoints(Points->GetCount(),Points->GetPoints(),Points->GetPointsInit());
  NpartInit=Points->CountPointsInit()*Layers;

  //-Compute flow velocity factor.
  if(InOutVel->GetFlowActive()){
    if(ZsurfMode!=InZsurf_Fixed)Run_Exceptioon(fun::PrintStr("Inlet/outlet zone %d: The use of flow velocity is only supported by fixed zsurf configuration.",IdZone));
    const unsigned nptok=InOutZsurf->ComputeActivePoints(Points->GetCount(),Points->GetPoints());
    InOutVel->ConfigFlowToVel(nptok);
  }

  InputCheck=(InOutZsurf->GetRemoveZsurf() || InputMode!=InInput_Free);
}

//==============================================================================
/// Returns true when AWAS-velocity is configured.
//==============================================================================
bool JSphInOutZone::Use_AwasVel()const{ 
  return(InOutVel && InOutVel->UseAwasVel());
}

//==============================================================================
/// Returns a byte with information about VelMode, VelProfile and RhopMode.
//==============================================================================
byte JSphInOutZone::GetConfigZone()const{
  const TpInVelProfile velprofile=InOutVel->GetVelProfile();
  if((velprofile&InVelP_MASK)!=velprofile)Run_Exceptioon("VelProfile value is invalid to code using mask.");
  if((VelMode   &InVelM_MASK)!=VelMode   )Run_Exceptioon("VelMode value is invalid to code using mask.");
  if((RhopMode  &InRhop_MASK)!=RhopMode  )Run_Exceptioon("RhopMode value is invalid to code using mask.");
  //if(byte(ZsurfMode) >3)Run_Exceptioon("ZsurfMode value is invalid to code in 2 bits.");
  byte ret=byte(velprofile);
  ret|=byte(VelMode);
  ret|=byte(RhopMode);
  if(InputCheck)ret|=byte(CheckInput_MASK);
  //-Checks coded configuration.
  TpInVelProfile vprof =GetConfigVelProfile(ret);
  TpInVelMode    vmode =GetConfigVelMode(ret);
  TpInRhopMode   rmode =GetConfigRhopMode(ret);
  bool         checkf=GetConfigCheckInputDG(ret);
  if(vprof!=velprofile || vmode!=VelMode || rmode!=RhopMode || checkf!=InputCheck)
    Run_Exceptioon("Coded configuration is not right.");
  return(ret);
}

//==============================================================================
/// Returns a byte with information about refilling mode and RemoveZsurf.
//==============================================================================
byte JSphInOutZone::GetConfigUpdate()const{
  byte ret=byte(RefillingMode==InRefill_Advanced  ? RefillAdvanced_MASK: 0);
  ret= ret|byte(RefillingMode==InRefill_SimpleFull? RefillSpFull_MASK  : 0);
  ret= ret|byte(InputMode==InInput_Remove         ? RemoveInput_MASK : 0);
  ret= ret|byte(InOutZsurf->GetRemoveZsurf()      ? RemoveZsurf_MASK : 0);
  ret= ret|byte(InputMode==InInput_Convert        ? ConvertInput_MASK: 0);
  return(ret);
}

//==============================================================================
/// Loads positon of inlet/outlet points.
//==============================================================================
unsigned JSphInOutZone::LoadInletPoints(tdouble3 *pos){
  const tdouble3 *ptpos=Points->GetPoints();
  const tdouble3 dir=Direction*(CSP.dp/2);
  for(unsigned cp=0;cp<NptInit;cp++)pos[cp]=ptpos[cp]+dir;
  return(NptInit);
}

//==============================================================================
/// Loads positon of initial inlet/outlet particles.
//==============================================================================
void JSphInOutZone::LoadInitialParticles(unsigned npartinit,tdouble3 *pos){
  const tdouble3 *ptpos=Points->GetPoints();
  const byte *ptok=Points->GetPointsInit();
  unsigned p=0;
  for(unsigned layer=0;layer<Layers;layer++){
    const tdouble3 dir=Direction*(-CSP.dp*layer);
    for(unsigned cp=0;cp<NptInit;cp++)if(ptok[cp]){
      if(p<npartinit)pos[p]=ptpos[cp]+dir;
      p++;
    }
  }
  if(p!=npartinit)Run_Exceptioon("Number of initial particles does not match.");
  //Log->Printf(" LoadInletParticles--> p:%u  NptInit:%u  Layers:%u",p,NptInit,Layers);
  Points->ResetPoints();//-Frees memory of points.
}

//==============================================================================
/// Returns velocity according profile configuration.
//==============================================================================
float JSphInOutZone::CalcVel(TpInVelProfile vprof,const tfloat4 &vdata,double posz){
  float vel=0;
  if(vprof==InVelP_Uniform)vel=vdata.x;
  else if(vprof==InVelP_Linear){
    const float m=vdata.x;
    const float b=vdata.y;
    vel=m*float(posz)+b;
  }
  else if(vprof==InVelP_Parabolic){
    const float a=vdata.x;
    const float b=vdata.y;
    const float c=vdata.z;
    vel=a*float(posz)*float(posz)+b*float(posz)+c;
  }
  return(vel);
}

//==============================================================================
/// Indicates if position is inside inlet zone using BoxLimitMin/Max.
//==============================================================================
bool JSphInOutZone::InZoneBox(const tfloat3 &ps)const{
  return(BoxLimitMin.x<=ps.x && ps.x<=BoxLimitMax.x && BoxLimitMin.y<=ps.y && ps.y<=BoxLimitMax.y && BoxLimitMin.z<=ps.z && ps.z<=BoxLimitMax.z);
}

//==============================================================================
/// Indicates if position is inside inlet zone, using plane and BoxLimitMin/Max.
//==============================================================================
bool JSphInOutZone::InZone(bool useboxlimit,const tfloat3 &ps)const{
  return((!useboxlimit || InZoneBox(ps)) && fgeo::PlanePoint(Plane,ps)<0);
}


