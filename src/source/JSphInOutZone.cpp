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
#include "JSphInOutGridData.h"
//#include "JSphCpu.h"
//#include "JSphMk.h"
#include "JXml.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "FunctionsGeo3d.h"
//#include "JMatrix4.h"
#include "JLinearValue.h"
//#include "JRangeFilter.h"
//#include "JDataArrays.h"
//#include "JVtkLib.h"
//#include "JSimpleNeigs.h"
//#include "JTimeControl.h"
#include "JDsGaugeSystem.h"
//#include "JNumexLib.h"

#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  //#include "JSphGpu_InOut_iker.h"
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
JSphInOutZone::JSphInOutZone(bool cpu,JLog2 *log,unsigned idzone,bool simulate2d
  ,double simulate2dposy,double dp,const tdouble3 &posmin,const tdouble3 &posmax
  ,float gravityz,const JXml *sxml,TiXmlElement* ele,const std::string &dirdatafile
  ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem)
 :Cpu(cpu),Log(log),IdZone(idzone),Simulate2D(simulate2d),Simulate2DPosY(simulate2dposy)
 ,Dp(dp),MapRealPosMin(posmin),MapRealPosMax(posmax),GravityZ(gravityz)
{
  ClassName="JSphInOutZone";
  Points=NULL;
  TimeInputVelData=NULL;
  InputVelGrid=NULL;
  TimeInputZsurf=NULL;
  PtzPos=NULL;
  #ifdef _WITHGPU
    PtzPosg=NULL; PtzAuxg=NULL; PtzAux=NULL;
  #endif

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
  delete TimeInputVelData; TimeInputVelData=NULL;
  delete InputVelGrid; InputVelGrid=NULL;
  ResetZVelGrid=false;
  TimeVelIdx0=TimeVelIdx1=UINT_MAX; 
  SaveVelProfile=true;
  delete TimeInputZsurf; TimeInputZsurf=NULL;
  Layers=0;
  InputMode=TIN_Free;
  InputCheck=false;
  RefillingMode=TFI_SimpleZsurf;
  for(unsigned c=0;c<10;c++)PtDom[c]=TDouble3(DBL_MAX);
  BoxLimitMin=BoxLimitMax=TFloat3(FLT_MAX);
  VelMode=MVEL_Fixed;
  VelProfile=PVEL_Constant;
  InputVel=InputVel2=InputVel3=0;
  InputVelPosz=InputVelPosz2=InputVelPosz3=0;
  VelMin=VelMax=0;
  RhopMode=MRHOP_Constant;
  ZsurfMode=ZSURF_Undefined;
  InputZsurf=InputZbottom=FLT_MAX;
  RemoveZsurf=false;
  SvVtkZsurf=false;
  ExternalVarInput=false;
  Direction=PtPlane=TDouble3(0);
  Plane=TPlane3f(0);
  NptInit=NpartInit=0;
  delete[] PtzPos;  PtzPos=NULL;
  #ifdef _WITHGPU
    if(PtzPosg)cudaFree(PtzPosg);  PtzPosg=NULL;
    if(PtzAuxg)cudaFree(PtzAuxg);  PtzAuxg=NULL;
    delete[] PtzAux;  PtzAux=NULL;
  #endif
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphInOutZone::GetConfig(std::vector<std::string> &lines)const{
  const byte tvel=(VelMin>=0 && VelMax>=0? 1: (VelMin<=0 && VelMax<=0? 2: (VelMin<0 && VelMax>=0 && VelMin!=-FLT_MAX && VelMax!=FLT_MAX? 3: 0))); //-0=Unknown, 1=Inlet, 2=Outlet, 3=Reverse.
  lines.push_back(fun::PrintStr("Zone type: %s",(tvel==1? "Inlet (velocity>=0)": (tvel==2? "Outlet (velocity<0)": (tvel==3? "Inlet & Outlet (velocity>0 and velocity<0)": "Unknown")))));
  Points->GetConfig(lines);
  if(Simulate2D)lines.push_back(fun::PrintStr("InOut Position(x,z): (%g,%g)",PtPlane.x,PtPlane.z));
  else          lines.push_back(fun::PrintStr("InOut Position: (%g,%g,%g)",PtPlane.x,PtPlane.y,PtPlane.z));
  if(Simulate2D)lines.push_back(fun::PrintStr("InOut Direction(x,z): (%g,%g)",Direction.x,Direction.z));
  else          lines.push_back(fun::PrintStr("InOut Direction: (%g,%g,%g)",Direction.x,Direction.y,Direction.z));
  lines.push_back(fun::PrintStr("InOut points: %u",NptInit));
  lines.push_back(fun::PrintStr("Initial InOut particles: %u",NpartInit));
  lines.push_back(fun::PrintStr("Layers: %u",Layers));
  lines.push_back(fun::PrintStr("Refilling mode: %s",(RefillingMode==TFI_SimpleFull? "Simple-Full": (RefillingMode==TFI_SimpleZsurf? "Simple-Zsurf": (RefillingMode==TFI_Advanced? "Advanced": "???")))));
  lines.push_back(fun::PrintStr("Input treatment: %s",(InputMode==TIN_Free? "Free (no changes)": (InputMode==TIN_Convert? "Convert fluid": (InputMode==TIN_Remove? "Remove fluid": "???")))));
  if(VelMode==MVEL_Fixed || VelMode==MVEL_Variable){
    lines.push_back(fun::PrintStr("Velocity mode: %s",(VelMode==MVEL_Fixed? "Fixed": "Variable")));
    if(VelMode==MVEL_Variable && !TimeInputVelData->GetFile().empty())lines.push_back(fun::PrintStr("  Velocity file: %s",TimeInputVelData->GetFile().c_str()));
    if(VelProfile==PVEL_Constant )lines.push_back(fun::PrintStr("  Velocity profile: Constant %g",InputVel));
    if(VelProfile==PVEL_Linear   )lines.push_back(fun::PrintStr("  Velocity profile: Linear %g(z=%g), %g(z=%g)",InputVel,InputVelPosz,InputVel2,InputVelPosz2));
    if(VelProfile==PVEL_Parabolic)lines.push_back(fun::PrintStr("  Velocity profile: Parabolic %g(z=%g), %g(z=%g), %g(z=%g)",InputVel,InputVelPosz,InputVel2,InputVelPosz2,InputVel3,InputVelPosz3));
  }
  else if(VelMode==MVEL_Extrapolated)lines.push_back("Velocity mode: Extrapolated");
  else if(VelMode==MVEL_Interpolated){
    lines.push_back("Velocity mode: Interpolated");
    lines.push_back(fun::PrintStr("  Velocity file: %s",InputVelGrid->GetFile().c_str()));
    lines.push_back(fun::PrintStr("  Reset Z velocity: %s",(ResetZVelGrid? "True": "False")));
  }
  else lines.push_back("???");
  lines.push_back(fun::PrintStr("Density mode: %s",(RhopMode==MRHOP_Constant? "Constant": (RhopMode==MRHOP_Hydrostatic? "Hydrostatic": (RhopMode==MRHOP_Extrapolated? "Extrapolated": "???")))));
  lines.push_back(fun::PrintStr("Z-Surface mode: %s",(ZsurfMode==ZSURF_Undefined? "Undefined": (ZsurfMode==ZSURF_Fixed? "Fixed": (ZsurfMode==ZSURF_Variable? "Variable": (ZsurfMode==ZSURF_Calculated? "Calculated": "???"))))));
  if(ZsurfMode==ZSURF_Variable && !TimeInputZsurf->GetFile().empty())lines.push_back(fun::PrintStr("Zsurf file: %s",TimeInputZsurf->GetFile().c_str()));
  if(ZsurfMode!=ZSURF_Undefined || RhopMode==MRHOP_Hydrostatic){
    if(InputZsurf!=FLT_MAX)  lines.push_back(fun::PrintStr("  Z-Surface value (initial): %g",InputZsurf));
    else                     lines.push_back(fun::PrintStr("  Z-Surface value (initial): Undefined"));
    if(InputZbottom!=FLT_MAX)lines.push_back(fun::PrintStr("  Z-Bottom value: %g",InputZbottom));
    else                     lines.push_back(fun::PrintStr("  Z-Bottom value (initial): Undefined"));
  }
  if(ZsurfMode!=ZSURF_Undefined){
    lines.push_back(fun::PrintStr("  Z-Surface Remove: %s",(RemoveZsurf? "True": "False")));
    lines.push_back(fun::PrintStr("  Z-Surface SaveVTK: %s",(SvVtkZsurf? "True": "False")));
  }
}

//==============================================================================
/// Checks options, shows warnings and throws exceptions.
//==============================================================================
void JSphInOutZone::CheckConfig()const{
  const byte tvel=(VelMin>=0 && VelMax>=0? 1: (VelMin<=0 && VelMax<=0? 2: (VelMin<0 && VelMax>=0 && VelMin!=-FLT_MAX && VelMax!=FLT_MAX? 3: 0))); //-0=Unknown, 1=Inlet, 2=Outlet, 3=Reverse.
  if(RefillingMode==TFI_SimpleFull && RemoveZsurf)Log->PrintfWarning("Inlet/outlet %d: If RemoveZsurf==true then RefillingMode==SimpleFull is equivalent to SimpleZsurf.",IdZone);
  if(RefillingMode==TFI_SimpleFull || RefillingMode==TFI_SimpleZsurf){
    if(tvel==0)Log->PrintfWarning("Inlet/outlet %d: Velocity limits are unknown and RefillingMode==SimpleFull/SimpleZsurf is invalid for reverse flows.",IdZone);
    if(tvel==3)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleFull/SimpleZsurf is invalid for reverse flows.",IdZone));
  }
  if(InputMode==TIN_Free    && tvel==3)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Free is invalid for reverse flows (inlet & outlet).",IdZone));
  if(InputMode==TIN_Remove  && tvel==3)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Remove is invalid for reverse flows (inlet & outlet).",IdZone));
  if(InputMode==TIN_Free    && tvel==2)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Free is invalid for outlet flows.",IdZone));
  if(InputMode==TIN_Remove  && tvel==2)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: InputTreatment==Remove is invalid for outlet flows.",IdZone));
  if(InputMode==TIN_Convert && tvel==1)Log->PrintfWarning("Inlet/outlet %d: InputTreatment==Convert is not recommended for inlet flows.",IdZone);
  if(InputMode==TIN_Free    && tvel==0)Log->PrintfWarning("Inlet/outlet %d: InputTreatment==Free is not recommended for undefined flows since it is invalid for outlet flows.",IdZone);
  if(InputMode==TIN_Remove  && tvel==0)Log->PrintfWarning("Inlet/outlet %d: InputTreatment==Remove is not recommended for undefined flows since it is invalid for outlet flows.",IdZone);
  if(RhopMode==MRHOP_Hydrostatic && (InputZsurf==FLT_MAX || InputZbottom==FLT_MAX))Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: Z-Surface or Z-Bottom is undefined for hydrostatic density calculation.",IdZone));
  if(VelMode==MVEL_Extrapolated)Log->PrintfWarning("Inlet/outlet %d: VelocityMode==Extrapolated is not recommended for inlet (or reverse) flows.",IdZone);
  if(ZsurfMode==ZSURF_Variable && RefillingMode==TFI_SimpleZsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleZsurf is invalid with ZSurface==Variable.",IdZone));
  if(ZsurfMode==ZSURF_Variable && RefillingMode==TFI_SimpleFull && RemoveZsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleFull is invalid with ZSurface==Variable and RemoveZsurf==true.",IdZone));
  if(ZsurfMode==ZSURF_Calculated && RefillingMode==TFI_SimpleZsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleZsurf is invalid with ZSurface==Calculated.",IdZone));
  if(ZsurfMode==ZSURF_Calculated && RefillingMode==TFI_SimpleFull && RemoveZsurf)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: RefillingMode==SimpleFull is invalid with ZSurface==Calculated and RemoveZsurf==true.",IdZone));
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
  for(unsigned c=1;c<unsigned(Simulate2D? 4: 8);c++){
    const tdouble3 ps=PtDom[c];
    pmin=MinValues(pmin,ps);
    pmax=MaxValues(pmax,ps);
  }
  BoxLimitMin=ToTFloat3(pmin);
  BoxLimitMax=ToTFloat3(pmax);
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
      case 0:  InputMode=TIN_Free;     break;
      case 1:  InputMode=TIN_Convert;  break;
      case 2:  InputMode=TIN_Remove;   break;
      default: sxml->ErrReadElement(ele,"inputtreatment",false);
    }
  }
  //-Loads RefillinMode.
  {
    const unsigned imode=sxml->ReadElementUnsigned(ele,"refilling","value");
    switch(imode){
      case 0:  RefillingMode=TFI_SimpleFull;   break;
      case 1:  RefillingMode=TFI_SimpleZsurf;  break;
      case 2:  RefillingMode=TFI_Advanced;     break;
      default: sxml->ErrReadElement(ele,"refilling",false);
    }
  }
  //-Loads layers value.
  unsigned layers=sxml->ReadElementUnsigned(ele,"layers","value");
  if(layers<1)sxml->ErrReadElement(ele,"layers",false,"Minumum number of layers is 1.");
  if(layers>250)sxml->ErrReadElement(ele,"layers",false,"Maximum number of layers is 250.");
  Layers=byte(layers);
  //-Creates inlet points.
  Points=new JSphInOutPoints(Log,Simulate2D,Simulate2DPosY,Layers,Dp,0,MapRealPosMin,MapRealPosMax);
  if(Simulate2D){
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
  TiXmlElement* xele=ele->FirstChildElement("imposevelocity");
  if(xele){
    const unsigned mode=sxml->GetAttributeUint(xele,"mode",true);
    switch(mode){
      case 0:  VelMode=MVEL_Fixed;         break;
      case 1:  VelMode=MVEL_Variable;      break;
      case 2:  VelMode=MVEL_Extrapolated;  break;
      case 3:  VelMode=MVEL_Interpolated;  break;
      default: sxml->ErrReadAtrib(xele,"mode",false,"Value is not valid.");
    }
    if(VelMode==MVEL_Fixed){
      const byte vel1=(sxml->ExistsElement(xele,"velocity" )? 1: 0);
      const byte vel2=(sxml->ExistsElement(xele,"velocity2")? 1: 0);
      const byte vel3=(sxml->ExistsElement(xele,"velocity3")? 1: 0);
      if(vel1+vel2+vel3>1)sxml->ErrReadElement(xele,"velocity",false,"Several definitions for velocity were found.");
      if(vel1 || vel1+vel2+vel3==0){
        VelProfile=PVEL_Constant;
        sxml->CheckAttributeNames(xele,"velocity","time v comment units_comment");
        InputVel=sxml->ReadElementFloat(xele,"velocity","v");
      }
      if(vel2){
        VelProfile=PVEL_Linear;
        sxml->CheckAttributeNames(xele,"velocity2","time v v2 z z2 comment units_comment");
        InputVel=sxml->ReadElementFloat(xele,"velocity2","v");
        InputVel2=sxml->ReadElementFloat(xele,"velocity2","v2");
        InputVelPosz=sxml->ReadElementFloat(xele,"velocity2","z");
        InputVelPosz2=sxml->ReadElementFloat(xele,"velocity2","z2");
      }
      if(vel3){
        VelProfile=PVEL_Parabolic;
        sxml->CheckAttributeNames(xele,"velocity3","time v v2 v3 z z2 z3 comment units_comment");
        InputVel =sxml->ReadElementFloat(xele,"velocity3","v");
        InputVel2=sxml->ReadElementFloat(xele,"velocity3","v2");
        InputVel3=sxml->ReadElementFloat(xele,"velocity3","v3");
        InputVelPosz =sxml->ReadElementFloat(xele,"velocity3","z");
        InputVelPosz2=sxml->ReadElementFloat(xele,"velocity3","z2");
        InputVelPosz3=sxml->ReadElementFloat(xele,"velocity3","z3");
      }
      sxml->CheckElementNames(xele,true,"*velocity *velocity2 *velocity3");
    }
    else if(VelMode==MVEL_Variable){
      const unsigned NTMIN=50;
      const byte vel1= (sxml->ExistsElement(xele,"velocitytimes" )? 1: 0);
      const byte vel2= (sxml->ExistsElement(xele,"velocitytimes2")? 1: 0);
      const byte vel3= (sxml->ExistsElement(xele,"velocitytimes3")? 1: 0);
      const byte vel1f=(sxml->ExistsElement(xele,"velocityfile"  )? 1: 0);
      const byte vel2f=(sxml->ExistsElement(xele,"velocityfile2" )? 1: 0);
      const byte vel3f=(sxml->ExistsElement(xele,"velocityfile3" )? 1: 0);
      if(vel1+vel2+vel3+vel1f+vel2f+vel3f>1)sxml->ErrReadElement(xele,"velocitytimes/velocityfile",false,"Several definitions for velocity were found.");
      if(vel1+vel2+vel3+vel1f+vel2f+vel3f==0)sxml->ErrReadElement(xele,"velocitytimes/velocityfile",false,"No definitions for variable velocity were found.");
      if(vel1 || vel1f){
        VelProfile=PVEL_Constant;
        TimeInputVelData=new JLinearValue(1);
        TiXmlElement* xlis=xele->FirstChildElement("velocitytimes");
        if(xlis){
          TimeInputVelData->SetSize(NTMIN);
          TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
          while(elet){
            sxml->CheckAttributeNames(elet,"time v");
            double t=sxml->GetAttributeDouble(elet,"time");
            double v=sxml->GetAttributeDouble(elet,"v");
            TimeInputVelData->AddTimeValue(t,v);
            elet=elet->NextSiblingElement("timevalue");
          }
        }
        else{
          TimeInputVelData->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"velocityfile","file"));
        }
      }
      if(vel2 || vel2f){
        VelProfile=PVEL_Linear;
        TimeInputVelData=new JLinearValue(4);
        TiXmlElement* xlis=xele->FirstChildElement("velocitytimes2");
        if(xlis){
          TimeInputVelData->SetSize(NTMIN);
          TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
          while(elet){
            sxml->CheckAttributeNames(elet,"time v v2 z z2");
            double t=sxml->GetAttributeDouble(elet,"time");
            double v=sxml->GetAttributeDouble(elet,"v");
            double v2=sxml->GetAttributeDouble(elet,"v2");
            double z=sxml->GetAttributeDouble(elet,"z");
            double z2=sxml->GetAttributeDouble(elet,"z2");
            TimeInputVelData->AddTimeValue(t,v,v2,z,z2);
            elet=elet->NextSiblingElement("timevalue");
          }
        }
        else{
          TimeInputVelData->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"velocityfile2","file"));
        }
      }
      if(vel3 || vel3f){
        VelProfile=PVEL_Parabolic;
        TimeInputVelData=new JLinearValue(6);
        TiXmlElement* xlis=xele->FirstChildElement("velocitytimes3");
        if(xlis){
          TimeInputVelData->SetSize(NTMIN);
          TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
          while(elet){
            sxml->CheckAttributeNames(elet,"time v v2 v3 z z2 z3");
            double t=sxml->GetAttributeDouble(elet,"time");
            double v=sxml->GetAttributeDouble(elet,"v");
            double v2=sxml->GetAttributeDouble(elet,"v2");
            double v3=sxml->GetAttributeDouble(elet,"v3");
            double z=sxml->GetAttributeDouble(elet,"z");
            double z2=sxml->GetAttributeDouble(elet,"z2");
            double z3=sxml->GetAttributeDouble(elet,"z3");
            TimeInputVelData->AddTimeValue(t,v,v2,v3,z,z2,z3);
            elet=elet->NextSiblingElement("timevalue");
          }
        }
        else{
          TimeInputVelData->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"velocityfile3","file"));
        }
      }
      if(!TimeInputVelData || TimeInputVelData->GetCount()<1)Run_Exceptioon("There are not inlet/outlet velocity values.");
      sxml->CheckElementNames(xele,true,"*velocitytimes *velocitytimes2 *velocitytimes3 velocityfile velocityfile2 velocityfile3");
    }
    else if(VelMode==MVEL_Interpolated){
      string checklist="gridveldata gridresetvelz gridposzero";
      sxml->CheckElementNames(xele,true,checklist);
      InputVelGrid=new JSphInOutGridData(Log);
      InputVelGrid->ConfigFromFile(dirdatafile+sxml->ReadElementStr(xele,"gridveldata","file"));
      ResetZVelGrid=sxml->ReadElementBool(xele,"gridresetvelz","value",true,false);
      double xmin=sxml->ReadElementDouble(xele,"gridposzero","x",true,0);
      double zmin=sxml->ReadElementDouble(xele,"gridposzero","z",true,0);
      InputVelGrid->SetPosMin(TDouble3(xmin,(Simulate2D? Simulate2DPosY: MapRealPosMin.y),zmin));
    }
    else if(VelMode!=MVEL_Extrapolated)Run_Exceptioon("Inlet/outlet velocity profile is unknown.");
  }
  else sxml->ErrReadElement(ele,"imposevelocity",true);

  //-Z-surface configuration.
  xele=ele->FirstChildElement("imposezsurf");
  if(xele){
    const unsigned mode=sxml->GetAttributeUint(xele,"mode",true);
    switch(mode){
      case 0:  
        ZsurfMode=ZSURF_Fixed;       
        sxml->CheckElementNames(xele,true,"zbottom zsurf remove");
      break;
      case 1:  
        ZsurfMode=ZSURF_Variable;    
        sxml->CheckElementNames(xele,true,"zbottom savevtk remove zsurftimes zsurffile");
      break;
      case 2:  
        ZsurfMode=ZSURF_Calculated;  
        sxml->CheckElementNames(xele,true,"zbottom zsurf savevtk remove");
      break;
      case 3:  
        ZsurfMode=ZSURF_WaveTheory;  
      break;
      default: sxml->ErrReadAtrib(xele,"mode",false,"Value is not valid.");
    }
    InputZbottom=sxml->ReadElementFloat(xele,"zbottom","value");
    if(ZsurfMode==ZSURF_Fixed || ZsurfMode==ZSURF_Calculated)InputZsurf=sxml->ReadElementFloat(xele,"zsurf","value");
    else if(ZsurfMode==ZSURF_Variable){
      const byte zlist=(sxml->ExistsElement(xele,"zsurftimes")? 1: 0);
      const byte zfile=(sxml->ExistsElement(xele,"zsurffile" )? 1: 0);
      if(zlist+zfile>1)sxml->ErrReadElement(xele,"zsurftimes/zsurffile",false,"Several definitions for zsurf were found.");
      if(zlist+zfile==0)sxml->ErrReadElement(xele,"zsurftimes/zsurffile",false,"No definitions for variable zsurf were found.");
      //-Loads list of time-values.
      TimeInputZsurf=new JLinearValue(1);
      TiXmlElement* xlis=xele->FirstChildElement("zsurftimes");
      if(xlis){
        const unsigned NTMIN=50;
        TimeInputZsurf->SetSize(NTMIN);
        TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
        while(elet){
          double t=sxml->GetAttributeDouble(elet,"time");
          double v=sxml->GetAttributeDouble(elet,"zsurf");
          TimeInputZsurf->AddTimeValue(t,v);
          elet=elet->NextSiblingElement("timevalue");
        }
      }
      else{
        TimeInputZsurf->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"zsurffile","file"));
      }
      //-Select firt zsurf value.
      InputZsurf=(float)TimeInputZsurf->GetValue(0);
    }
    RemoveZsurf=sxml->ReadElementBool(xele,"remove","value",true,false);
    if(ZsurfMode==ZSURF_Calculated || ZsurfMode==ZSURF_Variable)SvVtkZsurf=sxml->ReadElementBool(xele,"savevtk","value",true,false);
  }
  else if(ZsurfMode!=ZSURF_WaveTheory)ZsurfMode=ZSURF_Undefined;

  //-Rhop configuration.
  xele=ele->FirstChildElement("imposerhop");
  if(xele){
    const unsigned mode=sxml->GetAttributeUint(xele,"mode",true);
    switch(mode){
      case 0:  RhopMode=MRHOP_Constant;      break;
      case 1:  RhopMode=MRHOP_Hydrostatic;   break;
      case 2:  RhopMode=MRHOP_Extrapolated;  break;
      default: sxml->ErrReadAtrib(xele,"mode",false,"Value is not valid.");
    }
    //-Checks old configuration.
    if(sxml->ExistsElement(ele,"zsurf"  ))Run_ExceptioonFile(fun::PrintStr("Inlet/outlet zone %d: <zsurf> is not supported by current inlet/outlet version.",IdZone),sxml->ErrGetFileRow(ele,"zsurf"));
    if(sxml->ExistsElement(ele,"zbottom"))Run_ExceptioonFile(fun::PrintStr("Inlet/outlet zone %d: <zbottom> is not supported by current inlet/outlet version.",IdZone),sxml->ErrGetFileRow(ele,"zbottom"));
    //-Checks configuration.
    sxml->CheckElementNames(xele,true," ");
  }
  else{//-Default configuration for rhop.
    RhopMode=MRHOP_Constant;
  }

  //-Calcule initial number of new particles.
  if(ZsurfMode==ZSURF_Undefined || RefillingMode==TFI_SimpleFull)NpartInit=NptInit*Layers;
  else NpartInit=Points->GetCountZmax(InputZsurf)*Layers;

  //-Loads points in fluid domain to calculate zsurf.
  if(ZsurfMode==ZSURF_Calculated){
    PtzPos=new tfloat3[NptInit];
    const tdouble3 *ptpos=Points->GetPoints();
    const double dist=Dp/2 + GetDistPtzPos();//-It adds Dp/2 because the inlet limit is at +Dp/2.
    const tdouble3 dirmov=Direction*dist;
    for(unsigned p=0;p<NptInit;p++)PtzPos[p]=ToTFloat3(ptpos[p]+dirmov);
    #ifdef _WITHGPU
    if(!Cpu){
      fcuda::Malloc(&PtzPosg,NptInit);
      cudaMemcpy(PtzPosg,PtzPos,sizeof(float3)*NptInit,cudaMemcpyHostToDevice);
      fcuda::Malloc(&PtzAuxg,NptInit);
      PtzAux=new float[NptInit];
    }
    #endif
  }
  InputCheck=(RemoveZsurf || InputMode!=TIN_Free);

  //-Calculate velocity limits.
  CalculateVelMinMax(VelMin,VelMax);
}

//==============================================================================
/// Calculates minimum and maximum velocity according velocity configuration.
/// Returns -FLT_MAX and FLT_MAX when velocity is unknown.
//==============================================================================
void JSphInOutZone::CalculateVelMinMax(float &velmin,float &velmax)const{
  velmin=FLT_MAX;
  velmax=-FLT_MAX;
  if(VelMode==MVEL_Fixed){
    switch(VelProfile){
      case PVEL_Constant:
        velmin=velmax=InputVel;  
      break;
      case PVEL_Linear:    
        velmin=min(InputVel,InputVel2);
        velmax=max(InputVel,InputVel2);  
      break;
      case PVEL_Parabolic:    
        velmin=min(InputVel,min(InputVel2,InputVel3));
        velmax=max(InputVel,max(InputVel2,InputVel3));
      break;
      default: Run_Exceptioon("Velocity profile is unknown.");
    }
  }
  else if(VelMode==MVEL_Variable){
    const unsigned count=TimeInputVelData->GetCount();
    switch(VelProfile){
      case PVEL_Constant:
        for(unsigned c=0;c<count;c++){
          const float v=float(TimeInputVelData->GetValueByIdx(c));
          velmin=min(velmin,v);
          velmax=max(velmax,v);
        }
      break;
      case PVEL_Linear:    
        for(unsigned c=0;c<count;c++){
          const float v0=float(TimeInputVelData->GetValueByIdx(c,0));
          const float v1=float(TimeInputVelData->GetValueByIdx(c,1));
          velmin=min(velmin,min(v0,v1));
          velmax=max(velmax,max(v0,v1));
        }
      break;
      case PVEL_Parabolic:    
        for(unsigned c=0;c<count;c++){
          const float v0=float(TimeInputVelData->GetValueByIdx(c,0));
          const float v1=float(TimeInputVelData->GetValueByIdx(c,1));
          const float v2=float(TimeInputVelData->GetValueByIdx(c,2));
          velmin=min(velmin,min(v0,min(v1,v2)));
          velmax=max(velmax,max(v0,max(v1,v2)));
        }
      break;
      default: Run_Exceptioon("Velocity profile is unknown.");
    }
  }
  else if(VelMode==MVEL_Extrapolated || VelMode==MVEL_Interpolated){
    velmin=-FLT_MAX;
    velmax=FLT_MAX;
  }
  else Run_Exceptioon("Velocity mode is unknown.");
}

//==============================================================================
/// Returns a byte with information about VelMode, VelProfile and RhopMode.
//==============================================================================
byte JSphInOutZone::GetConfigZone()const{
  if((VelProfile&PVEL_MASK )!=VelProfile)Run_Exceptioon("VelProfile value is invalid to code using mask.");
  if((VelMode   &MVEL_MASK )!=VelMode   )Run_Exceptioon("VelMode value is invalid to code using mask.");
  if((RhopMode  &MRHOP_MASK)!=RhopMode  )Run_Exceptioon("RhopMode value is invalid to code using mask.");
  //if(byte(ZsurfMode) >3)Run_Exceptioon("ZsurfMode value is invalid to code in 2 bits.");
  byte ret=byte(VelProfile);
  ret|=byte(VelMode);
  ret|=byte(RhopMode);
  if(InputCheck)ret|=byte(CheckInput_MASK);
  //-Checks coded configuration.
  TpVelProfile vprof =GetConfigVelProfile(ret);
  TpVelMode    vmode =GetConfigVelMode(ret);
  TpRhopMode   rmode =GetConfigRhopMode(ret);
  bool         checkf=GetConfigCheckInputDG(ret);
  //TpZsurfMode  zmode=GetConfigZsurfMode(ret);
  //if(vprof!=VelProfile || vmode!=VelMode || rmode!=RhopMode || zmode!=ZsurfMode)Run_Exceptioon("Coded configuration is not right.");
  if(vprof!=VelProfile || vmode!=VelMode || rmode!=RhopMode || checkf!=InputCheck)Run_Exceptioon("Coded configuration is not right.");
  return(ret);
}

//==============================================================================
/// Returns a byte with information about refilling mode and RemoveZsurf.
//==============================================================================
byte JSphInOutZone::GetConfigUpdate()const{
  byte ret=byte(RefillingMode==TFI_Advanced  ? RefillAdvanced_MASK: 0);
  ret= ret|byte(RefillingMode==TFI_SimpleFull? RefillSpFull_MASK  : 0);
  ret= ret|byte(InputMode==TIN_Remove ? RemoveInput_MASK : 0);
  ret= ret|byte(GetRemoveZsurf()      ? RemoveZsurf_MASK : 0);
  ret= ret|byte(InputMode==TIN_Convert? ConvertInput_MASK: 0);
  return(ret);
}

//==============================================================================
/// Loads positon of inlet/outlet points.
//==============================================================================
unsigned JSphInOutZone::LoadInletPoints(tdouble3 *pos){
  const tdouble3 *ptpos=Points->GetPoints();
  const tdouble3 dir=Direction*(Dp/2);
  for(unsigned cp=0;cp<NptInit;cp++)pos[cp]=ptpos[cp]+dir;
  return(NptInit);
}

//==============================================================================
/// Loads positon of inlet/outlet particles.
//==============================================================================
void JSphInOutZone::LoadInletParticles(tdouble3 *pos){
  const tdouble3 *ptpos=Points->GetPoints();
  unsigned p=0;
  for(unsigned layer=0;layer<Layers;layer++){
    const tdouble3 dir=Direction*(-Dp*layer);
    for(unsigned cp=0;cp<NptInit;cp++)if(RefillingMode==TFI_SimpleFull || ZsurfMode==ZSURF_Undefined || float(ptpos[cp].z)<=InputZsurf){
      pos[p]=ptpos[cp]+dir;
      p++;
    }
  }
  //Log->Printf(" LoadInletParticles--> p:%u  NptInit:%u  Layers:%u",p,NptInit,Layers);
  Points->ResetPoints();//-Frees memory of points.
}

//==============================================================================
/// Returns velocity according profile configuration.
//==============================================================================
float JSphInOutZone::CalcVel(JSphInOutZone::TpVelProfile vprof,const tfloat4 &vdata,double posz){
  float vel=0;
  if(vprof==JSphInOutZone::PVEL_Constant)vel=vdata.x;
  else if(vprof==JSphInOutZone::PVEL_Linear){
    const float m=vdata.x;
    const float b=vdata.y;
    vel=m*float(posz)+b;
  }
  else if(vprof==JSphInOutZone::PVEL_Parabolic){
    const float a=vdata.x;
    const float b=vdata.y;
    const float c=vdata.z;
    vel=a*float(posz)*float(posz)+b*float(posz)+c;
  }
  return(vel);
}

//==============================================================================
/// Updates velocity coefficients for imposed velocity according timestep. 
/// Returns true when the data was changed.
//==============================================================================
bool JSphInOutZone::UpdateVelData(double timestep,bool full,tfloat4 &vel0,tfloat4 &vel1){
  bool modified=full;
  if(VelMode==MVEL_Fixed){
    if(full){
      if(VelProfile==PVEL_Constant){
        vel0=TFloat4(InputVel,0,0,0);  vel1=TFloat4(0);
        //const float v=InputVel;
      }
      else if(VelProfile==PVEL_Linear){
        const float m=(InputVel2-InputVel)/(InputVelPosz2-InputVelPosz);
        const float b=InputVel-m*InputVelPosz;
        vel0=TFloat4(m,b,0,0);  vel1=TFloat4(0);
        //const float v=m*float(vecpos[p].z)+b;
        //if(SaveVelProfile && Npt){
        //  SaveVelProfile=false;
        //  JSaveCsv sv(Log->GetDirOut()+fun::PrintStr("Inlet_%d_VelProfileLinear.csv",IdZone),true);
        //  sv.AddHead("time;z;vel");
        //  const double y0=vecpos[0].y;
        //  for(unsigned p=0;p<Npt;p++)if(y0==vecpos[p].y){
        //    sv.AddValuesf("%g;%g;%f",timestep,vecpos[p].z,vecvel[p]);
        //    sv.AddEndl();
        //  }
        //}
      }
      else if(VelProfile==PVEL_Parabolic){
        const tmatrix3f inv=fmath::InverseMatrix3x3(TMatrix3f(InputVelPosz*InputVelPosz,InputVelPosz,1,InputVelPosz2*InputVelPosz2,InputVelPosz2,1,InputVelPosz3*InputVelPosz3,InputVelPosz3,1));
        const float a=inv.a11*InputVel+inv.a12*InputVel2+inv.a13*InputVel3;
        const float b=inv.a21*InputVel+inv.a22*InputVel2+inv.a23*InputVel3;
        const float c=inv.a31*InputVel+inv.a32*InputVel2+inv.a33*InputVel3;
        vel0=TFloat4(a,b,c,0);  vel1=TFloat4(0);
        //const float v=a*float(vecpos[p].z)*float(vecpos[p].z)+b*float(vecpos[p].z)+c;
        //if(SaveVelProfile && Npt){
        //  SaveVelProfile=false;
        //  JSaveCsv sv(Log->GetDirOut()+fun::PrintStr("Inlet_%d_VelProfileParabolic.csv",IdZone),true);
        //  sv.AddHead("time;z;vel");
        //  const double y0=vecpos[0].y;
        //  for(unsigned p=0;p<Npt;p++)if(y0==vecpos[p].y){
        //    sv.AddValuesf("%g;%g;%f",timestep,vecpos[p].z,vecvel[p]);
        //    sv.AddEndl();
        //  }
        //}
      }
      else Run_Exceptioon("Inlet/outlet velocity profile is unknown.");
    }
  }
  else if(VelMode==MVEL_Variable){
    TimeInputVelData->FindTime(timestep);
    unsigned idx0=TimeVelIdx0,idx1=TimeVelIdx1;
    TimeVelIdx0=TimeInputVelData->GetPos();
    TimeVelIdx1=TimeInputVelData->GetPosNext();
    if(full || idx0!=TimeVelIdx0 || idx1!=TimeVelIdx1){
      modified=true;
      const float t =(float)TimeInputVelData->GetTimeByIdx(TimeVelIdx0);
      const float t2=(float)TimeInputVelData->GetTimeByIdx(TimeVelIdx1);
      if(VelProfile==PVEL_Constant){
        //TimeInputVelData->AddTimeValue(t,v);
        vel0=TFloat4((float)TimeInputVelData->GetValueByIdx(TimeVelIdx0),0,0,t);
        vel1=TFloat4((float)TimeInputVelData->GetValueByIdx(TimeVelIdx1),0,0,t2);
      }
      else if(VelProfile==PVEL_Linear){
        //TimeInputVelData->AddTimeValue(t,v,v2,z,z2);
        for(unsigned cc=0;cc<=1;cc++){
          const unsigned idx=(!cc? TimeVelIdx0: TimeVelIdx1);
          const float v =(float)TimeInputVelData->GetValueByIdx(idx,0);
          const float v2=(float)TimeInputVelData->GetValueByIdx(idx,1);
          const float z =(float)TimeInputVelData->GetValueByIdx(idx,2);
          const float z2=(float)TimeInputVelData->GetValueByIdx(idx,3);
          const float m=(v2-v)/(z2-z);
          const float b=v-m*z;
          if(!cc)vel0=TFloat4(m,b,0,t);
          else   vel1=TFloat4(m,b,0,t2);
        }
      }
      else if(VelProfile==PVEL_Parabolic){
        //TimeInputVelData->AddTimeValue(t,v,v2,v3,z,z2,z3);
        for(unsigned cc=0;cc<=1;cc++){
          const unsigned idx=(!cc? TimeVelIdx0: TimeVelIdx1);
          const float v =(float)TimeInputVelData->GetValueByIdx(idx,0);
          const float v2=(float)TimeInputVelData->GetValueByIdx(idx,1);
          const float v3=(float)TimeInputVelData->GetValueByIdx(idx,2);
          const float z =(float)TimeInputVelData->GetValueByIdx(idx,3);
          const float z2=(float)TimeInputVelData->GetValueByIdx(idx,4);
          const float z3=(float)TimeInputVelData->GetValueByIdx(idx,5);
          const tmatrix3f inv=fmath::InverseMatrix3x3(TMatrix3f(z*z,z,1,z2*z2,z2,1,z3*z3,z3,1));
          const float a=inv.a11*v+inv.a12*v2+inv.a13*v3;
          const float b=inv.a21*v+inv.a22*v2+inv.a23*v3;
          const float c=inv.a31*v+inv.a32*v2+inv.a33*v3;
          if(!cc)vel0=TFloat4(a,b,c,t);
          else   vel1=TFloat4(a,b,c,t2);
        }
      }
      else Run_Exceptioon("Inlet/outlet velocity profile is unknown.");
    }
  }
  else if(VelMode!=MVEL_Extrapolated && VelMode!=MVEL_Interpolated)Run_Exceptioon("Inlet/outlet velocity mode is unknown.");
  return(modified);
}

//==============================================================================
/// Updates zsurf according timestep. 
/// Returns true when the data was changed.
//==============================================================================
bool JSphInOutZone::UpdateZsurf(double timestep,bool full,float &zsurf){
  bool modified=full;
  if(ZsurfMode==ZSURF_Undefined || ZsurfMode==ZSURF_Fixed){
    if(full)zsurf=InputZsurf;
  }
  else if(ZsurfMode==ZSURF_Variable){
    const float zsurf1=(float)TimeInputZsurf->GetValue(timestep);
    //Log->Printf(" JSphInOutZone::UpdateZsurf---> zsurf1:%f  time:%f",zsurf1,timestep);
    modified|=(InputZsurf!=zsurf1);
    if(modified)zsurf=InputZsurf=zsurf1;
  }
  else if(ZsurfMode==ZSURF_Calculated){
    modified|=(InputZsurf!=zsurf);
    if(modified)zsurf=InputZsurf;
  }
  else Run_Exceptioon("Inlet/outlet zsurf mode is unknown.");
  return(modified);
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

//==============================================================================
/// Saves grid nodes in VTK file.
//==============================================================================
void JSphInOutZone::SaveVtkVelGrid(std::string filename)const{
  if(InputVelGrid)InputVelGrid->SaveDataVtk(filename,0);
}


//==============================================================================
/// Returns true when zone is MVEL_Variable, PVEL_Constant and ZSURF_Variable.
//==============================================================================
bool JSphInOutZone::IsVariableInput()const{
  return(VelMode==MVEL_Variable && VelProfile==PVEL_Constant && ZsurfMode==ZSURF_Variable);
}

//==============================================================================
/// Activate or deactivate use of external data (ExternalVarInput).
//==============================================================================
void JSphInOutZone::SetExternalVarInput(bool active){
  if(!IsVariableInput())Run_Exceptioon("Input configuration is invalid.");
  if(!ExternalVarInput && active)InitExtVarInput();
  ExternalVarInput=active;
}

//==============================================================================
/// Initialisation of TimeInputVelData and TimeInputZsurf.
//==============================================================================
void JSphInOutZone::InitExtVarInput(){
  if(!IsVariableInput())Run_Exceptioon("Input configuration is invalid.");
  TimeInputVelData->Reset();
  TimeInputVelData->AddTimeValue(0,0);
  TimeInputVelData->AddTimeValue(DBL_MAX,0);
  TimeInputZsurf->Reset();
  TimeInputZsurf->AddTimeValue(0,0);
  TimeInputZsurf->AddTimeValue(DBL_MAX,0);
}

//==============================================================================
/// Returns velocity configuration in TimeInputVelData.
/// Returns TFloat4(time1,vel1,time2,vel2).
//==============================================================================
tfloat4 JSphInOutZone::GetExtVarVelocity1()const{
  if(VelMode!=MVEL_Variable || VelProfile!=PVEL_Constant)Run_Exceptioon("Velocity mode or profile is invalid.");
  const unsigned nt=TimeInputVelData->GetCount();
  unsigned idx=0;
  if(nt>2)idx=TimeInputVelData->GetPos();
  const float t1=(idx<nt?   float(TimeInputVelData->GetTimeByIdx (idx))    : 0);
  const float v1=(idx<nt?   float(TimeInputVelData->GetValueByIdx(idx,0))  : 0);
  const float t2=(idx+1<nt? float(TimeInputVelData->GetTimeByIdx (idx+1))  : t1);
  const float v2=(idx+1<nt? float(TimeInputVelData->GetValueByIdx(idx+1,0)): v1);
  return(TFloat4(t1,v1,t2,v2));
}

//==============================================================================
/// Returns zsurf configuration in TimeInputZsurf.
/// Returns TFloat4(time1,z1,time2,z2).
//==============================================================================
tfloat4 JSphInOutZone::GetExtVarZsurf()const{
  if(ZsurfMode!=ZSURF_Variable)Run_Exceptioon("Zsurf mode is not Variable.");
  const unsigned nt=TimeInputZsurf->GetCount();
  unsigned idx=0;
  if(nt>2)idx=TimeInputZsurf->GetPos();
  const float t1=(idx<nt?   float(TimeInputZsurf->GetTimeByIdx (idx))    : 0);
  const float v1=(idx<nt?   float(TimeInputZsurf->GetValueByIdx(idx,0))  : 0);
  const float t2=(idx+1<nt? float(TimeInputZsurf->GetTimeByIdx (idx+1))  : t1);
  const float v2=(idx+1<nt? float(TimeInputZsurf->GetValueByIdx(idx+1,0)): v1);
  return(TFloat4(t1,v1,t2,v2));
}

//==============================================================================
/// Modify velocity configuration in TimeInputVelData.
//==============================================================================
void JSphInOutZone::SetExtVarVelocity1(double t1,double v1,double t2,double v2){
  if(!ExternalVarInput)Run_Exceptioon("ExternalVarInput is not activated.");
  TimeInputVelData->SetTimeValue(0,t1,v1);
  TimeInputVelData->SetTimeValue(1,t2,v2);
}

//==============================================================================
/// Modify zsurf configuration in TimeInputZsurf.
//==============================================================================
void JSphInOutZone::SetExtVarZsurf(double t1,double v1,double t2,double v2){
  if(!ExternalVarInput)Run_Exceptioon("ExternalVarInput is not activated.");
  TimeInputZsurf->SetTimeValue(0,t1,v1);
  TimeInputZsurf->SetTimeValue(1,t2,v2);
}


