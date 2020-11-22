//HEAD_DSPH
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

/// \file JSphInOutVelAwas.cpp \brief Implements the class \ref JSphInOutVelAwas.

#include "JSphInOutVelAwas.h"
#include "Functions.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JXml.h"
#include "JLinearValue.h"
#include "JDsGaugeSystem.h"

#include <cstring>
#include <cfloat>
#include <cmath>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSphInOutVelAwas
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutVelAwas::JSphInOutVelAwas(unsigned idzone,double inletx,tdouble3 inletdir
  ,float gravityz,const std::string &dirdatafile,JGaugeSystem *gaugesystem
  ,const JXml *sxml,TiXmlElement* ele)
  :Log(AppInfo.LogPtr()),IdZone(idzone),InletX(inletx),InletDir(inletdir),GravityZ(gravityz)
{
  ClassName="JSphInOutVelAwas";
  ZsurfTarget=NULL;
  StepData=NULL;
  Reset();
  ReadXml(sxml,ele,dirdatafile,gaugesystem);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOutVelAwas::~JSphInOutVelAwas(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInOutVelAwas::Reset(){
  //-Initial configuration variables.
  InletMode=false;
  StartAwas=InitDepth=CoefDepth=0;
  delete ZsurfTarget;  ZsurfTarget=NULL;
  GaugeX=GaugeXh=GaugeXdp=GaugeY=GaugeZmin=GaugeZmax=GaugeDpXml=GaugeDp=0;
  SaveData=0;
  GaugeSwl=NULL;
  //-Saves step data for CSV.
  CountStepData=0;
  delete[] StepData; StepData=NULL;
  //-Values for the current step.
  LastTimeStep=0;
  LastZgauge=LastZtarget=LastVelCorr=0;
}

//==============================================================================
/// Reads initial configuration in the XML node.
//==============================================================================
void JSphInOutVelAwas::ReadXml(const JXml *sxml,TiXmlElement* ele
  ,const std::string &dirdatafile,JGaugeSystem *gaugesystem)
{
  //-Checks element names.
  sxml->CheckElementNames(ele,true,"inletmode startawas depth zsurffile gaugex gaugey gaugezmin gaugezmax gaugedp savedata");
  //-Loads configuration.
  InletMode=sxml->ReadElementBool(ele,"inletmode","value");
  StartAwas=sxml->ReadElementDouble(ele,"startawas","value",true,0);
  InitDepth=sxml->ReadElementDouble(ele,"depth","value");
  CoefDepth=sqrt(-GravityZ/InitDepth);
  string filezsurf=sxml->ReadElementStr(ele,"zsurffile","file");
  ZsurfTarget=new JLinearValue(1);
  ZsurfTarget->LoadFile(dirdatafile+filezsurf);
  GaugeX=sxml->ReadElementDouble(ele,"gaugex","value",true,DBL_MAX);
  GaugeXh=sxml->ReadElementDouble(ele,"gaugex","valueh",true,DBL_MAX);
  GaugeXdp=sxml->ReadElementDouble(ele,"gaugex","valuedp",true,DBL_MAX);
  if((GaugeX!=DBL_MAX? 1: 0)+(GaugeXh!=DBL_MAX? 1: 0)+(GaugeXdp!=DBL_MAX? 1: 0)>1){
    sxml->ErrReadElement(sxml->GetFirstElement(ele,"gaugex",false),"gaugex",false,"GaugeX value must be defined only one time (using value, valueh or valuedp).");
  }
  GaugeY=sxml->ReadElementDouble(ele,"gaugey","value",gaugesystem->GetSimulate2D(),gaugesystem->GetSimulate2DPosY());
  GaugeZmin=sxml->ReadElementDouble(ele,"gaugezmin","value",true,-DBL_MAX);
  GaugeZmax=sxml->ReadElementDouble(ele,"gaugezmax","value",true,DBL_MAX);
  GaugeDpXml=sxml->ReadElementDouble(ele,"gaugedp","value",true,0.1);
  //-Save data configuration.
  SaveData=(byte)sxml->ReadElementUnsigned(ele,"savedata","value",true,0);
  if(SaveData>2)sxml->ErrReadElement(ele,"savedata",false,"Value out of range.");
  if(SaveData==2){
    CountStepData=0;
    StepData=new tfloat4[SizeStepData];
  }
  //-Defines gauge to compute water level.
  if(!InletDir.x || InletDir.y!=0 || InletDir.z!=0)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: Direction is invalid for AWAS. Only X-direction is valid.",IdZone));
  if(GaugeX==DBL_MAX){
    if(GaugeX==DBL_MAX)GaugeX=gaugesystem->GetDp()*5;
    if(GaugeXdp!=DBL_MAX)GaugeX=GaugeXdp*gaugesystem->GetDp();
    if(GaugeXh!=DBL_MAX)GaugeX=GaugeXh*gaugesystem->GetKernelH();
    GaugeX=InletX+(InletDir.x>0? GaugeX: -GaugeX);
  }
  GaugeDp=GaugeDpXml*gaugesystem->GetDp();
  if(GaugeZmin<gaugesystem->GetDomPosMin().z)GaugeZmin=gaugesystem->GetDomPosMin().z;
  if(GaugeZmax>gaugesystem->GetDomPosMax().z)GaugeZmax=gaugesystem->GetDomPosMax().z;
  unsigned gaugenmax=unsigned((GaugeZmax-GaugeZmin)/GaugeDp);
  GaugeZmax=GaugeZmin+GaugeDp*gaugenmax;
  const tdouble3 point0=TDouble3(GaugeX,GaugeY,GaugeZmin);
  const tdouble3 point2=TDouble3(GaugeX,GaugeY,GaugeZmax);
  GaugeSwl=gaugesystem->AddGaugeSwl(fun::PrintStr("AwasInlet%02u",IdZone),StartAwas,DBL_MAX,0,point0,point2,GaugeDp,0);
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphInOutVelAwas::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back("  AWAS configuration:");
  if(StartAwas)lines.push_back(fun::PrintStr("    StartAWAS.: %f [s]",StartAwas));
               lines.push_back(fun::PrintStr("    Depth.....: %g",InitDepth));
               lines.push_back(fun::PrintStr("    Zsurf file: %s",ZsurfTarget->GetFile().c_str()));
               lines.push_back(fun::PrintStr("    GaugePos..: (%g, %g, %g -- %g) GaugeDp:%g [m]",GaugeX,GaugeY,GaugeZmin,GaugeZmax,GaugeDp));
   if(SaveData)lines.push_back(fun::PrintStr("    SaveData..: %s",(SaveData==1? "PARTs": "Steps")));
}

//==============================================================================
/// Computes and returns AWAS correction for current timestep.
//==============================================================================
float JSphInOutVelAwas::GetVelCorr(double timestep){
  const bool active=GaugeSwl->GetResult().modified;
  //if(IdZone==1)Log->Printf("--> t:%g t:%g mod:%d",timestep,GaugeSwl->GetResult().timestep,(GaugeSwl->GetResult().modified?1:0));
  const double zgauge=GaugeSwl->GetResult().posswl.z;
  const double ztarget=ZsurfTarget->GetValue(timestep);
  LastVelCorr=(active? float((InletMode? zgauge-ztarget: ztarget-zgauge)*CoefDepth): 0);
  LastTimeStep=timestep;
  LastZgauge =float(zgauge);
  LastZtarget=float(ztarget);
  if(StepData){
    if(CountStepData>=SizeStepData)SaveCsvData();
    StepData[CountStepData++]=TFloat4(float(LastTimeStep),LastZgauge,LastZtarget,LastVelCorr);
  }
  return(LastVelCorr);
}

//==============================================================================
/// Saves CSV file with AWAS information.
//==============================================================================
void JSphInOutVelAwas::SaveCsvData(){
  if(SaveData==1 || (SaveData==2 && CountStepData)){
    const string file=AppInfo.GetDirOut()+fun::PrintStr("InOut_AWAS_i%d.csv",IdZone); 
    jcsv::JSaveCsv2 scsv(file,true,AppInfo.GetCsvSepComa());
    if(scsv.GetAppendMode())Log->AddFileInfo(AppInfo.GetDirOut()+"InOut_AWAS_i?.csv","Saves CSV with Inlet-AWAS information.");
    scsv.SetHead();
    scsv << "Time [s];Zgauge [m];ZsurfTarget [m]; VelCorr [m/s]" << jcsv::Endl();
    scsv.SetData();
    if(StepData)for(unsigned c=0;c<CountStepData;c++)scsv << StepData[c] << jcsv::Endl();
    else scsv << LastTimeStep << LastZgauge << LastZtarget << LastVelCorr << jcsv::Endl();
    scsv.SaveData(true);
  }
  CountStepData=0;
}

