//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSaveDt.cpp \brief Implements the class \ref JSaveDt.

#include "JSaveDt.h"
#include "JLog2.h"
#include "JXml.h"
#include "Functions.h"
#include "JSaveCsv2.h"
#include <cfloat>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

//##############################################################################
//# JSaveDt
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSaveDt::JSaveDt(JLog2* log):Log(log){
  ClassName="JSaveDt";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSaveDt::~JSaveDt(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSaveDt::Reset(){
  TimeStart=TimeFinish=TimeInterval=0;
  FullInfo=AllDt=false;
  Count=0;
  LastInterval=0;
  memset(&ValueNull,0,sizeof(StValue));
  LastDtf=LastDt1=LastDt2=ValueNull;
  LastAceMax=LastViscDtMax=LastVelMax=ValueNull;
  CountAllDts=0;
}

//==============================================================================
/// Configures object.
//==============================================================================
void JSaveDt::Config(JXml *sxml,const std::string &place,double timemax,double timeout){
  Reset();
  LoadXml(sxml,place);
  if(TimeFinish<=0)TimeFinish=DBL_MAX;
  if(TimeInterval<0)TimeInterval=timeout;
  SizeValuesSave=max(1u,min(GetSizeValues(),unsigned(timeout/TimeInterval)));
}


//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSaveDt::LoadXml(JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JSaveDt::ReadXml(JXml *sxml,TiXmlElement* ele){
  const char met[]="ReadXml";
  TimeStart=sxml->ReadElementFloat(ele,"start","value",true);
  TimeFinish=sxml->ReadElementFloat(ele,"finish","value",true,-1);
  TimeInterval=sxml->ReadElementFloat(ele,"interval","value",true,-1);
  AllDt=sxml->ReadElementBool(ele,"alldt","value",true,false);
  FullInfo=sxml->ReadElementBool(ele,"fullinfo","value",true,false);
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSaveDt::VisuConfig(std::string txhead,std::string txfoot){
  if(!txhead.empty())Log->Print(txhead);
  if(TimeFinish==DBL_MAX)Log->Printf("  Time    : (%f - END)",TimeStart);
  else Log->Printf("  Time    : (%f - %f)",TimeStart,TimeFinish);
  Log->Printf("  Interval: %f (group:%u)",TimeInterval,SizeValuesSave);
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Stores file buffer values.
/// Graba valores de buffer en fichero.
//==============================================================================
void JSaveDt::SaveFileValues(){
  const char met[]="SaveFileValues";
  const bool firstsv=FileDtInfo.empty();
  if(firstsv){
    FileDtInfo=Log->GetDirOut()+"DtInfo.csv";
    Log->AddFileInfo(FileDtInfo,"Saves statistical information about DT values.");
  }
  jcsv::JSaveCsv2 scsv(FileDtInfo,!firstsv,Log->GetCsvSepComa());
  //-Saves head.
  if(firstsv){
    scsv.SetHead();
    scsv << "Time;Values";
    scsv << "Dtf_mean;Dtf_min;Dtf_max";
    scsv << "Dt1_mean;Dt1_min;Dt1_max";
    scsv << "Dt2_mean;Dt2_min;Dt2_max";
    if(FullInfo){
      scsv << "AceMax_mean;AceMax_min;AceMax_max";
      scsv << "ViscDtMax_mean;ViscDtMax_min;ViscDtMax_max";
      scsv << "VelMax_mean;VelMax_min;VelMax_max";
    }
    scsv << jcsv::Endl();
  }
  //-Saves data.
  scsv.SetData();
  scsv << jcsv::Fmt(jcsv::TpDouble1,"%20.12E");
  for(unsigned c=0;c<Count;c++){
    StValue v;
    v=DtFinal[c];    scsv << v.tini << v.num << v.vmean << v.vmin << v.vmax;
    v=Dt1[c];        scsv << v.vmean << v.vmin << v.vmax;
    v=Dt2[c];        scsv << v.vmean << v.vmin << v.vmax;
    if(FullInfo){
      v=AceMax[c];     scsv << v.vmean << v.vmin << v.vmax;
      v=ViscDtMax[c];  scsv << v.vmean << v.vmin << v.vmax;
      v=VelMax[c];     scsv << v.vmean << v.vmin << v.vmax;
    }
    scsv << jcsv::Endl();
  }
  scsv.SaveData();
  Count=0;
}

//==============================================================================
/// Stores buffer values.
/// Graba valores de buffer.
//==============================================================================
void JSaveDt::SaveFileValuesEnd(){
  if(LastDtf.num)AddLastValues();
  if(Count)SaveFileValues();
}

//==============================================================================
/// Stores file buffer values.
/// Graba valores de buffer en fichero.
//==============================================================================
void JSaveDt::SaveFileAllDts(){
  const char met[]="SaveFileAllDts";
  const bool firstsv=FileDtAllInfo.empty();
  if(firstsv){
    FileDtAllInfo=Log->GetDirOut()+"DtAllInfo.csv";
    Log->AddFileInfo(FileDtAllInfo,"Saves DT values for each simulation step.");
  }
  jcsv::JSaveCsv2 scsv(FileDtAllInfo,!firstsv,Log->GetCsvSepComa());
  //-Saves head.
  if(firstsv){
    scsv.SetHead();
    scsv << "Time;Dtf" << jcsv::Endl();
  }
  //-Saves data.
  scsv.SetData();
  scsv << jcsv::Fmt(jcsv::TpDouble1,"%20.12E");
  for(unsigned c=0;c<CountAllDts;c++)scsv << AllDts[c].x << AllDts[c].y << jcsv::Endl();
  CountAllDts=0;
}

//==============================================================================
/// Saves indicated info for dt. If it matches with timestep.
/// Guarda info del dt inicado. Si coincide timestep lo sobre.
//==============================================================================
void JSaveDt::AddValueData(double timestep,double dt,StValue &value){
  if(!value.num){
    value.tini=timestep;
    value.vmean=value.vmin=value.vmax=dt;
    value.num=1;
  }
  else{
    if(value.vmin>dt)value.vmin=dt;
    if(value.vmax<dt)value.vmax=dt;
    value.vmean=(value.vmean*value.num+dt)/(value.num+1);
    value.num++;
  }
}

//==============================================================================
/// Saves buffer info.
/// Guarda info en buffer.
//==============================================================================
void JSaveDt::AddLastValues(){
  if(Count>=GetSizeValues())SaveFileValues();
  DtFinal[Count]=LastDtf;
  Dt1[Count]=LastDt1;
  Dt2[Count]=LastDt2;
  LastDtf=LastDt1=LastDt2=ValueNull;
  if(FullInfo){
    AceMax[Count]=LastAceMax;
    ViscDtMax[Count]=LastViscDtMax;
    VelMax[Count]=LastVelMax;
    LastAceMax=LastViscDtMax=LastVelMax=ValueNull;
  }
  Count++;
}

//==============================================================================
/// Saves indicated info for dt. If it matches with timestep.
/// Guarda info del dt inicado. Si coincide timestep lo sobre
//==============================================================================
void JSaveDt::AddValues(double timestep,double dtfinal,double dt1,double dt2,double acemax,double viscdtmax,double velmax){
  if(TimeStart<=timestep && timestep<=TimeFinish){
    unsigned interval=unsigned((timestep-TimeStart)/TimeInterval);
    if(LastInterval!=interval && LastDtf.num){
      AddLastValues();
      if(Count>=SizeValuesSave)SaveFileValues();
    }
    LastInterval=interval;
    AddValueData(timestep,dtfinal,LastDtf);
    AddValueData(timestep,dt1,LastDt1);
    AddValueData(timestep,dt2,LastDt2);
    if(FullInfo){
      AddValueData(timestep,acemax,LastAceMax);
      AddValueData(timestep,viscdtmax,LastViscDtMax);
      AddValueData(timestep,velmax,LastVelMax);
    }
    //-Management of AllDt.
    //-Gestion de AllDt.
    if(AllDt){
      if(CountAllDts>=SizeAllDts)SaveFileAllDts();
      AllDts[CountAllDts]=TDouble2(timestep,dtfinal);
      CountAllDts++;
    }
  }
  else if(timestep>TimeFinish && Count)SaveFileValuesEnd();
}

//==============================================================================
/// Saves current data in output files.
//==============================================================================
void JSaveDt::SaveData(){
  SaveFileValuesEnd();
  if(AllDt)SaveFileAllDts();
}


