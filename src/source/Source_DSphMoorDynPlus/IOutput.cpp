/*===================================================================================
<MOORDYNPLUS> Copyright (c) 2025
Ivan Martinez Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
Matt Hall (github.com/mattEhall)

This file is part of MoorDynPlus. MoorDynPlus is free software: you can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

Linking the MoorDynPlus library statically or dynamically with other modules is
making a combined work based on this library. Thus, the terms and conditions
of the GNU General Public License cover the whole combination. As a special
exception, the copyright holders of MoorDynPlus give you permission to dynamically
link this library with the program DualSPHysics to produce a combined model
featuring the capabilities of both DualSPHysics and MoorDynPlus. This exception
is strictly limited to linking between the compiled MoorDynPlus library and
DualSPHysics. It does not extend to other programs or the use of the MoorDynPlus
source code beyond the stipulations of the GPL. When the exception is used,
this paragraph must be included in the copyright notice.

MoorDynPlus is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.

You should have received a copy of the GNU General Public License along with
MoorDynPlus. If not, see <http://www.gnu.org/licenses/>.
===================================================================================*/

/// \file IOutput.cpp \brief Implements the class \ref IOutput.

#include "IOutput.h"
#include "FunMoorDynPlus.h"
#include "Functions.h"
#include "JSaveCsv2.h"

//==============================================================================
/// Constructor.
//==============================================================================
IOutProperties::IOutProperties() {
  Reset();
  ClassName="IOutProperties";
}

//==============================================================================
/// Destructor.
//==============================================================================
IOutProperties::~IOutProperties() {
  Reset();
  Lines.clear();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void IOutProperties::Reset() {
  LineCount=0;
  Units="";
}

//==============================================================================
/// Constructor.
//==============================================================================
IOutput::IOutput(double timemax,double dtout) {
  Reset();
  ClassName="IOutput";
  TimeMax=timemax;
  DtOut=dtout;
  IsFirst=true;
}

//==============================================================================
/// Destructor.
//==============================================================================
IOutput::~IOutput() {
  fmdp::FreeVector(OutProps);
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void IOutput::Reset() {
  LineCount=0;
  TypeCount=0;
  TimeStart=0;
  TimeMax=0;
  NextTime=0;
  FileName="";
  IsFirst=false;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void IOutput::LoadXml(JXml* sxml,const std::string& place,std::vector<ILine*> lines) {
  std::string function="LoadXml";

  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
  ReadXml(sxml,node->ToElement(),lines);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void IOutput::ReadXml(JXml* sxml,TiXmlElement* ele,std::vector<ILine*> lines) {
  LineCount=(unsigned)lines.size();

  TiXmlElement* eleOut;
  if     (sxml->ExistsElement(ele,"savedata"))eleOut=ele->FirstChildElement("savedata");
  else if(sxml->ExistsElement(ele,"output")  )eleOut=ele->FirstChildElement("output");//-For compatibility with old versions

  if(eleOut) {
    TiXmlElement* eleTime=eleOut->FirstChildElement("time");
    if(eleTime) {
      TimeStart=sxml->GetAttributeDouble(eleTime,"starTime",true,0);
      TimeMax=sxml->GetAttributeDouble(eleTime,"endTime",true,TimeMax);
      DtOut=sxml->GetAttributeDouble(eleTime,"dtOut",true,DtOut);
      if(!TimeMax)TimeMax=DBL_MAX;
      NextTime=TimeStart;
    }
    TypeCount=0;
    //-Tensions
    if(sxml->ExistsElement(eleOut,"tension")) {
      TiXmlElement* elet=eleOut->FirstChildElement("tension");
      if(sxml->GetAttributeBool(elet,"value",true,true)){
        OutProps.push_back(new IOutProperties());
        OutProps[TypeCount]->SetMagnitude(TpMagnitude::tension);
        SetLinesSelected(OutProps[TypeCount],lines);
        OutProps[TypeCount]->SetUnits("[N]");
        TypeCount++;
      }
    }
    //-Forces
    if(sxml->ExistsElement(eleOut,"force")) {
      TiXmlElement* elef=eleOut->FirstChildElement("force");
      if(sxml->GetAttributeBool(elef,"value",true,true)){
        OutProps.push_back(new IOutProperties());
        OutProps[TypeCount]->SetMagnitude(TpMagnitude::force);
        SetLinesSelected(OutProps[TypeCount],lines);
        OutProps[TypeCount]->SetUnits("[N]");
        TypeCount++;
      }
    }
    //-Velocities
    if(sxml->ExistsElement(eleOut,"velocity")) {
      TiXmlElement* elev=eleOut->FirstChildElement("velocity");
      if(sxml->GetAttributeBool(elev,"value",true,true)){
        OutProps.push_back(new IOutProperties());
        OutProps[TypeCount]->SetMagnitude(TpMagnitude::velocity);
        SetLinesSelected(OutProps[TypeCount],lines);
        OutProps[TypeCount]->SetUnits("[m/s]");
        TypeCount++;
      }
    }
    //-Positions
    if(sxml->ExistsElement(eleOut,"position")) {
      TiXmlElement* elep=eleOut->FirstChildElement("position");
      if(sxml->GetAttributeBool(elep,"value",true,true)){
        OutProps.push_back(new IOutProperties());
        OutProps[TypeCount]->SetMagnitude(TpMagnitude::position);
        SetLinesSelected(OutProps[TypeCount],lines);
        OutProps[TypeCount]->SetUnits("[m]");
        TypeCount++;
      }
    }
    TypeCount=(unsigned)OutProps.size();
    eleOut=eleOut->NextSiblingElement();
  }
}

//==============================================================================
/// Initializes output file
//==============================================================================
void IOutput::Setup(std::string dir){
  std::string function="Setup";
  DataDir=dir;
}

//==============================================================================
/// Checks the string passed and store in IOutProperties the lines selected
//==============================================================================
void IOutput::SetLinesSelected(IOutProperties* oProps,std::vector<ILine*> lines) {
  for(unsigned nl=0;nl<LineCount;nl++) lines[nl]->SetDtOut(DtOut); 
  oProps->SetNLines(LineCount);
  oProps->SetLines(lines);
}

//==============================================================================
/// Write the outputs for each line
//==============================================================================
void IOutput::SaveCsv(double timestep,double dt) {
  const double start=TimeStart;  //-Start time. | Tiempo de inicio.
  const double finish=TimeMax;   //-End time. | Tiempo de finalizacion.
  const double dtout=DtOut;      //-Save frequency. | Frecuencia de guardado.
  bool savedata=false;
  
  //if(start<=timestep && timestep<=finish) savedata=!(timestep<(floor((timestep-dt)/dtout)+1.0)*dtout);
  if(start<=timestep && timestep<=finish) savedata=NextTime<=timestep;

  if(savedata){
    const std::string path=DataDir+"/";
    if(!fun::DirExists(path))fun::MkdirPath(path);
    for(unsigned lf=0; lf<TypeCount; lf++){
      const IOutProperties* out_p=OutProps[lf];
      const TpMagnitude mag=out_p->GetMagnitude();
      if(mag==TpMagnitude::tension) SaveCsvTen(path,out_p,timestep,dt); //-Tensions
      if(mag==TpMagnitude::force)   SaveCsvFor(path,out_p,timestep,dt); //-Forces
      if(mag==TpMagnitude::velocity)SaveCsvVel(path,out_p,timestep,dt); //-Velocities
      if(mag==TpMagnitude::position)SaveCsvPos(path,out_p,timestep,dt); //-Position
    }
    //-Disable the first save. 
    if(IsFirst)IsFirst=false;
    //-Calculates the next time to save data
    NextTime=NextTime+dtout;
  }
}

//==============================================================================
/// Write the tensions for each line
//==============================================================================
void IOutput::SaveCsvTen(const std::string path,const IOutProperties* out_p,double timestep,double dt){
  std::string units="";
  const unsigned nlines=out_p->GetNLines();
  std::vector<ILine*> lines=out_p->GetLines();
  const std::string file=path+fun::PrintStr("MoorDynPlus_tension.csv");
  jcsv::JSaveCsv2 scsv(file,!IsFirst,0);
  if(!scsv.GetAppendMode()){
    //-Saves head.
    //-Guarda cabecera.
    scsv.SetHead();
    scsv << "time [s];dt [s]";
    for(unsigned l=0; l<nlines; l++) {
      units=out_p->GetUnits();
      const ILine* line=lines[l];
      scsv << fun::PrintStr("Anchor_L%u %s",line->GetId(),units.c_str());
      scsv << fun::PrintStr("Fairlead_L%u %s",line->GetId(),units.c_str());
    }
    scsv << jcsv::Endl();
  }
  //-Saves data.
  //-Guarda datos.
  scsv.SetData();
  scsv << timestep << dt;
  for(unsigned l=0; l<nlines; l++) {
    const ILine* line=lines[l];
    scsv << line->GetTensionOutput(0);
    scsv << line->GetTensionOutput(lines[l]->GetN());
  }
  scsv << jcsv::Endl();
  scsv.SaveData();
}

//==============================================================================
/// Write the forces for each line
//==============================================================================
void IOutput::SaveCsvFor(const std::string path,const IOutProperties* out_p,double timestep,double dt){
  std::string units="";
  const unsigned nlines=out_p->GetNLines();
  std::vector<ILine*> lines=out_p->GetLines();
  const std::string file=path+fun::PrintStr("MoorDynPlus_force.csv");
  jcsv::JSaveCsv2 scsv(file,!IsFirst,0);
  if(!scsv.GetAppendMode()){
    //-Saves head.
    //-Guarda cabecera.
    scsv.SetHead();
    scsv << "time [s];dt [s]";
    for(unsigned l=0; l<nlines; l++) {
      units=out_p->GetUnits();
      const ILine* line=lines[l];
      scsv << fun::PrintStr("Anchor.x_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Anchor.y_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Anchor.z_L%u %s",line->GetId(),units.c_str());
      scsv << fun::PrintStr("Fairlead.x_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Fairlead.y_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Fairlead.z_L%u %s",line->GetId(),units.c_str());
    }
    scsv << jcsv::Endl();
  }
  //-Saves data.
  //-Guarda datos.
  scsv.SetData();
  scsv << timestep << dt;
  for(unsigned l=0; l<nlines; l++) {
    const ILine* line=lines[l];
    tdouble3 out_val=line->GetForceOutput(0);
    scsv << out_val.x << out_val.y << out_val.z;
    out_val=line->GetForceOutput(line->GetN());
    scsv << out_val.x << out_val.y << out_val.z;
  }
  scsv << jcsv::Endl();
  scsv.SaveData();
}

//==============================================================================
/// Write the velocities for each line
//==============================================================================
void IOutput::SaveCsvVel(const std::string path,const IOutProperties* out_p,double timestep,double dt){
  std::string units="";
  const unsigned nlines=out_p->GetNLines();
  std::vector<ILine*> lines=out_p->GetLines();
  const std::string file=path+fun::PrintStr("MoorDynPlus_velocity.csv");
  jcsv::JSaveCsv2 scsv(file,!IsFirst,0);
  if(!scsv.GetAppendMode()){
    //-Saves head.
    //-Guarda cabecera.
    scsv.SetHead();
    scsv << "time [s];dt [s]";
    for(unsigned l=0; l<nlines; l++) {
      units=out_p->GetUnits();
      const ILine* line=lines[l];
      scsv << fun::PrintStr("Anchor.x_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Anchor.y_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Anchor.z_L%u %s",line->GetId(),units.c_str());
      scsv << fun::PrintStr("Fairlead.x_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Fairlead.y_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Fairlead.z_L%u %s",line->GetId(),units.c_str());
    }
    scsv << jcsv::Endl();
  }
  //-Saves data.
  //-Guarda datos.
  scsv.SetData();
  scsv << timestep << dt;
  for(unsigned l=0; l<nlines; l++) {
    const ILine* line=lines[l];
    tdouble3 out_val=line->GetVelocityOutput(0);
    scsv << out_val.x  << out_val.y << out_val.z;
    out_val=line->GetVelocityOutput(line->GetN());
    scsv << out_val.x << out_val.y << out_val.z;
  }
  scsv << jcsv::Endl();
  scsv.SaveData();
}

//==============================================================================
/// Write the positions for each line
//==============================================================================
void IOutput::SaveCsvPos(const std::string path,const IOutProperties* out_p,double timestep,double dt){
  std::string units="";
  const unsigned nlines=out_p->GetNLines();
  std::vector<ILine*> lines=out_p->GetLines();
  const std::string file=path+fun::PrintStr("MoorDynPlus_position.csv");
  jcsv::JSaveCsv2 scsv(file,!IsFirst,0);
  if(!scsv.GetAppendMode()){
    //-Saves head.
    //-Guarda cabecera.
    scsv.SetHead();
    scsv << "time [s];dt [s]";
    for(unsigned l=0; l<nlines; l++) {
      units=out_p->GetUnits();
      const ILine* line=lines[l];
      scsv << fun::PrintStr("Anchor.x_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Anchor.y_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Anchor.z_L%u %s",line->GetId(),units.c_str());
      scsv << fun::PrintStr("Fairlead.x_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Fairlead.y_L%u %s",line->GetId(),units.c_str())
        << fun::PrintStr("Fairlead.z_L%u %s",line->GetId(),units.c_str());
    }
    scsv << jcsv::Endl();
  }
  //-Saves data.
  //-Guarda datos.
  scsv.SetData();
  scsv << timestep << dt;
  for(unsigned l=0; l<nlines; l++) {
    const ILine* line=lines[l];
    tdouble3 out_val=line->GetPositionOutput(0);
    scsv << out_val.x << out_val.y  << out_val.z;
    out_val=line->GetPositionOutput(line->GetN());
    scsv << out_val.x  << out_val.y  << out_val.z;
  }
  scsv << jcsv::Endl();
  scsv.SaveData();
}