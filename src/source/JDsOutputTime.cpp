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

/// \file JDsOutputTime.cpp \brief Implements the class \ref JDsOutputTime.

#include "JDsOutputTime.h"
#include "Functions.h"
#include "JXml.h"
#include <cstring>
#include <cfloat>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JDsOutputTime::JDsOutputTime(const JDsOutputTime* tout){
  ClassName="JDsOutputTime";
  Reset();
  if(tout)CopyFrom(tout);
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsOutputTime::~JDsOutputTime(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputTime::Reset(){
  Times.clear();
  TimeBase=0;
  SpecialConfig=false;
  LastTimeInput=LastTimeOutput=-1;
}

//==============================================================================
/// Copy data from other object.
//==============================================================================
void JDsOutputTime::CopyFrom(const JDsOutputTime* tout){
  Reset();
  const unsigned nt=unsigned(tout->Times.size());
  for(unsigned c=0;c<nt;c++)Times.push_back(tout->Times[c]);
  TimeBase=tout->TimeBase;
  SpecialConfig=tout->SpecialConfig;
  LastTimeInput=tout->LastTimeInput;
  LastTimeOutput=tout->LastTimeOutput;
}

//==============================================================================
/// Configures object.
//==============================================================================
void JDsOutputTime::Config(double timeoutdef){
  Reset();
  if(timeoutdef<=0)Run_Exceptioon("Value timeout by default is invalid. It is not greater than zero.");
  AddTimeOut(0,timeoutdef);
  SpecialConfig=false;
}

//==============================================================================
/// Configures object.
//==============================================================================
void JDsOutputTime::Config(std::string filexml,const std::string &place,double timeoutdef){
  Reset();
  JXml xml; xml.LoadFile(filexml);
  LoadXml(&xml,place);
  //-It uses timeoutdef when there is not special configuration in XML the file.
  if(!GetCount())Config(timeoutdef);
  else SpecialConfig=true;
}

//==============================================================================
/// Checks and adds new value of timeout. (returns true if it is wrong).
//==============================================================================
bool JDsOutputTime::AddTimeOut(double t,double tout){
  bool err=(t<0 || tout<=0);
  if(!err && GetCount()){
    StTimeOut last=Times[GetCount()-1];
    err=(last.time>=t);
  }
  if(!err){
    StTimeOut to;
    to.time=t; to.tout=tout;
    Times.push_back(to);
  }
  return(err);
}

//==============================================================================
/// Reads configuration from the XML node.
//==============================================================================
void JDsOutputTime::ReadXml(const JXml *sxml,TiXmlElement* ele){
  TiXmlElement* elet=ele->FirstChildElement("tout"); 
  while(elet){
    double t=sxml->GetAttributeDouble(elet,"time");
    double v=sxml->GetAttributeDouble(elet,"timeout");
    if(AddTimeOut(t,v))sxml->ErrReadElement(elet,"tout",false);
    elet=elet->NextSiblingElement("tout");
  }
}

//==============================================================================
/// Loads configuration from XML object.
//==============================================================================
void JDsOutputTime::LoadXml(const JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNodeSimple(place);
  //if(!node)Run_Exception(string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Shows object configuration.
//==============================================================================
void JDsOutputTime::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())printf("%s\n",txhead.c_str());
  for(unsigned c=0;c<GetCount();c++){
    printf("  Time: %f  (tout:%f)\n",Times[c].time,Times[c].tout);
  }
  if(!txfoot.empty())printf("%s\n",txfoot.c_str());
}

//==============================================================================
/// Returns configuration information.
//==============================================================================
void JDsOutputTime::GetConfig(std::string txhead,std::string txfoot
  ,std::vector<std::string> &lines)const
{
  if(!txhead.empty())lines.push_back(txhead); 
  for(unsigned c=0;c<GetCount();c++){
    lines.push_back(fun::PrintStr("  Time: %f  (tout:%f)",Times[c].time,Times[c].tout)); 
  }
  if(!txfoot.empty())lines.push_back(txfoot); 
}

//==============================================================================
/// Returns next time to save PART file.
//==============================================================================
double JDsOutputTime::GetNextTime(double t){
  if(LastTimeInput>=0 && t==LastTimeInput)return(LastTimeOutput);
  LastTimeInput=t;
  //printf("\n++> GetNextTime(%f)\n",t);
  double nexttime=0;
  double tb=Times[TimeBase].time;
  double tnext=(TimeBase+1<GetCount()? Times[TimeBase+1].time: DBL_MAX);
  //printf("++>INI TimeBase::%d tb:%f \n",TimeBase,tb);
  //-Avanza hasta tb <= t < tnext.
  while(t>=tnext){
    TimeBase++;
    tb=Times[TimeBase].time;
    tnext=(TimeBase+1<GetCount()? Times[TimeBase+1].time: DBL_MAX);
  }
  //printf("++> TimeBase::%d tb:%f \n",TimeBase,tb);
  if(t<tb)nexttime=tb;
  else{
    const double tbout=Times[TimeBase].tout;
    unsigned nt=unsigned((t-tb)/tbout);
    nexttime=tb+tbout*nt;
    for(;nexttime<=t;nt++){
      nexttime=tb+tbout*nt;
      //printf("++> nexttime:%f dif:%f %20.12E\n",nexttime,nexttime-t,nexttime-t);
    }
    if(nexttime>tnext)nexttime=tnext;
  }
  //printf("++> nexttime:%f\n",nexttime);
  LastTimeOutput=nexttime;
  return(nexttime);
}


