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

/// \file JTimeOut.cpp \brief Implements the class \ref JTimeOut.

#include "JTimeOut.h"
#include "Functions.h"
#include "JXml.h"
#include "JLog2.h"
#include <cstring>
#include <cfloat>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JTimeOut::JTimeOut(){
  ClassName="JTimeOut";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JTimeOut::~JTimeOut(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JTimeOut::Reset(){
  Times.clear();
  TimeBase=0;
  SpecialConfig=false;
}

//==============================================================================
/// Configures object.
//==============================================================================
void JTimeOut::Config(double timeoutdef){
  Reset();
  if(timeoutdef<=0)RunException("Config","Value timeout by default is invalid.");
  AddTimeOut(0,timeoutdef);
  SpecialConfig=false;
}

//==============================================================================
/// Configures object.
//==============================================================================
void JTimeOut::Config(std::string filexml,const std::string &place,double timeoutdef){
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
bool JTimeOut::AddTimeOut(double t,double tout){
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
void JTimeOut::ReadXml(JXml *sxml,TiXmlElement* ele){
  const char met[]="ReadXml";
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
void JTimeOut::LoadXml(JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNode(place,false);
  //if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  if(node)ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JTimeOut::VisuConfig(JLog2 *log,std::string txhead,std::string txfoot){
  if(!txhead.empty()){
    if(log)log->Print(txhead); 
    else printf("%s\n",txhead.c_str());
  }
  for(unsigned c=0;c<GetCount();c++){
    if(log)log->Printf("  Time: %f  (tout:%f)",Times[c].time,Times[c].tout);
    else printf("  Time: %f  (tout:%f)\n",Times[c].time,Times[c].tout);
  }
  if(!txfoot.empty()){
    if(log)log->Print(txfoot);
    else printf("%s\n",txfoot.c_str());
  }
}

//==============================================================================
/// Returns next time to save PART file.
//==============================================================================
double JTimeOut::GetNextTime(double t){
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
  return(nexttime);
}


