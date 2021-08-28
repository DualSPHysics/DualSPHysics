//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2021 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsTimers.cpp \brief Implements the class \ref JDsTimers.

#include "JDsTimers.h"
#include "Functions.h"
#include "JLog2.h"

using namespace std;

//##############################################################################
//# JDsTimers
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsTimers::JDsTimers(std::string classname){
  ClassName=classname;
  List=new StDsTimer[TIMERSIZE];
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsTimers::~JDsTimers(){
  DestructorActive=true;
  Reset();
  delete[] List; List=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsTimers::Reset(){
  for(int c=0;c<TIMERSIZE;c++)ResetTimer(c);
  CtMax=0;
  SvTimers=false;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsTimers::ResetTimer(unsigned c){
  StDsTimer &t=List[c];
  t.active=false;
  t.timer.Reset();
  t.time=0;
  t.level=0;
  t.name="";
}

//==============================================================================
/// Add timer to list.
//==============================================================================
void JDsTimers::AddTimer(unsigned c,std::string name,unsigned level,bool active){
  if(c>=TIMERSIZE)Run_Exceptioon("Timer Id exceeds the allowed limit.");
  StDsTimer &t=List[c];
  if(!t.name.empty())Run_Exceptioon("Timer Id is already defined.");
  t.name=name;
  t.level=level;
  t.active=active;
  if(c>CtMax)CtMax=c;
}

//==============================================================================
/// Initialises the time accumulated by all timers.
//==============================================================================
void JDsTimers::ResetTimes(){
  for(unsigned c=0;c<=CtMax;c++)List[c].time=0; 
}  

//==============================================================================
/// Returns string with the name of timer and value (empty for inactive timers).
//==============================================================================
std::string JDsTimers::TimerToText(unsigned c,unsigned maxlen)const{
  string ret;
  const StDsTimer &t=List[c];
  if(t.active){
    ret=t.name;
    for(unsigned cv=0;cv<t.level;cv++)ret=string("  ")+ret;
    while(ret.length()<maxlen)ret+=".";
    ret=ret+": ";
    if(t.time)ret=ret+fun::DoubleStr(t.time/1000.,"%f")+" sec.";
    else      ret=ret+"0 sec.";
  }
  return(ret);
}

//==============================================================================
/// Shows active timers.
//==============================================================================
void JDsTimers::ShowTimes(std::string title,JLog2 *log,bool onlyfile)const{
  JLog2::TpMode_Out mode=(onlyfile? JLog2::Out_File: JLog2::Out_ScrFile);
  if(!title.empty())log->Print(title,mode);
  if(!SvTimers)log->Print("none",mode);
  else{
    const unsigned maxlen=33;
    for(unsigned c=0;c<=CtMax;c++){
      const string tx=TimerToText(c,maxlen);
      if(!tx.empty())log->Print(tx,mode);
    }
  }
  log->Print(" ");
}

//==============================================================================
/// Return string with names and values of active timers.
/// Devuelve string con nombres y valores de los timers activos.
//==============================================================================
void JDsTimers::GetTimersInfo(std::string &hinfo,std::string &dinfo)const{
  for(unsigned c=0;c<=CtMax;c++)if(List[c].active){
    hinfo=hinfo+";"+List[c].name;
    dinfo=dinfo+";"+fun::DoubleStr(List[c].time/1000.);
  }
}
