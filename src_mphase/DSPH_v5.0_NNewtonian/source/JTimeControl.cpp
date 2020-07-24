//HEAD_DSCODES
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

/// \file JTimeControl.cpp \brief Implements the class \ref JTimeControl.

#include "JTimeControl.h"
#include "Functions.h"
#include <cstring>

using namespace std;

//##############################################################################
//# JTimeControl
//##############################################################################
//==============================================================================
/// Constructor for periodic evaluation.
/// Constructor para evaluacion periodica.
//==============================================================================
JTimeControl::JTimeControl(double tout){
  ClassName="JTimeControl";
  ConfigPeriodic(0,tout);
}

//==============================================================================
/// Constructor for periodic evaluation using an initial time.
/// Constructor para evaluacion periodica con un instane inicial.
//==============================================================================
JTimeControl::JTimeControl(double tfirst,double tout){
  ClassName="JTimeControl";
  ConfigPeriodic(tfirst,tout);
}

//==============================================================================
/// Constructor.
//==============================================================================
JTimeControl::JTimeControl(unsigned ntimes,const double *vtimes){
  ClassName="JTimeControl";
  ConfigTimes(ntimes,vtimes);
}

//==============================================================================
/// Constructor.
//==============================================================================
JTimeControl::JTimeControl(const std::string &times){
  ClassName="JTimeControl";
  ConfigTimes(times);
}

////==============================================================================
///// Destructor.
////==============================================================================
//JTimeControl::~JTimeControl(){
//}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JTimeControl::Reset(){
  Active=true;
  Timer.Reset();
  NextTime=0;
  Periodic=false;
  FirstTime=TimeOut=0;
  TimeOutNum=0;
  Times.clear();
  IteStart=0;
  IteNum=NextIte=0;
  LastTime=0;
}

//==============================================================================
/// Configures for periodic evaluation.
/// Configura para evaluacion periodica.
//==============================================================================
void JTimeControl::ConfigPeriodic(double tfirst,double tout){
  Reset();
  Periodic=true;
  FirstTime=(tfirst<tout? tfirst: 0);
  TimeOut=tout;
  NextTime=CalcNextTime();
  Timer.Start();
}

//==============================================================================
/// Configures list of times.
/// Configura lista de tiempos.
//==============================================================================
void JTimeControl::ConfigTimes(unsigned ntimes,const double *vtimes){
  Reset();
  for(unsigned c=0;c<ntimes;c++)Times.push_back(vtimes[c]);
  PrepareTimes();
}

//==============================================================================
/// Loads list of times from a string using the separator ','.
/// Configura lista de tiempos a partir de string con valoes separados por ','.
//==============================================================================
void JTimeControl::ConfigTimes(std::string times){
  Reset();
  while(!times.empty()){
    std::string value=fun::StrSplit(",",times);
    if(!value.empty())Times.push_back(atof(value.c_str()));
  }
  PrepareTimes();
}

//==============================================================================
/// Prepares and sort Times[] to use it.
/// Prepara y ordena Times[] para empezar a usarlo.
//==============================================================================
void JTimeControl::PrepareTimes(){
  TimesSize=unsigned(Times.size());
  //-Sorts values from lowest to highest. | Ordena valores de menor a mayor.
  for(unsigned c=0;c<TimesSize;c++)for(unsigned c2=c+1;c2<TimesSize;c2++)if(Times[c]>Times[c2])swap(Times[c],Times[c2]);
  NextTime=CalcNextTime();
  Timer.Start();
}

//==============================================================================
/// Returns next NextTime.
/// Devuelve siguiente NextTime.
//==============================================================================
double JTimeControl::CalcNextTime(){
  double ret=0;
  if(Periodic){
    TimeOutNum++;
    ret=double(TimeOutNum)*TimeOut;
    if(FirstTime){
      TimeOutNum--;
      ret=FirstTime; 
      FirstTime=0;
    }
  }
  else{
    if(TimeOutNum<TimesSize){ 
      ret=Times[TimeOutNum];
      TimeOutNum++;
    }
    else Active=false;
  }
  //:printf("---------->CalcNextTime:%f \n",ret);
  return(ret);
}

//==============================================================================
/// Returns true when it changes to next NextTime.
/// Devuelve true si paso el siguiente NextTime.
//==============================================================================
bool JTimeControl::CheckRealTime(){
  bool ret=false;
  Timer.Stop();
  const double t=Timer.GetElapsedTimeD()/1000.0;//-En segundos.
  //:printf("---------->nite:%u t:%f \n",IteNum,t);
  if(t>NextTime){
    ret=true;
    LastTime=t;
    NextTime=CalcNextTime();
    //-Initializes IteNum.
    IteStart=t;
    IteNum=NextIte=0;
  }
  if(IteNum>10){
    const unsigned nite=unsigned(double(((NextTime-t)*IteNum)/(t-IteStart))/10.);
    //:printf("--> nite:%u  t:%f\n",nite,(NextTime-t));
    NextIte=IteNum+(nite<1000000000? nite: 1000000000); //-The maximum is 1000M.
    if(NextIte<IteNum){//-It has gone over and initializes IteSum. | Se ha pasado de vueltas e inicializa IteNum.
      IteStart=t;
      IteNum=0; NextIte=nite;
    }
  }
  return(ret);
}

//==============================================================================
/// Returns text with inoformation.
/// Devuelve texto "Progress 12.5% in 15.1 seconds. Estimated finish time 26-04-2016 12:21:23"
/// done es progreso actual en tanto por 1.
//==============================================================================
std::string JTimeControl::GetInfoFinish(double done){
  return(fun::PrintStr("Progress %.2f%% in %.1fs. Estimated finish time: ",done*100,LastTime)+fun::GetDateTimeAfter(int((LastTime/done)-LastTime)));
}


