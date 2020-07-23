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

/// \file JAppInfo.cpp \brief Implements the class \ref JAppInfo.

#include <cstdio>
#include "JAppInfo.h"
#include "Functions.h"

#ifdef JAppInfo_UseLog
  #include "JLog2.h"
#endif

using namespace std;

//##############################################################################
//# JAppInfo
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JAppInfo::JAppInfo(std::string name,std::string ver,std::string date){
  ClassName="JAppInfo";
  #ifdef JAppInfo_UseLog
    Log=NULL;
  #endif
  Reset();
  MainName=name; MainVer=ver; Date=date;
}

//==============================================================================
// Constructor with subname.
//==============================================================================
JAppInfo::JAppInfo(std::string name,std::string ver
  ,std::string subname,std::string subver,std::string date)
{
  ClassName="JAppInfo";
  #ifdef JAppInfo_UseLog
    Log=NULL;
  #endif
  Reset();
  MainName=name; MainVer=ver; Date=date;
  SubName=subname;
  SubVer=subver;
}

//==============================================================================
// Destructor.
//==============================================================================
JAppInfo::~JAppInfo(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JAppInfo::Reset(){
  //-Application information.
  MainName=MainNameExtra=MainVer=Date="";
  SubName=SubVer="";
  //-Execution paths.
  RunCommand=RunPath=ProgramPath="";
  //-Output configuration.
  CreateDirs=true;
  CsvSepComa=false;
  DirOut=DirDataOut="";
  //-Log definition.
  #ifdef JAppInfo_UseLog
    delete Log; Log=NULL;
  #endif
}

//==============================================================================
// Adds extra name.
//==============================================================================
void JAppInfo::AddNameExtra(std::string extra){
  if(!MainNameExtra.empty())MainNameExtra=MainNameExtra+"+";
  MainNameExtra=MainNameExtra+extra;
}

//==============================================================================
// Configures execution paths.
//==============================================================================
void JAppInfo::ConfigRunPaths(std::string runcommand){
  RunCommand=runcommand;
  RunPath=fun::GetCurrentDir();
  ProgramPath=fun::GetDirParent(fun::GetCanonicalPath(RunPath,RunCommand));
}

//==============================================================================
// Configures general output options.
//==============================================================================
void JAppInfo::ConfigOutput(bool createdirs,bool csvsepcoma,std::string dirout,std::string dirdataout){
  CreateDirs=createdirs;
  CsvSepComa=csvsepcoma;
  DirOut=fun::GetDirWithSlash(dirout);
  DirDataOut=(!dirdataout.empty()? fun::GetDirWithSlash(DirOut+dirdataout): DirOut);
}

#ifdef JAppInfo_UseLog
//==============================================================================
// Configures general output options.
// When fname is empty, log file is not created.
//==============================================================================
void JAppInfo::LogInit(std::string fname,bool mpirun,int mpirank,int mpilaunch){
  delete Log; Log=NULL;
  //-Creates directory for log file when it is necessary.
  if(AppInfo.GetCreateDirs() && !fname.empty()){
    fun::MkdirPath(fun::GetDirParent(fname));
  }
  //-Creates object for log.
  Log=new JLog2;
  Log->Init(fname);
}
#endif


//==============================================================================
// Returns short application name.
//==============================================================================
std::string JAppInfo::GetShortName()const{
  return(MainName+(!SubName.empty()? "-": "")+SubName);
}

//==============================================================================
// Returns full application name.
//==============================================================================
std::string JAppInfo::GetFullName()const{
  string extra=(!MainNameExtra.empty()? string(" [")+MainNameExtra+"]": "");
  string subnamever=string(" ")+(!SubName.empty()? SubVer+" ("+MainVer+")": MainVer);
  string date=string(" (")+Date+")"; 
  return(GetShortName()+extra+subnamever+date); 
}

//==============================================================================
// Create directory path when necessary.
//==============================================================================
int JAppInfo::MkdirPath(const std::string &dir)const{
  int ret=0;
  if(AppInfo.GetCreateDirs())ret=fun::MkdirPath(dir);
  return(ret);
}

//==============================================================================
// Create file path when necessary.
//==============================================================================
int JAppInfo::MkdirPathFile(const std::string &file)const{
  int ret=0;
  if(AppInfo.GetCreateDirs())ret=fun::MkdirPath(fun::GetDirParent(file));
  return(ret);
}







