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

//:#############################################################################
//:# Cambios:
//:# =========
//:# Clase para gestionar informacion general de la aplicacion y proporcionar 
//:# funcionalidades transverales. (21-05-2018)
//:# - Se incluye extraname en constructor y se elimina ConfigNameExtra(). (21-05-2018)
//:# - Nuevos metodos MkdirPath() y MkdirPathFile(). (21-05-2018)
//:# - Nuevos metodos ClearNameExtra() y AddNameExtra(). (23-05-2018)
//:# - Nuevos metodos GetMainName(), GetMainVer() y GetDate(). (07-03-2019)
//:# - El uso de JLog2 o no se define en JAppInfoDef.h. (17-06-2020)
//:#############################################################################

/// \file JAppInfo.h \brief Declares the class \ref JAppInfo and the global object AppInfo.

#ifndef _JAppInfo_
#define _JAppInfo_

#include "JAppInfoDef.h"
#include "JObject.h"

#include <string>

class JLog2;

//##############################################################################
//# JAppInfo
//##############################################################################
/// \brief Defines a global object with general information and resources for the whole application.

class JAppInfo : public JObject
{
private:
  //-Application information.
  std::string MainName;
  std::string MainNameExtra;
  std::string MainVer;
  std::string SubName;
  std::string SubVer;
  std::string Date;

  //-Execution paths.
  std::string RunCommand;
  std::string RunPath;
  std::string ProgramPath;

  //-Output configuration.
  bool CreateDirs;   ///<Creates full path for output files (true by default).
  bool CsvSepComa;   ///<Separator character in CSV files (0=semicolon, 1=coma).
  std::string DirOut;
  std::string DirDataOut;

  //-Log definition.
#ifdef JAppInfo_UseLog
  JLog2* Log;
#endif

public:
  JAppInfo(std::string name,std::string ver,std::string date);
  JAppInfo(std::string name,std::string ver,std::string subname,std::string subver,std::string date);
  ~JAppInfo();
  void Reset();
  void ClearNameExtra(){ MainNameExtra=""; };
  void AddNameExtra(std::string extra);

  void ConfigRunPaths(std::string runcommand);
  void ConfigOutput(bool createdirs,bool csvsepcoma,std::string dirout,std::string dirdataout="");

  void SetMainName(const std::string &mname){ MainName=mname; }

#ifdef JAppInfo_UseLog
  void LogInit(std::string fname,bool mpirun=false,int mpirank=0,int mpilaunch=0);
  JLog2* LogPtr(){ return(Log); }
  bool LogDefined()const{ return(Log!=NULL); }
#else
  bool LogDefined()const{ return(false); }
#endif

  std::string GetShortName()const;
  std::string GetFullName()const;

  std::string GetMainName()const{ return(MainName); }
  std::string GetMainVer()const{ return(MainVer); }
  std::string GetDate()const{ return(Date); }

  std::string GetRunCommand()const{ return(RunCommand); };
  std::string GetRunPath()const{ return(RunPath); };
  std::string GetProgramPath()const{ return(ProgramPath); };

  bool GetCreateDirs()const{ return(CreateDirs); };
  bool GetCsvSepComa()const{ return(CsvSepComa); };
  std::string GetDirOut()const{ return(DirOut); };
  std::string GetDirDataOut()const{ return(DirDataOut); };

  int MkdirPath(const std::string &dir)const;
  int MkdirPathFile(const std::string &file)const;

};

#endif

/// \brief Declares the global object AppInfo with general information and resources for the whole application.
extern JAppInfo AppInfo;

