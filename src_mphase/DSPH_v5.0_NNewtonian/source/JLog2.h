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
//:# - Soporte para MPI. (28-10-2011)
//:# - Nuevo metodo PrintDbg() que por defecto hace un fflush(stdout). (11-01-2012)
//:# - Traduccion de comentarios al ingles. (10-02-2012)
//:# - Nuevo metodo Printf() y PrintfDbg() usando argumentos variables. (10-03-2014)
//:# - Nuevo metodo GetDirOut() para obtener el directiorio de salida. (07-05-2014)
//:# - Nuevo metodo Print() para std::vector<std::string>. (01-02-2017)
//:# - Error coregido en Printf() y PrintfDbg() cuando el string era demasiado 
//:#   grande. (20-02-2017)
//:# - Nuevas funciones Printp() y Printfp() a las que se le puede anhadir un 
//:#   prefijo. (21-02-2017)
//:# - New attribute CsvSepComa to configure separator in CSV files. (24-10-2017)
//:# - Se incluye DirDataOut para facilitar su uso en distintos ambitos. (19-02-2017)
//:# - Funciones para gestion especial de warnings. (10-03-2018)
//:# - Permite usar un Mutex para sicronizar el acceso en multithreading ejecuciones. (22-08-2019)
//:# - Permite crear un log como referencia a otro para incluir un prefijo de forma automatica. (06-09-2019)
//:# - El uso de mutex o no se define en JLog2Def.h. (17-06-2020)
//:#############################################################################

/// \file JLog2.h \brief Declares the class \ref JLog2.

#ifndef _JLog2_
#define _JLog2_

#include "JObject.h"

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>

//##############################################################################
//# JLog2
//##############################################################################
/// \brief Manages the output of information in the file Run.out and on screen

class JLog2 : protected JObject
{
public:
  typedef enum{ Out_Default=4,Out_ScrFile=3,Out_File=2,Out_Screen=1,Out_None=0 }TpMode_Out;

  ///Structure to save file descriptions.
  typedef struct StrFileInfo{
    std::string file;
    std::string info;
    StrFileInfo(const std::string &xfile,const std::string &xinfo){ file=xfile; info=xinfo; }
  }StFileInfo;

protected:
  JLog2 *Parent;
  std::string ParentPrefix;

  std::string FileName;
  std::ofstream *Pf;
  bool Ok;
  bool MpiRun;
  int MpiRank,MpiLaunch;
  TpMode_Out ModeOutDef;

  std::vector<std::string> Warnings; ///<List of warnings.

  std::vector<StFileInfo> FileInfo; ///<List of file descriptions.

  //-General output configuration.
  std::string DirOut;      ///<Specifies the general output directory.

public:
  JLog2(TpMode_Out modeoutdef=Out_ScrFile);
  JLog2(JLog2 *parent,std::string prefix);
  ~JLog2();
  void Reset();
  void Init(std::string fname,bool mpirun=false,int mpirank=0,int mpilaunch=0);
  void SetModeOutDef(TpMode_Out modeoutdef){ ModeOutDef=modeoutdef; }
  void Print(const std::string &tx,TpMode_Out mode=Out_Default,bool flush=false);
  void Print(const std::vector<std::string> &lines,TpMode_Out mode=Out_Default,bool flush=false);
  void PrintDbg(const std::string &tx,TpMode_Out mode=Out_Default){ Print(tx,mode,true); }
  void PrintFile(const std::string &tx,bool flush=false){ Print(tx,Out_File,flush); }
  bool IsOk()const{ return(Ok); }
  int GetMpiRank()const{ return(MpiRun? MpiRank: -1); }
  std::string GetParentPrefix()const{ return(ParentPrefix); }

  //std::string GetDirOut()const{ return(DirOut); }

  void Printf(const char *format,...);
  void PrintfDbg(const char *format,...);
  //-Adding a prefix.
  void Printp(const std::string &prefix,const std::string &tx,TpMode_Out mode=Out_Default,bool flush=false){ Print(prefix+tx,mode,flush); }
  void Printp(const std::string &prefix,const std::vector<std::string> &lines,TpMode_Out mode=Out_Default,bool flush=false);
  void PrintpDbg(const std::string &prefix,const std::string &tx,TpMode_Out mode=Out_Default){ Printp(prefix,tx,mode,true); }
  void Printfp(const std::string &prefix,const char *format,...);
  void PrintfpDbg(const std::string &prefix,const char *format,...);

  //-Warning system.
  void AddWarning(const std::string &tx);
  void PrintWarning(const std::string &tx,TpMode_Out mode=Out_Default,bool flush=false);
  void PrintfWarning(const char *format,...);
  unsigned WarningCount()const{ return(unsigned(Warnings.size())); }
  void PrintWarningList(const std::string &txhead,const std::string &txfoot,TpMode_Out mode=Out_Default,bool flush=false);
  void PrintWarningList(TpMode_Out mode=Out_Default,bool flush=false);

  //-File description.
  void AddFileInfo(std::string fname,const std::string &finfo);
  unsigned FilesCount()const{ return(unsigned(FileInfo.size())); }
  void PrintFilesList(const std::string &txhead,const std::string &txfoot,TpMode_Out mode=Out_Default,bool flush=false);
  void PrintFilesList(TpMode_Out mode=Out_Default,bool flush=false);

};

#endif


