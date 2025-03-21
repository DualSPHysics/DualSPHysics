//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# Descripcion:
//:# =============
//:# Clase para grabar malla de datos.
//:# - Implementacion. (09-08-2020)
//:# - Permite grabar ficheros mbi4 en formato multidata. (24-10-2021)
//:#############################################################################

/// \file JMeshTDatasSave.h \brief Declares the class \ref JMeshTDatasSave.

#ifndef _JMeshTDatasSave_
#define _JMeshTDatasSave_

#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"
#include "JMeshDataDef.h"
#include <string>
#include <vector>
#include <fstream>

namespace jcsv{
  class JSaveCsv2;
}

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

class JMeshData;

//##############################################################################
//# JMeshTDatasSave
//##############################################################################
/// \brief Allows writing mesh data using different formats.

class JMeshTDatasSave : protected JObject
{
 private:
  JBinaryData* Data;      ///<Stores general information of the grid data (constant for each record).
  JBinaryData* DataTime;  ///<Belongs to Data and stores information of an instant.
  unsigned Ctime;         ///<Number of DataTimes.

  //-Management of variables.
  static const unsigned FormatVerDef=200817;  ///<Version of file format by default.
  unsigned FormatVer;    ///<Current format version.

  std::string AppName;   ///<Application Name.
  std::string SaveFile;  ///<File name to save data.

  StMeshPts Mesh;      ///<Mesh definition.
  unsigned Npt;        ///<Total number of points (Npt=Npt1*Npt2*Npt3).
  unsigned Npt12;      ///<Total number of points (Npt12=Npt1*Npt2).
  std::vector<StMeshData> Datdefs; ///<Data definitions.

  bool InitialSaved;   ///<Indicates if header information is saved.

 private:
  static void RunExceptioonStatic(const std::string& srcfile,int srcline
    ,const std::string& method
    ,const std::string& msg,const std::string& file="");

  static std::string GetNameDataTime(unsigned ctime);
  void CheckDataDefinitions(const JMeshData* mdat);

 public:
  JMeshTDatasSave();
  ~JMeshTDatasSave();
  void Reset();
  void ResetData();
  void ResetDataTime();

  llong GetAllocMemory()const;

  void Config(const std::string file,std::string appname,const JMeshData* mdat);
  void SaveInitial();

  void SaveDataTime(const JMeshData* mdat);
  void SaveDataTimes(unsigned n,JMeshData** vmdat);

  static void SaveCsvHead(const StMeshPts& m,std::string name,std::string tunits
    ,TpTypeData type,bool data12,jcsv::JSaveCsv2& scsv);
  static void SaveCsv(const JMeshData* mdat,unsigned cdata,jcsv::JSaveCsv2& scsv,bool svhead);
  
  static void SaveVtk(std::string file,int fnum,const JMeshData* mdat,bool svdata12);

  static void SaveVtkScheme(std::string file,StMeshPts m);

};



}

#endif


