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
//:# Clase para cargar grid de datos de velocidad y SWL.
//:# - Implementacion. (13-08-2020)
//:# - Error corregido en valist. (01-09-2020)
//:# - Soporta definicion de loops desde la carga de datos por fichero. (03-08-2022)
//:#############################################################################

/// \file JMeshTDatasLoad.h \brief Declares the class \ref JMeshTDatasLoad.

#ifndef _JMeshTDatasLoad_
#define _JMeshTDatasLoad_

#include "JObject.h"
#include "TypesDef.h"
#include "JMeshDataDef.h"
#include <string>
#include <vector>
#include <cfloat>

class JBinaryData;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

class JMeshData;

//##############################################################################
//# JMeshTDatasLoad
//##############################################################################
/// \brief Allows loading mesh data using different formats.

class JMeshTDatasLoad : protected JObject
{
 private:
  const std::string AppName;   ///<Application Name.

  static const unsigned FormatVerDef=200817;  ///<Version of file format by default.
  unsigned FormatVer;    ///<Current format version.

  std::string FileIn;    ///<Name of input file.

  const bool FilterTime; ///<Time filter (DBL_MAX=disabled).
  const double TimeMin;  ///<Time filter (DBL_MAX=disabled).
  const double TimeMax;  ///<Time filter (DBL_MAX=disabled).

  const double LoopTmaxRq;   ///<Requested final time of loop (DBL_MAX=disabled).
  const double LoopTbeginRq; ///<Requested begin time of loop (DBL_MAX=disabled).
  double LoopTmax;           ///<Final time of loop (DBL_MAX=disabled).
  double LoopTsub;           ///<Time interval to begin time of loop (0=disabled).
  


  const bool FilterPos;  ///<Position filter (DBL_MAX=disabled).
  const tdouble3 PosMin; ///<Position filter (DBL_MAX=disabled).
  const tdouble3 PosMax; ///<Position filter (DBL_MAX=disabled).

  const double SetTime;  ///<Time offset (applied after filters).
  const tdouble3 SetPos; ///<Position offset (applied after filters).

 private:
  void LoadDataDefs(const JBinaryData* bdat,std::string varlist
    ,std::vector<StMeshData>& datadef);
  void LoadFileBin(std::vector<JMeshData*>& vdata,std::string varlist="");

  void LoadFileCsv(std::vector<JMeshData*>& vdata);

  StMeshPts ComputeFilterPosRanges(StMeshPts mp,unsigned& c1ini,unsigned& c1fin
   ,unsigned& c2ini,unsigned& c2fin,unsigned& c3ini,unsigned& c3fin)const;
  void RunFilterPos(std::vector<JMeshData*>& vdata)const;

  void SetTimePos(std::vector<JMeshData*>& vdata)const;

 public:
  JMeshTDatasLoad(double tmin=DBL_MAX,double tmax=DBL_MAX
    ,double looptmax=DBL_MAX,double looptbegin=DBL_MAX
    ,tdouble3 pmin=TDouble3(DBL_MAX),tdouble3 pmax=TDouble3(DBL_MAX)
    ,double settime=0,tdouble3 setpos=TDouble3(0));
  ~JMeshTDatasLoad();
  void Reset();

  void LoadFile(const std::string file,std::vector<JMeshData*>& vdata,std::string varlist="");

  double GetLoopTmax()const{ return(LoopTmax); }
  double GetLoopTsub()const{ return(LoopTsub); }
};

}

#endif


