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

//#############################################################################
//# Cambios:
//# =========
//# - Clase para grabar informacion de los floatings con mayor frecuencia.
//#
//# Cambios:
//# =========
//# - Implementacion inicial. (19-10-2020)
//# - Nueva implementacion simplificada basada en JFtMotionSave. (23-07-2023)
//#############################################################################

/// \file JDsPartFloatSave.h \brief Declares the class \ref JDsPartFloatSave.

#ifndef _JDsPartFloatSave_
#define _JDsPartFloatSave_

#include "JObject.h"
#include "TypesDef.h"
#include "DualSphDef.h"
#include <string>

#ifdef _WITHGPU
#include <cuda_runtime_api.h>
#endif

class JLog2;
class JPartFloatInfoBi4Save;
class JPartFloatInfoBi4Data;
class JDsFtForcePoints;

//##############################################################################
//# JDsPartFloatSave
//##############################################################################
/// \brief Saves floating data with several frequencies.

class JDsPartFloatSave : protected JObject
{
private:
  JLog2* Log;

  const bool Cpu;
  const std::string AppName; ///<Application Name.
  const std::string Dir;     ///<Data Directory.
  const double TimeOut;      ///<Defined output period (for information only).
  const double TimeOut2;     ///<Defined second output period when TimeOut2>=0 (for information only).

  const word MkBoundFirst;   ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  const unsigned FtCount;    ///<Number of floating bodies.
  unsigned FptCount;         ///<Number of force points.

  JPartFloatInfoBi4Data* FtData; ///Store last floating data.
  
  JPartFloatInfoBi4Save* FtFile1; ///Main object to save floating data according to TimeOut as file ibi4.
  JPartFloatInfoBi4Save* FtFile2; ///Extra object to save floating data according to TimeOut2 as file ibi4.

public:
  JDsPartFloatSave(bool cpu,std::string appname,std::string dir,word mkboundfirst
    ,unsigned ftcount,double timeout,double timeout2);
  ~JDsPartFloatSave();
  void Reset();

  void ConfigFtData(unsigned ftcount,const StFloatingData* ftobjs);
  void ReconfigForcePoints(const JDsFtForcePoints* forcepoints);

  void SaveInitial();

  bool ExtraIsActive()const{ return(FtFile2!=NULL); };
  std::string GetFileName(bool mainfile)const;

  void SetFtData(int cpart,double timestep,unsigned nstep
    ,const StFloatingData* ftobjs,const JDsFtForcePoints* forcepoints);

  void SaveDataMain(int cpart,double timestep,unsigned nstep);
  
  void AddDataExtra(int cpart,double timestep,unsigned nstep);
  void SaveDataExtra();
};

#endif

