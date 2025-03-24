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
//:# Clase grabar posiciones de referencia de moving y floating bodies con 
//:# distintas frecuencias a partir de datos de particulas en CPU y GPU.
//:#
//:# Cambios:
//:# =========
//:# - Implementacion inicial. (06-07-2023)
//:#############################################################################

/// \file JDsPartMotionSave.h \brief Declares the class \ref JDsPartMotionSave.

#ifndef _JDsPartMotionSave_
#define _JDsPartMotionSave_

#include "JObject.h"
#include "TypesDef.h"
#include "DualSphDef.h"
#include "JPartMotionDef.h"
#include <string>
#include <vector>

#ifdef _WITHGPU
#include <cuda_runtime_api.h>
#endif

class JLog2;
class JPartMotRefBi4Save;

//##############################################################################
//# JDsPartMotionSave
//##############################################################################
/// \brief Allows writing motion reference data of moving and floating objects during simulation.

class JDsPartMotionSave : protected JObject
{
 private:
  JLog2* Log;

  const bool Cpu;
  const std::string AppName;  ///<Application Name.
  const std::string Dir;      ///<Data Directory.
  const double TimeOut;       ///<Defined output period (for information only).
  const double TimeOut2;      ///<Defined second output period when TimeOut2>=0 (for information only).

  const unsigned CaseNfixed;  ///<Number of fixed boundary particles. 
  const unsigned CaseNmoving; ///<Number of moving boundary particles. 
  const unsigned CaseNfloat;  ///<Number of floating boundary particles. 
  const word MkBoundFirst;    ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).

  const unsigned MkCount;     ///<Number of moving and floating bodies.
  const unsigned PsCount;     ///<Number of reference positions (PsCount=MkCount*3).

  unsigned MkMovingCount; ///<Number of moving bodies.
  unsigned MkFloatCount;  ///<Number of floating bodies.
  StMkMotionData* MkMotionData; ///<Information of moving/floating boundary Mk blocks [MkCount].

  unsigned IdpBegin;  ///<First idp of selection (first moving idp).
  unsigned IdpCount;  ///<Size of selection.

  unsigned LastStep;
  double   LastTimestep;

  //-Auxiliary memory.
  unsigned* IdpRef;  ///<Idp of selected particles [PsCount].
  tdouble3* PosRef;  ///<Position of selected particles [PsCount].

  //-Auxiliary memory on GPU.
#ifdef _WITHGPU
  unsigned* IdpRefg; ///<Idp of selected particles [PsCount].
  double3*  PosRefg; ///<Position of selected particles [PsCount*3].
#endif

  JPartMotRefBi4Save* MotData1; ///<Main object to save motion reference data according to TimeOut as file ibi4.
  JPartMotRefBi4Save* MotData2; ///<Extra object to save motion reference data according to TimeOut2 as file ibi4.

  bool UseUnits; ///<Activate the use of Units vector.

 private:
  //void ClearPartData();
  //static std::string GetNamePart(unsigned cpart);
  void ConfigAddMk(bool floating,word mktype,unsigned begin,unsigned mknp);

  void LoadPosRefCpu(double timestep,unsigned step,unsigned np
    ,const tdouble3* pos,const unsigned* ridpmot);

#ifdef _WITHGPU
  void ConfigDataGpu();
  void LoadPosRefGpu(double timestep,unsigned step,unsigned np
    ,const double2* posxy,const double* posz,const unsigned* ridpmot);
#endif

 public:
  JDsPartMotionSave(bool cpu,std::string appname,std::string dir
    ,unsigned casenfixed,unsigned casenmoving,unsigned casenfloat
    ,word mkboundfirst,unsigned mkcount,double timeout,double timeout2);
  ~JDsPartMotionSave();
  void Reset();

  void ConfigAddMovingMk(word mktype,unsigned begin,unsigned mknp){
    ConfigAddMk(false,mktype,begin,mknp);
  }
  void ConfigAddFloatingMk(word mktype,unsigned begin,unsigned mknp){
    ConfigAddMk(true,mktype,begin,mknp);
  }
  void ConfigMotionRefs(unsigned np,const tdouble3* pos,const unsigned* idp);

  void SaveInitial();

  bool ExtraIsActive()const{ return(MotData2!=NULL); };
  std::string GetFileName(bool mainfile)const;
  
  void SaveDataMainCpu(int cpart,double timestep,unsigned step
    ,unsigned np,const tdouble3* pos,const unsigned* ridpmot);
  void AddDataExtraCpu(int cpart,double timestep,unsigned step
    ,unsigned np,const tdouble3* pos,const unsigned* ridpmot);

  void SaveDataExtra();

#ifdef _WITHGPU
  void ResetDataGpu();
  void SaveDataMainGpu(int cpart,double timestep,unsigned step
    ,unsigned np,const double2* posxy,const double* posz,const unsigned* ridpmot);
  void AddDataExtraGpu(int cpart,double timestep,unsigned step
    ,unsigned np,const double2* posxy,const double* posz,const unsigned* ridpmot);
#endif

};

#endif


