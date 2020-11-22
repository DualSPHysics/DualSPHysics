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

//#############################################################################
//# Cambios:
//# =========
//# - Clase para grabar movimiento y otra informacion de los floatings con 
//#   mayor frecuencia. (19-10-2020)
//#############################################################################

/// \file JFtMotionSave.h \brief Declares the class \ref JFtMotionSave.

#ifndef _JFtMotionSave_
#define _JFtMotionSave_

#include <string>
#include "JObject.h"
#include "TypesDef.h"
#include "DualSphDef.h"
#include "JPartMotionDef.h"

#ifdef _WITHGPU
#include <cuda_runtime_api.h>
#endif

class JLog2;
class JPartFloatBi4Save;

//##############################################################################
//# JFtMotionSave
//##############################################################################
/// \brief Saves floating motion data with high frequency.

class JFtMotionSave : protected JObject
{
private:
  const double TimeOut;
  JLog2 *Log;

  unsigned FtCount;      ///<Number of floating bodies.
  word MkBoundFirst;     ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  unsigned CaseFtBegin;  ///<Id of first floating boundary particles. 
  unsigned CaseNfloat;   ///<Number of floating boundary particles. 

  StMkMotionData *FtMks; ///<Information of floating boundary Mk blocks [FtCount].
  
  JPartFloatBi4Save *FtData; ///Object to save data as file fbi4.

  unsigned Num;          ///<Number of saved data.
  double NextTimeOutput; ///<Next time to save data.

  //-Auxiliary memory.
  unsigned *IdpRef;      ///<Idp of selected particles [FtCount*3].
  tdouble3 *PosRef;      ///<Position of selected particles [FtCount*3].

  //-Auxiliary memory on GPU.
  unsigned *IdpRefg;     ///<Idp of selected particles [FtCount*3].
  double   *PosRefg;     ///<Position of selected particles [FtCount*3*3].

private:
  void ConfigPosRef(unsigned np,const tdouble3 *pos,const unsigned *idp);
  double GetNextTime(double t)const;
  void SaveFtData(double timestep,unsigned nstep,const StFloatingData *ftobjs);

public:
  JFtMotionSave(double tout);
  ~JFtMotionSave();
  void Reset();

  void Config(std::string appname,std::string dirout,word mkboundfirst
    ,unsigned ftcount,const StFloatingData *ftobjs,unsigned np
    ,const tdouble3 *pos,const unsigned *idp);

  double GetTimeOut()const{ return(TimeOut); }

  bool CheckTime(double timestep)const{ return(timestep>=NextTimeOutput); }

  void SaveFtDataCpu(double timestep,unsigned nstep,const StFloatingData *ftobjs
    ,unsigned np,const tdouble3 *pos,const unsigned *ftridp);

#ifdef _WITHGPU
  void SaveFtDataGpu(double timestep,unsigned nstep,const StFloatingData *ftobjs
    ,unsigned np,const double2 *posxyg,const double *poszg,const unsigned *ftridpg);
#endif

};

#endif

