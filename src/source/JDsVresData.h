//HEAD_DSPH
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

/// \file JDsVresData.h \brief Declares the class \ref JDsExtraData.

#ifndef _JDsVresData_
#define _JDsVresData_

#include "JObject.h"
#include "TypesDef.h"

class JLog2;
class JBinaryData;
class JRangeFilter;

//##############################################################################
//# JDsVResDataSave
//##############################################################################
/// \brief Saves extra data for restart option.

class JDsVResDataSave : protected JObject
{
protected:
  static const unsigned FormatVerDef=211030;   ///<Version de formato by default.
  JLog2* Log;

  const std::string AppName;   ///<Nombre de aplicacion. Application Name.
  const std::string Dir;       ///<Directorio de datos. Data Directory.
  unsigned CaseNbound;   ///<Number of boundary particles ( \ref Nfixed + \ref Nmoving + \ref Nfloat ).
  unsigned CaseNfloat;   ///<Number of floating boundary particles. 

  int SvParts; ///<Part interval (or list) for saving extra data for restart option (default=empty=disabled)
  JRangeFilter* FilterParts;

  int Cpart;
  JBinaryData* Data;

public:
  JDsVResDataSave(std::string appname,std::string dir,JLog2* log);
  ~JDsVResDataSave();
  void Reset();
  static std::string GetFileNamePart(int cpart);

  void Config(std::string svparts);

  bool CheckSave(int cpart)const;

  void InitPartData(int cpart,double timestep,int nstep);
  void AddArray(unsigned np,unsigned msize,const tdouble3* pos,const tfloat3* normals
  ,const tfloat3* velmot,const float* mass,const double* matarray);
  void SavePartData();
};


//##############################################################################
//# JDsVResDataLoad
//##############################################################################
/// \brief Loads extra data for restart option.

class JDsVResDataLoad : protected JObject
{
protected:
  static const unsigned FormatVerDef=211030;   ///<Version de formato by default.
  JLog2* Log;


  std::string Dir;   ///<Data Directory.
  int Cpart;
  std::string FileData;
  JBinaryData* Data;

public:
  JDsVResDataLoad(JLog2* log);
  ~JDsVResDataLoad();
  void Reset();
  static std::string GetFileNamePart(std::string dir,int cpart);

  static bool ExistsPartData(std::string dir,int cpart);
  void LoadPartData(std::string dir,int cpart);

  void LoadArray(unsigned np,unsigned msize,tfloat3* velmot,float* mass,double* matarray);
  tmatrix4d GetMat();
};

#endif


