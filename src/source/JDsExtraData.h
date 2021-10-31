//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2021 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsExtraData.h \brief Declares the class \ref JDsExtraData.

#ifndef _JDsExtraData_
#define _JDsExtraData_

#include "JObject.h"
#include "TypesDef.h"
#include "DualSphDef.h"

class JLog2;
class JBinaryData;
class JRangeFilter;

//##############################################################################
//# JDsExtraDataSave
//##############################################################################
/// \brief Saves extra data for restart option.

class JDsExtraDataSave : protected JObject
{
protected:
  static const unsigned FormatVerDef=211030;   ///<Version de formato by default.
  JLog2* Log;

  const std::string AppName;   ///<Nombre de aplicacion. Application Name.
  const std::string Dir;       ///<Directorio de datos. Data Directory.
  const unsigned CaseNbound;   ///<Number of boundary particles ( \ref Nfixed + \ref Nmoving + \ref Nfloat ).
  const unsigned CaseNfloat;   ///<Number of floating boundary particles. 

  int SvParts; ///<Part interval (or list) for saving extra data for restart option (default=empty=disabled)
  JRangeFilter *FilterParts;

  int Cpart;
  JBinaryData *Data;

public:
  JDsExtraDataSave(std::string appname,std::string dir
    ,unsigned casenbound,unsigned casenfloat,JLog2* log);
  ~JDsExtraDataSave();
  void Reset();
  static std::string GetFileNamePart(int cpart);

  void Config(std::string svparts);

  bool CheckSave(int cpart)const;

  void InitPartData(int cpart,double timestep,int nstep);
  void AddNormals(bool usenormalsft,unsigned np,unsigned npb
    ,const unsigned *idp,const typecode *code,const tfloat3 *boundnormal);
  void SavePartData();
};


//##############################################################################
//# JDsExtraDataLoad
//##############################################################################
/// \brief Loads extra data for restart option.

class JDsExtraDataLoad : protected JObject
{
protected:
  static const unsigned FormatVerDef=211030;   ///<Version de formato by default.
  JLog2* Log;

  const unsigned CaseNbound;   ///<Number of boundary particles ( \ref Nfixed + \ref Nmoving + \ref Nfloat ).
  const unsigned CaseNfloat;   ///<Number of floating boundary particles. 

  std::string Dir;   ///<Data Directory.
  int Cpart;
  std::string FileData;
  JBinaryData *Data;

public:
  JDsExtraDataLoad(unsigned casenbound,unsigned casenfloat,JLog2* log);
  ~JDsExtraDataLoad();
  void Reset();
  static std::string GetFileNamePart(std::string dir,int cpart);

  static bool ExistsPartData(std::string dir,int cpart);
  void LoadPartData(std::string dir,int cpart);

  bool LoadNormals(unsigned np,unsigned npb,const unsigned *idp,tfloat3 *boundnormal);

};

#endif


