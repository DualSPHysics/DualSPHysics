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
//:# Clase JPartMotRefBi4Save para grabar la informacion de referencia para el 
//:# calculo de movimiento de floating and moving bodies.
//:# =========
//:# - Implementacion. (21-07-2023)
//:#############################################################################

/// \file JPartMotRefBi4Save.h \brief Declares the class \ref JPartMotRefBi4Save.

#ifndef _JPartMotRefBi4Save_
#define _JPartMotRefBi4Save_

#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"
#include "JPartMotionDef.h"
#include <string>
#include <vector>


//##############################################################################
//# JPartMotRefBi4Save
//##############################################################################
/// \brief Allows writing information of floating objects during simulation.

class JPartMotRefBi4Save : protected JObject
{
 private:
  const std::string AppName; ///<Nombre de aplicacion. Application Name.
  const std::string Dir;     ///<Directorio de datos. Data Directory.

  JBinaryData* BdData;   ///<Almacena la informacion general de los datos (constante para cada PART). Stores general information of data (constant for each PART).
  JBinaryData* BdPart;   ///<Pertenece a Data y almacena informacion de un part (incluyendo datos de floatings). Belongs to data and stores information of a part (including data of floatings).

  //-Variables de gestion. Management of variables.
  static const unsigned FormatVerDef=230729;  ///<Version de formato by default. Version of format by default.
  unsigned FormatVer;    ///<Version de formato. Format version.

  bool MainFile;         ///<Main output frequency.
  double TimeOut;        ///<Defined output period (for information only).
  std::string FileFull;  ///<Full name of output file.

  word MkBoundFirst;     ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  unsigned MkMovingCount;///<Number of moving bodies.
  unsigned MkFloatCount; ///<Number of floating bodies.
  unsigned MkCount;      ///<Number of moving and floating bodies (MkMovingCount+MkFloatCount).
  unsigned PsCount;      ///<Number of reference positions (PsCount=MkCount*3).

  bool InitialSaved;     ///<Indica si se grabo la informacion de cabecera. Indicates if header information is recorded.
  unsigned PartCount;    ///<Number of saved PART data.

  //-Stores data of several steps.
  unsigned  DatSize;   ///<Number of allocated data items.
  unsigned  DatCount;  ///<Number of stored data items.
  int       DatPart;   ///<Part number of data.
  double*   DatTime;   ///<Step number [DataSize]
  unsigned* DatStep;   ///<Step number [DataSize]
  tdouble3* DatPosRef; ///<Position of selected particles [DataSize*PsCount].

 private:
  void ResetBdData();
  void ResetBdPart();
  void ResizeDat(unsigned size);

  static std::string GetNamePart(unsigned cpart);
  JBinaryData* MakeBdPartSingle(int cpart,double timestep,unsigned step
    ,unsigned npos,const tdouble3* posref);
  JBinaryData* MakeBdPartArray(int cpart);
  void SaveBdPart();

 public:
  JPartMotRefBi4Save(std::string appname,std::string dir);
  ~JPartMotRefBi4Save();
  void Reset();

  llong GetAllocMemory()const;
  static std::string GetFileNameDef(bool mainfile,std::string dir="");

  void Config(bool mainfile,double timeout,word mkboundfirst
    ,unsigned mvcount,unsigned ftcount,const StMkMotionData* mkmotiondata);
  void SaveInitial();

  void SavePart(int cpart,double timestep,unsigned step,unsigned npos,const tdouble3* posref);

  void AddDataPart(int cpart,double timestep,unsigned step,unsigned npos,const tdouble3* posref);
  void SaveStoredData();
};


#endif


