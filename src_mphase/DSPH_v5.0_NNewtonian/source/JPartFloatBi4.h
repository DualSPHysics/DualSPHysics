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
//:# Descripcion:
//:# =============
//:# Clase JPartFloatBi4Save para grabar la informacion de los floatings.
//:# Clase JPartFloatBi4Load para recuperar la informacion de los floatings.
//:# Algunas de sus funcionalidades son:
//:# - Graba datos de cabecera y de estado de floatings por PART.
//:# - Comprueba cabecera y recupera datos de estado de floatings por PART.
//:# Cambios:
//:# =========
//:# - Implementacion. (04-12-2014)
//:# - Implementacion independiente de StFloatingData en DualSphDef.h. (11-05-2015)
//:# - Incluye informacion de MkBoundFirst (FormatVer=180423). (23-04-2018)
//:# - En JPartFloatBi4Load se comprueba que los PARTs cargados sean consecutivos. (11-08-2019)
//:# - Se guarda tambien Massp de floatings. (10-03-2020)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Graba aceleracion lineal y angular. (19-10-2020)
//:# - Graba numero de step. (19-10-2020)
//:# - Permite graba posiciones de referencia para el calculo de movimiento. (19-10-2020)
//:# - Nuevos metodos y otras mejoras. (20-10-2020)
//:#############################################################################

/// \file JPartFloatBi4.h \brief Declares the classes \ref JPartFloatBi4Save and class \ref JPartFloatBi4Load.

#ifndef _JPartFloatBi4_
#define _JPartFloatBi4_

#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"
#include <string>
#include <vector>
#include <fstream>
#include <climits>


//##############################################################################
//# JPartFloatBi4Save
//##############################################################################
/// \brief Allows writing information of floating objects during simulation.

class JPartFloatBi4Save : protected JObject
{
 private:
  JBinaryData *Data;      ///<Almacena la informacion general de los datos (constante para cada PART). Stores general information of data (constant for each PART).
  JBinaryData *Part;      ///<Pertenece a Data y almacena informacion de un part (incluyendo datos de floatings). Belongs to data and stores information of a part (including data of floatings).

  //-Variables de gestion. Management of variables.
  static const unsigned FormatVerDef=180423;    ///<Version de formato by default. Version of format by default.
  unsigned FormatVer;    ///<Version de formato. Format version.

  bool InitialSaved;     ///<Indica si se grabo la informacion de cabecera. Indicates if header information is recorded.

  std::string AppName;   ///<Nombre de aplicacion. Application Name.
  std::string Dir;       ///<Directorio de datos. Data Directory.
  std::string FileFull;  ///<Full name of output file.
  word MkBoundFirst;     ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  unsigned FtCount;      ///<Numero de floatings. Number of floats.
  bool PosRefData;       ///<Saves position reference data.

  //-Datos constantes de floatings (head). Constant data of floats (head).
  word *HeadMkbound;
  unsigned *HeadBegin;
  unsigned *HeadCount;
  float *HeadMass;
  float *HeadMassp;
  float *HeadRadius;

  //-Datos variables de floatings (PARTs). Data variables of floats (parts).
  tdouble3 *PartCenter; ///<Center of the floating object [FtCount]. 
  tdouble3 *PartPosRef; ///<Reference positions of the floating object [FtCount*3].. 
  tfloat3  *PartVelLin; ///<Linear velocity of the floating object (units:m/s) [FtCount].
  tfloat3  *PartVelAng; ///<Angular velocity of the floating object (units:rad/s) [FtCount].
  tfloat3  *PartAceLin; ///<Linear acceleration of the floating object (units:m/s^2) [FtCount].
  tfloat3  *PartAceAng; ///<Angular acceleration of the floating object (units:rad/s^2) [FtCount].
  
  unsigned Cpart;    ///<Numero de PART. PART number.

 private:
  void ResizeFtData(unsigned ftcount);
  void ClearPartData();
  static std::string GetNamePart(unsigned cpart);

 public:
  JPartFloatBi4Save();
  ~JPartFloatBi4Save();
  void Reset();
  void ResetData();
  void ResetPart();

  long long GetAllocMemory()const;
  static std::string GetFileNamePart();

  //Grabacion de datos:
  //Data recording:
  //====================
  //-Configuracion de objeto. Object Configuration.
  void Config(std::string appname,const std::string &dir,word mkboundfirst
    ,unsigned ftcount,bool saveposref);
  void AddHeadData(unsigned cf,word mkbound,unsigned begin,unsigned count,float mass,float massp,float radius);
  void SaveInitial(std::string file="");

  //-Configuracion de parts.  Parts Configuration.
  void AddPartData(unsigned cf,const tdouble3 &center,const tfloat3 &fvellin
    ,const tfloat3 &fvelang,const tfloat3 &facelin,const tfloat3 &faceang);
  void AddPartDataPosRef(unsigned ftcount,const tdouble3 *posref);
  JBinaryData* AddPartFloat(unsigned cpart,unsigned step,double timestep,double demdtforce);

  ////-Grabacion de fichero. File recording.
  void SavePartFloat();
  void SavePartFloat(unsigned cpart,unsigned step,double timestep,double demdtforce){   
    AddPartFloat(cpart,step,timestep,demdtforce); SavePartFloat();
  }
};


//##############################################################################
//# JPartFloatBi4Load
//##############################################################################
/// \brief Allows reading information of floating objects saved during simulation.

class JPartFloatBi4Load : protected JObject
{
 private:
  static const unsigned FormatVerDef=180423;    ///<Version de formato by default. Version of format by default.
  unsigned FormatVer;    ///<Version de formato. Format version.
   
  std::string FileData;   ///<Name of input file with data.
  JBinaryData *Data;      ///<Almacena la informacion general de los datos (constante para cada PART). Stores general information of data (constant for each PART).

  word MkBoundFirst;      ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  bool PosRefData;        ///<Position reference data is available.

  unsigned FtCount;       ///<Numero de floatings. Floating number
  unsigned PartCount;     ///<Numero de PARTs. PARTs number
  JBinaryData *Part;      ///<Pertenece a Data y almacena informacion de un part (incluyendo datos de floatings). Belongs to data and stores information of a part (including data of floatings).
  unsigned FirstPart;     ///<Primer PART almacenado. First number of stored PART.

  //-Datos constantes de floatings (head). Constant data of floatings (head).
  word *HeadMkbound;
  unsigned *HeadBegin;
  unsigned *HeadCount;
  float *HeadMass;
  float *HeadMassp;
  float *HeadRadius;

  //-Informacion de PART. PART information.
  unsigned Cpart;
  unsigned Step;
  double TimeStep;
  double DemDtForce;
  bool AceData;         ///<Linear and angular acceleration is available.  
  //-Datos variables de floatings (PARTs). Data variables of floatings (parts).
  tdouble3 *PartCenter;
  tdouble3 *PartPosRef; ///<Reference positions of the floating object [FtCount*3]. 
  tfloat3  *PartVelLin; ///<Linear velocity of the floating object (units:m/s) [FtCount].
  tfloat3  *PartVelAng; ///<Angular velocity of the floating object (units:rad/s) [FtCount].
  tfloat3  *PartAceLin; ///<Linear acceleration of the floating object (units:m/s^2) [FtCount].
  tfloat3  *PartAceAng; ///<Angular acceleration of the floating object (units:rad/s^2) [FtCount].

 private:
  JBinaryDataArray* CheckArray(JBinaryData *bd,const std::string &name
    ,JBinaryDataDef::TpData type,unsigned count=UINT_MAX);
  void ResetPart();
  void ResizeFtData(unsigned ftcount);
  void CheckPartList()const;
  void CheckPart()const;
  void CheckFloating(unsigned cf)const;

 public:
  JPartFloatBi4Load();
  ~JPartFloatBi4Load();
  void Reset();
  static std::string GetFileNamePart();

  static std::string GetLoadFile(const std::string &dir,std::string filename="");
  void LoadFile(const std::string &dir,std::string filename="");
  void CheckHeadData(unsigned cf,word mkbound,unsigned begin,unsigned count,float mass,float massp);

  word GetMkBoundFirst()const{ return(MkBoundFirst); }
  bool GetPosRefData()const{ return(PosRefData); }

  unsigned GetFtCount()const{ return(FtCount); }
  unsigned GetCount()const{ return(PartCount); }
  unsigned GetFirstPart()const{ return(FirstPart); }

  word     GetHeadMkbound(unsigned cf)const{ CheckFloating(cf); return(HeadMkbound[cf]); }
  word     GetHeadMk     (unsigned cf)const{ return(MkBoundFirst+GetHeadMkbound(cf)); }
  unsigned GetHeadBegin  (unsigned cf)const{ CheckFloating(cf); return(HeadBegin [cf]); }
  unsigned GetHeadCount  (unsigned cf)const{ CheckFloating(cf); return(HeadCount [cf]); }
  float    GetHeadMass   (unsigned cf)const{ CheckFloating(cf); return(HeadMass  [cf]); }
  float    GetHeadMassp  (unsigned cf)const{ CheckFloating(cf); return(HeadMassp [cf]); }
  float    GetHeadRadius (unsigned cf)const{ CheckFloating(cf); return(HeadRadius[cf]); }

  void LoadPartItem(unsigned cp);
  void LoadPart(unsigned cpart);

  unsigned GetPartCpart     ()const{ CheckPart(); return(Cpart); }
  unsigned GetPartStep      ()const{ CheckPart(); return(Step); }
  double   GetPartTimeStep  ()const{ CheckPart(); return(TimeStep); }
  double   GetPartDemDtForce()const{ CheckPart(); return(DemDtForce); }
  bool     GetPartAceData   ()const{ return(AceData); }
  
  tdouble3 GetPartCenter(unsigned cf)const{ CheckFloating(cf); return(PartCenter[cf]); }
  tfloat3  GetPartVelLin(unsigned cf)const{ CheckFloating(cf); return(PartVelLin[cf]); }
  tfloat3  GetPartVelAng(unsigned cf)const{ CheckFloating(cf); return(PartVelAng[cf]); }
  tfloat3  GetPartAceLin(unsigned cf)const{ CheckFloating(cf); return(AceData? PartAceLin[cf]: TFloat3(0)); }
  tfloat3  GetPartAceAng(unsigned cf)const{ CheckFloating(cf); return(AceData? PartAceAng[cf]: TFloat3(0)); }

  const tdouble3* GetPartCenter()const{ CheckPart(); return(PartCenter); }
  const tfloat3*  GetPartVelLin()const{ CheckPart(); return(PartVelLin); }
  const tfloat3*  GetPartVelAng()const{ CheckPart(); return(PartVelAng); }
  const tfloat3*  GetPartAceLin()const{ CheckPart(); return(AceData? PartAceLin: NULL); }
  const tfloat3*  GetPartAceAng()const{ CheckPart(); return(AceData? PartAceAng: NULL); }
  const tdouble3* GetPartPosRef()const{ CheckPart(); return(PosRefData? PartPosRef: NULL); }
};


#endif


