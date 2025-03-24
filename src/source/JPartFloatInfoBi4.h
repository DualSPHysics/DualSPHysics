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
//:# Clase JPartFloatInfoBi4Save para grabar la informacion de los floatings.
//:# Clase JPartFloatInfoBi4Load para recuperar la informacion de los floatings.
//:# Algunas de sus funcionalidades son:
//:# - Graba datos de cabecera y de estado de floatings por PART.
//:# - Comprueba cabecera y recupera datos de estado de floatings por PART.
//:# Cambios:
//:# =========
//:# - Implementacion de JPartFloatBi4Save y JPartFloatBi4Load. (04-12-2014)
//:# - Implementacion independiente de StFloatingData en DualSphDef.h. (11-05-2015)
//:# - Incluye informacion de MkBoundFirst (FormatVer=180423). (23-04-2018)
//:# - En JPartFloatBi4Load se comprueba que los PARTs cargados sean consecutivos. (11-08-2019)
//:# - Se guarda tambien Massp de floatings. (10-03-2020)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Graba aceleracion lineal y angular. (19-10-2020)
//:# - Graba numero de step. (19-10-2020)
//:# - Permite graba posiciones de referencia para el calculo de movimiento. (19-10-2020)
//:# - Nuevos metodos y otras mejoras. (20-10-2020)
//:# - Graba mas valoras de fuerzas y aceleraciones. (30-03-2023)
//:# - Updated for JBinaryData with 64-bit size. (20-06-2023)
//:# - Nueva implementacion simplificada en JPartFloatInfoBi4Save y JPartFloatInfoBi4Load. (23-07-2023)
//:# - Error corregido en GetIdxCpt(), GetPartFptXXX() y GetPPartFptXXX(). (24-11-2024)
//:#############################################################################

/// \file JPartFloatInfoBi4.h \brief Declares the classes \ref JPartFloatInfoBi4Save and class \ref JPartFloatInfoBi4Load.

#ifndef _JPartFloatInfoBi4_
#define _JPartFloatInfoBi4_

#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"
#include <string>
#include <vector>


//##############################################################################
//# JPartFloatInfoBi4Data
//##############################################################################
/// \brief Stores information of floating objects.

class JPartFloatInfoBi4Data : protected JObject
{
 public:
  friend class JPartFloatInfoBi4Save;
  friend class JPartFloatInfoBi4Load;

  /// Structure with basic and constant floating data.
  typedef struct {
    word mkbound;   ///<MkBound of the floating object.
    unsigned begin; ///<First particle of the floating object.
    unsigned count; ///<Number of floating objects.
    float mass;     ///<Mass of the floating object (units:Kg).
    float massp;    ///<Mass of the particle of the floating object (units:Kg).
    float radius;   ///<Maximum distance between particles and center (units:m).
  }StFtDataCte;

 public:
  const unsigned FtCount;  ///<Number of floats.
  const unsigned FptCount; ///<Number of force points.
  const word MkBoundFirst; ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).

 private:
  StFtDataCte* CteData; ///<Basic and constant floating data [FtCount].

  //-Stores data of several PARTs and several steps by PART.
  unsigned  PartSize;  ///<Number of allocated data items.
  unsigned  PartCount; ///<Number of stored data items.

  //-Data variables of PARTs [PartSize].
  int*      PartNum;
  double*   PartTime;
  unsigned* PartStep;

  //-Data variables of floating bodies [PartSize*FtCount].
  tdouble3* PartCenter;  ///<Center of the floating object [PartSize*FtCount]. 
  tfloat3*  PartVelLin;  ///<Linear velocity of the floating object (units:m/s) [PartSize*FtCount].
  tfloat3*  PartVelAng;  ///<Angular velocity of the floating object (units:rad/s) [PartSize*FtCount].
  tfloat3*  PartAceLin;  ///<Linear acceleration of the floating object (units:m/s^2) [PartSize*FtCount].
  tfloat3*  PartAceAng;  ///<Angular acceleration of the floating object (units:rad/s^2) [PartSize*FtCount].
  //-Extra data [PartSize*FtCount].
  tfloat3*  PartExtForceLin;///<External linear forces (moorings and imposed forces) (units:N) [PartSize*FtCount].
  tfloat3*  PartExtForceAng;///<External angular forces (moorings and imposed forces) (units:N*m*rad) [PartSize*FtCount].
  tfloat3*  PartFluForceLin;///<Linear forces from fluid (sum in eq.48 at Dominguez et al 2022) (units:N) [PartSize*FtCount].
  tfloat3*  PartFluForceAng;///<Angular forces from fluid (sum in eq.49 at Dominguez et al 2022) (units:N*m*rad) [PartSize*FtCount].
  tfloat3*  PartPreAceLin;  ///<Linear acceleration before constraints (includes external forces and gravity) (units:m/s^2) [PartSize*FtCount].
  tfloat3*  PartPreAceAng;  ///<Angular acceleration before constraints (multiplied by rotated inertia tensor) (units:rad/s^2) [PartSize*FtCount].

  //-Data of force points [PartSize*FptCount].
  word*     FptMkbound;  ///<MkBound of floating body. [PartSize*FptCount].
  tdouble3* FptPos;      ///<Position. [PartSize*FptCount].
  tfloat3*  FptForce;    ///<Force. [PartSize*FptCount].

 private:
  void ClearFtData(bool clearctedata);

 public:
  JPartFloatInfoBi4Data(unsigned ftcount,unsigned fptcount,word mkboundfirst);
  ~JPartFloatInfoBi4Data();
  void Reset();

  llong GetAllocMemory()const;

  void ResizeData(unsigned size);

  void SetCteData(unsigned cf,word mkbound,unsigned begin
    ,unsigned count,float mass,float massp,float radius);
  void LoadCteData(const JPartFloatInfoBi4Data* ftdata);

  void SetPartData0(const JPartFloatInfoBi4Data* ftdata);
  void SetPartData0(unsigned cf,const tdouble3& center,const tfloat3& fvellin
    ,const tfloat3& fvelang,const tfloat3& facelin,const tfloat3& faceang
    ,const tfloat3& extforcelin,const tfloat3& extforceang
    ,const tfloat3& fluforcelin,const tfloat3& fluforceang
    ,const tfloat3& preacelin,const tfloat3& preaceang);
  void SetForcePoints0(unsigned npt,const word* mkbound
    ,const tdouble3* pos,const tfloat3* force);
  void SetTimeData0(int cpart,double timestep,unsigned step);

  void ClearData(){ PartCount=0; }
  void AddData(const JPartFloatInfoBi4Data* ftdata);

  unsigned GetFtCount()const{ return(FtCount); }
  unsigned GetPartCount()const{ return(PartCount); }

  unsigned GetIdxPart(unsigned cpart,bool unique)const;

  unsigned GetCf(unsigned cf)const;   ///<Use index of floating.
  unsigned GetIdx(unsigned idx)const; ///<Use index of PARTs.
  unsigned GetIdxCf(unsigned idx,unsigned cf)const{ return(GetIdx(idx)*FtCount+GetCf(cf)); }

  StFtDataCte GetHead(unsigned cf)const{ return(CteData[GetCf(cf)]); }
  word     GetHeadMkbound(unsigned cf)const{ return(GetHead(cf).mkbound); }
  word     GetHeadMk     (unsigned cf)const{ return(MkBoundFirst+GetHeadMkbound(cf)); }
  unsigned GetHeadBegin  (unsigned cf)const{ return(GetHead(cf).begin ); }
  unsigned GetHeadCount  (unsigned cf)const{ return(GetHead(cf).count ); }
  float    GetHeadMass   (unsigned cf)const{ return(GetHead(cf).mass  ); }
  float    GetHeadMassp  (unsigned cf)const{ return(GetHead(cf).massp ); }
  float    GetHeadRadius (unsigned cf)const{ return(GetHead(cf).radius); }

  unsigned GetPartCpart   (unsigned idx)const{ return(PartNum [GetIdx(idx)]); }
  unsigned GetPartStep    (unsigned idx)const{ return(PartStep[GetIdx(idx)]); }
  double   GetPartTimeStep(unsigned idx)const{ return(PartTime[GetIdx(idx)]); }

  tdouble3 GetPartCenter(unsigned idx,unsigned cf)const{ return(PartCenter[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartVelLin(unsigned idx,unsigned cf)const{ return(PartVelLin[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartVelAng(unsigned idx,unsigned cf)const{ return(PartVelAng[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartAceLin(unsigned idx,unsigned cf)const{ return(PartAceLin[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartAceAng(unsigned idx,unsigned cf)const{ return(PartAceAng[GetIdxCf(idx,cf)]); }

  tfloat3  GetPartExtForceLin(unsigned idx,unsigned cf)const{ return(PartExtForceLin[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartExtForceAng(unsigned idx,unsigned cf)const{ return(PartExtForceAng[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartFluForceLin(unsigned idx,unsigned cf)const{ return(PartFluForceLin[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartFluForceAng(unsigned idx,unsigned cf)const{ return(PartFluForceAng[GetIdxCf(idx,cf)]); }
  tfloat3  GetPartPreAceLin  (unsigned idx,unsigned cf)const{ return(PartPreAceLin  [GetIdxCf(idx,cf)]); }
  tfloat3  GetPartPreAceAng  (unsigned idx,unsigned cf)const{ return(PartPreAceAng  [GetIdxCf(idx,cf)]); }


  const tdouble3* GetPPartCenter(unsigned idx)const{ return(PartCenter+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartVelLin(unsigned idx)const{ return(PartVelLin+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartVelAng(unsigned idx)const{ return(PartVelAng+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartAceLin(unsigned idx)const{ return(PartAceLin+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartAceAng(unsigned idx)const{ return(PartAceAng+GetIdxCf(idx,0)); }

  const tfloat3*  GetPPartExtForceLin(unsigned idx)const{ return(PartExtForceLin+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartExtForceAng(unsigned idx)const{ return(PartExtForceAng+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartFluForceLin(unsigned idx)const{ return(PartFluForceLin+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartFluForceAng(unsigned idx)const{ return(PartFluForceAng+GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartPreAceLin  (unsigned idx)const{ return(PartPreAceLin  +GetIdxCf(idx,0)); }
  const tfloat3*  GetPPartPreAceAng  (unsigned idx)const{ return(PartPreAceAng  +GetIdxCf(idx,0)); }


  unsigned GetFptCount()const{ return(FptCount); }
  unsigned GetCpt(unsigned cpt)const;
  unsigned GetIdxCpt(unsigned idx,unsigned cpt)const{ return(GetIdx(idx)*FptCount+GetCpt(cpt)); }

  word     GetPartFptMkbound(unsigned idx,unsigned cpt)const{ return(FptMkbound[GetIdxCpt(idx,cpt)]); }
  tdouble3 GetPartFptPos    (unsigned idx,unsigned cpt)const{ return(FptPos    [GetIdxCpt(idx,cpt)]); }
  tfloat3  GetPartFptForce  (unsigned idx,unsigned cpt)const{ return(FptForce  [GetIdxCpt(idx,cpt)]); }


  const word*     GetPPartFptMkbound(unsigned idx)const{ return(FptMkbound+GetIdxCpt(idx,0)); }
  const tdouble3* GetPPartFptPos    (unsigned idx)const{ return(FptPos    +GetIdxCpt(idx,0)); }
  const tfloat3*  GetPPartFptForce  (unsigned idx)const{ return(FptForce  +GetIdxCpt(idx,0)); }

};


//##############################################################################
//# JPartFloatInfoBi4Save
//##############################################################################
/// \brief Allows writing information of floating objects during simulation.

class JPartFloatInfoBi4Save : protected JObject
{
 private:
  const std::string AppName; ///<Nombre de aplicacion. Application Name.
  const std::string Dir;     ///<Directorio de datos. Data Directory.

  JBinaryData* BdData;      ///<Almacena la informacion general de los datos (constante para cada PART). Stores general information of data (constant for each PART).
  JBinaryData* BdPart;      ///<Pertenece a Data y almacena informacion de un part (incluyendo datos de floatings). Belongs to data and stores information of a part (including data of floatings).

  //-Variables de gestion. Management of variables.
  static const unsigned FormatVerDef=230729;  ///<Version de formato by default. Version of format by default.
  unsigned FormatVer;    ///<Version de formato. Format version.

  bool MainFile;         ///<Main output frequency.
  double TimeOut;        ///<Defined output period (for information only).
  std::string FileFull;  ///<Full name of output file.

  word MkBoundFirst;     ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  unsigned FtCount;      ///<Numero de floatings. Number of floats.
  unsigned FptCount;     ///<Number of force points.

  bool InitialSaved;     ///<Indica si se grabo la informacion de cabecera. Indicates if header information is recorded.
  unsigned PartCount;    ///<Number of saved PART data.

  JPartFloatInfoBi4Data* FtData; ///<Store floating data.

 private:
  void ResetBdData();
  void ResetBdPart();
  static std::string GetNamePart(unsigned cpart);
  JBinaryData* MakeBdPartSingle(int cpart,double timestep,unsigned step
    ,const JPartFloatInfoBi4Data* ftdata);
  JBinaryData* MakeBdPartArray(const JPartFloatInfoBi4Data* ftdata);
  void SaveBdPart();

 public:
  JPartFloatInfoBi4Save(std::string appname,std::string dir);
  ~JPartFloatInfoBi4Save();
  void Reset();

  llong GetAllocMemory()const;
  static std::string GetFileNameDef(bool mainfile,std::string dir="");

  void Config(bool mainfile,double timeout,const JPartFloatInfoBi4Data* ftdata);
  void SaveInitial();

  void SavePart(int cpart,unsigned step,const JPartFloatInfoBi4Data* ftdata);
  
  void AddDataPart(int cpart,unsigned step,const JPartFloatInfoBi4Data* ftdata);
  void SaveStoredData();
};

//##############################################################################
//# JPartFloatInfoBi4Load
//##############################################################################
/// \brief Allows reading information of floating objects saved during simulation.

class JPartFloatInfoBi4Load : protected JObject
{
 private:
  static const unsigned FormatVerDef=230729;  ///<Version de formato by default. Version of format by default.
  unsigned FormatVer;    ///<Version de formato. Format version.
   
  std::string FileData;   ///<Name of input file with data.

  JPartFloatInfoBi4Data* FtData; ///<Store floating data.

  //-Loaded data in FtData.
  word MkBoundFirst;     ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  unsigned FtCount;      ///<Numero de floatings. Number of floats.
  unsigned FptCount;     ///<Number of force points.
  unsigned PartCount;    ///<Number of PARTs.
  unsigned FirstPart;    ///<Primer PART almacenado. First number of stored PART.

 private:
  const word*     GetArrayWord   (JBinaryData* bd,const std::string& name,unsigned count);
  const unsigned* GetArrayUint   (JBinaryData* bd,const std::string& name,unsigned count);
  const float*    GetArrayFloat  (JBinaryData* bd,const std::string& name,unsigned count);
  const double*   GetArrayDouble (JBinaryData* bd,const std::string& name,unsigned count);
  const tfloat3*  GetArrayFloat3 (JBinaryData* bd,const std::string& name,unsigned count);
  const tdouble3* GetArrayDouble3(JBinaryData* bd,const std::string& name,unsigned count);

 public:
  JPartFloatInfoBi4Load();
  ~JPartFloatInfoBi4Load();
  void Reset();

  static std::string GetFileNameDef(bool mainfile,std::string dir="");
  void LoadFile(std::string filename);

  void CheckHeadData(unsigned cf,word mkbound,unsigned begin
    ,unsigned count,float mass,float massp);

  word GetMkBoundFirst()const{ return(MkBoundFirst); }

  unsigned GetFtCount()const{ return(FtCount); }
  unsigned GetCount()const{ return(PartCount); }
  unsigned GetFirstPart()const{ return(FirstPart); }

  const JPartFloatInfoBi4Data* GetFtData()const{ return(FtData); }
};

#endif


