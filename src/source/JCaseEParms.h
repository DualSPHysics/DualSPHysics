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
//:# Cambios:
//:# =========
//:# - Nuevos metodos LoadFileXml() y SaveFileXml() para cargar o generar un
//:#   fichero xml de forma directa. (28-11-2010)
//:# - Nuevos metodos GetValueNumInt() y GetValueNumDouble() permiten leer un
//:#   valor entre varios dentro del atributo Value. (04-12-2011)
//:# - Traduccion de comentarios al ingles. (10-02-2012)
//:# - Nuevo metodo GetValueNumStr() para leer atributos string. (10-11-2012)
//:# - Nuevo metodo GetValueDouble3() para leer atributos tdouble3. (01-10-2015)
//:# - Se muestran unidades de los parametros. (15-12-2015)
//:# - Nueva definicion del dominio de simulacion. (22-03-2019)
//:# - Objeto JXml pasado como const para operaciones de lectura. (17-03-2020)  
//:# - Comprueba que los valores enteros y reales sean validos. (17-03-2020)  
//:# - Improved exception managment. (18-03-2020)
//:# - Cambio de nombre de J.SpaceEParms a J.CaseEParms. (28-06-2020)
//:#############################################################################

/// \file JCaseEParms.h \brief Declares the class \ref JCaseEParms.

#ifndef _JCaseEParms_
#define _JCaseEParms_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"

class JXml;
class TiXmlElement;

//##############################################################################
//# JCaseEParms
//##############################################################################
/// \brief Manages the info of execution parameters from the input XML file.

class JCaseEParms : protected JObject
{
public:
  /// Structure used to store information about each parameter.
  typedef struct{
    std::string key;
    std::string value;
    std::string comment;
    std::string unitscomment;
  }JCaseEParmsItem;

  /// Domain configuration mode.
  typedef enum{ 
    DC_Fixed=0    ///< Fixed value.
   ,DC_Default=1  ///< Default.
   ,DC_DefValue=2 ///< Default +/- value.
   ,DC_DefPrc=3   ///< Default +/- % of size.
  }TpDomainMode;

  /// Structure used to store information about default domain position.
  typedef struct{
    TpDomainMode mode;
    double value;
    std::string textmod;
  }JCaseEParmsPos;

private:
  typedef std::vector<JCaseEParmsItem> VecList;
  typedef std::vector<JCaseEParmsItem>::iterator VecListIt;
  
  VecList List;
  //-Defines simulation domain.
  std::string Posminx,Posminy,Posminz;
  std::string Posmaxx,Posmaxy,Posmaxz;

  int CheckPosValue(const std::string &value,bool isposmin,JCaseEParmsPos &ps)const;
  std::string ReadPosValue(const JXml *sxml,TiXmlElement* ele,const std::string &name,const std::string &subname)const;

  JCaseEParmsItem* GetItemPointer(const std::string &key);
  std::string GetValueNum(const std::string &key,int num);
  std::string GetValue(const std::string &key);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void WriteXml(JXml *sxml,TiXmlElement* lis)const;
public:
  JCaseEParms();
  ~JCaseEParms();
  void Reset();

  void Add(const std::string &key,const std::string &value,const std::string &comment,const std::string &unitscomment="");
  void SetValue(const std::string &key,const std::string &value);
  void SetComment(const std::string &key,const std::string &comment);
  bool Exists(const std::string &key){ return(GetItemPointer(key)!=NULL); }

  void SetPosmin(std::string x,std::string y,std::string z);
  void SetPosmax(std::string x,std::string y,std::string z);
  JCaseEParmsPos GetPosminValue(char key)const;
  JCaseEParmsPos GetPosmaxValue(char key)const;
  bool IsPosDefault()const{ return(Posminx=="default" && Posminy=="default" && Posminz=="default" && Posmaxx=="default" && Posmaxy=="default" && Posmaxz=="default"); }

  int GetValueNumInt(const std::string &key,int num,bool optional=false,int valdef=0);
  double GetValueNumDouble(const std::string &key,int num,bool optional=false,double valdef=0);
  float GetValueNumFloat(const std::string &key,int num,bool optional=false,float valdef=0){ return(float(GetValueNumDouble(key,num,optional,valdef))); }
  std::string GetValueNumStr(const std::string &key,int num,bool optional=false,std::string valdef="");
  
  int GetValueInt(const std::string &key,bool optional=false,int valdef=0){ return(GetValueNumInt(key,0,optional,valdef)); }
  double GetValueDouble(const std::string &key,bool optional=false,double valdef=0){ return(GetValueNumDouble(key,0,optional,valdef)); }
  float GetValueFloat(const std::string &key,bool optional=false,float valdef=0){ return(GetValueNumFloat(key,0,optional,valdef)); }
  std::string GetValueStr(const std::string &key,bool optional=false,std::string valdef=""){ return(GetValueNumStr(key,0,optional,valdef)); }
  tdouble3 GetValueDouble3(const std::string &key,bool optional=false,tdouble3 valdef=TDouble3(0)){ return(TDouble3(GetValueNumDouble(key,0,optional,valdef.x),GetValueNumDouble(key,1,optional,valdef.y),GetValueNumDouble(key,2,optional,valdef.z))); }

  unsigned Count()const{ return(unsigned(List.size())); }
  std::string ToString(unsigned pos)const;
  JCaseEParmsItem GetParm(unsigned pos)const;
  void LoadFileXml(const std::string &file,const std::string &path);
  void SaveFileXml(const std::string &file,const std::string &path,bool newfile=true)const;
  void LoadXml(const JXml *sxml,const std::string &place);
  void SaveXml(JXml *sxml,const std::string &place)const;
};

#endif


