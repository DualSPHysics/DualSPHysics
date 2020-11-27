//HEAD_DSPH
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

/// \file JLinearValue.h \brief Declares the class \ref JLinearValue.

#ifndef _JLinearValue_
#define _JLinearValue_

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase obtener valores interpolados linealmente de una lista. (05-01-2014)
//:# - Nuevo metodo LoadFile() para cargar valores de un fichero. (20-01-2014)
//:# - Se eliminaron includes innecesarios. (01-05-2014)
//:# - La funcion GetNewInterval() indica cuando se cambia de intervalo. (25-01-2017)
//:# - Puede gestionar varios valores para cada instante. (17-04-2017)
//:# - Nuevos metodos para calcular interpolacion de forma externa. (17-04-2017)
//:# - Nuevos metodos SetTimeValue() para modificar datos de la lista. (24-01-2018)
//:# - Nuevo constructor con input file. (13-02-2019)
//:# - Nuevo constructor de copias. (14-02-2019)
//:# - Permite uso de DBL_MAX como valor especial. (21-02-2019)
//:# - El uso de valores especiales se establece en el constructor. (25-02-2019)
//:# - Permite cargar datos directamente del XML. (26-02-2019)
//:# - Amplia SIZEMAX por defecto de 200k a 20M. (21-10-2019)
//:# - Error corregido cuando se cargaban valores none del XML. (08-02-2020)
//:# - Objeto JXml pasado como const para operaciones de lectura. (17-03-2020)  
//:# - Mejora la gestion de excepciones. (23-04-2020)
//:# - Nuevos metodos GetValue3ByIdx() y GetValue3d3d(). (23-04-2020)
//:# - Con la opcion OptionalValues los valores que falten del XML se leen como 
//:#   cero. No es compatible con la opcion SpecialValues. (24-04-2020)
//:# - Nuevos metodos RnSetValues(). (05-09-2020)
//:#############################################################################

#include "JObject.h"
#include "TypesDef.h"
#include "JLinearValueDef.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>

class JReadDatafile;
#ifdef JLinearValue_UseJXml
class JXml;
class TiXmlElement;
#endif

//##############################################################################
//# JLinearValue
//##############################################################################
/// \brief Implements a list of values according time using linear interpolation.

class JLinearValue : protected JObject
{
protected:
  static const unsigned SIZEMAX=20000000;
  static const unsigned SIZEINITIAL=500;

  std::string File;
  unsigned Size;
  unsigned Count;
  double *Times;
  double *Values;

  bool NewInterval;

  //-Variables with last search.
  double TimeStep;
  unsigned Position;
  unsigned PositionNext;
  double TimePre;
  double TimeNext;
  double TimeFactor;

public:
  const unsigned Nvalues;
  const bool SpecialValues;  //<Uses the special value DBL_MAX. Missing values in XML configuration are considered as DBL_MAX.
  const bool OptionalValues; //<Only in ReadXmlValues() when SpecialValues=false, missing values are considered as zero. 

  JLinearValue(unsigned nvalues=1,bool specialvalues=false,bool optionalvalues=false);
  JLinearValue(const std::string &inputfile,unsigned nvalues=1,bool specialvalues=false,bool optionalvalues=false);
  JLinearValue(const JLinearValue &obj);
  ~JLinearValue();
  void Reset();
  void CopyFrom(const JLinearValue &obj);
  unsigned GetAllocMemory()const;

  void SetSize(unsigned size);
  unsigned GetSize()const{ return(Size); }

  unsigned AddTimeValue(double time,double value);
  unsigned AddTimeValue(double time,double value,double value2){ const unsigned idx=AddTimeValue(time,value); SetValue(idx,1,value2); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3){ const unsigned idx=AddTimeValue(time,value,value2); SetValue(idx,2,value3); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3,double value4){ const unsigned idx=AddTimeValue(time,value,value2,value3); SetValue(idx,3,value4); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3,double value4,double value5){ const unsigned idx=AddTimeValue(time,value,value2,value3,value4); SetValue(idx,4,value5); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3,double value4,double value5,double value6){ const unsigned idx=AddTimeValue(time,value,value2,value3,value4,value5); SetValue(idx,5,value6); return(idx); }
  void SetValue(unsigned idx,unsigned cvalue,double value);

  void SetTimeValue(unsigned idx,double time,double value);
  void SetTimeValue(unsigned idx,double time,double value,double value2){  SetTimeValue(idx,time,value); SetValue(idx,1,value2);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3){  SetTimeValue(idx,time,value,value2); SetValue(idx,2,value3);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3,double value4){  SetTimeValue(idx,time,value,value2,value3); SetValue(idx,3,value4);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3,double value4,double value5){  SetTimeValue(idx,time,value,value2,value3,value4); SetValue(idx,4,value5);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3,double value4,double value5,double value6){  SetTimeValue(idx,time,value,value2,value3,value4,value5); SetValue(idx,5,value6);  }

  double ReadNextDouble(JReadDatafile &rdat,bool in_line=false);
  void LoadFile(std::string file);
  std::string GetFile()const{ return(File); };

  unsigned GetCount()const{ return(Count); }
  void VisuData();

  //-For simple use.
  double GetValue(double timestep,unsigned cvalue=0);
  float GetValuef(double timestep,unsigned cvalue=0);
  tdouble3 GetValue3d(double timestep);
  tfloat3  GetValue3f(double timestep);
  void     GetValue3d3d(double timestep,tdouble3 &v1,tdouble3 &v2);

  bool GetNewInterval()const{ return(NewInterval); }

  //-For external use.
  void FindTime(double timestep);
  unsigned GetPos()const{ return(Position); };
  unsigned GetPosNext()const{ return(PositionNext); };
  double GetTimeByIdx(unsigned idx)const;
  double GetValueByIdx(unsigned idx,unsigned cvalue=0)const;
  tdouble3 GetValue3ByIdx(unsigned idx,unsigned cvalue=0)const;

  //-Loads or saves on XML file.
#ifdef JLinearValue_UseJXml
  void ReadXmlValues(const JXml *sxml,TiXmlElement* ele,std::string name
    ,std::string subname,std::string attributes);
  TiXmlElement* WriteXmlValues(JXml *sxml,TiXmlElement* ele,std::string name
    ,std::string subname,std::string attributes)const;
#endif

  void RnSetValues(double v){ RnSetValues(0,v,0,v); }
  void RnSetValues(double t0,double v0,double t1,double v1);
};

#endif


