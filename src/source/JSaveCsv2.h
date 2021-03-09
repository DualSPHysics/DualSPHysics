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

//:NO_COMENTARIO
//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para Simplificar la creacion de ficheros CSV en debug. (23-03-2014)
//:# - Ahora al ejecutar Save() vacia el buffer. (06-02-2015)
//:# - Funciones AddValue() para datos triples. (26-05-2016)
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:# - New attribute CsvSepComa to configure separator in CSV files. (24-10-2017)
//:# - New improved code to simplify CSV file creation. (08-11-2017)
//:# - Error corregido: Reescribia fichero cuando ya estaba cerrado pos SaveData(). (02-03-2018)
//:# - Opcion para activar o no el uso de SetSeparators((). (12-03-2018)
//:# - Control de excepciones en el destructor. (21-03-2018)
//:# - Nuevo metodo GetAppendMode() para saber si se va ampliar un fichero existente. (06-07-2018)
//:# - Error corregido en SaveData() por el que siempre se grababa el head. (13-08-2018)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:# - Nuevo metodo AddHead3() y clase Head3 para cabeceras de datos triples. (02-11-2020)
//:# - Soporta tipo tmatrix4d. (04-02-2021)
//:#############################################################################

/// \file JSaveCsv2.h \brief Declares the class \ref JSaveCsv2.

#ifndef _JSaveCsv2_
#define _JSaveCsv2_

#include "TypesDef.h"
#include "JObject.h"

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>

/// Implements a set of classes to create CSV files.
namespace jcsv{

///Output formats.
typedef enum{ 
  TpSigned1=0,
  TpSigned2=1,
  TpSigned3=2,
  TpSigned4=3,
  TpUnsigned1=4,
  TpUnsigned2=5,
  TpUnsigned3=6,
  TpUnsigned4=7,
  TpLlong1=8,
  TpUllong1=9,
  TpFloat1=10,
  TpFloat2=11,
  TpFloat3=12,
  TpFloat4=13,
  TpDouble1=14,
  TpDouble2=15,
  TpDouble3=16,
  TpDouble4=17
}TpFormat; 


//##############################################################################
//# Fmt
//##############################################################################
/// \brief Class to set format for output (operator use).

class Fmt
{
public:
  const TpFormat TypeFmt;
  const std::string Format;
  Fmt(TpFormat typefmt):TypeFmt(typefmt),Format(""){};
  Fmt(TpFormat typefmt,const std::string &format):TypeFmt(typefmt),Format(format){};
};


//##############################################################################
//# Sep
//##############################################################################
/// \brief Class to add one or several field separators (operator use).

class Sep{
public:
  const unsigned Count;
  Sep(unsigned count=1):Count(count){};
};

//##############################################################################
//# AutoSepOn
//##############################################################################
/// \brief Class to enable auto separator change according configuration.

class AutoSepOn{
public:
  AutoSepOn(){};
};

//##############################################################################
//# AutoSepOff
//##############################################################################
/// \brief Class to disables auto separator change according configuration.

class AutoSepOff{
public:
  AutoSepOff(){};
};

//##############################################################################
//# Endl
//##############################################################################
/// \brief Class to add end of line (operator use).

class Endl{
public:
  Endl(){};
};

//##############################################################################
//# Head3  Eg: "cen [m];pos" -> "cen.x [m];cen.y [m];cen.z [m];pos.x;pos.y;pos.z"
//##############################################################################
/// \brief Class to add header for triple data (operator use).

class Head3
{
public:
  const std::string HeadText;
  Head3(std::string headtext):HeadText(headtext){};
};


//##############################################################################
//# JSaveCsv2
//##############################################################################
/// \brief Creates CSV files in a simple way.

class JSaveCsv2 : protected JObject
{
private:
  bool ExceptionThrown;    ///<An exception was thrown.
  std::string FileName;    ///<Name of file to store data.
  const bool App;          ///<Append mode enabled.
  const bool CsvSepComa;   ///<Separator character in CSV files (0=semicolon, 1=coma).

  std::fstream *Pf;
  bool AppendMode;         ///<File already exists and was opened to append data.
  bool FirstSaveData;

  bool AutoSepEnable; ///<Enable or disable auto separator change according configuration (def=true).

  static const unsigned SizeFmt=18;  ///<Number of different formats.
  std::string FmtDefault[SizeFmt];   ///<Default formats.
  std::string FmtCurrent[SizeFmt];   ///<Current formats.

  bool DataSelected;
  std::string Head;
  bool HeadLineEmpty;
  std::string Data;
  bool DataLineEmpty;

  void InitFmt();
  void AddStr(std::string tx);
  void AddSeparator(unsigned count);
  void AddEndl();
  void AddHead3(std::string headtx);

  void Save(const std::string &tx);
  void SetSeparators(std::string &tx)const;
  void OpenFile();
  //void R unException(const std::string &method,const std::string &msg){
  //  ExceptionThrown=true; JObject::R unException(method,msg);
  //}
  //void R unException(const std::string &method,const std::string &msg,const std::string &file){
  //  ExceptionThrown=true; JObject::R unException(method,msg,file);
  //}

public:
  JSaveCsv2(std::string fname,bool app,bool csvsepcoma);
  ~JSaveCsv2();
  void Reset();

  std::string GetFileName()const{ return(FileName); }
  bool GetAppendMode()const{ return(AppendMode); }

  const std::string& GetFmtDefault(TpFormat tfmt)const{ return(FmtDefault[tfmt]); };
  const std::string& GetFmtCurrent(TpFormat tfmt)const{ return(FmtCurrent[tfmt]); };
  void SetFmtCurrent(TpFormat tfmt,const std::string &fmt){ FmtCurrent[tfmt]=fmt; };

  void SetHead(){ DataSelected=false; }
  void SetData(){ DataSelected=true;  }

  std::string ToStr(const char *format,...)const;

  JSaveCsv2& operator <<(const Sep &obj);
  JSaveCsv2& operator <<(const AutoSepOn &obj);
  JSaveCsv2& operator <<(const AutoSepOff &obj);
  JSaveCsv2& operator <<(const Endl &obj);
  JSaveCsv2& operator <<(const Fmt &obj);
  JSaveCsv2& operator <<(const Head3 &obj);

  JSaveCsv2& operator <<(const std::string &v){ AddStr(v); return(*this); }

  JSaveCsv2& operator <<(char     v){ AddStr(ToStr(FmtCurrent[TpSigned1  ].c_str(),v)); return(*this); }
  JSaveCsv2& operator <<(short    v){ AddStr(ToStr(FmtCurrent[TpSigned1  ].c_str(),v)); return(*this); }
  JSaveCsv2& operator <<(int      v){ AddStr(ToStr(FmtCurrent[TpSigned1  ].c_str(),v)); return(*this); }
  JSaveCsv2& operator <<(byte     v){ AddStr(ToStr(FmtCurrent[TpUnsigned1].c_str(),v)); return(*this); }
  JSaveCsv2& operator <<(word     v){ AddStr(ToStr(FmtCurrent[TpUnsigned1].c_str(),v)); return(*this); }
  JSaveCsv2& operator <<(unsigned v){ AddStr(ToStr(FmtCurrent[TpUnsigned1].c_str(),v)); return(*this); }

  JSaveCsv2& operator <<(llong    v){ AddStr(ToStr(FmtCurrent[TpLlong1   ].c_str(),v)); return(*this); }
  JSaveCsv2& operator <<(ullong   v){ AddStr(ToStr(FmtCurrent[TpUllong1  ].c_str(),v)); return(*this); }

  JSaveCsv2& operator <<(float    v){ AddStr(ToStr(FmtCurrent[TpFloat1   ].c_str(),v)); return(*this); }
  JSaveCsv2& operator <<(double   v){ AddStr(ToStr(FmtCurrent[TpDouble1  ].c_str(),v)); return(*this); }

  JSaveCsv2& operator <<(const tint2   &v){ AddStr(ToStr(FmtCurrent[TpSigned2  ].c_str(),v.x,v.y));         return(*this); }
  JSaveCsv2& operator <<(const tint3   &v){ AddStr(ToStr(FmtCurrent[TpSigned3  ].c_str(),v.x,v.y,v.z));     return(*this); }
  JSaveCsv2& operator <<(const tint4   &v){ AddStr(ToStr(FmtCurrent[TpSigned4  ].c_str(),v.x,v.y,v.z,v.w)); return(*this); }

  JSaveCsv2& operator <<(const tuint2  &v){ AddStr(ToStr(FmtCurrent[TpUnsigned2].c_str(),v.x,v.y));         return(*this); }
  JSaveCsv2& operator <<(const tuint3  &v){ AddStr(ToStr(FmtCurrent[TpUnsigned3].c_str(),v.x,v.y,v.z));     return(*this); }
  JSaveCsv2& operator <<(const tuint4  &v){ AddStr(ToStr(FmtCurrent[TpUnsigned4].c_str(),v.x,v.y,v.z,v.w)); return(*this); }

  JSaveCsv2& operator <<(const tfloat2  &v){ AddStr(ToStr(FmtCurrent[TpFloat2  ].c_str(),v.x,v.y));         return(*this); }
  JSaveCsv2& operator <<(const tfloat3  &v){ AddStr(ToStr(FmtCurrent[TpFloat3  ].c_str(),v.x,v.y,v.z));     return(*this); }
  JSaveCsv2& operator <<(const tfloat4  &v){ AddStr(ToStr(FmtCurrent[TpFloat4  ].c_str(),v.x,v.y,v.z,v.w)); return(*this); }

  JSaveCsv2& operator <<(const tdouble2 &v){ AddStr(ToStr(FmtCurrent[TpDouble2 ].c_str(),v.x,v.y));         return(*this); }
  JSaveCsv2& operator <<(const tdouble3 &v){ AddStr(ToStr(FmtCurrent[TpDouble3 ].c_str(),v.x,v.y,v.z));     return(*this); }
  JSaveCsv2& operator <<(const tdouble4 &v){ AddStr(ToStr(FmtCurrent[TpDouble4 ].c_str(),v.x,v.y,v.z,v.w)); return(*this); }

  JSaveCsv2& operator <<(const tmatrix4d &v){ 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a11)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a12)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a13)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a14)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a21)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a22)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a23)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a24)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a31)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a32)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a33)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a34)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a41)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a42)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a43)); 
    AddStr(ToStr(FmtCurrent[TpDouble1].c_str(),v.a44)); 
    return(*this); 
  }

  void SaveData(bool closefile=false);
};

}

#endif


