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
//:# - Nuevas funciones para gestion de nombres de ficheros: GetFile(), 
//:#   GetFileNameSplit(). (10-08-2010)
//:# - Nuevas funciones para pasar de valores numericos a texto: UintStr(),
//:#   IntStr(). (17-12-2010)
//:# - Nueva funcion para ruats de ficheros: GetWithoutExtension(). (22-12-2010)
//:# - Funciones para convertir datos entre BigEndian y LittleEndian. (09-03-2011)
//:# - Agrupa funciones en namespace fun. (09-03-2011)
//:# - Nuevas funciones FileExists() y DirExists(). (10-03-2011)
//:# - Correccion en GetExtension() y GetWithoutExtension(), ahora busca la 
//:#   extension apartir del punto del ultimo fichero o directorio. (08-05-2011)
//:# - Nuevas funciones VarStr() para vectores de datos. (02-11-2011)
//:# - Funcion StrSplit() para extraer partes de un texto. (27-01-2012)
//:# - Traduccion de comentarios al ingles. (10-02-2012)
//:# - Error corregido en ReverseByteOrder(). (21-03-2012)
//:# - Nuevas funciones ResizeAlloc para redimensionar la cantidad de memoria
//:#   reservada conservando los datos. (22-03-2012)
//:# - Nuevas funciones GetDateTimeFormat() y GetTextRandomCode(). (29-09-2012)
//:# - Nueva funcion LongStr(). (05-04-2013)
//:# - Mejora de funcion StrSplit(). (06-06-2013)
//:# - Nuevas funciones StrSplitCount() y StrSplitValue(). (06-06-2013)
//:# - Algunas funciones nuevas para tipos double. (14-11-2013)
//:# - Nueva funcion StrWithoutChar(). (13-12-2013)
//:# - Nueva funcion ResizeAlloc para tfloat4 y tdouble2. (21-12-2013)
//:# - Nueva funcion PrintStr usando argumentos como el printf(). (10-03-2014)
//:# - Nuevos metodos VarStr para arrays de unsigned y word. (18-03-2014)
//:# - Nuevas funcion StrTrimRepeated(). (08-05-2014)
//:# - Nuevas funcion StrRepeat(). (03-10-2014)
//:# - Nuevas funciones GetFirstValueXXX(). (15-12-2014)
//:# - Remplaza long long por llong. (01-10-2015)
//:# - Nueva funcion StrOnlyChars(). (03-10-2015)
//:# - Nuevas funciones ResizeAlloc(). (10-10-2015)
//:# - Nuevas funciones IsInfinity() y IsNaN(). (18-10-2015)
//:# - Nuevas funciones StrTo valor numerico. (14-01-2016)
//:# - Nuevas funciones GetCanonicalPath() y GetPathLevels(). (07-01-2016)
//:# - Error coregido en PrintStr() cuando el string era demasiado grande. (20-02-2017)
//:# - Nueva funcion UintStrFill(). (03-08-2017)
//:# - Nuevas funciones: GetCurrentDir(), Mkdir(), MkdirPath(). (17-08-2017)
//:# - Nuevas funciones: StrAddSlashes(), StrStripSlashes(). (19-08-2017)
//:# - Nuevas funciones: StrReplace(), CpyFile(). (02-10-2017)
//:# - Error corregido en StrReplace(). Solo permitia 4 remplazos. (23-10-2017)
//:# - Nuevas funciones StrCsvSep() y PrintStrCsv(). (23-10-2017)
//:# - Error corregido en GetCanonicalPath(). No soportaba rutas absolutas de 
//:#   Windows. (24-10-2017)
//:# - Nuevas funciones: Float3Str(),Float3xStr(),Float3xRangeStr(),DoublexStr()
//:#   ,Double3xStr(),Double3xRangeStr(). (31-01-2018)
//:# - Error corregido en MkdirPath(). No soportaba rutas absolutas de Linux. (21-05-2018)
//:# - Nuevas funciones StrFileToVector(), StrVectorToFile() y StrFileError(). (11-07-2018)
//:# - Nuevas funciones StrTrimBegin() y StrTrimEnd(). (29-08-2018)
//:# - Nuevas funciones GetDateTimeFormatUTC(), GetWeekDay(), GetYearDay() y GetWeekNumber(). (24-10-2018)
//:# - Nuevas funciones Delay() y GetRuntime(). (22-11-2018)
//:# - Nuevas funcion VectorSplitStr().  (26-02-2019)
//:# - Nueva funcion StrIsNumber().  (21-03-2019)
//:# - Nueva funcion FileSize().  (23-03-2019)
//:# - Nuevas funciones IsEqual(), IsGtEqual(), IsLtEqual().  (08-04-2019)
//:# - Nueva funcion VectorSplitDouble().  (05-06-2019)
//:# - Nuevas funciones StrIsIntegerNumber() and StrIsRealNumber().  (13-06-2019)
//:# - Updates FileType() for files larger than 2GB.  (20-06-2019)
//:# - StrIsIntegerNumber() tambien valida enteros con parte decimal nula.  (12-08-2019)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:# - Output functions in JSON format.  (16-09-2019)
//:# - Nueva funcion VectorFind().  (25-11-2019)
//:# - Nueva funcion CompareVersions().  (08-12-2019)
//:# - Error corregido en StrStripSlashes().  (04-01-2020)
//:# - Nueva funcion VectorFindMask().  (17-03-2020)
//:# - Nueva funcion NaturalFmt().  (15-04-2020)
//:# - Nuevas funciones NewToTFloat3() y NewToTDouble3().  (09-05-2020)
//:# - Nueva funcion FileModifTime().  (02-06-2020)
//:# - Nueva funcion Length() para vectores.  (27-06-2020)
//:# - Nueva funcion VectorLower().  (09-09-2020)
//:# - Nueva funcion GetFirstTextBetween().  (20-09-2020)
//:#############################################################################

/// \file Functions.h \brief Declares basic/general functions for the entire application.

#ifndef _Functions_
#define _Functions_

#include <ctime>
#include <string>
#include <vector>
#include <sys/stat.h>
#include "TypesDef.h"
#include "RunExceptionDef.h"


/// Implements a set of basic/general functions.
namespace fun{
void RunExceptioonFun(const std::string &srcfile,int srcline,const std::string &fun
  ,const std::string &msg,const std::string &file="");

std::string GetDateTimeFormat(const char* format,int nseg=0);
inline std::string GetDateTime(){ return(GetDateTimeFormat("%d-%m-%Y %H:%M:%S",0)); }
inline std::string GetDateTimeAfter(int nseg){ return(GetDateTimeFormat("%d-%m-%Y %H:%M:%S",nseg)); }

std::string GetDateTimeFormatUTC(const char* format,int day,int month,int year,int hour=0,int min=0,int sec=0);
int GetWeekDay(int day,int month,int year);
int GetYearDay(int day,int month,int year);
int GetWeekNumber(int day,int month,int year);

void Delay(int ms);
double GetRuntime();

std::string GetHoursOfSeconds(double s);

std::string GetTextRandomCode(unsigned length);

std::string PrintStr(const char *format,...);

std::string PrintStrCsv(bool csvsepcoma,const char *format,...);
std::string StrCsvSep(bool csvsepcoma,const std::string &cad);

std::string NaturalFmt(double v,unsigned ndigits,bool removezeros);

std::string IntStrFill(int v,int vmax);
std::string UintStrFill(unsigned v,unsigned vmax,const char fillchar='0');
std::string LongStr(llong v);
std::string UlongStr(ullong v);
std::string UintStr(unsigned v,const char* fmt="%u");
std::string IntStr(int v);
std::string Int3Str(const tint3 &v);
std::string Uint3Str(const tuint3 &v);
/// Converts range of tint3 values to string.  
inline std::string Int3RangeStr(const tint3 &v,const tint3 &v2){ return(std::string("(")+Int3Str(v)+")-("+Int3Str(v2)+")"); }
/// Converts range of tuint3 values to string.  
inline std::string Uint3RangeStr(const tuint3 &v,const tuint3 &v2){ return(std::string("(")+Uint3Str(v)+")-("+Uint3Str(v2)+")"); }

std::string FloatStr(float v,const char* fmt="%f");
std::string FloatxStr(float v,const char* fmt="%f");
std::string Float3Str(const tfloat3 &v,const char* fmt="%f,%f,%f");
/// Converts real value to string with format g.
inline std::string Float3gStr(const tfloat3 &v){ return(Float3Str(v,"%g,%g,%g")); }
/// Converts real value to string (-FLT_MAX=MIN and FLT_MAX=MAX).
inline std::string Float3xStr(const tfloat3 &v,const char* fmt="%f"){ return(FloatxStr(v.x,fmt)+","+FloatxStr(v.y,fmt)+","+FloatxStr(v.z,fmt)); }
/// Converts range of tfloat3 values to string.  
inline std::string Float3gRangeStr(const tfloat3 &v,const tfloat3 &v2){ return(std::string("(")+Float3gStr(v)+")-("+Float3gStr(v2)+")"); }
/// Converts range of tfloat3 values to string (-FLT_MAX=MIN and FLT_MAX=MAX).
inline std::string Float3xRangeStr(const tfloat3 &v,const tfloat3 &v2,const char* fmt="%f"){ return(std::string("(")+Float3xStr(v,fmt)+")-("+Float3xStr(v2,fmt)+")"); }

std::string DoubleStr(double v,const char* fmt="%g");
std::string DoublexStr(double v,const char* fmt="%f");
std::string Double3Str(const tdouble3 &v,const char* fmt="%f,%f,%f");
/// Converts real values to string with format g.
inline std::string Double3gStr(const tdouble3 &v){ return(Double3Str(v,"%g,%g,%g")); }
/// Converts real values to string (-DBL_MAX=MIN and DBL_MAX=MAX).
inline std::string Double3xStr(const tdouble3 &v,const char* fmt="%f"){ return(DoublexStr(v.x,fmt)+","+DoublexStr(v.y,fmt)+","+DoublexStr(v.z,fmt)); }
/// Converts range of tdouble3 values to string.  
inline std::string Double3gRangeStr(const tdouble3 &v,const tdouble3 &v2){ return(std::string("(")+Double3gStr(v)+")-("+Double3gStr(v2)+")"); }
/// Converts range of tdouble3 values to string (-DBL_MAX=MIN and DBL_MAX=MAX).
inline std::string Double3xRangeStr(const tdouble3 &v,const tdouble3 &v2,const char* fmt="%f"){ return(std::string("(")+Double3xStr(v,fmt)+")-("+Double3xStr(v2,fmt)+")"); }

std::string Double4Str(const tdouble4 &v,const char* fmt="%f,%f,%f");
/// Converts range of tdouble4 values to string.  
inline std::string Double4gStr(const tdouble4 &v){ return(Double4Str(v,"%g,%g,%g,%g")); }

std::string VectorStr(const std::vector<std::string> &v);


bool StrIsIntegerNumber(const std::string &v);
bool StrIsRealNumber(const std::string &v);

int      StrToInt    (const std::string &v);
tint3    StrToInt3   (std::string v);

double   StrToDouble (const std::string &v);
tdouble3 StrToDouble3(std::string v);

inline byte     StrToByte   (const std::string &v){ return(byte(StrToInt(v)));          }
inline word     StrToWord   (const std::string &v){ return(word(StrToInt(v)));          }
inline unsigned StrToUint   (const std::string &v){ return(unsigned(StrToInt(v)));      }
inline tuint3   StrToUint3  (const std::string &v){ return(ToTUint3(StrToInt3(v)));     }
inline float    StrToFloat  (const std::string &v){ return(float(StrToDouble(v)));      }
inline tfloat3  StrToFloat3 (const std::string &v){ return(ToTFloat3(StrToDouble3(v))); }

std::string StrUpper(const std::string &cad);
std::string StrLower(const std::string &cad);
std::string StrTrim(const std::string &cad);
std::string StrTrimBegin(const std::string &cad);
std::string StrTrimEnd(const std::string &cad);
std::string StrTrimRepeated(const std::string &cad);
std::string StrWithoutChar(const std::string &cad,char let);
std::string StrRepeat(const std::string &cad,unsigned count);
std::string StrReplace(const std::string &cad,const std::string &key,const std::string &newcad);
std::string StrAddSlashes(const std::string &cad);
std::string StrStripSlashes(const std::string &cad);

bool StrOnlyChars(const std::string &cad,const std::string &chars);

int StrFileToVector(const std::string &file,std::vector<std::string> &lines);
int StrVectorToFile(const std::string &file,const std::vector<std::string> &lines);
std::string StrFileError(int error);

std::string StrSplit(const std::string mark,std::string &text);
unsigned StrSplitCount(const std::string mark,std::string text);
std::string StrSplitValue(const std::string mark,std::string text,unsigned value);
unsigned VectorSplitStr(const std::string mark,const std::string &text,std::vector<std::string> &vec);
unsigned VectorSplitInt(const std::string mark,const std::string &text,std::vector<int> &vec);
unsigned VectorSplitDouble(const std::string mark,const std::string &text,std::vector<double> &vec);
unsigned VectorSplitFloat(const std::string mark,const std::string &text,std::vector<float> &vec);
void     VectorLower(std::vector<std::string> &vec);
unsigned VectorFind(const std::string &key,const std::vector<std::string> &vec,unsigned first=0);
unsigned VectorFindMask(const std::string &keymask,const std::vector<std::string> &vec,unsigned first=0);

double GetFirstValueDouble(std::string tex,std::string pretex="");
double GetFirstValueDouble(std::string tex,std::string &resttex,std::string pretex);
int GetFirstValueInt(std::string tex,std::string pretex="");
int GetFirstValueInt(std::string tex,std::string &resttex,std::string pretex);
std::string GetFirstTextBetween(std::string tex,std::string &resttex,std::string pretex,std::string endtex);

int CompareVersions(std::string v1,std::string v2);

std::string VarStr(const std::string &name,const char *value);
std::string VarStr(const std::string &name,const std::string &value);
std::string VarStr(const std::string &name,float value);
std::string VarStr(const std::string &name,tfloat3 value);
std::string VarStr(const std::string &name,double value);
std::string VarStr(const std::string &name,tdouble3 value);
std::string VarStr(const std::string &name,bool value);
std::string VarStr(const std::string &name,int value);
std::string VarStr(const std::string &name,unsigned value);

std::string VarStr(const std::string &name,unsigned n,const int *values,std::string size="?");
std::string VarStr(const std::string &name,unsigned n,const unsigned *values,std::string size="?");
std::string VarStr(const std::string &name,unsigned n,const word *values,std::string size="?");
std::string VarStr(const std::string &name,unsigned n,const float *values,std::string size="?",const char *fmt="%f");
std::string VarStr(const std::string &name,unsigned n,const double *values,std::string size="?",const char *fmt="%f");
std::string VarStr(const std::string &name,unsigned n,const tdouble3 *values,std::string size="?",const char *fmt="%g");

std::string VarStr(const std::string &name,const std::vector<int> &values,std::string size="?");
std::string VarStr(const std::string &name,const std::vector<tdouble3> &values,std::string size="?",const char *fmt="%g");



void PrintVar(const std::string &name,const char *value,const std::string &post="");
void PrintVar(const std::string &name,const std::string &value,const std::string &post="");
void PrintVar(const std::string &name,float value,const std::string &post="");
void PrintVar(const std::string &name,double value,const std::string &post="");
void PrintVar(const std::string &name,tfloat3 value,const std::string &post="");
void PrintVar(const std::string &name,tdouble3 value,const std::string &post="");
void PrintVar(const std::string &name,bool value,const std::string &post="");
void PrintVar(const std::string &name,int value,const std::string &post="");
void PrintVar(const std::string &name,unsigned value,const std::string &post="");

//-Output functions in JSON format.
inline std::string JSONValue(const std::string &v){ return(std::string(" \"")+v+"\" "); }
inline std::string JSONValue(bool     v){ return(v? " true ": " false "); }
inline std::string JSONValue(int      v){ return(std::string(" ")+IntStr(v)+" "); }
inline std::string JSONValue(unsigned v){ return(std::string(" ")+UintStr(v)+" "); }
inline std::string JSONValue(double   v){ return(std::string(" ")+DoubleStr(v)+" "); }
inline std::string JSONPropertyValue(const std::string &name,const std::string &v){ return(std::string(" \"")+name+"\" : "+StrTrimBegin(v)); }
inline std::string JSONProperty(const std::string &name,const std::string &v){ return(JSONPropertyValue(name,JSONValue(v))); }
inline std::string JSONProperty(const std::string &name,bool               v){ return(JSONPropertyValue(name,JSONValue(v))); }
inline std::string JSONProperty(const std::string &name,int                v){ return(JSONPropertyValue(name,JSONValue(v))); }
inline std::string JSONProperty(const std::string &name,unsigned           v){ return(JSONPropertyValue(name,JSONValue(v))); }
inline std::string JSONProperty(const std::string &name,double             v){ return(JSONPropertyValue(name,JSONValue(v))); }
std::string JSONObject(const std::vector<std::string> &properties);
std::string JSONArray(const std::vector<std::string> &values);

int FileType(const std::string &name);
ullong FileModifTime(const std::string &name);
inline bool FileExists(const std::string &name){ return(FileType(name)==2); }
inline bool DirExists(const std::string &name){ return(FileType(name)==1); }
llong FileSize(const std::string &name);

std::string GetCurrentDir();
int Mkdir(const std::string &dirname);
int MkdirPath(std::string path);

std::string GetDirParent(const std::string &ruta);
std::string GetCanonicalPath(std::string pathbase,std::string path);
std::string GetPathLevels(std::string path,unsigned levels);
std::string GetFile(const std::string &ruta);
std::string GetDirWithSlash(const std::string &ruta);
std::string GetDirWithoutSlash(const std::string &ruta);
std::string GetExtension(const std::string &file);
std::string GetWithoutExtension(const std::string &ruta);
void GetFileNameSplit(const std::string &file,std::string &dir,std::string &fname,std::string &fext);
std::string AddExtension(const std::string &file,const std::string &ext);
std::string FileNameSec(std::string fname,unsigned fnumber);
std::string ShortFileName(const std::string &file,unsigned maxlen,bool withpoints=true);
std::string TextWithShortFileName(const std::string &txpre,const std::string &txpos,const std::string &file,unsigned maxlen);

bool FileMask(std::string text,std::string mask);
int CpyFile(std::string src,const std::string dest);

typedef enum{ BigEndian=1,LittleEndian=0 }TpByteOrder;
TpByteOrder GetByteOrder();
void ReverseByteOrder(llong *data,int count,llong *result);
void ReverseByteOrder(int *data,int count,int *result);
void ReverseByteOrder(short *data,int count,short *result);
inline void ReverseByteOrder(llong *data,int count){ ReverseByteOrder(data,count,data); }
inline void ReverseByteOrder(int *data,int count){ ReverseByteOrder(data,count,data); }
inline void ReverseByteOrder(short *data,int count){ ReverseByteOrder(data,count,data); }

byte*     ResizeAlloc(byte     *data,unsigned ndata,unsigned newsize);
word*     ResizeAlloc(word     *data,unsigned ndata,unsigned newsize);
unsigned* ResizeAlloc(unsigned *data,unsigned ndata,unsigned newsize);
tuint2*   ResizeAlloc(tuint2   *data,unsigned ndata,unsigned newsize);
tuint3*   ResizeAlloc(tuint3   *data,unsigned ndata,unsigned newsize);
tuint4*   ResizeAlloc(tuint4   *data,unsigned ndata,unsigned newsize);
int*      ResizeAlloc(int      *data,unsigned ndata,unsigned newsize);
tint2*    ResizeAlloc(tint2    *data,unsigned ndata,unsigned newsize);
tint3*    ResizeAlloc(tint3    *data,unsigned ndata,unsigned newsize);
float*    ResizeAlloc(float    *data,unsigned ndata,unsigned newsize);
tfloat2*  ResizeAlloc(tfloat2  *data,unsigned ndata,unsigned newsize);
tfloat3*  ResizeAlloc(tfloat3  *data,unsigned ndata,unsigned newsize);
tfloat4*  ResizeAlloc(tfloat4  *data,unsigned ndata,unsigned newsize);
double*   ResizeAlloc(double   *data,unsigned ndata,unsigned newsize);
tdouble2* ResizeAlloc(tdouble2 *data,unsigned ndata,unsigned newsize);
tdouble3* ResizeAlloc(tdouble3 *data,unsigned ndata,unsigned newsize);
tdouble4* ResizeAlloc(tdouble4 *data,unsigned ndata,unsigned newsize);

tfloat3*  NewToTFloat3 (const tdouble3* data,unsigned ndata);
tdouble3* NewToTDouble3(const tfloat3* data,unsigned ndata);

float  Length(const tfloat3 &v);
double Length(const tdouble3 &v);

bool IsInfinity(float v);
bool IsInfinity(double v);
bool IsNAN(float v);
bool IsNAN(double v);

bool IsEqual(float  v1,float  v2,float  tolerance);
bool IsEqual(double v1,double v2,double tolerance);
bool IsGtEqual(float  v1,float  v2,float  tolerance);
bool IsGtEqual(double v1,double v2,double tolerance);
bool IsLtEqual(float  v1,float  v2,float  tolerance);
bool IsLtEqual(double v1,double v2,double tolerance);

bool IsEqual(const tdouble3 &v1,const tdouble3 &v2,double tolerance);
bool IsEqual(const tdouble4 &v1,const tdouble4 &v2,double tolerance);

}

#endif


