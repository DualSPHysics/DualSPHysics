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
//:# Clase base para la definicion de cualquier formato binario de fichero.
//:# Algunas de sus funcionalidades son:
//:# - Permite guardar datos de todo tipo en formato binario.
//:# - Todos los datos contienen nombre y tipo para permitir su exploracion sin
//:#   sin conocer su estructura.
//:# - Los datos se pueden organizar en distintos niveles con forma de arbol.
//:# - Permite guardar datos basicos y arrays de tipos basicos y strings.
//:# - La lectura del contenido de arrays desde fichero puede ser de forma 
//:#   selectiva, solo los que se necesiten.
//:# - Los arrays se pueden redimensionar de forma automatica segun se vayan
//:#   introduciendo mas datos. 
//:# - Permite uso de punteros externos para reducir el consumo de memoria. 
//:# - Implementacion de constructor de copia y sobrecarga de asigancion. 
//:# - Permite la grabacion y lectuara de varios items consecutivos en un fichero. 
//:#   Esto permite la ampliacion de ficheros sin reescribir todo el contenido.
//:#
//:# Cambios:
//:# =========
//:# - Implementacion. (25-08-2013 <-> 23-11-2013)
//:# - Ahora el metodo SaveFileListApp() graba los datos del Parent al principio
//:#   del fichero. (12-01-2014)
//:# - Opcion en SaveFileXml() para grabar datos de arrays. (04-12-2014)
//:# - Nuevos metodos CheckCopyArrayData() y CopyArrayData(). (13-04-2020)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Lanza excepcion cuando el fichero no es soportado por ser mayor de 4GB. (25-10-2021)
//:# - Nuevo metodo CreateArrayFloat3() para facilitar uso eficiente. (31-10-2021)
//:# - Uso de size_t en lugar de unsigned para gestionar mas valores y mayor 
//:#   volumen de datos. (20-06-2023)
//:# - Soporta volumenes de datos superiores a 4GB (head.si64=true). (20-06-2023) 
//:# - Solo usa head.si64=true para grabar volumenes grandes de forma automatica y
//:#   guarda head.byteorder+10 para provocar una excepcion en codigos de lectura 
//:#   obsoletos que no reconocen head.si64. (20-06-2023)
//:# - Graba datos en partes de MaxSizeSi32 para evitar errores de escritura y 
//:#   comprueba lecturas y escrituras principales. (21-06-2023) 
//:# - Nuevos metodos CreateArrayTYPE() para facilitar un uso eficiente. (23-07-2023)
//:# - Nuevos metodos GetArrayType() y GetArrayTpSize() para obtener arrays comprobados. (28-07-2023)
//:#############################################################################

/// \file JBinaryData.h \brief Declares the class \ref JBinaryData.

#ifndef _JBinaryData_
#define _JBinaryData_

#include "JObject.h"
#include "TypesDef.h"
#include <string>
#include <vector>
#include <fstream>
#include <climits>
#include <cstdint>

class JBinaryData;

//##############################################################################
//# JBinaryDataDef
//##############################################################################
/// \brief Defines types to be used in \ref JBinaryData and \ref JBinaryDataArray
// Define tipos usados en \ref JBinaryData y \ref JBinaryDataArray

class JBinaryDataDef
{
 public:
  typedef enum{ DatNull=0,DatText=1,DatBool=2,DatChar=3,DatUchar=4
    ,DatShort=5,DatUshort=6,DatInt=7,DatUint=8,DatLlong=9,DatUllong=10
    ,DatFloat=11,DatDouble=12
    ,DatInt3=20,DatUint3=21,DatFloat3=22,DatDouble3=23 
  }TpData; 

  static void RunExceptioonStatic(const std::string& srcfile,int srcline
    ,const std::string& method
    ,const std::string& msg,const std::string& file="");

  static std::string TypeToStr(TpData type);
  static size_t SizeOfType(TpData type);
  static bool TypeIsTriple(TpData type);

  static const bool ForceSi64=false; //-False by default.
};


//##############################################################################
//# JBinaryDataArray
//##############################################################################
/// \brief Defines data arrays included in binary files.
// Define arrays de datos incluidos en ficheros binarios.

class JBinaryDataArray : protected JObject
{
  //friend class JBinaryData;
 private:
  JBinaryData* Parent;
  std::string Name;
  bool Hide;
  const JBinaryDataDef::TpData Type;
  size_t Count;         ///<Numero de elementos almacenados en pointer. Number of elements stored in pointer.
  size_t Size;          ///<Numero de elementos para los que hay memoria reservada. Number of elements for which there is reserved memory.
  void* Pointer;
  bool ExternalPointer; ///<Indica que el puntero es externo y no debe liberarse. Indicates that the pointer is external, and should not be released.
  size_t FileDataPox;   ///<Valor mayor o igual a cero indica la posicion de lectura en el fichero abierto en el ItemHead. Value greater than or equal to zero indicates the position of reading in the file opened in the ItemHead.
  size_t FileDataCount; ///<Numero de elemetos del array en fichero. Number of elements in the array in a file.
  size_t FileDataSize;  ///<Size de datos del array en fichero. Size of array data in file.

  void FreePointer(void* ptr)const;
  void* AllocPointer(size_t size)const;
  void CheckMemory(size_t count,bool resize);

  void        OutData(size_t& count,size_t size,const byte* ptr,byte* dat,size_t sdat)const;
  std::string OutStr (size_t& count,size_t size,const byte* ptr)const;
  /// Extrae unsigned de ptr. Extract ptr unsigned.
  unsigned    OutUint(size_t& count,size_t size,const byte* ptr)const{  
    unsigned v; OutData(count,size,ptr,(byte*)&v,sizeof(unsigned)); return(v);
  }  

 public:
  JBinaryDataArray(JBinaryData* parent,const std::string& name,JBinaryDataDef::TpData type);
  ~JBinaryDataArray();
  size_t GetAllocMemory()const;

  void SetName(const std::string& name);
  std::string GetName()const{ return(Name); };

  void SetHide(bool hide){ Hide=hide; }
  bool GetHide()const{ return(Hide); }

  JBinaryDataDef::TpData GetType()const{ return(Type); };
  size_t GetCount()const{ return(Count); };
  size_t GetSize()const{ return(Size); };
  const void* GetPointer()const{ return(Pointer); };

  bool PointerIsExternal()const{ return(ExternalPointer); };
  bool DataInPointer()const{ return(Pointer && Count); }
  bool DataInFile()const{ return(FileDataPox!=SIZE_MAX); }

  void FreeMemory();
  void AllocMemory(size_t size,bool savedata=false);
  void AllocMemoryCount(size_t count,bool clear);
  void ConfigExternalMemory(size_t size,void* pointer);

  void ReadData(size_t count,size_t size,std::ifstream* pf,bool resize);
  void AddData(size_t count,const void* data,bool resize);
  void SetData(size_t count,const void* data,bool externalpointer);

  const void* GetDataPointer()const;
  size_t GetDataCopy(size_t size,void* pointer)const;

  void AddText(const std::string& str,bool resize);
  void AddTexts(size_t count,const std::string* strs,bool resize);

  void ConfigFileData(size_t filepos,size_t datacount,size_t datasize);
  void ClearFileData();
  size_t GetFileDataCount()const{ return(FileDataCount); }
  size_t GetFileDataSize()const{ return(FileDataSize); }
  void ReadFileData(bool resize);
};

//##############################################################################
//# JBinaryData
//##############################################################################
/// \brief Defines any binary format of a file.
// Clase base para la definicion de cualquier formato binario de fichero.

class JBinaryData : protected JObject
{
 private:
  ///Structure that describes the header of binary format files.
  typedef struct{
    char titu[60];               ///<Title of the file eg: "#File JBinaryData".
    byte byteorder;              ///<1:BigEndian, 0:LittleEndian (+10 when si64==true).
    byte si64;                   ///<1:64-bit size for large elements, 0:32-bit size
    byte void2;                  ///<Not used.
    byte void3;                  ///<Not used.
  }StHeadFmtBin;//-sizeof(64)

  static const std::string CodeItemDef;
  static const std::string CodeValuesDef;
  static const std::string CodeArrayDef;
  static const unsigned MaxSizeSi32=4228102018u;

 public:

  ///Structure that describes the information of a value.
  typedef struct{
    std::string name;
    JBinaryDataDef::TpData type;
    std::string vtext;
    union{
      char vchar;
      byte vuchar;
      short vshort;
      unsigned vushort;
      int vint;
      unsigned vuint;
      llong vllong;
      ullong vullong;
      float vfloat;
      double vdouble;
      tint3 vint3;
      tuint3 vuint3;
      tfloat3 vfloat3;
      tdouble3 vdouble3;   //- Elemento de mayor tamanho utilizado para poner a Zero. Large item used to zero elements.
    };
  }StValue;

 private:
  std::string Name;      ///<Nombre de item. Name of item.
  bool HideAll;          ///<Ignora el item en determinados metodos como SaveData(). It ignores the item in certain functions as SaveData().
  bool HideValues;       ///<Ignora los Values en determinados metodos como SaveData(). It ignores the values in certain functions as SaveData().
  std::string FmtFloat;  ///<Formato para valores float, por defecto "%.7E". Format for float, by default " %.7E". 
  std::string FmtDouble; ///<Formato para valores double, por defecto "%.15E" Format for double, by default " %.15E".

  JBinaryData* Parent;
  std::vector<JBinaryData*> Items;
  std::vector<JBinaryDataArray*> Arrays;
  std::vector<StValue> Values;

  std::ifstream* FileStructure;

  //-Variables para cache de values. Variables to cache values.
  bool ValuesModif;
  byte* ValuesData;
  size_t ValuesSize;

  //-Gestion de Values. Management of Values.
  void ValuesCacheReset();
  size_t ChecksGetValue(const std::string& name,bool optional,JBinaryDataDef::TpData type)const;
  size_t ChecksSetValue(const std::string& name,JBinaryDataDef::TpData type);
  static void ResetValue(const std::string& name,JBinaryDataDef::TpData type,StValue& v);
  std::string ValueToXml(const StValue& v)const;

  void InData   (size_t& count,size_t size,byte* ptr,const byte* dat,size_t sdat)const;
  void InStr    (size_t& count,size_t size,byte* ptr,const std::string& cad)const;
  void InBool   (size_t& count,size_t size,byte* ptr,bool v)const{     int vv=(v? 1: 0); InInt(count,size,ptr,vv);                }  ///<Introduce bool en ptr. Introduces bool in ptr.
  void InChar   (size_t& count,size_t size,byte* ptr,char v)const{     InData(count,size,ptr,(byte*)&v,sizeof(char));     }  ///<Introduce char en ptr. Introduces char in ptr.
  void InUchar  (size_t& count,size_t size,byte* ptr,byte v)const{     InData(count,size,ptr,(byte*)&v,sizeof(byte));     }  ///<Introduce byte en ptr. Introduces byte in ptr.
  void InShort  (size_t& count,size_t size,byte* ptr,short v)const{    InData(count,size,ptr,(byte*)&v,sizeof(short));    }  ///<Introduce short en ptr. Introduces short in ptr.
  void InUshort (size_t& count,size_t size,byte* ptr,word v)const{     InData(count,size,ptr,(byte*)&v,sizeof(word));     }  ///<Introduce word en ptr. Introduces word in ptr.
  void InInt    (size_t& count,size_t size,byte* ptr,int v)const{      InData(count,size,ptr,(byte*)&v,sizeof(int));      }  ///<Introduce int en ptr. Introduces int in ptr.
  void InUint   (size_t& count,size_t size,byte* ptr,unsigned v)const{ InData(count,size,ptr,(byte*)&v,sizeof(unsigned)); }  ///<Introduce unsigned en ptr. Introduces unsigned in ptr
  void InLlong  (size_t& count,size_t size,byte* ptr,llong v)const{    InData(count,size,ptr,(byte*)&v,sizeof(llong));    }  ///<Introduce long long en ptr. Introduces long long in ptr
  void InUllong (size_t& count,size_t size,byte* ptr,ullong v)const{   InData(count,size,ptr,(byte*)&v,sizeof(ullong));   }  ///<Introduce unsigned long long en ptr. Introduces unsigned long in ptr
  void InFloat  (size_t& count,size_t size,byte* ptr,float v)const{    InData(count,size,ptr,(byte*)&v,sizeof(float));    }  ///<Introduce float en ptr. Introduces float in ptr
  void InDouble (size_t& count,size_t size,byte* ptr,double v)const{   InData(count,size,ptr,(byte*)&v,sizeof(double));   }  ///<Introduce double en ptr. Introduces double in ptr
  void InInt3   (size_t& count,size_t size,byte* ptr,tint3 v)const{    InData(count,size,ptr,(byte*)&v,sizeof(tint3));    }  ///<Introduce tint3 en ptr. Introduces tint3 in ptr
  void InUint3  (size_t& count,size_t size,byte* ptr,tuint3 v)const{   InData(count,size,ptr,(byte*)&v,sizeof(tuint3));   }  ///<Introduce tuint3 en ptr. Introduces tuint3 in ptr
  void InFloat3 (size_t& count,size_t size,byte* ptr,tfloat3 v)const{  InData(count,size,ptr,(byte*)&v,sizeof(tfloat3));  }  ///<Introduce tfloat3 en ptr. Introduces tfloat3 in ptr
  void InDouble3(size_t& count,size_t size,byte* ptr,tdouble3 v)const{ InData(count,size,ptr,(byte*)&v,sizeof(tdouble3)); }  ///<Introduce tdouble3 en ptr. Introduces tdouble3 in ptr

  void         OutData   (size_t& count,size_t size,const byte* ptr,byte* dat,size_t sdat)const;
  std::string  OutStr    (size_t& count,size_t size,const byte* ptr)const;
  bool         OutBool   (size_t& count,size_t size,const byte* ptr)const{  return(OutInt(count,size,ptr)!=0);  }  /// Extrae bool de ptr.
  char         OutChar   (size_t& count,size_t size,const byte* ptr)const{  char v;     OutData(count,size,ptr,(byte*)&v,sizeof(char));     return(v);  }  ///< Extrae char de ptr. Extracts char of ptr.
  byte         OutUchar  (size_t& count,size_t size,const byte* ptr)const{  byte v;     OutData(count,size,ptr,(byte*)&v,sizeof(byte));     return(v);  }  ///< Extrae byte de ptr. Extracts byte of ptr.
  short        OutShort  (size_t& count,size_t size,const byte* ptr)const{  short v;    OutData(count,size,ptr,(byte*)&v,sizeof(short));    return(v);  }  ///< Extrae short de ptr. Extracts short of ptr.
  word         OutUshort (size_t& count,size_t size,const byte* ptr)const{  word v;     OutData(count,size,ptr,(byte*)&v,sizeof(word));     return(v);  }  ///< Extrae word de ptr. Extracts word of ptr.
  int          OutInt    (size_t& count,size_t size,const byte* ptr)const{  int v;      OutData(count,size,ptr,(byte*)&v,sizeof(int));      return(v);  }  ///< Extrae int de ptr. Extracts int of ptr.
  unsigned     OutUint   (size_t& count,size_t size,const byte* ptr)const{  unsigned v; OutData(count,size,ptr,(byte*)&v,sizeof(unsigned)); return(v);  }  ///< Extrae unsigned de ptr. Extracts unsigned of ptr.
  llong        OutLlong  (size_t& count,size_t size,const byte* ptr)const{  llong v;    OutData(count,size,ptr,(byte*)&v,sizeof(llong));    return(v);  }  ///< Extrae long long de ptr. Extracts long long of ptr.
  ullong       OutUllong (size_t& count,size_t size,const byte* ptr)const{  ullong v;   OutData(count,size,ptr,(byte*)&v,sizeof(ullong));   return(v);  }  ///< Extrae unsigned long long de ptr. Extracts unsigned long long of ptr.
  float        OutFloat  (size_t& count,size_t size,const byte* ptr)const{  float v;    OutData(count,size,ptr,(byte*)&v,sizeof(float));    return(v);  }  ///< Extrae float de ptr. Extracts float of ptr.
  double       OutDouble (size_t& count,size_t size,const byte* ptr)const{  double v;   OutData(count,size,ptr,(byte*)&v,sizeof(double));   return(v);  }  ///< Extrae double de ptr. Extracts double of ptr.
  tint3        OutInt3   (size_t& count,size_t size,const byte* ptr)const{  tint3 v;    OutData(count,size,ptr,(byte*)&v,sizeof(tint3));    return(v);  }  ///< Extrae tint3 de ptr. Extracts tint3 of ptr.
  tuint3       OutUint3  (size_t& count,size_t size,const byte* ptr)const{  tuint3 v;   OutData(count,size,ptr,(byte*)&v,sizeof(tuint3));   return(v);  }  ///< Extrae tuint3 de ptr. Extracts tuint3 of ptr.
  tfloat3      OutFloat3 (size_t& count,size_t size,const byte* ptr)const{  tfloat3 v;  OutData(count,size,ptr,(byte*)&v,sizeof(tfloat3));  return(v);  }  ///< Extrae tfloat3 de ptr. Extracts tfloat3 of ptr.
  tdouble3     OutDouble3(size_t& count,size_t size,const byte* ptr)const{  tdouble3 v; OutData(count,size,ptr,(byte*)&v,sizeof(tdouble3)); return(v);  }  ///< Extrae tdouble3 de ptr. Extracts tdouble3 of ptr.

  void InValue(size_t& count,size_t size,byte* ptr,const StValue& v)const;
  void OutValue(size_t& count,size_t size,const byte* ptr);

  void InArrayBase(bool si64,size_t& count,size_t size,byte* ptr,const JBinaryDataArray* ar)const;
  void InArrayData(size_t& count,size_t size,byte* ptr,const JBinaryDataArray* ar)const;
  void InArray(bool si64,size_t& count,size_t size,byte* ptr,const JBinaryDataArray* ar)const;
  void InItemBase(size_t& count,size_t size,byte* ptr,bool all)const;
  void InItem(bool si64,size_t& count,size_t size,byte* ptr,bool all)const;

  JBinaryDataArray* OutArrayBase(bool si64,size_t& count,size_t size,const byte* ptr,size_t& countdata,size_t& sizedata);
  void OutArrayData(size_t& count,size_t size,const byte* ptr,JBinaryDataArray* ar,size_t countdata,size_t sizedata);
  void OutArray(bool si64,size_t& count,size_t size,const byte* ptr);
  JBinaryData* OutItemBase(size_t& count,size_t size,const byte* ptr,bool create,size_t& narrays,size_t& nitems,size_t& sizevalues);
  void OutItem(bool si64,size_t& count,size_t size,const byte* ptr,bool create);

  size_t GetSizeValues()const;
  void SaveValues(size_t& count,size_t size,byte* ptr)const;
  void ValuesCachePrepare(bool down);

  void WriteArrayData(std::fstream* pf,const JBinaryDataArray* ar)const;
  void WriteArray(bool si64,std::fstream* pf,size_t sbuf,byte* buf,const JBinaryDataArray* ar)const;
  void WriteItem(bool si64,std::fstream* pf,size_t sbuf,byte* buf,bool all)const;

  unsigned ReadUint(std::ifstream* pf)const;
  void ReadArrayData(std::ifstream* pf,JBinaryDataArray* ar,size_t countdata,size_t sizedata,bool loadarraysdata);
  void ReadArray(bool si64,std::ifstream* pf,size_t sbuf,byte* buf,bool loadarraysdata);
  void ReadItem(bool si64,std::ifstream* pf,size_t sbuf,byte* buf,bool create,bool loadarraysdata);

  JBinaryData::StHeadFmtBin MakeFileHead(bool si64,const std::string& filecode)const;
  size_t GetFileHead(std::ifstream* pf,JBinaryData::StHeadFmtBin& head)const;
  void CheckHead(const std::string& file,const StHeadFmtBin& head,const std::string& filecode)const;
  size_t CheckFileHead(const std::string& file,std::ifstream* pf,const std::string& filecode,bool& file_si64)const;
  size_t CheckFileListHead(const std::string& file,std::fstream* pf,const std::string& filecode,bool& file_si64)const;
  bool SaveFileData(std::fstream* pf,bool head,const std::string& filecode,bool memory,bool all)const;

  void WriteFileXmlArray(const std::string& tabs,std::ofstream* pf,bool svarrays,const JBinaryDataArray* ar)const;

  void WriteFileXml(const std::string& tabs,std::ofstream* pf,bool svarrays)const;

 public:
  JBinaryData(std::string name="JBinary_Data");
  JBinaryData(const JBinaryData& src);
  ~JBinaryData();
  JBinaryData& operator=(const JBinaryData& src);
  void Clear();
  size_t GetAllocMemory()const;

  void SetName(const std::string& name);
  std::string GetName()const{ return(Name); };
  
  void SetHide(bool hide){ HideAll=hide; }
  bool GetHide()const{ return(HideAll); }
  void SetHideValues(bool hide,bool down);
  bool GetHideValues()const{ return(HideValues); }
  void SetHideArrays(bool hide,bool down);
  void SetHideItems(bool hide,bool down);

  void SetFmtFloat(const std::string& fmt,bool down);
  void SetFmtDouble(const std::string& fmt,bool down);
  std::string GetFmtFloat()const{ return(FmtFloat); };
  std::string GetFmtDouble()const{ return(FmtDouble); };

  size_t GetSizeDataConst(bool si64,bool all)const;
  size_t SaveDataConst(bool si64,size_t size,byte* ptr,bool all)const;

  size_t GetSizeData(bool si64,bool all);
  size_t SaveData(bool si64,size_t size,byte* ptr,bool all);
  void LoadData(bool si64,size_t size,const byte* ptr);

  void SaveFile(const std::string& file,bool memory=false,bool all=true);
  void LoadFile(const std::string& file,const std::string& filecode="",bool memory=false);

  void SaveFileListApp(const std::string& file,const std::string& filecode,bool memory=false,bool all=true);
  void LoadFileListApp(const std::string& file,const std::string& filecode,bool memory=false);
  
  void OpenFileStructure(const std::string& file,const std::string& filecode="");
  void CloseFileStructure();
  std::ifstream* GetFileStructure()const;

  void SaveFileXml(std::string file,bool svarrays=false,const std::string& head=" fmt=\"JBinaryData\"")const;

  //-Gestion de items. Management of items.
  JBinaryData* GetParent(){ return(Parent); };
  JBinaryData* GetItemRoot();
  size_t GetItemsCount()const{ return(Items.size()); }
  size_t GetVisibleItemsCount()const;
  size_t GetItemIndex64(const std::string& name);
  JBinaryData* GetItem(const std::string& name);
  JBinaryData* GetItem(size_t index);
  JBinaryData* CreateItem(const std::string& name);
  void RemoveItem(const std::string& name);
  void RemoveItems();
  
  //-Gestion de Arrays. Management of Arrays.
  size_t GetArraysCount()const{ return(Arrays.size()); }
  size_t GetVisibleArraysCount()const;
  size_t GetArrayIndex64(const std::string& name)const;
  JBinaryDataArray* GetArray(const std::string& name);
  JBinaryDataArray* GetArray(size_t index);
  JBinaryDataArray* GetArrayType(const std::string& name,JBinaryDataDef::TpData type,std::string filerror="");
  JBinaryDataArray* GetArrayTpSize(const std::string& name,JBinaryDataDef::TpData type,size_t count,std::string filerror="");
  JBinaryDataArray* CreateArray(const std::string& name,JBinaryDataDef::TpData type);
  JBinaryDataArray* CreateArray(const std::string& name,JBinaryDataDef::TpData type,size_t count,const void* data,bool externalpointer);

  word*     CreateArrayUshort (const std::string& name,size_t count,bool clear);
  unsigned* CreateArrayUint   (const std::string& name,size_t count,bool clear);
  float*    CreateArrayFloat  (const std::string& name,size_t count,bool clear);
  double*   CreateArrayDouble (const std::string& name,size_t count,bool clear);
  tfloat3*  CreateArrayFloat3 (const std::string& name,size_t count,bool clear);
  tdouble3* CreateArrayDouble3(const std::string& name,size_t count,bool clear);

  void RemoveArray(const std::string& name);
  void RemoveArrays();
  JBinaryDataArray* CheckCopyArrayData(const std::string& name,size_t size,JBinaryDataDef::TpData type);
  void CopyArrayData(const std::string& name,size_t size,char*     ptr);
  void CopyArrayData(const std::string& name,size_t size,byte*     ptr);
  void CopyArrayData(const std::string& name,size_t size,short*    ptr);
  void CopyArrayData(const std::string& name,size_t size,word*     ptr);
  void CopyArrayData(const std::string& name,size_t size,int*      ptr);
  void CopyArrayData(const std::string& name,size_t size,unsigned* ptr);
  void CopyArrayData(const std::string& name,size_t size,llong*    ptr);
  void CopyArrayData(const std::string& name,size_t size,ullong*   ptr);
  void CopyArrayData(const std::string& name,size_t size,float*    ptr);
  void CopyArrayData(const std::string& name,size_t size,double*   ptr);
  void CopyArrayData(const std::string& name,size_t size,tint3*    ptr);
  void CopyArrayData(const std::string& name,size_t size,tuint3*   ptr);
  void CopyArrayData(const std::string& name,size_t size,tfloat3*  ptr);
  void CopyArrayData(const std::string& name,size_t size,tdouble3* ptr);

  //-Gestion de values. Management of values.
  size_t GetValuesCount()const{ return(Values.size()); }
  size_t GetValueIndex64(const std::string& name)const;
  std::string NameOfValue(size_t index)const;
  JBinaryDataDef::TpData TypeOfValue(const std::string& name)const;
  JBinaryDataDef::TpData TypeOfValue(size_t index)const;
  bool ExistsValue(const std::string& name)const;
  bool ExistsValue(const std::string& name,JBinaryDataDef::TpData type)const;
  void RemoveValue(const std::string& name);
  void RemoveValues();
  
  std::string GetvText   (const std::string& name,bool optional=false,std::string valdef="")const;
  bool        GetvBool   (const std::string& name,bool optional=false,bool valdef=false)const;
  char        GetvChar   (const std::string& name,bool optional=false,char valdef=0)const;
  byte        GetvUchar  (const std::string& name,bool optional=false,byte valdef=0)const;
  short       GetvShort  (const std::string& name,bool optional=false,short valdef=0)const;
  word        GetvUshort (const std::string& name,bool optional=false,word valdef=0)const;
  int         GetvInt    (const std::string& name,bool optional=false,int valdef=0)const;
  unsigned    GetvUint   (const std::string& name,bool optional=false,unsigned valdef=0)const;
  llong       GetvLlong  (const std::string& name,bool optional=false,llong valdef=0)const;
  ullong      GetvUllong (const std::string& name,bool optional=false,ullong valdef=0)const;
  float       GetvFloat  (const std::string& name,bool optional=false,float valdef=0)const;
  double      GetvDouble (const std::string& name,bool optional=false,double valdef=0)const;
  tint3       GetvInt3   (const std::string& name,bool optional=false,tint3 valdef=TInt3(0))const;
  tuint3      GetvUint3  (const std::string& name,bool optional=false,tuint3 valdef=TUint3(0))const;
  tfloat3     GetvFloat3 (const std::string& name,bool optional=false,tfloat3 valdef=TFloat3(0))const;
  tdouble3    GetvDouble3(const std::string& name,bool optional=false,tdouble3 valdef=TDouble3(0))const;

  void SetvText   (const std::string& name,const std::string& v);
  void SetvBool   (const std::string& name,bool v);
  void SetvChar   (const std::string& name,char v);
  void SetvUchar  (const std::string& name,byte v);
  void SetvShort  (const std::string& name,short v);
  void SetvUshort (const std::string& name,word v);
  void SetvInt    (const std::string& name,int v);
  void SetvUint   (const std::string& name,unsigned v);
  void SetvLlong  (const std::string& name,llong v);
  void SetvUllong (const std::string& name,ullong v);
  void SetvFloat  (const std::string& name,float v);
  void SetvDouble (const std::string& name,double v);
  void SetvInt3   (const std::string& name,tint3 v);
  void SetvUint3  (const std::string& name,tuint3 v);
  void SetvFloat3 (const std::string& name,tfloat3 v);
  void SetvDouble3(const std::string& name,tdouble3 v);
};

/*
Structure of file JBinaryData:
===================================
StHeadFmtBin header
uint size_item_def
- "ITEM"
- str name
- bool hide
- bool hidevalues
- str fmtfloat
- str fmtdouble
- uint num_arrays
- uint num_items
- uint size_values
  - "VALUES"
  - uint num_values
  - [value_0]
      ...
  - [value_n]
[array_0]    
uint size_array_def
- "ARRAY"
- str name
- bool hide
- int type
- (si64? ullong: uint) count  
- (si64? ullong: uint) size_data
  - [contenido de array]
[array_1]    
... 
[array_n]    
*/



#endif


