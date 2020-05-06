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
//:# Clase base para la definicion de cualquier formato binario de fichero.
//:# Algunas de sus funcionalidades son:
//:# - Permite guardar datos de todo tipo en formato binario.
//:# - Todos los datos contienen nombre y tipo para permitir su exploracion sin
//:#   sin conocer su estructura.
//:# - Los datos se pueden organizar en distintos niveles con forma de arbol.
//:# - Permite guardar datos basicos y arrays de tipos basicos y strings.
//:# - La lectura del contenido de arrays desde fichero puede ser de forma 
//:#   selectiva, sólo los que se necesiten.
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
//:#############################################################################

/// \file JBinaryData.h \brief Declares the class \ref JBinaryData.

#ifndef _JBinaryData_
#define _JBinaryData_

#include "JObject.h"
#include "TypesDef.h"
#include <string>
#include <vector>
#include <fstream>

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

  static std::string TypeToStr(TpData type);
  static size_t SizeOfType(TpData type);
  static bool TypeIsTriple(TpData type);
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
  unsigned Count;         ///<Numero de elementos almacenados en pointer. Number of elements stored in pointer.
  unsigned Size;          ///<Numero de elementos para los que hay memoria reservada. Number of elements for which there is reserved memory.
  void* Pointer;
  bool ExternalPointer;   ///<Indica que el puntero es externo y no debe liberarse. Indicates that the pointer is external, and should not be released.
  llong FileDataPos;      ///<Valor mayor o igual a cero indica la posicion de lectura en el fichero abierto en el ItemHead. Value greater than or equal to zero indicates the position of reading in the file opened in the ItemHead.
  unsigned FileDataCount; ///<Numero de elemetos del array en fichero. Number of elements in the array in a file.
  unsigned FileDataSize;  ///<Size de datos del array en fichero. Size of array data in file.

  void FreePointer(void* ptr)const;
  void* AllocPointer(unsigned size)const;
  void CheckMemory(unsigned count,bool resize);

  void        OutData(unsigned &count,unsigned size,const byte *ptr,byte *dat,unsigned sdat)const;
  std::string OutStr (unsigned &count,unsigned size,const byte *ptr)const;
  unsigned    OutUint(unsigned &count,unsigned size,const byte *ptr)const{  unsigned v; OutData(count,size,ptr,(byte*)&v,sizeof(unsigned)); return(v);  }  /// Extrae unsigned de ptr. Extract ptr unsigned.

 public:
  JBinaryDataArray(JBinaryData* parent,const std::string &name,JBinaryDataDef::TpData type);
  ~JBinaryDataArray();
  llong GetAllocMemory()const;

  void SetName(const std::string &name);
  std::string GetName()const{ return(Name); };

  void SetHide(bool hide){ Hide=hide; }
  bool GetHide()const{ return(Hide); }

  JBinaryDataDef::TpData GetType()const{ return(Type); };
  unsigned GetCount()const{ return(Count); };
  unsigned GetSize()const{ return(Size); };
  const void* GetPointer()const{ return(Pointer); };

  bool PointerIsExternal()const{ return(ExternalPointer); };
  bool DataInPointer()const{ return(Pointer&&Count); }
  bool DataInFile()const{ return(FileDataPos>=0); }

  void FreeMemory();
  void AllocMemory(unsigned size,bool savedata=false);
  void ConfigExternalMemory(unsigned size,void* pointer);

  void ReadData(unsigned count,unsigned size,std::ifstream *pf,bool resize);
  void AddData(unsigned count,const void* data,bool resize);
  void SetData(unsigned count,const void* data,bool externalpointer);

  const void* GetDataPointer()const;
  unsigned GetDataCopy(unsigned size,void* pointer)const;

  void AddText(const std::string &str,bool resize);
  void AddTexts(unsigned count,const std::string *strs,bool resize);

  void ConfigFileData(llong filepos,unsigned datacount,unsigned datasize);
  void ClearFileData();
  unsigned GetFileDataCount()const{ return(FileDataCount); }
  unsigned GetFileDataSize()const{ return(FileDataSize); }
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
    byte byteorder;              ///<1:BigEndian 0:LittleEndian.
    byte void1;                  ///<Not used.
    byte void2;                  ///<Not used.
    byte void3;                  ///<Not used.
  }StHeadFmtBin;//-sizeof(64)

  static const std::string CodeItemDef;
  static const std::string CodeValuesDef;
  static const std::string CodeArrayDef;

 public:

  ///Structure that describes the information of a value.
  typedef struct{
    std::string name;
    JBinaryDataDef::TpData type;
    std::string vtext;
    union{
      char vchar;
      unsigned char vuchar;
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

  std::ifstream *FileStructure;

  //-Variables para cache de values. Variables to cache values.
  bool ValuesModif;
  byte* ValuesData;
  unsigned ValuesSize;

  //-Gestion de Values. Management of Values.
  void ValuesCacheReset();
  int CheckGetValue(const std::string &name,bool optional,JBinaryDataDef::TpData type)const;
  int CheckSetValue(const std::string &name,JBinaryDataDef::TpData type);
  static void ResetValue(const std::string &name,JBinaryDataDef::TpData type,StValue &v);
  std::string ValueToXml(const StValue &v)const;

  void InData   (unsigned &count,unsigned size,byte *ptr,const byte *dat,unsigned sdat)const;
  void InStr    (unsigned &count,unsigned size,byte *ptr,const std::string &cad)const;
  void InBool   (unsigned &count,unsigned size,byte *ptr,bool v)const{            int vv=(v? 1: 0); InInt(count,size,ptr,vv);               }  ///<Introduce bool en ptr. Introduces bool in ptr.
  void InChar   (unsigned &count,unsigned size,byte *ptr,char v)const{            InData(count,size,ptr,(byte*)&v,sizeof(char));            }  ///<Introduce char en ptr. Introduces char in ptr.
  void InUchar  (unsigned &count,unsigned size,byte *ptr,unsigned char v)const{   InData(count,size,ptr,(byte*)&v,sizeof(unsigned char));   }  ///<Introduce unsigned char en ptr. Introduces unsigned char in ptr.
  void InShort  (unsigned &count,unsigned size,byte *ptr,short v)const{           InData(count,size,ptr,(byte*)&v,sizeof(short));           }  ///<Introduce short en ptr. Introduces short in ptr.
  void InUshort (unsigned &count,unsigned size,byte *ptr,unsigned short v)const{  InData(count,size,ptr,(byte*)&v,sizeof(unsigned short));  }  ///<Introduce unsigned short en ptr. Introduces unsigned short in ptr.
  void InInt    (unsigned &count,unsigned size,byte *ptr,int v)const{             InData(count,size,ptr,(byte*)&v,sizeof(int));             }  ///<Introduce int en ptr. Introduces int in ptr.
  void InUint   (unsigned &count,unsigned size,byte *ptr,unsigned v)const{        InData(count,size,ptr,(byte*)&v,sizeof(unsigned));        }  ///<Introduce unsigned en ptr. Introduces unsigned in ptr
  void InLlong  (unsigned &count,unsigned size,byte *ptr,llong v)const{           InData(count,size,ptr,(byte*)&v,sizeof(llong));           }  ///<Introduce long long en ptr. Introduces long long in ptr
  void InUllong (unsigned &count,unsigned size,byte *ptr,ullong v)const{          InData(count,size,ptr,(byte*)&v,sizeof(ullong));          }  ///<Introduce unsigned long long en ptr. Introduces unsigned long in ptr
  void InFloat  (unsigned &count,unsigned size,byte *ptr,float v)const{           InData(count,size,ptr,(byte*)&v,sizeof(float));           }  ///<Introduce float en ptr. Introduces float in ptr
  void InDouble (unsigned &count,unsigned size,byte *ptr,double v)const{          InData(count,size,ptr,(byte*)&v,sizeof(double));          }  ///<Introduce double en ptr. Introduces double in ptr
  void InInt3   (unsigned &count,unsigned size,byte *ptr,tint3 v)const{           InData(count,size,ptr,(byte*)&v,sizeof(tint3));           }  ///<Introduce tint3 en ptr. Introduces tint3 in ptr
  void InUint3  (unsigned &count,unsigned size,byte *ptr,tuint3 v)const{          InData(count,size,ptr,(byte*)&v,sizeof(tuint3));          }  ///<Introduce tuint3 en ptr. Introduces tuint3 in ptr
  void InFloat3 (unsigned &count,unsigned size,byte *ptr,tfloat3 v)const{         InData(count,size,ptr,(byte*)&v,sizeof(tfloat3));         }  ///<Introduce tfloat3 en ptr. Introduces tfloat3 in ptr
  void InDouble3(unsigned &count,unsigned size,byte *ptr,tdouble3 v)const{        InData(count,size,ptr,(byte*)&v,sizeof(tdouble3));        }  ///<Introduce tdouble3 en ptr. Introduces tdouble3 in ptr

  void           OutData   (unsigned &count,unsigned size,const byte *ptr,byte *dat,unsigned sdat)const;
  std::string    OutStr    (unsigned &count,unsigned size,const byte *ptr)const;
  bool           OutBool   (unsigned &count,unsigned size,const byte *ptr)const{  return(OutInt(count,size,ptr)!=0);  }  /// Extrae bool de ptr.
  char           OutChar   (unsigned &count,unsigned size,const byte *ptr)const{  char v;           OutData(count,size,ptr,(byte*)&v,sizeof(char));           return(v);  }  ///< Extrae char de ptr. Extracts char of ptr.
  unsigned char  OutUchar  (unsigned &count,unsigned size,const byte *ptr)const{  unsigned char v;  OutData(count,size,ptr,(byte*)&v,sizeof(unsigned char));  return(v);  }  ///< Extrae unsigned char de ptr. Extracts unsigned char of ptr.
  short          OutShort  (unsigned &count,unsigned size,const byte *ptr)const{  short v;          OutData(count,size,ptr,(byte*)&v,sizeof(short));          return(v);  }  ///< Extrae short de ptr. Extracts short of ptr.
  unsigned short OutUshort (unsigned &count,unsigned size,const byte *ptr)const{  unsigned short v; OutData(count,size,ptr,(byte*)&v,sizeof(unsigned short)); return(v);  }  ///< Extrae unsigned short de ptr. Extracts unsigned short of ptr.
  int            OutInt    (unsigned &count,unsigned size,const byte *ptr)const{  int v;            OutData(count,size,ptr,(byte*)&v,sizeof(int));            return(v);  }  ///< Extrae int de ptr. Extracts int of ptr.
  unsigned       OutUint   (unsigned &count,unsigned size,const byte *ptr)const{  unsigned v;       OutData(count,size,ptr,(byte*)&v,sizeof(unsigned));       return(v);  }  ///< Extrae unsigned de ptr. Extracts unsigned of ptr.
  llong          OutLlong  (unsigned &count,unsigned size,const byte *ptr)const{  llong v;          OutData(count,size,ptr,(byte*)&v,sizeof(llong));          return(v);  }  ///< Extrae long long de ptr. Extracts long long of ptr.
  ullong         OutUllong (unsigned &count,unsigned size,const byte *ptr)const{  ullong v;         OutData(count,size,ptr,(byte*)&v,sizeof(ullong));         return(v);  }  ///< Extrae unsigned long long de ptr. Extracts unsigned long long of ptr.
  float          OutFloat  (unsigned &count,unsigned size,const byte *ptr)const{  float v;          OutData(count,size,ptr,(byte*)&v,sizeof(float));          return(v);  }  ///< Extrae float de ptr. Extracts float of ptr.
  double         OutDouble (unsigned &count,unsigned size,const byte *ptr)const{  double v;         OutData(count,size,ptr,(byte*)&v,sizeof(double));         return(v);  }  ///< Extrae double de ptr. Extracts double of ptr.
  tint3          OutInt3   (unsigned &count,unsigned size,const byte *ptr)const{  tint3 v;          OutData(count,size,ptr,(byte*)&v,sizeof(tint3));          return(v);  }  ///< Extrae tint3 de ptr. Extracts tint3 of ptr.
  tuint3         OutUint3  (unsigned &count,unsigned size,const byte *ptr)const{  tuint3 v;         OutData(count,size,ptr,(byte*)&v,sizeof(tuint3));         return(v);  }  /// Extrae tuint3 de ptr. Extracts tuint3 of ptr.
  tfloat3        OutFloat3 (unsigned &count,unsigned size,const byte *ptr)const{  tfloat3 v;        OutData(count,size,ptr,(byte*)&v,sizeof(tfloat3));        return(v);  }  /// Extrae tfloat3 de ptr. Extracts tfloat3 of ptr.
  tdouble3       OutDouble3(unsigned &count,unsigned size,const byte *ptr)const{  tdouble3 v;       OutData(count,size,ptr,(byte*)&v,sizeof(tdouble3));       return(v);  }  /// Extrae tdouble3 de ptr. Extracts tdouble3 of ptr.

  void InValue(unsigned &count,unsigned size,byte *ptr,const StValue &v)const;
  void OutValue(unsigned &count,unsigned size,const byte *ptr);

  void InArrayBase(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const;
  void InArrayData(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const;
  void InArray(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const;
  void InItemBase(unsigned &count,unsigned size,byte *ptr,bool all)const;
  void InItem(unsigned &count,unsigned size,byte *ptr,bool all)const;

  JBinaryDataArray* OutArrayBase(unsigned &count,unsigned size,const byte *ptr,unsigned &countdata,unsigned &sizedata);
  void OutArrayData(unsigned &count,unsigned size,const byte *ptr,JBinaryDataArray *ar,unsigned countdata,unsigned sizedata);
  void OutArray(unsigned &count,unsigned size,const byte *ptr);
  JBinaryData* OutItemBase(unsigned &count,unsigned size,const byte *ptr,bool create,unsigned &narrays,unsigned &nitems,unsigned &sizevalues);
  void OutItem(unsigned &count,unsigned size,const byte *ptr,bool create);

  unsigned GetSizeValues()const;
  void SaveValues(unsigned &count,unsigned size,byte *ptr)const;
  void ValuesCachePrepare(bool down);

  void WriteArrayData(std::fstream *pf,const JBinaryDataArray *ar)const;
  void WriteArray(std::fstream *pf,unsigned sbuf,byte *buf,const JBinaryDataArray *ar)const;
  void WriteItem(std::fstream *pf,unsigned sbuf,byte *buf,bool all)const;

  unsigned ReadUint(std::ifstream *pf)const;
  void ReadArrayData(std::ifstream *pf,JBinaryDataArray *ar,unsigned countdata,unsigned sizedata,bool loadarraysdata);
  void ReadArray(std::ifstream *pf,unsigned sbuf,byte *buf,bool loadarraysdata);
  void ReadItem(std::ifstream *pf,unsigned sbuf,byte *buf,bool create,bool loadarraysdata);

  JBinaryData::StHeadFmtBin MakeFileHead(const std::string &filecode)const;
  unsigned GetFileHead(std::ifstream *pf,JBinaryData::StHeadFmtBin &head)const;
  void CheckHead(const std::string &file,const StHeadFmtBin &head,const std::string &filecode)const;
  unsigned CheckFileHead(const std::string &file,std::ifstream *pf,const std::string &filecode)const;
  unsigned CheckFileListHead(const std::string &file,std::fstream *pf,const std::string &filecode)const;
  void SaveFileData(std::fstream *pf,bool head,const std::string &filecode,bool memory,bool all)const;

  void WriteFileXmlArray(const std::string &tabs,std::ofstream* pf,bool svarrays,const JBinaryDataArray* ar)const;

  void WriteFileXml(const std::string &tabs,std::ofstream* pf,bool svarrays)const;

 public:
  JBinaryData(std::string name="JBinary_Data");
  JBinaryData(const JBinaryData &src);
  ~JBinaryData();
  JBinaryData& operator=(const JBinaryData &src);
  void Clear();
  llong GetAllocMemory()const;

  void SetName(const std::string &name);
  std::string GetName()const{ return(Name); };
  
  void SetHide(bool hide){ HideAll=hide; }
  bool GetHide()const{ return(HideAll); }
  void SetHideValues(bool hide,bool down);
  bool GetHideValues()const{ return(HideValues); }
  void SetHideArrays(bool hide,bool down);
  void SetHideItems(bool hide,bool down);

  void SetFmtFloat(const std::string &fmt,bool down);
  void SetFmtDouble(const std::string &fmt,bool down);
  std::string GetFmtFloat()const{ return(FmtFloat); };
  std::string GetFmtDouble()const{ return(FmtDouble); };

  unsigned GetSizeDataConst(bool all)const;
  unsigned SaveDataConst(unsigned size,byte* ptr,bool all)const;

  unsigned GetSizeData(bool all);
  unsigned SaveData(unsigned size,byte* ptr,bool all);
  void LoadData(unsigned size,const byte* ptr);

  void SaveFile(const std::string &file,bool memory=false,bool all=true);
  void LoadFile(const std::string &file,const std::string &filecode="",bool memory=false);

  void SaveFileListApp(const std::string &file,const std::string &filecode,bool memory=false,bool all=true);
  void LoadFileListApp(const std::string &file,const std::string &filecode,bool memory=false);
  
  void OpenFileStructure(const std::string &file,const std::string &filecode="");
  void CloseFileStructure();
  std::ifstream* GetFileStructure()const;

  void SaveFileXml(std::string file,bool svarrays=false,const std::string &head=" fmt=\"JBinaryData\"")const;

  //-Gestion de items. Management of items.
  JBinaryData* GetParent(){ return(Parent); };
  JBinaryData* GetItemRoot();
  unsigned GetItemsCount()const{ return(unsigned(Items.size())); }
  unsigned GetVisibleItemsCount()const;
  int GetItemIndex(const std::string &name);
  JBinaryData* GetItem(const std::string &name);
  JBinaryData* GetItem(unsigned index);
  JBinaryData* CreateItem(const std::string &name);
  void RemoveItem(const std::string &name);
  void RemoveItems();
  
  //-Gestion de Arrays. Management of Arrays.
  unsigned GetArraysCount()const{ return(unsigned(Arrays.size())); }
  unsigned GetVisibleArraysCount()const;
  int GetArrayIndex(const std::string &name)const;
  JBinaryDataArray* GetArray(const std::string &name);
  JBinaryDataArray* GetArray(unsigned index);
  JBinaryDataArray* CreateArray(const std::string &name,JBinaryDataDef::TpData type);
  JBinaryDataArray* CreateArray(const std::string &name,JBinaryDataDef::TpData type,unsigned count,const void *data,bool externalpointer);
  void RemoveArray(const std::string &name);
  void RemoveArrays();
  JBinaryDataArray* CheckCopyArrayData(const std::string &name,unsigned size,JBinaryDataDef::TpData type);
  void CopyArrayData(const std::string &name,unsigned size,char           *ptr);
  void CopyArrayData(const std::string &name,unsigned size,unsigned char  *ptr);
  void CopyArrayData(const std::string &name,unsigned size,short          *ptr);
  void CopyArrayData(const std::string &name,unsigned size,unsigned short *ptr);
  void CopyArrayData(const std::string &name,unsigned size,int            *ptr);
  void CopyArrayData(const std::string &name,unsigned size,unsigned       *ptr);
  void CopyArrayData(const std::string &name,unsigned size,llong          *ptr);
  void CopyArrayData(const std::string &name,unsigned size,ullong         *ptr);
  void CopyArrayData(const std::string &name,unsigned size,float          *ptr);
  void CopyArrayData(const std::string &name,unsigned size,double         *ptr);
  void CopyArrayData(const std::string &name,unsigned size,tint3          *ptr);
  void CopyArrayData(const std::string &name,unsigned size,tuint3         *ptr);
  void CopyArrayData(const std::string &name,unsigned size,tfloat3        *ptr);
  void CopyArrayData(const std::string &name,unsigned size,tdouble3       *ptr);

  //-Gestion de values. Management of values.
  unsigned GetValuesCount()const{ return(unsigned(Values.size())); }
  int GetValueIndex(const std::string &name)const;
  std::string NameOfValue(unsigned index)const;
  JBinaryDataDef::TpData TypeOfValue(const std::string &name)const;
  JBinaryDataDef::TpData TypeOfValue(unsigned index)const;
  bool ExistsValue(const std::string &name)const;
  bool ExistsValue(const std::string &name,JBinaryDataDef::TpData type)const;
  void RemoveValue(const std::string &name);
  void RemoveValues();
  
  std::string    GetvText   (const std::string &name,bool optional=false,std::string valdef="")const;
  bool           GetvBool   (const std::string &name,bool optional=false,bool valdef=false)const;
  char           GetvChar   (const std::string &name,bool optional=false,char valdef=0)const;
  unsigned char  GetvUchar  (const std::string &name,bool optional=false,unsigned char valdef=0)const;
  short          GetvShort  (const std::string &name,bool optional=false,short valdef=0)const;
  unsigned short GetvUshort (const std::string &name,bool optional=false,unsigned short valdef=0)const;
  int            GetvInt    (const std::string &name,bool optional=false,int valdef=0)const;
  unsigned       GetvUint   (const std::string &name,bool optional=false,unsigned valdef=0)const;
  llong          GetvLlong  (const std::string &name,bool optional=false,llong valdef=0)const;
  ullong         GetvUllong (const std::string &name,bool optional=false,ullong valdef=0)const;
  float          GetvFloat  (const std::string &name,bool optional=false,float valdef=0)const;
  double         GetvDouble (const std::string &name,bool optional=false,double valdef=0)const;
  tint3          GetvInt3   (const std::string &name,bool optional=false,tint3 valdef=TInt3(0))const;
  tuint3         GetvUint3  (const std::string &name,bool optional=false,tuint3 valdef=TUint3(0))const;
  tfloat3        GetvFloat3 (const std::string &name,bool optional=false,tfloat3 valdef=TFloat3(0))const;
  tdouble3       GetvDouble3(const std::string &name,bool optional=false,tdouble3 valdef=TDouble3(0))const;

  void SetvText   (const std::string &name,const std::string &v);
  void SetvBool   (const std::string &name,bool v);
  void SetvChar   (const std::string &name,char v);
  void SetvUchar  (const std::string &name,unsigned char v);
  void SetvShort  (const std::string &name,short v);
  void SetvUshort (const std::string &name,unsigned short v);
  void SetvInt    (const std::string &name,int v);
  void SetvUint   (const std::string &name,unsigned v);
  void SetvLlong  (const std::string &name,llong v);
  void SetvUllong (const std::string &name,ullong v);
  void SetvFloat  (const std::string &name,float v);
  void SetvDouble (const std::string &name,double v);
  void SetvInt3   (const std::string &name,tint3 v);
  void SetvUint3  (const std::string &name,tuint3 v);
  void SetvFloat3 (const std::string &name,tfloat3 v);
  void SetvDouble3(const std::string &name,tdouble3 v);
};

/*
Structure of file JBinaryData:
===================================
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
uint size_array_def [n]
- "ARRAY"
- str name
- bool hide
- int type
- uint count
- uint size_contenido
  - [contenido de array]
[array_1]    
... 
[array_n]    
*/



#endif


