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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestiona vector de arrays para intercambio de datos entre funciones o
//:#   generacion de ficheros de salida. (23-08-2019)
//:# - Nuevos metodos para crear arrays con memoria dinamica de forma automatica. (05-09-2019)
//:# - Mejora la gestion de excepciones. (14-09-2019)
//:# - Nuevos metodos GetDataCount(), GetArrayCte(), GetArrayDim()... (11-12-2019)
//:# - Nuevos metodos FilterApply(), FilterList(), SortData(), FilterSortList()... (18-01-2020)
//:# - Nuevos tipos int e int3. (19-02-2020)
//:# - Nuevos metodos CreateArrayPtrXXX(). (04-03-2020)
//:#############################################################################

/// \file JDataArrays.h \brief Declares the class \ref JDataArrays.

#ifndef _JDataArrays_
#define _JDataArrays_

#include "JObject.h"
#include "TypesDef.h"
#include <string>
#include <vector>
#include <climits>

//-Defines for normal exceptions for static methods.
#ifndef Run_ExceptioonSta
#define Run_ExceptioonSta(msg) RunExceptioonStatic(__FILE__,__LINE__,__func__,msg)
#endif
#ifndef Run_ExceptioonFileSta
#define Run_ExceptioonFileSta(msg,file) RunExceptioonStatic(__FILE__,__LINE__,__func__,msg,file)
#endif

//##############################################################################
//# JDataArrays
//##############################################################################
/// \brief Manage a vector of arrays of data.

class JDataArrays : protected JObject
{
protected:
  static void RunExceptioonStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,const std::string &msg,const std::string &file="");

public:
  typedef struct {
    std::string fullname; ///<Name + other information. "name:outputformat:units", e.g. "velocity:%f:m/s"
    std::string keyname;  ///<Name for identification.
    TpTypeData type;      ///<Type of values.
    unsigned count;       ///<Number of values in pointer.
    void *ptr;            ///<Data pointer [count].
    bool delptr;          ///<Automatically frees memory in destructor or Reset() method.
  }StDataArray;

protected:
  std::vector<StDataArray> Arrays;

  void FreeMemory(StDataArray &arr);
  void FreeMemory();

  template<class T> void TReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,T *data,T *aux)const;
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,byte     *data,byte     *aux)const{ TReindexData<byte>    (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,word     *data,word     *aux)const{ TReindexData<word>    (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,unsigned *data,unsigned *aux)const{ TReindexData<unsigned>(sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,int      *data,int      *aux)const{ TReindexData<int>     (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,float    *data,float    *aux)const{ TReindexData<float>   (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,double   *data,double   *aux)const{ TReindexData<double>  (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,tuint3   *data,tuint3   *aux)const{ TReindexData<tuint3>  (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,tint3    *data,tint3    *aux)const{ TReindexData<tint3>   (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,tfloat3  *data,tfloat3  *aux)const{ TReindexData<tfloat3> (sreindex,reindex,ndata,data,aux); }
  void ReindexData(unsigned sreindex,const unsigned *reindex,unsigned ndata,tdouble3 *data,tdouble3 *aux)const{ TReindexData<tdouble3>(sreindex,reindex,ndata,data,aux); }

  unsigned SortData(unsigned count,const unsigned *reindex);

public:
  JDataArrays();
  ~JDataArrays();
  void Reset();

  void CopyFrom(const JDataArrays &arr);

  unsigned AddArray(std::string fullname,TpTypeData type,unsigned count,void *ptr,bool delptr);
  unsigned AddArray(std::string fullname,unsigned count,const byte     *ptr,bool delptr=false){ return(AddArray(fullname,TypeUchar  ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const word     *ptr,bool delptr=false){ return(AddArray(fullname,TypeUshort ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const unsigned *ptr,bool delptr=false){ return(AddArray(fullname,TypeUint   ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const int      *ptr,bool delptr=false){ return(AddArray(fullname,TypeInt    ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const float    *ptr,bool delptr=false){ return(AddArray(fullname,TypeFloat  ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const double   *ptr,bool delptr=false){ return(AddArray(fullname,TypeDouble ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const tuint3   *ptr,bool delptr=false){ return(AddArray(fullname,TypeUint3  ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const tint3    *ptr,bool delptr=false){ return(AddArray(fullname,TypeInt3   ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const tfloat3  *ptr,bool delptr=false){ return(AddArray(fullname,TypeFloat3 ,count,(void*)ptr,delptr)); }
  unsigned AddArray(std::string fullname,unsigned count,const tdouble3 *ptr,bool delptr=false){ return(AddArray(fullname,TypeDouble3,count,(void*)ptr,delptr)); }

  unsigned CreateArrayByte   (std::string fullname,unsigned count,byte     value=0          ,bool delptr=true){ return(AddArray(fullname,count,NewArrayByte   (count,true,value),delptr)); }
  unsigned CreateArrayWord   (std::string fullname,unsigned count,word     value=0          ,bool delptr=true){ return(AddArray(fullname,count,NewArrayWord   (count,true,value),delptr)); }
  unsigned CreateArrayUint   (std::string fullname,unsigned count,unsigned value=0          ,bool delptr=true){ return(AddArray(fullname,count,NewArrayUint   (count,true,value),delptr)); }
  unsigned CreateArrayInt    (std::string fullname,unsigned count,int      value=0          ,bool delptr=true){ return(AddArray(fullname,count,NewArrayInt    (count,true,value),delptr)); }
  unsigned CreateArrayFloat  (std::string fullname,unsigned count,float    value=0          ,bool delptr=true){ return(AddArray(fullname,count,NewArrayFloat  (count,true,value),delptr)); }
  unsigned CreateArrayDouble (std::string fullname,unsigned count,double   value=0          ,bool delptr=true){ return(AddArray(fullname,count,NewArrayDouble (count,true,value),delptr)); }
  unsigned CreateArrayUint3  (std::string fullname,unsigned count,tuint3   value=TUint3(0)  ,bool delptr=true){ return(AddArray(fullname,count,NewArrayUint3  (count,true,value),delptr)); }
  unsigned CreateArrayInt3   (std::string fullname,unsigned count,tint3    value=TInt3(0)   ,bool delptr=true){ return(AddArray(fullname,count,NewArrayInt3   (count,true,value),delptr)); }
  unsigned CreateArrayFloat3 (std::string fullname,unsigned count,tfloat3  value=TFloat3(0) ,bool delptr=true){ return(AddArray(fullname,count,NewArrayFloat3 (count,true,value),delptr)); }
  unsigned CreateArrayDouble3(std::string fullname,unsigned count,tdouble3 value=TDouble3(0),bool delptr=true){ return(AddArray(fullname,count,NewArrayDouble3(count,true,value),delptr)); }

  byte*     CreateArrayPtrByte   (std::string fullname,unsigned count,byte     value=0          ,bool delptr=true){ return((byte*)    GetArrayByte   (CreateArrayByte   (fullname,count,value,delptr))); }
  word*     CreateArrayPtrWord   (std::string fullname,unsigned count,word     value=0          ,bool delptr=true){ return((word*)    GetArrayWord   (CreateArrayWord   (fullname,count,value,delptr))); }
  unsigned* CreateArrayPtrUint   (std::string fullname,unsigned count,unsigned value=0          ,bool delptr=true){ return((unsigned*)GetArrayUint   (CreateArrayUint   (fullname,count,value,delptr))); }
  int*      CreateArrayPtrInt    (std::string fullname,unsigned count,int      value=0          ,bool delptr=true){ return((int*)     GetArrayInt    (CreateArrayInt    (fullname,count,value,delptr))); }
  float*    CreateArrayPtrFloat  (std::string fullname,unsigned count,float    value=0          ,bool delptr=true){ return((float*)   GetArrayFloat  (CreateArrayFloat  (fullname,count,value,delptr))); }
  double*   CreateArrayPtrDouble (std::string fullname,unsigned count,double   value=0          ,bool delptr=true){ return((double*)  GetArrayDouble (CreateArrayDouble (fullname,count,value,delptr))); }
  tuint3*   CreateArrayPtrUint3  (std::string fullname,unsigned count,tuint3   value=TUint3(0)  ,bool delptr=true){ return((tuint3*)  GetArrayUint3  (CreateArrayUint3  (fullname,count,value,delptr))); }
  tint3*    CreateArrayPtrInt3   (std::string fullname,unsigned count,tint3    value=TInt3(0)   ,bool delptr=true){ return((tint3*)   GetArrayInt3   (CreateArrayInt3   (fullname,count,value,delptr))); }
  tfloat3*  CreateArrayPtrFloat3 (std::string fullname,unsigned count,tfloat3  value=TFloat3(0) ,bool delptr=true){ return((tfloat3*) GetArrayFloat3 (CreateArrayFloat3 (fullname,count,value,delptr))); }
  tdouble3* CreateArrayPtrDouble3(std::string fullname,unsigned count,tdouble3 value=TDouble3(0),bool delptr=true){ return((tdouble3*)GetArrayDouble3(CreateArrayDouble3(fullname,count,value,delptr))); }

  void DeleteArray(unsigned idx);
  void DeleteArray(std::string keyname);

  void EraseArray(std::string keyname);
  void MoveArray(unsigned idx,unsigned idx2);

  unsigned Count()const{ return(unsigned(Arrays.size())); }

  unsigned GetDataCount(bool minimum)const;
  unsigned GetDataCount()const;

  unsigned GetIdxName(const std::string &keyname)const;
  bool ExistsName(const std::string &keyname)const{ return(GetIdxName(keyname)!=UINT_MAX); }

  std::string CheckErrorArray(unsigned idx,TpTypeData type,unsigned count)const;
  std::string CheckErrorArray(const std::string &keyname,TpTypeData type,unsigned count)const;

  JDataArrays::StDataArray& GetArray(unsigned idx);
  JDataArrays::StDataArray& GetArray(const std::string &keyname);

  const JDataArrays::StDataArray& GetArrayCte(unsigned idx)const;
  const JDataArrays::StDataArray& GetArrayCte(const std::string &keyname)const;

  JDataArrays::StDataArray GetArrayData(unsigned idx)const;
  JDataArrays::StDataArray GetArrayData(const std::string &keyname)const;

  int GetArrayDim(unsigned idx)const;

  std::string GetArrayFmt(unsigned idx)const;
  static std::string GetFmtByType(TpTypeData type);

  std::string GetArrayUnits(unsigned idx)const;
  static std::string GetUnitsByName(std::string keyname);

  const void* GetArrayPtr(unsigned idx,TpTypeData type,unsigned count=0)const;
  const byte*     GetArrayByte   (unsigned idx,unsigned count=0)const{ return((byte*    )GetArrayPtr(idx,TypeUchar  ,count)); }
  const word*     GetArrayWord   (unsigned idx,unsigned count=0)const{ return((word*    )GetArrayPtr(idx,TypeUshort ,count)); }
  const unsigned* GetArrayUint   (unsigned idx,unsigned count=0)const{ return((unsigned*)GetArrayPtr(idx,TypeUint   ,count)); }
  const int*      GetArrayInt    (unsigned idx,unsigned count=0)const{ return((int*     )GetArrayPtr(idx,TypeInt    ,count)); }
  const float*    GetArrayFloat  (unsigned idx,unsigned count=0)const{ return((float*   )GetArrayPtr(idx,TypeFloat  ,count)); }
  const double*   GetArrayDouble (unsigned idx,unsigned count=0)const{ return((double*  )GetArrayPtr(idx,TypeDouble ,count)); }
  const tuint3*   GetArrayUint3  (unsigned idx,unsigned count=0)const{ return((tuint3*  )GetArrayPtr(idx,TypeUint3  ,count)); }
  const tint3*    GetArrayInt3   (unsigned idx,unsigned count=0)const{ return((tint3*   )GetArrayPtr(idx,TypeInt3   ,count)); }
  const tfloat3*  GetArrayFloat3 (unsigned idx,unsigned count=0)const{ return((tfloat3* )GetArrayPtr(idx,TypeFloat3 ,count)); }
  const tdouble3* GetArrayDouble3(unsigned idx,unsigned count=0)const{ return((tdouble3*)GetArrayPtr(idx,TypeDouble3,count)); }

  const void* GetArrayPtr(const std::string &keyname,TpTypeData type,unsigned count=0)const;
  const byte*     GetArrayByte   (const std::string &keyname,unsigned count=0)const{ return((byte*    )GetArrayPtr(keyname,TypeUchar  ,count)); }
  const word*     GetArrayWord   (const std::string &keyname,unsigned count=0)const{ return((word*    )GetArrayPtr(keyname,TypeUshort ,count)); }
  const unsigned* GetArrayUint   (const std::string &keyname,unsigned count=0)const{ return((unsigned*)GetArrayPtr(keyname,TypeUint   ,count)); }
  const int*      GetArrayInt    (const std::string &keyname,unsigned count=0)const{ return((int*     )GetArrayPtr(keyname,TypeInt    ,count)); }
  const float*    GetArrayFloat  (const std::string &keyname,unsigned count=0)const{ return((float*   )GetArrayPtr(keyname,TypeFloat  ,count)); }
  const double*   GetArrayDouble (const std::string &keyname,unsigned count=0)const{ return((double*  )GetArrayPtr(keyname,TypeDouble ,count)); }
  const tuint3*   GetArrayUint3  (const std::string &keyname,unsigned count=0)const{ return((tuint3*  )GetArrayPtr(keyname,TypeUint3  ,count)); }
  const tint3*    GetArrayInt3   (const std::string &keyname,unsigned count=0)const{ return((tint3*   )GetArrayPtr(keyname,TypeInt3   ,count)); }
  const tfloat3*  GetArrayFloat3 (const std::string &keyname,unsigned count=0)const{ return((tfloat3* )GetArrayPtr(keyname,TypeFloat3 ,count)); }
  const tdouble3* GetArrayDouble3(const std::string &keyname,unsigned count=0)const{ return((tdouble3*)GetArrayPtr(keyname,TypeDouble3,count)); }

  void Print()const;

public:
  static byte*     NewArrayByte   (unsigned count,bool defvalue=false,byte     value=0);
  static word*     NewArrayWord   (unsigned count,bool defvalue=false,word     value=0);
  static unsigned* NewArrayUint   (unsigned count,bool defvalue=false,unsigned value=0);
  static int*      NewArrayInt    (unsigned count,bool defvalue=false,int      value=0);
  static float*    NewArrayFloat  (unsigned count,bool defvalue=false,float    value=0);
  static double*   NewArrayDouble (unsigned count,bool defvalue=false,double   value=0);
  static tuint3*   NewArrayUint3  (unsigned count,bool defvalue=false,tuint3   value=TUint3(0));
  static tint3*    NewArrayInt3   (unsigned count,bool defvalue=false,tint3    value=TInt3(0));
  static tfloat3*  NewArrayFloat3 (unsigned count,bool defvalue=false,tfloat3  value=TFloat3(0));
  static tdouble3* NewArrayDouble3(unsigned count,bool defvalue=false,tdouble3 value=TDouble3(0));
  static unsigned* NewArraySeqUint(unsigned count,unsigned start=0,unsigned step=1);

  static void ToFloat3xyz(unsigned count,const tfloat4 *data,tfloat3 *dest);
  static void ToFloat1w  (unsigned count,const tfloat4 *data,float   *dest);
  static tfloat3*  NewArrayFloat3xyz(unsigned count,const tfloat4 *data);
  static float*    NewArrayFloat1w  (unsigned count,const tfloat4 *data);

  unsigned FilterApply(unsigned count,const byte *filter);
  unsigned FilterList(unsigned n,const unsigned *list);
  unsigned FilterSortList(unsigned n,const unsigned *list);

};

#endif


