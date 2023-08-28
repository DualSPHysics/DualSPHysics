//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2023 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JArraysGpu.h \brief Declares the class \ref JArraysGpu.

#ifndef _JArraysGpu_
#define _JArraysGpu_

//:#############################################################################
//:# Cambios:
//:# =========
//:# - New version 3 of JArraysGpu for single- and multi-GPU. (18-08-2023)
//:#   - Se usan objetos JArrayGpuXXX para gestionar la memoria GPU de JArraysGPU.
//:#   - Se mejora la reserva y liberacion de arrays de memoria GPU.
//:#   - Se permite aumentar el numero de arrays de memoria GPU automaticamente 
//:#     segun SetAutoAddArrays().
//:#   - Se permite redimensionar la memoria asignada a los arrays sin perdida 
//:#     de datos y los objetos JArrayGpuXXX se actualizan automaticamente.
//:#   - Los objetos JArrayGpuXXX liberan su reserva de memoria de forma 
//:#     automatica cuando se destruyen.
//:#   - Los objetos JArrayGpuXXX permiten la descarga de datos en CPU de forma
//:#     simple.
//:#############################################################################

#include "JObjectGpu.h"
#include "TypesDef.h"
#include "JArraysCpu.h"

class JLog2;
class JArraysGpu;
class JArrayXpu;

//#define DG_ARRXPU_PRINT  //-Active prints for debug.


#define AG_PTR(a) (a!=NULL? a->ptr(): NULL)
#define AG_CPTR(a) (a!=NULL? a->cptr(): NULL)

//##############################################################################
//# JArraysGpuList
//##############################################################################
class JArraysGpuList : protected JObjectGpu
{
public:
  static const unsigned MAX_ARRAYS=15;

protected:
  JArraysGpu* Head;

  bool AutoAddArrays;  ///<It allows the automatic creation of arrays in calls Reserve().

  unsigned   Count;                   ///<Number of pointers with allocated memory.  
  void*      Pointers  [MAX_ARRAYS];  ///<Pointers to GPU memory. 
  JArrayXpu* PointersAr[MAX_ARRAYS];  ///<Pointer to JArraysGpu using the pointer Pointers[]. 

  unsigned PointersUnused[MAX_ARRAYS]; ///<List of unused pointers (index to Pointers[]).
  unsigned CountUnused;                ///<Number of unused pointers.  
  unsigned CountUsedMax;

  unsigned ArraySize; ///<Number of values allocated by array.

  //-Variables for saving GPU data during resizing memory.
  unsigned SizeDataCpu;   ///<Allocated elements for GPU data.
  unsigned CountDataCpu;  ///<Number of saved data arrays.
  unsigned   SvCount;                  ///<Saves Count value. 
  void*      SvDataCpu[MAX_ARRAYS];    ///<Allocated CPU memory to save GPU data for reszing.  
  JArrayXpu* SvPointersAr[MAX_ARRAYS]; ///<Saves PointersAr in use. 

  void ResetDataCpu();
  void FreeMemory();

public:
  const unsigned ValueSize;  ///<Size of one value according to data type. 

public:
  JArraysGpuList(JArraysGpu* head,unsigned valuesize);
  ~JArraysGpuList();
  void Reset();

  llong GetAllocMemoryGpu()const{ return(llong(ValueSize)*ArraySize*Count); };

  void* AllocCpuMemory(unsigned size)const;
  void  FreeCpuMemory(void* pointer)const;

  void SetAutoAddArrays(bool active){ AutoAddArrays=active; }
  bool GetAutoAddArrays()const{ return(AutoAddArrays); }
  
  unsigned GetArrayCount()const{ return(Count); }
  unsigned GetArrayCountUsed()const{ return(Count-CountUnused); }
  unsigned GetArrayCountUsedMax()const{ return(CountUsedMax); }

  unsigned GetArraySize()const{ return(ArraySize); }

  void SetArrayCount(unsigned count);
  void SetArraySize(unsigned size);

  void Reserve(JArrayXpu* ar);
  void Free(JArrayXpu* ar);
  void SwapPtr(JArrayXpu* ar1,JArrayXpu* ar2);

  void SaveDataArray(unsigned savesizedata);
  void RestoreDataArray(unsigned newsize);
};


//##############################################################################
//# JArraysGpu
//##############################################################################
class JArraysGpu : protected JObjectGpu
{
  friend class JArrayXpu;

public:
  static const unsigned MAX_ASIZE=9;
  typedef enum{ 
     ASIZE_1B=0
    ,ASIZE_2B=1
    ,ASIZE_4B=2
    ,ASIZE_8B=3
    ,ASIZE_12B=4
    ,ASIZE_16B=5
    ,ASIZE_24B=6
    ,ASIZE_32B=7
    ,ASIZE_72B=8  // +1 -> MAX_ASIZE
  }TpASizeId; //-Arrays sizes for ArraysList.

  static unsigned GetTypeSize(TpASizeId asize){
    return(asize==ASIZE_1B?   1: (asize==ASIZE_2B?   2: 
          (asize==ASIZE_4B?   4: (asize==ASIZE_8B?   8: 
          (asize==ASIZE_12B? 12: (asize==ASIZE_16B? 16: 
          (asize==ASIZE_24B? 24: (asize==ASIZE_32B? 32: 
          (asize==ASIZE_72B? 72: 0)))))))));
  }

protected:
  JArraysGpuList* ArraysList[MAX_ASIZE];
  unsigned CountResizeDataCalls;  ///<Number of calls to SetDataArraySize().
  unsigned LastCountArrays;       ///<Number of arrays updated by GetArrayCountUpdated().

protected:
  JArraysGpuList* GetArraysList(TpASizeId asize)const{ 
    return(ArraysList[asize]);
  }
  
public:
  JLog2 *Log;
  const int GpuId;
  const int Id;

public:
  JArraysGpu(int gpuid,JLog2* log,int id=-1);
  ~JArraysGpu();
  void Reset();
  llong GetAllocMemoryGpu()const;
  unsigned GetCountResizeDataCalls()const{ return(CountResizeDataCalls); }
  void PrintAllocMemory(bool all)const;

  void SetAutoAddArrays(bool active);
  bool GetAutoAddArrays()const{ return(ArraysList[0]->GetAutoAddArrays()); }

  unsigned GetArrayCount(TpASizeId asize)const{     return(ArraysList[asize]->GetArrayCount()); }
  unsigned GetArrayCountUsed(TpASizeId asize)const{ return(ArraysList[asize]->GetArrayCountUsed()); }

  unsigned GetArrayCount()const;
  unsigned GetArrayCountUpdated();

  void AddArrayCount(TpASizeId asize,unsigned count);

  unsigned GetArraySize()const{ return(ArraysList[0]->GetArraySize()); }

  void SetArraySize(unsigned size);
  void SetDataArraySize(unsigned newsize,unsigned savesizedata);
};


//##############################################################################
//# JArraysGpu
//##############################################################################
/// \brief Defines the type of elements of the arrays managed in \ref JArraysGpu.
class JArrayXpu : protected JObjectGpu
{
  friend class JArraysGpuList;

public:
  const std::string Name;
  const TpTypeData Type;

protected:
  JArraysGpu*     Head;
  JArraysGpuList* ArraysList;

  unsigned PtrId;
  void* Ptr;
  void* DataCpu;

protected:
  void AssignReserve(void* ptr,unsigned ptrid);
  void ClearReserve();
  void PSwapPtr(JArrayXpu* ar);
  void PMemsetOffset(void* ptr_offset,unsigned offset,byte value,size_t size);

  void PCopyFrom(const JArrayXpu* src,size_t size);
  void PCopyFromOffset(void* dst_ptr,unsigned dst_offset,const JArrayXpu* src
    ,const void* src_ptr,unsigned src_offset,size_t size);
  void PCopyFromPointer(const void* src_ptr,size_t size);
  void PCopyFromPointerOffset(void* dst_ptr,unsigned dst_offset
    ,const void* src_ptr,unsigned src_offset,size_t size);

  void PCopyTo(void* dst_ptr,size_t size)const;
  void PCopyToOffset(const void* src_ptr,unsigned src_offset
    ,void* dst_ptr,size_t size)const;

  void PCopyToHost(JArraySpu* dst,size_t size)const;
  void PCopyToHostPointer(void* dst_ptr,size_t size)const;
  void PCopyToHostOffset(const void* src_ptr,unsigned src_offset
    ,JArraySpu* dst,void* dst_ptr,unsigned dst_offset,size_t size)const;
  void PCopyToHostPointerOffset(const void* src_ptr,unsigned src_offset
    ,void* dst_ptr,unsigned dst_offset,size_t size)const;

  void PCopyFromHost(const JArraySpu* src,size_t size);
  void PCopyFromHostPointer(const void* src_ptr,size_t size);
  void PCopyFromHostOffset(void* dst_ptr,unsigned dst_offset
    ,const JArraySpu* src,const void* src_ptr,unsigned src_offset,size_t size);
  void PCopyFromHostPointerOffset(void* dst_ptr,unsigned dst_offset
    ,const void* src_ptr,unsigned scr_offset,size_t size);


public:
  JArrayXpu(const std::string& name,TpTypeData type
    ,JArraysGpu::TpASizeId asize,JArraysGpu* head,bool reserve);
  ~JArrayXpu();
  std::string ObjectId()const;

  unsigned GetSize()const{ return(ArraysList->GetArraySize()); }
  unsigned GetValueSize()const{ return(ArraysList->ValueSize); }
  size_t Sizeof()const{ return(size_t(ArraysList->ValueSize)); }
  size_t Sizeof(size_t size)const{ return(Sizeof()*size); }

  bool Active()const{ return(Ptr!=NULL); }
  void Reserve();
  void Free();
  
  void CuMemset(byte value,size_t size);

  //-For data on CPU.
  void DataAlloc();
  void DataFree();
  void DataDown(){ DataDown(GetSize()); }
  void DataDown(unsigned ndata);
};


//##############################################################################
//# JArraysGpu
//##############################################################################
template <class T,class T2,TpTypeData ttype,JArraysGpu::TpASizeId tasize,class TCPU> 
  class JArrayXpuT : public JArrayXpu
{
public:
  JArrayXpuT(const std::string& name,JArraysGpu* head,bool reserve=false)
    :JArrayXpu(name,ttype,tasize,head,reserve){};
  //----------------------------------------------------------------------------
  void SwapPtr(JArrayXpuT* ar){ PSwapPtr((JArrayXpu*)ar); }
  //----------------------------------------------------------------------------
  T* ptr(){ return((T*)Ptr); }
  //----------------------------------------------------------------------------
  const T* cptr()const{ return((const T*)Ptr); }
  //----------------------------------------------------------------------------
  T* dat(){ return((T*)DataCpu); }
  //----------------------------------------------------------------------------
  const T* cdat()const{ return((const T*)DataCpu); }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  void CuMemsetOffset(unsigned offset,byte value,size_t size){
    PMemsetOffset(
      (void*)(ptr()? ptr()+offset: NULL)
      ,offset,value,size);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  void CuCopyFrom(const JArrayXpuT* src,size_t size){ 
    PCopyFrom((const JArrayXpu*)src,size); 
  }
  //----------------------------------------------------------------------------
  void CuCopyFrom(const T* src_ptr,size_t size){
    PCopyFromPointer((const void*)src_ptr,size); 
  }
  //----------------------------------------------------------------------------
  void CuCopyFrom2(const T2* src_ptr,size_t size){
    CuCopyFrom((const T*)src_ptr,size); 
  }
  //----------------------------------------------------------------------------
  void CuCopyFromOffset(unsigned dst_offset,const JArrayXpuT* src
    ,unsigned src_offset,size_t size)
  { 
    PCopyFromOffset(
      (void*)(ptr()? ptr()+dst_offset: NULL)
      ,dst_offset
      ,(const JArrayXpu*)src
      ,(const void*)(src && src->cptr()? src->cptr()+src_offset: NULL)
      ,src_offset,size); 
  }
  //----------------------------------------------------------------------------
  void CuCopyFromOffset(unsigned dst_offset,const T* src_ptr
    ,unsigned src_offset,size_t size)
  { 
    PCopyFromPointerOffset(
      (void*)(ptr()? ptr()+dst_offset: NULL)
      ,dst_offset
      ,(const void*)(src_ptr? src_ptr+src_offset: NULL)
      ,src_offset,size);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  void CuCopyTo(T* dst_ptr,size_t size)const{ 
    PCopyTo((void*)dst_ptr,size); 
  }
  //----------------------------------------------------------------------------
  void CuCopyToOffset(unsigned src_offset,T* dst_ptr
    ,unsigned dst_offset,size_t size)const
  { 
    PCopyToOffset(
      (const void*)(cptr()? cptr()+src_offset: NULL)
      ,src_offset
      ,(void*)(dst_ptr? dst_ptr+dst_offset: NULL)
      ,size);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  void CuCopyToHost(TCPU* dst,size_t size)const{ 
    PCopyToHost((JArraySpu*)dst,size);
  }
  //----------------------------------------------------------------------------
  void CuCopyToHost(T* dst_ptr,size_t size)const{ 
    PCopyToHostPointer((void*)dst_ptr,size);
  }
  //----------------------------------------------------------------------------
  void CuCopyToHostOffset(unsigned src_offset,TCPU* dst,unsigned dst_offset
    ,size_t size)const
  { 
    PCopyToHostOffset(
      (const void*)(cptr()? cptr()+src_offset: NULL)
      ,src_offset
      ,(JArraySpu*)dst
      ,(void*)(dst && dst->cptr()? dst->cptr()+dst_offset: NULL)
      ,dst_offset
      ,size);
  }
  //----------------------------------------------------------------------------
  void CuCopyToHostOffset(unsigned src_offset,T* dst_ptr,unsigned dst_offset
    ,size_t size)const
  { 
    PCopyToHostPointerOffset(
      (const void*)(cptr()? cptr()+src_offset: NULL)
      ,src_offset
      ,(void*)(dst_ptr? dst_ptr+dst_offset: NULL)
      ,dst_offset
      ,size);
  }
  //----------------------------------------------------------------------------
  //----------------------------------------------------------------------------
  void CuCopyFromHost(const TCPU* src,size_t size){ 
    PCopyFromHost((const JArraySpu*)src,size);
  }
  //----------------------------------------------------------------------------
  void CuCopyFromHost(const T* src_ptr,size_t size){ 
    PCopyFromHostPointer((const void*)src_ptr,size);
  }
  //----------------------------------------------------------------------------
  void CuCopyFromHost2(const T2* src_ptr,size_t size){ 
    CuCopyFromHost((const T*)src_ptr,size);
  }
  //----------------------------------------------------------------------------
  void CuCopyFromHostOffset(unsigned dst_offset,const TCPU* src,unsigned src_offset
    ,size_t size)
  { 
    PCopyFromHostOffset(
      (void*)(cptr()? cptr()+dst_offset: NULL)
      ,dst_offset
      ,(const JArraySpu*)src
      ,(const void*)(src && src->cptr()? src->cptr()+src_offset: NULL)
      ,src_offset
      ,size);
  }
  //----------------------------------------------------------------------------
  void CuCopyFromHostOffset(unsigned dst_offset,const T* src_ptr,unsigned src_offset
    ,size_t size)
  { 
    PCopyFromHostPointerOffset(
      (void*)(cptr()? cptr()+dst_offset: NULL)
      ,dst_offset
      ,(const void*)(src_ptr? src_ptr+src_offset: NULL)
      ,src_offset
      ,size);
  }

};

//##############################################################################
typedef JArrayXpuT<byte    ,byte    ,TypeUchar  ,JArraysGpu::ASIZE_1B ,JArraySpuByte>    JArrayXpuByte;
typedef JArrayXpuT<word    ,word    ,TypeUshort ,JArraysGpu::ASIZE_2B ,JArraySpuWord>    JArrayXpuWord;
typedef JArrayXpuT<ushort4 ,tword4  ,TypeNull   ,JArraysGpu::ASIZE_8B ,JArraySpuWord4>   JArrayXpuWord4;  // ***Use TypeNull***
typedef JArrayXpuT<int     ,int     ,TypeInt    ,JArraysGpu::ASIZE_4B ,JArraySpuInt>     JArrayXpuInt;
typedef JArrayXpuT<unsigned,unsigned,TypeUint   ,JArraysGpu::ASIZE_4B ,JArraySpuUint>    JArrayXpuUint;
typedef JArrayXpuT<float   ,float   ,TypeFloat  ,JArraysGpu::ASIZE_4B ,JArraySpuFloat>   JArrayXpuFloat;
typedef JArrayXpuT<float2  ,tfloat2 ,TypeFloat2 ,JArraysGpu::ASIZE_8B ,JArraySpuFloat2>  JArrayXpuFloat2;
typedef JArrayXpuT<float3  ,tfloat3 ,TypeFloat3 ,JArraysGpu::ASIZE_12B,JArraySpuFloat3>  JArrayXpuFloat3;
typedef JArrayXpuT<float4  ,tfloat4 ,TypeFloat4 ,JArraysGpu::ASIZE_16B,JArraySpuFloat4>  JArrayXpuFloat4;
typedef JArrayXpuT<double  ,double  ,TypeDouble ,JArraysGpu::ASIZE_8B ,JArraySpuDouble>  JArrayXpuDouble;
typedef JArrayXpuT<double2 ,double2 ,TypeDouble2,JArraysGpu::ASIZE_16B,JArraySpuDouble2> JArrayXpuDouble2;
typedef JArrayXpuT<double3 ,double3 ,TypeDouble3,JArraysGpu::ASIZE_24B,JArraySpuDouble3> JArrayXpuDouble3;
typedef JArrayXpuT<double4 ,double4 ,TypeDouble4,JArraysGpu::ASIZE_32B,JArraySpuDouble4> JArrayXpuDouble4;

typedef JArrayXpuT<tmatrix2d  ,tmatrix2d  ,TypeNull,JArraysGpu::ASIZE_32B,JArraySpuMatrix2d> JArrayXpuMatrix2d;       // ***Use TypeNull***
typedef JArrayXpuT<tmatrix3d  ,tmatrix3d  ,TypeNull,JArraysGpu::ASIZE_72B,JArraySpuMatrix3d> JArrayXpuMatrix3d;       // ***Use TypeNull***
typedef JArrayXpuT<tsymatrix3f,tsymatrix3f,TypeNull,JArraysGpu::ASIZE_24B,JArraySpuSymMatrix3f> JArrayXpuSymMatrix3f; // ***Use TypeNull***

//##############################################################################
typedef JArrayXpuByte     agbyte;
typedef JArrayXpuWord     agword;
typedef JArrayXpuWord4    agword4;
typedef JArrayXpuInt      agint;
typedef JArrayXpuUint     aguint;
typedef JArrayXpuFloat    agfloat;
typedef JArrayXpuFloat2   agfloat2;
typedef JArrayXpuFloat3   agfloat3;
typedef JArrayXpuFloat4   agfloat4;
typedef JArrayXpuDouble   agdouble;
typedef JArrayXpuDouble2  agdouble2;
typedef JArrayXpuDouble3  agdouble3;
typedef JArrayXpuDouble4  agdouble4;
typedef JArrayXpuMatrix2d agmatrix2d;
typedef JArrayXpuMatrix3d agmatrix3d;
typedef JArrayXpuSymMatrix3f agsymatrix3f;

//-Code type (typecode) size 2-bytes or 4-bytes is defined in DualSphDef.h
#ifdef CODE_SIZE4  
  typedef JArrayXpuUint JArrayXpuTypecode;
  typedef JArrayXpuUint agtypecode;
#else
  typedef JArrayXpuWord JArrayXpuTypecode;
  typedef JArrayXpuWord agtypecode;
#endif

#endif

