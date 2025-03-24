//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JArraysCpu.h \brief Declares the class \ref JArraysCpu.

#ifndef _JArraysCpu_
#define _JArraysCpu_

//:#############################################################################
//:# Cambios:
//:# =========
//:# - New version 3 of JArraysCpu equivalent to a JArraysGpu. (18-08-2023)
//:#   - Se usan objetos JArrayCpuXXX para gestionar la memoria CPU de JArraysCPU.
//:#   - Se mejora la reserva y liberacion de arrays de memoria CPU.
//:#   - Se permite aumentar el numero de arrays de memoria CPU automaticamente 
//:#     segun SetAutoAddArrays().
//:#   - Se permite redimensionar la memoria asignada a los arrays sin perdida
//:#     de datos y los objetos JArrayCpuXXX se actualizan automaticamente.
//:#   - Los objetos JArrayCpuXXX liberan su reserva de memoria de forma
//:#     automatica cuando se destruyen.
//:#############################################################################

#include "JObject.h"
#include "TypesDef.h"

class JLog2;
class JArraysCpu;
class JArrayCpu;

//#define DG_ARRSPU_PRINT  //-Active prints for debug.
  

#define AC_PTR(a)     (a!=NULL? a->ptr()   : NULL)
#define AC_PTRV(a,v)  (a!=NULL? a->ptr()+v : NULL)
#define AC_CPTR(a)    (a!=NULL? a->cptr()  : NULL)
#define AC_CPTRV(a,v) (a!=NULL? a->cptr()+v: NULL)

//##############################################################################
//# JArraysCpuList
//##############################################################################
class JArraysCpuList : protected JObject
{
public:
  static const unsigned MAX_ARRAYS=15;

protected:
  JArraysCpu* Head;

  bool AutoAddArrays;  ///<It allows the automatic creation of arrays in calls Reserve().

  unsigned   Count;                   ///<Number of pointers with allocated memory.  
  void*      Pointers  [MAX_ARRAYS];  ///<Pointers to CPU memory. 
  JArrayCpu* PointersAr[MAX_ARRAYS];  ///<Pointer to JArraysCpu using the pointer Pointers[]. 

  unsigned PointersUnused[MAX_ARRAYS]; ///<List of unused pointers (index to Pointers[]).
  unsigned CountUnused;                ///<Number of unused pointers.  
  unsigned CountUsedMax;

  unsigned ArraySize; ///<Number of values allocated by array.

  void FreeMemory();

public:
  const unsigned ValueSize;  ///<Size of one value according to data type. 

public:
  JArraysCpuList(JArraysCpu* head,unsigned valuesize);
  ~JArraysCpuList();
  void Reset();
  void PrintArraysInfo()const;

  llong GetAllocMemoryCpu()const{ return(llong(ValueSize)*ArraySize*Count); };

  void* AllocCpuMemory(unsigned size)const;
  void  FreeCpuMemory(void* pointer)const;

  void SetAutoAddArrays(bool active){ AutoAddArrays=active; }
  bool GetAutoAddArrays()const{ return(AutoAddArrays); }
  
  unsigned GetArrayCount()const{ return(Count); }
  unsigned GetArrayCountUsed()const{ return(Count-CountUnused); }
  unsigned GetArrayCountUsedMax()const{ return(CountUsedMax); }

  void SetArrayCount(unsigned count);
  void FreeUnusedArrays();

  unsigned GetArraySize()const{ return(ArraySize); }

  void SetArraySize(unsigned size);

  void Reserve(JArrayCpu* ar);
  void Free(JArrayCpu* ar);
  void SwapPtr(JArrayCpu* ar1,JArrayCpu* ar2);

  void SetDataArraySize(unsigned newsize,unsigned savesizedata);
};


//##############################################################################
//# JArraysCpu
//##############################################################################
class JArraysCpu : protected JObject
{
  friend class JArrayCpu;

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
  JArraysCpuList* ArraysList[MAX_ASIZE];
  unsigned CountResizeDataCalls;  ///<Number of calls to SetDataArraySize().
  unsigned LastCountArrays;       ///<Number of arrays updated by GetArrayCountUpdated().

protected:
  JArraysCpuList* GetArraysList(TpASizeId asize)const{ 
    return(ArraysList[asize]);
  }
  
public:
  JLog2* Log;
  const int Id;

public:
  JArraysCpu(JLog2* log,int id=-1);
  ~JArraysCpu();
  void Reset();
  llong GetAllocMemoryCpu()const;
  unsigned GetCountResizeDataCalls()const{ return(CountResizeDataCalls); }
  void PrintAllocMemory(bool emptyblocks,bool arrayinfo)const;

  void SetAutoAddArrays(bool active);
  bool GetAutoAddArrays()const{ return(ArraysList[0]->GetAutoAddArrays()); }

  unsigned GetArrayCount(TpASizeId asize)const{     return(ArraysList[asize]->GetArrayCount()); }
  unsigned GetArrayCountUsed(TpASizeId asize)const{ return(ArraysList[asize]->GetArrayCountUsed()); }

  unsigned GetArrayCount()const;
  unsigned GetArrayCountUpdated();

  void AddArrayCount(TpASizeId asize,unsigned count);
  void FreeUnusedArrays();

  unsigned GetArraySize()const{ return(ArraysList[0]->GetArraySize()); }

  void SetArraySize(unsigned size);
  void SetDataArraySize(unsigned newsize,unsigned savesizedata);
};


//##############################################################################
//# JArraysCpu
//##############################################################################
/// \brief Defines the type of elements of the arrays managed in \ref JArraysCpu.
class JArrayCpu : protected JObject
{
  friend class JArraysCpuList;

public:
  const std::string Name;
  const TpTypeData Type;

protected:
  JArraysCpu*     Head;
  JArraysCpuList* ArraysList;

  unsigned PtrId;
  void* Ptr;
  bool Locked;

protected:
  void AssignReserve(void* ptr,unsigned ptrid);
  void ClearReserve();
  void PSwapPtr(JArrayCpu* ar);
  void PMemsetOffset(void* ptr_offset,unsigned offset,byte value,size_t size);
  void PCopyFrom(const JArrayCpu* src,size_t size);
  void PCopyFromOffset(void* dst_ptr,unsigned dst_offset,const JArrayCpu* src
    ,const void* src_ptr,unsigned src_offset,size_t size);
  void PCopyFromPointer(const void* src_ptr,size_t size);
  void PCopyFromPointerOffset(void* dst_ptr,unsigned dst_offset
    ,const void* src_ptr,unsigned src_offset,size_t size);
  void PCopyTo(void* dst_ptr,size_t size)const;
  void PCopyToOffset(const void* src_ptr,unsigned src_offset
    ,void* dst_ptr,size_t size)const;

public:
  JArrayCpu(const std::string& name,TpTypeData type
    ,JArraysCpu::TpASizeId asize,JArraysCpu* head,bool reserve);
  ~JArrayCpu();
  std::string ObjectId()const;

  void* ptrvoid(){ return(Ptr); }
  const void* cptrvoid()const{ return(Ptr); }

  unsigned GetSize()const{ return(ArraysList->GetArraySize()); }
  unsigned GetValueSize()const{ return(ArraysList->ValueSize); }
  size_t Sizeof()const{ return(size_t(ArraysList->ValueSize)); }
  size_t Sizeof(size_t size)const{ return(Sizeof()*size); }

  bool Active()const{ return(Ptr!=NULL); }
  void Reserve();
  void Free();

  bool IsLocked()const{ return(Locked); };
  void LockPtr();
  void UnlockPtr();

  void Memset(byte value,size_t size);
};


//##############################################################################
//# JArraysCpu
//##############################################################################
template <class T,TpTypeData ttype,JArraysCpu::TpASizeId tasize> class JArrayCpuT : public JArrayCpu
{
public:
  JArrayCpuT(const std::string& name,JArraysCpu* head,bool reserve=false)
    :JArrayCpu(name,ttype,tasize,head,reserve){};
  //----------------------------------------------------------------------------
  void SwapPtr(JArrayCpuT* ar){ PSwapPtr((JArrayCpu*)ar); }
  //----------------------------------------------------------------------------
  T* ptr(){ return((T*)Ptr); }
  //----------------------------------------------------------------------------
  const T* cptr()const{ return((const T*)Ptr); }
  //----------------------------------------------------------------------------
  void MemsetOffset(unsigned offset,byte value,size_t size){
    PMemsetOffset(
      (void*)(ptr()? ptr()+offset: NULL)
      ,offset,value,size);
  }
  //----------------------------------------------------------------------------
  void CopyFrom(const JArrayCpuT* src,size_t size){ 
    PCopyFrom((const JArrayCpu*)src,size); 
  }
  //----------------------------------------------------------------------------
  void CopyFrom(const T* src_ptr,size_t size){
    PCopyFromPointer((const void*)src_ptr,size); 
  }
  //----------------------------------------------------------------------------
  void CopyFromOffset(unsigned dst_offset,const JArrayCpuT* src
    ,unsigned src_offset,size_t size)
  { 
    PCopyFromOffset(
      (void*)(ptr()? ptr()+dst_offset: NULL)
      ,dst_offset
      ,(const JArrayCpu*)src
      ,(const void*)(src && src->cptr()? src->cptr()+src_offset: NULL)
      ,src_offset,size); 
  }
  //----------------------------------------------------------------------------
  void CopyFromOffset(unsigned dst_offset,const T* src_ptr
    ,unsigned src_offset,size_t size)
  { 
    PCopyFromPointerOffset(
      (void*)(ptr()? ptr()+dst_offset: NULL)
      ,dst_offset
      ,(const void*)(src_ptr? src_ptr+src_offset: NULL)
      ,src_offset,size);
  }
  //----------------------------------------------------------------------------
  void CopyTo(T* dst_ptr,size_t size)const{ 
    PCopyTo((void*)dst_ptr,size); 
  }
  //----------------------------------------------------------------------------
  void CopyToOffset(unsigned src_offset,T* dst_ptr
    ,unsigned dst_offset,size_t size)const
  { 
    PCopyToOffset(
      (const void*)(cptr()? cptr()+src_offset: NULL)
      ,src_offset
      ,(void*)(dst_ptr? dst_ptr+dst_offset: NULL)
      ,size);
  }
};

//##############################################################################
typedef JArrayCpuT<byte    ,TypeUchar  ,JArraysCpu::ASIZE_1B>  JArrayCpuByte;
typedef JArrayCpuT<word    ,TypeUshort ,JArraysCpu::ASIZE_2B>  JArrayCpuWord;
typedef JArrayCpuT<tword4  ,TypeNull   ,JArraysCpu::ASIZE_8B>  JArrayCpuWord4;
typedef JArrayCpuT<int     ,TypeInt    ,JArraysCpu::ASIZE_4B>  JArrayCpuInt;
typedef JArrayCpuT<unsigned,TypeUint   ,JArraysCpu::ASIZE_4B>  JArrayCpuUint;
typedef JArrayCpuT<float   ,TypeFloat  ,JArraysCpu::ASIZE_4B>  JArrayCpuFloat;
typedef JArrayCpuT<tfloat2 ,TypeFloat2 ,JArraysCpu::ASIZE_8B>  JArrayCpuFloat2;
typedef JArrayCpuT<tfloat3 ,TypeFloat3 ,JArraysCpu::ASIZE_12B> JArrayCpuFloat3;
typedef JArrayCpuT<tfloat4 ,TypeFloat4 ,JArraysCpu::ASIZE_16B> JArrayCpuFloat4;
typedef JArrayCpuT<double  ,TypeDouble ,JArraysCpu::ASIZE_8B>  JArrayCpuDouble;
typedef JArrayCpuT<tdouble2,TypeDouble2,JArraysCpu::ASIZE_16B> JArrayCpuDouble2;
typedef JArrayCpuT<tdouble3,TypeDouble3,JArraysCpu::ASIZE_24B> JArrayCpuDouble3;
typedef JArrayCpuT<tdouble4,TypeDouble4,JArraysCpu::ASIZE_32B> JArrayCpuDouble4;

typedef JArrayCpuT<tmatrix2d  ,TypeNull,JArraysCpu::ASIZE_32B> JArrayCpuMatrix2d;    // ***Use TypeNull***
typedef JArrayCpuT<tmatrix3d  ,TypeNull,JArraysCpu::ASIZE_72B> JArrayCpuMatrix3d;    // ***Use TypeNull***
typedef JArrayCpuT<tsymatrix3f,TypeNull,JArraysCpu::ASIZE_24B> JArrayCpuSymMatrix3f; // ***Use TypeNull***

//##############################################################################
typedef JArrayCpuByte     acbyte;
typedef JArrayCpuWord     acword;
typedef JArrayCpuWord4    acword4;
typedef JArrayCpuInt      acint;
typedef JArrayCpuUint     acuint;
typedef JArrayCpuFloat    acfloat;
typedef JArrayCpuFloat2   acfloat2;
typedef JArrayCpuFloat3   acfloat3;
typedef JArrayCpuFloat4   acfloat4;
typedef JArrayCpuDouble   acdouble;
typedef JArrayCpuDouble2  acdouble2;
typedef JArrayCpuDouble3  acdouble3;
typedef JArrayCpuDouble4  acdouble4;
typedef JArrayCpuMatrix2d acmatrix2d;
typedef JArrayCpuMatrix3d acmatrix3d;
typedef JArrayCpuSymMatrix3f acsymatrix3f;

//-Code type (typecode) size 2-bytes or 4-bytes is defined in DualSphDef.h
#ifdef CODE_SIZE4  
  typedef JArrayCpuUint JArrayCpuTypecode;
  typedef JArrayCpuUint actypecode;
#else
  typedef JArrayCpuWord JArrayCpuTypecode;
  typedef JArrayCpuWord actypecode;
#endif

#endif

