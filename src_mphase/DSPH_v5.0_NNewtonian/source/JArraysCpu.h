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
//:# - Codigo creado a partir de JArraysGpu para usar con memoria CPU. (10-03-2014)
//:# - Remplaza long long por llong. (01-10-2015)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:#############################################################################

/// \file JArraysCpu.h \brief Declares the class \ref JArraysCpu.

#ifndef _JArraysCpu_
#define _JArraysCpu_

#include "JObject.h"
#include "DualSphDef.h"

//##############################################################################
//# JArraysCpuSize
//##############################################################################
/// \brief Defines the type of elements of the arrays managed in \ref JArraysCpu with a given size.

class JArraysCpuSize : protected JObject
{
protected:
  const unsigned ElementSize;
  unsigned ArraySize;

  static const unsigned MAXPOINTERS=30;
  void* Pointers[MAXPOINTERS];
  unsigned Count;
  unsigned CountUsed;

  unsigned CountMax,CountUsedMax;
  
  void* AllocPointer(unsigned size)const;
  void FreePointer(void* pointer)const;

  void FreeMemory();
  unsigned FindPointerUsed(void *pointer)const;

public:
  JArraysCpuSize(unsigned elementsize);
  ~JArraysCpuSize();
  void Reset();
  
  void SetArrayCount(unsigned count);
  unsigned GetArrayCount()const{ return(Count); }
  unsigned GetArrayCountUsed()const{ return(CountUsed); }
  
  unsigned GetArrayCountMax()const{ return(CountMax); }
  unsigned GetArrayCountUsedMax()const{ return(CountUsedMax); }

  void SetArraySize(unsigned size);
  unsigned GetArraySize()const{ return(ArraySize); }

  llong GetAllocMemoryCpu()const{ return((llong)(Count)*ElementSize*ArraySize); };

  void* Reserve();
  void Free(void *pointer);
};


//##############################################################################
//# JArraysCpu
//##############################################################################
/// \brief Defines the type of elements of the arrays managed in \ref JArraysCpu. 

class JArraysCpu : protected JObject
{
public:
  typedef enum{ SIZE_1B=1,SIZE_2B=2,SIZE_4B=4,SIZE_8B=8,SIZE_12B=12,SIZE_16B=16,SIZE_24B=24,SIZE_32B=32 }TpArraySize;  //-Tipos de arrays.

protected:
  JArraysCpuSize *Arrays1b;
  JArraysCpuSize *Arrays2b;
  JArraysCpuSize *Arrays4b;
  JArraysCpuSize *Arrays8b;
  JArraysCpuSize *Arrays12b;
  JArraysCpuSize *Arrays16b;
  JArraysCpuSize *Arrays24b;
  JArraysCpuSize *Arrays32b;
  
  JArraysCpuSize* GetArrays(TpArraySize tsize)const{ return(tsize==SIZE_32B? Arrays32b: (tsize==SIZE_24B? Arrays24b: (tsize==SIZE_16B? Arrays16b: (tsize==SIZE_12B? Arrays12b: (tsize==SIZE_8B? Arrays8b: (tsize==SIZE_4B? Arrays4b: (tsize==SIZE_2B? Arrays2b: Arrays1b))))))); }

public:
  JArraysCpu();
  ~JArraysCpu();
  void Reset();
  llong GetAllocMemoryCpu()const;
  
  void SetArrayCount(TpArraySize tsize,unsigned count){ GetArrays(tsize)->SetArrayCount(count); }
  void AddArrayCount(TpArraySize tsize,unsigned count=1){ SetArrayCount(tsize,GetArrayCount(tsize)+count); }
  unsigned GetArrayCount(TpArraySize tsize)const{ return(GetArrays(tsize)->GetArrayCount()); }
  unsigned GetArrayCountUsed(TpArraySize tsize)const{ return(GetArrays(tsize)->GetArrayCountUsed()); }

  void SetArraySize(unsigned size);
  unsigned GetArraySize()const{ return(Arrays1b->GetArraySize()); }

  byte*        ReserveByte(){       return((byte*)Arrays1b->Reserve());         }
  word*        ReserveWord(){       return((word*)Arrays2b->Reserve());         }
  unsigned*    ReserveUint(){       return((unsigned*)Arrays4b->Reserve());     }
  int*         ReserveInt(){        return((int*)Arrays4b->Reserve());          }
  float*       ReserveFloat(){      return((float*)Arrays4b->Reserve());        }
  tfloat3*     ReserveFloat3(){     return((tfloat3*)Arrays12b->Reserve());     }
  tfloat4*     ReserveFloat4(){     return((tfloat4*)Arrays16b->Reserve());     }
  double*      ReserveDouble(){     return((double*)Arrays8b->Reserve());       }
  tdouble2*    ReserveDouble2(){    return((tdouble2*)Arrays16b->Reserve());    }
  tdouble3*    ReserveDouble3(){    return((tdouble3*)Arrays24b->Reserve());    }
  tsymatrix3f* ReserveSymatrix3f(){ return((tsymatrix3f*)Arrays24b->Reserve()); }
#ifdef CODE_SIZE4
  typecode*    ReserveTypeCode(){   return(ReserveUint());                      }
#else
  typecode*    ReserveTypeCode(){   return(ReserveWord());                      }
#endif

  void Free(byte        *pointer){ Arrays1b->Free(pointer);  }
  void Free(word        *pointer){ Arrays2b->Free(pointer);  }
  void Free(unsigned    *pointer){ Arrays4b->Free(pointer);  }
  void Free(int         *pointer){ Arrays4b->Free(pointer);  }
  void Free(float       *pointer){ Arrays4b->Free(pointer);  }
  void Free(tfloat3     *pointer){ Arrays12b->Free(pointer); }
  void Free(tfloat4     *pointer){ Arrays16b->Free(pointer); }
  void Free(double      *pointer){ Arrays8b->Free(pointer);  }
  void Free(tdouble2    *pointer){ Arrays16b->Free(pointer); }
  void Free(tdouble3    *pointer){ Arrays24b->Free(pointer); }
  void Free(tsymatrix3f *pointer){ Arrays24b->Free(pointer); }
};


#endif


