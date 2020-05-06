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

/// \file JRadixSort.cpp \brief Implements the class \ref JRadixSort.

#include "JRadixSort.h"
#include <string>
#include <cstring>
#include <climits>
#include <algorithm>

using namespace std;


//==============================================================================
/// Constructor de objetos.
/// Object constructor.
//==============================================================================
JRadixSort::JRadixSort(bool useomp):UseOmp(useomp && CompiledOMP()){
  ClassName="JRadixSort";
  InitData32=NULL; InitData64=NULL;
  Data32=NULL; Data64=NULL;
  PrevData32=NULL; PrevData64=NULL;
  BeginKeys=NULL;
  Index=NULL; PrevIndex=NULL;
  Reset();
}

//==============================================================================
/// Destructor de objetos.
/// Object destructor.
//==============================================================================
JRadixSort::~JRadixSort(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Reinicializa el estado del objeto, recuperando la configuracion por defecto.
/// Resets the object's state, restoring default settings.
//==============================================================================
void JRadixSort::Reset(){
  if(Data32==InitData32)Data32=NULL;
  if(Data64==InitData64)Data64=NULL;
  if(PrevData32==InitData32)PrevData32=NULL;
  if(PrevData64==InitData64)PrevData64=NULL;
  InitData32=NULL; InitData64=NULL;
  delete[] Data32; Data32=NULL;
  delete[] Data64; Data64=NULL;
  delete[] PrevData32; PrevData32=NULL;
  delete[] PrevData64; PrevData64=NULL;
  Size=Nbits=Nkeys=0;
  delete[] BeginKeys; BeginKeys=NULL;
  delete[] Index; Index=NULL;
  delete[] PrevIndex; PrevIndex=NULL;
}

//==============================================================================
/// Indica si el codigo fue compilado con OpenMP.
/// Returns if the code was compiled with OpenMP.
//==============================================================================
bool JRadixSort::CompiledOMP(){
  #ifdef OMP_USE_RADIXSORT
    return(true);
  #else
    return(false);
  #endif
}

//==============================================================================
/// Devuelve el numero de bits necesarios para codificar el valor indicado.
/// Returns the number of bits needed to encode the specified value.
//==============================================================================
template<class T> unsigned JRadixSort::TBitsSize(T v,unsigned smax)const{
  unsigned nbits=1;
  for(;v>>nbits&&nbits<smax;nbits++);
  return(nbits);
}

//==============================================================================
/// Devuelve el numero de bits necesarios para codificar el valor indicado.
/// Returns the number of bits needed to encode the specified value.
//==============================================================================
unsigned JRadixSort::BitsSize(unsigned v)const{ return(TBitsSize<unsigned>(v,32)); };

//==============================================================================
/// Devuelve el numero de bits necesarios para codificar el valor indicado.
/// Returns the number of bits needed to encode the specified value.
//==============================================================================
unsigned JRadixSort::BitsSize(ullong v)const{ return(TBitsSize<ullong>(v,64)); };

//==============================================================================
/// Calcula Nbits para los datos indicados.
/// Computes Nbits for the indicated data.
//==============================================================================
template<class T> unsigned JRadixSort::TCalcNbits(unsigned size,const T *data)const{
  const int threads=omp_get_max_threads();
  T mxdata=0;

  if(!UseOmp || threads<2){//-Secuencial. //-Sequential
    T vmax=0;
    for(unsigned c=0;c<size;c++)vmax=max(vmax,data[c]);
    mxdata=vmax;
  }
  else{//-Con OpenMP. //-With OpenMP.
    //-Calcula bloques de ejecucion.
    //-Computes execution blocks.
    const int nk=int(size/OMPSIZE)+1;
    if(nk<0)Run_Exceptioon("Number of values is invalid.");
    const int rk=int(size%OMPSIZE);
    //-Calcula maximo de nk bloques con varios hilos.
    //-Calculate maximum of nk blocks with several threads.
    T *vmax=new T[threads*OMPSTRIDE];
    memset(vmax,0,sizeof(T)*threads*OMPSTRIDE);
    #ifdef OMP_USE_RADIXSORT
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      int th=omp_get_thread_num();
      T mx=vmax[OMPSTRIDE*th];
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      for(unsigned c2=c2ini;c2<c2fin;c2++)mx=max(mx,data[c2]);
      vmax[OMPSTRIDE*th]=mx;
    }
    //-Calcula reduce maximo de todos los hilos.
    //-Computes maximum for all threads.
    T mx=0;
    for(int t=0;t<threads;t++)mx=max(mx,vmax[OMPSTRIDE*t]);
    delete[] vmax; vmax=NULL;
    mxdata=mx;
  }
  //-Calcula nbits para el valor maximo.
  //-Computes nbits for the maximum value.
  return(BitsSize(mxdata));
}

//==============================================================================
/// Calcula Nbits para los datos indicados.
/// Computes Nbits for the indicated data.
//==============================================================================
unsigned JRadixSort::CalcNbits(unsigned size,const unsigned *data)const{ return(TCalcNbits<unsigned>(size,data)); }

//==============================================================================
/// Calcula Nbits para los datos indicados.
/// Computes Nbits for the indicated data.
//==============================================================================
unsigned JRadixSort::CalcNbits(unsigned size,const ullong *data)const{ return(TCalcNbits<ullong>(size,data)); }

//==============================================================================
/// Reserva la memoria necesaria.
/// Allocates the necessary memory.
//==============================================================================
void JRadixSort::AllocMemory(unsigned s){
  try{
    if(Type32)Data32=new unsigned[s];
    else Data64=new ullong[s];
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Cannot allocate the requested memory.");
  }
}

//==============================================================================
/// Contabiliza numero de valores para cada clave.
/// Counts number of values for each key.
//==============================================================================
template<class T> void JRadixSort::LoadBeginKeys(const T* data){
  const int threads=omp_get_max_threads();
  //-Reserva espacio para contadores de claves.
  //-Allocates space for key counters.
  Nkeys=unsigned((Nbits+(KEYSBITS-1))/KEYSBITS);
  delete[] BeginKeys; BeginKeys=NULL;
  BeginKeys=new unsigned[Nkeys*KEYSRANGE];

  //-Inicia proceso.
  //-Initialises process.
  if(!UseOmp || threads<2){//-Secuencial. //-Sequential
    unsigned *nkeys=new unsigned[Nkeys*KEYSRANGE];
    memset(nkeys,0,sizeof(unsigned)*Nkeys*KEYSRANGE);
    for(unsigned c2=0;c2<Size;c2++){
      const T v=data[c2];
      for(unsigned ck=0;ck<Nkeys;ck++){ 
        const unsigned k=((v>>(ck*KEYSBITS))&KEYSMASK);
        nkeys[ck*KEYSRANGE+k]++;
      } 
    }
    //-Carga valores en BeginKeys.
    //-Loads values in BeginKeys.
    for(unsigned ck=0;ck<Nkeys;ck++){
      BeginKeys[ck*KEYSRANGE]=0;
      for(unsigned c=1;c<KEYSRANGE;c++){
        BeginKeys[ck*KEYSRANGE+c]=BeginKeys[ck*KEYSRANGE+c-1]+nkeys[ck*KEYSRANGE+c-1];
      }
    }
    delete[] nkeys;
  }
  else{//-con OpenMP. //-with OpenMP.
    //-Calcula bloques de ejecucion.
    //-Computes execution blocks.
    const int nk=int(Size/OMPSIZE)+1;
    if(nk<0)Run_Exceptioon("Number of values is invalid.");
    const int rk=int(Size%OMPSIZE);
    //-Reserva memoria auxiliar para conteo.
    //-Allocates auxiliary memory for counting.
    const unsigned skeys=Nkeys*KEYSRANGE+100;
    unsigned *nkeys=new unsigned[skeys*threads];
    memset(nkeys,0,sizeof(unsigned)*skeys*threads);
    //-Realiza conteo con varios hilos.
    //-Performs count with several threads.
    #ifdef OMP_USE_RADIXSORT
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      int th=omp_get_thread_num();
      unsigned *n=nkeys+(skeys*th);
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      //printf(">> c2ini:%d  c2fin:%d\n",c2ini,c2fin);
      for(unsigned c2=c2ini;c2<c2fin;c2++){
        T v=data[c2];
        for(unsigned ck=0;ck<Nkeys;ck++){ 
          unsigned k=((v>>(ck*KEYSBITS))&KEYSMASK);
          n[ck*KEYSRANGE+k]++;
        } 
      }
    }
    //-Reduce conteo de todos los hilos.
    //-Reduced count for all threads.
    for(int t=1;t<threads;t++)for(unsigned ck=0;ck<Nkeys;ck++)for(unsigned c=0;c<KEYSRANGE;c++)nkeys[ck*KEYSRANGE+c]+=nkeys[t*skeys+ck*KEYSRANGE+c];
    //-Carga valores en BeginKeys.
    //-Loads values in BeginKeys.
    for(unsigned ck=0;ck<Nkeys;ck++){
      BeginKeys[ck*KEYSRANGE]=0;
      for(unsigned c=1;c<KEYSRANGE;c++){
        BeginKeys[ck*KEYSRANGE+c]=BeginKeys[ck*KEYSRANGE+c-1]+nkeys[ck*KEYSRANGE+c-1];
      }
    }
    //-Libera memoria auxiliar.
    //-Frees auxiliary memory.
    delete[] nkeys;
  }
}

//==============================================================================
/// Realiza un paso de ordenacion en funcion de 1 bit.
/// Performs a sorting step in an 1 bit function.
//==============================================================================
template<class T> void JRadixSort::SortStep(unsigned ck,const T* data,T* data2){
  unsigned p2[KEYSRANGE];
  memcpy(p2,BeginKeys+(ck*KEYSRANGE),sizeof(unsigned)*KEYSRANGE);
  const unsigned ckmov=ck*KEYSBITS;
  for(unsigned p=0;p<Size;p++){
    unsigned k=((data[p]>>ckmov)&KEYSMASK);
    data2[p2[k]]=data[p];
    p2[k]++;
  }
}

//==============================================================================
/// Realiza un paso de ordenacion en funcion de 1 bit.
/// Performs a sorting step in an 1 bit function.
//==============================================================================
template<class T> void JRadixSort::SortStepIndex(unsigned ck,const T* data,T* data2,const unsigned *index,unsigned *index2){
  unsigned p2[KEYSRANGE];
  memcpy(p2,BeginKeys+(ck*KEYSRANGE),sizeof(unsigned)*KEYSRANGE);
  const unsigned ckmov=ck*KEYSBITS;
  for(unsigned p=0;p<Size;p++){
    unsigned k=((data[p]>>ckmov)&KEYSMASK);
    unsigned pk=p2[k];
    data2[pk]=data[p];
    index2[pk]=index[p];
    p2[k]++;
  }
}

//==============================================================================
/// Crea e inicializa el vector Index[].
/// Creates and initializes the Index[] array.
//==============================================================================
void JRadixSort::IndexCreate(){
  const int threads=omp_get_max_threads();
  //-Reserva memoria.
  //-Allocates memeory.
  try{
    Index=new unsigned[Size];
    PrevIndex=new unsigned[Size];
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Cannot allocate the requested memory.");
  }

  //-Carga PrevIndex[] con valores consecutivos.
  //-Loads PrevIndex[] with consecutive values.
  if(!UseOmp || threads<2){//-Secuencial. //-Sequential.
    for(unsigned c2=0;c2<Size;c2++)PrevIndex[c2]=c2;
  }
  else{//-con OpenMP.
    const int nk=int(Size/OMPSIZE)+1;
    if(nk<0)Run_Exceptioon("Number of values is invalid.");
    const int rk=int(Size%OMPSIZE);
    //-Realiza proceso con varios hilos.
    //-Performs process with several threads.
    #ifdef OMP_USE_RADIXSORT
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      for(unsigned c2=c2ini;c2<c2fin;c2++)PrevIndex[c2]=c2;
    }
  }
}

//==============================================================================
/// Ordena valores de data.
/// Reorders data values.
//==============================================================================
void JRadixSort::Sort(bool makeindex,unsigned size,unsigned *data,unsigned nbits){
  Reset();
  Nbits=nbits; Size=size; 
  Type32=true; InitData32=data; PrevData32=data;
  AllocMemory(Size);
  if(makeindex)IndexCreate();
  LoadBeginKeys<unsigned>(PrevData32);
  for(unsigned ck=0;ck<Nkeys;ck++){
    if(makeindex){
      SortStepIndex(ck,PrevData32,Data32,PrevIndex,Index);
      swap(PrevIndex,Index);
    }
    else SortStep(ck,PrevData32,Data32);
    swap(PrevData32,Data32);
  }
  if(makeindex){ 
    swap(PrevIndex,Index);
    delete[] PrevIndex; PrevIndex=NULL;
  }
  //-Copia los datos en el puntero recibido como parametro.
  //-Copies data in the pointer received as a parameter.
  if(PrevData32!=InitData32)memcpy(InitData32,PrevData32,sizeof(unsigned)*Size);
}

//==============================================================================
/// Ordena valores de data.
/// Reorders data values.
//==============================================================================
void JRadixSort::Sort(bool makeindex,unsigned size,ullong *data,unsigned nbits){
  Reset();
  Nbits=nbits; Size=size; 
  Type32=false; InitData64=data; PrevData64=data;
  AllocMemory(Size);
  if(makeindex)IndexCreate();
  LoadBeginKeys<ullong>(PrevData64);
  for(unsigned ck=0;ck<Nkeys;ck++){
    if(makeindex){
      SortStepIndex(ck,PrevData64,Data64,PrevIndex,Index);
      swap(PrevIndex,Index);
    }
    else SortStep(ck,PrevData64,Data64);
    swap(PrevData64,Data64);
  }
  if(makeindex){ 
    swap(PrevIndex,Index);
    delete[] PrevIndex; PrevIndex=NULL;
  }
  //-Copia los datos en el puntero recibido como parametro.
  //-Copies data in the pointer received as a parameter.
  if(PrevData64!=InitData64)memcpy(InitData64,PrevData64,sizeof(ullong)*Size);
}

//==============================================================================
/// Crea indice para ordenacion pero sin modificar los datos pasados.
/// Creates sorting index without modifying the processed data.
//==============================================================================
void JRadixSort::MakeIndex(unsigned size,const unsigned *data,unsigned nbits){
  unsigned* auxdata=new unsigned[size];
  memcpy(auxdata,data,sizeof(unsigned)*size);
  Sort(true,size,auxdata,nbits);
  delete[] auxdata;
}

//==============================================================================
/// Crea indice para ordenacion pero sin modificar los datos pasados.
/// Creates sorting index without modifying the processed data.
//==============================================================================
void JRadixSort::MakeIndex(unsigned size,const ullong *data,unsigned nbits){
  ullong* auxdata=new ullong[size];
  memcpy(auxdata,data,sizeof(ullong)*size);
  Sort(true,size,auxdata,nbits);
  delete[] auxdata;
}

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
template<class T> void JRadixSort::TSortData(unsigned size,const T *data,T *result){
  const int threads=omp_get_max_threads();
  if(!Index)Run_Exceptioon("There is no index to sort data.");
  if(size!=Size)Run_Exceptioon("The size of data is invalid.");
  T *res=result;
  if(data==res){//-Reserva vector auxiliar para la ordenacion. //-Allocates auxiliary array for sorting.
    try{
      res=new T[size];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Cannot allocate the requested memory.");
    }
  }

  //-Reordena data[] en res[]
  //-Reorders data[] in res[]
  if(!UseOmp || threads<2){//-Secuencial. //-Sequential
    for(unsigned c2=0;c2<Size;c2++)res[c2]=data[Index[c2]]; 
  }
  else{//-con OpenMP. //-with OpenMP
    const int nk=int(Size/OMPSIZE)+1;
    if(nk<0)Run_Exceptioon("Number of values is invalid.");
    const int rk=int(Size%OMPSIZE);
    #ifdef OMP_USE_RADIXSORT
      #pragma omp parallel for schedule (static)
    #endif
    for(int c=0;c<nk;c++){
      const unsigned c2ini=OMPSIZE*c;
      const unsigned c2fin=c2ini+(c+1<nk? OMPSIZE: rk);
      for(unsigned c2=c2ini;c2<c2fin;c2++)res[c2]=data[Index[c2]]; 
    }
  }
  //-Coloca resultado y libera memoria.
  //-Copies result and frees memory.
  if(res!=result){
    memcpy(result,res,sizeof(T)*size);
    delete[] res;
  }
}

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const byte *data,byte *result){ TSortData<byte>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const word *data,word *result){ TSortData<word>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const unsigned *data,unsigned *result){ TSortData<unsigned>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const int *data,int *result){ TSortData<int>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const float *data,float *result){ TSortData<float>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const double *data,double *result){ TSortData<double>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const tuint2 *data,tuint2 *result){ TSortData<tuint2>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const tfloat2 *data,tfloat2 *result){ TSortData<tfloat2>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const tfloat3 *data,tfloat3 *result){ TSortData<tfloat3>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const tfloat4 *data,tfloat4 *result){ TSortData<tfloat4>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const tdouble3 *data,tdouble3 *result){ TSortData<tdouble3>(size,data,result); }

//==============================================================================
/// Ordena vector de datos en funcion del Index[] calculado previamente.
/// Reorders data arrays as a function of the previously calculated Index[].
//==============================================================================
void JRadixSort::SortData(unsigned size,const tdouble2 *data,tdouble2 *result){ TSortData<tdouble2>(size,data,result); }

//==============================================================================
/// Comprueba ordenacion de datos.
/// Checks sorting data.
//==============================================================================
void JRadixSort::DgCheckResult32()const{
  unsigned p=1;
  for(;p<Size&&InitData32[p-1]<=InitData32[p];p++);
  if(p!=Size)Run_Exceptioon("The order is not correct");
}
//==============================================================================
/// Comprueba ordenacion de datos.
/// Checks sorting data.
//==============================================================================
void JRadixSort::DgCheckResult64()const{
  unsigned p=1;
  for(;p<Size&&InitData64[p-1]<=InitData64[p];p++);
  if(p!=Size)Run_Exceptioon("The order is not correct");
}


