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

/// \file JArraysCpu.cpp \brief Implements the class \ref JArraysCpu.

#include "JArraysCpu.h"
#include "Functions.h"
#include "JLog2.h"
#include <climits>
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

//##############################################################################
//# JArraysCpuList
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArraysCpuList::JArraysCpuList(JArraysCpu* head,unsigned valuesize)
  :Head(head),ValueSize(valuesize)
{
  if(Head->Id>=0)ClassName=fun::PrintStr("JArraysCpuList_%dB[g%u]",valuesize,Head->Id);
  else ClassName=fun::PrintStr("JArraysCpuList_%dB",valuesize);
  Count=0;
  for(unsigned c=0;c<MAX_ARRAYS;c++){
    Pointers  [c]=NULL;
    PointersAr[c]=NULL;
    PointersUnused[c]=0;
  }
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JArraysCpuList::~JArraysCpuList(){
  Reset();
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JArraysCpuList::Reset(){
  AutoAddArrays=true;
  FreeMemory();
  CountUsedMax=0;
  ArraySize=0;
}

//==============================================================================
/// Frees allocated memory.
//==============================================================================
void JArraysCpuList::FreeMemory(){
  //-Remove allocations to JArraySpu objects.
  for(unsigned c=0;c<Count;c++)if(PointersAr[c])PointersAr[c]->Free();
  //-Free CPU memory in Pointers[].
  for(unsigned c=0;c<Count;c++)if(Pointers[c]){
    FreeCpuMemory(Pointers[c]);
    Pointers[c]=NULL;
  }
  Count=0;
  CountUnused=0;
}

//==============================================================================
/// Allocates CPU memory and returns pointer with allocated memory.
//==============================================================================
void* JArraysCpuList::AllocCpuMemory(unsigned size)const{
  void* pointer=NULL;
  try{
    switch(ValueSize){
      case 1:   pointer=new char[size];      break;
      case 2:   pointer=new word[size];      break;
      case 4:   pointer=new int[size];       break;
      case 8:   pointer=new double[size];    break;
      case 12:  pointer=new int[size*3];     break;
      case 16:  pointer=new int[size*4];     break;
      case 24:  pointer=new double[size*3];  break;
      case 32:  pointer=new double[size*4];  break;
      case 72:  pointer=new double[size*9];  break;
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Cannot allocate the requested memory.");
  }
  if(!pointer)Run_Exceptioon("The size of value is invalid.");
  return(pointer);
}

//==============================================================================
/// Frees CPU memory allocated to pointers.
//==============================================================================
void JArraysCpuList::FreeCpuMemory(void* pointer)const{
  switch(ValueSize){
    case 1:   delete[] ((char*)pointer);    pointer=NULL;   break;
    case 2:   delete[] ((word*)pointer);    pointer=NULL;   break;
    case 4:   delete[] ((int*)pointer);     pointer=NULL;   break;
    case 8:   delete[] ((double*)pointer);  pointer=NULL;   break;
    case 12:  delete[] ((int*)pointer);     pointer=NULL;   break;
    case 16:  delete[] ((int*)pointer);     pointer=NULL;   break;
    case 24:  delete[] ((double*)pointer);  pointer=NULL;   break;
    case 32:  delete[] ((double*)pointer);  pointer=NULL;   break;
    case 72:  delete[] ((double*)pointer);  pointer=NULL;   break;
  }
  if(pointer)Run_Exceptioon("The size of value is invalid.");
}

//==============================================================================
/// Increases the number of arrays stored.
//==============================================================================
void JArraysCpuList::SetArrayCount(unsigned count){
  if(count>MAX_ARRAYS)Run_Exceptioon(fun::PrintStr("Number of requested arrays exceeds the maximum (MAX_ARRAYS=%u).",MAX_ARRAYS));
  if(count<Count)Run_Exceptioon("Unable to reduce the number of allocated arrays.");
  if(ArraySize){
    if(count>Count){//-Generates new arrays.
      #ifdef DG_ARRSPU_PRINT
        Head->Log->Printf("ARRSPU_JArraysCpuList>  Allocate new arrays_%uB: %d-%d  (arraysize:%u)"
          ,ValueSize,Count,count-1,ArraySize);
      #endif
      for(unsigned c=Count;c<count;c++){
        Pointers[c]=AllocCpuMemory(ArraySize);
        PointersAr[c]=NULL;
        PointersUnused[CountUnused]=c;
        CountUnused++;
      }
      Count=count;
    }
  }
  else Count=CountUnused=count;
}

//==============================================================================
/// Changes the number of allocated elements in the arrays.
/// If there is any array in use raises an exception.
//==============================================================================
void JArraysCpuList::SetArraySize(unsigned size){
  if(Count!=CountUnused)Run_Exceptioon("Unable to change the dimension of the arrays because some of them are in use.");
  if(ArraySize!=size){
    #ifdef DG_ARRSPU_PRINT
      Head->Log->Printf("ARRSPU_JArraysCpuList>  Set ArraySize of %uB from %u to %u (count:%u)"
        ,ValueSize,ArraySize,size,Count);
    #endif
    ArraySize=size;
    const unsigned countpre=Count;
    FreeMemory();
    if(countpre)SetArrayCount(countpre);
  }
}

//==============================================================================
/// Requests an allocated array.
//==============================================================================
void JArraysCpuList::Reserve(JArraySpu* ar){
  if(!ArraySize)Run_Exceptioon(fun::PrintStr("There are no allocated memory for array \'%s\'.",ar->Name.c_str()));
  if(!CountUnused && AutoAddArrays)SetArrayCount(Count+1);
  if(!CountUnused)Run_Exceptioon(fun::PrintStr("There are no arrays available of %u bytes for array \'%s\'.",ValueSize,ar->Name.c_str()));
  //-Reserve last available array.
  CountUnused--;
  const unsigned ptrid=PointersUnused[CountUnused];
  if(PointersAr[ptrid])Run_Exceptioon(fun::PrintStr("Error: Array was not available for array \'%s\'.",ar->Name.c_str()));
  PointersAr[ptrid]=ar;
  ar->AssignReserve(Pointers[ptrid],ptrid);
  CountUsedMax=max(CountUsedMax,GetArrayCountUsed());
}

//==============================================================================
/// Frees an allocated array.
//==============================================================================
void JArraysCpuList::Free(JArraySpu* ar){
  const unsigned ptrid=ar->PtrId;
  if(ptrid>=Count || ar->Ptr!=Pointers[ptrid] || ar!=PointersAr[ptrid])
    Run_Exceptioon(fun::PrintStr("Error: Array was not reserved for array \'%s\'.",ar->Name.c_str()));
  PointersAr[ptrid]=NULL;
  PointersUnused[CountUnused]=ptrid;
  CountUnused++;
  ar->ClearReserve();
}  

//==============================================================================
/// Swap Ptr and PtrId data between objects.
//==============================================================================
void JArraysCpuList::SwapPtr(JArraySpu* ar1,JArraySpu* ar2){
  void* ptr1=ar1->Ptr;
  unsigned ptrid1=ar1->PtrId;
  void* ptr2=ar2->Ptr;
  unsigned ptrid2=ar2->PtrId;
  ar1->AssignReserve(ptr2,ptrid2);
  if(ptrid2<Count)PointersAr[ptrid2]=ar1;
  ar2->AssignReserve(ptr1,ptrid1);
  if(ptrid1<Count)PointersAr[ptrid1]=ar2;
}

//==============================================================================
/// Changes the number of elements in the arrays keeping the data  after resize
/// of CPU memory.
//==============================================================================
void JArraysCpuList::SetDataArraySize(unsigned newsize,unsigned savesizedata){
  if(savesizedata>ArraySize)Run_Exceptioon("Size of data to save is invalid.");
  //-Resize arrays and update Ptr in JArraySpuXXX.
  unsigned nsave=0;
  for(unsigned c=0;c<Count;c++){ 
    void* ptrnew=AllocCpuMemory(newsize);
    if(PointersAr[c]){
      if(savesizedata)memcpy(ptrnew,Pointers[c],size_t(ValueSize)*savesizedata);
      FreeCpuMemory(Pointers[c]);
      Pointers[c]=ptrnew;
      PointersAr[c]->AssignReserve(ptrnew,c);
    }
    else{
      FreeCpuMemory(Pointers[c]);
      Pointers[c]=ptrnew;
    }
  }
  ArraySize=newsize;
}


//##############################################################################
//# JArraysCpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArraysCpu::JArraysCpu(JLog2* log,int id):Log(log),Id(id){
  ClassName=(Id>=0? fun::PrintStr("JArraysCpu[c%u]",Id): string("JArraysCpu"));
  for(unsigned a=0;a<MAX_ASIZE;a++){
    ArraysList[a]=new JArraysCpuList(this,GetTypeSize(TpASizeId(a)));
  }
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JArraysCpu::~JArraysCpu(){
  for(unsigned a=0;a<MAX_ASIZE;a++){
    delete ArraysList[a]; ArraysList[a]=NULL;
  }
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JArraysCpu::Reset(){
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->Reset(); 
  CountResizeDataCalls=0;
}
 
//==============================================================================
/// Devuelve la cantidad de memoria reservada.
/// Returns amount of allocated memory.
//==============================================================================
llong JArraysCpu::GetAllocMemoryCpu()const{ 
  llong m=0;
  for(unsigned a=0;a<MAX_ASIZE;a++)m+=ArraysList[a]->GetAllocMemoryCpu();
  return(m);
}
 
//==============================================================================
/// Muestra la memoria reservada de cada tipo de array.
/// Prints the allocated memory by each array.
//==============================================================================
void JArraysCpu::PrintAllocMemory(bool all)const{ 
  Log->Printf("Allocated memory for %s:",ClassName.c_str());
  for(unsigned a=0;a<MAX_ASIZE;a++){
    const JArraysCpuList* ar=ArraysList[a];
    if(all || ar->GetArrayCount() || ar->GetArrayCountUsedMax()){
      Log->Printf("  %2dB x [%u] x %2d = %7.2f MiB (used:%d maxused:%d)",ar->ValueSize,ar->GetArraySize()
        ,ar->GetArrayCount(),double(ar->GetAllocMemoryCpu())/MEBIBYTE
        ,ar->GetArrayCountUsed(),ar->GetArrayCountUsedMax());
    }
  }
  Log->Printf("  Total = %.2f MiB",double(GetAllocMemoryCpu())/MEBIBYTE);
}

//==============================================================================
/// Activa/desactiva la creacion automanita de arrays al llamar Reserve().
/// Active/deactivate the automatic creation of new arrays calling Reserve().
//==============================================================================
void JArraysCpu::SetAutoAddArrays(bool active){
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetAutoAddArrays(active);
}

//==============================================================================
/// Asigna memoria mas arrays de la lista indicada (ArraysList).
/// Add new array allocations to the indicated ArraysList.
//==============================================================================
void JArraysCpu::AddArrayCount(TpASizeId asize,unsigned count){
  JArraysCpuList* arlis=ArraysList[asize];
  arlis->SetArrayCount(arlis->GetArrayCount()+count);
}

//==============================================================================
/// Cambia el numero de elementos de los arrays.
/// Si hay algun array en uso lanza una excepcion.
/// Changes the number of elements in the arrays.
/// If there is any array in use raises an exception.
//==============================================================================
void JArraysCpu::SetArraySize(unsigned size){ 
  //-Frees memory.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetArraySize(0); 
  //-Allocates memory.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetArraySize(size); 
}

//==============================================================================
/// Changes the number of elements in the arrays keeping the data after resize
/// of CPU memory.
//==============================================================================
void JArraysCpu::SetDataArraySize(unsigned newsize,unsigned savesizedata){
  #ifdef DG_ARRSPU_PRINT
    Log->Printf("ARRSPU_JArraysCpu>  SetDataArraySize( newsize:%d savesizedata:%d )",newsize,savesizedata);
    PrintAllocMemory(false);
  #endif
  if(savesizedata>GetArraySize() || savesizedata>newsize)Run_Exceptioon("Size of data to save is higher than the current size or the new size.");
  //-Resize allocated CPU memory of arrays.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetDataArraySize(newsize,savesizedata); 
  CountResizeDataCalls++;
}


//##############################################################################
//# JArraySpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArraySpu::JArraySpu(const std::string& name,TpTypeData type
  ,JArraysCpu::TpASizeId asize,JArraysCpu* head,bool reserve)
  :Name(name),Type(type),Head(head),ArraysList(NULL)
{
  ArraysList=Head->GetArraysList(asize);
  PtrId=UINT_MAX;
  Ptr=NULL;
  #ifdef DG_ARRSPU_PRINT
    Head->Log->Printf("ARRSPU_JArraySpu>  Create array \'%s\'",Name.c_str());
  #endif
  if(reserve)Reserve();
}
//==============================================================================
/// Destructor.
//==============================================================================
JArraySpu::~JArraySpu(){
  Free();
  #ifdef DG_ARRSPU_PRINT
    Head->Log->Printf("ARRSPU_JArraySpu>  Destroy array \'%s\'",Name.c_str());
  #endif
  Head=NULL;
  ArraysList=NULL;
}
//==============================================================================
/// Returns the object identification for exception message and debug.
//==============================================================================
std::string JArraySpu::ObjectId()const{
  string id=fun::PrintStr("JArraysCpu_%dB.",(ArraysList? ArraysList->ValueSize: 0))+Name;
  if(Head && Head->Id>=0)id=id+fun::PrintStr("[c%u]",Head->Id);
  return(id);
}
//==============================================================================
/// Assigns the reserve by JArraysCpuList object.
//==============================================================================
void JArraySpu::AssignReserve(void* ptr,unsigned ptrid){
  Ptr=ptr;
  PtrId=ptrid;
}
//==============================================================================
/// Clears the reserve by JArraysCpuList object.
//==============================================================================
void JArraySpu::ClearReserve(){
  Ptr=NULL;
  PtrId=UINT_MAX;
}
//==============================================================================
/// Swap Ptr and PtrId data between objects.
//==============================================================================
void JArraySpu::PSwapPtr(JArraySpu* ar){
  if(this!=ar){
    if(ar==NULL)Run_Exceptioon("Invalid array for swap.");
    if(ArraysList!=ar->ArraysList)
      Run_Exceptioon(fun::PrintStr("Swap between Array_%uB \'%s\' and Array_%uB \'%s\' is invalid because they belong to different array lists."
        ,ArraysList->ValueSize,Name.c_str(),ar->ArraysList->ValueSize,ar->Name.c_str()));
    ArraysList->SwapPtr(this,ar);
  }
}
//==============================================================================
/// Run memset with Ptr using offset.
//==============================================================================
void JArraySpu::PMemsetOffset(void* ptr_offset,unsigned offset,byte value
  ,size_t size)
{
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size+offset>GetSize())Run_Exceptioon("Invalid size.");
  memset(ptr_offset,value,Sizeof()*size);
}

//==============================================================================
/// Copy data from src array.
//==============================================================================
void JArraySpu::PCopyFrom(const JArraySpu* src,size_t size){
  if(this!=src){
    if(!Active() || src==NULL || !src->Active())Run_Exceptioon("Invalid arrays or pointers.");
    if(GetValueSize()!=src->GetValueSize())Run_Exceptioon("Size element does not match.");
    if(size>GetSize() || size>src->GetSize())Run_Exceptioon("Invalid size.");
    memcpy(Ptr,src->Ptr,Sizeof()*size);
  }
}
//==============================================================================
/// Copy data from src array using offsets.
//==============================================================================
void JArraySpu::PCopyFromOffset(void* dst_ptr,unsigned dst_offset
  ,const JArraySpu* src,const void* src_ptr,unsigned src_offset,size_t size)
{
  if(!Active() || src==NULL || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(GetValueSize()!=src->GetValueSize())Run_Exceptioon("Size element does not match.");
  if(size+dst_offset>GetSize() || size+src_offset>src->GetSize())Run_Exceptioon("Invalid offset.");
  memcpy(dst_ptr,src_ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data from src pointer.
//==============================================================================
void JArraySpu::PCopyFromPointer(const void* src_ptr,size_t size){
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  memcpy(Ptr,src_ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data from src pointer using offsets.
//==============================================================================
void JArraySpu::PCopyFromPointerOffset(void* dst_ptr,unsigned dst_offset
  ,const void* src_ptr,unsigned src_offset,size_t size)
{
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+dst_offset>GetSize())Run_Exceptioon("Invalid offset.");
  memcpy(dst_ptr,src_ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data to dst pointer.
//==============================================================================
void JArraySpu::PCopyTo(void* dst_ptr,size_t size)const{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  memcpy(dst_ptr,Ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data to dst pointer using offsets.
//==============================================================================
void JArraySpu::PCopyToOffset(const void* src_ptr,unsigned src_offset
  ,void* dst_ptr,size_t size)const
{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+src_offset>GetSize())Run_Exceptioon("Invalid offset.");
  memcpy(dst_ptr,src_ptr,Sizeof()*size);
}

//==============================================================================
/// Requests an allocated array.
//==============================================================================
void JArraySpu::Reserve(){
  if(!Ptr){
    ArraysList->Reserve(this);
    #ifdef DG_ARRSPU_PRINT
      Head->Log->Printf("ARRSPU_JArraySpu>  Reserve array_%uB \'%s\' ptr:%d_%p"
        ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
    #endif
  }
}
//==============================================================================
/// Frees an allocated array.
//==============================================================================
void JArraySpu::Free(){
  if(Ptr){
    #ifdef DG_ARRSPU_PRINT
      Head->Log->Printf("ARRSPU_JArraySpu>  Free array_%uB \'%s\' ptr:%d_%p"
        ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
    #endif
    ArraysList->Free(this);
  }
}
//==============================================================================
/// Run memset with Ptr.
//==============================================================================
void JArraySpu::Memset(byte value,size_t size){
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  memset(Ptr,value,Sizeof()*size);
}

