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
/// Muestra informacion de los arrays con memoria reservada.
/// Prints information on arrays with reserved memory.
//==============================================================================
void JArraysCpuList::PrintArraysInfo()const{
  for(unsigned c=0;c<Count;c++)if(PointersAr[c]){
    const JArrayCpu* ar=PointersAr[c];
    Head->Log->Printf("        \"%s\" => ptrid:%2u  ptr:%p"
      ,ar->Name.c_str(),ar->PtrId,ar->Ptr);
    if(ar->PtrId!=c)Run_Exceptioon("PtrId value does not match.");
    if(ar->Ptr!=Pointers[c])Run_Exceptioon("Ptr value does not match.");
  }
}

//==============================================================================
/// Frees allocated memory.
//==============================================================================
void JArraysCpuList::FreeMemory(){
  //-Remove allocations to JArrayCpu objects.
  for(unsigned c=0;c<Count;c++)if(PointersAr[c]){
    if(PointersAr[c]->IsLocked())Run_Exceptioon(fun::PrintStr(
      "Pointer of array %s is locked.",PointersAr[c]->ObjectId().c_str()));
    PointersAr[c]->Free();
  }
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
/// Frees memory of unused allocated arrays.
//==============================================================================
void JArraysCpuList::FreeUnusedArrays(){
  if(GetArrayCountUsed()){
    unsigned nused=0;
    //-Groups arrays in use at the beginning of the list.
    for(unsigned c=0;c<Count;c++){
      if(PointersAr[c]){//-In use.
        if(c!=nused){
          Pointers[nused]=Pointers[c];
          PointersAr[nused]=PointersAr[c];
          PointersAr[nused]->AssignReserve(Pointers[nused],nused);
          Pointers[c]=NULL;
          PointersAr[c]=NULL;
        }
        nused++;
      }
      else{//-Not in use.
        FreeCpuMemory(Pointers[c]);
        Pointers[c]=NULL;
      }
    }
    //-Update variables on unused arrays and Count.
    CountUnused=0;
    for(unsigned c=0;c<Count;c++)PointersUnused[c]=0;
    Count=nused;
  }
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
void JArraysCpuList::Reserve(JArrayCpu* ar){
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
void JArraysCpuList::Free(JArrayCpu* ar){
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
void JArraysCpuList::SwapPtr(JArrayCpu* ar1,JArrayCpu* ar2){
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
  //-Resize arrays and update Ptr in JArrayCpuXXX.
  unsigned nsave=0;
  for(unsigned c=0;c<Count;c++){ 
    void* ptrnew=AllocCpuMemory(newsize);
    if(PointersAr[c]){
      if(PointersAr[c]->IsLocked())Run_Exceptioon(fun::PrintStr(
        "Pointer of array %s is locked.",PointersAr[c]->ObjectId().c_str()));
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
  LastCountArrays=0;
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
void JArraysCpu::PrintAllocMemory(bool emptyblocks,bool arrayinfo)const{ 
  Log->Printf("Allocated memory for %s:",ClassName.c_str());
  for(unsigned a=0;a<MAX_ASIZE;a++){
    const JArraysCpuList* ar=ArraysList[a];
    if(emptyblocks || ar->GetArrayCount() || ar->GetArrayCountUsedMax()){
      Log->Printf("  %2dB x [%u] x %2d = %7.2f MiB (used:%d maxused:%d)",ar->ValueSize,ar->GetArraySize()
        ,ar->GetArrayCount(),double(ar->GetAllocMemoryCpu())/MEBIBYTE
        ,ar->GetArrayCountUsed(),ar->GetArrayCountUsedMax());
      if(arrayinfo)ar->PrintArraysInfo();
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
/// Return the total number of allocated arrays.
//==============================================================================
unsigned JArraysCpu::GetArrayCount()const{
  unsigned n=0;
  for(unsigned a=0;a<MAX_ASIZE;a++)n+=ArraysList[a]->GetArrayCount();
  return(n);
}

//==============================================================================
/// Return the total number of allocated arrays when it was changed.
//==============================================================================
unsigned JArraysCpu::GetArrayCountUpdated(){
  unsigned n=GetArrayCount();
  if(LastCountArrays!=n)LastCountArrays=n;
  else n=0;
  return(n);
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
/// Frees memory of unused allocated arrays.
//==============================================================================
void JArraysCpu::FreeUnusedArrays(){
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->FreeUnusedArrays(); 
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
    PrintAllocMemory(false,false);
  #endif
  if(savesizedata>GetArraySize() || savesizedata>newsize)Run_Exceptioon("Size of data to save is higher than the current size or the new size.");
  //-Resize allocated CPU memory of arrays.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetDataArraySize(newsize,savesizedata); 
  CountResizeDataCalls++;
}


//##############################################################################
//# JArrayCpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArrayCpu::JArrayCpu(const std::string& name,TpTypeData type
  ,JArraysCpu::TpASizeId asize,JArraysCpu* head,bool reserve)
  :Name(name),Type(type),Head(head),ArraysList(NULL)
{
  ArraysList=Head->GetArraysList(asize);
  PtrId=UINT_MAX;
  Ptr=NULL;
  Locked=false;
  #ifdef DG_ARRSPU_PRINT
    Head->Log->Printf("ARRSPU_JArrayCpu>  Create array \'%s\'",Name.c_str());
  #endif
  if(reserve)Reserve();
}
//==============================================================================
/// Destructor.
//==============================================================================
JArrayCpu::~JArrayCpu(){
  Free();
  #ifdef DG_ARRSPU_PRINT
    Head->Log->Printf("ARRSPU_JArrayCpu>  Destroy array \'%s\'",Name.c_str());
  #endif
  Head=NULL;
  ArraysList=NULL;
}
//==============================================================================
/// Returns the object identification for exception messages and debug.
//==============================================================================
std::string JArrayCpu::ObjectId()const{
  string id=fun::PrintStr("JArraysCpu_%dB__",(ArraysList? ArraysList->ValueSize: 0))+Name;
  if(Head && Head->Id>=0)id=id+fun::PrintStr("[c%u]",Head->Id);
  return(id);
}
//==============================================================================
/// Assigns the reserve by JArraysCpuList object.
//==============================================================================
void JArrayCpu::AssignReserve(void* ptr,unsigned ptrid){
  if(Locked && Ptr!=ptr)Run_Exceptioon("Pointer is locked.");
  Ptr=ptr;
  PtrId=ptrid;
}
//==============================================================================
/// Clears the reserve by JArraysCpuList object.
//==============================================================================
void JArrayCpu::ClearReserve(){
  if(Locked)Run_Exceptioon("Pointer is locked.");
  Ptr=NULL;
  PtrId=UINT_MAX;
}
//==============================================================================
/// Swap Ptr and PtrId data between objects.
//==============================================================================
void JArrayCpu::PSwapPtr(JArrayCpu* ar){
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
void JArrayCpu::PMemsetOffset(void* ptr_offset,unsigned offset,byte value
  ,size_t size)
{
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size+offset>GetSize())Run_Exceptioon("Invalid size.");
  memset(ptr_offset,value,Sizeof()*size);
}

//==============================================================================
/// Copy data from src array.
//==============================================================================
void JArrayCpu::PCopyFrom(const JArrayCpu* src,size_t size){
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
void JArrayCpu::PCopyFromOffset(void* dst_ptr,unsigned dst_offset
  ,const JArrayCpu* src,const void* src_ptr,unsigned src_offset,size_t size)
{
  if(!Active() || src==NULL || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(GetValueSize()!=src->GetValueSize())Run_Exceptioon("Size element does not match.");
  if(size+dst_offset>GetSize() || size+src_offset>src->GetSize())Run_Exceptioon("Invalid offset.");
  memcpy(dst_ptr,src_ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data from src pointer.
//==============================================================================
void JArrayCpu::PCopyFromPointer(const void* src_ptr,size_t size){
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  memcpy(Ptr,src_ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data from src pointer using offsets.
//==============================================================================
void JArrayCpu::PCopyFromPointerOffset(void* dst_ptr,unsigned dst_offset
  ,const void* src_ptr,unsigned src_offset,size_t size)
{
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+dst_offset>GetSize())Run_Exceptioon("Invalid offset.");
  memcpy(dst_ptr,src_ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data to dst pointer.
//==============================================================================
void JArrayCpu::PCopyTo(void* dst_ptr,size_t size)const{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  memcpy(dst_ptr,Ptr,Sizeof()*size);
}
//==============================================================================
/// Copy data to dst pointer using offsets.
//==============================================================================
void JArrayCpu::PCopyToOffset(const void* src_ptr,unsigned src_offset
  ,void* dst_ptr,size_t size)const
{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+src_offset>GetSize())Run_Exceptioon("Invalid offset.");
  memcpy(dst_ptr,src_ptr,Sizeof()*size);
}

//==============================================================================
/// Requests an allocated array.
//==============================================================================
void JArrayCpu::Reserve(){
  if(!Ptr){
    ArraysList->Reserve(this);
    #ifdef DG_ARRSPU_PRINT
      Head->Log->Printf("ARRSPU_JArrayCpu>  Reserve array_%uB \'%s\' ptr:%d_%p"
        ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
    #endif
  }
}
//==============================================================================
/// Frees an allocated array.
//==============================================================================
void JArrayCpu::Free(){
  if(Ptr){
    if(Locked)Run_Exceptioon("Pointer is locked.");
    #ifdef DG_ARRSPU_PRINT
      Head->Log->Printf("ARRSPU_JArrayCpu>  Free array_%uB \'%s\' ptr:%d_%p"
        ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
    #endif
    ArraysList->Free(this);
  }
}
//==============================================================================
/// Locks memory pointer in use.
//==============================================================================
void JArrayCpu::LockPtr(){
  if(!Active())Run_Exceptioon("Invalid pointer.");
  #ifdef DG_ARRSPU_PRINT
    Head->Log->Printf("ARRSPU_JArrayCpu>  Lock array_%uB \'%s\' ptr:%d_%p"
      ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
  #endif
  Locked=true;
}
//==============================================================================
/// Release pointer lock.
//==============================================================================
void JArrayCpu::UnlockPtr(){
  #ifdef DG_ARRSPU_PRINT
    Head->Log->Printf("ARRSPU_JArrayCpu>  Unlock array_%uB \'%s\' ptr:%d_%p"
      ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
  #endif
  Locked=false;
}
//==============================================================================
/// Run memset with Ptr.
//==============================================================================
void JArrayCpu::Memset(byte value,size_t size){
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  memset(Ptr,value,Sizeof()*size);
}

