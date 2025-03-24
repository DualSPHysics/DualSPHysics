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

/// \file JArraysGpu.cpp \brief Implements the class \ref JArraysGpu.

#include "JArraysGpu.h"
#include "Functions.h"
#include "JLog2.h"
#include <climits>
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

//##############################################################################
//# JArraysGpuList
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArraysGpuList::JArraysGpuList(JArraysGpu* head,unsigned valuesize)
  :Head(head),ValueSize(valuesize)
{
  if(Head->Id>=0)ClassName=fun::PrintStr("JArraysGpuList_%dB[g%u]",valuesize,Head->Id);
  else ClassName=fun::PrintStr("JArraysGpuList_%dB",valuesize);
  Count=0;
  for(unsigned c=0;c<MAX_ARRAYS;c++){
    Pointers  [c]=NULL;
    PointersAr[c]=NULL;
    PointersUnused[c]=0;
  }
  //-Variables for saving GPU data during resizing memory.
  for(unsigned c=0;c<MAX_ARRAYS;c++){
    SvDataCpu[c]=NULL;
    SvPointersAr[c]=NULL;
  }
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JArraysGpuList::~JArraysGpuList(){
  Reset();
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JArraysGpuList::Reset(){
  AutoAddArrays=true;
  ResetDataCpu();
  FreeMemory();
  CountUsedMax=0;
  ArraySize=0;
}
 
//==============================================================================
/// Muestra informacion de los arrays con memoria reservada.
/// Prints information on arrays with reserved memory.
//==============================================================================
void JArraysGpuList::PrintArraysInfo()const{
  for(unsigned c=0;c<Count;c++)if(PointersAr[c]){
    const JArrayGpu* ar=PointersAr[c];
    Head->Log->Printf("        \"%s\" => ptrid:%2u  ptr:%p"
      ,ar->Name.c_str(),ar->PtrId,ar->Ptr);
    if(ar->PtrId!=c)Run_Exceptioon("PtrId value does not match.");
    if(ar->Ptr!=Pointers[c])Run_Exceptioon("Ptr value does not match.");
  }
}

//==============================================================================
/// Frees CPU memory allocated and initialise variables for saving GPU data.
//==============================================================================
void JArraysGpuList::ResetDataCpu(){
  for(unsigned c=0;c<MAX_ARRAYS;c++){
    FreeCpuMemory(SvDataCpu[c]);
    SvDataCpu[c]=NULL;
    SvPointersAr[c]=NULL;
  }
  SizeDataCpu=CountDataCpu=SvCount=0;
}

//==============================================================================
/// Frees allocated memory.
//==============================================================================
void JArraysGpuList::FreeMemory(){
  //-Remove allocations to JArrayGpu objects.
  for(unsigned c=0;c<Count;c++)if(PointersAr[c]){
    if(PointersAr[c]->IsLocked())Run_Exceptioon(fun::PrintStr(
      "Pointer of array %s is locked.",PointersAr[c]->ObjectId().c_str()));
    PointersAr[c]->Free();
  }
  //-Free GPU memory in Pointers[].
  for(unsigned c=0;c<Count;c++)if(Pointers[c]){ 
    cudaFree(Pointers[c]); Pointers[c]=NULL;
    Check_CudaErroor("Failed to free GPU memory.");
  }
  Count=0;
  CountUnused=0;
}

//==============================================================================
/// Allocates CPU memory and returns pointer with allocated memory.
//==============================================================================
void* JArraysGpuList::AllocCpuMemory(unsigned size)const{
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
void JArraysGpuList::FreeCpuMemory(void* pointer)const{
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
void JArraysGpuList::SetArrayCount(unsigned count){
  if(count>MAX_ARRAYS)Run_Exceptioon(fun::PrintStr("Number of requested arrays exceeds the maximum (MAX_ARRAYS=%u).",MAX_ARRAYS));
  if(count<Count)Run_Exceptioon("Unable to reduce the number of allocated arrays.");
  if(ArraySize){
    if(count>Count){//-Generates new arrays.
      #ifdef DG_ARRXPU_PRINT
        Head->Log->Printf("ARRXPU_JArraysGpuList>  Allocate new arrays_%uB: %d-%d  (arraysize:%u)"
          ,ValueSize,Count,count-1,ArraySize);
      #endif
      for(unsigned c=Count;c<count;c++){
        cudaMalloc((void**)(Pointers+c),size_t(ValueSize)*ArraySize);
        Check_CudaErroor("Failed to allocate GPU memory.");
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
void JArraysGpuList::FreeUnusedArrays(){
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
        cudaFree(Pointers[c]);
        Pointers[c]=NULL;
        Check_CudaErroor("Failed to free GPU memory.");
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
void JArraysGpuList::SetArraySize(unsigned size){
  if(Count!=CountUnused)Run_Exceptioon("Unable to change the dimension of the arrays because some of them are in use.");
  if(ArraySize!=size){
    #ifdef DG_ARRXPU_PRINT
      Head->Log->Printf("ARRXPU_JArraysGpuList>  Set ArraySize of %uB from %u to %u (count:%u)"
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
void JArraysGpuList::Reserve(JArrayGpu* ar){
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
void JArraysGpuList::Free(JArrayGpu* ar){
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
void JArraysGpuList::SwapPtr(JArrayGpu* ar1,JArrayGpu* ar2){
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
/// Saves current GPU data on CPU memory buffers and free GPU memory.
//==============================================================================
void JArraysGpuList::SaveDataArray(unsigned savesizedata){
  if(SizeDataCpu || CountDataCpu || SvCount)
    Run_Exceptioon("Error: Variables for saving GPU data are in use.");
  SizeDataCpu=savesizedata;
  if(SizeDataCpu>ArraySize)Run_Exceptioon("Size of data to save is invalid.");
  //-Copy GPU data in use to CPU memory.
  unsigned nsave=0;
  for(unsigned c=0;c<Count;c++){ 
    if(PointersAr[c]){
      if(PointersAr[c]->IsLocked())Run_Exceptioon(fun::PrintStr(
        "Pointer of array %s is locked.",PointersAr[c]->ObjectId().c_str()));
      SvPointersAr[nsave]=PointersAr[c];
      SvDataCpu[nsave]=AllocCpuMemory(SizeDataCpu);
      if(SizeDataCpu){
        cudaMemcpy(SvDataCpu[nsave],Pointers[c],size_t(ValueSize)*SizeDataCpu,cudaMemcpyDeviceToHost);
        Check_CudaErroor("Failed to save GPU data to CPU buffer.");
      }
      nsave++;
    }
  }
  CountDataCpu=nsave;
  SvCount=Count;
  //-Free pointers in use and GPU memory.
  FreeMemory();
}

//==============================================================================
/// Allocate GPU memory and recover data from CPU memory.
//==============================================================================
void JArraysGpuList::RestoreDataArray(unsigned newsize){
  if(Count)Run_Exceptioon("Error: Object should be empty.");
  SetArraySize(newsize);
  if(SizeDataCpu>ArraySize)Run_Exceptioon("Size of data to load is higher than allocated memory.");
  //-Allocates GPU memory for previous number of arrays.
  SetArrayCount(SvCount);
  if(CountDataCpu>Count)Run_Exceptioon("Error: Number of arrays to load is higher than allocated arrays.");
  //-Allocates GPU memory in SvPointersAr.
  for(unsigned cs=0;cs<CountDataCpu;cs++){
    SvPointersAr[cs]->Reserve();
    if(SizeDataCpu){
      cudaMemcpy(SvPointersAr[cs]->Ptr,SvDataCpu[cs],size_t(ValueSize)*SizeDataCpu,cudaMemcpyHostToDevice);
      Check_CudaErroor("Failed to restore GPU data from CPU buffer.");
    }
  }
  //-Free CPU data buffers.
  ResetDataCpu();
}


//##############################################################################
//# JArraysGpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArraysGpu::JArraysGpu(int gpuid,JLog2* log,int id):Log(log),GpuId(gpuid),Id(id){
  ClassName=(Id>=0? fun::PrintStr("JArraysGpu[g%u]",Id): string("JArraysGpu"));
  for(unsigned a=0;a<MAX_ASIZE;a++){
    ArraysList[a]=new JArraysGpuList(this,GetTypeSize(TpASizeId(a)));
  }
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JArraysGpu::~JArraysGpu(){
  for(unsigned a=0;a<MAX_ASIZE;a++){
    delete ArraysList[a]; ArraysList[a]=NULL;
  }
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JArraysGpu::Reset(){
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->Reset(); 
  CountResizeDataCalls=0;
  LastCountArrays=0;
}
 
//==============================================================================
/// Devuelve la cantidad de memoria reservada.
/// Returns amount of allocated memory.
//==============================================================================
llong JArraysGpu::GetAllocMemoryGpu()const{ 
  llong m=0;
  for(unsigned a=0;a<MAX_ASIZE;a++)m+=ArraysList[a]->GetAllocMemoryGpu();
  return(m);
}
 
//==============================================================================
/// Muestra la memoria reservada de cada tipo de array.
/// Prints the allocated memory by each array.
//==============================================================================
void JArraysGpu::PrintAllocMemory(bool emptyblocks,bool arrayinfo)const{ 
  Log->Printf("Allocated memory for %s:",ClassName.c_str());
  for(unsigned a=0;a<MAX_ASIZE;a++){
    const JArraysGpuList* ar=ArraysList[a];
    if(emptyblocks || ar->GetArrayCount() || ar->GetArrayCountUsedMax()){
      Log->Printf("  %2dB x [%u] x %2d = %7.2f MiB (used:%d maxused:%d)",ar->ValueSize,ar->GetArraySize()
        ,ar->GetArrayCount(),double(ar->GetAllocMemoryGpu())/MEBIBYTE
        ,ar->GetArrayCountUsed(),ar->GetArrayCountUsedMax());
      if(arrayinfo)ar->PrintArraysInfo();
    }
  }
  Log->Printf("  Total = %.2f MiB",double(GetAllocMemoryGpu())/MEBIBYTE);
}

//==============================================================================
/// Activa/desactiva la creacion automanita de arrays al llamar Reserve().
/// Active/deactivate the automatic creation of new arrays calling Reserve().
//==============================================================================
void JArraysGpu::SetAutoAddArrays(bool active){
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetAutoAddArrays(active);
}

//==============================================================================
/// Return the total number of allocated arrays.
//==============================================================================
unsigned JArraysGpu::GetArrayCount()const{
  unsigned n=0;
  for(unsigned a=0;a<MAX_ASIZE;a++)n+=ArraysList[a]->GetArrayCount();
  return(n);
}

//==============================================================================
/// Return the total number of allocated arrays when it was changed.
//==============================================================================
unsigned JArraysGpu::GetArrayCountUpdated(){
  unsigned n=GetArrayCount();
  if(LastCountArrays!=n)LastCountArrays=n;
  else n=0;
  return(n);
}

//==============================================================================
/// Asigna memoria mas arrays de la lista indicada (ArraysList).
/// Add new array allocations to the indicated ArraysList.
//==============================================================================
void JArraysGpu::AddArrayCount(TpASizeId asize,unsigned count){
  JArraysGpuList* arlis=ArraysList[asize];
  arlis->SetArrayCount(arlis->GetArrayCount()+count);
}

//==============================================================================
/// Frees memory of unused allocated arrays.
//==============================================================================
void JArraysGpu::FreeUnusedArrays(){
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->FreeUnusedArrays(); 
}

//==============================================================================
/// Cambia el numero de elementos de los arrays.
/// Si hay algun array en uso lanza una excepcion.
/// Changes the number of elements in the arrays.
/// If there is any array in use raises an exception.
//==============================================================================
void JArraysGpu::SetArraySize(unsigned size){ 
  //-Frees memory.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetArraySize(0); 
  //-Allocates memory.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SetArraySize(size); 
}

//==============================================================================
/// Cambia el numero de elementos de los arrays guardando previamente los datos 
/// en CPU para despues recupearlos.
/// Changes the number of elements in the arrays keeping the data to recover 
/// after resize of GPU memory.
//==============================================================================
void JArraysGpu::SetDataArraySize(unsigned newsize,unsigned savesizedata){
  #ifdef DG_ARRXPU_PRINT
    Log->Printf("ARRXPU_JArraysGpu>  SetDataArraySize( newsize:%d savesizedata:%d )",newsize,savesizedata);
    PrintAllocMemory(false,false);
  #endif
  if(savesizedata>GetArraySize() || savesizedata>newsize)Run_Exceptioon("Size of data to save is higher than the current size or the new size.");
  //-Saves data on CPU memory and free GPU memory.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->SaveDataArray(savesizedata); 
  //-Allocate new GPU memory size and recover data from CPU memory.
  for(unsigned a=0;a<MAX_ASIZE;a++)ArraysList[a]->RestoreDataArray(newsize);
  CountResizeDataCalls++;
}


//##############################################################################
//# JArrayGpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArrayGpu::JArrayGpu(const std::string& name,TpTypeData type
  ,JArraysGpu::TpASizeId asize,JArraysGpu* head,bool reserve)
  :Name(name),Type(type),Head(head),ArraysList(NULL)
{
  ArraysList=Head->GetArraysList(asize);
  PtrId=UINT_MAX;
  Ptr=NULL;
  Locked=false;
  DataCpu=NULL;
  #ifdef DG_ARRXPU_PRINT
    Head->Log->Printf("ARRXPU_JArrayGpu>  Create array \'%s\'",Name.c_str());
  #endif
  if(reserve)Reserve();
}
//==============================================================================
/// Destructor.
//==============================================================================
JArrayGpu::~JArrayGpu(){
  Free();
  #ifdef DG_ARRXPU_PRINT
    Head->Log->Printf("ARRXPU_JArrayGpu>  Destroy array \'%s\'",Name.c_str());
  #endif
  Head=NULL;
  ArraysList=NULL;
}
//==============================================================================
/// Returns the object identification for exception messages and debug.
//==============================================================================
std::string JArrayGpu::ObjectId()const{
  string id=fun::PrintStr("JArraysGpu_%dB__",(ArraysList? ArraysList->ValueSize: 0))+Name;
  if(Head && Head->Id>=0)id=id+fun::PrintStr("[g%u]",Head->Id);
  return(id);
}
//==============================================================================
/// Assigns the reserve by JArraysGpuList object.
//==============================================================================
void JArrayGpu::AssignReserve(void* ptr,unsigned ptrid){
  if(Locked && Ptr!=ptr)Run_Exceptioon("Pointer is locked.");
  Ptr=ptr;
  PtrId=ptrid;
}
//==============================================================================
/// Clears the reserve by JArraysGpuList object.
//==============================================================================
void JArrayGpu::ClearReserve(){
  if(Locked)Run_Exceptioon("Pointer is locked.");
  Ptr=NULL;
  PtrId=UINT_MAX;
}
//==============================================================================
/// Swap Ptr and PtrId data between objects.
//==============================================================================
void JArrayGpu::PSwapPtr(JArrayGpu* ar){
  if(this!=ar){
    if(ar==NULL)Run_Exceptioon("Invalid array for swap.");
    if(ArraysList!=ar->ArraysList)
      Run_Exceptioon(fun::PrintStr("Swap between Array_%uB \'%s\' and Array_%uB \'%s\' is invalid because they belong to different array lists."
        ,ArraysList->ValueSize,Name.c_str(),ar->ArraysList->ValueSize,ar->Name.c_str()));
    if(DataCpu)DataFree();
    if(ar->DataCpu)ar->DataFree();
    ArraysList->SwapPtr(this,ar);
  }
}
//==============================================================================
/// Run cudaMemset with Ptr using offset.
//==============================================================================
void JArrayGpu::PMemsetOffset(void* ptr_offset,unsigned offset,byte value
  ,size_t size)
{
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size+offset>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemset(ptr_offset,value,Sizeof()*size);
}
//==============================================================================
/// Run cudaMemsetAsync with Ptr using offset.
//==============================================================================
void JArrayGpu::PMemsetAsyncOffset(void* ptr_offset,unsigned offset,byte value
  ,size_t size,cudaStream_t stm)
{
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size+offset>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemsetAsync(ptr_offset,value,Sizeof()*size,stm);
}

//==============================================================================
/// Copy data from src array.
//==============================================================================
void JArrayGpu::PCopyFrom(const JArrayGpu* src,size_t size){
  if(this!=src){
    if(!Active() || src==NULL || !src->Active())Run_Exceptioon("Invalid arrays or pointers.");
    if(GetValueSize()!=src->GetValueSize())Run_Exceptioon("Size element does not match.");
    if(size>GetSize() || size>src->GetSize())Run_Exceptioon("Invalid size.");
    cudaMemcpy(Ptr,src->Ptr,Sizeof()*size,cudaMemcpyDeviceToDevice);
  }
}
//==============================================================================
/// Copy data from src array using offsets.
//==============================================================================
void JArrayGpu::PCopyFromOffset(void* dst_ptr,unsigned dst_offset
  ,const JArrayGpu* src,const void* src_ptr,unsigned src_offset,size_t size)
{
  if(!Active() || src==NULL || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(GetValueSize()!=src->GetValueSize())Run_Exceptioon("Size element does not match.");
  if(size+dst_offset>GetSize() || size+src_offset>src->GetSize())Run_Exceptioon("Invalid offset.");
  cudaMemcpy(dst_ptr,src_ptr,Sizeof()*size,cudaMemcpyDeviceToDevice);
}
//==============================================================================
/// Copy data from src pointer.
//==============================================================================
void JArrayGpu::PCopyFromPointer(const void* src_ptr,size_t size){
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(Ptr,src_ptr,Sizeof()*size,cudaMemcpyDeviceToDevice);
}
//==============================================================================
/// Copy data from src pointer using offsets.
//==============================================================================
void JArrayGpu::PCopyFromPointerOffset(void* dst_ptr,unsigned dst_offset
  ,const void* src_ptr,unsigned src_offset,size_t size)
{
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+dst_offset>GetSize())Run_Exceptioon("Invalid offset.");
  cudaMemcpy(dst_ptr,src_ptr,Sizeof()*size,cudaMemcpyDeviceToDevice);
}

//==============================================================================
/// Copy data to dst pointer.
//==============================================================================
void JArrayGpu::PCopyTo(void* dst_ptr,size_t size)const{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(dst_ptr,Ptr,Sizeof()*size,cudaMemcpyDeviceToDevice);
}
//==============================================================================
/// Copy data to dst pointer using offsets.
//==============================================================================
void JArrayGpu::PCopyToOffset(const void* src_ptr,unsigned src_offset
  ,void* dst_ptr,size_t size)const
{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+src_offset>GetSize())Run_Exceptioon("Invalid offset.");
  cudaMemcpy(dst_ptr,src_ptr,Sizeof()*size,cudaMemcpyDeviceToDevice);
}

//==============================================================================
/// Copy gpu data to cpu array (dst).
//==============================================================================
void JArrayGpu::PCopyToHost(JArrayCpu* dst,size_t size)const{
  if(!Active() || dst==NULL || !dst->Active())Run_Exceptioon("Invalid arrays or pointers.");
  if(GetValueSize()!=dst->GetValueSize())Run_Exceptioon("Size element does not match.");
  if(size>GetSize() || size>dst->GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(dst->ptrvoid(),Ptr,Sizeof()*size,cudaMemcpyDeviceToHost);
}
//==============================================================================
/// Copy gpu data to cpu pointer (dst_ptr).
//==============================================================================
void JArrayGpu::PCopyToHostPointer(void* dst_ptr,size_t size)const{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(dst_ptr,Ptr,Sizeof()*size,cudaMemcpyDeviceToHost);
}
//==============================================================================
/// Copy gpu data to cpu array (dst) using offset.
//==============================================================================
void JArrayGpu::PCopyToHostOffset(const void* src_ptr,unsigned src_offset
  ,JArrayCpu* dst,void* dst_ptr,unsigned dst_offset,size_t size)const
{
  if(!Active() || dst==NULL || !dst->Active())Run_Exceptioon("Invalid arrays or pointers.");
  if(GetValueSize()!=dst->GetValueSize())Run_Exceptioon("Size element does not match.");
  if(size+src_offset>GetSize() || size+dst_offset>dst->GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(dst_ptr,src_ptr,Sizeof()*size,cudaMemcpyDeviceToHost);
}
//==============================================================================
/// Copy gpu data to cpu pointer (dst_ptr) using offset.
//==============================================================================
void JArrayGpu::PCopyToHostPointerOffset(const void* src_ptr,unsigned src_offset
  ,void* dst_ptr,unsigned dst_offset,size_t size)const
{
  if(!Active() || dst_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+src_offset>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(dst_ptr,src_ptr,Sizeof()*size,cudaMemcpyDeviceToHost);
}

//==============================================================================
/// Copy from cpu array (src) to GPU.
//==============================================================================
void JArrayGpu::PCopyFromHost(const JArrayCpu* src,size_t size){
  if(!Active() || src==NULL || !src->Active())Run_Exceptioon("Invalid arrays or pointers.");
  if(GetValueSize()!=src->GetValueSize())Run_Exceptioon("Size element does not match.");
  if(size>GetSize() || size>src->GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(Ptr,src->cptrvoid(),Sizeof()*size,cudaMemcpyHostToDevice);
}
//==============================================================================
/// Copy from cpu pointer (src_ptr) to GPU.
//==============================================================================
void JArrayGpu::PCopyFromHostPointer(const void* src_ptr,size_t size){
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(Ptr,src_ptr,Sizeof()*size,cudaMemcpyHostToDevice);
}
//==============================================================================
/// Copy from cpu array (src) to GPU using offset.
//==============================================================================
void JArrayGpu::PCopyFromHostOffset(void* dst_ptr,unsigned dst_offset
  ,const JArrayCpu* src,const void* src_ptr,unsigned src_offset,size_t size)
{
  if(!Active() || src==NULL || !src->Active())Run_Exceptioon("Invalid arrays or pointers.");
  if(GetValueSize()!=src->GetValueSize())Run_Exceptioon("Size element does not match.");
  if(size+dst_offset>GetSize() || size+src_offset>src->GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(dst_ptr,src_ptr,Sizeof()*size,cudaMemcpyHostToDevice);
}
//==============================================================================
/// Copy gpu data to cpu pointer (dst_ptr) using offset.
//==============================================================================
void JArrayGpu::PCopyFromHostPointerOffset(void* dst_ptr,unsigned dst_offset
  ,const void* src_ptr,unsigned scr_offset,size_t size)
{
  if(!Active() || src_ptr==NULL)Run_Exceptioon("Invalid arrays or pointers.");
  if(size+dst_offset>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemcpy(dst_ptr,src_ptr,Sizeof()*size,cudaMemcpyHostToDevice);
}


//==============================================================================
/// Requests an allocated array (and frees CPU memory).
//==============================================================================
void JArrayGpu::Reserve(){
  if(!Ptr){
    if(DataCpu)DataFree();
    ArraysList->Reserve(this);
    #ifdef DG_ARRXPU_PRINT
      Head->Log->Printf("ARRXPU_JArrayGpu>  Reserve array_%uB \'%s\' ptr:%d_%p"
        ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
    #endif
  }
}
//==============================================================================
/// Frees an allocated array (and frees CPU memory).
//==============================================================================
void JArrayGpu::Free(){
  if(DataCpu)DataFree();
  if(Ptr){
    if(Locked)Run_Exceptioon("Pointer is locked.");
    #ifdef DG_ARRXPU_PRINT
      Head->Log->Printf("ARRXPU_JArrayGpu>  Free array_%uB \'%s\' ptr:%d_%p"
        ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
    #endif
    ArraysList->Free(this);
  }
}
//==============================================================================
/// Locks memory pointer in use.
//==============================================================================
void JArrayGpu::LockPtr(){
  if(!Active())Run_Exceptioon("Invalid pointer.");
  #ifdef DG_ARRXPU_PRINT
    Head->Log->Printf("ARRSXPU_JArrayGpu>  Lock array_%uB \'%s\' ptr:%d_%p"
      ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
  #endif
  Locked=true;
}
//==============================================================================
/// Release pointer lock.
//==============================================================================
void JArrayGpu::UnlockPtr(){
  #ifdef DG_ARRXPU_PRINT
    Head->Log->Printf("ARRXPU_JArrayGpu>  Unlock array_%uB \'%s\' ptr:%d_%p"
      ,ArraysList->ValueSize,Name.c_str(),PtrId,Ptr);
  #endif
  Locked=false;
}
//==============================================================================
/// Run cudaMemset with Ptr.
//==============================================================================
void JArrayGpu::CuMemset(byte value,size_t size){
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemset(Ptr,value,Sizeof()*size);
}
//==============================================================================
/// Run cudaMemsetAsync with Ptr.
//==============================================================================
void JArrayGpu::CuMemsetAsync(byte value,size_t size,cudaStream_t stm){
  if(!Active())Run_Exceptioon("Invalid pointer.");
  if(size>GetSize())Run_Exceptioon("Invalid size.");
  cudaMemsetAsync(Ptr,value,Sizeof()*size,stm);
}


//==============================================================================
/// Allocates CPU memory for data copy.
//==============================================================================
void JArrayGpu::DataAlloc(){
  if(!Ptr)Run_Exceptioon(fun::PrintStr("Array_%uB \'%s\' is not active."
    ,ArraysList->ValueSize,Name.c_str()));
  if(!DataCpu)DataCpu=ArraysList->AllocCpuMemory(GetSize());
}
//==============================================================================
/// Frees allocated CPU memory for data.
//==============================================================================
void JArrayGpu::DataFree(){
  if(DataCpu){
    ArraysList->FreeCpuMemory(DataCpu);
    DataCpu=NULL;
  }
}
//==============================================================================
/// Allocates CPU memory and copy data from GPU memory.
//==============================================================================
void JArrayGpu::DataDown(unsigned ndata){
  if(ndata>GetSize())Run_Exceptioon(fun::PrintStr("Array_%uB \'%s\' does not have the requested amount of data."
    ,ArraysList->ValueSize,Name.c_str()));
  DataAlloc();
  if(ndata)cudaMemcpy(DataCpu,Ptr,size_t(GetValueSize())*ndata,cudaMemcpyDeviceToHost);
}

