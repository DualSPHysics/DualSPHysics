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

/// \file JArraysCpu.cpp \brief Implements the class \ref JArraysCpu.

#include "JArraysCpu.h"
#include "Functions.h"
#include <cstdio>
#include <algorithm>

using namespace std;

//##############################################################################
//# JArraysCpuSize
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArraysCpuSize::JArraysCpuSize(unsigned elementsize):ElementSize(elementsize){
  ClassName="JArraysCpuSize";
  for(unsigned c=0;c<MAXPOINTERS;c++)Pointers[c]=NULL;
  Count=0;
  CountMax=CountUsedMax=0;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JArraysCpuSize::~JArraysCpuSize(){
  DestructorActive=true;
  Reset();
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JArraysCpuSize::Reset(){
  FreeMemory();
  ArraySize=0;
}

//==============================================================================
/// Libera memoria reservada.
/// Frees allocated memory.
//==============================================================================
void JArraysCpuSize::FreeMemory(){
  for(unsigned c=0;c<Count;c++)if(Pointers[c]){ FreePointer(Pointers[c]); Pointers[c]=NULL; }
  CountUsed=Count=0;
}

//==============================================================================
/// Reserva memoria y devuelve puntero con memoria asignada.
/// Allocates memory and returns pointers with allocated memory.
//==============================================================================
void* JArraysCpuSize::AllocPointer(unsigned size)const{
  void* pointer=NULL;
  try{
    switch(ElementSize){
      case 1:   pointer=new char[size];      break;
      case 2:   pointer=new word[size];      break;
      case 4:   pointer=new int[size];       break;
      case 8:   pointer=new double[size];    break;
      case 12:  pointer=new int[size*3];     break;
      case 16:  pointer=new int[size*4];     break;
      case 24:  pointer=new double[size*3];  break;
      case 32:  pointer=new double[size*4];  break;
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Cannot allocate the requested memory.");
  }
  if(!pointer)Run_Exceptioon("The elementsize value is invalid.");
  return(pointer);
}

//==============================================================================
/// Libera la memoria asignada del puntero.
/// Frees memory allocated to pointers.
//==============================================================================
void JArraysCpuSize::FreePointer(void* pointer)const{
  switch(ElementSize){
    case 1:   delete[] ((char*)pointer);    pointer=NULL;   break;
    case 2:   delete[] ((word*)pointer);    pointer=NULL;   break;
    case 4:   delete[] ((int*)pointer);     pointer=NULL;   break;
    case 8:   delete[] ((double*)pointer);  pointer=NULL;   break;
    case 12:  delete[] ((int*)pointer);     pointer=NULL;   break;
    case 16:  delete[] ((int*)pointer);     pointer=NULL;   break;
    case 24:  delete[] ((double*)pointer);  pointer=NULL;   break;
    case 32:  delete[] ((double*)pointer);  pointer=NULL;   break;
  }
  if(pointer)Run_Exceptioon("The elementsize value is invalid.");
}

//==============================================================================
/// Cambia el numero de arrays almacenados. Asignando nuevos arrays o liberando
/// los de los actuales sin uso. 
/// Si count es inferior al numero de los que estan en uso lanza una excepcion.
/// Changes the number of arrays stored. Assigns or releases new arrays if
/// the current are unused.
/// If the count is less than the number of those in use raises an exception.
//==============================================================================
void JArraysCpuSize::SetArrayCount(unsigned count){
  if(count>MAXPOINTERS)Run_Exceptioon("Number of requested arrays exceeds the maximum.");
  if(count<CountUsed)Run_Exceptioon("Unable to free arrays in use.");
  if(ArraySize){
    if(Count<count){//-Genera nuevos arrays. //-Generates new arrays.
      for(unsigned c=Count;c<count;c++)Pointers[c]=AllocPointer(ArraySize);
    }
    if(Count>count){//-Libera arrays. //-Frees arrays.
      for(unsigned c=count;c<Count;c++){ FreePointer(Pointers[c]); Pointers[c]=NULL; }
    }
  }
  Count=count;
  CountMax=max(CountMax,Count);
}

//==============================================================================
/// Cambia el numero de elementos de los arrays.
/// Si hay algun array en uso lanza una excepcion.
/// Changes the number of elements in the arrays.
/// If there is any array in use raises an exception.
//==============================================================================
void JArraysCpuSize::SetArraySize(unsigned size){
  if(CountUsed)Run_Exceptioon("Unable to change the dimension of the arrays because some are in use.");
  if(ArraySize!=size){
    ArraySize=size;
    unsigned count=Count;
    FreeMemory();
    if(count)SetArrayCount(count);
  }
}

//==============================================================================
/// Solicita la reserva de un array.
/// Requests allocating an array.
//==============================================================================
void* JArraysCpuSize::Reserve(){
  if(CountUsed==Count||!ArraySize)Run_Exceptioon(fun::PrintStr("There are no arrays available with %u bytes.",ElementSize));
  CountUsed++;
  CountUsedMax=max(CountUsedMax,CountUsed);
  return(Pointers[CountUsed-1]);
}

//==============================================================================
/// Devuelve la posicion del puntero indicado. Si no existe devuelve MAXPOINTERS.
/// Returns the position of indicated pointer. If it doesn't exist returns MAXPOINTERS.
//==============================================================================
unsigned JArraysCpuSize::FindPointerUsed(void *pointer)const{
  unsigned pos=0;
  for(;pos<CountUsed&&Pointers[pos]!=pointer;pos++);
  return(pos>=CountUsed? MAXPOINTERS: pos);
}

//==============================================================================
/// Libera la reserva de un array.
/// Frees an allocated array.
//==============================================================================
void JArraysCpuSize::Free(void *pointer){
  if(pointer){
    unsigned pos=FindPointerUsed(pointer);
    if(pos==MAXPOINTERS)Run_Exceptioon("The pointer indicated was not reserved.");
    if(pos+1<CountUsed){
      void *aux=Pointers[CountUsed-1]; Pointers[CountUsed-1]=Pointers[pos]; Pointers[pos]=aux;
    }
    CountUsed--;
  }
}  


//##############################################################################
//# JArraysCpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JArraysCpu::JArraysCpu(){
  ClassName="JArraysCpu";
  Arrays1b=new JArraysCpuSize(1);
  Arrays2b=new JArraysCpuSize(2);
  Arrays4b=new JArraysCpuSize(4);
  Arrays8b=new JArraysCpuSize(8);
  Arrays12b=new JArraysCpuSize(12);
  Arrays16b=new JArraysCpuSize(16);
  Arrays24b=new JArraysCpuSize(24);
  Arrays32b=new JArraysCpuSize(32);
}

//==============================================================================
/// Destructor.
//==============================================================================
JArraysCpu::~JArraysCpu(){
  DestructorActive=true;
  delete Arrays1b;
  delete Arrays2b;
  delete Arrays4b;
  delete Arrays8b;
  delete Arrays12b;
  delete Arrays16b;
  delete Arrays24b;
  delete Arrays32b;
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JArraysCpu::Reset(){
  Arrays1b->Reset(); 
  Arrays2b->Reset(); 
  Arrays4b->Reset(); 
  Arrays8b->Reset(); 
  Arrays12b->Reset();
  Arrays16b->Reset();
  Arrays24b->Reset();
  Arrays32b->Reset();
}
 
//==============================================================================
/// Devuelve la cantidad de memoria reservada.
/// Returns amount of allocated memory.
//==============================================================================
llong JArraysCpu::GetAllocMemoryCpu()const{ 
  llong m=Arrays1b->GetAllocMemoryCpu();
  m+=Arrays2b->GetAllocMemoryCpu();
  m+=Arrays4b->GetAllocMemoryCpu();
  m+=Arrays8b->GetAllocMemoryCpu();
  m+=Arrays12b->GetAllocMemoryCpu();
  m+=Arrays16b->GetAllocMemoryCpu();
  m+=Arrays24b->GetAllocMemoryCpu();
  m+=Arrays32b->GetAllocMemoryCpu();
  return(m);
}

//==============================================================================
/// Cambia el numero de elementos de los arrays.
/// Si hay algun array en uso lanza una excepcion.
/// Changes the number of elements in the arrays.
/// If there is any array in use raises an exception.
//==============================================================================
void JArraysCpu::SetArraySize(unsigned size){ 
  //-Frees memory.
  Arrays1b->SetArraySize(0); 
  Arrays2b->SetArraySize(0); 
  Arrays4b->SetArraySize(0); 
  Arrays8b->SetArraySize(0); 
  Arrays12b->SetArraySize(0);
  Arrays16b->SetArraySize(0);
  Arrays24b->SetArraySize(0);
  Arrays32b->SetArraySize(0);
  //-Allocates memory.
  Arrays1b->SetArraySize(size); 
  Arrays2b->SetArraySize(size); 
  Arrays4b->SetArraySize(size); 
  Arrays8b->SetArraySize(size); 
  Arrays12b->SetArraySize(size);
  Arrays16b->SetArraySize(size);
  Arrays24b->SetArraySize(size);
  Arrays32b->SetArraySize(size);
}


