//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2019 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDataArrays.cpp \brief Implements the class \ref JDataArrays.

#include "JDataArrays.h"
#include "JException.h"
#include "Functions.h"
#include <cstring>
#include <cstdio>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JDataArrays::JDataArrays(){
  ClassName="JDataArrays";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDataArrays::~JDataArrays(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Throws exception related to a file from a static method.
//==============================================================================
void JDataArrays::RunExceptioonStatic(const std::string &srcfile,int srcline
  ,const std::string &method
  ,const std::string &msg,const std::string &file)
{
  throw JException(srcfile,srcline,"JDataArrays",method,msg,file);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDataArrays::Reset(){
  FreeMemory();
  Arrays.clear();
}

//==============================================================================
/// Copy data from other JDataArrays but delptr is changed to false.
//==============================================================================
void JDataArrays::CopyFrom(const JDataArrays &arr){
  Reset();
  const unsigned na=arr.Count();
  for(unsigned ca=0;ca<na;ca++){
    StDataArray ar=arr.Arrays[ca];
    ar.delptr=false;
    Arrays.push_back(ar);
  }
}

//==============================================================================
/// Frees dynamic memory of array pointer using delete[] when delptr=true.
//==============================================================================
void JDataArrays::FreeMemory(StDataArray &arr){
  if(arr.delptr){
    void *ptr=arr.ptr;
    switch(arr.type){
      case TypeUchar:    delete[] (byte    *)ptr;  break;
      case TypeUshort:   delete[] (word    *)ptr;  break;
      case TypeUint:     delete[] (unsigned*)ptr;  break;
      case TypeFloat:    delete[] (float   *)ptr;  break;
      case TypeDouble:   delete[] (double  *)ptr;  break;
      case TypeUint3:    delete[] (tuint3  *)ptr;  break;
      case TypeFloat3:   delete[] (tfloat3 *)ptr;  break;
      case TypeDouble3:  delete[] (tdouble3*)ptr;  break;
      default: Run_Exceptioon(fun::PrintStr("Type of pointer \'%s\' is invalid.",TypeToStr(arr.type)));
    }
    arr.ptr=NULL;
    arr.count=0;
  }
}

//==============================================================================
/// Frees pointers with delptr=true using delete[].
//==============================================================================
void JDataArrays::FreeMemory(){
  const unsigned na=Count();
  for(unsigned c=0;c<na;c++)if(Arrays[c].delptr)FreeMemory(Arrays[c]);
}

//==============================================================================
/// Returns idx of array with given name.
//==============================================================================
unsigned JDataArrays::GetIdxName(const std::string &keyname)const{ 
  unsigned idx=0;
  const unsigned na=Count();
  for(;idx<na && Arrays[idx].keyname!=keyname;idx++);
  return(idx>=na? UINT_MAX: idx);
}

//==============================================================================
/// Returns idx of array with given name.
//==============================================================================
std::string JDataArrays::CheckErrorArray(const std::string &keyname,TpTypeData type,unsigned count)const{ 
  string err;
  unsigned idx=GetIdxName(keyname);
  if(idx==UINT_MAX)err=fun::PrintStr("Array \'%s\' is missing.",keyname.c_str());
  else err=CheckErrorArray(idx,type,count);
  return(err);
}

//==============================================================================
/// Returns idx of array with given name.
//==============================================================================
std::string JDataArrays::CheckErrorArray(unsigned idx,TpTypeData type,unsigned count)const{ 
  string err;
  if(idx>=Count())err=fun::PrintStr("Array %u is missing.",idx);
  else{
    const string keyname=Arrays[idx].keyname;
    if(Arrays[idx].type!=type)err=fun::PrintStr("Type of array \'%s\' is not %s.",keyname.c_str(),TypeToStr(type));
    else if(count && Arrays[idx].count!=count)err=fun::PrintStr("Size of array \'%s\' is not %u.",keyname.c_str(),count);
    else if(count && !Arrays[idx].ptr)err=fun::PrintStr("Array \'%s\' without data.",keyname.c_str());
  }
  return(err);
}

//==============================================================================
/// Adds generic array defined by the user.
//==============================================================================
unsigned JDataArrays::AddArray(std::string fullname,TpTypeData type
  ,unsigned count,void *ptr,bool delptr)
{
  StDataArray ar;
  ar.fullname=fullname;
  ar.keyname=fun::StrSplit(":",fullname);
  ar.type=type;
  ar.count=count;
  ar.ptr=ptr;
  ar.delptr=delptr;
  if(ar.keyname.empty())Run_Exceptioon(fun::PrintStr("Array name \'%s\' is invalid.",ar.keyname.c_str()));
  if(ExistsName(ar.keyname))Run_Exceptioon(fun::PrintStr("Array \'%s\' already exists.",ar.keyname.c_str()));
  Arrays.push_back(ar);
  return(Count()-1);
}

//==============================================================================
/// Removes array from vector and dynamic memory is freed when delptr=true.
//==============================================================================
void JDataArrays::DeleteArray(unsigned idx){
  if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
  FreeMemory(Arrays[idx]);
  Arrays.erase(Arrays.begin()+idx);
}

//==============================================================================
/// Removes array from vector and dynamic memory is freed when delptr=true.
//==============================================================================
void JDataArrays::DeleteArray(std::string keyname){
  const unsigned idx=GetIdxName(keyname);
  if(idx==UINT_MAX)Run_Exceptioon(fun::PrintStr("Array \'%s\' is missing.",keyname.c_str()));
  else DeleteArray(idx);
}

//==============================================================================
/// Removes array from vector, but dynamic memory never is freed.
//==============================================================================
void JDataArrays::EraseArray(std::string keyname){
  const unsigned idx=GetIdxName(keyname);
  if(idx==UINT_MAX)Run_Exceptioon(fun::PrintStr("Array \'%s\' is missing.",keyname.c_str()));
  Arrays.erase(Arrays.begin()+idx);
}

//==============================================================================
/// Moves array to other position.
//==============================================================================
void JDataArrays::MoveArray(unsigned idx,unsigned idx2){
  if(idx!=idx2){
    if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
    StDataArray ar=Arrays[idx];
    Arrays.erase(Arrays.begin()+idx);
    if(idx2>idx)idx2--;
    if(idx2>=Count())Arrays.push_back(ar);
    else Arrays.insert(Arrays.begin()+idx2,ar);
  }
}

//==============================================================================
/// Returns reference to requested array by idx.
//==============================================================================
JDataArrays::StDataArray& JDataArrays::GetArray(unsigned idx){
  if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
  return(Arrays[idx]);
}

//==============================================================================
/// Returns reference to requested array by name.
//==============================================================================
JDataArrays::StDataArray& JDataArrays::GetArray(const std::string &keyname){
  const unsigned idx=GetIdxName(keyname);
  if(idx==UINT_MAX)Run_Exceptioon(fun::PrintStr("Array \'%s\' is missing.",keyname.c_str()));
  return(Arrays[idx]);
}

//==============================================================================
/// Returns data of requested array by idx.
//==============================================================================
JDataArrays::StDataArray JDataArrays::GetArrayData(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
  return(Arrays[idx]);
}

//==============================================================================
/// Returns data of requested array by name.
//==============================================================================
JDataArrays::StDataArray JDataArrays::GetArrayData(const std::string &keyname)const{
  const unsigned idx=GetIdxName(keyname);
  if(idx==UINT_MAX)Run_Exceptioon(fun::PrintStr("Array \'%s\' is missing.",keyname.c_str()));
  return(Arrays[idx]);
}

//==============================================================================
/// Returns the pointer to data. Throw exception when idx, type or count is wrong.
//==============================================================================
const void* JDataArrays::GetArrayPtr(unsigned idx,TpTypeData type,unsigned count)const{
  string err=CheckErrorArray(idx,type,count);
  if(!err.empty())Run_Exceptioon(err);
  return(Arrays[idx].ptr);
}

//==============================================================================
/// Returns the pointer to data. Throw exception when keyname, type or count is wrong.
//==============================================================================
const void* JDataArrays::GetArrayPtr(const std::string &keyname,TpTypeData type,unsigned count)const{
  string err=CheckErrorArray(keyname,type,count);
  if(!err.empty())Run_Exceptioon(err);
  return(Arrays[GetIdxName(keyname)].ptr);
}

//==============================================================================
/// Print list of arrays.
//==============================================================================
void JDataArrays::Print()const{
  for(unsigned ca=0;ca<Count();ca++){
    StDataArray ar=Arrays[ca];
    printf("[%02d] \'%s\' (%s)  %s[%u] %s\n",ca,ar.keyname.c_str(),ar.fullname.c_str(),TypeToStr(ar.type),ar.count,(ar.delptr? "(delptr)": ""));
  }
}

//==============================================================================
/// Returns dynamic pointer with byte array. (this pointer must be deleted)
//==============================================================================
byte* JDataArrays::NewArrayByte(unsigned count,bool defvalue,byte value){
  try{
    byte *v=new byte[count];
    if(count && defvalue)memset(v,value,sizeof(byte)*count);
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with word array. (this pointer must be deleted)
//==============================================================================
word* JDataArrays::NewArrayWord(unsigned count,bool defvalue,word value){
  try{
    word *v=new word[count];
    if(count && defvalue){
      if(!value)memset(v,0,sizeof(word)*count);
      else for(unsigned c=0;c<count;c++)v[c]=value;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with unsigned array. (this pointer must be deleted)
//==============================================================================
unsigned* JDataArrays::NewArrayUint(unsigned count,bool defvalue,unsigned value){
  try{
    unsigned *v=new unsigned[count];
    if(count && defvalue){
      if(!value)memset(v,0,sizeof(unsigned)*count);
      else for(unsigned c=0;c<count;c++)v[c]=value;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with float array. (this pointer must be deleted)
//==============================================================================
float* JDataArrays::NewArrayFloat(unsigned count,bool defvalue,float value){
  try{
    float *v=new float[count];
    if(count && defvalue){
      if(!value)memset(v,0,sizeof(float)*count);
      else for(unsigned c=0;c<count;c++)v[c]=value;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with double array. (this pointer must be deleted)
//==============================================================================
double* JDataArrays::NewArrayDouble(unsigned count,bool defvalue,double value){
  try{
    double *v=new double[count];
    if(count && defvalue){
      if(!value)memset(v,0,sizeof(double)*count);
      else for(unsigned c=0;c<count;c++)v[c]=value;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with tuint3 array. (this pointer must be deleted)
//==============================================================================
tuint3* JDataArrays::NewArrayUint3(unsigned count,bool defvalue,tuint3 value){
  try{
    tuint3 *v=new tuint3[count];
    if(count && defvalue){
      if(value==TUint3(0))memset(v,0,sizeof(tuint3)*count);
      else for(unsigned c=0;c<count;c++)v[c]=value;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with tfloat3 array. (this pointer must be deleted)
//==============================================================================
tfloat3* JDataArrays::NewArrayFloat3(unsigned count,bool defvalue,tfloat3 value){
  try{
    tfloat3 *v=new tfloat3[count];
    if(count && defvalue){
      if(value==TFloat3(0))memset(v,0,sizeof(tfloat3)*count);
      else for(unsigned c=0;c<count;c++)v[c]=value;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with double array. (this pointer must be deleted)
//==============================================================================
tdouble3* JDataArrays::NewArrayDouble3(unsigned count,bool defvalue,tdouble3 value){
  try{
    tdouble3 *v=new tdouble3[count];
    if(count && defvalue){
      if(value==TDouble3(0))memset(v,0,sizeof(tdouble3)*count);
      else for(unsigned c=0;c<count;c++)v[c]=value;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Returns dynamic pointer with a sequence. (this pointer must be deleted)
//==============================================================================
unsigned* JDataArrays::NewArraySeqUint(unsigned count,unsigned start,unsigned step){
  try{
    unsigned *v=new unsigned[count];
    unsigned vv=start;
    for(unsigned c=0;c<count;c++){
      v[c]=vv; vv+=step;
    }
    return(v);
  }
  catch(const std::bad_alloc){
    Run_ExceptioonSta(fun::PrintStr("Could not allocate the requested memory (size=%u).",count));
  }
}

//==============================================================================
/// Obtains x,y,z values from tfloat4 values.
//==============================================================================
void JDataArrays::ToFloat3xyz(unsigned count,const tfloat4 *data,tfloat3 *dest){
  for(unsigned c=0;c<count;c++){
    const tfloat4 v=data[c];
    dest[c]=TFloat3(v.x,v.y,v.z);
  }
}

//==============================================================================
/// Obtains w values from tfloat4 values.
//==============================================================================
void JDataArrays::ToFloat1w(unsigned count,const tfloat4 *data,float *dest){
  for(unsigned c=0;c<count;c++)dest[c]=data[c].w;
}

//==============================================================================
/// Returns dynamic pointer with x,y,z values from tfloat4 values. (this pointer must be deleted)
//==============================================================================
tfloat3* JDataArrays::NewArrayFloat3xyz(unsigned count,const tfloat4 *data){
  tfloat3 *v=NewArrayFloat3(count,false);
  ToFloat3xyz(count,data,v);
  return(v);
}

//==============================================================================
/// Returns dynamic pointer with x,y,z values from tfloat4 values. (this pointer must be deleted)
//==============================================================================
float* JDataArrays::NewArrayFloat1w(unsigned count,const tfloat4 *data){
  float *v=NewArrayFloat(count,false);
  ToFloat1w(count,data,v);
  return(v);
}



