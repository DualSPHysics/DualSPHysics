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

/// \file JDataArrays.cpp \brief Implements the class \ref JDataArrays.

#include "JDataArrays.h"
#include "JException.h"
#include "Functions.h"
#include <cstring>
#include <cstdio>

using namespace std;

//##############################################################################
//# JDataArrays
//##############################################################################
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
      case TypeInt:      delete[] (int     *)ptr;  break;
      case TypeFloat:    delete[] (float   *)ptr;  break;
      case TypeDouble:   delete[] (double  *)ptr;  break;
      case TypeUint3:    delete[] (tuint3  *)ptr;  break;
      case TypeInt3:     delete[] (tint3   *)ptr;  break;
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
/// Returns maximum or minimum number of data values.
//==============================================================================
unsigned JDataArrays::GetDataCount(bool minimum)const{ 
  unsigned smin=0,smax=0;
  const unsigned na=Count();
  if(na)smin=smax=Arrays[0].count;
  if(minimum)for(unsigned c=1;c<na;c++)if(smin>Arrays[c].count)smin=Arrays[c].count;
  else       for(unsigned c=1;c<na;c++)if(smax<Arrays[c].count)smax=Arrays[c].count;
  return(minimum? smin: smax);
}

//==============================================================================
/// Returns number of data values and throws exception when they do not match.
//==============================================================================
unsigned JDataArrays::GetDataCount()const{ 
  unsigned smin=GetDataCount(true );
  unsigned smax=GetDataCount(false);
  if(smin!=smax)Run_Exceptioon("The minimum data count does not match the maximum data count.");
  return(smin);
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
/// Returns constant reference to requested array by idx.
//==============================================================================
const JDataArrays::StDataArray& JDataArrays::GetArrayCte(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
  return(Arrays[idx]);
}

//==============================================================================
/// Returns constant reference to requested array by name.
//==============================================================================
const JDataArrays::StDataArray& JDataArrays::GetArrayCte(const std::string &keyname)const{
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
/// Returns dimension of type of requested array by idx.
//==============================================================================
int JDataArrays::GetArrayDim(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
  return(DimOfType(Arrays[idx].type));
}

//==============================================================================
/// Returns units of requested array by idx.
//==============================================================================
std::string JDataArrays::GetArrayFmt(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
  string fmt=fun::StrSplitValue(":",Arrays[idx].fullname,1);
  if(fmt.empty())fmt=GetFmtByType(Arrays[idx].type);
  return(fmt);
}

//==============================================================================
/// Returns output format according type of array.
//==============================================================================
std::string JDataArrays::GetFmtByType(TpTypeData type){
  string fmt;
  switch(type){
    case TypeUchar:  
    case TypeUshort:
    case TypeUint:     fmt="%u";       break;
    case TypeInt:      fmt="%d";       break;
    case TypeFloat:    fmt="%15.7E";   break;
    case TypeDouble:   fmt="%20.12E";  break;
    case TypeUint3:    fmt="%u";       break;
    case TypeInt3:     fmt="%d";       break;
    case TypeFloat3:   fmt="%15.7E";   break;
    case TypeDouble3:  fmt="%20.12E";  break;
  }
  if(fmt.empty())Run_ExceptioonSta(fun::PrintStr("Type \'%s\' without output-format.",TypeToStr(type)));
  return(fmt);
}

//==============================================================================
/// Returns units of requested array by idx.
//==============================================================================
std::string JDataArrays::GetArrayUnits(unsigned idx)const{
  if(idx>=Count())Run_Exceptioon("Array idx is invalid.");
  string units=fun::StrSplitValue(":",Arrays[idx].fullname,2);
  if(units.empty())units=GetUnitsByName(Arrays[idx].keyname);
  if(units=="NONE")units="";
  return(units);
}

//==============================================================================
/// Returns units according name of array.
//==============================================================================
std::string JDataArrays::GetUnitsByName(std::string keyname){
  const string var=fun::StrLower(keyname);
  if(var=="pos")return(" [m]");
  else if(var=="vel")return(" [m/s]");
  else if(var=="rhop")return(" [kg/m^3]");
  else if(var=="mass")return(" [kg]");
  else if(var=="press")return(" [Pa]");
  else if(var=="vol")return(" [m^3]");
  else if(var=="ace")return(" [m/s^2]");
  else if(var=="vor")return(" [1/s]");
  else if(var=="height")return(" [m]");
  return("");
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
/// Returns dynamic pointer with int array. (this pointer must be deleted)
//==============================================================================
int* JDataArrays::NewArrayInt(unsigned count,bool defvalue,int value){
  try{
    int *v=new int[count];
    if(count && defvalue){
      if(!value)memset(v,0,sizeof(int)*count);
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
/// Returns dynamic pointer with tuint3 array. (this pointer must be deleted)
//==============================================================================
tint3* JDataArrays::NewArrayInt3(unsigned count,bool defvalue,tint3 value){
  try{
    tint3 *v=new tint3[count];
    if(count && defvalue){
      if(value==TInt3(0))memset(v,0,sizeof(tint3)*count);
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

//==============================================================================
/// Move array values according a reindex array.
//==============================================================================
template<class T> void JDataArrays::TReindexData(unsigned sreindex,const unsigned *reindex
  ,unsigned ndata,T *data,T *aux)const
{
  if(aux)memcpy(aux,data,sizeof(T)*ndata); //-Copy current data in auxiliary memory.
  else aux=data; //-Auxiliary memory is not used.
  for(unsigned p=0;p<sreindex;p++){
    const unsigned p0=reindex[p];
    if(p0>=ndata || p>=ndata)Run_Exceptioon("Value number is invalid.");
    if(p!=p0)data[p]=aux[p0];
  }
}

//==============================================================================
/// Apply filter array[count] (1:selected, 0:discarded).
//==============================================================================
unsigned JDataArrays::FilterApply(unsigned count,const byte *filter){
  //-Create index to compact filtered data.
  unsigned *reindex=NewArrayUint(count);
  unsigned nsel=0;
  for(unsigned p=0;p<count;p++)if(filter[p]!=0){
    reindex[nsel++]=p;
  }
  //-Compact filtered data.
  const unsigned na=Count();
  for(unsigned c=0;c<na;c++){
    StDataArray &arr=Arrays[c];
    if(arr.count!=count)Run_Exceptioon(fun::PrintStr("Number of values of \'%s\' does not match size of filter.",arr.keyname.c_str()));
    switch(arr.type){
      case TypeUchar:    ReindexData(nsel,reindex,arr.count,(byte    *)arr.ptr,NULL);  break;
      case TypeUshort:   ReindexData(nsel,reindex,arr.count,(word    *)arr.ptr,NULL);  break;
      case TypeUint:     ReindexData(nsel,reindex,arr.count,(unsigned*)arr.ptr,NULL);  break;
      case TypeInt:      ReindexData(nsel,reindex,arr.count,(int     *)arr.ptr,NULL);  break;
      case TypeFloat:    ReindexData(nsel,reindex,arr.count,(float   *)arr.ptr,NULL);  break;
      case TypeDouble:   ReindexData(nsel,reindex,arr.count,(double  *)arr.ptr,NULL);  break;
      case TypeUint3:    ReindexData(nsel,reindex,arr.count,(tuint3  *)arr.ptr,NULL);  break;
      case TypeInt3:     ReindexData(nsel,reindex,arr.count,(tint3   *)arr.ptr,NULL);  break;
      case TypeFloat3:   ReindexData(nsel,reindex,arr.count,(tfloat3 *)arr.ptr,NULL);  break;
      case TypeDouble3:  ReindexData(nsel,reindex,arr.count,(tdouble3*)arr.ptr,NULL);  break;
      default: Run_Exceptioon(fun::PrintStr("Type of pointer \'%s\' is invalid.",TypeToStr(arr.type)));
    }
    arr.count=nsel;
  }
  //-Free memory.
  delete[] reindex; reindex=NULL;
  return(nsel);
}

//==============================================================================
/// Sort and cut data.
//==============================================================================
unsigned JDataArrays::SortData(unsigned count,const unsigned *reindex){
  const unsigned na=Count();
  //printf("00> count:%u \n",count);
  //-Create auxiliary memory to sort data.
  unsigned maxsize4b=0;
  for(unsigned c=0;c<na;c++){
    const unsigned s4b=(SizeOfType(Arrays[c].type)+3)/4;
    //printf("00b> [%s].%s s4b:%u \n",Arrays[c].keyname.c_str(),TypeToStr(Arrays[c].type),s4b);
    maxsize4b=(maxsize4b>=s4b? maxsize4b: s4b);
  }
  unsigned* aux=NewArrayUint(maxsize4b*GetDataCount(false));
  //printf("01> maxsize4b:%u \n",maxsize4b);
  //-Sort and cut each array data.
  for(unsigned c=0;c<na;c++){
    StDataArray &arr=Arrays[c];
    //printf("01b> [%s].%s count:%u \n",arr.keyname.c_str(),TypeToStr(arr.type),arr.count);
    switch(arr.type){
      case TypeUchar:    ReindexData(count,reindex,arr.count,(byte    *)arr.ptr,(byte    *)aux);  break;
      case TypeUshort:   ReindexData(count,reindex,arr.count,(word    *)arr.ptr,(word    *)aux);  break;
      case TypeUint:     ReindexData(count,reindex,arr.count,(unsigned*)arr.ptr,(unsigned*)aux);  break;
      case TypeInt:      ReindexData(count,reindex,arr.count,(int     *)arr.ptr,(int     *)aux);  break;
      case TypeFloat:    ReindexData(count,reindex,arr.count,(float   *)arr.ptr,(float   *)aux);  break;
      case TypeDouble:   ReindexData(count,reindex,arr.count,(double  *)arr.ptr,(double  *)aux);  break;
      case TypeUint3:    ReindexData(count,reindex,arr.count,(tuint3  *)arr.ptr,(tuint3  *)aux);  break;
      case TypeInt3:     ReindexData(count,reindex,arr.count,(tint3   *)arr.ptr,(tint3   *)aux);  break;
      case TypeFloat3:   ReindexData(count,reindex,arr.count,(tfloat3 *)arr.ptr,(tfloat3 *)aux);  break;
      case TypeDouble3:  ReindexData(count,reindex,arr.count,(tdouble3*)arr.ptr,(tdouble3*)aux);  break;
      default: Run_Exceptioon(fun::PrintStr("Type of pointer \'%s\' is invalid.",TypeToStr(arr.type)));
    }
    arr.count=count;
  }
  //-Free memory.
  delete[] aux; aux=NULL;
  return(count);
}

//==============================================================================
/// Filter list of values according its memory position.
//==============================================================================
unsigned JDataArrays::FilterList(unsigned n,const unsigned *list){
  const unsigned count=GetDataCount(false);
  if(count!=GetDataCount(true))Run_Exceptioon("All arrays must have the same number of values.");
  byte *filter=NewArrayByte(count,true);
  for(unsigned c=0;c<n;c++)if(list[c]<count)filter[list[c]]=1;
  const unsigned nfinal=FilterApply(count,filter);
  delete[] filter; filter=NULL;
  return(nfinal);
}

//==============================================================================
/// Sort and filter list of values according its memory position.
//==============================================================================
unsigned JDataArrays::FilterSortList(unsigned n,const unsigned *list){
  const unsigned count=GetDataCount(false);
  if(count!=GetDataCount(true))Run_Exceptioon("All arrays must have the same number of values.");
  if(n>count)Run_Exceptioon("Size of list is higher than number of values.");
  SortData(n,list);
  return(n);
}


