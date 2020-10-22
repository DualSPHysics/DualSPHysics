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

/// \file JRangeFilter.cpp \brief Implements the class \ref JRangeFilter.

#include "JRangeFilter.h"
#include "Functions.h"
#include <cstdlib>
#include <cstring>
#include <climits>
#include <algorithm>

//##############################################################################
//# JRangeFilter
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JRangeFilter::JRangeFilter(std::string filter){
  ClassName="JRangeFilter";
  Ranges=NULL; FastValue=NULL;
  Reset();
  Config(filter);
}

//==============================================================================
/// Constructor.
//==============================================================================
JRangeFilter::JRangeFilter(const std::vector<unsigned> &values){
  ClassName="JRangeFilter";
  Ranges=NULL; FastValue=NULL;
  Reset();
  Config(values);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JRangeFilter::Reset(){
  ResizeRanges(0);
  ValueMin=1; ValueMax=0;
  delete[] FastValue; FastValue=NULL;
}

//==============================================================================
/// Allocates memory for ranges of values.
//==============================================================================
void JRangeFilter::ResizeRanges(unsigned size){
  unsigned* rg=NULL;
  if(size)rg=new unsigned[size*2];
  if(size&&Size)memcpy(rg,Ranges,sizeof(unsigned)*2*Size);
  delete[] Ranges;
  Ranges=rg;
  Size=size;
  Count=std::min(Count,Size);
}

//==============================================================================
/// Checks whether a new value passes the filter.
//==============================================================================
bool JRangeFilter::CheckNewValue(unsigned v)const{
  bool ret=false;
  for(unsigned c=0;c<Count&&!ret;c++)ret=(Ranges[c<<1]<=v&&v<=Ranges[(c<<1)+1]);
  return(ret);
}

//==============================================================================
/// Adds a value.
//==============================================================================
void JRangeFilter::AddValue(unsigned v){
  if(!CheckNewValue(v)){
    if(Count+1>Size)ResizeRanges(Size+50);
    Ranges[Count<<1]=v;
    Ranges[(Count<<1)+1]=v;
    Count++;
  }
}

//==============================================================================
/// Adds a range of values.
//==============================================================================
void JRangeFilter::AddRange(unsigned v,unsigned v2){
  if(v==v2)AddValue(v);
  else if(v<v2){
    int len;
    do{
      len=v2-v+1;
      for(int c=0;c<int(Count)&&v<=v2;c++){
        unsigned r=Ranges[c<<1],r2=Ranges[(c<<1)+1];
        if(r<=v&&v<=r2)v=r2+1;
        if(r<=v2&&v2<=r2)v2=r-1;
      }
    }while(len!=v2-v+1&&v<=v2);
    if(v<=v2){
      if(Count+1>Size)ResizeRanges(Size+50);
      Ranges[Count<<1]=v;
      Ranges[(Count<<1)+1]=v2;
      Count++;
    }
  }
}

//==============================================================================
/// Adds a range of values with step.
//==============================================================================
void JRangeFilter::AddRangeStep(unsigned v,unsigned v2,unsigned step){ 
  if(step==1)AddRange(v,v2);
  else if(step>0)for(unsigned c=v;v<=c && c<=v2;c=(c+step>c? c+step: v2+1))AddValue(c);
}

//==============================================================================
/// Reorders ranges of values.
//==============================================================================
void JRangeFilter::SortRanges(){
  for(int c=0;c<int(Count)-1;c++){
    unsigned r=Ranges[c<<1];
    for(int c2=c+1;c2<int(Count);c2++)if(r>Ranges[c2<<1]){
      Ranges[c<<1]=Ranges[c2<<1];  Ranges[c2<<1]=r; 
      r=Ranges[(c<<1)+1];
      Ranges[(c<<1)+1]=Ranges[(c2<<1)+1];  Ranges[(c2<<1)+1]=r; 
      r=Ranges[c<<1];
    }
  }
}

//==============================================================================
/// Merges ranges of values.
//==============================================================================
void JRangeFilter::JoinRanges(){
  if(Count){
    unsigned n=0;
    unsigned r=Ranges[0],r2=Ranges[1];
    for(unsigned c=1;c<Count;c++){
      unsigned x=Ranges[c<<1],x2=Ranges[(c<<1)+1];
      if(r2+1==x)r2=x2;
      else{
        Ranges[n<<1]=r; Ranges[(n<<1)+1]=r2; n++;
        r=x; r2=x2;
      }
    }
    Ranges[n<<1]=r; Ranges[(n<<1)+1]=r2;
    Count=n+1;
  }
}

//==============================================================================
/// Returns ranges of values in text format.
//==============================================================================
std::string JRangeFilter::ToString()const{
  std::string tx="";
  for(unsigned c=0;c<Count;c++){
    unsigned r=Ranges[c<<1],r2=Ranges[(c<<1)+1];
    if(!tx.empty())tx=tx+',';
    if(r==r2)tx=tx+fun::UintStr(r);
    else tx=tx+fun::UintStr(r)+"-"+fun::UintStr(r2);
  } 
  return(tx);
}

//==============================================================================
/// Configures the given filter.
//==============================================================================
void JRangeFilter::Config(std::string filter){
  Reset();
  while(!filter.empty()){
    //printf("filter --> [%s]\n",filter.c_str());
    std::string tx=fun::StrSplit(",",filter);
    if(!tx.empty()){
      //printf("tx --> [%s]\n",tx.c_str());
      int pos=int(tx.find("-"));
      std::string tx1=(pos<0? tx: (pos>0? tx.substr(0,pos): ""));
      std::string tx2=(pos>0? tx.substr(pos+1): "");
      std::string tx3="";
      pos=int(tx2.find(":"));
      if(pos>0){
        tx3=(pos>0? tx2.substr(pos+1): "");
        tx2=(pos>0? tx2.substr(0,pos): "");
      }
      //printf("tx:[%s]-[%s]-[%s]\n",tx1.c_str(),tx2.c_str(),tx3.c_str());
      if(!tx1.empty()){
        if(tx2.empty())AddValue(atoi(tx1.c_str()));
        else if(tx3.empty())AddRange(atoi(tx1.c_str()),atoi(tx2.c_str()));
        else AddRangeStep(atoi(tx1.c_str()),atoi(tx2.c_str()),atoi(tx3.c_str()));
      }
    }
  }
  Prepare();
}

//==============================================================================
/// Configures the given filter.
//==============================================================================
void JRangeFilter::Config(const std::vector<unsigned> &values){
  Reset();
  const unsigned nv=unsigned(values.size());
  if(nv){
    unsigned v1=values[0],v2=values[0];
    for(unsigned cv=1;cv<nv;cv++){
      const unsigned v=values[cv];
      if(v==v2 || v==v2+1)v2=v;
      else{
        AddRange(v1,v2);
        v1=v2=v;
      }
    }
    AddRange(v1,v2);
  }
  Prepare();
}

//==============================================================================
/// Prepares data to use.
//==============================================================================
void JRangeFilter::Prepare(){
  SortRanges();
  JoinRanges();
  if(Count){
    ValueMin=Ranges[0]; ValueMax=Ranges[((int(Count)-1)<<1)+1];
    if(ValueMax-ValueMin<1000 && ValueMax-ValueMin>1 && Count>1){
      FastValue=new byte[ValueMax-ValueMin+1];
      memset(FastValue,0,sizeof(byte)*(ValueMax-ValueMin+1));
      for(unsigned c=0;c<Count;c++){
        unsigned r=Ranges[c<<1],r2=Ranges[(c<<1)+1];
        for(;r<=r2;r++)FastValue[r-ValueMin]=1;
      }
    }
  }
}

//==============================================================================
/// Checks whether a value passes the filter.
//==============================================================================
bool JRangeFilter::CheckValue(unsigned v)const{
  return(ValueMin<=v&&v<=ValueMax&&( Count==1||(FastValue&&FastValue[v-ValueMin])||(!FastValue&&CheckNewValue(v)) ));
}

//==============================================================================
/// Returns the first valid value.
/// Returns UINT_MAX if there is not.
//==============================================================================
unsigned JRangeFilter::GetFirstValue()const{
  return(Count? Ranges[0]: UINT_MAX);
}

//==============================================================================
/// Returns the next valid value to the given one.
/// Returns UINT_MAX if there is not.
//==============================================================================
unsigned JRangeFilter::GetNextValue(unsigned v)const{
  unsigned ret=UINT_MAX;
  for(unsigned c=0;c<Count && ret==UINT_MAX;c++){
    unsigned r=Ranges[c<<1],r2=Ranges[(c<<1)+1];
    if(v<r)ret=r;
    else if(v<r2)ret=v+1;
  }
  return(ret);
}

//==============================================================================
/// Returns vector with values.
//==============================================================================
void JRangeFilter::GetValues(std::vector<unsigned> &values)const{
  unsigned v=GetFirstValue();
  while(v!=UINT_MAX){
    values.push_back(v);
    v=GetNextValue(v);
  }
}


