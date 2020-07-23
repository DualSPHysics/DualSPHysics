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

/// \file JMeanValues.cpp \brief Implements the class \ref JMeanValue and class \ref JMeanMoving.

#include "JMeanValues.h"
#include <cmath>

using namespace std;


//##############################################################################
//# JMeanMoving
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMeanMoving::JMeanMoving(unsigned size){
  ClassName="JMeanMoving";
  Values=NULL; Weights=NULL;
  InitSimple(size);
}

//==============================================================================
/// Destructor.
//==============================================================================
JMeanMoving::~JMeanMoving(){
  DestructorActive=true;
  Reset();
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JMeanMoving::Reset(){
  delete[] Values;  Values=NULL;
  delete[] Weights; Weights=NULL;
  SizeValues=NextValue=0;
  ValuesFull=false;
}
 
//==============================================================================
/// Configures type and size of mean.
//==============================================================================
void JMeanMoving::Init(unsigned size,bool weighted){ 
  Reset();
  if(size){
    SizeValues=size;
    try{
      Values=new double[SizeValues];
      if(weighted)Weights=new double[SizeValues];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Cannot allocate the requested memory.");
    }
  }
}
 
//==============================================================================
/// Configures size of simple mean.
//==============================================================================
void JMeanMoving::InitSimple(unsigned size){ 
  Init(size,false);
}
 
//==============================================================================
/// Configures size of weighted linear mean.
//==============================================================================
void JMeanMoving::InitWeightedLinear(unsigned size){ 
  Init(size,true);
  double sum=0;
  for(unsigned c=0;c<SizeValues;c++){
    Weights[c]=double(c+1);
    sum+=Weights[c];
  }
  for(unsigned c=0;c<SizeValues;c++)Weights[c]=Weights[c]/sum;
}
 
//==============================================================================
/// Configures size of weighted exponential mean.
//==============================================================================
void JMeanMoving::InitWeightedExponential(unsigned size,float fac){ 
  Init(size,true);
  double sum=0;
  for(unsigned c=0;c<SizeValues;c++){
    Weights[c]=(exp(double(c+1)*fac/SizeValues)-1.);
    sum+=Weights[c];
  }
  for(unsigned c=0;c<SizeValues;c++)Weights[c]=Weights[c]/sum;
}
 
//==============================================================================
/// Adds new value to mean.
//==============================================================================
void JMeanMoving::AddValue(double v){ 
  if(SizeValues){
    Values[NextValue]=v;
    NextValue=(NextValue+1)%SizeValues;
    if(!NextValue)ValuesFull=true;
  }
}  

//==============================================================================
/// Returns simple mean of values.
//==============================================================================
double JMeanMoving::GetSimpleMean()const{
  double sum=0;
  const unsigned n=(ValuesFull? SizeValues: NextValue);
  for(unsigned c=0;c<n;c++)sum+=Values[c];
  return(n? sum/n: 0);
}

//==============================================================================
/// Returns weighted mean of values.
//==============================================================================
double JMeanMoving::GetWeightedMean()const{
  double sum=0;
  const unsigned n=(ValuesFull? SizeValues: NextValue);
  if(n && Weights){
    const unsigned ini=(ValuesFull? NextValue: 0);
    for(unsigned c=0;c<n;c++)sum+=(Values[(ini+c)%SizeValues]*Weights[c]);
    if(!ValuesFull){ 
      double wsum=0;
      for(unsigned c=0;c<n;c++)wsum+=Weights[c];
      sum=sum/wsum;
    }
  }
  return(sum);
}


