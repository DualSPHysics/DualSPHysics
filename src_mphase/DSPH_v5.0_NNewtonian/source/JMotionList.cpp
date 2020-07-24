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

/// \file JMotionList.cpp \brief Implements the class \ref JMotionList.

#include "JMotionList.h"

using namespace std;

//##############################################################################
//# JMotionListData
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JMotionListData::JMotionListData(){
  ClassName="JMotionListData";
  Active=false;
  TypeSimple=true;
  MvSimple1=MvSimple2=TDouble3(0);
  MvMatrix1=MvMatrix2=TMatrix4d(); //-Indentity matrix.
  VelSimple=AceSimple=TDouble3(0);
}

//==============================================================================
// Anota unico movimiento simple.
//==============================================================================
void JMotionListData::Sp_Movedt(const tdouble3 &mvsimple,double dt){
  Active=true;
  TypeSimple=true;
  MvSimple1=mvsimple;
  MvSimple2=TDouble3(0);
  VelSimple=mvsimple/dt;
  AceSimple=TDouble3(0);
}

//==============================================================================
// Anota unico movimiento con matriz.
//==============================================================================
void JMotionListData::Sp_Movedt(const tmatrix4d &mvmatrix,double dt){
  Active=true;
  TypeSimple=false;
  MvMatrix1=mvmatrix;
  MvMatrix2=TMatrix4d(); //-Indentity matrix.
}

//==============================================================================
// Anota primer movimiento simple.
//==============================================================================
void JMotionListData::Ace2_Move1dt(const tdouble3 &mvsimple){
  Active=true;
  TypeSimple=true;
  MvSimple1=mvsimple;
  MvSimple2=TDouble3(0);
}

//==============================================================================
// Anota segundo movimiento simple.
//==============================================================================
void JMotionListData::Ace2_Move2dt(const tdouble3 &mvsimple){
  if(!Active)Ace2_Move1dt(TDouble3(0));
  if(TypeSimple)MvSimple2=mvsimple;
  else{
    MvMatrix2.a14=mvsimple.x;
    MvMatrix2.a24=mvsimple.y; 
    MvMatrix2.a34=mvsimple.z;
  }
}

//==============================================================================
// Anota primer movimiento con matriz.
//==============================================================================
void JMotionListData::Ace2_Move1dt(const tmatrix4d &mvmatrix){
  Active=true;
  TypeSimple=false;
  MvMatrix1=mvmatrix;
  MvMatrix2=TMatrix4d(); //-Indentity matrix.
}

//==============================================================================
// Anota segundo movimiento con matriz.
//==============================================================================
void JMotionListData::Ace2_Move2dt(const tmatrix4d &mvmatrix){
  if(!Active)Ace2_Move1dt(TMatrix4d());
  if(TypeSimple){
    TypeSimple=false;
    MvMatrix1=TMatrix4d(); //-Indentity matrix.
    MvMatrix1.a14=MvSimple1.x;
    MvMatrix1.a24=MvSimple1.y; 
    MvMatrix1.a34=MvSimple1.z;
  }
  MvMatrix2=mvmatrix;
}

//==============================================================================
// Procesa fin de movimiento.
//==============================================================================
void JMotionListData::Ace2_PosMotion(double dt){
  if(Active && TypeSimple){
    VelSimple=(MvSimple2/(dt+dt));
    AceSimple=((MvSimple2-MvSimple1-MvSimple1)/(dt*dt));
  }
}

//==============================================================================
/// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JMotionListData::GetData(bool &typesimple,tdouble3 &simplemov,tdouble3 &simplevel
  ,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2)const
{
  if(Active){
    if(typesimple=TypeSimple){
      simplemov=MvSimple1;
      simplevel=VelSimple;
      simpleace=AceSimple;
    }
    else{
      matmov=MvMatrix1;
      matmov2=MvMatrix2;
    }
  }
  return(Active);
}

//==============================================================================
/// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JMotionListData::GetData(bool &typesimple,tdouble3 &simplemov,tmatrix4d &matmov)const
{
  if(Active){
    if(typesimple=TypeSimple){
      simplemov=MvSimple1;
    }
    else{
      matmov=MvMatrix1;
    }
  }
  return(Active);
}


//##############################################################################
//# JMotionList
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JMotionList::JMotionList(int nref):Nref(unsigned(nref)){
  ClassName="JMotionList";
  MotionData=new JMotionListData[Nref];
  TimeStep=0;
}

//==============================================================================
// Destructor.
//==============================================================================
JMotionList::~JMotionList(){
  DestructorActive=true;
  delete[] MotionData; MotionData=NULL;
}

//==============================================================================
/// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JMotionList::GetData(unsigned ref,bool &typesimple,tdouble3 &simplemov
  ,tdouble3 &simplevel,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2)const
{
  if(ref>=Nref)Run_Exceptioon("Reference is invalid.");
  return(MotionData[ref].GetData(typesimple,simplemov,simplevel,simpleace,matmov,matmov2));
}

//==============================================================================
/// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JMotionList::GetData(unsigned ref,bool &typesimple,tdouble3 &simplemov,tmatrix4d &matmov)const
{
  if(ref>=Nref)Run_Exceptioon("Reference is invalid.");
  return(MotionData[ref].GetData(typesimple,simplemov,matmov));
}



