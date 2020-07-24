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

/// \file JMotionPos.cpp \brief Implements the class \ref JMotionPos.

#include "JMotionPos.h"

//##############################################################################
//# JMotionPos
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JMotionPos::JMotionPos(){
  PosSimple=TDouble3(0);
  TypeSimple=true;
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JMotionPos::Reset(){
  PosSimple=TDouble3(0);
  PosMatrix.SetIdentity();
  TypeSimple=true;
}

//==============================================================================
// Aplica movimiento lineal.
//==============================================================================
void JMotionPos::Move(const tdouble3 &dis){
  if(TypeSimple)PosSimple=TDouble3(PosSimple.x+dis.x,PosSimple.y+dis.y,PosSimple.z+dis.z);
  else PosMatrix.Mul(JMatrix4d::MatrixMov(dis));
/*
  else{
//  printf("\nMove:(%g,%g,%g)\n",dis.x,dis.y,dis.z);
//    PosMatrix.Print("PreMove");
    PosMatrix.Mul(JMatrix4::MatrixMov(dis));
//    PosMatrix.Print("PostMove");
  }*/
}

//==============================================================================
// Aplica rotacion.
//==============================================================================
void JMotionPos::Rotate(double ang,const tdouble3 &axisp1,const tdouble3 &axisp2){
  if(TypeSimple)ToMatrix();
//  JMatrix4d m=JMatrix4d::MatrixRot(ang,axisp1,axisp2);
//  m.Print("m");
  PosMatrix.Mul(JMatrix4d::MatrixRot(ang,axisp1,axisp2));
/*
  JMatrix4d m=JMatrix4d::MatrixRot(ang,axisp1,axisp2);
  PosMatrix=m;
  TypeSimple=false;
/*
  if(TypeSimple&&modpos.TypeSimple){
    PosSimple=TDouble3(PosSimple.x+modpos.PosSimple.x,PosSimple.y+modpos.PosSimple.y,PosSimple.z+modpos.PosSimple.z);
  }
  else{
    //---->PDTE
  }
*/
}

//==============================================================================
// Aplica movimiento.
//==============================================================================
void JMotionPos::MoveMix(const JMotionPos &modpos){
  if(TypeSimple&&!modpos.TypeSimple)ToMatrix();
  if(modpos.TypeSimple)Move(modpos.PosSimple);
  else PosMatrix.Mul(modpos.PosMatrix);//<-Usando MulPre() da problemas...
}

//==============================================================================
// Convierte objeto a tipo matriz.
//==============================================================================
void JMotionPos::ToMatrix(){
  if(TypeSimple){
    PosMatrix=JMatrix4d::MatrixMov(PosSimple);
//  printf("PosSimple:(%g,%g,%g)\n",PosSimple.x,PosSimple.y,PosSimple.z);
//  PosMatrix.Print("\nToMatrix");
    TypeSimple=false;
  }
}

//==============================================================================
// Devuelve punto modificado al aplicarle el desplazamiento de PosSimple/PosMatrix
//==============================================================================
tdouble3 JMotionPos::PointMove(const tdouble3 &p) const{
  return(TypeSimple? TDouble3(p.x+PosSimple.x,p.y+PosSimple.y,p.z+PosSimple.z): PosMatrix.MulPoint(p));
}

