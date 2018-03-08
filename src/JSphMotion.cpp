//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

#include "JSphMotion.h"
#include "JMotion.h"
#include "JXml.h"

using namespace std;
//==============================================================================
/// Constructor.
//==============================================================================
JSphMotion::JSphMotion(){
  Mot=new JMotion();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphMotion::~JSphMotion(){
  delete Mot;     Mot=NULL;
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JSphMotion::Reset(){
  Mot->Reset();
}

//==============================================================================
/// Initialization of configuration and returns number of moving objects.
//==============================================================================
unsigned JSphMotion::Init(JXml *jxml,const std::string &path,const std::string &dirdata){
  Mot->ReadXml(dirdata,jxml,path,false);
  Mot->Prepare();
  int numobjects=Mot->GetMaxRef()+1;
  return(unsigned(numobjects));
}

//==============================================================================
/// Returns number of moving objects.
//==============================================================================
unsigned JSphMotion::GetNumObjects()const{
  const int numobjects=(Mot? Mot->GetMaxRef()+1: 0);
  return(unsigned(numobjects));
}

//==============================================================================
/// Processes next time interval and returns true if there are active motions.
//==============================================================================
bool JSphMotion::ProcesTime(TpMotionMode mode,double timestep,double dt){
  if(mode==MOMT_Simple)return(Mot->ProcesTimeSimple(timestep,dt));
  if(mode==MOMT_Ace2dt)return(Mot->ProcesTimeAce(timestep,dt));
  return(false);
}

//==============================================================================
/// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JSphMotion::ProcesTimeGetData(unsigned ref,bool &typesimple,tdouble3 &simplemov,tdouble3 &simplevel,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2)const{
  return(Mot->ProcesTimeGetData(ref,typesimple,simplemov,simplevel,simpleace,matmov,matmov2));
}

////==============================================================================
///// Processes next time interval and returns true if there are active motions.
////==============================================================================
//bool JSphMotion::ProcesTime(double timestep,double dt){
//  return(Mot->ProcesTime(timestep,dt));
//}

////==============================================================================
///// Reinicia movimientos en objeto Mot para iniciar el calculo desde el principio.
////==============================================================================
//void JSphMotion::ResetTime(double timestep){
//  Mot->ResetTime(timestep);
//}

////==============================================================================
///// Returns the number of performed movements.
////==============================================================================
//unsigned JSphMotion::GetMovCount()const{
//  return(Mot->GetMovCount());
//}

////==============================================================================
///// Returns data of the motion of an object.
////==============================================================================
//bool JSphMotion::GetMov(unsigned mov,unsigned &ref,tfloat3 &mvsimple,tmatrix4f &mvmatrix)const{
//  JMatrix4d aux;
//  tdouble3 mvsimpled;
//  bool ret=Mot->GetMov(mov,ref,mvsimpled,aux);
//  mvsimple=ToTFloat3(mvsimpled);
//  mvmatrix=aux.GetMatrix4f();
//  return(ret);
//}

////==============================================================================
///// Returns data of the motion of an object.
////==============================================================================
//bool JSphMotion::GetMov(unsigned mov,unsigned &ref,tdouble3 &mvsimple,tmatrix4d &mvmatrix)const{
//  JMatrix4d aux;
//  bool ret=Mot->GetMov(mov,ref,mvsimple,aux);
//  mvmatrix=aux.GetMatrix4d();
//  return(ret);
//}

////==============================================================================
///// Returns the number of finished movements.
////==============================================================================
//unsigned JSphMotion::GetStopCount()const{
//  return(Mot->GetStopCount());
//}

////==============================================================================
///// Returns the reference of the stopped object.
////==============================================================================
//unsigned JSphMotion::GetStopRef(unsigned mov)const{
//  return(Mot->GetStopRef(mov));
//}


