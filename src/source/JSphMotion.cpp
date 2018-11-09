//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphMotion.cpp \brief Implements the class \ref JSphMotion.

#include "JSphMotion.h"
#include "JSpaceParts.h"
#include "JMotion.h"
#include "JXml.h"
#include <climits>

using namespace std;

//##############################################################################
//# JSphMotion
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphMotion::JSphMotion(){
  ClassName="JSphMotion";
  ObjBegin=NULL; ObjMkBound=NULL;
  Mot=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphMotion::~JSphMotion(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JSphMotion::Reset(){
  TimeMod=0;
  ObjCount=0;
  delete[] ObjBegin;   ObjBegin=NULL;
  delete[] ObjMkBound; ObjMkBound=NULL;
  delete Mot; Mot=NULL;
  ActiveMotion=false;
}

//==============================================================================
/// Configures moving objects.
//==============================================================================
void JSphMotion::ConfigObjects(const JSpaceParts *parts){
  const char met[]="ConfigObjects";
  ObjCount=parts->CountBlocks(TpPartMoving);
  if(ObjCount>CODE_MKRANGEMAX)RunException(met,"The number of mobile objects exceeds the maximum.");
  //-Prepares memory.
  ObjBegin=new unsigned[ObjCount+1];
  ObjMkBound=new word[ObjCount];
  memset(ObjBegin,0,sizeof(unsigned)*(ObjCount+1));
  memset(ObjMkBound,0,sizeof(word)*ObjCount);
  //-Loads configuration.
  unsigned cmot=0;
  for(unsigned c=0;c<parts->CountBlocks();c++){
    const JSpacePartBlock &block=parts->GetBlock(c);
    if(block.Type==TpPartMoving){
      if(cmot>=ObjCount)RunException(met,"The number of mobile objects exceeds the expected maximum.");
      //:printf("block[%2d]=%d -> %d\n",c,block.GetBegin(),block.GetCount());
      ObjBegin[cmot]=block.GetBegin();
      ObjBegin[cmot+1]=ObjBegin[cmot]+block.GetCount();
      ObjMkBound[cmot]=block.GetMkType();
      cmot++;
    }
  }
  if(cmot!=ObjCount)RunException(met,"The number of mobile objects is invalid.");
}

//==============================================================================
/// Initialisation of configuration for moving objects.
//==============================================================================
void JSphMotion::Init(const JSpaceParts *parts,JXml *jxml,const std::string &path,const std::string &dirdata){
  const char met[]="Init";
  Reset();
  //-Configures moving objects.
  ConfigObjects(parts);
  //-Configures predefined motion object.
  Mot=new JMotion();
  Mot->ReadXml(dirdata,jxml,path,false);
  Mot->Prepare();
  if(ObjCount!=unsigned(Mot->GetMaxRef()+1))RunException(met,"The number of mobile objects do not match the predefined motions in XML file.");
}

//==============================================================================
/// Returns MkBound of requested moving object.
//==============================================================================
word JSphMotion::GetObjMkBound(unsigned idx)const{
  if(idx>=ObjCount)RunException("GetObjMkBound","Moving object does not exist.");
  return(ObjMkBound[idx]);
}

//==============================================================================
/// Returns first particle of requested moving object.
//==============================================================================
unsigned JSphMotion::GetObjBegin(unsigned idx)const{
  if(idx>=ObjCount)RunException("GetObjBegin","Moving object does not exist.");
  return(ObjBegin[idx]);
}

//==============================================================================
/// Returns number of particles of requested moving object.
//==============================================================================
unsigned JSphMotion::GetObjSize(unsigned idx)const{
  if(idx>=ObjCount)RunException("GetObjSize","Moving object does not exist.");
  return(ObjBegin[idx+1]-ObjBegin[idx]);
}

//==============================================================================
/// Returns Idx of requested moving object accoring to MkBound.
/// Returns UINT_MAX when it was not found.
//==============================================================================
unsigned JSphMotion::GetObjIdxByMkBound(word mkbound)const{
  unsigned idx=0;
  for(;idx<GetNumObjects() && GetObjMkBound(idx)!=mkbound;idx++);
  return(idx<GetNumObjects()? idx: UINT_MAX);
}

//==============================================================================
/// Processes next time interval and returns true if there are active motions.
//==============================================================================
bool JSphMotion::ProcesTime(TpMotionMode mode,double timestep,double dt){
  ActiveMotion=false;
  if(mode==MOMT_Simple)ActiveMotion=Mot->ProcesTimeSimple(timestep+TimeMod,dt);
  if(mode==MOMT_Ace2dt)ActiveMotion=Mot->ProcesTimeAce(timestep+TimeMod,dt);
  return(ActiveMotion);
}

//==============================================================================
/// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JSphMotion::ProcesTimeGetData(unsigned ref,bool &typesimple,tdouble3 &simplemov
  ,tdouble3 &simplevel,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2
  ,unsigned &nparts,unsigned &idbegin)const
{
  const bool active=Mot->ProcesTimeGetData(ref,typesimple,simplemov,simplevel,simpleace,matmov,matmov2);
  if(active){
    idbegin=(ref<ObjCount? ObjBegin[ref]: 0);
    nparts=(ref<ObjCount? ObjBegin[ref+1]-idbegin: 0);
  }
  return(active);
}

//==============================================================================
/// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JSphMotion::ProcesTimeGetData(unsigned ref,word &mkbound
  ,bool &typesimple,tdouble3 &simplemov,tmatrix4d &matmov)const
{
  mkbound=GetObjMkBound(ref);
  const bool active=Mot->ProcesTimeGetData(ref,typesimple,simplemov,matmov);
  return(active);
}



