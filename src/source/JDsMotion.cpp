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

/// \file JDsMotion.cpp \brief Implements the class \ref JDsMotion.

#include "JDsMotion.h"
#include "JCaseParts.h"
#include "JMotion.h"
#include "JXml.h"
#include <climits>

using namespace std;

//##############################################################################
//# JDsMotion
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsMotion::JDsMotion(bool simulate2d):Simulate2D(simulate2d){
  ClassName="JDsMotion";
  ObjMotion=NULL;
  Mot=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsMotion::~JDsMotion(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JDsMotion::Reset(){
  TimeMod=0;
  ObjCount=0;
  delete[] ObjMotion;  ObjMotion=NULL;
  memset(&MotionNull,0,sizeof(StMotionData));
  MotionNull.ref=USHRT_MAX; MotionNull.type=MOTT_None;
  delete Mot; Mot=NULL;
  ActiveMotion=false;
  LastDt=0;
}

//==============================================================================
/// Configures moving objects.
//==============================================================================
void JDsMotion::ConfigObjects(const JCaseParts *parts){
  ObjCount=parts->CountBlocks(TpPartMoving);
  if(ObjCount>CODE_MKRANGEMAX)Run_Exceptioon("The number of mobile objects exceeds the maximum.");
  //-Prepares memory.
  ObjMotion=new StMotionData[ObjCount];
  memset(ObjMotion,0,sizeof(StMotionData)*ObjCount);
  //-Loads configuration.
  unsigned cmot=0;
  for(unsigned c=0;c<parts->CountBlocks();c++){
    const JCasePartBlock &block=parts->GetBlock(c);
    if(block.Type==TpPartMoving){
      if(cmot>=ObjCount)Run_Exceptioon("The number of mobile objects exceeds the expected maximum.");
      //:printf("block[%2d]=%d -> %d\n",c,block.GetBegin(),block.GetCount());
      ObjMotion[cmot].ref=word(cmot);
      ObjMotion[cmot].mkbound=block.GetMkType();
      ObjMotion[cmot].idbegin=block.GetBegin();
      ObjMotion[cmot].count  =block.GetCount();
      ObjMotion[cmot].type   =MOTT_None;
      cmot++;
    }
  }
  if(cmot!=ObjCount)Run_Exceptioon("The number of mobile objects is invalid.");
}

//==============================================================================
/// Initialisation of configuration for moving objects.
//==============================================================================
void JDsMotion::Init(const JCaseParts *parts,JXml *jxml,const std::string &path
  ,const std::string &dirdata)
{
  Reset();
  //-Configures moving objects.
  ConfigObjects(parts);
  //-Configures predefined motion object.
  Mot=new JMotion();
  Mot->ReadXml(dirdata,jxml,path,false);
  Mot->Prepare();
  if(ObjCount!=unsigned(Mot->GetMaxRef()+1))
    Run_Exceptioon("The number of mobile objects do not match the predefined motions in XML file.");
}

//==============================================================================
/// Returns Idx of requested moving object accoring to MkBound.
/// Returns UINT_MAX when it was not found.
//==============================================================================
unsigned JDsMotion::GetObjIdxByMkBound(word mkbound)const{
  unsigned idx=0;
  for(;idx<ObjCount && ObjMotion[idx].mkbound!=mkbound;idx++);
  return(idx<ObjCount? idx: UINT_MAX);
}

//==============================================================================
/// Processes next time interval and returns true if there are active motions.
//==============================================================================
bool JDsMotion::ProcesTime(TpMotionMode mode,double timestep,double dt){
  ActiveMotion=false; LastDt=dt;
  if(mode==MOMT_Simple)ActiveMotion=Mot->ProcesTimeSimple(timestep+TimeMod,dt);
  if(mode==MOMT_Ace2dt)ActiveMotion=Mot->ProcesTimeAce(timestep+TimeMod,dt);
  //-Load motion data in ObjMotion[].
  if(true || ActiveMotion)for(unsigned ref=0;ref<ObjCount;ref++){
    StMotionData &m=ObjMotion[ref];
    bool typelinear;
    bool active=Mot->ProcesTimeGetData(ref,typelinear,m.linmov,m.linvel,m.linace,m.matmov,m.matmov2);
    if(active){
      m.type=(typelinear? MOTT_Linear: MOTT_Matrix);
      if(Simulate2D && typelinear==MOTT_Linear)m.linmov.y=m.linvel.y=m.linace.y=0;
    }
    else m.type=MOTT_None;
  }
  return(ActiveMotion);
}

//==============================================================================
/// Defines no motion for indicated object.
//==============================================================================
const StMotionData& JDsMotion::GetMotionData(unsigned idx)const{
  if(idx<GetNumObjects())return(ObjMotion[idx]);
  return(MotionNull);
}

//==============================================================================
/// Defines motion for indicated object.
//==============================================================================
void JDsMotion::SetMotionData(const StMotionData& d){
  if(d.ref>=GetNumObjects())Run_Exceptioon("Moving object does not exist.");
  StMotionData &m=ObjMotion[d.ref];
  m.type=d.type;
  if(m.type==MOTT_Linear){
    m.linmov=d.linmov;
    m.linvel=d.linvel;
  }
  if(m.type==MOTT_Matrix)m.matmov=d.matmov;
  ActiveMotion=(m.type!=MOTT_None);
}

//==============================================================================
/// Defines motion with acceleration for indicated object.
//==============================================================================
void JDsMotion::SetMotionDataAce(const StMotionData& d){
  SetMotionData(d);
  StMotionData &m=ObjMotion[d.ref];
  m.linace=d.linace;
  m.matmov2=d.matmov2;
}

//==============================================================================
/// Defines no motion for indicated object.
//==============================================================================
void JDsMotion::SetMotionDataNone(unsigned idx){
  if(idx<GetNumObjects())ObjMotion[idx].type=MOTT_None;
}

//==============================================================================
/// Defines linear motion for indicated object.
//==============================================================================
void JDsMotion::SetMotionDataLin(unsigned idx,const tdouble3 &linmov){
  if(idx<GetNumObjects()){
    //printf("JDsMotion::SetMotionData-> idx:%d tp:%d  motion:(%g,%g,%g)\n",idx,tpmov,simplemov.x,simplemov.y,simplemov.z);
    ObjMotion[idx].type=MOTT_Linear;
    ObjMotion[idx].linmov=linmov;
    ObjMotion[idx].linvel=(LastDt? linmov/LastDt: TDouble3(0));
    ObjMotion[idx].linace=TDouble3(0);
    ActiveMotion=true;
  }
}

//==============================================================================
/// Defines matrix motion for indicated object.
//==============================================================================
void JDsMotion::SetMotionDataMat(unsigned idx,const tmatrix4d &matmov){
  if(idx<GetNumObjects()){
    //printf("JDsMotion::SetMotionData-> idx:%d tp:%d  motion:(%g,%g,%g)\n",idx,tpmov,simplemov.x,simplemov.y,simplemov.z);
    ObjMotion[idx].type=MOTT_Matrix;
    ObjMotion[idx].matmov=matmov;
    ObjMotion[idx].matmov2=TMatrix4d(0);
    ActiveMotion=true;
  }
}


