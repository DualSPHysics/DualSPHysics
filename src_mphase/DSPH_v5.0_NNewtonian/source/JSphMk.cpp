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

/// \file JSphMk.cpp \brief Implements the class \ref JSphMk.

#include "JSphMk.h"
#include "Functions.h"
#include "JCaseParts.h"
#include "JPartDataHead.h"
#include <algorithm>
#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JSphMkBlock
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphMkBlock::JSphMkBlock(TpParticles type,unsigned mktype,unsigned mk,typecode code,unsigned begin,unsigned count)
  :Bound(IsBound(type)),Type(type),MkType(mktype),Mk(mk),Code(code),Begin(begin),Count(count)
{
  ClassName="JSphMkBlock";
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphMkBlock::Reset(){
  PosDefined=false;
  PosMin=PosMax=TDouble3(DBL_MAX);
}


//##############################################################################
//# JSphMk
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphMk::JSphMk(){
  ClassName="JSphMk";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphMk::~JSphMk(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphMk::Reset(){
  MkBoundFirst=MkFluidFirst=0;
  for(unsigned c=0;c<MkList.size();c++)delete MkList[c];
  MkList.clear();
  MkListSize=MkListFixed=MkListMoving=MkListFloat=MkListBound=MkListFluid=0;
  CodeNewFluid=0;
}

//==============================================================================
/// Load MK information of particles.
//==============================================================================
void JSphMk::Config(const JCaseParts *parts){
  MkBoundFirst=parts->GetMkBoundFirst();
  MkFluidFirst=parts->GetMkFluidFirst();
  MkListSize=parts->CountBlocks();
  MkListFixed=parts->CountBlocks(TpPartFixed);
  MkListMoving=parts->CountBlocks(TpPartMoving);
  MkListFloat=parts->CountBlocks(TpPartFloating);
  MkListFluid=parts->CountBlocks(TpPartFluid);
  MkListBound=MkListFixed+MkListMoving+MkListFloat;
  //:Log->Printf("__MkInfo> MkListSize  :%u",MkListSize);
  //:Log->Printf("__MkInfo> MkListFixed :%u",MkListFixed);
  //:Log->Printf("__MkInfo> MkListMoving:%u",MkListMoving);
  //:Log->Printf("__MkInfo> MkListFloat :%u",MkListFloat);
  //:Log->Printf("__MkInfo> MkListFluid :%u",MkListFluid);
  //:Log->Printf("__MkInfo> MkListBound :%u",MkListBound);
  //-Checks number of objects for each type.
  if(MkListFixed >CODE_MKRANGEMAX)Run_Exceptioon("The number of fixed particle blocks exceeds the maximum.");
  if(MkListMoving>CODE_MKRANGEMAX)Run_Exceptioon("The number of moving particle blocks exceeds the maximum.");
  if(MkListFloat >CODE_MKRANGEMAX)Run_Exceptioon("The number of floating particle blocks exceeds the maximum.");
  if(MkListFluid >CODE_MKRANGEMAX)Run_Exceptioon("The number of fluid particle blocks exceeds the maximum.");
  //-Gets info for each block of particles.
  for(unsigned c=0;c<MkListSize;c++){
    const JCasePartBlock &block=parts->GetBlock(c);
    const bool bound=(block.Bound);
    typecode code=0;
    switch(block.Type){
      case TpPartFixed:     code=CodeSetType(0,block.Type,c);                           break;
      case TpPartMoving:    code=CodeSetType(0,block.Type,c-MkListFixed);               break;
      case TpPartFloating:  code=CodeSetType(0,block.Type,c-MkListFixed-MkListMoving);  break;
      case TpPartFluid:     code=CodeSetType(0,block.Type,c-MkListBound);               break;
    }
    JSphMkBlock* pmk=new JSphMkBlock(block.Type,block.GetMkType(),block.GetMk(),code,block.GetBegin(),block.GetCount());
    MkList.push_back(pmk);
  }
  //-Checks number of fluid blocks.
  CodeNewFluid=CodeSetType(0,TpPartFluid,MkListFluid);
  if(CODE_GetTypeValue(CodeNewFluid)>=CODE_GetTypeValue(CODE_TYPE_FLUID_LIMITFREE))Run_Exceptioon("There are not free fluid codes for new particles created during the simulation.");
}

//==============================================================================
/// Returns number of blocks with a give type.
//==============================================================================
unsigned JSphMk::CountBlockType(TpParticles type)const{
  switch(type){
    case TpPartFixed:     return(MkListFixed);   break;
    case TpPartMoving:    return(MkListMoving);  break;
    case TpPartFloating:  return(MkListFloat);   break;
    case TpPartFluid:     return(MkListFluid);   break;
  }
  return(0);
}

//==============================================================================
/// Returns the first block in MkList according to a given type.
//==============================================================================
unsigned JSphMk::GetFirstBlockType(TpParticles type)const{
  switch(type){
    case TpPartFixed:     return(0);                          break;
    case TpPartMoving:    return(MkListFixed);                break;
    case TpPartFloating:  return(MkListFixed+MkListMoving);   break;
    case TpPartFluid:     return(MkListBound);                break;
  }
  return(MkListSize);
}

//==============================================================================
/// Returns the block in MkList according to a given Id.
//==============================================================================
unsigned JSphMk::GetMkBlockById(unsigned id)const{
  unsigned c=0;
  for(;c<MkListSize && id>=(MkList[c]->Begin+MkList[c]->Count);c++);
  return(c);
}

//==============================================================================
/// Returns the particle code according to a given Id.
//==============================================================================
typecode JSphMk::GetCodeById(unsigned id)const{
  const unsigned cmk=GetMkBlockById(id);
  if(cmk>=Size())Run_Exceptioon(fun::PrintStr("Mk block of particle (idp=%u) was not found.",id));
  return(MkList[cmk]->Code);
}

//==============================================================================
/// Returns the block in MkList according to a given MK.
//==============================================================================
unsigned JSphMk::GetMkBlockByMk(word mk)const{
  unsigned c=0;
  for(;c<MkListSize && unsigned(mk)!=MkList[c]->Mk;c++);
  return(c);
}

//==============================================================================
/// Returns the block in MkList according to a given bound MK.
//==============================================================================
unsigned JSphMk::GetMkBlockByMkBound(word mkbound)const{
  unsigned c=0;
  for(;c<MkListBound && unsigned(mkbound)!=MkList[c]->MkType;c++);
  return(c<MkListBound? c: MkListSize);
}

//==============================================================================
/// Returns the block in MkList according to a given fluid MK.
//==============================================================================
unsigned JSphMk::GetMkBlockByMkFluid(word mkfluid)const{
  unsigned c=MkListBound;
  for(;c<MkListSize && unsigned(mkfluid)!=MkList[c]->MkType;c++);
  return(c);
}

//==============================================================================
/// Returns the block in MkList according to a given Code.
/// Returns UINT_MAX if number of block is invalid.
//==============================================================================
unsigned JSphMk::GetMkBlockByCode(typecode code)const{
  unsigned ret=UINT_MAX;
  const unsigned type=CODE_GetType(code);
  const unsigned cblock=CODE_GetTypeValue(code);
  switch(type){
    case CODE_TYPE_FIXED:    if(cblock<MkListFixed) ret=cblock;                           break;
    case CODE_TYPE_MOVING:   if(cblock<MkListMoving)ret=cblock+MkListFixed;               break;
    case CODE_TYPE_FLOATING: if(cblock<MkListFloat) ret=cblock+MkListFixed+MkListMoving;  break;
    case CODE_TYPE_FLUID:    if(cblock<MkListFluid) ret=cblock+MkListBound;               break;
    default: Run_Exceptioon("Type of particle is invalid.");
  }
  return(ret);
}

//==============================================================================
/// Returns the code of a particle according to the given parameters.
//==============================================================================
typecode JSphMk::CodeSetType(typecode code,TpParticles type,unsigned value)const{
  //-Chooses type.
  typecode tp;
  if(type==TpPartFixed)tp=CODE_TYPE_FIXED;
  else if(type==TpPartMoving)tp=CODE_TYPE_MOVING;
  else if(type==TpPartFloating)tp=CODE_TYPE_FLOATING;
  else if(type==TpPartFluid)tp=CODE_TYPE_FLUID;
  else Run_Exceptioon("Type of particle is invalid.");
  //-Checks the value.
  typecode v=typecode(value&CODE_MASKVALUE);
  if(unsigned(v)!=value)Run_Exceptioon("The value is invalid.");
  //-Returns the new code.
  return(code&(~CODE_MASKTYPEVALUE)|tp|v);
}

////==============================================================================
///// Calculates domain limits for each Mk value.
////==============================================================================
//void JSphMk::ComputeMkDomains(bool bound,const std::vector<unsigned> &mklist,unsigned np,const tdouble3 *pos,const typecode *code){
//  for(unsigned c=0;c<unsigned(mklist.size());c++){
//    const unsigned cmk=(bound? GetMkBlockByMkBound(mklist[c]): GetMkBlockByMkFluid(mklist[c]));
//    if(cmk<Size() && !MkList[cmk]->GetPosDefined()){
//      const typecode rcode=MkList[cmk]->Code;
//      tdouble3 pmin=TDouble3(DBL_MAX),pmax=TDouble3(-DBL_MAX);
//      //-Calculates minimum and maximum position. 
//      for(unsigned p=0;p<np;p++)if(code[p]==rcode){
//        const tdouble3 ps=pos[p];
//        if(pmin.x>ps.x)pmin.x=ps.x;
//        if(pmin.y>ps.y)pmin.y=ps.y;
//        if(pmin.z>ps.z)pmin.z=ps.z;
//        if(pmax.x<ps.x)pmax.x=ps.x;
//        if(pmax.y<ps.y)pmax.y=ps.y;
//        if(pmax.z<ps.z)pmax.z=ps.z;
//      }
//      if(pmin<=pmax)MkList[cmk]->SetPosMinMax(pmin,pmax);
//      //if(pmin<=pmax)printf("------> mkfluid[%u]:%u pos:%s\n",c,mkfluidlist[c],fun::Double3gRangeStr(pmin,pmax).c_str());
//    }
//  }
//}

//==============================================================================
/// Calculates domain limits for each Mk value.
//==============================================================================
void JSphMk::ComputeMkDomains(unsigned np,const tdouble3 *pos,const typecode *code){
  //-Allocates memory and initializes dommain limits.
  tdouble3 *pmin=new tdouble3[MkListSize];
  tdouble3 *pmax=new tdouble3[MkListSize];
  for(unsigned c=0;c<MkListSize;c++){
    pmin[c]=TDouble3(DBL_MAX);
    pmax[c]=TDouble3(-DBL_MAX);
  }
  //-Calculates minimum and maximum position. 
  for(unsigned p=0;p<np;p++){
    const unsigned c=GetMkBlockByCode(code[p]);
    const tdouble3 ps=pos[p];
    if(pmin[c].x>ps.x)pmin[c].x=ps.x;
    if(pmin[c].y>ps.y)pmin[c].y=ps.y;
    if(pmin[c].z>ps.z)pmin[c].z=ps.z;
    if(pmax[c].x<ps.x)pmax[c].x=ps.x;
    if(pmax[c].y<ps.y)pmax[c].y=ps.y;
    if(pmax[c].z<ps.z)pmax[c].z=ps.z;
  }
  //-Configures minimum and maximum position for each MK. 
  for(unsigned c=0;c<MkListSize;c++)MkList[c]->SetPosMinMax(pmin[c],pmax[c]);
  //-Frees memory. 
  delete[] pmin;
  delete[] pmax;
}

//==============================================================================
/// Configures particle blocks in a JPartDataHead object.
//==============================================================================
void JSphMk::ConfigPartDataHead(JPartDataHead *parthead)const{
  for(unsigned c=0;c<Size();c++){
    const JSphMkBlock* pmbk=MkList[c];
    TpParticles type;
    switch(CODE_GetType(pmbk->Code)){
      case CODE_TYPE_FIXED:     type=TpPartFixed;     break;
      case CODE_TYPE_MOVING:    type=TpPartMoving;    break;
      case CODE_TYPE_FLOATING:  type=TpPartFloating;  break;
      case CODE_TYPE_FLUID:     type=TpPartFluid;     break;
      default: Run_Exceptioon("Type of particle block is invalid.");
    }
    parthead->ConfigParticles(type,pmbk->Mk,pmbk->MkType,pmbk->Begin,pmbk->Count);
  }
}

