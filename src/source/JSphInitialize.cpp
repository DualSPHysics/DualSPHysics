//HEAD_DSPH
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

/// \file JSphInitialize.cpp \brief Implements the class \ref JSphInitialize.

#include "JSphInitialize.h"
#include "JSpaceProperties.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JRangeFilter.h"
#include "JXml.h"

using namespace std;

//##############################################################################
//# JSphInitializeOp_FluidVel
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInitializeOp_FluidVel::Reset(){
  VelType=TVEL_Constant;
  MkFluid="";
  Direction=TFloat3(0);
  Vel1=Vel2=Vel3=0;
  Posz1=Posz2=Posz3=0;
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JSphInitializeOp_FluidVel::ReadXml(JXml *sxml,TiXmlElement* xele){
  MkFluid=sxml->GetAttributeStr(xele,"mkfluid",true);
  Direction=sxml->ReadElementFloat3(xele,"direction");
  const byte vel1=(sxml->ExistsElement(xele,"velocity" )? 1: 0);
  const byte vel2=(sxml->ExistsElement(xele,"velocity2")? 1: 0);
  const byte vel3=(sxml->ExistsElement(xele,"velocity3")? 1: 0);
  if(vel1+vel2+vel3>1)sxml->ErrReadElement(xele,"velocity",false,"Several definitions for velocity were found.");
  if(vel1 || vel1+vel2+vel3==0){
    VelType=TVEL_Constant;
    Vel1=sxml->ReadElementFloat(xele,"velocity","v");
  }
  if(vel2){
    VelType=TVEL_Linear;
    Vel1 =sxml->ReadElementFloat(xele,"velocity2","v");
    Vel2 =sxml->ReadElementFloat(xele,"velocity2","v2");
    Posz1=sxml->ReadElementFloat(xele,"velocity2","z");
    Posz2=sxml->ReadElementFloat(xele,"velocity2","z2");
  }
  if(vel3){
    VelType=TVEL_Parabolic;
    Vel1 =sxml->ReadElementFloat(xele,"velocity3","v");
    Vel2 =sxml->ReadElementFloat(xele,"velocity3","v2");
    Vel3 =sxml->ReadElementFloat(xele,"velocity3","v3");
    Posz1=sxml->ReadElementFloat(xele,"velocity3","z");
    Posz2=sxml->ReadElementFloat(xele,"velocity3","z2");
    Posz3=sxml->ReadElementFloat(xele,"velocity3","z3");
  }
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JSphInitializeOp_FluidVel::Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop){
  const char met[]="Run";
  const tfloat3 dir=fmath::VecUnitary(Direction);
  float m2=0,b2=0;
  float a3=0,b3=0,c3=3;
  float v=0;
  if(VelType==TVEL_Constant)v=Vel1;
  else if(VelType==TVEL_Linear){
    m2=(Vel2-Vel1)/(Posz2-Posz1);
    b2=Vel1-m2*Posz1;
    //const float v=m2*float(pos[p].z)+b2;
  }
  else if(VelType==TVEL_Parabolic){
    const tmatrix3f inv=fmath::InverseMatrix3x3(TMatrix3f(Posz1*Posz1,Posz1,1,Posz2*Posz2,Posz2,1,Posz3*Posz3,Posz3,1));
    a3=inv.a11*Vel1+inv.a12*Vel2+inv.a13*Vel3;
    b3=inv.a21*Vel1+inv.a22*Vel2+inv.a23*Vel3;
    c3=inv.a31*Vel1+inv.a32*Vel2+inv.a33*Vel3;
    //const float v=a3*float(pos[p].z)*float(pos[p].z)+b3*float(pos[p].z)+c3;
  }
  else RunException(met,"Velocity profile is unknown.");
  JRangeFilter rg(MkFluid);
  bool all=(MkFluid.empty());
  for(unsigned p=npb;p<np;p++)if(all||rg.CheckValue(mktype[p])){
    float v1=v;
    if(VelType==TVEL_Linear)v1=m2*float(pos[p].z)+b2;
    else if(VelType==TVEL_Parabolic)v1=a3*float(pos[p].z)*float(pos[p].z)+b3*float(pos[p].z)+c3;
    velrhop[p].x=dir.x*v1;
    velrhop[p].y=dir.y*v1;
    velrhop[p].z=dir.z*v1;
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JSphInitializeOp_FluidVel::GetConfig(std::vector<std::string> &lines)const{
  const char met[]="GetConfig";
  lines.push_back(fun::PrintStr("  Operation: %s",ClassName.substr(17).c_str()));
  lines.push_back(fun::PrintStr("  MkFluid: %s",(MkFluid.empty()? "ALL": MkFluid.c_str())));
  if(VelType==TVEL_Constant)lines.push_back(fun::PrintStr("  Constant velocity: %g",Vel1));
  else if(VelType==TVEL_Linear)lines.push_back(fun::PrintStr("  Linear velocity: %g(z=%g), %g(z=%g)",Vel1,Posz1,Vel2,Posz2));
  else if(VelType==TVEL_Parabolic)lines.push_back(fun::PrintStr("  Parabolic velocity: %g(z=%g), %g(z=%g), %g(z=%g)",Vel1,Posz1,Vel2,Posz2,Vel3,Posz3));
  else RunException(met,"  Velocity profile is unknown.");
}


//##############################################################################
//# JSphInitialize
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInitialize::JSphInitialize(const std::string &file){
  ClassName="JSphInitialize";
  Reset();
  LoadFileXml(file,"case.execution.special.initialize");
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInitialize::~JSphInitialize(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInitialize::Reset(){
  for(unsigned c=0;c<Count();c++)delete Opes[c];
  Opes.clear();
}

//==============================================================================
/// Loads data in XML format from a file.
//==============================================================================
void JSphInitialize::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path);
}

//==============================================================================
/// Loads particles information from the object XML.
//==============================================================================
void JSphInitialize::LoadXml(JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNode(place,false);
  //if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  if(node)ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads particles information in XML format.
//==============================================================================
void JSphInitialize::ReadXml(JXml *sxml,TiXmlElement* lis){
  const char met[]="ReadXml";
  //-Loads fluidvelocity elements.
  TiXmlElement* ele=lis->FirstChildElement("fluidvelocity"); 
  while(ele){
    string cmd=ele->Value();
    if(cmd.length() && cmd[0]!='_'){
      if(cmd=="fluidvelocity"){
        JSphInitializeOp_FluidVel *ope=new JSphInitializeOp_FluidVel(sxml,ele);
        Opes.push_back(ope);
      }
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JSphInitialize::Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop){
  for(unsigned c=0;c<Count();c++){
    Opes[c]->Run(np,npb,pos,idp,mktype,velrhop);
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JSphInitialize::GetConfig(std::vector<std::string> &lines)const{
  for(unsigned c=0;c<Count();c++){
    lines.push_back(fun::PrintStr("Initialize_%u",c));
    Opes[c]->GetConfig(lines);
  }
}
