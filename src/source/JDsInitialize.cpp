//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsInitialize.cpp \brief Implements the class \ref JDsInitialize.

#include "JDsInitialize.h"
#include "DualSphDef.h"
#include "JCaseProperties.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "FunGeo3d.h"
#include "JRangeFilter.h"
#include "JAppInfo.h"
#include "JXml.h"

#include <climits>
#include <cfloat>
#include <algorithm>

using namespace std;

//##############################################################################
//# JDsInitializeOp
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsInitializeOp::Reset(){
  OnlyPos=false;
  OnlyPosMin=OnlyPosMax=TDouble3(0);
  NpUpdated=NpTotal=0;
}

//==============================================================================
/// Reads onlypos filter information.
//==============================================================================
void JDsInitializeOp::ReadXmlOnlyPos(const JXml* sxml,TiXmlElement* ele){
  OnlyPos=false;
  OnlyPosMin=OnlyPosMax=TDouble3(0);
  ele=ele->FirstChildElement("onlypos"); 
  if(ele && sxml->CheckElementActive(ele)){
    OnlyPos=true;
    OnlyPosMin=TDouble3(-DBL_MAX);
    OnlyPosMax=TDouble3(DBL_MAX);
    //-Minimum position.
    TiXmlElement* elepos=sxml->GetFirstElement(ele,"posmin",true); 
    if(elepos){
      OnlyPosMin.x=sxml->GetAttributeDouble(elepos,"x",true,-DBL_MAX);
      OnlyPosMin.y=sxml->GetAttributeDouble(elepos,"y",true,-DBL_MAX);
      OnlyPosMin.z=sxml->GetAttributeDouble(elepos,"z",true,-DBL_MAX);
    }
    //-Maximum position.
    elepos=sxml->GetFirstElement(ele,"posmax",true); 
    if(elepos){
      OnlyPosMax.x=sxml->GetAttributeDouble(elepos,"x",true,DBL_MAX);
      OnlyPosMax.y=sxml->GetAttributeDouble(elepos,"y",true,DBL_MAX);
      OnlyPosMax.z=sxml->GetAttributeDouble(elepos,"z",true,DBL_MAX);
    }
  }
}

//==============================================================================
/// Calculates domain limits of MkType particles and returns number of particles.
//==============================================================================
unsigned JDsInitializeOp::ComputeDomainMk(bool bound,word mktp,unsigned np
  ,const word* mktype,const unsigned* idp,const tdouble3* pos
  ,tdouble3& posmin,tdouble3& posmax)const
{
  tdouble3 pmin=TDouble3(DBL_MAX),pmax=TDouble3(-DBL_MAX);
  unsigned n=0;
  for(unsigned p=0;p<np;p++)if(mktype[p]==mktp && bound==(idp[p]<InitCt.nbound)){
    const tdouble3 ps=pos[p];
    if(!OnlyPos || (OnlyPosMin<=ps && ps<=OnlyPosMax)){
      pmin=MinValues(pmin,ps);
      pmax=MaxValues(pmax,ps);
      n++;
    }
  }
  posmin=pmin;
  posmax=pmax;
  return(n);
}

//==============================================================================
/// Returns string with information about updated particles.
//==============================================================================
std::string JDsInitializeOp::GetConfigNp()const{
  string tx;
  if(NpUpdated!=NpTotal)tx=fun::PrintStr("(%u / %u particles)",NpUpdated,NpTotal);
  else                  tx=fun::PrintStr("(%u particles)",NpUpdated);
  return(tx);
}

//==============================================================================
/// Returns string with Mk configuration.
//==============================================================================
std::string JDsInitializeOp::GetConfigMkBound(std::string mktype)const{
  const string tmk=(mktype.empty()? "ALL": mktype.c_str());
  return(fun::PrintStr("  MkBound..: %s %s",tmk.c_str(),GetConfigNp().c_str()));
}

//==============================================================================
/// Returns string with Mk configuration.
//==============================================================================
std::string JDsInitializeOp::GetConfigMkFluid(std::string mktype)const{
  const string tmk=(mktype.empty()? "ALL": mktype.c_str());
  return(fun::PrintStr("  MkFluid..: %s %s",tmk.c_str(),GetConfigNp().c_str()));
}

//==============================================================================
/// Returns string with OnlyPos configuration.
//==============================================================================
std::string JDsInitializeOp::GetConfigOnlyPos()const{
  return(fun::PrintStr("  OnlyPos..: %s",fun::Double3xRangeStr(OnlyPosMin,OnlyPosMax,"%g").c_str()));
}


//##############################################################################
//# JDsInitializeOp_FluidVel
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsInitializeOp_FluidVel::Reset(){
  JDsInitializeOp::Reset();
  VelType=TVEL_Constant;
  MkFluid="";
  Direction=TFloat3(0);
  Vel1=Vel2=Vel3=0;
  Posz1=Posz2=Posz3=0;
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JDsInitializeOp_FluidVel::ReadXml(const JXml* sxml,TiXmlElement* xele){
  sxml->CheckElementNames(xele,true,"onlypos direction velocity velocity2 velocity3");
  ReadXmlOnlyPos(sxml,xele);
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
void JDsInitializeOp_FluidVel::Run(unsigned np,unsigned npb,const tdouble3* pos
  ,const unsigned* idp,const word* mktype,tfloat4* velrho,tfloat3* boundnor)
{
  const tfloat3 dir=fgeo::VecUnitary(Direction);
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
  else Run_Exceptioon("Velocity profile is unknown.");
  //-Updates selected particles.
  JRangeFilter rg(MkFluid);
  const bool all=(MkFluid.empty());
  for(unsigned p=npb;p<np;p++)if((all || rg.CheckValue(mktype[p])) && CheckPos(p,pos)){
    float v1=v;
    if(VelType==TVEL_Linear)v1=m2*float(pos[p].z)+b2;
    else if(VelType==TVEL_Parabolic)v1=a3*float(pos[p].z)*float(pos[p].z)+b3*float(pos[p].z)+c3;
    velrho[p].x=dir.x*v1;
    velrho[p].y=dir.y*v1;
    velrho[p].z=dir.z*v1;
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitializeOp_FluidVel::GetConfig(std::vector<std::string>& lines)const{
  lines.push_back(fun::PrintStr("  Operation: %s",ClassName.substr(BaseNameSize).c_str()));
  lines.push_back(GetConfigMkFluid(MkFluid));
  if(OnlyPos)lines.push_back(GetConfigOnlyPos());
  if(VelType==TVEL_Constant)lines.push_back(fun::PrintStr("  Constant velocity: %g",Vel1));
  else if(VelType==TVEL_Linear)lines.push_back(fun::PrintStr("  Linear velocity: %g(z=%g), %g(z=%g)",Vel1,Posz1,Vel2,Posz2));
  else if(VelType==TVEL_Parabolic)lines.push_back(fun::PrintStr("  Parabolic velocity: %g(z=%g), %g(z=%g), %g(z=%g)",Vel1,Posz1,Vel2,Posz2,Vel3,Posz3));
  else Run_Exceptioon("Velocity profile is unknown.");
}


//##############################################################################
//# JDsInitializeOp_BoundNormalSet
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsInitializeOp_BoundNormalSet::Reset(){
  JDsInitializeOp::Reset();
  MkBound=""; 
  Normal=TFloat3(0);
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JDsInitializeOp_BoundNormalSet::ReadXml(const JXml* sxml,TiXmlElement* xele){
  sxml->CheckElementNames(xele,true,"onlypos normal");
  ReadXmlOnlyPos(sxml,xele);
  MkBound=sxml->GetAttributeStr(xele,"mkbound",true);
  Normal=sxml->ReadElementFloat3(xele,"normal");
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalSet::Run(unsigned np,unsigned npb
  ,const tdouble3* pos,const unsigned* idp,const word* mktype,tfloat4* velrho
  ,tfloat3* boundnor)
{
  JRangeFilter rg(MkBound);
  const bool all=(MkBound.empty());
  for(unsigned p=0;p<np;p++)if(idp[p]<InitCt.nbound && (all || rg.CheckValue(mktype[p])) && CheckPos(p,pos)){
    boundnor[p]=Normal;
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalSet::GetConfig(std::vector<std::string>& lines)const{
  lines.push_back(fun::PrintStr("  Operation: %s",ClassName.substr(BaseNameSize).c_str()));
  lines.push_back(GetConfigMkBound(MkBound));
  if(OnlyPos)lines.push_back(GetConfigOnlyPos());
  lines.push_back(fun::PrintStr("  Normal...: (%g,%g,%g)",Normal.x,Normal.y,Normal.z));
}


//##############################################################################
//# JDsInitialize
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsInitialize::JDsInitialize(bool sim2d,double sim2dy,tdouble3 posmin,tdouble3 posmax
  ,double dp,float kernelh,const std::string& dirdatafile
  ,unsigned nbound,bool boundnormals):BoundNormals(boundnormals)
  ,InitCt(JDsInitializeOp::StrInitCt(sim2d,sim2dy,posmin,posmax,dp,kernelh,nbound,dirdatafile))
{
  ClassName="JDsInitialize";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsInitialize::~JDsInitialize(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsInitialize::Reset(){
  for(unsigned c=0;c<Count();c++)delete Opes[c];
  Opes.clear();
}

//==============================================================================
/// Loads data in XML format from a file.
//==============================================================================
void JDsInitialize::LoadFileXml(const std::string& file,const std::string& path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path);
}

//==============================================================================
/// Loads particles information from the object XML.
//==============================================================================
void JDsInitialize::LoadXml(const JXml* sxml,const std::string& place){
  //Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place);
  //if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads particles information in XML format.
//==============================================================================
void JDsInitialize::ReadXml(const JXml* sxml,TiXmlElement* lis){
  //-Loads elements.
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    string cmd=ele->Value();
    if(cmd.length() && cmd[0]!='_' && sxml->CheckElementActive(ele)){
      //printf("-----------> [%s]\n",cmd.c_str());
      if(cmd=="fluidvelocity"){ 
        JDsInitializeOp_FluidVel* ope=new JDsInitializeOp_FluidVel(sxml,ele,InitCt);
        Opes.push_back(ope); 
      }
      else if(cmd=="boundnormal_set"){
        if(BoundNormals){ 
          JDsInitializeOp_BoundNormalSet* ope=new JDsInitializeOp_BoundNormalSet(sxml,ele,InitCt);
          Opes.push_back(ope);
        }
      }
      else if(cmd=="boundnormal_plane" || cmd=="boundnormal_sphere" 
        || cmd=="boundnormal_cylinder" || cmd=="boundnormal_parts")Run_ExceptioonFile(fun::PrintStr(
          "The <initialize><%s> option is not supported by the current version. Use the <casedef><normals> section of the XML and GenCase v5.4.350 (or higher) for additional configuration options."
          ,cmd.c_str()),sxml->ErrGetFileRow(ele));
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JDsInitialize::Run(unsigned np,unsigned npb,const tdouble3* pos
  ,const unsigned* idp,const word* mktype,tfloat4* velrho,tfloat3* boundnor)
{
  for(unsigned c=0;c<Count();c++){
    Opes[c]->Run(np,npb,pos,idp,mktype,velrho,boundnor);
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitialize::GetConfig(std::vector<std::string>& lines)const{
  for(unsigned c=0;c<Count();c++){
    lines.push_back(fun::PrintStr("Initialize_%u",c));
    Opes[c]->GetConfig(lines);
  }
}
