//HEAD_DSPH
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

/// \file JDsInitialize.cpp \brief Implements the class \ref JDsInitialize.

#include "JDsInitialize.h"
#include "JCaseProperties.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "FunGeo3d.h"
#include "JRangeFilter.h"
#include "JAppInfo.h"
#include "JXml.h"

#include <climits>
#include <cfloat>

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
void JDsInitializeOp::ReadXmlOnlyPos(const JXml *sxml,TiXmlElement* ele){
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
  ,const word *mktype,const unsigned *idp,const tdouble3 *pos
  ,tdouble3 &posmin,tdouble3 &posmax)const
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
void JDsInitializeOp_FluidVel::ReadXml(const JXml *sxml,TiXmlElement* xele){
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
void JDsInitializeOp_FluidVel::Run(unsigned np,unsigned npb,const tdouble3 *pos
  ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)
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
    velrhop[p].x=dir.x*v1;
    velrhop[p].y=dir.y*v1;
    velrhop[p].z=dir.z*v1;
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitializeOp_FluidVel::GetConfig(std::vector<std::string> &lines)const{
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
void JDsInitializeOp_BoundNormalSet::ReadXml(const JXml *sxml,TiXmlElement* xele){
  sxml->CheckElementNames(xele,true,"onlypos normal");
  ReadXmlOnlyPos(sxml,xele);
  MkBound=sxml->GetAttributeStr(xele,"mkbound",true);
  Normal=sxml->ReadElementFloat3(xele,"normal");
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalSet::Run(unsigned np,unsigned npb,const tdouble3 *pos
  ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)
{
  JRangeFilter rg(MkBound);
  const bool all=(MkBound.empty());
  for(unsigned p=0;p<np;p++)if(idp[p]<InitCt.nbound && (all || rg.CheckValue(mktype[p])) && CheckPos(p,pos)){
    boundnormal[p]=Normal;
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalSet::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back(fun::PrintStr("  Operation: %s",ClassName.substr(BaseNameSize).c_str()));
  lines.push_back(GetConfigMkBound(MkBound));
  if(OnlyPos)lines.push_back(GetConfigOnlyPos());
  lines.push_back(fun::PrintStr("  Normal...: (%g,%g,%g)",Normal.x,Normal.y,Normal.z));
}


//##############################################################################
//# JDsInitializeOp_BoundNormalPlane
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsInitializeOp_BoundNormalPlane::Reset(){
  JDsInitializeOp::Reset();
  MkBound="";
  PointAuto=false; 
  LimitDist=0; 
  Point=Normal=TFloat3(0); 
  MaxDisteH=0;
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JDsInitializeOp_BoundNormalPlane::ReadXml(const JXml *sxml,TiXmlElement* xele){
  sxml->CheckElementNames(xele,true,"onlypos point limitdist normal maxdisth");
  ReadXmlOnlyPos(sxml,xele);
  MkBound=sxml->GetAttributeStr(xele,"mkbound",true);
  PointAuto=sxml->ReadElementBool(xele,"point","auto",true,false);
  if(!PointAuto)Point=sxml->ReadElementFloat3(xele,"point");
  LimitDist=sxml->ReadElementFloat(xele,"limitdist","vdp",true,0.5f);
  Normal=sxml->ReadElementFloat3(xele,"normal");
  MaxDisteH=sxml->ReadElementFloat(xele,"maxdisth","v",true,2.f);
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalPlane::Run(unsigned np,unsigned npb,const tdouble3 *pos
  ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)
{
  const double maxdist=(MaxDisteH>0? InitCt.kernelh*MaxDisteH: DBL_MAX);
  const double limitdis=InitCt.dp*LimitDist;
  const tdouble3 nor=fgeo::VecUnitary(ToTDouble3(Normal));
  //-Define direction of normal (undefined, top, bottom, left, right...).
  byte nordir=0;
       if(fun::IsEqual(nor,TDouble3( 0, 0, 1),0.0000001))nordir=1;
  else if(fun::IsEqual(nor,TDouble3( 0, 0,-1),0.0000001))nordir=2;
  else if(fun::IsEqual(nor,TDouble3(-1, 0, 0),0.0000001))nordir=3;
  else if(fun::IsEqual(nor,TDouble3( 1, 0, 0),0.0000001))nordir=4;
  else if(fun::IsEqual(nor,TDouble3( 0,-1, 0),0.0000001))nordir=5;
  else if(fun::IsEqual(nor,TDouble3( 0, 1, 0),0.0000001))nordir=6;
  //-Processes particles.
  JRangeFilter rg(MkBound);
  unsigned v=rg.GetFirstValue();
  while(v!=UINT_MAX){
    const word mktp=word(v);
    tdouble3 point=ToTDouble3(Point);
    //-Calculates point when it is undefined.
    if(PointAuto){
      //-Calculates domain limits.
      tdouble3 pmin,pmax;
      const unsigned n=ComputeDomainMk(true,mktp,np,mktype,idp,pos,pmin,pmax);
      if(n){
        const tdouble3 pmed=(pmin+pmax)/2.;
        //-Calculate point in boundary limit plane.
        if(nordir){//-Defines point according to an orthogonal normal.
          point=pmed;
          switch(nordir){
            case 1:  point.z=pmax.z+limitdis;  break;
            case 2:  point.z=pmin.z-limitdis;  break;
            case 3:  point.x=pmin.x-limitdis;  break;
            case 4:  point.x=pmax.x+limitdis;  break;
            case 5:  point.y=pmin.y-limitdis;  break;
            case 6:  point.y=pmax.y+limitdis;  break;
          }
        }
        else{//-Defines point according to an arbitrary normal.
          const tplane3d pla=fgeo::PlanePtVec(pmed,nor);
          double dismax=-DBL_MAX;
          for(unsigned p=0;p<np;p++)if(mktype[p]==mktp && idp[p]<InitCt.nbound){
            const double dist=fgeo::PlaneDistSign(pla,pos[p]);
            if(dist>dismax)dismax=dist;
          }
          point=pmed+(nor*(dismax+limitdis));
        }
      }
    }
    //-Compute normals according to calculated plane.
    const tplane3d pla=fgeo::PlanePtVec(point,nor);
    for(unsigned p=0;p<np;p++)if(mktype[p]==mktp && idp[p]<InitCt.nbound && CheckPos(p,pos)){
      const tdouble3 ps=pos[p];
      const tdouble3 psb=fgeo::PlaneOrthogonalPoint(ps,pla);
      boundnormal[p]=(fgeo::PointsDist(ps,psb)<maxdist? ToTFloat3(psb-ps): TFloat3(0));
    }
    //-Next mktp value.
    v=rg.GetNextValue(v);
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalPlane::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back(fun::PrintStr("  Operation: %s",ClassName.substr(BaseNameSize).c_str()));
  lines.push_back(GetConfigMkBound(MkBound));
  if(OnlyPos)lines.push_back(GetConfigOnlyPos());
  if(PointAuto)      lines.push_back("  Point....: Automatic");
  else lines.push_back(fun::PrintStr("  Point....: (%g,%g,%g)",Point.x,Point.y,Point.z));
  lines.push_back(fun::PrintStr("  Normal...: (%g,%g,%g)",Normal.x,Normal.y,Normal.z));
  lines.push_back(fun::PrintStr("  MaxDistH.: %g",MaxDisteH));
}


//##############################################################################
//# JDsInitializeOp_BoundNormalSphere
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsInitializeOp_BoundNormalSphere::Reset(){
  JDsInitializeOp::Reset();
  MkBound=""; 
  Center=TFloat3(0); 
  MaxDisteH=Radius=0; 
  Inside=true;
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JDsInitializeOp_BoundNormalSphere::ReadXml(const JXml *sxml,TiXmlElement* xele){
  sxml->CheckElementNames(xele,true,"onlypos center radius inside maxdisth");
  ReadXmlOnlyPos(sxml,xele);
  MkBound=sxml->GetAttributeStr(xele,"mkbound",true);
  Center=sxml->ReadElementFloat3(xele,"center");
  Radius=sxml->ReadElementFloat(xele,"radius","v");
  Inside=sxml->ReadElementBool(xele,"inside","v");
  MaxDisteH=sxml->ReadElementFloat(xele,"maxdisth","v",true,2.f);
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalSphere::Run(unsigned np,unsigned npb,const tdouble3 *pos
  ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)
{
  const tdouble3 pcen=ToTDouble3(Center);
  const double ra=double(Radius);
  const double maxdist=(MaxDisteH>0? InitCt.kernelh*MaxDisteH: DBL_MAX);
  //-Processes particles.
  JRangeFilter rg(MkBound);
  const bool all=(MkBound.empty());
  for(unsigned p=0;p<np;p++)if(idp[p]<InitCt.nbound && (all || rg.CheckValue(mktype[p])) && CheckPos(p,pos)){
    const tdouble3 ps=pos[p];
    const double dissurf=ra-fgeo::PointsDist(pcen,ps); //-Distance to surface of sphere (inside:+dis, outside:-dis).
    const bool ps_in =(dissurf>0);//-Particle inside the sphere.
    const bool ps_out=(dissurf<0);//-Particle outside the sphere.
    tdouble3 psb=TDouble3(DBL_MAX);
    if((Inside  && ps_in  && dissurf<=maxdist) || (!Inside && ps_out && fabs(dissurf)<=maxdist)){
      psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra); //-Normal to surface limit.
    }
    boundnormal[p]=(psb.x!=DBL_MAX? ToTFloat3(psb-ps): TFloat3(0));
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalSphere::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back(fun::PrintStr("  Operation: %s",ClassName.substr(BaseNameSize).c_str()));
  lines.push_back(GetConfigMkBound(MkBound));
  if(OnlyPos)lines.push_back(GetConfigOnlyPos());
  lines.push_back(fun::PrintStr("  Center...: (%g,%g,%g)",Center.x,Center.y,Center.z));
  lines.push_back(fun::PrintStr("  Radius...: %g",Radius));
  lines.push_back(fun::PrintStr("  Inside...: %s",(Inside? "true": "false")));
  lines.push_back(fun::PrintStr("  MaxDistH.: %g",MaxDisteH));
}


//##############################################################################
//# JDsInitializeOp_BoundNormalCylinder
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsInitializeOp_BoundNormalCylinder::Reset(){
  JDsInitializeOp::Reset();
  MkBound="";
  Center1=Center2=TFloat3(0);
  MaxDisteH=Radius=0;
  Inside=true;
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JDsInitializeOp_BoundNormalCylinder::ReadXml(const JXml *sxml,TiXmlElement* xele){
  sxml->CheckElementNames(xele,true,"onlypos center1 center2 radius inside maxdisth");
  ReadXmlOnlyPos(sxml,xele);
  MkBound=sxml->GetAttributeStr(xele,"mkbound",true);
  Center1=sxml->ReadElementFloat3(xele,"center1");
  Center2=sxml->ReadElementFloat3(xele,"center2");
  Radius=sxml->ReadElementFloat(xele,"radius","v");
  Inside=sxml->ReadElementBool(xele,"inside","v");
  MaxDisteH=sxml->ReadElementFloat(xele,"maxdisth","v",true,2.f);
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalCylinder::Run(unsigned np,unsigned npb,const tdouble3 *pos
  ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)
{
  const tdouble3 cen1=ToTDouble3(Center1);
  const tdouble3 cen2=ToTDouble3(Center2);
  const double ra=double(Radius);
  const double maxdist=(MaxDisteH>0? InitCt.kernelh*MaxDisteH: DBL_MAX);
  const double tolerance=0.01*InitCt.kernelh;
  const tdouble3 vbot=fgeo::VecUnitary(cen1-cen2);
  const tdouble3 vtop=fgeo::VecUnitary(cen2-cen1);
  const tplane3d platop=fgeo::PlanePtVec(cen2,vbot);
  const tplane3d plabot=fgeo::PlanePtVec(cen1,vtop);
  //-Processes particles.
  JRangeFilter rg(MkBound);
  const bool all=(MkBound.empty());
  for(unsigned p=0;p<np;p++)if(idp[p]<InitCt.nbound && (all || rg.CheckValue(mktype[p])) && CheckPos(p,pos)){
    const tdouble3 ps=pos[p];
    const tdouble3 pcen=fgeo::LineOrthogonalPoint(pos[p],cen1,cen2);
    const double disbody=ra-fgeo::PointsDist(pcen,ps);    //-Distance to boundary limit of body.
    const double distop =fgeo::PlaneDistSign(platop,ps);  //-Distance to boundary limit of top.
    const double disbot =fgeo::PlaneDistSign(plabot,ps);  //-Distance to boundary limit of bottom.
    const bool ps_in =(disbody>0 && distop>0 && disbot>0);//-Particle inside the cylinder.
    const bool ps_out=(disbody<0 || distop<0 || disbot<0);//-Particle outside the cylinder.
    tdouble3 psb=TDouble3(DBL_MAX);
    byte sel=0;
    if(Inside && ps_in && (disbody<=maxdist || distop<=maxdist || disbot<=maxdist)){
      if(disbot<=distop){
        if(fun::IsEqual(disbody,disbot,tolerance))psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra)+vbot*disbot; //-Normal to bottom border.
        else if(disbody<disbot)psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra);//-Normal to body limit.
        else psb=ps+vbot*disbot; //-Normal to bottom limit.
      }
      else{
        if(fun::IsEqual(disbody,distop,tolerance))psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra)+vtop*distop; //-Normal to top border.
        else if(disbody<distop)psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra);//-Normal to body limit.
        else psb=ps+vtop*distop; //-Normal to top limit.
      } 
    }
    if(!Inside && ps_out && ((disbody<0 && fabs(disbody)<=maxdist) || (distop<0 && fabs(distop)<=maxdist) || (disbot<0 && fabs(disbot)<=maxdist))){
      const double adisbody=fabs(disbody),adistop=fabs(distop),adisbot=fabs(disbot);
      if(disbody<0){
             if(disbot<0)psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra)+vbot*disbot; //-Normal to body-bottom border.
        else if(distop<0)psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra)+vtop*distop; //-Normal to body-top border.
        else psb=pcen+(fgeo::VecUnitary(ps-pcen)*ra);//-Normal to body limit. 
      }
      else if(disbot<0)psb=ps+vbot*disbot; //-Normal to bottom limit.
      else if(distop<0)psb=ps+vtop*distop; //-Normal to top limit.
    }
    boundnormal[p]=(psb.x!=DBL_MAX? ToTFloat3(psb-ps): TFloat3(0));
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitializeOp_BoundNormalCylinder::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back(fun::PrintStr("  Operation: %s",ClassName.substr(BaseNameSize).c_str()));
  lines.push_back(GetConfigMkBound(MkBound));
  if(OnlyPos)lines.push_back(GetConfigOnlyPos());
  lines.push_back(fun::PrintStr("  Centers..: (%g,%g,%g)-(%g,%g,%g)",Center1.x,Center1.y,Center1.z,Center2.x,Center2.y,Center2.z));
  lines.push_back(fun::PrintStr("  Radius...: %g",Radius));
  lines.push_back(fun::PrintStr("  Inside...: %s",(Inside? "true": "false")));
  lines.push_back(fun::PrintStr("  MaxDistH.: %g",MaxDisteH));
}


//##############################################################################
//# JDsInitialize
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsInitialize::JDsInitialize(const JXml *sxml,const std::string &place
  ,const std::string &dirdatafile,float kernelh,float dp,unsigned nbound
  ,bool boundnormals):BoundNormals(boundnormals)
  ,InitCt(JDsInitializeOp::StrInitCt(kernelh,dp,nbound,dirdatafile))
{
  ClassName="JDsInitialize";
  Reset();
  LoadXml(sxml,place);
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
void JDsInitialize::LoadFileXml(const std::string &file,const std::string &path){
  JXml jxml;
  jxml.LoadFile(file);
  LoadXml(&jxml,path);
}

//==============================================================================
/// Loads particles information from the object XML.
//==============================================================================
void JDsInitialize::LoadXml(const JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place);
  //if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads particles information in XML format.
//==============================================================================
void JDsInitialize::ReadXml(const JXml *sxml,TiXmlElement* lis){
  //-Loads elements.
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    string cmd=ele->Value();
    if(cmd.length() && cmd[0]!='_' && sxml->CheckElementActive(ele)){
      //printf("-----------> [%s]\n",cmd.c_str());
      if(cmd=="fluidvelocity"){ 
        JDsInitializeOp_FluidVel *ope=new JDsInitializeOp_FluidVel(sxml,ele,InitCt);
        Opes.push_back(ope); 
      }
      else if(cmd=="boundnormal_set"     ){ if(BoundNormals){ JDsInitializeOp_BoundNormalSet      *ope=new JDsInitializeOp_BoundNormalSet     (sxml,ele,InitCt); Opes.push_back(ope); } }
      else if(cmd=="boundnormal_plane"   ){ if(BoundNormals){ JDsInitializeOp_BoundNormalPlane    *ope=new JDsInitializeOp_BoundNormalPlane   (sxml,ele,InitCt); Opes.push_back(ope); } }
      else if(cmd=="boundnormal_sphere"  ){ if(BoundNormals){ JDsInitializeOp_BoundNormalSphere   *ope=new JDsInitializeOp_BoundNormalSphere  (sxml,ele,InitCt); Opes.push_back(ope); } }
      else if(cmd=="boundnormal_cylinder"){ if(BoundNormals){ JDsInitializeOp_BoundNormalCylinder *ope=new JDsInitializeOp_BoundNormalCylinder(sxml,ele,InitCt); Opes.push_back(ope); } }
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
}

//==============================================================================
/// Initializes data of particles according XML configuration.
//==============================================================================
void JDsInitialize::Run(unsigned np,unsigned npb,const tdouble3 *pos
  ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)
{
  for(unsigned c=0;c<Count();c++){
    Opes[c]->Run(np,npb,pos,idp,mktype,velrhop,boundnormal);
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsInitialize::GetConfig(std::vector<std::string> &lines)const{
  for(unsigned c=0;c<Count();c++){
    lines.push_back(fun::PrintStr("Initialize_%u",c));
    Opes[c]->GetConfig(lines);
  }
}
