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

/// \file JDsOutputParts.cpp \brief Implements the class \ref JDsOutputParts.

#include "JDsOutputParts.h"
#include "Functions.h"
#include "JXml.h"
#include "FunGeo3d.h"
#include "JSpVtkShape.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "JSphMk.h"
#include <cstring>
#include <cstdlib>
#include <cfloat>

#ifdef _WITHGPU
#include "JSphGpu_ker.h"
#include "FunctionsBasic_iker.h"
#endif

using namespace std;

//##############################################################################
//# JDsOutputPartsOp_Init
//##############################################################################
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Init::GetConfig(std::vector<std::string>& lines)const{
  lines.push_back(fun::PrintStr("  Initial selection: %s (lv: %d)",(SelAll? "ALL": "NONE"),BitResult));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Init::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    const byte psel=sel[p];
    sel[p]=(SelAll? psel|resmask: psel&(~resmask));
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_Init::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  const byte resmask=(1<<BitResult);
  cusph::ComputeOutputPartsInit(resmask,SelAll,np,0,sel);
}
#endif


//##############################################################################
//# JDsOutputPartsOp_GroupFin
//##############################################################################
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_GroupFin::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode: %s, lv: %d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_GroupFin::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const byte resprev=(1<<BitResultPrev);
  const int n=int(np);
  //printf("---> inv:%d  cmband:%d",(Inverse?1:0),(CombineAnd?1:0));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    byte psel=sel[p];
    bool r=((psel&resprev)!=0);  //-lv=0
    bool ok=((psel&resmask)!=0); //-lv+1
    if(Inverse)ok=!ok;
    r=(CombineAnd? r&&ok: r||ok);
    psel=psel&(~resprev);
    if(r)psel=psel|resprev;
    sel[p]=psel;
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_GroupFin::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  const byte resmask=(1<<BitResult);
  const byte resprev=(1<<BitResultPrev);
  cusph::ComputeOutputPartsGroup(resmask,resprev,CombineAnd,Inverse,np,0,sel);
}
#endif


//##############################################################################
//# JDsOutputPartsOp_Pos
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputPartsOp_Pos::Reset(){
  PosMinInit=PosMin=TDouble3(-DBL_MAX); 
  PosMaxInit=PosMax=TDouble3(DBL_MAX);
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Pos::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  tdouble3 posmin=TDouble3(-DBL_MAX);
  tdouble3 posmax=TDouble3(DBL_MAX);
  //-Minimum position.
  {
    const TiXmlElement* elepos=sxml->GetFirstElement(ele,"posmin",true); 
    if(elepos){
      posmin.x=sxml->GetAttributeDouble(elepos,"x",true,-DBL_MAX);
      posmin.y=sxml->GetAttributeDouble(elepos,"y",true,-DBL_MAX);
      posmin.z=sxml->GetAttributeDouble(elepos,"z",true,-DBL_MAX);
    }
  }
  //-Maximum position.
  {
    const TiXmlElement* elepos=sxml->GetFirstElement(ele,"posmax",true); 
    if(elepos){
      posmax.x=sxml->GetAttributeDouble(elepos,"x",true,DBL_MAX);
      posmax.y=sxml->GetAttributeDouble(elepos,"y",true,DBL_MAX);
      posmax.z=sxml->GetAttributeDouble(elepos,"z",true,DBL_MAX);
    }
  }
  PosMinInit=PosMin=posmin;
  PosMaxInit=PosMax=posmax;
}
//==============================================================================
/// Update position of filter according to floating center.
//==============================================================================
void JDsOutputPartsOp_Pos::UpdateFtPos(const tdouble3& ftcenter){
  const tdouble3 mv=(ftcenter-FtFollowCen0);
  if(PosMin.x!=-DBL_MAX)PosMin.x=PosMinInit.x+mv.x;
  if(PosMin.y!=-DBL_MAX)PosMin.y=PosMinInit.y+mv.y;
  if(PosMin.z!=-DBL_MAX)PosMin.z=PosMinInit.z+mv.z;
  if(PosMax.x!= DBL_MAX)PosMax.x=PosMaxInit.x+mv.x;
  if(PosMax.y!= DBL_MAX)PosMax.y=PosMaxInit.y+mv.y;
  if(PosMax.z!= DBL_MAX)PosMax.z=PosMaxInit.z+mv.z;
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Pos::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode: %s, lv: %d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    LimitPoints: %s",fun::Double3xRangeStr(PosMin,PosMax).c_str()));
  if(FtFollowId!=UINT_MAX)lines.push_back(fun::PrintStr("    Following floating mkbound: %d",FtFollowMkb));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Pos::SaveVtkConfig(double size,JSpVtkShape* ss)const{
  const double size05=size/2;
  const bool defx=(PosMin.x!=-DBL_MAX && PosMax.x!=DBL_MAX);
  const bool defy=(PosMin.y!=-DBL_MAX && PosMax.y!=DBL_MAX);
  const bool defz=(PosMin.z!=-DBL_MAX && PosMax.z!=DBL_MAX);
  const bool undefx=(PosMin.x==-DBL_MAX && PosMax.x==DBL_MAX);
  const bool undefy=(PosMin.y==-DBL_MAX && PosMax.y==DBL_MAX);
  const bool undefz=(PosMin.z==-DBL_MAX && PosMax.z==DBL_MAX);
  const word id=word(Id);
  if(defx && defy && defz){
    tdouble3 p1=PosMin,p2=PosMax;
    ss->AddBoxSize(p1,p2-p1,id);
  }
  else if(defx && undefy && defz){
    tdouble3 p1=PosMin,p2=PosMax;
    p1.y=-size05; p2.y=size05;
    ss->AddBoxSize(p1,p2-p1,id);
  }
  else if(defx && defy && undefz){
    tdouble3 p1=PosMin,p2=PosMax;
    p1.z=-size05; p2.z=size05;
    ss->AddBoxSize(p1,p2-p1,id);
  }
  else if(undefx && defy && defz){
    tdouble3 p1=PosMin,p2=PosMax;
    p1.x=-size05; p2.x=size05;
    ss->AddBoxSize(p1,p2-p1,id);
  }
  else{
    tdouble3 pmin=PosMin,pmax=PosMax;
    if(!defx){
      if(pmin.x!=-DBL_MAX)pmax.x=pmin.x+size;
      else if(pmax.x!=DBL_MAX)pmin.x=pmax.x-size;
      else{ pmin.x=-size05; pmax.x=size05; }
    }
    if(!defy){
      if(pmin.y!=-DBL_MAX)pmax.y=pmin.y+size;
      else if(pmax.y!=DBL_MAX)pmin.y=pmax.y-size;
      else{ pmin.y=-size05; pmax.y=size05; }
    }
    if(!defz){
      if(pmin.z!=-DBL_MAX)pmax.z=pmin.z+size;
      else if(pmax.z!=DBL_MAX)pmin.z=pmax.z-size;
      else{ pmin.z=-size05; pmax.z=size05; }
    }
    const tdouble3 pm=(pmin+pmax)/2;
    if(PosMin.x!=-DBL_MAX){
      const tdouble3 p1=pmin;
      const tdouble3 p2=TDouble3(pmin.x,pmin.y,pmax.z);
      const tdouble3 p3=TDouble3(pmin.x,pmax.y,pmax.z);
      const tdouble3 p4=TDouble3(pmin.x,pmax.y,pmin.z);
      const tdouble3 pm1=TDouble3(pmin.x,pm.y,pm.z);
      const tdouble3 pm2=TDouble3(pmin.x+size05,pm.y,pm.z);
      ss->AddQuad(p1,p2,p3,p4,id);
      ss->AddQuadWire(p1,p2,p3,p4,id);
      ss->AddLine(pm1,pm2,id);
    }
    if(PosMin.y!=-DBL_MAX){
      const tdouble3 p1=pmin;
      const tdouble3 p2=TDouble3(pmin.x,pmin.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmin.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmin.y,pmin.z);
      const tdouble3 pm1=TDouble3(pm.x,pmin.y,pm.z);
      const tdouble3 pm2=TDouble3(pm.x,pmin.y+size05,pm.z);
      ss->AddQuad(p1,p2,p3,p4,id);
      ss->AddQuadWire(p1,p2,p3,p4,id);
      ss->AddLine(pm1,pm2,id);
    }
    if(PosMin.z!=-DBL_MAX){
      const tdouble3 p1=pmin;
      const tdouble3 p2=TDouble3(pmin.x,pmax.y,pmin.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmin.z);
      const tdouble3 p4=TDouble3(pmax.x,pmin.y,pmin.z);
      const tdouble3 pm1=TDouble3(pm.x,pm.y,pmin.z);
      const tdouble3 pm2=TDouble3(pm.x,pm.y,pmin.z+size05);
      ss->AddQuad(p1,p2,p3,p4,id);
      ss->AddQuadWire(p1,p2,p3,p4,id);
      ss->AddLine(pm1,pm2,id);
    }
    if(PosMax.x!=DBL_MAX){
      const tdouble3 p1=TDouble3(pmax.x,pmin.y,pmin.z);
      const tdouble3 p2=TDouble3(pmax.x,pmin.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmax.y,pmin.z);
      const tdouble3 pm1=TDouble3(pmax.x,pm.y,pm.z);
      const tdouble3 pm2=TDouble3(pmax.x-size05,pm.y,pm.z);
      ss->AddQuad(p1,p2,p3,p4,id);
      ss->AddQuadWire(p1,p2,p3,p4,id);
      ss->AddLine(pm1,pm2,id);
    }
    if(PosMax.y!=DBL_MAX){
      const tdouble3 p1=TDouble3(pmin.x,pmax.y,pmin.z);
      const tdouble3 p2=TDouble3(pmin.x,pmax.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmax.y,pmin.z);
      const tdouble3 pm1=TDouble3(pm.x,pmax.y,pm.z);
      const tdouble3 pm2=TDouble3(pm.x,pmax.y-size05,pm.z);
      ss->AddQuad(p1,p2,p3,p4,id);
      ss->AddQuadWire(p1,p2,p3,p4,id);
      ss->AddLine(pm1,pm2,id);
    }
    if(PosMax.z!=DBL_MAX){
      const tdouble3 p1=TDouble3(pmin.x,pmin.y,pmax.z);
      const tdouble3 p2=TDouble3(pmin.x,pmax.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmin.y,pmax.z);
      const tdouble3 pm1=TDouble3(pm.x,pm.y,pmax.z);
      const tdouble3 pm2=TDouble3(pm.x,pm.y,pmax.z-size05);
      ss->AddQuad(p1,p2,p3,p4,id);
      ss->AddQuadWire(p1,p2,p3,p4,id);
      ss->AddLine(pm1,pm2,id);
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Pos::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==CombineAnd){ //if((r && CombineAnd) || (!r && !CombineAnd)){
      const tdouble3 ps=pos[p];
      bool ok=(PosMin.x<=ps.x && ps.x <=PosMax.x &&
               PosMin.z<=ps.z && ps.z <=PosMax.z &&
               PosMin.y<=ps.y && ps.y <=PosMax.y);
      if(Inverse)ok=!ok;
      r=(CombineAnd? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_Pos::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  const byte resmask=(1<<BitResult);
  const double3 pmin=Double3(PosMin);
  const double3 pmax=Double3(PosMax);
  cusph::ComputeOutputPartsPos(resmask,CombineAnd,Inverse,pmin,pmax
    ,np,0,posxy,posz,sel);
}
#endif


//##############################################################################
//# JDsOutputPartsOp_Plane
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputPartsOp_Plane::Reset(){
  PointInit=Point=Vector=TDouble3(0); 
  Distance=0;
  Plane=TPlane3d(0);
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Plane::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  PointInit=Point=sxml->ReadElementDouble3(ele,"point");
  Vector=sxml->ReadElementDouble3(ele,"vector");
  Distance=sxml->ReadElementDouble(ele,"distance","v",true,DBL_MAX);
  Plane=fgeo::PlanePtVec(Point,Vector);
}
//==============================================================================
/// Update position of filter according to floating center.
//==============================================================================
void JDsOutputPartsOp_Plane::UpdateFtPos(const tdouble3& ftcenter){
  const tdouble3 mv=(ftcenter-FtFollowCen0);
  Point=PointInit+mv;
  Plane=fgeo::PlanePtVec(Point,Vector);
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Plane::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode: %s, lv: %d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Point: (%g,%g,%g)  Vector: (%g,%g,%g)  Dist: %g"
    ,Point.x,Point.y,Point.z,Vector.x,Vector.y,Vector.z,Distance));
  if(FtFollowId!=UINT_MAX)lines.push_back(fun::PrintStr("    Following floating mkbound: %d",FtFollowMkb));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Plane::SaveVtkConfig(double size,JSpVtkShape* ss)const{
  const word id=word(Id);
  tdouble3 vec=fgeo::VecUnitary(Vector);
  ss->AddQuadOrtho(Point,vec,size,id);
  ss->AddQuadOrthoWire(Point,vec,size,id);
  const tdouble3 p2=Point+(vec*(Distance<DBL_MAX? Distance: size/2));
  ss->AddLine(Point,p2,id);
  if(Distance<DBL_MAX){
    ss->AddQuadOrtho(p2,vec*-1.,size,id);
    ss->AddQuadOrthoWire(p2,vec*-1.,size,id);
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Plane::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const bool cmband=CombineAnd;
  const bool inverse=Inverse;
  const double maxdist=Distance;
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==cmband){ //if((r && cmband) || (!r && !cmband)){
      const tdouble3 ps=pos[p];
      const double dist=fgeo::PlanePoint(Plane,ps);
      //bool ok=(dist>=0);
      bool ok=(dist>=0 && dist<=maxdist);
      if(inverse)ok=!ok;
      r=(cmband? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_Plane::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const double4 plane=Double4(Plane);
  const float maxdist=float(Distance);
  cusph::ComputeOutputPartsPlane(resmask,CombineAnd,Inverse,plane,maxdist
    ,np,0,posxy,posz,sel);
}
#endif


//##############################################################################
//# JDsOutputPartsOp_Sphere
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputPartsOp_Sphere::Reset(){
  CentreInit=Centre=TDouble3(0); 
  Radius=0;
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Sphere::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  if(FtFollowId==UINT_MAX)CentreInit=Centre=sxml->ReadElementDouble3(ele,"centre");
  else CentreInit=Centre=sxml->ReadElementDouble3(ele,"centre",true,FtFollowCen0);
  Radius=sxml->ReadElementFloat(ele,"radius","v");
}
//==============================================================================
/// Update position of filter according to floating center.
//==============================================================================
void JDsOutputPartsOp_Sphere::UpdateFtPos(const tdouble3& ftcenter){
  const tdouble3 mv=(ftcenter-FtFollowCen0);
  Centre=CentreInit+mv;
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Sphere::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode: %s, lv: %d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Centre: (%g,%g,%g)  Radius: %g",Centre.x,Centre.y,Centre.z,Radius));
  if(FtFollowId!=UINT_MAX)lines.push_back(fun::PrintStr("    Following floating mkbound: %d",FtFollowMkb));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Sphere::SaveVtkConfig(double size,JSpVtkShape* ss)const{
  ss->AddSphere(Centre,Radius,word(Id));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Sphere::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const float radius2=float(Radius*Radius);
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==CombineAnd){ //if((r && CombineAnd) || (!r && !CombineAnd)){
      const tdouble3 ps=pos[p];
      const float dx=float(Centre.x-ps.x);
      const float dy=float(Centre.y-ps.y);
      const float dz=float(Centre.z-ps.z);
      bool ok=(dx*dx+dy*dy+dz*dz <= radius2);
      if(Inverse)ok=!ok;
      r=(CombineAnd? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_Sphere::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  const byte resmask=(1<<BitResult);
  const double3 pcen=Double3(Centre);
  const float radius2=float(Radius*Radius);
  cusph::ComputeOutputPartsSphere(resmask,CombineAnd,Inverse,pcen,radius2
    ,np,0,posxy,posz,sel);
}
#endif


//##############################################################################
//# JDsOutputPartsOp_Cylinder
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputPartsOp_Cylinder::Reset(){
  Point1Init=Point2Init=TDouble3(0); 
  Point1=Point2=TDouble3(0); 
  Radius=0;
  Plane=TPlane3d(0);
  Distance=0;
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Cylinder::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  Point1Init=Point1=sxml->ReadElementDouble3(ele,"point1");
  Point2Init=Point2=sxml->ReadElementDouble3(ele,"point2");
  Radius=sxml->ReadElementDouble(ele,"radius","v",true,DBL_MAX);
  //Plane=fgeo::PlaneNormalized(fgeo::PlanePtVec(Point,Vector));
  const tdouble3 vector=Point2-Point1;
  Plane=fgeo::PlanePtVec(Point1,vector);
  Distance=fgeo::PointDist(vector);
}
//==============================================================================
/// Update position of filter according to floating center.
//==============================================================================
void JDsOutputPartsOp_Cylinder::UpdateFtPos(const tdouble3& ftcenter){
  const tdouble3 mv=(ftcenter-FtFollowCen0);
  Point1=Point1Init+mv;
  Point2=Point2Init+mv;
  const tdouble3 vector=Point2-Point1;
  Plane=fgeo::PlanePtVec(Point1,vector);
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Cylinder::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode: %s, lv: %d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Point1: (%g,%g,%g)  Point2: (%g,%g,%g)  Radius: %g"
    ,Point1.x,Point1.y,Point1.z,Point2.x,Point2.y,Point2.z,Radius));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Cylinder::SaveVtkConfig(double size,JSpVtkShape* ss)const{
  ss->AddCylinder(Point1,Point2,Radius,word(Id));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Cylinder::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const bool isvertical=(Point1.x==Point2.x && Point1.y==Point2.y);
  const byte resmask=(1<<BitResult);
  const float radius2=float(Radius*Radius);
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==CombineAnd){ //if((r && CombineAnd) || (!r && !CombineAnd)){
      const tdouble3 ps=pos[p];

      const double dist=fgeo::PlanePoint(Plane,ps);
      bool ok=(dist>=0 && dist<=Distance);
      if(ok && isvertical){
        const double dx=(ps.x-Point1.x);
        const double dy=(ps.y-Point1.y);
        ok=(dx*dx+dy*dy <= radius2);
      }
      if(ok && !isvertical){
        //fgeo::LinePointDist(ps,Point1,Point2)
        const double ar=fgeo::TriangleArea(ps,Point1,Point2);
        const double dis=(ar*2)/Distance;
        ok=(dis<=Radius);
      }
      if(Inverse)ok=!ok;
      r=(CombineAnd? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_Cylinder::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  const byte resmask=(1<<BitResult);
  const double3 pcen1=Double3(Point1);
  const double3 pcen2=Double3(Point2);
  const double4 plane=Double4(Plane);
  const float maxdist=float(Distance);
  const float radius=float(Radius);
  cusph::ComputeOutputPartsCylinder(resmask,CombineAnd,Inverse
    ,plane,maxdist,pcen1,pcen2,radius,np,0,posxy,posz,sel);
}
#endif


//##############################################################################
//# JDsOutputPartsOp_Type
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputPartsOp_Type::Reset(){
  Types=0;
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Type::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  string txtypes=fun::StrLower(sxml->ReadElementStr(ele,"type","v"));
  vector<string> vtypes;
  unsigned n=fun::VectorSplitStr(",",txtypes,vtypes);
  for(unsigned c=0;c<n;c++){
    const string tx=fun::StrTrim(vtypes[c]);
    //printf("--->[%s]\n",tx.c_str());
    if(tx=="fixed")Types=(Types|1);
    else if(tx=="moving")Types=(Types|2);
    else if(tx=="floating")Types=(Types|4);
    else if(tx=="fluid")Types=(Types|8);
    else if(tx=="bound")Types=(Types|7);
    else sxml->ErrReadElement(ele,"type",false,"List of types is invalid.");
  }
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Type::GetConfig(std::vector<std::string>& lines)const{
  string txtypes;
  if(Types&7)txtypes=txtypes+(txtypes.empty()? "": ",")+"bound";
  else{
    if(Types&1)txtypes=txtypes+(txtypes.empty()? "": ",")+"fixed";
    if(Types&2)txtypes=txtypes+(txtypes.empty()? "": ",")+"moving";
    if(Types&4)txtypes=txtypes+(txtypes.empty()? "": ",")+"floating";
  }
  if(Types&8)txtypes=txtypes+(txtypes.empty()? "": ",")+"fluid";
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode: %s, lv: %d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Types: %s",txtypes.c_str()));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Type::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==CombineAnd){ //if((r && CombineAnd) || (!r && !CombineAnd)){
      const typecode codetp=CODE_GetType(code[p]);
      bool ok=((Types&1 && codetp==CODE_TYPE_FIXED) ||
               (Types&2 && codetp==CODE_TYPE_MOVING) ||
               (Types&4 && codetp==CODE_TYPE_FLOATING) ||
               (Types&8 && codetp==CODE_TYPE_FLUID) );
      if(Inverse)ok=!ok;
      r=(CombineAnd? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_Type::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  const byte resmask=(1<<BitResult);
  cusph::ComputeOutputPartsType(resmask,CombineAnd,Inverse,Types,np,0,code,sel);
}
#endif


//##############################################################################
//# JDsOutputPartsOp_Mk
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputPartsOp_Mk::Reset(){
  MkDef="";
  Mk1=Mk2=0;
  MkCode1=MkCode2=0;
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Mk::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  unsigned n=0;
  const bool usemk=sxml->ExistsElement(ele,"mk");
  const bool usemkb=sxml->ExistsElement(ele,"mkbound");
  const bool usemkf=sxml->ExistsElement(ele,"mkfluid");
  if(usemk)n++;
  if(usemkb)n++;
  if(usemkf)n++;
  if(n>1)Run_ExceptioonFile("Only one definition of mk is allowed (mk, mkbound or mkfluid).",sxml->ErrGetFileRow(ele));
  if(!n )Run_ExceptioonFile("No definition of mk, mkbound or mkfluid was found.",sxml->ErrGetFileRow(ele));
  string txmk;
  if(usemk )txmk=sxml->ReadElementStr(ele,"mk","v",true);
  if(usemkb)txmk=sxml->ReadElementStr(ele,"mkbound","v",true);
  if(usemkf)txmk=sxml->ReadElementStr(ele,"mkfluid","v",true);
  string txmk1,txmk2;
  if(fun::StrSplitCount("-",txmk)==1)txmk1=txmk2=txmk;
  else{
    txmk1=fun::StrSplitValue("-",txmk,0);
    txmk2=fun::StrSplitValue("-",txmk,1);
  }
  Mk1=atoi(txmk1.c_str());
  Mk2=atoi(txmk2.c_str());
  unsigned cbk1,cbk2;
  if(usemk ){
    MkDef=(Mk1!=Mk2? fun::PrintStr("mk: %d-%d",Mk1,Mk2): fun::PrintStr("mk: %d",Mk1));
    cbk1=MkInfo->GetMkBlockByMk(word(Mk1));
    cbk2=MkInfo->GetMkBlockByMk(word(Mk2));
    if(cbk1>=MkInfo->Size())Run_ExceptioonFile(fun::PrintStr("Mk value %d is invalid.",Mk1),sxml->ErrGetFileRow(ele));
    if(cbk2>=MkInfo->Size())Run_ExceptioonFile(fun::PrintStr("Mk value %d is invalid.",Mk2),sxml->ErrGetFileRow(ele));
  }
  if(usemkb){
    MkDef=(Mk1!=Mk2? fun::PrintStr("mkbound: %d-%d",Mk1,Mk2): fun::PrintStr("mkbound: %d",Mk1));
    cbk1=MkInfo->GetMkBlockByMkBound(word(Mk1));
    cbk2=MkInfo->GetMkBlockByMkBound(word(Mk2));
    if(cbk1>=MkInfo->Size())Run_ExceptioonFile(fun::PrintStr("Mkbound value %d is invalid.",Mk1),sxml->ErrGetFileRow(ele));
    if(cbk2>=MkInfo->Size())Run_ExceptioonFile(fun::PrintStr("Mkbound value %d is invalid.",Mk2),sxml->ErrGetFileRow(ele));
  }
  if(usemkf){
    MkDef=(Mk1!=Mk2? fun::PrintStr("mkfluid: %d-%d",Mk1,Mk2): fun::PrintStr("mkfluid: %d",Mk1));
    cbk1=MkInfo->GetMkBlockByMkFluid(word(Mk1));
    cbk2=MkInfo->GetMkBlockByMkFluid(word(Mk2));
    if(cbk1>=MkInfo->Size())Run_ExceptioonFile(fun::PrintStr("Mkfluid value %d is invalid.",Mk1),sxml->ErrGetFileRow(ele));
    if(cbk2>=MkInfo->Size())Run_ExceptioonFile(fun::PrintStr("Mkfluid value %d is invalid.",Mk2),sxml->ErrGetFileRow(ele));
  }
  if(MkInfo->Mkblock(cbk1)->Type!=MkInfo->Mkblock(cbk2)->Type)Run_ExceptioonFile("The values of the mk range include different types of particles.",sxml->ErrGetFileRow(ele));
  MkCode1=MkInfo->Mkblock(cbk1)->Code;
  MkCode2=MkInfo->Mkblock(cbk2)->Code;
  printf("--->  mk1[%s]:%d  mk2[%s]:%d\n",txmk1.c_str(),Mk1,txmk2.c_str(),Mk2);
  printf("--->  cbk1:%d  cbk2:%d  mksize:%d\n",cbk1,cbk2,MkInfo->Size());
  printf("--->  code1:%d  code2:%d\n",MkCode1,MkCode2);
  //Run_Exceptioon("aaaa");
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Mk::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode: %s, lv: %d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Selected %s",MkDef.c_str()));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Mk::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  const byte resmask=(1<<BitResult);
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<n;p++){
    byte psel=sel[p];
    bool r=((psel&resmask)!=0);
    if(r==CombineAnd){ //if((r && CombineAnd) || (!r && !CombineAnd)){
      const typecode codetp=CODE_GetTypeAndValue(code[p]);
      bool ok=(MkCode1<=codetp && codetp<=MkCode2);
      if(Inverse)ok=!ok;
      r=(CombineAnd? r&&ok: r||ok);
      psel=psel&(~resmask);
      if(r)psel=psel|resmask;
      sel[p]=psel;
    }
  }
}
#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputPartsOp_Mk::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  const byte resmask=(1<<BitResult);
  cusph::ComputeOutputPartsMk(resmask,CombineAnd,Inverse,MkCode1,MkCode2,np,0,code,sel);
}
#endif


//##############################################################################
//# JDsOutputParts
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsOutputParts::JDsOutputParts(bool cpu,double vtksize):Cpu(cpu),VtkSize(vtksize){
  ClassName="JDsOutputParts";
  Reset();
}
//==============================================================================
/// Destructor.
//==============================================================================
JDsOutputParts::~JDsOutputParts(){
  DestructorActive=true;
  Reset();
}
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputParts::Reset(){
  for(unsigned c=0;c<Count();c++)delete List[c];
  List.clear();
}
//==============================================================================
/// Reads configuration from the XML node.
//==============================================================================
void JDsOutputParts::ReadXml(const JXml* sxml,const TiXmlElement* ele,byte level
  ,const JSphMk* mkinfo,unsigned ftcount,const StFloatingData* ftobjs)
{
  const string presel=fun::StrUpper(sxml->GetAttributeStr(ele,"preselection",false));
  if(presel!="ALL" && presel!="NONE")sxml->ErrReadAtrib(ele,"preselection",false,"");
  byte result=level;
  JDsOutputPartsOp_Init* opeini=new JDsOutputPartsOp_Init(Count(),presel=="ALL",result,sxml,ele);
  List.push_back(opeini);
  const TiXmlElement* ele2=ele->FirstChildElement(); 
  while(ele2){
    const string cmd=ele2->Value();
    if(cmd.length() && cmd[0]!='_' && cmd!="ignore" && sxml->CheckElementActive(ele2)){
      //-Loads common configuration.
      const string oper=fun::StrUpper(sxml->GetAttributeStr(ele2,"operation",false));
      if(oper!="ADD" && oper!="DEL" && oper!="CONFIRM")sxml->ErrReadAtrib(ele2,"operation",false,"");
      const bool cmband=(oper=="CONFIRM" || oper=="DEL");
      bool inverse=sxml->GetAttributeBool(ele2,"inverse",true,false);
      if(oper=="DEL")inverse=!inverse;
      //-Loads following floating configuration.
      unsigned ftid=UINT_MAX;
      tdouble3 ftcen=TDouble3(0);
      unsigned ftmkb=UINT_MAX;
      const TiXmlElement* eleft=sxml->GetFirstElement(ele2,"ftfollow",true);
      if(eleft){
        ftmkb=sxml->ReadElementUnsigned(ele2,"ftfollow","mkbound",true,UINT_MAX);
        const unsigned ftmk =sxml->ReadElementUnsigned(ele2,"ftfollow","mk",true,UINT_MAX);
        if(ftmkb!=UINT_MAX && ftmk!=UINT_MAX)Run_ExceptioonFile("Only one definition of mk is allowed (mk or mkbound).",sxml->ErrGetFileRow(eleft));
        if(ftmkb==UINT_MAX && ftmk==UINT_MAX)ftmkb=sxml->ReadElementUnsigned(ele2,"ftfollow","mkbound",false);
        if(ftmkb==UINT_MAX)ftmkb=ftmk-unsigned(mkinfo->GetMkBoundFirst());
        for(unsigned cf=0;cf<ftcount && ftid==UINT_MAX;cf++)if(ftobjs[cf].mkbound==ftmkb)ftid=cf;
        if(ftid==UINT_MAX)Run_ExceptioonFile("No floating bodies with the specified mkbound (or mk).",sxml->ErrGetFileRow(eleft));
        ftcen=ftobjs[ftid].center;
      }
      //-Loads filters.
           if(cmd=="filterpos"     )List.push_back(new JDsOutputPartsOp_Pos     (Count(),inverse,cmband,result,ftid,ftmkb,ftcen,sxml,ele2));
      else if(cmd=="filterplane"   )List.push_back(new JDsOutputPartsOp_Plane   (Count(),inverse,cmband,result,ftid,ftmkb,ftcen,sxml,ele2));
      else if(cmd=="filtersphere"  )List.push_back(new JDsOutputPartsOp_Sphere  (Count(),inverse,cmband,result,ftid,ftmkb,ftcen,sxml,ele2));
      else if(cmd=="filtercylinder")List.push_back(new JDsOutputPartsOp_Cylinder(Count(),inverse,cmband,result,ftid,ftmkb,ftcen,sxml,ele2));
      else if(cmd=="filtertype"    )List.push_back(new JDsOutputPartsOp_Type    (Count(),inverse,cmband,result,sxml,ele2));
      else if(cmd=="filtermk"      )List.push_back(new JDsOutputPartsOp_Mk      (Count(),inverse,cmband,result,sxml,ele2,mkinfo));
      else if(cmd=="filtergroup"   ){
        if(level+1>7)Run_ExceptioonFile("No more than 7 nesting levels are allowed.",sxml->ErrGetFileRow(ele2));
        ReadXml(sxml,ele2,level+1,mkinfo,ftcount,ftobjs);
        JDsOutputPartsOp_GroupFin* gfin=new JDsOutputPartsOp_GroupFin(Count(),presel=="ALL",inverse,cmband,level+1,result,sxml,ele);
        List.push_back(gfin);
      }
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele2=ele2->NextSiblingElement();
  }
  if(level==0){
    Nparts=sxml->ReadElementInt(ele,"ignore","nparts",true,0);
    SaveVtkConfig(VtkSize);
  }
}
//==============================================================================
/// Loads configuration from XML object.
//==============================================================================
void JDsOutputParts::LoadXml(const JXml* sxml,const std::string& place
  ,const JSphMk* mkinfo,unsigned ftcount,const StFloatingData* ftobjs)
{
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place,false);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement(),0,mkinfo,ftcount,ftobjs);
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputParts::GetConfig(std::string txhead,std::string txfoot
  ,std::vector<std::string>& lines)const
{
  if(!txhead.empty())lines.push_back(txhead);
  if(Nparts>0)lines.push_back(fun::PrintStr("  Ignore every %d PARTs",Nparts));
  for(unsigned c=0;c<Count();c++)List[c]->GetConfig(lines);
  if(!txfoot.empty())lines.push_back(txfoot); 
}
//==============================================================================
/// Check PART number to apply filters.
//==============================================================================
bool JDsOutputParts::CheckFilters(int cpart)const{
  return(Nparts<=0 || (cpart%Nparts)!=0);
}

//==============================================================================
/// Update position of filters according to floating data.
//==============================================================================
void JDsOutputParts::UpdateFtPos(unsigned ftcount,const StFloatingData* ftobjs)
{
  const unsigned nf=Count();
  for(unsigned cf=0;cf<nf;cf++)if(List[cf]->FtFollowId!=UINT_MAX){
    const unsigned ftid=List[cf]->FtFollowId;
    List[cf]->UpdateFtPos(ftobjs[ftid].center);
  }
}

//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputParts::SaveVtkConfig(double size)const{
  if(Count()){
    JSpVtkShape ss;
    for(unsigned c=0;c<Count();c++)List[c]->SaveVtkConfig(size,&ss);
    string filevtk=AppInfo.GetDirOut()+"CfgOutputParts_Scheme.vtk";
    ss.SaveVtk(filevtk,(Count()>1? "num": ""));
    AppInfo.LogPtr()->AddFileInfo(filevtk,"Saves VTK file with OutputParts filters configurations.");
    ss.SaveVtk(AppInfo.GetDirOut()+"CfgOutputParts_Scheme_dg01.vtk",(Count()>1? "num": ""));
    ss.SaveVtk(AppInfo.GetDirOut()+"CfgOutputParts_Scheme_dg02.vtk",(Count()>1? "num": ""));
    ss.SaveVtk(AppInfo.GetDirOut()+"CfgOutputParts_Scheme_dg03.vtk",(Count()>1? "num": ""));
  }
}

//==============================================================================
/// Compute filters for particles.
//==============================================================================
void JDsOutputParts::ComputeFilterCpu(unsigned np,const unsigned* idp
  ,const tdouble3* pos,const typecode* code,byte* sel)const
{
  memset(sel,0,sizeof(byte)*np);
  for(unsigned cf=0;cf<Count();cf++)List[cf]->ComputeFilterCpu(np,idp,pos,code,sel);
}

#ifdef _WITHGPU
//==============================================================================
/// Compute filters for particles on GPU.
//==============================================================================
void JDsOutputParts::ComputeFilterGpu(unsigned np,const double2* posxy
  ,const double*posz,const typecode* code,byte*sel)const
{
  cudaMemset(sel,0,sizeof(byte)*np);
  for(unsigned cf=0;cf<Count();cf++)List[cf]->ComputeFilterGpu(np,posxy,posz,code,sel);
}
#endif

