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

/// \file JDsOutputParts.cpp \brief Implements the class \ref JDsOutputParts.

#include "JDsOutputParts.h"
#include "Functions.h"
#include "JXml.h"
#include "FunGeo3d.h"
#include "JVtkLib.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include <cstring>
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
  lines.push_back(fun::PrintStr("  Initial selection: %s (lv:%d)",(SelAll? "ALL": "NONE"),BitResult));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Init::ComputeFilterCpu(unsigned np,const tdouble3* pos
  ,const typecode* code,byte* sel)const
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

}
#endif


//##############################################################################
//# JDsOutputPartsOp_Pos
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsOutputPartsOp_Pos::Reset(){
  PosMin=TDouble3(-DBL_MAX); 
  PosMax=TDouble3(DBL_MAX);
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
  PosMin=posmin;
  PosMax=posmax;
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Pos::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode:%s, lv:%d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    LimitPoints: %s",fun::Double3xRangeStr(PosMin,PosMax).c_str()));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Pos::SaveVtkConfig(double size,JVtkLib* sh)const{
  const double size05=size/2;
  const bool defx=(PosMin.x!=-DBL_MAX && PosMax.x!=DBL_MAX);
  const bool defy=(PosMin.y!=-DBL_MAX && PosMax.y!=DBL_MAX);
  const bool defz=(PosMin.z!=-DBL_MAX && PosMax.z!=DBL_MAX);
  const bool undefx=(PosMin.x==-DBL_MAX && PosMax.x==DBL_MAX);
  const bool undefy=(PosMin.y==-DBL_MAX && PosMax.y==DBL_MAX);
  const bool undefz=(PosMin.z==-DBL_MAX && PosMax.z==DBL_MAX);
  sh->SetShapeWireMode(false);
  if(defx && defy && defz){
    tdouble3 p1=PosMin,p2=PosMax;
    sh->AddShapeBoxSize(p1,p2-p1,int(Id));
  }
  else if(defx && undefy && defz){
    tdouble3 p1=PosMin,p2=PosMax;
    p1.y=-size05; p2.y=size05;
    sh->AddShapeBoxSize(p1,p2-p1,int(Id));
  }
  else if(defx && defy && undefz){
    tdouble3 p1=PosMin,p2=PosMax;
    p1.z=-size05; p2.z=size05;
    sh->AddShapeBoxSize(p1,p2-p1,int(Id));
  }
  else if(undefx && defy && defz){
    tdouble3 p1=PosMin,p2=PosMax;
    p1.x=-size05; p2.x=size05;
    sh->AddShapeBoxSize(p1,p2-p1,int(Id));
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
      sh->AddShapeQuad(p1,p2,p3,p4,int(Id));
      sh->AddShapeQuadWire(p1,p2,p3,p4,int(Id));
      sh->AddShapeLine(pm1,pm2,int(Id));
    }
    if(PosMin.y!=-DBL_MAX){
      const tdouble3 p1=pmin;
      const tdouble3 p2=TDouble3(pmin.x,pmin.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmin.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmin.y,pmin.z);
      const tdouble3 pm1=TDouble3(pm.x,pmin.y,pm.z);
      const tdouble3 pm2=TDouble3(pm.x,pmin.y+size05,pm.z);
      sh->AddShapeQuad(p1,p2,p3,p4,int(Id));
      sh->AddShapeQuadWire(p1,p2,p3,p4,int(Id));
      sh->AddShapeLine(pm1,pm2,int(Id));
    }
    if(PosMin.z!=-DBL_MAX){
      const tdouble3 p1=pmin;
      const tdouble3 p2=TDouble3(pmin.x,pmax.y,pmin.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmin.z);
      const tdouble3 p4=TDouble3(pmax.x,pmin.y,pmin.z);
      const tdouble3 pm1=TDouble3(pm.x,pm.y,pmin.z);
      const tdouble3 pm2=TDouble3(pm.x,pm.y,pmin.z+size05);
      sh->AddShapeQuad(p1,p2,p3,p4,int(Id));
      sh->AddShapeQuadWire(p1,p2,p3,p4,int(Id));
      sh->AddShapeLine(pm1,pm2,int(Id));
    }
    if(PosMax.x!=DBL_MAX){
      const tdouble3 p1=TDouble3(pmax.x,pmin.y,pmin.z);
      const tdouble3 p2=TDouble3(pmax.x,pmin.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmax.y,pmin.z);
      const tdouble3 pm1=TDouble3(pmax.x,pm.y,pm.z);
      const tdouble3 pm2=TDouble3(pmax.x-size05,pm.y,pm.z);
      sh->AddShapeQuad(p1,p2,p3,p4,int(Id));
      sh->AddShapeQuadWire(p1,p2,p3,p4,int(Id));
      sh->AddShapeLine(pm1,pm2,int(Id));
    }
    if(PosMax.y!=DBL_MAX){
      const tdouble3 p1=TDouble3(pmin.x,pmax.y,pmin.z);
      const tdouble3 p2=TDouble3(pmin.x,pmax.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmax.y,pmin.z);
      const tdouble3 pm1=TDouble3(pm.x,pmax.y,pm.z);
      const tdouble3 pm2=TDouble3(pm.x,pmax.y-size05,pm.z);
      sh->AddShapeQuad(p1,p2,p3,p4,int(Id));
      sh->AddShapeQuadWire(p1,p2,p3,p4,int(Id));
      sh->AddShapeLine(pm1,pm2,int(Id));
    }
    if(PosMax.z!=DBL_MAX){
      const tdouble3 p1=TDouble3(pmin.x,pmin.y,pmax.z);
      const tdouble3 p2=TDouble3(pmin.x,pmax.y,pmax.z);
      const tdouble3 p3=TDouble3(pmax.x,pmax.y,pmax.z);
      const tdouble3 p4=TDouble3(pmax.x,pmin.y,pmax.z);
      const tdouble3 pm1=TDouble3(pm.x,pm.y,pmax.z);
      const tdouble3 pm2=TDouble3(pm.x,pm.y,pmax.z-size05);
      sh->AddShapeQuad(p1,p2,p3,p4,int(Id));
      sh->AddShapeQuadWire(p1,p2,p3,p4,int(Id));
      sh->AddShapeLine(pm1,pm2,int(Id));
    }
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Pos::ComputeFilterCpu(unsigned np,const tdouble3* pos
  ,const typecode* code,byte* sel)const
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
  Point=Vector=TDouble3(0); 
  Distance=0;
  Plane=TPlane3d(0);
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Plane::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  Point=sxml->ReadElementDouble3(ele,"point");
  Vector=sxml->ReadElementDouble3(ele,"vector");
  Distance=sxml->ReadElementDouble(ele,"distance","v",true,DBL_MAX);
  Plane=fgeo::PlanePtVec(Point,Vector);
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Plane::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode:%s, lv:%d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Point:(%g,%g,%g)  Vector:(%g,%g,%g)  Dist:%g"
    ,Point.x,Point.y,Point.z,Vector.x,Vector.y,Vector.z,Distance));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Plane::SaveVtkConfig(double size,JVtkLib* sh)const{
  tdouble3 vec=fgeo::VecUnitary(Vector);
  sh->AddShapeQuad(Point,vec,size,int(Id));
  sh->AddShapeQuadWire(Point,vec,size,int(Id));
  const tdouble3 p2=Point+(vec*(Distance<DBL_MAX? Distance: size/2));
  sh->AddShapeLine(Point,p2,int(Id));
  if(Distance<DBL_MAX){
    sh->AddShapeQuad(p2,vec*-1.,size,int(Id));
    sh->AddShapeQuadWire(p2,vec*-1.,size,int(Id));
  }
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Plane::ComputeFilterCpu(unsigned np,const tdouble3* pos
  ,const typecode* code,byte* sel)const
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
  Centre=TDouble3(0); 
  Radius=0;
}
//==============================================================================
/// Read extra configuration from XML.
//==============================================================================
void JDsOutputPartsOp_Sphere::ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){
  Reset();
  Centre=sxml->ReadElementDouble3(ele,"centre");
  Radius=sxml->ReadElementFloat(ele,"radius","v");
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Sphere::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode:%s, lv:%d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Centre:(%g,%g,%g)  Radius:%g",Centre.x,Centre.y,Centre.z,Radius));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Sphere::SaveVtkConfig(double size,JVtkLib* sh)const{
  sh->AddShapeSphere(Centre,Radius,28,int(Id));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Sphere::ComputeFilterCpu(unsigned np,const tdouble3* pos
  ,const typecode* code,byte* sel)const
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
  Point1=sxml->ReadElementDouble3(ele,"point1");
  Point2=sxml->ReadElementDouble3(ele,"point2");
  Radius=sxml->ReadElementDouble(ele,"radius","v",true,DBL_MAX);
  //Plane=fgeo::PlaneNormalized(fgeo::PlanePtVec(Point,Vector));
  const tdouble3 vector=Point2-Point1;
  Plane=fgeo::PlanePtVec(Point1,vector);
  Distance=fgeo::PointDist(vector);
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputPartsOp_Cylinder::GetConfig(std::vector<std::string>& lines)const{
  const string mode=(CombineAnd? (Inverse? "Del": "Confirm"): "Add");
  lines.push_back(fun::PrintStr("  Filter_%s (mode:%s, lv:%d): ",GetNameType(Type).c_str(),mode.c_str(),BitResult));
  lines.push_back(fun::PrintStr("    Point1:(%g,%g,%g)  Point2:(%g,%g,%g)  Radius:%g"
    ,Point1.x,Point1.y,Point1.z,Point2.x,Point2.y,Point2.z,Radius));
}
//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputPartsOp_Cylinder::SaveVtkConfig(double size,JVtkLib* sh)const{
  sh->AddShapeCylinder(Point1,Point2,Radius,28,int(Id));
}
//==============================================================================
/// Compute filter for particles.
//==============================================================================
void JDsOutputPartsOp_Cylinder::ComputeFilterCpu(unsigned np,const tdouble3* pos
  ,const typecode* code,byte* sel)const
{
  const bool isvertical=(Point1.x==Point2.x && Point1.y==Point2.y);
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

      const double dist=fgeo::PlanePoint(Plane,ps);
      bool ok=(dist>=0 && dist<=Distance);
      if(ok && isvertical){
        const double dx=(ps.x-Point1.x);
        const double dy=(ps.y-Point1.y);
        ok=(dx*dx+dy*dy<=Radius);
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
void JDsOutputParts::ReadXml(const JXml* sxml,const TiXmlElement* ele){
  const string presel=fun::StrUpper(sxml->GetAttributeStr(ele,"preselection",false));
  if(presel!="ALL" && presel!="NONE")sxml->ErrReadAtrib(ele,"preselection",false,"");
  byte result=0;
  JDsOutputPartsOp_Init* opeini=new JDsOutputPartsOp_Init(Count(),presel=="ALL",result,sxml,ele);
  List.push_back(opeini);
  const TiXmlElement* ele2=ele->FirstChildElement(); 
  while(ele2){
    const string cmd=ele2->Value();
    if(cmd.length() && cmd[0]!='_' && sxml->CheckElementActive(ele2)){
      const string oper=fun::StrUpper(sxml->GetAttributeStr(ele2,"operation",false));
      if(oper!="ADD" && oper!="DEL" && oper!="CONFIRM")sxml->ErrReadAtrib(ele2,"operation",false,"");
      const bool cmband=(oper=="CONFIRM" || oper=="DEL");
      bool inverse=sxml->GetAttributeBool(ele2,"inverse",true,false);
      if(oper=="DEL")inverse=!inverse;
      if(cmd=="filterpos")List.push_back(new JDsOutputPartsOp_Pos(Count(),inverse,cmband,result,sxml,ele2));
      else if(cmd=="filterplane")List.push_back(new JDsOutputPartsOp_Plane(Count(),inverse,cmband,result,sxml,ele2));
      else if(cmd=="filtersphere")List.push_back(new JDsOutputPartsOp_Sphere(Count(),inverse,cmband,result,sxml,ele2));
      else if(cmd=="filtercylinder")List.push_back(new JDsOutputPartsOp_Cylinder(Count(),inverse,cmband,result,sxml,ele2));
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele2=ele2->NextSiblingElement();
  }
  SaveVtkConfig(VtkSize);
}
//==============================================================================
/// Loads configuration from XML object.
//==============================================================================
void JDsOutputParts::LoadXml(const JXml* sxml,const std::string& place){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place,false);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}
//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsOutputParts::GetConfig(std::string txhead,std::string txfoot
  ,std::vector<std::string>& lines)const
{
  if(!txhead.empty())lines.push_back(txhead); 
  for(unsigned c=0;c<Count();c++)List[c]->GetConfig(lines);
  if(!txfoot.empty())lines.push_back(txfoot); 
}

//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsOutputParts::SaveVtkConfig(double size)const{
  if(Count()){
    JVtkLib sh;
    for(unsigned c=0;c<Count();c++)List[c]->SaveVtkConfig(size,&sh);
    string filevtk=AppInfo.GetDirOut()+"CfgOutputParts_Scheme.vtk";
    sh.SaveShapeVtk(filevtk,(Count()>1? "num": ""));
    AppInfo.LogPtr()->AddFileInfo(filevtk,"Saves VTK file with OutputParts filters configurations.");
    sh.SaveShapeVtk(AppInfo.GetDirOut()+"CfgOutputParts_Scheme_dg01.vtk",(Count()>1? "num": ""));
    sh.SaveShapeVtk(AppInfo.GetDirOut()+"CfgOutputParts_Scheme_dg02.vtk",(Count()>1? "num": ""));
    sh.SaveShapeVtk(AppInfo.GetDirOut()+"CfgOutputParts_Scheme_dg03.vtk",(Count()>1? "num": ""));
  }
}

//==============================================================================
/// Compute filters for particles.
//==============================================================================
void JDsOutputParts::ComputeFilterCpu(unsigned np,const tdouble3* pos
  ,const typecode* code,byte* sel)const
{
  memset(sel,0,sizeof(byte)*np);
  for(unsigned cf=0;cf<Count();cf++)List[cf]->ComputeFilterCpu(np,pos,code,sel);
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

