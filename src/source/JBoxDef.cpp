//HEAD_DSCODES
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

/// \file JBoxDef.cpp \brief Implements the class \ref JBoxDef.

#include "JBoxDef.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "FunGeo3d.h"
#include "JMatrix4.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cfloat>

using namespace std;

//##############################################################################
//# JBoxDef
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JBoxDef::JBoxDef(bool is2d,double posy2d){
  ClassName="JBoxDef";
  Reset();
  Is2D=is2d;
  Posy2D=posy2d;
}

//==============================================================================
/// Constructor for simple mode.
//==============================================================================
JBoxDef::JBoxDef(const tdouble3& pmin,const tdouble3& pmax,bool is2d)
  :JBoxDef(is2d,(is2d ? pmin.y : 0))
{
  SetPos(pmin,pmax);
}

//==============================================================================
/// Constructor for non-simple mode.
//==============================================================================
JBoxDef::JBoxDef(const tdouble3& pmin,const tdouble3& vx,const tdouble3& vy
  ,const tdouble3& vz,bool is2d):JBoxDef(is2d,(is2d ? pmin.y : 0))
{
  SetPos(pmin,vx,vy,vz);
}

//==============================================================================
/// Constructor for non-simple mode (simple + transformation).
//==============================================================================
JBoxDef::JBoxDef(const tdouble3& pmin,const tdouble3& pmax,bool is2d
  ,const tdouble3& mov,const tdouble3& rot,const tdouble3& rotcen)
  :JBoxDef(is2d,(is2d ? pmin.y : 0))
{
  SetPos(pmin,pmax);
  MovRotate(mov,rot,rotcen);
}

//==============================================================================
/// Constructor for copy.
//==============================================================================
JBoxDef::JBoxDef(const JBoxDef& src):JBoxDef(src.Is2D,src.Posy2D){
  *this=src;
}

//==============================================================================
/// Overload assignment operator.
//==============================================================================
JBoxDef& JBoxDef::operator=(const JBoxDef& src){
  if(this!=&src){
    Is2D=src.Is2D;
    Posy2D=src.Posy2D;
    Simple=src.Simple;
    PosMin=src.PosMin;
    PosMax=src.PosMax;
    Vx=src.Vx;
    Vy=src.Vy;
    Vz=src.Vz;
  }
  return(*this);
}

//==============================================================================
/// Destructor.
//==============================================================================
JBoxDef::~JBoxDef(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JBoxDef::Reset(){
  Is2D=false;
  Posy2D=0;
  Simple=true;
  PosMin=PosMax=TDouble3(0,Posy2D,0);
  Vx=Vy=Vz=TDouble3(0);
}

//==============================================================================
/// Set minimum and maximum postion as simple mode.
//==============================================================================
void JBoxDef::SetPos(const tdouble3& pmin,const tdouble3& pmax){
  Simple=true;
  PosMin=pmin;
  PosMax=pmax;
  Vx=TDouble3(PosMax.x-PosMin.x,0,0);
  Vy=TDouble3(0,PosMax.y-PosMin.y,0);
  Vz=TDouble3(0,0,PosMax.z-PosMin.z);
  if(!(PosMin<=PosMax))Run_Exceptioon("Minimum limit is higher than the maximum limit.");
  if(Is2D && (PosMin.y!=Posy2D || PosMax.y!=Posy2D))
    Run_Exceptioon("Y value is invalid for current 2D configuration.");
}

//==============================================================================
/// Set minimum postion and vectors as non-simple mode.
//==============================================================================
void JBoxDef::SetPos(const tdouble3& pmin,const tdouble3& vx,const tdouble3& vy
  ,const tdouble3& vz)
{
  Simple=false;
  PosMin=pmin;
  Vx=vx;
  Vy=vy;
  Vz=vz;
  PosMax=PosMin + (vx + vy + vz);
  if(Is2D && (PosMin.y!=Posy2D || PosMax.y!=Posy2D))
    Run_Exceptioon("Y value is invalid for current 2D configuration.");
}

//==============================================================================
/// Applies translation and rotation and converts to non-simple mode.
//==============================================================================
void JBoxDef::MovRotate(const tdouble3& mov,const tdouble3& rot
  ,const tdouble3& cen)
{
  if(mov!=TDouble3(0) || rot!=TDouble3(0)){
    //-Converts to non-simple mode.
    Simple=false;
    //-Computes transformation matrix.
    JMatrix4d mat;
    if(mov!=TDouble3(0))mat.Move(mov);
    if(rot!=TDouble3(0)){
      if(cen!=TDouble3(0))mat.Move(cen);
      mat.Rotate(rot);
      if(cen!=TDouble3(0))mat.Move(cen*-1.);
    }
    //-Checks 2D limitations.
    if(Is2D && mov.y)
      Run_Exceptioon("Y-translation is invalid for current 2D configuration.");
    if(Is2D && (rot.x || rot.z))
      Run_Exceptioon("Only Y-rotation is valid for current 2D configuration.");
    //-Applies transformation matrix.
    PosMin=mat.MulPoint(PosMin);
    PosMax=mat.MulPoint(PosMax);
    Vx=mat.MulNormal(Vx);
    Vy=mat.MulNormal(Vy);
    Vz=mat.MulNormal(Vz);
  }
}

//==============================================================================
/// Applies matrix transformation and converts to non-simple mode.
//==============================================================================
void JBoxDef::MovRotate(const tmatrix4d& mat){
  if(mat!=TMatrix4d()){
    //-Converts to non-simple mode.
    Simple=false;
    //-Applies transformation matrix.
    JMatrix4d jmat(mat);
    PosMin=jmat.MulPoint(PosMin);
    PosMax=jmat.MulPoint(PosMax);
    Vx=jmat.MulNormal(Vx);
    Vy=jmat.MulNormal(Vy);
    Vz=jmat.MulNormal(Vz);
    //-Checks invalid tranformations.
    if(Is2D && (PosMin.y!=PosMax.y || PosMin.y!=Posy2D))
      Run_Exceptioon("Y translation or rotation is invalid for current 2D configuration.");
  }
}

//==============================================================================
/// Extend or reduce limits around the box domain.
//==============================================================================
void JBoxDef::Resize(double resize){
  if(resize){
    if(Simple){
      const tdouble3 rs=TDouble3(resize,(Is2D? 0: resize),resize);
      PosMin=PosMin-rs;
      PosMax=PosMax+rs;
      Vx=TDouble3(PosMax.x-PosMin.x,0,0);
      Vy=TDouble3(0,PosMax.y-PosMin.y,0);
      Vz=TDouble3(0,0,PosMax.z-PosMin.z);
    }
    else{
      const tdouble3 rs=fgeo::VecModule(Vx,resize) 
        + (Is2D? TDouble3(0): fgeo::VecModule(Vy,resize))
        + fgeo::VecModule(Vz,resize);
      //-Set minimum position.
      PosMin=PosMin-rs;
      //-Resize vectors.
      const double d1=fgeo::PointDist(Vx);
      Vx=Vx*((d1+resize+resize)/d1);
      if(!Is2D){
        const double d2=fgeo::PointDist(Vy);
        Vy=Vy*((d2+resize+resize)/d2);
      }
      const double d3=fgeo::PointDist(Vz);
      Vz=Vz*((d3+resize+resize)/d3);
      //-Set maximum position.
      PosMax=PosMin+(Vx+Vy+Vz);
    }
  }
}

//==============================================================================
/// Returns minimum valid size.
//==============================================================================
double JBoxDef::GetMinSize()const{
  double smin=0;
  if(Simple){
    const tdouble3 size=PosMax-PosMin;
    smin=min(size.x,size.z);
    if(!Is2D)smin=min(smin,size.y);
  }
  else{
    const double dx=fgeo::PointDist(Vx);
    const double dz=fgeo::PointDist(Vz);
    const double dy=(Is2D? 0: fgeo::PointDist(Vy));
    smin=min(dx,dz);
    if(!Is2D)smin=min(smin,dy);
  }
  return(smin);
}

//==============================================================================
/// Returns box definition as string.
//==============================================================================
std::string JBoxDef::GetDomainStr()const{
  string tx;
  if(Simple)tx=fun::Double3gRangeStr(PosMin,PosMax);
  else tx=fun::PrintStr("(%g,%g,%g)-[(%g,%g,%g)+(%g,%g,%g)+(%g,%g,%g)]"
    ,PosMin.x,PosMin.y,PosMin.z ,Vx.x,Vx.y,Vx.z ,Vy.x,Vy.y,Vy.z ,Vz.x,Vz.y,Vz.z);
  return(tx);
}

//==============================================================================
/// Returns 8 points of definition.
//==============================================================================
void JBoxDef::GetBoxPoints(std::vector<tdouble3>& points)const{
  points.resize(8);
  GetBoxPoints(points.data());
}

//==============================================================================
/// Returns 8 points of definition.
//==============================================================================
void JBoxDef::GetBoxPoints(tdouble3* points)const{
  if(Simple){
    points[0]=PosMin;
    points[1]=TDouble3(PosMax.x,PosMin.y,PosMin.z);
    points[2]=TDouble3(PosMax.x,PosMax.y,PosMin.z);
    points[3]=TDouble3(PosMin.x,PosMax.y,PosMin.z);
    points[4]=TDouble3(PosMin.x,PosMin.y,PosMax.z);
    points[5]=TDouble3(PosMax.x,PosMin.y,PosMax.z);
    points[6]=TDouble3(PosMax.x,PosMax.y,PosMax.z);
    points[7]=TDouble3(PosMin.x,PosMax.y,PosMax.z);
  }
  else{
    points[0]=PosMin;
    points[1]=PosMin + Vx;
    points[2]=PosMin + Vx + Vy;
    points[3]=PosMin + Vy;
    points[4]=points[0] + Vz;
    points[5]=points[1] + Vz;
    points[6]=points[2] + Vz;
    points[7]=points[3] + Vz;
  }
}

//==============================================================================
/// Returns minimum position of the simple box that includes it.
//==============================================================================
tdouble3 JBoxDef::OutBoxMin()const{
  tdouble3 pt=PosMin;
  if(!Simple){
    tdouble3 points[8];
    GetBoxPoints(points);
    for(unsigned c=1;c<8;c++)pt=MinValues(pt,points[c]);
  }
  return(pt);
}

//==============================================================================
/// Returns maximum position of the simple box that includes it.
//==============================================================================
tdouble3 JBoxDef::OutBoxMax()const{
  tdouble3 pt=PosMax;
  if(!Simple){
    tdouble3 points[8];
    GetBoxPoints(points);
    for(unsigned c=1;c<8;c++)pt=MaxValues(pt,points[c]);
  }
  return(pt);
}

//==============================================================================
/// Returns true when box is fully inside.
//==============================================================================
bool JBoxDef::Inside(const JBoxDef& box)const{
  bool ret=false;
  if(Simple){
    if(box.Simple)ret=(box.PosMin>=PosMin && box.PosMax<=PosMax);
    else{
      const tdouble3 bmin=box.OutBoxMin();
      const tdouble3 bmax=box.OutBoxMax();
      ret=(bmin>=PosMin && bmax<=PosMax);
    }
  }
  else{
    //-Get out-box points.
    tdouble3 points[8];
    box.GetBoxPoints(points);
    //-Get size values and planes.
    const double vdis[3]={fgeo::PointDist(Vx),
                          fgeo::PointDist(Vz),
                          fgeo::PointDist(Vy)};
    const tplane3d vpla[3]={fgeo::PlaneAxisDist(PosMin,Vx),
                            fgeo::PlaneAxisDist(PosMin,Vz),
        (Is2D? TPlane3d(0): fgeo::PlaneAxisDist(PosMin,Vy))};
    const unsigned na=(Is2D? 2: 3);
    //-Check if pointe of box are inside.
    ret=true;
    for(unsigned a=0;a<na && ret;a++)for(unsigned c=0;c<8 && ret;c++){
      const double d=fgeo::PlanePoint(vpla[a],points[c]); //-Computes distance.
      if(d<0 || d>vdis[a])ret=false;
    }
  }
  return(ret);
}

//==============================================================================
/// Loads planes to check positions inside box in vplanes[] and return the 
/// number of planes (3 for 3-D boxes and 2 for 2-D boxes). 
//==============================================================================
unsigned JBoxDef::GetInsidePlanes(tplane3d* vplanes)const{
  //-Get size values.
  const double dx=fgeo::PointDist(Vx);
  const double dz=fgeo::PointDist(Vz);
  const double dy=fgeo::PointDist(Vy);
  vplanes[0]=fgeo::PlaneAxisDist(PosMin,Vx,dx);
  vplanes[1]=fgeo::PlaneAxisDist(PosMin,Vz,dz);
  vplanes[2]=(Is2D? TPlane3d(0): fgeo::PlaneAxisDist(PosMin,Vy,dy));
  return(Is2D? 2: 3);
}

