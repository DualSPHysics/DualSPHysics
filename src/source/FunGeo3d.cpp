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

/// \file FunGeo3d.cpp \brief Implements geometry functions for 3D.

#include "FunGeo3d.h"
#include "FunctionsMath.h"
#include "Functions.h"

#include <cstdio>
#include <climits>
#include <algorithm>

namespace fgeo{

//==============================================================================
/// Devuelve proyeccion ortogonal del punto en la linea (pr1,pr2).
/// Returns orthogonal projection of the point in the line (pr1,pr2).
//==============================================================================
tdouble3 LineOrthogonalPoint(const tdouble3 &pt,const tdouble3 &pr1,const tdouble3 &pr2){
  return(PlaneLineIntersec(PlanePtVec(pt,pr2-pr1),pr1,pr2));
}

//==============================================================================
/// Devuelve proyeccion ortogonal del punto en la linea (pr1,pr2).
/// Returns orthogonal projection of the point in the line (pr1,pr2).
//==============================================================================
tfloat3 LineOrthogonalPoint(const tfloat3 &pt,const tfloat3 &pr1,const tfloat3 &pr2){
  return(PlaneLineIntersec(PlanePtVec(pt,pr2-pr1),pr1,pr2));
}

  
//==============================================================================
/// Devuelve el plano formado por 3 puntos. La normal es (a,b,c)
/// Returns the plane defined by 3 points. The normal is (a,b,c)
//==============================================================================
tplane3d Plane3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3){ 
  tplane3d plano={0,0,0,0};
  if(p1!=p2 && p1!=p3 && p2!=p3){
    const tdouble3 v1=TDouble3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
    const tdouble3 v2=TDouble3(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
    const tdouble3 v=ProductVec(v1,v2); //-Vector normal del plano. //-Normal plane vector
    plano.a=v.x; plano.b=v.y; plano.c=v.z; 
    plano.d=-((v.x*p1.x)+(v.y*p1.y)+(v.z*p1.z));
  }
  return(plano);
}

//==============================================================================
/// Devuelve el plano formado por 3 puntos. La normal es (a,b,c)
/// Returns the plane defined by 3 points. The normal is (a,b,c)
//==============================================================================
tplane3f Plane3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3){ 
  tplane3f plano={0,0,0,0};
  if(p1!=p2 && p1!=p3 && p2!=p3){
    const tfloat3 v1=TFloat3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
    const tfloat3 v2=TFloat3(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
    const tfloat3 v=ProductVec(v1,v2); //-Vector normal del plano. //-Normal plane vector
    plano.a=v.x; plano.b=v.y; plano.c=v.z; 
    plano.d=-((v.x*p1.x)+(v.y*p1.y)+(v.z*p1.z));
  }
  return(plano);
}


//==============================================================================
/// Devuelve true cuando todos los puntos estan en el plano.
/// Returns true when all points are in the plane.
//==============================================================================
bool PlanePointsIn(const tplane3d &pla,unsigned np,const tdouble3 *vpt,double tolerance){
  bool ret=true;
  for(unsigned p=0;p<np && ret;p++)ret=(fgeo::PlaneDist(pla,vpt[p])<=tolerance);
  return(ret);
}

//==============================================================================
/// Devuelve true cuando todos los puntos estan en el plano.
/// Returns true when all points are in the plane.
//==============================================================================
bool PlanePointsIn(const tplane3f &pla,unsigned np,const tfloat3 *vpt,float tolerance){
  bool ret=true;
  for(unsigned p=0;p<np && ret;p++)ret=(fgeo::PlaneDist(pla,vpt[p])<=tolerance);
  return(ret);
}


//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si usando
/// la Regla de Cramer (para sistemas compatibles determinados).
/// Returns point of intersection of three non-parallel planes using
/// Cramer's rule (for determinants of compatible systems).
//==============================================================================
tdouble3 PlanesIntersec(const tplane3d &pla1,const tplane3d &pla2,const tplane3d &pla3){
  tdouble3 res=TDouble3(0);
  const double dm=fmath::Determinant3x3(TMatrix3d(pla1.a,pla1.b,pla1.c,pla2.a,pla2.b,pla2.c,pla3.a,pla3.b,pla3.c));
  if(dm){
    const double dx=fmath::Determinant3x3(TMatrix3d(-pla1.d,pla1.b,pla1.c,-pla2.d,pla2.b,pla2.c,-pla3.d,pla3.b,pla3.c));
    const double dy=fmath::Determinant3x3(TMatrix3d(pla1.a,-pla1.d,pla1.c,pla2.a,-pla2.d,pla2.c,pla3.a,-pla3.d,pla3.c));
    const double dz=fmath::Determinant3x3(TMatrix3d(pla1.a,pla1.b,-pla1.d,pla2.a,pla2.b,-pla2.d,pla3.a,pla3.b,-pla3.d));
    res=TDouble3(dx/dm,dy/dm,dz/dm);
  }
  return(res);
}

//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si usando
/// la Regla de Cramer (para sistemas compatibles determinados).
/// Returns point of intersection of three non-parallel planes using
/// Cramer's rule (for determinants of compatible systems).
//==============================================================================
tfloat3 PlanesIntersec(const tplane3f &pla1,const tplane3f &pla2,const tplane3f &pla3){
  tfloat3 res=TFloat3(0);
  const float dm=fmath::Determinant3x3(TMatrix3f(pla1.a,pla1.b,pla1.c,pla2.a,pla2.b,pla2.c,pla3.a,pla3.b,pla3.c));
  if(dm){
    const float dx=fmath::Determinant3x3(TMatrix3f(-pla1.d,pla1.b,pla1.c,-pla2.d,pla2.b,pla2.c,-pla3.d,pla3.b,pla3.c));
    const float dy=fmath::Determinant3x3(TMatrix3f(pla1.a,-pla1.d,pla1.c,pla2.a,-pla2.d,pla2.c,pla3.a,-pla3.d,pla3.c));
    const float dz=fmath::Determinant3x3(TMatrix3f(pla1.a,pla1.b,-pla1.d,pla2.a,pla2.b,-pla2.d,pla3.a,pla3.b,-pla3.d));
    res=TFloat3(dx/dm,dy/dm,dz/dm);
  }
  return(res);
}


//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tdouble3 PlaneLineIntersec(const tplane3d &pla,const tdouble3 &pt1,const tdouble3 &pt2){
  const tdouble3 v=(pt2-pt1);
  const double x1=pt1.x;
  const double x2=v.x;
  const double y1=pt1.y;
  const double y2=v.y;
  const double z1=pt1.z;
  const double z2=v.z;
  const double t=-(pla.a*x1+pla.b*y1+pla.c*z1+pla.d)/(pla.a*x2+pla.b*y2+pla.c*z2);
  return(TDouble3(x1+x2*t,y1+y2*t,z1+z2*t));
}

//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tfloat3 PlaneLineIntersec(const tplane3f &pla,const tfloat3 &pt1,const tfloat3 &pt2){
  const tfloat3 v=(pt2-pt1);
  const float x1=pt1.x;
  const float x2=v.x;
  const float y1=pt1.y;
  const float y2=v.y;
  const float z1=pt1.z;
  const float z2=v.z;
  const float t=-(pla.a*x1+pla.b*y1+pla.c*z1+pla.d)/(pla.a*x2+pla.b*y2+pla.c*z2);
  return(TFloat3(x1+x2*t,y1+y2*t,z1+z2*t));
}


//==============================================================================
/// Devuelve distancia maxima entre plano y 4 puntos.
/// Returns maximum distance between plane and 4 points.
//==============================================================================
double PlaneDistMax(const tplane3d &pla,const tdouble3 &p1,const tdouble3 &p2
  ,const tdouble3 &p3,const tdouble3 &p4)
{ 
  const double d1=PlaneDistSign(pla,p1);
  const double d2=PlaneDistSign(pla,p2);
  const double d3=PlaneDistSign(pla,p3);
  const double d4=PlaneDistSign(pla,p4);
  const double d=std::max(std::max(d1,d2),std::max(d3,d4));
  return(fabs(d));
}

//==============================================================================
/// Devuelve planos y distancias para delimitar un dominio en forma de caja.
/// Returns planes and distnaces to limit a box domain.
//==============================================================================
void PlanesDomain(const tdouble3 &pt,const tdouble3 &vx,const tdouble3 &vy
  ,const tdouble3 &vz,tplane3d &plax,tplane3d &play,tplane3d &plaz,tdouble3 &pladist)
{
  const tdouble3 p1=pt;
  const tdouble3 p2=p1+vx;
  const tdouble3 p3=p2+vz;
  const tdouble3 p4=p1+vz;
  const tdouble3 p5=p1+vy;
  const tdouble3 p6=p2+vy;
  const tdouble3 p7=p3+vy;
  const tdouble3 p8=p4+vy;
  const tplane3d plx=Plane3Pt(p1,p8,p4);
  const double dx=PlaneDistMax(plx,p2,p3,p6,p7)*sqrt(plx.a*plx.a+plx.b*plx.b+plx.c*plx.c);
  const tplane3d ply=Plane3Pt(p1,p3,p2);
  const double dy=PlaneDistMax(ply,p5,p6,p7,p8)*sqrt(ply.a*ply.a+ply.b*ply.b+ply.c*ply.c);
  const tplane3d plz=Plane3Pt(p1,p2,p6);
  const double dz=PlaneDistMax(plz,p3,p4,p7,p8)*sqrt(plz.a*plz.a+plz.b*plz.b+plz.c*plz.c);
  plax=plx; play=ply; plaz=plz; pladist=TDouble3(dx,dy,dz);
}

//==============================================================================
/// Comprueba si el punto esta dentro del dominio definido.
/// Checks the point is inside the defined domain.
//==============================================================================
bool PlanesDomainCheck(const tdouble3 &pt,const tplane3d &plax,const tplane3d &play
  ,const tplane3d &plaz,const tdouble3 &pladist)
{
  const double dx=PlanePoint(plax,pt);
  const double dy=PlanePoint(play,pt);
  const double dz=PlanePoint(plaz,pt);
  return(dx>=0 && dx<=pladist.x && dy>=0 && dy<=pladist.y && dz>=0 && dz<=pladist.z);
}


//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
/// Returns the area of a triangle formed by 3 points.
//==============================================================================
double TriangleArea(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3){
  //Se obtienen los vectores del triangulo.
  //Obtains the triangle vectors.
  double PQx=p2.x-p1.x;
  double PQy=p2.y-p1.y;
  double PQz=p2.z-p1.z;
  double PRx=p3.x-p1.x;
  double PRy=p3.y-p1.y;
  double PRz=p3.z-p1.z;
  //Se hace el producto cruz.
  //Computes the cross product.
  double Vi=PQy*PRz-PRy*PQz;
  double Vj=-(PQx*PRz-PRx*PQz);
  double Vk=PQx*PRy-PRx*PQy;
  //Se obtiene el area del triangulo que es igual a la mitad de la magnitud del vector resultante.
  //Obtains the triangle area that equals half the magnitude of the resulting vector.
  return(double(.5)*sqrt(Vi*Vi+Vj*Vj+Vk*Vk));
}

//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
/// Returns the area of a triangle formed by 3 points.
//==============================================================================
float TriangleArea(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3){
  //Se obtienen los vectores del triangulo.
  //Obtains the triangle vectors.
  float PQx=p2.x-p1.x;
  float PQy=p2.y-p1.y;
  float PQz=p2.z-p1.z;
  float PRx=p3.x-p1.x;
  float PRy=p3.y-p1.y;
  float PRz=p3.z-p1.z;
  //Se hace el producto cruz.
  //Computes the cross product.
  float Vi=PQy*PRz-PRy*PQz;
  float Vj=-(PQx*PRz-PRx*PQz);
  float Vk=PQx*PRy-PRx*PQy;
  //Se obtiene el area del triangulo que es igual a la mitad de la magnitud del vector resultante.
  //Obtains the triangle area that equals half the magnitude of the resulting vector.
  return(float(.5)*sqrt(Vi*Vi+Vj*Vj+Vk*Vk));
}


}


