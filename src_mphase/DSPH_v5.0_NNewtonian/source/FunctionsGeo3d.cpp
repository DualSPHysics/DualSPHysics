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

/// \file FunctionsGeo3d.cpp \brief Implements geometry functions for 3D.

#include "FunctionsGeo3d.h"
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
/// Devuelve proyeccion ortogonal del punto en la linea (pr1,pr2).
/// Returns orthogonal projection of the point in the line (pr1,pr2).
//==============================================================================
double LineOrthogonalPointFromPr1(const tdouble3 &pt,const tdouble3 &pr1,const tdouble3 &pr2){
  const tdouble3 pr=LineOrthogonalPoint(pt,pr1,pr2);
  const double vx=fabs(pr2.x-pr1.x);
  const double vy=fabs(pr2.y-pr1.y);
  const double vz=fabs(pr2.z-pr1.z);
  double t=0;
  if(vx>0 || vy>0 || vz>0){
    if(vx>=vy && vx>=vz)t=(pr.x-pr1.x)/(pr2.x-pr1.x);
    else if(vy>=vx && vy>=vz)t=(pr.y-pr1.y)/(pr2.y-pr1.y);
    else t=(pr.z-pr1.z)/(pr2.z-pr1.z);
  }
  return(t);
}

//==============================================================================
/// Devuelve proyeccion ortogonal del punto en la linea (pr1,pr2).
/// Returns orthogonal projection of the point in the line (pr1,pr2).
//==============================================================================
float LineOrthogonalPointFromPr1(const tfloat3 &pt,const tfloat3 &pr1,const tfloat3 &pr2){
  const tfloat3 pr=LineOrthogonalPoint(pt,pr1,pr2);
  const float vx=fabs(pr2.x-pr1.x);
  const float vy=fabs(pr2.y-pr1.y);
  const float vz=fabs(pr2.z-pr1.z);
  float t=0;
  if(vx>0 || vy>0 || vz>0){
    if(vx>=vy && vx>=vz)t=(pr.x-pr1.x)/(pr2.x-pr1.x);
    else if(vy>=vx && vy>=vz)t=(pr.y-pr1.y)/(pr2.y-pr1.y);
    else t=(pr.z-pr1.z)/(pr2.z-pr1.z);
  }
  return(t);
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
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal planes.
//==============================================================================
void TriangleNormalPlanes(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3
  ,double openingdist,tplane3d &pla1,tplane3d &pla2,tplane3d &pla3)
{
  //-Calcula normal unitaria.
  //-Computes unit normal.
  const tdouble3 nor=VecUnitary(TriangleNormal(p1,p2,p3));
  //-Calcula planos normales a triangulo (p1,p2,p3).
  //-Computes normal triangle planes (p1,p2,p3).
  const tdouble3 p1n=p1+nor;
  const tdouble3 p2n=p2+nor;
  tplane3d plad1=Plane3Pt(p1,p2,p1n);
  tplane3d plad2=Plane3Pt(p2,p3,p2n);
  tplane3d plad3=Plane3Pt(p3,p1,p1n);
  //-Abre planos normales.
  //-Opens normal planes.
  if(openingdist){
    double md=openingdist/sqrt(plad1.a*plad1.a+plad1.b*plad1.b+plad1.c*plad1.c);
    tdouble3 v=TDouble3(plad1.a*md,plad1.b*md,plad1.c*md);
    plad1=Plane3Pt(p1+v,p2+v,p1n+v);
    md=openingdist/sqrt(plad2.a*plad2.a+plad2.b*plad2.b+plad2.c*plad2.c);
    v=TDouble3(plad2.a*md,plad2.b*md,plad2.c*md);
    plad2=Plane3Pt(p2+v,p3+v,p2n+v);
    md=openingdist/sqrt(plad3.a*plad3.a+plad3.b*plad3.b+plad3.c*plad3.c);
    v=TDouble3(plad3.a*md,plad3.b*md,plad3.c*md);
    plad3=Plane3Pt(p3+v,p1+v,p1n+v);
  }
  //-Guarda planos normales.
  //-Stores normal planes.
  pla1=plad1; pla2=plad2; pla3=plad3;
}

//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal planes.
//==============================================================================
void TriangleNormalPlanes(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3
  ,float openingdist,tplane3f &pla1,tplane3f &pla2,tplane3f &pla3)
{
  //-Calcula normal unitaria.
  //-Computes unit normal.
  const tfloat3 nor=VecUnitary(TriangleNormal(p1,p2,p3));
  //-Calcula planos normales a triangulo (p1,p2,p3).
  //-Computes normal triangle planes (p1,p2,p3).
  const tfloat3 p1n=p1+nor;
  const tfloat3 p2n=p2+nor;
  tplane3f plad1=Plane3Pt(p1,p2,p1n);
  tplane3f plad2=Plane3Pt(p2,p3,p2n);
  tplane3f plad3=Plane3Pt(p3,p1,p1n);
  //-Abre planos normales.
  //-Opens normal planes.
  if(openingdist){
    float md=openingdist/sqrt(plad1.a*plad1.a+plad1.b*plad1.b+plad1.c*plad1.c);
    tfloat3 v=TFloat3(plad1.a*md,plad1.b*md,plad1.c*md);
    plad1=Plane3Pt(p1+v,p2+v,p1n+v);
    md=openingdist/sqrt(plad2.a*plad2.a+plad2.b*plad2.b+plad2.c*plad2.c);
    v=TFloat3(plad2.a*md,plad2.b*md,plad2.c*md);
    plad2=Plane3Pt(p2+v,p3+v,p2n+v);
    md=openingdist/sqrt(plad3.a*plad3.a+plad3.b*plad3.b+plad3.c*plad3.c);
    v=TFloat3(plad3.a*md,plad3.b*md,plad3.c*md);
    plad3=Plane3Pt(p3+v,p1+v,p1n+v);
  }
  //-Guarda planos normales.
  //-Stores normal planes.
  pla1=plad1; pla2=plad2; pla3=plad3;
}


//==============================================================================
/// Devuelve true cuando el punto esta dentro de 3 planos.
/// Returns true when the point is inside 3 planes.
//==============================================================================
bool TriangleInside(const tdouble3 &pt,const tplane3d &pla1,const tplane3d &pla2,const tplane3d &pla3){
  const double p1=PlanePoint(pla1,pt);
  const double p2=PlanePoint(pla2,pt);
  const double p3=PlanePoint(pla3,pt);
  //printf("--->pt(%f,%f,%f)=  %f  %f  %f\n",pt.x,pt.y,pt.z,p1,p2,p3);
  return((p1<=0 && p2<=0 && p3<=0) || (p1>=0 && p2>=0 && p3>=0));
}

//==============================================================================
/// Devuelve true cuando el punto esta dentro de 3 planos.
/// Returns true when the point is inside 3 planes.
//==============================================================================
bool TriangleInside(const tfloat3 &pt,const tplane3f &pla1,const tplane3f &pla2,const tplane3f &pla3){
  const double p1=PlanePoint(pla1,pt);
  const double p2=PlanePoint(pla2,pt);
  const double p3=PlanePoint(pla3,pt);
  return((p1<=0 && p2<=0 && p3<=0) || (p1>=0 && p2>=0 && p3>=0));
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


//==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
/// Starting from a triangle formed by 3 points returns the points that form
/// a triangle more or less open according to openingdist.
//==============================================================================
void TriangleOpen(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3
  ,double openingdist,tdouble3 &pt1,tdouble3 &pt2,tdouble3 &pt3)
{
  //-Calcula plano de triangulo.
  //-Computes plane of the triangle.
  const tplane3d pla=Plane3Pt(p1,p2,p3);
  //-Calcula planos normales abiertos.
  //-Computes open normal planes.
  tplane3d pla1n;
  tplane3d pla2n;
  tplane3d pla3n;
  TriangleNormalPlanes(p1,p2,p3,openingdist,pla1n,pla2n,pla3n);
  //-Calcula interseccion entre planos.
  //-Computes the intersection of three planes
  pt1=PlanesIntersec(pla,pla1n,pla3n);
  pt2=PlanesIntersec(pla,pla1n,pla2n);
  pt3=PlanesIntersec(pla,pla2n,pla3n);
}

//==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
/// Starting from a triangle formed by 3 points returns the points that form
/// a triangle more or less open according to openingdist.
//==============================================================================
void TriangleOpen(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3
  ,float openingdist,tfloat3 &pt1,tfloat3 &pt2,tfloat3 &pt3)
{
  //-Calcula plano de triangulo.
  //-Computes plane of the triangle.
  const tplane3f pla=Plane3Pt(p1,p2,p3);
  //-Calcula planos normales abiertos.
  //-Computes open normal planes.
  tplane3f pla1n;
  tplane3f pla2n;
  tplane3f pla3n;
  TriangleNormalPlanes(p1,p2,p3,openingdist,pla1n,pla2n,pla3n);
  //-Calcula interseccion entre planos.
  //-Computes the intersection of three planes
  pt1=PlanesIntersec(pla,pla1n,pla3n);
  pt2=PlanesIntersec(pla,pla1n,pla2n);
  pt3=PlanesIntersec(pla,pla2n,pla3n);
}


//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal planes.
//==============================================================================
void PolygonNormalPlanes(const std::vector<tdouble3> &vpt
  ,double openingdist,std::vector<tplane3d> &vpla)
{
  vpla.clear();
  const unsigned np=unsigned(vpt.size());
  if(np>2){
    //-Calcula normal unitaria.
    //-Computes unit normal.
    const tdouble3 nor=VecUnitary(TriangleNormal(vpt[0],vpt[1],vpt[2]));
    //-Calcula planos normales.
    //-Computes normal planes.
    for(unsigned p=0;p<np;p++){
      const unsigned p2=(p+1<np? p+1: 0);
      vpla.push_back(Plane3Pt(vpt[p],vpt[p2],vpt[p]+nor));
    }
    //-Abre planos normales.
    //-Opens normal planes.
    if(openingdist)for(unsigned p=0;p<np;p++){
      const tplane3d &pla=vpla[p];
      double md=openingdist/sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c);
      tdouble3 v=TDouble3(pla.a*md,pla.b*md,pla.c*md);
      const unsigned p2=(p+1<np? p+1: 0);
      vpla[p]=Plane3Pt(vpt[p]+v,vpt[p2]+v,vpt[p]+nor+v);
    }
  }
}

//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal planes.
//==============================================================================
void PolygonNormalPlanes(const std::vector<tfloat3> &vpt
  ,float openingdist,std::vector<tplane3f> &vpla)
{
  vpla.clear();
  const unsigned np=unsigned(vpt.size());
  if(np>2){
    //-Calcula normal unitaria.
    //-Computes unit normal.
    const tfloat3 nor=VecUnitary(TriangleNormal(vpt[0],vpt[1],vpt[2]));
    //-Calcula planos normales.
    //-Computes normal planes.
    for(unsigned p=0;p<np;p++){
      const unsigned p2=(p+1<np? p+1: 0);
      vpla.push_back(Plane3Pt(vpt[p],vpt[p2],vpt[p]+nor));
    }
    //-Abre planos normales.
    //-Opens normal planes.
    if(openingdist)for(unsigned p=0;p<np;p++){
      const tplane3f &pla=vpla[p];
      float md=openingdist/sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c);
      tfloat3 v=TFloat3(pla.a*md,pla.b*md,pla.c*md);
      const unsigned p2=(p+1<np? p+1: 0);
      vpla[p]=Plane3Pt(vpt[p]+v,vpt[p2]+v,vpt[p]+nor+v);
    }
  }
}


//==============================================================================
/// Indica si un poligon es concavo.
/// Returns true when the polygon is concave.
//==============================================================================
bool PolygonIsConcave(unsigned np,const tdouble3* vp,double tolerance){
  bool concave=false;
  tdouble3 vec0;
  for(unsigned p=0;p<np && !concave;p++){
    const tdouble3 p0=vp[(p+np-1)%np];
    const tdouble3 p1=vp[p%np];
    const tdouble3 p2=vp[(p+1)%np];
    tdouble3 vec=fgeo::VecUnitary(fgeo::ProductVec(p0-p1,p2-p1));
    //printf("pvec[%d](%g,%g,%g):: (%g,%g,%g)\n",p,p1.x,p1.y,p1.z,vec.x,vec.y,vec.z);
    if(!p)vec0=vec;
    else{
      if(!fun::IsEqual(vec0,vec,tolerance))concave=true;
      //if(concave)printf("concave\n");
    }
  }
  return(concave);
}

//==============================================================================
/// Devuelve poligono (convexo o concavo) como un conjunto de triangulos.
/// Returns polygon (convex or concave) as a set of triangles.
//==============================================================================
void PolygonConcave(unsigned np,const tdouble3* vpt,std::vector<tint3> &vtri,double tolerance){
  //printf("\nPolygonConcave\n");
  vtri.clear();
  unsigned ntri0=0;
  if(np==3)vtri.push_back(TInt3(0,1,2));
  else if(np>3){
    //-Allocate memory.
    int* vref=new int[np];
    const tdouble3 vec0=fgeo::VecUnitary(fgeo::ProductVec(vpt[np-1]-vpt[0],vpt[1]-vpt[0]));
    for(int mode=0;mode<2;mode++){
      //printf("\nMode_%d\n",mode);
      if(mode==1)ntri0=unsigned(vtri.size());
      //-Initialise polygon data.
      for(unsigned p=0;p<np;p++)vref[p]=int(p);
      unsigned n=np;
      while(n){
        //-Look for triangle without other vertices inside.
        if(n>3){
          double trisize=DBL_MAX;
          tint3 tri=TInt3(0);
          unsigned tricp1=UINT_MAX;
          for(unsigned p=0;p<n;p++){
            const unsigned cp0=(p+n-1)%n;
            const unsigned cp1=p%n;
            const unsigned cp2=(p+1)%n;
            const tdouble3 p0=vpt[vref[cp0]];
            const tdouble3 p1=vpt[vref[cp1]];
            const tdouble3 p2=vpt[vref[cp2]];
            const double tsize=TrianglePerimeter(p0,p1,p2);
            //printf("Checks triangle_%d [%d,%d,%d] tsize:%f\n",p,vref[cp0],vref[cp1],vref[cp2],tsize);
            if(tsize<trisize){
              const tdouble3 vec=fgeo::VecUnitary(fgeo::ProductVec(p0-p1,p2-p1));
              //printf("  Checks triangle_%d vec:(%f,%f,%f)\n",p,vec.x,vec.y,vec.z);
              const bool eq0=fun::IsEqual(vec0,vec,tolerance);
              if((!mode && eq0) || (mode && !eq0)){
                tplane3d pla1,pla2,pla3;
                fgeo::TriangleNormalPlanes(p0,p1,p2,0,pla1,pla2,pla3);
                bool out=true;
                for(unsigned p2=0;p2<n-3 && out;p2++){
                  const tdouble3 pt=vpt[vref[(p2+p+2)%n]];
                  out=!fgeo::TriangleInside(pt,pla1,pla2,pla3);
                  //if(!out)printf("    Point (%g,%g,%g) in triangle[%d,%d,%d]\n",pt.x,pt.y,pt.z,vref[cp0],vref[cp1],vref[cp2]);
                }
                if(out){
                  trisize=tsize;
                  tri=TInt3(vref[cp0],vref[cp1],vref[cp2]);
                  tricp1=cp1;
                }
              }
            }
          }
          if(tricp1!=UINT_MAX){
            //printf("Removes triangle [%d,%d,%d]\n",tri.x,tri.y,tri.z);
            vtri.push_back(tri);
            for(unsigned cp=tricp1;cp+1<n;cp++)vref[cp]=vref[cp+1];
            n--;
          }
          else n=0;
        }
        else{
          vtri.push_back(TInt3(vref[0],vref[1],vref[2]));
          n=0;
        }
      }
    }
    //-Select complete division (using all points) with lower area.
    unsigned ntri=unsigned(vtri.size());
    if(ntri>ntri0){
      const unsigned ntri1=ntri-ntri0;
      if(ntri0!=ntri1){//-Select complete division (using more triangles).
        if(ntri0>ntri1)vtri.erase(vtri.begin()+ntri0,vtri.begin()+ntri);
        else vtri.erase(vtri.begin(),vtri.begin()+ntri0);
      }
      else{//-Select division with lower area.
        double area0=0;
        for(unsigned t=0;t<ntri0;t++){
          const tint3 tt=vtri[t];
          area0+=TriangleArea(vpt[tt.x],vpt[tt.y],vpt[tt.z]);
        }
        double area1=0;
        const unsigned ntri=unsigned(vtri.size());
        for(unsigned t=ntri0;t<ntri;t++){
          const tint3 tt=vtri[t];
          area1+=TriangleArea(vpt[tt.x],vpt[tt.y],vpt[tt.z]);
        }
        //printf("--> Triangle area: %f  =?  %f\n",area0,area1);
        if(area0<=area1)vtri.erase(vtri.begin()+ntri0,vtri.begin()+ntri);
        else vtri.erase(vtri.begin(),vtri.begin()+ntri0);
      }
    }
    //-Frees memory.
    delete[] vref; vref=NULL;
  }
}


//==============================================================================
/// Devuelve true cuando el punto esta dentro de los planos.
/// Returns true when the point is inside planes.
//==============================================================================
bool PolygonInside(const tdouble3 &pt,const std::vector<tplane3d> &vpla){
  int neg=0,pos=0;
  const unsigned npa=unsigned(vpla.size());
  for(unsigned p=0;p<npa;p++){
    const double r=PlanePoint(vpla[p],pt);
    if(r<=0)neg++;
    if(r>=0)pos++;
  }
  return(neg==npa || pos==npa);
}

//==============================================================================
/// Devuelve true cuando el punto esta dentro de los planos.
/// Returns true when the point is inside planes.
//==============================================================================
bool PolygonInside(const tfloat3 &pt,const std::vector<tplane3f> &vpla){
  int neg=0,pos=0;
  const unsigned npa=unsigned(vpla.size());
  for(unsigned p=0;p<npa;p++){
    const float r=PlanePoint(vpla[p],pt);
    if(r<=0)neg++;
    if(r>=0)pos++;
  }
  return(neg==npa || pos==npa);
}


//==============================================================================
/// Devuelve true cuando el punto esta dentro de los planos.
/// Returns true when the point is inside planes.
//==============================================================================
bool PolygonInside(const tdouble3 &pt,unsigned npla,const tplane3d *vpla){
  int neg=0,pos=0;
  for(unsigned p=0;p<npla;p++){
    const double r=PlanePoint(vpla[p],pt);
    if(r<=0)neg++;
    if(r>=0)pos++;
  }
  return(neg==npla || pos==npla);
}

//==============================================================================
/// Devuelve true cuando el punto esta dentro de los planos.
/// Returns true when the point is inside planes.
//==============================================================================
bool PolygonInside(const tfloat3 &pt,unsigned npla,const tplane3f *vpla){
  int neg=0,pos=0;
  for(unsigned p=0;p<npla;p++){
    const float r=PlanePoint(vpla[p],pt);
    if(r<=0)neg++;
    if(r>=0)pos++;
  }
  return(neg==npla || pos==npla);
}













}


