//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file FunctionsMath.cpp \brief Implements basic/general math functions.

#include "FunctionsMath.h"
#include <cstdio>
#include <cstdlib>
#include <cmath>

namespace fmath{

//==============================================================================
/// Devuelve el plano formado por 3 puntos. La normal es (a,b,c)
/// Returns the plane defined by 3 points. The normal is (a,b,c)
//==============================================================================
tdouble4 Plane3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3){ 
  tdouble4 plano={0,0,0,0};
  if(p1!=p2 && p1!=p3 && p2!=p3){
    const tdouble3 v1=TDouble3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
    const tdouble3 v2=TDouble3(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
    const tdouble3 v=ProductVec(v1,v2); //-Vector normal del plano. //-Normal plane vector
    plano.x=v.x; plano.y=v.y; plano.z=v.z; 
    plano.w=-((v.x*p1.x)+(v.y*p1.y)+(v.z*p1.z));
  }
  return(plano);
}

//==============================================================================
/// Devuelve el plano formado por 3 puntos. La normal es (a,b,c)
/// Returns the plane defined by 3 points. The normal is (a,b,c)
//==============================================================================
tfloat4 Plane3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3){ 
  tfloat4 plano={0,0,0,0};
  if(p1!=p2 && p1!=p3 && p2!=p3){
    const tfloat3 v1=TFloat3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z);
    const tfloat3 v2=TFloat3(p3.x-p1.x,p3.y-p1.y,p3.z-p1.z);
    const tfloat3 v=ProductVec(v1,v2); //-Vector normal del plano. //-Normal plane vector
    plano.x=v.x; plano.y=v.y; plano.z=v.z; 
    plano.w=-((v.x*p1.x)+(v.y*p1.y)+(v.z*p1.z));
  }
  return(plano);
}


//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal planes.
//==============================================================================
void NormalPlanes3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3,double openingdist,tdouble4 &pla1,tdouble4 &pla2,tdouble4 &pla3){
  //-Calcula normal unitaria.
  //-Computes unit normal.
  const tdouble3 nor=VecUnitary(NormalTriangle(p1,p2,p3));
  //-Calcula planos normales a triangulo (p1,p2,p3).
  //-Computes normal triangle planes (p1,p2,p3).
  const tdouble3 p1n=p1+nor;
  const tdouble3 p2n=p2+nor;
  tdouble4 plad1=Plane3Pt(p1,p2,p1n);
  tdouble4 plad2=Plane3Pt(p2,p3,p2n);
  tdouble4 plad3=Plane3Pt(p3,p1,p1n);
  //-Abre planos normales.
  //-Opens normal planes.
  if(openingdist){
    double md=openingdist/sqrt(plad1.x*plad1.x+plad1.y*plad1.y+plad1.z*plad1.z);
    tdouble3 v=TDouble3(plad1.x*md,plad1.y*md,plad1.z*md);
    plad1=Plane3Pt(p1+v,p2+v,p1n+v);
    md=openingdist/sqrt(plad2.x*plad2.x+plad2.y*plad2.y+plad2.z*plad2.z);
    v=TDouble3(plad2.x*md,plad2.y*md,plad2.z*md);
    plad2=Plane3Pt(p2+v,p3+v,p2n+v);
    md=openingdist/sqrt(plad3.x*plad3.x+plad3.y*plad3.y+plad3.z*plad3.z);
    v=TDouble3(plad3.x*md,plad3.y*md,plad3.z*md);
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
void NormalPlanes3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat4 &pla1,tfloat4 &pla2,tfloat4 &pla3){
  //-Calcula normal unitaria.
  //-Computes unit normal.
  const tfloat3 nor=VecUnitary(NormalTriangle(p1,p2,p3));
  //-Calcula planos normales a triangulo (p1,p2,p3).
  //-Computes normal triangle planes (p1,p2,p3).
  const tfloat3 p1n=p1+nor;
  const tfloat3 p2n=p2+nor;
  tfloat4 plad1=Plane3Pt(p1,p2,p1n);
  tfloat4 plad2=Plane3Pt(p2,p3,p2n);
  tfloat4 plad3=Plane3Pt(p3,p1,p1n);
  //-Abre planos normales.
  //-Opens normal planes.
  if(openingdist){
    float md=openingdist/sqrt(plad1.x*plad1.x+plad1.y*plad1.y+plad1.z*plad1.z);
    tfloat3 v=TFloat3(plad1.x*md,plad1.y*md,plad1.z*md);
    plad1=Plane3Pt(p1+v,p2+v,p1n+v);
    md=openingdist/sqrt(plad2.x*plad2.x+plad2.y*plad2.y+plad2.z*plad2.z);
    v=TFloat3(plad2.x*md,plad2.y*md,plad2.z*md);
    plad2=Plane3Pt(p2+v,p3+v,p2n+v);
    md=openingdist/sqrt(plad3.x*plad3.x+plad3.y*plad3.y+plad3.z*plad3.z);
    v=TFloat3(plad3.x*md,plad3.y*md,plad3.z*md);
    plad3=Plane3Pt(p3+v,p1+v,p1n+v);
  }
  //-Guarda planos normales.
  //-Stores normal planes.
  pla1=plad1; pla2=plad2; pla3=plad3;
}


//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si usando
/// la Regla de Cramer (para sistemas compatibles determinados).
/// Returns point of intersection of three non-parallel planes using
/// Cramer's rule (for determinants of compatible systems).
//==============================================================================
tdouble3 Intersec3Planes(const tdouble4 &pla1,const tdouble4 &pla2,const tdouble4 &pla3){
  tdouble3 res=TDouble3(0);
  const double dm=Determinant3x3(TMatrix3d(pla1.x,pla1.y,pla1.z,pla2.x,pla2.y,pla2.z,pla3.x,pla3.y,pla3.z));
  if(dm){
    const double dx=Determinant3x3(TMatrix3d(-pla1.w,pla1.y,pla1.z,-pla2.w,pla2.y,pla2.z,-pla3.w,pla3.y,pla3.z));
    const double dy=Determinant3x3(TMatrix3d(pla1.x,-pla1.w,pla1.z,pla2.x,-pla2.w,pla2.z,pla3.x,-pla3.w,pla3.z));
    const double dz=Determinant3x3(TMatrix3d(pla1.x,pla1.y,-pla1.w,pla2.x,pla2.y,-pla2.w,pla3.x,pla3.y,-pla3.w));
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
tfloat3 Intersec3Planes(const tfloat4 &pla1,const tfloat4 &pla2,const tfloat4 &pla3){
  tfloat3 res=TFloat3(0);
  const float dm=Determinant3x3(TMatrix3f(pla1.x,pla1.y,pla1.z,pla2.x,pla2.y,pla2.z,pla3.x,pla3.y,pla3.z));
  if(dm){
    const float dx=Determinant3x3(TMatrix3f(-pla1.w,pla1.y,pla1.z,-pla2.w,pla2.y,pla2.z,-pla3.w,pla3.y,pla3.z));
    const float dy=Determinant3x3(TMatrix3f(pla1.x,-pla1.w,pla1.z,pla2.x,-pla2.w,pla2.z,pla3.x,-pla3.w,pla3.z));
    const float dz=Determinant3x3(TMatrix3f(pla1.x,pla1.y,-pla1.w,pla2.x,pla2.y,-pla2.w,pla3.x,pla3.y,-pla3.w));
    res=TFloat3(dx/dm,dy/dm,dz/dm);
  }
  return(res);
}


//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tdouble3 IntersecPlaneLine(const tdouble4 &pla,const tdouble3 &pt1,const tdouble3 &pt2){
  const tdouble3 v=(pt2-pt1);
  const double x1=pt1.x;
  const double x2=v.x;
  const double y1=pt1.y;
  const double y2=v.y;
  const double z1=pt1.z;
  const double z2=v.z;
  const double t=-(pla.x*x1+pla.y*y1+pla.z*z1+pla.w)/(pla.x*x2+pla.y*y2+pla.z*z2);
  return(TDouble3(x1+x2*t,y1+y2*t,z1+z2*t));
}

//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tfloat3 IntersecPlaneLine(const tfloat4 &pla,const tfloat3 &pt1,const tfloat3 &pt2){
  const tfloat3 v=(pt2-pt1);
  const float x1=pt1.x;
  const float x2=v.x;
  const float y1=pt1.y;
  const float y2=v.y;
  const float z1=pt1.z;
  const float z2=v.z;
  const float t=-(pla.x*x1+pla.y*y1+pla.z*z1+pla.w)/(pla.x*x2+pla.y*y2+pla.z*z2);
  return(TFloat3(x1+x2*t,y1+y2*t,z1+z2*t));
}


//==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
/// Starting from a triangle formed by 3 points returns the points that form
/// a triangle more or less open according to openingdist.
//==============================================================================
void OpenTriangle3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3,double openingdist,tdouble3 &pt1,tdouble3 &pt2,tdouble3 &pt3){
  //-Calcula plano de triangulo.
  //-Computes plane of the triangle.
  const tdouble4 pla=Plane3Pt(p1,p2,p3);
  //-Calcula planos normales abiertos.
  //-Computes open normal planes.
  tdouble4 pla1n;
  tdouble4 pla2n;
  tdouble4 pla3n;
  NormalPlanes3Pt(p1,p2,p3,openingdist,pla1n,pla2n,pla3n);
  //-Calcula interseccion entre planos.
  //-Computes the intersection of three planes
  pt1=Intersec3Planes(pla,pla1n,pla3n);
  pt2=Intersec3Planes(pla,pla1n,pla2n);
  pt3=Intersec3Planes(pla,pla2n,pla3n);
}

//==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
/// Starting from a triangle formed by 3 points returns the points that form
/// a triangle more or less open according to openingdist.
//==============================================================================
void OpenTriangle3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat3 &pt1,tfloat3 &pt2,tfloat3 &pt3){
  //-Calcula plano de triangulo.
  //-Computes plane of the triangle.
  const tfloat4 pla=Plane3Pt(p1,p2,p3);
  //-Calcula planos normales abiertos.
  //-Computes open normal planes.
  tfloat4 pla1n;
  tfloat4 pla2n;
  tfloat4 pla3n;
  NormalPlanes3Pt(p1,p2,p3,openingdist,pla1n,pla2n,pla3n);
  //-Calcula interseccion entre planos.
  //-Computes the intersection of three planes
  pt1=Intersec3Planes(pla,pla1n,pla3n);
  pt2=Intersec3Planes(pla,pla1n,pla2n);
  pt3=Intersec3Planes(pla,pla2n,pla3n);
}

//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
/// Returns the area of a triangle formed by 3 points.
//==============================================================================
double AreaTriangle(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3){
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
float AreaTriangle(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3){
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
/// Devuelve posicion en la recta para un valor de X o DBL_MAX para posiciones no validas.
/// Returns position on the rect for a X value or DBL_MAX for invalid positions.
//==============================================================================
tdouble3 RectPosX(const StRect3d &r,double x){
  if(r.v.x==0)return(TDouble3(DBL_MAX));
  const double f=(x-r.p.x)/r.v.x;
  return(TDouble3(x,r.p1.y+f*r.v.y,r.p1.z+f*r.v.z));
}

//==============================================================================
/// Devuelve posicion en la recta para un valor de Y o DBL_MAX para posiciones no validas.
/// Returns position on the rect for a Y value or DBL_MAX for invalid positions.
//==============================================================================
tdouble3 RectPosY(const StRect3d &r,double y){
  if(r.v.y==0)return(TDouble3(DBL_MAX));
  const double f=(y-r.p.y)/r.v.y;
  return(TDouble3(r.p1.x+f*r.v.x,y,r.p1.z+f*r.v.z));
}

//==============================================================================
/// Devuelve posicion en la recta para un valor de Z o DBL_MAX para posiciones no validas.
/// Returns position on the rect for a Z value or DBL_MAX for invalid positions.
//==============================================================================
tdouble3 RectPosZ(const StRect3d &r,double z){
  if(r.v.z==0)return(TDouble3(DBL_MAX));
  const double f=(z-r.p.z)/r.v.z;
  return(TDouble3(r.p1.x+f*r.v.x,r.p1.y+f*r.v.y,z));
}






}


