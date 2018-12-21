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

//:#############################################################################
//:# Descripcion:
//:# =============
//:# Conjunto de funciones tipicas de geometria y demas.
//:#
//:# Cambios:
//:# =========
//:# - Implementacion. (19-03-2013)
//:# - Metodos para calcular area de un triangulo. (01-04-2013)
//:# - Nuevos metodos para interpolacion lineal y bilineal. (08-05-2013)
//:# - Nuevas funciones trigonometricas. (20-08-2015)
//:# - Nuevas funcion DistLine(). (15-03-2016)
//:# - Nuevas funciones PointPlane() y PlanePtVec(). (23-03-2016)
//:# - Nuevas funciones. (05-04-2016)
//:# - Nuevas funciones para matrices. (24-01-2017)
//:# - En el calculo de la matriz inversa puedes pasarle el determinante. (08-02-2017)
//:# - Nuevas funciones IntersecPlaneLine(). (08-09-2016)
//:# - Nuevas funciones MulMatrix3x3(), TrasMatrix3x3() y RotMatrix3x3(). (29-11-2017)
//:# - Nueva funcion VecOrthogonal(). (10-08-2018)
//:# - Nueva funciones Rect3d2pt(), RectPosX(), RectPosY(), RectPosZ(). (21-08-2018)
//:# - Nuevas funciones VecOrthogonal2(). (05-10-2018)
//:#############################################################################

/// \file FunctionsMath.h \brief Declares basic/general math functions.

#ifndef _FunctionsMath_
#define _FunctionsMath_

#include "TypesDef.h"
#include <cstdlib>
#include <cmath>
#include <cfloat>

/// Implements a set of basic/general math functions.
namespace fmath{

//==============================================================================
/// Devuelve la interpolacion lineal de dos valores.
/// Returns the linear interpolation value.
//==============================================================================
inline double InterpolationLinear(double x,double x0,double x1,double v0,double v1){
  const double fx=(x-x0)/(x1-x0);
  return(fx*(v1-v0)+v0);
}

//==============================================================================
/// Devuelve la interpolacion lineal de dos valores.
/// Returns the linear interpolation value.
//==============================================================================
inline float InterpolationLinear(float x,float x0,float x1,float v0,float v1){
  const float fx=(x-x0)/(x1-x0);
  return(fx*(v1-v0)+v0);
}

//==============================================================================
/// Devuelve la interpolacion bilineal de cuatro valores que forman un cuadrado.
/// Returns the bilinear interpolation of four values that form a square.
//==============================================================================
inline double InterpolationBilinear(double x,double y,double px,double py,double dx,double dy,double vxy,double vxyy,double vxxy,double vxxyy){
  double vy0=InterpolationLinear(x,px,px+dx,vxy,vxxy);
  double vy1=InterpolationLinear(x,px,px+dx,vxyy,vxxyy);
  return(InterpolationLinear(y,py,py+dy,vy0,vy1));
}


//==============================================================================
/// Devuelve el producto escalar de 2 vectores.
/// Returns the scalar product of two vectors.
//==============================================================================
inline double ProductScalar(tdouble3 v1,tdouble3 v2){
  return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}

//==============================================================================
/// Devuelve el producto escalar de 2 vectores.
/// Returns the scalar product of two vectors.
//==============================================================================
inline float ProductScalar(tfloat3 v1,tfloat3 v2){
  return(v1.x*v2.x + v1.y*v2.y + v1.z*v2.z);
}


//==============================================================================
/// Devuelve el producto vectorial de 2 vectores.
/// Returns the vectorial product of two vectors.
//==============================================================================
inline tdouble3 ProductVec(const tdouble3 &v1,const tdouble3 &v2){
  tdouble3 r;
  r.x=v1.y*v2.z - v1.z*v2.y;
  r.y=v1.z*v2.x - v1.x*v2.z;
  r.z=v1.x*v2.y - v1.y*v2.x;
  return(r);
}

//==============================================================================
/// Devuelve el producto vectorial de 2 vectores.
/// Returns the vectorial product of two vectors.
//==============================================================================
inline tfloat3 ProductVec(const tfloat3 &v1,const tfloat3 &v2){
  tfloat3 r;
  r.x=v1.y*v2.z - v1.z*v2.y;
  r.y=v1.z*v2.x - v1.x*v2.z;
  r.z=v1.x*v2.y - v1.y*v2.x;
  return(r);
}


//==============================================================================
/// Resuelve punto en el plano.
/// Solves point in the plane.
//==============================================================================
inline double PointPlane(const tdouble4 &pla,const tdouble3 &pt){ 
  return(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w);
}

//==============================================================================
/// Resuelve punto en el plano.
/// Solves point in the plane.
//==============================================================================
inline float PointPlane(const tfloat4 &pla,const tfloat3 &pt){ 
  return(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w);
}


//==============================================================================
/// Devuelve la distancia entre un punto y un plano con signo.
/// Returns the distance between a point and a plane with sign.
//==============================================================================
inline double DistPlaneSign(const tdouble4 &pla,const tdouble3 &pt){
  return(PointPlane(pla,pt)/sqrt(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z));
}

//==============================================================================
/// Devuelve la distancia entre un punto y un plano con signo.
/// Returns the distance between a point and a plane with sign.
//==============================================================================
inline float DistPlaneSign(const tfloat4 &pla,const tfloat3 &pt){ 
  return(PointPlane(pla,pt)/sqrt(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z));
}


//==============================================================================
/// Devuelve la distancia entre un punto y un plano.
/// Returns the distance between a point and a plane.
//==============================================================================
inline double DistPlane(const tdouble4 &pla,const tdouble3 &pt){ 
  return(fabs(DistPlaneSign(pla,pt)));
}

//==============================================================================
/// Devuelve la distancia entre un punto y un plano.
/// Returns the distance between a point and a plane.
//==============================================================================
inline float DistPlane(const tfloat4 &pla,const tfloat3 &pt){ 
  return(fabs(DistPlaneSign(pla,pt)));
}


//==============================================================================
/// Devuelve la distancia entre dos puntos.
/// Returns the distance between two points.
//==============================================================================
inline double DistPoints(const tdouble3 &p1,const tdouble3 &p2){
  const tdouble3 v=p1-p2;
  return(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
}

//==============================================================================
/// Devuelve la distancia entre dos puntos.
/// Returns the distance between two points.
//==============================================================================
inline float DistPoints(const tfloat3 &p1,const tfloat3 &p2){
  const tfloat3 v=p1-p2;
  return(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
}


//==============================================================================
/// Devuelve la distancia al (0,0,0).
/// Returns the distance from (0,0,0).
//==============================================================================
inline double DistPoint(const tdouble3 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}

//==============================================================================
/// Devuelve la distancia al (0,0,0).
/// Returns the distance from (0,0,0).
//==============================================================================
inline float DistPoint(const tfloat3 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}


//==============================================================================
/// Devuelve vector unitario del vector (0,0,0)->p1.
/// Returns a unit vector of the vector (0,0,0)->p1.
//==============================================================================
inline tdouble3 VecUnitary(const tdouble3 &p1){
  return(p1/TDouble3(DistPoint(p1)));
}

//==============================================================================
/// Devuelve vector unitario del vector (0,0,0)->p1.
/// Returns a unit vector of the vector (0,0,0)->p1.
//==============================================================================
inline tfloat3 VecUnitary(const tfloat3 &p1){
  return(p1/TFloat3(DistPoint(p1)));
}


//==============================================================================
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//==============================================================================
inline tdouble3 VecUnitarySafe(const tdouble3 &v){
  const double m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  return(m? TDouble3(v.x/m,v.y/m,v.z/m): v);
}

//==============================================================================
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//==============================================================================
inline tfloat3 VecUnitarySafe(const tfloat3 &v){
  const float m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  return(m? TFloat3(v.x/m,v.y/m,v.z/m): v);
}


//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tdouble3 VecOrthogonal2(const tdouble3 &v,double module,bool first=true){
  tdouble3 r=TDouble3(0.f,0.f,0.f);
  if(v.x)r=(first? TDouble3(-v.z,0,v.x): TDouble3(-v.y,v.x,0));       //-When a!=0 in (a,b,c) => (-b,a,0)*y+(-c,0,a)*z 
  else if(v.y)r=(first? TDouble3(0,-v.z,v.y): TDouble3(v.y,-v.x,0));  //-When b!=0 in (a,b,c) => (0,-c,b)*z+(b,-a,0)*x 
  else if(v.z)r=(first? TDouble3(v.z,0,-v.x): TDouble3(0,v.z,-v.y));  //-When z!=0 in (a,b,c) => (c,0,-a)*x+(0,c,-b)*y 
  double m=sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
  if(m){
    m=module/m;
    r=TDouble3(r.x*m,r.y*m,r.z*m);
  }
  return(r);
}

//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tfloat3 VecOrthogonal2(const tfloat3 &v,float module,bool first=true){
  tfloat3 r=TFloat3(0.f,0.f,0.f);
  if(v.x)r=(first? TFloat3(-v.z,0,v.x): TFloat3(-v.y,v.x,0));       //-When a!=0 in (a,b,c) => (-b,a,0)*y+(-c,0,a)*z 
  else if(v.y)r=(first? TFloat3(0,-v.z,v.y): TFloat3(v.y,-v.x,0));  //-When b!=0 in (a,b,c) => (0,-c,b)*z+(b,-a,0)*x 
  else if(v.z)r=(first? TFloat3(v.z,0,-v.x): TFloat3(0,v.z,-v.y));  //-When z!=0 in (a,b,c) => (c,0,-a)*x+(0,c,-b)*y 
  float m=sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
  if(m){
    m=module/m;
    r=TFloat3(r.x*m,r.y*m,r.z*m);
  }
  return(r);
}


//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tdouble3 VecOrthogonal(const tdouble3 &v,double module){
  tdouble3 r=TDouble3(0.f,0.f,0.f);
  if(v.x)r=TDouble3((-v.y-v.z)/v.x,1,1);
  else if(v.y)r=TDouble3(1,(-v.x-v.z)/v.y,1);
  else if(v.z)r=TDouble3(1,1,(-v.x-v.y)/v.z);
  double m=sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
  if(m){
    m=module/m;
    r=TDouble3(r.x*m,r.y*m,r.z*m);
  }
  return(r);
}

//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tfloat3 VecOrthogonal(const tfloat3 &v,float module){
  tfloat3 r=TFloat3(0.f,0.f,0.f);
  if(v.x)r=TFloat3((-v.y-v.z)/v.x,1,1);
  else if(v.y)r=TFloat3(1,(-v.x-v.z)/v.y,1);
  else if(v.z)r=TFloat3(1,1,(-v.x-v.y)/v.z);
  float m=sqrt(r.x*r.x+r.y*r.y+r.z*r.z);
  if(m){
    m=module/m;
    r=TFloat3(r.x*m,r.y*m,r.z*m);
  }
  return(r);
}


//==============================================================================
/// Devuelve la normal de un triangulo.
/// Returns the normal of a triangle.
//==============================================================================
inline tdouble3 NormalTriangle(const tdouble3& p1,const tdouble3& p2,const tdouble3& p3){
  return(ProductVec(p1-p2,p2-p3));
}

//==============================================================================
/// Devuelve la normal de un triangulo.
/// Returns the normal of a triangle.
//==============================================================================
inline tfloat3 NormalTriangle(const tfloat3& p1,const tfloat3& p2,const tfloat3& p3){
  return(ProductVec(p1-p2,p2-p3));
}


//==============================================================================
/// Calcula el determinante de una matriz de 3x3.
/// Returns the determinant of a 3x3 matrix.
//==============================================================================
inline double Determinant3x3(const tmatrix3d &d){
  return(d.a11 * d.a22 * d.a33 + d.a12 * d.a23 * d.a31 + d.a13 * d.a21 * d.a32 - d.a31 * d.a22 * d.a13 - d.a32 * d.a23 * d.a11 - d.a33 * d.a21 * d.a12);
}

//==============================================================================
/// Calcula el determinante de una matriz de 3x3.
/// Returns the determinant of a 3x3 matrix.
//==============================================================================
inline float Determinant3x3(const tmatrix3f &d){
  return(d.a11 * d.a22 * d.a33 + d.a12 * d.a23 * d.a31 + d.a13 * d.a21 * d.a32 - d.a31 * d.a22 * d.a13 - d.a32 * d.a23 * d.a11 - d.a33 * d.a21 * d.a12);
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//==============================================================================
inline tmatrix3f InverseMatrix3x3(const tmatrix3f &d,const float det){
  tmatrix3f inv;
  if(det){
    inv.a11= (d.a22*d.a33-d.a23*d.a32)/det;
    inv.a12=-(d.a12*d.a33-d.a13*d.a32)/det;
    inv.a13= (d.a12*d.a23-d.a13*d.a22)/det;
    inv.a21=-(d.a21*d.a33-d.a23*d.a31)/det;
    inv.a22= (d.a11*d.a33-d.a13*d.a31)/det;
    inv.a23=-(d.a11*d.a23-d.a13*d.a21)/det;
    inv.a31= (d.a21*d.a32-d.a22*d.a31)/det;
    inv.a32=-(d.a11*d.a32-d.a12*d.a31)/det;
    inv.a33= (d.a11*d.a22-d.a12*d.a21)/det;
  }
  else inv=TMatrix3f(0);
  return(inv);
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//==============================================================================
inline tmatrix3f InverseMatrix3x3(const tmatrix3f &d){
  return(InverseMatrix3x3(d,Determinant3x3(d)));
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//==============================================================================
inline tmatrix3d InverseMatrix3x3(const tmatrix3d &d,const double det){
  tmatrix3d inv;
  if(det){
    inv.a11= (d.a22*d.a33-d.a23*d.a32)/det;
    inv.a12=-(d.a12*d.a33-d.a13*d.a32)/det;
    inv.a13= (d.a12*d.a23-d.a13*d.a22)/det;
    inv.a21=-(d.a21*d.a33-d.a23*d.a31)/det;
    inv.a22= (d.a11*d.a33-d.a13*d.a31)/det;
    inv.a23=-(d.a11*d.a23-d.a13*d.a21)/det;
    inv.a31= (d.a21*d.a32-d.a22*d.a31)/det;
    inv.a32=-(d.a11*d.a32-d.a12*d.a31)/det;
    inv.a33= (d.a11*d.a22-d.a12*d.a21)/det;
  }
  else inv=TMatrix3d(0);
  return(inv);
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//==============================================================================
inline tmatrix3d InverseMatrix3x3(const tmatrix3d &d){
  return(InverseMatrix3x3(d,Determinant3x3(d)));
}

//==============================================================================
/// Calcula el determinante de una matriz de 4x4.
/// Returns the determinant of a 4x4 matrix.
//==============================================================================
inline double Determinant4x4(const tmatrix4d &d){
  return(d.a14*d.a23*d.a32*d.a41 - d.a13*d.a24*d.a32*d.a41-
         d.a14*d.a22*d.a33*d.a41 + d.a12*d.a24*d.a33*d.a41+
         d.a13*d.a22*d.a34*d.a41 - d.a12*d.a23*d.a34*d.a41-
         d.a14*d.a23*d.a31*d.a42 + d.a13*d.a24*d.a31*d.a42+
         d.a14*d.a21*d.a33*d.a42 - d.a11*d.a24*d.a33*d.a42-
         d.a13*d.a21*d.a34*d.a42 + d.a11*d.a23*d.a34*d.a42+
         d.a14*d.a22*d.a31*d.a43 - d.a12*d.a24*d.a31*d.a43-
         d.a14*d.a21*d.a32*d.a43 + d.a11*d.a24*d.a32*d.a43+
         d.a12*d.a21*d.a34*d.a43 - d.a11*d.a22*d.a34*d.a43-
         d.a13*d.a22*d.a31*d.a44 + d.a12*d.a23*d.a31*d.a44+
         d.a13*d.a21*d.a32*d.a44 - d.a11*d.a23*d.a32*d.a44-
         d.a12*d.a21*d.a33*d.a44 + d.a11*d.a22*d.a33*d.a44);
}

//==============================================================================
/// Calcula el determinante de una matriz de 4x4.
/// Returns the determinant of a 4x4 matrix.
//==============================================================================
inline float Determinant4x4(const tmatrix4f &d){
  return(d.a14*d.a23*d.a32*d.a41 - d.a13*d.a24*d.a32*d.a41-
         d.a14*d.a22*d.a33*d.a41 + d.a12*d.a24*d.a33*d.a41+
         d.a13*d.a22*d.a34*d.a41 - d.a12*d.a23*d.a34*d.a41-
         d.a14*d.a23*d.a31*d.a42 + d.a13*d.a24*d.a31*d.a42+
         d.a14*d.a21*d.a33*d.a42 - d.a11*d.a24*d.a33*d.a42-
         d.a13*d.a21*d.a34*d.a42 + d.a11*d.a23*d.a34*d.a42+
         d.a14*d.a22*d.a31*d.a43 - d.a12*d.a24*d.a31*d.a43-
         d.a14*d.a21*d.a32*d.a43 + d.a11*d.a24*d.a32*d.a43+
         d.a12*d.a21*d.a34*d.a43 - d.a11*d.a22*d.a34*d.a43-
         d.a13*d.a22*d.a31*d.a44 + d.a12*d.a23*d.a31*d.a44+
         d.a13*d.a21*d.a32*d.a44 - d.a11*d.a23*d.a32*d.a44-
         d.a12*d.a21*d.a33*d.a44 + d.a11*d.a22*d.a33*d.a44);
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 4x4.
/// Returns the inverse matrix of a 4x4 matrix.
//==============================================================================
inline tmatrix4f InverseMatrix4x4(const tmatrix4f &d,const float det){
  tmatrix4f inv;
  if(det){
    inv.a11=(d.a22*(d.a33*d.a44-d.a34*d.a43)+d.a23*(d.a34*d.a42-d.a32*d.a44)+d.a24*(d.a32*d.a43-d.a33*d.a42))/det;
    inv.a21=(d.a21*(d.a34*d.a43-d.a33*d.a44)+d.a23*(d.a31*d.a44-d.a34*d.a41)+d.a24*(d.a33*d.a41-d.a31*d.a43))/det;
    inv.a31=(d.a21*(d.a32*d.a44-d.a34*d.a42)+d.a22*(d.a34*d.a41-d.a31*d.a44)+d.a24*(d.a31*d.a42-d.a32*d.a41))/det;
    inv.a41=(d.a21*(d.a33*d.a42-d.a32*d.a43)+d.a22*(d.a31*d.a43-d.a33*d.a41)+d.a23*(d.a32*d.a41-d.a31*d.a42))/det;
    inv.a12=(d.a12*(d.a34*d.a43-d.a33*d.a44)+d.a13*(d.a32*d.a44-d.a34*d.a42)+d.a14*(d.a33*d.a42-d.a32*d.a43))/det;
    inv.a22=(d.a11*(d.a33*d.a44-d.a34*d.a43)+d.a13*(d.a34*d.a41-d.a31*d.a44)+d.a14*(d.a31*d.a43-d.a33*d.a41))/det;
    inv.a32=(d.a11*(d.a34*d.a42-d.a32*d.a44)+d.a12*(d.a31*d.a44-d.a34*d.a41)+d.a14*(d.a32*d.a41-d.a31*d.a42))/det;
    inv.a42=(d.a11*(d.a32*d.a43-d.a33*d.a42)+d.a12*(d.a33*d.a41-d.a31*d.a43)+d.a13*(d.a31*d.a42-d.a32*d.a41))/det;
    inv.a13=(d.a12*(d.a23*d.a44-d.a24*d.a43)+d.a13*(d.a24*d.a42-d.a22*d.a44)+d.a14*(d.a22*d.a43-d.a23*d.a42))/det;
    inv.a23=(d.a11*(d.a24*d.a43-d.a23*d.a44)+d.a13*(d.a21*d.a44-d.a24*d.a41)+d.a14*(d.a23*d.a41-d.a21*d.a43))/det;
    inv.a33=(d.a11*(d.a22*d.a44-d.a24*d.a42)+d.a12*(d.a24*d.a41-d.a21*d.a44)+d.a14*(d.a21*d.a42-d.a22*d.a41))/det;
    inv.a43=(d.a11*(d.a23*d.a42-d.a22*d.a43)+d.a12*(d.a21*d.a43-d.a23*d.a41)+d.a13*(d.a22*d.a41-d.a21*d.a42))/det;
    inv.a14=(d.a12*(d.a24*d.a33-d.a23*d.a34)+d.a13*(d.a22*d.a34-d.a24*d.a32)+d.a14*(d.a23*d.a32-d.a22*d.a33))/det;
    inv.a24=(d.a11*(d.a23*d.a34-d.a24*d.a33)+d.a13*(d.a24*d.a31-d.a21*d.a34)+d.a14*(d.a21*d.a33-d.a23*d.a31))/det;
    inv.a34=(d.a11*(d.a24*d.a32-d.a22*d.a34)+d.a12*(d.a21*d.a34-d.a24*d.a31)+d.a14*(d.a22*d.a31-d.a21*d.a32))/det;
    inv.a44=(d.a11*(d.a22*d.a33-d.a23*d.a32)+d.a12*(d.a23*d.a31-d.a21*d.a33)+d.a13*(d.a21*d.a32-d.a22*d.a31))/det;
  }
  else inv=TMatrix4f(0);
  return(inv);
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 4x4.
/// Returns the inverse matrix of a 4x4 matrix.
//==============================================================================
inline tmatrix4f InverseMatrix4x4(const tmatrix4f &d){
  return(InverseMatrix4x4(d,Determinant4x4(d)));
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 4x4.
/// Returns the inverse matrix of a 4x4 matrix.
//==============================================================================
inline tmatrix4d InverseMatrix4x4(const tmatrix4d &d,const double det){
  tmatrix4d inv;
  if(det){
    inv.a11=(d.a22*(d.a33*d.a44-d.a34*d.a43)+d.a23*(d.a34*d.a42-d.a32*d.a44)+d.a24*(d.a32*d.a43-d.a33*d.a42))/det;
    inv.a21=(d.a21*(d.a34*d.a43-d.a33*d.a44)+d.a23*(d.a31*d.a44-d.a34*d.a41)+d.a24*(d.a33*d.a41-d.a31*d.a43))/det;
    inv.a31=(d.a21*(d.a32*d.a44-d.a34*d.a42)+d.a22*(d.a34*d.a41-d.a31*d.a44)+d.a24*(d.a31*d.a42-d.a32*d.a41))/det;
    inv.a41=(d.a21*(d.a33*d.a42-d.a32*d.a43)+d.a22*(d.a31*d.a43-d.a33*d.a41)+d.a23*(d.a32*d.a41-d.a31*d.a42))/det;
    inv.a12=(d.a12*(d.a34*d.a43-d.a33*d.a44)+d.a13*(d.a32*d.a44-d.a34*d.a42)+d.a14*(d.a33*d.a42-d.a32*d.a43))/det;
    inv.a22=(d.a11*(d.a33*d.a44-d.a34*d.a43)+d.a13*(d.a34*d.a41-d.a31*d.a44)+d.a14*(d.a31*d.a43-d.a33*d.a41))/det;
    inv.a32=(d.a11*(d.a34*d.a42-d.a32*d.a44)+d.a12*(d.a31*d.a44-d.a34*d.a41)+d.a14*(d.a32*d.a41-d.a31*d.a42))/det;
    inv.a42=(d.a11*(d.a32*d.a43-d.a33*d.a42)+d.a12*(d.a33*d.a41-d.a31*d.a43)+d.a13*(d.a31*d.a42-d.a32*d.a41))/det;
    inv.a13=(d.a12*(d.a23*d.a44-d.a24*d.a43)+d.a13*(d.a24*d.a42-d.a22*d.a44)+d.a14*(d.a22*d.a43-d.a23*d.a42))/det;
    inv.a23=(d.a11*(d.a24*d.a43-d.a23*d.a44)+d.a13*(d.a21*d.a44-d.a24*d.a41)+d.a14*(d.a23*d.a41-d.a21*d.a43))/det;
    inv.a33=(d.a11*(d.a22*d.a44-d.a24*d.a42)+d.a12*(d.a24*d.a41-d.a21*d.a44)+d.a14*(d.a21*d.a42-d.a22*d.a41))/det;
    inv.a43=(d.a11*(d.a23*d.a42-d.a22*d.a43)+d.a12*(d.a21*d.a43-d.a23*d.a41)+d.a13*(d.a22*d.a41-d.a21*d.a42))/det;
    inv.a14=(d.a12*(d.a24*d.a33-d.a23*d.a34)+d.a13*(d.a22*d.a34-d.a24*d.a32)+d.a14*(d.a23*d.a32-d.a22*d.a33))/det;
    inv.a24=(d.a11*(d.a23*d.a34-d.a24*d.a33)+d.a13*(d.a24*d.a31-d.a21*d.a34)+d.a14*(d.a21*d.a33-d.a23*d.a31))/det;
    inv.a34=(d.a11*(d.a24*d.a32-d.a22*d.a34)+d.a12*(d.a21*d.a34-d.a24*d.a31)+d.a14*(d.a22*d.a31-d.a21*d.a32))/det;
    inv.a44=(d.a11*(d.a22*d.a33-d.a23*d.a32)+d.a12*(d.a23*d.a31-d.a21*d.a33)+d.a13*(d.a21*d.a32-d.a22*d.a31))/det;
  }
  else inv=TMatrix4d(0);
  return(inv);
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 4x4.
/// Returns the inverse matrix of a 4x4 matrix.
//==============================================================================
inline tmatrix4d InverseMatrix4x4(const tmatrix4d &d){
  return(InverseMatrix4x4(d,Determinant4x4(d)));
}


//==============================================================================
/// Devuelve producto de 2 matrices de 3x3.
/// Returns the product of 2 matrices of 3x3.
//==============================================================================
inline tmatrix3f MulMatrix3x3(const tmatrix3f &a,const tmatrix3f &b){
  return(TMatrix3f(
    a.a11*b.a11 + a.a12*b.a21 + a.a13*b.a31, a.a11*b.a12 + a.a12*b.a22 + a.a13*b.a32, a.a11*b.a13 + a.a12*b.a23 + a.a13*b.a33,
    a.a21*b.a11 + a.a22*b.a21 + a.a23*b.a31, a.a21*b.a12 + a.a22*b.a22 + a.a23*b.a32, a.a21*b.a13 + a.a22*b.a23 + a.a23*b.a33,
    a.a31*b.a11 + a.a32*b.a21 + a.a33*b.a31, a.a31*b.a12 + a.a32*b.a22 + a.a33*b.a32, a.a31*b.a13 + a.a32*b.a23 + a.a33*b.a33
  ));
}

//==============================================================================
/// Devuelve traspuesta de matriz 3x3.
/// Returns the transpose from matrix 3x3.
//==============================================================================
inline tmatrix3f TrasMatrix3x3(const tmatrix3f &a){
  return(TMatrix3f(
    a.a11, a.a21, a.a31,
    a.a12, a.a22, a.a32,
    a.a13, a.a23, a.a33
  ));
}

//==============================================================================
/// Devuelve la matriz de rotacion.
/// Returns the rotation matrix.
//==============================================================================
inline tmatrix3f RotMatrix3x3(const tfloat3 &ang){
  const float cosx=cos(ang.x),cosy=cos(ang.y),cosz=cos(ang.z);
  const float sinx=sin(ang.x),siny=sin(ang.y),sinz=sin(ang.z);
  return(TMatrix3f(
     cosy*cosz,                   -cosy*sinz,                    siny,
     sinx*siny*cosz + cosx*sinz,  -sinx*siny*sinz + cosx*cosz,  -sinx*cosy,
    -cosx*siny*cosz + sinx*sinz,   cosx*siny*sinz + sinx*cosz,   cosx*cosy
  ));
}


//==============================================================================
/// Devuelve proyeccion ortogonal del punto en el plano.
/// Returns orthogonal projection of the point in the plane.
//==============================================================================
inline tdouble3 PtOrthogonal(const tdouble3 &pt,const tdouble4 &pla){
  const double t=-(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w)/(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z);
  return(TDouble3(pt.x+pla.x*t,pt.y+pla.y*t,pt.z+pla.z*t));
}

//==============================================================================
/// Devuelve proyeccion ortogonal del punto en el plano.
/// Returns orthogonal projection of the point in the plane.
//==============================================================================
inline tfloat3 PtOrthogonal(const tfloat3 &pt,const tfloat4 &pla){
  const float t=-(pla.x*pt.x+pla.y*pt.y+pla.z*pt.z+pla.w)/(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z);
  return(TFloat3(pt.x+pla.x*t,pt.y+pla.y*t,pt.z+pla.z*t));
}


//==============================================================================
/// Devuelve el plano formado por 3 puntos.
/// Returns the plane defined by 3 points.
//==============================================================================
tdouble4 Plane3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3);

//==============================================================================
/// Devuelve el plano formado por 3 puntos.
/// Returns the plane defined by 3 points.
//==============================================================================
tfloat4 Plane3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3);


//==============================================================================
/// Devuelve el plano formado por un punto y un vector.
/// Returns the plane defined by a point and a vector.
//==============================================================================
inline tdouble4 PlanePtVec(const tdouble3 &pt,const tdouble3 &vec){
  const tdouble3 v=VecUnitary(vec);//-No es necesario pero asi el modulo del vector no afecta al resultado de PointPlane().
  return(TDouble4(v.x,v.y,v.z,-v.x*pt.x-v.y*pt.y-v.z*pt.z));
}

//==============================================================================
/// Devuelve el plano formado por un punto y un vector.
/// Returns the plane defined by a point and a vector.
//==============================================================================
inline tfloat4 PlanePtVec(const tfloat3 &pt,const tfloat3 &vec){
  const tfloat3 v=VecUnitary(vec);//-No es necesario pero asi el modulo del vector no afecta al resultado de PointPlane().
  return(TFloat4(v.x,v.y,v.z,-v.x*pt.x-v.y*pt.y-v.z*pt.z));
}


//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal planes.
//==============================================================================
void NormalPlanes3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3,double openingdist,tdouble4 &pla1,tdouble4 &pla2,tdouble4 &pla3);

//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Los calculos internos se hacen con double precision.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal levels.
/// The internal computation is performed with double precision.
//==============================================================================
inline void NormalPlanes3Pt_dbl(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat4 &pla1,tfloat4 &pla2,tfloat4 &pla3){
  tdouble4 plad1,plad2,plad3;
  NormalPlanes3Pt(ToTDouble3(p1),ToTDouble3(p2),ToTDouble3(p3),double(openingdist),plad1,plad2,plad3);
  pla1=ToTFloat4(plad1); pla2=ToTFloat4(plad2); pla3=ToTFloat4(plad3);
}

//==============================================================================
/// Devuelve los tres planos normales que limitan un triangulo formado por 3 puntos.
/// Con openingdist puedes abrir o cerrar los planos normales.
/// Returns the three normal planes which bound a triangle formed by 3 points.
/// With openingdist you can open or close normal planes.
//==============================================================================
void NormalPlanes3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat4 &pla1,tfloat4 &pla2,tfloat4 &pla3);


//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si.
/// Returns intersection of three planes not parallel to each other.
//==============================================================================
tdouble3 Intersec3Planes(const tdouble4 &pla1,const tdouble4 &pla2,const tdouble4 &pla3);

//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si.
/// Returns intersection of three planes not parallel to each other.
//==============================================================================
tfloat3 Intersec3Planes(const tfloat4 &pla1,const tfloat4 &pla2,const tfloat4 &pla3);


//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tdouble3 IntersecPlaneLine(const tdouble4 &pla,const tdouble3 &pt1,const tdouble3 &pt2);

//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tfloat3 IntersecPlaneLine(const tfloat4 &pla,const tfloat3 &pt1,const tfloat3 &pt2);


//==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
/// Starting from a triangle formed by 3 points returns the points that form
/// a triangle more or less open according to openingdist.
//==============================================================================
void OpenTriangle3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3,double openingdist,tdouble3 &pt1,tdouble3 &pt2,tdouble3 &pt3);

//==============================================================================
/// A partir de un triangulo formado por 3 puntos devuelve los puntos que forman
/// un triangulo mas o menos abierto segun openingdist.
/// Starting from a triangle formed by 3 points returns the points that form
/// a triangle more or less open according to openingdist.
//==============================================================================
void OpenTriangle3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3,float openingdist,tfloat3 &pt1,tfloat3 &pt2,tfloat3 &pt3);

//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
/// Returns the area of a triangle formed by 3 points.
//==============================================================================
double AreaTriangle(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3);

//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
/// Returns the area of a triangle formed by 3 points.
//==============================================================================
float AreaTriangle(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3);


//==============================================================================
/// Devuelve la distancia entre un punto y una recta entre dos puntos.
/// Returns the distance between a point and a line between two points.
//==============================================================================
inline double DistLine(const tdouble3 &pt,const tdouble3 &pr1,const tdouble3 &pr2){
  double ar=AreaTriangle(pt,pr1,pr2);
  double dis=DistPoints(pr1,pr2);
  return((ar*2)/dis);
}

//==============================================================================
/// Devuelve la distancia entre un punto y una recta entre dos puntos.
/// Returns the distance between a point and a line between two points.
//==============================================================================
inline float DistLine(const tfloat3 &pt,const tfloat3 &pr1,const tfloat3 &pr2){
  float ar=AreaTriangle(pt,pr1,pr2);
  float dis=DistPoints(pr1,pr2);
  return((ar*2)/dis);
}


//==============================================================================
/// Devuelve el angulo en grados que forman dos vectores.
/// Returns angle in degrees between two vectors.
//==============================================================================
inline double AngleVector(const tdouble3 &v1,const tdouble3 &v2){
  return(acos(ProductScalar(v1,v2)/(DistPoint(v1)*DistPoint(v2)))*TODEG);
}

//==============================================================================
/// Devuelve el angulo en grados que forman dos vectores.
/// Returns angle in degrees between two vectors.
//==============================================================================
inline float AngleVector(const tfloat3 &v1,const tfloat3 &v2){
  return(float(acos(ProductScalar(v1,v2)/(DistPoint(v1)*DistPoint(v2)))*TODEG));
}


//==============================================================================
/// Devuelve el angulo en grados que forman dos planos.
/// Returns angle in degrees between two planes.
//==============================================================================
inline double AnglePlanes(tdouble4 v1,tdouble4 v2){
  return(AngleVector(TDouble3(v1.x,v1.y,v1.z),TDouble3(v2.x,v2.y,v2.z)));
}

//==============================================================================
/// Devuelve el angulo en grados que forman dos planos.
/// Returns angle in degrees between two planes.
//==============================================================================
inline float AnglePlanes(tfloat4 v1,tfloat4 v2){
  return(AngleVector(TFloat3(v1.x,v1.y,v1.z),TFloat3(v2.x,v2.y,v2.z)));
}


///Structure with two points to define a rect.
typedef struct{
  tdouble3 p; ///<Point of rect.
  tdouble3 v; ///<Vector of rect.
  tdouble3 p1; ///<Vector of rect.
  tdouble3 p2; ///<Vector of rect.
}StRect3d;

//==============================================================================
/// Devuelve recta definida por 2 puntos.
/// Returns rect defined by 2 points.
//==============================================================================
inline StRect3d Rect3d2pt(tdouble3 p1,tdouble3 p2){
  StRect3d r={p1,VecUnitary(p2-p1),p1,p2};
  return(r);
}

//==============================================================================
/// Devuelve posicion en la recta para un valor de X o DBL_MAX para posiciones no validas.
/// Returns position on the rect for a X value or DBL_MAX for invalid positions.
//==============================================================================
tdouble3 RectPosX(const StRect3d &r,double x);

//==============================================================================
/// Devuelve posicion en la recta para un valor de Y o DBL_MAX para posiciones no validas.
/// Returns position on the rect for a Y value or DBL_MAX for invalid positions.
//==============================================================================
tdouble3 RectPosY(const StRect3d &r,double y);

//==============================================================================
/// Devuelve posicion en la recta para un valor de Z o DBL_MAX para posiciones no validas.
/// Returns position on the rect for a Z value or DBL_MAX for invalid positions.
//==============================================================================
tdouble3 RectPosZ(const StRect3d &r,double z);


//==============================================================================
/// Devuelve normal eliminando error de precision en double.
/// Returns normal removing the error of precision in double.
//==============================================================================
inline tdouble3 CorrectNormal(tdouble3 n){
  if(fabs(n.x)<DBL_EPSILON*10)n.x=0;
  if(fabs(n.y)<DBL_EPSILON*10)n.y=0;
  if(fabs(n.z)<DBL_EPSILON*10)n.z=0;
  return(VecUnitary(n));
}

//==============================================================================
/// Devuelve normal eliminando error de precision en float.
/// Returns normal removing the error of precision in float.
//==============================================================================
inline tfloat3 CorrectNormal(tfloat3 n){
  if(fabs(n.x)<FLT_EPSILON*10)n.x=0;
  if(fabs(n.y)<FLT_EPSILON*10)n.y=0;
  if(fabs(n.z)<FLT_EPSILON*10)n.z=0;
  return(VecUnitary(n));
}


//==============================================================================
/// Returns cotangent of angle in radians.
//==============================================================================
inline double cot(double z){ return(1.0 / tan(z)); }

//==============================================================================
/// Returns hyperbolic cotangent of angle in radians.
//==============================================================================
inline double coth(double z){ return(cosh(z) / sinh(z)); }

//==============================================================================
/// Returns secant of angle in radians.
//==============================================================================
inline double sec(double z){ return(1.0 / cos(z)); }

//==============================================================================
/// Returns cosecant of input angle in radians.
//==============================================================================
inline double csc(double z){ return(1.0 / sin(z)); }

}

#endif


