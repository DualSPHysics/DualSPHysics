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

//:#############################################################################
//:# Descripcion:
//:# =============
//:# Conjunto de funciones tipicas de geometria 3D.
//:#
//:# Cambios:
//:# =========
//:# - Implementacion a partir de FunctionsMath. (08-02-2019)
//:# - Cambio de nombres:
//:#     PointPlane() ----------> PlanePoint()
//:#     DistPlaneSign() -------> PlaneDistSign()
//:#     DistPlane() -----------> PlaneDist()
//:#     DistPoints() ----------> PointsDist()
//:#     DistPoint() -----------> PointDist()
//:#     VecUnitary() ----------> none
//:#     VecUnitarySafe() ------> VecUnitary()
//:#     CorrectNormal() -------> NormalCorrect()
//:#     NormalTriangle() ------> TriangleNormal()
//:#     PtOrthogonal() --------> PlaneOrthogonalPoint()
//:#     NormalPlanes3Pt() -----> TriangleNormalPlanes()
//:#     NormalPlanes3Pt_dbl() -> TriangleNormalPlanes_dbl()
//:#     Intersec3Planes() -----> PlanesIntersec()
//:#     IntersecPlaneLine() ---> PlaneLineIntersec()  
//:#     OpenTriangle3Pt() -----> TriangleOpen()
//:#     AreaTriangle() --------> TriangleArea()
//:#     DistLine() ------------> LinePointDist()
//:#     AngleVector() ---------> VectorsAngle()
//:#     AnglePlanes() ---------> PlanesAngle()
//:#     struct StRect3d -------> struct tline3d
//:#     Rect3d2pt() -----------> TLine3d2Pt()
//:#     RectPosX() ------------> LinePointX()
//:#     RectPosY() ------------> LinePointY()
//:#     RectPosZ() ------------> LinePointZ()
//:# - Nuevas funciones TriangleInside(), PolygonNormalPlanes(), PolygonInside(). (08-02-2019)
//:# - Nuevas funciones LineOrthogonalPoint(), LineOrthogonalPointFromPr1(). (29-05-2019)
//:# - Nuevas funciones LineNearestPoint(). (04-06-2019)
//:# - Nuevas funciones PlanesDomain() y PlanesDomainCheck(). (26-11-2019)
//:# - Nuevas funciones: TrianglePerimeter(), PolygonIsConcave() y PolygonConcave(). (28-11-2019)
//:# - Nuevas funciones: PointInMinMax(). (01-12-2019)
//:# - Error corregido en PolygonConcave(). (05-12-2019)
//:# - Nuevas funciones: PlanePointsIn(). (17-12-2019)
//:# - Nuevas funciones: PlaneTriangleIntersec() y PlanePolygonIntersec(). (17-12-2019)
//:# - Nuevas funciones: PointsLower(), PointsSortLower(). (11-02-2020)
//:# - Nuevas funciones: LineMerge(). (11-02-2020)
//:# - Nuevas funciones: PlaneAxisDist(). (26-08-2020)
//:# - Nuevas funciones: PlaneNormalized(). (04-10-2020)
//:# - Nuevas funciones: PointsDist2(). (06-10-2020)
//:# - Nuevas funciones: DomainsIntersection(). (01-11-2020)
//:# - Nuevas funciones: VectorsUnitaryAngle(). (03-11-2020)
//:# - Controla error de acos() en funciones VectorsAngle() y VectorsUnitaryAngle(). (03-11-2020)
//:# - Cambio menor en PolygonConcave(). (20-11-2020)
//:# - Mueve funciones de poligonos a FunGeo3dPolygon. (21-11-2020)
//:# - Mueve funciones de triangulos a FunGeo3dTriangle. (21-11-2020)
//:# - Nuevas funciones: VecBounce(). (05-05-2021)
//:#############################################################################

/// \file FunGeo3d.h \brief Declares geometry functions for 3D.

#ifndef _FunGeo3d_
#define _FunGeo3d_

#include "TypesDef.h"
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <vector>

/// Implements a set of geometry functions for 3D.
namespace fgeo{

double TriangleArea(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3);
float TriangleArea(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3);


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
/// Devuelve true cuando pmin <= pt <= pmax.
/// Returns true when pmin <= pt <= pmax.
//==============================================================================
inline bool PointInMinMax(const tdouble3 &pt,const tdouble3 &pmin,const tdouble3 &pmax){
  return(pmin.x<=pt.x && pmin.y<=pt.y && pmin.z<=pt.z && pt.x<=pmax.x && pt.y<=pmax.x && pt.z<=pmax.z);
}

//==============================================================================
/// Devuelve true cuando pmin <= pt <= pmax.
/// Returns true when pmin <= pt <= pmax.
//==============================================================================
inline bool PointInMinMax(const tfloat3 &pt,const tfloat3 &pmin,const tfloat3 &pmax){
  return(pmin.x<=pt.x && pmin.y<=pt.y && pmin.z<=pt.z && pt.x<=pmax.x && pt.y<=pmax.x && pt.z<=pmax.z);
}


//==============================================================================
/// Devuelve true cuando parte de un dominio esta dentro de otro dominio.
/// Returns true when part of domain is inside other domain.
//==============================================================================
inline bool DomainsIntersection(const tdouble3 &pmin1,const tdouble3 &pmax1
  ,const tdouble3 &pmin2,const tdouble3 &pmax2)
{
  return(((pmin1.x<=pmin2.x && pmin2.x<=pmax1.x) || (pmin2.x<=pmin1.x && pmin1.x<=pmax2.x))
      && ((pmin1.y<=pmin2.y && pmin2.y<=pmax1.y) || (pmin2.y<=pmin1.y && pmin1.y<=pmax2.y))
      && ((pmin1.z<=pmin2.z && pmin2.z<=pmax1.z) || (pmin2.z<=pmin1.z && pmin1.z<=pmax2.z)));
}

//==============================================================================
/// Devuelve true cuando parte de un dominio esta dentro de otro dominio.
/// Returns true when part of domain is inside other domain.
//==============================================================================
inline bool DomainsIntersection(const tfloat3 &pmin1,const tfloat3 &pmax1
  ,const tfloat3 &pmin2,const tfloat3 &pmax2)
{
  return(((pmin1.x<=pmin2.x && pmin2.x<=pmax1.x) || (pmin2.x<=pmin1.x && pmin1.x<=pmax2.x))
      && ((pmin1.y<=pmin2.y && pmin2.y<=pmax1.y) || (pmin2.y<=pmin1.y && pmin1.y<=pmax2.y))
      && ((pmin1.z<=pmin2.z && pmin2.z<=pmax1.z) || (pmin2.z<=pmin1.z && pmin1.z<=pmax2.z)));
}


//==============================================================================
/// Devuelve verdadero cuando el punto a es menor que el b segun (z,y,x).
/// Returns true when point a is lower than point b according to (z,y,x).
//==============================================================================
inline bool PointsLower(const tfloat3 &a,const tfloat3 b){
  return((a.z<b.z) || (a.z==b.z && a.y<=b.y) || (a.z==b.z && a.y==b.y && a.x<=b.x));
}
//==============================================================================
/// Devuelve verdadero cuando el punto a es menor que el b segun (z,y,x).
/// Returns true when point a is lower than point b according to (z,y,x).
//==============================================================================
inline bool PointsLower(const tdouble3 &a,const tdouble3 b){
  return((a.z<b.z) || (a.z==b.z && a.y<=b.y) || (a.z==b.z && a.y==b.y && a.x<=b.x));
}


//==============================================================================
/// Reordena puntos de menor a mayor segun (z,y,x).
/// Reorders points from lowest to highest according to (z,y,x).
//==============================================================================
inline bool PointsSortLower(tfloat3 &a,tfloat3 &b){
  const bool sort=PointsLower(b,a);
  if(sort){ tfloat3 p=a; a=b; b=p; }
  return(sort);
}
//==============================================================================
/// Reordena puntos de menor a mayor segun (z,y,x).
/// Reorders points from lowest to highest according to (z,y,x).
//==============================================================================
inline bool PointsSortLower(tdouble3 &a,tdouble3 &b){
  const bool sort=PointsLower(b,a);
  if(sort){ tdouble3 p=a; a=b; b=p; }
  return(sort);
}


//==============================================================================
/// Devuelve la distancia entre dos puntos.
/// Returns the distance between two points.
//==============================================================================
inline double PointsDist(const tdouble3 &p1,const tdouble3 &p2){
  const tdouble3 v=p1-p2;
  return(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
}

//==============================================================================
/// Devuelve la distancia entre dos puntos.
/// Returns the distance between two points.
//==============================================================================
inline float PointsDist(const tfloat3 &p1,const tfloat3 &p2){
  const tfloat3 v=p1-p2;
  return(sqrt(v.x*v.x+v.y*v.y+v.z*v.z));
}


//==============================================================================
/// Devuelve la distancia al (0,0,0).
/// Returns the distance from (0,0,0).
//==============================================================================
inline double PointDist(const tdouble3 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}

//==============================================================================
/// Devuelve la distancia al (0,0,0).
/// Returns the distance from (0,0,0).
//==============================================================================
inline float PointDist(const tfloat3 &p1){
  return(sqrt(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z));
}


//==============================================================================
/// Devuelve la distancia^2 entre dos puntos.
/// Returns the distance^2 between two points.
//==============================================================================
inline double PointsDist2(const tdouble3 &p1,const tdouble3 &p2){
  const tdouble3 v=p1-p2;
  return(v.x*v.x+v.y*v.y+v.z*v.z);
}

//==============================================================================
/// Devuelve la distancia^2 entre dos puntos.
/// Returns the distance^2 between two points.
//==============================================================================
inline float PointsDist2(const tfloat3 &p1,const tfloat3 &p2){
  const tfloat3 v=p1-p2;
  return(v.x*v.x+v.y*v.y+v.z*v.z);
}


//==============================================================================
/// Devuelve la distancia^2 al (0,0,0).
/// Returns the distance^2 from (0,0,0).
//==============================================================================
inline double PointDist2(const tdouble3 &p1){
  return(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
}

//==============================================================================
/// Devuelve la distancia^2 al (0,0,0).
/// Returns the distance^2 from (0,0,0).
//==============================================================================
inline float PointDist2(const tfloat3 &p1){
  return(p1.x*p1.x+p1.y*p1.y+p1.z*p1.z);
}


//==============================================================================
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//==============================================================================
inline tdouble3 VecUnitary(const tdouble3 &v){
  const double m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  return(m? TDouble3(v.x/m,v.y/m,v.z/m): v);
}

//==============================================================================
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//==============================================================================
inline tfloat3 VecUnitary(const tfloat3 &v){
  const float m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  return(m? TFloat3(v.x/m,v.y/m,v.z/m): v);
}


//==============================================================================
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//==============================================================================
inline tdouble3 VecModule(const tdouble3 &v,double module){
  const double m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  const double m2=(m? module/m: 0);
  return(m? TDouble3(v.x*m2,v.y*m2,v.z*m2): v);
}

//==============================================================================
/// Devuelve vector unitario valido del vector or (0,0,0).
/// Returns a valid unit vector of the vector or (0,0,0).
//==============================================================================
inline tfloat3 VecModule(const tfloat3 &v,float module){
  const float m=sqrt(v.x*v.x+v.y*v.y+v.z*v.z);
  const float m2=(m? module/m: 0);
  return(m? TFloat3(v.x*m2,v.y*m2,v.z*m2): v);
}


//==============================================================================
/// Devuelve punto eliminando error de precision en double.
/// Returns point removing the error of precision in double.
//==============================================================================
inline tdouble3 PointCorrect(const tdouble3 &v,int m=10){
  return(TDouble3((fabs(v.x)<DBL_EPSILON*m? 0: v.x)
                 ,(fabs(v.y)<DBL_EPSILON*m? 0: v.y)
                 ,(fabs(v.z)<DBL_EPSILON*m? 0: v.z)));
}

//==============================================================================
/// Devuelve punto eliminando error de precision en double.
/// Returns point removing the error of precision in double.
//==============================================================================
inline tfloat3 PointCorrect(const tfloat3 &v,int m=10){
  return(TFloat3((fabs(v.x)<FLT_EPSILON*m? 0: v.x)
                ,(fabs(v.y)<FLT_EPSILON*m? 0: v.y)
                ,(fabs(v.z)<FLT_EPSILON*m? 0: v.z)));
}


//==============================================================================
/// Devuelve normal eliminando error de precision en double.
/// Returns normal removing the error of precision in double.
//==============================================================================
inline tdouble3 NormalCorrect(tdouble3 n,int m=10){
  return(VecUnitary(PointCorrect(n,m)));
}

//==============================================================================
/// Devuelve normal eliminando error de precision en float.
/// Returns normal removing the error of precision in float.
//==============================================================================
inline tfloat3 NormalCorrect(tfloat3 n,int m=10){
  return(VecUnitary(PointCorrect(n,m)));
}


//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tdouble3 VecOrthogonal2(const tdouble3 &v,double module,bool first=true){
  tdouble3 r=TDouble3(0);
  if(v.x)r=(first? TDouble3(-v.z,0,v.x): TDouble3(-v.y,v.x,0));       //-When a!=0 in (a,b,c) => (-b,a,0)*y+(-c,0,a)*z 
  else if(v.y)r=(first? TDouble3(0,-v.z,v.y): TDouble3(v.y,-v.x,0));  //-When b!=0 in (a,b,c) => (0,-c,b)*z+(b,-a,0)*x 
  else if(v.z)r=(first? TDouble3(v.z,0,-v.x): TDouble3(0,v.z,-v.y));  //-When z!=0 in (a,b,c) => (c,0,-a)*x+(0,c,-b)*y 
  return(VecModule(r,module));
}

//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tfloat3 VecOrthogonal2(const tfloat3 &v,float module,bool first=true){
  tfloat3 r=TFloat3(0);
  if(v.x)r=(first? TFloat3(-v.z,0,v.x): TFloat3(-v.y,v.x,0));       //-When a!=0 in (a,b,c) => (-b,a,0)*y+(-c,0,a)*z 
  else if(v.y)r=(first? TFloat3(0,-v.z,v.y): TFloat3(v.y,-v.x,0));  //-When b!=0 in (a,b,c) => (0,-c,b)*z+(b,-a,0)*x 
  else if(v.z)r=(first? TFloat3(v.z,0,-v.x): TFloat3(0,v.z,-v.y));  //-When z!=0 in (a,b,c) => (c,0,-a)*x+(0,c,-b)*y 
  return(VecModule(r,module));
}

//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tdouble3 VecOrthogonal(const tdouble3 &v,double module){
  tdouble3 r=TDouble3(0);
  if(v.x)r=TDouble3((-v.y-v.z)/v.x,1,1);
  else if(v.y)r=TDouble3(1,(-v.x-v.z)/v.y,1);
  else if(v.z)r=TDouble3(1,1,(-v.x-v.y)/v.z);
  return(VecModule(r,module));
}

//==============================================================================
// Devuelve un vector ortogonal al dado con el modulo indicado.
// Returns an orthogonal vector with indicated module.
//==============================================================================
inline tfloat3 VecOrthogonal(const tfloat3 &v,float module){
  tfloat3 r=TFloat3(0);
  if(v.x)r=TFloat3((-v.y-v.z)/v.x,1,1);
  else if(v.y)r=TFloat3(1,(-v.x-v.z)/v.y,1);
  else if(v.z)r=TFloat3(1,1,(-v.x-v.y)/v.z);
  return(VecModule(r,module));
}


//==============================================================================
/// Devuelve vector de rebote a una normal.
/// Returns bounce vector to a normal.
//==============================================================================
inline tfloat3 VecBounce(const tfloat3 &vec,const tfloat3 &normal){
  const tfloat3 u=normal*(fgeo::ProductScalar(vec,normal)/fgeo::ProductScalar(normal,normal));
  return(vec-u-u);
}

//==============================================================================
/// Devuelve vector de rebote a una normal.
/// Returns bounce vector to a normal.
//==============================================================================
inline tdouble3 VecBounce(const tdouble3 &vec,const tdouble3 &normal){
  const tdouble3 u=normal*(fgeo::ProductScalar(vec,normal)/fgeo::ProductScalar(normal,normal));
  return(vec-u-u);
}

//==============================================================================
/// Devuelve el angulo en grados que forman dos vectores.
/// Returns angle in degrees between two vectors.
//==============================================================================
inline double VectorsAngle(const tdouble3 &v1,const tdouble3 &v2){
  const double v=ProductScalar(v1,v2)/(PointDist(v1)*PointDist(v2));
  return(acos(v<-1.? -1.: (v>1.? 1.: v))*TODEG);
}

//==============================================================================
/// Devuelve el angulo en grados que forman dos vectores.
/// Returns angle in degrees between two vectors.
//==============================================================================
inline float VectorsAngle(const tfloat3 &v1,const tfloat3 &v2){
  const float v=ProductScalar(v1,v2)/(PointDist(v1)*PointDist(v2));
  return(float(acos(v<-1.? -1.: (v>1.? 1.: v))*TODEG));
}

//==============================================================================
/// Devuelve el angulo en grados que forman dos vectores unitarios.
/// Returns angle in degrees between two unitary vectors.
//==============================================================================
inline double VectorsUnitaryAngle(const tdouble3 &v1,const tdouble3 &v2){
  const double v=ProductScalar(v1,v2);
  return(acos(v<-1.? -1.: (v>1.? 1.: v))*TODEG);
}

//==============================================================================
/// Devuelve el angulo en grados que forman dos vectores unitarios.
/// Returns angle in degrees between two unitary vectors.
//==============================================================================
inline float VectorsUnitaryAngle(const tfloat3 &v1,const tfloat3 &v2){
  const float v=ProductScalar(v1,v2);
  return(float(acos(v<-1.? -1.: (v>1.? 1.: v))*TODEG));
}


//==============================================================================
/// Devuelve recta definida por 2 puntos.
/// Returns rect defined by 2 points.
//==============================================================================
inline tline3d TLine3d2Pt(tdouble3 p1,tdouble3 p2){
  tline3d r={p1,VecUnitary(p2-p1)};
  return(r);
}

//==============================================================================
/// Devuelve recta definida por un puntos y un vector.
/// Returns rect defined by one point and one vector.
//==============================================================================
inline tline3d TLine3dPtVec(tdouble3 pt,tdouble3 vec){
  tline3d r={pt,VecUnitary(vec)};
  return(r);
}


//==============================================================================
/// Devuelve la distancia entre un punto y una recta entre dos puntos.
/// Returns the distance between a point and a line between two points.
//==============================================================================
inline double LinePointDist(const tdouble3 &pt,const tdouble3 &pr1,const tdouble3 &pr2){
  double ar=TriangleArea(pt,pr1,pr2);
  double dis=PointsDist(pr1,pr2);
  return((ar*2)/dis);
}

//==============================================================================
/// Devuelve la distancia entre un punto y una recta entre dos puntos.
/// Returns the distance between a point and a line between two points.
//==============================================================================
inline float LinePointDist(const tfloat3 &pt,const tfloat3 &pr1,const tfloat3 &pr2){
  float ar=TriangleArea(pt,pr1,pr2);
  float dis=PointsDist(pr1,pr2);
  return((ar*2)/dis);
}


//==============================================================================
/// Devuelve proyeccion ortogonal del punto en la linea (pr1,pr2).
/// Returns orthogonal projection of the point in the line (pr1,pr2).
//==============================================================================
tdouble3 LineOrthogonalPoint(const tdouble3 &pt,const tdouble3 &pr1,const tdouble3 &pr2);

//==============================================================================
/// Devuelve proyeccion ortogonal del punto en la linea (pr1,pr2).
/// Returns orthogonal projection of the point in the line (pr1,pr2).
//==============================================================================
tfloat3 LineOrthogonalPoint(const tfloat3 &pt,const tfloat3 &pr1,const tfloat3 &pr2);


//==============================================================================
/// Devuelve el plano formado por 3 puntos.
/// Returns the plane defined by 3 points.
//==============================================================================
tplane3d Plane3Pt(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3);

//==============================================================================
/// Devuelve el plano formado por 3 puntos.
/// Returns the plane defined by 3 points.
//==============================================================================
tplane3f Plane3Pt(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3);


//==============================================================================
/// Devuelve el plano formado por un punto y un vector.
/// Returns the plane defined by a point and a vector.
//==============================================================================
inline tplane3d PlanePtVec(const tdouble3 &pt,const tdouble3 &vec){
  const tdouble3 v=VecUnitary(vec);//-No es necesario pero asi el modulo del vector no afecta al resultado de PointPlane().
  return(TPlane3d(v.x,v.y,v.z,-v.x*pt.x-v.y*pt.y-v.z*pt.z));
}

//==============================================================================
/// Devuelve el plano formado por un punto y un vector.
/// Returns the plane defined by a point and a vector.
//==============================================================================
inline tplane3f PlanePtVec(const tfloat3 &pt,const tfloat3 &vec){
  const tfloat3 v=VecUnitary(vec);//-No es necesario pero asi el modulo del vector no afecta al resultado de PointPlane().
  return(TPlane3f(v.x,v.y,v.z,-v.x*pt.x-v.y*pt.y-v.z*pt.z));
}


//==============================================================================
/// Devuelve el plano normalizado (A^2+B^2+C^2=1).
/// Returns the normalized plane (A^2+B^2+C^2=1).
//==============================================================================
inline tplane3d PlaneNormalized(const tplane3d &pla){
  const double dis=sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c);
  return(TPlane3d(pla.a/dis,pla.b/dis,pla.c/dis,pla.d/dis));
}

//==============================================================================
/// Devuelve el plano normalizado (A^2+B^2+C^2=1).
/// Returns the normalized plane (A^2+B^2+C^2=1).
//==============================================================================
inline tplane3f PlaneNormalized(const tplane3f &pla){
  const float dis=sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c);
  return(TPlane3f(pla.a/dis,pla.b/dis,pla.c/dis,pla.d/dis));
}


//==============================================================================
/// Devuelve el plano para calcular distancias desde pt en un eje arbitrario.
/// Returns the plane to calculate distance from pt on one arbitrary axis.
//==============================================================================
inline tplane3d PlaneAxisDist(const tdouble3 &pt,const tdouble3 &vec,double distdp=1){
  const tplane3d pla=PlanePtVec(pt,vec);
  const double dis=sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c)*distdp;
  return(TPlane3d(pla.a/dis,pla.b/dis,pla.c/dis,pla.d/dis));
}

//==============================================================================
/// Devuelve el plano para calcular distancias desde pt en un eje arbitrario.
/// Returns the plane to calculate distance from pt on one arbitrary axis.
//==============================================================================
inline tplane3f PlaneAxisDist(const tfloat3 &pt,const tfloat3 &vec,float distdp=1){
  const tplane3f pla=PlanePtVec(pt,vec);
  const float dis=sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c)*distdp;
  return(TPlane3f(pla.a/dis,pla.b/dis,pla.c/dis,pla.d/dis));
}


//==============================================================================
/// Resuelve punto en el plano.
/// Solves point in the plane.
//==============================================================================
inline double PlanePoint(const tplane3d &pla,const tdouble3 &pt){ 
  return(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d);
}

//==============================================================================
/// Resuelve punto en el plano.
/// Solves point in the plane.
//==============================================================================
inline float PlanePoint(const tplane3f &pla,const tfloat3 &pt){ 
  return(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d);
}


//==============================================================================
/// Devuelve la distancia entre un punto y un plano con signo.
/// Returns the distance between a point and a plane with sign.
//==============================================================================
inline double PlaneDistSign(const tplane3d &pla,const tdouble3 &pt){
  return(PlanePoint(pla,pt)/sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c));
}

//==============================================================================
/// Devuelve la distancia entre un punto y un plano con signo.
/// Returns the distance between a point and a plane with sign.
//==============================================================================
inline float PlaneDistSign(const tplane3f &pla,const tfloat3 &pt){ 
  return(PlanePoint(pla,pt)/sqrt(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c));
}


//==============================================================================
/// Devuelve la distancia entre un punto y un plano.
/// Returns the distance between a point and a plane.
//==============================================================================
inline double PlaneDist(const tplane3d &pla,const tdouble3 &pt){ 
  return(fabs(PlaneDistSign(pla,pt)));
}

//==============================================================================
/// Devuelve la distancia entre un punto y un plano.
/// Returns the distance between a point and a plane.
//==============================================================================
inline float PlaneDist(const tplane3f &pla,const tfloat3 &pt){ 
  return(fabs(PlaneDistSign(pla,pt)));
}


//==============================================================================
/// Devuelve el angulo en grados que forman dos planos.
/// Returns angle in degrees between two planes.
//==============================================================================
inline double PlanesAngle(tplane3d v1,tplane3d v2){
  return(VectorsAngle(TDouble3(v1.a,v1.b,v1.c),TDouble3(v2.a,v2.b,v2.c)));
}

//==============================================================================
/// Devuelve el angulo en grados que forman dos planos.
/// Returns angle in degrees between two planes.
//==============================================================================
inline float PlanesAngle(tplane3f v1,tplane3f v2){
  return(VectorsAngle(TFloat3(v1.a,v1.b,v1.c),TFloat3(v2.a,v2.b,v2.c)));
}


//==============================================================================
/// Devuelve proyeccion ortogonal del punto en el plano.
/// Returns orthogonal projection of the point in the plane.
//==============================================================================
inline tdouble3 PlaneOrthogonalPoint(const tdouble3 &pt,const tplane3d &pla){
  const double t=-(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d)/(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c);
  return(TDouble3(pt.x+pla.a*t,pt.y+pla.b*t,pt.z+pla.c*t));
}

//==============================================================================
/// Devuelve proyeccion ortogonal del punto en el plano.
/// Returns orthogonal projection of the point in the plane.
//==============================================================================
inline tfloat3 PlaneOrthogonalPoint(const tfloat3 &pt,const tplane3f &pla){
  const float t=-(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d)/(pla.a*pla.a+pla.b*pla.b+pla.c*pla.c);
  return(TFloat3(pt.x+pla.a*t,pt.y+pla.b*t,pt.z+pla.c*t));
}


//==============================================================================
/// Devuelve true cuando todos los puntos estan en el plano.
/// Returns true when all points are in the plane.
//==============================================================================
bool PlanePointsIn(const tplane3d &pla,unsigned np,const tdouble3 *vpt,double tolerance);

//==============================================================================
/// Devuelve true cuando todos los puntos estan en el plano.
/// Returns true when all points are in the plane.
//==============================================================================
bool PlanePointsIn(const tplane3f &pla,unsigned np,const tfloat3 *vpt,float tolerance);


//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si.
/// Returns intersection of three planes not parallel to each other.
//==============================================================================
tdouble3 PlanesIntersec(const tplane3d &pla1,const tplane3d &pla2,const tplane3d &pla3);

//==============================================================================
/// Devuelve punto de interseccion entre 3 planos no paralelos entre si.
/// Returns intersection of three planes not parallel to each other.
//==============================================================================
tfloat3 PlanesIntersec(const tplane3f &pla1,const tplane3f &pla2,const tplane3f &pla3);


//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tdouble3 PlaneLineIntersec(const tplane3d &pla,const tdouble3 &pt1,const tdouble3 &pt2);

//==============================================================================
/// Devuelve punto de interseccion entre un plano y una linea.
/// Returns intersection of a plane and a line.
//==============================================================================
tfloat3 PlaneLineIntersec(const tplane3f &pla,const tfloat3 &pt1,const tfloat3 &pt2);


//==============================================================================
/// Devuelve planos y distancias para delimitar un dominio en forma de caja.
/// Returns planes and distnaces to limit a box domain.
//==============================================================================
void PlanesDomain(const tdouble3 &pt,const tdouble3 &vx,const tdouble3 &vy
  ,const tdouble3 &vz,tplane3d &plax,tplane3d &play,tplane3d &plaz,tdouble3 &pladist);

//==============================================================================
/// Comprueba si el punto esta dentro del dominio definido.
/// Checks the point is inside the defined domain.
//==============================================================================
bool PlanesDomainCheck(const tdouble3 &pt,const tplane3d &plax,const tplane3d &play
  ,const tplane3d &plaz,const tdouble3 &pladist);


//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
/// Returns the area of a triangle formed by 3 points.
//==============================================================================
double TriangleArea(const tdouble3 &p1,const tdouble3 &p2,const tdouble3 &p3);

//==============================================================================
/// Devuelve el area de un triangulo formado por 3 puntos.
/// Returns the area of a triangle formed by 3 points.
//==============================================================================
float TriangleArea(const tfloat3 &p1,const tfloat3 &p2,const tfloat3 &p3);


}

#endif


