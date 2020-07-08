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

/// \file FunctionsMath_iker.h \brief Implements basic/general math functions for the GPU executions.

#include "TypesDef.h"
#include <cuda_runtime_api.h>

namespace cumath{

//------------------------------------------------------------------------------
/// Resuelve punto en el plano.
/// Solves point in the plane.
//------------------------------------------------------------------------------
__device__ double PointPlane(const float4 &pla,const double3 &pt){ 
  return(pt.x*pla.x+pt.y*pla.y+pt.z*pla.z+pla.w);
}

//------------------------------------------------------------------------------
/// Resuelve punto en el plano.
/// Solves point in the plane.
//------------------------------------------------------------------------------
__device__ float PointPlane(const float4 &pla,float px,float py,float pz){ 
  return(pla.x*px+pla.y*py+pla.z*pz+pla.w);
}

//------------------------------------------------------------------------------
/// Returns the distance between a point and a plane.
/// Devuelve la distancia entre un punto y un plano.
//------------------------------------------------------------------------------
__device__ double DistPlaneSign(const float4 &pla,const double3 &pt){
  return(PointPlane(pla,pt)/sqrt(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z));
}

//------------------------------------------------------------------------------
/// Returns the distance between a point and a plane.
/// Devuelve la distancia entre un punto y un plano.
//------------------------------------------------------------------------------
__device__ float KerDistPlaneSign(const float4 &pla,float px,float py,float pz){ 
  return(PointPlane(pla,px,py,pz)/sqrt(pla.x*pla.x+pla.y*pla.y+pla.z*pla.z));
}

//------------------------------------------------------------------------------
/// Returns the distance between a point and a plane.
/// Devuelve la distancia entre un punto y un plano.
//------------------------------------------------------------------------------
__device__ double DistPlane(const float4 &pla,const double3 &pt){ 
  return(fabs(DistPlaneSign(pla,pt)));
}

//------------------------------------------------------------------------------
/// Initializes matrix to zero.
/// Inicializa matriz a cero.
//------------------------------------------------------------------------------
__device__ void Tmatrix3fReset(tmatrix3f &m){ 
  m.a11=m.a12=m.a13=m.a21=m.a22=m.a23=m.a31=m.a32=m.a33=0; 
}

//------------------------------------------------------------------------------
/// Initializes matrix to zero.
/// Inicializa matriz a cero.
//------------------------------------------------------------------------------
__device__ void Tmatrix3dReset(tmatrix3d &m){ 
  m.a11=m.a12=m.a13=m.a21=m.a22=m.a23=m.a31=m.a32=m.a33=0; 
}

//------------------------------------------------------------------------------
/// Initializes matrix to zero.
/// Inicializa matriz a cero.
//------------------------------------------------------------------------------
__device__ void Tmatrix4fReset(tmatrix4f &m){ 
  m.a11=m.a12=m.a13=m.a14=m.a21=m.a22=m.a23=m.a24=m.a31=m.a32=m.a33=m.a34=m.a41=m.a42=m.a43=m.a44=0; 
}

//------------------------------------------------------------------------------
/// Initializes matrix to zero.
/// Inicializa matriz a cero.
//------------------------------------------------------------------------------
__device__ void Tmatrix4dReset(tmatrix4d &m){ 
  m.a11=m.a12=m.a13=m.a14=m.a21=m.a22=m.a23=m.a24=m.a31=m.a32=m.a33=m.a34=m.a41=m.a42=m.a43=m.a44=0; 
}

//------------------------------------------------------------------------------
/// Calcula el determinante de una matriz de 3x3.
/// Returns the determinant of a 3x3 matrix.
//------------------------------------------------------------------------------
__device__ float Determinant3x3(const tmatrix3f &d){
  return(d.a11*d.a22*d.a33 + d.a12*d.a23*d.a31 + d.a13*d.a21*d.a32 - d.a31*d.a22*d.a13 - d.a32*d.a23*d.a11 - d.a33*d.a21*d.a12);
}

//------------------------------------------------------------------------------
/// Calcula el determinante de una matriz de 3x3.
/// Returns the determinant of a 3x3 matrix.
//------------------------------------------------------------------------------
__device__ double Determinant3x3dbl(const tmatrix3f &d){
  return(double(d.a11)*double(d.a22)*double(d.a33) + double(d.a12)*double(d.a23)*double(d.a31) + double(d.a13)*double(d.a21)*double(d.a32) - double(d.a31)*double(d.a22)*double(d.a13) - double(d.a32)*double(d.a23)*double(d.a11) - double(d.a33)*double(d.a21)*double(d.a12));
}

//------------------------------------------------------------------------------
/// Calcula el determinante de una matriz de 3x3.
/// Returns the determinant of a 3x3 matrix.
//------------------------------------------------------------------------------
__device__ double Determinant3x3(const tmatrix3d &d){
  return(d.a11*d.a22*d.a33 + d.a12*d.a23*d.a31 + d.a13*d.a21*d.a32 - d.a31*d.a22*d.a13 - d.a32*d.a23*d.a11 - d.a33*d.a21*d.a12);
}

//------------------------------------------------------------------------------
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//------------------------------------------------------------------------------
__device__ tmatrix3f InverseMatrix3x3(const tmatrix3f &d,const float det){
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
  else Tmatrix3fReset(inv);
  return(inv);
}

//------------------------------------------------------------------------------
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//------------------------------------------------------------------------------
__device__ tmatrix3f InverseMatrix3x3dbl(const tmatrix3f &d,const double det){
  tmatrix3f inv;
  if(det){
    inv.a11=float( (double(d.a22)*double(d.a33)-double(d.a23)*double(d.a32))/det);
    inv.a12=float(-(double(d.a12)*double(d.a33)-double(d.a13)*double(d.a32))/det);
    inv.a13=float( (double(d.a12)*double(d.a23)-double(d.a13)*double(d.a22))/det);
    inv.a21=float(-(double(d.a21)*double(d.a33)-double(d.a23)*double(d.a31))/det);
    inv.a22=float( (double(d.a11)*double(d.a33)-double(d.a13)*double(d.a31))/det);
    inv.a23=float(-(double(d.a11)*double(d.a23)-double(d.a13)*double(d.a21))/det);
    inv.a31=float( (double(d.a21)*double(d.a32)-double(d.a22)*double(d.a31))/det);
    inv.a32=float(-(double(d.a11)*double(d.a32)-double(d.a12)*double(d.a31))/det);
    inv.a33=float( (double(d.a11)*double(d.a22)-double(d.a12)*double(d.a21))/det);
  }
  else Tmatrix3fReset(inv);
  return(inv);
}

//------------------------------------------------------------------------------
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//------------------------------------------------------------------------------
__device__ tmatrix3d InverseMatrix3x3(const tmatrix3d &d,const double det){
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
  else Tmatrix3dReset(inv);
  return(inv);
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//==============================================================================
__device__ tmatrix3f InverseMatrix3x3(const tmatrix3f &d){
  return(InverseMatrix3x3(d,Determinant3x3(d)));
}

//==============================================================================
/// Devuelve la matriz inversa de una matriz de 3x3.
/// Returns the inverse matrix of a 3x3 matrix.
//==============================================================================
__device__ tmatrix3d InverseMatrix3x3(const tmatrix3d &d){
  return(InverseMatrix3x3(d,Determinant3x3(d)));
}

//------------------------------------------------------------------------------
/// Calcula el determinante de una matriz de 4x4.
/// Returns the determinant of a 4x4 matrix.
//------------------------------------------------------------------------------
__device__ float Determinant4x4(const tmatrix4f &d){
  return(d.a14*d.a23*d.a32*d.a41 - d.a13*d.a24*d.a32*d.a41 -
         d.a14*d.a22*d.a33*d.a41 + d.a12*d.a24*d.a33*d.a41 +
         d.a13*d.a22*d.a34*d.a41 - d.a12*d.a23*d.a34*d.a41 -
         d.a14*d.a23*d.a31*d.a42 + d.a13*d.a24*d.a31*d.a42 +
         d.a14*d.a21*d.a33*d.a42 - d.a11*d.a24*d.a33*d.a42 -
         d.a13*d.a21*d.a34*d.a42 + d.a11*d.a23*d.a34*d.a42 +
         d.a14*d.a22*d.a31*d.a43 - d.a12*d.a24*d.a31*d.a43 -
         d.a14*d.a21*d.a32*d.a43 + d.a11*d.a24*d.a32*d.a43 +
         d.a12*d.a21*d.a34*d.a43 - d.a11*d.a22*d.a34*d.a43 -
         d.a13*d.a22*d.a31*d.a44 + d.a12*d.a23*d.a31*d.a44 +
         d.a13*d.a21*d.a32*d.a44 - d.a11*d.a23*d.a32*d.a44 -
         d.a12*d.a21*d.a33*d.a44 + d.a11*d.a22*d.a33*d.a44);
}

//------------------------------------------------------------------------------
/// Calcula el determinante de una matriz de 4x4.
/// Returns the determinant of a 4x4 matrix.
//------------------------------------------------------------------------------
__device__ double Determinant4x4dbl(const tmatrix4f &d){
  return(double(d.a14)*double(d.a23)*double(d.a32)*double(d.a41) - double(d.a13)*double(d.a24)*double(d.a32)*double(d.a41) -
         double(d.a14)*double(d.a22)*double(d.a33)*double(d.a41) + double(d.a12)*double(d.a24)*double(d.a33)*double(d.a41) +
         double(d.a13)*double(d.a22)*double(d.a34)*double(d.a41) - double(d.a12)*double(d.a23)*double(d.a34)*double(d.a41) -
         double(d.a14)*double(d.a23)*double(d.a31)*double(d.a42) + double(d.a13)*double(d.a24)*double(d.a31)*double(d.a42) +
         double(d.a14)*double(d.a21)*double(d.a33)*double(d.a42) - double(d.a11)*double(d.a24)*double(d.a33)*double(d.a42) -
         double(d.a13)*double(d.a21)*double(d.a34)*double(d.a42) + double(d.a11)*double(d.a23)*double(d.a34)*double(d.a42) +
         double(d.a14)*double(d.a22)*double(d.a31)*double(d.a43) - double(d.a12)*double(d.a24)*double(d.a31)*double(d.a43) -
         double(d.a14)*double(d.a21)*double(d.a32)*double(d.a43) + double(d.a11)*double(d.a24)*double(d.a32)*double(d.a43) +
         double(d.a12)*double(d.a21)*double(d.a34)*double(d.a43) - double(d.a11)*double(d.a22)*double(d.a34)*double(d.a43) -
         double(d.a13)*double(d.a22)*double(d.a31)*double(d.a44) + double(d.a12)*double(d.a23)*double(d.a31)*double(d.a44) +
         double(d.a13)*double(d.a21)*double(d.a32)*double(d.a44) - double(d.a11)*double(d.a23)*double(d.a32)*double(d.a44) -
         double(d.a12)*double(d.a21)*double(d.a33)*double(d.a44) + double(d.a11)*double(d.a22)*double(d.a33)*double(d.a44));
}

//------------------------------------------------------------------------------
/// Calcula el determinante de una matriz de 4x4.
/// Returns the determinant of a 4x4 matrix.
//------------------------------------------------------------------------------
__device__ double Determinant4x4(const tmatrix4d &d){
  return(d.a14*d.a23*d.a32*d.a41 - d.a13*d.a24*d.a32*d.a41 -
         d.a14*d.a22*d.a33*d.a41 + d.a12*d.a24*d.a33*d.a41 +
         d.a13*d.a22*d.a34*d.a41 - d.a12*d.a23*d.a34*d.a41 -
         d.a14*d.a23*d.a31*d.a42 + d.a13*d.a24*d.a31*d.a42 +
         d.a14*d.a21*d.a33*d.a42 - d.a11*d.a24*d.a33*d.a42 -
         d.a13*d.a21*d.a34*d.a42 + d.a11*d.a23*d.a34*d.a42 +
         d.a14*d.a22*d.a31*d.a43 - d.a12*d.a24*d.a31*d.a43 -
         d.a14*d.a21*d.a32*d.a43 + d.a11*d.a24*d.a32*d.a43 +
         d.a12*d.a21*d.a34*d.a43 - d.a11*d.a22*d.a34*d.a43 -
         d.a13*d.a22*d.a31*d.a44 + d.a12*d.a23*d.a31*d.a44 +
         d.a13*d.a21*d.a32*d.a44 - d.a11*d.a23*d.a32*d.a44 -
         d.a12*d.a21*d.a33*d.a44 + d.a11*d.a22*d.a33*d.a44);
}

//------------------------------------------------------------------------------
/// Devuelve la matriz inversa de una matriz de 4x4.
/// Returns the inverse matrix of a 4x4 matrix.
//------------------------------------------------------------------------------
__device__ tmatrix4f InverseMatrix4x4(const tmatrix4f &d,const float det){
  tmatrix4f inv;
  if(det){
    inv.a11=(d.a22*(d.a33*d.a44-d.a34*d.a43) + d.a23*(d.a34*d.a42-d.a32*d.a44) + d.a24*(d.a32*d.a43-d.a33*d.a42)) /det;
    inv.a21=(d.a21*(d.a34*d.a43-d.a33*d.a44) + d.a23*(d.a31*d.a44-d.a34*d.a41) + d.a24*(d.a33*d.a41-d.a31*d.a43)) /det;
    inv.a31=(d.a21*(d.a32*d.a44-d.a34*d.a42) + d.a22*(d.a34*d.a41-d.a31*d.a44) + d.a24*(d.a31*d.a42-d.a32*d.a41)) /det;
    inv.a41=(d.a21*(d.a33*d.a42-d.a32*d.a43) + d.a22*(d.a31*d.a43-d.a33*d.a41) + d.a23*(d.a32*d.a41-d.a31*d.a42)) /det;
    inv.a12=(d.a12*(d.a34*d.a43-d.a33*d.a44) + d.a13*(d.a32*d.a44-d.a34*d.a42) + d.a14*(d.a33*d.a42-d.a32*d.a43)) /det;
    inv.a22=(d.a11*(d.a33*d.a44-d.a34*d.a43) + d.a13*(d.a34*d.a41-d.a31*d.a44) + d.a14*(d.a31*d.a43-d.a33*d.a41)) /det;
    inv.a32=(d.a11*(d.a34*d.a42-d.a32*d.a44) + d.a12*(d.a31*d.a44-d.a34*d.a41) + d.a14*(d.a32*d.a41-d.a31*d.a42)) /det;
    inv.a42=(d.a11*(d.a32*d.a43-d.a33*d.a42) + d.a12*(d.a33*d.a41-d.a31*d.a43) + d.a13*(d.a31*d.a42-d.a32*d.a41)) /det;
    inv.a13=(d.a12*(d.a23*d.a44-d.a24*d.a43) + d.a13*(d.a24*d.a42-d.a22*d.a44) + d.a14*(d.a22*d.a43-d.a23*d.a42)) /det;
    inv.a23=(d.a11*(d.a24*d.a43-d.a23*d.a44) + d.a13*(d.a21*d.a44-d.a24*d.a41) + d.a14*(d.a23*d.a41-d.a21*d.a43)) /det;
    inv.a33=(d.a11*(d.a22*d.a44-d.a24*d.a42) + d.a12*(d.a24*d.a41-d.a21*d.a44) + d.a14*(d.a21*d.a42-d.a22*d.a41)) /det;
    inv.a43=(d.a11*(d.a23*d.a42-d.a22*d.a43) + d.a12*(d.a21*d.a43-d.a23*d.a41) + d.a13*(d.a22*d.a41-d.a21*d.a42)) /det;
    inv.a14=(d.a12*(d.a24*d.a33-d.a23*d.a34) + d.a13*(d.a22*d.a34-d.a24*d.a32) + d.a14*(d.a23*d.a32-d.a22*d.a33)) /det;
    inv.a24=(d.a11*(d.a23*d.a34-d.a24*d.a33) + d.a13*(d.a24*d.a31-d.a21*d.a34) + d.a14*(d.a21*d.a33-d.a23*d.a31)) /det;
    inv.a34=(d.a11*(d.a24*d.a32-d.a22*d.a34) + d.a12*(d.a21*d.a34-d.a24*d.a31) + d.a14*(d.a22*d.a31-d.a21*d.a32)) /det;
    inv.a44=(d.a11*(d.a22*d.a33-d.a23*d.a32) + d.a12*(d.a23*d.a31-d.a21*d.a33) + d.a13*(d.a21*d.a32-d.a22*d.a31)) /det;
  }
  else Tmatrix4fReset(inv);
  return(inv);
}

//------------------------------------------------------------------------------
/// Devuelve la matriz inversa de una matriz de 4x4.
/// Returns the inverse matrix of a 4x4 matrix.
//------------------------------------------------------------------------------
__device__ tmatrix4f InverseMatrix4x4dbl(const tmatrix4f &d,const double det){
  tmatrix4f inv;
  if(det){
    inv.a11=(double(d.a22)*(double(d.a33)*double(d.a44)-double(d.a34)*double(d.a43)) + double(d.a23)*(double(d.a34)*double(d.a42)-double(d.a32)*double(d.a44)) + double(d.a24)*(double(d.a32)*double(d.a43)-double(d.a33)*double(d.a42))) /det;
    inv.a21=(double(d.a21)*(double(d.a34)*double(d.a43)-double(d.a33)*double(d.a44)) + double(d.a23)*(double(d.a31)*double(d.a44)-double(d.a34)*double(d.a41)) + double(d.a24)*(double(d.a33)*double(d.a41)-double(d.a31)*double(d.a43))) /det;
    inv.a31=(double(d.a21)*(double(d.a32)*double(d.a44)-double(d.a34)*double(d.a42)) + double(d.a22)*(double(d.a34)*double(d.a41)-double(d.a31)*double(d.a44)) + double(d.a24)*(double(d.a31)*double(d.a42)-double(d.a32)*double(d.a41))) /det;
    inv.a41=(double(d.a21)*(double(d.a33)*double(d.a42)-double(d.a32)*double(d.a43)) + double(d.a22)*(double(d.a31)*double(d.a43)-double(d.a33)*double(d.a41)) + double(d.a23)*(double(d.a32)*double(d.a41)-double(d.a31)*double(d.a42))) /det;
    inv.a12=(double(d.a12)*(double(d.a34)*double(d.a43)-double(d.a33)*double(d.a44)) + double(d.a13)*(double(d.a32)*double(d.a44)-double(d.a34)*double(d.a42)) + double(d.a14)*(double(d.a33)*double(d.a42)-double(d.a32)*double(d.a43))) /det;
    inv.a22=(double(d.a11)*(double(d.a33)*double(d.a44)-double(d.a34)*double(d.a43)) + double(d.a13)*(double(d.a34)*double(d.a41)-double(d.a31)*double(d.a44)) + double(d.a14)*(double(d.a31)*double(d.a43)-double(d.a33)*double(d.a41))) /det;
    inv.a32=(double(d.a11)*(double(d.a34)*double(d.a42)-double(d.a32)*double(d.a44)) + double(d.a12)*(double(d.a31)*double(d.a44)-double(d.a34)*double(d.a41)) + double(d.a14)*(double(d.a32)*double(d.a41)-double(d.a31)*double(d.a42))) /det;
    inv.a42=(double(d.a11)*(double(d.a32)*double(d.a43)-double(d.a33)*double(d.a42)) + double(d.a12)*(double(d.a33)*double(d.a41)-double(d.a31)*double(d.a43)) + double(d.a13)*(double(d.a31)*double(d.a42)-double(d.a32)*double(d.a41))) /det;
    inv.a13=(double(d.a12)*(double(d.a23)*double(d.a44)-double(d.a24)*double(d.a43)) + double(d.a13)*(double(d.a24)*double(d.a42)-double(d.a22)*double(d.a44)) + double(d.a14)*(double(d.a22)*double(d.a43)-double(d.a23)*double(d.a42))) /det;
    inv.a23=(double(d.a11)*(double(d.a24)*double(d.a43)-double(d.a23)*double(d.a44)) + double(d.a13)*(double(d.a21)*double(d.a44)-double(d.a24)*double(d.a41)) + double(d.a14)*(double(d.a23)*double(d.a41)-double(d.a21)*double(d.a43))) /det;
    inv.a33=(double(d.a11)*(double(d.a22)*double(d.a44)-double(d.a24)*double(d.a42)) + double(d.a12)*(double(d.a24)*double(d.a41)-double(d.a21)*double(d.a44)) + double(d.a14)*(double(d.a21)*double(d.a42)-double(d.a22)*double(d.a41))) /det;
    inv.a43=(double(d.a11)*(double(d.a23)*double(d.a42)-double(d.a22)*double(d.a43)) + double(d.a12)*(double(d.a21)*double(d.a43)-double(d.a23)*double(d.a41)) + double(d.a13)*(double(d.a22)*double(d.a41)-double(d.a21)*double(d.a42))) /det;
    inv.a14=(double(d.a12)*(double(d.a24)*double(d.a33)-double(d.a23)*double(d.a34)) + double(d.a13)*(double(d.a22)*double(d.a34)-double(d.a24)*double(d.a32)) + double(d.a14)*(double(d.a23)*double(d.a32)-double(d.a22)*double(d.a33))) /det;
    inv.a24=(double(d.a11)*(double(d.a23)*double(d.a34)-double(d.a24)*double(d.a33)) + double(d.a13)*(double(d.a24)*double(d.a31)-double(d.a21)*double(d.a34)) + double(d.a14)*(double(d.a21)*double(d.a33)-double(d.a23)*double(d.a31))) /det;
    inv.a34=(double(d.a11)*(double(d.a24)*double(d.a32)-double(d.a22)*double(d.a34)) + double(d.a12)*(double(d.a21)*double(d.a34)-double(d.a24)*double(d.a31)) + double(d.a14)*(double(d.a22)*double(d.a31)-double(d.a21)*double(d.a32))) /det;
    inv.a44=(double(d.a11)*(double(d.a22)*double(d.a33)-double(d.a23)*double(d.a32)) + double(d.a12)*(double(d.a23)*double(d.a31)-double(d.a21)*double(d.a33)) + double(d.a13)*(double(d.a21)*double(d.a32)-double(d.a22)*double(d.a31))) /det;
  }
  else Tmatrix4fReset(inv);
  return(inv);
}

//------------------------------------------------------------------------------
/// Devuelve la matriz inversa de una matriz de 4x4.
/// Returns the inverse matrix of a 4x4 matrix.
//------------------------------------------------------------------------------
__device__ tmatrix4d InverseMatrix4x4(const tmatrix4d &d,const double det){
  tmatrix4d inv;
  if(det){
    inv.a11=(d.a22*(d.a33*d.a44-d.a34*d.a43) + d.a23*(d.a34*d.a42-d.a32*d.a44) + d.a24*(d.a32*d.a43-d.a33*d.a42)) /det;
    inv.a21=(d.a21*(d.a34*d.a43-d.a33*d.a44) + d.a23*(d.a31*d.a44-d.a34*d.a41) + d.a24*(d.a33*d.a41-d.a31*d.a43)) /det;
    inv.a31=(d.a21*(d.a32*d.a44-d.a34*d.a42) + d.a22*(d.a34*d.a41-d.a31*d.a44) + d.a24*(d.a31*d.a42-d.a32*d.a41)) /det;
    inv.a41=(d.a21*(d.a33*d.a42-d.a32*d.a43) + d.a22*(d.a31*d.a43-d.a33*d.a41) + d.a23*(d.a32*d.a41-d.a31*d.a42)) /det;
    inv.a12=(d.a12*(d.a34*d.a43-d.a33*d.a44) + d.a13*(d.a32*d.a44-d.a34*d.a42) + d.a14*(d.a33*d.a42-d.a32*d.a43)) /det;
    inv.a22=(d.a11*(d.a33*d.a44-d.a34*d.a43) + d.a13*(d.a34*d.a41-d.a31*d.a44) + d.a14*(d.a31*d.a43-d.a33*d.a41)) /det;
    inv.a32=(d.a11*(d.a34*d.a42-d.a32*d.a44) + d.a12*(d.a31*d.a44-d.a34*d.a41) + d.a14*(d.a32*d.a41-d.a31*d.a42)) /det;
    inv.a42=(d.a11*(d.a32*d.a43-d.a33*d.a42) + d.a12*(d.a33*d.a41-d.a31*d.a43) + d.a13*(d.a31*d.a42-d.a32*d.a41)) /det;
    inv.a13=(d.a12*(d.a23*d.a44-d.a24*d.a43) + d.a13*(d.a24*d.a42-d.a22*d.a44) + d.a14*(d.a22*d.a43-d.a23*d.a42)) /det;
    inv.a23=(d.a11*(d.a24*d.a43-d.a23*d.a44) + d.a13*(d.a21*d.a44-d.a24*d.a41) + d.a14*(d.a23*d.a41-d.a21*d.a43)) /det;
    inv.a33=(d.a11*(d.a22*d.a44-d.a24*d.a42) + d.a12*(d.a24*d.a41-d.a21*d.a44) + d.a14*(d.a21*d.a42-d.a22*d.a41)) /det;
    inv.a43=(d.a11*(d.a23*d.a42-d.a22*d.a43) + d.a12*(d.a21*d.a43-d.a23*d.a41) + d.a13*(d.a22*d.a41-d.a21*d.a42)) /det;
    inv.a14=(d.a12*(d.a24*d.a33-d.a23*d.a34) + d.a13*(d.a22*d.a34-d.a24*d.a32) + d.a14*(d.a23*d.a32-d.a22*d.a33)) /det;
    inv.a24=(d.a11*(d.a23*d.a34-d.a24*d.a33) + d.a13*(d.a24*d.a31-d.a21*d.a34) + d.a14*(d.a21*d.a33-d.a23*d.a31)) /det;
    inv.a34=(d.a11*(d.a24*d.a32-d.a22*d.a34) + d.a12*(d.a21*d.a34-d.a24*d.a31) + d.a14*(d.a22*d.a31-d.a21*d.a32)) /det;
    inv.a44=(d.a11*(d.a22*d.a33-d.a23*d.a32) + d.a12*(d.a23*d.a31-d.a21*d.a33) + d.a13*(d.a21*d.a32-d.a22*d.a31)) /det;
  }
  else Tmatrix4dReset(inv);
  return(inv);
}

//==============================================================================
/// Devuelve producto de 2 matrices de 3x3.
/// Returns the product of 2 matrices of 3x3.
//==============================================================================
__device__ tmatrix3f MulMatrix3x3(const tmatrix3f &a,const tmatrix3f &b){
  tmatrix3f ret;
  ret.a11=a.a11*b.a11 + a.a12*b.a21 + a.a13*b.a31;
  ret.a12=a.a11*b.a12 + a.a12*b.a22 + a.a13*b.a32;
  ret.a13=a.a11*b.a13 + a.a12*b.a23 + a.a13*b.a33;
  ret.a21=a.a21*b.a11 + a.a22*b.a21 + a.a23*b.a31;
  ret.a22=a.a21*b.a12 + a.a22*b.a22 + a.a23*b.a32;
  ret.a23=a.a21*b.a13 + a.a22*b.a23 + a.a23*b.a33;
  ret.a31=a.a31*b.a11 + a.a32*b.a21 + a.a33*b.a31;
  ret.a32=a.a31*b.a12 + a.a32*b.a22 + a.a33*b.a32;
  ret.a33=a.a31*b.a13 + a.a32*b.a23 + a.a33*b.a33;
  return(ret);
}

//==============================================================================
/// Devuelve producto de 2 matrices de 3x3.
/// Returns the product of 2 matrices of 3x3.
//==============================================================================
__device__ tmatrix3d MulMatrix3x3(const tmatrix3d &a,const tmatrix3d &b){
  tmatrix3d ret;
  ret.a11=a.a11*b.a11 + a.a12*b.a21 + a.a13*b.a31;
  ret.a12=a.a11*b.a12 + a.a12*b.a22 + a.a13*b.a32;
  ret.a13=a.a11*b.a13 + a.a12*b.a23 + a.a13*b.a33;
  ret.a21=a.a21*b.a11 + a.a22*b.a21 + a.a23*b.a31;
  ret.a22=a.a21*b.a12 + a.a22*b.a22 + a.a23*b.a32;
  ret.a23=a.a21*b.a13 + a.a22*b.a23 + a.a23*b.a33;
  ret.a31=a.a31*b.a11 + a.a32*b.a21 + a.a33*b.a31;
  ret.a32=a.a31*b.a12 + a.a32*b.a22 + a.a33*b.a32;
  ret.a33=a.a31*b.a13 + a.a32*b.a23 + a.a33*b.a33;
  return(ret);
}

//==============================================================================
/// Devuelve traspuesta de matriz 3x3.
/// Returns the transpose from matrix 3x3.
//==============================================================================
__device__ tmatrix3f TrasMatrix3x3(const tmatrix3f &a){
  tmatrix3f ret;
  ret.a11=a.a11;  ret.a12=a.a21;  ret.a13=a.a31;
  ret.a21=a.a12;  ret.a22=a.a22;  ret.a23=a.a32;
  ret.a31=a.a13;  ret.a32=a.a23;  ret.a33=a.a33;
  return(ret);
}

//==============================================================================
/// Devuelve traspuesta de matriz 3x3.
/// Returns the transpose from matrix 3x3.
//==============================================================================
__device__ tmatrix3d TrasMatrix3x3(const tmatrix3d &a){
  tmatrix3d ret;
  ret.a11=a.a11;  ret.a12=a.a21;  ret.a13=a.a31;
  ret.a21=a.a12;  ret.a22=a.a22;  ret.a23=a.a32;
  ret.a31=a.a13;  ret.a32=a.a23;  ret.a33=a.a33;
  return(ret);
}

//==============================================================================
/// Devuelve la matriz de rotacion.
/// Returns the rotation matrix.
//==============================================================================
__device__ tmatrix3f RotMatrix3x3(const float3 &ang){
  const float cosx=cos(ang.x),cosy=cos(ang.y),cosz=cos(ang.z);
  const float sinx=sin(ang.x),siny=sin(ang.y),sinz=sin(ang.z);
  tmatrix3f ret;
  ret.a11= cosy*cosz;  
  ret.a12=-cosy*sinz;  
  ret.a13= siny;
  ret.a21= sinx*siny*cosz + cosx*sinz;  
  ret.a22=-sinx*siny*sinz + cosx*cosz;  
  ret.a23=-sinx*cosy;
  ret.a31=-cosx*siny*cosz + sinx*sinz;
  ret.a32= cosx*siny*sinz + sinx*cosz;
  ret.a33= cosx*cosy;
  return(ret);
}

//==============================================================================
/// Devuelve la matriz de rotacion.
/// Returns the rotation matrix.
//==============================================================================
__device__ tmatrix3d RotMatrix3x3(const double3 &ang){
  const double cosx=cos(ang.x),cosy=cos(ang.y),cosz=cos(ang.z);
  const double sinx=sin(ang.x),siny=sin(ang.y),sinz=sin(ang.z);
  tmatrix3d ret;
  ret.a11= cosy*cosz;  
  ret.a12=-cosy*sinz;  
  ret.a13= siny;
  ret.a21= sinx*siny*cosz + cosx*sinz;  
  ret.a22=-sinx*siny*sinz + cosx*cosz;  
  ret.a23=-sinx*cosy;
  ret.a31=-cosx*siny*cosz + sinx*sinz;
  ret.a32= cosx*siny*sinz + sinx*cosz;
  ret.a33= cosx*cosy;
  return(ret);
}



}


