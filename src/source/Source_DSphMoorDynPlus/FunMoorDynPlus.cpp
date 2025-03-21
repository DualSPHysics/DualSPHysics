/*===================================================================================
<MOORDYNPLUS> Copyright (c) 2025 by Ivan Martinez-Estevez

Ivan Martinez-Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
Matt Hall (github.com/mattEhall)

This file is part of MoorDynPlus. MoorDynPlus is free software: you can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

Linking the MoorDynPlus library statically or dynamically with other modules is
making a combined work based on this library. Thus, the terms and conditions
of the GNU General Public License cover the whole combination. As a special
exception, the copyright holders of MoorDynPlus give you permission to dynamically
link this library with the program DualSPHysics to produce a combined model
featuring the capabilities of both DualSPHysics and MoorDynPlus. This exception
is strictly limited to linking between the compiled MoorDynPlus library and
DualSPHysics. It does not extend to other programs or the use of the MoorDynPlus
source code beyond the stipulations of the GPL. When the exception is used,
this paragraph must be included in the copyright notice.

MoorDynPlus is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.

You should have received a copy of the GNU General Public License along with
MoorDynPlus. If not, see <http://www.gnu.org/licenses/>.
===================================================================================*/

/// \file FunMoorDynPlus.cpp \brief Implements the common functions for MoorDynPlus.

#include "FunMoorDynPlus.h"
#include "Functions.h"
#include <cstdlib>
#include <cstring>
#include <cmath>

namespace fmdp{
  //=====================================================================
  /// 2D double array destruction functions.
  //=====================================================================
  void FreeArray2D(double** v,unsigned sizex){
    for(unsigned c=0;c<sizex;c++) delete[] v[c];
    delete[] v;
  }

  //==============================================================================
  /// 3D double array destruction functions.
  //==============================================================================
  void FreeArray3D(double*** v,unsigned sizex) {
    for(unsigned c=0;c<sizex;c++) {
      unsigned sizey=sizeof(v[c])/sizeof(double);
      for(unsigned f=0;f<sizey;f++)delete[] v[c][f];
      delete[] v[c];
    }
    delete[] v;
  }

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double** GetArray2D(const unsigned sizex,const unsigned sizey) {
    double** array2d=NULL;
    if(sizex) {
      try {
        array2d=new double*[sizex];
      }
      catch (std::bad_alloc const&){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      for(unsigned c=0;c<sizex;c++) {
        array2d[c]=new double[sizey];
        for(unsigned e=0;e<sizey;e++) {
          array2d[c][e]=0;
        }
      }
    }
    return array2d;
  }

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double*** GetArray3D(const unsigned sizex,const unsigned sizey, const unsigned sizez) {
    double*** array3d=NULL;
    if(sizex) {
      try {
        array3d=new double**[sizex];
      }
      catch (std::bad_alloc const&){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      for(unsigned c=0;c<sizex;c++) {  
        array3d[c]=new double*[sizey];
        for(unsigned e=0;e<sizey;e++) {
          array3d[c][e]=new double[sizez];
          for(unsigned f=0; f<sizey; f++) { array3d[c][e][f]=0; }
        }
      }
    }
    return array3d;
  }

  //==============================================================================
  /// Return a Pointer to Pointer of int
  //==============================================================================
  int* GetArrayInt2D(const unsigned size) {
    int* pToInt=NULL;
    if(size) {
      try {
        pToInt=new int[size];
      }
      catch (std::bad_alloc const&){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      memset(pToInt,0,sizeof(int)*size);
    }
    return pToInt;
  }

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double* GetArray(const unsigned size) {
    double* array=NULL;
    if(size) {
      try {
        array=new double[size];
      }
      catch (std::bad_alloc const&){ fun::Run_ExceptioonFun(fun::PrintStr("Could not allocate the requested memory.")); }
      memset(array,0,sizeof(double)*size);
    }
    return array;
  }

  //==============================================================================
  /// Check for Not-A-Number
  //==============================================================================
  bool CheckForNan(const double* v,const unsigned size) {
    bool err=false;
    for(unsigned i=0; i<size && !err; i++)if(fun::IsNAN(v[i]))err=true;
    return (err);
  }

  //==============================================================================
  /// Convert rad/s to Hz
  //==============================================================================
  double RadToHertz(const double r){
    return r*(1/(2*PI));
  }

  //==============================================================================
  /// Function that returns a row in a 3x3 matrix
  //==============================================================================
  inline tdouble3 GetRowMatrix3x3(const unsigned& row,const tmatrix3d& mat){
    tdouble3 rowv=TDouble3(0);
    if(row==0)rowv=TDouble3(mat.a11,mat.a12,mat.a13);
    if(row==1)rowv=TDouble3(mat.a21,mat.a22,mat.a23);
    if(row==2)rowv=TDouble3(mat.a31,mat.a32,mat.a33);
    return rowv;
  }

  //==============================================================================
  /// Update the states of the nodes
  //==============================================================================
  void UpdateStates(const unsigned offset,const double* X,double* Xd,const tmatrix3d& mat,const tdouble3& fnet) {
    // calculate RHS constant (premultiplying force vector by inverse of mass matrix  ... i.e. rhs=S*Forces) 
    for(unsigned i=0;i<3;i++) {
      const tdouble3 RHSiI=GetRowMatrix3x3(i,mat)*fnet;
      Xd[offset+i]=X[i]; // dxdt=V  (velocities)
      Xd[       i]=RHSiI.x+RHSiI.y+RHSiI.z; // dVdt=RHS*A  (accelerations)
    }
  }
};