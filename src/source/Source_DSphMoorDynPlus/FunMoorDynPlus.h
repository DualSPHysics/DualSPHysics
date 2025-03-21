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

/// \file FunMoorDynPlus.cpp \brief Defines the common functions for MoorDynPlus.

#ifndef _FunMoorDynPlus_
#define _FunMoorDynPlus_

#include "TypesDef.h"

#include <string>
#include <vector>
#include <sstream>

namespace fmdp{
  //==============================================================================
  /// Converts numeric values to string and returns it
  //==============================================================================
  template <class T> std::string ToString(T value){
    std::string ret="";
    std::stringstream ss;
    ss<<value;
    ret=(ss.str());
    return ret;
  }  
 
  //==============================================================================
  /// Frees memory vector
  //==============================================================================
  template <class T> void FreeVector(std::vector<T*> v){
    for(unsigned i=0;i<v.size();i++){delete v[i];}
    v.clear();
  }

  //=====================================================================
  /// 2D double array destruction functions.
  //=====================================================================
  void FreeArray2D(double** v,unsigned sizex);

  //==============================================================================
  /// 3D double array destruction functions.
  //==============================================================================
  void FreeArray3D(double*** v,unsigned sizex);

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double** GetArray2D(const unsigned sizex,const unsigned sizey);

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double*** GetArray3D(const unsigned sizex,const unsigned sizey,const unsigned sizez);

  //==============================================================================
  /// Return a Pointer to Pointer of int
  //==============================================================================
  int* GetArrayInt2D(const unsigned size);

  //==============================================================================
  /// Return a Pointer to Pointer of doubles
  //==============================================================================
  double* GetArray(const unsigned size);

  //==============================================================================
  /// Check for Not-A-Number
  //==============================================================================
  bool CheckForNan(const double* v,const unsigned size);

  //==============================================================================
  /// Convert rad/s to Hz
  //==============================================================================
  double RadToHertz(const double r);

  //==============================================================================
  /// Update the states of the nodes
  //==============================================================================
  void UpdateStates(const unsigned offset,const double* X,double* Xd,const tmatrix3d& mat,const tdouble3& fnet);
}
#endif // !_FunMoorDynPlus_
