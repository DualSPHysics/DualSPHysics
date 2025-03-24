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

/// \file JTrianglesMesh.h \brief Declares the class \ref JTrianglesMesh.

#ifndef _JTrianglesMesh_
#define _JTrianglesMesh_

#include "JObject.h"
#include "TypesDef.h"
#include <string>

//##############################################################################
//# JTrianglesMesh
//##############################################################################
/// \brief Manages the information of mesh of triangles.

class JTrianglesMesh : protected JObject
{
protected:
  unsigned TriCount;    ///<Number of triangles. 
  tuint3*  Triangles;   ///<Positions of vertices. [TriCount]
  unsigned VertsCount;  ///<Number of vertices.
  tfloat3* Verts;       ///<Positions of vertices [VertsCount].

  tfloat3 PosMin;      ///<Minimum position of triangles domain.
  tfloat3 PosMax;      ///<Maximum position of triangles domain.

protected:
  void AllocTriangles(unsigned tricount,unsigned vertscount);

public:
  JTrianglesMesh();
  ~JTrianglesMesh();
  void Reset();
  void LoadData(unsigned ntri,const tuint3* vtri
    ,unsigned nverts,const tfloat3* verts);
  void LoadData(const std::vector<tuint3>& vtri
    ,const std::vector<tfloat3>& verts);

  void TransformMove(const tfloat3& move);
  void TransformMatrix(const tmatrix4d& mat);
  void InvertNormals();

  void ComputeDomain();

  unsigned       GetTriCount  ()const{ return(TriCount);   } 
  const tuint3*  GetTriangles ()const{ return(Triangles);  } 
  unsigned       GetVertsCount()const{ return(VertsCount); } 
  const tfloat3* GetVerts     ()const{ return(Verts);      } 

  tfloat3 GetPosMin ()const{ return(PosMin);  }
  tfloat3 GetPosMax ()const{ return(PosMax);  }

};

#endif

