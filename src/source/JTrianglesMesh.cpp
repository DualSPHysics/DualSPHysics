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

/// \file JTrianglesMesh.cpp \brief Implements the class \ref JTrianglesMesh.

#include "JTrianglesMesh.h"
#include "Functions.h"

#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstring>
#include <cmath>

using namespace std;

//##############################################################################
//# JTrianglesMesh
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JTrianglesMesh::JTrianglesMesh(){
  ClassName="JTrianglesMesh";
  Triangles=NULL;
  Verts=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JTrianglesMesh::~JTrianglesMesh(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JTrianglesMesh::Reset(){
  AllocTriangles(0,0);
  PosMin=PosMax=TFloat3(FLT_MAX);
}

//==============================================================================
/// Allocates memory for triangles and vertices.
//==============================================================================
void JTrianglesMesh::AllocTriangles(unsigned tricount,unsigned vertscount)
{
  //-Frees allocated memory.
  TriCount=VertsCount=0;
  delete[] Triangles;   Triangles=NULL;
  delete[] Verts;       Verts=NULL;
  if(tricount && vertscount){
    TriCount=tricount;
    VertsCount=vertscount;
    try{
      Triangles=new tuint3[TriCount];
      Verts=new tfloat3[VertsCount];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon(fun::PrintStr(
        "Could not allocate the requested memory (triangles=%u, nverts=%u)."
        ,TriCount,VertsCount));
    }
  }
}

//==============================================================================
/// Replace basic data pointers.
//==============================================================================
void JTrianglesMesh::LoadData(unsigned ntri,const tuint3* vtri,unsigned nverts
  ,const tfloat3* verts)
{
  Reset();
  AllocTriangles(ntri,nverts);
  memcpy(Triangles,vtri ,sizeof(tuint3 )*ntri);
  memcpy(Verts    ,verts,sizeof(tfloat3)*nverts);
}

//==============================================================================
/// Loads input triangles from vectors.
//==============================================================================
void JTrianglesMesh::LoadData(const std::vector<tuint3>& vtri
  ,const std::vector<tfloat3>& verts)
{
  const unsigned ntri  =unsigned(vtri .size());
  const unsigned nverts=unsigned(verts.size());
  LoadData(ntri,vtri.data(),nverts,verts.data());
 }

//==============================================================================
/// Applies displacement transformation.
//==============================================================================
void JTrianglesMesh::TransformMove(const tfloat3& move){
  for(unsigned c=0;c<VertsCount;c++)Verts[c]=Verts[c]+move;
  //-Reset calculated info.
  PosMin=PosMax=TFloat3(FLT_MAX);
}

//==============================================================================
/// Applies matrix transformation.
//==============================================================================
void JTrianglesMesh::TransformMatrix(const tmatrix4d& m){
  //-Implements JMatrix4d::MulPoint()
  if(m!=TMatrix4d())for(unsigned c=0;c<VertsCount;c++){
    const tdouble3 ps=ToTDouble3(Verts[c]);
    Verts[c]=TFloat3(
       float(m.a11*ps.x + m.a12*ps.y + m.a13*ps.z + m.a14)
      ,float(m.a21*ps.x + m.a22*ps.y + m.a23*ps.z + m.a24)
      ,float(m.a31*ps.x + m.a32*ps.y + m.a33*ps.z + m.a34)
    );
  }
  //-Reset calculated info.
  PosMin=PosMax=TFloat3(FLT_MAX);
}

//==============================================================================
/// Inverts normals of triangles.
//==============================================================================
void JTrianglesMesh::InvertNormals(){
  for(unsigned ct=0;ct<TriCount;ct++){
    const tuint3 t=Triangles[ct];
    Triangles[ct]=TUint3(t.x,t.z,t.y);
  }
}

//==============================================================================
/// Compute domain limits of vertices.
//==============================================================================
void JTrianglesMesh::ComputeDomain(){
  PosMin=PosMax=TFloat3(FLT_MAX);
  if(VertsCount){
    //-Compute limits of triangles vertices.
    tfloat3 pmin=Verts[0],pmax=Verts[0];
    for(unsigned c=1;c<VertsCount;c++){
      pmin=MinValues(pmin,Verts[c]);
      pmax=MaxValues(pmax,Verts[c]);
    }
    PosMin=pmin; PosMax=pmax;
  }
}

