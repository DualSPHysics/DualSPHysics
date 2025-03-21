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

/// \file JSpVtkShape.h \brief Declares the class \ref JSpVtkShape.

#ifndef _JSpVtkShape_
#define _JSpVtkShape_

#include "JObject.h"
#include "TypesDef.h"
#include <string>
#include <climits>
#include <vector>

//##############################################################################
//# JSpVtkShape
//##############################################################################
/// \brief Creates VTK files with shapes.

class JSpVtkShape : protected JObject
{
protected:
  typedef struct StrCell{
    byte tp;
    word np;
    word v;
  }StCell;

protected:
  unsigned NpTp[3];
  unsigned NcTp[3];
  std::vector<tfloat3> Ptos;
  std::vector<StCell> Cells;

protected:
  unsigned CellCount()const{ return(unsigned(Cells.size())); }
  unsigned PtCount()const{ return(unsigned(Ptos.size())); }
  void SavePtos(std::ofstream& pf)const;
  void SaveVerts(std::ofstream& pf)const;
  void SaveLines(std::ofstream& pf)const;
  void SavePoly(std::ofstream& pf)const;
  void SaveData(const std::string& vname,std::ofstream& pf)const;

public:
  JSpVtkShape();
  ~JSpVtkShape();
  void Reset();

  void AddPoint(const tfloat3& p1,word v=0);
  void AddPoint(const tdouble3& p1,word v=0);
  void AddPoints(unsigned npoints,const tfloat3* points,word v=0);
  void AddPoints(unsigned npoints,const tdouble3* points,word v=0);

  void AddLine(const tfloat3& p1,const tfloat3& p2,word v=0);
  void AddLine(const tdouble3& p1,const tdouble3& p2,word v=0);
  void AddLines(unsigned npoints,const tfloat3* points,word v=0);
  void AddLines(unsigned npoints,const tdouble3* points,word v=0);

  void AddTriangle(const tfloat3& p1,const tfloat3& p2,const tfloat3& p3,word v=0);
  void AddTriangle(const tdouble3& p1,const tdouble3& p2,const tdouble3& p3,word v=0);
  void AddTriangleWire(const tfloat3& p1,const tfloat3& p2,const tfloat3& p3,word v=0);
  void AddTriangleWire(const tdouble3& p1,const tdouble3& p2,const tdouble3& p3,word v=0);

  void AddQuad(const tfloat3& p1,const tfloat3& p2,const tfloat3& p3
    ,const tfloat3& p4,word v=0);
  void AddQuad(const tdouble3& p1,const tdouble3& p2,const tdouble3& p3
    ,const tdouble3& p4,word v=0);
  void AddQuadWire(const tfloat3& p1,const tfloat3& p2,const tfloat3& p3
    ,const tfloat3& p4,word v=0);
  void AddQuadWire(const tdouble3& p1,const tdouble3& p2,const tdouble3& p3
    ,const tdouble3& p4,word v=0);
  void AddQuadOrtho(const tfloat3& p1,const tfloat3& v1,float sidesize,word v=0);
  void AddQuadOrtho(const tdouble3& p1,const tdouble3& v1,double sidesize,word v=0);
  void AddQuadOrthoWire(const tfloat3& p1,const tfloat3& v1,float sidesize,word v=0);
  void AddQuadOrthoWire(const tdouble3& p1,const tdouble3& v1,double sidesize,word v=0);

  void AddPolygon(unsigned np,const tfloat3* vp,word v=0);

  void AddBoxSizeVec(const tfloat3& p0,const tfloat3& sizex
    ,const tfloat3& sizey,const tfloat3& sizez,word v=0);
  void AddBoxSizeVec(const tdouble3& p0,const tdouble3& sizex
    ,const tdouble3& sizey,const tdouble3& sizez,word v=0);
  void AddBoxSizeVecWire(const tfloat3& p0,const tfloat3& sizex
    ,const tfloat3& sizey,const tfloat3& sizez,word v=0);
  void AddBoxSizeVecWire(const tdouble3& p0,const tdouble3& sizex
    ,const tdouble3& sizey,const tdouble3& sizez,word v=0);

  void AddBoxSize(const tfloat3& p0,const tfloat3& sizexyz,word v=0);
  void AddBoxSize(const tdouble3& p0,const tdouble3& sizexyz,word v=0);
  void AddBoxSizeWire(const tfloat3& p0,const tfloat3& sizexyz,word v=0);
  void AddBoxSizeWire(const tdouble3& p0,const tdouble3& sizexyz,word v=0);

  void AddBoxFront(const tfloat3& p0,const tfloat3& px,const tfloat3& pxz
    ,const tfloat3& pz,const tfloat3& py,const tfloat3& pyx
    ,const tfloat3& pyxz,const tfloat3& pyz,word v=0);
  void AddBoxFront(const tdouble3& p0,const tdouble3& px,const tdouble3& pxz
    ,const tdouble3& pz,const tdouble3& py,const tdouble3& pyx
    ,const tdouble3& pyxz,const tdouble3& pyz,word v=0);

  void AddBoxes(unsigned nbox,const tfloat3* vbox,float sizemin);
  void AddBoxes(unsigned nbox,const tdouble3* vbox,double sizemin);

  void AddCylinder(const tfloat3& cen1,const tfloat3& cen2
    ,float radius,word v=0);
  void AddCylinder(const tdouble3& cen1,const tdouble3& cen2
    ,double radius,word v=0);
  void AddCylinderWire(const tfloat3& cen1,const tfloat3& cen2
    ,float radius,word v=0);
  void AddCylinderWire(const tdouble3& cen1,const tdouble3& cen2
    ,double radius,word v=0);

  void AddSphere(const tfloat3& pcen,float radius,word v=0);
  void AddSphere(const tdouble3& pcen,double radius,word v=0);
  void AddSphereWire(const tfloat3& pcen,float radius,word v=0);
  void AddSphereWire(const tdouble3& pcen,double radius,word v=0);

  void AddSpring(const tfloat3& point1,const tfloat3& point2
    ,float restlength,float scalesize,float radius,float revlength
    ,word v=0);
  void AddSpring(const tdouble3& point1,const tdouble3& point2
    ,double restlength,double scalesize,double radius,double revlength
    ,word v=0);

  void AddCross(const tfloat3& pcen,float size,word v=0);
  void AddCross(const tdouble3& pcen,double size,word v=0);

  void AddCircle(const tfloat3& pcen,const tfloat3& vec,float radius,word v=0);
  void AddCircle(const tdouble3& pcen,const tdouble3& vec,double radius,word v=0);


  void SaveVtk(std::string file,std::string vname="")const;
};

#endif

