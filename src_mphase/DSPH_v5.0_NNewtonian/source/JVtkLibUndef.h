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

/// \file JVtkLibUndef.h \brief Declares the empty class \ref JVtkLib.

#ifndef _JVtkLibUndef_
#define _JVtkLibUndef_

//##############################################################################
//# JVtkLib
//##############################################################################
/// \brief Saves VTK files with particle data and shapes.

#ifdef DISABLE_VTKLIB
class JVtkLib : protected JObject
{
public:
  /// Modes to define the normals (for CHRONO coupling).
  typedef enum{ NorNULL,NorOriginal,NorInvert,NorTwoFace }TpModeNormal; 

private:
  JShapeVtk *Shapes;

public:
  JVtkLib(){}
  ~JVtkLib(){}
  void Reset(){}

  /// Returns true when this feature is available.
  static bool Available(){ return(false); }

  static void RunExceptioonStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,const std::string &msg,const std::string &file=""){}

  //==============================================================================
  // Different functions to create special VTK files.
  //==============================================================================
  static void SaveVtkData(std::string fname,const JDataArrays &arrays,std::string posfield,bool createpath=true){}
  static void SaveVtkCells(const std::string &fname,const tfloat3 &posmin,const tuint3 &cells,float scell,bool createpath=true){}
  static void SaveVtkBoxes(const std::string &fname,unsigned nbox,const tfloat3 *vbox,float sizemin=0,bool createpath=true){}
  static void SaveVtkBoxes(const std::string &fname,unsigned nbox,const tdouble3 *vbox,float sizemin=0,bool createpath=true){}

  //==============================================================================
  // Functions to create VTK files with shapes.
  //==============================================================================
  void SaveShapeVtk(std::string file,std::string varname,bool createpath=true){}
  void AddShapePoint(const tfloat3 &pt,int value){}
  void AddShapePoint(const tdouble3 &pt,int value){}
  void AddShapePoints(unsigned np,const tfloat3 *vp,int value){}
  void AddShapePoints(unsigned np,const tdouble3 *vp,int value){}
  void AddShapeLine(const tfloat3  &pt1,const tfloat3  &pt2,int value){}
  void AddShapeLine(const tdouble3 &pt1,const tdouble3 &pt2,int value){}
  void AddShapePolyLine(unsigned np,const tfloat3  *vp,int value){}
  void AddShapePolyLine(unsigned np,const tdouble3 *vp,int value){}
  void AddShapeTriangle(const tfloat3 &pt1,const tfloat3 &pt2,const tfloat3 &pt3,int value){}
  void AddShapeTriangle(const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,int value){}
  void AddShapeQuad(const tfloat3  &pt1,const tfloat3  &pt2,const tfloat3  &pt3,const tfloat3  &pt4,int value){}
  void AddShapeQuad(const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,const tdouble3 &pt4,int value){}
  void AddShapeQuad(const tfloat3  &pt,const tfloat3  &vec,float  size,int value){}
  void AddShapeQuad(const tdouble3 &pt,const tdouble3 &vec,double size,int value){}
  void AddShapeQuadWire(const tfloat3  &pt1,const tfloat3  &pt2,const tfloat3  &pt3,const tfloat3  &pt4,int value){}
  void AddShapeQuadWire(const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,const tdouble3 &pt4,int value){}
  void AddShapeQuadWire(const tfloat3  &pt,const tfloat3  &vec,float  size,int value){}
  void AddShapeQuadWire(const tdouble3 &pt,const tdouble3 &vec,double size,int value){}
  void AddShapeBox(const tfloat3  &pt1,const tfloat3  &vx,const tfloat3  &vy,const tfloat3  &vz,int value){}
  void AddShapeBox(const tdouble3 &pt1,const tdouble3 &vx,const tdouble3 &vy,const tdouble3 &vz,int value){}
  void AddShapeBoxSize(const tfloat3  &pt1,const tfloat3  &size,int value){}
  void AddShapeBoxSize(const tdouble3 &pt1,const tdouble3 &size,int value){}
  void AddShapeBoxFront(const tfloat3  &p,const tfloat3  &px,const tfloat3  &pxz,const tfloat3  &pz
    ,const tfloat3  &py,const tfloat3  &pyx,const tfloat3  &pyxz,const tfloat3  &pyz,int value){}
  void AddShapeBoxFront(const tdouble3 &p,const tdouble3 &px,const tdouble3 &pxz,const tdouble3 &pz
    ,const tdouble3 &py,const tdouble3 &pyx,const tdouble3 &pyxz,const tdouble3 &pyz,int value){}
  void AddShapeSphere(const tfloat3  &p,float  radius,int nside,int value){}
  void AddShapeSphere(const tdouble3 &p,double radius,int nside,int value){}
  void AddShapeCylinder(const tfloat3  &p1,const tfloat3  &p2,float  radius,int nside,int value,unsigned maskfaceshide=0){}
  void AddShapeCylinder(const tdouble3 &p1,const tdouble3 &p2,double radius,int nside,int value,unsigned maskfaceshide=0){}
  void AddShapeCross(const tfloat3  &pt,float  radius,int value){}
  void AddShapeCross(const tdouble3 &pt,double radius,int value){}
  void AddShapeSpring(const tfloat3 &p1,const tfloat3 &p2,float restlength,float scalesize 
    ,float cornersout,float cornersin,float radius,float revlength,int nsides,int value){}
  void AddShapeSpring(const tdouble3 &p1,const tdouble3 &p2,double restlength,double scalesize 
    ,double cornersout,double cornersin,double radius,double revlength,int nsides,int value){}

  //==============================================================================
  // Functions to create OBJ files starting from VTK files (for CHRONO coupling).
  //==============================================================================
  static void* CreateMkShapes(const std::vector<std::string> &vtkfiles){ return(NULL); }
  static void DeleteMkShapes(void* ptr_vtksimple){}
  static unsigned CreateOBJsByMk(void* ptr_vtksimple,std::string filein,std::string filesout
    ,const std::vector<unsigned> &mkbounds,unsigned mkboundfirst,TpModeNormal normalmode){ return(0); }

};
#endif

#endif

