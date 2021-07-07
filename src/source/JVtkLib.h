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

//#############################################################################
//# Cambios:
//# =========
//# - Clase para generacion de ficheros VTK de particulas o formas. (10-12-2019)
//# - Permite compilar sin libreria de VTK. (13-12-2019)
//# - Nuevas funciones AddShapePolyLine(). (23-12-2019)
//# - Parametro creatpath que por defecto es true. (27-12-2019)
//# - Nuevas funciones AddShapePoint() y AddShapePoints(). (04-08-2020)
//# - Nuevas funciones AddShapeTriangle(). (26-08-2020)
//# - La funcion CreateOBJsByMk() devuelve en numero de faces creadas. (v5.0.158 / 18-10-2020)
//#############################################################################

/// \file JVtkLib.h \brief Declares the class \ref JVtkLib.

#ifndef _JVtkLib_
#define _JVtkLib_

#include "TypesDef.h"
#include "JObject.h"
#include "JDataArrays.h"
#include <string>
#include <cstring>
#include <string>
#include <vector>
#include "JVtkLibDef.h"      //Defines DISABLE_VTKLIB to compile without VTK library and more options.

class JShapeVtk;

//##############################################################################
//# JVtkLib
//##############################################################################
/// \brief Saves VTK files with particle data and shapes.

#ifdef DISABLE_VTKLIB
#include "JVtkLibUndef.h"
#else
class JVtkLib : protected JObject
{
public:
  /// Modes to define the normals (for CHRONO coupling).
  typedef enum{ NorNULL,NorOriginal,NorInvert,NorTwoFace }TpModeNormal; 

private:
  JShapeVtk *Shapes;

public:
  JVtkLib();
  ~JVtkLib();
  void Reset();

  /// Returns true when this feature is available.
  static bool Available(){ return(true); }

  //==============================================================================
  /// Throws exception related to a file from a static method.
  //==============================================================================
  static void RunExceptioonStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,const std::string &msg,const std::string &file="");


  //==============================================================================
  // Different functions to create special VTK files.
  //==============================================================================
  /// Stores particle data in VTK file.
  static void SaveVtkData(std::string fname,const JDataArrays &arrays,std::string posfield,bool createpath=true);

  /// Generates a VTK file with map cells.
  static void SaveVtkCells(const std::string &fname,const tfloat3 &posmin,const tuint3 &cells,float scell,bool createpath=true);

  /// Generates a VTK file with boxes.
  static void SaveVtkBoxes(const std::string &fname,unsigned nbox,const tfloat3 *vbox,float sizemin=0,bool createpath=true);
  /// Generates a VTK file with boxes.
  static void SaveVtkBoxes(const std::string &fname,unsigned nbox,const tdouble3 *vbox,float sizemin=0,bool createpath=true);

  //==============================================================================
  // Functions to create VTK files with shapes.
  //==============================================================================
  /// Generates a VTK file with shapes.
  void SaveShapeVtk(std::string file,std::string varname,bool createpath=true);

  /// Adds shape point.
  void AddShapePoint(const tfloat3 &pt,int value);
  /// Adds shape point.
  void AddShapePoint(const tdouble3 &pt,int value);

  /// Adds shape set of points.
  void AddShapePoints(unsigned np,const tfloat3 *vp,int value);
  /// Adds shape set of points.
  void AddShapePoints(unsigned np,const tdouble3 *vp,int value);

  /// Adds shape line.
  void AddShapeLine(const tfloat3  &pt1,const tfloat3  &pt2,int value);
  /// Adds shape line.
  void AddShapeLine(const tdouble3 &pt1,const tdouble3 &pt2,int value);

  /// Adds shape polyline.
  void AddShapePolyLine(unsigned np,const tfloat3  *vp,int value);
  /// Adds shape polyline.
  void AddShapePolyLine(unsigned np,const tdouble3 *vp,int value);

  /// Adds shape triangle using 3 points.
  void AddShapeTriangle(const tfloat3 &pt1,const tfloat3 &pt2,const tfloat3 &pt3,int value);
  /// Adds shape triangle using 3 points.
  void AddShapeTriangle(const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,int value);

  /// Adds shape quad using 4 points.
  void AddShapeQuad(const tfloat3  &pt1,const tfloat3  &pt2,const tfloat3  &pt3,const tfloat3  &pt4,int value);
  /// Adds shape quad using 4 points.
  void AddShapeQuad(const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,const tdouble3 &pt4,int value);

  /// Adds shape quad using an orthogonal vector.
  void AddShapeQuad(const tfloat3  &pt,const tfloat3  &vec,float  size,int value);
  /// Adds shape quad using an orthogonal vector.
  void AddShapeQuad(const tdouble3 &pt,const tdouble3 &vec,double size,int value);

  /// Adds shape lines of quad using 4 points.
  void AddShapeQuadWire(const tfloat3  &pt1,const tfloat3  &pt2,const tfloat3  &pt3,const tfloat3  &pt4,int value);
  /// Adds shape lines of quad using 4 points.
  void AddShapeQuadWire(const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,const tdouble3 &pt4,int value);

  /// Adds shape lines of quad using an orthogonal vector.
  void AddShapeQuadWire(const tfloat3  &pt,const tfloat3  &vec,float  size,int value);
  /// Adds shape lines of quad using an orthogonal vector.
  void AddShapeQuadWire(const tdouble3 &pt,const tdouble3 &vec,double size,int value);

  /// Adds shape box (1pt + 3x vec3).
  void AddShapeBox(const tfloat3  &pt1,const tfloat3  &vx,const tfloat3  &vy,const tfloat3  &vz,int value);
  /// Adds shape box (1pt + 3x vec3).
  void AddShapeBox(const tdouble3 &pt1,const tdouble3 &vx,const tdouble3 &vy,const tdouble3 &vz,int value);

  /// Adds shape box (1pt + size).
  void AddShapeBoxSize(const tfloat3  &pt1,const tfloat3  &size,int value);
  /// Adds shape box (1pt + 3x vec3).
  void AddShapeBoxSize(const tdouble3 &pt1,const tdouble3 &size,int value);

  /// Adds shape box (4pt front + 4pt back).
  void AddShapeBoxFront(const tfloat3  &p,const tfloat3  &px,const tfloat3  &pxz,const tfloat3  &pz
    ,const tfloat3  &py,const tfloat3  &pyx,const tfloat3  &pyxz,const tfloat3  &pyz,int value);
  /// Adds shape box (4pt front + 4pt back).
  void AddShapeBoxFront(const tdouble3 &p,const tdouble3 &px,const tdouble3 &pxz,const tdouble3 &pz
    ,const tdouble3 &py,const tdouble3 &pyx,const tdouble3 &pyxz,const tdouble3 &pyz,int value);

  /// Adds shape sphere using quads.
  void AddShapeSphere(const tfloat3  &p,float  radius,int nside,int value);
  /// Adds shape sphere using quads.
  void AddShapeSphere(const tdouble3 &p,double radius,int nside,int value);

  /// Adds shape cylinder.
  void AddShapeCylinder(const tfloat3  &p1,const tfloat3  &p2,float  radius,int nside,int value,unsigned maskfaceshide=0);
  /// Adds shape cylinder.
  void AddShapeCylinder(const tdouble3 &p1,const tdouble3 &p2,double radius,int nside,int value,unsigned maskfaceshide=0);

  /// Adds lines to create a cross.
  void AddShapeCross(const tfloat3  &pt,float  radius,int value);
  /// Adds lines to create a cross.
  void AddShapeCross(const tdouble3 &pt,double radius,int value);

  /// Adds spring using lines.
  /// \param cornersout Size of corner.
  /// \param cornersin  Size of corner (inside).
  /// \param radius     Spring radius.
  /// \param revlength  Length for each revolution.
  /// \param nsides     Number of sections for each revolution.
  void AddShapeSpring(const tfloat3 &p1,const tfloat3 &p2,float restlength,float scalesize 
    ,float cornersout,float cornersin,float radius,float revlength,int nsides,int value);
  /// Adds spring using lines.
  void AddShapeSpring(const tdouble3 &p1,const tdouble3 &p2,double restlength,double scalesize 
    ,double cornersout,double cornersin,double radius,double revlength,int nsides,int value);


  //==============================================================================
  // Functions to create OBJ files starting from VTK files (for CHRONO coupling).
  //==============================================================================
  /// Creates object with geometry (triangles and quads) and mk data from VTK files.
  static void* CreateMkShapes(const std::vector<std::string> &vtkfiles);

  /// Frees object with geometry and mk data from VTK files.
  static void DeleteMkShapes(void* ptr_vtksimple);

  /// Creates OBJ file with MK geometry in VTK file. Returns number of created shapes.
  static unsigned CreateOBJsByMk(void* ptr_vtksimple,std::string filein,std::string filesout
    ,const std::vector<unsigned> &mkbounds,unsigned mkboundfirst,TpModeNormal normalmode);


  //==============================================================================
  // Functions to compute normals from final particles (for mDBC and under development).
  //==============================================================================
  static void ComputeNormalsPartCells(bool data2d,double data2dposy,double dp
    ,tdouble3 mapposmin,tdouble3 mapposmax,double dist,std::string dirout
    ,unsigned nsel,const unsigned *partsel,unsigned np,const tdouble3 *pos
    ,tfloat3 *boundnormal);


};
#endif

#endif

