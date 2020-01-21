//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2016, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para gestionar la creacion de los puntos inlet. (25-01-2017)
//:#############################################################################

/// \file JSphInOutPoints.h \brief Declares the class \ref JSphInOutPoints.

#ifndef _JSphInOutPoints_
#define _JSphInOutPoints_

#include <string>
#include <vector>
#include "JObject.h"
#include "JMatrix4.h"
#include "DualSphDef.h"

#define DBG_INOUT_PTINIT 0   ///<JSphInOut: Saves VTK files (CfgInOut_PtInit.vtk and CfgInOut_PtInitZ.vtk) with initial inout points (0/1).
#define DBG_INOUT_PARTINIT 0 ///<JSphInOut: Saves VTK files (CfgInOut_InletIni_XXXX.vtk) with initial inout particles (0/1).


class JSphMk;
class JSphMkBlock;
class JXml;
class TiXmlElement;
class JLog2;
class JSphPartsInit;

//##############################################################################
//# XML format in JSphInOut_fmt.xml.
//##############################################################################

//##############################################################################
//# JSphInOutPoints
//##############################################################################
/// \brief Defines the position of inlet points.
class JSphInOutPoints : protected JObject
{
public:
  /// Defines structure to store particle information. 
  typedef struct StrPartData{
    const JSphMk* mkinfo;
    unsigned np;
    const tdouble3 *pos;
    const typecode *code;

    StrPartData(){ Clear(); }
    StrPartData(const JSphMk* vmkinfo,unsigned vnp,const tdouble3 *vpos,const typecode *vcode){
      mkinfo=vmkinfo; np=vnp; pos=vpos; code=vcode;
    }
    void Clear(){ 
      mkinfo=NULL; np=0; pos=NULL; code=NULL;
    }
  }StPartData;

private:
  JLog2 *Log;

  //-Basic simulation parameters.
  const bool Simulate2D;        ///<Indicates 2D simulation.
  const double Simulate2DPosY;  ///<Y value in 2D simulations.
  const byte Layers;            ///<Number of inlet particle layers.
  const double Dp;              ///<Distance between particles.
  const tdouble3 MapRealPosMin;
  const tdouble3 MapRealPosMax;
  const double InitialMove;     ///<Initial movement applied to points according Direction.

  std::vector<std::string> ConfigInfo;

  tdouble3 Direction;    ///<Inflow direction.

  unsigned Size;      ///<Size of allocated memory for Points[].
  unsigned Count;     ///<Number of valid points.
  tdouble3* Points;   ///<Position of points [Size].

  //-Domain data. PtDom[8] is the reference point in inout plane.
  tdouble3 PtDom[10]; 

  void ResizeMemory(unsigned newnpt);
  JMatrix4d ReadRotate2D(JXml *sxml,TiXmlElement* ele,const tdouble3 &pt);
  JMatrix4d ReadRotate3D(JXml *sxml,TiXmlElement* ele);

  tdouble3 DirectionFromStr(const std::string &strdir)const;
  std::string CheckParticlesDirection(const JSphMkBlock *pmk,const tdouble3 &dir)const;


  void Create2d3d_Particles(JXml *sxml,TiXmlElement* ele,const JSphPartsInit *partsdata);
  void Create2d_Line(JXml *sxml,TiXmlElement* ele);
  void Create3d_Box(JXml *sxml,TiXmlElement* ele);
  void Create3d_Circle(JXml *sxml,TiXmlElement* ele);
  void CheckPoints(const std::string &xmlrow);
  void ComputeDomainFromPoints();

public:
  JSphInOutPoints(JLog2 *log,bool simulate2d,double simulate2dposy,byte layers
    ,double dp,double initialmove,tdouble3 posmin,tdouble3 posmax);
  ~JSphInOutPoints();
  void Reset();
  void ResetPoints();

  void CreatePoints(JXml *sxml,TiXmlElement* ele,const JSphPartsInit *partsdata);

  void GetConfig(std::vector<std::string> &lines)const;

  tdouble3 GetDirection()const{ return(Direction); }
  unsigned GetCount()const{ return(Count); }
  tdouble3* GetPoints()const{ return(Points); }
  unsigned GetCountZmax(float zsurf)const;

  const tdouble3* GetPtDomain()const{ return(PtDom); };

};


#endif


