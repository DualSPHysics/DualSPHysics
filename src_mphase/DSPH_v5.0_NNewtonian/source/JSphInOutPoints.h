//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Objeto JXml pasado como const para operaciones de lectura. (18-03-2020)  
//:#############################################################################

/// \file JSphInOutPoints.h \brief Declares the class \ref JSphInOutPoints.

#ifndef _JSphInOutPoints_
#define _JSphInOutPoints_

#include <string>
#include <vector>
#include "JObject.h"
#include "JMatrix4.h"
#include "DualSphDef.h"

class JSphMk;
class JSphMkBlock;
class JXml;
class TiXmlElement;
class JLog2;
class JDsPartsInit;

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

  unsigned Size;        ///<Size of allocated memory for Points[].
  unsigned Count;       ///<Number of valid points.
  tdouble3 *Points;     ///<Position of points [Size].
  byte     *PointsInit; ///<Indicates an initial valid (z<zsurf) point [Size].

  //-Domain data. PtDom[8] is the reference point in inout plane.
  tdouble3 PtDom[10]; 

  tdouble3 ZonePosMin;
  tdouble3 ZonePosMax;

  void ResizeMemory(unsigned newnpt);
  JMatrix4d ReadRotate2D(const JXml *sxml,TiXmlElement* ele,const tdouble3 &pt);
  JMatrix4d ReadRotate3D(const JXml *sxml,TiXmlElement* ele);

  tdouble3 DirectionFromStr(const std::string &strdir)const;
  std::string CheckParticlesDirection(const JSphMkBlock *pmk,const tdouble3 &dir)const;


  void Create2d3d_Particles(const JXml *sxml,TiXmlElement* ele,const JDsPartsInit *partsdata);
  void Create2d_Line(const JXml *sxml,TiXmlElement* ele);
  void Create3d_Box(const JXml *sxml,TiXmlElement* ele);
  void Create3d_Circle(const JXml *sxml,TiXmlElement* ele);
  void CheckPoints(const std::string &xmlrow);
  void ComputeDomainLimits(tdouble3 &posmin,tdouble3 &posmax)const;
  void ComputeDomainFromPoints();

public:
  JSphInOutPoints(bool simulate2d,double simulate2dposy,byte layers
    ,double dp,double initialmove,tdouble3 posmin,tdouble3 posmax);
  ~JSphInOutPoints();
  void Reset();
  void ResetPoints();

  void CreatePoints(const JXml *sxml,TiXmlElement* ele,const JDsPartsInit *partsdata);

  void GetConfig(std::vector<std::string> &lines)const;

  tdouble3 GetDirection()const{ return(Direction); }
  unsigned GetCount()const{ return(Count); }
  tdouble3* GetPoints()const{ return(Points); }
  byte* GetPointsInit()const{ return(PointsInit); }

  void SetPointsInit(bool active);
  unsigned CountPointsInit()const;

  const tdouble3* GetPtDomain()const{ return(PtDom); };
  void GetPtDomain(std::vector<tdouble3> &ptdom)const;

  tdouble3 GetZonePosMin()const{ return(ZonePosMin); }
  tdouble3 GetZonePosMax()const{ return(ZonePosMax); }

};


#endif


