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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Class to define rotated box for Variable Resolution. (20-05-2024)
//:# - Adds rotation with center definition. (04-06-2024)
//:#############################################################################

/// \file JBoxDef.h \brief Declares the class \ref JBoxDef.

#ifndef _JBoxDef_
#define _JBoxDef_

#include "JObject.h"
#include "TypesDef.h"
#include <string>
#include <vector>

//##############################################################################
//# JBoxDef
//##############################################################################
/// \brief Class to manage a box domain.
class JBoxDef : protected JObject 
{
protected:
  bool   Is2D;      ///< 2-D domain.
  double Posy2D;    ///< Y-position for 2-D domain.

  bool Simple;      ///< Use simple box definition (posmin & posmax).
  tdouble3 PosMin;  ///< Minimum limits of domain. 
  tdouble3 PosMax;  ///< Maximum limits of domain (Simple=true) or final point. 
  tdouble3 Vx;      ///< Vector in X with size in X (Simple=false). 
  tdouble3 Vy;      ///< Vector in Y with size in Y (Simple=false). 
  tdouble3 Vz;      ///< Vector in Z with size in Z (Simple=false). 

public:
  JBoxDef(const JBoxDef& src);
  JBoxDef(bool is2d=false,double posy2d=0);
  JBoxDef(const tdouble3& pmin,const tdouble3& pmax,bool is2d=false);
  JBoxDef(const tdouble3& pmin,const tdouble3& vx,const tdouble3& vy
    ,const tdouble3& vz,bool is2d=false);
  JBoxDef(const tdouble3& pmin,const tdouble3& pmax,bool is2d
    ,const tdouble3& mov,const tdouble3& rot,const tdouble3& rotcen);
  ~JBoxDef();
  JBoxDef& operator=(const JBoxDef& src);

  void Reset();

  void SetPos(const tdouble3& pmin,const tdouble3& pmax);
  void SetPos(const tdouble3& pmin,const tdouble3& vx
    ,const tdouble3& vy,const tdouble3& vz);
  void MovRotate(const tdouble3& mov,const tdouble3& rot
    ,const tdouble3& cen);
  void MovRotate(const tmatrix4d& mat);

  void Resize(double resize);

  bool     GetIs2D  ()const{ return(Is2D); }
  double   GetPosy2D()const{ return(Posy2D); }


  bool     IsSimple ()const{ return(Simple); }
  tdouble3 GetPosMin()const{ return(PosMin); }
  tdouble3 GetPosMax()const{ return(PosMax); }
  tdouble3 GetVx    ()const{ return(Vx); }
  tdouble3 GetVy    ()const{ return(Vy); }
  tdouble3 GetVz    ()const{ return(Vz); }

  double GetMinSize()const;

  std::string GetDomainStr()const;

  void GetBoxPoints(std::vector<tdouble3>& points)const;
  void GetBoxPoints(tdouble3* points)const;

  tdouble3 OutBoxMin()const;
  tdouble3 OutBoxMax()const;

  bool Inside(const JBoxDef& box)const;

  unsigned GetInsidePlanes(tplane3d* vplanes)const;

};

#endif


