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
//:# - Clase para almacenar informacion inicial de particulas usada en configuraciones automaticas. (09-10-2018)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Cambio de nombre de J.SphPartsInit a J.DsPartsInit. (28-06-2020)
//:#############################################################################

/// \file JDsPartsInit.h \brief Declares the class \ref JDsPartsInit.

#ifndef _JDsPartsInit_
#define _JDsPartsInit_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"

class JSphMk;

//##############################################################################
//# JDsPartsInit
//##############################################################################
/// \brief Stores initial particle information for automatic configurations.
class JDsPartsInit : protected JObject
{
private:
  const JSphMk* MkInfo;

  unsigned Np;
  tdouble3 *Pos;
  typecode *Code;

  void FreeParticleData();
  void LoadParticleData(unsigned np,const tdouble3 *pos,const typecode *code);

public:
  const bool Simulate2D;         ///<Indicates 2D simulation.
  const double Simulate2DPosY;   ///<Y value in 2D simulations.
  const double Dp;               ///<Distance between particles.

  JDsPartsInit(bool simulate2d,double simulate2dposy,double dp
    ,const JSphMk* mkinfo,unsigned np,const tdouble3 *pos,const typecode *code);
  ~JDsPartsInit();
  void Reset();

  unsigned GetNp()const{ return(Np); }
  const typecode* GetCode()const{ return(Code); }
  const tdouble3* GetPos()const{ return(Pos); }
  const JSphMk* GetMkInfo()const{ return(MkInfo); }
};


#endif


