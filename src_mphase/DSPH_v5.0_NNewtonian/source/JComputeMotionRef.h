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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase simple para calcular puntos de referencia para el calculo de 
//:#   movimiento de objetos a partir de particulas. (19-10-2020)
//:#############################################################################

/// \file JComputeMotionRef.h \brief Declares the class \ref JComputeMotionRef.

#ifndef _JComputeMotionRef_
#define _JComputeMotionRef_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"
#include "JMatrix4.h"
#include "JPartMotionDef.h"

//##############################################################################
//# JComputeMotionRef
//##############################################################################
/// \brief Computes reference data to computes motion starting from particles.

class JComputeMotionRef : protected JObject
{
private:
  word MkBoundFirst;      ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  std::vector<StMkMotionData> Mks; ///<Information of fixed and moving boundary Mk blocks.

  unsigned CaseNfixed;   ///<Number of fixed boundary particles. 
  unsigned CaseNmoving;  ///<Number of moving boundary particles. 
  unsigned CaseNfloat;   ///<Number of floating boundary particles. 
  unsigned CaseNbound;   ///<Number of boundary particles (fixed+moving+float). 

  unsigned *PartRidp;    ///<To look for particles by Id. [CaseNbound]

private:
  unsigned* ComputeRidp(unsigned np,const unsigned *idp);

public:
  JComputeMotionRef();
  ~JComputeMotionRef();
  void Reset();

  unsigned CountMkBlock()const{ return(unsigned(Mks.size())); }
  unsigned IdxMkBlock(word mkbound)const;
  const StMkMotionData& GetMkBlock(unsigned idx)const{ return(Mks[idx]); }

  void AddMkBlock(word mk,word mkbound,unsigned begin,unsigned np);
  void AddMkBlocks(word mkboundfirst,const std::vector<StMkMotionData> &mks);
  void AddMkBlocks(word mkboundfirst,unsigned nmk,const StMkMotionData *mks);

  void ComputeRefPoints(unsigned casenfixed,unsigned casenmoving
    ,unsigned casenfloat,unsigned np,const unsigned *idp,const tdouble3 *posd
    ,const tfloat3 *posf,const unsigned *ridp);

  void GetMotionRefData(StMkMotionData &v)const;
  void GetMotionRef(std::vector<StMkMotionData> &mks)const;
  void GetMotionRef(unsigned nmk,StMkMotionData *mks)const;

};

#endif


