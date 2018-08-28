//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Creacion de clase para gestionar informacion relativa a MK de particulas. (25-01-2018)
//:# - Nuevos metodos CountBlockType() y GetFirstBlockType(). (14-08-2018)
//:#############################################################################

/// \file JSphMk.h \brief Declares the class \ref JSphMk.

#ifndef _JSphMk_
#define _JSphMk_

#include <string>
#include <vector>
#include "JObject.h"
#include "Types.h"

class JSpaceParts;
class JPartDataHead;

//##############################################################################
//# JSphMkBlock
//##############################################################################
/// \brief Manages the info of each block of particles.
class JSphMkBlock : public JObject
{
private:
  bool PosDefined;
  tdouble3 PosMin;
  tdouble3 PosMax;

public:
  const bool Bound;       ///<Indicates whether a particle is boundary or not.
  const TpParticles Type; ///<Type of particle.
  const unsigned MkType;  ///<Label of block fluid or bound.
  const unsigned Mk;      ///<Absolute label.
  const typecode Code;
  const unsigned Begin;   ///<Id of the first particle of the block.
  const unsigned Count;   ///<Number of particles.

  JSphMkBlock(TpParticles type,unsigned mktype,unsigned mk,typecode code,unsigned begin,unsigned count);
  void Reset();

  bool GetPosDefined()const{ return(PosDefined); }
  tdouble3 GetPosMin()const{ return(PosMin); }
  tdouble3 GetPosMax()const{ return(PosMax); }
  void SetPosMinMax(const tdouble3 &pmin,const tdouble3 &pmax){ PosDefined=true; PosMin=pmin; PosMax=pmax; }
};


//##############################################################################
//# JSphMk
//##############################################################################
/// \brief Manages the info of particles according Mk.

class JSphMk : protected JObject
{
private:
  word MkBoundFirst;     ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  word MkFluidFirst;     ///<First Mk for fluid blocks (Mk=MkFluid+MkFluidFirst).

  std::vector<JSphMkBlock*> MkList;

  unsigned MkListSize;       ///<Total number of Mk blocks.
  unsigned MkListFixed;      ///<Number of Mk blocks of fixed type.
  unsigned MkListMoving;     ///<Number of Mk blocks of moving type.
  unsigned MkListFloat;      ///<Number of Mk blocks of floating type.
  unsigned MkListBound;      ///<Number of Mk blocks of boundary types. MkListBound=MkListFixed+MkListMoving+MkListFloat
  unsigned MkListFluid;      ///<Number of Mk blocks of fluid type.

  typecode CodeNewFluid;     ///<Code for new fluid particles created during the simulation.

public:
  JSphMk();
  ~JSphMk();
  void Reset();
  void Config(const JSpaceParts *parts);

  unsigned Size()const{ return(MkListSize); }
  const JSphMkBlock* Mkblock(unsigned c)const{ return(MkList[c]); }

  unsigned CountBlockType(TpParticles type)const;
  unsigned GetFirstBlockType(TpParticles type)const;

  typecode GetCodeNewFluid()const{ return(CodeNewFluid); }

  unsigned GetMkBlockById(unsigned id)const;
  typecode GetCodeById(unsigned id)const;

  word GetMkBoundFirst()const{ return(MkBoundFirst); }
  word GetMkFluidFirst()const{ return(MkFluidFirst); }

  unsigned GetMkBlockByMk(word mk)const;
  unsigned GetMkBlockByMkBound(word mkbound)const;
  unsigned GetMkBlockByMkFluid(word mkfluid)const;

  unsigned GetMkBlockByCode(typecode code)const;

  typecode CodeSetType(typecode code,TpParticles type,unsigned value)const;

  //void ComputeMkDomains(bool bound,const std::vector<unsigned> &mklist,unsigned np,const tdouble3 *pos,const typecode *code);
  void ComputeMkDomains(unsigned np,const tdouble3 *pos,const typecode *code);

  void ConfigPartDataHead(JPartDataHead *parthead)const;
};



#endif


