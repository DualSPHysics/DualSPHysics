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

/// \file JSphFlexibleStructure.h \brief Declares the class \ref JSphFlexibleStructure.

#ifndef _JSphFlexibleStructureUndef_
#define _JSphFlexibleStructureUndef_

#ifdef DISABLE_FLEXSTRUCT

#include "JObject.h"
#include "DualSphDef.h"

class JXml;
class JSphMk;

//##############################################################################
//# JSphFlexibleStructureBody
//##############################################################################
/// \brief Manages the info of a single flexible structure.
class JSphFlexibleStructureBody : protected JObject
{
public:
  const unsigned IdBody;
  const word MkBound;

  JSphFlexibleStructureBody(unsigned idbody,word mkbound,float particlevolume
    ,float density,double youngmod,double poissratio,TpConstitModel constitmodel
    ,float hgfactor):IdBody(0),MkBound(0){};
  ~JSphFlexibleStructureBody(){};
  void ConfigBoundCode(typecode boundcode){};

  typecode GetBoundCode()const{ return(NULL); }

  float GetParticleVolume()const{ return(0); };
  float GetDensity()const{ return(0); };
  float GetYoungMod()const{ return(0); };
  float GetPoissRatio()const{ return(0); };
  float GetHgFactor()const{ return(0); };

  float GetParticleMass()const{ return(0); };
  tmatrix6f GetConstitMatrix()const{ return(TMatrix6f(0)); };

};

//##############################################################################
//# JSphFlexibleStructure
//##############################################################################
/// \brief Manages the info of flexible structures.
class JSphFlexibleStructure : protected JObject
{
public:

  JSphFlexibleStructure(bool simulate2d,double dp,JXml *sxml,const std::string &place,const JSphMk *mkinfo){};
  ~JSphFlexibleStructure() {};

  unsigned GetCount()const{ return(0); }
  const JSphFlexibleStructureBody* GetMkBody(unsigned idx)const{ return(NULL); }

  void ConfigCode(unsigned npb,typecode *code){};

};
#endif
#endif