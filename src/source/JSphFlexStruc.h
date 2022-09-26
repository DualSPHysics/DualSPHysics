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

/// \file JSphFlexStruc.h \brief Declares the class \ref JSphFlexStruc.

#ifndef _JSphFlexStruc_
#define _JSphFlexStruc_

#include "JObject.h"
#include "DualSphDef.h"

#include <cmath>

class JLog2;
class JXml;
class TiXmlElement;
class JSphMk;

//##############################################################################
//# JSphFlexStrucBody
//##############################################################################
/// \brief Manages the info of a single flexible structure.
class JSphFlexStrucBody : protected JObject
{
private:
  JLog2* Log;

  //-Selection of particles
  typecode BoundCode;             ///<Code to select boundary particles.
  typecode ClampCode;             ///<Code to select clamping particles.

  //-Body parameters
  float ParticleVolume;           ///<Initial particle volume.
  float Density;                  ///<Initial particle density.
  float YoungMod;                 ///<Young's modulus.
  float PoissRatio;               ///<Poisson ratio.
  TpConstitModel ConstitModel;    ///<Constitutive model.
  float HgFactor;                 ///<Hourglass correction factor.

  tmatrix6f ConstitMatrix;        ///<Constitutive matrix.

  void Reset();

public:
  const unsigned IdBody;          ///<Flexible structure ID.
  const word MkBound;             ///<MkBound of flexible structure.
  const word MkClamp;             ///<MkBound of clamp.

  JSphFlexStrucBody(unsigned idbody,word mkbound,word mkclamp,float particlevolume,float density,double youngmod,double poissratio,TpConstitModel constitmodel,float hgfactor);
  ~JSphFlexStrucBody();

  void ConfigBoundCode(typecode boundcode);
  void ConfigClampCode(typecode clampcode);

  void GetConfig(std::vector<std::string> &lines)const;

  typecode GetBoundCode()const{ return(BoundCode); }
  typecode GetClampCode()const{ return(ClampCode); }

  float GetParticleVolume()const{ return ParticleVolume; };
  float GetDensity()const{ return Density; };
  float GetYoungMod()const{ return YoungMod; };
  float GetPoissRatio()const{ return PoissRatio; };
  float GetHgFactor()const{ return HgFactor; };

  float GetParticleMass()const{ return GetParticleVolume()*GetDensity(); };
  tmatrix6f GetConstitMatrix()const{ return ConstitMatrix; };
  float GetSoundSpeed()const{ return(float(sqrt(GetYoungMod()*(1.0-GetPoissRatio())/(GetDensity()*(1.0+GetPoissRatio())*(1.0-2.0*GetPoissRatio()))))); }

};

//##############################################################################
//# JSphFlexStruc
//##############################################################################
/// \brief Manages the info of flexible structures.
class JSphFlexStruc : protected JObject
{
private:
  JLog2* Log;
  const bool Simulate2D;
  const double Dp;                          ///<Initial distance between particles [m].

  std::vector<JSphFlexStrucBody*> List;     ///<List of flexible structure bodies.

  void Reset();
  bool ExistMk(word mkbound)const;
  void LoadXml(const JXml *sxml,const std::string &place);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void ConfigBoundCode(const JSphMk *mkinfo);
  void ConfigClampCode(const JSphMk *mkinfo);

public:

  JSphFlexStruc(bool simulate2d,double dp,JXml *sxml,const std::string &place,const JSphMk *mkinfo);
  ~JSphFlexStruc();

  void VisuConfig(std::string txhead,std::string txfoot);

  unsigned GetCount()const{ return(unsigned(List.size())); }
  const JSphFlexStrucBody* GetBody(unsigned idx)const{ return(idx<GetCount()?List[idx]:NULL); }

  void ConfigCode(unsigned npb,typecode *code);
  double GetInitialSoundSpeed();
};

#endif
