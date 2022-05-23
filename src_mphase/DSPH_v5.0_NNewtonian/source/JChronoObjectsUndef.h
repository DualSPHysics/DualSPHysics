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

/// \file JChronoObjects.h \brief Declares the empty class \ref JChronoObjects.

#ifndef _JChronoObjectsUndef_
#define _JChronoObjectsUndef_

//##############################################################################
//# JChronoObjects
//##############################################################################
/// \brief Manages Chrono objects.

#ifdef DISABLE_CHRONO
class JChronoObjects : protected JObject
{
protected:
public:
  JChronoObjects(const std::string &dirdata,const std::string &casename
    ,const JXml *sxml,const std::string &place,double dp,word mkboundfirst
    ,bool simulate2d){}
  ~JChronoObjects(){};
  void Reset(){};
  static bool Available(){ return(false); }

  bool UseDataDVI(word mkbound)const{ return(false); };
  bool GetUseCollision()const{ return(false); };
  unsigned GetCollisionShapes()const{ return(0); }

  bool ConfigBodyFloating(word mkbound,double mass,const tdouble3 &center
    ,const tmatrix3d &inertia,const tint3 &translationfree,const tint3 &rotationfree
    ,const tfloat3 &linvelini,const tfloat3 &angvelini){ return(false); };

  void ConfigDataBodyFloating(word mkbound,float kfric,float sfric,float restitu,float young,float poisson){};
  void ConfigDataBodyMoving  (word mkbound,float kfric,float sfric,float restitu,float young,float poisson){};
  void ConfigDataBodyFixed   (word mkbound,float kfric,float sfric,float restitu,float young,float poisson){};

  void Init(const JSphMk* mkinfo){};
  void VisuConfig(std::string txhead, std::string txfoot)const{};

  bool GetWithMotion()const{ return(false); }

  void SetFtData(word mkbound,const tfloat3 &face,const tfloat3 &fomegaace){};
  void SetFtDataVel(word mkbound,const tfloat3 &vlin,const tfloat3 &vang){};
  void GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const{};

  void SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double dt){};

  void RunChrono(unsigned nstep, double timestep, double dt, bool predictor){};

  void SavePart(int part){};
};
#endif

#endif


