//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

//:NO_COMENTARIO
//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestiona la inizializacion inicial de las particulas (01-02-2017)
//:#############################################################################

/// \file JSphInitialize.h \brief Declares the class \ref JSphInitialize.

#ifndef _JSphInitialize_
#define _JSphInitialize_

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "JObject.h"
#include "TypesDef.h"

class JXml;
class TiXmlElement;

//##############################################################################
//# XML format in _FmtXML_Initialize.xml.
//##############################################################################

//##############################################################################
//# JSphInitializeOp
//##############################################################################
/// \brief Base clase for initialization of particle data.
class JSphInitializeOp : public JObject
{
public:
  ///<Types of initializations.
  typedef enum{ IT_FluidVel=1 }TpInitialize; 

public:
  const TpInitialize Type;   ///<Type of particle.

  JSphInitializeOp(TpInitialize type,const char* name):Type(type){ 
    ClassName=std::string("JSphInitializeOp_")+name;
  } 
  virtual ~JSphInitializeOp(){ DestructorActive=true; }
  virtual void ReadXml(JXml *sxml,TiXmlElement* ele)=0;
  virtual void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop)=0;
  virtual void GetConfig(std::vector<std::string> &lines)const=0;
};

//##############################################################################
//# JSphInitializeOp_FluidVel
//##############################################################################
/// Initializes velocity of fluid particles.
class JSphInitializeOp_FluidVel : public JSphInitializeOp
{
private:
  ///Controls profile of imposed velocity.
  typedef enum{ 
    TVEL_Constant=0,    ///<Velocity profile uniform.
    TVEL_Linear=1,      ///<Velocity profile linear.
    TVEL_Parabolic=2    ///<Velocity profile parabolic.
  }TpVelocity;

private:
  TpVelocity VelType;  ///<Type of velocity.
  std::string MkFluid;
  tfloat3 Direction;
  float Vel1,Vel2,Vel3;
  float Posz1,Posz2,Posz3;

public:
  JSphInitializeOp_FluidVel(JXml *sxml,TiXmlElement* ele):JSphInitializeOp(IT_FluidVel,"FluidVel"){ 
    Reset();
    ReadXml(sxml,ele); 
  }
  void Reset();
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop);
  void GetConfig(std::vector<std::string> &lines)const;
};  


//##############################################################################
//# JSphInitialize
//##############################################################################
/// \brief Manages the info of particles from the input XML file.
class JSphInitialize  : protected JObject
{
private:
  std::vector<JSphInitializeOp*> Opes;

  void LoadFileXml(const std::string &file,const std::string &path);
  void LoadXml(JXml *sxml,const std::string &place);
  void ReadXml(JXml *sxml,TiXmlElement* lis);

public:
  JSphInitialize(const std::string &file);
  ~JSphInitialize();
  void Reset();
  unsigned Count()const{ return(unsigned(Opes.size())); }

  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop);
  void GetConfig(std::vector<std::string> &lines)const;

};

#endif


