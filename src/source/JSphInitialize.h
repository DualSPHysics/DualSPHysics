//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2019 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Implementa normales para cilindors (body & tank) y esfera (tank). (31-01-2020)
//:# - Limita calculo de normales con MaxDisteH. (31-01-2020)
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
  typedef enum{ 
    IT_FluidVel=1,
    IT_BoundNormalSet=2,       //<vs_mddbc>
    IT_BoundNormalPlane=3,     //<vs_mddbc>
    IT_BoundNormalSphere=4,    //<vs_mddbc>
    IT_BoundNormalCylinder=5,  //<vs_mddbc>
  }TpInitialize; 

public:
  const TpInitialize Type;   ///<Type of particle.

  JSphInitializeOp(TpInitialize type,const char* name):Type(type){ 
    ClassName=std::string("JSphInitializeOp_")+name;
  } 
  virtual ~JSphInitializeOp(){ DestructorActive=true; }
  virtual void ReadXml(JXml *sxml,TiXmlElement* ele)=0;
  virtual void Run(unsigned np,unsigned npb,unsigned nbound,const tdouble3 *pos
    ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)=0;
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
  JSphInitializeOp_FluidVel(JXml *sxml,TiXmlElement* ele):JSphInitializeOp(IT_FluidVel,"FluidVel"){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,unsigned nbound,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//<vs_mddbc_ini>
//##############################################################################
//# JSphInitializeOp_BoundNormalSet
//##############################################################################
/// Initializes normals of boundary particles.
class JSphInitializeOp_BoundNormalSet : public JSphInitializeOp
{
private:
  std::string MkBound;
  tfloat3 Normal;
public:
  JSphInitializeOp_BoundNormalSet(JXml *sxml,TiXmlElement* ele)
    :JSphInitializeOp(IT_BoundNormalSet,"BoundNormalSet"){ Reset(); ReadXml(sxml,ele); }
  void Reset(){ MkBound=""; Normal=TFloat3(0); }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,unsigned nbound,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JSphInitializeOp_BoundNormalPlane
//##############################################################################
/// Initializes normals of boundary particles.
class JSphInitializeOp_BoundNormalPlane : public JSphInitializeOp
{
private:
  const float H;
  std::string MkBound;
  tfloat3 Point;
  tfloat3 Normal;
  float MaxDisteH;  ///<Maximum distance to boundary limit. It uses H*distanceh (default=2).
public:
  JSphInitializeOp_BoundNormalPlane(JXml *sxml,TiXmlElement* ele,float h)
    :JSphInitializeOp(IT_BoundNormalPlane,"BoundNormalPlane"),H(h){ Reset(); ReadXml(sxml,ele); }
  void Reset(){ MkBound=""; Point=Normal=TFloat3(0); MaxDisteH=0; }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,unsigned nbound,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JSphInitializeOp_BoundNormalSphere
//##############################################################################
/// Initializes normals of boundary particles.
class JSphInitializeOp_BoundNormalSphere : public JSphInitializeOp
{
private:
  const float H;
  std::string MkBound;
  tfloat3 Center;
  float Radius;
  bool Inside;      ///<Boundary particles inside the sphere.
  float MaxDisteH;  ///<Maximum distance to boundary limit. It uses H*distanceh (default=2).
public:
  JSphInitializeOp_BoundNormalSphere(JXml *sxml,TiXmlElement* ele,float h)
    :JSphInitializeOp(IT_BoundNormalSphere,"BoundNormalSphere"),H(h){ Reset(); ReadXml(sxml,ele); }
  void Reset(){ MkBound=""; Center=TFloat3(0); MaxDisteH=Radius=0; Inside=true; }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,unsigned nbound,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JSphInitializeOp_BoundNormalCylinder
//##############################################################################
/// Initializes normals of boundary particles.
class JSphInitializeOp_BoundNormalCylinder : public JSphInitializeOp
{
private:
  const float H;
  std::string MkBound;
  tfloat3 Center1;
  tfloat3 Center2;
  float Radius;
  bool Inside;      ///<Boundary particles inside the cylinder.
  float MaxDisteH;  ///<Maximum distance to boundary limit. It uses H*distanceh (default=2).
public:
  JSphInitializeOp_BoundNormalCylinder(JXml *sxml,TiXmlElement* ele,float h)
    :JSphInitializeOp(IT_BoundNormalCylinder,"BoundNormalCylinder"),H(h){ Reset(); ReadXml(sxml,ele); }
  void Reset(){ MkBound=""; Center1=Center2=TFloat3(0); MaxDisteH=Radius=0; Inside=true; }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,unsigned nbound,const tdouble3 *pos,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  
//<vs_mddbc_end>


//##############################################################################
//# JSphInitialize
//##############################################################################
/// \brief Manages the info of particles from the input XML file.
class JSphInitialize  : protected JObject
{
private:
  const bool BoundNormals;
  const float H;
  std::vector<JSphInitializeOp*> Opes;

  void LoadFileXml(const std::string &file,const std::string &path);
  void LoadXml(JXml *sxml,const std::string &place);
  void ReadXml(JXml *sxml,TiXmlElement* lis);

public:
  JSphInitialize(const std::string &file,float h,bool boundnormals);
  ~JSphInitialize();
  void Reset();
  unsigned Count()const{ return(unsigned(Opes.size())); }

  void Run(unsigned np,unsigned npb,unsigned nbound,const tdouble3 *pos
    ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;

};

#endif


