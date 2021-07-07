//HEAD_DSPH
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

//:NO_COMENTARIO
//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestiona la inizializacion inicial de las particulas (01-02-2017)
//:# - Implementa normales para cilindors (body & tank) y esfera (tank). (31-01-2020)
//:# - Limita calculo de normales con MaxDisteH. (31-01-2020)
//:# - Objeto JXml pasado como const para operaciones de lectura. (18-03-2020)  
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (19-03-2020)  
//:# - Opcion para calcular boundary limit de forma automatica. (19-05-2020)  
//:# - Cambio de nombre de J.SphInitialize a J.DsInitialize. (28-06-2020)
//:# - Error corregido al obtener nombre de operacion a partir de la clase. (02-07-2020)
//:# - New filter onlypos according to particle position. (25-07-2020)
//:# - Nueva opcion IT_BoundNormalParts para calcular normales a partir de 
//:#   particulas finales (en desarrollo y solo para 2D). (07-07-2020)
//:# - Se permite la configuracion por parametros de ejecucion de IT_BoundNormalPlane 
//:#   y IT_BoundNormalParts. (07-07-2020)
//:#############################################################################

/// \file JDsInitialize.h \brief Declares the class \ref JDsInitialize.

#ifndef _JDsInitialize_
#define _JDsInitialize_

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
//# JDsInitializeOp
//##############################################################################
/// \brief Base clase for initialization of particle data.
class JDsInitializeOp : public JObject
{
public:
  ///Types of initializations.
  typedef enum{ 
    IT_FluidVel=10
   ,IT_BoundNormalSet=30
   ,IT_BoundNormalPlane=31
   ,IT_BoundNormalSphere=32
   ,IT_BoundNormalCylinder=33
   ,IT_BoundNormalParts=34
  }TpInitialize; 

  ///Structure with constant values needed for initialization tasks.
  typedef struct StrInitCt{
    bool simulate2d;         ///<Toggles 2D simulation.
    double simulate2dposy;   ///<Y value in 2D simulations.
    tdouble3 maprealposmin;  ///<Real lower limit simulation.
    tdouble3 maprealposmax;  ///<Real upper limit simulation.
    double dp;               ///<Initial distance between particles [m].
    float kernelh;           ///<The smoothing length of SPH kernel [m].
    unsigned nbound;         ///<Initial number of boundary particles (fixed+moving+floating).
    std::string dirdatafile; ///<Directory to data files.

    StrInitCt(bool sim2d,double sim2dy,tdouble3 posmin,tdouble3 posmax
      ,double dp_,float kernelh_,unsigned nbound_,std::string dirdatafile_)
    {
      simulate2d=sim2d; simulate2dposy=sim2dy;
      maprealposmin=posmin; maprealposmax=posmax; dp=dp_;
      kernelh=kernelh_; nbound=nbound_; dirdatafile=dirdatafile_;
    }
  }StInitCt;

protected:
  bool OnlyPos;         ///<Activate filter according to position.
  tdouble3 OnlyPosMin;  ///<Minimum positon for filtering.
  tdouble3 OnlyPosMax;  ///<Maximum positon for filtering.
  unsigned NpUpdated;   ///<Number of updated particles.
  unsigned NpTotal;     ///<Total number of particles.
public:
  const TpInitialize Type;   ///<Type of particle.
  const StInitCt InitCt;     ///<Constant values needed for initialization tasks.
  const unsigned BaseNameSize;

public:
  JDsInitializeOp(TpInitialize type,const char* name,StInitCt initct)
    :Type(type),InitCt(initct),BaseNameSize(unsigned(std::string("JDsInitializeOp_").size()))
  { 
    ClassName=std::string("JDsInitializeOp_")+name;
    Reset();
  } 
  virtual ~JDsInitializeOp(){ DestructorActive=true; }
  void Reset();
  void ReadXmlOnlyPos(const JXml *sxml,TiXmlElement* ele);
  virtual void ReadXml(const JXml *sxml,TiXmlElement* ele)=0;
  virtual void Run(unsigned np,unsigned npb,const tdouble3 *pos
    ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal)=0;
  virtual void GetConfig(std::vector<std::string> &lines)const=0;
  unsigned ComputeDomainMk(bool bound,word mktp,unsigned np,const word *mktype
    ,const unsigned *idp,const tdouble3 *pos,tdouble3 &posmin,tdouble3 &posmax)const;
  inline bool CheckPos(unsigned p,const tdouble3 *pos){
    NpTotal++;
    const bool sel=(!OnlyPos || (OnlyPosMin<=pos[p] && pos[p]<=OnlyPosMax));
    if(sel)NpUpdated++;
    return(sel);
  }
  std::string GetConfigNp()const;
  std::string GetConfigMkBound(std::string mktype)const;
  std::string GetConfigMkFluid(std::string mktype)const;
  std::string GetConfigOnlyPos()const;
};

//##############################################################################
//# JDsInitializeOp_FluidVel
//##############################################################################
/// Initializes velocity of fluid particles.
class JDsInitializeOp_FluidVel : public JDsInitializeOp
{
private:
  ///Controls profile of imposed velocity.
  typedef enum{ 
    TVEL_Constant=0   ///<Velocity profile uniform.
   ,TVEL_Linear=1     ///<Velocity profile linear.
   ,TVEL_Parabolic=2  ///<Velocity profile parabolic.
  }TpVelocity;
private:
  TpVelocity VelType;  ///<Type of velocity.
  std::string MkFluid;
  tfloat3 Direction;
  float Vel1,Vel2,Vel3;
  float Posz1,Posz2,Posz3;
public:
  JDsInitializeOp_FluidVel(const JXml *sxml,TiXmlElement* ele,StInitCt initct)
    :JDsInitializeOp(IT_FluidVel,"FluidVel",initct){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp
    ,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JDsInitializeOp_BoundNormalSet
//##############################################################################
/// Initializes normals of boundary particles.
class JDsInitializeOp_BoundNormalSet : public JDsInitializeOp
{
private:
  std::string MkBound;
  tfloat3 Normal;
public:
  JDsInitializeOp_BoundNormalSet(const JXml *sxml,TiXmlElement* ele,StInitCt initct)
    :JDsInitializeOp(IT_BoundNormalSet,"BoundNormalSet",initct){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp
    ,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JDsInitializeOp_BoundNormalPlane
//##############################################################################
/// Initializes normals of boundary particles.
class JDsInitializeOp_BoundNormalPlane : public JDsInitializeOp
{
private:
  std::string MkBound;
  bool PointAuto;   ///<Point is calculated automatically accoding to normal configuration.
  float LimitDist;  ///<Minimun distance (Dp*vdp) between particles and boundary limit to calculate the point (default=0.5).
  tfloat3 Point;
  tfloat3 Normal;
  float MaxDisteH;  ///<Maximum distance to boundary limit. It uses H*distanceh (default=2).
public:
  JDsInitializeOp_BoundNormalPlane(const JXml *sxml,TiXmlElement* ele,StInitCt initct)
    :JDsInitializeOp(IT_BoundNormalPlane,"BoundNormalPlane",initct){ Reset(); ReadXml(sxml,ele); }
  JDsInitializeOp_BoundNormalPlane(const std::string &eparm,StInitCt initct)
    :JDsInitializeOp(IT_BoundNormalPlane,"BoundNormalPlane",initct){ Reset(); ReadKeyvals(eparm); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void ReadKeyvals(const std::string &eparm);
  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp
    ,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JDsInitializeOp_BoundNormalSphere
//##############################################################################
/// Initializes normals of boundary particles.
class JDsInitializeOp_BoundNormalSphere : public JDsInitializeOp
{
private:
  std::string MkBound;
  tfloat3 Center;
  float Radius;
  bool Inside;      ///<Boundary particles inside the sphere.
  float MaxDisteH;  ///<Maximum distance to boundary limit. It uses H*distanceh (default=2).
public:
  JDsInitializeOp_BoundNormalSphere(const JXml *sxml,TiXmlElement* ele,StInitCt initct)
    :JDsInitializeOp(IT_BoundNormalSphere,"BoundNormalSphere",initct){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp
    ,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JDsInitializeOp_BoundNormalCylinder
//##############################################################################
/// Initializes normals of boundary particles.
class JDsInitializeOp_BoundNormalCylinder : public JDsInitializeOp
{
private:
  std::string MkBound;
  tfloat3 Center1;
  tfloat3 Center2;
  float Radius;
  bool Inside;      ///<Boundary particles inside the cylinder.
  float MaxDisteH;  ///<Maximum distance to boundary limit. It uses H*distanceh (default=2).
public:
  JDsInitializeOp_BoundNormalCylinder(const JXml *sxml,TiXmlElement* ele,StInitCt initct)
    :JDsInitializeOp(IT_BoundNormalCylinder,"BoundNormalCylinder",initct){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp
    ,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  

//##############################################################################
//# JDsInitializeOp_BoundNormalParts
//##############################################################################
/// Initializes normals of boundary particles.
class JDsInitializeOp_BoundNormalParts : public JDsInitializeOp
{
private:
  std::string MkBound;
  float MaxDisteH;  ///<Maximum distance to boundary limit. It uses H*distanceh (default=2).
public:
  JDsInitializeOp_BoundNormalParts(const JXml *sxml,TiXmlElement* ele,StInitCt initct)
    :JDsInitializeOp(IT_BoundNormalParts,"BoundNormalParts",initct){ Reset(); ReadXml(sxml,ele); }
  JDsInitializeOp_BoundNormalParts(const std::string &eparm,StInitCt initct)
    :JDsInitializeOp(IT_BoundNormalPlane,"BoundNormalParts",initct){ Reset(); ReadKeyvals(eparm); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void ReadKeyvals(const std::string &eparm);
  void Run(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp
    ,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;
};  


//##############################################################################
//# JDsInitialize
//##############################################################################
/// \brief Manages the info of particles from the input XML file.
class JDsInitialize  : protected JObject
{
private:
  const bool BoundNormals;
  const JDsInitializeOp::StInitCt InitCt;  ///<Constant values needed for initialization tasks.
  std::vector<JDsInitializeOp*> Opes;

  void LoadFileXml(const std::string &file,const std::string &path);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);

public:
  JDsInitialize(bool sim2d,double sim2dy,tdouble3 posmin,tdouble3 posmax
    ,double dp,float kernelh,const std::string &dirdatafile
    ,unsigned nbound,bool boundnormals);
  ~JDsInitialize();
  void Reset();

  void LoadXml(const JXml *sxml,const std::string &place);
  void LoadExecParms(const std::vector<std::string> &execparms);

  unsigned Count()const{ return(unsigned(Opes.size())); }

  void Run(unsigned np,unsigned npb,const tdouble3 *pos
    ,const unsigned *idp,const word *mktype,tfloat4 *velrhop,tfloat3 *boundnormal);
  void GetConfig(std::vector<std::string> &lines)const;

};

#endif


