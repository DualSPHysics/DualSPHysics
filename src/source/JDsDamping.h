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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para gestionar zonas de amortiguamiento. (13-01-2014)
//:# - Modificacion para usar la funcion exponencial de S.J. Lind et al. 2012. (02-09-2014)
//:# - Muestra vector normal del plano para facilitar su comprension. (26-11-2014)
//:# - Se puede aplicar un factor de amortiguacion para cada componente. (26-11-2014)
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:# - Mejora configuracion XML para definir el dominio del damping. (21-12-2018)
//:# - Genera VTK con esquema de configuraciones. (21-12-2018)
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (18-03-2020)  
//:# - Cambio de nombre de J.Damping a J.DsDamping. (28-06-2020)
//:# - Mejora de codigo para la implementacion de nuevas opcioens de damping. (14-10-2021)
//:#############################################################################

/// \file JDsDamping.h \brief Declares the class \ref JDsDamping.

#ifndef _JDsDamping_
#define _JDsDamping_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

class JXml;
class TiXmlElement;
class JLog2;
class JVtkLib;

//##############################################################################
//# XML format in _FmtXML_Damping.xml.
//##############################################################################

//##############################################################################
//# JDsDampingOp
//##############################################################################
/// \brief Base clase for damping configurations.
class JDsDampingOp : public JObject
{
public:
  ///Types of damping configurations.
  typedef enum{ 
    DA_Plane=1
   ,DA_Box=2
   ,DA_Cylinder=3
  }TpDamping; 

  ///Returns damping type as string.
  static std::string GetNameType(TpDamping type){
    switch(type){
      case DA_Plane:    return("Plane");
      case DA_Box:      return("Box");
      case DA_Cylinder: return("Cylinder");
    }
    return("???");
  }

protected:
  float OverLimit;    ///<Distance after limit with maximum reduction. | Distancia despues de limite con reduccion maxima.
  float ReduMax;      ///<Percentage of maximum reduction. | Porcentaje de reduccion maxima.
  tfloat3 Factorxyz;  ///<Factor applied on each axis. | Factor de aplicacion en cada eje.
public:
  const unsigned Id;    ///<Id of damping.
  const TpDamping Type; ///<Type of damping.

public:
  JDsDampingOp(unsigned id,TpDamping type):Id(id),Type(type)
  { 
    ClassName=std::string("JDsDampingOp_")+GetNameType(type);
    Reset();
  } 
  virtual ~JDsDampingOp(){ DestructorActive=true; }
  void Reset(){ OverLimit=ReduMax=0; Factorxyz=TFloat3(0); }

  float   GetOverLimit()const{ return(OverLimit); };
  float   GetRedumax  ()const{ return(ReduMax); };
  tfloat3 GetFactorxyz()const{ return(Factorxyz); };

  virtual void ReadXml(const JXml *sxml,TiXmlElement* ele)=0;
  virtual void SaveVtkConfig(double dp,JVtkLib *sh)const=0;
  virtual void GetConfig(std::vector<std::string> &lines)const=0;

  virtual void ComputeDampingCpu(double dt,unsigned n,unsigned pini
    ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const=0;

#ifdef _WITHGPU
  virtual void ComputeDampingGpu(double dt,unsigned n,unsigned pini
    ,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop)const=0;
#endif
};


//##############################################################################
//# JDsDampingOp_Plane
//##############################################################################
/// Damping according to distance to plane with or without domain limit.
class JDsDampingOp_Plane : public JDsDampingOp
{
private:
  tdouble3 LimitMin;  ///<Minimal reduction position. | Posicion de reduccion minima.
  tdouble3 LimitMax;  ///<Miximum reduction position. | Posicion de reduccion maxima.

  bool UseDomain;     ///<Indicates use of domain planes. | Indica uso de planos del dominio.
  double DomzMin;     ///<Domain definition - Z minimum. | Definicion de dominio - Z minima.
  double DomzMax;     ///<Domain definition - Z maximum. | Definicion de dominio - Z maxima.
  tdouble2 DomPt0;    ///<Domain point.
  tdouble2 DomPt1;    ///<Domain point.
  tdouble2 DomPt2;    ///<Domain point.
  tdouble2 DomPt3;    ///<Domain point.
  tplane3d DomPla0;   ///<Domain definition - plane 0.
  tplane3d DomPla1;   ///<Domain definition - plane 1.
  tplane3d DomPla2;   ///<Domain definition - plane 2.
  tplane3d DomPla3;   ///<Domain definition - plane 3.

  tplane3d Plane;     ///<Plane at the limitmin point. | Plano en el punto limitmin.
  float Dist;         ///<Distance between limitmin and limitmax points. | Distancia entre puntos limitmin y limitmax.

private:
  void ComputeDomPlanes();

public:
  JDsDampingOp_Plane(unsigned id,const JXml *sxml,TiXmlElement* ele)
    :JDsDampingOp(id,DA_Plane){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void SaveVtkConfig(double dp,JVtkLib *sh)const;
  void GetConfig(std::vector<std::string> &lines)const;

  void ComputeDampingCpu(double dt,unsigned n,unsigned pini
    ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const;

#ifdef _WITHGPU
  void ComputeDampingGpu(double dt,unsigned n,unsigned pini
    ,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop)const;
#endif
};  


//##############################################################################
//# JDsDampingOp_Box
//##############################################################################
/// Damping according to distance to box limits.
class JDsDampingOp_Box : public JDsDampingOp
{
public:
///Direction mode.
  typedef enum{ 
    BDIR_Void=0
   ,BDIR_Top=1
   ,BDIR_Bottom=2
   ,BDIR_Left=4
   ,BDIR_Right=8
   ,BDIR_Front=16
   ,BDIR_Back=32
   ,BDIR_All=63
   ,BDIR_Error=64
  }TpDirections;  

  static std::string GetBoxDirText(TpDirections bdir);
  static TpDirections GetBoxDir(std::string txdir);


private:
  TpDirections Directions; ///<Define the activated directions.
  tdouble3 LimitMin1;  ///<Initial box point for minimal reduction.
  tdouble3 LimitMin2;  ///<Final box point for minimal reduction.
  tdouble3 LimitMax1;  ///<Initial box point for Miximum reduction.
  tdouble3 LimitMax2;  ///<Final box point for Miximum reduction.

  tdouble3 LimitOver1; ///<Initial box point for overlimit domain.
  tdouble3 LimitOver2; ///<Final box point for overlimit domain.

  tdouble3 BoxSize1;
  tdouble3 BoxSize2;

  //tplane3d Plane;     ///<Plane at the limitmin point. | Plano en el punto limitmin.
  //float Dist;         ///<Distance between limitmin and limitmax points. | Distancia entre puntos limitmin y limitmax.

private:
  TpDirections ReadXmlDirections(const JXml *sxml,TiXmlElement* ele)const;

public:
  JDsDampingOp_Box(unsigned id,const JXml *sxml,TiXmlElement* ele)
    :JDsDampingOp(id,DA_Box){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void SaveVtkConfig(double dp,JVtkLib *sh)const;
  void GetConfig(std::vector<std::string> &lines)const;

  void ComputeDampingCpu(double dt,unsigned n,unsigned pini
    ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const;

#ifdef _WITHGPU
  void ComputeDampingGpu(double dt,unsigned n,unsigned pini
    ,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop)const;
#endif
};  


//##############################################################################
//# JDsDampingOp_Cylinder
//##############################################################################
/// Damping according to distance to cylinder limits.
class JDsDampingOp_Cylinder : public JDsDampingOp
{
private:
  tdouble3 Point1;   ///<Point for axis definition.
  tdouble3 Point2;   ///<Point for axis definition.
  double LimitMin;  ///<Radius for minimal reduction.
  double LimitMax;  ///<Radius for maximum reduction.

  //tplane3d Plane;     ///<Plane at the limitmin point. | Plano en el punto limitmin.
  //float Dist;         ///<Distance between limitmin and limitmax points. | Distancia entre puntos limitmin y limitmax.

public:
  JDsDampingOp_Cylinder(unsigned id,const JXml *sxml,TiXmlElement* ele)
    :JDsDampingOp(id,DA_Cylinder){ Reset(); ReadXml(sxml,ele); }
  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void SaveVtkConfig(double dp,JVtkLib *sh)const;
  void GetConfig(std::vector<std::string> &lines)const;

  void ComputeDampingCpu(double dt,unsigned n,unsigned pini
    ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const;

#ifdef _WITHGPU
  void ComputeDampingGpu(double dt,unsigned n,unsigned pini
    ,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop)const;
#endif
};  


//##############################################################################
//# JDsDamping
//##############################################################################
/// \brief Manages the info of damping zones.

class JDsDamping : protected JObject
{
private:
  JLog2* Log;
  const double Dp;      ///<Initial distance between particles [m].
  std::vector<JDsDampingOp*> List;

  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void SaveVtkConfig(double dp)const;

public:
  JDsDamping(double dp);
  ~JDsDamping();
  void Reset();

  unsigned Count()const{ return(unsigned(List.size())); }

  void LoadXml(const JXml *sxml,const std::string &place);
  void VisuConfig(std::string txhead,std::string txfoot);

  void ComputeDampingCpu(double timestep,double dt,unsigned n,unsigned pini
    ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const;

#ifdef _WITHGPU
  void ComputeDampingGpu(double timestep,double dt,unsigned n,unsigned pini
    ,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop)const;
#endif
};


#endif


