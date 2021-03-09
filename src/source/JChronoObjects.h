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
//:# - Hace de interface con la libreria dsphchrono.dll. (03-05-2016).
//:# - Permite compilar sin librerias de CHRONO. (13-12-2019).
//:# - Permitidos SMooth Contacts (SMC) y ejecuciones Multi Core (17-01-2020).
//:# - Nuevo objeto link_coulomdamping (21-02-2020).
//:# - Nuevo objeto link_lulley (27-02-2020).
//:# - Se anade un bloque collision para activar o desactivar las colisiones en chrono (28-02-2020).
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (19-03-2020). 
//:#   y damping. (04-10-2020).
//:# - Contabiliza y muestra el numero de caras para las colisiones. (18-10-2020).
//:# - Permite imponer el valor de kfric de un objeto sobre otros. (30-10-2020).
//:# - Permite escalar fuerzas usando la etiqueta <scaleforce> (10-03-2021).
//:#############################################################################

/// \file JChronoObjects.h \brief Declares the class \ref JChronoObjects.

#ifndef _JChronoObjects_
#define _JChronoObjects_

#include "JObject.h"
#include "DualSphDef.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

//#define DISABLE_CHRONO     ///<It allows compile without ChronoLib library.

class JLog2;
class JSphMk;
class JXml;
class JVtkLib;
class TiXmlElement;
class JChronoData;
class JChValues;
class JChBody;
class JChBodyFloating;
class JChLink;
class DSPHChronoLib;
class JLinearValue;

//##############################################################################
//# JChronoObjects
//##############################################################################
/// \brief Manages Chrono objects.

#ifdef DISABLE_CHRONO
#include "JChronoObjectsUndef.h"
#else
class JChronoObjects : protected JObject
{
protected:
  JLog2 *Log;
  std::string DirData;
  std::string CaseName;
  const double Dp;
  const word MkBoundFirst;
  const bool Simulate2D;

  bool UseVariableCoeff; ///<Indicates the use of variable coefficients
  std::vector<JLinearValue*> StiffnessV; ///<For variable stiffness
  std::vector<JLinearValue*> DampingV;   ///<For variable damping

  int OmpThreads;     ///<Max number of OpenMP threads in execution on CPU host (minimum 1).
  const bool UseDVI;  ///<Uses Differential Variational Inequality (DVI) method.
  bool UseChronoSMC;  ///<Uses Smooth Contacts for collisions.
  bool UseCollision;  ///<Activates collisions between chrono objects.
  
  bool WithMotion;   ///<Some Chrono object with geometry is a moving object.

  void* Ptr_VtkSimple_AutoActual;
  void* Ptr_VtkSimple_AutoDp;

  JChronoData *ChronoDataXml; ///<Chrono data loaded from XML. | Datos de Chrono cargados del XML.
  DSPHChronoLib *ChronoLib;   ///<Object for the integrarion with Chrono Engine library. | Objeto para integracion con libreria de Chrono Engine.

  float CollisionDp;        ///<Allowed collision overlap according Dp (default=0.5).
  double SchemeScale;       ///<Scale value to create initial scheme of configuration.
  unsigned CollisionShapes; ///<Number of shapes for collisions created from VTK files.

  double SaveDataTime;  ///<Saves CSV with data exchange (0=all steps, <0:none).
  double NextTime;
  double LastTimeOk;

  void LoadPtrAutoActual(const JXml *sxml,std::string xmlrow);
  void LoadPtrAutoDp(const JXml *sxml,std::string xmlrow);
  
  unsigned CreateObjFiles(std::string idname,const std::vector<unsigned> &mkbounds
    ,std::string datadir,std::string mfile,byte normalmode,std::string diroutobj,std::string xmlrow);

  void LoadXml(const JXml *sxml, const std::string &place);
  void ReadXml(const JXml *sxml, TiXmlElement* lis);
  void ReadXmlValues(const JXml *sxml,TiXmlElement* lis,JChValues* values);
  std::string ReadXmlModelFile(const JXml *sxml,TiXmlElement* ele)const;
  void ConfigMovingBodies(const JSphMk* mkinfo);
  void ConfigOmp();

  void CheckParams(const JChBody *body)const;
  void VisuValues(const JChValues *values)const;
  void VisuBody(const JChBody *body)const;
  void VisuLink(const JChLink *link)const;
  void SaveVtkScheme_Spring(JVtkLib *sh,word mk,word mk1,tdouble3 pt0,tdouble3 pt1
    ,double restlength,double radius,double revlength,int nside)const;
  void SaveVtkScheme()const;

  void SetVariableCoeff(const double timestep);
  void ReadCoeffs(JChLink *link,const JXml *sxml,TiXmlElement* ele);

public:
  JChronoObjects(const std::string &dirdata,const std::string &casename
    ,const JXml *sxml,const std::string &place,double dp,word mkboundfirst
    ,bool simulate2d);
  ~JChronoObjects();
  void Reset();
  static bool Available(){ return(true); }

  bool UseDataDVI(word mkbound)const;
  bool GetUseCollision()const{ return(UseCollision); }
  unsigned GetCollisionShapes()const{ return(CollisionShapes); }

  bool ConfigBodyFloating(word mkbound,double mass,const tdouble3 &center
    ,const tmatrix3d &inertia,const tint3 &translationfree,const tint3 &rotationfree
    ,const tfloat3 &linvelini,const tfloat3 &angvelini);

  void ConfigDataBodyFloating(word mkbound,float kfric,float sfric,float restitu,float young,float poisson);
  void ConfigDataBodyMoving  (word mkbound,float kfric,float sfric,float restitu,float young,float poisson);
  void ConfigDataBodyFixed   (word mkbound,float kfric,float sfric,float restitu,float young,float poisson);

  void Init(const JSphMk* mkinfo);
  void VisuConfig(std::string txhead, std::string txfoot)const;

  bool GetWithMotion()const{ return(WithMotion); }

  void SetFtData(word mkbound,const tfloat3 &face,const tfloat3 &fomegaace);
  void SetFtDataVel(word mkbound,const tfloat3 &vlin,const tfloat3 &vang);
  void GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const;

  void SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double stepdt);

  void RunChrono(unsigned nstep,double timestep,double dt,bool predictor);

  void SavePart(int part);

  void ReadScaleForces(const JXml *sxml,TiXmlElement* lis);
};
#endif

#endif

