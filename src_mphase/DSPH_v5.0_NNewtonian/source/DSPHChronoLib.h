/*
 <DUALSPHYSICS>  Copyright (c) 2019, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

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
//:# - Recibe parametros del objeto en una structuctura. (v0.01 / 03-05-2016)
//:# - Codigo adaptado al estilo de DualSPHysics. (v0.02 / 03-05-2016)
//:# - Uso de JChronoData para configurar y gestionar varios objetos con 
//:#   con multiples uniones. (v0.03 / 10-05-2016)
//:# - Nuevo atributo DataDir en JChronoData. (v0.04 / 10-05-2016)
//:# - Cambio de solver e integrador para evitar valores NaN. (v0.05 / 23-05-2016)
//:# - Procesa todos los objetos de una vez. (v0.06 / 26-05-2016)
//:# - Se cambio el nombre de fomegavel por fomegaace. (v0.07 / 27-05-2016)
//:# - Acoplamiento con Chrono4.0.0. (v1.01 / 13-12-2019)
//:# - Permite ejecuciones Multi Core utilizando el modulo MODULE_PARALLEL. (v2.01 / 20-12-2019)
//:# - Metodo de contacto SMC (SMooth Contact) (v3.01 / 14-01-2020)
//:# - Envio desde DSPH los parametros Young's Modulos y Poisson Ratio para SMC. (v3.02 /16-01-2020)
//:# - Se separan en dos clases los modos de ejecucion Single Core y Multi Core. (v3.03 / 07-02-2020)
//:# - Nuevo objeto link_pulley. (v3.04 / 24-02-2020).
//:# - Nuevo objeto link_coulombdamping. (v3.05 / 27-02-2020)
//:# - Permitidas las ejecuciones Multi Core con SMC en Linux. (v3.06 / 17-04-2020)
//:# - Funciones ApplyInitialVel y ApplyImposedVel para aplicar velocidades
//:#   externas. (v3.07 / 29-06-2020)
//:# - Simulaciones con objetos FEA. (v4.01 / 14-07-2020)
//:# - Elementos flexibles con formas circulares 2D (v4.02 / 14-07-2020)
//:# - Elementos flexibles con formas rectangulares 2D (v4.03 / 28-09-2020)
//:# - Permite la simulacion de ChLinks con coeficientes variables de stiffness 
//:#   y damping. (v4.04 / 04-10-2020)
//:# - Activadas las colisiones entre objetos flexibles y objetos rigidos usando
//:#   la clase chrono/fea/ChContactSurfaceMesh (v4.05 / 08-10-2020)
//:# - Imponer valor de friction de un material sobre otros (v4.06 / 30-10-2020) 
//:# - Permite seleccionar el modo de colision de los objetos flexibles (v4.07 / 06-11-2020) 
//:# - Nuevo link ChLinkPointRotFrame para conectar FEA nodos a objetos rigidos (v4.08 / 09-11-2020)
//:# - Posibilidad de multiplicar las fuerzas por coeficientes introducidos por
//:#   por el usuario para escalar las fuerzas en simulaciones 2D (v4.09 / 15-12-2020)
//:# - Uso de coeficientes de friccion dinamico (Kfric) y estatico (Sfric) (v4.10 / 01-03-2021)
//:# - Se permite escalar las fuerzas de manera individual para cada uno de los
//:#   objetos de chrono (v4.11 / 10-03-2021)
//:#############################################################################

/// \file DSPHChronoLib.h \brief Declares the class \ref DSPHChronoLib which is the interface between DualSPHysics and Chrono.

#ifndef DSPHCHRONOLIB_H
#define DSPHCHRONOLIB_H

#include "TypesDef.h"
#include "JChronoData.h"
#include <string>
#include <iostream>
#include <memory>

//-Forward declarations to avoid including chrono classes.
namespace chrono {
  class ChSystem;
  class ChSystemParallel;
  class ChMaterialSurface;
  class ChBody;
  //<vs_chronoo_fea_ini>
  namespace fea{
    class ChMesh;
  }//<vs_chronoo_fea_end>
};

class BuildersFEA;	//<vs_chronoo_fea>

//##############################################################################
//# DSPHChronoLib
//##############################################################################
/// \brief Defines the class interface between DualSPHysics and Chrono.
class DSPHChronoLib {
public:
  //-States of execution.
  typedef enum { RSTATE_Init,RSTATE_Loading,RSTATE_Results }TpRunState;
  const std::string version;       ///<DualSPHysics version
  const std::string DsphChVersion; ///<Interface version

protected:
  const std::string ClassName;

  //-Chrono physical system.
  std::string DirOut;
  bool Simulate2D;          ///<True for 2D Simulations.
  BuildersFEA *BuildList;   ///<Array of BuilderFEA.
  bool UseOmp;              ///<Indicates if use of ChronoEngine_Parallel module is enabled.
  int OmpThreads;           ///<Threads number used by OpenMP.
  unsigned SolverIndex;     ///<Indicates the index of chrono solver enum.
  unsigned MaxIter;         ///<Indicates the maximun number of iterations for the solver.
  bool DG;                  ///<Used for Debug.
  bool UseSMC;              ///<True if it is using SMC (SMooth Contacts) 
  bool UseFEA;              ///<True if the use of Finite Element Analysis is enabled.
  JChronoData ChData;
  TpRunState RunState;
  double CollisionCoef;

  /// Constructor
  DSPHChronoLib(const JChronoData &chdata);

  /// Initialisation of variables
  void Reset();

  /// Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  virtual void SaveForcesHead(){};
  
  /// Establishes the variable coefficients to the link objects.
  virtual void SetVariableCoeff(){};
  
  /// Returns a valid unit vector of the vector or (0,0,0).
  static tdouble3 VecUnitarySafe(const tdouble3 &v);

  /// Adds the material properties to a object to enable collisions
  void ConfigSurfaceBody(const JChBody &body,chrono::ChBody *chbody);
   
  /// Adds the initial velocity.
  void ApplyInitialVel(const JChBody &body,chrono::ChBody *chbody);

  /// Adds the imposed velocity.
  void ApplyImposedVel(const JChBodyFloating &body,chrono::ChBody *chbody);

public:

  /// Loads data for bodies and configures objects.
  virtual void Config(std::string dirout,bool svdata,bool simulate2d){};

  /// Loads inertia for bodies.
  virtual void Config_Inertia(){};

  /// Compute a single timestep for each floating and moving body.
  virtual bool RunChrono(double timestep,double dt,bool predictor)=0;
 
  /// Saves forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  virtual void SaveForces(){};

  /// Obtains positions of Spring link.
  virtual bool GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const=0;

  /// Obtains RestLength of Spring link.
  virtual double GetSpringLinkRestLength(const std::string &linkname)const=0;

  /// Modifies RestLength of Spring link.
  virtual void SetSpringLinkRestLength(const std::string &linkname,double restlength)const{};

  /// Obtains center of body.
  virtual bool GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const=0;

  /// Returns pointer to ChronoData object.
  const JChronoData* GetChronoData(){ return(&ChData); }

  /// Loads floating data to calculate coupling with Chrono.
  bool SetFtData(word mkbound,const tfloat3 &face,const tfloat3 &fomegaace);

  /// Loads imposed velocity for floating to calculate coupling with Chrono.
  bool SetFtDataVel(word mkbound,const tfloat3 &vlin,const tfloat3 &vang); 
  
  /// Obtains floating data from coupling with Chrono.
  bool GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const;

  /// Loads motion data to calculate coupling with Chrono.
  bool SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double stepdt);
};

//##############################################################################
//# DSPHChronoLibSC
//##############################################################################
/// \brief Defines the class for single-core executions.
class DSPHChronoLibSC : public DSPHChronoLib {
private:
  chrono::ChSystem *MphysicalSystem; ///<Pointer to Chrono System

  /// Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  void SaveForcesHead();
  
  /// Establishes the variable coefficients to the link objects.
  void SetVariableCoeff();
  
  /// Configures floating bodies
  void ConfigFloating(const JChBody* body);

  /// Configures moving bodies
  void ConfigMoving(const JChBody* body);

  /// Configures fixed bodies
  void ConfigFixed(const JChBody* body);

  
public:
  /// Constructor
  DSPHChronoLibSC(const JChronoData &chdata);

  /// Destructor
  ~DSPHChronoLibSC();

  /// Loads data for bodies and configures objects.
  void Config(std::string dirout,bool svdata,bool simulate2d);

  /// Loads inertia for bodies.
  void Config_Inertia();

  /// Compute a single timestep for each floating and moving body.
  bool RunChrono(double timestep,double dt,bool predictor);

  /// Saves forces for each body and link (ChronoLink_forces.csv,ChronoBody_forces.csv).
  void SaveForces();

  /// Obtains positions of Spring link.
  bool GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const;

  /// Obtains RestLength of Spring link.
  double GetSpringLinkRestLength(const std::string &linkname)const;

  /// Modifies RestLength of Spring link.
  void SetSpringLinkRestLength(const std::string &linkname,double restlength)const;

  /// Obtains center of body.
  bool GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const;
};


#ifndef DISABLE_CHRONO_OMP
//##############################################################################
//# DSPHChronoLibMC
//##############################################################################
/// \brief Defines the class for multi-core executions.
class DSPHChronoLibMC : public DSPHChronoLib {
private:
  chrono::ChSystemParallel *MphysicalSystem;  ///<Pointer to Chrono System 

  /// Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  void SaveForcesHead();

  /// Establishes the variable coefficients to the link objects.
  void SetVariableCoeff();

  /// Configures floating bodies
  void ConfigFloating(const JChBody* body);

  /// Configures moving bodies
  void ConfigMoving(const JChBody* body);

  /// Configures fixed bodies
  void ConfigFixed(const JChBody* body);

public:
  /// Constructor
  DSPHChronoLibMC(const JChronoData &chdata);

  /// Destructor
  ~DSPHChronoLibMC();

  /// Loads data for bodies and configures objects.
  void Config(std::string dirout,bool svdata,bool simulate2d);

  /// Loads inertia for bodies.
  void Config_Inertia();

  /// Compute a single timestep for each floating and moving body.
  bool RunChrono(double timestep,double dt,bool predictor);
  
  /// Saves forces for each body and link (ChronoLink_forces.csv,ChronoBody_forces.csv).
  void SaveForces();

  /// Obtains positions of Spring link.
  bool GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const;

  /// Obtains RestLength of Spring link.
  double GetSpringLinkRestLength(const std::string &linkname)const;

  /// Modifies RestLength of Spring link.
  void SetSpringLinkRestLength(const std::string &linkname,double restlength)const;

  /// Obtains center of body.
  bool GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const;
};
#endif //!DISABLE_CHRONO_OMP

#endif //!DSPHCHRONOLIB_H
