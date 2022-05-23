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
//:# - Acoplamiento con la version Chrono-6.0.0 (v6.01 / 10-03-2021)
//:# - Se permite la simulacion de vigas de tipo Euler (v6.02 / 17-03-2021).
//:# - Creacion de vigas manualmente a partir de una lista de nodos (v6.03 / 29-03-2021).
//:# - Creacion de vigas manualmente entre dos puntos a una distancia de Dp (v6.04 / 05-04-2021).
//:# - Calculo de rotationes de cada nodo FEA  (v6.05 / 12-04-2021).
//:# - Velocidades iniciales lineales y angulares para FEA (v6.06 / 19-04-2021).
//:# - Funcion para aplicar un desplazamiento dado a cada nodo (v6.07 / 08-05-2021).
//:# - FEA elementos como cuerpos flotantes (v6.08 / 15-05-2021).
//:# - Comprueba si las velocidades lineales y angulares fueron inicializadas para los FEA (v6.09 / 17-05-2021).
//:# - Agrega un multiplicador para la rigidez axial para FEA (v6.10 / 18-06-2021).
//:# - Resuelto error en la deteccion decolisiones para FEA (v6.11/ 18-06-2021).
//:# - Renombra BuilderFEA por FeaBeam (v6.12/ 29-06-2021).
//:# - Construye vigas con nodos compartidos y se utiliza una unica malla (v6.13/ 29-06-2021).
//:# - Colisiones mesh vs mesh habilitadas para FEA (v6.14 / 01-07-2021).
//:# - Link_pointframe permite asignar directamente nodos a un objeto rigido (v6.15 / 01-07-2021).
//:# - Permite el uso de varios hilos con OpenMP para resolver colisiones (v6.16 / 20-07-2021).
//:# - Habilita el uso de FtPause con Chrono (v6.17 / 21-07-2021).
//:# - Se mantiene independencia entre la masa del objeto SPH y de la viga FEA (v6.18 / 22-09-2021)
//:# - Permite la opcion modified Newton matrix para evaluarlo y ensamblarlo en cada paso (HHT) (v6.19 / 19-11-2021).
//:# - Almacena los nodos de cada beam como copias en vez de punteros (v6.20 / 22-11-2021).
//:# - Permite configurar los parametros para el integrador HHT usando una estructura StHHT (v6.21 / 03-01-2022).
//:# - Acomplamiento con chrono-4.0.0 debido a los problemas de rendimiento de v6 (v4.21 / 07-04-2022).
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
  class ChMaterialSurface;
  class ChBody;
};

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
  const std::string ChronoVersion; ///<Chrono version

protected:
  const std::string ClassName;

  //-Chrono physical system.
  std::string DirOut;
  bool Simulate2D;          ///<True for 2D Simulations.
  bool UseOmp;              ///<Indicates if use of ChronoEngine_Multicore module is enabled.
  int OmpThreads;           ///<Threads number used by OpenMP.
  unsigned SolverIx;        ///<Indicates the index of chrono solver enum.
  unsigned TimeStepperIx;   ///<Indicates the index of chrono timestepper enum.
  bool UseSMC;              ///<True if it is using SMC (SMooth Contacts) 
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

  /// Establishes the variable coefficients to the link objects.
  virtual void SetVariableCoeff();

  /// Saves header for forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  void SaveForcesHead();
  
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
#endif //!DSPHCHRONOLIB_H
