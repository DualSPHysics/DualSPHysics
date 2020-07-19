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

//#############################################################################
//# Cambios:
//# =========
//# - Recibe parametros del objeto en una structuctura. (03-05-2016)
//# - Codigo adaptado al estilo de DualSPHysics. (03-05-2016)
//# - Uso de JChronoData para configurar y gestionar varios objetos con 
//#   con multiples uniones. (10-05-2016)
//# - Nuevo atributo DataDir en JChronoData. (10-05-2016)
//# - Cambio de solver e integrador para evitar valores NaN. (23-05-2016)
//# - Procesa todos los objetos de una vez. (26-05-2016)
//# - Se cambio el nombre de fomegavel por fomegaace. (27-05-2016)
//# - Funciones ApplyInitialVel y ApplyImposedVel para aplicar velocidades
//#  externas (29-06-2020).
//#############################################################################

/// \file DSPHChronoLib.h \brief Declares the class \ref DSPHChronoLib which is the interface between DualSPHysics and Chrono.

#ifndef DSPHCHRONOLIB_H
#define DSPHCHRONOLIB_H

#include "TypesDef.h"
#include "JChronoData.h"
#include <string>
#include <iostream>
#include <memory>
//ChApi const ChVector<double> VNULL(0., 0., 0.);

//-Forward declarations to avoid including chrono classes.
namespace chrono {
  class ChSystem;
  class ChSystemParallel;     //<chrono_multicore>
  class ChMaterialSurface;    //<chrono_contacts>
  class ChBody;
  namespace fea{
    class ChMesh;
  }
};

//##############################################################################
//# DSPHChronoLib
//##############################################################################
/// \brief Defines the class interface between DualSPHysics and Chrono.
class DSPHChronoLib {
public:
  //-States of execution.
  typedef enum { RSTATE_Init, RSTATE_Loading, RSTATE_Results }TpRunState;
  const std::string version;  ///<Number version
protected:
  const std::string ClassName; 

  //-Chrono physical system.
  std::string DirOut;
  bool Simulate2D;          ///<True for 2D Simulations.
  bool UseOmp;              ///<Indicates if use of ChronoEngine_Parallel module is enabled
  int OmpThreads;           ///<Threads number used by OpenMP.
  unsigned SolverIndex;     ///<Indicates the index of chrono solver enum.
  unsigned MaxIter;         ///<Indicates the maximun number of iterations for the solver.
  bool DG;                  ///<Used for Debug
  bool UseSMC;              ///<True if it is using SMC (SMooth Contacts) (chrono_contacts)
  bool UseFEA;              ///<True if the use of Finite Element Analysis is enabled
  JChronoData ChData; 
  TpRunState RunState;
  virtual void SaveForcesHead(){};
  static tdouble3 VecUnitarySafe(const tdouble3 &v);
  DSPHChronoLib(const JChronoData &chdata);

public:

  ///Initialize floating body.
  virtual void Config(std::string dirout, bool svdata, bool simulate2d){};

  ///Initialize floating body.
  virtual void Config_Inertia(){};

  ///Compute a single timestep for each floating and moving body.
  virtual bool RunChrono(double timestep, double dt, bool predictor)=0;

  ///Saves forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
  virtual void SaveForces(){};

  ///Obtains positions of Spring link.
  virtual bool GetSpringLinkPositions(const std::string &linkname, tdouble3 &p1, tdouble3 &p2)const=0;

  /// Obtains RestLength of Spring link.
  virtual double GetSpringLinkRestLength(const std::string &linkname)const=0;

  /// Modifies RestLength of Spring link.
  virtual void SetSpringLinkRestLength(const std::string &linkname, double restlength)const{};

  ///Obtains center of body.
  virtual bool GetBodyCenter(const std::string &bodyname, tdouble3 &pcen)const=0;

  ///Returns pointer to ChronoData object.
  const JChronoData* GetChronoData(){ return(&ChData); }

  ///Loads floating data to calculate coupling with Chrono.
  bool SetFtData(word mkbound, const tfloat3 &face, const tfloat3 &fomegaace);

  ///Loads imposed velocity for floating to calculate coupling with Chrono.    //<vs_fttvel>
  bool SetFtDataVel(word mkbound,const tfloat3 &vlin,const tfloat3 &vang);   //<vs_fttvel>
  
  ///Obtains floating data from coupling with Chrono.
  bool GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const;

  ///Loads motion data to calculate coupling with Chrono.
  bool SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double stepdt);
  
  /// Adds the material properties to chrono object
  void ConfigMaterial(const JChBody &body,chrono::ChBody *chbody);
    
  /// Adds the initial velocity
  void ApplyInitialVel(const JChBody &body,chrono::ChBody *chbody);

  /// Adds the imposed velocity.
  void ApplyImposedVel(const JChBodyFloating &body,chrono::ChBody *chbody);
};


//##############################################################################
//# DSPHChronoLibSC
//##############################################################################
/// \brief Defines the class for single-core executions.
class DSPHChronoLibSC:public DSPHChronoLib {
private:
  chrono::ChSystem *MphysicalSystem;
  void SaveForcesHead();

public:
  DSPHChronoLibSC(const JChronoData &chdata);
  ~DSPHChronoLibSC();

  ///Initialize floating body.
  void Config(std::string dirout,bool svdata,bool simulate2d);

  ///Initialize floating body.
  void Config_Inertia();

  ///Compute a single timestep for each floating and moving body.
  bool RunChrono(double timestep,double dt,bool predictor);

  ///Saves forces for each body and link (ChronoLink_forces.csv,ChronoBody_forces.csv).
  void SaveForces();

  ///Obtains positions of Spring link.
  bool GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const;

  /// Obtains RestLength of Spring link.
  double GetSpringLinkRestLength(const std::string &linkname)const;

  /// Modifies RestLength of Spring link.
  void SetSpringLinkRestLength(const std::string &linkname,double restlength)const;

  ///Obtains center of body.
  bool GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const;
};


#ifndef DISABLE_CHRONO_OMP
//##############################################################################
//# DSPHChronoLibMC
//##############################################################################
/// \brief Defines the class for multi-core executions.
class DSPHChronoLibMC:public DSPHChronoLib {
private:
  chrono::ChSystemParallel *MphysicalSystem;  //<chrono_multicore> //TODO: SMC multicore
  void SaveForcesHead();

public:
  DSPHChronoLibMC(const JChronoData &chdata);
  ~DSPHChronoLibMC();

  ///Initialize floating body.
  void Config(std::string dirout,bool svdata,bool simulate2d);

  ///Initialize floating body.
  void Config_Inertia();

  ///Compute a single timestep for each floating and moving body.
  bool RunChrono(double timestep,double dt,bool predictor);

  ///Saves forces for each body and link (ChronoLink_forces.csv,ChronoBody_forces.csv).
  void SaveForces();

  ///Obtains positions of Spring link.
  bool GetSpringLinkPositions(const std::string &linkname,tdouble3 &p1,tdouble3 &p2)const;

  /// Obtains RestLength of Spring link.
  double GetSpringLinkRestLength(const std::string &linkname)const;
 
  /// Modifies RestLength of Spring link.
  void SetSpringLinkRestLength(const std::string &linkname,double restlength)const;

  ///Obtains center of body.
  bool GetBodyCenter(const std::string &bodyname,tdouble3 &pcen)const;
};
#endif //!DISABLE_CHRONO_OMP

#endif //!DSPHCHRONOLIB_H