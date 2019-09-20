#ifndef DSPHCHRONOLIB_H
#define DSPHCHRONOLIB_H

//#############################################################################
//# Cambios:
//# =========
//# - Recibe parametros del objeto en una structuctura. (v0.001 / 03-05-2016)
//# - Codigo adaptado al estilo de DualSPHysics. (v0.001 / 03-05-2016)
//# - Uso de JChronoData para configurar y gestionar varios objetos con 
//#   con multiples uniones. (v0.002 / 10-05-2016)
//# - Nuevo atributo DataDir en JChronoData. (v0.003 / 10-05-2016)
//# - Cambio de solver e integrador para evitar valores NaN. (v0.004 / 23-05-2016)
//# - Procesa todos los objetos de una vez. (v0.005 / 26-05-2016)
//# - Se cambio el nombre de fomegavel por fomegaace. (v0.006 / 27-05-2016)
//#############################################################################

#include "TypesDef.h"
#include "JChronoData.h"
#include <string>

//ChApi const ChVector<double> VNULL(0., 0., 0.);


//-Forward declarations to avoid including chrono classes.
namespace chrono{
  class ChSystemNSC;
};


class DSPHChronoLib{
 public:
  //-States of execution.
  typedef enum{ RSTATE_Init,RSTATE_Loading,RSTATE_Results }TpRunState; 
  const std::string version;

 private:
  //-Chrono physical system.
  chrono::ChSystemNSC *mphysicalSystem;

  std::string dirOut;
  bool Simulate2D;
  JChronoData chData;
  TpRunState RunState;

  void SaveForcesHead();
  static tdouble3 VecUnitarySafe(const tdouble3 &v);

 public:


  DSPHChronoLib();
  ~DSPHChronoLib();

  ///Initialize floating body.
  void Config(std::string dirout,bool svdata,bool simulate2d,const JChronoData &chdata);

  ///Initialize floating body.
  void Config_Inertia();

  ///Returns pointer to ChronoData object.
  const JChronoData* GetChronoData(){ return(&chData); }

  ///Loads floating data to calculate coupling with Chrono.
  bool SetFtData(word mkbound,const tfloat3 &face,const tfloat3 &fomegaace);

  ///Obtains floating data from coupling with Chrono.
  bool GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const;

  ///Loads motion data to calculate coupling with Chrono.
  bool SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double stepdt);

  ///Compute a single timestep for each floating and moving body.
  bool RunChrono(double timestep,double dt,bool predictor);

  ///Saves forces for each body and link (ChronoLink_forces.csv, ChronoBody_forces.csv).
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

#endif
