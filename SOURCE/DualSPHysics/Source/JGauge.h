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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para medir velocidad en el fluido durante la simulacion.
//:# - Se comento la parte donde se aplica la correccion del kernel en el 
//:#   calculo de la velocidad. (04-02-2015)
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:#############################################################################

/// \file JGauge.h \brief Declares the class \ref JGauge.

#ifndef _JGauge_
#define _JGauge_

#include <string>
#include "JObject.h"
#include "Types.h"
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

class JLog2;
class JSaveCsv;

//##############################################################################
//# JGaugeBase
//##############################################################################
/// \brief Defines the common part of different gauge types.

class JGaugeBase : protected JObject
{
public:
  /// Structure with the constants of simulation for gauges.
  typedef struct{
    bool simulate2d;
    TpCellOrder cellorder;
    float massfluid;
    double dp;
    float dosh;
    float scell;
    int hdiv;
    tdouble3 domposmin;
    tdouble3 domrealposmin;
    tdouble3 domrealposmax;
  }StCtes;
  
  static StCtes SetCtes(bool simulate2d,TpCellOrder cellorder,float massfluid,double dp,float dosh,float scell,int hdiv,tdouble3 domposmin,tdouble3 domrealposmin,tdouble3 domrealposmax){
    StCtes ctes;
    ctes.simulate2d=simulate2d;
    ctes.cellorder=cellorder;
    ctes.massfluid=massfluid;
    ctes.dp=dp;
    ctes.dosh=dosh;
    ctes.scell=scell;
    ctes.hdiv=hdiv;
    ctes.domposmin=domposmin;
    ctes.domrealposmin=domrealposmin;
    ctes.domrealposmax=domrealposmax;
    return(ctes);
  }


protected:
  JLog2* Log;

  //-Constants during simulation.
  bool Simulate2D;
  TpCellOrder CellOrder;
  float MassFluid;
  double Dp;
  float Dosh,H;
  tdouble3 DomPosMin; ///<Simulation limits plus border 2h (using periodic conditions). | Limites de simulacion + borde 2h si hay condiciones periodicas.
  tdouble3 DomRealPosMin,DomRealPosMax; ///<Real simulation limits according to DomCelIni/End (without borders for periodic conditions). | Limites reales de simulacion segun DomCelIni/Fin (sin bordes de condiciones periodicas).
  float Scell;
  int Hdiv;
  float Awen;
  float Fourh2;

  JGaugeBase(JLog2* log);
  void Reset();
  void ConfigConstants(StCtes ctes);

  bool PointIsOut(double px,double py,double pz)const{ return(px!=px || py!=py || pz!=pz || px<DomRealPosMin.x || py<DomRealPosMin.y || pz<DomRealPosMin.z || px>=DomRealPosMax.x || py>=DomRealPosMax.y || pz>=DomRealPosMax.z); }
  bool PointIsOut(double px,double py)const{ return(px!=px || py!=py || px<DomRealPosMin.x || py<DomRealPosMin.y || px>=DomRealPosMax.x || py>=DomRealPosMax.y); }

 #ifdef _WITHGPU
  void RunExceptionCuda(const std::string &method,const std::string &msg,cudaError_t error);
  void CheckCudaError(const std::string &method,const std::string &msg);
 #endif

public:
  double GetDp()const{ return(Dp); }
};

//##############################################################################
//# JGaugeVel
//##############################################################################
/// \brief Calculates velocity in fluid domain.
class JGaugeVel : public JGaugeBase
{
protected:
  //-Results of measurement. | Resultados de medicion.
  double TimeStep;
  bool PtOut; ///<Indicates if the point is outside the real simulation domain. | Indica si el punto esta fuera del dominio real de simulacion.
  tdouble3 PtPos;
  tfloat3 PtVel;

  void Reset();
  void SaveResult();

public:
  JGaugeVel(JLog2* log);
  ~JGaugeVel();

  void Init(StCtes ctes);

  tfloat3 CalculeCpu(double timestep,tdouble3 ptpos,tuint3 ncells,tuint3 cellmin
    ,const unsigned *begincell,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop);

 #ifdef _WITHGPU
  tfloat3 CalculeGpu(double timestep,tdouble3 ptpos,tuint3 ncells,tuint3 cellmin,const int2 *beginendcell
    ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,float3 *aux);
 #endif
};

//##############################################################################
//# JGaugeZsurf
//##############################################################################
/// \brief Calculates surface in fluid domain.
class JGaugeZsurf : public JGaugeBase
{
protected:
  double GaugeZmin;   ///<Minimum position in Z to measure free-surface water, it must be in water (def=domain limits).
  double GaugeZmax;   ///<Maximum position in Z to measure free-surface water (def=domain limits).
  double GaugeDp;     ///<Distance between measuring positions in Z. | Espacio entre posiciones de medicion en Z.
  float MassLimit;    ///<Mass value to detect the free surface. | Valor de masa para detectar la superfice libre.

  //-Results of measurement. | Resultados de medicion.
  double TimeStep;
  bool PtOut;      ///<Indicates if the point is outside the real simulation domain. | Indica si el punto esta fuera del dominio real de simulacion.
  tdouble3 PtPos;  ///<Measurement position (X,Y,Zsurf previa). | Posicion de medicion (X,Y,Zsurf previa).
  double PtzMin;   ///<Lower point used to calculate the free surface. | Punto inferior utilizado para calcular la superficie libre.
  double Zsurf;    ///<Z position of free surface. | Posicion Z de superficie libre.

  void Reset();
  void SaveResult();
  void GetLimitsFreeSurface(double px,double py,double zsurf,tuint3 ncells,tuint3 cellmin,unsigned &czmin,unsigned &czmax)const;
  float CalculeMass(tdouble3 ptpos,tuint3 ncells,tuint3 cellmin,const unsigned *begincell
    ,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)const;

public:
  JGaugeZsurf(JLog2* log);
  ~JGaugeZsurf();

  void Init(StCtes ctes,double gaugezmin,double gaugezmax,double gaugedp,float masslimit);

  double CalculeCpu(double timestep,double px,double py,double zsurf0,tuint3 ncells,tuint3 cellmin
    ,const unsigned *begincell,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop);

#ifdef _WITHGPU
  double CalculeGpu(double timestep,double px,double py,double zsurf0,tuint3 ncells,tuint3 cellmin
    ,const int2 *beginendcell,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,double *aux);
 #endif

  double GetTimeStep()const{ return(TimeStep); }
  bool GetPtOut()const{ return(PtOut); }
  tdouble3 GetPtPos()const{ return(PtPos); }
  double GetPtzMin()const{ return(PtzMin); }
};


#endif


