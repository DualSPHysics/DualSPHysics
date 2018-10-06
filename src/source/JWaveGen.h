//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Incluye la aplicacion de TimeMod. (23-04-2018)
//:#############################################################################

/// \file JWaveGen.h \brief Declares the class \ref JWaveGen.

#ifndef _JWaveGen_
#define _JWaveGen_

#include "Types.h"
#include <string>

class JXml;
class JLog2;
class JWavePaddles;
class JGaugeSystem;
class JSphMk;

//##############################################################################
//# JWaveGen
//##############################################################################
/// \brief Implements wave generation for regular and irregular waves.

class JWaveGen
{
private:
  const bool UseOmp;
  const bool UseGpu;
  double TimeMod;       ///<Modifies the timestep for paddle motion.
  JWavePaddles* WavPad; 
  
  //-Vars. auxiliares cargadas tras Init().
  //-Auxiliary variables loaded after Init().
  bool Use_Awas;    //-Use of AWAS-Zsurf.
  unsigned Count;

public:

  //==============================================================================
  /// Constructor.
  //==============================================================================
  JWaveGen(bool useomp,bool usegpu,JLog2* log,std::string dirdata,JXml *sxml,const std::string &place);

  //==============================================================================
  /// Destructor.
  //==============================================================================
  ~JWaveGen();

  //==============================================================================
  /// Configura paddle con datos de las particulas.
  /// Set paddle with the particle data.
  //==============================================================================
  bool ConfigPaddle(word mkbound,word paddleid,unsigned idbegin,unsigned np);

  //==============================================================================
  /// Prepara movimiento de paddles.
  /// Prepares paddle movement.
  //==============================================================================
  void Init(JGaugeSystem *gaugesystem,const JSphMk *mkinfo,double timemax,tfloat3 gravity);

  //==============================================================================
  /// Adjust motion paddles for the first simulation instant.
  //==============================================================================
  void SetTimeMod(double timemod){ TimeMod=timemod; };

  //==============================================================================
  /// Shows object configuration using Log.
  //==============================================================================
  void VisuConfig(std::string txhead,std::string txfoot);

  //==============================================================================
  /// Devuelve si es un movimiento lineal y los datos de movimiento para el intervalo indicado.
  /// Returns if it is a linear motion and the motion data for the specified interval.
  //==============================================================================
  bool GetMotion(bool svdata,unsigned cp,double timestep,double dt
    ,tdouble3 &simplemov,tdouble3 &simplevel,tmatrix4d &matmov
    ,unsigned &np,unsigned &idbegin);

  //==============================================================================
  /// Devuelve si es un movimiento lineal y los datos de movimiento para el intervalo indicado.
  /// Returns if it is a linear motion and the motion data for the specified interval.
  //==============================================================================
  bool GetMotionAce(bool svdata,unsigned cp,double timestep,double dt
    ,tdouble3 &simplemov,tdouble3 &simplevel,tdouble3 &simpleace
    ,tmatrix4d &matmov,tmatrix4d &matmov2,unsigned &np,unsigned &idbegin);

  //==============================================================================
  /// Devuelve Mkbound de paddle.
  /// Returns Mkbound of paddle.
  //==============================================================================
  word GetPaddleMkbound(unsigned cp)const;


  unsigned GetCount()const{ return(Count); }
  bool UseAwas()const{ return(Use_Awas); } 

};


#endif


