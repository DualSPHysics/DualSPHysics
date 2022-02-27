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
//:# - Incluye la aplicacion de TimeMod. (23-04-2018)
//:# - Ahora devuelve los datos de motion como StMotionData. (11-09-2019)
//:# - Permite compilar sin libreria WaveGen. (10-12-2019)
//:# - Usa el valor de gravity de la simuacion. (03-02-2020)
//:# - Permite configurar varios <savemotion>. (04-02-2020)
//:# - <savemotion> usa por defecto TimeMax y TimePart de la simulacion. (04-02-2020)
//:# - Comprueba opcion active en elementos de primer, segundo nivel y AWAS. (19-03-2020)  
//:# - Devuelve movimiento MOTT_None para tiempos fuera del intervalo de generacion. (03-02-2021)
//:# - Indica tipos de generacion de olas configurados. (23-07-2021)
//:#############################################################################

/// \file JWaveGen.h \brief Declares the class \ref JWaveGen.

#ifndef _JWaveGen_
#define _JWaveGen_

#include "DualSphDef.h"
#include <string>

//#define DISABLE_WAVEGEN     ///<It allows compile without WaveGen library.

class JXml;
class JLog2;
class JWavePaddles;
class JGaugeSystem;
class JSphMk;

//##############################################################################
//# JWaveGen
//##############################################################################
/// \brief Implements wave generation for regular and irregular waves.

#ifdef DISABLE_WAVEGEN
#include "JWaveGenUndef.h"
#else
class JWaveGen
{
private:
  const bool UseOmp;
  const bool UseGpu;
  StMotionData MotionNull; ///<Motion data for null movement.

  double TimeMod;          ///<Modifies the timestep for paddle motion.
  JWavePaddles* WavPad; 

  //-Auxiliary variables loaded after Init().
  unsigned Count;

  bool Use_Awas;       ///<Use of AWAS-Zsurf.
  bool Waves_Regular;  ///<Regular waves configured.
  bool Waves_Spectrum; ///<Irregular waves configured.
  bool Waves_File;     ///<Waves from external file configured.
  bool Waves_Solitary; ///<Solitary waves configured.

public:

  //==============================================================================
  /// Constructor.
  //==============================================================================
  JWaveGen(bool useomp,bool usegpu,JLog2 *log,std::string dirdata,const JXml *sxml
    ,const std::string &place,tdouble3 gravity3);

  //==============================================================================
  /// Destructor.
  //==============================================================================
  ~JWaveGen();

  //==============================================================================
  /// Returns true when this feature is available.
  //==============================================================================
  static bool Available(){ return(true); }

  //==============================================================================
  /// Configura paddle con datos de las particulas.
  /// Set paddle with the particle data.
  //==============================================================================
  bool ConfigPaddle(word mkbound,word motionref,unsigned idbegin,unsigned np);

  //==============================================================================
  /// Prepara movimiento de paddles.
  /// Prepares paddle movement.
  //==============================================================================
  void Init(JGaugeSystem *gaugesystem,const JSphMk *mkinfo,double timemax,double timepart);

  //==============================================================================
  /// Adjust motion paddles for the first simulation instant.
  //==============================================================================
  void SetTimeMod(double timemod){ TimeMod=timemod; };

  //==============================================================================
  /// Shows object configuration using Log.
  //==============================================================================
  void VisuConfig(std::string txhead,std::string txfoot);

  //==============================================================================
  /// Devuelve datos de movimiento nulo.
  /// Returns null motion data.
  //==============================================================================
  const StMotionData& GetMotionNull()const{ return(MotionNull); };

  //==============================================================================
  /// Devuelve datos del movimiento lineal o matricial para el intervalo indicado.
  /// Returns linear or matrix motion data for the specified interval.
  //==============================================================================
  const StMotionData& GetMotion(bool svdata,unsigned cp,double timestep,double dt);

  //==============================================================================
  /// Devuelve datos del movimiento lineal o matricial para el intervalo indicado.
  /// Returns linear or matrix motion data for the specified interval.
  //==============================================================================
  const StMotionData& GetMotionAce(bool svdata,unsigned cp,double timestep,double dt);

  //==============================================================================
  /// Devuelve Mkbound de paddle.
  /// Returns Mkbound of paddle.
  //==============================================================================
  word GetPaddleMkbound(unsigned cp)const;

  unsigned GetCount()const{ return(Count); }

  bool UseAwas()      const{ return(Use_Awas);       } 
  bool WavesRegular() const{ return(Waves_Regular);  } 
  bool WavesSpectrum()const{ return(Waves_Spectrum); } 
  bool WavesFile()    const{ return(Waves_File);     } 
  bool WavesSolitary()const{ return(Waves_Solitary); } 

};
#endif

#endif


