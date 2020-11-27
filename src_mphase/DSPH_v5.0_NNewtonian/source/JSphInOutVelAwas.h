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

/// \file JSphInOutVelAwas.h \brief Declares the class \ref JSphInOutVelAwas.

#ifndef _JSphInOutVelAwas_
#define _JSphInOutVelAwas_

#include "JObject.h"
#include "DualSphDef.h"
#include "JSphInOutDef.h"
#include <string>
#include <vector>

class JXml;
class TiXmlElement;
class JLog2;
class JLinearValue;
class JGaugeSystem;
class JGaugeSwl; 


//##############################################################################
//# JSphInOutVelAwas
//##############################################################################
/// \brief AWAS implementation to calculate velocity correction according to SWL.
class JSphInOutVelAwas : protected JObject
{
private:
  JLog2 *Log;

  const unsigned IdZone;
  const double InletX;       ///<Inlet limit in X.
  const tdouble3 InletDir;   ///<Inlet direction in X.
  const float GravityZ;

  //-Initial configuration variables.
  bool InletMode;            ///<Applies Inlet correction or Outlet correction.
  double StartAwas;          ///<Time to start AWAS correction (def=start+ramp*waveperiod).
  double InitDepth;          ///<Initial depth.
  double CoefDepth;          ///<Coefficient from initial depth. CoefDepth=sqrt(-GravityZ/InitDepth)
  JLinearValue *ZsurfTarget; ///<Zsurf target to compute AWAS correction.
  double GaugeX;             ///<Position in X from piston to measure free-surface water (def=5*Dp).
  double GaugeXh;            ///<Position in X from piston to measure free-surface water (according H value).
  double GaugeXdp;           ///<Position in X from piston to measure free-surface water (according Dp value).
  double GaugeY;             ///<Position in Y to measure free-surface water.
  double GaugeZmin;          ///<Minimum position in Z to measure free-surface water, it must be in water (def=domain limits).
  double GaugeZmax;          ///<Maximum position in Z to measure free-surface water (def=domain limits).
  double GaugeDpXml;         ///<Resolution to measure free-surface water, it uses Dp*gaugedp (def=0.1).
  double GaugeDp;            ///<Gauge resolution. GaugeDp=Dp*GaugeDpXml
  byte SaveData;             ///<Saves CSV with AWAS information. 0:none, 1:PART data, 2:all steps.

  JGaugeSwl *GaugeSwl;       ///<Gauge object to measure water level in front of the inlet/outlet zone.

  //-Saves step data for CSV.
  static const unsigned SizeStepData=200;
  unsigned CountStepData;
  tfloat4 *StepData;         ///<Saves values of each step. [SizeSaveData]

  //-Values for the current step.
  double LastTimeStep;
  float LastZgauge;
  float LastZtarget;
  float LastVelCorr;

private:
  void ReadXml(const JXml *sxml,TiXmlElement* lis,const std::string &dirdatafile
    ,JGaugeSystem *gaugesystem);

public:
  JSphInOutVelAwas(unsigned idzone,double inletx,tdouble3 inletdir
    ,float gravityz,const std::string &dirdatafile,JGaugeSystem *gaugesystem
    ,const JXml *sxml,TiXmlElement* ele);
  ~JSphInOutVelAwas();
  void Reset();
  void GetConfig(std::vector<std::string> &lines)const;
  float GetVelCorr(double timestep);
  void SaveCsvData();
};

#endif


