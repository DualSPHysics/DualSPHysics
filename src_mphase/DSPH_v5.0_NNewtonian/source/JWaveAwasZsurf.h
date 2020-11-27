//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2020 by Dr Jose M. Dominguez
 All rights reserved.

 DualSPHysics is an international collaboration between:
 - EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 - School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
 the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the distribution.
 * Neither the name of the DualSPHysics nor the names of its contributors may be used to endorse or promote products derived 
   from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
 SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//#############################################################################
//# Cambios:
//# =========
//# - Clase para la generacion de oleaje con absorcion activa (Didier et al 2012).
//# - Se escriben las unidades en las cabeceras de los ficheros CSV. (26-04-2018)
//#############################################################################

#ifndef _JWaveAwasZsurf_
#define _JWaveAwasZsurf_

#include "JObject.h"
#include "DualSphDef.h"
#include "JMeanValues.h"
#include "JSaveCsv2.h"

#include <string>
#include <vector>

class JXml;
class TiXmlElement;
class JLog2;
class JSphMk;
class JGaugeSystem;
class JGaugeSwl;

//##############################################################################
//# XML format in _FmtXML_WavePaddlesAwas.xml.
//##############################################################################

//##############################################################################
//# JWaveAwasZsurf
//##############################################################################
class JWaveAwasZsurf : protected JObject
{
private:
  JLog2 *Log;
  JGaugeSwl* GaugeSwl;   //-Objeto para medir la posicion de la superficie libre.
  const word MkBound;    //-Mk-bound del piston asociado.
  float Scell;           ///<Cell size: KernelSize/ScellDiv (KernelSize or KernelSize/2).

  double StartAwas;      //-Time to start AWAS correction (def=start+ramp*waveperiod).
  double Swl;            //-Still water level (free-surface water).
  bool Ele2ndOrder;      //-Calcules wave elevation using 2nd order theory.
  double GaugeX;         //-Position in X from piston to measure free-surface water (def=5*Dp).
  double GaugeXh;        //-Position in X from piston to measure free-surface water (according H value).
  double GaugeXdp;       //-Position in X from piston to measure free-surface water (according Dp value).
  double GaugeY;         //-Position in Y to measure free-surface water.
  double GaugeZmin;      //-Minimum position in Z to measure free-surface water, it must be in water (def=domain limits).
  double GaugeZmax;      //-Maximum position in Z to measure free-surface water (def=domain limits).
  double GaugeDpXml;     //-Resolution to measure free-surface water, it uses Dp*gaugedp (def=0.1).

  double CoefMassLimit;  //-Coefficient to calculate mass of free-surface (def=0.5 on 3D and 0.4 on 2D).
  byte SaveData;         //-Saves CSV with information 1:by part, 2:more info 3:by step (def=0).
  double LimitAceFactor; //-Factor to limit maximum value of acceleration, with 0 disabled (def=5).
  double CorrCoefStroke; //-Coef. de stroke maximo permitido (Drift Correction).
  double CorrCoefPeriod; //-Coef. de periodo para correccion (Drift Correction).
  double CorrPowerFunc;  //-Pontecia de funcion para recuperacion (Drift Correction).

  double CorrLimit;      //-Limite de stroke para activar correccion (CorrLimit==0 sin correccion).
  double CorrTime;       //-Tiempo para aplicar correccion.
  double CorrStart;      //-Instante en que se inicia la correccion (CorrStart!=0 correccion en curso).
  bool CorrPositive;     //-Signo de la correcion iniciada.

  double Xpos0;          //-Posicion X inicial del limite del piston con el fluido.

  //-Constants for waves.
  double GravityZ;
  double Depth;         //-Water depth.

  double GaugeDp;        //-GaugeDp=Dp*GaugeDpXml
  float MassLimit;       //-MassLimit=CoefMassLimit*MassFluid

  //-Calculated values during simulation.
  bool SvDataStep;   //-Saves information of the current step.
  double Ztarget;    //-Value of expected free-surface.
  double Zsurf;      //-Value of free-surface in fluid.

  double Xpiston;    //-Last value of piston position in space.
  double Upiston;    //-Last value of piston velocity.

  unsigned NumSave;
  jcsv::JSaveCsv2 *FileSurf;     //-Saves CSV with free-surface information (AWASFreeSurf.csv).
  jcsv::JSaveCsv2 *FileData;     //-Saves CSV with information (AWASData.csv).
  jcsv::JSaveCsv2 *FileSurfStep; //-Saves CSV with free-surface information (AWASFreeSurfStep_XXXX.csv).

  //-Variables for stadistical information of zsurf variation.
  double TimeStepM1,ZsurfM1;
  static const bool SaveVarInfo=false; //-Para usar descomentar codigo.
  unsigned VarNum;
  tdouble3 VarZsurf,VarDt,VarZsurfDt; //-Valores (mean,min,max)
  void SetStats(tdouble3 &var,double v){ var=TDouble3((var.x*VarNum)/(VarNum+1)+v/(VarNum+1),(var.y<=v? var.y: v),(var.z>=v? var.z: v)); }

  //-Saves infor per step.
  //static const unsigned SurfStepSize=1u;       //-Desactiva la grabacion de info por step.
  static const unsigned SurfStepSize=1000u;  //-Activa la grabacion de info por step.
  tdouble4 SurfStepData[SurfStepSize]; //-{timestep,xpos,zsurf,zsurftarget}
  tdouble2 SurfStepData2[SurfStepSize]; //-{LimAcePre,LimAce}
  unsigned SurfStepCount;
  unsigned SurfStepFile;

  //-Limits extreme movements.
  double LimMotionAceMax;  //-Aceleracion maxima segun periodo y amplitud de oscilacion regular.
  double LimAceMax;        //-Aceleracion maxima permitida.

  double LimVel;     //-Velocidad del ultimo movimiento.
  double LimAcePre;  //-Aceleracion del ultimo movimiento (sin corregir).
  double LimAce;     //-Aceleracion del ultimo movimiento (valor aplicado).
  JMeanValue LimAceStats; //-Datos estadisticos de aceleracion.

  void Reset();
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void RunAwas(double timestep,double ztarget,bool svdata);
  void UpdateGaugePoints();
  
  void CreateFileSurf(bool spectrum,double waveperiod,double waveheight,double wavelength,double amplitude,double initphase);
  void SaveInfoSurf();
  double GetDriftCorrection(double timecorr)const;

  void CreateFileSurfStep(unsigned numfile);
  void SaveFileSurfStep();
  void SaveFileSurfStep(double timestep,double xpos,double zsurf,double zsurftarget);

public:
  JWaveAwasZsurf(word mkbound,const JXml *sxml,TiXmlElement* lis);
  ~JWaveAwasZsurf(){ DestructorActive=true; Reset(); }
  void Init(JGaugeSystem *gaugesystem,const JSphMk *mkinfo,tdouble3 pistondir
    ,double gravityz,double depth,double wavestart,double waveheight,double waveperiod
    ,double initphase,double timeramp,double wavelength,double amplitude,double maxacepiston);
  void GetConfig(std::vector<std::string> &lines)const;
  bool GetEle2ndOrder()const{ return(Ele2ndOrder); }
  double GetGaugeX()const{ return(Xpiston+GaugeX); } //-Devuelve posicion de medida para calculo de elevacion.
  double GetGaugeXDist()const{ return(GetGaugeX()-Xpos0); } //-Devuelve distancia de la posicion de medida a la posicion inicial del piston.

  double GetPosCorrection(double timestep,double ztarget,bool svdata,double dt,double xmov,bool dgvoid=false);

  double GetDepth()const{ return(Depth); }
};


#endif

