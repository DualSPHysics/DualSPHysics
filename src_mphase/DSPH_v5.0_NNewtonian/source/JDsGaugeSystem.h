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
//:# - Clase para gestionar la medicion de distintas magnitudes de forma 
//:#   automatica y simple. (12-02-2017)
//:# - Error corregido cargando <default><output>. (03-03-2017)
//:# - Nueva opcion para calcular fuerzas sobre fixed o moving boundary. (20-11-2018)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:# - Objeto JXml pasado como const para operaciones de lectura. (18-03-2020)  
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (18-03-2020)  
//:# - Cambio de nombre de fichero J.GaugeSystem a J.DsGaugeSystem. (28-06-2020)
//:# - Nuevos metodos CalculeLastInputXXX(). (27-08-2020)
//:#############################################################################

/// \file JDsGaugeSystem.h \brief Declares the class \ref JGaugeSystem.

#ifndef _JDsGaugeSystem_
#define _JDsGaugeSystem_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#include "JDsGaugeItem.h"

class JXml;
class TiXmlElement;
class JLog2;
class JSphMk;

//##############################################################################
//# XML format in _FmtXML_Gauges.xml.
//##############################################################################

//##############################################################################
//# JGaugeSystem
//##############################################################################
/// \brief Manages the list of configured \ref gauges.

class JGaugeSystem : protected JObject
{
public:
  ///Structure to save execution input on CPU.
  typedef struct StrInputCpu{
    bool ready;
    double timestep;
    StDivDataCpu dvd;
    unsigned npbok;
    unsigned npb;
    unsigned np;
    const tdouble3 *pos;
    const typecode *code;
    const unsigned *idp;
    const tfloat4 *velrhop;
    StrInputCpu(){ ready=false; npbok=npb=np=0; }
    StrInputCpu(double timestep_,const StDivDataCpu &dvd_
      ,unsigned npbok_,unsigned npb_,unsigned np_,const tdouble3 *pos_
      ,const typecode *code_,const unsigned *idp_,const tfloat4 *velrhop_)
    {
      ready=true;  timestep=timestep_;
      dvd=dvd_;  npbok=npbok_;  npb=npb_;   np=np_;  
      pos=pos_;  code=code_;  idp=idp_;  velrhop=velrhop_;
    }
  }StInputCpu;


 #ifdef _WITHGPU
  ///Structure to save execution input on GPU.
  typedef struct StrInputGpu{
    bool ready;
    double timestep;
    StDivDataGpu dvd;
    unsigned npbok;
    unsigned npb;
    unsigned np;
    const double2 *posxy;
    const double *posz;
    const typecode *code;
    const unsigned *idp;
    const float4 *velrhop;
    StrInputGpu(){ ready=false; npbok=npb=np=0; }
    StrInputGpu(double timestep_,const StDivDataGpu &dvd_
      ,unsigned npbok_,unsigned npb_,unsigned np_,const double2 *posxy_,const double *posz_
      ,const typecode *code_,const unsigned *idp_,const float4 *velrhop_)
    {
      ready=true;  timestep=timestep_;
      dvd=dvd_;  npbok=npbok_;  npb=npb_;   np=np_;  
      posxy=posxy_;  posz=posz_;  code=code_;  idp=idp_;  velrhop=velrhop_;
    }
  }StInputGpu;
 #endif

private:
  JLog2* Log;
  const bool Cpu;

  bool Configured;
  StCteSph CSP;           ///<Structure with main SPH constants values and configurations.
  bool Symmetry;          ///<Use of symmetry in plane y=0.
  double TimeMax;         ///<Total time to simulate [s].
  double TimePart;        ///<Time of output data [s].
  tdouble3 DomPosMin;     ///<Lower limit of simulation + edge (KernelSize) if periodic conditions. DomPosMin=Map_PosMin+(DomCelIni*Scell); | Limite inferior de simulacion + borde (KernelSize) si hay condiciones periodicas. 
  tdouble3 DomPosMax;     ///<Upper limit of simulation + edge (KernelSize) if periodic conditions. DomPosMax=min(Map_PosMax,Map_PosMin+(DomCelFin*Scell)); | Limite inferior de simulacion + borde (KernelSize) si hay condiciones periodicas. 
  float Scell;            ///<Cell size: KernelSize/ScellDiv (KernelSize or KernelSize/2).
  int ScellDiv;           ///<Value to divide KernelSize (1 or 2).


  JGaugeItem::StDefault CfgDefault; ///<Default configuration.
  
  std::vector<JGaugeItem*> Gauges;

  //-Variables for GPU.
 #ifdef _WITHGPU
  float3* AuxMemoryg;  ///<Auxiliary allocated memory on GPU [1].
 #endif

  //-Variables for saving input state.
  StInputCpu InputCpu;
 #ifdef _WITHGPU
  StInputGpu InputGpu;
 #endif

  void ResetCfgDefault();
  void LoadLinePoints(double coefdp,const tdouble3 &point1,const tdouble3 &point2,std::vector<tdouble3> &points,const std::string &ref)const;
  void LoadLinePoints(unsigned count,const tdouble3 &point1,const tdouble3 &point2,std::vector<tdouble3> &points,const std::string &ref)const;
  void LoadPoints(JXml *sxml,TiXmlElement* lis,std::vector<tdouble3> &points)const;
  JGaugeItem::StDefault ReadXmlCommon(const JXml *sxml,TiXmlElement* ele)const;
  void ReadXml(const JXml *sxml,TiXmlElement* ele,const JSphMk* mkinfo);

public:
  JGaugeSystem(bool cpu);
  ~JGaugeSystem();
  void Reset();

  void Config(const StCteSph &csp,bool symmetry,double timemax,double timepart
    ,tdouble3 posmin,tdouble3 posmax,float scell,int scelldiv);

  void LoadXml(const JXml *sxml,const std::string &place,const JSphMk* mkinfo);
  void VisuConfig(std::string txhead,std::string txfoot);

  bool GetSimulate2D()const{ return(CSP.simulate2d); };
  double GetSimulate2DPosY()const{ return(CSP.simulate2dposy); };
  double GetDp()const{ return(CSP.dp); }
  tdouble3 GetDomPosMin()const{ return(DomPosMin); }
  tdouble3 GetDomPosMax()const{ return(DomPosMax); }
  float GetMassFluid()const{ return(CSP.massfluid); }
  //float GetMassBound()const{ return(MassBound); }
  float GetKernelH()const{ return(CSP.kernelh); }
  float GetScell()const{ return(Scell); }

  void LoadLinePoints(double coefdp,const tdouble3 &point1,const tdouble3 &point2,std::vector<tdouble3> &points)const{ LoadLinePoints(coefdp,point1,point2,points,""); }
  void LoadLinePoints(unsigned count,const tdouble3 &point1,const tdouble3 &point2,std::vector<tdouble3> &points)const{ LoadLinePoints(count,point1,point2,points,""); }

  JGaugeVelocity* AddGaugeVel  (std::string name,double computestart,double computeend,double computedt
    ,const tdouble3 &point);
  JGaugeSwl*      AddGaugeSwl  (std::string name,double computestart,double computeend,double computedt
    ,tdouble3 point0,tdouble3 point2,double pointdp,float masslimit=0);
  JGaugeMaxZ*     AddGaugeMaxZ (std::string name,double computestart,double computeend,double computedt
    ,tdouble3 point0,double height,float distlimit);
  JGaugeForce*    AddGaugeForce(std::string name,double computestart,double computeend,double computedt
    ,const JSphMk* mkinfo,word mkbound);

  void SaveVtkInitPoints()const;

  unsigned GetCount()const{ return(unsigned(Gauges.size())); }
  unsigned GetGaugeIdx(const std::string &name)const;
  JGaugeItem* GetGauge(unsigned c)const;

  void CalculeCpu(double timestep,const StDivDataCpu &dvd
    ,unsigned npbok,unsigned npb,unsigned np,const tdouble3 *pos
    ,const typecode *code,const unsigned *idp,const tfloat4 *velrhop
    ,bool saveinput=false);

  void CalculeLastInputCpu(std::string gaugename);

 #ifdef _WITHGPU
  void CalculeGpu(double timestep,const StDivDataGpu &dvd
    ,unsigned npbok,unsigned npb,unsigned np,const double2 *posxy,const double *posz
    ,const typecode *code,const unsigned *idp,const float4 *velrhop
    ,bool saveinput=false);

  void CalculeLastInputGpu(std::string gaugename);
 #endif

  void CalculeLastInput(std::string gaugename);

  void SaveResults(unsigned cpart);
};


#endif


