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

/// \file JSphInOutVel.h \brief Declares the class \ref JSphInOutVel.

#ifndef _JSphInOutVel_
#define _JSphInOutVel_

#include "JObject.h"
#include "DualSphDef.h"
#include "JSphInOutDef.h"

#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
  #include "JSphTimersGpu.h"
#endif

#include <string>
#include <vector>
 
class JXml;
class TiXmlElement;
class JLog2;
class JLinearValue;
class JSphInOutVelAwas;
class JSphInOutGridData;
class JGaugeSystem;


//##############################################################################
//# JSphInOutVel
//##############################################################################
/// \brief Manages the velocity data for inlet/outlet conditions.

class JSphInOutVel : protected JObject
{
protected:
  JLog2 *Log;
  const bool Cpu;
  const StCteSph CSP; ///<Structure with main SPH constants values and configurations.
  const unsigned IdZone;
  const tdouble3 Direction;   ///<Inflow direction.
  const tdouble3 PtPlane;     ///<Position to create inlet plane.
  const tdouble3 ZonePosMin;  ///<Minimum position of inlet points.
  const tdouble3 ZonePosMax;  ///<Maximum position of inlet points.

  TpInVelMode VelMode;        ///<Inflow velocity mode (fixed, variable, extrapolated or interpolated).
  TpInVelProfile VelProfile;  ///<Inflow velocity profile (uniform, linear or parabolic).
  TpInBehaviour VelBehaviour; ///<Behaviour of inlet/outlet according to the velocity.

  float VelMin;        ///<Minimum input velocity or -FLT_MAX (when it is unknown).
  float VelMax;        ///<Miximum input velocity or +FLT_MAX (when it is unknown).

  double CurrentTime;     ///<Timestep of current velocity data.
  tfloat4 CurrentCoefs0;  ///<Current velocity coefficients for fiexed or variable velocity.
  tfloat4 CurrentCoefs1;  ///<Current velocity coefficients for fiexed or variable velocity.


  //-Velocity definition for Fixed and Variable velocity.
  float InputVel;             ///<Input velocity (used when VelMode==InVelM_Fixed).
  float InputVel2;            ///<2nd input velocity (used when VelMode==InVelM_Fixed).
  float InputVel3;            ///<3rd input velocity (used when VelMode==InVelM_Fixed).
  float InputVelPosz;         ///<1st input velocity Z (used when VelMode==InVelM_Fixed).
  float InputVelPosz2;        ///<2nd input velocity Z (used when VelMode==InVelM_Fixed).
  float InputVelPosz3;        ///<3rd input velocity Z (used when VelMode==InVelM_Fixed).

  JLinearValue *InputTimeVel; ///<Input velocity in time (for VelMode==InVelM_Variable).
  unsigned TimeVelIdx0;
  unsigned TimeVelIdx1;
  
  //-Velocity according flow [l/s].
  bool FlowActive;
  float FlowRatio;       ///<It is used to compute particle volume (ratio*dp^3). 
  unsigned FlowPointsOk; ///<Initial active inlet points (below zsurf). 
  float FlowToVel;       ///<Factor to compute velocity [m/s] starting from flow [l/s]. 

  //unsigned TimeVelIdx0,TimeVelIdx1; ///<Interval of time for velocity variable.
  bool SaveVelProfile;              ///<Indicate when input velocity profile must be saved in CSV file.

  JSphInOutVelAwas* AwasVel;     ///<AWAS object to compute velocity.

  JSphInOutGridData *InputVelGrid;  ///<Input velocity is interpolated in time from a data grid (for VelMode==InVelM_Interpolated).

protected:
  void CalculateVelMinMax(float &velmin,float &velmax)const;
  void DefineBehaviour();
  void ComputeInitialVel();
  void UpdateVelVariable(double timestep);


public:
  JSphInOutVel(bool cpu,unsigned idzone,const StCteSph &csp,tdouble3 direction
    ,tdouble3 ptplane,tdouble3 zoneposmin,tdouble3 zoneposmax);
  ~JSphInOutVel();
  void Reset();
  TpInVelMode ReadXml(const JXml *sxml,TiXmlElement *lis,const std::string &dirdatafile
    ,JGaugeSystem *gaugesystem,double maprealposminy);

  bool GetFlowActive()const{ return(FlowActive); }
  void ConfigFlowToVel(unsigned initnptok);

  void GetConfig(std::vector<std::string> &lines)const;

  float GetVelMin()const{ return(VelMin); }
  float GetVelMax()const{ return(VelMax); }
  TpInBehaviour GetInletBehaviour()const{ return(VelBehaviour); }
  std::string GetInletBehaviourName()const;
  bool UseAwasVel()const{ return(AwasVel!=NULL); }

  TpInVelProfile GetVelProfile()const{ return(VelProfile); }
  bool UseCoefficients()const{ return(VelMode==InVelM_Fixed || VelMode==InVelM_Variable); }
  tfloat4 GetCurrentCoefs0()const{ return(CurrentCoefs0); }
  tfloat4 GetCurrentCoefs1()const{ return(CurrentCoefs1); }

  void UpdateVel(double timestep);

  void UpdateVelInterpolateCpu(double timestep,unsigned nplist,const int *plist
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop);

#ifdef _WITHGPU
  void UpdateVelInterpolateGpu(double timestep,unsigned nplist,const int *plist
    ,const double2 *posxyg,const double *poszg,const typecode *codeg
    ,const unsigned *idpg,float4 *velrhopg);
#endif


  void SaveAwasVelCsv();
  void SaveVtkVelGrid();

};

#endif


