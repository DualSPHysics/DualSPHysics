//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

//:NO_COMENTARIO
//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para gestionar nuevas condiciones inlet/outlet. (08-04-2017)
//:# - Se divide en dos ficheros independientes. (26-04-2020)
//:#############################################################################

/// \file JSphInOutZone.h \brief Declares the class \ref JSphInOutZone.

#ifndef _JSphInOutZone_
#define _JSphInOutZone_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
  #include "JSphTimersGpu.h"
#endif

class JXml;
class TiXmlElement;
class JLog2;
class JLinearValue;
class JSphInOutPoints;
class JSphInOutGridData;
class JWaveTheoryReg;
class JDsPartsInit;
class JGaugeSystem;
class JGaugeSwl;

//##############################################################################
//# XML format in _FmtXML_InOut.xml.
//##############################################################################

//##############################################################################
//# JSphInOutZone
//##############################################################################
/// \brief Manages one inlet zone.
class JSphInOutZone : protected JObject
{
public:

  ///Defines the treatment of fluid particles entering a inlet/outlet zone.
  typedef enum{ 
     TIN_Free=0     ///<Free input mode (no changes).
    ,TIN_Convert=1  ///<Convert fluid to inlet/outlet.
    ,TIN_Remove=2   ///<Remove fluid.
  }TpInput;

  ///Defines refilling modes.
  typedef enum{ 
     TFI_SimpleFull=0   ///<It creates new inout particles when inout particles is moved to normal fluid.
    ,TFI_SimpleZsurf=1  ///<It only creates new inout particles below zsurf.
    ,TFI_Advanced=2     ///<Advanced for reverse flows (slow).
  }TpRefilling;

  ///Controls imposed velocity.
  typedef enum{ 
    MVEL_Fixed=0x00,        ///<Imposed fixed velocity (xxxx 00xx).
    MVEL_Variable=0x04,     ///<Imposed variable velocity (xxxx 01xx).
    MVEL_Extrapolated=0x08, ///<Extrapolated from ghost nodes (xxxx 10xx).
    MVEL_Interpolated=0x0C, ///<Interpolated velocity (xxxx 11xx).
    MVEL_MASK=0x0C          ///<Mask to obtain value (0000 1100).
  }TpVelMode;   

  ///Controls profile of imposed velocity.
  typedef enum{ 
    PVEL_Constant=0x00,    ///<Imposed velocity profile constant (xxxx xx00).
    PVEL_Linear=0x01,      ///<Imposed velocity profile linear (xxxx xx01).
    PVEL_Parabolic=0x02,   ///<Imposed velocity profile parabolic (xxxx xx10).
    PVEL_MASK=0x03         ///<Mask to obtain value (0000 0011).
  }TpVelProfile;

  ///Controls imposed rhop.
  typedef enum{ 
    MRHOP_Constant=0x00,     ///<Imposed rhop profile constant (xx00 xxxx).
    MRHOP_Hydrostatic=0x10,  ///<Imposed rhop profile hydrostatic (xx01 xxxx).
    MRHOP_Extrapolated=0x20, ///<Extrapolated from ghost nodes (xx10 xxxx).
    MRHOP_MASK=0x30          ///<Mask to obtain value (0011 0000).
  }TpRhopMode;   

  ///Controls imposed Z surface.
  typedef enum{ 
    ZSURF_Undefined=0,    ///<Imposed undefined zsurf.
    ZSURF_Fixed=1,        ///<Imposed fixed zsurf.
    ZSURF_Variable=2,     ///<Imposed variable zsurf.
    ZSURF_Calculated=3,   ///<Zsurf is calculated from fluid domain.
    ZSURF_WaveTheory=4    ///<Zsurf is calculated from wave theory.
  }TpZsurfMode;   

  static const unsigned CheckInput_MASK  =0x80;  ///<Mask to obtain value (1000 0000).

  static const unsigned RefillAdvanced_MASK=0x01;  ///<Mask to obtain value (0000 0001).
  static const unsigned RefillSpFull_MASK  =0x02;  ///<Mask to obtain value (0000 0010).
  static const unsigned RemoveInput_MASK   =0x04;  ///<Mask to obtain value (0000 0100).
  static const unsigned RemoveZsurf_MASK   =0x08;  ///<Mask to obtain value (0000 1000).
  static const unsigned ConvertInput_MASK  =0x10;  ///<Mask to obtain value (0001 0000).

private:
  const bool Cpu;
  JLog2 *Log;
  const unsigned IdZone;
  const bool Simulate2D;         ///<Indicates 2D simulation.
  const double Simulate2DPosY;   ///<Y value in 2D simulations.
  const double Dp;               ///<Distance between particles.
  const tdouble3 MapRealPosMin;
  const tdouble3 MapRealPosMax;
  const float GravityZ;;

  //-Configuration parameters.
  JSphInOutPoints *Points;  ///<Definition of inlet points.
  byte Layers;              ///<Number of inlet particle layers.
  TpInput InputMode;        ///<Defines the treatment of fluid particles entering a inlet/outlet zone (default=0).
  bool InputCheck;          ///<Checks fluid input for RemoveZsurf=true or InputMode!=TIN_Free.
  TpRefilling RefillingMode;  ///<Refilling mode. 0:Simple full, 1:Simple below zsurf, 2:Advanced for reverse flows (very slow) (default=1).

  //-Domain data from Points object. PtDom[8] is the reference point in inout plane.
  tdouble3 PtDom[10];
  tfloat3 BoxLimitMin;
  tfloat3 BoxLimitMax;
  
  TpVelMode VelMode;        ///<Inflow velocity mode (fixed or variable).
  TpVelProfile VelProfile;  ///<Inflow velocity profile (constant, linear or parabolic).
  float InputVel;           ///<Input velocity (used when VelMode==MVEL_Fixed).
  float InputVel2;          ///<2nd input velocity (used when VelMode==MVEL_Fixed).
  float InputVel3;          ///<3rd input velocity (used when VelMode==MVEL_Fixed).
  float InputVelPosz;       ///<1st input velocity Z (used when VelMode==MVEL_Fixed).
  float InputVelPosz2;      ///<2nd input velocity Z (used when VelMode==MVEL_Fixed).
  float InputVelPosz3;      ///<3rd input velocity Z (used when VelMode==MVEL_Fixed).

  JLinearValue *TimeInputVelData;   ///<Input velocity in time (for VelMode==MVEL_Variable).
  unsigned TimeVelIdx0,TimeVelIdx1; ///<Interval of time for velocity variable.
  bool SaveVelProfile;              ///<Indicate when input velocity profile must be saved in CSV file.

  JSphInOutGridData *InputVelGrid;  ///<Input velocity is interpolated in time from a data grid (for VelMode==MVEL_Interpolated).
  bool ResetZVelGrid;               ///<Reset Z velocity after interaction (default=false).

  float VelMin;        ///<Minimum input velocity or -FLT_MAX (when it is unknown).
  float VelMax;        ///<Miximum input velocity or +FLT_MAX (when it is unknown).

  TpZsurfMode ZsurfMode;            ///<Inflow zsurf mode (fixed or variable).
  float InputZbottom;               ///<Bottom level of water (it is constant).
  float InputZsurf;                 ///<Input zsurf. It is FLT_MAX when ZsurfMode==ZSURF_Undefined and it changes when ZsurfMode==ZSURF_Variable.
  JLinearValue *TimeInputZsurf;     ///<Input zsurf in time (for ZsurfMode==ZSURF_Variable).
  bool RemoveZsurf;                 ///<Removes particles above the Zsurf limit (default=false)
  bool SvVtkZsurf;                  ///<Creates VTK files with Zsurf for each PART (default=false)

  bool ExternalVarInput;            ///<Receives velocity and zsurf values from external sources. (default=false)

  JWaveTheoryReg *WaveReg;

  TpRhopMode RhopMode;          ///<Inflow rhop mode (fixed or extrapolated).

  tdouble3 Direction;    ///<Inflow direction.
  tdouble3 PtPlane;      ///<Position to create Plane.
  tplane3f Plane;        ///<Inflow plane.
  unsigned NptInit;      ///<Number of inout points at the begining (first layer).
  unsigned NpartInit;    ///<Number of inout particles at the begining.

  tfloat3 *PtzPos;       ///<Positions in fluid domain to calculate zsurf [NptInit].
#ifdef _WITHGPU
  float3 *PtzPosg;       ///<Positions in fluid domain to calculate zsurf on GPU memory [NptInit].
  float  *PtzAuxg;       ///<Memory auxiliar to recude maximum zsurf on GPU memory [NptInit].
  float  *PtzAux;        ///<Memory auxiliar to recude maximum zsurf on CPU memory [NptInit].
#endif

  void ReadXml(const JXml *sxml,TiXmlElement* lis,const std::string &dirdatafile
    ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem);
  void CalculateVelMinMax(float &velmin,float &velmax)const;
  void LoadDomain();

public:
  JSphInOutZone(bool cpu,JLog2 *log,unsigned idzone,bool simulate2d,double simulate2dposy
    ,double dp,const tdouble3 &posmin,const tdouble3 &posmax,float gravityz
    ,const JXml *sxml,TiXmlElement* ele,const std::string &dirdatafile
    ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem);
  ~JSphInOutZone();
  void Reset();
  void GetConfig(std::vector<std::string> &lines)const;
  void CheckConfig()const;

  static TpVelProfile GetConfigVelProfile  (byte cfg){ return((TpVelProfile)(cfg&PVEL_MASK));  }
  static TpVelMode    GetConfigVelMode     (byte cfg){ return((TpVelMode)   (cfg&MVEL_MASK));  }
  static TpRhopMode   GetConfigRhopMode    (byte cfg){ return((TpRhopMode)  (cfg&MRHOP_MASK)); }
  static bool         GetConfigCheckInputDG(byte cfg){ return((cfg&CheckInput_MASK)!=0); }
  //static TpZsurfMode  GetConfigZsurfMode (byte cfg){ return((TpZsurfMode) ((cfg>>6)&3)); }
  byte GetConfigZone()const;
  byte GetConfigUpdate()const;

  unsigned LoadInletPoints(tdouble3 *pos);
  void LoadInletParticles(tdouble3 *pos);
  static float CalcVel(JSphInOutZone::TpVelProfile vprof,const tfloat4 &vdata,double posz);
  bool UpdateVelData(double timestep,bool full,tfloat4 &vel0,tfloat4 &vel1);
  bool UpdateZsurf(double timestep,bool full,float &zsurf);

  unsigned GetIdZone()const{ return(IdZone); }
  byte GetLayers()const{ return(Layers); }
  tdouble3 GetDirection()const{ return(Direction); }
  tplane3f GetPlane()const{ return(Plane); }

  const tdouble3* GetPtDomain()const{ return(PtDom); };
  tfloat3 GetBoxLimitMin()const{ return(BoxLimitMin); };
  tfloat3 GetBoxLimitMax()const{ return(BoxLimitMax); };
  inline bool InZoneBox(const tfloat3 &ps)const;
  bool InZone(bool useboxlimit,const tfloat3 &ps)const;

  float GetInputZbottom()const{ return(InputZbottom); }
  unsigned GetNptInit()const{ return(NptInit); }
  unsigned GetNpartInit()const{ return(NpartInit); }

  TpRhopMode   GetRhopMode()const{ return(RhopMode); }
  TpVelMode    GetVelMode()const{ return(VelMode); }
  TpVelProfile GetVelProfile()const{ return(VelProfile); }
  TpZsurfMode  GetZsurfMode()const{ return(ZsurfMode); }

  bool GetRefillAdvanced()const{ return(RefillingMode==TFI_Advanced); }
  bool GetExtrapolatedData()const{ return(RhopMode==MRHOP_Extrapolated || VelMode==MVEL_Extrapolated); }
  bool GetNoExtrapolatedData()const{ return(RhopMode!=MRHOP_Extrapolated || VelMode!=MVEL_Extrapolated); }
  bool GetVariableVel()const{ return(VelMode==MVEL_Variable); }
  bool GetInterpolatedVel()const{ return(VelMode==MVEL_Interpolated); }
  bool GetVariableZsurf()const{ return(ZsurfMode==ZSURF_Variable); }
  bool GetCalculatedZsurf()const{ return(ZsurfMode==ZSURF_Calculated); }
  bool GetRemoveZsurf()const{ return(RemoveZsurf); }

  bool GetSvVtkZsurf()const{ return(SvVtkZsurf); }

  double GetDistPtzPos()const{ return(PtzPos? Dp*2.: 0); }
  unsigned GetCountPtzPos()const{ return(PtzPos? NptInit: 0); }
  const tfloat3* GetPtzPos()const{ return(PtzPos); }
#ifdef _WITHGPU
  const float3* GetPtzPosg()const{ return(PtzPosg); }
  float* GetPtzAuxg()const{ return(PtzAuxg); }
  float* GetPtzAux()const{ return(PtzAux); }
#endif
  float GetInputZsurf()const{ return(InputZsurf); }
  void SetInputZsurf(float zsurf){ if(ZsurfMode==ZSURF_Calculated)InputZsurf=zsurf; }

  JSphInOutGridData* GetInputVelGrid()const{ return(InputVelGrid); }
  bool GetResetZVelGrid()const{ return(ResetZVelGrid); }

  void SaveVtkVelGrid(std::string filename)const;

  bool IsVariableInput()const;
  void SetExternalVarInput(bool active);
  bool GetExternalVarInput()const{ return(ExternalVarInput); };
  void InitExtVarInput();
  tfloat4 GetExtVarVelocity1()const;
  tfloat4 GetExtVarZsurf()const;
  void SetExtVarVelocity1(double t1,double v1,double t2,double v2);
  void SetExtVarZsurf(double t1,double z1,double t2,double z2);
};


#endif


