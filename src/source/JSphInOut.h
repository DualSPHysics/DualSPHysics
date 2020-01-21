//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2016, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:#############################################################################

/// \file JSphInOut.h \brief Declares the class \ref JSphInOut.

#ifndef _JSphInOut_
#define _JSphInOut_

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
class JSphCpu;
class JWaveTheoryReg;
class JSphMk;
class JSphPartsInit;
class JSphGpuSingle;

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

  bool RevNewnpPerSec;   ///<Indicate when NewnpPerSec should be calculate.
  double NewnpPerSec;    ///<Number of new fluid particles per second.

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

  void ReadXml(JXml *sxml,TiXmlElement* lis,const std::string &dirdatafile,const JSphPartsInit *partsdata);
  void CalculateVelMinMax(float &velmin,float &velmax)const;
  void LoadDomain();

public:
  JSphInOutZone(bool cpu,JLog2 *log,unsigned idzone,bool simulate2d,double simulate2dposy
    ,double dp,const tdouble3 &posmin,const tdouble3 &posmax,JXml *sxml,TiXmlElement* ele
    ,const std::string &dirdatafile,const JSphPartsInit *partsdata);
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

  unsigned CalcResizeNp(double timestep,double timeinterval);
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
  inline bool InZone(bool useboxlimit,const tfloat3 &ps)const;

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

//##############################################################################
//# JSphInOut
//##############################################################################
/// \brief Manages the info of inlet conditions.
class JSphInOut : protected JObject
{
private:
  JLog2 *Log;
  const std::string XmlFile;
  const std::string XmlPath;
  static const unsigned MaxZones=CODE_TYPE_FLUID_INOUTNUM;  ///<Maximum number of inout configurations.

  //-Basic simulation parameters.
  bool Stable;
  bool Simulate2D;        ///<Indicates 2D simulation.
  double Simulate2DPosY;  ///<Y value in 2D simulations.
  std::string DirDataFile;

  byte PeriActive;
  float RhopZero;
  float CteB;
  float Gamma;
  float GravityZ;
  float CoefHydro;        ///<Constant to calculate hydrostatic rhop. CoefHydro=RhopZero*(-GravityZ)/CteB
  double Dp;              ///<Distance between particles.
  tdouble3 MapRealPosMin;
  tdouble3 MapRealPosMax;
  typecode CodeNewPart;   ///<Code for new fluid particles created starting from inlet particles.

  bool ReuseIds;         ///<Id of particles excluded values ​​are reused.
  double ResizeTime;     ///<Time to calculate number of new particles.
  float DetermLimit;     ///<Limit for determinant. Use 1e-3 for first_order or 1e+3 for zeroth_order (default=1e+3).
  byte ExtrapolateMode;  ///<Calculation mode for rhop and velocity extrapolation from ghost nodes 1:fast-single, 2:single, 3:double (default=1).

  bool UseBoxLimit;     ///<In/out process is only applied to InOut zones delimited by BoxLimit (default=true).
  tfloat3 FreeCentre;   ///<Centre of zone where InOut is not applied (default=centre of simulation domain).
  tfloat3 FreeLimitMin; ///<Minimum limit where InOut is not applied.
  tfloat3 FreeLimitMax; ///<Maximum limit where InOut is not applied.

  ullong NewNpTotal;   ///<Total number of fluid particles created by inlet.
  unsigned NewNpPart;  ///<Number of fluid particles created by inlet since last PART.
  unsigned CurrentNp;  ///<Current number of inout particles (including periodic particles). 

  bool RefillAdvanced;       ///<Indicates if some inlet configuration uses advanced refilling method.
  bool VariableVel;          ///<Indicates if some inlet configuration uses variable velocity.
  bool VariableZsurf;        ///<Indicates if some inlet configuration uses variable zsurf.
  bool CalculatedZsurf;      ///<Indicates if some inlet configuration uses calculated zsurf.
  bool ExtrapolatedData;     ///<Indicates if some inlet configuration uses extrapolated data.
  bool NoExtrapolatedData;   ///<Indicates if any inlet configuration uses extrapolated data.
  bool InterpolatedVel;      ///<Indicates if some inlet configuration uses interpolated velocity.
  std::vector<JSphInOutZone*> List; ///<List of inlet/outlet configurations.
  unsigned ListSize;  ///<Number of inlet/outlet configurations.

  //-Data to manage all inlet/outlet conditions.
  tplane3f *Planes;    ///<Planes for inlet/outlet zones [ListSize].
  byte     *CfgZone;   ///<Information about VelMode, VelProfile, RhopMode and ConvertFluid. [ListSize].
  byte     *CfgUpdate; ///<Information about refilling mode, InputMode and RemoveZsurf. [ListSize].
  float    *Width;     ///<Zone width [ListSize].
  tfloat3  *DirData;   ///<Inflow direction [ListSize].
  tfloat4  *VelData;   ///<Velocity coefficients for imposed velocity [ListSize*2].
  float    *Zbottom;   ///<Zbottom for hydrostatic rhop (it is constant) [ListSize].
  float    *Zsurf;     ///<Zsurf (it can be variable) [ListSize].
  
#ifdef _WITHGPU
  //-Data to manage all inlet/outlet conditions on GPU.
  float4 *Planesg;    ///<Planes for inlet/outlet zones [ListSize].
  float2 *BoxLimitg;  ///<BoxLimitMin/Max for inlet/outlet zones [ListSize*(xmin,xmax),ListSize*(ymin,ymax),ListSize*(zmin,zmax)]=[3*ListSize].
  byte   *CfgZoneg;   ///<Information about VelMode, VelProfile, RhopMode and ConvertFluid. [ListSize].
  byte   *CfgUpdateg; ///<Information about refilling mode, InputMode and RemoveZsurf. [ListSize].
  float  *Widthg;     ///<Zone width [ListSize].
  float3 *DirDatag;   ///<Inflow direction [ListSize].
  float  *Zsurfg;     ///<Zsurf (it can be variable) [ListSize].
#endif

  //-Data to refill inlet/outlet zone.
  unsigned PtCount;  ///<Number of points.
  byte     *PtZone;  ///<Zone for each point [PtCount]. Data is constant after configuration.
  tdouble3 *PtPos;   ///<Position of points [PtCount]. Data is constant after configuration.
  float    *PtAuxDist;

#ifdef _WITHGPU
  //-Data to refill inlet/outlet zone on GPU.
  byte    *PtZoneg;   ///<Zone for each point [PtCount]. Data is constant after configuration.
  double2 *PtPosxyg;  ///<Position of points [PtCount]. Data is constant after configuration.
  double  *PtPoszg;   ///<Position of points [PtCount]. Data is constant after configuration.
  float   *PtAuxDistg;
#endif

  void LoadXmlInit(JXml *sxml,const std::string &place);
  void LoadFileXml(const std::string &file,const std::string &path,const JSphPartsInit *partsdata);
  void LoadXml(JXml *sxml,const std::string &place,const JSphPartsInit *partsdata);
  void ReadXml(JXml *sxml,TiXmlElement* ele,const JSphPartsInit *partsdata);

  void ComputeFreeDomain();
  void SaveVtkDomains();
  void SaveVtkVelGrid();

  void AllocateMemory(unsigned listsize);
  void FreeMemory();

  void AllocatePtMemory(unsigned ptcount);
  void FreePtMemory();
#ifdef _WITHGPU
  void AllocatePtMemoryGpu(unsigned ptcount);
  void FreePtMemoryGpu();
#endif

  bool UpdateVelData(double timestep,bool full);
  bool UpdateZsurf(double timestep,bool full);

#ifdef _WITHGPU
  void AllocateMemoryGpu(unsigned listsize);
  void FreeMemoryGpu();
#endif

public:
  const bool Cpu;

  std::vector<unsigned> MkFluidList;
    
  JSphInOut(bool cpu,JLog2 *log,std::string xmlfile,JXml *sxml,std::string xmlpath,const std::string &dirdatafile);
  ~JSphInOut();
  void Reset();

  unsigned Config(double timestep,bool stable,bool simulate2d,double simulate2dposy
    ,byte periactive,float rhopzero,float cteb,float gamma,tfloat3 gravity,double dp
    ,tdouble3 posmin,tdouble3 posmax,typecode codenewpart,const JSphPartsInit *partsdata);
    
  void LoadInitPartsData(unsigned idpfirst,unsigned npart,unsigned* idp,typecode* code,tdouble3* pos,tfloat4* velrhop);
  void InitCheckProximity(unsigned np,unsigned newnp,float scell,const tdouble3* pos,const unsigned *idp,typecode *code);


//-Specific code for CPU.
  unsigned CreateListSimpleCpu(unsigned nstep,unsigned npf,unsigned pini
    ,const typecode *code,int *inoutpart);
  unsigned CreateListCpu(unsigned nstep,unsigned npf,unsigned pini
    ,const tdouble3 *pos,const unsigned *idp,typecode *code,int *inoutpart);
  void UpdateDataCpu(float timestep,bool full,unsigned inoutcount,const int *inoutpart
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop);

  void InterpolateVelCpu(float timestep,unsigned inoutcount,const int *inoutpart
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop);
  void InterpolateResetZVelCpu(unsigned inoutcount,const int *inoutpart
    ,const typecode *code,tfloat4 *velrhop);


  void CheckPartsIzone(std::string key,unsigned nstep,unsigned inoutcount,const int *inoutpart,typecode *code,unsigned *idp);
  unsigned ComputeStepCpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpart
    ,const JSphCpu *sphcpu,unsigned idnext,unsigned sizenp,unsigned np
    ,tdouble3 *pos,unsigned *dcell,typecode *code,unsigned *idp,tfloat4 *velrhop,byte *newizone);
  unsigned ComputeStepFillingCpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpart
    ,const JSphCpu *sphcpu,unsigned idnext,unsigned sizenp,unsigned np
    ,tdouble3 *pos,unsigned *dcell,typecode *code,unsigned *idp,tfloat4 *velrhop
    ,float *prodist,tdouble3 *propos);

  void UpdateVelrhopM1Cpu(unsigned inoutcount,const int *inoutpart
    ,const tfloat4 *velrhop,tfloat4 *velrhopm1);

  void ClearInteractionVarsCpu(unsigned npf,unsigned pini,const typecode *code
    ,tfloat3 *ace,float *ar,tfloat4 *shiftposfs);


//-Specific code for GPU.
#ifdef _WITHGPU
  unsigned CreateListSimpleGpu(unsigned nstep,unsigned npf,unsigned pini
    ,const typecode *codeg,unsigned size,int *inoutpartg);
  unsigned CreateListGpu(unsigned nstep,unsigned npf,unsigned pini
    ,const double2 *posxyg,const double *poszg,typecode *codeg,unsigned size,int *inoutpartg);
  void UpdateDataGpu(float timestep,bool full,unsigned inoutcount,const int *inoutpartg
    ,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg,float4 *velrhopg);

  void InterpolateVelGpu(float timestep,unsigned inoutcount,const int *inoutpartg
    ,double2 *posxyg,double *poszg,typecode *codeg,unsigned *idpg,float4 *velrhopg);
  void InterpolateResetZVelGpu(unsigned inoutcount,const int *inoutpartg
    ,typecode *codeg,float4 *velrhopg);

  unsigned ComputeStepGpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpartg
    ,unsigned idnext,unsigned sizenp,unsigned np,double2 *posxyg,double *poszg
    ,unsigned *dcellg,typecode *codeg,unsigned *idpg,float4 *velrhopg,byte *newizoneg,const JSphGpuSingle *gp);
  unsigned ComputeStepFillingGpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpartg
    ,unsigned idnext,unsigned sizenp,unsigned np
    ,double2 *posxyg,double *poszg,unsigned *dcellg,typecode *codeg,unsigned *idpg,float4 *velrhopg
    ,float *prodistg,double2 *proposxyg,double *proposzg,TimersGpu timers);

  void UpdateVelrhopM1Gpu(unsigned inoutcount,const int *inoutpartg
    ,const float4 *velrhopg,float4 *velrhopm1g);

  void ClearInteractionVarsGpu(unsigned npf,unsigned pini,const typecode *codeg
    ,float3 *aceg,float *arg,float *viscdtg,float4 *shiftposfsg);
#endif


  void UpdateZsurf(double timestep){ UpdateZsurf(timestep,false); }

  void VisuConfig(std::string txhead,std::string txfoot)const;
  void SaveVtkZsurf(unsigned part);

  unsigned GetCount()const{ return(unsigned(List.size())); }
  const JSphInOutZone* GetZone(unsigned ci)const{ return(ci<ListSize? List[ci]: NULL); }

  bool GetNoExtrapolatedData()const{ return(NoExtrapolatedData); }
  bool GetExtrapolatedData()const{ return(ExtrapolatedData); }
  bool GetInterpolatedVel()const{ return(InterpolatedVel); }
  bool GetVariableZsurf()const{ return(VariableZsurf); }
  bool GetCalculatedZsurf()const{ return(CalculatedZsurf); }
  bool GetRefillAdvanced()const{ return(RefillAdvanced); }

  double GetDistPtzPos(unsigned ci)const{     return(ci<ListSize? List[ci]->GetDistPtzPos()  : 0); }
  unsigned GetCountPtzPos(unsigned ci)const{  return(ci<ListSize? List[ci]->GetCountPtzPos() : 0); }
  const tfloat3* GetPtzPos(unsigned ci)const{ return(ci<ListSize? List[ci]->GetPtzPos()      : NULL); }
#ifdef _WITHGPU
  const float3* GetPtzPosg(unsigned ci)const{     return(ci<ListSize? List[ci]->GetPtzPosg()    : NULL); }
  float* GetPtzAuxg(unsigned ci)const{ return(ci<ListSize? List[ci]->GetPtzAuxg(): NULL); }
  float* GetPtzAux (unsigned ci)const{ return(ci<ListSize? List[ci]->GetPtzAux() : NULL); }
#endif
  void SetInputZsurf(unsigned ci,float zsurf){ if(ci<ListSize)List[ci]->SetInputZsurf(zsurf); }
  float GetZbottom(unsigned ci)const{ return(ci<ListSize? List[ci]->GetInputZbottom(): 0); }

  unsigned CalcResizeNp(double timestep)const;

  void AddNewNp(unsigned newnp){ NewNpPart+=newnp; NewNpTotal+=newnp; }
  void ClearNewNpPart(){ NewNpPart=0; }
  unsigned GetNewNpPart()const{ return(NewNpPart); }
  ullong GetNewNpTotal()const{ return(NewNpTotal); }

  void SetCurrentNp(unsigned n){ CurrentNp=n; }
  unsigned GetCurrentNp()const{ return(CurrentNp); }

  bool GetReuseIds()const{ return(ReuseIds); }
  float GetDetermLimit()const{ return(DetermLimit); };
  byte GetExtrapolateMode()const{ return(ExtrapolateMode); };

  const tplane3f* GetPlanes() const{ return(Planes);  };
  const byte*     GetCfgZone()const{ return(CfgZone); };
  const float*    GetWidth()  const{ return(Width);   };
  const tfloat3*  GetDirData()const{ return(DirData); };
#ifdef _WITHGPU
  const float4*  GetPlanesg() const{ return(Planesg);  };
  const byte*    GetCfgZoneg()const{ return(CfgZoneg); };
  const float*   GetWidthg()  const{ return(Widthg);   };
  const float3*  GetDirDatag()const{ return(DirDatag); };
  byte  GetExtrapRhopMask()const{ return(byte(JSphInOutZone::MRHOP_Extrapolated)); };
  byte  GetExtrapVelMask() const{ return(byte(JSphInOutZone::MVEL_Extrapolated)); };
#endif

};



#endif


