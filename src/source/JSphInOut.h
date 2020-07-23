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
//:# - Objeto JXml pasado como const para operaciones de lectura. (18-03-2020)  
//:# - Uso de JNumexLib para evaluar expresiones del XML. (18-03-2020)  
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (18-03-2020)  
//:# - Memory resizing configuration was improved. (02-07-2020)
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

#define DBG_INOUT_PTINIT 0   ///<JSphInOut: Saves VTK files (CfgInOut_PtInit.vtk and CfgInOut_PtInitZ.vtk) with initial inout points (0/1).
#define DBG_INOUT_PARTINIT 0 ///<JSphInOut: Saves VTK files (CfgInOut_InletIni_XXXX.vtk) with initial inout particles (0/1).

class JXml;
class TiXmlElement;
class JLog2;
//class JLinearValue;
class JSphInOutZone;
//class JSphInOutPoints;
//class JSphInOutGridData;
class JSphCpu;
class JSphGpuSingle;
//class JWaveTheoryReg;
//class JSphMk;
class JDsPartsInit;
class JGaugeSystem;
//class JGaugeSwl;
class JNumexLib;

//##############################################################################
//# XML format in _FmtXML_InOut.xml.
//##############################################################################

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

  bool ReuseIds;          ///<Id of particles excluded values ​​are reused.
  float MemoryResize0;    ///<Multipier to calculate the minimum initial resize (default=2).
  float MemoryResize1;    ///<Multipier to calculate the following resizes (default=4).
  unsigned NpResizePlus0; ///<Extra number of particles for initial resize memory.
  unsigned NpResizePlus1; ///<Extra number of particles for resize memory (0 does not allow resizing).

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

  void LoadXmlInit(const JXml *sxml,const std::string &place);
  void LoadFileXml(const std::string &file,const std::string &path
    ,JNumexLib *nuxlib,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem);
  void LoadXml(const JXml *sxml,const std::string &place
    ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem);
  void ReadXml(const JXml *sxml,TiXmlElement* ele
    ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem);

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
    ,tdouble3 posmin,tdouble3 posmax,typecode codenewpart,const JDsPartsInit *partsdata
    ,JGaugeSystem *gaugesystem,JNumexLib *nuxlib);
    
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
#endif


  void UpdateZsurf(double timestep){ UpdateZsurf(timestep,false); }

  void VisuConfig(std::string txhead,std::string txfoot)const;
  void SavePartFiles(unsigned part);
  void SaveVtkZsurf(unsigned part);

  unsigned GetCount()const{ return(unsigned(List.size())); }
  const JSphInOutZone* GetZone(unsigned ci)const{ return(ci<ListSize? List[ci]: NULL); }

  unsigned GetNpResizePlus0()const{ return(NpResizePlus0); }
  unsigned GetNpResizePlus1()const{ return(NpResizePlus1); }

  bool GetNoExtrapolatedData()const{ return(NoExtrapolatedData); }
  bool GetExtrapolatedData()const{ return(ExtrapolatedData); }
  bool GetInterpolatedVel()const{ return(InterpolatedVel); }
  bool GetVariableZsurf()const{ return(VariableZsurf); }
  bool GetCalculatedZsurf()const{ return(CalculatedZsurf); }
  bool GetRefillAdvanced()const{ return(RefillAdvanced); }

  double GetDistPtzPos(unsigned ci)const;
  unsigned GetCountPtzPos(unsigned ci)const;
  const tfloat3* GetPtzPos(unsigned ci)const;
#ifdef _WITHGPU
  const float3* GetPtzPosg(unsigned ci)const;
  float* GetPtzAuxg(unsigned ci)const;
  float* GetPtzAux (unsigned ci)const;
#endif
  void SetInputZsurf(unsigned ci,float zsurf);
  float GetZbottom(unsigned ci)const;

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
  byte  GetExtrapRhopMask()const;
  byte  GetExtrapVelMask() const;
#endif

};



#endif


