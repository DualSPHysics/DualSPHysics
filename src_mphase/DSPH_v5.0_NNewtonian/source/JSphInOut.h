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
#include "JSphInOutDef.h"
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
  #include "JSphTimersGpu.h"
#endif

#define DBG_INOUT_PTINIT 0   ///<JSphInOut: Saves VTK files (CfgInOut_PtInit.vtk and CfgInOut_PtInitZ.vtk) with initial inout points (0/1).
#define DBG_INOUT_PARTINIT 0 ///<JSphInOut: Saves VTK files (CfgInOut_InletIni_XXXX.vtk) with initial inout particles (0/1).

class JXml;
class TiXmlElement;
class JLog2;
class JSphInOutZone;
class JSphCpu;
class JSphGpuSingle;
class JDsPartsInit;
class JGaugeSystem;
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
  std::string DirDataFile;
  byte PeriActive;
  float CoefHydro;        ///<Constant to calculate hydrostatic rhop. CoefHydro=RhopZero*(-GravityZ)/CteB
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

  bool UseRefillAdvanced;   ///<Indicates if some inlet configuration uses advanced refilling method.
  bool UseZsurfNonUniform;  ///<Indicates if some inlet configuration uses Non-uniform zsurf. 
  bool UseAnalyticalData;   ///<Indicates if some inlet configuration uses analytical solution for density (Constant or Hydrostatic) or velocity (Fixed or Variable).
  bool UseExtrapolatedData; ///<Indicates if some inlet configuration uses extrapolated data for density or velocity.
  bool UseInterpolatedVel;  ///<Indicates if some inlet configuration uses interpolated velocity.

  bool VariableZsurf;         ///<Indicates if some inlet configuration uses variable zsurf.
  bool CalculatedZsurf;       ///<Indicates if some inlet configuration uses calculated zsurf.
  std::vector<JSphInOutZone*> List; ///<List of inlet/outlet configurations.
  unsigned ListSize;  ///<Number of inlet/outlet configurations.

  //-Data to manage all inlet/outlet conditions.
  tplane3f *Planes;    ///<Planes for inlet/outlet zones [ListSize].
  byte     *CfgZone;   ///<Information about VelMode, VelProfile, RhopMode and ConvertFluid. [ListSize].
  byte     *CfgUpdate; ///<Information about refilling mode, InputMode and RemoveZsurf. [ListSize].
  float    *Width;     ///<Zone width [ListSize].
  tfloat3  *DirData;   ///<Inflow direction [ListSize].
  tfloat3  *DirVel;    ///<Inflow direction for velocity (use FLT_MAX to ignore it). [ListSize].
  tfloat4  *VelData;   ///<Velocity coefficients for imposed velocity [ListSize*2].
  float    *Zsurf;     ///<Zsurf (it can be variable) [ListSize].
  
#ifdef _WITHGPU
  //-Data to manage all inlet/outlet conditions on GPU.
  float4 *Planesg;    ///<Planes for inlet/outlet zones [ListSize].
  float2 *BoxLimitg;  ///<BoxLimitMin/Max for inlet/outlet zones [ListSize*(xmin,xmax),ListSize*(ymin,ymax),ListSize*(zmin,zmax)]=[3*ListSize].
  byte   *CfgZoneg;   ///<Information about VelMode, VelProfile, RhopMode and ConvertFluid. [ListSize].
  byte   *CfgUpdateg; ///<Information about refilling mode, InputMode and RemoveZsurf. [ListSize].
  float  *Widthg;     ///<Zone width [ListSize].
  float3 *DirDatag;   ///<Inflow direction [ListSize].
  float3 *DirVelg;    ///<Inflow direction for velocity (use FLT_MAX to ignore it). [ListSize].
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

  void AllocateMemoryGpu(unsigned listsize);
  void FreeMemoryGpu();
#endif

public:
  const bool Cpu;
  const StCteSph CSP; ///<Structure with main SPH constants values and configurations.
  unsigned Nstep;     ///<Number of step in execution (for debug).

  std::vector<unsigned> MkFluidList;
    
  JSphInOut(bool cpu,const StCteSph &csp,std::string xmlfile,JXml *sxml
    ,std::string xmlpath,const std::string &dirdatafile);
  ~JSphInOut();
  void Reset();

  unsigned Config(double timestep,bool stable,byte periactive
    ,tdouble3 posmin,tdouble3 posmax,typecode codenewpart,const JDsPartsInit *partsdata
    ,JGaugeSystem *gaugesystem,JNumexLib *nuxlib);
    
  void LoadInitPartsData(unsigned idpfirst,unsigned npart,unsigned* idp,typecode* code,tdouble3* pos,tfloat4* velrhop);
  void InitCheckProximity(unsigned np,unsigned newnp,float scell,const tdouble3* pos,const unsigned *idp,typecode *code);


//-Specific code for CPU.
  unsigned CreateListSimpleCpu(unsigned npf,unsigned pini
    ,const typecode *code,int *inoutpart);
  unsigned CreateListCpu(unsigned npf,unsigned pini
    ,const tdouble3 *pos,const unsigned *idp,typecode *code,int *inoutpart);

  void SetAnalyticalDataCpu(float timestep,unsigned inoutcount,const int *inoutpart
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp,const float *zsurfpart
    ,tfloat4 *velrhop);

  void InterpolateVelCpu(float timestep,unsigned inoutcount,const int *inoutpart
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop);


  void CheckPartsIzone(std::string key,unsigned nstep,unsigned inoutcount,const int *inoutpart,typecode *code,unsigned *idp);
  unsigned ComputeStepCpu(unsigned inoutcount,int *inoutpart
    ,const JSphCpu *sphcpu,unsigned idnext,unsigned sizenp,unsigned np
    ,tdouble3 *pos,unsigned *dcell,typecode *code,unsigned *idp,const byte *zsurfok
    ,tfloat4 *velrhop,byte *newizone);
  unsigned ComputeStepFillingCpu(unsigned inoutcount,int *inoutpart
    ,const JSphCpu *sphcpu,unsigned idnext,unsigned sizenp,unsigned np
    ,tdouble3 *pos,unsigned *dcell,typecode *code,unsigned *idp,tfloat4 *velrhop
    ,const byte *zsurfok,float *prodist,tdouble3 *propos);

  void UpdateVelrhopM1Cpu(unsigned inoutcount,const int *inoutpart
    ,const tfloat4 *velrhop,tfloat4 *velrhopm1);


//-Specific code for GPU.
#ifdef _WITHGPU
  unsigned CreateListSimpleGpu(unsigned npf,unsigned pini
    ,const typecode *codeg,unsigned size,int *inoutpartg);
  unsigned CreateListGpu(unsigned npf,unsigned pini
    ,const double2 *posxyg,const double *poszg,typecode *codeg,unsigned size,int *inoutpartg);

  void SetAnalyticalDataGpu(float timestep,unsigned inoutcount,const int *inoutpartg
    ,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg
    ,const float *zsurfpart,float4 *velrhopg);

  void InterpolateVelGpu(float timestep,unsigned inoutcount,const int *inoutpartg
    ,const double2 *posxyg,const double *poszg,const typecode *codeg
    ,const unsigned *idpg,float4 *velrhopg);

  unsigned ComputeStepGpu(unsigned inoutcount,int *inoutpartg
    ,unsigned idnext,unsigned sizenp,unsigned np,double2 *posxyg,double *poszg
    ,unsigned *dcellg,typecode *codeg,unsigned *idpg,const byte *zsurfok
    ,float4 *velrhopg,byte *newizoneg,const JSphGpuSingle *gp);
  unsigned ComputeStepFillingGpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpartg
    ,unsigned idnext,unsigned sizenp,unsigned np
    ,double2 *posxyg,double *poszg,unsigned *dcellg,typecode *codeg,unsigned *idpg,float4 *velrhopg
    ,const byte* zsurfokg,float *prodistg,double2 *proposxyg,double *proposzg,TimersGpu timers);

  void UpdateVelrhopM1Gpu(unsigned inoutcount,const int *inoutpartg
    ,const float4 *velrhopg,float4 *velrhopm1g);
#endif

  void UpdateVelData(double timestep);
  void UpdateZsurfData(double timestep,bool full);

  void VisuConfig(std::string txhead,std::string txfoot)const;
  void SavePartFiles(unsigned part);
  void SaveVtkZsurf(unsigned part);

  unsigned GetCount()const{ return(unsigned(List.size())); }
  const JSphInOutZone* GetZone(unsigned ci)const{ return(ci<ListSize? List[ci]: NULL); }

  unsigned GetNpResizePlus0()const{ return(NpResizePlus0); }
  unsigned GetNpResizePlus1()const{ return(NpResizePlus1); }

  bool Use_RefillAdvanced()const{ return(UseRefillAdvanced); }
  bool Use_ZsurfNonUniform()const{ return(UseZsurfNonUniform); }
  bool Use_AnalyticalData()const{ return(UseAnalyticalData); }
  bool Use_ExtrapolatedData()const{ return(UseExtrapolatedData); }
  bool Use_InterpolatedVel()const{ return(UseInterpolatedVel); }

  bool GetVariableZsurf()const{ return(VariableZsurf); }
  bool GetCalculatedZsurf()const{ return(CalculatedZsurf); }

  void AddNewNp(unsigned newnp){ NewNpPart+=newnp; NewNpTotal+=newnp; }
  void ClearNewNpPart(){ NewNpPart=0; }
  unsigned GetNewNpPart()const{ return(NewNpPart); }
  ullong GetNewNpTotal()const{ return(NewNpTotal); }

  void SetCurrentNp(unsigned n){ CurrentNp=n; }
  unsigned GetCurrentNp()const{ return(CurrentNp); }

  bool GetReuseIds()const{ return(ReuseIds); }
  float GetDetermLimit()const{ return(DetermLimit); };
  byte GetExtrapolateMode()const{ return(ExtrapolateMode); };

  const tplane3f* GetPlanes ()const{ return(Planes);  };
  const byte*     GetCfgZone()const{ return(CfgZone); };
  const float*    GetWidth  ()const{ return(Width);   };
  const tfloat3*  GetDirData()const{ return(DirData); };
  const tfloat3*  GetDirVel ()const{ return(DirVel);  };
#ifdef _WITHGPU
  const float4*  GetPlanesg ()const{ return(Planesg);  };
  const byte*    GetCfgZoneg()const{ return(CfgZoneg); };
  const float*   GetWidthg  ()const{ return(Widthg);   };
  const float3*  GetDirDatag()const{ return(DirDatag); };
  const float3*  GetDirVelg ()const{ return(DirVelg);  };
  byte  GetExtrapRhopMask()const;
  byte  GetExtrapVelMask() const;
#endif

};



#endif


