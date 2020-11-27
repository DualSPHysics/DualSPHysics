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
//:# - BoxLimit en Y se amplia en 2D para evitar problemas en comparaciones 
//:#   float/double. (18-07-2020)
//:#############################################################################

/// \file JSphInOutZone.h \brief Declares the class \ref JSphInOutZone.

#ifndef _JSphInOutZone_
#define _JSphInOutZone_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#include "JSphInOutDef.h"

#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
  #include "JSphTimersGpu.h"
#endif

class JXml;
class TiXmlElement;
class JLog2;
class JLinearValue;
class JSphInOutPoints;
class JSphInOutVel;
class JSphInOutZsurf;
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
  
  static const unsigned CheckInput_MASK  =0x80;  ///<Mask to obtain value (1000 0000).

  static const unsigned RefillAdvanced_MASK=0x01;  ///<Mask to obtain value (0000 0001).
  static const unsigned RefillSpFull_MASK  =0x02;  ///<Mask to obtain value (0000 0010).
  static const unsigned RemoveInput_MASK   =0x04;  ///<Mask to obtain value (0000 0100).
  static const unsigned RemoveZsurf_MASK   =0x08;  ///<Mask to obtain value (0000 1000).
  static const unsigned ConvertInput_MASK  =0x10;  ///<Mask to obtain value (0001 0000).

private:
  const bool Cpu;
  const StCteSph CSP; ///<Structure with main SPH constants values and configurations.
  JLog2 *Log;
  const unsigned IdZone;
  const tdouble3 MapRealPosMin;
  const tdouble3 MapRealPosMax;

  //-Configuration parameters.
  JSphInOutPoints *Points;      ///<Definition of inlet points.
  byte Layers;                  ///<Number of inlet particle layers.
  TpInInput InputMode;          ///<Defines the treatment of fluid particles entering a inlet/outlet zone (default=0).
  bool InputCheck;              ///<Checks fluid input for RemoveZsurf=true or InputMode!=InInput_Free.
  TpInRefilling RefillingMode;  ///<Refilling mode. 0:Simple full, 1:Simple below zsurf, 2:Advanced for reverse flows (very slow) (default=1).

  //-Domain data from Points object. The PtDom[8] is the reference point in inout plane.
  static const unsigned PTDOMSIZE=10;
  tdouble3 PtDom[PTDOMSIZE];
  tfloat3 BoxLimitMin;
  tfloat3 BoxLimitMax;

  tdouble3 Direction;      ///<Inflow direction.
  tdouble3 PtPlane;        ///<Position to create inlet plane.
  tplane3f Plane;          ///<Inlet plane.
  unsigned NptInit;        ///<Number of inout points at the begining (first layer).
  unsigned NpartInit;      ///<Number of inout particles at the begining.

  //-Configuration of velocity, zsurf and density.
  TpInVelMode VelMode;         ///<Inflow velocity mode (fixed or variable).
  JSphInOutVel *InOutVel;      ///<Manages Velocity configuration.

  TpInZsurfMode ZsurfMode;     ///<Inflow zsurf mode (fixed, variable...).
  JSphInOutZsurf *InOutZsurf;  ///<Manages Zsurf configuration.

  TpInRhopMode RhopMode;       ///<Inflow rhop mode (fixed or extrapolated).
  
  bool ExternalVarInput;       ///<Receives velocity and zsurf values from external sources. (default=false)

private:
  void ReadXml(const JXml *sxml,TiXmlElement* lis,const std::string &dirdatafile
    ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem);
  void LoadDomain();

public:
  JSphInOutZone(bool cpu,unsigned idzone,const StCteSph &csp
    ,const tdouble3 &posmin,const tdouble3 &posmax
    ,const JXml *sxml,TiXmlElement* ele,const std::string &dirdatafile
    ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem);
  ~JSphInOutZone();
  void Reset();
  void GetConfig(std::vector<std::string> &lines)const;
  void CheckConfig()const;

  static TpInVelProfile GetConfigVelProfile  (byte cfg){ return((TpInVelProfile)(cfg&InVelP_MASK));  }
  static TpInVelMode    GetConfigVelMode     (byte cfg){ return((TpInVelMode)   (cfg&InVelM_MASK));  }
  static TpInRhopMode   GetConfigRhopMode    (byte cfg){ return((TpInRhopMode)  (cfg&InRhop_MASK)); }
  static bool           GetConfigCheckInputDG(byte cfg){ return((cfg&CheckInput_MASK)!=0); }
  byte GetConfigZone()const;
  byte GetConfigUpdate()const;

  unsigned LoadInletPoints(tdouble3 *pos);
  void LoadInitialParticles(unsigned npartinit,tdouble3 *pos);
  static float CalcVel(TpInVelProfile vprof,const tfloat4 &vdata,double posz);

  unsigned GetIdZone()const{ return(IdZone); }
  byte GetLayers()const{ return(Layers); }
  tdouble3 GetDirection()const{ return(Direction); }
  tplane3f GetPlane()const{ return(Plane); }

  const tdouble3* GetPtDomain()const{ return(PtDom); };
  tfloat3 GetBoxLimitMin()const{ return(BoxLimitMin); };
  tfloat3 GetBoxLimitMax()const{ return(BoxLimitMax); };
  inline bool InZoneBox(const tfloat3 &ps)const;
  bool InZone(bool useboxlimit,const tfloat3 &ps)const;

  unsigned GetNptInit()const{ return(NptInit); }
  unsigned GetNpartInit()const{ return(NpartInit); }

  TpInRhopMode   GetRhopMode()const{ return(RhopMode); }
  TpInVelMode    GetVelMode()const{ return(VelMode); }
  TpInZsurfMode  GetZsurfMode()const{ return(ZsurfMode); }

  bool Use_AnalyticalData()const{ 
    return(RhopMode==InRhop_Constant || RhopMode==InRhop_Hydrostatic
    || VelMode==InVelM_Fixed || VelMode==InVelM_Variable); 
  }
  bool Use_ExtrapolatedData()const{ 
    return(RhopMode==InRhop_Extrapolated || VelMode==InVelM_Extrapolated); 
  }
  bool Use_InterpolatedVel()const{ return(VelMode==InVelM_Interpolated); }
  bool Use_VelFixed()const{ return(VelMode==InVelM_Fixed); }
  bool Use_VelVariable()const{ return(VelMode==InVelM_Variable); }
  bool Use_RefillAdvanced()const{ return(RefillingMode==InRefill_Advanced); }

  bool GetVariableZsurf()const{ return(ZsurfMode==InZsurf_Variable); }
  bool GetCalculatedZsurf()const{ return(ZsurfMode==InZsurf_Calculated); }


  JSphInOutVel*   GetInOutVel  ()const{ return(InOutVel); }  
  JSphInOutZsurf* GetInOutZsurf()const{ return(InOutZsurf); }  

};


#endif


