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
//:# - Clase para gestionar zonas de shifting. (24-11-2019)
//:# - Improved exception managment. (19-03-2020)  
//:# - Objeto JXml pasado como const para operaciones de lectura. (19-03-2020)  
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (19-03-2020)  
//:# - Nuevo metodo GetConfigInfo(). (03-06-2020)  
//:# - Cambio de nombre de J.Shifting a J.SphShifting. (28-06-2020)
//:#############################################################################

/// \file JSphShifting.h \brief Declares the class \ref JSphShifting.

#ifndef _JSphShifting_
#define _JSphShifting_

#include "JObject.h"
#include "DualSphDef.h"
#include "JMatrix4.h"
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

#include <string>
#include <vector>

class JLog2;
class JXml;
class TiXmlElement;

//##############################################################################
//# XML format in _FmtXML_Shifting.xml.
//##############################################################################

//##############################################################################
//# JSphShifting
//##############################################################################
/// \brief Manages the info of damping zones.

class JSphShiftingZone : protected JObject
{
private:
  tdouble3 PosRef;    ///<Reference position of shifting domain.
  tdouble3 Vecx;      ///<Vector to define shifting domain.
  tdouble3 Vecy;      ///<Vector to define shifting domain.
  tdouble3 Vecz;      ///<Vector to define shifting domain.

  bool UsePosMax;     ///<Use parallel box domain definition.
  tdouble3 PosMax;    ///<Maximum position to define parallel box domain.

  tplane3d DomPlax;   ///<Plane x for plane-domain definition.
  tplane3d DomPlay;   ///<Plane y for plane-domain definition.
  tplane3d DomPlaz;   ///<Plane z for plane-domain definition.
  tdouble3 DomPladis; ///<Distance coefficients for plane-domain definition.

private:
  void PrepareZone();

public:
  const unsigned Id;
  JSphShiftingZone(unsigned id,const tdouble3 &posref,const tdouble3 &vx,const tdouble3 &vy,const tdouble3 &vz);
  ~JSphShiftingZone();
  void Reset();

  tdouble3 GetVecx()const{ return(Vecx); }
  tdouble3 GetVecy()const{ return(Vecy); }
  tdouble3 GetVecz()const{ return(Vecz); }

  bool GetUsePosMax()const{ return(UsePosMax); }
  
  tdouble3 GetPosMin()const{ return(PosRef); }
  tdouble3 GetPosMax()const{ return(PosMax); }

  tplane3d GetDomPlax()const{ return(DomPlax); }
  tplane3d GetDomPlay()const{ return(DomPlay); }
  tplane3d GetDomPlaz()const{ return(DomPlaz); }
  tdouble3 GetDomPladis()const{ return(DomPladis); }
};

//##############################################################################
//# JSphShifting
//##############################################################################
/// \brief Manages the info of damping zones.

class JSphShifting : protected JObject
{
private:
  JLog2* Log;
  const bool Simulate2D;
  const double Dp;       ///<Initial distance between particles [m].
  float KernelH;         ///<The smoothing length of SPH kernel [m].

  TpShifting ShiftMode;  ///<Type of Shifting: None, NoBound, NoFixed, Full.
  float ShiftCoef;       ///<Coefficient for shifting computation.
  float ShiftTFS;        ///<Threshold to detect free surface. Typically 1.5 for 2D and 2.75 for 3D (def=0).

  unsigned ZonesXml;     ///<Number of shifting zones defined by XML.
  unsigned ZonesPosmax;  ///<Number of zones defined by position min-max.
  std::vector<JSphShiftingZone*> Zones;

  JMatrix4d ReadRotate3D(const JXml *sxml,TiXmlElement* ele);
  void ReadXml(const JXml *sxml,TiXmlElement* ele);

  template<bool first,bool dbl> void InitCpuPosMax(unsigned n,unsigned pini
    ,const tdouble3& pmin1,const tdouble3& pmax1,const tdouble3& pmin2,const tdouble3& pmax2
    ,const tdouble3* pos,tfloat4* shiftposfs)const;
  template<bool first,bool dbl> void InitCpuPlanes(unsigned n,unsigned pini
    ,const tplane3d& plax1,const tplane3d& play1,const tplane3d& plaz1,const tdouble3& pladis1
    ,const tplane3d& plax2,const tplane3d& play2,const tplane3d& plaz2,const tdouble3& pladis2
    ,const tdouble3* pos,tfloat4* shiftposfs)const;


public:
  JSphShifting(bool simulate2d,double dp,float kernelh);
  ~JSphShifting();
  void Reset();

  void ConfigBasic(TpShifting shiftmode,float shiftcoef=-2.f,float shifttfs=0);
  void LoadXml(const JXml *sxml,const std::string &place);
  void AddZone(bool fromxml,const tdouble3 &posref,const tdouble3 &vx,const tdouble3 &vy,const tdouble3 &vz);

  void VisuConfig(std::string txhead="",std::string txfoot="");
  std::string GetConfigInfo()const;
  void SaveVtkConfig()const;

  unsigned GetCount()const{ return(unsigned(Zones.size())); }

  std::string GetShiftingModeStr()const;
  TpShifting  GetShiftMode()const{ return(ShiftMode); }
  float       GetShiftCoef()const{ return(ShiftCoef); }
  float       GetShiftTFS ()const{ return(ShiftTFS); }

  void InitCpu(unsigned n,unsigned pini,const tdouble3* pos,tfloat4* shiftposfs)const;
  void RunCpu(unsigned n,unsigned pini,double dt,const tfloat4* velrhop,tfloat4* shiftposfs)const;

#ifdef _WITHGPU
  void InitGpu(unsigned n,unsigned pini,const double2* posxy,const double* posz,float4* shiftposfs,cudaStream_t stm=NULL)const;
  void RunGpu(unsigned n,unsigned pini,double dt,const float4* velrhop,float4* shiftposfs,cudaStream_t stm=NULL)const;
#endif

};


#endif


