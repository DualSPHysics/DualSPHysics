//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphInOutZsurf.h \brief Declares the class \ref JSphInOutZsurf.

#ifndef _JSphInOutZsurf_
#define _JSphInOutZsurf_

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
class JGaugeMesh; //<vs_meeshdat>

namespace jmsh{       //<vs_meeshdat>
  class JMeshTDatas;  //<vs_meeshdat>
}                     //<vs_meeshdat>

//##############################################################################
//# JSphInOutZsurf
//##############################################################################
/// \brief Manages the Zsurf data for inlet/outlet conditions.

class JSphInOutZsurf : protected JObject
{
protected:
  JLog2* Log;
  const bool Cpu;
  const StCteSph CSP; ///<Structure with main SPH constants values and configurations.
  const unsigned IdZone;
  const tdouble3 Direction;   ///<Inflow direction.
  const tdouble3 ZonePosMin;  ///<Minimum position of inlet points.
  const tdouble3 ZonePosMax;  ///<Maximum position of inlet points.
  bool OldCode;

  TpInZsurfMode ZsurfMode; ///<Inflow zsurf mode (fixed, variable...).
  bool UniformZsurf;       ///<Zsurf is the same for all inlet points.

  float InputZbottom;      ///<Bottom level of water (it is constant).
  float InputZsurf;        ///<Input uniform zsurf. It is FLT_MAX when ZsurfMode==InZsurf_Undefined or UniformZsurf==false.
  bool RemoveZsurf;        ///<Removes particles above the Zsurf limit (default=false).
  bool SvVtkZsurf;         ///<Creates VTK files with Zsurf for each PART (default=false).
  float Zsurf0;            ///<Initial Z-surface value for ZsurfMode==InZsurf_Calculated.
  double ZsurfMin;         ///<Minimum Z-surface to reduce the calculation area (default=minimum Z inlet).
  float ZsurfFit;          ///<Set calculated zsurf (default -0.5*dp).
  JLinearValue* InputTime; ///<Input zsurf in time (for ZsurfMode==InZsurf_Variable).

  //<vs_meeshdat_ini>
  //-Zsurf variable from mesh-data file.
  std::string MeshFile;    ///<File with mesh-data.
  tdouble3 MeshSetPos;     ///<Position offset (applied after filters).
  double MeshInitialTime;  ///<Defines initial time of data.
  double MeshLoopTbegRq;   ///<Requested begin time of loop (DBL_MAX=disabled).
  double MeshLoopTmaxRq;   ///<Requested final time of loop (DBL_MAX=disabled).
  double MeshLoopTbeg;     ///<Begin time of loop (DBL_MAX=disabled).
  double MeshLoopTmax;     ///<Final time of loop (DBL_MAX=disabled).

  jmsh::JMeshTDatas* InputMesh; ///<Input zsurf in time from mesh data (for ZsurfMode==InZsurf_Variable).
  //<vs_meeshdat_end>

  StZsurfResult ZsurfResults; ///<Stores zsurf results to send. 

  //-Zsurf line definition.
  tdouble3 GaugePt0;  ///<Reference position to measure Zsurf (at CSP.kernelsize distance).
  tdouble3 GaugePtz;  ///<Maximum Z position to measure Zsurf (at CSP.kernelsize distance).
  tdouble3 GaugePtx;  ///<Maximum horizontal position to measure Zsurf (at CSP.kernelsize distance).
  double Dptz;        ///<Distance between gauge points in Z.
  double Dptx;        ///<Distance between gauge points in horizontal direction.
  unsigned Nptz;      ///<Number of points in Z.
  unsigned Nptx;      ///<Number of points in horizontal direction.
  tplane3d PlaDisx;   ///<Reference plane to compute horizontal distance.

  JGaugeSwl*  GaugeSwl;  ///<Gauge object to compute zsurf value.
  JGaugeMesh* GaugeMesh; ///<Gauge object to compute zsurf values. //<vs_meeshdat>


  //-Zsurf data when UniformZsurf==false (calculated or from mesh data file).
  unsigned TimeCount;    ///<Number times with zsurf data.
  double*  Times;       ///<Times with data (minimum allocated size is 2). [TimeCount]
  float*   TimesZsurf;  ///<Zsurf data for different times and positions on CPU (minimum size is 2*Nptx). [Nptx*TimeCount]
  float*   TimesZsurfg; ///<Zsurf data for different times and positions on GPU (minimum size is 2*Nptx). [Nptx*TimeCount]

  double CurrentTime;      ///<Timestep of zsurf in Zsurf line.
  bool   CurrentExternal;  ///<Pointers CurrentZsurf and CurrentZsurfg are external from GaugeMesh object.
  float* CurrentZsurf;    ///<Zsurf data for timestep==CurrentTime on CPU. [Nptx]
  float* CurrentZsurfg;   ///<Zsurf data for timestep==CurrentTime on GPU. [Nptx]

  //-Variables with last search in Times[].
  double TimeStep;
  unsigned TimePosition;
  unsigned TimePositionNext;
  double TimePre;
  double TimeNext;
  double TimeFactor;

protected:
  void ComputeZsurfLine(bool forgauge,bool forceuniform);
  bool ConfigGaugeZsurf(JGaugeSystem* gaugesystem);

  void ConfigZsurfResults();
  void MeshToUniformData(jmsh::JMeshTDatas* mdatas);     //<vs_meeshdat>
  void MeshToZsurfLinePoints(jmsh::JMeshTDatas* mdatas); //<vs_meeshdat>

  void ResetTimes();
  void FindTime(double timestep);
  void InterpolateZsurfTime(double timestep,bool full);

public:
  JSphInOutZsurf(bool cpu,unsigned idzone,const StCteSph& csp
    ,tdouble3 direction,tdouble3 zoneposmin,tdouble3 zoneposmax);
  ~JSphInOutZsurf();
  void Reset();
  TpInZsurfMode ReadXml(const JXml* sxml,TiXmlElement* lis
    ,const std::string& dirdatafile,JGaugeSystem* gaugesystem);

  void GetConfig(std::vector<std::string>& lines)const;

  TpInZsurfMode GetZsurfMode()const{ return(ZsurfMode); }
  bool GetUniformZsurf()const{ return(UniformZsurf); }

  unsigned GetIdZone()const{ return(IdZone); }
  bool GetRemoveZsurf()const{ return(RemoveZsurf); }
  bool GetSvVtkZsurf()const{ return(SvVtkZsurf); }

  unsigned ComputeActivePoints(unsigned npt,const tdouble3* ptpos)const;
  void SetInitialPoints(unsigned npt,const tdouble3* ptpos,byte* ptok)const;

  float UpdateZsurf(double timestep);
  float GetCurrentZsurfUniform()const{ return(InputZsurf); }
  float GetCurrentZsurfNonUniform(const tdouble3& ps)const;

  tplane3d GetPlaDisx()const{ return(PlaDisx); }
  unsigned GetNptx()const{ return(Nptx); }
  const float* GetCurrentZsurfg()const{ return(CurrentZsurfg); }

  const StZsurfResult& GetZsurfResults()const;

//<vs_meeshdat_ini>
  //-Set uniform Zsurf.
  void RnSetZsurfUniform(double time0,double zsurf0,double time1,double zsurf1);
  //-Set non-uniform Zsurf.
  StRnZsurfData RnGetZsurfPtr(double time0,double time1);
//<vs_meeshdat_end>
};

#endif


