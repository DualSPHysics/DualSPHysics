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

/// \file JSphCfgRun.h \brief Declares the class \ref JSphCfgRun.

#ifndef _JSphCfgRun_
#define _JSphCfgRun_

#include "JCfgRunBase.h"
#include "DualSphDef.h"

//##############################################################################
//# JSphCfgRun
//##############################################################################
/// \brief Defines the class responsible for collecting the execution parameters by command line.

class JSphCfgRun : public JCfgRunBase
{
protected:
  int DirsDef;

public:
  bool Cpu;
  bool Gpu;
  int GpuId;
  bool GpuFree;
  bool Stable;
  int SvPosDouble;  ///<Saves particle position using double precision (default=0)

  int OmpThreads;

  bool CellDomFixed;    ///<The Cell domain is fixed according maximum domain size.
  TpCellMode CellMode;  ///<Cell division mode.
  int TBoundary;        ///<Boundary method: 0:None, 1:DBC (by default), 2:mDBC (SlipMode: 1:DBC vel=0)
  int SlipMode;         ///<Slip mode for mDBC: 0:None, 1:DBC vel=0, 2:No-slip, 3:Free slip (default=1).
  int MdbcFastSingle;   ///<Matrix calculations are done in single precision (default=1). 
  float MdbcThreshold;  ///<Kernel support limit to apply mDBC correction (default=0).
  std::vector<std::string> InitParms;

  TpStep TStep;
  int VerletSteps;
  TpKernel TKernel;
  TpVisco TVisco;
  float Visco;
  float ViscoBoundFactor;
  double TimeMax,TimePart;
  int TDensity;   ///<Density Diffusion Term 0:None, 1:Molteni, 2:Fourtakas, 3:Fourtakas(full) (default=0)
  float DDTValue; ///<Value used with Density Diffusion Term (default=0.1)
  int Shifting;   ///<Shifting mode -1:no defined, 0:none, 1:nobound, 2:nofixed, 3:full
  bool SvRes,SvTimers,SvDomainVtk;
  bool Sv_Binx,Sv_Info,Sv_Csv,Sv_Vtk;
  std::string CaseName,RunName,DirOut,DirDataOut;
  std::string PartBeginDir;
  unsigned PartBegin,PartBeginFirst;
  float FtPause;
  bool RhopOutModif;              ///<Indicates whether \ref RhopOutMin or RhopOutMax is changed.
  float RhopOutMin,RhopOutMax;    ///<Limits for \ref RhopOut density correction.

  byte DomainMode;  ///<Domain configuration 0:No configured, 2:Fixed
  tdouble3 DomainFixedMin,DomainFixedMax;

  int NstepsBreak;  ///<Maximum number of steps allowed (debug).
  bool SvAllSteps;  ///<Saves a PART for each step (debug).

  unsigned PipsMode;   ///<Defines mode of PIPS calculation (0:No computed (default), 1:Computed, 2:computed and save detail).
  unsigned PipsSteps;  ///<Number of steps per interval to compute PIPS (100 by default).

public:
  JSphCfgRun();
  void Reset();
  void VisuInfo()const;
  void VisuConfig()const;
  void LoadOpts(std::string *optlis,int optn,int lv,const std::string &file);
  void ValidaCfg(){}
};

#endif


