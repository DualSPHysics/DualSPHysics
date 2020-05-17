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

/// \file JCfgRun.h \brief Declares the class \ref JCfgRun.

#ifndef _JCfgRun_
#define _JCfgRun_

#include "DualSphDef.h"
#include "Functions.h"
#include "JObject.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

//##############################################################################
//# JCfgRun
//##############################################################################
/// \brief Defines the class responsible for collecting the execution parameters by command line.

class JCfgRun : protected JObject
{
public:
protected:
  bool SvDef;
  int DirsDef;
  void LoadDsphConfig(std::string path);

  static void LoadDouble3(std::string txopt,double def,tdouble3 &v1);
  static void LoadFloat3(std::string txopt,float def,tfloat3 &v1);
  static void LoadDouble6(std::string txopt,double def,tdouble3 &v1,tdouble3 &v2);
  static void LoadFloat6(std::string txopt,float def,tfloat3 &v1,tfloat3 &v2);

  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull,std::string &txopt1,std::string &txopt2,std::string &txopt3,std::string &txopt4)const;
  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull,std::string &txopt1,std::string &txopt2,std::string &txopt3)const{
    std::string tx4; SplitsOpts(opt,txword,txoptfull,txopt1,txopt2,txopt3,tx4);
  }
  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull,std::string &txopt1,std::string &txopt2)const{
    std::string tx3,tx4; SplitsOpts(opt,txword,txoptfull,txopt1,txopt2,tx3,tx4);
  }
  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull)const{
    std::string tx1,tx2,tx3,tx4; SplitsOpts(opt,txword,txoptfull,tx1,tx2,tx3,tx4);
  }

public:
  bool PrintInfo;

  bool Cpu;
  bool Gpu;
  int GpuId;
  bool GpuFree;
  bool Stable;
  int SvPosDouble;  ///<Saves particle position using double precision (default=0)

  int OmpThreads;

  TpCellMode  CellMode;
  int TBoundary;        ///<Boundary method: 0:None, 1:DBC (by default), 2:mDBC (SlipMode: 1:DBC vel=0)
  int SlipMode;         ///<Slip mode for mDBC: 0:None, 1:DBC vel=0, 2:No-slip, 3:Free slip (default=1).
  float MdbcThreshold;  ///<Kernel support limit to apply mDBC correction (default=0).
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

  //-General configuration from DsphConfig.xml
  bool CreateDirs;   ///<Creates full path for output files (true by default).
  bool CsvSepComa;   ///<Separator character in CSV files (false=semicolon, true=coma).

public:
  JCfgRun();
  void Reset();
  void VisuInfo()const;
  void VisuConfig()const;
  void LoadArgv(int argc,char** argv);
  void LoadFile(std::string fname,int lv);
  void LoadOpts(std::string *optlis,int optn,int lv,std::string file);
  void ErrorParm(const std::string &opt,int optc,int lv,const std::string &file)const;
};

#endif


