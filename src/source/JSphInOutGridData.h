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
//:# - Clase para gestionar la creacion de los puntos inlet. (25-01-2017)
//:#############################################################################

/// \file JSphInOutGrid.h \brief Declares the class \ref JSphInOutGrid.

#ifndef _JSphInOutGrid_
#define _JSphInOutGrid_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

class JLog2;

//##############################################################################
//# JSphInOutPointsParticles
//##############################################################################
/// \brief Stores data points in time.
class JSphInOutGridDataTime
{
private:
  double Time;
  float *Velx;
  float *Velz;

  inline void AllocData(bool usevelz);
  void ResetInit();

public:
  const unsigned Nx;
  const unsigned Nz;
  const unsigned Npt;
  JSphInOutGridDataTime(unsigned nx,unsigned nz);
  JSphInOutGridDataTime(unsigned nx,unsigned nz,double time,const float *velx,const float *velz);
  ~JSphInOutGridDataTime();

  void SetData(double time,const float *velx,const float *velz);
  void CopyFrom(double time,const JSphInOutGridDataTime *gdt){ SetData(time,gdt->Velx,gdt->Velz); }
  void Interpolate(double time,const JSphInOutGridDataTime *gdt,const JSphInOutGridDataTime *gdt2);

  double GetTime()const{ return(Time); }
  const float* GetVelx()const{ return(Velx); };
  const float* GetVelz()const{ return(Velz); };
};


//##############################################################################
//# JSphInOutGridData
//##############################################################################
/// \brief Defines object to manage interpolation grid points.
class JSphInOutGridData : protected JObject
{
private:
  JLog2 *Log;
  std::string File;

  unsigned Nx;
  unsigned Nz;
  unsigned Npt;
  double Dpx;
  double Dpz;
  bool UseVelz;

  tdouble3 PosMin;
  tdouble3 PosMax;

  std::vector<JSphInOutGridDataTime*> DataTimes;

  void LoadDataCsv(const std::string &filename);
  void LoadDataBin(const std::string &filename);

  unsigned SelCt;
  JSphInOutGridDataTime* SelData;

 #ifdef _WITHGPU
  double TimeGpu;
  unsigned CtVel0;  //-DataTime on GPU memory (CtVel0=SelCt).
  unsigned CtVel1;  //-DataTime on GPU memory (CtVel1=SelCt+1).
  unsigned CtSelVel;//-DataTime on GPU memory.

  float* Velx0g;    //-X velocity data on GPU at CtVel0.
  float* Velx1g;    //-X velocity data on GPU at CtVel1.
  float* SelVelxg;  //-X velocity data on GPU at requested time.

  float* Velz0g;    //-Z velocity data on GPU at CtVel0.
  float* Velz1g;    //-Z velocity data on GPU at CtVel1.
  float* SelVelzg;  //-Z velocity data on GPU at requested time.

  void AllocateMemoryGpu();
  void FreeMemoryGpu();
  void ComputeTimeGpu(double t);
 #endif

  void ComputeTime(double t);
  void SaveVtk(const JSphInOutGridDataTime *gdt,std::string filename)const;

public:
  const static unsigned FmtVersion=1;

  JSphInOutGridData();
  ~JSphInOutGridData();
  void Reset();
  void ConfigFromFile(const std::string &filename);
  void ConfigGridData(unsigned nx,unsigned nz,double dpx,double dpz,bool usevelz);
  void SetPosMin(const tdouble3 &posmin);

  void AddDataTime(double time,unsigned npt,const float *velx,const float *velz);

  unsigned GetNx()const{ return(Nx); }
  unsigned GetNz()const{ return(Nz); }
  double GetDpx()const{ return(Dpx); }
  double GetDpz()const{ return(Dpz); }
  bool GetUseVelz()const{ return(UseVelz); }

  unsigned CountTimes()const{ return(unsigned(DataTimes.size())); }
  std::string GetFile()const{ return(File); };

  void InterpolateVelCpu(double time,unsigned izone,unsigned np,const int *plist
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop
    ,float velcorr);
  void InterpolateZVelCpu(double time,unsigned izone,unsigned np,const int *plist
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop
    ,float velcorr);

#ifdef _WITHGPU
  void InterpolateZVelGpu(double time,unsigned izone,unsigned np,const int *plist
    ,const double2 *posxyg,const double *poszg,const typecode *codeg
    ,const unsigned *idpg,float4 *velrhopg,float velcorr);
#endif



  void SaveDataCsv(std::string filename)const;
  void SaveDataBin(std::string filename)const;

  void SaveDataVtk(std::string filename,int ctime=-1)const;
  void SaveDataVtkTime(std::string filename,double tmax,double dt);

  //void SaveVtkGrid(std::string filename,tfloat3 pos0)const;

};


#endif


