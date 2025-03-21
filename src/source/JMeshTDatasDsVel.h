//HEAD_DSCODES
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

//:#############################################################################
//:# Descripcion:
//:# =============
//:# Clase derivada de JMeshTDatas para interpolacion de datos para el inlet.
//:# - Implementacion. (24-08-2020)
//:#############################################################################

/// \file JMeshTDatasDsVel.h \brief Declares the class \ref JMeshTDatasDsVel.

#ifndef _JMeshTDatasDsVel_
#define _JMeshTDatasDsVel_

#include "JMeshTDatas.h"
#include "JMeshDataDef.h"
#include "DualSphDef.h"
#include "JSphInOutDef.h"

#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

#include <string>
#include <vector>
#include <fstream>
#include <cfloat>

class JDataArrays;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshTDatasDsVel
//##############################################################################
/// \brief Allows loading and processing mesh data for inlet/outlet code.

class JMeshTDatasDsVel : public JMeshTDatas
{
 public:
   const bool Cpu;
 protected:
   tdouble3 SetPos;     ///<Position offset (applied after filters).
   double InitialTime;  ///<Defines initial time of data.
   double LoopTbeginRq; ///<Requested begin time of loop (DBL_MAX=disabled).
   double LoopTmaxRq;   ///<Requested final time of loop (DBL_MAX=disabled).
   bool VelMagnitude;   ///<Only magnitude of velocity is used.
   tfloat3 VelDir;      ///<Direction vector for velocity.
   bool VelReverse;     ///<Reverses velocity data (v=-v).
   tdouble3 SetVelMul;  ///<Multiply value to velocity data (v*=v2).
   tdouble3 SetVelAdd;  ///<Add value to velocity data (v+=v2).

   //-Variables for GPU calculations.
#ifdef _WITHGPU
   unsigned GpuCtime1;
   unsigned GpuCtime2;

   float*  Vel1Data1g;
   float*  Vel1Data2g;
   float*  Vel1DataTg;

   float3* Vel3Data1g;
   float3* Vel3Data2g;
   float3* Vel3DataTg;
#endif

 protected:
#ifdef _WITHGPU
  void ResetGpu();
  void FreeMemoryGpu();
  void AllocMemoryGpu();
  void IntpComputeTimeGpu(double timestep);
#endif

 public:
  JMeshTDatasDsVel(std::string appname,bool cpu);
  ~JMeshTDatasDsVel();
  void Reset();
  void ConfigVel(std::string file,tdouble3 setpos
    ,double initialtime,double looptmax,double looptbegin
    ,bool magnitude,tfloat3 veldir,bool reverse
    ,tdouble3 velmul,tdouble3 veladd);

  double GetInitialTime()const{ return(InitialTime); }
  double GetLoopTbeginRq()const{ return(LoopTbeginRq); };
  double GetLoopTmaxRq()const{ return(LoopTmaxRq); };

  tdouble3 GetSetPos()const{ return(SetPos); }
  bool GetVelMagnitude()const{ return(VelMagnitude); }
  bool GetVelReverse()const{ return(VelReverse); }
  tdouble3 GetSetVelMul()const{ return(SetVelMul); }
  tdouble3 GetSetVelAdd()const{ return(SetVelAdd); }

  void InterpolateInOutVelCpu(double timestep,unsigned izone,unsigned np,const int* plist
    ,const tdouble3* pos,const typecode* code,const unsigned* idp
    ,tfloat4* velrhop,float velcorr);

  void InterpolateRzVelCpu(double t,byte zoneid,unsigned np,unsigned pini
    ,const byte* rzid,const float* rzfactor,const tdouble3* pos,tfloat4* velrhop);

#ifdef _WITHGPU
  void InterpolateInOutVelGpu(double timestep,unsigned izone
    ,unsigned np,const int* plistg,const double2* posxyg
    ,const double* poszg,const typecode* codeg
    ,const unsigned* idpg,float4* velrhopg,float velcorr);

#endif

  StRnVelData RnGetVelPtr(double time0,double time1);

};



}


#endif


