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
//:# - Contabiliza numero de interacciones entre particulas para calcular PIPS
//:#   (Particle Interactions Per Second). 03-06-2020
//:#############################################################################

/// \file JDsPips.h \brief Declares the class \ref JDsPips.

#ifndef _JDsPips_
#define _JDsPips_

#include "JObject.h"
#include "DualSphDef.h"
#include "JCellDivDataCpu.h"

#ifdef _WITHGPU
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>
#endif

#include <string>
#include <vector>

class JLog2;

//##############################################################################
//# JDsPips
//##############################################################################
/// \brief Count particle interactions to compute PIPS.

class JDsPips : protected JObject
{
public:
///Structure with the information of an interval.
typedef struct{
  unsigned nstep;  ///<Number of simulation step.
  double tstep;    ///<Physical time of simulation.
  double tsim;     ///<Runtime of simulation.
  ullong pirf;     ///<Real particle interactions between fluid particles with other particles.
  ullong pirb;     ///<Real particle interactions between bound particles with other particles.
  ullong picf;     ///<Checked particle interactions between fluid particles with other particles.
  ullong picb;     ///<Checked particle interactions between bound particles with other particles.
}StPipsInfo;

protected:
  JLog2* Log;

  unsigned NextNstep;
  std::vector<StPipsInfo> Data;
  unsigned NewData;

  unsigned SizeResultAux;
  ullong* ResultAux;  //-To copy final result from GPU memory. [SizeResultAux]


  double GetGPIs(unsigned cdata)const;
  double GetGPIsType(unsigned cdata,bool fluid)const;
  double GetCheckGPIsType(unsigned cdata,bool fluid)const;
  tdouble2 GetTotalPIs()const;

public:
  const bool Cpu;
  const unsigned StepsNum;  ///<Number of steps per interval to compute PIPS.
  const bool SvData;        ///<Store and save all data.
  const unsigned Ntimes;    ///<Interaction number per step (Verlet:1, Symplectic:2).

public:
  JDsPips(bool cpu,unsigned stepsnum,bool svdata,unsigned ntimes);
  ~JDsPips();
  long long GetAllocMemory()const;

  void SaveData();

  double GetGPIPS(double tsim)const;
  std::string GetTotalPIsInfo()const;

  bool CheckRun(unsigned nstep)const{ return(nstep>=NextNstep); }

  void ComputeCpu(unsigned nstep,double tstep,double tsim
    ,const StCteSph &csp,int ompthreads
    ,unsigned np,unsigned npb,unsigned npbok
    ,const StDivDataCpu &dvd,const unsigned *dcell,const tdouble3 *pos);

#ifdef _WITHGPU
  void ComputeGpu(unsigned nstep,double tstep,double tsim
    ,unsigned np,unsigned npb,unsigned npbok
    ,const StDivDataGpu &dvd,const unsigned *dcell,const float4 *poscell
    ,unsigned sauxmem,unsigned *auxmem);
#endif

};

#endif


