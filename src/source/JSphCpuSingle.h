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

/// \file JSphCpuSingle.h \brief Declares the class \ref JSphCpuSingle.

#ifndef _JSphCpuSingle_
#define _JSphCpuSingle_

#include "DualSphDef.h"
#include "JSphCpu.h"
#include <string>

class JCellDivCpuSingle;

//##############################################################################
//# JSphCpuSingle
//##############################################################################
/// \brief Defines the attributes and functions used only in Single-Core implementation.

class JSphCpuSingle : public JSphCpu
{
protected:
  JCellDivCpuSingle* CellDivSingle;

  llong GetAllocMemoryCpu()const;
  void UpdateMaxValues();
  void LoadConfig(const JSphCfgRun* cfg);
  void ConfigDomain();

  void ResizeParticlesSizeData(unsigned ndatacpu,unsigned newsize,unsigned minsize,float oversize,bool updatedivide);
  unsigned PeriodicMakeList(unsigned np,unsigned pini,bool stable,unsigned nmax,tdouble3 perinc,const tdouble3* pos,const typecode* code,unsigned* listp)const;
  void PeriodicDuplicatePos(unsigned pnew,unsigned pcopy,bool inverse,double dx,double dy,double dz,tuint3 cellmax,tdouble3* pos,unsigned* dcell)const;
  void PeriodicDuplicateVerlet(unsigned np,unsigned pini,tuint3 cellmax,tdouble3 perinc,const unsigned* listp
    ,unsigned* idp,typecode* code,unsigned* dcell,tdouble3* pos,tfloat4* velrho,tsymatrix3f* spstau,tfloat4* velrhom1)const;
  void PeriodicDuplicateSymplectic(unsigned np,unsigned pini,tuint3 cellmax,tdouble3 perinc
    ,const unsigned* listp,unsigned* idp,typecode* code,unsigned* dcell,tdouble3* pos
    ,tfloat4* velrho,tsymatrix3f* spstau,tdouble3* pospre,tfloat4* velrhopre)const;
  void PeriodicDuplicateNormals(unsigned np,unsigned pini,tuint3 cellmax
    ,tdouble3 perinc,const unsigned* listp,tfloat3* normals
    ,tfloat3* motionvel,tfloat3* motionace)const;
  void PeriodicSaveParent(unsigned np,unsigned pini,const unsigned* listp,unsigned* periparent)const;

  void PeriodicIgnore(unsigned np,typecode* code)const;
  void RunPeriodic();

  void RunCellDivide(bool updateperiodic);
  void AbortBoundOut();
  void SaveFluidOut();
  
  void MdbcBoundCorrection(TpInterStep interstep);
  void PreLoopProcedure(TpInterStep interstep);  //<vs_advshift>
  void ComputeFSParticles();                     //<vs_advshift>
  void ComputeUmbrellaRegion();                  //<vs_advshift>

  void Interaction_Forces(TpInterStep interstep);

  double ComputeAceMax()const;
  template<bool checkcode> double ComputeAceMaxSeq(unsigned np,const tfloat3* ace,const typecode* code)const;
  template<bool checkcode> double ComputeAceMaxOmp(unsigned np,const tfloat3* ace,const typecode* code)const;
  
  void RunInitialDDTRamp(); //<vs_ddramp>

  double ComputeStep(){ return(TStep==STEP_Verlet? ComputeStep_Ver(): ComputeStep_Sym()); }
  double ComputeStep_Ver();
  double ComputeStep_Sym();

  void RunFloating(double dt,bool predictor);
  inline tfloat3 FtPeriodicDist(const tdouble3& pos
    ,const tdouble3& center,float radius)const;
  void FtPartsSumAce(const tdouble3* posc,const tfloat3* acec
    ,const unsigned* ridpmot,tfloat6* acelinang)const;
  void FtPartsUpdate(double dt,bool updatenormals
    ,const tfloat6* fto_vellinang,const tdouble3* fto_center
    ,const unsigned* ridpmot,tdouble3* posc,tfloat4* velrhoc
    ,unsigned* dcellc,typecode* codec,tfloat3* boundnorc
    ,tfloat3* motionvelc,tfloat3* motionacec)const;

  void RunFirstGaugeSystem(double timestep);
  void RunGaugeSystem(double timestep);

  void ComputePips(bool run);
  
  void SaveData();
  void SaveExtraData();
  void FinishRun(bool stop);

  //<vs_flexstruc_ini>
  void FlexStrucInit();
  void UpdateFlexStrucGeometry();
  //<vs_flexstruc_end>

public:
  JSphCpuSingle();
  ~JSphCpuSingle();
  void Run(std::string appname,const JSphCfgRun* cfg,JLog2* log);

//-Code for InOut in JSphCpuSingle_InOut.cpp
//--------------------------------------------
protected:
  void InOutInit(double timestepini);
  void InOutIgnoreFluidDef(const std::vector<unsigned>& mkfluidlist,typecode* code);
  void InOutCheckProximity(unsigned newnp);
  void InOutComputeStep(double stepdt);
  void InOutUpdatePartsData(double timestepnew);
  void InOutExtrapolateData(unsigned inoutcount,const int* inoutpart);
};

#endif


