//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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

#include "Types.h"
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
  void LoadConfig(JCfgRun *cfg);
  void ConfigDomain();

  void ResizeParticlesSize(unsigned newsize,float oversize,bool updatedivide);
  unsigned PeriodicMakeList(unsigned np,unsigned pini,bool stable,unsigned nmax,tdouble3 perinc,const tdouble3 *pos,const typecode *code,unsigned *listp)const;
  void PeriodicDuplicatePos(unsigned pnew,unsigned pcopy,bool inverse,double dx,double dy,double dz,tuint3 cellmax,tdouble3 *pos,unsigned *dcell)const;
  void PeriodicDuplicateVerlet(unsigned np,unsigned pini,tuint3 cellmax,tdouble3 perinc,const unsigned *listp
    ,unsigned *idp,typecode *code,unsigned *dcell,tdouble3 *pos,tfloat4 *velrhop,tsymatrix3f *spstau,tfloat4 *velrhopm1)const;
  void PeriodicDuplicateSymplectic(unsigned np,unsigned pini,tuint3 cellmax,tdouble3 perinc,const unsigned *listp
    ,unsigned *idp,typecode *code,unsigned *dcell,tdouble3 *pos,tfloat4 *velrhop,tsymatrix3f *spstau,tdouble3 *pospre,tfloat4 *velrhoppre)const;
  void RunPeriodic();

  void RunCellDivide(bool updateperiodic);
  void AbortBoundOut();

  inline void GetInteractionCells(unsigned rcell
    ,int hdiv,const tint4 &nc,const tint3 &cellzero
    ,int &cxini,int &cxfin,int &yini,int &yfin,int &zini,int &zfin)const;

  void Interaction_Forces(TpInterStep tinterstep);
  
  double ComputeAceMax(unsigned np,const tfloat3* ace,const typecode *code)const;
  template<bool checkperiodic> double ComputeAceMaxSeq(unsigned np,const tfloat3* ace,const typecode *code)const;
  template<bool checkperiodic> double ComputeAceMaxOmp(unsigned np,const tfloat3* ace,const typecode *code)const;
  
  double ComputeStep(){ return(TStep==STEP_Verlet? ComputeStep_Ver(): ComputeStep_Sym()); }
  double ComputeStep_Ver();
  double ComputeStep_Sym();

  inline tfloat3 FtPeriodicDist(const tdouble3 &pos,const tdouble3 &center,float radius)const;
  void FtCalcForcesSum(unsigned cf,tfloat3 &face,tfloat3 &fomegaace)const;
  void FtCalcForces(StFtoForces *ftoforces)const;
  void FtCalcForcesRes(double dt,const StFtoForces *ftoforces,StFtoForcesRes *ftoforcesres)const;
  void RunFloating(double dt,bool predictor);
  void RunGaugeSystem(double timestep);
  
  void SaveData();
  void FinishRun(bool stop);

public:
  JSphCpuSingle();
  ~JSphCpuSingle();
  void Run(std::string appname,JCfgRun *cfg,JLog2 *log);

};

#endif


