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

/// \file JSphCpuSingle_VRes.h \brief Declares the class \ref JSphCpuSingle_VRes.

#ifndef _JSphCpuSingle_VRes_
#define _JSphCpuSingle_VRes_

#include "JSphCpuSingle.h"
#include "JCellDivCpuSingle.h"
#include "JSphVRes.h"
#include "JCaseVRes.h"
#include "JSphCpu.h"
#include "JSphCpu_VRes.h"
#include "JSphVResDef.h"

class JSphVRes;
class JXml;

///Collects parameters for particle interaction on CPU.
inline stinterparmscb StInterparmscb(unsigned np,unsigned npb,unsigned npbok
  ,StDivDataCpu divdata,const unsigned *dcell ,const tdouble3 *pos
  ,const tfloat4 *velrhop,const unsigned *idp,const typecode *code,StCteSph csp){
  stinterparmscb  d={np,npb,npbok,np-npb
    ,divdata,dcell,pos,velrhop,idp,code,csp};
  return(d);
}

//##############################################################################
//# JSphGpuSingle_VRes
//##############################################################################
/// \brief Defines the attributes and functions used only in Single-CPU Variable resolution implementation.
class JSphCpuSingle_VRes  :  public JSphCpuSingle
{
protected:
	TpVresOrder   VResOrder;
  TpVresMethod  VResMethod;
  JSphVRes*     VRes;
  float         VResThreshold;

  double SymplecticDtPre1;

  void LoadVResConfigParameters(const JSphCfgRun* cfg);
  void VisuConfigVRes();
  void ComputeUmbrellaRegionVRes();
  void ComputeFSParticlesVRes();
  void PreLoopProcedureVRes(TpInterStep interstep);
  void BufferShifting();

  JMatrix4d CalcMotionMoving(const StMotionData m,double dt);
  JMatrix4d CalcMotionFloating(const StFloatingData m,double dt);

public:
	JSphCpuSingle_VRes();
	~JSphCpuSingle_VRes();
  void Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
    ,unsigned vrescount,unsigned vresid);
  double Init2();
  void AddWarningVRes();
  double ComputeStepVRes();
  void Finish(double dt1);
  void Finish2();
  void VResInit(const JSphCfgRun *cfg, JCaseVRes casemultires, unsigned id);

  stinterparmscb GetVResParms();

  void CallRunCellDivide();
  void BufferInit(stinterparmscb *parms);
  
  double  getSymplecticDtPre(){return SymplecticDtPre;};
  double  getSymplecticDtPre1(){return SymplecticDtPre1;};
  void    setSymplecticDtPre(double dt1){SymplecticDtPre=dt1;};  
  int     getNStep(){return Nstep;};
  int     getNStepsBreak(){return NstepsBreak;};
  double  GetTimeStep(){return TimeStep;};

  void BufferExtrapolateData(stinterparmscb *parms);
  void ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,stinterparmscb *parms);

  JMatrix4d CalcVelMotion(unsigned trackingmk,double dt);
};




#endif