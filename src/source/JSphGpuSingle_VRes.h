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

/// \file JSphGpuSingle_VRes.h \brief Declares the class \ref JSphGpuSingle_VRes.


#ifndef _JSphGpuSingle_VRes_
#define _JSphGpuSingle_VRes_



#include "JSphGpuSingle.h"
#include "JCellDivGpuSingle.h"
#include "JSphGpu_VRes_iker.h"
#include "JSphGpu_cte.h"
#include "JSphVRes.h"
#include "JCaseVRes.h"
#include "JSphVResDef.h"



class JSphVRes;
class JXml;


//##############################################################################
//# JSphGpuSingle_VRes
//##############################################################################
/// \brief Defines the attributes and functions used only in Single-GPU Variable resolution implementation.
class JSphGpuSingle_VRes  :  public JSphGpuSingle
{
private:
  JSphVRes*     VRes;
  TpVresOrder   VResOrder;
  TpVresMethod  VResMethod;
  bool          VResFastSingle;
  float         VResThreshold;

  StCteInteraction CTE;
  double SymplecticDtPre1;

  void LoadVResConfigParameters(const JSphCfgRun* cfg);
  void VisuConfigVRes();
  void ComputeUmbrellaRegionVRes();
  void ComputeFSParticlesVRes();
  void PreLoopProcedureVRes(TpInterStep interstep);

  JMatrix4d CalcMotionMoving(const StMotionData m,double dt);
  JMatrix4d CalcMotionFloating(const StFloatingData m,double dt);

public:
	JSphGpuSingle_VRes();
	~JSphGpuSingle_VRes();
  void Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
    ,unsigned vrescount,unsigned vresid);
  void VResInit(const JSphCfgRun *cfg, JCaseVRes casemultires, unsigned id);
  StInterParmsbg GetVResParms();
  void BufferInit(StInterParmsbg *parms);
  double Init2();

  void AddWarningVRes();

  double ComputeStepVRes();
  
  void Finish(double dt1);
  void Finish2();
  void BufferShifting();
  
  double  getSymplecticDtPre(){return SymplecticDtPre;};
  double  getSymplecticDtPre1(){return SymplecticDtPre1;};
  void    setSymplecticDtPre(double dt1){SymplecticDtPre=dt1;};
  double  GetTimeStep(){return TimeStep;};
  int     getNStep(){return Nstep;};
  int     getNStepsBreak(){return NstepsBreak;};


  void BufferExtrapolateData(StInterParmsbg *parms);
  void ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,StInterParmsbg *parms);
  void CallRunCellDivide();
  
  JMatrix4d CalcVelMotion(unsigned trackingmk,double dt);

    




};




#endif