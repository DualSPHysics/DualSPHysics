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

/// \file JSphVResDriver.h \brief Declares and implements the class \ref JSphVResDriver.

#ifndef _JSphVResDriver_
#define _JSphVResDriver_

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>

#include "JObject.h"
#include "JLog2.h"
#include "JSphVRes.h"
#include "JSphCfgRun.h"
#include "JCaseVRes.h"
#include "JMatrix4.h"
#include "JVResDataLimits.h"
#ifdef _WITHGPU
  #include "JSphGpu_VRes_iker.h"
#endif
//##############################################################################
//# JSphVResDriver
//##############################################################################
/// \brief Defines the attributes and functions used only in Variable resolution implementation.
template<typename T, typename TP>
class JSphVResDriver : protected JObject {
public:
  void Run(std::string appname,JSphCfgRun* cfg,JLog2* log);

private:
  std::vector<T> VResObj; 
  std::vector<TP> Parms;   
  JCaseVRes CaseVRes;
  unsigned VResCount  = 1;
  std::vector<std::vector<JMatrix4d>> CalcBufferMotion(JCaseVRes& casemultires, double dt);
  JMatrix4d GetBufferMotion(const JCaseVRes_Box* box, double dt, int mkbound = 0);
  void LoadCaseVRes(const JSphCfgRun* cfg);
  void InitProc(std::string appname, JSphCfgRun* cfg, JLog2* log);
  void LoadParms();
  void UpdateParms();
  double SynchTimeStep();
  void BufferExtrapolateData();
  void ComputeStepBuffer(double dt,JCaseVRes& casemultires);
  void RunVRes(JCaseVRes& casemultires);
};

// using namespace std;

//==============================================================================
/// Updates case configuration in VRes simulations.
//==============================================================================   
template<typename T, typename TP>
void JSphVResDriver<T, TP>::InitProc(std::string appname, JSphCfgRun* cfg, JLog2* log) {
    const std::string casename=cfg->CaseName;
    const std::string begindir=cfg->PartBeginDir;

    for(unsigned i=0;i<VResCount;i++){
      const std::string vrname =fun::PrintStr("_vres%02u", i);
      cfg->CaseName=casename+vrname;
      if (!begindir.empty()) cfg->PartBeginDir=begindir+vrname;
      cfg->DirDataOut = std::string("data") + vrname;
      VResObj[i].Init(appname, cfg, log, VResCount, i);
    }
    cfg->CaseName     =casename;
    cfg->PartBeginDir =begindir;
    cfg->DirDataOut   ="data";
}

//==============================================================================
/// Initial loading of parameters for VRes Coupling.
//============================================================================== 
template<typename T,typename TP>
void JSphVResDriver<T, TP>::LoadParms(){
  for(unsigned i=0;i<VResCount;i++) Parms.push_back(VResObj[i].GetVResParms());
}

//==============================================================================
/// Updates parameters for VRes Coupling.
//==============================================================================  
template<typename T,typename TP>
void JSphVResDriver<T, TP>::UpdateParms(){
  for(unsigned i=0;i<VResCount;i++) Parms[i]=VResObj[i].GetVResParms();
}

//==============================================================================
template<typename T,typename TP>
void JSphVResDriver<T, TP>::BufferExtrapolateData() {
  for(auto &vrobj : VResObj) vrobj.CallRunCellDivide();
  UpdateParms();
  for(auto &vrobj : VResObj) vrobj.BufferExtrapolateData(Parms.data());
  for(auto &vrobj : VResObj) vrobj.CallRunCellDivide();
}

//==============================================================================
/// Compute Step for buffer particles and obtain interpolation.
//==============================================================================   
template<typename T,typename TP>
void JSphVResDriver<T, TP>::ComputeStepBuffer(double dt, JCaseVRes& casemultires) {
  UpdateParms();
  std::vector<std::vector<JMatrix4d>> mat=CalcBufferMotion(casemultires, dt);
  for(unsigned i=0;i<VResCount;i++) {
    UpdateParms();
    VResObj[i].ComputeStepBuffer(dt,mat[i],Parms.data());
  }
  for(auto &vrobj : VResObj) vrobj.CallRunCellDivide();
  UpdateParms();
  for(auto &vrobj : VResObj) vrobj.BufferExtrapolateData(Parms.data());
}

//==============================================================================
/// Obtain movements information from VRes simulations
//==============================================================================   
template<typename T,typename TP>
std::vector<std::vector<JMatrix4d>> JSphVResDriver<T, TP>::CalcBufferMotion(JCaseVRes& casemultires,double dt) {
  std::vector<std::vector<JMatrix4d>> mat;
  unsigned count = casemultires.Count();
  mat.resize(count);

  for(unsigned i=0;i<count;i++){
    const JCaseVRes_Box* box = casemultires.GetZoneBox(i);
    if(box->Parent){
      mat[i].push_back(GetBufferMotion(box, dt));
    }
    unsigned subcount =box->Count();
    for (unsigned j=0;j<subcount;j++) mat[i].push_back(GetBufferMotion(box->GetSubZoneBox(j),dt));
  }
  return mat;
}

//==============================================================================
/// Obtain movements information from VRes simulations
//==============================================================================  
template<typename T,typename TP>
JMatrix4d JSphVResDriver<T, TP>::GetBufferMotion(const JCaseVRes_Box* box,double dt,int mkbound) {
  unsigned idzone = box->Id;
  bool IsTracking = box->TrackingIsActive();
  JMatrix4d mat;
  if(IsTracking){
    if(mkbound==0) mkbound=box->GetTrackingMkBound();
    // Check if idzone is within bounds
    if(idzone<VResCount) {
      mat=VResObj[idzone].CalcVelMotion(mkbound, dt);
    }
    unsigned subcount = box->Count();
    for (unsigned j=0;j<subcount;j++) mat=GetBufferMotion(box->GetSubZoneBox(j), dt);
  }
  return mat;
}

//==============================================================================
/// Load VRes configuration from xml.
//==============================================================================  
template<typename T,typename TP>
void JSphVResDriver<T, TP>::LoadCaseVRes(const JSphCfgRun* cfg) {
  const std::string filexml = cfg->CaseName + fun::PrintStr("_vres%02u", 0) + ".xml";
  CaseVRes.LoadFileXmlRun(filexml, "case.execution.vres");
  //-Saves data of VR zones for post-processing.
  JVResDataLimits vrlim;
  vrlim.LoadData(&CaseVRes);
  vrlim.SaveFile(cfg->DirOut);
}

//==============================================================================
/// Synch time-step between VRes simulations.
//==============================================================================  
template<typename T,typename TP>
double JSphVResDriver<T, TP>::SynchTimeStep(){
  std::vector<double> dtsync(VResCount);

  for(unsigned i=0;i<VResCount;i++) dtsync[i]=VResObj[i].getSymplecticDtPre1();
  double dtmin = *std::min_element(dtsync.begin(), dtsync.end());
  for(auto &vrobj : VResObj) vrobj.setSymplecticDtPre(dtmin);
  return(dtmin);
}

//==============================================================================
/// Main Loop for VRes simulation.
//==============================================================================  
template<typename T,typename TP>
void JSphVResDriver<T, TP>::RunVRes(JCaseVRes& casemultires) {
  LoadParms();
  
  for(auto &vrobj : VResObj) vrobj.BufferInit(Parms.data());
  BufferExtrapolateData();
  ComputeStepBuffer(0.0,casemultires);

  double TimeMax=0;
  double TimeStep=VResObj[0].GetTimeStep();
  double stepdt=0;

  for(unsigned i=0;i<VResCount;i++) TimeMax=VResObj[i].Init2();
  
  SynchTimeStep();
  while(TimeStep<TimeMax){
    for(auto &vrobj : VResObj) stepdt=vrobj.ComputeStepVRes();
    stepdt=VResObj[0].getSymplecticDtPre();
    const double dtmin=SynchTimeStep();
    ComputeStepBuffer(stepdt, casemultires);
    TimeStep += stepdt;
    for(auto &vrobj : VResObj) vrobj.Finish(stepdt);
    if (VResObj[0].getNStepsBreak() && VResObj[0].getNStep() >= VResObj[0].getNStepsBreak()) break;
  }  
  for(auto &vrobj : VResObj) vrobj.Finish2();
}

//==============================================================================
/// Perform VRes simulation.
//============================================================================== 
template<typename T,typename TP>
void JSphVResDriver<T, TP>::Run(std::string appname,JSphCfgRun* cfg,JLog2* log) {
  LoadCaseVRes(cfg);
  VResCount = CaseVRes.Count();
  VResObj.resize(VResCount);
  InitProc(appname, cfg, log);
  for(unsigned i=0;i<VResCount;i++) VResObj[i].VResInit(cfg,CaseVRes,i);
  RunVRes(CaseVRes);
}
#endif 
