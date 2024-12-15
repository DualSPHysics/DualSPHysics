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



template<typename T, typename TP>
class JSphVResDriver : protected JObject {
public:
  void Run(std::string appname,JSphCfgRun* cfg,JLog2* log);

private:
  std::vector<T> VResObj;    
  JCaseVRes CaseVRes;
  unsigned VResCount = 1;
  std::vector<std::vector<JMatrix4d>> CalcBufferMotion(JCaseVRes& casemultires, double dt);
  JMatrix4d GetBufferMotion(const JCaseVRes_Box* box, double dt, int mkbound = 0);
  void LoadCaseVRes(const JSphCfgRun* cfg);
  void InitProc(std::string appname, JSphCfgRun* cfg, JLog2* log);
  void SynchTimeStep();
  void BufferExtrapolateData(TP* parms);
  void ComputeStepBuffer(double dt,TP* parms,JCaseVRes& casemultires);
  void RunNormal(JCaseVRes& casemultires);
};

using namespace std;

//==============================================================================
template<typename T, typename TP>
void JSphVResDriver<T, TP>::InitProc(std::string appname, JSphCfgRun* cfg, JLog2* log) {
    const string casename=cfg->CaseName;
    const string begindir=cfg->PartBeginDir;

    for(unsigned i=0;i<VResCount;i++){
      const string vrname =fun::PrintStr("_vres%02u", i);
      cfg->CaseName=casename+vrname;
      if (!begindir.empty()) cfg->PartBeginDir=begindir+vrname;
      cfg->DirDataOut = string("data") + vrname;
      VResObj[i].Init(appname, cfg, log, VResCount, i);
    }
    cfg->CaseName     =casename;
    cfg->PartBeginDir =begindir;
    cfg->DirDataOut   ="data";
}

//==============================================================================
template<typename T, typename TP>
void JSphVResDriver<T, TP>::BufferExtrapolateData(TP* parms) {
    for (unsigned i = 0; i < VResCount; i++) VResObj[i].CallRunCellDivide();
    for (unsigned i = 0; i < VResCount; i++) parms[i] = VResObj[i].getParms();
    for (unsigned i = 0; i < VResCount; i++) VResObj[i].BufferExtrapolateData(parms);
    for (unsigned i = 0; i < VResCount; i++) VResObj[i].CallRunCellDivide();
}

template<typename T, typename TP>
void JSphVResDriver<T, TP>::ComputeStepBuffer(double dt, TP* parms, JCaseVRes& casemultires) {
    for (unsigned i = 0; i < VResCount; i++) parms[i] = VResObj[i].getParms();
    std::vector<std::vector<JMatrix4d>> mat = CalcBufferMotion(casemultires, dt);
    for (unsigned i = 0; i < VResCount; i++) {
        for (unsigned j = 0; j < VResCount; j++) parms[j] = VResObj[j].getParms();
        VResObj[i].ComputeStepBuffer(dt, mat[i], parms);
    }
    for (unsigned i = 0; i < VResCount; i++) VResObj[i].CallRunCellDivide();
    for (unsigned i = 0; i < VResCount; i++) parms[i] = VResObj[i].getParms();
    for (unsigned i = 0; i < VResCount; i++) VResObj[i].BufferExtrapolateData(parms);
}

template<typename T, typename TP>
std::vector<std::vector<JMatrix4d>> JSphVResDriver<T, TP>::CalcBufferMotion(JCaseVRes& casemultires, double dt) {
    std::vector<std::vector<JMatrix4d>> mat;
    unsigned count = casemultires.Count();
    mat.resize(count);

    for (unsigned i = 0; i < count; i++) {
        const JCaseVRes_Box* box = casemultires.GetZoneBox(i);
        if (box->Parent) {
            mat[i].push_back(GetBufferMotion(box, dt));
        }

        unsigned subcount = box->Count();
        for (unsigned j = 0; j < subcount; j++) mat[i].push_back(GetBufferMotion(box->GetSubZoneBox(j), dt));
    }
    return mat;
}

template<typename T, typename TP>
JMatrix4d JSphVResDriver<T, TP>::GetBufferMotion(const JCaseVRes_Box* box, double dt, int mkbound) {
    unsigned idzone = box->Id;
    bool IsTracking = box->TrackingIsActive();
    JMatrix4d mat;
    if (IsTracking) {
        if (mkbound == 0) mkbound = box->GetTrackingMkBound();
        // Check if idzone is within bounds
        if (idzone < VResCount) {
            mat = VResObj[idzone].CalcVelMotion(mkbound, dt);
        }
        unsigned subcount = box->Count();
        for (unsigned j = 0; j < subcount; j++) mat = GetBufferMotion(box->GetSubZoneBox(j), dt);
    }
    return mat;
}

//==============================================================================
template<typename T, typename TP>
void JSphVResDriver<T, TP>::LoadCaseVRes(const JSphCfgRun* cfg) {
    const string filexml = cfg->CaseName + fun::PrintStr("_vres%02u", 0) + ".xml";
    CaseVRes.LoadFileXmlRun(filexml, "case.execution.vres");
    //-Saves data of VR zones for post-processing.
    JVResDataLimits vrlim;
    vrlim.LoadData(&CaseVRes);
    vrlim.SaveFile(cfg->DirOut);
}

//==============================================================================
template<typename T, typename TP>
void JSphVResDriver<T, TP>::Run(std::string appname, JSphCfgRun* cfg, JLog2* log) {
  LoadCaseVRes(cfg);
  VResCount = CaseVRes.Count();
  // VResObj= new T[VResCount];
  VResObj.resize(VResCount);
  InitProc(appname, cfg, log);
  for (unsigned i = 0; i < VResCount; i++) VResObj[i].InitMultires(cfg, CaseVRes, i);
  RunNormal(CaseVRes);
  // for (unsigned i = 0; i < VResCount; i++){
  //   delete VResObj[i];
  // }
  // delete[] VResObj;
}

template<typename T, typename TP>
void JSphVResDriver<T, TP>::SynchTimeStep(){
  std::vector<double> dtsync(VResCount);

  for(unsigned i=0;i<VResCount;i++) dtsync[i]=VResObj[i].getSymplecticDtPre1();
  double dtmin = *std::min_element(dtsync.begin(), dtsync.end());
  for(unsigned i=0;i<VResCount;i++) VResObj[i].setSymplecticDtPre(dtmin);
}

template<typename T, typename TP>
void JSphVResDriver<T, TP>::RunNormal(JCaseVRes& casemultires) {
  std::vector<TP> parms;
  if (VResCount > 1) {
    for (unsigned i = 0; i < VResCount; i++) parms.push_back(VResObj[i].getParms());
    for (unsigned i = 0; i < VResCount; i++) parms[i] = VResObj[i].getParms();
    for (unsigned i = 0; i < VResCount; i++) VResObj[i].BufferInit(parms.data());
    BufferExtrapolateData(parms.data());
  }
    
  ComputeStepBuffer(0.0, parms.data(), casemultires);

  double TimeMax = 0;

  for(unsigned i=0;i<VResCount;i++)TimeMax = VResObj[i].Init2();

  double TimeStep=VResObj[0].GetTimeStep();
  double stepdt=0;

  SynchTimeStep();
  while(TimeStep<TimeMax){
    for(auto &vrobj : VResObj) stepdt=vrobj.ComputeStep_SymB();
    SynchTimeStep();
    ComputeStepBuffer(stepdt, parms.data(), casemultires);
    TimeStep += stepdt;
    for(unsigned i=0; i<VResCount;i++) VResObj[i].Finish(stepdt);
    if (VResObj[0].getNStepsBreak() && VResObj[0].getNStep() >= VResObj[0].getNStepsBreak()) break;
  }
  
  for (unsigned i = 0; i < VResCount; i++) VResObj[i].Finish2();
}



#endif /* JSPHGPUWRAPPER_H_ */
