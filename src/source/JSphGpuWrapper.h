#ifndef JSPHGPUWRAPPER_H_
#define JSPHGPUWRAPPER_H_

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "JObject.h"
#include "JLog2.h"
#include "JSphBuffer.h"
#include "JSphCfgRun.h"
#include "JCaseVRes.h"
#include "JMatrix4.h"
#include "JVResDataLimits.h"



template<typename T, typename TP>
class JSphGpuWrapper : protected JObject {
public:
    void Run(std::string appname, JSphCfgRun* cfg, JLog2* log);
    void BufferExtrapolateData(TP* parms);
    void ComputeStepBuffer(double dt, TP* parms, JCaseVRes& casemultires);
    void RunNormal(JCaseVRes& casemultires);

private:
    T* SubDomains;
    JCaseVRes CaseVRes;
    unsigned Num_SubDomains = 1;
    std::vector<std::vector<JMatrix4d>> CalcBufferMotion(JCaseVRes& casemultires, double dt);
    JMatrix4d GetBufferMotion(const JCaseVRes_Box* box, double dt, int mkbound = 0);
    void AllocateProc() { SubDomains = new T[Num_SubDomains]; }
    void LoadCaseVRes(const JSphCfgRun* cfg);
    void InitProc(std::string appname, JSphCfgRun* cfg, JLog2* log);
};

using namespace std;

//==============================================================================
template<typename T, typename TP>
void JSphGpuWrapper<T, TP>::InitProc(std::string appname, JSphCfgRun* cfg, JLog2* log) {
    const string casename = cfg->CaseName;
    const string begindir = cfg->PartBeginDir;
    for (unsigned i = 0; i < Num_SubDomains; i++) {
        const string vrname = fun::PrintStr("_vres%02u", i);
        cfg->CaseName = casename + vrname;
        if (!begindir.empty()) cfg->PartBeginDir = begindir + vrname;
        cfg->DirDataOut = string("data") + vrname;
        SubDomains[i].Init(appname, cfg, log, Num_SubDomains, i);
    }
    cfg->CaseName = casename;
    cfg->PartBeginDir = begindir;
    cfg->DirDataOut = "data";
}

//==============================================================================
template<typename T, typename TP>
void JSphGpuWrapper<T, TP>::BufferExtrapolateData(TP* parms) {
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].CallRunCellDivide();
    for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].BufferExtrapolateData(parms);
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].CallRunCellDivide();
}

template<typename T, typename TP>
void JSphGpuWrapper<T, TP>::ComputeStepBuffer(double dt, TP* parms, JCaseVRes& casemultires) {
    for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
    std::vector<std::vector<JMatrix4d>> mat = CalcBufferMotion(casemultires, dt);
    for (unsigned i = 0; i < Num_SubDomains; i++) {
        for (unsigned j = 0; j < Num_SubDomains; j++) parms[j] = SubDomains[j].getParms();
        SubDomains[i].ComputeStepBuffer(dt, mat[i], parms);
    }
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].CallRunCellDivide();
    for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].BufferExtrapolateData(parms);
}

template<typename T, typename TP>
std::vector<std::vector<JMatrix4d>> JSphGpuWrapper<T, TP>::CalcBufferMotion(JCaseVRes& casemultires, double dt) {
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
JMatrix4d JSphGpuWrapper<T, TP>::GetBufferMotion(const JCaseVRes_Box* box, double dt, int mkbound) {
    unsigned idzone = box->Id;
    bool IsTracking = box->TrackingIsActive();
    JMatrix4d mat;
    if (IsTracking) {
        if (mkbound == 0) mkbound = box->GetTrackingMkBound();
        // Check if idzone is within bounds
        if (idzone < Num_SubDomains) {
            mat = SubDomains[idzone].CalcVelMotion(mkbound, dt);
        }
        unsigned subcount = box->Count();
        for (unsigned j = 0; j < subcount; j++) mat = GetBufferMotion(box->GetSubZoneBox(j), dt);
    }
    return mat;
}

//==============================================================================
template<typename T, typename TP>
void JSphGpuWrapper<T, TP>::LoadCaseVRes(const JSphCfgRun* cfg) {
    const string filexml = cfg->CaseName + fun::PrintStr("_vres%02u", 0) + ".xml";
    CaseVRes.LoadFileXmlRun(filexml, "case.execution.vres");
    //-Saves data of VR zones for post-processing.
    JVResDataLimits vrlim;
    vrlim.LoadData(&CaseVRes);
    vrlim.SaveFile(cfg->DirOut);
}

//==============================================================================
template<typename T, typename TP>
void JSphGpuWrapper<T, TP>::Run(std::string appname, JSphCfgRun* cfg, JLog2* log) {
    LoadCaseVRes(cfg);
    Num_SubDomains = CaseVRes.Count();
    AllocateProc();
    InitProc(appname, cfg, log);
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].InitMultires(cfg, CaseVRes, i);
    RunNormal(CaseVRes);
}

template<typename T, typename TP>
void JSphGpuWrapper<T, TP>::RunNormal(JCaseVRes& casemultires) {
    std::vector<TP> parms;
    if (Num_SubDomains > 1) {
        for (unsigned i = 0; i < Num_SubDomains; i++) parms.push_back(SubDomains[i].getParms());
        for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
        for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].BufferInit(parms.data());
        BufferExtrapolateData(parms.data());
    }
    if (Num_SubDomains > 1) ComputeStepBuffer(0.0, parms.data(), casemultires);

    double TimeMax = 0;

    for (unsigned i = 0; i < Num_SubDomains; i++) {
        TimeMax = SubDomains[i].Init2();
    }
    double TimeStep = SubDomains[0].GetTimeStep();
    double dt;
    double stepdt = 0;
    double* dtsync = new double[Num_SubDomains];
    double dtr = 0;
    for (unsigned i = 0; i < Num_SubDomains; i++) dtsync[i] = SubDomains[i].getSymplecticDtPre();
    double dtmin = 1000.0;
#undef min
    for (unsigned i = 0; i < Num_SubDomains; i++) dtmin = std::min(dtmin, dtsync[i]);
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].setSymplecticDtPre(dtmin);
    while (TimeStep < TimeMax) {
        for (unsigned i = 0; i < Num_SubDomains; i++) dt = SubDomains[i].ComputeStep_SymB();
        stepdt = dt;
        for (unsigned i = 0; i < Num_SubDomains; i++) dtsync[i] = SubDomains[i].getSymplecticDtPre1();
        dtmin = 1000.0;
#undef min
        for (unsigned i = 0; i < Num_SubDomains; i++) dtmin = std::min(dtmin, dtsync[i]);
        for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].setSymplecticDtPre(dtmin);
        if (Num_SubDomains > 1) {
            ComputeStepBuffer(stepdt, parms.data(), casemultires);
        }
        TimeStep += stepdt;
        for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].Finish(stepdt);
        if (SubDomains[0].getNStepsBreak() && SubDomains[0].getNStep() >= SubDomains[0].getNStepsBreak()) break;
    }
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].Finish2();

    delete[] dtsync;
}



#endif /* JSPHGPUWRAPPER_H_ */
