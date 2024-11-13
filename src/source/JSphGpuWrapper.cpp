#include "JSphGpuWrapper.h"
#include "JVResDataLimits.h"

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <limits>

#include "JSphBuffer.h"
#include "JSphGpu.h"
#include "JSphGpuSingle.h"
#include "JSphGpuSingle_Mr.h"
#include "JSphGpu_Buffer_iker.h"

using namespace std;

//==============================================================================
/// ???
//==============================================================================
void JSphGpuWrapper::InitProc(std::string appname,JSphCfgRun* cfg,JLog2* log){
  const string casename=cfg->CaseName;
  const string begindir=cfg->PartBeginDir;
  for(unsigned i=0;i<Num_SubDomains;i++){
    const string vrname=fun::PrintStr("_vres%02u",i);
    cfg->CaseName=casename+vrname;
    if(!begindir.empty())cfg->PartBeginDir=begindir+vrname;
    cfg->DirDataOut=string("data")+vrname;
    SubDomains[i].Init(appname,cfg,log,Num_SubDomains,i);
  }
  cfg->CaseName=casename;
  cfg->PartBeginDir=begindir;
  cfg->DirDataOut="data";
}

//==============================================================================
/// ???
//==============================================================================
void JSphGpuWrapper::BufferExtrapolateData(StInterParmsbg *parms) {
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].CallRunCellDivide();
    for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].BufferExtrapolateData(parms);
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].CallRunCellDivide();
}
void JSphGpuWrapper::ComputeStepBuffer(double dt, StInterParmsbg *parms,JCaseVRes& casemultires) {
    
    for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
    std::vector<std::vector<JMatrix4d>> mat=CalcBufferMotion(casemultires,dt);
	for(unsigned i=0;i<Num_SubDomains;i++){
        for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
        SubDomains[i].ComputeStepBuffer(dt,mat[i],parms);
    }
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].CallRunCellDivide();
    for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].BufferExtrapolateData(parms);
    
}

std::vector<std::vector<JMatrix4d>> JSphGpuWrapper::CalcBufferMotion(JCaseVRes& casemultires,double dt){
    // tdouble3* velmot= new tdouble3[NProc];
    std::vector<std::vector<JMatrix4d>> mat;
    
    //get the number of multi-resolution zones
    unsigned count=casemultires.Count();
    mat.resize(count);

    for (unsigned i = 0; i<count; i++){
        const JCaseVRes_Box* box=casemultires.GetZoneBox(i);
        if(box->Parent){           
            mat[i].push_back(GetBufferMotion(box,dt));
        }

        unsigned subcount=box->Count();
        for(unsigned j=0; j<subcount; j++) mat[i].push_back(GetBufferMotion(box->GetSubZoneBox(j),dt));   
    }   
    return(mat);
}

JMatrix4d JSphGpuWrapper::GetBufferMotion(const JCaseVRes_Box* box,double dt, int mkbound){
    //First check if the tracking is enable for the current zone, if not return (0,0,0)
    unsigned idzone=box->Id;
    bool IsTracking=box->TrackingIsActive();
    JMatrix4d mat;
    if (IsTracking){
        if(mkbound==0) mkbound=box->GetTrackingMkBound();
        //first check if in this domain there is the motion associated to the mkbound and if not return (0,0,0)
        mat=SubDomains[idzone].CalcVelMotion(mkbound,dt);
        //if not present in this domain let's try in the child domains by exploring the subdomains
        unsigned subcount=box->Count();
        for(unsigned j=0; j<subcount; j++) mat=GetBufferMotion(box->GetSubZoneBox(j),dt);
    }
    return(mat);
}

//==============================================================================
/// ???
//==============================================================================
void JSphGpuWrapper::LoadCaseVRes(const JSphCfgRun* cfg){
  const string filexml=cfg->CaseName+fun::PrintStr("_vres%02u",0)+".xml";
  CaseVRes.LoadFileXmlRun(filexml,"case.execution.vres");
  //-Saves data of VR zones for post-processing.
  JVResDataLimits vrlim;
  vrlim.LoadData(&CaseVRes);
  vrlim.SaveFile(cfg->DirOut);
}

//==============================================================================
/// ???
//==============================================================================
void JSphGpuWrapper::Run(std::string appname,JSphCfgRun* cfg,JLog2* log){
  LoadCaseVRes(cfg);
  Num_SubDomains=CaseVRes.Count();
  AllocateProc();
  InitProc(appname,cfg,log);
  for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].InitMultires(cfg,CaseVRes,i); 
  RunNormal(CaseVRes);
}

void JSphGpuWrapper::RunNormal(JCaseVRes& casemultires){
    std::vector<StInterParmsbg> parms;
     if (Num_SubDomains > 1) {
        for (unsigned i = 0; i < Num_SubDomains; i++) parms.push_back(SubDomains[i].getParms());
        for (unsigned i = 0; i < Num_SubDomains; i++) parms[i] = SubDomains[i].getParms();
        for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].BufferInit(parms.data());
        BufferExtrapolateData(parms.data());
    }
    if (Num_SubDomains > 1) ComputeStepBuffer(0.0, parms.data(),casemultires);

    double TimeMax = 0;

    for (unsigned i = 0; i < Num_SubDomains; i++) {
        TimeMax = SubDomains[i].Init2();
    }
    double TimeStep = SubDomains[0].GetTimeStep();
    double dt;
    double stepdt = 0;
    double *dtsync = new double[Num_SubDomains];
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
            ComputeStepBuffer(stepdt, parms.data(),casemultires);
        }
        TimeStep += stepdt;
        for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].Finish(stepdt);
        if(SubDomains[0].getNStepsBreak() && SubDomains[0].getNStep()>=SubDomains[0].getNStepsBreak())break;

    }
    for (unsigned i = 0; i < Num_SubDomains; i++) SubDomains[i].Finish2();

}
