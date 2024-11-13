#ifndef JSPHGPUWRAPPER_H_
#define JSPHGPUWRAPPER_H_

#include <string>

#include "JLog2.h"
#include "JSphBuffer.h"
#include "JSphCfgRun.h"
#include "JSphGpu.h"
#include "JSphGpuSingle_Mr.h"

class JSphGpuSingle_Mr;

class JSphGpuWrapper : protected JObject {
   public:
    void Run(std::string appname, JSphCfgRun *cfg, JLog2 *log);
    void BufferExtrapolateData(StInterParmsbg *parms);
    void BufferExtrapolateDataInter(unsigned nzone,StInterParmsbg *parms);
    void ComputeStepBuffer(double dt, StInterParmsbg *parms,JCaseVRes& casemultires);
    void RunNormal(JCaseVRes& casemultires);


   private:
    JSphGpuSingle_Mr *SubDomains;



    JCaseVRes CaseVRes;

    unsigned Num_SubDomains = 1;
    std::vector<std::vector<JMatrix4d>>  CalcBufferMotion(JCaseVRes& casemultires,double dt);
    JMatrix4d GetBufferMotion(const JCaseVRes_Box* box,double dt,int mkbound=0);
    void AllocateProc() { SubDomains = new JSphGpuSingle_Mr[Num_SubDomains]; };
    void LoadCaseVRes(const JSphCfgRun* cfg);
    void InitProc(std::string appname,JSphCfgRun* cfg,JLog2* log);

};

#endif /* JSPHGPUWRAPPER_H_ */