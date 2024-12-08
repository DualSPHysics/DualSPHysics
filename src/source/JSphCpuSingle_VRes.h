#ifndef _JSphCpuSingle_VRes_
#define _JSphCpuSingle_VRes_



#include "JSphCpuSingle.h"
#include "JCellDivCpuSingle.h"
// #include "JSphGpu_Buffer_iker.h"
// #include "JSphGpu_cte.h"
#include "JSphVRes.h"
#include "JCaseVRes.h"
#include "JSphCpu.h"



class JSphVRes;
class JXml;




///Collects parameters for particle interaction on CPU.
inline stinterparmscb StInterparmscb(unsigned np,unsigned npb,unsigned npbok
  ,StDivDataCpu divdata,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat4 *velrhop,const unsigned *idp,const typecode *code, StCteSph csp
)
{
  stinterparmscb  d={np,npb,npbok,(np-npb)
    ,divdata,dcell
    ,pos,velrhop,idp,code,csp
  };
  return(d);
}


class JSphCpuSingle_VRes  :  public JSphCpuSingle
{
protected:
	JSphVRes* Multires;

//   StCteInteraction CTE;
  double SymplecticDtPre1;



  bool MRfastsingle=true;
  unsigned MROrder=1;
  float MRThreshold;




public:
	JSphCpuSingle_VRes();
	~JSphCpuSingle_VRes();
  void Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
    ,unsigned vrescount,unsigned vresid);
  void InitMultires(const JSphCfgRun *cfg, JCaseVRes casemultires, unsigned id);
  stinterparmscb getParms();
  void CallRunCellDivide();
  void BufferInit(stinterparmscb *parms);
  double Init2();
  double ComputeStep_SymB();
  double getSymplecticDtPre(){return SymplecticDtPre;};
      double getSymplecticDtPre1(){return SymplecticDtPre1;};
   void setSymplecticDtPre(double dt1){SymplecticDtPre=dt1;};
    void Finish(double dt1);
    void Finish2();
    void BufferShifting();
  int getNStep(){return Nstep;};
  int getNStepsBreak(){return NstepsBreak;};


  void BufferExtrapolateData(stinterparmscb *parms);
  void ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,stinterparmscb *parms);
  void ComputeUmbrellaRegionVRes();
  void ComputeFSParticlesVRes();
  void PreLoopProcedureVRes(TpInterStep interstep);


  double GetTimeStep(){return TimeStep;};
  void Interaction_ForcesB(TpInterStep interstep);

  JMatrix4d CalcVelMotion(unsigned trackingmk,double dt);
  JMatrix4d CalcMotionMoving(const StMotionData m,double dt);
  JMatrix4d CalcMotionFloating(const StFloatingData m,double dt);
    




};




#endif