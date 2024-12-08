#ifndef _JSphGpuSingle_VRes_
#define _JSphGpuSingle_VRes_



#include "JSphGpuSingle.h"
#include "JCellDivGpuSingle.h"
#include "JSphGpu_Buffer_iker.h"
#include "JSphGpu_cte.h"
#include "JSphVRes.h"
#include "JCaseVRes.h"
#include "JSphVResDef.h"



class JSphVRes;
class JXml;



class JSphGpuSingle_VRes  :  public JSphGpuSingle
{
protected:
	JSphVRes* Multires;

  


  void CollectCTEdata();

  

private:
  TpVresOrder VResOrder;
  TpVresMethod VResMethod;

  void LoadVResConfigParameters(const JSphCfgRun* cfg);
  void ComputeUmbrellaRegionVRes();
  void ComputeFSParticlesVRes();
  void PreLoopProcedureVRes(TpInterStep interstep);

  JMatrix4d CalcMotionMoving(const StMotionData m,double dt);
  JMatrix4d CalcMotionFloating(const StFloatingData m,double dt);


  bool MRfastsingle=true;
  unsigned MROrder=1;
  float MRThreshold;

  StCteInteraction CTE;
  double SymplecticDtPre1;


public:
	JSphGpuSingle_VRes();
	~JSphGpuSingle_VRes();
  void Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
    ,unsigned vrescount,unsigned vresid);
  void InitMultires(const JSphCfgRun *cfg, JCaseVRes casemultires, unsigned id);
  StInterParmsbg getParms();
  void CallRunCellDivide();
  void BufferInit(StInterParmsbg *parms);
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


  void BufferExtrapolateData(StInterParmsbg *parms);
  void ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,StInterParmsbg *parms);
  


  double GetTimeStep(){return TimeStep;};
  void Interaction_ForcesB(TpInterStep interstep);

  JMatrix4d CalcVelMotion(unsigned trackingmk,double dt);

    




};




#endif