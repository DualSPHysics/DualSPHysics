#ifndef _JSphGpuSingle_Mr_
#define _JSphGpuSingle_Mr_



#include "JSphGpuSingle.h"
#include "JCellDivGpuSingle.h"
#include "JSphGpu_Buffer_iker.h"
#include "JSphGpu_cte.h"
#include "JSphBuffer.h"
#include "JCaseVRes.h"
#include "JSphVResDef.h"



class JSphBuffer;
class JXml;



class JSphGpuSingle_Mr  :  public JSphGpuSingle
{
protected:
	JSphBuffer* Multires;

  StCteInteraction CTE;
  double SymplecticDtPre1;


  void CollectCTEdata();

  bool MRfastsingle=true;
  unsigned MROrder=1;
  float MRThreshold;

private:
  TpVresOrder VResOrder;
  TpVresMethod VResMethod;

  void LoadVResConfigParameters(const JSphCfgRun* cfg);


public:
	JSphGpuSingle_Mr();
	~JSphGpuSingle_Mr();
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