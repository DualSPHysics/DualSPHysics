#include "JSphCpuSingle_VRes.h"
#include "JSphCpuSingle.h"
#include "JSphCpu_VRes.h"
#include "JSphVRes.h"
#include "JTimeControl.h"
#include "JDsOutputTime.h"
#include "JCellDivCpuSingle.h"
#include "JSphInOut.h"
#include "JXml.h"
#include "JDsMotion.h"
#include "JDsPartMotionSave.h"
#include "JDsPartFloatSave.h"
#include "JSphShifting.h"
#include "JSphShiftingAdv.h"
#include "JSph.h"



using namespace std;
JSphCpuSingle_VRes::JSphCpuSingle_VRes():JSphCpuSingle(){
  ClassName="JSphGpuSingle";
  MRfastsingle=true;
  MROrder=0;
  MRThreshold=100;
  Multires=NULL;


}



JSphCpuSingle_VRes::~JSphCpuSingle_VRes(){
  DestructorActive=true;
  delete Multires; Multires=NULL;
}

//==============================================================================
/// ???
//==============================================================================
void JSphCpuSingle_VRes::Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
  ,unsigned vrescount,unsigned vresid)
{
  if(!cfg || !log)return;
  AppName=appname; Log=log; CfgRun=cfg;
  VResCount=vrescount;
  VResId=vresid; 

  


  if(cfg->MRFastSingle>=0)MRfastsingle=(cfg->MRFastSingle>0);
  if(cfg->MROrder>=0)MROrder=cfg->MROrder;
  if(cfg->MRThreshold >=0)MRThreshold=cfg->MRThreshold;

  //-Creates array system for particles.
  Arrays_Cpu=new JArraysCpu(Log);

  //-Configure timers.
  Timersc->Config(cfg->SvTimers);
  Timersc->TmStart(TMC_Init);

  //-Load parameters and input data.
  LoadConfig(cfg);
  LoadCaseParticles();
  VisuConfig();
  ConfigDomain();
  ConfigRunMode();
  VisuParticleSummary();

  //-Initialisation of execution variables.
  InitRunCpu();
  RunFirstGaugeSystem(TimeStep);
  if(InOut)InOutInit(TimeStepIni);
  FreePartsInit();
  PrintAllocMemory(GetAllocMemoryCpu());
  UpdateMaxValues();
  SaveData();
  Arrays_Cpu->FreeUnusedArrays();
  Timersc->ResetTimes();
  Timersc->TmStop(TMC_Init);
  if(Log->WarningCount())Log->PrintWarningList("\n[WARNINGS]","");
  PartNstep=0; Part++;
}

double JSphCpuSingle_VRes::Init2(){
	 //-Main Loop.
	JTimeControl tc("30,60,300,600");//-Shows information at 0.5, 1, 5 y 10 minutes (before first PART).
  bool minfluidstopped=false;
  TimerSim.Start();
  TimerPart.Start();
  Log->Print(string("\n[Initialising simulation (")+RunCode+")  "+fun::GetDateTime()+"]");
  if(DsPips)ComputePips(true);
  PrintHeadPart();
	return(TimeMax);
}


double JSphCpuSingle_VRes::ComputeStep_SymB(){
  const double dt=SymplecticDtPre;
  if(CaseNmoving)CalcMotion(dt);               //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
  InterStep=INTERSTEP_SymPredictor;
  DemDtForce=dt*0.5f;
  MdbcBoundCorrection(INTERSTEP_SymPredictor);  //-Mdbc correction
  PreInteraction_Forces(INTERSTEP_SymPredictor);//-Allocating temporary arrays.
  PreLoopProcedureVRes(INTERSTEP_SymPredictor);     //-Calculate variables for interaction forces (Shifting,DDT,etc...).                          //-For DEM interaction.
  Interaction_Forces(INTERSTEP_SymPredictor);  //-Interaction.
  const double dt_p=DtVariable(false);         //-Calculate dt of predictor step.
  if(Shifting)RunShifting(dt*.5);              //-Shifting.
  ComputeSymplecticPre(dt);                    //-Apply Symplectic-Predictor to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt*.5,true);       //-Control of floating bodies.
  PosInteraction_Forces();                     //-Free memory used for interaction.
  //-Corrector
  //-----------
  InterStep=INTERSTEP_SymCorrector;
  DemDtForce=dt;                               //-For DEM interaction.
  RunCellDivide(true);
  MdbcBoundCorrection(INTERSTEP_SymCorrector);  //-Mdbc correction
  PreInteraction_Forces(INTERSTEP_SymCorrector);//-Allocating temporary arrays.
  PreLoopProcedureVRes(INTERSTEP_SymCorrector);     //-Calculate variables for interaction forces (Shifting,DDT,etc...).
  Interaction_Forces(INTERSTEP_SymCorrector);  //-Interaction.
  const double dt_c=DtVariable(true);          //-Calculate dt of corrector step.
  if(Shifting)RunShifting(dt);                 //-Shifting.
  ComputeSymplecticCorr(dt);                   //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);         //-Control of floating bodies.
  PosInteraction_Forces();                     //-Free memory used for interaction.
  if(Damping)RunDamping(dt);                   //-Applies Damping.
  if(RelaxZones)RunRelaxZone(dt);              //-Generate waves using RZ.
  SymplecticDtPre1=min(dt_p,dt_c);              //-Calculate dt for next ComputeStep.
  return(dt);
}


void JSphCpuSingle_VRes::Finish(double dt1){
	    RunGaugeSystem(TimeStep+dt1);
	    if(CaseNmoving)RunMotion(dt1);
	    if(InOut)InOutComputeStep(dt1);
	    // else RunCellDivide(true);
	    TimeStep+=dt1;
	    LastDt=dt1;
      Nstep++;
    //-Save extra PART data.
    if(TimeOutExtraCheck(TimeStep)){
      if(PartFloatSave)PartFloatSave->AddDataExtra(Part,TimeStep,Nstep);
      if(PartMotionSave)PartMotionSave->AddDataExtraCpu(Part,TimeStep,Nstep,Np
        ,Pos_c->cptr(),RidpMot);
      TimeOutExtraUpdate(TimeStep);
    }
	    if(TimeStep>=TimePartNext ){
        // MultiRes->SaveFluxesData(DirDataOut,Part);
	    	if(VRES_DG_SAVEVTK){
          DgSaveVtkParticlesCpuMRBuffer("Compute_Debug_Buffer.vtk",Part,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL);
	        DgSaveVtkParticlesCpuMR("Compute_Debug_Fluid.vtk",Part,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL,NULL);
        }
        if(VRES_DG_SAVEVTK)Multires->SaveNormals(DirDataOut+"NormalsBuffer_",Part);
        if(VRES_DG_SAVEVTK)Multires->SaveVtkDomains(DirDataOut+"Domains",Part,Simulate2D);
        if(VRES_DG_SAVEVTK)Multires->SaveVResData(Part,TimeStep,Nstep);

	      SaveData();
	      Part++;
	      PartNstep=Nstep;
	      TimeStepM1=TimeStep;
        TimePartNext=(SvAllSteps? TimeStep: OutputTime->GetNextTime(TimeStep));
	      TimerPart.Start();
	    }
	    UpdateMaxValues();
}

void JSphCpuSingle_VRes::Finish2(){
	TimerSim.Stop(); TimerTot.Stop();

	  //-End of Simulation.
	  //--------------------
	  FinishRun(false);

}

void JSphCpuSingle_VRes::InitMultires(const JSphCfgRun* cfg, JCaseVRes casemultires, unsigned id){
   Multires=new JSphVRes(Cpu,CSP,casemultires,id,AppName,DirDataOut,PartBegin,PartBeginDir);
   Multires->Config();
}

stinterparmscb JSphCpuSingle_VRes::getParms()
{
	 stinterparmscb parms=StInterparmscb(Np,Npb,NpbOk
	    ,DivData,Dcell_c->cptr()
	    ,Pos_c->cptr(),Velrho_c->cptr(),Idp_c->cptr(),Code_c->cptr(),CSP
	  );

	return(parms);
}



void JSphCpuSingle_VRes::BufferInit(stinterparmscb *parms){

  for(unsigned i=0;i<Multires->GetCount();i++){
      acint bufferpartc("-",Arrays_Cpu,true);
      unsigned buffercountpre;

  //		DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
      buffercountpre=Multires->CreateListCpuInit(Np,0,Pos_c->cptr(),Idp_c->cptr(),Code_c->ptr(),bufferpartc.ptr(),i);
      // std::cout << buffercountpre << std::endl;

    }

    RunCellDivide(true);
    if(VRES_DG_SAVEVTK)DgSaveVtkParticlesCpuMRBuffer("Compute_step",Nstep,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL);

}

void JSphCpuSingle_VRes::BufferExtrapolateData(stinterparmscb *parms){

	for(unsigned i=0;i<Multires->GetCount();i++){
		
   acint bufferpartc("-",Arrays_Cpu,true);
      unsigned buffercountpre;

  //		DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
      buffercountpre=Multires->CreateListCpuInit(Np,0,Pos_c->cptr(),Idp_c->cptr(),Code_c->ptr(),bufferpartc.ptr(),i);
		
		unsigned id=Multires->GetZone(i)->getZone()-1;


	 fvres::Interaction_BufferExtrap(buffercountpre,bufferpartc.ptr(),parms[id],Pos_c->ptr(),Velrho_c->ptr(),Code_c->ptr(),2);
    // DgSaveVtkParticlesCpuMRBuffer("Debug_Buffer_CpuInit.vtk",Nstep,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL);
	  // DgSaveVtkParticlesCpuMR("Debug_Multi_CpuInit.vtk",Nstep,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL,NULL);
		

		
	}
}

void JSphCpuSingle_VRes::ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,stinterparmscb *parms){
	// // TmgStart(Timers,TMG_SuBuffer);
    Multires->UpdateMatMov(mat);
	  Multires->MoveBufferZone(dt,mat);
	// // TmgStart(Timers,TMG_SuBuffer);

	// // Multires->MoveBufferZone(dt,velmot);
	
	for(unsigned i=0;i<Multires->GetCount();i++){
    
    acint bufferpartc("-",Arrays_Cpu,true);
    unsigned buffercountpre;

    buffercountpre=Multires->CreateListCpuInit(Np,0,Pos_c->cptr(),Idp_c->cptr(),Code_c->ptr(),bufferpartc.ptr(),i);

    StrDataVresCpu vresdata=Multires->GetZoneFluxInfoCpu(i);
			
  //   // cusphvres::MoveBufferZone(nini,ntot,posxy,posz,dt,velmot[i]);


		unsigned id=Multires->GetZone(i)->getZone()-1;
		
    fvres::Interaction_BufferExtrapFlux(vresdata.ntot,vresdata.nini,parms[id],vresdata.points,vresdata.normals,vresdata.velmot,vresdata.mass,2,Dp,dt,100);

		fvres::CheckMassFlux(vresdata.ntot,vresdata.nini,CSP,DivData
      ,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
      ,vresdata.points,vresdata.normals,vresdata.mass);
    
    unsigned newnp=Multires->ComputeStepCpu(buffercountpre,bufferpartc.ptr(),Code_c->ptr(),Pos_c->cptr(),i);

    if(newnp){
      if(!CheckCpuParticlesSize(Np+newnp)){
        const unsigned ndatacpu=Np;
        ResizeParticlesSizeData(ndatacpu,Np+newnp,Np+newnp, 0.2,true);
        CellDivSingle->SetIncreaseNp(newnp);
      }
      Multires->CreateNewPart(IdMax+1,Dcell_c->ptr(),Code_c->ptr(),Pos_c->ptr(),Idp_c->ptr(),Velrho_c->ptr()
        ,this,Np,i);

    }

    if(newnp){      
      if(SpsTauRho2_c)  SpsTauRho2_c->MemsetOffset(Np,0,newnp);
      if(BoundNor_c)    BoundNor_c->MemsetOffset(Np,0,newnp);
      if(FSType_c)      FSType_c->MemsetOffset(Np,0,newnp);
      if(ShiftVel_c)    ShiftVel_c->MemsetOffset(Np,0,newnp);
      Np+=newnp; 
      TotalNp+=newnp;
      IdMax=unsigned(TotalNp-1);
		} 
	}
}

void JSphCpuSingle_VRes::BufferShifting(){
  // StrGeomVresGpu& vresgdata=Multires->GetGeomInfoVres();
	// cusphvres::BufferShiftingGpu(Np,Npb,Posxy_g->ptr(),Posz_g->ptr(),ShiftVel_g->ptr(),Code_g->ptr(),vresgdata,NULL);
}



void JSphCpuSingle_VRes::CallRunCellDivide(){
    RunCellDivide(true);
}


//==============================================================================
/// PreLoop for additional models computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphCpuSingle_VRes::PreLoopProcedureVRes(TpInterStep interstep){
  const bool runshift=(ShiftingAdv && interstep==INTERSTEP_SymPredictor && Nstep!=0);
  if(runshift){
    Timersc->TmStart(TMC_SuShifting);
    ComputeFSParticlesVRes();
    ComputeUmbrellaRegionVRes();

    StrGeomVresCpu vresdata=Multires->GetGeomInfoVresCpu();
    fvres::CallCorrectShiftBuff(Np,Npb,CSP,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,ShiftVel_c->ptr(),FSType_c->ptr(),vresdata);

    PreLoopInteraction_ct(DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
      ,Velrho_c->cptr(),FSType_c->ptr(),ShiftVel_c->ptr(),FSNormal_c->ptr()
      ,FSMinDist_c->ptr());
    ComputeShiftingVel(Simulate2D,ShiftingAdv->GetShiftCoef()
      ,ShiftingAdv->GetAleActive(),SymplecticDtPre,FSType_c->ptr()
      ,FSNormal_c->ptr(),FSMinDist_c->ptr(),ShiftVel_c->ptr());
    fvres::BufferShiftingCpu(Np,Npb,Pos_c->cptr(),ShiftVel_c->ptr(),Code_c->ptr(),vresdata);
    //-Updates pre-loop variables in periodic particles.
    if(PeriParent_c){
      const unsigned* periparent=PeriParent_c->ptr();
      unsigned* fstype  =FSType_c->ptr();
      tfloat4*  shiftvel=ShiftVel_c->ptr();
      for(unsigned p=Npb;p<Np;p++)if(periparent[p]!=UINT_MAX){
        fstype[p]  =fstype[periparent[p]];
        shiftvel[p]=shiftvel[periparent[p]];
      }
    }
    //-Saves VTK for debug.
    if(1 && runshift && TimeStep+LastDt>=TimePartNext)DgSaveVtkParticlesCpu("Compute_FreeSurface_",Part,0,Np,Pos_c->cptr()
      ,Code_c->cptr(),FSType_c->cptr(),ShiftVel_c->cptr(),FSNormal_c->cptr());
    Timersc->TmStop(TMC_SuShifting);
  }
}


//==============================================================================
/// Compute free-surface particles and their normals.
//==============================================================================
void JSphCpuSingle_VRes::ComputeFSParticlesVRes(){
  StrGeomVresCpu vresdata=Multires->GetGeomInfoVresCpu();

  acuint fspart("-",Arrays_Cpu,true);
  fvres::CallComputeFSNormals(Np,Npb,CSP,DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,Velrho_c->cptr(),FSType_c->ptr(),FSNormal_c->ptr(),fspart.ptr(),vresdata);

}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
//==============================================================================
void JSphCpuSingle_VRes::ComputeUmbrellaRegionVRes(){
  StrGeomVresCpu vresdata=Multires->GetGeomInfoVresCpu();

  acuint fspart("-",Arrays_Cpu,true);
  fvres::CallScanUmbrellaRegion(Np,Npb,CSP,DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,FSNormal_c->cptr(),fspart.ptr(),FSType_c->ptr(),vresdata);
}


//==============================================================================
/// Interaction for force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphCpuSingle_VRes::Interaction_ForcesB(TpInterStep interstep){
  //-Boundary correction for mDBC.
  const bool runmdbc=(TBoundary==BC_MDBC && (MdbcCorrector || interstep!=INTERSTEP_SymCorrector));
  const bool mdbc2=(runmdbc && SlipMode>=SLIP_NoSlip); //<vs_m2dbc>
  

  tfloat3* dengradcorr=NULL;

  Timersc->TmStart(TMC_CfForces);
  //-Interaction of Fluid-Fluid/Bound & Bound-Fluid (forces and DEM). | Interaccion Fluid-Fluid/Bound & Bound-Fluid (forces and DEM).
  const stinterparmsc parms=StInterparmsc(Np,Npb,NpbOk
    ,DivData,Dcell_c->cptr()
    ,Pos_c->cptr(),Velrho_c->cptr(),Idp_c->cptr(),Code_c->cptr(),Press_c->cptr()
    ,AC_CPTR(BoundMode_c),AC_CPTR(TangenVel_c),AC_CPTR(MotionVel_c) //<vs_m2dbc>
    ,AC_CPTR(BoundNor_c),AC_PTR(NoPenShift_c) //<vs_m2dbcNP>
    ,dengradcorr
    ,Ar_c->ptr(),Ace_c->ptr(),AC_PTR(Delta_c)
    ,ShiftingMode,AC_PTR(ShiftPosfs_c)
    ,AC_PTR(SpsTauRho2_c),AC_PTR(Sps2Strain_c)
    ,AC_PTR(FSType_c),AC_PTR(ShiftVel_c) //<vs_advshift>
  );
  StInterResultc res;
  res.viscdt=0;
  JSphCpu::Interaction_Forces_ct(parms,res);

  //-For 2-D simulations zero the 2nd component. | Para simulaciones 2D anula siempre la 2nd componente.
  if(Simulate2D){
    tfloat3* acec=Ace_c->ptr();
    const int ini=int(Npb),fin=int(Np),npf=int(Np-Npb);
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTELIGHT)
    #endif
    for(int p=ini;p<fin;p++)acec[p].y=0;
  }

  //-Add Delta-SPH correction to Ar_c[].
  if(AC_CPTR(Delta_c)){
    const float* deltac=Delta_c->cptr();
    float* arc=Ar_c->ptr();
    const int ini=int(Npb),fin=int(Np),npf=int(Np-Npb);
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTELIGHT)
    #endif
    for(int p=ini;p<fin;p++)if(deltac[p]!=FLT_MAX)arc[p]+=deltac[p];
  }

  //-Calculates maximum value of ViscDt.
  ViscDtMax=res.viscdt;
  //-Calculates maximum value of Ace (periodic particles are ignored).
  AceMax=ComputeAceMax();

  InterNum++;
  Timersc->TmStop(TMC_CfForces);
}



JMatrix4d JSphCpuSingle_VRes::CalcVelMotion(unsigned trackingmk,double dt){
  //Check if the mkbound is in one of the motion objects
  JMatrix4d mat;
  if(CaseNmoving){
     const unsigned nref=DsMotion->GetNumObjects();
      for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      const unsigned mk=m.mkbound;
      if(trackingmk==mk) return(CalcMotionMoving(m,dt));
      }
  }
  //Check if the mkbound is in one of the floating objects
  if(CaseNfloat){
    for(unsigned cf=0;cf<FtCount;cf++){
      //-Get Floating object values.
      const StFloatingData fobj=FtObjs[cf];
      const unsigned mk=fobj.mkbound;
      if(trackingmk==mk) return(CalcMotionFloating(fobj,dt));
  }
  }

  return (mat);

	

}



JMatrix4d JSphCpuSingle_VRes::CalcMotionMoving(const StMotionData m,double dt){
  JMatrix4d mat;
  if(m.type==MOTT_Linear){
    mat.Move(TDouble3(dt*m.linvel.x,dt*m.linvel.y,dt*m.linvel.z));
    mat.IsMovMatrix();
  } else if (m.type==MOTT_Matrix) {
    mat=m.matmov;
  }
  return(mat);
  
}

JMatrix4d JSphCpuSingle_VRes::CalcMotionFloating(const StFloatingData m,double dt){
  JMatrix4d mat;
  const tfloat3  fvel   =m.fvel;
  const tfloat3  fomega =m.fomega;
  const tdouble3 fcenter=m.center;
  const tdouble3 dang=(ToTDouble3(fomega)*dt)*TODEG;
  // const tdouble3 cen2=(PeriActive? UpdatePeriodicPos(fcenter): fcenter);

  mat.Move(fcenter);
  mat.Rotate(dang);
  mat.Move(TDouble3(-m.center.x,-m.center.y,-m.center.z));
  mat.Move(TDouble3(dt*fvel.x,dt*fvel.y,dt*fvel.z)); 
  return(mat);
  
}