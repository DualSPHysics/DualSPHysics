#include "JSphGpuSingle_VRes.h"
#include "JSphGpuSingle.h"
#include "JSphVRes.h"
#include "JSphGpu_VRes_iker.h"
#include "JTimeControl.h"
#include "JDsOutputTime.h"
#include "JCellDivGpuSingle.h"
#include "JSphInOut.h"
#include "JXml.h"
#include "JDsMotion.h"
#include "JDsPartMotionSave.h"
#include "JDsPartFloatSave.h"
#include "JSphShifting.h"
#include "JSphShiftingAdv.h"
#include "JSphGpu_preloop_iker.h"


//==============================================================================
/// Constructor.
//==============================================================================
using namespace std;
JSphGpuSingle_VRes::JSphGpuSingle_VRes():JSphGpuSingle(){
  ClassName="JSphGpuSingle";
  MRfastsingle=false;
  MROrder=0;
  MRThreshold=100;
  Multires=NULL;
  VResOrder=VrOrder_1st;


}


//==============================================================================
/// Destructor.
//==============================================================================
JSphGpuSingle_VRes::~JSphGpuSingle_VRes(){
  DestructorActive=true;
  delete Multires; Multires=NULL;
}


//==============================================================================
/// Load VRes configuration.
//==============================================================================
void JSphGpuSingle_VRes::LoadVResConfigParameters(const JSphCfgRun* cfg){
  if(cfg->MROrder>=0){
    switch(cfg->MROrder){
      case 0:   VResOrder=VrOrder_0th;
      case 1:   VResOrder=VrOrder_1st;
      case 2:   VResOrder=VrOrder_2nd;
      default:  Run_Exceptioon("Variable resolution reconstruction method is not valid.");
    }
  }
}

//==============================================================================
/// Initialises VRes simulation.
//==============================================================================
void JSphGpuSingle_VRes::Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
  ,unsigned vrescount,unsigned vresid)
{
  if(!cfg || !log)return;
  AppName=appname; Log=log; CfgRun=cfg;
  VResCount=vrescount;
  VResId=vresid; 

  //-Selection of GPU.
  const int gpuid=SelecDevice(cfg->GpuId);
  LoadVResConfigParameters(cfg);


  if(cfg->MRFastSingle>=0)MRfastsingle=(cfg->MRFastSingle>0);
  if(cfg->MROrder>=0)MROrder=cfg->MROrder;
  if(cfg->MRThreshold >=0)MRThreshold=cfg->MRThreshold;


  //-Creates array system for particles.
  Arrays_Cpu=new JArraysCpu(Log);
  Arrays_Gpu=new JArraysGpu(gpuid,Log);

  //-Configures timers.
  Timersg->Config(cfg->SvTimers);
  Timersg->TmStart(TMG_Init,false);

  //-Load parameters and input data.
  LoadConfig(cfg);
  LoadCaseParticles();
  VisuConfig();
  ConfigDomain();
  ConfigRunMode();
  VisuParticleSummary();
  
  //-Store Interaction parameters for vres.
  GetConstantData(CTE);


  //-Initialisation of execution variables.
  InitRunGpu();
  RunFirstGaugeSystem(TimeStep);
  if(InOut)InOutInit(TimeStepIni);
  FreePartsInit();
  PrintAllocMemory(GetAllocMemoryCpu(),GetAllocMemoryGpu());
  UpdateMaxValues();
  SaveData(); 
  Arrays_Cpu->FreeUnusedArrays();
  Arrays_Gpu->FreeUnusedArrays();
  Timersg->ResetTimes();
  Timersg->TmStop(TMG_Init,true);
  if(Log->WarningCount())Log->PrintWarningList("\n[WARNINGS]","");
  PartNstep=-1; Part++;
}

//==============================================================================
/// Complete initialization of VRes simulation.
//==============================================================================
double JSphGpuSingle_VRes::Init2(){
  cusph::CteInteractionUp(&CTE);

	//-Main Loop.
	//------------
	JTimeControl tc("30,60,300,600");//-Shows information at 0.5, 1, 5 y 10 minutes (before first PART).
	bool partoutstop=false;
	TimerSim.Start();
	TimerPart.Start();
	Log->Print(string("\n[Initialising simulation (")+RunCode+")  "+fun::GetDateTime()+"]");
	if(DsPips)ComputePips(true);
	PrintHeadPart();
	return(TimeMax);
}

//==============================================================================
/// Compute time step of VRes simulation.
//==============================================================================
double JSphGpuSingle_VRes::ComputeStep_SymB(){

  cusph::CteInteractionUp(&CTE);
  const double dt=SymplecticDtPre;
  if(CaseNmoving)CalcMotion(dt);               //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
  InterStep=INTERSTEP_SymPredictor;
  DemDtForce=dt*0.5f;                           //-For DEM interaction.
  MdbcBoundCorrection(InterStep);               //-Mdbc correction
  PreInteraction_Forces(InterStep);             //-Allocating temporary arrays.
  PreLoopProcedureVRes(InterStep);              //-Calculate variables for interaction forces (Shifting,DDT,etc...).  
  Interaction_ForcesB(InterStep);               //-Interaction.
  const double dt_p=DtVariable(false);          //-Calculate dt of predictor step.
  if(Shifting)RunShifting(dt*.5);               //-Shifting.
  ComputeSymplecticPre(dt);                     //-Apply Symplectic-Predictor to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt*.5,true);        //-Control of floating bodies.
  PosInteraction_Forces();                      //-Free memory used for interaction.
  //-Corrector
  //-----------
  InterStep=INTERSTEP_SymCorrector;
  DemDtForce=dt;                                //-For DEM interaction.
  RunCellDivide(true);
  MdbcBoundCorrection(InterStep);               //-Mdbc correction
  PreInteraction_Forces(InterStep);             //-Allocating temporary arrays.
  PreLoopProcedureVRes(InterStep);              //-Calculate variables for interaction forces (Shifting,DDT,etc...).
  Interaction_ForcesB(InterStep);               //-Interaction.
  const double dt_c=DtVariable(true);           //-Calculate dt of corrector step.
  if(Shifting)RunShifting(dt);                  //-Shifting.
  ComputeSymplecticCorr(dt);                    //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);          //-Control of floating bodies.
  PosInteraction_Forces();                      //-Free memory used for interaction.
  if(Damping)RunDamping(dt);                    //-Aplies Damping.
  if(RelaxZones)RunRelaxZone(dt);               //-Generate waves using RZ.

  SymplecticDtPre1=min(dt_p,dt_c);              //-Calculate dt for next ComputeStep.
  return(dt);
}

//==============================================================================
/// Complete time step of VRes simulation.
//==============================================================================
void JSphGpuSingle_VRes::Finish(double dt1){
	cusph::CteInteractionUp(&CTE);
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
    if(PartMotionSave)PartMotionSave->AddDataExtraGpu(Part,TimeStep,Nstep,Np
      ,Posxy_g->cptr(),Posz_g->cptr(),RidpMotg);
    TimeOutExtraUpdate(TimeStep);
  }
	if(TimeStep>=TimePartNext ){
    if(VRES_DG_SAVEVTK)DgSaveVtkParticlesGpuMultiRes("Compute_",Part,0,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->cptr(),NULL,NULL);
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

void JSphGpuSingle_VRes::Finish2(){
	TimerSim.Stop(); TimerTot.Stop();
	//-End of Simulation.
	//--------------------
	FinishRun(false);
}

void JSphGpuSingle_VRes::InitMultires(const JSphCfgRun* cfg, JCaseVRes casemultires, unsigned id){
   Multires=new JSphVRes(Cpu,CSP,casemultires,id,AppName,DirDataOut,PartBegin,PartBeginDir);
   Multires->Config();
}

StInterParmsbg JSphGpuSingle_VRes::getParms(){
	StrInterParmsbg parms=StrInterParmsbg(Simulate2D,TKernel
		    ,DivData,CTE,Map_PosMin,Dcell_g->cptr()
		    ,Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr()
        ,Velrho_g->cptr(),Idp_g->cptr(),Code_g->cptr()
		    ,NULL,NULL);
 return(parms);
}


void JSphGpuSingle_VRes::CollectCTEdata(){
  GetConstantData(CTE);
}

void JSphGpuSingle_VRes::BufferInit(StInterParmsbg *parms){
  cusph::CteInteractionUp(&CTE);


	for(unsigned i=0;i<Multires->GetCount();i++){
		agint bufferpartg("-",Arrays_Gpu,true);
		unsigned buffercountpre;

//		DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
		buffercountpre=Multires->CreateListGpuInit(Np,0,Posxy_g->cptr()
      ,Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize,bufferpartg.ptr(),i);
	}

  RunCellDivide(true);
  if(VRES_DG_SAVEVTK)DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->cptr());

}


void JSphGpuSingle_VRes::BufferExtrapolateData(StInterParmsbg *parms){
  const unsigned count=Multires->GetCount();
	for(unsigned i=0;i<count;i++){
    
    //-Compute list of buffer particles.	
    agint bufferpartg("-",Arrays_Gpu,true);
		const unsigned buffercountpre=Multires->CreateListGpu(Np-Npb,Npb
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize
      ,bufferpartg.ptr(),i);
	
    //-Update constant memory and perform interpolation.  
		unsigned id=Multires->GetZone(i)->getZone()-1;
		cusph::CteInteractionUp(&parms[id].cte);
		cusphvres::Interaction_BufferExtrap(buffercountpre,bufferpartg.ptr(),parms[id]
      ,Posxy_g->cptr(),Posz_g->cptr(),Velrho_g->ptr(),Code_g->ptr(),true,VrOrder_1st,100);
	}
}

//==============================================================================
/// ComputeStep over buffer regions:
/// - If buffer particle is moved to fluid zone then it changes to fluid particle.
/// - If fluid particle is moved to buffer zone then it changes to buffer particle.
/// - If buffer particle is moved out the domain then it changes to ignore particle.
/// - Create buffer particle based on eulerian flux on the accumulations points.
/// - Update position of buffer regions and compute motion velocity.
//==============================================================================
void JSphGpuSingle_VRes::ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,StInterParmsbg *parms){
    
  //-Compute movement for VRes regions
  Multires->UpdateMatMov(mat);
	Multires->MoveBufferZone(dt,mat);

	//- ComputeStep for each buffer region.
	for(unsigned i=0;i<Multires->GetCount();i++){
    
    cusph::CteInteractionUp(&CTE);                            //-Update constant memory.

    //-Compute list of buffer particles.
		agint bufferpartg("-",Arrays_Gpu,true);       
		const unsigned buffercountpre=Multires->CreateListGpu(Np-Npb,Npb
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize
      ,bufferpartg.ptr(),i);

    //-Compute eulerian flux on the mass accumulation points.
  	unsigned id=Multires->GetZone(i)->getZone()-1;      
		cusph::CteInteractionUp(&parms[id].cte);                  //-Update constant memory.

    StrDataVresGpu vresdata=Multires->GetZoneFluxInfoGpu(i);  //-Retrieve buffer parameters.			

		
    cusphvres::Interaction_BufferExtrapFlux(parms[id]
      ,vresdata,CTE.dp,dt,true,VrOrder_1st,100);

		
    Multires->ComputeStepGpu(buffercountpre,bufferpartg.ptr(),IdMax+1,Posxy_g->ptr()
      ,Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr(),Idp_g->ptr(),Velrho_g->ptr(),NULL,this,i);

		
    cusph::CteInteractionUp(&CTE);                            //-Update constant memory.
    
    int *newpart=NULL;
		cudaMalloc((void**)&newpart, sizeof(int)*(vresdata.ntot+1));
		cudaMemset(newpart,0,vresdata.ntot*sizeof(int));

		
    cusphvres::CheckMassFlux(vresdata.ntot,vresdata.nini,DivData,Map_PosMin
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),PosCell_g->cptr()
      ,vresdata.ptposxy,vresdata.ptposz,vresdata.normals,vresdata.mass);
		
    
    unsigned newnp=cusphvres::BufferListCreate(false,vresdata.ntot,vresdata.nini
      ,vresdata.ntot, vresdata.mass, newpart, CTE.massf);
		
    if(newnp){			
      if(!CheckGpuParticlesSize(Np+newnp)){
        const unsigned ndatacpu=0,ndatagpu=Np;
        ResizeParticlesSizeData(ndatacpu,ndatagpu,Np+newnp,Np+newnp, 0.2,true);
        CellDivSingle->SetIncreaseNp(newnp);
			}  
      cusphvres::BufferCreateNewPart(PeriActive,newnp,vresdata.nini,newpart,Np,IdMax+1,vresdata.normals,CTE.dp,Posxy_g->ptr(),Posz_g->ptr()
        ,vresdata.ptposxy,vresdata.ptposz,Dcell_g->ptr(),Code_g->ptr(),Idp_g->ptr(),Velrho_g->ptr(),i);
      
      //-Updates basic arrays.
      if(SpsTauRho2_g)SpsTauRho2_g->CuMemsetOffset(Np,0,newnp);
      if(BoundNor_g)BoundNor_g->CuMemsetOffset(Np,0,newnp);
      if(FSType_g)FSType_g->CuMemsetOffset(Np,0,newnp);
      if(ShiftVel_g)ShiftVel_g->CuMemsetOffset(Np,0,newnp);
      Np+=newnp; 
      TotalNp+=newnp;
      IdMax=unsigned(TotalNp-1);
    }      
    cudaFree(newpart);  newpart=NULL; 
	}
}

void JSphGpuSingle_VRes::BufferShifting(){
  StrGeomVresGpu vresgdata=Multires->GetGeomInfoVresGpu();
	cusphvres::BufferShiftingGpu(Np,Npb,Posxy_g->ptr(),Posz_g->ptr(),ShiftVel_g->ptr(),Code_g->ptr(),vresgdata,NULL);
}



void JSphGpuSingle_VRes::CallRunCellDivide(){
  	cusph::CteInteractionUp(&CTE);
    RunCellDivide(true);
}


//==============================================================================
/// PreLoop for additional models computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================

//==============================================================================
/// PreLoop for additional models computation.
//==============================================================================
void JSphGpuSingle_VRes::PreLoopProcedureVRes(TpInterStep interstep){
  const bool runshift=(ShiftingAdv && interstep==INTERSTEP_SymPredictor && Nstep!=0);
  if(runshift){
    Timersg->TmStart(TMG_SuShifting,false);
    ComputeFSParticlesVRes();
    ComputeUmbrellaRegionVRes();
    const unsigned bsfluid=BlockSizes.forcesfluid;
     StrGeomVresGpu vresgdata=Multires->GetGeomInfoVresGpu();

  if(runshift)cusphvres::PreLoopInteraction(TKernel,Simulate2D,runshift,bsfluid,Np-Npb,Npb,DivData
    ,Posxy_g->cptr(),Posz_g->cptr(),Dcell_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr(),Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr()
    ,FSType_g->ptr(),FSNormal_g->ptr(),FSMinDist_g->ptr(),vresgdata,NULL);

      
    cusph::ComputeShiftingVel(bsfluid,Np-Npb,Npb,Simulate2D,ShiftingAdv->GetShiftCoef()
      ,ShiftingAdv->GetAleActive(),float(SymplecticDtPre),FSType_g->cptr()
      ,FSNormal_g->cptr(),FSMinDist_g->cptr(),ShiftVel_g->ptr(),NULL);
    //-Updates pre-loop variables in periodic particles.
    BufferShifting();
    if(PeriParent_g){
      cusph::PeriPreLoopCorr(Np,0,PeriParent_g->cptr(),FSType_g->ptr()
        ,ShiftVel_g->ptr());
    }

    //-Saves VTK for debug.
    if(1 && runshift && TimeStep+LastDt>=TimePartNext){
		  DgSaveVtkParticlesGpu("Compute_FreeSurface_",Part,0,Np,Posxy_g->cptr()
        ,Posz_g->cptr(),Code_g->cptr(),FSType_g->cptr(),ShiftVel_g->cptr()
        ,FSNormal_g->cptr());
	  }
    
    Timersg->TmStop(TMG_SuShifting,true);
  }
}


//==============================================================================
/// Compute free-surface particles and their normals.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle_VRes::ComputeFSParticlesVRes(){
  unsigned bsfluid=BlockSizes.forcesfluid;
  StrGeomVresGpu vresgdata=Multires->GetGeomInfoVresGpu();

  aguint    inoutpartg("-",Arrays_Gpu,true);
  cusphvres::ComputeFSNormals(TKernel,Simulate2D,bsfluid,Npb,Np-Npb,DivData
    ,Dcell_g->cptr(),Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr()
    ,Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr(),FSType_g->ptr(),FSNormal_g->ptr()
    ,inoutpartg.ptr(),vresgdata,NULL);

}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle_VRes::ComputeUmbrellaRegionVRes(){
  unsigned bsfluid=BlockSizes.forcesfluid;
  StrGeomVresGpu vresgdata=Multires->GetGeomInfoVresGpu();

  aguint    inoutpartg("-",Arrays_Gpu,true);
  cusphvres::ComputeUmbrellaRegion(TKernel,Simulate2D,bsfluid,Npb,Np-Npb,DivData
    ,Dcell_g->cptr(),Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr()
    ,Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr(),FSType_g->ptr(),FSNormal_g->ptr()
    ,inoutpartg.ptr(),vresgdata,NULL);

}


//==============================================================================
/// Interaction for force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle_VRes::Interaction_ForcesB(TpInterStep interstep){
  const bool runmdbc=(TBoundary==BC_MDBC && (MdbcCorrector || interstep!=INTERSTEP_SymCorrector));
  const bool mdbc2=(runmdbc && SlipMode>=SLIP_NoSlip);
  InterStep=interstep;

  float3* dengradcorr=NULL;

  Timersg->TmStart(TMG_CfForces,true);
  unsigned bsfluid=BlockSizes.forcesfluid;
  unsigned bsbound=BlockSizes.forcesbound;

  //<ShiftingAdvanced>
  bool shiftadv=ShiftingAdv!=NULL;
  bool corrector= InterStep==INTERSTEP_SymCorrector;
  bool aleform=shiftadv ? ShiftingAdv->GetAleActive() : false;
  bool ncpress=shiftadv ? ShiftingAdv->GetNcPress()   : false;
  // if(interstep==INTERSTEP_SymPredictor && ShiftingAdv!=NULL)DgSaveVtkParticlesGpu("Compute_FreeSurface_",Nstep,0,Np,Posxy_g->cptr()
  //       ,Posz_g->cptr(),Code_g->cptr(),FSType_g->cptr(),ShiftVel_g->cptr(),FSNormal_g->cptr());

  //-Interaction Fluid-Fluid/Bound & Bound-Fluid.
  const StInterParmsg parms=StrInterParmsg(Simulate2D
    ,TKernel,FtMode
    ,TVisco,TDensity,ShiftingMode,TMdbc2 //<vs_m2dbcNP>
    ,shiftadv,corrector,aleform,ncpress //<vs_advshift>
    ,Visco*ViscoBoundFactor,Visco
    ,bsbound,bsfluid,Np,Npb,NpbOk
    ,0,Nstep,DivData,Dcell_g->cptr()
    ,Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr()
    ,Velrho_g->cptr(),Idp_g->cptr(),Code_g->cptr()
    ,AG_CPTR(BoundMode_g),AG_CPTR(TangenVel_g),AG_CPTR(MotionVel_g)//<vs_m2dbc>
    ,AG_CPTR(BoundNor_g),AG_PTR(NoPenShift_g) //<vs_m2dbcNP>
    ,FtoMasspg,AG_CPTR(SpsTauRho2_g),dengradcorr
    ,ViscDt_g->ptr(),Ar_g->ptr(),Ace_g->ptr(),AG_PTR(Delta_g)
    ,AG_PTR(Sps2Strain_g)
    ,AG_PTR(ShiftPosfs_g)
    ,AG_PTR(FSType_g),AG_CPTR(ShiftVel_g) //<vs_advshift>
    ,NULL,NULL);
  cusph::Interaction_Forces(parms);


  //-Interaction DEM Floating-Bound & Floating-Floating. //(DEM)
  if(UseDEM)cusph::Interaction_ForcesDem(BlockSizes.forcesdem,CaseNfloat
    ,DivData,Dcell_g->cptr(),RidpMotg+CaseNmoving,DemDatag,FtoMasspg
    ,float(DemDtForce),PosCell_g->cptr(),Velrho_g->cptr(),Code_g->cptr()
    ,Idp_g->cptr(),ViscDt_g->ptr(),Ace_g->ptr(),NULL);

  //-For 2D simulations always overrides the 2nd component (Y axis).
  //-Para simulaciones 2D anula siempre la 2nd componente.
  if(Simulate2D)cusph::Resety(Np-Npb,Npb,Ace_g->ptr());

  //-Computes Tau for Laminar+SPS.
  if(TVisco==VISCO_LaminarSPS)cusph::ComputeSpsTau(Np,Npb,SpsSmag,SpsBlin
    ,Velrho_g->cptr(),Sps2Strain_g->cptr(),SpsTauRho2_g->ptr());
  
  //-Add Delta-SPH correction to Ar_g[].
  if(AG_CPTR(Delta_g))cusph::AddDelta(Np-Npb,Delta_g->cptr()+Npb,Ar_g->ptr()+Npb);

  cudaDeviceSynchronize();
  Check_CudaErroor("Failed while executing kernels of interaction.");

  //-Calculates maximum value of ViscDt.
  if(Np)ViscDtMax=cusph::ReduMaxFloat(Np,0,ViscDt_g->ptr(),CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np)));
  //-Calculates maximum value of Ace (periodic particles are ignored). ViscDtg is used like auxiliary memory.
  AceMax=ComputeAceMax(ViscDt_g->ptr()); 

  InterNum++;
  Timersg->TmStop(TMG_CfForces,true);
  Check_CudaErroor("Failed in reduction of viscdt.");
}



JMatrix4d JSphGpuSingle_VRes::CalcVelMotion(unsigned trackingmk,double dt){
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



JMatrix4d JSphGpuSingle_VRes::CalcMotionMoving(const StMotionData m,double dt){
  JMatrix4d mat;
  if(m.type==MOTT_Linear){
    mat.Move(TDouble3(dt*m.linvel.x,dt*m.linvel.y,dt*m.linvel.z));
    mat.IsMovMatrix();
  } else if (m.type==MOTT_Matrix) {
    mat=m.matmov;
  }
  return(mat);
  
}

JMatrix4d JSphGpuSingle_VRes::CalcMotionFloating(const StFloatingData m,double dt){
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