#include "JSphGpuSingle_Mr.h"
#include "JSphGpuSingle.h"
#include "JSphBuffer.h"
#include "JSphGpu_Buffer_iker.h"
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
#include "JSphBuffer.h"
#include "JSphGpu_preloop_iker.h"



using namespace std;
JSphGpuSingle_Mr::JSphGpuSingle_Mr():JSphGpuSingle(){
  ClassName="JSphGpuSingle";
  MRfastsingle=true;
  MROrder=0;
  MRThreshold=100;
  Multires=NULL;


}



JSphGpuSingle_Mr::~JSphGpuSingle_Mr(){
  DestructorActive=true;
  delete Multires; Multires=NULL;
}

//==============================================================================
/// ???
//==============================================================================
void JSphGpuSingle_Mr::Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
  ,unsigned vrescount,unsigned vresid)
{
  if(!cfg || !log)return;
  AppName=appname; Log=log; CfgRun=cfg;
  VResCount=vrescount;
  VResId=vresid; 

  //-Selection of GPU.
  const int gpuid=SelecDevice(cfg->GpuId);


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
  CollectCTEdata();

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

double JSphGpuSingle_Mr::Init2(){
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


double JSphGpuSingle_Mr::ComputeStep_SymB(){

  cusph::CteInteractionUp(&CTE);
  const double dt=SymplecticDtPre;
  if(CaseNmoving)CalcMotion(dt);               //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
   DemDtForce=dt*0.5f;                          //-For DEM interaction.
  MdbcBoundCorrection(INTERSTEP_SymPredictor);  //-Mdbc correction
  PreInteraction_Forces(INTERSTEP_SymPredictor);//-Allocating temporary arrays.
  PreLoopProcedureVRes(INTERSTEP_SymPredictor); //-Calculate variables for interaction forces (Shifting,DDT,etc...).  
  Interaction_ForcesB(INTERSTEP_SymPredictor);  //-Interaction.
  const double dt_p=DtVariable(false);          //-Calculate dt of predictor step.
  if(Shifting)RunShifting(dt*.5);               //-Shifting.
  ComputeSymplecticPre(dt);                     //-Apply Symplectic-Predictor to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt*.5,true);        //-Control of floating bodies.
  PosInteraction_Forces();                      //-Free memory used for interaction.
  //-Corrector
  //-----------
  DemDtForce=dt;                                //-For DEM interaction.
  RunCellDivide(true);
  MdbcBoundCorrection(INTERSTEP_SymCorrector);  //-Mdbc correction
  PreInteraction_Forces(INTERSTEP_SymCorrector);//-Allocating temporary arrays.
  PreLoopProcedureVRes(INTERSTEP_SymCorrector); //-Calculate variables for interaction forces (Shifting,DDT,etc...).
  Interaction_ForcesB(INTERSTEP_SymCorrector);  //-Interaction.
  const double dt_c=DtVariable(true);           //-Calculate dt of corrector step.
  if(Shifting)RunShifting(dt);                  //-Shifting.
  ComputeSymplecticCorr(dt);                    //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);          //-Control of floating bodies.
  PosInteraction_Forces();                      //-Free memory used for interaction.
  if(Damping)RunDamping(dt);                    //-Aplies Damping.
  if(RelaxZones)RunRelaxZone(dt);               //-Generate waves using RZ.

  SymplecticDtPre1=min(dt_p,dt_c);            //-Calculate dt for next ComputeStep.
  return(dt);
}


void JSphGpuSingle_Mr::Finish(double dt1){
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
        // MultiRes->SaveFluxesData(DirDataOut,Part);
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

void JSphGpuSingle_Mr::Finish2(){
	TimerSim.Stop(); TimerTot.Stop();

	  //-End of Simulation.
	  //--------------------
	  FinishRun(false);

}

void JSphGpuSingle_Mr::InitMultires(const JSphCfgRun* cfg, JCaseVRes casemultires, unsigned id){
  // unsigned numzone=casemultires.Count();

  // const JCaseVRes_Box* Zone= casemultires.GetZoneBox(id);

  // unsigned numsubzone=Zone->Count();

  // tdouble3 innerdomainmin= Zone->GetBuffMax();
  // tdouble3 innerdomainmax=Zone->GetBuffMin();


  // JXml xml; xml.LoadFile(FileXml);
   Multires=new JSphBuffer(Cpu,CSP,casemultires,id,AppName,DirDataOut,PartBegin,PartBeginDir);
   Multires->Config();
}

StInterParmsbg JSphGpuSingle_Mr::getParms(){
	StrInterParmsbg parms=StrInterParmsbg(Simulate2D,TKernel
		    ,DivData,CTE,Map_PosMin,Dcell_g->cptr()
		    ,Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr()
        ,Velrho_g->cptr(),Idp_g->cptr(),Code_g->cptr()
		    ,NULL,NULL);
 return(parms);
}


void JSphGpuSingle_Mr::CollectCTEdata(){
  GetConstantData(CTE);
}

void JSphGpuSingle_Mr::BufferInit(StInterParmsbg *parms){
    cusph::CteInteractionUp(&CTE);

	for(unsigned i=0;i<Multires->GetCount();i++){
		int *bufferpart=NULL;
		cudaMalloc((void**)&bufferpart, sizeof(int)*(Np+1));
		unsigned buffercountpre=0;

//		DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg);
		buffercountpre=Multires->CreateListGpuInit(Np,0,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize,bufferpart,i);
		// std::cout << buffercountpre << std::endl;
		cudaFree(bufferpart);

	}

  RunCellDivide(true);
    if(VRES_DG_SAVEVTK)DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->cptr());

}

void JSphGpuSingle_Mr::BufferExtrapolateData(StInterParmsbg *parms){

	for(unsigned i=0;i<Multires->GetCount();i++){
		
    agint bufferpartg("-",Arrays_Gpu,true);
		unsigned buffercountpre;

		buffercountpre=Multires->CreateListGpu(Np-Npb,Npb,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize,bufferpartg.ptr(),i);

		
		unsigned id=Multires->GetZone(i)->getZone()-1;


		cusph::CteInteractionUp(&parms[id].cte);
		cusphbuffer::Interaction_BufferExtrap(buffercountpre,bufferpartg.ptr(),parms[id],Posxy_g->cptr(),Posz_g->cptr(),Velrho_g->ptr(),NULL,Code_g->ptr(),true,1,100);
		// if(Nstep==0)  DgSaveVtkParticlesGpuMultiRes("Compute_Init_",Nstep,0,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg,rcond);

		double t1=TimeStep+LastDt;
		double t2=TimePartNext;
		if(t1>=t2 ){
  // if(VRES_DG_SAVEVTK)DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->cptr());
		}


		
	}
	// cudaFree(rcond);
}

void JSphGpuSingle_Mr::ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,StInterParmsbg *parms){
	// TmgStart(Timers,TMG_SuBuffer);
    Multires->UpdateMatMov(mat);
	  Multires->MoveBufferZone(dt,mat);
	// TmgStart(Timers,TMG_SuBuffer);

	// Multires->MoveBufferZone(dt,velmot);
	
	for(unsigned i=0;i<Multires->GetCount();i++){
    
    cusph::CteInteractionUp(&CTE);

		agint bufferpartg("-",Arrays_Gpu,true);
		
    
    unsigned buffercountpre=0;

		buffercountpre=Multires->CreateListGpu(Np-Npb,Npb,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize,bufferpartg.ptr(),i);

    StrDataVresGpu vresdata=Multires->GetZoneFluxInfoGpu(i);
			
    // cusphbuffer::MoveBufferZone(nini,ntot,posxy,posz,dt,velmot[i]);


		unsigned id=Multires->GetZone(i)->getZone()-1;


		cusph::CteInteractionUp(&parms[id].cte);
		
    cusphbuffer::Interaction_BufferExtrapFlux(vresdata.ntot,vresdata.nini,parms[id],vresdata.ptposxy,vresdata.ptposz,vresdata.normals,vresdata.mass,CTE.dp,dt,vresdata.velmot,true,1,100);

		Multires->ComputeStepGpu(buffercountpre,bufferpartg.ptr(),IdMax+1,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr()
        ,Code_g->ptr(),Idp_g->ptr(),Velrho_g->ptr(),NULL,this,i);

		
    
    int *newpart=NULL;
			cudaMalloc((void**)&newpart, sizeof(int)*(vresdata.ntot+1));
			cudaMemset(newpart,0,vresdata.ntot*sizeof(int));

		cusph::CteInteractionUp(&CTE);
		unsigned newnp=cusphbuffer::BufferListCreate(false,vresdata.ntot,vresdata.nini,vresdata.ntot, vresdata.mass, newpart, CTE.massf);
		if(newnp){			
      if(!CheckGpuParticlesSize(Np+newnp)){
        const unsigned ndatacpu=0,ndatagpu=Np;
        ResizeParticlesSizeData(ndatacpu,ndatagpu,Np+newnp,Np+newnp, 0.2,true);
        CellDivSingle->SetIncreaseNp(newnp);
			}
      
      cusphbuffer::BufferCreateNewPart(PeriActive,newnp,vresdata.nini,newpart,Np,IdMax+1,vresdata.normals,CTE.dp,Posxy_g->ptr(),Posz_g->ptr()
          ,vresdata.ptposxy,vresdata.ptposz,Dcell_g->ptr(),Code_g->ptr(),Idp_g->ptr(),Velrho_g->ptr(),i);
    }
      
        
    if(newnp){      
      if(SpsTauRho2_g)SpsTauRho2_g->CuMemsetOffset(Np,0,newnp);
      if(BoundNor_g)BoundNor_g->CuMemsetOffset(Np,0,newnp);
      // if(MotionAce_g)MotionAce_g->CuMemsetOffset(Np,0,newnp);
      if(FSType_g)FSType_g->CuMemsetOffset(Np,0,newnp);
      if(ShiftVel_g)ShiftVel_g->CuMemsetOffset(Np,0,newnp);
      Np+=newnp; 
      TotalNp+=newnp;
      IdMax=unsigned(TotalNp-1);
		} 
    cudaFree(newpart);
          newpart=NULL;

	}
  // if(VRES_DG_SAVEVTK)DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->cptr());
}

void JSphGpuSingle_Mr::BufferShifting(){
  StrGeomVresGpu* vresgdata=Multires->GetGeomInfoVres();
	cusphbuffer::BufferShiftingGpu(Np,Npb,Posxy_g->ptr(),Posz_g->ptr(),ShiftVel_g->ptr(),Code_g->ptr(),vresgdata,NULL);
}



void JSphGpuSingle_Mr::CallRunCellDivide(){
  	cusph::CteInteractionUp(&CTE);
    RunCellDivide(true);
}


//==============================================================================
/// PreLoop for additional models computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle_Mr::PreLoopProcedureVRes(TpInterStep interstep){
  bool runshift=(interstep==INTERSTEP_SymPredictor && Nstep!=0 && ShiftingAdv!=NULL);
  if(runshift){
    ComputeFSParticlesVRes();
    ComputeUmbrellaRegionVRes();
  }
  
  unsigned bsfluid=BlockSizes.forcesfluid;
  StrGeomVresGpu* vresgdata=Multires->GetGeomInfoVres();

  if(runshift)cusphbuffer::PreLoopInteraction(TKernel,Simulate2D,runshift,false,bsfluid,Np-Npb,Npb,DivData
    ,Posxy_g->cptr(),Posz_g->cptr(),Dcell_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr(),Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr()
    ,FSType_g->ptr(),FSNormal_g->ptr(),FSMinDist_g->ptr(),vresgdata,NULL);
  
  if(runshift)cusph::ComputeShiftingVel(bsfluid,Np-Npb,Npb,Simulate2D,ShiftVel_g->ptr(),FSType_g->cptr(),FSNormal_g->cptr()
    ,FSMinDist_g->cptr(),SymplecticDtPre,ShiftingAdv->GetShiftCoef(),ShiftingAdv->GetAleActive(),NULL);

  if(runshift )BufferShifting();
  //-Updates preloop variables in periodic particles.
  if(PeriParent_g){
    cusph::PeriPreLoopCorr(Np,0,PeriParent_g->cptr(),FSType_g->ptr(),ShiftVel_g->ptr());
  }

  double t1=TimeStep+LastDt;
	double t2=TimePartNext;
  if(t1>=t2 ){
		if(interstep==INTERSTEP_SymPredictor && ShiftingAdv!=NULL)DgSaveVtkParticlesGpu("Compute_FreeSurface_",Part,0,Np,Posxy_g->cptr()
        ,Posz_g->cptr(),Code_g->cptr(),FSType_g->cptr(),ShiftVel_g->cptr(),FSNormal_g->cptr());
	}
  
}


//==============================================================================
/// Compute free-surface particles and their normals.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle_Mr::ComputeFSParticlesVRes(){
  unsigned bsfluid=BlockSizes.forcesfluid;
  StrGeomVresGpu* vresgdata=Multires->GetGeomInfoVres();

  aguint    inoutpartg("-",Arrays_Gpu,true);
  cusphbuffer::ComputeFSNormals(TKernel,Simulate2D,Symmetry,bsfluid,Npb,Np-Npb,DivData
    ,Dcell_g->cptr(),Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr()
    ,Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr(),FSType_g->ptr(),FSNormal_g->ptr()
    ,inoutpartg.ptr(),vresgdata,NULL);

}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle_Mr::ComputeUmbrellaRegionVRes(){
  unsigned bsfluid=BlockSizes.forcesfluid;
  StrGeomVresGpu* vresgdata=Multires->GetGeomInfoVres();

  aguint    inoutpartg("-",Arrays_Gpu,true);
  cusphbuffer::ComputeUmbrellaRegion(TKernel,Simulate2D,Symmetry,bsfluid,Npb,Np-Npb,DivData
    ,Dcell_g->cptr(),Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr()
    ,Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr(),FSType_g->ptr(),FSNormal_g->ptr()
    ,inoutpartg.ptr(),vresgdata,NULL);

}


//==============================================================================
/// Interaction for force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle_Mr::Interaction_ForcesB(TpInterStep interstep){
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
    ,Symmetry //<vs_syymmetry>
    ,TKernel,FtMode
    ,TVisco,TDensity,ShiftingMode,mdbc2 //<vs_m2dbc>
    ,shiftadv,corrector,aleform,ncpress         //<ShiftingAdvanced>
    ,Visco*ViscoBoundFactor,Visco
    ,bsbound,bsfluid,Np,Npb,NpbOk
    ,0,Nstep,DivData,Dcell_g->cptr()
    ,Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr()
    ,Velrho_g->cptr(),Idp_g->cptr(),Code_g->cptr()
    ,AG_CPTR(BoundMode_g),AG_CPTR(TangenVel_g),AG_CPTR(MotionVel_g) //<vs_m2dbc>
    ,FtoMasspg,AG_CPTR(SpsTauRho2_g),dengradcorr
    ,ViscDt_g->ptr(),Ar_g->ptr(),Ace_g->ptr(),AG_PTR(Delta_g)
    ,AG_PTR(Sps2Strain_g)
    ,AG_PTR(ShiftPosfs_g)
    ,AG_PTR(FSType_g),AG_CPTR(ShiftVel_g)           //<ShiftingAdvanced>
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



JMatrix4d JSphGpuSingle_Mr::CalcVelMotion(unsigned trackingmk,double dt){
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



JMatrix4d JSphGpuSingle_Mr::CalcMotionMoving(const StMotionData m,double dt){
  JMatrix4d mat;
  if(m.type==MOTT_Linear){
    mat.Move(TDouble3(dt*m.linvel.x,dt*m.linvel.y,dt*m.linvel.z));
    mat.IsMovMatrix();
  } else if (m.type==MOTT_Matrix) {
    mat=m.matmov;
  }
  return(mat);
  
}

JMatrix4d JSphGpuSingle_Mr::CalcMotionFloating(const StFloatingData m,double dt){
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