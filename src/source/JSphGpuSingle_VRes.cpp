//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpuSingle_VRes.cpp \brief Implements the class \ref JSphGpuSingle_VRes.

#include "JSphGpuSingle_VRes.h"
#include "JSphGpuSingle.h"
#include "JSphVRes.h"
#include "JSphGpu.h"
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

using namespace std;
//==============================================================================
/// Constructor.
//==============================================================================
JSphGpuSingle_VRes::JSphGpuSingle_VRes():JSphGpuSingle(){
  ClassName="JSphGpuSingle";
  VResFastSingle=true;
  VResThreshold=100;
  VRes=NULL;
  VResOrder=VrOrder_1st;
  VResMethod=VrMethod_Liu;
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphGpuSingle_VRes::~JSphGpuSingle_VRes(){
  DestructorActive=true;
  delete VRes; VRes=NULL;
}

//==============================================================================
/// Load VRes configuration.
//==============================================================================
void JSphGpuSingle_VRes::LoadVResConfigParameters(const JSphCfgRun* cfg){
  if(cfg->VResOrder>=0){
    switch(cfg->VResOrder){
      case 0:   VResOrder=VrOrder_0th;  break;
      case 1:   VResOrder=VrOrder_1st;  break;
      case 2:   VResOrder=VrOrder_2nd;  break;
      default:  Run_Exceptioon("Variable resolution reconstruction order is not valid.");
    }
  }
  if(cfg->VResMethod>=0){
    switch(cfg->VResMethod){
      case 0:   VResMethod=VrMethod_Liu;  break;
      case 1:   VResMethod=VrMethod_Mls;  break;
      default:  Run_Exceptioon("Variable resolution reconstruction method is not valid.");
    }
  }
  if(cfg->MRFastSingle>=0)  VResFastSingle=cfg->MRFastSingle>0;
  if(cfg->VResThreshold>=0) VResThreshold=cfg->VResThreshold;
}

//==============================================================================
/// Print VRes configuration.
//==============================================================================
void JSphGpuSingle_VRes::VisuConfigVRes(){
  Log->Print("\nVRes Configuration:");
  Log->Printf(" Interpolation Method: %s",(VResMethod==VrMethod_Liu? "Liu-Liu Correction": "Moving Least Square"));
  Log->Printf(" Interpolation Order: %s" ,(VResOrder==VrOrder_2nd? "2nd": (VResOrder==VrOrder_1st? "1st": "0th")));
  Log->Printf(" Inter Threshold: %g" ,VResThreshold);
  Log->Printf(" Inter Precision: %s" ,(VResFastSingle ? "Single": "Double"));
}

//==============================================================================
/// Initialize VRes object.
//==============================================================================
void JSphGpuSingle_VRes::VResInit(const JSphCfgRun* cfg,JCaseVRes casemultires,unsigned id){
   VRes=new JSphVRes(Cpu,CSP,casemultires,id,MapRealPosMin
    ,MapRealPosMax,AppName,DirDataOut,PartBegin,PartBeginDir);
   VRes->Config();
}

//==============================================================================
/// Add warning for VRes execution.
//==============================================================================
void JSphGpuSingle_VRes::AddWarningVRes(){
   if(Moorings)Log->PrintWarning("The use of MoorDynPlus with Variable resolution is an experimental feature.");
   if(UseChrono)Log->PrintWarning("The use of CHRONO with Variable resolution is an experimental feature.");
   if(FlexStruc)Log->PrintWarning("The use of FlexStruct with Variable resolution is an experimental feature.");
   if(TKernel==KERNEL_Cubic) Run_Exceptioon("Cubic Kernel is not available with Variable resolution execution.");
   if(TStep==STEP_Verlet) Run_Exceptioon("Verlet Time Integration is not available with Variable resolution execution.");
   if(Shifting)Log->PrintWarning("Shifting==4 is recommended with Variable resolution execution.");
   #ifdef AVAILABLE_DIVCLEAN
   if(DivCleaning) Run_Exceptioon("Divergence cleaning is not available with Variable resolution execution.");
   #endif
}

//==============================================================================
/// Return parameters for VRes coupling.
//==============================================================================
StInterParmsbg JSphGpuSingle_VRes::GetVResParms(){
	StrInterParmsbg parms=StrInterParmsbg(Simulate2D,TKernel
		    ,DivData,CTE,Map_PosMin,Dcell_g->cptr()
		    ,Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr()
        ,Velrho_g->cptr(),Idp_g->cptr(),Code_g->cptr()
		    ,NULL,NULL);
 return(parms);
}

//==============================================================================
/// Initial definition of buffer particles and reordering.
//==============================================================================
void JSphGpuSingle_VRes::BufferInit(StInterParmsbg *parms){
  cusph::CteInteractionUp(&CTE);

  acfloat3* boundnormc=new acfloat3("-",Arrays_Cpu,false);
  boundnormc->Reserve();
  if(TBoundary==BC_MDBC)BoundNor_g->CuCopyToHost(boundnormc,Np);
  VRes->CheckNormals(TBoundary,Npb,0,AuxPos_c->cptr()
    ,Idp_c->cptr(),Code_c->ptr(),AC_CPTR(boundnormc),VResId);
  boundnormc->Free();

	for(unsigned i=0;i<VRes->GetCount();i++){

		agint bufferpartg("-",Arrays_Gpu,true);
		const unsigned buffercountpre=VRes->CreateListGpuInit(Np,0,Posxy_g->cptr()
      ,Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize,bufferpartg.ptr(),i);   
	}

  RunCellDivide(true);
  if(VRES_DG_SAVEVTK)DgSaveVtkParticlesGpuMultiRes("Compute_step",Nstep,0,Np
    ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),Idp_g->cptr(),Velrho_g->cptr());

}

//==============================================================================
/// Perform interpolation of buffer particles over coupled VRes simulation.
//==============================================================================
void JSphGpuSingle_VRes::BufferExtrapolateData(StInterParmsbg *parms){
  const unsigned count=VRes->GetCount();
	for(unsigned i=0;i<count;i++){
    
    //-Compute list of buffer particles.	
    agint bufferpartg("-",Arrays_Gpu,true);
		const unsigned buffercountpre=VRes->CreateListGpu(Np-Npb,Npb
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize
      ,bufferpartg.ptr(),i);
	
    //-Update constant memory and perform interpolation.  
		unsigned id=VRes->GetZone(i)->getZone()-1;
		cusph::CteInteractionUp(&parms[id].cte);
		cusphvres::Interaction_BufferExtrap(buffercountpre,bufferpartg.ptr()
      ,parms[id],Posxy_g->cptr(),Posz_g->cptr(),Velrho_g->ptr(),Code_g->ptr()
      ,VResFastSingle,VResOrder,VResMethod,VResThreshold);
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
  cusph::CteInteractionUp(&CTE);                            //-Update constant memory.

  //-Compute movement for VRes regions
  VRes->UpdateMatMov(mat);
	VRes->MoveBufferZone(dt,mat);

	//- ComputeStep for each buffer region.
	for(unsigned i=0;i<VRes->GetCount();i++){
    
    StrDataVresGpu vresdata=VRes->GetZoneFluxInfoGpu(i);  //-Retrieve buffer parameters.			 

    //-Compute eulerian flux on the mass accumulation points.
  	unsigned id=VRes->GetZone(i)->getZone()-1;      
		cusph::CteInteractionUp(&parms[id].cte);                  //-Update constant memory.
    cusphvres::Interaction_BufferExtrapFlux(parms[id]
      ,vresdata,CTE.dp,dt,VResFastSingle,VResOrder,VResMethod,VResThreshold);

		
    //-Compute step on buffer particles.
		agint bufferpartg("-",Arrays_Gpu,true);       
		const unsigned buffercountpre=VRes->CreateListGpu(Np-Npb,Npb
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->ptr(),GpuParticlesSize
      ,bufferpartg.ptr(),i);
    VRes->ComputeStepGpu(buffercountpre,bufferpartg.ptr(),IdMax+1,Posxy_g->ptr()
      ,Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr(),Idp_g->ptr(),Velrho_g->ptr(),NULL,this,i);

		
    cusph::CteInteractionUp(&CTE);                            //-Update constant memory.
		
    cusphvres::CheckMassFlux(vresdata.ntot,vresdata.nini,DivData,Map_PosMin
      ,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr(),PosCell_g->cptr()
      ,vresdata.ptposxy,vresdata.ptposz,vresdata.normals,vresdata.mass);
		
    int *newpart=NULL;
		cudaMalloc((void**)&newpart, sizeof(int)*(vresdata.ntot+1));
		cudaMemset(newpart,0,vresdata.ntot*sizeof(int));


    unsigned newnp = VRes->NewPartListCreate(newpart,i);
		
    if(newnp){			
      if(!CheckGpuParticlesSize(Np+newnp)){ //-Check memory allocation and resize.
        const unsigned ndatacpu=0,ndatagpu=Np;
        ResizeParticlesSizeData(ndatacpu,ndatagpu,Np+newnp,Np+newnp,0.2f,true);
        CellDivSingle->SetIncreaseNp(newnp);
			}  
      VRes->CreateNewPartGpu(Np,newnp,IdMax+1,Posxy_g->ptr(),Posz_g->ptr()
        ,Dcell_g->ptr(),Code_g->ptr(),Idp_g->ptr(),Velrho_g->ptr(),newpart,i);
      
      //-Updates basic arrays.
      if(SpsTauRho2_g)  SpsTauRho2_g->CuMemsetOffset(Np,0,newnp);
      if(BoundNor_g)    BoundNor_g  ->CuMemsetOffset(Np,0,newnp);
      if(FSType_g)      FSType_g    ->CuMemsetOffset(Np,0,newnp);
      if(ShiftVel_g)    ShiftVel_g  ->CuMemsetOffset(Np,0,newnp);
      #ifdef AVAILABLE_DIVCLEAN
      if(PsiClean_g)PsiClean_g->CuMemsetOffset(Np,0,newnp);
      #endif
      Np+=newnp; 
      TotalNp+=newnp;
      IdMax=unsigned(TotalNp-1);
    }      
    cudaFree(newpart);  newpart=NULL; 
	}
}

//==============================================================================
/// Remove the normal component to the interface of the shifting for buffer particles.
//==============================================================================
void JSphGpuSingle_VRes::BufferShifting(){
  StrGeomVresGpu vresgdata=VRes->GetGeomInfoVresGpu();
	cusphvres::BufferShiftingGpu(Np,Npb,Posxy_g->ptr(),Posz_g->ptr(),ShiftVel_g->ptr(),Code_g->ptr(),vresgdata,NULL);
}

//==============================================================================
/// Wrapper for call particle sorting procedure in VResDriver.
//==============================================================================
void JSphGpuSingle_VRes::CallRunCellDivide(){
  	cusph::CteInteractionUp(&CTE);
    RunCellDivide(true);
}

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
     StrGeomVresGpu vresgdata=VRes->GetGeomInfoVresGpu();

  if(runshift)cusphvres::PreLoopInteraction(TKernel,Simulate2D,runshift,bsfluid,Np-Npb,Npb,DivData
    ,Posxy_g->cptr(),Posz_g->cptr(),Dcell_g->cptr(),PosCell_g->cptr()
    ,Velrho_g->cptr(),Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr()
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
    if(0 && runshift && TimeStep+LastDt>=TimePartNext){
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
  StrGeomVresGpu vresgdata=VRes->GetGeomInfoVresGpu();

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
  StrGeomVresGpu vresgdata=VRes->GetGeomInfoVresGpu();

  aguint    inoutpartg("-",Arrays_Gpu,true);
  cusphvres::ComputeUmbrellaRegion(TKernel,Simulate2D,bsfluid,Npb,Np-Npb,DivData
    ,Dcell_g->cptr(),Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr()
    ,Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr(),FSType_g->ptr(),FSNormal_g->ptr()
    ,inoutpartg.ptr(),vresgdata,NULL);
}

//==============================================================================
/// Return movement of an object for VRes tracking.
//==============================================================================
JMatrix4d JSphGpuSingle_VRes::CalcVelMotion(unsigned trackingmk,double dt){
  //-Check if the mkbound a motion object.
  JMatrix4d mat;
  if(CaseNmoving){
    const unsigned nref=DsMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      const unsigned mk=m.mkbound;
      if(trackingmk==mk) return(CalcMotionMoving(m,dt));
    }
  }
  //-Check if the mkbound is a floating object.
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

//==============================================================================
/// Return movement if mkbound is a moving object.
//==============================================================================
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

//==============================================================================
/// Return movement if mkbound is a floating object.
//==============================================================================
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

  AbortNoNormals=false;

  //-Selection of GPU.
  const int gpuid=SelecDevice(cfg->GpuId);
  LoadVResConfigParameters(cfg);

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
  VisuConfigVRes();
  ConfigDomain();
  ConfigRunMode();
  VisuParticleSummary();
  
  //-Store Interaction parameters for vres.
  GetConstantData(CTE);


  //-Initialisation of execution variables.
  InitRunGpu();
  RunFirstGaugeSystem(TimeStep);
  if(InOut)InOutInit(TimeStepIni);
  if(FlexStruc)FlexStrucInit(); //<vs_flexstruc>
  FreePartsInit();
  PrintAllocMemory(GetAllocMemoryCpu(),GetAllocMemoryGpu());
  UpdateMaxValues();
  SaveData(); 
  AddWarningVRes();
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
SymplecticDtPre1=SymplecticDtPre;

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
double JSphGpuSingle_VRes::ComputeStepVRes(){

  cusph::CteInteractionUp(&CTE);
  const double dt=SymplecticDtPre;
  if(CaseNmoving)CalcMotion(dt);                //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
  InterStep=INTERSTEP_SymPredictor;
  DemDtForce=dt*0.5f;                           //-For DEM interaction.
  MdbcBoundCorrection(InterStep);               //-Mdbc correction
  PreInteraction_Forces(InterStep);             //-Allocating temporary arrays.
  PreLoopProcedureVRes(InterStep);              //-Calculate variables for interaction forces (Shifting,DDT,etc...).  
  Interaction_Forces(InterStep);                //-Interaction.
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
  if(FlexStruc)UpdateFlexStrucGeometry();       //-Update the geometric information for each flexible structure particle.
  MdbcBoundCorrection(InterStep);               //-Mdbc correction
  PreInteraction_Forces(InterStep);             //-Allocating temporary arrays.
  PreLoopProcedureVRes(InterStep);              //-Calculate variables for interaction forces (Shifting,DDT,etc...).
  Interaction_Forces(InterStep);                //-Interaction.
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
	else RunCellDivide(true);
  if(FlexStruc)UpdateFlexStrucGeometry();
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
    if(VRES_DG_SAVEVTK)VRes->SaveNormals(DirDataOut+"NormalsBuffer_",Part);
    if(VRES_DG_SAVEVTK)VRes->SaveVtkDomains(DirDataOut+"Domains",Part,Simulate2D);
    if(VRES_DG_SAVEVTK)VRes->SaveVResData(Part,TimeStep,Nstep);

    SaveData();
    Part++;
    PartNstep=Nstep;
    TimeStepM1=TimeStep;
    TimePartNext=(SvAllSteps? TimeStep: OutputTime->GetNextTime(TimeStep));
    TimerPart.Start();
	}
	UpdateMaxValues();
}

//==============================================================================
/// Complete VRes simulation.
//==============================================================================
void JSphGpuSingle_VRes::Finish2(){
	TimerSim.Stop(); TimerTot.Stop();
	//-End of Simulation.
	//--------------------
	FinishRun(false);
}