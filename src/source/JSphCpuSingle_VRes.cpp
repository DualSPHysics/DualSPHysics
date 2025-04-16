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

/// \file JSphCpuSingle_VRes.cpp \brief Implements the class \ref JSphCpuSingle_VRes.

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
//==============================================================================
/// Constructor.
//==============================================================================
JSphCpuSingle_VRes::JSphCpuSingle_VRes():JSphCpuSingle(){
  ClassName="JSphGpuSingle";
  VResThreshold=100;
  VRes=NULL;
  VResOrder=VrOrder_1st;
  VResMethod=VrMethod_Liu;
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphCpuSingle_VRes::~JSphCpuSingle_VRes(){
  DestructorActive=true;
  delete VRes; VRes=NULL;
}

//==============================================================================
/// Load VRes configuration.
//==============================================================================
void JSphCpuSingle_VRes::LoadVResConfigParameters(const JSphCfgRun* cfg){
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
  if(cfg->VResThreshold>=0) VResThreshold=cfg->VResThreshold;
}

//==============================================================================
/// Print VRes configuration.
//==============================================================================
void JSphCpuSingle_VRes::VisuConfigVRes(){
  Log->Print("\nVRes Configuration:");
  Log->Printf(" Interpolation Method: %s",(VResMethod==VrMethod_Liu? "Liu-Liu Correction": "Moving Least Square"));
  Log->Printf(" Interpolation Order: %s" ,(VResOrder==VrOrder_2nd? "2nd": (VResOrder==VrOrder_1st? "1st": "0th")));
  Log->Printf(" Inter Threshold: %g" ,VResThreshold);
}

//==============================================================================
/// Initialize VRes object.
//==============================================================================
void JSphCpuSingle_VRes::VResInit(const JSphCfgRun* cfg, JCaseVRes casemultires, unsigned id){
   VRes=new JSphVRes(Cpu,CSP,casemultires,id,MapRealPosMin
   ,MapRealPosMax,AppName,DirDataOut,PartBegin,PartBeginDir);
   VRes->Config();
}

//==============================================================================
/// Add warning for VRes execution.
//==============================================================================
void JSphCpuSingle_VRes::AddWarningVRes(){
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
stinterparmscb JSphCpuSingle_VRes::GetVResParms(){
	stinterparmscb parms=StInterparmscb(Np,Npb,NpbOk,DivData
    ,Dcell_c->cptr(),Pos_c->cptr(),Velrho_c->cptr()
    ,Idp_c->cptr(),Code_c->cptr(),CSP);  
	return(parms);
}

//==============================================================================
/// Initial definition of buffer particles and reordering.
//==============================================================================
void JSphCpuSingle_VRes::BufferInit(stinterparmscb *parms){
  VRes->CheckNormals(TBoundary,Npb,0,Pos_c->cptr()
    ,Idp_c->cptr(),Code_c->ptr(),AC_CPTR(BoundNor_c),VResId);
  for(unsigned i=0;i<VRes->GetCount();i++){
    acint bufferpartc("-",Arrays_Cpu,true);
    unsigned buffercountpre=VRes->CreateListCpuInit(Np,0,Pos_c->cptr()
      ,Idp_c->cptr(),Code_c->ptr(),bufferpartc.ptr(),i);
    }
    RunCellDivide(true);
    if(VRES_DG_SAVEVTK)DgSaveVtkParticlesCpuMRBuffer("Compute_step",Nstep,0,Np
      ,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL);

}

//==============================================================================
/// Perform interpolation of buffer particles over coupled VRes simulation:
//==============================================================================
void JSphCpuSingle_VRes::BufferExtrapolateData(stinterparmscb *parms){
  const unsigned count=VRes->GetCount();
	for(unsigned i=0;i<count;i++){
		
    //-Compute list of buffer particles.	
    acint bufferpartc("-",Arrays_Cpu,true);
    unsigned buffercountpre=VRes->CreateListCpuInit(Np,0,Pos_c->cptr()
      ,Idp_c->cptr(),Code_c->ptr(),bufferpartc.ptr(),i);

    //-Update constant memory and perform interpolation.  
		unsigned id=VRes->GetZone(i)->getZone()-1;
    fvres::Interaction_BufferExtrap(buffercountpre,bufferpartc.ptr()
      ,parms[id],Pos_c->ptr(),Velrho_c->ptr(),Code_c->ptr()
      ,VResOrder,VResMethod,VResThreshold);
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
void JSphCpuSingle_VRes::ComputeStepBuffer(double dt,std::vector<JMatrix4d> mat,stinterparmscb *parms){
  //-Compute movement for VRes regions
  VRes->UpdateMatMov(mat);
	VRes->MoveBufferZone(dt,mat);
	
  //- ComputeStep for each buffer region.
	for(unsigned i=0;i<VRes->GetCount();i++){

    StrDataVresCpu vresdata=VRes->GetZoneFluxInfoCpu(i);    //-Retrieve buffer parameters

    //-Compute eulerian flux on the mass accumulation points.
		unsigned id=VRes->GetZone(i)->getZone()-1;	
    fvres::Interaction_BufferExtrapFlux(parms[id],vresdata,Dp,dt
      ,VResOrder,VResMethod,VResThreshold);
    
    //-Compute step on buffer particles.
    acint bufferpartc("-",Arrays_Cpu,true);
    unsigned buffercountpre=VRes->CreateListCpuInit(Np,0,Pos_c->cptr()
      ,Idp_c->cptr(),Code_c->ptr(),bufferpartc.ptr(),i);

		fvres::CheckMassFlux(vresdata.ntot,vresdata.nini,CSP,DivData
      ,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
      ,vresdata.points,vresdata.normals,vresdata.mass);
    
    unsigned newnp=VRes->ComputeStepCpu(buffercountpre,bufferpartc.ptr(),Code_c->ptr(),Pos_c->cptr(),i);

    if(newnp){
      if(!CheckCpuParticlesSize(Np+newnp)){
        const unsigned ndatacpu=Np;
        ResizeParticlesSizeData(ndatacpu,Np+newnp,Np+newnp,0.2f,true);
        CellDivSingle->SetIncreaseNp(newnp);
      }
      VRes->CreateNewPart(IdMax+1,Dcell_c->ptr(),Code_c->ptr()
        ,Pos_c->ptr(),Idp_c->ptr(),Velrho_c->ptr(),this,Np,i);

    }

    if(newnp){      
      if(SpsTauRho2_c)  SpsTauRho2_c->MemsetOffset(Np,0,newnp);
      if(BoundNor_c)    BoundNor_c->MemsetOffset(Np,0,newnp);
      if(FSType_c)      FSType_c->MemsetOffset(Np,0,newnp);
      if(ShiftVel_c)    ShiftVel_c->MemsetOffset(Np,0,newnp);
      #ifdef AVAILABLE_DIVCLEAN
      if(PsiClean_c)PsiClean_c->MemsetOffset(Np,0,newnp);
      #endif
      Np+=newnp; 
      TotalNp+=newnp;
      IdMax=unsigned(TotalNp-1);
		} 
	}
}

void JSphCpuSingle_VRes::BufferShifting(){
  // StrGeomVresGpu& vresgdata=VRes->GetGeomInfoVres();
	// cusphvres::BufferShiftingGpu(Np,Npb,Posxy_g->ptr(),Posz_g->ptr(),ShiftVel_g->ptr(),Code_g->ptr(),vresgdata,NULL);
}

//==============================================================================
/// Wrapper for call particle sorting procedure in VResDriver.
//==============================================================================
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

    StrGeomVresCpu vresdata=VRes->GetGeomInfoVresCpu();
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
    if(0 && runshift && TimeStep+LastDt>=TimePartNext)DgSaveVtkParticlesCpu("Compute_FreeSurface_",Part,0,Np,Pos_c->cptr()
      ,Code_c->cptr(),FSType_c->cptr(),ShiftVel_c->cptr(),FSNormal_c->cptr());
    Timersc->TmStop(TMC_SuShifting);
  }
}

//==============================================================================
/// Compute free-surface particles and their normals.
//==============================================================================
void JSphCpuSingle_VRes::ComputeFSParticlesVRes(){
  StrGeomVresCpu vresdata=VRes->GetGeomInfoVresCpu();

  acuint fspart("-",Arrays_Cpu,true);
  fvres::CallComputeFSNormals(Np,Npb,CSP,DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,Velrho_c->cptr(),FSType_c->ptr(),FSNormal_c->ptr(),fspart.ptr(),vresdata);
}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
//==============================================================================
void JSphCpuSingle_VRes::ComputeUmbrellaRegionVRes(){
  StrGeomVresCpu vresdata=VRes->GetGeomInfoVresCpu();

  acuint fspart("-",Arrays_Cpu,true);
  fvres::CallScanUmbrellaRegion(Np,Npb,CSP,DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,FSNormal_c->cptr(),fspart.ptr(),FSType_c->ptr(),vresdata);
}

//==============================================================================
/// Return movement of an object for VRes tracking.
//==============================================================================
JMatrix4d JSphCpuSingle_VRes::CalcVelMotion(unsigned trackingmk,double dt){
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

//==============================================================================
/// Return movement if mkbound is a floating object.
//==============================================================================
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

//==============================================================================
/// Initialises VRes simulation.
//==============================================================================
void JSphCpuSingle_VRes::Init(std::string appname,const JSphCfgRun* cfg,JLog2* log
  ,unsigned vrescount,unsigned vresid)
{
  if(!cfg || !log)return;
  AppName=appname; Log=log; CfgRun=cfg;
  VResCount=vrescount;
  VResId=vresid; 
  AbortNoNormals=false;

  LoadVResConfigParameters(cfg);

  //-Creates array system for particles.
  Arrays_Cpu=new JArraysCpu(Log);

  //-Configure timers.
  Timersc->Config(cfg->SvTimers);
  Timersc->TmStart(TMC_Init);

  //-Load parameters and input data.
  LoadConfig(cfg);
  LoadCaseParticles();
  VisuConfig();
  VisuConfigVRes();
  ConfigDomain();
  ConfigRunMode();
  VisuParticleSummary();

  //-Initialisation of execution variables.
  InitRunCpu();
  RunFirstGaugeSystem(TimeStep);
  if(InOut)InOutInit(TimeStepIni);
  if(FlexStruc)FlexStrucInit(); //<vs_flexstruc>
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

//==============================================================================
/// Complete initialization of VRes simulation.
//==============================================================================
double JSphCpuSingle_VRes::Init2(){
  SymplecticDtPre1=SymplecticDtPre;

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

//==============================================================================
/// Compute time step of VRes simulation.
//==============================================================================
double JSphCpuSingle_VRes::ComputeStepVRes(){
  const double dt=SymplecticDtPre;
  if(CaseNmoving)CalcMotion(dt);          //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
  InterStep=INTERSTEP_SymPredictor;
  DemDtForce=dt*0.5f;
  MdbcBoundCorrection(InterStep);         //-Mdbc correction
  PreInteraction_Forces(InterStep);       //-Allocating temporary arrays.
  PreLoopProcedureVRes(InterStep);        //-Calculate variables for interaction forces (Shifting,DDT,etc...).                          //-For DEM interaction.
  Interaction_Forces(InterStep);          //-Interaction.
  const double dt_p=DtVariable(false);    //-Calculate dt of predictor step.
  if(Shifting)RunShifting(dt*.5);         //-Shifting.
  ComputeSymplecticPre(dt);               //-Apply Symplectic-Predictor to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt*.5,true);  //-Control of floating bodies.
  PosInteraction_Forces();                //-Free memory used for interaction.
  //-Corrector
  //-----------
  InterStep=INTERSTEP_SymCorrector;
  DemDtForce=dt;                          //-For DEM interaction.
  RunCellDivide(true);
  if(FlexStruc)UpdateFlexStrucGeometry(); //-Update the geometric information for each flexible structure particle.
  MdbcBoundCorrection(InterStep);         //-Mdbc correction
  PreInteraction_Forces(InterStep);       //-Allocating temporary arrays.
  PreLoopProcedureVRes(InterStep);        //-Calculate variables for interaction forces (Shifting,DDT,etc...).
  Interaction_Forces(InterStep);          //-Interaction.
  const double dt_c=DtVariable(true);     //-Calculate dt of corrector step.
  if(Shifting)RunShifting(dt);            //-Shifting.
  ComputeSymplecticCorr(dt);              //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);    //-Control of floating bodies.
  PosInteraction_Forces();                //-Free memory used for interaction.
  if(Damping)RunDamping(dt);              //-Applies Damping.
  if(RelaxZones)RunRelaxZone(dt);         //-Generate waves using RZ.
  SymplecticDtPre1=min(dt_p,dt_c);        //-Calculate dt for next ComputeStep.
  return(dt);
}

//==============================================================================
/// Complete time step of VRes simulation.
//==============================================================================
void JSphCpuSingle_VRes::Finish(double dt1){
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
    if(PartMotionSave)PartMotionSave->AddDataExtraCpu(Part,TimeStep,Nstep,Np
      ,Pos_c->cptr(),RidpMot);
    TimeOutExtraUpdate(TimeStep);
  }
	if(TimeStep>=TimePartNext ){
	  if(VRES_DG_SAVEVTK){
      DgSaveVtkParticlesCpuMRBuffer("Compute_Debug_Buffer.vtk",Part,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL);
	    DgSaveVtkParticlesCpuMR("Compute_Debug_Fluid.vtk",Part,0,Np,Pos_c->cptr(),Code_c->cptr(),Idp_c->cptr(),Velrho_c->cptr(),NULL,NULL,NULL,NULL);
    }
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
void JSphCpuSingle_VRes::Finish2(){
	TimerSim.Stop(); TimerTot.Stop();
	//-End of Simulation.
	//--------------------
	FinishRun(false);

}