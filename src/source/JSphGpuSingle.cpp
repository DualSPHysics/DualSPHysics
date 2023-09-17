//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2023 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphGpuSingle.cpp \brief Implements the class \ref JSphGpuSingle.

#include "JSphGpuSingle.h"
#include "JCellDivGpuSingle.h"
#include "JSphMk.h"
#include "JPartsLoad4.h"
#include "Functions.h"
#include "JXml.h"
#include "JDsMotion.h"
#include "JDsViscoInput.h"
#include "JWaveGen.h"
#include "JMLPistons.h"
#include "JRelaxZones.h"
#include "JChronoObjects.h"
#include "JDsMooredFloatings.h"
#include "JDsFtForcePoints.h"
#include "JDsOutputTime.h"
#include "JTimeControl.h"
#include "JSphGpu_ker.h"
#include "JSphGpuSimple_ker.h"
#include "JDsGaugeSystem.h"
#include "JSphInOut.h"
#include "JDsPartMotionSave.h"
#include "JDsPartFloatSave.h"
#include "JLinearValue.h"
#include "JDataArrays.h"
#include "JDebugSphGpu.h"
#include "JSphShifting.h"
#include "JDsPips.h"
#include "JDsExtraData.h"
#include "FunctionsCuda.h"
#include "JDsOutputParts.h" //<vs_outpaarts>

#include <climits>

using namespace std;
//==============================================================================
/// Constructor.
//==============================================================================
JSphGpuSingle::JSphGpuSingle():JSphGpu(false){
  ClassName="JSphGpuSingle";
  CellDivSingle=NULL;
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphGpuSingle::~JSphGpuSingle(){
  DestructorActive=true;
  delete CellDivSingle; CellDivSingle=NULL;
}

//==============================================================================
/// Returns the memory allocated to the CPU.
/// Devuelve la memoria reservada en CPU.
//==============================================================================
llong JSphGpuSingle::GetAllocMemoryCpu()const{  
  llong s=JSphGpu::GetAllocMemoryCpu();
  //-Allocated in other objects.
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemoryCpu();
  return(s);
}

//==============================================================================
/// Returns the memory allocated to the GPU.
/// Devuelve la memoria reservada en GPU.
//==============================================================================
llong JSphGpuSingle::GetAllocMemoryGpu()const{  
  llong s=JSphGpu::GetAllocMemoryGpu();
  //-Allocated in other objects.
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemoryGpu();
  return(s);
}

//==============================================================================
/// Returns the GPU memory allocated or used for particles
/// Devuelve la memoria GPU reservada o usada para particulas.
//==============================================================================
llong JSphGpuSingle::GetMemoryGpuNp()const{
  llong s=JSphGpu::GetAllocMemoryGpu();
  //-Allocated in other objects.
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemoryGpuNp();
  return(s);
}

//==============================================================================
/// Returns the GPU memory allocated or used for cells.
/// Devuelve la memoria GPU reservada o usada para celdas.
//==============================================================================
llong JSphGpuSingle::GetMemoryGpuNct()const{
  return(CellDivSingle? CellDivSingle->GetAllocMemoryGpuNct(): 0);
}

//==============================================================================
/// Updates the maximum values of memory, particles and cells.
/// Actualiza los valores maximos de memory, particles y cells.
//==============================================================================
void JSphGpuSingle::UpdateMaxValues(){
  const llong mcpu=GetAllocMemoryCpu();
  const llong mgpu=GetAllocMemoryGpu();
  MaxNumbers.memcpu=max(MaxNumbers.memcpu,mcpu);
  if(mgpu>MaxNumbers.memgpu){
     MaxNumbers.memgpu=mgpu;
     MaxNumbers.memgpunct=GetMemoryGpuNct();
  }
  MaxNumbers.particles=max(MaxNumbers.particles,Np);
  if(CellDivSingle)MaxNumbers.cells=max(MaxNumbers.cells,CellDivSingle->GetNct());
  //-Prints allocated arrays when it is changed.
  if(Arrays_Cpu->GetArrayCountUpdated() || Arrays_Gpu->GetArrayCountUpdated()){
    PrintAllocMemory(GetAllocMemoryCpu(),GetAllocMemoryGpu());
  }
}

//==============================================================================
/// Loads the configuration of the execution.
/// Carga la configuracion de ejecucion.
//==============================================================================
void JSphGpuSingle::LoadConfig(const JSphCfgRun* cfg){
  //-Loads general configuration.
  JSph::LoadConfig(cfg);
  //-Checks compatibility of selected options.
  Log->Print("**Special case configuration is loaded");
}

//==============================================================================
/// Configuration of the current domain.
/// Configuracion del dominio actual.
//==============================================================================
void JSphGpuSingle::ConfigDomain(){
  //-Configure cell map division (defines ScellDiv, Scell, Map_Cells). 
  ConfigCellDivision();
  //-Computes the number of particles.
  Np=PartsLoaded->GetCount();
  Npb=CaseNpb;
  NpbOk=Npb;

  //-Allocates fixed CPU memory for moving & floating particles.
  AllocCpuMemoryFixed();
  //-Allocates CPU memory for particles.
  AllocCpuMemoryParticles(Np);

  //-Copies particle data.
  AuxPos_c->CopyFrom(PartsLoaded->GetPos(),Np);
  Idp_c   ->CopyFrom(PartsLoaded->GetIdp(),Np);
  Velrho_c->CopyFrom(PartsLoaded->GetVelRho(),Np);

  //-Computes radius of floating bodies.
  if(CaseNfloat && PeriActive!=0 && !PartBegin)
    CalcFloatingRadius(Np,AuxPos_c->cptr(),Idp_c->cptr());

  //-Configures Multi-Layer Pistons according to particles.
  if(MLPistons)MLPistons->PreparePiston(Dp,Np,Idp_c->cptr(),AuxPos_c->cptr());

  //-Loads Code of the particles.
  LoadCodeParticles(Np,Idp_c->cptr(),Code_c->ptr());

  //-Load normals for boundary particles (fixed and moving).
  acfloat3 boundnorc("boundnor",Arrays_Cpu,UseNormals);
  if(UseNormals)LoadBoundNormals(Np,Npb,Idp_c->cptr(),Code_c->cptr(),boundnorc.ptr());

  //-Creates PartsInit object with initial particle data for automatic configurations.
  CreatePartsInit(Np,AuxPos_c->cptr(),Code_c->cptr());

  //-Runs initialization operations from XML.
  RunInitialize(Np,Npb,AuxPos_c->cptr(),Idp_c->cptr(),Code_c->cptr()
    ,Velrho_c->ptr(),boundnorc.ptr());
  if(UseNormals)ConfigBoundNormals(Np,Npb,AuxPos_c->cptr(),Idp_c->cptr(),boundnorc.ptr());

  //-Computes MK domain for boundary and fluid particles.
  MkInfo->ComputeMkDomains(Np,AuxPos_c->cptr(),Code_c->cptr());

  //-Sets local domain of the simulation within Map_Cells and computes DomCellCode.
  //-Establece dominio de simulacion local dentro de Map_Cells y calcula DomCellCode.
  SelecDomain(TUint3(0,0,0),Map_Cells);
  //-Computes inital cell of the particles and checks if there are unexpected excluded particles.
  //-Calcula celda inicial de particulas y comprueba si hay excluidas inesperadas.
  LoadDcellParticles(Np,Code_c->cptr(),AuxPos_c->cptr(),Dcell_c->ptr());

  //-Allocates fixed GPU memory for moving & floating particles.
  AllocGpuMemoryFixed();
  //-Allocates GPU memory for particles.
  AllocGpuMemoryParticles(Np);

  //-Uploads particle data on the GPU.
  Pos3ToPos21(Np,AuxPos_c->cptr(),Posxy_c->ptr(),Posz_c->ptr());
  ParticlesDataUp(Np,boundnorc.cptr());
  boundnorc.Free();

  //-Uploads constants on the GPU.
  ConstantDataUp();

  //-Creates object for Celldiv on the GPU and selects a valid cellmode.
  //-Crea objeto para divide en GPU y selecciona un cellmode valido.
  CellDivSingle=new JCellDivGpuSingle(Stable,FtCount!=0,PeriActive,KernelSize2,PosCellSize
    ,CellDomFixed,CellMode,Scell,Map_PosMin,Map_PosMax,Map_Cells,CaseNbound,CaseNfixed,CaseNpb,DirOut);
  CellDivSingle->DefineDomain(DomCellCode,DomCelIni,DomCelFin,DomPosMin,DomPosMax);
  ConfigCellDiv((JCellDivGpu*)CellDivSingle);

  ConfigBlockSizes(false,PeriActive!=0);
  ConfigSaveData(0,1,"",Np,AuxPos_c->cptr(),Idp_c->cptr());

  //-Reorders particles according to cells.
  BoundChanged=true;
  RunCellDivide(true);
}

//==============================================================================
/// Resizes the allocated memory for particles on the CPU and the GPU saving 
/// the current data (ndatacpu and ndatagpu) and measures the time spent 
/// using TMG_SuResizeNp.
/// At the end updates the division.
///
/// Redimensiona el espacio reservado para particulas en CPU y GPU manteniendo 
/// los datos actuales (ndatacpu y ndatagpu) y mide el tiempo con TMG_SuResizeNp.
/// Al terminar actualiza el divide.
//==============================================================================
void JSphGpuSingle::ResizeParticlesSizeData(unsigned ndatacpu,unsigned ndatagpu
  ,unsigned newsize,unsigned minsize,float oversize,bool updatedivide)
{
  Timersg->TmStart(TMG_SuResizeNp,false);
  newsize=newsize+(oversize>0? unsigned(oversize*newsize): 0);
  //-Free up almost all of the CPU allocated memory when saving is not required.
  if(!ndatacpu)Arrays_Cpu->SetDataArraySize(1,0);
  //-Free GPU memory for NL.
  CellDivSingle->FreeMemoryGpu();
  DivData=DivDataGpuNull();
  //-Resize GPU memory for particles saving current data.
  ResizeGpuMemoryParticlesData(ndatagpu,newsize,minsize);
  //-Resize CPU memory for particles saving current data when is required.
  CpuParticlesSize=GpuParticlesSize;
  Arrays_Cpu->SetDataArraySize(CpuParticlesSize,ndatacpu);
  Timersg->TmStop(TMG_SuResizeNp,true);
  if(updatedivide)RunCellDivide(true);
}

//==============================================================================
/// Creates duplicate particles for periodic conditions.
/// Creates new periodic particles and marks the old ones to be ignored.
/// The new particles are lccated from the value of Np, first the NpbPer for 
/// boundaries and then the NpfPer for the fluids. The Np output also contains 
/// the new periodic particles.
///
/// Crea particulas duplicadas de condiciones periodicas.
/// Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
/// Las nuevas periodicas se situan a partir del Np de entrada, primero las NpbPer
/// de contorno y despues las NpfPer fluidas. El Np de salida contiene tambien las
/// nuevas periodicas.
//==============================================================================
void JSphGpuSingle::RunPeriodic(){
  Timersg->TmStart(TMG_SuPeriodic,false);
  //-Stores the current number of periodic particles.
  //-Guarda numero de periodicas actuales.
  NpfPerM1=NpfPer;
  NpbPerM1=NpbPer;
  //-Marks current periodic particles to be ignored.
  //-Marca periodicas actuales para ignorar.
  cusph::PeriodicIgnore(Np,Code_g->ptr());
  //-Creates new periodic particles.
  //-Crea las nuevas periodicas.
  aguint listpg("listp",Arrays_Gpu,true);
  const unsigned npb0=Npb;
  const unsigned npf0=Np-Npb;
  const unsigned np0=Np;
  NpbPer=NpfPer=0;
  BoundChanged=true;
  for(unsigned ctype=0;ctype<2;ctype++){//-0:bound, 1:fluid+floating.
    //-Calculate range of particles to be examined (bound or fluid).
    //-Calcula rango de particulas a examinar (bound o fluid).
    const unsigned pini=(ctype? npb0: 0);
    const unsigned num= (ctype? npf0: npb0);
    //-Search for periodic particles in each direction (X, Y, or Z).
    //-Busca periodicas en cada eje (X, Y e Z).
    for(unsigned cper=0;cper<3;cper++)if((cper==0 && PeriX) || (cper==1 && PeriY) || (cper==2 && PeriZ)){
      const tdouble3 perinc=(cper==0? PeriXinc: (cper==1? PeriYinc: PeriZinc));
      //-First search in the list of new periodic particles and then in the initial list of particles (this is needed for periodic particles in more than one direction).
      //-Primero busca en la lista de periodicas nuevas y despues en la lista inicial de particulas (necesario para periodicas en mas de un eje).
      for(unsigned cblock=0;cblock<2;cblock++){//-0:new periodic particles, 1:original periodic particles
        const unsigned nper=(ctype? NpfPer: NpbPer);  //-number of new periodic particles for the type currently computed (bound or fluid). | Numero de periodicas nuevas del tipo a procesar.
        const unsigned pini2=(cblock? pini: Np-nper);
        const unsigned num2= (cblock? num:  nper);
        //-Repeats search if the available memory was insufficient and had to be increased.
        //-Repite la busqueda si la memoria disponible resulto insuficiente y hubo que aumentarla.
        bool run=true;
        while(run && num2){
          //-Reserve memory to create list of periodic particles.
          //-Reserva memoria para crear lista de particulas periodicas.
          listpg.Reserve();
          const unsigned nmax=GpuParticlesSize-1; //-Maximum number of particles that can be included in the list. | Numero maximo de particulas que caben en la lista.
          //-Generates list of new periodic particles.
          if(Np>=0x80000000)Run_Exceptioon("The number of particles is too big.");//-Because the last bit is used to mark the reason the new periodical is created. | Porque el ultimo bit se usa para marcar el sentido en que se crea la nueva periodica. 
          const unsigned count=cusph::PeriodicMakeList(num2,pini2,Stable
            ,nmax,Map_PosMin,Map_PosMax,perinc,Posxy_g->cptr()
            ,Posz_g->cptr(),Code_g->cptr(),listpg.ptr());
          //-Resizes the allocated memory for the particles if there is not sufficient space and repeats the serach process.
          //-Redimensiona memoria para particulas si no hay espacio suficiente y repite el proceso de busqueda.
          if(count>nmax || !CheckGpuParticlesSize(count+Np)){
            listpg.Free(); //-Avoids unnecessary copying of its data during resizing.
            Timersg->TmStop(TMG_SuPeriodic,true);
            const unsigned ndatacpu=0,ndatagpu=Np;
            ResizeParticlesSizeData(ndatacpu,ndatagpu,Np+count,Np+count,PERIODIC_OVERMEMORYNP,false);
            Timersg->TmStart(TMG_SuPeriodic,false);
          }
          else{
            run=false;
            //-Create new periodic particles duplicating the particles from the list
            //-Crea nuevas particulas periodicas duplicando las particulas de la lista.
            if(TStep==STEP_Verlet){
              cusph::PeriodicDuplicateVerlet(count,Np,DomCells,perinc,listpg.cptr()
                ,Idp_g->ptr(),Code_g->ptr(),Dcell_g->ptr(),Posxy_g->ptr(),Posz_g->ptr()
                ,Velrho_g->ptr(),AG_PTR(SpsTau_g),VelrhoM1_g->ptr());
            }
            if(TStep==STEP_Symplectic){
              if(PosxyPre_g->Active()!=PoszPre_g->Active() || PoszPre_g->Active()!=VelrhoPre_g->Active())
                Run_Exceptioon("Symplectic data is invalid.");
              cusph::PeriodicDuplicateSymplectic(count,Np,DomCells,perinc,listpg.cptr()
                ,Idp_g->ptr(),Code_g->ptr(),Dcell_g->ptr(),Posxy_g->ptr(),Posz_g->ptr()
                ,Velrho_g->ptr(),AG_PTR(SpsTau_g),PosxyPre_g->ptr(),PoszPre_g->ptr(),VelrhoPre_g->ptr());
            }
            if(UseNormals){
              cusph::PeriodicDuplicateNormals(count,Np,listpg.cptr()
                ,BoundNor_g->ptr(),AG_PTR(MotionVel_g));
            }
            //-Update the total number of particles.
            Np+=count;
            //-Update number of new periodic particles.
            if(!ctype)NpbPer+=count;
            else NpfPer+=count;
          }
        }
      }
    }
  }
  Timersg->TmStop(TMG_SuPeriodic,true);
  Check_CudaErroor("Failed in creation of periodic particles.");
}

//==============================================================================
/// Executes divide of particles in cells.
/// Ejecuta divide de particulas en celdas.
//==============================================================================
void JSphGpuSingle::RunCellDivide(bool updateperiodic){
  //JDebugSphGpu::SaveVtk("_DG_Divide_Pre.vtk",Nstep,0,Np,"all",this);
  DivData=DivDataGpuNull();
  //-Creates new periodic particles and marks the old ones to be ignored.
  //-Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
  if(updateperiodic && PeriActive)RunPeriodic();

  //-Initiates Divide process.
  CellDivSingle->Divide(Npb,Np-Npb-NpbPer-NpfPer,NpbPer,NpfPer
    ,BoundChanged,Dcell_g->cptr(),Code_g->cptr(),Posxy_g->cptr()
    ,Posz_g->cptr(),Idp_g->cptr(),Timersg);
  DivData=CellDivSingle->GetCellDivData();

  //-Sorts particle data. | Ordena datos de particulas.
  Timersg->TmStart(TMG_NlSortData,false);
  {
    //-Creates auxiliary arrays for sorted data.
    aguint     idpg   ("-",Arrays_Gpu,true);
    agtypecode codeg  ("-",Arrays_Gpu,true);
    aguint     dcellg ("-",Arrays_Gpu,true);
    agdouble2  posxyg ("-",Arrays_Gpu,true);
    agdouble   poszg  ("-",Arrays_Gpu,true);
    agfloat4   velrhog("-",Arrays_Gpu,true);
    //-Sorts data in auxiliary arrays.
    CellDivSingle->SortBasicArrays(Idp_g->cptr(),Code_g->cptr()
      ,Dcell_g->cptr(),Posxy_g->cptr(),Posz_g->cptr(),Velrho_g->cptr()
      ,idpg.ptr(),codeg.ptr(),dcellg.ptr()
      ,posxyg.ptr(),poszg.ptr(),velrhog.ptr());
    //-Swap pointers with data and auxiliary arrays are automatically freed.
    Idp_g   ->SwapPtr(&idpg);
    Code_g  ->SwapPtr(&codeg);
    Dcell_g ->SwapPtr(&dcellg);
    Posxy_g ->SwapPtr(&posxyg);
    Posz_g  ->SwapPtr(&poszg);
    Velrho_g->SwapPtr(&velrhog);
  }
  if(TStep==STEP_Verlet){
    agfloat4 auxg("-",Arrays_Gpu,true);
    CellDivSingle->SortDataArrays(VelrhoM1_g->cptr(),auxg.ptr());
    VelrhoM1_g->SwapPtr(&auxg);
  }
  else if(TStep==STEP_Symplectic && (PosxyPre_g->Active() || PoszPre_g->Active() || VelrhoPre_g->Active())){ //-In reality, only necessary in the corrector not the predictor step??? | En realidad solo es necesario en el divide del corrector, no en el predictor??? 
    if(!PosxyPre_g->Active() || !PoszPre_g->Active() || !VelrhoPre_g->Active())
      Run_Exceptioon("Symplectic data is invalid.") ;
    agdouble2 posxyg ("-",Arrays_Gpu,true);
    agdouble  poszg  ("-",Arrays_Gpu,true);
    agfloat4  velrhog("-",Arrays_Gpu,true);
    CellDivSingle->SortDataArrays(PosxyPre_g->cptr(),PoszPre_g->cptr()
      ,VelrhoPre_g->cptr(),posxyg.ptr(),poszg.ptr(),velrhog.ptr());
    PosxyPre_g ->SwapPtr(&posxyg);
    PoszPre_g  ->SwapPtr(&poszg);
    VelrhoPre_g->SwapPtr(&velrhog);
  }
  if(TVisco==VISCO_LaminarSPS){
    agsymatrix3f spstaug("-",Arrays_Gpu,true);
    CellDivSingle->SortDataArrays(SpsTau_g->cptr(),spstaug.ptr());
    SpsTau_g->SwapPtr(&spstaug);
  }
  if(UseNormals){
    agfloat3 auxg("-",Arrays_Gpu,true);
    CellDivSingle->SortDataArrays(BoundNor_g->cptr(),auxg.ptr());
    BoundNor_g->SwapPtr(&auxg);
    if(MotionVel_g){
      CellDivSingle->SortDataArrays(MotionVel_g->cptr(),auxg.ptr());
      MotionVel_g->SwapPtr(&auxg);
    }
  }

  //-Collect divide data. | Recupera datos del divide.
  Np=CellDivSingle->GetNpFinal();
  Npb=CellDivSingle->GetNpbFinal();
  NpbOk=Npb-CellDivSingle->GetNpbIgnore();

  //-Update PosCell_g[] according to current position of particles.
  cusphs::UpdatePosCell(Np,Map_PosMin,PosCellSize,Posxy_g->cptr()
    ,Posz_g->cptr(),PosCell_g->ptr(),NULL);

  //-Manages excluded particles fixed, moving and floating before aborting the execution.
  if(CellDivSingle->GetNpbOut())AbortBoundOut();

  //-Collects index of moving and floating particles.
  if(CaseNmoving || CaseNfloat){
    const unsigned np=(CaseNmoving? Npb: 0) + (CaseNfloat? Np-Npb: 0);
    const unsigned pini=(CaseNmoving? 0: Npb);
    cusph::CalcRidp(PeriActive!=0,np,pini,CaseNfixed,CaseNbound
      ,Code_g->cptr(),Idp_g->cptr(),RidpMotg);
  }
  Timersg->TmStop(TMG_NlSortData,true);

  //-Control of excluded particles (only fluid because excluded boundary were checked before).
  if(CellDivSingle->GetNpfOut()){
    Timersg->TmStart(TMG_NlOutCheck,false);
    SaveFluidOut();
    Timersg->TmStop(TMG_NlOutCheck,true);
  }
  BoundChanged=false;
}

//==============================================================================
/// Manages excluded particles fixed, moving and floating before aborting the execution.
/// Gestiona particulas excluidas fixed, moving y floating antes de abortar la ejecucion.
//==============================================================================
void JSphGpuSingle::AbortBoundOut(){
  const unsigned nboundout=CellDivSingle->GetNpbOut();
  //-Get data of excluded boundary particles.
  unsigned nfilter=0;
  ParticlesDataDown(nboundout,Np,true,false,NULL,nfilter);
  //-Shows excluded particles information and aborts execution.
  JSph::AbortBoundOut(Log,nboundout,Idp_c->cptr(),AuxPos_c->cptr()
    ,AuxVel_c->cptr(),AuxRho_c->cptr(),Code_c->cptr());
}

//==============================================================================
/// Manages the excluded fluid particles.
/// Gestiona las particulas fluid excluidas.
//==============================================================================
void JSphGpuSingle::SaveFluidOut(){
  const unsigned npfout=CellDivSingle->GetNpfOut();
  //-Get data of excluded fluid particles.
  unsigned nfilter=0;
  ParticlesDataDown(npfout,Np,true,false,NULL,nfilter);
  //-Stores new excluded particles until recordering next PART.
  AddParticlesOut(npfout,Idp_c->cptr(),AuxPos_c->cptr()
    ,AuxVel_c->cptr(),AuxRho_c->cptr(),Code_c->cptr());
}

//==============================================================================
/// Interaction for force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle::Interaction_Forces(TpInterStep interstep){
  //-Boundary correction for mDBC.
  if(TBoundary==BC_MDBC && (MdbcCorrector || interstep!=INTERSTEP_SymCorrector)){
    MdbcBoundCorrection(); 
  }
  InterStep=interstep;
  PreInteraction_Forces();

  float3* dengradcorr=NULL;

  Timersg->TmStart(TMG_CfForces,true);
  const bool lamsps=(TVisco==VISCO_LaminarSPS);
  unsigned bsfluid=BlockSizes.forcesfluid;
  unsigned bsbound=BlockSizes.forcesbound;

  //-Interaction Fluid-Fluid/Bound & Bound-Fluid.
  const StInterParmsg parms=StrInterParmsg(Simulate2D
    ,Symmetry //<vs_syymmetry>
    ,TKernel,FtMode
    ,lamsps,TDensity,ShiftingMode
    ,Visco*ViscoBoundFactor,Visco
    ,bsbound,bsfluid,Np,Npb,NpbOk
    ,0,Nstep,DivData,Dcell_g->cptr()
    ,Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr()
    ,Velrho_g->cptr(),Idp_g->cptr(),Code_g->cptr()
    ,FtoMasspg,AG_CPTR(SpsTau_g),dengradcorr
    ,ViscDt_g->ptr(),Ar_g->ptr(),Ace_g->ptr(),AG_PTR(Delta_g)
    ,AG_PTR(SpsGradvel_g)
    ,AG_PTR(ShiftPosfs_g)
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
  if(lamsps)cusph::ComputeSpsTau(Np,Npb,SpsSmag,SpsBlin,Velrho_g->cptr()
    ,SpsGradvel_g->cptr(),SpsTau_g->ptr());

  //-Add Delta-SPH correction to Ar_g[].
  if(AG_CPTR(Delta_g))cusph::AddDelta(Np-Npb,Delta_g->cptr()+Npb,Ar_g->ptr()+Npb);

  cudaDeviceSynchronize();
  Check_CudaErroor("Failed while executing kernels of interaction.");

  //-Calculates maximum value of ViscDt.
  if(Np)ViscDtMax=cusph::ReduMaxFloat(Np,0,ViscDt_g->ptr(),CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np)));
  //-Calculates maximum value of Ace (periodic particles are ignored). ViscDtg is used like auxiliary memory.
  AceMax=ComputeAceMax(ViscDt_g->ptr()); 

  Timersg->TmStop(TMG_CfForces,true);
  Check_CudaErroor("Failed in reduction of viscdt.");
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphGpuSingle::MdbcBoundCorrection(){
  Timersg->TmStart(TMG_CfPreForces,false);
  const unsigned n=(UseNormalsFt? Np: NpbOk);
  cusph::Interaction_MdbcCorrection(TKernel,Simulate2D,SlipMode,MdbcFastSingle
    ,n,CaseNbound,MdbcThreshold,DivData,Map_PosMin,Posxy_g->cptr(),Posz_g->cptr()
    ,PosCell_g->cptr(),Code_g->cptr(),Idp_g->cptr(),BoundNor_g->cptr()
    ,AG_CPTR(MotionVel_g),Velrho_g->ptr());
  Timersg->TmStop(TMG_CfPreForces,true);
}


//==============================================================================
/// Returns the maximum value of  (ace.x^2 + ace.y^2 + ace.z^2) from Acec[].
/// Devuelve valor maximo de (ace.x^2 + ace.y^2 + ace.z^2) a partir de Acec[].
//==============================================================================
double JSphGpuSingle::ComputeAceMax(float* auxmemg){
  const bool check=(PeriActive!=0 || InOut!=NULL);
  const unsigned npf=Np-Npb;
  float acemax=0;
  if(check)cusph::ComputeAceMod(npf,Code_g->cptr()+Npb,Ace_g->cptr()+Npb,auxmemg);
  else     cusph::ComputeAceMod(npf,Ace_g->cptr()+Npb,auxmemg);      
  if(npf)acemax=cusph::ReduMaxFloat(npf,0,auxmemg,CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(npf)));
  return(sqrt(double(acemax)));
}

//<vs_ddramp_ini>
//==============================================================================
/// Applies initial DDT ramp.
//==============================================================================
void JSphGpuSingle::RunInitialDDTRamp(){
  if(TimeStep<DDTRamp.x){
    if((Nstep%10)==0){//-DDTkh value is updated every 10 calculation steps.
      if(TimeStep<=DDTRamp.y)DDTkh=KernelSize * float(DDTRamp.z);
      else{
        const double tt=TimeStep-DDTRamp.y;
        const double tr=DDTRamp.x-DDTRamp.y;
        DDTkh=KernelSize * float(((tr-tt)/tr)*(DDTRamp.z-DDTValue)+DDTValue);
      }
      ConstantDataUp(); //-Updates value in constant memory of GPU.
    }
  }
  else{
    if(DDTkh!=DDTkhCte){
      CSP.ddtkh=DDTkh=DDTkhCte;
      ConstantDataUp();
    }
    DDTRamp.x=0;
  }
}//<vs_ddramp_end>

//==============================================================================
/// Particle interaction and update of particle data according to
/// the computed forces using the Verlet time stepping scheme.
///
/// Realiza interaccion y actualizacion de particulas segun las fuerzas 
/// calculadas en la interaccion usando Verlet.
//==============================================================================
double JSphGpuSingle::ComputeStep_Ver(){
  Interaction_Forces(INTERSTEP_Verlet);  //-Interaction.
  const double dt=DtVariable(true);      //-Calculate new dt.
  if(CaseNmoving)CalcMotion(dt);         //-Calculate motion for moving bodies.
  DemDtForce=dt;                         //-For DEM interaction.
  if(Shifting)RunShifting(dt);           //-Shifting.
  ComputeVerlet(dt);                     //-Update particles using Verlet (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);   //-Control of floating bodies.
  PosInteraction_Forces();               //-Free memory used for interaction.
  if(Damping)RunDamping(dt);             //-Aplies Damping.
  if(RelaxZones)RunRelaxZone(dt);        //-Generate waves using RZ.
  return(dt);
}

//==============================================================================
/// Particle interaction and update of particle data according to
/// the computed forces using the Symplectic time stepping scheme.
///
/// Realiza interaccion y actualizacion de particulas segun las fuerzas 
/// calculadas en la interaccion usando Symplectic.
//==============================================================================
double JSphGpuSingle::ComputeStep_Sym(){
  const double dt=SymplecticDtPre;
  if(CaseNmoving)CalcMotion(dt);               //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
  DemDtForce=dt*0.5f;                          //-For DEM interaction.
  Interaction_Forces(INTERSTEP_SymPredictor);  //-Interaction.
  const double dt_p=DtVariable(false);         //-Calculate dt of predictor step.
  if(Shifting)RunShifting(dt*.5);              //-Shifting.
  ComputeSymplecticPre(dt);                    //-Apply Symplectic-Predictor to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt*.5,true);       //-Control of floating bodies.
  PosInteraction_Forces();                     //-Free memory used for interaction.
  //-Corrector
  //-----------
  DemDtForce=dt;                               //-For DEM interaction.
  RunCellDivide(true);
  Interaction_Forces(INTERSTEP_SymCorrector);  //-Interaction.
  const double dt_c=DtVariable(true);          //-Calculate dt of corrector step.
  if(Shifting)RunShifting(dt);                 //-Shifting.
  ComputeSymplecticCorr(dt);                   //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);         //-Control of floating bodies.
  PosInteraction_Forces();                     //-Free memory used for interaction.
  if(Damping)RunDamping(dt);                   //-Aplies Damping.
  if(RelaxZones)RunRelaxZone(dt);              //-Generate waves using RZ.

  SymplecticDtPre=min(dt_p,dt_c);              //-Calculate dt for next ComputeStep.
  return(dt);
}

//==============================================================================
/// Updates information in FtObjs[] copying data from GPU.
/// Actualiza informacion en FtObjs[] copiando los datos en GPU.
//==============================================================================
void JSphGpuSingle::UpdateFtObjs(){
  if(FtCount && FtObjsOutdated){
    tdouble3* fcen=FtoAuxDouble6;
    tfloat3*  fang=FtoAuxFloat15;
    tfloat3*  fvellin=fang+FtCount;
    tfloat3*  fvelang=fvellin+FtCount;
    tfloat3*  facelin=fvelang+FtCount;
    tfloat3*  faceang=facelin+FtCount;
    cudaMemcpy(fcen,FtoCenterg,sizeof(double3)*FtCount,cudaMemcpyDeviceToHost);
    cudaMemcpy(fang,FtoAnglesg,sizeof(float3) *FtCount,cudaMemcpyDeviceToHost);
    cudaMemcpy(fvellin,FtoVelAceg,sizeof(float3)*FtCount*4,cudaMemcpyDeviceToHost);
    for(unsigned cf=0;cf<FtCount;cf++){
      FtObjs[cf].center=fcen[cf];
      FtObjs[cf].angles=fang[cf];
      FtObjs[cf].fvel  =fvellin[cf];
      FtObjs[cf].fomega=fvelang[cf];
      FtObjs[cf].facelin=facelin[cf];
      FtObjs[cf].faceang=faceang[cf];
    }
  }
  FtObjsOutdated=false;
}

//==============================================================================
/// Applies imposed velocity.
/// Aplica velocidad predefinida.
//==============================================================================
void JSphGpuSingle::FtApplyImposedVel(float3* ftoforcesresg)const{
  tfloat3* ftoforcesresc=NULL;
  for(unsigned cf=0;cf<FtCount;cf++)if(!FtObjs[cf].usechrono && (FtLinearVel[cf]!=NULL || FtAngularVel[cf]!=NULL)){
    const tfloat3 v1=(FtLinearVel [cf]!=NULL? FtLinearVel [cf]->GetValue3f(TimeStep): TFloat3(FLT_MAX));
    const tfloat3 v2=(FtAngularVel[cf]!=NULL? FtAngularVel[cf]->GetValue3f(TimeStep): TFloat3(FLT_MAX));
    if(!ftoforcesresc && (v1!=TFloat3(FLT_MAX) || v2!=TFloat3(FLT_MAX))){
      //-Copies data on GPU memory to CPU memory.
      ftoforcesresc=FtoAuxFloat15;
      cudaMemcpy(ftoforcesresc,ftoforcesresg,sizeof(tfloat3)*FtCount*2,cudaMemcpyDeviceToHost);
    }
    unsigned cfpos=cf*2+1;
    if(v1.x!=FLT_MAX)ftoforcesresc[cfpos].x=v1.x;
    if(v1.y!=FLT_MAX)ftoforcesresc[cfpos].y=v1.y;
    if(v1.z!=FLT_MAX)ftoforcesresc[cfpos].z=v1.z;
    cfpos--;
    if(v2.x!=FLT_MAX)ftoforcesresc[cfpos].x=v2.x;
    if(v2.y!=FLT_MAX)ftoforcesresc[cfpos].y=v2.y;
    if(v2.z!=FLT_MAX)ftoforcesresc[cfpos].z=v2.z;
  }
  //-Updates data on GPU memory.
  if(ftoforcesresc!=NULL){
    cudaMemcpy(ftoforcesresg,ftoforcesresc,sizeof(tfloat3)*FtCount*2,cudaMemcpyHostToDevice);
  }
}

//==============================================================================
/// Process floating objects.
/// Procesa floating objects.
//==============================================================================
void JSphGpuSingle::RunFloating(double dt,bool predictor){
  Timersg->TmStart(TMG_SuFloating,false);
  const bool saveftmot=(!predictor && (PartFloatSave && TimeOutExtraCheck(TimeStep+dt)) );
  const bool saveftvalues=(!predictor && (TimeStep+dt>=TimePartNext || saveftmot));
  const bool ftpaused=(TimeStep<FtPause);

  if(!ftpaused || saveftvalues){
    //-Adds external forces (ForcePoints, Moorings, external file) to FtoForces[].
    if(ForcePoints!=NULL || FtLinearForce!=NULL){
      StFtoForces* ftoforces=(StFtoForces*)FtoAuxFloat15;
      memset(ftoforces,0,sizeof(StFtoForces)*FtCount);
      //-Loads sum of linear and angular forces from ForcePoints and Moorings.
      if(ForcePoints)ForcePoints->GetFtForcesSum(ftoforces);
      //-Adds the external forces.
      if(FtLinearForce!=NULL){
        for(unsigned cf=0;cf<FtCount;cf++){
          ftoforces[cf].face     =ftoforces[cf].face     +GetFtExternalForceLin(cf,TimeStep);
          ftoforces[cf].fomegaace=ftoforces[cf].fomegaace+GetFtExternalForceAng(cf,TimeStep);
        }
      }
      //-Copies data to GPU memory.
      cudaMemcpy(FtoForcesg,ftoforces,sizeof(StFtoForces)*FtCount,cudaMemcpyHostToDevice);
      //-Saves sum of external forces applied to floating body.
      if(saveftvalues)for(unsigned cf=0;cf<FtCount;cf++){
        FtObjs[cf].extforcelin=ftoforces[cf].face;
        FtObjs[cf].extforceang=ftoforces[cf].fomegaace;
      }
    }
    else{
      //-Initialises forces of floatings when no external forces are applied.
      cudaMemset(FtoForcesg,0,sizeof(StFtoForces)*FtCount);
      //-Saves sum of external forces applied to floating body.
      if(saveftvalues)for(unsigned cf=0;cf<FtCount;cf++){
        FtObjs[cf].extforcelin=TFloat3(0);
        FtObjs[cf].extforceang=TFloat3(0);
      }
    }

    //-Calculate forces summation (face,fomegaace) starting from floating particles and add in FtoForcesg[].
    cusph::FtCalcForcesSum(PeriActive!=0,FtCount,FtoDatpg,FtoCenterg,RidpMotg
      ,Posxy_g->cptr(),Posz_g->cptr(),Ace_g->cptr(),FtoForcesg);
    //-Saves sum of fluid forces applied to floating body.
    if(saveftvalues){
      StFtoForces* ftoforces=(StFtoForces*)FtoAuxFloat15;
      cudaMemcpy(ftoforces,FtoForcesg,sizeof(StFtoForces)*FtCount,cudaMemcpyDeviceToHost);
      for(unsigned cf=0;cf<FtCount;cf++){
        FtObjs[cf].fluforcelin=ftoforces[cf].face-FtObjs[cf].extforcelin;
        FtObjs[cf].fluforceang=ftoforces[cf].fomegaace-FtObjs[cf].extforceang;
      }    
    }

    //-Computes final acceleration from particles and from external forces in FtoForcesg[].
    cusph::FtCalcForces(FtCount,Gravity,FtoMassg,FtoAnglesg,FtoInertiaini8g,FtoInertiaini1g,FtoForcesg);
    //-Saves acceleration before constraints (includes external forces, gravity and rotated inertia tensor)
    if(saveftvalues){
      StFtoForces* ftoforces=(StFtoForces*)FtoAuxFloat15;
      cudaMemcpy(ftoforces,FtoForcesg,sizeof(StFtoForces)*FtCount,cudaMemcpyDeviceToHost);
      for(unsigned cf=0;cf<FtCount;cf++){
        FtObjs[cf].preacelin=ftoforces[cf].face;
        FtObjs[cf].preaceang=ftoforces[cf].fomegaace;
      }    
    }
  }

  if(!ftpaused){//-Operator >= is used because when FtPause=0 in symplectic-predictor, code would not enter here. | Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
    //-Calculate data to update floatings / Calcula datos para actualizar floatings.
    cusph::FtCalcForcesRes(FtCount,Simulate2D,dt,FtoVelAceg,FtoCenterg,FtoForcesg,FtoForcesResg,FtoCenterResg);
    //-Applies imposed velocity.
    if(FtLinearVel!=NULL)FtApplyImposedVel(FtoForcesResg);
    //-Applies motion constraints.
    if(FtConstraints)cusph::FtApplyConstraints(FtCount,FtoConstraintsg,FtoForcesg,FtoForcesResg);
    
    //-Saves face and fomegace for debug.
    if(SaveFtAce){
      StFtoForces* ftoforces=(StFtoForces*)FtoAuxFloat15;
      cudaMemcpy(ftoforces,FtoForcesg,sizeof(tfloat3)*FtCount*2,cudaMemcpyDeviceToHost);
      SaveFtAceFun(dt,predictor,ftoforces);
    }

    //-Run floating with Chrono library.
    if(ChronoObjects){      
      Timersg->TmStop(TMG_SuFloating,true);
      Timersg->TmStart(TMG_SuChrono,false);
      //-Export data / Exporta datos.
      tfloat3* ftoforces=FtoAuxFloat15;
      cudaMemcpy(ftoforces,FtoForcesg,sizeof(tfloat3)*FtCount*2,cudaMemcpyDeviceToHost);
      for(unsigned cf=0;cf<FtCount;cf++)if(FtObjs[cf].usechrono){
        ChronoObjects->SetFtData(FtObjs[cf].mkbound,ftoforces[cf*2],ftoforces[cf*2+1]);
      }
      //-Applies the external velocities to each floating body of Chrono.
      if(FtLinearVel!=NULL)ChronoFtApplyImposedVel();
      //-Calculate data using Chrono / Calcula datos usando Chrono.
      ChronoObjects->RunChrono(Nstep,TimeStep,dt,predictor);
      //-Load calculated data by Chrono / Carga datos calculados por Chrono.
      tdouble3* ftocenter=FtoAuxDouble6;
      cudaMemcpy(ftocenter,FtoCenterResg,sizeof(tdouble3)*FtCount  ,cudaMemcpyDeviceToHost);//-Necesario para cargar datos de floatings sin chrono.
      cudaMemcpy(ftoforces,FtoForcesResg,sizeof(tfloat3) *FtCount*2,cudaMemcpyDeviceToHost);//-Necesario para cargar datos de floatings sin chrono.
      for(unsigned cf=0;cf<FtCount;cf++)if(FtObjs[cf].usechrono)ChronoObjects->GetFtData(FtObjs[cf].mkbound,ftocenter[cf],ftoforces[cf*2+1],ftoforces[cf*2]);
      cudaMemcpy(FtoCenterResg,ftocenter,sizeof(tdouble3)*FtCount  ,cudaMemcpyHostToDevice);
      cudaMemcpy(FtoForcesResg,ftoforces,sizeof(float3)  *FtCount*2,cudaMemcpyHostToDevice);
      Timersg->TmStop(TMG_SuChrono,false);
      Timersg->TmStart(TMG_SuFloating,false);
    }

    //-Apply movement around floating objects / Aplica movimiento sobre floatings.
    cusph::FtUpdate(PeriActive!=0,predictor,FtCount,dt,FtoDatpg,FtoForcesResg,FtoCenterResg
      ,RidpMotg,FtoCenterg,FtoAnglesg,FtoVelAceg,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr()
      ,Velrho_g->ptr(),Code_g->ptr());

    //-Stores floating data.
    if(!predictor){
      FtObjsOutdated=true;
      //-Updates floating normals for mDBC.
      if(UseNormalsFt){
        tdouble3* fcen=FtoAuxDouble6;
        tfloat3*  fang=FtoAuxFloat15;
        cudaMemcpy(fcen,FtoCenterg,sizeof(double3)*FtCount,cudaMemcpyDeviceToHost);
        cudaMemcpy(fang,FtoAnglesg,sizeof(float3) *FtCount,cudaMemcpyDeviceToHost);
        for(unsigned cf=0;cf<FtCount;cf++){
          const StFloatingData fobj=FtObjs[cf];
          FtObjs[cf].center=fcen[cf];
          FtObjs[cf].angles=fang[cf];
          const tdouble3 dang=ToTDouble3(FtObjs[cf].angles-fobj.angles)*TODEG;
          const tdouble3 cen=FtObjs[cf].center;
          JMatrix4d mat;
          mat.Move(cen);
          mat.Rotate(dang);
          mat.Move(fobj.center*-1);
          cusph::FtNormalsUpdate(fobj.count,fobj.begin-CaseNfixed,mat.GetMatrix()
            ,RidpMotg,AG_PTR(BoundNor_g));
        }
      }
    }
  }
  //-Saves current floating data for output.
  if(saveftvalues){
    UpdateFtObjs(); //-Updates floating information on CPU memory.
    if(PartFloatSave)PartFloatSave->SetFtData(Part,TimeStep+dt,Nstep+1,FtObjs,ForcePoints);
  }
  Timersg->TmStop(TMG_SuFloating,true);

  //-Update data of points in FtForces and calculates motion data of affected floatings.
  if(!predictor && ForcePoints){
    Timersg->TmStart(TMG_SuMoorings,false);
    UpdateFtObjs(); //-Updates floating information on CPU memory.
    ForcePoints->UpdatePoints(TimeStep,dt,ftpaused,FtObjs);
    if(Moorings)Moorings->ComputeForces(Nstep,TimeStep,dt,ForcePoints);
    ForcePoints->ComputeForcesSum();
    Timersg->TmStop(TMG_SuMoorings,true);
  }
}

//==============================================================================
/// Runs calculations in configured gauges.
/// Ejecuta calculos en las posiciones de medida configuradas.
//==============================================================================
void JSphGpuSingle::RunGaugeSystem(double timestep,bool saveinput){
  if(!Nstep || GaugeSystem->GetCount()){
    Timersg->TmStart(TMG_SuGauges,false);
    //const bool svpart=(TimeStep>=TimePartNext);
    GaugeSystem->CalculeGpu(timestep,DivData
      ,NpbOk,Npb,Np,Posxy_g->cptr(),Posz_g->cptr(),Code_g->cptr()
      ,Idp_g->cptr(),Velrho_g->cptr(),saveinput);
    Timersg->TmStop(TMG_SuGauges,true);
  }
}

//==============================================================================
/// Compute PIPS information of current particles.
/// Calcula datos de PIPS de particulas actuales.
//==============================================================================
void JSphGpuSingle::ComputePips(bool run){
  if(run || DsPips->CheckRun(Nstep)){
    TimerSim.Stop();
    const double timesim=TimerSim.GetElapsedTimeD()/1000.;
    const unsigned sauxmemg=Arrays_Gpu->GetArraySize();
    aguint auxmemg("auxmemg",Arrays_Gpu,true);
    DsPips->ComputeGpu(Nstep,TimeStep,timesim,Np,Npb,NpbOk
      ,DivData,Dcell_g->cptr(),PosCell_g->cptr(),sauxmemg,auxmemg.ptr());
  }
}

//==============================================================================
/// Initialises execution of simulation.
/// Inicia ejecucion de simulacion.
//==============================================================================
void JSphGpuSingle::Run(std::string appname,const JSphCfgRun* cfg,JLog2* log){
  if(!cfg || !log)return;
  AppName=appname; Log=log; CfgRun=cfg;

  //-Selection of GPU.
  const int gpuid=SelecDevice(cfg->GpuId);

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

  //-Initialisation of execution variables.
  InitRunGpu();
  RunGaugeSystem(TimeStep,true);
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

  //-Main Loop.
  //------------
  JTimeControl tc("30,60,300,600");//-Shows information at 0.5, 1, 5 y 10 minutes (before first PART).
  bool minfluidstopped=false;
  TimerSim.Start();
  TimerPart.Start();
  Log->Print(string("\n[Initialising simulation (")+RunCode+")  "+fun::GetDateTime()+"]");
  if(DsPips)ComputePips(true);
  PrintHeadPart();
  while(TimeStep<TimeMax){
    InterStep=(TStep==STEP_Symplectic? INTERSTEP_SymPredictor: INTERSTEP_Verlet);
    if(ViscoTime)Visco=ViscoTime->GetVisco(float(TimeStep));
    if(DDTRamp.x)RunInitialDDTRamp(); //<vs_ddramp>
    const double stepdt=ComputeStep();
    RunGaugeSystem(TimeStep+stepdt);
    if(CaseNmoving)RunMotion(stepdt);
    if(InOut)InOutComputeStep(stepdt);
    else RunCellDivide(true);
    TimeStep+=stepdt;
    LastDt=stepdt;
    Nstep++;
    //-Save extra PART data.
    if(TimeOutExtraCheck(TimeStep)){
      if(PartFloatSave)PartFloatSave->AddDataExtra(Part,TimeStep,Nstep);
      if(PartMotionSave)PartMotionSave->AddDataExtraGpu(Part,TimeStep,Nstep,Np
        ,Posxy_g->cptr(),Posz_g->cptr(),RidpMotg);
      TimeOutExtraUpdate(TimeStep);
    }
    //-Check minimum fluid allowed.
    const unsigned npfnormal=Np-NpbPer-NpfPer-CaseNbound;
    minfluidstopped=(npfnormal<NpfMinimum || !Np);
    //-Save main PART data.
    if(TimeStep>=TimePartNext || minfluidstopped){
      if(minfluidstopped){
        Log->PrintWarning(fun::PrintStr("The minimum number of fluid particles (%s) was reached.",KINT(NpfMinimum)));
        TimeMax=TimeStep;
      }
      SaveData();
      Part++;
      PartNstep=Nstep;
      TimeStepM1=TimeStep;
      TimePartNext=(SvAllSteps? TimeStep: OutputTime->GetNextTime(TimeStep));
      TimerPart.Start();
    }
    UpdateMaxValues();
    const bool laststep=(TimeStep>=TimeMax || (NstepsBreak && Nstep>=NstepsBreak));
    if(DsPips)ComputePips(laststep);
    if(Part<=PartIni+1 && tc.CheckTime())Log->Print(string("  ")
      +tc.GetInfoFinish((TimeStep-TimeStepIni)/(TimeMax-TimeStepIni)));
    //-Terminates the simulation according to NstepsBreak (for debugging).
    if(NstepsBreak && Nstep>=NstepsBreak)break;
  }
  TimerSim.Stop(); TimerTot.Stop();

  //-End of Simulation.
  //--------------------
  FinishRun(minfluidstopped);
}

//==============================================================================
/// Generates files with output data.
/// Genera los ficheros de salida de datos.
//==============================================================================
void JSphGpuSingle::SaveData(){
  //-Retrieve floating object data from the GPU. | Recupera datos de floatings en GPU.
  if(FtCount){
    Timersg->TmStart(TMG_SuDownData,false);
    UpdateFtObjs();
    Timersg->TmStop(TMG_SuDownData,true);
  }
  const bool save=(SvData!=SDAT_None && SvData!=SDAT_Info);
  const unsigned npnormal=Np-NpbPer-NpfPer; //-Subtracts the periodic particles if they exist. | Resta las periodicas si las hubiera.
  unsigned npsave=npnormal;
  //-Retrieves particle data from the GPU. | Recupera datos de particulas en GPU.
  if(save){
    Timersg->TmStart(TMG_SuDownData,false);
    //-Prepare filter for output particles data. //<vs_outpaarts>
    agbyte filterg("filterg",Arrays_Gpu,false);
    //<vs_outpaarts_ini>
    const bool svextra=(SvExtraDataBi4 && SvExtraDataBi4->CheckSave(Part));
    if(OutputParts && OutputParts->CheckFilters(Part) && !svextra){
      filterg.Reserve();
      OutputParts->UpdateFtPos(FtCount,FtObjs);
      OutputParts->ComputeFilterGpu(Np,Posxy_g->cptr(),Posz_g->cptr()
        ,Code_g->cptr(),filterg.ptr());
    }//<vs_outpaarts_end>
    //-Obtain output particles data.
    unsigned npfilterdel=0;
    const unsigned npsel=ParticlesDataDown(Np,0,false,PeriActive!=0
      ,filterg.ptr(),npfilterdel);
    if(npsel+npfilterdel!=npnormal)Run_Exceptioon("The number of particles is invalid.");
    npsave=npsel;
    Timersg->TmStop(TMG_SuDownData,true);
  }

  Timersg->TmStart(TMG_SuSavePart,false);
  //-Saves main motion reference data from particles (moving and floating bodies).
  if(PartMotionSave){
    PartMotionSave->SaveDataMainGpu(Part,TimeStep,Nstep,Np,Posxy_g->cptr()
      ,Posz_g->cptr(),RidpMotg);
    PartMotionSave->SaveDataExtra();
  }

  //-Collects additional information.
  StInfoPartPlus infoplus;
  if(SvData&SDAT_Info){
    infoplus.gpudata=true;
    infoplus.SetNct(CellDivSingle->GetNct(),CellDivSingle->GetSizeNct());
    infoplus.SetNp(Np,GpuParticlesSize,npnormal,npsave);
    if(InOut)infoplus.npnew=InOut->GetNewNpPart();
    infoplus.SetNpExtra(NpbOk,Npb-NpbOk,Np-Npb,NpbPer,NpfPer);
    infoplus.memorycpualloc=GetAllocMemoryCpu();
    infoplus.memorynctalloc=infoplus.memorynctused=GetMemoryGpuNct();
    infoplus.memorynpalloc=infoplus.memorynpused=GetMemoryGpuNp();
    TimerSim.Stop();
    infoplus.timesim=TimerSim.GetElapsedTimeD()/1000.;
  }
  //-Obtains current domain limits.
  const tdouble6 vdom=CellDivSingle->GetDomainLimitsMinMax();
  //-Stores particle data. | Graba datos de particulas.
  JDataArrays arrays;
  AddBasicArrays(arrays,npsave,AuxPos_c->cptr(),Idp_c->cptr(),AuxVel_c->cptr(),AuxRho_c->cptr());
  JSph::SaveData(npsave,arrays,1,&vdom,infoplus);
  //-Save VTK file with current boundary normals (for debug).
  if(UseNormals && SvNormals)SaveVtkNormalsGpu(DirVtkOut+"Normals.vtk",Part
    ,npsave,Npb,Posxy_g->cptr(),Posz_g->cptr(),Idp_g->cptr(),BoundNor_g->cptr());
  //-Save extra data.
  if(SvExtraDataBi4)SaveExtraData();
  Timersg->TmStop(TMG_SuSavePart,true);
}

//==============================================================================
/// Displays and stores final summary of the execution.
/// Muestra y graba resumen final de ejecucion.
//==============================================================================
void JSphGpuSingle::SaveExtraData(){
  const bool svextra=(BoundNor_g!=NULL);
  if(svextra && SvExtraDataBi4->CheckSave(Part)){
    SvExtraDataBi4->InitPartData(Part,TimeStep,Nstep);
    //-Saves normals of mDBC.
    if(BoundNor_g){
      acfloat3* nor_c=AuxVel_c;
      const unsigned nsize=(UseNormalsFt? Np: Npb);
      Idp_g->CuCopyToHost(Idp_c,nsize);
      BoundNor_g->CuCopyToHost(nor_c,nsize);
      if(PeriActive)Code_g->CuCopyToHost(Code_c,nsize);
      SvExtraDataBi4->AddNormals(UseNormalsFt,Np,Npb,Idp_c->cptr()
        ,(PeriActive? Code_c->cptr(): NULL),nor_c->cptr());
    }
    //-Saves file.
    SvExtraDataBi4->SavePartData();
  }
}

//==============================================================================
/// Displays and stores final summary of the execution.
/// Muestra y graba resumen final de ejecucion.
//==============================================================================
void JSphGpuSingle::FinishRun(bool stop){
  const float tsim=TimerSim.GetElapsedTimeF()/1000.f;
  const float ttot=TimerTot.GetElapsedTimeF()/1000.f;
  JSph::ShowResume(stop,tsim,ttot,true,"");
  Log->Print(" ");
  string hinfo,dinfo;
  if(SvTimers){
    Timersg->ShowTimes("[GPU Timers]",Log);
    Timersg->GetTimersInfo(hinfo,dinfo);
  }
  if(SvRes)SaveRes(tsim,ttot,hinfo,dinfo);
  if(SvData&SDAT_Info)SaveRunPartsCsvFinal();
  Log->PrintFilesList();
  Log->PrintWarningList();
  VisuRefs();
}

