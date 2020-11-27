//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "JArraysGpu.h"
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
#include "JSphGpu_NN_ker.h"  //<vs_non-Newtonian>
#include "JDsGaugeSystem.h"
#include "JSphInOut.h"
#include "JFtMotionSave.h" //<vs_ftmottionsv>  
#include "JLinearValue.h"
#include "JDataArrays.h"
#include "JDebugSphGpu.h"
#include "JSphShifting.h"
#include "JDsPips.h"

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
  llong s=CellDivSingle->GetAllocMemoryGpuNct();
  return(CellDivSingle->GetAllocMemoryGpuNct());
}

//==============================================================================
/// Updates the maximum values of memory, particles and cells.
/// Actualiza los valores maximos de memory, particles y cells.
//==============================================================================
void JSphGpuSingle::UpdateMaxValues(){
  const llong mcpu=GetAllocMemoryCpu();
  const llong mgpu=GetAllocMemoryGpu();
  MaxNumbers.memcpu=max(MaxNumbers.memcpu,mcpu);
  MaxNumbers.memgpu=max(MaxNumbers.memgpu,mgpu);
  MaxNumbers.particles=max(MaxNumbers.particles,Np);
  if(CellDivSingle)MaxNumbers.cells=max(MaxNumbers.cells,CellDivSingle->GetNct());
}

//==============================================================================
/// Loads the configuration of the execution.
/// Carga la configuracion de ejecucion.
//==============================================================================
void JSphGpuSingle::LoadConfig(JSphCfgRun *cfg){
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
  //-Computes the number of particles.
  Np=PartsLoaded->GetCount(); Npb=CaseNpb; NpbOk=Npb;
  //-Allocates memory for arrays with fixed size (motion and floating bodies).
  AllocGpuMemoryFixed();
  //-Allocates GPU memory for particles.
  AllocGpuMemoryParticles(Np,0);
  //-Allocates memory on the CPU.
  AllocCpuMemoryFixed();
  AllocCpuMemoryParticles(Np);

  //-Copies particle data.
  memcpy(AuxPos,PartsLoaded->GetPos(),sizeof(tdouble3)*Np);
  memcpy(Idp,PartsLoaded->GetIdp(),sizeof(unsigned)*Np);
  memcpy(Velrhop,PartsLoaded->GetVelRhop(),sizeof(tfloat4)*Np);

  //-Computes radius of floating bodies.
  if(CaseNfloat && PeriActive!=0 && !PartBegin)CalcFloatingRadius(Np,AuxPos,Idp);
  //-Configures floating motion data storage with high frequency. //<vs_ftmottionsv>  
  if(FtMotSave)ConfigFtMotionSave(Np,AuxPos,Idp);                 //<vs_ftmottionsv>  

  //-Configures Multi-Layer Pistons according particles. | Configura pistones Multi-Layer segun particulas.
  if(MLPistons)MLPistons->PreparePiston(Dp,Np,Idp,AuxPos);

  //-Loads Code of the particles.
  LoadCodeParticles(Np,Idp,Code);
  //<vs_non-Newtonian>
  LoadMultiphaseData(Np, Idp, Code, Velrhop, AuxNN);

  //-Load normals for boundary particles (fixed and moving).
  tfloat3 *boundnormal=NULL;
  if(UseNormals){
    boundnormal=new tfloat3[Np];
    memset(boundnormal,0,sizeof(tfloat3)*Np);
    LoadBoundNormals(Np,Npb,Idp,Code,boundnormal);
  }

  //-Creates PartsInit object with initial particle data for automatic configurations.
  CreatePartsInit(Np,AuxPos,Code);

  //-Runs initialization operations from XML.
  RunInitialize(Np,Npb,AuxPos,Idp,Code,Velrhop,boundnormal);
  if(UseNormals)ConfigBoundNormals(Np,Npb,AuxPos,Idp,boundnormal);

  //-Computes MK domain for boundary and fluid particles.
  MkInfo->ComputeMkDomains(Np,AuxPos,Code);

  //-Configure cell division.
  ConfigCellDivision();
  //-Sets local domain of the simulation within Map_Cells and computes DomCellCode.
  //-Establece dominio de simulacion local dentro de Map_Cells y calcula DomCellCode.
  SelecDomain(TUint3(0,0,0),Map_Cells);
  //-Computes inital cell of the particles and checks if there are unexpected excluded particles.
  //-Calcula celda inicial de particulas y comprueba si hay excluidas inesperadas.
  LoadDcellParticles(Np,Code,AuxPos,Dcell);

  //-Uploads particle data on the GPU.
  ReserveBasicArraysGpu();
  for(unsigned p=0;p<Np;p++){ Posxy[p]=TDouble2(AuxPos[p].x,AuxPos[p].y); Posz[p]=AuxPos[p].z; }
  ParticlesDataUp(Np,boundnormal);
  delete[] boundnormal; boundnormal=NULL;
  //-Uploads constants on the GPU.
  ConstantDataUp();

  //-Creates object for Celldiv on the GPU and selects a valid cellmode.
  //-Crea objeto para divide en GPU y selecciona un cellmode valido.
  CellDivSingle=new JCellDivGpuSingle(Stable,FtCount!=0,PeriActive,KernelSize2,PosCellSize
    ,CellMode,Scell,Map_PosMin,Map_PosMax,Map_Cells,CaseNbound,CaseNfixed,CaseNpb,DirOut);
  CellDivSingle->DefineDomain(DomCellCode,DomCelIni,DomCelFin,DomPosMin,DomPosMax);
  ConfigCellDiv((JCellDivGpu*)CellDivSingle);

  ConfigBlockSizes(false,PeriActive!=0);
  ConfigSaveData(0,1,"");

  //-Reorders particles according to cells.
  //-Reordena particulas por celda.
  BoundChanged=true;
  RunCellDivide(true);
}

//==============================================================================
/// Resizes the allocated space for particles on the CPU and the GPU measuring
/// the time spent with TMG_SuResizeNp. At the end updates the division.
///
/// Redimensiona el espacio reservado para particulas en CPU y GPU midiendo el
/// tiempo consumido con TMG_SuResizeNp. Al terminar actualiza el divide.
//==============================================================================
void JSphGpuSingle::ResizeParticlesSize(unsigned newsize,float oversize,bool updatedivide){
  TmgStart(Timers,TMG_SuResizeNp);
  newsize+=(oversize>0? unsigned(oversize*newsize): 0);
  FreeCpuMemoryParticles();
  CellDivSingle->FreeMemoryGpu();
  DivData=DivDataGpuNull();
  ResizeGpuMemoryParticles(newsize);
  AllocCpuMemoryParticles(newsize);
  TmgStop(Timers,TMG_SuResizeNp);
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
  TmgStart(Timers,TMG_SuPeriodic);
  //-Stores the current number of periodic particles.
  //-Guarda numero de periodicas actuales.
  NpfPerM1=NpfPer;
  NpbPerM1=NpbPer;
  //-Marks current periodic particles to be ignored.
  //-Marca periodicas actuales para ignorar.
  cusph::PeriodicIgnore(Np,Codeg);
  //-Creates new periodic particles.
  //-Crea las nuevas periodicas.
  const unsigned npb0=Npb;
  const unsigned npf0=Np-Npb;
  const unsigned np0=Np;
  NpbPer=NpfPer=0;
  BoundChanged=true;
  for(unsigned ctype=0;ctype<2;ctype++){//-0:bound, 1:fluid+floating.
    //-Computes the particle range to be examined (bound and fluid).
    //-Calcula rango de particulas a examinar (bound o fluid).
    const unsigned pini=(ctype? npb0: 0);
    const unsigned num= (ctype? npf0: npb0);
    //-Searches for periodic zones on each axis (X, Y and Z).
    //-Busca periodicas en cada eje (X, Y e Z).
    for(unsigned cper=0;cper<3;cper++)if((cper==0 && PeriX) || (cper==1 && PeriY) || (cper==2 && PeriZ)){
      tdouble3 perinc=(cper==0? PeriXinc: (cper==1? PeriYinc: PeriZinc));
      //-First searches in the list of new periodic particles and then in the initial particle list (necessary for periodic zones in more than one axis)
      //-Primero busca en la lista de periodicas nuevas y despues en la lista inicial de particulas (necesario para periodicas en mas de un eje).
      for(unsigned cblock=0;cblock<2;cblock++){//-0:new periodic particles, 1:original periodic particles
        const unsigned nper=(ctype? NpfPer: NpbPer);  //-number of new periodic particles for the type currently computed (bound or fluid). | Numero de periodicas nuevas del tipo a procesar.
        const unsigned pini2=(cblock? pini: Np-nper);
        const unsigned num2= (cblock? num:  nper);
        //-Repeats search if the available memory was insufficient and had to be increased.
        //-Repite la busqueda si la memoria disponible resulto insuficiente y hubo que aumentarla.
        bool run=true;
        while(run && num2){
          //-Allocates memory to create the periodic particle list.
          //-Reserva memoria para crear lista de particulas periodicas.
          unsigned* listpg=ArraysGpu->ReserveUint();
          unsigned nmax=GpuParticlesSize-1; //-Maximum number of particles that can be included in the list. | Numero maximo de particulas que caben en la lista.
          //-Generates list of new periodic particles
          if(Np>=0x80000000)Run_Exceptioon("The number of particles is too big.");//-Because the last bit is used to mark the reason the new periodical is created. | Porque el ultimo bit se usa para marcar el sentido en que se crea la nueva periodica. 
          unsigned count=cusph::PeriodicMakeList(num2,pini2,Stable,nmax,Map_PosMin,Map_PosMax,perinc,Posxyg,Poszg,Codeg,listpg);
          //-Resizes the memory size for the particles if there is not sufficient space and repeats the serach process.
          //-Redimensiona memoria para particulas si no hay espacio suficiente y repite el proceso de busqueda.
          if(count>nmax || !CheckGpuParticlesSize(count+Np)){
            ArraysGpu->Free(listpg); listpg=NULL;
            TmgStop(Timers,TMG_SuPeriodic);
            ResizeParticlesSize(Np+count,PERIODIC_OVERMEMORYNP,false);
            TmgStart(Timers,TMG_SuPeriodic);
          }
          else{
            run=false;
            //-Create new periodic particles duplicating the particles from the list
            //-Crea nuevas particulas periodicas duplicando las particulas de la lista.
            if(TStep==STEP_Verlet) {
              if(!MultiPhase)     cusph::PeriodicDuplicateVerlet(count,Np,DomCells,perinc,listpg,Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,SpsTaug,VelrhopM1g);  //<vs_non-Newtonian>-pdte
              else                cusphNN::PeriodicDuplicateVerlet(count,Np,DomCells,perinc,listpg,Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,AuxNNg,VelrhopM1g);  //<vs_non-Newtonian>-pdte
            }
            if(TStep==STEP_Symplectic){
              if((PosxyPreg || PoszPreg || VelrhopPreg) && (!PosxyPreg || !PoszPreg || !VelrhopPreg))Run_Exceptioon("Symplectic data is invalid.");
              if(!MultiPhase)  cusph::PeriodicDuplicateSymplectic(count,Np,DomCells,perinc,listpg,Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,SpsTaug,PosxyPreg,PoszPreg,VelrhopPreg);
              else              cusphNN::PeriodicDuplicateSymplectic(count,Np,DomCells,perinc,listpg,Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,AuxNNg,PosxyPreg,PoszPreg,VelrhopPreg);  //<vs_non-Newtonian>
            }
            if(UseNormals)cusph::PeriodicDuplicateNormals(count,Np,listpg,BoundNormalg,MotionVelg);

            //-Frees memory and updates the particle number.
            //-Libera lista y actualiza numero de particulas.
            ArraysGpu->Free(listpg); listpg=NULL;
            Np+=count;
            //-Updated number of new periodic particles.
            //-Actualiza numero de periodicas nuevas.
            if(!ctype)NpbPer+=count;
            else NpfPer+=count;
          }
        }
      }
    }
  }
  Check_CudaErroor("Failed in creation of periodic particles.");
  TmgStop(Timers,TMG_SuPeriodic);
}

//==============================================================================
/// Executes divide of particles in cells.
/// Ejecuta divide de particulas en celdas.
//==============================================================================
void JSphGpuSingle::RunCellDivide(bool updateperiodic){
  DivData=DivDataGpuNull();
  //-Creates new periodic particles and marks the old ones to be ignored.
  //-Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
  if(updateperiodic && PeriActive)RunPeriodic();

  //-Initiates Divide.
  CellDivSingle->Divide(Npb,Np-Npb-NpbPer-NpfPer,NpbPer,NpfPer,BoundChanged,Dcellg,Codeg,Timers,Posxyg,Poszg,Idpg);
  DivData=CellDivSingle->GetCellDivData();

  //-Sorts particle data. | Ordena datos de particulas.
  TmgStart(Timers,TMG_NlSortData);
  {
    unsigned* idpg=ArraysGpu->ReserveUint();
    typecode* codeg=ArraysGpu->ReserveTypeCode();
    unsigned* dcellg=ArraysGpu->ReserveUint();
    double2*  posxyg=ArraysGpu->ReserveDouble2();
    double*   poszg=ArraysGpu->ReserveDouble();
    float4*   velrhopg=ArraysGpu->ReserveFloat4();
    CellDivSingle->SortBasicArrays(Idpg,Codeg,Dcellg,Posxyg,Poszg,Velrhopg,idpg,codeg,dcellg,posxyg,poszg,velrhopg);
    swap(Idpg,idpg);           ArraysGpu->Free(idpg);
    swap(Codeg,codeg);         ArraysGpu->Free(codeg);
    swap(Dcellg,dcellg);       ArraysGpu->Free(dcellg);
    swap(Posxyg,posxyg);       ArraysGpu->Free(posxyg);
    swap(Poszg,poszg);         ArraysGpu->Free(poszg);
    swap(Velrhopg,velrhopg);   ArraysGpu->Free(velrhopg);
  }
  if(TStep==STEP_Verlet){
    float4* velrhopg=ArraysGpu->ReserveFloat4();
    CellDivSingle->SortDataArrays(VelrhopM1g,velrhopg);
    swap(VelrhopM1g,velrhopg);   ArraysGpu->Free(velrhopg);
  }
  else if(TStep==STEP_Symplectic && (PosxyPreg || PoszPreg || VelrhopPreg)){ //-In reality, only necessary in the corrector not the predictor step??? | En realidad solo es necesario en el divide del corrector, no en el predictor??? 
    if(!PosxyPreg || !PoszPreg || !VelrhopPreg)Run_Exceptioon("Symplectic data is invalid.") ;
    double2* posxyg=ArraysGpu->ReserveDouble2();
    double* poszg=ArraysGpu->ReserveDouble();
    float4* velrhopg=ArraysGpu->ReserveFloat4();
    CellDivSingle->SortDataArrays(PosxyPreg,PoszPreg,VelrhopPreg,posxyg,poszg,velrhopg);
    swap(PosxyPreg,posxyg);      ArraysGpu->Free(posxyg);
    swap(PoszPreg,poszg);        ArraysGpu->Free(poszg);
    swap(VelrhopPreg,velrhopg);  ArraysGpu->Free(velrhopg);
  }
  if(!MultiPhase && TVisco==VISCO_LaminarSPS) { //<vs_non-Newtonian>
    tsymatrix3f *spstaug=ArraysGpu->ReserveSymatrix3f();
    CellDivSingle->SortDataArrays(SpsTaug,spstaug);
    swap(SpsTaug,spstaug);  ArraysGpu->Free(spstaug);
  }
  if(MultiPhase) {  //<vs_non-Newtonian>
    float *auxg=ArraysGpu->ReserveFloat();
    CellDivSingle->SortDataArrays(AuxNNg,auxg);
    swap(AuxNNg,auxg);  ArraysGpu->Free(auxg);

    //tsymatrix3f *spstaug = ArraysGpu->ReserveSymatrix3f();
    //CellDivSingle->SortDataArrays(SpsTaug, spstaug);
    //swap(SpsTaug, spstaug);  ArraysGpu->Free(spstaug);

    //tsymatrix3f *d_tensor = ArraysGpu->ReserveSymatrix3f();
    //CellDivSingle->SortDataArrays(D_tensorg, d_tensor);
    //swap(D_tensorg, d_tensor);  ArraysGpu->Free(d_tensor);
  }
  if(UseNormals){
    float3* boundnormalg=ArraysGpu->ReserveFloat3();
    CellDivSingle->SortDataArrays(BoundNormalg,boundnormalg);
    swap(BoundNormalg,boundnormalg); ArraysGpu->Free(boundnormalg);
    if(MotionVelg){
      float3* motionvelg=ArraysGpu->ReserveFloat3();
      CellDivSingle->SortDataArrays(MotionVelg,motionvelg);
      swap(MotionVelg,motionvelg); ArraysGpu->Free(motionvelg);
    }
  }

  //-Collect divide data. | Recupera datos del divide.
  Np=CellDivSingle->GetNpFinal();
  Npb=CellDivSingle->GetNpbFinal();
  NpbOk=Npb-CellDivSingle->GetNpbIgnore();

  //-Update PosCellg[] according to current position of particles.
  cusphs::UpdatePosCell(Np,Map_PosMin,PosCellSize,Posxyg,Poszg,PosCellg,NULL);

  //-Manages excluded particles fixed, moving and floating before aborting the execution.
  if(CellDivSingle->GetNpbOut())AbortBoundOut();

  //-Collect position of floating particles. | Recupera posiciones de floatings.
  if(CaseNfloat)cusph::CalcRidp(PeriActive!=0,Np-Npb,Npb,CaseNpb,CaseNpb+CaseNfloat,Codeg,Idpg,FtRidpg);
  TmgStop(Timers,TMG_NlSortData);

  //-Control of excluded particles (only fluid because excluded boundary are checked before).
  //-Gestion de particulas excluidas (solo fluid porque las boundary excluidas se comprueban antes).
  TmgStart(Timers,TMG_NlOutCheck);
  unsigned npfout=CellDivSingle->GetNpfOut();
  if(npfout){
    ParticlesDataDown(npfout,Np,true,false);
    AddParticlesOut(npfout,Idp,AuxPos,AuxVel,AuxRhop,Code);
  }
  TmgStop(Timers,TMG_NlOutCheck);
  BoundChanged=false;
}

//------------------------------------------------------------------------------
/// Manages excluded particles fixed, moving and floating before aborting the execution.
/// Gestiona particulas excluidas fixed, moving y floating antes de abortar la ejecucion.
//------------------------------------------------------------------------------
void JSphGpuSingle::AbortBoundOut(){
  const unsigned nboundout=CellDivSingle->GetNpbOut();
  //-Get data of excluded boundary particles.
  ParticlesDataDown(nboundout,Np,true,false);
  //-Shows excluded particles information and aborts execution.
  JSph::AbortBoundOut(Log,nboundout,Idp,AuxPos,AuxVel,AuxRhop,Code);
}

//==============================================================================
/// Interaction for force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle::Interaction_Forces(TpInterStep interstep){
  if(TBoundary==BC_MDBC && (MdbcCorrector || interstep!=INTERSTEP_SymCorrector))MdbcBoundCorrection(); //-Boundary correction for mDBC.
  InterStep=interstep;
  PreInteraction_Forces();
  TmgStart(Timers,TMG_CfForces);

  const bool lamsps=(TVisco==VISCO_LaminarSPS);
  unsigned bsfluid=BlockSizes.forcesfluid;
  unsigned bsbound=BlockSizes.forcesbound;

  //-Interaction Fluid-Fluid/Bound & Bound-Fluid.
  const StInterParmsg parms=StrInterParmsg(Simulate2D
    ,Symmetry //<vs_syymmetry>
    ,TKernel,FtMode
    ,lamsps, MultiPhase, TVisco, TVelGrad, TDensity, ShiftingMode    
    ,Visco*ViscoBoundFactor,Visco
    ,bsbound,bsfluid,Np,Npb,NpbOk
    ,0,DivData,Dcellg
    ,Posxyg,Poszg,PosCellg,Velrhopg,Idpg,Codeg
    ,FtoMasspg
    ,ViscDtg, ViscEtaDtg,Arg,Aceg,Deltag
    ,Visco_etag, SpsTaug, SpsGradvelg, D_tensorg, AuxNNg
    ,ShiftPosfsg
    ,NULL,NULL);

  if(!MultiPhase)cusph::Interaction_Forces(parms);
  else         cusphNN::Interaction_ForcesNN(parms);  //<vs_non-Newtonian>

  //-Interaction DEM Floating-Bound & Floating-Floating. //(DEM)
  if(UseDEM)cusph::Interaction_ForcesDem(BlockSizes.forcesdem,CaseNfloat
    ,DivData,Dcellg,FtRidpg,DemDatag,FtoMasspg,float(DemDtForce)
    ,PosCellg,Velrhopg,Codeg,Idpg,ViscDtg,Aceg,NULL);

  //-For 2D simulations always overrides the 2nd component (Y axis).
  //-Para simulaciones 2D anula siempre la 2nd componente.
  if(Simulate2D)cusph::Resety(Np-Npb,Npb,Aceg);

  //-Computes Tau for Laminar+SPS.
  if(!MultiPhase && lamsps)cusph::ComputeSpsTau(Np,Npb,SpsSmag,SpsBlin,Velrhopg,SpsGradvelg,SpsTaug); //<vs_non-Newtonian>

  //-Applies DDT.
  if(Deltag)cusph::AddDelta(Np-Npb,Deltag+Npb,Arg+Npb);//-Adds the Delta-SPH correction for the density. | Anhade correccion de Delta-SPH a Arg[]. 
  Check_CudaErroor("Failed while executing kernels of interaction.");

  //-Calculates maximum value of ViscDt.
  if(Np)ViscDtMax=cusph::ReduMaxFloat(Np,0,ViscDtg,CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np)));
  if(Np && MultiPhase)ViscEtaDtMax = cusph::ReduMaxFloat(Np, 0, ViscEtaDtg, CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np)));
  //-Calculates maximum value of Ace (periodic particles are ignored). ViscDtg is used like auxiliary memory.
  AceMax=ComputeAceMax(ViscDtg); 

  Check_CudaErroor("Failed in reduction of viscdt.");
  TmgStop(Timers,TMG_CfForces);
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphGpuSingle::MdbcBoundCorrection(){
  TmgStart(Timers,TMG_CfPreForces);
  unsigned n=NpbOk;
  cusph::Interaction_MdbcCorrection(TKernel,Simulate2D,SlipMode,MdbcFastSingle
    ,n,CaseNbound,MdbcThreshold,DivData,Map_PosMin,Posxyg,Poszg,PosCellg,Codeg
    ,Idpg,BoundNormalg,MotionVelg,Velrhopg);
  TmgStop(Timers,TMG_CfPreForces);
}

//==============================================================================
/// Returns the maximum value of  (ace.x^2 + ace.y^2 + ace.z^2) from Acec[].
/// Devuelve valor maximo de (ace.x^2 + ace.y^2 + ace.z^2) a partir de Acec[].
//==============================================================================
double JSphGpuSingle::ComputeAceMax(float *auxmem){
  const bool check=(PeriActive!=0 || InOut!=NULL);
  float acemax=0;
  const unsigned npf=Np-Npb;
  if(!check)cusph::ComputeAceMod(npf,Aceg+Npb,auxmem);//-Without periodic conditions. | Sin condiciones periodicas.
  else cusph::ComputeAceMod(npf,Codeg+Npb,Aceg+Npb,auxmem);//-With periodic conditions ignores the periodic particles. | Con condiciones periodicas ignora las particulas periodicas.
  if(npf)acemax=cusph::ReduMaxFloat(npf,0,auxmem,CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(npf)));
  return(sqrt(double(acemax)));
}

//==============================================================================
/// Particle interaction and update of particle data according to
/// the computed forces using the Verlet time stepping scheme.
///
/// Realiza interaccion y actualizacion de particulas segun las fuerzas 
/// calculadas en la interaccion usando Verlet.
//==============================================================================
double JSphGpuSingle::ComputeStep_Ver(){
  if(BoundCorr)BoundCorrectionData();      //-Apply BoundCorrection.
  Interaction_Forces(INTERSTEP_Verlet);    //-Interaction.
  const double dt=DtVariable(true);        //-Calculate new dt.
  if(CaseNmoving)CalcMotion(dt);           //-Calculate motion for moving bodies.
  DemDtForce=dt;                           //(DEM)
  if(Shifting)RunShifting(dt);             //-Shifting.
  ComputeVerlet(dt);                       //-Update particles using Verlet (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);     //-Control of floating bodies.
  PosInteraction_Forces();                 //-Free memory used for interaction.
  if(Damping)RunDamping(dt,Np,Npb,Posxyg,Poszg,Codeg,Velrhopg); //-Aplies Damping.
  if(RelaxZones)RunRelaxZone(dt);          //-Generate waves using RZ.
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
  DemDtForce=dt*0.5f;                          //(DEM)
  if(BoundCorr)BoundCorrectionData();          //-Apply BoundCorrection.
  Interaction_Forces(INTERSTEP_SymPredictor);  //-Interaction.
  const double ddt_p=DtVariable(false);        //-Calculate dt of predictor step.
  if(Shifting)RunShifting(dt*.5);              //-Shifting.
  ComputeSymplecticPre(dt);                    //-Apply Symplectic-Predictor to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt*.5,true);       //-Control of floating bodies.
  PosInteraction_Forces();                     //-Free memory used for interaction.
  //-Corrector
  //-----------
  DemDtForce=dt;                               //(DEM)
  RunCellDivide(true);
  Interaction_Forces(INTERSTEP_SymCorrector);  //-Interaction.
  const double ddt_c=DtVariable(true);         //-Calculate dt of corrector step.
  if(Shifting)RunShifting(dt);                 //-Shifting.
  ComputeSymplecticCorr(dt);                   //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);         //-Control of floating bodies.
  PosInteraction_Forces();                     //-Free memory used for interaction.
  if(Damping)RunDamping(dt,Np,Npb,Posxyg,Poszg,Codeg,Velrhopg); //-Aplies Damping.
  if(RelaxZones)RunRelaxZone(dt);              //-Generate waves using RZ.

  SymplecticDtPre=min(ddt_p,ddt_c);            //-Calculate dt for next ComputeStep.
  return(dt);
}

//==============================================================================
/// Updates information in FtObjs[] copying data from GPU.
/// Actualiza informacion en FtObjs[] copiando los datos en GPU.
//==============================================================================
void JSphGpuSingle::UpdateFtObjs(){
  if(FtCount && FtObjsOutdated){
    tdouble3 *fcen=FtoAuxDouble6;
    tfloat3  *fang=FtoAuxFloat15;
    tfloat3  *fvellin=fang+FtCount;
    tfloat3  *fvelang=fvellin+FtCount;
    tfloat3  *facelin=fvelang+FtCount;
    tfloat3  *faceang=facelin+FtCount;
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
void JSphGpuSingle::FtApplyImposedVel(float3 *ftoforcesresg)const{
  tfloat3 *ftoforcesresc=NULL;
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
/// Copies the external forces to FtoExtForcesg array for GPU.
/// Copia las fuerzas externas al array FtoExtForcesg para GPU.
//==============================================================================
void JSphGpuSingle::FtCopyExternalForces(){
  tfloat3 *ftoextforces=FtoAuxFloat15;
  memset(ftoextforces,0,sizeof(tfloat3)*FtCount*2);
  if(FtLinearForce!=NULL && FtAngularForce!=NULL)for(unsigned cf=0;cf<FtCount;cf++){
    //-Adds external linear forces.
    if(FtLinearForce[cf]!=NULL){
      const tfloat3 flin=FtLinearForce[cf]->GetValue3f(TimeStep);
      if(flin.x!=FLT_MAX)ftoextforces[cf*2].x=flin.x;
      if(flin.y!=FLT_MAX)ftoextforces[cf*2].y=flin.y;
      if(flin.z!=FLT_MAX)ftoextforces[cf*2].z=flin.z;
    }
    //-Adds external angular forces.
    if(FtAngularForce[cf]!=NULL){
      const tfloat3 fang=FtAngularForce[cf]->GetValue3f(TimeStep);
      if(fang.x!=FLT_MAX)ftoextforces[cf*2+1].x=fang.x;
      if(fang.y!=FLT_MAX)ftoextforces[cf*2+1].y=fang.y;
      if(fang.z!=FLT_MAX)ftoextforces[cf*2+1].z=fang.z;
    }
  }
  cudaMemcpy(FtoExtForcesg,ftoextforces,sizeof(tfloat3)*FtCount*2,cudaMemcpyHostToDevice);
}

//==============================================================================
/// Process floating objects.
/// Procesa floating objects.
//==============================================================================
void JSphGpuSingle::RunFloating(double dt,bool predictor){
  if(TimeStep>=FtPause){//-Operator >= is used because when FtPause=0 in symplectic-predictor, code would not enter here. | Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
    TmgStart(Timers,TMG_SuFloating);
    //-Initialises forces of floatings.
    cudaMemset(FtoForcesg,0,sizeof(StFtoForces)*FtCount);

    //-Adds accelerations from ForcePoints and Moorings.
    if(ForcePoints){
      //-Initialises forces of floatings.
      StFtoForces *ftoforces=(StFtoForces *)FtoAuxFloat15;
      memset(ftoforces,0,sizeof(StFtoForces)*FtCount);
      ForcePoints->GetFtMotionData(ftoforces);
      //-Copies data to GPU memory.
      cudaMemcpy(FtoForcesg,ftoforces,sizeof(StFtoForces)*FtCount,cudaMemcpyHostToDevice);
    }

    //-Adds acceleration from particles and from external forces to FtoForces[].
    //-Copies the external forces to FtoExtForcesg array.
    if(FtLinearForce!=NULL)FtCopyExternalForces();
    //-Calculate forces summation (face,fomegaace) starting from floating particles in ftoforcessum[].
    cusph::FtCalcForcesSum(PeriActive!=0,FtCount,FtoDatpg,FtoCenterg,FtRidpg,Posxyg,Poszg,Aceg,FtoForcesSumg);
    //-Adds acceleration from particles and from external forces to FtoForcesg[].
    cusph::FtCalcForces(FtCount,Gravity,FtoMassg,FtoAnglesg,FtoInertiaini8g,FtoInertiaini1g,FtoForcesSumg,FtoForcesg,FtoExtForcesg);
    
    //-Calculate data to update floatings / Calcula datos para actualizar floatings.
    cusph::FtCalcForcesRes(FtCount,Simulate2D,dt,FtoVelAceg,FtoCenterg,FtoForcesg,FtoForcesResg,FtoCenterResg);
    //-Applies imposed velocity.
    if(FtLinearVel!=NULL)FtApplyImposedVel(FtoForcesResg);
    //-Applies motion constraints.
    if(FtConstraints)cusph::FtApplyConstraints(FtCount,FtoConstraintsg,FtoForcesg,FtoForcesResg);
    
    //-Saves face and fomegace for debug.
    if(SaveFtAce){
      StFtoForces *ftoforces=(StFtoForces *)FtoAuxFloat15;
      cudaMemcpy(ftoforces,FtoForcesg,sizeof(tfloat3)*FtCount*2,cudaMemcpyDeviceToHost);
      SaveFtAceFun(dt,predictor,ftoforces);
    }

    //-Run floating with Chrono library.
    if(ChronoObjects){      
      TmgStop(Timers,TMG_SuFloating);
      TmgStart(Timers,TMG_SuChrono);
      //-Export data / Exporta datos.
      tfloat3* ftoforces=FtoAuxFloat15;
      cudaMemcpy(ftoforces,FtoForcesg,sizeof(tfloat3)*FtCount*2,cudaMemcpyDeviceToHost);
      for(unsigned cf=0;cf<FtCount;cf++)if(FtObjs[cf].usechrono)
        ChronoObjects->SetFtData(FtObjs[cf].mkbound,ftoforces[cf*2],ftoforces[cf*2+1]);
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
      TmgStop(Timers,TMG_SuChrono);
      TmgStart(Timers,TMG_SuFloating);
    }

    //-Apply movement around floating objects / Aplica movimiento sobre floatings.
    cusph::FtUpdate(PeriActive!=0,predictor,FtCount,dt,FtoDatpg,FtoForcesResg,FtoCenterResg
      ,FtRidpg,FtoCenterg,FtoAnglesg,FtoVelAceg,Posxyg,Poszg,Dcellg,Velrhopg,Codeg);

    //-Stores floating data.
    if(!predictor){
      FtObjsOutdated=true;
      //<vs_ftmottionsv_ini>
      if(FtMotSave && FtMotSave->CheckTime(TimeStep+dt)){
        UpdateFtObjs(); //-Updates floating information on CPU memory.
        FtMotSave->SaveFtDataGpu(TimeStep+dt,Nstep+1,FtObjs,Np,Posxyg,Poszg,FtRidpg);
      }
      //<vs_ftmottionsv_end>
    }

    //-Update data of points in FtForces and calculates motion data of affected floatings.
    if(!predictor && ForcePoints){
      UpdateFtObjs(); //-Updates floating information on CPU memory.
      ForcePoints->UpdatePoints(TimeStep,dt,FtObjs);
      if(Moorings)Moorings->ComputeForces(Nstep,TimeStep,dt,ForcePoints);
      ForcePoints->ComputeFtMotion();
    }
    TmgStop(Timers,TMG_SuFloating);
  }
}

//==============================================================================
/// Runs calculations in configured gauges.
/// Ejecuta calculos en las posiciones de medida configuradas.
//==============================================================================
void JSphGpuSingle::RunGaugeSystem(double timestep,bool saveinput){
  //const bool svpart=(TimeStep>=TimePartNext);
  GaugeSystem->CalculeGpu(timestep,DivData
    ,NpbOk,Npb,Np,Posxyg,Poszg,Codeg,Idpg,Velrhopg,saveinput);
}

 //==============================================================================
/// Compute PIPS information of current particles.
/// Calcula datos de PIPS de particulas actuales.
//==============================================================================
void JSphGpuSingle::ComputePips(bool run){
  if(run || DsPips->CheckRun(Nstep)){
    TimerSim.Stop();
    const double timesim=TimerSim.GetElapsedTimeD()/1000.;
    const unsigned sauxmemg=ArraysGpu->GetArraySize();
    unsigned* auxmemg=ArraysGpu->ReserveUint();
    DsPips->ComputeGpu(Nstep,TimeStep,timesim,Np,Npb,NpbOk
      ,DivData,Dcellg,PosCellg,sauxmemg,auxmemg);
    ArraysGpu->Free(auxmemg);
  }
}

//==============================================================================
/// Initialises execution of simulation.
/// Inicia ejecucion de simulacion.
//==============================================================================
void JSphGpuSingle::Run(std::string appname,JSphCfgRun *cfg,JLog2 *log){
  if(!cfg||!log)return;
  AppName=appname; Log=log;

  //-Selection of GPU.
  //-------------------
  SelecDevice(cfg->GpuId);

  //-Configures timers.
  //-------------------
  TmgCreation(Timers,cfg->SvTimers);
  TmgStart(Timers,TMG_Init);

  //-Load parameters and values of input. | Carga de parametros y datos de entrada.
  //--------------------------------------------------------------------------------
  LoadConfig(cfg);
  LoadCaseParticles();
  VisuConfig();
  ConfigDomain();
  ConfigRunMode("Single-Gpu");
  VisuParticleSummary();

  //-Initialisation of execution variables. | Inicializacion de variables de ejecucion.
  //------------------------------------------------------------------------------------
  InitRunGpu();
  RunGaugeSystem(TimeStep,true);
  if(InOut)InOutInit(TimeStepIni);
  FreePartsInit();
  UpdateMaxValues();
  PrintAllocMemory(GetAllocMemoryCpu(),GetAllocMemoryGpu());
  SaveData(); 
  TmgResetValues(Timers);
  TmgStop(Timers,TMG_Init);
  if(Log->WarningCount())Log->PrintWarningList("\n[WARNINGS]","");
  PartNstep=-1; Part++;

  //-Main Loop.
  //------------
  JTimeControl tc("30,60,300,600");//-Shows information at 0.5, 1, 5 y 10 minutes (before first PART).
  bool partoutstop=false;
  TimerSim.Start();
  TimerPart.Start();
  Log->Print(string("\n[Initialising simulation (")+RunCode+")  "+fun::GetDateTime()+"]");
  if(DsPips)ComputePips(true);
  PrintHeadPart();
  while(TimeStep<TimeMax){
    InterStep=(TStep==STEP_Symplectic? INTERSTEP_SymPredictor: INTERSTEP_Verlet);
    if(ViscoTime)Visco=ViscoTime->GetVisco(float(TimeStep));
    double stepdt=ComputeStep();
    RunGaugeSystem(TimeStep+stepdt);
    if(CaseNmoving)RunMotion(stepdt);
    if(InOut)InOutComputeStep(stepdt);
    else RunCellDivide(true);
    TimeStep+=stepdt;
    LastDt=stepdt;
    partoutstop=(Np<NpMinimum || !Np);
    if(TimeStep>=TimePartNext || partoutstop){
      if(partoutstop){
        Log->PrintWarning("Particles OUT limit reached...");
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
    Nstep++;
    const bool laststep=(TimeStep>=TimeMax || (NstepsBreak && Nstep>=NstepsBreak));
    if(DsPips)ComputePips(laststep);
    if(Part<=PartIni+1 && tc.CheckTime())Log->Print(string("  ")+tc.GetInfoFinish((TimeStep-TimeStepIni)/(TimeMax-TimeStepIni)));
    if(NstepsBreak && Nstep>=NstepsBreak)break; //-For debugging.
  }
  TimerSim.Stop(); TimerTot.Stop();

  //-End of Simulation.
  //--------------------
  FinishRun(partoutstop);
}

//==============================================================================
/// Generates files with output data.
/// Genera los ficheros de salida de datos.
//==============================================================================
void JSphGpuSingle::SaveData(){
  const bool save=(SvData!=SDAT_None && SvData!=SDAT_Info);
  const unsigned npsave=Np-NpbPer-NpfPer; //-Subtracts the periodic particles if they exist. | Resta las periodicas si las hubiera.
  //-Retrieves particle data from the GPU. | Recupera datos de particulas en GPU.
  if(save){
    TmgStart(Timers,TMG_SuDownData);
    unsigned npnormal=ParticlesDataDown(Np,0,false,PeriActive!=0);
    if(npnormal!=npsave)Run_Exceptioon("The number of particles is invalid.");
    TmgStop(Timers,TMG_SuDownData);
  }
  //-Retrieve floating object data from the GPU. | Recupera datos de floatings en GPU.
  if(FtCount){
    TmgStart(Timers,TMG_SuDownData);
    UpdateFtObjs();
    TmgStop(Timers,TMG_SuDownData);
  }
  //-Collects additional information. | Reune informacion adicional.
  TmgStart(Timers,TMG_SuSavePart);
  StInfoPartPlus infoplus;
  memset(&infoplus,0,sizeof(StInfoPartPlus));
  if(SvData&SDAT_Info){
    infoplus.nct=CellDivSingle->GetNct();
    infoplus.npbin=NpbOk;
    infoplus.npbout=Npb-NpbOk;
    infoplus.npf=Np-Npb;
    infoplus.npbper=NpbPer;
    infoplus.npfper=NpfPer;
    infoplus.newnp=(InOut? InOut->GetNewNpPart(): 0);
    infoplus.memorycpualloc=this->GetAllocMemoryCpu();
    infoplus.gpudata=true;
    infoplus.memorynctalloc=infoplus.memorynctused=GetMemoryGpuNct();
    infoplus.memorynpalloc=infoplus.memorynpused=GetMemoryGpuNp();
    TimerSim.Stop();
    infoplus.timesim=TimerSim.GetElapsedTimeD()/1000.;
  }
  //-Obtains current domain limits.
  const tdouble3 vdom[2]={CellDivSingle->GetDomainLimits(true),CellDivSingle->GetDomainLimits(false)};
  //-Stores particle data. | Graba datos de particulas.
  JDataArrays arrays;
  AddBasicArrays(arrays,npsave,AuxPos,Idp,AuxVel,AuxRhop);
  JSph::SaveData(npsave,arrays,1,vdom,&infoplus);
  if(UseNormals && SvNormals)SaveVtkNormalsGpu("normals/Normals.vtk",Part,npsave,Npb,Posxyg,Poszg,Idpg,BoundNormalg);
  TmgStop(Timers,TMG_SuSavePart);
}

//==============================================================================
/// Displays and stores final summary of the execution.
/// Muestra y graba resumen final de ejecucion.
//==============================================================================
void JSphGpuSingle::FinishRun(bool stop){
  float tsim=TimerSim.GetElapsedTimeF()/1000.f,ttot=TimerTot.GetElapsedTimeF()/1000.f;
  JSph::ShowResume(stop,tsim,ttot,true,"");
  Log->Print(" ");
  string hinfo,dinfo;
  if(SvTimers){
    ShowTimers();
    GetTimersInfo(hinfo,dinfo);
    Log->Print(" ");
  }
  if(SvRes)SaveRes(tsim,ttot,hinfo,dinfo);
  Log->PrintFilesList();
  Log->PrintWarningList();
}

