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

/// \file JSphGpuSingle.cpp \brief Implements the class \ref JSphGpuSingle.

#include "JSphGpuSingle.h"
#include "JCellDivGpuSingle.h"
#include "JSphGpu_preloop_iker.h" //<vs_advshift>
#include "JSphGpu_mdbc_iker.h"
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
#include "JSphFlexStruc.h"  //<vs_flexstruc>
#include "JDsPartMotionSave.h"
#include "JDsPartFloatSave.h"
#include "JLinearValue.h"
#include "JDataArrays.h"
#include "JDebugSphGpu.h"
#include "JSphShifting.h"
#include "JSphShiftingAdv.h" //<vs_advshift>
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

  //-Copies particle data from input file.
  AuxPos_c->CopyFrom(PartsLoaded->GetPos(),Np);
  Idp_c   ->CopyFrom(PartsLoaded->GetIdp(),Np);
  Velrho_c->CopyFrom(PartsLoaded->GetVelRho(),Np);
  acfloat3 boundnorc("boundnor",Arrays_Cpu,UseNormals);
  if(UseNormals){
    boundnorc.Memset(0,Np);
    if(PartsLoaded->GetBoundNor())boundnorc.CopyFrom(PartsLoaded->GetBoundNor(),CaseNbound);
    else if(AbortNoNormals)Run_ExceptioonFile(
      "No normal data for mDBC in the input file.",PartsLoaded->GetFileLoaded());
    else Log->PrintWarning(fun::PrintStr("No normal data for mDBC in the input file (%s)."
      ,PartsLoaded->GetFileLoaded().c_str()));
  }

  //-Computes radius of floating bodies.
  if(CaseNfloat && PeriActive!=0 && !PartBegin)
    CalcFloatingRadius(Np,AuxPos_c->cptr(),Idp_c->cptr());

  //-Configures Multi-Layer Pistons according to particles.
  if(MLPistons)MLPistons->PreparePiston(Dp,Np,Idp_c->cptr(),AuxPos_c->cptr());

  //-Loads Code of the particles.
  LoadCodeParticles(Np,Idp_c->cptr(),Code_c->ptr());

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
  //-Computes initial cell of the particles and checks if there are unexpected excluded particles.
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
  if(PeriParent_g)PeriParent_g->CuMemset(255,Np);
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
          //-Resizes the allocated memory for the particles if there is not sufficient space and repeats the search process.
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
                ,Velrho_g->ptr(),AG_PTR(SpsTauRho2_g),VelrhoM1_g->ptr());
            }
            if(TStep==STEP_Symplectic){
              if(PosxyPre_g->Active()!=PoszPre_g->Active() || PoszPre_g->Active()!=VelrhoPre_g->Active())
                Run_Exceptioon("Symplectic data is invalid.");
              cusph::PeriodicDuplicateSymplectic(count,Np,DomCells,perinc,listpg.cptr()
                ,Idp_g->ptr(),Code_g->ptr(),Dcell_g->ptr(),Posxy_g->ptr(),Posz_g->ptr()
                ,Velrho_g->ptr(),AG_PTR(SpsTauRho2_g),PosxyPre_g->ptr(),PoszPre_g->ptr(),VelrhoPre_g->ptr());
            }
            if(UseNormals){
              cusph::PeriodicDuplicateNormals(count,Np,listpg.cptr()
                ,BoundNor_g->ptr(),AG_PTR(MotionVel_g),AG_PTR(MotionAce_g)); //<vs_m2dbc>
            }
            if(PeriParent_g){
              cusph::PeriodicSaveParent(count,Np,listpg.cptr(),PeriParent_g->ptr());
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
    CellDivSingle->SortDataArrays(SpsTauRho2_g->cptr(),spstaug.ptr());
    SpsTauRho2_g->SwapPtr(&spstaug);
  }
  if(UseNormals){
    if(SlipMode<SLIP_NoSlip){
      agfloat3 bnorg("-",Arrays_Gpu,true);
      CellDivSingle->SortDataArrays(BoundNor_g->cptr(),bnorg.ptr());
      BoundNor_g->SwapPtr(&bnorg);
    }
    else{//<vs_m2dbc_ini>
      agfloat3 bnorg("-",Arrays_Gpu,true);
      agfloat3 mvelg("-",Arrays_Gpu,true);
      agfloat3 maceg("-",Arrays_Gpu,true);
      CellDivSingle->SortDataArrays(BoundNor_g->cptr(),MotionVel_g->cptr(),MotionAce_g->cptr()
        ,bnorg.ptr(),mvelg.ptr(),maceg.ptr());
      BoundNor_g ->SwapPtr(&bnorg);
      MotionVel_g->SwapPtr(&mvelg);
      MotionAce_g->SwapPtr(&maceg);
    }//<vs_m2dbc_end>
  }
  if(ShiftingAdv!=NULL){//<vs_advshift_ini>
      aguint      fstypeg     ("-",Arrays_Gpu,true);
      agfloat4    shiftvelg   ("-",Arrays_Gpu,true);
      CellDivSingle->SortDataArrays(FSType_g->cptr(),ShiftVel_g->cptr(),fstypeg.ptr(),shiftvelg.ptr());
      FSType_g  ->SwapPtr(&fstypeg);
      ShiftVel_g->SwapPtr(&shiftvelg);
  }//<vs_advshift_end>
  if(PeriParent_g){
    aguint auxg("-",Arrays_Gpu,true);
    aguint periparentg("-",Arrays_Gpu,true);
    CellDivSingle->SortArrayPeriParent(auxg.ptr(),PeriParent_g->cptr(),periparentg.ptr());
    PeriParent_g->SwapPtr(&periparentg);
  }
  if(FlexStruc&&FlexStrucRidpg)CellDivSingle->UpdateIndices(CaseNflexstruc,FlexStrucRidpg); //<vs_flexstruc>

  #ifdef AVAILABLE_DIVCLEAN
  if(PsiClean_g){//<vs_divclean>
    agfloat   psiclean("-",Arrays_Gpu,true);
    CellDivSingle->SortDataArrays(PsiClean_g->cptr(),psiclean.ptr());
    PsiClean_g->SwapPtr(&psiclean);
  }
  #endif

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
  //JDebugSphGpu::SaveVtk("_DG_Divide_End.vtk",Nstep,0,Np,"idp,vel,rho",this);
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

//<vs_advshift_ini>
//==============================================================================
/// PreLoop for additional models computation.
//==============================================================================
void JSphGpuSingle::PreLoopProcedure(TpInterStep interstep){
  const bool runshift=(ShiftingAdv && interstep==INTERSTEP_SymPredictor && Nstep!=0);
  if(runshift){
    Timersg->TmStart(TMG_SuShifting,false);
    ComputeFSParticles();
    ComputeUmbrellaRegion();
    const unsigned bsfluid=BlockSizes.forcesfluid;
    cusph::PreLoopInteraction(TKernel,Simulate2D,runshift,bsfluid,Np-Npb
      ,Npb,DivData,Dcell_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr()
      ,Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr(),FSType_g->ptr()
      ,FSNormal_g->ptr(),FSMinDist_g->ptr(),NULL);
    cusph::ComputeShiftingVel(bsfluid,Np-Npb,Npb,Simulate2D,ShiftingAdv->GetShiftCoef()
      ,ShiftingAdv->GetAleActive(),float(SymplecticDtPre),FSType_g->cptr()
      ,FSNormal_g->cptr(),FSMinDist_g->cptr(),ShiftVel_g->ptr(),NULL);
    //-Updates pre-loop variables in periodic particles.
    if(PeriParent_g){
      cusph::PeriPreLoopCorr(Np,0,PeriParent_g->cptr(),FSType_g->ptr()
        ,ShiftVel_g->ptr());
    }
    //-Saves VTK for debug.
    if(0 && TimeStep+LastDt>=TimePartNext){
		  DgSaveVtkParticlesGpu("Compute_FreeSurface_",Part,0,Np,Posxy_g->cptr()
        ,Posz_g->cptr(),Code_g->cptr(),FSType_g->cptr(),ShiftVel_g->cptr()
        ,FSNormal_g->cptr());
	  }
    Timersg->TmStop(TMG_SuShifting,true);
  }
}

//==============================================================================
/// Compute free-surface particles and their normals.
//==============================================================================
void JSphGpuSingle::ComputeFSParticles(){
  const unsigned bsfluid=BlockSizes.forcesfluid;
  aguint fspartg("-",Arrays_Gpu,true);
  cusph::ComputeFSNormals(TKernel,Simulate2D,bsfluid,Npb,Np-Npb,DivData
    ,Dcell_g->cptr(),Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr(),Velrho_g->cptr()
    ,Code_g->cptr(),FtoMasspg,ShiftVel_g->ptr(),FSType_g->ptr(),FSNormal_g->ptr()
    ,fspartg.ptr(),NULL);
}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
//==============================================================================
void JSphGpuSingle::ComputeUmbrellaRegion(){
  const unsigned bsfluid=BlockSizes.forcesfluid;
  aguint fspartg("-",Arrays_Gpu,true);
  cusph::ComputeUmbrellaRegion(TKernel,Simulate2D,bsfluid,Npb,Np-Npb,DivData
    ,Dcell_g->cptr(),PosCell_g->cptr(),Code_g->cptr(),FSNormal_g->cptr()
    ,fspartg.ptr(),FSType_g->ptr(),NULL);
}
//<vs_advshift_end>

//==============================================================================
/// Interaction for force computation.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphGpuSingle::Interaction_Forces(TpInterStep interstep){
  float3* dengradcorr=NULL;

  Timersg->TmStart(TMG_CfForces,true);
  const unsigned bsfluid=BlockSizes.forcesfluid;
  const unsigned bsbound=BlockSizes.forcesbound;

  //<vs_advshift_ini>
  const bool shiftadv=(ShiftingAdv!=NULL);
  const bool corrector=(InterStep==INTERSTEP_SymCorrector);
  const bool aleform=(shiftadv? ShiftingAdv->GetAleActive(): false);
  const bool ncpress=(shiftadv? ShiftingAdv->GetNcPress()  : false);
  //<vs_advshift_end>

  //-Interaction of Fluid-Fluid/Bound & Bound-Fluid (forces).
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
    ,AG_CPTR(PsiClean_g),AG_PTR(PsiCleanRhs_g),AG_PTR(CsPsiClean_g),DivCleanKp,DivCleaning
    ,NULL,NULL);
  cusph::Interaction_Forces(parms);

  //-DEM interaction of Floating-Bound & Floating-Floating.
  if(UseDEM)cusph::Interaction_ForcesDem(BlockSizes.forcesdem,CaseNfloat
    ,DivData,Dcell_g->cptr(),RidpMotg+CaseNmoving,DemDatag,FtoMasspg
    ,float(DemDtForce),PosCell_g->cptr(),Velrho_g->cptr(),Code_g->cptr()
    ,Idp_g->cptr(),ViscDt_g->ptr(),Ace_g->ptr(),NULL);

  //<vs_flexstruc_ini>
  //-Interaction flexible structure-flexible structure.
  if(FlexStruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    const StInterParmsFlexStrucg parmsfs=StrInterParmsFlexStrucg(Simulate2D,TKernel,TVisco,TMdbc2
      ,Visco*ViscoBoundFactor,CaseNflexstruc,DivData,Dcell_g->cptr()
      ,PosCell_g->cptr(),Velrho_g->cptr(),Code_g->cptr()
      ,AG_CPTR(BoundMode_g),AG_CPTR(TangenVel_g)
      ,FlexStrucDatag,FlexStrucRidpg,PosCell0g,NumPairsg,PairIdxg,KerCorrg,BoundNor0g,DefGradg,AG_PTR(BoundNor_g),FlexStrucDtg,Ace_g->ptr(),NULL);
    cusph::Interaction_ForcesFlexStruc(parmsfs);
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end>

  //-For 2-D simulations zero the 2nd component.
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

  #ifdef AVAILABLE_DIVCLEAN
  if(DivCleaning){
    if(Np)CsPsiCleanMax=cusph::ReduMaxFloat(Np,0,CsPsiClean_g->ptr(),CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np)));
  }
  #endif
  //<vs_flexstruc_ini>
  //-Calculates maximum value of FlexStrucDt.
  if(CaseNflexstruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    FlexStrucDtMax=cusph::ReduMaxFloat(CaseNflexstruc,0,FlexStrucDtg,CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(CaseNflexstruc)));
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end> 

  InterNum++;
  Timersg->TmStop(TMG_CfForces,true);
  Check_CudaErroor("Failed in reduction of viscdt.");
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphGpuSingle::MdbcBoundCorrection(TpInterStep interstep){
  const bool runmdbc=(TBoundary==BC_MDBC
    && (MdbcCorrector || interstep!=INTERSTEP_SymCorrector));
  if(runmdbc){
    Timersg->TmStart(TMG_CfPreMDBC,false);
    if(SlipMode==SLIP_Vel0){
      const unsigned n=(UseNormalsFt? Np: NpbOk);
      cusph::Interaction_MdbcCorrection(TKernel,Simulate2D,n,CaseNbound
        ,DivData,Map_PosMin,Posxy_g->cptr(),Posz_g->cptr(),PosCell_g->cptr()
        ,Code_g->cptr(),Idp_g->cptr(),BoundNor_g->cptr(),Velrho_g->ptr());
    }
    else{
  //else if(SlipMode==SLIP_NoSlip){ //<vs_m2dbc_ini>
      const unsigned n=(UseNormalsFt? Np: Npb);
      BoundMode_g->Reserve();     //-BoundOnOff_g is freed in PosInteraction_Forces().
      BoundMode_g->CuMemset(0,n); //-BoundMode_g[]=0=BMODE_DBC
      TangenVel_g->Reserve();     //-TangenVel_g is freed in PosInteraction_Forces().
      cusph::Interaction_Mdbc2Correction(TKernel,Simulate2D,SlipMode,n,CaseNbound
        ,Gravity,DivData,Map_PosMin,Posxy_g->cptr(),Posz_g->cptr()
        ,PosCell_g->cptr(),Code_g->cptr(),Idp_g->cptr(),BoundNor_g->cptr()
        ,MotionVel_g->cptr(),MotionAce_g->cptr(),Velrho_g->ptr()
        ,BoundMode_g->ptr(),TangenVel_g->ptr());
    } //<vs_m2dbc_end>
    //else Run_Exceptioon("Error: SlipMode is invalid.");
    Timersg->TmStop(TMG_CfPreMDBC,true);
  }
}


//==============================================================================
/// Returns the maximum value of  (ace.x^2 + ace.y^2 + ace.z^2) from Acec[].
/// Devuelve valor maximo de (ace.x^2 + ace.y^2 + ace.z^2) a partir de Acec[].
//==============================================================================
double JSphGpuSingle::ComputeAceMax(float* auxmemg){
  const bool check=(PeriActive!=0 || InOut!=NULL);
  const unsigned pini=(CaseNflexstruc? 0: Npb);
  float acemax=0;
  if(check)cusph::ComputeAceMod(Np-pini,Code_g->cptr()+pini,Ace_g->cptr()+pini,auxmemg);
  else     cusph::ComputeAceMod(Np-pini,Ace_g->cptr()+pini,auxmemg);      
  if(Np-pini)acemax=cusph::ReduMaxFloat(Np-pini,0,auxmemg,CellDivSingle->GetAuxMem(cusph::ReduMaxFloatSize(Np-pini)));
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
  InterStep=INTERSTEP_Verlet;
  MdbcBoundCorrection(InterStep);        //-mDBC correction
  PreInteraction_Forces(InterStep);      //-Allocating temporary arrays.
  Interaction_Forces(InterStep);         //-Interaction.
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
  if(CaseNmoving)CalcMotion(dt);          //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
  InterStep=INTERSTEP_SymPredictor;
  DemDtForce=dt*0.5f;                     //-For DEM interaction.
  MdbcBoundCorrection(InterStep);         //-mDBC correction
  PreInteraction_Forces(InterStep);       //-Allocating temporary arrays.
  PreLoopProcedure(InterStep);            //-Pre-calculation for advanced shifting and other formulations. //<vs_advshift>
  Interaction_Forces(InterStep);          //-Interaction.
  const double dt_p=DtVariable(false);    //-Calculate dt of predictor step.
  if(Shifting)RunShifting(dt*.5);         //-Standard shifting.
  ComputeSymplecticPre(dt);               //-Apply Symplectic-Predictor to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt*.5,true);  //-Control of floating bodies.
  PosInteraction_Forces();                //-Free memory used for interaction.
  //-Corrector
  //-----------
  InterStep=INTERSTEP_SymCorrector;
  DemDtForce=dt;                          //-For DEM interaction.
  RunCellDivide(true);                    //-Rearrange particles in cells.
  if(FlexStruc)UpdateFlexStrucGeometry(); //-Update the geometric information for each flexible structure particle.
  MdbcBoundCorrection(InterStep);         //-mDBC correction.
  PreInteraction_Forces(InterStep);       //-Allocating temporary arrays.
  PreLoopProcedure(InterStep);            //-Pre-calculation for advanced shifting and other formulations. //<vs_advshift>
  Interaction_Forces(InterStep);          //-Interaction.
  const double dt_c=DtVariable(true);     //-Calculate dt of corrector step.
  if(Shifting)RunShifting(dt);            //-Shifting.
  ComputeSymplecticCorr(dt);              //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);    //-Control of floating bodies.
  PosInteraction_Forces();                //-Free memory used for interaction.
  if(Damping)RunDamping(dt);              //-Aplies Damping.
  if(RelaxZones)RunRelaxZone(dt);         //-Generate waves using RZ.

  SymplecticDtPre=min(dt_p,dt_c);         //-Calculate dt for next ComputeStep.
  return(dt);
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
    //-Computes sum of linear and angular acceleration of floating particles.
    cusph::FtPartsSumAce(PeriActive!=0,FtCount,FtoDatpg,FtoCenterg,RidpMotg
      ,Posxy_g->cptr(),Posz_g->cptr(),Ace_g->cptr(),FtoAceg);
    cudaMemcpy(Fto_AceLinAng,FtoAceg,sizeof(tfloat6)*FtCount,cudaMemcpyDeviceToHost);
    //-Compute new linear and angular acceleration, velocity and center to update floatings.
    FtComputeAceVel(dt,predictor,saveftvalues,Fto_AceLinAng,Fto_VelLinAng,Fto_Center);
  }

  //-Run floatings with Chrono library (change computed center, linear and angular velocity).
  if(!ftpaused && ChronoObjects){
    Timersg->TmStop(TMG_SuFloating,true);
    Timersg->TmStart(TMG_SuChrono,false);
    FtComputeChrono(dt,predictor,Fto_AceLinAng,Fto_VelLinAng,Fto_Center);
    Timersg->TmStop(TMG_SuChrono,false);
    Timersg->TmStart(TMG_SuFloating,false);
  }

  if(!ftpaused){//-Operator >= is used because when FtPause=0 in symplectic-predictor, code would not enter here. | Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
    //-Apply displacement and new velocity to floating particles.
    const bool updatenormals=(!predictor && UseNormalsFt);
    tmatrix4d matc=TMatrix4d();
    for(unsigned cf=0;cf<FtCount;cf++){
      //-Get Floating object values.
      const StFloatingData fobj=FtObjs[cf];
      const float fradius=fobj.radius;
      const unsigned fpini=fobj.begin-CaseNfixed;
      const unsigned fnp=fobj.count;
      const tfloat3  fvel  =Fto_VelLinAng[cf].getlo();
      const tfloat3  fomega=Fto_VelLinAng[cf].gethi();
      const tdouble3 fcenter=Fto_Center[cf];
      if(updatenormals){
        //-Compute matrix to update floating normals for mDBC.
        JMatrix4d mat;
        const tdouble3 dang=(ToTDouble3(fomega)*dt)*TODEG;
        const tdouble3 cen2=(PeriActive? UpdatePeriodicPos(fcenter): fcenter);
        mat.Move(cen2);
        mat.Rotate(dang);
        mat.Move(fobj.center*-1);
        matc=mat.GetMatrix();
      }
      //-Run CUDA kernel.
      cudaStream_t stm=(NStmFloatings? StmFloatings[cf%NStmFloatings]: NULL);
      cusph::FtPartsUpdate(PeriActive!=0,dt,updatenormals
        ,fnp,fpini,fradius,matc,fvel,fomega,fcenter
        ,RidpMotg,Posxy_g->ptr(),Posz_g->ptr(),Velrho_g->ptr()
        ,Dcell_g->ptr(),Code_g->ptr(),AG_PTR(BoundNor_g)
        ,AG_PTR(MotionVel_g),AG_PTR(MotionAce_g),stm);
    }
    if(NStmFloatings)cudaDeviceSynchronize();

    //-Update floating data (FtObjs[]) for next step.
    if(!predictor){
      FtUpdateFloatings(dt,Fto_VelLinAng,Fto_Center);
      //-Update center in GPU memory (FtoCenterg[]).
      for(unsigned cf=0;cf<FtCount;cf++)FtoCenterc[cf]=FtObjs[cf].center;
      cudaMemcpy(FtoCenterg,FtoCenterc,sizeof(tdouble3)*FtCount,cudaMemcpyHostToDevice);
    }
  }

  //-Saves current floating data for output.
  if(saveftvalues){
    if(PartFloatSave)PartFloatSave->SetFtData(Part,TimeStep+dt,Nstep+1,FtObjs,ForcePoints);
  }
  Timersg->TmStop(TMG_SuFloating,true);

  //-Update data of points in FtForces and calculates motion data of affected floatings.
  if(!predictor && ForcePoints){
    Timersg->TmStart(TMG_SuMoorings,false);
    ForcePoints->UpdatePoints(TimeStep,dt,ftpaused,FtObjs);
    if(Moorings)Moorings->ComputeForces(Nstep,TimeStep,dt,ForcePoints);
    ForcePoints->ComputeForcesSum();
    Timersg->TmStop(TMG_SuMoorings,true);
  }
}

//==============================================================================
/// Runs first calculations in configured gauges.
/// Ejecuta primeros calculos en las posiciones de medida configuradas.
//==============================================================================
void JSphGpuSingle::RunFirstGaugeSystem(double timestep){
  GaugeSystem->ConfigArraysGpu(0,Posxy_g,Posz_g,Code_g,Idp_g,Velrho_g);  
  GaugeSystem->CalculeGpu(0,timestep,DivData,NpbOk,Npb,Np,true);
}

//==============================================================================
/// Runs calculations in configured gauges.
/// Ejecuta calculos en las posiciones de medida configuradas.
//==============================================================================
void JSphGpuSingle::RunGaugeSystem(double timestep){
  if(GaugeSystem->GetCount()){
    Timersg->TmStart(TMG_SuGauges,false);
    //const bool svpart=(TimeStep>=TimePartNext);
    GaugeSystem->CalculeGpu(0,timestep,DivData,NpbOk,Npb,Np,false);
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
    const double timesim=TimerSim.GetSecs();
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
  RunFirstGaugeSystem(TimeStep);
  if(InOut)InOutInit(TimeStepIni);
  if(FlexStruc)FlexStrucInit(); //<vs_flexstruc>
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
    if(FlexStruc)UpdateFlexStrucGeometry();
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
    infoplus.timesim=TimerSim.GetSecs();
  }
  else infoplus.SetBasic(Np,npnormal,Np-Npb,CellDivSingle->GetNct());

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
  const float tsim=(float)TimerSim.GetSecs();
  const float ttot=(float)TimerTot.GetSecs();
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

//<vs_flexstruc_ini>
//==============================================================================
/// Initialises the associated arrays for the flexible structures.
/// Inicializa las matrices asociadas para las estructuras flexibles.
//==============================================================================
void JSphGpuSingle::FlexStrucInit(){
  //-Start timer and print info.
  Timersg->TmStart(TMG_SuFlexStruc,false);
  Log->Print("\nInitialising Flexible Structures...");
  //-Allocate array.
  size_t m=0;
  m=sizeof(StFlexStrucData)*FlexStrucCount; cudaMalloc((void**)&FlexStrucDatag,m); MemGpuFixed+=m;
  //-Get flexible structure data for each body and copy to GPU
  vector<StFlexStrucData> flexstrucdata(FlexStrucCount);
  for(unsigned c=0;c<FlexStrucCount;c++){
    std::vector<typecode> clampcode=FlexStruc->GetBody(c)->GetClampCode();
    flexstrucdata[c].nc=unsigned(clampcode.size());
    std::copy(clampcode.begin(),clampcode.end(),flexstrucdata[c].clampcode);
    flexstrucdata[c].vol0=FlexStruc->GetBody(c)->GetParticleVolume();
    flexstrucdata[c].rho0=FlexStruc->GetBody(c)->GetDensity();
    flexstrucdata[c].youngmod=FlexStruc->GetBody(c)->GetYoungMod();
    flexstrucdata[c].poissratio=FlexStruc->GetBody(c)->GetPoissRatio();
    flexstrucdata[c].hgfactor=FlexStruc->GetBody(c)->GetHgFactor();
    flexstrucdata[c].cmat=FlexStruc->GetBody(c)->GetConstitMatrix();
  }
  cudaMemcpy(FlexStrucDatag,flexstrucdata.data(),sizeof(StFlexStrucData)*FlexStrucCount,cudaMemcpyHostToDevice);
  //-Configure code for flexible structures.
  cudaMemcpy(Code_c->ptr(),Code_g->ptr(),sizeof(typecode)*Npb,cudaMemcpyDeviceToHost);
  FlexStruc->ConfigCode(Npb,Code_c->ptr());
  cudaMemcpy(Code_g->ptr(),Code_c->ptr(),sizeof(typecode)*Npb,cudaMemcpyHostToDevice);
  cusph::SetFlexStrucClampCodes(Npb,PosCell_g->cptr(),FlexStrucDatag,Code_g->ptr());
  cudaMemcpy(Code_c->ptr(),Code_g->cptr(),sizeof(typecode)*Npb,cudaMemcpyDeviceToHost);
  //-Count number of flexible structure particles.
  CaseNflexstruc=cusph::CountFlexStrucParts(Npb,Code_g->cptr());
  //-Allocate arrays.
  m=sizeof(unsigned)*CaseNflexstruc;              cudaMalloc((void**)&FlexStrucRidpg,   m);  MemGpuFixed+=m;
  m=sizeof(float4)*CaseNflexstruc;                cudaMalloc((void**)&PosCell0g,        m);  MemGpuFixed+=m;
  m=sizeof(unsigned)*CaseNflexstruc;              cudaMalloc((void**)&NumPairsg,        m);  MemGpuFixed+=m;
  m=sizeof(unsigned*)*CaseNflexstruc;             cudaMalloc((void**)&PairIdxg,         m);  MemGpuFixed+=m;
  m=sizeof(tmatrix3f)*CaseNflexstruc;             cudaMalloc((void**)&KerCorrg,         m);  MemGpuFixed+=m;
  m=sizeof(tmatrix3f)*CaseNflexstruc;             cudaMalloc((void**)&DefGradg,         m);  MemGpuFixed+=m;
  if(UseNormals)m=sizeof(float3)*CaseNflexstruc;  cudaMalloc((void**)&BoundNor0g,       m);  MemGpuFixed+=m;
  m=sizeof(float)*CaseNflexstruc;                 cudaMalloc((void**)&FlexStrucDtg,     m);  MemGpuFixed+=m;
  //-Calculate array for indexing into flexible structure particles.
  cusph::CalcFlexStrucRidp(Npb,Code_g->cptr(),FlexStrucRidpg);
  //-Copy current position and normals into initial position and normals.
  cusph::GatherToFlexStrucArray(CaseNflexstruc,FlexStrucRidpg,PosCell_g->cptr(),PosCell0g);
  if(UseNormals)cusph::GatherToFlexStrucArray(CaseNflexstruc,FlexStrucRidpg,BoundNor_g->cptr(),BoundNor0g);
  //-Get number of particle pairs for each flexible structure particle.
  const unsigned numpairstot=cusph::CountFlexStrucPairs(CaseNflexstruc,PosCell0g,NumPairsg);
  //-Download number of pairs for each particle.
  vector<unsigned> numpairs(CaseNflexstruc);
  cudaMemcpy(numpairs.data(),NumPairsg,sizeof(unsigned)*CaseNflexstruc,cudaMemcpyDeviceToHost);
  //-Allocate memory for raw buffer for storing pair indices and set the pointers to the indices.
  m=sizeof(unsigned)*numpairstot; cudaMalloc((void**)&PairIdxBufferg,m);  MemGpuFixed+=m;
  unsigned* offset=PairIdxBufferg;
  vector<unsigned*> pairidx(CaseNflexstruc);
  for(unsigned p=0;p<CaseNflexstruc;p++){
    pairidx[p]=offset;
    offset+=numpairs[p];
  }
  cudaMemcpy(PairIdxg,pairidx.data(),sizeof(unsigned*)*CaseNflexstruc,cudaMemcpyHostToDevice);
  //-Set the indices for each particle pair.
  cusph::SetFlexStrucPairs(CaseNflexstruc,PosCell0g,PairIdxg);
  //-Interaction parameters.
  const StInterParmsFlexStrucg parmsfs=StrInterParmsFlexStrucg(Simulate2D,TKernel,TVisco,TMdbc2
    ,Visco*ViscoBoundFactor,CaseNflexstruc,DivData,Dcell_g->cptr()
    ,PosCell_g->cptr(),Velrho_g->cptr(),Code_g->cptr()
    ,AG_CPTR(BoundMode_g),AG_CPTR(TangenVel_g)
    ,FlexStrucDatag,FlexStrucRidpg,PosCell0g,NumPairsg,PairIdxg,KerCorrg,BoundNor0g,DefGradg,AG_PTR(BoundNor_g),FlexStrucDtg,Ace_g->ptr(),NULL);
  //-Calculate kernel correction and update geometry.
  cusph::CalcFlexStrucKerCorr(parmsfs);
  cusph::UpdateFlexStrucGeometry(parmsfs);
  Log->Print("");
  Timersg->TmStop(TMG_SuFlexStruc,false);
}

//==============================================================================
/// Updates the geometric information for each flexible structure particle.
/// Actualiza la información geométrica de cada partícula de estructura flexible.
//==============================================================================
void JSphGpuSingle::UpdateFlexStrucGeometry(){
  Timersg->TmStart(TMG_SuFlexStruc,false);
  const StInterParmsFlexStrucg parmsfs=StrInterParmsFlexStrucg(Simulate2D,TKernel,TVisco,TMdbc2
    ,Visco*ViscoBoundFactor,CaseNflexstruc,DivData,Dcell_g->cptr()
    ,PosCell_g->cptr(),Velrho_g->cptr(),Code_g->cptr()
    ,AG_CPTR(BoundMode_g),AG_CPTR(TangenVel_g)
    ,FlexStrucDatag,FlexStrucRidpg,PosCell0g,NumPairsg,PairIdxg,KerCorrg,BoundNor0g,DefGradg,AG_PTR(BoundNor_g),FlexStrucDtg,Ace_g->ptr(),NULL);
  cusph::UpdateFlexStrucGeometry(parmsfs);
  Timersg->TmStop(TMG_SuFlexStruc,false);
}
//<vs_flexstruc_end>

