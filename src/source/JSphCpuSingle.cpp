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

/// \file JSphCpuSingle.cpp \brief Implements the class \ref JSphCpuSingle.

#include "JSphCpuSingle.h"
#include "JCellDivCpuSingle.h"
#include "JSphMk.h"
#include "JPartsLoad4.h"
#include "Functions.h"
#include "FunctionsMath.h"
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
#include "JDsGaugeSystem.h"
#include "JSphInOut.h"
#include "JSphFlexStruc.h"  //<vs_flexstruc>
#include "JDsPartMotionSave.h"
#include "JDsPartFloatSave.h"
#include "JLinearValue.h"
#include "JDataArrays.h"
#include "JDebugSphCpu.h"
#include "JSphShifting.h"
#include "JSphShiftingAdv.h" //<vs_advshift>
#include "JDsPips.h"
#include "JDsExtraData.h"
#include "JDsOutputParts.h" //<vs_outpaarts>

#include <climits>

using namespace std;
//==============================================================================
/// Constructor.
//==============================================================================
JSphCpuSingle::JSphCpuSingle():JSphCpu(false){
  ClassName="JSphCpuSingle";
  CellDivSingle=NULL;
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphCpuSingle::~JSphCpuSingle(){
  DestructorActive=true;
  delete CellDivSingle; CellDivSingle=NULL;
}

//==============================================================================
/// Return memory reserved in CPU.
/// Devuelve la memoria reservada en cpu.
//==============================================================================
llong JSphCpuSingle::GetAllocMemoryCpu()const{  
  llong s=JSphCpu::GetAllocMemoryCpu();
  //-Allocated in other objects.
  if(CellDivSingle)s+=CellDivSingle->GetAllocMemory();
  return(s);
}

//==============================================================================
/// Update maximum values of memory, particles & cells.
/// Actualiza los valores maximos de memory, particles y cells.
//==============================================================================
void JSphCpuSingle::UpdateMaxValues(){
  const llong mcpu=GetAllocMemoryCpu();
  MaxNumbers.memcpu=max(MaxNumbers.memcpu,mcpu);
  MaxNumbers.memgpu=MaxNumbers.memgpunct=0;
  MaxNumbers.particles=max(MaxNumbers.particles,Np);
  if(CellDivSingle)MaxNumbers.cells=max(MaxNumbers.cells,CellDivSingle->GetNct());
  //-Prints allocated arrays when it is changed.
  if(Arrays_Cpu->GetArrayCountUpdated()){
    PrintAllocMemory(GetAllocMemoryCpu());
  }
}

//==============================================================================
/// Load the execution configuration.
/// Carga la configuracion de ejecucion.
//==============================================================================
void JSphCpuSingle::LoadConfig(const JSphCfgRun* cfg){
  //-Load OpenMP configuraction. | Carga configuracion de OpenMP.
  ConfigOmp(cfg);
  //-Load basic general configuraction. | Carga configuracion basica general.
  JSph::LoadConfig(cfg);
  //-Checks compatibility of selected options.
  Log->Print("**Special case configuration is loaded");
}

//==============================================================================
/// Configuration of current domain.
/// Configuracion del dominio actual.
//==============================================================================
void JSphCpuSingle::ConfigDomain(){
  //-Configure cell map division (defines ScellDiv, Scell, Map_Cells). 
  ConfigCellDivision();
  //-Calculate number of particles. | Calcula numero de particulas.
  Np=PartsLoaded->GetCount();
  Npb=CaseNpb;
  NpbOk=Npb;

  //-Allocates fixed CPU memory for moving & floating particles.
  AllocCpuMemoryFixed();
  //-Allocates CPU memory for particles.
  AllocCpuMemoryParticles(Np);

  //-Copies particle data from input file.
  Pos_c   ->CopyFrom(PartsLoaded->GetPos(),Np);
  Idp_c   ->CopyFrom(PartsLoaded->GetIdp(),Np);
  Velrho_c->CopyFrom(PartsLoaded->GetVelRho(),Np);
  if(UseNormals){
    BoundNor_c->Memset(0,Np);
    if(PartsLoaded->GetBoundNor())BoundNor_c->CopyFrom(PartsLoaded->GetBoundNor(),CaseNbound);
    else if(AbortNoNormals)Run_ExceptioonFile(
      "No normal data for mDBC in the input file.",PartsLoaded->GetFileLoaded());
    else Log->PrintWarning(fun::PrintStr("No normal data for mDBC in the input file (%s)."
      ,PartsLoaded->GetFileLoaded().c_str()));
  }

  //-Computes radius of floating bodies.
  if(CaseNfloat && PeriActive!=0 && !PartBegin)
    CalcFloatingRadius(Np,Pos_c->cptr(),Idp_c->cptr());

  //-Configures Multi-Layer Pistons according to particles.
  if(MLPistons)MLPistons->PreparePiston(Dp,Np,Idp_c->cptr(),Pos_c->cptr());

  //-Load particle code. | Carga code de particulas.
  LoadCodeParticles(Np,Idp_c->cptr(),Code_c->ptr());

  //-Creates PartsInit object with initial particle data for automatic configurations.
  CreatePartsInit(Np,Pos_c->cptr(),Code_c->cptr());

  //-Runs initialization operations from XML.
  RunInitialize(Np,Npb,Pos_c->cptr(),Idp_c->cptr(),Code_c->cptr()
    ,Velrho_c->ptr(),AC_PTR(BoundNor_c));
  if(UseNormals)ConfigBoundNormals(Np,Npb,Pos_c->cptr(),Idp_c->cptr(),BoundNor_c->ptr());

  //-Computes MK domain for boundary and fluid particles.
  MkInfo->ComputeMkDomains(Np,Pos_c->cptr(),Code_c->cptr());

  //-Sets local domain of the simulation within Map_Cells and computes DomCellCode.
  //-Establece dominio de simulacion local dentro de Map_Cells y calcula DomCellCode.
  SelecDomain(TUint3(0,0,0),Map_Cells);
  //-Computes initial cell of the particles and checks if there are unexpected excluded particles.
  //-Calcula celda inicial de particulas y comprueba si hay excluidas inesperadas.
  LoadDcellParticles(Np,Code_c->cptr(),Pos_c->cptr(),Dcell_c->ptr());

  //-Creates object for Celldiv on the CPU and selects a valid cellmode.
  //-Crea objeto para divide en CPU y selecciona un cellmode valido.
  CellDivSingle=new JCellDivCpuSingle(Stable,FtCount!=0,PeriActive,CellDomFixed,CellMode
    ,Scell,Map_PosMin,Map_PosMax,Map_Cells,CaseNbound,CaseNfixed,CaseNpb,DirOut);
  CellDivSingle->DefineDomain(DomCellCode,DomCelIni,DomCelFin,DomPosMin,DomPosMax);
  ConfigCellDiv((JCellDivCpu*)CellDivSingle);

  ConfigSaveData(0,1,"",Np,Pos_c->cptr(),Idp_c->cptr());

  //-Reorders particles according to cells.
  BoundChanged=true;
  RunCellDivide(true);
}

//==============================================================================
/// Resizes the allocated memory for particles on the CPU saving the current 
/// data (ndatacpu) and measures the time spent using TMC_SuResizeNp. 
/// At the end updates the division.
///
/// Redimensiona el espacio reservado para particulas en CPU manteniendo los 
/// datos actuales (ndatacpu) y mide el tiempo consumido con TMC_SuResizeNp.
/// Al terminar actualiza el divide.
//==============================================================================
void JSphCpuSingle::ResizeParticlesSizeData(unsigned ndatacpu
  ,unsigned newsize,unsigned minsize,float oversize,bool updatedivide)
{
  Timersc->TmStart(TMC_SuResizeNp);
  newsize+=(oversize>0? unsigned(oversize*newsize): 0);
  //-Resize CPU memory for particles saving current data.
  ResizeCpuMemoryParticlesData(ndatacpu,newsize,minsize);
  Timersc->TmStop(TMC_SuResizeNp);
  if(updatedivide)RunCellDivide(true);
}

//==============================================================================
/// Create list of new periodic particles to duplicate.
/// With stable activated reordered list of periodic particles.
///
/// Crea lista de nuevas particulas periodicas a duplicar.
/// Con stable activado reordena lista de periodicas.
//==============================================================================
unsigned JSphCpuSingle::PeriodicMakeList(unsigned n,unsigned pini,bool stable
  ,unsigned nmax,tdouble3 perinc,const tdouble3* pos,const typecode* code
  ,unsigned* listp)const
{
  unsigned count=0;
  if(n){
    //-Initialize size of list lsph to zero. | Inicializa tamanho de lista lspg a cero.
    listp[nmax]=0;
    for(unsigned p=0;p<n;p++){
      const unsigned p2=p+pini;
      //-Keep normal or periodic particles. | Se queda con particulas normales o periodicas.
      if(CODE_GetSpecialValue(code[p2])<=CODE_PERIODIC){
        //-Get particle position. | Obtiene posicion de particula.
        const tdouble3 ps=pos[p2];
        tdouble3 ps2=ps+perinc;
        if(Map_PosMin<=ps2 && ps2<Map_PosMax){
          unsigned cp=listp[nmax]; listp[nmax]++; if(cp<nmax)listp[cp]=p2;
        }
        ps2=ps-perinc;
        if(Map_PosMin<=ps2 && ps2<Map_PosMax){
          unsigned cp=listp[nmax]; listp[nmax]++; if(cp<nmax)listp[cp]=(p2|0x80000000);
        }
      }
    }
    count=listp[nmax];
    //-Reorder list if it is valid and stability is activated. | Reordena lista si es valida y stable esta activado.
    if(stable && count && count<=nmax){
      //-Don't make mistake because at the moment the list is not created using OpenMP. | No hace falta porque de momento no se crea la lista usando OpenMP.
    }
  }
  return(count);
}

//==============================================================================
/// Duplicate the indicated particle position applying displacement.
/// Duplicated particles are considered to be always valid and are inside
/// of the domain.
/// This kernel works for single-cpu & multi-cpu because the computations are done  
/// starting from domposmin.
/// It is controlled that the coordinates of the cell do not exceed the maximum.
///
/// Duplica la posicion de la particula indicada aplicandole un desplazamiento.
/// Las particulas duplicadas se considera que siempre son validas y estan dentro
/// del dominio.
/// Este kernel vale para single-cpu y multi-cpu porque los calculos se hacen 
/// a partir de domposmin.
/// Se controla que las coordendas de celda no sobrepasen el maximo.
//==============================================================================
void JSphCpuSingle::PeriodicDuplicatePos(unsigned pnew,unsigned pcopy,bool inverse
  ,double dx,double dy,double dz,tuint3 cellmax,tdouble3* pos,unsigned* dcell)const
{
  //-Get pos of particle to be duplicated. | Obtiene pos de particula a duplicar.
  tdouble3 ps=pos[pcopy];
  //-Apply displacement. | Aplica desplazamiento.
  ps.x+=(inverse? -dx: dx);
  ps.y+=(inverse? -dy: dy);
  ps.z+=(inverse? -dz: dz);
  //-Calculate coordinates of cell inside of domain. | Calcula coordendas de celda dentro de dominio.
  unsigned cx=unsigned((ps.x-DomPosMin.x)/Scell);
  unsigned cy=unsigned((ps.y-DomPosMin.y)/Scell);
  unsigned cz=unsigned((ps.z-DomPosMin.z)/Scell);
  //-Adjust coordinates of cell is they exceed maximum. | Ajusta las coordendas de celda si sobrepasan el maximo.
  cx=(cx<=cellmax.x? cx: cellmax.x);
  cy=(cy<=cellmax.y? cy: cellmax.y);
  cz=(cz<=cellmax.z? cz: cellmax.z);
  //-Record position and cell of new particles. |  Graba posicion y celda de nuevas particulas.
  pos[pnew]=ps;
  dcell[pnew]=DCEL_Cell(DomCellCode,cx,cy,cz);
}

//==============================================================================
/// Create periodic particles starting from a list of the particles to duplicate.
/// Assume that all the particles are valid.
/// This kernel works for single-cpu & multi-cpu because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-cpu y multi-cpu porque usa domposmin. 
//==============================================================================
void JSphCpuSingle::PeriodicDuplicateVerlet(unsigned np,unsigned pini
  ,tuint3 cellmax,tdouble3 perinc,const unsigned* listp,unsigned* idp
  ,typecode* code,unsigned* dcell,tdouble3* pos,tfloat4* velrho
  ,tsymatrix3f* spstau,tfloat4* velrhom1)const
{
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=0;p<n;p++){
    const unsigned pnew=unsigned(p)+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Adjust position and cell of new particle. | Ajusta posicion y celda de nueva particula.
    PeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000),perinc.x,perinc.y,perinc.z,cellmax,pos,dcell);
    //-Copy the rest of the values. | Copia el resto de datos.
    idp     [pnew]=idp[pcopy];
    code    [pnew]=CODE_SetPeriodic(code[pcopy]);
    velrho  [pnew]=velrho[pcopy];
    velrhom1[pnew]=velrhom1[pcopy];
    if(spstau)spstau[pnew]=spstau[pcopy];
  }
}

//==============================================================================
/// Create periodic particles starting from a list of the particles to duplicate.
/// Assume that all the particles are valid.
/// This kernel works for single-cpu & multi-cpu because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-cpu y multi-cpu porque usa domposmin. 
//==============================================================================
void JSphCpuSingle::PeriodicDuplicateSymplectic(unsigned np,unsigned pini
  ,tuint3 cellmax,tdouble3 perinc,const unsigned* listp,unsigned* idp
  ,typecode* code,unsigned* dcell,tdouble3* pos,tfloat4* velrho
  ,tsymatrix3f* spstau,tdouble3* pospre,tfloat4* velrhopre)const
{
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=0;p<n;p++){
    const unsigned pnew=unsigned(p)+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    //-Adjust position and cell of new particle. | Ajusta posicion y celda de nueva particula.
    PeriodicDuplicatePos(pnew,pcopy,(rp>=0x80000000),perinc.x,perinc.y,perinc.z,cellmax,pos,dcell);
    //-Copy the rest of the values. | Copia el resto de datos.
    idp   [pnew]=idp[pcopy];
    code  [pnew]=CODE_SetPeriodic(code[pcopy]);
    velrho[pnew]=velrho[pcopy];
    if(pospre)   pospre   [pnew]=pospre[pcopy];
    if(velrhopre)velrhopre[pnew]=velrhopre[pcopy];
    if(spstau)   spstau   [pnew]=spstau[pcopy];
  }
}

//==============================================================================
/// Create periodic particles starting from a list of the particles to duplicate.
/// Assume that all the particles are valid.
/// This kernel works for single-cpu & multi-cpu because it uses domposmin.
///
/// Crea particulas periodicas a partir de una lista con las particulas a duplicar.
/// Se presupone que todas las particulas son validas.
/// Este kernel vale para single-cpu y multi-cpu porque usa domposmin. 
//==============================================================================
void JSphCpuSingle::PeriodicDuplicateNormals(unsigned np,unsigned pini
  ,tuint3 cellmax,tdouble3 perinc,const unsigned* listp,tfloat3* normals
  ,tfloat3* motionvel,tfloat3* motionace)const
{
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=0;p<n;p++){
    const unsigned pnew=unsigned(p)+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    normals[pnew]=normals[pcopy];
    if(motionvel)motionvel[pnew]=motionvel[pcopy]; //<vs_m2dbc>
    if(motionace)motionace[pnew]=motionace[pcopy]; //<vs_m2dbc>
  }
}

//==============================================================================
/// Saves particle index to access to the parent of periodic particles.
//==============================================================================
void JSphCpuSingle::PeriodicSaveParent(unsigned np,unsigned pini
  ,const unsigned* listp,unsigned* periparent)const
{
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=0;p<n;p++){
    const unsigned pnew=unsigned(p)+pini;
    const unsigned rp=listp[p];
    const unsigned pcopy=(rp&0x7FFFFFFF);
    periparent[pnew]=pcopy;
  }
}

//==============================================================================
/// Marks current periodic particles to be ignored.
/// Marca periodicas actuales para ignorar.
//==============================================================================
void JSphCpuSingle::PeriodicIgnore(unsigned np,typecode* code)const{
  for(unsigned p=0;p<np;p++){
    const typecode rcode=code[p];
    if(CODE_IsPeriodic(rcode))code[p]=CODE_SetOutIgnore(rcode);
  }
}

//==============================================================================
/// Create duplicate particles for periodic conditions.
/// Create new periodic particles and mark the old ones to be ignored.
/// New periodic particles are created from Np of the beginning, first the NpbPer
/// of the boundary and then the NpfPer fluid ones. The Np of the those leaving contains also the
/// new periodic ones.
///
/// Crea particulas duplicadas de condiciones periodicas.
/// Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
/// Las nuevas periodicas se situan a partir del Np de entrada, primero las NpbPer
/// de contorno y despues las NpfPer fluidas. El Np de salida contiene tambien las
/// nuevas periodicas.
//==============================================================================
void JSphCpuSingle::RunPeriodic(){
  Timersc->TmStart(TMC_SuPeriodic);
  if(PeriParent_c)PeriParent_c->Memset(255,Np);
  //-Stores the current number of periodic particles.
  //-Guarda numero de periodicas actuales.
  NpfPerM1=NpfPer;
  NpbPerM1=NpbPer;
  //-Marks current periodic particles to be ignored.
  //-Marca periodicas actuales para ignorar.
  PeriodicIgnore(Np,Code_c->ptr());
  //-Creates new periodic particles.
  //-Crea las nuevas periodicas.
  acuint listp("listp",Arrays_Cpu,true);
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
      for(unsigned cblock=0;cblock<2;cblock++){//-0:new periodic, 1:original particles. | 0:periodicas nuevas, 1:particulas originales
        const unsigned nper=(ctype? NpfPer: NpbPer); //-Number of new periodic particles of type to be processed. | Numero de periodicas nuevas del tipo a procesar.
        const unsigned pini2=(cblock? pini: Np-nper);
        const unsigned num2= (cblock? num:  nper);
        //-Repeats search if the available memory was insufficient and had to be increased.
        //-Repite la busqueda si la memoria disponible resulto insuficiente y hubo que aumentarla.
        bool run=true;
        while(run && num2){
          //-Reserve memory to create list of periodic particles.
          //-Reserva memoria para crear lista de particulas periodicas.
          listp.Reserve();
          const unsigned nmax=CpuParticlesSize-1; //-Maximmum number of particles that fit in the list. | Numero maximo de particulas que caben en la lista.
          //-Generates list of new periodic particles.
          if(Np>=0x80000000)Run_Exceptioon("The number of particles is too big.");//-Because the last bit is used to mark the direction in which a new periodic particle is created. | Porque el ultimo bit se usa para marcar el sentido en que se crea la nueva periodica.
          const unsigned count=PeriodicMakeList(num2,pini2,Stable
            ,nmax,perinc,Pos_c->cptr(),Code_c->cptr(),listp.ptr());
          //-Resizes the allocated memory for the particles if there is not sufficient space and repeats the search process.
          //-Redimensiona memoria para particulas si no hay espacio suficiente y repite el proceso de busqueda.
          if(count>nmax || !CheckCpuParticlesSize(count+Np)){
            listp.Free(); //-Avoids unnecessary copying of its data during resizing.
            Timersc->TmStop(TMC_SuPeriodic);
            const unsigned ndatacpu=Np;
            ResizeParticlesSizeData(ndatacpu,Np+count,Np+count,PERIODIC_OVERMEMORYNP,false);
            Timersc->TmStart(TMC_SuPeriodic);
          }
          else{
            run=false;
            //-Create new periodic particles duplicating the particles from the list
            //-Crea nuevas particulas periodicas duplicando las particulas de la lista.
            if(TStep==STEP_Verlet){
              PeriodicDuplicateVerlet(count,Np,DomCells,perinc,listp.cptr()
                ,Idp_c->ptr(),Code_c->ptr(),Dcell_c->ptr(),Pos_c->ptr(),Velrho_c->ptr()
                ,AC_PTR(SpsTauRho2_c),VelrhoM1_c->ptr());
            }
            if(TStep==STEP_Symplectic){
              if(PosPre_c->Active()!=VelrhoPre_c->Active())
                Run_Exceptioon("Symplectic data is invalid.");
              PeriodicDuplicateSymplectic(count,Np,DomCells,perinc,listp.cptr()
                ,Idp_c->ptr(),Code_c->ptr(),Dcell_c->ptr(),Pos_c->ptr(),Velrho_c->ptr()
                ,AC_PTR(SpsTauRho2_c),PosPre_c->ptr(),VelrhoPre_c->ptr());
            }
            if(UseNormals){
              PeriodicDuplicateNormals(count,Np,DomCells,perinc,listp.cptr()
                ,BoundNor_c->ptr(),AC_PTR(MotionVel_c),AC_PTR(MotionAce_c)); //<vs_m2dbc>
            }
            if(PeriParent_c){
              PeriodicSaveParent(count,Np,listp.cptr(),PeriParent_c->ptr());
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
  Timersc->TmStop(TMC_SuPeriodic);
}

//==============================================================================
/// Executes divide of particles in cells.
/// Ejecuta divide de particulas en celdas.
//==============================================================================
void JSphCpuSingle::RunCellDivide(bool updateperiodic){
  DivData=DivDataCpuNull();
  //-Creates new periodic particles and marks the old ones to be ignored.
  //-Crea nuevas particulas periodicas y marca las viejas para ignorarlas.
  if(updateperiodic && PeriActive)RunPeriodic();

  //-Initiates Divide process.
  CellDivSingle->Divide(Npb,Np-Npb-NpbPer-NpfPer,NpbPer,NpfPer
    ,BoundChanged,Dcell_c->cptr(),Code_c->cptr(),Idp_c->cptr()
    ,Pos_c->cptr(),Timersc);
  DivData=CellDivSingle->GetCellDivData();

  //-Sorts particle data. | Ordena datos de particulas.
  Timersc->TmStart(TMC_NlSortData);
  CellDivSingle->SortArray(Idp_c->ptr());
  CellDivSingle->SortArray(Code_c->ptr());
  CellDivSingle->SortArray(Dcell_c->ptr());
  CellDivSingle->SortArray(Pos_c->ptr());
  CellDivSingle->SortArray(Velrho_c->ptr());
  if(TStep==STEP_Verlet){
    CellDivSingle->SortArray(VelrhoM1_c->ptr());
  }
  else if(TStep==STEP_Symplectic && (PosPre_c->Active() || VelrhoPre_c->Active())){//-In reality, this is only necessary in divide for corrector, not in predictor??? | En realidad solo es necesario en el divide del corrector, no en el predictor???
    if(!PosPre_c->Active() || !VelrhoPre_c->Active())
      Run_Exceptioon("Symplectic data is invalid.") ;
    CellDivSingle->SortArray(PosPre_c->ptr());
    CellDivSingle->SortArray(VelrhoPre_c->ptr());
  }
  if(TVisco==VISCO_LaminarSPS){
    CellDivSingle->SortArray(SpsTauRho2_c->ptr());
  }
  if(UseNormals){
    CellDivSingle->SortArray(BoundNor_c->ptr());
    if(MotionVel_c || MotionAce_c){ //<vs_m2dbc_ini>
      CellDivSingle->SortArray(MotionVel_c->ptr());
      CellDivSingle->SortArray(MotionAce_c->ptr());
    } //<vs_m2dbc_end>
  }
  if(ShiftingAdv){ //<vs_advshift_ini>
    CellDivSingle->SortArray(FSType_c->ptr());
    CellDivSingle->SortArray(ShiftVel_c->ptr());
  } //<vs_advshift_end>
  if(PeriParent_c){
    CellDivSingle->SortArrayPeriParent(PeriParent_c->ptr());
  }
  if(FlexStruc&&FlexStrucRidpc)CellDivSingle->UpdateIndices(CaseNflexstruc,FlexStrucRidpc); //<vs_flexstruc>

  #ifdef AVAILABLE_DIVCLEAN
    if(DivCleaning){
      CellDivSingle->SortArray(PsiClean_c->ptr());
    }
  #endif

  //-Collect divide data. | Recupera datos del divide.
  Np=CellDivSingle->GetNpFinal();
  Npb=CellDivSingle->GetNpbFinal();
  NpbOk=Npb-CellDivSingle->GetNpbIgnore();

  //-Manages excluded particles fixed, moving and floating before aborting the execution.
  if(CellDivSingle->GetNpbOut())AbortBoundOut();

  //-Collects index of moving and floating particles.
  if(CaseNmoving || CaseNfloat){
    const unsigned np=(CaseNmoving? Npb: 0) + (CaseNfloat? Np-Npb: 0);
    const unsigned pini=(CaseNmoving? 0: Npb);
    CalcRidp(PeriActive!=0,np,pini,CaseNfixed,CaseNbound
      ,Code_c->cptr(),Idp_c->cptr(),RidpMot);
  }
  Timersc->TmStop(TMC_NlSortData);

  //-Control of excluded particles (only fluid because excluded boundary were checked before).
  if(CellDivSingle->GetNpfOut()){
    Timersc->TmStart(TMC_NlOutCheck);
    SaveFluidOut();
    Timersc->TmStop(TMC_NlOutCheck);
  }
  BoundChanged=false;
}

//==============================================================================
/// Manages excluded particles fixed, moving and floating before aborting the execution.
/// Gestiona particulas excluidas fixed, moving y floating antes de abortar la ejecucion.
//==============================================================================
void JSphCpuSingle::AbortBoundOut(){
  const unsigned nboundout=CellDivSingle->GetNpbOut();
  //-Get data of excluded boundary particles.
  acuint     idp("idp",Arrays_Cpu,true);
  acdouble3  pos("pos",Arrays_Cpu,true);
  acfloat3   vel("vel",Arrays_Cpu,true);
  acfloat    rho("rho",Arrays_Cpu,true);
  actypecode cod("cod",Arrays_Cpu,true);
  unsigned nfilter=0;
  GetParticlesData(nboundout,Np,false,idp.ptr(),pos.ptr()
    ,vel.ptr(),rho.ptr(),cod.ptr(),NULL,nfilter);
  //-Shows excluded particles information and aborts execution.
  JSph::AbortBoundOut(Log,nboundout,idp.cptr(),pos.cptr(),vel.cptr()
    ,rho.cptr(),cod.cptr());
}

//==============================================================================
/// Manages the excluded fluid particles.
/// Gestiona las particulas fluid excluidas.
//==============================================================================
void JSphCpuSingle::SaveFluidOut(){
  const unsigned npfout=CellDivSingle->GetNpfOut();
  //-Get data of excluded fluid particles.
  acuint     idp("idp",Arrays_Cpu,true);
  acdouble3  pos("pos",Arrays_Cpu,true);
  acfloat3   vel("vel",Arrays_Cpu,true);
  acfloat    rho("rho",Arrays_Cpu,true);
  actypecode cod("cod",Arrays_Cpu,true);
  unsigned nfilter=0;
  GetParticlesData(npfout,Np,false,idp.ptr(),pos.ptr()
    ,vel.ptr(),rho.ptr(),cod.ptr(),NULL,nfilter);
  //-Stores new excluded particles until recordering next PART.
  AddParticlesOut(npfout,idp.cptr(),pos.cptr(),vel.cptr()
    ,rho.cptr(),cod.cptr());
}

//<vs_advshift_ini>
//==============================================================================
/// PreLoop for additional models computation.
//==============================================================================
void JSphCpuSingle::PreLoopProcedure(TpInterStep interstep){
  const bool runshift=(ShiftingAdv && interstep==INTERSTEP_SymPredictor && Nstep!=0);
  if(runshift){
    Timersc->TmStart(TMC_SuShifting);
    ComputeFSParticles();
    ComputeUmbrellaRegion();
    PreLoopInteraction_ct(DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
      ,Velrho_c->cptr(),FSType_c->ptr(),ShiftVel_c->ptr(),FSNormal_c->ptr()
      ,FSMinDist_c->ptr());
    ComputeShiftingVel(Simulate2D,ShiftingAdv->GetShiftCoef()
      ,ShiftingAdv->GetAleActive(),SymplecticDtPre,FSType_c->ptr()
      ,FSNormal_c->ptr(),FSMinDist_c->ptr(),ShiftVel_c->ptr());
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
    if(0)DgSaveVtkParticlesCpu("Compute_FreeSurface_",Part,0,Np,Pos_c->cptr()
      ,Code_c->cptr(),FSType_c->cptr(),ShiftVel_c->cptr(),FSNormal_c->cptr());
    Timersc->TmStop(TMC_SuShifting);
  }
}

//==============================================================================
/// Compute free-surface particles and their normals.
//==============================================================================
void JSphCpuSingle::ComputeFSParticles(){
  acuint fspart("-",Arrays_Cpu,true);
  CallComputeFSNormals(DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,Velrho_c->cptr(),FSType_c->ptr(),FSNormal_c->ptr(),fspart.ptr());
}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
//==============================================================================
void JSphCpuSingle::ComputeUmbrellaRegion(){
  acuint fspart("-",Arrays_Cpu,true);
  CallScanUmbrellaRegion(DivData,Dcell_c->cptr(),Pos_c->cptr(),Code_c->cptr()
    ,FSNormal_c->cptr(),fspart.ptr(),FSType_c->ptr());
}
//<vs_advshift_end>

//==============================================================================
/// Interaction to calculate forces.
/// Interaccion para el calculo de fuerzas.
//==============================================================================
void JSphCpuSingle::Interaction_Forces(TpInterStep interstep){
  tfloat3* dengradcorr=NULL;

  Timersc->TmStart(TMC_CfForces);
  //-Interaction of Fluid-Fluid/Bound & Bound-Fluid (forces and DEM).
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
    ,AC_CPTR(PsiClean_c),AC_PTR(PsiCleanRhs_c)    //<vs_divclean>
    ,DivCleanKp,DivCleaning                       //<vs_divclean>
  );
  StInterResultc res;
  res.viscdt=0;
  JSphCpu::Interaction_Forces_ct(parms,res);

  //<vs_flexstruc_ini>
  //-Interaction flexible structure-flexible structure.
  if(FlexStruc){
    Timersc->TmStart(TMC_SuFlexStruc);
    JSphCpu::Interaction_ForcesFlexStruc(FlexStrucDtMax);
    Timersc->TmStop(TMC_SuFlexStruc);
  }
  //<vs_flexstruc_end>

  //-For 2-D simulations zero the 2nd component.
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

  #ifdef AVAILABLE_DIVCLEAN
  CsPsiCleanMax=res.cspsiclean;
  #endif
  //-Calculates maximum value of Ace (periodic particles are ignored).
  AceMax=ComputeAceMax();

  InterNum++;
  Timersc->TmStop(TMC_CfForces);
}

//==============================================================================
/// Calculates extrapolated data on boundary particles from fluid domain for mDBC.
/// Calcula datos extrapolados en el contorno para mDBC.
//==============================================================================
void JSphCpuSingle::MdbcBoundCorrection(TpInterStep interstep){
  const bool runmdbc=(TBoundary==BC_MDBC 
    && (MdbcCorrector || interstep!=INTERSTEP_SymCorrector));
  if(runmdbc){
    Timersc->TmStart(TMC_CfPreMDBC);
    if(SlipMode==SLIP_Vel0){
      Interaction_MdbcCorrection(DivData,Pos_c->cptr(),Code_c->cptr()
        ,Idp_c->cptr(),BoundNor_c->cptr(),Velrho_c->ptr());
    }
    else{ //if(SlipMode==SLIP_NoSlip){ //<vs_m2dbc_ini>
      const unsigned nmode=(UseNormalsFt? Np: Npb);
      BoundMode_c->Reserve();       //-BoundMode_c is freed in PosInteraction_Forces().
      BoundMode_c->Memset(0,nmode); //-BoundMode_c[]=0=BMODE_DBC
      TangenVel_c->Reserve();       //-TangenVel_c is freed in PosInteraction_Forces().
      Interaction_Mdbc2Correction(DivData,Pos_c->cptr(),Code_c->cptr()
        ,Idp_c->cptr(),BoundNor_c->cptr(),MotionVel_c->cptr(),MotionAce_c->cptr()
        ,Velrho_c->ptr(),BoundMode_c->ptr(),TangenVel_c->ptr());
    } //<vs_m2dbc_end>
   // else Run_Exceptioon("Error: SlipMode is invalid.");
    Timersc->TmStop(TMC_CfPreMDBC);
  }
}


//==============================================================================
/// Returns maximum value of ace (modulus), periodic and inout particles must be ignored.
/// Devuelve el valor maximo de ace (modulo), se deben ignorar las particulas periodicas e inout.
//==============================================================================
double JSphCpuSingle::ComputeAceMax()const{
  const bool check=(PeriActive!=0 || InOut!=NULL);
  const unsigned pini=(CaseNflexstruc? 0: Npb);
  if(check)return(ComputeAceMaxOmp<true >(Np-pini,Ace_c->cptr()+pini,Code_c->cptr()+pini));
  else     return(ComputeAceMaxOmp<false>(Np-pini,Ace_c->cptr()+pini,Code_c->cptr()+pini));
}

//==============================================================================
/// Returns maximum value of ace (modulus), periodic particles must be ignored.
/// Devuelve el valor maximo de ace (modulo), se deben ignorar las particulas periodicas.
//==============================================================================
template<bool checkcode> double JSphCpuSingle::ComputeAceMaxSeq(unsigned np
  ,const tfloat3* ace,const typecode* code)const
{
  float acemax=0;
  const int n=int(np);
  for(int p=0;p<n;p++){
    const typecode cod=(checkcode? code[p]: 0);
    const tfloat3 a=(!checkcode || (CODE_IsNormal(cod) && !CODE_IsFluidInout(cod))? ace[p]: TFloat3(0));
    const float a2=a.x*a.x+a.y*a.y+a.z*a.z;
    acemax=max(acemax,a2);
  }
  return(sqrt(double(acemax)));
}

//==============================================================================
/// Returns maximum value of ace (modulus) using OpenMP, periodic particles must be ignored.
/// Devuelve el valor maximo de ace (modulo) using OpenMP, se deben ignorar las particulas periodicas.
//==============================================================================
template<bool checkcode> double JSphCpuSingle::ComputeAceMaxOmp(unsigned np
  ,const tfloat3* ace,const typecode* code)const
{
  double acemax=0;
  #ifdef OMP_USE
    if(np>OMP_LIMIT_COMPUTELIGHT){
      const int n=int(np);
      if(n<0)Run_Exceptioon("Number of values is too big.");
      float amax=0;
      #pragma omp parallel 
      {
        float amax2=0;
        #pragma omp for nowait
        for(int p=0;p<n;++p){
          const typecode cod=(checkcode? code[p]: 0);
          const tfloat3 a=(!checkcode || (CODE_IsNormal(cod) && !CODE_IsFluidInout(cod))? ace[p]: TFloat3(0));
          const float a2=a.x*a.x+a.y*a.y+a.z*a.z;
          if(amax2<a2)amax2=a2;
        }
        #pragma omp critical 
        {
          if(amax<amax2)amax=amax2;
        }
      }
      //-Saves result.
      acemax=sqrt(double(amax));
    }
    else if(np)acemax=ComputeAceMaxSeq<checkcode>(np,ace,code);
  #else
    if(np)acemax=ComputeAceMaxSeq<checkcode>(np,ace,code);
  #endif
  return(acemax);
}

//<vs_ddramp_ini>
//==============================================================================
/// Applies initial DDT ramp.
//==============================================================================
void JSphCpuSingle::RunInitialDDTRamp(){
  if(TimeStep<DDTRamp.x){
    if((Nstep%10)==0){//-DDTkh value is updated every 10 calculation steps.
      if(TimeStep<=DDTRamp.y)DDTkh=KernelSize * float(DDTRamp.z);
      else{
        const double tt=TimeStep-DDTRamp.y;
        const double tr=DDTRamp.x-DDTRamp.y;
        DDTkh=KernelSize * float(((tr-tt)/tr)*(DDTRamp.z-DDTValue)+DDTValue);
      }
    }
  }
  else{
    if(DDTkh!=DDTkhCte)CSP.ddtkh=DDTkh=DDTkhCte;
    DDTRamp.x=0;
  }
}//<vs_ddramp_end>

//==============================================================================
/// Perform interactions and updates of particles according to forces 
/// calculated in the interaction using Verlet.
///
/// Realiza interaccion y actualizacion de particulas segun las fuerzas 
/// calculadas en la interaccion usando Verlet.
//==============================================================================
double JSphCpuSingle::ComputeStep_Ver(){
  InterStep=INTERSTEP_Verlet;
  MdbcBoundCorrection(InterStep);        //-mDBC correction
  PreInteraction_Forces(InterStep);      //-Allocating temporary arrays.
  Interaction_Forces(InterStep);         //-Interaction.
  const double dt=DtVariable(true);      //-Calculate new dt.
  if(CaseNmoving)CalcMotion(dt);         //-Calculate motion for moving bodies.
  DemDtForce=dt;                         //-For DEM interaction.
  if(Shifting)RunShifting(dt);           //-Shifting.
  ComputeVerlet(dt);                     //-Update particles using Verlet.
  if(CaseNfloat)RunFloating(dt,false);   //-Control of floating bodies.
  PosInteraction_Forces();               //-Free memory used for interaction.
  if(Damping)RunDamping(dt);             //-Applies Damping.
  if(RelaxZones)RunRelaxZone(dt);        //-Generate waves using RZ.
  return(dt);
}

//==============================================================================
/// Perform interactions and updates of particles according to forces 
/// calculated in the interaction using Symplectic.
///
/// Realiza interaccion y actualizacion de particulas segun las fuerzas 
/// calculadas en la interaccion usando Symplectic.
//==============================================================================
double JSphCpuSingle::ComputeStep_Sym(){
  const double dt=SymplecticDtPre;
  if(CaseNmoving)CalcMotion(dt);          //-Calculate motion for moving bodies.
  //-Predictor
  //-----------
  InterStep=INTERSTEP_SymPredictor;
  DemDtForce=dt*0.5f;
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
  if(Shifting)RunShifting(dt);            //-Standard shifting.
  ComputeSymplecticCorr(dt);              //-Apply Symplectic-Corrector to particles (periodic particles become invalid).
  if(CaseNfloat)RunFloating(dt,false);    //-Control of floating bodies.
  PosInteraction_Forces();                //-Free memory used for interaction.
  if(Damping)RunDamping(dt);              //-Applies Damping.
  if(RelaxZones)RunRelaxZone(dt);         //-Generate waves using RZ.
  SymplecticDtPre=min(dt_p,dt_c);         //-Calculate dt for next ComputeStep.
  return(dt);
}

//==============================================================================
/// Process floating objects
/// Procesa floating objects.
//==============================================================================
void JSphCpuSingle::RunFloating(double dt,bool predictor){
  Timersc->TmStart(TMC_SuFloating);
  const bool saveftmot=(!predictor && (PartFloatSave && TimeOutExtraCheck(TimeStep+dt)) );
  const bool saveftvalues=(!predictor && (TimeStep+dt>=TimePartNext || saveftmot));
  const bool ftpaused=(TimeStep<FtPause);

  if(!ftpaused || saveftvalues){
    //-Computes sum of linear and angular acceleration of floating particles.
    FtPartsSumAce(Pos_c->cptr(),Ace_c->cptr(),RidpMot,Fto_AceLinAng);
    //-Compute new linear and angular acceleration, velocity and center to update floatings.
    FtComputeAceVel(dt,predictor,saveftvalues,Fto_AceLinAng,Fto_VelLinAng,Fto_Center);
  }

  //-Run floatings with Chrono library (change computed center, linear and angular velocity).
  if(!ftpaused && ChronoObjects){
    Timersc->TmStop(TMC_SuFloating);
    Timersc->TmStart(TMC_SuChrono);
    FtComputeChrono(dt,predictor,Fto_AceLinAng,Fto_VelLinAng,Fto_Center);
    Timersc->TmStop(TMC_SuChrono);
    Timersc->TmStart(TMC_SuFloating);
  }

  if(!ftpaused){//-Operator >= is used because when FtPause=0 in symplectic-predictor, code would not enter here. | Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
    //-Apply displacement and new velocity to floating particles.
    const bool updatenormals=(!predictor && UseNormalsFt);
    FtPartsUpdate(dt,updatenormals,Fto_VelLinAng,Fto_Center,RidpMot
      ,Pos_c->ptr(),Velrho_c->ptr(),Dcell_c->ptr(),Code_c->ptr()
      ,AC_PTR(BoundNor_c),AC_PTR(MotionVel_c),AC_PTR(MotionAce_c));

    //-Update floating data (FtObjs[]) for next step.
    if(!predictor)FtUpdateFloatings(dt,Fto_VelLinAng,Fto_Center);
  }

  //-Saves current floating data for output.
  if(saveftvalues){
    if(PartFloatSave)PartFloatSave->SetFtData(Part,TimeStep+dt,Nstep+1,FtObjs,ForcePoints);
  }
  Timersc->TmStop(TMC_SuFloating);

  //-Update data of points in FtForces and calculates motion data of affected floatings.
  if(!predictor && ForcePoints){
    Timersc->TmStart(TMC_SuMoorings);
    ForcePoints->UpdatePoints(TimeStep,dt,ftpaused,FtObjs);
    if(Moorings)Moorings->ComputeForces(Nstep,TimeStep,dt,ForcePoints);
    ForcePoints->ComputeForcesSum();
    Timersc->TmStop(TMC_SuMoorings);
  }
}

//==============================================================================
/// Calculate distance between floating particles & centre according to periodic conditions.
/// Calcula distancia entre pariculas floatin y centro segun condiciones periodicas.
//==============================================================================
tfloat3 JSphCpuSingle::FtPeriodicDist(const tdouble3& pos,const tdouble3& center
  ,float radius)const
{
  tdouble3 distd=(pos-center);
  if(PeriX && fabs(distd.x)>radius){
    if(distd.x>0)distd=distd+PeriXinc;
    else distd=distd-PeriXinc;
  }
  if(PeriY && fabs(distd.y)>radius){
    if(distd.y>0)distd=distd+PeriYinc;
    else distd=distd-PeriYinc;
  }
  if(PeriZ && fabs(distd.z)>radius){
    if(distd.z>0)distd=distd+PeriZinc;
    else distd=distd-PeriZinc;
  }
  return(ToTFloat3(distd));
}

//==============================================================================
/// Calculate summation of linear and angular acceleration of floating particles.
/// Calcula suma de aceleracion lineal y angular de las particulas floating.
//==============================================================================
void JSphCpuSingle::FtPartsSumAce(const tdouble3* posc,const tfloat3* acec
  ,const unsigned* ridpmot,tfloat6* acelinang)const
{
  const int ftcount=int(FtCount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int cf=0;cf<ftcount;cf++){
    const StFloatingData fobj=FtObjs[cf];
    const unsigned fpini=fobj.begin-CaseNfixed;
    const unsigned fpfin=fpini+fobj.count;
    const float    fradius=fobj.radius;
    const tdouble3 fcenter=fobj.center;
    tfloat3 acelin=TFloat3(0);
    tfloat3 aceang=TFloat3(0);
    for(unsigned fp=fpini;fp<fpfin;fp++){
      const unsigned p=ridpmot[fp];
      const tfloat3 acep=acec[p];
      acelin=acelin+acep;
      const tfloat3 dist=(PeriActive? FtPeriodicDist(posc[p],fcenter,fradius): ToTFloat3(posc[p]-fcenter)); 
      aceang.x+= acep.z*dist.y - acep.y*dist.z;
      aceang.y+= acep.x*dist.z - acep.z*dist.x;
      aceang.z+= acep.y*dist.x - acep.x*dist.y;
    }
    acelinang[cf]=TFloat6(acelin,aceang);
  }
}

//==============================================================================
/// Apply displacement and new velocity to floating particles.
//==============================================================================
void JSphCpuSingle::FtPartsUpdate(double dt,bool updatenormals
  ,const tfloat6* fto_vellinang,const tdouble3* fto_center,const unsigned* ridpmot
  ,tdouble3* posc,tfloat4* velrhoc,unsigned* dcellc,typecode* codec
  ,tfloat3* boundnorc,tfloat3* motionvelc,tfloat3* motionacec)const
{
  const int ftcount=int(FtCount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int cf=0;cf<ftcount;cf++){
    //-Get Floating object values.
    const StFloatingData fobj=FtObjs[cf];
    const float fradius=fobj.radius;
    const unsigned fpini=fobj.begin-CaseNfixed;
    const unsigned fpfin=fpini+fobj.count;
    const tfloat3  fvel   =fto_vellinang[cf].getlo();
    const tfloat3  fomega =fto_vellinang[cf].gethi();
    const tdouble3 fcenter=fto_center[cf];
    //-Compute matrix to update floating normals for mDBC.
    JMatrix4d mat;
    if(updatenormals){
      const tdouble3 dang=(ToTDouble3(fomega)*dt)*TODEG;
      const tdouble3 cen2=(PeriActive? UpdatePeriodicPos(fcenter): fcenter);
      mat.Move(cen2);
      mat.Rotate(dang);
      mat.Move(fobj.center*-1);
    }
    //-Updates floating particles.
    for(unsigned fp=fpini;fp<fpfin;fp++){
      const int p=ridpmot[fp];
      if(p!=UINT_MAX){
        tfloat4 vr=velrhoc[p];
        //-Computes displacement and updates position.
        const double dx=dt*double(vr.x);
        const double dy=dt*double(vr.y);
        const double dz=dt*double(vr.z);
        UpdatePos(posc[p],dx,dy,dz,false,p,posc,dcellc,codec);
        //-Computes and updates velocity.
        const tfloat3 dist=(PeriActive? FtPeriodicDist(posc[p],fcenter,fradius): ToTFloat3(posc[p]-fcenter)); 
        vr.x=fvel.x+(fomega.y*dist.z - fomega.z*dist.y);
        vr.y=fvel.y+(fomega.z*dist.x - fomega.x*dist.z);
        vr.z=fvel.z+(fomega.x*dist.y - fomega.y*dist.x);
        velrhoc[p]=vr;
        //-Updates motionvel and motionace for mDBC no-slip.
        if(TMdbc2>=MDBC2_Std){ //<vs_m2dbc_ini>
          const tfloat3 mvel0=motionvelc[p];
          motionacec[p]=TFloat3(float((double(vr.x)-mvel0.x)/dt),
                                float((double(vr.y)-mvel0.y)/dt),
                                float((double(vr.z)-mvel0.z)/dt));
          motionvelc[p]=TFloat3(vr.x,vr.y,vr.z);
        } //<vs_m2dbc_end>
        //-Updates floating normals for mDBC.
        if(updatenormals){
          boundnorc[p]=ToTFloat3(mat.MulNormal(ToTDouble3(boundnorc[p])));
        }
      }
    }
  }
}

//==============================================================================
/// Runs first calculations in configured gauges.
/// Ejecuta primeros calculos en las posiciones de medida configuradas.
//==============================================================================
void JSphCpuSingle::RunFirstGaugeSystem(double timestep){
  GaugeSystem->ConfigArraysCpu(Pos_c,Code_c,Idp_c,Velrho_c);
  GaugeSystem->CalculeCpu(timestep,DivData,NpbOk,Npb,Np,true);
}

//==============================================================================
/// Runs calculations in configured gauges.
/// Ejecuta calculos en las posiciones de medida configuradas.
//==============================================================================
void JSphCpuSingle::RunGaugeSystem(double timestep){
  if(GaugeSystem->GetCount()){
    Timersc->TmStart(TMC_SuGauges);
    //const bool svpart=(TimeStep>=TimePartNext);
    GaugeSystem->CalculeCpu(timestep,DivData,NpbOk,Npb,Np,false);
    Timersc->TmStop(TMC_SuGauges);
  }
}

//==============================================================================
/// Compute PIPS information of current particles.
/// Calcula datos de PIPS de particulas actuales.
//==============================================================================
void JSphCpuSingle::ComputePips(bool run){
  if(run || DsPips->CheckRun(Nstep)){
    TimerSim.Stop();
    const double timesim=TimerSim.GetSecs();
    DsPips->ComputeCpu(Nstep,TimeStep,timesim,CSP,OmpThreads,Np,Npb
      ,NpbOk,DivData,Dcell_c->cptr(),Pos_c->cptr());
  }
}

//==============================================================================
/// Initialises execution of simulation.
/// Inicia ejecucion de simulacion.
//==============================================================================
void JSphCpuSingle::Run(std::string appname,const JSphCfgRun* cfg,JLog2* log){
  if(!cfg || !log)return;
  AppName=appname; Log=log; CfgRun=cfg;

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
      if(PartMotionSave)PartMotionSave->AddDataExtraCpu(Part,TimeStep,Nstep,Np
        ,Pos_c->cptr(),RidpMot);
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
void JSphCpuSingle::SaveData(){
  const bool save=(SvData!=SDAT_None && SvData!=SDAT_Info);
  const unsigned npnormal=Np-NpbPer-NpfPer; //-Subtracts the periodic particles if they exist. | Resta las periodicas si las hubiera.
  unsigned npsave=npnormal;
  Timersc->TmStart(TMC_SuSavePart);
  //-Collect particle values in original order. | Recupera datos de particulas en orden original.
  acuint    svidp("svidp",Arrays_Cpu,save);
  acdouble3 svpos("svpos",Arrays_Cpu,save);
  acfloat3  svvel("svvel",Arrays_Cpu,save);
  acfloat   svrho("svrho",Arrays_Cpu,save);
  if(save){
    //-Prepare filter for output particles data. //<vs_outpaarts>
    acbyte filter("filter",Arrays_Cpu,false);
    //<vs_outpaarts_ini>
    const bool svextra=(SvExtraDataBi4 && SvExtraDataBi4->CheckSave(Part));
    if(OutputParts && OutputParts->CheckFilters(Part) && !svextra){
      filter.Reserve();
      OutputParts->UpdateFtPos(FtCount,FtObjs);
      OutputParts->ComputeFilterCpu(Np,Idp_c->cptr(),Pos_c->cptr()
        ,Code_c->cptr(),filter.ptr());
    }//<vs_outpaarts_end>
    //-Obtain output particles data.
    unsigned npfilterdel=0;
    const unsigned npsel=GetParticlesData(Np,0,PeriActive!=0
      ,svidp.ptr(),svpos.ptr(),svvel.ptr(),svrho.ptr(),NULL
      ,filter.cptr(),npfilterdel);
    if(npsel+npfilterdel!=npnormal)Run_Exceptioon("The number of particles is invalid.");
    npsave=npsel;
  }
  //-Saves main motion reference data from particles (moving and floating bodies).
  if(PartMotionSave){
    PartMotionSave->SaveDataMainCpu(Part,TimeStep,Nstep,Np,Pos_c->cptr()
      ,RidpMot);
    PartMotionSave->SaveDataExtra();
  }

  //-Collects additional information.
  StInfoPartPlus infoplus;
  if(SvData&SDAT_Info){
    infoplus.gpudata=false;
    infoplus.SetNct(CellDivSingle->GetNct(),CellDivSingle->GetSizeNct());
    infoplus.SetNp(Np,CpuParticlesSize,npnormal,npsave);
    if(InOut)infoplus.npnew=InOut->GetNewNpPart();
    infoplus.SetNpExtra(NpbOk,Npb-NpbOk,Np-Npb,NpbPer,NpfPer);
    infoplus.memorycpualloc=GetAllocMemoryCpu();
    TimerSim.Stop();
    infoplus.timesim=TimerSim.GetSecs();
  }
  else infoplus.SetBasic(Np,npnormal,Np-Npb,CellDivSingle->GetNct());

  //-Obtains current domain limits.
  const tdouble6 vdom=CellDivSingle->GetDomainLimitsMinMax();
  //-Stores particle data. | Graba datos de particulas.
  JDataArrays arrays;
  AddBasicArrays(arrays,npsave,svpos.cptr(),svidp.cptr(),svvel.cptr(),svrho.cptr());
  JSph::SaveData(npsave,arrays,1,&vdom,infoplus);
  //-Save VTK file with current boundary normals (for debug).
  if(UseNormals && SvNormals)SaveVtkNormals(DirVtkOut+"Normals.vtk",Part
    ,npsave,Npb,Pos_c->cptr(),Idp_c->cptr(),BoundNor_c->cptr(),1.f);
  //-Save extra data.
  if(SvExtraDataBi4)SaveExtraData();
  Timersc->TmStop(TMC_SuSavePart);
}

//==============================================================================
/// Displays and stores final summary of the execution.
/// Muestra y graba resumen final de ejecucion.
//==============================================================================
void JSphCpuSingle::SaveExtraData(){
  const bool svextra=(BoundNor_c!=NULL);
  if(svextra && SvExtraDataBi4->CheckSave(Part)){
    SvExtraDataBi4->InitPartData(Part,TimeStep,Nstep);
    //-Saves normals of mDBC.
    if(BoundNor_c){
      SvExtraDataBi4->AddNormals(UseNormalsFt,Np,Npb,Idp_c->cptr()
        ,(PeriActive? Code_c->cptr(): NULL),BoundNor_c->cptr());
    }
    //-Saves file.
    SvExtraDataBi4->SavePartData();
  }
}

//==============================================================================
/// Displays and stores final summary of the execution.
/// Muestra y graba resumen final de ejecucion.
//==============================================================================
void JSphCpuSingle::FinishRun(bool stop){
  const float tsim=(float)TimerSim.GetSecs();
  const float ttot=(float)TimerTot.GetSecs();
  JSph::ShowResume(stop,tsim,ttot,true,"");
  Log->Print(" ");
  string hinfo,dinfo;
  if(SvTimers){
    Timersc->ShowTimes("[CPU Timers]",Log);
    Timersc->GetTimersInfo(hinfo,dinfo);
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
void JSphCpuSingle::FlexStrucInit(){
  //-Start timer and print info.
  Timersc->TmStart(TMC_SuFlexStruc);
  Log->Print("\nInitialising Flexible Structures...");
  //-Allocate array.
  FlexStrucDatac=new StFlexStrucData[FlexStrucCount]; MemCpuFixed+=(sizeof(StFlexStrucData)*FlexStrucCount);
  //-Get flexible structure data for each body and copy to GPU
  for(unsigned c=0;c<FlexStrucCount;c++){
    std::vector<typecode> clampcode=FlexStruc->GetBody(c)->GetClampCode();
    FlexStrucDatac[c].nc=unsigned(clampcode.size());
    std::copy(clampcode.begin(),clampcode.end(),FlexStrucDatac[c].clampcode);
    FlexStrucDatac[c].vol0=FlexStruc->GetBody(c)->GetParticleVolume();
    FlexStrucDatac[c].rho0=FlexStruc->GetBody(c)->GetDensity();
    FlexStrucDatac[c].youngmod=FlexStruc->GetBody(c)->GetYoungMod();
    FlexStrucDatac[c].poissratio=FlexStruc->GetBody(c)->GetPoissRatio();
    FlexStrucDatac[c].hgfactor=FlexStruc->GetBody(c)->GetHgFactor();
    FlexStrucDatac[c].cmat=FlexStruc->GetBody(c)->GetConstitMatrix();
  }
  //-Configure code for flexible structures.
  FlexStruc->ConfigCode(Npb,Code_c->ptr());
  JSphCpu::SetFlexStrucClampCodes(Npb,Pos_c->cptr(),FlexStrucDatac,Code_c->ptr());
  //-Count number of flexible structure particles.
  CaseNflexstruc=JSphCpu::CountFlexStrucParts(Npb,Code_c->cptr());
  //-Allocate arrays.
  FlexStrucRidpc          =new unsigned[CaseNflexstruc];  MemCpuFixed+=(sizeof(unsigned)*CaseNflexstruc);
  Pos0c                   =new tdouble3[CaseNflexstruc];  MemCpuFixed+=(sizeof(tdouble3)*CaseNflexstruc);
  NumPairsc               =new unsigned[CaseNflexstruc];  MemCpuFixed+=(sizeof(unsigned)*CaseNflexstruc);
  PairIdxc                =new unsigned*[CaseNflexstruc]; MemCpuFixed+=(sizeof(unsigned*)*CaseNflexstruc);
  KerCorrc                =new tmatrix3f[CaseNflexstruc]; MemCpuFixed+=(sizeof(tmatrix3f)*CaseNflexstruc);
  DefGradc                =new tmatrix3f[CaseNflexstruc]; MemCpuFixed+=(sizeof(tmatrix3f)*CaseNflexstruc);
  if(UseNormals)BoundNor0c=new tfloat3[CaseNflexstruc];   MemCpuFixed+=(sizeof(tfloat3)*CaseNflexstruc);
  //-Calculate array for indexing into flexible structure particles.
  JSphCpu::CalcFlexStrucRidp(Npb,Code_c->cptr(),FlexStrucRidpc);
  //-Copy current position and normals into initial position and normals.
  JSphCpu::GatherToFlexStrucArray(CaseNflexstruc,FlexStrucRidpc,Pos_c->cptr(),Pos0c);
  if(UseNormals)JSphCpu::GatherToFlexStrucArray(CaseNflexstruc,FlexStrucRidpc,BoundNor_c->cptr(),BoundNor0c);
  //-Get number of particle pairs for each flexible structure particle.
  const unsigned numpairstot=JSphCpu::CountFlexStrucPairs(CaseNflexstruc,Pos0c,NumPairsc);
  //-Allocate memory for raw buffer for storing pair indices and set the pointers to the indices.
  PairIdxBufferc=new unsigned[numpairstot]; MemCpuFixed+=(sizeof(unsigned)*numpairstot);
  unsigned* offset=PairIdxBufferc;
  vector<unsigned*> pairidx(CaseNflexstruc);
  for(unsigned p=0;p<CaseNflexstruc;p++){
    PairIdxc[p]=offset;
    offset+=NumPairsc[p];
  }
  //-Set the indices for each particle pair.
  JSphCpu::SetFlexStrucPairs(CaseNflexstruc,Pos0c,PairIdxc);
  //-Calculate kernel correction and update geometry.
  JSphCpu::CalcFlexStrucKerCorr();
  JSphCpu::UpdateFlexStrucGeometry();
  Log->Print("");
  Timersc->TmStop(TMC_SuFlexStruc);
}

//==============================================================================
/// Updates the geometric information for each flexible structure particle.
/// Actualiza la información geométrica de cada partícula de estructura flexible.
//==============================================================================
void JSphCpuSingle::UpdateFlexStrucGeometry(){
  Timersc->TmStart(TMC_SuFlexStruc);
  JSphCpu::UpdateFlexStrucGeometry();
  Timersc->TmStop(TMC_SuFlexStruc);
}
//<vs_flexstruc_end>
