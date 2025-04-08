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

/// \file JSphGpu.cpp \brief Implements the class \ref JSphGpu.

#include "JSphGpu.h"
#include "JException.h"
#include "JSphGpu_ker.h"
#include "JSphGpu_mdbc_iker.h"
#include "JSphGpuSimple_ker.h"
#include "JCellDivGpu.h"
#include "JDsGpuInfo.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include "JDsMotion.h"
#include "JDsFixedDt.h"
#include "JDsSaveDt.h"
#include "JDsOutputTime.h"
#include "JWaveGen.h"
#include "JMLPistons.h"
#include "JRelaxZones.h"
#include "JChronoObjects.h"
#include "JDsFtForcePoints.h"
#include "JDsDamping.h"
#include "JDsAccInput.h"
#include "JXml.h"
#include "JDsGaugeSystem.h"
#include "JSphInOut.h"
#include "JSphShifting.h"
#include "JDataArrays.h"
#include "JSpVtkData.h"

#include <climits>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JSphGpu::JSphGpu(bool withmpi):JSph(1,withmpi),DivAxis(MGDIV_None){
  ClassName="JSphGpu";
  CellDiv=NULL;
  Arrays_Cpu=NULL;
  Arrays_Gpu=NULL;
  GpuInfo=new JDsGpuInfo;
  Timersg=new JDsTimersGpu;
  InitVars();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphGpu::~JSphGpu(){
  DestructorActive=true;
  FreeCpuMemoryParticles();
  FreeGpuMemoryParticles();
  FreeGpuMemoryFixed();
  FreeCpuMemoryFixed();
  delete Arrays_Cpu; Arrays_Cpu=NULL;
  delete Arrays_Gpu; Arrays_Gpu=NULL;
  delete GpuInfo;    GpuInfo=NULL;
  delete Timersg;    Timersg=NULL;
  // cudaDeviceReset();
}

//==============================================================================
/// Throws exception related to a CUDA error from a static method.
//==============================================================================
void JSphGpu::RunExceptioonCudaStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,cudaError_t cuerr,std::string msg)
{
  msg=msg+fun::PrintStr(" (CUDA error %d (%s)).\n",cuerr,cudaGetErrorString(cuerr));
  throw JException(srcfile,srcline,"JSphGpu",method,msg,"");
}

//==============================================================================
/// Checks CUDA error and throws exception from a static method.
//==============================================================================
void JSphGpu::CheckCudaErroorStatic(const std::string& srcfile,int srcline
  ,const std::string& method,std::string msg)
{
  const cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess)RunExceptioonCudaStatic(srcfile,srcline,method,cuerr,msg);
}

//==============================================================================
/// Throws exception related to a CUDA error.
//==============================================================================
void JSphGpu::RunExceptioonCuda(const std::string& srcfile,int srcline
  ,const std::string& classname,const std::string& method
  ,cudaError_t cuerr,std::string msg)const
{
  msg=msg+fun::PrintStr(" (CUDA error %d (%s)).\n",cuerr,cudaGetErrorString(cuerr));
  throw JException(srcfile,srcline,classname,method,msg,"");
}

//==============================================================================
/// Checks CUDA error and throws exception.
/// Comprueba error de CUDA y lanza excepcion si lo hubiera.
//==============================================================================
void JSphGpu::CheckCudaErroor(const std::string& srcfile,int srcline
  ,const std::string& classname,const std::string& method
  ,std::string msg)const
{
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess)RunExceptioonCuda(srcfile,srcline,classname,method,cuerr,msg);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphGpu::InitVars(){
  RunMode="";
  memset(&BlockSizes,0,sizeof(StBlockSizes));
  BlockSizesStr="";

  DivData=DivDataGpuNull();

  Np=Npb=NpbOk=0;
  NpfPer=0;
  NpbPer=0;
  NpfPerM1=0;
  NpbPerM1=0;

  Idp_c=NULL;
  Code_c=NULL;
  Dcell_c=NULL;
  Posxy_c=NULL;
  Posz_c=NULL;
  Velrho_c=NULL;

  AuxPos_c=NULL;
  AuxVel_c=NULL;
  AuxRho_c=NULL;
  FreeCpuMemoryParticles();

  Idp_g=NULL;
  Code_g=NULL;
  Dcell_g=NULL;
  Posxy_g=NULL;
  Posz_g=NULL;
  PosCell_g=NULL;
  Velrho_g=NULL;

  PeriParent_g=NULL;

  BoundNor_g=NULL;     //-mDBC
  MotionVel_g=NULL;    //-mDBC2  //<vs_m2dbc>
  MotionAce_g=NULL;    //-mDBC2  //<vs_m2dbc>
  BoundMode_g=NULL;    //-mDBC2  //<vs_m2dbc>
  TangenVel_g=NULL;    //-mDBC2  //<vs_m2dbc>

  VelrhoM1_g=NULL;     //-Verlet
  PosxyPre_g=NULL;     //-Symplectic
  PoszPre_g=NULL;      //-Symplectic
  VelrhoPre_g=NULL;    //-Symplectic

  ViscDt_g=NULL; 
  Ace_g=NULL;
  Ar_g=NULL;
  Delta_g=NULL;
  NoPenShift_g = NULL; //<vs_m2dbcNP>

  SpsTauRho2_g=NULL;   //-Laminar+SPS.
  Sps2Strain_g=NULL;   //-Laminar+SPS.

  ShiftPosfs_g=NULL;   //-Shifting.

  ShiftVel_g=NULL;     //-ShiftingAdvanced //<vs_advshift>
  FSType_g=NULL;       //-ShiftingAdvanced //<vs_advshift>
  FSMinDist_g=NULL;    //-ShiftingAdvanced //<vs_advshift>
  FSNormal_g=NULL;     //-ShiftingAdvanced //<vs_advshift>

  //<vs_flexstruc_ini>
  FlexStrucDatag=NULL;
  FlexStrucRidpg=NULL;
  PosCell0g=NULL;
  NumPairsg=NULL;
  PairIdxBufferg=NULL;
  PairIdxg=NULL;
  KerCorrg=NULL;
  DefGradg=NULL;
  BoundNor0g=NULL;
  FlexStrucDtg=NULL;
  FlexStrucDtMax=0;
  //<vs_flexstruc_end>

  PsiClean_g=NULL;       //<vs_divclean>
  PsiCleanPre_g=NULL;    //<vs_divclean>
  PsiCleanRhs_g=NULL;    //<vs_divclean>
  CsPsiClean_g=NULL;     //<vs_divclean>
  CsPsiCleanMax=0;       //<vs_divclean>
  
  FreeGpuMemoryParticles();

  //-Variables for moving and floating bodies (GPU memory).
  RidpMotg=NULL;
  FtoMasspg=NULL;
  FtoDatpg=NULL;
  FtoCenterg=NULL;
  FtoAceg=NULL;
  NStmFloatings=0;
  for(unsigned c=0;c<MaxNStmFloatings;c++)StmFloatings[c]=NULL;
  //-DEM.
  DemDatag=NULL;
  FreeGpuMemoryFixed();

  //-Variables for floating bodies (CPU memory).
  FtoCenterc=NULL;
  FreeCpuMemoryFixed();
}

//==============================================================================
/// Frees fixed memory on CPU for moving and floating bodies.
/// Libera memoria fija en CPU para moving y floating.
//==============================================================================
void JSphGpu::FreeCpuMemoryFixed(){
  MemCpuFixed=0;
  delete[] FtoCenterc; FtoCenterc=NULL;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphGpu::AllocCpuMemoryFixed(){
  MemCpuFixed=0;
  try{
    //-Allocates memory for floating bodies.
    if(CaseNfloat){
      FtoCenterc=new tdouble3[FtCount];  MemCpuFixed+=(sizeof(tdouble3)*FtCount);
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Frees fixed memory on the GPU for moving and floating bodies.
/// Libera memoria fija en GPU para moving y floating.
//==============================================================================
void JSphGpu::FreeGpuMemoryFixed(){
  MemGpuFixed=0;
  //-Memory for moving and floating particles.
  if(RidpMotg)  cudaFree(RidpMotg);    RidpMotg=NULL;
  //-Memory for floating bodies.
  if(FtoMasspg) cudaFree(FtoMasspg);   FtoMasspg=NULL;
  if(FtoDatpg)  cudaFree(FtoDatpg);    FtoDatpg=NULL;
  if(FtoCenterg)cudaFree(FtoCenterg);  FtoCenterg=NULL;
  if(FtoAceg)   cudaFree(FtoAceg);     FtoAceg=NULL;
  //-Memory for DEM coefficients.
  if(DemDatag)  cudaFree(DemDatag);    DemDatag=NULL;
  //<vs_flexstruc_ini>
  if(FlexStrucDatag)    cudaFree(FlexStrucDatag);     FlexStrucDatag=NULL;
  if(FlexStrucRidpg)    cudaFree(FlexStrucRidpg);     FlexStrucRidpg=NULL;
  if(PosCell0g)         cudaFree(PosCell0g);          PosCell0g=NULL;
  if(NumPairsg)         cudaFree(NumPairsg);          NumPairsg=NULL;
  if(PairIdxBufferg)    cudaFree(PairIdxBufferg);     PairIdxBufferg=NULL;
  if(PairIdxg)          cudaFree(PairIdxg);           PairIdxg=NULL;
  if(KerCorrg)          cudaFree(KerCorrg);           KerCorrg=NULL;
  if(DefGradg)          cudaFree(DefGradg);           DefGradg=NULL;
  if(BoundNor0g)        cudaFree(BoundNor0g);         BoundNor0g=NULL;
  if(FlexStrucDtg)      cudaFree(FlexStrucDtg);       FlexStrucDtg=NULL;
  //<vs_flexstruc_end>
  //-Frees streams for floating bodies.
  for(unsigned c=0;c<MaxNStmFloatings;c++){
    if(StmFloatings[c])cudaStreamDestroy(StmFloatings[c]);
    StmFloatings[c]=NULL;
  }
  NStmFloatings=0;
}

//==============================================================================
/// Allocates memory for arrays with fixed size (motion and floating bodies).
//==============================================================================
void JSphGpu::AllocGpuMemoryFixed(){
  MemGpuFixed=0;
  //-Allocates memory for moving and floating particles.
  if(CaseNmoving || CaseNfloat){
    MemGpuFixed+=fcuda::Malloc(&RidpMotg,CaseNmoving+CaseNfloat);
  }
  //-Allocates memory for floating bodies.
  if(CaseNfloat){
    MemGpuFixed+=fcuda::Malloc(&FtoMasspg ,FtCount);
    MemGpuFixed+=fcuda::Malloc(&FtoDatpg  ,FtCount);
    MemGpuFixed+=fcuda::Malloc(&FtoCenterg,FtCount);
    MemGpuFixed+=fcuda::Malloc(&FtoAceg   ,FtCount*2);
  }
  //-Allocates streams for floating bodies.
  NStmFloatings=(FtCount>1? min(MaxNStmFloatings,FtCount): 0);
  for(unsigned c=0;c<NStmFloatings;c++){
    cudaStreamCreate(StmFloatings+c);
  }
  //-Allocates memory for DEM coefficients.
  if(UseDEM){
    MemGpuFixed+=fcuda::Malloc(&DemDatag,DemDataSize);
  }
  Check_CudaErroor("Failed GPU memory allocation.");
}

//==============================================================================
/// Frees CPU memory for the particles.
/// Libera memoria en CPU para particulas.
//==============================================================================
void JSphGpu::FreeCpuMemoryParticles(){
  //-Free array objects.
  delete Idp_c;      Idp_c=NULL;
  delete Code_c;     Code_c=NULL;
  delete Dcell_c;    Dcell_c=NULL;
  delete Posxy_c;    Posxy_c=NULL;
  delete Posz_c;     Posz_c=NULL;
  delete Velrho_c;   Velrho_c=NULL;

  delete AuxPos_c;   AuxPos_c=NULL;
  delete AuxVel_c;   AuxVel_c=NULL;
  delete AuxRho_c;   AuxRho_c=NULL;
  
  //-Free CPU memory for array objects.
  CpuParticlesSize=0;
  if(Arrays_Cpu)Arrays_Cpu->Reset();
}

//==============================================================================
/// Allocates memory for the main particle data.
/// Reserva memoria para datos principales de particulas.
//==============================================================================
void JSphGpu::AllocCpuMemoryParticles(unsigned np){
  FreeCpuMemoryParticles();
  //-Calculate number of partices to allocate memory.
  if(np>0)np=np+PARTICLES_OVERMEMORY_MIN;
  CpuParticlesSize=np;
  //-Set size of arrays.
  Arrays_Cpu->SetArraySize(CpuParticlesSize);
  //-Create arrays for basic particle data.
  Idp_c   =new acuint    ("Idpc"   ,Arrays_Cpu,true);
  Code_c  =new actypecode("Codec"  ,Arrays_Cpu,true);
  Dcell_c =new acuint    ("Dcellc" ,Arrays_Cpu,true);
  Posxy_c =new acdouble2 ("Posxyc" ,Arrays_Cpu,true);
  Posz_c  =new acdouble  ("Poszc"  ,Arrays_Cpu,true);
  Velrho_c=new acfloat4  ("Velrhoc",Arrays_Cpu,true);

  AuxPos_c=new acdouble3 ("AuxPosc",Arrays_Cpu,true);
  AuxVel_c=new acfloat3  ("AuxVelc",Arrays_Cpu,true);
  AuxRho_c=new acfloat   ("AuxRhoc",Arrays_Cpu,true);
}

//==============================================================================
/// Frees GPU memory for the particles.
/// Libera memoria en GPU para particulas.
//==============================================================================
void JSphGpu::FreeGpuMemoryParticles(){
  //-Free array objects.
  delete Idp_g;         Idp_g=NULL;
  delete Code_g;        Code_g=NULL;
  delete Dcell_g;       Dcell_g=NULL;
  delete Posxy_g;       Posxy_g=NULL;
  delete Posz_g;        Posz_g=NULL;
  delete PosCell_g;     PosCell_g=NULL;
  delete Velrho_g;      Velrho_g=NULL;

  delete PeriParent_g;  PeriParent_g=NULL;

  delete BoundNor_g;    BoundNor_g=NULL;    //-mDBC
  delete MotionVel_g;   MotionVel_g=NULL;   //-mDBC2  //<vs_m2dbc>
  delete MotionAce_g;   MotionAce_g=NULL;   //-mDBC2  //<vs_m2dbc>
  delete BoundMode_g;   BoundMode_g=NULL;   //-mDBC2  //<vs_m2dbc>
  delete TangenVel_g;   TangenVel_g=NULL;   //-mDBC2  //<vs_m2dbc>

  delete VelrhoM1_g;    VelrhoM1_g=NULL;    //-Verlet
  delete PosxyPre_g;    PosxyPre_g=NULL;    //-Symplectic
  delete PoszPre_g;     PoszPre_g=NULL;     //-Symplectic
  delete VelrhoPre_g;   VelrhoPre_g=NULL;   //-Symplectic

  delete ViscDt_g;      ViscDt_g=NULL; 
  delete Ace_g;         Ace_g=NULL;
  delete Ar_g;          Ar_g=NULL;
  delete Delta_g;       Delta_g=NULL;
  delete NoPenShift_g;  NoPenShift_g=NULL;  //-NoPenetration //<vs_m2dbcNP>

  delete SpsTauRho2_g;  SpsTauRho2_g=NULL;  //-Laminar+SPS.
  delete Sps2Strain_g;  Sps2Strain_g=NULL;  //-Laminar+SPS.

  delete ShiftPosfs_g;  ShiftPosfs_g=NULL;  //-Shifting.

  delete ShiftVel_g;    ShiftVel_g=NULL;    //-ShiftingAdvanced //<vs_advshift>
  delete FSType_g;      FSType_g=NULL;      //-ShiftingAdvanced //<vs_advshift>
  delete FSMinDist_g;   FSMinDist_g=NULL;   //-ShiftingAdvanced //<vs_advshift>
  delete FSNormal_g;    FSNormal_g=NULL;    //-ShiftingAdvanced //<vs_advshift>

  delete PsiClean_g;    PsiClean_g=NULL;    //<vs_divclean>
  delete PsiCleanPre_g; PsiCleanPre_g=NULL; //<vs_divclean>
  delete PsiCleanRhs_g; PsiCleanRhs_g=NULL; //<vs_divclean>
  delete CsPsiClean_g;  CsPsiClean_g=NULL;  //<vs_divclean>

  //-Free GPU memory for array objects.
  GpuParticlesSize=0;
  if(Arrays_Gpu)Arrays_Gpu->Reset();
}

//==============================================================================
/// Allocates GPU memory for the particles.
/// Reserva memoria en Gpu para las particulas. 
//==============================================================================
void JSphGpu::AllocGpuMemoryParticles(unsigned np){
  FreeGpuMemoryParticles();
  //-Calculate number of partices to allocate memory.
  GpuParticlesSize=np+PARTICLES_OVERMEMORY_MIN;
  //-Checks maximum number.
  if(GpuParticlesSize>=unsigned(INT_MAX))
    Run_Exceptioon("Particle number for allocation is too large.");
  //-Set size of arrays.
  Arrays_Gpu->SetArraySize(GpuParticlesSize);
  //-Create arrays for basic particle data.
  Idp_g    =new aguint    ("Idpg"    ,Arrays_Gpu,true);
  Code_g   =new agtypecode("Codeg"   ,Arrays_Gpu,true);
  Dcell_g  =new aguint    ("Dcellg"  ,Arrays_Gpu,true);
  Posxy_g  =new agdouble2 ("Posxyg"  ,Arrays_Gpu,true);
  Posz_g   =new agdouble  ("Poszg"   ,Arrays_Gpu,true);
  PosCell_g=new agfloat4  ("PosCellg",Arrays_Gpu,true);
  Velrho_g =new agfloat4  ("Velrhog" ,Arrays_Gpu,true);
  //-Arrays for mDBC.
  if(UseNormals){
    BoundNor_g=new agfloat3("BoundNorg",Arrays_Gpu,true);
    if(SlipMode>=SLIP_NoSlip){ //<vs_m2dbc_ini>
      MotionVel_g=new agfloat3("MotionVelg",Arrays_Gpu,true);
      MotionAce_g=new agfloat3("MotionAceg",Arrays_Gpu,true);
      BoundMode_g=new agbyte  ("BoundModeg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
      TangenVel_g=new agfloat3("TangenVelg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
    } //<vs_m2dbc_end>
  }
  //-Arrays for Verlet.
  if(TStep==STEP_Verlet){
    VelrhoM1_g=new agfloat4  ("VelrhoM1g",Arrays_Gpu,true);
  }
  //-Arrays for Symplectic.
  if(TStep==STEP_Symplectic){
    PosxyPre_g =new agdouble2("PosxyPreg" ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
    PoszPre_g  =new agdouble ("PoszPreg"  ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
    VelrhoPre_g=new agfloat4 ("VelrhoPreg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  }
  //-Arrays for forces computation.
  ViscDt_g    =new agfloat ("ViscDtg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  Ace_g       =new agfloat3("Aceg"   ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
  Ar_g        =new agfloat ("Arg"    ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
  Delta_g     =new agfloat ("Deltag" ,Arrays_Gpu,false); //-NO INITIAL MEMORY.
  //-Arrays for No Pentration.
  if(TMdbc2==MDBC2_NoPen){
    NoPenShift_g = new agfloat4("NoPenShiftg", Arrays_Gpu, false); //-NO INITIAL MEMORY.
  }
  //-Arrays for Laminar+SPS.
  if(TVisco==VISCO_LaminarSPS){
    SpsTauRho2_g=new agsymatrix3f("SpsTauRho2g",Arrays_Gpu,true);
    Sps2Strain_g=new agsymatrix3f("Sps2Straing",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  }
  //-Arrays for Shifting.
  ShiftPosfs_g=new agfloat4("ShiftPosfsg",Arrays_Gpu,false); //-NO INITIAL MEMORY.
  //-Arrays for Advanced shifting. //<vs_advshift_ini>
  if(ShiftingAdv!=NULL){
    ShiftVel_g    =new agfloat4 ("ShiftVel_g" ,Arrays_Gpu,true); 
    FSType_g      =new aguint   ("FSType_g"   ,Arrays_Gpu,true);
    FSNormal_g    =new agfloat3 ("FSNormal_g" ,Arrays_Gpu,false);  //-NO INITIAL MEMORY.
    FSMinDist_g   =new agfloat  ("FSMinDist_g",Arrays_Gpu,false);  //-NO INITIAL MEMORY.
    if(PeriActive && !PeriParent_g){
      PeriParent_g=new aguint(  "PeriParentg" ,Arrays_Gpu,true);
    }
  }//<vs_advshift_end>

  //<vs_divclean_ini>
  if(DivCleaning){
    PsiClean_g    =new agfloat("PhiClean_g"   ,Arrays_Gpu,true);
    PsiCleanPre_g =new agfloat("PhiCleanPre_g",Arrays_Gpu,false);
    PsiCleanRhs_g =new agfloat("PhiCleanRhs_g",Arrays_Gpu,false);
    CsPsiClean_g  =new agfloat("CsPhiClean_g" ,Arrays_Gpu,false); 
  }
  //<vs_divclean_end>

  //-Check CUDA errors.
  Check_CudaErroor("Failed GPU memory allocation.");
}

//==============================================================================
/// Resizes GPU memory for particles saving current data (ndatagpu).
//==============================================================================
void JSphGpu::ResizeGpuMemoryParticlesData(unsigned ndatagpu
  ,unsigned npnew,unsigned npmin)
{
  npnew=npnew+PARTICLES_OVERMEMORY_MIN;
  //-Shows GPU memory allocation.
  const llong memgpuparticles=Arrays_Gpu->GetAllocMemoryGpu();
  const double mbparticle=(double(memgpuparticles)/MEBIBYTE)/GpuParticlesSize; //-MiB per particle.
  const string txover=(npmin>1? fun::PrintStr(" (over-allocation: %.2fX)",double(npnew)/npmin): "");
  Log->Printf("**JSphGpu: Requesting GPU memory for %s particles%s: %.1f MiB (%u times)."
    ,KINT(npnew),txover.c_str(),mbparticle*npnew,Arrays_Gpu->GetCountResizeDataCalls()+1);
  //-Resizes GPU memory allocation.
  Arrays_Gpu->SetDataArraySize(npnew,ndatagpu);
  //-Updates values.
  GpuParticlesSize=npnew;
}

//==============================================================================
/// Returns the allocated memory on the CPU.
/// Devuelve la memoria reservada en CPU.
//==============================================================================
llong JSphGpu::GetAllocMemoryCpu()const{  
  llong s=JSph::GetAllocMemoryCpu();
  //-Allocated for particle arrays on CPU.
  if(Arrays_Cpu)s+=Arrays_Cpu->GetAllocMemoryCpu();
  //-Allocated in AllocCpuMemoryFixed().
  s+=MemCpuFixed;
  //-Allocated in other objects.
  if(MLPistons)s+=MLPistons->GetAllocMemoryCpu();
  return(s);
}

//==============================================================================
/// Returns the allocated memory in the GPU.
/// Devuelve la memoria reservada en GPU.
//==============================================================================
llong JSphGpu::GetAllocMemoryGpu()const{  
  llong s=0;
  //-Allocated for particle arrays on GPU.
  if(Arrays_Gpu)s+=Arrays_Gpu->GetAllocMemoryGpu();
  //-Allocated in AllocGpuMemoryFixed().
  s+=MemGpuFixed;
  //-Allocated in their objects.
  if(MLPistons)s+=MLPistons->GetAllocMemoryGpu();
  return(s);
}

//==============================================================================
/// Displays the allocated memory.
/// Visualiza la memoria reservada.
//==============================================================================
void JSphGpu::PrintAllocMemory(llong mcpu,llong mgpu)const{
  if(!Nstep)Log->Print("");
  Log->Printf("%s allocated memory in CPU: %s (%.2f MiB)  (%u particle arrays)"
    ,(!Nstep? "Initial": "**Updated"),KINT(mcpu),double(mcpu)/MEBIBYTE,Arrays_Cpu->GetArrayCount());
  Log->Printf("%s allocated memory in GPU: %s (%.2f MiB)  (%u particle arrays)"
    ,(!Nstep? "Initial": "**Updated"),KINT(mgpu),double(mgpu)/MEBIBYTE,Arrays_Gpu->GetArrayCount());
  Arrays_Cpu->GetArrayCountUpdated();
  Arrays_Gpu->GetArrayCountUpdated();
}

//==============================================================================
/// Loads constant data for GPU in ctes.
/// Carga constantes para GPU en ctes.
//==============================================================================
void JSphGpu::GetConstantData(StCteInteraction& ctes)const{
  memset(&ctes,0,sizeof(StCteInteraction));
  //-Set mass particle values.
  SetCtegMass(ctes,MassBound,MassFluid);
  //-Set distance values depending on h and dp.
  SetCtegKsize(ctes,KernelH,KernelSize,KernelSize2,PosCellSize
    ,Eta2,float(Dp),Scell,MovLimit);
  //-Wendland constants are always computed since this kernel is used in some parts where other kernels are not defined (e.g. mDBC, inlet/outlet...).
  SetCtegKerWendland(ctes,KWend);
  //-Copies constants for other kernels.
  if(TKernel==KERNEL_Cubic)SetCtegKerCubic(ctes,KCubic);
  //-Set values for density and pressure.
  SetCtegRho(ctes,RhopZero,1.f/RhopZero,Gamma,float(Cs0),CteB);
  //-Set values for DDT.
  SetCtegDdt(ctes,DDTkh,DDTgz);
  //-Set values on formulation options.
  SetCtegOpts(ctes,unsigned(TBoundary));
  //-Set values for open periodic boundaries.
  SetCtegPeriodic(ctes,PeriActive,PeriXinc,PeriYinc,PeriZinc);
  //-Set values for map definition.
  SetCtegMap(ctes,MapRealPosMin,MapRealSize);
  //-Set values depending on the assigned domain (can change).
  //-MGDIV_Z is necessary to avoid errors in KerGetInteraction_Cells().
  SetCtegDomain(ctes,MGDIV_Z,DomCellCode,DomPosMin);
}

//==============================================================================
/// Uploads constants to the GPU.
/// Copia constantes a la GPU.
//==============================================================================
void JSphGpu::ConstantDataUp(){
  StCteInteraction ctes;
  GetConstantData(ctes);
  //-Copy constants to GPU memory.
  cusph::CteInteractionUp(&ctes);
  Check_CudaErroor("Failed copying constants to GPU.");
}

//==============================================================================
/// Uploads particle data to the GPU.
/// Sube datos de particulas a la GPU.
//==============================================================================
void JSphGpu::ParticlesDataUp(unsigned n,const tfloat3* boundnor){
  Idp_g   ->CuCopyFromHost(Idp_c,n);
  Code_g  ->CuCopyFromHost(Code_c,n);
  Dcell_g ->CuCopyFromHost(Dcell_c,n);
  Posxy_g ->CuCopyFromHost(Posxy_c,n);
  Posz_g  ->CuCopyFromHost(Posz_c,n);
  Velrho_g->CuCopyFromHost(Velrho_c,n);
  if(UseNormals){
    BoundNor_g->CuCopyFromHost2(boundnor,n);
  }
  Check_CudaErroor("Failed copying data to GPU.");
}

//==============================================================================
/// Recovers particle data from the GPU and returns the particle number that
/// are less than n if the paeriodic particles are removed.
/// - code: Recovers data of Codeg.
/// - onlynormal: Only retains the normal particles, removes the periodic ones.
///
/// Recupera datos de particulas de la GPU y devuelve el numero de particulas que
/// sera menor que n si se eliminaron las periodicas.
/// - code: Recupera datos de Codeg.
/// - onlynormal: Solo se queda con las normales, elimina las particulas periodicas.
//==============================================================================
unsigned JSphGpu::ParticlesDataDown(unsigned n,unsigned pini,bool code
  ,bool onlynormal,const byte* filterg,unsigned& npfilterdel)
{
  unsigned num=n;
  //-Copy data from GPU memory.
  Idp_g   ->CuCopyToHostOffset(pini,Idp_c  ,0,n);
  Posxy_g ->CuCopyToHostOffset(pini,Posxy_c,0,n);
  Posz_g  ->CuCopyToHostOffset(pini,Posz_c,0,n);
  Velrho_g->CuCopyToHostOffset(pini,Velrho_c,0,n);
  if(code || onlynormal){
    Code_g->CuCopyToHostOffset(pini,Code_c,0,n);
  }
  //-Creates simple CPU pointers.
  unsigned* idpc   =Idp_c->ptr(); 
  tdouble2* posxyc =Posxy_c->ptr(); 
  double*   poszc  =Posz_c->ptr(); 
  tfloat4*  velrhoc=Velrho_c->ptr(); 
  typecode* codec  =Code_c->ptr();

  //-Obtain filter data on CPU memory. //<vs_outpaarts_ini>
  byte* filter=NULL; 
  if(filterg){
    filter=(byte*)AuxRho_c->ptr();
    cudaMemcpy(filter,filterg+pini,sizeof(byte)*n,cudaMemcpyDeviceToHost);
  }//<vs_outpaarts_end>
  Check_CudaErroor("Failed copying data from GPU.");
  
  //-Eliminate non-normal particles (periodic & others).
  const bool usefilter=(filter!=NULL);
  unsigned nfilter=0;
  if(onlynormal || usefilter){
    unsigned ndel=0;
    for(unsigned p=0;p<n;p++){
      const bool isnormal=(!onlynormal || CODE_IsNormal(codec[p]));
      const bool selected=(isnormal && (!usefilter || (filter[p]&1)));
      if(ndel && selected){
        const unsigned p2=p-ndel;
        idpc   [p2]=idpc   [p];
        posxyc [p2]=posxyc [p];
        poszc  [p2]=poszc  [p];
        velrhoc[p2]=velrhoc[p];
        codec  [p2]=codec  [p];
      }
      if(!selected){
        ndel++;
        if(isnormal)nfilter++;//-Normal particles removed by filters.
      }
    }
    num-=ndel;
  }
  npfilterdel=nfilter;
  filter=NULL;
  //-Converts data to a simple format in AuxPos_c ,AuxVel_c and AuxRhop_c.
  Pos21Vel4ToPos3Vel31(num,posxyc,poszc,velrhoc
    ,AuxPos_c->ptr(),AuxVel_c->ptr(),AuxRho_c->ptr());
  return(num);
}

//==============================================================================
/// Initialises CUDA device and returns GPU id.
/// Inicializa dispositivo CUDA y devuelve id de GPU.
//==============================================================================
int JSphGpu::SelecDevice(int gpuid){
  //-Shows information on available GPUs.
  JDsGpuInfo::ShowGpusInfo(Log);
  //-GPU selection.
  Log->Print("[GPU Hardware]");
  GpuInfo->SelectGpu(gpuid);
  GpuInfo->ShowSelectGpusInfo(Log);
  cudaSetDevice(GpuInfo->GetGpuId());
  return(GpuInfo->GetGpuId());
}

//==============================================================================
/// Configures BlockSize for main interaction CUDA kernels.
//==============================================================================
void JSphGpu::ConfigBlockSizes(bool usezone,bool useperi){
  Log->Print(" ");
  BlockSizesStr="";
  if(CellMode==CELLMODE_Full || CellMode==CELLMODE_Half){
    BlockSizes.forcesbound=BlockSizes.forcesfluid=BlockSizes.forcesdem=BSIZE_FORCES;
    //-Collects kernel information.
    StKerInfo kerinfo;
    memset(&kerinfo,0,sizeof(StKerInfo));
    StDivDataGpu divdatag;
    memset(&divdatag,0,sizeof(StDivDataGpu));
    #ifndef DISABLE_BSMODES
      const StInterParmsg parms=StrInterParmsg(Simulate2D
        ,TKernel,FtMode
        ,TVisco,TDensity,ShiftingMode,TMdbc2 //<vs_m2dbc>
        ,false,false,false,false //<vs_advshift>
        ,0,0,0,0,100,0,0
        ,0,0,divdatag,NULL
        ,NULL,NULL,NULL
        ,NULL,NULL,NULL
        ,NULL,NULL,NULL //<vs_m2dbc>
        ,NULL //<vs_m2dbc> //SHABA4
        ,NULL,NULL,NULL
        ,NULL,NULL,NULL,NULL
        ,NULL
        ,NULL
        ,NULL // SHABA
        ,NULL,NULL     //<vs_advshift>
        ,NULL,NULL,NULL,0,false //<vs_divclean>
        ,NULL,&kerinfo);
      cusph::Interaction_Forces(parms);
      if(UseDEM)cusph::Interaction_ForcesDem(BlockSizes.forcesdem,CaseNfloat,divdatag,NULL,NULL,NULL,NULL,0,NULL,NULL,NULL,NULL,NULL,NULL,&kerinfo);
    #endif
    //Log->Printf("====> bound -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesbound_rg,kerinfo.forcesbound_bs,kerinfo.forcesbound_bsmax);
    //Log->Printf("====> fluid -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesfluid_rg,kerinfo.forcesfluid_bs,kerinfo.forcesfluid_bsmax);
    //Log->Printf("====> dem   -> rg:%d  bs:%d  bsmax:%d",kerinfo.forcesdem_rg  ,kerinfo.forcesdem_bs  ,kerinfo.forcesdem_bsmax);
    Log->Printf("BlockSize calculation mode: Fixed.");
    const string txb=fun::PrintStr("BsForcesBound=%d ",BlockSizes.forcesbound)+(kerinfo.forcesbound_rg? fun::PrintStr("(%d regs)",kerinfo.forcesbound_rg): string("(? regs)"));
    const string txf=fun::PrintStr("BsForcesFluid=%d ",BlockSizes.forcesfluid)+(kerinfo.forcesfluid_rg? fun::PrintStr("(%d regs)",kerinfo.forcesfluid_rg): string("(? regs)"));
    const string txd=fun::PrintStr("BsForcesDem=%d "  ,BlockSizes.forcesdem)  +(kerinfo.forcesdem_rg  ? fun::PrintStr("(%d regs)",kerinfo.forcesdem_rg  ): string("(? regs)"));
    Log->Print(string("  ")+txb);
    Log->Print(string("  ")+txf);
    if(UseDEM)Log->Print(string("  ")+txd);
    if(!BlockSizesStr.empty())BlockSizesStr=BlockSizesStr+" - ";
    BlockSizesStr=BlockSizesStr+txb+" - "+txf;
    if(UseDEM)BlockSizesStr=BlockSizesStr+" - "+txd;
  }
  else Run_Exceptioon("CellMode unrecognised.");
  Log->Print(" ");
}

//==============================================================================
/// Configures execution mode in the GPU.
/// Configura modo de ejecucion en GPU.
//==============================================================================
void JSphGpu::ConfigRunMode(){
  Hardware=GpuInfo->GetHardware();
  //-Defines RunMode.
  RunMode="";
  if(Stable)RunMode=RunMode+(!RunMode.empty()? " - ": "")+"Stable";
  RunMode=RunMode+(!RunMode.empty()? " - ": "")+"Pos-Cell";
  RunMode=RunMode+(!RunMode.empty()? " - ": "")+"Single-GPU";
  //-Shows RunMode.
  Log->Print(" ");
  Log->Print(fun::VarStr("RunMode",RunMode));
  Log->Print(" ");
}

//==============================================================================
/// Adjusts variables of particles of floating bodies for GPU execution.
//==============================================================================
void JSphGpu::InitFloatingsGpu(float* ftomasspg,float4* ftodatag
  ,double3* ftocenterg,float4* demdatag)const
{
  //-Copies massp values to GPU.
  {
    float* massp=new float[FtCount];
    for(unsigned cf=0;cf<FtCount;cf++)massp[cf]=FtObjs[cf].massp;
    cudaMemcpy(ftomasspg,massp,sizeof(float)*FtCount,cudaMemcpyHostToDevice);
    delete[] massp; massp=NULL;
  }
  //-Copies floating values to GPU.
  {
    typedef struct{
      unsigned pini;
      unsigned np;
      float radius;
      float massp;
    }stdata;
    stdata* datp=new stdata[FtCount];
    tdouble3* centerc=new tdouble3[FtCount];
    for(unsigned cf=0;cf<FtCount;cf++){
      const StFloatingData& fobj=FtObjs[cf];
      datp[cf].pini=fobj.begin-CaseNfixed;
      datp[cf].np=fobj.count;
      datp[cf].radius=fobj.radius;
      datp[cf].massp=fobj.massp;
      centerc[cf]=fobj.center;
    }
    cudaMemcpy(ftodatag,datp,sizeof(float4)*FtCount,cudaMemcpyHostToDevice);
    cudaMemcpy(ftocenterg,centerc,sizeof(double3)*FtCount,cudaMemcpyHostToDevice);
    delete[] datp;    datp=NULL;
    delete[] centerc; centerc=NULL;
  }
  //-Copies data object for DEM to GPU.
  if(UseDEM){
    float4* ddata=new float4[DemDataSize];
    for(unsigned c=0;c<DemDataSize;c++){
      ddata[c].x=DemData[c].mass;
      ddata[c].y=DemData[c].tau;
      ddata[c].z=DemData[c].kfric;
      ddata[c].w=DemData[c].restitu;
    }
    cudaMemcpy(demdatag,ddata,sizeof(float4)*DemDataSize,cudaMemcpyHostToDevice);
    delete[] ddata; ddata=NULL;
  }
}

//==============================================================================
/// Initialises arrays and variables for the execution.
/// Inicializa vectores y variables para la ejecucion.
//==============================================================================
void JSphGpu::InitRunGpu(){
  unsigned nfilter=0;
  ParticlesDataDown(Np,0,false,false,NULL,nfilter);
  InitRun(Np,Idp_c->cptr(),AuxPos_c->cptr());
  if(CaseNfloat)InitFloatingsGpu(FtoMasspg,FtoDatpg,FtoCenterg,DemDatag);
  if(TStep==STEP_Verlet)VelrhoM1_g->CuCopyFrom(Velrho_g,Np);
  if(TVisco==VISCO_LaminarSPS)SpsTauRho2_g->CuMemset(0,Np);
  if(MotionVel_g)MotionVel_g->CuMemset(0,Np); //<vs_m2dbc>
  if(MotionAce_g)MotionAce_g->CuMemset(0,Np); //<vs_m2dbc>
  if(ShiftVel_g)ShiftVel_g->CuMemset(0,Np);   //<vs_advshift>
  if(ShiftVel_g)FSType_g->CuMemset(3,Np);     //<vs_advshift>
  if(PsiClean_g)PsiClean_g->CuMemset(0,Np);   //<vs_divclean>
  Check_CudaErroor("Failed initializing variables for execution.");
}

//==============================================================================
/// Prepares variables for interaction (pre-loop and forces).
/// Prepara variables para interaccion (pre-loop y fuerzas).
//==============================================================================
void JSphGpu::PreInteraction_Forces(TpInterStep interstep){
  Timersg->TmStart(TMG_CfPreForces,false);
  //-Assign memory.
  ViscDt_g->Reserve();
  Ar_g->Reserve();
  Ace_g->Reserve();
  if(DDTArray)Delta_g->Reserve();
  if(Shifting)ShiftPosfs_g->Reserve();
  if(TVisco==VISCO_LaminarSPS)Sps2Strain_g->Reserve();
  if(TMdbc2==MDBC2_NoPen)NoPenShift_g->Reserve(); // SHABA
  if(ShiftingAdv){ //<vs_advshift_ini>
    FSMinDist_g->Reserve();
    FSNormal_g->Reserve();
  } //<vs_advshift_end>

  //<vs_divclean_ini>
  if(DivCleaning){
    PsiCleanRhs_g->Reserve(); 
    CsPsiClean_g->Reserve();
    PsiCleanRhs_g->CuMemset(0,Np);
    CsPsiClean_g->CuMemset(0,Np);
    CsPsiCleanMax=0;
  }  
  //vs_divclean_end>

  //-Initialise arrays.
  const unsigned npf=Np-Npb;
  ViscDt_g->CuMemset(0,Np);                                         //ViscDtg[]=0
  Ar_g->CuMemset(0,Np);                                             //Arg[]=0
  Ace_g->CuMemset(0,Np);                                            //Aceg[]=(0)
  if(AG_CPTR(Delta_g))Delta_g->CuMemset(0,Np);                      //Deltag[]=0
  if(AG_CPTR(Sps2Strain_g))Sps2Strain_g->CuMemsetOffset(Npb,0,npf); //Sps2Straing[]=(0)
  if(AG_CPTR(NoPenShift_g))NoPenShift_g->CuMemset(0,Np);            //NoPenShiftg[]=(0) //<vs_m2dbcNP
  
  //-Select particles for shifting.
  if(AC_CPTR(ShiftPosfs_g))Shifting->InitGpu(npf,Npb,Posxy_g->cptr()
    ,Posz_g->cptr(),ShiftPosfs_g->ptr());

  //<vs_advshift_ini>
  if(AC_CPTR(ShiftVel_g) && interstep==INTERSTEP_SymPredictor)
    ShiftVel_g->CuMemset(0,Np);
  if(AC_CPTR(FSMinDist_g))FSMinDist_g->CuMemset(0,Np);
  if(AC_CPTR(FSNormal_g))FSNormal_g->CuMemset(0,Np);
  //<vs_advshift_end>

  //-Adds variable acceleration from input configuration.
  if(AccInput)AccInput->RunGpu(TimeStep,Gravity,npf,Npb,Code_g->cptr()
    ,Posxy_g->cptr(),Posz_g->cptr(),Velrho_g->cptr(),Ace_g->ptr());

  //-Computes VelMax: Includes the particles from floating bodies and does not affect the periodic conditions.
  //-Calcula VelMax: Se incluyen las particulas floatings y no afecta el uso de condiciones periodicas.
  const unsigned pini=(DtAllParticles? 0: Npb);
  cusph::ComputeVelMod(Np-pini,Velrho_g->cptr()+pini,ViscDt_g->ptr());
  const float velmax=cusph::ReduMaxFloat(Np-pini,0,ViscDt_g->ptr()
    ,CellDiv->GetAuxMem(cusph::ReduMaxFloatSize(Np-pini)));
  VelMax=sqrt(velmax);

  ViscDt_g->CuMemset(0,Np); //ViscDtg[]=0
  ViscDtMax=0;
  //<vs_flexstruc_ini>
  if(FlexStruc){
    cudaMemset(FlexStrucDtg,0,sizeof(float)*CaseNflexstruc);  //FlexStrucDtg[]=0
    FlexStrucDtMax=0;
  }
  //<vs_flexstruc_end>
  Timersg->TmStop(TMG_CfPreForces,true);
  Check_CudaErroor("Failed calculating VelMax.");
}

//==============================================================================
/// Frees memory allocated from Arrays_Gpu in PreInteraction_Forces().
/// Libera memoria asignada de Arrays_Gpu en PreInteraction_Forces().
//==============================================================================
void JSphGpu::PosInteraction_Forces(){
  ViscDt_g->Free();
  Ar_g->Free();
  Ace_g->Free();
  Delta_g->Free();
  ShiftPosfs_g->Free();
  if(Sps2Strain_g)Sps2Strain_g->Free();
  if(BoundMode_g)BoundMode_g->Free(); //-Reserved in MdbcBoundCorrection(). //<vs_m2dbc>
  if(TangenVel_g)TangenVel_g->Free(); //-Reserved in MdbcBoundCorrection(). //<vs_m2dbc>
  if(NoPenShift_g)NoPenShift_g->Free(); //<vs_m2dbcNP>
  if(FSMinDist_g)FSMinDist_g->Free(); //<vs_advshift>
  if(FSNormal_g)FSNormal_g->Free();   //<vs_advshift>
  if(PsiCleanRhs_g)PsiCleanRhs_g->Free();  //<vs_divclean>
  if(CsPsiClean_g)CsPsiClean_g->Free() ;   //<vs_divclean>
}

//==============================================================================
/// Updates particles according to forces and dt using Verlet.
/// Actualizacion de particulas segun fuerzas y dt usando Verlet.
//==============================================================================
void JSphGpu::ComputeVerlet(double dt){  //pdtedom
  Timersg->TmStart(TMG_SuComputeStep,false);
  const bool shift=(ShiftingMode!=SHIFT_None);
  const bool inout=(InOut!=NULL);
  const float3* indirvel=(inout? InOut->GetDirVelg(): NULL);
  const byte* boundmode=AG_CPTR(BoundMode_g); //<vs_m2dbc>
  VerletStep++;
  //-Allocates memory to compute the displacement.
  agdouble2 movxyg("movxyg",Arrays_Gpu,true);
  agdouble movzg("movzg",Arrays_Gpu,true);
  //-Computes displacement, velocity and density.
  if(VerletStep<VerletSteps){
    const double twodt=dt+dt;
    cusphs::ComputeStepVerlet(WithFloating,shift,inout,TMdbc2,Np,Npb
      ,Velrho_g->cptr(),VelrhoM1_g->cptr(),boundmode,Ar_g->cptr()
      ,Ace_g->cptr(),ShiftPosfs_g->cptr(),indirvel,AG_CPTR(NoPenShift_g),dt,twodt
      ,RhopZero,RhopOutMin,RhopOutMax,Gravity,Code_g->ptr()
      ,movxyg.ptr(),movzg.ptr(),VelrhoM1_g->ptr(),NULL);
  }
  else{
    cusphs::ComputeStepVerlet(WithFloating,shift,inout,TMdbc2,Np,Npb
      ,Velrho_g->cptr(),Velrho_g->cptr(),boundmode,Ar_g->cptr()
      ,Ace_g->cptr(),ShiftPosfs_g->cptr(),indirvel,AG_CPTR(NoPenShift_g),dt,dt
      ,RhopZero,RhopOutMin,RhopOutMax,Gravity,Code_g->ptr()
      ,movxyg.ptr(),movzg.ptr(),VelrhoM1_g->ptr(),NULL);
    VerletStep=0;
  }
  //-Applies displacement to non-periodic fluid particles.
  cusph::ComputeStepPos(PeriActive,WithFloating,Np,Npb,movxyg.cptr(),movzg.cptr()
    ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    #ifndef MDBC2_KEEPVEL
      if(MotionVel_g)cusphs::CopyMotionVelFlexStruc(CaseNflexstruc,Code_g->cptr(),FlexStrucRidpg,MotionVel_g->cptr(),Velrho_g->ptr());
    #endif
    cusphs::ComputeStepFlexStrucSemiImplicitEuler(CaseNflexstruc,Velrho_g->cptr(),Code_g->cptr(),FlexStrucRidpg,Ace_g->cptr(),dt,Gravity,movxyg.ptr(),movzg.ptr(),VelrhoM1_g->ptr(),NULL);
    cusph::ComputeStepPosFlexStruc(CaseNflexstruc,FlexStrucRidpg,Posxy_g->cptr(),Posz_g->cptr(),movxyg.cptr(),movzg.cptr(),Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
    if(!cusph::FlexStrucStepIsValid(Npb,Code_g->cptr()))
      Run_Exceptioon("Issue with FlexStruc particle update step (either naturally moved out of the domain or blew up).");
    BoundChanged=true;
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end>

  //-The new values are calculated in VelrhoM1_g.
  Velrho_g->SwapPtr(VelrhoM1_g); //-Exchanges Velrho_g and VelrhoM1_g.

  Timersg->TmStop(TMG_SuComputeStep,true);
}

//==============================================================================
/// Updates particles according to forces and dt using Symplectic-Predictor.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Predictor.
//==============================================================================
void JSphGpu::ComputeSymplecticPre(double dt){
  Timersg->TmStart(TMG_SuComputeStep,false);
  const bool shift=false; //(ShiftingMode!=SHIFT_None); //-We strongly recommend running the shifting correction only for the corrector. If you want to re-enable shifting in the predictor, change the value here to "true".
  const bool inout=(InOut!=NULL);
  //-Assign memory to PRE variables.
  PosxyPre_g->Reserve();
  PoszPre_g->Reserve();
  VelrhoPre_g->Reserve();
  if(DivCleaning)PsiCleanPre_g->Reserve();  //<vs_divclean>
  //-Move current data to PRE variables for calculating the new data.
  PosxyPre_g->SwapPtr(Posxy_g);   //- PosxyPre_g[]   <= Posxy_g[]
  PoszPre_g->SwapPtr(Posz_g);     //- PoszPre_g[]    <= Posz_g[]
  VelrhoPre_g->SwapPtr(Velrho_g); //- VelrhoPre_g[]  <= Velrho_g[]
  if(DivCleaning)PsiCleanPre_g->SwapPtr(PsiClean_g); //-PhiCleanPre[]<= PhiClean_g[] //<vs_divclean>
  //-Allocate memory to compute the diplacement.
  agdouble2 movxyg("movxyg",Arrays_Gpu,true);
  agdouble  movzg("movzg",Arrays_Gpu,true);
  //-Compute displacement, velocity and density.
  const double dt05=dt*.5;
  const float3* indirvel=(InOut? InOut->GetDirVelg(): NULL);
  const byte* boundmode=AG_CPTR(BoundMode_g); //<vs_m2dbc>
  cusphs::ComputeStepSymplecticPre(WithFloating,shift,inout,TMdbc2,Np,Npb
    ,VelrhoPre_g->cptr(),boundmode,Ar_g->cptr(),Ace_g->cptr(),ShiftPosfs_g->cptr()
    ,indirvel,dt05,RhopZero,RhopOutMin,RhopOutMax,Gravity
    ,Code_g->ptr(),movxyg.ptr(),movzg.ptr(),Velrho_g->ptr()
    ,AG_PTR(PsiClean_g),AG_CPTR(PsiCleanPre_g),AG_CPTR(PsiCleanRhs_g),DivCleaning,NULL);  //<vs_divclean>

  //-Applies displacement to non-periodic fluid particles.
  cusph::ComputeStepPos2(PeriActive,WithFloating,Np,Npb
    ,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr()
    ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
  //-Copy previous position of the boundary particles.
  Posxy_g->CuCopyFrom(PosxyPre_g,Npb);
  Posz_g->CuCopyFrom(PoszPre_g,Npb);

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    #ifndef MDBC2_KEEPVEL
      if(MotionVel_g)cusphs::CopyMotionVelFlexStruc(CaseNflexstruc,Code_g->cptr(),FlexStrucRidpg,MotionVel_g->cptr(),VelrhoPre_g->ptr());
    #endif
    cusphs::ComputeStepFlexStrucSymplecticPre(CaseNflexstruc,VelrhoPre_g->cptr(),Code_g->cptr(),FlexStrucRidpg,Ace_g->cptr(),dt05,Gravity,movxyg.ptr(),movzg.ptr(),Velrho_g->ptr(),NULL);
    cusph::ComputeStepPosFlexStruc(CaseNflexstruc,FlexStrucRidpg,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr(),Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
    if(!cusph::FlexStrucStepIsValid(Npb,Code_g->cptr()))
      Run_Exceptioon("Issue with FlexStruc particle update step (either naturally moved out of the domain or blew up).");
    BoundChanged=true;
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end>

  Timersg->TmStop(TMG_SuComputeStep,false);
}

//==============================================================================
/// Updates particles according to forces and dt using Symplectic-Corrector.
/// Actualizacion de particulas segun fuerzas y dt usando Symplectic-Corrector.
//==============================================================================
void JSphGpu::ComputeSymplecticCorr(double dt){
  Timersg->TmStart(TMG_SuComputeStep,false);
  const bool shift=(ShiftingMode!=SHIFT_None);
  const bool inout=(InOut!=NULL);
  const bool shiftadv=(ShiftingAdv!=NULL); //<vs_advshift>
  //-Allocate memory to compute the diplacement.
  agdouble2 movxyg("movxyg",Arrays_Gpu,true);
  agdouble  movzg("movzg",Arrays_Gpu,true);
  //-Computes displacement, velocity and density.
  const double dt05=dt*.5;
  const float3* indirvel=(InOut? InOut->GetDirVelg(): NULL);
  const byte*   boundmode=AG_CPTR(BoundMode_g); //<vs_m2dbc>
  const float4* shiftvel=AG_CPTR(ShiftVel_g); //<vs_advshift>
  cusphs::ComputeStepSymplecticCor(WithFloating,shift,shiftadv,inout,TMdbc2,Np,Npb
    ,VelrhoPre_g->cptr(),boundmode,Ar_g->cptr(),Ace_g->cptr(),ShiftPosfs_g->cptr()
    ,indirvel,AG_CPTR(NoPenShift_g),shiftvel,dt05,dt,RhopZero,RhopOutMin,RhopOutMax,Gravity // SHABA
    ,Code_g->ptr(),movxyg.ptr(),movzg.ptr(),Velrho_g->ptr()
    ,AG_PTR(PsiClean_g),AG_CPTR(PsiCleanPre_g),AG_CPTR(PsiCleanRhs_g),DivCleaning,NULL);  //<vs_divclean>

  //-Applies displacement to non-periodic fluid particles.
  cusph::ComputeStepPos2(PeriActive,WithFloating,Np,Npb
    ,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr()
    ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());

  //<vs_flexstruc_ini>
  //-Computes new position and velocity for flexible structures.
  if(CaseNflexstruc){
    Timersg->TmStart(TMG_SuFlexStruc,false);
    cusphs::ComputeStepFlexStrucSymplecticCor(CaseNflexstruc,VelrhoPre_g->cptr(),Code_g->cptr(),FlexStrucRidpg,Ace_g->cptr(),dt05,dt,Gravity,movxyg.ptr(),movzg.ptr(),Velrho_g->ptr(),NULL);
    cusph::ComputeStepPosFlexStruc(CaseNflexstruc,FlexStrucRidpg,PosxyPre_g->cptr(),PoszPre_g->cptr(),movxyg.cptr(),movzg.cptr(),Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Code_g->ptr());
    if(!cusph::FlexStrucStepIsValid(Npb,Code_g->cptr()))
      Run_Exceptioon("Issue with FlexStruc particle update step (either naturally moved out of the domain or blew up).");
    BoundChanged=true;
    Timersg->TmStop(TMG_SuFlexStruc,false);
  }
  //<vs_flexstruc_end>

  //-Free memory assigned to PRE variables in ComputeSymplecticPre().
  PosxyPre_g->Free();
  PoszPre_g->Free();
  VelrhoPre_g->Free();
  if(DivCleaning)PsiCleanPre_g->Free();  //<vs_divclean>
  Timersg->TmStop(TMG_SuComputeStep,true);
}

//==============================================================================
/// Computes the variable Dt.
/// Calcula un Dt variable.
//==============================================================================
double JSphGpu::DtVariable(bool final){
  //-dt1 depends on force per unit mass.
  const double acemax=sqrt(double(AceMax));
  const double dt1=(AceMax? (sqrt(double(KernelH)/AceMax)): DBL_MAX); 
  //-dt2 combines the Courant and the viscous time-step controls.
  const double dt2=double(KernelH)/(max(Cs0,VelMax*10.)+double(KernelH)*ViscDtMax);
  //-dt3 uses the maximum speed of sound across all structure particles.
  const double dt3=(FlexStrucDtMax? double(KernelH)/FlexStrucDtMax: DBL_MAX); //<vs_flexstruc>
  //-dt new value of time step.
  #ifdef AVAILABLE_DIVCLEAN
    //-dt4 uses the maximum local speed of sound across all fluid particles.
    const double dt4=(DivCleaning? double(KernelH)/CsPsiCleanMax: DBL_MAX); //<vs_divclean>
    double dt=CFLnumber*min(dt1,min(dt2,min(dt3,dt4)));
  #else
  double dt=CFLnumber*min(dt1,min(dt2,dt3));
  #endif
  //-Use dt value defined by FixedDt object.
  if(FixedDt)dt=FixedDt->GetDt(TimeStep,dt);
  //-Check that the dt is valid.
  if(fun::IsNAN(dt) || fun::IsInfinity(dt)){
    if(FlexStruc)Run_Exceptioon(fun::PrintStr(
    "The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f, FlexStrucDtMax=%f) is NaN or infinity at nstep=%u."
    ,dt,AceMax,VelMax,ViscDtMax,FlexStrucDtMax,Nstep));
    else         Run_Exceptioon(fun::PrintStr(
    "The computed Dt=%f (from AceMax=%f, VelMax=%f, ViscDtMax=%f) is NaN or infinity at nstep=%u."
    ,dt,AceMax,VelMax,ViscDtMax,Nstep));
  }
  //-Correct dt with minimum allowed value.
  if(dt<double(DtMin)){ 
    dt=double(DtMin); DtModif++;
    if(DtModif>=DtModifWrn){
      Log->PrintfWarning("%d DTs adjusted to DtMin (t:%g, nstep:%u)",DtModif,TimeStep,Nstep);
      DtModifWrn*=10;
    }
  }
  //-Saves information about dt.
  if(final){
    if(PartDtMin>dt)PartDtMin=dt;
    if(PartDtMax<dt)PartDtMax=dt;
    //-Saves detailed information about dt in SaveDt object.
    if(SaveDt)SaveDt->AddValues(TimeStep,dt,dt1*CFLnumber,dt2*CFLnumber
      ,dt3*CFLnumber,AceMax,ViscDtMax,FlexStrucDtMax,VelMax);
  }
  return(dt);
}

//==============================================================================
/// Computes final shifting distance for the particle position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphGpu::RunShifting(double dt){
  Timersg->TmStart(TMG_SuShifting,false);
  Shifting->RunGpu(Np-Npb,Npb,dt,Velrho_g->cptr(),ShiftPosfs_g->ptr());
  Timersg->TmStop(TMG_SuShifting,true);
}

//==============================================================================
/// Calculates predefined movement of boundary particles.
/// Calcula movimiento predefinido de boundary particles.
//==============================================================================
void JSphGpu::CalcMotion(double stepdt){
  Timersg->TmStart(TMG_SuMotion,false);
  JSph::CalcMotion(stepdt);
  Timersg->TmStop(TMG_SuMotion,false);
}

//==============================================================================
/// Processes boundary particle movement.
/// Procesa movimiento de boundary particles.
//==============================================================================
void JSphGpu::RunMotion(double stepdt){
  Timersg->TmStart(TMG_SuMotion,false);
  BoundChanged=false;
  //-Add motion from automatic wave generation.
  if(WaveGen)CalcMotionWaveGen(stepdt);
  //-Process particles motion.
  if(DsMotion->GetActiveMotion()){
    BoundChanged=true;
    const unsigned nref=DsMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      if(m.type==MOTT_Linear){//-Linear movement.
        //Log->PrintfDbg("%u] t:%g  dt:%g  mov:(%g,%g,%g)  vel:(%g,%g,%g)"
        //  ,Nstep,TimeStep,stepdt,m.linmov.x,m.linmov.y,m.linmov.z,m.linvel.x,m.linvel.y,m.linvel.z);
        cusph::MoveLinBound(PeriActive,m.count,m.idbegin-CaseNfixed,m.linmov
          ,ToTFloat3(m.linvel),RidpMotg,Posxy_g->ptr(),Posz_g->ptr()
          ,Dcell_g->ptr(),Velrho_g->ptr(),Code_g->ptr());
      }
      if(m.type==MOTT_Matrix){//-Matrix movement (for rotations).
        cusph::MoveMatBound(PeriActive,Simulate2D,m.count,m.idbegin-CaseNfixed
          ,m.matmov,stepdt,RidpMotg,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr()
          ,Velrho_g->ptr(),Code_g->ptr(),AG_PTR(BoundNor_g));
      }      
    }
  }
  //-Management of Multi-Layer Pistons.
  if(MLPistons){
    BoundChanged=true;
    if(MLPistons->GetPiston1dCount()){//-Process motion for pistons 1D.
      MLPistons->CalculateMotion1d(TimeStep+MLPistons->GetTimeMod()+stepdt);
      cusph::MovePiston1d(PeriActive!=0,CaseNmoving,0,Dp,MLPistons->GetPoszMin()
        ,MLPistons->GetPoszCount(),MLPistons->GetPistonIdGpu()
        ,MLPistons->GetMovxGpu(),MLPistons->GetVelxGpu(),RidpMotg
        ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Velrho_g->ptr()
        ,Code_g->ptr());
    }
    for(unsigned cp=0;cp<MLPistons->GetPiston2dCount();cp++){//-Process motion for pistons 2D.
      const JMLPistons::StMotionInfoPiston2D mot=MLPistons->CalculateMotion2d(cp
        ,TimeStep+MLPistons->GetTimeMod()+stepdt);
      cusph::MovePiston2d(PeriActive!=0,mot.np,mot.idbegin-CaseNfixed,Dp
        ,mot.posymin,mot.poszmin,mot.poszcount,mot.movyz,mot.velyz,RidpMotg
        ,Posxy_g->ptr(),Posz_g->ptr(),Dcell_g->ptr(),Velrho_g->ptr()
        ,Code_g->ptr());
    }
  }
  //<vs_m2dbc_ini>
  //-Copy motion velocity and compute acceleration of moving particles.
  if(MotionVel_g){
    cusph::CopyMotionVelAce(CaseNmoving,stepdt,RidpMotg,Velrho_g->cptr()
      ,MotionVel_g->ptr(),MotionAce_g->ptr());
  } //<vs_m2dbc_end>
  Timersg->TmStop(TMG_SuMotion,true);
}

//==============================================================================
/// Applies RelaxZone to selected particles.
/// Aplica RelaxZone a las particulas indicadas.
//==============================================================================
void JSphGpu::RunRelaxZone(double dt){
  Timersg->TmStart(TMG_SuMotion,false);
  agbyte  rzid("rzid",Arrays_Gpu,false);
  agfloat rzfactor("rzfactor",Arrays_Gpu,false);
  RelaxZones->SetFluidVelGpu(TimeStep,dt,Np-Npb,Npb,(const tdouble2*)Posxy_g->cptr()
    ,Posz_g->cptr(),Idp_g->cptr(),(tfloat4*)Velrho_g->ptr(),rzid.ptr(),rzfactor.ptr());
  Timersg->TmStop(TMG_SuMotion,true);
}

//==============================================================================
/// Applies Damping to indicated particles.
/// Aplica Damping a las particulas indicadas.
//==============================================================================
void JSphGpu::RunDamping(double dt){
  const typecode* codeptr=(CaseNfloat || PeriActive? Code_g->cptr(): NULL);
  Damping->ComputeDampingGpu(TimeStep,dt,Np-Npb,Npb,Posxy_g->cptr(),Posz_g->cptr()
    ,codeptr,Velrho_g->ptr());
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::SaveVtkNormalsGpu(std::string filename,int numfile,unsigned np,unsigned npb
  ,const double2* posxyg,const double* poszg,const unsigned* idpg,const float3* boundnorg)
{
  //-Allocates memory (simple way for debug method).
  const unsigned n=(UseNormalsFt? np: npb);
  tdouble3* pos=fcuda::ToHostPosd3(0,n,posxyg,poszg);
  unsigned* idp=fcuda::ToHostUint(0,n,idpg);
  tfloat3*  nor=fcuda::ToHostFloat3(0,n,boundnorg);
  //-Generates VTK file.
  SaveVtkNormals(filename,numfile,np,npb,pos,idp,nor,1.f);
  //-Frees memory.
  delete[] pos;
  delete[] idp;
  delete[] nor;
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const double2* posxyg,const double* poszg
  ,const typecode* codeg,const unsigned* idpg,const float4* velrhog,const float3* fsnormal)const
{
  const unsigned np=pfin-pini;
  //-Copy data to CPU memory.
  tdouble3* posh=fcuda::ToHostPosd3(pini,np,posxyg,poszg);
  #ifdef CODE_SIZE4
    typecode* codeh=fcuda::ToHostUint(pini,np,codeg);
  #else
    typecode* codeh=fcuda::ToHostWord(pini,np,codeg);
  #endif
  unsigned* idph=fcuda::ToHostUint(pini,np,idpg);
  tfloat4*  velrhoh=fcuda::ToHostFloat4(pini,np,velrhog);
  tfloat3*  fsnormalh=fcuda::ToHostFloat3(pini,np,fsnormal);
  //-Creates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,np,posh,codeh,idph,velrhoh,fsnormalh);
  //-Deallocates memory.
  delete[] posh;  
  delete[] codeh;
  delete[] idph;
  delete[] velrhoh;
  delete[] fsnormalh;
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const double2* posxyg,const double* poszg
  ,const typecode* codeg,const unsigned* idpg,const float4* velrhog)const
{
  const unsigned np=pfin-pini;
  //-Copy data to CPU memory.
  tdouble3* posh=fcuda::ToHostPosd3(pini,np,posxyg,poszg);
  #ifdef CODE_SIZE4
    typecode* codeh=fcuda::ToHostUint(pini,np,codeg);
  #else
    typecode* codeh=fcuda::ToHostWord(pini,np,codeg);
  #endif
  unsigned* idph=fcuda::ToHostUint(pini,np,idpg);
  tfloat4*  velrhoh=fcuda::ToHostFloat4(pini,np,velrhog);
  //-Creates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,np,posh,codeh,idph,velrhoh);
  //-Deallocates memory.
  delete[] posh;  
  delete[] codeh;
  delete[] idph;
  delete[] velrhoh;
}

//==============================================================================
/// Saves VTK file with particle data (debug).
/// Graba fichero vtk con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,unsigned cellcode,const double2* posxyg
  ,const double* poszg,const unsigned* idpg,const unsigned* dcelg
  ,const typecode* codeg,const float4* velrhog,const float4* velrhom1g
  ,const float3* aceg)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  //-Loads position.
  tfloat3* pos=new tfloat3[n];
  {
    tdouble2* pxy=new tdouble2[n];
    double* pz=new double[n];
    cudaMemcpy(pxy,posxyg+pini,sizeof(double2)*n,cudaMemcpyDeviceToHost);
    cudaMemcpy(pz,poszg+pini,sizeof(double)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)pos[p]=TFloat3(float(pxy[p].x),float(pxy[p].y),float(pz[p]));
    delete[] pxy;
    delete[] pz;
  }
  //-Loads idp.
  unsigned* idp=NULL;
  if(idpg){
    idp=new unsigned[n];
    cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  }
  //-Loads dcel.
  tuint3* dcel=NULL;
  if(dcelg){
    dcel=new tuint3[n];
    unsigned* aux=new unsigned[n];
    cudaMemcpy(aux,dcelg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++)dcel[p]=TUint3(unsigned(DCEL_Cellx(cellcode,aux[p])),unsigned(DCEL_Celly(cellcode,aux[p])),unsigned(DCEL_Cellz(cellcode,aux[p])));
    delete[] aux;
  }
  //-Loads vel and rho.
  tfloat3* vel=NULL;
  float* rho=NULL;
  if(velrhog){
    vel=new tfloat3[n];
    rho=new float[n];
    tfloat4* aux=new tfloat4[n];
    cudaMemcpy(aux,velrhog+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ vel[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rho[p]=aux[p].w; }
    delete[] aux;
  }
  //-Loads velm1 and rhom1.
  tfloat3* velm1=NULL;
  float* rhom1=NULL;
  if(velrhom1g){
    velm1=new tfloat3[n];
    rhom1=new float[n];
    tfloat4* aux=new tfloat4[n];
    cudaMemcpy(aux,velrhom1g+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ velm1[p]=TFloat3(aux[p].x,aux[p].y,aux[p].z); rhom1[p]=aux[p].w; }
    delete[] aux;
  }
  //-Loads ace.
  tfloat3* ace=NULL;
  if(aceg){
    ace=new tfloat3[n];
    cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  }
  //-Loads type.
  byte* type=NULL;
  if(codeg){
    type=new byte[n];
    typecode* aux=new typecode[n];
    cudaMemcpy(aux,codeg+pini,sizeof(typecode)*n,cudaMemcpyDeviceToHost);
    for(unsigned p=0;p<n;p++){ 
      const typecode cod=aux[p];
      byte tp=99;
      if(CODE_IsFixed(cod))tp=0;
      else if(CODE_IsMoving(cod))tp=1;
      else if(CODE_IsFloating(cod))tp=2;
      else if(CODE_IsFluid(cod))tp=3;
      if(CODE_IsNormal(cod))tp+=0;
      else if(CODE_IsPeriodic(cod))tp+=10;
      else if(CODE_IsOutIgnore(cod))tp+=20;
      else tp+=30;
      type[p]=tp;
    }
    delete[] aux;
  }
  //-Saves VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos",n,pos,true);
  if(idp)   arrays.AddArray("Idp"  ,n,idp  ,true);
  if(dcel)  arrays.AddArray("Dcel" ,n,dcel ,true);
  if(vel)   arrays.AddArray("Vel"  ,n,vel  ,true);
  if(rho)   arrays.AddArray("Rhop" ,n,rho  ,true);
  if(velm1) arrays.AddArray("Velm1",n,velm1,true);
  if(rhom1) arrays.AddArray("Rhom1",n,rhom1,true);
  if(ace)   arrays.AddArray("Ace"  ,n,ace  ,true);
  if(type)  arrays.AddArray("Typex",n,type ,true);
  JSpVtkData::Save(fun::FileNameSec(filename,numfile),arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,bool idp,bool vel,bool rho,bool code)
{
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  if(Log->GetMpiRank()>=0)filename=string("p")+fun::IntStr(Log->GetMpiRank())+"_"+filename;
  filename=DirDataOut+filename;
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned nfilter=0;
  ParticlesDataDown(n,pini,code,false,NULL,nfilter);
  tfloat3* pos=new tfloat3[n];
  for(unsigned p=0;p<n;p++)pos[p]=ToTFloat3(AuxPos_c->cptr()[p]);
  byte* type=new byte[n];
  for(unsigned p=0;p<n;p++){
    const typecode cod=Code_c->cptr()[p];
    byte tp=99;
    if(CODE_IsFixed(cod))tp=0;
    else if(CODE_IsMoving(cod))tp=1;
    else if(CODE_IsFloating(cod))tp=2;
    else if(CODE_IsFluid(cod))tp=3;
    if(CODE_IsNormal(cod))tp+=0;
    else if(CODE_IsPeriodic(cod))tp+=10;
    else if(CODE_IsOutIgnore(cod))tp+=20;
    else tp+=30;
    type[p]=tp;
  }
  //-Saves VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos",n,pos,true);
  if(idp) arrays.AddArray("Idp"  ,n,Idp_c->cptr()   ,false);
  if(type)arrays.AddArray("Typex",n,type            ,true);
  if(vel) arrays.AddArray("Vel"  ,n,AuxVel_c->cptr(),false);
  if(rho) arrays.AddArray("Rho"  ,n,AuxRho_c->cptr(),false);
#ifdef CODE_SIZE4
  if(code)arrays.AddArray("Code" ,n,(unsigned*)Code_c->cptr(),false);
#else
  if(code)arrays.AddArray("Code" ,n,(word*)    Code_c->cptr(),false);
#endif
  JSpVtkData::Save(filename,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Save VTK file with particle data (debug).
/// Graba fichero VTK con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveVtkParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const float3* posg,const byte* checkg
  ,const unsigned* idpg,const float3* velg,const float* rhog)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  tfloat3*  pos=new tfloat3[n];
  byte*     check=NULL;
  unsigned* idp=NULL;
  tfloat3*  vel=NULL;
  float*    rho=NULL;
  if(checkg)check=new byte[n];
  if(idpg)  idp=new unsigned[n];
  if(velg)  vel=new tfloat3[n];
  if(rhog)  rho=new float[n];
  //-Copies data from GPU.
  cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(checkg)cudaMemcpy(check,checkg+pini,sizeof(byte)*n,cudaMemcpyDeviceToHost);
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhog)cudaMemcpy(rho,rhog+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  //-Generates VTK file.
  DgSaveVtkParticlesCpu(filename,numfile,0,n,pos,check,idp,vel,rho);
  //-Frees memory.
  delete[] pos;
  delete[] check;
  delete[] idp;
  delete[] vel;
  delete[] rho;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,std::string head,const float3* posg
  ,const unsigned* idpg,const float3* velg,const float* rhog,const float* arg
  ,const float3* aceg,const float3* vcorrg)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned* idp=NULL;   if(idpg)idp=new unsigned[n];
  tfloat3*  pos=NULL;   if(posg)pos=new tfloat3[n];
  tfloat3*  vel=NULL;   if(velg)vel=new tfloat3[n];
  float*    rho=NULL;   if(rhog)rho=new float[n];
  float*    ar=NULL;    if(arg) ar=new float[n];
  tfloat3*  ace=NULL;   if(aceg)ace=new tfloat3[n];
  tfloat3*  vcorr=NULL; if(vcorrg)vcorr=new tfloat3[n];
  //-Copies data from GPU.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhog)cudaMemcpy(rho,rhog+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(arg)cudaMemcpy(ar,arg+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(aceg)cudaMemcpy(ace,aceg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(vcorrg)cudaMemcpy(vcorr,vcorrg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying data from GPU.");
  //-Generates CSV file.
  DgSaveCsvParticlesCpu(filename,numfile,0,n,head,pos,idp,vel,rho,ar,ace,vcorr);
  //-Frees memory.
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rho;
  delete[] ar;
  delete[] ace;
  delete[] vcorr;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticlesGpu2(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,std::string head,const float3* posg
  ,const unsigned* idpg,const float3* velg,const float* rhog
  ,const float4* pospresg,const float4* velrhog)
{
  //-Allocates memory.
  const unsigned n=pfin-pini;
  unsigned* idp=NULL;     if(idpg)idp=new unsigned[n];
  tfloat3*  pos=NULL;     if(posg)pos=new tfloat3[n];
  tfloat3*  vel=NULL;     if(velg)vel=new tfloat3[n];
  float*    rho=NULL;     if(rhog)rho=new float[n];
  tfloat4*  pospres=NULL; if(pospresg)pospres=new tfloat4[n];
  tfloat4*  velrho=NULL;  if(velrhog)velrho=new tfloat4[n];
  //-Copies data from GPU.
  if(idpg)cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*n,cudaMemcpyDeviceToHost);
  if(posg)cudaMemcpy(pos,posg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(velg)cudaMemcpy(vel,velg+pini,sizeof(float3)*n,cudaMemcpyDeviceToHost);
  if(rhog)cudaMemcpy(rho,rhog+pini,sizeof(float)*n,cudaMemcpyDeviceToHost);
  if(pospresg)cudaMemcpy(pospres,pospresg+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  if(velrhog)cudaMemcpy(velrho,velrhog+pini,sizeof(float4)*n,cudaMemcpyDeviceToHost);
  Check_CudaErroor("Failed copying data from GPU.");
  //-Generates CSV file.
  DgSaveCsvParticles2(filename,numfile,0,n,head,pos,idp,vel,rho,pospres,velrho);
  //-Frees memory.
  delete[] idp;
  delete[] pos;
  delete[] vel;
  delete[] rho;
  delete[] pospres;
  delete[] velrho;
}

//==============================================================================
/// Save CSV file with particle data (debug).
/// Graba fichero CSV con datos de las particulas (debug).
//==============================================================================
void JSphGpu::DgSaveCsvParticles2(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,std::string head,const tfloat3* pos
  ,const unsigned* idp,const tfloat3* vel,const float* rho
  ,const tfloat4* pospres,const tfloat4* velrho)
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Generates CSV file.
  ofstream pf;
  pf.open(filename.c_str());
  if(pf){
    if(!head.empty())pf << head << endl;
    pf << "Num";
    if(idp)pf << ";Idp";
    if(pos)pf << ";PosX;PosY;PosZ";
    if(vel)pf << ";VelX;VelY;VelZ";
    if(rho)pf << ";Rhop";
    if(pospres)pf << ";Px;Py;Pz;Pp";
    if(velrho) pf << ";Vx;Vy;Vz;Vr";
    pf << endl;
    const char fmt1[]="%f"; //="%24.16f";
    const char fmt3[]="%f;%f;%f"; //="%24.16f;%24.16f;%24.16f";
    for(unsigned p=pini;p<pfin;p++){
      pf << fun::UintStr(p-pini);
      if(idp)pf << ";" << fun::UintStr(idp[p]);
      if(pos)pf << ";" << fun::Float3Str(pos[p],fmt3);
      if(vel)pf << ";" << fun::Float3Str(vel[p],fmt3);
      if(rho)pf << ";" << fun::FloatStr(rho[p],fmt1);
      if(pospres)pf << ";" <<  fun::FloatStr(pospres[p].x,fmt1) << ";" << fun::FloatStr(pospres[p].y,fmt1) << ";" << fun::FloatStr(pospres[p].z,fmt1) << ";" << fun::FloatStr(pospres[p].w,fmt1);
      if(velrho)pf << ";" <<  fun::FloatStr(velrho[p].x,fmt1) << ";" << fun::FloatStr(velrho[p].y,fmt1) << ";" << fun::FloatStr(velrho[p].z,fmt1) << ";" << fun::FloatStr(velrho[p].w,fmt1);
      pf << endl;
    }
    if(pf.fail())Run_ExceptioonFile("Failed writing to file.",filename);
    pf.close();
  }
  else Run_ExceptioonFile("File could not be opened.",filename);
}

const unsigned JSphGpu::MaxNStmFloatings;


#ifdef _WITHMR //<vs_vrres_ini>
void JSphGpu::DgSaveVtkParticlesGpuMultiRes(std::string filename,int numfile,unsigned pini,unsigned pfin
 ,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg
 ,const float4 *velrhopg,const double* rcond,const float4 *shiftposg)const
{
  const unsigned np=pfin-pini;
  //-Copy data to CPU memory.
  tdouble3 *posh=fcuda::ToHostPosd3(pini,np,posxyg,poszg);
  #ifdef CODE_SIZE4
    typecode *codeh=fcuda::ToHostUint(pini,np,codeg);
  #else
    typecode *codeh=fcuda::ToHostWord(pini,np,codeg);
  #endif
  unsigned *idph=fcuda::ToHostUint(pini,np,idpg);
  tfloat4  *velrhoph=fcuda::ToHostFloat4(pini,np,velrhopg);
  tfloat4  *shiftpos=NULL;
  if(shiftposg)shiftpos=fcuda::ToHostFloat4(pini,np,shiftposg);
  double   *rcondh=NULL;
  if(rcond)rcondh=fcuda::ToHostDouble(pini,np,rcond);
  //-Creates VTK file.
  DgSaveVtkParticlesCpuMR(filename+"Fluid.vtk",numfile,0,np,posh,codeh,idph,velrhoph,NULL,NULL,NULL,shiftpos);
  DgSaveVtkParticlesCpuMRBuffer(filename+"Buffer.vtk",numfile,0,np,posh,codeh,idph,velrhoph,NULL,rcondh,shiftpos);
  DgSaveVtkParticlesCpuMRBound(filename+"Bound.vtk",numfile,0,np,posh,codeh,idph,velrhoph,NULL);
  //-Deallocates memory.
  delete[] posh;
  delete[] codeh;
  delete[] idph;
  delete[] velrhoph;
  delete[] rcondh;
}
#endif  //<vs_vrres_end>



