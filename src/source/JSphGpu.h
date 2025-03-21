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

/// \file JSphGpu.h \brief Declares the class \ref JSphGpu.

#ifndef _JSphGpu_
#define _JSphGpu_

#include "DualSphDef.h"
#include "JDsTimersGpu.h"
#include "JCellDivDataGpu.h"
#include "JSphGpu_cte.h"
#include "JSph.h"
#include "JArraysCpu.h"
#include "JArraysGpu.h"
#include <string>


class JDsGpuInfo;
class JDsPartsOut;
class JCellDivGpu;

//##############################################################################
//# JSphGpu
//##############################################################################
/// \brief Defines the attributes and functions used only in GPU simulations.

class JSphGpu : public JSph
{
  friend class JDebugSphGpu;

protected:
  static void RunExceptioonCudaStatic(const std::string& srcfile,int srcline
    ,const std::string& method
    ,cudaError_t cuerr,std::string msg);
  static void CheckCudaErroorStatic(const std::string& srcfile,int srcline
    ,const std::string& method
    ,std::string msg);
  void RunExceptioonCuda(const std::string& srcfile,int srcline
    ,const std::string& classname,const std::string& method
    ,cudaError_t cuerr,std::string msg)const;
  void CheckCudaErroor(const std::string& srcfile,int srcline
    ,const std::string& classname,const std::string& method
    ,std::string msg)const;

private:
  JCellDivGpu* CellDiv;

public:
  ///Structure that stores the block size to be used in each interaction kernel during GPU execution.
  typedef struct {
    unsigned forcesfluid;
    unsigned forcesbound;
    unsigned forcesdem;
  }StBlockSizes;

protected:
  StBlockSizes BlockSizes;    ///<Stores configuration of BlockSizes. | Almacena configuracion de BlockSizes.
  std::string BlockSizesStr;  ///<Stores configuration of BlockSizes in text form. | Almacena configuracion de BlockSizes en texto.

  JDsGpuInfo* GpuInfo;    ///<Main information of selected GPU.

  const TpMgDivMode DivAxis;  ///<Axis used in current division. It is used to sort particle data. MGDIV_Z is used for single GPU.

  StDivDataGpu DivData;   ///<Current data of cell division for neighborhood search on GPU.

  //-Number of particles in the domain.
  //-Numero de particulas del dominio.
  unsigned Np;        ///<Total number of particles (including duplicate periodic particles). | Numero total de particulas (incluidas las duplicadas periodicas). 
  unsigned Npb;       ///<Number of boundary particles (including periodic boundaries). | Numero de particulas contorno (incluidas las contorno periodicas). 
  unsigned NpbOk;     ///<Number of boundary particles interacting the fluid (including the periodic boundaries). | Numero de particulas contorno cerca del fluido (incluidas las contorno periodicas). 

  unsigned NpfPer;    ///<Number of periodic particles (fluid-floating). | Numero de particulas fluidas-floating periodicas. 
  unsigned NpbPer;    ///<Number of periodic boundary particles. | Numero de particulas contorno periodicas. 
  unsigned NpfPerM1;  ///<Number of fluid-floating periodic particles (previous values). | Numero de particulas fluidas-floating periodicas (valores anteriores). 
  unsigned NpbPerM1;  ///<Number of periodic boundary particles (previous values). | Numero de particulas contorno periodicas (valores anteriores).

  bool BoundChanged;  ///<Indicates if a selected boundary particle has changed since the last time step. | Indica si el contorno seleccionado a cambiado desde el ultimo divide.

  //-CPU memory allocated.
  unsigned CpuParticlesSize; ///<Number of particles for which CPU memory was allocated.
  llong MemCpuFixed;         ///<Allocated memory in AllocCpuMemoryFixed. | Mermoria reservada en AllocCpuMemoryFixed. 

  //-List of particle arrays on CPU [CpuParticlesSize=GpuParticlesSize].
  JArraysCpu* Arrays_Cpu;

  //-Execution Variables for particles [CpuParticlesSize].
  acuint*     Idp_c;    ///<Identifier of particle.
  actypecode* Code_c;   ///<Indicator of group of particles & other special markers.
  acuint*     Dcell_c;  ///<Cells inside DomCells coded with DomCellCode.
  acdouble2*  Posxy_c;  ///<Position (x,y).
  acdouble*   Posz_c;   ///<Position (z).
  acfloat4*   Velrho_c; ///<Velocity (x,y,z) and density (w).

  //-Auxiliary variables for the conversion [CpuParticlesSize].
  acdouble3*  AuxPos_c;
  acfloat3*   AuxVel_c; 
  acfloat*    AuxRho_c;


  unsigned GpuParticlesSize;  ///<Number of particles for which GPU memory was allocated. | Numero de particulas para las cuales se reservo memoria en gpu.
  llong MemGpuFixed;          ///<Allocated memory in AllocGpuMemoryFixed. | Memoria reservada en AllocGpuMemoryFixed. 

  unsigned* RidpMotg;  ///<Particle index according to Idp (only for moving and floating particles and updated after RunCellDivide) [CaseNmoving+CaseNfloat]. 

  //-List of particle arrays on GPU [GpuParticlesSize=CpuParticlesSize].
  JArraysGpu* Arrays_Gpu;
  
  //-Execution Variables for particles on GPU [GpuParticlesSize].
  aguint*     Idp_g;     ///<Identifier of particle.
  agtypecode* Code_g;    ///<Indicator of group of particles & other special markers.
  aguint*     Dcell_g;   ///<Cells inside DomCells coded with DomCellCode.
  agdouble2*  Posxy_g;
  agdouble*   Posz_g;
  agfloat4*   PosCell_g; ///<Relative position and cell coordinates for particle interaction {posx,posy,posz,cellxyz}
  agfloat4*   Velrho_g;

  aguint*     PeriParent_g; ///<Particle index to access to the parent of periodic particles (Opt).

  //-Variables for mDBC (Opt).
  agfloat3*   BoundNor_g;   ///<Normal (x,y,z) pointing from boundary particles to ghost nodes (Opt).
  agfloat3*   MotionVel_g;  ///<Velocity of a moving boundary particle (Opt).                  //<vs_m2dbc>
  agfloat3*   MotionAce_g;  ///<Acceleration of a moving boundary (Opt).                       //<vs_m2dbc>
  agbyte*     BoundMode_g;  ///<Boundary particle on off switch to multiply massp2 (Opt,Null). //<vs_m2dbc>
  agfloat3*   TangenVel_g;  ///<Velocity tangent to boundary (Opt,Null).                       //<vs_m2dbc>
    
  //-Variables for compute step VERLET (Opt).
  agfloat4*   VelrhoM1_g;  ///<Verlet: in order to keep previous values (Opt).

  //-Variables for compute step SYMPLECTIC (Opt,Null).
  agdouble2*  PosxyPre_g;  ///<Sympletic: in order to keep predictor values (Opt,Null).
  agdouble*   PoszPre_g;   ///<Sympletic: in order to keep predictor values (Opt,Null).
  agfloat4*   VelrhoPre_g; ///<Sympletic: in order to keep predictor values (Opt,Null).

  //-Variables for floating bodies (GPU memory).
  float*    FtoMasspg;   ///<Mass of the particle for each floating body [FtCount] in GPU (used in interaction forces).
  float4*   FtoDatpg;    ///<Constant data of floatings {pini_u,np_u,radius_f,massp_f} [FtCount] //__device__ int __float_as_int(float x) //__device__ float __int_as_float(int x).
  double3*  FtoCenterg;  ///<Auxiliary to store the center in floaging calculations. [FtCount].   
  float3*   FtoAceg;     ///<Auxiliary to compute sum of linear and angular acceleration of floating bodies {acelin_f3,aceang_f3} [2*FtCount].
  static const unsigned MaxNStmFloatings=18;
  cudaStream_t StmFloatings[MaxNStmFloatings];
  unsigned NStmFloatings;
  //-Variables for floating bodies (CPU memory).
  tdouble3* FtoCenterc;   ///<Auxiliary CPU Memory to swap floating center with GPU [FtCount].

  //-Variables for DEM.
  float4*   DemDatag;  ///<Data of the object {mass, (1-poisson^2)/young, kfric, restitu} in GPU [DemObjsSize].

  //<vs_flexstruc_ini>
  //-Variables for flexible structures.
  StFlexStrucData* FlexStrucDatag;  ///<Data for each individual flexible structure body [FlexStruc->GetCount()]
  unsigned* FlexStrucRidpg;         ///<Identifier to access the particles of the flexible structures [CaseNflexstruc].
  float4* PosCell0g;                ///<Relative initial position and cell coordinates {posx,posy,posz,cellxyz} [CaseNflexstruc].
  unsigned* NumPairsg;              ///<Number of initial neighbours [CaseNflexstruc].
  unsigned* PairIdxBufferg;         ///<Raw buffer to particle indices [NumPairsTot].
  unsigned** PairIdxg;              ///<List of indices to each initial neighbour [CaseNflexstruc].
  tmatrix3f* KerCorrg;              ///<Kernel correction [CaseNflexstruc].
  tmatrix3f* DefGradg;              ///<Deformation gradient tensor [CaseNflexstruc].
  float3* BoundNor0g;               ///<Initial boundary normals for mDBC [CaseNflexstruc].
  float* FlexStrucDtg;              ///<Structural speed of sound for each flexible structure particle [CaseNflexstruc].
  float FlexStrucDtMax;             ///<Maximum value of FlexStrucDt computed in Interaction_ForcesFlexStruc().
  //<vs_flexstruc_end>

  //-Variables for computing forces (Null).
  agfloat*  ViscDt_g;     ///< (Null).
  agfloat3* Ace_g;        ///<Sum of interaction acceleration (Null).
  agfloat*  Ar_g;         ///<Sum of density variation (Null). 
  agfloat*  Delta_g;      ///<Sum of Delta-SPH value when DELTA_DynamicExt (Null).
  agfloat4* ShiftPosfs_g; ///<Particle displacement and free surface detection for Shifting (Null).
  agfloat4* NoPenShift_g; ///<Particle velocity correction to prevent particle penetrating boundary (Null). //<vs_m2dbcNP>

  //<vs_advshift_ini>
  //-Variable for advanced shifting formulation.
  agfloat4* ShiftVel_g;   ///<Shifting Velocity vector for advanced shifting.
  aguint*   FSType_g;     ///<Free-surface identification.
  agfloat*  FSMinDist_g;  ///<Distance from the Free-Surface (needed for advanced shifting).
  agfloat3* FSNormal_g;   ///<Normals of Free-Surface particles (needed for advanced shifting).
  //<vs_advshift_end>

  //<vs_divclean_ini>
  agfloat*  PsiClean_g;       ///<Scalar value for divergence cleaning model.
  agfloat*  PsiCleanPre_g;    ///<Sympletic: in order to keep predictor values.
  agfloat*  PsiCleanRhs_g;    ///<Sum of Phi value for divergence cleaning.
  agfloat*  CsPsiClean_g;     ///<Local Speed of sound.
  float     CsPsiCleanMax;    ///<Max value of the local speed of sound for stability condition.
  //<vs_divclean_end>

  double VelMax;      ///<Maximum value of Vel[] sqrt(vel.x^2 + vel.y^2 + vel.z^2) computed in PreInteraction_Forces().
  double AceMax;      ///<Maximum value of Ace[] (ace.x^2 + ace.y^2 + ace.z^2) computed in Interaction_Forces().
  float ViscDtMax;    ///<Maximum value of ViscDt computed in Interaction_Forces().

  //-Variables for Laminar+SPS viscosity (Opt) & (Opt,Null).  
  agsymatrix3f* SpsTauRho2_g; ///<SPS sub-particle stress tensor divided by rho^2 (tau/rho^2) (Opt).
  agsymatrix3f* Sps2Strain_g; ///<Two times strain tensor for SPS (2S^ij) (Opt,Null).

  JDsTimersGpu* Timersg;  ///<Manages timers for GPU execution.

  void InitVars();

  void FreeCpuMemoryFixed();
  void AllocCpuMemoryFixed();
  void FreeGpuMemoryFixed();
  void AllocGpuMemoryFixed();
  void FreeCpuMemoryParticles();
  void AllocCpuMemoryParticles(unsigned np);
  void FreeGpuMemoryParticles();
  void AllocGpuMemoryParticles(unsigned np);

  void ResizeGpuMemoryParticlesData(unsigned ndatagpu,unsigned np,unsigned npmin);
  bool CheckGpuParticlesSize(unsigned requirednp)const{
    return(requirednp+PARTICLES_OVERMEMORY_MIN<=GpuParticlesSize);
  }

  llong GetAllocMemoryCpu()const;
  llong GetAllocMemoryGpu()const;
  void PrintAllocMemory(llong mcpu,llong mgpu)const;

  void GetConstantData(StCteInteraction& ctes)const;
  void ConstantDataUp();
  void ParticlesDataUp(unsigned n,const tfloat3* boundnor);
  unsigned ParticlesDataDown(unsigned n,unsigned pini,bool code
    ,bool onlynormal,const byte* filterg,unsigned& npfilterdel);
  
  int SelecDevice(int gpuid);
  void ConfigBlockSizes(bool usezone,bool useperi);

  void ConfigRunMode();
  void ConfigCellDiv(JCellDivGpu* celldiv){ CellDiv=celldiv; }
  void InitFloatingsGpu(float* ftomasspg,float4* ftodatag
    ,double3* ftocenterg,float4* demdatag)const;
  void InitRunGpu();

  void PreInteraction_Forces(TpInterStep interstep);
  void PosInteraction_Forces();
  
  void ComputeVerlet(double dt);
  void ComputeSymplecticPre(double dt);
  void ComputeSymplecticCorr(double dt);
  double DtVariable(bool final);

  void RunShifting(double dt);

  void CalcMotion(double stepdt);
  void RunMotion(double stepdt);
  void RunRelaxZone(double dt);
  void RunDamping(double dt);

  void SaveVtkNormalsGpu(std::string filename,int numfile,unsigned np,unsigned npb
    ,const double2* posxyg,const double* poszg,const unsigned* idpg,const float3* boundnorg);

public:
  JSphGpu(bool withmpi);
  ~JSphGpu();

//-Functions for debug.
//----------------------
public:
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,const double2* posxyg,const double* poszg,const typecode* codeg,const unsigned* idpg
    ,const float4* velrhog,const float3* fsnormal)const;
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,const double2* posxyg,const double* poszg,const typecode* codeg,const unsigned* idpg
    ,const float4* velrhog)const;
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,unsigned cellcode,const double2* posxyg,const double* poszg,const unsigned* idpg
    ,const unsigned* dcelg,const typecode* codeg,const float4* velrhog
    ,const float4* velrhom1g,const float3* aceg);
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,bool idp,bool vel,bool rho,bool code);
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,const float3* posg,const byte* checkg=NULL,const unsigned* idpg=NULL
    ,const float3* velg=NULL,const float* rhog=NULL);
  void DgSaveCsvParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,std::string head,const float3* posg=NULL,const unsigned* idpg=NULL
    ,const float3* velg=NULL,const float* rhog=NULL,const float* arg=NULL
    ,const float3* aceg=NULL,const float3* vcorrg=NULL);
  void DgSaveCsvParticlesGpu2(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,std::string head,const float3* posg=NULL,const unsigned* idpg=NULL
    ,const float3* velg=NULL,const float* rhog=NULL,const float4* pospres=NULL
    ,const float4* velrho=NULL);
  void DgSaveCsvParticles2(std::string filename,int numfile,unsigned pini,unsigned pfin
    ,std::string head,const tfloat3* pos=NULL,const unsigned* idp=NULL
    ,const tfloat3* vel=NULL,const float* rho=NULL,const tfloat4* pospres=NULL
    ,const tfloat4* velrho=NULL);
    
#ifdef _WITHMR //<vs_vrres_ini>
  void DgSaveVtkParticlesGpuMultiRes(std::string filename,int numfile,unsigned pini,unsigned pfin,const double2 *posxyg
  ,const double *poszg,const typecode *codeg,const unsigned *idpg,const float4 *velrhopg, const double* rcond=NULL,const float4 *shiftposg=NULL)const;
#endif        //<vs_vrres_end>
};

#endif


