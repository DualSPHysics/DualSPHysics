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

/// \file JSphGpu.h \brief Declares the class \ref JSphGpu.

#ifndef _JSphGpu_
#define _JSphGpu_

#include "DualSphDef.h"
#include "JSphTimersGpu.h"
#include "JCellDivDataGpu.h"
#include "JSph.h"

#include <string>


class JDsPartsOut;
class JArraysGpu;
class JCellDivGpu;

//##############################################################################
//# JSphGpu
//##############################################################################
/// \brief Defines the attributes and functions used only in GPU simulations.

class JSphGpu : public JSph
{
  friend class JDebugSphGpu;

protected:
  static void RunExceptioonCudaStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,cudaError_t cuerr,std::string msg);
  static void CheckCudaErroorStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,std::string msg);
  void RunExceptioonCuda(const std::string &srcfile,int srcline
    ,const std::string &classname,const std::string &method
    ,cudaError_t cuerr,std::string msg)const;
  void CheckCudaErroor(const std::string &srcfile,int srcline
    ,const std::string &classname,const std::string &method
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
  StBlockSizes BlockSizes;        ///<Stores configuration of BlockSizes. | Almacena configuracion de BlockSizes.
  std::string BlockSizesStr;      ///<Stores configuration of BlockSizes in text form. | Almacena configuracion de BlockSizes en texto.

  //-Variables with information for the GPU hardware.
  //-Variables con informacion del hardware GPU.
  int GpuSelect;          ///<GPU Selection (-1:no selection). | Gpu seleccionada (-1:sin seleccion).
  std::string GpuName;    ///<Name of the selected GPU.
  size_t GpuGlobalMem;    ///<Size of global memory in bytes.
  unsigned GpuSharedMem;  ///<Size of shared memory for each block in bytes.
  unsigned GpuCompute;    ///<Compute capability: 10,11,12,20... 

  const TpMgDivMode DivAxis;  ///<Axis used in current division. It is used to sort particle data. MGDIV_Z is used for single GPU.

  StDivDataGpu DivData;   ///<Current data of cell division for neighborhood search on GPU.

  //-Number of particles in the domain.
  //-Numero de particulas del dominio.
  unsigned Np;        ///<Total number of particles (including duplicate periodic particles). | Numero total de particulas (incluidas las duplicadas periodicas). 
  unsigned Npb;       ///<Number of boundary particles (including periodic boundaries). | Numero de particulas contorno (incluidas las contorno periodicas). 
  unsigned NpbOk;     ///<Number of boundary particles interacting the fluid (including the periodic bounaries). | Numero de particulas contorno cerca del fluido (incluidas las contorno periodicas). 

  unsigned NpfPer;    ///<Number of periodic particles (fluid-floating). | Numero de particulas fluidas-floating periodicas. 
  unsigned NpbPer;    ///<Number of periodic boundary particles. | Numero de particulas contorno periodicas. 
  unsigned NpfPerM1;  ///<Number of fluid-floating periodic particles (previous values). | Numero de particulas fluidas-floating periodicas (valores anteriores). 
  unsigned NpbPerM1;  ///<Number of periodic boundary particles (previous values). | Numero de particulas contorno periodicas (valores anteriores).

  bool BoundChanged;  ///<Indicates if a selected boundary particle has changed since the last time step. | Indica si el contorno seleccionado a cambiado desde el ultimo divide.

  unsigned CpuParticlesSize; ///<Number of particles for which CPU memory was allocated. | Numero de particulas para las cuales se reservo memoria en cpu. 
  llong MemCpuFixed;         ///<Allocated memory in AllocCpuMemoryFixed. | Mermoria reservada en AllocCpuMemoryFixed. 
  llong MemCpuParticles;     ///<Allocated CPU memory for arrays with particle data. | Mermoria reservada para vectores de datos de particulas. 

  //-Variables holding particle data for the execution (size=ParticlesSize).
  //-Variables con datos de las particulas para ejecucion (size=ParticlesSize).
  unsigned *Idp;      ///<Identifier of particle | Identificador de particula.
  typecode *Code;     ///<Indicator of group of particles & other special markers. | Indica el grupo de las particulas y otras marcas especiales.
  unsigned *Dcell;    ///<Cells inside DomCells coded with DomCellCode. | Celda dentro de DomCells codificada con DomCellCode.
  tdouble2 *Posxy;
  double *Posz;
  tfloat4 *Velrhop;

  //-Auxiliary variables for the conversion (size=ParticlesSize).
  //-Variables auxiliares para conversion (size=ParticlesSize).
  tdouble3 *AuxPos;
  tfloat3 *AuxVel; 
  float *AuxRhop;

  unsigned GpuParticlesAllocs;///<Number of allocations.
  unsigned GpuParticlesSize;  ///<Number of particles for which GPU memory was allocated. | Numero de particulas para las cuales se reservo memoria en gpu.
  llong MemGpuParticles;      ///<Allocated GPU memory for arrays with particle data. | Mermoria reservada para vectores de datos de particulas.
  llong MemGpuFixed;          ///<Allocated memory in AllocGpuMemoryFixed. | Memoria reservada en AllocGpuMemoryFixed. 

  //-Particle position according to the identifier for the motion.
  //-Posicion de particula segun id para motion.
  unsigned *RidpMoveg;  ///<Only for moving boundary particles [CaseNmoving] and when CaseNmoving!=0 | Solo para boundary moving particles [CaseNmoving] y cuando CaseNmoving!=0 

  //-List of particle arrays on GPU. | Lista de arrays en GPU para particulas.
  JArraysGpu* ArraysGpu;
  //-Variables holding particle data for the execution (size=ParticlesSize).
  //-Variables con datos de las particulas para ejecucion (size=ParticlesSize).
  unsigned *Idpg;   ///<Identifier of particle | Identificador de particula.
  typecode *Codeg;  ///<Indicator of group of particles & other special markers. | Indica el grupo de las particulas y otras marcas especiales.
  unsigned *Dcellg; ///<Cells inside DomCells coded with DomCellCode. | Celda dentro de DomCells codificada con DomCellCode.
  double2 *Posxyg;
  double *Poszg;
  float4 *PosCellg; ///<Relative position and cell coordiantes for particle interaction {posx,posy,posz,cellxyz}
  float4 *Velrhopg;

  float3 *BoundNormalg;  ///<Normal (x,y,z) pointing from boundary particles to ghost nodes.
  float3 *MotionVelg;    ///<Velocity of a moving boundary particle.
    
  //-Variables for compute step: VERLET.
  float4 *VelrhopM1g;  ///<Verlet: in order to keep previous values. | Verlet: para guardar valores anteriores.

  //-Variables for compute step: SYMPLECTIC.
  double2 *PosxyPreg;  ///<Sympletic: in order to keep previous values. | Sympletic: para guardar valores en predictor.
  double *PoszPreg;
  float4 *VelrhopPreg;

  //-Variables for floating bodies.
  unsigned *FtRidpg;      ///<Identifier to access to the particles of the floating object [CaseNfloat].
  float *FtoMasspg;       ///<Mass of the particle for each floating body [FtCount] in GPU (used in interaction forces).

  float4 *FtoDatpg;        ///<Constant data of floatings {pini_u,np_u,radius_f,massp_f} [FtCount] //__device__ int __float_as_int(float x) //__device__ float __int_as_float(int x).
  float  *FtoMassg;        ///<Constant data of floatings (mass_f) [FtCount] 
  byte   *FtoConstraintsg; ///<Constant value to define motion constraints.
  float3 *FtoForcesSumg;   ///<Stores forces summation for the floating bodies {sumface_f3,sumfomegaace_f3}[FtCount]. | Almacena sumatorio de fuerzas de floatings {sumface_f3,sumfomegaace_f3} [FtCount]. 
  float3 *FtoForcesg;      ///<Stores forces for the floating bodies {face_f3,fomegaace_f3} equivalent to JSphCpu::FtoForces [FtCount]. | Almacena fuerzas de floatings {face_f3,fomegaace_f3} equivalente a JSphCpu::FtoForces [FtCount]. 
  float3 *FtoForcesResg;   ///<Stores data to update floatings {fomegares_f3,fvelres_f3} equivalent to JSphCpu::FtoForcesRes. [FtCount]. | Almacena datos para actualizar floatings {fomegares_f3,fvelres_f3} equivalente a JSphCpu::FtoForcesRes. [FtCount].
  double3 *FtoCenterResg;  ///<Stores centre to update floatings. [Ftcount]. | Almacena centro para actualizar floatings. [FtCount]. 
  float3  *FtoExtForcesg;  ///<Stores the external forces to sum of each floating body. [Ftcount].| Almacena las fuerzas externas para sumar a cada objeto flotante. [FtCount]. 

  tdouble3 *FtoAuxDouble6; ///<Memory to swap floating data with GPU. [2*FtCount]. | Memoria para intercambiar datos de floatings con GPU. [2*FtCount].
  tfloat3  *FtoAuxFloat15; ///<Memory to swap floating data with GPU. [5*FtCount]. | Memoria para intercambiar datos de floatings con GPU. [5*FtCount].

  double3 *FtoCenterg;      ///<Maintains centre of floating bodies [Ftcount].   | Mantiene centro de floating. [FtCount].   
  float3  *FtoAnglesg;      ///<Maintains rotation angles from center (angle xz, angle yz, angle xy) (units:Rad) [FtCount].   
  float3  *FtoVelAceg;      ///<Maintains velocity and acceleration (linear and angular) of floating bodies (vellin,velang,acelin,aceang)  [FtCount*4].
  //float3  *FtoVelg;         ///<Maintains velocity of floating bodies [FtCount]. | Mantiene vel de floating. [FtCount].
  //float3  *FtoOmegag;       ///<Maintains omega of floating bodies [FtCount].    | Mantiene omega de floating. [FtCount].
  float4  *FtoInertiaini8g; ///<Initial state inertia tensor in world coordinates (computed or user-given) (a11,...,a21,a22,...,a32) [Ftcount*2].
  float   *FtoInertiaini1g; ///<Initial state inertia tensor in world coordinates (computed or user-given) (a33) [Ftcount].

  bool FtObjsOutdated;      ///<FtObjs[] was not updated with new GPU values.

  //-Variables for DEM. (DEM)
  float4 *DemDatag;       ///<Data of the object {mass, (1-poisson^2)/young, kfric, restitu} in GPU [DemObjsSize].

  //-Variables for computing forces
  float *ViscDtg;
  float3 *Aceg;      ///<Accumulates acceleration of the particles. | Acumula fuerzas de interaccion.
  float *Arg; 
  float *Deltag;     ///<Accumulates adjustment of Delta-SPH with DELTA_DynamicExt. | Acumula ajuste de Delta-SPH con DELTA_DynamicExt.

  float4 *ShiftPosfsg;  ///<Particle displacement and free surface detection for Shifting.

  double VelMax;      ///<Maximum value of Vel[] sqrt(vel.x^2 + vel.y^2 + vel.z^2) computed in PreInteraction_Forces().
  double AceMax;      ///<Maximum value of Ace[] (ace.x^2 + ace.y^2 + ace.z^2) computed in Interaction_Forces().
  float ViscDtMax;    ///<Maximum value of ViscDt computed in Interaction_Forces().

  //-Variables for Laminar+SPS viscosity.  
  tsymatrix3f *SpsTaug;       ///<SPS sub-particle stress tensor.
  tsymatrix3f *SpsGradvelg;   ///<Velocity gradients.

  TimersGpu Timers;  ///<Declares an array with timers for CPU (type structure \ref StSphTimerGpu).

  void InitVars();

  void FreeCpuMemoryFixed();
  void AllocCpuMemoryFixed();
  void FreeGpuMemoryFixed();
  void AllocGpuMemoryFixed();
  void FreeCpuMemoryParticles();
  void AllocCpuMemoryParticles(unsigned np);
  void FreeGpuMemoryParticles();
  void AllocGpuMemoryParticles(unsigned np,float over);

  void ResizeGpuMemoryParticles(unsigned np);
  void ReserveBasicArraysGpu();

  bool CheckGpuParticlesSize(unsigned requirednp){ return(requirednp+PARTICLES_OVERMEMORY_MIN<=GpuParticlesSize); }

  template<class T> T* TSaveArrayGpu(unsigned np,const T *datasrc)const;
  word*        SaveArrayGpu(unsigned np,const word        *datasrc)const{ return(TSaveArrayGpu<word>       (np,datasrc)); }
  unsigned*    SaveArrayGpu(unsigned np,const unsigned    *datasrc)const{ return(TSaveArrayGpu<unsigned>   (np,datasrc)); }
  int*         SaveArrayGpu(unsigned np,const int         *datasrc)const{ return(TSaveArrayGpu<int>        (np,datasrc)); }
  float*       SaveArrayGpu(unsigned np,const float       *datasrc)const{ return(TSaveArrayGpu<float>      (np,datasrc)); }
  float3*      SaveArrayGpu(unsigned np,const float3      *datasrc)const{ return(TSaveArrayGpu<float3>     (np,datasrc)); }
  float4*      SaveArrayGpu(unsigned np,const float4      *datasrc)const{ return(TSaveArrayGpu<float4>     (np,datasrc)); }
  double*      SaveArrayGpu(unsigned np,const double      *datasrc)const{ return(TSaveArrayGpu<double>     (np,datasrc)); }
  double2*     SaveArrayGpu(unsigned np,const double2     *datasrc)const{ return(TSaveArrayGpu<double2>    (np,datasrc)); }
  tsymatrix3f* SaveArrayGpu(unsigned np,const tsymatrix3f *datasrc)const{ return(TSaveArrayGpu<tsymatrix3f>(np,datasrc)); }
  template<class T> void TRestoreArrayGpu(unsigned np,T *data,T *datanew)const;
  void RestoreArrayGpu(unsigned np,word        *data,word        *datanew)const{ TRestoreArrayGpu<word>       (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,unsigned    *data,unsigned    *datanew)const{ TRestoreArrayGpu<unsigned>   (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,int         *data,int         *datanew)const{ TRestoreArrayGpu<int>        (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,float       *data,float       *datanew)const{ TRestoreArrayGpu<float>      (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,float3      *data,float3      *datanew)const{ TRestoreArrayGpu<float3>     (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,float4      *data,float4      *datanew)const{ TRestoreArrayGpu<float4>     (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,double      *data,double      *datanew)const{ TRestoreArrayGpu<double>     (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,double2     *data,double2     *datanew)const{ TRestoreArrayGpu<double2>    (np,data,datanew); }
  void RestoreArrayGpu(unsigned np,tsymatrix3f *data,tsymatrix3f *datanew)const{ TRestoreArrayGpu<tsymatrix3f>(np,data,datanew); }

  llong GetAllocMemoryCpu()const;
  llong GetAllocMemoryGpu()const;
  void PrintAllocMemory(llong mcpu,llong mgpu)const;

  void ConstantDataUp();
  void ParticlesDataUp(unsigned n,const tfloat3 *boundnormal);
  unsigned ParticlesDataDown(unsigned n,unsigned pini,bool code,bool onlynormal);
  
  void SelecDevice(int gpuid);
  void ConfigBlockSizes(bool usezone,bool useperi);

  void ConfigRunMode(std::string preinfo);
  void ConfigCellDiv(JCellDivGpu* celldiv){ CellDiv=celldiv; }
  void InitFloating();
  void InitRunGpu();

  void PreInteractionVars_Forces(unsigned np,unsigned npb);
  void PreInteraction_Forces();
  void PosInteraction_Forces();
  
  void ComputeVerlet(double dt);
  void ComputeSymplecticPre(double dt);
  void ComputeSymplecticCorr(double dt);
  double DtVariable(bool final);

  void RunShifting(double dt);

  void CalcMotion(double stepdt);
  void RunMotion(double stepdt);
  void RunRelaxZone(double dt);
  void RunDamping(double dt,unsigned np,unsigned npb,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop);

  void SaveVtkNormalsGpu(std::string filename,int numfile,unsigned np,unsigned npb
    ,const double2 *posxyg,const double *poszg,const unsigned *idpg,const float3 *boundnormalg);

  void ShowTimers(bool onlyfile=false);
  void GetTimersInfo(std::string &hinfo,std::string &dinfo)const;
  unsigned TimerGetCount()const{ return(TmgGetCount()); }
  bool TimerIsActive(unsigned ct)const{ return(TmgIsActive(Timers,(CsTypeTimerGPU)ct)); }
  float TimerGetValue(unsigned ct)const{ return(TmgGetValue(Timers,(CsTypeTimerGPU)ct)); }
  const double* TimerGetPtrValue(unsigned ct)const{ return(TmgGetPtrValue(Timers,(CsTypeTimerGPU)ct)); }
  std::string TimerGetName(unsigned ct)const{ return(TmgGetName((CsTypeTimerGPU)ct)); }
  std::string TimerToText(unsigned ct)const{ return(JSph::TimerToText(TimerGetName(ct),TimerGetValue(ct))); }

public:
  JSphGpu(bool withmpi);
  ~JSphGpu();

//-Functions for debug.
//----------------------
public:
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg,const float4 *velrhopg)const;
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,unsigned cellcode,const double2 *posxyg,const double *poszg,const unsigned *idpg,const unsigned *dcelg,const typecode *codeg,const float4 *velrhopg,const float4 *velrhopm1g,const float3 *aceg);
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,bool idp,bool vel,bool rhop,bool code);
  void DgSaveVtkParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,const float3 *posg,const byte *checkg=NULL,const unsigned *idpg=NULL,const float3 *velg=NULL,const float *rhopg=NULL);
  void DgSaveCsvParticlesGpu(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const float3 *posg=NULL,const unsigned *idpg=NULL,const float3 *velg=NULL,const float *rhopg=NULL,const float *arg=NULL,const float3 *aceg=NULL,const float3 *vcorrg=NULL);
  void DgSaveCsvParticlesGpu2(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const float3 *posg=NULL,const unsigned *idpg=NULL,const float3 *velg=NULL,const float *rhopg=NULL,const float4 *pospres=NULL,const float4 *velrhop=NULL);
  void DgSaveCsvParticles2(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const tfloat3 *pos=NULL,const unsigned *idp=NULL,const tfloat3 *vel=NULL,const float *rhop=NULL,const tfloat4 *pospres=NULL,const tfloat4 *velrhop=NULL);
};

#endif


