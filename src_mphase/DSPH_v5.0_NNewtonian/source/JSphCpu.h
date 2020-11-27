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

/// \file JSphCpu.h \brief Declares the class \ref JSphCpu.

#ifndef _JSphCpu_
#define _JSphCpu_

#include "DualSphDef.h"
#include "JSphTimersCpu.h"
#include "JCellDivDataCpu.h"
#include "JSph.h"
#include <string>


///Structure with the parameters for particle interaction on CPU.
typedef struct{
  unsigned np,npb,npbok,npf; // npf=np-npb
  StDivDataCpu divdata;
  const unsigned *dcell;
  const tdouble3 *pos;
  const tfloat4 *velrhop;
  const unsigned *idp;
  const typecode *code;
  const float *press;
  float* ar;
  tfloat3 *ace;
  float *delta;
  TpShifting shiftmode;
  tfloat4 *shiftposfs;
  tsymatrix3f *spstau;
  tsymatrix3f *spsgradvel;
  //<vs_non-Newtonian_ini>
  float *visco_eta; 
  tsymatrix3f *d_tensor;
  float *auxnn;          
  //<vs_non-Newtonian_end>
}stinterparmsc;

///Collects parameters for particle interaction on CPU.
inline stinterparmsc StInterparmsc(unsigned np,unsigned npb,unsigned npbok
  ,StDivDataCpu divdata,const unsigned *dcell
  ,const tdouble3 *pos,const tfloat4 *velrhop,const unsigned *idp,const typecode *code
  ,const float *press
  ,float* ar,tfloat3 *ace,float *delta
  ,TpShifting shiftmode,tfloat4 *shiftposfs
  ,tsymatrix3f *spstau,tsymatrix3f *spsgradvel
  ,float* visco_eta,tsymatrix3f *d_tensor,float *auxnn  //<vs_non-Newtonian>
)
{
  stinterparmsc d={np,npb,npbok,(np-npb)
    ,divdata,dcell
    ,pos,velrhop,idp,code
    ,press
    ,ar,ace,delta
    ,shiftmode,shiftposfs
    ,spstau,spsgradvel
    ,visco_eta,d_tensor,auxnn  //<vs_non-Newtonian>
  };
  return(d);
}


///Structure to collect interaction results.
typedef struct{
  float viscdt;
  float viscetadt;  //<vs_non-Newtonian>
}StInterResultc;


class JDsPartsOut;
class JArraysCpu;
class JCellDivCpu;

//##############################################################################
//# JSphCpu
//##############################################################################
/// \brief Defines the attributes and functions to be used only in CPU simulations.

class JSphCpu : public JSph
{
private:
  JCellDivCpu* CellDiv;

protected:
  int OmpThreads;       ///<Max number of OpenMP threads in execution on CPU host (minimum 1). | Numero maximo de hilos OpenMP en ejecucion por host en CPU (minimo 1).

  StDivDataCpu DivData; ///<Current data of cell division for neighborhood search on CPU.

  //-Number of particles in domain | Numero de particulas del dominio.
  unsigned Np;        ///<Total number of particles (including periodic duplicates). | Numero total de particulas (incluidas las duplicadas periodicas).
  unsigned Npb;       ///<Total number of boundary particles (including periodic boundaries). | Numero de particulas contorno (incluidas las contorno periodicas).
  unsigned NpbOk;     ///<Total number of boundary particles near fluid (including periodic duplicates). | Numero de particulas contorno cerca del fluido (incluidas las contorno periodicas).

  unsigned NpfPer;    ///<Number of periodic floating-fluid particles. | Numero de particulas fluidas-floating periodicas.
  unsigned NpbPer;    ///<Number of periodic boundary particles. | Numero de particulas contorno periodicas.
  unsigned NpfPerM1;  ///<Number of periodic floating-fluid particles (previous values). | Numero de particulas fluidas-floating periodicas (valores anteriores).
  unsigned NpbPerM1;  ///<Number of periodic boundary particles (previous values). | Numero de particulas contorno periodicas (valores anteriores).

  bool BoundChanged;  ///<Indicates if selected boundary has changed since last call of divide. | Indica si el contorno seleccionado a cambiado desde el ultimo divide.

  unsigned CpuParticlesSize;  ///<Number of particles with reserved memory on the CPU. | Numero de particulas para las cuales se reservo memoria en cpu.
  llong MemCpuParticles;      ///<Memory reserved for particles' vectors. | Mermoria reservada para vectores de datos de particulas.
  llong MemCpuFixed;          ///<Memory reserved in AllocMemoryFixed. | Mermoria reservada en AllocMemoryFixed.

  //-Particle Position according to id. | Posicion de particula segun id.
  unsigned *RidpMove; ///<Only for moving boundary particles [CaseNmoving] and when CaseNmoving!=0 | Solo para boundary moving particles [CaseNmoving] y cuando CaseNmoving!=0 

  //-List of particle arrays on CPU. | Lista de arrays en CPU para particulas.
  JArraysCpu* ArraysCpu;

  //-Execution Variables for particles (size=ParticlesSize). | Variables con datos de las particulas para ejecucion (size=ParticlesSize).
  unsigned *Idpc;    ///<Identifier of particle | Identificador de particula.
  typecode *Codec;   ///<Indicator of group of particles & other special markers. | Indica el grupo de las particulas y otras marcas especiales.
  unsigned *Dcellc;  ///<Cells inside DomCells coded with DomCellCode. | Celda dentro de DomCells codificada con DomCellCode.
  tdouble3 *Posc;
  tfloat4 *Velrhopc;

  tfloat3 *BoundNormalc;  ///<Normal (x,y,z) pointing from boundary particles to ghost nodes.
  tfloat3 *MotionVelc;    ///<Velocity of a moving boundary particle.
    
  //-Variables for compute step: VERLET. | Vars. para compute step: VERLET.
  tfloat4 *VelrhopM1c;  ///<Verlet: in order to keep previous values. | Verlet: para guardar valores anteriores.

  //-Variables for compute step: SYMPLECTIC. | Vars. para compute step: SYMPLECTIC.
  tdouble3 *PosPrec;    ///<Sympletic: in order to keep previous values. | Sympletic: para guardar valores en predictor.
  tfloat4 *VelrhopPrec;

  //-Variables for floating bodies.
  unsigned *FtRidp;             ///<Identifier to access to the particles of the floating object [CaseNfloat].
  StFtoForces *FtoForces;       ///<Stores forces of floatings [FtCount].
  StFtoForcesRes *FtoForcesRes; ///<Stores data to update floatings [FtCount].

  //-Variables for computation of forces | Vars. para computo de fuerzas.
  tfloat3 *Acec;         ///<Sum of interaction forces | Acumula fuerzas de interaccion
  float *Arc; 
  float *Deltac;         ///<Adjusted sum with Delta-SPH with DELTA_DynamicExt | Acumula ajuste de Delta-SPH con DELTA_DynamicExt

  tfloat4 *ShiftPosfsc;    ///<Particle displacement and free surface detection for Shifting.

  double VelMax;        ///<Maximum value of Vel[] sqrt(vel.x^2 + vel.y^2 + vel.z^2) computed in PreInteraction_Forces().
  double AceMax;        ///<Maximum value of Ace[] sqrt(ace.x^2 + ace.y^2 + ace.z^2) computed in Interaction_Forces().
  float ViscDtMax;      ///<Max value of ViscDt calculated in Interaction_Forces().

  //-Variables for computing forces. | Vars. derivadas para computo de fuerzas.
  float *Pressc;       ///<Pressure computed starting from density for interaction. Press[]=fsph::ComputePress(Rhop,CSP)

  //-Variables for Laminar+SPS viscosity.  
  tsymatrix3f *SpsTauc;       ///<SPS sub-particle stress tensor.
  tsymatrix3f *SpsGradvelc;   ///<Velocity gradients.

  TimersCpu Timers;


  void InitVars();

  void FreeCpuMemoryFixed();
  void AllocCpuMemoryFixed();
  void FreeCpuMemoryParticles();
  void AllocCpuMemoryParticles(unsigned np,float over);

  void ResizeCpuMemoryParticles(unsigned np);
  void ReserveBasicArraysCpu();

  bool CheckCpuParticlesSize(unsigned requirednp){ return(requirednp+PARTICLES_OVERMEMORY_MIN<=CpuParticlesSize); }

  template<class T> T* TSaveArrayCpu(unsigned np,const T *datasrc)const;
  word*        SaveArrayCpu(unsigned np,const word        *datasrc)const{ return(TSaveArrayCpu<word>       (np,datasrc)); }
  unsigned*    SaveArrayCpu(unsigned np,const unsigned    *datasrc)const{ return(TSaveArrayCpu<unsigned>   (np,datasrc)); }
  int*         SaveArrayCpu(unsigned np,const int         *datasrc)const{ return(TSaveArrayCpu<int>        (np,datasrc)); }
  float*       SaveArrayCpu(unsigned np,const float       *datasrc)const{ return(TSaveArrayCpu<float>      (np,datasrc)); }
  tfloat3*     SaveArrayCpu(unsigned np,const tfloat3     *datasrc)const{ return(TSaveArrayCpu<tfloat3>    (np,datasrc)); }
  tfloat4*     SaveArrayCpu(unsigned np,const tfloat4     *datasrc)const{ return(TSaveArrayCpu<tfloat4>    (np,datasrc)); }
  double*      SaveArrayCpu(unsigned np,const double      *datasrc)const{ return(TSaveArrayCpu<double>     (np,datasrc)); }
  tdouble3*    SaveArrayCpu(unsigned np,const tdouble3    *datasrc)const{ return(TSaveArrayCpu<tdouble3>   (np,datasrc)); }
  tsymatrix3f* SaveArrayCpu(unsigned np,const tsymatrix3f *datasrc)const{ return(TSaveArrayCpu<tsymatrix3f>(np,datasrc)); }
  template<class T> void TRestoreArrayCpu(unsigned np,T *data,T *datanew)const;
  void RestoreArrayCpu(unsigned np,word        *data,word        *datanew)const{ TRestoreArrayCpu<word>       (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,unsigned    *data,unsigned    *datanew)const{ TRestoreArrayCpu<unsigned>   (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,int         *data,int         *datanew)const{ TRestoreArrayCpu<int>        (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,float       *data,float       *datanew)const{ TRestoreArrayCpu<float>      (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,tfloat3     *data,tfloat3     *datanew)const{ TRestoreArrayCpu<tfloat3>    (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,tfloat4     *data,tfloat4     *datanew)const{ TRestoreArrayCpu<tfloat4>    (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,double      *data,double      *datanew)const{ TRestoreArrayCpu<double>     (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,tdouble3    *data,tdouble3    *datanew)const{ TRestoreArrayCpu<tdouble3>   (np,data,datanew); }
  void RestoreArrayCpu(unsigned np,tsymatrix3f *data,tsymatrix3f *datanew)const{ TRestoreArrayCpu<tsymatrix3f>(np,data,datanew); }

  llong GetAllocMemoryCpu()const;
  void PrintAllocMemory(llong mcpu)const;

  unsigned GetParticlesData(unsigned n,unsigned pini,bool onlynormal
    ,unsigned *idp,tdouble3 *pos,tfloat3 *vel,float *rhop,typecode *code);
  void ConfigOmp(const JSphCfgRun *cfg);

  void ConfigRunMode(const JSphCfgRun *cfg,std::string preinfo="");
  void ConfigCellDiv(JCellDivCpu* celldiv){ CellDiv=celldiv; }
  void InitFloating();
  void InitRunCpu();

  float CalcVelMaxSeq(unsigned np,const tfloat4* velrhop)const;
  float CalcVelMaxOmp(unsigned np,const tfloat4* velrhop)const;

  void PreInteractionVars_Forces(unsigned np,unsigned npb);
  void PreInteraction_Forces();
  void PosInteraction_Forces();

  template<TpKernel tker,TpFtMode ftmode> void InteractionForcesBound
    (unsigned n,unsigned pini,StDivDataCpu divdata,const unsigned *dcell
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *id
    ,float &viscdt,float *ar)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift> 
    void InteractionForcesFluid(unsigned n,unsigned pini,bool boundp2,float visco
    ,StDivDataCpu divdata,const unsigned *dcell
    ,const tsymatrix3f* tau,tsymatrix3f* gradvel
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
    ,const float *press
    ,float &viscdt,float *ar,tfloat3 *ace,float *delta
    ,TpShifting shiftmode,tfloat4 *shiftposfs)const;

  void InteractionForcesDEM(unsigned nfloat,StDivDataCpu divdata,const unsigned *dcell
    ,const unsigned *ftridp,const StDemData* demobjs
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
    ,float &viscdt,tfloat3 *ace)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift> 
    void Interaction_ForcesCpuT(const stinterparmsc &t,StInterResultc &res)const;
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity> void Interaction_Forces_ct5(const stinterparmsc &t,StInterResultc &res)const;
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco> void Interaction_Forces_ct4(const stinterparmsc &t,StInterResultc &res)const;
  template<TpKernel tker,TpFtMode ftmode> void Interaction_Forces_ct3(const stinterparmsc &t,StInterResultc &res)const;
  template<TpKernel tker> void Interaction_Forces_ct2(const stinterparmsc &t,StInterResultc &res)const;
  void Interaction_Forces_ct(const stinterparmsc &t,StInterResultc &res)const;

  template<TpKernel tker,bool sim2d,TpSlipMode tslip> void InteractionMdbcCorrectionT2
    (unsigned n,StDivDataCpu divdata,float determlimit,float mdbcthreshold
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp
    ,const tfloat3 *boundnormal,const tfloat3 *motionvel,tfloat4 *velrhop);
  template<TpKernel tker> void Interaction_MdbcCorrectionT(TpSlipMode slipmode,const StDivDataCpu &divdata
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp
    ,const tfloat3 *boundnormal,const tfloat3 *motionvel,tfloat4 *velrhop);
  void Interaction_MdbcCorrection(TpSlipMode slipmode,const StDivDataCpu &divdata
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp
    ,const tfloat3 *boundnormal,const tfloat3 *motionvel,tfloat4 *velrhop);

  void ComputeSpsTau(unsigned n,unsigned pini,const tfloat4 *velrhop,const tsymatrix3f *gradvel,tsymatrix3f *tau)const;

  void ComputeVerletVarsFluid(bool shift,const tfloat3 *indirvel,const tfloat4 *velrhop1,const tfloat4 *velrhop2,double dt,double dt2,tdouble3 *pos,unsigned *cell,typecode *code,tfloat4 *velrhopnew)const;
  void ComputeVelrhopBound(const tfloat4* velrhopold,double armul,tfloat4* velrhopnew)const;

  void ComputeVerlet(double dt);
  void ComputeSymplecticPre(double dt);
  void ComputeSymplecticCorr(double dt);
  double DtVariable(bool final);

  void RunShifting(double dt);

  void CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini,unsigned idfin
    ,const typecode *code,const unsigned *idp,unsigned *ridp)const;
  void MoveLinBound(unsigned np,unsigned ini,const tdouble3 &mvpos,const tfloat3 &mvvel
    ,const unsigned *ridp,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code)const;
  void MoveMatBound(unsigned np,unsigned ini,tmatrix4d m,double dt,const unsigned *ridpmv
    ,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code,tfloat3 *boundnormal)const;
  void CopyMotionVel(unsigned nmoving,const unsigned *ridp,const tfloat4 *velrhop,tfloat3 *motionvel)const;
  void CalcMotion(double stepdt);
  void RunMotion(double stepdt);
  void RunRelaxZone(double dt);
  void RunDamping(double dt,unsigned np,unsigned npb,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const;

  void MovePiston1d(unsigned np,unsigned ini,double poszmin,unsigned poszcount
    ,const byte *pistonid,const double* movx,const double* velx
    ,const unsigned *ridpmv,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code)const;
  void MovePiston2d(unsigned np,unsigned ini
    ,double posymin,double poszmin,unsigned poszcount,const double* movx,const double* velx
    ,const unsigned *ridpmv,tdouble3 *pos,unsigned *dcell,tfloat4 *velrhop,typecode *code)const;

  void ShowTimers(bool onlyfile=false);
  void GetTimersInfo(std::string &hinfo,std::string &dinfo)const;
  unsigned TimerGetCount()const{ return(TmcGetCount()); }
  bool TimerIsActive(unsigned ct)const{ return(TmcIsActive(Timers,(CsTypeTimerCPU)ct)); }
  float TimerGetValue(unsigned ct)const{ return(TmcGetValue(Timers,(CsTypeTimerCPU)ct)); }
  const double* TimerGetPtrValue(unsigned ct)const{ return(TmcGetPtrValue(Timers,(CsTypeTimerCPU)ct)); }
  std::string TimerGetName(unsigned ct)const{ return(TmcGetName((CsTypeTimerCPU)ct)); }
  std::string TimerToText(unsigned ct)const{ return(JSph::TimerToText(TimerGetName(ct),TimerGetValue(ct))); }

public:
  JSphCpu(bool withmpi);
  ~JSphCpu();

  void UpdatePos(tdouble3 pos0,double dx,double dy,double dz,bool outrhop,unsigned p,tdouble3 *pos,unsigned *cell,typecode *code)const;

//-Code for InOut in JSphCpu_InOut.cpp
//--------------------------------------
protected:
  tdouble3 Interaction_PosNoPeriodic(tdouble3 posp1)const;

  template<bool sim2d,TpKernel tker> void InteractionInOutExtrap_Double
    (unsigned inoutcount,const int *inoutpart,const byte *cfgzone
    ,const tplane3f *planes,const float* width,const tfloat3 *dirdata,float determlimit
    ,StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,const typecode *code
    ,const unsigned *idp,tfloat4 *velrhop);
  
  template<bool sim2d,TpKernel tker> void InteractionInOutExtrap_Single
    (unsigned inoutcount,const int *inoutpart,const byte *cfgzone
    ,const tplane3f *planes,const float* width,const tfloat3 *dirdata,float determlimit
    ,StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,const typecode *code
    ,const unsigned *idp,tfloat4 *velrhop);
  
  template<TpKernel tker> inline void Interaction_InOutExtrapT
    (byte doublemode,unsigned inoutcount,const int *inoutpart
    ,const byte *cfgzone,const tplane3f *planes
    ,const float* width,const tfloat3 *dirdata,float determlimit
    ,const unsigned *dcell,const tdouble3 *pos,const typecode *code
    ,const unsigned *idp,tfloat4 *velrhop);

  void Interaction_InOutExtrap(byte doublemode,unsigned inoutcount,const int *inoutpart
    ,const byte *cfgzone,const tplane3f *planes
    ,const float* width,const tfloat3 *dirdata,float determlimit
    ,const unsigned *dcell,const tdouble3 *pos,const typecode *code
    ,const unsigned *idp,tfloat4 *velrhop);

  float Interaction_InOutZsurf(unsigned nptz,const tfloat3 *ptzpos,float maxdist,float zbottom
    ,const StDivDataCpu &divdata,const tdouble3 *pos,const typecode *code);


  template<bool sim2d,TpKernel tker> void InteractionBoundCorr_Double
    (unsigned npb,typecode boundcode,tplane3f plane,tfloat3 direction,float determlimit
    ,const StDivDataCpu &dvd,const tdouble3 *pos,const typecode *code
    ,const unsigned *idp,tfloat4 *velrhop);

  template<bool sim2d,TpKernel tker> void InteractionBoundCorr_Single
    (unsigned npb,typecode boundcode,tplane3f plane,tfloat3 direction,float determlimit
    ,const StDivDataCpu &dvd,const tdouble3 *pos,const typecode *code
    ,const unsigned *idp,tfloat4 *velrhop);

  template<TpKernel tker> inline void Interaction_BoundCorrT
    (byte doublemode,typecode boundcode
    ,tplane3f plane,tfloat3 direction,float determlimit
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp
    ,tfloat4 *velrhop);

  void Interaction_BoundCorr(byte doublemode,typecode boundcode
    ,tplane3f plane,tfloat3 direction,float determlimit
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp
    ,tfloat4 *velrhop);


//<vs_non-Newtonian_ini>
// NN-MultiPhase: Variables and methods.  
//------------------------------------------
protected:
  float ViscEtaDtMax;      ///<Max value of ViscDt calculated in Interaction_Forces(). 

  //-Extra variables for viscosity.  
  float *Visco_etac;       ///<Effective viscosity.  
  tsymatrix3f *D_tensorc;  ///<Deformation tensor. 
  float *AuxNN;            ///<Auxilary. 


  unsigned GetParticlesData(unsigned n,unsigned pini,bool onlynormal
    ,unsigned *idp,tdouble3 *pos,tfloat3 *vel,float *rhop,typecode *code,float *aux_n);

  //<vs_non-Newtonian> pressure
  void ComputePress_NN(unsigned np,unsigned npb);
  //functions for tensors
  void GetVelocityGradients_FDA(float rr2,float drx,float dry,float drz
    ,float dvx,float dvy,float dvz,tmatrix3f &dvelp1,float &div_vel)const;
  void GetStrainRateTensor(const tmatrix3f &dvelp1,float div_vel
    ,float &I_D,float &II_D,float &J1_D,float &J2_D,float &div_D_tensor
    ,float &D_tensor_magn,tmatrix3f &D_tensor)const;
  void GetEta_Effective(const typecode ppx,float tau_yield,float D_tensor_magn
    ,float visco,float m_NN,float n_NN,float &visco_etap1)const;
  void GetStressTensor(const tmatrix3f &D_tensor,float visco_etap1
    ,float &I_t,float &II_t,float &J1_t,float &J2_t,float &tau_tensor_magn,tmatrix3f &tau_tensor)const;

  //SPH Symetric
  void GetVelocityGradients_SPH_tsym(float massp2,const tfloat4 &velrhop2
    ,float dvx,float dvy,float dvz,float frx,float fry,float frz
    ,tsymatrix3f &gradvelp1)const;
  void GetStrainRateTensor_tsym(const tsymatrix3f &dvelp1
    ,float &I_D,float &II_D,float &J1_D,float &J2_D,float &div_D_tensor
    ,float &D_tensor_magn,tsymatrix3f &D_tensor)const;
  void GetStressTensor_sym(const tsymatrix3f &D_tensorp1,float visco_etap1
    ,float &I_t,float &II_t,float &J1_t,float &J2_t,float &tau_tensor_magn,tsymatrix3f &tau_tensorp1)const;

  //BCs
  template<TpKernel tker,TpFtMode ftmode> void InteractionForcesBound_NN_FDA
  (unsigned n,unsigned pini,StDivDataCpu divdata,const unsigned *dcell
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *id
    ,float &viscdt,float *ar)const;
  template<TpKernel tker,TpFtMode ftmode> void InteractionForcesBound_NN_SPH
  (unsigned n,unsigned pini,StDivDataCpu divdata,const unsigned *dcell
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *id
    ,float &viscdt,float *ar)const;

  //SPH
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
  void InteractionForcesFluid_NN_SPH_ConsEq(unsigned n,unsigned pini,bool boundp2,float visco
    ,StDivDataCpu divdata,const unsigned *dcell
    ,float *visco_eta,const tsymatrix3f* tau,float *auxnn
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
    ,tfloat3 *ace)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
  void InteractionForcesFluid_NN_SPH_Morris
  (unsigned n,unsigned pini,bool boundp2,float visco
    ,StDivDataCpu divdata,const unsigned *dcell
    ,float *visco_eta,const tsymatrix3f* tau,tsymatrix3f* gradvel,float *auxnn
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
    ,tfloat3 *ace)const;

  template<TpFtMode ftmode,TpVisco tvisco> void InteractionForcesFluid_NN_SPH_Visco_Stress_tensor
  (unsigned n,unsigned pinit,float visco,float *visco_eta
    ,tsymatrix3f* tau,const tsymatrix3f* D_tensor,float *auxnn
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp)const;

  template<TpFtMode ftmode,TpVisco tvisco> void InteractionForcesFluid_NN_SPH_Visco_eta
  (unsigned n,unsigned pinit,float visco,float *visco_eta,const tfloat4 *velrhop
    ,const tsymatrix3f* gradvel,tsymatrix3f* D_tensor,float *auxnn,float &viscetadt
    ,const tdouble3 *pos,const typecode *code,const unsigned *idp)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
  void InteractionForcesFluid_NN_SPH_PressGrad
  (unsigned n,unsigned pini,bool boundp2,float visco
    ,StDivDataCpu divdata,const unsigned *dcell
    ,tsymatrix3f* gradvel
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
    ,const float *press
    ,float &viscdt,float *ar,tfloat3 *ace,float *delta
    ,TpShifting shiftmode,tfloat4 *shiftposfs)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
  void InteractionForcesFluid_NN_FDA_All(unsigned n,unsigned pini,bool boundp2,float visco,float *visco_eta
    ,StDivDataCpu divdata,const unsigned *dcell
    ,const tsymatrix3f* tau,tsymatrix3f* gradvel
    ,const tdouble3 *pos,const tfloat4 *velrhop,const typecode *code,const unsigned *idp
    ,const float *press
    ,float &viscdt,float &viscetadt,float *ar,tfloat3 *ace,float *delta
    ,TpShifting shiftmode,tfloat4 *shiftposfs)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
  void Interaction_ForcesCpuT_NN_FDA(const stinterparmsc &t,StInterResultc &res)const;

  //FDA or SPH
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift>
  void Interaction_ForcesCpuT_NN_SPH(const stinterparmsc &t,StInterResultc &res)const;
  //<vs_non-Newtonian_end>
  
};

#endif


