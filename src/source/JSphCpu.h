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

/// \file JSphCpu.h \brief Declares the class \ref JSphCpu.

#ifndef _JSphCpu_
#define _JSphCpu_

#include "DualSphDef.h"
#include "JDsTimersCpu.h"
#include "JCellDivDataCpu.h"
#include "JSph.h"
#include "JArraysCpu.h"
#include <string>

typedef struct{
  unsigned np,npb,npbok,npf; // npf=np-npb
  StDivDataCpu divdata;
  const unsigned *dcell;
  const tdouble3 *pos;
  const tfloat4 *velrhop;
  const unsigned *idp;
  const typecode *code;
   StCteSph csp;
}stinterparmscb;


///Structure with the parameters for particle interaction on CPU.
typedef struct{
  unsigned np,npb,npbok,npf; // npf=np-npb
  StDivDataCpu divdata;
  const unsigned* dcell;
  const tdouble3* pos;
  const tfloat4*  velrho;
  const unsigned* idp;
  const typecode* code;
  const float*    press;
  const byte*     boundmode;    //<vs_m2dbc>
  const tfloat3*  tangenvel;    //<vs_m2dbc>
  const tfloat3*  motionvel;    //<vs_m2dbc>
  const tfloat3*  dengradcorr;
  float*   ar;
  tfloat3* ace;
  float*   delta;
  TpShifting shiftmode;
  tfloat4*   shiftposfs;
  tsymatrix3f* spstaurho2;
  tsymatrix3f* sps2strain;
  unsigned*     fstype;         //<AdvancedShifting>
  tfloat4*      shiftvel;       //<AdvancedShifting>
  tmatrix3d*    lcorr;          //<AdvancedShifting>        
  float*        fstresh;        //<AdvancedShifting> 
  tfloat3*      presssym;       //<AdvancedShifting> 
  tfloat3*      pressasym;      //<AdvancedShifting>     
}stinterparmsc;

///Collects parameters for particle interaction on CPU.
inline stinterparmsc StInterparmsc(unsigned np,unsigned npb,unsigned npbok
  ,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const tfloat4* velrho,const unsigned* idp
  ,const typecode* code,const float* press
  ,const byte*    boundmode    //<vs_m2dbc>
  ,const tfloat3* tangenvel    //<vs_m2dbc>
  ,const tfloat3* motionvel    //<vs_m2dbc>
  ,const tfloat3* dengradcorr
  ,float* ar,tfloat3* ace,float* delta
  ,TpShifting shiftmode,tfloat4* shiftposfs
  ,tsymatrix3f* spstaurho2,tsymatrix3f* sps2strain
  ,unsigned* fstype,tfloat4* shiftvel,tmatrix3d* lcorr    //<AdvancedShifting
  ,float* fstresh,tfloat3* presssym,tfloat3* pressasym       //<AdvancedShifting
)
{
  stinterparmsc d={np,npb,npbok,(np-npb)
    ,divdata,dcell
    ,pos,velrho,idp
    ,code,press
    ,boundmode,tangenvel,motionvel //<vs_m2dbc>
    ,dengradcorr
    ,ar,ace,delta
    ,shiftmode,shiftposfs
    ,spstaurho2,sps2strain
    ,fstype,shiftvel,lcorr        //<AdvancedShifting
    ,fstresh,presssym,pressasym   //<AdvancedShifting      
  };
  return(d);
}


///Structure to collect interaction results.
typedef struct{
  float viscdt;
}StInterResultc;


class JDsPartsOut;
class JCellDivCpu;

//##############################################################################
//# JSphCpu
//##############################################################################
/// \brief Defines the attributes and functions to be used only in CPU simulations.

class JSphCpu : public JSph
{
  friend class JDebugSphCpu;

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

  //-CPU memory allocated.
  unsigned CpuParticlesSize;  ///<Number of particles with reserved memory on the CPU.
  llong MemCpuFixed;          ///<Memory reserved in AllocMemoryFixed. | Mermoria reservada en AllocMemoryFixed.

  unsigned* RidpMot; ///<Particle index according to Idp (only for moving and floating particles and updated after RunCellDivide) [CaseNmoving+CaseNfloat]. 

  //-List of particle arrays on CPU [CpuParticlesSize].
  JArraysCpu* Arrays_Cpu;

  //-Execution Variables for particles [CpuParticlesSize].
  acuint*     Idp_c;     ///<Identifier of particle.
  actypecode* Code_c;    ///<Indicator of group of particles & other special markers.
  acuint*     Dcell_c;   ///<Cells inside DomCells coded with DomCellCode.
  acdouble3*  Pos_c;
  acfloat4*   Velrho_c;

  //-Variables for mDBC (Opt).
  acfloat3*   BoundNor_c;   ///<Normal (x,y,z) pointing from boundary particles to ghost nodes (Opt).
  acfloat3*   MotionVel_c;  ///<Velocity of a moving boundary particle (Opt).                  //<vs_m2dbc>
  acfloat3*   MotionAce_c;  ///<Acceleration of a moving boundary (Opt).                       //<vs_m2dbc>
  acbyte*     BoundMode_c;  ///<Boundary particle on off switch to multiply massp2 (Opt,Null). //<vs_m2dbc>
  acfloat3*   TangenVel_c;  ///<Velocity tangent to boundary (Opt,Null).                       //<vs_m2dbc>
    
  //-Variables for compute step VERLET (Opt).
  acfloat4*   VelrhoM1_c;   ///<Verlet: in order to keep previous values (Opt).

  //-Variables for compute step SYMPLECTIC (Opt,Null).
  acdouble3*  PosPre_c;     ///<Sympletic: in order to keep predictor values (Opt,Null).
  acfloat4*   VelrhoPre_c;  ///<Sympletic: in order to keep predictor values (Opt,Null).

  //-Variables for computing forces (Null).
  acfloat3*   Ace_c;        ///<Sum of interaction acceleration (Null).
  acfloat*    Ar_c;         ///<Sum of density variation (Null). 
  acfloat*    Press_c;      ///<Pressure computed starting from density for interaction (Null). Press[]=fsph::ComputePress(Rho,CSP)
  acfloat*    Delta_c;      ///<Sum of Delta-SPH value when DELTA_DynamicExt (Null).
  acfloat4*   ShiftPosfs_c; ///<Particle displacement and free surface detection for Shifting (Null).

  //-Variable for advanced shifting formulation.
  acfloat4*   ShiftVel_c;       ///<Shifting Velocity vector for advanced shifting.
  acuint*     FSType_c;         ///<Free-surface identification.
  acfloat*    FSMinDist_c;      ///<Distance from the Free-Surface (needed for advanced shifting).
  acfloat3*   FSNormal_c;       ///<Normals of Free-Surface particles (needed for advanced shifting).
  acfloat*    FSTresh_c;        ///<Divergence of position needed to identify probably free-surface particles).
  acmatrix3d* LCorr_c;          ///<Correction matrix needed for non-conservative pressure formulation (only in Cpu).
  acfloat3*   PressSym_c;        ///<Array to store symmetric part of the pressure gradient;
  acfloat3*   PressAsym_c;       ///<Array to store asymmetric part of the pressure gradient;

  acuint*   PeriParent_c;     ///<Particle index to access to the parent of periodic particles (Opt). //<ShiftingAdvanced>

  double VelMax;        ///<Maximum value of Vel[] sqrt(vel.x^2 + vel.y^2 + vel.z^2) computed in PreInteraction_Forces().
  double AceMax;        ///<Maximum value of Ace[] sqrt(ace.x^2 + ace.y^2 + ace.z^2) computed in Interaction_Forces().
  float ViscDtMax;      ///<Max value of ViscDt calculated in Interaction_Forces().

  //-Variables for Laminar+SPS viscosity (Opt) & (Opt,Null).  
  acsymatrix3f* SpsTauRho2_c; ///<SPS sub-particle stress tensor divided by rho^2 (tau/rho^2) (Opt).
  acsymatrix3f* Sps2Strain_c; ///<Two times strain tensor for SPS (2S^ij) (Opt,Null).

  JDsTimersCpu* Timersc;  ///<Manages timers for CPU execution.

  void InitVars();

  void FreeCpuMemoryFixed();
  void AllocCpuMemoryFixed();
  void FreeCpuMemoryParticles();
  void AllocCpuMemoryParticles(unsigned np);

  void ResizeCpuMemoryParticlesData(unsigned ndatacpu,unsigned np,unsigned npmin);
  bool CheckCpuParticlesSize(unsigned requirednp)const{ 
    return(requirednp+PARTICLES_OVERMEMORY_MIN<=CpuParticlesSize);
  }

  llong GetAllocMemoryCpu()const;
  void PrintAllocMemory(llong mcpu)const;

  unsigned GetParticlesData(unsigned n,unsigned pini,bool onlynormal
    ,unsigned* idp,tdouble3* pos,tfloat3* vel,float* rho,typecode* code
    ,const byte* filter,unsigned& npfilterdel);
  void ConfigOmp(const JSphCfgRun* cfg);

  void ConfigRunMode();
  void ConfigCellDiv(JCellDivCpu* celldiv){ CellDiv=celldiv; }
  void InitRunCpu();

  float CalcVelMaxSeq(unsigned np,const tfloat4* velrho)const;
  float CalcVelMaxOmp(unsigned np,const tfloat4* velrho)const;

  void PreInteraction_Forces(TpInterStep interstep);
  void PosInteraction_Forces();

  template<TpKernel tker,TpFtMode ftmode> void InteractionForcesBound
    (unsigned n,unsigned pini,StDivDataCpu divdata,const unsigned* dcell
    ,const tdouble3* pos,const tfloat4* velrho,const typecode* code,const unsigned* id
    ,float& viscdt,float* ar)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity
    ,bool shift,bool mdbc2
    ,bool shiftadv,bool aleform,bool ncpress> void InteractionForcesFluid       //>AdvancedShifting>
    (unsigned n,unsigned pinit,bool boundp2,float visco
    ,StDivDataCpu divdata,const unsigned* dcell
    ,const tsymatrix3f* tau,tsymatrix3f* gradvel
    ,const tdouble3* pos,const tfloat4* velrho,const typecode* code
    ,const unsigned* idp,const float* press,const tfloat3* dengradcorr
    ,const byte* boundmode,const tfloat3* tangenvel,const tfloat3* motionvel //<vs_m2dbc>
    ,float& viscdt,float* ar,tfloat3* ace,float* delta
    ,TpShifting shiftmode,tfloat4* shiftposfs
    ,unsigned* fstype,tfloat4* shiftvel,tmatrix3d* lcorr              //<AdvancedShifting>
    ,float* fstresh,tfloat3* presssym,tfloat3* pressasym)const;       //<AdvancedShifting>

  void InteractionForcesDEM(unsigned nfloat,StDivDataCpu divdata,const unsigned* dcell
    ,const unsigned* ftridp,const StDemData* demobjs
    ,const tdouble3* pos,const tfloat4* velrho,const typecode* code
    ,const unsigned* idp,float& viscdt,tfloat3* ace)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool mdbc2
    ,bool shiftadv,bool aleform,bool ncpress>         //<AdvancedShifting>
    void Interaction_ForcesCpuT(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift,bool mdbc2>
    void Interaction_Forces_ct6(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity> 
    void Interaction_Forces_ct5(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco> 
    void Interaction_Forces_ct4(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker,TpFtMode ftmode>
    void Interaction_Forces_ct3(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker>
    void Interaction_Forces_ct2(const stinterparmsc& t,StInterResultc& res)const;
  void Interaction_Forces_ct(const stinterparmsc& t,StInterResultc& res)const;


  //------------------------------------------
  //-mDBC implementation in JSphCpu_mdbc.cpp
  //------------------------------------------
  template<TpKernel tker,bool sim2d> void InteractionMdbcCorrectionT2(unsigned n
    ,StDivDataCpu divdata,const tdouble3* pos,const typecode* code
    ,const unsigned* idp,const tfloat3* boundnor,tfloat4* velrho);
  template<TpKernel tker> void Interaction_MdbcCorrectionT(const StDivDataCpu& divdata
    ,const tdouble3* pos,const typecode* code,const unsigned* idp
    ,const tfloat3* boundnor,tfloat4* velrho);
  void Interaction_MdbcCorrection(const StDivDataCpu& divdata
    ,const tdouble3* pos,const typecode* code,const unsigned* idp
    ,const tfloat3* boundnor,tfloat4* velrho);
  //------------------------------------------

  //<vs_m2dbc_ini>
  //------------------------------------------
  //-mDBC2 implementation in JSphCpu_mdbc.cpp
  //------------------------------------------
  float Mdbc2PressClone(bool sim2d,const float rhoghost,tfloat3 bnormalp1
    ,const tfloat3 gravity,const tfloat3 motacep1,const tfloat3 dpos)const;
  float Mdbc2InfNorm3x3(tmatrix3d mat)const;
  float Mdbc2InfNorm4x4(tmatrix4d mat)const;
  tfloat3 Mdbc2TangenVel(const tfloat3& boundnor,const tfloat3& velfinal)const;
  template<TpKernel tker,bool sim2d> void InteractionMdbc2CorrectionT2
    (unsigned n,const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
    ,const unsigned* idp,const tfloat3* boundnor,const tfloat3* motionvel
    ,const tfloat3* motionace,tfloat4* velrho,byte* boundmode,tfloat3* tangenvel);
  template<TpKernel tker> void Interaction_Mdbc2CorrectionT
    (const StDivDataCpu &divdata,const tdouble3* pos,const typecode* code
    ,const unsigned* idp,const tfloat3* boundnor,const tfloat3* motionvel
    ,const tfloat3* motionace,tfloat4* velrho,byte* boundmode,tfloat3* tangenvel);
  void Interaction_Mdbc2Correction(const StDivDataCpu& divdata
    ,const tdouble3* pos,const typecode* code,const unsigned* idp
    ,const tfloat3* boundnor,const tfloat3* motionvel,const tfloat3* motionace
    ,tfloat4* velrho,byte* boundmode,tfloat3* tangenvel);
  void CopyMotionVelAce(unsigned nmoving,double dt,const unsigned* ridpmot
    ,const tfloat4* velrho,tfloat3* motionvel,tfloat3* motionace)const;
  //------------------------------------------
  //<vs_m2dbc_end>

  //<ShiftingAdvanced_ini>
  //------------------------------------------
  //-Shifting Advanced implementation in JSphCpu_preloop.cpp
  //------------------------------------------


  unsigned CountFreeSurfaceParticles(unsigned npf,unsigned pini
  ,const unsigned* fstype,unsigned* listp)const;

  template<TpKernel tker,bool sim2d> void InteractionComputeFSNormals
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat3* fsnormal,unsigned* listp)const;
  template<TpKernel tker,bool sim2d> void CallComputeFSNormalsT1
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat3* fsnormal,unsigned* listp)const;
  void CallComputeFSNormals(const StDivDataCpu& divdata
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,const tfloat4* velrho,unsigned* fstype,tfloat3* fsnormal,unsigned* listp)const;

  template<TpKernel tker,bool sim2d> void InteractionCallScanUmbrellaRegion
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,const tfloat3* fsnormal,unsigned* listp)const;
  template<TpKernel tker,bool sim2d> void CallScanUmbrellaRegionT1
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,const tfloat3* fsnormal,unsigned* listp)const;
  void CallScanUmbrellaRegion(const StDivDataCpu& divdata
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,const tfloat4* velrho,unsigned* fstype,const tfloat3* fsnormal,unsigned* listp)const;


  template<TpKernel tker,bool sim2d,bool shiftadv> void PreLoopInteraction
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat4* shiftvel,tfloat3* fsnormal,float* fsmindist)const;
  template<TpKernel tker,bool sim2d,bool shiftadv> void PreLoopInteraction_ct2
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat4* shiftvel,tfloat3* fsnormal,float* fsmindist)const;
  template<TpKernel tker,bool sim2d> void PreLoopInteraction_ct1
  (unsigned n,unsigned pinit,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat4* shiftvel,tfloat3* fsnormal,float* fsmindist)const;
  void PreLoopInteraction_ct(const StDivDataCpu& divdata
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,const tfloat4* velrho,unsigned* fstype,tfloat4* shiftvel,tfloat3* fsnormal,float* fsmindist)const;

  void ComputeShiftingVel(bool simulate2d,tfloat4* shiftvel
    ,const unsigned* fstype,const tfloat3* fsnormal,const float* fsmindist,double dt,float shiftcoef,bool ale)const;

  void ComputeFsType(unsigned n,unsigned pini,unsigned* fstype
  ,const float* fstresh,bool sim2d)const;

  //------------------------------------------
  //<ShiftingAdvanced_end>


  void ComputeSpsTau(unsigned n,unsigned pini,const tfloat4* velrho
    ,const tsymatrix3f* sps2strain,tsymatrix3f* tau_rho2)const;

  void ComputeVerletVarsFluid(bool shift,const tfloat3* indirvel
    ,const tfloat4* velrho1,const tfloat4* velrho2,const byte* boundmode
    ,double dt,double dt2,const float* ar,const tfloat3* ace,const tfloat4* shiftposfs
    ,tdouble3* pos,unsigned* cell,typecode* code,tfloat4* velrhonew)const;
  void ComputeVelrhoBound(const tfloat4* velrhoold,const byte* boundmode
    ,const float* ar,double armul,tfloat4* velrhonew)const;
  void ComputeVerlet(double dt);

  void ComputeSymplecticPre(double dt);
  void ComputeSymplecticCorr(double dt);

  double DtVariable(bool final);

  void RunShifting(double dt);

  void CalcRidp(bool periactive,unsigned np,unsigned pini,unsigned idini,unsigned idfin
    ,const typecode* code,const unsigned* idp,unsigned* ridp)const;
  void MoveLinBound(unsigned np,unsigned ini,const tdouble3& mvpos,const tfloat3& mvvel
    ,const unsigned* ridpmot,tdouble3* pos,unsigned* dcell,tfloat4* velrho,typecode* code)const;
  void MoveMatBound(unsigned np,unsigned ini,tmatrix4d m,double dt,const unsigned* ridpmot
    ,tdouble3* pos,unsigned* dcell,tfloat4* velrho,typecode* code,tfloat3* boundnor)const;
  void CalcMotion(double stepdt);
  void RunMotion(double stepdt);
  void RunRelaxZone(double dt);
  void RunDamping(double dt);

  void MovePiston1d(unsigned np,unsigned ini,double poszmin,unsigned poszcount
    ,const byte* pistonid,const double* movx,const double* velx
    ,const unsigned* ridpmot,tdouble3* pos,unsigned* dcell,tfloat4* velrho,typecode* code)const;
  void MovePiston2d(unsigned np,unsigned ini
    ,double posymin,double poszmin,unsigned poszcount,const double* movx,const double* velx
    ,const unsigned* ridpmot,tdouble3* pos,unsigned* dcell,tfloat4* velrho,typecode* code)const;

public:
  JSphCpu(bool withmpi);
  ~JSphCpu();

  void UpdatePos(tdouble3 pos0,double dx,double dy,double dz,bool outrho
    ,unsigned p,tdouble3* pos,unsigned* cell,typecode* code)const;

//-Code for InOut in JSphCpu_InOut.cpp
//--------------------------------------
protected:
  tdouble3 Interaction_PosNoPeriodic(tdouble3 posp1)const;

  template<bool sim2d,TpKernel tker> void InteractionInOutExtrap_Double
    (unsigned inoutcount,const int* inoutpart,const byte* cfgzone
    ,const tplane3f* planes,const float* width,const tfloat3* dirdata,float determlimit
    ,StDivDataCpu dvd,const unsigned* dcell,const tdouble3* pos,const typecode* code
    ,const unsigned* idp,tfloat4* velrho);
  
  template<bool sim2d,TpKernel tker> void InteractionInOutExtrap_Single
    (unsigned inoutcount,const int* inoutpart,const byte* cfgzone
    ,const tplane3f* planes,const float* width,const tfloat3* dirdata,float determlimit
    ,StDivDataCpu dvd,const unsigned* dcell,const tdouble3* pos,const typecode* code
    ,const unsigned* idp,tfloat4* velrho);
  
  template<TpKernel tker> inline void Interaction_InOutExtrapT
    (byte doublemode,unsigned inoutcount,const int* inoutpart
    ,const byte* cfgzone,const tplane3f* planes
    ,const float* width,const tfloat3* dirdata,float determlimit
    ,const unsigned* dcell,const tdouble3* pos,const typecode* code
    ,const unsigned* idp,tfloat4* velrho);

  void Interaction_InOutExtrap(byte doublemode,unsigned inoutcount,const int* inoutpart
    ,const byte* cfgzone,const tplane3f* planes
    ,const float* width,const tfloat3* dirdata,float determlimit
    ,const unsigned* dcell,const tdouble3* pos,const typecode* code
    ,const unsigned* idp,tfloat4* velrho);

  float Interaction_InOutZsurf(unsigned nptz,const tfloat3* ptzpos,float maxdist,float zbottom
    ,const StDivDataCpu& divdata,const tdouble3* pos,const typecode* code);

  //-Code for VRes in JSphCpu_Buffer.cpp
  //--------------------------------------

  template<bool sim2d,TpKernel tker,unsigned order> void InteractionBufferExtrap(unsigned bufferpartcount,const int *bufferpart,
			StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,
			const typecode *code,const unsigned *idp,const tfloat4 *velrhop,const StCteSph csp,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb);

  template<TpKernel tker> void Interaction_BufferExtrapT(unsigned bufferpartcount,
	const int *bufferpart,const stinterparmscb &t,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb,unsigned order);

  void Interaction_BufferExtrap(unsigned bufferpartcount,const int *bufferpart
     ,const stinterparmscb &t,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb,unsigned order);

  template<bool sim2d,TpKernel tker,unsigned order> void InteractionBufferExtrapFlux(const unsigned n,const int pini,
		StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,
		const typecode *code,const unsigned *idp,const tfloat4 *velrhop,const StCteSph csp,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold);

  template<TpKernel tker> void Interaction_BufferExtrapFluxT(const unsigned n,const int pini
  ,const stinterparmscb &t,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold);

  void Interaction_BufferExtrapFlux(const unsigned n,const int pini
  ,const stinterparmscb &t,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold);

};

#endif


