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
  const tfloat3*  dengradcorr;
  float*   ar;
  tfloat3* ace;
  float*   delta;
  TpShifting shiftmode;
  tfloat4*   shiftposfs;
  tsymatrix3f* spstaurho2;
  tsymatrix3f* sps2strain;
}stinterparmsc;

///Collects parameters for particle interaction on CPU.
inline stinterparmsc StInterparmsc(unsigned np,unsigned npb,unsigned npbok
  ,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const tfloat4* velrho,const unsigned* idp
  ,const typecode* code,const float* press
  ,const tfloat3* dengradcorr
  ,float* ar,tfloat3* ace,float* delta
  ,TpShifting shiftmode,tfloat4* shiftposfs
  ,tsymatrix3f* spstaurho2,tsymatrix3f* sps2strain
)
{
  stinterparmsc d={np,npb,npbok,(np-npb)
    ,divdata,dcell
    ,pos,velrho,idp
    ,code,press
    ,dengradcorr
    ,ar,ace,delta
    ,shiftmode,shiftposfs
    ,spstaurho2,sps2strain
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
  acfloat3*   MotionVel_c;  ///<Velocity of a moving boundary particle (Opt).
    
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

  //<vs_flexstruc_ini>
  //-Variables for flexible structures.
  unsigned NumPairsTot;             ///<Total number of pairs across all flexible structure bodies.
  StFlexStrucData *FlexStrucDatac;  ///<Data for each individual flexible structure body [FlexStruc->GetCount()]
  unsigned *FlexStrucRidpc;         ///<Identifier to access to the particles of the flexible structures [CaseNflexstruc].
  tdouble3 *Pos0c;                  ///<Initial particle positions [CaseNflexstruc].
  unsigned *NumPairsc;              ///<Number of initial neighbours [CaseNflexstruc].
  unsigned *PairIdxBufferc;         ///<Raw buffer to particle indices [NumPairsTot].
  unsigned **PairIdxc;              ///<List of indices to each initial neighbour [CaseNflexstruc].
  tmatrix3f *KerCorrc;              ///<Kernel correction [CaseNflexstruc].
  tmatrix3f *DefGradc;              ///<Deformation gradient tensor [CaseNflexstruc].
  float FlexStrucDtMax;             ///<Maximum value of FlexStrucDt computed in Interaction_ForcesFlexStruc().
  //<vs_flexstruc_end>

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

  void PreInteraction_Forces();
  void PosInteraction_Forces();

  template<TpKernel tker,TpFtMode ftmode> void InteractionForcesBound
    (unsigned n,unsigned pini,StDivDataCpu divdata,const unsigned* dcell
    ,const tdouble3* pos,const tfloat4* velrho,const typecode* code,const unsigned* id
    ,float& viscdt,float* ar)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift> 
    void InteractionForcesFluid(unsigned n,unsigned pini,bool boundp2,float visco
    ,StDivDataCpu divdata,const unsigned* dcell
    ,const tsymatrix3f* tau,tsymatrix3f* gradvel
    ,const tdouble3* pos,const tfloat4* velrho,const typecode* code,const unsigned* idp
    ,const float* press,const tfloat3* dengradcorr
    ,float& viscdt,float* ar,tfloat3* ace,float* delta
    ,TpShifting shiftmode,tfloat4* shiftposfs)const;

  void InteractionForcesDEM(unsigned nfloat,StDivDataCpu divdata,const unsigned* dcell
    ,const unsigned* ftridp,const StDemData* demobjs
    ,const tdouble3* pos,const tfloat4* velrho,const typecode* code,const unsigned* idp
    ,float& viscdt,tfloat3* ace)const;

  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity,bool shift> 
    void Interaction_ForcesCpuT(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco,TpDensity tdensity> void Interaction_Forces_ct5(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker,TpFtMode ftmode,TpVisco tvisco> void Interaction_Forces_ct4(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker,TpFtMode ftmode> void Interaction_Forces_ct3(const stinterparmsc& t,StInterResultc& res)const;
  template<TpKernel tker> void Interaction_Forces_ct2(const stinterparmsc& t,StInterResultc& res)const;
  void Interaction_Forces_ct(const stinterparmsc& t,StInterResultc& res)const;

  template<TpKernel tker,bool sim2d,TpSlipMode tslip> void InteractionMdbcCorrectionT2
    (unsigned n,StDivDataCpu divdata,float determlimit,float mdbcthreshold
    ,const tdouble3* pos,const typecode* code,const unsigned* idp
    ,const tfloat3* boundnor,const tfloat3* motionvel,tfloat4* velrho);
  template<TpKernel tker> void Interaction_MdbcCorrectionT(TpSlipMode slipmode,const StDivDataCpu& divdata
    ,const tdouble3* pos,const typecode* code,const unsigned* idp
    ,const tfloat3* boundnor,const tfloat3* motionvel,tfloat4* velrho);
  void Interaction_MdbcCorrection(TpSlipMode slipmode,const StDivDataCpu& divdata
    ,const tdouble3* pos,const typecode* code,const unsigned* idp
    ,const tfloat3* boundnor,const tfloat3* motionvel,tfloat4* velrho);

  void ComputeSpsTau(unsigned n,unsigned pini,const tfloat4* velrho
    ,const tsymatrix3f* sps2strain,tsymatrix3f* tau_rho2)const;

  void ComputeVerletVarsFluid(bool shift,const tfloat3* indirvel
    ,const tfloat4* velrho1,const tfloat4* velrho2,double dt,double dt2
    ,const float* ar,const tfloat3* ace,const tfloat4* shiftposfs 
    ,tdouble3* pos,unsigned* cell,typecode* code,tfloat4* velrhonew)const;
  void ComputeVelrhoBound(const tfloat4* velrhoold,const float* ar
    ,double armul,tfloat4* velrhonew)const;
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
  void CopyMotionVel(unsigned nmoving,const unsigned* ridpmot,const tfloat4* velrho,tfloat3* motionvel)const;
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

  //<vs_flexstruc_ini>
  void SetClampCodes(unsigned np,const tdouble3 *pos,const StFlexStrucData *flexstrucdata,typecode *code)const;
  bool FlexStrucHasNormals(unsigned npb,const typecode *code,const tfloat3 *boundnormals)const;
  unsigned CountFlexStrucParts(unsigned npb,const typecode *code)const;
  void CalcFlexStrucRidp(unsigned npb,const typecode *code,unsigned *flexstrucridp)const;
  void GatherToFlexStrucArray(unsigned npfs,const unsigned *flexstrucridp,const tdouble3 *fullarray,tdouble3 *flexstrucarray)const;
  unsigned CountFlexStrucPairs(unsigned np,const tdouble3 *pos0,unsigned *numpairs)const;
  void SetFlexStrucPairs(unsigned np,const tdouble3 *pos0,unsigned **pairidx)const;
  template<TpKernel tker,bool simulate2d> void CalcFlexStrucKerCorr(unsigned np,const typecode *code,const StFlexStrucData *flexstrucdata
      ,const unsigned *flexstrucridp,const tdouble3 *pos0,const unsigned *numpairs,const unsigned *const *pairidx
      ,tmatrix3f *kercorr)const;
  template<TpKernel tker,bool simulate2d> void CalcFlexStrucKerCorrT()const;
  template<TpKernel tker> void CalcFlexStrucKerCorr_ct0()const;
  void CalcFlexStrucKerCorr()const;
  template<TpKernel tker,bool simulate2d> void ComputeDefGradFlexStruc(unsigned np,const tdouble3 *pos,const typecode *code
      ,const StFlexStrucData *flexstrucdata,const unsigned *flexstrucridp
      ,const tdouble3 *pos0,const unsigned *numpairs,const unsigned *const *pairidx,const tmatrix3f *kercorr
      ,tmatrix3f *defgrad)const;
  inline tmatrix3f ComputePK1StressFlexStruc(const tmatrix3f &defgrad,const tmatrix6f &cmat)const;
  template<TpKernel tker,bool simulate2d,bool lamsps> void InteractionForcesFlexStruc(unsigned np,float visco
      ,StDivDataCpu divdata,const unsigned *dcell
      ,const tdouble3 *pos,const tfloat4 *velrhop,const float *press,const typecode *code
      ,const StFlexStrucData *flexstrucdata,const unsigned *flexstrucridp
      ,const tdouble3 *pos0,const unsigned *numpairs,const unsigned *const *pairidx,const tmatrix3f *kercorr,const tmatrix3f *defgrad
      ,float &flexstrucdt,tfloat3 *ace)const;
  template<TpKernel tker,bool simulate2d,bool lamsps> void Interaction_ForcesFlexStrucT(float &flexstrucdtmax)const;
  template<TpKernel tker,bool simulate2d> void Interaction_ForcesFlexStruc_ct1(float &flexstrucdtmax)const;
  template<TpKernel tker> void Interaction_ForcesFlexStruc_ct0(float &flexstrucdtmax)const;
  void Interaction_ForcesFlexStruc(float &flexstrucdtmax)const;
  void ComputeSemiImplicitEulerFlexStruc(double dt,tdouble3 *pos,unsigned *dcell,typecode *code)const;
  void ComputeSymplecticPreFlexStruc(double dtm,tdouble3 *pos,unsigned *dcell,typecode *code)const;
  void ComputeSymplecticCorrFlexStruc(double dtm,double dt,tdouble3 *pos,unsigned *dcell,typecode *code)const;
  //<vs_flexstruc_end>

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

};

#endif


