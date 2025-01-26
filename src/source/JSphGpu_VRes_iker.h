/*
 * JSphGpu_Buffer_iker.h
 *
 *  Created on: Apr 25, 2022
 *      Author: francesco
 */

#ifndef _JSphGpu_VRes_iker_
#define _JSphGpu_VRes_iker_
#include "DualSphDef.h"
#include "JCellDivDataGpu.h"
#include "JSphGpu_ker.h"
#include <cuda_runtime_api.h>
#include "JSphVRes.h"
#include "JSphVResDef.h"

typedef struct StrInterParmsbg{
  //-Configuration options.
  bool simulate2d;
  TpKernel tkernel;
  StDivDataGpu divdatag;
  StCteInteraction cte;
  tdouble3 mapposmin;
  //-Input data arrays.
  const unsigned *dcell;
  const double2 *posxy;
  const double *posz;
  const float4 *poscell;
  const float4 *velrho;
  const unsigned *idp;
  const typecode *code;
  //-Other values and objects.
  cudaStream_t stm;
  StKerInfo *kerinfo;

  ///Structure constructor.
  StrInterParmsbg(
     bool simulate2d_
    ,TpKernel tkernel_
    ,const StDivDataGpu &divdatag_,const StCteInteraction cte_,const tdouble3 mapposmin_
	,const unsigned *dcell_
    ,const double2 *posxy_,const double *posz_,const float4 *poscell_
    ,const float4 *velrhop_,const unsigned *idp_,const typecode *code_
    ,cudaStream_t stm_
    ,StKerInfo *kerinfo_)
  {
    //-Configuration options.
    this->simulate2d=simulate2d_;
    this->tkernel=tkernel_;
    this->divdatag=divdatag_;
    this->cte=cte_;
    this->mapposmin=mapposmin_;
    //-Input data arrays.
    this->dcell=dcell_;
    this->posxy=posxy_; this->posz=posz_; this->poscell=poscell_;
    this->velrho=velrhop_; this->idp=idp_; this->code=code_;
    //-Other values and objects.
    this->stm=stm_;
    this->kerinfo=kerinfo_;
  }

}StInterParmsbg;

namespace cusphvres{

unsigned BufferCreateList(bool stable,unsigned n,unsigned pini
  ,const tdouble3 boxlimitmininner,const tdouble3 boxlimitmaxinner
  ,const tdouble3 boxlimitminouter,const tdouble3 boxlimitmaxouter
  ,const bool inner,const double2 *posxy,const double *posz,typecode *code
  ,unsigned *listp,tmatrix4f* mat,bool tracking,unsigned nzone);

unsigned BufferCreateListInit(bool stable,unsigned n,unsigned pini
  ,const tdouble3 boxlimitmininner,const tdouble3 boxlimitmaxinner
  ,const tdouble3 boxlimitminouter,const tdouble3 boxlimitmaxouter
  ,const tdouble3 boxlimitminmid,const tdouble3 boxlimitmaxmid
  ,const bool inner,const double2 *posxy,const double *posz
  ,typecode *code,unsigned *listp,tmatrix4f mat,bool tracking,unsigned nzone);


unsigned BufferCheckNormals(TpBoundary tboundary,bool stable,unsigned n,unsigned pini
  ,const tdouble3 boxlimitmininner,const tdouble3 boxlimitmaxinner
  ,const tdouble3 boxlimitminouter,const tdouble3 boxlimitmaxouter
  ,const bool inner,const double2 *posxy,const double *posz,const float3* boundnor
  ,typecode *code,unsigned *listp,tmatrix4f* mat,bool tracking,unsigned nzone);

void Interaction_BufferExtrap(unsigned bufferpartcount,const int *bufferpart
  ,const StInterParmsbg &t,const double2 *posxyb,const double *poszb,float4* velrhop
  ,typecode *code1,bool fastsingle,const TpVresOrder order,TpVresMethod vrmethod,float mrthreshold);

void Interaction_BufferExtrapFlux(const StInterParmsbg &t,StrDataVresGpu &vres
  ,double dp,double dt,bool fastsingle,const TpVresOrder vrorder,TpVresMethod vrmethod,float mrthreshold);

void CheckMassFlux(unsigned n,unsigned pini
  ,const StDivDataGpu& dvd,const tdouble3& mapposmin,const double2* posxy
  ,const double* posz,const typecode *code,const float4* poscell,const double2 *posxyb
  ,const double *poszb,float3 *normals,float *fluxes);

void BufferComputeStep(unsigned n,int *inoutpart,const double2 *posxy,const double *posz,typecode *code
    ,const tdouble3 boxlimitmininner,const tdouble3 boxlimitmaxinner,const tdouble3 boxlimitminouter
    ,const tdouble3 boxlimitmaxouter,const bool inner,tmatrix4f* mat,bool tracking,unsigned nzone);

void ComputeFluxesBuffer(unsigned n,double3 *normals,double *fluxes,double* rhs,double dp,double dt,bool simulate2d);

void CopySolutionBuffer(unsigned n,int *bufferpart,float4 *velrhop,double* rhs,bool simulate2d);

void CreateNewPart(unsigned newnp,unsigned pini,int *newpart
  ,unsigned np,unsigned idnext,double2 *posxy,double *posz
  ,unsigned *dcell,typecode *code,unsigned *idp,float4 *velrhop
  ,const float3 *normals,const float dp,double2 *posxyb,double *poszb
  ,unsigned nzone);

unsigned NewPartListCreate(unsigned n,unsigned pini,unsigned nmax
  ,float *fluxes,int *bufferpart,double massf);

void MoveBufferZone(unsigned pini,unsigned ntot, double2 *posxy,double *posz,float3* normals,float3* velflux,double dt,tmatrix4d mat,int zone);

void BufferShiftingGpu(unsigned np,unsigned npb,const double2 *posxy,const double *posz
  ,float4 *shiftpos,typecode *code,StrGeomVresGpu& vresgdata,cudaStream_t stm);

void ComputeFSNormals(TpKernel tkernel,bool simulate2d,unsigned bsfluid,unsigned fluidini,unsigned fluidnum
    ,StDivDataGpu& dvd,const unsigned* dcell,const double2* posxy,const double* posz
    ,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp,float4* shiftposfs
    ,unsigned* fstype,float3* fsnormal,unsigned* listp,StrGeomVresGpu& vresgdata,cudaStream_t stm);
    
void ComputeUmbrellaRegion(TpKernel tkernel,bool simulate2d,unsigned bsfluid,unsigned fluidini,unsigned fluidnum
    ,StDivDataGpu& dvd,const unsigned* dcell,const double2* posxy,const double* posz
    ,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp,float4* shiftposfs
    ,unsigned* fstype,float3* fsnormal,unsigned* listp,StrGeomVresGpu& vresgdata,cudaStream_t stm);

void PreLoopInteraction(TpKernel tkernel,bool simulate2d,bool shiftadv
    ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd,const double2* posxy,const double* posz
    ,const unsigned* dcell,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp
    ,float4* shiftvel,unsigned* fstype,float3* fsnormal,float* fsmindist,StrGeomVresGpu& vresgdata,cudaStream_t stm);

}







#endif 
