//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphVRes.h \brief Declares the class \ref JSphVRes.

#ifndef _JSphVResh_
#define _JSphVResh_

#include <string>
#include <vector>

#include "DualSphDef.h"
#include "JCellDivDataCpu.h"
#include "JObject.h"
#include "JSphVResZone.h"
#include "JCaseVRes.h"
#include "JDsMotion.h"
#include "JDsMotionDef.h"
#include "JMatrix4.h"

#ifdef _WITHGPU
#include <cuda_runtime_api.h>

#include "JSphGpu.h"
#endif

#define VRES_DG_SAVEVTK 0

class JXml;
class TiXmlElement;
class JLog2;
class JSphInOutZone;
class JSphCpu;
class JSphGpuSingle;
class JDsPartsInit;
class JGaugeSystem;
class JNumexLib;
class JCaseVRes;
class JDsVResDataSave;
class JDsVResDataLoad;


	
typedef struct StrDataVresCpu{
  //-Info Vres zone data.
  unsigned ntot;
  unsigned nini;

  //-Arrays with points data.
  tdouble3*	points;
  tfloat3* 	normals;
  tfloat3*  velmot;
  float*   	mass;
    //-Methods.
  StrDataVresCpu(const unsigned npoints,const unsigned nini,tdouble3* points
    ,tfloat3* normals,tfloat3* velmot,float* mass){ 
    this->ntot=npoints; this->nini=nini;  this->points=points;	
    this->normals=normals;	this->velmot=velmot; this->mass=mass;}
}StrDataVresCpu;


typedef struct StrGeomVresCpu{
  //-Info Vres zone data.
  tdouble3*    boxlimitmin;
  tdouble3*    boxlimitmax;
  tdouble3*    boxdommin;
  tdouble3*    boxdommax;
  bool*       tracking;
  bool*       inner;
  JMatrix4d*  matmov; 
  //-Methods.
  StrGeomVresCpu(tdouble3* boxlimitmin,tdouble3* boxlimitmax,tdouble3* boxdommin
    ,tdouble3* boxdommax,bool* tracking,bool* inner,JMatrix4d* matmov){ 
    	this->boxlimitmin=boxlimitmin;  this->boxlimitmax=boxlimitmax;  
      this->boxdommin=boxdommin;  this->boxdommax=boxdommax;
      this->tracking=tracking;   this->inner=inner;	this->matmov=matmov;}
}StrGeomVresCpu;
  

#ifdef _WITHGPU
  typedef struct StrDataVresGpu{
  //-Info Vres zone data.
  unsigned ntot;
  unsigned nini;
  //-Arrays with points data.
  double2*		ptposxy;
  double*     ptposz;
  float3* 	  normals;
  float3*  	  velmot;
  float*   	  mass;
  //-Methods.
  StrDataVresGpu(const unsigned npoints,const unsigned nini,double2* ptposxy,double* ptposz
    ,float3* normals,float3* velmot,float* mass){ 
    	this->ntot=npoints;     this->ptposxy=ptposxy;  this->ptposz=ptposz;	
      this->normals=normals;	this->velmot=velmot;    this->mass=mass;
      this->nini=nini;}
}StrDataVresGpu;
typedef struct StrGeomVresGpu{
  //-Info Vres zone data.
  double3*    boxlimitmin;
  double3*    boxlimitmax;
  double3*    boxdommin;
  double3*    boxdommax;
  bool*       tracking;
  bool*       inner;
  tmatrix4f*  matmov; 
  //-Methods.
  StrGeomVresGpu(double3* boxlimitmin,double3* boxlimitmax,double3* boxdommin,double3* boxdommax,
    bool* tracking,bool* inner,tmatrix4f* matmov){ 
    	this->boxlimitmin=boxlimitmin;  this->boxlimitmax=boxlimitmax;  
      this->boxdommin=boxdommin;  this->boxdommax=boxdommax;
      this->tracking=tracking;   this->inner=inner;	this->matmov=matmov;}
}StrGeomVresGpu;
#endif

class JSphVRes : protected JObject{  
private:

  const bool Cpu;
  const StCteSph CSP;     ///< Structure with main SPH constants values and configurations.
  unsigned ZoneId;
  JCaseVRes VresZone;			///< Pointer to the data structure that load vres configuration from XML file.
    
	JLog2 *Log;
  std::string AppName;
  std::string DirDataOut;
  unsigned PartBegin;
  std::string PartBeginDir;
  JDsVResDataSave* SvVResDataBi4;

  std::vector<JSphVResZone *> List;   ///<List of buffer zones.        
  unsigned ListSize;                    ///<Number of buffer zones.

  const tdouble3 MapRealPosMin;
  const tdouble3 MapRealPosMax;

  tdouble3*	 BoxLimitMin;
  tdouble3*	 BoxLimitMax;
  tdouble3*	 BoxDomMin;
  tdouble3*	 BoxDomMax;
  float*		 Width;
	bool*		   Inner;
	bool*			 Tracking;
	unsigned*	 NIni;
	unsigned*	 NPoints;
  JMatrix4d* Matmov;

	#ifdef _WITHGPU
		double3*   BoxLimitMing;
		double3*	 BoxLimitMaxg;
    double3*   BoxDomMing;
		double3*	 BoxDomMaxg;
		float*	   Widthg;
		bool*		   Innerg;
		bool*		   Trackingg;
		unsigned*	 NInig;
		unsigned*	 NPointsg;
    tmatrix4f* Matmovg;
  #endif

	unsigned	PtCount;        ///> Total number of interface points;
  tdouble3* PtPointsIni;    ///> Initial position of interface points.
  tfloat3*  PtNormalsIni;   ///> Initial definition of interface points normals.
  tdouble3* PtPoints;       ///> Position of interface points.
  tfloat3*  PtNormals;      ///> Normals of interface points.
  tfloat3*  PtVelMot;       ///> Motion Velocity of interface points.
  float*    PtMass;         ///> Mass accumulated of interface points.

  #ifdef _WITHGPU
		double2*  PtPosxyg;			///> Position of interface points.
    double*   PtPoszg;			///> Position of interface points.
		float3*   PtNormalsg;		///> Normals of interface points.
		float3*   PtVelMotg;		///> Motion Velocity of interface points.
		float*    PtMassg;			///> Mass accumulated of interface points.
  #endif
    

  void AllocateMemory(unsigned listsize);
  void FreeMemory();

  void AllocatePtMemory(unsigned ptcount);
  void FreePtMemory();

#ifdef _WITHGPU
  void AllocatePtMemoryGpu(unsigned ptcount);
  void FreePtMemoryGpu();

  void AllocateMemoryGpu(unsigned listsize);
  void FreeMemoryGpu();
#endif

void CreateZones();
void GetQuadPoints2d(tdouble3 pmin,tdouble3 pmax,tdouble3* vpt)const;
tdouble3 MovePoint(tdouble3 oldpos,const tmatrix4d& mat);

public: 
  JSphVRes(bool cpu, const StCteSph &csp,const JCaseVRes vreszone,unsigned zoneid
    ,tdouble3 maprealposmin,tdouble3 maprealposmax
    ,std::string appname,std::string dirdataout,unsigned partbegin,std::string partbegindir);
  ~JSphVRes();
  void Reset();
  void Config();

  StrDataVresCpu GetZoneFluxInfoCpu(unsigned nzone);
  StrGeomVresCpu GetGeomInfoVresCpu();
#ifdef _WITHGPU
  StrDataVresGpu GetZoneFluxInfoGpu(unsigned nzone);
  StrGeomVresGpu GetGeomInfoVresGpu();
#endif

  void SaveVResData(int part,double timestep,int nstep);
  void LoadVResData();

  unsigned GetCount()const{ return(unsigned(List.size())); }
  const JSphVResZone* GetZone(unsigned ci)const{ return(ci<ListSize? List[ci]: NULL); }

  void SaveVtkNormals(std::string filename,int numfile,unsigned np
    ,const tdouble3* pos,const tfloat3* boundnor,const tfloat3* velflux,const float* flux)const;
  void SaveNormals(std::string filename,int numfile);
  void SaveVtkDomains(std::string filename,int numfile,bool is2d)const;
  void GetRotMatrix(const JBoxDef& boxdef,JMatrix4d& mat,const tdouble3 posmin);

  void UpdateMatMov(std::vector<JMatrix4d> mat);
	void MoveBufferZone(double dt,std::vector<JMatrix4d> mat);

  unsigned CreateListCpu(unsigned npf,unsigned pini,const tdouble3 *pos
    ,const unsigned *idp,typecode *code,int *inoutpart,unsigned nzone);
  unsigned CreateListCpuInit(unsigned npf, unsigned pini,const tdouble3 *pos
    ,const unsigned *idp,typecode *code,int *inoutpart, unsigned nzone);
  unsigned CheckNormals(TpBoundary tboundary,unsigned npf,unsigned pini
    ,const tdouble2 *posxy,const double* posz,const unsigned *idp,typecode *code
    ,const tfloat3* boundnor,int *inoutpart,unsigned nzone);
  void CheckNormals(TpBoundary tboundary,unsigned np,unsigned pini
    ,const tdouble3* pos,const unsigned *idp,typecode *code
    ,const tfloat3* boundnor,unsigned vresid);
  unsigned ComputeStepCpu(unsigned bufferpartcount,int *bufferpart
		,typecode *code,const tdouble3 *pos,unsigned nzone);
  void CreateNewPart(const unsigned idnext,unsigned *dcell,typecode *code
    ,tdouble3 *pos,unsigned *idp,tfloat4 *velrhop,const JSphCpu *sphcpu
    ,unsigned np,unsigned nzone); 
  unsigned NewPartListCreate(int* newpart,unsigned idzone);
  
#ifdef _WITHGPU
  unsigned CreateListGpuInit(unsigned npf,unsigned pini,const double2 *posxyg
    ,const double *poszg,typecode *codeg,unsigned size,int *inoutpartg,unsigned nzone);
  unsigned CreateListGpu(unsigned npf,unsigned pini,const double2 *posxyg
    ,const double *poszg,typecode *codeg,unsigned size,int *inoutpartg,unsigned nzone);
  // unsigned CheckNormals(TpBoundary tboundary,unsigned npf,unsigned pini,const double2 *posxyg
  //   ,const double *poszg,const float3* boundnor,typecode *codeg
  //   ,unsigned size,int *inoutpartg,unsigned nzone);
  void ComputeStepGpu(unsigned bufferpartcount,int *bufferpart,unsigned idnext
    ,double2 *posxyg,double *poszg,unsigned *dcellg,typecode *codeg,unsigned *idpg
    ,float4 *velrhopg,byte *newizoneg, const JSphGpuSingle *gp,unsigned nzone);
  void CreateNewPartGpu(unsigned np,unsigned newnp,unsigned idnext,double2 *posxy
  ,double *posz,unsigned *dcell,typecode *code,unsigned *idp,float4 *velrhop
  ,int *newpart,unsigned nzone);
#endif
};

#endif 
