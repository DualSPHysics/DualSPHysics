#ifndef JSPHBUFFER_H_
#define JSPHBUFFER_H_

#include <string>
#include <vector>

#include "DualSphDef.h"
#include "JCellDivDataCpu.h"
#include "JObject.h"
#include "JSphBufferZone.h"
#include "JCaseVRes.h"
#include "JDsMotion.h"
#include "JDsMotionDef.h"
#include "JMatrix4.h"

#ifdef _WITHGPU
#include <cuda_runtime_api.h>

#include "JSphGpu.h"
#endif

#define VRES_DG_SAVEVTK 1

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


	// typedef struct StrDataVresCpu{
  //   //-Info Vres zone data.
  //   const tfloat3*	boxlimitmin;
	// 	const tfloat3*	boxlimitmax; 
	// 	const bool*			inner;
	// 	const bool*			tracking;

  //   //-Arrays with points data.
  //   tfloat3*	points;
  //   tfloat3* 	normals;
  //   tfloat3*  velmot;
  //   float*   	mass;
  //   //-Methods.
  //   StrDataVresCpu(const tfloat3* boxlimitmin,const tfloat3*boxlimitmax,const bool* inner,
	// 		const bool* tracking,tfloat3* points,tfloat3* normals,tfloat3* velmot,float* mass){ 
  //     this->boxlimitmin=boxlimitmin;	this->boxlimitmax=boxlimitmax;
	// 		this->inner=inner;	this->tracking=tracking;
	// 		this->points=points;	this->normals=normals;	this->velmot=velmot; this->mass=mass;
  //   }
  // }StrDataVresCpu;
  

// #ifdef _WITHGPU

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
      this->nini=nini;
    }
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

    	this->boxlimitmin=boxlimitmin;  this->boxlimitmax=boxlimitmax;  this->boxdommin=boxdommin;  this->boxdommax=boxdommax;
      this->tracking=tracking;   this->inner=inner;	
      this->matmov=matmov;
    }
  }StrGeomVresGpu;

	
// #endif

class JSphBuffer : protected JObject {
  
private:
    
	JLog2 *Log;
  const std::string XmlFile;
  const std::string XmlPath;
  static const unsigned MaxZones = CODE_TYPE_FLUID_INOUTNUM;  ///< Maximum number of buffer zones.
  std::string AppName;
  std::string DirDataOut;
  unsigned PartBegin;
  std::string PartBeginDir;
  JDsVResDataSave* SvVResDataBi4;

	JCaseVRes VresZone;				///< Pointer to the data structure that load vres configuration from XML file.

  std::vector<JSphBufferZone *> List;   ///<List of buffer zones.        
  unsigned ListSize;                    ///<Number of buffer zones.
  unsigned ZoneId;

  StrGeomVresGpu* GeomInfo;
    

  std::vector<unsigned> ListNum;
  // const std::string Datafile;

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

	unsigned	PtCount;
  tdouble3* PtPointsIni;    ///> Initial position of points.
  tfloat3*  PtNormalsIni;   ///> Initial definition of points normals.

  tdouble3* PtPoints;       ///> Position of points.
  tfloat3*  PtNormals;      ///> Normals of points.
  tfloat3*  PtVelMot;       ///> Motion Velocity of points.
  float*    PtMass;         ///> Mass accumulated of points.

  #ifdef _WITHGPU
		double2*  PtPosxyg;			///> Position of points.
    double*   PtPoszg;			///> Position of points.
		float3*   PtNormalsg;			///> Normals of points.
		float3*   PtVelMotg;			///> Motion Velocity of points.
		float*    PtMassg;				///> Mass accumulated of points.
  #endif
    

    unsigned ParticlesIni=0;

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

#ifdef _WITHGPU
    double3 *BoxLimitMinInnerg = NULL;
    double3 *BoxLimitMaxInnerg = NULL;
    double3 *BoxLimitMinOuterg = NULL;
    double3 *BoxLimitMaxOuterg = NULL;
    double3 *BoxLimitMinMidg = NULL;
    double3 *BoxLimitMaxMidg = NULL;
    // bool *Innerg = NULL;
    double2 *Posxyg = NULL;
    double *Poszg = NULL;
    double3 *Normalsg = NULL;
    double *Fluxesg = NULL;
    double3 *Origing = NULL;
    double3 *BoxSizeg = NULL;
#endif

void CreateZones();
void GetQuadPoints2d(tdouble3 pmin,tdouble3 pmax,tdouble3* vpt)const;

public:


  
  const bool Cpu;
  const StCteSph CSP;  ///< Structure with main SPH constants values and configurations.
  JSphBuffer(bool cpu, const StCteSph &csp,const JCaseVRes vreszone,unsigned zoneid
      ,std::string appname,std::string dirdataout,unsigned partbegin,std::string partbegindir);
  ~JSphBuffer();
  void Reset();
  void Config();
  StrDataVresGpu GetZoneFluxInfoGpu(unsigned nzone);
  StrGeomVresGpu* GetGeomInfoVres(){return GeomInfo;};
  void SaveVResData(int part,double timestep,int nstep);
  void LoadVResData();

	


  unsigned GetCount()const{ return(unsigned(List.size())); }
  const JSphBufferZone* GetZone(unsigned ci)const{ return(ci<ListSize? List[ci]: NULL); }


  // void SaveNormals(std::string filename,int numfile);
  // void CalcMotion(const StMotionData m,double dt);
  // void CalcMotionFloating(const StFloatingData m,double dt);
  void SaveVtkNormals(std::string filename,int numfile,unsigned np
  ,const tdouble3* pos,const tfloat3* boundnor,const tfloat3* velflux,const float* flux)const;
      void SaveNormals(std::string filename,int numfile);
  void SaveVtkDomains(std::string filename,int numfile,bool is2d)const;
  void GetRotMatrix(const JBoxDef& boxdef,JMatrix4d& mat,const tdouble3 posmin);

  void UpdateMatMov(std::vector<JMatrix4d> mat);
	void MoveBufferZone(double dt,std::vector<JMatrix4d> mat);



  unsigned CreateListCpu(unsigned npf, unsigned pini, const tdouble3 *pos, const unsigned *idp, typecode *code, int *inoutpart, unsigned nzone);
  unsigned CreateListCpuInit(unsigned npf, unsigned pini, const tdouble3 *pos, const unsigned *idp, typecode *code, int *inoutpart, unsigned nzone);
  unsigned ComputeStepCpu(unsigned bufferpartcount, int *bufferpart,const unsigned idnext, unsigned *dcell, typecode *code, tdouble3 *pos, unsigned *idp,
                            tfloat4 *velrhop, const JSphCpu *sphcpu, byte *newizone, unsigned np, double dt, StDivDataCpu divdata, StCteSph CSP, unsigned i, const tdouble3 *normals, const tdouble3 *points, const tdouble4 *vel);
#ifdef _WITHGPU
  unsigned CreateListGpuInit(unsigned npf, unsigned pini, const double2 *posxyg, const double *poszg, typecode *codeg, unsigned size, int *inoutpartg, unsigned nzone);
  unsigned CreateListGpu(unsigned npf, unsigned pini, const double2 *posxyg, const double *poszg, typecode *codeg, unsigned size, int *inoutpartg, unsigned nzone);
  void ComputeStepGpu(unsigned bufferpartcount, int *bufferpart, unsigned idnext, double2 *posxyg, double *poszg, unsigned *dcellg, typecode *codeg, unsigned *idpg, float4 *velrhopg, byte *newizoneg, const JSphGpuSingle *gp, unsigned nzone);

  double2 *getPosxy() { return Posxyg; };
  double *getPosz() { return Poszg; };
  double3 *getNormals() { return Normalsg; };
  double *getFluxes() { return Fluxesg; };
  double3 *getOrigin() { return Origing; };
  double3 *getBoxSize() { return BoxSizeg; };
  bool *getInner(){ return Innerg;};
#endif





};

#endif /* JSPHBUFFER_H_ */
