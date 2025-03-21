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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para medir magnitudes fisicas durante la simulacion. (12-02-2018)
//:# - Se escriben las unidades en las cabeceras de los ficheros CSV. (26-04-2018)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:# - Cambio de nombre de fichero J.GaugeItem a J.DsGaugeItem. (28-06-2020)
//:# - Permite carga de ficheros mbi4 en formato multidata. (24-10-2021)
//:#############################################################################

/// \file JDsGaugeItem.h \brief Declares the class \ref JGaugeItem.

#ifndef _JDsGaugeItem_
#define _JDsGaugeItem_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#include "JSaveCsv2.h"
#include "JArraysCpu.h"
#include "JCellDivDataCpu.h"
#include "JMeshDataDef.h" //<vs_meeshdat>

#ifdef _WITHGPU
#include "JArraysGpu.h"
#include "JCellDivDataGpu.h"
#include <cuda_runtime_api.h>
#endif

class JLog2;
namespace jmsh{          //<vs_meeshdat>
  class JMeshData;       //<vs_meeshdat>
  class JMeshTDatasSave; //<vs_meeshdat>
}                        //<vs_meeshdat>

//##############################################################################
//# JGaugeItem
//##############################################################################
/// \brief Defines the common part of different gauge types.

class JGaugeItem : protected JObject
{
protected:
  static const int MAXGPUS=1;

protected:
 #ifdef _WITHGPU
  void RunExceptioonCuda(const std::string& srcfile,int srcline
    ,const std::string& classname,const std::string& method
    ,cudaError_t cuerr,std::string msg)const;
  void CheckCudaErroor(const std::string& srcfile,int srcline
    ,const std::string& classname,const std::string& method
    ,std::string msg)const;
 #endif

public:

  ///Types of gauges.
  typedef enum{ 
     GAUGE_Vel
    ,GAUGE_Swl
    ,GAUGE_MaxZ
    ,GAUGE_Mesh   //<vs_meeshdat>
    ,GAUGE_Force
  }TpGauge;

  ///Structure with data to gauges execution on CPU.
  typedef struct StrDataCpu{
    //-Divide data.
    bool divstate;   ///<Divide state is ready.
    double timestep;
    StDivDataCpu dvd;
    unsigned npbok;
    unsigned npb;
    unsigned np;
    //-Arrays with particle data.
    const acdouble3*  pos_c;
    const actypecode* code_c;
    const acuint*     idp_c;
    const acfloat4*   velrho_c;
    //-Methods.
    StrDataCpu(){ 
      divstate=false;
      npbok=npb=np=0;
      ConfigArrays(NULL,NULL,NULL,NULL);
    }
    void ConfigArrays(const acdouble3* pos,const actypecode* code
      ,const acuint* idp,const acfloat4* velrho)
    {
      this->pos_c=pos;  this->code_c=code;  
      this->idp_c=idp;  this->velrho_c=velrho;
    }
    void SetDivState(double timestep,const StDivDataCpu& dvd
      ,unsigned npbok,unsigned npb,unsigned np)
    {
      this->divstate=true;
      this->timestep=timestep;   this->dvd=dvd;
      this->npbok=npbok;  this->npb=npb;  this->np=np;  
    }
    void ClearDivState(){ this->divstate=false; }
  }StDataCpu;

 #ifdef _WITHGPU
  ///Structure with data to gauges execution on GPU.
  typedef struct StrDataGpu{
    //-Divide data.
    bool divstate;    ///<Divide state is ready.
    double timestep;
    StDivDataGpu dvd;
    unsigned npbok;
    unsigned npb;
    unsigned np;
    //-Arrays with particle data.
    const agdouble2*  posxy_g;
    const agdouble*   posz_g;
    const agtypecode* code_g;
    const aguint*     idp_g;
    const agfloat4*   velrho_g;
    //-Methods.
    StrDataGpu(){ 
      divstate=false;
      npbok=npb=np=0;
      ConfigArrays(NULL,NULL,NULL,NULL,NULL);
    }
    void ConfigArrays(const agdouble2* posxy,const agdouble* posz
      ,const agtypecode* code,const aguint* idp,const agfloat4* velrho)
    {
      this->posxy_g=posxy;  this->posz_g=posz;
      this->code_g=code;    this->idp_g=idp;
      this->velrho_g=velrho;
    }
    void SetDivState(double timestep,const StDivDataGpu& dvd
      ,unsigned npbok,unsigned npb,unsigned np)
    {
      this->divstate=true;
      this->timestep=timestep;   this->dvd=dvd;
      this->npbok=npbok;  this->npb=npb;  this->np=np;  
    }
    void ClearDivState(){ this->divstate=false; }
  }StDataGpu;
 #endif

  ///Structure with default configuration for JGaugeItem objects.
  typedef struct{
    bool savevtkpart;
    double computedt;
    double computestart;
    double computeend;
    bool output;
    double outputdt;
    double outputstart;
    double outputend;
  }StDefault;

protected:
  JLog2* Log;
  const bool Cpu;
  const int GpuCount;  ///<Number of GPUs (units) in use.
  std::string FileInfo;

  //-Constant values for calculation (they are constant).
  StCteSph CSP;        ///<Structure with main SPH constants values and configurations.
  float Scell;         ///<Cell size: KernelSize/ScellDiv (KernelSize or KernelSize/2).
  int ScellDiv;        ///<Value to divide KernelSize (1 or 2).
  tdouble3 MapPosMin;  ///<Lower limit of simulation + edge (KernelSize) if periodic conditions. MapPosMin=MapRealPosMin-KernelSize(in periodic axis) | Limite inferior de simulacion + borde (KernelSize) si hay condiciones periodicas.
  tdouble3 DomPosMin;  ///<Lower limit of simulation + edge (KernelSize) if periodic conditions. DomPosMin=Map_PosMin+(DomCelIni*Scell); | Limite inferior de simulacion + borde (KernelSize) si hay condiciones periodicas. 
  tdouble3 DomPosMax;  ///<Upper limit of simulation + edge (KernelSize) if periodic conditions. DomPosMax=min(Map_PosMax,Map_PosMin+(DomCelFin*Scell)); | Limite inferior de simulacion + borde (KernelSize) si hay condiciones periodicas. 

  //-Configuration on limits of calculation area.
  bool   FixedDomMCel; ///<Fixed calculation area.
  tuint3 DomMCelIni0;  ///<Initial first cell within the Map defining calculation area.
  tuint3 DomMCelFin0;  ///<Initial last cell within the Map defining calculation area.

  //-Configuration variables.
  bool SaveVtkPart; //-Creates VTK files for each PART.
  double ComputeDt;
  double ComputeStart;
  double ComputeEnd;
  double ComputeNext;
  bool OutputSave;
  double OutputDt;
  double OutputStart;
  double OutputEnd;
  double OutputNext;

  //-Results of measurement.
  double TimeStep;

  //-Variables to store the results in buffer.
  const unsigned OutSize; ///<Maximum number of results in buffer.
  unsigned OutCount;      ///<Number of stored results in buffer.
  std::string OutFile;

  JGaugeItem(TpGauge type,unsigned idx,std::string name,int gpucount
    ,unsigned outsize=200);
  void Reset();
  void SetTimeStep(double timestep);

  bool PointIsOut(double px,double py,double pz,const tdouble3& domposmin,const tdouble3& domposmax)const{ 
    return(px!=px || py!=py || pz!=pz || 
           px<domposmin.x  || py<domposmin.y  || pz<domposmin.z || 
           px>=domposmax.x || py>=domposmax.y || pz>=domposmax.z);
  }
  bool PointIsOut(double px,double py,double pz)const{ 
    return(px!=px || py!=py || pz!=pz || 
           px<DomPosMin.x  || py<DomPosMin.y  || pz<DomPosMin.z || 
           px>=DomPosMax.x || py>=DomPosMax.y || pz>=DomPosMax.z);
  }
  bool PointIsOut(double px,double py)const{ 
    return(px!=px || py!=py || px<DomPosMin.x || py<DomPosMin.y || 
           px>=DomPosMax.x  || py>=DomPosMax.y); 
  }

  static std::string GetNameType(TpGauge type);

  tuint3 CalcMCelFromPos(const tdouble3& ps)const;
  void   CalcMCelIniFinFromPos(const tdouble3& ps,tuint3& mcelini,tuint3& mcelfin)const;
  void   CalcMCelIniFinFromPos(const tdouble3& psmin,const tdouble3& psmax,tuint3& mcelini,tuint3& mcelfin)const;

  virtual void ClearResult()=0;
  virtual void StoreResult()=0;

public:
  const TpGauge Type;
  const unsigned Idx;
  const std::string Name;

  void Config(const StCteSph& csp,float scell,int scelldiv
    ,tdouble3 mapposmin,tdouble3 domposmin,tdouble3 domposmax);
  virtual void ConfigDomMCel(bool fixed)=0;
  void SetSaveVtkPart(bool save){ SaveVtkPart=save; }
  void ConfigComputeTiming(double start,double end,double dt);
  void ConfigOutputTiming(bool save,double start,double end,double dt);

  void GetConfig(std::vector<std::string>& lines)const;

  std::string GetResultsFile(bool dirgauges,const std::string& fext,const std::string& subname="")const;
  std::string GetResultsFileCsv(const std::string& subname="")const;
  std::string GetResultsFileVtk(const std::string& subname="")const;
  virtual void SaveResults()=0;
  virtual void SaveVtkResult(unsigned cpart)=0;
  virtual unsigned GetPointDef(std::vector<tfloat3>& points)const=0;
  virtual void SaveVtkScheme()const{};

  void SaveResults(unsigned cpart);

  double GetComputeDt()const{ return(ComputeDt); }
  double GetComputeStart()const{ return(ComputeStart); }
  double GetComputeEnd()const{ return(ComputeEnd); }

  double GetOutputDt()const{ return(OutputDt); }
  double GetOutputStart()const{ return(OutputStart); }
  double GetOutputEnd()const{ return(OutputEnd); }

  bool Update(double timestep)const{ return(timestep>=ComputeNext && ComputeStart<=timestep && timestep<=ComputeEnd); }
  bool Output(double timestep)const{ return(OutputSave && timestep>=OutputNext && OutputStart<=timestep && timestep<=OutputEnd); }

  virtual void CalculeCpu(const StDataCpu& datacpu)=0;

 #ifdef _WITHGPU
  virtual bool AllocatedGpuMemory(int id)const=0;
  virtual void FreeGpuMemory(int id)=0;
  virtual void AllocGpuMemory(int id)=0;
  virtual void CalculeGpu(const StDataGpu& datagpu)=0;
 #endif
};


//##############################################################################
//# JGaugeVelocity
//##############################################################################
/// \brief Calculates velocity in fluid domain.
class JGaugeVelocity : public JGaugeItem
{
private:
 #ifdef _WITHGPU
  ///Structure with auxiliary memory for execution on GPU.
  typedef struct StrGaugeVelDataGpu{
    bool GpuMemory;     ///<Indicates when GPU memory is allocated.
    float3* Resultg;    ///<Stores final result from GPU [1].
    //-Methods.
    StrGaugeVelDataGpu(){
      GpuMemory=false;
      Resultg=NULL;
    }
  }StGaugeVelDataGpu;
 #endif

public:
  ///Structure with result of JGaugeVelocity object.
  typedef struct StrGaugeVelRes{
    double timestep;
    tfloat3 point;
    tfloat3 vel;
    bool modified;
    StrGaugeVelRes(){ Reset(); }
    void Reset(){
      Set(0,TFloat3(0),TFloat3(0));
      modified=false;
    }
    void Set(double t,const tfloat3& pt,const tfloat3& vel){
      timestep=t; point=pt; this->vel=vel; modified=true;
    }
  }StGaugeVelRes;

protected:
  //-Definition.
  tdouble3 Point;

  //-Auxiliary variables for GPU execution.
 #ifdef _WITHGPU
  StGaugeVelDataGpu AuxDataGpu[MAXGPUS];
 #endif

  //-Result variables.
  StGaugeVelRes Result;      ///<Result of the last measure.
  std::vector<StGaugeVelRes> OutBuff; ///<Results in buffer.

  void Reset();
 #ifdef _WITHGPU
  void ResetGpuMemory();
 #endif

  void ClearResult(){ Result.Reset(); }
  void StoreResult();

public:
  JGaugeVelocity(unsigned idx,std::string name,tdouble3 point
    ,int gpucount);
  ~JGaugeVelocity();
  void ConfigDomMCel(bool fixed);

  void SaveResults();
  void SaveVtkResult(unsigned cpart);
  unsigned GetPointDef(std::vector<tfloat3>& points)const;

  tdouble3 GetPoint()const{ return(Point); }
  const StGaugeVelRes& GetResult()const{ return(Result); }

  void SetPoint(const tdouble3& point){ ClearResult(); Point=point; }

  template<TpKernel tker> void CalculeCpuT(const StDataCpu& datacpu);
  void CalculeCpu(const StDataCpu& datacpu);

 #ifdef _WITHGPU
  bool AllocatedGpuMemory(int id)const;
  void FreeGpuMemory(int id);
  void AllocGpuMemory(int id);
  void CalculeGpu(const StDataGpu& datagpu);
 #endif
};


//##############################################################################
//# JGaugeSwl
//##############################################################################
/// \brief Calculates Surface Water Level in fluid domain.
class JGaugeSwl : public JGaugeItem
{
private:
 #ifdef _WITHGPU
  ///Structure with auxiliary memory for execution on GPU.
  typedef struct StrGaugeSwlDataGpu{
    bool GpuMemory;     ///<Indicates when GPU memory is allocated.
    float3* Resultg;    ///<Stores final result from GPU [1].
    //-Methods.
    StrGaugeSwlDataGpu(){ 
      GpuMemory=false;
      Resultg=NULL;
    }
  }StGaugeSwlDataGpu;
 #endif

public:
  ///Structure with result of JGaugeVelocity object.
  typedef struct StrGaugeSwlRes{
    double timestep;
    tfloat3 point0;
    tfloat3 point2;
    tfloat3 posswl;
    bool modified;
    StrGaugeSwlRes(){ Reset(); }
    void Reset(){
      Set(0,TFloat3(0),TFloat3(0),TFloat3(0));
      modified=false;
    }
    void Set(double t,const tfloat3& pt0,const tfloat3& pt2,const tfloat3& ps){
      timestep=t; point0=pt0; point2=pt2; posswl=ps; modified=true;
    }
  }StGaugeSwlRes;

protected:
  //-Definition.
  tdouble3 Point0;   ///<Initial point.
  tdouble3 Point2;   ///<Final point.
  double PointDp;    ///<Distance between points.
  float MassLimit; 
  //-Auxiliary variables.
  unsigned PointNp;  ///<Number of points.
  tdouble3 PointDir; ///<Vector to compute points.

  //-Auxiliary variables for GPU execution.
 #ifdef _WITHGPU
  StGaugeSwlDataGpu AuxDataGpu[MAXGPUS];
 #endif

  //-Result variables.
  StGaugeSwlRes Result;      ///<Result of the last measure.
  std::vector<StGaugeSwlRes> OutBuff; ///<Results in buffer.

  void Reset();
 #ifdef _WITHGPU
  void ResetGpuMemory();
 #endif

  void ClearResult(){ Result.Reset(); }
  void StoreResult();
  template<TpKernel tker> float CalculeMassCpu(const tdouble3& ptpos,const StDivDataCpu& dvd
    ,const tdouble3* pos,const typecode* code,const tfloat4* velrho)const;

public:
  JGaugeSwl(unsigned idx,std::string name,tdouble3 point0,tdouble3 point2
    ,double pointdp,float masslimit,int gpucount);
  ~JGaugeSwl();
  void ConfigDomMCel(bool fixed);

  void SaveResults();
  void SaveVtkResult(unsigned cpart);
  unsigned GetPointDef(std::vector<tfloat3>& points)const;

  tdouble3 GetPoint0()const{ return(Point0); }
  tdouble3 GetPoint2()const{ return(Point2); }
  double GetPointDp()const{ return(PointDp); }
  float GetMassLimit()const{ return(MassLimit); }
  const StGaugeSwlRes& GetResult()const{ return(Result); }

  void SetPoints(const tdouble3& point0,const tdouble3& point2,double pointdp=0);

  template<TpKernel tker> void CalculeCpuT(const StDataCpu& datacpu);
  void CalculeCpu(const StDataCpu& datacpu);

 #ifdef _WITHGPU
  bool AllocatedGpuMemory(int id)const;
  void FreeGpuMemory(int id);
  void AllocGpuMemory(int id);
  void CalculeGpu(const StDataGpu& datagpu);
 #endif
};


//##############################################################################
//# JGaugeMaxZ
//##############################################################################
/// \brief Calculates maximum z of fluid at distance of a vertical line.
class JGaugeMaxZ : public JGaugeItem
{
private:
 #ifdef _WITHGPU
  ///Structure with auxiliary memory for execution on GPU.
  typedef struct StrGaugeMaxzDataGpu{
    bool GpuMemory;     ///<Indicates when GPU memory is allocated.
    float3* Resultg;    ///<Stores final result from GPU [1].
    //-Methods.
    StrGaugeMaxzDataGpu(){ 
      GpuMemory=false;
      Resultg=NULL;
    }
  }StGaugeMaxzDataGpu;
 #endif

public:
  ///Structure with result of JGaugeMaxZ object.
  typedef struct StrGaugeMaxzRes{
    double timestep;
    tfloat3 point0;
    float zmax;
    bool modified;
    StrGaugeMaxzRes(){ Reset(); }
    void Reset(){
      Set(0,TFloat3(0),0);
      modified=false;
    }
    void Set(double t,const tfloat3& pt,float z){
      timestep=t; point0=pt; zmax=z; modified=true;
    }
  }StGaugeMaxzRes;

protected:
  //-Definition.
  tdouble3 Point0;
  double Height;
  float DistLimit;

  //-Auxiliary variables for GPU execution.
 #ifdef _WITHGPU
  StGaugeMaxzDataGpu AuxDataGpu[MAXGPUS];
 #endif

  //-Result variables.
  StGaugeMaxzRes Result;      ///<Result of the last measure.
  std::vector<StGaugeMaxzRes> OutBuff; ///<Results in buffer.

  void Reset();
 #ifdef _WITHGPU
  void ResetGpuMemory();
 #endif

  void ClearResult(){ Result.Reset(); }
  void StoreResult();
  void GetInteractionCellsMaxZ(const tdouble3& pos,const tint4& nc
    ,const tdouble3& domposmin,const tint3& cellzero
    ,int& cxini,int& cxfin,int& yini,int& yfin,int& zini,int& zfin)const;

public:
  JGaugeMaxZ(unsigned idx,std::string name,tdouble3 point0,double height
    ,float distlimit,int gpucount);
  ~JGaugeMaxZ();
  void ConfigDomMCel(bool fixed);

  void SaveResults();
  void SaveVtkResult(unsigned cpart);
  unsigned GetPointDef(std::vector<tfloat3>& points)const;

  tdouble3 GetPoint0()const{ return(Point0); }
  double GetHeight()const{ return(Height); }
  float GetDistLimit()const{ return(DistLimit); }
  const StGaugeMaxzRes& GetResult()const{ return(Result); }

  void SetPoint0   (const tdouble3& point0){ ClearResult(); Point0=point0; }
  void SetHeight   (double height){          ClearResult(); Height=height; }
  void SetDistLimit(float distlimit){        ClearResult(); DistLimit=distlimit; }

   void CalculeCpu(const StDataCpu& datacpu);

 #ifdef _WITHGPU
  bool AllocatedGpuMemory(int id)const;
  void FreeGpuMemory(int id);
  void AllocGpuMemory(int id);
  void CalculeGpu(const StDataGpu& datagpu);
 #endif
};


//<vs_meeshdat_ini>
//##############################################################################
//# JGaugeMesh
//##############################################################################
/// \brief Calculates Velocity and Surface Water Level in a grid of positions.
class JGaugeMesh : public JGaugeItem
{
private:
 #ifdef _WITHGPU
  ///Structure with auxiliary memory for execution on GPU.
  typedef struct StrGaugeMeshDataGpu{
    bool GpuMemory;     ///<Indicates when GPU memory is allocated.
    //float3* Resultg;  ///<Stores final result from GPU [1].
    float*  DataRhopg;  ///<Stores on GPU memory density. [GridPts.npt1*GridPts.npt2*GridPts.npt3]
    float3* DataVxyzg;  ///<Stores on GPU memory velocity in X, Y, Z. [GridPts.npt1*GridPts.npt2*GridPts.npt3]
    float*  DataVdirg;  ///<Stores on GPU memory velocity in requested direction. [GridPts.npt1*GridPts.npt2*GridPts.npt3]
    float*  DataZsurfg; ///<Stores on GPU memory Z surface water. [GridPts.npt1*GridPts.npt2]
    float*  DataMassg;  ///<Stores on GPU memory mass value. [GridPts.npt1*GridPts.npt2*GridPts.npt3]
    //-Methods.
    StrGaugeMeshDataGpu(){ 
      GpuMemory=false;
      //Resultg=NULL;
      DataRhopg=NULL;
      DataVxyzg=NULL;
      DataVdirg=NULL;
      DataZsurfg=NULL;
      DataMassg=NULL;
    }
  }StGaugeMeshDataGpu;
 #endif

public:
  ///Structure with basic information about configuration.
  typedef struct{
    tdouble3 ptref;      ///<Initial measurement position.
    tdouble3 ptend;      ///<Final measurement position.
    tdouble3 vec1;       ///<First axis vector to define the measurement grid.
    tdouble3 vec2;       ///<Second axis vector to define the measurement grid.
    tdouble3 vec3;       ///<Third axis vector to define the measurement grid (used for swl calculation).
    tdouble3 dispt;      ///<Distance between measurement points.
    tuint4   npt;        ///<Number of positions.
    tfloat3 dirdat;      ///<Direction vector for computed linear velocity or other variables.
    std::string outdata; ///<Output data selection.
  }StInfo;

  ///Structure with result of JGaugeMesh object.
  typedef struct StrMeshRes{
    double timestep;
    jmsh::JMeshData* meshdat;
    bool modified;
    StrMeshRes(){ Reset(); }
    void Reset(){
      Set(0,NULL);
      modified=false;
    }
    void Set(double t,jmsh::JMeshData* meshdat){
      timestep=t; this->meshdat=meshdat; modified=true;
    }
  }StMeshRes;

protected:
  //-Definition.
  std::string OutDataList;   ///<Output data selection.
  bool ComputeVelxyz;        ///<Stores velocity in X, Y, Z.
  bool ComputeVeldir;        ///<Stores velocity in requested direction.
  bool ComputeRhop;          ///<Stores density.
  bool ComputeZsurf;         ///<Stores Z surface water.
  bool SaveCsv;              ///<Saves data in CSV file.
  bool SaveBin;              ///<Saves data in binary file.
  jmsh::StMeshBasic MeshBas; ///<Basic mesh configuration.
  jmsh::StMeshPts MeshPts;   ///<Mesh configuration for calculations. 
  float KcLimit;             ///<Minimum value of sum_wab_vol to apply the Kernel Correction. FLT_MAX to disable (default=0.5)
  float KcDummy;             ///<Dummy value for non-corrected values. FLT_MAX to disable (default=0)
  float MassLimit;

  //-Auxiliary variables.
  jmsh::JMeshData* MeshDat;
  float* MassDatCpu;   ///<Auxiliar memory to compute mass. [MeshPts.npt]

  //-Auxiliary variables for GPU execution.
 #ifdef _WITHGPU
  StGaugeMeshDataGpu AuxDataGpu[MAXGPUS];
 #endif

  //-Result variables.
  StrMeshRes Result;      ///<Result of the last measure.
  std::vector<StrMeshRes> OutBuff; ///<Results in buffer.

  jmsh::JMeshTDatasSave* MeshDataSave;  ///<Saves data in binary file.

  void Reset();
 #ifdef _WITHGPU
  void ResetGpuMemory();
 #endif

  void SetMeshData(const jmsh::StMeshBasic& meshbas,std::string outdata);
  void ClearResult(){ Result.Reset(); }
  void StoreResult();

public:
  JGaugeMesh(unsigned idx,std::string name,const jmsh::StMeshBasic& meshbas
    ,std::string outdata,unsigned tfmt,unsigned buffersize
    ,float kclimit,float kcdummy,float masslimit,int gpucount);
  ~JGaugeMesh();
  void ConfigDomMCel(bool fixed);

  static std::string CorrectDataList(std::string datalist);
  void ConfigDataList(std::string datalist);
  std::string GetDataList()const;

  StInfo GetInfo()const;

  void SaveResults();
  void SaveVtkResult(unsigned cpart);
  unsigned GetPointDef(std::vector<tfloat3>& points)const;
  void SaveVtkScheme()const;

  jmsh::StMeshBasic GetMeshBas()const{ return(MeshBas); }
  jmsh::StMeshPts GetMesh()const{ return(MeshPts); }
  float GetMassLimit()const{ return(MassLimit); }
  float GetKcLimit()const{ return(KcLimit); }
  float GetKcDummy()const{ return(KcDummy); }

  const StMeshRes& GetResult()const{ return(Result); }

  const float* GetPtrDataZsurf()const;
 #ifdef _WITHGPU
  const float* GetPtrDataZsurfg(int id)const;
 #endif

  template<TpKernel tker> void CalculeCpuT(const StDataCpu& datacpu);
  void CalculeCpu(const StDataCpu& datacpu);

 #ifdef _WITHGPU
  bool AllocatedGpuMemory(int id)const;
  void FreeGpuMemory(int id);
  void AllocGpuMemory(int id);
  void CalculeGpu(const StDataGpu& datagpu);
 #endif

};
//<vs_meeshdat_end>


//##############################################################################
//# JGaugeForce
//##############################################################################
/// \brief Calculates force sumation on selected particles (using only fluid particles).
class JGaugeForce : public JGaugeItem
{
private:
 #ifdef _WITHGPU
  ///Structure with auxiliary memory for execution on GPU.
  typedef struct StrGaugeForceDataGpu{
    bool GpuMemory;     ///<Indicates when GPU memory is allocated.
    float3* Resultg;    ///<Stores final result from GPU [1].
    float3* PartAceg;   ///<Ace of particles [Count].
    float3* AuxSumg;    ///<Used for Ace reduction [size reduction].
    //-Methods.
    StrGaugeForceDataGpu(){ 
      GpuMemory=false;
      Resultg=NULL;
      PartAceg=NULL;
      AuxSumg=NULL;
    }
  }StGaugeForceDataGpu;
 #endif

public:
  ///Structure with result of JGaugeForce object.
  typedef struct StrGaugeForceRes{
    double timestep;
    word mkbound;
    tfloat3 force;
    bool modified;
    StrGaugeForceRes(){ Reset(); }
    void Reset(){
      Set(0,TFloat3(0));
      modified=false;
    }
    void Set(double t,const tfloat3& forceres){
      timestep=t; force=forceres; modified=true;
    }
  }StGaugeForceRes;

protected:
  //-Definition.
  word MkBound;
  TpParticles TypeParts;
  unsigned IdBegin;
  unsigned Count;
  typecode Code;
  tfloat3 InitialCenter;

  //-Auxiliary variables.
  tfloat3* PartAcec;

  //-Auxiliary variables for GPU execution.
 #ifdef _WITHGPU
  StGaugeForceDataGpu AuxDataGpu[MAXGPUS];
 #endif

  //-Result variables.
  StGaugeForceRes Result;      ///<Result of the last measure.
  std::vector<StGaugeForceRes> OutBuff; ///<Results in buffer.

  void Reset();
 #ifdef _WITHGPU
  void ResetGpuMemory();
 #endif

  void ClearResult(){ Result.Reset(); }
  void StoreResult();

public:
  JGaugeForce(unsigned idx,std::string name,word mkbound
    ,TpParticles typeparts,unsigned idbegin,unsigned count
    ,typecode code,tfloat3 center,int gpucount);
  ~JGaugeForce();
  void ConfigDomMCel(bool fixed);

  void SaveResults();
  void SaveVtkResult(unsigned cpart);
  unsigned GetPointDef(std::vector<tfloat3>& points)const;

  word        GetMkBound()  const{ return(MkBound); }
  TpParticles GetTypeParts()const{ return(TypeParts); }
  unsigned    GetIdBegin()  const{ return(IdBegin); }
  unsigned    GetCount()    const{ return(Count); }
  typecode    GetCode()     const{ return(Code); }
  tfloat3     GetInitialCenter()const{ return(InitialCenter); }
  const StGaugeForceRes& GetResult()const{ return(Result); }

  template<TpKernel tker> void CalculeCpuT(const StDataCpu& datacpu);
  void CalculeCpu(const StDataCpu& datacpu);

 #ifdef _WITHGPU
  bool AllocatedGpuMemory(int id)const;
  void FreeGpuMemory(int id);
  void AllocGpuMemory(int id);
  void CalculeGpu(const StDataGpu& datagpu);
 #endif
};


#endif

