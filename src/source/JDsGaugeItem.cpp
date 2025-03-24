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

/// \file JDsGaugeItem.cpp \brief Implements the class \ref JGauge.

#include "JDsGaugeItem.h"
#include "JCellSearch_inline.h"
#include "FunSphKernel.h"
#include "FunSphEos.h"
#include "JException.h"
#include "JLog2.h"
#include "JSaveCsv2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JDataArrays.h"
#include "JSpVtkData.h"
#include "JMeshData.h"        //<vs_meeshdat>
#include "JMeshTDatasSave.h"  //<vs_meeshdat>
#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JDsGauge_ker.h"
  #include "JReduSum_ker.h"
#endif
#include <cmath>
#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JGaugeItem
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeItem::JGaugeItem(TpGauge type,unsigned idx,std::string name,int gpucount
  ,unsigned outsize):Log(AppInfo.LogPtr()),Type(type),Idx(idx),Name(name)
  ,Cpu(gpucount==0),GpuCount(gpucount),OutSize(outsize)
{
  ClassName="JGaugeItem";
  Reset();
}

#ifdef _WITHGPU
//==============================================================================
/// Throws exception related to a CUDA error.
//==============================================================================
void JGaugeItem::RunExceptioonCuda(const std::string& srcfile,int srcline
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
void JGaugeItem::CheckCudaErroor(const std::string& srcfile,int srcline
  ,const std::string& classname,const std::string& method
  ,std::string msg)const
{
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess)RunExceptioonCuda(srcfile,srcline,classname,method,cuerr,msg);
}
#endif

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeItem::Reset(){
  Config(CteSphNull(),0,0,TDouble3(0),TDouble3(0),TDouble3(0));
  FixedDomMCel=false;
  DomMCelIni0=DomMCelFin0=TUint3(0);
  SaveVtkPart=false;
  ConfigComputeTiming(0,0,0);
  ConfigOutputTiming(false,0,0,0);
  TimeStep=0;
  OutCount=0;
  OutFile="";
}

//==============================================================================
/// Configures object.
//==============================================================================
void JGaugeItem::Config(const StCteSph& csp,float scell,int scelldiv
  ,tdouble3 mapposmin,tdouble3 domposmin,tdouble3 domposmax)
{
  CSP=csp;
  Scell=scell;
  ScellDiv=scelldiv;
  MapPosMin=mapposmin;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
}

//==============================================================================
/// Configures compute timing.
//==============================================================================
void JGaugeItem::ConfigComputeTiming(double start,double end,double dt){
  ComputeDt=dt;
  ComputeStart=start;
  ComputeEnd=end;
  ComputeNext=0;
}

//==============================================================================
/// Configures output timing.
//==============================================================================
void JGaugeItem::ConfigOutputTiming(bool save,double start,double end,double dt){
  OutputSave=save;
  OutputDt=dt;
  OutputStart=start;
  OutputEnd=end;
  OutputNext=0;
}

//==============================================================================
/// Returns type in string.
//==============================================================================
std::string JGaugeItem::GetNameType(TpGauge type){
  switch(type){
    case GAUGE_Vel:     return("Vel");
    case GAUGE_Swl:     return("SWL");
    case GAUGE_MaxZ:    return("MaxZ");
    case GAUGE_Mesh:    return("Mesh");  //<vs_meeshdat>
    case GAUGE_Force:   return("Force");
  }
  return("???");
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JGaugeItem::GetConfig(std::vector<std::string>& lines)const{
  lines.push_back(fun::PrintStr("Type.......: %s"
    ,JGaugeItem::GetNameType(Type).c_str()));
  const string cpend=fun::DoublexStr(ComputeEnd,"%g");
  lines.push_back(fun::PrintStr("Compute....: %g - %s   dt:%g"
    ,ComputeStart,cpend.c_str(),ComputeDt));
  const string ouend=fun::DoublexStr(OutputEnd,"%g");
  lines.push_back(fun::PrintStr("Output.....:%s%g - %s   dt:%g"
    ,(OutputSave? " ": " (disabled)  "),OutputStart ,ouend.c_str(),OutputDt));
  lines.push_back(fun::PrintStr("SaveVtkPart: %s",(SaveVtkPart? "True": "False")));
  if(Type!=GAUGE_Force){
    //-Computes map-cell max according to final map-cell.
    const tuint3 dcelmax=MinValues(DomMCelFin0,DomMCelFin0-TUint3(1));
    lines.push_back(fun::PrintStr("DomainCells: %s   mode:%s"
      ,fun::Uint3RangeStr(DomMCelIni0,dcelmax).c_str()
      ,(FixedDomMCel? "fixed": "unfixed")));
  }
  if(Type==GAUGE_Vel){
    const JGaugeVelocity* gau=(JGaugeVelocity*)this;
    lines.push_back(fun::PrintStr("Point......: (%g,%g,%g)",gau->GetPoint().x,gau->GetPoint().y,gau->GetPoint().z));
  }
  else if(Type==GAUGE_Swl){
    const JGaugeSwl* gau=(JGaugeSwl*)this;
    lines.push_back(fun::PrintStr("MassLimit..: %g",gau->GetMassLimit()));
    lines.push_back(fun::PrintStr("GaugePoints: %s",fun::Double3gRangeStr(gau->GetPoint0(),gau->GetPoint2()).c_str()));
    lines.push_back(fun::PrintStr("PointDp....: %g",gau->GetPointDp()));
  }
  else if(Type==GAUGE_MaxZ){
    const JGaugeMaxZ* gau=(JGaugeMaxZ*)this;
    lines.push_back(fun::PrintStr("Point0.....: (%g,%g,%g)   Height:%g",gau->GetPoint0().x,gau->GetPoint0().y,gau->GetPoint0().z,gau->GetHeight()));
    lines.push_back(fun::PrintStr("DistLimit..: %g",gau->GetDistLimit()));
  }
  else if(Type==GAUGE_Mesh){  //<vs_meeshdat_ini>
    const JGaugeMesh* gau=(JGaugeMesh*)this;
    lines.push_back(fun::PrintStr("MassLimit..: %g",gau->GetMassLimit()));
    JGaugeMesh::StInfo v=gau->GetInfo();
    lines.push_back(fun::PrintStr("GaugePoints: %u  (%u x %u x %u)",v.npt.w,v.npt.x,v.npt.y,v.npt.z));
    lines.push_back(fun::PrintStr("GaugePos...: %s",fun::Double3gRangeStr(v.ptref,v.ptend).c_str()));
    lines.push_back(fun::PrintStr("Vectors....: (%s) - (%s) - (%s)",fun::Double3gStr(v.vec1).c_str(),fun::Double3gStr(v.vec2).c_str(),fun::Double3gStr(v.vec3).c_str()));
    lines.push_back(fun::PrintStr("PointDp....: (%s)",fun::Double3gStr(v.dispt).c_str()));
    lines.push_back(fun::PrintStr("DirDat.....: (%s)",fun::Float3gStr(v.dirdat).c_str()));
    lines.push_back(fun::PrintStr("Compute....: %s",v.outdata.c_str()));
    string tkcorr="none";
    if(gau->GetKcLimit()!=FLT_MAX){
      tkcorr=fun::PrintStr("ksupport >= %g",gau->GetKcLimit());
      if(gau->GetKcDummy()!=FLT_MAX)tkcorr=tkcorr+fun::PrintStr(" (dummy:%g)",gau->GetKcDummy());
      else tkcorr=tkcorr+" (dummy:none)";
    }
    lines.push_back(fun::PrintStr("KernelCorr.: %s",tkcorr.c_str()));
  }  //<vs_meeshdat_end>
  else if(Type==GAUGE_Force){
    const JGaugeForce* gau=(JGaugeForce*)this;
    lines.push_back(fun::PrintStr("MkBound.....: %u (%s particles)",gau->GetMkBound(),TpPartGetStrCode(gau->GetTypeParts())));
    lines.push_back(fun::PrintStr("Particles id: %u - %u",gau->GetIdBegin(),gau->GetIdBegin()+gau->GetCount()-1));
  }
  else Run_Exceptioon("Type unknown.");
}

//==============================================================================
/// Returns filename for output results files.
//==============================================================================
std::string JGaugeItem::GetResultsFile(bool dirgauges,const std::string& fext
  ,const std::string& subname)const
{
  //const string fdir=(dirgauges? AppInfo.GetDirDataOut(): AppInfo.GetDirOut());
  const string fdir=(dirgauges? AppInfo.GetDirOut()+"gaugesvtks/": AppInfo.GetDirOut());
  return(fdir+"Gauges"+GetNameType(Type)+"_"+Name+subname+"."+fun::StrLower(fext));
}

//==============================================================================
/// Returns filename for output results for CSV files.
//==============================================================================
std::string JGaugeItem::GetResultsFileCsv(const std::string& subname)const{
  return(GetResultsFile(false,"csv",subname));
}

//==============================================================================
/// Returns filename for output results for VTK files.
//==============================================================================
std::string JGaugeItem::GetResultsFileVtk(const std::string& subname)const{
  return(GetResultsFile(true,"vtk",subname));
}

//==============================================================================
/// Updates TimeStep and calculates ComputeNext.
//==============================================================================
void JGaugeItem::SetTimeStep(double timestep){
  if(ComputeDt){
    unsigned nt=unsigned(timestep/ComputeDt);
    ComputeNext=ComputeDt*nt;
    if(ComputeNext<=timestep)ComputeNext=ComputeDt*(nt+1);
  }
  TimeStep=timestep;
}

//==============================================================================
/// Saves results in VTK and/or CSV file.
//==============================================================================
void JGaugeItem::SaveResults(unsigned cpart){
  SaveResults();
  if(SaveVtkPart)SaveVtkResult(cpart);
}

//==============================================================================
/// Returns MapCell of position. Note that the minimum cell is always zero.
//==============================================================================
tuint3 JGaugeItem::CalcMCelFromPos(const tdouble3& ps)const{
  const double dx=ps.x-MapPosMin.x;
  const double dy=ps.y-MapPosMin.y;
  const double dz=ps.z-MapPosMin.z;
  const unsigned cx=(dx>=0? unsigned(dx/Scell): 0);
  const unsigned cy=(dy>=0? unsigned(dy/Scell): 0);
  const unsigned cz=(dz>=0? unsigned(dz/Scell): 0);
  return(TUint3(cx,cy,cz));
}

//==============================================================================
/// Returns MapCell domain of one position.
//==============================================================================
void JGaugeItem::CalcMCelIniFinFromPos(const tdouble3& ps,tuint3& mcelini
  ,tuint3& mcelfin)const
{
  const unsigned ncel=ScellDiv;
  const tuint3 mcel=CalcMCelFromPos(ps);
  const unsigned cxi=(mcel.x>ncel? mcel.x-ncel: 0);
  const unsigned cyi=(mcel.y>ncel? mcel.y-ncel: 0);
  const unsigned czi=(mcel.z>ncel? mcel.z-ncel: 0);
  mcelini=TUint3(cxi,cyi,czi);
  mcelfin=mcel+TUint3(ncel+1);
}

//==============================================================================
/// Returns MapCell domain of two positions.
//==============================================================================
void JGaugeItem::CalcMCelIniFinFromPos(const tdouble3& psmin
  ,const tdouble3& psmax,tuint3& mcelini,tuint3& mcelfin)const
{
  const unsigned ncel=ScellDiv;
  const tuint3 mcel=CalcMCelFromPos(psmin);
  const unsigned cxi=(mcel.x>ncel? mcel.x-ncel: 0);
  const unsigned cyi=(mcel.y>ncel? mcel.y-ncel: 0);
  const unsigned czi=(mcel.z>ncel? mcel.z-ncel: 0);
  mcelini=TUint3(cxi,cyi,czi);
  const tuint3 mcel2=CalcMCelFromPos(psmax);
  mcelfin=mcel2+TUint3(ncel+1);
}


//##############################################################################
//# JGaugeVelocity
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeVelocity::JGaugeVelocity(unsigned idx,std::string name,tdouble3 point
  ,int gpucount):JGaugeItem(GAUGE_Vel,idx,name,gpucount)
{
  ClassName="JGaugeVel";
  FileInfo=string("Saves velocity data measured from fluid particles (by ")+ClassName+").";
  Reset();
  SetPoint(point);
}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeVelocity::~JGaugeVelocity(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeVelocity::Reset(){
  SetPoint(TDouble3(0));
  #ifdef _WITHGPU
    ResetGpuMemory();
  #endif
  JGaugeItem::Reset();
}

#ifdef _WITHGPU
//==============================================================================
/// Reset GPU memory.
//==============================================================================
void JGaugeVelocity::ResetGpuMemory(){
  for(int g=0;g<GpuCount;g++){
    if(AuxDataGpu[g].GpuMemory)
      Run_Exceptioon(fun::PrintStr("Auxiliary GPU memory for unit %d is not free.",g));
    AuxDataGpu[g]=StrGaugeVelDataGpu();
  }
}
#endif

//==============================================================================
/// Configuration of initial limits of calculation area.
//==============================================================================
void JGaugeVelocity::ConfigDomMCel(bool fixed){
  if(!fixed)Run_Exceptioon("Unfixed calculation domain is invalid for Vel-Gauge.");
  FixedDomMCel=fixed;
  CalcMCelIniFinFromPos(Point,DomMCelIni0,DomMCelFin0);
}

//==============================================================================
/// Record the last measure result.
//==============================================================================
void JGaugeVelocity::StoreResult(){
  if(OutputSave){
    //-Allocates memory.
    while(unsigned(OutBuff.size())<OutSize)OutBuff.push_back(StrGaugeVelRes());
    //-Empty buffer.
    if(OutCount+1>=OutSize)SaveResults();
    //-Stores last results.
    OutBuff[OutCount]=Result;
    OutCount++;
    //-Updates OutputNext.
    if(OutputDt){
      const unsigned nt=unsigned(TimeStep/OutputDt);
      OutputNext=OutputDt*nt;
      if(OutputNext<=TimeStep)OutputNext=OutputDt*(nt+1);
    }
  }
}

//==============================================================================
/// Saves stored results in CSV file.
//==============================================================================
void JGaugeVelocity::SaveResults(){
  if(OutCount){
    const bool first=OutFile.empty();
    if(first){
      OutFile=GetResultsFileCsv();
      Log->AddFileInfo(OutFile,FileInfo);
    }
    jcsv::JSaveCsv2 scsv(OutFile,!first,AppInfo.GetCsvSepComa());
    //-Saves head.
    if(first){
      scsv.SetHead();
      scsv << "time [s];velx [m/s];vely [m/s];velz [m/s];posx [m];posy [m];posz [m]" << jcsv::Endl();
    }
    //-Saves data.
    scsv.SetData();
    scsv << jcsv::Fmt(jcsv::TpFloat1,"%g") << jcsv::Fmt(jcsv::TpFloat3,"%g;%g;%g");
    for(unsigned c=0;c<OutCount;c++){
      scsv << OutBuff[c].timestep << OutBuff[c].vel << OutBuff[c].point << jcsv::Endl();
    }
    OutCount=0;
  }
}

//==============================================================================
/// Saves last result in VTK file.
//==============================================================================
void JGaugeVelocity::SaveVtkResult(unsigned cpart){
  //-Prepares data.
  JDataArrays arrays;
  arrays.AddArray("Pos",1,&(Result.point),false);
  arrays.AddArray("Vel",1,&(Result.vel),false);
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  JSpVtkData::Save(fun::FileNameSec(GetResultsFileVtk(),cpart),arrays,"Pos");
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeVelocity::GetPointDef(std::vector<tfloat3>& points)const{
  points.push_back(ToTFloat3(Point));
  return(1);
}

//==============================================================================
/// Calculates velocity at indicated points (on CPU).
//==============================================================================
template<TpKernel tker> void JGaugeVelocity::CalculeCpuT(const StDataCpu& datacpu)
{
  //-Prepare input data.
  const double& timestep =datacpu.timestep;
  const StDivDataCpu& dvd=datacpu.dvd;
  const tdouble3* pos    =datacpu.pos_c->cptr();
  const typecode* code   =datacpu.code_c->cptr();
  //const unsigned* idp  =datacpu.idp_c->cptr();
  const tfloat4*  velrho =datacpu.velrho_c->cptr();
  SetTimeStep(timestep);
  //-Start measure.
  tfloat3 ptvel=TFloat3(0);
  const bool ptout=PointIsOut(Point.x,Point.y,Point.z);//-Verify that the point is within domain boundaries. | Comprueba que el punto este dentro de limites del dominio.
  if(!ptout){
    //-Auxiliary variables.
    double sumwab=0;
    tdouble3 sumvel=TDouble3(0);
    //-Search for fluid neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(Point,false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float rr2=nsearch::Distance2(Point,pos[p2]);
        //-Interaction with real neighbouring particles.
        if(rr2<=CSP.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
          float wab=fsph::GetKernel_Wab<tker>(CSP,rr2);
          const tfloat4 velrho2=velrho[p2];
          wab*=CSP.massfluid/velrho2.w;
          sumwab+=wab;
          sumvel.x+=wab*velrho2.x;
          sumvel.y+=wab*velrho2.y;
          sumvel.z+=wab*velrho2.z;
        }
      }
    }
    //-Applies kernel correction.
    //if(sumwab!=0){
    //  sumvel.x/=sumwab;
    //  sumvel.y/=sumwab;
    //  sumvel.z/=sumwab;
    //  PtVel=ToTFloat3(sumvel);
    //}
    //-Stores result. | Guarda resultado.
    ptvel=ToTFloat3(sumvel);
  }
  //-Stores result. | Guarda resultado.
  Result.Set(timestep,ToTFloat3(Point),ptvel);
  //Log->Printf("------> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}
//==============================================================================
/// Calculates velocity at indicated points (on CPU).
//==============================================================================
void JGaugeVelocity::CalculeCpu(const StDataCpu& datacpu){
  switch(CSP.tkernel){
    case KERNEL_Cubic:       //Kernel Wendland is used since Cubic is not available.
    case KERNEL_Wendland:    CalculeCpuT<KERNEL_Wendland>  (datacpu);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Check if GPU memory was allocated.
//==============================================================================
bool JGaugeVelocity::AllocatedGpuMemory(int id)const{
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  return(AuxDataGpu[id].GpuMemory);
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JGaugeVelocity::FreeGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeVelDataGpu& aug=AuxDataGpu[id];
  aug.GpuMemory=false;
  if(aug.Resultg)cudaFree(aug.Resultg); aug.Resultg=NULL;
}

//==============================================================================
/// Allocates GPU memory.
//==============================================================================
void JGaugeVelocity::AllocGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeVelDataGpu& aug=AuxDataGpu[id];
  if(aug.GpuMemory)FreeGpuMemory(id);
  fcuda::Malloc(&aug.Resultg,1);
  aug.GpuMemory=true; 
}

//==============================================================================
/// Calculates velocity at indicated points (on GPU).
//==============================================================================
void JGaugeVelocity::CalculeGpu(const StDataGpu& datagpu){
  const int id=0;
  const int indomain=2;
  const tdouble3 domposmin=DomPosMin,domposmax=DomPosMax;
  //-Run GPU calculation.
  if(indomain==1 || indomain==2){
    //-Prepare input data.
    const StDivDataGpu& dvd=datagpu.dvd;
    const double2*  posxy  =datagpu.posxy_g->cptr();
    const double*   posz   =datagpu.posz_g->cptr();
    const typecode* code   =datagpu.code_g->cptr();
    //const unsigned* idp  =datagpu.idp_g->cptr();
    const float4*   velrho =datagpu.velrho_g->cptr();
    //-Start measure.
    StGaugeVelDataGpu& aug=AuxDataGpu[id];
    if(!aug.GpuMemory)AllocGpuMemory(id);
    tfloat3 ptvel=TFloat3(0);
    //const bool ptout=PointIsOut(Point.x,Point.y,Point.z,domposmin,domposmax);
    //if(!ptout){//-Verify that the point is within domain limits.
    {
      cugauge::Interaction_GaugeVel(CSP,dvd,Point,posxy,posz,code,velrho,aug.Resultg);
      cudaMemcpy(&ptvel,aug.Resultg,sizeof(float3),cudaMemcpyDeviceToHost);
      Check_CudaErroor("Failed in velocity calculation.");
    }
    //-Stores calculated result.
    if(indomain==2){
      const double& timestep=datagpu.timestep;
      SetTimeStep(timestep);
      Result.Set(timestep,ToTFloat3(Point),ptvel);
      if(Output(timestep))StoreResult();
      //Log->Printf("------> t:%f",TimeStep);
    }
  }
}

#endif


//##############################################################################
//# JGaugeSwl
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeSwl::JGaugeSwl(unsigned idx,std::string name,tdouble3 point0
  ,tdouble3 point2,double pointdp,float masslimit,int gpucount)
  :JGaugeItem(GAUGE_Swl,idx,name,gpucount)
{
  ClassName="JGaugeSwl";
  FileInfo=string("Saves SWL data measured from fluid particles (by ")+ClassName+").";
  Reset();
  SetPoints(point0,point2,pointdp);
  MassLimit=masslimit;

}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeSwl::~JGaugeSwl(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeSwl::Reset(){
  PointDp=0;
  SetPoints(TDouble3(0),TDouble3(0),0);
  MassLimit=0;
  #ifdef _WITHGPU
    ResetGpuMemory();
  #endif
  JGaugeItem::Reset();
}

#ifdef _WITHGPU
//==============================================================================
/// Reset GPU memory.
//==============================================================================
void JGaugeSwl::ResetGpuMemory(){
  for(int g=0;g<GpuCount;g++){
    if(AuxDataGpu[g].GpuMemory)
      Run_Exceptioon(fun::PrintStr("Auxiliary GPU memory for unit %d is not free.",g));
    AuxDataGpu[g]=StrGaugeSwlDataGpu();
  }
}
#endif

//==============================================================================
/// Changes points definition.
//==============================================================================
void JGaugeSwl::SetPoints(const tdouble3& point0,const tdouble3& point2
  ,double pointdp)
{
  //if(PointDp<=0)Run_Exceptioon(fun::PrintStr("The value of PointDp is <= zero in gauge \'%s\'.",Name.c_str()));
  Point0=point0;
  Point2=point2;
  if(pointdp)PointDp=pointdp;
  const double dis=fgeo::PointsDist(Point0,Point2);
  if(dis>0 && PointDp>0){
    PointNp=unsigned(dis/PointDp);
    if(dis-(PointDp*PointNp)>=PointDp*0.1)PointNp++;
    if(PointNp<1)PointNp++;
    const double dp=dis/PointNp;
    //printf("------> PointNp:%d dp:%f\n",PointNp,dp);
    PointDir=fgeo::VecUnitary(Point2-Point0)*dp;
  }
  else{
    PointNp=0;
    PointDir=TDouble3(0);
  }
  //-Updates limits of calculation area.
  if(Scell){
    if(FixedDomMCel)Run_Exceptioon("Set points is invalid for fixed calculation domain.");
    ConfigDomMCel(FixedDomMCel);
  }
}

//==============================================================================
/// Configuration of initial limits of calculation area.
//==============================================================================
void JGaugeSwl::ConfigDomMCel(bool fixed){
  FixedDomMCel=fixed;
  const tdouble3 psmin=MinValues(Point0,Point2);
  const tdouble3 psmax=MaxValues(Point0,Point2);
  CalcMCelIniFinFromPos(psmin,psmax,DomMCelIni0,DomMCelFin0);
}

//==============================================================================
/// Record the last measure result.
//==============================================================================
void JGaugeSwl::StoreResult(){
  if(OutputSave){
    //-Allocates memory.
    while(unsigned(OutBuff.size())<OutSize)OutBuff.push_back(StrGaugeSwlRes());
    //-Empty buffer.
    if(OutCount+1>=OutSize)SaveResults();
    //-Stores last results.
    OutBuff[OutCount]=Result;
    OutCount++;
    //-Updates OutputNext.
    if(OutputDt){
      const unsigned nt=unsigned(TimeStep/OutputDt);
      OutputNext=OutputDt*nt;
      if(OutputNext<=TimeStep)OutputNext=OutputDt*(nt+1);
    }
  }
}

//==============================================================================
/// Saves stored results in CSV file.
//==============================================================================
void JGaugeSwl::SaveResults(){
  if(OutCount){
    const bool first=OutFile.empty();
    if(first){
      OutFile=GetResultsFileCsv();
      Log->AddFileInfo(OutFile,FileInfo);
    }
    jcsv::JSaveCsv2 scsv(OutFile,!first,AppInfo.GetCsvSepComa());
    //-Saves head.
    if(first){
      //-Head of values.
      scsv.SetHead();
      scsv <<"time [s];swlx [m];swly [m];swlz [m];pos0x [m];pos0y [m];pos0z [m];pos2x [m];pos2y [m];pos2z [m]" << jcsv::Endl();
    }
    //-Saves data.
    scsv.SetData();
    scsv << jcsv::Fmt(jcsv::TpFloat1,"%g") << jcsv::Fmt(jcsv::TpFloat3,"%g;%g;%g");
    for(unsigned c=0;c<OutCount;c++){
      scsv << OutBuff[c].timestep << OutBuff[c].posswl << OutBuff[c].point0 << OutBuff[c].point2 << jcsv::Endl();
    }
    OutCount=0;
  }
}

//==============================================================================
/// Saves last result in VTK file.
//==============================================================================
void JGaugeSwl::SaveVtkResult(unsigned cpart){
  JDataArrays arrays;
  arrays.AddArray("Pos",1,&(Result.posswl),false);
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  JSpVtkData::Save(fun::FileNameSec(GetResultsFileVtk(),cpart),arrays,"Pos");
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeSwl::GetPointDef(std::vector<tfloat3>& points)const{
  tdouble3 ps=Point0;
  for(unsigned p=0;p<=PointNp;p++){
    points.push_back(ToTFloat3(ps));
    ps=ps+PointDir;
  }
  return(PointNp+1);
}

//==============================================================================
/// Returns the interpolated mass value at the indicated point. The point must 
/// belong to the cell domain.
///
/// Devuelve valor de masa interpolado en el punto indicado. El punto debe 
/// pertenecer al dominio de celdas.
//==============================================================================
template<TpKernel tker> float JGaugeSwl::CalculeMassCpu(const tdouble3& ptpos
  ,const StDivDataCpu& dvd,const tdouble3* pos,const typecode* code
  ,const tfloat4* velrho)const
{
  //-Auxiliary variables.
  double sumwab=0;
  double summass=0;
  //-Search for fluid neighbours in adjacent cells.
  const StNgSearch ngs=nsearch::Init(ptpos,false,dvd);
  for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
    const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
    for(unsigned p2=pif.x;p2<pif.y;p2++){
      const float rr2=nsearch::Distance2(ptpos,pos[p2]);
      //-Interaction with real neighbouring particles.
      if(rr2<=CSP.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
        float wab=fsph::GetKernel_Wab<tker>(CSP,rr2);
        wab*=CSP.massfluid/velrho[p2].w;
        sumwab+=wab;
        summass+=wab*CSP.massfluid;
      }
    }
  }
  //-Applies kernel correction.
  //if(sumwab!=0)summass/=sumwab;
  return(float(summass));
}

//==============================================================================
/// Calculates surface water level at indicated points (on CPU).
//==============================================================================
template<TpKernel tker> void JGaugeSwl::CalculeCpuT(const StDataCpu& datacpu)
{
  //-Prepare input data.
  const double& timestep =datacpu.timestep;
  const StDivDataCpu& dvd=datacpu.dvd;
  const tdouble3* pos    =datacpu.pos_c->cptr();
  const typecode* code   =datacpu.code_c->cptr();
  //const unsigned* idp  =datacpu.idp_c->cptr();
  const tfloat4*  velrho =datacpu.velrho_c->cptr();
  SetTimeStep(timestep);
  //-Look for change of fluid to empty. | Busca paso de fluido a vacio.
  tdouble3 ptsurf=TDouble3(DBL_MAX);
  float mpre=0;
  tdouble3 ptpos=Point0;
  for(unsigned cp=0;cp<=PointNp;cp++){
    const float mass=CalculeMassCpu<tker>(ptpos,dvd,pos,code,velrho);
    if(mass>MassLimit)mpre=mass;
    if(mass<MassLimit && mpre){
      const float fxm1=(MassLimit-mpre)/(mass-mpre)-1;
      ptsurf=ptpos+(PointDir*double(fxm1));
      cp=PointNp+1;
    }
    ptpos=ptpos+PointDir;
  }
  if(ptsurf.x==DBL_MAX)ptsurf=Point0+(PointDir*(mpre? PointNp: 0));
  //-Stores result. | Guarda resultado.
  Result.Set(timestep,ToTFloat3(Point0),ToTFloat3(Point2),ToTFloat3(ptsurf));
  //Log->Printf("---JGaugeSwl::CalculeCpuT---> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}

//==============================================================================
/// Calculates surface water level at indicated points (on CPU).
//==============================================================================
void JGaugeSwl::CalculeCpu(const StDataCpu& datacpu){
  switch(CSP.tkernel){
    case KERNEL_Cubic:       //Kernel Wendland is used since Cubic is not available.
    case KERNEL_Wendland:    CalculeCpuT<KERNEL_Wendland>  (datacpu);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Check if GPU memory was allocated.
//==============================================================================
bool JGaugeSwl::AllocatedGpuMemory(int id)const{
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  return(AuxDataGpu[id].GpuMemory);
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JGaugeSwl::FreeGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeSwlDataGpu& aug=AuxDataGpu[id];
  aug.GpuMemory=false;
  if(aug.Resultg)cudaFree(aug.Resultg); aug.Resultg=NULL;
}

//==============================================================================
/// Allocates GPU memory.
//==============================================================================
void JGaugeSwl::AllocGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeSwlDataGpu& aug=AuxDataGpu[id];
  if(aug.GpuMemory)FreeGpuMemory(id);
  fcuda::Malloc(&aug.Resultg,1);
  aug.GpuMemory=true; 
}

//==============================================================================
/// Calculates surface water level at indicated points (on GPU).
//==============================================================================
void JGaugeSwl::CalculeGpu(const StDataGpu& datagpu){
  const int id=0;
  const int indomain=2;
  //-Run GPU calculation.
  if(indomain==1 || indomain==2){
    //-Prepare input data.
    const StDivDataGpu& dvd=datagpu.dvd;
    const double2*  posxy  =datagpu.posxy_g->cptr();
    const double*   posz   =datagpu.posz_g->cptr();
    const typecode* code   =datagpu.code_g->cptr();
    //const unsigned* idp  =datagpu.idp_g->cptr();
    const float4*   velrho =datagpu.velrho_g->cptr();
    //-Start measure.
    StGaugeSwlDataGpu& aug=AuxDataGpu[id];
    if(!aug.GpuMemory)AllocGpuMemory(id);
    if(indomain==2){
      cugauge::Interaction_GaugeSwl(CSP,dvd,Point0,PointDir,PointNp,MassLimit
        ,posxy,posz,code,velrho,aug.Resultg);
      tfloat3 ptsurf=TFloat3(0);
      cudaMemcpy(&ptsurf,aug.Resultg,sizeof(float3),cudaMemcpyDeviceToHost);
      Check_CudaErroor("Failed in Swl calculation.");
      //-Stores calculated result.
      const double& timestep=datagpu.timestep;
      SetTimeStep(timestep);
      Result.Set(timestep,ToTFloat3(Point0),ToTFloat3(Point2),ptsurf);
      if(Output(timestep))StoreResult();
      //Log->Printf("------> t:%f",TimeStep);
    }
  }
}

#endif


//##############################################################################
//# JGaugeMaxZ
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeMaxZ::JGaugeMaxZ(unsigned idx,std::string name,tdouble3 point0
  ,double height,float distlimit,int gpucount)
  :JGaugeItem(GAUGE_MaxZ,idx,name,gpucount)
{
  ClassName="JGaugeMaxZ";
  FileInfo=string("Saves maximum Z position of fluid particles (by ")+ClassName+").";
  Reset();
  SetPoint0(point0);
  SetHeight(height);
  SetDistLimit(distlimit);
}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeMaxZ::~JGaugeMaxZ(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeMaxZ::Reset(){
  SetPoint0(TDouble3(0));
  SetHeight(0);
  SetDistLimit(0);
  #ifdef _WITHGPU
    ResetGpuMemory();
  #endif
  JGaugeItem::Reset();
}

#ifdef _WITHGPU
//==============================================================================
/// Reset GPU memory.
//==============================================================================
void JGaugeMaxZ::ResetGpuMemory(){
  for(int g=0;g<GpuCount;g++){
    if(AuxDataGpu[g].GpuMemory)
      Run_Exceptioon(fun::PrintStr("Auxiliary GPU memory for unit %d is not free.",g));
    AuxDataGpu[g]=StrGaugeMaxzDataGpu();
  }
}
#endif

//==============================================================================
/// Configuration of initial limits of calculation area.
//==============================================================================
void JGaugeMaxZ::ConfigDomMCel(bool fixed){
  if(!fixed)Run_Exceptioon("Unfixed calculation domain is invalid for MaxZ-Gauge.");
  FixedDomMCel=fixed;
  const unsigned ncel=unsigned(ceil(DistLimit/Scell));
  tuint3 mcel0=CalcMCelFromPos(Point0);
  const unsigned cxi=(mcel0.x>ncel? mcel0.x-ncel: 0);
  const unsigned cyi=(mcel0.y>ncel? mcel0.y-ncel: 0);
  const unsigned czi=mcel0.z;
  const unsigned cxf=mcel0.x+(ncel+1);
  const unsigned cyf=mcel0.y+(ncel+1);
  const double   dz =(Point0.z+Height)-MapPosMin.z;
  const unsigned czf=(dz>=0? unsigned(dz/Scell): 0)+1;
  DomMCelIni0=TUint3(cxi,cyi,czi);
  DomMCelFin0=TUint3(cxf,cyf,czf);
}

//==============================================================================
/// Record the last measure result.
//==============================================================================
void JGaugeMaxZ::StoreResult(){
  if(OutputSave){
    //-Allocates memory.
    while(unsigned(OutBuff.size())<OutSize)OutBuff.push_back(StrGaugeMaxzRes());
    //-Empty buffer.
    if(OutCount+1>=OutSize)SaveResults();
    //-Stores last results.
    OutBuff[OutCount]=Result;
    OutCount++;
    //-Updates OutputNext.
    if(OutputDt){
      const unsigned nt=unsigned(TimeStep/OutputDt);
      OutputNext=OutputDt*nt;
      if(OutputNext<=TimeStep)OutputNext=OutputDt*(nt+1);
    }
  }
}

//==============================================================================
/// Saves stored results in CSV file.
//==============================================================================
void JGaugeMaxZ::SaveResults(){
  if(OutCount){
    const bool first=OutFile.empty();
    if(first){
      OutFile=GetResultsFileCsv();
      Log->AddFileInfo(OutFile,FileInfo);
    }
    jcsv::JSaveCsv2 scsv(OutFile,!first,AppInfo.GetCsvSepComa());
    //-Saves head.
    if(first){
      scsv.SetHead();
      scsv << "time [s];zmax [m];posx [m];posy [m];posz [m]" << jcsv::Endl();
    }
    //-Saves data.
    scsv.SetData();
    scsv << jcsv::Fmt(jcsv::TpFloat1,"%g") << jcsv::Fmt(jcsv::TpFloat3,"%g;%g;%g");
    for(unsigned c=0;c<OutCount;c++){
      scsv << OutBuff[c].timestep << OutBuff[c].zmax << OutBuff[c].point0 << jcsv::Endl();
    }
    OutCount=0;
  }
}

//==============================================================================
/// Saves last result in VTK file.
//==============================================================================
void JGaugeMaxZ::SaveVtkResult(unsigned cpart){
  //-Prepares data.
  const tfloat3 pt0=Result.point0;
  const tfloat3 ptz=TFloat3(pt0.x,pt0.y,Result.zmax);
  const float height=ptz.z-pt0.z;
  //Log->Printf("---->ptz:(%g,%g,%g)  h:%g",ptz.x,ptz.y,ptz.z,height);
  JDataArrays arrays;
  arrays.AddArray("Pos",1,&ptz,false);
  arrays.AddArray("Height",1,&height,false);
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  JSpVtkData::Save(fun::FileNameSec(GetResultsFileVtk(),cpart),arrays,"Pos");
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeMaxZ::GetPointDef(std::vector<tfloat3>& points)const{
  const float dp0=Scell/2;
  unsigned np=unsigned(Height/dp0);
  if(Height-(dp0*np)>=dp0*0.1)np++;
  if(np<1)np++;
  const double dp=Height/np;
  const tdouble3 dir=TDouble3(0,0,dp);
  tdouble3 ps=Point0;
  for(unsigned p=0;p<=np;p++){
    points.push_back(ToTFloat3(ps));
    ps=ps+dir;
  }
  return(np+1);
}

//==============================================================================
/// Return cell limits for interaction starting from position.
/// Devuelve limites de celdas para interaccion a partir de posicion.
//==============================================================================
void JGaugeMaxZ::GetInteractionCellsMaxZ(const tdouble3& pos,const tint4& nc
  ,const tdouble3& domposmin,const tint3& cellzero
  ,int& cxini,int& cxfin,int& yini,int& yfin,int& zini,int& zfin)const
{
  //-Get cell coordinates of position pos.
  const int cx=int((pos.x-domposmin.x)/Scell)-cellzero.x;
  const int cy=int((pos.y-domposmin.y)/Scell)-cellzero.y;
  const int cz=int((pos.z-domposmin.z)/Scell)-cellzero.z;
  const int cz2=int((pos.z+Height-domposmin.z)/Scell)-cellzero.z;
  //-Calculates range of cells to check.
  const int ncel=int(ceil(DistLimit/Scell));
  cxini=max(cx-ncel,0);
  cxfin=min(cx+ncel+1,nc.x);
  yini=max(cy-ncel,0);
  yfin=min(cy+ncel+1,nc.y);
  zini=max(cz-ncel,0);
  zfin=min(cz2+ncel+1,nc.z);
}

//==============================================================================
/// Calculates maximum z of fluid at distance of a vertical line (on CPU).
//==============================================================================
void JGaugeMaxZ::CalculeCpu(const StDataCpu& datacpu){
  //Log->Printf("JGaugeMaxZ----> timestep:%g  (%d)",timestep,(DG?1:0));
  //-Prepare input data.
  const double& timestep =datacpu.timestep;
  const StDivDataCpu& dvd=datacpu.dvd;
  const tdouble3* pos    =datacpu.pos_c->cptr();
  const typecode* code   =datacpu.code_c->cptr();
  //const unsigned* idp  =datacpu.idp_c->cptr();
  //const tfloat4*  velrho=datacpu.velrho_c->cptr();
  SetTimeStep(timestep);
  //-Compute auxiliary constants.
  const float maxdist2=DistLimit*DistLimit;
  //-Obtain limits of interaction.
  int cxini,cxfin,yini,yfin,zini,zfin;
  GetInteractionCellsMaxZ(Point0,dvd.nc,dvd.domposmin,dvd.cellzero
    ,cxini,cxfin,yini,yfin,zini,zfin);
  //-Start measure.
  unsigned pmax=UINT_MAX;
  float zmax=-FLT_MAX;
  //-Search for neighbours in adjacent cells. | Busqueda de vecinos en celdas adyacentes.
  if(cxini<cxfin)for(int z=zfin-1;z>=zini && pmax==UINT_MAX;z--){
    const int zmod=(dvd.nc.w)*z+dvd.cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
    for(int y=yini;y<yfin;y++){
      int ymod=zmod+dvd.nc.x*y;
      const unsigned pini=dvd.begincell[cxini+ymod];
      const unsigned pfin=dvd.begincell[cxfin+ymod];

      //-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
      //---------------------------------------------------------------------------------------------
      for(unsigned p2=pini;p2<pfin;p2++){
        if(pos[p2].z>zmax){
          const float drx=float(Point0.x-pos[p2].x);
          const float dry=float(Point0.y-pos[p2].y);
          const float rr2=drx*drx+dry*dry;
          if(rr2<=maxdist2 && CODE_IsFluid(code[p2])){//-Only with fluid particles.
            zmax=float(pos[p2].z);
            pmax=p2;
          }
        }
      }
    }
  }
  //-Stores result. | Guarda resultado.
  Result.Set(timestep,ToTFloat3(Point0),max(zmax,float(DomPosMin.z)));
  //Log->Printf("JGaugeMaxZ> zmax:%g",Result.zmax);
  if(Output(timestep))StoreResult();
}

#ifdef _WITHGPU

//==============================================================================
/// Check if GPU memory was allocated.
//==============================================================================
bool JGaugeMaxZ::AllocatedGpuMemory(int id)const{
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  return(AuxDataGpu[id].GpuMemory);
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JGaugeMaxZ::FreeGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeMaxzDataGpu& aug=AuxDataGpu[id];
  aug.GpuMemory=false;
  if(aug.Resultg)cudaFree(aug.Resultg); aug.Resultg=NULL;
}

//==============================================================================
/// Allocates GPU memory.
//==============================================================================
void JGaugeMaxZ::AllocGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeMaxzDataGpu& aug=AuxDataGpu[id];
  if(aug.GpuMemory)FreeGpuMemory(id);
  fcuda::Malloc(&aug.Resultg,1);
  aug.GpuMemory=true; 
}

//==============================================================================
/// Calculates maximum z of fluid at distance of a vertical line (on GPU).
//==============================================================================
void JGaugeMaxZ::CalculeGpu(const StDataGpu& datagpu){
  const int id=0;
  const int indomain=2;
  //-Run GPU calculation.
  if(indomain==1 || indomain==2){
    //-Prepare input data.
    const StDivDataGpu& dvd=datagpu.dvd;
    const double2*  posxy  =datagpu.posxy_g->cptr();
    const double*   posz   =datagpu.posz_g->cptr();
    const typecode* code   =datagpu.code_g->cptr();
    //const unsigned* idp    =datagpu.idp_g->cptr();
    //const float4*   velrho =datagpu.velrho_g->cptr();
    //-Start measure.
    StGaugeMaxzDataGpu& aug=AuxDataGpu[id];
    if(!aug.GpuMemory)AllocGpuMemory(id);
    tfloat3 ptsurf=TFloat3(0);
    //-Compute auxiliary constants.
    const float maxdist2=DistLimit*DistLimit;
    //-Obtain limits of interaction.
    const tint4 nc=TInt4(dvd.nc.x,dvd.nc.y,dvd.nc.z,dvd.nc.w);
    const tdouble3 domposmin=TDouble3(dvd.domposmin.x,dvd.domposmin.y,dvd.domposmin.z);
    const tint3 cellzero=TInt3(dvd.cellzero.x,dvd.cellzero.y,dvd.cellzero.z);
    int cxini,cxfin,yini,yfin,zini,zfin;
    {
      GetInteractionCellsMaxZ(Point0,nc,domposmin,cellzero
        ,cxini,cxfin,yini,yfin,zini,zfin);
      //-Start measure.
      cugauge::Interaction_GaugeMaxz(Point0,maxdist2,dvd
        ,cxini,cxfin,yini,yfin,zini,zfin,posxy,posz,code,aug.Resultg);
    }
    cudaMemcpy(&ptsurf,aug.Resultg,sizeof(float3),cudaMemcpyDeviceToHost);
    Check_CudaErroor("Failed in MaxZ calculation.");
    //-Stores calculated result.
    if(indomain==2){
      const double& timestep=datagpu.timestep;
      SetTimeStep(timestep);
      Result.Set(timestep,ToTFloat3(Point0),ptsurf.z);
      if(Output(timestep))StoreResult();
      //Log->Printf("------> t:%f",TimeStep);
    }
  }
}

#endif


//<vs_meeshdat_ini>
//##############################################################################
//# JGaugeMesh
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeMesh::JGaugeMesh(unsigned idx,std::string name,const jmsh::StMeshBasic& meshbas
  ,std::string outdata,unsigned tfmt,unsigned buffersize
  ,float kclimit,float kcdummy,float masslimit,int gpucount)
  :JGaugeItem(GAUGE_Mesh,idx,name,gpucount,(buffersize<1? 30: buffersize))
{
  ClassName="JGaugeMesh";
  FileInfo=string("Saves mesh data measured from fluid particles (by ")+ClassName+").";
  MeshDat=NULL;
  MassDatCpu=NULL;
  MeshDataSave=NULL; 
  Reset();
  SetMeshData(meshbas,outdata);
  KcLimit=kclimit;
  KcDummy=kcdummy;
  MassLimit=masslimit;
  SaveBin=(tfmt&jmsh::TpFmtBin)!=0;
  SaveCsv=(tfmt&jmsh::TpFmtCsv)!=0;
}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeMesh::~JGaugeMesh(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeMesh::Reset(){
  OutDataList="";
  ComputeRhop=ComputeVelxyz=ComputeVeldir=ComputeZsurf=false;
  SaveCsv=SaveBin=false;
  memset(&MeshBas,0,sizeof(jmsh::StMeshBasic));
  memset(&MeshPts,0,sizeof(jmsh::StMeshPts));
  KcLimit=KcDummy=0;
  MassLimit=0;
  delete   MeshDat;    MeshDat=NULL;
  delete[] MassDatCpu; MassDatCpu=NULL;
  #ifdef _WITHGPU
    ResetGpuMemory();
  #endif
  Result.Reset();
  for(unsigned c=0;c<unsigned(OutBuff.size());c++)delete OutBuff[c].meshdat;
  OutBuff.clear();
  delete MeshDataSave; MeshDataSave=NULL; 
  JGaugeItem::Reset();
}

#ifdef _WITHGPU
//==============================================================================
/// Reset GPU memory.
//==============================================================================
void JGaugeMesh::ResetGpuMemory(){
  for(int g=0;g<GpuCount;g++){
    if(AuxDataGpu[g].GpuMemory)
      Run_Exceptioon(fun::PrintStr("Auxiliary GPU memory for unit %d is not free.",g));
    AuxDataGpu[g]=StrGaugeMeshDataGpu();
  }
}
#endif

//==============================================================================
/// Configures mesh positions for calculation.
//==============================================================================
void JGaugeMesh::SetMeshData(const jmsh::StMeshBasic& meshbas,std::string outdata){
  //-Free memory.
  delete   MeshDat; MeshDat=NULL;
  delete[] MassDatCpu; MassDatCpu=NULL;
  //-Configures mesh definition.
  MeshBas=meshbas;
  //printf("DDD_000 v3:(%f,%f,%f) dis3:%f/%f\n",MeshBas.vec3.x,MeshBas.vec3.y,MeshBas.vec3.z,MeshBas.dis3,MeshBas.dispt3);
  MeshPts=jmsh::JMeshData::MakeMeshPt(meshbas);
  //jmsh::JMeshData::PrintMeshPts(MeshPts);
  //-Configures output data.
  ConfigDataList(outdata);
  //-Allocates memory on CPU.
  MeshDat=new jmsh::JMeshData();
  MeshDat->ConfigMesh(MeshPts,0,OutDataList);
  if(Cpu && ComputeZsurf){
    try{
      MassDatCpu=new float[MeshPts.npt];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Returns corrected datalist or "error" when it is invalid.
//==============================================================================
std::string JGaugeMesh::CorrectDataList(std::string datalist){
  bool err=false;
  string ret;
  datalist=fun::StrLower(fun::StrWithoutChar(datalist,' '));
  vector<string> vstr;
  unsigned nstr=fun::VectorSplitStr(",",datalist,vstr);
  for(unsigned c=0;c<nstr && !err;c++){
         if(vstr[c]=="vel"   )ret=ret+(ret.empty()? "": ",")+"Vel";
    else if(vstr[c]=="veldir")ret=ret+(ret.empty()? "": ",")+"VelDir";
    else if(vstr[c]=="rhop"  )ret=ret+(ret.empty()? "": ",")+"Rhop";
    else if(vstr[c]=="zsurf" )ret=ret+(ret.empty()? "": ",")+"Zsurf";
    else err=true;
  }
  return(err? string("error"): ret);
}

//==============================================================================
/// Configure the output data selection.
//==============================================================================
void JGaugeMesh::ConfigDataList(std::string datalist){
  ComputeRhop=ComputeVelxyz=ComputeVeldir=ComputeZsurf=false;
  datalist=fun::StrLower(fun::StrWithoutChar(datalist,' '));
  vector<string> vstr;
  unsigned nstr=fun::VectorSplitStr(",",datalist,vstr);
  for(unsigned c=0;c<nstr;c++){
         if(vstr[c]=="vel"   )ComputeVelxyz=true;
    else if(vstr[c]=="veldir")ComputeVeldir=true;
    else if(vstr[c]=="rhop"  )ComputeRhop=true;
    else if(vstr[c]=="zsurf" )ComputeZsurf=true;
    else Run_Exceptioon("Output data name is unknown.");
  }
  OutDataList=GetDataList();
}

//==============================================================================
/// Returns the output data selection as string.
//==============================================================================
std::string JGaugeMesh::GetDataList()const{
  string ret;
  if(ComputeVelxyz)ret=ret+(ret.empty()? "": ",")+"Vel";
  if(ComputeVeldir)ret=ret+(ret.empty()? "": ",")+"VelDir";
  if(ComputeRhop  )ret=ret+(ret.empty()? "": ",")+"Rhop";
  if(ComputeZsurf )ret=ret+(ret.empty()? "": ",")+"Zsurf";
  return(ret);
}

//==============================================================================
/// Returns structure with basic information about configuration.
//==============================================================================
JGaugeMesh::StInfo JGaugeMesh::GetInfo()const{
  StInfo v={TDouble3(0),TDouble3(0),TDouble3(0),TDouble3(0),TDouble3(0),TDouble3(0),TUint4(0),TFloat3(0),""};
  const jmsh::StMeshBasic& mb=MeshBas;
  const jmsh::StMeshPts& m=MeshPts;
  v.ptref=m.ptref;
  v.ptend=m.ptref+(m.vdp1*m.npt1)+(m.vdp2*m.npt2)+(m.vdp3*m.npt3);
  v.vec1=mb.vec1;
  v.vec2=mb.vec2;
  v.vec3=mb.vec3;
  v.dispt=TDouble3(mb.dispt1,mb.dispt2,mb.dispt3);
  v.npt=TUint4(m.npt1,m.npt2,m.npt3,m.npt);
  v.dirdat=m.dirdat;
  v.outdata=OutDataList;
  return(v);
}

//==============================================================================
/// Configuration of initial limits of calculation area.
//==============================================================================
void JGaugeMesh::ConfigDomMCel(bool fixed){
  if(!fixed)Run_Exceptioon("Unfixed calculation domain is invalid for Mesh-Gauge.");
  FixedDomMCel=fixed;
  const jmsh::StMeshPts& m=MeshPts;
  const tdouble3 ptref=m.ptref;
  const tdouble3 ptend=m.ptref+(m.vdp1*m.npt1)+(m.vdp2*m.npt2)+(m.vdp3*m.npt3);
  const tdouble3 psmin=MinValues(ptref,ptend);
  const tdouble3 psmax=MaxValues(ptref,ptend);
  CalcMCelIniFinFromPos(psmin,psmax,DomMCelIni0,DomMCelFin0);
}

//==============================================================================
/// Record the last measure result.
//==============================================================================
void JGaugeMesh::StoreResult(){
  if(OutputSave){
    //-Allocates memory for OutSize records.
    while(unsigned(OutBuff.size())<OutSize){
      StMeshRes res=StrMeshRes();
      res.meshdat=new jmsh::JMeshData();
      OutBuff.push_back(res);
    }
    //-Empty buffer.
    if(OutCount+1>=OutSize)SaveResults();
    //-Stores last results.
    jmsh::JMeshData* meshdat=OutBuff[OutCount].meshdat;
    if(!meshdat->GetNpt()){
      meshdat->ConfigMesh(MeshPts,0,OutDataList);
    }
    meshdat->CopyDataFrom(Result.meshdat);
    OutBuff[OutCount].Set(Result.timestep,meshdat);
    OutCount++;
    //-Updates OutputNext.
    if(OutputDt){
      const unsigned nt=unsigned(TimeStep/OutputDt);
      OutputNext=OutputDt*nt;
      if(OutputNext<=TimeStep)OutputNext=OutputDt*(nt+1);
    }
  }
}

//==============================================================================
/// Saves stored results in CSV file.
//==============================================================================
void JGaugeMesh::SaveResults(){
  if(OutCount){
    if(SaveCsv){
      bool first=OutFile.empty();
      if(first){
        OutFile=GetResultsFileCsv("XXX");
        Log->AddFileInfo(OutFile,FileInfo);
      }
      const JDataArrays* arr=OutBuff[0].meshdat->GetArrays();
      const unsigned na=arr->Count();
      for(unsigned ca=0;ca<na;ca++){
        const JDataArrays::StDataArray& ar=arr->GetArrayCte(ca);
        string file=GetResultsFileCsv(string("_")+ar.keyname);
        jcsv::JSaveCsv2 scsv(file,!first,AppInfo.GetCsvSepComa());
        for(unsigned c=0;c<OutCount;c++){
          jmsh::JMeshTDatasSave::SaveCsv(OutBuff[c].meshdat,ca,scsv,first && !c);
        }
      }
    }
    if(SaveBin){
      if(!MeshDataSave){
        const string file=GetResultsFile(false,"mbi4");
        Log->AddFileInfo(file,FileInfo);
        MeshDataSave=new jmsh::JMeshTDatasSave();
        MeshDataSave->Config(file,AppInfo.GetFullName(),OutBuff[0].meshdat);
      }
      if(OutCount>1){//-MultiData version.
        //printf("---> OutCount:%u\n",OutCount);
        jmsh::JMeshData** vmeshdat=new jmsh::JMeshData*[OutCount];
        for(unsigned c=0;c<OutCount;c++)vmeshdat[c]=OutBuff[c].meshdat;
        MeshDataSave->SaveDataTimes(OutCount,vmeshdat);
        delete[] vmeshdat; vmeshdat=NULL;
      }
      else{//-SingleData version.
        for(unsigned c=0;c<OutCount;c++){
          MeshDataSave->SaveDataTime(OutBuff[c].meshdat);
        }
      }
    }
    OutCount=0;
  }
}

//==============================================================================
/// Saves last result in VTK file.
//==============================================================================
void JGaugeMesh::SaveVtkResult(unsigned cpart){
  //-Saves VTK with npt size data.
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  jmsh::JMeshTDatasSave::SaveVtk(GetResultsFileVtk(),int(cpart),MeshDat,false);
  //-Saves VTK with zsurf data.
  if(ComputeZsurf){
    Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk("_Zsurf"),UINT_MAX),FileInfo);
    jmsh::JMeshTDatasSave::SaveVtk(GetResultsFileVtk("_Zsurf"),int(cpart),MeshDat,true);
  }
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeMesh::GetPointDef(std::vector<tfloat3>& points)const{
  //printf("=> GetPointDef  MeshPts.npt:%u \n",MeshPts.npt);
  const unsigned np=MeshPts.npt;
  tfloat3* pos=new tfloat3[np];
  MeshDat->GetPosf(np,pos);
  for(unsigned p=0;p<np;p++)points.push_back(pos[p]);
  delete[] pos; pos=NULL;
  return(np);
}

//==============================================================================
/// Creates VTK file with the scheme of gauge configuration.
//==============================================================================
void JGaugeMesh::SaveVtkScheme()const{
  const string filevtk=AppInfo.GetDirOut()+"CfgGaugeMesh_"+Name+".vtk";
  Log->AddFileInfo(filevtk,"Saves VTK file with scheme of Mesh gauge.");
  jmsh::JMeshTDatasSave::SaveVtkScheme(filevtk,MeshPts);
}

//==============================================================================
/// Returns internal pointer to save Zsurf results on CPU.
//==============================================================================
const float* JGaugeMesh::GetPtrDataZsurf()const{
  return(MeshDat->GetVarFloat("Zsurf"));
}

#ifdef _WITHGPU
//==============================================================================
/// Returns internal pointer to save Zsurf results on GPU.
//==============================================================================
const float* JGaugeMesh::GetPtrDataZsurfg(int id)const{
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  return(AuxDataGpu[id].DataZsurfg);
}
#endif

//==============================================================================
/// Calculates velocity and swl at indicated points (on CPU).
//==============================================================================
template<TpKernel tker> void JGaugeMesh::CalculeCpuT(const StDataCpu& datacpu){
  //-Prepare input data.
  const double& timestep =datacpu.timestep;
  const StDivDataCpu& dvd=datacpu.dvd;
  const tdouble3* pos    =datacpu.pos_c->cptr();
  const typecode* code   =datacpu.code_c->cptr();
  //const unsigned* idp  =datacpu.idp_c->cptr();
  const tfloat4*  velrho =datacpu.velrho_c->cptr();
  SetTimeStep(timestep);
  //-Start measure.
  //MeshDat->ClearData();  //-It is not necessary.
  tfloat3* datavxyz =MeshDat->GetVarFloat3("Vel");
  float*   datavdir =MeshDat->GetVarFloat("VelDir");
  float*   datarhop =MeshDat->GetVarFloat("Rhop");
  float*   datazsurf=MeshDat->GetVarFloat("Zsurf");
  //-Obtains basic data for calculations.
  const jmsh::StMeshPts& m=MeshPts;
  const bool vel=(ComputeVelxyz || ComputeVeldir);
  const unsigned npt1=m.npt1;
  const unsigned npt2=m.npt2;
  const unsigned npt3=m.npt3;
  const unsigned npt12=npt1*npt2;
  const tdouble3 ptref=m.ptref;
  const tdouble3 vdp1=m.vdp1;
  const tdouble3 vdp2=m.vdp2;
  const tdouble3 vdp3=m.vdp3;
  const tfloat3 vdir=m.dirdat;
  //-Computes velocity and mass.
  const int npt=int(m.npt);
  for(int cp=0;cp<npt;cp++){
    const unsigned cp3=unsigned(cp)/npt12;
    const unsigned cp3r=unsigned(cp)-cp3*npt12;
    const unsigned cp2=unsigned(cp3r)/npt1;
    const unsigned cp1=unsigned(cp3r)-cp2*npt1;
    //const tdouble3 pt2=ptref+(vdp2*cp2);
    const tdouble3 ptpos=ptref+(vdp1*cp1)+(vdp2*cp2)+(vdp3*cp3);
    float sumwab=0;
    float summass=0;
    float sumrhop=0;
    tfloat3 sumvel=TFloat3(0);
    //-Search for fluid neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(ptpos,false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float rr2=nsearch::Distance2(ptpos,pos[p2]);
        //-Interaction with real neighbouring particles.
        if(rr2<=CSP.kernelsize2 && CODE_IsFluid(code[p2])){
          float wab=fsph::GetKernel_Wab<tker>(CSP,rr2);
          const tfloat4 velrho2=velrho[p2];
          wab*=CSP.massfluid/velrho2.w;
          sumwab+=wab;
          summass+=wab*CSP.massfluid;
          sumrhop+=wab*velrho2.w;
          if(vel){
            sumvel.x+=wab*velrho2.x;
            sumvel.y+=wab*velrho2.y;
            sumvel.z+=wab*velrho2.z;
          }
        }
      }
    }
    //-Applies kernel correction.
    if(KcLimit!=FLT_MAX){
      if(sumwab>=KcLimit){
        sumvel=sumvel/sumwab;
        sumrhop/=sumwab;
      }
      else if(KcDummy!=FLT_MAX){
        sumvel=TFloat3(KcDummy);
        sumrhop=KcDummy;
      }
    }
    //-Stores results.
    if(vel){
      if(datavxyz)datavxyz[cp]=sumvel;
      if(datavdir)datavdir[cp]=(sumvel.x*vdir.x + sumvel.y*vdir.y + sumvel.z*vdir.z);
    }
    if(datarhop)datarhop[cp]=sumrhop;
    if(MassDatCpu)MassDatCpu[cp]=summass;
  }
  //-Computes surface water level.
  if(ComputeZsurf){
    for(unsigned cp2=0;cp2<npt2;cp2++)for(unsigned cp1=0;cp1<npt1;cp1++){
      float masspre=0;
      unsigned cpsurf=0;
      float    fsurf=0;
      unsigned cp=cp1+cp2*npt1;
      //const bool DG=(cp1==1 && cp2==0);
      for(unsigned cp3=0;cp3<npt3 && !cpsurf;cp3++){
        const float mass=MassDatCpu[cp];
        //if(DG)Log->Printf("==> mass[%u,%u,%u=%u]:%f  Masslimit:%f  z:%f",cp1,cp2,cp3,cp,mass,MassLimit,   ptref.z+vdp3.z*cp3);
        if(mass>MassLimit)masspre=mass;
        if(mass<MassLimit && masspre){
          //if(DG)Log->Printf("==> mass<ML %f<%f  masspre:%f  cp3:%u  z:%f",mass,MassLimit,masspre,cp3,   ptref.z+vdp3.z*cp3);
          fsurf=(MassLimit-masspre)/(mass-masspre);
          cpsurf=cp3;
        }
        cp+=npt12;
      }
      float zsurf=float(ptref.z+(vdp1.z*cp1)+(vdp2.z*cp2));  //-Minimum zsurf.
      if(cpsurf==0 && masspre)zsurf+=float(vdp3.z*(npt3-1)); //-Maximum zsurf.
      if(cpsurf){
        zsurf+=float((vdp3.z*(cpsurf-1))+(vdp3.z*fsurf));    //-Found zsurf.
      }
      datazsurf[cp1+cp2*npt1]=zsurf;
      //if(DG)Log->Printf("==> zsurf[%u,%u]:%f  cpsurf:%d",cp1,cp2,zsurf,cpsurf);
    }
  }
  //-Stores result. | Guarda resultado.
  MeshDat->SetTimeStep(timestep);
  Result.Set(timestep,MeshDat);
  //Log->Printf("------> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}

//==============================================================================
/// Calculates data at indicated mesh points (on CPU).
//==============================================================================
void JGaugeMesh::CalculeCpu(const StDataCpu& datacpu){
  switch(CSP.tkernel){
    case KERNEL_Cubic:       //Kernel Wendland is used since Cubic is not available.
    case KERNEL_Wendland:    CalculeCpuT<KERNEL_Wendland>  (datacpu);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Check if GPU memory was allocated.
//==============================================================================
bool JGaugeMesh::AllocatedGpuMemory(int id)const{
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  return(AuxDataGpu[id].GpuMemory);
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JGaugeMesh::FreeGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeMeshDataGpu& aug=AuxDataGpu[id];
  aug.GpuMemory=false;
  //if(aug.Resultg   )cudaFree(aug.Resultg   ); aug.Resultg=NULL;
  if(aug.DataRhopg )cudaFree(aug.DataRhopg ); aug.DataRhopg=NULL;
  if(aug.DataVxyzg )cudaFree(aug.DataVxyzg ); aug.DataVxyzg=NULL;
  if(aug.DataVdirg )cudaFree(aug.DataVdirg ); aug.DataVdirg=NULL;
  if(aug.DataZsurfg)cudaFree(aug.DataZsurfg); aug.DataZsurfg=NULL;
  if(aug.DataMassg )cudaFree(aug.DataMassg ); aug.DataMassg=NULL;
}

//==============================================================================
/// Allocates GPU memory.
//==============================================================================
void JGaugeMesh::AllocGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeMeshDataGpu& aug=AuxDataGpu[id];
  if(aug.GpuMemory)FreeGpuMemory(id);
  //fcuda::Malloc(&aug.Resultg,1);
  const unsigned npt12=MeshPts.npt1*MeshPts.npt2;
  const unsigned npt  =MeshPts.npt;
  if(ComputeRhop  )fcuda::Malloc(&aug.DataRhopg ,npt);
  if(ComputeVelxyz)fcuda::Malloc(&aug.DataVxyzg ,npt);
  if(ComputeVeldir)fcuda::Malloc(&aug.DataVdirg ,npt);
  if(ComputeZsurf )fcuda::Malloc(&aug.DataZsurfg,npt12);
  if(ComputeZsurf )fcuda::Malloc(&aug.DataMassg ,npt);
  aug.GpuMemory=true; 
}

//==============================================================================
/// Calculates data at indicated mesh points (on GPU).
//==============================================================================
void JGaugeMesh::CalculeGpu(const StDataGpu& datagpu){
  const int id=0;
  const int indomain=2;
  //-Run GPU calculation.
  if(indomain==2){
    //-Prepare input data.
    const StDivDataGpu& dvd=datagpu.dvd;
    const double2*  posxy  =datagpu.posxy_g->cptr();
    const double*   posz   =datagpu.posz_g->cptr();
    const typecode* code   =datagpu.code_g->cptr();
    //const unsigned* idp  =datagpu.idp_g->cptr();
    const float4*   velrho =datagpu.velrho_g->cptr();
    //-Prepare output arrays.
    tfloat3* datavxyz =MeshDat->GetVarFloat3("Vel");
    float*   datavdir =MeshDat->GetVarFloat("VelDir");
    float*   datarhop =MeshDat->GetVarFloat("Rhop");
    float*   datazsurf=MeshDat->GetVarFloat("Zsurf");
    //-Start measure.
    StGaugeMeshDataGpu& aug=AuxDataGpu[id];
    if(!aug.GpuMemory)AllocGpuMemory(id);
    cugauge::ComputeGaugeMesh(CSP,dvd,MeshPts,KcLimit,KcDummy,posxy,posz,code
      ,velrho,aug.DataVxyzg,aug.DataVdirg,aug.DataRhopg,aug.DataMassg);
    if(ComputeZsurf)cugauge::ComputeGaugeMeshZsurf(MassLimit,MeshPts
      ,aug.DataMassg,aug.DataZsurfg);
    //-Copy data from GPU memory.
    const unsigned npt12=MeshPts.npt1*MeshPts.npt2;
    const unsigned npt=MeshPts.npt;
    if(aug.DataVxyzg )cudaMemcpy(datavxyz ,aug.DataVxyzg ,sizeof(float3)*npt  ,cudaMemcpyDeviceToHost);
    if(aug.DataVdirg )cudaMemcpy(datavdir ,aug.DataVdirg ,sizeof(float )*npt  ,cudaMemcpyDeviceToHost);
    if(aug.DataRhopg )cudaMemcpy(datarhop ,aug.DataRhopg ,sizeof(float )*npt  ,cudaMemcpyDeviceToHost);
    if(aug.DataZsurfg)cudaMemcpy(datazsurf,aug.DataZsurfg,sizeof(float )*npt12,cudaMemcpyDeviceToHost);
    Check_CudaErroor("Failed in GaugeMesh calculation.");
    //-Stores calculated result.
    {
      const double& timestep=datagpu.timestep;
      SetTimeStep(timestep);
      MeshDat->SetTimeStep(timestep);
      Result.Set(timestep,MeshDat);
      if(Output(timestep))StoreResult();
      //Log->Printf("------> t:%f",TimeStep);
    }
  }
}

#endif
//<vs_meeshdat_end>


//##############################################################################
//# JGaugeForce
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeForce::JGaugeForce(unsigned idx,std::string name,word mkbound
  ,TpParticles typeparts,unsigned idbegin,unsigned count,typecode code
  ,tfloat3 center,int gpucount)
  :JGaugeItem(GAUGE_Force,idx,name,gpucount)
{
  ClassName="JGaugeForce";
  FileInfo=string("Saves Force data measured from boundary particles (by ")+ClassName+").";
  PartAcec=NULL; 
  Reset();
  MkBound=mkbound;
  TypeParts=typeparts;
  IdBegin=idbegin;
  Count=count;
  Code=code;
  InitialCenter=center;
  //-Allocates memory to calculate acceleration in selected particles on CPU.
  if(Cpu)PartAcec=new tfloat3[Count];
}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeForce::~JGaugeForce(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeForce::Reset(){
  MkBound=0;
  TypeParts=TpPartUnknown;
  IdBegin=Count=0;
  Code=0;
  InitialCenter=TFloat3(0);
  delete[] PartAcec; PartAcec=NULL;
  #ifdef _WITHGPU
    ResetGpuMemory();
  #endif
  JGaugeItem::Reset();
}

#ifdef _WITHGPU
//==============================================================================
/// Reset GPU memory.
//==============================================================================
void JGaugeForce::ResetGpuMemory(){
  for(int g=0;g<GpuCount;g++){
    if(AuxDataGpu[g].GpuMemory)
      Run_Exceptioon(fun::PrintStr("Auxiliary GPU memory for unit %d is not free.",g));
    AuxDataGpu[g]=StrGaugeForceDataGpu();
  }
}
#endif

//==============================================================================
/// Configuration of initial limits of calculation area.
//==============================================================================
void JGaugeForce::ConfigDomMCel(bool fixed){
  if(!fixed)Run_Exceptioon("Unfixed calculation domain is invalid for Force-Gauge.");
  FixedDomMCel=fixed;
  DomMCelIni0=TUint3(0);
  DomMCelFin0=TUint3(UINT_MAX);
}

//==============================================================================
/// Record the last measure result.
//==============================================================================
void JGaugeForce::StoreResult(){
  if(OutputSave){
    //-Allocates memory.
    while(unsigned(OutBuff.size())<OutSize)OutBuff.push_back(StrGaugeForceRes());
    //-Empty buffer.
    if(OutCount+1>=OutSize)SaveResults();
    //-Stores last results.
    OutBuff[OutCount]=Result;
    OutCount++;
    //-Updates OutputNext.
    if(OutputDt){
      const unsigned nt=unsigned(TimeStep/OutputDt);
      OutputNext=OutputDt*nt;
      if(OutputNext<=TimeStep)OutputNext=OutputDt*(nt+1);
    }
  }
}

//==============================================================================
/// Saves stored results in CSV file.
//==============================================================================
void JGaugeForce::SaveResults(){
  if(OutCount){
    const bool first=OutFile.empty();
    if(first){
      OutFile=GetResultsFileCsv();
      Log->AddFileInfo(OutFile,FileInfo);
    }
    jcsv::JSaveCsv2 scsv(OutFile,!first,AppInfo.GetCsvSepComa());
    //-Saves head.
    if(first){
      //-Head of values.
      scsv.SetHead();
      scsv <<"time [s];force [N];forcex [N];forcey [N];forcez [N]" << jcsv::Endl();
    }
    //-Saves data.
    scsv.SetData();
    scsv << jcsv::Fmt(jcsv::TpFloat1,"%g") << jcsv::Fmt(jcsv::TpFloat3,"%g;%g;%g");
    for(unsigned c=0;c<OutCount;c++){
      scsv << OutBuff[c].timestep << fgeo::PointDist(OutBuff[c].force) << OutBuff[c].force << jcsv::Endl();
    }
    OutCount=0;
  }
}

//==============================================================================
/// Saves last result in VTK file.
//==============================================================================
void JGaugeForce::SaveVtkResult(unsigned cpart){
  JDataArrays arrays;
  arrays.AddArray("Pos",1,&InitialCenter,false);
  arrays.AddArray("Force",1,&(Result.force),false);
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  JSpVtkData::Save(fun::FileNameSec(GetResultsFileVtk(),cpart),arrays,"Pos");
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeForce::GetPointDef(std::vector<tfloat3>& points)const{
  points.push_back(InitialCenter);
  return(1);
}

//==============================================================================
/// Calculates force sumation on selected fixed or moving particles using only fluid particles (on CPU).
/// Ignores periodic boundary particles to avoid race condition problems.
//==============================================================================
template<TpKernel tker> void JGaugeForce::CalculeCpuT(const StDataCpu& datacpu){
  //-Prepare input data.
  const double& timestep =datacpu.timestep;
  const StDivDataCpu& dvd=datacpu.dvd;
  const tdouble3* pos    =datacpu.pos_c->cptr();
  const typecode* code   =datacpu.code_c->cptr();
  const unsigned* idp    =datacpu.idp_c->cptr();
  const tfloat4*  velrho =datacpu.velrho_c->cptr();
  const unsigned  npbok  =datacpu.npbok;
  const unsigned  np     =datacpu.np;
  SetTimeStep(timestep);
  //-Computes acceleration in selected boundary particles.
  memset(PartAcec,0,sizeof(tfloat3)*Count);
  const int n=int(TypeParts==TpPartFixed || TypeParts==TpPartMoving? npbok: np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<n;p1++)if(CODE_GetTypeAndValue(code[p1])==Code && CODE_IsNormal(code[p1])){//-Ignores periodic boundaries.
    const tdouble3 pos1=pos[p1];
    const tfloat4 velrho1=velrho[p1];
    const float press1=fsph::ComputePress(velrho1.w,CSP);
    //-Auxiliary variables.
    tfloat3 ace=TFloat3(0);

    //-Search for fluid neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(pos1,false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const tfloat4 dr=nsearch::Distances(pos1,pos[p2]);
        //-Interaction with real neighboring particles. | Interaccion con particulas vecinas reales.
        if(dr.w<=CSP.kernelsize2 && dr.w>=ALMOSTZERO && CODE_IsFluid(code[p2])){
          const float fac=fsph::GetKernel_Fac<tker>(CSP,dr.w);
          const float frx=fac*dr.x;
          const float fry=fac*dr.y;
          const float frz=fac*dr.z;
          //-Velocity derivative (Momentum equation).
          const float mass2=CSP.massfluid;
          const tfloat4 velrho2=velrho[p2];
          const float press2=fsph::ComputePress(velrho2.w,CSP);
          //-Velocity derivative (Momentum equation).
          const float prs=(press1+press2)/(velrho1.w*velrho2.w)
            + (tker==KERNEL_Cubic? fsph::GetKernelCubic_Tensil(CSP,dr.w,velrho1.w,press1,velrho2.w,press2): 0);
          const float p_vpm1=-prs*mass2;
          ace.x+=p_vpm1*frx;  ace.y+=p_vpm1*fry;  ace.z+=p_vpm1*frz;
          //-Aceleration from viscosity not included.
        }
      }
    }
    //-Saves ace.
    PartAcec[idp[p1]-IdBegin]=ace;
  }
  //-Computes total ace.
  tfloat3 acesum=TFloat3(0);
  for(unsigned p=0;p<Count;p++)acesum=acesum+PartAcec[p];

  //-Stores result. | Guarda resultado.
  Result.Set(timestep,acesum*CSP.massbound);
  //Log->Printf("------> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}

//==============================================================================
/// Calculates force sumation on selected fixed or moving particles using only fluid particles (on CPU).
/// Ignores periodic boundary particles to avoid race condition problems.
//==============================================================================
void JGaugeForce::CalculeCpu(const StDataCpu& datacpu){
  switch(CSP.tkernel){
    case KERNEL_Cubic:       //Kernel Wendland is used since Cubic is not available.
    case KERNEL_Wendland:    CalculeCpuT<KERNEL_Wendland>  (datacpu);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Check if GPU memory was allocated.
//==============================================================================
bool JGaugeForce::AllocatedGpuMemory(int id)const{
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  return(AuxDataGpu[id].GpuMemory);
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JGaugeForce::FreeGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeForceDataGpu& aug=AuxDataGpu[id];
  aug.GpuMemory=false;
  if(aug.Resultg )cudaFree(aug.Resultg ); aug.Resultg=NULL;
  if(aug.PartAceg)cudaFree(aug.PartAceg); aug.PartAceg=NULL;
  if(aug.AuxSumg )cudaFree(aug.AuxSumg ); aug.AuxSumg=NULL;
}

//==============================================================================
/// Allocates GPU memory.
//==============================================================================
void JGaugeForce::AllocGpuMemory(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  StGaugeForceDataGpu& aug=AuxDataGpu[id];
  if(aug.GpuMemory)FreeGpuMemory(id);
  fcuda::Malloc(&aug.Resultg,1);
  fcuda::Malloc(&aug.PartAceg,Count);
  const unsigned saux=curedus::GetAuxSize_ReduSumFloat3(Count);
  fcuda::Malloc(&aug.AuxSumg,saux);
  aug.GpuMemory=true; 
}

//==============================================================================
/// Calculates force sumation on selected fixed or moving particles using only fluid particles (on GPU).
/// Ignores periodic boundary particles.
//==============================================================================
void JGaugeForce::CalculeGpu(const StDataGpu& datagpu){
  const int id=0;
  const int indomain=2;
  //-Run GPU calculation.
  if(indomain==1 || indomain==2){
    //-Prepare input data.
    const StDivDataGpu& dvd=datagpu.dvd;
    const double2*  posxy  =datagpu.posxy_g->cptr();
    const double*   posz   =datagpu.posz_g->cptr();
    const typecode* code   =datagpu.code_g->cptr();
    const unsigned* idp    =datagpu.idp_g->cptr();
    const float4*   velrho =datagpu.velrho_g->cptr();
    const unsigned  npbok  =datagpu.npbok;
    const unsigned  np     =datagpu.np;
    //-Start measure.
    StGaugeForceDataGpu& aug=AuxDataGpu[id];
    if(!aug.GpuMemory)AllocGpuMemory(id);

    //-Initializes acceleration array to zero.
    cudaMemset(aug.PartAceg,0,sizeof(float3)*Count);
    const int n=int(TypeParts==TpPartFixed || TypeParts==TpPartMoving? npbok: np);
    //-Computes acceleration in selected boundary particles.
    cugauge::Interaction_GaugeForce(CSP,dvd,n,IdBegin,Code
      ,posxy,posz,code,idp,velrho,aug.PartAceg);

    //-Computes total ace.
    tfloat3 acesum=TFloat3(0);
    {//-Computes total ace on GPU.
      const float3 result=curedus::ReduSumFloat3(Count,0,aug.PartAceg,aug.AuxSumg);
      acesum=TFloat3(result.x,result.y,result.z);
      Check_CudaErroor("Failed in Force calculation.");
    }
    //-Stores calculated result.
    if(indomain==2){
      const double& timestep=datagpu.timestep;
      SetTimeStep(timestep);
      Result.Set(timestep,acesum*CSP.massbound);
      if(Output(timestep))StoreResult();
      //Log->Printf("------> t:%f",TimeStep);
    }
  }
}

#endif

