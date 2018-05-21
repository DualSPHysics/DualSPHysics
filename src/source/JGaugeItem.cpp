//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JGauge.cpp \brief Implements the class \ref JGauge.

#include "JGaugeItem.h"
#include "JLog2.h"
#include "JSaveCsv2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JFormatFiles2.h"
#ifdef _WITHGPU
  #include "JGauge_ker.h"
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
JGaugeItem::JGaugeItem(TpGauge type,unsigned idx,std::string name,JLog2* log)
  :Type(type),Idx(idx),Name(name),Log(log)
{
  ClassName="JGaugeItem";
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeItem::Reset(){
  Config(false,TDouble3(0),TDouble3(0),0,0,0,0);
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
void JGaugeItem::Config(bool simulate2d,tdouble3 domposmin,tdouble3 domposmax
  ,float scell,int hdiv,float h,float massfluid)
{
  Simulate2D=simulate2d;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
  Scell=scell;
  Hdiv=hdiv;
  H=h;
  Fourh2=float(h*h*4); 
  Awen=(h? float(Simulate2D? 0.557/(h*h): 0.41778/(h*h*h)): 0);
  MassFluid=massfluid;
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
    case GAUGE_Vel:  return("Vel");
    case GAUGE_Swl:  return("SWL");
    case GAUGE_MaxZ: return("MaxZ");
  }
  return("???");
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JGaugeItem::GetConfig(std::vector<std::string> &lines)const{
  const char met[]="GetConfig";
  lines.push_back(fun::PrintStr("Type.......: %s",JGaugeItem::GetNameType(Type).c_str()));
  const string cpend=fun::DoublexStr(ComputeEnd,"%g");
  lines.push_back(fun::PrintStr("Compute....: %g - %s   dt:%g",ComputeStart,cpend.c_str(),ComputeDt));
  const string ouend=fun::DoublexStr(OutputEnd,"%g");
  lines.push_back(fun::PrintStr("Output.....:%s%g - %s   dt:%g",(OutputSave? " ": " (disabled)  "),OutputStart ,ouend.c_str(),OutputDt));
  lines.push_back(fun::PrintStr("SaveVtkPart: %s",(SaveVtkPart? "True": "False")));
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
  else RunException(met,"Type unknown.");
}

//==============================================================================
/// Returns filename for output results for CSV files.
//==============================================================================
std::string JGaugeItem::GetResultsFileCsv()const{
  return(AppInfo.GetDirOut()+"Gauges"+GetNameType(Type)+"_"+Name)+".csv";
}

//==============================================================================
/// Returns filename for output results for VTK files.
//==============================================================================
std::string JGaugeItem::GetResultsFileVtk()const{
  return(AppInfo.GetDirDataOut()+"Gauges"+GetNameType(Type)+"_"+Name)+".vtk";
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
/// Return cell limits for interaction starting from position.
/// Devuelve limites de celdas para interaccion a partir de posicion.
//==============================================================================
void JGaugeItem::GetInteractionCells(const tdouble3 &pos,const tint4 &nc,const tint3 &cellzero
  ,int &cxini,int &cxfin,int &yini,int &yfin,int &zini,int &zfin)const
{
  //-Get cell coordinates of position pos.
  const int cx=int((pos.x-DomPosMin.x)/Scell)-cellzero.x;
  const int cy=int((pos.y-DomPosMin.y)/Scell)-cellzero.y;
  const int cz=int((pos.z-DomPosMin.z)/Scell)-cellzero.z;
  //-code for hdiv 1 or 2 but not zero. | Codigo para hdiv 1 o 2 pero no cero.
  cxini=cx-min(cx,Hdiv);
  cxfin=cx+min(nc.x-cx-1,Hdiv)+1;
  yini=cy-min(cy,Hdiv);
  yfin=cy+min(nc.y-cy-1,Hdiv)+1;
  zini=cz-min(cz,Hdiv);
  zfin=cz+min(nc.z-cz-1,Hdiv)+1;
}

#ifdef _WITHGPU
//==============================================================================
/// Throws exception for Cuda error.
/// Lanza excepcion por un error Cuda.
//==============================================================================
void JGaugeItem::RunExceptionCuda(const std::string &method,const std::string &msg,cudaError_t error){
  std::string tx=fun::PrintStr("%s (CUDA error: %s).\n",msg.c_str(),cudaGetErrorString(error)); 
  Log->Print(GetExceptionText(method,tx));
  RunException(method,msg);
}

//==============================================================================
/// Check error and throw exception if there was one. 
/// Comprueba error y lanza excepcion si lo hubiera.
//==============================================================================
void JGaugeItem::CheckCudaError(const std::string &method,const std::string &msg){
  cudaError_t err=cudaGetLastError();
  if(err!=cudaSuccess)RunExceptionCuda(method,msg,err);
}
#endif

//##############################################################################
//# JGaugeVelocity
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeVelocity::JGaugeVelocity(unsigned idx,std::string name,tdouble3 point,JLog2* log)
  :JGaugeItem(GAUGE_Vel,idx,name,log)
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
  JGaugeItem::Reset();
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
  std::vector<JFormatFiles2::StScalarData> fields;
  fields.push_back(JFormatFiles2::DefineField("Vel",JFormatFiles2::Float32,3,&(Result.vel)));
  //-Saves VTK file.
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  JFormatFiles2::SaveVtk(fun::FileNameSec(GetResultsFileVtk(),cpart),1,&(Result.point),fields);
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeVelocity::GetPointDef(std::vector<tfloat3> &points)const{
  points.push_back(ToTFloat3(Point));
  return(1);
}

//==============================================================================
/// Calculates velocity at indicated points (on CPU).
//==============================================================================
void JGaugeVelocity::CalculeCpu(double timestep,tuint3 ncells,tuint3 cellmin,const unsigned *begincell
  ,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)
{
  SetTimeStep(timestep);
  //-Start measure.
  tfloat3 ptvel=TFloat3(0);
  const bool ptout=PointIsOut(Point.x,Point.y,Point.z);//-Verify that the point is within domain boundaries. | Comprueba que el punto este dentro de limites del dominio.
  if(!ptout){
    const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
    const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
    const unsigned cellfluid=nc.w*nc.z+1;
    //-Obtain limits of interaction. | Obtiene limites de interaccion.
    int cxini,cxfin,yini,yfin,zini,zfin;
    GetInteractionCells(Point,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin);

    //-Auxiliary variables.
    double sumwab=0;
    tdouble3 sumvel=TDouble3(0);

    //-Search for neighbors in adjacent cells.
    //-Busqueda de vecinos en celdas adyacentes.
    if(cxini<cxfin)for(int z=zini;z<zfin;z++){
      const int zmod=(nc.w)*z+cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        const unsigned pini=begincell[cxini+ymod];
        const unsigned pfin=begincell[cxfin+ymod];

        //-Interaction with Fluid/Floating | Interaccion con varias Fluid/Floating.
        //--------------------------------------------------------------------------
        for(unsigned p2=pini;p2<pfin;p2++){
          const float drx=float(Point.x-pos[p2].x);
          const float dry=float(Point.y-pos[p2].y);
          const float drz=float(Point.z-pos[p2].z);
          const float rr2=(drx*drx+dry*dry+drz*drz);
          //-Interaction with real neighboring particles.
          //-Interaccion con particulas vecinas reales.
          if(rr2<=Fourh2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
            float wab;
            {//-Wendland kernel.
              const float qq=sqrt(rr2)/H;
              const float wqq=2.f*qq+1.f;
              const float wqq1=1.f-0.5f*qq;
              const float wqq2=wqq1*wqq1;
              wab=Awen*wqq*wqq2*wqq2;
            }
            wab*=MassFluid/velrhop[p2].w;
            sumwab+=wab;
            sumvel.x+=wab*velrhop[p2].x;
            sumvel.y+=wab*velrhop[p2].y;
            sumvel.z+=wab*velrhop[p2].z;
          }
        }
      }
    }
    //-Aplies kernel correction.
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

#ifdef _WITHGPU
//==============================================================================
/// Calculates velocity at indicated points (on GPU).
//==============================================================================
void JGaugeVelocity::CalculeGpu(double timestep,tuint3 ncells,tuint3 cellmin,const int2 *beginendcell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,float3 *aux)
{
  SetTimeStep(timestep);
  //-Start measure.
  tfloat3 ptvel=TFloat3(0);
  const bool ptout=PointIsOut(Point.x,Point.y,Point.z);//-Verify that the point is within domain boundaries. | Comprueba que el punto este dentro de limites del dominio.
  if(!ptout){
    cugauge::Interaction_GaugeVel(Point,Awen,Hdiv,ncells,cellmin,beginendcell,posxy,posz,code,velrhop,aux,DomPosMin,Scell,Fourh2,H,MassFluid);
    cudaMemcpy(&ptvel,aux,sizeof(float3),cudaMemcpyDeviceToHost);
    CheckCudaError("CalculeGpu","Failed in velocity calculation.");
  }
  //-Stores result. | Guarda resultado.
  Result.Set(timestep,ToTFloat3(Point),ptvel);
  //Log->Printf("------> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}
#endif


//##############################################################################
//# JGaugeSwl
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeSwl::JGaugeSwl(unsigned idx,std::string name,tdouble3 point0,tdouble3 point2,double pointdp,float masslimit,JLog2* log)
  :JGaugeItem(GAUGE_Swl,idx,name,log)
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
  SetPoints(TDouble3(0),TDouble3(0),0);
  MassLimit=0;
  JGaugeItem::Reset();
}

//==============================================================================
/// Changes points definition.
//==============================================================================
void JGaugeSwl::SetPoints(const tdouble3 &point0,const tdouble3 &point2,double pointdp){
  //if(PointDp<=0)RunException("SetPoints",fun::PrintStr("The value of PointDp is <= zero in gauge \'%s\'.",Name.c_str()));
  Point0=point0;
  Point2=point2;
  PointDp=pointdp;
  const double dis=fmath::DistPoints(Point0,Point2);
  if(dis>0 && PointDp>0){
    PointNp=unsigned(dis/PointDp);
    if(dis-(PointDp*PointNp)>=PointDp*0.1)PointNp++;
    if(PointNp<1)PointNp++;
    const double dp=dis/PointNp;
    PointDir=fmath::VecUnitary(Point2-Point0)*dp;
  }
  else{
    PointNp=0;
    PointDir=TDouble3(0);
  } 
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
  std::vector<JFormatFiles2::StScalarData> fields;
  //-Saves VTK file.
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  JFormatFiles2::SaveVtk(fun::FileNameSec(GetResultsFileVtk(),cpart),1,&(Result.posswl),fields);
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeSwl::GetPointDef(std::vector<tfloat3> &points)const{
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
float JGaugeSwl::CalculeMassCpu(const tdouble3 &ptpos,const tint4 &nc
  ,const tint3 &cellzero,unsigned cellfluid,const unsigned *begincell
  ,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)const
{
  //-Obtain limits of interaction. | Obtiene limites de interaccion.
  int cxini,cxfin,yini,yfin,zini,zfin;
  GetInteractionCells(ptpos,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin);

  //-Auxiliary variables.
  double sumwab=0;
  double summass=0;

  //-Search for neighbors in adjacent cells.
  //-Busqueda de vecinos en celdas adyacentes.
  if(cxini<cxfin)for(int z=zini;z<zfin;z++){
    const int zmod=(nc.w)*z+cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
    for(int y=yini;y<yfin;y++){
      int ymod=zmod+nc.x*y;
      const unsigned pini=begincell[cxini+ymod];
      const unsigned pfin=begincell[cxfin+ymod];

      //-Interaction with Fluid/Floating | Interaccion con varias Fluid/Floating.
      //--------------------------------------------------------------------------
      for(unsigned p2=pini;p2<pfin;p2++){
        const float drx=float(ptpos.x-pos[p2].x);
        const float dry=float(ptpos.y-pos[p2].y);
        const float drz=float(ptpos.z-pos[p2].z);
        const float rr2=(drx*drx+dry*dry+drz*drz);
        //-Interaction with real neighboring particles.
        //-Interaccion con particulas vecinas reales.
        if(rr2<=Fourh2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2])){
          float wab;
          {//-Wendland kernel.
            const float qq=sqrt(rr2)/H;
            const float wqq=2.f*qq+1.f;
            const float wqq1=1.f-0.5f*qq;
            const float wqq2=wqq1*wqq1;
            wab=Awen*wqq*wqq2*wqq2;
          }
          //:Log->Printf("----> p2:%u  wab:%f  vol:%f",p2,wab,MassFluid/velrhop[p2].w);
          wab*=MassFluid/velrhop[p2].w;
          sumwab+=wab;
          summass+=wab*MassFluid;
        }
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
void JGaugeSwl::CalculeCpu(double timestep,tuint3 ncells,tuint3 cellmin,const unsigned *begincell
  ,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)
{
  SetTimeStep(timestep);
  const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
  const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
  const unsigned cellfluid=nc.w*nc.z+1;

  //-Look for change of fluid to empty. | Busca paso de fluido a vacio.
  tdouble3 ptsurf=TDouble3(DBL_MAX);
  float mpre=0;
  tdouble3 ptpos=Point0;
  for(unsigned cp=0;cp<=PointNp;cp++){
    const float mass=CalculeMassCpu(ptpos,nc,cellzero,cellfluid,begincell,pos,code,velrhop);
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
  //Log->Printf("------> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}

#ifdef _WITHGPU
//==============================================================================
/// Calculates surface water level at indicated points (on GPU).
//==============================================================================
void JGaugeSwl::CalculeGpu(double timestep,tuint3 ncells,tuint3 cellmin,const int2 *beginendcell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,float3 *aux)
{
  SetTimeStep(timestep);
  //-Start measure.
  cugauge::Interaction_GaugeSwl(Point0,PointDir,PointNp,MassLimit,Awen,Hdiv,ncells,cellmin,beginendcell,posxy,posz,code,velrhop,DomPosMin,Scell,Fourh2,H,MassFluid,aux);
  tfloat3 ptsurf=TFloat3(0);
  cudaMemcpy(&ptsurf,aux,sizeof(float3),cudaMemcpyDeviceToHost);
  CheckCudaError("CalculeGpu","Failed in Swl calculation.");
  //-Stores result. | Guarda resultado.
  Result.Set(timestep,ToTFloat3(Point0),ToTFloat3(Point2),ptsurf);
  //Log->Printf("------> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}
#endif


//##############################################################################
//# JGaugeMaxZ
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeMaxZ::JGaugeMaxZ(unsigned idx,std::string name,tdouble3 point0,double height,float distlimit,JLog2* log)
  :JGaugeItem(GAUGE_MaxZ,idx,name,log)
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
  JGaugeItem::Reset();
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
  std::vector<JFormatFiles2::StScalarData> fields;
  fields.push_back(JFormatFiles2::DefineField("Height",JFormatFiles2::Float32,1,&height));
  //-Saves VTK file.
  Log->AddFileInfo(fun::FileNameSec(GetResultsFileVtk(),UINT_MAX),FileInfo);
  JFormatFiles2::SaveVtk(fun::FileNameSec(GetResultsFileVtk(),cpart),1,&ptz,fields);
}

//==============================================================================
/// Loads and returns number definition points.
//==============================================================================
unsigned JGaugeMaxZ::GetPointDef(std::vector<tfloat3> &points)const{
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
void JGaugeMaxZ::GetInteractionCellsMaxZ(const tdouble3 &pos,const tint4 &nc,const tint3 &cellzero
  ,int &cxini,int &cxfin,int &yini,int &yfin,int &zini,int &zfin)const
{
  //-Get cell coordinates of position pos.
  const int cx=int((pos.x-DomPosMin.x)/Scell)-cellzero.x;
  const int cy=int((pos.y-DomPosMin.y)/Scell)-cellzero.y;
  const int cz=int((pos.z-DomPosMin.z)/Scell)-cellzero.z;
  const int cz2=int((pos.z+Height-DomPosMin.z)/Scell)-cellzero.z;
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
void JGaugeMaxZ::CalculeCpu(double timestep,tuint3 ncells,tuint3 cellmin,const unsigned *begincell
  ,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)
{
  //Log->Printf("JGaugeMaxZ----> timestep:%g  (%d)",timestep,(DG?1:0));
  SetTimeStep(timestep);
  //-Compute auxiliary constants.
  const float maxdist2=DistLimit*DistLimit;
  //-Obtain limits of interaction.
  const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
  const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
  const unsigned cellfluid=nc.w*nc.z+1;
  int cxini,cxfin,yini,yfin,zini,zfin;
  GetInteractionCellsMaxZ(Point0,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin);
  //-Start measure.
  unsigned pmax=UINT_MAX;
  float zmax=-FLT_MAX;
  //-Search for neighbours in adjacent cells. | Busqueda de vecinos en celdas adyacentes.
  if(cxini<cxfin)for(int z=zfin-1;z>=zini && pmax==UINT_MAX;z--){
    const int zmod=(nc.w)*z+cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
    for(int y=yini;y<yfin;y++){
      int ymod=zmod+nc.x*y;
      const unsigned pini=begincell[cxini+ymod];
      const unsigned pfin=begincell[cxfin+ymod];

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
/// Calculates maximum z of fluid at distance of a vertical line (on GPU).
//==============================================================================
void JGaugeMaxZ::CalculeGpu(double timestep,tuint3 ncells,tuint3 cellmin,const int2 *beginendcell
  ,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,float3 *aux)
{
  SetTimeStep(timestep);
  //-Compute auxiliary constants.
  const float maxdist2=DistLimit*DistLimit;
  //-Obtain limits of interaction.
  const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
  const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
  const unsigned cellfluid=nc.w*nc.z+1;
  int cxini,cxfin,yini,yfin,zini,zfin;
  GetInteractionCellsMaxZ(Point0,nc,cellzero,cxini,cxfin,yini,yfin,zini,zfin);
  //-Start measure.
  cugauge::Interaction_GaugeMaxz(Point0,maxdist2,cxini,cxfin,yini,yfin,zini,zfin,cugauge::Int4(nc),cellfluid,beginendcell,posxy,posz,code,aux);
  tfloat3 ptsurf=TFloat3(0);
  cudaMemcpy(&ptsurf,aux,sizeof(float3),cudaMemcpyDeviceToHost);
  CheckCudaError("CalculeGpu","Failed in MaxZ calculation.");
  //-Stores result. | Guarda resultado.
  Result.Set(timestep,ToTFloat3(Point0),ptsurf.z);
  //Log->Printf("------> t:%f",TimeStep);
  if(Output(timestep))StoreResult();
}
#endif

