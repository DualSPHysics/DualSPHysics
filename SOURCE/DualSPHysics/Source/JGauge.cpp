//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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

#include "JGauge.h"
#include "JLog2.h"
#include "JSaveCsv.h"
#include "Functions.h"
#ifdef _WITHGPU
  #include "JGauge_ker.h"
#endif
#include <cmath>
#include <cfloat>

using std::string;
using std::min;
using std::max;

//##############################################################################
//# JGaugeBase
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeBase::JGaugeBase(JLog2* log):Log(log){
  ClassName="JGaugeBase";
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeBase::Reset(){
  ConfigConstants(SetCtes(false,ORDER_XYZ,0,0,0,0,0,TDouble3(0),TDouble3(0),TDouble3(0)));
}

//==============================================================================
/// Configures constants during simulation.
//==============================================================================
void JGaugeBase::ConfigConstants(StCtes ctes){
  Simulate2D=ctes.simulate2d;
  CellOrder=ctes.cellorder;
  MassFluid=ctes.massfluid;
  Dp=ctes.dp;
  Dosh=ctes.dosh;
  Scell=ctes.scell;
  Hdiv=ctes.hdiv;
  DomPosMin=ctes.domposmin;
  DomRealPosMin=ctes.domrealposmin;
  DomRealPosMax=ctes.domrealposmax;
  if(CellOrder!=ORDER_XYZ)RunException("ConfigConstants","Only order XYZ is valid for now...");
  //-Calculation of constants. | Calculo de constantes.
  if(Dosh){
    double h=double(Dosh)/2;
    H=float(h);
    Awen=float(Simulate2D? 0.557/(h*h): 0.41778/(h*h*h));
    Fourh2=float(h*h*4); 
  }
  else H=Awen=Fourh2=0;
}

#ifdef _WITHGPU
//==============================================================================
/// Throws exception for Cuda error.
/// Lanza excepcion por un error Cuda.
//==============================================================================
void JGaugeBase::RunExceptionCuda(const std::string &method,const std::string &msg,cudaError_t error){
  std::string tx=fun::PrintStr("%s (CUDA error: %s).\n",msg.c_str(),cudaGetErrorString(error)); 
  Log->Print(GetExceptionText(method,tx));
  RunException(method,msg);
}

//==============================================================================
/// Check error and throw exception if there was one. 
/// Comprueba error y lanza excepcion si lo hubiera.
//==============================================================================
void JGaugeBase::CheckCudaError(const std::string &method,const std::string &msg){
  cudaError_t err=cudaGetLastError();
  if(err!=cudaSuccess)RunExceptionCuda(method,msg,err);
}
#endif

//##############################################################################
//# JGaugeVel
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeVel::JGaugeVel(JLog2* log):JGaugeBase(log){
  ClassName="JGaugeVel";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeVel::~JGaugeVel(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeVel::Reset(){
  JGaugeBase::Reset();
  TimeStep=0;
  PtOut=true;
  PtPos=TDouble3(0);
  PtVel=TFloat3(0);
}

//==============================================================================
/// Configures object.
//==============================================================================
void JGaugeVel::Init(StCtes ctes){
  Reset();
  ConfigConstants(ctes);
}

//==============================================================================
/// Record the last measurement result.
/// Graba ultimo resultado de medicion.
//==============================================================================
void JGaugeVel::SaveResult(){
  JSaveCsv scsv(Log->GetDirOut()+"_JGaugeVel.csv",true,Log->GetCsvSepComa());
  scsv.AddHead("Time;Velx;Velz");
  scsv.AddValuesf("%g;%g;%g",TimeStep,PtVel.x,PtVel.z);
  scsv.AddEndl();
}

//==============================================================================
/// Calculates velocity in indicated point.
//==============================================================================
tfloat3 JGaugeVel::CalculeCpu(double timestep,tdouble3 ptpos,tuint3 ncells,tuint3 cellmin
  ,const unsigned *begincell,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)
{
  TimeStep=timestep;
  PtPos=ptpos;
  PtVel=TFloat3(0);
  PtOut=PointIsOut(ptpos.x,ptpos.y,ptpos.z);//-Verify that the point is within real domain boundaries. | Comprueba que el punto este dentro de limites del dominio real.
  if(!PtOut){
    const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
    const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
    const unsigned cellfluid=nc.w*nc.z+1;
    const double px=ptpos.x;
    const double py=ptpos.y;
    const double pz=ptpos.z;

    int cx=int((px-DomPosMin.x)/Scell)-cellzero.x;
    int cy=int((py-DomPosMin.y)/Scell)-cellzero.y;
    int cz=int((pz-DomPosMin.z)/Scell)-cellzero.z;
    
    //-Code for hdiv 1 or 2. The result is always within the range [0,nc-1+1].
    //-Codigo para hdiv 1 o 2. El resultado siempre esta dentro del intervalo [0,nc-1+1].
    const int cxini=cx-min(cx,Hdiv);
    const int cxfin=cx+min(nc.x-cx-1,Hdiv)+1;
    const int yini=cy-min(cy,Hdiv);
    const int yfin=cy+min(nc.y-cy-1,Hdiv)+1;
    const int zini=cz-min(cz,Hdiv);
    const int zfin=cz+min(nc.z-cz-1,Hdiv)+1;

    //-Auxiliary variables.
    double sumwab=0;
    tdouble3 sumvel=TDouble3(0);

    //-Search for neighbors in adjacent cells.
    //-Busqueda de vecinos en celdas adyacentes.
    for(int z=zini;z<zfin;z++){
      const int zmod=(nc.w)*z+cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
      for(int y=yini;y<yfin;y++){
        int ymod=zmod+nc.x*y;
        const unsigned pini=begincell[cxini+ymod];
        const unsigned pfin=begincell[cxfin+ymod];

        //-Interaction with Fluid/Floating | Interaccion con varias Fluid/Floating.
        //--------------------------------------------------------------------------
        for(unsigned p2=pini;p2<pfin;p2++){
          const float drx=float(px-pos[p2].x);
          const float dry=float(py-pos[p2].y);
          const float drz=float(pz-pos[p2].z);
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
    PtVel=ToTFloat3(sumvel);
  }
  //SaveResult();
  return(PtVel);
}

#ifdef _WITHGPU
//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
tfloat3 JGaugeVel::CalculeGpu(double timestep,tdouble3 ptpos,tuint3 ncells,tuint3 cellmin
  ,const int2 *beginendcell,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,float3 *aux)
{
  TimeStep=timestep;
  PtPos=ptpos;
  PtVel=TFloat3(0);
  PtOut=PointIsOut(ptpos.x,ptpos.y,ptpos.z);//-Verify that the point is within real domain boundaries. | Comprueba que el punto este dentro de limites del dominio real.
  if(!PtOut){
    cugauge::Interaction_GaugeVel(ptpos,Awen,Hdiv,ncells,cellmin,beginendcell,posxy,posz,code,velrhop,aux,DomPosMin,Scell,Fourh2,H,MassFluid);
    //-Stores result. | Guarda resultado.
    cudaMemcpy(&PtVel,aux,sizeof(float3),cudaMemcpyDeviceToHost);
    CheckCudaError("CalculeGpu","Failed in velocity interpolation.");
  }
  //SaveResult();
  return(PtVel);
}
#endif


//##############################################################################
//# JGaugeZsurf
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeZsurf::JGaugeZsurf(JLog2* log):JGaugeBase(log){
  ClassName="JGaugeZsurf";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeZsurf::~JGaugeZsurf(){
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeZsurf::Reset(){
  JGaugeBase::Reset();
  GaugeZmin=GaugeZmax=GaugeDp=0;
  MassLimit=0;
  TimeStep=0;
  PtOut=true;
  PtPos=TDouble3(0);
  PtzMin=Zsurf=0;
}

//==============================================================================
/// Configures object.
//==============================================================================
void JGaugeZsurf::Init(StCtes ctes,double gaugezmin,double gaugezmax,double gaugedp,float masslimit){
  Reset();
  ConfigConstants(ctes);
  GaugeZmin=gaugezmin;
  GaugeZmax=gaugezmax;
  GaugeDp=gaugedp;
  MassLimit=masslimit;
}

//==============================================================================
/// Record the last measurement result.
/// Graba ultimo resultado de medicion.
//==============================================================================
void JGaugeZsurf::SaveResult(){
  JSaveCsv scsv(Log->GetDirOut()+"_JGaugeZsurf.csv",true,Log->GetCsvSepComa());
  scsv.AddHead("Time;Zsurf");
  scsv.AddValuesf("%g;%g",TimeStep,Zsurf);
  scsv.AddEndl();
}

//==============================================================================
/// Calculates measurement limits according cell domain and last measurement.
/// Calcula limites de medicion segun dominio en celdas y ultima medicion.
//==============================================================================
void JGaugeZsurf::GetLimitsFreeSurface(double px,double py,double zsurf,tuint3 ncells,tuint3 cellmin,unsigned &czmin,unsigned &czmax)const{
  double zmin=DomPosMin.z+Scell*cellmin.z+Dp;
  double zmax=zmin+Scell*ncells.z-Dp;
  if(zmin<zsurf-Dosh)zmin=zsurf-Dosh;//-Increases zmin according to last measurement minus Dosh. | Aumenta zmin segun ultima medicion menos Dosh.
  zmin=max(zmin,GaugeZmin);
  zmax=min(zmax,GaugeZmax);
  czmin=unsigned((zmin-GaugeZmin)/GaugeDp);
  czmax=unsigned((zmax-GaugeZmin)/GaugeDp);
}

//==============================================================================
/// Returns the interpolated mass value at the indicated point. The point must 
/// belong to the cell domain.
///
/// Devuelve valor de masa interpolado en el punto indicado. El punto debe 
/// pertenecer al dominio de celdas.
//==============================================================================
float JGaugeZsurf::CalculeMass(tdouble3 ptpos,tuint3 ncells,tuint3 cellmin,const unsigned *begincell
  ,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)const
{
  const tint4 nc=TInt4(int(ncells.x),int(ncells.y),int(ncells.z),int(ncells.x*ncells.y));
  const tint3 cellzero=TInt3(cellmin.x,cellmin.y,cellmin.z);
  const unsigned cellfluid=nc.w*nc.z+1;
  const double px=ptpos.x;
  const double py=ptpos.y;
  const double pz=ptpos.z;

  int cx=int((px-DomPosMin.x)/Scell)-cellzero.x;
  int cy=int((py-DomPosMin.y)/Scell)-cellzero.y;
  int cz=int((pz-DomPosMin.z)/Scell)-cellzero.z;

  //-Code for hdiv 1 or 2. The result is always within the range [0,nc-1+1].
  //-Codigo para hdiv 1 o 2. El resultado siempre esta dentro del intervalo [0,nc-1+1].
  const int cxini=cx-min(cx,Hdiv);
  const int cxfin=cx+min(nc.x-cx-1,Hdiv)+1;
  const int yini=cy-min(cy,Hdiv);
  const int yfin=cy+min(nc.y-cy-1,Hdiv)+1;
  const int zini=cz-min(cz,Hdiv);
  const int zfin=cz+min(nc.z-cz-1,Hdiv)+1;

  //-Auxiliary variables.
  double sumwab=0;
  double summass=0;

  //-Search for neighbors in adjacent cells.
  //-Busqueda de vecinos en celdas adyacentes.
  for(int z=zini;z<zfin;z++){
    const int zmod=(nc.w)*z+cellfluid; //-Sum from start of fluid cells. | Le suma donde empiezan las celdas de fluido.
    for(int y=yini;y<yfin;y++){
      int ymod=zmod+nc.x*y;
      const unsigned pini=begincell[cxini+ymod];
      const unsigned pfin=begincell[cxfin+ymod];

      //-Interaction with Fluid/Floating | Interaccion con varias Fluid/Floating.
      //--------------------------------------------------------------------------
      for(unsigned p2=pini;p2<pfin;p2++){
        const float drx=float(px-pos[p2].x);
        const float dry=float(py-pos[p2].y);
        const float drz=float(pz-pos[p2].z);
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
  //-Aplies kernel correction.
  //if(sumwab!=0)summass/=sumwab;
  return(float(summass));
}

//==============================================================================
/// Calculates free-surface water in indicated point (x,y) on CPU. 
/// The zsurf0 must be the last zsurf value or z minimum.
//==============================================================================
double JGaugeZsurf::CalculeCpu(double timestep,double px,double py,double zsurf0,tuint3 ncells,tuint3 cellmin
  ,const unsigned *begincell,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop)
{
  TimeStep=timestep;
  PtPos=TDouble3(px,py,zsurf0);
  PtzMin=GaugeZmin;
  Zsurf=0; 
  PtOut=PointIsOut(px,py);//-Verify that the point is within real domain boundaries. | Comprueba que el punto este dentro de limites del dominio real.
  if(!PtOut){
    //-Prepares zsurf calculation. | Prepara calculo de zsurf.
    tdouble3 ps=TDouble3(px,py,0);
    unsigned czmin,czmax;
    GetLimitsFreeSurface(px,py,zsurf0,ncells,cellmin,czmin,czmax);
    PtzMin=GaugeZmin+GaugeDp*czmin;
    //-Look for change of fluid to empty. | Busca paso de fluido a vacio.
    double zsurf=DBL_MAX;
    double zpre=GaugeZmin;
    float mpre=MassLimit;
    for(unsigned cz=czmin;cz<=czmax && zsurf==DBL_MAX;cz++){
      ps.z=GaugeZmin+GaugeDp*cz;
      float mass=CalculeMass(ps,ncells,cellmin,begincell,pos,code,velrhop);
      if(mass<MassLimit){
        const double fx=(MassLimit-mpre)/(mass-mpre);
        zsurf=fx*(ps.z-zpre)+zpre;
      }
      zpre=ps.z;
      mpre=mass;
    }
    if(zsurf==DBL_MAX)zsurf=GaugeZmax;
    //-Stores result. | Guarda resultado.
    Zsurf=zsurf;
  }
  //SaveResult();
  return(Zsurf);
}

#ifdef _WITHGPU
//==============================================================================
/// Calculates free-surface water in indicated point (x,y) on CPU. 
/// The zsurf0 must be the last zsurf value or z minimum.
//==============================================================================
double JGaugeZsurf::CalculeGpu(double timestep,double px,double py,double zsurf0,tuint3 ncells,tuint3 cellmin
  ,const int2 *beginendcell,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,double *aux)
{
  TimeStep=timestep;
  PtPos=TDouble3(px,py,zsurf0);
  PtzMin=GaugeZmin;
  Zsurf=0;
  PtOut=PointIsOut(px,py);//-Verify that the point is within real domain boundaries. | Comprueba que el punto este dentro de limites del dominio real.
  if(!PtOut){
    //-Prepares zsurf calculation. | Prepara calculo de zsurf.
    unsigned czmin,czmax;
    GetLimitsFreeSurface(px,py,zsurf0,ncells,cellmin,czmin,czmax);
    //:Log->Printf("WD-> zsurf0:%f",zsurf0);
    //:Log->Printf("WD-> czmin-max:%u - %u (%u)  pzmin-max:%f - %f",czmin,czmax,(czmax-czmin+1),GaugeZmin+GaugeDp*czmin,GaugeZmin+GaugeDp*czmax);
    PtzMin=GaugeZmin+GaugeDp*czmin;
    //-Look for change of fluid to empty. | Busca paso de fluido a vacio.
    cugauge::Interaction_GaugeZsurf(px,py,GaugeZmin,GaugeDp,MassLimit,czmin,czmax,Awen,Hdiv,ncells,cellmin,beginendcell,posxy,posz,code,velrhop,aux,DomPosMin,Scell,Fourh2,H,MassFluid);
    double zsurf=0;
    cudaMemcpy(&zsurf,aux,sizeof(double),cudaMemcpyDeviceToHost);
    CheckCudaError("CalculeGpu","Failed in free-surface calculation.");
    if(zsurf==DBL_MAX)zsurf=GaugeZmax;
    //-Stores result. | Guarda resultado.
    Zsurf=zsurf;
    //:Log->Printf("WD-> Zsurf:%f\n",Zsurf);
  }
  //SaveResult();
  return(Zsurf);
}
#endif


