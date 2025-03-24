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

/// \file JSphVRes.cpp \brief Implements the class \ref JSphVRes.

#include "JSphVRes.h"
#include "JSphCpu.h"
#include "JXml.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JSimpleNeigs.h"
#include "JTimeControl.h"
#include "JDsGaugeSystem.h"
#include "JCaseVRes.h"
#include "JSpVtkData.h"
#include "JSpVtkShape.h"
#include "JDataArrays.h"
#include "JDsVresData.h"
#include "JBoxDef.h"
#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JSphGpu_VRes_iker.h"
  #include "JDebugSphGpu.h"
#endif

#include <fstream>
#include <cstring>
#include <cfloat>
#include <climits>
#include <algorithm>
#include <fstream>
#include <iomanip>


using namespace std;

JSphVRes::JSphVRes(bool cpu,const StCteSph &csp,const JCaseVRes vreszone,unsigned zoneid
  ,tdouble3 maprealposmin,tdouble3 maprealposmax
  ,std::string appname,std::string dirdataout,unsigned partbegin,std::string partbegindir)
  : Log(AppInfo.LogPtr()),Cpu(cpu),CSP(csp),VresZone(vreszone),ZoneId(zoneid)
  ,MapRealPosMin(maprealposmin),MapRealPosMax(maprealposmax)
  ,AppName(appname),DirDataOut(dirdataout),PartBegin(partbegin),PartBeginDir(partbegindir)
{
  ClassName = "JSphVRes";
  BoxLimitMin=NULL;  BoxLimitMax=NULL; BoxDomMin=NULL;  BoxDomMax=NULL; 
  Width=NULL; Inner=NULL;  Tracking=NULL;
  NIni=NULL;  NPoints=NULL; Matmov=NULL;
  PtPointsIni=NULL;   PtNormalsIni=NULL; PtPoints=NULL; PtNormals=NULL;
  PtVelMot=NULL;  PtMass=NULL;
#ifdef _WITHGPU
  BoxLimitMing=NULL;  BoxLimitMaxg=NULL; BoxDomMing=NULL;  BoxDomMaxg=NULL; 
  Widthg=NULL; Innerg=NULL;  Trackingg=NULL;
  NInig=NULL;  NPointsg=NULL;  Matmovg=NULL; 
  PtPosxyg=NULL; PtPoszg=NULL; PtNormalsg=NULL; PtVelMotg=NULL;
  PtMassg=NULL;
#endif

  Reset();
  
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphVRes::~JSphVRes()
{
  DestructorActive = true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphVRes::Reset()
{

  // ZoneId=0;
  SvVResDataBi4=NULL;
  #ifdef _WITHGPU
  // GeomInfo=NULL;
  #endif
  ListSize=0;
  PtCount=0;
  for(int c=0;c<List.size();c++)delete List[c];
  List.clear();
  FreeMemory();
  FreePtMemory();
}

//==============================================================================
/// Allocates memory for vres zone configurations.
//==============================================================================
void JSphVRes::AllocateMemory(unsigned listsize){
  ListSize=listsize;
  try{
    BoxLimitMin   =new tdouble3 [ListSize];
    BoxLimitMax   =new tdouble3 [ListSize];
    BoxDomMin     =new tdouble3 [ListSize];
    BoxDomMax     =new tdouble3 [ListSize];
    Width         =new float    [ListSize];
    Inner         =new bool     [ListSize];
    Tracking      =new bool     [ListSize];
    NIni          =new unsigned [ListSize];
    NPoints       =new unsigned [ListSize];
    Matmov        =new JMatrix4d[ListSize];
    {
      memset(BoxLimitMin  ,0,sizeof(tdouble3)*ListSize);
      memset(BoxLimitMax  ,0,sizeof(tdouble3)*ListSize);
      memset(BoxDomMin    ,0,sizeof(tdouble3)*ListSize);
      memset(BoxDomMax    ,0,sizeof(tdouble3)*ListSize);
      memset(Width        ,0,sizeof(float   )*ListSize);
      memset(Inner        ,0,sizeof(bool    )*ListSize);
      memset(Tracking     ,0,sizeof(bool    )*ListSize);
      memset(NIni         ,0,sizeof(unsigned)*ListSize);
      memset(NPoints      ,0,sizeof(unsigned)*ListSize);
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
  #ifdef _WITHGPU
    if(!Cpu)AllocateMemoryGpu(ListSize);
  #endif
}

void JSphVRes::FreeMemory(){
  ListSize = 0;
  delete[]BoxLimitMin;  BoxLimitMin=NULL;
  delete[]BoxLimitMax;  BoxLimitMax=NULL;
  delete[]BoxDomMin;    BoxDomMin=NULL;
  delete[]BoxDomMax;    BoxDomMax=NULL;
  delete[]Width;        Width=NULL;
  delete[]Inner;        Inner=NULL;
  delete[]Tracking;     Tracking=NULL;
  delete[]NIni;         NIni=NULL;
  delete[]NPoints;      NPoints=NULL;
  delete[]Matmov;       Matmov=NULL;
  #ifdef _WITHGPU
    if(!Cpu)FreeMemoryGpu();
  #endif
}

//==============================================================================
/// Allocates memory for reference points.
//==============================================================================
void JSphVRes::AllocatePtMemory(unsigned ptcount){
  PtCount=ptcount;
  try{
    PtPointsIni     =new tdouble3   [ptcount];
    PtNormalsIni    =new tfloat3    [ptcount];
    PtPoints        =new tdouble3   [ptcount];
    PtNormals       =new tfloat3    [ptcount];
    PtVelMot        =new tfloat3    [ptcount];
    PtMass          =new float      [ptcount];
    {
      memset(PtPointsIni    ,0,sizeof(tdouble3)*PtCount);
      memset(PtNormalsIni   ,0,sizeof(tfloat3 )*PtCount);
      memset(PtPoints       ,0,sizeof(tdouble3)*PtCount);
      memset(PtNormals      ,0,sizeof(tfloat3 )*PtCount);
      memset(PtVelMot       ,0,sizeof(tfloat3 )*PtCount);
      memset(PtMass         ,0,sizeof(float   )*PtCount);
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
  #ifdef _WITHGPU
    if(!Cpu)AllocatePtMemoryGpu(PtCount);
  #endif
}

//==============================================================================
/// Frees allocated memory for interface points and auxiliary memory on GPU.
//==============================================================================
void JSphVRes::FreePtMemory(){
  delete[] PtPointsIni;     PtPointsIni=NULL;
  delete[] PtNormalsIni;    PtNormalsIni=NULL;
  delete[] PtPoints;        PtPoints=NULL;
  delete[] PtNormals;       PtNormals=NULL;
  delete[] PtVelMot;        PtVelMot=NULL;
  delete[] PtMass;          PtMass=NULL;
  #ifdef _WITHGPU
    if(!Cpu) FreePtMemoryGpu();
  #endif

}

#ifdef _WITHGPU
//==============================================================================
/// Allocates memory on GPU.
//==============================================================================
void JSphVRes::AllocateMemoryGpu(unsigned listsize){
 try{
    fcuda::Malloc(&BoxLimitMing   ,ListSize);
    fcuda::Malloc(&BoxLimitMaxg   ,ListSize);
    fcuda::Malloc(&BoxDomMing     ,ListSize);
    fcuda::Malloc(&BoxDomMaxg     ,ListSize);
    fcuda::Malloc(&Widthg         ,ListSize);
    fcuda::Malloc(&Innerg         ,ListSize);
    fcuda::Malloc(&Trackingg      ,ListSize);
    fcuda::Malloc(&NInig          ,ListSize);
    fcuda::Malloc(&NPointsg       ,ListSize);
    fcuda::Malloc(&Matmovg        ,ListSize);
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Frees allocated memory on GPU.
//==============================================================================
void JSphVRes::FreeMemoryGpu(){
  if(BoxLimitMing)    cudaFree(BoxLimitMing);   BoxLimitMing=NULL;
  if(BoxLimitMaxg)    cudaFree(BoxLimitMaxg);   BoxLimitMaxg=NULL;
  if(BoxDomMing)      cudaFree(BoxDomMing);     BoxDomMing=NULL;
  if(BoxDomMaxg)      cudaFree(BoxDomMaxg);     BoxDomMaxg=NULL;
  if(Widthg)          cudaFree(Widthg);         Widthg=NULL;
  if(Innerg)          cudaFree(Innerg);         Innerg=NULL;
  if(Trackingg)       cudaFree(Trackingg);      Trackingg=NULL;
  if(NInig)           cudaFree(NInig);          NInig=NULL;
  if(NPointsg)        cudaFree(NPointsg);       NPointsg=NULL; 
  if(Matmovg)         cudaFree(Matmovg);        Matmovg=NULL;
}

//==============================================================================
/// Allocates memory for interface points on GPU.
//==============================================================================
void JSphVRes::AllocatePtMemoryGpu(unsigned ptcount){
  fcuda::Malloc(&PtPosxyg   ,ptcount);
  fcuda::Malloc(&PtPoszg    ,ptcount);
  fcuda::Malloc(&PtNormalsg ,ptcount);
  fcuda::Malloc(&PtVelMotg  ,ptcount);
  fcuda::Malloc(&PtMassg    ,ptcount);
}

//==============================================================================
/// Frees allocated memory for reference points and auxiliary memory on GPU.
//==============================================================================
void JSphVRes::FreePtMemoryGpu(){
  if(PtPosxyg)     cudaFree(PtPosxyg);      PtPosxyg=NULL;
  if(PtPoszg)      cudaFree(PtPoszg);       PtPoszg=NULL;
  if(PtNormalsg)   cudaFree(PtNormalsg);    PtNormalsg=NULL;
  if(PtVelMotg)    cudaFree(PtVelMotg);     PtVelMotg=NULL;
  if(PtMassg)      cudaFree(PtMassg);       PtMassg=NULL;
}
#endif

//==============================================================================
/// Store vres configuration from XML.
//==============================================================================
void JSphVRes::CreateZones(){

  const JCaseVRes_Box* zone= VresZone.GetZoneBox(ZoneId);

  bool parent=(zone->Parent ? true : false);
  unsigned SubZones_Count=zone->Count();
  if(parent){
    JMatrix4d mat; mat.SetIdentity();
    bool isSimple=zone->GetBuffBox().IsSimple();
    if(!isSimple)GetRotMatrix(zone->GetPtBox(),mat,zone->GetPtBox().GetPosMin());
    unsigned zoneid           =zone->Parent->Id;
    tdouble3 boxlimitmininner =zone->GetPtBox().GetPosMin();
    tdouble3 boxlimitmaxinner =zone->GetPtBox().GetPosMax();
    tdouble3 boxlimitminouter =zone->GetBuffBox().GetPosMin(); 
    tdouble3 boxlimitmaxouter =zone->GetBuffBox().GetPosMax(); 
    tdouble3 boxlimitminmid   =zone->GetFixedBox().GetPosMin(); 
    tdouble3 boxlimitmaxmid   =zone->GetFixedBox().GetPosMax();
    bool trackingisactive     =zone->TrackingIsActive();
    unsigned trackingmk       =zone->GetTrackingMkBound();
    if(!isSimple)trackingisactive=true;
    JSphVResZone *zo = new JSphVResZone(Cpu, CSP,true ,zoneid, boxlimitmininner, boxlimitmaxinner
      ,boxlimitminouter, boxlimitmaxouter,boxlimitminmid,boxlimitmaxmid,MapRealPosMin,MapRealPosMax,trackingisactive,trackingmk,mat,isSimple);
      List.push_back(zo);
  }
  for(unsigned i=0; i<SubZones_Count; i++){
    const JCaseVRes_Box* subzone=zone->GetSubZoneBox(i);
    JMatrix4d mat; mat.SetIdentity();
    bool isSimple=subzone->GetParentBuffEndBox().IsSimple();
    if(!isSimple)GetRotMatrix(subzone->GetParentBuffIniBox(),mat,subzone->GetParentBuffIniBox().GetPosMin());
    unsigned zoneid             =subzone->Id;
    tdouble3 boxlimitmininner   =subzone->GetParentBuffEndBox().GetPosMin();
    tdouble3 boxlimitmaxinner   =subzone->GetParentBuffEndBox().GetPosMax();
    tdouble3 boxlimitminouter   =subzone->GetParentBuffIniBox().GetPosMin(); 
    tdouble3 boxlimitmaxouter   =subzone->GetParentBuffIniBox().GetPosMax();
    tdouble3 boxlimitminmid     =subzone->GetParentFixedBox().GetPosMin(); 
    tdouble3 boxlimitmaxmid     =subzone->GetParentFixedBox().GetPosMax(); 
    bool trackingisactive       =subzone->TrackingIsActive();
    unsigned trackingmk         =subzone->GetTrackingMkBound();
    if(!isSimple)trackingisactive=true;
    JSphVResZone *zo = new JSphVResZone(Cpu, CSP,false ,zoneid, boxlimitmininner, boxlimitmaxinner
      ,boxlimitminouter, boxlimitmaxouter,boxlimitminmid,boxlimitmaxmid,MapRealPosMin,MapRealPosMax,trackingisactive,trackingmk,mat,isSimple);
      List.push_back(zo);
  }
  parent=NULL;
}

//==============================================================================
/// Configures basic parameter of the simulation and prepares execution.
//==============================================================================
void JSphVRes::Config()
{
  //-Read XML file for vres zone configurations.
  CreateZones();

  //-Create objects for saving VRes data
  SvVResDataBi4 = new JDsVResDataSave(AppName,DirDataOut,Log);

  //-Allocates memory for vres zone configurations.
  AllocateMemory(GetCount());
  //-Prepares data for vres zone configurations.
  for(unsigned ci=0;ci<ListSize;ci++){
    BoxLimitMin   [ci]=List[ci]->GetBoxLimitMin();
    BoxLimitMax   [ci]=List[ci]->GetBoxLimitMax();
    BoxDomMin     [ci]=List[ci]->GetBoxDomMin();
    BoxDomMax     [ci]=List[ci]->GetBoxDomMax();
    Width         [ci]=List[ci]->GetWidth();
    Inner         [ci]=List[ci]->GetInner();;
    Tracking      [ci]=List[ci]->TrackingIsActive();
    NPoints       [ci]=List[ci]->getCount();
    Matmov        [ci]=JMatrix4d(List[ci]->GetMat());
  }

  unsigned npt=0;
  for(unsigned ci=0;ci<ListSize;ci++){
    NIni[ci]=npt;
    npt+=NPoints[ci];
  }

  //-Upload vres zone configurations on the gpu memory.
#ifdef _WITHGPU
  if(!Cpu){
    cudaMemcpy(BoxLimitMing ,BoxLimitMin  ,sizeof(double3)  *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(BoxLimitMaxg ,BoxLimitMax  ,sizeof(double3)  *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(BoxDomMing   ,BoxDomMin    ,sizeof(double3)  *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(BoxDomMaxg   ,BoxDomMax    ,sizeof(double3)  *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(Widthg       ,Width        ,sizeof(float)    *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(Innerg       ,Inner        ,sizeof(bool)     *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(Trackingg    ,Tracking     ,sizeof(bool)     *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(NInig        ,NIni         ,sizeof(unsigned) *ListSize,cudaMemcpyHostToDevice);
    cudaMemcpy(NPointsg     ,NPoints      ,sizeof(unsigned) *ListSize,cudaMemcpyHostToDevice);

    tmatrix4f* matc=  new tmatrix4f[ListSize];
    for(unsigned ci=0; ci<ListSize; ci++){
      tmatrix4f mat_new=Matmov[ci].GetMatrix4f();
      matc[ci]=mat_new;
    }
    cudaMemcpy(Matmovg,matc,sizeof(tmatrix4f)*ListSize,cudaMemcpyHostToDevice);
    delete[] matc;  matc=NULL;
  }
#endif
  

  //-Allocate memory on the cpu for interface points.
  AllocatePtMemory(npt);

  for(unsigned ci=0;ci<ListSize;ci++){
    npt=0;
    for(unsigned i=0;i<NPoints[ci];i++){
      PtPointsIni[NIni[ci]+i]  =List[ci]->getPoints()[i];
      PtNormalsIni[NIni[ci]+i] =TFloat3(static_cast<float>(List[ci]->getNormals()[i].x)
        , static_cast<float>(List[ci]->getNormals()[i].y), static_cast<float>(List[ci]->getNormals()[i].z));
    }
  }
  
  memcpy(PtPoints  ,PtPointsIni  ,sizeof(tdouble3)*PtCount);
  memcpy(PtNormals ,PtNormalsIni ,sizeof(tfloat3) *PtCount);

  //-Load data for restart.
  if(PartBegin)LoadVResData();

  //-Upload interface points info on the gpu memory.
#ifdef _WITHGPU
  if(!Cpu){
    tdouble2* pxy=new tdouble2[PtCount];
    double*   pz =new double[PtCount];
    for(unsigned c=0;c<PtCount;c++){
    pxy[c]=TDouble2(PtPoints[c].x,PtPoints[c].y);
    pz[c]=PtPoints[c].z;
    }
    cudaMemcpy(PtPosxyg     ,pxy          ,sizeof(double2)  *PtCount ,cudaMemcpyHostToDevice);
    cudaMemcpy(PtPoszg      ,pz           ,sizeof(double)   *PtCount ,cudaMemcpyHostToDevice);
    cudaMemcpy(PtNormalsg   ,PtNormals    ,sizeof(float3)   *PtCount ,cudaMemcpyHostToDevice);
    cudaMemcpy(PtVelMotg    ,PtVelMot     ,sizeof(float3)   *PtCount ,cudaMemcpyHostToDevice);
    cudaMemcpy(PtMassg      ,PtMass       ,sizeof(float)    *PtCount ,cudaMemcpyHostToDevice);
    delete[] pxy; pxy=NULL;
    delete[] pz;  pz=NULL;
  } 
#endif    
}

//========================================================================================
/// Obtain basic arrays for interface points for Cpu computations.
//=========================================================================================
StrDataVresCpu JSphVRes::GetZoneFluxInfoCpu(unsigned nzone){
  unsigned nini=NIni[nzone];
  unsigned npoints=NPoints[nzone];

  tdouble3*  points  =PtPoints   ;
  tfloat3*   normals =PtNormals  ;
  tfloat3*   velmot  =PtVelMot   ;
  float*     mass    =PtMass     ;
  return(StrDataVresCpu(npoints,nini,points,normals,velmot,mass));
}

StrGeomVresCpu JSphVRes::GetGeomInfoVresCpu(){
    return(StrGeomVresCpu(BoxLimitMin,BoxLimitMax,BoxDomMin,
      BoxDomMax,Tracking,Inner,Matmov));
  };

#ifdef _WITHGPU
//========================================================================================
/// Obtain basic arrays for interface points for Gpu computations.
//=========================================================================================
StrDataVresGpu JSphVRes::GetZoneFluxInfoGpu(unsigned nzone){
  unsigned nini=NIni[nzone];
  unsigned npoints=NPoints[nzone];

  double2*  posxyg  =PtPosxyg   ;
  double*   poszg   =PtPoszg    ;
  float3*   normals =PtNormalsg ;
  float3*   velmot  =PtVelMotg  ;
  float*    mass    =PtMassg    ;
  return(StrDataVresGpu(npoints,nini,posxyg,poszg,normals,velmot,mass));
}

  StrGeomVresGpu JSphVRes::GetGeomInfoVresGpu(){
    return(StrGeomVresGpu(BoxLimitMing,BoxLimitMaxg,BoxDomMing,
      BoxDomMaxg,Trackingg,Innerg,Matmovg));
  };
#endif

//========================================================================================
/// Update the movement of VRes region by matrix representation on Cpu and Gpu.
//=========================================================================================
void JSphVRes::UpdateMatMov(std::vector<JMatrix4d> mat){
  for(unsigned ci=0; ci<ListSize; ci++){
    JMatrix4d mat_new=mat[ci];    
    mat_new.Mul(Matmov[ci]);
    Matmov[ci]=mat_new;
  }
#ifdef _WITHGPU
  if(!Cpu){
    tmatrix4f* matc=  new tmatrix4f[ListSize];
    for(unsigned ci=0; ci<ListSize; ci++){
      tmatrix4f mat_1=Matmov[ci].GetMatrix4f();
      matc[ci]=mat_1;
    }
    cudaMemcpy(Matmovg,matc,sizeof(tmatrix4f)*ListSize,cudaMemcpyHostToDevice);
    delete[] matc;  matc=NULL;
  }
#endif
}

//========================================================================================
/// Save interface points arrays.
//=========================================================================================
void JSphVRes::SaveVResData(int part,double timestep,int nstep){

  SvVResDataBi4->InitPartData(part,timestep,nstep);

  unsigned msize=16*ListSize;
  double* matarray = new double[msize];
  for(unsigned i=0;i<ListSize;i++){
      tmatrix4d mat = Matmov[i].GetMatrix4d();
      double *ptr = (double*)&mat;
      for(unsigned j=0;j<16;j++){
        matarray[i*16+j]=ptr[j];
      }
  }

  #ifdef _WITHGPU
  if(!Cpu){
    cudaMemcpy(PtVelMot   ,PtVelMotg     ,sizeof(float3)   *PtCount ,cudaMemcpyDeviceToHost);
    cudaMemcpy(PtMass     ,PtMassg       ,sizeof(float)    *PtCount ,cudaMemcpyDeviceToHost);
  }
  #endif

  SvVResDataBi4->AddArray(PtCount,msize,PtPointsIni,PtNormalsIni,PtVelMot,PtMass,matarray);
  SvVResDataBi4->SavePartData();
  delete  []  matarray;   matarray=NULL;
}

//========================================================================================
/// Load interface points arrays for restaring simulation.
//=========================================================================================
void JSphVRes::LoadVResData(){

  JDsVResDataLoad edat(Log);
  edat.LoadPartData(PartBeginDir,int(PartBegin));
  

  unsigned msize=16*ListSize;
  double* matarray  = new double[msize];
  tdouble3* pos     = new tdouble3[PtCount];
  tfloat3*  normals = new tfloat3[PtCount];
  tfloat3*  velmot  = new tfloat3[PtCount];
  float*    mass    = new float  [PtCount];

  edat.LoadArray(PtCount,msize,velmot,mass,matarray);

  memcpy(PtPoints   ,pos      ,sizeof(tdouble3) *PtCount);
  memcpy(PtNormals  ,normals  ,sizeof(tfloat3)  *PtCount);
  memcpy(PtVelMot   ,normals  ,sizeof(tfloat3)  *PtCount);
  memcpy(PtMass     ,normals  ,sizeof(float)    *PtCount);

  for(unsigned i=0;i<ListSize;i++){
    tmatrix4d mat = reinterpret_cast<tmatrix4d*>(matarray)[i]; 
    Matmov[i]=JMatrix4d(mat);  
  }


  delete  []  pos     ;   pos     =NULL;
  delete  []  normals ;   normals =NULL;
  delete  []  velmot  ;   velmot  =NULL;
  delete  []  mass    ;   mass    =NULL;
  delete  []  matarray;   matarray=NULL;
}




unsigned JSphVRes::CreateListCpu(unsigned npf,unsigned pini,const tdouble3 *pos
  ,const unsigned *idp,typecode *code,int *inoutpart,unsigned nzone)
{
  unsigned count=0;
  const unsigned pfin=pini+npf;
  for(unsigned p=pini;p<pfin;p++){
    const typecode rcode=code[p];
    if(CODE_IsNormal(rcode) && CODE_IsFluid(rcode)){//-It includes only normal fluid particles (no periodic).
      if(CODE_IsFluidBuffer(rcode)){//-Particles already selected as Buffer.
        inoutpart[count]=p; count++;
      }
      else{//-Fluid particles no buffer.
        tdouble3 ps=pos[p];
        if(Tracking[nzone]) ps=MovePoint(ps,Matmov[nzone].GetMatrix4d());
        byte zone=255;
        if(List[nzone]->InZoneBox(ps))zone=byte(0);
          if(zone!=255){//-Particulas fluid que pasan a in/out.
            code[p]=CODE_ToFluidBuffer(rcode,zone)|CODE_TYPE_FLUID_BUFFERNUM; //-Adds 16 to indicate new particle in zone.
            inoutpart[count]=p; count++;
          }
         // if(List[nzone]->is_Out(pos[p]))code[p]=CODE_SetOutIgnore(rcode);
        }
      }
    }
  return(count);
}


unsigned JSphVRes::CreateListCpuInit(unsigned npf,unsigned pini,const tdouble3 *pos
  ,const unsigned *idp,typecode *code,int *inoutpart,unsigned nzone)
{
  unsigned count=0;
  const unsigned pfin=pini+npf;
  for(unsigned p=pini;p<pfin;p++){
    const typecode rcode=code[p];
    if(CODE_IsNormal(rcode) && CODE_IsFluid(rcode)){//-It includes only normal fluid particles (no periodic).
      if(CODE_IsFluidBuffer(rcode)){//-Particles already selected as Buffer.
        inoutpart[count]=p; count++;
      }else{//-Fluid particles no buffer.
        tdouble3 ps=pos[p];
        if(Tracking[nzone]) ps=MovePoint(ps,Matmov[nzone].GetMatrix4d());
        byte zone=255;
        if(List[nzone]->InZoneBox(ps))zone=byte(0);
        if(zone!=255){//-Particulas fluid que pasan a in/out.
          code[p]=CODE_ToFluidBuffer(rcode,zone)|CODE_TYPE_FLUID_BUFFERNUM; //-Adds 16 to indicate new particle in zone.
          inoutpart[count]=p; count++;
        }
        if(List[nzone]->is_Out(ps)){
          if(!List[nzone]->InZoneBoxMid(ps))code[p]=CODE_TYPE_FLUID_FIXED;
          else code[p]=CODE_SetOutIgnore(rcode);
        }
      }
    }
  }
  return(count);
}



void JSphVRes::CheckNormals(TpBoundary tboundary,unsigned np
  ,unsigned pini,const tdouble3 *pos,const unsigned *idp,typecode *code
  ,const tfloat3* boundnor,unsigned vresid)
{
  unsigned nerr=0;
  byte* errpart=new byte[np];
  memset(errpart,0,sizeof(byte)*np);
  const unsigned nzone=static_cast<unsigned>(List.size());
  const unsigned pfin=pini+np;
  for(unsigned p=pini;p<pfin;p++){
    const typecode rcode=code[p];
    tdouble3 ps=pos[p]; 
    for(unsigned zone=0;zone<nzone;zone++){
      if(Tracking[zone]) ps=MovePoint(ps,Matmov[zone].GetMatrix4d());
      if(List[zone]->InZoneBox(ps)){
        if(tboundary==BC_DBC){
          errpart[p]=1; nerr++;
        }else{
          if(boundnor[p]==TFloat3(0)){
            errpart[p]=1; nerr++;
          } 
        }
      }
    }
  }
  if(nerr>0){
    tfloat3* vpos=new tfloat3[nerr];
    byte*    vtype=new byte[nerr];
    unsigned pp=0;
    for(unsigned p=0;p<np;p++)if(errpart[p]){
      vpos[pp]=ToTFloat3(pos[p]);
      vtype[pp]=errpart[p];
      pp++;
    }
    JDataArrays arrays;
    arrays.AddArray("Pos",nerr,vpos,false);
    arrays.AddArray("ErrorType",nerr,vtype,false);
    const std::string vrname =std::string("Cfg")+fun::PrintStr("_vres%02u", vresid)+std::string("_ErrorParticles.vtk");
    const string filevtk=AppInfo.GetDirOut()+vrname;
    JSpVtkData::Save(filevtk,arrays,"Pos");
    Log->AddFileInfo(filevtk,"Saves boundary particles inside VRes buffer region.");
    delete[] vpos;  vpos=NULL;
    delete[] vtype; vtype=NULL;
  }
  if(nerr>0) Run_Exceptioon("Boundary particles in the buffer region are not allowed. Check VTK file Cfg_vres%02_ErrorParticles.vtk with excluded particles.");
}

#ifdef _WITHGPU
unsigned JSphVRes::CreateListGpuInit(unsigned npf,unsigned pini,const double2 *posxyg
  ,const double *poszg,typecode *codeg,unsigned size,int *inoutpartg,unsigned nzone)
{
  unsigned count = 0;
  tmatrix4f mat=Matmov[nzone].GetMatrix4f();
  tmatrix4f matc=mat;

  const bool Stable=false;
  count = cusphvres::BufferCreateListInit(Stable,npf,pini
    ,List[nzone]->GetBoxLimitMinInner(),List[nzone]->GetBoxLimitMaxInner()
    ,List[nzone]->GetBoxLimitMinOuter(),List[nzone]->GetBoxLimitMaxOuter()
    ,List[nzone]->GetBoxLimitMinMid(),List[nzone]->GetBoxLimitMaxMid()
    ,List[nzone]->GetInner(),posxyg,poszg,codeg,(unsigned *)inoutpartg
    ,matc,Tracking[nzone],nzone);
  return (count);
}
unsigned JSphVRes::CreateListGpu(unsigned npf,unsigned pini,const double2 *posxyg
  ,const double *poszg,typecode *codeg,unsigned size,int *inoutpartg,unsigned nzone)
{
  unsigned count = 0;
  const bool Stable=false;
  count = cusphvres::BufferCreateList(Stable,npf,pini
    ,List[nzone]->GetBoxLimitMinInner(),List[nzone]->GetBoxLimitMaxInner()
    ,List[nzone]->GetBoxLimitMinOuter(),List[nzone]->GetBoxLimitMaxOuter()
    ,List[nzone]->GetInner(),posxyg,poszg,codeg,(unsigned *)inoutpartg
    ,Matmovg,Tracking[nzone],nzone);
  return (count);
}


// unsigned JSphVRes::CheckNormals(TpBoundary tboundary,unsigned npf,unsigned pini
//   ,const double2 *posxyg,const double *poszg,const float3* boundnor
//   ,typecode *codeg,unsigned size,int *inoutpartg,unsigned nzone)
// {
//   unsigned count = 0;
//   bool Stable = false;
//   count = cusphvres::BufferCheckNormals(tboundary,Stable,npf,pini,List[nzone]->BoxLimitMinInner
//     ,List[nzone]->BoxLimitMaxInner,List[nzone]->BoxLimitMinOuter
//     ,List[nzone]->BoxLimitMaxOuter, List[nzone]->Inner,posxyg,poszg,boundnor
//     ,codeg,(unsigned *)inoutpartg,Matmovg,Tracking[nzone],nzone);

//   return (count);
// }
#endif

tdouble3 JSphVRes::MovePoint(tdouble3 oldpos,const tmatrix4d& mat){
  tdouble3 newpos=TDouble3(0);
  oldpos.x-=mat.a14;  oldpos.y-=mat.a24;  oldpos.z-=mat.a34;
  newpos.x=oldpos.x*mat.a11+oldpos.y*mat.a21+oldpos.z*mat.a31;
  newpos.y=oldpos.x*mat.a12+oldpos.y*mat.a22+oldpos.z*mat.a32;
  newpos.z=oldpos.x*mat.a13+oldpos.y*mat.a23+oldpos.z*mat.a33;
  return(newpos);
}

//==============================================================================
/// ComputeStep over buffer regions:
/// - If buffer particle is moved to fluid zone then it changes to fluid particle.
/// - If fluid particle is moved to buffer zone then it changes to buffer particle.
/// - If buffer particle is moved out the domain then it changes to ignore particle.
/// - Compute number of new particle to be created.
//==============================================================================
unsigned JSphVRes::ComputeStepCpu(unsigned bufferpartcount,int *bufferpart
	,typecode *code,const tdouble3 *pos,unsigned nzone){
  //-Updates code according to particle position and define new particles to create.
	const int ncp=int(bufferpartcount);
	#ifdef OMP_USE
	  #pragma omp parallel for schedule (static)
	#endif
	for(int cp=0;cp<ncp;cp++){
		typecode cod=0;
	  byte newiz=255;
	  const unsigned p=(unsigned)bufferpart[cp];
	  const typecode rcode=code[p];
	  const unsigned izone0=CODE_GetIzoneFluidBuffer(rcode);
	  tdouble3 ps=pos[p]; 
    if(Tracking[nzone]) ps=MovePoint(ps,Matmov[nzone].GetMatrix4d());

	  if(izone0>=16){//-Normal fluid particle in buffer zone.
		  cod= rcode^0x10 ; //-Converts to buffer particle or not.
		  code[p]=cod;
	  }else{//-Previous inout fluid particle.
			const bool isnormal=List[nzone]->is_Normal(ps);
			if(List[nzone]->is_Out(ps)){
			  cod=CODE_SetOutIgnore(rcode); //-Particle is moved out domain.
				code[p]=cod;
			}else if(isnormal){
				cod=CODE_TYPE_FLUID; //-Particle become normal;
			  code[p]=cod;
			}
		}
	}
  //-Compute number of new particles.
	unsigned newcp=0;
  const unsigned nini=NIni[nzone];
  const unsigned npoints=NPoints[nzone];
  for(unsigned i=nini;i<npoints;i++){
    float mass=PtMass[i];
    if(mass>CSP.massfluid) newcp++;
  }
	return(newcp);
}

//==============================================================================
/// Create new buffer particles at Dp*0.5 in the normal direction to the interface.
//==============================================================================
void JSphVRes::CreateNewPart(const unsigned idnext,unsigned *dcell,typecode *code,tdouble3 *pos,unsigned *idp,
	tfloat4 *velrhop,const JSphCpu *sphcpu,unsigned np,unsigned nzone){    
  unsigned newcp=0;
  const unsigned nini=NIni[nzone];
  const unsigned npoints=NPoints[nzone];
  for(unsigned i=nini;i<npoints;i++){
	  float mass=PtMass[i];
    if(mass>CSP.massfluid){
      tdouble3 rpos=PtPoints[i];
      tfloat3  normal=PtNormals[i];
      PtMass[i]-=CSP.massfluid;
      const double dis=0.5*CSP.dp;
      rpos.x-=(double)dis*normal.x;
      rpos.y-=(double)dis*normal.y;
      rpos.z-=(double)dis*normal.z;
      const unsigned p2=np+newcp;
      code[p2]=CODE_ToFluidBuffer(CODE_TYPE_FLUID,byte(nzone));
			sphcpu->UpdatePos(rpos,0,0,0,false,p2,pos,dcell,code);
		  idp[p2]=idnext+newcp;
		  velrhop[p2]=TFloat4(0,0,0,1000);
		  newcp++;
    }
  }
}

#ifdef _WITHGPU
//==============================================================================
/// ComputeStep over buffer regions:
/// - If buffer particle is moved to fluid zone then it changes to fluid particle.
/// - If fluid particle is moved to buffer zone then it changes to buffer particle.
/// - If buffer particle is moved out the domain then it changes to ignore particle.
//==============================================================================
void JSphVRes::ComputeStepGpu(unsigned bufferpartcount,int *bufferpart,unsigned idnext,double2 *posxyg,double *poszg,unsigned *dcellg,typecode *codeg
    ,unsigned *idpg,float4 *velrhopg,byte *newizoneg,const JSphGpuSingle *gp,unsigned nzone)
{
  //-Checks particle position.
  cusphvres::BufferComputeStep(bufferpartcount,bufferpart,posxyg,poszg,codeg
    ,List[nzone]->GetBoxLimitMinInner(),List[nzone]->GetBoxLimitMaxInner()
    ,List[nzone]->GetBoxLimitMinOuter(),List[nzone]->GetBoxLimitMaxOuter()
    ,List[nzone]->GetInner(),Matmovg,Tracking[nzone],nzone);

}

//==============================================================================
/// Compute list of new buffer particles.
//==============================================================================
unsigned JSphVRes::NewPartListCreate(int* newpart,unsigned idzone){
  unsigned newnp=0;
  newnp = cusphvres::NewPartListCreate(NPoints[idzone],NIni[idzone]
    ,NPoints[idzone],PtMassg,newpart,CSP.massfluid);

  return(newnp);
}

//==============================================================================
/// Create new buffer particles at Dp*0.5 in the normal direction to the interface.
//==============================================================================
void JSphVRes::CreateNewPartGpu(unsigned np,unsigned newnp,unsigned idnext,double2 *posxy
  ,double *posz,unsigned *dcell,typecode *code,unsigned *idp,float4 *velrhop
  ,int *newpart,unsigned nzone)
{
  cusphvres::CreateNewPart(newnp,NIni[nzone],newpart,np,idnext
    ,posxy,posz,dcell,code,idp,velrhop,PtNormalsg, static_cast<float>(CSP.dp),PtPosxyg,PtPoszg,nzone);
}
#endif

//========================================================================================
/// Move position and orient normal for accumulation points on the VRes interface.
//=========================================================================================
void JSphVRes::MoveBufferZone(double dt,std::vector<JMatrix4d> mat){
  for(int i=0; i<List.size(); i++){
    unsigned ntot=NPoints[i];
    unsigned nini=NIni[i];
    tmatrix4d mat_i=mat[i].GetMatrix();
    if(Tracking[i]){
      if(Cpu){
        for(int p=nini;p<int(ntot);p++){
          const tdouble3 posp=PtPoints[p];
          const tdouble3 normalp=ToTDouble3(PtNormals[p]);
          const tdouble3 pospnew=mat[i].MulPoint(posp);
          const tfloat3 normalpnew=ToTFloat3(mat[i].MulNormal(normalp));
          if(dt>0) PtVelMot[p]=ToTFloat3(pospnew-posp)/ static_cast<float>(dt);
          PtPoints[p]=pospnew;
          PtNormals[p]=normalpnew;
        }
      }else{
        #ifdef _WITHGPU
        cusphvres::MoveBufferZone(NIni[i],ntot,PtPosxyg,PtPoszg,PtNormalsg,PtVelMotg,dt,mat_i,i);
        #endif
      }
    }
  }
}

//==============================================================================
// Saves VTK file with accumulation points info.
//==============================================================================
void JSphVRes::SaveNormals(std::string filename,int numfile){
  #ifdef _WITHGPU
  if(!Cpu){
    tdouble3 *posh=fcuda::ToHostPosd3(0,PtCount,PtPosxyg,PtPoszg);
    memcpy(PtPoints,posh,sizeof(tdouble3)*PtCount);
    delete[] posh;

    cudaMemcpy(PtNormals, PtNormalsg, sizeof(float3) * PtCount, cudaMemcpyDeviceToHost);
    cudaMemcpy(PtVelMot,  PtVelMotg,  sizeof(float3) * PtCount, cudaMemcpyDeviceToHost);
    cudaMemcpy(PtMass,    PtMassg,    sizeof(float)  * PtCount, cudaMemcpyDeviceToHost);
  }
  #endif

  SaveVtkNormals(filename,numfile,PtCount,PtPoints,PtNormals,PtVelMot,PtMass);
}

//==============================================================================
// Saves VTK file with accumulation points info.
//==============================================================================
void JSphVRes::SaveVtkNormals(std::string filename,int numfile,unsigned np
  ,const tdouble3* pos,const tfloat3* boundnor,const tfloat3* velflux,const float* flux)const
{
    if(fun::GetExtension(filename).empty())filename=fun::AddExtension(filename,"vtk");
    if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
    //-Find floating particles.
       //-Allocate memory for boundary particles.
    const unsigned npsel=np;
    JDataArrays arrays;
    tdouble3* vpos =arrays.CreateArrayPtrDouble3("Pos",npsel);
    tfloat3*  vnor =arrays.CreateArrayPtrFloat3 ("Normal",npsel);
    tfloat3* vflux =arrays.CreateArrayPtrFloat3("VelFlux",npsel);
    float* vmass =arrays.CreateArrayPtrFloat("Flux",npsel);

    //-Loads data of fixed and moving particles.
    memcpy(vpos,pos,sizeof(tdouble3)*npsel);
    memcpy(vnor,boundnor,sizeof(tfloat3)*npsel);
    memcpy(vflux,velflux,sizeof(tfloat3)*npsel);
    memcpy(vmass,flux,sizeof(float)*npsel);
    //-Saves VTK file.
    JSpVtkData::Save(filename,arrays,"Pos");
    //-Frees memory.
    arrays.Reset();
}


void JSphVRes::GetQuadPoints2d(tdouble3 pmin,tdouble3 pmax,tdouble3* vpt)const{
  const tdouble3 s=pmax-pmin;
  vpt[0]=pmin;
  vpt[1]=TDouble3(pmin.x+s.x,pmin.y,pmin.z);
  vpt[2]=TDouble3(pmin.x+s.x,pmin.y,pmin.z+s.z);
  vpt[3]=TDouble3(pmin.x    ,pmin.y,pmin.z+s.z);
}

//==============================================================================
// Saves VTK file with domain of zones.
//==============================================================================
void JSphVRes::SaveVtkDomains(std::string filename,int numfile,bool is2d)const{
  string file=fun::FileNameSec(filename,numfile);
  const unsigned nz=ListSize;
  JSpVtkShape sh;
  if(is2d){
    for(unsigned id=0;id<nz;id++){
      tdouble3 pt[4];
      GetQuadPoints2d(BoxLimitMin[id],BoxLimitMax[id],pt);
      sh.AddQuadWire(Matmov[id].MulPoint(pt[0]),Matmov[id].MulPoint(pt[1])
        ,Matmov[id].MulPoint(pt[2]),Matmov[id].MulPoint(pt[3]),id);
    }
    sh.SaveVtk(file+".vtk","Zone");
  }
  else{
    for(unsigned id=0;id<nz;id++){
      const tdouble3 p0=Matmov[id].MulPoint(BoxLimitMin[id]),s0=Matmov[id].MulPoint(BoxLimitMax[id])-p0;
      const tdouble3 size=BoxLimitMax[id]-BoxLimitMin[id];
      const tdouble3 vx=Matmov[id].MulNormal(TDouble3(size.x,0,0));
      const tdouble3 vy=Matmov[id].MulNormal(TDouble3(0,size.y,0));
      const tdouble3 vz=Matmov[id].MulNormal(TDouble3(0,0,size.z));

      sh.AddBoxSizeVec(p0,vx,vy,vz,id);
    }     
    sh.SaveVtk(file+".vtk","Zone");
  }
}

void JSphVRes::GetRotMatrix(const JBoxDef& boxdef,JMatrix4d& mat,const tdouble3 posmin){
  tmatrix4d rotmat=TMatrix4d();
  tdouble3 vx=boxdef.GetVx();
  tdouble3 vy=boxdef.GetVy();
  tdouble3 vz=boxdef.GetVz();


  double norm_vx=sqrt(vx.x*vx.x+vx.y*vx.y+vx.z*vx.z);
  double norm_vy=sqrt(vy.x*vy.x+vy.y*vy.y+vy.z*vy.z);
  double norm_vz=sqrt(vz.x*vz.x+vz.y*vz.y+vz.z*vz.z);
  rotmat.a11=vx.x/norm_vx; rotmat.a21=vx.y/norm_vx; rotmat.a31=vx.z/norm_vx;
  rotmat.a13=vz.x/norm_vz; rotmat.a23=vz.y/norm_vz; rotmat.a33=vz.z/norm_vz;
  if(!boxdef.GetIs2D()) {
    rotmat.a12=vy.x/norm_vy; rotmat.a22=vy.y/norm_vy; rotmat.a32=vy.z/norm_vy;
  }
  if(boxdef.GetIs2D()){
    rotmat.a21=0.0;
    rotmat.a23=0.0;
  }
  rotmat.a14=(1.0-rotmat.a11)*posmin.x-rotmat.a12*posmin.y-rotmat.a13*posmin.z;
  rotmat.a24=-rotmat.a21*posmin.x+(1.0-rotmat.a22)*posmin.y-rotmat.a23*posmin.z;
  rotmat.a34=-rotmat.a31*posmin.x-rotmat.a32*posmin.y+(1.0-rotmat.a33)*posmin.z;

  JMatrix4d rotmatj(rotmat);
  mat.MulPre(rotmatj);
}
