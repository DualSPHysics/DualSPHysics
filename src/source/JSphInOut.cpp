//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphInOut.cpp \brief Implements the class \ref JSphInOut.

#include "JSphInOut.h"
#include "JSphInOutZone.h"
#include "JSphInOutGridData.h"
#include "JSphCpu.h"
#include "JXml.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsGeo3d.h"
#include "JVtkLib.h"
#include "JSimpleNeigs.h"
#include "JTimeControl.h"
#include "JDsGaugeSystem.h"
#include "JNumexLib.h"

#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JSphGpu_InOut_iker.h"
  #include "JDebugSphGpu.h"
#endif

#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSphInOut
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOut::JSphInOut(bool cpu,JLog2 *log,std::string xmlfile,JXml *sxml,std::string xmlpath,const std::string &dirdatafile)
  :Cpu(cpu),Log(log),XmlFile(xmlfile),XmlPath(xmlpath),DirDataFile(dirdatafile)
{
  ClassName="JSphInOut";
  Planes=NULL;
  CfgZone=NULL; CfgUpdate=NULL;  Width=NULL;  DirData=NULL;  VelData=NULL;  Zbottom=NULL; Zsurf=NULL;
  PtZone=NULL;  PtPos=NULL;  PtAuxDist=NULL;
  #ifdef _WITHGPU
    Planesg=NULL; BoxLimitg=NULL; CfgZoneg=NULL; CfgUpdateg=NULL; Widthg=NULL; DirDatag=NULL; Zsurfg=NULL;
    PtZoneg=NULL; PtPosxyg=NULL; PtPoszg=NULL; PtAuxDistg=NULL;
  #endif
  Reset();
  //-Loads basic configuration.
  LoadXmlInit(sxml,xmlpath);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOut::~JSphInOut(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInOut::Reset(){

  Stable=Simulate2D=false;
  Simulate2DPosY=0;
  PeriActive=0;
  RhopZero=CteB=Gamma=GravityZ=CoefHydro=0;
  Dp=0;
  MapRealPosMin=MapRealPosMax=TDouble3(0);
  CodeNewPart=0;

  ReuseIds=false;
  MemoryResize0=MemoryResize1=0;
  NpResizePlus0=NpResizePlus1=0;

  DetermLimit=0;
  ExtrapolateMode=0;

  UseBoxLimit=true;
  FreeCentre=TFloat3(FLT_MAX);
  FreeLimitMin=TFloat3(-FLT_MAX);
  FreeLimitMax=TFloat3(FLT_MAX);

  NewNpTotal=0;
  NewNpPart=0;
  CurrentNp=0;

  RefillAdvanced=false;
  VariableVel=false;
  VariableZsurf=false;
  CalculatedZsurf=false;
  ExtrapolatedData=NoExtrapolatedData=false;
  InterpolatedVel=false;
  for(int c=0;c<List.size();c++)delete List[c];
  List.clear();
  FreeMemory();

  FreePtMemory();
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSphInOut::LoadXmlInit(const JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNodeSimple(place,true);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  TiXmlElement* ele=node->ToElement();
  //-Checks old configuration.
  if(sxml->ExistsElement(ele,"userefilling"))Run_ExceptioonFile("Inlet/outlet: <userefilling> is not supported by current inlet/outlet version.",sxml->ErrGetFileRow(ele,"userefilling"));
  //-Checks element names.
  sxml->CheckElementNames(ele,true,"memoryresize useboxlimit determlimit extrapolatemode *inoutzone");
  //-Loads general configuration.
  ReuseIds=false; //sxml->GetAttributeBool(ele,"reuseids",true,false);//-Since ReuseIds=true is not implemented.
  //-Loads configuration for memoryresize.
  MemoryResize0=sxml->ReadElementFloat(ele,"memoryresize","size0",true,2.f);
  MemoryResize1=sxml->ReadElementFloat(ele,"memoryresize","size" ,true,4.f);
  if(MemoryResize0<0)Run_ExceptioonFile("Value of memoryresize.size0 lower than zero is invalid.",sxml->ErrGetFileRow(ele,"memoryresize"));
  if(MemoryResize1<0)Run_ExceptioonFile("Value of memoryresize.size lower than zero is invalid.",sxml->ErrGetFileRow(ele,"memoryresize"));
  //-Loads value determlimit.
  DetermLimit=sxml->ReadElementFloat(ele,"determlimit","value",true,1e+3f);
  //-Loads ExtrapolateMode.
  ExtrapolateMode=sxml->ReadElementInt(ele,"extrapolatemode","value",true,1);
  if(ExtrapolateMode>3)ExtrapolateMode=3;
  if(ExtrapolateMode<1)ExtrapolateMode=1;
  if(ExtrapolateMode<2 && Cpu)ExtrapolateMode=2;
  //-Loads UseBoxLimit.
  UseBoxLimit=sxml->ReadElementBool(ele,"useboxlimit","value",true,true);
  {
    TiXmlElement* ele2=ele->FirstChildElement("useboxlimit"); 
    if(ele2 && sxml->ExistsElement(ele2,"freecentre"))FreeCentre=sxml->ReadElementFloat3(ele2,"freecentre");
    else FreeCentre=TFloat3(FLT_MAX);
  }
  //-Loads MkFluidList.
  MkFluidList.clear();
  TiXmlElement* ele2=ele->FirstChildElement("inoutzone"); 
  while(ele2){
    if(sxml->CheckElementActive(ele2)){
      TiXmlElement* zone=ele2->FirstChildElement("zone2d");
      if(!zone)zone=ele2->FirstChildElement("zone3d");
      if(zone){
        const unsigned mkfluid=sxml->ReadElementUnsigned(zone,"particles","mkfluid",true,UINT_MAX);
        if(mkfluid!=UINT_MAX){
          unsigned c=0;
          for(;c<unsigned(MkFluidList.size()) && mkfluid!=MkFluidList[c];c++);
          if(c<unsigned(MkFluidList.size()))Run_Exceptioon(fun::PrintStr("Mkfluid=%u is used in several <inoutzone> definitions.",mkfluid));
          MkFluidList.push_back(mkfluid);
        }
      }
    }
    ele2=ele2->NextSiblingElement("inoutzone");
  }
}

//==============================================================================
/// Loads data of a file in XML format.
//==============================================================================
void JSphInOut::LoadFileXml(const std::string &file,const std::string &path
  ,JNumexLib *nuxlib,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem)
{
  JXml jxml;
  jxml.LoadFile(file);
  jxml.SetNuxLib(nuxlib); //-Enables the use of NuxLib in XML configuration.
  LoadXml(&jxml,path,partsdata,gaugesystem);
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSphInOut::LoadXml(const JXml *sxml,const std::string &place
  ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem)
{
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement(),partsdata,gaugesystem);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JSphInOut::ReadXml(const JXml *sxml,TiXmlElement* lis
  ,const JDsPartsInit *partsdata,JGaugeSystem *gaugesystem)
{
  //-Loads inflow elements.
  const unsigned idmax=CODE_MASKTYPEVALUE-CODE_TYPE_FLUID_INOUT;
  if(idmax-CODE_TYPE_FLUID_INOUTNUM!=MaxZones-1)Run_Exceptioon("Maximum number of inlet/outlet zones is invalid.");
  TiXmlElement* ele=lis->FirstChildElement("inoutzone"); 
  while(ele){
    if(sxml->CheckElementActive(ele)){
      const unsigned id=GetCount();
      if(id>idmax)Run_Exceptioon("Maximum number of inlet/outlet zones has been reached.");
      JSphInOutZone* iet=new JSphInOutZone(Cpu,Log,id,Simulate2D,Simulate2DPosY,Dp
        ,MapRealPosMin,MapRealPosMax,GravityZ,sxml,ele,DirDataFile,partsdata,gaugesystem);
      List.push_back(iet);
    }
    ele=ele->NextSiblingElement("inoutzone");
  }
}

//==============================================================================
/// Allocates memory for inlet/outlet configurations.
//==============================================================================
void JSphInOut::AllocateMemory(unsigned listsize){
  ListSize=listsize;
  try{
    const unsigned size=32;
    Planes   =new tplane3f[size];
    CfgZone  =new byte    [size];
    CfgUpdate=new byte    [size];
    Width    =new float   [size];
    DirData  =new tfloat3 [size];
    VelData  =new tfloat4 [size*2];
    Zbottom  =new float   [size];
    Zsurf    =new float   [size];
    if(1){
      memset(Planes   ,255,sizeof(tplane3f)*size);
      memset(CfgZone  ,255,sizeof(byte    )*size);
      memset(CfgUpdate,255,sizeof(byte    )*size);
      memset(Width    ,255,sizeof(float   )*size);
      memset(DirData  ,255,sizeof(tfloat3 )*size);
      memset(VelData  ,255,sizeof(tfloat4 )*size*2);
      memset(Zbottom  ,255,sizeof(float   )*size);
      memset(Zsurf    ,255,sizeof(float   )*size);
    }
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
  #ifdef _WITHGPU
    if(!Cpu)AllocateMemoryGpu(ListSize);
  #endif
}

//==============================================================================
/// Frees allocated memory.
//==============================================================================
void JSphInOut::FreeMemory(){
  ListSize=0;
  delete[] Planes;    Planes=NULL;
  delete[] CfgZone;   CfgZone=NULL;
  delete[] CfgUpdate; CfgUpdate=NULL;
  delete[] Width;     Width=NULL;
  delete[] DirData;   DirData=NULL;
  delete[] VelData;   VelData=NULL;
  delete[] Zbottom;   Zbottom=NULL;
  delete[] Zsurf;     Zsurf=NULL;
  #ifdef _WITHGPU
    if(!Cpu)FreeMemoryGpu();
  #endif
}

//==============================================================================
/// Allocates memory for reference points.
//==============================================================================
void JSphInOut::AllocatePtMemory(unsigned ptcount){
  PtCount=ptcount;
  try{
    PtZone   =new byte    [ptcount];
    PtPos    =new tdouble3[ptcount];
    PtAuxDist=new float   [ptcount];
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
  #ifdef _WITHGPU
    if(!Cpu)AllocatePtMemoryGpu(PtCount);
  #endif
}

//==============================================================================
/// Frees allocated memory for reference points and auxiliary memory.
//==============================================================================
void JSphInOut::FreePtMemory(){
  PtCount=0;
  delete[] PtZone;    PtZone=NULL;
  delete[] PtPos;     PtPos=NULL;
  delete[] PtAuxDist; PtAuxDist=NULL;
  #ifdef _WITHGPU
    if(!Cpu)FreePtMemoryGpu();
  #endif
}

#ifdef _WITHGPU
//==============================================================================
/// Allocates memory for reference points on GPU.
//==============================================================================
void JSphInOut::AllocatePtMemoryGpu(unsigned ptcount){
  fcuda::Malloc(&PtZoneg,ptcount);
  fcuda::Malloc(&PtPosxyg,ptcount);
  fcuda::Malloc(&PtPoszg,ptcount);
  fcuda::Malloc(&PtAuxDistg,ptcount);
}

//==============================================================================
/// Frees allocated memory for reference points and auxiliary memory on GPU.
//==============================================================================
void JSphInOut::FreePtMemoryGpu(){
  if(PtZoneg)   cudaFree(PtZoneg);    PtZoneg=NULL;
  if(PtPosxyg)  cudaFree(PtPosxyg);   PtPosxyg=NULL;
  if(PtPoszg)   cudaFree(PtPoszg);    PtPoszg=NULL;
  if(PtAuxDistg)cudaFree(PtAuxDistg); PtAuxDistg=NULL;
}

//==============================================================================
/// Allocates memory on GPU for inlet/outlet configurations.
//==============================================================================
void JSphInOut::AllocateMemoryGpu(unsigned listsize){
  try{
    fcuda::Malloc(&Planesg   ,ListSize);
    fcuda::Malloc(&BoxLimitg ,ListSize*3);
    fcuda::Malloc(&CfgZoneg  ,ListSize);
    fcuda::Malloc(&CfgUpdateg,ListSize);
    fcuda::Malloc(&Widthg    ,ListSize);
    fcuda::Malloc(&DirDatag  ,ListSize);
    fcuda::Malloc(&Zsurfg    ,ListSize);
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
}

//==============================================================================
/// Frees allocated memory on GPU.
//==============================================================================
void JSphInOut::FreeMemoryGpu(){
  if(Planesg)   cudaFree(Planesg);     Planesg=NULL;
  if(BoxLimitg) cudaFree(BoxLimitg);   BoxLimitg=NULL;
  if(CfgZoneg)  cudaFree(CfgZoneg);    CfgZoneg=NULL;
  if(CfgUpdateg)cudaFree(CfgUpdateg);  CfgUpdateg=NULL;
  if(Widthg)    cudaFree(Widthg);      Widthg=NULL;
  if(DirDatag)  cudaFree(DirDatag);    DirDatag=NULL;
  if(Zsurfg)    cudaFree(Zsurfg);      Zsurfg=NULL;
}
#endif


//==============================================================================
/// Saves domain zones in VTK file.
//==============================================================================
void JSphInOut::ComputeFreeDomain(){
  //-Recalculates FreeCentre starting from domain simulation.
  if(FreeCentre==TFloat3(FLT_MAX))FreeCentre=ToTFloat3((MapRealPosMax+MapRealPosMin)/2);
  //Log->Printf("----> FreeCentre:(%g,%g,%g)",FreeCentre.x,FreeCentre.y,FreeCentre.z);
  //-Calculates domain starting from FreeCentre for each axis.
  tfloat3 pmin=ToTFloat3(MapRealPosMin),pmax=ToTFloat3(MapRealPosMax);
  const unsigned nzone=GetCount();
  tfloat3 *zolimit=new tfloat3[nzone];
  byte *sel=new byte[nzone];
  tfloat3 *dmin=new tfloat3[nzone+1];
  tfloat3 *dmax=new tfloat3[nzone+1];
  for(unsigned ci=0;ci<nzone;ci++){ 
    sel[ci]=0;
    const tfloat3 boxmin=List[ci]->GetBoxLimitMin();
    const tfloat3 boxmax=List[ci]->GetBoxLimitMax();
    zolimit[ci].x=(FreeCentre.x<=boxmin.x? boxmin.x: (FreeCentre.x>=boxmax.x? boxmax.x: FLT_MAX));
    zolimit[ci].y=(FreeCentre.y<=boxmin.y? boxmin.y: (FreeCentre.y>=boxmax.y? boxmax.y: FLT_MAX));
    zolimit[ci].z=(FreeCentre.z<=boxmin.z? boxmin.z: (FreeCentre.z>=boxmax.z? boxmax.z: FLT_MAX));
    if(zolimit[ci].x==FLT_MAX && zolimit[ci].z==FLT_MAX && (Simulate2D || zolimit[ci].y==FLT_MAX))
      Run_Exceptioon(fun::PrintStr("FreeCentre position (%g,%g,%g) is within the inout zone %d. FreeCentre must be changed in XML file.",FreeCentre.x,FreeCentre.y,FreeCentre.z,ci));
  }
  //-Look for the best solution for FreeLimitMin/Max.
  const byte nsel=(Simulate2D? 2: 3);
  float bestsize=-FLT_MAX;
  tfloat3 bestmin,bestmax;
  sel[0]=0;
  dmin[0]=pmin;
  dmax[0]=pmax;
  unsigned ci=0;
  bool run=true;
  //unsigned cd=0;
  while(run){
    if(sel[ci]<nsel){
      sel[ci]++;
      pmin=dmin[ci]; pmax=dmax[ci];
      bool ok=false;
      if(sel[ci]==1 && zolimit[ci].x!=FLT_MAX){//-Adjust in X.
        if(zolimit[ci].x<=FreeCentre.x)pmin.x=max(pmin.x,zolimit[ci].x);
        else                           pmax.x=min(pmax.x,zolimit[ci].x);
        ok=true;
      }
      if(sel[ci]==2 && zolimit[ci].z!=FLT_MAX){//-Adjust in Z.
        if(zolimit[ci].z<=FreeCentre.z)pmin.z=max(pmin.z,zolimit[ci].z);
        else                           pmax.z=min(pmax.z,zolimit[ci].z);
        ok=true;
      }
      if(sel[ci]==3 && zolimit[ci].y!=FLT_MAX){//-Adjust in Y.
        if(zolimit[ci].y<=FreeCentre.y)pmin.y=max(pmin.y,zolimit[ci].y);
        else                           pmax.y=min(pmax.y,zolimit[ci].y);
        ok=true;
      }
      if(ok){
        const tfloat3 ss=pmax-pmin;
        const float size=(Simulate2D? ss.x*ss.z: ss.x*ss.y*ss.z);
        if(size>bestsize){
          if(ci+1==nzone){//-Last zone was used.
            //if(1){
            //  string tx;
            //  for(unsigned cc=0;cc<nzone;cc++)tx=tx+fun::PrintStr("-[%u]",sel[cc]);
            //  Log->Printf("-----> %u >---%s: %f",cd,tx.c_str(),size);cd++;
            //  JVtkLib sh;
            //  sh.AddShapeBoxSize(pmin,pmax-pmin,cd);
            //  sh.SaveShapeVtk(Log->GetDirOut()+fun::FileNameSec("CfgInOut_DG_FreeDomain.vtk",cd),"");
            //}
            bestsize=size;
            bestmin=pmin; bestmax=pmax;
          }
          else{//-Use next zone.
            ci++;
            sel[ci]=0;
            dmin[ci]=pmin; dmax[ci]=pmax;
          }
        }
      }
    }
    else{
      sel[ci]=0;
      if(ci>0)ci--;    //-Go to previous zone
      else run=false;  //-or finish.
    }
  }
  //-Saves best solution.
  FreeLimitMin=bestmin;
  FreeLimitMax=bestmax;
  //SaveVtkDomains();
  //Run_Exceptioon("Stop");
}

//==============================================================================
/// Saves domain zones in VTK file.
//==============================================================================
void JSphInOut::SaveVtkDomains(){
  //-InOut real domains.
  {
    JVtkLib sh;
    for(unsigned ci=0;ci<GetCount();ci++){
      const JSphInOutZone *izone=List[ci];
      const tdouble3* ptdom=izone->GetPtDomain();
      if(Simulate2D)sh.AddShapeQuad(ptdom[0],ptdom[1],ptdom[2],ptdom[3],ci);
      else sh.AddShapeBoxFront(ptdom[0],ptdom[1],ptdom[2],ptdom[3],ptdom[4],ptdom[5],ptdom[6],ptdom[7],ci);
      sh.AddShapeLine(ptdom[8],ptdom[9],ci); //-Normal line.
    }
    const string filevtk=AppInfo.GetDirOut()+"CfgInOut_DomainReal.vtk";
    sh.SaveShapeVtk(filevtk,"izone");
    Log->AddFileInfo(filevtk,"Saves real domain of InOut configurations.");
  }
  //-InOut box domains.
  {
    JVtkLib sh;
    for(unsigned ci=0;ci<GetCount();ci++){
      tfloat3 boxmin=List[ci]->GetBoxLimitMin();
      tfloat3 boxmax=List[ci]->GetBoxLimitMax();
      if(Simulate2D){
        boxmin.y=boxmax.y=float(Simulate2DPosY);
        const tfloat3 pt1=TFloat3(boxmax.x,boxmin.y,boxmin.z);
        const tfloat3 pt2=TFloat3(boxmin.x,boxmax.y,boxmax.z);
        sh.AddShapeQuad(boxmin,pt1,boxmax,pt2,ci);
      }
      else sh.AddShapeBoxSize(boxmin,boxmax-boxmin,ci);
    }
    //-Draws FreeCentre.
    {
      tfloat3 pc0=FreeCentre-TFloat3(float(Dp/2));
      tfloat3 pc2=FreeCentre+TFloat3(float(Dp/2));
      tfloat3 pc1=TFloat3(pc2.x,pc0.y,pc0.z);
      tfloat3 pc3=TFloat3(pc0.x,pc2.y,pc2.z);
      if(Simulate2D){
        pc0.y=pc1.y=pc2.y=pc3.y=float(Simulate2DPosY);
        sh.AddShapeQuad(pc0,pc1,pc2,pc3,GetCount());
      }
      else sh.AddShapeBoxSize(pc0,pc2-pc0,GetCount());
    }
    //-Draws FreeLimitMin/Max.
    {
      tfloat3 pc0=MaxValues(FreeLimitMin,ToTFloat3(MapRealPosMin));
      tfloat3 pc2=MinValues(FreeLimitMax,ToTFloat3(MapRealPosMax));
      tfloat3 pc1=TFloat3(pc2.x,pc0.y,pc0.z);
      tfloat3 pc3=TFloat3(pc0.x,pc2.y,pc2.z);
      if(Simulate2D){
        pc0.y=pc1.y=pc2.y=pc3.y=float(Simulate2DPosY);
        sh.AddShapeQuad(pc0,pc1,pc2,pc3,GetCount());
      }
      else sh.AddShapeBoxSize(pc0,pc2-pc0,GetCount());
    }
    const string filevtk=AppInfo.GetDirOut()+"CfgInOut_DomainBox.vtk";
    sh.SaveShapeVtk(filevtk,"izone");
    Log->AddFileInfo(filevtk,"Saves box domain of InOut configurations.");
  }
}

//==============================================================================
/// Saves VTK of InputVelGrid nodes.
//==============================================================================
void JSphInOut::SaveVtkVelGrid(){
  for(unsigned ci=0;ci<GetCount();ci++)if(List[ci]->GetInterpolatedVel()){
    const string filevtk=AppInfo.GetDirOut()+fun::PrintStr("CfgInOut_VelGrid_i%d.vtk",ci);
    List[ci]->SaveVtkVelGrid(filevtk);
    Log->AddFileInfo(filevtk,"Saves VTK file with InputVelGrid nodes (by JSphInOut).");
  }
}

//==============================================================================
/// Configures basic parameter of the simulation and prepares execution and 
/// returns number of initial inlet particles.
//==============================================================================
unsigned JSphInOut::Config(double timestep,bool stable,bool simulate2d,double simulate2dposy
  ,byte periactive,float rhopzero,float cteb,float gamma,tfloat3 gravity,double dp
  ,tdouble3 posmin,tdouble3 posmax,typecode codenewpart,const JDsPartsInit *partsdata
  ,JGaugeSystem *gaugesystem,JNumexLib *nuxlib)
{
  Stable=stable;
  Simulate2D=simulate2d;
  Simulate2DPosY=simulate2dposy;
  PeriActive=periactive;
  RhopZero=rhopzero;
  CteB=cteb;
  Gamma=gamma;
  GravityZ=gravity.z;
  if(gravity.x!=0 || gravity.y!=0)Log->PrintfWarning("Only gravity in Z (0,0,%g) is used in inlet/outlet code (e.g.: hydrostatic density or water elevation calculation).",GravityZ);
  CoefHydro=RhopZero*(-GravityZ)/CteB;
  Dp=dp;
  MapRealPosMin=posmin; MapRealPosMax=posmax;
  CodeNewPart=codenewpart;
  //-Loads Xml configuration.
  LoadFileXml(XmlFile,XmlPath,nuxlib,partsdata,gaugesystem);

  //-Calculates and saves domain zones.
  ComputeFreeDomain();
  SaveVtkDomains();
  //-Saves VTK of InputVelGrid nodes.
  SaveVtkVelGrid();

  //-Allocates memory for inlet/outlet configurations.
  AllocateMemory(GetCount());
  //-Prepares data for inlet/outlet configurations.
  for(unsigned ci=0;ci<ListSize;ci++){
    Planes   [ci]=List[ci]->GetPlane();
    CfgZone  [ci]=List[ci]->GetConfigZone();
    CfgUpdate[ci]=List[ci]->GetConfigUpdate();
    Width    [ci]=float(Dp*List[ci]->GetLayers());
    DirData  [ci]=ToTFloat3(List[ci]->GetDirection());
    Zbottom  [ci]=List[ci]->GetInputZbottom();
  }
  UpdateVelData(timestep,true);
  UpdateZsurf(timestep,true);

  #ifdef _WITHGPU
    if(INOUT_RefillAdvanced_MASK!=JSphInOutZone::RefillAdvanced_MASK)Run_Exceptioon("RefillAdvanced mask does not match.");
    if(INOUT_RefillSpFull_MASK  !=JSphInOutZone::RefillSpFull_MASK  )Run_Exceptioon("RefillSpFull mask does not match.");
    if(INOUT_RemoveInput_MASK   !=JSphInOutZone::RemoveInput_MASK   )Run_Exceptioon("RemoveInput mask does not match.");
    if(INOUT_RemoveZsurf_MASK   !=JSphInOutZone::RemoveZsurf_MASK   )Run_Exceptioon("RemoveZsurf mask does not match.");
    if(INOUT_ConvertInput_MASK  !=JSphInOutZone::ConvertInput_MASK  )Run_Exceptioon("ConvertInput mask does not match.");
    if(!Cpu){
      //-Copies data to GPU memory.
      cudaMemcpy(Planesg   ,Planes   ,sizeof(float4)*ListSize,cudaMemcpyHostToDevice);
      cudaMemcpy(CfgZoneg  ,CfgZone  ,sizeof(byte)  *ListSize,cudaMemcpyHostToDevice);
      cudaMemcpy(CfgUpdateg,CfgUpdate,sizeof(byte)  *ListSize,cudaMemcpyHostToDevice);
      cudaMemcpy(Widthg    ,Width    ,sizeof(float) *ListSize,cudaMemcpyHostToDevice);
      cudaMemcpy(DirDatag  ,DirData  ,sizeof(float3)*ListSize,cudaMemcpyHostToDevice);
      //cudaMemcpy(Zsurfg  ,Zsurf  ,sizeof(float) *ListSize,cudaMemcpyHostToDevice); //It is done in UpdateZsurf().
      //-Copies data for BoxLimitg to GPU memory.
      if(UseBoxLimit){
        tfloat2* boxlimit=new tfloat2[ListSize*3];
        for(unsigned ci=0;ci<ListSize;ci++){
          const tfloat3 boxmin=List[ci]->GetBoxLimitMin();
          const tfloat3 boxmax=List[ci]->GetBoxLimitMax();
          boxlimit[ci]=TFloat2(boxmin.x,boxmax.x);
          boxlimit[ListSize+ci]=TFloat2(boxmin.y,boxmax.y);
          boxlimit[ListSize*2+ci]=TFloat2(boxmin.z,boxmax.z);
        }
        cudaMemcpy(BoxLimitg,boxlimit,sizeof(float2)*ListSize*3,cudaMemcpyHostToDevice);
        delete[] boxlimit; boxlimit=NULL;
      }
    }
  #endif

  //-Calculates number of particles for resize memory (NpResizePlus0 and NpResizePlus1).
  {
    unsigned npfull=0;
    for(unsigned ci=0;ci<ListSize;ci++)npfull+=List[ci]->GetNptInit()*List[ci]->GetLayers();
    NpResizePlus0=unsigned(MemoryResize0*npfull);
    NpResizePlus1=unsigned(MemoryResize1*npfull);
    //-NpResizePlus1 is only zero when MemoryResize1 is also zero, since MemoryResize1=0 does not allow resizing.
    if(!NpResizePlus1 && MemoryResize1>0)NpResizePlus1=1;
  }

  //-Calculates total number of inout points and total number of inout particles.
  unsigned npt=0,npart=0;
  for(unsigned ci=0;ci<ListSize;ci++){
    npt+=List[ci]->GetNptInit();
    npart+=List[ci]->GetNpartInit();
  }

  //-Checks if some inlet configuration uses extrapolated data or variable velocity.
  RefillAdvanced=false;
  VariableVel=false;
  VariableZsurf=false;
  CalculatedZsurf=false;
  ExtrapolatedData=false;
  NoExtrapolatedData=false;
  InterpolatedVel=false;
  for(unsigned ci=0;ci<GetCount();ci++){
    const JSphInOutZone* zo=List[ci];
    if(zo->GetRefillAdvanced())RefillAdvanced=true; 
    if(zo->GetInterpolatedVel())InterpolatedVel=true;
    if(zo->GetExtrapolatedData())ExtrapolatedData=true;
    if(zo->GetNoExtrapolatedData())NoExtrapolatedData=true;
    if(zo->GetVariableVel())VariableVel=true;
    if(zo->GetVariableZsurf())VariableZsurf=true;
    if(zo->GetCalculatedZsurf())CalculatedZsurf=true;
  }

  //-Prepares data of points to refilling or calculated zsurf.
  if(RefillAdvanced || CalculatedZsurf){
    AllocatePtMemory(npt);
    npt=0;
    for(unsigned ci=0;ci<ListSize;ci++){
      const unsigned n=List[ci]->LoadInletPoints(PtPos+npt);
      memset(PtZone+npt,byte(List[ci]->GetIdZone()),sizeof(byte)*n);
      npt+=n;
    }
    #ifdef _WITHGPU
      if(!Cpu){
        //-Allocates and prepares auxiliary memory.
        tdouble2 *pxy=new tdouble2[PtCount];
        double   *pz =new double  [PtCount];
        for(unsigned c=0;c<PtCount;c++){
          pxy[c]=TDouble2(PtPos[c].x,PtPos[c].y);
          pz[c]=PtPos[c].z;
        }
        //-Copies data to GPU memory.
        cudaMemcpy(PtZoneg ,PtZone,sizeof(byte)   *PtCount,cudaMemcpyHostToDevice);
        cudaMemcpy(PtPosxyg,pxy   ,sizeof(double2)*PtCount,cudaMemcpyHostToDevice);
        cudaMemcpy(PtPoszg ,pz    ,sizeof(double) *PtCount,cudaMemcpyHostToDevice);
        //-Frees auxiliary memory.
        delete[] pxy; pxy=NULL;
        delete[] pz;  pz=NULL;
      }
    #endif

    //-Creates VTK file.
    if(DBG_INOUT_PTINIT){
      JDataArrays arrays;
      arrays.AddArray("Pos",PtCount,PtPos,false);
      if(PtZone)arrays.AddArray("PtZone",PtCount,PtZone,false);
      const string filevtk=AppInfo.GetDirOut()+"CfgInOut_PtInit.vtk";
      JVtkLib::SaveVtkData(filevtk,arrays,"Pos");
      Log->AddFileInfo(filevtk,"Saves initial InOut points for DEBUG (by JSphInOut).");
    }
    //-Creates VTK file.
    if(DBG_INOUT_PTINIT){
      unsigned npz=0;
      for(unsigned ci=0;ci<ListSize;ci++)npz+=List[ci]->GetCountPtzPos();
      if(npz){
        tfloat3 *pos =new tfloat3[npz];
        byte    *zone=new byte[npz];
        npz=0;
        for(unsigned ci=0;ci<ListSize;ci++){
          unsigned n=List[ci]->GetCountPtzPos();
          memcpy(pos+npz,List[ci]->GetPtzPos(),sizeof(tfloat3)*n);
          memset(zone+npz,byte(ci),sizeof(byte)*n);
          npz+=n;
        }
        JDataArrays arrays;
        arrays.AddArray("Pos",npz,pos,false);
        if(zone)arrays.AddArray("Zone",npz,zone,false);
        const string filevtk=AppInfo.GetDirOut()+"CfgInOut_PtInitZ.vtk";
        JVtkLib::SaveVtkData(filevtk,arrays,"Pos");
        Log->AddFileInfo(filevtk,"Saves initial InOut points for DEBUG (by JSphInOut).");
        //-Frees memory.
        delete[] pos;  pos=NULL;
        delete[] zone; zone=NULL;
      }
    }
  }
  return(npart);
}

//==============================================================================
/// Loads basic data (pos,idp,code,velrhop=0) for initial inout particles.
//==============================================================================
void JSphInOut::LoadInitPartsData(unsigned idpfirst,unsigned nparttot
  ,unsigned* idp,typecode* code,tdouble3* pos,tfloat4* velrhop)
{
  //Log->Printf(" LoadInitPartsData--> nparttot:%u",nparttot);
  unsigned npart=0;
  for(unsigned ci=0;ci<GetCount();ci++){
    const unsigned np=List[ci]->GetNpartInit();
    //Log->Printf(" LoadInitPartsData--> np:%u  npart:%u",np,npart);
    if(npart+np>nparttot)Run_Exceptioon("Number of initial inlet/outlet particles is invalid.");
    List[ci]->LoadInletParticles(pos+npart);
    for(unsigned cp=0;cp<np;cp++){
      const unsigned p=npart+cp;
      idp[p]=idpfirst+p;
      code[p]=typecode(CODE_TYPE_FLUID_INOUT)+ci;
      velrhop[p]=TFloat4(0,0,0,1000);
    }
    npart+=np;
  }
  //Log->Printf(" LoadInitPartsData--> npart:%u",npart);
  if(npart!=nparttot)Run_Exceptioon("Number of initial inlet/outlet particles is invalid.");
}

//==============================================================================
/// Checks proximity of inout particles to other particles and excludes fluid 
/// particles near the inout particles.
///
/// Comprueba proximidad de particulas inout con otras particulas y excluye 
/// particulas fluid cerca de particulas inout.
//==============================================================================
void JSphInOut::InitCheckProximity(unsigned np,unsigned newnp,float scell
  ,const tdouble3* pos,const unsigned *idp,typecode *code)
{
  //-Look for nearby particles.
  const double disterror=Dp*0.8;
  JSimpleNeigs neigs(np,pos,scell);
  byte* errpart=new byte[np];
  memset(errpart,0,sizeof(byte)*np);
  const unsigned pini=np-newnp;
  JTimeControl tc(5,60);
  for(unsigned p=pini;p<np;p++){//-Only inout particles.
    const unsigned n=neigs.NearbyPositions(pos[p],p,disterror);
    const unsigned *selpos=neigs.GetSelectPos();
    for(unsigned cp=0;cp<n;cp++)errpart[selpos[cp]]=1;
    if(tc.CheckTime())Log->Print(string("  ")+tc.GetInfoFinish(double(p-pini)/double(np-pini)));
  }
  //-Obtain number and type of nearby particles.
  unsigned nfluid=0,nfluidinout=0,nbound=0;
  for(unsigned p=0;p<np;p++)if(errpart[p]){
    const typecode cod=code[p];
    if(CODE_IsNormal(cod)){
      if(CODE_IsFluid(cod)){
        if(CODE_IsFluidNotInout(cod)){ //-Normal fluid.
          errpart[p]=1;
          nfluid++;
        }
        else{ //-Inout fluid.
          errpart[p]=2;
          nfluidinout++;
        }
      }
      else{ //-Boundary.
        errpart[p]=3;
        nbound++;
      } 
    }
    else errpart[p]=0; //-Ignores non-normal particles.
  }
  const unsigned nerr=nfluid+nfluidinout+nbound;
  //-Saves VTK file with nearby particles.
  if(nerr>0){
    const unsigned n=nfluid+nfluidinout+nbound;
    tfloat3* vpos=new tfloat3[n];
    byte* vtype=new byte[n];
    unsigned pp=0;
    for(unsigned p=0;p<np;p++)if(errpart[p]){
      vpos[pp]=ToTFloat3(pos[p]);
      vtype[pp]=errpart[p];
      pp++;
    }
    JDataArrays arrays;
    arrays.AddArray("Pos",n,vpos,false);
    arrays.AddArray("ErrorType",n,vtype,false);
    const string filevtk=AppInfo.GetDirOut()+(n>nfluid? "CfgInOut_ErrorParticles.vtk": "CfgInOut_ExcludedParticles.vtk");
    JVtkLib::SaveVtkData(filevtk,arrays,"Pos");
    if(nerr>nfluid)Log->AddFileInfo(filevtk,"Saves error fluid and boundary particles too close to inout particles.");
    else Log->AddFileInfo(filevtk,"Saves excluded fluid particles too close to inout particles.");
    delete[] vpos;  vpos=NULL;
    delete[] vtype; vtype=NULL;
  }
  //-Checks errors and remove nearby fluid particles.
  if(nerr>nfluid)Run_Exceptioon("There are inout fluid or boundary particles too close to inout particles. Check VTK file CfgInOut_ErrorParticles.vtk with excluded particles.");
  else{
    if(nfluid)Log->PrintfWarning("%u fluid particles were excluded since they are too close to inout particles. Check VTK file CfgInOut_ExcludedParticles.vtk",nfluid);
    //-Mark fluid particles to ignore.
    for(unsigned p=0;p<np;p++)if(errpart[p]==1){
      code[p]=CODE_SetOutIgnore(code[p]); //-Mark fluid particles to ignore.
    }
  }
  //-Free memory.
  delete[] errpart; errpart=NULL;
}

//==============================================================================
/// Creates list with current inout particles (normal and periodic).
//==============================================================================
unsigned JSphInOut::CreateListSimpleCpu(unsigned nstep,unsigned npf,unsigned pini
  ,const typecode *code,int *inoutpart)
{
  unsigned count=0;
  if(ListSize){
    const unsigned pfin=pini+npf;
    for(unsigned p=pini;p<pfin;p++){
      const typecode rcode=code[p];
      if(CODE_IsNotOut(rcode) && CODE_IsFluidInout(rcode)){//-It includes normal and periodic particles.
        inoutpart[count]=p; count++;
      }
    }
    if(0){ //DG_INOUT
      Log->Printf("AAA_000 ListSize:%u",ListSize);
      Log->Printf("AAA_000 Planes[0]:(%f,%f,%f,%f)",Planes[0].a,Planes[0].b,Planes[0].c,Planes[0].d);
      const float xini=-5.f;
      const float dp=0.01f;
      const unsigned np=1000;
      tfloat3 vpos[np];
      float vdis0[np];
      float vdis1[np];
      for(unsigned p=0;p<np;p++){
        vpos[p]=TFloat3(xini+dp*p,0,0);
        const float dis0=fgeo::PlanePoint(Planes[0],vpos[p]);
        vdis0[p]=(dis0<0? -1.f: (dis0>0? 1.f: 0));
        const float dis1=fgeo::PlanePoint(Planes[1],vpos[p]);
        vdis1[p]=(dis1<0? -1.f: (dis1>0? 1.f: 0));
      }
      //-Generates VTK file.
      JDataArrays arrays;
      arrays.AddArray("Pos",np,vpos,false);
      if(vdis0)arrays.AddArray("vdis0",np,vdis0,false);
      if(vdis1)arrays.AddArray("vdis1",np,vdis1,false);
      JVtkLib::SaveVtkData(AppInfo.GetDirOut()+"_Planes.vtk",arrays,"Pos");
      //-Old Style...
      //JFormatFiles2::StScalarData fields[5];
      //unsigned nfields=0;
      //if(vdis0){ fields[nfields]=JFormatFiles2::DefineField("vdis0",JFormatFiles2::Float32,1,vdis0);    nfields++; }
      //if(vdis1){ fields[nfields]=JFormatFiles2::DefineField("vdis1",JFormatFiles2::Float32,1,vdis1);    nfields++; }
      ////string fname=DirOut+fun::FileNameSec("DgParts.vtk",numfile);
      //JFormatFiles2::SaveVtk(AppInfo.GetDirOut()+"_Planes.vtk",np,vpos,nfields,fields);
    }
  }
  //Log->Printf("%u> -------->CreateListXXX>> InOutcount:%u",nstep,count);
  return(count);
}

//==============================================================================
/// Creates list with current inout particles and normal (no periodic) fluid in 
/// inlet/outlet zones (update its code).
//==============================================================================
unsigned JSphInOut::CreateListCpu(unsigned nstep,unsigned npf,unsigned pini
  ,const tdouble3 *pos,const unsigned *idp,typecode *code,int *inoutpart)
{
  unsigned count=0;
  if(ListSize){
    const byte chkinputmask=byte(JSphInOutZone::CheckInput_MASK);
    const bool checkfreelimit=(UseBoxLimit && ListSize>2);
    const unsigned pfin=pini+npf;
    for(unsigned p=pini;p<pfin;p++){
      const typecode rcode=code[p];
      if(CODE_IsNormal(rcode) && CODE_IsFluid(rcode)){//-It includes only normal fluid particles (no periodic).
        if(CODE_IsFluidInout(rcode)){//-Particles already selected as InOut.
          inoutpart[count]=p; count++;
        }
        else{//-Fluid particles no inout.
          const tfloat3 ps=ToTFloat3(pos[p]);
          if(!checkfreelimit || ps.x<=FreeLimitMin.x || FreeLimitMax.x<=ps.x || ps.z<=FreeLimitMin.z || FreeLimitMax.z<=ps.z || ps.y<=FreeLimitMin.y || FreeLimitMax.y<=ps.y){
            byte zone=255;
            for(unsigned cp=0;cp<ListSize && zone==255;cp++)
              if((CfgZone[cp]&chkinputmask)!=0 && List[cp]->InZone(UseBoxLimit,ps))zone=byte(cp);
            if(zone!=255){//-Particulas fluid que pasan a in/out.
              code[p]=CODE_ToFluidInout(rcode,zone)|CODE_TYPE_FLUID_INOUTNUM; //-Adds 16 to indicate new particle in zone.
              inoutpart[count]=p; count++;
            }
          }
        }
      }
    }
//    Log->Printf("==>> nold:%d  nnew:%d",nold,nnew);
    if(0){ //DG_INOUT
      Log->Printf("AAA_000 ListSize:%u",ListSize);
      Log->Printf("AAA_000 Planes[0]:(%f,%f,%f,%f)",Planes[0].a,Planes[0].b,Planes[0].c,Planes[0].d);
      const float xini=-5.f;
      const float dp=0.01f;
      const unsigned np=1000;
      tfloat3 vpos[np];
      float vdis0[np];
      float vdis1[np];
      for(unsigned p=0;p<np;p++){
        vpos[p]=TFloat3(xini+dp*p,0,0);
        const float dis0=fgeo::PlanePoint(Planes[0],vpos[p]);
        vdis0[p]=(dis0<0? -1.f: (dis0>0? 1.f: 0));
        const float dis1=fgeo::PlanePoint(Planes[1],vpos[p]);
        vdis1[p]=(dis1<0? -1.f: (dis1>0? 1.f: 0));
      }
      //-Generates VTK file.
      JDataArrays arrays;
      arrays.AddArray("Pos",np,vpos,false);
      if(vdis0)arrays.AddArray("vdis0",np,vdis0,false);
      if(vdis1)arrays.AddArray("vdis1",np,vdis1,false);
      JVtkLib::SaveVtkData(AppInfo.GetDirOut()+"_Planes.vtk",arrays,"Pos");
      //-Old Style...
      //JFormatFiles2::StScalarData fields[5];
      //unsigned nfields=0;
      //if(vdis0){ fields[nfields]=JFormatFiles2::DefineField("vdis0",JFormatFiles2::Float32,1,vdis0);    nfields++; }
      //if(vdis1){ fields[nfields]=JFormatFiles2::DefineField("vdis1",JFormatFiles2::Float32,1,vdis1);    nfields++; }
      ////string fname=DirOut+fun::FileNameSec("DgParts.vtk",numfile);
      //JFormatFiles2::SaveVtk(AppInfo.GetDirOut()+"_Planes.vtk",np,vpos,nfields,fields);
    }
  }
  //Log->Printf("%u> -------->CreateListXXX>> InOutcount:%u",nstep,count);
  return(count);
}

#ifdef _WITHGPU
//==============================================================================
/// Creates list with current inout particles (normal and periodic).
//==============================================================================
unsigned JSphInOut::CreateListSimpleGpu(unsigned nstep,unsigned npf,unsigned pini
  ,const typecode *codeg,unsigned size,int *inoutpartg)
{
  unsigned count=0;
  if(ListSize){
    if(npf+2>=size)Run_Exceptioon("GPU memory allocated is not enough.");
    count=cusphinout::InOutCreateListSimple(Stable,npf,pini,codeg,(unsigned*)inoutpartg);
  }
  //Log->Printf("%u> -------->CreateListXXX>> InOutcount:%u",nstep,count);
  return(count);
}
//==============================================================================
/// Creates list with current inout particles and normal (no periodic) fluid in 
/// inlet/outlet zones (update its code).
//==============================================================================
unsigned JSphInOut::CreateListGpu(unsigned nstep,unsigned npf,unsigned pini
  ,const double2 *posxyg,const double *poszg,typecode *codeg,unsigned size,int *inoutpartg)
{
  unsigned count=0;
  if(ListSize){
    if(npf+2>=size)Run_Exceptioon("GPU memory allocated is not enough.");
    const tfloat3 freemin=(UseBoxLimit? FreeLimitMin: TFloat3(FLT_MAX));
    const tfloat3 freemax=(UseBoxLimit? FreeLimitMax: TFloat3(-FLT_MAX));
    const byte chkinputmask=byte(JSphInOutZone::CheckInput_MASK);
    count=cusphinout::InOutCreateList(Stable,npf,pini,chkinputmask,byte(ListSize)
      ,CfgZoneg,Planesg,freemin,freemax,(UseBoxLimit? BoxLimitg: NULL),posxyg,poszg,codeg,(unsigned*)inoutpartg);
  }
  //Log->Printf("%u> -------->CreateListXXX>> InOutcount:%u",nstep,count);
  return(count);
}
#endif

//==============================================================================
/// Updates velocity coefficients for imposed velocity according timestep. 
/// Returns true when the data was changed.
//==============================================================================
bool JSphInOut::UpdateVelData(double timestep,bool full){
  bool modified=full;
  for(unsigned ci=0;ci<ListSize;ci++)if(full || List[ci]->GetVariableVel()){
    tfloat4 vel0,vel1;
    if(List[ci]->UpdateVelData(timestep,full,vel0,vel1)){
      modified=true;
      VelData[ci*2]=vel0;
      VelData[ci*2+1]=vel1;
    }
  }
  return(modified);
}

//==============================================================================
/// Updates zsurf according timestep. 
/// Returns true when the data was changed.
//==============================================================================
bool JSphInOut::UpdateZsurf(double timestep,bool full){
  bool modified=full;
  for(unsigned ci=0;ci<ListSize;ci++)if(full || List[ci]->GetVariableZsurf() || List[ci]->GetCalculatedZsurf()){
    modified|=List[ci]->UpdateZsurf(timestep,full,Zsurf[ci]);
  }
  #ifdef _WITHGPU
    if(modified && !Cpu)cudaMemcpy(Zsurfg,Zsurf,sizeof(float)*ListSize,cudaMemcpyHostToDevice);
  #endif
  return(modified);
}

//==============================================================================
/// Updates velocity and rhop of inlet/outlet particles when it is not extrapolated. 
/// Actualiza velocidad y densidad de particulas inlet/outlet cuando no es extrapolada.
//==============================================================================
void JSphInOut::UpdateDataCpu(float timestep,bool full,unsigned inoutcount,const int *inoutpart
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop)
{
  const bool modifvel=UpdateVelData(timestep,full);
  //const bool modifzsurf=UpdateZsurf(timestep,full);
  const int ncp=int(inoutcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int cp=0;cp<ncp;cp++){
    const unsigned p=(unsigned)inoutpart[cp];
    const unsigned izone=CODE_GetIzoneFluidInout(code[p]);
    const byte cfg=CfgZone[izone];
    const bool refillspfull=(CfgUpdate[izone]&JSphInOutZone::RefillSpFull_MASK)!=0;
    const double posz=pos[p].z;
    tfloat4 rvelrhop=velrhop[p];
    //-Compute rhop value.
    const JSphInOutZone::TpRhopMode rmode=JSphInOutZone::GetConfigRhopMode(cfg);
    if(rmode==JSphInOutZone::MRHOP_Constant)rvelrhop.w=RhopZero;
    if(rmode==JSphInOutZone::MRHOP_Hydrostatic){
      const float depth=float(double(Zsurf[izone])-posz);
      const float rh=1.f+CoefHydro*depth;     //rh=1.+rhop0*(-gravity.z)*(Dp*ptdata.GetDepth(p))/vCteB;
      rvelrhop.w=RhopZero*pow(rh,1.f/Gamma);  //rhop[id]=rhop0*pow(rh,(1./gamma));
    }
    //-Compute velocity value.
    const JSphInOutZone::TpVelMode    vmode=JSphInOutZone::GetConfigVelMode(cfg);
    const JSphInOutZone::TpVelProfile vprof=JSphInOutZone::GetConfigVelProfile(cfg);
    float vel=0;
    if(!refillspfull || posz<=Zsurf[izone]){//-It is necessary for RefillingMode==Simple-Full
      if(vmode==JSphInOutZone::MVEL_Fixed){
        vel=JSphInOutZone::CalcVel(vprof,VelData[izone*2],posz);
      }
      else if(vmode==JSphInOutZone::MVEL_Variable){
        const float vel1=JSphInOutZone::CalcVel(vprof,VelData[izone*2],posz);
        const float vel2=JSphInOutZone::CalcVel(vprof,VelData[izone*2+1],posz);
        const float time1=VelData[izone*2].w;
        const float time2=VelData[izone*2+1].w;
        if(timestep<=time1 || time1==time2)vel=vel1;
        else if(timestep>=time2)vel=vel2;
        else vel=(timestep-time1)/(time2-time1)*(vel2-vel1)+vel1;
      }
    }
    if(vmode!=JSphInOutZone::MVEL_Extrapolated){
      rvelrhop.x=vel*DirData[izone].x;
      rvelrhop.y=vel*DirData[izone].y;
      rvelrhop.z=vel*DirData[izone].z;
    }
    velrhop[p]=rvelrhop;
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Updates velocity and rhop of inlet/outlet particles when it is not extrapolated. 
/// Actualiza velocidad y densidad de particulas inlet/outlet cuando no es extrapolada.
//==============================================================================
void JSphInOut::UpdateDataGpu(float timestep,bool full,unsigned inoutcount,const int *inoutpartg
  ,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg,float4 *velrhopg)
{
  const bool modifvel=UpdateVelData(timestep,full);
  //const bool modifzsurf=UpdateZsurf(timestep,full);
  for(unsigned izone=0;izone<ListSize;izone++){
    const byte refillspfull=((CfgUpdate[izone]&JSphInOutZone::RefillSpFull_MASK)!=0? 1: 0);
    const byte cfg=CfgZone[izone];
    const JSphInOutZone::TpRhopMode   rmode=JSphInOutZone::GetConfigRhopMode(cfg);
    const JSphInOutZone::TpVelMode    vmode=JSphInOutZone::GetConfigVelMode(cfg);
    const JSphInOutZone::TpVelProfile vprof=JSphInOutZone::GetConfigVelProfile(cfg);
    const byte brmode=(rmode==JSphInOutZone::MRHOP_Constant? 0: (rmode==JSphInOutZone::MRHOP_Hydrostatic? 1: (rmode==JSphInOutZone::MRHOP_Extrapolated? 2: 99)));
    const byte bvmode=(vmode==JSphInOutZone::MVEL_Fixed?     0: (vmode==JSphInOutZone::MVEL_Variable?     1: (vmode==JSphInOutZone::MVEL_Extrapolated?  2: 99)));
    const byte bvprof=(vprof==JSphInOutZone::PVEL_Constant?  0: (vprof==JSphInOutZone::PVEL_Linear?       1: (vprof==JSphInOutZone::PVEL_Parabolic?     2: 99)));
    cusphinout::InOutUpdateData(inoutcount,(unsigned*)inoutpartg
      ,byte(izone),brmode,bvmode,bvprof,refillspfull
      ,timestep,Zsurf[izone],VelData[izone*2],VelData[izone*2+1],DirData[izone]
      ,CoefHydro,RhopZero,Gamma,codeg,poszg,velrhopg);
  }
}
#endif

//==============================================================================
/// Interpolate velocity of inlet/outlet particles from data in InputVelGrid object.
/// Interpola velocidad de particulas inlet/outlet a partir de datos en el objeto InputVelGrid.
//==============================================================================
void JSphInOut::InterpolateVelCpu(float timestep,unsigned inoutcount,const int *inoutpart
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop)
{
  for(unsigned ci=0;ci<GetCount();ci++)if(List[ci]->GetInterpolatedVel()){
    JSphInOutZone* zo=List[ci];
    float velcorr=0;
    JSphInOutGridData* gd=zo->GetInputVelGrid();
    if(gd->GetNx()==1)gd->InterpolateZVelCpu(timestep,byte(ci),inoutcount,inoutpart,pos,code,idp,velrhop,velcorr);
    else gd->InterpolateVelCpu(timestep,byte(ci),inoutcount,inoutpart,pos,code,idp,velrhop,velcorr);
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Interpolate velocity of inlet/outlet particles from data in InputVelGrid object.
/// Interpola velocidad de particulas inlet/outlet a partir de datos en el objeto InputVelGrid.
//==============================================================================
void JSphInOut::InterpolateVelGpu(float timestep,unsigned inoutcount,const int *inoutpartg
  ,double2 *posxyg,double *poszg,typecode *codeg,unsigned *idpg,float4 *velrhopg)
{
  for(unsigned ci=0;ci<GetCount();ci++)if(List[ci]->GetInterpolatedVel()){
    JSphInOutZone* zo=List[ci];
    float velcorr=0;
    JSphInOutGridData* gd=zo->GetInputVelGrid();
    if(gd->GetNx()==1)gd->InterpolateZVelGpu(timestep,byte(ci),inoutcount,inoutpartg,posxyg,poszg,codeg,idpg,velrhopg,velcorr);
    else Run_Exceptioon("GPU code was not implemented for nx>1.");
  }
}
#endif

//==============================================================================
/// Removes interpolated Z velocity of inlet/outlet particles.
/// Elimina velocidad interpolada en Z de particulas inlet/outlet.
//==============================================================================
void JSphInOut::InterpolateResetZVelCpu(unsigned inoutcount,const int *inoutpart
  ,const typecode *code,tfloat4 *velrhop)
{
  for(unsigned ci=0;ci<GetCount();ci++)if(List[ci]->GetResetZVelGrid()){
    const JSphInOutGridData* gd=List[ci]->GetInputVelGrid();
    if(gd->GetUseVelz()){
      const int n=int(inoutcount);
      #ifdef OMP_USE
        #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
      #endif
      for(int cp=0;cp<n;cp++){
        const unsigned p=inoutpart[cp];
        if(ci==CODE_GetIzoneFluidInout(code[p]))velrhop[p].z=0;
      }
    }
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Removes interpolated Z velocity of inlet/outlet particles.
/// Elimina velocidad interpolada en Z de particulas inlet/outlet.
//==============================================================================
void JSphInOut::InterpolateResetZVelGpu(unsigned inoutcount,const int *inoutpartg
  ,typecode *codeg,float4 *velrhopg)
{
  for(unsigned ci=0;ci<GetCount();ci++)if(List[ci]->GetResetZVelGrid()){
    const JSphInOutGridData* gd=List[ci]->GetInputVelGrid();
    if(gd->GetUseVelz()){
      cusphinout::InOutInterpolateResetZVel(ci,inoutcount,inoutpartg,codeg,velrhopg);
    }
  }
}
#endif

//==============================================================================
/// Checks izone code in list of inout particles.
//==============================================================================
void JSphInOut::CheckPartsIzone(std::string key,unsigned nstep
  ,unsigned inoutcount,const int *inoutpart,typecode *code,unsigned *idp)
{
  for(unsigned c=0;c<inoutcount;c++){
    const unsigned p=inoutpart[c];
    if(CODE_IsFluidInout(code[p])){
      unsigned izone=CODE_GetIzoneFluidInout(code[p]);
      if(izone>=ListSize)Run_Exceptioon(fun::PrintStr("%d> [%s] Value izone %d is invalid of cp=%d idp[%d]=%d.",nstep,key.c_str(),izone,c,p,idp[p]));
    }
  }
}


//==============================================================================
/// ComputeStep over inout particles:
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new inout particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
unsigned JSphInOut::ComputeStepCpu(unsigned nstep,double dt,unsigned inoutcount
  ,int *inoutpart,const JSphCpu *sphcpu,unsigned idnext,unsigned sizenp,unsigned np
  ,tdouble3 *pos,unsigned *dcell,typecode *code,unsigned *idp,tfloat4 *velrhop,byte *newizone)
{
  //-Updates code according to particle position and define new particles to create.
  const int ncp=int(inoutcount);
  //Log->Printf("%d>==>> ComputeStepCpu  ncp:%d",nstep,ncp);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int cp=0;cp<ncp;cp++){
    typecode cod=0;
    byte newiz=255;
    const unsigned p=(unsigned)inoutpart[cp];
    const typecode rcode=code[p];
    const unsigned izone0=CODE_GetIzoneFluidInout(rcode);
    const unsigned izone=(izone0&0xf); //-Substract 16 to obtain the actual zone (0-15).
if(izone>=ListSize)Run_Exceptioon(fun::PrintStr("%d>> Value izone %d is invalid of idp[%d]=%d.",nstep,izone,p,idp[p]));
    const byte cfupdate=CfgUpdate[izone];
    const bool refilladvan=(cfupdate&JSphInOutZone::RefillAdvanced_MASK)!=0;
    const bool refillsfull=(cfupdate&JSphInOutZone::RefillSpFull_MASK  )!=0;
    const bool removeinput=(cfupdate&JSphInOutZone::RemoveInput_MASK   )!=0;
    const bool removezsurf=(cfupdate&JSphInOutZone::RemoveZsurf_MASK   )!=0;
    const bool converinput=(cfupdate&JSphInOutZone::ConvertInput_MASK  )!=0;
    const tfloat3 ps=ToTFloat3(pos[p]);
    if(izone0>=16){//-Normal fluid particle in zone inlet/outlet.
      if(removeinput || (removezsurf && ps.z>Zsurf[izone]))cod=CODE_SetOutPos(rcode); //-Normal fluid particle in zone inlet/outlet is removed.
      else cod=(converinput? rcode^0x10: CodeNewPart); //-Converts to inout particle or not.
    }
    else{//-Previous inout fluid particle.
      const float displane=-fgeo::PlaneDistSign(Planes[izone],ps);
      if(displane>Width[izone] || (removezsurf && ps.z>Zsurf[izone])){
        cod=CODE_SetOutIgnore(rcode); //-Particle is moved out domain.
      }
      else if(displane<0){
        cod=CodeNewPart;//-Inout particle changes to fluid particle.
        if(!refilladvan && (refillsfull || ps.z<=Zsurf[izone]))newiz=byte(izone); //-A new particle is created.
      }
    }
    newizone[cp]=newiz;
    if(cod!=0)code[p]=cod;
  }

  //-Create list for new inlet particles to create.
  unsigned inoutcount2=inoutcount;
  for(int cp=0;cp<ncp;cp++)if(newizone[cp]<16){
    if(inoutcount2<sizenp)inoutpart[inoutcount2]=cp;
    inoutcount2++;
  }
  if(inoutcount2>=sizenp)Run_Exceptioon("Allocated memory is not enough for new particles inlet.");

  //-Creates new inlet particles to replace the particles moved to fluid domain.
  const int newnp=int(inoutcount2-inoutcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int cp=0;cp<newnp;cp++){
    const unsigned cp0=(unsigned)inoutpart[inoutcount+cp];
    const unsigned p=(unsigned)inoutpart[cp0];
    const unsigned izone=newizone[cp0];
    const double dis=Width[izone];
    tdouble3 rpos=pos[p];
    rpos.x-=dis*DirData[izone].x;
    rpos.y-=dis*DirData[izone].y;
    rpos.z-=dis*DirData[izone].z;
    const unsigned p2=np+cp;
    code[p2]=CODE_ToFluidInout(CodeNewPart,izone);
    sphcpu->UpdatePos(rpos,0,0,0,false,p2,pos,dcell,code);
    idp[p2]=idnext+cp;
    velrhop[p2]=TFloat4(0,0,0,1000);
  }
  //-Returns number of new inlet particles.
  return(unsigned(newnp));
}

#ifdef _WITHGPU
//==============================================================================
/// ComputeStep over inlet/outlet particles:
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new in/out particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
unsigned JSphInOut::ComputeStepGpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpartg
  ,unsigned idnext,unsigned sizenp,unsigned np,double2 *posxyg,double *poszg
  ,unsigned *dcellg,typecode *codeg,unsigned *idpg,float4 *velrhopg,byte *newizoneg,const JSphGpuSingle *gp)
{
  //-Checks particle position.
  cusphinout::InOutComputeStep(inoutcount,inoutpartg,Planesg,Widthg,CfgUpdateg,Zsurfg
    ,CodeNewPart,posxyg,poszg,codeg,newizoneg);
  //-Create list for new inlet particles to create.
  const unsigned newnp=cusphinout::InOutListCreate(Stable,inoutcount,sizenp-1,newizoneg,inoutpartg);
  if(inoutcount+newnp>=sizenp)Run_Exceptioon("Allocated memory is not enough for new particles inlet.");
  //-Creates new inlet particles to replace the particles moved to fluid domain.
  cusphinout::InOutCreateNewInlet(PeriActive,newnp,(unsigned*)inoutpartg,inoutcount,newizoneg,np,idnext
    ,CodeNewPart,DirDatag,Widthg,posxyg,poszg,dcellg,codeg,idpg,velrhopg);
  //-Returns number of new inlet particles.
  return(unsigned(newnp));
}
#endif

//==============================================================================
/// ComputeStep over inout particles and filling inout domain:
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new inout particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
unsigned JSphInOut::ComputeStepFillingCpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpart
  ,const JSphCpu *sphcpu,unsigned idnext,unsigned sizenp,unsigned np
  ,tdouble3 *pos,unsigned *dcell,typecode *code,unsigned *idp,tfloat4 *velrhop
  ,float *prodist,tdouble3 *propos)
{
  //-Updates position of particles and computes projection data to filling mode.
  const int ncp=int(inoutcount);
  memset(prodist,0,sizeof(float)*ncp);
  memset(propos,0,sizeof(tdouble3)*ncp);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int cp=0;cp<ncp;cp++){
    const unsigned p=(unsigned)inoutpart[cp];
    const typecode rcode=code[p];
    if(CODE_IsNotOut(rcode) && CODE_IsFluidInout(rcode)){
      const unsigned izone=CODE_GetIzoneFluidInout(rcode);
      if((CfgUpdate[izone]&JSphInOutZone::RefillAdvanced_MASK)!=0){
        const tdouble3 rpos=pos[p];
        const tplane3f rplanes=Planes[izone];
        //-Compute distance to plane.
        const double v1=rpos.x*rplanes.a + rpos.y*rplanes.b + rpos.z*rplanes.c + rplanes.d;
        const double v2=rplanes.a*rplanes.a+rplanes.b*rplanes.b+rplanes.c*rplanes.c;
        prodist[cp]=-float(v1/sqrt(v2));//-Equivalent to fgeo::PlaneDistSign().
        //-Calculates point on plane.
        const double t=-v1/v2;
        propos[cp]=TDouble3(rpos.x+t*rplanes.a,rpos.y+t*rplanes.b,rpos.z+t*rplanes.c);
      }
    }
  }

  //-Compute maximum distance to create points in each PtPos.
  //const bool checkzsurf=(VariableZsurf || CalculatedZsurf);
  const float dpmin=float(Dp*1);
  const float dpmin2=dpmin*dpmin;
  const int npt=int(PtCount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static)
  #endif
  for(int cpt=0;cpt<npt;cpt++){
    float distmax=FLT_MAX;
    const byte izone=PtZone[cpt];
    if((CfgUpdate[izone]&JSphInOutZone::RefillAdvanced_MASK)!=0){
      const tdouble3 ps=PtPos[cpt];
      if(float(ps.z)<=Zsurf[izone]){
        distmax=0;
        for(int cp=0;cp<ncp;cp++){
          const tfloat3 dis=ToTFloat3(ps-propos[cp]);
          if(dis.x<=dpmin && dis.y<=dpmin && dis.z<=dpmin){//-particle near to ptpoint (approx.)
            const float dist2=(dis.x*dis.x+dis.y*dis.y+dis.z*dis.z);
            if(dist2<dpmin2){//-particle near to ptpoint.
              const float dmax=prodist[cp]+sqrt(dpmin2-dist2);
              distmax=max(distmax,dmax);
            }
          }
        }
      }
    }
    PtAuxDist[cpt]=(distmax==0? float(Dp): distmax);
  }

  //-Creates new inout particles.
  unsigned newnp=0;
//  for(int c=0;c<nc;c++)if(PtAuxDist[c]<Width[PtZone[c]]*0.99f){//-The 0.99 value avoids the creation of new particle to replace a removed particle in same step.
  for(int cpt=0;cpt<npt;cpt++)if(PtAuxDist[cpt]<Width[PtZone[cpt]]){
    const unsigned p=np+newnp;
    if(p<sizenp){
      const byte izone=PtZone[cpt];
      code[p]=CODE_ToFluidInout(CodeNewPart,izone);
      const double dis=PtAuxDist[cpt];
      tdouble3 rpos=PtPos[cpt];
      rpos.x-=dis*DirData[izone].x;
      rpos.y-=dis*DirData[izone].y;
      rpos.z-=dis*DirData[izone].z;
      sphcpu->UpdatePos(rpos,0,0,0,false,p,pos,dcell,code);
      idp[p]=idnext+newnp;
      velrhop[p]=TFloat4(0,0,0,1000);
    }
    newnp++;
  }
  if(np+newnp>=sizenp)Run_Exceptioon("Allocated memory is not enough for new particles inlet.");

  //-Returns number of new inlet particles.
  //Log->Printf("%u> -------->ComputeStepFillingXXX>> NewInletCount:%u",nstep,newnp);
  return(newnp);
}

#ifdef _WITHGPU
//==============================================================================
/// ComputeStep over inlet/outlet particles:
/// - If particle is moved to fluid zone then it changes to fluid particle and 
///   it creates a new in/out particle.
/// - If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
unsigned JSphInOut::ComputeStepFillingGpu(unsigned nstep,double dt,unsigned inoutcount,int *inoutpartg
  ,unsigned idnext,unsigned sizenp,unsigned np
  ,double2 *posxyg,double *poszg,unsigned *dcellg,typecode *codeg,unsigned *idpg,float4 *velrhopg
  ,float *prodistg,double2 *proposxyg,double *proposzg,TimersGpu timers)
{
  //-Computes projection data to filling mode.
  cusphinout::InOutFillProjection(inoutcount,(unsigned *)inoutpartg,CfgUpdateg,Planesg,posxyg,poszg
    ,codeg,prodistg,proposxyg,proposzg);

  //-Create list of selected ptpoints and its distance to create new inlet/outlet particles.
  const float dpmin=float(Dp*1);
  const float dpmin2=dpmin*dpmin;
  const unsigned newnp=cusphinout::InOutFillListCreate(Stable,PtCount,PtPosxyg,PtPoszg
    ,PtZoneg,CfgUpdateg,Zsurfg,Widthg,inoutcount,prodistg,proposxyg,proposzg
    ,dpmin,dpmin2,float(Dp),PtAuxDistg,sizenp-1,(unsigned*)inoutpartg);

  //-Creates new inlet/outlet particles to fill inlet/outlet domain.
  cusphinout::InOutFillCreate(PeriActive,newnp,(unsigned *)inoutpartg,PtPosxyg,PtPoszg,PtZoneg,PtAuxDistg
    ,np,idnext,CodeNewPart,DirDatag,posxyg,poszg,dcellg,codeg,idpg,velrhopg);
  return(newnp);
}
#endif

//==============================================================================
/// Updates velocity and rhop for M1 variable when Verlet is used. 
/// Actualiza velocidad y densidad de varible M1 cuando se usa Verlet.
//==============================================================================
void JSphInOut::UpdateVelrhopM1Cpu(unsigned inoutcount,const int *inoutpart
  ,const tfloat4 *velrhop,tfloat4 *velrhopm1)
{
  const int ncp=int(inoutcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(ncp>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int cp=0;cp<ncp;cp++){
    const unsigned p=(unsigned)inoutpart[cp];
    velrhopm1[p]=velrhop[p];
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Updates velocity and rhop for M1 variable when Verlet is used. 
/// Actualiza velocidad y densidad de varible M1 cuando se usa Verlet.
//==============================================================================
void JSphInOut::UpdateVelrhopM1Gpu(unsigned inoutcount,const int *inoutpartg
  ,const float4 *velrhopg,float4 *velrhopm1g)
{
  cusphinout::InOutUpdateVelrhopM1(inoutcount,inoutpartg,velrhopg,velrhopm1g);
}
#endif

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSphInOut::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  Log->Printf("MemoryResizeNp: +%u particles (initial: +%u particles)",GetNpResizePlus1(),GetNpResizePlus0());
  Log->Printf("UseBoxLimit: %s",(UseBoxLimit? "True": "False"));
  if(UseBoxLimit){
    Log->Printf("  FreeLimits:%s",fun::Float3xRangeStr(FreeLimitMin,FreeLimitMax,"%g").c_str());
    Log->Printf("  FreeCentre:(%s)",fun::Float3gStr(FreeCentre).c_str());
  }
  Log->Printf("DetermLimit.: %g %s",DetermLimit,(DetermLimit==1e-3f? "(1st order)": (DetermLimit==1e+3f? "(0th order)": " ")));
  Log->Printf("ExtrapolateMode: %s",(ExtrapolateMode==1? "FastSingle": (ExtrapolateMode==2? "Single": (ExtrapolateMode==3? "Double": "???"))));
  for(unsigned ci=0;ci<GetCount();ci++){
    JSphInOutZone *izone=List[ci];
    Log->Printf("InOut_%u",izone->GetIdZone());
    std::vector<std::string> lines;
    izone->GetConfig(lines);
    for(unsigned i=0;i<unsigned(lines.size());i++)Log->Print(string("  ")+lines[i]);
    izone->CheckConfig();
  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Saves VTK and CSV files per PART.
//==============================================================================
void JSphInOut::SavePartFiles(unsigned part){
  //-Creates VTK file with Zsurf.
  SaveVtkZsurf(part);
}

//==============================================================================
/// Creates VTK files with Zsurf.
//==============================================================================
void JSphInOut::SaveVtkZsurf(unsigned part){
  JVtkLib sh;
  bool usesh=false;
  for(unsigned ci=0;ci<GetCount();ci++){
    const JSphInOutZone *izone=List[ci];
    if(izone->GetSvVtkZsurf()){
      const float zsurf=izone->GetInputZsurf();
      //const tdouble3* ptdom=izone->GetPtDomain();
      tfloat3 boxmin=List[ci]->GetBoxLimitMin();
      tfloat3 boxmax=List[ci]->GetBoxLimitMax();
      if(Simulate2D){
        const float py=float(Simulate2DPosY);
        const tfloat3 pt1=TFloat3(boxmin.x,py,zsurf);
        const tfloat3 pt2=TFloat3(boxmax.x,py,zsurf);
        sh.AddShapeLine(pt1,pt2,ci);
      }
      else{
        boxmin.z=boxmax.z=zsurf;
        const tfloat3 pt1=TFloat3(boxmax.x,boxmin.y,boxmin.z);
        const tfloat3 pt2=TFloat3(boxmin.x,boxmax.y,boxmax.z);
        sh.AddShapeQuad(boxmin,pt1,boxmax,pt2,ci);
      }
      usesh=true;
    }
  }
  if(usesh){
    const string filevtk=AppInfo.GetDirDataOut()+"InOut_Zsurf.vtk";
    sh.SaveShapeVtk(fun::FileNameSec(filevtk,part),"izone");
    Log->AddFileInfo(filevtk,"Saves VTK files with Zsurf (by JSphInOut).");
  }
}

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
double JSphInOut::GetDistPtzPos(unsigned ci)const{
  return(ci<ListSize? List[ci]->GetDistPtzPos(): 0); 
}

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
unsigned JSphInOut::GetCountPtzPos(unsigned ci)const{  
  return(ci<ListSize? List[ci]->GetCountPtzPos(): 0); 
}

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
const tfloat3* JSphInOut::GetPtzPos(unsigned ci)const{ 
  return(ci<ListSize? List[ci]->GetPtzPos(): NULL); 
}

#ifdef _WITHGPU
//==============================================================================
/// Returns data from requested zone.
//==============================================================================
const float3* JSphInOut::GetPtzPosg(unsigned ci)const{     
  return(ci<ListSize? List[ci]->GetPtzPosg(): NULL); 
}

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
float* JSphInOut::GetPtzAuxg(unsigned ci)const{ 
  return(ci<ListSize? List[ci]->GetPtzAuxg(): NULL); 
}

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
float* JSphInOut::GetPtzAux (unsigned ci)const{ 
  return(ci<ListSize? List[ci]->GetPtzAux(): NULL); 
}
#endif

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
void JSphInOut::SetInputZsurf(unsigned ci,float zsurf){ 
  if(ci<ListSize)List[ci]->SetInputZsurf(zsurf); 
}

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
float JSphInOut::GetZbottom(unsigned ci)const{ 
  return(ci<ListSize? List[ci]->GetInputZbottom(): 0); 
}

#ifdef _WITHGPU
//==============================================================================
/// Returns data from requested zone.
//==============================================================================
byte JSphInOut::GetExtrapRhopMask()const{ 
  return(byte(JSphInOutZone::MRHOP_Extrapolated)); 
}

//==============================================================================
/// Returns data from requested zone.
//==============================================================================
byte JSphInOut::GetExtrapVelMask()const{ 
  return(byte(JSphInOutZone::MVEL_Extrapolated)); 
}
#endif


