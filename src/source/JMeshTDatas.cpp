//HEAD_DSCODES
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

/// \file JMeshTDatas.cpp \brief Implements the classes JMeshTDatas.

#include "JMeshTDatas.h"
#include "Functions.h"

#include "JMeshTDatasLoad.h"
#include "JMeshTDatasSave.h"
#include "JMeshData.h"
#include "JDataArrays.h"
#include "JSaveCsv2.h"
#include "FunGeo3d.h"
#include "JTimeControl.h"


#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <climits>
#include <algorithm>

using namespace std;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshTDatas
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMeshTDatas::JMeshTDatas(std::string appname,bool printtime)
  :AppName(appname),PrintTime(printtime)
{
  ClassName="JMeshTDatas";
  Times=NULL;
  OutPos=NULL;
  PosCell1=NULL;  PosCell2=NULL;  PosCell3=NULL;  PosFdis=NULL;
  Arrays1=NULL;  Arrays2=NULL;
  FrData=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JMeshTDatas::~JMeshTDatas(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JMeshTDatas::Reset(){
  FileIn="";
  ClearDatas();
  LoopTmax=DBL_MAX; 
  LoopTsub=0;
  LoopTbeg=DBL_MAX; 
}

////==============================================================================
///// Devuelve la memoria reservada.
///// Returns allocated memory
////==============================================================================
//llong JMeshTDatas::GetAllocMemory()const{  
//  return(0);
//}

//==============================================================================
/// Reset vector Datas and interpolation variables.
//==============================================================================
void JMeshTDatas::ClearDatas(){
  const unsigned sdatas=GetCount();
  for(unsigned c=0;c<sdatas;c++)delete Datas[c];
  Datas.clear();
  IntpReset();
}

//==============================================================================
/// Initialise data about times.
//==============================================================================
void JMeshTDatas::ResetTimes(){
  SizeTimes=0;
  delete[] Times; Times=NULL;
  TimeStep=TimePre=TimeNext=TimeFactor=0;
  Position=PositionNext=UINT_MAX;
}

//==============================================================================
/// Prepares times.
//==============================================================================
void JMeshTDatas::PrepareTimes(){
  ResetTimes();
  SizeTimes=GetCount();
  Times=new double[SizeTimes];
  for(unsigned ct=0;ct<SizeTimes;ct++)Times[ct]=Datas[ct]->GetTimeStep();
}

//==============================================================================
/// Reset vector vdata.
//==============================================================================
void JMeshTDatas::IntpReset(){
  ResetTimes();
  //-Variables for interpolation I.
  InPrepared=false;
  memset(&InMesh,0,sizeof(StMeshPts));
  Dir1=Dir2=Dir3=false;
  InMode=InMode12=0;
  Npt1=Npt2=Npt3=Npt12=Npt=0;
  Pla1=Pla2=Pla3=TPlane3d(0);
  Dp=TDouble3(0);
  Cmax=TUint3(0);
  //-Variables for interpolation II.
  InConfigured=IntpNone;
  IntpFreeMemory();
  //-Variables for interpolation III.
  FrTimeStep=FrTimeFactor=DBL_MAX;
  FrCtime=0;
  FrMode12=false;
  IntpAllocFrMemory(0,TypeNull);
  FrMode=99;
  FrCmax1=FrCmax2=FrCmax3=0;
  FrPla1=FrPla2=FrPla3=TPlane3d(0);
  FrNum1=FrNum2=0;
}

//==============================================================================
/// Compute plane to calculate distance to ptref on one arbitrary axis.
//==============================================================================
tplane3d JMeshTDatas::ComputeAxisPlane(const tdouble3& ps,const tdouble3& vdp)const{
  //tplane3d pla1=fgeo::PlanePtVec(ps,vdp);
  //double dis=sqrt(pla1.a*pla1.a+pla1.b*pla1.b+pla1.c*pla1.c)*fgeo::PointDist(vdp);
  //return(TPlane3d(pla1.a/dis,pla1.b/dis,pla1.c/dis,pla1.d/dis));
  return(fgeo::PlaneAxisDist(ps,vdp,fgeo::PointDist(vdp)));
}

//==============================================================================
/// Prepares basic variables for interpolation.
//==============================================================================
void JMeshTDatas::IntpPrepare(){
  //printf("==> IntpPrepare>> \n");
  InMesh=GetMeshPt();
  const StMeshPts& m=InMesh;
  Dir1=(m.npt1>1);
  Dir2=(m.npt2>1);
  Dir3=(m.npt3>1);

  InMode=InMode12=0;
  if(m.npt==1)                    { InMode=10; InMode12=10; } //-Dirs:(0,0,0)
  else if( Dir1 && !Dir2 && !Dir3){ InMode=1;  InMode12=1;  } //-Dirs:(1,0,0)
  else if(!Dir1 &&  Dir2 && !Dir3){ InMode=2;  InMode12=2;  } //-Dirs:(0,1,0)
  else if(!Dir1 && !Dir2 &&  Dir3){ InMode=3;  InMode12=3;  } //-Dirs:(0,0,1)
  else if(!Dir1 &&  Dir2 &&  Dir3){ InMode=13; InMode12=2;  } //-Dirs:(0,1,1)
  else if( Dir1 && !Dir2 &&  Dir3){ InMode=15; InMode12=1;  } //-Dirs:(1,0,1)
  else if( Dir1 &&  Dir2 && !Dir3){ InMode=16; InMode12=1;  } //-Dirs:(1,1,0)
  else if( Dir1 &&  Dir2 &&  Dir3){ InMode=17; InMode12=16; } //-Dirs:(1,1,1)
  else Run_Exceptioon("Mode interpolation is invalid.");

  Npt1=m.npt1;
  Npt2=m.npt2;
  Npt3=m.npt3;
  Npt12=Npt1*Npt2;
  Npt=m.npt;
  Pla1=(Dir1? ComputeAxisPlane(m.ptref,m.vdp1): TPlane3d(0));
  Pla2=(Dir2? ComputeAxisPlane(m.ptref,m.vdp2): TPlane3d(0));
  Pla3=(Dir3? ComputeAxisPlane(m.ptref,m.vdp3): TPlane3d(0));
  Dp.x=(Dir1? fgeo::PointDist(m.vdp1): 0);
  Dp.y=(Dir2? fgeo::PointDist(m.vdp2): 0);
  Dp.z=(Dir3? fgeo::PointDist(m.vdp3): 0);
  Cmax=TUint3(m.npt1-1,m.npt2-1,m.npt3-1);
  InPrepared=true;
}

//==============================================================================
/// Free allocated memory for interpolation.
//==============================================================================
void JMeshTDatas::IntpFreeMemory(){
  InConfigured=IntpNone;
  Npos12=NposTot=0;
  memset(&OutMesh,0,sizeof(StMeshPts));
  delete[] OutPos;    OutPos=NULL;
  delete[] PosCell1;  PosCell1=NULL;
  delete[] PosCell2;  PosCell2=NULL;
  delete[] PosCell3;  PosCell3=NULL;
  delete[] PosFdis;   PosFdis=NULL;
  Ctime1=UINT_MAX;
  delete Arrays1;  Arrays1=NULL;
  Ctime2=UINT_MAX;
  delete Arrays2;  Arrays2=NULL;
}

//==============================================================================
/// Define configuration and allocate dynamic memory for interpolation.
//==============================================================================
void JMeshTDatas::IntpAllocMemory(TpIntpMode mode,unsigned npos12,unsigned npostot){
  IntpFreeMemory();
  InConfigured=mode;
  Npos12=npos12;
  NposTot=npostot;
  //printf("==> IntpAllocMemory> Npos12:%u  NposTot:%u \n",Npos12,NposTot);
  //-Allocates basic memory for interpolation.
  try{
    if(InConfigured==IntpPos)OutPos=new tdouble3 [NposTot];
    if(Dir1)PosCell1=new tuint2 [NposTot];
    if(Dir2)PosCell2=new tuint2 [NposTot];
    if(Dir3)PosCell3=new tuint2 [NposTot];
    PosFdis=new tfloat3[NposTot];
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
  //-Allocates and prepare auxiliary arrays.
  Arrays1=new JDataArrays();
  Arrays2=new JDataArrays();
  IntpPrepareArray(Npos12,NposTot,Arrays1);
  IntpPrepareArray(Npos12,NposTot,Arrays2);
}

//==============================================================================
/// Prepare auxiliary arrays for position-interpolation.
//==============================================================================
void JMeshTDatas::IntpPrepareArray(unsigned npos12,unsigned npostot
  ,JDataArrays* arrays)const
{
  if(!GetCount())Run_Exceptioon("No data available for interplation configuration.");
  const bool resetdata=false;
  const JDataArrays* ars0=Datas[0]->GetArrays();
  const unsigned na=ars0->Count();
  for(unsigned ca=0;ca<na;ca++){
    const JDataArrays::StDataArray& ar0=ars0->GetArrayCte(ca);
    const bool dat12=(ar0.tag==JMeshData::TAGDATA12);
    const unsigned size=(dat12? npos12: npostot);
    const unsigned idx=arrays->Count();
    switch(ar0.type){
      case TypeFloat:    arrays->CreateArrayFloat  (ar0.fullname,size,resetdata);  break;
      case TypeFloat3:   arrays->CreateArrayFloat3 (ar0.fullname,size,resetdata);  break;
      case TypeDouble:   arrays->CreateArrayDouble (ar0.fullname,size,resetdata);  break;
      case TypeDouble3:  arrays->CreateArrayDouble3(ar0.fullname,size,resetdata);  break;
      default: Run_Exceptioon(fun::PrintStr("Type of data \'%s\' is invalid.",ar0.keyname.c_str()));
    }
    arrays->GetArray(idx).tag=ar0.tag;
  }
}

//==============================================================================
/// Computes cell and factor position according to given data.
//==============================================================================
float JMeshTDatas::IntpComputePosCellDist(const tdouble3& ps,const tplane3d& pla
  ,unsigned cmax,tuint2& cell)const
{
  float fdis=0;
  double d=(ps.x*pla.a + ps.y*pla.b + ps.z*pla.c + pla.d);
  d=max(d,0.);
  const unsigned c=unsigned(d);
  if(c>=cmax)cell=TUint2(cmax);
  else{
    cell=TUint2(c,c+1);
    fdis=float(d-c);
  }
  return(fdis);
}

//==============================================================================
/// Select previous and following times with data.
//==============================================================================
void JMeshTDatas::FindTime(const double timestep){
  if(!SizeTimes)PrepareTimes();
  if(!SizeTimes)Run_Exceptioon("There are not times.");
  if(timestep!=TimeStep || Position==UINT_MAX){
    const double timesteploop=(timestep>=LoopTmax? fmod((timestep-LoopTbeg),LoopTsub)+LoopTbeg: timestep);
    const double timestep2=(timesteploop>=LoopTmax? timesteploop-LoopTsub: timesteploop);
    unsigned pos=(Position==UINT_MAX? 0: Position);
    unsigned posnext=pos;
    double tpre=Times[pos];
    double tnext=tpre;
    if(SizeTimes>1){
      while(tpre>=timestep2 && pos>0){//-Retrocede.
        pos--;
        tpre=Times[pos];
      }
      posnext=(pos+1<SizeTimes? pos+1: pos);
      tnext=Times[posnext];
      while(tnext<timestep2 && posnext+1<SizeTimes){//-Avanza.
        posnext++;
        tnext=Times[posnext];
      }
      if(posnext-pos>1){
        pos=posnext-1;
        tpre=Times[pos];
      }
    }
    TimeStep=timestep;
    Position=pos;
    PositionNext=posnext;
    TimePre=tpre;
    TimeNext=tnext;
    const double tdif=TimeNext-TimePre;
    TimeFactor=(tdif? (timestep2-TimePre)/tdif: 0);
    if(TimeFactor<0.)TimeFactor=0.;
    if(TimeFactor>1.)TimeFactor=1.;
  }
}

//==============================================================================
/// Compute interpolated values in position according to ctime.
//==============================================================================
template<class T> void JMeshTDatas::IntpComputePosVar(unsigned ctime,bool data12
  ,const T* data,unsigned npos,T* res)
{
  const byte inmode=(data12? InMode12: InMode);
  //printf("===> IntpComputePosVar> ct:%u data12:%d npos:%u inmode:%d\n",ctime,(data12?1:0),npos,inmode);
  //printf("===> Dir:%d,%d,%d\n",(Dir1?1:0),(Dir2?1:0),(Dir3?1:0));
  //if(ctime>0)throw "aaaa";
  switch(inmode){
    case 1: //-Dirs:(1,0,0)
      for(unsigned cp=0;cp<npos;cp++){
        const float fdis=PosFdis[cp].x;
        const tuint2 cell=PosCell1[cp];
        const T v0=data[cell.x];
        const T v1=data[cell.y];
        res[cp]=(v1-v0)*fdis+v0;
        //printf("   ctime:%d> fdis[%d]:%f  v0_1:%f_%f  v:%f\n",ctime,cp,fdis,*((float*)&v0),*((float*)&v1),*((float*)(res+cp)));
      }
    break;
    case 2: //-Dirs:(0,1,0)
      for(unsigned cp=0;cp<npos;cp++){
        const float fdis=PosFdis[cp].y;
        const tuint2 cell=PosCell2[cp];
        const T v0=data[cell.x];
        const T v1=data[cell.y];
        res[cp]=(v1-v0)*fdis+v0;
        //printf("   ctime:%d> fdis[%d]:%f  v0_1:%f_%f  v:%f\n",ctime,cp,fdis,v0,v1,res[cp]);
      }
    break;
    case 3: //-Dirs:(0,0,1)
      if(data12 && InMode==3)for(unsigned cp=0;cp<npos;cp++)res[cp]=data[0]; //-Solves special case.
      else for(unsigned cp=0;cp<npos;cp++){
        const float fdis=PosFdis[cp].z;
        const tuint2 cell=PosCell3[cp];
        const T v0=data[cell.x];
        const T v1=data[cell.y];
        res[cp]=(v1-v0)*fdis+v0;
        //printf("   ctime:%d> fdis[%d]:%f  v0_1:%f_%f  v:%f\n",ctime,cp,fdis,v0,v1,res[cp]);
      }
    break;
    case 13: //-Dirs:(0,1,1)
      for(unsigned cp=0;cp<npos;cp++){
        const tfloat3 fdis=PosFdis[cp];
        const tuint2 cell2=PosCell2[cp];
        const tuint2 cell3=PosCell3[cp];
        const T v000=data[cell2.x + Npt12*cell3.x];
        const T v010=data[cell2.y + Npt12*cell3.x];
        const T v001=data[cell2.x + Npt12*cell3.y];
        const T v011=data[cell2.y + Npt12*cell3.y];
        const T v0x0=(v010-v000)*fdis.y+v000;
        const T v0x1=(v011-v001)*fdis.y+v001;
        res[cp]=(v0x1-v0x0)*fdis.z+v0x0;
        //printf("   ctime:%d> fdis[%d]:%f  v0_1:%f_%f  v:%f\n",ctime,cp,fdis,v0,v1,res[cp]);
      }
    break;
    case 15: //-Dirs:(1,0,1)
      for(unsigned cp=0;cp<npos;cp++){
        const tfloat3 fdis=PosFdis[cp];
        const tuint2 cell1=PosCell1[cp];
        const tuint2 cell3=PosCell3[cp];
        const T v000=data[cell1.x + Npt12*cell3.x];
        const T v100=data[cell1.y + Npt12*cell3.x];
        const T v001=data[cell1.x + Npt12*cell3.y];
        const T v101=data[cell1.y + Npt12*cell3.y];
        const T vx00=(v100-v000)*fdis.x+v000;
        const T vx01=(v101-v001)*fdis.x+v001;
        res[cp]=(vx01-vx00)*fdis.z+vx00;
        //printf("   ctime:%d> fdis[%d]:%f  v0_1:%f_%f  v:%f\n",ctime,cp,fdis,v0,v1,res[cp]);
      }
    break;
    case 16: //-Dirs:(1,1,0)
      for(unsigned cp=0;cp<npos;cp++){
        const tfloat3 fdis=PosFdis[cp];
        const tuint2 cell1=PosCell1[cp];
        const tuint2 cell2=PosCell2[cp];
        const T v000=data[cell1.x + Npt1*cell2.x];
        const T v100=data[cell1.y + Npt1*cell2.x];
        const T v010=data[cell1.x + Npt1*cell2.y];
        const T v110=data[cell1.y + Npt1*cell2.y];
        const T vx00=(v100-v000)*fdis.x+v000;
        const T vx10=(v110-v010)*fdis.x+v010;
        res[cp]=(vx10-vx00)*fdis.y+vx00;
        //printf("   ctime:%d> fdis[%d]:%f  v0_1:%f_%f  v:%f\n",ctime,cp,fdis,v0,v1,res[cp]);
      }
    break;
    case 17: //-Dirs:(1,1,1)
      for(unsigned cp=0;cp<npos;cp++){
        const tfloat3 fdis=PosFdis[cp];
        const tuint2 cell1=PosCell1[cp];
        const tuint2 cell2=PosCell2[cp];
        const tuint2 cell3=PosCell3[cp];
        T res3x,res3y;
        for(unsigned c3=0;c3<2;c3++){
          const unsigned mod3=Npt12*(c3? cell3.y: cell3.x);
          const T v000=data[cell1.x + Npt1*cell2.x + mod3];
          const T v100=data[cell1.y + Npt1*cell2.x + mod3];
          const T v010=data[cell1.x + Npt1*cell2.y + mod3];
          const T v110=data[cell1.y + Npt1*cell2.y + mod3];
          const T vx00=(v100-v000)*fdis.x+v000;
          const T vx10=(v110-v010)*fdis.x+v010;
          if(c3)res3y=(vx10-vx00)*fdis.y+vx00;
          else  res3x=(vx10-vx00)*fdis.y+vx00;
        }
        res[cp]=(res3y-res3x)*fdis.z+res3x;
        //printf("ct:%d [%d]  ce1:%u-%u  ce2:%u-%u  ce3:%u-%u\n",ctime,cp,cell1.x,cell1.y,cell2.x,cell2.y,cell3.x,cell3.y);
      }
    break;
    case 10:{ //-Dirs:(1,1,1)
      const T v=data[0];
      for(unsigned cp=0;cp<npos;cp++)res[cp]=v;
    }break;
    default: Run_Exceptioon("Mode interpolation is invalid.");
  }
}

//==============================================================================
/// Compute position-interpolation of arrays for time interval ctime.
//==============================================================================
void JMeshTDatas::IntpComputePosArray(unsigned ctime,JDataArrays* ars){
  //printf("===>  IntpComputePosArray> ctime:%u\n",ctime);
  const JDataArrays* ardata=Datas[ctime]->GetArrays();
  const unsigned na=ars->Count();
  for(unsigned ca=0;ca<na;ca++){
    const void* ptrdata=ardata->GetArrayCte(ca).ptr;
    const JDataArrays::StDataArray& ar=ars->GetArray(ca);
    const bool data12=(ar.tag==JMeshData::TAGDATA12);
    //if(ar.tag==JMeshData::TAGDATA12)Run_Exceptioon("TAGDATA12: Not implemented yet!!");
    //printf("======> ca:%u [%s]\n",ca,ar.keyname.c_str());
    switch(ar.type){
      case TypeFloat:   IntpComputePosVar(ctime,data12,(const float*   )ptrdata,ar.count,(float*   )ar.ptr);  break;
      case TypeFloat3:  IntpComputePosVar(ctime,data12,(const tfloat3* )ptrdata,ar.count,(tfloat3* )ar.ptr);  break;
      case TypeDouble:  IntpComputePosVar(ctime,data12,(const double*  )ptrdata,ar.count,(double*  )ar.ptr);  break;
      case TypeDouble3: IntpComputePosVar(ctime,data12,(const tdouble3*)ptrdata,ar.count,(tdouble3*)ar.ptr);  break;
      default: Run_Exceptioon(fun::PrintStr("Type of data \'%s\' is invalid.",ar.keyname.c_str()));
    }
  }
  ////-Save VTK with interpolated data. Only for debug.
  //for(unsigned cf=0;cf<2;cf++){
  //  const int seltag=(cf? JMeshData::TAGDATA12: 0);
  //  JDataArrays arrays;
  //  arrays.CopyFrom(*ars);
  //  const unsigned na=ars->Count();
  //  for(unsigned ca=0;ca<na;ca++){
  //    const JDataArrays::StDataArray& ar=ars->GetArray(ca);
  //    if(ar.tag!=seltag)arrays.DeleteArray(ar.keyname);
  //  }
  //  if(arrays.Count()){
  //    const unsigned np=arrays.GetDataCount();
  //    tfloat3* pos=arrays.CreateArrayPtrFloat3("Pos",np,false);
  //    if(InConfigured==IntpMesh)JMeshData::GetMeshPos(OutMesh,np,pos);
  //    else{
  //      if(np>NposTot || !OutPos)Run_Exceptioon("OutPos[] is invalid.");
  //      for(unsigned p=0;p<np;p++)pos[p]=ToTFloat3(OutPos[p]);
  //    }
  //    if(cf && arrays.GetArrayCte(0).keyname=="Zsurf"){
  //      const float* zs=arrays.GetArrayFloat("Zsurf");
  //      for(unsigned p=0;p<np;p++)pos[p].z=zs[p];
  //    }
  //    const string file=string("DgOut/_DG_InterpolatedData")+(cf? "_Zsurf.vtk": ".vtk");
  //    JSpVtkData::Save(fun::FileNameSec(file,ctime),arrays,"Pos");
  //  }
  //}
}

//==============================================================================
/// Compute interpolated values according to positions in time t.
//==============================================================================
void JMeshTDatas::IntpComputePosIntervals(double t){
  //printf("=> IntpComputePosIntervals>> NposTot:%u Npos12:%u\n",NposTot,Npos12);
  //-Select previous and following times with data.
  FindTime(t);
  //-Computes interpolation data in position for previous and next times.
  //printf("=> t:%f  ctimes:%u-%u  times:%f-%f  tf:%f\n",t,Position,PositionNext,TimePre,TimeNext,TimeFactor);
  if(Ctime1!=Position || Ctime2!=PositionNext){
    //-Swap data when it is possible.
    if(Ctime2==Position || Ctime1==PositionNext){
      swap(Ctime1,Ctime2);
      swap(Arrays1,Arrays2);
    }
    if(Ctime1!=Position){
      IntpComputePosArray(Position,Arrays1);
      Ctime1=Position;
    }
    if(Ctime2!=PositionNext){
      if(Ctime1==PositionNext){//-Copy from datatime1.
        const unsigned na=Arrays2->Count();
        for(unsigned ca=0;ca<na;ca++){
          const JDataArrays::StDataArray& ar2=Arrays2->GetArray(ca);
          const size_t sizetot=TypeSize(ar2.type)*ar2.count;
          const void* ptr1=(const void*)Arrays1->GetArrayCte(ca).ptr;
          memcpy(ar2.ptr,ptr1,sizetot);
        }
      }
      else IntpComputePosArray(PositionNext,Arrays2);
      Ctime2=PositionNext;
    }
  }
}

//==============================================================================
/// Returns variable number after checking requested variable and interpolation 
/// configuration. Throws exceptions in case of error.
//==============================================================================
unsigned JMeshTDatas::IntpCheckVarName(const std::string& varname,unsigned size)const{
  if(InConfigured!=IntpMesh && InConfigured!=IntpPos)Run_Exceptioon("Fixed interplation configuration is missing.");
  const unsigned varnum=Arrays1->GetIdxName(varname);
  if(varnum==UINT_MAX)Run_Exceptioon(fun::PrintStr("The variable \'%s\' is missing.",varname.c_str()));
  if(size!=Arrays1->GetArrayCte(varnum).count)Run_Exceptioon(fun::PrintStr("Size of variable \'%s\' does not match.",varname.c_str()));
  return(varnum);
}

//==============================================================================
/// Checks requested variable and interpolation configuration. 
/// Throws exceptions in case of error.
//==============================================================================
unsigned JMeshTDatas::IntpCheckVarNum(unsigned varnum,unsigned size)const{
  if(InConfigured!=IntpMesh && InConfigured!=IntpPos)Run_Exceptioon("Fixed interplation configuration is missing.");
  if(varnum>=Arrays1->Count())Run_Exceptioon("Number of variable is invalid.");
  if(size!=Arrays1->GetArrayCte(varnum).count)Run_Exceptioon("Size of requested variable does not match.");
  return(varnum);
}

//==============================================================================
/// Computes interpolation on requeseted variable.
//==============================================================================
template<class T,class T1> void JMeshTDatas::IntpComputeVarT(double timestep
  ,unsigned varnum,T* result)
{
  IntpComputePosIntervals(timestep);
  if(TimeFactor==0 || TimeFactor==1.f){//-Copies data from previous or next time-data. 
    const JDataArrays::StDataArray& ar=(TimeFactor==0? Arrays1: Arrays2)->GetArrayCte(varnum);
    const T* pdat=(const T*)ar.ptr;
    memcpy(result,(const T*)ar.ptr,sizeof(T)*ar.count);
  }
  else{//-Interpolates result between previous next time-data. 
    const JDataArrays::StDataArray& ar1=Arrays1->GetArrayCte(varnum);
    const T* pv1=(const T*)ar1.ptr;
    const T* pv2=(const T*)Arrays2->GetArrayCte(varnum).ptr;
    for(unsigned c=0;c<ar1.count;c++){
      const T v1=pv1[c];
      result[c]=((pv2[c]-v1)*T1(TimeFactor)+v1);
    }
  }
}
//==============================================================================
void JMeshTDatas::IntpComputeVar(double timestep,unsigned varnum,unsigned size,float*    result){  IntpComputeVarT<float   ,float >(timestep,IntpCheckVarNum(varnum,size),result);  }
void JMeshTDatas::IntpComputeVar(double timestep,unsigned varnum,unsigned size,tfloat3*  result){  IntpComputeVarT<tfloat3 ,float >(timestep,IntpCheckVarNum(varnum,size),result);  }
void JMeshTDatas::IntpComputeVar(double timestep,unsigned varnum,unsigned size,double*   result){  IntpComputeVarT<double  ,double>(timestep,IntpCheckVarNum(varnum,size),result);  }
void JMeshTDatas::IntpComputeVar(double timestep,unsigned varnum,unsigned size,tdouble3* result){  IntpComputeVarT<tdouble3,double>(timestep,IntpCheckVarNum(varnum,size),result);  }
//==============================================================================
void JMeshTDatas::IntpComputeVar(double timestep,const std::string& varname,unsigned size,float*    result){  IntpComputeVarT<float   ,float >(timestep,IntpCheckVarName(varname,size),result);  }
void JMeshTDatas::IntpComputeVar(double timestep,const std::string& varname,unsigned size,tfloat3*  result){  IntpComputeVarT<tfloat3 ,float >(timestep,IntpCheckVarName(varname,size),result);  }
void JMeshTDatas::IntpComputeVar(double timestep,const std::string& varname,unsigned size,double*   result){  IntpComputeVarT<double  ,double>(timestep,IntpCheckVarName(varname,size),result);  }
void JMeshTDatas::IntpComputeVar(double timestep,const std::string& varname,unsigned size,tdouble3* result){  IntpComputeVarT<tdouble3,double>(timestep,IntpCheckVarName(varname,size),result);  }


//==============================================================================
/// Loads data from a file (.mbi4 or .csv).
//==============================================================================
void JMeshTDatas::LoadFile(const std::string file,std::string varlist
  ,double tmin,double tmax,double looptmax,double looptbegin
  ,tdouble3 pmin,tdouble3 pmax,double settime,tdouble3 setpos)
{
  Reset();
  FileIn=file;
  JMeshTDatasLoad fload(tmin,tmax,looptmax,looptbegin,pmin,pmax,settime,setpos);
  fload.LoadFile(file,Datas,varlist);
  LoopTmax=fload.GetLoopTmax();
  LoopTsub=fload.GetLoopTsub();
  LoopTbeg=LoopTmax-LoopTsub;
}

//==============================================================================
/// Add mesh Data-time.
//==============================================================================
void JMeshTDatas::AddData(JMeshData* mdat){
  if(InPrepared || InConfigured!=IntpNone)IntpReset();
  Datas.push_back(mdat);
}

//==============================================================================
/// Reverse selected variables (v=-v).
//==============================================================================
void JMeshTDatas::ReverseData(std::string varlist){
  const unsigned ndata=GetCount();
  for(unsigned c=0;c<ndata;c++)Datas[c]->ReverseData(varlist);
}

//==============================================================================
/// Set data of selected variables "all,name1,name1.x,name1.y,name2.z,name3" (v*=v2).
//==============================================================================
void JMeshTDatas::SetMulData(std::string varlist,double v2){
  const unsigned ndata=GetCount();
  for(unsigned c=0;c<ndata;c++)Datas[c]->SetMulData(varlist,v2);
}

//==============================================================================
/// Set data of selected variables "all,name1,name1.x,name1.y,name2.z,name3" (v+=v2).
//==============================================================================
void JMeshTDatas::SetAddData(std::string varlist,double v2){
  const unsigned ndata=GetCount();
  for(unsigned c=0;c<ndata;c++)Datas[c]->SetAddData(varlist,v2);
}

//==============================================================================
/// Returns requested time.
//==============================================================================
double JMeshTDatas::GetTime(unsigned ct)const{
  const unsigned ntimes=GetCount();
  if(ct>=GetCount())Run_Exceptioon("No data available for requested time.");
  return(Datas[ct]->GetTimeStep());
}

//==============================================================================
/// Save times in vtimes[] and returns number of times.
//==============================================================================
unsigned JMeshTDatas::GetTimes(std::vector<double>& vtimes)const{
  const unsigned ntimes=GetCount();
  vtimes.clear();
  vtimes.reserve(ntimes);
  for(unsigned ct=0;ct<ntimes;ct++)vtimes.push_back(Datas[ct]->GetTimeStep());
  return(ntimes);
}

//==============================================================================
/// Returns current Mesh data definition.
//==============================================================================
StMeshPts JMeshTDatas::GetMeshPt()const{
  StMeshPts m;
  if(GetCount())m=Datas[0]->GetMeshPt();
  else memset(&m,0,sizeof(StMeshPts));
  return(m);
}

//==============================================================================
/// Configures interpolation for positions defined by a fixed mesh.
//==============================================================================
void JMeshTDatas::IntpConfigFixedMesh(const StMeshPts& m){
  //printf("==> IntpConfigFixedMesh>\n");
  if(!GetCount())Run_Exceptioon("No data available for interpolation configuration.");
  //-Prepares configuration.
  IntpReset();
  IntpPrepare();
  IntpAllocMemory(IntpMesh,m.npt1*m.npt2,m.npt);
  //-Copy output mesh definition.
  OutMesh=m;
  //JMeshData::PrintMeshPts("==> OutMesh:",OutMesh);
  //-Compute position info for interpolation (PosCell1[],PosCell2[],PosCell3[],PosFdis[]).
  unsigned cp=0;
  for(unsigned cp3=0;cp3<m.npt3;cp3++){
    const tdouble3 pt3=m.ptref+(m.vdp3*cp3);
    for(unsigned cp2=0;cp2<m.npt2;cp2++){
      const tdouble3 pt2=pt3+(m.vdp2*cp2);
      for(unsigned cp1=0;cp1<m.npt1;cp1++,cp++){
        const tdouble3 ps=pt2+(m.vdp1*cp1);
        tfloat3 fdis=TFloat3(0);
        if(Dir1)fdis.x=IntpComputePosCellDist(ps,Pla1,Cmax.x,PosCell1[cp]);
        if(Dir2)fdis.y=IntpComputePosCellDist(ps,Pla2,Cmax.y,PosCell2[cp]);
        if(Dir3)fdis.z=IntpComputePosCellDist(ps,Pla3,Cmax.z,PosCell3[cp]);
        PosFdis[cp]=fdis;
      }
    }
  }
}

//==============================================================================
/// Configures interpolation for a fixed list of positions.
//==============================================================================
void JMeshTDatas::IntpConfigFixedPos(unsigned npos,const tdouble3* pos){
  //printf("==> IntpConfigPos> npos:%u\n",npos);
  if(!GetCount())Run_Exceptioon("No data available for interpolation configuration.");
  //-Prepares configuration.
  IntpReset();
  IntpPrepare();
  IntpAllocMemory(IntpPos,npos,npos);
  //-Copy output pos data.
  memcpy(OutPos,pos,sizeof(tdouble3)*npos);
  //-Compute position info for interpolation (PosCell1[],PosCell2[],PosCell3[],PosFdis[]).
  //printf("==> IntpConfigPos> NposTot:%u\n",NposTot);
  for(unsigned cp=0;cp<NposTot;cp++){
    tfloat3 fdis=TFloat3(0);
    if(Dir1)fdis.x=IntpComputePosCellDist(pos[cp],Pla1,Cmax.x,PosCell1[cp]);
    if(Dir2)fdis.y=IntpComputePosCellDist(pos[cp],Pla2,Cmax.y,PosCell2[cp]);
    if(Dir3)fdis.z=IntpComputePosCellDist(pos[cp],Pla3,Cmax.z,PosCell3[cp]);
    PosFdis[cp]=fdis;
    //printf("==> IntpConfigPos> [%u]  ce1:%u-%u  ce2:%u-%u  ce3:%u-%u\n",cp,PosCell1[cp].x,PosCell1[cp].y,PosCell2[cp].x,PosCell2[cp].y,PosCell3[cp].x,PosCell3[cp].y);
  }
}

//==============================================================================
/// Configures interpolation for a variable list of positions.
//==============================================================================
void JMeshTDatas::IntpConfigFreePos(){
  //printf("==> IntpConfigFreePos>\n");
  if(!GetCount())Run_Exceptioon("No data available for free-points interpolation configuration.");
  if(Datas[0]->GetVarCount()==0)Run_Exceptioon("No variables for free-points interpolation configuration.");
  if(Datas[0]->GetVarCount()>1)Run_Exceptioon("Only one variable is allowed for free-points interpolation configuration.");
  //-Prepares configuration.
  IntpReset();
  IntpPrepare();
  InConfigured=IntpVar;
  const JDataArrays::StDataArray& ar=Datas[0]->GetArrays()->GetArrayCte(0);
  FrMode12=(ar.tag==JMeshData::TAGDATA12);
  IntpAllocFrMemory(ar.count,ar.type);
  //-Configures FrMode and parameters for interpolation.
  const byte inmode=(FrMode12? InMode12: InMode);
  //printf("==> IntpConfigFreePos> inmode:%u\n",inmode);
  switch(inmode){
    case 1: //-Dirs:(1,0,0) 
      FrMode=1;  FrCmax1=Cmax.x;  FrPla1=Pla1;
    break;      
    case 2: //-Dirs:(0,1,0)  
      FrMode=1;  FrCmax1=Cmax.y;  FrPla1=Pla2;
    break;      
    case 3: //-Dirs:(0,0,1)
      if(FrMode12 && InMode==3)FrMode=0;
      else{ FrMode=1;  FrCmax1=Cmax.z;  FrPla1=Pla3; }
    break;
    case 13: //-Dirs:(0,1,1)
      FrMode=2;  FrCmax1=Cmax.y;  FrPla1=Pla2;  FrCmax2=Cmax.z;  FrPla2=Pla3;  FrNum1=Npt12;
    break;
    case 15: //-Dirs:(1,0,1)
      FrMode=2;  FrCmax1=Cmax.x;  FrPla1=Pla1;  FrCmax2=Cmax.z;  FrPla2=Pla3;  FrNum1=Npt12;
    break;
    case 16: //-Dirs:(1,1,0)
      FrMode=2;  FrCmax1=Cmax.x;  FrPla1=Pla1;  FrCmax2=Cmax.y;  FrPla2=Pla2;  FrNum1=Npt1;
    break;
    case 17: //-Dirs:(1,1,1)
      FrMode=3;  FrCmax1=Cmax.x;  FrPla1=Pla1;  FrCmax2=Cmax.y;  FrPla2=Pla2;
      FrCmax3=Cmax.z;  FrPla3=Pla3;  FrNum1=Npt1;
    break;
    case 10: //-Dirs:(1,1,1)
        FrMode=0;
    break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
}

//==============================================================================
/// Define configuration and allocate dynamic memory for interpolation.
//==============================================================================
void JMeshTDatas::IntpAllocFrMemory(unsigned size,TpTypeData type){
  //-Free allocated memory.
  if(FrData){
    switch(FrType){
      case TypeFloat:   delete[] (float*)  FrData;  break;
      case TypeFloat3:  delete[] (tfloat3*)FrData;  break;
      default: Run_Exceptioon(fun::PrintStr("Type \'%s\' is invalid.",TypeToStr(FrType)));
    }
  }
  FrSize=0;
  FrType=TypeNull;
  FrData=NULL;
  //-Allocate memory.
  if(size){
    FrSize=size;
    FrType=type;
    try{
      switch(FrType){
        case TypeFloat:   FrData=new float  [FrSize];  break;
        case TypeFloat3:  FrData=new tfloat3[FrSize];  break;
        default: Run_Exceptioon(fun::PrintStr("Type \'%s\' is invalid.",TypeToStr(FrType)));
      }
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Computes interpolation on requeseted time when it is necessary.
//==============================================================================
template<class T> void JMeshTDatas::IntpComputeTime(double timestep){
  FindTime(timestep);
  if(FrCtime!=Position || FrTimeFactor!=TimeFactor){
    const float tfactor=float(TimeFactor);
    T* res=(T*)FrData;
    const T* dat1=(T*)Datas[Position    ]->GetArrays()->GetArrayCte(0).ptr;
    const T* dat2=(T*)Datas[PositionNext]->GetArrays()->GetArrayCte(0).ptr;
    if(Position==PositionNext || TimeFactor==0)memcpy(res,dat1,sizeof(T)*FrSize);
    else if(TimeFactor>=1)memcpy(res,dat2,sizeof(T)*FrSize);
    else for(unsigned c=0;c<FrSize;c++){
      res[c]=(dat2[c]-dat1[c])*tfactor + dat1[c];
    }
    FrCtime=Position;
    FrTimeFactor=TimeFactor;
  }
  FrTimeStep=timestep;
}

//==============================================================================
/// Compute interpolated value in position for FrMode==1.
//==============================================================================
template<class T> T JMeshTDatas::IntpComputeVarFr1(const tdouble3& ps)const{
  const T* tdat=(T*)FrData;
  tuint2 cell;
  const float fdis=IntpComputePosCellDist(ps,FrPla1,FrCmax1,cell);
  const T v0=tdat[cell.x];
  const T v1=tdat[cell.y];
  return((v1-v0)*fdis+v0);
}

//==============================================================================
/// Compute interpolated value in position for FrMode==2.
//==============================================================================
template<class T> T JMeshTDatas::IntpComputeVarFr2(const tdouble3& ps)const{
  const T* tdat=(T*)FrData;
  tuint2 cell1,cell2;
  const float fdis1=IntpComputePosCellDist(ps,FrPla1,FrCmax1,cell1);
  const float fdis2=IntpComputePosCellDist(ps,FrPla2,FrCmax2,cell2);
  const T v000=tdat[cell1.x + FrNum1*cell2.x];
  const T v010=tdat[cell1.y + FrNum1*cell2.x];
  const T v001=tdat[cell1.x + FrNum1*cell2.y];
  const T v011=tdat[cell1.y + FrNum1*cell2.y];
  const T v0x0=(v010-v000)*fdis1+v000;
  const T v0x1=(v011-v001)*fdis1+v001;
  return((v0x1-v0x0)*fdis2+v0x0);
}

//==============================================================================
/// Compute interpolated value in position for FrMode==3.
//==============================================================================
template<class T> T JMeshTDatas::IntpComputeVarFr3(const tdouble3& ps)const{
  const T* tdat=(T*)FrData;
  tuint2 cell1,cell2,cell3;
  const float fdis1=IntpComputePosCellDist(ps,FrPla1,FrCmax1,cell1);
  const float fdis2=IntpComputePosCellDist(ps,FrPla2,FrCmax2,cell2);
  const float fdis3=IntpComputePosCellDist(ps,FrPla3,FrCmax3,cell3);
  T resx,resy;
  {
    const unsigned mod3=Npt12*cell3.x;
    const T v000=tdat[cell1.x + Npt1*cell2.x + mod3];
    const T v100=tdat[cell1.y + Npt1*cell2.x + mod3];
    const T v010=tdat[cell1.x + Npt1*cell2.y + mod3];
    const T v110=tdat[cell1.y + Npt1*cell2.y + mod3];
    const T vx00=(v100-v000)*fdis1+v000;
    const T vx10=(v110-v010)*fdis1+v010;
    resx=(vx10-vx00)*fdis2+vx00;
  }
  {
    const unsigned mod3=Npt12*cell3.y;
    const T v000=tdat[cell1.x + Npt1*cell2.x + mod3];
    const T v100=tdat[cell1.y + Npt1*cell2.x + mod3];
    const T v010=tdat[cell1.x + Npt1*cell2.y + mod3];
    const T v110=tdat[cell1.y + Npt1*cell2.y + mod3];
    const T vx00=(v100-v000)*fdis1+v000;
    const T vx10=(v110-v010)*fdis1+v010;
    resy=(vx10-vx00)*fdis2+vx00;
  }
  return((resy-resx)*fdis3+resx);
}

//==============================================================================
/// Computes interpolation on requeseted time when it is necessary.
//==============================================================================
void JMeshTDatas::IntpComputeTime_f(double timestep){
  IntpComputeTime<float>(timestep);
}

//==============================================================================
/// Compute interpolated value in position for FrMode==1.
//==============================================================================
float JMeshTDatas::IntpComputeVarFr1_f(const tdouble3& ps){
  return(IntpComputeVarFr1<float>(ps));
}

//==============================================================================
/// Compute interpolated value in position for FrMode==2.
//==============================================================================
float JMeshTDatas::IntpComputeVarFr2_f(const tdouble3& ps){
  return(IntpComputeVarFr2<float>(ps));
}

//==============================================================================
/// Compute interpolated value in position for FrMode==3.
//==============================================================================
float JMeshTDatas::IntpComputeVarFr3_f(const tdouble3& ps){
  return(IntpComputeVarFr3<float>(ps));
}

//==============================================================================
/// Computes interpolation on requeseted time when it is necessary.
//==============================================================================
void JMeshTDatas::IntpComputeTime_f3(double timestep){
  IntpComputeTime<tfloat3>(timestep);
}

//==============================================================================
/// Compute interpolated value in position for FrMode==1.
//==============================================================================
tfloat3 JMeshTDatas::IntpComputeVarFr1_f3(const tdouble3& ps){
  return(IntpComputeVarFr1<tfloat3>(ps));
}

//==============================================================================
/// Compute interpolated value in position for FrMode==2.
//==============================================================================
tfloat3 JMeshTDatas::IntpComputeVarFr2_f3(const tdouble3& ps){
  return(IntpComputeVarFr2<tfloat3>(ps));
}

//==============================================================================
/// Compute interpolated value in position for FrMode==3.
//==============================================================================
tfloat3 JMeshTDatas::IntpComputeVarFr3_f3(const tdouble3& ps){
  return(IntpComputeVarFr3<tfloat3>(ps));
}

//==============================================================================
/// Computes interpolation of data on requeseted time.
//==============================================================================
void JMeshTDatas::IntpComputeInTime(double t){
  if(InConfigured!=IntpVar)Run_Exceptioon("Free interplation configuration is missing.");
  //-Select previous and following times with data.
  if(FrTimeStep!=t){
    if(FrType==TypeFloat) IntpComputeTime<float>(t);
    if(FrType==TypeFloat3)IntpComputeTime<tfloat3>(t);
  }
}

//==============================================================================
/// Computes interpolation on variable with type float.
//==============================================================================
void JMeshTDatas::IntpCompute(unsigned np,const tdouble3* pos,float* res)const{
  if(InConfigured!=IntpVar || FrType!=TypeFloat)Run_Exceptioon("Free interplation configuration is missing or data type does not match.");
  //-Interpolate in positions.
  const float* tdat=(float*)FrData;
  switch(FrMode){
    case 0: for(unsigned cp=0;cp<np;cp++)res[cp]=tdat[0];                            break;
    case 1: for(unsigned cp=0;cp<np;cp++)res[cp]=IntpComputeVarFr1<float>(pos[cp]);  break;
    case 2: for(unsigned cp=0;cp<np;cp++)res[cp]=IntpComputeVarFr2<float>(pos[cp]);  break;
    case 3: for(unsigned cp=0;cp<np;cp++)res[cp]=IntpComputeVarFr3<float>(pos[cp]);  break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
}

//==============================================================================
/// Computes interpolation point on variable with type float.
//==============================================================================
float JMeshTDatas::IntpComputeFloat(tdouble3 pos)const{
  if(InConfigured!=IntpVar || FrType!=TypeFloat)Run_Exceptioon("Free interplation configuration is missing or data type does not match.");
  //-Interpolate in positions.
  const float* tdat=(float*)FrData;
  float res=0;
  switch(FrMode){
    case 0: res=tdat[0];                        break;
    case 1: res=IntpComputeVarFr1<float>(pos);  break;
    case 2: res=IntpComputeVarFr2<float>(pos);  break;
    case 3: res=IntpComputeVarFr3<float>(pos);  break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
  return(res);
}

//==============================================================================
/// Computes interpolation on variable with type float.
//==============================================================================
void JMeshTDatas::IntpCompute(unsigned np,const tdouble3* pos,tfloat3* res)const{
  if(InConfigured!=IntpVar || FrType!=TypeFloat3)Run_Exceptioon("Free interplation configuration is missing or data type does not match.");
  //-Interpolate in positions.
  const tfloat3* tdat=(tfloat3*)FrData;
  switch(FrMode){
    case 0: for(unsigned cp=0;cp<np;cp++)res[cp]=tdat[0];                              break;
    case 1: for(unsigned cp=0;cp<np;cp++)res[cp]=IntpComputeVarFr1<tfloat3>(pos[cp]);  break;
    case 2: for(unsigned cp=0;cp<np;cp++)res[cp]=IntpComputeVarFr2<tfloat3>(pos[cp]);  break;
    case 3: for(unsigned cp=0;cp<np;cp++)res[cp]=IntpComputeVarFr3<tfloat3>(pos[cp]);  break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
}

//==============================================================================
/// Computes interpolation point on variable with type float.
//==============================================================================
tfloat3 JMeshTDatas::IntpComputeFloat3(tdouble3 pos)const{
  if(InConfigured!=IntpVar || FrType!=TypeFloat3)Run_Exceptioon("Free interplation configuration is missing or data type does not match.");
  //-Interpolate in positions.
  const tfloat3* tdat=(tfloat3*)FrData;
  tfloat3 res;
  switch(FrMode){
    case 0: res=tdat[0];                          break;
    case 1: res=IntpComputeVarFr1<tfloat3>(pos);  break;
    case 2: res=IntpComputeVarFr2<tfloat3>(pos);  break;
    case 3: res=IntpComputeVarFr3<tfloat3>(pos);  break;
    default: Run_Exceptioon("Interpolation mode is invalid.");
  }
  return(res);
}




//==============================================================================
/// Creates VTK file with the scheme of grid configuration.
//==============================================================================
void JMeshTDatas::SaveVtkScheme(std::string file)const{
  file=fun::GetWithoutExtension(file)+"_Scheme.vtk";
  if(GetCount())JMeshTDatasSave::SaveVtkScheme(file,Datas[0]->GetMeshPt());
}

//==============================================================================
/// Saves VTK files for each datatime.
//==============================================================================
void JMeshTDatas::SaveVtks(std::string file)const{
  const unsigned sdatas=GetCount();
  if(sdatas){
    file=fun::GetWithoutExtension(file);
    //-Get file name with first name of data with tag==TAGDATA12.
    string file12;
    if(sdatas){
      string namevar12=Datas[0]->GetVarList(false,true);
      namevar12=fun::StrSplitValue(",",namevar12,0);
      if(!namevar12.empty())file12=file+"_"+namevar12;
    }
    //-Saves VTK files.
    JTimeControl tc(5,60);
    for(unsigned c=0;c<sdatas;c++){
      JMeshTDatasSave::SaveVtk(file,int(c),Datas[c],false);
      if(!file12.empty())JMeshTDatasSave::SaveVtk(file12,int(c),Datas[c],true);
      if(PrintTime && tc.CheckTime())printf("  %s\n",tc.GetInfoFinish(double(c)/double(sdatas+1)).c_str());
    }
  }
}

//==============================================================================
/// Saves CSV file with all datatimes.
//==============================================================================
void JMeshTDatas::SaveCsv(std::string file,bool csvsepcoma)const{
  const unsigned sdatas=GetCount();
  if(sdatas){
    file=fun::GetWithoutExtension(file);
    const unsigned na=Datas[0]->GetVarCount();
    JTimeControl tc(5,60);
    unsigned tcnum=0,tctot=na*sdatas;
    for(unsigned ca=0;ca<na;ca++){
      const string filecsv=file+"_"+Datas[0]->GetVarName(ca)+".csv";;
      jcsv::JSaveCsv2 scsv(filecsv,false,csvsepcoma);
      for(unsigned c=0;c<sdatas;c++){
        JMeshTDatasSave::SaveCsv(Datas[c],ca,scsv,!c);
        if(scsv.GetDataSize()>18102800)scsv.SaveData(false);
        tcnum++;
        if(PrintTime && tc.CheckTime())printf("  %s\n",tc.GetInfoFinish(double(tcnum)/double(tctot+1)).c_str());
      }
      scsv.SaveData(true);
    }
  }
}

//==============================================================================
/// Saves binary file with all datatimes.
//==============================================================================
void JMeshTDatas::SaveBin(std::string file)const{
  file=fun::GetWithoutExtension(file)+".mbi4";
  const unsigned sdatas=GetCount();
  if(sdatas){
    JMeshTDatasSave fsave;
    fsave.Config(file,AppName,Datas[0]);
    for(unsigned c=0;c<sdatas;c++)fsave.SaveDataTime(Datas[c]);
  }
}

//==============================================================================
/// Saves example files of mesh-data (CSV and binary).
//==============================================================================
void JMeshTDatas::SaveTemplate(std::string filename,std::string appname,bool csvsepcoma){
  JMeshTDatas mdatas(appname);
  //-Creates mesh definition.
  StMeshPts m;
  memset(&m,0,sizeof(StMeshPts));
  m.npt1  =3;
  m.npt2  =2;
  m.npt3  =3;
  m.ptref =TDouble3(0,0,0);
  m.vdp1  =TDouble3(0.05,0,0);
  m.vdp2  =TDouble3(0,0.05,0);
  m.vdp3  =TDouble3(0,0,0.02);
  m.dirdat=fgeo::VecUnitary(TFloat3(1,0,0));
  m.npt=m.npt1*m.npt2*m.npt3;
  //-Creates mesh data for 2 instants.
  for(unsigned ct=0;ct<2;ct++){
    JMeshData* mdat=new JMeshData();
    mdat->ConfigMesh(m,0.1*ct,"Vel,VelDir,Rhop,Zsurf");
    {//-Density array.
      float* vrhop=mdat->GetVarFloat("Rhop");
      for(unsigned cp3=0,cp=0;cp3<m.npt3;cp3++)for(unsigned cp2=0;cp2<m.npt2;cp2++)for(unsigned cp1=0;cp1<m.npt1;cp1++){
        const float v=1000.f+(m.npt3-cp3-1);
        vrhop[cp++]=(!ct? v: 1000.f-v);
      }
    }
    {//-Vel and VelDir arrays.
      float* vveldir=mdat->GetVarFloat("VelDir");
      tfloat3* vvel =mdat->GetVarFloat3("Vel");
      for(unsigned cp3=0,cp=0;cp3<m.npt3;cp3++)for(unsigned cp2=0;cp2<m.npt2;cp2++)for(unsigned cp1=0;cp1<m.npt1;cp1++){
        const tfloat3 v=fgeo::VecUnitary(TFloat3(0.1f*cp1,0.2f*cp2,0.3f*cp3)*2.f)*(!ct? 1.f: -1.f);
        vvel[cp]=v;
        vveldir[cp++]=(v.x*m.dirdat.x + v.y*m.dirdat.y + v.z*m.dirdat.z);
      }
    }
    {//-Zsurf array.
      float* vzsurf=mdat->GetVarFloat("Zsurf");
      for(unsigned cp2=0,cp=0;cp2<m.npt2;cp2++)for(unsigned cp1=0;cp1<m.npt1;cp1++){
        vzsurf[cp++]=0.02f*m.npt3+0.02f*(!ct? cp1*cp2: cp1*cp1+cp2*cp2);
      }
    }
    mdatas.Datas.push_back(mdat);
  }
  //-Saves data to binary and CSV files.
  filename=fun::GetWithoutExtension(filename);
  mdatas.SaveBin(filename+".mbi4");
  mdatas.SaveCsv(filename+".csv",csvsepcoma);
}

}
