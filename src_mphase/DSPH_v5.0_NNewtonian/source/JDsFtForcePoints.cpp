//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2020 by Dr Jose M. Dominguez
 All rights reserved.

 DualSPHysics is an international collaboration between:
 - EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 - School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
 the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the distribution.
 * Neither the name of the DualSPHysics nor the names of its contributors may be used to endorse or promote products derived 
   from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
 SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/// \file JDsFtForcePoints.cpp \brief Implements the class \ref JDsFtForcePoints.

#include "JDsFtForcePoints.h"
#include "JLog2.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JDataArrays.h"
#include "JVtkLib.h"
#include "JSaveCsv2.h"
#include "JAppInfo.h"
#include "JSphMk.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <cfloat>
#include <algorithm>

using namespace std;

//##############################################################################
//# JDsMooredFloatingLink
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsFtForcePoints::JDsFtForcePoints(bool iscpu,double dp,unsigned ftcount)
  :Log(AppInfo.LogPtr()),Cpu(iscpu),Dp(dp),FtCount(word(ftcount))
{
  ClassName="JDsFtForcePoints";
  FtRadius=NULL;  FtMass=NULL;
  FtCountPt=NULL;  FtBeginPt=NULL;
  PtId=NULL;  PtFtid=NULL;  PartDist=NULL;
  PtPos=NULL;  PtVel=NULL;  PtForce=NULL;
  SelFtIndex=NULL;  SelFtid=NULL;  SelFtCenter=NULL;  
  SelFtAce=NULL;  SelFtOmega=NULL;
  Reset();
  if(FtCount>=USHRT_MAX)Run_Exceptioon("Number of floatings is invalid.");
}

//==============================================================================
/// Destructor
//==============================================================================
JDsFtForcePoints::~JDsFtForcePoints(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsFtForcePoints::Reset(){
  SvCsvPoints=SvVtkPoints=false;
  ConfigPeri(0,false,false,false,TDouble3(0),TDouble3(0),TDouble3(0));
  AllocMemoryFt(0);
  AllocMemorySelFt(0);
  ResizeMemoryPt(0);
  TimeStep=0;
}

//==============================================================================
/// Allocates memory according FtCount.
//==============================================================================
void JDsFtForcePoints::AllocMemoryFt(word ftcount){
  delete[] FtRadius;   FtRadius=NULL;
  delete[] FtMass;     FtMass=NULL;
  delete[] FtCountPt;  FtCountPt=NULL;
  delete[] FtBeginPt;  FtBeginPt=NULL;
  if(ftcount){
    try{
      FtRadius =new float[ftcount];
      FtMass   =new float[ftcount];
      FtCountPt=new word[ftcount];
      FtBeginPt=new word[ftcount];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
    memset(FtRadius ,0,sizeof(float)*ftcount);
    memset(FtMass   ,0,sizeof(float)*ftcount);
    memset(FtCountPt,0,sizeof(word)*ftcount);
    memset(FtBeginPt,0,sizeof(word)*ftcount);
  }
}

//==============================================================================
/// Allocates memory for SelXXX arrays.
//==============================================================================
void JDsFtForcePoints::AllocMemorySelFt(word selftcount){
  delete[] SelFtIndex;  SelFtIndex=NULL;
  delete[] SelFtid;     SelFtid=NULL;
  delete[] SelFtCenter; SelFtCenter=NULL;
  delete[] SelFtAce;    SelFtAce=NULL;
  delete[] SelFtOmega;  SelFtOmega=NULL;
  if(selftcount){
    try{
      SelFtIndex =new word[FtCount];
      SelFtid    =new word[selftcount];
      SelFtCenter=new tdouble3[selftcount];
      SelFtAce   =new tfloat3[selftcount];
      SelFtOmega =new tfloat3[selftcount];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
    memset(SelFtIndex ,0,sizeof(word)*FtCount);
    memset(SelFtid    ,0,sizeof(word)*selftcount);
    memset(SelFtCenter,0,sizeof(tdouble3)*selftcount);
    memset(SelFtAce   ,0,sizeof(tfloat3)*selftcount);
    memset(SelFtOmega ,0,sizeof(tfloat3)*selftcount);
  }
  SelFtCount=selftcount;
}

//==============================================================================
/// Allocates memory according PtCount.
//==============================================================================
void JDsFtForcePoints::ResizeMemoryPt(word ptcount){
  if(!ptcount){
    PtCount=0;
    delete[] PtId;     PtId=NULL;
    delete[] PtFtid;   PtFtid=NULL;
    delete[] PartDist; PartDist=NULL;
    delete[] PtPos;    PtPos=NULL;
    delete[] PtVel;    PtVel=NULL;
    delete[] PtForce;  PtForce=NULL;
  }
  if(ptcount){
    try{
      PtId    =fun::ResizeAlloc(PtId    ,PtCount,ptcount);
      PtFtid  =fun::ResizeAlloc(PtFtid  ,PtCount,ptcount);
      PartDist=fun::ResizeAlloc(PartDist,PtCount,ptcount);
      PtPos   =fun::ResizeAlloc(PtPos   ,PtCount,ptcount);
      PtVel   =fun::ResizeAlloc(PtVel   ,PtCount,ptcount);
      PtForce =fun::ResizeAlloc(PtForce ,PtCount,ptcount);
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
  }
  PtCount=ptcount;
}

//==============================================================================
/// Returns the allocated memory in CPU.
//==============================================================================
llong JDsFtForcePoints::GetAllocMemory()const{  
  llong s=0;
  //Allocated in AllocMemoryFt().
  if(FtRadius) s+=sizeof(float)*FtCount;
  if(FtMass)   s+=sizeof(float)*FtCount;
  if(FtCountPt)s+=sizeof(word)*FtCount;
  if(FtBeginPt)s+=sizeof(word)*FtCount;
  //Allocated in AllocMemorySelFt().
  if(SelFtIndex) s+=sizeof(word)*FtCount;
  if(SelFtid)    s+=sizeof(word)*SelFtCount;
  if(SelFtCenter)s+=sizeof(tdouble3)*SelFtCount;
  if(SelFtAce)   s+=sizeof(tfloat3)*SelFtCount;
  if(SelFtOmega) s+=sizeof(tfloat3)*SelFtCount;
  //Allocated in ResizeMemoryPt().
  if(PtId)    s+=sizeof(word)*PtCount;
  if(PtFtid)  s+=sizeof(word)*PtCount;
  if(PartDist)s+=sizeof(float)*PtCount;
  if(PtPos)   s+=sizeof(tdouble3)*PtCount;
  if(PtVel)   s+=sizeof(tfloat3)*PtCount;
  if(PtForce) s+=sizeof(tfloat3)*PtCount;
  return(s);
}

//==============================================================================
/// Adds point and returns ptid of the new point.
//==============================================================================
word JDsFtForcePoints::AddPoint(unsigned ftid,const tdouble3 &pos){
  if(unsigned(ftid)>=unsigned(FtCount))Run_Exceptioon("Id of floating is invalid.");
  word c=PtCount;
  ResizeMemoryPt(PtCount+1);
  PtId[c]=c;
  PtFtid[c]=ftid;
  PartDist[c]=0;
  PtPos[c]=pos;
  PtVel[c]=TFloat3(0);
  PtForce[c]=TFloat3(0);
  return(c);
}

//==============================================================================
/// Returns idx of requested ptid.
//==============================================================================
word JDsFtForcePoints::GetIdx(word ptid)const{
  word c;
  for(c=0;c<PtCount && PtId[c]!=ptid;c++);
  if(c>=PtCount)Run_Exceptioon("Value ptid is invalid.");
  return(c);
}

//==============================================================================
/// Configures periodic parameters.
//==============================================================================
void JDsFtForcePoints::ConfigPeri(byte periactive,bool perix,bool periy,bool periz
  ,tdouble3 perixinc,tdouble3 periyinc,tdouble3 perizinc)
{
  PeriActive=periactive;
  PeriX=perix;
  PeriY=periy;
  PeriZ=periz;
  PeriXinc=perixinc;
  PeriYinc=periyinc;
  PeriZinc=perizinc;
}

//==============================================================================
/// Configures object for execution.
//==============================================================================
void JDsFtForcePoints::Config(unsigned ftcount,const StFloatingData *ftdata
  ,byte periactive,bool perix,bool periy,bool periz
  ,tdouble3 perixinc,tdouble3 periyinc,tdouble3 perizinc)
{
  if(ftcount!=unsigned(FtCount))Run_Exceptioon("Number of floating is invalid.");
  //-Configures periodic parameters.
  ConfigPeri(periactive,perix,periy,periz,perixinc,periyinc,perizinc);
  //-Stores floating data.
  AllocMemoryFt(FtCount);
  for(word cf=0;cf<FtCount;cf++){
    FtRadius[cf]=ftdata[cf].radius;
    FtMass  [cf]=ftdata[cf].mass;
  }
  const int n=int(PtCount);
  //-Sorts points according ftid and ptid.
  for(int c=0;c<n-1;c++)for(int c2=c+1;c2<n;c2++)if(PtFtid[c]>PtFtid[c2] || (PtFtid[c]==PtFtid[c2] && PtId[c]>PtId[c2])){
    const word      id=PtId  [c]; PtId  [c]=PtId  [c2]; PtId  [c2]=id;
    const word    ftid=PtFtid[c]; PtFtid[c]=PtFtid[c2]; PtFtid[c2]=ftid;
    const tdouble3 pos=PtPos [c]; PtPos [c]=PtPos [c2]; PtPos [c2]=pos;
  }
  //-Prepares PtCountFt[] and PtBeginFt[].
  memset(FtBeginPt,0,sizeof(word)*FtCount);
  for(int c=0;c<n;c++)FtCountPt[PtFtid[c]]++;
  FtBeginPt[0]=0;
  for(word cf=1;cf<FtCount;cf++)FtBeginPt[cf]=FtBeginPt[cf-1]+FtCountPt[cf-1];
  //-Prepares SelFtIndex[] with index to selected floatings.
  unsigned selft=0;
  for(word cf=0;cf<FtCount;cf++)if(FtCountPt[cf])selft++;
  AllocMemorySelFt(selft);
  selft=0;
  for(word cf=0;cf<FtCount;cf++){
    if(FtCountPt[cf]){
      SelFtIndex[cf]=selft;
      SelFtid[selft]=cf;
      selft++;
    }
    else SelFtIndex[cf]=USHRT_MAX;
  }
}

//==============================================================================
/// Checks force points according positions of floating particles.
//==============================================================================
void JDsFtForcePoints::CheckPoints(const JSphMk *mkinfo
  ,unsigned np,const unsigned *idp,const tdouble3 *pos)
{
  for(word c=0;c<PtCount;c++){
    //-Selects floating data.
    const unsigned ftid=PtFtid[c]+mkinfo->GetFirstBlockType(TpPartFloating);
    if(ftid>=mkinfo->Size())Run_Exceptioon("Floating body is missing.");
    const JSphMkBlock* mkb=mkinfo->Mkblock(ftid);
    const unsigned idini=mkb->Begin;
    const unsigned idfin=idini+mkb->Count;
    //-Look for the nearest particle.
    const tdouble3 ptpos=PtPos[c];
    tdouble3 partpos=TDouble3(0);
    double partdist=DBL_MAX;
    for(unsigned p=0;p<np;p++){
      const unsigned id=idp[p];
      if(idini<=id && id<idfin){
        double dis=fgeo::PointsDist(ptpos,pos[p]);
        //Log->Printf("idp[%u]:%u   dis:%g (%g)",p,idp[p],dis,maxdist);
        if(dis<partdist){ partdist=dis; partpos=pos[p]; }
      }
    }
    PartDist[c]=float(partdist);
  }
}

//==============================================================================
/// Shows force points configuration.
//==============================================================================
void JDsFtForcePoints::VisuConfig(std::string txhead,std::string txfoot
  ,unsigned ftcount,const StFloatingData *ftdata)const
{
  if(!txhead.empty())Log->Print(txhead);
  for(word c=0;c<PtCount;c++){
    Log->Printf("  Point_%u(%g,%g,%g) in Floating_%u with mkbound=%u. Distance to particle: %g (%g x Dp)"
      ,PtId[c],PtPos[c].x,PtPos[c].y,PtPos[c].z,PtFtid[c],ftdata[PtFtid[c]].mkbound,PartDist[c],PartDist[c]/Dp);
  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Calculate distance between floating particles & centre according to periodic conditions.
/// Calcula distancia entre pariculas floatin y centro segun condiciones periodicas.
//==============================================================================
tfloat3 JDsFtForcePoints::FtPeriodicDist(const tdouble3 &pos,const tdouble3 &center,float radius)const{
  tdouble3 distd=(pos-center);
  while(PeriX && fabs(distd.x)>radius){
    if(distd.x>0)distd=distd+PeriXinc;
    else distd=distd-PeriXinc;
  }
  while(PeriY && fabs(distd.y)>radius){
    if(distd.y>0)distd=distd+PeriYinc;
    else distd=distd-PeriYinc;
  }
  while(PeriZ && fabs(distd.z)>radius){
    if(distd.z>0)distd=distd+PeriZinc;
    else distd=distd-PeriZinc;
  }
  return(ToTFloat3(distd));
}

//==============================================================================
/// Updates position and velocity of points.
//==============================================================================
void JDsFtForcePoints::UpdatePoints(double timestep,double dt,const StFloatingData *ftdata){
  TimeStep=timestep+dt;
  for(word c=0;c<PtCount;c++){
    const word ftid=PtFtid[c];
    tdouble3 pos=PtPos[c];
    tfloat3 vel=PtVel[c];
    pos.x+=dt*double(vel.x);
    pos.y+=dt*double(vel.y);
    pos.z+=dt*double(vel.z);
    const tfloat3 fvel=ftdata[ftid].fvel;
    const tfloat3 fomega=ftdata[ftid].fomega;
    const tdouble3 fcenter=ftdata[ftid].center;
    SelFtCenter[SelFtIndex[ftid]]=fcenter;
    const tfloat3 dist=(PeriActive? FtPeriodicDist(pos,fcenter,FtRadius[ftid]): ToTFloat3(pos-fcenter)); 
    vel.x=fvel.x+(fomega.y*dist.z-fomega.z*dist.y);
    vel.y=fvel.y+(fomega.z*dist.x-fomega.x*dist.z);
    vel.z=fvel.z+(fomega.x*dist.y-fomega.y*dist.x);
    PtPos[c]=pos; 
    PtVel[c]=vel; 
    PtForce[c]=TFloat3(0);
  }
}

//==============================================================================
/// Computes motion data for floatings (SelFtAce[] and SelFtOmega[]).
//==============================================================================
void JDsFtForcePoints::ComputeFtMotion(){
  for(word cs=0;cs<SelFtCount;cs++){
    const word cf=SelFtid[cs];
    const tdouble3 fcenter=SelFtCenter[cs];
    tdouble3 facetot=TDouble3(0);
    tdouble3 fomegatot=TDouble3(0);
    const word cpini=FtBeginPt[cf];
    const word cpfin=cpini+FtCountPt[cf];
    for(word cp=cpini;cp<cpfin;cp++){
      const tdouble3 face=ToTDouble3(PtForce[cp])/TDouble3(FtMass[cf]); //-Acceleration in this point.
      facetot=facetot+face;
      const tfloat3 dist=(PeriActive? FtPeriodicDist(PtPos[cp],fcenter,FtRadius[cf]): ToTFloat3(PtPos[cp]-fcenter)); 
      fomegatot.x+= face.z*dist.y - face.y*dist.z;
      fomegatot.y+= face.x*dist.z - face.z*dist.x;
      fomegatot.z+= face.y*dist.x - face.x*dist.y;
    }
    SelFtAce[cs]=ToTFloat3(facetot);
    SelFtOmega[cs]=ToTFloat3(fomegatot);
  }
}

//==============================================================================
/// Stores motion data for floatings in ftoforces[].
//==============================================================================
void JDsFtForcePoints::GetFtMotionData(StFtoForces *ftoforces)const{
  for(word cs=0;cs<SelFtCount;cs++){
    const word cf=SelFtid[cs];
    ftoforces[cf].face=SelFtAce[cs];
    ftoforces[cf].fomegaace=SelFtOmega[cs];
  }
}

//==============================================================================
/// Saves VTK with force points.
//==============================================================================
void JDsFtForcePoints::SaveVtkPoints(unsigned numfile)const{
  JDataArrays arrays;
  arrays.AddArray("Pos",PtCount,PtPos,false);
  if(PtFtid) arrays.AddArray("FtId" ,PtCount,PtFtid ,false);
  if(PtVel)  arrays.AddArray("Vel"  ,PtCount,PtVel  ,false);
  if(PtForce)arrays.AddArray("Force",PtCount,PtForce,false);
  const string files=AppInfo.GetDirOut()+"MooringsVtk/FtForcesPoints.vtk";
  Log->AddFileInfo(fun::FileNameSec(files,UINT_MAX),"Saves VTK file with force points (Moordyn coupling).");
  JVtkLib::SaveVtkData(fun::FileNameSec(files,numfile),arrays,"Pos");
}

//==============================================================================
/// Saves CSV with force points.
//==============================================================================
void JDsFtForcePoints::SaveCsvPoints(unsigned numfile)const{
  if(PtCount)Log->AddFileInfo("FtForcesPoints_ft????_pt??.csv","Saves CSV file with force points (Moordyn coupling).");
  for(word cp=0;cp<PtCount;cp++){
    const string file=AppInfo.GetDirOut()+fun::PrintStr("FtForcePoints_ft%04d_pt%02u.csv",PtFtid[cp],cp);
    jcsv::JSaveCsv2 scsv(file,true,AppInfo.GetCsvSepComa());
    //-Saves head.
    scsv.SetHead();
    scsv << "Part;Time [s];PosX [m];PosY [m];PosZ [m];ForceX [N];ForceY [N];ForceZ [N];VelX [m/s];VelY [m/s];VelZ [m/s]" << jcsv::Endl();
    //-Saves data.
    scsv.SetData();
    scsv << jcsv::Fmt(jcsv::TpDouble1,"%g") << jcsv::Fmt(jcsv::TpFloat3,"%g;%g;%g") << jcsv::Fmt(jcsv::TpDouble3,"%g;%g;%g");
    scsv << numfile << TimeStep << PtPos[cp] << PtForce[cp] << PtVel[cp] << jcsv::Endl();
    scsv.SaveData(true);
  }
}

//==============================================================================
/// Saves VTK and CSV files with force points.
//==============================================================================
void JDsFtForcePoints::SaveData(unsigned numfile)const{
  if(SvVtkPoints)SaveVtkPoints(numfile);
  if(SvCsvPoints)SaveCsvPoints(numfile);
}




