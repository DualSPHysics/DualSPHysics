//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2025 by Dr Jose M. Dominguez
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

/// \file JDsMooredFloatings.cpp \brief Implements the class \ref JDsMooredFloatings.

#include "JDsMooredFloatings.h"
#include "JXml.h"
#include "JLog2.h"
#include "JSpVtkShape.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JSaveCsv2.h"
#include "JDsFtForcePoints.h"
#include "JAppInfo.h"
#include "DSphMoorDynPlus.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JDsMooredFloating
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsMooredFloating::JDsMooredFloating(word mkbound)
  :Log(AppInfo.LogPtr()),FloatingMkBound(mkbound)
{
  ClassName="JDsMooredFloating";
  Reset();
}

//==============================================================================
/// Destructor
//==============================================================================
JDsMooredFloating::~JDsMooredFloating(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsMooredFloating::Reset(){
  FtIdx=FtId=UINT_MAX;
  Fairleads.clear();
}

//==============================================================================
/// Configura mooring con datos de floating asociado.
//==============================================================================
void JDsMooredFloating::ConfigIds(unsigned ftidx,unsigned ftid){
  if(FtIdx!=UINT_MAX || FtId!=UINT_MAX)Run_Exceptioon("Floating has already configured.");
  FtIdx=ftidx;
  FtId=ftid;
}

//==============================================================================
/// Adds new link point.
//==============================================================================
void JDsMooredFloating::AddFairlead(unsigned fairnum,const tdouble3& linkpos
  ,word ptid)
{
  Fairleads.push_back(StrLinkData(FtId,fairnum,linkpos,ptid));
}

//==============================================================================
/// Shows mooring configuration.
//==============================================================================
void JDsMooredFloating::VisuConfig()const{
  Log->Printf("  Floating_%u (MkBound:%u)",FtId,FloatingMkBound);
  for(unsigned c=0;c<Count();c++){
    const tdouble3 pos=Fairleads[c].linkpos;
    Log->Printf("    Link_%u:  LinkPos:(%g,%g,%g)  Point_%u",Fairleads[c].fairnum,pos.x,pos.y,pos.z,Fairleads[c].ptid);
  }
}

//==============================================================================
/// Returns the requested fairlead data.
//==============================================================================
JDsMooredFloating::StLinkData JDsMooredFloating::GetFairlead(unsigned fairnum)const{
  if(fairnum>=Count())Run_Exceptioon("The requested fairlead is invalid.");
  return(Fairleads[fairnum]);
}


//##############################################################################
//# JDsMooredFloatings
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsMooredFloatings::JDsMooredFloatings(std::string dircase,std::string casename
  ,tfloat3 gravity,double timemax,double dtout)
  :Log(AppInfo.LogPtr()),DirCase(dircase),CaseName(casename),Gravity(gravity)
  ,TimeMax(timemax),DtOut(dtout)
{
  ClassName="JDsMooredFloatings";
  MoorDynPlusReady=false;
  FairArrays=false;
  FairNftm=0;
  FairFtmNum=NULL;
  FairleadPos=NULL;
  FairleadVel=NULL;
  FairleadForce=NULL;
  MoorDynPlus_LogInit(Log);
  Reset();
}
//==============================================================================
/// Destructor
//==============================================================================
JDsMooredFloatings::~JDsMooredFloatings(){
  DestructorActive=true;
  Reset();
}
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsMooredFloatings::Reset(){
  FileLines="";
  MdpDir="";
  StartTime=StartRamp=StartEnd=0;
  SvVtkLines=SvVtkMoorings=SvCsvPoints=SvVtkPoints=false;
  for(unsigned c=0;c<Count();c++)delete Floatings[c];
  Floatings.clear();
  if(MoorDynPlusReady){
    if(MoorDynPlus_LinesClose())Run_Exceptioon("Error releasing moorings in MoorDynPlus library.");
  }
  MoorDynPlusReady=false;
  FreeFairMemory();  //-Frees link data arrays.
  //SaveDataTime=NextTime=0;
  //LastTimeOk=-1;
}

//==============================================================================
/// Returns idx of floating with requested mkbound.
//==============================================================================
unsigned JDsMooredFloatings::GetFloatingByMkbound(word mkbound)const{
  unsigned c=0;
  for(c=0;c<Count() && Floatings[c]->FloatingMkBound!=mkbound;c++);
  return(c>=Count()? UINT_MAX: c);
}

//==============================================================================
/// Reads list of mooredfloatings in the XML node.
//==============================================================================
void JDsMooredFloatings::ReadXml(const JXml* sxml,TiXmlElement* lis){
  sxml->CheckElementNames(lis,true,"moordynplus moordyn start savevtk_lines savevtk_moorings savecsv_points savevtk_points mooredfloatings");
  //-Loads configuration file for MoorDynPlus solver.
  if(sxml->CheckElementActive(lis,"moordyn")) Run_Exceptioon(fun::PrintStr("Option <moordyn> is deprecated, use <moordynplus> instead."));
  if(sxml->CheckElementActive(lis,"moordynplus")){
    FileLines=sxml->ReadElementStr(lis,"moordynplus","file",true);
    if(FileLines.empty()){
      FileLines=fun::StrReplace(FileLines,"[CaseName]",CaseName);
      MdpDir=AppInfo.GetDirOut()+"moordynplus_data";
      fun::Mkdir(MdpDir);
      MdpDir=MdpDir+"/";
      //IME//if(fun::CpyFile(DirData+FileLines,MdpDir+"lines.txt"))Run_ExceptioonFile("Error: File could not be created.",MdpDir+"lines.txt");
    }
  }
  //-Loads initial configuration.
  StartTime=sxml->ReadElementDouble(lis,"start","value",true,0);
  StartRamp=sxml->ReadElementDouble(lis,"start","ramptime",true,0);
  StartEnd=(StartRamp>0? StartTime+StartRamp: StartTime);
  //-Loads configuration to save VTK and CSV files.
  SvVtkLines=sxml->ReadElementBool(lis,"savevtk_lines","value",true,true);
  SvVtkMoorings=sxml->ReadElementBool(lis,"savevtk_moorings","value",true,true);
  SvCsvPoints=sxml->ReadElementBool(lis,"savecsv_points","value",true,true);
  SvVtkPoints=sxml->ReadElementBool(lis,"savevtk_points","value",true,false);
  //-Loads floatings with moorings.
  TiXmlElement* mlis=lis->FirstChildElement("mooredfloatings");
  if(mlis){
    TiXmlElement* ele=mlis->FirstChildElement("floating"); 
    while(ele){
      if(sxml->CheckElementActive(ele)){
        const word floatingmkbound=sxml->GetAttributeWord(ele,"mkbound");
        if(GetFloatingByMkbound(floatingmkbound)!=UINT_MAX)Run_Exceptioon(fun::PrintStr("Floating mkbound=%d is already configured.",floatingmkbound));
        JDsMooredFloating* mo=new JDsMooredFloating(floatingmkbound);
        Floatings.push_back(mo);
      }
      ele=ele->NextSiblingElement("floating");
    }
  }
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JDsMooredFloatings::LoadXml(const JXml* sxml,const std::string& place){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Asocia moorings con floatings y los reordena.
//==============================================================================
void JDsMooredFloatings::ConfigFloatings(unsigned ftcount
  ,const StFloatingData* ftdata)
{
  //-Ordena floatings amarrados segun MkBound de floating.
  const unsigned nftm=Count();
  for(unsigned c=0;c<nftm-1;c++)for(unsigned c2=c+1;c2<nftm;c2++)if(Floatings[c]->FloatingMkBound>Floatings[c2]->FloatingMkBound){
    JDsMooredFloating* mo=Floatings[c]; Floatings[c]=Floatings[c2]; Floatings[c2]=mo;
  }
  //-Asigna id general de floating y numero de floating amarrado.
  for(unsigned cfm=0;cfm<nftm;cfm++){
    JDsMooredFloating* mo=Floatings[cfm];
    unsigned cf;
    for(cf=0;cf<ftcount && ftdata[cf].mkbound!=mo->FloatingMkBound;cf++);
    if(cf>=ftcount)Run_Exceptioon("Floating body used in mooring not found.");
    mo->ConfigIds(cfm,cf);
  }
}

//==============================================================================
/// Configures object.
//==============================================================================
void JDsMooredFloatings::Config(unsigned ftcount,const StFloatingData* ftdata
  ,JDsFtForcePoints* forcepoints)
{
  //-Asocia moorings con floatings y los reordena.
  ConfigFloatings(ftcount,ftdata);
  const unsigned nftm=Count();
  //-Checks errors.
  if(nftm<1)Run_Exceptioon("There are not moored floatings.");
  //if(nftm>1)Run_Exceptioon("MoorDynPlus only supports one moored floating.");
  
  //-Initilizes MoorDynPlus moorings.
  {
    //-Prepares data to initialize moorings.
    unsigned* ftmkb=new unsigned[nftm];
    tdouble3* ftvellin=new tdouble3[nftm];
    tdouble3* ftvelang=new tdouble3[nftm];
    for(unsigned cfm=0;cfm<nftm;cfm++){
      const unsigned ftid=Floatings[cfm]->GetFtId();
      ftmkb   [cfm]=Floatings[cfm]->FloatingMkBound;
      ftvellin[cfm]=ToTDouble3(ftdata[ftid].fvel);
      ftvelang[cfm]=ToTDouble3(ftdata[ftid].fomega);
    }
    //-Initilizes MoorDynPlus.
    string filexml=FileLines;
    string nodexml="moordynplus";
    if(filexml.empty()){
      filexml=DirCase+CaseName+".xml";
      nodexml="case.execution.special.moorings.moordynplus";
    }
    if(Gravity.z==0)Run_Exceptioon("Gravity.z equal to zero is not allowed.");
    if(Gravity.x || Gravity.y)Log->PrintfWarning("Gravity.x or Gravity.y are not zero but only gravity.z=%f is used for MoorDynPlus.",fabs(Gravity.z));
    if(MoorDynPlus_LinesInit(filexml,nodexml,MdpDir,nftm,ftmkb,ftvellin,ftvelang,Gravity,TimeMax,DtOut))
      Run_Exceptioon("Error initializing moorings in MoorDynPlus library.");
    MoorDynPlusReady=true;
    //-Free memory.
    delete[] ftmkb;    ftmkb   =NULL;
    delete[] ftvellin; ftvellin=NULL;
    delete[] ftvelang; ftvelang=NULL;
  }

  //-Adds link positions for each floating.
  for(unsigned cfm=0;cfm<nftm;cfm++){
    const unsigned ftid=Floatings[cfm]->GetFtId();
    const word     ftmkb=Floatings[cfm]->FloatingMkBound;
    const unsigned nfairleads=MoorDynPlus_FairsCount(cfm);
    for(unsigned ck=0;ck<nfairleads;ck++){
      const unsigned nodes=MoorDynPlus_SegsCount(cfm,ck);
      const tdouble3 pos=MoorDynPlus_GetNodePosLink(cfm,ck);
      const word ptid=forcepoints->AddPoint(ftid,ftmkb,pos);
      Floatings[cfm]->AddFairlead(ck,pos,ptid);
    }
  }

  //-Configures savedata in JDsFtForcePoints object.
  forcepoints->SetSaveData(SvCsvPoints,SvVtkPoints);
}

//==============================================================================
/// Shows moorings configuration.
//==============================================================================
void JDsMooredFloatings::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  Log->Printf("  Start...........: %g [s]  (ramp time: %g [s])",StartTime,StartRamp);
  Log->Printf("  Output directory: %s",MdpDir.c_str());
  for(unsigned cfm=0;cfm<Count();cfm++)Floatings[cfm]->VisuConfig();
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Saves VTK with moorings full or only lines.
//==============================================================================
void JDsMooredFloatings::SaveVtkMoorings(unsigned numfile,bool svlines)const{
  JSpVtkShape ss;
  unsigned vpossize=0;
  tfloat3* vpos=NULL;
  const unsigned nlines=MoorDynPlus_LinesCount();
  for(unsigned cl=0;cl<nlines;cl++){
    const unsigned nodes=MoorDynPlus_SegsCount(cl)+1;
    if(nodes>vpossize){
      vpossize=nodes;
      vpos=fun::ResizeAlloc(vpos,0,vpossize);
    }
    for(unsigned cn=0;cn<nodes;cn++)vpos[cn]=ToTFloat3(MoorDynPlus_GetNodePos(cl,cn));
    ss.AddLines(nodes,vpos,word(cl));
    if(svlines)ss.AddPoints(nodes,vpos,word(cl));
    const float dist=(nodes>=2? fgeo::PointsDist(vpos[0],vpos[1])/10: 1);
    if(!svlines)for(unsigned cn=0;cn<nodes;cn++)ss.AddSphere(vpos[cn],dist,word(cl));
  }
  delete[] vpos; vpos=NULL;
  //-Genera fichero VTK.
  if(svlines){
    Log->AddFileInfo(AppInfo.GetDirOut()+"MooringsVtk/MooringsLines_????.vtk","Saves VTK file with moorings.");
    const string file=AppInfo.GetDirOut()+fun::FileNameSec("MooringsVtk/MooringsLines.vtk",numfile);
    ss.SaveVtk(file,"Line");
  }
  else{
    Log->AddFileInfo(AppInfo.GetDirOut()+"MooringsVtk/Moorings_????.vtk","Saves VTK file with moorings.");
    const string file=AppInfo.GetDirOut()+fun::FileNameSec("MooringsVtk/Moorings.vtk",numfile);
    ss.SaveVtk(file,"Line");
  }
}

//==============================================================================
/// Saves data of moorings.
//==============================================================================
void JDsMooredFloatings::SaveData(unsigned numfile)const{
  if(SvVtkLines)   SaveVtkMoorings(numfile,true); 
  if(SvVtkMoorings)SaveVtkMoorings(numfile,false); 
}

//==============================================================================
/// Allocates memory for fairlead link data.
//==============================================================================
void JDsMooredFloatings::AllocFairMemory(){
  if(FairArrays)FreeFairMemory();
  FairNftm=Count();
  if(FairNftm){
    FairFtmNum   =new unsigned[FairNftm];
    FairleadPos  =new double**[FairNftm];
    FairleadVel  =new double**[FairNftm];
    FairleadForce=new double**[FairNftm];
    for(unsigned cf=0;cf<FairNftm;cf++){
      const JDsMooredFloating* mo=Floatings[cf];
      const unsigned nfairs=mo->Count();
      FairFtmNum   [cf]=nfairs;
      if(nfairs){
        FairleadPos  [cf]=new double*[nfairs];
        FairleadVel  [cf]=new double*[nfairs];
        FairleadForce[cf]=new double*[nfairs];
        for(unsigned cfa=0;cfa<nfairs;cfa++){
          FairleadPos  [cf][cfa]=new double[3];
          FairleadVel  [cf][cfa]=new double[3];
          FairleadForce[cf][cfa]=new double[3];
        }
      }
    }
  }
  FairArrays=true;
}

//==============================================================================
/// Frees memory for fairlead link data.
//==============================================================================
void JDsMooredFloatings::FreeFairMemory(){
  for(unsigned cf=0;cf<FairNftm;cf++){
    const unsigned nfairs=FairFtmNum[cf];
    for(unsigned cfa=0;cfa<nfairs;cfa++){
      delete[] FairleadPos  [cf][cfa];  FairleadPos  [cf][cfa]=NULL;
      delete[] FairleadVel  [cf][cfa];  FairleadVel  [cf][cfa]=NULL;
      delete[] FairleadForce[cf][cfa];  FairleadForce[cf][cfa]=NULL;
    }
    delete[] FairleadPos  [cf];  FairleadPos  [cf]=NULL;
    delete[] FairleadVel  [cf];  FairleadVel  [cf]=NULL;
    delete[] FairleadForce[cf];  FairleadForce[cf]=NULL;
  }
  delete[] FairFtmNum;     FairFtmNum   =NULL;
  delete[] FairleadPos;    FairleadPos  =NULL;
  delete[] FairleadVel;    FairleadVel  =NULL;
  delete[] FairleadForce;  FairleadForce=NULL;
  FairArrays=false;
  FairNftm=0;
}

//==============================================================================
/// Calcula acelaracion y momento angular para aplicar al floating.
//==============================================================================
void JDsMooredFloatings::ComputeForces(unsigned nstep,double timestep,double dt
  ,JDsFtForcePoints* forcepoints)
{
  //Log->Printf("Moorings> timestep:%f",timestep);
  //-Allocates memory for fairlead link data.
  if(!FairArrays)AllocFairMemory();

  //-Loads position and velocity data for MoorDynPlus calculation.
  for(unsigned cf=0;cf<FairNftm;cf++){
    const JDsMooredFloating* mo=Floatings[cf];
    const unsigned nfairs=FairFtmNum[cf];
    for(unsigned cfa=0;cfa<nfairs;cfa++){
      JDsMooredFloating::StLinkData fairlead=mo->GetFairlead(cfa);
      const word     ptid=fairlead.ptid;
      const tdouble3 pos=forcepoints->GetPos(ptid);
      const tfloat3  vel=forcepoints->GetVel(ptid);
      FairleadPos[cf][cfa][0]=pos.x;  FairleadPos[cf][cfa][1]=pos.y;  FairleadPos[cf][cfa][2]=pos.z;
      FairleadVel[cf][cfa][0]=vel.x;  FairleadVel[cf][cfa][1]=vel.y;  FairleadVel[cf][cfa][2]=vel.z;
    }
  }

  //-Computes forces on lines.
  if(MoorDynPlus_FairleadsCalc(FairNftm,FairleadPos,FairleadVel,FairleadForce,timestep,dt))
    Run_Exceptioon("Error calculating forces by MoorDynPlus.");
  
  //-Computes ramp factor according start configuration.
  double framp=1.;
  if(timestep<StartEnd){
    framp=(timestep<StartTime || StartRamp==0? 0.: 1.-(StartEnd-timestep)/StartRamp);
    framp=framp*framp;
  }

  //-Updates forces in JDsFtForcePoints object.
  for(unsigned cf=0;cf<FairNftm;cf++){
    const JDsMooredFloating* mo=Floatings[cf];
    const unsigned nfairs=FairFtmNum[cf];
    for(unsigned cfa=0;cfa<nfairs;cfa++){
      const word ptid=mo->GetFairlead(cfa).ptid;
      tfloat3 force=ToTFloat3(TDouble3(FairleadForce[cf][cfa][0],FairleadForce[cf][cfa][1],FairleadForce[cf][cfa][2]));
      if(timestep<StartEnd)force=force*float(framp);
      forcepoints->SetForce(ptid,force);
    }
  }
}

