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

/// \file JDsPartFloatSave.cpp \brief Implements the class \ref JDsPartFloatSave.

#include "JDsPartFloatSave.h"
#include "JPartFloatInfoBi4.h"
#include "JDsFtForcePoints.h"
#include "Functions.h"
#include "JAppInfo.h"
#include "JLog2.h"

#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JDsPartFloatSave
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsPartFloatSave::JDsPartFloatSave(bool cpu,std::string appname,std::string dir
  ,word mkboundfirst,unsigned ftcount,double timeout,double timeout2)
  :Log(AppInfo.LogPtr()),Cpu(cpu),AppName(appname),Dir(dir)
  ,MkBoundFirst(mkboundfirst),FtCount(ftcount)
  ,TimeOut(timeout),TimeOut2(timeout2)
{
  ClassName="JDsPartFloatSave";
  FtData=NULL;
  FtFile1=NULL;
  FtFile2=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsPartFloatSave::~JDsPartFloatSave(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsPartFloatSave::Reset(){
  FptCount=0;
  delete FtData;   FtData=NULL;
  delete FtFile1;  FtFile1=NULL;
  delete FtFile2;  FtFile2=NULL;
}

//==============================================================================
/// Initial configuration of floating data.
//==============================================================================
void JDsPartFloatSave::ConfigFtData(unsigned ftcount
  ,const StFloatingData* ftobjs)
{
  Reset();
  if(!FtCount)Run_Exceptioon("No floating bodies available.");
  if(FtCount!=ftcount)Run_Exceptioon("Number of floating bodies is invalid.");
  //-Loads constant data.
  FtData=new JPartFloatInfoBi4Data(FtCount,FptCount,MkBoundFirst);
  for(unsigned cf=0;cf<FtCount;cf++){
    const StFloatingData& v=ftobjs[cf];
    FtData->SetCteData(cf,v.mkbound,v.begin,v.count,v.mass,v.massp,v.radius);
  }
  //-Creates main JPartFloatBi4Save object for output data.
  FtFile1=new JPartFloatInfoBi4Save(AppName,Dir);
  FtFile1->Config(true,TimeOut,FtData);
  //-Creates extra JPartFloatBi4Save object for output data.
  if(TimeOut2>=0){
    FtFile2=new JPartFloatInfoBi4Save(AppName,Dir);
    FtFile2->Config(false,TimeOut2,FtData);
  }
  FtData->ResizeData(1);
}

//==============================================================================
/// Reconfigures PartFloatSave to include information from ForcePoints.
//==============================================================================
void JDsPartFloatSave::ReconfigForcePoints(const JDsFtForcePoints* forcepoints){
  //Log->Printf("\nAAA_000  fp:%p",forcepoints);
  if(FptCount==0 && forcepoints && forcepoints->GetPtCount()){
    FptCount=forcepoints->GetPtCount();
    //-Reconfigures FtData.
    if(FtData){
      JPartFloatInfoBi4Data* ftdata0=FtData;
      if(ftdata0->GetPartCount()!=1)Run_Exceptioon("Reconfiguration is invalid.");
      FtData=new JPartFloatInfoBi4Data(FtCount,FptCount,MkBoundFirst);
      FtData->LoadCteData(ftdata0);
      FtData->ResizeData(1);
      //-Copy data from ftdata0 and add forcepoints data.
      FtData->SetPartData0(ftdata0);
      FtData->SetForcePoints0(FptCount,forcepoints->GetPtMkBound()
        ,forcepoints->GetPtPos(),forcepoints->GetPtForce());
      FtData->SetTimeData0(ftdata0->GetPartCpart(0),ftdata0->GetPartTimeStep(0),ftdata0->GetPartStep(0));
      delete ftdata0; ftdata0=NULL;
    }
    //-Reconfigures FtFile1 and FtFile2.
    delete FtFile1; FtFile1=NULL;
    delete FtFile2; FtFile2=NULL;

    //-Creates main JPartFloatBi4Save object for output data.
    FtFile1=new JPartFloatInfoBi4Save(AppName,Dir);
    FtFile1->Config(true,TimeOut,FtData);
    //-Creates extra JPartFloatBi4Save object for output data.
    if(TimeOut2>=0){
      FtFile2=new JPartFloatInfoBi4Save(AppName,Dir);
      FtFile2->Config(false,TimeOut2,FtData);
    }

    SaveInitial();
    if(ExtraIsActive()){
      AddDataExtra(FtData->GetPartCpart(0),FtData->GetPartTimeStep(0),FtData->GetPartStep(0));
    }
  }
}

//==============================================================================
/// Initial files recording with info from Data.
//==============================================================================
void JDsPartFloatSave::SaveInitial(){
  if(FtFile1)FtFile1->SaveInitial();
  if(FtFile2)FtFile2->SaveInitial();
}

//==============================================================================
/// Returns file name BI4.
//==============================================================================
std::string JDsPartFloatSave::GetFileName(bool mainfile)const{
  string fname;
  if(mainfile  && FtFile1)fname=JPartFloatInfoBi4Save::GetFileNameDef(true);
  if(!mainfile && FtFile2)fname=JPartFloatInfoBi4Save::GetFileNameDef(false);
  return(fname);
}

//==============================================================================
/// Saves floating data from simulation.
//==============================================================================
void JDsPartFloatSave::SetFtData(int cpart,double timestep,unsigned nstep
  ,const StFloatingData* ftobjs,const JDsFtForcePoints* forcepoints)
{
  //Log->Printf("**JDsPartFloatSave::SetFtData>> nstep:%u  t:%g  ftn:%u fcn:%u\n"
  //  ,nstep,timestep,FtCount,(forcepoints? forcepoints->GetPtCount(): 0));
  //-Stores data from ftobjs.
  for(unsigned cf=0;cf<FtCount;cf++){
    const StFloatingData& v=ftobjs[cf];
    FtData->SetPartData0(cf,v.center,v.fvel,v.fomega,v.facelin,v.faceang
     ,v.extforcelin,v.extforceang,v.fluforcelin,v.fluforceang,v.preacelin,v.preaceang);
  }
  //-Stores data from force points.
  const unsigned fptcount=(forcepoints? forcepoints->GetPtCount(): 0);
  if(FptCount!=fptcount)Run_Exceptioon("Number of force points does not match.");
  if(FptCount)FtData->SetForcePoints0(FptCount,forcepoints->GetPtMkBound()
    ,forcepoints->GetPtPos(),forcepoints->GetPtForce());
  //-Stores time data.
  FtData->SetTimeData0(cpart,timestep,nstep);
}

//==============================================================================
/// Saves floating data in file bi4 for main output.
//==============================================================================
void JDsPartFloatSave::SaveDataMain(int cpart,double timestep,unsigned nstep){
  //Log->Printf("**JDsPartFloatSave::SaveData_Main>> nstep:%u  t:%g\n",nstep,timestep);
  FtFile1->SavePart(cpart,nstep,FtData);
}

//==============================================================================
/// Add floating data for extra output.
//==============================================================================
void JDsPartFloatSave::AddDataExtra(int cpart,double timestep,unsigned nstep){
  //Log->Printf("**JDsPartFloatSave::SaveData_Extra>> nstep:%u  t:%g\n",nstep,timestep);
  FtFile2->AddDataPart(cpart,nstep,FtData);
}

//==============================================================================
/// Saves stored EXTRA floating data in BI4 file.
//==============================================================================
void JDsPartFloatSave::SaveDataExtra(){
  if(FtFile2)FtFile2->SaveStoredData();
}

