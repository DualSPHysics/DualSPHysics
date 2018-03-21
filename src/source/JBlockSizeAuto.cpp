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

/// \file JBlockSizeAuto.cpp \brief Implements the class \ref JBlockSizeAuto.

#include "JBlockSizeAuto.h"
#include "Functions.h"
#include "JLog2.h"
#include "JSaveCsv2.h"
#include <climits>
#include <cfloat>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;


//##############################################################################
//# JBlockSizeAutoSize
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JBlockSizeAutoKer::JBlockSizeAutoKer(JLog2 *log,std::string name,int bsmin,int bsnum,int bsinc,int bsdef):Log(log),Name(name),BsDef(bsdef),BsMin(bsmin),BsInc(bsinc),BsNum(bsnum){
  ClassName="JBlockSizeAutoKer";
  BsActive=NULL;
  Times=NULL;
  OverMean=NULL;
  MeanTot=NULL;  
  MeanExp=NULL;
  AllocateMemory(BsNum);
  InfoData=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JBlockSizeAutoKer::~JBlockSizeAutoKer(){
  DestructorActive=true;
  if(SAVEINFO)SaveFileInfoData();
  Reset();
  AllocateMemory(0);
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JBlockSizeAutoKer::Reset(){
  Nrun=0;
  BsSel=BsDef;
  for(int c=0;c<BsNum;c++)BsActive[c]=true;
  NumActive=BsNum;
  memset(Times,0,sizeof(float)*BsNum);
  memset(OverMean,0,sizeof(float)*BsNum);
  for(int c=0;c<BsNum;c++){
    MeanTot[c].Reset();
    MeanExp[c].InitWeightedExponential(MEANDEPTH,4);
  }
  InfoDataSaved=false;
  InfoDataSizeLine=InfoDataLines=InfoDataCount=0;
  AllocateInfoData(0);
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JBlockSizeAutoKer::AllocateMemory(unsigned size){
  delete[] BsActive;   BsActive=NULL;
  delete[] Times;      Times=NULL;
  delete[] OverMean;   OverMean=NULL;
  delete[] MeanTot;    MeanTot=NULL;
  delete[] MeanExp;    MeanExp=NULL;
  if(size){
    try{
      BsActive  =new bool[size];
      Times     =new float[size];
      OverMean  =new float[size];
      MeanTot   =new JMeanValue[size];
      MeanExp   =new JMeanMoving[size];
    }
    catch(const std::bad_alloc){
      RunException("AllocateMemory","Cannot allocate the requested memory.");
    }
  }
}
 
//==============================================================================
/// Assigns memory for InfoData.
//==============================================================================
void JBlockSizeAutoKer::AllocateInfoData(unsigned nlines){
  delete[] InfoData; InfoData=NULL;
  InfoDataCount=0;
  InfoDataLines=nlines;
  if(InfoDataLines){
    try{
      InfoData=new float[InfoDataSizeLine*InfoDataLines];
    }
    catch(const std::bad_alloc){
      RunException("AllocateInfoData","Cannot allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Stores CSV file and clears InfoData[].
//==============================================================================
void JBlockSizeAutoKer::SaveFileInfoData(){ 
  const char* met="SaveFileInfoData";
  const bool firstsv=!InfoDataSaved;
  const string file=Log->GetDirOut()+Name+".csv";
  if(firstsv)Log->AddFileInfo(file,"Saves information from the automatic block size calculation.");
  jcsv::JSaveCsv2 scsv(file,!firstsv,Log->GetCsvSepComa());
  //-Saves head.
  if(firstsv){
    scsv.SetHead();
    scsv << "Step;Time;Bs";
    for(int ct=0;ct<BsNum;ct++)scsv << fun::PrintStr("Time_%d",BsMin+BsInc*ct);
    for(int ct=0;ct<BsNum;ct++)scsv << fun::PrintStr("Mexp_%d%%",BsMin+BsInc*ct);
    for(int ct=0;ct<BsNum;ct++)scsv << fun::PrintStr("Mtot_%d%%",BsMin+BsInc*ct);
    scsv << jcsv::Endl();
  }
  //-Saves lines in InfoData[].
  scsv.SetData();
  unsigned cpos=0;
  //Log->Printf("%s] SaveInfoData lines:%u",Name.c_str(),InfoDataCount);
  for(unsigned c=0;c<InfoDataCount;c++){
    scsv << fun::PrintStr("%g;%f;%g",InfoData[cpos],InfoData[cpos+1],InfoData[cpos+2]); cpos+=3;
    //Log->Printf("--[%s]--",tx.c_str());
    unsigned nv=BsNum*3;
    for(unsigned cv=0;cv<nv;cv++){ 
      scsv << InfoData[cpos]; cpos++;
    }
    scsv << jcsv::Endl();
  }
  scsv.SaveData();
  InfoDataSaved=true;
  //-Clears buffer.
  InfoDataCount=0;
}

//==============================================================================
/// Stores statistical data in InfoData[].
//==============================================================================
void JBlockSizeAutoKer::SaveInfoData(unsigned nstep,float timestep){ 
  //-Creates buffer for statistical data.
  if(!InfoData){
    //Step;Time;Bs;Time_X;Mexp_X;Mtot_X
    InfoDataSizeLine=unsigned(3+3*BsNum);
    AllocateInfoData(1000);
    //Log->Printf("%s] Allocate",Name.c_str());
  }
  //-Stores CSV file and clears buffer.
  if(InfoDataCount+1>=InfoDataLines)SaveFileInfoData();
  //-Stores new line of data.
  float *dat=InfoData+(InfoDataSizeLine*InfoDataCount);
  unsigned cpos=0;
  dat[cpos]=float(nstep);    cpos++;
  dat[cpos]=float(timestep); cpos++;
  dat[cpos]=float(BsSel);    cpos++;
  for(int ct=0;ct<BsNum;ct++){
    dat[cpos]=Times[ct];     cpos++;
  }
  for(int ct=0;ct<BsNum;ct++){
    dat[cpos]=OverMean[ct];  cpos++;
  }
  for(int ct=0;ct<BsNum;ct++){
    dat[cpos]=float(MeanTot[ct].GetMean());  cpos++;
  }
  if(cpos!=InfoDataSizeLine)RunException("SaveInfoData","Error in InfoDataSizeLine...");
  InfoDataCount++;
}

//==============================================================================
/// Processes calculated times.
//==============================================================================
void JBlockSizeAutoKer::ProcessTimes(double timestep,unsigned nstep){
  //-Computes minimum time.
  int ctmin=0;
  float tmin=FLT_MAX;
  for(int ct=0;ct<BsNum;ct++)if(BsActive[ct]){
    const float t=Times[ct];
    if(t==FLT_MAX){ 
      //Log->Printf("-------------->> bs %u DESCARTADO.",GetBs(ct));
      BsActive[ct]=false;
      MeanTot[ct].Reset();
      MeanExp[ct].Reset();
      NumActive--;
    }
    if(tmin>t){ 
      tmin=t;
      ctmin=ct;
    }
  }
  /// Adds values of overhead to compute mean values and computes minimum overmean.
  int ctmin2=0;
  float tmin2=FLT_MAX;
  for(int ct=0;ct<BsNum;ct++){
    if(BsActive[ct]){
      //printf("------> Times[%u]:%f\n",ct,Times[ct]);
      const double overhead=Times[ct]/tmin-1;
      MeanTot[ct].AddValue(overhead);
      //MeanSimple[ct].AddValue(overhead);
      MeanExp[ct].AddValue(overhead);
      OverMean[ct]=float(MeanExp[ct].GetWeightedMean());
      if(tmin2>OverMean[ct]){
        tmin2=OverMean[ct];
        ctmin2=ct;
      }
    }
    else Times[ct]=OverMean[ct]=0;
  }
  //-Selects optimum Blocksize.
  BsSel=GetBs(ctmin2);
  //-Stores statistical data.
  if(SAVEINFO)SaveInfoData(nstep,float(timestep));
  //-Discards the lowest sizes of block.
  if(Nrun>=REMOVESTART && NumActive>REMOVELIMIT){
    int nremove=int(float(NumActive)*REMOVEPRC/100);
    if(!nremove)nremove=1;
    for(int cr=0;cr<nremove;cr++){
      int ctmax=0;
      float tmax=-FLT_MAX;
      for(int ct=0;ct<BsNum;ct++)if(BsActive[ct] && tmax<float(MeanTot[ct].GetMean())){
        tmax=float(MeanTot[ct].GetMean());
        ctmax=ct;
      }
      BsActive[ctmax]=false;
      MeanTot[ctmax].Reset();
      MeanExp[ctmax].Reset();
      NumActive--;
    }
  }
  Nrun++;
}


//##############################################################################
//# JBlockSizeAuto
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JBlockSizeAuto::JBlockSizeAuto(JLog2 *log,unsigned steps):Log(log),StepsInterval(steps){
  ClassName="JBlockSizeAuto";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JBlockSizeAuto::~JBlockSizeAuto(){
  DestructorActive=true;
  Reset();
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JBlockSizeAuto::Reset(){
  for(unsigned c=0;c<GetKernelsCount();c++)delete Kernels[c];
  Kernels.clear();
}
 
//==============================================================================
/// Returns index of kernel by name (UINT_MAX when it does not exist).
//==============================================================================
unsigned JBlockSizeAuto::GetKernelByName(std::string name)const{ 
  unsigned ret=UINT_MAX;
  for(unsigned c=0;c<GetKernelsCount() && ret==UINT_MAX;c++)if(Kernels[c]->Name==name)ret=c;
  return(ret);
}
 
//==============================================================================
/// Adds new kernel.
//==============================================================================
void JBlockSizeAuto::AddKernel(std::string name,int bsmin,int bsnum,int bsinc,int bsdefault){ 
  if(GetKernelByName(name)!=UINT_MAX)RunException("AddKernel","The name is already in use.");
  JBlockSizeAutoKer *ker=new JBlockSizeAutoKer(Log,name,bsmin,bsnum,bsinc,bsdefault);
  Kernels.push_back(ker);
}
 
//==============================================================================
/// Processes calculated times.
//==============================================================================
void JBlockSizeAuto::ProcessTimes(double timestep,unsigned nstep){ 
  for(unsigned c=0;c<GetKernelsCount();c++)Kernels[c]->ProcessTimes(timestep,nstep);
}


