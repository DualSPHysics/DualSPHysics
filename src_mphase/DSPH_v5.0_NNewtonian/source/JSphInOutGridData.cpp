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

/// \file JSphInOutGrid.cpp \brief Implements the class \ref JSphInOutGrid.

#include "JSphInOutGridData.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JReadDatafile.h"
#include "JSaveCsv2.h"
#include "JDataArrays.h"
#include "JVtkLib.h"
#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JSphGpu_InOut_iker.h"
#endif

#include <cfloat>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSphInOutGridDataTime
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutGridDataTime::JSphInOutGridDataTime(unsigned nx,unsigned nz)
 :Nx(nx),Nz(nz),Npt(Nx*Nz)
{
  ResetInit();
}

//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutGridDataTime::JSphInOutGridDataTime(unsigned nx,unsigned nz,double time,const float *velx,const float *velz)
 :Nx(nx),Nz(nz),Npt(Nx*Nz)
{
  ResetInit();
  SetData(time,velx,velz);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOutGridDataTime::~JSphInOutGridDataTime(){
  delete[] Velx; Velx=NULL;
  delete[] Velz; Velz=NULL;
}

//==============================================================================
/// Initialisation of variables first time.
//==============================================================================
void JSphInOutGridDataTime::ResetInit(){
  Time=-1;
  Velx=NULL;
  Velz=NULL;
}

//==============================================================================
/// Allocate memory for data.
//==============================================================================
void JSphInOutGridDataTime::AllocData(bool usevelz){
  if(!Velx)Velx=new float[Npt];
  if(usevelz && !Velz)Velz=new float[Npt];
}

//==============================================================================
/// Set data.
//==============================================================================
void JSphInOutGridDataTime::SetData(double time,const float *velx,const float *velz){
  AllocData(velz!=NULL);
  Time=time;
  memcpy(Velx,velx,sizeof(float)*Npt);
  if(velz)memcpy(Velz,velz,sizeof(float)*Npt);
}

//==============================================================================
/// Interpolate data between gdt and gdt2.
//==============================================================================
void JSphInOutGridDataTime::Interpolate(double time,const JSphInOutGridDataTime *gdt,const JSphInOutGridDataTime *gdt2){
  AllocData(gdt->Velz!=NULL);
  Time=time;
  const double fx=((time-gdt->Time)/(gdt2->Time-gdt->Time));
  const float* v0=gdt->GetVelx();
  const float* v1=gdt2->GetVelx();
  for(unsigned p=0;p<Npt;p++)Velx[p]=float(fx*(v1[p]-v0[p])+v0[p]);
  if(Velz){
    v0=gdt->GetVelz();
    v1=gdt2->GetVelz();
    for(unsigned p=0;p<Npt;p++)Velz[p]=float(fx*(v1[p]-v0[p])+v0[p]);
  }
}


//##############################################################################
//# JSphInOutGrid
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutGridData::JSphInOutGridData():Log(AppInfo.LogPtr()){
  ClassName="JSphInOutGridData";
  SelData=NULL;
  #ifdef _WITHGPU
    Velx0g=Velx1g=SelVelxg=NULL;
    Velz0g=Velz1g=SelVelzg=NULL;
  #endif
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOutGridData::~JSphInOutGridData(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInOutGridData::Reset(){
  File="";
  const unsigned nt=CountTimes();
  for(unsigned ct=0;ct<nt;ct++){ delete DataTimes[ct]; DataTimes[ct]=NULL; }
  DataTimes.clear();
  Nx=Nz=Npt=0;
  Dpx=Dpz=0;
  UseVelz=false;
  PosMin=PosMax=TDouble3(0);
  SelCt=0;
  delete SelData; SelData=NULL;
  #ifdef _WITHGPU
    TimeGpu=-1;
    CtVel0=CtVel1=CtSelVel=UINT_MAX;
    FreeMemoryGpu();
  #endif
}

#ifdef _WITHGPU
//==============================================================================
/// Allocate GPU memory for data.
//==============================================================================
void JSphInOutGridData::AllocateMemoryGpu(){
  if(Npt){
    fcuda::Malloc(&Velx0g,Npt);
    fcuda::Malloc(&Velx1g,Npt);
    fcuda::Malloc(&SelVelxg,Npt);
    if(UseVelz){
      fcuda::Malloc(&Velz0g,Npt);
      fcuda::Malloc(&Velz1g,Npt);
      fcuda::Malloc(&SelVelzg,Npt);
    }
  }
}

//==============================================================================
/// Free GPU memory for data.
//==============================================================================
void JSphInOutGridData::FreeMemoryGpu(){
  if(Velx0g)  cudaFree(Velx0g);    Velx0g  =NULL;
  if(Velx1g)  cudaFree(Velx1g);    Velx1g  =NULL;
  if(SelVelxg)cudaFree(SelVelxg);  SelVelxg=NULL;
  if(Velz0g)  cudaFree(Velz0g);    Velz0g  =NULL;
  if(Velz1g)  cudaFree(Velz1g);    Velz1g  =NULL;
  if(SelVelzg)cudaFree(SelVelzg);  SelVelzg=NULL;
}
#endif

//==============================================================================
/// Sets origin of grid.
//==============================================================================
void JSphInOutGridData::SetPosMin(const tdouble3 &posmin){
  PosMin=posmin;
  PosMax.x=PosMin.x+Dpx*(Nx-1); 
  PosMax.y=PosMin.y;
  PosMax.z=PosMin.z+Dpz*(Nz-1); 
}

//==============================================================================
/// Configures and load data from CSV or BIN file.
//==============================================================================
void JSphInOutGridData::ConfigFromFile(const std::string &filename){
  Reset();
  string ext=fun::StrUpper(fun::GetExtension(filename));
  if(ext=="CSV")LoadDataCsv(filename);
  else if(ext=="BIN")LoadDataBin(filename);
  else Run_ExceptioonFile("Unknown file extension.",filename);
  File=filename;
}

//==============================================================================
/// Configures and load data from CSV file.
//==============================================================================
void JSphInOutGridData::LoadDataCsv(const std::string &filename){
  JReadDatafile rdat;
  rdat.LoadFile(filename);
  //-Load and check fmtversion.
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  //printf("==> rows: %u\n",rows);
  if(rows<2)Run_ExceptioonFile("Number of rows is invalid.",filename);
  if(rdat.ReadNextValue()!="fmtversion")Run_ExceptioonFile("fmtversion is missing.",filename);
  rdat.SetReadLine(1);
  const unsigned fmtver=rdat.ReadNextUnsigned(true);
  if(fmtver!=FmtVersion)Run_ExceptioonFile(fun::PrintStr("fmtversion value (%u) is invalid. The expected format version is %u.",fmtver,FmtVersion),filename);
  //-Load head data according FmtVersion.
  if(rows<5)Run_ExceptioonFile("Number of data rows is invalid.",filename);
  rdat.SetReadLine(0); rdat.ReadNextValue(true);
  bool headerror=false;
  if(rdat.ReadNextValue(true)!="grid_dpx")headerror=true;
  if(rdat.ReadNextValue(true)!="grid_dpz")headerror=true;
  if(rdat.ReadNextValue(true)!="grid_nx") headerror=true;
  if(rdat.ReadNextValue(true)!="grid_nz") headerror=true;
  if(rdat.ReadNextValue(true)!="vars")    headerror=true;
  if(headerror)Run_ExceptioonFile("First line is invalid.",filename);
  rdat.ReadNextValue();
  const double dpx=rdat.ReadNextDouble(true);
  const double dpz=rdat.ReadNextDouble(true);
  const unsigned nx=rdat.ReadNextUnsigned(true);
  const unsigned nz=rdat.ReadNextUnsigned(true);
  const string vars=rdat.ReadNextValue(true);
  bool usevelz=false;
  if(vars=="velx")usevelz=false;
  else if(vars=="velx velz")usevelz=true;
  else Run_ExceptioonFile("Head value \'vars\' is invalid.",filename);
  ConfigGridData(nx,nz,dpx,dpz,usevelz);
  //-Load values data according FmtVersion.
  rdat.SetReadLine(4);
  const unsigned npt=nx*nz;
  float *velx=new float[npt];
  float *velz=(usevelz? new float[npt]: NULL);
  for(unsigned cr=4;cr<rows;cr++){
    const double time=rdat.ReadNextDouble();
    for(unsigned p=0;p<npt;p++)velx[p]=rdat.ReadNextFloat(true);
    if(usevelz)for(unsigned p=0;p<npt;p++)velz[p]=rdat.ReadNextFloat(true);
    AddDataTime(time,npt,velx,velz);
  }
  //SaveDataCsv("dg.csv");
  //-Frees memory.
  delete[] velx; velx=NULL;
  delete[] velz; velz=NULL;
}

//==============================================================================
/// Configures and load data from BIN file.
//==============================================================================
void JSphInOutGridData::LoadDataBin(const std::string &filename){
  Run_Exceptioon("NOT IMPLEMENTED...");
}

//==============================================================================
/// Configures grid parameters.
//==============================================================================
void JSphInOutGridData::ConfigGridData(unsigned nx,unsigned nz,double dpx,double dpz,bool usevelz){
  Reset();
  if(dpx<=0 || dpz<=0 || !nx || !nz || nx>30 || nz>1000)Run_Exceptioon("Grid configuration is invalid.");
  Nx=nx; Nz=nz; Npt=Nx*Nz; 
  Dpx=dpx; Dpz=dpz;
  UseVelz=usevelz;
  SetPosMin(PosMin);
}

//==============================================================================
/// Adds data for another time.
//==============================================================================
void JSphInOutGridData::AddDataTime(double time,unsigned npt,const float *velx,const float *velz){
  if(CountTimes() && DataTimes[CountTimes()-1]->GetTime()>=time)Run_Exceptioon("New time of data is not higher than previous one.");
  if(npt!=Npt)Run_Exceptioon("The number of points does not match.");
  JSphInOutGridDataTime *gdt=new JSphInOutGridDataTime(Nx,Nz,time,velx,(UseVelz? velz: NULL));
  DataTimes.push_back(gdt);
}

//==============================================================================
/// Saves DataTimes in CSV file.
//==============================================================================
void JSphInOutGridData::SaveDataCsv(std::string filename)const{
  filename=fun::GetWithoutExtension(filename)+".csv";
  jcsv::JSaveCsv2 scsv(filename,false,AppInfo.GetCsvSepComa());
  //-Saves head in CSV file.
  scsv.SetHead();
  scsv << "fmtversion;grid_dpx;grid_dpz;grid_nx;grid_nz;vars" << jcsv::Endl();
  const string vars=(UseVelz? "velx velz": "velx");
  scsv << FmtVersion << Dpx << Dpz << Nx << Nz << vars << jcsv::Endl();
  scsv << jcsv::Endl();
  scsv << "time";
  for(unsigned cx=0;cx<Nx;cx++)for(unsigned cz=0;cz<Nz;cz++)scsv << fun::PrintStr("vx_x%u_z%u",cx,cz);
  if(UseVelz)for(unsigned cx=0;cx<Nx;cx++)for(unsigned cz=0;cz<Nz;cz++)scsv << fun::PrintStr("vz_x%u_z%u",cx,cz);
  scsv << jcsv::Endl();
  scsv.SaveData(true);
  //-Saves data in CSV file.
  scsv.SetData();
  const unsigned nt=CountTimes();
  for(unsigned ct=0;ct<nt;ct++){
    const JSphInOutGridDataTime *gdt=DataTimes[ct];
    const float *velx=gdt->GetVelx();
    const float *velz=gdt->GetVelz();
    scsv << gdt->GetTime();
    for(unsigned c=0;c<Npt;c++)scsv << velx[c];
    if(UseVelz)for(unsigned c=0;c<Npt;c++)scsv << velz[c];
    scsv << jcsv::Endl();
  }
  //-Saves file.
  scsv.SaveData(true);
}

//==============================================================================
/// Loads values for time t in SelData object.
//==============================================================================
void JSphInOutGridData::ComputeTime(double t){
  if(!SelData)SelData=new JSphInOutGridDataTime(Nx,Nz);
  if(SelData->GetTime()!=t){
    const unsigned nt=CountTimes();
    if(SelCt>=nt)SelCt=nt-1;
    while(SelCt+1<nt && t>=DataTimes[SelCt+1]->GetTime())SelCt++;
    while(SelCt>0 && t<DataTimes[SelCt-1]->GetTime())SelCt--;
    const double seltime=DataTimes[SelCt]->GetTime();
    if(t<=seltime || (t>=seltime && SelCt+1>=nt)){//-Copy data from DataTimes[SelCt].
      SelData->CopyFrom(t,DataTimes[SelCt]);
    }
    else{ //-Interpolate data between DataTimes[SelCt] and DataTimes[SelCt+1].
      SelData->Interpolate(t,DataTimes[SelCt],DataTimes[SelCt+1]);
    }
  }
}

//==============================================================================
/// Interpolate velocity in time and position of selected partiles in a list.
//==============================================================================
void JSphInOutGridData::InterpolateVelCpu(double time,unsigned izone,unsigned np,const int *plist
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop,float velcorr)
{
  ComputeTime(time);
  const float *velx=SelData->GetVelx();
  const float *velz=SelData->GetVelz();
  const int nx1=Nx-1;
  const int nz1=Nz-1;
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int cp=0;cp<n;cp++){
    const unsigned p=plist[cp];
    if(izone==CODE_GetIzoneFluidInout(code[p])){
      const double px=pos[p].x-PosMin.x;
      const double pz=pos[p].z-PosMin.z;
      int cx=int(px/Dpx);
      cx=max(cx,0);
      cx=min(cx,nx1);
      const double fx=(px/Dpx-cx);  //const double fx=(px-Dpx*cx)/Dpx;
      int cz=int(pz/Dpz);
      cz=max(cz,0);
      cz=min(cz,nz1);
      const double fz=(pz/Dpz-cz);  //const double fz=(pz-Dpz*cz)/Dpz;
      //-Interpolation in Z.
      const unsigned cp=Nz*cx+cz;
      const float v00=velx[cp];
      const float v01=(cz<nz1? velx[cp+1]:    v00);
      const float v10=(cx<nx1? velx[cp+Nz]:   v00);
      const float v11=(cx<nx1? (cz<nz1? velx[cp+Nz+1]: v10): v01);
      const float v0=float(fz*(v01-v00)+v00);
      const float v1=float(fz*(v11-v10)+v10);
      const float v=float(fx*(v1-v0)+v0);
      velrhop[p]=TFloat4(v-velcorr,0,0,velrhop[p].w);
      if(UseVelz){
        const float v00=velz[cp];
        const float v01=(cz<nz1? velz[cp+1]:    v00);
        const float v10=(cx<nx1? velz[cp+Nz]:   v00);
        const float v11=(cx<nx1? (cz<nz1? velz[cp+Nz+1]: v10): v01);
        const float v0=float(fz*(v01-v00)+v00);
        const float v1=float(fz*(v11-v10)+v10);
        const float v=float(fx*(v1-v0)+v0);
        velrhop[p].z=v;
      }
    }
  }
}

//==============================================================================
/// Interpolate velocity in time and Z-position of selected partiles in a list.
//==============================================================================
void JSphInOutGridData::InterpolateZVelCpu(double time,unsigned izone,unsigned np,const int *plist
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,tfloat4 *velrhop,float velcorr)
{
  ComputeTime(time);
  const float *velx=SelData->GetVelx();
  const float *velz=SelData->GetVelz();
  const int nx1=0;
  const int nz1=Nz-1;
  const int n=int(np);
  //#ifdef OMP_USE
  //  #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  //#endif
  for(int cp=0;cp<n;cp++){
    const unsigned p=plist[cp];
    if(izone==CODE_GetIzoneFluidInout(code[p])){
      const double pz=pos[p].z-PosMin.z;
      int cz=int(pz/Dpz);
      cz=max(cz,0);
      cz=min(cz,nz1);
      const double fz=(pz/Dpz-cz);  //const double fz=(pz-Dpz*cz)/Dpz;
      //-Interpolation in Z.
      const unsigned cp=cz;
      const float v00=velx[cp];
      const float v01=(cz<nz1? velx[cp+1]: v00);
      const float v=float(fz*(v01-v00)+v00);
      velrhop[p]=TFloat4(v-velcorr,0,0,velrhop[p].w);
      if(UseVelz){
        const float v00=velz[cp];
        const float v01=(cz<nz1? velz[cp+1]:    v00);
        const float v=float(fz*(v01-v00)+v00);
        velrhop[p].z=v;
      }
    }
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Loads values for time t in SelData object.
//==============================================================================
void JSphInOutGridData::ComputeTimeGpu(double t){
  if(!Velx0g)AllocateMemoryGpu();
  if(TimeGpu!=t){
    const unsigned nt=CountTimes();
    if(SelCt>=nt)SelCt=nt-1;
    while(SelCt+1<nt && t>=DataTimes[SelCt+1]->GetTime())SelCt++;
    while(SelCt>0 && t<DataTimes[SelCt-1]->GetTime())SelCt--;
    const double seltime=DataTimes[SelCt]->GetTime();
    //-Swap data between variables for time t0 and t1.
    if(CtVel0!=SelCt && CtVel1==SelCt){
      swap(CtVel0,CtVel1);
      swap(Velx0g,Velx1g);
      swap(Velz0g,Velz1g);
    }
    //-Updata data for time t0.
    if(CtVel0!=SelCt){
      CtVel0=SelCt;
      cudaMemcpy(Velx0g,DataTimes[SelCt]->GetVelx(),sizeof(float)*Npt,cudaMemcpyHostToDevice);
      if(UseVelz)cudaMemcpy(Velz0g,DataTimes[SelCt]->GetVelz(),sizeof(float)*Npt,cudaMemcpyHostToDevice);
    }
    //-Updata data for time t1.
    unsigned selct1=(SelCt+1<nt? SelCt+1: SelCt);
    if(CtVel1!=selct1){
      CtVel1=selct1;
      cudaMemcpy(Velx1g,DataTimes[selct1]->GetVelx(),sizeof(float)*Npt,cudaMemcpyHostToDevice);
      if(UseVelz)cudaMemcpy(Velz1g,DataTimes[selct1]->GetVelz(),sizeof(float)*Npt,cudaMemcpyHostToDevice);
    }
    //-Updata data for the requested time (t).
    if(t<=seltime || (t>=seltime && SelCt+1>=nt)){//-Copy data from DataTimes[SelCt].
      if(CtSelVel!=SelCt){
        CtSelVel=SelCt;
        cudaMemcpy(SelVelxg,Velx0g,sizeof(float)*Npt,cudaMemcpyDeviceToDevice);
        if(UseVelz)cudaMemcpy(SelVelzg,Velz0g,sizeof(float)*Npt,cudaMemcpyDeviceToDevice);
      }
    }
    else{ //-Interpolate data between DataTimes[SelCt] and DataTimes[SelCt+1].
      CtSelVel=UINT_MAX;
      const double t0=DataTimes[CtVel0]->GetTime();
      const double t1=DataTimes[CtVel1]->GetTime();
      cusphinout::InOutInterpolateTime(Npt,t,t0,t1,Velx0g,Velx1g,SelVelxg,Velz0g,Velz1g,SelVelzg);
    }
    TimeGpu=t;
  }
}

//==============================================================================
/// Interpolate velocity in time and Z-position of selected partiles in a list.
//==============================================================================
void JSphInOutGridData::InterpolateZVelGpu(double time,unsigned izone,unsigned np,const int *plist
  ,const double2 *posxyg,const double *poszg,const typecode *codeg,const unsigned *idpg
  ,float4 *velrhopg,float velcorr)
{
  ComputeTimeGpu(time);
  cusphinout::InOutInterpolateZVel(izone,PosMin.z,Dpz,Nz-1,SelVelxg,SelVelzg,np,plist,poszg,codeg,velrhopg,velcorr);
}
#endif

//==============================================================================
/// Saves DataTime object in VTK file.
//==============================================================================
void JSphInOutGridData::SaveVtk(const JSphInOutGridDataTime *gdt,std::string filename)const{
  const tfloat3 pos0=ToTFloat3(PosMin);
  //-Allocates memory.
  tfloat3 *pos=new tfloat3[Npt];
  tfloat3 *vel=new tfloat3[Npt];
  //-Computes position.
  unsigned p=0;
  for(unsigned cx=0;cx<Nx;cx++)for(unsigned cz=0;cz<Nz;cz++,p++){
    pos[p]=pos0+TFloat3(float(Dpx*cx),0,float(Dpz*cz));
  }
  //-Computes velocity.
  const float *velx=gdt->GetVelx();
  const float *velz=gdt->GetVelz();
  if(UseVelz)for(unsigned p=0;p<Npt;p++)vel[p]=TFloat3(velx[p],0,velz[p]);
  else       for(unsigned p=0;p<Npt;p++)vel[p]=TFloat3(velx[p],0,0);
  //-Saves VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos",Npt,pos,true);
  arrays.AddArray("Vel",Npt,vel,true);
  JVtkLib::SaveVtkData(filename,arrays,"Pos");
  arrays.Reset();
  ////-Old style...
  //std::vector<JFormatFiles2::StScalarData> fields;
  //fields.push_back(JFormatFiles2::DefineField("Vel",JFormatFiles2::Float32,3,vel));
  //JFormatFiles2::SaveVtk(filename,Npt,pos,fields);
  ////-Frees memory.
  //delete[] pos; pos=NULL;
  //delete[] vel; vel=NULL;
}

//==============================================================================
/// Saves DataTimes in VTK file.
//==============================================================================
void JSphInOutGridData::SaveDataVtk(std::string filename,int ctime)const{
  const bool onefile=(ctime!=-1);
  filename=fun::GetWithoutExtension(filename)+".vtk";
  const unsigned ctini=(!onefile? 0: unsigned(ctime));
  const unsigned ctfin=(!onefile? CountTimes(): unsigned(ctime)+1);
  if(ctini>=CountTimes())Run_Exceptioon("Number of DataTime is invalid.");
  if(onefile)SaveVtk(DataTimes[ctini],filename);
  else for(unsigned ct=ctini;ct<ctfin;ct++)SaveVtk(DataTimes[ct],fun::FileNameSec(filename,ct));
}

//==============================================================================
/// Saves interpolated values at time t in a VTK file.
//==============================================================================
void JSphInOutGridData::SaveDataVtkTime(std::string filename,double tmax,double dt){
  filename=fun::GetWithoutExtension(filename)+".vtk";
  const unsigned nt=max(unsigned(ceil(tmax/dt))+1,1u);
  for(unsigned ct=0;ct<nt;ct++){
    //-Calculates data for each time.
    ComputeTime(dt*ct);
    SaveVtk(SelData,fun::FileNameSec(filename,ct));
  }
}

////==============================================================================
///// Saves grid nodes in VTK file.
////==============================================================================
//void JSphInOutGridData::SaveVtkGrid(std::string filename,tfloat3 pos0)const{
//  std::vector<JFormatFiles2::StShapeData> shapes;
//  tfloat3 pmin=pos0;
//  tfloat3 pmax=pos0+TFloat3(float(Dpx*(Nx-1)),0,float(Dpz*(Nz-1)));
//  //-Vertical lines.
//  for(unsigned c=0;c<Nx;c++){
//    const float px=pmin.x+float(Dpx*c);
//    shapes.push_back(JFormatFiles2::DefineShape_Line(TFloat3(px,pmin.y,pmin.z),TFloat3(px,pmin.y,pmax.z),0,0));
//  }
//  //-Horizontal lines.
//  for(unsigned c=0;c<Nx;c++){
//    const float pz=pmin.z+float(Dpz*c);
//    shapes.push_back(JFormatFiles2::DefineShape_Line(TFloat3(px,pmin.y,pmin.z),TFloat3(px,pmin.y,pmax.z),0,0));
//  }
//
//
//    for(unsigned ci=0;ci<GetCount();ci++){
//      const JSphInOutZone *izone=List[ci];
//      const tdouble3* ptdom=izone->GetPtDomain();
//      if(Simulate2D){
//        shapes.push_back(JFormatFiles2::DefineShape_Quad(ptdom[0],ptdom[1],ptdom[2],ptdom[3],ci,0));
//
//
//
//  const bool onefile=(ctime!=-1);
//  filename=fun::GetWithoutExtension(filename)+".vtk";
//  const unsigned ctini=(!onefile? 0: unsigned(ctime));
//  const unsigned ctfin=(!onefile? CountTimes(): unsigned(ctime)+1);
//  if(ctini>=CountTimes())Run_Exceptioon("Number of DataTime is invalid.");
//  if(onefile)SaveVtk(DataTimes[ctini],filename,pos0);
//  else for(unsigned ct=ctini;ct<ctfin;ct++)SaveVtk(DataTimes[ct],fun::FileNameSec(filename,ct),pos0);
//}



