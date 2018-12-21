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

/// \file JSph.cpp \brief Implements the class \ref JSph.

#include "JSph.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "JPartDataHead.h"
#include "JSphMk.h"
#include "JSphPartsInit.h"
#include "JPartsLoad4.h"
#include "JSphMotion.h"
#include "JXml.h"
#include "JSpaceCtes.h"
#include "JSpaceEParms.h"
#include "JSpaceParts.h"
#include "JFormatFiles2.h"
#include "JSphDtFixed.h"
#include "JSaveDt.h"
#include "JTimeOut.h"
#include "JSphVisco.h"
#include "JGaugeSystem.h"
#include "JWaveGen.h"
#include "JSphAccInput.h"
#include "JPartDataBi4.h"
#include "JPartOutBi4Save.h"
#include "JPartFloatBi4.h"
#include "JPartsOut.h"
#include "JDamping.h"
#include "JSphInitialize.h"
#include <climits>

//using namespace std;
using std::string;
using std::ofstream;
using std::endl;

//==============================================================================
/// Constructor.
//==============================================================================
JSph::JSph(bool cpu,bool withmpi):Cpu(cpu),WithMpi(withmpi){
  ClassName="JSph";
  DgNum=0;
  DataBi4=NULL;
  DataOutBi4=NULL;
  DataFloatBi4=NULL;
  PartsOut=NULL;
  Log=NULL;
  ViscoTime=NULL;
  DtFixed=NULL;
  SaveDt=NULL;
  TimeOut=NULL;
  MkInfo=NULL;
  PartsInit=NULL;
  SphMotion=NULL;
  FtObjs=NULL;
  DemData=NULL;
  GaugeSystem=NULL;
  WaveGen=NULL;
  Damping=NULL;
  AccInput=NULL;
  PartsLoaded=NULL;
  InitVars();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSph::~JSph(){
  DestructorActive=true;
  delete DataBi4;       DataBi4=NULL;
  delete DataOutBi4;    DataOutBi4=NULL;
  delete DataFloatBi4;  DataFloatBi4=NULL;
  delete PartsOut;      PartsOut=NULL;
  delete ViscoTime;     ViscoTime=NULL;
  delete DtFixed;       DtFixed=NULL;
  delete SaveDt;        SaveDt=NULL;
  delete TimeOut;       TimeOut=NULL;
  delete MkInfo;        MkInfo=NULL;
  delete PartsInit;     PartsInit=NULL;
  delete SphMotion;     SphMotion=NULL;
  AllocMemoryFloating(0);
  delete[] DemData;     DemData=NULL;
  delete GaugeSystem;   GaugeSystem=NULL;
  delete WaveGen;       WaveGen=NULL;
  delete Damping;       Damping=NULL;
  delete AccInput;      AccInput=NULL; 
  delete PartsLoaded;   PartsLoaded=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSph::InitVars(){
  ClearCfgDomain();
  OutPosCount=OutRhopCount=OutMoveCount=0;
  Simulate2D=false;
  Simulate2DPosY=0;
  Stable=false;
  Psingle=true;
  SvDouble=false;
  RunCode=CalcRunCode();
  RunTimeDate="";
  CaseName=""; DirCase=""; RunName="";
  DirOut=""; 
  DirDataOut=""; 
  FileXml="";
  TStep=STEP_None;
  InterStep=INTERSTEP_None;
  VerletSteps=40;
  TKernel=KERNEL_Wendland;
  Awen=Bwen=Agau=Bgau=0;
  memset(&CubicCte,0,sizeof(StCubicCte));
  TVisco=VISCO_None;
  TDeltaSph=DELTA_None; DeltaSph=0;
  TShifting=SHIFT_None; ShiftCoef=ShiftTFS=0;
  Visco=0; ViscoBoundFactor=1;
  UseDEM=false;  //(DEM)
  delete[] DemData; DemData=NULL;  //(DEM)
  RhopOut=true; RhopOutMin=700; RhopOutMax=1300;
  TimeMax=TimePart=0;
  DtIni=DtMin=0; CoefDtMin=0; DtAllParticles=false;
  PartsOutMax=0;
  NpMinimum=0;
  PartsOutWrn=1; PartsOutTotWrn=10;

  SvData=byte(SDAT_Binx)|byte(SDAT_Info);
  SvRes=false;
  SvTimers=false;
  SvDomainVtk=false;

  H=CteB=Gamma=RhopZero=CFLnumber=0;
  Dp=0;
  Cs0=0;
  Delta2H=0;
  MassFluid=MassBound=0;
  Gravity=TFloat3(0);
  Dosh=H2=Fourh2=Eta2=0;
  SpsSmag=SpsBlin=0;

  CasePosMin=CasePosMax=TDouble3(0);
  CaseNp=CaseNbound=CaseNfixed=CaseNmoving=CaseNfloat=CaseNfluid=CaseNpb=0;

  PeriActive=0; PeriX=PeriY=PeriZ=false;
  PeriXinc=PeriYinc=PeriZinc=TDouble3(0);

  PartBeginDir=""; 
  PartBegin=PartBeginFirst=0;
  PartBeginTimeStep=0; 
  PartBeginTotalNp=0;

  FtCount=0;
  FtPause=0;
  FtMode=FTMODE_None;
  WithFloating=false;

  AllocMemoryFloating(0);

  CellMode=CELLMODE_None;
  Hdiv=0;
  Scell=0;
  MovLimit=0;

  Map_PosMin=Map_PosMax=Map_Size=TDouble3(0);
  Map_Cells=TUint3(0);
  MapRealPosMin=MapRealPosMax=MapRealSize=TDouble3(0);

  DomCelIni=DomCelFin=TUint3(0);
  DomCells=TUint3(0);
  DomPosMin=DomPosMax=DomSize=TDouble3(0);
  DomRealPosMin=DomRealPosMax=TDouble3(0);
  DomCellCode=0;

  NpDynamic=ReuseIds=false;
  TotalNp=0; IdMax=0;

  DtModif=0;
  DtModifWrn=1;
  PartDtMin=DBL_MAX; PartDtMax=-DBL_MAX;

  MaxMemoryCpu=MaxMemoryGpu=MaxParticles=MaxCells=0;

  PartIni=Part=0; 
  Nstep=0; PartNstep=-1;
  PartOut=0;

  TimeStepIni=0;
  TimeStep=TimeStepM1=0;
  TimePartNext=0;

  VerletStep=0;
  SymplecticDtPre=0;
  DemDtForce=0;  //(DEM)
}

//==============================================================================
/// Generates a random code to identify the file of the results of the execution.
//==============================================================================
std::string JSph::CalcRunCode()const{
  srand((unsigned)time(NULL));
  const unsigned len=8;
  char code[len+1];
  for(unsigned c=0;c<len;c++){
    char let=char(float(rand())/float(RAND_MAX)*36);
    code[c]=(let<10? let+48: let+87);
  } 
  code[len]=0;
  return(code);
}

//==============================================================================
/// Sets the configuration of the domain limits by default.
//==============================================================================
void JSph::ClearCfgDomain(){
  CfgDomainParticlesMin=CfgDomainParticlesMax=TDouble3(DBL_MAX);
  CfgDomainParticlesPrcMin=CfgDomainParticlesPrcMax=TDouble3(DBL_MAX);
  CfgDomainFixedMin=CfgDomainFixedMax=TDouble3(DBL_MAX);
}

//==============================================================================
/// Sets the configuration of the domain limits using given values.
//==============================================================================
void JSph::ConfigDomainFixed(tdouble3 vmin,tdouble3 vmax){
  ClearCfgDomain();
  CfgDomainFixedMin=vmin; CfgDomainFixedMax=vmax;
}

//==============================================================================
/// Sets the configuration of the domain limits using given values.
//==============================================================================
void JSph::ConfigDomainFixedValue(std::string key,double v){
  const char met[]="ConfigDomainFixedValue";
  const string keyend=(key.size()>=4? key.substr(key.size()-4,4): "");
       if(keyend=="Xmin")CfgDomainFixedMin.x=v;
  else if(keyend=="Ymin")CfgDomainFixedMin.y=v;
  else if(keyend=="Zmin")CfgDomainFixedMin.z=v;
  else if(keyend=="Xmax")CfgDomainFixedMax.x=v;
  else if(keyend=="Ymax")CfgDomainFixedMax.y=v;
  else if(keyend=="Zmax")CfgDomainFixedMax.z=v;
  else RunException(met,"Key for limit is invalid.");
}

//==============================================================================
/// Sets the configuration of the domain limits using positions of particles.
//==============================================================================
void JSph::ConfigDomainParticles(tdouble3 vmin,tdouble3 vmax){
  CfgDomainParticlesMin=vmin; CfgDomainParticlesMax=vmax;
}

//==============================================================================
/// Sets the configuration of the domain limits using positions of particles.
//==============================================================================
void JSph::ConfigDomainParticlesValue(std::string key,double v){
  const char met[]="ConfigDomainParticlesValue";
  const string keyend=(key.size()>=4? key.substr(key.size()-4,4): "");
       if(keyend=="Xmin")CfgDomainParticlesMin.x=v;
  else if(keyend=="Ymin")CfgDomainParticlesMin.y=v;
  else if(keyend=="Zmin")CfgDomainParticlesMin.z=v;
  else if(keyend=="Xmax")CfgDomainParticlesMax.x=v;
  else if(keyend=="Ymax")CfgDomainParticlesMax.y=v;
  else if(keyend=="Zmax")CfgDomainParticlesMax.z=v;
  else RunException(met,"Key for limit is invalid.");
}

//==============================================================================
/// Sets the configuration of the domain limits using positions plus a percentage.
//==============================================================================
void JSph::ConfigDomainParticlesPrc(tdouble3 vmin,tdouble3 vmax){
  CfgDomainParticlesPrcMin=vmin; CfgDomainParticlesPrcMax=vmax;
}

//==============================================================================
/// Sets the configuration of the domain limits using positions plus a percentage.
//==============================================================================
void JSph::ConfigDomainParticlesPrcValue(std::string key,double v){
  const char met[]="ConfigDomainParticlesValue";
  const string keyend=(key.size()>=4? key.substr(key.size()-4,4): "");
       if(keyend=="Xmin")CfgDomainParticlesPrcMin.x=v;
  else if(keyend=="Ymin")CfgDomainParticlesPrcMin.y=v;
  else if(keyend=="Zmin")CfgDomainParticlesPrcMin.z=v;
  else if(keyend=="Xmax")CfgDomainParticlesPrcMax.x=v;
  else if(keyend=="Ymax")CfgDomainParticlesPrcMax.y=v;
  else if(keyend=="Zmax")CfgDomainParticlesPrcMax.z=v;
  else RunException(met,"Key for limit is invalid.");
}

//==============================================================================
/// Allocates memory of floating objectcs.
//==============================================================================
void JSph::AllocMemoryFloating(unsigned ftcount){
  delete[] FtObjs; FtObjs=NULL;
  if(ftcount)FtObjs=new StFloatingData[ftcount];
}

//==============================================================================
/// Returns the allocated memory in CPU.
//==============================================================================
llong JSph::GetAllocMemoryCpu()const{  
  //-Allocated in AllocMemoryCase().
  llong s=0;
  //-Allocated in AllocMemoryFloating().
  if(FtObjs)s+=sizeof(StFloatingData)*FtCount;
  //-Allocated in other objects.
  if(PartsOut)s+=PartsOut->GetAllocMemory();
  if(ViscoTime)s+=ViscoTime->GetAllocMemory();
  if(DtFixed)s+=DtFixed->GetAllocMemory();
  if(AccInput)s+=AccInput->GetAllocMemory();
  if(PartsLoaded)s+=PartsLoaded->GetAllocMemory();
  return(s);
}

//==============================================================================
/// Loads the configuration of the execution.
//==============================================================================
void JSph::LoadConfig(const JCfgRun *cfg){
  const char* met="LoadConfig";
  TimerTot.Start();
  Stable=cfg->Stable;
  Psingle=true; SvDouble=false; //-Options by default.
  DirOut=fun::GetDirWithSlash(cfg->DirOut);
  DirDataOut=(!cfg->DirDataOut.empty()? fun::GetDirWithSlash(DirOut+cfg->DirDataOut): DirOut);
  CaseName=cfg->CaseName; 
  DirCase=fun::GetDirWithSlash(fun::GetDirParent(CaseName));
  CaseName=CaseName.substr(DirCase.length());
  if(!CaseName.length())RunException(met,"Name of the case for execution was not indicated.");
  RunName=(cfg->RunName.length()? cfg->RunName: CaseName);
  FileXml=DirCase+CaseName+".xml";
  PartBeginDir=cfg->PartBeginDir; PartBegin=cfg->PartBegin; PartBeginFirst=cfg->PartBeginFirst;

  //-Output options:
  CsvSepComa=cfg->CsvSepComa;
  SvData=byte(SDAT_None); 
  if(cfg->Sv_Csv&&!WithMpi)SvData|=byte(SDAT_Csv);
  if(cfg->Sv_Binx)SvData|=byte(SDAT_Binx);
  if(cfg->Sv_Info)SvData|=byte(SDAT_Info);
  if(cfg->Sv_Vtk)SvData|=byte(SDAT_Vtk);

  SvRes=cfg->SvRes;
  SvTimers=cfg->SvTimers;
  SvDomainVtk=cfg->SvDomainVtk;

  printf("\n");
  RunTimeDate=fun::GetDateTime();
  Log->Printf("[Initialising %s  %s]",ClassName.c_str(),RunTimeDate.c_str());

  const string runpath=AppInfo.GetRunPath();
  Log->Printf("ProgramFile=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,AppInfo.GetRunCommand()),3).c_str());
  Log->Printf("ExecutionDir=\"%s\"",fun::GetPathLevels(runpath,3).c_str());
  Log->Printf("XmlFile=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,FileXml),3).c_str());
  Log->Printf("OutputDir=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,DirOut),3).c_str());
  Log->Printf("OutputDataDir=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,DirDataOut),3).c_str());
  
  //-Creates Output directories when it is necessary.
  if(AppInfo.GetCreateDirs()){
    fun::MkdirPath(DirOut);
    if(DirOut!=DirDataOut)fun::MkdirPath(DirDataOut);
  }

  if(PartBegin){
    Log->Print(fun::VarStr("PartBegin",PartBegin));
    Log->Print(fun::VarStr("PartBeginDir",PartBeginDir));
    Log->Print(fun::VarStr("PartBeginFirst",PartBeginFirst));
  }

  LoadCaseConfig();

  //-Aplies configuration using command line.
  if(cfg->PosDouble==0){      Psingle=true;  SvDouble=false; }
  else if(cfg->PosDouble==1){ Psingle=false; SvDouble=false; }
  else if(cfg->PosDouble==2){ Psingle=false; SvDouble=true;  }
  if(cfg->TStep)TStep=cfg->TStep;
  if(cfg->VerletSteps>=0)VerletSteps=cfg->VerletSteps;
  if(cfg->TKernel)TKernel=cfg->TKernel;
  if(cfg->TVisco){ TVisco=cfg->TVisco; Visco=cfg->Visco; }
  if(cfg->ViscoBoundFactor>=0)ViscoBoundFactor=cfg->ViscoBoundFactor;
  if(cfg->DeltaSph>=0){
    DeltaSph=cfg->DeltaSph;
    TDeltaSph=(DeltaSph? DELTA_Dynamic: DELTA_None);
  }
  if(TDeltaSph==DELTA_Dynamic && Cpu)TDeltaSph=DELTA_DynamicExt; //-It is necessary because the interaction is divided in two steps: fluid-fluid/float and fluid-bound.

  if(cfg->Shifting>=0){
    switch(cfg->Shifting){
      case 0:  TShifting=SHIFT_None;     break;
      case 1:  TShifting=SHIFT_NoBound;  break;
      case 2:  TShifting=SHIFT_NoFixed;  break;
      case 3:  TShifting=SHIFT_Full;     break;
      default: RunException(met,"Shifting mode is not valid.");
    }
    if(TShifting!=SHIFT_None){
      ShiftCoef=-2; ShiftTFS=0;
    }
    else ShiftCoef=ShiftTFS=0;
  }

  if(cfg->FtPause>=0)FtPause=cfg->FtPause;
  if(cfg->TimeMax>0)TimeMax=cfg->TimeMax;
  //-Configuration of JTimeOut with TimePart.
  TimeOut=new JTimeOut();
  if(cfg->TimePart>=0){
    TimePart=cfg->TimePart;
    TimeOut->Config(TimePart);
  }
  else TimeOut->Config(FileXml,"case.execution.special.timeout",TimePart);

  CellMode=cfg->CellMode;
  if(cfg->DomainMode==1){
    ConfigDomainParticles(cfg->DomainParticlesMin,cfg->DomainParticlesMax);
    ConfigDomainParticlesPrc(cfg->DomainParticlesPrcMin,cfg->DomainParticlesPrcMax);
  }
  else if(cfg->DomainMode==2)ConfigDomainFixed(cfg->DomainFixedMin,cfg->DomainFixedMax);
  if(cfg->RhopOutModif){
    RhopOutMin=cfg->RhopOutMin; RhopOutMax=cfg->RhopOutMax;
  }
  RhopOut=(RhopOutMin<RhopOutMax);
  if(!RhopOut){ RhopOutMin=-FLT_MAX; RhopOutMax=FLT_MAX; }
}

//==============================================================================
/// Loads the case configuration to be executed.
//==============================================================================
void JSph::LoadCaseConfig(){
  const char* met="LoadCaseConfig";
  if(!fun::FileExists(FileXml))RunException(met,"Case configuration was not found.",FileXml);
  JXml xml; xml.LoadFile(FileXml);
  JSpaceCtes ctes;     ctes.LoadXmlRun(&xml,"case.execution.constants");
  JSpaceEParms eparms; eparms.LoadXml(&xml,"case.execution.parameters");
  JSpaceParts parts;   parts.LoadXml(&xml,"case.execution.particles");

  //-Execution parameters.
  switch(eparms.GetValueInt("PosDouble",true,0)){
    case 0:  Psingle=true;  SvDouble=false;  break;
    case 1:  Psingle=false; SvDouble=false;  break;
    case 2:  Psingle=false; SvDouble=true;   break;
    default: RunException(met,"PosDouble value is not valid.");
  }
  switch(eparms.GetValueInt("RigidAlgorithm",true,1)){ //(DEM)
    case 1:  UseDEM=false;    break;
    case 2:  UseDEM=true;     break;
    default: RunException(met,"Rigid algorithm is not valid.");
  }
  switch(eparms.GetValueInt("StepAlgorithm",true,1)){
    case 1:  TStep=STEP_Verlet;      break;
    case 2:  TStep=STEP_Symplectic;  break;
    default: RunException(met,"Step algorithm is not valid.");
  }
  VerletSteps=eparms.GetValueInt("VerletSteps",true,40);
  switch(eparms.GetValueInt("Kernel",true,2)){
    case 1:  TKernel=KERNEL_Cubic;     break;
    case 2:  TKernel=KERNEL_Wendland;  break;
    case 3:  TKernel=KERNEL_Gaussian;  break;
    default: RunException(met,"Kernel choice is not valid.");
  }
  switch(eparms.GetValueInt("ViscoTreatment",true,1)){
    case 1:  TVisco=VISCO_Artificial;  break;
    case 2:  TVisco=VISCO_LaminarSPS;  break;
    default: RunException(met,"Viscosity treatment is not valid.");
  }
  Visco=eparms.GetValueFloat("Visco");
  ViscoBoundFactor=eparms.GetValueFloat("ViscoBoundFactor",true,1.f);
  string filevisco=eparms.GetValueStr("ViscoTime",true);
  if(!filevisco.empty()){
    ViscoTime=new JSphVisco();
    ViscoTime->LoadFile(DirCase+filevisco);
  }
  DeltaSph=eparms.GetValueFloat("DeltaSPH",true,0);
  TDeltaSph=(DeltaSph? DELTA_Dynamic: DELTA_None);

  switch(eparms.GetValueInt("Shifting",true,0)){
    case 0:  TShifting=SHIFT_None;     break;
    case 1:  TShifting=SHIFT_NoBound;  break;
    case 2:  TShifting=SHIFT_NoFixed;  break;
    case 3:  TShifting=SHIFT_Full;     break;
    default: RunException(met,"Shifting mode is not valid.");
  }
  if(TShifting!=SHIFT_None){
    ShiftCoef=eparms.GetValueFloat("ShiftCoef",true,-2);
    if(ShiftCoef==0)TShifting=SHIFT_None;
    else ShiftTFS=eparms.GetValueFloat("ShiftTFS",true,0);
  }

  FtPause=eparms.GetValueFloat("FtPause",true,0);
  TimeMax=eparms.GetValueDouble("TimeMax");
  TimePart=eparms.GetValueDouble("TimeOut");

  DtIni=eparms.GetValueDouble("DtIni",true,0);
  DtMin=eparms.GetValueDouble("DtMin",true,0);
  CoefDtMin=eparms.GetValueFloat("CoefDtMin",true,0.05f);
  DtAllParticles=(eparms.GetValueInt("DtAllParticles",true,0)==1);

  string filedtfixed=eparms.GetValueStr("DtFixed",true);
  if(!filedtfixed.empty()){
    DtFixed=new JSphDtFixed();
    DtFixed->LoadFile(DirCase+filedtfixed);
  }
  if(eparms.Exists("RhopOutMin"))RhopOutMin=eparms.GetValueFloat("RhopOutMin");
  if(eparms.Exists("RhopOutMax"))RhopOutMax=eparms.GetValueFloat("RhopOutMax");
  PartsOutMax=eparms.GetValueFloat("PartsOutMax",true,1);

  //-Configuration of periodic boundaries.
  {
    PeriX=PeriY=PeriZ=false;
    PeriXinc=PeriYinc=PeriZinc=TDouble3(0);
    if(eparms.Exists("XPeriodicIncY")){ PeriXinc.y=eparms.GetValueDouble("XPeriodicIncY"); PeriX=true; }
    if(eparms.Exists("XPeriodicIncZ")){ PeriXinc.z=eparms.GetValueDouble("XPeriodicIncZ"); PeriX=true; }
    if(eparms.Exists("YPeriodicIncX")){ PeriYinc.x=eparms.GetValueDouble("YPeriodicIncX"); PeriY=true; }
    if(eparms.Exists("YPeriodicIncZ")){ PeriYinc.z=eparms.GetValueDouble("YPeriodicIncZ"); PeriY=true; }
    if(eparms.Exists("ZPeriodicIncX")){ PeriZinc.x=eparms.GetValueDouble("ZPeriodicIncX"); PeriZ=true; }
    if(eparms.Exists("ZPeriodicIncY")){ PeriZinc.y=eparms.GetValueDouble("ZPeriodicIncY"); PeriZ=true; }
    if(eparms.Exists("XYPeriodic")){ PeriX=true;  PeriY=true;  PeriZ=false;  PeriXinc=PeriYinc=PeriZinc=TDouble3(0); }
    if(eparms.Exists("XZPeriodic")){ PeriX=true;  PeriY=false; PeriZ=true;   PeriXinc=PeriYinc=PeriZinc=TDouble3(0); }
    if(eparms.Exists("YZPeriodic")){ PeriX=false; PeriY=true;  PeriZ=true;   PeriXinc=PeriYinc=PeriZinc=TDouble3(0); }
    PeriActive=DefPeriActive(PeriX,PeriY,PeriZ);
  }

  //-Configuration of domain size.
  float incz=eparms.GetValueFloat("IncZ",true,0.f);
  if(incz){
    ClearCfgDomain();
    CfgDomainParticlesPrcMax.z=incz;
  }
  string key;
  if(eparms.Exists(key="DomainParticles"))ConfigDomainParticles(TDouble3(eparms.GetValueNumDouble(key,0),eparms.GetValueNumDouble(key,1),eparms.GetValueNumDouble(key,2)),TDouble3(eparms.GetValueNumDouble(key,3),eparms.GetValueNumDouble(key,4),eparms.GetValueNumDouble(key,5)));
  if(eparms.Exists(key="DomainParticlesXmin"))ConfigDomainParticlesValue(key,-eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesYmin"))ConfigDomainParticlesValue(key,-eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesZmin"))ConfigDomainParticlesValue(key,-eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesXmax"))ConfigDomainParticlesValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesYmax"))ConfigDomainParticlesValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesZmax"))ConfigDomainParticlesValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesPrc"))ConfigDomainParticlesPrc(TDouble3(eparms.GetValueNumDouble(key,0),eparms.GetValueNumDouble(key,1),eparms.GetValueNumDouble(key,2)),TDouble3(eparms.GetValueNumDouble(key,3),eparms.GetValueNumDouble(key,4),eparms.GetValueNumDouble(key,5)));
  if(eparms.Exists(key="DomainParticlesPrcXmin"))ConfigDomainParticlesPrcValue(key,-eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesPrcYmin"))ConfigDomainParticlesPrcValue(key,-eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesPrcZmin"))ConfigDomainParticlesPrcValue(key,-eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesPrcXmax"))ConfigDomainParticlesPrcValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesPrcYmax"))ConfigDomainParticlesPrcValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainParticlesPrcZmax"))ConfigDomainParticlesPrcValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainFixed"))ConfigDomainFixed(TDouble3(eparms.GetValueNumDouble(key,0),eparms.GetValueNumDouble(key,1),eparms.GetValueNumDouble(key,2)),TDouble3(eparms.GetValueNumDouble(key,3),eparms.GetValueNumDouble(key,4),eparms.GetValueNumDouble(key,5)));
  if(eparms.Exists(key="DomainFixedXmin"))ConfigDomainFixedValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainFixedYmin"))ConfigDomainFixedValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainFixedZmin"))ConfigDomainFixedValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainFixedXmax"))ConfigDomainFixedValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainFixedYmax"))ConfigDomainFixedValue(key,eparms.GetValueDouble(key));
  if(eparms.Exists(key="DomainFixedZmax"))ConfigDomainFixedValue(key,eparms.GetValueDouble(key));

  //-Predefined constantes.
  if(ctes.GetEps()!=0)Log->PrintWarning("Eps value is not used (this correction is deprecated).");
  H=(float)ctes.GetH();
  CteB=(float)ctes.GetB();
  Gamma=(float)ctes.GetGamma();
  RhopZero=(float)ctes.GetRhop0();
  CFLnumber=(float)ctes.GetCFLnumber();
  Dp=ctes.GetDp();
  Gravity=ToTFloat3(ctes.GetGravity());
  MassFluid=(float)ctes.GetMassFluid();
  MassBound=(float)ctes.GetMassBound();

  //-Particle data.
  CaseNp=parts.Count();
  CaseNfixed=parts.Count(TpPartFixed);
  CaseNmoving=parts.Count(TpPartMoving);
  CaseNfloat=parts.Count(TpPartFloating);
  CaseNfluid=parts.Count(TpPartFluid);
  CaseNbound=CaseNp-CaseNfluid;
  CaseNpb=CaseNbound-CaseNfloat;
  
  NpDynamic=ReuseIds=false;
  TotalNp=CaseNp; IdMax=CaseNp-1;

  //-Loads and configures MK of particles.
  MkInfo=new JSphMk();
  MkInfo->Config(&parts);

  //-Configuration of GaugeSystem.
  GaugeSystem=new JGaugeSystem(Cpu,Log);

  //-Configuration of AccInput.
  if(xml.GetNode("case.execution.special.accinputs",false)){
    AccInput=new JSphAccInput(Log,DirCase,&xml,"case.execution.special.accinputs");
  }

  //-Loads and configures moving objects.
  if(parts.CountBlocks(TpPartMoving)>0){
    SphMotion=new JSphMotion();
    SphMotion->Init(&parts,&xml,"case.execution.motion",DirCase);
  }

  //-Configuration of WaveGen.
  if(xml.GetNode("case.execution.special.wavepaddles",false)){
    bool useomp=false,usegpu=false;
    #ifdef OMP_USE_WAVEGEN
      useomp=(omp_get_max_threads()>1);
    #endif
    WaveGen=new JWaveGen(useomp,!Cpu,Log,DirCase,&xml,"case.execution.special.wavepaddles");
    if(SphMotion)for(unsigned c=0;c<SphMotion->GetNumObjects();c++){
      WaveGen->ConfigPaddle(SphMotion->GetObjMkBound(c),c,SphMotion->GetObjBegin(c),SphMotion->GetObjSize(c));
    }
  }

  //-Configuration of damping zones.
  if(xml.GetNode("case.execution.special.damping",false)){
    Damping=new JDamping(Dp,Log);
    Damping->LoadXml(&xml,"case.execution.special.damping");
  }

  //-Loads floating objects.
  FtCount=parts.CountBlocks(TpPartFloating);
  if(FtCount){
    if(FtCount>CODE_MKRANGEMAX)RunException(met,"The number of floating objects exceeds the maximum.");
    AllocMemoryFloating(FtCount);
    unsigned cobj=0;
    for(unsigned c=0;c<parts.CountBlocks()&&cobj<FtCount;c++){
      const JSpacePartBlock &block=parts.GetBlock(c);
      if(block.Type==TpPartFloating){
        const JSpacePartBlock_Floating &fblock=(const JSpacePartBlock_Floating &)block;
        StFloatingData* fobj=FtObjs+cobj;
        fobj->mkbound=fblock.GetMkType();
        fobj->begin=fblock.GetBegin();
        fobj->count=fblock.GetCount();
        fobj->mass=(float)fblock.GetMassbody();
        fobj->massp=fobj->mass/fobj->count;
        fobj->radius=0;
        fobj->center=fblock.GetCenter();
        fobj->angles=TFloat3(0);
        fobj->fvel=ToTFloat3(fblock.GetVelini());
        fobj->fomega=ToTFloat3(fblock.GetOmegaini());
        fobj->inertiaini=ToTMatrix3f(fblock.GetInertia());
        cobj++;
      }
    }
  }
  else{
    if(UseDEM){
      UseDEM=false;
      Log->PrintWarning("The use of DEM was disabled because there are no floating objects...");
    }
  }
  WithFloating=(FtCount>0);
  FtMode=(WithFloating? FTMODE_Sph: FTMODE_None);
  if(UseDEM)FtMode=FTMODE_Ext;

  //-Loads DEM data for boundary objects. (DEM)
  if(UseDEM){
    if(UseDEM){
      DemData=new StDemData[DemDataSize];
      memset(DemData,0,sizeof(StDemData)*DemDataSize);
    }
    for(unsigned c=0;c<parts.CountBlocks();c++){
      const JSpacePartBlock &block=parts.GetBlock(c);
      if(IsBound(block.Type)){
        const word mkbound=block.GetMkType();
        const unsigned cmk=MkInfo->GetMkBlockByMkBound(mkbound);
        if(cmk>=MkInfo->Size())RunException(met,fun::PrintStr("Error loading boundary objects. Mkbound=%u is unknown.",mkbound));
        const StDemData data=LoadDemData(UseDEM,UseDEM,&block);
        if(UseDEM){
          const unsigned tav=CODE_GetTypeAndValue(MkInfo->Mkblock(cmk)->Code); //:Log->Printf("___> tav[%u]:%u",cmk,tav);
          DemData[tav]=data;
        }
      }
    }
  }

  NpMinimum=CaseNp-unsigned(PartsOutMax*CaseNfluid);
  Log->Print("**Basic case configuration is loaded");
}

//==============================================================================
/// Loads coefficients used for DEM or Chrono objects.
//==============================================================================
StDemData JSph::LoadDemData(bool basicdata,bool extradata,const JSpacePartBlock* block)const{
  const char met[]="LoadDemData";
  const word mk=block->GetMk();
  const word mkbound=block->GetMkType();
  StDemData data;
  memset(&data,0,sizeof(StDemData));
  if(basicdata){
    if(!block->ExistsSubValue("Kfric","value"))RunException(met,fun::PrintStr("Object mk=%u (mkbound=%u) - Value of Kfric is invalid.",mk,mkbound));
    if(!block->ExistsSubValue("Restitution_Coefficient","value"))RunException(met,fun::PrintStr("Object mk=%u (mkbound=%u) - Value of Restitution_Coefficient is invalid.",mk,mkbound));
    data.kfric=block->GetSubValueFloat("Kfric","value",true,0);
    data.restitu=block->GetSubValueFloat("Restitution_Coefficient","value",true,0);
    if(block->ExistsValue("Restitution_Coefficient_User"))data.restitu=block->GetValueFloat("Restitution_Coefficient_User");
  }
  if(extradata){
    data.massp=MassBound;
    if(block->Type==TpPartFloating){
      const JSpacePartBlock_Floating *fblock=(const JSpacePartBlock_Floating *)block;
      data.mass=(float)fblock->GetMassbody();
      data.massp=(float)(fblock->GetMassbody()/fblock->GetCount());
    }
    if(!block->ExistsSubValue("Young_Modulus","value"))RunException(met,fun::PrintStr("Object mk=%u (mkbound=%u) - Value of Young_Modulus is invalid.",mk,mkbound));
    if(!block->ExistsSubValue("PoissonRatio","value"))RunException(met,fun::PrintStr("Object mk=%u (mkbound=%u) - Value of PoissonRatio is invalid.",mk,mkbound));
    data.young=block->GetSubValueFloat("Young_Modulus","value",true,0);
    data.poisson=block->GetSubValueFloat("PoissonRatio","value",true,0);
    data.tau=(data.young? (1-data.poisson*data.poisson)/data.young: 0);
  }
  return(data);
}

//==============================================================================
/// Shows coefficients used for DEM objects.
//==============================================================================
void JSph::VisuDemCoefficients()const{
  //-Gets info for each block of particles.
  Log->Printf("Coefficients for DEM:");
  for(unsigned c=0;c<MkInfo->Size();c++){
    const JSphMkBlock *pmk=MkInfo->Mkblock(c);
    const typecode code=pmk->Code;
    const typecode type=CODE_GetType(code);
    const unsigned tav=CODE_GetTypeAndValue(code);
    if(type==CODE_TYPE_FIXED || type==CODE_TYPE_MOVING || type==CODE_TYPE_FLOATING){
      Log->Printf("  Object %s  mkbound:%u  mk:%u",(type==CODE_TYPE_FIXED? "Fixed": (type==CODE_TYPE_MOVING? "Moving": "Floating")),pmk->MkType,pmk->Mk);
      Log->Printf("    Young_Modulus: %g",DemData[tav].young);
      Log->Printf("    PoissonRatio.: %g",DemData[tav].poisson);
      Log->Printf("    Kfric........: %g",DemData[tav].kfric);
      Log->Printf("    Restitution..: %g",DemData[tav].restitu);
    }
  }
}

//==============================================================================
/// Loads the code of a particle group and flags the last "nout" 
/// particles as excluded. 
///
/// Carga el codigo de grupo de las particulas y marca las nout ultimas
/// particulas como excluidas.
//==============================================================================
void JSph::LoadCodeParticles(unsigned np,const unsigned *idp,typecode *code)const{
  const char met[]="LoadCodeParticles"; 
  //-Assigns code to each group of particles.
  for(unsigned p=0;p<np;p++)code[p]=MkInfo->GetCodeById(idp[p]);
}

//==============================================================================
/// Sets DBL_MAX values by indicated values.
//==============================================================================
void JSph::PrepareCfgDomainValues(tdouble3 &v,tdouble3 vdef)const{
  if(v.x==DBL_MAX)v.x=vdef.x;
  if(v.y==DBL_MAX)v.y=vdef.y;
  if(v.z==DBL_MAX)v.z=vdef.z;
}

//==============================================================================
/// Resizes limits of the map according to case configuration.
//==============================================================================
void JSph::ResizeMapLimits(){
  Log->Print(string("MapRealPos(border)=")+fun::Double3gRangeStr(MapRealPosMin,MapRealPosMax));
  tdouble3 dmin=MapRealPosMin,dmax=MapRealPosMax;
  //-Sets Y configuration when it is a 2-D simulation.
  if(Simulate2D){
    CfgDomainParticlesMin.y=CfgDomainParticlesMax.y=DBL_MAX;
    CfgDomainParticlesPrcMin.y=CfgDomainParticlesPrcMax.y=DBL_MAX;
    CfgDomainFixedMin.y=CfgDomainFixedMax.y=DBL_MAX;
  }
  //-Configuration according particles domain.
  PrepareCfgDomainValues(CfgDomainParticlesMin);
  PrepareCfgDomainValues(CfgDomainParticlesMax);
  PrepareCfgDomainValues(CfgDomainParticlesPrcMin);
  PrepareCfgDomainValues(CfgDomainParticlesPrcMax);
  const tdouble3 dif=dmax-dmin;
  dmin=dmin-dif*CfgDomainParticlesPrcMin;
  dmax=dmax+dif*CfgDomainParticlesPrcMax;
  dmin=dmin-CfgDomainParticlesMin;
  dmax=dmax+CfgDomainParticlesMax;
  //-Fixed domain configuration.
  PrepareCfgDomainValues(CfgDomainFixedMin,dmin);
  PrepareCfgDomainValues(CfgDomainFixedMax,dmax);
  dmin=CfgDomainFixedMin; 
  dmax=CfgDomainFixedMax; 
  //-Checks domain limits.
  if(dmin.x>MapRealPosMin.x||dmin.y>MapRealPosMin.y||dmin.z>MapRealPosMin.z||dmax.x<MapRealPosMax.x||dmax.y<MapRealPosMax.y||dmax.z<MapRealPosMax.z)
    RunException("ResizeMapLimits",fun::PrintStr("Domain limits %s are not valid.",fun::Double3gRangeStr(dmin,dmax).c_str()));
  //-Periodic domain configuration.
  if(!PeriX){ MapRealPosMin.x=dmin.x; MapRealPosMax.x=dmax.x; }
  if(!PeriY){ MapRealPosMin.y=dmin.y; MapRealPosMax.y=dmax.y; }
  if(!PeriZ){ MapRealPosMin.z=dmin.z; MapRealPosMax.z=dmax.z; }
}

//==============================================================================
/// Configures value of constants.
//==============================================================================
void JSph::ConfigConstants(bool simulate2d){
  const char* met="ConfigConstants";
  //-Computation of constants.
  const double h=H;
  Delta2H=float(h*2*DeltaSph);
  Cs0=sqrt(double(Gamma)*double(CteB)/double(RhopZero));
  if(!DtIni)DtIni=h/Cs0;
  if(!DtMin)DtMin=(h/Cs0)*CoefDtMin;
  Dosh=float(h*2); 
  H2=float(h*h);
  Fourh2=float(h*h*4); 
  Eta2=float((h*0.1)*(h*0.1));
  if(simulate2d){
    if(TKernel==KERNEL_Wendland){
      Awen=float(0.557/(h*h));
      Bwen=float(-2.7852/(h*h*h));
    }
    else if(TKernel==KERNEL_Gaussian){
      const double a1=4./PI;
      const double a2=a1/(h*h);
      const double aa=a1/(h*h*h);
      Agau=float(a2);
      Bgau=float(-8.*aa);
    }
    else if(TKernel==KERNEL_Cubic){
      const double a1=10./(PI*7.);
      const double a2=a1/(h*h);
      const double aa=a1/(h*h*h);
      const double deltap=1./1.5;
      const double wdeltap=a2*(1.-1.5*deltap*deltap+0.75*deltap*deltap*deltap);
      CubicCte.od_wdeltap=float(1./wdeltap);
      CubicCte.a1=float(a1);
      CubicCte.a2=float(a2);
      CubicCte.aa=float(aa);
      CubicCte.a24=float(0.25*a2);
      CubicCte.c1=float(-3.*aa);
      CubicCte.d1=float(9.*aa/4.);
      CubicCte.c2=float(-3.*aa/4.);
    }
  }
  else{
    if(TKernel==KERNEL_Wendland){
      Awen=float(0.41778/(h*h*h));
      Bwen=float(-2.08891/(h*h*h*h));
    }
    else if(TKernel==KERNEL_Gaussian){
      const double a1=8./5.5683;
      const double a2=a1/(h*h*h);
      const double aa=a1/(h*h*h*h); 
      Agau=float(a2);
      Bgau=float(-8.*aa);
    }
    else if(TKernel==KERNEL_Cubic){
      const double a1=1./PI;
      const double a2=a1/(h*h*h);
      const double aa=a1/(h*h*h*h);
      const double deltap=1./1.5;
      const double wdeltap=a2*(1.-1.5*deltap*deltap+0.75*deltap*deltap*deltap);
      CubicCte.od_wdeltap=float(1./wdeltap);
      CubicCte.a1=float(a1);
      CubicCte.a2=float(a2);
      CubicCte.aa=float(aa);
      CubicCte.a24=float(0.25*a2);
      CubicCte.c1=float(-3.*aa);
      CubicCte.d1=float(9.*aa/4.);
      CubicCte.c2=float(-3.*aa/4.);
    }
  }
  //-Constants for Laminar viscosity + SPS turbulence model.
  if(TVisco==VISCO_LaminarSPS){  
    double dp_sps=(Simulate2D? sqrt(Dp*Dp*2.)/2.: sqrt(Dp*Dp*3.)/3.);  
    SpsSmag=float(pow((0.12*dp_sps),2));
    SpsBlin=float((2./3.)*0.0066*dp_sps*dp_sps); 
  }
  VisuConfig();
}

//==============================================================================
/// Prints out configuration of the case.
//==============================================================================
void JSph::VisuConfig()const{
  const char* met="VisuConfig";
  Log->Print(Simulate2D? "**2D-Simulation parameters:": "**3D-Simulation parameters:");
  Log->Print(fun::VarStr("CaseName",CaseName));
  Log->Print(fun::VarStr("RunName",RunName));
  if(Simulate2D)Log->Print(fun::VarStr("Simulate2DPosY",Simulate2DPosY));
  Log->Print(fun::VarStr("PosDouble",GetPosDoubleName(Psingle,SvDouble)));
  Log->Print(fun::VarStr("SvTimers",SvTimers));
  Log->Print(fun::VarStr("StepAlgorithm",GetStepName(TStep)));
  if(TStep==STEP_None)RunException(met,"StepAlgorithm value is invalid.");
  if(TStep==STEP_Verlet)Log->Print(fun::VarStr("VerletSteps",VerletSteps));
  Log->Print(fun::VarStr("Kernel",GetKernelName(TKernel)));
  Log->Print(fun::VarStr("Viscosity",GetViscoName(TVisco)));
  Log->Print(fun::VarStr("Visco",Visco));
  Log->Print(fun::VarStr("ViscoBoundFactor",ViscoBoundFactor));
  if(ViscoTime)Log->Print(fun::VarStr("ViscoTime",ViscoTime->GetFile()));
  Log->Print(fun::VarStr("DeltaSph",GetDeltaSphName(TDeltaSph)));
  if(TDeltaSph!=DELTA_None)Log->Print(fun::VarStr("DeltaSphValue",DeltaSph));
  Log->Print(fun::VarStr("Shifting",GetShiftingName(TShifting)));
  if(TShifting!=SHIFT_None){
    Log->Print(fun::VarStr("ShiftCoef",ShiftCoef));
    if(ShiftTFS)Log->Print(fun::VarStr("ShiftTFS",ShiftTFS));
  }
  string rigidalgorithm=(!FtCount? "None": (UseDEM? "SPH+DCDEM": "SPH"));
  Log->Print(fun::VarStr("RigidAlgorithm",rigidalgorithm));
  Log->Print(fun::VarStr("FloatingCount",FtCount));
  if(FtCount)Log->Print(fun::VarStr("FtPause",FtPause));
  Log->Print(fun::VarStr("CaseNp",CaseNp));
  Log->Print(fun::VarStr("CaseNbound",CaseNbound));
  Log->Print(fun::VarStr("CaseNfixed",CaseNfixed));
  Log->Print(fun::VarStr("CaseNmoving",CaseNmoving));
  Log->Print(fun::VarStr("CaseNfloat",CaseNfloat));
  Log->Print(fun::VarStr("CaseNfluid",CaseNfluid));
  Log->Print(fun::VarStr("PeriodicActive",TpPeriName(TpPeri(PeriActive))));
  if(PeriX)Log->Print(fun::VarStr("PeriodicXinc",PeriXinc));
  if(PeriY)Log->Print(fun::VarStr("PeriodicYinc",PeriYinc));
  if(PeriZ)Log->Print(fun::VarStr("PeriodicZinc",PeriZinc));
  Log->Print(fun::VarStr("Dx",Dp));
  Log->Print(fun::VarStr("H",H));
  Log->Print(fun::VarStr("CoefficientH",H/(Dp*sqrt(Simulate2D? 2.f: 3.f))));
  Log->Print(fun::VarStr("CteB",CteB));
  Log->Print(fun::VarStr("Gamma",Gamma));
  Log->Print(fun::VarStr("RhopZero",RhopZero));
  Log->Print(fun::VarStr("Cs0",Cs0));
  Log->Print(fun::VarStr("CFLnumber",CFLnumber));
  Log->Print(fun::VarStr("DtIni",DtIni));
  Log->Print(fun::VarStr("DtMin",DtMin));
  Log->Print(fun::VarStr("DtAllParticles",DtAllParticles));
  if(DtFixed)Log->Print(fun::VarStr("DtFixed",DtFixed->GetFile()));
  Log->Print(fun::VarStr("MassFluid",MassFluid));
  Log->Print(fun::VarStr("MassBound",MassBound));
  if(TKernel==KERNEL_Wendland){
    Log->Print(fun::VarStr("Awen (Wendland)",Awen));
    Log->Print(fun::VarStr("Bwen (Wendland)",Bwen));
  }
  else if(TKernel==KERNEL_Gaussian){
    Log->Print(fun::VarStr("Agau (Gaussian)",Agau));
    Log->Print(fun::VarStr("Bgau (Gaussian)",Bgau));
  }
  else if(TKernel==KERNEL_Cubic){
    Log->Print(fun::VarStr("CubicCte.a1",CubicCte.a1));
    Log->Print(fun::VarStr("CubicCte.aa",CubicCte.aa));
    Log->Print(fun::VarStr("CubicCte.a24",CubicCte.a24));
    Log->Print(fun::VarStr("CubicCte.c1",CubicCte.c1));
    Log->Print(fun::VarStr("CubicCte.c2",CubicCte.c2));
    Log->Print(fun::VarStr("CubicCte.d1",CubicCte.d1));
    Log->Print(fun::VarStr("CubicCte.od_wdeltap",CubicCte.od_wdeltap));
  }
  if(TVisco==VISCO_LaminarSPS){     
    Log->Print(fun::VarStr("SpsSmag",SpsSmag));
    Log->Print(fun::VarStr("SpsBlin",SpsBlin));
  }
  if(UseDEM)VisuDemCoefficients();
  if(CaseNfloat)Log->Print(fun::VarStr("FtPause",FtPause));
  Log->Print(fun::VarStr("TimeMax",TimeMax));
  Log->Print(fun::VarStr("TimePart",TimePart));
  Log->Print(fun::VarStr("Gravity",Gravity));
  Log->Print(fun::VarStr("NpMinimum",NpMinimum));
  Log->Print(fun::VarStr("RhopOut",RhopOut));
  if(RhopOut){
    Log->Print(fun::VarStr("RhopOutMin",RhopOutMin));
    Log->Print(fun::VarStr("RhopOutMax",RhopOutMax));
  }
  if(CteB==0)RunException(met,"Constant \'b\' cannot be zero.\n\'b\' is zero when fluid height is zero (or fluid particles were not created)");
}

//==============================================================================
/// Shows particle and MK blocks summary.
//==============================================================================
void JSph::VisuParticleSummary()const{
  JXml xml; xml.LoadFile(FileXml);
  JSpaceParts parts; 
  parts.LoadXml(&xml,"case.execution.particles");
  std::vector<std::string> summary;
  parts.GetParticleSummary(summary);
  Log->Print(summary);
  Log->Print(" ");
}

//==============================================================================
/// Computes cell particles and checks if there are more particles
/// excluded than expected.
///
/// Calcula celda de las particulas y comprueba que no existan mas particulas
/// excluidas de las previstas.
//==============================================================================
void JSph::LoadDcellParticles(unsigned n,const typecode *code,const tdouble3 *pos,unsigned *dcell)const{
  const char met[]="LoadDcellParticles";
  for(unsigned p=0;p<n;p++){
    typecode codeout=CODE_GetSpecialValue(code[p]);
    if(codeout<CODE_OUTIGNORE){
      const tdouble3 ps=pos[p];
      if(ps>=DomRealPosMin && ps<DomRealPosMax){//-Particle in.
        const double dx=ps.x-DomPosMin.x;
        const double dy=ps.y-DomPosMin.y;
        const double dz=ps.z-DomPosMin.z;
        unsigned cx=unsigned(dx/Scell),cy=unsigned(dy/Scell),cz=unsigned(dz/Scell);
        dcell[p]=PC__Cell(DomCellCode,cx,cy,cz);
      }
      else{//-Particle out.
        RunException(met,"Found new particles out."); //-There can not be new particles excluded. | No puede haber nuevas particulas excluidas.
        dcell[p]=PC__CodeOut;
      }
    }
    else dcell[p]=PC__CodeOut;
  }
}

//==============================================================================
/// Initializes data of particles according XML configuration.
///
/// Inicializa datos de las particulas a partir de la configuracion en el XML.
//==============================================================================
void JSph::RunInitialize(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp,const typecode *code,tfloat4 *velrhop){
  const char met[]="RunInitialize";
  if(!PartBegin){
    JSphInitialize init(FileXml);
    if(init.Count()){
      //-Creates array with mktype value.
      word *mktype=new word[np];
      for(unsigned p=0;p<np;p++){
        const unsigned cmk=MkInfo->GetMkBlockByCode(code[p]);
        mktype[p]=(cmk<MkInfo->Size()? word(MkInfo->Mkblock(cmk)->MkType): USHRT_MAX);
      }
      init.Run(np,npb,pos,idp,mktype,velrhop);
      init.GetConfig(InitializeInfo);
      //-Frees memory.
      delete[] mktype; mktype=NULL;
    }
  }
}

//==============================================================================
/// Creates PartsInit object with initial particle data for automatic 
/// configurations.
///
/// Crea el objeto PartsInit con los datos iniciales de las particulas para 
/// configuraciones automaticas.
//==============================================================================
void JSph::CreatePartsInit(unsigned np,const tdouble3 *pos,const typecode *code){
  PartsInit=new JSphPartsInit(Simulate2D,Simulate2DPosY,Dp,MkInfo,np,pos,code);
}

//==============================================================================
/// Free memory of PartsInit.
///
/// Libera memoria de PartsInit.
//==============================================================================
void JSph::FreePartsInit(){
  delete PartsInit; PartsInit=NULL;
}

//==============================================================================
/// Configures cell division.
//==============================================================================
void JSph::ConfigCellDivision(){
  if(CellMode!=CELLMODE_2H && CellMode!=CELLMODE_H)RunException("ConfigCellDivision","The CellMode is invalid.");
  Hdiv=(CellMode==CELLMODE_2H? 1: 2);
  Scell=Dosh/Hdiv;
  MovLimit=Scell*0.9f;
  Map_Cells=TUint3(unsigned(ceil(Map_Size.x/Scell)),unsigned(ceil(Map_Size.y/Scell)),unsigned(ceil(Map_Size.z/Scell)));
  //-Prints configuration.
  Log->Print(fun::VarStr("CellMode",string(GetNameCellMode(CellMode))));
  Log->Print(fun::VarStr("Hdiv",Hdiv));
  Log->Print(string("MapCells=(")+fun::Uint3Str(Map_Cells)+")");
  //-Creates VTK file with map cells.
  if(SaveMapCellsVtkSize()<1024*1024*10)SaveMapCellsVtk(Scell);
  else Log->PrintWarning("File CfgInit_MapCells.vtk was not created because number of cells is too high.");
}

//==============================================================================
/// Sets local domain of simulation within Map_Cells and computes DomCellCode.
/// Establece dominio local de simulacion dentro de Map_Cells y calcula DomCellCode.
//==============================================================================
void JSph::SelecDomain(tuint3 celini,tuint3 celfin){
  const char met[]="SelecDomain";
  DomCelIni=celini;
  DomCelFin=celfin;
  DomCells=DomCelFin-DomCelIni;
  if(DomCelIni.x>=Map_Cells.x || DomCelIni.y>=Map_Cells.y || DomCelIni.z>=Map_Cells.z )RunException(met,"DomCelIni is invalid.");
  if(DomCelFin.x>Map_Cells.x || DomCelFin.y>Map_Cells.y || DomCelFin.z>Map_Cells.z )RunException(met,"DomCelFin is invalid.");
  if(DomCells.x<1 || DomCells.y<1 || DomCells.z<1 )RunException(met,"The domain of cells is invalid.");
  //-Computes local domain limits.
  DomPosMin.x=Map_PosMin.x+(DomCelIni.x*Scell);
  DomPosMin.y=Map_PosMin.y+(DomCelIni.y*Scell);
  DomPosMin.z=Map_PosMin.z+(DomCelIni.z*Scell);
  DomPosMax.x=Map_PosMin.x+(DomCelFin.x*Scell);
  DomPosMax.y=Map_PosMin.y+(DomCelFin.y*Scell);
  DomPosMax.z=Map_PosMin.z+(DomCelFin.z*Scell);
  //-Adjusts final limits.
  if(DomPosMax.x>Map_PosMax.x)DomPosMax.x=Map_PosMax.x;
  if(DomPosMax.y>Map_PosMax.y)DomPosMax.y=Map_PosMax.y;
  if(DomPosMax.z>Map_PosMax.z)DomPosMax.z=Map_PosMax.z;
  //-Computes actual limits of local domain.
  //-Calcula limites reales del dominio local.
  DomRealPosMin=DomPosMin;
  DomRealPosMax=DomPosMax;
  if(DomRealPosMax.x>MapRealPosMax.x)DomRealPosMax.x=MapRealPosMax.x;
  if(DomRealPosMax.y>MapRealPosMax.y)DomRealPosMax.y=MapRealPosMax.y;
  if(DomRealPosMax.z>MapRealPosMax.z)DomRealPosMax.z=MapRealPosMax.z;
  if(DomRealPosMin.x<MapRealPosMin.x)DomRealPosMin.x=MapRealPosMin.x;
  if(DomRealPosMin.y<MapRealPosMin.y)DomRealPosMin.y=MapRealPosMin.y;
  if(DomRealPosMin.z<MapRealPosMin.z)DomRealPosMin.z=MapRealPosMin.z;
  //-Computes cofification of cells for the selected domain.
  //-Calcula codificacion de celdas para el dominio seleccionado.
  DomCellCode=CalcCellCode(DomCells+TUint3(1));
  if(!DomCellCode)RunException(met,string("Failed to select a valid CellCode for ")+fun::UintStr(DomCells.x)+"x"+fun::UintStr(DomCells.y)+"x"+fun::UintStr(DomCells.z)+" cells (CellMode="+GetNameCellMode(CellMode)+").");
  //-Prints configurantion.
  Log->Print(string("DomCells=(")+fun::Uint3Str(DomCells)+")");
  Log->Print(fun::VarStr("DomCellCode",fun::UintStr(PC__GetSx(DomCellCode))+"_"+fun::UintStr(PC__GetSy(DomCellCode))+"_"+fun::UintStr(PC__GetSz(DomCellCode))));
}

//==============================================================================
/// Selects an adequate code for cell configuration.
/// Selecciona un codigo adecuado para la codificion de celda.
//==============================================================================
unsigned JSph::CalcCellCode(tuint3 ncells){
  unsigned sxmin=2; for(;ncells.x>>sxmin;sxmin++);
  unsigned symin=2; for(;ncells.y>>symin;symin++);
  unsigned szmin=2; for(;ncells.z>>szmin;szmin++);
  unsigned smin=sxmin+symin+szmin;
  unsigned ccode=0;
  if(smin<=32){
    unsigned sx=sxmin,sy=symin,sz=szmin;
    unsigned rest=32-smin;
    while(rest){
      if(rest){ sx++; rest--; }
      if(rest){ sy++; rest--; }
      if(rest){ sz++; rest--; }
    }
    ccode=PC__GetCode(sx,sy,sz);
  }
  return(ccode);
}

//==============================================================================
/// Computes maximum distance between particles and center of floating.
/// Calcula distancia maxima entre particulas y centro de cada floating.
//==============================================================================
void JSph::CalcFloatingRadius(unsigned np,const tdouble3 *pos,const unsigned *idp){
  const char met[]="CalcFloatingsRadius";
  const float overradius=1.2f; //-Percentage of ration increase. | Porcentaje de incremento de radio. 
  unsigned *ridp=new unsigned[CaseNfloat];
  //-Assigns values UINT_MAX. 
  memset(ridp,255,sizeof(unsigned)*CaseNfloat); 
  //-Computes position according to id assuming that all particles are not periodic.
  //-Calcula posicion segun id suponiendo que todas las particulas son normales (no periodicas).
  const unsigned idini=CaseNpb,idfin=CaseNpb+CaseNfloat;
  for(unsigned p=0;p<np;p++){
    const unsigned id=idp[p];
    if(idini<=id && id<idfin)ridp[id-idini]=p;
  }
  //-Checks that all floating particles are located.  
  //-Comprueba que todas las particulas floating estan localizadas.
  for(unsigned fp=0;fp<CaseNfloat;fp++){
    if(ridp[fp]==UINT_MAX)RunException(met,"There are floating particles not found.");
  }
  //-Calculates maximum distance between particles and center of the floating (all are valid).  
  //-Calcula distancia maxima entre particulas y centro de floating (todas son validas).
  float radiusmax=0;
  for(unsigned cf=0;cf<FtCount;cf++){
    StFloatingData *fobj=FtObjs+cf;
    const unsigned fpini=fobj->begin-CaseNpb;
    const unsigned fpfin=fpini+fobj->count;
    const tdouble3 fcen=fobj->center;
    double r2max=0;
    for(unsigned fp=fpini;fp<fpfin;fp++){
      const int p=ridp[fp];
      const double dx=fcen.x-pos[p].x,dy=fcen.y-pos[p].y,dz=fcen.z-pos[p].z;
      double r2=dx*dx+dy*dy+dz*dz;
      if(r2max<r2)r2max=r2;
    }
    fobj->radius=float(sqrt(r2max)*overradius);
    if(radiusmax<fobj->radius)radiusmax=fobj->radius;
  }
  //-Deallocate of memory.  
  delete[] ridp; ridp=NULL;
  //-Checks maximum radius < dimensions of the periodic domain.
  //-Comprueba que el radio maximo sea menor que las dimensiones del dominio periodico.
  if(PeriX && fabs(PeriXinc.x)<=radiusmax)RunException(met,fun::PrintStr("The floating radius (%g) is too large for periodic distance in X (%g).",radiusmax,abs(PeriXinc.x)));
  if(PeriY && fabs(PeriYinc.y)<=radiusmax)RunException(met,fun::PrintStr("The floating radius (%g) is too large for periodic distance in Y (%g).",radiusmax,abs(PeriYinc.y)));
  if(PeriZ && fabs(PeriZinc.z)<=radiusmax)RunException(met,fun::PrintStr("The floating radius (%g) is too large for periodic distance in Z (%g).",radiusmax,abs(PeriZinc.z)));
}

//==============================================================================
/// Returns the corrected position after applying periodic conditions.
/// Devuelve la posicion corregida tras aplicar condiciones periodicas.
//==============================================================================
tdouble3 JSph::UpdatePeriodicPos(tdouble3 ps)const{
  double dx=ps.x-MapRealPosMin.x;
  double dy=ps.y-MapRealPosMin.y;
  double dz=ps.z-MapRealPosMin.z;
  const bool out=(dx!=dx || dy!=dy || dz!=dz || dx<0 || dy<0 || dz<0 || dx>=MapRealSize.x || dy>=MapRealSize.y || dz>=MapRealSize.z);
  //-Adjusts position according to periodic conditions and checks again domain limtis.
  //-Ajusta posicion segun condiciones periodicas y vuelve a comprobar los limites del dominio.
  if(PeriActive && out){
    if(PeriX){
      if(dx<0)             { dx-=PeriXinc.x; dy-=PeriXinc.y; dz-=PeriXinc.z; }
      if(dx>=MapRealSize.x){ dx+=PeriXinc.x; dy+=PeriXinc.y; dz+=PeriXinc.z; }
    }
    if(PeriY){
      if(dy<0)             { dx-=PeriYinc.x; dy-=PeriYinc.y; dz-=PeriYinc.z; }
      if(dy>=MapRealSize.y){ dx+=PeriYinc.x; dy+=PeriYinc.y; dz+=PeriYinc.z; }
    }
    if(PeriZ){
      if(dz<0)             { dx-=PeriZinc.x; dy-=PeriZinc.y; dz-=PeriZinc.z; }
      if(dz>=MapRealSize.z){ dx+=PeriZinc.x; dy+=PeriZinc.y; dz+=PeriZinc.z; }
    }
    ps=TDouble3(dx,dy,dz)+MapRealPosMin;
  }
  return(ps);
}

//==============================================================================
/// Sets configuration for recordering of particles.
/// Establece configuracion para grabacion de particulas.
//==============================================================================
void JSph::RestartCheckData(){
  const char met[]="RestartCheckData";
  if(PartBegin){
    //-Loads particle blocks information.
    JPartDataHead parthead;
    parthead.LoadFile(PartBeginDir);
    const string filehead=fun::GetDirWithSlash(PartBeginDir)+parthead.GetFileName();
    //-Checks particle blocks information.
    const unsigned nmk=MkInfo->Size();
    if(nmk!=parthead.MkBlockCount())RunException(met,"Number of Mk blocks is invalid.",filehead);
    for(unsigned c=0;c<nmk;c++){
      const JSphMkBlock* pmk=MkInfo->Mkblock(c);
      const JPartDataHeadMkBlock& mbk=parthead.Mkblock(c);
      if(pmk->Type!=mbk.Type)RunException(met,fun::PrintStr("Type of Mk block %u does not match.",c),filehead);
      if(pmk->Mk!=mbk.Mk)RunException(met,fun::PrintStr("Mk value of Mk block %u does not match.",c),filehead);
      if(pmk->MkType!=mbk.MkType)RunException(met,fun::PrintStr("MkType value of Mk block %u does not match.",c),filehead);
      if(pmk->Count!=mbk.Count)RunException(met,fun::PrintStr("Count value of Mk block %u does not match.",c),filehead);
    }
  }
}

//==============================================================================
/// Load particles of case and process.
/// Carga particulas del caso a procesar.
//==============================================================================
void JSph::LoadCaseParticles(){
  Log->Print("Loading initial state of particles...");
  PartsLoaded=new JPartsLoad4(Cpu);
  PartsLoaded->LoadParticles(DirCase,CaseName,PartBegin,PartBeginDir);
  PartsLoaded->CheckConfig(CaseNp,CaseNfixed,CaseNmoving,CaseNfloat,CaseNfluid,TpPeri(PeriActive));
  if(PartBegin)RestartCheckData();
  Log->Printf("Loaded particles: %u",PartsLoaded->GetCount());

  //-Collect information of loaded particles.
  //-Recupera informacion de las particulas cargadas.
  Simulate2D=PartsLoaded->GetSimulate2D();
  Simulate2DPosY=PartsLoaded->GetSimulate2DPosY();
  if(Simulate2D && PeriY)RunException("LoadCaseParticles","Cannot use periodic conditions in Y with 2D simulations");
  CasePosMin=PartsLoaded->GetCasePosMin();
  CasePosMax=PartsLoaded->GetCasePosMax();

  //-Computes actual limits of simulation.
  //-Calcula limites reales de la simulacion.
  if(PartsLoaded->MapSizeLoaded())PartsLoaded->GetMapSize(MapRealPosMin,MapRealPosMax);
  else{
    PartsLoaded->CalculeLimits(double(H)*BORDER_MAP,Dp/2.,PeriX,PeriY,PeriZ,MapRealPosMin,MapRealPosMax);
    ResizeMapLimits();
  }
  Log->Print(string("MapRealPos(final)=")+fun::Double3gRangeStr(MapRealPosMin,MapRealPosMax));
  MapRealSize=MapRealPosMax-MapRealPosMin;
  Log->Print("**Initial state of particles is loaded");

  //-Configure limits of periodic axes. 
  //-Configura limites de ejes periodicos.
  if(PeriX)PeriXinc.x=-MapRealSize.x;
  if(PeriY)PeriYinc.y=-MapRealSize.y;
  if(PeriZ)PeriZinc.z=-MapRealSize.z;
  //-Computes simulation limits with periodic boundaries.
  //-Calcula limites de simulacion con bordes periodicos.
  Map_PosMin=MapRealPosMin; Map_PosMax=MapRealPosMax;
  float dosh=float(H*2);
  if(PeriX){ Map_PosMin.x=Map_PosMin.x-dosh;  Map_PosMax.x=Map_PosMax.x+dosh; }
  if(PeriY){ Map_PosMin.y=Map_PosMin.y-dosh;  Map_PosMax.y=Map_PosMax.y+dosh; }
  if(PeriZ){ Map_PosMin.z=Map_PosMin.z-dosh;  Map_PosMax.z=Map_PosMax.z+dosh; }
  Map_Size=Map_PosMax-Map_PosMin;
  //-Saves initial domain in a VTK file (CasePosMin/Max, MapRealPosMin/Max and Map_PosMin/Max).
  SaveInitialDomainVtk();
}

//==============================================================================
/// Initialisation of variables and objects for execution.
/// Inicializa variables y objetos para la ejecucion.
//==============================================================================
void JSph::InitRun(unsigned np,const unsigned *idp,const tdouble3 *pos){
  const char met[]="InitRun";
  VerletStep=0;
  if(TStep==STEP_Symplectic)SymplecticDtPre=DtIni;
  if(UseDEM)DemDtForce=DtIni; //(DEM)
  
  //-Loads data from other simulation to restart it.
  if(PartBegin){
    PartBeginTimeStep=PartsLoaded->GetPartBeginTimeStep();
    PartBeginTotalNp=PartsLoaded->GetPartBeginTotalNp();
    if(TStep==STEP_Symplectic && PartsLoaded->GetSymplecticDtPre())SymplecticDtPre=PartsLoaded->GetSymplecticDtPre();
    if(UseDEM  && PartsLoaded->GetDemDtForce())DemDtForce=PartsLoaded->GetDemDtForce(); //(DEM)
  }
  //-Free memory of PartsLoaded.
  delete PartsLoaded; PartsLoaded=NULL;

  //-Adjust paramaters to start.
  PartIni=PartBeginFirst;
  TimeStepIni=(!PartIni? 0: PartBeginTimeStep);
  //-Adjust motion for the instant of the loaded PART.
  if(SphMotion){
    SphMotion->SetTimeMod(!PartIni? PartBeginTimeStep: 0);
    SphMotion->ProcesTime(JSphMotion::MOMT_Simple,0,TimeStepIni);
  }
  //-Adjust motion paddles for the instant of the loaded PART.
  if(WaveGen)WaveGen->SetTimeMod(!PartIni? PartBeginTimeStep: 0);

  //-Uses Inlet information from PART read.
  if(PartBeginTimeStep && PartBeginTotalNp){
    TotalNp=PartBeginTotalNp;
    IdMax=unsigned(TotalNp-1);
  }

  //-Shows Initialize configuration.
  if(InitializeInfo.size()){
    Log->Print("Initialization configuration:");
    Log->Print(InitializeInfo);
    Log->Print(" ");
  }

  //-Process Special configurations in XML.
  JXml xml; xml.LoadFile(FileXml);

  //-Configuration of GaugeSystem.
  GaugeSystem->Config(Simulate2D,Simulate2DPosY,TimeMax,TimePart,Dp,DomPosMin,DomPosMax,Scell,Hdiv,H,MassFluid,MassBound,CteB,Gamma,RhopZero);
  if(xml.GetNode("case.execution.special.gauges",false))GaugeSystem->LoadXml(&xml,"case.execution.special.gauges",MkInfo);

  //-Prepares WaveGen configuration.
  if(WaveGen){
    Log->Print("Wave paddles configuration:");
    WaveGen->Init(GaugeSystem,MkInfo,TimeMax,Gravity);
    WaveGen->VisuConfig(""," ");
  }

  //-Prepares Damping configuration.
  if(Damping){
    Damping->VisuConfig("Damping configuration:"," ");
  }

  //-Prepares AccInput configuration.
  if(AccInput){
    Log->Print("AccInput configuration:");
    AccInput->Init(TimeMax);
    AccInput->VisuConfig(""," ");
  }

  //-Configuration of SaveDt.
  if(xml.GetNode("case.execution.special.savedt",false)){
    SaveDt=new JSaveDt(Log);
    SaveDt->Config(&xml,"case.execution.special.savedt",TimeMax,TimePart);
    SaveDt->VisuConfig("SaveDt configuration:"," ");
  }

  //-Shows configuration of JGaugeSystem.
  if(GaugeSystem->GetCount())GaugeSystem->VisuConfig("GaugeSystem configuration:"," ");

  //-Shows configuration of JTimeOut.
  if(TimeOut->UseSpecialConfig())TimeOut->VisuConfig(Log,"TimeOut configuration:"," ");

  Part=PartIni; Nstep=0; PartNstep=0; PartOut=0;
  TimeStep=TimeStepIni; TimeStepM1=TimeStep;
  if(DtFixed)DtIni=DtFixed->GetDt(TimeStep,DtIni);
  TimePartNext=TimeOut->GetNextTime(TimeStep);
}

//==============================================================================
/// Display a message with reserved memory for the basic data of particles.
/// Muestra un mensaje con la memoria reservada para los datos basicos de las particulas.
//==============================================================================
void JSph::PrintSizeNp(unsigned np,llong size)const{
  Log->Printf("**Requested %s memory for %u particles: %.1f MB.",(Cpu? "cpu": "gpu"),np,double(size)/(1024*1024));
}

//==============================================================================
/// Display headers of PARTs
/// Visualiza cabeceras de PARTs.
//==============================================================================
void JSph::PrintHeadPart(){
  Log->Print("PART       PartTime      TotalSteps    Steps    Time/Sec   Finish time        ");
  Log->Print("=========  ============  ============  =======  =========  ===================");
  fflush(stdout);
}

//==============================================================================
/// Sets configuration for recordering of particles.
/// Establece configuracion para grabacion de particulas.
//==============================================================================
void JSph::ConfigSaveData(unsigned piece,unsigned pieces,std::string div){
  const char met[]="ConfigSaveData";
  //-Stores basic information of simulation data.
  JPartDataHead parthead;
  parthead.ConfigBasic(RunCode,AppName,CaseName,CasePosMin,CasePosMax
    ,Simulate2D,Simulate2DPosY,pieces,PartBeginFirst);
  MkInfo->ConfigPartDataHead(&parthead);
  parthead.ConfigCtes(Dp,H,CteB,RhopZero,Gamma,MassBound,MassFluid,Gravity);
  parthead.ConfigSimNp(NpDynamic,ReuseIds);
  parthead.ConfigSimMap(MapRealPosMin,MapRealPosMax);
  parthead.ConfigSimPeri(TpPeriFromPeriActive(PeriActive),PeriXinc,PeriYinc,PeriZinc);
  switch(TVisco){
    case VISCO_None:        parthead.ConfigVisco(JPartDataHead::VISCO_None      ,Visco,ViscoBoundFactor);  break;
    case VISCO_Artificial:  parthead.ConfigVisco(JPartDataHead::VISCO_Artificial,Visco,ViscoBoundFactor);  break;
    case VISCO_LaminarSPS:  parthead.ConfigVisco(JPartDataHead::VISCO_LaminarSPS,Visco,ViscoBoundFactor);  break;
    default: RunException(met,"Viscosity type is unknown.");
  }
  if(SvData&SDAT_Binx){
    Log->AddFileInfo(DirDataOut+"Part_Head.ibi4","Binary file with basic information of simulation data.");
    parthead.SaveFile(DirDataOut);
  }
  //-Configures object to store particles and information.  
  //-Configura objeto para grabacion de particulas e informacion.
  if(SvData&SDAT_Info || SvData&SDAT_Binx){
    DataBi4=new JPartDataBi4();
    DataBi4->Config(piece,pieces,DirDataOut,&parthead);
    if(div.empty())DataBi4->ConfigSimDiv(JPartDataBi4::DIV_None);
    else if(div=="X")DataBi4->ConfigSimDiv(JPartDataBi4::DIV_X);
    else if(div=="Y")DataBi4->ConfigSimDiv(JPartDataBi4::DIV_Y);
    else if(div=="Z")DataBi4->ConfigSimDiv(JPartDataBi4::DIV_Z);
    else RunException(met,"The division configuration is invalid.");
    if(SvData&SDAT_Binx)Log->AddFileInfo(DirDataOut+"Part_????.bi4","Binary file with particle data in different instants.");
    if(SvData&SDAT_Info)Log->AddFileInfo(DirDataOut+"PartInfo.ibi4","Binary file with execution information for each instant (input for PartInfo program).");
  }
  //-Configures object to store excluded particles.  
  //-Configura objeto para grabacion de particulas excluidas.
  if(SvData&SDAT_Binx){
    DataOutBi4=new JPartOutBi4Save();
    DataOutBi4->ConfigBasic(piece,pieces,RunCode,AppName,Simulate2D,DirDataOut);
    DataOutBi4->ConfigParticles(CaseNp,CaseNfixed,CaseNmoving,CaseNfloat,CaseNfluid);
    DataOutBi4->ConfigLimits(MapRealPosMin,MapRealPosMax,(RhopOut? RhopOutMin: 0),(RhopOut? RhopOutMax: 0));
    DataOutBi4->SaveInitial();
    Log->AddFileInfo(DirDataOut+"PartOut_???.obi4","Binary file with particles excluded during simulation (input for PartVtkOut program).");
  }
  //-Configures object to store data of floatings.
  //-Configura objeto para grabacion de datos de floatings.
  if(SvData&SDAT_Binx && FtCount){
    DataFloatBi4=new JPartFloatBi4Save();
    DataFloatBi4->Config(AppName,DirDataOut,MkInfo->GetMkBoundFirst(),FtCount);
    for(unsigned cf=0;cf<FtCount;cf++){
      DataFloatBi4->AddHeadData(cf,FtObjs[cf].mkbound,FtObjs[cf].begin,FtObjs[cf].count,FtObjs[cf].mass,FtObjs[cf].radius);
    }
    DataFloatBi4->SaveInitial();
    Log->AddFileInfo(DirDataOut+"PartFloat.fbi4","Binary file with floating body information for each instant (input for FloatingInfo program).");
  }
  //-Creates object to store excluded particles until recordering. 
  //-Crea objeto para almacenar las particulas excluidas hasta su grabacion.
  PartsOut=new JPartsOut();
}

//==============================================================================
/// Stores new excluded particles until recordering next PART.
/// Almacena nuevas particulas excluidas hasta la grabacion del proximo PART.
//==============================================================================
void JSph::AddParticlesOut(unsigned nout,const unsigned *idp,const tdouble3 *pos
  ,const tfloat3 *vel,const float *rhop,const typecode *code)
{
  PartsOut->AddParticles(nout,idp,pos,vel,rhop,code);
}

//==============================================================================
/// Manages excluded particles fixed, moving and floating before aborting the execution.
/// Gestiona particulas excluidas fixed, moving y floating antes de abortar la ejecucion.
//==============================================================================
void JSph::AbortBoundOut(unsigned nout,const unsigned *idp,const tdouble3 *pos,const tfloat3 *vel,const float *rhop,const typecode *code){
  //-Prepares data of excluded boundary particles.
  byte* type=new byte[nout];
  byte* motive=new byte[nout];
  unsigned outfixed=0,outmoving=0,outfloat=0;
  unsigned outpos=0,outrhop=0,outmove=0;
  for(unsigned p=0;p<nout;p++){
    //-Checks type of particle.
    switch(CODE_GetType(code[p])){
    case CODE_TYPE_FIXED:     type[p]=0;  outfixed++;   break;
    case CODE_TYPE_MOVING:    type[p]=1;  outmoving++;  break; 
    case CODE_TYPE_FLOATING:  type[p]=2;  outfloat++;   break; 
    default:                  type[p]=99;               break; 
    }
    //-Checks reason for exclusion.
    switch(CODE_GetSpecialValue(code[p])){
    case CODE_OUTPOS:   motive[p]=1; outpos++;   break;
    case CODE_OUTRHOP:  motive[p]=2; outrhop++;  break; 
    case CODE_OUTMOVE:  motive[p]=3; outmove++;  break; 
    default:            motive[p]=0;             break; 
    }
  }
  //-Shows excluded particles information.
  Log->Print(" ");
  Log->Print("*** ERROR: Some boundary particle was excluded. ***");
  Log->Printf("TimeStep: %f  (Nstep: %u)",TimeStep,Nstep);
  unsigned npunknown=nout-outfixed-outmoving-outfloat;
  if(!npunknown)Log->Printf("Total boundary: %u  (fixed=%u  moving=%u  floating=%u)",nout,outfixed,outmoving,outfloat);
  else Log->Printf("Total boundary: %u  (fixed=%u  moving=%u  floating=%u  UNKNOWN=%u)",nout,outfixed,outmoving,outfloat,npunknown);
  npunknown=nout-outpos-outrhop-outmove;
  if(!npunknown)Log->Printf("Excluded for: position=%u  rhop=%u  velocity=%u",outpos,outrhop,outmove);
  else Log->Printf("Excluded for: position=%u  rhop=%u  velocity=%u  UNKNOWN=%u",outpos,outrhop,outmove,npunknown);
  Log->Print(" ");
  //-Creates VTK file.
  std::vector<JFormatFiles2::StScalarData> fields;
  fields.push_back(JFormatFiles2::DefineField("Idp"   ,JFormatFiles2::UInt32 ,1,idp));
  fields.push_back(JFormatFiles2::DefineField("Vel"   ,JFormatFiles2::Float32,3,vel));
  fields.push_back(JFormatFiles2::DefineField("Rhop"  ,JFormatFiles2::Float32,1,rhop));
  fields.push_back(JFormatFiles2::DefineField("Type"  ,JFormatFiles2::UChar8 ,1,type));
  fields.push_back(JFormatFiles2::DefineField("Motive",JFormatFiles2::UChar8 ,1,motive));
  const string file=DirOut+"Error_BoundaryOut.vtk";
  Log->AddFileInfo(file,"Saves the excluded boundary particles.");
  JFormatFiles2::SaveVtk(file,nout,pos,fields);
  //-Aborts execution.
  RunException("AbortBoundOut","Fixed, moving or floating particles were excluded. Check VTK file Error_BoundaryOut.vtk with excluded particles.");
}

//==============================================================================
/// Returns dynamic memory pointer with data transformed in tfloat3.
/// THE POINTER MUST BE RELEASED AFTER USING IT.
///
/// Devuelve puntero de memoria dinamica con los datos transformados en tfloat3.
/// EL PUNTERO DEBE SER LIBERADO DESPUES DE USARLO.
//==============================================================================
tfloat3* JSph::GetPointerDataFloat3(unsigned n,const tdouble3* v)const{
  tfloat3* v2=new tfloat3[n];
  for(unsigned c=0;c<n;c++)v2[c]=ToTFloat3(v[c]);
  return(v2);
}

//==============================================================================
/// Stores files of particle data.
/// Graba los ficheros de datos de particulas.
//==============================================================================
void JSph::SavePartData(unsigned npok,unsigned nout,const unsigned *idp,const tdouble3 *pos,const tfloat3 *vel,const float *rhop,unsigned ndom,const tdouble3 *vdom,const StInfoPartPlus *infoplus){
  //-Stores particle data and/or information in bi4 format.
  //-Graba datos de particulas y/o informacion en formato bi4.
  if(DataBi4){
    tfloat3* posf3=NULL;
    TimerPart.Stop();
    JBinaryData* bdpart=DataBi4->AddPartInfo(Part,TimeStep,npok,nout,Nstep,TimerPart.GetElapsedTimeD()/1000.,vdom[0],vdom[1],TotalNp);
    if(TStep==STEP_Symplectic)bdpart->SetvDouble("SymplecticDtPre",SymplecticDtPre);
    if(UseDEM)bdpart->SetvDouble("DemDtForce",DemDtForce); //(DEM)
    if(infoplus && SvData&SDAT_Info){
      bdpart->SetvDouble("dtmean",(!Nstep? 0: (TimeStep-TimeStepM1)/(Nstep-PartNstep)));
      bdpart->SetvDouble("dtmin",(!Nstep? 0: PartDtMin));
      bdpart->SetvDouble("dtmax",(!Nstep? 0: PartDtMax));
      if(DtFixed)bdpart->SetvDouble("dterror",DtFixed->GetDtError(true));
      bdpart->SetvDouble("timesim",infoplus->timesim);
      bdpart->SetvUint("nct",infoplus->nct);
      bdpart->SetvUint("npbin",infoplus->npbin);
      bdpart->SetvUint("npbout",infoplus->npbout);
      bdpart->SetvUint("npf",infoplus->npf);
      bdpart->SetvUint("npbper",infoplus->npbper);
      bdpart->SetvUint("npfper",infoplus->npfper);
      bdpart->SetvUint("newnp",infoplus->newnp);
      bdpart->SetvLlong("cpualloc",infoplus->memorycpualloc);
      if(infoplus->gpudata){
        bdpart->SetvLlong("nctalloc",infoplus->memorynctalloc);
        bdpart->SetvLlong("nctused",infoplus->memorynctused);
        bdpart->SetvLlong("npalloc",infoplus->memorynpalloc);
        bdpart->SetvLlong("npused",infoplus->memorynpused);
      }
    }
    if(SvData&SDAT_Binx){
      if(SvDouble)DataBi4->AddPartData(npok,idp,pos,vel,rhop);
      else{
        posf3=GetPointerDataFloat3(npok,pos);
        DataBi4->AddPartData(npok,idp,posf3,vel,rhop);
      }
      float *press=NULL;
      if(0){//-Example saving a new array (Pressure) in files BI4.
        press=new float[npok];
        for(unsigned p=0;p<npok;p++)press[p]=(idp[p]>=CaseNbound? CteB*(pow(rhop[p]/RhopZero,Gamma)-1.0f): 0.f);
        DataBi4->AddPartData("Pressure",npok,press);
      }
      DataBi4->SaveFilePart();
      delete[] press; press=NULL;//-Memory must to be deallocated after saving file because DataBi4 uses this memory space.
    }
    if(SvData&SDAT_Info)DataBi4->SaveFileInfo();
    delete[] posf3;
  }

  //-Stores VTK nd/or CSV files.
  if((SvData&SDAT_Csv) || (SvData&SDAT_Vtk)){
    //-Generates array with posf3 and type of particle.
    tfloat3* posf3=GetPointerDataFloat3(npok,pos);
    byte *type=new byte[npok];
    for(unsigned p=0;p<npok;p++){
      const unsigned id=idp[p];
      type[p]=(id>=CaseNbound? 3: (id<CaseNfixed? 0: (id<CaseNpb? 1: 2)));
    }
    //-Defines fields to be stored.
    JFormatFiles2::StScalarData fields[8];
    unsigned nfields=0;
    if(idp){   fields[nfields]=JFormatFiles2::DefineField("Idp" ,JFormatFiles2::UInt32 ,1,idp);   nfields++; }
    if(vel){   fields[nfields]=JFormatFiles2::DefineField("Vel" ,JFormatFiles2::Float32,3,vel);   nfields++; }
    if(rhop){  fields[nfields]=JFormatFiles2::DefineField("Rhop",JFormatFiles2::Float32,1,rhop);  nfields++; }
    if(type){  fields[nfields]=JFormatFiles2::DefineField("Type",JFormatFiles2::UChar8 ,1,type);  nfields++; }
    if(SvData&SDAT_Vtk)JFormatFiles2::SaveVtk(DirDataOut+fun::FileNameSec("PartVtk.vtk",Part),npok,posf3,nfields,fields);
    if(SvData&SDAT_Csv)JFormatFiles2::SaveCsv(DirDataOut+fun::FileNameSec("PartCsv.csv",Part),CsvSepComa,npok,posf3,nfields,fields);
    //-Deallocate of memory.
    delete[] posf3;
    delete[] type; 
  }

  //-Stores data of excluded particles.
  if(DataOutBi4 && PartsOut->GetCount()){
    DataOutBi4->SavePartOut(SvDouble,Part,TimeStep,PartsOut->GetCount(),PartsOut->GetIdpOut(),NULL,PartsOut->GetPosOut(),PartsOut->GetVelOut(),PartsOut->GetRhopOut(),PartsOut->GetMotiveOut());
  }

  //-Stores data of floating bodies.
  if(DataFloatBi4){
    for(unsigned cf=0;cf<FtCount;cf++)DataFloatBi4->AddPartData(cf,FtObjs[cf].center,FtObjs[cf].fvel,FtObjs[cf].fomega);
    DataFloatBi4->SavePartFloat(Part,TimeStep,(UseDEM? DemDtForce: 0));
  }

  //-Empties stock of excluded particles.
  //-Vacia almacen de particulas excluidas.
  PartsOut->Clear();
}

//==============================================================================
/// Generates data output files.
/// Genera los ficheros de salida de datos.
//==============================================================================
void JSph::SaveData(unsigned npok,const unsigned *idp,const tdouble3 *pos,const tfloat3 *vel,const float *rhop
  ,unsigned ndom,const tdouble3 *vdom,const StInfoPartPlus *infoplus)
{
  const char met[]="SaveData";
  string suffixpartx=fun::PrintStr("_%04d",Part);

  //-Counts new excluded particles.
  //-Contabiliza nuevas particulas excluidas.
  const unsigned noutpos=PartsOut->GetOutPosCount(),noutrhop=PartsOut->GetOutRhopCount(),noutmove=PartsOut->GetOutMoveCount();
  const unsigned nout=noutpos+noutrhop+noutmove;
  if(nout!=PartsOut->GetCount())RunException(met,"Excluded particles with unknown reason.");
  AddOutCount(noutpos,noutrhop,noutmove);

  //-Stores data files of particles.
  SavePartData(npok,nout,idp,pos,vel,rhop,ndom,vdom,infoplus);

  //-Reinitialises limits of dt. | Reinicia limites de dt.
  PartDtMin=DBL_MAX; PartDtMax=-DBL_MAX;

  //-Computation of time.
  if(Part>PartIni||Nstep){
    TimerPart.Stop();
    double tpart=TimerPart.GetElapsedTimeD()/1000;
    double tseg=tpart/(TimeStep-TimeStepM1);
    TimerSim.Stop();
    double tcalc=TimerSim.GetElapsedTimeD()/1000;
    double tleft=(tcalc/(TimeStep-TimeStepIni))*(TimeMax-TimeStep);
    Log->Printf("Part%s  %12.6f  %12d  %7d  %9.2f  %14s",suffixpartx.c_str(),TimeStep,(Nstep+1),Nstep-PartNstep,tseg,fun::GetDateTimeAfter(int(tleft)).c_str());
  }
  else Log->Printf("Part%s        %u particles successfully stored",suffixpartx.c_str(),npok);   
  
  //-Shows info of the excluded particles.
  if(nout){
    PartOut+=nout;
    Log->Printf("  Particles out: %u  (total: %u)",nout,PartOut);
  }

  //-Cheks number of excluded particles.
  if(nout){
    //-Cheks number of excluded particles in one PART.
    if(PartsOutWrn<=100 && nout>=float(infoplus->npf)*(float(PartsOutWrn)/100.f)){
      Log->PrintfWarning("More than %d%% of current fluid particles were excluded in one PART (t:%g, nstep:%u)",PartsOutWrn,TimeStep,Nstep);
      if(PartsOutWrn==1)PartsOutWrn=2;
      else if(PartsOutWrn==2)PartsOutWrn=5;
      else if(PartsOutWrn==5)PartsOutWrn=10;
      else PartsOutWrn+=10;
    }
    //-Cheks number of total excluded particles.
    const unsigned noutt=GetOutPosCount()+GetOutRhopCount()+GetOutMoveCount();
    if(PartsOutTotWrn<=100 && noutt>=float(TotalNp)*(float(PartsOutTotWrn)/100.f)){
      Log->PrintfWarning("More than %d%% of particles were excluded (t:%g, nstep:%u)",PartsOutTotWrn,TimeStep,Nstep);
      PartsOutTotWrn+=10;
    }
  }  

  if(SvDomainVtk)SaveDomainVtk(ndom,vdom);
  if(SaveDt)SaveDt->SaveData();
  if(GaugeSystem)GaugeSystem->SaveResults(Part);
}

//==============================================================================
/// Generates VTK file with domain of the particles.
/// Genera fichero VTK con el dominio de las particulas.
//==============================================================================
void JSph::SaveDomainVtk(unsigned ndom,const tdouble3 *vdom)const{ 
  if(vdom){
    string fname=fun::FileNameSec("Domain.vtk",Part);
    tfloat3 *vdomf3=new tfloat3[ndom*2];
    for(unsigned c=0;c<ndom*2;c++)vdomf3[c]=ToTFloat3(vdom[c]);
    JFormatFiles2::SaveVtkBoxes(DirDataOut+fname,ndom,vdomf3,H*0.5f);
    delete[] vdomf3;
  }
}

//==============================================================================
/// Saves initial domain of simulation in a VTK file (CasePosMin/Max, 
/// MapRealPosMin/Max and Map_PosMin/Max).
///
/// Graba dominio inicial de simulacion en fichero VTK (CasePosMin/Max, 
/// MapRealPosMin/Max and Map_PosMin/Max).
//==============================================================================
void JSph::SaveInitialDomainVtk()const{
  const unsigned nbox=(MapRealPosMin!=Map_PosMin || MapRealPosMax!=Map_PosMax? 3: 2);
  tfloat3 *vdomf3=new tfloat3[nbox*2];
  vdomf3[0]=ToTFloat3(CasePosMin);
  vdomf3[1]=ToTFloat3(CasePosMax);
  vdomf3[2]=ToTFloat3(MapRealPosMin);
  vdomf3[3]=ToTFloat3(MapRealPosMax);
  if(nbox==3){
    vdomf3[4]=ToTFloat3(Map_PosMin);
    vdomf3[5]=ToTFloat3(Map_PosMax);
  }
  const string file=DirOut+"CfgInit_Domain.vtk";
  Log->AddFileInfo(file,"Saves the limits of the case and the simulation domain limits.");
  JFormatFiles2::SaveVtkBoxes(file,nbox,vdomf3,0);
  delete[] vdomf3;
}

//==============================================================================
/// Returns size of VTK file with map cells.
/// Devuelve tamaño de fichero VTK con las celdas del mapa.
//==============================================================================
unsigned JSph::SaveMapCellsVtkSize()const{
  const tuint3 cells=Map_Cells;
  unsigned nlin=cells.x+cells.z+2;//-Back lines.
  if(!Simulate2D){
    nlin+=cells.x+cells.y+2;//-Bottom lines.
    nlin+=cells.y+cells.z+2;//-Left lines.
  }
  const unsigned slin=sizeof(tfloat3)*2+sizeof(int)*4; //-Size per line is 40 bytes.
  return(nlin*slin);
}

//==============================================================================
/// Generates VTK file with map cells.
/// Genera fichero VTK con las celdas del mapa.
//==============================================================================
void JSph::SaveMapCellsVtk(float scell)const{
  const tuint3 cells=Map_Cells;
  tdouble3 pmin=MapRealPosMin;
  tdouble3 pmax=pmin+TDouble3(scell*cells.x,scell*cells.y,scell*cells.z);
  if(Simulate2D)pmin.y=pmax.y=Simulate2DPosY;
  //-Creates lines.
  std::vector<JFormatFiles2::StShapeData> shapes;
  //-Back lines.
  tdouble3 p0=TDouble3(pmin.x,pmax.y,pmin.z),p1=TDouble3(pmin.x,pmax.y,pmax.z);
  for(unsigned cx=0;cx<=cells.x;cx++)shapes.push_back(JFormatFiles2::DefineShape_Line(p0+TDouble3(scell*cx,0,0),p1+TDouble3(scell*cx,0,0),0,0));
  p1=TDouble3(pmax.x,pmax.y,pmin.z);
  for(unsigned cz=0;cz<=cells.z;cz++)shapes.push_back(JFormatFiles2::DefineShape_Line(p0+TDouble3(0,0,scell*cz),p1+TDouble3(0,0,scell*cz),0,0));
  if(!Simulate2D){
    //-Bottom lines.
    p0=TDouble3(pmin.x,pmin.y,pmin.z),p1=TDouble3(pmax.x,pmin.y,pmin.z);
    for(unsigned cy=0;cy<=cells.y;cy++)shapes.push_back(JFormatFiles2::DefineShape_Line(p0+TDouble3(0,scell*cy,0),p1+TDouble3(0,scell*cy,0),1,0));
    p1=TDouble3(pmin.x,pmax.y,pmin.z);
    for(unsigned cx=0;cx<=cells.x;cx++)shapes.push_back(JFormatFiles2::DefineShape_Line(p0+TDouble3(scell*cx,0,0),p1+TDouble3(scell*cx,0,0),1,0));
    //-Left lines.
    p0=TDouble3(pmin.x,pmin.y,pmin.z),p1=TDouble3(pmin.x,pmax.y,pmin.z);
    for(unsigned cz=0;cz<=cells.z;cz++)shapes.push_back(JFormatFiles2::DefineShape_Line(p0+TDouble3(0,0,scell*cz),p1+TDouble3(0,0,scell*cz),2,0));
    p1=TDouble3(pmin.x,pmin.y,pmax.z);
    for(unsigned cy=0;cy<=cells.y;cy++)shapes.push_back(JFormatFiles2::DefineShape_Line(p0+TDouble3(0,scell*cy,0),p1+TDouble3(0,scell*cy,0),2,0));
  }
  const string file=DirOut+"CfgInit_MapCells.vtk";
  Log->AddFileInfo(file,"Saves the cell division of the simulation domain.");
  JFormatFiles2::SaveVtkShapes(file,"axis","",shapes);
}

//==============================================================================
/// Adds basic information of resume to hinfo & dinfo.
/// Añade la informacion basica de resumen a hinfo y dinfo.
//==============================================================================
void JSph::GetResInfo(float tsim,float ttot,const std::string &headplus,const std::string &detplus,std::string &hinfo,std::string &dinfo){
  hinfo=hinfo+"#RunName;RunCode;DateTime;Np;TSimul;TSeg;TTotal;MemCpu;MemGpu;Steps;PartFiles;PartsOut;MaxParticles;MaxCells;Hw;StepAlgo;Kernel;Viscosity;ViscoValue;DeltaSPH;TMax;Nbound;Nfixed;H;RhopOut;PartsRhopOut;PartsVelOut;CellMode"+headplus;
  dinfo=dinfo+ RunName+ ";"+ RunCode+ ";"+ RunTimeDate+ ";"+ fun::UintStr(CaseNp);
  dinfo=dinfo+ ";"+ fun::FloatStr(tsim)+ ";"+ fun::FloatStr(tsim/float(TimeStep))+ ";"+ fun::FloatStr(ttot);
  dinfo=dinfo+ ";"+ fun::LongStr(MaxMemoryCpu)+ ";"+ fun::LongStr(MaxMemoryGpu);
  const unsigned nout=GetOutPosCount()+GetOutRhopCount()+GetOutMoveCount();
  dinfo=dinfo+ ";"+ fun::IntStr(Nstep)+ ";"+ fun::IntStr(Part)+ ";"+ fun::UintStr(nout);
  dinfo=dinfo+ ";"+ fun::UintStr(MaxParticles)+ ";"+ fun::UintStr(MaxCells);
  dinfo=dinfo+ ";"+ Hardware+ ";"+ GetStepName(TStep)+ ";"+ GetKernelName(TKernel)+ ";"+ GetViscoName(TVisco)+ ";"+ fun::FloatStr(Visco);
  dinfo=dinfo+ ";"+ fun::FloatStr(DeltaSph,"%G")+ ";"+ fun::FloatStr(float(TimeMax));
  dinfo=dinfo+ ";"+ fun::UintStr(CaseNbound)+ ";"+ fun::UintStr(CaseNfixed)+ ";"+ fun::FloatStr(H);
  std::string rhopcad;
  if(RhopOut)rhopcad=fun::PrintStr("(%G-%G)",RhopOutMin,RhopOutMax); else rhopcad="None";
  dinfo=dinfo+ ";"+ rhopcad+ ";"+ fun::UintStr(GetOutRhopCount())+ ";"+ fun::UintStr(GetOutMoveCount())+ ";"+ GetNameCellMode(CellMode)+ detplus;
}

//==============================================================================
/// Generates file Run.csv with resume of execution.
/// Genera fichero Run.csv con resumen de ejecucion.
//==============================================================================
void JSph::SaveRes(float tsim,float ttot,const std::string &headplus,const std::string &detplus){
  const char* met="SaveRes";
  const string fname=DirOut+"Run.csv";
  Log->AddFileInfo(fname,"One line CSV file with execution parameters and other simulation data.");
  ofstream pf;
  pf.open(fname.c_str());
  if(pf){
    string hinfo,dinfo;
    GetResInfo(tsim,ttot,headplus,detplus,hinfo,dinfo);
    pf << fun::StrCsvSep(CsvSepComa,hinfo) << endl << fun::StrCsvSep(CsvSepComa,dinfo) << endl;
    if(pf.fail())RunException(met,"Failed writing to file.",fname);
    pf.close();
  }
  else RunException(met,"File could not be opened.",fname);
}

//==============================================================================
/// Shows resume of execution.
/// Muestra resumen de ejecucion.
//==============================================================================
void JSph::ShowResume(bool stop,float tsim,float ttot,bool all,std::string infoplus){
  Log->Printf("\n[Simulation %s  %s]",(stop? "INTERRUPTED": "finished"),fun::GetDateTime().c_str());
  Log->Printf("Particles of simulation (initial): %u",CaseNp);
  if(NpDynamic)Log->Printf("Particles of simulation (total)..: %llu",TotalNp);
  if(all){
    Log->Printf("DTs adjusted to DtMin............: %d",DtModif);
    const unsigned nout=GetOutPosCount()+GetOutRhopCount()+GetOutMoveCount();
    Log->Printf("Excluded particles...............: %d",nout);
    if(GetOutRhopCount())Log->Printf("Excluded particles due to RhopOut: %u",GetOutRhopCount());
    if(GetOutMoveCount())Log->Printf("Excluded particles due to Velocity: %u",GetOutMoveCount());
  }
  Log->Printf("Total Runtime....................: %f sec.",ttot);
  Log->Printf("Simulation Runtime...............: %f sec.",tsim);
  if(all){
    float tseg=tsim/float(TimeStep);
    float nstepseg=float(Nstep)/tsim;
    Log->Printf("Time per second of simulation....: %f sec.",tseg);
    Log->Printf("Steps per second.................: %f",nstepseg);
    Log->Printf("Steps of simulation..............: %d",Nstep);
    Log->Printf("PART files.......................: %d",Part-PartIni);
    while(!infoplus.empty()){
      string lin=fun::StrSplit("#",infoplus);
      if(!lin.empty()){
        string tex=fun::StrSplit("=",lin);
        string val=fun::StrSplit("=",lin);
        while(tex.size()<33)tex=tex+".";
        Log->Print(tex+": "+val);
      }
    }
  }
  Log->Printf("Maximum number of particles......: %u",MaxParticles);
  Log->Printf("Maximum number of cells..........: %u",MaxCells);
  Log->Printf("CPU Memory.......................: %lld (%.2f MB)",MaxMemoryCpu,double(MaxMemoryCpu)/(1024*1024));
  if(MaxMemoryGpu)Log->Printf("GPU Memory.......................: %lld (%.2f MB)",MaxMemoryGpu,double(MaxMemoryGpu)/(1024*1024));
}

//==============================================================================
/// Returns text about PosDouble configuration.
/// Devuelve texto sobre la configuracion de PosDouble.
//==============================================================================
std::string JSph::GetPosDoubleName(bool psingle,bool svdouble){
  string tx;
  if(psingle && !svdouble)tx="0: Uses and stores in single precision";
  else if(!psingle && !svdouble)tx="1: Uses double and stores in single precision";
  else if(!psingle && svdouble)tx="2: Uses and stores in double precision";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns the name of the time algorithm in text format.
/// Devuelve el nombre del algoritmo en texto.
//==============================================================================
std::string JSph::GetStepName(TpStep tstep){
  string tx;
  if(tstep==STEP_Verlet)tx="Verlet";
  else if(tstep==STEP_Symplectic)tx="Symplectic";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns the name of the kernel in text format.
/// Devuelve el nombre del kernel en texto.
//==============================================================================
std::string JSph::GetKernelName(TpKernel tkernel){
  string tx;
  if(tkernel==KERNEL_Cubic)tx="Cubic";
  else if(tkernel==KERNEL_Wendland)tx="Wendland";
  else if(tkernel==KERNEL_Gaussian)tx="Gaussian";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns value of viscosity in text format.
/// Devuelve el nombre de la viscosidad en texto.
//==============================================================================
std::string JSph::GetViscoName(TpVisco tvisco){
  string tx;
  if(tvisco==VISCO_Artificial)tx="Artificial";
  else if(tvisco==VISCO_LaminarSPS)tx="Laminar+SPS";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns value of DeltaSPH in text format.
/// Devuelve el valor de DeltaSPH en texto.
//==============================================================================
std::string JSph::GetDeltaSphName(TpDeltaSph tdelta){
  string tx;
  if(tdelta==DELTA_None)tx="None";
  else if(tdelta==DELTA_Dynamic)tx="Dynamic";
  else if(tdelta==DELTA_DynamicExt)tx="DynamicExt";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns value of Shifting in text format.
/// Devuelve el valor de Shifting en texto.
//==============================================================================
std::string JSph::GetShiftingName(TpShifting tshift){
  string tx;
  if(tshift==SHIFT_None)tx="None";
  else if(tshift==SHIFT_NoBound)tx="NoBound";
  else if(tshift==SHIFT_NoFixed)tx="NoFixed";
  else if(tshift==SHIFT_Full)tx="Full";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns string with the name of timer and value.
/// Devuelve string con el nombre del temporizador y su valor.
//==============================================================================
std::string JSph::TimerToText(const std::string &name,float value){
  string ret=name;
  while(ret.length()<33)ret+=".";
  return(ret+": "+fun::FloatStr(value/1000)+" sec.");
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSph::DgSaveVtkParticlesCpu(std::string filename,int numfile,unsigned pini,unsigned pfin
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,const tfloat4 *velrhop
  ,const tfloat3 *ace)const
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Allocates memory.
  const unsigned np=pfin-pini;
  tfloat3 *xpos=new tfloat3[np];
  tfloat3 *xvel=new tfloat3[np];
  tfloat3 *xace=(ace? new tfloat3[np]: NULL);
  float *xrhop=new float[np];
  byte *xtype=new byte[np];
  byte *xkind=new byte[np];
  for(unsigned p=0;p<np;p++){
    xpos[p]=ToTFloat3(pos[p+pini]);
    tfloat4 vr=velrhop[p+pini];
    xvel[p]=TFloat3(vr.x,vr.y,vr.z);
    if(xace)xace[p]=ace[p+pini];
    xrhop[p]=vr.w;
    typecode t=CODE_GetType(code[p+pini]);
    xtype[p]=(t==CODE_TYPE_FIXED? 0: (t==CODE_TYPE_MOVING? 1: (t==CODE_TYPE_FLOATING? 2: 3)));
    typecode k=CODE_GetSpecialValue(code[p+pini]);
    xkind[p]=(k==CODE_NORMAL? 0: (k==CODE_PERIODIC? 1: (k==CODE_OUTIGNORE? 2: 3)));
  }
  //-Generates VTK file.
  JFormatFiles2::StScalarData fields[10];
  unsigned nfields=0;
  if(idp){   fields[nfields]=JFormatFiles2::DefineField("Idp" ,JFormatFiles2::UInt32 ,1,idp+pini); nfields++; }
  if(xtype){ fields[nfields]=JFormatFiles2::DefineField("Type",JFormatFiles2::UChar8 ,1,xtype);    nfields++; }
  if(xkind){ fields[nfields]=JFormatFiles2::DefineField("Kind",JFormatFiles2::UChar8 ,1,xkind);    nfields++; }
  if(xvel){  fields[nfields]=JFormatFiles2::DefineField("Vel" ,JFormatFiles2::Float32,3,xvel);     nfields++; }
  if(xrhop){ fields[nfields]=JFormatFiles2::DefineField("Rhop",JFormatFiles2::Float32,1,xrhop);    nfields++; }
  if(xace){  fields[nfields]=JFormatFiles2::DefineField("Ace" ,JFormatFiles2::Float32,3,xace);     nfields++; }
  //string fname=DirOut+fun::FileNameSec("DgParts.vtk",numfile);
  JFormatFiles2::SaveVtk(filename,np,xpos,nfields,fields);
  //-Deallocates memory.
  delete[] xpos;
  delete[] xtype;
  delete[] xkind;
  delete[] xvel;
  delete[] xrhop;
  delete[] xace;
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSph::DgSaveVtkParticlesCpu(std::string filename,int numfile,unsigned pini,unsigned pfin,const tfloat3 *pos,const byte *check,const unsigned *idp,const tfloat3 *vel,const float *rhop){
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Reserva memoria basica.
  const unsigned n=pfin-pini;
  unsigned *num=new unsigned[n];
  for(unsigned p=0;p<n;p++)num[p]=p;
  //-Generates VTK file.
  JFormatFiles2::StScalarData fields[10];
  unsigned nfields=0;
  if(idp){   fields[nfields]=JFormatFiles2::DefineField("Idp"  ,JFormatFiles2::UInt32 ,1,idp+pini);   nfields++; }
  if(vel){   fields[nfields]=JFormatFiles2::DefineField("Vel"  ,JFormatFiles2::Float32,3,vel+pini);   nfields++; }
  if(rhop){  fields[nfields]=JFormatFiles2::DefineField("Rhop" ,JFormatFiles2::Float32,1,rhop+pini);  nfields++; }
  if(check){ fields[nfields]=JFormatFiles2::DefineField("Check",JFormatFiles2::UChar8 ,1,check+pini); nfields++; }
  if(num){   fields[nfields]=JFormatFiles2::DefineField("Num"  ,JFormatFiles2::UInt32 ,1,num);        nfields++; }
  //-Generates VTK file.
  JFormatFiles2::SaveVtk(filename,n,pos+pini,nfields,fields);
  //-Deallocates memory.
  delete[] num;
}

//==============================================================================
/// Saves CSV file with particle data (degug).
/// Graba fichero CSV con datos de las particulas (degug).
//==============================================================================
void JSph::DgSaveCsvParticlesCpu(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string head,const tfloat3 *pos,const unsigned *idp,const tfloat3 *vel,const float *rhop,const float *ar,const tfloat3 *ace,const tfloat3 *vcorr){
  const char met[]="DgSaveCsvParticlesCpu";
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirOut+filename;
  //-Generates CSV file.
  ofstream pf;
  pf.open(filename.c_str());
  if(pf){
    if(!head.empty())pf << head << endl;
    pf << "Num";
    if(idp)pf << ";Idp";
    if(pos)pf << ";PosX;PosY;PosZ";
    if(vel)pf << ";VelX;VelY;VelZ";
    if(rhop)pf << ";Rhop";
    if(ar)pf << ";Ar";
    if(ace)pf << ";AceX;AceY;AceZ";
    if(vcorr)pf << ";VcorrX;VcorrY;VcorrZ";
    pf << endl;
    const char fmt1[]="%f"; //="%24.16f";
    const char fmt3[]="%f;%f;%f"; //="%24.16f;%24.16f;%24.16f";
    for(unsigned p=pini;p<pfin;p++){
      pf << fun::UintStr(p-pini);
      if(idp)pf << ";" << fun::UintStr(idp[p]);
      if(pos)pf << ";" << fun::Float3Str(pos[p],fmt3);
      if(vel)pf << ";" << fun::Float3Str(vel[p],fmt3);
      if(rhop)pf << ";" << fun::FloatStr(rhop[p],fmt1);
      if(ar)pf << ";" << fun::FloatStr(ar[p],fmt1);
      if(ace)pf << ";" << fun::Float3Str(ace[p],fmt3);
      if(vcorr)pf << ";" << fun::Float3Str(vcorr[p],fmt3);
      pf << endl;
    }
    if(pf.fail())RunException(met,"Failed writing to file.",filename);
    pf.close();
  }
  else RunException(met,"File could not be opened.",filename);
}


