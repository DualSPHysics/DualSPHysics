//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "FunGeo3d.h"
#include "FunSphKernelsCfg.h"
#include "FunSphKernel.h"
#include "JPartDataHead.h"
#include "JSphMk.h"
#include "JDsPartsInit.h"
#include "JPartsLoad4.h"
#include "JDsMotion.h"
#include "JXml.h"
#include "JCaseCtes.h"
#include "JCaseEParms.h"
#include "JCaseParts.h"
#include "JDsFixedDt.h"
#include "JDsSaveDt.h"
#include "JDsOutputTime.h"
#include "JDsViscoInput.h"
#include "JDsGaugeSystem.h"
#include "JWaveGen.h"
#include "JMLPistons.h"
#include "JRelaxZones.h"
#include "JChronoObjects.h"
#include "JDsMooredFloatings.h"
#include "JDsFtForcePoints.h"
#include "JDsAccInput.h"
#include "JPartDataBi4.h"
#include "JPartOutBi4Save.h"
#include "JPartFloatBi4.h"
#include "JDsPartsOut.h"
#include "JSphShifting.h"
#include "JDsDamping.h"
#include "JDsInitialize.h"
#include "JSphInOut.h"
#include "JSphBoundCorr.h"
#include "JFtMotionSave.h"  //<vs_ftmottionsv>
#include "JDsPips.h"
#include "JLinearValue.h"
#include "JPartNormalData.h"
#include "JNormalsMarrone.h"
#include "JDataArrays.h"
#include "JOutputCsv.h"
#include "JVtkLib.h"
#include "JNumexLib.h"
#include "JCaseUserVars.h"
#include <algorithm>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JSph::JSph(bool cpu,bool mgpu,bool withmpi):Cpu(cpu),Mgpu(mgpu),WithMpi(withmpi)
{
  ClassName="JSph";
  DgNum=0;
  DataBi4=NULL;
  DataOutBi4=NULL;
  DataFloatBi4=NULL;
  PartsOut=NULL;
  Log=NULL;
  CfgRun=NULL;
  ViscoTime=NULL;
  FixedDt=NULL;
  SaveDt=NULL;
  OutputTime=NULL;
  MkInfo=NULL;
  PartsInit=NULL;
  DsMotion=NULL;
  FtObjs=NULL;
  FtLinearVel=NULL;
  FtAngularVel=NULL;
  FtLinearForce=NULL;
  FtAngularForce=NULL;
  DemData=NULL;
  GaugeSystem=NULL;
  WaveGen=NULL;
  MLPistons=NULL;
  RelaxZones=NULL;
  ChronoObjects=NULL;
  Moorings=NULL;
  ForcePoints=NULL;
  Shifting=NULL;
  Damping=NULL;
  AccInput=NULL;
  PartsLoaded=NULL;
  InOut=NULL;
  BoundCorr=NULL;
  FtMotSave=NULL; //<vs_ftmottionsv>
  DsPips=NULL;
  NuxLib=NULL;
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
  delete FixedDt;       FixedDt=NULL;
  delete SaveDt;        SaveDt=NULL;
  delete OutputTime;    OutputTime=NULL;
  delete MkInfo;        MkInfo=NULL;
  delete PartsInit;     PartsInit=NULL;
  delete DsMotion;      DsMotion=NULL;
  AllocMemoryFloating(0,false);
  delete[] DemData;     DemData=NULL;
  delete GaugeSystem;   GaugeSystem=NULL;
  delete WaveGen;       WaveGen=NULL;
  delete MLPistons;     MLPistons=NULL;
  delete RelaxZones;    RelaxZones=NULL;
  delete ChronoObjects; ChronoObjects=NULL;
  delete Moorings;      Moorings=NULL;
  delete ForcePoints;   ForcePoints=NULL;
  delete Shifting;      Shifting=NULL;
  delete Damping;       Damping=NULL;
  delete AccInput;      AccInput=NULL; 
  delete PartsLoaded;   PartsLoaded=NULL;
  delete InOut;         InOut=NULL;
  delete BoundCorr;     BoundCorr=NULL;
  delete FtMotSave;     FtMotSave=NULL;   //<vs_ftmottionsv>
  delete DsPips;        DsPips=NULL;
  delete NuxLib;        NuxLib=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSph::InitVars(){
  ClearCfgDomain();
  OutPosCount=OutRhopCount=OutMoveCount=0;
  Simulate2D=false;
  Simulate2DPosY=0;
  Symmetry=false;
  Stable=false;
  SvPosDouble=false;
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
  KCubic ={0,0,0,0,0,0,0,0};
  KWend  ={0,0};
  TVisco=VISCO_None;
  TDensity=DDT_None; DDTValue=0; DDTArray=false;
  ShiftingMode=(Shifting? Shifting->GetShiftMode(): SHIFT_None);
  Visco=0; ViscoBoundFactor=1;
  TBoundary=BC_DBC;
  SlipMode=SLIP_Vel0;
  MdbcCorrector=false;
  MdbcFastSingle=true;
  MdbcThreshold=0;
  UseNormals=false;
  UseNormalsFt=false;
  SvNormals=false;
  UseDEM=false;  //(DEM)
  delete[] DemData; DemData=NULL;  //(DEM)
  UseChrono=false;
  RhopOut=true; RhopOutMin=700; RhopOutMax=1300;
  TimeMax=TimePart=0;
  NstepsBreak=0;
  SvAllSteps=false;
  TerminateMt=0;
  DtIni=DtMin=0; CoefDtMin=0; DtAllParticles=false;
  PartsOutMax=0;
  NpMinimum=0;
  PartsOutWrn=1; PartsOutTotWrn=10;

  SvData=byte(SDAT_Binx)|byte(SDAT_Info);
  SvRes=false;
  SvTimers=false;
  SvDomainVtk=false;

  KernelH=CteB=Gamma=RhopZero=CFLnumber=0;
  Dp=0;
  MassFluid=MassBound=0;
  Gravity=TFloat3(0);

  KernelSize=KernelSize2=0;
  Cs0=0;
  Eta2=0;
  SpsSmag=SpsBlin=0;
  DDTkh=DDTgz=0;
  memset(&CSP,0,sizeof(StCteSph));

  CasePosMin=CasePosMax=TDouble3(0);
  CaseNp=CaseNbound=CaseNfixed=CaseNmoving=CaseNfloat=CaseNfluid=CaseNpb=0;

  PeriActive=0; PeriX=PeriY=PeriZ=false;
  PeriXinc=PeriYinc=PeriZinc=TDouble3(0);

  PartBeginDir=""; 
  PartBegin=PartBeginFirst=0;
  PartBeginTimeStep=0; 
  PartBeginTotalNp=0;

  WrnPartsOut=true;

  FtCount=0;
  FtPause=0;
  FtMode=FTMODE_Sph;
  FtConstraints=false;
  FtIgnoreRadius=false;
  WithFloating=false;

  AllocMemoryFloating(0,false);

  CellDomFixed=false;
  CellMode=CELLMODE_None;
  ScellDiv=0;
  Scell=0;
  MovLimit=0;

  PosCellSize=0;

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

  PartIni=Part=0; 
  Nstep=0; PartNstep=-1;
  PartOut=0;

  TimeStepIni=0;
  TimeStep=TimeStepM1=0;
  TimePartNext=0;
  LastDt=0;

  VerletStep=0;
  SymplecticDtPre=0;
  DemDtForce=0;  //(DEM)
  MaxNumbers.Clear();

  SaveFtAce=false;
}

//==============================================================================
/// Saves the linear and angular acceleration of each floating object in a csv file
//==============================================================================
void JSph::SaveFtAceFun(double dt,bool predictor,StFtoForces *ftoforces){
  const unsigned nstep=Nstep;
  const double timestep=TimeStep;
  const bool savedata=!(timestep<(floor((timestep-dt)/TimePart)+1.0)*TimePart);
  if(savedata && !predictor){
    for(unsigned cf=0;cf<FtCount;cf++){
      const string file=AppInfo.GetDirOut()+fun::PrintStr("FloatingAce_mkbound_%u.csv",FtObjs[cf].mkbound);
      jcsv::JSaveCsv2 scsv(file,true,AppInfo.GetCsvSepComa());
      if(!scsv.GetAppendMode()){
        Log->AddFileInfo("FloatingAce_mkbound_XX.csv", "Saves information of acceleration used to move the floating object.");
        //-Saves head.
        scsv.SetHead();
        scsv << "nstep;time [s];dt [s];predictor;face.x [m/s^2];face.y [m/s^2];face.z [m/s^2]";
        scsv << "fomegaace.x [rad/s^2];fomegaace.y [rad/s^2];fomegaace.z [rad/s^2]" << jcsv::Endl();
      }
      //-Saves data.
      scsv.SetData();
      scsv << nstep << timestep << dt << (predictor?"True":"False");
      scsv << ftoforces[cf].face;
      scsv << ftoforces[cf].fomegaace;
      scsv << jcsv::Endl();
      scsv.SaveData();
    }
  }
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
  const string keyend=(key.size()>=4? key.substr(key.size()-4,4): "");
       if(keyend=="Xmin")CfgDomainFixedMin.x=v;
  else if(keyend=="Ymin")CfgDomainFixedMin.y=v;
  else if(keyend=="Zmin")CfgDomainFixedMin.z=v;
  else if(keyend=="Xmax")CfgDomainFixedMax.x=v;
  else if(keyend=="Ymax")CfgDomainFixedMax.y=v;
  else if(keyend=="Zmax")CfgDomainFixedMax.z=v;
  else Run_Exceptioon("Key for limit is invalid.");
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
  const string keyend=(key.size()>=4? key.substr(key.size()-4,4): "");
       if(keyend=="Xmin")CfgDomainParticlesMin.x=v;
  else if(keyend=="Ymin")CfgDomainParticlesMin.y=v;
  else if(keyend=="Zmin")CfgDomainParticlesMin.z=v;
  else if(keyend=="Xmax")CfgDomainParticlesMax.x=v;
  else if(keyend=="Ymax")CfgDomainParticlesMax.y=v;
  else if(keyend=="Zmax")CfgDomainParticlesMax.z=v;
  else Run_Exceptioon("Key for limit is invalid.");
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
  const string keyend=(key.size()>=4? key.substr(key.size()-4,4): "");
       if(keyend=="Xmin")CfgDomainParticlesPrcMin.x=v;
  else if(keyend=="Ymin")CfgDomainParticlesPrcMin.y=v;
  else if(keyend=="Zmin")CfgDomainParticlesPrcMin.z=v;
  else if(keyend=="Xmax")CfgDomainParticlesPrcMax.x=v;
  else if(keyend=="Ymax")CfgDomainParticlesPrcMax.y=v;
  else if(keyend=="Zmax")CfgDomainParticlesPrcMax.z=v;
  else Run_Exceptioon("Key for limit is invalid.");
}

//==============================================================================
/// Loads the case configuration to be executed.
//==============================================================================
void JSph::ConfigDomainResize(std::string key,const JCaseEParms *eparms){
  const char axis=fun::StrLower(key)[0];
  if(axis!='x' && axis!='y' && axis!='z')Run_Exceptioon("Axis value is invalid.");
  if(key.substr(1,3)!="min" && key.substr(1,3)!="max")Run_Exceptioon("Key value is invalid.");
  if(key.substr(1,3)=="min"){
    JCaseEParms::JCaseEParmsPos ps=eparms->GetPosminValue(axis);
    switch(ps.mode){
      case JCaseEParms::DC_Fixed:     ConfigDomainFixedValue(string("DomainFixed")+key,ps.value);                     break;
      case JCaseEParms::DC_DefValue:  ConfigDomainParticlesValue(string("DomainParticles")+key,ps.value);             break;
      case JCaseEParms::DC_DefPrc:    ConfigDomainParticlesPrcValue(string("DomainParticlesPrc")+key,ps.value/100);   break;
    }
  }
  else{
    JCaseEParms::JCaseEParmsPos ps=eparms->GetPosmaxValue(axis);
    switch(ps.mode){
      case JCaseEParms::DC_Fixed:     ConfigDomainFixedValue(string("DomainFixed")+key,ps.value);                     break;
      case JCaseEParms::DC_DefValue:  ConfigDomainParticlesValue(string("DomainParticles")+key,ps.value);             break;
      case JCaseEParms::DC_DefPrc:    ConfigDomainParticlesPrcValue(string("DomainParticlesPrc")+key,ps.value/100);   break;
    }
  }
}

//==============================================================================
/// Allocates memory of floating objectcs.
//==============================================================================
void JSph::AllocMemoryFloating(unsigned ftcount,bool imposedvel,bool addedforce){
  delete[] FtObjs; FtObjs=NULL;
  if(FtLinearVel){
    for(unsigned c=0;c<FtCount;c++)delete FtLinearVel[c];
    delete[] FtLinearVel; FtLinearVel=NULL;
  }
  if(FtAngularVel){
    for(unsigned c=0;c<FtCount;c++)delete FtAngularVel[c];
    delete[] FtAngularVel; FtAngularVel=NULL;
  } 
  if(FtLinearForce){
    for(unsigned c=0;c<FtCount;c++)delete FtLinearForce[c];
    delete[] FtLinearForce; FtLinearForce=NULL;
  }
  if(FtAngularForce){
    for(unsigned c=0;c<FtCount;c++)delete FtAngularForce[c];
    delete[] FtAngularForce; FtAngularForce=NULL;
  }
  if(ftcount){
    FtObjs=new StFloatingData[ftcount];
    if(imposedvel){
      FtLinearVel   =new JLinearValue*[ftcount];
      FtAngularVel  =new JLinearValue*[ftcount];
      for(unsigned c=0;c<ftcount;c++)FtLinearVel[c]=FtAngularVel[c]=NULL;
    } 
    if(addedforce){
      FtLinearForce =new JLinearValue*[ftcount];
      FtAngularForce=new JLinearValue*[ftcount];
      for(unsigned c=0;c<ftcount;c++)FtLinearForce[c]=FtAngularForce[c]=NULL;
    }
  }
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
  if(FixedDt)s+=FixedDt->GetAllocMemory();
  if(AccInput)s+=AccInput->GetAllocMemory();
  if(PartsLoaded)s+=PartsLoaded->GetAllocMemory();
  return(s);
}

//==============================================================================
/// Loads the configuration of the execution.
//==============================================================================
void JSph::LoadConfig(const JSphCfgRun *cfg){
  TimerTot.Start();

  //-Loads basic configuration from execution parameters.
  //------------------------------------------------------
  Stable=cfg->Stable;
  SvPosDouble=false; //-Options by default.
  DirOut=fun::GetDirWithSlash(cfg->DirOut);
  DirDataOut=(!cfg->DirDataOut.empty()? fun::GetDirWithSlash(DirOut+cfg->DirDataOut): DirOut);
  CaseName=cfg->CaseName; 
  DirCase=fun::GetDirWithSlash(fun::GetDirParent(CaseName));
  CaseName=CaseName.substr(DirCase.length());
  if(!CaseName.length())Run_Exceptioon("Name of the case for execution was not indicated.");
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
  if(!JVtkLib::Available())Log->PrintWarning("Code for VTK format files is not included in the current compilation, so no output VTK files will be created.");
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

  //-Loads case configuration from XML and command line.
  LoadCaseConfig(cfg);

  //-PIPS configuration.
  if(cfg->PipsMode)DsPips=new JDsPips(Cpu,cfg->PipsSteps,(cfg->PipsMode==2),(TStep==STEP_Symplectic? 2: 1));
}

//==============================================================================
/// Loads kernel selection to compute kernel values.
//==============================================================================
void JSph::LoadKernelSelection(const JSphCfgRun *cfg,const JXml *xml){
  //-Load kernel selection from execution parameters from XML.
  JCaseEParms eparms;
  eparms.LoadXml(xml,"case.execution.parameters");
  switch(eparms.GetValueInt("Kernel",true,2)){
    case 1:  TKernel=KERNEL_Cubic;     break;
    case 2:  TKernel=KERNEL_Wendland;  break;
    default: Run_Exceptioon("Kernel choice is not valid.");
  }
  //-Load kernel selection from execution parameters from commands.
  if(cfg->TKernel)TKernel=cfg->TKernel;
}

//==============================================================================
/// Loads predefined constans from XML.
//==============================================================================
void JSph::LoadConfigCtes(const JXml *xml){
  JCaseCtes ctes;
  ctes.LoadXmlRun(xml,"case.execution.constants");

  Simulate2D=ctes.GetData2D();
  Simulate2DPosY=ctes.GetData2DPosY();
  KernelH=(float)ctes.GetH();
  CteB=(float)ctes.GetB();
  Gamma=(float)ctes.GetGamma();
  RhopZero=(float)ctes.GetRhop0();
  CFLnumber=(float)ctes.GetCFLnumber();
  Dp=ctes.GetDp();
  MassFluid=(float)ctes.GetMassFluid();
  MassBound=(float)ctes.GetMassBound();
  Gravity=ToTFloat3(ctes.GetGravity());
  if(ctes.GetEps()!=0)Log->PrintWarning("Eps value is not used (this correction was removed).");
}

//==============================================================================
/// Loads execution parameters from XML.
//==============================================================================
void JSph::LoadConfigParameters(const JXml *xml){
  JCaseEParms eparms;
  eparms.LoadXml(xml,"case.execution.parameters");
  if(eparms.Exists("FtSaveAce"))SaveFtAce=(eparms.GetValueInt("FtSaveAce",true,0)!=0); //-For Debug.
  if(eparms.GetValueDouble("FtSaveMotion",true,-1.)>=0)FtMotSave=new JFtMotionSave(eparms.GetValueDouble("FtSaveMotion")); //<vs_ftmottionsv>
  if(eparms.Exists("PosDouble")){
    Log->PrintWarning("The parameter \'PosDouble\' is deprecated.");
    SvPosDouble=(eparms.GetValueInt("PosDouble")==2);
  }
  if(eparms.Exists("SavePosDouble"))SvPosDouble=(eparms.GetValueInt("SavePosDouble",true,0)!=0);
  switch(eparms.GetValueInt("RigidAlgorithm",true,1)){ //(DEM)
    case 0:  FtMode=FTMODE_Ext;                   break;
    case 1:  FtMode=FTMODE_Sph;                   break;
    case 2:  FtMode=FTMODE_Ext;  UseDEM=true;     break;
    case 3:  FtMode=FTMODE_Ext;  UseChrono=true;  break;
    default: Run_Exceptioon("Rigid algorithm is not valid.");
  }
  switch(eparms.GetValueInt("StepAlgorithm",true,1)){
    case 1:  TStep=STEP_Verlet;      break;
    case 2:  TStep=STEP_Symplectic;  break;
    default: Run_Exceptioon("Step algorithm is not valid.");
  }
  VerletSteps=eparms.GetValueInt("VerletSteps",true,40);
  switch(eparms.GetValueInt("ViscoTreatment",true,1)){
    case 1:  TVisco=VISCO_Artificial;  break;
    case 2:  TVisco=VISCO_LaminarSPS;  break;
    default: Run_Exceptioon("Viscosity treatment is not valid.");
  }
  Visco=eparms.GetValueFloat("Visco");
  ViscoBoundFactor=eparms.GetValueFloat("ViscoBoundFactor",true,1.f);
  string filevisco=eparms.GetValueStr("ViscoTime",true);
  if(!filevisco.empty()){
    ViscoTime=new JDsViscoInput();
    ViscoTime->LoadFile(DirCase+filevisco);
  }

  //-Boundary configuration.
  switch(eparms.GetValueInt("Boundary",true,1)){
    case 1:  TBoundary=BC_DBC;      break;
    case 2:  TBoundary=BC_MDBC;     break;
    default: Run_Exceptioon("Boundary Condition method is not valid.");
  }
  if(TBoundary==BC_MDBC){
    UseNormals=true;
    switch(eparms.GetValueInt("SlipMode",true,1)){
      case 1:  SlipMode=SLIP_Vel0;      break;
      case 2:  SlipMode=SLIP_NoSlip;    break;
      case 3:  SlipMode=SLIP_FreeSlip;  break;
      default: Run_Exceptioon("Slip mode is not valid.");
    }
    MdbcCorrector=(eparms.GetValueInt("MDBCCorrector",true,0)!=0);
    MdbcFastSingle=(eparms.GetValueInt("MDBCFastSingle",true,1)!=0);
    if(Cpu)MdbcFastSingle=false;
  } 

  //-Density Diffusion Term configuration.
  if(eparms.Exists("DeltaSPH")){
    if(eparms.Exists("DensityDT"))Run_Exceptioon("The parameters \'DeltaSPH\' and \'DensityDT\' cannot be combined. Only \'DensityDT\' should be used.");
    float v=eparms.GetValueFloat("DeltaSPH");
    if(v>0){
      TDensity=DDT_DDT;
      DDTValue=v;
    }
    Log->PrintWarning("The parameter \'DeltaSPH\' is deprecated. It should be replaced by \'DensityDT\'.");
  }
  else{
    const int tddt=eparms.GetValueInt("DensityDT",true,0);
    switch(tddt){
      case 0:  TDensity=DDT_None;     break;
      case 1:  TDensity=DDT_DDT;      break;
      case 2:  TDensity=DDT_DDT2;     break;
      case 3:  TDensity=DDT_DDT2Full; break;
      default: Run_Exceptioon("Density Diffusion Term mode is not valid.");
    }
    DDTValue=eparms.GetValueFloat("DensityDTvalue",true,0.1f);
  }
  DDTArray=false;

  //-Old shifting configuration.
  if(eparms.Exists("Shifting")){
    TpShifting shiftmode=SHIFT_None;
    float shiftcoef=0,shifttfs=0;
    switch(eparms.GetValueInt("Shifting",true,0)){
      case 0:  shiftmode=SHIFT_None;     break;
      case 1:  shiftmode=SHIFT_NoBound;  break;
      case 2:  shiftmode=SHIFT_NoFixed;  break;
      case 3:  shiftmode=SHIFT_Full;     break;
      default: Run_ExceptioonFile("Shifting mode in <execution><parameters> is not valid.",FileXml);
    }
    if(shiftmode!=SHIFT_None){
      shiftcoef=eparms.GetValueFloat("ShiftCoef",true,-2);
      if(shiftcoef==0)shiftmode=SHIFT_None;
      else shifttfs=eparms.GetValueFloat("ShiftTFS",true,0);
    }
    Shifting=new JSphShifting(Simulate2D,Dp,KernelH);
    Shifting->ConfigBasic(shiftmode,shiftcoef,shifttfs);
  }

  WrnPartsOut=(eparms.GetValueInt("WrnPartsOut",true,1)!=0);
  FtPause=eparms.GetValueFloat("FtPause",true,0);
  FtIgnoreRadius=(eparms.GetValueInt("FtIgnoreRadius",true,0)!=0);

  TimeMax=eparms.GetValueDouble("TimeMax");
  TimePart=eparms.GetValueDouble("TimeOut");

  DtIni=max(0.0,eparms.GetValueDouble("DtIni",true,0));
  DtMin=max(0.0,eparms.GetValueDouble("DtMin",true,0));
  CoefDtMin=eparms.GetValueFloat("CoefDtMin",true,0.05f);
  DtAllParticles=(eparms.GetValueInt("DtAllParticles",true,0)==1);

  double valuefixeddt=max(0.0,eparms.GetValueDouble("DtFixed",true,0));
  string filefixeddt=eparms.GetValueStr("DtFixedFile",true);
  if(fun::StrUpper(filefixeddt)=="NONE")filefixeddt="";
  if(valuefixeddt && !filefixeddt.empty())Run_Exceptioon("The parameters \'DtFixed\' and \'DtFixedFile\' cannot be used at the same time.");
  if(!filefixeddt.empty()){
    FixedDt=new JDsFixedDt();
    FixedDt->LoadFile(DirCase+filefixeddt);
  }
  else if(valuefixeddt)FixedDt=new JDsFixedDt(valuefixeddt);

  if(eparms.Exists("RhopOutMin"))RhopOutMin=eparms.GetValueFloat("RhopOutMin");
  if(eparms.Exists("RhopOutMax"))RhopOutMax=eparms.GetValueFloat("RhopOutMax");
  PartsOutMax=eparms.GetValueFloat("PartsOutMax",true,1);

  //-Configuration of symmetry calculation.            //<vs_syymmetry>
  Symmetry=(eparms.GetValueInt("Symmetry",true,0)!=0); //<vs_syymmetry>

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
    if(Simulate2D && PeriY)Run_Exceptioon("Cannot use periodic conditions in Y with 2D simulations");
  }

  //-Configuration of domain size.
  bool resizeold=false;
  float incz=eparms.GetValueFloat("IncZ",true,0.f);
  if(incz){
    ClearCfgDomain();
    CfgDomainParticlesPrcMax.z=incz;
    resizeold=true;
  }
  string key;
  if(eparms.Exists(key="DomainFixed")){ ConfigDomainFixed(TDouble3(eparms.GetValueNumDouble(key,0),eparms.GetValueNumDouble(key,1),eparms.GetValueNumDouble(key,2)),TDouble3(eparms.GetValueNumDouble(key,3),eparms.GetValueNumDouble(key,4),eparms.GetValueNumDouble(key,5))); resizeold=true; }
  if(eparms.Exists(key="DomainFixedXmin")){ ConfigDomainFixedValue(key,eparms.GetValueDouble(key)); resizeold=true; }
  if(eparms.Exists(key="DomainFixedYmin")){ ConfigDomainFixedValue(key,eparms.GetValueDouble(key)); resizeold=true; }
  if(eparms.Exists(key="DomainFixedZmin")){ ConfigDomainFixedValue(key,eparms.GetValueDouble(key)); resizeold=true; }
  if(eparms.Exists(key="DomainFixedXmax")){ ConfigDomainFixedValue(key,eparms.GetValueDouble(key)); resizeold=true; }
  if(eparms.Exists(key="DomainFixedYmax")){ ConfigDomainFixedValue(key,eparms.GetValueDouble(key)); resizeold=true; }
  if(eparms.Exists(key="DomainFixedZmax")){ ConfigDomainFixedValue(key,eparms.GetValueDouble(key)); resizeold=true; }
  if(!eparms.IsPosDefault() && resizeold)Run_ExceptioonFile("Combination of <simulationdomain> with IncZ or DomainFixedXXX in <parameters> section of XML is not allowed.",FileXml);
  if(resizeold)Log->PrintWarning("The options IncZ and DomainFixedXXXX are deprecated.");
  ConfigDomainResize("Xmin",&eparms);
  ConfigDomainResize("Ymin",&eparms);
  ConfigDomainResize("Zmin",&eparms);
  ConfigDomainResize("Xmax",&eparms);
  ConfigDomainResize("Ymax",&eparms);
  ConfigDomainResize("Zmax",&eparms);
}

//==============================================================================
/// Loads the case configuration to be executed.
//==============================================================================
void JSph::LoadConfigCommands(const JSphCfgRun *cfg){
  //-Aplies configuration using command line.
  if(cfg->SvPosDouble>=0)SvPosDouble=(cfg->SvPosDouble!=0);
  if(cfg->TBoundary){
    TBoundary=BC_DBC;
    SlipMode=SLIP_Vel0;
    MdbcCorrector=false;
    MdbcFastSingle=true;
    MdbcThreshold=0;
    switch(cfg->TBoundary){
      case 1:  TBoundary=BC_DBC;   break;
      case 2:  TBoundary=BC_MDBC;  break;
      default: Run_Exceptioon("Boundary method is not valid.");
    }
    if(TBoundary==BC_MDBC)switch(cfg->SlipMode){
      case 1:  SlipMode=SLIP_Vel0;      break;
      case 2:  SlipMode=SLIP_NoSlip;    break;
      case 3:  SlipMode=SLIP_FreeSlip;  break;
      default: Run_Exceptioon("Slip mode for mDBC is not valid.");
    }
    UseNormals=(TBoundary==BC_MDBC);
  }
  if(TBoundary==BC_MDBC){
    if(cfg->MdbcThreshold >=0)MdbcThreshold=cfg->MdbcThreshold;
    if(cfg->MdbcFastSingle>=0)MdbcFastSingle=(cfg->MdbcFastSingle>0);
    if(SlipMode!=SLIP_Vel0)Run_Exceptioon("Only the slip mode velocity=0 is allowed with mDBC conditions."); //SHABA
    if(Cpu)MdbcFastSingle=false;
  }
    
  if(cfg->TStep)TStep=cfg->TStep;
  if(cfg->VerletSteps>=0)VerletSteps=cfg->VerletSteps;
  if(cfg->TVisco){ TVisco=cfg->TVisco; Visco=cfg->Visco; }
  if(cfg->ViscoBoundFactor>=0)ViscoBoundFactor=cfg->ViscoBoundFactor;

  //-Density Diffusion Term configuration.
  if(cfg->TDensity>=0){
    switch(cfg->TDensity){
      case 0:  TDensity=DDT_None;     break;
      case 1:  TDensity=DDT_DDT;      break;
      case 2:  TDensity=DDT_DDT2;     break;
      case 3:  TDensity=DDT_DDT2Full; break;
      default: Run_Exceptioon("Density Diffusion Term mode is not valid.");
    }
    if(cfg->DDTValue>=0)DDTValue=cfg->DDTValue;
    else if(DDTValue==0)DDTValue=0.1f;
  }
  if(TDensity==DDT_None)DDTValue=0;
  else if(cfg->DDTValue>=0)DDTValue=cfg->DDTValue;
  DDTArray=(TDensity!=DDT_None && Cpu); //-It is necessary because the interaction is divided in two steps: fluid-fluid/float and fluid-bound.
  //-Checks warnings of DDT according to gravity.
  if(Gravity.z==0){
    if(TDensity==DDT_DDT2)Log->PrintWarning("The option DDT:2 is equal to DDT:1 when gravity.z is zero.");
    if(TDensity==DDT_DDT2Full)Log->PrintWarning("The option DDT:3 is equal to DDT:1 (but applied to the entire fluid) when gravity.z is zero.");
  }
  if((Gravity.x!=0 || Gravity.y!=0) && (TDensity==DDT_DDT2 || TDensity==DDT_DDT2Full)){
    Log->PrintWarning("Gravity.x or Gravity.y is not zero, but only gravity.z is used in DDT:2 or DDT:3 calculations.");
  }
  
  //-Shifting configuration.
  if(cfg->Shifting>=0){
    TpShifting shiftmode=SHIFT_None;
    switch(cfg->Shifting){
      case 0:  shiftmode=SHIFT_None;     break;
      case 1:  shiftmode=SHIFT_NoBound;  break;
      case 2:  shiftmode=SHIFT_NoFixed;  break;
      case 3:  shiftmode=SHIFT_Full;     break;
      default: Run_Exceptioon("Shifting mode is not valid.");
    }
    if(!Shifting)Shifting=new JSphShifting(Simulate2D,Dp,KernelH);
    Shifting->ConfigBasic(shiftmode);
  }

  if(cfg->FtPause>=0)FtPause=cfg->FtPause;
  if(cfg->TimeMax>0)TimeMax=cfg->TimeMax;
  NstepsBreak=cfg->NstepsBreak;
  if(NstepsBreak)Log->PrintfWarning("The execution will be cancelled after %d simulation steps.",NstepsBreak);
  SvAllSteps=cfg->SvAllSteps;
  //-Configuration of JDsOutputTime with TimePart.
  OutputTime=new JDsOutputTime();
  if(cfg->TimePart>=0){
    TimePart=cfg->TimePart;
    OutputTime->Config(TimePart);
  }
  else OutputTime->Config(FileXml,"case.execution.special.timeout",TimePart);

  //-Configuration of domain limits.
  CellDomFixed=cfg->CellDomFixed;
  CellMode=cfg->CellMode;
  if(cfg->DomainMode==2)ConfigDomainFixed(cfg->DomainFixedMin,cfg->DomainFixedMax);

  //-Configuration of density limits.
  if(cfg->RhopOutModif){
    RhopOutMin=cfg->RhopOutMin; RhopOutMax=cfg->RhopOutMax;
  }
  RhopOut=(RhopOutMin<RhopOutMax);
  if(!RhopOut){ RhopOutMin=-FLT_MAX; RhopOutMax=FLT_MAX; }
  if(RhopZero<RhopOutMin || RhopZero>RhopOutMax)
    Run_Exceptioon(fun::PrintStr("The reference density value %f is outside the defined limits [%f,%f].",RhopZero,RhopOutMin,RhopOutMax));
}

//==============================================================================
/// Creates and load NuxLib object to evaluate user-defined expressions.
//==============================================================================
void JSph::LoadConfigVars(const JXml *xml){
  if(!JNumexLib::Available())Log->PrintWarning("Code for JNumex libary is not included in the current compilation, so user-defined expresions in XML file are not evaluated.");
  else{
    NuxLib=new JNumexLib();
    //-Loads user variables from XML.
    JCaseUserVars uvars;
    uvars.LoadXml(xml,"case.execution.uservars",true);
    for(unsigned c=0;c<uvars.CountVars();c++){
      const JCaseUserVars::StVar v=uvars.GetVar(c);
      if(v.isnum)NuxLib->CreateVar(v.name,true,false,v.valuenum);
      else       NuxLib->CreateVar(v.name,true,false,v.valuestr);
    }
    //-Defines basic constant values on NuxLib object.
    bool rep=false;
    rep|=NuxLib->CreateVar("CaseName"  ,true,true,CaseName);
    rep|=NuxLib->CreateVar("Data2D"    ,true,true,Simulate2D);
    rep|=NuxLib->CreateVar("Data2DPosy",true,true,Simulate2DPosY);
    rep|=NuxLib->CreateVar("H"         ,true,true,KernelH);
    rep|=NuxLib->CreateVar("KernelSize",true,true,KernelSize);
    rep|=NuxLib->CreateVar("B"         ,true,true,CteB);
    rep|=NuxLib->CreateVar("Gamma"     ,true,true,Gamma);
    rep|=NuxLib->CreateVar("Rhop0"     ,true,true,RhopZero);
    rep|=NuxLib->CreateVar("Dp"        ,true,true,Dp);
    rep|=NuxLib->CreateVar("Gravity"   ,true,true,Gravity);
    rep|=NuxLib->CreateVar("MassFluid" ,true,true,MassFluid);
    rep|=NuxLib->CreateVar("MassBound" ,true,true,MassBound);
    if(rep)Log->PrintWarning("Some user-defined variable from XML was replaced by values of constants.");
    //-Shows list of available variables.
    //Log->Printf("List of available variables for user expressions: %s\n",NuxLib->ListVarsToStr().c_str());
    Log->Printf("XML-Vars (uservars + ctes): %s",NuxLib->ListVarsToStr().c_str());
  }
}

//==============================================================================
/// Loads execution parameters in NuxLib object to evaluate user-defined expressions.
//==============================================================================
void JSph::LoadConfigVarsExec(){
  if(NuxLib){
    const unsigned nv=NuxLib->CountVars();
    bool rep=false;
    rep|=NuxLib->CreateVar("TimeMax",true,true,TimeMax);
    rep|=NuxLib->CreateVar("TimeOut",true,true,TimePart);
    if(rep)Log->PrintWarning("Some user-defined variable from XML was replaced by values of parameters.");
    //-Shows list of available variables.
    //Log->Printf("List of available variables for user expressions: %s\n",NuxLib->ListVarsToStr().c_str());
    Log->Printf("XML-Vars (parameters): %s",NuxLib->ListVarsToStr(nv).c_str());
  }
}

//==============================================================================
/// Loads the case configuration to be executed.
//==============================================================================
void JSph::LoadCaseConfig(const JSphCfgRun *cfg){
  if(!fun::FileExists(FileXml))Run_ExceptioonFile("Case configuration was not found.",FileXml);
  JXml xml; xml.LoadFile(FileXml);
  //-Shows pre-processing application generating the XML file.
  Log->Printf("XML-App: %s",xml.GetAttributeStr(xml.GetNodeError("case")->ToElement(),"app",true,"unknown").c_str());

  //-Loads kernel selection to compute kernel values.
  LoadKernelSelection(cfg,&xml);
  //-Loads predefined constants.
  LoadConfigCtes(&xml);
  //-Configures value of calculated constants and loads CSP structure.
  ConfigConstants1(Simulate2D);

  //-Define variables on NuxLib object.
  LoadConfigVars(&xml);
  //-Enables the use of NuxLib in XML configuration.
  xml.SetNuxLib(NuxLib);

  //-Execution parameters from XML.
  LoadConfigParameters(&xml);
  //-Execution parameters from commands.
  LoadConfigCommands(cfg);
  //-Configures other constants according to formulation options and loads more values in CSP structure.
  ConfigConstants2();

  //-Define variables of execution parameters and shows list of available variables.
  LoadConfigVarsExec();

  //-Particle data.
  JCaseParts parts;
  parts.LoadXml(&xml,"case.execution.particles");
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
  GaugeSystem=new JGaugeSystem(Cpu);

  //-Configuration of AccInput.
  if(xml.GetNodeSimple("case.execution.special.accinputs",true)){
    AccInput=new JDsAccInput(DirCase,&xml,"case.execution.special.accinputs");
  }

  //-Configuration of ChronoObjects.
  if(UseChrono){
    if(xml.GetNodeSimple("case.execution.special.chrono",true)){
      if(!JChronoObjects::Available())Run_Exceptioon("DSPHChronoLib to use Chrono is not included in the current compilation.");
      ChronoObjects=new JChronoObjects(DirCase,CaseName,&xml,"case.execution.special.chrono",Dp,parts.GetMkBoundFirst(),Gravity,Simulate2D,FtPause);
    }
    else Run_ExceptioonFile("Chrono configuration in XML file is missing.",FileXml);
  }
  else if(xml.GetNodeSimple("case.execution.special.chrono",true))Log->PrintfWarning("The use of Chrono is disabled (RigidAlgorithm is not 3) although XML file includes the Chrono configuration.");

  //-Loads and configures moving objects.
  if(parts.CountBlocks(TpPartMoving)>0){
    DsMotion=new JDsMotion(Simulate2D);
    DsMotion->Init(&parts,&xml,"case.execution.motion",DirCase);
  }

  //-Configuration of WaveGen.
  if(xml.GetNodeSimple("case.execution.special.wavepaddles",true)){
    if(!JWaveGen::Available())Run_Exceptioon("Code for Wave-Paddles is not included in the current compilation.");
    bool useomp=false,usegpu=false;
    #ifdef OMP_USE_WAVEGEN
      useomp=(omp_get_max_threads()>1);
    #endif
    WaveGen=new JWaveGen(useomp,!Cpu,Log,DirCase,&xml,"case.execution.special.wavepaddles",ToTDouble3(Gravity));
    if(DsMotion)for(unsigned ref=0;ref<DsMotion->GetNumObjects();ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      WaveGen->ConfigPaddle(m.mkbound,ref,m.idbegin,m.count);
    }
  }

  //-Configuration of MLPistons.
  if(xml.GetNodeSimple("case.execution.special.mlayerpistons",true)){
    if(!JMLPistons::Available())Run_Exceptioon("Code for Multi-Layer Pistons is not included in the current compilation.");
    bool useomp=false,usegpu=false;
    MLPistons=new JMLPistons(!Cpu,Log,DirCase);
    MLPistons->LoadXml(&xml,"case.execution.special.mlayerpistons");  
    if(DsMotion)for(unsigned c=0;c<DsMotion->GetNumObjects();c++){
      const StMotionData& m=DsMotion->GetMotionData(c);
      MLPistons->ConfigPiston(m.mkbound,c,m.idbegin,m.count,TimeMax);
    }
    MLPistons->CheckPistons();
  }

  //-Configuration of RelaxZones.
  if(xml.GetNodeSimple("case.execution.special.relaxationzones",true)){
    if(!JRelaxZones::Available())Run_Exceptioon("Code for Relaxation Zones is not included in the current compilation.");
    bool useomp=false,usegpu=false;
    #ifdef OMP_USE_WAVEGEN
      useomp=(omp_get_max_threads()>1);
    #endif
    RelaxZones=new JRelaxZones(useomp,!Cpu,Log,DirCase,CaseNfloat>0,CaseNbound,ToTDouble3(Gravity));
    RelaxZones->LoadXml(&xml,"case.execution.special.relaxationzones");
  }

  //-Configuration of Shifting with zones.
  if(Shifting && xml.GetNodeSimple("case.execution.special.shifting",true))Run_ExceptioonFile("Shifting is defined several times (in <special><shifting> and <execution><parameters>).",FileXml);
  if(xml.GetNodeSimple("case.execution.special.shifting",true)){
    Shifting=new JSphShifting(Simulate2D,Dp,KernelH);
    Shifting->LoadXml(&xml,"case.execution.special.shifting");
  }
  if(Shifting && !Shifting->GetShiftMode()){ delete Shifting; Shifting=NULL; }
  ShiftingMode=(Shifting? Shifting->GetShiftMode(): SHIFT_None);

  //-Configuration of damping zones.
  if(xml.GetNodeSimple("case.execution.special.damping",true)){
    Damping=new JDsDamping(Dp);
    Damping->LoadXml(&xml,"case.execution.special.damping");
  }

  //-Loads floating objects.
  FtCount=parts.CountBlocks(TpPartFloating);
  if(FtCount){
    FtConstraints=false;
    if(FtCount>CODE_MKRANGEMAX)Run_Exceptioon("The number of floating objects exceeds the maximum.");
    AllocMemoryFloating(FtCount,parts.UseImposedFtVel(),parts.UseAddedFtForce());
    unsigned cobj=0;
    for(unsigned c=0;c<parts.CountBlocks()&&cobj<FtCount;c++){
      const JCasePartBlock &block=parts.GetBlock(c);
      if(block.Type==TpPartFloating){
        const JCasePartBlock_Floating &fblock=(const JCasePartBlock_Floating &)block;
        StFloatingData* fobj=FtObjs+cobj;
        fobj->mkbound=fblock.GetMkType();
        fobj->begin=fblock.GetBegin();
        fobj->count=fblock.GetCount();
        fobj->mass=(float)fblock.GetMassbody();
        fobj->massp=(float)fblock.GetMasspart();
        fobj->radius=0;
        fobj->constraints=ComputeConstraintsValue(fblock.GetTranslationFree(),fblock.GetRotationFree());
        if(fobj->constraints!=FTCON_Free)FtConstraints=true;
        if(fblock.GetLinearVel()){
          FtLinearVel[cobj]=new JLinearValue(*fblock.GetLinearVel());
          if(!FtLinearVel[cobj]->GetFile().empty())FtLinearVel[cobj]->LoadFile(DirCase+FtLinearVel[cobj]->GetFile());
        }
        if(fblock.GetAngularVel()){
          FtAngularVel[cobj]=new JLinearValue(*fblock.GetAngularVel());
          if(!FtAngularVel[cobj]->GetFile().empty())FtAngularVel[cobj]->LoadFile(DirCase+FtAngularVel[cobj]->GetFile());
        } 
        if(fblock.GetLinearForce()){
          FtLinearForce[cobj]=new JLinearValue(*fblock.GetLinearForce());
          if(!FtLinearForce[cobj]->GetFile().empty())FtLinearForce[cobj]->LoadFile(DirCase+FtLinearForce[cobj]->GetFile());
        } 
        if(fblock.GetAngularForce()){
          FtAngularForce[cobj]=new JLinearValue(*fblock.GetAngularForce());
          if(!FtAngularForce[cobj]->GetFile().empty())FtAngularForce[cobj]->LoadFile(DirCase+FtAngularForce[cobj]->GetFile());
        } 
        fobj->center=fblock.GetCenter();
        fobj->angles=TFloat3(0);
        fobj->fvel=ToTFloat3(fblock.GetLinearVelini());
        fobj->fomega=ToTFloat3(fblock.GetAngularVelini());
        fobj->facelin=fobj->faceang=TFloat3(0);
        tmatrix3d inertiaini=fblock.GetInertia();
        //-Set to zero the values Ixx and Izz in the inertia matrix when it is a 2D-Simulation.
        //if(Simulate2D)inertiaini=TMatrix3d(0,0,0,0,inertiaini.a22,0,0,0,0);
        fobj->inertiaini=ToTMatrix3f(inertiaini);
        //-Chrono configuration.
        fobj->usechrono=(ChronoObjects && ChronoObjects->ConfigBodyFloating(fblock.GetMkType()
          ,fblock.GetMassbody(),fblock.GetCenter(),inertiaini,fblock.GetTranslationFree()
          ,fblock.GetRotationFree(),fobj->fvel,fobj->fomega));
        cobj++;
      }
    }
  }
  else{
    if(UseDEM){
      UseDEM=false;
      Log->PrintWarning("The use of DEM was disabled because there are no floating objects...");
    }
    if(UseChrono){
      UseChrono=false;
      Log->PrintWarning("The use of CHRONO was disabled because there are no floating objects...");
      delete ChronoObjects; ChronoObjects=NULL;
    }
  }
  WithFloating=(FtCount>0);
  if(!WithFloating)FtMode=FTMODE_None; //-Disables floatings when there are not floating particles.
  if(UseChrono && PeriActive!=0)Log->PrintfWarning("The use of Chrono with periodic limits is only recommended when moving and floating objects do not move beyond those periodic limits.");

  //-Loads DEM and DVI data for boundary objects.
  if(UseDEM || UseChrono){
    if(UseDEM){
      DemData=new StDemData[DemDataSize];
      memset(DemData,0,sizeof(StDemData)*DemDataSize);
    }
    for(unsigned c=0;c<parts.CountBlocks();c++){
      const JCasePartBlock &block=parts.GetBlock(c);
      if(IsBound(block.Type)){
        const word mkbound=block.GetMkType();
        const unsigned cmk=MkInfo->GetMkBlockByMkBound(mkbound);
        if(cmk>=MkInfo->Size())Run_Exceptioon(fun::PrintStr("Error loading boundary objects. Mkbound=%u is unknown.",mkbound));
        //-Loads data for collisions using DEM.
        if(UseDEM){
          const unsigned tav=CODE_GetTypeAndValue(MkInfo->Mkblock(cmk)->Code);
          DemData[tav]=LoadDemData(true,&block);
        }
        //-Loads data for collisions using Chrono.
        if(ChronoObjects && ChronoObjects->UseDataDVI(mkbound)){
          const StDemData data=LoadDemData(false,&block);
          if(block.Type==TpPartFloating)ChronoObjects->ConfigDataBodyFloating(mkbound,data.kfric,data.sfric,data.restitu,data.young,data.poisson);
          if(block.Type==TpPartMoving)  ChronoObjects->ConfigDataBodyMoving  (mkbound,data.kfric,data.sfric,data.restitu,data.young,data.poisson);
          if(block.Type==TpPartFixed)   ChronoObjects->ConfigDataBodyFixed   (mkbound,data.kfric,data.sfric,data.restitu,data.young,data.poisson);
        }
      }
    }
  }

  //-Configuration of Inlet/Outlet.
  if(xml.GetNodeSimple("case.execution.special.inout",true)){
    InOut=new JSphInOut(Cpu,CSP,FileXml,&xml,"case.execution.special.inout",DirCase);
    NpDynamic=true;
    ReuseIds=InOut->GetReuseIds();
    if(ReuseIds)Run_Exceptioon("Inlet/Outlet with ReuseIds is not a valid option for now...");
  }
  
  //-Configuration of boundary extrapolated correction.
  if(xml.GetNodeSimple("case.execution.special.boundextrap",false))Run_Exceptioon("The XML section 'boundextrap' is obsolete.");
  if(xml.GetNodeSimple("case.execution.special.boundcorr",true)){
    BoundCorr=new JSphBoundCorr(Cpu,Dp,&xml,"case.execution.special.boundcorr",MkInfo);
  }
 
  //-Configuration of Moorings object.
  if(xml.GetNodeSimple("case.execution.special.moorings",true)){
    if(WithFloating){
      if(!AVAILABLE_MOORDYN)Run_Exceptioon("Code for moorings and MoorDyn+ coupling is not included in the current compilation.");
      Moorings=new JDsMooredFloatings(DirCase,CaseName,Gravity,TimeMax,TimePart);
      Moorings->LoadXml(&xml,"case.execution.special.moorings");
    }
    else Log->PrintWarning("The use of Moorings was disabled because there are no floating objects...");
  }

  //-Configuration of ForcePoints object.
  if(Moorings || xml.GetNodeSimple("case.execution.special.forcepoints",true)){
    if(WithFloating){
      ForcePoints=new JDsFtForcePoints(Cpu,Dp,FtCount);
      //FtForces->LoadXml(&xml,"case.execution.special.forcepoints");
    }
    else Log->PrintWarning("The use of impossed force to floatings was disabled because there are no floating objects...");
  }

  //-Checks invalid options for symmetry. //<vs_syymmetry_ini>
  if(Symmetry){
    if(Simulate2D)  Run_Exceptioon("Symmetry is not allowed with 2-D simulations.");
    if(PeriY)       Run_Exceptioon("Symmetry is not allowed with periodic conditions in axis Y.");
    if(WithFloating)Run_Exceptioon("Symmetry is not allowed with floating bodies.");
    if(UseChrono)   Run_Exceptioon("Symmetry is not allowed with Chrono objects.");
    if(BoundCorr)   Run_Exceptioon("Symmetry is not allowed with BoundCor.");
    if(TVisco!=VISCO_Artificial)Run_Exceptioon("Symmetry is only allowed with Artificial viscosity.");
  } //<vs_syymmetry_end>

  NpMinimum=CaseNp-unsigned(PartsOutMax*CaseNfluid);
  Log->Print("**Basic case configuration is loaded");
}

//==============================================================================
/// Loads coefficients used for DEM or Chrono objects.
//==============================================================================
StDemData JSph::LoadDemData(bool checkdata,const JCasePartBlock* block)const{
  const word mk=block->GetMk();
  const word mkbound=block->GetMkType();
  StDemData data;
  memset(&data,0,sizeof(StDemData));
  //-Checks necessary values for collisions using DEM or Chrono.
  if(checkdata){
    const string objdesc=fun::PrintStr("Object mk=%u (mkbound=%u)",mk,mkbound);
    if(!block->ExistsSubValue("Kfric","value"))Run_Exceptioon(objdesc+" - Value of Kfric is invalid.");
    //if(!block->ExistsSubValue("Sfric","value"))Run_Exceptioon(objdesc+" - Value of Sfric is invalid.");
    if(!block->ExistsSubValue("Restitution_Coefficient","value"))Run_Exceptioon(objdesc+" - Value of Restitution_Coefficient is invalid.");
    if(!block->ExistsSubValue("Young_Modulus","value"))Run_Exceptioon(objdesc+" - Value of Young_Modulus is invalid.");
    if(!block->ExistsSubValue("PoissonRatio","value"))Run_Exceptioon(objdesc+" - Value of PoissonRatio is invalid.");
  }
  //-Loads necessary values for collisions using DEM or Chrono.
  data.kfric=block->GetSubValueFloat("Kfric","value",true,FLT_MAX);
  //data.sfric=block->GetSubValueFloat("Sfric","value",true,FLT_MAX); 
  data.restitu=block->GetSubValueFloat("Restitution_Coefficient","value",true,FLT_MAX);
  data.young=block->GetSubValueFloat("Young_Modulus","value",true,FLT_MAX);
  data.poisson=block->GetSubValueFloat("PoissonRatio","value",true,FLT_MAX);
  if(block->ExistsValue("Kfric_User"))data.kfric=block->GetValueFloat("Kfric_User");
  //if(block->ExistsValue("Sfric_User"))data.sfric=block->GetValueFloat("Sfric_User");
  if(block->ExistsValue("Restitution_Coefficient_User"))data.restitu=block->GetValueFloat("Restitution_Coefficient_User");

  //-Checks and initialises the friction coefficients when one of them was not defined.
  data.sfric=FLT_MAX; //-Disable sfric.
  if     (data.kfric==FLT_MAX && data.sfric!=FLT_MAX)data.kfric=data.sfric; 
  else if(data.kfric!=FLT_MAX && data.sfric==FLT_MAX)data.sfric=data.kfric;

  //-Loads necessary values for DEM.
  data.massp=MassBound;
  if(block->Type==TpPartFloating){
    const JCasePartBlock_Floating *fblock=(const JCasePartBlock_Floating *)block;
    data.mass=(float)fblock->GetMassbody();
    data.massp=(float)fblock->GetMasspart();
  }
  data.tau=(data.young? (1-data.poisson*data.poisson)/data.young: 0);
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
  //-Assigns code to each group of particles.
  for(unsigned p=0;p<np;p++)code[p]=MkInfo->GetCodeById(idp[p]);
}

//==============================================================================
/// Load normals for boundary particles (fixed and moving).
//==============================================================================
void JSph::LoadBoundNormals(unsigned np,unsigned npb,const unsigned *idp
  ,const typecode *code,tfloat3 *boundnormal)
{
  memset(boundnormal,0,sizeof(tfloat3)*np);
  string filenordata=JNormalsMarrone::GetNormalDataFile(DirCase+CaseName);
  if(fun::FileExists(filenordata)){
    //-Load or compute final normals.
    const tdouble3 *pnor=NULL;
    unsigned pnorsize=0;
    JPartNormalData nd;
    JNormalsMarrone nmarrone;
    nd.LoadFile(DirCase+CaseName);
    UseNormalsFt=nd.GetFtSupport();
    if(nd.GetCountNormals()){
      //-Compute Marrone normals starting from normal data in NBI4 file.
      nd.Reset();
      const bool savevtknor=false; //-Saves normals calculated starting from NBI4 file.
      nmarrone.RunCase(DirCase+CaseName,DirOut,savevtknor);
      pnor=nmarrone.GetPartNor();
      pnorsize=nmarrone.GetPartNorSize();
    }
    else{
      //-Loads final normals in NBI4 file.
      pnor=nd.GetPartNormals();
      pnorsize=nd.GetNbound();
    }
    //-Applies final normals. Loads normals from boundary particle to boundary limit.
    if(pnorsize){
      if(pnorsize<npb)Run_ExceptioonFile("The number of final normals does not match fixed and moving particles.",filenordata);
      for(unsigned p=0;p<npb;p++)boundnormal[p]=ToTFloat3(pnor[p]);  //-For fixed and moving particles.
      Log->Printf("NormalDataFile=\"%s\"",filenordata.c_str());
    }
    //for(unsigned p=0;p<npb;p++)boundnormal[p]=ToTFloat3(pnor[p]*2.);  //-Normal from boundary particle to ghost node (full distance).
    ////-Removes normal of unselected boundaries.
    //if(BoundCorr && !BoundCorr->GetMkBoundList().empty()){
    //  JRangeFilter rg(BoundCorr->GetMkBoundList());
    //  const unsigned nc=MkInfo->Size();
    //  for(unsigned c=0;c<nc;c++){
    //    const JSphMkBlock* bk=MkInfo->Mkblock(c);
    //    if(!rg.CheckValue(bk->MkType) && (bk->Type==TpPartFixed || bk->Type==TpPartMoving)){
    //      for(unsigned p=0;p<bk->Count;p++)boundnormal[p+bk->Begin]=TFloat3(0);
    //    }
    //  }
    //}
  }
  else Log->Print("**File with normal data not found.");
}

//==============================================================================
/// Config normals to point from boundary particle to ghost node (full distance).
//==============================================================================
void JSph::ConfigBoundNormals(unsigned np,unsigned npb,const tdouble3 *pos
  ,const unsigned *idp,tfloat3 *boundnormal)
{
  //-Checks current normals.
  bool ftnor=false;
  for(unsigned p=0;p<np;p++)if(idp[p]<CaseNbound && boundnormal[p]!=TFloat3(0)){
    if(idp[p]>=CaseNpb){
      if(!UseNormalsFt)boundnormal[p]=TFloat3(0);
    }
  }
  UseNormalsFt=ftnor;
  //-Saves normals from boundary particles to boundary limit.
  const string file1="CfgInit_Normals.vtk";
  Log->AddFileInfo(DirOut+file1,"Saves VTK file with initial normals (from boundary particles to boundary limit).");
  SaveVtkNormals(file1,-1,np,npb,pos,idp,boundnormal);
  //-Config normals.
  unsigned nerr=0,nerrft=0;
  for(unsigned p=0;p<np;p++)if(idp[p]<CaseNbound){
    if(boundnormal[p]==TFloat3(0)){
      if(idp[p]<CaseNpb)nerr++;
      else nerrft++;
    }
    boundnormal[p]=(boundnormal[p]*2.f);
  }
  //-Saves normals from boundary particles to ghost node.
  const string file2="CfgInit_NormalsGhost.vtk";
  Log->AddFileInfo(DirOut+file2,"Saves VTK file with initial normals (from boundary particles to ghost node).");
  SaveVtkNormals(file2,-1,np,npb,pos,idp,boundnormal);
  if(nerr>0)Log->PrintfWarning("There are %u of %u fixed or moving boundary particles without normal data.",nerr,npb);
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
    Run_Exceptioon(fun::PrintStr("Domain limits %s are not valid.",fun::Double3gRangeStr(dmin,dmax).c_str()));
  //-Periodic domain configuration.
  if(!PeriX){ MapRealPosMin.x=dmin.x; MapRealPosMax.x=dmax.x; }
  if(!PeriY){ MapRealPosMin.y=dmin.y; MapRealPosMax.y=dmax.y; }
  if(!PeriZ){ MapRealPosMin.z=dmin.z; MapRealPosMax.z=dmax.z; }
  //-Symmetry domain configuration. //<vs_syymmetry>
  if(Symmetry)MapRealPosMin.y=0;    //<vs_syymmetry>
}

//==============================================================================
/// Configures value of main constants and loads CSP structure.
//==============================================================================
void JSph::ConfigConstants1(bool simulate2d){
  //-Computation of constants.
  const float kernelk=fsph::GetKernelFactor(TKernel);
  const double h=KernelH;
  const double kh=h*kernelk;
  KernelSize=float(kh); 
  KernelSize2=KernelSize*KernelSize;
  //-Computes kernel constants.
  switch(TKernel){
    case KERNEL_Cubic:       KCubic =fsph::GetKernelCubic_Ctes     (simulate2d,h);  break;
    case KERNEL_Wendland:    KWend  =fsph::GetKernelWendland_Ctes  (simulate2d,h);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
  //-Wendland constants are always computed since this kernel is used in some parts where other kernels are not defined (e.g. Gauge code does not support Cubic Spline).
  if(TKernel!=KERNEL_Wendland)KWend=fsph::GetKernelWendland_Ctes  (simulate2d,h);
#ifdef DISABLE_KERNELS_EXTRA
  if(TKernel!=KERNEL_Wendland)Run_Exceptioon("Only Wendland kernel was compiled in this binary distribution.");
#endif
  //-Computes other constants.
  Cs0=sqrt(double(Gamma)*double(CteB)/double(RhopZero));
  Eta2=float((h*0.1)*(h*0.1));
  //-Loads main SPH constants and configurations in CSP.
  memset(&CSP,0,sizeof(StCteSph));
  CSP.simulate2d    =Simulate2D;
  CSP.simulate2dposy=Simulate2DPosY;
  CSP.tkernel       =TKernel;
  CSP.kcubic        =KCubic;
  CSP.kwend         =KWend;
  CSP.kernelh       =KernelH;
  CSP.cteb          =CteB;
  CSP.gamma         =Gamma;
  CSP.rhopzero      =RhopZero;
  CSP.dp            =Dp;
  CSP.massfluid     =MassFluid;
  CSP.massbound     =MassBound;
  CSP.gravity       =Gravity;
  CSP.kernelsize    =KernelSize;
  CSP.kernelsize2   =KernelSize2;
  CSP.cs0           =Cs0;
  CSP.eta2          =Eta2;
}

//==============================================================================
/// Configures other constants and loads more values in CSP structure.
//==============================================================================
void JSph::ConfigConstants2(){
  //-Constants for Laminar viscosity + SPS turbulence model.
  if(TVisco==VISCO_LaminarSPS){
    const double dp_sps=(Simulate2D? sqrt(Dp*Dp*2.)/2.: sqrt(Dp*Dp*3.)/3.);  
    SpsSmag=float(pow((0.12*dp_sps),2));
    SpsBlin=float((2./3.)*0.0066*dp_sps*dp_sps); 
  }
  //-Constants for DDT.
  DDTkh=KernelSize*DDTValue;
  DDTgz=float(double(RhopZero)*double(fabs(Gravity.z))/double(CteB));
  //-Constants for Dt.
  if(!DtIni)DtIni=KernelH/Cs0;
  if(!DtMin)DtMin=(KernelH/Cs0)*CoefDtMin;

  //-Loads main SPH constants and configurations in CSP.
  CSP.spssmag       =SpsSmag;
  CSP.spsblin       =SpsBlin;
  CSP.ddtkh         =DDTkh;
  CSP.ddtgz         =DDTgz;
}


//==============================================================================
/// Prints out configuration of the case.
//==============================================================================
void JSph::VisuConfig(){
  const string sep=" - ";
  Log->Print(Simulate2D? "**2D-Simulation parameters:": "**3D-Simulation parameters:");
  Log->Print(fun::VarStr("CaseName",CaseName));
  ConfigInfo=CaseName;
  Log->Print(fun::VarStr("RunName",RunName));
  //-Simulation dimension.
  if(Simulate2D){
    Log->Print(fun::VarStr("Simulate2DPosY",Simulate2DPosY));
    ConfigInfo=ConfigInfo+sep+"2D";
  }
  else ConfigInfo=ConfigInfo+sep+"3D";
  //-Symmetry. 
  Log->Print(fun::VarStr("Symmetry",Symmetry));  //<vs_syymmetry>
  if(Symmetry)ConfigInfo=ConfigInfo+sep+"Symmetry";
  //-SavePosDouble. 
  Log->Print(fun::VarStr("SavePosDouble",SvPosDouble));
  if(SvPosDouble)ConfigInfo=ConfigInfo+sep+"SvPosDouble";
  //-Other configurations. 
  Log->Print(fun::VarStr("SaveFtAce",SaveFtAce));
  if(FtMotSave)Log->Printf("SaveFtMotion=%s  (tout:%g)",(FtMotSave? "True": "False"),FtMotSave->GetTimeOut()); //<vs_ftmottionsv>
  Log->Print(fun::VarStr("SvTimers",SvTimers));
  if(DsPips)Log->Print(fun::VarStr("PIPS-steps",DsPips->StepsNum));
  //-Boundary. 
  Log->Print(fun::VarStr("Boundary",GetBoundName(TBoundary)));
  ConfigInfo=ConfigInfo+sep+GetBoundName(TBoundary);
  if(TBoundary==BC_MDBC){
    Log->Print(fun::VarStr("  SlipMode",GetSlipName(SlipMode)));
    Log->Print(fun::VarStr("  mDBC-Corrector",MdbcCorrector));
    Log->Print(fun::VarStr("  mDBC-FastSingle",MdbcFastSingle));
    Log->Print(fun::VarStr("  mDBC-Threshold",MdbcThreshold));
    ConfigInfo=ConfigInfo+"("+GetSlipName(SlipMode);
    if(MdbcCorrector)ConfigInfo=ConfigInfo+" - Corrector";
    if(MdbcFastSingle)ConfigInfo=ConfigInfo+" - FastSingle";
    if(MdbcThreshold>0)ConfigInfo=ConfigInfo+fun::PrintStr(" - Threshold=%g",MdbcThreshold);
    ConfigInfo=ConfigInfo+")";
  }
  //-StepAlgorithm. 
  Log->Print(fun::VarStr("StepAlgorithm",GetStepName(TStep)));
  ConfigInfo=ConfigInfo+sep+GetStepName(TStep);
  if(TStep==STEP_None)Run_Exceptioon("StepAlgorithm value is invalid.");
  if(TStep==STEP_Verlet){
    Log->Print(fun::VarStr("  VerletSteps",VerletSteps));
    ConfigInfo=ConfigInfo+sep+fun::PrintStr("(%d)",VerletSteps);
  }
  //-Kernel.
  std::vector<std::string> lines;
  fsph::GetKernelConfig(CSP,lines);
  Log->Print(lines);
  ConfigInfo=ConfigInfo+sep+fsph::GetKernelName(TKernel);

  //-Viscosity.
  Log->Print(fun::VarStr("Viscosity",GetViscoName(TVisco)));
  Log->Print(fun::VarStr("  Visco",Visco));
  Log->Print(fun::VarStr("  ViscoBoundFactor",ViscoBoundFactor));
  if(ViscoTime)Log->Print(fun::VarStr("ViscoTime",ViscoTime->GetFile()));
  ConfigInfo=ConfigInfo+sep+"Visco_"+GetViscoName(TVisco)+fun::PrintStr("(%gb%g)",Visco,ViscoBoundFactor);
  //-DensityDiffusion.
  Log->Print(fun::VarStr("DensityDiffusion",GetDDTName(TDensity)));
  ConfigInfo=ConfigInfo+sep+fun::PrintStr("DDT%d",int(TDensity));
  if(TDensity!=DDT_None){
    Log->Print(fun::VarStr("  DensityDiffusionValue",DDTValue));
    //Log->Print(fun::VarStr("DensityDiffusionArray",DDTArray));
    ConfigInfo=ConfigInfo+fun::PrintStr("(%g)",DDTValue);
  }
  if(TDensity==DDT_DDT2Full && KernelH/Dp>1.5)Log->PrintWarning("It is advised that selected DDT: \'Fourtakas et al 2019 (full)\' is used with several boundary layers of particles when h/dp>1.5 (2h <= layers*Dp)");
  //-Shifting.
  if(Shifting){
    Shifting->VisuConfig();
    ConfigInfo=ConfigInfo+sep+Shifting->GetConfigInfo();
  }
  else Log->Print(fun::VarStr("Shifting","None"));
  //-RigidAlgorithm.
  string rigidalgorithm=(!FtCount? "None": (FtMode==FTMODE_Sph? "SPH": "Collision-Free"));
  if(UseDEM)rigidalgorithm="SPH+DCDEM";
  if(UseChrono)rigidalgorithm="SPH+CHRONO";
  Log->Print(fun::VarStr("RigidAlgorithm",rigidalgorithm));
  if(FtCount){
    if(UseChrono)ConfigInfo=ConfigInfo+sep+"Ft-Chrono";
    else if(UseDEM)ConfigInfo=ConfigInfo+sep+"Ft-DCDEM";
    else ConfigInfo=ConfigInfo+sep+"Ft-SPH";
  }
  //-Moorings.
  if(Moorings)ConfigInfo=ConfigInfo+sep+"MoorDyn+";
  //-Other configurations. 
  Log->Print(fun::VarStr("FloatingCount",FtCount));
  if(FtCount)Log->Print(fun::VarStr("FtPause",FtPause));
  if(FtCount)Log->Print(fun::VarStr("FtConstraints",FtConstraints));
  if(FtCount)Log->Print(fun::VarStr("FtIgnoreRadius",FtIgnoreRadius));
  Log->Print(fun::VarStr("CaseNp",CaseNp));
  Log->Print(fun::VarStr("CaseNbound",CaseNbound));
  Log->Print(fun::VarStr("CaseNfixed",CaseNfixed));
  Log->Print(fun::VarStr("CaseNmoving",CaseNmoving));
  Log->Print(fun::VarStr("CaseNfloat",CaseNfloat));
  Log->Print(fun::VarStr("CaseNfluid",CaseNfluid));
  //-Periodic boundaries.
  Log->Print(fun::VarStr("PeriodicActive",TpPeriName(TpPeri(PeriActive))));
  if(PeriX)Log->Print(fun::VarStr("PeriodicXinc",PeriXinc));
  if(PeriY)Log->Print(fun::VarStr("PeriodicYinc",PeriYinc));
  if(PeriZ)Log->Print(fun::VarStr("PeriodicZinc",PeriZinc));
  if(PeriActive)ConfigInfo=ConfigInfo+sep+"Periodic_"+TpPeriName(TpPeri(PeriActive));
  //-Other configurations. 
  Log->Print(fun::VarStr("Dp",Dp));
  const float coefh=float(KernelH/(Dp*sqrt(Simulate2D? 2.f: 3.f)));
  Log->Printf("KernelH=%f  (CoefficientH=%g; H/Dp=%g)",KernelH,coefh,KernelH/Dp);
  const float kernelk=fsph::GetKernelFactor(TKernel);
  if(kernelk!=2.f)Log->Print(fun::VarStr("KernelK",kernelk));
  Log->Print(fun::VarStr("KernelSize",KernelSize));
  Log->Print(fun::VarStr("CteB",CteB));
  Log->Print(fun::VarStr("Gamma",Gamma));
  Log->Print(fun::VarStr("RhopZero",RhopZero));
  Log->Print(fun::VarStr("Cs0",Cs0));
  Log->Print(fun::VarStr("CFLnumber",CFLnumber));
  Log->Print(fun::VarStr("DtIni",DtIni));
  Log->Print(fun::VarStr("DtMin",DtMin));
  Log->Print(fun::VarStr("DtAllParticles",DtAllParticles));
  if(FixedDt){
    if(FixedDt->GetFixedValue())Log->Print(fun::VarStr("FixedDt",FixedDt->GetFixedValue()));
    else Log->Print(fun::VarStr("FixedDtFile",FixedDt->GetFile()));
  }
  Log->Print(fun::VarStr("MassFluid",MassFluid));
  Log->Print(fun::VarStr("MassBound",MassBound));
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
  //-RhopOut limits.
  Log->Print(fun::VarStr("RhopOut",RhopOut));
  if(RhopOut){
    Log->Print(fun::VarStr("RhopOutMin",RhopOutMin));
    Log->Print(fun::VarStr("RhopOutMax",RhopOutMax));
    ConfigInfo=ConfigInfo+sep+fun::PrintStr("RhopOut(%G-%G)",RhopOutMin,RhopOutMax);
  }
  Log->Print(fun::VarStr("WrnPartsOut",WrnPartsOut));
  //-Other configurations. 
  if(CteB==0)Run_Exceptioon("Constant \'b\' cannot be zero.\n\'b\' is zero when fluid height is zero (or fluid particles were not created)");
}


//==============================================================================
/// Prints out references of the case.
//==============================================================================
void JSph::VisuRefs(){
  Log->Print("[References]");
/////////////|---------1---------2---------3---------4---------5---------6---------7--------X8
  Log->Print("- Official solver reference DualSPHysics v5.0: J.M. Dominguez, G. Fourtakas,");
  Log->Print("    C. Altomare, R.B. Canelas, A. Tafuni, O. Garcia-Feal, I. Martinez-Estevez,"); 
  Log->Print("    A. Mokos, R. Vacondio, A.J.C. Crespo, B.D. Rogers, P.K. Stansby, M. Gomez-Gesteira."); 
  Log->Print("    2021. DualSPHysics: from fluid dynamics to multiphysics problems.");
  Log->Print("    Computational Particle Mechanics. doi: https://doi.org/10.1007/s40571-021-00404-2");
  Log->Print("");
  //-Code implementation:
  Log->Print("- Optimised CPU multi-core and GPU implementation (Dominguez et al., 2013  https://doi.org/10.1016/j.cpc.2012.10.015)");
  //-Boundary conditions:
  if(TBoundary==BC_DBC )Log->Print("- Dynamic boundary conditions (Crespo et al., 2007  https://doi.org/10.3970/cmc.2007.005.173)");
  if(TBoundary==BC_MDBC)Log->Print("- Modified Dynamic boundary conditions (English et al., 2021  https://doi.org/10.1007/s40571-021-00403-3)");
  //-Density diffusion Term:
  if(TDensity==DDT_DDT     )Log->Print("- Density diffusion Term: Molteni (Molteni and Colagrossi, 2009  https://doi.org/10.1016/j.cpc.2008.12.004)");
  if(TDensity==DDT_DDT2    )Log->Print("- Density diffusion Term: Fourtakas (Fourtakas et al., 2019  https://doi.org/10.1016/j.compfluid.2019.06.009)");
  if(TDensity==DDT_DDT2Full)Log->Print("- Density diffusion Term: Fourtakas (Fourtakas et al., 2019  https://doi.org/10.1016/j.compfluid.2019.06.009)");
  //-Viscosity:
  if(TVisco==VISCO_Artificial)Log->Print("- Viscosity: Artificial (Monaghan, 1992  https://doi.org/10.1146/annurev.aa.30.090192.002551)");
  if(TVisco==VISCO_LaminarSPS)Log->Print("- Viscosity: Laminar + SPS turbulence model (Dalrymple and Rogers, 2006  https://doi.org/10.1016/j.coastaleng.2005.10.004)");
  if(ViscoBoundFactor!=1     )Log->Print("- Viscosity: ViscoBoundFactor coefficient (Barreiro et al., 2014  https://doi.org/10.1371/journal.pone.0111031)");
  //-Kernel fuctions:
  if(TKernel==KERNEL_Cubic   )Log->Print("- Kernel: Cubic Spline (Monaghan, 1992  https://doi.org/10.1146/annurev.aa.30.090192.002551)");
  if(TKernel==KERNEL_Wendland)Log->Print("- Kernel: Quintic Wendland (Wendland, 1995  https://doi.org/10.1007/BF02123482)");
  //-Time integration scheme: 
  if(TStep==STEP_Verlet    )Log->Print("- Time integration scheme: Verlet (Verlet, 1967  https://doi.org/10.1103/PhysRev.159.98)");
  if(TStep==STEP_Symplectic)Log->Print("- Time integration scheme: Symplectic (Leimkhuler, 1996  https://doi.org/10.1007/978-3-319-16375-8_1)");
  //-Other features:
  if(PeriActive!=0)Log->Print("- Periodic open boundaries (Gomez-Gesteira et al., 2012  https://doi.org/10.1016/j.cageo.2012.02.029)");
  if(InOut        )Log->Print("- Inflow-outflow boundary conditions (Tafuni et al., 2018  https://doi.org/10.1016/j.cma.2018.08.004)");
  if(WithFloating )Log->Print("- Floating objects (Canelas et al., 2015  https://doi.org/10.1002/fld.4031)");
  if(UseDEM       )Log->Print("- Coupling SPH-DCDEM (Canelas et al., 2017  https://doi.org/10.1061/(ASCE)HY.1943-7900.0001331)");
  if(UseChrono    )Log->Print("- Coupling with Project Chrono (Canelas et al., 2018  https://doi.org/10.1016/j.apor.2018.04.015)");
  if(Moorings     )Log->Print("- Coupling with MoorDyn+ (Dominguez et al., 2019  https://doi.org/10.1016/j.coastaleng.2019.103560)");
  if(Shifting     )Log->Print("- Shifting algorithm (Lind et al., 2012  https://doi.org/10.1016/j.jcp.2011.10.027)");
  if(AccInput     )Log->Print("- External imposed forces (Longshaw and Rogers, 2015  https://doi.org/10.1016/j.advengsoft.2015.01.008)");
  //-Wave generation:
  const bool damp=(Damping!=NULL);
  const bool awas=(WaveGen && WaveGen->UseAwas());
  const bool lonw=(WaveGen && !WaveGen->WavesSolitary());
  const bool solw=(WaveGen && WaveGen->WavesSolitary());
  const bool inow=(InOut && InOut->Use_AwasVel());
  if(lonw        )Log->Print("- Long-crested wave generation (Altomare et al., 2017  https://doi.org/10.1016/j.coastaleng.2017.06.004)");
  if(solw        )Log->Print("- Solitary wave generation (Dominguez et al., 2019  https://doi.org/10.1080/21664250.2018.1560682)");
  if(MLPistons   )Log->Print("- Multi-layer Piston for wave generation (Altomare et al., 2015  https://doi.org/10.1142/S0578563415500242)");
  if(RelaxZones  )Log->Print("- Relaxation Zone for wave generation (Altomare et al., 2018  https://doi.org/10.1016/j.apor.2018.09.013)");
  if(inow        )Log->Print("- Wave generation and absorption using open boundaries (Verbrugghe et al., 2019  https://doi.org/10.1016/j.cpc.2019.02.003)");
  if(awas && damp)Log->Print("- AWAS & Damping: Active and passive wave absorption (Altomare et al., 2017  https://doi.org/10.1016/j.coastaleng.2017.06.004)");
  else if(awas)   Log->Print("- AWAS: Active Wave Absorption System (Altomare et al., 2017  https://doi.org/10.1016/j.coastaleng.2017.06.004)");
  else if(damp)   Log->Print("- Damping system (Altomare et al., 2017  https://doi.org/10.1016/j.coastaleng.2017.06.004)");
}

//==============================================================================
/// Shows particle and MK blocks summary.
//==============================================================================
void JSph::VisuParticleSummary()const{
  JXml xml; xml.LoadFile(FileXml);
  JCaseParts parts; 
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
        Run_Exceptioon("Found new particles out."); //-There can not be new particles excluded. | No puede haber nuevas particulas excluidas.
        dcell[p]=PC__CodeMapOut;
      }
    }
    else dcell[p]=PC__CodeMapOut;
  }
}

//==============================================================================
/// Initializes data of particles according XML configuration.
///
/// Inicializa datos de las particulas a partir de la configuracion en el XML.
//==============================================================================
void JSph::RunInitialize(unsigned np,unsigned npb,const tdouble3 *pos,const unsigned *idp
  ,const typecode *code,tfloat4 *velrhop,tfloat3 *boundnormal)
{
  if(!PartBegin){
    JDsInitialize init(Simulate2D,Simulate2DPosY,MapRealPosMin,MapRealPosMax,Dp,KernelH,DirCase,CaseNbound,boundnormal!=NULL);
    //-Loads configuration from XML.
    {
      JXml xml; xml.LoadFile(FileXml);
      xml.SetNuxLib(NuxLib); //-Enables the use of NuxLib in XML configuration.
      if(xml.GetNodeSimple("case.execution.special.initialize",true)){
        init.LoadXml(&xml,"case.execution.special.initialize");
      }
    }
    //-Loads configuration from execution parameters.
    init.LoadExecParms(CfgRun->InitParms);
    //-Executes initialize tasks.
    if(init.Count()){
      //-Creates array with mktype value.
      word *mktype=new word[np];
      for(unsigned p=0;p<np;p++){
        const unsigned cmk=MkInfo->GetMkBlockByCode(code[p]);
        mktype[p]=(cmk<MkInfo->Size()? word(MkInfo->Mkblock(cmk)->MkType): USHRT_MAX);
      }
      init.Run(np,npb,pos,idp,mktype,velrhop,boundnormal);
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
  PartsInit=new JDsPartsInit(Simulate2D,Simulate2DPosY,Dp,MkInfo,np,pos,code);
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
  if(CellMode!=CELLMODE_Full && CellMode!=CELLMODE_Half)Run_Exceptioon("The CellMode is invalid.");
  ScellDiv=(CellMode==CELLMODE_Full? 1: 2);
  Scell=KernelSize/ScellDiv;
  MovLimit=Scell*0.9f;
  Map_Cells=TUint3(unsigned(ceil(Map_Size.x/Scell)),unsigned(ceil(Map_Size.y/Scell)),unsigned(ceil(Map_Size.z/Scell)));
  //-Prints configuration.
  Log->Print(fun::VarStr("CellMode",string(GetNameCellMode(CellMode))));
  Log->Print(fun::VarStr("ScellDiv",ScellDiv));
  Log->Print(string("MapCells=(")+fun::Uint3Str(Map_Cells)+")");
  Log->Print(fun::VarStr("CellDomFixed",CellDomFixed));
  //-Creates VTK file with map cells.
  if(SaveMapCellsVtkSize()<1024*1024*10)SaveMapCellsVtk(Scell);
  else Log->PrintWarning("File CfgInit_MapCells.vtk was not created because number of cells is too high.");
}

//==============================================================================
/// Sets local domain of simulation within Map_Cells and computes DomCellCode.
/// Establece dominio local de simulacion dentro de Map_Cells y calcula DomCellCode.
//==============================================================================
void JSph::SelecDomain(tuint3 celini,tuint3 celfin){
  DomCelIni=celini;
  DomCelFin=celfin;
  DomCells=DomCelFin-DomCelIni;
  if(DomCelIni.x>=Map_Cells.x || DomCelIni.y>=Map_Cells.y || DomCelIni.z>=Map_Cells.z )Run_Exceptioon("DomCelIni is invalid.");
  if(DomCelFin.x>Map_Cells.x || DomCelFin.y>Map_Cells.y || DomCelFin.z>Map_Cells.z )Run_Exceptioon("DomCelFin is invalid.");
  if(DomCells.x<1 || DomCells.y<1 || DomCells.z<1 )Run_Exceptioon("The domain of cells is invalid.");
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
  if(!DomCellCode)Run_Exceptioon(string("Failed to select a valid CellCode for ")+fun::UintStr(DomCells.x)+"x"+fun::UintStr(DomCells.y)+"x"+fun::UintStr(DomCells.z)+" cells (CellMode="+GetNameCellMode(CellMode)+").");
  //-Prints configurantion.
  Log->Print(string("DomCells=(")+fun::Uint3Str(DomCells)+")");
  Log->Print(fun::VarStr("DomCellCode",fun::UintStr(PC__GetSx(DomCellCode))+"_"+fun::UintStr(PC__GetSy(DomCellCode))+"_"+fun::UintStr(PC__GetSz(DomCellCode))));
  //-Checks dimension for PosCell use.
  if(!Cpu)ConfigPosCellGpu();
}

//==============================================================================
/// Calculates the distribution of 32 bits between the X, Y and Z axes.
/// Calcula el reparto de 32 bits entre los ejes X, Y y Z.
//==============================================================================
tuint3 JSph::CalcCellDistribution(tuint3 ncells){
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
    sxmin=sx; symin=sy; szmin=sz;
  }
  return(TUint3(sxmin,symin,szmin));
}

//==============================================================================
/// Selects an adequate code for cell configuration.
/// Selecciona un codigo adecuado para la codificion de celda.
//==============================================================================
unsigned JSph::CalcCellCode(tuint3 ncells){
  const tuint3 scells=CalcCellDistribution(ncells);
  unsigned ccode=0;
  if(scells.x+scells.y+scells.z==32)ccode=PC__GetCode(scells.x,scells.y,scells.z);
  return(ccode);
}

//==============================================================================
/// Checks static configuration for PosCell on GPU.
/// Comprueba la configuracion estatica para PosCell en GPU.
//==============================================================================
void JSph::ConfigPosCellGpu(){
  //-Checks PosCellCode configuration is valid.
  const unsigned bz=CEL_MOVY;
  const unsigned by=CEL_MOVX-bz;
  const unsigned bx=32-by-bz;
  const unsigned nx=1<<bx,ny=1<<by,nz=1<<bz;
  Log->Print(fun::VarStr("PosCellCode",fun::PrintStr("%d_%d_%d (%d,%d,%d)",bx,by,bz,nx,ny,nz)));
  if(bx+by+bz!=32)Run_Exceptioon("Number of bits for cells is not 32.");
  const unsigned maskx=(UINT_MAX>>(by+bz))<<(by+bz);
  const unsigned maskz=(UINT_MAX<<(bx+by))>>(bx+by);
  const unsigned masky=(UINT_MAX&(~(maskx|maskz)));
  //Log->Printf("===> X(%02d):[%u]==[%u] cells:",bx,maskx,CEL_X);
  //Log->Printf("===> Y(%02d):[%u]==[%u]",by,masky,CEL_Y);
  //Log->Printf("===> Z(%02d):[%u]==[%u]",bz,maskz,CEL_Z);
  if(maskx!=CEL_X)Run_Exceptioon("Mask for cell-X is wrong.");
  if(masky!=CEL_Y)Run_Exceptioon("Mask for cell-Y is wrong.");
  if(maskz!=CEL_Z)Run_Exceptioon("Mask for cell-Z is wrong.");
  //-Config PosCellSize and check PosCellCode is enough for current simulation.
  PosCellSize=KernelSize;//-It is KernelSize by default. 
  int nks=1;
  tuint3 nposcells=TUint3(unsigned(ceil(Map_Size.x/PosCellSize)),unsigned(ceil(Map_Size.y/PosCellSize)),unsigned(ceil(Map_Size.z/PosCellSize)));
  const bool adapt_poscellsize=true;
  if(adapt_poscellsize){
    while(nposcells.x>nx || nposcells.y>ny || nposcells.z>nz){
      nks=(nks==1? 2: nks+2);
      PosCellSize=KernelSize*nks;
      nposcells=TUint3(unsigned(ceil(Map_Size.x/PosCellSize)),unsigned(ceil(Map_Size.y/PosCellSize)),unsigned(ceil(Map_Size.z/PosCellSize)));
    }
  }
  //-Shows PosCellSize value.
  Log->Printf("PosCellSize=%f  (%d x KernelSize)",PosCellSize,nks);
  //-Checks PosCellCode is enough for current simulation.
  if(nposcells.x>nx || nposcells.y>ny || nposcells.z>nz){
    Log->Printf("\n*** Attention ***");
    string tx1="The static cell configuration for PosCell on GPU approach is invalid for the current domain";
    string tx2="since the maximum number of cells for each axis is";
    Log->Printf("%s (%u x %u x %u cells), %s %u x %u x %u cells.",tx1.c_str()
      ,Map_Cells.x,Map_Cells.y,Map_Cells.z,tx2.c_str(),nx,ny,nz);
    const tuint3 scells=CalcCellDistribution(Map_Cells);
    if(scells.x+scells.y+scells.z>32)Run_Exceptioon("The number of cells is too large for a 32-bit PosCell configuration. The number of cells should be reduced.");
    Log->Printf("\nThe current configuration can be changed by the user by modifying the DualSphDef.h file and compiling the program again. ");
    Log->Printf("Replace the following code in DualSphDef.h:");
    Log->Printf("  //#define CEL_CONFIG_USER");
    Log->Printf("  #ifdef CEL_CONFIG_USER");
    Log->Printf("     ...");
    Log->Printf("  #else\n");
    Log->Printf("With:");
    Log->Printf("  #define CEL_CONFIG_USER");
    Log->Printf("  #ifdef CEL_CONFIG_USER");
    const unsigned cex0=(1<<(scells.y+scells.z));
    const unsigned cey0=(1<<(scells.z));
    const unsigned cez0=(1);
    unsigned cex=cex0,cey=cey0,cez=cez0;
    for(unsigned c=1;c<scells.x;c++)cex=(cex<<1)|cex0;
    for(unsigned c=1;c<scells.y;c++)cey=(cey<<1)|cey0;
    for(unsigned c=1;c<scells.z;c++)cez=(cez<<1)|cez0;
    Log->Printf("    #define CEL1_X 0x%08x  //-Mask of bits for cell X: %u bits for %u cells",cex,scells.x,1<<scells.x);
    Log->Printf("    #define CEL1_Y 0x%08x  //-Mask of bits for cell Y: %u bits for %u cells",cey,scells.y,1<<scells.y);
    Log->Printf("    #define CEL1_Z 0x%08x  //-Mask of bits for cell Z: %u bits for %u cells",cez,scells.z,1<<scells.z);
    Log->Printf("    #define CEL1_MOVX %7u  //-Displacement to obaint X cell.",scells.y+scells.z);
    Log->Printf("    #define CEL1_MOVY %7u  //-Displacement to obaint Y cell.",scells.z);
    Log->Printf("  #else\n");
    Run_Exceptioon("Current configuration for PosCell on GPU is invalid. More information above.");
  }
}

//==============================================================================
/// Computes maximum distance between particles and center of floating.
/// Calcula distancia maxima entre particulas y centro de cada floating.
//==============================================================================
void JSph::CalcFloatingRadius(unsigned np,const tdouble3 *pos,const unsigned *idp){
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
    if(ridp[fp]==UINT_MAX)Run_Exceptioon("There are floating particles not found.");
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
  const string errtex="The floating body radius (%g [m]) is too large for periodic distance in %c (%g [m]). If the floating body crosses the periodical limits, the simulation may be incorrect. *** If you want to avoid this initial verification please use \'FtIgnoreRadius\' in the XML file (parameters section).";
  if(PeriX && fabs(PeriXinc.x)<=radiusmax){
    const string tx=fun::PrintStr(errtex.c_str(),radiusmax,'X',abs(PeriXinc.x));
    if(FtIgnoreRadius)Log->PrintWarning(tx); else Run_Exceptioon(tx);
  }
  if(PeriY && fabs(PeriYinc.y)<=radiusmax){
    const string tx=fun::PrintStr(errtex.c_str(),radiusmax,'Y',abs(PeriYinc.y));
    if(FtIgnoreRadius)Log->PrintWarning(tx); else Run_Exceptioon(tx);
  }
  //Log->Printf("\nPeriYinc.y:%f <= %f:radiusmax",fabs(PeriYinc.y),radiusmax);
  if(PeriZ && fabs(PeriZinc.z)<=radiusmax){
    const string tx=fun::PrintStr(errtex.c_str(),radiusmax,'Z',abs(PeriZinc.z));
    if(FtIgnoreRadius)Log->PrintWarning(tx); else Run_Exceptioon(tx);
  }
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
  if(PartBegin){
    //-Loads particle blocks information.
    JPartDataHead parthead;
    parthead.LoadFile(PartBeginDir);
    const string filehead=fun::GetDirWithSlash(PartBeginDir)+parthead.GetFileName();
    //-Checks particle blocks information.
    const unsigned nmk=MkInfo->Size();
    if(nmk!=parthead.MkBlockCount())Run_ExceptioonFile("Number of Mk blocks is invalid.",filehead);
    for(unsigned c=0;c<nmk;c++){
      const JSphMkBlock* pmk=MkInfo->Mkblock(c);
      const JPartDataHeadMkBlock& mbk=parthead.Mkblock(c);
      if(pmk->Type!=mbk.Type)Run_ExceptioonFile(fun::PrintStr("Type of Mk block %u does not match.",c),filehead);
      if(pmk->Mk!=mbk.Mk)Run_ExceptioonFile(fun::PrintStr("Mk value of Mk block %u does not match.",c),filehead);
      if(pmk->MkType!=mbk.MkType)Run_ExceptioonFile(fun::PrintStr("MkType value of Mk block %u does not match.",c),filehead);
      if(pmk->Count!=mbk.Count)Run_ExceptioonFile(fun::PrintStr("Count value of Mk block %u does not match.",c),filehead);
    }
  }
}

//==============================================================================
/// Checks the initial density of fluid particles. If some particle is out of 
/// limits throws an exception.
/// Comprueba la densidad inicial de las particulas fluido. Si alguna particula 
/// esta fuera de los limites, lanza una excepcion.
//==============================================================================
void JSph::CheckRhopLimits(){
  const tfloat4 *velrhop=PartsLoaded->GetVelRhop(); ///<Velocity and density of each particle
  const unsigned *idp   =PartsLoaded->GetIdp();     ///<Identifier of each particle
  const unsigned n=PartsLoaded->GetCount();
  //-Checks the initial density of each fluid particle
  for(unsigned p=0;p<n;p++)if(idp[p]>=CaseNbound){
    if(velrhop[p].w<RhopOutMin || RhopOutMax<velrhop[p].w)
      Run_Exceptioon("Initial fluid density is out of limits. *** To change the limits modify the value of \'RhopOutMin\' and \'RhopOutMax\' in the XML file (parameters section).");
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
  PartsLoaded->CheckConfig(CaseNp,CaseNfixed,CaseNmoving,CaseNfloat,CaseNfluid,Simulate2D,Simulate2DPosY,TpPeri(PeriActive));

  if(PartBegin)RestartCheckData();
  Log->Printf("Loaded particles: %u",PartsLoaded->GetCount());

  //-Checks if the initial density of fluid particles is out of limits.
  //-Comprueba si la densidad inicial de las particulas fluido esta fuera de los limites.
  CheckRhopLimits();

  //-Collect information of loaded particles.
  //-Recupera informacion de las particulas cargadas.
  CasePosMin=PartsLoaded->GetCasePosMin();
  CasePosMax=PartsLoaded->GetCasePosMax();

  //-Computes actual limits of simulation.
  //-Calcula limites reales de la simulacion.
  if(PartsLoaded->MapSizeLoaded())PartsLoaded->GetMapSize(MapRealPosMin,MapRealPosMax);
  else{
    PartsLoaded->CalculeLimits(double(KernelH)*BORDER_MAP,Dp/2.,PeriX,PeriY,PeriZ,MapRealPosMin,MapRealPosMax);
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
  if(PeriX){ Map_PosMin.x=Map_PosMin.x-KernelSize;  Map_PosMax.x=Map_PosMax.x+KernelSize; }
  if(PeriY){ Map_PosMin.y=Map_PosMin.y-KernelSize;  Map_PosMax.y=Map_PosMax.y+KernelSize; }
  if(PeriZ){ Map_PosMin.z=Map_PosMin.z-KernelSize;  Map_PosMax.z=Map_PosMax.z+KernelSize; }
  Map_Size=Map_PosMax-Map_PosMin;
  //-Saves initial domain in a VTK file (CasePosMin/Max, MapRealPosMin/Max and Map_PosMin/Max).
  SaveInitialDomainVtk();
}

//==============================================================================
/// Initialisation of variables and objects for execution.
/// Inicializa variables y objetos para la ejecucion.
//==============================================================================
void JSph::InitRun(unsigned np,const unsigned *idp,const tdouble3 *pos){
  InterStep=(TStep==STEP_Symplectic? INTERSTEP_SymPredictor: INTERSTEP_Verlet);
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
  if(DsMotion){
    DsMotion->SetTimeMod(!PartIni? PartBeginTimeStep: 0);
    DsMotion->ProcesTime(JDsMotion::MOMT_Simple,0,TimeStepIni);
  }
  //-Adjust motion paddles for the instant of the loaded PART.
  if(WaveGen)WaveGen->SetTimeMod(!PartIni? PartBeginTimeStep: 0);
  //-Adjust Multi-layer pistons for the instant of the loaded PART.
  if(MLPistons)MLPistons->SetTimeMod(!PartIni? PartBeginTimeStep: 0);

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
  xml.SetNuxLib(NuxLib); //-Enables the use of NuxLib in XML configuration.

  //-Configuration of GaugeSystem.
  GaugeSystem->Config(CSP,Symmetry,TimeMax,TimePart,DomPosMin,DomPosMax,Scell,ScellDiv);
  if(xml.GetNodeSimple("case.execution.special.gauges",true))
    GaugeSystem->LoadXml(&xml,"case.execution.special.gauges",MkInfo);

  //-Prepares WaveGen configuration.
  if(WaveGen){
    Log->Print("Wave paddles configuration:");
    WaveGen->Init(GaugeSystem,MkInfo,TimeMax,TimePart);
    WaveGen->VisuConfig(""," ");
  }

  //-Prepares MLPistons configuration.
  if(MLPistons){
    Log->Printf("Multi-Layer Pistons configuration:");
    MLPistons->VisuConfig(""," ");
  }

  //-Prepares RelaxZones configuration.
  if(RelaxZones){
    Log->Print("Relaxation Zones configuration:");
    RelaxZones->Init(DirCase,TimeMax,Dp);
    RelaxZones->VisuConfig(""," ");
  }

  //-Prepares ChronoObjects configuration.
  if(ChronoObjects){
    Log->Print("Chrono Objects configuration:");
    if(PartBegin)Run_Exceptioon("Simulation restart not allowed when Chrono is used.");
    ChronoObjects->Init(MkInfo);
    ChronoObjects->VisuConfig(""," ");
  }

  //-Prepares Moorings configuration.
  if(Moorings){
    Log->Print("Moorings configuration:");
    if(!PartBegin)Moorings->Config(FtCount,FtObjs,ForcePoints);
    else Run_Exceptioon("Simulation restart not allowed when moorings are used.");
    Moorings->VisuConfig(""," ");
  }

  //-Prepares ForcePoints configuration.
  if(ForcePoints){
    Log->Printf("ForcePoints configuration:");
    ForcePoints->Config(FtCount,FtObjs,PeriActive,PeriX,PeriY,PeriZ,PeriXinc,PeriYinc,PeriZinc);
    if(!PartBegin)ForcePoints->CheckPoints(MkInfo,np,idp,pos);
    else Run_Exceptioon("Simulation restart not allowed when FtForces is used.");
    ForcePoints->VisuConfig(""," ",FtCount,FtObjs);
  }

  //-Prepares Damping configuration.
  if(Damping){
    Damping->VisuConfig("Damping configuration:"," ");
  }

  //-Prepares AccInput configuration.
  if(AccInput){
    Log->Print("AccInput configuration:");
    AccInput->VisuConfig(""," ");
    AccInput->Init(MkInfo);
  }

  //-Configuration of SaveDt.
  if(xml.GetNodeSimple("case.execution.special.savedt",true)){
    SaveDt=new JDsSaveDt();
    SaveDt->Config(&xml,"case.execution.special.savedt",TimeMax,TimePart);
    SaveDt->VisuConfig("SaveDt configuration:"," ");
  }

  //-Prepares BoundCorr configuration.
  if(BoundCorr){
    Log->Print("BoundCorr configuration:");
    if(PartBegin)Run_Exceptioon("Simulation restart not allowed when BoundCorr is used.");
    BoundCorr->RunAutoConfig(PartsInit);
    BoundCorr->VisuConfig(""," ");
  }

  //-Shows configuration of JGaugeSystem.
  if(GaugeSystem->GetCount())GaugeSystem->VisuConfig("GaugeSystem configuration:"," ");

  //-Shows configuration of JDsOutputTime.
  if(OutputTime->UseSpecialConfig()){
    vector<string> lines;
    OutputTime->GetConfig("TimeOut configuration:"," ",lines);
    Log->Print(lines);
  }

  Part=PartIni; Nstep=0; PartNstep=0; PartOut=0;
  TimeStep=TimeStepIni; TimeStepM1=TimeStep;
  if(FixedDt)DtIni=FixedDt->GetDt(TimeStep,DtIni);
  TimePartNext=(SvAllSteps? TimeStep: OutputTime->GetNextTime(TimeStep));
}

//==============================================================================
/// Calculates predefined movement of boundary particles.
/// Calcula movimiento predefinido de boundary particles.
//==============================================================================
bool JSph::CalcMotion(double stepdt){
  const bool motsim=true;
  const JDsMotion::TpMotionMode mode=(motsim? JDsMotion::MOMT_Simple: JDsMotion::MOMT_Ace2dt);
  DsMotion->ProcesTime(mode,TimeStep,stepdt);
  const bool active=DsMotion->GetActiveMotion();
  if(ChronoObjects && ChronoObjects->GetWithMotion() && active){
    const unsigned nref=DsMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=DsMotion->GetMotionData(ref);
      if(m.type!=MOTT_None)ChronoObjects->SetMovingData(m.mkbound,m.type==MOTT_Linear,m.linmov,m.matmov,stepdt);
    }
  }
  return(active);
}

//==============================================================================
/// Add motion from automatic wave generation.
/// Anhade movimiento de paddles calculado por generacion automatica de olas.
//==============================================================================
void JSph::CalcMotionWaveGen(double stepdt){
  const bool motsim=true;
  if(WaveGen){
    const bool svdata=(TimeStep+stepdt>=TimePartNext);
    for(unsigned c=0;c<WaveGen->GetCount();c++){
      const StMotionData m=(motsim? WaveGen->GetMotion(svdata,c,TimeStep,stepdt): WaveGen->GetMotionAce(svdata,c,TimeStep,stepdt));
      //Log->Printf("%u> t:%f  tp:%d  mx:%f  SetMotionData-WaveGen",Nstep,TimeStep,m.type,m.linmov.x);
      if(m.type!=MOTT_None){
        if(motsim)DsMotion->SetMotionData   (m);
        else      DsMotion->SetMotionDataAce(m);
      }
    }
  }
}

//==============================================================================
/// Applies the external velocities to each floating body of Chrono.
/// Aplica las velocidades externas a cada objeto flotante de Chrono.
//==============================================================================
void JSph::ChronoFtApplyImposedVel(){
  //-Applies imposed velocity.
  for(unsigned cf=0;cf<FtCount;cf++)if(FtObjs[cf].usechrono){
    const tfloat3 v1=(FtLinearVel [cf]!=NULL? FtLinearVel [cf]->GetValue3f(TimeStep): TFloat3(FLT_MAX));
    const tfloat3 v2=(FtAngularVel[cf]!=NULL? FtAngularVel[cf]->GetValue3f(TimeStep): TFloat3(FLT_MAX));
    ChronoObjects->SetFtDataVel(FtObjs[cf].mkbound,v1,v2);
  } 
}

//==============================================================================
/// Display a message with reserved memory for the basic data of particles.
/// Muestra un mensaje con la memoria reservada para los datos basicos de las particulas.
//==============================================================================
void JSph::PrintSizeNp(unsigned np,llong size,unsigned allocs)const{
  const double s=double(size)/(1024*1024);
  if(Cpu)Log->Printf("**Requested CPU memory for %u particles: %.1f MB.",np,s);
  else   Log->Printf("**Requested GPU memory for %u particles: %.1f MB (%u times).",np,s,allocs);
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
  //-Stores basic information of simulation data.
  JPartDataHead parthead;
  parthead.ConfigBasic(RunCode,AppName,CaseName,CasePosMin,CasePosMax
    ,Simulate2D,Simulate2DPosY,pieces,PartBeginFirst);
  MkInfo->ConfigPartDataHead(&parthead);
  parthead.ConfigCtes(Dp,KernelH,CteB,RhopZero,Gamma,MassBound,MassFluid,Gravity);
  parthead.ConfigSimNp(NpDynamic,ReuseIds);
  parthead.ConfigSimMap(MapRealPosMin,MapRealPosMax);
  parthead.ConfigSimPeri(TpPeriFromPeriActive(PeriActive),PeriXinc,PeriYinc,PeriZinc);
  parthead.ConfigSymmetry(Symmetry); //<vs_syymmetry>
  switch(TVisco){
    case VISCO_None:        parthead.ConfigVisco(JPartDataHead::VISCO_None      ,Visco,ViscoBoundFactor);  break;
    case VISCO_Artificial:  parthead.ConfigVisco(JPartDataHead::VISCO_Artificial,Visco,ViscoBoundFactor);  break;
    case VISCO_LaminarSPS:  parthead.ConfigVisco(JPartDataHead::VISCO_LaminarSPS,Visco,ViscoBoundFactor);  break;
    default: Run_Exceptioon("Viscosity type is unknown.");
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
    else Run_Exceptioon("The division configuration is invalid.");
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
    DataFloatBi4->Config(AppName,DirDataOut,MkInfo->GetMkBoundFirst(),FtCount,false);
    for(unsigned cf=0;cf<FtCount;cf++){
      const StFloatingData &ft=FtObjs[cf];
      DataFloatBi4->AddHeadData(cf,ft.mkbound,ft.begin,ft.count,ft.mass,ft.massp,ft.radius);
    }
    DataFloatBi4->SaveInitial();
    Log->AddFileInfo(DirDataOut+"PartFloat.fbi4","Binary file with floating body information for each instant (input for FloatingInfo program).");
  }
  //-Creates object to store excluded particles until recordering. 
  //-Crea objeto para almacenar las particulas excluidas hasta su grabacion.
  PartsOut=new JDsPartsOut();
}

//<vs_ftmottionsv_ini>  
//==============================================================================
/// Configures object to store floating motion data with high frequency.
//==============================================================================
void JSph::ConfigFtMotionSave(unsigned np,const tdouble3 *pos,const unsigned *idp){
  if(FtCount){
    FtMotSave->Config(AppName,DirDataOut,MkInfo->GetMkBoundFirst(),FtCount,FtObjs,np,pos,idp);
  }
  else{ 
    delete FtMotSave; FtMotSave=NULL; 
  }
}
//<vs_ftmottionsv_end>

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
void JSph::AbortBoundOut(JLog2 *log,unsigned nout,const unsigned *idp,const tdouble3 *pos
  ,const tfloat3 *vel,const float *rhop,const typecode *code)
{
  //-Prepares data of excluded boundary particles.
  byte* type=new byte[nout];
  byte* motive=new byte[nout];
  unsigned outfixed=0,outmoving=0,outfloat=0;
  unsigned outpos=0,outrhop=0,outmove=0;
  bool outxmin=false,outymin=false,outzmin=false;
  bool outxmax=false,outymax=false,outzmax=false;
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
    //-Checks out-position limits.
    if(CODE_GetSpecialValue(code[p])==CODE_OUTPOS){
      const tdouble3 rpos=pos[p];
      //-Check limits of real domain. | Comprueba limites del dominio reales.
      const double dx=rpos.x-MapRealPosMin.x;
      const double dy=rpos.y-MapRealPosMin.y;
      const double dz=rpos.z-MapRealPosMin.z;
      if(!PeriX && dx<0)outxmin=true;
      if(!PeriX && dx>=MapRealSize.x)outxmax=true;
      if(!PeriY && dy<0)outymin=true;
      if(!PeriY && dy>=MapRealSize.y)outymax=true;
      if(!PeriZ && dz<0)outzmin=true;
      if(!PeriZ && dz>=MapRealSize.z)outzmax=true;
    }
  }
  //-Shows excluded particles information.
  log->Print(" ");
  log->Print("*** ERROR: Some boundary particle was excluded. ***");
  if(UseChrono && PeriActive!=0)log->Print("*** Maybe some Chrono object went beyond the periodic limits. Be careful when combining the use of Chrono with periodic limits.");
  log->Printf("TimeStep: %f  (Nstep: %u)",TimeStep,Nstep);
  unsigned npunknown=nout-outfixed-outmoving-outfloat;
  if(!npunknown)log->Printf("Total boundary: %u  (fixed=%u  moving=%u  floating=%u)",nout,outfixed,outmoving,outfloat);
  else log->Printf("Total boundary: %u  (fixed=%u  moving=%u  floating=%u  UNKNOWN=%u)",nout,outfixed,outmoving,outfloat,npunknown);
  npunknown=nout-outpos-outrhop-outmove;
  if(!npunknown)log->Printf("Excluded for: position=%u  rhop=%u  velocity=%u",outpos,outrhop,outmove);
  else log->Printf("Excluded for: position=%u  rhop=%u  velocity=%u  UNKNOWN=%u",outpos,outrhop,outmove,npunknown);
  if(outxmin)log->Print("Some boundary particle exceeded the -X limit (left limit) of the simulation domain.");
  if(outxmax)log->Print("Some boundary particle exceeded the +X limit (right limit) of the simulation domain.");
  if(outymin)log->Print("Some boundary particle exceeded the -Y limit (front limit) of the simulation domain.");
  if(outymax)log->Print("Some boundary particle exceeded the +Y limit (back limit) of the simulation domain.");
  if(outzmin)log->Print("Some boundary particle exceeded the -Z limit (bottom limit) of the simulation domain.");
  if(outzmax)log->Print("Some boundary particle exceeded the +Z limit (top limit) of the simulation domain.");
  log->Print(" ");
  //-Creates VTK file.
  if(JVtkLib::Available()){
    JDataArrays arrays;
    arrays.AddArray("Pos"   ,nout,pos   ,false);
    arrays.AddArray("Idp"   ,nout,idp   ,false);
    arrays.AddArray("Vel"   ,nout,vel   ,false);
    arrays.AddArray("Rhop"  ,nout,rhop  ,false);
    arrays.AddArray("Type"  ,nout,type  ,false);
    arrays.AddArray("Motive",nout,motive,false);
    const string file=DirOut+"Error_BoundaryOut.vtk";
    log->AddFileInfo(file,"Saves the excluded boundary particles.");
    JVtkLib::SaveVtkData(file,arrays,"Pos");
  }
  //-Aborts execution.
  Run_Exceptioon("Fixed, moving or floating particles were excluded. Check VTK file Error_BoundaryOut.vtk with excluded particles.");
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
/// Adds basic data arrays in object JDataArrays.
//==============================================================================
void JSph::AddBasicArrays(JDataArrays &arrays,unsigned np,const tdouble3 *pos
  ,const unsigned *idp,const tfloat3 *vel,const float *rhop)const
{
  arrays.AddArray("Pos" ,np,pos);
  arrays.AddArray("Idp" ,np,idp);
  arrays.AddArray("Vel" ,np,vel);
  arrays.AddArray("Rhop",np,rhop);
}

//==============================================================================
/// Stores files of particle data.
/// Graba los ficheros de datos de particulas.
//==============================================================================
void JSph::SavePartData(unsigned npok,unsigned nout,const JDataArrays& arrays
  ,unsigned ndom,const tdouble3 *vdom,const StInfoPartPlus *infoplus)
{
  //-Stores particle data and/or information in bi4 format.
  //-Graba datos de particulas y/o informacion en formato bi4.
  if(DataBi4){
    tfloat3* posf3=NULL;
    TimerPart.Stop();
    tdouble3 domainmin=vdom[0];
    tdouble3 domainmax=vdom[1];
    for(unsigned c=1;c<ndom;c++){
      domainmin=MinValues(domainmin,vdom[c*2  ]);
      domainmax=MaxValues(domainmax,vdom[c*2+1]);
    }
    JBinaryData* bdpart=DataBi4->AddPartInfo(Part,TimeStep,npok,nout,Nstep,TimerPart.GetElapsedTimeD()/1000.,domainmin,domainmax,TotalNp);
    if(TStep==STEP_Symplectic)bdpart->SetvDouble("SymplecticDtPre",SymplecticDtPre);
    if(UseDEM)bdpart->SetvDouble("DemDtForce",DemDtForce); //(DEM)
    if(infoplus && SvData&SDAT_Info){
      bdpart->SetvDouble("dtmean",(!Nstep? 0: (TimeStep-TimeStepM1)/(Nstep-PartNstep)));
      bdpart->SetvDouble("dtmin",(!Nstep? 0: PartDtMin));
      bdpart->SetvDouble("dtmax",(!Nstep? 0: PartDtMax));
      if(FixedDt)bdpart->SetvDouble("dterror",FixedDt->GetDtError(true));
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
      if(ndom>1){
        bdpart->SetvUint("subdomain_count",ndom);
        for(unsigned c=0;c<ndom;c++){
          bdpart->SetvDouble3(fun::PrintStr("subdomainmin_%02u",c),vdom[c*2  ]);
          bdpart->SetvDouble3(fun::PrintStr("subdomainmax_%02u",c),vdom[c*2+1]);
        }
      }
    }
    if(SvData&SDAT_Binx){
      string err;
      if(!(err=arrays.CheckErrorArray("Pos" ,TypeDouble3,npok)).empty())Run_Exceptioon(err);
      if(!(err=arrays.CheckErrorArray("Idp" ,TypeUint   ,npok)).empty())Run_Exceptioon(err);
      if(!(err=arrays.CheckErrorArray("Vel" ,TypeFloat3 ,npok)).empty())Run_Exceptioon(err);
      if(!(err=arrays.CheckErrorArray("Rhop",TypeFloat  ,npok)).empty())Run_Exceptioon(err);
      const tdouble3 *pos =arrays.GetArrayDouble3("Pos");
      const unsigned *idp =arrays.GetArrayUint   ("Idp");
      const tfloat3  *vel =arrays.GetArrayFloat3 ("Vel");
      const float    *rhop=arrays.GetArrayFloat  ("Rhop");
      if(SvPosDouble){
        DataBi4->AddPartData(npok,idp,pos,vel,rhop);
      }
      else{
        posf3=GetPointerDataFloat3(npok,pos);
        DataBi4->AddPartData(npok,idp,posf3,vel,rhop);
      }
      //-Adds other arrays.
      const string arrignore=":Pos:Idp:Vel:Rhop:";
      for(unsigned ca=0;ca<arrays.Count();ca++){
        const JDataArrays::StDataArray arr=arrays.GetArrayData(ca);
        if(int(arrignore.find(string(":")+arr.keyname+":"))<0){//-Ignore main arrays.
          DataBi4->AddPartData(arr.keyname,npok,arr.ptr,arr.type);
        }
      }
      DataBi4->SaveFilePart();
    }
    if(SvData&SDAT_Info)DataBi4->SaveFileInfo();
    delete[] posf3;
  }

  //-Stores VTK nd/or CSV files.
  if((SvData&SDAT_Csv) || (SvData&SDAT_Vtk)){
    JDataArrays arrays2;
    arrays2.CopyFrom(arrays);

    string err;
    if(!(err=arrays2.CheckErrorArray("Pos" ,TypeDouble3,npok)).empty())Run_Exceptioon(err);
    if(!(err=arrays2.CheckErrorArray("Idp" ,TypeUint   ,npok)).empty())Run_Exceptioon(err);
    const tdouble3 *pos =arrays2.GetArrayDouble3("Pos");
    const unsigned *idp =arrays2.GetArrayUint   ("Idp");
    //-Generates array with posf3 and type of particle.
    tfloat3* posf3=GetPointerDataFloat3(npok,pos);
    byte *type=new byte[npok];
    for(unsigned p=0;p<npok;p++){
      const unsigned id=idp[p];
      type[p]=(id>=CaseNbound? 3: (id<CaseNfixed? 0: (id<CaseNpb? 1: 2)));
    }
    arrays2.DeleteArray("Pos");
    arrays2.AddArray("Pos",npok,posf3);
    arrays2.MoveArray(arrays2.Count()-1,0);
    arrays2.AddArray("Type",npok,type);
    arrays2.MoveArray(arrays2.Count()-1,4);
    //-Defines fields to be stored.
    if(SvData&SDAT_Vtk){
      JVtkLib::SaveVtkData(DirDataOut+fun::FileNameSec("PartVtk.vtk",Part),arrays2,"Pos");
    }
    if(SvData&SDAT_Csv){ 
      JOutputCsv ocsv(AppInfo.GetCsvSepComa());
      ocsv.SaveCsv(DirDataOut+fun::FileNameSec("PartCsv.csv",Part),arrays2);
    }
    //-Deallocate of memory.
    delete[] posf3;
    delete[] type; 
  }

  //-Stores data of excluded particles.
  if(DataOutBi4 && PartsOut->GetCount()){
    DataOutBi4->SavePartOut(SvPosDouble,Part,TimeStep,PartsOut->GetCount(),PartsOut->GetIdpOut(),NULL,PartsOut->GetPosOut(),PartsOut->GetVelOut(),PartsOut->GetRhopOut(),PartsOut->GetMotiveOut());
  }

  //-Stores data of floating bodies.
  if(DataFloatBi4){
    for(unsigned cf=0;cf<FtCount;cf++){
      const StFloatingData *v=FtObjs+cf;
      DataFloatBi4->AddPartData(cf,v->center,v->fvel,v->fomega,v->facelin,v->faceang);
    }
    DataFloatBi4->SavePartFloat(Part,Nstep,TimeStep,(UseDEM? DemDtForce: 0));
  }

  //-Empties stock of excluded particles.
  //-Vacia almacen de particulas excluidas.
  PartsOut->Clear();
}

//==============================================================================
/// Generates data output files.
/// Genera los ficheros de salida de datos.
//==============================================================================
void JSph::SaveData(unsigned npok,const JDataArrays& arrays
  ,unsigned ndom,const tdouble3 *vdom,const StInfoPartPlus *infoplus)
{
  string suffixpartx=fun::PrintStr("_%04d",Part);

  //-Counts new excluded particles.
  //-Contabiliza nuevas particulas excluidas.
  const unsigned noutpos=PartsOut->GetOutPosCount(),noutrhop=PartsOut->GetOutRhopCount(),noutmove=PartsOut->GetOutMoveCount();
  const unsigned nout=noutpos+noutrhop+noutmove;
  if(nout!=PartsOut->GetCount())Run_Exceptioon("Excluded particles with unknown reason.");
  AddOutCount(noutpos,noutrhop,noutmove);

  //-Stores data files of particles.
  SavePartData(npok,nout,arrays,ndom,vdom,infoplus);

  //-Reinitialises limits of dt. | Reinicia limites de dt.
  PartDtMin=DBL_MAX; PartDtMax=-DBL_MAX;

  //-Computation of time.
  if(Part>PartIni || Nstep){
    TimerPart.Stop();
    double tpart=TimerPart.GetElapsedTimeD()/1000;
    double tseg=tpart/(TimeStep-TimeStepM1);
    TimerSim.Stop();
    double tcalc=TimerSim.GetElapsedTimeD()/1000;
    double tleft=(tcalc/(TimeStep-TimeStepIni))*(TimeMax-TimeStep);
    Log->Printf("Part%s  %12.6f  %12d  %7d  %9.2f  %14s",suffixpartx.c_str(),TimeStep,(Nstep+1),Nstep-PartNstep,tseg,fun::GetDateTimeAfter(int(tleft)).c_str());
  }
  else Log->Printf("Part%s        %u particles successfully stored",suffixpartx.c_str(),npok);   
  
  //-Shows info of the new inlet particles.
  bool printnp=true;
  if(InOut && InOut->GetNewNpPart()){
    Log->Printf("  Particles new: %u (total new: %llu)  -  Current np: %u",InOut->GetNewNpPart(),InOut->GetNewNpTotal(),npok);
    InOut->ClearNewNpPart();
    printnp=false;
  }
  
  //-Shows info of the excluded particles.
  if(nout){
    PartOut+=nout;
    if(printnp)Log->Printf("  Particles out: %u  (total out: %u)  -  Current np: %u",nout,PartOut,npok);
    else       Log->Printf("  Particles out: %u  (total out: %u)",nout,PartOut);
  }

  //-Cheks number of excluded particles.
  if(WrnPartsOut && nout){
    //-Cheks number of excluded particles in one PART.
    if(PartsOutWrn<=100 && nout>=float(max(CaseNfluid,infoplus->npf))*(float(PartsOutWrn)/100.f)){
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
  if(ChronoObjects)ChronoObjects->SavePart(Part);
  if(Moorings)Moorings->SaveData(Part);
  if(ForcePoints)ForcePoints->SaveData(Part);
  if(BoundCorr && BoundCorr->GetUseMotion())BoundCorr->SaveData(Part);
  if(DsPips && DsPips->SvData && SvData!=SDAT_None)DsPips->SaveData();

  //-Checks request for simulation termination.
  CheckTermination();
}

//==============================================================================
/// Checks request for simulation termination.
/// Comprueba solicitud de terminar la simulacion.
//==============================================================================
void JSph::CheckTermination(){
  const string file=DirOut+"TERMINATE";
  const ullong tmodif=fun::FileModifTime(file);
  if(tmodif && tmodif!=TerminateMt){
    double tmax=0;
    ifstream pf;
    pf.open(file.c_str(),ios::binary);
    if(pf){
      pf.seekg(0,ios::end);
      const unsigned fsize=min(127u,unsigned(pf.tellg()));
      if(fsize){
        char buff[128];
        pf.seekg(0,ios::beg);
        pf.read(buff,fsize);
        buff[fsize]='\0';
        tmax=atof(buff);
      }
      pf.close();
    }
    if(tmax<TimeStep)tmax=TimeStep;
    Log->PrintfWarning("TERMINATE file has updated TimeMax from %gs to %gs (current time: %fs).",TimeMax,tmax,TimeStep);
    TimeMax=tmax;
  }
  TerminateMt=tmodif;
}

//==============================================================================
/// Generates VTK file with domain of the particles.
/// Genera fichero VTK con el dominio de las particulas.
//==============================================================================
void JSph::SaveDomainVtk(unsigned ndom,const tdouble3 *vdom)const{ 
  if(vdom){
    string fname=fun::FileNameSec("Domain.vtk",Part);
    JVtkLib::SaveVtkBoxes(DirDataOut+fname,ndom,vdom,KernelH*0.5f);
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
  JVtkLib::SaveVtkBoxes(file,nbox,vdomf3,0);
  delete[] vdomf3;
}

//==============================================================================
/// Returns size of VTK file with map cells.
/// Devuelve tamanho de fichero VTK con las celdas del mapa.
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
  tdouble3 pmin=Map_PosMin; //MapRealPosMin;
  tdouble3 pmax=pmin+TDouble3(scell*cells.x,scell*cells.y,scell*cells.z);
  if(Simulate2D)pmin.y=pmax.y=Simulate2DPosY;
  //-Creates lines.
  JVtkLib sh;
  //-Back lines.
  tdouble3 p0=TDouble3(pmin.x,pmax.y,pmin.z),p1=TDouble3(pmin.x,pmax.y,pmax.z);
  for(unsigned cx=0;cx<=cells.x;cx++)sh.AddShapeLine(p0+TDouble3(scell*cx,0,0),p1+TDouble3(scell*cx,0,0),0);
  p1=TDouble3(pmax.x,pmax.y,pmin.z);
  for(unsigned cz=0;cz<=cells.z;cz++)sh.AddShapeLine(p0+TDouble3(0,0,scell*cz),p1+TDouble3(0,0,scell*cz),0);
  if(!Simulate2D){
    //-Bottom lines.
    p0=TDouble3(pmin.x,pmin.y,pmin.z),p1=TDouble3(pmax.x,pmin.y,pmin.z);
    for(unsigned cy=0;cy<=cells.y;cy++)sh.AddShapeLine(p0+TDouble3(0,scell*cy,0),p1+TDouble3(0,scell*cy,0),1);
    p1=TDouble3(pmin.x,pmax.y,pmin.z);
    for(unsigned cx=0;cx<=cells.x;cx++)sh.AddShapeLine(p0+TDouble3(scell*cx,0,0),p1+TDouble3(scell*cx,0,0),1);
    //-Left lines.
    p0=TDouble3(pmin.x,pmin.y,pmin.z),p1=TDouble3(pmin.x,pmax.y,pmin.z);
    for(unsigned cz=0;cz<=cells.z;cz++)sh.AddShapeLine(p0+TDouble3(0,0,scell*cz),p1+TDouble3(0,0,scell*cz),2);
    p1=TDouble3(pmin.x,pmin.y,pmax.z);
    for(unsigned cy=0;cy<=cells.y;cy++)sh.AddShapeLine(p0+TDouble3(0,scell*cy,0),p1+TDouble3(0,scell*cy,0),2);
  }
  const string file=DirOut+"CfgInit_MapCells.vtk";
  Log->AddFileInfo(file,"Saves the cell division of the simulation domain.");
  sh.SaveShapeVtk(file,"axis");
}

//==============================================================================
/// Saves VTK file with normals of particles (degug).
/// Only normal (non-periodic) particles are allowed.
/// Graba fichero VTK con normales de las particulas (degug).
/// Solo se permiten particulas normales (no periodicas).
//==============================================================================
void JSph::SaveVtkNormals(std::string filename,int numfile,unsigned np,unsigned npb
  ,const tdouble3 *pos,const unsigned *idp,const tfloat3 *boundnormal)const
{
  if(JVtkLib::Available()){
    if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
    filename=DirOut+filename;
    //-Find floating particles.
    unsigned nfloat=0;
    unsigned* ftidx=NULL;
    //-Allocate memory for boundary particles.
    const unsigned npsel=npb+nfloat;
    JDataArrays arrays;
    tdouble3* vpos =arrays.CreateArrayPtrDouble3("Pos",npsel);
    unsigned* vidp =arrays.CreateArrayPtrUint   ("Idp",npsel);
    word*     vmk  =arrays.CreateArrayPtrWord   ("Mk" ,npsel);
    tfloat3*  vnor =arrays.CreateArrayPtrFloat3 ("Normal",npsel);
    float*    vsnor=arrays.CreateArrayPtrFloat  ("NormalSize",npsel);
    //-Loads data of fixed and moving particles.
    memcpy(vpos,pos,sizeof(tdouble3)*npb);
    memcpy(vidp,idp,sizeof(unsigned)*npb);
    MkInfo->GetMkByIds(npb,idp,vmk);
    memcpy(vnor,boundnormal,sizeof(tfloat3)*npb);
    //-Computes normalsize.
    for(unsigned p=0;p<npsel;p++)vsnor[p]=fgeo::PointDist(vnor[p]);
    //-Saves VTK file.
    JVtkLib::SaveVtkData(filename,arrays,"Pos");
    //-Frees memory.
    arrays.Reset();
  }
}

//==============================================================================
/// Adds basic information of resume to hinfo & dinfo.
/// Anhade la informacion basica de resumen a hinfo y dinfo.
//==============================================================================
void JSph::GetResInfo(float tsim,float ttot,std::string headplus,std::string detplus
  ,std::string &hinfo,std::string &dinfo)const
{
  hinfo=hinfo+"#RunName;Rcode-VersionInfo;DateTime;Np;TSimul;TSeg;TTotal;MemCpu;MemGpu";
  dinfo=dinfo+ RunName+ ";"+ RunCode+ "-"+ AppName+ ";"+ RunTimeDate+ ";"+ fun::UintStr(CaseNp);
  dinfo=dinfo+ ";"+ fun::FloatStr(tsim)+ ";"+ fun::FloatStr(tsim/float(TimeStep))+ ";"+ fun::FloatStr(ttot);
  dinfo=dinfo+ ";"+ fun::LongStr(MaxNumbers.memcpu)+ ";"+ fun::LongStr(MaxNumbers.memgpu);
  hinfo=hinfo+";Steps;GPIPS;PhysicalTime;PartFiles;PartsOut;MaxParticles;MaxCells";
  const unsigned nout=GetOutPosCount()+GetOutRhopCount()+GetOutMoveCount();
  const string gpips=(DsPips? fun::DoublexStr(DsPips->GetGPIPS(tsim),"%.10f"): "");
  dinfo=dinfo+ ";"+ fun::IntStr(Nstep)+ ";"+ gpips+ ";"+ fun::DoublexStr(TimeStep) + ";"+ fun::IntStr(Part)+ ";"+ fun::UintStr(nout);
  dinfo=dinfo+ ";"+ fun::UintStr(MaxNumbers.particles)+ ";"+ fun::UintStr(MaxNumbers.cells);
  hinfo=hinfo+";Hardware;RunMode;Configuration";
  dinfo=dinfo+ ";"+ Hardware+ ";"+ "Cells"+GetNameCellMode(CellMode) + " - " + RunMode+ ";"+ ConfigInfo;
  hinfo=hinfo+";Nbound;Nfixed;Dp;H";
  dinfo=dinfo+ ";"+ fun::UintStr(CaseNbound)+ ";"+ fun::UintStr(CaseNfixed);
  dinfo=dinfo+ ";"+ fun::FloatStr(float(Dp))+ ";"+ fun::FloatStr(KernelH);
  hinfo=hinfo+";PartsOutRhop;PartsOutVel";
  dinfo=dinfo+ ";"+ fun::UintStr(GetOutRhopCount())+ ";"+ fun::UintStr(GetOutMoveCount());
  hinfo=hinfo+ headplus;
  dinfo=dinfo+ detplus;
}

//==============================================================================
/// Generates file Run.csv with resume of execution.
/// Genera fichero Run.csv con resumen de ejecucion.
//==============================================================================
void JSph::SaveRes(float tsim,float ttot,const std::string &headplus,const std::string &detplus){
  const string fname=DirOut+"Run.csv";
  Log->AddFileInfo(fname,"One line CSV file with execution parameters and other simulation data.");
  ofstream pf;
  pf.open(fname.c_str());
  if(pf){
    string hinfo,dinfo;
    GetResInfo(tsim,ttot,headplus,detplus,hinfo,dinfo);
    pf << fun::StrCsvSep(CsvSepComa,hinfo) << endl << fun::StrCsvSep(CsvSepComa,dinfo) << endl;
    if(pf.fail())Run_ExceptioonFile("Failed writing to file.",fname);
    pf.close();
  }
  else Run_ExceptioonFile("File could not be opened.",fname);
}

//==============================================================================
/// Shows resume of execution.
/// Muestra resumen de ejecucion.
//==============================================================================
void JSph::ShowResume(bool stop,float tsim,float ttot,bool all,std::string infoplus){
  if(DsPips && DsPips->SvData)DsPips->SaveData();
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
    Log->Printf("Runtime per physical second......: %f sec.",tseg);
    //Log->Printf("Time per second of simulation....: %f sec.",tseg);
    Log->Printf("Steps per second.................: %f",nstepseg);
    Log->Printf("Steps of simulation..............: %d",Nstep);
    if(DsPips){
      Log->Printf("Particle Interactions Per Second.: %.8f GPIPS",DsPips->GetGPIPS(tsim));
      Log->Printf("Total particle interactions (f+b): %s",DsPips->GetTotalPIsInfo().c_str());
    }
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
  Log->Printf("Maximum number of particles......: %u",MaxNumbers.particles);
  Log->Printf("Maximum number of cells..........: %u",MaxNumbers.cells);
  Log->Printf("CPU Memory.......................: %lld (%.2f MB)",MaxNumbers.memcpu,double(MaxNumbers.memcpu)/(1024*1024));
  if(MaxNumbers.memgpu)Log->Printf("GPU Memory.......................: %lld (%.2f MB)",MaxNumbers.memgpu,double(MaxNumbers.memgpu)/(1024*1024));
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
/// Returns name of boundary method in text format.
/// Devuelve nombre de condicioens de contorno en texto.
//==============================================================================
std::string JSph::GetBoundName(TpBoundary tboundary){
  string tx;
  if(tboundary==BC_DBC)tx="DBC";
  else if(tboundary==BC_MDBC)tx="mDBC";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns name of Slip mode for mDBC in text format.
/// Devuelve nombre de modo Slip para mDBC en texto.
//==============================================================================
std::string JSph::GetSlipName(TpSlipMode tslip){
  string tx;
  if(tslip==SLIP_Vel0)tx="DBC vel=0";
  else if(tslip==SLIP_NoSlip)tx="No-slip";
  else if(tslip==SLIP_FreeSlip)tx="Free slip";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns name of Density Diffusion Term in text format.
/// Devuelve nombre del Density Diffusion Term en texto.
//==============================================================================
std::string JSph::GetDDTName(TpDensity tdensity)const{
  string tx;
  if(tdensity==DDT_None)tx="None";
  else if(tdensity==DDT_DDT)tx="Molteni and Colagrossi 2009";
  else if(tdensity==DDT_DDT2)tx="Fourtakas et al 2019 (inner)";
  else if(tdensity==DDT_DDT2Full)tx="Fourtakas et al 2019 (full)";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns Density Diffusion Term configuration in text format.
/// Devuelve configuracion de Density Diffusion Term en texto.
//==============================================================================
std::string JSph::GetDDTConfig()const{
  string tx="None";
  if(TDensity!=DDT_None)tx=GetDDTName(TDensity)+fun::PrintStr("(%g)",DDTValue);
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
    {//-For InOut particles.
      typecode tv=CODE_GetTypeAndValue(code[p+pini]);
      if(tv>=CODE_TYPE_FLUID_INOUT)xkind[p]=byte(tv-CODE_TYPE_FLUID_INOUT+10);
    }
  }
  //-Generates VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos" ,np,xpos,true);
  if(idp)  arrays.AddArray("Idp" ,np,idp+pini,false);
  if(xtype)arrays.AddArray("Type",np,xtype,true);
  if(xkind)arrays.AddArray("Kind",np,xkind,true);
  if(xvel) arrays.AddArray("Vel" ,np,xvel ,true);
  if(xrhop)arrays.AddArray("Rhop",np,xrhop,true);
  if(xace) arrays.AddArray("Ace" ,np,xace ,true);
  JVtkLib::SaveVtkData(filename,arrays,"Pos");
  arrays.Reset();
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
  JDataArrays arrays;
  arrays.AddArray("Pos" ,n,pos+pini,false);
  if(idp)  arrays.AddArray("Idp"  ,n,idp+pini,false);
  if(vel)  arrays.AddArray("Vel"  ,n,vel+pini,false);
  if(rhop) arrays.AddArray("Rhop" ,n,rhop+pini,false);
  if(check)arrays.AddArray("Check",n,check+pini,false);
  if(num)  arrays.AddArray("Num"  ,n,num,true);
  JVtkLib::SaveVtkData(filename,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Saves CSV file with particle data (degug).
/// Graba fichero CSV con datos de las particulas (degug).
//==============================================================================
void JSph::DgSaveCsvParticlesCpu(std::string filename,int numfile,unsigned pini,unsigned pfin
  ,std::string head,const tfloat3 *pos,const unsigned *idp,const tfloat3 *vel
  ,const float *rhop,const float *ar,const tfloat3 *ace,const tfloat3 *vcorr)
{
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
    if(pf.fail())Run_ExceptioonFile("Failed writing to file.",filename);
    pf.close();
  }
  else Run_ExceptioonFile("File could not be opened.",filename);
}


