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
#include "FunctionsGeo3d.h"
#include "JPartDataHead.h"
#include "JSphMk.h"
#include "JSphPartsInit.h"
#include "JPartsLoad4.h"
#include "JSphMotion.h"
#include "JXml.h"
#include "JSpaceCtes.h"
#include "JSpaceEParms.h"
#include "JSpaceParts.h"
#include "JSphDtFixed.h"
#include "JSaveDt.h"
#include "JTimeOut.h"
#include "JSphVisco.h"
#include "JGaugeSystem.h"
#include "JWaveGen.h"
#include "JMLPistons.h"     //<vs_mlapiston>
#include "JRelaxZones.h"    //<vs_rzone>
#include "JChronoObjects.h" //<vs_chroono>
#include "JMooredFloatings.h"  //<vs_moordyyn>
#include "JSphFtForcePoints.h" //<vs_moordyyn>
#include "JSphAccInput.h"
#include "JPartDataBi4.h"
#include "JPartOutBi4Save.h"
#include "JPartFloatBi4.h"
#include "JPartsOut.h"
#include "JShifting.h"
#include "JDamping.h"
#include "JSphInitialize.h"
#include "JSphInOut.h"       //<vs_innlet> 
#include "JSphBoundCorr.h"   //<vs_innlet> 
#include "JLinearValue.h"
#include "JPartNormalData.h" //<vs_mddbc>
#include "JNormalsMarrone.h" //<vs_mddbc>
#include "JDataArrays.h"
#include "JOutputCsv.h"
#include "JVtkLib.h"
#include "JNumexLib.h"
#include "JSpaceUserVars.h"

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
  ViscoTime=NULL;
  DtFixed=NULL;
  SaveDt=NULL;
  TimeOut=NULL;
  MkInfo=NULL;
  PartsInit=NULL;
  SphMotion=NULL;
  FtObjs=NULL;
  FtLinearVel=NULL;  //<vs_fttvel>
  FtAngularVel=NULL; //<vs_fttvel>
  DemData=NULL;
  GaugeSystem=NULL;
  WaveGen=NULL;
  MLPistons=NULL;     //<vs_mlapiston>
  RelaxZones=NULL;    //<vs_rzone>
  ChronoObjects=NULL; //<vs_chroono>
  Moorings=NULL;    //<vs_moordyyn>
  ForcePoints=NULL; //<vs_moordyyn>
  Shifting=NULL;
  Damping=NULL;
  AccInput=NULL;
  PartsLoaded=NULL;
  InOut=NULL;       //<vs_innlet>
  BoundCorr=NULL;   //<vs_innlet>
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
  delete DtFixed;       DtFixed=NULL;
  delete SaveDt;        SaveDt=NULL;
  delete TimeOut;       TimeOut=NULL;
  delete MkInfo;        MkInfo=NULL;
  delete PartsInit;     PartsInit=NULL;
  delete SphMotion;     SphMotion=NULL;
  AllocMemoryFloating(0,false);
  delete[] DemData;     DemData=NULL;
  delete GaugeSystem;   GaugeSystem=NULL;
  delete WaveGen;       WaveGen=NULL;
  delete MLPistons;     MLPistons=NULL;     //<vs_mlapiston>
  delete RelaxZones;    RelaxZones=NULL;    //<vs_rzone>
  delete ChronoObjects; ChronoObjects=NULL; //<vs_chroono>
  delete Moorings;      Moorings=NULL;      //<vs_moordyyn>
  delete ForcePoints;   ForcePoints=NULL;   //<vs_moordyyn>
  delete Shifting;      Shifting=NULL;
  delete Damping;       Damping=NULL;
  delete AccInput;      AccInput=NULL; 
  delete PartsLoaded;   PartsLoaded=NULL;
  delete InOut;         InOut=NULL;       //<vs_innlet>
  delete BoundCorr;     BoundCorr=NULL;   //<vs_innlet>
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
  Awen=Bwen=Agau=Bgau=0;
  Awc6=Bwc6=0;  //<vs_praticalsskq>
  memset(&CubicCte,0,sizeof(StCubicCte));
  TVisco=VISCO_None;
  TDensity=DDT_None; DDTValue=0; DDTArray=false;
  ShiftingMode=(Shifting? Shifting->GetShiftMode(): SHIFT_None);
  Visco=0; ViscoBoundFactor=1;
  TBoundary=BC_DBC;
  SlipMode=SLIP_Vel0;  //<vs_mddbc>
  MdbcCorrector=false; //<vs_mddbc>
  UseNormals=false;    //<vs_mddbc>
  UseNormalsFt=false;  //<vs_mddbc>
  UseDEM=false;  //(DEM)
  delete[] DemData; DemData=NULL;  //(DEM)
  UseChrono=false; //<vs_chroono>
  RhopOut=true; RhopOutMin=700; RhopOutMax=1300;
  TimeMax=TimePart=0;
  NstepsBreak=0;
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
  DDT2h=DDTgz=0;
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

  WrnPartsOut=true;

  FtCount=0;
  FtPause=0;
  FtMode=FTMODE_None;
  FtConstraints=false;
  FtIgnoreRadius=false;
  WithFloating=false;

  AllocMemoryFloating(0,false);

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
/// Loads the case configuration to be executed.
//==============================================================================
void JSph::ConfigDomainResize(std::string key,const JSpaceEParms *eparms){
  const char met[]="ConfigDomainResize";
  const char axis=fun::StrLower(key)[0];
  if(axis!='x' && axis!='y' && axis!='z')RunException(met,"Axis value is invalid.");
  if(key.substr(1,3)!="min" && key.substr(1,3)!="max")RunException(met,"Key value is invalid.");
  if(key.substr(1,3)=="min"){
    JSpaceEParms::JSpaceEParmsPos ps=eparms->GetPosminValue(axis);
    switch(ps.mode){
      case JSpaceEParms::DC_Fixed:     ConfigDomainFixedValue(string("DomainFixed")+key,ps.value);                     break;
      case JSpaceEParms::DC_DefValue:  ConfigDomainParticlesValue(string("DomainParticles")+key,ps.value);             break;
      case JSpaceEParms::DC_DefPrc:    ConfigDomainParticlesPrcValue(string("DomainParticlesPrc")+key,ps.value/100);   break;
    }
  }
  else{
    JSpaceEParms::JSpaceEParmsPos ps=eparms->GetPosmaxValue(axis);
    switch(ps.mode){
      case JSpaceEParms::DC_Fixed:     ConfigDomainFixedValue(string("DomainFixed")+key,ps.value);                     break;
      case JSpaceEParms::DC_DefValue:  ConfigDomainParticlesValue(string("DomainParticles")+key,ps.value);             break;
      case JSpaceEParms::DC_DefPrc:    ConfigDomainParticlesPrcValue(string("DomainParticlesPrc")+key,ps.value/100);   break;
    }
  }
}

//==============================================================================
/// Allocates memory of floating objectcs.
//==============================================================================
void JSph::AllocMemoryFloating(unsigned ftcount,bool imposedvel){
  delete[] FtObjs; FtObjs=NULL;
  if(FtLinearVel){ //<vs_fttvel_ini>
    for(unsigned c=0;c<FtCount;c++)delete FtLinearVel[c];
    delete[] FtLinearVel; FtLinearVel=NULL;
  }
  if(FtAngularVel){
    for(unsigned c=0;c<FtCount;c++)delete FtAngularVel[c];
    delete[] FtAngularVel; FtAngularVel=NULL;
  } //<vs_fttvel_end>
  if(ftcount){
    FtObjs=new StFloatingData[ftcount];
    if(imposedvel){ //<vs_fttvel_ini>
      FtLinearVel =new JLinearValue*[ftcount];
      FtAngularVel=new JLinearValue*[ftcount];
      for(unsigned c=0;c<ftcount;c++)FtLinearVel[c]=FtAngularVel[c]=NULL;
    } //<vs_fttvel_end>
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
  if(DtFixed)s+=DtFixed->GetAllocMemory();
  if(AccInput)s+=AccInput->GetAllocMemory();
  if(PartsLoaded)s+=PartsLoaded->GetAllocMemory();
  return(s);
}

//==============================================================================
/// Loads the configuration of the execution.
//==============================================================================
void JSph::LoadConfig(const JCfgRun *cfg){
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
}

//==============================================================================
/// Loads predefined constans from XML.
//==============================================================================
void JSph::LoadConfigCtes(const JXml *xml){
  JSpaceCtes ctes;
  ctes.LoadXmlRun(xml,"case.execution.constants");

  Simulate2D=ctes.GetData2D();
  Simulate2DPosY=ctes.GetData2DPosY();
  H=(float)ctes.GetH();
  CteB=(float)ctes.GetB();
  Gamma=(float)ctes.GetGamma();
  RhopZero=(float)ctes.GetRhop0();
  CFLnumber=(float)ctes.GetCFLnumber();
  Dp=ctes.GetDp();
  Gravity=ToTFloat3(ctes.GetGravity());
  MassFluid=(float)ctes.GetMassFluid();
  MassBound=(float)ctes.GetMassBound();
  if(ctes.GetEps()!=0)Log->PrintWarning("Eps value is not used (this correction is deprecated).");
}

//==============================================================================
/// Loads execution parameters from XML.
//==============================================================================
void JSph::LoadConfigParameters(const JXml *xml){
  JSpaceEParms eparms;
  eparms.LoadXml(xml,"case.execution.parameters");

  if(eparms.Exists("FtSaveAce"))SaveFtAce=(eparms.GetValueInt("FtSaveAce",true,0)!=0);

  if(eparms.Exists("PosDouble")){
    Log->PrintWarning("The parameter \'PosDouble\' is deprecated.");
    SvPosDouble=(eparms.GetValueInt("PosDouble")==2);
  }
  if(eparms.Exists("SavePosDouble"))SvPosDouble=(eparms.GetValueInt("SavePosDouble",true,0)!=0);
  switch(eparms.GetValueInt("RigidAlgorithm",true,1)){ //(DEM)
    case 1:  UseDEM=false;    break;
    case 2:  UseDEM=true;     break;
    case 3:  UseChrono=true;  break;  //<vs_chroono>
    default: Run_Exceptioon("Rigid algorithm is not valid.");
  }
  switch(eparms.GetValueInt("StepAlgorithm",true,1)){
    case 1:  TStep=STEP_Verlet;      break;
    case 2:  TStep=STEP_Symplectic;  break;
    default: Run_Exceptioon("Step algorithm is not valid.");
  }
  VerletSteps=eparms.GetValueInt("VerletSteps",true,40);
  switch(eparms.GetValueInt("Kernel",true,2)){
    case 1:  TKernel=KERNEL_Cubic;     break;
    case 2:  TKernel=KERNEL_Wendland;  break;
    case 3:  TKernel=KERNEL_Gaussian;  break;
#ifdef PRASS3_WENDLANDC6                        //<vs_praticalsskq>
    case 4:  TKernel=KERNEL_WendlandC6; break;  //<vs_praticalsskq>
#endif                                          //<vs_praticalsskq>
    default: Run_Exceptioon("Kernel choice is not valid.");
  }
  switch(eparms.GetValueInt("ViscoTreatment",true,1)){
    case 1:  TVisco=VISCO_Artificial;  break;
    case 2:  TVisco=VISCO_LaminarSPS;  break;
    default: Run_Exceptioon("Viscosity treatment is not valid.");
  }
  Visco=eparms.GetValueFloat("Visco");
  ViscoBoundFactor=eparms.GetValueFloat("ViscoBoundFactor",true,1.f);
  string filevisco=eparms.GetValueStr("ViscoTime",true);
  if(!filevisco.empty()){
    ViscoTime=new JSphVisco();
    ViscoTime->LoadFile(DirCase+filevisco);
  }

  //-Boundary configuration.
  switch(eparms.GetValueInt("Boundary",true,1)){
    case 1:  TBoundary=BC_DBC;      break;
    case 2:  TBoundary=BC_MDBC;     break;   //<vs_mddbc>
    default: Run_Exceptioon("Boundary Condition method is not valid.");
  }
  if(TBoundary==BC_MDBC){ //<vs_mddbc_ini>
    UseNormals=true;
    switch(eparms.GetValueInt("SlipMode",true,1)){
      case 1:  SlipMode=SLIP_Vel0;      break;
      case 2:  SlipMode=SLIP_NoSlip;    break;
      case 3:  SlipMode=SLIP_FreeSlip;  break;
      default: Run_Exceptioon("Slip mode is not valid.");
    }
    MdbcCorrector=(eparms.GetValueInt("MDBCCorrector",true,0)!=0);
  } 
  if(SlipMode==SLIP_FreeSlip)Run_Exceptioon("SlipMode=\'Free slip\' is not available for now.");
  //<vs_mddbc_end>

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
      case 2:  TDensity=DDT_DDT2;     break;  //<vs_dtt2>
      case 3:  TDensity=DDT_DDT2Full; break;  //<vs_dtt2>
      default: 
        if(tddt==2 || tddt==3)Run_Exceptioon("Density Diffusion Term \'Fourtakas et al 2019\' is not available in this version.");
        else Run_Exceptioon("Density Diffusion Term mode is not valid.");
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
    Shifting=new JShifting(Simulate2D,Dp,H,Log);
    Shifting->ConfigBasic(shiftmode,shiftcoef,shifttfs);
  }

  WrnPartsOut=(eparms.GetValueInt("WrnPartsOut",true,1)!=0);
  FtPause=eparms.GetValueFloat("FtPause",true,0);
  FtIgnoreRadius=(eparms.GetValueInt("FtIgnoreRadius",true,0)!=0);

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
void JSph::LoadConfigCommands(const JCfgRun *cfg){
  //-Aplies configuration using command line.
  if(cfg->SvPosDouble>=0)SvPosDouble=(cfg->SvPosDouble!=0);
  if(cfg->TBoundary){
    TBoundary=BC_DBC;  SlipMode=SLIP_Vel0;  MdbcCorrector=false;
    switch(cfg->TBoundary){
      case 1:  TBoundary=BC_DBC;                          break;
      case 2:  TBoundary=BC_MDBC;  SlipMode=SLIP_Vel0;    break;
      case 3:  TBoundary=BC_MDBC;  SlipMode=SLIP_NoSlip;  break;
      default: Run_Exceptioon("Boundary method is not valid.");
    }
    UseNormals=(TBoundary==BC_MDBC);
  }
  if(TBoundary==BC_MDBC && SlipMode!=SLIP_Vel0)Run_Exceptioon("Only the slip mode velocity=0 is allowed with mDBC conditions.");
  if(cfg->TStep)TStep=cfg->TStep;
  if(cfg->VerletSteps>=0)VerletSteps=cfg->VerletSteps;
  if(cfg->TKernel)TKernel=cfg->TKernel;
  if(cfg->TVisco){ TVisco=cfg->TVisco; Visco=cfg->Visco; }
  if(cfg->ViscoBoundFactor>=0)ViscoBoundFactor=cfg->ViscoBoundFactor;

  //-Density Diffusion Term configuration.
  if(cfg->TDensity>=0){
    switch(cfg->TDensity){
      case 0:  TDensity=DDT_None;     break;
      case 1:  TDensity=DDT_DDT;      break;
      case 2:  TDensity=DDT_DDT2;     break;  //<vs_dtt2>
      case 3:  TDensity=DDT_DDT2Full; break;  //<vs_dtt2>
      default: 
        if(cfg->TDensity==2 || cfg->TDensity==3)Run_Exceptioon("Density Diffusion Term \'Fourtakas et al 2019\' is not available in this version.");
        else Run_Exceptioon("Density Diffusion Term mode is not valid.");
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
    if(!Shifting)Shifting=new JShifting(Simulate2D,Dp,H,Log);
    Shifting->ConfigBasic(shiftmode);
  }

  if(cfg->FtPause>=0)FtPause=cfg->FtPause;
  if(cfg->TimeMax>0)TimeMax=cfg->TimeMax;
  NstepsBreak=cfg->NstepsBreak;
  if(NstepsBreak)Log->PrintfWarning("The execution will be cancelled after %d simulation steps.",NstepsBreak);
  //-Configuration of JTimeOut with TimePart.
  TimeOut=new JTimeOut();
  if(cfg->TimePart>=0){
    TimePart=cfg->TimePart;
    TimeOut->Config(TimePart);
  }
  else TimeOut->Config(FileXml,"case.execution.special.timeout",TimePart);

  CellMode=cfg->CellMode;
  if(cfg->DomainMode==2)ConfigDomainFixed(cfg->DomainFixedMin,cfg->DomainFixedMax);
  if(cfg->RhopOutModif){
    RhopOutMin=cfg->RhopOutMin; RhopOutMax=cfg->RhopOutMax;
  }
  RhopOut=(RhopOutMin<RhopOutMax);
  if(!RhopOut){ RhopOutMin=-FLT_MAX; RhopOutMax=FLT_MAX; }
}

//==============================================================================
/// Creates and load NuxLib object to evaluate user-defined expressions.
//==============================================================================
void JSph::LoadConfigVars(const JXml *xml){
  if(!JNumexLib::Available())Log->PrintWarning("Code for JNumex libary is not included in the current compilation, so user-defined expresions in XML file are not evaluated.");
  else{
    NuxLib=new JNumexLib();
    //-Loads user variables from XML.
    JSpaceUserVars uvars;
    uvars.LoadXml(xml,"case.execution.uservars",true);
    for(unsigned c=0;c<uvars.CountVars();c++){
      const JSpaceUserVars::StVar v=uvars.GetVar(c);
      if(v.isnum)NuxLib->CreateVar(v.name,true,false,v.valuenum);
      else       NuxLib->CreateVar(v.name,true,false,v.valuestr);
    }
    //-Defines basic constant values on NuxLib object.
    bool rep=false;
    rep|=NuxLib->CreateVar("CaseName"  ,true,true,CaseName);
    rep|=NuxLib->CreateVar("Data2D"    ,true,true,Simulate2D);
    rep|=NuxLib->CreateVar("Data2DPosy",true,true,Simulate2DPosY);
    rep|=NuxLib->CreateVar("H"         ,true,true,H);
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
void JSph::LoadCaseConfig(const JCfgRun *cfg){
  if(!fun::FileExists(FileXml))Run_ExceptioonFile("Case configuration was not found.",FileXml);
  JXml xml; xml.LoadFile(FileXml);

  //-Predefined constantes.
  LoadConfigCtes(&xml);

  //-Define variables on NuxLib object.
  LoadConfigVars(&xml);
  //-Enables the use of NuxLib in XML configuration.
  xml.SetNuxLib(NuxLib);

  //-Execution parameters from XML.
  LoadConfigParameters(&xml);

  //-Execution parameters from commands.
  LoadConfigCommands(cfg);

  //-Define variables of execution parameters and shows list of available variables.
  LoadConfigVarsExec();

  //-Particle data.
  JSpaceParts parts;
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
  GaugeSystem=new JGaugeSystem(Cpu,Log);

  //-Configuration of AccInput.
  if(xml.GetNodeSimple("case.execution.special.accinputs",true)){
    AccInput=new JSphAccInput(Log,DirCase,&xml,"case.execution.special.accinputs");
  }

  //-Configuration of ChronoObjects. //<vs_chroono_ini>
  if(UseChrono){
    if(xml.GetNodeSimple("case.execution.special.chrono",true)){
      if(!JChronoObjects::Available())Run_Exceptioon("DSPHChronoLib to use Chrono is not included in the current compilation.");
      ChronoObjects=new JChronoObjects(Log,DirCase,CaseName,&xml,"case.execution.special.chrono",Dp,parts.GetMkBoundFirst());
    }
    else Run_ExceptioonFile("Chrono configuration in XML file is missing.",FileXml);
  }
  else if(xml.GetNodeSimple("case.execution.special.chrono",true))Log->PrintfWarning("The use of Chrono is disabled (RigidAlgorithm is not 3) although XML file includes the Chrono configuration.");
  //<vs_chroono_end>

  //-Loads and configures moving objects.
  if(parts.CountBlocks(TpPartMoving)>0){
    SphMotion=new JSphMotion(Simulate2D);
    SphMotion->Init(&parts,&xml,"case.execution.motion",DirCase);
  }

  //-Configuration of WaveGen.
  if(xml.GetNodeSimple("case.execution.special.wavepaddles",true)){
    if(!JWaveGen::Available())Run_Exceptioon("Code for Wave-Paddles is not included in the current compilation.");
    bool useomp=false,usegpu=false;
    #ifdef OMP_USE_WAVEGEN
      useomp=(omp_get_max_threads()>1);
    #endif
    WaveGen=new JWaveGen(useomp,!Cpu,Log,DirCase,&xml,"case.execution.special.wavepaddles",ToTDouble3(Gravity));
    if(SphMotion)for(unsigned ref=0;ref<SphMotion->GetNumObjects();ref++){
      const StMotionData& m=SphMotion->GetMotionData(ref);
      WaveGen->ConfigPaddle(m.mkbound,ref,m.idbegin,m.count);
    }
  }

  //-Configuration of MLPistons.  //<vs_mlapiston_ini>
  if(xml.GetNodeSimple("case.execution.special.mlayerpistons",true)){
    if(!JMLPistons::Available())Run_Exceptioon("Code for Multi-Layer Pistons is not included in the current compilation.");
    bool useomp=false,usegpu=false;
    MLPistons=new JMLPistons(!Cpu,Log,DirCase);
    MLPistons->LoadXml(&xml,"case.execution.special.mlayerpistons");  
    if(SphMotion)for(unsigned c=0;c<SphMotion->GetNumObjects();c++){
      const StMotionData& m=SphMotion->GetMotionData(c);
      MLPistons->ConfigPiston(m.mkbound,c,m.idbegin,m.count,TimeMax);
    }
    MLPistons->CheckPistons();
  }//<vs_mlapiston_end>

  //-Configuration of RelaxZones.  //<vs_rzone_ini>
  if(xml.GetNodeSimple("case.execution.special.relaxationzones",true)){
    if(!JRelaxZones::Available())Run_Exceptioon("Code for Relaxation Zones is not included in the current compilation.");
    bool useomp=false,usegpu=false;
    #ifdef OMP_USE_WAVEGEN
      useomp=(omp_get_max_threads()>1);
    #endif
    RelaxZones=new JRelaxZones(useomp,!Cpu,Log,DirCase,CaseNfloat>0,CaseNbound,ToTDouble3(Gravity));
    RelaxZones->LoadXml(&xml,"case.execution.special.relaxationzones");
  }//<vs_rzone_end>

  //-Configuration of Shifting with zones.
  if(Shifting && xml.GetNodeSimple("case.execution.special.shifting",true))Run_ExceptioonFile("Shifting is defined several times (in <special><shifting> and <execution><parameters>).",FileXml);
  if(xml.GetNodeSimple("case.execution.special.shifting",true)){
    Shifting=new JShifting(Simulate2D,Dp,H,Log);
    Shifting->LoadXml(&xml,"case.execution.special.shifting");
  }
  if(Shifting && !Shifting->GetShiftMode()){ delete Shifting; Shifting=NULL; }
  ShiftingMode=(Shifting? Shifting->GetShiftMode(): SHIFT_None);

  //-Configuration of damping zones.
  if(xml.GetNodeSimple("case.execution.special.damping",true)){
    Damping=new JDamping(Dp,Log);
    Damping->LoadXml(&xml,"case.execution.special.damping");
  }

  //-Loads floating objects.
  FtCount=parts.CountBlocks(TpPartFloating);
  if(FtCount){
    FtConstraints=false;
    if(FtCount>CODE_MKRANGEMAX)Run_Exceptioon("The number of floating objects exceeds the maximum.");
    AllocMemoryFloating(FtCount,parts.UseImposedFtVel());
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
        fobj->massp=(float)fblock.GetMasspart();
        fobj->radius=0;
        fobj->constraints=ComputeConstraintsValue(fblock.GetTranslationFree(),fblock.GetRotationFree());
        if(fobj->constraints!=FTCON_Free)FtConstraints=true;
        if(fblock.GetLinearVel()){ //<vs_fttvel_ini>
          FtLinearVel[cobj]=new JLinearValue(*fblock.GetLinearVel());
          if(!FtLinearVel[cobj]->GetFile().empty())FtLinearVel[cobj]->LoadFile(DirCase+FtLinearVel[cobj]->GetFile());
        }
        if(fblock.GetAngularVel()){
          FtAngularVel[cobj]=new JLinearValue(*fblock.GetAngularVel());
          if(!FtAngularVel[cobj]->GetFile().empty())FtAngularVel[cobj]->LoadFile(DirCase+FtAngularVel[cobj]->GetFile());
        } //<vs_fttvel_end>
        fobj->center=fblock.GetCenter();
        fobj->angles=TFloat3(0);
        fobj->fvel=ToTFloat3(fblock.GetLinearVelini());
        fobj->fomega=ToTFloat3(fblock.GetAngularVelini());
        fobj->inertiaini=ToTMatrix3f(fblock.GetInertia());
        //-Chrono configuration. //<vs_chroono_ini> 
        fobj->usechrono=(ChronoObjects && ChronoObjects->ConfigBodyFloating(fblock.GetMkType()
          ,fblock.GetMassbody(),fblock.GetCenter(),fblock.GetInertia()
          ,fblock.GetTranslationFree(),fblock.GetRotationFree(),fobj->fvel,fobj->fomega));
        //<vs_chroono_end>
        cobj++;
      }
    }
  }
  else{
    if(UseDEM){
      UseDEM=false;
      Log->PrintWarning("The use of DEM was disabled because there are no floating objects...");
    }
    if(UseChrono){ //<vs_chroono_ini> 
      UseChrono=false;
      Log->PrintWarning("The use of CHRONO was disabled because there are no floating objects...");
      delete ChronoObjects; ChronoObjects=NULL;
    } //<vs_chroono_end>
  }
  WithFloating=(FtCount>0);
  FtMode=(WithFloating? FTMODE_Sph: FTMODE_None);
  if(UseDEM)FtMode=FTMODE_Ext;
  if(UseChrono)FtMode=FTMODE_Ext; //<vs_chroono>
  if(UseChrono && PeriActive!=0)Log->PrintfWarning("The use of Chrono with periodic limits is only recommended when moving and floating objects do not move beyond those periodic limits."); //<vs_chroono>

  //-Loads DEM and DVI data for boundary objects.   //<vs_chroono>
  if(UseDEM || UseChrono){/*                        //<vs_chroono>
  //-Loads DEM data for boundary objects. (DEM)
  if(UseDEM){
  */                                                //<vs_chroono>
    if(UseDEM){
      DemData=new StDemData[DemDataSize];
      memset(DemData,0,sizeof(StDemData)*DemDataSize);
    }
    for(unsigned c=0;c<parts.CountBlocks();c++){
      const JSpacePartBlock &block=parts.GetBlock(c);
      if(IsBound(block.Type)){
        const word mkbound=block.GetMkType();
        const unsigned cmk=MkInfo->GetMkBlockByMkBound(mkbound);
        if(cmk>=MkInfo->Size())Run_Exceptioon(fun::PrintStr("Error loading boundary objects. Mkbound=%u is unknown.",mkbound));
        //-Loads data for collisions using DEM.
        if(UseDEM){
          const unsigned tav=CODE_GetTypeAndValue(MkInfo->Mkblock(cmk)->Code);
          DemData[tav]=LoadDemData(true,&block);
        }
        //-Loads data for collisions using Chrono. //<vs_chroono_ini>
        if(ChronoObjects && ChronoObjects->UseDataDVI(mkbound)){
          const StDemData data=LoadDemData(false,&block);
          if(block.Type==TpPartFloating)ChronoObjects->ConfigDataBodyFloating(mkbound,data.kfric,data.restitu,data.young,data.poisson);
          if(block.Type==TpPartMoving)  ChronoObjects->ConfigDataBodyMoving  (mkbound,data.kfric,data.restitu,data.young,data.poisson);
          if(block.Type==TpPartFixed)   ChronoObjects->ConfigDataBodyFixed   (mkbound,data.kfric,data.restitu,data.young,data.poisson);
        } //<vs_chroono_end>
      }
    }
  }

  //-Configuration of Inlet/Outlet.  //<vs_innlet_ini> 
  if(xml.GetNodeSimple("case.execution.special.inout",true)){
    InOut=new JSphInOut(Cpu,Log,FileXml,&xml,"case.execution.special.inout",DirCase);
    NpDynamic=true;
    ReuseIds=InOut->GetReuseIds();
    if(ReuseIds)Run_Exceptioon("Inlet/Outlet with ReuseIds is not a valid option for now...");
  }
  
  //-Configuration of boundary extrapolated correction.
  if(xml.GetNodeSimple("case.execution.special.boundextrap",false))Run_Exceptioon("The XML section 'boundextrap' is obsolete.");
  if(xml.GetNodeSimple("case.execution.special.boundcorr",true)){
    BoundCorr=new JSphBoundCorr(Cpu,Dp,Log,&xml,"case.execution.special.boundcorr",MkInfo);
  } //<vs_innlet_end> 
 
  //-Configuration of Moorings object. //<vs_moordyyn_ini>
  if(xml.GetNodeSimple("case.execution.special.moorings",true)){
    if(WithFloating){
      if(!AVAILABLE_MOORDYN)Run_Exceptioon("Code for moorings and MoorDyn+ coupling is not included in the current compilation.");
      Moorings=new JMooredFloatings(Log,DirCase,CaseName,Gravity);
      Moorings->LoadXml(&xml,"case.execution.special.moorings");
    }
    else Log->PrintWarning("The use of Moorings was disabled because there are no floating objects...");
  }

  //-Configuration of FtForces object.
  if(Moorings || xml.GetNodeSimple("case.execution.special.forcepoints",true)){
    if(WithFloating){
      ForcePoints=new JSphFtForcePoints(Log,Cpu,Dp,FtCount);
      //FtForces->LoadXml(&xml,"case.execution.special.forcepoints");
    }
    else Log->PrintWarning("The use of impossed force to floatings was disabled because there are no floating objects...");
  } //<vs_moordyyn_end>

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
StDemData JSph::LoadDemData(bool checkdata,const JSpacePartBlock* block)const{
  const word mk=block->GetMk();
  const word mkbound=block->GetMkType();
  StDemData data;
  memset(&data,0,sizeof(StDemData));
  //-Checks necessary values for collisions using DEM or Chrono.
  if(checkdata){
    const string objdesc=fun::PrintStr("Object mk=%u (mkbound=%u)",mk,mkbound);
    if(!block->ExistsSubValue("Kfric","value"))Run_Exceptioon(objdesc+" - Value of Kfric is invalid.");
    if(!block->ExistsSubValue("Restitution_Coefficient","value"))Run_Exceptioon(objdesc+" - Value of Restitution_Coefficient is invalid.");
    if(!block->ExistsSubValue("Young_Modulus","value"))Run_Exceptioon(objdesc+" - Value of Young_Modulus is invalid.");
    if(!block->ExistsSubValue("PoissonRatio","value"))Run_Exceptioon(objdesc+" - Value of PoissonRatio is invalid.");
  }
  //-Loads necessary values for collisions using DEM or Chrono.
  data.kfric=block->GetSubValueFloat("Kfric","value",true,FLT_MAX);
  data.restitu=block->GetSubValueFloat("Restitution_Coefficient","value",true,FLT_MAX);
  if(block->ExistsValue("Restitution_Coefficient_User"))data.restitu=block->GetValueFloat("Restitution_Coefficient_User");
  data.young=block->GetSubValueFloat("Young_Modulus","value",true,FLT_MAX);
  data.poisson=block->GetSubValueFloat("PoissonRatio","value",true,FLT_MAX);
  //-Loads necessary values for DEM.
  data.massp=MassBound;
  if(block->Type==TpPartFloating){
    const JSpacePartBlock_Floating *fblock=(const JSpacePartBlock_Floating *)block;
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
  const char met[]="LoadCodeParticles"; 
  //-Assigns code to each group of particles.
  for(unsigned p=0;p<np;p++)code[p]=MkInfo->GetCodeById(idp[p]);
}

//<vs_mddbc_ini>
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
      nmarrone.RunCase(DirCase+CaseName,true);
      pnor=nmarrone.GetPartNor();
      pnorsize=nmarrone.GetPartNorSize();
    }
    else{
      //-Loads final normals in NBI4 file.
      pnor=nd.GetPartNormals();
      if(pnor==NULL)Run_ExceptioonFile("Normal data and final normals are missing in NBI4 file.",filenordata);
      pnorsize=nd.GetNbound();
    }
    //-Applies final normals. Loads normals from boundary particle to boundary limit.
    if(pnorsize<npb)Run_ExceptioonFile("The number of final normals does not match fixed and moving particles.",filenordata);
    for(unsigned p=0;p<npb;p++)boundnormal[p]=ToTFloat3(pnor[p]);  //-For fixed and moving particles.
    Log->Printf("NormalDataFile=\"%s\"",filenordata.c_str());
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
//<vs_mddbc_end>

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
  //-Symmetry domain configuration. //<vs_syymmetry>
  if(Symmetry)MapRealPosMin.y=0;    //<vs_syymmetry>
}

//==============================================================================
/// Calculate Wendland constants.
//==============================================================================
void JSph::WendlandConstants(bool simulate2d,double h,float &awen,float &bwen){
  if(simulate2d){
    awen=float(0.557/(h*h));
    bwen=float(-2.7852/(h*h*h));
  }
  else{
    awen=float(0.41778/(h*h*h));
    bwen=float(-2.08891/(h*h*h*h));
  }
}

//==============================================================================
/// Configures value of constants.
//==============================================================================
void JSph::ConfigConstants(bool simulate2d){
  const char* met="ConfigConstants";
  //-Computation of constants.
  const double h=H;
  DDT2h=float(h*2*DDTValue);
  DDTgz=float(double(RhopZero)*double(fabs(Gravity.z))/double(CteB));
  const double cs0=sqrt(double(Gamma)*double(CteB)/double(RhopZero));
  Cs0=cs0;
  if(!DtIni)DtIni=h/cs0;
  if(!DtMin)DtMin=(h/cs0)*CoefDtMin;
  Dosh=float(h*2); 
  H2=float(h*h);
  Fourh2=float(h*h*4); 
  Eta2=float((h*0.1)*(h*0.1));
  //-Wendland constants are always computed since this kernel is used in some parts where other kernels are not defined (e.g. mDBC, inlet/outlet, boundcorr...).
  WendlandConstants(simulate2d,H,Awen,Bwen);
  //-Computes constants for other kernels.
  if(TKernel==KERNEL_Cubic){
    if(simulate2d){
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
    else{
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
  else if(TKernel==KERNEL_Gaussian){
    if(simulate2d){
      const double a1=4./PI;
      const double a2=a1/(h*h);
      const double aa=a1/(h*h*h);
      Agau=float(a2);
      Bgau=float(-8.*aa);
    }
    else{
      const double a1=8./5.5683;
      const double a2=a1/(h*h*h);
      const double aa=a1/(h*h*h*h); 
      Agau=float(a2);
      Bgau=float(-8.*aa);
    }
  }
  else if(TKernel==KERNEL_WendlandC6){  //<vs_praticalsskq_ini>
    if(simulate2d){
      Awc6=float(78./(28.*PI*h*h));
      Bwc6=float(Awc6/h);
    }
    else{
      Awc6=float(1365./(512.*PI*h*h*h));
      Bwc6=float(Awc6/h);
    }
  }  //<vs_praticalsskq_end>
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
  Log->Print(Simulate2D? "**2D-Simulation parameters:": "**3D-Simulation parameters:");
  Log->Print(fun::VarStr("CaseName",CaseName));
  Log->Print(fun::VarStr("RunName",RunName));
  if(Simulate2D)Log->Print(fun::VarStr("Simulate2DPosY",Simulate2DPosY));
  Log->Print(fun::VarStr("Symmetry",Symmetry));  //<vs_syymmetry>
  Log->Print(fun::VarStr("SavePosDouble",SvPosDouble));
  Log->Print(fun::VarStr("SaveFtAce",SaveFtAce));
  Log->Print(fun::VarStr("SvTimers",SvTimers));
  Log->Print(fun::VarStr("Boundary",GetBoundName(TBoundary)));
  if(TBoundary==BC_MDBC){ //<vs_mddbc_ini>
    Log->Print(fun::VarStr("  SlipMode",GetSlipName(SlipMode)));
    Log->Print(fun::VarStr("  mDBC-Corrector",MdbcCorrector));
  } //<vs_mddbc_end>
  Log->Print(fun::VarStr("StepAlgorithm",GetStepName(TStep)));
  if(TStep==STEP_None)Run_Exceptioon("StepAlgorithm value is invalid.");
  if(TStep==STEP_Verlet)Log->Print(fun::VarStr("  VerletSteps",VerletSteps));
  Log->Print(fun::VarStr("Kernel",GetKernelName(TKernel)));
  Log->Print(fun::VarStr("Viscosity",GetViscoName(TVisco)));
  Log->Print(fun::VarStr("  Visco",Visco));
  Log->Print(fun::VarStr("  ViscoBoundFactor",ViscoBoundFactor));
  if(ViscoTime)Log->Print(fun::VarStr("ViscoTime",ViscoTime->GetFile()));
  Log->Print(fun::VarStr("DensityDiffusion",GetDDTName(TDensity)));
  if(TDensity!=DDT_None){
    Log->Print(fun::VarStr("  DensityDiffusionValue",DDTValue));
    //Log->Print(fun::VarStr("DensityDiffusionArray",DDTArray));
  }
  if(TDensity==DDT_DDT2Full && H/Dp>1.5)Log->PrintWarning("The selected DDT \'(Fourtakas et al 2019 (full)\' needs several boundary layers when h/dp>1.5");  //<vs_dtt2>
  if(Shifting)Shifting->VisuConfig();
  else Log->Print(fun::VarStr("Shifting","None"));
  string rigidalgorithm=(!FtCount? "None": (UseDEM? "SPH+DCDEM": "SPH"));
  if(UseChrono)rigidalgorithm="SPH+CHRONO"; //<vs_chroono>
  Log->Print(fun::VarStr("RigidAlgorithm",rigidalgorithm));
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
  else if(TKernel==KERNEL_Cubic){
    Log->Print(fun::VarStr("CubicCte.a1",CubicCte.a1));
    Log->Print(fun::VarStr("CubicCte.aa",CubicCte.aa));
    Log->Print(fun::VarStr("CubicCte.a24",CubicCte.a24));
    Log->Print(fun::VarStr("CubicCte.c1",CubicCte.c1));
    Log->Print(fun::VarStr("CubicCte.c2",CubicCte.c2));
    Log->Print(fun::VarStr("CubicCte.d1",CubicCte.d1));
    Log->Print(fun::VarStr("CubicCte.od_wdeltap",CubicCte.od_wdeltap));
  }
  else if(TKernel==KERNEL_Gaussian){
    Log->Print(fun::VarStr("Agau (Gaussian)",Agau));
    Log->Print(fun::VarStr("Bgau (Gaussian)",Bgau));
  }
  else if(TKernel==KERNEL_WendlandC6){ //<vs_praticalsskq_ini>
    Log->Print(fun::VarStr("Awc6 (WendlandC6)",Awc6));
    Log->Print(fun::VarStr("Bwc6 (WendlandC6)",Bwc6));
  } //<vs_praticalsskq_end>
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
  Log->Print(fun::VarStr("WrnPartsOut",WrnPartsOut));
  if(CteB==0)Run_Exceptioon("Constant \'b\' cannot be zero.\n\'b\' is zero when fluid height is zero (or fluid particles were not created)");
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
    JXml xml; xml.LoadFile(FileXml);
    xml.SetNuxLib(NuxLib); //-Enables the use of NuxLib in XML configuration.
    if(xml.GetNodeSimple("case.execution.special.initialize",true)){
      JSphInitialize init(&xml,"case.execution.special.initialize",H,boundnormal!=NULL);
      if(init.Count()){
        //-Creates array with mktype value.
        word *mktype=new word[np];
        for(unsigned p=0;p<np;p++){
          const unsigned cmk=MkInfo->GetMkBlockByCode(code[p]);
          mktype[p]=(cmk<MkInfo->Size()? word(MkInfo->Mkblock(cmk)->MkType): USHRT_MAX);
        }
        init.Run(np,npb,CaseNbound,pos,idp,mktype,velrhop,boundnormal);
        init.GetConfig(InitializeInfo);
        //-Frees memory.
        delete[] mktype; mktype=NULL;
      }
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
  if(!Cpu){
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
    if(Map_Cells.x>nx || Map_Cells.y>ny || Map_Cells.z>nz)Run_Exceptioon("PosCellCode is invalid for MapCells.");
  }
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
  const string errtex="The floating body radius (%g [m]) is too large for periodic distance in %c (%g [m]). If the floating body crosses the periodical limits, the simulation may be incorrect.";
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
  PartsLoaded->CheckConfig(CaseNp,CaseNfixed,CaseNmoving,CaseNfloat,CaseNfluid,Simulate2D,Simulate2DPosY,TpPeri(PeriActive));
  if(PartBegin)RestartCheckData();
  Log->Printf("Loaded particles: %u",PartsLoaded->GetCount());

  //-Collect information of loaded particles.
  //-Recupera informacion de las particulas cargadas.
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
  if(SphMotion){
    SphMotion->SetTimeMod(!PartIni? PartBeginTimeStep: 0);
    SphMotion->ProcesTime(JSphMotion::MOMT_Simple,0,TimeStepIni);
  }
  //-Adjust motion paddles for the instant of the loaded PART.
  if(WaveGen)WaveGen->SetTimeMod(!PartIni? PartBeginTimeStep: 0);
  //-Adjust Multi-layer pistons for the instant of the loaded PART.   //<vs_mlapiston>
  if(MLPistons)MLPistons->SetTimeMod(!PartIni? PartBeginTimeStep: 0);     //<vs_mlapiston>

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
  GaugeSystem->Config(Simulate2D,Simulate2DPosY,Symmetry,TimeMax,TimePart
    ,Dp,DomPosMin,DomPosMax,Scell,Hdiv,H,MassFluid,MassBound,float(Cs0),CteB,Gamma,RhopZero);
  if(xml.GetNodeSimple("case.execution.special.gauges",true))
    GaugeSystem->LoadXml(&xml,"case.execution.special.gauges",MkInfo);

  //-Prepares WaveGen configuration.
  if(WaveGen){
    Log->Print("Wave paddles configuration:");
    WaveGen->Init(GaugeSystem,MkInfo,TimeMax,TimePart);
    WaveGen->VisuConfig(""," ");
  }

  //-Prepares MLPistons configuration.  //<vs_mlapiston_ini>
  if(MLPistons){
    Log->Printf("Multi-Layer Pistons configuration:");
    //if(MLPistons->UseAwasZsurf())MLPistons->InitAwas(Gravity,Simulate2D,CellOrder,MassFluid,Dp,Dosh,Scell,Hdiv,DomPosMin,DomRealPosMin,DomRealPosMax); ->awaspdte
    MLPistons->VisuConfig(""," ");
  }  //<vs_mlapiston_end>

  //-Prepares RelaxZones configuration.  //<vs_rzone_ini>
  if(RelaxZones){
    Log->Print("Relaxation Zones configuration:");
    RelaxZones->Init(DirCase,TimeMax,Dp);
    RelaxZones->VisuConfig(""," ");
  }  //<vs_rzone_end>

  //-Prepares ChronoObjects configuration.  //<vs_chroono_ini>
  if(ChronoObjects){
    Log->Print("Chrono Objects configuration:");
    if(PartBegin)RunException(met,"Simulation restart not allowed when Chrono is used.");
    ChronoObjects->Init(Simulate2D,MkInfo);
    ChronoObjects->VisuConfig(""," ");
  }  //<vs_chroono_end>

  //-Prepares Moorings configuration.  //<vs_moordyyn_ini>
  if(Moorings){
    Log->Print("Moorings configuration:");
    if(!PartBegin)Moorings->Config(FtCount,FtObjs,ForcePoints);
    else RunException(met,"Simulation restart not allowed when moorings are used.");
    Moorings->VisuConfig(""," ");
  }  //<vs_moordyyn_end>

  //-Prepares ForcePoints configuration.  //<vs_moordyyn_ini>
  if(ForcePoints){
    Log->Printf("ForcePoints configuration:");
    ForcePoints->Config(FtCount,FtObjs,PeriActive,PeriX,PeriY,PeriZ,PeriXinc,PeriYinc,PeriZinc);
    if(!PartBegin)ForcePoints->CheckPoints(MkInfo,np,idp,pos);
    else RunException(met,"Simulation restart not allowed when FtForces is used.");
    ForcePoints->VisuConfig(""," ",FtCount,FtObjs);
  }  //<vs_moordyyn_end>

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
  if(xml.GetNodeSimple("case.execution.special.savedt",true)){
    SaveDt=new JSaveDt(Log);
    SaveDt->Config(&xml,"case.execution.special.savedt",TimeMax,TimePart);
    SaveDt->VisuConfig("SaveDt configuration:"," ");
  }

  //-Prepares BoundCorr configuration.  //<vs_innlet_ini>
  if(BoundCorr){
    Log->Print("BoundCorr configuration:");
    if(PartBegin)RunException(met,"Simulation restart not allowed when BoundCorr is used.");
    BoundCorr->RunAutoConfig(PartsInit);
    BoundCorr->VisuConfig(""," ");
  }//<vs_innlet_end>

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
/// Calculates predefined movement of boundary particles.
/// Calcula movimiento predefinido de boundary particles.
//==============================================================================
bool JSph::CalcMotion(double stepdt){
  const bool motsim=true;
  const JSphMotion::TpMotionMode mode=(motsim? JSphMotion::MOMT_Simple: JSphMotion::MOMT_Ace2dt);
  SphMotion->ProcesTime(mode,TimeStep,stepdt);
  const bool active=SphMotion->GetActiveMotion();
  if(ChronoObjects && ChronoObjects->GetWithMotion() && active){ //<vs_chroono_ini> 
    const unsigned nref=SphMotion->GetNumObjects();
    for(unsigned ref=0;ref<nref;ref++){
      const StMotionData& m=SphMotion->GetMotionData(ref);
      if(m.type!=MOTT_None)ChronoObjects->SetMovingData(m.mkbound,m.type==MOTT_Linear,m.linmov,m.matmov,stepdt);
    }
  } //<vs_chroono_end>
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
      if(motsim)SphMotion->SetMotionData   (WaveGen->GetMotion   (svdata,c,TimeStep,stepdt));
      else      SphMotion->SetMotionDataAce(WaveGen->GetMotionAce(svdata,c,TimeStep,stepdt));
    }
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
  parthead.ConfigSymmetry(Symmetry); //<vs_syymmetry>
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
      const StFloatingData &ft=FtObjs[cf];
      DataFloatBi4->AddHeadData(cf,ft.mkbound,ft.begin,ft.count,ft.mass,ft.massp,ft.radius);
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
  if(UseChrono && PeriActive!=0)log->Print("*** Maybe some Chrono object went beyond the periodic limits. Be careful when combining the use of Chrono with periodic limits."); //<vs_chroono>
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
  const char met[]="SavePartData";
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
      if(!(err=arrays.CheckErrorArray("Pos" ,TypeDouble3,npok)).empty())RunException(met,err);
      if(!(err=arrays.CheckErrorArray("Idp" ,TypeUint   ,npok)).empty())RunException(met,err);
      if(!(err=arrays.CheckErrorArray("Vel" ,TypeFloat3 ,npok)).empty())RunException(met,err);
      if(!(err=arrays.CheckErrorArray("Rhop",TypeFloat  ,npok)).empty())RunException(met,err);
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
    if(!(err=arrays2.CheckErrorArray("Pos" ,TypeDouble3,npok)).empty())RunException(met,err);
    if(!(err=arrays2.CheckErrorArray("Idp" ,TypeUint   ,npok)).empty())RunException(met,err);
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
void JSph::SaveData(unsigned npok,const JDataArrays& arrays
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
  
  //-Shows info of the new inlet particles.  //<vs_innlet_ini> 
  bool printnp=true;
  if(InOut && InOut->GetNewNpPart()){
    Log->Printf("  Particles new: %u (total new: %llu)  -  Current np: %u",InOut->GetNewNpPart(),InOut->GetNewNpTotal(),npok);
    InOut->ClearNewNpPart();
    printnp=false;
  }  //<vs_innlet_end> 
  
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
  if(ChronoObjects)ChronoObjects->SavePart(Part); //<vs_chroono>
  if(Moorings)Moorings->SaveData(Part);        //<vs_moordyyn>
  if(ForcePoints)ForcePoints->SaveData(Part);  //<vs_moordyyn>
  if(BoundCorr && BoundCorr->GetUseMotion())BoundCorr->SaveData(Part);  //<vs_innlet>
}

//==============================================================================
/// Generates VTK file with domain of the particles.
/// Genera fichero VTK con el dominio de las particulas.
//==============================================================================
void JSph::SaveDomainVtk(unsigned ndom,const tdouble3 *vdom)const{ 
  if(vdom){
    string fname=fun::FileNameSec("Domain.vtk",Part);
    JVtkLib::SaveVtkBoxes(DirDataOut+fname,ndom,vdom,H*0.5f);
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
  tdouble3 pmin=MapRealPosMin;
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

//<vs_mddbc_ini>
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
} //<vs_mddbc_end>

//==============================================================================
/// Adds basic information of resume to hinfo & dinfo.
/// Anhade la informacion basica de resumen a hinfo y dinfo.
//==============================================================================
void JSph::GetResInfo(float tsim,float ttot,const std::string &headplus,const std::string &detplus,std::string &hinfo,std::string &dinfo){
  hinfo=hinfo+"#RunName;RunCode;DateTime;Np;TSimul;TSeg;TTotal;MemCpu;MemGpu;Steps;PartFiles;PartsOut;MaxParticles;MaxCells;Hw;StepAlgo;Kernel;Viscosity;ViscoValue;DensityCorrection;TMax;Nbound;Nfixed;H;RhopOut;PartsRhopOut;PartsVelOut;CellMode"+headplus;
  dinfo=dinfo+ RunName+ ";"+ RunCode+ ";"+ RunTimeDate+ ";"+ fun::UintStr(CaseNp);
  dinfo=dinfo+ ";"+ fun::FloatStr(tsim)+ ";"+ fun::FloatStr(tsim/float(TimeStep))+ ";"+ fun::FloatStr(ttot);
  dinfo=dinfo+ ";"+ fun::LongStr(MaxNumbers.memcpu)+ ";"+ fun::LongStr(MaxNumbers.memgpu);
  const unsigned nout=GetOutPosCount()+GetOutRhopCount()+GetOutMoveCount();
  dinfo=dinfo+ ";"+ fun::IntStr(Nstep)+ ";"+ fun::IntStr(Part)+ ";"+ fun::UintStr(nout);
  dinfo=dinfo+ ";"+ fun::UintStr(MaxNumbers.particles)+ ";"+ fun::UintStr(MaxNumbers.cells);
  dinfo=dinfo+ ";"+ Hardware+ ";"+ GetStepName(TStep)+ ";"+ GetKernelName(TKernel)+ ";"+ GetViscoName(TVisco)+ ";"+ fun::FloatStr(Visco);
  dinfo=dinfo+ ";"+ GetDDTConfig()+ ";"+ fun::FloatStr(float(TimeMax));
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
/// Returns the name of the kernel in text format.
/// Devuelve el nombre del kernel en texto.
//==============================================================================
std::string JSph::GetKernelName(TpKernel tkernel){
  string tx;
  if(tkernel==KERNEL_Cubic)tx="Cubic";
  else if(tkernel==KERNEL_Wendland)tx="Wendland";
  else if(tkernel==KERNEL_Gaussian)tx="Gaussian";
  else if(tkernel==KERNEL_WendlandC6)tx="WendlandC6"; //<vs_praticalsskq>
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
  else if(tboundary==BC_MDBC)tx="mDBC"; //<vs_mddbc>
  else tx="???";
  return(tx);
}

//<vs_mddbc_ini>
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
} //<vs_mddbc_end>

//==============================================================================
/// Returns name of Density Diffusion Term in text format.
/// Devuelve nombre del Density Diffusion Term en texto.
//==============================================================================
std::string JSph::GetDDTName(TpDensity tdensity){
  string tx;
  if(tdensity==DDT_None)tx="None";
  else if(tdensity==DDT_DDT)tx="Molteni and Colagrossi 2009";
  else if(tdensity==DDT_DDT2)tx="Fourtakas et al 2019 (inner)";     //<vs_dtt2>
  else if(tdensity==DDT_DDT2Full)tx="Fourtakas et al 2019 (full)";  //<vs_dtt2>
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
    {//-For InOut particles.  //<vs_innlet_ini> 
      typecode tv=CODE_GetTypeAndValue(code[p+pini]);
      if(tv>=CODE_TYPE_FLUID_INOUT)xkind[p]=byte(tv-CODE_TYPE_FLUID_INOUT+10);
    }  //<vs_innlet_end> 
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


