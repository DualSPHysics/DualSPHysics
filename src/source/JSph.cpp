//HEAD_DSPH
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

/// \file JSph.cpp \brief Implements the class \ref JSph.

#include "JSph.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsMath.h"
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
#include "JDsDcell.h"
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
#include "JDsPartMotionSave.h"
#include "JDsPartFloatSave.h"
#include "JPartFloatInfoBi4.h"
#include "JDsExtraData.h"
#include "JDsPartsOut.h"
#include "JSphShifting.h"
#include "JSphShiftingAdv.h" //<vs_advshift>
#include "JDsDamping.h"
#include "JDsInitialize.h"
#include "JSphInOut.h"
#include "JSphFlexStruc.h"  //<vs_flexstruc>
#include "JDsOutputParts.h" //<vs_outpaarts>
#include "JDsPips.h"
#include "JLinearValue.h"
#include "JDataArrays.h"
#include "JOutputCsv.h"
#include "JSpVtkData.h"
#include "JSpVtkShape.h"
#include <algorithm>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JSph::JSph(int gpucount,bool withmpi):Cpu(gpucount==0)
  ,Mgpu(gpucount>1),GpuCount(gpucount),WithMpi(withmpi)
{
  ClassName="JSph";
  DgNum=0;
  DataBi4=NULL;
  DataOutBi4=NULL;
  PartMotionSave=NULL;
  PartFloatSave=NULL;
  SvExtraDataBi4=NULL;
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
  ShiftingAdv=NULL;   //<vs_advshift>
  Damping=NULL;
  AccInput=NULL;
  PartsLoaded=NULL;
  InOut=NULL;
  FlexStruc=NULL;   //<vs_flexstruc>
  OutputParts=NULL; //<vs_outpaarts>
  DsPips=NULL;
  //-Auxiliary variables for floating bodies calculations.
  Fto_ExForceLinAng=NULL;
  Fto_AceLinAng=NULL;
  Fto_VelLinAng=NULL;
  Fto_Center=NULL;

  InitVars();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSph::~JSph(){
  DestructorActive=true;
  delete DataBi4;        DataBi4=NULL;
  delete DataOutBi4;     DataOutBi4=NULL;
  delete PartMotionSave; PartMotionSave=NULL;
  delete PartFloatSave;  PartFloatSave=NULL;
  delete SvExtraDataBi4; SvExtraDataBi4=NULL;
  delete PartsOut;       PartsOut=NULL;
  delete ViscoTime;      ViscoTime=NULL;
  delete FixedDt;        FixedDt=NULL;
  delete SaveDt;         SaveDt=NULL;
  delete OutputTime;     OutputTime=NULL;
  delete MkInfo;         MkInfo=NULL;
  delete PartsInit;      PartsInit=NULL;
  delete DsMotion;       DsMotion=NULL;
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
  delete ShiftingAdv;   ShiftingAdv=NULL; //<vs_advshift>
  delete Damping;       Damping=NULL;
  delete AccInput;      AccInput=NULL; 
  delete PartsLoaded;   PartsLoaded=NULL;
  delete InOut;         InOut=NULL;
  delete FlexStruc;     FlexStruc=NULL;   //<vs_flexstruc>
  delete DsPips;        DsPips=NULL;
  delete OutputParts;   OutputParts=NULL; //<vs_outpaarts>
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSph::InitVars(){
  ClearCfgDomain();
  OutPosCount=OutRhoCount=OutMovCount=0;
  Simulate2D=false;
  Simulate2DPosY=0;
  Stable=false;
  SvPosDouble=false;
  RunCode=CalcRunCode();
  RunTimeDate="";
  CaseName=""; DirCase=""; RunName="";
  DirOut="";
  DirDataOut="";
  DirVtkOut="";
  FileXml="";
  TStep=STEP_None;
  InterStep=INTERSTEP_None;
  InterNum=0;
  VerletSteps=40;
  TKernel=KERNEL_Wendland;
  KCubic ={0,0,0,0,0,0,0,0};
  KWend  ={0,0};
  TVisco=VISCO_None;
  TDensity=DDT_None; DDTValue=0; DDTArray=false;
  DDTRamp=TDouble3(0); //<vs_ddramp>
  ShiftingMode=(Shifting? Shifting->GetShiftMode(): SHIFT_None);
  Visco=0; ViscoBoundFactor=1;
  TBoundary=BC_DBC;
  SlipMode=SLIP_None;
  MdbcCorrector=false;
  UseNormals=false;
  UseNormalsFt=false;
  SvNormals=false;
  AbortNoNormals=true;
  NoPenetration=false;                //<vs_m2dbcNP>
  TMdbc2=MDBC2_None;                  //<vs_m2dbcNP>  
  UseDEM=false;  //(DEM)
  delete[] DemData; DemData=NULL;  //(DEM)
  UseChrono=false;
  UseWavegenAWAS=false;
  RhopOut=true; RhopOutMin=700; RhopOutMax=1300;
  TimeMax=TimePart=0;
  TimePartExtra=-1;
  TimePartExtraNext=DBL_MAX;
  NstepsBreak=0;
  SvAllSteps=false;
  NoRtimes=false;
  TerminateMt=0;
  TerminateTimeMax=DBL_MAX;
  DtIni=0;
  DtMin=0;
  CoefDtMin=0;
  DtAllParticles=false;
  MinFluidStop=0;
  NpfMinimum=0;
  PartsOutWrn=1; PartsOutTotWrn=10;

  SvData=byte(SDAT_Binx)|byte(SDAT_Info);
  SvExtraParts="";
  SvRes=false;
  SvTimers=false;
  SvDomainVtk=false;

  KernelH=CteB=Gamma=RhopZero=0;
  CFLnumber=0;
  Dp=0;
  MassFluid=MassBound=0;
  Gravity=TFloat3(0);

  KernelSize=KernelSize2=0;
  Cs0=0;
  Eta2=0;
  SpsSmag=SpsBlin=0;
  DDTkhCte=DDTkh=DDTgz=0;
  memset(&CSP,0,sizeof(StCteSph));

  CasePosMin=CasePosMax=TDouble3(0);
  CaseNp=CaseNbound=CaseNfixed=CaseNmoving=CaseNfloat=CaseNfluid=CaseNpb=0;
  CaseNflexstruc=0; //<vs_flexstruc>

  PeriActive=0; PeriX=PeriY=PeriZ=false;
  PeriXinc=PeriYinc=PeriZinc=TDouble3(0);

  PartBeginDir=""; 
  PartBegin=PartBeginFirst=0;
  PartBeginTimeStep=0; 
  PartBeginTotalNp=0;
  RestartChrono=false;

  WrnPartsOut=true;

  FtCount=0;
  FtPause=0;
  RigidMode=FTRIGID_Sph;
  FtMode=FTMODE_Sph;
  FtConstraints=false;
  FtIgnoreRadius=false;
  WithFloating=false;

  FlexStrucCount=0; //<vs_flexstruc>
  FlexStrucCs0=0;   //<vs_flexstruc>

  DivCleaning=false;  //<vs_divclean>
  DivCleanKp=0;       //<vs_divclean>

  AllocMemoryFloating(0,false);

  CellMode=CELLMODE_None;
  CellDomFixed=false;
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
  PartDtModif=0;

  PartIni=Part=0; 
  Nstep=0; PartNstep=0;
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
/// Saves the linear and angular acceleration of each floating object in 
/// a csv file (for debug).
//==============================================================================
void JSph::SaveFtAceFun(double dt,bool predictor,const StFloatingData* ftobjs){
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
        scsv << "nstep;time [s];dt [s];predictor;acelin.x [m/s^2];acelin.y [m/s^2];acelin.z [m/s^2]";
        scsv << "aceang.x [rad/s^2];aceang.y [rad/s^2];aceang.z [rad/s^2]" << jcsv::Endl();
      }
      //-Saves data.
      scsv.SetData();
      scsv << nstep << timestep << dt << (predictor?"True":"False");
      scsv << ftobjs[cf].preacelin;
      scsv << ftobjs[cf].preaceang;
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
void JSph::ConfigDomainResize(std::string key,const JCaseEParms* eparms){
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
  //-Free allocated memory for floatings.
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
  //-Auxiliary variables for floating bodies calculations.
  delete[] Fto_ExForceLinAng; Fto_ExForceLinAng=NULL;
  delete[] Fto_AceLinAng;     Fto_AceLinAng=NULL;
  delete[] Fto_VelLinAng;     Fto_VelLinAng=NULL;
  delete[] Fto_Center;        Fto_Center=NULL;
  
  //-Allocate memory for floatings.
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
    //-Auxiliary variables for floating bodies calculations.
    Fto_ExForceLinAng =new tfloat6 [ftcount];
    Fto_AceLinAng     =new tfloat6 [ftcount];
    Fto_VelLinAng     =new tfloat6 [ftcount];
    Fto_Center        =new tdouble3[ftcount];
  }
}

//==============================================================================
/// Returns the allocated memory in CPU.
//==============================================================================
llong JSph::GetAllocMemoryCpu()const{  
  //-Allocated in AllocMemoryCase().
  llong s=0;
  //-Allocated in AllocMemoryFloating().
  if(FtCount){
    s+=(sizeof(StFloatingData)+sizeof(tfloat6)*3+sizeof(tdouble3))*FtCount;
  }
  //-Allocated in other objects.
  if(PartsOut)s+=PartsOut->GetAllocMemory();
  if(ViscoTime)s+=ViscoTime->GetAllocMemory();
  if(FixedDt)s+=FixedDt->GetAllocMemory();
  if(AccInput)s+=AccInput->GetAllocMemory();
  if(PartsLoaded)s+=PartsLoaded->GetAllocMemory();
  return(s);
}

//==============================================================================
/// Returns list of available features in current compilation.
//==============================================================================
std::string JSph::GetFeatureList(){
  string list;
  if(AVAILABLE_CHRONO     )list=list+", Project Chrono coupling";
  if(AVAILABLE_MOORDYNPLUS)list=list+", MoorDynPlus coupling";
  if(AVAILABLE_WAVEGEN    )list=list+", Wave generation";
  if(AVAILABLE_FLEXSTRUC  )list=list+", FlexStruc";    //<vs_flexstruc>
                           list=list+", mDBC no-slip"; //<vs_m2dbc>
  #ifdef CODE_SIZE4
                           list=list+", MkWord";
  #endif
  if(AVAILABLE_GPU        )list=list+", GPU execution";
  #ifdef OMP_USE
                           list=list+", OpenMP execution";
  #endif
  if(list.size()>2)list=list.substr(2);
  return(list);
}

//==============================================================================
/// Loads the configuration of the execution.
//==============================================================================
void JSph::LoadConfig(const JSphCfgRun* cfg){
  TimerTot.Start();

  //-Loads basic configuration from execution parameters.
  //------------------------------------------------------
  Stable=cfg->Stable;
  SvPosDouble=false; //-Options by default.
  SvExtraParts="";   //-Options by default.
  DirOut=fun::GetDirWithSlash(cfg->DirOut);
  DirDataOut=(!cfg->DirDataOut.empty()? fun::GetDirWithSlash(DirOut+cfg->DirDataOut): DirOut);
  DirVtkOut=fun::GetDirWithSlash(DirOut+"vtks");
  CaseName=cfg->CaseName; 
  DirCase=fun::GetDirWithSlash(fun::GetDirParent(CaseName));
  CaseName=CaseName.substr(DirCase.length());
  if(CaseName.empty())Run_Exceptioon(
    "Name of the case for execution was not indicated.");
  RunName=(cfg->RunName.length()? cfg->RunName: CaseName);
  FileXml=DirCase+CaseName+".xml";
  PartBeginDir=cfg->PartBeginDir; 
  PartBegin=cfg->PartBegin; 
  PartBeginFirst=cfg->PartBeginFirst;
  RestartChrono=cfg->RestartChrono;

  //-Output options:
  CsvSepComa=cfg->CsvSepComa;
  SvData=byte(SDAT_None); 
  if(cfg->Sv_Csv&&!WithMpi)SvData|=byte(SDAT_Csv);
  if(cfg->Sv_Binx)SvData|=byte(SDAT_Binx);
  if(cfg->Sv_Info)SvData|=byte(SDAT_Info);
  if(cfg->Sv_Vtk)SvData|=byte(SDAT_Vtk);
  SvNormals=cfg->SvNormals;
  SvRes=cfg->SvRes;
  SvTimers=cfg->SvTimers;
  SvDomainVtk=cfg->SvDomainVtk;

  printf("\n");
  RunTimeDate=fun::GetDateTime();
  Log->Printf("[Initialising %s  %s]",ClassName.c_str(),RunTimeDate.c_str());
  const string runpath=AppInfo.GetRunPath();
  Log->Printf("ProgramFile=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,AppInfo.GetRunCommand()),3).c_str());
  Log->Printf("Available features: %s.",GetFeatureList().c_str());
  Log->Printf("ExecutionDir=\"%s\"",fun::GetPathLevels(runpath,3).c_str());
  Log->Printf("XmlFile=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,FileXml),3).c_str());
  Log->Printf("OutputDir=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,DirOut),3).c_str());
  Log->Printf("OutputDataDir=\"%s\"",fun::GetPathLevels(fun::GetCanonicalPath(runpath,DirDataOut),3).c_str());
  
  //-Creates Output directories when it is necessary.
  if(AppInfo.GetCreateDirs()){
    fun::MkdirPath(DirOut);
    if(DirOut!=DirDataOut)fun::MkdirPath(DirDataOut);
    if(SvData&SDAT_Vtk || SvNormals)fun::MkdirPath(DirVtkOut);
  }

  if(PartBegin){
    Log->Print(fun::VarStr("PartBegin",PartBegin));
    Log->Print(fun::VarStr("PartBeginDir",PartBeginDir));
    Log->Print(fun::VarStr("PartBeginFirst",PartBeginFirst));
  }

  //-Loads case configuration from XML and command line.
  LoadCaseConfig(cfg);

  //-PIPS configuration.
  if(cfg->PipsMode){
    const bool svdata=(cfg->PipsMode==2);
    const int ntimes=(TStep==STEP_Symplectic? 2: 1);
    DsPips=new JDsPips(Cpu,cfg->PipsSteps,svdata,ntimes);
  }
}

//==============================================================================
/// Loads kernel selection to compute kernel values.
//==============================================================================
void JSph::LoadKernelSelection(const JSphCfgRun* cfg,const JXml* cxml){
  //-Load kernel selection from execution parameters from XML.
  JCaseEParms eparms;
  eparms.LoadXml(cxml,"case.execution.parameters");
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
void JSph::LoadConfigCtes(const JXml* cxml){
  JCaseCtes ctes;
  ctes.LoadXmlRun(cxml,"case.execution.constants");

  Simulate2D=ctes.GetData2D();
  Simulate2DPosY=ctes.GetData2DPosY();
  KernelH=(float)ctes.GetH();
  CteB=(float)ctes.GetB();
  Gamma=(float)ctes.GetGamma();
  RhopZero=(float)ctes.GetRhop0();
  CFLnumber=ctes.GetCFLnumber();
  Dp=ctes.GetDp();
  MassFluid=(float)ctes.GetMassFluid();
  MassBound=(float)ctes.GetMassBound();
  Gravity=ToTFloat3(ctes.GetGravity());
  if(ctes.GetEps()!=0)Log->PrintWarning("Eps value is not used (this correction was removed).");
}

//==============================================================================
/// Loads execution parameters from XML.
//==============================================================================
void JSph::LoadConfigParameters(const JXml* cxml){
  JCaseEParms eparms;
  eparms.LoadXml(cxml,"case.execution.parameters");
  if(eparms.Exists("FtSaveAce"))SaveFtAce=(eparms.GetValueInt("FtSaveAce",true,0)!=0); //-For Debug.
  if(eparms.Exists("PosDouble")){
    Log->PrintWarning("The parameter \'PosDouble\' is deprecated.");
    SvPosDouble=(eparms.GetValueInt("PosDouble")==2);
  }
  if(eparms.Exists("SavePosDouble"))SvPosDouble=(eparms.GetValueInt("SavePosDouble",true,0)!=0);
  if(eparms.Exists("SaveExtraParts"))SvExtraParts=eparms.GetValue("SaveExtraParts");
  switch(eparms.GetValueInt("RigidAlgorithm",true,1)){ //(DEM)
    case 0:  RigidMode=FTRIGID_Free;    FtMode=FTMODE_Ext;                   break;
    case 1:  RigidMode=FTRIGID_Sph;     FtMode=FTMODE_Sph;                   break;
    case 2:  RigidMode=FTRIGID_Dem;     FtMode=FTMODE_Ext;  UseDEM=true;     break;
    case 3:  RigidMode=FTRIGID_Chrono;  FtMode=FTMODE_Ext;  UseChrono=true;  break;
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
    case 3:  TVisco=VISCO_Laminar;     break;
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
    //<vs_m2dbcNP_ini>
    if(SlipMode>=SLIP_NoSlip) NoPenetration=eparms.GetValueBool("NoPenetration",true,false);   
    if(SlipMode>=SLIP_NoSlip) TMdbc2=(NoPenetration? MDBC2_NoPen: MDBC2_Std);
    //<vs_m2dbcNP_end>
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
    switch(eparms.GetValueInt("Shifting",true,0)){
      case 0:  shiftmode=SHIFT_None;     break;
      case 1:  shiftmode=SHIFT_NoBound;  break;
      case 2:  shiftmode=SHIFT_NoFixed;  break;
      case 3:  shiftmode=SHIFT_Full;     break;
      case 4:  shiftmode=SHIFT_FS;       break; //<vs_advshift>
      default: Run_ExceptioonFile("Shifting mode in <execution><parameters> is not valid.",FileXml);
    }
    if(shiftmode<=SHIFT_Full){
      const float shiftcoef=eparms.GetValueFloat("ShiftCoef",true,-2);
      const float shifttfs=eparms.GetValueFloat("ShiftTFS",true,0);
      Shifting=new JSphShifting(Simulate2D,Dp,KernelH);
      Shifting->ConfigBasic(shiftmode,shiftcoef,shifttfs);
    }
    if(shiftmode==SHIFT_FS){ //<vs_advshift_ini>
      const float shiftcoef=eparms.GetValueFloat("ShiftAdvCoef",true,-0.01f);
      const bool  aleform=eparms.GetValueBool("ShiftAdvALE",true,false);
      const bool  ncpress=eparms.GetValueBool("ShiftAdvNCPress",true,false);
      ShiftingAdv=new JSphShiftingAdv(Simulate2D,Dp,KernelH);
      ShiftingAdv->ConfigBasic(shiftcoef,aleform,ncpress);
    } //<vs_advshift_end>
  }

  //<vs_divclean_ini>
  //Divergence cleaning configuration.
  if(eparms.Exists("DivCleanKp")){
    #ifndef AVAILABLE_DIVCLEAN
      Run_Exceptioon("Divergence cleaning is not available in the current compilation.");
    #else
      DivCleanKp=(eparms.GetValueFloat("DivCleanKp",true,0));
      DivCleaning=DivCleanKp>0;
      if(DivCleaning && TStep==STEP_Verlet)Run_Exceptioon("Divergence cleaning is not available with Verlet time integrator.");
    #endif
  }
  //<vs_divclean_ini>

  WrnPartsOut=(eparms.GetValueInt("WrnPartsOut",true,1)!=0);
  FtPause=eparms.GetValueFloat("FtPause",true,0);
  FtIgnoreRadius=(eparms.GetValueInt("FtIgnoreRadius",true,0)!=0);

  TimeMax=eparms.GetValueDouble("TimeMax");
  TimePart=eparms.GetValueDouble("TimeOut");

  TimePartExtra=eparms.GetValueDouble("TimeOutExtra",true,-1);
  if(TimePartExtra>=0)TimePartExtraNext=0;

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

  if(eparms.Exists("PartsOutMax"))Log->PrintWarning("The XML option \'PartsOutMax\' is deprecated. Use \'MinFluidStop\' option.");
  MinFluidStop=eparms.GetValueFloat("MinFluidStop",true,0);

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
/// Computes next time for extra output data.
//==============================================================================
void JSph::TimeOutExtraUpdate(double timestep){
  if(TimePartExtra>=0){
    double tnext=timestep;
    if(TimePartExtra){
      unsigned ct=unsigned(timestep/TimePartExtra);
      tnext=TimePartExtra*ct;
      for(;tnext<=timestep;ct++)tnext=TimePartExtra*ct;
    }
    TimePartExtraNext=tnext;
  }
}

//==============================================================================
/// Loads the case configuration to be executed.
//==============================================================================
void JSph::LoadConfigCommands(const JSphCfgRun* cfg){
  //-Aplies configuration using command line.
  if(cfg->SvPosDouble>=0)SvPosDouble=(cfg->SvPosDouble!=0);
  if(cfg->SvExtraParts!="undefined")SvExtraParts=cfg->SvExtraParts;
  if(cfg->TBoundary>=0){
    TBoundary=BC_DBC;
    SlipMode=SLIP_None;
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
  }
 // if(TBoundary==BC_MDBC){
  //  if(SlipMode!=SLIP_Vel0 && SlipMode!=SLIP_NoSlip && SlipMode != SLIP_FreeSlip)Run_Exceptioon(
  //    "Only the slip modes velocity=0 and no-slip are allowed with mDBC conditions.");
  //}
  MdbcCorrector=(SlipMode>=SLIP_NoSlip);
  UseNormals=(TBoundary==BC_MDBC);
  if(cfg->TBoundary>=0){
    NoPenetration=false;
    TMdbc2=MDBC2_None;
    if(SlipMode>=SLIP_NoSlip)NoPenetration=cfg->NoPenetration;
    if(SlipMode>=SLIP_NoSlip) TMdbc2=(NoPenetration? MDBC2_NoPen: MDBC2_Std);
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
  if(TDensity!=DDT_None && cfg->DDTValueTRamp>0 && cfg->DDTValueMax>0){ //<vs_ddramp_ini>
    DDTRamp=TDouble3(cfg->DDTValueTRamp,cfg->DDTValueTMax,cfg->DDTValueMax);
  } //<vs_ddramp_end>
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
      case 4:  shiftmode=SHIFT_FS;       break;
      default: Run_Exceptioon("Shifting mode is not valid.");
    }
    if(shiftmode<=shiftmode){
      delete ShiftingAdv; ShiftingAdv=NULL; //<vs_advshift>
      if(!Shifting)Shifting=new JSphShifting(Simulate2D,Dp,KernelH);
      Shifting->ConfigBasic(shiftmode);
    }
    if(shiftmode==SHIFT_FS){ //<vs_advshift_ini>
      delete Shifting; Shifting=NULL;
      if(!ShiftingAdv)ShiftingAdv=new JSphShiftingAdv(Simulate2D,Dp,KernelH);
      ShiftingAdv->ConfigBasic(ShiftingAdv->GetShiftCoef(),cfg->ShiftAdvALE,cfg->ShiftAdvNCP);
    } //<vs_advshift_end>
  }

  if(cfg->CFLnumber>0)CFLnumber=cfg->CFLnumber;
  if(cfg->FtPause>=0)FtPause=cfg->FtPause;
  if(cfg->TimeMax>0)TimeMax=cfg->TimeMax;
  NstepsBreak=cfg->NstepsBreak;
  if(NstepsBreak)Log->PrintfWarning("The execution will be cancelled after %d simulation steps.",NstepsBreak);
  SvAllSteps=cfg->SvAllSteps;
  NoRtimes=cfg->NoRtimes;
  //-Configuration of JDsOutputTime with TimePart.
  OutputTime=new JDsOutputTime();
  if(cfg->TimePart>=0){
    TimePart=cfg->TimePart;
    OutputTime->Config(TimePart);
  }
  else OutputTime->Config(FileXml,"case.execution.special.timeout",TimePart);
  //-Configuration of extra output time.
  if(cfg->TimePartExtra!=DBL_MAX){
    TimePartExtra=-1;
    TimePartExtraNext=DBL_MAX;
    if(cfg->TimePartExtra>=0){
      TimePartExtra=cfg->TimePartExtra;
      TimePartExtraNext=0;
    }
  }

  //-Configuration of domain limits.
  CellMode=cfg->CellMode;
  CellDomFixed=cfg->CellDomFixed;
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
/// Loads the case configuration to be executed.
//==============================================================================
void JSph::LoadCaseConfig(const JSphCfgRun* cfg){
  if(!fun::FileExists(FileXml))
    Run_ExceptioonFile("Case configuration was not found.",FileXml);
  JXml xml;
  xml.LoadFile(FileXml);
  const JXml* cxml=&xml;
  //-Shows pre-processing application generating the XML file.
  Log->Printf("XML-App: %s",cxml->GetAttributeStr(cxml->GetNodeError("case")->ToElement()
    ,"app",true,"unknown").c_str());

  //-Loads kernel selection to compute kernel values.
  LoadKernelSelection(cfg,cxml);
  //-Loads predefined constants.
  LoadConfigCtes(cxml);
  //-Configures value of calculated constants and loads CSP structure.
  ConfigConstants1(Simulate2D);

  //-Execution parameters from XML.
  LoadConfigParameters(cxml);
  //-Execution parameters from commands.
  LoadConfigCommands(cfg);
  //-Configures other constants according to formulation options and loads more values in CSP structure.
  ConfigConstants2();

  //-Particle data.
  JCaseParts parts;
  parts.LoadXml(cxml,"case.execution.particles");
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
  GaugeSystem=new JGaugeSystem(GpuCount);

  //-Configuration of AccInput.
  if(xml.GetNodeSimple("case.execution.special.accinputs",true)){
    AccInput=new JDsAccInput(DirCase,&xml,"case.execution.special.accinputs");
  }

  //-Configuration of ChronoObjects.
  if(xml.GetNodeSimple("case.execution.special.chrono",true)){
    if(!JChronoObjects::Available())Run_Exceptioon("DSPHChronoLib to use Chrono is not included in the current compilation.");
    ChronoObjects=new JChronoObjects(DirCase,CaseName,&xml,"case.execution.special.chrono",Dp,parts.GetMkBoundFirst(),Gravity,Simulate2D,FtPause);
    UseChrono=true;
    if(RigidMode==FTRIGID_Chrono && !ChronoObjects->GetUseCollision())
      Log->PrintfWarning("Collisions are configured in terms of Chrono (RigidAlgorithm=3) but the collisions are disabled in Chrono configuration.");
  }
  else if(RigidMode==FTRIGID_Chrono)Run_ExceptioonFile("Chrono configuration in XML file for defined rigid-algorithm is missing.",FileXml);

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
      WaveGen->ConfigPaddleParts(m.mkbound,ref,m.idbegin,m.count);
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
    RelaxZones=new JRelaxZones(useomp,!Cpu,Log,AppName,DirCase,CaseNfloat>0,CaseNbound,ToTDouble3(Gravity));
    RelaxZones->LoadXml(&xml,"case.execution.special.relaxationzones");
  }

  //-Configuration of Shifting with zones.
  if(xml.GetNodeSimple("case.execution.special.shifting",true)){
    if(ShiftingAdv)Run_ExceptioonFile("Advanced shifting is not compatible with standard shifting configuration in <special><shifting>.",FileXml); //<vs_advshift>
    if(Shifting)Run_ExceptioonFile("Shifting is defined several times (in <special><shifting> and <execution><parameters>).",FileXml);
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
      const JCasePartBlock& block=parts.GetBlock(c);
      if(block.Type==TpPartFloating){
        const JCasePartBlock_Floating& fblock=(const JCasePartBlock_Floating&)block;
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
        //-Additional data.
        fobj->extforcelin=TFloat3(0);
        fobj->extforceang=TFloat3(0);
        fobj->fluforcelin=TFloat3(0);
        fobj->fluforceang=TFloat3(0);
        fobj->preacelin=TFloat3(0);
        fobj->preaceang=TFloat3(0);
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
      const JCasePartBlock& block=parts.GetBlock(c);
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

  //<vs_flexstruc_ini>
  //-Configuration of flexible structures.
  if(xml.GetNodeSimple("case.execution.special.flexstrucs",true)){
    #ifndef AVAILABLE_FLEXSTRUC
      Run_Exceptioon("FlexStruc feature is not available in the current compilation.");
    #else
      if(!AVAILABLE_FLEXSTRUC)Run_Exceptioon("FlexStruc feature is not available in the current compilation.");
    #endif
    FlexStruc=new JSphFlexStruc(Simulate2D,Dp,&xml,"case.execution.special.flexstrucs",MkInfo);
    FlexStrucCount=FlexStruc->GetCount();
    FlexStrucCs0=FlexStruc->GetInitialSoundSpeed();
    if(!CaseNfluid||FlexStrucCs0>Cs0){
      DtIni=KernelH/FlexStrucCs0;
      DtMin=(KernelH/FlexStrucCs0)*CoefDtMin;
    }
  }
  //<vs_flexstruc_end>

  //-Configuration of Inlet/Outlet.
  if(xml.GetNodeSimple(JSphInOut::XmlPath,true)){
    InOut=new JSphInOut(Cpu,CSP,cxml,DirCase);
    NpDynamic=true;
    ReuseIds=InOut->GetReuseIds();
    if(ReuseIds)Run_Exceptioon("Inlet/Outlet with ReuseIds is not a valid option for now...");
  }
  
  //-Configuration of boundary extrapolated correction.
  if(xml.GetNodeSimple("case.execution.special.boundextrap",false))
    Run_Exceptioon("The XML section 'boundextrap' is obsolete. Use mDBC boundary conditions.");
  if(xml.GetNodeSimple("case.execution.special.boundcorr",false))
    Run_Exceptioon("The XML section 'boundcorr' is obsolete. Use mDBC boundary conditions.");
 
  //-Configuration of Moorings object.
  if(xml.GetNodeSimple("case.execution.special.moorings",true)){
    if(WithFloating){
      if(!AVAILABLE_MOORDYNPLUS)Run_Exceptioon(
        "Code for moorings and MoorDynPlus coupling is not included in the current compilation.");
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

  //-Configuration of OutputParts object. //<vs_outpaarts_ini>
  if(xml.GetNodeSimple("case.execution.special.outputparts",true)){
    OutputParts=new JDsOutputParts(Cpu,Dp*10);
    OutputParts->LoadXml(&xml,"case.execution.special.outputparts"
      ,MkInfo,FtCount,FtObjs);
  } //<vs_outpaarts_end>

  //-Checks invalid options for advanced shifting. //<vs_advshift_ini>
  if(ShiftingAdv){
    if(TStep==STEP_Verlet)Run_Exceptioon("Advanced shifting is not allowed with Verlet.");
  } //<vs_advshift_end>
  
  //<vs_flexstruc_ini>
  //-Checks invalid options for flexible structures.
  if(FlexStruc){
    if(PartBegin)   Run_Exceptioon("Simulation restart not allowed when flexible structures are used.");
    if(PeriActive)  Run_Exceptioon("Flexible structures are not allowed with periodic conditions.");
    if(TBoundary==BC_MDBC&&SlipMode==SLIP_Vel0)Run_Exceptioon("Flexible structures are not allowed with original mDBC (SlipMode=1).");
    if(TBoundary==BC_MDBC&&TStep==STEP_Verlet) Run_Exceptioon("Flexible structures with mDBC are not allowed with Verlet.");
  }
  //<vs_flexstruc_end>

  //-Defines NpfMinimum according to CaseNfluid. It is updated later to add initial inlet fluid particles.
  NpfMinimum=unsigned(MinFluidStop*CaseNfluid);

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
    const JCasePartBlock_Floating* fblock=(const JCasePartBlock_Floating* )block;
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
    const JSphMkBlock* pmk=MkInfo->Mkblock(c);
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
void JSph::LoadCodeParticles(unsigned np,const unsigned* idp,typecode* code)const{
  //-Assigns code to each group of particles.
  for(unsigned p=0;p<np;p++)code[p]=MkInfo->GetCodeById(idp[p]);
}

//==============================================================================
/// Config normals to point from boundary particle to ghost node (full distance).
//==============================================================================
void JSph::ConfigBoundNormals(unsigned np,unsigned npb,const tdouble3* pos
  ,const unsigned* idp,tfloat3* boundnor)
{
  //-Checks current normals.
  if(!PartBegin){
    bool ftnor=false;
    for(unsigned p=0;p<np && !ftnor;p++)
      ftnor=(idp[p]<CaseNbound && idp[p]>=CaseNpb && boundnor[p]!=TFloat3(0));
    UseNormalsFt=ftnor;
  }

  //-Loads normals for restart mode.
  if(PartBegin){
    if(JDsExtraDataLoad::ExistsPartData(PartBeginDir,int(PartBegin))){
      JDsExtraDataLoad edat(CaseNbound,CaseNfloat,Log);
      edat.LoadPartData(PartBeginDir,int(PartBegin));
      UseNormalsFt=edat.LoadNormals(np,npb,idp,boundnor);
    }
    else Run_Exceptioon(fun::PrintStr("No extra data available to restart at PART_%04d with mDBC.",PartBegin).c_str());
  }

  //-Saves normals from boundary particles to boundary limit.
  const string file1=DirOut+"CfgInit_Normals.vtk";
  Log->AddFileInfo(file1,"Saves VTK file with initial normals (from boundary particles to boundary limit).");
  SaveVtkNormals(file1,-1,np,npb,pos,idp,boundnor,(PartBegin? 0.5f: 1.f));
  //-Counts the null normals.
  unsigned nerr=0,nerrft=0;
  for(unsigned p=0;p<np;p++)if(idp[p]<CaseNbound){
    if(boundnor[p]==TFloat3(0)){
      if(idp[p]<CaseNpb)nerr++;
      else nerrft++;
    }
    if(!PartBegin)boundnor[p]=(boundnor[p]*2.f);
  }
  //-Saves normals from boundary particles to ghost node.
  const string file2=DirOut+"CfgInit_NormalsGhost.vtk";
  Log->AddFileInfo(file2,"Saves VTK file with initial normals (from boundary particles to ghost node).");
  SaveVtkNormals(file2,-1,np,npb,pos,idp,boundnor,1.f);
  if(nerr  >0)Log->PrintfWarning("There are %u of %u fixed or moving boundary particles without normal data.",nerr,npb);
  if(nerrft>0)Log->PrintfWarning("There are %u of %u floating particles without normal data.",nerrft,CaseNfloat);
  if(TBoundary==BC_MDBC && nerr==npb && nerrft==CaseNfloat){
    if(AbortNoNormals)Run_Exceptioon("No valid normal vectors for using mDBC.");
    else Log->PrintfWarning("No valid normal vectors for using mDBC.");
  }
  if(UseNormalsFt && (!UseChrono || !ChronoObjects->GetUseCollision()))
    Log->PrintWarning("When mDBC is applied to floating bodies, their collisions should be solved using Chrono (RigidAlgorithm=3).");
}

//==============================================================================
/// Sets DBL_MAX values by indicated values.
//==============================================================================
void JSph::PrepareCfgDomainValues(tdouble3& v,tdouble3 vdef)const{
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
    SpsSmag=float((2./3.)*pow((0.12*dp_sps),2));
    SpsBlin=float((2./3.)*0.0066*dp_sps*dp_sps); 
  }
  //-Constants for DDT.
  DDTkhCte=DDTkh=KernelSize*DDTValue;
  DDTgz=float(double(RhopZero)*double(fabs(Gravity.z))/double(CteB));
  //-Constants for Dt.
  if(!DtIni)DtIni=KernelH/Cs0;
  if(!DtMin)DtMin=(KernelH/Cs0)*CoefDtMin;

  //-Loads main SPH constants and configurations in CSP.
  CSP.spssmag       =SpsSmag;
  CSP.spsblin       =SpsBlin;
  CSP.ddtkhcte      =DDTkhCte;
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
  //-SavePosDouble. 
  Log->Print(fun::VarStr("SavePosDouble",SvPosDouble));
  if(SvPosDouble)ConfigInfo=ConfigInfo+sep+"SvPosDouble";
  //-SaveExtraData. 
  Log->Print(fun::VarStr("SvExtraParts",SvExtraParts));
  if(!SvExtraParts.empty())ConfigInfo=ConfigInfo+sep+"SvExtraParts";
  //-Other configurations. 
  Log->Print(fun::VarStr("SaveFtAce",SaveFtAce));
  Log->Print(fun::VarStr("SvTimers",SvTimers));
  if(DsPips)Log->Print(fun::VarStr("PIPS-steps",DsPips->StepsNum));
  //-Boundary. 
  Log->Print(fun::VarStr("Boundary",GetBoundName(TBoundary)));
  ConfigInfo=ConfigInfo+sep+GetBoundName(TBoundary);
  if(TBoundary==BC_MDBC){
    Log->Print(fun::VarStr("  SlipMode",GetSlipName(SlipMode)));
    Log->Print(fun::VarStr("  mDBC-Corrector",MdbcCorrector));
    Log->Print(fun::VarStr("  No Penetration",NoPenetration));
    ConfigInfo=ConfigInfo+"("+GetSlipName(SlipMode);
    if(MdbcCorrector)ConfigInfo=ConfigInfo+" - Corrector";
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
  if(ViscoBoundFactor!=1.f && SlipMode>=SLIP_NoSlip)
    Log->PrintWarning("ViscoBoundFactor should be 1.0 when mDBC no-slip or free slip is used.");
  if(ViscoTime)Log->Print(fun::VarStr("ViscoTime",ViscoTime->GetFile()));
  ConfigInfo=ConfigInfo+sep+"Visco_"+GetViscoName(TVisco)+fun::PrintStr("(%gb%g)",Visco,ViscoBoundFactor);
  //-DensityDiffusion.
  Log->Print(fun::VarStr("DensityDiffusion",GetDDTName(TDensity)));
  ConfigInfo=ConfigInfo+sep+fun::PrintStr("DDT%d",int(TDensity));
  if(TDensity!=DDT_None){
    Log->Print(fun::VarStr("  DensityDiffusionValue",DDTValue));
    //Log->Print(fun::VarStr("DensityDiffusionArray",DDTArray));
    string cinfo=fun::PrintStr("(%g)",DDTValue);
    if(DDTRamp.x){ //<vs_ddramp_ini>
      Log->Print(fun::VarStr("  DensityDiffusionTRamp"   ,DDTRamp.x));
      Log->Print(fun::VarStr("  DensityDiffusionTMax"    ,DDTRamp.y));
      Log->Print(fun::VarStr("  DensityDiffusionValueMax",DDTRamp.z));
      cinfo=fun::PrintStr("(%g:t%g:t%g:v%g)",DDTValue,DDTRamp.x,DDTRamp.y,DDTRamp.z);
    } //<vs_ddramp_end>
    ConfigInfo=ConfigInfo+cinfo;
  }
  if(TDensity==DDT_DDT2Full && KernelH/Dp>1.5)Log->PrintWarning(
    "It is advised that selected DDT: \'Fourtakas et al 2019 (full)\' is used with several boundary layers of particles when h/dp>1.5 (2h <= layers*Dp)");
  #ifdef AVAILABLE_DIVCLEAN
  if(DivCleaning){
    Log->Print(fun::VarStr("Div Cleaning","true"));
    Log->Print(fun::VarStr(" Div Cleaning Kp",DivCleanKp));
  }else{
    Log->Print(fun::VarStr("Div Cleaning","false"));
  }
  #endif
  //-Shifting.
  if(Shifting){
    Shifting->VisuConfig();
    ConfigInfo=ConfigInfo+sep+Shifting->GetConfigInfo();
  }
  else if(ShiftingAdv){ //<vs_advshift_ini>
    ShiftingAdv->VisuConfig();
    ConfigInfo=ConfigInfo+sep+ShiftingAdv->GetConfigInfo();
  } //<vs_advshift_end>
  else Log->Print(fun::VarStr("Shifting","None"));
  //-RigidAlgorithm.
  string rigidalgorithm=(!FtCount? "None": GetNameRigidMode(RigidMode));
  Log->Print(fun::VarStr("RigidAlgorithm",rigidalgorithm));
  if(FtCount)ConfigInfo=ConfigInfo+sep+"Ft-"+rigidalgorithm;
  //-Moorings.
  if(Moorings)ConfigInfo=ConfigInfo+sep+"MoorDynPlus";
  //-Other configurations. 
  Log->Print(fun::VarStr("FloatingCount",FtCount));
  if(FtCount)Log->Print(fun::VarStr("FtPause",FtPause));
  if(FtCount)Log->Print(fun::VarStr("FtConstraints",FtConstraints));
  if(FtCount)Log->Print(fun::VarStr("FtIgnoreRadius",FtIgnoreRadius));
  Log->Print(fun::VarKStr("CaseNp",CaseNp));
  Log->Print(fun::VarKStr("CaseNbound",CaseNbound));
  Log->Print(fun::VarKStr("CaseNfixed",CaseNfixed));
  Log->Print(fun::VarKStr("CaseNmoving",CaseNmoving));
  Log->Print(fun::VarKStr("CaseNfloat",CaseNfloat));
  Log->Print(fun::VarKStr("CaseNfluid",CaseNfluid));
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
  Log->Printf("TimeMax=%g",TimeMax);
  Log->Printf("TimePart=%g",TimePart);
  if(TimePartExtra>=0)Log->Printf("TimePartExtra=%g",TimePartExtra);
  Log->Print(fun::VarStr("Gravity",Gravity));
  //-RhopOut limits.
  Log->Print(fun::VarStr("RhopOut",RhopOut));
  if(RhopOut){
    Log->Print(fun::VarStr("RhopOutMin",RhopOutMin));
    Log->Print(fun::VarStr("RhopOutMax",RhopOutMax));
    ConfigInfo=ConfigInfo+sep+fun::PrintStr("RhopOut(%G-%G)",RhopOutMin,RhopOutMax);
  }
  Log->Print(fun::VarStr("WrnPartsOut",WrnPartsOut));
  Log->Print(fun::VarStr("MinFluidStop",MinFluidStop));
  Log->Print(fun::VarKStr("NpfMinimum",NpfMinimum));
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
  Log->Print("    2022. DualSPHysics: from fluid dynamics to multiphysics problems.");
  Log->Print("    Computational Particle Mechanics, 9:867-895. doi: https://doi.org/10.1007/s40571-021-00404-2");
  Log->Print("");
  //-Code implementation:
  if(Mgpu)Log->Print("- New multi-GPU implementation for Smoothed Particle Hydrodynamics on heterogeneous clusters (Dominguez et al., 2013  https://doi.org/10.1016/j.cpc.2013.03.008)");
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
  if(TVisco==VISCO_Laminar || TVisco==VISCO_LaminarSPS)Log->Print("- Viscosity: Laminar + SPS turbulence model (Dalrymple and Rogers, 2006  https://doi.org/10.1016/j.coastaleng.2005.10.004)");
  if(ViscoBoundFactor!=1     )Log->Print("- Viscosity: ViscoBoundFactor coefficient (Barreiro et al., 2014  https://doi.org/10.1371/journal.pone.0111031)");
  //-Kernel functions:
  if(TKernel==KERNEL_Cubic   )Log->Print("- Kernel: Cubic Spline (Monaghan, 1992  https://doi.org/10.1146/annurev.aa.30.090192.002551)");
  if(TKernel==KERNEL_Wendland)Log->Print("- Kernel: Quintic Wendland (Wendland, 1995  https://doi.org/10.1007/BF02123482)");
  //-Time integration scheme: 
  if(TStep==STEP_Verlet    )Log->Print("- Time integration scheme: Verlet (Verlet, 1967  https://doi.org/10.1103/PhysRev.159.98)");
  if(TStep==STEP_Symplectic)Log->Print("- Time integration scheme: Symplectic (Leimkhuler, 1996  https://doi.org/10.1007/978-3-319-16375-8_1)");
  //-Other features:
  if(PeriActive!=0)Log->Print("- Periodic open boundaries (Gomez-Gesteira et al., 2012  https://doi.org/10.1016/j.cageo.2012.02.029)");
  //<vs_meeshdat_ini>
  if(InOut && InOut->Use_InterpolatedVel())Log->Print("- Inflow-outflow with interpolated mesh data (Ruffini et al., 2023  https://doi.org/10.1016/j.oceaneng.2022.113400)");
  //<vs_meeshdat_end>
  if(InOut        )Log->Print("- Inflow-outflow boundary conditions (Tafuni et al., 2018  https://doi.org/10.1016/j.cma.2018.08.004)");
  if(WithFloating )Log->Print("- Floating objects (Canelas et al., 2015  https://doi.org/10.1002/fld.4031)");
  if(UseDEM       )Log->Print("- Coupling SPH-DCDEM (Canelas et al., 2017  https://doi.org/10.1061/(ASCE)HY.1943-7900.0001331)");
  if(UseChrono    )Log->Print("- New coupling with Project Chrono (Martinez-Estevez et al., 2023  https://doi.org/10.1016/j.cpc.2022.108581)");
  if(Moorings     )Log->Print("- Coupling with MoorDynPlus (Dominguez et al., 2019  https://doi.org/10.1016/j.coastaleng.2019.103560)");
  if(Shifting     )Log->Print("- Shifting algorithm (Lind et al., 2012  https://doi.org/10.1016/j.jcp.2011.10.027)");
  if(AccInput     )Log->Print("- External imposed forces (Longshaw and Rogers, 2015  https://doi.org/10.1016/j.advengsoft.2015.01.008)");
  if(FlexStruc    )Log->Print("- Flexible fluid-structure interaction (O'Connor and Rogers, 2021  https://doi.org/10.1016/j.jfluidstructs.2021.103312)"); //<vs_flexstruc>
  //-Wave generation:
  const bool damp=(Damping!=NULL);
  const bool awas=(WaveGen && WaveGen->UseAwas());
  const bool lonw=(WaveGen && !WaveGen->WavesSolitary());
  const bool solw=(WaveGen && WaveGen->WavesSolitary());
  const bool focw=(WaveGen && WaveGen->WavesFocused());
  const bool inow=(InOut && InOut->Use_AwasVel());
  if(lonw        )Log->Print("- Long-crested wave generation (Altomare et al., 2017  https://doi.org/10.1016/j.coastaleng.2017.06.004)");
  if(solw        )Log->Print("- Solitary wave generation (Dominguez et al., 2019  https://doi.org/10.1080/21664250.2018.1560682)");
  if(focw        )Log->Print("- Focused waves theory (Whittaker et al., 2017  https://doi.org/10.1016/j.coastaleng.2016.12.001)");
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
void JSph::LoadDcellParticles(unsigned n,const typecode* code,const tdouble3* pos
  ,unsigned* dcell)const
{
  for(unsigned p=0;p<n;p++){
    typecode codeout=CODE_GetSpecialValue(code[p]);
    if(codeout<CODE_OUTIGNORE){
      const tdouble3 ps=pos[p];
      if(ps>=DomRealPosMin && ps<DomRealPosMax){//-Particle in.
        const double dx=ps.x-DomPosMin.x;
        const double dy=ps.y-DomPosMin.y;
        const double dz=ps.z-DomPosMin.z;
        const unsigned cx=unsigned(dx/Scell);
        const unsigned cy=unsigned(dy/Scell);
        const unsigned cz=unsigned(dz/Scell);
        dcell[p]=DCEL_Cell(DomCellCode,cx,cy,cz);
      }
      else{//-Particle out.
        Run_Exceptioon("Found new particles out."); //-There can not be new particles excluded. | No puede haber nuevas particulas excluidas.
        dcell[p]=DCEL_CodeMapOut;
      }
    }
    else dcell[p]=DCEL_CodeMapOut;
  }
}

//==============================================================================
/// Initializes data of particles according XML configuration.
///
/// Inicializa datos de las particulas a partir de la configuracion en el XML.
//==============================================================================
void JSph::RunInitialize(unsigned np,unsigned npb,const tdouble3* pos
  ,const unsigned* idp,const typecode* code,tfloat4* velrho,tfloat3* boundnor)
{
  if(!PartBegin){
    JDsInitialize init(Simulate2D,Simulate2DPosY,MapRealPosMin,MapRealPosMax
      ,Dp,KernelH,DirCase,CaseNbound,boundnor!=NULL);
    //-Loads configuration from XML.
    {
      JXml xml;
      xml.LoadFile(FileXml);
      if(xml.GetNodeSimple("case.execution.special.initialize",true)){
        init.LoadXml(&xml,"case.execution.special.initialize");
      }
    }
    //-Executes initialize tasks.
    if(init.Count()){
      //-Creates array with mktype value.
      word* mktype=new word[np];
      for(unsigned p=0;p<np;p++){
        const unsigned cmk=MkInfo->GetMkBlockByCode(code[p]);
        mktype[p]=(cmk<MkInfo->Size()? word(MkInfo->Mkblock(cmk)->MkType): USHRT_MAX);
      }
      init.Run(np,npb,pos,idp,mktype,velrho,boundnor);
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
void JSph::CreatePartsInit(unsigned np,const tdouble3* pos,const typecode* code){
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
/// Configure cell map division (defines ScellDiv, Scell, Map_Cells). 
//==============================================================================
void JSph::ConfigCellDivision(){
  if(CellMode!=CELLMODE_Full && CellMode!=CELLMODE_Half)
    Run_Exceptioon("The CellMode is invalid.");
  ScellDiv=(CellMode==CELLMODE_Full? 1: 2);
  Scell=KernelSize/ScellDiv;
  MovLimit=Scell*0.9f;
  Map_Cells=TUint3(unsigned(ceil(Map_Size.x/Scell)),
                   unsigned(ceil(Map_Size.y/Scell)),
                   unsigned(ceil(Map_Size.z/Scell)));
  //-Prints configuration.
  Log->Print(fun::VarStr("CellMode",string(GetNameCellMode(CellMode))));
  Log->Print(fun::VarStr("ScellDiv",ScellDiv));
  const ullong ncells=ullong(Map_Cells.x)*ullong(Map_Cells.y)*ullong(Map_Cells.z);
  Log->Printf("MapCells=(%s)  (%s cells)",fun::Uint3Str(Map_Cells).c_str(),KINT(ncells));
  Log->Print(fun::VarStr("CellDomFixed",CellDomFixed));
  //-Creates VTK file with map cells.
  if(SaveMapCellsVtkSize()<10*MEBIBYTE)SaveMapCellsVtk(Scell);
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
  if(DomCelIni.x>=Map_Cells.x || DomCelIni.y>=Map_Cells.y || DomCelIni.z>=Map_Cells.z )
    Run_Exceptioon("DomCelIni is invalid.");
  if(DomCelFin.x>Map_Cells.x || DomCelFin.y>Map_Cells.y || DomCelFin.z>Map_Cells.z )
    Run_Exceptioon("DomCelFin is invalid.");
  if(DomCells.x<1 || DomCells.y<1 || DomCells.z<1 )
    Run_Exceptioon("The domain of cells is invalid.");
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
  //-Computes codification of cells for the selected domain.
  //-Calcula codificacion de celdas para el dominio seleccionado.
  DomCellCode=JDsDcell::CalcCellCode(DomCells+TUint3(1));
  if(!DomCellCode){
    Run_Exceptioon(fun::PrintStr("Failed to select a valid CellCode for %ux%ux%u cells (CellMode=%s)."
      ,DomCells.x,DomCells.y,DomCells.z,GetNameCellMode(CellMode)));
  }
  //-Prints configurantion.
  const ullong ncells=ullong(DomCells.x)*ullong(DomCells.y)*ullong(DomCells.z);
  Log->Printf("DomCells=(%s)  (%s cells)",fun::Uint3Str(DomCells).c_str(),KINT(ncells));
  Log->Print(fun::VarStr("DomCellCode",JDsDcell::DcellCodeStr(DomCellCode)));
}

//==============================================================================
/// Defines PosCellSize for PosCell on GPU according to Map_Size and KernelSize.
/// Define PosCellSize para PosCell en GPU segun Map_Size y KernelSize.
//==============================================================================
void JSph::ConfigPosCellGpu(){
  //-Checks PosCellCode configuration is valid.
  const unsigned bz=PSCEL_MOVY;
  const unsigned by=PSCEL_MOVX-bz;
  const unsigned bx=32-by-bz;
  const unsigned nx=1<<bx,ny=1<<by,nz=1<<bz;
  Log->Print(fun::VarStr("PosCellCode",fun::PrintStr("%d_%d_%d (%d,%d,%d)",bx,by,bz,nx,ny,nz)));
  if(bx+by+bz!=32)Run_Exceptioon("Number of bits for cells is not 32.");
  const unsigned maskx=(UINT_MAX>>(by+bz))<<(by+bz);
  const unsigned maskz=(UINT_MAX<<(bx+by))>>(bx+by);
  const unsigned masky=(UINT_MAX&(~(maskx|maskz)));
  //Log->Printf("===> X(%02d):[%u]==[%u] cells:",bx,maskx,PSCEL_X);
  //Log->Printf("===> Y(%02d):[%u]==[%u]",by,masky,PSCEL_Y);
  //Log->Printf("===> Z(%02d):[%u]==[%u]",bz,maskz,PSCEL_Z);
  if(maskx!=PSCEL_X)Run_Exceptioon("Mask for cell-X is wrong.");
  if(masky!=PSCEL_Y)Run_Exceptioon("Mask for cell-Y is wrong.");
  if(maskz!=PSCEL_Z)Run_Exceptioon("Mask for cell-Z is wrong.");
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
    const tuint3 scells=JDsDcell::CalcCellDistribution(Map_Cells,32);
    if(scells==TUint3(0) || scells.x+scells.y+scells.z>32)
      Run_Exceptioon("The number of cells is too large for a 32-bit PosCell configuration. The number of cells should be reduced.");
    Log->Printf("\nThe current configuration can be changed by the user by modifying the DualSphDef.h file and compiling the program again. ");
    Log->Printf("Replace the following code in DualSphDef.h:");
    Log->Printf("  //#define PSCEL_CONFIG_USER");
    Log->Printf("  #ifdef PSCEL_CONFIG_USER");
    Log->Printf("     ...");
    Log->Printf("  #else\n");
    Log->Printf("With:");
    Log->Printf("  #define PSCEL_CONFIG_USER");
    Log->Printf("  #ifdef PSCEL_CONFIG_USER");
    const unsigned cex0=(1<<(scells.y+scells.z));
    const unsigned cey0=(1<<(scells.z));
    const unsigned cez0=(1);
    unsigned cex=cex0,cey=cey0,cez=cez0;
    for(unsigned c=1;c<scells.x;c++)cex=(cex<<1)|cex0;
    for(unsigned c=1;c<scells.y;c++)cey=(cey<<1)|cey0;
    for(unsigned c=1;c<scells.z;c++)cez=(cez<<1)|cez0;
    Log->Printf("    #define PSCEL1_X 0x%08x  //-Mask of bits for cell X: %u bits for %u cells",cex,scells.x,1<<scells.x);
    Log->Printf("    #define PSCEL1_Y 0x%08x  //-Mask of bits for cell Y: %u bits for %u cells",cey,scells.y,1<<scells.y);
    Log->Printf("    #define PSCEL1_Z 0x%08x  //-Mask of bits for cell Z: %u bits for %u cells",cez,scells.z,1<<scells.z);
    Log->Printf("    #define PSCEL1_MOVX %7u  //-Displacement to obaint X cell.",scells.y+scells.z);
    Log->Printf("    #define PSCEL1_MOVY %7u  //-Displacement to obaint Y cell.",scells.z);
    Log->Printf("  #else\n");
    Run_Exceptioon("Current configuration for PosCell on GPU is invalid. More information above.");
  }
}

//==============================================================================
/// Computes maximum distance between particles and center of floating.
/// Calcula distancia maxima entre particulas y centro de cada floating.
//==============================================================================
void JSph::CalcFloatingRadius(unsigned np,const tdouble3* pos,const unsigned* idp){
  const float overradius=1.2f; //-Percentage of ration increase. | Porcentaje de incremento de radio. 
  unsigned* ridp=new unsigned[CaseNfloat];
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
    StFloatingData* fobj=FtObjs+cf;
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
void JSph::RestartCheckData(bool loadpsingle){
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
    //-Checks warnings.
    if(loadpsingle)Log->PrintWarning("Restarting from single precision data generates some minor differences in the results.");
  }
}

//==============================================================================
/// Checks the initial density of fluid particles. If some particle is out of 
/// limits throws an exception.
/// Comprueba la densidad inicial de las particulas fluido. Si alguna particula 
/// esta fuera de los limites, lanza una excepcion.
//==============================================================================
void JSph::CheckRhoLimits(){
  const tfloat4*  velrho=PartsLoaded->GetVelRho(); ///<Velocity and density of each particle
  const unsigned* idp   =PartsLoaded->GetIdp();    ///<Identifier of each particle
  const unsigned n=PartsLoaded->GetCount();
  //-Checks the initial density of each fluid particle
  for(unsigned p=0;p<n;p++)if(idp[p]>=CaseNbound){
    if(velrho[p].w<RhopOutMin || RhopOutMax<velrho[p].w)
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
  PartsLoaded->CheckConfig(CaseNp,CaseNfixed,CaseNmoving,CaseNfloat,CaseNfluid
    ,Simulate2D,Simulate2DPosY,TpPeri(PeriActive));

  if(PartBegin)RestartCheckData(PartsLoaded->GetPosSingle());
  Log->Printf("Loaded particles: %s",KINT(PartsLoaded->GetCount()));

  //-Checks if the initial density of fluid particles is out of limits.
  //-Comprueba si la densidad inicial de las particulas fluido esta fuera de los limites.
  CheckRhoLimits();

  //-Collect information of loaded particles.
  //-Recupera informacion de las particulas cargadas.
  CasePosMin=PartsLoaded->GetCasePosMin();
  CasePosMax=PartsLoaded->GetCasePosMax();
  if(PartsLoaded->GetCount()==0)
    Run_Exceptioon("The initial condition does not include any initial particle. Some particle is required to start the simulation.");

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
  //-Defines PosCellSize for PosCell on GPU according to Map_Size and KernelSize.
  if(!Cpu)ConfigPosCellGpu();
  //-Saves initial domain in CfgInit_Domain.vtk (CasePosMin/Max, MapRealPosMin/Max and Map_PosMin/Max).
  SaveInitialDomainVtk();
}

//==============================================================================
/// Initialisation of variables and objects for execution.
/// Inicializa variables y objetos para la ejecucion.
//==============================================================================
void JSph::InitRun(unsigned np,const unsigned* idp,const tdouble3* pos){
  InterStep=(TStep==STEP_Symplectic? INTERSTEP_SymPredictor: INTERSTEP_Verlet);
  VerletStep=0;
  if(FixedDt)DtIni=FixedDt->GetDt(TimeStep,DtIni);
  if(TStep==STEP_Symplectic)SymplecticDtPre=DtIni;
  if(UseDEM)DemDtForce=DtIni; //(DEM)
  
  //-Loads data from other simulation to restart it.
  if(PartBegin){
    PartBeginTimeStep=PartsLoaded->GetPartBeginTimeStep();
    PartBeginTotalNp=PartsLoaded->GetPartBeginTotalNp();
    if(TStep==STEP_Symplectic && PartsLoaded->GetSymplecticDtPre())SymplecticDtPre=PartsLoaded->GetSymplecticDtPre();
    if(UseDEM && PartsLoaded->GetDemDtForce())DemDtForce=PartsLoaded->GetDemDtForce(); //(DEM)
  }
  //-Free memory of PartsLoaded.
  delete PartsLoaded; PartsLoaded=NULL;

  //-Adjust parameters to start.
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
  if(MLPistons){
    MLPistons->SetTimeMod(!PartIni? PartBeginTimeStep: 0);
    MLPistons->CalculateMontionInit(TimeStepIni);
  }

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
  JXml xml;
  xml.LoadFile(FileXml);

  //-Configuration of GaugeSystem.
  GaugeSystem->ConfigCtes(CSP,TimeMax,TimePart,Scell,ScellDiv
    ,Map_PosMin,Map_PosMin,Map_PosMax);
  if(xml.GetNodeSimple("case.execution.special.gauges",true))
    GaugeSystem->LoadXml(&xml,"case.execution.special.gauges",MkInfo);

  //-Prepares WaveGen configuration.
  if(WaveGen){
    Log->Print("Wave paddles configuration:");
    WavesInit(GaugeSystem,MkInfo,TimeMax,TimePart);
    UseWavegenAWAS=WaveGen->UseAwas();
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
    if(PartBegin){
      if(!RestartChrono)Run_Exceptioon("Simulation restart not allowed when Chrono is used.");
      else Log->PrintWarning("Be careful as restart mode is not fully supported by Chrono.");
    }
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

  //-Reconfigures PartFloatSave to include information from ForcePoints.
  if(ForcePoints && ForcePoints->GetPtCount() && PartFloatSave){
    PartFloatSave->ReconfigForcePoints(ForcePoints);
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

  //<vs_flexstruc_ini>
  //-Prepares flexible structures configuration.
  if(FlexStruc){
    Log->Printf("Flexible Structures configuration:");
    FlexStruc->VisuConfig(""," ");
  }
  //<vs_flexstruc_end>

  //-Shows configuration of JGaugeSystem.
  if(GaugeSystem->GetCount())GaugeSystem->VisuConfig("GaugeSystem configuration:"," ");

  //-Shows configuration of JDsOutputTime.
  if(OutputTime->UseSpecialConfig()){
    vector<string> lines;
    OutputTime->GetConfig("TimeOut configuration:"," ",lines);
    Log->Print(lines);
  }

  //-Shows configuration of JDsOutputParts. //<vs_outpaarts_ini>
  if(OutputParts){
    vector<string> lines;
    OutputParts->GetConfig("Output particles filter configuration:"," ",lines);
    Log->Print(lines);
  } //<vs_outpaarts_end>

  //-Adjust variables of floating body particles for restart.
  if(CaseNfloat && PartBegin)InitFloatingsRestart();

  Part=PartIni; Nstep=0; PartNstep=0; PartOut=0;
  TimeStep=TimeStepIni;
  TimeStepM1=TimeStep;
  TimePartNext=(SvAllSteps? TimeStep: OutputTime->GetNextTime(TimeStep));
}

//==============================================================================
/// Adjust variables of floating body particles.
/// Ajusta variables de particulas floating body.
//==============================================================================
void JSph::InitFloatingsRestart(){
  JPartFloatInfoBi4Load ftfile;
  const string filedata=JPartFloatInfoBi4Load::GetFileNameDef(true,PartBeginDir);
  ftfile.LoadFile(filedata);
  //-Checks if the constant data match.
  for(unsigned cf=0;cf<FtCount;cf++){
    const StFloatingData& ft=FtObjs[cf];
    ftfile.CheckHeadData(cf,ft.mkbound,ft.begin,ft.count,ft.mass,ft.massp);
  }
  //-Load PART data. | Carga datos de PART.
  const JPartFloatInfoBi4Data* ftdata=ftfile.GetFtData();
  const unsigned idx=ftdata->GetIdxPart(PartBegin,true);
  for(unsigned cf=0;cf<FtCount;cf++){
    FtObjs[cf].center=ftdata->GetPartCenter(idx,cf);
    FtObjs[cf].fvel  =ftdata->GetPartVelLin(idx,cf);
    FtObjs[cf].fomega=ftdata->GetPartVelAng(idx,cf);
    FtObjs[cf].radius=ftdata->GetHeadRadius(cf);
  }
}

//==============================================================================
/// Returns linear forces from external file according to timestep.
//==============================================================================
void JSph::WavesInit(JGaugeSystem* gaugesystem,const JSphMk* mkinfo
  ,double timemax,double timepart)
{
  StWvgDimensions wdims;
  wdims.dp=GaugeSystem->GetDp();
  wdims.scell=GaugeSystem->GetScell();
  wdims.kernelh=GaugeSystem->GetKernelH();
  wdims.massfluid=GaugeSystem->GetMassFluid();
  wdims.domposmin=GaugeSystem->GetDomPosMin();
  wdims.domposmax=GaugeSystem->GetDomPosMax();
  wdims.padposmin=wdims.padposmax=TDouble3(DBL_MAX);
  const bool ignore_missingmk=(false); //-Enabled for VRes simulations.
  for(unsigned cp=0;cp<WaveGen->GetCount();cp++){
    const word mkbound=WaveGen->GetPaddleMkbound(cp);
    //-Obtains initial limits of paddle for AWAS configuration.
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(mkbound);
    const bool missingmk=(cmk>=mkinfo->Size());
    if(!missingmk){
      wdims.padposmin=mkinfo->Mkblock(cmk)->GetPosMin();
      wdims.padposmax=mkinfo->Mkblock(cmk)->GetPosMax();
    }
    else{
      wdims.padposmin=wdims.padposmax=TDouble3(DBL_MAX);
      if(ignore_missingmk)Log->PrintfWarning(
        "Wave paddle is ignored since particles with mkbound=%d are missing.",mkbound);
    }
    if(!ignore_missingmk || !missingmk){
      WaveGen->InitPaddle(cp,timemax,timepart,wdims);
      //-Configures gauge according AWAS configuration.
      if(WaveGen->PaddleUseAwas(cp)){
        const string gname=fun::PrintStr("AwasMkb%02u",mkbound);
        const double coefmassdef=(Simulate2D? 0.4: 0.5);
        double masslimit,tstart,gdp;
        tdouble3 point0,point2;
        WaveGen->PaddleGetAwasInfo(cp,coefmassdef,wdims.massfluid,masslimit
          ,tstart,gdp,point0,point2);

        //-Creates gauge for AWAS.
        JGaugeSwl* gswl=gaugesystem->AddGaugeSwl(gname,tstart,DBL_MAX,0
          ,false,point0,point2,gdp,float(masslimit));
        WaveGen->PaddleGaugeInit(cp,(void*)gswl);
      }
    }
  }
}

//==============================================================================
/// Loads the last gauge results.
//==============================================================================
void JSph::WavesLoadLastGaugeResults(){
  for(unsigned cp=0;cp<WaveGen->GetCount();cp++){
    const JGaugeSwl* gswl=(const JGaugeSwl*)WaveGen->PaddleGaugeGet(cp);
    if(gswl){
      WaveGen->LoadLastGaugeResults(cp,gswl->GetResult().timestep
        ,gswl->GetResult().posswl,gswl->GetResult().point0);
    }
  }
}

//==============================================================================
/// Updates measurement limits according cell size and last measurement.
//==============================================================================
void JSph::WavesUpdateGaugePoints(){
  for(unsigned cp=0;cp<WaveGen->GetCount();cp++){
    JGaugeSwl* gswl=(JGaugeSwl*)WaveGen->PaddleGaugeGet(cp);
    if(gswl){
      tdouble3 point0=TDouble3(0),point2=TDouble3(0);
      WaveGen->GetUpdateGaugePoints(cp,gswl->GetResult().modified,point0,point2);
      gswl->SetPoints(point0,point2);
    }
  }
}

//==============================================================================
/// Returns linear forces from external file according to timestep.
//==============================================================================
tfloat3 JSph::GetFtExternalForceLin(unsigned cf,double timestep)const{
  return(FtLinearForce[cf]!=NULL? FtLinearForce[cf]->GetValue3f(timestep): TFloat3(0));
}

//==============================================================================
/// Returns angular forces from external file according to timestep.
//==============================================================================
tfloat3 JSph::GetFtExternalForceAng(unsigned cf,double timestep)const{
  return(FtAngularForce[cf]!=NULL? FtAngularForce[cf]->GetValue3f(timestep): TFloat3(0));
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
    if(WaveGen->UseAwas())WavesLoadLastGaugeResults();
    const bool svdata=(TimeStep+stepdt>=TimePartNext);
    for(unsigned c=0;c<WaveGen->GetCount();c++){
      const StMotionData m=(motsim? WaveGen->GetMotion(svdata,c,TimeStep,stepdt):
                                    WaveGen->GetMotionAce(svdata,c,TimeStep,stepdt));
      //Log->PrintfDbg("%u] CalcMotionWaveGen> t:%f  tp:%d  mx:%f",Nstep,TimeStep,m.type,m.linmov.x);
      if(m.type!=MOTT_None){
        if(motsim)DsMotion->SetMotionData   (m);
        else      DsMotion->SetMotionDataAce(m);
      }
    }
    if(WaveGen->UseAwas())WavesUpdateGaugePoints();
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
/// Applies imposed velocity to floating velocity.
/// Aplica velocidad predefinida a la velocidad del floating.
//==============================================================================
void JSph::FtApplyImposedVel(double timestep,int cf,tfloat3& vellin
  ,tfloat3& velang)const
{
  if(!FtObjs[cf].usechrono && (FtLinearVel[cf]!=NULL || FtAngularVel[cf]!=NULL)){
    if(FtLinearVel[cf]!=NULL){
      const tfloat3 v=FtLinearVel[cf]->GetValue3f(timestep);
      if(v.x!=FLT_MAX)vellin.x=v.x;
      if(v.y!=FLT_MAX)vellin.y=v.y;
      if(v.z!=FLT_MAX)vellin.z=v.z;
    }
    if(FtAngularVel[cf]!=NULL){
      const tfloat3 v=FtAngularVel[cf]->GetValue3f(timestep);
      if(v.x!=FLT_MAX)velang.x=v.x;
      if(v.y!=FLT_MAX)velang.y=v.y;
      if(v.z!=FLT_MAX)velang.z=v.z;
    }
  }
}

//==============================================================================
/// Compute new linear and angular acceleration, velocity and center of floating
/// bodies in variables: fto_acelinang[], fto_vellinang[], fto_center[].
//==============================================================================
void JSph::FtComputeAceVel(double dt,bool predictor,bool saveftvalues
  ,tfloat6* fto_acelinang,tfloat6* fto_vellinang,tdouble3* fto_center)
{
  const bool ftpaused=(TimeStep<FtPause);
  if(!ftpaused || saveftvalues){
    //-Get external forces from ForcePoints (moorings).
    memset(Fto_ExForceLinAng,0,sizeof(tfloat6)*FtCount); 
    if(ForcePoints)ForcePoints->GetFtForcesSum(Fto_ExForceLinAng);

    for(unsigned cf=0;cf<FtCount;cf++){
      const StFloatingData& fobj=FtObjs[cf];

      //-Compute total external force (moorings + external file).
      tfloat3 eforcelin=Fto_ExForceLinAng[cf].getlo();
      tfloat3 eforceang=Fto_ExForceLinAng[cf].gethi();
      //-Add predefined forces from external file.
      if(FtLinearForce!=NULL){
        eforcelin=eforcelin+GetFtExternalForceLin(cf,TimeStep);
        eforceang=eforceang+GetFtExternalForceAng(cf,TimeStep);
      }

      //-Compute summation of linear and angular forces starting from acceleration of particles.
      const float fmassp=fobj.massp;
      const tfloat3 fforcelin=fto_acelinang[cf].getlo() * fmassp;
      const tfloat3 fforceang=fto_acelinang[cf].gethi() * fmassp;

      //-Computes total linear acceleration of fluid force + external force + weight due to gravity.
      const float fmass=fobj.mass;
      tfloat3 acelin=(fforcelin + eforcelin + (Gravity*fmass)) / fmass;

      //-Computes total angular acceleration with rotated inertia tensor.
      tfloat3 aceang;
      {
        const tfloat3 fang=fobj.angles;
        tmatrix3f inert=fobj.inertiaini;
        //-Compute a cumulative rotation matrix.
        const tmatrix3f frot=fmath::RotMatrix3x3(fang);
        //-Compute the inertia tensor by rotating the initial tensor to the current orientation I=(R*I_0)*R^T.
        inert=fmath::MulMatrix3x3(fmath::MulMatrix3x3(frot,inert),fmath::TrasMatrix3x3(frot));
        //-Calculates the inverse of the inertia matrix to compute the I^-1 * L= W
        const tmatrix3f invinert=fmath::InverseMatrix3x3(inert);
        //-Calculate total angular acceleration using invinert.
        const tfloat3 forceang=fforceang + eforceang; //-Fluid angular force + external angular force.
        aceang.x=(forceang.x*invinert.a11+forceang.y*invinert.a12+forceang.z*invinert.a13);
        aceang.y=(forceang.x*invinert.a21+forceang.y*invinert.a22+forceang.z*invinert.a23);
        aceang.z=(forceang.x*invinert.a31+forceang.y*invinert.a32+forceang.z*invinert.a33);
      }

      //-Saves computed floating data for FloatingInfo file.
      if(saveftvalues){
        FtObjs[cf].fluforcelin=fforcelin;
        FtObjs[cf].fluforceang=fforceang;
        FtObjs[cf].extforcelin=eforcelin;
        FtObjs[cf].extforceang=eforceang;
        FtObjs[cf].preacelin=acelin;
        FtObjs[cf].preaceang=aceang;
      }

      //-Calculate data to update floatings.
      if(!ftpaused){//-Operator >= is used because when FtPause=0 in symplectic-predictor, code would not enter here. | Se usa >= pq si FtPause es cero en symplectic-predictor no entraria.
        //-Compute new center.
        const tfloat3 fvel0=fobj.fvel;
        tdouble3 fcenter=fobj.center;
        fcenter.x+=                dt*fvel0.x;
        fcenter.y+=(Simulate2D? 0: dt*fvel0.y);
        fcenter.z+=                dt*fvel0.z;
        //-Compute new linear velocity.
        tfloat3 fvel;
        fvel.x=                float(dt*acelin.x + fvel0.x);
        fvel.y=(Simulate2D? 0: float(dt*acelin.y + fvel0.y));
        fvel.z=                float(dt*acelin.z + fvel0.z);
        //-Compute new angular velocity.
        tfloat3 fomega=fobj.fomega;
        fomega.x=(Simulate2D? 0: float(dt*aceang.x + fomega.x));
        fomega.y=                float(dt*aceang.y + fomega.y);
        fomega.z=(Simulate2D? 0: float(dt*aceang.z + fomega.z));

        //-Applies imposed velocity.
        if(FtLinearVel!=NULL)FtApplyImposedVel(TimeStep,cf,fvel,fomega);
        //-Applies motion constraints to computed velocity and acceleration.
        if(fobj.constraints!=FTCON_Free){
          ApplyConstraints(fobj.constraints,acelin,aceang);
          ApplyConstraints(fobj.constraints,fvel,fomega);
        }

        //-Store data to update floating.
        fto_acelinang[cf]=TFloat6(acelin,aceang);
        fto_vellinang[cf]=TFloat6(fvel,fomega);
        fto_center[cf]=fcenter;
      }
    }
    //-Saves linear and angular acceleration of floatings (for debug).
    if(SaveFtAce)SaveFtAceFun(dt,predictor,FtObjs);
  }
}

//==============================================================================
/// Run floating bodiess with Chrono library starting from computed acceleration
/// in fto_acelinang[] and change computed center, linear and angular velocity.
//==============================================================================
void JSph::FtComputeChrono(double dt,bool predictor
  ,const tfloat6* fto_acelinang,tfloat6* fto_vellinang,tdouble3* fto_center)
{
  //-Export data / Exporta datos.
  for(unsigned cf=0;cf<FtCount;cf++)if(FtObjs[cf].usechrono){
    const tfloat6 acelinang=fto_acelinang[cf];
    ChronoObjects->SetFtData(FtObjs[cf].mkbound,acelinang.getlo(),acelinang.gethi());
  }
  //-Applies the external velocities to each floating body of Chrono.
  if(FtLinearVel!=NULL)ChronoFtApplyImposedVel();
  //-Calculate data using Chrono / Calcula datos usando Chrono.
  ChronoObjects->RunChrono(Nstep,TimeStep,dt,predictor);
  //-Load calculated data by Chrono / Carga datos calculados por Chrono.
  for(unsigned cf=0;cf<FtCount;cf++)if(FtObjs[cf].usechrono){
    tfloat3 vellin=TFloat3(0),velang=TFloat3(0);
    ChronoObjects->GetFtData(FtObjs[cf].mkbound,fto_center[cf],vellin,velang);
    fto_vellinang[cf]=TFloat6(vellin,velang);
  }
}

//==============================================================================
/// Update floating data (FtObjs[]) for next step according to new velocity 
/// and new center. 
//==============================================================================
void JSph::FtUpdateFloatings(double dt,const tfloat6* fto_vellinang
  ,const tdouble3* fto_center)
{
  for(unsigned cf=0;cf<FtCount;cf++){
    const tfloat3  fvel   =fto_vellinang[cf].getlo();
    const tfloat3  fomega =fto_vellinang[cf].gethi();
    const tdouble3 fcenter=fto_center[cf];
    //-Stores floating data.
    FtObjs[cf].center=(PeriActive? UpdatePeriodicPos(fcenter): fcenter);
    FtObjs[cf].angles=ToTFloat3(ToTDouble3(FtObjs[cf].angles)+ToTDouble3(fomega)*dt);
    FtObjs[cf].facelin=(fvel  -FtObjs[cf].fvel  )/float(dt);
    FtObjs[cf].faceang=(fomega-FtObjs[cf].fomega)/float(dt);
    FtObjs[cf].fvel=fvel;
    FtObjs[cf].fomega=fomega;
  }
}



//==============================================================================
/// Display a message with reserved memory for the basic data of particles.
/// Muestra un mensaje con la memoria reservada para los datos basicos de las particulas.
//==============================================================================
void JSph::PrintSizeNp(unsigned np,llong size,unsigned allocs)const{
  const double s=double(size)/MEBIBYTE;
  if(Cpu)Log->Printf("**Requested CPU memory for %s particles: %.1f MiB.",KINT(np),s);
  else   Log->Printf("**Requested GPU memory for %s particles: %.1f MiB (%u times).",KINT(np),s,allocs);
}

//==============================================================================
/// Display headers of PARTs
/// Visualiza cabeceras de PARTs.
//==============================================================================
void JSph::PrintHeadPart(){
  Log->Print("PART   PartTime   TotalSteps   Steps    Particles    Cells        Time/Sec   Finish time        ");
  Log->Print("=====  =========  ===========  =======  ===========  ===========  =========  ===================");
  fflush(stdout);
}

//==============================================================================
/// Sets configuration for recordering of particles.
/// Establece configuracion para grabacion de particulas.
//==============================================================================
void JSph::ConfigSaveData(unsigned piece,unsigned pieces,std::string div
  ,unsigned np,const tdouble3* pos,const unsigned* idp)
{
  //-Stores basic information of simulation data.
  JPartDataHead parthead;
  parthead.ConfigBasic(RunCode,AppName,CaseName,CasePosMin,CasePosMax
    ,Simulate2D,Simulate2DPosY,pieces,PartBeginFirst);
  MkInfo->ConfigPartDataHead(&parthead);
  parthead.ConfigCtes(Dp,KernelH,CteB,RhopZero,Gamma,MassBound,MassFluid,Gravity);
  parthead.ConfigSimNp(NpDynamic,ReuseIds);
  parthead.ConfigSimMap(MapRealPosMin,MapRealPosMax);
  parthead.ConfigSimPeri(TpPeriFromPeriActive(PeriActive),PeriXinc,PeriYinc,PeriZinc);
  parthead.ConfigVisco(TVisco,Visco,ViscoBoundFactor);
  if(SvData&SDAT_Binx){
    Log->AddFileInfo(DirDataOut+"Part_Head.ibi4","Binary file with basic information of simulation data.");
    parthead.SaveFile(DirDataOut);
  }
  //-Configures object to store particles and information.  
  //-Configura objeto para grabacion de particulas e informacion.
  if(SvData&SDAT_Info || SvData&SDAT_Binx){
    DataBi4=new JPartDataBi4();
    DataBi4->Config(NoRtimes,piece,pieces,DirDataOut,&parthead);
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
    DataOutBi4->ConfigBasic(NoRtimes,piece,pieces,RunCode,AppName,Simulate2D,DirDataOut);
    DataOutBi4->ConfigParticles(CaseNp,CaseNfixed,CaseNmoving,CaseNfloat,CaseNfluid);
    DataOutBi4->ConfigLimits(MapRealPosMin,MapRealPosMax,(RhopOut? RhopOutMin: 0),(RhopOut? RhopOutMax: 0));
    DataOutBi4->SaveInitial();
    Log->AddFileInfo(DirDataOut+"PartOut_???.obi4","Binary file with particles excluded during simulation (input for PartVtkOut program).");
  }
  //-Configures object to store motion reference data.
  //-Configura objeto para grabacion de datos de referencia de movimiento.
  if(SvData&SDAT_Binx && (CaseNmoving || CaseNfloat)){
    const word mkboundfirst=MkInfo->GetMkBoundFirst();
    const unsigned countmk=MkInfo->CountBlockType(TpPartMoving)+MkInfo->CountBlockType(TpPartFloating);
    PartMotionSave=new JDsPartMotionSave(Cpu||Mgpu,AppName,DirDataOut,CaseNfixed,CaseNmoving,CaseNfloat,mkboundfirst,countmk,TimePart,TimePartExtra);
    for(unsigned cmk=0;cmk<MkInfo->Size();cmk++){
      const JSphMkBlock* pmk=MkInfo->Mkblock(cmk);
      if(pmk->Type==TpPartMoving  )PartMotionSave->ConfigAddMovingMk  (pmk->MkType,pmk->Begin,pmk->Count);
      if(pmk->Type==TpPartFloating)PartMotionSave->ConfigAddFloatingMk(pmk->MkType,pmk->Begin,pmk->Count);
    }
    PartMotionSave->ConfigMotionRefs(np,pos,idp);
    PartMotionSave->SaveInitial();
    Log->AddFileInfo(DirDataOut+PartMotionSave->GetFileName(true)
      ,"Main binary file with motion reference data of moving and floating objects during simulation.");
    if(PartMotionSave->ExtraIsActive()){
      PartMotionSave->AddDataExtraCpu(Part,TimeStep,Nstep,np,NULL,NULL);
      Log->AddFileInfo(DirDataOut+PartMotionSave->GetFileName(false)
        ,"Extra binary file with motion reference data of moving and floating objects during simulation.");
    }
  }
  //-Configures object to store floating data.
  //-Configura objeto para grabacion de datos de objetos flotantes.
  if(SvData&SDAT_Binx && CaseNfloat){
    const word mkboundfirst=MkInfo->GetMkBoundFirst();
    PartFloatSave=new JDsPartFloatSave(Cpu,AppName,DirDataOut,mkboundfirst
      ,FtCount,TimePart,TimePartExtra);
    PartFloatSave->ConfigFtData(FtCount,FtObjs);
    PartFloatSave->SetFtData(Part,TimeStep,Nstep,FtObjs,NULL);
    PartFloatSave->SaveInitial();
    Log->AddFileInfo(DirDataOut+PartFloatSave->GetFileName(true)
      ,"Main binary file with floating data during simulation.");
    if(PartFloatSave->ExtraIsActive()){
      PartFloatSave->AddDataExtra(Part,TimeStep,Nstep);
      Log->AddFileInfo(DirDataOut+PartFloatSave->GetFileName(false)
        ,"Extra binary file with floating data during simulation.");
    }
  }
  TimeOutExtraUpdate(TimeStep);
  //-Configures object to store extra data.
  if(SvData&SDAT_Binx && !SvExtraParts.empty()){
    if(TBoundary==BC_MDBC || !SvPosDouble){
      SvExtraDataBi4=new JDsExtraDataSave(AppName,DirDataOut,CaseNbound,CaseNfloat,Log);
      SvExtraDataBi4->Config(SvExtraParts);
    }
  }
  //-Creates object to store excluded particles until recordering. 
  //-Crea objeto para almacenar las particulas excluidas hasta su grabacion.
  PartsOut=new JDsPartsOut();
}

//==============================================================================
/// Stores new excluded particles until recordering next PART.
/// Almacena nuevas particulas excluidas hasta la grabacion del proximo PART.
//==============================================================================
void JSph::AddParticlesOut(unsigned nout,const unsigned* idp,const tdouble3* pos
  ,const tfloat3* vel,const float* rho,const typecode* code)
{
  PartsOut->AddParticles(nout,idp,pos,vel,rho,code);
}

//==============================================================================
/// Manages excluded particles fixed, moving and floating before aborting the execution.
/// Gestiona particulas excluidas fixed, moving y floating antes de abortar la ejecucion.
//==============================================================================
void JSph::AbortBoundOut(JLog2* log,unsigned nout,const unsigned* idp
  ,const tdouble3* pos,const tfloat3* vel,const float* rho,const typecode* code)
{
  //-Prepares data of excluded boundary particles.
  byte* type=new byte[nout];
  byte* motive=new byte[nout];
  unsigned outfixed=0,outmoving=0,outfloat=0;
  unsigned outpos=0,outrho=0,outmov=0;
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
      case CODE_OUTPOS:  motive[p]=1; outpos++;  break;
      case CODE_OUTRHO:  motive[p]=2; outrho++;  break; 
      case CODE_OUTMOV:  motive[p]=3; outmov++;  break; 
      default:           motive[p]=0;            break; 
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
  npunknown=nout-outpos-outrho-outmov;
  if(!npunknown)log->Printf("Excluded for: position=%u  rho=%u  velocity=%u",outpos,outrho,outmov);
  else log->Printf("Excluded for: position=%u  rho=%u  velocity=%u  UNKNOWN=%u",outpos,outrho,outmov,npunknown);
  if(outxmin)log->Print("Some boundary particle exceeded the -X limit (left limit) of the simulation domain.");
  if(outxmax)log->Print("Some boundary particle exceeded the +X limit (right limit) of the simulation domain.");
  if(outymin)log->Print("Some boundary particle exceeded the -Y limit (front limit) of the simulation domain.");
  if(outymax)log->Print("Some boundary particle exceeded the +Y limit (back limit) of the simulation domain.");
  if(outzmin)log->Print("Some boundary particle exceeded the -Z limit (bottom limit) of the simulation domain.");
  if(outzmax)log->Print("Some boundary particle exceeded the +Z limit (top limit) of the simulation domain.");
  log->Print(" ");
  //-Creates VTK file.
  {
    JDataArrays arrays;
    arrays.AddArray("Pos"   ,nout,pos   ,false);
    arrays.AddArray("Idp"   ,nout,idp   ,false);
    arrays.AddArray("Vel"   ,nout,vel   ,false);
    arrays.AddArray("Rho"   ,nout,rho   ,false);
    arrays.AddArray("Type"  ,nout,type  ,false);
    arrays.AddArray("Motive",nout,motive,false);
    const string file=DirOut+"Error_BoundaryOut.vtk";
    log->AddFileInfo(file,"Saves the excluded boundary particles.");
    JSpVtkData::Save(file,arrays,"Pos");
  }
  //-Aborts execution.
  Run_Exceptioon("Fixed, moving or floating particles were excluded. Check VTK file Error_BoundaryOut.vtk with excluded particles.");
}

//==============================================================================
/// Loads array tdouble3 with position in tdouble2+double.
//==============================================================================
void JSph::Pos21ToPos3(unsigned n,const tdouble2* posxy,const double* posz
  ,tdouble3* pos)const
{
  for(unsigned p=0;p<n;p++){
    pos[p]=TDouble3(posxy[p].x,posxy[p].y,posz[p]);
  }
}

//==============================================================================
/// Loads arrays tdouble2+double with position in tdouble3.
//==============================================================================
void JSph::Pos3ToPos21(unsigned n,const tdouble3* pos,tdouble2* posxy
  ,double* posz)const
{
  for(unsigned p=0;p<n;p++){
    posxy[p]=TDouble2(pos[p].x,pos[p].y);
    posz [p]=pos[p].z;
  }
}

//==============================================================================
/// Loads arrays pos3,vel3,rho from posxy,posz,velrho.
//==============================================================================
void JSph::Pos21Vel4ToPos3Vel31(unsigned n,const tdouble2* posxy
  ,const double* posz,const tfloat4* velrho,tdouble3* pos,tfloat3* vel
  ,float* rho)const
{
  for(unsigned p=0;p<n;p++){
    pos[p]=TDouble3(posxy[p].x,posxy[p].y,posz[p]);
    const tfloat4 vr=velrho[p];
    vel[p]=TFloat3(vr.x,vr.y,vr.z);
    rho[p]=vr.w;
  }
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
void JSph::AddBasicArrays(JDataArrays& arrays,unsigned np,const tdouble3* pos
  ,const unsigned* idp,const tfloat3* vel,const float* rho)const
{
  arrays.AddArray("Pos",np,pos);
  arrays.AddArray("Idp",np,idp);
  arrays.AddArray("Vel",np,vel);
  arrays.AddArray("Rho",np,rho);
}

//==============================================================================
/// Saves RunPARTs.csv with extra information of simulated PARTs.
//==============================================================================
void JSph::SaveRunPartsCsv(const StInfoPartPlus& infoplus,double tpart
  ,double tsim)const
{
  const string file=AppInfo.GetDirOut()+"RunPARTs.csv";
  jcsv::JSaveCsv2 scsv(file,true,false);
  if(!scsv.GetAppendMode()){
    Log->AddFileInfo("RunPARTs.csv","Saves extra information of simulated PARTs.");
    //-Saves head.
    scsv.SetHead();
    scsv << "Part;TimeStep [s];Steps";
    scsv << "DTsMin;PartRuntime [s]";
    scsv << "NpSave;NpSim;NpNew;NpOut";
    scsv << "NctSim;NpAlloc [X];NctAlloc [X]";
    scsv << "SimRuntime [s];NpbSim;NpfSim;NpNormal";
    scsv << "NpOutPos;NpOutRho;NpOutMov";
    scsv << "DtMin [s];DtMax [s]";
    scsv << "MemCPU [MiB];MemGPU [MiB];MemGPU_Cells [MiB]";
    scsv << "NpAlloc;NctAlloc";
    scsv << jcsv::Endl();
  }
  //-Saves data.
  scsv.SetData();
  scsv << jcsv::AutoSepOff();
  const int partnsteps=(Nstep-PartNstep);
  scsv << Part << fun::RealStr(TimeStep) << KINT(partnsteps);
  scsv << KINT(DtModif-PartDtModif) << tpart;
  const unsigned nout=(infoplus.noutpos+infoplus.noutrho+infoplus.noutmov);
  scsv << KINT(infoplus.npsave) << KINT(infoplus.npsim) << KINT(infoplus.npnew) << KINT(nout);
  scsv << KINT(infoplus.nct) << double(infoplus.npsize)/infoplus.npsim << double(infoplus.nctsize)/infoplus.nct;
  scsv << tsim;
  scsv << KINT(infoplus.npbin+infoplus.npbout) << KINT(infoplus.npf) << KINT(infoplus.npnormal);
  scsv << KINT(infoplus.noutpos) << KINT(infoplus.noutrho) << KINT(infoplus.noutmov);
  //scsv << jcsv::Fmt(jcsv::TpDouble1,"%g") << jcsv::Fmt(jcsv::TpDouble1,"%20.12E");
  scsv << fun::RealStr(partnsteps? PartDtMin: 0) << fun::RealStr(partnsteps? PartDtMax: 0);
  const double memgpunp =(!Cpu? double(infoplus.memorynpalloc)/MEBIBYTE: 0);
  const double memgpunct=(!Cpu? double(infoplus.memorynctalloc)/MEBIBYTE: 0);
  scsv << double(infoplus.memorycpualloc)/MEBIBYTE << memgpunp+memgpunct << memgpunct;
  scsv << KINT(infoplus.npsize) << KINT(infoplus.nctsize);
  scsv << jcsv::Endl();
  scsv.SaveData(true);
}

//==============================================================================
/// Saves RunPARTs.csv with extra information of simulated PARTs.
//==============================================================================
void JSph::SaveRunPartsCsvFinal()const{
  const string file=AppInfo.GetDirOut()+"RunPARTs.csv";
  jcsv::JSaveCsv2 scsv(file,true,false);
  if(scsv.GetAppendMode()){
    scsv.SetData();
    scsv << jcsv::AutoSepOff();
    scsv << jcsv::Endl();
    scsv << "# Part:  Number of output file with particles data." << jcsv::Endl();
    scsv << "# TimeStep [s]:  Physical time of simulation." << jcsv::Endl();
    scsv << "# Steps:  Number of calculation steps since the previous PART." << jcsv::Endl();
    scsv << "# DTsMin:  Number of times the Dt is updated with the minimum allowed value (can be 2 times per step when Symplectic is used)." << jcsv::Endl();
    scsv << "# PartRuntime [s]:  Runtime of last PART simulation." << jcsv::Endl();
    scsv << "# NpSave:  Number of selected particles to save." << jcsv::Endl();
    scsv << "# NpSim:  Number of particles used in simulation (normal + periodic particles)." << jcsv::Endl();
    scsv << "# NpNew:  Number of new fluid particles created by inlet conditions since the previous PART." << jcsv::Endl();
    scsv << "# NpOut:  Number of excluded fluid particles since the previous PART." << jcsv::Endl();
    scsv << "# NctSim:  Number of cells used in simulation for particle search." << jcsv::Endl();
    scsv << "# NpAlloc [X]:  Over-allocation memory for particles NpAlloc/NpSim." << jcsv::Endl();
    scsv << "# NctAlloc [X]:  Over-allocation memory for cells NctAlloc/NctSim." << jcsv::Endl();
    scsv << "# SimRuntime [s]:  Total runtime of simulation (initialisation tasks not included)." << jcsv::Endl();
    scsv << "# NpbSim:  Number of fixed and moving particles (includes periodic particles)." << jcsv::Endl();
    scsv << "# NpfSim:  Number of fluid and floating particles (includes periodic particles)." << jcsv::Endl();
    scsv << "# NpNormal:  Total number of particles without periodic particles." << jcsv::Endl();
    scsv << "# NpOutPos:  Number of excluded particles due to invalid position." << jcsv::Endl();
    scsv << "# NpOutRho:  Number of excluded particles due to invalid density." << jcsv::Endl();
    scsv << "# NpOutMov:  Number of excluded particles due to invalid movement." << jcsv::Endl();
    scsv << "# DtMin [s]:  Minimum value of Dt since the previous PART." << jcsv::Endl();
    scsv << "# DtMax [s]:  Maximum value of Dt since the previous PART." << jcsv::Endl();
    scsv << "# MemCPU [MiB]:  Current CPU memory allocated (only includes the most memory intensive objects)." << jcsv::Endl();
    scsv << "# MemGPU [MiB]:  Current GPU memory allocated (only includes the most memory intensive objects)." << jcsv::Endl();
    scsv << "# MemGPU_Cells [MiB]:  Current GPU memory allocated for cells according to NctAlloc." << jcsv::Endl();
    scsv << "# NpAlloc:  Number of supported particles with memory allocated." << jcsv::Endl();
    scsv << "# NctAlloc:  Number of supported cells with memory allocated." << jcsv::Endl();
  }
}
//==============================================================================
/// Stores files of particle data.
/// Graba los ficheros de datos de particulas.
//==============================================================================
void JSph::SavePartData(unsigned npsave,unsigned nout,const JDataArrays& arrays
  ,unsigned ndom,const tdouble6* vdom,const StInfoPartPlus& infoplus)
{
  TimerPart.Stop();
  TimerSim.Stop();
  //-Saves RunPARTs.csv
  if(SvData&SDAT_Info){
    SaveRunPartsCsv(infoplus,TimerPart.GetSecs(),TimerSim.GetSecs());
  }

  //-Stores particle data and other information in bi4 format.
  if(DataBi4){
    tfloat3* posf3=NULL;
    tdouble3 domainmin=vdom[0].getlo();
    tdouble3 domainmax=vdom[0].gethi();
    for(unsigned c=1;c<ndom;c++){
      domainmin=MinValues(domainmin,vdom[c].getlo());
      domainmax=MaxValues(domainmax,vdom[c].gethi());
    }
    JBinaryData* bdpart=DataBi4->AddPartInfo(Part,TimeStep,npsave,nout,Nstep
      ,TimerPart.GetSecs(),domainmin,domainmax,TotalNp);
    if(TStep==STEP_Symplectic)bdpart->SetvDouble("SymplecticDtPre",SymplecticDtPre);
    if(UseDEM)bdpart->SetvDouble("DemDtForce",DemDtForce); //(DEM)
    if(SvData&SDAT_Info){
      bdpart->SetvDouble("dtmean",(!Nstep? 0: (TimeStep-TimeStepM1)/(Nstep-PartNstep)));
      bdpart->SetvDouble("dtmin",(!Nstep? 0: PartDtMin));
      bdpart->SetvDouble("dtmax",(!Nstep? 0: PartDtMax));
      if(FixedDt)bdpart->SetvDouble("dterror",FixedDt->GetDtError(true));
      bdpart->SetvDouble("timesim",infoplus.timesim);
      bdpart->SetvUint("nct"     ,infoplus.nct);
      bdpart->SetvUint("nctsize" ,infoplus.nctsize);
      bdpart->SetvUint("npsim"   ,infoplus.npsim);
      bdpart->SetvUint("npsize"  ,infoplus.npsize);
      bdpart->SetvUint("npnormal",infoplus.npnormal);
      bdpart->SetvUint("npsave"  ,infoplus.npsave);
      bdpart->SetvUint("npnew"   ,infoplus.npnew);
      bdpart->SetvUint("noutpos" ,infoplus.noutpos);
      bdpart->SetvUint("noutrho" ,infoplus.noutrho);
      bdpart->SetvUint("noutmov" ,infoplus.noutmov);
      bdpart->SetvUint("npbin"   ,infoplus.npbin);
      bdpart->SetvUint("npbout"  ,infoplus.npbout);
      bdpart->SetvUint("npf"     ,infoplus.npf);
      bdpart->SetvUint("npbper"  ,infoplus.npbper);
      bdpart->SetvUint("npfper"  ,infoplus.npfper);
      bdpart->SetvLlong("cpualloc",infoplus.memorycpualloc);
      if(infoplus.gpudata){
        bdpart->SetvLlong("nctalloc",infoplus.memorynctalloc);
        bdpart->SetvLlong("nctused" ,infoplus.memorynctused);
        bdpart->SetvLlong("npalloc" ,infoplus.memorynpalloc);
        bdpart->SetvLlong("npused"  ,infoplus.memorynpused);
      }
      if(ndom>1){
        bdpart->SetvUint("subdomain_count",ndom);
        for(unsigned c=0;c<ndom;c++){
          bdpart->SetvDouble3(fun::PrintStr("subdomainmin_%02u",c),vdom[c].getlo());
          bdpart->SetvDouble3(fun::PrintStr("subdomainmax_%02u",c),vdom[c].gethi());
        }
      }
    }
    if(SvData&SDAT_Binx){
      string err;
      if(!(err=arrays.CheckErrorArray("Pos",TypeDouble3,npsave)).empty())Run_Exceptioon(err);
      if(!(err=arrays.CheckErrorArray("Idp",TypeUint   ,npsave)).empty())Run_Exceptioon(err);
      if(!(err=arrays.CheckErrorArray("Vel",TypeFloat3 ,npsave)).empty())Run_Exceptioon(err);
      if(!(err=arrays.CheckErrorArray("Rho",TypeFloat  ,npsave)).empty())Run_Exceptioon(err);
      const tdouble3* pos=arrays.GetArrayDouble3("Pos");
      const unsigned* idp=arrays.GetArrayUint   ("Idp");
      const tfloat3*  vel=arrays.GetArrayFloat3 ("Vel");
      const float*    rho=arrays.GetArrayFloat  ("Rho");
      if(SvPosDouble || (SvExtraDataBi4 && SvExtraDataBi4->CheckSave(Part))){
        DataBi4->AddPartData(npsave,idp,pos,vel,rho);
      }
      else{
        posf3=GetPointerDataFloat3(npsave,pos);
        DataBi4->AddPartData(npsave,idp,posf3,vel,rho);
      }
      //-Adds other arrays.
      const string arrignore=":Pos:Idp:Vel:Rho:";
      for(unsigned ca=0;ca<arrays.Count();ca++){
        const JDataArrays::StDataArray arr=arrays.GetArrayData(ca);
        if(int(arrignore.find(string(":")+arr.keyname+":"))<0){//-Ignore main arrays.
          DataBi4->AddPartData(arr.keyname,npsave,arr.ptr,arr.type);
        }
      }
      DataBi4->SaveFilePart();
    }
    if(SvData&SDAT_Info)DataBi4->SaveFileInfo();
    delete[] posf3;
  }

  //-Stores VTK and/or CSV files.
  if((SvData&SDAT_Csv) || (SvData&SDAT_Vtk)){
    JDataArrays arrays2;
    arrays2.CopyFrom(arrays);

    string err;
    if(!(err=arrays2.CheckErrorArray("Pos" ,TypeDouble3,npsave)).empty())Run_Exceptioon(err);
    if(!(err=arrays2.CheckErrorArray("Idp" ,TypeUint   ,npsave)).empty())Run_Exceptioon(err);
    const tdouble3* pos =arrays2.GetArrayDouble3("Pos");
    const unsigned* idp =arrays2.GetArrayUint   ("Idp");
    //-Generates array with posf3 and type of particle.
    tfloat3* posf3=GetPointerDataFloat3(npsave,pos);
    byte*    type=new byte[npsave];
    for(unsigned p=0;p<npsave;p++){
      const unsigned id=idp[p];
      type[p]=(id>=CaseNbound? 3: (id<CaseNfixed? 0: (id<CaseNpb? 1: 2)));
    }
    arrays2.DeleteArray("Pos");
    arrays2.AddArray("Pos",npsave,posf3);
    arrays2.MoveArray(arrays2.Count()-1,0);
    arrays2.AddArray("Type",npsave,type);
    arrays2.MoveArray(arrays2.Count()-1,4);
    //-Defines fields to be stored.
    if(SvData&SDAT_Vtk){
      JSpVtkData::Save(DirVtkOut+fun::FileNameSec("PartVtk.vtk",Part),arrays2,"Pos");
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
    DataOutBi4->SavePartOut(SvPosDouble,Part,TimeStep,PartsOut->GetCount()
      ,PartsOut->GetIdpOut(),NULL,PartsOut->GetPosOut(),PartsOut->GetVelOut()
      ,PartsOut->GetRhoOut(),PartsOut->GetMotiveOut());
  }

  //-Saves data of floating bodies (and force points) collected in RunFloating().
  if(PartFloatSave){
    PartFloatSave->SaveDataMain(Part,TimeStep,Nstep);
    PartFloatSave->SaveDataExtra();
  }

  //-Empties stock of excluded particles.
  //-Vacia almacen de particulas excluidas.
  PartsOut->Clear();
}

//==============================================================================
/// Generates data output files.
/// Genera los ficheros de salida de datos.
//==============================================================================
void JSph::SaveData(unsigned npsave,const JDataArrays& arrays
  ,unsigned ndom,const tdouble6* vdom,StInfoPartPlus infoplus)
{
  string suffixpartx=fun::PrintStr("_%04d",Part);

  //-Counts new excluded particles.
  const unsigned noutpos=PartsOut->GetOutPosCount();
  const unsigned noutrho=PartsOut->GetOutRhoCount();
  const unsigned noutmov=PartsOut->GetOutMovCount();
  const unsigned nout=noutpos+noutrho+noutmov;
  if(nout!=PartsOut->GetCount())Run_Exceptioon("Excluded particles with unknown reason.");
  AddOutCount(noutpos,noutrho,noutmov);
  infoplus.SetNout(noutpos,noutrho,noutmov);

  //-Stores data files of particles.
  SavePartData(npsave,nout,arrays,ndom,vdom,infoplus);

  //-Reinitialises limits of dt. | Reinicia limites de dt.
  PartDtMin=DBL_MAX; PartDtMax=-DBL_MAX;
  PartDtModif=DtModif;

  //-Computation of time.
  if(Part>PartIni || Nstep){
    TimerPart.Stop();
    const double tpart=TimerPart.GetSecs();
    const double tseg=tpart/(TimeStep-TimeStepM1);
    TimerSim.Stop();
    const double tcalc=TimerSim.GetSecs();
    const double tleft=(tcalc/(TimeStep-TimeStepIni))*(TimeMax-TimeStep);
    const string xparttime=fun::PrintStr((TimeStep>=100? "%9.4f": (TimeStep>=10? "%9.5f": "%9.6f")),TimeStep);
    string xtseg=fun::PrintStr("%9.2f",tseg);
    if(xtseg.size()>9)xtseg=fun::PrintStr("%6.3e",tseg);
    Log->Printf("%05d  %s  %11s  %7s  %11s  %11s  %s  %14s"
      ,Part,xparttime.c_str(),KINT(Nstep),KINT(Nstep-PartNstep)
      ,KINT(infoplus.npsim),KINT(infoplus.nct),xtseg.c_str()
      ,fun::GetDateTimeAfter(int(tleft)).c_str());
  }
  else Log->Printf("Part%s        %s (%.1f%%) particles successfully stored"
    ,suffixpartx.c_str(),KINT(npsave),double(npsave)/infoplus.npnormal*100.);   
  
  //-Shows info of the new inlet particles.
  bool printnp=true;
  if(InOut && InOut->GetNewNpPart()){
    Log->Printf("  Particles new: %s (total new: %s)  -  Current np: %s"
      ,KINT(InOut->GetNewNpPart()),KINT(InOut->GetNewNpTotal()),KINT(infoplus.npsim));
    InOut->ClearNewNpPart();
    printnp=false;
  }
  
  //-Shows info of the excluded particles.
  if(nout){
    PartOut+=nout;
    if(printnp)Log->Printf("  Particles out: %s  (total out: %s)  -  Current np: %s",KINT(nout),KINT(PartOut),KINT(infoplus.npsim));
    else       Log->Printf("  Particles out: %s  (total out: %s)",KINT(nout),KINT(PartOut));
  }

  //-Cheks number of excluded particles.
  if(WrnPartsOut && nout){
    //-Cheks number of excluded particles in one PART.
    if(PartsOutWrn<=100 && nout>=float(max(CaseNfluid,infoplus.npf))*(float(PartsOutWrn)/100.f)){
      Log->PrintfWarning("More than %d%% of current fluid particles were excluded in one PART (t:%g, nstep:%u)"
        ,PartsOutWrn,TimeStep,Nstep);
      if(PartsOutWrn==1)PartsOutWrn=2;
      else if(PartsOutWrn==2)PartsOutWrn=5;
      else if(PartsOutWrn==5)PartsOutWrn=10;
      else PartsOutWrn+=10;
    }
    //-Cheks number of total excluded particles.
    const unsigned noutt=GetOutTotCount();
    if(PartsOutTotWrn<=100 && noutt>=float(TotalNp)*(float(PartsOutTotWrn)/100.f)){
      Log->PrintfWarning("More than %d%% of particles were excluded (t:%g, nstep:%u)"
        ,PartsOutTotWrn,TimeStep,Nstep);
      PartsOutTotWrn+=10;
    }
  }  

  if(SvDomainVtk)SaveDomainVtk(ndom,vdom);
  if(SaveDt)SaveDt->SaveData();
  if(GaugeSystem)GaugeSystem->SaveResults(Part);
  if(ChronoObjects)ChronoObjects->SavePart(Part);
  if(Moorings)Moorings->SaveData(Part);
  if(ForcePoints)ForcePoints->SaveData(Part);
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
    Log->PrintfWarning(
      "TERMINATE file has updated TimeMax from %gs to %gs (current time: %fs)."
      ,TimeMax,tmax,TimeStep);
    TerminateTimeMax=tmax;
    if(!Mgpu)TimeMax=tmax;
  }
  TerminateMt=tmodif;
}

//==============================================================================
/// Generates VTK file with domain of the particles.
/// Genera fichero VTK con el dominio de las particulas.
//==============================================================================
void JSph::SaveDomainVtk(unsigned ndom,const tdouble6* vdom)const{ 
  if(vdom){
    const string file=DirDataOut+fun::FileNameSec("Domain.vtk",Part);
    JSpVtkShape ss;
    ss.AddBoxes(ndom,(const tdouble3*)vdom,KernelH*0.5f);
    ss.SaveVtk(file,"Box");
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
  tfloat3* vdomf3=new tfloat3[nbox*2];
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
  JSpVtkShape ss;
  ss.AddBoxes(nbox,vdomf3,0);
  ss.SaveVtk(file,"Box");
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
  JSpVtkShape ss;
  //-Back lines.
  tdouble3 p0=TDouble3(pmin.x,pmax.y,pmin.z),p1=TDouble3(pmin.x,pmax.y,pmax.z);
  for(unsigned cx=0;cx<=cells.x;cx++)ss.AddLine(p0+TDouble3(scell*cx,0,0),p1+TDouble3(scell*cx,0,0),0);
  p1=TDouble3(pmax.x,pmax.y,pmin.z);
  for(unsigned cz=0;cz<=cells.z;cz++)ss.AddLine(p0+TDouble3(0,0,scell*cz),p1+TDouble3(0,0,scell*cz),0);
  if(!Simulate2D){
    //-Bottom lines.
    p0=TDouble3(pmin.x,pmin.y,pmin.z),p1=TDouble3(pmax.x,pmin.y,pmin.z);
    for(unsigned cy=0;cy<=cells.y;cy++)ss.AddLine(p0+TDouble3(0,scell*cy,0),p1+TDouble3(0,scell*cy,0),1);
    p1=TDouble3(pmin.x,pmax.y,pmin.z);
    for(unsigned cx=0;cx<=cells.x;cx++)ss.AddLine(p0+TDouble3(scell*cx,0,0),p1+TDouble3(scell*cx,0,0),1);
    //-Left lines.
    p0=TDouble3(pmin.x,pmin.y,pmin.z),p1=TDouble3(pmin.x,pmax.y,pmin.z);
    for(unsigned cz=0;cz<=cells.z;cz++)ss.AddLine(p0+TDouble3(0,0,scell*cz),p1+TDouble3(0,0,scell*cz),2);
    p1=TDouble3(pmin.x,pmin.y,pmax.z);
    for(unsigned cy=0;cy<=cells.y;cy++)ss.AddLine(p0+TDouble3(0,scell*cy,0),p1+TDouble3(0,scell*cy,0),2);
  }
  const string file=DirOut+"CfgInit_MapCells.vtk";
  Log->AddFileInfo(file,"Saves the cell division of the simulation domain.");
  ss.SaveVtk(file,"Axis");
}

//==============================================================================
/// Saves VTK file with normals of particles (degug).
/// Only normal (non-periodic) particles are allowed.
/// Graba fichero VTK con normales de las particulas (degug).
/// Solo se permiten particulas normales (no periodicas).
//==============================================================================
void JSph::SaveVtkNormals(std::string filename,int numfile,unsigned np,unsigned npb
  ,const tdouble3* pos,const unsigned* idp,const tfloat3* boundnor,float resize)const
{
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  //-Find floating particles.
  unsigned nfloat=0;
  unsigned* ftidx=NULL;
  if(UseNormalsFt){
    const unsigned size=min(CaseNfloat,np);
    ftidx=new unsigned[size];
    for(unsigned p=npb;p<np;p++)if(idp[p]<CaseNbound)ftidx[nfloat++]=p;
    if(nfloat>CaseNfloat)Run_Exceptioon("More floating particles were found than expected.");
  }
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
  memcpy(vnor,boundnor,sizeof(tfloat3)*npb);
  //-Loads data of floating particles.
  for(unsigned p=0;p<nfloat;p++){
    const unsigned idx=ftidx[p];
    vpos[npb+p]=pos[idx];
    vidp[npb+p]=idp[idx];
    vmk [npb+p]=MkInfo->GetMkById(idp[idx]);
    vnor[npb+p]=boundnor[idx];
  }
  delete[] ftidx; ftidx=NULL;
  //-Computes normalsize.
  if(resize!=1.f)for(unsigned p=0;p<npsel;p++)vnor[p]=vnor[p]*resize;
  for(unsigned p=0;p<npsel;p++)vsnor[p]=fgeo::PointDist(vnor[p]);
  //-Saves VTK file.
  JSpVtkData::Save(filename,arrays,"Pos");
  //-Frees memory.
  arrays.Reset();
}

//==============================================================================
/// Adds basic information of resume to hinfo & dinfo.
/// Anhade la informacion basica de resumen a hinfo y dinfo.
//==============================================================================
void JSph::GetResInfo(float tsim,float ttot,std::string headplus,std::string detplus
  ,std::string& hinfo,std::string& dinfo)const
{
  hinfo=hinfo+"#RunName;Rcode-VersionInfo;DateTime;Np;TSimul;TSeg;TTotal;MemCpu;MemGpu;MemGpuCells";
  dinfo=dinfo+ RunName+ ";"+ RunCode+ "-"+ AppName+ ";"+ RunTimeDate+ ";"+ KINT(CaseNp);
  dinfo=dinfo+ ";"+ fun::FloatStr(tsim)+ ";"+ fun::FloatStr(tsim/float(TimeStep))+ ";"+ fun::FloatStr(ttot);
  dinfo=dinfo+ ";"+ KINT(MaxNumbers.memcpu);
  dinfo=dinfo+ ";"+ KINT(MaxNumbers.memgpu)+ ";"+ KINT(MaxNumbers.memgpunct);
  hinfo=hinfo+";Steps;GPIPS;PhysicalTime;PartFiles;PartsOut;MaxParticles;MaxCells";
  const unsigned nout=GetOutTotCount();
  const string gpips=(DsPips? fun::DoublexStr(DsPips->GetGPIPS(tsim),"%.10f"): "");
  dinfo=dinfo+ ";"+ KINT(Nstep)+ ";"+ gpips+ ";"+ fun::DoublexStr(TimeStep);
  dinfo=dinfo+ ";"+ fun::IntStr(Part)+ ";"+ KINT(nout);
  dinfo=dinfo+ ";"+ KINT(MaxNumbers.particles)+ ";"+ KINT(MaxNumbers.cells);
  hinfo=hinfo+";Hardware;RunMode;Configuration";
  dinfo=dinfo+ ";"+ Hardware+ ";"+ "Cells"+GetNameCellMode(CellMode) + " - " + RunMode+ ";"+ ConfigInfo;
  hinfo=hinfo+";Nbound;Nfixed;Dp;H";
  dinfo=dinfo+ ";"+ KINT(CaseNbound)+ ";"+ KINT(CaseNfixed);
  dinfo=dinfo+ ";"+ fun::FloatStr(float(Dp))+ ";"+ fun::FloatStr(KernelH);
  hinfo=hinfo+";PartsOutRho;PartsOutVel";
  dinfo=dinfo+ ";"+ KINT(GetOutRhoCount())+ ";"+ KINT(GetOutMovCount());
  hinfo=hinfo+ headplus;
  dinfo=dinfo+ detplus;
}

//==============================================================================
/// Generates file Run.csv with resume of execution.
/// Genera fichero Run.csv con resumen de ejecucion.
//==============================================================================
void JSph::SaveRes(float tsim,float ttot,const std::string& headplus
  ,const std::string& detplus)
{
  const string fname=DirOut+"Run.csv";
  Log->AddFileInfo(fname,"One line CSV file with execution parameters and other simulation data.");
  ofstream pf;
  pf.open(fname.c_str());
  if(pf){
    string hinfo,dinfo;
    GetResInfo(tsim,ttot,headplus,detplus,hinfo,dinfo);
    pf << hinfo << endl;
    pf << dinfo << endl;
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
  Log->Printf("Particles of simulation (initial): %s",KINT(CaseNp));
  if(NpDynamic)Log->Printf("Particles of simulation (total)..: %s",KINT(TotalNp));
  if(all){
    Log->Printf("DTs adjusted to DtMin............: %s",KINT(DtModif));
    const unsigned nout=GetOutTotCount();
    Log->Printf("Excluded particles...............: %s",KINT(nout));
    if(GetOutRhoCount())Log->Printf("Excluded particles due to Density: %s",KINT(GetOutRhoCount()));
    if(GetOutMovCount())Log->Printf("Excluded particles due to Velocity: %s",KINT(GetOutMovCount()));
  }
  Log->Printf("Total Runtime....................: %f sec.",ttot);
  Log->Printf("Simulation Runtime...............: %f sec.",tsim);
  if(all){
    float tseg=tsim/float(TimeStep);
    float nstepseg=float(Nstep)/tsim;
    Log->Printf("Runtime per physical second......: %f sec.",tseg);
    //Log->Printf("Time per second of simulation....: %f sec.",tseg);
    Log->Printf("Steps per second.................: %f",nstepseg);
    Log->Printf("Steps of simulation..............: %s",KINT(Nstep));
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
  Log->Printf("Maximum number of particles......: %s",KINT(MaxNumbers.particles));
  Log->Printf("Maximum number of cells..........: %s",KINT(MaxNumbers.cells));
  Log->Printf("CPU Memory.......................: %s (%.2f MiB)"
    ,KINT(MaxNumbers.memcpu),double(MaxNumbers.memcpu)/MEBIBYTE);
  if(MaxNumbers.memgpu)   Log->Printf("GPU Memory.......................: %s (%.2f MiB)"
    ,KINT(MaxNumbers.memgpu),double(MaxNumbers.memgpu)/MEBIBYTE);
  if(MaxNumbers.memgpunct)Log->Printf("GPU Memory (for cells)...........: %s (%.2f MiB)"
    ,KINT(MaxNumbers.memgpunct),double(MaxNumbers.memgpunct)/MEBIBYTE);
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
  else if(tslip==SLIP_FreeSlip)tx="Free-slip";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Returns name of Density Diffusion Term in text format.
/// Devuelve nombre del Density Diffusion Term en texto.
//==============================================================================
std::string JSph::GetDDTName(TpDensity tdensity)const{
  string tx;
       if(tdensity==DDT_None    )tx="None";
  else if(tdensity==DDT_DDT     )tx="Molteni and Colagrossi 2009";
  else if(tdensity==DDT_DDT2    )tx="Fourtakas et al 2019 (inner)";
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
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSph::DgSaveVtkParticlesCpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const tdouble3* pos,const typecode* code
  ,const unsigned* idp,const tfloat4* velrho,const tfloat3* ace)const
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  if(fun::GetExtension(filename).empty())filename=fun::AddExtension(filename,"vtk");
  //-Allocates memory.
  const unsigned np=pfin-pini;
  tfloat3* xpos=new tfloat3[np];
  tfloat3* xvel=new tfloat3[np];
  tfloat3* xace=(ace? new tfloat3[np]: NULL);
  float*   xrho=new float[np];
  byte*    xtype=new byte[np];
  byte*    xkind=new byte[np];
  for(unsigned p=0;p<np;p++){
    xpos[p]=ToTFloat3(pos[p+pini]);
    tfloat4 vr=velrho[p+pini];
    xvel[p]=TFloat3(vr.x,vr.y,vr.z);
    if(xace)xace[p]=ace[p+pini];
    xrho[p]=vr.w;
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
  if(xrho) arrays.AddArray("Rhop",np,xrho ,true);
  if(xace) arrays.AddArray("Ace" ,np,xace ,true);
  JSpVtkData::Save(filename,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Saves VTK file with particle data (degug).
/// Graba fichero VTK con datos de las particulas (degug).
//==============================================================================
void JSph::DgSaveVtkParticlesCpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,const tfloat3* pos,const byte* check
  ,const unsigned* idp,const tfloat3* vel,const float* rho)
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Reserva memoria basica.
  const unsigned n=pfin-pini;
  unsigned* num=new unsigned[n];
  for(unsigned p=0;p<n;p++)num[p]=p;
  //-Generates VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos" ,n,pos+pini,false);
  if(idp)  arrays.AddArray("Idp"  ,n,idp+pini,false);
  if(vel)  arrays.AddArray("Vel"  ,n,vel+pini,false);
  if(rho)  arrays.AddArray("Rho"  ,n,rho+pini,false);
  if(check)arrays.AddArray("Check",n,check+pini,false);
  if(num)  arrays.AddArray("Num"  ,n,num,true);
  JSpVtkData::Save(filename,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Saves CSV file with particle data (degug).
/// Graba fichero CSV con datos de las particulas (degug).
//==============================================================================
void JSph::DgSaveCsvParticlesCpu(std::string filename,int numfile
  ,unsigned pini,unsigned pfin,std::string head,const tfloat3* pos
  ,const unsigned* idp,const tfloat3* vel,const float* rho,const float* ar
  ,const tfloat3* ace,const tfloat3* vcorr)
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
    if(rho)pf << ";Rho";
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
      if(rho)pf << ";" << fun::FloatStr(rho[p],fmt1);
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


#ifdef _WITHMR //<vs_vrres_ini>
void JSph::DgSaveVtkParticlesCpuMR(std::string filename,int numfile,unsigned pini,unsigned pfin
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,const tfloat4 *velrhop
  ,const tfloat3 *ace,const float *arc,const float *pou,const tfloat4 *shiftpos)const
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Allocates memory.
  const unsigned np=pfin-pini;
  unsigned count=0;
  for(unsigned p=0;p<np;p++){
      if(CODE_IsFluid(code[p]) && !CODE_IsFluidBuffer(code[p])&&CODE_IsNormal(code[p])&&!CODE_IsFluidFixed(code[p]))count++;
  }
  tfloat3 *xpos=new tfloat3[count];
  tfloat3 *xvel=new tfloat3[count];
  // tfloat3 *xace=(ace? new tfloat3[count]: NULL);
  // float *xarc=(arc?new float[count]: NULL);
  // float *xpou=(pou? new float[count]:NULL);
  float *xrhop=new float[count];
  unsigned short *xcode= new  unsigned short[count];
  // tfloat3 *xshift=(shiftpos? new tfloat3[count]:NULL);

  unsigned index=0;
  for(unsigned p=0;p<np;p++){
      if(CODE_IsFluid(code[p]) && !CODE_IsFluidBuffer(code[p])&&CODE_IsNormal(code[p])&&!CODE_IsFluidFixed(code[p])){
	xpos[index]=ToTFloat3(pos[p]);
	tfloat4 vr=velrhop[p];
	xvel[index]=TFloat3(vr.x,vr.y,vr.z);
	// tfloat4 sh=(shiftpos?shiftpos[p]: TFloat4(0.0,0.0,0.0,0.0));
	// if(xshift)xshift[index]=TFloat3(sh.x,sh.y,sh.z);
	// if(xace)xace[index]=ace[p];
	xrhop[index]=vr.w;
	xcode[index]=code[p];
	// if(xarc)xarc[index]=arc[p];
	// if(xpou)xpou[index]=pou[p];

	index++;
      }
  }
  //-Generates VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos" ,count,xpos,true);
  if(idp)  arrays.AddArray("Idp" ,count,idp,false);
  if(xcode)arrays.AddArray("Code" ,count,xcode,true);
  if(xvel) arrays.AddArray("Vel" ,count,xvel ,true);
  if(xrhop)arrays.AddArray("Rhop",count,xrhop,true);
  // if(xace) arrays.AddArray("Ace" ,count,xace ,true);
  // if(xarc) arrays.AddArray("Arc" ,count,xarc ,true);
  // if(xpou) arrays.AddArray("Pou" ,count,xpou ,true);
  // if(xshift) arrays.AddArray("Shift" ,count,xshift ,true);
//   delete[] xpos;
//  delete[] xvel;
//  delete[] xrhop;
//  if(xace) delete[] xace;
//  if(xarc) delete[] xarc;
//  if(xpou) delete[] xpou;

  JSpVtkData::Save(filename,arrays,"Pos");
  arrays.Reset();

 
}



void JSph::DgSaveVtkParticlesCpuMRBuffer(std::string filename,int numfile,unsigned pini,unsigned pfin
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,const tfloat4 *velrhop
  ,const tfloat3 *ace,const double *rcond,const tfloat4 *shiftpos)const
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Allocates memory.
  const unsigned np=pfin-pini;
  unsigned count=0;
  for(unsigned p=0;p<np;p++){
      if(CODE_IsFluidBuffer(code[p]))count++;
  }
  tfloat3 *xpos=new tfloat3[count];
  tfloat3 *xvel=new tfloat3[count];
  // tfloat3 *xace=(ace? new tfloat3[count]: NULL);
  // double  *xrcond=(rcond? new double[count]: NULL);
  float *xrhop=new float[count];
  // tfloat3 *xshift=(shiftpos ? new tfloat3[count]:NULL);

  unsigned index=0;
  for(unsigned p=0;p<np;p++){
      if(CODE_IsFluidBuffer(code[p])){
	xpos[index]=ToTFloat3(pos[p]);
	tfloat4 vr=velrhop[p];
	xvel[index]=TFloat3(vr.x,vr.y,vr.z);
	tfloat4 sh=(shiftpos?shiftpos[p]: TFloat4(0.0,0.0,0.0,0.0));
	// if(xshift)xshift[index]=TFloat3(sh.x,sh.y,sh.z);
	// if(xace)xace[index]=ace[p];
	// if(xrcond)xrcond[index]=rcond[p];
	xrhop[index]=vr.w;
	index++;
      }
  }
  //-Generates VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos" ,count,xpos,true);
  if(idp)  arrays.AddArray("Idp" ,count,idp,false);
  if(xvel) arrays.AddArray("Vel" ,count,xvel ,true);
  if(xrhop)arrays.AddArray("Rhop",count,xrhop,true);
  // if(xace) arrays.AddArray("Ace" ,count,xace ,true);
  // if(xrcond) arrays.AddArray("Rcond" ,count,xrcond ,true);
  // if(xshift) arrays.AddArray("Shift" ,count,xshift ,true);
  JSpVtkData::Save(filename,arrays,"Pos");
  arrays.Reset();

//  delete[] xpos;
//    delete[] xvel;
//    delete[] xrhop;
//    if(xace) delete[] xace;
//  if(xrcond) delete [] xrcond;
}


void JSph::DgSaveVtkParticlesCpuMRBound(std::string filename,int numfile,unsigned pini,unsigned pfin
  ,const tdouble3 *pos,const typecode *code,const unsigned *idp,const tfloat4 *velrhop
  ,const tfloat3 *ace)const
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)filename=string("p")+fun::IntStr(mpirank)+"_"+filename;
  if(numfile>=0)filename=fun::FileNameSec(filename,numfile);
  filename=DirDataOut+filename;
  //-Allocates memory.
  const unsigned np=pfin-pini;
  unsigned count=0;
  for(unsigned p=0;p<np;p++){
      if(!CODE_IsFluid(code[p]))count++;
  }
  tfloat3 *xpos=new tfloat3[count];
  tfloat3 *xvel=new tfloat3[count];
  tfloat3 *xace=(ace? new tfloat3[count]: NULL);
  float *xrhop=new float[count];

  unsigned index=0;
  for(unsigned p=0;p<np;p++){
      if(!CODE_IsFluid(code[p])){
	xpos[index]=ToTFloat3(pos[p]);
	tfloat4 vr=velrhop[p];
	xvel[index]=TFloat3(vr.x,vr.y,vr.z);
	if(xace)xace[index]=ace[p];
	xrhop[index]=vr.w;
	index++;
      }
  }
  //-Generates VTK file.
  JDataArrays arrays;
  arrays.AddArray("Pos" ,count,xpos,true);
  if(idp)  arrays.AddArray("Idp" ,count,idp,false);
  if(xvel) arrays.AddArray("Vel" ,count,xvel ,true);
  if(xrhop)arrays.AddArray("Rhop",count,xrhop,true);
  if(xace) arrays.AddArray("Ace" ,count,xace ,true);
  JSpVtkData::Save(filename,arrays,"Pos");
  arrays.Reset();

//  delete[] xpos;
//    delete[] xvel;
//    delete[] xrhop;
//    if(xace) delete[] xace;
}
#endif        //<vs_vrres_end>


