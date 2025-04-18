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

/// \file JSphCfgRun.cpp \brief Implements the class \ref JSphCfgRun.

#include "JSphCfgRun.h"
#include "JAppInfo.h"
#include "JDsphConfig.h"
#include <cfloat>

using namespace std;

//##############################################################################
//# JSphCfgRun
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphCfgRun::JSphCfgRun():JCfgRunBase(){
  ClassName="JSphCfgRun";
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphCfgRun::Reset(){
  PrintInfo=false;
  DirsDef=0;
  Cpu=true;
  Gpu=false;
  GpuId=-1;
  Stable=false;
  SvPosDouble=-1;
  SvExtraParts="undefined";
  OmpThreads=0;
  SvTimers=true;
  CellMode=CELLMODE_Full;
  CellDomFixed=false;
  TBoundary=-1;
  SlipMode=-1;
  NoPenetration=false;
  DomainMode=0;
  DomainFixedMin=DomainFixedMax=TDouble3(0);
  TStep=STEP_None;
  VerletSteps=-1;
  TKernel=KERNEL_None;
  TVisco=VISCO_None;
  Visco=0;
  ViscoBoundFactor=-1;
  TimeMax=-1;
  TimePart=-1;
  TimePartExtra=DBL_MAX;
  TDensity=-1;
  DDTValue=-1;
  DDTValueTRamp=DDTValueTMax=DDTValueMax=0;  //<vs_ddramp>
  Shifting=-1;
  ShiftAdvALE=false; //<vs_advshift>
  ShiftAdvNCP=false; //<vs_advshift>
  Sv_Binx=true; 
  Sv_Info=true;
  Sv_Vtk=false;
  Sv_Csv=false;
  SvNormals=false; 
  SvRes=true; 
  SvDomainVtk=false;
  CaseName="";
  RunName="";
  DirOut="";
  DirDataOut="data"; 
  PartBegin=0;
  PartBeginFirst=0;
  PartBeginDir="";
  RestartChrono=false;
  CFLnumber=-1;
  RhopOutModif=false;
  RhopOutMin=700;
  RhopOutMax=1300;
  FtPause=-1;
  NstepsBreak=0;
  SvAllSteps=false;
  NoRtimes=true;
  PipsMode=0; PipsSteps=100;
  CreateDirs=true;
  CsvSepComa=false;
  #ifdef _WITHMR //<vs_vrres_ini>
    VRes=-1;
    VResOrder=-1;
    VResMethod=-1;
    MRFastSingle=true;
    VResThreshold=-1;
  #endif         //<vs_vrres_end>

}

//==============================================================================
/// Shows information about execution parameters.
//==============================================================================
void JSphCfgRun::VisuInfo()const{
  if(!FeatureList.empty())printf("Available features: %s.\n",FeatureList.c_str());
  printf("\n");
/////////|---------1---------2---------3---------4---------5---------6---------7--------X8
  printf("Information about execution parameters:\n");
  printf("\n");
  printf("  DualSPHysics [name_case [dir_out]] [options]\n");
  printf("\n");
  printf("  General options:\n");
  printf("    -h          Shows information about parameters\n");
  printf("    -ver        Shows version information\n");
  printf("    -info       Shows version features in JSON format\n");
  printf("    -opt <file> Loads a file configuration\n");
  printf("\n");

  printf("  Execution options for CPU:\n");
  printf("    -cpu        Execution on CPU (option by default)\n");
#ifdef OMP_USE
  printf("    -ompthreads:<int>  Only for CPU execution, indicates the number of threads\n");
  printf("                by host for parallel execution, this takes the number of cores\n");
  printf("                of the device by default (or using zero value)\n");
#endif
  printf("\n");

#ifdef _WITHGPU
  printf("  Execution options for Single-GPU:\n");
  printf("    -gpu[:id]   Execution on GPU and optional id of the device\n");
  printf("\n");
#endif

  printf("  General execution options:\n");
  printf("    -stable     The result is always the same but the execution is slower\n");
  printf("    -saveposdouble:<0/1>  Saves position using double precision (default=0)\n");
  printf("    -svextraparts:<int>  PART interval for saving extra data (default=0)\n");
  printf("    -svextraparts:<list> List of PARTs for saving extra data (default=0)\n");
  printf("\n");
  printf("    -cellmode:<mode>  Specifies the cell division mode\n");
  printf("        full      Lowest and the least expensive in memory (by default)\n");
  printf("        half      Fastest and the most expensive in memory\n");
  printf("    -cellfixed:<0/1>  Cell domain is fixed according maximum domain size\n");
  printf("                      (default=0)\n");
  printf("\n");

  printf("  Formulation options:\n");
  printf("    -dbc                    Dynamic Boundary Condition DBC (by default)\n");
  printf("    -mdbc                   Modified Dynamic Boundary Condition mDBC (vel=0 mode)\n");
  printf("    -mdbc_noslip[:nopen]    Modified Dynamic Boundary Condition mDBC (no-slip mode)\n");
  printf("    -mdbc_freeslip[:nopen]  Modified Dynamic Boundary Condition mDBC (free-slip mode)\n");
/////////|---------1---------2---------3---------4---------5---------6---------7--------X8
  printf("\n");
  printf("    -initnorpla:<inlinecfg>  Initialize definition for <boundnormal_plane>\n");
  printf("    -initnorpart:<inlinecfg> Initialize definition for <boundnormal_parts>\n");
  printf("\n");
  printf("    -symplectic      Symplectic algorithm as time step algorithm\n");
  printf("    -verlet[:steps]  Verlet algorithm as time step algorithm and number of\n");
  printf("                     time steps to switch equations\n");
  printf("\n");
  printf("    -wendland        Wendland kernel (by default)\n");
#ifndef DISABLE_KERNELS_EXTRA
  printf("    -cubic           Cubic spline kernel\n");
#endif
  printf("\n");
  printf("    -viscoart:<float>          Artificial viscosity [0-1]\n");
  printf("    -viscolam:<float>          Laminar viscosity [order of 1E-6]\n");  
  printf("    -viscolamsps:<float>       Laminar+SPS viscosity [order of 1E-6]\n");  
  printf("    -viscoboundfactor:<float>  Multiplies the viscosity value of boundary\n");
  printf("\n");
  printf("    -ddt:<mode> Specifies the Density Diffusion Term to correct density\n");
  printf("        none       Not used (by default)\n");
  printf("        1          Diffusion term by Molteni and Colagrossi 2009\n");
  printf("        2          Diffusion term by Fourtakas et al 2019 (inner fluid particles)\n");
  printf("        3          Diffusion term by Fourtakas et al 2019 (all fluid particles)\n");
  printf("    -ddtvalue:<float> Constant for DDT (0.1 by default)\n");
  printf("    -ddtramp:tramp:tmax:maxvalue  Total time of DDT ramp and time for maxvalue\n"); //<vs_ddramp>
  printf("\n");
  printf("    -shifting:<mode> Set Shifting correction (with default parameters)\n");
  printf("        none       Shifting is disabled (by default)\n");
  printf("        nobound    Shifting is not applied near boundary\n");
  printf("        nofixed    Shifting is not applied near fixed boundary\n");
  printf("        full       Shifting is always applied\n");
  printf("        fulladv    Advanced shifting for free-surface (mode:ale:ncp)\n"); //<vs_advshift>
  printf("\n");

  printf("  Simulation options:\n");
  printf("    -name <string>      Specifies path and name of the case \n");
  printf("    -runname <string>   Specifies name for case execution\n");
  printf("    -dirout <dir>       Specifies the general output directory \n");
  printf("    -dirdataout <dir>   Specifies the output subdirectory for binary data\n");
  printf("                        files (default=data) \n");
  printf("\n");
  printf("    -partbegin:begin[:first] dir \n");
  printf("     Specifies the beginning of the simulation starting from a given PART\n");
  printf("     (begin) and located in the directory (dir), (first) indicates the\n");
  printf("     number of the first PART to be generated\n");
  printf("    -restartchrono:<0/1>    Allows restart with Chrono active (default=0)\n");
  printf("\n");
  printf("    -tmax:<float>   Maximum time of simulation\n");
  printf("    -tout:<float>   Time between output files\n");
  printf("    -toutx:<float>  Time between extra output files on motion and floatings\n");
  printf("\n");
  printf("    -cfl:<float> CFL number coefficient to multiply dt\n");
  printf("    -ftpause:<float> Time to start floating bodies movement. By default 0\n");
  printf("    -rhopout:min:max Excludes fluid particles out of these density limits\n");
  printf("    -domain_fixed:xmin:ymin:zmin:xmax:ymax:zmax    The domain is fixed\n");
  printf("     with the specified values\n");
  printf("\n");

  printf("  Output options:\n");
  printf("    -sv:[formats,...] Specifies the output formats.\n");
  printf("        none    No particles files are generated\n");
  printf("        binx    Binary files (by default)\n");
  printf("        info    Information about execution in .ibi4 format (by default)\n");
  printf("        vtk     VTK files\n");
  printf("        csv     CSV files\n");
  printf("    -svnormals:<0/1> Saves normal vector of boundary particles (default=0)\n");
  printf("    -svres:<0/1>     Generates file that summarises the execution process\n");
  printf("    -svtimers:<0/1>  Obtains timing for each individual process\n");
  printf("    -svdomainvtk:<0/1>  Generates VTK file with domain limits\n");
/////////|---------1---------2---------3---------4---------5---------6---------7--------X8
  printf("    -svpips:<mode>:n  Compute PIPS of simulation each n steps (100 by default),\n");
  printf("       mode options: 0=disabled (by default), 1=no save details, 2=save details\n");
  printf("\n");
  printf("    -createdirs:<0/1> Creates full path for output files\n");
  printf("                      (value by default is read from DsphConfig.xml or 1)\n");
  printf("    -csvsep:<0/1>     Separator character in CSV files (0=semicolon, 1=coma)\n");
  printf("                      (value by default is read from DsphConfig.xml or 0)\n");
  printf("\n");

  printf("  Debug options:\n");
  printf("    -nsteps:<uint>  Maximum number of steps allowed (activates nortimes)\n");
  printf("    -svsteps:<0/1>  Saves a PART for each step (activates nortimes)\n");
  printf("    -nortimes:<0/1> Removes execution dependent values from bi4 files\n");
  printf("\n");

  printf("  Variable resolution options:\n");
  printf("    -vres                     Execution with variable resolution algorithm\n");
  printf("    -vres_fast:<0/1>          Precision of variable resolution interpolation (0=double, 1=single)\n");
  printf("    -vres_order:<uint>        Maximum order vres interpolation (0=0th, 1:1st, 2:2nd)\n");
  printf("    -vres_threshold:<float>   Threshold for vres interpolation\n");
  printf("\n");

  printf("  Examples:\n");
  printf("    DualSPHysics case out_case -sv:binx,csv \n");
  printf("    DualSPHysics -name case -dirout out_case -sv:binx,csv \n");
}

//==============================================================================
/// Shows current configuration.
//==============================================================================
void JSphCfgRun::VisuConfig()const{
  printf("\nConfiguration of execution:\n");
  string ln="\n";
  fun::PrintVar("  CaseName",CaseName,ln);
  fun::PrintVar("  RunName",RunName,ln);
  fun::PrintVar("  DirOut",DirOut,ln);
  fun::PrintVar("  DirDataOut",DirDataOut,ln);
  fun::PrintVar("  PartBegin",PartBegin,ln);
  fun::PrintVar("  PartBeginFirst",PartBeginFirst,ln);
  fun::PrintVar("  PartBeginDir",PartBeginDir,ln);
  fun::PrintVar("  Cpu",Cpu,ln);
  printf("  %s  %s\n",fun::VarStr("Gpu",Gpu).c_str(),fun::VarStr("GpuId",GpuId).c_str());
  fun::PrintVar("  Stable",Stable,ln);
  fun::PrintVar("  SvPosDouble",SvPosDouble,ln);
  fun::PrintVar("  OmpThreads",OmpThreads,ln);
  fun::PrintVar("  CellMode",GetNameCellMode(CellMode),ln);
  fun::PrintVar("  TStep",TStep,ln);
  fun::PrintVar("  VerletSteps",VerletSteps,ln);
  fun::PrintVar("  TKernel",TKernel,ln);
  fun::PrintVar("  TVisco",TVisco,ln);
  fun::PrintVar("  Visco",Visco,ln);
  fun::PrintVar("  ViscoBoundFactor",ViscoBoundFactor,ln);
  fun::PrintVar("  TDensity",TDensity,ln);
  fun::PrintVar("  DDTValue",DDTValue,ln);
  fun::PrintVar("  Shifting",Shifting,ln);
  fun::PrintVar("  SvRes",SvRes,ln);
  fun::PrintVar("  SvTimers",SvTimers,ln);
  fun::PrintVar("  SvDomainVtk",SvDomainVtk,ln);
  fun::PrintVar("  Sv_Binx",Sv_Binx,ln);
  fun::PrintVar("  Sv_Info",Sv_Info,ln);
  fun::PrintVar("  Sv_Vtk",Sv_Vtk,ln);
  fun::PrintVar("  Sv_Csv",Sv_Csv,ln);
  fun::PrintVar("  RhopOutModif",RhopOutModif,ln);
  if(RhopOutModif){
    fun::PrintVar("  RhopOutMin",RhopOutMin,ln);
    fun::PrintVar("  RhopOutMax",RhopOutMax,ln);
  }
  fun::PrintVar("  TimeMax",TimeMax,ln);
  fun::PrintVar("  TimePart",TimePart,ln);
  if(DomainMode==2){
    fun::PrintVar("  DomainFixedMin",DomainFixedMin,ln);
    fun::PrintVar("  DomainFixedMax",DomainFixedMax,ln);
  }
  fun::PrintVar("  FtPause",FtPause,ln);
}

//==============================================================================
/// Loads execution parameters.
//==============================================================================
void JSphCfgRun::LoadOpts(const std::string* optlis,int optn,int lv
  ,const std::string& file)
{
  if(lv>=10)Run_Exceptioon("No more than 10 levels of recursive configuration.");
  for(int c=0;c<optn;c++){
    const string opt=optlis[c];
    if(opt[0]!='-' && opt[0]!='#'){
      if(!DirsDef){ CaseName=opt; DirsDef++; }
      else if(DirsDef==1){ DirOut=opt; DirsDef++; }
      else ErrorParm(opt,c,lv,file);
    }
    else if(opt[0]=='-'){
      //-Splits options in txoptfull, txopt1, txopt2, txopt3 and txopt4.
      string txword,txoptfull,txopt1,txopt2,txopt3,txopt4;
      SplitsOpts(opt,txword,txoptfull,txopt1,txopt2,txopt3,txopt4);
      //-Checks keywords in commands.
      if(txword=="CPU"){
        Cpu=Gpu=false; 
        Cpu=true;
      }
      else if(txword=="GPU"){ 
        Cpu=Gpu=false; 
        Gpu=true;
        if(txoptfull!="")GpuId=atoi(txoptfull.c_str()); 
      }

      else if(txword=="STABLE")Stable=OptIsEnabled(txoptfull);
      else if(txword=="SAVEPOSDOUBLE"){
        const int v=(txoptfull!=""? atoi(txoptfull.c_str()): 1);
        SvPosDouble=(!v? 0: 1);
      }
      else if(txword=="SVEXTRAPARTS")SvExtraParts=txoptfull;
#ifdef OMP_USE
      else if(txword=="OMPTHREADS"){ 
        OmpThreads=atoi(txoptfull.c_str());
        if(OmpThreads<0)OmpThreads=0;
      } 
#endif
      else if(txword=="CELLMODE"){
        bool ok=true;
        if(!txoptfull.empty()){
          txoptfull=fun::StrUpper(txoptfull);
          if(txoptfull=="HALF" || txoptfull=="H")CellMode=CELLMODE_Half;
          else if(txoptfull=="FULL" || txoptfull=="2H")CellMode=CELLMODE_Full;
          else ok=false;
        }
        else ok=false;
        if(!ok)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="CELLFIXED")CellDomFixed=OptIsEnabled(txoptfull);
      else if(txword=="DBC"){
        TBoundary=int(BC_DBC);
        SlipMode=int(SLIP_None);
      }
      else if(txword=="MDBC"){
        TBoundary=int(BC_MDBC);
        SlipMode=int(SLIP_Vel0);
      }
      else if(txword=="MDBC_NOSLIP"){
        TBoundary=int(BC_MDBC);
        SlipMode=int(SLIP_NoSlip);
       if(txoptfull!="")NoPenetration=OptIsEnabled(txoptfull);
      }
      else if(txword=="MDBC_FREESLIP"){
        TBoundary=int(BC_MDBC);
        SlipMode=int(SLIP_FreeSlip);
        if(txoptfull!="")NoPenetration=OptIsEnabled(txoptfull);
      }
      else if(txword=="SYMPLECTIC")TStep=STEP_Symplectic;
      else if(txword=="VERLET"){
        TStep=STEP_Verlet; 
        if(txoptfull!="")VerletSteps=atoi(txoptfull.c_str()); 
      }
      else if(txword=="WENDLAND")TKernel=KERNEL_Wendland;
      else if(txword=="CUBIC")TKernel=KERNEL_Cubic;
      else if(txword=="VISCOART"){ 
        Visco=float(atof(txoptfull.c_str())); 
        if(Visco>10)ErrorParm(opt,c,lv,file);
        TVisco=VISCO_Artificial;
      }
      else if(txword=="VISCOLAM"){ 
        Visco=float(atof(txoptfull.c_str())); 
        if(Visco>0.001f)ErrorParm(opt,c,lv,file);
        TVisco=VISCO_Laminar;
      }
      else if(txword=="VISCOLAMSPS"){ 
        Visco=float(atof(txoptfull.c_str())); 
        if(Visco>0.001f)ErrorParm(opt,c,lv,file);
        TVisco=VISCO_LaminarSPS;
      }
      else if(txword=="VISCOBOUNDFACTOR"){ 
        ViscoBoundFactor=float(atof(txoptfull.c_str())); 
        if(ViscoBoundFactor<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="DDT"){
        const int v=atoi(txoptfull.c_str());
        if(v>=0 && v<=3)TDensity=TpDensity(v);
        else ErrorParm(opt,c,lv,file);
      }
      else if(txword=="DDTVALUE"){
        DDTValue=float(atof(txoptfull.c_str())); 
        if(DDTValue<0 || DDTValue>1)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="DDTRAMP"){ //<vs_ddramp_ini>
        tdouble3 v;
        LoadDouble3(txoptfull,0,v);
        DDTValueTRamp=v.x;
        DDTValueTMax=v.y;
        DDTValueMax=v.z;
        if(DDTValueTRamp<=0 || DDTValueMax<=0)DDTValueTRamp=DDTValueTMax=DDTValueMax=0;
        if(DDTValueTMax>DDTValueTRamp)DDTValueTMax=DDTValueTRamp;
      } //<vs_ddramp_end>
      else if(txword=="SHIFTING"){
        const string tx=fun::StrUpper(txopt1);
        if(tx=="NONE")Shifting=0;
        else if(tx=="NOBOUND")Shifting=1;
        else if(tx=="NOFIXED")Shifting=2;
        else if(tx=="FULL"   )Shifting=3;
        else if(tx=="FULLADV"){ //<vs_advshift_ini>
          Shifting=4;
          if(!txopt2.empty())ShiftAdvALE=OptIsEnabled(txopt2);
          if(!txopt3.empty())ShiftAdvNCP=OptIsEnabled(txopt3);
        } //<vs_advshift_end>
        else ErrorParm(opt,c,lv,file);
      }
      else if(txword=="SVNORMALS")SvNormals=OptIsEnabled(txoptfull);
      else if(txword=="SVRES")SvRes=OptIsEnabled(txoptfull);
      else if(txword=="SVTIMERS")SvTimers=OptIsEnabled(txoptfull);
      else if(txword=="SVDOMAINVTK")SvDomainVtk=OptIsEnabled(txoptfull);
      else if(txword=="SV"){
        string txop=fun::StrUpper(txoptfull);
        while(!txop.empty()){
          string op=fun::StrSplit(",",txop);
          if(op=="NONE")Sv_Binx=Sv_Info=Sv_Csv=Sv_Vtk=false;
          else if(op=="BINX" || op=="BIN")Sv_Binx=true;
          else if(op=="INFO" || op=="INF")Sv_Info=true;
          else if(op=="VTK")Sv_Vtk=true;
          else if(op=="CSV")Sv_Csv=true;
          else ErrorParm(opt,c,lv,file);
        }
      }
      else if(txword=="CREATEDIRS")CreateDirs=OptIsEnabled(txoptfull);
      else if(txword=="CSVSEP")CsvSepComa=OptIsEnabled(txoptfull);
      else if(txword=="NAME" && c+1<optn){ CaseName=optlis[c+1]; c++; }
      else if(txword=="RUNNAME" && c+1<optn){ RunName=optlis[c+1]; c++; }
      else if(txword=="DIROUT" && c+1<optn){ DirOut=optlis[c+1]; c++; }
      else if(txword=="DIRDATAOUT" && c+1<optn){ 
        DirDataOut=optlis[c+1]; c++;
        if(DirDataOut==".")DirDataOut="";
      }
      else if(txword=="PARTBEGIN" && c+1<optn){ 
        int v1=atoi(txopt1.c_str());
        int v2=atoi(txopt2.c_str());
        if(v1<0 || v2<0)ErrorParm(opt,c,lv,file);
        else{
          PartBegin=unsigned(v1);
          PartBeginFirst=(txopt2.empty()? PartBegin: unsigned(v2));
        }
        PartBeginDir=optlis[c+1]; c++; 
      }
      else if(txword=="RESTARTCHRONO")RestartChrono=OptIsEnabled(txoptfull);
      else if(txword=="RHOPOUT"){ 
        RhopOutMin=float(atof(txopt1.c_str())); 
        RhopOutMax=float(atof(txopt2.c_str())); 
        RhopOutModif=true;
      }
      else if(txword=="CFL"){ 
        CFLnumber=atof(txoptfull.c_str()); 
        if(CFLnumber<=0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="FTPAUSE"){ 
        FtPause=float(atof(txoptfull.c_str())); 
        if(FtPause<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="TMAX"){ 
        TimeMax=atof(txoptfull.c_str()); 
        if(TimeMax<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="TOUT"){ 
        TimePart=atof(txoptfull.c_str()); 
        if(TimePart<0)ErrorParm(opt,c,lv,file);
      }
      else if(txword=="TOUTX"){ 
        TimePartExtra=atof(txoptfull.c_str()); 
      }
      else if(txword=="DOMAIN_FIXED"){
        LoadDouble6(txoptfull,0,DomainFixedMin,DomainFixedMax);
        DomainMode=2;
      }
      else if(txword=="NSTEPS"){
        NstepsBreak=atoi(txoptfull.c_str()); 
        if(NstepsBreak)NoRtimes=true;
      }
      else if(txword=="SVSTEPS"){
        SvAllSteps=OptIsEnabled(txoptfull);
        if(SvAllSteps)NoRtimes=true;
      }
      else if(txword=="NORTIMES")NoRtimes=OptIsEnabled(txoptfull);
      else if(txword=="SVPIPS"){
        PipsMode=(unsigned)atoi(txopt1.c_str());
        if(PipsMode>2)ErrorParm(opt,c,lv,file);
        if(!txopt2.empty())PipsSteps=(unsigned)atoi(txopt2.c_str());
      }
#ifdef _WITHMR //<vs_vrres_ini>
      else if(txword=="VRES")VRes=(OptIsEnabled(txoptfull)? 1: 0);
      else if(txword=="VRES_FAST")MRFastSingle=(txoptfull!=""? atoi(txoptfull.c_str()): 1);
      else if(txword=="VRES_ORDER")VResOrder=(txoptfull!=""? atoi(txoptfull.c_str()): -1);
      else if(txword=="VRES_METHOD")VResMethod=(txoptfull!=""? atoi(txoptfull.c_str()): -1);
      else if(txword=="VRES_THRESHOLD")VResThreshold=float(txoptfull!=""?atof(txoptfull.c_str()): -1); 
      
#endif         //<vs_vrres_end>
      else if(txword=="OPT" && c+1<optn){ LoadFile(optlis[c+1],lv+1); c++; }
      else if(txword=="H" || txword=="HELP" || txword=="?")PrintInfo=true;
      else ErrorParm(opt,c,lv,file);
    }
  }
#ifdef _WITHMR //<vs_vrres_ini>
  if(!lv && VRes<0){
    VRes=(fun::FileExists(CaseName+"_vres00.bi4")? 1: 0);
  }
#endif         //<vs_vrres_end>
}

