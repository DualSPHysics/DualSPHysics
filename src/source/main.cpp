/*
 <DUALSPHYSICS>  Copyright (c) 2025, 
 Dr Jose M. Dominguez Alonso, Prof. Alejandro Crespo, 
 Dr Georgios Fourtakas, Prof. Benedict Rogers, 
 Dr Renato Vacondio, Dr Corrado Altomare, 
 Dr Angelo Tafuni, Dr Orlando Garcia Feal, 
 Dr Ivan Martinez Estevez, Dr Joseph O'Connor, 
 Dr Aaron English, Dr Francesco Ricci,
 Prof. Moncho Gomez Gesteira, Prof. Peter Stansby

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 

 Free download source files and documentation from DualSPHysics website (http://dual.sphysics.org).
 Full list of developers and contributors is available at https://dual.sphysics.org/developers
 Properly cite the use of DualSPHysics (Dominguez et al., 2022) (https://dual.sphysics.org/references).
*/

/** \mainpage DualSPHysics Documentation
\section main-des Description
DualSPHysics is based on the Smoothed Particle Hydrodynamics <a href="http://www.sphysics.org">SPHysics code.</a> \n
The package is a set of C++ and CUDA codes. \n
DualSPHysics is developed to deal with real-life engineering problems <a href="http://www.youtube.com/user/DualSPHysics">DualSPHysics animations.</a> \n

EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain. \n
School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.
\section compile_sec Project files
Please download source files and documentation from <a href="http://dual.sphysics.org">DualSPHysics website.</a> \n
\author <a href="http://dual.sphysics.org/index.php/developers">DualSPHysics Developers.</a> 
\version 5.4.355
\date 07-04-2025
\copyright GNU Lesser General Public License <a href="http://www.gnu.org/licenses/">GNU licenses.</a>
*/

/// \file main.cpp \brief Main file of the project that executes the code on CPU or GPU.

#include <string>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include "JAppInfo.h"
#include "JLog2.h"
#include "JException.h"
#include "JSphCfgRun.h"
#include "JSphCpuSingle.h"
#ifdef _WITHGPU
  #include "JSphGpuSingle.h"
#endif
#ifdef _WITHMR
  #include "JSphVResDriver.h"
  #include "JSphCpuSingle_VRes.h"
  #ifdef _WITHGPU
    #include "JSphGpuSingle_VRes.h"
  #endif
#endif

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

using namespace std;

JAppInfo AppInfo("DualSPHysics5","v5.4.355","08-04-2025");
//JAppInfo AppInfo("DualSPHysics5","v5.4.353","ShiftAdv","v0.32","11-12-2024"); //-for user versions.
//JAppInfo AppInfo("DualSPHysics5","v5.0.???","UserVersion","v1.0","??-??-????"); //-for user versions.

//==============================================================================
/// LGPL License.
//==============================================================================
std::string getlicense_lgpl(const std::string& name,bool simple){
  std::string tx=(simple? "": "\n");
  tx=tx+"\n <"+fun::StrUpper(name)+"> Copyright (c) 2025 by"; 
  tx=tx+"\n Dr Jose M. Dominguez Alonso, Prof. Alejandro Crespo,";
  tx=tx+"\n Dr Georgios Fourtakas, Prof. Benedict Rogers,";
  tx=tx+"\n Dr Renato Vacondio, Dr Corrado Altomare,";
  tx=tx+"\n Dr Angelo Tafuni, Dr Orlando Garcia Feal,";
  tx=tx+"\n Dr Ivan Martinez Estevez, Dr Joseph O'Connor,";
  tx=tx+"\n Dr Aaron English, Dr Francesco Ricci,";
  tx=tx+"\n Prof. Moncho Gomez Gesteira, Prof. Peter Stansby\n";
  if(!simple){
    tx=tx+"\n EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo";
    tx=tx+"\n School of Mechanical, Aerospace and Civil Engineering, University of Manchester\n";
    tx=tx+"\n DualSPHysics is free software: you can redistribute it and/or"; 
    tx=tx+"\n modify it under the terms of the GNU Lesser General Public License";
    tx=tx+"\n as published by the Free Software Foundation, either version 2.1 of"; 
    tx=tx+"\n the License, or (at your option) any later version.\n"; 
    tx=tx+"\n DualSPHysics is distributed in the hope that it will be useful,"; 
    tx=tx+"\n but WITHOUT ANY WARRANTY; without even the implied warranty of"; 
    tx=tx+"\n MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"; 
    tx=tx+"\n GNU Lesser General Public License for more details.\n";
    tx=tx+"\n You should have received a copy of the GNU Lesser General Public License"; 
    tx=tx+"\n along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.\n";
    tx=tx+"\n Free download source files and documentation from DualSPHysics website (http://dual.sphysics.org).";
    tx=tx+"\n Full list of developers and contributors is available at https://dual.sphysics.org/developers";
    tx=tx+"\n Properly cite the use of DualSPHysics (Dominguez et al., 2022) (https://dual.sphysics.org/references).\n\n";
  }
  return(tx);
}

//==============================================================================
///  Shows program version and JSON information and finishes the execution.
//==============================================================================
bool ShowsVersionInfo(int argc,char** argv,bool show=true){
  const string option=fun::StrLower(argc==2? argv[1]: "");
  bool finish=true;
  if(fun::StrRemoveAfter(option,":")=="-ver"){
    const string vtex=JCfgRunBase::VerText(AppInfo.GetFullName(),option);
    if(show){
      printf("%s\n",vtex.c_str());
      if(vtex==AppInfo.GetFullName())printf("%s",getlicense_lgpl(AppInfo.GetShortName(),true).c_str());
    }
  }
  else if(option=="-info"){
    //-Defines the features included in the program.
    std::vector<std::string> features;
    features.push_back(fun::JSONProperty("CPU",true));
    features.push_back(fun::JSONProperty("GPU",AVAILABLE_GPU));
    features.push_back(fun::JSONProperty("MultiGPU",AVAILABLE_MGPU));
    features.push_back(fun::JSONProperty("CHRONO_Coupling",AVAILABLE_CHRONO));
    features.push_back(fun::JSONProperty("MoorDynPlus_Coupling",AVAILABLE_MOORDYNPLUS));
    features.push_back(fun::JSONProperty("WaveGen",AVAILABLE_WAVEGEN));
    features.push_back(fun::JSONProperty("FlexStruc",AVAILABLE_FLEXSTRUC)); //<vs_flexstruc>
    features.push_back(fun::JSONProperty("DDT_Fourtakas",true));
    //-Defines main information about the version program.
    std::vector<std::string> info;
    info.push_back(fun::JSONProperty("ShortName",AppInfo.GetShortName()));
    info.push_back(fun::JSONProperty("FullName" ,AppInfo.GetFullName()));
    info.push_back(fun::JSONProperty("Version"  ,AppInfo.GetMainVer()));
    info.push_back(fun::JSONProperty("Date"     ,AppInfo.GetDate()));
    info.push_back(fun::JSONPropertyValue("Features",fun::JSONObject(features)));
    if(show)printf("%s\n",fun::JSONObject(info).c_str());
  }
  else finish=false;
  return(finish);
}

//==============================================================================
/// Print exception message on standard error output and log file.
//==============================================================================
void PrintExceptionLog(const std::string& prefix,const std::string& text
  ,JLog2* log)
{
  const bool prt=(text.empty() || text[0]!='#');
  const string tx=(prt? prefix+text: text.substr(1));
  //if(prt)printf("%s\n",tx.c_str());
  fflush(stdout);
  if(prt)cerr << tx << endl;
  if(log && log->IsOk())log->PrintFile(tx,true);
}

//==============================================================================
//==============================================================================
int main(int argc, char** argv){
  int errcode=1;

  //AppInfo.AddNameExtra("SaveFtAce");
  #ifdef CODE_SIZE4
    AppInfo.AddNameExtra("MK65k");
  #endif

  AppInfo.ConfigRunPaths(argv[0]);
  const std::string appname=AppInfo.GetFullName();
  const std::string license=getlicense_lgpl(AppInfo.GetShortName(),false);
  if(ShowsVersionInfo(argc,argv))return(errcode);
  printf("%s",license.c_str());
  const std::string appnamesub=fun::StrFillEnd("","=",unsigned(appname.size())+1);
  printf("\n%s\n%s\n",appname.c_str(),appnamesub.c_str());
  JLog2* log=NULL;
  JSphCfgRun cfg;
  cfg.SetFeatureList(JSph::GetFeatureList());
  try{
    cfg.LoadArgv(argc,argv);
    //cfg.VisuConfig();
    if(!cfg.PrintInfo){
      AppInfo.ConfigOutput(cfg.CreateDirs,cfg.CsvSepComa,cfg.DirOut,cfg.DirDataOut);
      AppInfo.LogInit(AppInfo.GetDirOut()+"/Run.out");
      log=AppInfo.LogPtr();
      log->AddFileInfo(cfg.DirOut+"/Run.out","Log file of the simulation.");
      log->Print(license,JLog2::Out_File);
      log->Print(appname,JLog2::Out_File);
      log->Print(appnamesub,JLog2::Out_File);
      //-SPH Execution.
      #ifndef _WITHGPU
        cfg.Cpu=true;
      #endif
      if(cfg.Cpu){
        #if defined(_WITHMR)               //<vs_vrres_ini>
        if(cfg.VRes){
          JSphVResDriver<JSphCpuSingle_VRes, stinterparmscb> sph;
          sph.Run(appname,&cfg,log);
        }else{
          JSphCpuSingle sph;
          sph.Run(appname,&cfg,log);        //<vs_vrres_end>
        }
        #else
          JSphCpuSingle sph;
          sph.Run(appname,&cfg,log);
        #endif
      }else{
      #if defined(_WITHGPU) && defined(_WITHMR)               //<vs_vrres_ini>
      if(cfg.VRes){                             
        JSphVResDriver<JSphGpuSingle_VRes, StInterParmsbg> sph;
        sph.Run(appname,&cfg,log);
      }else{
        JSphGpuSingle sph;
        sph.Run(appname,&cfg,log);
      }
      #elif defined(_WITHGPU)                                 //<vs_vrres_end>
      else{
        JSphGpuSingle sph;
        sph.Run(appname,&cfg,log);
      }
      #endif
      }
    }
    errcode=0;
  }
  catch(const char* cad){
    PrintExceptionLog("\n*** Exception(chr): ",cad,log);
  }
  catch(const string& e){
    PrintExceptionLog("\n*** Exception(str): ",e,log);
  }
  catch (const JException& e){
    if(log && log->IsOk())log->PrintFile(e.what());
  }
  catch (const exception& e){
    PrintExceptionLog("\n*** Exception(exc): ",e.what(),log);
  }
  catch(...){
    PrintExceptionLog("","\n*** Attention: Unknown exception...",log);
  }
  //-Finished execution.
  if(log && log->IsOk())log->PrintFile(fun::PrintStr("\nFinished execution (code=%d).\n",errcode),true);
  printf("\nFinished execution (code=%d).\n",errcode);
  fflush(stdout);
  return(errcode);
}


