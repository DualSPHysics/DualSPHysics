/*
 <DUALSPHYSICS>  Copyright (c) 2018, 
 Dr Jose M. Dominguez, Dr Alejandro Crespo, Prof. Moncho Gomez Gesteira, Orlando Garcia-Feal, Dr Anxo Barreiro,
 Dr Benedict Rogers, Dr Georgios Fourtakas, Dr Athanasios Mokos, Prof. Peter Stansby,
 Dr Renato Vacondio, Dr Ricardo Canelas, Dr Stephen Longshaw, Dr Corrado Altomare.

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
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
\version 4.2.050
\date 21-05-2018
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
#include "JCfgRun.h"
#include "JException.h"
#include "JSphCpuSingle.h"
#ifdef _WITHGPU
  #include "JSphGpuSingle.h"
#endif

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

using namespace std;


JAppInfo AppInfo("DualSPHysics4","v4.2.112","21-12-2018");
//JAppInfo AppInfo("DualSPHysics4","v4.2.???","UserVersion","v1.0","??-??-????"); //-for user versions.

//==============================================================================
/// LGPL License.
//==============================================================================
std::string getlicense_lgpl(const std::string &name,bool simple){
  std::string tx=(simple? "": "\n");
  tx=tx+"\n <"+fun::StrUpper(name)+"> Copyright (C) 2018 by"; 
  tx=tx+"\n Dr Jose M. Dominguez, Dr Alejandro Crespo,";
  tx=tx+"\n Prof. Moncho Gomez Gesteira, Orlando Garcia-Feal, Dr Anxo Barreiro";
  tx=tx+"\n Dr Benedict Rogers, Dr Georgios Fourtakas, Dr Athanasios Mokos,";
  tx=tx+"\n Prof. Peter Stansby, Dr Renato Vacondio, Dr Ricardo Canelas,";
  tx=tx+"\n Dr Stephen Longshaw, Dr Corrado Altomare.\n";
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
    tx=tx+"\n along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.\n\n";
  }
  return(tx);
}

//==============================================================================
///  Shows program version and finishes the execution.
//==============================================================================
bool ShowsVersion(int argc,char** argv){
  string ver=(argc==2? argv[1]: "");
  bool ret=(ver=="-ver" || ver=="-VER");
  if(ret){
    printf("%s\n",AppInfo.GetFullName().c_str());
    printf("%s",getlicense_lgpl(AppInfo.GetShortName(),true).c_str());
  }
  return(ret);
}

//==============================================================================
//==============================================================================
int main(int argc, char** argv){
  int errcode=1;
#ifdef CODE_SIZE4
  AppInfo.AddNameExtra("MK65k");
#endif
  AppInfo.ConfigRunPaths(argv[0]);
  if(ShowsVersion(argc,argv))return(errcode);
  std::string license=getlicense_lgpl(AppInfo.GetShortName(),false);
  printf("%s",license.c_str());
  std::string appname=AppInfo.GetFullName();
  std::string appnamesub;
  for(unsigned c=0;c<=unsigned(appname.size());c++)appnamesub=appnamesub+"=";
  printf("\n%s\n%s\n",appname.c_str(),appnamesub.c_str());
  JLog2 *log=NULL;
  JCfgRun cfg;
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
        JSphCpuSingle sph;
        sph.Run(appname,&cfg,log);
      }
      #ifdef _WITHGPU
      else{
        JSphGpuSingle sph;
        sph.Run(appname,&cfg,log);
      }
      #endif
    }
    errcode=0;
  }
  catch(const char *cad){
    string tx=string("\n*** Exception: ")+cad+"\n";
    if(log && log->IsOk())log->Print(tx); else printf("%s",tx.c_str());
  }
  catch(const string &e){
    string tx=string("\n*** Exception: ")+e+"\n";
    if(log && log->IsOk())log->Print(tx); else printf("%s",tx.c_str());
  }
  catch (const exception &e){
    string tx=string("\n*** ")+e.what()+"\n";
    if(log && log->IsOk())log->Print(tx); else printf("%s",tx.c_str());
  }
  catch(...){
    printf("\n*** Attention: Unknown exception...\n");
  }
  return(errcode);
}


