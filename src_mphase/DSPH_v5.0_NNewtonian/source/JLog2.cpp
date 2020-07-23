//HEAD_DSCODES
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

/// \file JLog2.cpp \brief Implements the class \ref JLog2.

#include "JLog2.h"
#include "Functions.h"
#include <stdarg.h>
#include <algorithm>

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

using std::string;
using std::ofstream;
using std::endl;

//##############################################################################
//# JLog2
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JLog2::JLog2(TpMode_Out modeoutdef):ModeOutDef(modeoutdef){
  ClassName="JLog2";
  Pf=NULL;
  Reset();
}

//==============================================================================
/// Constructor.
//==============================================================================
JLog2::JLog2(JLog2 *parent,std::string prefix):ModeOutDef(parent->ModeOutDef){
  ClassName="JLog2";
  Pf=NULL;
  Reset();
  Parent=parent;
  ParentPrefix=prefix;
}

//==============================================================================
/// Destructor.
//==============================================================================
JLog2::~JLog2(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JLog2::Reset(){
  Parent=NULL;
  ParentPrefix="";
  FileName="";
  Ok=false;
  MpiRun=false;
  MpiRank=0; MpiLaunch=0;
  if(Pf){
    if(Pf->is_open())Pf->close();
    delete Pf; Pf=NULL;
  }
  Warnings.clear();
  FileInfo.clear();
  DirOut="";
}

//==============================================================================
/// Initialisation of log file.
//==============================================================================
void JLog2::Init(std::string fname,bool mpirun,int mpirank,int mpilaunch){
  if(Parent)Run_Exceptioon("Method is not available for dependent log.");
  Reset();
  MpiRun=mpirun; MpiRank=mpirank; MpiLaunch=mpilaunch;
  if(!fname.empty()){//-When file for log is defined.
    if(MpiRun){
      string ext=fun::GetExtension(fname);
      fname=fun::GetWithoutExtension(fname);
      fname=fname+"_"+fun::IntStrFill(mpirank,MpiLaunch-1);
      if(!ext.empty())fname=fname+"."+ext;
    }
    FileName=fname;
    DirOut=fun::GetDirWithSlash(fun::GetDirParent(FileName));
    if(ModeOutDef&Out_File){
      Pf=new ofstream; 
      Pf->open(fname.c_str());
      if(Pf)Ok=true;
      else Run_ExceptioonFile("Cannot open the file.",fname);
    }
  }
}
  
//==============================================================================
/// Visualises and/or stores information of the execution.
//==============================================================================
void JLog2::Print(const std::string &tx,TpMode_Out mode,bool flush){
  if(Parent){ Parent->Print(ParentPrefix+tx,mode,flush); return; }
  if(mode==Out_Default)mode=ModeOutDef;
  if(mode&Out_Screen){
    if(MpiRun){
      int pos=0;
      for(;pos<int(tx.length())&&tx[pos]=='\n';pos++)printf("%d>\n",MpiRank);
      printf("%d>%s\n",MpiRank,tx.substr(pos).c_str());
    }
    else printf("%s\n",tx.c_str());
  }
  if((mode&Out_File) && Pf)(*Pf) << tx << endl;
  if(flush)fflush(stdout);
}
  
//==============================================================================
/// Visualises and/or stores information of the execution.
//==============================================================================
void JLog2::Print(const std::vector<std::string> &lines,TpMode_Out mode,bool flush){
  for(unsigned c=0;c<unsigned(lines.size());c++)Print(lines[c],mode,false);
  if(flush)fflush(stdout);
}
  
//==============================================================================
/// Visualises and/or stores information of the execution.
//==============================================================================
void JLog2::Printf(const char *format,...){
  const int SIZE=1024;
  char buffer[SIZE+1];
  va_list args;
  va_start(args,format);
  int size=vsnprintf(buffer,SIZE,format,args);
  if(size>=0 && size<SIZE)Print(buffer);
  else{
    int rsize=-1;
    int size2=SIZE+SIZE*2;
    for(int c=0;c<10 && rsize<0;c++,size2+=SIZE*2){
      char *buff2=new char[size2+1];
      rsize=vsnprintf(buff2,size2,format,args);
      if(rsize>=0)Print(buff2);
      delete[] buff2;
    }
    if(rsize<0)Print("[***ERROR: Text is too long***]",Out_Default,true);
  }
  va_end(args);
}
  
//==============================================================================
/// Visualises and/or stores information of the execution.
//==============================================================================
void JLog2::PrintfDbg(const char *format,...){
  const int SIZE=1024;
  char buffer[SIZE+1];
  va_list args;
  va_start(args,format);
  int size=vsnprintf(buffer,SIZE,format,args);
  if(size>=0 && size<SIZE)Print(buffer,Out_Default,true);
  else{
    int rsize=-1;
    int size2=SIZE+SIZE*2;
    for(int c=0;c<10 && rsize<0;c++,size2+=SIZE*2){
      char *buff2=new char[size2+1];
      rsize=vsnprintf(buff2,size2,format,args);
      if(rsize>=0)Print(buff2,Out_Default,true);
      delete[] buff2;
    }
    if(rsize<0)Print("[***ERROR: Text is too long***]",Out_Default,true);
  }
  va_end(args);
}
  
//==============================================================================
/// Visualises and/or stores information of the execution adding a prefix.
//==============================================================================
void JLog2::Printp(const std::string &prefix,const std::vector<std::string> &lines,JLog2::TpMode_Out mode,bool flush){
  for(unsigned c=0;c<unsigned(lines.size());c++)Printp(prefix,lines[c],mode,false);
  if(flush)fflush(stdout);
}
  
//==============================================================================
/// Visualises and/or stores information of the execution adding a prefix.
//==============================================================================
void JLog2::Printfp(const std::string &prefix,const char *format,...){
  const int SIZE=1024;
  char buffer[SIZE+1];
  va_list args;
  va_start(args,format);
  int size=vsnprintf(buffer,SIZE,format,args);
  if(size>=0 && size<SIZE)Printp(prefix,buffer);
  else{
    int rsize=-1;
    int size2=SIZE+SIZE*2;
    for(int c=0;c<10 && rsize<0;c++,size2+=SIZE*2){
      char *buff2=new char[size2+1];
      rsize=vsnprintf(buff2,size2,format,args);
      if(rsize>=0)Printp(prefix,buff2);
      delete[] buff2;
    }
    if(rsize<0)Printp(prefix,"[***ERROR: Text is too long***]",JLog2::Out_Default,true);
  }
  va_end(args);
}
  
//==============================================================================
/// Visualises and/or stores information of the execution adding a prefix.
//==============================================================================
void JLog2::PrintfpDbg(const std::string &prefix,const char *format,...){
  const int SIZE=1024;
  char buffer[SIZE+1];
  va_list args;
  va_start(args,format);
  int size=vsnprintf(buffer,SIZE,format,args);
  if(size>=0 && size<SIZE)Printp(prefix,buffer,JLog2::Out_Default,true);
  else{
    int rsize=-1;
    int size2=SIZE+SIZE*2;
    for(int c=0;c<10 && rsize<0;c++,size2+=SIZE*2){
      char *buff2=new char[size2+1];
      rsize=vsnprintf(buff2,size2,format,args);
      if(rsize>=0)Printp(prefix,buff2,JLog2::Out_Default,true);
      delete[] buff2;
    }
    if(rsize<0)Printp(prefix,"[***ERROR: Text is too long***]",JLog2::Out_Default,true);
  }
  va_end(args);
}
  
//==============================================================================
/// Adds warning to warning list.
//==============================================================================
void JLog2::AddWarning(const std::string &tx){
  if(Parent){ Parent->AddWarning(tx); return; }
  Warnings.push_back(tx);
}
  
//==============================================================================
/// Visualises and stores warning.
//==============================================================================
void JLog2::PrintWarning(const std::string &tx,TpMode_Out mode,bool flush){
  AddWarning(tx);
  Print(string("\n*** WARNING: ")+tx+"\n",mode,flush);
}
  
//==============================================================================
/// Visualises and stores warning.
//==============================================================================
void JLog2::PrintfWarning(const char *format,...){
  const int SIZE=1024;
  char buffer[SIZE+1];
  va_list args;
  va_start(args,format);
  int size=vsnprintf(buffer,SIZE,format,args);
  if(size>=0 && size<SIZE)PrintWarning(buffer);
  else{
    int rsize=-1;
    int size2=SIZE+SIZE*2;
    for(int c=0;c<10 && rsize<0;c++,size2+=SIZE*2){
      char *buff2=new char[size2+1];
      rsize=vsnprintf(buff2,size2,format,args);
      if(rsize>=0)PrintWarning(buff2);
      delete[] buff2;
    }
    if(rsize<0)Print("[***ERROR: Text is too long***]",Out_Default,true);
  }
  va_end(args);
}
  
//==============================================================================
/// Visualises list of warnings.
//==============================================================================
void JLog2::PrintWarningList(const std::string &txhead,const std::string &txfoot,TpMode_Out mode,bool flush){
  if(Parent){ Parent->PrintWarningList(txhead,txfoot,mode,flush); return; }
  const unsigned nw=WarningCount();
  //Print(fun::PrintStr("[WARNINGS #:%u]",nw),mode,flush);
  if(!txhead.empty())Print(txhead,mode,flush);
  string fmt=fun::PrintStr("%%0%dd. ",std::max(1u,unsigned(fun::UintStr(nw).size())));
  for(unsigned c=0;c<nw;c++){
    string pref=fun::PrintStr(fmt.c_str(),c+1);
    Print(pref+Warnings[c],mode,flush);
  }
  if(!txfoot.empty())Print(txfoot,mode,flush);
}
  
//==============================================================================
/// Visualises list of warnings.
//==============================================================================
void JLog2::PrintWarningList(TpMode_Out mode,bool flush){
  const unsigned nw=WarningCount();
  //if(nw)PrintWarningList(fun::PrintStr("[WARNINGS #:%u]",nw)," ",mode,flush);
  if(nw)PrintWarningList("[WARNINGS]"," ",mode,flush);
}

//==============================================================================
/// Adds file description.
//==============================================================================
void JLog2::AddFileInfo(std::string fname,const std::string &finfo){
  if(Parent){ Parent->AddFileInfo(fname,finfo); return; }
  //-Removes execution path.
  if(int(fname.find(DirOut))==0)fname=fname.substr(DirOut.size());
  //-Checks if it already exists.
  const unsigned nf=FilesCount();
  unsigned cf=0;
  for(;cf<nf && FileInfo[cf].file!=fname;cf++);
  //-Adds new file description.
  if(cf==nf)FileInfo.push_back(StrFileInfo(fname,finfo));
}

//==============================================================================
/// Sorts file descriptions.
//==============================================================================
bool SortFilesList(JLog2::StFileInfo a,JLog2::StFileInfo b){
  return(a.file<b.file); 
}

//==============================================================================
/// Visualises list of file descriptions.
//==============================================================================
void JLog2::PrintFilesList(const std::string &txhead,const std::string &txfoot,TpMode_Out mode,bool flush){
  if(Parent){ Parent->PrintFilesList(txhead,txfoot,mode,flush); return; }
  const unsigned nf=FilesCount();
  if(!txhead.empty())Print(txhead,mode,flush);
//  string fmt=fun::PrintStr("%%0%dd. ",std::max(1u,unsigned(fun::UintStr(nw).size())));
  //-Compute maximum size of filenames.
  unsigned sfile=0;
  for(unsigned c=0;c<nf;c++)sfile=std::max(sfile,unsigned(FileInfo[c].file.size()));
  //-Sorts refined list.
  std::sort(FileInfo.begin(),FileInfo.end(),SortFilesList);
  //-Prints file list.
  for(unsigned cf=0;cf<nf;cf++){
    const string file=FileInfo[cf].file;
    Print(string("- ")+file+string(sfile-unsigned(file.size()),'.')+": "+FileInfo[cf].info,mode,flush);
  }
  if(!txfoot.empty())Print(txfoot,mode,flush);
}
  
//==============================================================================
/// Visualises list of file descriptions.
//==============================================================================
void JLog2::PrintFilesList(TpMode_Out mode,bool flush){
  const unsigned nf=FilesCount();
  if(nf)PrintFilesList("[Output files]"," ",mode,flush);
}


