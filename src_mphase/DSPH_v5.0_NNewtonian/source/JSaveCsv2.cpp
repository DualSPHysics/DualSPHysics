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

/// \file JSaveCsv2.cpp \brief Implements the class \ref JSaveCsv2.

#include "JSaveCsv2.h"
#include "Functions.h"

#include <cstring>
#include <cstdlib>
#include <stdarg.h>

using namespace std;

namespace jcsv{

//##############################################################################
//# JSaveCsv2
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSaveCsv2::JSaveCsv2(std::string fname,bool app,bool csvsepcoma):ExceptionThrown(false),App(app),CsvSepComa(csvsepcoma){
  ClassName="JSaveCsv2";
  Pf=NULL;
  Reset();
  FileName=fname;
  OpenFile();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSaveCsv2::~JSaveCsv2(){ 
  DestructorActive=true;
  if(!ExceptionThrown)SaveData(true);
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JSaveCsv2::Reset(){
  delete Pf; Pf=NULL;
  AppendMode=false;
  FileName="";
  FirstSaveData=true;
  AutoSepEnable=true;
  InitFmt();
  SetData();
  Head="";
  HeadLineEmpty=true;
  Data="";
  DataLineEmpty=true;
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JSaveCsv2::OpenFile(){
  if(Pf==NULL){
    Pf=new fstream();
    const bool fexists=fun::FileExists(FileName);
    if(App && fexists)Pf->open(FileName.c_str(),ios::binary|ios::out|ios::in|ios::app);
    else Pf->open(FileName.c_str(),ios::binary|ios::out);
    if(!(*Pf)){
      ExceptionThrown=true;
      Run_ExceptioonFile("File could not be opened.",FileName);
    }
    if(App && fexists){
      AppendMode=true;
      Pf->seekp(0,Pf->end);
    }
  }
}

//==============================================================================
/// Initialization of output formats.
//==============================================================================
void JSaveCsv2::InitFmt(){
  FmtCurrent[TpSigned1]  =FmtDefault[TpSigned1]  ="%d";
  FmtCurrent[TpSigned2]  =FmtDefault[TpSigned2]  ="%d;%d";
  FmtCurrent[TpSigned3]  =FmtDefault[TpSigned3]  ="%d;%d;%d";
  FmtCurrent[TpSigned4]  =FmtDefault[TpSigned4]  ="%d;%d;%d;%d";
  FmtCurrent[TpUnsigned1]=FmtDefault[TpUnsigned1]="%u";
  FmtCurrent[TpUnsigned2]=FmtDefault[TpUnsigned2]="%u;%u";
  FmtCurrent[TpUnsigned3]=FmtDefault[TpUnsigned3]="%u;%u;%u";
  FmtCurrent[TpUnsigned4]=FmtDefault[TpUnsigned4]="%u;%u;%u;%u";
  FmtCurrent[TpLlong1]   =FmtDefault[TpLlong1 ]  ="%lld";
  FmtCurrent[TpUllong1]  =FmtDefault[TpUllong1]  ="%llu";
  FmtCurrent[TpFloat1]   =FmtDefault[TpFloat1]   ="%f";
  FmtCurrent[TpFloat2]   =FmtDefault[TpFloat2]   ="%f;%f";
  FmtCurrent[TpFloat3]   =FmtDefault[TpFloat3]   ="%f;%f;%f";
  FmtCurrent[TpFloat4]   =FmtDefault[TpFloat4]   ="%f;%f;%f;%f";
  FmtCurrent[TpDouble1]  =FmtDefault[TpDouble1]  ="%f";
  FmtCurrent[TpDouble2]  =FmtDefault[TpDouble2]  ="%f;%f";
  FmtCurrent[TpDouble3]  =FmtDefault[TpDouble3]  ="%f;%f;%f";
  FmtCurrent[TpDouble4]  =FmtDefault[TpDouble4]  ="%f;%f;%f;%f";
}
 
//==============================================================================
/// Adds string to data or head.
//==============================================================================
void JSaveCsv2::AddStr(std::string tx){
  if(AutoSepEnable)SetSeparators(tx);
  const bool jump=(!tx.empty() && tx[0]=='\n');
  if(DataSelected){
    if(!(Data.empty() || Data[Data.size()-1]=='\n' || jump))Data.append(";");
    Data.append(tx);
    //printf("Data:[%s]\n",Data.c_str());
  }
  else{
    if(!(Head.empty() || Head[Head.size()-1]=='\n' || jump))Head.append(";");
    Head.append(tx);
    //printf("Head:[%s]\n",Head.c_str());
  }
}
 
//==============================================================================
/// Adds one or several field separators.
//==============================================================================
void JSaveCsv2::AddSeparator(unsigned count){
  if(DataSelected)for(unsigned c=0;c<count;c++)Data=Data+";";
  else for(unsigned c=0;c<count;c++)Head=Head+";";
}
 
//==============================================================================
/// Adds end of line.
//==============================================================================
void JSaveCsv2::AddEndl(){
  if(DataSelected)Data=Data+"\n";
  else Head=Head+"\n";
}
 
//==============================================================================
/// Add header for triple data.
//==============================================================================
void JSaveCsv2::AddHead3(std::string headtx){
  if(!DataSelected){
    const string mark=(fun::StrSplitCount(";",headtx)>=fun::StrSplitCount(",",headtx)? ";": ",");
    vector<string> vec;
    const unsigned nv=fun::VectorSplitStr(mark,headtx,vec);
    for(unsigned cv=0;cv<nv;cv++)if(!vec[cv].empty()){
      string txend=vec[cv];
      string txini=fun::StrSplit(" ",txend);
      if(!txini.empty()){
        AddStr(txini+".x "+txend);
        AddStr(txini+".y "+txend);
        AddStr(txini+".z "+txend);
      }
    }
  }
}
  
//==============================================================================
/// Returns string using the same parameters used in printf().
//==============================================================================
std::string JSaveCsv2::ToStr(const char *format,...)const{
  std::string ret;
  const unsigned SIZE=1024;
  char buffer[SIZE+1];
  va_list args;
  va_start(args, format);
  int size=vsnprintf(buffer,SIZE,format,args);
  if(size>=0 && size<SIZE)ret=buffer;
  else{
    int rsize=-1;
    int size2=SIZE+SIZE*2;
    for(int c=0;c<10 && rsize<0;c++,size2+=SIZE*2){
      char *buff2=new char[size2+1];
      rsize=vsnprintf(buff2,size2,format,args);
      if(rsize>=0)ret=buff2;
      delete[] buff2;
    }
    if(rsize<0)Run_Exceptioon("Error in JSaveCsv2::ToStr(): Output text is too long.");
  }
  va_end(args);
  return(ret);
}

//==============================================================================
/// Operator: Adds one or several field separators.
//==============================================================================
JSaveCsv2& JSaveCsv2::operator <<(const Sep &obj){
  AddSeparator(obj.Count);
  return(*this);
}

//==============================================================================
/// Operator: Enable auto separator change according configuration.
//==============================================================================
JSaveCsv2& JSaveCsv2::operator <<(const AutoSepOn &obj){
  AutoSepEnable=true;
  return(*this);
}

//==============================================================================
/// Operator: Disable auto separator change according configuration.
//==============================================================================
JSaveCsv2& JSaveCsv2::operator <<(const AutoSepOff &obj){
  AutoSepEnable=false;
  return(*this);
}

//==============================================================================
/// Operator: Adds end of line.
//==============================================================================
JSaveCsv2& JSaveCsv2::operator <<(const Endl &obj){
  AddEndl();
  return(*this);
}

//==============================================================================
/// Operator: Sets format for output.
//==============================================================================
JSaveCsv2& JSaveCsv2::operator <<(const Fmt &obj){
  FmtCurrent[obj.TypeFmt]=(obj.Format.empty()? FmtDefault[obj.TypeFmt]: obj.Format);
  return(*this);
}

//==============================================================================
/// Operator: Add header for triple data.
//==============================================================================
JSaveCsv2& JSaveCsv2::operator <<(const Head3 &obj){
  AddHead3(obj.HeadText);
  return(*this);
}

//==============================================================================
/// Write string to file checking for line breaks.
/// Graba string en fichero comprobando los saltos de linea.
//==============================================================================
void JSaveCsv2::Save(const std::string &tx){
  Pf->write(tx.c_str(),tx.size());
  Pf->flush();
  //fflush(NULL);//-Vacia todos los bufers
  if(Pf->fail())Run_ExceptioonFile("File writing failure.",FileName);
}

//==============================================================================
/// Sets separators according configuration.
/// Cambia separadores segun configuracion.
//==============================================================================
void JSaveCsv2::SetSeparators(std::string &tx)const{
  const char sep0=(CsvSepComa? ';': ',');
  const char sep1=(CsvSepComa? ',': ';');
  const unsigned size=unsigned(tx.size());
  for(unsigned c=0;c<size;c++)if(tx[c]==sep0)tx[c]=sep1;
}

//==============================================================================
/// Writes data in file.
/// Graba datos en fichero.
//==============================================================================
void JSaveCsv2::SaveData(bool closefile){
  if(FirstSaveData || !Data.empty()){
    if(Pf==NULL)OpenFile();
    if(FirstSaveData && !AppendMode)Save(Head);
    Save(Data);
    Data="";
    FirstSaveData=false;
  }
  if(Pf!=NULL && closefile){
    Pf->close();
    delete Pf; Pf=NULL;
  }
}

}

