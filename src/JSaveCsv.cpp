//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSaveCsv.cpp \brief Implements the class \ref JSaveCsv.

#include "JSaveCsv.h"
#include "Functions.h"

#include <cstring>
#include <cstdlib>
#include <stdarg.h>

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

using namespace std;

//##############################################################################
//# JSaveCsv
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSaveCsv::JSaveCsv(std::string fname,bool app,bool csvsepcoma):CsvSepComa(csvsepcoma){
  ClassName="JSaveCsv";
  CsvSep[0]=(csvsepcoma? ',': ';');
  CsvSep[1]='\0';
  Pf=NULL;
  Reset();
  FileName=fname;
  App=app;
}

//==============================================================================
/// Destructor.
//==============================================================================
JSaveCsv::~JSaveCsv(){ 
  if(!FileError)SaveData(true);
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JSaveCsv::Reset(){
  delete Pf; Pf=NULL;
  FileError=false;
  FileName="";
  App=false;
  Head="";
  Data="";
}

//==============================================================================
/// Write string to file checking for line breaks.
/// Graba string en fichero comprobando los saltos de linea.
//==============================================================================
void JSaveCsv::Save(const std::string &tx){
  Pf->write(tx.c_str(),tx.size());
  Pf->flush();
  //fflush(NULL);//-Vacia todos los bufers
}

//==============================================================================
/// Writes data in file.
/// Graba datos en fichero.
//==============================================================================
void JSaveCsv::SaveData(bool closefile){
  const char met[]="SaveData";
  if(Pf==NULL){
    Pf=new fstream();
    const bool fexists=fun::FileExists(FileName);
    if(App && fexists)Pf->open(FileName.c_str(),ios::binary|ios::out|ios::in|ios::app);
    else Pf->open(FileName.c_str(),ios::binary|ios::out);
    if(!(*Pf)){
      FileError=true;
      RunException(met,"File could not be opened.",FileName);
    }
    if(App && fexists)Pf->seekp(0,Pf->end);
    else{
      Save(Head);
      if(Pf->fail())RunException(met,"File writing failure.",FileName);
    }
  }
  if(!FileError){
    Save(Data);
    Data="";
    if(Pf->fail())RunException(met,"File writing failure.",FileName);
    if(closefile){
      Pf->close();
      delete Pf; Pf=NULL;
    }
  }
}

//==============================================================================
/// Adds heads.
/// Añade texto de cabecera.
//==============================================================================
void JSaveCsv::AddHead(const std::string &head,bool endl){
  string tx=(!Head.empty() && Head[Head.size()-1]!='\n'? ";": "");
  tx=tx+(endl? head+"\n": head);
  Head=Head+fun::StrCsvSep(CsvSepComa,tx);
}

//==============================================================================
/// Adds data.
/// Añade texto en detalle.
//==============================================================================
void JSaveCsv::AddData(const std::string &data,bool endl){
  Data=Data+(!Data.empty() && Data[Data.size()-1]!='\n'? CsvSep: "")+(endl? data+"\n": data);
}

//==============================================================================
/// Adds data using a format.
/// Añade valores en detalle con formato.
//==============================================================================
void JSaveCsv::AddValuesf(const char *format,...){
  const std::string format2=fun::StrCsvSep(CsvSepComa,format);
  const char *formatok=format2.c_str();
  const unsigned SIZE=1024;
  char buffer[SIZE+1];
  va_list args;
  va_start(args, format);
  unsigned size=vsnprintf(buffer,SIZE,formatok,args);
  if(size<SIZE)AddData(buffer,false);
  else{
    char *buff=new char[size+1];
    vsnprintf(buff,size,formatok,args);
    AddData(buff,false);
    delete[] buff;
  }
  va_end(args);
}

//==============================================================================
/// Adds value with string type.
/// Añade valor de tipo string en detalle.
//==============================================================================
void JSaveCsv::AddValue(const std::string &v){
  Data=Data+(!Data.empty() && Data[Data.size()-1]!='\n'? CsvSep: "")+fun::StrCsvSep(CsvSepComa,v);
}

//==============================================================================
/// Adds value with usnigned type.
/// Añade valor de tipo unsigned en detalle.
//==============================================================================
void JSaveCsv::AddValue(unsigned v){
  Data=Data+(!Data.empty() && Data[Data.size()-1]!='\n'? CsvSep: "")+fun::UintStr(v);
}

//==============================================================================
/// Adds value with int type.
/// Añade valor de tipo int en detalle.
//==============================================================================
void JSaveCsv::AddValue(int v){
  Data=Data+(!Data.empty() && Data[Data.size()-1]!='\n'? CsvSep: "")+fun::IntStr(v);
}

//==============================================================================
/// Adds value with float type.
/// Añade valor de tipo float en detalle.
//==============================================================================
void JSaveCsv::AddValue(float v){
  Data=Data+(!Data.empty() && Data[Data.size()-1]!='\n'? CsvSep: "")+fun::FloatStr(v,"%15.7E");
}

//==============================================================================
/// Adds value with double type.
/// Añade valor de tipo double en detalle.
//==============================================================================
void JSaveCsv::AddValue(double v){
  Data=Data+(!Data.empty() && Data[Data.size()-1]!='\n'? CsvSep: "")+fun::DoubleStr(v,"%20.12E");
}

//==============================================================================
/// Adds end of line.
/// Añade final de linea.
//==============================================================================
void JSaveCsv::AddEndl(){
  Data=Data+"\n";
}


