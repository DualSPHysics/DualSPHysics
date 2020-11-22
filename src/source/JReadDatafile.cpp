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

/// \file JReadDatafile.cpp \brief Implements the class \ref JReadDatafile.

#include "JReadDatafile.h"
#include "Functions.h"
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JReadDatafile::JReadDatafile(){
  ClassName="JReadDatafile";
  Data=NULL; LineBegin=NULL;
  Reset();
}
//==============================================================================
/// Destructor.
//==============================================================================
JReadDatafile::~JReadDatafile(){
  DestructorActive=true;
  Reset();
}
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JReadDatafile::Reset(){
  delete[] Data; Data=NULL;
  delete[] LineBegin; LineBegin=NULL;
  SizeFile=Size=0;
  File=""; Sep=";";
  LineCount=RemLineCount=0;
  ResetReadLine();
}

//==============================================================================
/// Load file data.
//==============================================================================
void JReadDatafile::LoadFile(const std::string &file,unsigned maxsize){
  Reset();
  File=file;
  ifstream pf;
  pf.open(file.c_str(),ios::binary);
  if(pf){
    pf.seekg(0,ios::end);
    SizeFile=unsigned(pf.tellg())+1;
    if(SizeFile>maxsize)Run_ExceptioonFile(fun::PrintStr("File exceeds the maximum size allowed (%u bytes).",maxsize),file);
    Data=new char[SizeFile];
    pf.seekg(0,ios::beg);
    pf.read(Data,SizeFile-1);
    Data[SizeFile-1]='\n';
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",file);
  Size=SizeFile;
  //-Prepares to first use.
  ProcessLines();
  SetReadLine(0);
}

//==============================================================================
/// Replaces several spaces and tabulations by a tabulation.
/// Removes spaces and tabulations at begin and end of line.
//==============================================================================
void JReadDatafile::ProcessSpaces(){
  unsigned nsp=0;
  unsigned c2=0;
  bool begin=true;
  for(unsigned c=0;c<Size;c++){
    char let=Data[c];
    if(let==' ' || let=='\t')nsp++;
    else{
      if(let=='\n'){ nsp=0; begin=true; }
      if(nsp){ 
        if(!begin){ Data[c2]='\t'; c2++; }
        nsp=0; 
      }
      if(let!='\n')begin=false;
      if(c!=c2)Data[c2]=let;
      c2++;
    }
  }
  //if(nsp){ Data[c2]='\t'; c2++; nsp=0; }
  //printf("++> Remove spaces+tabs: %u -> %u\n",Size,c2);
  Size=c2;
}

//==============================================================================
/// Load begin and end of each line and counts lines.
//==============================================================================
void JReadDatafile::ProcessLines(){
  //-Remove character \r and counts lines.
  LineCount=0;
  unsigned c2=0;
  for(unsigned c=0;c<Size;c++){
    char let=Data[c];
    if(let=='\n')LineCount++;
    if(c!=c2)Data[c2]=let;
    if(let!='\r')c2++;
  }
  //printf("++> Remove \\r: %u -> %u\n",Size,c2);
  Size=c2;
  //-Remove empty lines at end.
  while(Size>1 && Data[Size-2]=='\n'){ Size--; LineCount--; }
  //printf("++> Remove empty lines: %u\n",c2-Size);
  //-Allocate memory.
  delete[] LineBegin; LineBegin=NULL;
  LineBegin=new unsigned[LineCount+1];
  //-Prepares lines and looks for separator.
  bool run=true;
  while(run){
    run=false;
    //-Load LineBegin[].
    unsigned lin=0;
    LineBegin[0]=0;
    for(unsigned c=0;c<Size;c++){
      char let=Data[c];
      if(let=='\n'){
        lin++; LineBegin[lin]=c+1;
      }
    }
    if(lin!=LineCount)Run_Exceptioon("Error counting lines.");
    LineBegin[LineCount]=Size;
    //-Counts remark lines.
    RemLineCount=0;
    for(int c=0;c<LineCount;c++)if(Data[LineBegin[c]]=='#')RemLineCount++;
    //printf("++> RemLineCount: %u\n",RemLineCount);
    //-Determines the separator.
    {
      unsigned sep0=0,sep1=0,sep2=0,sep3=0;
      unsigned nlin=20;
      for(int c=0;c<LineCount && nlin;c++)if(Data[LineBegin[c]]!='#'){
        nlin--;
        const unsigned pini=LineBegin[c];
        const unsigned pfin=LineBegin[c+1];
        for(unsigned p=pini;p<pfin;p++){
          switch(Data[p]){
            case ' ':   sep0++;  break;
            case '\t':  sep1++;  break;
            case ';':   sep2++;  break;
            case ',':   sep3++;  break;
          }
        }
      }
      //printf("++> Sep: [ ]:%u [\\t]:%u [,]:%u [;]:%u\n",sep0,sep1,sep2,sep3);
      sep1+=sep0;
      if(sep1>=sep2 && sep1>=sep3)Sep="\t";
      else if(sep2>=sep1 && sep2>=sep3)Sep=";";
      else Sep=",";
      if(Sep=="\t" && sep0){
        ProcessSpaces();
        run=true;
      }
    }
  }
  //printf("++> LineCount:%u(%u)  Size:%u -> %u\n",LineCount,RemLineCount,SizeFile,Size);
}


//==============================================================================
/// Removes character.
//==============================================================================
void JReadDatafile::RemoveChar(char let){
  unsigned c2=0;
  for(unsigned c=0;c<Size;c++){
    const char v=Data[c];
    if(v!=let){
      if(c!=c2)Data[c2]=v;
      c2++;
    }
  }
  //printf("++> Remove char[%u:%c]: %u -> %u\n",let,let,Size,c2);
  Size=c2;
  //-Prepares to first use.
  ProcessLines();
  SetReadLine(0);
}


//==============================================================================
/// Restarts reading position.
//==============================================================================
void JReadDatafile::ResetReadLine(){ 
  ReadLin=ReadLinValue=-1;
  ReadLine="";
  ReadValue="";
}

//==============================================================================
/// Configures reading position. (line=0...n-1)
//==============================================================================
void JReadDatafile::SetReadLine(int line){ 
  ResetReadLine();
  ReadLin=line;
  ReadLine=GetLine(ReadLin);
  ReadLinValue=-1;
}

//==============================================================================
/// Returns selected line.
//==============================================================================
std::string JReadDatafile::GetLine(int line)const{
  if(line<0 || line>=LineCount)Run_ExceptioonFile("Line number is invalid.",File);
  const unsigned ini=LineBegin[line];
  const unsigned len=LineBegin[line+1]-ini-1;//-Con -1 se quita el salto de linea.
  string tex;
  if(len){
    tex.resize(len);
    memcpy((char*)tex.c_str(),Data+ini,len);
  }
  return(tex);
}

//==============================================================================
/// Returns next value in the current line or next line if in_line is false.
//==============================================================================
std::string JReadDatafile::ReadNextValue(bool in_line){
  if(in_line && ReadLine.empty())Run_ExceptioonFile(fun::PrintStr("Value %d does not exist in line %d.",ReadLinValue+1,ReadLin+1),File);
  while(ReadLine.empty() || ReadLine[0]=='#')SetReadLine(ReadLin+1);
  //printf("==>> ReadLine:[%s] Sep:[%u]\n",ReadLine.c_str(),Sep[0]);
  ReadValue=fun::StrSplit(Sep,ReadLine); 
  ReadLinValue++;
  return(ReadValue);
}

//==============================================================================
/// Returns next bool in the current line or next line if in_line is false.
//==============================================================================
bool JReadDatafile::ReadNextBool(bool in_line){
  bool ret;
  const string value=fun::StrLower(ReadNextValue(in_line));
  if(value=="1" || value=="true")ret=true;
  else if(value=="0" || value=="false")ret=false;
  else if(fun::StrIsRealNumber(value))ret=(atof(value.c_str())!=0);
  else ret=true;
  return(ret);
}

//==============================================================================
/// Returns next double in the current line or next line if in_line is false.
//==============================================================================
double JReadDatafile::ReadNextDouble(bool in_line){
  string value=ReadNextValue(in_line);
  return(atof(ReadValue.c_str()));
}

//==============================================================================
/// Returns next tdouble3 in the current line or next line if in_line is false.
//==============================================================================
tdouble3 JReadDatafile::ReadNextDouble3(bool in_line){
  tdouble3 v;
  v.x=ReadNextDouble(in_line);
  v.y=ReadNextDouble(true);
  v.z=ReadNextDouble(true);
  return(v);
}

//==============================================================================
/// Returns next int in the current line or next line if in_line is false.
//==============================================================================
int JReadDatafile::ReadNextInt(bool in_line){
  string value=ReadNextValue(in_line);
  return(atoi(ReadValue.c_str()));
}

//==============================================================================
/// Returns next tint3 in the current line or next line if in_line is false.
//==============================================================================
tint3 JReadDatafile::ReadNextInt3(bool in_line){
  tint3 v;
  v.x=ReadNextInt(in_line);
  v.y=ReadNextInt(true);
  v.z=ReadNextInt(true);
  return(v);
}

//==============================================================================
/// Returns next tuint3 in the current line or next line if in_line is false.
//==============================================================================
tuint3 JReadDatafile::ReadNextUnsigned3(bool in_line){
  tuint3 v;
  v.x=ReadNextUnsigned(in_line);
  v.y=ReadNextUnsigned(true);
  v.z=ReadNextUnsigned(true);
  return(v);
}

//==============================================================================
/// Returns line (and position inside the line) after firstline where key appears.
//==============================================================================
tint2 JReadDatafile::Find(std::string key,int firstline)const{
  tint2 v=TInt2(-1);
  for(int c=firstline;c<LineCount && v.x<0;c++){
    const string line=GetLine(c);
    const int pos=int(line.find(key));
    if(pos>=0)v=TInt2(c,pos);
  }
  return(v);
}

//==============================================================================
/// Returns string value after key string.
//==============================================================================
std::string JReadDatafile::FindValueStr(std::string key,bool optional,std::string valdef)const{
  string ret=valdef;
  tint2 res=Find(key);
  if(res.x<0 && !optional)Run_ExceptioonFile(fun::PrintStr("The KEY \'%s\' was not found.",key.c_str()),File);
  if(res.x>=0){
    //printf("--->Line[%s] res.y:%d \n",GetLine(res.x).c_str(),res.y);
    string value=GetLine(res.x).substr(res.y+key.size());
    //printf("--->1[%s]\n",value.c_str());
    if(value.size()>1 && value[0]=='\"' && value[value.size()-1]=='\"')value=value.substr(1,value.size()-2);
    //printf("--->2[%s]\n",value.c_str());
    ret=value;
  }
  //printf("----->RetStr[%s]\n",ret.c_str());
  return(ret);
}

//==============================================================================
/// Returns double value after key string.
//==============================================================================
double JReadDatafile::FindValueDbl(std::string key,bool optional,double valdef)const{
  double ret=valdef;
  tint2 res=Find(key);
  if(res.x<0 && !optional)Run_ExceptioonFile(fun::PrintStr("The KEY \'%s\' was not found.",key.c_str()),File);
  if(res.x>=0){
    string value=GetLine(res.x).substr(res.y+key.size());
    ret=atof(value.c_str());
  }
  //printf("----->RetDouble[%f]\n",ret);
  return(ret);
}



