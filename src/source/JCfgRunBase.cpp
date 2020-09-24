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

/// \file JCfgRunBase.cpp \brief Implements the class \ref JCfgRunBase.

#include "JCfgRunBase.h"
#include "JAppInfo.h"
#include "JDsphConfig.h"

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCfgRunBase::JCfgRunBase(bool noparms):NoParms(noparms){
  ClassName="JCfgRunBase";
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JCfgRunBase::Reset(){
  PrintInfo=false; ParmDef=0;
  CreateDirs=true;
  CsvSepComa=false;
}

//==============================================================================
/// Load configuration from DsphConfig.xml.
//==============================================================================
void JCfgRunBase::LoadDsphConfig(std::string path){
  JDsphConfig dsphconfig;
  dsphconfig.Init(path);
  if(!dsphconfig.GetFileCfg().empty())printf("LoadDsphConfig> %s\n",fun::GetPathLevels(dsphconfig.GetFileCfg(),3).c_str());
  if(dsphconfig.GetCreateDirs()!=-1)CreateDirs=(dsphconfig.GetCreateDirs()==1);
  if(dsphconfig.GetCsvSeparator()!=-1)CsvSepComa=(dsphconfig.GetCsvSeparator()==1);
}

//==============================================================================
/// Loads execution parameters from the command line.
//==============================================================================
void JCfgRunBase::LoadArgv(int argc,char** argv){
  Reset();
  //-Loads configuration from DsphConfig.xml.
  LoadDsphConfig(AppInfo.GetProgramPath());
  //-Loads execution parameters.
  const int MAXOPTS=100;
  string *optlis=new string[MAXOPTS];
  int optn=0;
  for(int c=0;c<argc-1;c++){
    string tex=fun::StrTrim(argv[c+1]);
    int pos=int(tex.find(" "));
    if(pos>0){
      while(pos>0){
        bool divide=((tex[0]=='-' || tex[0]=='#') || (pos+2<tex.size() && ((tex[pos+1]=='-' && tex[pos+2]!=' ') || tex[pos+1]=='#')));
        //printf("  tex[%s]  pos:%d  divide=%d\n",tex.c_str(),pos,(divide? 1: 0));
        if(divide){
          if(optn>=MAXOPTS)Run_Exceptioon("Has exceeded the maximum configuration options.");
          optlis[optn]=tex.substr(0,pos); optn++;
          tex=fun::StrTrim(tex.substr(pos+1)); //-StrTrim() removes spaces between options.
          pos=int(tex.find(" "));
        }
        else pos=int(tex.find(" ",pos+1));
      }
    }
    if(optn>=MAXOPTS)Run_Exceptioon("Has exceeded the maximum configuration options.");
    if(!tex.empty()){//-Ignores empty parameters.
      optlis[optn]=tex; optn++;
    }
  }
  //for(int c=0;c<optn;c++)printf("[%d]=[%s]\n",c,optlis[c].c_str());
  if(optn || NoParms)LoadOpts(optlis,optn,0,"");
  delete[] optlis;
  if(!optn && !NoParms)PrintInfo=true;
  if(!PrintInfo){ //-Configuracion por defecto
    //VisuConfig(); 
  }
  else VisuInfo();
}

//==============================================================================
/// Loads execution parameters from a text file.
//==============================================================================
void JCfgRunBase::LoadFile(std::string fname,int lv){
  //printf("\nFile:[%s] lv:%d\n",fname.c_str(),lv);
  const int MAXOPTS=50;
  int optn=0;
  std::string *optlis=new std::string[MAXOPTS];
  std::ifstream pf;
  pf.open(fname.c_str());
  if(pf){
    while(!pf.eof()&&optn<MAXOPTS){
      std::string tex;  pf >> tex;
      if(tex!=""){
        if(optn<MAXOPTS)optlis[optn]=tex;
        optn++;
      }
    } 
    if(!pf.eof() && pf.fail())Run_ExceptioonFile("Error reading data from the file.",fname);
    pf.close();
  }
  else Run_ExceptioonFile("The file can not be opened.",fname);
  if(optn>=MAXOPTS)Run_ExceptioonFile(fun::PrintStr("File with too many lines (Maximum=%d)",MAXOPTS),fname);
  if(optn>0)LoadOpts(optlis,optn,lv,fname);
  delete[] optlis;
}

//==============================================================================
/// Generates error of unknown parameter.
//==============================================================================
void JCfgRunBase::ErrorParm(const std::string &opt,int optc,int lv,const std::string &file)const{
  std::string tx=fun::PrintStr("Parameter \"%s\" unrecognised or invalid. ",opt.c_str());
  tx=tx+fun::PrintStr("(Level cfg:%d, Parameter:%d)",lv,optc);
  Run_ExceptioonFile(tx,file);
}


//==============================================================================
/// Loads nv values float using command options. Returns number of loaded 
/// values.
//==============================================================================
unsigned JCfgRunBase::LoadFloats(std::string txopt,float def,unsigned nv
  ,std::vector<float> &vv)
{
  //printf("txopt=[%s]\n",txopt.c_str());
  vv.clear();
  for(unsigned c=0;c<nv;c++)vv.push_back(def);
  unsigned num=0;
  std::string aux=txopt;
  for(;!aux.empty() && num<nv;num++){
    std::string txv=fun::StrSplit(":",aux);
    if(!txv.empty())vv[num]=float(atof(txv.c_str()));
  }
  return(num);
}

//==============================================================================
/// Loads nv values double using command options. Returns number of loaded 
/// values.
//==============================================================================
unsigned JCfgRunBase::LoadDoubles(std::string txopt,double def,unsigned nv
  ,std::vector<double> &vv)
{
  //printf("txopt=[%s]\n",txopt.c_str());
  vv.clear();
  for(unsigned c=0;c<nv;c++)vv.push_back(def);
  unsigned num=0;
  std::string aux=txopt;
  for(;!aux.empty() && num<nv;num++){
    std::string txv=fun::StrSplit(":",aux);
    if(!txv.empty())vv[num]=atof(txv.c_str());
  }
  return(num);
}

//==============================================================================
/// Loads nv values tfloat3 using command options. Returns number of loaded 
/// values.
//==============================================================================
unsigned JCfgRunBase::LoadFloat3(std::string txopt,float def,unsigned nv,tfloat3 *v){
  //printf("txopt=[%s]\n",txopt.c_str());
  unsigned num=0;
  float *values=(float*)v;
  for(unsigned c=0;c<nv*3;c++)values[c]=def;
  string ttx=txopt;
  for(unsigned tc=0;ttx!="" && tc<nv*3;tc++){
    int tpos=int(ttx.find(":"));
    string ttxopt=(tpos>0? ttx.substr(0,tpos): ttx);
    string ttxopt2;
    if(tpos>0)ttxopt2=ttx.substr(tpos+1);
    values[tc]=float(atof(ttxopt.c_str()));
    ttx=ttxopt2;
    num++;
  }
  return(num);
}

//==============================================================================
/// Loads nv values tdouble3 using command options. Returns number of loaded 
/// values.
//==============================================================================
unsigned JCfgRunBase::LoadDouble3(std::string txopt,double def,unsigned nv,tdouble3 *v){
  //printf("txopt=[%s]\n",txopt.c_str());
  unsigned num=0;
  double *values=(double*)v;
  for(unsigned c=0;c<nv*3;c++)values[c]=def;
  string ttx=txopt;
  for(unsigned tc=0;ttx!="" && tc<nv*3;tc++){
    int tpos=int(ttx.find(":"));
    string ttxopt=(tpos>0? ttx.substr(0,tpos): ttx);
    string ttxopt2;
    if(tpos>0)ttxopt2=ttx.substr(tpos+1);
    values[tc]=atof(ttxopt.c_str());
    ttx=ttxopt2;
    num++;
  }
  return(num);
}

//==============================================================================
/// Load 1 value tfloat3 using command options.
//==============================================================================
void JCfgRunBase::LoadFloat3(std::string txopt,float def,tfloat3 &v1){
  //printf("txopt=[%s]\n",txopt.c_str());
  tdouble3 vd;
  LoadDouble3(txopt,def,vd);
  v1=ToTFloat3(vd);
}

//==============================================================================
/// Load 2 values tfloat3 using command options.
//==============================================================================
void JCfgRunBase::LoadFloat6(std::string txopt,float def,tfloat3 &v1,tfloat3 &v2){
  tdouble3 vd1,vd2;
  LoadDouble6(txopt,def,vd1,vd2);
  v1=ToTFloat3(vd1);
  v2=ToTFloat3(vd2);
}

//==============================================================================
/// Load 1 value tdouble3 using command options.
//==============================================================================
void JCfgRunBase::LoadDouble3(std::string txopt,double def,tdouble3 &v1){
  //printf("txopt=[%s]\n",txopt.c_str());
  double values[3]={def,def,def};
  string ttx=txopt;
  for(int tc=0;ttx!="" && tc<3;tc++){
    int tpos=int(ttx.find(":"));
    string ttxopt=(tpos>0? ttx.substr(0,tpos): ttx);
    string ttxopt2;
    if(tpos>0)ttxopt2=ttx.substr(tpos+1);
    values[tc]=atof(ttxopt.c_str());
    ttx=ttxopt2;
  } 
  v1=TDouble3(values[0],values[1],values[2]);
}

//==============================================================================
/// Load 2 values tdouble3 using command options.
//==============================================================================
void JCfgRunBase::LoadDouble6(std::string txopt,double def,tdouble3 &v1,tdouble3 &v2){
  //printf("txopt=[%s]\n",txopt.c_str());
  double values[6]={def,def,def,def,def,def};
  string ttx=txopt;
  for(int tc=0;ttx!="" && tc<6;tc++){
    int tpos=int(ttx.find(":"));
    string ttxopt=(tpos>0? ttx.substr(0,tpos): ttx);
    string ttxopt2;
    if(tpos>0)ttxopt2=ttx.substr(tpos+1);
    values[tc]=atof(ttxopt.c_str());
    ttx=ttxopt2;
  } 
  v1=TDouble3(values[0],values[1],values[2]);
  v2=TDouble3(values[3],values[4],values[5]);
}

//==============================================================================
/// Splits options in txoptfull, txopt, txopt2, txopt3 and txopt4.
//==============================================================================
void JCfgRunBase::SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull
  ,std::string &txopt1,std::string &txopt2,std::string &txopt3,std::string &txopt4)const
{
  txword=txoptfull=txopt1=txopt2=txopt3=txopt4="";
  string tx=opt.substr(1);
  int pos=int(tx.find("#"));
  if(pos>0)tx=tx.substr(0,pos);
  pos=int(tx.find(":"));
  txword=fun::StrUpper(pos>0? tx.substr(0,pos): tx);
  if(pos>=0)txopt1=tx.substr(pos+1);
  txoptfull=txopt1;
  tx=txopt1;
  pos=int(tx.find(":"));
  txopt1=(pos>=0? tx.substr(0,pos): tx);
  if(pos>=0)txopt2=tx.substr(pos+1);
  tx=txopt2;
  pos=int(tx.find(":"));
  txopt2=(pos>=0? tx.substr(0,pos): tx);
  if(pos>=0)txopt3=tx.substr(pos+1);
  tx=txopt3;
  pos=int(tx.find(":"));
  txopt3=(pos>=0? tx.substr(0,pos): tx);
  if(pos>=0)txopt4=tx.substr(pos+1);
}


