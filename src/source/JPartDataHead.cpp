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

/// \file JPartDataHead.cpp \brief Implements the class \ref JPartDataHead.

#include "JPartDataHead.h"
//#include "JSpaceParts.h"
#include "JBinaryData.h"
#include "Functions.h"
#include <algorithm>
#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JPartDataHead
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartDataHead::JPartDataHead(){
  ClassName="JPartDataHead";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartDataHead::~JPartDataHead(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartDataHead::Reset(){
  DirData="";
  //-General variables.
  ConfigBasic("","","",TDouble3(0),TDouble3(0),false,0);
  ConfigCtes(0,0,0,0,0,0,0,TFloat3(0));
  ConfigSimPeri(false,false,false,false,false,false,TDouble3(0),TDouble3(0),TDouble3(0));
  ConfigSplitting(false);
  //-Mk blocks data.
  MkList.clear();
  UptadeMkNumbers();
}

//==============================================================================
/// Configuracion de variables basicas.
/// Configuration of basic variables.
//==============================================================================
void JPartDataHead::ConfigBasic(std::string runcode,std::string appname
  ,std::string casename,tdouble3 caseposmin,tdouble3 caseposmax
  ,bool data2d,double data2dposy)
{
  AppName=appname;
  Date=fun::GetDateTime();
  RunCode=runcode;
  CaseName=casename;
  Data2d=data2d;
  Data2dPosY=data2dposy;
  CasePosMin=caseposmin;
  CasePosMax=caseposmax;
}

//==============================================================================
/// Configuracion de informacion de cada bloque de particulas.
/// Setting information of each particle block.
//==============================================================================
void JPartDataHead::ConfigParticles(TpParticles type,unsigned mk,unsigned mktype,unsigned begin,unsigned count){
  const char met[]="ConfigParticles";
  if(MkListSize){
    const JPartDataHeadMkBlock mbk0=MkList[MkListSize-1];
    if(mbk0.Type>type)RunException(met,"Type of new particles block is invalid.");
    if(mbk0.Begin+mbk0.Count!=begin)RunException(met,"Begin of new particles block is invalid.");
    if(GetMkBlockByMk(mk)!=UINT_MAX)RunException(met,"Mk is already defined.");
    if(IsBound(type) && GetMkBlockByMkBound(mktype)!=UINT_MAX)RunException(met,"MkBound is already defined.");
    if(IsFluid(type) && GetMkBlockByMkFluid(mktype)!=UINT_MAX)RunException(met,"MkFluid is already defined.");
  }
  else if(begin)RunException(met,"Begin of first particles block is invalid.");
  if(!count)RunException(met,"Count of new particles block is invalid.");
  MkList.push_back(JPartDataHeadMkBlock(type,mk,mktype,begin,count));
  UptadeMkNumbers();
}

//==============================================================================
/// Actuliza numeros de bloques y particulas.
/// Updates numbers of Mk blocks and numbers of particles.
//==============================================================================
void JPartDataHead::UptadeMkNumbers(){
  MkListSize=unsigned(MkList.size());
  MkListFixed=MkListMoving=MkListFloat=MkListBound=MkListFluid=0;
  CaseNp=CaseNfixed=CaseNmoving=CaseNfloat=CaseNfluid=0;
  for(unsigned c=0;c<MkListSize;c++){
    const JPartDataHeadMkBlock& mbk=Mkblock(c);
    switch(mbk.Type){
      case TpPartFixed:     MkListFixed++;  CaseNfixed+=mbk.Count;   break;
      case TpPartMoving:    MkListMoving++; CaseNmoving+=mbk.Count;  break;
      case TpPartFloating:  MkListFloat++;  CaseNfloat+=mbk.Count;   break;
      case TpPartFluid:     MkListFluid++;  CaseNfluid+=mbk.Count;   break;
      default: RunException("UptadeMkNumbers","Type is unknown.");
    }
  }
  MkListBound=MkListFixed+MkListMoving+MkListFloat;
  CaseNp=CaseNfixed+CaseNmoving+CaseNfloat+CaseNfluid;
}

//==============================================================================
/// Configuracion de constantes.
/// Configuration of constants.
//==============================================================================
void JPartDataHead::ConfigCtes(double dp,double h,double b,double rhop0,double gamma
  ,double massbound,double massfluid,tfloat3 gravity)
{
  Dp=dp; H=h; B=b;
  RhopZero=rhop0;
  Gamma=gamma;
  MassBound=massbound;
  MassFluid=massfluid;
  Gravity=gravity;
}

//==============================================================================
/// Configuracion relativa al Np.
/// Sconfiguration related to the Np.
//==============================================================================
void JPartDataHead::ConfigSimNp(bool npdynamic,bool reuseids){
  NpDynamic=npdynamic;
  ReuseIds=reuseids;
}

//==============================================================================
/// Configuracion de variables de simulacion: map limits.
/// Configuration of variables of simulation: map limits.
//==============================================================================
void JPartDataHead::ConfigSimMap(tdouble3 mapposmin,tdouble3 mapposmax){
  MapPosMin=mapposmin;
  MapPosMax=mapposmax;
}

//==============================================================================
/// Configuracion de variables de condiciones periodicas.
/// Configuration of variables of periodic conditions.
//==============================================================================
void JPartDataHead::ConfigSimPeri(bool peri_x,bool peri_y,bool peri_z
  ,bool peri_xy,bool peri_xz,bool peri_yz
  ,tdouble3 perixinc,tdouble3 periyinc,tdouble3 perizinc)
{
  const char met[]="ConfigSimPeri";
  PeriActive=PERI_None;
  if(peri_xy)PeriActive=PERI_XY;
  else if(peri_xz)PeriActive=PERI_XZ;
  else if(peri_yz)PeriActive=PERI_YZ;
  else if(peri_x)PeriActive=PERI_X;
  else if(peri_y)PeriActive=PERI_Y;
  else if(peri_z)PeriActive=PERI_Z;
  PeriXinc=perixinc;
  PeriYinc=periyinc;
  PeriZinc=perizinc;
}

//==============================================================================
/// Configuracion uso de Splitting.
/// Configuration used for Splitting.
//==============================================================================
void JPartDataHead::ConfigSplitting(bool splitting){
  Splitting=splitting;
}

//==============================================================================
/// Saves binary file Part_Head.ibi4.
//==============================================================================
void JPartDataHead::SaveFile(std::string dir){
  DirData=fun::GetDirWithSlash(dir);
  JBinaryData bdat(ClassName);
  bdat.SetvUint("FmtVersion",FmtVersion);
  //-Saves general variables.
  bdat.SetvText("AppName",AppName);
  bdat.SetvText("Date",Date);
  bdat.SetvText("RunCode",RunCode);
  bdat.SetvText("CaseName",CaseName);
  bdat.SetvBool("Data2d",Data2d);
  bdat.SetvDouble("Data2dPosY",Data2dPosY);

  bdat.SetvDouble3("CasePosMin",CasePosMin);
  bdat.SetvDouble3("CasePosMax",CasePosMax);

  bdat.SetvBool("NpDynamic",NpDynamic);
  bdat.SetvBool("ReuseIds",ReuseIds);

  bdat.SetvDouble3("MapPosMin",MapPosMin);
  bdat.SetvDouble3("MapPosMax",MapPosMax);

  bdat.SetvInt("PeriActive",int(PeriActive));
  bdat.SetvDouble3("PeriXinc",PeriXinc);
  bdat.SetvDouble3("PeriYinc",PeriYinc);
  bdat.SetvDouble3("PeriZinc",PeriZinc);

  bdat.SetvBool("Splitting",Splitting);

  bdat.SetvDouble("Dp",Dp);
  bdat.SetvDouble("H",H);
  bdat.SetvDouble("B",B);
  bdat.SetvDouble("Gamma",Gamma);
  bdat.SetvDouble("RhopZero",RhopZero);
  bdat.SetvDouble("MassBound",MassBound);
  bdat.SetvDouble("MassFluid",MassFluid);
  bdat.SetvFloat3("Gravity",Gravity);

  //-Number of particles (is optional).
  bdat.SetvUllong("CaseNp",CaseNp);
  bdat.SetvUllong("CaseNfixed",CaseNfixed);
  bdat.SetvUllong("CaseNmoving",CaseNmoving);
  bdat.SetvUllong("CaseNfloat",CaseNfloat);
  bdat.SetvUllong("CaseNfluid",CaseNfluid);

  //-Saves information of particle blocks.
  JBinaryData* bdatamk=bdat.CreateItem("MkBlocks");
  bdatamk->SetvUint("Count",MkListSize);
  for(unsigned c=0;c<MkListSize;c++){
    const JPartDataHeadMkBlock& mbk=MkList[c];
    JBinaryData* item=bdatamk->CreateItem(fun::PrintStr("MkBlock_%03u",c));
    item->SetvText("Type"  ,TpPartGetStrCode(mbk.Type));
    item->SetvUint("Mk"    ,mbk.Mk);
    item->SetvUint("MkType",mbk.MkType);
    item->SetvUint("Count" ,mbk.Count);
  }
  bdat.SaveFile(DirData+GetFileName(),false,true);
}

//==============================================================================
/// Loads binary file Part_Head.ibi4.
//==============================================================================
void JPartDataHead::LoadFile(std::string dir){
  const char met[]="LoadFile";
  Reset();
  DirData=fun::GetDirWithSlash(dir);
  const string file=DirData+GetFileName();
  JBinaryData bdat;
  bdat.LoadFile(file,ClassName);
  FmtVersion=bdat.GetvUint("FmtVersion");
  //-Loads general variables.
  AppName   =bdat.GetvText("AppName");
  Date      =bdat.GetvText("Date");
  RunCode   =bdat.GetvText("RunCode");
  CaseName  =bdat.GetvText("CaseName");
  Data2d    =bdat.GetvBool("Data2d");
  Data2dPosY=bdat.GetvDouble("Data2dPosY");

  CasePosMin=bdat.GetvDouble3("CasePosMin");
  CasePosMax=bdat.GetvDouble3("CasePosMax");

  NpDynamic=bdat.GetvBool("NpDynamic");
  ReuseIds=bdat.GetvBool("ReuseIds");

  MapPosMin=bdat.GetvDouble3("MapPosMin");
  MapPosMax=bdat.GetvDouble3("MapPosMax");

  PeriActive=(TpPeri)bdat.GetvInt("PeriActive");
  PeriXinc=bdat.GetvDouble3("PeriXinc");
  PeriYinc=bdat.GetvDouble3("PeriYinc");
  PeriZinc=bdat.GetvDouble3("PeriZinc");

  Splitting=bdat.GetvBool("Splitting");

  Dp=bdat.GetvDouble("Dp");
  H=bdat.GetvDouble("H");
  B=bdat.GetvDouble("B");
  Gamma=bdat.GetvDouble("Gamma");
  RhopZero=bdat.GetvDouble("RhopZero");
  MassBound=bdat.GetvDouble("MassBound");
  MassFluid=bdat.GetvDouble("MassFluid");
  Gravity=bdat.GetvFloat3("Gravity");

  //-Loads information of particle blocks.
  JBinaryData* bdatamk=bdat.GetItem("MkBlocks");
  if(!bdatamk)RunException(met,"Item \'MkBlocks\' is missing.",file);
  const unsigned nmk=bdatamk->GetvUint("Count");
  unsigned begin=0;
  for(unsigned c=0;c<nmk;c++){
    const JBinaryData* item=bdatamk->GetItem(fun::PrintStr("MkBlock_%03u",c));
    if(!item)RunException(met,fun::PrintStr("Item \'MkBlock_%03u\' is missing.",c),file);
    const TpParticles type=TpPartGetType(item->GetvText("Type"));
    const unsigned mk=item->GetvUint("Mk");
    const unsigned mktype=item->GetvUint("MkType");
    const unsigned count=item->GetvUint("Count");
    ConfigParticles(type,mk,mktype,begin,count);
    begin+=count;
  }

}

//==============================================================================
/// Returns the block in MkList according to a given MK.
//==============================================================================
unsigned JPartDataHead::GetMkBlockByMk(unsigned mk)const{
  unsigned c=0;
  for(;c<MkListSize && mk!=MkList[c].Mk;c++);
  return(c<MkListSize? c: UINT_MAX);
}

//==============================================================================
/// Returns the block in MkList according to a given bound MK.
//==============================================================================
unsigned JPartDataHead::GetMkBlockByMkBound(unsigned mkbound)const{
  unsigned c=0;
  for(;c<MkListBound && mkbound!=MkList[c].MkType;c++);
  return(c<MkListBound? c: UINT_MAX);
}

//==============================================================================
/// Returns the block in MkList according to a given fluid MK.
//==============================================================================
unsigned JPartDataHead::GetMkBlockByMkFluid(unsigned mkfluid)const{
  unsigned c=MkListBound;
  for(;c<MkListSize && mkfluid!=MkList[c].MkType;c++);
  return(c<MkListSize? c: UINT_MAX);
}

//==============================================================================
/// Returns the block in MkList according to a given Id.
//==============================================================================
unsigned JPartDataHead::GetMkBlockById(unsigned id)const{
  unsigned c=0;
  for(;c<MkListSize && id>=(MkList[c].Begin+MkList[c].Count);c++);
  return(c);
}

