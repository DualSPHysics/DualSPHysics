//HEAD_DSCODES
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

/// \file JPartFloatInfoBi4.cpp \brief Implements the classes JPartFloatInfoBi4Save and class JPartFloatInfoBi4Load.

#include "JPartFloatInfoBi4.h"
#include "Functions.h"
#include <cstring>

using namespace std;

//##############################################################################
//# JPartFloatInfoBi4Data
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartFloatInfoBi4Data::JPartFloatInfoBi4Data(unsigned ftcount,unsigned fptcount
  ,word mkboundfirst):FtCount(ftcount),FptCount(fptcount),MkBoundFirst(mkboundfirst)
{
  ClassName="JPartFloatInfoBi4Data";
  //-Basic and constant floating data.
  CteData=new StFtDataCte[FtCount];
  //-Data variables of PARTs.
  PartNum=NULL;
  PartTime=NULL;
  PartStep=NULL;
  //-Data variables of floating bodies.
  PartCenter=NULL;
  PartVelLin=NULL;      PartVelAng=NULL;
  PartAceLin=NULL;      PartAceAng=NULL;
  PartExtForceLin=NULL; PartExtForceAng=NULL;
  PartFluForceLin=NULL; PartFluForceAng=NULL;
  PartPreAceLin=NULL;   PartPreAceAng=NULL;
  //-Data of force points.
  FptMkbound=NULL;
  FptPos=NULL;
  FptForce=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartFloatInfoBi4Data::~JPartFloatInfoBi4Data(){
  DestructorActive=true;
  Reset();
  delete[] CteData; CteData=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartFloatInfoBi4Data::Reset(){
  ResizeData(0);
}

//==============================================================================
/// Resizes memory space for data items.
//==============================================================================
void JPartFloatInfoBi4Data::ResizeData(unsigned size){
  //-Free memory.
  if(!size){
    PartSize=PartCount=0;
    //-Data variables of PARTs.
    delete[] PartNum;     PartNum=NULL;
    delete[] PartTime;    PartTime=NULL;
    delete[] PartStep;    PartStep=NULL;
    //-Data variables of floating bodies.
    delete[] PartCenter;  PartCenter=NULL;
    delete[] PartVelLin;  PartVelLin=NULL;
    delete[] PartVelAng;  PartVelAng=NULL;
    delete[] PartAceLin;  PartAceLin=NULL;
    delete[] PartAceAng;  PartAceAng=NULL;
    //-Extra data.
    delete[] PartExtForceLin; PartExtForceLin=NULL;
    delete[] PartExtForceAng; PartExtForceAng=NULL;
    delete[] PartFluForceLin; PartFluForceLin=NULL;
    delete[] PartFluForceAng; PartFluForceAng=NULL;
    delete[] PartPreAceLin;   PartPreAceLin=NULL;
    delete[] PartPreAceAng;   PartPreAceAng=NULL;
    //-Data of force points.
    delete[] FptMkbound;  FptMkbound=NULL;
    delete[] FptPos;      FptPos=NULL;
    delete[] FptForce;    FptForce=NULL;
  }
  else{
    if(PartCount>size)Run_Exceptioon("No memory available for stored data.");
    PartSize=size;
    try{
      //-Data variables of PARTs.
      PartNum =fun::ResizeAlloc(PartNum ,PartCount,PartSize);
      PartTime=fun::ResizeAlloc(PartTime,PartCount,PartSize);
      PartStep=fun::ResizeAlloc(PartStep,PartCount,PartSize);
      //-Data variables of floating bodies.
      const unsigned scount=FtCount*PartCount;
      const unsigned ssize =FtCount*PartSize;
      PartCenter=fun::ResizeAlloc(PartCenter,scount,ssize);
      PartVelLin=fun::ResizeAlloc(PartVelLin,scount,ssize);
      PartVelAng=fun::ResizeAlloc(PartVelAng,scount,ssize);
      PartAceLin=fun::ResizeAlloc(PartAceLin,scount,ssize);
      PartAceAng=fun::ResizeAlloc(PartAceAng,scount,ssize);
      //-Extra data.
      PartExtForceLin=fun::ResizeAlloc(PartExtForceLin,scount,ssize);
      PartExtForceAng=fun::ResizeAlloc(PartExtForceAng,scount,ssize);
      PartFluForceLin=fun::ResizeAlloc(PartFluForceLin,scount,ssize);
      PartFluForceAng=fun::ResizeAlloc(PartFluForceAng,scount,ssize);
      PartPreAceLin  =fun::ResizeAlloc(PartPreAceLin  ,scount,ssize);
      PartPreAceAng  =fun::ResizeAlloc(PartPreAceAng  ,scount,ssize);
      //-Data of force points.
      if(FptCount){
        const unsigned spcount=FptCount*PartCount;
        const unsigned spsize =FptCount*PartSize;
        FptMkbound=fun::ResizeAlloc(FptMkbound,spcount,spsize);
        FptPos    =fun::ResizeAlloc(FptPos    ,spcount,spsize);
        FptForce  =fun::ResizeAlloc(FptForce  ,spcount,spsize);
      }
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Clears floating data.
//==============================================================================
void JPartFloatInfoBi4Data::ClearFtData(bool clearctedata){
  if(clearctedata && CteData)memset(CteData,0,sizeof(StFtDataCte)*FtCount);
  if(PartCenter)memset(PartCenter,0,sizeof(tdouble3)*FtCount);
  if(PartVelLin)memset(PartVelLin,0,sizeof(tfloat3) *FtCount);
  if(PartVelAng)memset(PartVelAng,0,sizeof(tfloat3) *FtCount);
  if(PartAceLin)memset(PartAceLin,0,sizeof(tfloat3) *FtCount);
  if(PartAceAng)memset(PartAceAng,0,sizeof(tfloat3) *FtCount);
  if(PartExtForceLin)memset(PartExtForceLin,0,sizeof(tfloat3)*FtCount);
  if(PartExtForceAng)memset(PartExtForceAng,0,sizeof(tfloat3)*FtCount);
  if(PartFluForceLin)memset(PartFluForceLin,0,sizeof(tfloat3)*FtCount);
  if(PartFluForceAng)memset(PartFluForceAng,0,sizeof(tfloat3)*FtCount);
  if(PartPreAceLin  )memset(PartPreAceLin  ,0,sizeof(tfloat3)*FtCount);
  if(PartPreAceAng  )memset(PartPreAceAng  ,0,sizeof(tfloat3)*FtCount);
}

//==============================================================================
/// Devuelve la memoria reservada.
/// Returns allocated memory
//==============================================================================
llong JPartFloatInfoBi4Data::GetAllocMemory()const{  
  llong s=0;
  //-Basic and constant floating data [FtCount].
  if(CteData)s=s+sizeof(StFtDataCte)*FtCount;
  //-Data variables of PARTs [PartSize].
  if(PartNum) s=s+sizeof(int)     *PartSize;
  if(PartTime)s=s+sizeof(double)  *PartSize;
  if(PartStep)s=s+sizeof(unsigned)*PartSize;
  //-Data variables of floating bodies [PartSize*FtCount].
  const unsigned ssize=FtCount*PartSize;
  if(PartCenter) s=s+sizeof(tdouble3)*ssize;
  if(PartVelLin) s=s+sizeof(tfloat3) *ssize;
  if(PartVelAng) s=s+sizeof(tfloat3) *ssize;
  if(PartAceLin) s=s+sizeof(tfloat3) *ssize;
  if(PartAceAng) s=s+sizeof(tfloat3) *ssize;
  //-Extra data [PartSize*FtCount].
  if(PartExtForceLin)s=s+sizeof(tfloat3)*ssize;
  if(PartExtForceAng)s=s+sizeof(tfloat3)*ssize;
  if(PartFluForceLin)s=s+sizeof(tfloat3)*ssize;
  if(PartFluForceAng)s=s+sizeof(tfloat3)*ssize;
  if(PartPreAceLin)  s=s+sizeof(tfloat3)*ssize;
  if(PartPreAceAng)  s=s+sizeof(tfloat3)*ssize;
  //-Data of force points [PartSize*FptSize].
  const unsigned spsize=FptCount*PartSize;
  if(FptMkbound) s=s+sizeof(word)    *spsize;
  if(FptPos)     s=s+sizeof(tdouble3)*spsize;
  if(FptForce)   s=s+sizeof(tfloat3) *spsize;
  return(s);
}

//==============================================================================
/// Set basic floating data.
//==============================================================================
void JPartFloatInfoBi4Data::SetCteData(unsigned cf,word mkbound,unsigned begin
  ,unsigned count,float mass,float massp,float radius)
{
  if(cf>=FtCount)Run_Exceptioon("Index of floating body is invalid.");
  StFtDataCte& v=CteData[cf];
  v.mkbound=mkbound;
  v.begin=begin;
  v.count=count;
  v.mass=mass;
  v.massp=massp;
  v.radius=radius;
}

//==============================================================================
/// Load basic floating data from ftdata object.
//==============================================================================
void JPartFloatInfoBi4Data::LoadCteData(const JPartFloatInfoBi4Data* ftdata){
  if(ftdata->FtCount!=FtCount)Run_Exceptioon("Number of floating bodies does not match.");
  //if(ftdata->FptCount!=FptCount)Run_Exceptioon("Number of force points does not match.");
  for(unsigned cf=0;cf<FtCount;cf++)CteData[cf]=ftdata->CteData[cf];
}

//==============================================================================
/// Set simulation data of floating bodies.
//==============================================================================
void JPartFloatInfoBi4Data::SetPartData0(const JPartFloatInfoBi4Data* ftdata)
{
  if(!PartSize)Run_Exceptioon("No allocated memory for data.");
  if(ftdata->GetFtCount()!=FtCount)Run_Exceptioon("Number of floating bodies does not match.");
  if(ftdata->GetPartCount()!=1)Run_Exceptioon("Source data is invalid.");
  for(unsigned cf=0;cf<FtCount;cf++){
    PartCenter[cf]=ftdata->PartCenter[cf];
    PartVelLin[cf]=ftdata->PartVelLin[cf];
    PartVelAng[cf]=ftdata->PartVelAng[cf];
    PartAceLin[cf]=ftdata->PartAceLin[cf];
    PartAceAng[cf]=ftdata->PartAceAng[cf];
    PartExtForceLin[cf]=ftdata->PartExtForceLin[cf];
    PartExtForceAng[cf]=ftdata->PartExtForceAng[cf];
    PartFluForceLin[cf]=ftdata->PartFluForceLin[cf];
    PartFluForceAng[cf]=ftdata->PartFluForceAng[cf];
    PartPreAceLin  [cf]=ftdata->PartPreAceLin  [cf];
    PartPreAceAng  [cf]=ftdata->PartPreAceAng  [cf];
  }
}

//==============================================================================
/// Set simulation data of floating bodies.
//==============================================================================
void JPartFloatInfoBi4Data::SetPartData0(unsigned cf,const tdouble3& center
  ,const tfloat3& fvellin,const tfloat3& fvelang
  ,const tfloat3& facelin,const tfloat3& faceang
  ,const tfloat3& extforcelin,const tfloat3& extforceang
  ,const tfloat3& fluforcelin,const tfloat3& fluforceang
  ,const tfloat3& preacelin,const tfloat3& preaceang)
{
  if(!PartSize)Run_Exceptioon("No allocated memory for data.");
  if(cf>=FtCount)Run_Exceptioon("Index of floating body is invalid.");
  PartCenter[cf]=center;
  PartVelLin[cf]=fvellin;
  PartVelAng[cf]=fvelang;
  PartAceLin[cf]=facelin;
  PartAceAng[cf]=faceang;
  PartExtForceLin[cf]=extforcelin;
  PartExtForceAng[cf]=extforceang;
  PartFluForceLin[cf]=fluforcelin;
  PartFluForceAng[cf]=fluforceang;
  PartPreAceLin  [cf]=preacelin;
  PartPreAceAng  [cf]=preaceang;
}

//==============================================================================
/// Adds data of force points to new part.
//==============================================================================
void JPartFloatInfoBi4Data::SetForcePoints0(unsigned npt
  ,const word* mkbound,const tdouble3* pos,const tfloat3* force)
{
  if(!PartSize)Run_Exceptioon("No allocated memory for data.");
  if(npt!=FptCount)Run_Exceptioon("Number of force points does not match.");
  memcpy(FptMkbound,mkbound,sizeof(word    )*npt);
  memcpy(FptPos    ,pos    ,sizeof(tdouble3)*npt);
  memcpy(FptForce  ,force  ,sizeof(tfloat3 )*npt);
}

//==============================================================================
/// Set time of data.
//==============================================================================
void JPartFloatInfoBi4Data::SetTimeData0(int cpart,double timestep,unsigned step){
  if(!PartSize)Run_Exceptioon("No allocated memory for data.");
  PartNum [0]=cpart;
  PartTime[0]=timestep;
  PartStep[0]=step;
  PartCount=1;
}

//==============================================================================
/// Load floating data from ftdata object.
//==============================================================================
void JPartFloatInfoBi4Data::AddData(const JPartFloatInfoBi4Data* ftdata){
  if(ftdata->FtCount!=FtCount)Run_Exceptioon("Number of floating bodies does not match.");
  if(ftdata->FptCount!=FptCount)Run_Exceptioon("Number of force points does not match.");
  const unsigned count2=ftdata->PartCount;
  if(count2){
    if(PartCount+count2>PartSize)ResizeData(max(PartCount+count2,max(PartSize*2,10u)));
    //-Data variables of PARTs.
    memcpy(PartNum +PartCount,ftdata->PartNum ,sizeof(int)     *count2);
    memcpy(PartTime+PartCount,ftdata->PartTime,sizeof(double)  *count2);
    memcpy(PartStep+PartCount,ftdata->PartStep,sizeof(unsigned)*count2);
    //-Data variables of floating bodies.
    const unsigned fcount1=PartCount*FtCount;
    const unsigned fcount2=ftdata->PartCount*FtCount;
    memcpy(PartCenter+fcount1,ftdata->PartCenter,sizeof(tdouble3)*fcount2);
    memcpy(PartVelLin+fcount1,ftdata->PartVelLin,sizeof(tfloat3)*fcount2);
    memcpy(PartVelAng+fcount1,ftdata->PartVelAng,sizeof(tfloat3)*fcount2);
    memcpy(PartAceLin+fcount1,ftdata->PartAceLin,sizeof(tfloat3)*fcount2);
    memcpy(PartAceAng+fcount1,ftdata->PartAceAng,sizeof(tfloat3)*fcount2);
    //-Extra data.
    memcpy(PartExtForceLin+fcount1,ftdata->PartExtForceLin,sizeof(tfloat3)*fcount2);
    memcpy(PartExtForceAng+fcount1,ftdata->PartExtForceAng,sizeof(tfloat3)*fcount2);
    memcpy(PartFluForceLin+fcount1,ftdata->PartFluForceLin,sizeof(tfloat3)*fcount2);
    memcpy(PartFluForceAng+fcount1,ftdata->PartFluForceAng,sizeof(tfloat3)*fcount2);
    memcpy(PartPreAceLin  +fcount1,ftdata->PartPreAceLin  ,sizeof(tfloat3)*fcount2);
    memcpy(PartPreAceAng  +fcount1,ftdata->PartPreAceAng  ,sizeof(tfloat3)*fcount2);
    //-Data of force points
    const unsigned pcount1=PartCount*FptCount;
    const unsigned pcount2=ftdata->PartCount*FptCount;
    memcpy(FptMkbound+pcount1,ftdata->FptMkbound,sizeof(word)    *pcount2);
    memcpy(FptPos    +pcount1,ftdata->FptPos    ,sizeof(tdouble3)*pcount2);
    memcpy(FptForce  +pcount1,ftdata->FptForce  ,sizeof(tfloat3) *pcount2);
    PartCount+=ftdata->PartCount;
  }
}

//==============================================================================
/// Devuelve el indice del numero de PART solicitado (comprueba si es unico).
/// Verifies that the number of floating is valid (checks if it is unique).
//==============================================================================
unsigned JPartFloatInfoBi4Data::GetIdxPart(unsigned cpart,bool unique)const{
  unsigned fidx=0;
  unsigned nidx=0;
  for(unsigned idx=0;idx<PartCount;idx++){
    if(PartNum[idx]==cpart){
      if(!nidx)fidx=idx;
      nidx++;
    }
  }
  if(!nidx || (nidx>1 && unique))fidx=UINT_MAX;
  return(fidx);
}

//==============================================================================
/// Comprueba que el numero de floating sea valido.
/// Verifies that the number of floating is valid.
//==============================================================================
unsigned JPartFloatInfoBi4Data::GetCf(unsigned cf)const{
  if(cf>=FtCount)Run_Exceptioon("Number of floating is invalid.");
  return(cf);
}

//==============================================================================
/// Comprueba que el numero de elemento sea valido.
/// Verifies that the number of item is valid.
//==============================================================================
unsigned JPartFloatInfoBi4Data::GetIdx(unsigned idx)const{
  if(idx>=PartCount)Run_Exceptioon("Index is invalid.");
  return(idx);
}

//==============================================================================
/// Comprueba que el numero de puntos de fuerza sea valido.
/// Verifies that the number of force points is valid.
//==============================================================================
unsigned JPartFloatInfoBi4Data::GetCpt(unsigned cpt)const{
  if(cpt>=FptCount)Run_Exceptioon("Number of force points is invalid.");
  return(cpt);
}


//##############################################################################
//# JPartFloatInfoBi4Save
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartFloatInfoBi4Save::JPartFloatInfoBi4Save(std::string appname,std::string dir)
  :AppName(appname),Dir(fun::GetDirWithSlash(dir))
{
  ClassName="JPartFloatInfoBi4Save";
  BdData=NULL;
  BdPart=NULL;
  FtData=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartFloatInfoBi4Save::~JPartFloatInfoBi4Save(){
  DestructorActive=true;
  Reset();
  delete BdData; BdData=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartFloatInfoBi4Save::Reset(){
  FormatVer=FormatVerDef;
  MainFile=false;
  TimeOut=-1;
  FileFull="";
  MkBoundFirst=0;
  FtCount=0;
  FptCount=0;
  InitialSaved=false;
  ResetBdData();
  delete FtData; FtData=NULL;
}

//==============================================================================
/// Elimina informacion de Data.
/// Deletes information from data. 
//==============================================================================
void JPartFloatInfoBi4Save::ResetBdData(){
  delete BdData; 
  BdData=new JBinaryData("JPartFloatInfoBi4");
  BdPart=BdData->CreateItem("Part");
  PartCount=0;
}

//==============================================================================
/// Elimina informacion de PARTs.
/// Deletes information from PARTs.
//==============================================================================
void JPartFloatInfoBi4Save::ResetBdPart(){
  BdPart->Clear();
}

//==============================================================================
/// Devuelve la memoria reservada.
/// Returns allocated memory
//==============================================================================
llong JPartFloatInfoBi4Save::GetAllocMemory()const{  
  llong s=0;
  if(BdData)s+=BdData->GetAllocMemory();
  if(FtData)s+=FtData->GetAllocMemory();
  return(s);
}

//==============================================================================
/// Devuelve nombre de part segun su numero.
/// Returns name of part according to their number.
//==============================================================================
std::string JPartFloatInfoBi4Save::GetNamePart(unsigned cpart){
  return(fun::PrintStr("PART_%04u",cpart));
}

//==============================================================================
/// Devuelve nombre de fichero bi4 por defecto.
/// Returns the filename bi4 by default.
//==============================================================================
std::string JPartFloatInfoBi4Save::GetFileNameDef(bool mainfile,std::string dir){
  if(!dir.empty())dir=fun::GetDirWithSlash(dir);
  return(dir+(mainfile? "PartFloatInfo.ibi4": "PartFloatInfo2.ibi4"));
}

//==============================================================================
/// Configuracion de datos de cabecera.
/// Configuration of header data.
//==============================================================================
void JPartFloatInfoBi4Save::Config(bool mainfile,double timeout
  ,const JPartFloatInfoBi4Data* ftdata)
{
  Reset();
  MainFile=mainfile;
  TimeOut=timeout;
  FileFull=Dir+GetFileNameDef(MainFile);
  MkBoundFirst=ftdata->MkBoundFirst;
  FtCount=ftdata->FtCount;
  FptCount=ftdata->FptCount;
  if(!FtCount)Run_Exceptioon("No floating bodies.");
  //-Save head data in Data.
  BdData->SetvText("AppName",AppName);
  BdData->SetvUint("FormatVer",FormatVer);
  BdData->SetvBool("MainFile",MainFile);
  BdData->SetvDouble("TimeOut",TimeOut);
  BdData->SetvUshort("MkBoundFirst",MkBoundFirst);
  BdData->SetvUint("FtCount",FtCount);
  BdData->SetvUint("FptCount",FptCount);
  //-Creates head arrays.
  word*     pmkb =BdData->CreateArrayUshort ("MkBound",FtCount,true);
  unsigned* pbeg =BdData->CreateArrayUint   ("Beginp" ,FtCount,true);
  unsigned* pnp  =BdData->CreateArrayUint   ("Countp" ,FtCount,true);
  float*    pmass=BdData->CreateArrayFloat  ("Mass"   ,FtCount,true);
  float*    massp=BdData->CreateArrayFloat  ("Massp"  ,FtCount,true);
  float*    pradi=BdData->CreateArrayFloat  ("Radius" ,FtCount,true);
  //-Loads head arrays.
  for(unsigned cf=0;cf<FtCount;cf++){
    const JPartFloatInfoBi4Data::StFtDataCte v=ftdata->GetHead(cf);
    pmkb [cf]=v.mkbound;
    pbeg [cf]=v.begin;
    pnp  [cf]=v.count;
    pmass[cf]=v.mass;
    massp[cf]=v.massp;
    pradi[cf]=v.radius;
  }
  BdPart->SetHide(true);
  //-Creates FtData object.
  FtData=new JPartFloatInfoBi4Data(FtCount,FptCount,MkBoundFirst);
  FtData->LoadCteData(ftdata);
}

//==============================================================================
/// Grabacion inicial de fichero con info de Data.
/// Initial file recording with info from Data.
//==============================================================================
void JPartFloatInfoBi4Save::SaveInitial(){
  if(!InitialSaved){
    BdData->SaveFile(FileFull,true,false);
    InitialSaved=true;
  }
}

//==============================================================================
/// Adds data to BdPart object.
//==============================================================================
JBinaryData* JPartFloatInfoBi4Save::MakeBdPartSingle(int cpart,double timestep
  ,unsigned step,const JPartFloatInfoBi4Data* ftdata)
{
  //-Check step of ftdata.
  //printf("*** MakeBdPartSingle> PartSize:%u  step:%u\n\n",ftdata->PartSize,ftdata->PartStep[0]);
  if(ftdata->PartCount!=1 || step!=ftdata->PartStep[0])
    Run_Exceptioon("Step number does not match the current data or data invalid.");
  //-Configures item Part.
  BdPart->Clear();
  BdPart->SetName(GetNamePart(cpart));
  BdPart->SetvUint("Cpart",cpart);
  BdPart->SetvDouble("TimeStep",timestep);
  BdPart->SetvUint("Step",step);
  //-Crea array con datos de floatings. Create array with floatings data.
  const bool external=true;
  BdPart->CreateArray("center" ,JBinaryDataDef::DatDouble3,FtCount,ftdata->PartCenter,external);
  BdPart->CreateArray("fvel"   ,JBinaryDataDef::DatFloat3 ,FtCount,ftdata->PartVelLin,external);
  BdPart->CreateArray("fomega" ,JBinaryDataDef::DatFloat3 ,FtCount,ftdata->PartVelAng,external);
  BdPart->CreateArray("facelin",JBinaryDataDef::DatFloat3 ,FtCount,ftdata->PartAceLin,external);
  BdPart->CreateArray("faceang",JBinaryDataDef::DatFloat3 ,FtCount,ftdata->PartAceAng,external);
  BdPart->CreateArray("extforcelin",JBinaryDataDef::DatFloat3,FtCount,ftdata->PartExtForceLin,external);
  BdPart->CreateArray("extforceang",JBinaryDataDef::DatFloat3,FtCount,ftdata->PartExtForceAng,external);
  BdPart->CreateArray("fluforcelin",JBinaryDataDef::DatFloat3,FtCount,ftdata->PartFluForceLin,external);
  BdPart->CreateArray("fluforceang",JBinaryDataDef::DatFloat3,FtCount,ftdata->PartFluForceAng,external);
  BdPart->CreateArray("preacelin"  ,JBinaryDataDef::DatFloat3,FtCount,ftdata->PartPreAceLin  ,external);
  BdPart->CreateArray("preaceang"  ,JBinaryDataDef::DatFloat3,FtCount,ftdata->PartPreAceAng  ,external);
  //-Crea arrays con datos de force points. Create arrays with force points data.
  if(ftdata->FptCount){
    BdPart->CreateArray("FptMkbound",JBinaryDataDef::DatUshort ,ftdata->FptCount,ftdata->FptMkbound,external);
    BdPart->CreateArray("FptPos"    ,JBinaryDataDef::DatDouble3,ftdata->FptCount,ftdata->FptPos    ,external);
    BdPart->CreateArray("FptForce"  ,JBinaryDataDef::DatFloat3 ,ftdata->FptCount,ftdata->FptForce  ,external);
  }
  return(BdPart);
}

//==============================================================================
/// Adds data to BdPart object.
//==============================================================================
JBinaryData* JPartFloatInfoBi4Save::MakeBdPartArray(const JPartFloatInfoBi4Data* ftdata)
{
  const unsigned count=ftdata->PartCount;
  const int cpart=ftdata->PartNum[0];
  //-Configures item Part.
  BdPart->Clear();
  BdPart->SetName(GetNamePart(cpart));
  BdPart->SetvUint("Cpart",cpart);
  BdPart->SetvUint("Count",count);
  const bool external=true;
  BdPart->CreateArray("TimeStep",JBinaryDataDef::DatDouble,count,ftdata->PartTime,external);
  BdPart->CreateArray("Step"    ,JBinaryDataDef::DatUint,count,ftdata->PartStep,external);
  //-Crea array con datos de floatings. Create array with floatings data.
  const unsigned scount=FtCount*count;
  BdPart->CreateArray("center" ,JBinaryDataDef::DatDouble3,scount,ftdata->PartCenter,external);
  BdPart->CreateArray("fvel"   ,JBinaryDataDef::DatFloat3 ,scount,ftdata->PartVelLin,external);
  BdPart->CreateArray("fomega" ,JBinaryDataDef::DatFloat3 ,scount,ftdata->PartVelAng,external);
  BdPart->CreateArray("facelin",JBinaryDataDef::DatFloat3 ,scount,ftdata->PartAceLin,external);
  BdPart->CreateArray("faceang",JBinaryDataDef::DatFloat3 ,scount,ftdata->PartAceAng,external);
  BdPart->CreateArray("extforcelin",JBinaryDataDef::DatFloat3,scount,ftdata->PartExtForceLin,external);
  BdPart->CreateArray("extforceang",JBinaryDataDef::DatFloat3,scount,ftdata->PartExtForceAng,external);
  BdPart->CreateArray("fluforcelin",JBinaryDataDef::DatFloat3,scount,ftdata->PartFluForceLin,external);
  BdPart->CreateArray("fluforceang",JBinaryDataDef::DatFloat3,scount,ftdata->PartFluForceAng,external);
  BdPart->CreateArray("preacelin"  ,JBinaryDataDef::DatFloat3,scount,ftdata->PartPreAceLin  ,external);
  BdPart->CreateArray("preaceang"  ,JBinaryDataDef::DatFloat3,scount,ftdata->PartPreAceAng  ,external);
  //-Crea arrays con datos de force points. Create arrays with force points data.
  if(ftdata->FptCount){
    const unsigned spcount=ftdata->FptCount*count;
    BdPart->CreateArray("FptMkbound",JBinaryDataDef::DatUshort ,spcount,ftdata->FptMkbound,external);
    BdPart->CreateArray("FptPos"    ,JBinaryDataDef::DatDouble3,spcount,ftdata->FptPos    ,external);
    BdPart->CreateArray("FptForce"  ,JBinaryDataDef::DatFloat3 ,spcount,ftdata->FptForce  ,external);
  }
  return(BdPart);
}

//==============================================================================
/// Graba BdPart data en fichero bi4.
/// Saves BdPart data in file bi4.
//==============================================================================
void JPartFloatInfoBi4Save::SaveBdPart(){
  if(!InitialSaved)SaveInitial();
  BdPart->SaveFileListApp(FileFull,"JPartFloatInfoBi4",true,true);
  BdPart->RemoveArrays();
  FtData->ClearData();
  PartCount++;
}

//==============================================================================
/// Graba datos almacenados en fichero bi4.
/// Saves stored data in file bi4.
//==============================================================================
void JPartFloatInfoBi4Save::SaveStoredData(){
  if(!InitialSaved)SaveInitial();
  if(FtData->PartCount){
    if(FtData->PartCount==1)MakeBdPartSingle(FtData->PartNum[0],FtData->PartTime[0],FtData->PartStep[0],FtData);
    else MakeBdPartArray(FtData);
    SaveBdPart();
  }
}

//==============================================================================
/// Graba particulas excluidas del PART.
/// Records particles excluded from PART.
//==============================================================================
void JPartFloatInfoBi4Save::AddDataPart(int cpart,unsigned step
  ,const JPartFloatInfoBi4Data* ftdata)
{
  //-Check step of ftdata.
  if(ftdata->PartSize!=1 || step!=ftdata->PartStep[0])
    Run_Exceptioon("Step number does not match the current data or data invalid.");
  //-Add data from ftdata.
  const unsigned count=FtData->PartCount;
  if(count && FtData->PartNum[0]!=cpart)SaveStoredData();
  FtData->AddData(ftdata);
}

//==============================================================================
/// Graba particulas excluidas del PART.
/// Records particles excluded from PART.
//==============================================================================
void JPartFloatInfoBi4Save::SavePart(int cpart,unsigned step
  ,const JPartFloatInfoBi4Data* ftdata)
{
  if(!InitialSaved)SaveInitial();
  //-Check step of ftdata.
  if(ftdata->PartSize!=1 || step!=ftdata->PartStep[0])
    Run_Exceptioon("Step number does not match the current data or data invalid.");
  //-Configures item Part.
  MakeBdPartSingle(cpart,ftdata->PartTime[0],step,ftdata);
  //-Save file and clear Part object.
  SaveBdPart();
}

//##############################################################################
//# JPartFloatInfoBi4Load
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartFloatInfoBi4Load::JPartFloatInfoBi4Load(){
  ClassName="JPartFloatInfoBi4Load";
  FtData=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartFloatInfoBi4Load::~JPartFloatInfoBi4Load(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartFloatInfoBi4Load::Reset(){
  FormatVer=0;
  FileData="";
  delete FtData; FtData=NULL;
  MkBoundFirst=0;
  FtCount=0;
  FptCount=0;
  PartCount=0;
  FirstPart=0;
}

//==============================================================================
/// Devuelve nombre de fichero bi4 por defecto.
/// Returns the filename bi4 by default.
//==============================================================================
std::string JPartFloatInfoBi4Load::GetFileNameDef(bool mainfile,std::string dir){
  if(!dir.empty())dir=fun::GetDirWithSlash(dir);
  return(dir+(mainfile? "PartFloatInfo.ibi4": "PartFloatInfo2.ibi4"));
}

//==============================================================================
/// Checks and returns data pointer of requested array.
//==============================================================================
const word* JPartFloatInfoBi4Load::GetArrayWord(JBinaryData* bd
  ,const std::string& name,unsigned count)
{
  JBinaryDataArray* ar=bd->GetArrayTpSize(name,JBinaryDataDef::DatUshort,count,FileData);
  return((const word*)ar->GetDataPointer() );
}

//==============================================================================
/// Checks and returns data pointer of requested array.
//==============================================================================
const unsigned* JPartFloatInfoBi4Load::GetArrayUint(JBinaryData* bd
  ,const std::string& name,unsigned count)
{
  JBinaryDataArray* ar=bd->GetArrayTpSize(name,JBinaryDataDef::DatUint,count,FileData);
  return((const unsigned*)ar->GetDataPointer() );
}

//==============================================================================
/// Checks and returns data pointer of requested array.
//==============================================================================
const float* JPartFloatInfoBi4Load::GetArrayFloat(JBinaryData* bd
  ,const std::string& name,unsigned count)
{
  JBinaryDataArray* ar=bd->GetArrayTpSize(name,JBinaryDataDef::DatFloat,count,FileData);
  return((const float*)ar->GetDataPointer() );
}

//==============================================================================
/// Checks and returns data pointer of requested array.
//==============================================================================
const double* JPartFloatInfoBi4Load::GetArrayDouble(JBinaryData* bd
  ,const std::string& name,unsigned count)
{
  JBinaryDataArray* ar=bd->GetArrayTpSize(name,JBinaryDataDef::DatDouble,count,FileData);
  return((const double*)ar->GetDataPointer() );
}

//==============================================================================
/// Checks and returns data pointer of requested array.
//==============================================================================
const tfloat3* JPartFloatInfoBi4Load::GetArrayFloat3(JBinaryData* bd
  ,const std::string& name,unsigned count)
{
  JBinaryDataArray* ar=bd->GetArrayTpSize(name,JBinaryDataDef::DatFloat3,count,FileData);
  return((const tfloat3*)ar->GetDataPointer() );
}

//==============================================================================
/// Checks and returns data pointer of requested array.
//==============================================================================
const tdouble3* JPartFloatInfoBi4Load::GetArrayDouble3(JBinaryData* bd
  ,const std::string& name,unsigned count)
{
  JBinaryDataArray* ar=bd->GetArrayTpSize(name,JBinaryDataDef::DatDouble3,count,FileData);
  return((const tdouble3*)ar->GetDataPointer() );
}

//==============================================================================
/// Carga datos de fichero y comprueba cabecera.
/// Loads data from file and verifies header.
//==============================================================================
void JPartFloatInfoBi4Load::LoadFile(std::string filename){
  Reset();
  FileData=filename;
  JBinaryData* bd=new JBinaryData("JPartFloatInfoBi4");
  bd->LoadFileListApp(FileData,"JPartFloatInfoBi4");
  //-Load header numbers.
  JBinaryData* head=bd->GetItem("LS0000_JPartFloatInfoBi4");
  if(!head)Run_ExceptioonFile("The head item is missing.",FileData);
  FormatVer=head->GetvUint("FormatVer",true,0);
  if(FormatVer<FormatVerDef)Run_ExceptioonFile(fun::PrintStr("The data format version \'%u\' is not valid. Version \'%u\' required.",FormatVer,FormatVerDef),FileData);
  MkBoundFirst=head->GetvUshort("MkBoundFirst");
  FtCount=head->GetvUint("FtCount");
  FptCount=head->GetvUint("FptCount");
  //-Create object for data.
  FtData=new JPartFloatInfoBi4Data(FtCount,FptCount,MkBoundFirst);
  //-Load header data.
  {
    const word*     pmkb =GetArrayWord (head,"MkBound",FtCount);
    const unsigned* pbeg =GetArrayUint (head,"Beginp",FtCount);
    const unsigned* pnp  =GetArrayUint (head,"Countp",FtCount);
    const float*    pmass=GetArrayFloat(head,"Mass",FtCount);
    const float*    massp=GetArrayFloat(head,"Massp",FtCount);
    const float*    pradi=GetArrayFloat(head,"Radius",FtCount);
    for(unsigned cf=0;cf<FtCount;cf++){
      FtData->SetCteData(cf,pmkb[cf],pbeg[cf],pnp[cf]
        ,pmass[cf],massp[cf],pradi[cf]);
    }
  }
  //-Check PARTs and total data items.
  const unsigned nitems=unsigned(bd->GetItemsCount());
  unsigned ndata=0;
  unsigned cpart0=0;
  for(unsigned cp=1;cp<nitems;cp++){
    JBinaryData* item=bd->GetItem(cp);
    const unsigned cpart=item->GetvUint("Cpart");
    if(cp>1 && cpart0+1!=cpart)Run_ExceptioonFile("Loaded data is corrupted. The data could have been modified by several simultaneous executions.",FileData);
    cpart0=cpart;
    //const string tname=fun::PrintStr("LS%04u_PART_%04u",cp,cp-1);
    //if(item->GetName()!=tname)Run_ExceptioonFile("Name of item is invalid.",FileData);
    ndata+=item->GetvUint("Count",true,1);
  }
  //-Loads PARTs with total data items.
  FtData->ResizeData(ndata);
  unsigned count1=0;
  for(unsigned cp=1;cp<nitems;cp++){
    JBinaryData* item=bd->GetItem(cp);
    const unsigned cpart=item->GetvUint("Cpart");
    const unsigned count2=item->GetvUint("Count",true,1);
    //-Data variables of PARTs.
    if(count2==1){
      FtData->PartNum [count1]=cpart;
      FtData->PartTime[count1]=item->GetvDouble("TimeStep");
      FtData->PartStep[count1]=item->GetvUint  ("Step");
    }
    else{
      for(unsigned c=0;c<count2;c++)FtData->PartNum[count1+c]=cpart;
      item->CopyArrayData("TimeStep",count2,FtData->PartTime+count1);
      item->CopyArrayData("Step"    ,count2,FtData->PartStep+count1);
    }
    //-Data variables of floating bodies.
    const unsigned fcount1=FtCount*count1;
    const unsigned fcount2=FtCount*count2;
    item->CopyArrayData("center"     ,fcount2,FtData->PartCenter     +fcount1);
    item->CopyArrayData("fvel"       ,fcount2,FtData->PartVelLin     +fcount1);
    item->CopyArrayData("fomega"     ,fcount2,FtData->PartVelAng     +fcount1);
    item->CopyArrayData("facelin"    ,fcount2,FtData->PartAceLin     +fcount1);
    item->CopyArrayData("faceang"    ,fcount2,FtData->PartAceAng     +fcount1);
    item->CopyArrayData("extforcelin",fcount2,FtData->PartExtForceLin+fcount1);
    item->CopyArrayData("extforceang",fcount2,FtData->PartExtForceAng+fcount1);
    item->CopyArrayData("fluforcelin",fcount2,FtData->PartFluForceLin+fcount1);
    item->CopyArrayData("fluforceang",fcount2,FtData->PartFluForceAng+fcount1);
    item->CopyArrayData("preacelin"  ,fcount2,FtData->PartPreAceLin  +fcount1);
    item->CopyArrayData("preaceang"  ,fcount2,FtData->PartPreAceAng  +fcount1);
    //-Data of force points.
    if(FptCount){
      const unsigned fpcount1=FptCount*count1;
      const unsigned fpcount2=FptCount*count2;
      item->CopyArrayData("FptMkbound",fpcount2,FtData->FptMkbound+fpcount1);
      item->CopyArrayData("FptPos"    ,fpcount2,FtData->FptPos    +fpcount1);
      item->CopyArrayData("FptForce"  ,fpcount2,FtData->FptForce  +fpcount1);
    }
    count1+=count2;
  }
  if(count1!=ndata)Run_ExceptioonFile("Allocated memory does not match with loaded data.",FileData);
  FtData->PartCount=ndata;
  PartCount=ndata;
  if(PartCount)FirstPart=FtData->PartNum[0];
  delete bd; bd=NULL;
}

//==============================================================================
/// Carga datos de fichero sin comprobar cabecera.
/// Load file data without checking header.
//==============================================================================
void JPartFloatInfoBi4Load::CheckHeadData(unsigned cf,word mkbound,unsigned begin
  ,unsigned count,float mass,float massp)
{
  if(cf>=FtCount)Run_Exceptioon("Number of floating is invalid.");
  const JPartFloatInfoBi4Data::StFtDataCte fth=FtData->GetHead(cf);
  if(fth.mkbound!=mkbound)Run_Exceptioon("The mkbound does not match.");
  if(fth.begin  !=begin)  Run_Exceptioon("The begin does not match.");
  if(fth.count  !=count)  Run_Exceptioon("The count does not match.");
  if(fth.mass   !=mass)   Run_Exceptioon("The mass does not match.");
  if(fth.massp  !=massp)  Run_Exceptioon("The massp does not match.");
}

