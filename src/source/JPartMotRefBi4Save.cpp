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

/// \file JPartMotRefBi4Save.cpp \brief Implements the class JPartMotRefBi4Save.

#include "JPartMotRefBi4Save.h"
#include "Functions.h"
#include <cstring>

using namespace std;

//##############################################################################
//# JPartMotRefBi4Save
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartMotRefBi4Save::JPartMotRefBi4Save(std::string appname,std::string dir)
  :AppName(appname),Dir(fun::GetDirWithSlash(dir))
{
  ClassName="JPartMotRefBi4Save";
  BdData=NULL; BdPart=NULL;
  DatTime=NULL;
  DatStep=NULL;
  DatPosRef=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartMotRefBi4Save::~JPartMotRefBi4Save(){
  DestructorActive=true;
  Reset();
  delete BdData; BdData=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartMotRefBi4Save::Reset(){
  FormatVer=FormatVerDef;
  MainFile=false;
  TimeOut=-1;
  FileFull="";
  MkBoundFirst=0;
  MkMovingCount=MkFloatCount=0;
  MkCount=0;
  PsCount=0;
  InitialSaved=false;
  PartCount=0;
  ResetBdData();
  ResizeDat(0);
  DatPart=-1;
}

//==============================================================================
/// Elimina informacion de Data.
/// Deletes information from data. 
//==============================================================================
void JPartMotRefBi4Save::ResetBdData(){
  delete BdData; 
  BdData=new JBinaryData("JPartMotRefBi4");
  BdPart=BdData->CreateItem("Part");
  PartCount=0;
}

//==============================================================================
/// Elimina informacion de PARTs.
/// Deletes information from PARTs.
//==============================================================================
void JPartMotRefBi4Save::ResetBdPart(){
  BdPart->Clear();
}

//==============================================================================
/// Resizes memory space for data items.
//==============================================================================
void JPartMotRefBi4Save::ResizeDat(unsigned size){
  //-Free memory.
  if(!size){
    DatSize=DatCount=0;
    delete[] DatTime;   DatTime=NULL;
    delete[] DatStep;   DatStep=NULL;
    delete[] DatPosRef; DatPosRef=NULL;
  }
  else{
    if(DatCount>size)Run_Exceptioon("No memory available for stored data.");
    DatSize=size;
    try{
      DatTime  =fun::ResizeAlloc(DatTime  ,DatCount,DatSize);
      DatStep  =fun::ResizeAlloc(DatStep  ,DatCount,DatSize);
      DatPosRef=fun::ResizeAlloc(DatPosRef,PsCount*DatCount,PsCount*DatSize);
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
  }
}

//==============================================================================
/// Devuelve la memoria reservada.
/// Returns allocated memory
//==============================================================================
llong JPartMotRefBi4Save::GetAllocMemory()const{  
  llong s=0;
  s+=BdData->GetAllocMemory();
  s+=sizeof(double)  *DatSize;
  s+=sizeof(unsigned)*DatSize;
  s+=sizeof(tdouble3)*DatSize*PsCount;
  return(s);
}

//==============================================================================
/// Devuelve nombre de part segun su numero.
/// Returns name of part according to their number.
//==============================================================================
std::string JPartMotRefBi4Save::GetNamePart(unsigned cpart){
  return(fun::PrintStr("PART_%04u",cpart));
}

//==============================================================================
/// Devuelve nombre de fichero PART segun los parametros indicados.
/// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JPartMotRefBi4Save::GetFileNameDef(bool mainfile,std::string dir){
  if(!dir.empty())dir=fun::GetDirWithSlash(dir);
  return(dir+(mainfile? "PartMotionRef.ibi4": "PartMotionRef2.ibi4"));
}

//==============================================================================
/// Configuracion de datos de cabecera.
/// Configuration of header data.
//==============================================================================
void JPartMotRefBi4Save::Config(bool mainfile,double timeout,word mkboundfirst
  ,unsigned mvcount,unsigned ftcount,const StMkMotionData* mkmotiondata)
{
  Reset();
  MainFile=mainfile;
  TimeOut=timeout;
  FileFull=Dir+GetFileNameDef(MainFile);
  MkBoundFirst=mkboundfirst;
  MkMovingCount=mvcount;
  MkFloatCount=ftcount;
  MkCount=MkMovingCount+MkFloatCount;
  if(!MkCount)Run_Exceptioon("No moving or floating bodies.");
  PsCount=MkCount*3;
  //-Save head data in BdData.
  BdData->SetvText("AppName",AppName);
  BdData->SetvUint("FormatVer",FormatVer);
  BdData->SetvBool("MainFile",MainFile);
  BdData->SetvDouble("TimeOut",TimeOut);
  BdData->SetvUshort("MkBoundFirst",MkBoundFirst);
  BdData->SetvUint("MkMovingCount",MkMovingCount);
  BdData->SetvUint("MkFloatCount",MkFloatCount);
  BdData->SetvUint("MkCount",MkCount);
  //-Creates head arrays.
  word*     pmkb=BdData->CreateArrayUshort ("MkBound",MkCount,true);
  unsigned* pnid=BdData->CreateArrayUint   ("Nid",MkCount,true);
  unsigned* pid =BdData->CreateArrayUint   ("Id" ,MkCount*3,true);
  tdouble3* pps =BdData->CreateArrayDouble3("Ps" ,MkCount*3,true);
  double*   pdis=BdData->CreateArrayDouble ("Dis",MkCount*3,true);
  //-Loads head arrays.
  for(unsigned cmk=0;cmk<MkCount;cmk++){
    const StMkMotionData& v=mkmotiondata[cmk];
	pmkb[cmk]=v.mkbound;
	pnid[cmk]=v.nid;
    for(unsigned ci=0;ci<3;ci++){
	  const unsigned cmki=cmk*3+ci;
	  pid [cmki]=v.id [ci];
	  pps [cmki]=v.ps [ci];
	  pdis[cmki]=v.dis[ci];
    }
  }
  BdPart->SetHide(true);
  //-Allocates memory for data.
  ResizeDat((MainFile? 0: 10));
}

//==============================================================================
/// Grabacion inicial de fichero con info de Data.
/// Initial file recording with info from Data.
//==============================================================================
void JPartMotRefBi4Save::SaveInitial(){
  if(!InitialSaved){
    BdData->SaveFile(FileFull,true,false);
    InitialSaved=true;
  }
}

//==============================================================================
/// Adds data to BdPart object.
//==============================================================================
JBinaryData* JPartMotRefBi4Save::MakeBdPartSingle(int cpart,double timestep
  ,unsigned step,unsigned npos,const tdouble3* posref)
{
  //-Configures item Part.
  BdPart->Clear();
  BdPart->SetName(GetNamePart(cpart));
  BdPart->SetvUint("Cpart",cpart);
  BdPart->SetvDouble("TimeStep",timestep);
  BdPart->SetvUint("Step",step);
  BdPart->CreateArray("PosRef",JBinaryDataDef::DatDouble3,npos,posref,false);
  return(BdPart);
}

//==============================================================================
/// Adds data to BdPart object.
//==============================================================================
JBinaryData* JPartMotRefBi4Save::MakeBdPartArray(int cpart){
  //-Configures item Part.
  BdPart->Clear();
  BdPart->SetName(GetNamePart(cpart));
  BdPart->SetvUint("Cpart",cpart);
  BdPart->SetvUint("Count",DatCount);
  const bool external=true;
  BdPart->CreateArray("TimeStep",JBinaryDataDef::DatDouble,DatCount,DatTime,external);
  BdPart->CreateArray("Step"    ,JBinaryDataDef::DatUint,DatCount,DatStep,external);
  BdPart->CreateArray("PosRef"  ,JBinaryDataDef::DatDouble3,DatCount*PsCount,DatPosRef,external);
  return(BdPart);
}

//==============================================================================
/// Graba BdPart data en fichero bi4.
/// Saves BdPart data in file bi4.
//==============================================================================
void JPartMotRefBi4Save::SaveBdPart(){
  if(!InitialSaved)SaveInitial();
  BdPart->SaveFileListApp(FileFull,"JPartMotRefBi4",true,true);
  BdPart->RemoveArrays();
  PartCount++;
}

//==============================================================================
/// Graba datos almacenados en fichero bi4.
/// Saves stored data in file bi4.
//==============================================================================
void JPartMotRefBi4Save::SaveStoredData(){
  if(!InitialSaved)SaveInitial();
  if(DatCount){
    if(DatCount==1)MakeBdPartSingle(DatPart,DatTime[0],DatStep[0],PsCount,DatPosRef);
    else MakeBdPartArray(DatPart);
    SaveBdPart();
  }
  DatCount=0; DatPart=-1;
}

//==============================================================================
/// Graba posiciones de referencia del PART en fichero.
/// Saves reference positions from PART in file.
//==============================================================================
void JPartMotRefBi4Save::AddDataPart(int cpart,double timestep,unsigned step
  ,unsigned npos,const tdouble3* posref)
{   
  if(npos!=PsCount)Run_Exceptioon("Size of data does not match.");
  if(DatCount && DatPart!=cpart)SaveStoredData();
  if(DatCount>=DatSize)ResizeDat(max(DatSize*2,10u));
  DatTime  [DatCount]=timestep;
  DatStep  [DatCount]=step;
  memcpy(DatPosRef+(PsCount*DatCount),posref,sizeof(tdouble3)*PsCount);
  DatPart=cpart;
  DatCount++;
}

//==============================================================================
/// Graba posiciones de referencia del PART en fichero.
/// Saves reference positions from PART in file.
//==============================================================================
void JPartMotRefBi4Save::SavePart(int cpart,double timestep,unsigned step
  ,unsigned npos,const tdouble3* posref)
{   
  if(DatCount)Run_Exceptioon("There is data to be saved.");
  MakeBdPartSingle(cpart,timestep,step,npos,posref);
  SaveBdPart();
}


