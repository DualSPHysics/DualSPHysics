//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JPartFloatBi4.cpp \brief Implements the classes JPartFloatBi4Save and class JPartFloatBi4Load.

#include "JPartFloatBi4.h"
#include "Functions.h"
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>

using namespace std;

//##############################################################################
//# JPartFloatBi4Save
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartFloatBi4Save::JPartFloatBi4Save(){
  ClassName="JPartFloatBi4Save";
  Data=NULL;
  HeadMkbound=NULL; HeadBegin=NULL; HeadCount=NULL; HeadMass=NULL; HeadRadius=NULL;
  PartCenter=NULL; PartFvel=NULL; PartFomega=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartFloatBi4Save::~JPartFloatBi4Save(){
  DestructorActive=true;
  Reset();
  delete Data; Data=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartFloatBi4Save::Reset(){
  ResizeFtData(0);
  ResetData();
  FormatVer=FormatVerDef;
  Dir="";
  MkBoundFirst=0;
  InitialSaved=false;
}

//==============================================================================
/// Elimina informacion de Data.
/// Deletes information from data. 
//==============================================================================
void JPartFloatBi4Save::ResetData(){
  delete Data; 
  Data=new JBinaryData("JPartFloatBi4");
  Part=Data->CreateItem("Part");
  Cpart=0;
}

//==============================================================================
/// Elimina informacion de PARTs.
/// Deletes information from PARTs.
//==============================================================================
void JPartFloatBi4Save::ResetPart(){
  Part->Clear();
}

//==============================================================================
/// Devuelve la memoria reservada.
/// Returns allocated memory
//==============================================================================
long long JPartFloatBi4Save::GetAllocMemory()const{  
  long long s=0;
  s+=Data->GetAllocMemory();
  if(HeadMkbound)s=s+sizeof(word)    *FtCount;
  if(HeadBegin)  s=s+sizeof(unsigned)*FtCount;
  if(HeadCount)  s=s+sizeof(unsigned)*FtCount;
  if(HeadMass)   s=s+sizeof(float)   *FtCount;
  if(HeadRadius) s=s+sizeof(float)   *FtCount;
  if(PartCenter) s=s+sizeof(tdouble3)*FtCount;
  if(PartFvel)   s=s+sizeof(tfloat3) *FtCount;
  if(PartFomega) s=s+sizeof(tfloat3) *FtCount;
  return(s);
}

//==============================================================================
/// Redimensiona memoria para datos de floatings.
/// Resize memory for floating data.
//==============================================================================
void JPartFloatBi4Save::ResizeFtData(unsigned ftcount){
  FtCount=ftcount;
  //-Libera memoria. Free memory
  delete[] HeadMkbound; HeadMkbound=NULL;
  delete[] HeadBegin;   HeadBegin=NULL;
  delete[] HeadCount;   HeadCount=NULL;
  delete[] HeadMass;    HeadMass=NULL;
  delete[] HeadRadius;  HeadRadius=NULL;
  delete[] PartCenter;  PartCenter=NULL;
  delete[] PartFvel;    PartFvel=NULL;
  delete[] PartFomega;  PartFomega=NULL;
  //-Asigna memoria. Assign memory.
  if(FtCount){
    HeadMkbound=new word    [FtCount];
    HeadBegin  =new unsigned[FtCount];
    HeadCount  =new unsigned[FtCount];
    HeadMass   =new float   [FtCount];
    HeadRadius =new float   [FtCount];
    PartCenter =new tdouble3[FtCount];
    PartFvel   =new tfloat3 [FtCount];
    PartFomega =new tfloat3 [FtCount];
    memset(HeadMkbound,0,sizeof(word)    *FtCount);
    memset(HeadBegin  ,0,sizeof(unsigned)*FtCount);
    memset(HeadCount  ,0,sizeof(unsigned)*FtCount);
    memset(HeadMass   ,0,sizeof(float)   *FtCount);
    memset(HeadRadius ,0,sizeof(float)   *FtCount);
    ClearPartData();
  }
}

//==============================================================================
/// Vacia datos de floatings por PART.
/// Clears data of floatings by PART.
//==============================================================================
void JPartFloatBi4Save::ClearPartData(){
  if(PartCenter)memset(PartCenter ,0,sizeof(tdouble3)*FtCount);
  if(PartFvel)  memset(PartFvel   ,0,sizeof(tfloat3) *FtCount);
  if(PartFomega)memset(PartFomega ,0,sizeof(tfloat3) *FtCount);
}

//==============================================================================
/// Devuelve nombre de part segun su numero.
/// Returns name of part according to their number.
//==============================================================================
std::string JPartFloatBi4Save::GetNamePart(unsigned cpart){
  return(fun::PrintStr("PART_%04u",cpart));
}

//==============================================================================
/// Devuelve nombre de fichero PART segun los parametros indicados.
/// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JPartFloatBi4Save::GetFileNamePart(){
  return("PartFloat.fbi4");
}

//==============================================================================
/// Configuracion de datos de cabecera.
/// Configuration of header data.
//==============================================================================
void JPartFloatBi4Save::Config(std::string appname,const std::string &dir,word mkboundfirst,unsigned ftcount){
  Reset();
  AppName=appname;
  Dir=fun::GetDirWithSlash(dir);
  MkBoundFirst=mkboundfirst;
  ResizeFtData(ftcount);
}

//==============================================================================
/// Añade datos de cabecera de floatings.
/// Adds data to floating header.
//==============================================================================
void JPartFloatBi4Save::AddHeadData(unsigned cf,word mkbound,unsigned begin,unsigned count,float mass,float radius){
  if(cf>=FtCount)RunException("AddHeadData","Number of floating is invalid.");
  HeadMkbound[cf]=mkbound;
  HeadBegin[cf]=begin;
  HeadCount[cf]=count;
  HeadMass[cf]=mass;
  HeadRadius[cf]=radius;
}

//==============================================================================
/// Grabacion inicial de fichero con info de Data.
/// Initial file recording with info from Data.
//==============================================================================
void JPartFloatBi4Save::SaveInitial(){
  if(!InitialSaved){
    Data->SetvText("AppName",AppName);
    Data->SetvUint("FormatVer",FormatVer);
    Data->SetvUshort("MkBoundFirst",MkBoundFirst);
    Data->SetvUint("FtCount",FtCount);
    Data->CreateArray("mkbound",JBinaryDataDef::DatUshort,FtCount,HeadMkbound,false);
    Data->CreateArray("begin",JBinaryDataDef::DatUint,FtCount,HeadBegin,false);
    Data->CreateArray("count",JBinaryDataDef::DatUint,FtCount,HeadCount,false);
    Data->CreateArray("mass",JBinaryDataDef::DatFloat,FtCount,HeadMass,false);
    Data->CreateArray("radius",JBinaryDataDef::DatFloat,FtCount,HeadRadius,false);
    Part->SetHide(true);
    Data->SaveFile(Dir+GetFileNamePart(),true,false);
    InitialSaved=true;
  }
}

//==============================================================================
/// Añade datos de particulas de de nuevo part.
/// Adds data of particles to new part.
//==============================================================================
void JPartFloatBi4Save::AddPartData(unsigned cf,const tdouble3 &center,const tfloat3 &fvel,const tfloat3 &fomega){
  if(cf>=FtCount)RunException("AddPartData","Number of floating is invalid.");
  PartCenter[cf]=center;
  PartFvel[cf]=fvel;
  PartFomega[cf]=fomega;
}


//==============================================================================
/// Añade datos de particulas de de nuevo part.
/// Adds data of particles to new part.
//==============================================================================
JBinaryData* JPartFloatBi4Save::AddPartFloat(unsigned cpart,double timestep,double demdtforce){
  const char met[]="AddPartFloat";
  //-Configura item Part. Configures item Part.
  Part->Clear();
  Cpart=cpart;
  Part->SetName(GetNamePart(cpart));
  Part->SetvUint("Cpart",cpart);
  Part->SetvDouble("TimeStep",timestep);
  Part->SetvDouble("DemDtForce",demdtforce);
  //-Crea array con datos de floatings. Create array with floatings data.
  Part->CreateArray("center",JBinaryDataDef::DatDouble3,FtCount,PartCenter,false);
  Part->CreateArray("fvel",JBinaryDataDef::DatFloat3,FtCount,PartFvel,false);
  Part->CreateArray("fomega",JBinaryDataDef::DatFloat3,FtCount,PartFomega,false);
  ClearPartData();
  return(Part);
}

//==============================================================================
/// Graba particulas excluidas del PART.
/// Records particles excluded from PART.
//==============================================================================
void JPartFloatBi4Save::SavePartFloat(){
  if(!InitialSaved)SaveInitial();
  Part->SaveFileListApp(Dir+GetFileNamePart(),"JPartFloatBi4",true,true);
  Part->RemoveArrays();
}



//##############################################################################
//# JPartFloatBi4Load
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartFloatBi4Load::JPartFloatBi4Load(){
  ClassName="JPartFloatBi4Load";
  Data=NULL;
  HeadMkbound=NULL; HeadBegin=NULL; HeadCount=NULL; HeadMass=NULL; HeadRadius=NULL;
  PartCenter=NULL; PartFvel=NULL; PartFomega=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartFloatBi4Load::~JPartFloatBi4Load(){
  DestructorActive=true;
  Reset();
  delete Data; Data=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartFloatBi4Load::Reset(){
  FormatVer=0;
  delete Data; 
  Data=new JBinaryData("JPartFloatBi4");
  MkBoundFirst=0;
  FtCount=PartCount=0;
  ResizeFtData(0);
  Part=NULL;
  ResetPart();
}

//==============================================================================
/// Elimina informacion de PART.
/// Removes information from PART.
//==============================================================================
void JPartFloatBi4Load::ResetPart(){
  Part=NULL;
  TimeStep=DemDtForce=0;
}

//==============================================================================
/// Redimensiona memoria para datos de floatings.
/// Resize memory for floating data.
//==============================================================================
void JPartFloatBi4Load::ResizeFtData(unsigned ftcount){
  FtCount=ftcount;
  //-Libera memoria. Frees memory
  delete[] HeadMkbound; HeadMkbound=NULL;
  delete[] HeadBegin;   HeadBegin=NULL;
  delete[] HeadCount;   HeadCount=NULL;
  delete[] HeadMass;    HeadMass=NULL;
  delete[] HeadRadius;  HeadRadius=NULL;
  delete[] PartCenter;  PartCenter=NULL;
  delete[] PartFvel;    PartFvel=NULL;
  delete[] PartFomega;  PartFomega=NULL;
  //-Asigna memoria. Asign memory
  if(FtCount){
    HeadMkbound=new word    [FtCount];
    HeadBegin  =new unsigned[FtCount];
    HeadCount  =new unsigned[FtCount];
    HeadMass   =new float   [FtCount];
    HeadRadius =new float   [FtCount];
    PartCenter =new tdouble3[FtCount];
    PartFvel   =new tfloat3 [FtCount];
    PartFomega =new tfloat3 [FtCount];
    memset(HeadMkbound,0,sizeof(word)    *FtCount);
    memset(HeadBegin  ,0,sizeof(unsigned)*FtCount);
    memset(HeadCount  ,0,sizeof(unsigned)*FtCount);
    memset(HeadMass   ,0,sizeof(float)   *FtCount);
    memset(HeadRadius ,0,sizeof(float)   *FtCount);
    memset(PartCenter ,0,sizeof(tdouble3)*FtCount);
    memset(PartFvel   ,0,sizeof(tfloat3) *FtCount);
    memset(PartFomega ,0,sizeof(tfloat3) *FtCount);
  }
}

//==============================================================================
// Devuelve nombre de fichero PART segun los parametros indicados.
// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JPartFloatBi4Load::GetFileNamePart(){
  return("PartFloat.fbi4");
}

//==============================================================================
/// Carga datos de fichero sin comprobar cabecera.
/// Load file data without checking header.
//==============================================================================
JBinaryDataArray* JPartFloatBi4Load::CheckArray(JBinaryData *bd,const std::string &name,JBinaryDataDef::TpData type){
  const char met[]="CheckArray";
  JBinaryDataArray *ar=bd->GetArray(name);
  if(!ar)RunException(met,string("The array ")+name+" is missing.");
  if(ar->GetType()!=type)RunException(met,string("The type of array ")+name+" does not match.");
  if(ar->GetCount()!=FtCount)RunException(met,string("The size of array ")+name+" does not match.");
  return(ar);
}

//==============================================================================
/// Carga datos de fichero y comprueba cabecera.
/// Loads data from file and verifies header.
//==============================================================================
void JPartFloatBi4Load::LoadFile(const std::string &dir){
  const char met[]="LoadFile";
  string file=fun::GetDirWithSlash(dir)+GetFileNamePart();
  Reset();
  Data->LoadFileListApp(file,"JPartFloatBi4");
  JBinaryData *head=Data->GetItem("LS0000_JPartFloatBi4");
  if(!head)RunException(met,"The head item is missing.",file);
  FormatVer=head->GetvUint("FormatVer",true,0);
  if(FormatVer<FormatVerDef)RunException(met,fun::PrintStr("The data format version \'%u\' is not valid. Version \'%u\' required.",FormatVer,FormatVerDef),file);
  MkBoundFirst=head->GetvUshort("MkBoundFirst",true,0);
  FtCount=head->GetvUint("FtCount",true,0);
  PartCount=Data->GetItemsCount()-1;
  FirstPart=(Data->GetItemsCount()>1? Data->GetItem(1)->GetvUint("Cpart",true,0): 0);
  //-Carga datos constantes de floatings (head). Load constant data of floatings (head).
  ResizeFtData(FtCount);
  {//-Loads array mkbound.
    JBinaryDataArray *ar=CheckArray(head,"mkbound",JBinaryDataDef::DatUshort);
    memcpy(HeadMkbound,(const word *)ar->GetDataPointer(),sizeof(word)*FtCount);
  }
  {//-Loads array begin.
    JBinaryDataArray *ar=CheckArray(head,"begin",JBinaryDataDef::DatUint);
    memcpy(HeadBegin,(const unsigned *)ar->GetDataPointer(),sizeof(unsigned)*FtCount);
  }
  {//-Loads array count.
    JBinaryDataArray *ar=CheckArray(head,"count",JBinaryDataDef::DatUint);
    memcpy(HeadCount,(const unsigned *)ar->GetDataPointer(),sizeof(unsigned)*FtCount);
  }
  {//-Loads array mass.
    JBinaryDataArray *ar=CheckArray(head,"mass",JBinaryDataDef::DatFloat);
    memcpy(HeadMass,(const float *)ar->GetDataPointer(),sizeof(float)*FtCount);
  }
  {//-Loads array radius.
    JBinaryDataArray *ar=CheckArray(head,"radius",JBinaryDataDef::DatFloat);
    memcpy(HeadRadius,(const float *)ar->GetDataPointer(),sizeof(float)*FtCount);
  }
  //head->SaveFileXml("probando.xml",true);
  //Data->SaveFileXml("probando.xml",true);
}

//==============================================================================
/// Carga datos de fichero sin comprobar cabecera.
/// Load file data without checking header.
//==============================================================================
void JPartFloatBi4Load::CheckHeadData(unsigned cf,word mkbound,unsigned begin,unsigned count,float mass){
  const char met[]="CheckHeadData";
  if(cf>=FtCount)RunException(met,"Number of floating is invalid.");
  if(HeadMkbound[cf]!=mkbound)RunException(met,"The mkbound does not match.");
  if(HeadBegin  [cf]!=begin)  RunException(met,"The begin does not match.");
  if(HeadCount  [cf]!=count)  RunException(met,"The count does not match.");
  if(HeadMass   [cf]!=mass)   RunException(met,"The mass does not match.");
}

//==============================================================================
/// Selecciona el PART indicado y devuelve false en caso de error.
/// Selects the indicated PART and returns false in case of error.
//==============================================================================
void JPartFloatBi4Load::LoadPart(unsigned cpart){
  const char met[]="LoadPart";
  ResetPart();
  if(!Data)RunException(met,"No loaded data.");
  string partname=fun::PrintStr("PART_%04u",cpart);
  unsigned spartname=unsigned(partname.size());
  const unsigned count=Data->GetItemsCount();
  for(unsigned c=1;c<Data->GetItemsCount() && !Part;c++){
    string name=Data->GetItem(c)->GetName();
    unsigned sname=unsigned(name.size());
    if(sname>spartname && name.substr(sname-spartname)==partname)Part=Data->GetItem(c);
  }
  if(Part){
    TimeStep=Part->GetvDouble("TimeStep");
    DemDtForce=Part->GetvDouble("DemDtForce");
    {//-Loads array center.
      JBinaryDataArray *ar=CheckArray(Part,"center",JBinaryDataDef::DatDouble3);
      memcpy(PartCenter,(const tdouble3 *)ar->GetDataPointer(),sizeof(tdouble3)*FtCount);
    }
    {//-Loads array fvel.
      JBinaryDataArray *ar=CheckArray(Part,"fvel",JBinaryDataDef::DatFloat3);
      memcpy(PartFvel,(const tfloat3 *)ar->GetDataPointer(),sizeof(tfloat3)*FtCount);
    }
    {//-Loads array fomega.
      JBinaryDataArray *ar=CheckArray(Part,"fomega",JBinaryDataDef::DatFloat3);
      memcpy(PartFomega,(const tfloat3 *)ar->GetDataPointer(),sizeof(tfloat3)*FtCount);
    }
  }
  else RunException(met,"PART not found.");
}

//==============================================================================
/// Selecciona el PART indicado y devuelve false en caso de error.
/// Select the indicated PART and returns false on error.
//==============================================================================
void JPartFloatBi4Load::CheckPart()const{
  if(!Part)RunException("CheckPart","PART not found.");
}

//==============================================================================
/// Comprueba que el numero de floating sea valido.
/// Verifies that the number of floating is valid.
//==============================================================================
void JPartFloatBi4Load::CheckFloating(unsigned cf)const{
  if(cf>=FtCount)RunException("CheckFloating","Number of floating is invalid.");
}


