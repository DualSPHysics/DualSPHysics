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

/// \file JPartFloatBi4.cpp \brief Implements the classes JPartFloatBi4Save and class JPartFloatBi4Load.

#include "JPartFloatBi4.h"
#include "Functions.h"
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <climits>
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
  HeadMkbound=NULL; HeadBegin=NULL; HeadCount=NULL; 
  HeadMass=NULL; HeadMassp=NULL; HeadRadius=NULL;
  PartCenter=NULL; PartPosRef=NULL; PartFVelLin=NULL; PartFVelAng=NULL;
  PartFAceLin=NULL; PartFAceAng=NULL;
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
  FileFull="";
  MkBoundFirst=0;
  PosRefData=false;
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
  if(HeadMassp)  s=s+sizeof(float)   *FtCount;
  if(HeadRadius) s=s+sizeof(float)   *FtCount;
  if(PartCenter) s=s+sizeof(tdouble3)*FtCount;
  if(PartPosRef) s=s+sizeof(tdouble3)*FtCount*3;
  if(PartFVelLin)s=s+sizeof(tfloat3) *FtCount;
  if(PartFVelAng)s=s+sizeof(tfloat3) *FtCount;
  if(PartFAceLin)s=s+sizeof(tfloat3) *FtCount;
  if(PartFAceAng)s=s+sizeof(tfloat3) *FtCount;
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
  delete[] HeadMassp;   HeadMassp=NULL;
  delete[] HeadRadius;  HeadRadius=NULL;
  delete[] PartCenter;  PartCenter=NULL;
  delete[] PartPosRef;  PartPosRef=NULL;
  delete[] PartFVelLin; PartFVelLin=NULL;
  delete[] PartFVelAng; PartFVelAng=NULL;
  delete[] PartFAceLin; PartFAceLin=NULL;
  delete[] PartFAceAng; PartFAceAng=NULL;
  //-Asigna memoria. Assign memory.
  if(FtCount){
    HeadMkbound=new word    [FtCount];
    HeadBegin  =new unsigned[FtCount];
    HeadCount  =new unsigned[FtCount];
    HeadMass   =new float   [FtCount];
    HeadMassp  =new float   [FtCount];
    HeadRadius =new float   [FtCount];
    PartCenter =new tdouble3[FtCount];
    if(PosRefData)PartPosRef=new tdouble3[FtCount*3];
    PartFVelLin=new tfloat3 [FtCount];
    PartFVelAng=new tfloat3 [FtCount];
    PartFAceLin=new tfloat3 [FtCount];
    PartFAceAng=new tfloat3 [FtCount];
    memset(HeadMkbound,0,sizeof(word)    *FtCount);
    memset(HeadBegin  ,0,sizeof(unsigned)*FtCount);
    memset(HeadCount  ,0,sizeof(unsigned)*FtCount);
    memset(HeadMass   ,0,sizeof(float)   *FtCount);
    memset(HeadMassp  ,0,sizeof(float)   *FtCount);
    memset(HeadRadius ,0,sizeof(float)   *FtCount);
    ClearPartData();
  }
}

//==============================================================================
/// Vacia datos de floatings por PART.
/// Clears data of floatings by PART.
//==============================================================================
void JPartFloatBi4Save::ClearPartData(){
  if(PartCenter )memset(PartCenter ,0,sizeof(tdouble3)*FtCount);
  if(PartPosRef )memset(PartPosRef ,0,sizeof(tdouble3)*FtCount*3);
  if(PartFVelLin)memset(PartFVelLin,0,sizeof(tfloat3) *FtCount);
  if(PartFVelAng)memset(PartFVelAng,0,sizeof(tfloat3) *FtCount);
  if(PartFAceLin)memset(PartFAceLin,0,sizeof(tfloat3) *FtCount);
  if(PartFAceAng)memset(PartFAceAng,0,sizeof(tfloat3) *FtCount);
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
void JPartFloatBi4Save::Config(std::string appname,const std::string &dir
  ,word mkboundfirst,unsigned ftcount,bool saveposref){
  Reset();
  AppName=appname;
  Dir=fun::GetDirWithSlash(dir);
  MkBoundFirst=mkboundfirst;
  PosRefData=saveposref;
  ResizeFtData(ftcount);
}

//==============================================================================
/// Anhade datos de cabecera de floatings.
/// Adds data to floating header.
//==============================================================================
void JPartFloatBi4Save::AddHeadData(unsigned cf,word mkbound,unsigned begin
  ,unsigned count,float mass,float massp,float radius)
{
  if(cf>=FtCount)Run_Exceptioon("Number of floating is invalid.");
  HeadMkbound[cf]=mkbound;
  HeadBegin  [cf]=begin;
  HeadCount  [cf]=count;
  HeadMass   [cf]=mass;
  HeadMassp  [cf]=massp;
  HeadRadius [cf]=radius;
}

//==============================================================================
/// Grabacion inicial de fichero con info de Data.
/// Initial file recording with info from Data.
//==============================================================================
void JPartFloatBi4Save::SaveInitial(std::string file){
  if(!InitialSaved){
    Data->SetvText("AppName",AppName);
    Data->SetvUint("FormatVer",FormatVer);
    Data->SetvUshort("MkBoundFirst",MkBoundFirst);
    Data->SetvBool("PosRefData",PosRefData);
    Data->SetvUint("FtCount",FtCount);
    Data->CreateArray("mkbound",JBinaryDataDef::DatUshort,FtCount,HeadMkbound,false);
    Data->CreateArray("begin"  ,JBinaryDataDef::DatUint  ,FtCount,HeadBegin  ,false);
    Data->CreateArray("count"  ,JBinaryDataDef::DatUint  ,FtCount,HeadCount  ,false);
    Data->CreateArray("mass"   ,JBinaryDataDef::DatFloat ,FtCount,HeadMass   ,false);
    Data->CreateArray("massp"  ,JBinaryDataDef::DatFloat ,FtCount,HeadMassp  ,false);
    Data->CreateArray("radius" ,JBinaryDataDef::DatFloat ,FtCount,HeadRadius ,false);
    Part->SetHide(true);
    if(file.empty())FileFull=Dir+GetFileNamePart();
    else FileFull=Dir+file;
    Data->SaveFile(FileFull,true,false);
    InitialSaved=true;
  }
}

//==============================================================================
/// Anhade datos de particulas de de nuevo part.
/// Adds data of particles to new part.
//==============================================================================
void JPartFloatBi4Save::AddPartData(unsigned cf,const tdouble3 &center
  ,const tfloat3 &fvellin,const tfloat3 &fvelang
  ,const tfloat3 &facelin,const tfloat3 &faceang)
{
  if(cf>=FtCount)Run_Exceptioon("Number of floating is invalid.");
  PartCenter[cf]=center;
  PartFVelLin[cf]=fvellin;
  PartFVelAng[cf]=fvelang;
  PartFAceLin[cf]=facelin;
  PartFAceAng[cf]=faceang;
}

//==============================================================================
/// Anhade datos de posicion de particulas de de nuevo part.
/// Adds postion reference data of particles to new part.
//==============================================================================
void JPartFloatBi4Save::AddPartDataPosRef(unsigned ftcount,const tdouble3 *posref)
{
  if(ftcount!=FtCount)Run_Exceptioon("Number of floating is invalid.");
  memcpy(PartPosRef,posref,sizeof(tdouble3)*ftcount*3);
}


//==============================================================================
/// Anhade datos de particulas de de nuevo part.
/// Adds data of particles to new part.
//==============================================================================
JBinaryData* JPartFloatBi4Save::AddPartFloat(unsigned cpart,unsigned step
  ,double timestep,double demdtforce)
{
  //-Configura item Part. Configures item Part.
  Part->Clear();
  Cpart=cpart;
  Part->SetName(GetNamePart(cpart));
  Part->SetvUint("Cpart",cpart);
  Part->SetvUint("Step",step);
  Part->SetvDouble("TimeStep",timestep);
  Part->SetvDouble("DemDtForce",demdtforce);
  //-Crea array con datos de floatings. Create array with floatings data.
  Part->CreateArray("center" ,JBinaryDataDef::DatDouble3,FtCount,PartCenter,false);
  Part->CreateArray("fvel"   ,JBinaryDataDef::DatFloat3 ,FtCount,PartFVelLin,false);
  Part->CreateArray("fomega" ,JBinaryDataDef::DatFloat3 ,FtCount,PartFVelAng,false);
  Part->CreateArray("facelin",JBinaryDataDef::DatFloat3 ,FtCount,PartFAceLin,false);
  Part->CreateArray("faceang",JBinaryDataDef::DatFloat3 ,FtCount,PartFAceAng,false);
  if(PartPosRef)Part->CreateArray("posref" ,JBinaryDataDef::DatDouble3,FtCount*3,PartPosRef,false);
  ClearPartData();
  return(Part);
}

//==============================================================================
/// Graba particulas excluidas del PART.
/// Records particles excluded from PART.
//==============================================================================
void JPartFloatBi4Save::SavePartFloat(){
  if(!InitialSaved)SaveInitial();
  Part->SaveFileListApp(FileFull,"JPartFloatBi4",true,true);
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
  HeadMkbound=NULL; HeadBegin=NULL; HeadCount=NULL; 
  HeadMass=NULL; HeadMassp=NULL; HeadRadius=NULL;
  PartCenter=NULL; PartPosRef=NULL;
  PartFVelLin=NULL; PartFVelAng=NULL;
  PartFAceLin=NULL; PartFAceAng=NULL;
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
  FileData="";
  delete Data; 
  Data=new JBinaryData("JPartFloatBi4");
  MkBoundFirst=0;
  PosRefData=false;
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
  Step=UINT_MAX;
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
  delete[] HeadMassp;   HeadMassp=NULL;
  delete[] HeadRadius;  HeadRadius=NULL;
  delete[] PartCenter;  PartCenter=NULL;
  delete[] PartPosRef;  PartPosRef=NULL;
  delete[] PartFVelLin; PartFVelLin=NULL;
  delete[] PartFVelAng; PartFVelAng=NULL;
  delete[] PartFAceLin; PartFAceLin=NULL;
  delete[] PartFAceAng; PartFAceAng=NULL;
  //-Asigna memoria. Asign memory
  if(FtCount){
    HeadMkbound=new word    [FtCount];
    HeadBegin  =new unsigned[FtCount];
    HeadCount  =new unsigned[FtCount];
    HeadMass   =new float   [FtCount];
    HeadMassp  =new float   [FtCount];
    HeadRadius =new float   [FtCount];
    PartCenter =new tdouble3[FtCount];
    if(PosRefData)PartPosRef=new tdouble3[FtCount*3];
    PartFVelLin=new tfloat3 [FtCount];
    PartFVelAng=new tfloat3 [FtCount];
    PartFAceLin=new tfloat3 [FtCount];
    PartFAceAng=new tfloat3 [FtCount];
    memset(HeadMkbound,0,sizeof(word)    *FtCount);
    memset(HeadBegin  ,0,sizeof(unsigned)*FtCount);
    memset(HeadCount  ,0,sizeof(unsigned)*FtCount);
    memset(HeadMass   ,0,sizeof(float)   *FtCount);
    memset(HeadMassp  ,0,sizeof(float)   *FtCount);
    memset(HeadRadius ,0,sizeof(float)   *FtCount);
    memset(PartCenter ,0,sizeof(tdouble3)*FtCount);
    if(PartPosRef)memset(PartPosRef,0,sizeof(tdouble3)*FtCount*3);
    memset(PartFVelLin,0,sizeof(tfloat3) *FtCount);
    memset(PartFVelAng,0,sizeof(tfloat3) *FtCount);
    memset(PartFAceLin,0,sizeof(tfloat3) *FtCount);
    memset(PartFAceAng,0,sizeof(tfloat3) *FtCount);
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
JBinaryDataArray* JPartFloatBi4Load::CheckArray(JBinaryData *bd,const std::string &name
  ,JBinaryDataDef::TpData type)
{
  JBinaryDataArray *ar=bd->GetArray(name);
  if(!ar)Run_Exceptioon(string("The array ")+name+" is missing.");
  if(ar->GetType()!=type)Run_Exceptioon(string("The type of array ")+name+" does not match.");
  if(ar->GetCount()!=FtCount)Run_Exceptioon(string("The size of array ")+name+" does not match.");
  return(ar);
}

//==============================================================================
/// Comprueba lista de PARTs cargados.
/// Check list of loaded PARTs.
//==============================================================================
void JPartFloatBi4Load::CheckPartList()const{
  unsigned nitem=Data->GetItemsCount();
  if(nitem>1){
    unsigned cpart=Data->GetItem(1)->GetvUint("Cpart");
    for(unsigned c=2;c<nitem;c++){
      const unsigned cpart2=Data->GetItem(c)->GetvUint("Cpart");
      if(cpart2!=cpart+1)Run_ExceptioonFile("Loaded data is corrupted. The data could have been modified by several simultaneous executions.",FileData);
      cpart=cpart2;
    }
  }
}

//==============================================================================
/// Carga datos de fichero y comprueba cabecera.
/// Loads data from file and verifies header.
//==============================================================================
void JPartFloatBi4Load::LoadFile(const std::string &dir){
  Reset();
  FileData=fun::GetDirWithSlash(dir)+GetFileNamePart();
  Data->LoadFileListApp(FileData,"JPartFloatBi4");
  JBinaryData *head=Data->GetItem("LS0000_JPartFloatBi4");
  if(!head)Run_ExceptioonFile("The head item is missing.",FileData);
  FormatVer=head->GetvUint("FormatVer",true,0);
  if(FormatVer<FormatVerDef)Run_ExceptioonFile(fun::PrintStr("The data format version \'%u\' is not valid. Version \'%u\' required.",FormatVer,FormatVerDef),FileData);
  MkBoundFirst=head->GetvUshort("MkBoundFirst",true,0);
  PosRefData=head->GetvBool("PosRefData",true,false);
  FtCount=head->GetvUint("FtCount",true,0);
  PartCount=Data->GetItemsCount()-1;
  FirstPart=(PartCount? Data->GetItem(1)->GetvUint("Cpart",true,0): 0);
  CheckPartList();

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
  {//-Loads array massp.
    if(head->GetArray("massp")!=NULL){
      JBinaryDataArray *ar=CheckArray(head,"massp",JBinaryDataDef::DatFloat);
      memcpy(HeadMassp,(const float *)ar->GetDataPointer(),sizeof(float)*FtCount);
    }
    else{
      for(unsigned cf=0;cf<FtCount;cf++)HeadMassp[cf]=HeadMass[cf]/HeadCount[cf];
    }
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
void JPartFloatBi4Load::CheckHeadData(unsigned cf,word mkbound,unsigned begin
  ,unsigned count,float mass,float massp)
{
  if(cf>=FtCount)Run_Exceptioon("Number of floating is invalid.");
  if(HeadMkbound[cf]!=mkbound)Run_Exceptioon("The mkbound does not match.");
  if(HeadBegin  [cf]!=begin)  Run_Exceptioon("The begin does not match.");
  if(HeadCount  [cf]!=count)  Run_Exceptioon("The count does not match.");
  if(HeadMass   [cf]!=mass)   Run_Exceptioon("The mass does not match.");
  if(HeadMassp  [cf]!=massp)  Run_Exceptioon("The massp does not match.");
}

//==============================================================================
/// Selecciona el PART indicado y devuelve false en caso de error.
/// Selects the indicated PART and returns false in case of error.
//==============================================================================
void JPartFloatBi4Load::LoadPart(unsigned cpart){
  ResetPart();
  if(!Data)Run_Exceptioon("No loaded data.");
  string partname=fun::PrintStr("PART_%04u",cpart);
  unsigned spartname=unsigned(partname.size());
  const unsigned count=Data->GetItemsCount();
  for(unsigned c=1;c<Data->GetItemsCount() && !Part;c++){
    string name=Data->GetItem(c)->GetName();
    unsigned sname=unsigned(name.size());
    if(sname>spartname && name.substr(sname-spartname)==partname)Part=Data->GetItem(c);
  }
  if(Part){
    Step=Part->GetvUint("Step",true,UINT_MAX);
    TimeStep=Part->GetvDouble("TimeStep");
    DemDtForce=Part->GetvDouble("DemDtForce");
    {//-Loads array center.
      JBinaryDataArray *ar=CheckArray(Part,"center",JBinaryDataDef::DatDouble3);
      memcpy(PartCenter,(const tdouble3 *)ar->GetDataPointer(),sizeof(tdouble3)*FtCount);
    }
    {//-Loads array fvel.
      JBinaryDataArray *ar=CheckArray(Part,"fvel",JBinaryDataDef::DatFloat3);
      memcpy(PartFVelLin,(const tfloat3 *)ar->GetDataPointer(),sizeof(tfloat3)*FtCount);
    }
    {//-Loads array fomega.
      JBinaryDataArray *ar=CheckArray(Part,"fomega",JBinaryDataDef::DatFloat3);
      memcpy(PartFVelAng,(const tfloat3 *)ar->GetDataPointer(),sizeof(tfloat3)*FtCount);
    }
    FAceData=(Part->GetArray("facelin") && Part->GetArray("faceang"));
    if(FAceData){
      {//-Loads array facelin when it is available.
        JBinaryDataArray *ar=CheckArray(Part,"facelin",JBinaryDataDef::DatFloat3);
        memcpy(PartFAceLin,(const tfloat3 *)ar->GetDataPointer(),sizeof(tfloat3)*FtCount);
      }
      {//-Loads array faceang when it is available.
        JBinaryDataArray *ar=CheckArray(Part,"faceang",JBinaryDataDef::DatFloat3);
        memcpy(PartFAceAng,(const tfloat3 *)ar->GetDataPointer(),sizeof(tfloat3)*FtCount);
      }
    }
  }
  else Run_Exceptioon("PART not found.");
}

//==============================================================================
/// Selecciona el PART indicado y devuelve false en caso de error.
/// Select the indicated PART and returns false on error.
//==============================================================================
void JPartFloatBi4Load::CheckPart()const{
  if(!Part)Run_Exceptioon("PART not found.");
}

//==============================================================================
/// Comprueba que el numero de floating sea valido.
/// Verifies that the number of floating is valid.
//==============================================================================
void JPartFloatBi4Load::CheckFloating(unsigned cf)const{
  if(cf>=FtCount)Run_Exceptioon("Number of floating is invalid.");
}


