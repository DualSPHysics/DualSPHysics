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

/// \file JPartOutBi4Save.cpp \brief Implements the class \ref JPartOutBi4Save.

#include "JPartOutBi4Save.h"
#include "Functions.h"
#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

using namespace std;

//##############################################################################
//# JPartOutBi4Save
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartOutBi4Save::JPartOutBi4Save(){
  ClassName="JPartOutBi4Save";
  Data=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartOutBi4Save::~JPartOutBi4Save(){
  DestructorActive=true;
  delete Data; Data=NULL;
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartOutBi4Save::Reset(){
  ResetData();
  Dir="";
  Block=0;
  Piece=0;
  Npiece=1;
  BlockNoutMin=BLOCKNOUTMIN;
  BlockNoutMax=BLOCKNOUTMAX;
  InitialSaved=false;
  BlockNout=0;
}

//==============================================================================
/// Elimina informacion de Data.
//==============================================================================
void JPartOutBi4Save::ResetData(){
  delete Data; 
  Data=new JBinaryData("JPartOutBi4");
  Part=Data->CreateItem("Part");
  Cpart=0;
}

//==============================================================================
/// Elimina informacion de PARTs.
/// Deletes information from PARTs.
//==============================================================================
void JPartOutBi4Save::ResetPart(){
  Part->Clear();
}

//==============================================================================
/// Devuelve la memoria reservada.
/// Returns allocated memory
//==============================================================================
long long JPartOutBi4Save::GetAllocMemory()const{  
  long long s=0;
  s+=Data->GetAllocMemory();
  return(s);
}

//==============================================================================
/// Devuelve nombre de fichero PART segun los parametros indicados.
/// Returns the filename PART according to the specified parameters.
//==============================================================================
std::string JPartOutBi4Save::GetFileNamePart(unsigned block,unsigned piece,unsigned npiece){
  string fname="PartOut";
  char cad[32];
  if(npiece>1){
    sprintf(cad,"_p%02d",piece);
    fname=fname+cad;
  }
  sprintf(cad,"_%03u.obi4",block);
  return(fname+cad);
}

//==============================================================================
/// Configuracion de variables basicas.
/// Configuration of basic variables.
//==============================================================================
void JPartOutBi4Save::ConfigBasic(unsigned piece,unsigned npiece,std::string runcode,std::string appname,bool data2d,const std::string &dir){
  Reset();
  Piece=piece; Npiece=npiece;
  Dir=fun::GetDirWithSlash(dir);
  Data->SetvUint("Piece",Piece);
  Data->SetvUint("Npiece",Npiece);
  Data->SetvText("RunCode",runcode);
  Data->SetvText("Date",fun::GetDateTime());
  Data->SetvText("AppName",appname);
  Data->SetvBool("Data2d",data2d);
  Data->SetvUint("FmtVersion",FmtVersion);
  Data->SetvUint("Block",Block);
  ConfigLimits(TDouble3(0),TDouble3(0),0,0);
}

//==============================================================================
/// Configuracion de numero de particulas.
/// Configuration of number of particles.
//==============================================================================
void JPartOutBi4Save::ConfigParticles(ullong casenp,ullong casenfixed,ullong casenmoving,ullong casenfloat,ullong casenfluid){
  if(casenp!=casenfixed+casenmoving+casenfloat+casenfluid)Run_Exceptioon("Error in the number of particles.");
  Data->SetvUllong("CaseNp",casenp);
  Data->SetvUllong("CaseNfixed",casenfixed);
  Data->SetvUllong("CaseNmoving",casenmoving);
  Data->SetvUllong("CaseNfloat",casenfloat);
  Data->SetvUllong("CaseNfluid",casenfluid);
}

//==============================================================================
/// Configuracion de limites para exclusion de particulas.
/// Setting of limits for exclusion of particles.
//==============================================================================
void JPartOutBi4Save::ConfigLimits(const tdouble3 &mapposmin,const tdouble3 &mapposmax,float rhopmin,float rhopmax){
  Data->SetvDouble3("MapPosMin",mapposmin);
  Data->SetvDouble3("MapPosMax",mapposmax);
  Data->SetvFloat("RhopMin",rhopmin);
  Data->SetvFloat("RhopMax",rhopmax);
}

//==============================================================================
/// Grabacion inicial de fichero con info de Data.
/// Initial recording of file with Data info.
//==============================================================================
void JPartOutBi4Save::SaveInitial(){
  Part->SetHide(true);
  Data->SetvUint("Block",Block);
  Data->SaveFile(Dir+GetFileNamePart(Block,Piece,Npiece),true,false);
  InitialSaved=true;
}

//==============================================================================
/// Devuelve nombre de part segun su numero.
/// Returns name of part according to their number.
//==============================================================================
std::string JPartOutBi4Save::GetNamePart(unsigned cpart){
  char cad[64];
  sprintf(cad,"PART_%04u",cpart);
  return(cad);
}

//==============================================================================
/// Anhade datos de particulas de de nuevo part.
/// Adds data of particles to new part.
//==============================================================================
JBinaryData* JPartOutBi4Save::AddPartOut(unsigned cpart,double timestep,unsigned nout
  ,const unsigned *idp,const ullong *idpd,const tfloat3 *pos,const tdouble3 *posd
  ,const tfloat3 *vel,const float *rhop,const byte *motive)
{
  if(!idp && !idpd)Run_Exceptioon("The id of particles is invalid.");
  if(!pos && !posd)Run_Exceptioon("The position of particles is invalid.");
  //-Configura item Part. Configures item Part.
  Part->Clear();
  Cpart=cpart;
  Part->SetName(GetNamePart(cpart));
  Part->SetvUint("Cpart",cpart);
  Part->SetvDouble("TimeStep",timestep);
  Part->SetvUint("Nout",nout);
  //-Crea array con particulas excluidas. Creates array with excluded particles.
  if(idpd)Part->CreateArray("Idpd",JBinaryDataDef::DatUllong,nout,idpd,true);
  else    Part->CreateArray("Idp" ,JBinaryDataDef::DatUint,nout,idp,true);
  if(posd)Part->CreateArray("Posd",JBinaryDataDef::DatDouble3,nout,posd,true);
  else    Part->CreateArray("Pos" ,JBinaryDataDef::DatFloat3,nout,pos,true);
  Part->CreateArray("Vel",JBinaryDataDef::DatFloat3,nout,vel,true);
  Part->CreateArray("Rhop",JBinaryDataDef::DatFloat,nout,rhop,true);
  Part->CreateArray("Motive",JBinaryDataDef::DatUchar,nout,motive,true);
  return(Part);
}

//==============================================================================
/// Graba particulas excluidas del PART.
/// Records particles excluded from the PART.
//==============================================================================
void JPartOutBi4Save::SavePartOut(){
  if(!InitialSaved)SaveInitial();
  unsigned nout=Part->GetvUint("Nout");
  if(nout){
    if(BlockNout>=BlockNoutMin && BlockNout+nout>BlockNoutMax){//-Cambio de bloque, graba en otro fichero. Change of block, writes in another file.
      BlockNout=0;
      Block++;
      SaveInitial();
    }
    Part->SaveFileListApp(Dir+GetFileNamePart(Block,Piece,Npiece),"JPartOutBi4",true,true);
    BlockNout+=nout;
    Part->RemoveArrays();
  }
}

//==============================================================================
/// Graba particulas excluidas del PART.
/// Records particles excluded from the PART.
//==============================================================================
void JPartOutBi4Save::SavePartOut(bool posdouble,unsigned cpart,double timestep,unsigned nout
  ,const unsigned *idp,const tfloat3 *posf,const tdouble3 *posd,const tfloat3 *vel
  ,const float *rhop,const byte *motive)
{
  if(!posf && !posd)Run_Exceptioon("The position of particles is invalid.");
  if(posdouble){
    if(posd==NULL){
      tdouble3 *xpos=new tdouble3[nout];
      for(unsigned c=0;c<nout;c++)xpos[c]=ToTDouble3(posf[c]);
      SavePartOut(cpart,timestep,nout,idp,xpos,vel,rhop,motive);
      delete[] xpos; xpos=NULL;
    }
    else SavePartOut(cpart,timestep,nout,idp,posd,vel,rhop,motive);
  }
  else{
    if(posf==NULL){
      tfloat3 *xpos=new tfloat3[nout];
      for(unsigned c=0;c<nout;c++)xpos[c]=ToTFloat3(posd[c]);
      SavePartOut(cpart,timestep,nout,idp,xpos,vel,rhop,motive);
      delete[] xpos; xpos=NULL;
    }
    else SavePartOut(cpart,timestep,nout,idp,posf,vel,rhop,motive);
  }
}



