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

#include "JDebugSphCpu.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JException.h"
#include "Functions.h"
#include "JSpVtkData.h"
#include "JOutputCsv.h"
#include "JDataArrays.h"

#include "JSphCpuSingle.h"

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

using namespace std;

//##############################################################################
//# JDebugSphCpu
//##############################################################################
//==============================================================================
/// Throws exception related to a file from a static method.
//==============================================================================
void JDebugSphCpu::RunExceptioonStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,const std::string& msg,const std::string& file)
{
  throw JException(srcfile,srcline,"JDebugSphCpu",method,msg,file);
}

//==============================================================================
/// Returns dynamic pointer with code-type of particles.
/// (this pointer must be deleted)
//==============================================================================
byte* JDebugSphCpu::GetCodeType(unsigned n,const typecode* code){
  byte* codetype=JDataArrays::NewArrayByte(n,false);
  for(unsigned p=0;p<n;p++){
    const typecode type=CODE_GetType(code[p]);
    codetype[p]=(type==CODE_TYPE_FIXED? 0: (type==CODE_TYPE_MOVING? 1: (type==CODE_TYPE_FLOATING? 2: (type==CODE_TYPE_FLUID? 3: 99))));
  }
  return(codetype);
}

//==============================================================================
/// Returns dynamic pointer with code-typevalue of particles.
/// (this pointer must be deleted)
//==============================================================================
typecode* JDebugSphCpu::GetCodeTypeValue(unsigned n,const typecode* code){
  #ifdef CODE_SIZE4
    typecode* codetval=JDataArrays::NewArrayUint(n,false);
  #else
    typecode* codetval=JDataArrays::NewArrayWord(n,false);
  #endif
  for(unsigned c=0;c<n;c++)codetval[c]=CODE_GetTypeValue(code[c]);
  return(codetval);
}

//==============================================================================
/// Returns dynamic pointer with cell coordinates of particles. (this pointer 
/// must be deleted).
//==============================================================================
tuint3* JDebugSphCpu::GetCell3(unsigned n,const unsigned* dcell
  ,unsigned cellcode){

  tuint3* cell3=JDataArrays::NewArrayUint3(n,false);
  for(unsigned c=0;c<n;c++){
    const unsigned dcel=dcell[c];
    cell3[c]=TUint3(unsigned(DCEL_Cellx(cellcode,dcel)),
                    unsigned(DCEL_Celly(cellcode,dcel)),
                    unsigned(DCEL_Cellz(cellcode,dcel)));
  }
  return(cell3);
}

//==============================================================================
/// Returns dynamic pointer with position as tfloat3. (this pointer must be deleted)
//==============================================================================
tfloat3* JDebugSphCpu::GetPosf3(unsigned n,const tdouble3* pos){
  tfloat3* posf=JDataArrays::NewArrayFloat3(n,false);
  for(unsigned c=0;c<n;c++)posf[c]=ToTFloat3(pos[c]);
  return(posf);
}

//==============================================================================
/// Checks list of variables and returns the unknown variable.
//==============================================================================
std::string JDebugSphCpu::PrepareVars(const std::string& vlist){
  return(string(",")+fun::StrLower(vlist)+",");
}

//==============================================================================
/// Checks list of variables and returns the unknown variable.
//==============================================================================
std::string JDebugSphCpu::CheckVars(std::string vlist){
  string allvars= ",all,idp,seq,cell,code,vel,rho,ace,ar";
  allvars=allvars+",boundnor,motionvel,motionace,spstaurho2,sps2strain,";
  vlist=fun::StrLower(vlist);
  while(!vlist.empty()){
    string var=fun::StrSplit(",",vlist);
    if(!var.empty() && !FindVar(var,allvars))return(var);
  }
  return("");
}

//==============================================================================
/// Returns a file name.
//============================================================================== 
std::string JDebugSphCpu::GetFileName(std::string filename,int numfile){
  const string fext=fun::GetExtension(filename);
  string file=AppInfo.GetDirOut()+fun::GetWithoutExtension(filename)+"."+fext;
  if(numfile>=0)file=fun::FileNameSec(file,numfile);
  return(file);
}

//==============================================================================
/// Applies special filters to data arrays.
//============================================================================== 
void JDebugSphCpu::RunUserFilters(JDataArrays& arrays){
}

//==============================================================================
/// Loads particle data in object ffdata
/// Carga datos de las particulas en el objeto ffdata.
//==============================================================================
void JDebugSphCpu::LoadParticlesData(const JSphCpuSingle* cp,unsigned pini
  ,unsigned pfin,std::string vars,JDataArrays* arrays,std::string file)
{
  //if(DG)AppInfo.LogPtr()->Printf("file [%s]",file.c_str());
  if(pini>=pfin)Run_ExceptioonFileSta(fun::PrintStr("Invalid data (pini:%u >= pfin:%u)."
    ,pini,pfin),file);
  const unsigned n=pfin-pini;
  //-Checks list of variables.
  vars=PrepareVars(vars);
  string errvar=CheckVars(vars);
  if(!errvar.empty())Run_ExceptioonFileSta(fun::PrintStr("The variable \'%s\' is unknown.",errvar.c_str()),file);
  const bool all=FindVar("all",vars);
  //-Loads data in arrays object.
  arrays->Reset();
  if(all || FindVar("idp",vars)){
    const unsigned* idp=AC_CPTRV(cp->Idp_c,pini);
    arrays->AddArray("Idp",n,JDataArrays::NewArrayCpyUint(n,idp),true);
  }
  if(all || FindVar("seq",vars)){
    arrays->AddArray("Seq",n,JDataArrays::NewArraySeqUint(n,pini,1),true);
  }
  //-Loads dcell.
  if(all || FindVar("cell",vars)){
    const unsigned* dcell=AC_CPTRV(cp->Dcell_c,pini);
    arrays->AddArray("Cell",n,GetCell3(n,dcell,cp->DomCellCode),true);
  }
  //-Loads vel and rhop.
  if(all || FindVar("vel",vars) || FindVar("rho",vars)){
    const tfloat4* velrho=AC_CPTRV(cp->Velrho_c,pini);
    if(all || FindVar("vel",vars))arrays->AddArray("Vel",n,JDataArrays::NewArrayFloat3xyz(n,velrho),true);
    if(all || FindVar("rho",vars))arrays->AddArray("Rho",n,JDataArrays::NewArrayFloat1w  (n,velrho),true);
  }
  //-Loads code.
  if(all || FindVar("code",vars)){
    const typecode* code=AC_CPTRV(cp->Code_c,pini);
    #ifdef CODE_SIZE4
      arrays->AddArray("Code",n,JDataArrays::NewArrayCpyUint(n,code),true);
    #else
      arrays->AddArray("Code",n,JDataArrays::NewArrayCpyWord(n,code),true);
    #endif
    arrays->AddArray("Type",n,GetCodeType(n,code),true);
    byte* typesp=new byte[n];
    for(unsigned p=0;p<n;p++)typesp[p]=byte(CODE_GetSpecialByte(code[p]));
    arrays->AddArray("TypeSp",n,typesp,true);
    arrays->AddArray("TypeValue",n,GetCodeTypeValue(n,code),true);
  }
  //-Loads pos.
  {
    const tdouble3* pos=AC_CPTRV(cp->Pos_c,pini);
    arrays->AddArray("Pos",n,JDataArrays::NewArrayCpyDouble3(n,pos),true);
  }
  //-Loads boundnor.
  if(all || FindVar("boundnor",vars)){
    const tfloat3* ptr=AC_CPTRV(cp->BoundNor_c,pini);
    if(ptr)arrays->AddArray("BoundNor",n,JDataArrays::NewArrayCpyFloat3(n,ptr),true);
    else if(!all)Run_ExceptioonFileSta("The variable BoundNorc is NULL.",file);
  }
  //-Loads motionvel.
  if(all || FindVar("motionvel",vars)){
    const tfloat3* ptr=AC_CPTRV(cp->MotionVel_c,pini);
    if(ptr)arrays->AddArray("MotionVel",n,JDataArrays::NewArrayCpyFloat3(n,ptr),true);
    else if(!all)Run_ExceptioonFileSta("The variable MotionVelc is NULL.",file);
  }
  //-Loads motionace.
  if(all || FindVar("motionace",vars)){
    const tfloat3* ptr=AC_CPTRV(cp->MotionAce_c,pini);
    if(ptr)arrays->AddArray("MotionAce",n,JDataArrays::NewArrayCpyFloat3(n,ptr),true);
    else if(!all)Run_ExceptioonFileSta("The variable MotionAcec is NULL.",file);
  }
  //-Loads ace.
  if(all || FindVar("ace",vars)){
    const tfloat3* ptr=AC_CPTRV(cp->Ace_c,pini);
    if(ptr)arrays->AddArray("Ace",n,JDataArrays::NewArrayCpyFloat3(n,ptr),true);
    else if(!all)Run_ExceptioonFileSta("The variable Acec is NULL.",file);
  }
  //-Loads ar.
  if(all || FindVar("ar",vars)){
    const float* ptr=AC_CPTRV(cp->Ar_c,pini);
    if(ptr)arrays->AddArray("Ar",n,JDataArrays::NewArrayCpyFloat(n,ptr),true);
    else if(!all)Run_ExceptioonFileSta("The variable Arc is NULL.",file);
  }
  //-Loads spstau.
  if(all || FindVar("spstaurho2",vars)){
    const tsymatrix3f* ptr=AC_CPTRV(cp->SpsTauRho2_c,pini);
    if(ptr){
      const unsigned* idpc=arrays->GetArrayUint("Idp",n);
      float* ptr_xx=arrays->CreateArrayPtrFloat("SpsTauRho2_xx",n);
      float* ptr_xy=arrays->CreateArrayPtrFloat("SpsTauRho2_xy",n);
      float* ptr_xz=arrays->CreateArrayPtrFloat("SpsTauRho2_xz",n);
      float* ptr_yy=arrays->CreateArrayPtrFloat("SpsTauRho2_yy",n);
      float* ptr_yz=arrays->CreateArrayPtrFloat("SpsTauRho2_yz",n);
      float* ptr_zz=arrays->CreateArrayPtrFloat("SpsTauRho2_zz",n);
      const unsigned casenbound=cp->CaseNbound;
      for(unsigned p=0;p<n;p++){
        const tsymatrix3f t=ptr[p];
        if(idpc[p]>=casenbound){
          ptr_xx[p]=t.xx;
          ptr_xy[p]=t.xy;
          ptr_xz[p]=t.xz;
          ptr_yy[p]=t.yy;
          ptr_yz[p]=t.yz;
          ptr_zz[p]=t.zz;
        }
      }
    }
    else if(!all)Run_ExceptioonFileSta("The variable SpsTau_g is NULL.",file);
  }
  //-Loads spsgrad.
  if(all || FindVar("sps2strain",vars)){
    const tsymatrix3f* ptr=AC_CPTRV(cp->Sps2Strain_c,pini);
    if(ptr){
      const unsigned* idpc=arrays->GetArrayUint("Idp",n);
      float* ptr_xx=arrays->CreateArrayPtrFloat("Sps2Strain_xx",n);
      float* ptr_xy=arrays->CreateArrayPtrFloat("Sps2Strain_xy",n);
      float* ptr_xz=arrays->CreateArrayPtrFloat("Sps2Strain_xz",n);
      float* ptr_yy=arrays->CreateArrayPtrFloat("Sps2Strain_yy",n);
      float* ptr_yz=arrays->CreateArrayPtrFloat("Sps2Strain_yz",n);
      float* ptr_zz=arrays->CreateArrayPtrFloat("Sps2Strain_zz",n);
      const unsigned casenbound=cp->CaseNbound;
      for(unsigned p=0;p<n;p++){
        const tsymatrix3f t=ptr[p];
        if(idpc[p]>=casenbound){
          ptr_xx[p]=t.xx;
          ptr_xy[p]=t.xy;
          ptr_xz[p]=t.xz;
          ptr_yy[p]=t.yy;
          ptr_yz[p]=t.yz;
          ptr_zz[p]=t.zz;
        }
      }
    }
    else if(!all)Run_ExceptioonFileSta("The variable SpsGradvel_g is NULL.",file);
  }
}

//==============================================================================
/// Stores data in VTK format.
//============================================================================== 
void JDebugSphCpu::SaveVtk(std::string filename,int numfile,unsigned pini
  ,unsigned pfin,std::string vars,const JSphCpuSingle* cp)
{
  const string file=GetFileName(filename,numfile);
  //-Loads particle data.
  JDataArrays arrays;
  LoadParticlesData(cp,pini,pfin,vars,&arrays,file);
  //-Run user filter.
  RunUserFilters(arrays);
  //-Saves VTK file.
  JSpVtkData::Save(file,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Stores data in CSV format.
//============================================================================== 
void JDebugSphCpu::SaveCsv(std::string filename,int numfile,unsigned pini
  ,unsigned pfin,std::string vars,const JSphCpuSingle* cp)
{
  const string file=GetFileName(filename,numfile);
  //-Loads particle data.
  JDataArrays arrays;
  LoadParticlesData(cp,pini,pfin,vars,&arrays,file);
  //-Sort by Idp when it exists.
  if(arrays.ExistsName("Idp"))arrays.SortDataBy("Idp");
  //-Run user filter.
  RunUserFilters(arrays);
  //-Saves CSV file.
  JOutputCsv ocsv(AppInfo.GetCsvSepComa());
  ocsv.SaveCsv(file,arrays);
  arrays.Reset();
}

