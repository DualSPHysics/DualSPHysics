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

#include "JDebugSphGpu.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JException.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include "JSpVtkData.h"
#include "JOutputCsv.h"
#include "JDataArrays.h"

#include "JSphGpuSingle.h"

#pragma warning(disable : 4996) //Cancels sprintf() deprecated.

using namespace std;

//##############################################################################
//# JDebugSphGpu
//##############################################################################
//==============================================================================
/// Throws exception related to a file from a static method.
//==============================================================================
void JDebugSphGpu::RunExceptioonStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,const std::string& msg,const std::string& file)
{
  throw JException(srcfile,srcline,"JDebugSphGpu",method,msg,file);
}
//==============================================================================
/// Throws exception related to a CUDA error from a static method.
//==============================================================================
void JDebugSphGpu::RunExceptioonCudaStatic(const std::string& srcfile,int srcline
  ,const std::string& method
  ,cudaError_t cuerr,std::string msg)
{
  msg=msg+fun::PrintStr(" (CUDA error %d (%s)).\n",cuerr,cudaGetErrorString(cuerr));
  throw JException(srcfile,srcline,"JDebugSphGpu",method,msg,"");
}
//==============================================================================
/// Checks CUDA error and throws exception from a static method.
//==============================================================================
void JDebugSphGpu::CheckCudaErroorStatic(const std::string& srcfile,int srcline
  ,const std::string& method,std::string msg)
{
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess)RunExceptioonCudaStatic(srcfile,srcline,method,cuerr,msg);
}

//==============================================================================
/// Returns dynamic pointer with code-type of particles. (this pointer must be deleted)
//==============================================================================
byte* JDebugSphGpu::GetCodeType(unsigned n,const typecode* code){
  byte* codetype=JDataArrays::NewArrayByte(n,false);
  for(unsigned p=0;p<n;p++){
    const typecode type=CODE_GetType(code[p]);
    codetype[p]=(type==CODE_TYPE_FIXED? 0: (type==CODE_TYPE_MOVING? 1: (type==CODE_TYPE_FLOATING? 2: (type==CODE_TYPE_FLUID? 3: 99))));
  }
  return(codetype);
}

//==============================================================================
/// Returns dynamic pointer with code-typevalue of particles. (this pointer must be deleted)
//==============================================================================
typecode* JDebugSphGpu::GetCodeTypeValue(unsigned n,const typecode* code){
  #ifdef CODE_SIZE4
    typecode* codetval=JDataArrays::NewArrayUint(n,false);
  #else
    typecode* codetval=JDataArrays::NewArrayWord(n,false);
  #endif
  for(unsigned c=0;c<n;c++)codetval[c]=CODE_GetTypeValue(code[c]);
  return(codetval);
}

//==============================================================================
/// Returns dynamic pointer with cell coordinates of particles. (this pointer must be deleted)
//==============================================================================
tuint3* JDebugSphGpu::GetCell3(unsigned n,const unsigned* dcell,unsigned cellcode){
  tuint3* cell3=JDataArrays::NewArrayUint3(n,false);
  for(unsigned c=0;c<n;c++){
    const unsigned dcel=dcell[c];
    cell3[c]=TUint3(unsigned(DCEL_Cellx(cellcode,dcel)),unsigned(DCEL_Celly(cellcode,dcel)),unsigned(DCEL_Cellz(cellcode,dcel)));
  }
  return(cell3);
}

//==============================================================================
/// Returns dynamic pointer with position as tfloat3. (this pointer must be deleted)
//==============================================================================
tfloat3* JDebugSphGpu::GetPosf3(unsigned n,const tdouble3* pos){
  tfloat3* posf=JDataArrays::NewArrayFloat3(n,false);
  for(unsigned c=0;c<n;c++)posf[c]=ToTFloat3(pos[c]);
  return(posf);
}

//==============================================================================
/// Returns dynamic pointer with position as tfloat3. (this pointer must be deleted)
//==============================================================================
tfloat3* JDebugSphGpu::GetPosf3(unsigned n,const tdouble2* posxy
  ,const double* posz)
{
  tfloat3* posf=JDataArrays::NewArrayFloat3(n,false);
  for(unsigned c=0;c<n;c++)posf[c]=TFloat3(float(posxy[c].x),float(posxy[c].y),float(posz[c]));
  return(posf);
}

//==============================================================================
/// Returns dynamic pointer with position as tdouble3. (this pointer must be deleted)
//==============================================================================
tdouble3* JDebugSphGpu::GetPosd3(unsigned n,const tdouble2* posxy
  ,const double* posz)
{
  tdouble3* posd=JDataArrays::NewArrayDouble3(n,false);
  for(unsigned c=0;c<n;c++)posd[c]=TDouble3(posxy[c].x,posxy[c].y,posz[c]);
  return(posd);
}

//==============================================================================
/// Returns dynamic pointer with relative position from PosCell. (this pointer must be deleted)
//==============================================================================
tfloat3* JDebugSphGpu::GetPosCell_Pos(unsigned n,const tfloat4* poscell){
  tfloat3* posf=JDataArrays::NewArrayFloat3(n,false);
  for(unsigned c=0;c<n;c++){
    const tfloat4 ps=poscell[c];
    posf[c]=TFloat3(ps.x,ps.y,ps.z);
  }
  return(posf);
}

//==============================================================================
/// Returns dynamic pointer with cell coordinates from PosCell. (this pointer must be deleted)
//==============================================================================
tuint3* JDebugSphGpu::GetPosCell_Cell(unsigned n,const tfloat4* poscell){
  const tuint4* poscellu=(const tuint4*)poscell;
  tuint3* cell=JDataArrays::NewArrayUint3(n,false);
  for(unsigned c=0;c<n;c++){
    const unsigned cellu=poscellu[c].w;
    cell[c]=TUint3(PSCEL_GetX(cellu),PSCEL_GetY(cellu),PSCEL_GetZ(cellu));
  }
  return(cell);
}

//==============================================================================
/// Checks list of variables and returns the unknown variable.
//==============================================================================
std::string JDebugSphGpu::PrepareVars(const std::string& vlist){
  return(string(",")+fun::StrLower(vlist)+",");
}

//==============================================================================
/// Checks list of variables and returns the unknown variable.
//==============================================================================
std::string JDebugSphGpu::CheckVars(std::string vlist){
  string allvars= ",all,idp,gid,seq,cell,code,vel,rho,ace,ar,viscdt,poscell";
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
std::string JDebugSphGpu::GetFileName(std::string filename,int numfile,int gid){
  const string fext=fun::GetExtension(filename);
  string file=AppInfo.GetDirOut()+fun::GetWithoutExtension(filename)
    +(gid>=0? fun::PrintStr("_g%d",gid): "")+"."+fext;
  if(numfile>=0)file=fun::FileNameSec(file,numfile);
  return(file);
}

//==============================================================================
/// Applies special filters to data arrays.
//============================================================================== 
void JDebugSphGpu::RunUserFilters(JDataArrays& arrays){
}

//==============================================================================
/// Loads particle data in object ffdata
/// Carga datos de las particulas en el objeto ffdata.
//==============================================================================
void JDebugSphGpu::LoadParticlesData(const JSphGpuSingle* gp,unsigned pini
  ,unsigned pfin,std::string vars,JDataArrays* arrays,std::string file)
{
  //const bool DG=(file.find("_FcHaloRig")>=0);
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
    const unsigned* idp=fcuda::ToHostUint(pini,n,gp->Idp_g->cptr());
    arrays->AddArray("Idp",n,idp,true);
  }
  if(all || FindVar("seq",vars)){
    arrays->AddArray("Seq",n,JDataArrays::NewArraySeqUint(n,pini,1),true);
  }
  //-Loads dcell.
  if(all || FindVar("cell",vars)){
    const unsigned* dcell=fcuda::ToHostUint(pini,n,gp->Dcell_g->cptr());
    arrays->AddArray("Cell",n,GetCell3(n,dcell,gp->DomCellCode),true);
    delete[] dcell; dcell=NULL;
  }
  //-Loads vel and rhop.
  if(all || FindVar("vel",vars) || FindVar("rho",vars)){
    const tfloat4* velrho=fcuda::ToHostFloat4(pini,n,gp->Velrho_g->cptr());
    if(all || FindVar("vel",vars)) arrays->AddArray("Vel",n,JDataArrays::NewArrayFloat3xyz(n,velrho),true);
    if(all || FindVar("rho",vars))arrays->AddArray("Rho",n,JDataArrays::NewArrayFloat1w(n,velrho),true);
    delete[] velrho; velrho=NULL;
  }
  //-Loads code.
  if(all || FindVar("code",vars)){
    #ifdef CODE_SIZE4
      const typecode* code=fcuda::ToHostUint(pini,n,gp->Code_g->cptr());
    #else
      const typecode* code=fcuda::ToHostWord(pini,n,gp->Code_g->cptr());
    #endif
    arrays->AddArray("Code",n,code,true);
    arrays->AddArray("Type",n,GetCodeType(n,code),true);
    byte* typesp=new byte[n];
    for(unsigned p=0;p<n;p++)typesp[p]=byte(CODE_GetSpecialByte(code[p]));
    arrays->AddArray("TypeSp",n,typesp,true);
    arrays->AddArray("TypeValue",n,GetCodeTypeValue(n,code),true);
  }
  //-Loads pos.
  {
    const tdouble2* posxy=fcuda::ToHostDouble2(pini,n,gp->Posxy_g->cptr());
    const double*   posz =fcuda::ToHostDouble (pini,n,gp->Posz_g->cptr());
    if(FindVar("posd",vars))arrays->AddArray("Pos",n,GetPosd3(n,posxy,posz),true);
    else                    arrays->AddArray("Pos",n,GetPosf3(n,posxy,posz),true);
    delete[] posxy; posxy=NULL;
    delete[] posz;  posz=NULL;
  }
  //-Loads poscell. //<vs_tesmode_ini> 
  if(all || FindVar("poscell",vars)){
    const float4* ptrg=gp->PosCell_g->cptr();
    if(ptrg){
      const tfloat4* poscell=fcuda::ToHostFloat4(pini,n,ptrg);
      arrays->AddArray("poscell_Pos" ,n,GetPosCell_Pos (n,poscell),true);
      arrays->AddArray("poscell_Cell",n,GetPosCell_Cell(n,poscell),true);
      delete[] poscell; poscell=NULL;
    }
    else if(!all)Run_ExceptioonFileSta("The variable PosCellg is NULL.",file);
  } //<vs_tesmode_end> 
  //-Loads boundnor.
  if(all || FindVar("boundnor",vars)){
    const float3* ptrg=AG_CPTR(gp->BoundNor_g);
    if(ptrg)arrays->AddArray("BoundNor",n,fcuda::ToHostFloat3(pini,n,ptrg),true);
    else if(!all)Run_ExceptioonFileSta("The variable BoundNorg is NULL.",file);
  }
  //-Loads motionvel.
  if(all || FindVar("motionvel",vars)){
    const float3* ptrg=AG_CPTR(gp->MotionVel_g);
    if(ptrg)arrays->AddArray("MotionVel",n,fcuda::ToHostFloat3(pini,n,ptrg),true);
    else if(!all)Run_ExceptioonFileSta("The variable MotionVelg is NULL.",file);
  }
  //-Loads motionace.
  if(all || FindVar("motionace",vars)){
    const float3* ptrg=AG_CPTR(gp->MotionAce_g);
    if(ptrg)arrays->AddArray("MotionAce",n,fcuda::ToHostFloat3(pini,n,ptrg),true);
    else if(!all)Run_ExceptioonFileSta("The variable MotionAceg is NULL.",file);
  }
  //-Loads ace.
  if(all || FindVar("ace",vars)){
    const float3* ptrg=gp->Ace_g->cptr();
    if(ptrg)arrays->AddArray("Ace",n,fcuda::ToHostFloat3(pini,n,ptrg),true);
    else if(!all)Run_ExceptioonFileSta("The variable Aceg is NULL.",file);
  }
  //-Loads ar.
  if(all || FindVar("ar",vars)){
    const float* ptrg=gp->Ar_g->cptr();
    if(ptrg)arrays->AddArray("Ar",n,fcuda::ToHostFloat(pini,n,ptrg),true);
    else if(!all)Run_ExceptioonFileSta("The variable Arg is NULL.",file);
  }
  //-Loads viscdt.
  if(all || FindVar("viscdt",vars)){
    const float* ptrg=AG_CPTR(gp->ViscDt_g);
    if(ptrg)arrays->AddArray("ViscDt",n,fcuda::ToHostFloat(pini,n,ptrg),true);
    else if(!all)Run_ExceptioonFileSta("The variable ViscDtg is NULL.",file);
  }
  //-Loads spstau.
  if(all || FindVar("spstaurho2",vars)){
    const float* ptrg=(float*)AG_CPTR(gp->SpsTauRho2_g);
    if(ptrg){
      const unsigned* idpc=arrays->GetArrayUint("Idp",n);
      const tsymatrix3f* spstauc=(const tsymatrix3f*)fcuda::ToHostFloat(pini,n*6,(float*)ptrg);
      float* ptr_xx=arrays->CreateArrayPtrFloat("SpsTauRho2_xx",n);
      float* ptr_xy=arrays->CreateArrayPtrFloat("SpsTauRho2_xy",n);
      float* ptr_xz=arrays->CreateArrayPtrFloat("SpsTauRho2_xz",n);
      float* ptr_yy=arrays->CreateArrayPtrFloat("SpsTauRho2_yy",n);
      float* ptr_yz=arrays->CreateArrayPtrFloat("SpsTauRho2_yz",n);
      float* ptr_zz=arrays->CreateArrayPtrFloat("SpsTauRho2_zz",n);
      const unsigned casenbound=gp->CaseNbound;
      for(unsigned p=0;p<n;p++){
        const tsymatrix3f t=spstauc[p];
        if(idpc[p]>=casenbound){
          ptr_xx[p]=t.xx;
          ptr_xy[p]=t.xy;
          ptr_xz[p]=t.xz;
          ptr_yy[p]=t.yy;
          ptr_yz[p]=t.yz;
          ptr_zz[p]=t.zz;
        }
      }
      delete[] (float*)spstauc;
    }
    else if(!all)Run_ExceptioonFileSta("The variable SpsTau_g is NULL.",file);
  }
  //-Loads spsgrad.
  if(all || FindVar("sps2strain",vars)){
    const float* ptrg=(float*)AG_CPTR(gp->Sps2Strain_g);
    if(ptrg){
      const unsigned* idpc=arrays->GetArrayUint("Idp",n);
      const tsymatrix3f* spsgradc=(const tsymatrix3f*)fcuda::ToHostFloat(pini,n*6,(float*)ptrg);
      float* ptr_xx=arrays->CreateArrayPtrFloat("Sps2Strain_xx",n);
      float* ptr_xy=arrays->CreateArrayPtrFloat("Sps2Strain_xy",n);
      float* ptr_xz=arrays->CreateArrayPtrFloat("Sps2Strain_xz",n);
      float* ptr_yy=arrays->CreateArrayPtrFloat("Sps2Strain_yy",n);
      float* ptr_yz=arrays->CreateArrayPtrFloat("Sps2Strain_yz",n);
      float* ptr_zz=arrays->CreateArrayPtrFloat("Sps2Strain_zz",n);
      const unsigned casenbound=gp->CaseNbound;
      for(unsigned p=0;p<n;p++){
        const tsymatrix3f t=spsgradc[p];
        if(idpc[p]>=casenbound){
          ptr_xx[p]=t.xx;
          ptr_xy[p]=t.xy;
          ptr_xz[p]=t.xz;
          ptr_yy[p]=t.yy;
          ptr_yz[p]=t.yz;
          ptr_zz[p]=t.zz;
        }
      }
      delete[] (float*)spsgradc;
    }
    else if(!all)Run_ExceptioonFileSta("The variable SpsGradvel_g is NULL.",file);
  }
}

//==============================================================================
/// Stores data in VTK format.
//============================================================================== 
void JDebugSphGpu::SaveVtk(std::string filename,int numfile,unsigned pini
  ,unsigned pfin,std::string vars,const JSphGpuSingle* gp)
{
  const string file=GetFileName(filename,numfile);
  //-Loads particle data.
  JDataArrays arrays;
  LoadParticlesData(gp,pini,pfin,vars,&arrays,file);
  //-Run user filter.
  RunUserFilters(arrays);
  //-Saves VTK file.
  JSpVtkData::Save(file,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Stores data in CSV format.
//============================================================================== 
void JDebugSphGpu::SaveCsv(std::string filename,int numfile,unsigned pini
  ,unsigned pfin,std::string vars,const JSphGpuSingle* gp)
{
  const string file=GetFileName(filename,numfile);
  //-Loads particle data.
  JDataArrays arrays;
  LoadParticlesData(gp,pini,pfin,vars,&arrays,file);
  //-Sort by Idp when it exists.
  if(arrays.ExistsName("Idp"))arrays.SortDataBy("Idp");
  //-Run user filter.
  RunUserFilters(arrays);
  //-Saves CSV file.
  JOutputCsv ocsv(AppInfo.GetCsvSepComa());
  ocsv.SaveCsv(file,arrays);
  arrays.Reset();
}


