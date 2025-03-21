//HEAD_DSPH
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

/// \file JCellDivGpu.cpp \brief Implements the class \ref JCellDivGpu.

#include "JCellDivGpu.h"
#include "JCellDivGpu_ker.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunctionsCuda.h"

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCellDivGpu::JCellDivGpu(bool stable,bool floating,byte periactive
  ,float kernelsize2,float poscellsize
  ,bool celldomfixed,TpCellMode cellmode,float scell
  ,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
  ,unsigned casenbound,unsigned casenfixed,unsigned casenpb
  ,std::string dirout)
  :Log(AppInfo.LogPtr()),Stable(stable),Floating(floating),PeriActive(periactive)
  ,CellDomFixed(celldomfixed),CellMode(cellmode)
  ,ScellDiv(cellmode==CELLMODE_Full? 1: (cellmode==CELLMODE_Half? 2: 0))
  ,Scell(scell),OvScell(1.f/scell),KernelSize2(kernelsize2),PosCellSize(poscellsize)
  ,Map_PosMin(mapposmin),Map_PosMax(mapposmax),Map_PosDif(mapposmax-mapposmin)
  ,Map_Cells(mapcells),CaseNbound(casenbound),CaseNfixed(casenfixed)
  ,CaseNpb(casenpb),DirOut(dirout)
  ,OverMemoryNp(CELLDIV_OVERMEMORYNP)
  ,OverMemoryCells(CELLDIV_OVERMEMORYCELLS)
  ,OverMemoryNCells(CELLDIV_OVERMEMORYNCELLS)
{
  ClassName="JCellDivGpu";
  CellPart=NULL;  SortPart=NULL;  AuxMem=NULL;
  BeginEndCell=NULL;
  SortPart2=NULL; SortIdx=NULL; //<vs_flexstruc>
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JCellDivGpu::~JCellDivGpu(){
  DestructorActive=true;
  //Log->Printf("---> DivideFull:%u/%u",NdivFull,Ndiv); //:del:
  Reset();
}
 
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JCellDivGpu::Reset(){
  SizeNp=SizeAuxMem=SizeNct=0;
  IncreaseNp=0;
  FreeMemoryAll();
  Ndiv=NdivFull=0;
  Nptot=Npb1=Npf1=Npb2=Npf2=0;
  MemAllocGpuNp=MemAllocGpuNct=0;
  MemAllocGpuNpTimes=MemAllocGpuNctTimes=0;
  NpbOut=NpfOut=NpbOutIgnore=NpfOutIgnore=0;
  NpFinal=NpbFinal=0;
  NpbIgnore=0;
  CellDomainMin=TUint3(1);
  CellDomainMax=TUint3(0);
  Ncx=Ncy=Ncz=Nsheet=Nct=0;
  Nctt=0;
  BoundLimitOk=BoundDivideOk=false;
  BoundLimitCellMin=BoundLimitCellMax=TUint3(0);
  BoundDivideCellMin=BoundDivideCellMax=TUint3(0);
  DivideFull=false;
}

//==============================================================================
/// Free memory reserved for cells.
/// Libera memoria reservada para celdas.
//==============================================================================
void JCellDivGpu::FreeMemoryNct(){
  cudaFree(BeginEndCell); BeginEndCell=NULL;
  MemAllocGpuNct=0;
  BoundDivideOk=false;
}

//==============================================================================
/// Free the basic memory allocated for particles and cells.
/// Libera memoria basica reservada para particulas y celdas.
//==============================================================================
void JCellDivGpu::FreeMemoryAll(){
  FreeMemoryNct();
  cudaFree(CellPart);   CellPart=NULL;
  cudaFree(SortPart);   SortPart=NULL;
  cudaFree(AuxMem);     AuxMem=NULL;
  cudaFree(SortPart2);  SortPart2=NULL; //<vs_flexstruc>
  cudaFree(SortIdx);    SortIdx=NULL;   //<vs_flexstruc>
  MemAllocGpuNp=0;
  BoundDivideOk=false;
}

//==============================================================================
/// Assigns basic memory according to the particle number.
/// Asigna memoria basica segun numero de particulas.
//==============================================================================
void JCellDivGpu::AllocMemoryNp(ullong np,ullong npmin){
  FreeMemoryAll();
  np=np+PARTICLES_OVERMEMORY_MIN;
  SizeNp=unsigned(np);
  //-Checks particle number.
  //-Comprueba numero de particulas.
  if(np!=SizeNp)Run_Exceptioon(string("Failed GPU memory allocation for ")+fun::UlongStr(np)+" particles.");
  //-Allocates memory for the particles.
  //-Reserva memoria para particulas.
  MemAllocGpuNp=0;
  MemAllocGpuNp+=fcuda::Malloc(&CellPart,SizeNp);
  MemAllocGpuNp+=fcuda::Malloc(&SortPart,SizeNp);
  SizeAuxMem=cudiv::LimitsPosSize(SizeNp);
  MemAllocGpuNp+=fcuda::Malloc(&AuxMem,SizeAuxMem);
  MemAllocGpuNp+=fcuda::Malloc(&SortPart2,SizeNp);  //<vs_flexstruc>
  MemAllocGpuNp+=fcuda::Malloc(&SortIdx,SizeNp);    //<vs_flexstruc>
  MemAllocGpuNpTimes++;
  //-Checks allocated memory.
  //-Comprueba reserva de memoria.
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess){
    Run_ExceptioonCuda(cuerr,fun::PrintStr(
      "Failed GPU memory allocation of %.1f MiB for %u particles."
      ,double(MemAllocGpuNp)/MEBIBYTE,SizeNp));
  }
  //-Show requested memory.
  const string txover=(npmin>1? fun::PrintStr(" (over-allocation: %.2fX)",double(SizeNp)/npmin): "");
  Log->Printf("**CellDiv: Requested gpu memory for %s particles%s: %.1f MiB (%u times)."
    ,KINT(SizeNp),txover.c_str(),double(MemAllocGpuNp)/MEBIBYTE,MemAllocGpuNpTimes);
}

//==============================================================================
/// Assigns memory according to the cell number.
/// Asigna memoria segun numero de celdas.
//==============================================================================
void JCellDivGpu::AllocMemoryNct(ullong nct,ullong nctmin){
  FreeMemoryNct();
  SizeNct=unsigned(nct);
  const unsigned nctt=unsigned(SizeBeginEndCell(SizeNct));
  //-Checks number of cells.
  if(ullong(nctt)!=SizeBeginEndCell(nct))Run_Exceptioon(
    fun::PrintStr("Failed GPU memory allocation for %s cells.",KINT(nct)));
  //-Allocates memory for cells.
  MemAllocGpuNct=fcuda::Malloc(&BeginEndCell,nctt);
  MemAllocGpuNctTimes++;
  //-Checks allocated memory.
  const cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess){
    Run_ExceptioonCuda(cuerr,fun::PrintStr(
      "Failed GPU memory allocation of %.1f MiB for %s cells."
      ,double(MemAllocGpuNct)/MEBIBYTE,KINT(SizeNct)));
  }
  //-Shows requested memory.
  const string txover=(nctmin>1? fun::PrintStr(" (over-allocation: %.2fX)",double(SizeNct)/nctmin): "");
  Log->Printf("**CellDiv: Requested GPU memory for %s cells%s: %.1f MiB (%u times)."
    ,KINT(SizeNct),txover.c_str(),double(MemAllocGpuNct)/MEBIBYTE,MemAllocGpuNctTimes);
}

//==============================================================================
/// Checks allocated memory for the indicated number of particles.
/// If the allocated memory is not sufficient, reserve the required memory.
///
/// Comprueba la reserva de memoria para el numero indicado de particulas. 
/// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivGpu::CheckMemoryNp(unsigned npmin){
  if(SizeNp<npmin+PARTICLES_OVERMEMORY_MIN){
    AllocMemoryNp(ullong(npmin)+ullong(OverMemoryNp*npmin)+IncreaseNp,npmin);
    IncreaseNp=0;
  }
  else if(!CellPart){
    AllocMemoryNp(SizeNp+IncreaseNp,npmin);
    IncreaseNp=0;
  }
}

//==============================================================================
/// Checks allocated memory for the indicated number of cells.
/// If the allocated memory is not sufficient, reserve the required memory.
///
/// Comprueba la reserva de memoria para el numero indicado de celdas. 
/// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivGpu::CheckMemoryNct(unsigned nctmin){
  if(SizeNct<nctmin){
    unsigned nctnew=nctmin;
    if(OverMemoryCells>0){
      const ullong nct1=ullong(Ncx+OverMemoryCells)*ullong(Ncy+OverMemoryCells)*ullong(Ncz+OverMemoryCells);
      const ullong nct2=ullong(nctmin)+OverMemoryNCells;
      const ullong nct3=(nct1>nct2? nct1: nct2);
      const ullong nct4=min(DomCellsMax,nct3);
      if(SizeBeginEndCell(nct4)>=UINT_MAX)Run_Exceptioon("The number of cells is too big.");
      nctnew=unsigned(nct4);
    }
    AllocMemoryNct(nctnew,nctmin);
  }
  else if(!BeginEndCell)AllocMemoryNct(SizeNct,nctmin);  
}

//==============================================================================
/// Defines the domain to be used by the simulation.
/// Define el dominio de simulacion a usar.
//==============================================================================
void JCellDivGpu::DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin
  ,tdouble3 domposmin,tdouble3 domposmax)
{
  DomCellCode=cellcode;
  DomCelIni=domcelini;
  DomCelFin=domcelfin;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
  DomCells=DomCelFin-DomCelIni;
  DomCellsMax=ullong(DomCells.x)*ullong(DomCells.y)*ullong(DomCells.z);
  //Log->Printf("-----> MaxDomCells:%s",KINT(MaxDomCells));
}

/*:/==============================================================================
// Devuelve coordenadas de celda a partir de una posicion.
//==============================================================================
//tuint3 JCellDivGpu::GetMapCell(const tfloat3& pos)const{ pdtecell
//  float dx=pos.x-MapPosMin.x,dy=pos.y-MapPosMin.y,dz=pos.z-MapPosMin.z;
//  unsigned cx=unsigned(dx*OvScell),cy=unsigned(dy*OvScell),cz=unsigned(dz*OvScell);
//  return(TUint3(cx,cy,cz));
//}:*/

//==============================================================================
/// Computes the maximum and minimum position for the indicated range of boundary particles.
/// The excluded particles are already marked in code[].
///
/// Calcula posiciones minimas y maximas del rango de particulas Bound indicado.
/// En code[] ya estan marcadas las particulas excluidas.
//==============================================================================
void JCellDivGpu::CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2
  ,unsigned pini2,const unsigned* dcellg,const typecode* codeg
  ,tuint3& cellmin,tuint3& cellmax)
{
  tuint3 cmin,cmax;
  cudiv::LimitsCell(n,pini,DomCellCode,dcellg,codeg,(unsigned*)AuxMem,cmin,cmax);
  cellmin=(cmin.x>cmax.x? DomCells: cmin);
  cellmax=(cmin.x>cmax.x? TUint3(0): cmax);
  if(n2){
    cudiv::LimitsCell(n2,pini2,DomCellCode,dcellg,codeg,(unsigned*)AuxMem,cmin,cmax);
    cmin=(cmin.x>cmax.x? DomCells: cmin);
    cmax=(cmin.x>cmax.x? TUint3(0): cmax);
    cellmin=MinValues(cellmin,cmin);
    cellmax=MaxValues(cellmax,cmax);
  }
  //:Log->Printf("CalcDomainBound> cell:(%s)-(%s)",fun::Uint3Str(cellmin).c_str(),fun::Uint3Str(cellmax).c_str());
}

//==============================================================================
/// Computes the maximum and minimum position for the indicated range of fluid particles.
/// Ignores the excluded particles already marked in code[]
///
/// Calcula posiciones minimas y maximas del rango de particulas Fluid indicado.
/// Ignora particulas excluidas que ya estan marcadas en code[].
//==============================================================================
void JCellDivGpu::CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2
  ,unsigned pini2,const unsigned* dcellg,const typecode* codeg
  ,tuint3& cellmin,tuint3& cellmax)
{
  tuint3 cmin,cmax;
  cudiv::LimitsCell(n,pini,DomCellCode,dcellg,codeg,(unsigned*)AuxMem,cmin,cmax);
  cellmin=(cmin.x>cmax.x? DomCells: cmin);
  cellmax=(cmin.x>cmax.x? TUint3(0): cmax);
  if(n2){
    cudiv::LimitsCell(n2,pini2,DomCellCode,dcellg,codeg,(unsigned*)AuxMem,cmin,cmax);
    cmin=(cmin.x>cmax.x? DomCells: cmin);
    cmax=(cmin.x>cmax.x? TUint3(0): cmax);
    cellmin=MinValues(cellmin,cmin);
    cellmax=MaxValues(cellmax,cmax);
  }
  //:Log->Printf("CalcDomainFluid> cell:(%s)-(%s)",fun::Uint3Str(cellmin).c_str(),fun::Uint3Str(cellmax).c_str());
}

//==============================================================================
/// Returns first and last particle of the indicated cell.
/// Devuelve principo y final de la celda indicada.
//==============================================================================
void JCellDivGpu::CellBeginEnd(unsigned cell,unsigned ndata,unsigned* data)const{
  cudaMemcpy(data,BeginEndCell+cell,sizeof(int)*ndata,cudaMemcpyDeviceToHost);
}

//==============================================================================
/// Returns first and last particle of the indicated cell.
/// Devuelve principo y final de la celda indicada.
//==============================================================================
int2 JCellDivGpu::CellBeginEnd(unsigned cell)const{
  int2 v;
  cudaMemcpy(&v,BeginEndCell+cell,sizeof(int2),cudaMemcpyDeviceToHost);
  return(v);
}

/*://==============================================================================
// Devuelve un nombre de fichero formado por los datos indicados.
//==============================================================================
std::string JCellDivGpu::GetFileName(std::string name,std::string ext,int num)const{
  int r=Log->GetMpiRank();
  if(r>=0)name=string("p")+fun::IntStr(r)+"_"+name;
  if(num>=0){
    char cad[64];
    sprintf(cad,"%04d",num);
    name=name+cad;
  }
  return(DirOut+name+ext);
}:*/

//==============================================================================
/// Reorders basic arrays according to SortPart.
/// Ordena arrays basicos segun SortPart. 
//==============================================================================
void JCellDivGpu::SortBasicArrays(const unsigned* idp,const typecode* code
  ,const unsigned* dcell,const double2* posxy,const double* posz
  ,const float4* velrhop,unsigned* idp2,typecode* code2,unsigned* dcell2
  ,double2* posxy2,double* posz2,float4* velrhop2)
{
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,idp,code,dcell,posxy,posz,velrhop
    ,idp2,code2,dcell2,posxy2,posz2,velrhop2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for type float4).
/// Ordena arrays de datos segun SortPart (para tipo float4).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float4* a,float4* a2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for float values).
/// Ordena arrays de datos segun SortPart (para valores float).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float* a,const float* b,float* a2,float* b2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,b,a2,b2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for double and double2 values).
/// Ordena arrays de datos segun SortPart (para valores double y double2).
//==============================================================================
void JCellDivGpu::SortDataArrays(const double2* a,const double* b,const float4* c
  ,double2* a2,double* b2,float4* c2)
{
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,b,c,a2,b2,c2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for 3x float3).
/// Ordena arrays de datos segun SortPart (para 3x float3).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float3* a,const float3* b,const float3* c
  ,float3* a2,float3* b2,float3* c2)
{
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,b,c,a2,b2,c2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for type tsymatrix3f).
/// Ordena arrays de datos segun SortPart (para tipo tsymatrix3f).
//==============================================================================
void JCellDivGpu::SortDataArrays(const tsymatrix3f* a,tsymatrix3f* a2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for type float3).
/// Ordena arrays de datos segun SortPart (para tipo float3).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float3* a,float3* a2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for float values).
/// Ordena arrays de datos segun SortPart (para valores float).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float* a, float* a2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for uint,float4 values).
/// Ordena arrays de datos segun SortPart (para valores uint,float4).
//==============================================================================
void JCellDivGpu::SortDataArrays(const unsigned* a,const float4* b
  ,unsigned* a2,float4* b2)
{
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,b,a2,b2);
}

//==============================================================================
/// Reorders PeriParent references.
//==============================================================================
void JCellDivGpu::SortArrayPeriParent(unsigned* aux,const unsigned* a
  ,unsigned* a2)
{
  //const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortArrayPeriParent(Nptot,SortPart,aux,a,a2);
}

//==============================================================================
/// Returns a pointer with the auxiliary memory allocated in the GPU, only
/// used as intermediate in some tasks, in order to use it in other tasks.
/// This memoery is resized according to the particle number thus its
/// size and direction can vary.
///
/// Devuelve un puntero con la memoria auxiliar reservada en GPU, que solo se usa
/// como almacenamiento intermedio durante ciertos procesos. Asi es posible
/// aprovechar esta memoria para otros usos.
/// Esta memoria se redimensiona segun el numero de particulas por lo que su
/// tamanho y direccion pueden variar.
//==============================================================================
float* JCellDivGpu::GetAuxMem(unsigned size){
  //:printf("GetAuxMem> size:%u  SizeAuxMem:%u\n",size,SizeAuxMem);
  if(size>SizeAuxMem)Run_Exceptioon("The requested memory is not available.");
  return(AuxMem);
}

//==============================================================================
/// Returns the actual limits of the domain.
/// Devuelve limites actuales del dominio.
//==============================================================================
tdouble6 JCellDivGpu::GetDomainLimitsMinMax(unsigned slicecellmin)const{
  tuint3 celmin=GetCellDomainMin(),celmax=GetCellDomainMax();
  if(celmin.x>celmax.x)celmin.x=celmax.x=0; else celmax.x++;
  if(celmin.y>celmax.y)celmin.y=celmax.y=0; else celmax.y++;
  if(celmin.z>celmax.z)celmin.z=celmax.z=slicecellmin; else celmax.z++;
  double scell=double(Scell);
  tdouble3 pmin=DomPosMin+TDouble3(scell*celmin.x,scell*celmin.y,scell*celmin.z);
  tdouble3 pmax=DomPosMin+TDouble3(scell*celmax.x,scell*celmax.y,scell*celmax.z);
  return(TDouble6(pmin,pmax));
}

//==============================================================================
/// Returns the actual limits of the domain.
/// Devuelve limites actuales del dominio.
//==============================================================================
tdouble3 JCellDivGpu::GetDomainLimits(bool limitmin,unsigned slicecellmin)const{
  const tdouble6 limits=GetDomainLimitsMinMax(slicecellmin);
  return(limitmin? limits.getlo(): limits.gethi());
}

//<vs_flexstruc_ini>
//==============================================================================
/// Sorts and updates the indices of the flexible structure particles.
/// Ordena y actualiza los índices de las partículas de estructura flexible.
//==============================================================================
void JCellDivGpu::UpdateIndices(unsigned n,unsigned* idx){
  if(DivideFull){
    cudaMemcpy(SortPart2,SortPart,sizeof(unsigned)*NpbFinal,cudaMemcpyDeviceToDevice);
    cudiv::SortIndices(SortPart2,SortIdx,NpbFinal,Stable);
    cudiv::UpdateIndices(n,SortIdx,idx);
  }
}
//<vs_flexstruc_end>

/*:
////==============================================================================
//// Devuelve rango de particulas en el rango de celdas indicadas.
////==============================================================================
//uint2 JCellDivGpu::GetRangeParticlesCells(bool fluid,unsigned celini,unsigned celfin)const{
//  if(fluid){ celini+=BoxFluid; celfin+=BoxFluid; }
//  unsigned pmin=UINT_MAX,pmax=0;
//  if(celini<celfin){
//    bool memorynew=false;
//    unsigned* auxg=NULL;
//    unsigned size=cudiv::GetRangeParticlesCellsSizeAux(celini,celfin);
//    if(size<=SizeAuxMem)auxg=(unsigned*)AuxMem;
//    else{
//      memorynew=true;
//      cudaMalloc((void**)&auxg,sizeof(unsigned)*size);
//    } 
//    cudiv::GetRangeParticlesCells(celini,celfin,BeginEndCell,auxg,pmin,pmax,Log);
//    if(memorynew)cudaFree(auxg);
//  }
//  uint2 rg; rg.x=pmin; rg.y=pmax;
//  return(rg);
//}

////==============================================================================
//// Devuelve numero de particulas en el rango de celdas indicadas.
////==============================================================================
//unsigned JCellDivGpu::GetParticlesCells(unsigned celini,unsigned celfin){
//  unsigned count=0;
//  if(celini<celfin){
//    bool memorynew=false;
//    unsigned* auxg=NULL;
//    unsigned size=cudiv::GetParticlesCellsSizeAux(celini,celfin);
//    if(size<=SizeAuxMem)auxg=(unsigned*)AuxMem;
//    else{
//      memorynew=true;
//      cudaMalloc((void**)&auxg,sizeof(unsigned)*size);
//    } 
//    count=cudiv::GetParticlesCells(celini,celfin,BeginEndCell,auxg,Log);
//    if(memorynew)cudaFree(auxg);
//  }
//  return(count);
//}


//==============================================================================
// Graba fichero vtk con el rango de particulas indicado.
//==============================================================================
void JCellDivGpu::DgSaveVktRange(std::string file,unsigned pini,unsigned pfin
  ,const unsigned* idpg,const float3* posg)const
{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)file=string("p")+fun::IntStr(mpirank)+"_"+file;
  file=DirOut+file;
  unsigned np=pfin-pini;
  tfloat3* pos=new tfloat3[np];
  unsigned* idp=new unsigned[np];
  cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*np,cudaMemcpyDeviceToHost);
  cudaMemcpy(pos,posg+pini,sizeof(float3)*np,cudaMemcpyDeviceToHost);
  JFormatFiles2::ParticlesToVtk(file,pfin-pini,pos,NULL,NULL,NULL,NULL,idp,NULL,NULL,NULL,NULL);
  delete[] pos;
  delete[] idp;
}:*/


