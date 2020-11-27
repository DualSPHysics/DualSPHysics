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

/// \file JCellDivGpu.cpp \brief Implements the class \ref JCellDivGpu.

#include "JCellDivGpu.h"
#include "JCellDivGpu_ker.h"
#include "JAppInfo.h"
#include "Functions.h"

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCellDivGpu::JCellDivGpu(bool stable,bool floating,byte periactive,float kernelsize2,float poscellsize
  ,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
  ,unsigned casenbound,unsigned casenfixed,unsigned casenpb,std::string dirout
  ,bool allocfullnct,float overmemorynp,word overmemorycells)
  :Log(AppInfo.LogPtr()),Stable(stable),Floating(floating),PeriActive(periactive)
  ,CellMode(cellmode),ScellDiv(cellmode==CELLMODE_Full? 1: (cellmode==CELLMODE_Half? 2: 0))
  ,Scell(scell),OvScell(1.f/scell),KernelSize2(kernelsize2),PosCellSize(poscellsize)
  ,Map_PosMin(mapposmin),Map_PosMax(mapposmax),Map_PosDif(mapposmax-mapposmin)
  ,Map_Cells(mapcells),CaseNbound(casenbound),CaseNfixed(casenfixed),CaseNpb(casenpb)
  ,DirOut(dirout),AllocFullNct(allocfullnct),OverMemoryNp(overmemorynp),OverMemoryCells(overmemorycells)
{
  ClassName="JCellDivGpu";
  CellPart=NULL;  SortPart=NULL;  AuxMem=NULL;
  BeginEndCell=NULL;
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
  cudaFree(CellPart);  CellPart=NULL;
  cudaFree(SortPart);  SortPart=NULL;
  cudaFree(AuxMem);    AuxMem=NULL; 
  MemAllocGpuNp=0;
  BoundDivideOk=false;
}

//==============================================================================
/// Assigns basic memory according to the particle number.
/// Asigna memoria basica segun numero de particulas.
//==============================================================================
void JCellDivGpu::AllocMemoryNp(ullong np){
  FreeMemoryAll();
  np=np+PARTICLES_OVERMEMORY_MIN;
  SizeNp=unsigned(np);
  //-Checks particle number.
  //-Comprueba numero de particulas.
  if(np!=SizeNp)Run_Exceptioon(string("Failed GPU memory allocation for ")+fun::UlongStr(np)+" particles.");
  //-Allocates memory for the particles.
  //-Reserva memoria para particulas.
  MemAllocGpuNp=0;
  size_t m=sizeof(unsigned)*SizeNp;
  cudaMalloc((void**)&CellPart,m); MemAllocGpuNp+=m;
  cudaMalloc((void**)&SortPart,m); MemAllocGpuNp+=m;
  SizeAuxMem=cudiv::LimitsPosSize(SizeNp);
  m=sizeof(float)*SizeAuxMem;
  cudaMalloc((void**)&AuxMem,m);   MemAllocGpuNp+=m;
  //-Checks allocated memory.
  //-Comprueba reserva de memoria.
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess){
    Run_ExceptioonCuda(cuerr,fun::PrintStr("Failed CPU memory allocation of %.1f MB for %u particles.",double(MemAllocGpuNp)/(1024*1024),SizeNp));
  }
  //-Displays the requested memory.
  //-Muestra la memoria solicitada.
  Log->Printf("**CellDiv: Requested gpu memory for %u particles: %.1f MB.",SizeNp,double(MemAllocGpuNp)/(1024*1024));
}

//==============================================================================
/// Assigns memory according to the cell number.
/// Asigna memoria segun numero de celdas.
//==============================================================================
void JCellDivGpu::AllocMemoryNct(ullong nct){
  FreeMemoryNct();
  SizeNct=unsigned(nct);
  //-Checks cell number.
  //-Comprueba numero de celdas.
  if(nct!=SizeNct)Run_Exceptioon(string("Failed GPU memory allocation for ")+fun::UlongStr(nct)+" cells.");
  //-Allocates memory for cells.
  //-Reserva memoria para celdas.
  MemAllocGpuNct=0;
  size_t m=sizeof(int2)*SizeBeginEndCell(SizeNct);
  cudaMalloc((void**)&BeginEndCell,m); MemAllocGpuNct+=m;
  //-Checks allocated memory.
  //-Comprueba reserva de memoria.
  cudaError_t cuerr=cudaGetLastError();
  if(cuerr!=cudaSuccess){
    Run_ExceptioonCuda(cuerr,fun::PrintStr("Failed GPU memory allocation of %.1f MB for %u cells.",double(MemAllocGpuNct)/(1024*1024),SizeNct));
  }
  //-Displays requested memory.
  //-Muestra la memoria solicitada.
  Log->Printf("**CellDiv: Requested gpu memory for %u cells (CellMode=%s): %.1f MB.",SizeNct,GetNameCellMode(CellMode),double(MemAllocGpuNct)/(1024*1024));
}

//==============================================================================
/// Checks allocated memory for the indicated number of particles.
/// If the allocated memeory is not sufficient, reserve the required memory.
///
/// Comprueba la reserva de memoria para el numero indicado de particulas. 
/// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivGpu::CheckMemoryNp(unsigned npmin){
  if(SizeNp<npmin+PARTICLES_OVERMEMORY_MIN){
    AllocMemoryNp(ullong(npmin)+ullong(OverMemoryNp*npmin)+IncreaseNp);
    IncreaseNp=0;
  }
  else if(!CellPart){
    AllocMemoryNp(SizeNp+IncreaseNp);
    IncreaseNp=0;
  }
}

//==============================================================================
/// Checks allocated memory for the indicated number of cells.
/// If the allocated memeory is not sufficient, reserve the required memory.
///
/// Comprueba la reserva de memoria para el numero indicado de celdas. 
/// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivGpu::CheckMemoryNct(unsigned nctmin){
  if(SizeNct<nctmin){
    unsigned overnct=0;
    if(OverMemoryCells>0){
      ullong nct=ullong(Ncx+OverMemoryCells)*ullong(Ncy+OverMemoryCells)*ullong(Ncz+OverMemoryCells);
      ullong nctt=SizeBeginEndCell(nct);
      if(nctt!=unsigned(nctt))Run_Exceptioon("The number of cells is too big.");
      overnct=unsigned(nct);
    }
    AllocMemoryNct(nctmin>overnct? nctmin: overnct);
  }
  else if(!BeginEndCell)AllocMemoryNct(SizeNct);  
}

//==============================================================================
/// Defines the domain to be used by the simulation.
/// Define el dominio de simulacion a usar.
//==============================================================================
void JCellDivGpu::DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin,tdouble3 domposmin,tdouble3 domposmax){
  DomCellCode=cellcode;
  DomCelIni=domcelini;
  DomCelFin=domcelfin;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
  DomCells=DomCelFin-DomCelIni;
}

/*:/==============================================================================
// Devuelve coordenadas de celda a partir de una posicion.
//==============================================================================
//tuint3 JCellDivGpu::GetMapCell(const tfloat3 &pos)const{ pdtecell
//  float dx=pos.x-MapPosMin.x,dy=pos.y-MapPosMin.y,dz=pos.z-MapPosMin.z;
//  unsigned cx=unsigned(dx*OvScell),cy=unsigned(dy*OvScell),cz=unsigned(dz*OvScell);
//  return(TUint3(cx,cy,cz));
//}:*/

//==============================================================================
/// Computes the maximum and minimum postion for the indicated range of boundary particles.
/// The excluded particles are already marked in code[].
///
/// Calcula posiciones minimas y maximas del rango de particulas Bound indicado.
/// En code[] ya estan marcadas las particulas excluidas.
//==============================================================================
void JCellDivGpu::CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2,unsigned pini2
  ,const unsigned* dcellg,const typecode *codeg,tuint3 &cellmin,tuint3 &cellmax)
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
/// Computes the maximum and minimum postion for the indicated range of fluid particles.
/// Ignores the excluded particles already marked in code[]
///
/// Calcula posiciones minimas y maximas del rango de particulas Fluid indicado.
/// Ignora particulas excluidas que ya estan marcadas en code[].
//==============================================================================
void JCellDivGpu::CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2,unsigned pini2
  ,const unsigned* dcellg,const typecode *codeg,tuint3 &cellmin,tuint3 &cellmax)
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
void JCellDivGpu::SortBasicArrays(const unsigned *idp,const typecode *code,const unsigned *dcell,const double2 *posxy,const double *posz,const float4 *velrhop
  ,unsigned *idp2,typecode *code2,unsigned *dcell2,double2 *posxy2,double *posz2,float4 *velrhop2)
{
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,idp,code,dcell,posxy,posz,velrhop,idp2,code2,dcell2,posxy2,posz2,velrhop2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for type float4).
/// Ordena arrays de datos segun SortPart (para tipo float4).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float4 *a,float4 *a2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for float values).
/// Ordena arrays de datos segun SortPart (para valores float).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float *a,const float *b,float *a2,float *b2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,b,a2,b2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for double and double2 values).
/// Ordena arrays de datos segun SortPart (para valores double y double2).
//==============================================================================
void JCellDivGpu::SortDataArrays(const double2 *a,const double *b,const float4 *c,double2 *a2,double *b2,float4 *c2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,b,c,a2,b2,c2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for type tsymatrix3f).
/// Ordena arrays de datos segun SortPart (para tipo tsymatrix3f).
//==============================================================================
void JCellDivGpu::SortDataArrays(const tsymatrix3f *a,tsymatrix3f *a2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for type float3).
/// Ordena arrays de datos segun SortPart (para tipo float3).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float3 *a,float3 *a2){
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
}

//==============================================================================
/// Reorders data arrays according to SortPart (for float values).
/// Ordena arrays de datos segun SortPart (para valores float).
//==============================================================================
void JCellDivGpu::SortDataArrays(const float *a, float *a2) {
  const unsigned pini=(DivideFull? 0: NpbFinal);
  cudiv::SortDataParticles(Nptot,pini,SortPart,a,a2);
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
tdouble3 JCellDivGpu::GetDomainLimits(bool limitmin,unsigned slicecellmin)const{
  tuint3 celmin=GetCellDomainMin(),celmax=GetCellDomainMax();
  if(celmin.x>celmax.x)celmin.x=celmax.x=0; else celmax.x++;
  if(celmin.y>celmax.y)celmin.y=celmax.y=0; else celmax.y++;
  if(celmin.z>celmax.z)celmin.z=celmax.z=slicecellmin; else celmax.z++;
  double scell=double(Scell);
  tdouble3 pmin=DomPosMin+TDouble3(scell*celmin.x,scell*celmin.y,scell*celmin.z);
  tdouble3 pmax=DomPosMin+TDouble3(scell*celmax.x,scell*celmax.y,scell*celmax.z);
  return(limitmin? pmin: pmax);
}

/*:
////==============================================================================
//// Devuelve rango de particulas en el rango de celdas indicadas.
////==============================================================================
//uint2 JCellDivGpu::GetRangeParticlesCells(bool fluid,unsigned celini,unsigned celfin)const{
//  if(fluid){ celini+=BoxFluid; celfin+=BoxFluid; }
//  unsigned pmin=UINT_MAX,pmax=0;
//  if(celini<celfin){
//    bool memorynew=false;
//    unsigned *auxg=NULL;
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
//    unsigned *auxg=NULL;
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
void JCellDivGpu::DgSaveVktRange(std::string file,unsigned pini,unsigned pfin,const unsigned *idpg,const float3 *posg)const{
  int mpirank=Log->GetMpiRank();
  if(mpirank>=0)file=string("p")+fun::IntStr(mpirank)+"_"+file;
  file=DirOut+file;
  unsigned np=pfin-pini;
  tfloat3 *pos=new tfloat3[np];
  unsigned *idp=new unsigned[np];
  cudaMemcpy(idp,idpg+pini,sizeof(unsigned)*np,cudaMemcpyDeviceToHost);
  cudaMemcpy(pos,posg+pini,sizeof(float3)*np,cudaMemcpyDeviceToHost);
  JFormatFiles2::ParticlesToVtk(file,pfin-pini,pos,NULL,NULL,NULL,NULL,idp,NULL,NULL,NULL,NULL);
  delete[] pos;
  delete[] idp;
}:*/


