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

/// \file JCellDivCpu.cpp \brief Implements the class \ref JCellDivCpu.

#include "JCellDivCpu.h"
#include "JAppInfo.h"
#include "Functions.h"
#include <cfloat>
#include <climits>
#include <numeric>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCellDivCpu::JCellDivCpu(bool stable,bool floating,byte periactive
  ,bool celldomfixed,TpCellMode cellmode,float scell
  ,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
  ,unsigned casenbound,unsigned casenfixed,unsigned casenpb
  ,std::string dirout)
  :Log(AppInfo.LogPtr()),Stable(stable),Floating(floating),PeriActive(periactive)
  ,CellDomFixed(celldomfixed),CellMode(cellmode)
  ,ScellDiv(cellmode==CELLMODE_Full? 1: (cellmode==CELLMODE_Half? 2: 0))
  ,Scell(scell),OvScell(1.f/scell)
  ,Map_PosMin(mapposmin),Map_PosMax(mapposmax),Map_PosDif(mapposmax-mapposmin)
  ,Map_Cells(mapcells),CaseNbound(casenbound),CaseNfixed(casenfixed)
  ,CaseNpb(casenpb),DirOut(dirout)
  ,OverMemoryNp(CELLDIV_OVERMEMORYNP)
  ,OverMemoryCells(CELLDIV_OVERMEMORYCELLS)
  ,OverMemoryNCells(CELLDIV_OVERMEMORYNCELLS)
{
  ClassName="JCellDivCpu";
  CellPart=NULL;    SortPart=NULL;
  PartsInCell=NULL; BeginCell=NULL;
  VSort=NULL;
  SortIdx=NULL; //<vs_flexstruc>
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JCellDivCpu::~JCellDivCpu(){
  DestructorActive=true;
  //Log->Printf("---> DivideFull:%u/%u",NdivFull,Ndiv); //:del:
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JCellDivCpu::Reset(){
  SizeNp=SizeNct=0;
  IncreaseNp=0;
  FreeMemoryAll();
  Ndiv=NdivFull=0;
  Nptot=Npb1=Npf1=Npb2=Npf2=0;
  MemAllocNp=MemAllocNct=0;
  MemAllocNpTimes=MemAllocNctTimes=0;
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
void JCellDivCpu::FreeMemoryNct(){
  delete[] PartsInCell;   PartsInCell=NULL;
  delete[] BeginCell;     BeginCell=NULL; 
  MemAllocNct=0;
  BoundDivideOk=false;
}

//==============================================================================
/// Free memory reserved for particles.
/// Libera memoria reservada para particulas.
//==============================================================================
void JCellDivCpu::FreeMemoryNp(){
  delete[] CellPart;    CellPart=NULL;
  delete[] SortPart;    SortPart=NULL;
  delete[] VSort;       SetMemoryVSort(NULL);
  delete[] SortIdx;     SortIdx=NULL;   //<vs_flexstruc>
  MemAllocNp=0;
  BoundDivideOk=false;
}

//==============================================================================
/// Free memory reserved for particles and cells.
/// Libera memoria reservada para particulas y celdas.
//==============================================================================
void JCellDivCpu::FreeMemoryAll(){
  FreeMemoryNct();
  FreeMemoryNp();
}

//==============================================================================
/// Adjust buffers to reorder particle values.
/// Ajusta buffers para reordenar datos de particulas.
//==============================================================================
void JCellDivCpu::SetMemoryVSort(byte* vsort){
  VSort=vsort;
  VSortInt=(int*)VSort;        VSortWord=(word*)VSort;
  VSortFloat=(float*)VSort;    VSortFloat3=(tfloat3*)VSort;
  VSortFloat4=(tfloat4*)VSort; VSortDouble3=(tdouble3*)VSort;
  VSortSymmatrix3f=(tsymatrix3f*)VSort;
}

//==============================================================================
/// Assign memory according to number of particles. 
/// Asigna memoria segun numero de particulas. 
//==============================================================================
void JCellDivCpu::AllocMemoryNp(ullong np,ullong npmin){
  FreeMemoryNp();
  np=np+PARTICLES_OVERMEMORY_MIN;
  SizeNp=unsigned(np);
  //-Check number of particles | Comprueba numero de particulas.
  if(np!=SizeNp)Run_Exceptioon(string("Failed memory allocation for ")+fun::UlongStr(np)+" particles.");
  //-Reserve memory for particles | Reserva memoria para particulas.
  MemAllocNp=0;
  try{
    CellPart=new unsigned[SizeNp];                      MemAllocNp+=sizeof(unsigned)*SizeNp;
    SortPart=new unsigned[SizeNp];                      MemAllocNp+=sizeof(unsigned)*SizeNp;
    SetMemoryVSort(new byte[sizeof(tdouble3)*SizeNp]);  MemAllocNp+=sizeof(tdouble3)*SizeNp;
    SortIdx=new unsigned[SizeNp];                       MemAllocNp+=sizeof(unsigned)*SizeNp;  //<vs_flexstruc>
    MemAllocNpTimes++;
  }
  catch(const std::bad_alloc){
    Run_Exceptioon(fun::PrintStr("Failed CPU memory allocation of %.1f MB for %u particles."
      ,double(MemAllocNp)/MEBIBYTE,SizeNp));
  }
  //-Show requested memory.
  const string txover=(npmin>1? fun::PrintStr(" (over-allocation: %.2fX)",double(SizeNp)/npmin): "");
  Log->Printf("**CellDiv: Requested cpu memory for %s particles%s: %.1f MiB (%u times)."
    ,KINT(SizeNp),txover.c_str(),double(MemAllocNp)/MEBIBYTE,MemAllocNpTimes);
}

//==============================================================================
/// Assign memory according to number of cells. 
/// Asigna memoria segun numero de celdas. 
//==============================================================================
void JCellDivCpu::AllocMemoryNct(ullong nct,ullong nctmin){
  FreeMemoryNct();
  SizeNct=unsigned(nct);
  const unsigned nctt=unsigned(SizeBeginCell(SizeNct));
  //-Checks number of cells.
  if(ullong(nctt)!=SizeBeginCell(nct))Run_Exceptioon(
    fun::PrintStr("Failed CPU memory allocation for %s cells.",KINT(nct)));
  //-Allocates memory for cells.
  MemAllocNct=0;
  try{
    PartsInCell=new unsigned[nctt-1];  MemAllocNct+=sizeof(unsigned)*(nctt-1);
    BeginCell  =new unsigned[nctt];    MemAllocNct+=sizeof(unsigned)*(nctt);
    MemAllocNctTimes++;
  }
  catch(const std::bad_alloc){
    Run_Exceptioon(fun::PrintStr("Failed CPU memory allocation of %.1f MiB for %s cells."
      ,double(MemAllocNct)/MEBIBYTE,KINT(SizeNct)));
  }
  //-Shows requested memory.
  const string txover=(nctmin>1? fun::PrintStr(" (over-allocation: %.2fX)",double(SizeNct)/nctmin): "");
  Log->Printf("**CellDiv: Requested CPU memory for %s cells%s: %.1f MiB (%u times)."
    ,KINT(SizeNct),txover.c_str(),double(MemAllocNct)/MEBIBYTE,MemAllocNctTimes);
}

//==============================================================================
/// Check reserved memory for the indicated number of particles. 
/// If there is insufficient memory or it is not reserved, then reserve the requested memory.
///
/// Comprueba la reserva de memoria para el numero indicado de particulas. 
/// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivCpu::CheckMemoryNp(unsigned npmin){
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
/// Check reserved memory for the indicated number of cells. 
/// If there is insufficient memory or it is not reserved, then reserve the requested memory.
///
/// Comprueba la reserva de memoria para el numero indicado de celdas. 
/// Si no es suficiente o no hay reserva, entonces reserva la memoria requerida.
//==============================================================================
void JCellDivCpu::CheckMemoryNct(unsigned nctmin){
  if(SizeNct<nctmin){
    unsigned nctnew=nctmin;
    if(OverMemoryCells>0){
      const ullong nct1=ullong(Ncx+OverMemoryCells)*ullong(Ncy+OverMemoryCells)*ullong(Ncz+OverMemoryCells);
      const ullong nct2=ullong(nctmin)+OverMemoryNCells;
      const ullong nct3=(nct1>nct2? nct1: nct2);
      const ullong nct4=min(DomCellsMax,nct3);
      if(SizeBeginCell(nct4)>=UINT_MAX)Run_Exceptioon("The number of cells is too big.");
      nctnew=unsigned(nct4);
    }
    AllocMemoryNct(nctnew,nctmin);
  }
  else if(!BeginCell)AllocMemoryNct(SizeNct,nctmin);  
}

//==============================================================================
/// Define simulation domain to use.
/// Define el dominio de simulacion a usar.
//==============================================================================
void JCellDivCpu::DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin
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

/*:
//==============================================================================
// Devuelve coordenadas de celda a partir de una posicion.
//==============================================================================
//tuint3 JCellDivCpu::GetMapCell(const tfloat3& pos)const{
//  float dx=pos.x-MapPosMin.x,dy=pos.y-MapPosMin.y,dz=pos.z-MapPosMin.z;
//  unsigned cx=unsigned(dx*OvScell),cy=unsigned(dy*OvScell),cz=unsigned(dz*OvScell);
//  return(TUint3(cx,cy,cz));
//}:*/

//==============================================================================
/// Computes minimum and maximum cells for valid particles.
/// The excluded particles are already marked in code[]
/// In case of there being no valid particles, the minimum is set to be greater than the maximum.
///
/// Calcula celda minima y maxima de las particulas validas.
/// Las particulas excluidas ya estan marcadas en code[].
/// En caso de no haber ninguna particula valida el minimo sera mayor que el maximo.
//==============================================================================
void JCellDivCpu::LimitsCellBound(unsigned n,unsigned pini,const unsigned* dcellc
  ,const typecode* codec,tuint3& cellmin,tuint3& cellmax)const
{
  tuint3 cmin=TUint3(UINT_MAX);
  tuint3 cmax=TUint3(0);
  const unsigned pfin=pini+n;
  for(unsigned p=pini;p<pfin;p++){
    const unsigned rcell=dcellc[p];
    const unsigned cx=DCEL_Cellx(DomCellCode,rcell);
    const unsigned cy=DCEL_Celly(DomCellCode,rcell);
    const unsigned cz=DCEL_Cellz(DomCellCode,rcell);
    const typecode rcode=codec[p];
    const typecode rcodsp=CODE_GetSpecialValue(rcode);
    if(rcodsp<CODE_OUTIGNORE){ //-Particle not excluded | Particula no excluida.
      if(cmin.x>cx)cmin.x=cx;
      if(cmin.y>cy)cmin.y=cy;
      if(cmin.z>cz)cmin.z=cz;
      if(cmax.x<cx)cmax.x=cx;
      if(cmax.y<cy)cmax.y=cy;
      if(cmax.z<cz)cmax.z=cz;
    }
  }
  cellmin=cmin;
  cellmax=cmax;
}

//==============================================================================
/// Computes the maximum and minimum position for the indicated range of boundary particles.
/// The excluded particles are already marked in code[].
///
/// Calcula posiciones minimas y maximas del rango de particulas Bound indicado.
/// En code[] ya estan marcadas las particulas excluidas.
//==============================================================================
void JCellDivCpu::CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2,unsigned pini2
  ,const unsigned* dcellc,const typecode* codec,tuint3& cellmin,tuint3& cellmax)
{
  tuint3 cmin,cmax;
  LimitsCellBound(n,pini,dcellc,codec,cmin,cmax);
  cellmin=(cmin.x>cmax.x? DomCells: cmin);
  cellmax=(cmin.x>cmax.x? TUint3(0): cmax);
  if(n2){
    LimitsCellBound(n2,pini2,dcellc,codec,cmin,cmax);
    cmin=(cmin.x>cmax.x? DomCells: cmin);
    cmax=(cmin.x>cmax.x? TUint3(0): cmax);
    cellmin=MinValues(cellmin,cmin);
    cellmax=MaxValues(cellmax,cmax);
  }
}

//==============================================================================
/// Computes minimum and maximum cell for valid particles.
/// The excluded particles are already marked in code[].
/// In case of having no valid particles the minimum value igreater than the maximum.
///
/// Calcula celda minima y maxima de las particulas validas.
/// Las particulas excluidas ya estan marcadas en code[].
/// En caso de no haber ninguna particula valida el minimo sera mayor que el maximo.
//==============================================================================
void JCellDivCpu::LimitsCellFluid(unsigned n,unsigned pini,const unsigned* dcellc
  ,const typecode* codec,tuint3& cellmin,tuint3& cellmax)const
{
  tuint3 cmin=TUint3(UINT_MAX);
  tuint3 cmax=TUint3(0);
  const unsigned pfin=pini+n;
  for(unsigned p=pini;p<pfin;p++){
    const unsigned rcell=dcellc[p];
    const unsigned cx=DCEL_Cellx(DomCellCode,rcell);
    const unsigned cy=DCEL_Celly(DomCellCode,rcell);
    const unsigned cz=DCEL_Cellz(DomCellCode,rcell);
    if(CODE_GetSpecialValue(codec[p])<CODE_OUTIGNORE){ //-Particle not excluded | Particula no excluida.
      if(cmin.x>cx)cmin.x=cx;
      if(cmin.y>cy)cmin.y=cy;
      if(cmin.z>cz)cmin.z=cz;
      if(cmax.x<cx)cmax.x=cx;
      if(cmax.y<cy)cmax.y=cy;
      if(cmax.z<cz)cmax.z=cz;
    }
  }
  cellmin=cmin;
  cellmax=cmax;
}

//==============================================================================
/// Computes the maximum and minimum position for the indicated range of fluid particles.
/// Ignores the excluded particles already marked in code[]
///
/// Calcula posiciones minimas y maximas del rango de particulas Fluid indicado.
/// Ignora particulas excluidas que ya estan marcadas en code[].
//==============================================================================
void JCellDivCpu::CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2,unsigned pini2
  ,const unsigned* dcellc,const typecode* codec,tuint3& cellmin,tuint3& cellmax)
{
  tuint3 cmin,cmax;
  LimitsCellFluid(n,pini,dcellc,codec,cmin,cmax);
  cellmin=(cmin.x>cmax.x? DomCells: cmin);
  cellmax=(cmin.x>cmax.x? TUint3(0): cmax);
  if(n2){
    LimitsCellFluid(n2,pini2,dcellc,codec,cmin,cmax);
    cmin=(cmin.x>cmax.x? DomCells: cmin);
    cmax=(cmin.x>cmax.x? TUint3(0): cmax);
    cellmin=MinValues(cellmin,cmin);
    cellmax=MaxValues(cellmax,cmax);
  }
  //:Log->Printf("CalcDomainFluid> cell:(%s)-(%s)",fun::Uint3Str(cellmin).c_str(),fun::Uint3Str(cellmax).c_str());
}

//==============================================================================
/// Reorder values of all particles (for type word).
/// Reordena datos de todas las particulas (para tipo word).
//==============================================================================
void JCellDivCpu::SortArray(word* vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortWord[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortWord+ini,sizeof(word)*(n-ini));
}

//==============================================================================
/// Reorder values of all particles (for type unsigned).
/// Reordena datos de todas las particulas (para tipo unsigned).
//==============================================================================
void JCellDivCpu::SortArray(unsigned* vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortInt[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortInt+ini,sizeof(unsigned)*(n-ini));
}

//==============================================================================
/// Reorder values of all particles (for type float).
/// Reordena datos de todas las particulas (para tipo float).
//==============================================================================
void JCellDivCpu::SortArray(float* vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortFloat[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortFloat+ini,sizeof(float)*(n-ini));
}

//==============================================================================
/// Reorder values of all particles (for type tdouble3).
/// Reordena datos de todas las particulas (para tipo tdouble3).
//==============================================================================
void JCellDivCpu::SortArray(tdouble3* vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortDouble3[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortDouble3+ini,sizeof(tdouble3)*(n-ini));
}

//==============================================================================
/// Reorder values of all particles (for type tfloat3).
/// Reordena datos de todas las particulas (para tipo tfloat3).
//==============================================================================
void JCellDivCpu::SortArray(tfloat3* vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortFloat3[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortFloat3+ini,sizeof(tfloat3)*(n-ini));
}

//==============================================================================
/// Reorder values of all particles (for type tfloat4).
/// Reordena datos de todas las particulas (para tipo tfloat4).
//==============================================================================
void JCellDivCpu::SortArray(tfloat4* vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortFloat4[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortFloat4+ini,sizeof(tfloat4)*(n-ini));
}

//==============================================================================
/// Reorder values of all particles (for type tsymatrix3f).
/// Reordena datos de todas las particulas (para tipo tsymatrix3f).
//==============================================================================
void JCellDivCpu::SortArray(tsymatrix3f* vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortSymmatrix3f[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortSymmatrix3f+ini,sizeof(tsymatrix3f)*(n-ini));
}

//==============================================================================
/// Reorder PeriParent references.
//==============================================================================
void JCellDivCpu::SortArrayPeriParent(unsigned* vec){
  memset(VSortInt,255,sizeof(unsigned)*Nptot);
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  unsigned* resortpart=(unsigned*)(VSortInt+Nptot);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)resortpart[SortPart[p]]=p;
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++){
    const unsigned pp=vec[SortPart[p]];
    VSortInt[p]=(pp!=UINT_MAX? resortpart[pp]: pp);
  }
  memcpy(vec+ini,VSortInt+ini,sizeof(unsigned)*(n-ini));
}

//==============================================================================
/// Return current limites of domain.
/// Devuelve limites actuales del dominio.
//==============================================================================
tdouble6 JCellDivCpu::GetDomainLimitsMinMax(unsigned slicecellmin)const{
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
/// Return current limites of domain.
/// Devuelve limites actuales del dominio.
//==============================================================================
tdouble3 JCellDivCpu::GetDomainLimits(bool limitmin,unsigned slicecellmin)const{
  const tdouble6 limits=GetDomainLimitsMinMax(slicecellmin);
  return(limitmin? limits.getlo(): limits.gethi());
}

//==============================================================================
/// Returns cell division data for neighborhood search.
/// Devuelve datis de division en celdas para busqueda de vecinos.
//==============================================================================
StDivDataCpu JCellDivCpu::GetCellDivData()const{
  return(MakeDivDataCpu(ScellDiv,GetNcells(),GetCellDomainMin(),GetBeginCell()
    ,Scell,DomCellCode,DomPosMin));
}

//<vs_flexstruc_ini>
//==============================================================================
/// Sorts and updates the indices of the flexible structure particles.
/// Ordena y actualiza los índices de las partículas de estructura flexible.
//==============================================================================
void JCellDivCpu::UpdateIndices(unsigned n,unsigned* idx){
  if(DivideFull){
    SortIndices(SortPart,SortIdx,NpbFinal,Stable);
    UpdateIndices(n,SortIdx,idx);
  }
}

//==============================================================================
/// Sorts the indices of the flexible structure particles.
/// Ordena los índices de las partículas de estructura flexible.
//==============================================================================
void JCellDivCpu::SortIndices(unsigned* sortpart,unsigned* sortidx,unsigned np,bool stable){
  iota(sortidx,sortidx+np,0);
  if(stable)stable_sort(sortidx,sortidx+np,[sortpart](unsigned i1,unsigned i2){return sortpart[i1]<sortpart[i2];});
  else             sort(sortidx,sortidx+np,[sortpart](unsigned i1,unsigned i2){return sortpart[i1]<sortpart[i2];});
}

//==============================================================================
/// Updates the indices of the flexible structure particles.
/// Actualiza los índices de las partículas de estructura flexible..
//==============================================================================
void JCellDivCpu::UpdateIndices(unsigned np,const unsigned* sortidx,unsigned* idx){
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=0;p<n;p++)if(p<n)idx[p]=sortidx[idx[p]];
}
//<vs_flexstruc_end>

/*:
////==============================================================================
//// Indica si la celda esta vacia o no.
////==============================================================================
//bool JCellDivCpu::CellNoEmpty(unsigned box,byte kind)const{
//#ifdef DBG_JCellDivCpu
//  if(box>=Nct)Run_Exceptioon("Celda no valida.");
//#endif
//  if(kind==2)box+=BoxFluid;
//  return(BeginCell[box]<BeginCell[box+1]);
//}

////==============================================================================
//// Devuelve la primera particula de la celda solicitada.
////==============================================================================
//unsigned JCellDivCpu::CellBegin(unsigned box,byte kind)const{
//#ifdef DBG_JCellDivCpu
//  if(box>Nct)Run_Exceptioon("Celda no valida.");
//#endif
//  return(BeginCell[(kind==1? box: box+BoxFluid)]);
//}

////==============================================================================
//// Devuelve el numero de particulas de la celda solicitada.
////==============================================================================
//unsigned JCellDivCpu::CellSize(unsigned box,byte kind)const{
//#ifdef DBG_JCellDivCpu
//  if(box>Nct)Run_Exceptioon("Celda no valida.");
//#endif
//  if(kind==2)box+=BoxFluid;
//  return(BeginCell[box+1]-BeginCell[box]);
//}
:*/


