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

/// \file JCellDivCpu.cpp \brief Implements the class \ref JCellDivCpu.

#include "JCellDivCpu.h"
#include "JAppInfo.h"
#include "Functions.h"
#include <cfloat>
#include <climits>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCellDivCpu::JCellDivCpu(bool stable,bool floating,byte periactive
  ,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
  ,unsigned casenbound,unsigned casenfixed,unsigned casenpb,std::string dirout
  ,bool allocfullnct,float overmemorynp,word overmemorycells)
  :Log(AppInfo.LogPtr()),Stable(stable),Floating(floating),PeriActive(periactive)
  ,CellMode(cellmode),ScellDiv(cellmode==CELLMODE_Full? 1: (cellmode==CELLMODE_Half? 2: 0)),Scell(scell)
  ,OvScell(1.f/scell),Map_PosMin(mapposmin),Map_PosMax(mapposmax),Map_PosDif(mapposmax-mapposmin)
  ,Map_Cells(mapcells),CaseNbound(casenbound),CaseNfixed(casenfixed),CaseNpb(casenpb)
  ,DirOut(dirout),AllocFullNct(allocfullnct),OverMemoryNp(overmemorynp),OverMemoryCells(overmemorycells)
{
  ClassName="JCellDivCpu";
  CellPart=NULL;    SortPart=NULL;
  PartsInCell=NULL; BeginCell=NULL;
  VSort=NULL;
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
void JCellDivCpu::SetMemoryVSort(byte *vsort){
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
void JCellDivCpu::AllocMemoryNp(ullong np){
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
  }
  catch(const std::bad_alloc){
    Run_Exceptioon(fun::PrintStr("Failed CPU memory allocation of %.1f MB for %u particles.",double(MemAllocNp)/(1024*1024),SizeNp));
  }
  //-Show requested memory | Muestra la memoria solicitada.
  Log->Printf("**CellDiv: Requested cpu memory for %u particles: %.1f MB.",SizeNp,double(MemAllocNp)/(1024*1024));
}

//==============================================================================
/// Assign memory according to number of cells. 
/// Asigna memoria segun numero de celdas. 
//==============================================================================
void JCellDivCpu::AllocMemoryNct(ullong nct){
  FreeMemoryNct();
  SizeNct=unsigned(nct);
  //-Check number of cells | Comprueba numero de celdas.
  if(nct!=SizeNct)Run_Exceptioon(string("Failed GPU memory allocation for ")+fun::UlongStr(nct)+" cells.");
  //-Reserve memory for cells | Reserva memoria para celdas.
  MemAllocNct=0;
  const unsigned nc=(unsigned)SizeBeginCell(nct);
  try{
    PartsInCell=new unsigned[nc-1];  MemAllocNct+=sizeof(unsigned)*(nc-1);
    BeginCell=new unsigned[nc];      MemAllocNct+=sizeof(unsigned)*(nc);
  }
  catch(const std::bad_alloc){
    Run_Exceptioon(fun::PrintStr("Failed CPU memory allocation of %.1f MB for %u cells.",double(MemAllocNct)/(1024*1024),SizeNct));
  }
  //-Show requested memory | Muestra la memoria solicitada.
  Log->Printf("**CellDiv: Requested cpu memory for %u cells (CellMode=%s): %.1f MB.",SizeNct,GetNameCellMode(CellMode),double(MemAllocNct)/(1024*1024));
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
    AllocMemoryNp(ullong(npmin)+ullong(OverMemoryNp*npmin)+IncreaseNp);
    IncreaseNp=0;
  }
  else if(!CellPart){
    AllocMemoryNp(SizeNp+IncreaseNp);
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
    unsigned overnct=0;
    if(OverMemoryCells>0){
      ullong nct=ullong(Ncx+OverMemoryCells)*ullong(Ncy+OverMemoryCells)*ullong(Ncz+OverMemoryCells);
      ullong nctt=SizeBeginCell(nct);
      if(nctt!=unsigned(nctt))Run_Exceptioon("The number of cells is too big.");
      overnct=unsigned(nct);
    }
    AllocMemoryNct(nctmin>overnct? nctmin: overnct);
  }
  else if(!BeginCell)AllocMemoryNct(SizeNct);  
}

//==============================================================================
/// Define simulation domain to use.
/// Define el dominio de simulacion a usar.
//==============================================================================
void JCellDivCpu::DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin,tdouble3 domposmin,tdouble3 domposmax){
  DomCellCode=cellcode;
  DomCelIni=domcelini;
  DomCelFin=domcelfin;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
  DomCells=DomCelFin-DomCelIni;
}

/*:
//==============================================================================
// Devuelve coordenadas de celda a partir de una posicion.
//==============================================================================
//tuint3 JCellDivCpu::GetMapCell(const tfloat3 &pos)const{
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
  ,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax)const
{
  tuint3 cmin=TUint3(UINT_MAX);
  tuint3 cmax=TUint3(0);
  const unsigned pfin=pini+n;
  for(unsigned p=pini;p<pfin;p++){
    const unsigned rcell=dcellc[p];
    const unsigned cx=PC__Cellx(DomCellCode,rcell);
    const unsigned cy=PC__Celly(DomCellCode,rcell);
    const unsigned cz=PC__Cellz(DomCellCode,rcell);
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
/// Computes the maximum and minimum postion for the indicated range of boundary particles.
/// The excluded particles are already marked in code[].
///
/// Calcula posiciones minimas y maximas del rango de particulas Bound indicado.
/// En code[] ya estan marcadas las particulas excluidas.
//==============================================================================
void JCellDivCpu::CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2,unsigned pini2
  ,const unsigned* dcellc,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax)
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
/// Computes minimun and maximum cell for valid particles.
/// The excluded particles are already marked in code[].
/// In case of having no valid particles the minimum value igreater than the maximum.
///
/// Calcula celda minima y maxima de las particulas validas.
/// Las particulas excluidas ya estan marcadas en code[].
/// En caso de no haber ninguna particula valida el minimo sera mayor que el maximo.
//==============================================================================
void JCellDivCpu::LimitsCellFluid(unsigned n,unsigned pini,const unsigned* dcellc
  ,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax)const
{
  tuint3 cmin=TUint3(UINT_MAX);
  tuint3 cmax=TUint3(0);
  const unsigned pfin=pini+n;
  for(unsigned p=pini;p<pfin;p++){
    const unsigned rcell=dcellc[p];
    const unsigned cx=PC__Cellx(DomCellCode,rcell);
    const unsigned cy=PC__Celly(DomCellCode,rcell);
    const unsigned cz=PC__Cellz(DomCellCode,rcell);
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
/// Computes the maximum and minimum postion for the indicated range of fluid particles.
/// Ignores the excluded particles already marked in code[]
///
/// Calcula posiciones minimas y maximas del rango de particulas Fluid indicado.
/// Ignora particulas excluidas que ya estan marcadas en code[].
//==============================================================================
void JCellDivCpu::CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2,unsigned pini2
  ,const unsigned* dcellc,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax)
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
void JCellDivCpu::SortArray(word *vec){
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
void JCellDivCpu::SortArray(unsigned *vec){
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
void JCellDivCpu::SortArray(float *vec){
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
void JCellDivCpu::SortArray(tdouble3 *vec){
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
void JCellDivCpu::SortArray(tfloat3 *vec){
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
void JCellDivCpu::SortArray(tfloat4 *vec){
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
void JCellDivCpu::SortArray(tsymatrix3f *vec){
  const int n=int(Nptot);
  const int ini=(DivideFull? 0: int(NpbFinal));
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ini;p<n;p++)VSortSymmatrix3f[p]=vec[SortPart[p]];
  memcpy(vec+ini,VSortSymmatrix3f+ini,sizeof(tsymatrix3f)*(n-ini));
}

//==============================================================================
/// Return current limites of domain.
/// Devuelve limites actuales del dominio.
//==============================================================================
tdouble3 JCellDivCpu::GetDomainLimits(bool limitmin,unsigned slicecellmin)const{
  tuint3 celmin=GetCellDomainMin(),celmax=GetCellDomainMax();
  if(celmin.x>celmax.x)celmin.x=celmax.x=0; else celmax.x++;
  if(celmin.y>celmax.y)celmin.y=celmax.y=0; else celmax.y++;
  if(celmin.z>celmax.z)celmin.z=celmax.z=slicecellmin; else celmax.z++;
  double scell=double(Scell);
  tdouble3 pmin=DomPosMin+TDouble3(scell*celmin.x,scell*celmin.y,scell*celmin.z);
  tdouble3 pmax=DomPosMin+TDouble3(scell*celmax.x,scell*celmax.y,scell*celmax.z);
  return(limitmin? pmin: pmax);
}

//==============================================================================
/// Returns cell division data for neighborhood search.
/// Devuelve datis de division en celdas para busqueda de vecinos.
//==============================================================================
StDivDataCpu JCellDivCpu::GetCellDivData()const{
  return(MakeDivDataCpu(ScellDiv,GetNcells(),GetCellDomainMin(),GetBeginCell()
    ,Scell,DomCellCode,DomPosMin));
}

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


