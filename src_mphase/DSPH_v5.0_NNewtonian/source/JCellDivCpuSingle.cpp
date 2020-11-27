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

/// \file JCellDivCpuSingle.cpp \brief Implements the class \ref JCellDivCpuSingle.

#include "JCellDivCpuSingle.h"
#include "Functions.h"
#include <climits>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCellDivCpuSingle::JCellDivCpuSingle(bool stable,bool floating,byte periactive
  ,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
  ,unsigned casenbound,unsigned casenfixed,unsigned casenpb,std::string dirout)
  :JCellDivCpu(stable,floating,periactive,cellmode,scell,mapposmin,mapposmax,mapcells
  ,casenbound,casenfixed,casenpb,dirout)
{
  ClassName="JCellDivCpuSingle";
}

//==============================================================================
/// Calculate limits of domain in cells adjusting to fluid (CellDomainMin/Max). 
/// Calcula limites del dominio en celdas ajustando al fluido (CellDomainMin/Max). 
//==============================================================================
void JCellDivCpuSingle::CalcCellDomain(const unsigned *dcellc,const typecode *codec){
  //-Calculate boundary domain | Calcula dominio del contorno.
  tuint3 celbmin,celbmax;
  if(!BoundLimitOk){
    CalcCellDomainBound(Npb1,0,Npb2,Npb1+Npf1,dcellc,codec,celbmin,celbmax);
    BoundLimitOk=true; BoundLimitCellMin=celbmin; BoundLimitCellMax=celbmax;
  } 
  else{ celbmin=BoundLimitCellMin; celbmax=BoundLimitCellMax; }
  //Log->Printf("----->CalcCellDomain> BoundLimitCellMin/Max2:%s",fun::Uint3RangeStr(BoundLimitCellMin,BoundLimitCellMax).c_str());
  //-Calculate fluid domain | Calcula dominio del fluido.
  tuint3 celfmin,celfmax;
  CalcCellDomainFluid(Npf1,Npb1,Npf2,Npb1+Npf1+Npb2,dcellc,codec,celfmin,celfmax);
  //Log->Printf("----->CalcCellDomain> celfmin/max:%s",fun::Uint3RangeStr(celfmin,celfmax).c_str());
  //-Computes the domain adjusting to the boundary and the fluid (with KernelSize halo).
  //-Calcula dominio ajustando al contorno y al fluido (con halo de KernelSize). 
  MergeMapCellBoundFluid(celbmin,celbmax,celfmin,celfmax,CellDomainMin,CellDomainMax);
}

//==============================================================================
/// Combines cell limits of boundary and fluid with map limits.
/// If UseFluidDomain=TRUE, uses fluid domain plus KernelSize if there is a boundary;
/// if not, uses the fluid and boundary domain
/// If the domain is null CellDomainMin=CellDomainMax=(0,0,0).
///
/// Combina limite de celdas de contorno y fluido con limites de mapa.
/// Con UseFluidDomain=TRUE se queda con el dominio del fluido mas KernelSize 
/// si hay contorno, en caso contrario se queda con el dominio que incluya 
/// fluido y contorno.
/// En caso de que el dominio sea nulo CellDomainMin=CellDomainMax=(0,0,0).
//==============================================================================
void JCellDivCpuSingle::MergeMapCellBoundFluid(const tuint3 &celbmin,const tuint3 &celbmax,const tuint3 &celfmin,const tuint3 &celfmax,tuint3 &celmin,tuint3 &celmax)const{
  const unsigned scelldiv=unsigned(ScellDiv);
  celmin=TUint3(max(min(celbmin.x,celfmin.x),(celfmin.x>=scelldiv? celfmin.x-scelldiv: 0))
               ,max(min(celbmin.y,celfmin.y),(celfmin.y>=scelldiv? celfmin.y-scelldiv: 0))
               ,max(min(celbmin.z,celfmin.z),(celfmin.z>=scelldiv? celfmin.z-scelldiv: 0)));
  celmax=TUint3(min(max(celbmax.x,celfmax.x),celfmax.x+scelldiv)
               ,min(max(celbmax.y,celfmax.y),celfmax.y+scelldiv)
               ,min(max(celbmax.z,celfmax.z),celfmax.z+scelldiv));
  if(celmax.x>=DomCells.x)celmax.x=DomCells.x-1;
  if(celmax.y>=DomCells.y)celmax.y=DomCells.y-1;
  if(celmax.z>=DomCells.z)celmax.z=DomCells.z-1;
  if(celmin.x>celmax.x || celmin.y>celmax.y || celmin.z>celmax.z)celmin=celmax=TUint3(0,0,0);
}

//==============================================================================
/// Computes number of cells starting from (CellDomainMIN/Max).
/// Obtains location of special cells.
///
/// Calcula numero de celdas a partir de (CellDomainMin/Max). 
/// Obtiene localizacion de celdas especiales.
//==============================================================================
void JCellDivCpuSingle::PrepareNct(){
  //-Computes number of cells.
  Ncx=CellDomainMax.x-CellDomainMin.x+1;
  Ncy=CellDomainMax.y-CellDomainMin.y+1;
  Ncz=CellDomainMax.z-CellDomainMin.z+1;
  //:printf("======  ncx:%u ncy:%u ncz:%u\n",Ncx,Ncy,Ncz);
  Nsheet=Ncx*Ncy; Nct=Nsheet*Ncz; Nctt=SizeBeginCell(Nct);
  if(Nctt!=unsigned(Nctt))Run_Exceptioon("The number of cells is too big.");
  BoxBoundIgnore=Nct; 
  BoxFluid=BoxBoundIgnore+1; 
  BoxBoundOut=BoxFluid+Nct; 
  BoxFluidOut=BoxBoundOut+1; 
  BoxBoundOutIgnore=BoxFluidOut+1;
  BoxFluidOutIgnore=BoxBoundOutIgnore+1;
  //:Log->Printf("--->PrepareNct> BoxIgnore:%u BoxFluid:%u BoxBoundOut:%u BoxFluidOut:%u",BoxIgnore,BoxFluid,BoxBoundOut,BoxFluidOut);
  //:Log->Printf("--->PrepareNct> BoxBoundOutIgnore:%u BoxFluidOutIgnore:%u",BoxBoundOutIgnore,BoxFluidOutIgnore);
}

//==============================================================================
/// Computes cell of each boundary and fluid particle (cellpart[]) starting from its cell in 
/// the map. all the excluded particles were already marked in code[].
/// Excluded particles bound (fixed and moving) and floating are moved to BoxBoundOut.
/// Account for particles for cell (partsincell[]).
///
/// Calcula celda de cada particula bound y fluid (cellpart[]) a partir de su celda en
/// mapa. Todas las particulas excluidas ya fueron marcadas en code[].
/// Las particulas excluidas de tipo bound (fixed and moving) and floating se mueven a BoxBoundOut.
/// Contabiliza particulas por celda (partsincell[]).
//==============================================================================
void JCellDivCpuSingle::PreSortFull(unsigned np,const unsigned *dcellc,const typecode *codec
  ,unsigned* cellpart,unsigned* partsincell)const
{
  memset(partsincell,0,sizeof(unsigned)*(Nctt-1));
  for(unsigned p=0;p<np;p++){
    //-Computes cell according position.
    const unsigned rcell=dcellc[p];
    const unsigned cx=PC__Cellx(DomCellCode,rcell)-CellDomainMin.x;
    const unsigned cy=PC__Celly(DomCellCode,rcell)-CellDomainMin.y;
    const unsigned cz=PC__Cellz(DomCellCode,rcell)-CellDomainMin.z;
    const unsigned cellsort=cx+cy*Ncx+cz*Nsheet;
    //-Checks particle code.
    const typecode rcode=codec[p];
    const typecode codetype=CODE_GetType(rcode);
    const typecode codeout=CODE_GetSpecialValue(rcode);
    //-Assigns box.
    unsigned box;
    if(codetype<CODE_TYPE_FLOATING){//-Bound particles (except floating) | Particulas bound (excepto floating).
      box=(codeout<CODE_OUTIGNORE?   ((cx<Ncx && cy<Ncy && cz<Ncz)? cellsort: BoxBoundIgnore):   (codeout==CODE_OUTIGNORE? BoxBoundOutIgnore: BoxBoundOut));
    }
    else{//-Fluid and floating particles | Particulas fluid y floating.
      box=(codeout<=CODE_OUTIGNORE?   (codeout<CODE_OUTIGNORE? BoxFluid+cellsort: BoxFluidOutIgnore):   (codetype==CODE_TYPE_FLOATING? BoxBoundOut: BoxFluidOut));
    }
    cellpart[p]=box;
    partsincell[box]++;
  }
}

//==============================================================================
/// Computes cell of each fluid particle (cellpart[]) starting from its cell in 
/// the map. all the excluded particles were already marked in code[].
/// Excluded particles floating are moved to BoxBoundOut.
/// Account for particles for cell (partsincell[]).
///
/// Calcula celda de cada particula fluid (cellpart[]) a partir de su celda en
/// mapa. Todas las particulas excluidas ya fueron marcadas en code[].
/// Las particulas excluidas de tipo floating se mueven a BoxBoundOut.
/// Contabiliza particulas por celda (partsincell[]).
//==============================================================================
void JCellDivCpuSingle::PreSortFluid(unsigned np,unsigned pini,const unsigned *dcellc
  ,const typecode *codec,unsigned* cellpart,unsigned* partsincell)const
{
  memset(partsincell+BoxFluid,0,sizeof(unsigned)*(Nctt-1-BoxFluid));
  const unsigned pfin=pini+np;
  for(unsigned p=pini;p<pfin;p++){
    //-Computes cell according position.
    const unsigned rcell=dcellc[p];
    const unsigned cx=PC__Cellx(DomCellCode,rcell)-CellDomainMin.x;
    const unsigned cy=PC__Celly(DomCellCode,rcell)-CellDomainMin.y;
    const unsigned cz=PC__Cellz(DomCellCode,rcell)-CellDomainMin.z;
    const unsigned cellsortfluid=BoxFluid+cx+cy*Ncx+cz*Nsheet;
    //-Checks particle code.
    const typecode rcode=codec[p];
    const typecode codetype=CODE_GetType(rcode);
    const typecode codeout=CODE_GetSpecialValue(rcode);
    //-Assigns box.
    const unsigned box=(codeout<=CODE_OUTIGNORE?   (codeout<CODE_OUTIGNORE? cellsortfluid: BoxFluidOutIgnore):   (codetype==CODE_TYPE_FLOATING? BoxBoundOut: BoxFluidOut));
    cellpart[p]=box;
    partsincell[box]++;
  }
}

//==============================================================================
/// Calculate SortPart[] (where the particle is that must go in stated position).
/// If there are no excluded boundary particles, no problem exists.
///
/// Calcula SortPart[] (donde esta la particula que deberia ir en dicha posicion).
/// Si hay particulas de contorno excluidas no hay ningun problema.
//==============================================================================
void JCellDivCpuSingle::MakeSortFull(const unsigned* cellpart,unsigned* begincell,unsigned* partsincell,unsigned* sortpart)const{
  //-Adjust initial position of cells | Ajusta posiciones iniciales de celdas.
  begincell[0]=0;
  for(unsigned box=0;box<Nctt-1;box++)begincell[box+1]=begincell[box]+partsincell[box];
  //-Put particles in their boxes | Coloca las particulas en sus cajas.
  memset(partsincell,0,sizeof(unsigned)*(Nctt-1));
  for(unsigned p=0;p<Nptot;p++){
    unsigned box=cellpart[p];
    sortpart[begincell[box]+partsincell[box]]=p;
    partsincell[box]++;
  }
}

//==============================================================================
/// Calculate SortPart[] (where the particle is that must go in stated position).
/// In this case, there are mp excluded boundary particles because an exception is generated.
///
/// Calcula SortPart[] (donde esta la particula que deberia ir en dicha posicion).
/// En este caso nunca hay particulas bound excluidas pq se genera excepcion.
//==============================================================================
void JCellDivCpuSingle::MakeSortFluid(unsigned np,unsigned pini,const unsigned* cellpart,unsigned* begincell,unsigned* partsincell,unsigned* sortpart)const{
  //-Adjust initial position of cells | Ajusta posiciones iniciales de celdas.
  for(unsigned box=BoxFluid;box<Nctt-1;box++)begincell[box+1]=begincell[box]+partsincell[box];
  //-Put particles in their boxes | Coloca las particulas en sus cajas.
  memset(partsincell+BoxFluid,0,sizeof(unsigned)*(Nctt-1-BoxFluid));
  const unsigned pfin=pini+np;
  for(unsigned p=pini;p<pfin;p++){
    unsigned box=cellpart[p];
    sortpart[begincell[box]+partsincell[box]]=p;
    partsincell[box]++;
  }
}

//==============================================================================
/// Computes cell of each particle (CellPart[]) from dcell[], all the excluded 
/// particles have been marked  in code[].
/// Computes SortPart[] (where the particle is that must go in stated position).
///
/// Calcula celda de cada particula (CellPart[]) a partir de cell[], todas las
/// particulas excluidas ya fueron marcadas en code[].
/// Calcula SortPart[] (donde esta la particula que deberia ir en dicha posicion).
//==============================================================================
void JCellDivCpuSingle::PreSort(const unsigned* dcellc,const typecode *codec){
  //-Load SortPart[] with the current particle in the data vectors where the particle is that must go in stated position.
  //-Load BeginCell[] with first particle of each cell.
  //-Carga SortPart[] con la p actual en los vectores de datos donde esta la particula que deberia ir en dicha posicion.
  //-Carga BeginCell[] con primera particula de cada celda.
  if(DivideFull){
    PreSortFull(Nptot,dcellc,codec,CellPart,PartsInCell);
    MakeSortFull(CellPart,BeginCell,PartsInCell,SortPart);
  }
  else{
    PreSortFluid(Npf1,Npb1,dcellc,codec,CellPart,PartsInCell);
    MakeSortFluid(Npf1,Npb1,CellPart,BeginCell,PartsInCell,SortPart);
  }
  SortArray(CellPart); //-Order values of CellPart[] | Ordena valores de CellPart[].
}

//==============================================================================
/// Initial process of Divide: Calculate limits of domain and calculate new position
/// for each particle (SortPart).
/// The value np includes periodic bound & fluid particles (npbper & npfper).
/// Floating particles are treated as if they are fluid (as such they are to be excluded & effectively 
/// ignored), but in case of some excluding floating particles an exception is generated 
/// in CalcCellDomainFluid();
///
/// Inicia proceso de Divide: Calcula limites de dominio y calcula nueva posicion
/// para cada particula (SortPart).
/// El valor np incluye las periodicas bound y fluid (npbper y npfper).
/// Las floating se tratan como si fuesen fluido (tanto al ser excluidas como 
/// ignoradas), pero en caso de haber alguna floating excluida se genera una 
/// excepcion en CalcCellDomainFluid();
//==============================================================================
void JCellDivCpuSingle::Divide(unsigned npb1,unsigned npf1,unsigned npb2,unsigned npf2,bool boundchanged
  ,const unsigned *dcellc,const typecode* codec,const unsigned* idpc,const tdouble3* posc,TimersCpu timers)
{
  DivideFull=false;
  TmcStart(timers,TMC_NlLimits);

  //-Establish number of particles. | Establece numero de particulas.
  Npb1=npb1; Npf1=npf1; Npb2=npb2; Npf2=npf2;
  Nptot=Npb1+Npf1+Npb2+Npf2;
  NpbOut=NpfOut=NpbOutIgnore=NpfOutIgnore=0;
  NpFinal=NpbFinal=0;
  NpbIgnore=0;
  //:printf("---> Npb1:%u  Npf1:%u  Npb2:%u  Npf2:%u\n",Npb1,Npf1,Npb2,Npf2);

  //-Checks if the allocated memory is esufficient for Nptot.
  //-Comprueba si hay memoria reservada y si es suficiente para Nptot.
  CheckMemoryNp(Nptot);

  //-If the position of the boundary changes or there are periodic conditions it is necessary to recalculate the limits & reorder all the particles. 
  //-Si la posicion del contorno cambia o hay condiciones periodicas es necesario recalcular limites y reordenar todas las particulas. 
  if(boundchanged || PeriActive){
    BoundLimitOk=BoundDivideOk=false;
    BoundLimitCellMin=BoundLimitCellMax=TUint3(0);
    BoundDivideCellMin=BoundDivideCellMax=TUint3(0);
  }

  //-Calculate domain limits | Calcula limites del dominio.
  CalcCellDomain(dcellc,codec);
  //-Calculate number of cells for divide and check reservation of memory for cells.
  //-Calcula numero de celdas para el divide y comprueba reserva de memoria para celdas.
  PrepareNct();
  //-Check is there is memory reserved and if it is sufficient for Nptot.
  //-Comprueba si hay memoria reservada y si es suficiente para Nptot.
  CheckMemoryNct(Nct);
  TmcStop(timers,TMC_NlLimits);

  //-Determines if the divide affects all the particles.
  //-BoundDivideOk becomes false when the allocation memory changes for particles or cells.
  //-Determina si el divide afecta a todas las particulas.
  //-BoundDivideOk se vuelve false al reservar o liberar memoria para particulas o celdas.
  if(!BoundDivideOk || BoundDivideCellMin!=CellDomainMin || BoundDivideCellMax!=CellDomainMax){
    DivideFull=true;
    BoundDivideOk=true; BoundDivideCellMin=CellDomainMin; BoundDivideCellMax=CellDomainMax;
  }
  else DivideFull=false;

  //-Computes CellPart[] and SortPart[] (where the particle is that must go in stated position).
  //-Calcula CellPart[] y SortPart[] (donde esta la particula que deberia ir en dicha posicion).
  TmcStart(timers,TMC_NlMakeSort);
  PreSort(dcellc,codec);

  //-Calculate number of particles. | Calcula numeros de particulas.
  NpbIgnore=CellSize(BoxBoundIgnore);
  NpbOut=CellSize(BoxBoundOut);
  NpfOut=CellSize(BoxFluidOut);
  NpbOutIgnore=CellSize(BoxBoundOutIgnore);
  NpfOutIgnore=CellSize(BoxFluidOutIgnore);
  //:printf("---> Nct:%u  BoxBoundOut:%u  SizeBeginEndCell:%u\n",Nct,BoxBoundOut,SizeBeginEndCell(Nct));
  //:printf("---> NpbIgnore:%u  NpbOut:%u  NpfOut:%u  NpfOutIgnore:%u\n",NpbIgnore,NpbOut,NpfOut,NpfOutIgnore);
  NpFinal=Nptot-NpbOut-NpfOut-NpbOutIgnore-NpfOutIgnore;
  NpbFinal=Npb1+Npb2-NpbOutIgnore;
  if(NpbOut!=0 && DivideFull)NpbFinal=UINT_MAX; //-NpbOut can contain excluded particles fixed, moving and also floating.

  Ndiv++;
  if(DivideFull)NdivFull++;
  TmcStop(timers,TMC_NlMakeSort);
}


