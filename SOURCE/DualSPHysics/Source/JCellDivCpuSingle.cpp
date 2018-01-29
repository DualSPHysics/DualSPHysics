//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCellDivCpuSingle::JCellDivCpuSingle(bool stable,bool floating,byte periactive,TpCellOrder cellorder
  ,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
  ,unsigned casenbound,unsigned casenfixed,unsigned casenpb,JLog2 *log,std::string dirout)
  :JCellDivCpu(stable,floating,periactive,cellorder,cellmode,scell,mapposmin,mapposmax,mapcells
  ,casenbound,casenfixed,casenpb,log,dirout)
{
  ClassName="JCellDivCpuSingle";
}

//==============================================================================
/// Calculate limits of domain in cells adjusting to fluid (CellDomainMin/Max). 
/// If there are some particles excluded of type boundary (including floating) generate
/// an exception ands show its infor.
///
/// Calcula limites del dominio en celdas ajustando al fluido (CellDomainMin/Max). 
/// Si hay alguna particula excluida de tipo boundary (incluidas floating) genera
/// excepcion y muestra su info.
//==============================================================================
void JCellDivCpuSingle::CalcCellDomain(const unsigned *dcellc,const typecode *codec,const unsigned* idpc,const tdouble3* posc){
  //-Calculate boundary domain | Calcula dominio del contorno.
  tuint3 celbmin,celbmax;
  if(!BoundLimitOk){
    CalcCellDomainBound(Npb1,0,Npb2,Npb1+Npf1,dcellc,codec,idpc,posc,celbmin,celbmax);
    BoundLimitOk=true; BoundLimitCellMin=celbmin; BoundLimitCellMax=celbmax;
  } 
  else{ celbmin=BoundLimitCellMin; celbmax=BoundLimitCellMax; }
  //-Calculate fluid domain | Calcula dominio del fluido.
  tuint3 celfmin,celfmax;
  CalcCellDomainFluid(Npf1,Npb1,Npf2,Npb1+Npf1+Npb2,dcellc,codec,idpc,posc,celfmin,celfmax);
  //-Calculate domain adjusting to boundary and fluid (with halo of 2h) | Calcula dominio ajustando al contorno y al fluido (con halo de 2h). 
  MergeMapCellBoundFluid(celbmin,celbmax,celfmin,celfmax,CellDomainMin,CellDomainMax);
}

//==============================================================================
/// Combine limits of boundary and fluid cells with limits of map.
/// With UseFluidDomain=TRUE staying with the fluid domain pluss 2h if there is a 
/// a boundary, in the contrary case, stay with domain that includes fluid and
/// boundary.
/// In case that the domain is null CellDomainMin=CellDomainMax=(0,0,0).
///
/// Combina limite de celdas de contorno y fluido con limites de mapa.
/// Con UseFluidDomain=TRUE se queda con el dominio del fluido mas 2h si hay 
/// contorno, en caso contrario se queda con el dominio que incluya fluido y
/// contorno.
/// En caso de que el dominio sea nulo CellDomainMin=CellDomainMax=(0,0,0).
//==============================================================================
void JCellDivCpuSingle::MergeMapCellBoundFluid(const tuint3 &celbmin,const tuint3 &celbmax,const tuint3 &celfmin,const tuint3 &celfmax,tuint3 &celmin,tuint3 &celmax)const{
  celmin=TUint3(max(min(celbmin.x,celfmin.x),(celfmin.x>=Hdiv? celfmin.x-Hdiv: 0)),max(min(celbmin.y,celfmin.y),(celfmin.y>=Hdiv? celfmin.y-Hdiv: 0)),max(min(celbmin.z,celfmin.z),(celfmin.z>=Hdiv? celfmin.z-Hdiv: 0)));
  celmax=TUint3(min(max(celbmax.x,celfmax.x),celfmax.x+Hdiv),min(max(celbmax.y,celfmax.y),celfmax.y+Hdiv),min(max(celbmax.z,celfmax.z),celfmax.z+Hdiv));
  if(celmax.x>=DomCells.x)celmax.x=DomCells.x-1;
  if(celmax.y>=DomCells.y)celmax.y=DomCells.y-1;
  if(celmax.z>=DomCells.z)celmax.z=DomCells.z-1;
  if(celmin.x>celmax.x||celmin.y>celmax.y||celmin.z>celmax.z){ celmin=celmax=TUint3(0,0,0); }
}

//==============================================================================
/// Calculate number of cells starting from (CellDomainMin/Max). 
/// Get location of special cells.
///
/// Calcula numero de celdas a partir de (CellDomainMin/Max). 
/// Obtiene localizacion de celdas especiales.
//==============================================================================
void JCellDivCpuSingle::PrepareNct(){
  //-Calculate number of cells | Calcula numero de celdas.
  Ncx=CellDomainMax.x-CellDomainMin.x+1;
  Ncy=CellDomainMax.y-CellDomainMin.y+1;
  Ncz=CellDomainMax.z-CellDomainMin.z+1;
  //:printf("======  ncx:%u ncy:%u ncz:%u\n",Ncx,Ncy,Ncz);
  Nsheet=Ncx*Ncy; Nct=Nsheet*Ncz; Nctt=SizeBeginCell(Nct);
  if(Nctt!=unsigned(Nctt))RunException("PrepareNct","The number of cells is too big.");
  BoxIgnore=Nct; 
  BoxFluid=BoxIgnore+1; 
  BoxBoundOut=BoxFluid+Nct; 
  BoxFluidOut=BoxBoundOut+1; 
  BoxBoundOutIgnore=BoxFluidOut+1;
  BoxFluidOutIgnore=BoxBoundOutIgnore+1;
}

//==============================================================================
/// Calculate cell of each boundary and fluid particle (cellpart[]) starting from its cell in 
/// the map,all the excluded particles are already going to be marked in code[].
/// It cannot have particles excluded bound or fluid because in such a case 
/// an exception is already generaed in LimitsCellBound/Fluid().
/// Account for particles for cell (partsincell[]).
///
/// Calcula celda de cada particula bound y fluid (cellpart[]) a partir de su celda en
/// mapa, todas las particulas excluidas ya fueron marcadas en code[].
/// No puede haber particulas bound o floating excluidas porque en tal caso ya se 
/// genera una excepcion en LimitsCellBound/Fluid().
/// Contabiliza particulas por celda (partsincell[]).
//==============================================================================
void JCellDivCpuSingle::PreSortFull(unsigned np,const unsigned *dcellc,const typecode *codec,unsigned* cellpart,unsigned* partsincell)const{
  memset(partsincell,0,sizeof(unsigned)*(Nctt-1));
  for(unsigned p=0;p<np;p++){
    unsigned rcell=dcellc[p];
    unsigned cx=PC__Cellx(DomCellCode,rcell)-CellDomainMin.x;
    unsigned cy=PC__Celly(DomCellCode,rcell)-CellDomainMin.y;
    unsigned cz=PC__Cellz(DomCellCode,rcell)-CellDomainMin.z;
    const unsigned cellsort=cx+cy*Ncx+cz*Nsheet;
    const typecode rcode=codec[p];
    const bool xbound=(CODE_GetType(rcode)<CODE_TYPE_FLOATING);
    const typecode codeout=CODE_GetSpecialValue(rcode);
    unsigned box;
    if(xbound){//-Bound particles (except floating) | Particulas bound (excepto floating).
      box=(codeout<CODE_OUTIGNORE? ((cx<Ncx && cy<Ncy && cz<Ncz)? cellsort: BoxIgnore): (codeout==CODE_OUTIGNORE? BoxBoundOutIgnore: BoxBoundOut));
    }
    else{//-Fluid particles | Particulas fluid.
      box=(codeout<CODE_OUTIGNORE? BoxFluid+cellsort: (codeout==CODE_OUTIGNORE? BoxFluidOutIgnore: BoxFluidOut));
    }
    cellpart[p]=box;
    partsincell[box]++;
  }
}

//==============================================================================
/// Calculate cell of each boundary and fluid particle (cellpart[]) starting from its cell in 
/// the map, all the excluded particles are already going to be marked in code[].
/// It cannot have particles excluded bound or floating because in such a case 
/// an exception is already generaed in LimitsCellBound/Fluid().
/// Account for particles for cell (partsincell[]).
///
/// Calcula celda de cada particula bound y fluid (cellpart[]) a partir de su celda en
/// mapa, todas las particulas excluidas ya fueron marcadas en code[].
/// No puede haber particulas bound o floating excluidas porque en tal caso ya se 
/// genera una excepcion en LimitsCellBound/Fluid().
/// Contabiliza particulas por celda (partsincell[]).
//==============================================================================
void JCellDivCpuSingle::PreSortFluid(unsigned np,unsigned pini,const unsigned *dcellc,const typecode *codec,unsigned* cellpart,unsigned* partsincell)const{
  memset(partsincell+BoxFluid,0,sizeof(unsigned)*(Nctt-1-BoxFluid));
  const unsigned pfin=pini+np;
  for(unsigned p=pini;p<pfin;p++){
    unsigned rcell=dcellc[p];
    unsigned cx=PC__Cellx(DomCellCode,rcell)-CellDomainMin.x;
    unsigned cy=PC__Celly(DomCellCode,rcell)-CellDomainMin.y;
    unsigned cz=PC__Cellz(DomCellCode,rcell)-CellDomainMin.z;
    const unsigned cellsort=BoxFluid+cx+cy*Ncx+cz*Nsheet;
    const typecode codeout=CODE_GetSpecialValue(codec[p]);
    unsigned box=(codeout<CODE_OUTIGNORE? cellsort: (codeout==CODE_OUTIGNORE? BoxFluidOutIgnore: BoxFluidOut));
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
/// Calculate cell of each particle (CellPart[]) starting from cell[], all the
/// excluded particles will already be marked in code[].
/// Calculate SortPart[] (where the particle is that must go in stated position).
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
  const char met[]="Divide";
  DivideFull=false;
  TmcStart(timers,TMC_NlLimits);

  //-Establish number of particles. | Establece numero de particulas.
  Npb1=npb1; Npf1=npf1; Npb2=npb2; Npf2=npf2;
  Nptot=Npb1+Npf1+Npb2+Npf2;
  NpbOut=NpfOut=NpbOutIgnore=NpfOutIgnore=0;
  NpFinal=NpbFinal=0;
  NpfOutRhop=NpfOutMove=NpbIgnore=0;
  //:printf("---> Npb1:%u  Npf1:%u  Npb2:%u  Npf2:%u\n",Npb1,Npf1,Npb2,Npf2);

  //-Check if there is memory reserved and if it sufficient for Nptot | Comprueba si hay memoria reservada y si es suficiente para Nptot.
  CheckMemoryNp(Nptot);

  //-If the position of the boundary changes or there are periodic conditions it is necessary to recalculate the limits & reorder all the particles. 
  //-Si la posicion del contorno cambia o hay condiciones periodicas es necesario recalcular limites y reordenar todas las particulas. 
  if(boundchanged || PeriActive){
    BoundLimitOk=BoundDivideOk=false;
    BoundLimitCellMin=BoundLimitCellMax=TUint3(0);
    BoundDivideCellMin=BoundDivideCellMax=TUint3(0);
  }

  //-Calculate domain limits | Calcula limites del dominio.
  CalcCellDomain(dcellc,codec,idpc,posc);
  //-Calculate number of cells for divide and check reservation of memory for cells.
  //-Calcula numero de celdas para el divide y comprueba reserva de memoria para celdas.
  PrepareNct();
  //-Check is there is memory reserved and if it is sufficient for Nptot.
  //-Comprueba si hay memoria reservada y si es suficiente para Nptot.
  CheckMemoryNct(Nct);
  TmcStop(timers,TMC_NlLimits);

  //-Determine if divide affects all the particles.
  //-BoundDivideOk returns false in order to reserve or free memory for particles or cells.
  //-Determina si el divide afecta a todas las particulas.
  //-BoundDivideOk se vuelve false al reservar o liberar memoria para particulas o celdas.
  if(!BoundDivideOk || BoundDivideCellMin!=CellDomainMin || BoundDivideCellMax!=CellDomainMax){
    DivideFull=true;
    BoundDivideOk=true; BoundDivideCellMin=CellDomainMin; BoundDivideCellMax=CellDomainMax;
  }
  else DivideFull=false;

  //- Calculate CellPart[] and SortPart[] ((where the particle is that must go in stated position)).
  //- Calcula CellPart[] y SortPart[] (donde esta la particula que deberia ir en dicha posicion).
  TmcStart(timers,TMC_NlMakeSort);
  PreSort(dcellc,codec);

  //-Calculate number of particles. | Calcula numeros de particulas.
  NpbIgnore=CellSize(BoxIgnore);
  NpbOut=CellSize(BoxBoundOut);
  if(NpbOut)RunException(met,"There cannot be excluded boundary particles.");
  NpfOut=CellSize(BoxFluidOut);
  NpbOutIgnore=CellSize(BoxBoundOutIgnore);
  NpfOutIgnore=CellSize(BoxFluidOutIgnore);
  //:printf("---> Nct:%u  BoxBoundOut:%u  SizeBeginEndCell:%u\n",Nct,BoxBoundOut,SizeBeginEndCell(Nct));
  //:printf("---> NpbIgnore:%u  NpbOut:%u  NpfOut:%u  NpfOutIgnore:%u\n",NpbIgnore,NpbOut,NpfOut,NpfOutIgnore);
  NpFinal=Nptot-NpbOut-NpfOut-NpbOutIgnore-NpfOutIgnore;
  NpbFinal=Npb1+Npb2-NpbOut-NpbOutIgnore;

  Ndiv++;
  if(DivideFull)NdivFull++;
  TmcStop(timers,TMC_NlMakeSort);
}


