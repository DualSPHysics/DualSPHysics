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

/// \file JCellDivGpuSingle.cpp \brief Implements the class \ref JCellDivGpuSingle.

#include "JCellDivGpuSingle.h"
#include "JCellDivGpuSingle_ker.h"
#include "Functions.h"
#include <climits>

using namespace std;

//==============================================================================
/// Constructor.
//==============================================================================
JCellDivGpuSingle::JCellDivGpuSingle(bool stable,bool floating,byte periactive,float kernelsize2,float poscellsize
  ,TpCellMode cellmode,float scell,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
  ,unsigned casenbound,unsigned casenfixed,unsigned casenpb,std::string dirout)
  :JCellDivGpu(stable,floating,periactive,kernelsize2,poscellsize,cellmode,scell
  ,mapposmin,mapposmax,mapcells,casenbound,casenfixed,casenpb,dirout)
{
  ClassName="JCellDivGpuSingle";
}

//==============================================================================
/// Computes cell domains adjusting to the fluid (CellDomainMin/Max). 
/// Calcula limites del dominio en celdas ajustando al fluido (CellDomainMin/Max). 
//==============================================================================
void JCellDivGpuSingle::CalcCellDomain(const unsigned *dcellg,const typecode *codeg){
  //-Calculate boundary domain | Calcula dominio del contorno.
  tuint3 celbmin,celbmax;
  if(!BoundLimitOk){
    CalcCellDomainBound(Npb1,0,Npb2,Npb1+Npf1,dcellg,codeg,celbmin,celbmax);
    BoundLimitOk=true; BoundLimitCellMin=celbmin; BoundLimitCellMax=celbmax;
  } 
  else{ celbmin=BoundLimitCellMin; celbmax=BoundLimitCellMax; }
  //Log->Printf("----->CalcCellDomain> BoundLimitCellMin/Max:%s",fun::Uint3RangeStr(BoundLimitCellMin,BoundLimitCellMax).c_str());
  //-Calculate fluid domain | Calcula dominio del fluido.
  tuint3 celfmin,celfmax;
  CalcCellDomainFluid(Npf1,Npb1,Npf2,Npb1+Npf1+Npb2,dcellg,codeg,celfmin,celfmax);
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
void JCellDivGpuSingle::MergeMapCellBoundFluid(const tuint3 &celbmin,const tuint3 &celbmax,const tuint3 &celfmin,const tuint3 &celfmax,tuint3 &celmin,tuint3 &celmax)const{
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
  if(celmin.x>celmax.x||celmin.y>celmax.y||celmin.z>celmax.z)celmin=celmax=TUint3(0,0,0);
}

//==============================================================================
/// Computes number of cells starting from (CellDomainMIN/Max).
/// Obtains location of special cells.
///
/// Calcula numero de celdas a partir de (CellDomainMin/Max). 
/// Obtiene localizacion de celdas especiales.
//==============================================================================
void JCellDivGpuSingle::PrepareNct(){
  //-Computes number of cells.
  Ncx=CellDomainMax.x-CellDomainMin.x+1;
  Ncy=CellDomainMax.y-CellDomainMin.y+1;
  Ncz=CellDomainMax.z-CellDomainMin.z+1;
  //Log->Printf("======  ncx:%u ncy:%u ncz:%u\n",Ncx,Ncy,Ncz);
  Nsheet=Ncx*Ncy; Nct=Nsheet*Ncz; Nctt=SizeBeginEndCell(Nct);
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
/// Computes cell of each particle (CellPart[]) from dcell[], all the excluded 
/// particles have been marked  in code[].
/// Assigns consecutive values to SortPart[].
///
/// Calcula celda de cada particula (CellPart[]) a partir de dcell[], todas las
/// particulas excluidas ya fueron marcadas en code[].
/// Asigna valores consecutivos a SortPart[].
//==============================================================================
void JCellDivGpuSingle::PreSort(const unsigned *dcellg,const typecode *codeg){
  if(DivideFull)cudiv::PreSortFull(Nptot,DomCellCode,dcellg,codeg,CellDomainMin,TUint3(Ncx,Ncy,Ncz),CellPart,SortPart);
  else cudiv::PreSortFluid(Npf1,Npb1,DomCellCode,dcellg,codeg,CellDomainMin,TUint3(Ncx,Ncy,Ncz),CellPart,SortPart);
}

//==============================================================================
/// Initial processing of Divide: Calculte the limits of the domain and
/// compute the new position of each particle (SortPart).
/// The value for np includes periodic boundary and fluid particles (npbper and npfper).
/// The floating bodies are treated as fluids (both to be ignored as excluded).
///
/// Inicia proceso de Divide: Calcula limites de dominio y calcula nueva posicion
/// para cada particula (SortPart).
/// El valor np incluye las periodicas bound y fluid (npbper y npfper).
/// Las floating se tratan como si fuesen fluido (tanto al ser excluidas como ignoradas).
//==============================================================================
void JCellDivGpuSingle::Divide(unsigned npb1,unsigned npf1,unsigned npb2,unsigned npf2,bool boundchanged
  ,const unsigned *dcellg,const typecode *codeg,TimersGpu timers,const double2 *posxy,const double *posz,const unsigned *idp)
{
  DivideFull=false;
  TmgStart(timers,TMG_NlLimits);

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

  //-If the boundary postion changes or there are periodic conditions it is necessary to recalculate the limits and reorder every particle.
  //-Si la posicion del contorno cambia o hay condiciones periodicas es necesario recalcular limites y reordenar todas las particulas. 
  if(boundchanged || PeriActive){
    BoundLimitOk=BoundDivideOk=false;
    BoundLimitCellMin=BoundLimitCellMax=TUint3(0);
    BoundDivideCellMin=BoundDivideCellMax=TUint3(0);
  }

  //-Calculate domain limits | Calcula limites del dominio.
  CalcCellDomain(dcellg,codeg);
  //-Calculate number of cells for divide and check reservation of memory for cells.
  //-Calcula numero de celdas para el divide y comprueba reserva de memoria para celdas.
  PrepareNct();
  //-Checks if the allocated memory is sufficient for Nptot.
  //-Comprueba si hay memoria reservada y si es suficiente para Nptot.
  CheckMemoryNct(Nct);
  TmgStop(timers,TMG_NlLimits);

  //-Determines if the divide affects all the particles.
  //-BoundDivideOk becomes false when the allocation memory changes for particles or cells.
  //-Determina si el divide afecta a todas las particulas.
  //-BoundDivideOk se vuelve false al reservar o liberar memoria para particulas o celdas.
  if(!BoundDivideOk || BoundDivideCellMin!=CellDomainMin || BoundDivideCellMax!=CellDomainMax){
    DivideFull=true;
    BoundDivideOk=true; BoundDivideCellMin=CellDomainMin; BoundDivideCellMax=CellDomainMax;
  }
  else DivideFull=false;

  //-Computes CellPart[] and assigns consecutive values to SortPart[].
  //-Calcula CellPart[] y asigna valores consecutivos a SortPart[].
  TmgStart(timers,TMG_NlPreSort);
  PreSort(dcellg,codeg);
  TmgStop(timers,TMG_NlPreSort);

  //-Sorts CellPart and SortPart as a function of the cell.
  //-Ordena CellPart y SortPart en funcion de la celda.
  TmgStart(timers,TMG_NlRadixSort);
  if(DivideFull)cudiv::Sort(CellPart,SortPart,Nptot,Stable);
  else cudiv::Sort(CellPart+Npb1,SortPart+Npb1,Nptot-Npb1,Stable);
  TmgStop(timers,TMG_NlRadixSort);

  //-Computes the initial and the last paeticle of each cell (BeginEndCell).
  //-Calcula particula inicial y final de cada celda (BeginEndCell).
  TmgStart(timers,TMG_NlCellBegin);
  cudiv::CalcBeginEndCell(DivideFull,Nptot,Npb1,unsigned(SizeBeginEndCell(Nct)),BoxFluid,CellPart,BeginEndCell);

  //-Calculate number of particles. | Calcula numeros de particulas.
  NpbIgnore=CellSize(BoxBoundIgnore);
  unsigned beginendcell[8];
  CellBeginEnd(BoxBoundOut,8,beginendcell);
  NpbOut=beginendcell[1]-beginendcell[0];
  NpfOut=beginendcell[3]-beginendcell[2];
  NpbOutIgnore=beginendcell[5]-beginendcell[4];
  NpfOutIgnore=beginendcell[7]-beginendcell[6];
  //:printf("---> Nct:%u  BoxBoundOut:%u  SizeBeginEndCell:%u\n",Nct,BoxBoundOut,SizeBeginEndCell(Nct));
  //:printf("---> NpbIgnore:%u  NpbOut:%u  NpfOut:%u  NpfOutIgnore:%u\n",NpbIgnore,NpbOut,NpfOut,NpfOutIgnore);
  NpFinal=Nptot-NpbOut-NpfOut-NpbOutIgnore-NpfOutIgnore;
  NpbFinal=Npb1+Npb2-NpbOutIgnore;
  if(NpbOut!=0 && DivideFull)NpbFinal=UINT_MAX; //-NpbOut can contain excluded particles fixed, moving and also floating.

  Ndiv++;
  if(DivideFull)NdivFull++;
  Check_CudaErroor("Error in NL construction.");
  TmgStop(timers,TMG_NlCellBegin);
}

//==============================================================================
/// Returns cell division data for neighborhood search.
/// Devuelve datos de division en celdas para busqueda de vecinos.
//==============================================================================
StDivDataGpu JCellDivGpuSingle::GetCellDivData()const{
  return(MakeDivDataGpu(ScellDiv,GetNcells(),GetCellDomainMin(),GetBeginCell()
    ,GetScell(),DomCellCode,DomPosMin,KernelSize2,PosCellSize));
}


