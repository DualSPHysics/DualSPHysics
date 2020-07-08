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

/// \file JCellDivGpuSingle_ker.cu \brief Implements functions and CUDA kernels to compute operations of the Neighbour List.

#include "JCellDivGpuSingle_ker.h"
#include "DualSphDef.h"
#include <float.h>
#include "JLog2.h"

namespace cudiv{
#include "FunctionsBasic_iker.h"

//------------------------------------------------------------------------------
/// Processes bound and fluid particles that may be mixed.
/// Computes cell of each boundary and fluid particle (cellpart[]) starting from its cell in 
/// the map. all the excluded particles were already marked in code[].
/// Excluded particles bound (fixed and moving) and floating are moved to BoxBoundOut.
/// Assigns consecutive values to SortPart[].
///
/// Calcula celda de cada particula bound y fluid (cellpart[]) a partir de su celda en
/// mapa. Todas las particulas excluidas ya fueron marcadas en code[].
/// Las particulas excluidas de tipo bound (fixed and moving) and floating se mueven a BoxBoundOut.
/// Asigna valores consecutivos a SortPart[].
//------------------------------------------------------------------------------
__global__ void KerPreSortFull(unsigned np,unsigned cellcode,const unsigned *dcell
  ,const typecode *code,uint3 cellzero,uint3 ncells,unsigned *cellpart,unsigned *sortpart)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<np){
    sortpart[p]=p;
    const unsigned nsheet=ncells.x*ncells.y;
    const unsigned cellboundignore=nsheet*ncells.z;          //-cellboundignore==nct
    const unsigned cellfluid=cellboundignore+1;
    const unsigned cellboundout=cellfluid+cellboundignore;   //-For fixed, moving and floating.
    const unsigned cellfluidout=cellboundout+1;              //-For fluid.
    const unsigned cellboundoutignore=cellfluidout+1;        //-For bound.
    const unsigned cellfluidoutignore=cellboundoutignore+1;  //-For fluid and floatings.
    //-Computes cell according position.
    const unsigned rcell=dcell[p];
    const unsigned cx=PC__Cellx(cellcode,rcell)-cellzero.x;
    const unsigned cy=PC__Celly(cellcode,rcell)-cellzero.y;
    const unsigned cz=PC__Cellz(cellcode,rcell)-cellzero.z;
    const unsigned cellsort=cx+cy*ncells.x+cz*nsheet;
    //-Checks particle code.
    const typecode rcode=code[p];
    const typecode codetype=CODE_GetType(rcode);
    const typecode codeout=CODE_GetSpecialValue(rcode);
    //-Assigns box.
    if(codetype<CODE_TYPE_FLOATING){//-Bound particles (except floating) | Particulas bound (excepto floating).
      cellpart[p]=(codeout<CODE_OUTIGNORE?   ((cx<ncells.x && cy<ncells.y && cz<ncells.z)? cellsort: cellboundignore):   (codeout==CODE_OUTIGNORE? cellboundoutignore: cellboundout));
    }
    else{//-Fluid and floating particles | Particulas fluid y floating.
      cellpart[p]=(codeout<=CODE_OUTIGNORE?   (codeout<CODE_OUTIGNORE? cellfluid+cellsort: cellfluidoutignore):   (codetype==CODE_TYPE_FLOATING? cellboundout: cellfluidout));
    }
  }
}

//==============================================================================
/// Processes bound and fluid particles that may be mixed.
/// Computes cell of each boundary and fluid particle (cellpart[]) starting from its cell in 
/// the map. all the excluded particles were already marked in code[].
/// Excluded particles bound (fixed and moving) and floating are moved to BoxBoundOut.
/// Assigns consecutive values to SortPart[].
///
/// Calcula celda de cada particula bound y fluid (cellpart[]) a partir de su celda en
/// mapa. Todas las particulas excluidas ya fueron marcadas en code[].
/// Las particulas excluidas de tipo bound (fixed and moving) and floating se mueven a BoxBoundOut.
/// Asigna valores consecutivos a SortPart[].
//==============================================================================
void PreSortFull(unsigned np,unsigned cellcode,const unsigned *dcell,const typecode *code
  ,tuint3 cellmin,tuint3 ncells,unsigned *cellpart,unsigned *sortpart)
{
  if(np){
    const uint3 cellzero=make_uint3(cellmin.x,cellmin.y,cellmin.z);
    const uint3 xncells=make_uint3(ncells.x,ncells.y,ncells.z);
    dim3 sgrid=GetSimpleGridSize(np,DIVBSIZE);
    KerPreSortFull <<<sgrid,DIVBSIZE>>> (np,cellcode,dcell,code,cellzero,xncells,cellpart,sortpart);
  }
}

//==============================================================================
/// Processes only fluid particles.
/// Computes cell of each fluid particle (cellpart[]) starting from its cell in 
/// the map. all the excluded particles were already marked in code[].
/// Excluded particles floating are moved to BoxBoundOut.
/// Assigns consecutive values to SortPart[].
///
/// Calcula celda de cada particula fluid (cellpart[]) a partir de su celda en
/// mapa. Todas las particulas excluidas ya fueron marcadas en code[].
/// Las particulas excluidas de tipo floating se mueven a BoxBoundOut.
/// Asigna valores consecutivos a SortPart[].
//==============================================================================
__global__ void KerPreSortFluid(unsigned n,unsigned pini,unsigned cellcode
  ,const unsigned *dcell,const typecode *code,uint3 cellzero,uint3 ncells
  ,unsigned *cellpart,unsigned *sortpart)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Particle number.
  if(p<n){
    p+=pini;
    sortpart[p]=p;
    const unsigned nsheet=ncells.x*ncells.y;
    const unsigned cellfluid=nsheet*ncells.z+1;
    const unsigned cellfluidout=cellfluid+cellfluid;   //-For fluid.
    const unsigned cellboundout=cellfluidout-1;        //-For floatings.
    const unsigned cellfluidoutignore=cellfluidout+2;  //-For fluid and floatings.
    //-Computes cell according position.
    const unsigned rcell=dcell[p];
    const unsigned cx=PC__Cellx(cellcode,rcell)-cellzero.x;
    const unsigned cy=PC__Celly(cellcode,rcell)-cellzero.y;
    const unsigned cz=PC__Cellz(cellcode,rcell)-cellzero.z;
    const unsigned cellsortfluid=cellfluid+cx+cy*ncells.x+cz*nsheet;
    //-Checks particle code.
    const typecode rcode=code[p];
    const typecode codetype=CODE_GetType(rcode);
    const typecode codeout=CODE_GetSpecialValue(rcode);
    //-Assigns box.
    cellpart[p]=(codeout<=CODE_OUTIGNORE?   (codeout<CODE_OUTIGNORE? cellsortfluid: cellfluidoutignore):   (codetype==CODE_TYPE_FLOATING? cellboundout: cellfluidout));
  }
}

//==============================================================================
/// Processes only fluid particles.
/// Computes cell of each fluid particle (cellpart[]) starting from its cell in 
/// the map. all the excluded particles were already marked in code[].
/// Excluded particles floating are moved to BoxBoundOut.
/// Assigns consecutive values to SortPart[].
///
/// Calcula celda de cada particula fluid (cellpart[]) a partir de su celda en
/// mapa. Todas las particulas excluidas ya fueron marcadas en code[].
/// Las particulas excluidas de tipo floating se mueven a BoxBoundOut.
/// Asigna valores consecutivos a SortPart[].
//==============================================================================
void PreSortFluid(unsigned npf,unsigned pini,unsigned cellcode,const unsigned *dcell,const typecode *code
  ,tuint3 cellmin,tuint3 ncells,unsigned *cellpart,unsigned *sortpart)
{
  if(npf){
    const uint3 cellzero=make_uint3(cellmin.x,cellmin.y,cellmin.z);
    const uint3 xncells=make_uint3(ncells.x,ncells.y,ncells.z);
    dim3 sgrid=GetSimpleGridSize(npf,DIVBSIZE);
    KerPreSortFluid <<<sgrid,DIVBSIZE>>> (npf,pini,cellcode,dcell,code,cellzero,xncells,cellpart,sortpart);
  }
}


}


