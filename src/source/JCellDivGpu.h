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

/// \file JCellDivGpu.h \brief Declares the class \ref JCellDivGpu.

#ifndef _JCellDivGpu_
#define _JCellDivGpu_

#include "DualSphDef.h"
#include "JObjectGpu.h"
#include "JSphTimersGpu.h"
#include "JLog2.h"
#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>

//##############################################################################
//# JCellDivGpu
//##############################################################################
/// \brief Defines the class responsible of generating the Neighbour List in GPU.

class JCellDivGpu : protected JObjectGpu
{
protected:
  const bool Stable;
  const bool Floating;
  const byte PeriActive;

  const TpCellMode CellMode;  ///<Cell division mode.
  const int ScellDiv;         ///<Value to divide KernelSize (1 or 2).
  const float Scell;          ///<Cell size: KernelSize/ScellDiv (KernelSize or KernelSize/2).
  const float OvScell;        ///<OvScell=1/Scell
  const float KernelSize2;    ///<Maximum interaction distance squared (KernelSize^2).
  const float PosCellSize;    ///<Size of cells used for coding PosCell (it is usually KernelSize).
  const tdouble3 Map_PosMin,Map_PosMax,Map_PosDif;
  const tuint3 Map_Cells;
  const unsigned CaseNbound,CaseNfixed,CaseNpb;
  JLog2 *Log;
  std::string DirOut;

  bool AllocFullNct;     ///<Resserve memory for max number of cells of domain (DomCells). | Reserva memoria para el numero maximo de celdas del dominio (DomCells).
  float OverMemoryNp;    ///<Percentage that is added to the memory reserved for Np. (def=0) | Porcentaje que se anhade a la reserva de memoria de Np. (def=0).
  word OverMemoryCells;  ///<Cell number that is incremented in each dimension to reserve memory. | Numero celdas que se incrementa en cada dimension reservar memoria. (def=0).

  //-Variables to define the domain.
  unsigned DomCellCode;  ///<Key for codifying cell of position. | Clave para la codificacion de la celda de posicion.
  tuint3 DomCelIni,DomCelFin;
  tdouble3 DomPosMin,DomPosMax;
  tuint3 DomCells;

  //-Variables with allocated memory as a function of the number of particles in CPU.
  //-Memoria reservada en funcion de particulas en GPU.
  unsigned SizeNp;
  unsigned *CellPart;
  unsigned *SortPart;
  unsigned SizeAuxMem;
  float *AuxMem;

  unsigned IncreaseNp; ///<Possible number of particles to be created in the near future.

  //-Variables with allocated memory as a function of the number of cells in GPU.
  //-Memoria reservada en funcion de celdas en GPU.
  unsigned SizeNct;
  int2 *BeginEndCell;  ///<Contains the first and final particle of each cell. | Contiene el principio y final de cada celda. 
  // BeginEndCell=[BoundOk(nct),BoundIgnore(1),Fluid(nct),BoundOut(1),FluidOut(1),BoundOutIgnore(1),FluidOutIgnore(1)]

  ullong MemAllocGpuNp;  ///<GPU memory reserved for particles. | Mermoria GPU reservada para particulas.
  ullong MemAllocGpuNct; ///<GPU memory reserved for cells. | Mermoria GPU reservada para celdas.

  unsigned Ndiv,NdivFull;

  //-Number of particles by type to initialise in divide.
  //-Numero de particulas por tipo al iniciar el divide.
  unsigned Npb1;
  unsigned Npf1;
  unsigned Npb2;
  unsigned Npf2;

  unsigned Nptot;  ///<Total number of particles included that are excluded at the end of divide. | Numero total de particulas incluidas las que se excluyeron al terminar el divide.
  unsigned NpbOut,NpfOut,NpbOutIgnore,NpfOutIgnore;
  
  unsigned NpFinal,NpbFinal;
  unsigned NpbIgnore;

  const bool CellDomFixed; ///<The CellDomainMin-Max is fixed according maximum domain size.
  tuint3 CellDomainMin;    ///<Lower domain limit in cells inside of DomCells. | Limite inferior del dominio en celdas dentro de DomCells.
  tuint3 CellDomainMax;    ///<Upper domain limit in cells inside of DomCells. | Limite superior del dominio en celdas dentro de DomCells.
  unsigned Ncx,Ncy,Ncz,Nsheet,Nct;
  ullong Nctt;          ///<Total number of special cells included  Nctt=SizeBeginEndCell(). | Numero total de celdas incluyendo las especiales Nctt=SizeBeginEndCell().
  unsigned BoxBoundIgnore,BoxFluid,BoxBoundOut,BoxFluidOut,BoxBoundOutIgnore,BoxFluidOutIgnore;

  bool BoundLimitOk;    ///<Indicate that the boundary limits are already calculated in BoundLimitCellMin & BoundLimitCellMax. | Indica que los limites del contorno ya estan calculados en BoundLimitCellMin y BoundLimitCellMax.
  tuint3 BoundLimitCellMin,BoundLimitCellMax;

  bool BoundDivideOk;   ///<Indicate that the limits of boundaries used in  previous divide will go in BoundDivideCellMin & BoundDivideCellMax. | Indica que los limites del contorno utilizados en el divide previo fueron BoundDivideCellMin y BoundDivideCellMax.
  tuint3 BoundDivideCellMin,BoundDivideCellMax;

  bool DivideFull;      ///<Indicate that divie is applied to fluid & boundary (not only to fluid). | Indica que el divide se aplico a fluido y contorno (no solo al fluido).

  void Reset();

  //-Management of allocated dynamic memory.
  //-Gestion de reserva dinamica de memoria.
  void FreeMemoryNct();
  void FreeMemoryAll();
  void AllocMemoryNp(ullong np);
  void AllocMemoryNct(ullong nct);
  void CheckMemoryNp(unsigned npmin);
  void CheckMemoryNct(unsigned nctmin);

  ullong SizeBeginEndCell(ullong nct)const{ return((nct*2)+5); } //-[BoundOk(nct),BoundIgnore(1),Fluid(nct),BoundOut(1),FluidOut(1),BoundOutIgnore(1),FluidOutIgnore(1)]

  ullong GetAllocMemoryCpu()const{ return(0); }
  ullong GetAllocMemoryGpuNp()const{ return(MemAllocGpuNp); };
  ullong GetAllocMemoryGpuNct()const{ return(MemAllocGpuNct); };
  ullong GetAllocMemoryGpu()const{ return(GetAllocMemoryGpuNp()+GetAllocMemoryGpuNct()); };

  //:tuint3 GetMapCell(const tfloat3 &pos)const;
  void CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellg,const typecode *codeg,tuint3 &cellmin,tuint3 &cellmax);
  void CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellg,const typecode *codeg,tuint3 &cellmin,tuint3 &cellmax);

  void CellBeginEnd(unsigned cell,unsigned ndata,unsigned* data)const;
  int2 CellBeginEnd(unsigned cell)const;
  unsigned CellSize(unsigned cell)const{ int2 v=CellBeginEnd(cell); return(unsigned(v.y-v.x)); }

public:
  JCellDivGpu(bool stable,bool floating,byte periactive
    ,float kernelsize2,float poscellsize
    ,bool celldomfixed,TpCellMode cellmode,float scell
    ,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
    ,unsigned casenbound,unsigned casenfixed,unsigned casenpb,std::string dirout
    ,bool allocfullnct=true,float overmemorynp=CELLDIV_OVERMEMORYNP,word overmemorycells=CELLDIV_OVERMEMORYCELLS);
  ~JCellDivGpu();
  void FreeMemoryGpu(){ FreeMemoryAll(); }

  void DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin,tdouble3 domposmin,tdouble3 domposmax);

  void SortBasicArrays(const unsigned *idp,const typecode *code,const unsigned *dcell,const double2 *posxy,const double *posz,const float4 *velrhop,unsigned *idp2,typecode *code2,unsigned *dcell2,double2 *posxy2,double *posz2,float4 *velrhop2);
  void SortDataArrays(const float4 *a,float4 *a2);
  void SortDataArrays(const float *a,const float *b,float *a2,float *b2);
  void SortDataArrays(const double2 *a,const double *b,const float4 *c,double2 *a2,double *b2,float4 *c2);
  void SortDataArrays(const tsymatrix3f *a,tsymatrix3f *a2);
  void SortDataArrays(const float3 *a,float3 *a2);
  void SortDataArrays(const float *a,float *a2);

  float* GetAuxMem(unsigned size);

  TpCellMode GetCellMode()const{ return(CellMode); }
  int GetScellDiv()const{ return(ScellDiv); }
  float GetScell()const{ return(Scell); }
//:  tuint3 GetDomCells()const{ return(DomCells); };
//:  unsigned GetCellCode()const{ return(DomCellCode); };

  unsigned GetNct()const{ return(Nct); }
  unsigned GetNcx()const{ return(Ncx); }
  unsigned GetNcy()const{ return(Ncy); }
  unsigned GetNcz()const{ return(Ncz); }
  tuint3 GetNcells()const{ return(TUint3(Ncx,Ncy,Ncz)); }
  unsigned GetBoxFluid()const{ return(BoxFluid); }

  tuint3 GetCellDomainMin()const{ return(CellDomainMin); }
  tuint3 GetCellDomainMax()const{ return(CellDomainMax); }
  tdouble3 GetDomainLimits(bool limitmin,unsigned slicecellmin=0)const;

  unsigned GetNpFinal()const{ return(NpFinal); }
  unsigned GetNpbFinal()const{ return(NpbFinal); }
  unsigned GetNpbIgnore()const{ return(NpbIgnore); }
  unsigned GetNpbOut()const{ return(NpbOut); }
  unsigned GetNpfOut()const{ return(NpfOut); }
  unsigned GetNpbOutIgnore()const{ return(NpbOutIgnore); }
  unsigned GetNpfOutIgnore()const{ return(NpfOutIgnore); }

  //:const unsigned* GetCellPart()const{ return(CellPart); }
  const int2* GetBeginCell()const{ return(BeginEndCell); }

  void SetIncreaseNp(unsigned increasenp){ IncreaseNp=increasenp; }

  //:uint2 GetRangeParticlesCells(bool fluid,unsigned celini,unsigned celfin)const;
  //:unsigned GetParticlesCells(unsigned celini,unsigned celfin);


  //:tuint3 DgGetCell(const tfloat3 &pos)const;
  //:void DgSaveVktIdx(std::string file,unsigned np,const unsigned *idx,const unsigned *idp,const tfloat3 *pos)const;
  //:void DgSaveVktRange(std::string file,unsigned pini,unsigned pfin,const unsigned *idpg,const float3 *posg)const;
};

#endif


