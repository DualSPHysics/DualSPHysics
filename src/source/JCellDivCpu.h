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

/// \file JCellDivCpu.h \brief Declares the class \ref JCellDivCpu.

#ifndef _JCellDivCpu_
#define _JCellDivCpu_

#include "DualSphDef.h"
#include "JObject.h"
#include "JSphTimersCpu.h"
#include "JCellDivDataCpu.h"
#include "JLog2.h"
#include <cmath>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>

//#define DBG_JCellDivCpu 1 //:DEL:

//##############################################################################
//# JCellDivCpu
//##############################################################################
/// \brief Defines the class responsible of generating the Neighbour List in CPU.

class JCellDivCpu : protected JObject
{
protected:
  const bool Stable;
  const bool Floating;
  const byte PeriActive;
  const TpCellMode CellMode;  ///<Cell division mode.
  const int ScellDiv;         ///<Value to divide KernelSize (1 or 2).
  const float Scell;          ///<Cell size: KernelSize/ScellDiv (KernelSize or KernelSize/2).
  const float OvScell;        ///<OvScell=1/Scell
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
  //-Memoria reservada en funcion de particulas en CPU.
  unsigned SizeNp;
  unsigned *CellPart;
  unsigned *SortPart;

  unsigned IncreaseNp; ///<Possible number of particles to be created in the near future.

  //-Variables with allocated memory as a function of the number of cells in GPU.
  //-Memoria reservada en funcion de celdas en GPU.
  unsigned SizeNct;
  unsigned *PartsInCell;
  unsigned *BeginCell;   ///<Get first value of each cell. | Contiene el principio de cada celda. 
  // BeginCell=[BoundOk(nct),BoundIgnore(1),Fluid(nct),BoundOut(1),FluidOut(1),BoundOutIgnore(1),FluidOutIgnore(1),END)]

  //-Variables to reorder particles. | Variables para reordenar particulas.
  byte        *VSort;            ///<Memory to reorder particles. | Memoria para reordenar particulas. [sizeof(tdouble3)*Np]
  word        *VSortWord;        ///<To order word vectors (write to VSort). | Para ordenar vectores word (apunta a VSort).
  int         *VSortInt;         ///<To order vectors int (write to VSort). | Para ordenar vectores int (apunta a VSort).
  float       *VSortFloat;       ///<To order vectors float (write to VSort). | Para ordenar vectores float (apunta a VSort).
  tfloat3     *VSortFloat3;      ///<To order vectors tfloat3 (write to VSort). | Para ordenar vectores tfloat3 (apunta a VSort).
  tfloat4     *VSortFloat4;      ///<To order vectors tfloat4 (write to VSort). | Para ordenar vectores tfloat4 (apunta a VSort).
  tdouble3    *VSortDouble3;     ///<To order vectors tdouble3 (write to VSort). | Para ordenar vectores tdouble3 (apunta a VSort).
  tsymatrix3f *VSortSymmatrix3f; ///<To order vectors tsymatrix3f (write to VSort). | Para ordenar vectores tsymatrix3f (apunta a VSort).

  llong MemAllocNp;  ///<Memory reserved for particles. | Mermoria reservada para particulas.
  llong MemAllocNct; ///<Memory reserved for cells. | Mermoria reservada para celdas.

  unsigned Ndiv,NdivFull;

  //-Number of particles by type to initialise in divide.
  //-Numero de particulas por tipo al iniciar el divide.
  unsigned Npb1;
  unsigned Npf1;
  unsigned Npb2;
  unsigned Npf2;

  unsigned Nptot;  ///<Total number of particles included that are excluded at the end of divide | Numero total de particulas incluidas las que se excluyeron al terminar el divide.
  unsigned NpbOut,NpfOut,NpbOutIgnore,NpfOutIgnore;
  
  unsigned NpFinal,NpbFinal;
  unsigned NpbIgnore;

  const bool CellDomFixed; ///<The CellDomainMin-Max is fixed according maximum domain size.
  tuint3 CellDomainMin;    ///<Lower domain limit in cells inside of DomCells. | Limite inferior del dominio en celdas dentro de DomCells.
  tuint3 CellDomainMax;    ///<Upper domain limit in cells inside of DomCells. | Limite superior del dominio en celdas dentro de DomCells.
  unsigned Ncx,Ncy,Ncz;
  unsigned Nsheet;      ///<Nsheet=Ncx*Ncy
  unsigned Nct;         ///<Nct=Ncx*Ncy*Ncz
  ullong Nctt;          ///<Total number of special cells included  Nctt=SizeBeginCell(). | Numero total de celdas incluyendo las especiales Nctt=SizeBeginCell().
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
  void FreeMemoryNp();
  void FreeMemoryAll();
  void SetMemoryVSort(byte *vsort);
  void AllocMemoryNp(ullong np);
  void AllocMemoryNct(ullong nct);
  void CheckMemoryNp(unsigned npmin);
  void CheckMemoryNct(unsigned nctmin);

  ullong SizeBeginCell(ullong nct)const{ return((nct*2)+5+1); } //-[BoundOk(nct),BoundIgnore(1),Fluid(nct),BoundOut(1),FluidOut(1),BoundOutIgnore(1),FluidOutIgnore(1),END(1)]

  ullong GetAllocMemoryNp()const{ return(MemAllocNp); };
  ullong GetAllocMemoryNct()const{ return(MemAllocNct); };
  ullong GetAllocMemory()const{ return(GetAllocMemoryNp()+GetAllocMemoryNct()); };

  //tuint3 GetMapCell(const tfloat3 &pos)const;
  void LimitsCellBound(unsigned n,unsigned pini,const unsigned* dcellc,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax)const;
  void CalcCellDomainBound(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellc,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax);
  void LimitsCellFluid(unsigned n,unsigned pini,const unsigned* dcellc,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax)const;
  void CalcCellDomainFluid(unsigned n,unsigned pini,unsigned n2,unsigned pini2,const unsigned* dcellc,const typecode *codec,tuint3 &cellmin,tuint3 &cellmax);

  unsigned CellSize(unsigned box)const{ return(BeginCell[box+1]-BeginCell[box]); }

public:
  JCellDivCpu(bool stable,bool floating,byte periactive
    ,bool celldomfixed,TpCellMode cellmode,float scell
    ,tdouble3 mapposmin,tdouble3 mapposmax,tuint3 mapcells
    ,unsigned casenbound,unsigned casenfixed,unsigned casenpb,std::string dirout
    ,bool allocfullnct=true,float overmemorynp=CELLDIV_OVERMEMORYNP,word overmemorycells=CELLDIV_OVERMEMORYCELLS);
  ~JCellDivCpu();

  void DefineDomain(unsigned cellcode,tuint3 domcelini,tuint3 domcelfin,tdouble3 domposmin,tdouble3 domposmax);

  void SortArray(word *vec);
  void SortArray(unsigned *vec);
  void SortArray(float *vec);
  void SortArray(tdouble3 *vec);
  void SortArray(tfloat3 *vec);
  void SortArray(tfloat4 *vec);
  void SortArray(tsymatrix3f *vec);

  TpCellMode GetCellMode()const{ return(CellMode); }
  int GetScellDiv()const{ return(ScellDiv); }
  float GetScell()const{ return(Scell); }

  unsigned GetNct()const{ return(Nct); }
  unsigned GetNcx()const{ return(Ncx); }
  unsigned GetNcy()const{ return(Ncy); }
  unsigned GetNcz()const{ return(Ncz); }
  tuint3 GetNcells()const{ return(TUint3(Ncx,Ncy,Ncz)); }
  unsigned GetBoxFluid()const{ return(BoxFluid); }

  tuint3 GetCellDomainMin()const{ return(CellDomainMin); }
  tuint3 GetCellDomainMax()const{ return(CellDomainMax); }
  tdouble3 GetDomainLimits(bool limitmin,unsigned slicecellmin=0)const;

  StDivDataCpu GetCellDivData()const;

  unsigned GetNpFinal()const{ return(NpFinal); }
  unsigned GetNpbFinal()const{ return(NpbFinal); }
  unsigned GetNpbIgnore()const{ return(NpbIgnore); }
  unsigned GetNpbOut()const{ return(NpbOut); }
  unsigned GetNpfOut()const{ return(NpfOut); }
  unsigned GetNpbOutIgnore()const{ return(NpbOutIgnore); }
  unsigned GetNpfOutIgnore()const{ return(NpfOutIgnore); }

  //:const unsigned* GetCellPart()const{ return(CellPart); }
  const unsigned* GetBeginCell()const{ return(BeginCell); }

  void SetIncreaseNp(unsigned increasenp){ IncreaseNp=increasenp; }

  //:bool CellNoEmpty(unsigned box,byte kind)const;
  //:unsigned CellBegin(unsigned box,byte kind)const;
  //:unsigned CellSize(unsigned box,byte kind)const;
};

#endif


