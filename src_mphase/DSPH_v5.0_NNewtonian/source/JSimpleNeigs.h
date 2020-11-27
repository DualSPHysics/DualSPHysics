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

/// \file JSimpleNeigs.h \brief Declares the class \ref JSimpleNeigs.

#ifndef _JSimpleNeigs_
#define _JSimpleNeigs_

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para buscar posiciones cercanas de forma simple con un rendimiento
//:#   aceptable. (15-09-2018)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Nuevo metodo NearbyPositionsLt(). (30-10-2020)
//:#############################################################################

#include "JObject.h"
#include "TypesDef.h"
#include <climits>
#include <vector>

//##############################################################################
//# JSimpleNeigs
//##############################################################################
/// \brief Implements a class to search for nearby positions.

class JSimpleNeigs : protected JObject
{
protected:
  const unsigned Np;
  const tdouble3 *Pos;
  const double Scell;

  tdouble3 PosMin;
  tdouble3 PosMax;
  int Ncx,Ncy,Ncz,Nsheet,Nct;

  unsigned *PosInCell; //-Positions in cells [Np].
  unsigned *BeginCell; //-First positions in each cell [Nct+1].

  unsigned CountSelect;
  unsigned SizeSelect;
  unsigned *SelectPos; //-Selected positions [SizeSelect].

  void DefineMapCells();
  void CreateMapCells();

  tint3 GetCell3(const tdouble3 &ps)const{ return(TInt3(int((ps.x-PosMin.x)/Scell),int((ps.y-PosMin.y)/Scell),int((ps.z-PosMin.z)/Scell))); }
  unsigned GetCell(const tint3 &cel)const{ return(cel.x<0 || cel.y<0 || cel.z<0 || cel.x>=Ncx || cel.y>=Ncy || cel.z>=Ncz? UINT_MAX: unsigned(cel.x+cel.y*Ncx+cel.z*Nsheet)); }
  void GetNearbyCells(const tdouble3 &ps,double dist,tint3 &celmin,tint3 &celmax)const;
  void SelectAdd(unsigned p);

public:
  JSimpleNeigs(unsigned np,const tdouble3* pos,double scell);
  ~JSimpleNeigs();
  void Reset();
  unsigned GetAllocMemory()const;

  unsigned NearbyPositions(const tdouble3 &ps,unsigned pignore,double dist);
  unsigned NearbyPositionsLt(const tdouble3 &ps,unsigned pignore,double dist,std::vector<unsigned> &vsel)const;

  unsigned GetCountSelect()const{ return(CountSelect); }
  const unsigned* GetSelectPos()const{ return(SelectPos); }

};

#endif


