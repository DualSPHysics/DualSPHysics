//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2021 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDsDcell.cpp \brief Implements the class \ref JDsDcell.

#include "JDsDcell.h"
#include "JDsDcellDef.h"
#include "Functions.h"

using namespace std;

//==============================================================================
/// Calculates number of bits to store the indicated value.
//==============================================================================
unsigned JDsDcell::CalcBitsValue(unsigned v,unsigned minbits){
  minbits=(minbits<1||minbits>30? 1: minbits);
  unsigned n=minbits; for(;v>>n;n++);
  return(n);
}

//==============================================================================
/// Calculates the distribution of bits between the X, Y and Z values.
/// When maxbits==0 returns the minimum number of bits.
/// When maxbits!=0 returns 0 if maxbits is not enough.
//==============================================================================
tuint3 JDsDcell::CalcCellDistribution(const tuint3 &ncells,const unsigned maxbits){
  unsigned sxmin=CalcBitsValue(ncells.x,2);
  unsigned symin=CalcBitsValue(ncells.y,2);
  unsigned szmin=CalcBitsValue(ncells.z,2);
  unsigned smin=sxmin+symin+szmin;
  if(smin<=maxbits){
    unsigned sx=sxmin,sy=symin,sz=szmin;
    unsigned rest=maxbits-smin;
    while(rest){
      if(rest){ sx++; rest--; }
      if(rest){ sy++; rest--; }
      if(rest){ sz++; rest--; }
    }
    sxmin=sx; symin=sy; szmin=sz;
  }
  return(maxbits && sxmin+symin+szmin>maxbits? TUint3(0): TUint3(sxmin,symin,szmin));
}

//==============================================================================
/// Selects an adequate code for cell configuration or 0 for invalid cell number.
//==============================================================================
unsigned JDsDcell::CalcCellCode(const tuint3 &ncells){
  const tuint3 nbits=CalcCellDistribution(ncells,31);
  unsigned dcelcode=(nbits==TUint3(0)? 0: DCEL_GetCode(nbits.x,nbits.y,nbits.z));
  //printf("==>CalcCellCode> ncells:(%u,%u,%u) sb:(%u,%u,%u) dcc:%X\n",ncells.x,ncells.y,ncells.z,nbits.x,nbits.y,nbits.z,dcelcode);
  return(dcelcode);
}

//==============================================================================
/// Returns true when cellcode is invalid for the number of cells.
//==============================================================================
bool JDsDcell::InvalidCellCode(unsigned dcc,const tuint3 &ncells){
  return(dcc==0
      || DCEL_MaxCellx(dcc)<ncells.x 
      || DCEL_MaxCelly(dcc)<ncells.y 
      || DCEL_MaxCellz(dcc)<ncells.z);
}

//==============================================================================
/// Returns dcellcode as string with bits distribution.
//==============================================================================
std::string JDsDcell::DcellCodeStr(unsigned dcc){
  return(fun::PrintStr("1+%u_%u_%u",DCEL_GetSx(dcc)-1,DCEL_GetSy(dcc),DCEL_GetSz(dcc)));
}

//==============================================================================
/// Prints information on dcell code (for debug).
//==============================================================================
void JDsDcell::PrintDcellCodeInfo(unsigned dcc){
  const unsigned sx=DCEL_GetSx(dcc);
  const unsigned sy=DCEL_GetSy(dcc);
  const unsigned sz=DCEL_GetSz(dcc);
  const unsigned ss=sx+sy+sz;
  const unsigned vx=DCEL_MaxCellx(dcc);
  const unsigned vy=DCEL_MaxCelly(dcc);
  const unsigned vz=DCEL_MaxCellz(dcc);
  printf("[%X]-> bits(%u,%u,%u):%u  maxvalue:(%u,%u,%u)\n",dcc,sx,sy,sz,ss,vx,vy,vz);
}
