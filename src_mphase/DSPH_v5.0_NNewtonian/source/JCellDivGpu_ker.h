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

/// \file JCellDivGpu_ker.h \brief Declares functions and CUDA kernels to compute operations of the Neighbour List.

#ifndef _JCellDivGpu_ker_
#define _JCellDivGpu_ker_

#include "DualSphDef.h"
#include <cuda_runtime_api.h>

//:#define DG_LimitsCell //-En LimitsCell() comprueba que el resultado sea correcto.
//:#define DG_LimitsPos //-En LimitsPos() comprueba que el resultado sea correcto.
//:#define DG_GetRangeParticlesCells //-En GetParticlesCellRange() comprueba que el resultado sea correcto.

class JLog2;

#define DIVBSIZE 256

/// Implements a set of functions and CUDA kernels to compute operations of the Neighbour List.
namespace cudiv{

inline tfloat3 ToTFloat3(const float3& v){ return(TFloat3(v.x,v.y,v.z)); }

void Sort(unsigned* keys,unsigned* values,unsigned size,bool stable);

void ReduPosLimits(unsigned nblocks,float *aux,tfloat3 &pmin,tfloat3 &pmax);


inline unsigned LimitsPosSize(unsigned ndata){ ndata=(ndata>DIVBSIZE? ndata: DIVBSIZE); unsigned n=6,s=((ndata/DIVBSIZE)+1); return((s*n + ((s/DIVBSIZE)+1)*n) + DIVBSIZE); }


void LimitsCell(unsigned np,unsigned pini,unsigned cellcode,const unsigned *dcell
  ,const typecode *code,unsigned *aux,tuint3 &celmin,tuint3 &celmax);
void CalcBeginEndCell(bool full,unsigned np,unsigned npb,unsigned sizebegcell
  ,unsigned cellfluid,const unsigned *cellpart,int2 *begcell);

void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const unsigned *idp,const typecode *code,const unsigned *dcell,const double2 *posxy,const double *posz,const float4 *velrhop,unsigned *idp2,typecode *code2,unsigned *dcell2,double2 *posxy2,double *posz2,float4 *velrhop2);
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float4 *a,float4 *a2);
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float *a,const float *b,float *a2,float *b2);
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const double2 *a,const double *b,const float4 *c,double2 *a2,double *b2,float4 *c2);
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const tsymatrix3f *a,tsymatrix3f *a2);
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float3 *a,float3 *a2);
void SortDataParticles(unsigned np,unsigned pini,const unsigned *sortpart,const float *a,float *a2);

//:inline unsigned GetRangeParticlesCellsSizeAux(unsigned celini,unsigned celfin){ unsigned n=2,s=(((celfin-celini)/DIVBSIZE)+1); return((s*n + ((s/DIVBSIZE)+1)*n) + DIVBSIZE); } 
//:void GetRangeParticlesCells(unsigned celini,unsigned celfin,const int2 *begcell,unsigned *aux,unsigned &pmin,unsigned &pmax);

//:inline unsigned GetParticlesCellsSizeAux(unsigned celini,unsigned celfin){ unsigned n=1,s=(((celfin-celini)/DIVBSIZE)+1); return((s*n + ((s/DIVBSIZE)+1)*n) + DIVBSIZE); }  
//:unsigned GetParticlesCells(unsigned celini,unsigned celfin,const int2 *begcell,unsigned *aux);

}

#endif


