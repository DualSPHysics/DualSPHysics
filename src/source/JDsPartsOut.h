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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Almacena temporalmente las particulas excluidas hasta que se graben en
//:#   el siguiente part. (08-11-2011)
//:# - Las funcion GetAllocMemory() devuelve long long. (05-04-2013)
//:# - Pos pasa a ser tdouble3 en lugar de tfloat3. (24-11-2013)
//:# - Se incluye el motivo de exclusion. (20-03-2018)
//:# - Mejoras para compatibilidad con Multi-GPU. (10-09-2019)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Cambio de nombre de J.PartsOut a J.DsPartsOut. (28-06-2020)
//:#############################################################################

/// \file JDsPartsOut.h \brief Declares the class \ref JDsPartsOut.

#ifndef _JDsPartsOut_
#define _JDsPartsOut_

#include "DualSphDef.h"
#include "JObject.h"
#include <cstring>

//##############################################################################
//# JDsPartsOut
//##############################################################################
/// \brief Stores excluded particles at each instant until writing the output file. 

class JDsPartsOut : protected JObject
{
protected:
  const unsigned SizeUnit;
  unsigned Size;
  unsigned Count;
  
  unsigned OutPosCount,OutRhopCount,OutMoveCount;

  //-Normal CPU memory pointers.
  unsigned *Idp;
  tdouble3 *Pos;
  tfloat3 *Vel;
  float *Rhop;
  byte *Motive; ///<Motives for exclusion. 1:position, 2:rhop, 3:velocity.

  unsigned MemAllocs;     ///<Number of allocations.
  llong MemCpuParticles;  ///<Allocated normal CPU memory.


  void AllocMemory(unsigned size,bool reset);
  void AddData(unsigned np,const typecode* code);

public:
  JDsPartsOut(unsigned sizeunit=1024);
  ~JDsPartsOut();
  void Reset();

  unsigned GetMemAllocs()const{ return(MemAllocs); }
  llong GetAllocMemory()const{ return(MemCpuParticles); }

  void AddParticles(unsigned np,const unsigned* idp,const tdouble3* pos
    ,const tfloat3* vel,const float* rhop,const typecode* code);

  unsigned GetSize()const{ return(Size); }
  unsigned GetCount()const{ return(Count); }

  unsigned GetOutPosCount()const{ return(OutPosCount); }
  unsigned GetOutRhopCount()const{ return(OutRhopCount); }
  unsigned GetOutMoveCount()const{ return(OutMoveCount); }

  const unsigned* GetIdpOut(){ return(Idp); }
  const tdouble3* GetPosOut(){ return(Pos); }
  const tfloat3* GetVelOut(){ return(Vel); }
  const float* GetRhopOut(){ return(Rhop); }
  const byte* GetMotiveOut(){ return(Motive); }

  void Clear(){ Count=0; OutPosCount=OutRhopCount=OutMoveCount=0; };
};

#endif


