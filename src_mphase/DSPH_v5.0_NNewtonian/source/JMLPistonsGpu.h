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

//#############################################################################
//# Cambios:
//# =========
//# - Clase para uso de GPU por parte de JMLPistons. (24-05-2018)
//#############################################################################

/// \file JMLPistonsGpu.h \brief Declares the class \ref JMLPistonsGpu.

#ifndef _JMLPistonsGpu_
#define _JMLPistonsGpu_

#include "TypesDef.h"
#include "JObject.h"


//##############################################################################
//# JMLPistonsGpu
//##############################################################################
/// \brief Manages the GPU computations for Multi-Layer Pistons (JMLPistons class).

class JMLPistonsGpu : protected JObject
{
private:
  byte *PistonIdg;     ///<Data of PistonId on GPU.
  double *MovVelg;     ///<Data of MovVel on GPU.
  llong MemGpuFixed;  

public:
  JMLPistonsGpu();
  ~JMLPistonsGpu();
  void FreeMemoryGpu();
  inline llong GetAllocMemoryGpu()const{ return(MemGpuFixed); }

  void PreparePiston1d(unsigned sizepistonid,const byte *pistonid,unsigned sizemovvel);
  void CopyMovVel(unsigned sizemovvel,const double *movvel);

  inline const byte* GetPistonIdg()const{ return(PistonIdg);}
  inline const double* GetMovVelg()const{ return(MovVelg);}
};

//##############################################################################
//# JMLPiston2DGpu
//##############################################################################
/// \brief Manages the GPU computations for Multi-Layer Pistons (JMLPiston2D class).

class JMLPiston2DGpu : protected JObject
{
private:
  unsigned Size;
  double *MovVelyzg;  ///<Data of MovVelyz on GPU. [Size]

public:
  JMLPiston2DGpu();
  ~JMLPiston2DGpu();
  void FreeMemoryGpu();
  inline llong GetAllocMemoryGpu()const{ return(sizeof(double)*Size); }

  void AllocMemoryGpu(unsigned size);
  void CopyMovVelyz(const double *movvelyz);

  inline const double* GetMovVelyzg()const{ return(MovVelyzg); }

};

#endif

