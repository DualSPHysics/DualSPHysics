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
//# - Clase para uso de GPU por parte de JWaveSpectrum. (13-03-2018)
//#############################################################################

/// \file JWaveSpectrumGpu.h \brief Declares the class \ref JWaveSpectrumGpu.

#ifndef _JWaveSpectrumGPU_
#define _JWaveSpectrumGPU_

#include "TypesDef.h"
#include "JObject.h"


//##############################################################################
//# JWaveSpectrumGpu
//##############################################################################
/// \brief Manages the GPU computations for irregular wave generation (JWaveSpectrum class).

class JWaveSpectrumGpu : protected JObject
{
private:

  llong MemGpuFixed;
  tdouble4 *Order2CoefsEtag;  ///<Coefficients on GPU for elevation of each wave combination {dnm,dkl,aagnm,bbgnm} [SizeWaveCoefs].
  double *Order2CoefsDnmg;    ///<Coefficients on GPU for position of each wave combination {dnm} [SizeWaveCoefs].
  tdouble2 *Order2CoefsPosg;  ///<Coefficients on GPU for position of each wave combination {aaf1,bbf1} [SizeWaveCoefs].
  double *Order2Auxg;         ///<Auxiliary memory for reductions. [SizeWaveCoefs].

public:
  JWaveSpectrumGpu();
  ~JWaveSpectrumGpu(){ DestructorActive=true; FreeMemoryGpu(); }

  void AllocMemoryGpu(unsigned sizewavecoefs);
  void FreeMemoryGpu();
  void CopyCoefs(unsigned sizewavecoefs,const tdouble4 *d4,const double *d1,const tdouble2 *d2);

  double CalcPosition(double time,unsigned sizewavecoefs);
  double CalcElevation(double time,double x,unsigned sizewavecoefs);

  llong GetMemGpuFixed()const{ return(MemGpuFixed); }
};

#endif

