//HEAD_DSPH
/*
<DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

This file is part of DualSPHysics.

DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.

DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.
*/

/// \file JSphGpu_NN_ker.h \brief Declares functions and CUDA kernels for NN feature.

#ifndef _JSphGpu_NN_ker_
#define _JSphGpu_NN_ker_

//#include "Types.h" //depriciated
#include "DualSphDef.h"
#include <cuda_runtime_api.h>

/// Implements a set of functions and CUDA kernels for InOut feature.
namespace cusphNN {

  void CteInteractionUp_NN(unsigned phasecount,const StPhaseCte *phasecte,const StPhaseArray *phasearray);
  void PeriodicDuplicateVerlet(unsigned n,unsigned pini,tuint3 domcells,tdouble3 perinc
    ,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
    ,double2 *posxy,double *posz,float4 *velrhop,float *auxnn,float4 *velrhopm1);
  void PeriodicDuplicateSymplectic(unsigned n,unsigned pini
    ,tuint3 domcells,tdouble3 perinc,const unsigned *listp,unsigned *idp,typecode *code,unsigned *dcell
    ,double2 *posxy,double *posz,float4 *velrhop,float *auxnn,double2 *posxypre,double *poszpre,float4 *velrhoppre);

  //-Kernels for the force calculation.
  void Interaction_ForcesNN(const StInterParmsg &t);

}//end of file
#endif