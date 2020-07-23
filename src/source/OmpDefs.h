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

/// \file OmpDefs.h \brief Defines constants to use of OpenMP.

#ifndef _OmpDefs_
#define _OmpDefs_

#include "TypesDef.h"

#define OMP_USE  ///<Enables/Disables OpenMP.
#ifdef OMP_USE
  //#define OMP_USE_RADIXSORT ///<Enables/disables OpenMP in JRadixSort.
  #define OMP_USE_WAVEGEN    ///<Enables/disables OpenMP in JWaveGen.
#endif

#ifdef OMP_USE
  #include <omp.h>  //-Active also in config. properties -> C/C++ -> Lenguage -> OpenMp.
#else
  #define omp_get_thread_num() 0
  #define omp_get_max_threads() 1
#endif

#define OMP_MAXTHREADS 64  
#define OMP_STRIDE 200
#define OMP_LIMIT_COMPUTESTEP 25000
#define OMP_LIMIT_COMPUTEMEDIUM 10000
#define OMP_LIMIT_COMPUTELIGHT 100000
#define OMP_LIMIT_PREINTERACTION 100000
#define OMP_LIMIT_TRIANGLESCELLS 3000
#define OMP_LIMIT_LIGHT 100000

#endif


