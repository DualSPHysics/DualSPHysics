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
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

/// \file JWaveOrder2_ker.h \brief Declares functions and CUDA kernels for second order irregular wave generation.

#ifndef _JWaveOrder2_ker_
#define _JWaveOrder2_ker_

#include <cuda_runtime_api.h>

#define WAVEBSIZE 256

/// Implements a set of functions and CUDA kernels for second order irregular wave generation.
namespace cuwave2{

//-Kernels for JWaveSpectrum.
unsigned GetSizeAux(unsigned n);
double CalcPosition(double time,unsigned n,const double *dnm,const double2 *coefx,double *aux);
double CalcElevation(double time,double x,unsigned n,const double4 *coefe,double *aux);

//:double CalcPosition1(double time,unsigned n,const double *dnm,const double2 *coefx,double *aux);
//:double CalcPosition2(double time,unsigned n,const double *dnm,const double2 *coefx,double *aux);

}

#endif


