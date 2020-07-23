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
//:# - Nuevas funciones para ReduSumFloat3. (20-11-2018)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

/// \file JReduSum_ker.h \brief Declares functions and CUDA kernels for reduction using using the sum.

#ifndef _JReduSum_ker_
#define _JReduSum_ker_

#include <cuda_runtime_api.h>

#ifndef REDUBSIZE
#define REDUBSIZE 256
#endif

//#define DG_curedus_Print          //-Displays information on the screen when it is active.
//#define DG_curedus_ReduSumDouble  //-In ReduSumDouble() checks that the result is correct.
//#define DG_curedus_ReduSumFloat   //-In ReduSumFloat() checks that the result is correct.
//#define DG_curedus_ReduSumUint    //-In ReduSumUint() checks that the result is correct.

/// Implements a set of functions and CUDA kernels for reduction using using the sum.
namespace curedus{

//-Kernels for ReduSumDouble.
unsigned GetAuxSize_ReduSumDouble(unsigned ndata);

double ReduSumDouble(unsigned ndata,unsigned inidata,const double* data,double* resu);
void ReduSumDoubleAsyn(unsigned ndata,unsigned inidata,const double* data,double* resu,cudaStream_t stm);
void ReduSumDoubleAsyn(unsigned ndata,unsigned inidata,const double* data,double* resu,double *pim1_sum,cudaStream_t stm);
double DgReduSumDouble(unsigned ndata,unsigned inidata,const double* data);

//-Kernels for ReduSumFloat.
unsigned GetAuxSize_ReduSumFloat(unsigned ndata);

float ReduSumFloat(unsigned ndata,unsigned inidata,const float* data,float* resu);
void ReduSumFloatAsyn(unsigned ndata,unsigned inidata,const float* data,float* resu,cudaStream_t stm);
void ReduSumFloatAsyn(unsigned ndata,unsigned inidata,const float* data,float* resu,float *pim1_sum,cudaStream_t stm);
float DgReduSumFloat(unsigned ndata,unsigned inidata,const float* data);

//-Kernels for ReduSumUint.
unsigned GetAuxSize_ReduSumUint(unsigned ndata);

unsigned ReduSumUint(unsigned ndata,unsigned inidata,const unsigned* data,unsigned* resu);
void ReduSumUintAsyn(unsigned ndata,unsigned inidata,const unsigned* data,unsigned* resu,cudaStream_t stm);
void ReduSumUintAsyn(unsigned ndata,unsigned inidata,const unsigned* data,unsigned* resu,unsigned *pim1_sum,cudaStream_t stm);
unsigned DgReduSumUint(unsigned ndata,unsigned inidata,const unsigned* data);

//-Kernels for ReduSumFloat3.
unsigned GetAuxSize_ReduSumFloat3(unsigned ndata);

float3 ReduSumFloat3(unsigned ndata,unsigned inidata,const float3* data,float3* resu);
void ReduSumFloat3Asyn(unsigned ndata,unsigned inidata,const float3* data,float3* resu,cudaStream_t stm);
void ReduSumFloat3Asyn(unsigned ndata,unsigned inidata,const float3* data,float3* resu,float3 *pim1_sum,cudaStream_t stm);
float3 DgReduSumFloat3(unsigned ndata,unsigned inidata,const float3* data);


}

#endif


