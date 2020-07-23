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

/// \file JWaveSpectrumGpu.cpp \brief Implements the class \ref JWaveSpectrumGpu.

#include "JWaveSpectrumGpu.h"
#ifdef _WITHGPU
  #include "JWaveOrder2_ker.h"
#endif

using namespace std;

//##############################################################################
//# JWaveSpectrumGpu
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JWaveSpectrumGpu::JWaveSpectrumGpu(){
  ClassName="JWaveSpectrumGpu";
  MemGpuFixed=0; Order2CoefsEtag=NULL; Order2CoefsDnmg=NULL; Order2CoefsPosg=NULL; Order2Auxg=NULL;
}

//==============================================================================
/// Frees GPU memory.
//==============================================================================
void JWaveSpectrumGpu::FreeMemoryGpu(){
  MemGpuFixed=0;
  #ifdef _WITHGPU
    if(Order2CoefsEtag)cudaFree(Order2CoefsEtag);   Order2CoefsEtag=NULL;
    if(Order2CoefsDnmg)cudaFree(Order2CoefsDnmg);   Order2CoefsDnmg=NULL;
    if(Order2CoefsPosg)cudaFree(Order2CoefsPosg);   Order2CoefsPosg=NULL;
    if(Order2Auxg)     cudaFree(Order2Auxg);        Order2Auxg=NULL;
  #endif 
}

//==============================================================================
/// Allocates GPU memory.
//==============================================================================
void JWaveSpectrumGpu::AllocMemoryGpu(unsigned sizewavecoefs){
  FreeMemoryGpu();
  #ifdef _WITHGPU
    MemGpuFixed=0;
    size_t m=sizeof(double4)*sizewavecoefs;
    cudaMalloc((void**)&Order2CoefsEtag,m);     MemGpuFixed+=m;
    m=sizeof(double)*sizewavecoefs;
    cudaMalloc((void**)&Order2CoefsDnmg,m);     MemGpuFixed+=m;
    m=sizeof(double2)*sizewavecoefs;
    cudaMalloc((void**)&Order2CoefsPosg,m);     MemGpuFixed+=m;
    m=sizeof(double)*cuwave2::GetSizeAux(sizewavecoefs);
    //m=sizeof(double)*SizeWaveCoefs*2+512; 
    cudaMalloc((void**)&Order2Auxg,m);          MemGpuFixed+=m;
  #endif
}

//==============================================================================
/// Copy coefficients to GPU memory.
//==============================================================================
void JWaveSpectrumGpu::CopyCoefs(unsigned sizewavecoefs,const tdouble4 *d4,const double *d1,const tdouble2 *d2){
  #ifdef _WITHGPU
    cudaMemcpy(Order2CoefsEtag,d4,sizeof(double4)*sizewavecoefs,cudaMemcpyHostToDevice);
    cudaMemcpy(Order2CoefsDnmg,d1,sizeof(double) *sizewavecoefs,cudaMemcpyHostToDevice);
    cudaMemcpy(Order2CoefsPosg,d2,sizeof(double2)*sizewavecoefs,cudaMemcpyHostToDevice);
  #endif
}

//==============================================================================
/// Returns paddle position using 2nd order wave theory.
//==============================================================================
double JWaveSpectrumGpu::CalcPosition(double time,unsigned sizewavecoefs){
  #ifdef _WITHGPU
    return(cuwave2::CalcPosition(time,sizewavecoefs,Order2CoefsDnmg,(double2*)Order2CoefsPosg,Order2Auxg));
  #else
    return(0);
  #endif
}

//==============================================================================
/// Returns surface elevation using 2nd order wave theory.
//==============================================================================
double JWaveSpectrumGpu::CalcElevation(double time,double x,unsigned sizewavecoefs){
  #ifdef _WITHGPU
    return(cuwave2::CalcElevation(time,x,sizewavecoefs,(double4*)Order2CoefsEtag,Order2Auxg));
  #else
    return(0);
  #endif
}




