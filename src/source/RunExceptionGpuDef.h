//HEAD_DSCODES
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

//-Defines for CUDA exceptions.
#ifndef Run_ExceptioonCuda
#define Run_ExceptioonCuda(cuerr,msg) RunExceptioonCuda(__FILE__,__LINE__,ClassName,__func__,cuerr,msg)
#endif
#ifndef Run_ExceptioonCudaSta
#define Run_ExceptioonCudaSta(cuerr,msg) RunExceptioonCudaStatic(__FILE__,__LINE__,__func__,cuerr,msg)
#endif
#ifndef Check_CudaErroor
#define Check_CudaErroor(msg) CheckCudaErroor(__FILE__,__LINE__,ClassName,__func__,msg)
#endif
#ifndef Check_CudaErroorSta
#define Check_CudaErroorSta(msg) CheckCudaErroorStatic(__FILE__,__LINE__,__func__,msg)
#endif
#ifndef Check_CudaErroorFun
#define Check_CudaErroorFun(msg) CheckCudaErroorFun(__FILE__,__LINE__,__func__,msg)
#endif

