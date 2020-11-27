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

//-Defines for normal exceptions in general functions.
#ifndef Run_ExceptioonFun
#define Run_ExceptioonFun(msg) RunExceptioonFun(__FILE__,__LINE__,__func__,msg)
#endif
#ifndef Run_ExceptioonFileFun
#define Run_ExceptioonFileFun(msg,file) RunExceptioonFun(__FILE__,__LINE__,__func__,msg,file)
#endif

//-Defines for normal exceptions in static methods.
#ifndef Run_ExceptioonSta
#define Run_ExceptioonSta(msg) RunExceptioonStatic(__FILE__,__LINE__,__func__,msg)
#endif
#ifndef Run_ExceptioonFileSta
#define Run_ExceptioonFileSta(msg,file) RunExceptioonStatic(__FILE__,__LINE__,__func__,msg,file)
#endif

//-Defines for normal exceptions in class methods.
#ifndef Run_Exceptioon
#define Run_Exceptioon(msg) RunExceptioon(__FILE__,__LINE__,ClassName,__func__,msg)
#endif
#ifndef Run_ExceptioonFile
#define Run_ExceptioonFile(msg,file) RunExceptioon(__FILE__,__LINE__,ClassName,__func__,msg,file)
#endif

