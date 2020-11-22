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

//::NO_COMENTARIO
//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para simplificar tareas de debug en GPU. (15-03-2017)
//:# - Simplificacion de clase usando JDataArrays. (05-09-2019)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:# - Cambios para soportar CODE_SIZE4. (07-04-2020)
//:# - Nueva variable TypeSp con el valor de special del code. (02-08-2020)
//:#############################################################################

/// \file JDebugSphGpu.h \brief Declares the class \ref JDebugSphGpu.

#ifndef _JDebugSphGpu_
#define _JDebugSphGpu_

#include "TypesDef.h"
#include "DualSphDef.h"

#include <string>
#include <cstring>
#include <cuda_runtime_api.h>


class JDataArrays;
class JSphGpuSingle;
class JSphMgpuNodeUnit;
class JSphDomain;

//##############################################################################
//# JDebugSphGpu
//##############################################################################
/// \brief Provides functions to store particle data in formats VTK, CSV, ASCII.

class JDebugSphGpu
{
protected:
  static void RunExceptioonStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,const std::string &msg,const std::string &file="");
  static void RunExceptioonCudaStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,cudaError_t cuerr,std::string msg);
  static void CheckCudaErroorStatic(const std::string &srcfile,int srcline
    ,const std::string &method
    ,std::string msg);

public:

  static byte*     GetCodeType     (unsigned n,const typecode *code);
  static typecode* GetCodeTypeValue(unsigned n,const typecode *code);
  static tuint3*   GetCell3(unsigned n,const unsigned *dcell,unsigned cellcode);
  static tfloat3*  GetPosf3(unsigned n,const tdouble3 *pos);
  static tfloat3*  GetPosf3(unsigned n,const tdouble2 *posxy,const double *posz);
  static tdouble3* GetPosd3(unsigned n,const tdouble2 *posxy,const double *posz);
  static tfloat3*  GetPosCell_Pos (unsigned n,const tfloat4 *poscell);
  static tuint3*   GetPosCell_Cell(unsigned n,const tfloat4 *poscell);

  static std::string PrepareVars(const std::string &vlist);
  static bool FindVar(const std::string &var,const std::string &vlist){ return(int(vlist.find(std::string(",")+var+","))>=0); }
  static std::string CheckVars(std::string vlist);

  static std::string GetFileName(std::string filename,int numfile,int gid=-1);

  static void LoadParticlesData(const JSphGpuSingle *gp,unsigned pini,unsigned pfin,std::string vars,JDataArrays *arrays,std::string file="");
  static void SaveVtk(std::string filename,int numfile,unsigned pini,unsigned pfin,std::string vars,const JSphGpuSingle *gp);

};


#endif


