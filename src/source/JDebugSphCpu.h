//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Clase para simplificar tareas de debug en CPU. (10-09-2024)
//:#############################################################################

/// \file JDebugSphCpu.h \brief Declares the class \ref JDebugSphCpu.

#ifndef _JDebugSphCpu_
#define _JDebugSphCpu_

#include "TypesDef.h"
#include "DualSphDef.h"

#include <string>
#include <cstring>

class JDataArrays;
class JSphCpuSingle;

//##############################################################################
//# JDebugSphCpu
//##############################################################################
/// \brief Provides functions to store particle data in formats VTK, CSV, ASCII.

class JDebugSphCpu
{
protected:
  static void RunExceptioonStatic(const std::string& srcfile,int srcline
    ,const std::string& method
    ,const std::string& msg,const std::string& file="");

public:

  static byte*     GetCodeType     (unsigned n,const typecode* code);
  static typecode* GetCodeTypeValue(unsigned n,const typecode* code);
  static tuint3*   GetCell3(unsigned n,const unsigned* dcell,unsigned cellcode);
  static tfloat3*  GetPosf3(unsigned n,const tdouble3* pos);

  static std::string PrepareVars(const std::string& vlist);
  static bool FindVar(const std::string& var,const std::string& vlist){
    return(int(vlist.find(std::string(",")+var+","))>=0);
  }
  static std::string CheckVars(std::string vlist);

  static std::string GetFileName(std::string filename,int numfile);

  static void RunUserFilters(JDataArrays& arrays);

  static void LoadParticlesData(const JSphCpuSingle* cp,unsigned pini
    ,unsigned pfin,std::string vars,JDataArrays* arrays,std::string file="");
  static void SaveVtk(std::string filename,int numfile,unsigned pini
    ,unsigned pfin,std::string vars,const JSphCpuSingle* cp);
  static void SaveCsv(std::string filename,int numfile,unsigned pini
    ,unsigned pfin,std::string vars,const JSphCpuSingle* gp);

};


#endif


