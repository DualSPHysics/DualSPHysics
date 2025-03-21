//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics.

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.
*/

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Coupling para release v5.0. (23-12-2019)
//:# - Permite compilar sin libreria DSphMoorDynPlus. (23-12-2019)
//:#############################################################################

/// \file DSphMoorDynPlus.h \brief Declares basic interface functions between MoorDynPlus and DualSPHysics.

#ifndef _DSphMoorDynPlus_
#define _DSphMoorDynPlus_

#ifndef DISABLE_DSPH
#include "DualSphDef.h"
#endif // !DISABLE_DSPH

#include "TypesDef.h"
#include "JLog2.h"
#include "JXml.h"
#include <string>

//#define DISABLE_MOORDYNPLUS     ///<It allows compile without DSphMoorDynPlus library.

#ifdef DISABLE_MOORDYNPLUS
#include "DSphMoorDynPlusUndef.h"
#else
//==============================================================================
/// Initializes MoorDynPlus and returns true in case of error.
//==============================================================================
bool MoorDynPlus_LinesInit(const std::string filexml,const std::string nodexml
  ,const std::string dirout,const unsigned numFts,const unsigned ftmkbound[]
  ,const tdouble3 vellin[],const tdouble3 velang[],const tfloat3 gravity
  ,const double tmax,const double dtout);

//==============================================================================
/// Deallocates the variables used by MoorDynPlus (returns true in case of error).
//==============================================================================
bool MoorDynPlus_LinesClose();

//==============================================================================
/// Force calculation of moorings by MoorDynPlus (returns true in case of error).
//==============================================================================
bool MoorDynPlus_FairleadsCalc(const unsigned numFts,double*** fairpos,double*** fairvel
  ,double*** fairforce,double t,double dt);

//==============================================================================
/// Returns the tension at the fairlead of a given line. (line=0...)
//==============================================================================
double MoorDynPlus_GetFairTen(const unsigned line);

//==============================================================================
/// Returns number of fairlead connections (connections between floating and lines).
//==============================================================================
unsigned MoorDynPlus_FairsCount(const unsigned ftid);

//==============================================================================
/// Returns number of lines.
//==============================================================================
unsigned MoorDynPlus_LinesCount();

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDynPlus_SegsCount(const unsigned line);

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDynPlus_SegsCount(const unsigned ftid,const unsigned line);

//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDynPlus_GetNodePosLink(const unsigned ftid,const unsigned line);


//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDynPlus_GetNodePos(const unsigned line,const unsigned node);

//==============================================================================
/// Returns number of Moorings created
//==============================================================================
unsigned MoorDynPlus_MooringsCount();

//==============================================================================
/// Returns the mkbound of the Mooring
//==============================================================================
unsigned MoorDynPlus_GetMooringReference(const unsigned ftid);

//==============================================================================
/// Sends to MoorDynPlus the object to manage Logs from DualSPHysics
//==============================================================================
void MoorDynPlus_LogInit(JLog2* log);

#endif
#endif


