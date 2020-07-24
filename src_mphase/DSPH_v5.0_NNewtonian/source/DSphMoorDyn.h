//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Permite compilar sin libreria DSphMoorDyn. (23-12-2019)
//:#############################################################################

/// \file DSphMoorDyn.h \brief Declares basic interface functions between MoorDyn+ and DualSPHysics.

#ifndef _DSphMoorDyn_
#define _DSphMoorDyn_

#ifndef DISABLE_DSPH
#include "DualSphDef.h"
#endif // !DISABLE_DSPH

#include "TypesDef.h"
#include "JLog2.h"
#include <string>

//#define DISABLE_MOORDYN     ///<It allows compile without DSphMoorDyn library.

#ifdef DISABLE_MOORDYN
#include "DSphMoorDynUndef.h"
#else
//==============================================================================
/// Initializes MoorDyn and returns true in case of error.
//==============================================================================
bool MoorDyn_LinesInit(const std::string filexml,const std::string nodexml
  ,const std::string dirout,const unsigned numFts,const unsigned ftmkbound[]
  ,const tdouble3 vellin[],const tdouble3 velang[],const tfloat3 gravity);

//==============================================================================
/// Deallocates the variables used by MoorDyn (returns true in case of error).
//==============================================================================
bool MoorDyn_LinesClose();

//==============================================================================
/// Force calculation of moorings by MoorDyn (returns true in case of error).
//==============================================================================
bool MoorDyn_FairleadsCalc(const unsigned numFts, double*** fairpos, double*** fairvel, double*** fairforce, double t, double dt);

//==============================================================================
/// Returns the tension at the fairlead of a given line. (line=0...)
//==============================================================================
double MoorDyn_GetFairTen(const unsigned line);

//==============================================================================
/// Returns number of fairlead connections (connections between floating and lines).
//==============================================================================
unsigned MoorDyn_FairsCount(const unsigned ftid);

//==============================================================================
/// Returns number of lines.
//==============================================================================
unsigned MoorDyn_LinesCount();

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDyn_SegsCount(const unsigned line);

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDyn_SegsCount(const unsigned ftid, const unsigned line);

//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDyn_GetNodePosLink(const unsigned ftid, const unsigned line);


//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDyn_GetNodePos(const unsigned line, const unsigned node);

//==============================================================================
/// Returns number of Moorings created
//==============================================================================
unsigned MoorDyn_MooringsCount();

//==============================================================================
/// Returns the mkbound of the Mooring
//==============================================================================
unsigned MoorDyn_GetMooringReference(const unsigned ftid);

//==============================================================================
/// Sends to MoorDyn the DualSPHysics log to store and print messages
//==============================================================================
void MoorDyn_LogInit(JLog2 * log);
#endif

#endif


