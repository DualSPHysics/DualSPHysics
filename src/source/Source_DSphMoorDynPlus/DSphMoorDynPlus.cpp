//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics.

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.
*/

/// \file DSphMoorDynPlus.cpp \brief Implements basic interface functions between MoorDynPlus and DualSPHysics.

#include "DSphMoorDynPlus.h"
#include "MoorDynPlus.h"

MoorDynPlus moordynplus;

//==============================================================================
/// Initializes MoorDynPlus (returns true in case of error).
//==============================================================================
bool MoorDynPlus_LinesInit(const std::string filexml,const std::string nodexml
                      ,const std::string dirout,const unsigned numFt
                      ,const unsigned ftmkbound[],const tdouble3 vellin[]
                      ,const tdouble3 velang[],const tfloat3 gravity
                      ,const double tmax,const double dtout) {
  // double pos[6]={center.x,center.y,center.z,angles.x,angles.y,angles.z};
  // double vel[6]={vellin.x,vellin.y,vellin.z,velang.x,velang.y,velang.z};
  double pos[6]={ 0,0,0,0,0,0 };
  double vel[6]={ 0,0,0,0,0,0 };
  return (moordynplus.LinesInit(pos,vel,filexml,nodexml,dirout.c_str(),gravity,tmax,dtout,ftmkbound,numFt) != 0);
}

//==============================================================================
/// Deallocates the variables used by MoorDynPlus (returns true in case of error).
//==============================================================================
bool MoorDynPlus_LinesClose() {
  return (moordynplus.LinesClose() != 0);
}

//==============================================================================
/// Force calculation of moorings by MoorDynPlus (returns true in case of error).
//==============================================================================
bool MoorDynPlus_FairleadsCalc(const unsigned numFts,double*** fairpos,double*** fairvel,double*** fairforce,double t,double dt) {
  return (moordynplus.FairleadsCalc(numFts,fairpos,fairvel,fairforce,&t,&dt)!=0);
}

//==============================================================================
/// Returns the tension at the fairlead of a given line. (line=0...)
//==============================================================================
double MoorDynPlus_GetFairTen(unsigned line) {
  return (moordynplus.GetFairTen(line));
}

//==============================================================================
/// Returns number of fairlead connections.
//==============================================================================
unsigned MoorDynPlus_FairsCount(const unsigned ftid) {
  return (moordynplus.GetNFairs(ftid));
}

//==============================================================================
/// Returns number of lines.
//==============================================================================
unsigned MoorDynPlus_LinesCount() {
  return (moordynplus.GetNLines());
}

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDynPlus_SegsCount(const unsigned line) {
  return (static_cast<unsigned>(moordynplus.GetSegsCount(line)));
}

//==============================================================================
/// Returns number of segments in indicated line. (line=0...)
//==============================================================================
unsigned MoorDynPlus_SegsCount(const unsigned ftid,const unsigned line) {
  return (moordynplus.GetSegsCount(ftid,line));
}

//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDynPlus_GetNodePos(const unsigned line,const unsigned node) {
  return (moordynplus.GetNodePos(line,node));
}

//==============================================================================
/// Returns position of node. (line=0... node=0...)
//==============================================================================
tdouble3 MoorDynPlus_GetNodePosLink(const unsigned ftid,const unsigned line) {
  return (moordynplus.GetNodePosLink(ftid,line));
}

//==============================================================================
/// Returns number of Moorings created
//==============================================================================
unsigned MoorDynPlus_MooringsCount() {
  return (moordynplus.GetNBodies());
}

//==============================================================================
/// Returns the mkbound of the Mooring
//==============================================================================
unsigned MoorDynPlus_GetMooringReference(const unsigned ftid) {
  return (moordynplus.GetBodyReference(ftid));
}

//==============================================================================
/// Sends to MoorDynPlus the object to manage Logs from DualSPHysics
//==============================================================================
void MoorDynPlus_LogInit(JLog2* log) {
  moordynplus.LogInit(log);
}
