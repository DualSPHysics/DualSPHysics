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

/// \file DSphMoorDynUndef.h \brief Declares empty interface functions between MoorDyn+ and DualSPHysics.

#ifndef _DSphMoorDynUndef_
#define _DSphMoorDynUndef_

#ifdef DISABLE_MOORDYN
bool MoorDyn_LinesInit(const std::string filexml,const std::string nodexml
  ,const std::string dirout,const unsigned numFts,const unsigned ftmkbound[]
  ,const tdouble3 vellin[],const tdouble3 velang[],const tfloat3 gravity){ return(true); }
bool MoorDyn_LinesClose(){ return(true); }
bool MoorDyn_FairleadsCalc(const unsigned numFts, double*** fairpos, double*** fairvel, double*** fairforce, double t, double dt){ return(true); }
double MoorDyn_GetFairTen(const unsigned line){ return(0); }
unsigned MoorDyn_FairsCount(const unsigned ftid){ return(0); }
unsigned MoorDyn_LinesCount(){ return(0); }
unsigned MoorDyn_SegsCount(const unsigned line){ return(0); }
unsigned MoorDyn_SegsCount(const unsigned ftid, const unsigned line){ return(0); }
tdouble3 MoorDyn_GetNodePosLink(const unsigned ftid, const unsigned line){ return(TDouble3(0)); }
tdouble3 MoorDyn_GetNodePos(const unsigned line, const unsigned node){ return(TDouble3(0)); }
unsigned MoorDyn_MooringsCount(){ return(0); }
unsigned MoorDyn_GetMooringReference(const unsigned ftid){ return(0); }
void MoorDyn_LogInit(JLog2 * log){}
#endif

#endif


