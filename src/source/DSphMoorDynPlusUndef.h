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

/// \file DSphMoorDynPlusUndef.h \brief Declares empty interface functions between MoorDynPlus and DualSPHysics.

#ifndef _DSphMoorDynPlusUndef_
#define _DSphMoorDynPlusUndef_

#ifdef DISABLE_MOORDYNPLUS
bool MoorDynPlus_LinesInit(const std::string filexml,const std::string nodexml
  ,const std::string dirout,const unsigned numFts,const unsigned ftmkbound[]
  ,const tdouble3 vellin[],const tdouble3 velang[],const tfloat3 gravity
  ,const double tmax,const double dtout){ return(true); }
bool MoorDynPlus_LinesClose(){ return(true); }
bool MoorDynPlus_FairleadsCalc(const unsigned numFts,double*** fairpos,double*** fairvel,double*** fairforce,double t,double dt){ return(true); }
double MoorDynPlus_GetFairTen(const unsigned line){ return(0); }
unsigned MoorDynPlus_FairsCount(const unsigned ftid){ return(0); }
unsigned MoorDynPlus_LinesCount(){ return(0); }
unsigned MoorDynPlus_SegsCount(const unsigned line){ return(0); }
unsigned MoorDynPlus_SegsCount(const unsigned ftid, const unsigned line){ return(0); }
tdouble3 MoorDynPlus_GetNodePosLink(const unsigned ftid, const unsigned line){ return(TDouble3(0)); }
tdouble3 MoorDynPlus_GetNodePos(const unsigned line, const unsigned node){ return(TDouble3(0)); }
unsigned MoorDynPlus_MooringsCount(){ return(0); }
unsigned MoorDynPlus_GetMooringReference(const unsigned ftid){ return(0); }
void MoorDynPlus_LogInit(JLog2* log){}
#endif
#endif


