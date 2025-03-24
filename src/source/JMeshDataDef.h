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
//:# - Definiciones relacionadas con JMeshData clases. (03-08-2020)
//:#############################################################################

/// \file JMeshDataDef.h \brief Defines basic/general definitions and functions for JMeshData clases.

#ifndef _JMeshDataDef_
#define _JMeshDataDef_

#include "TypesDef.h"
#include <string>

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

///Structure with data definition.
typedef struct{
  std::string name;
  TpTypeData type;
  unsigned size;
  int tag;
}StMeshData;

/// Types of output format.
typedef enum{ 
  TpFmtNone=0      ///<None.
 ,TpFmtBin=1       ///<Binary format.
 ,TpFmtCsv=2       ///<CSV format.
}TpFormat; 

/// Returns the string key of output format type.
inline const char* GetTpFormatKey(const TpFormat t){
  switch(t){
    case TpFmtBin:  return("bin");
    case TpFmtCsv:  return("csv");
  }
  return("???");
}

/// Returns the string name of output format type.
inline const char* GetTpFormatName(const TpFormat t){
  switch(t){
    case TpFmtNone:  return("None");
    case TpFmtBin:   return("Binary");
    case TpFmtCsv:   return("CSV");
  }
  return("???");
}

///Structure with basic mesh configuration.
typedef struct{
  tdouble3 ptref;   ///<Initial mesh position.
  tdouble3 vec1;    ///<First axis vector to define the mesh.
  tdouble3 vec2;    ///<Second axis vector to define the mesh.
  tdouble3 vec3;    ///<Third axis vector to define the mesh grid (used for swl calculation).
  double dis1;      ///<Length for first axis.
  double dis2;      ///<Length for second axis.
  double dis3;      ///<Length for third axis.
  double dispt1;    ///<Distance between points for first axis.
  double dispt2;    ///<Distance between points for second axis.
  double dispt3;    ///<Distance between points for third axis.
  tfloat3 dirdat;   ///<Direction vector for computed linear velocity or other variables.
}StMeshBasic;

///Structure with mesh configuration for calculation.
typedef struct{
  tdouble3 ptref;   ///<Initial mesh position.
  tdouble3 vdp1;    ///<First axis vector with distance between points.
  tdouble3 vdp2;    ///<Second axis vector with distance between points.
  tdouble3 vdp3;    ///<Third axis vector with distance between points (used for swl calculation).
  unsigned npt1;    ///<Number of points for first axis.
  unsigned npt2;    ///<Number of points for second axis.
  unsigned npt3;    ///<Number of points for third axis.
  unsigned npt;     ///<Total number of points.
  tfloat3 dirdat;   ///<Direction vector for computed linear velocity or other variables.
}StMeshPts;

/// Returns StMeshPts with zeros.
inline StMeshPts MakeMeshPts(){
  StMeshPts v={{0,0,0},{0,0,0},{0,0,0},{0,0,0},0,0,0,0,{0,0,0}};
  return(v);
}

/// Returns StMeshPts with data.
inline StMeshPts MakeMeshPts(const tdouble3& ptref,unsigned npt1,const tdouble3& vdp1
  ,unsigned npt2,const tdouble3& vdp2,unsigned npt3,const tdouble3& vdp3
  ,const tfloat3& dirdat)
{
  StMeshPts v={ptref,vdp1,vdp2,vdp3,npt1,npt2,npt3,npt1*npt2*npt3,dirdat};
  return(v);
}

/// Returns StMeshPts with data in dimension 3.
inline StMeshPts MakeMeshPts3(const tdouble3& ptref,unsigned npt3,const tdouble3& vdp3
  ,const tfloat3& dirdat)
{
  StMeshPts v={ptref,{0,0,0},{0,0,0},vdp3,1,1,npt3,npt3,dirdat};
  return(v);
}

/// Returns true when some axis vector is valid.
inline bool CheckMeshVectors(const StMeshPts& g){
  return(g.vdp1!=TDouble3(0) || g.vdp2!=TDouble3(0) || g.vdp3!=TDouble3(0));
}


///Structure with basic data to load mesh-data velocity.
typedef struct StrMeshVelCfg{
  std::string file;     ///<File with mesh-data.
  tdouble3 setpos;      ///<Position offset (applied after filters).
  double initialtime;   ///<Defines initial time of data.
  bool velmagnitude;    ///<Uses velocity magnitude with Direction.
  tfloat3 veldir;       ///<Direction vector for computed linear velocity.
  bool velreverse;      ///<Reverses velocity data (v=-v).
  void Clear(){ 
    file=""; setpos=TDouble3(0); initialtime=0; 
    velmagnitude=false; veldir=TFloat3(0); velreverse=false;
  }
}StMeshVelCfg;

///Structure with basic data to load mesh-data density.
typedef struct StrMeshRhoCfg{
  std::string file;     ///<File with mesh-data.
  tdouble3 setpos;      ///<Position offset (applied after filters).
  double initialtime;   ///<Defines initial time of data.
  void Clear(){ 
    file=""; setpos=TDouble3(0); initialtime=0; 
  }
}StMeshRhoCfg;


///Structure with extended data to load mesh-data velocity.
typedef struct StrMeshVelExtCfg{
  std::string file;     ///<File with mesh-data.
  tdouble3 setpos;      ///<Applies a position offset to the mesh data (default=0,0,0).
  double initialtime;   ///<Defines the initial time of mesh data (default=0).
  double looptbeg;      ///<Defines the loop start (default=0).
  double looptmax;      ///<Defines the loop end. (default=0).
  bool velmagnitude;    ///<Velocity data using magnitude or 3 components (default=true).
  bool velreverse;      ///<Reverse velocity values (v=-v) (default=false).
  tdouble3 velmul;      ///<Multiply fixed value to velocity data (default=1,1,1).
  tdouble3 veladd;      ///<Add fixed value to velocity data (default=0,0,0).
  void Clear(){ 
    file=""; setpos=TDouble3(0); initialtime=0; 
    looptbeg=looptmax=0;
    velmagnitude=false; velreverse=false;
    velmul=TDouble3(1); veladd=TDouble3(0);
  }
}StMeshVelExtCfg;

}


#endif


