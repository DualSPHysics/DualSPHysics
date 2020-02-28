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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Se crea para unificar definiciones de condiciones periodicas. (27-04-2018)
//:#############################################################################

/// \file JPeriodicDef.h \brief Defines type of periodic conditions and basic functionalities.

#ifndef _JPeriodicDef_
#define _JPeriodicDef_

#include "TypesDef.h"
#include <string>

/// Types of periodic conditions.
typedef enum{ 
  PERI_None=0,      ///<No periodic conditions.
  PERI_X=1,         ///<Periodic conditions on axis X.
  PERI_Y=2,         ///<Periodic conditions on axis Y.
  PERI_Z=4,         ///<Periodic conditions on axis Z.
  PERI_XY=3,        ///<Periodic conditions on axis X and Y.
  PERI_XZ=5,        ///<Periodic conditions on axis X and Z.
  PERI_YZ=6,        ///<Periodic conditions on axis Y and Z.
  PERI_Unknown=96   ///<Unknown periodic conditions.
}TpPeri; 

/// Returns the string name of periodic condition type.
inline const char* TpPeriName(TpPeri tperi){
  switch(tperi){
  case PERI_None:    return("None");
  case PERI_X:       return("Axis-X");
  case PERI_Y:       return("Axis-Y");
  case PERI_Z:       return("Axis-Z");
  case PERI_XY:      return("Axes-XY");
  case PERI_XZ:      return("Axes-XZ");
  case PERI_YZ:      return("Axes-YZ");
  case PERI_Unknown: return("Unknown");
  }
  return("???");
}

#define PERI_AxisX(periactive) (periactive&1)
#define PERI_AxisY(periactive) (periactive&2)
#define PERI_AxisZ(periactive) (periactive&4)

#define PERI_Axis_X(periactive) ((periactive&1)!=0)
#define PERI_Axis_Y(periactive) ((periactive&2)!=0)
#define PERI_Axis_Z(periactive) ((periactive&4)!=0)

/// Returns PeriActive value according periodic axes.
inline byte DefPeriActive(bool perix,bool periy,bool periz){ return(byte((perix? 1: 0)+(periy? 2: 0)+(periz? 4: 0))); }

/// Returns PeriActive value according periodic axes.
inline byte DefPeriActive(TpPeri tperi){ return(tperi==PERI_Unknown? 0: byte(tperi)); }

/// Returns TpPeri value according periactive value.
inline TpPeri TpPeriFromPeriActive(byte periactive){
  return(periactive==DefPeriActive(PERI_Axis_X(periactive),PERI_Axis_Y(periactive),PERI_Axis_Z(periactive))? TpPeri(periactive): PERI_Unknown);
}


/// Structure with Periodic information.
typedef struct StrPeriodic{
  TpPeri PeriMode;
  tdouble3 PeriXinc;   ///<Value that is added at the outer limit to modify coordinates.
  tdouble3 PeriYinc;   ///<Value that is added at the outer limit to modify coordinates.
  tdouble3 PeriZinc;   ///<Value that is added at the outer limit to modify coordinates.

  //-Constructor by default.
  StrPeriodic(){  Reset(); }
  //-Constructor.
  StrPeriodic(byte periactive,tdouble3 xinc,tdouble3 yinc,tdouble3 zinc){
    PeriMode=TpPeriFromPeriActive(periactive);
    PeriXinc=xinc;
    PeriYinc=yinc;
    PeriZinc=zinc;
  }
  //-Reset values.
  void Reset(){
    PeriMode=PERI_Unknown;
    PeriXinc=PeriYinc=PeriZinc=TDouble3(0);
  }
}StPeriodic;



#endif


