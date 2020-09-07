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

/// \file JSphInOutDef.h \brief Defines basic/general definitions and functions for Inlet/Outlet code.

#ifndef _JSphInOutDef_
#define _JSphInOutDef_

//#include "TypesDef.h"


///Defines the treatment of fluid particles entering a inlet/outlet zone.
typedef enum{ 
   InInput_Free=0     ///<Free input mode (no changes).
  ,InInput_Convert=1  ///<Convert fluid to inlet/outlet.
  ,InInput_Remove=2   ///<Remove fluid.
}TpInInput;

///Defines refilling modes.
typedef enum{ 
   InRefill_SimpleFull=0   ///<It creates new inout particles when inout particles is moved to normal fluid.
  ,InRefill_SimpleZsurf=1  ///<It only creates new inout particles below zsurf.
  ,InRefill_Advanced=2     ///<Advanced for reverse flows (slow).
}TpInRefilling;

///Controls imposed velocity.
typedef enum{ 
   InVelM_Fixed=0x00        ///<Imposed fixed velocity (xxxx 00xx).
  ,InVelM_Variable=0x04     ///<Imposed variable velocity (xxxx 01xx).
  ,InVelM_Extrapolated=0x08 ///<Extrapolated from ghost nodes (xxxx 10xx).
  ,InVelM_Interpolated=0x0C ///<Interpolated velocity (xxxx 11xx).
  ,InVelM_MASK=0x0C          ///<Mask to obtain value (0000 1100).
}TpInVelMode;   

/// Returns the name of imposed velocity mode.
inline const char* TpInVelModeText(const TpInVelMode t){
  switch(t){
    case InVelM_Fixed:        return("Fixed");
    case InVelM_Variable:     return("Variable");
    case InVelM_Extrapolated: return("Extrapolated");
    case InVelM_Interpolated: return("Interpolated");
  }
  return("???");
}

///Controls profile of imposed velocity.
typedef enum{ 
   InVelP_Uniform=0x00     ///<Imposed velocity profile uniform (xxxx xx00).
  ,InVelP_Linear=0x01      ///<Imposed velocity profile linear (xxxx xx01).
  ,InVelP_Parabolic=0x02   ///<Imposed velocity profile parabolic (xxxx xx10).
  ,InVelP_MASK=0x03        ///<Mask to obtain value (0000 0011).
}TpInVelProfile;

///Behaviour of inlet/outlet according to the velocity.
typedef enum{ 
   InVelB_Unknown=0  ///<Unknown behaviour.
  ,InVelB_Inlet=1    ///<Inlet behaviour (velocity>=0).
  ,InVelB_Outlet=2   ///<Outlet behaviour (velocity<0).
  ,InVelB_Reverse=3  ///<Inlet & Outlet behaviour (velocity>0 and velocity<0).
}TpInBehaviour;

///Controls imposed rhop.
typedef enum{ 
   InRhop_Constant=0x00      ///<Imposed rhop profile constant (xx00 xxxx).
  ,InRhop_Hydrostatic=0x10   ///<Imposed rhop profile hydrostatic (xx01 xxxx).
  ,InRhop_Extrapolated=0x20  ///<Extrapolated from ghost nodes (xx10 xxxx).
  ,InRhop_MASK=0x30          ///<Mask to obtain value (0011 0000).
}TpInRhopMode;   

///Controls imposed Z surface.
typedef enum{ 
  InZsurf_Undefined=0,    ///<Imposed undefined zsurf.
  InZsurf_Fixed=1,        ///<Imposed fixed zsurf.
  InZsurf_Variable=2,     ///<Imposed variable zsurf.
  InZsurf_Calculated=3,   ///<Zsurf is calculated from fluid domain.
}TpInZsurfMode;   

/// Returns the name of zsurf mode.
inline const char* TpInZsurfModeText(const TpInZsurfMode t){
  switch(t){
    case InZsurf_Undefined:  return("Undefined");
    case InZsurf_Fixed:      return("Fixed");
    case InZsurf_Variable:   return("Variable");
    case InZsurf_Calculated: return("Calculated");
  }
  return("???");
}

///Structure with zsurf results.
typedef struct{
  unsigned npt;        ///<Number of zsurf values (uniform=1).
  const float *zsurf;  ///<Pointer to zsurf values [npt].
  tdouble3 pt;         ///<Reference point of gauge points.
  tdouble3 vdp;        ///<Vector to compute gauge positions.
  tdouble3 direction;  ///<Inflow direction.
}StZsurfResult;


#endif


