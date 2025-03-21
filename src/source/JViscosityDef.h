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

/// \file JViscosityDef.h \brief Defines SPH viscosity formulations and basic methods.

#ifndef _JViscosityDef_
#define _JViscosityDef_

///Types of SPH viscosity treatment.
typedef enum{ 
  VISCO_Laminar=3,     ///<Laminar viscosity.
  VISCO_LaminarSPS=2,  ///<Laminar viscosity and Sub-Partice Scale Turbulence.
  VISCO_Artificial=1,  ///<Artificial viscosity.
  VISCO_None=0 
}TpVisco;

///Returns the name of the viscosity formulation.
inline const char* GetViscoName(TpVisco tvisco){
  switch(tvisco){
    case VISCO_Laminar:    return("Laminar");
    case VISCO_LaminarSPS: return("Laminar+SPS");
    case VISCO_Artificial: return("Artificial");
  }
  return("???");
}

#endif

