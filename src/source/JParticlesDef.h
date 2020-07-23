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
//:# - Se crea para unificar definiciones de particulas. (23-03-2018)
//:#############################################################################

/// \file JParticlesDef.h \brief Defines type of particles and basic methods for particles.

#ifndef _JParticlesDef_
#define _JParticlesDef_

#include <string>

/// Types of particles. Order must be maintained although specific values can be changed.
/// Tipos de particulas. Se tiene que mantener el orden aunque se pueden cambiar los valores concretos.
typedef enum{ 
  TpPartFixed=0       ///<Fixed boundary particles.
 ,TpPartMoving=1      ///<Moving boundary particles.
 ,TpPartFloating=2    ///<Floating boundary particles.
 ,TpPartFluid=3       ///<Fluid particles.
 ,TpPartUnknown=9     ///<Unknown or undefined type of particles.
}TpParticles; 


#define TPPARTICLES_COUNT 4   ///<Number of valid particle types.

/// Returns if particle type is boundary or not.
inline bool IsBound(const TpParticles type){ return(type<TpPartFluid); }

/// Returns if particle type is fluid or not.
inline bool IsFluid(const TpParticles type){ return(type==TpPartFluid); }


/// Returns the string code of particle type.
inline const char* TpPartGetStrCode(const TpParticles type){
  switch(type){
    case TpPartFixed:     return("Fixed");
    case TpPartMoving:    return("Moving");
    case TpPartFloating:  return("Floating");
    case TpPartFluid:     return("Fluid");
  }
  return("???");
}

/// Returns particle type according the string code of particle type.
inline TpParticles TpPartGetType(std::string strcode){
  if(strcode=="Fixed")   return(TpPartFixed);
  if(strcode=="Moving")  return(TpPartMoving);
  if(strcode=="Floating")return(TpPartFloating);
  if(strcode=="Fluid")   return(TpPartFluid);
  return(TpPartUnknown);
}



/////Types of particles (OLD).
//typedef enum{ 
//  PART_BoundFx=1,          ///<Fixed boundary particles.
//  PART_BoundMv=2,          ///<Moving boundary particles.
//  PART_BoundFx_BoundMv=3,  ///<Both fixed and moving boundary particles.
//  PART_BoundFt=4,          ///<Floating boundary particles.
//  PART_Fluid=8,            ///<Fluid particles.
//  PART_BoundFt_Fluid=12    ///<Both floating and fluid particles.
//}TpParticle;


#endif


