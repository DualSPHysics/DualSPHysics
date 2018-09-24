//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Para los ficheros de datos usa ruta absoluta si el nombre contiene alguna
//:#   barra de directorio. (11-09-2013)
//:# - Incluye la gestion de objetos moving. (23-04-2018)
//:# - Nuevo metodo GetObjIdxByMkBound(). (09-08-2018)
//:# - Nuevos metodos GetActiveMotion() y ProcesTimeGetData() simple. (19-09-2018)
//:#############################################################################

/// \file JSphMotion.h \brief Declares the class \ref JSphMotion.

#ifndef _JSphMotion_
#define _JSphMotion_

#include "JObject.h"
#include "Types.h"
#include <string>

class JMotion;
class JXml;
class JSpaceParts;

//##############################################################################
//# JSphMotion
//##############################################################################
/// \brief Provides the displacement for moving particles during a time interval.  

class JSphMotion : protected JObject
{
public:

  ///Controls the output of information on the screen and/or log.
  typedef enum{ 
    MOMT_Simple=0,  ///<Simple mode for only forward.
    MOMT_Ace2dt=1,  ///<Calculates acceleration using one dt in the future (always from the beginning).
  }TpMotionMode;   

private:
  double TimeMod;       ///<Modifies the timestep for motion | Modificador del TimeStep para Motion.
  unsigned ObjCount;    ///<Number of moving objects.
  unsigned *ObjBegin;   ///<Initial particle of each moving object. [ObjCount+1]
  word     *ObjMkBound; ///<MkBound of each moving object. [ObjCount]

  JMotion *Mot;
  bool ActiveMotion;    ///<Indicates active motions after executing ProcesTime().
  void ConfigObjects(const JSpaceParts *parts);

public:
  JSphMotion();
  ~JSphMotion();
  void Reset();
  void Init(const JSpaceParts *parts,JXml *jxml,const std::string &path,const std::string &dirdata);

  unsigned GetNumObjects()const{ return(ObjCount); };
  word GetObjMkBound(unsigned idx)const;
  unsigned GetObjBegin(unsigned idx)const;
  unsigned GetObjSize(unsigned idx)const;

  unsigned GetObjIdxByMkBound(word mkbound)const;

  void SetTimeMod(double timemod){ TimeMod=timemod; };
  bool ProcesTime(TpMotionMode mode,double timestep,double dt);
  bool GetActiveMotion()const{ return(ActiveMotion); }
  bool ProcesTimeGetData(unsigned ref,bool &typesimple,tdouble3 &simplemov
    ,tdouble3 &simplevel,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2
    ,unsigned &nparts,unsigned &idbegin)const;
  bool ProcesTimeGetData(unsigned ref,word &mkbound
    ,bool &typesimple,tdouble3 &simplemov,tmatrix4d &matmov)const;
};

#endif


