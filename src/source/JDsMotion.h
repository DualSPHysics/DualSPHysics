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
//:# - Para los ficheros de datos usa ruta absoluta si el nombre contiene alguna
//:#   barra de directorio. (11-09-2013)
//:# - Incluye la gestion de objetos moving. (23-04-2018)
//:# - Nuevo metodo GetObjIdxByMkBound(). (09-08-2018)
//:# - Nuevos metodos GetActiveMotion() y ProcesTimeGetData() simple. (19-09-2018)
//:# - Nuevo metodo SetMotionData(). (15-01-2019)
//:# - Mejora interface usando la estructura StMotionData. (12-09-2019)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Cambio de nombre de J.SphMotion a J.DsMotion. (28-06-2020)
//:# - En ProcesTime() siempre se calcula el movimiento para todos los objetos. (03-02-2021)
//:#############################################################################

/// \file JDsMotion.h \brief Declares the class \ref JDsMotion.

#ifndef _JDsMotion_
#define _JDsMotion_

#include "JObject.h"
#include "DualSphDef.h"
#include <string>

class JMotion;
class JXml;
class JCaseParts;

//##############################################################################
//# JDsMotion
//##############################################################################
/// \brief Provides the displacement for moving particles during a time interval.  

class JDsMotion : protected JObject
{
public:

  ///Controls the output of information on the screen and/or log.
  typedef enum{ 
    MOMT_Simple=0,  ///<Simple mode for only forward.
    MOMT_Ace2dt=1   ///<Calculates acceleration using one dt in the future (always from the beginning).
  }TpMotionMode;   

private:
  const bool Simulate2D;   ///<Toggles 2D simulation (cancels motion in Y axis).

  double TimeMod;          ///<Modifies the timestep for motion | Modificador del TimeStep para Motion.
  unsigned ObjCount;       ///<Number of moving objects.
  StMotionData *ObjMotion; ///<Motion data of moving objects.
  StMotionData MotionNull;

  //unsigned *ObjBegin;   ///<Initial particle of each moving object. [ObjCount+1]
  //word     *ObjMkBound; ///<MkBound of each moving object. [ObjCount]

  //byte      *ObjTpmov;     ///<Type of motion (0:none, 1:linear, 2:matrix, 3:ignore). [ObjCount]
  //tdouble3  *ObjLinMov;    ///<Linear motion. [ObjCount]
  //tmatrix4d *ObjMatMov;    ///<Matrix motion. [ObjCount]

  JMotion *Mot;
  bool ActiveMotion;    ///<Indicates active motions after executing ProcesTime().
  double LastDt;        ///<Dt used in last call to ProcesTime().
  void ConfigObjects(const JCaseParts *parts);

public:
  JDsMotion(bool simulate2d);
  ~JDsMotion();
  void Reset();
  void Init(const JCaseParts *parts,JXml *jxml,const std::string &path,const std::string &dirdata);

  unsigned GetNumObjects()const{ return(ObjCount); };

  unsigned GetObjIdxByMkBound(word mkbound)const;

  void SetTimeMod(double timemod){ TimeMod=timemod; };
  bool ProcesTime(TpMotionMode mode,double timestep,double dt);
  bool GetActiveMotion()const{ return(ActiveMotion); }

  const StMotionData& GetMotionData(unsigned ref)const;

  void SetMotionData   (const StMotionData& d);
  void SetMotionDataAce(const StMotionData& d);

  void SetMotionDataNone(unsigned idx);
  void SetMotionDataLin (unsigned idx,const tdouble3 &linmov);
  void SetMotionDataMat (unsigned idx,const tmatrix4d &matmov);

};

#endif


