//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2016, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

//:NO_COMENTARIO
//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase para extrapolar densidad en un determinado contorno. (03-05-2017)
//:# - Simple configuration for parallel boxes of boundary particles. (13-06-2018)
//:# - Saves VTK file (CfgBoundCorr_Limit.vtk) with LimitPos and Direction 
//:#   configuration. (13-06-2018)
//:#############################################################################

/// \file JSphBoundCorr.h \brief Declares the class \ref JSphBoundCorr.

#ifndef _JSphBoundCorr_
#define _JSphBoundCorr_

#include <string>
#include <vector>
#include "JObject.h"
#include "Types.h"

class JXml;
class TiXmlElement;
class JLog2;
class JLinearValue;
class JSphCpu;
class JSphMk;

//##############################################################################
//# XML format in _FmtXML_BoundCorr.xml.
//##############################################################################

//##############################################################################
//# JSphBoundCorrZone
//##############################################################################
/// \brief Manages one configuration for boundary extrapolated correction.
class JSphBoundCorrZone : protected JObject
{
public:
///Direction mode.
typedef enum{ 
    DIR_None=0,
    DIR_Top=1,
    DIR_Bottom=2,
    DIR_Left=3,
    DIR_Right=4,
    DIR_Front=5,
    DIR_Back=6
}TpDirection;  

private:
  JLog2 *Log;

  //-Selection of particles
  typecode BoundCode;      ///<Code to select boundary particles.

  //-Configuration parameters.
  TpDirection AutoDir; ///<Direction configuration for automatic definition.
  tdouble3 LimitPos;   ///<Limit between boundary and fluid.
  tdouble3 Direction;  ///<Direction to fluid particles.
  tfloat4 Plane;       ///<Plane in limit.

  void Reset();

public:
  const unsigned IdZone;
  const word MkBound;

  JSphBoundCorrZone(JLog2 *log,unsigned idzone,word mkbound
    ,TpDirection autodir,tdouble3 limitpos,tdouble3 direction);
  ~JSphBoundCorrZone();
  void ConfigBoundCode(typecode boundcode);
  void ConfigAutoLimit(double halfdp,tdouble3 pmin,tdouble3 pmax);

  void GetConfig(std::vector<std::string> &lines)const;

  TpDirection GetAutoDir()const{ return(AutoDir); }
  tdouble3 GetLimitPos()const{ return(LimitPos); }
  tdouble3 GetDirection()const{ return(Direction); }
  tfloat4 GetPlane()const{ return(Plane); }
  typecode GetBoundCode()const{ return(BoundCode); }
};

//##############################################################################
//# JSphBoundCorr
//##############################################################################
/// \brief Manages configurations for boundary extrapolated correction.
class JSphBoundCorr : protected JObject
{
private:
  JLog2 *Log;

  float DetermLimit;   ///<Limit for determinant. Use 1e-3 for first_order or 1e+3 for zeroth_order (default=1e+3).

  std::vector<JSphBoundCorrZone*> List; ///<List of configurations.

  void Reset();
  bool ExistMk(word mkbound)const;
  void LoadXml(JXml *sxml,const std::string &place);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void UpdateMkCode(const JSphMk *mkinfo);
  void SaveVtkConfig(double dp)const;

public:
  JSphBoundCorr(JLog2 *log,JXml *sxml,const std::string &place,const JSphMk *mkinfo);
  ~JSphBoundCorr();

  void RunAutoConfig(double dp,const JSphMk *mkinfo);

  void VisuConfig(std::string txhead,std::string txfoot)const;
  unsigned GetCount()const{ return(unsigned(List.size())); };

  float GetDetermLimit()const{ return(DetermLimit); };

  const JSphBoundCorrZone* GetMkZone(unsigned idx)const{ return(idx<GetCount()? List[idx]: NULL); }

};



#endif


