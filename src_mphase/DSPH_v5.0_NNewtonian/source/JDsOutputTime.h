//HEAD_DSPH
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
//:# - Gestiona el uso tiempo de grabacion de PARTs variable. (26-03-2016)
//:# - Ahora el constructor permite clonar otro objeto JDsOutputTime. (23-08-2019)
//:# - GetNextTime() guarda entrada y salida para evitar calculos con llamadas 
//:#   consecutivas iguales. (29-08-2019)
//:# - Objeto JXml pasado como const para operaciones de lectura. (18-03-2020)  
//:# - Improved exception managment. (18-03-2020)
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (19-03-2020)  
//:# - Cambio de nombre de J.TimeOut a J.DsOutputTime. (28-06-2020)
//:#############################################################################

/// \file JDsOutputTime.h \brief Declares the class \ref JDsOutputTime.

#ifndef _JDsTimeOut_
#define _JDsTimeOut_

#include "JObject.h"
#include <string>
#include <vector>

class JXml;
class TiXmlElement;
class JLog2;

//##############################################################################
//# XML format in _FmtXML_TimeOut.xml.
//##############################################################################

//##############################################################################
//# JDsOutputTime
//##############################################################################
/// \brief Manage the use of variable output time to save PARTs.

class JDsOutputTime : protected JObject
{
protected:
  ///Structure used to store information about timeout.
  typedef struct {
    double time;
    double tout;
  }StTimeOut;

  std::vector<StTimeOut> Times;  ///<List values for timeout.
  unsigned TimeBase;

  bool SpecialConfig; ///<Configuration loaded from XML file in special section.

  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void LoadXml(const JXml *sxml,const std::string &place);
  unsigned GetCount()const{ return(unsigned(Times.size())); }
  bool AddTimeOut(double t,double tout);
  void CopyFrom(const JDsOutputTime* tout);

  double LastTimeInput;   ///<Saves the last value used with GetNextTime().
  double LastTimeOutput;  ///<Saves the last value returned by GetNextTime().

public:
  JDsOutputTime(const JDsOutputTime* tout=NULL);
  ~JDsOutputTime();
  void Reset();
  void Config(double timeoutdef);
  void Config(std::string filexml,const std::string &place,double timeoutdef);
  bool UseSpecialConfig()const{ return(SpecialConfig); }
  void VisuConfig(std::string txhead,std::string txfoot);
  double GetNextTime(double t);
};

#endif


