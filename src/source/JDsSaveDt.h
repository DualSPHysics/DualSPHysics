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
//:# - Clase para generar fichero CSV con info del dt. (08-01-2015)
//:# - Permite grabar mas valores (csoundmax,acemax...). (14-01-2015)
//:# - Se elimina el valor csound. (22-01-2015)
//:# - Error corregido cuando TimeFinish=0. (22-01-2018)
//:# - Error corregido al generar excepcion por error en fichero. (22-01-2018)
//:# - Se escriben las unidades en las cabeceras de los ficheros CSV. (26-04-2018)
//:# - Improved exception managment. (19-03-2020)  
//:# - Objeto JXml pasado como const para operaciones de lectura. (19-03-2020)  
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (19-03-2020)  
//:# - Cambio de nombre de J.SaveDt a J.DsSaveDt. (28-06-2020)
//:#############################################################################

/// \file JDsSaveDt.h \brief Declares the class \ref JDsSaveDt.

#ifndef _JDsSaveDt_
#define _JDsSaveDt_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"

class JXml;
class TiXmlElement;
class JLog2;

//##############################################################################
//# XML format in _FmtXML_SaveDt.xml.
//##############################################################################

//##############################################################################
//# JDsSaveDt
//##############################################################################
/// \brief Manages the information of dt.

class JDsSaveDt : protected JObject
{
public:

/// Structure with dt information.
  typedef struct {
    double tini;
    unsigned num;
    double vmean;
    double vmin;
    double vmax;
  }StValue;

private:
  JLog2* Log;
  std::string FileDtInfo;
  std::string FileDtAllInfo;
  double TimeStart;    ///<Time from which information about the DT begins to be collected. | Instante a partir del cual se empieza a recopilar informacion del dt.
  double TimeFinish;   ///<Time from which dt information is not collected. | Instante a partir del cual se deja de recopilar informacion del dt.
  double TimeInterval; ///<Time lapse every time dt information is saved. | Cada cuanto se guarda info del dt.
  bool FullInfo;       ///<Saves AceMax, ViscDtMax and VelMax.
  bool AllDt;
  unsigned SizeValuesSave;

  StValue ValueNull;

  unsigned Count;                       ///<Number of stored intervals. | Numero de intervalos almacenados.
  static const unsigned SizeValues=100; ///<Maximum number of intervals to be buffered. | Numero maximo de intervalos a almacenar en buffer.
  StValue DtFinal[SizeValues];          ///<Resultant minimum Dt [SizeValues]. | Dt minimo resultante [SizeValues].
  StValue Dt1[SizeValues];              ///<Dt1 [SizeValues].
  StValue Dt2[SizeValues];              ///<Dt2 [SizeValues].
  StValue AceMax[SizeValues];           ///<AceMax [SizeValues].
  StValue ViscDtMax[SizeValues];        ///<ViscDtMax [SizeValues].
  StValue VelMax[SizeValues];           ///<VelMax [SizeValues].

  unsigned GetSizeValues()const{ return(SizeValues); }

  static const unsigned SizeAllDts=1000;
  tdouble2 AllDts[SizeAllDts];           
  unsigned CountAllDts;

  unsigned LastInterval;
  StValue LastDtf,LastDt1,LastDt2;
  StValue LastAceMax,LastViscDtMax,LastVelMax;

  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void LoadXml(const JXml *sxml,const std::string &place);
  void SaveFileValues();
  void SaveFileValuesEnd();
  void SaveFileAllDts();
  void AddValueData(double timestep,double dt,StValue &value);
  void AddLastValues();

public:
  JDsSaveDt();
  ~JDsSaveDt();
  void Reset();
  void Config(const JXml *sxml,const std::string &place,double timemax,double timeout);
  void VisuConfig(std::string txhead,std::string txfoot);
  void AddValues(double timestep,double dtfinal,double dt1,double dt2,double acemax,double viscdtmax,double velmax);
  bool GetFullInfo()const{ return(FullInfo); }
  void SaveData();
};


#endif


