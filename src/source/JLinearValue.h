//HEAD_DSPH
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

/// \file JLinearValue.h \brief Declares the class \ref JLinearValue.

#ifndef _JLinearValue_
#define _JLinearValue_

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase obtener valores interpolados linealmente de una lista. (05-01-2014)
//:# - Nuevo metodo LoadFile() para cargar valores de un fichero. (20-01-2014)
//:# - Se eliminaron includes innecesarios. (01-05-2014)
//:# - La funcion GetNewInterval() indica cuando se cambia de intervalo. (25-01-2017)
//:# - Puede gestionar varios valores para cada instante. (17-04-2017)
//:# - Nuevos metodos para calcular interpolacion de forma externa. (17-04-2017)
//:# - Nuevos metodos SetTimeValue() para modificar datos de la lista. (24-01-2018)
//:#############################################################################

#include "JObject.h"
#include "TypesDef.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>


//##############################################################################
//# JLinearValue
//##############################################################################
/// \brief Implements a list of values according time using linear interpolation.

class JLinearValue : protected JObject
{
protected:
  static const unsigned SIZEMAX=200000;
  static const unsigned SIZEINITIAL=500;

  std::string File;
  unsigned Size;
  unsigned Count;
  double *Times;
  double *Values;

  bool NewInterval;

  //-Variables with last search.
  double TimeStep;
  unsigned Position;
  unsigned PositionNext;
  double TimePre;
  double TimeNext;
  double TimeFactor;

public:
  const unsigned Nvalues;

  JLinearValue(unsigned nvalues=1);
  ~JLinearValue();
  void Reset();
  unsigned GetAllocMemory()const;

  void SetSize(unsigned size);
  unsigned GetSize()const{ return(Size); }

  unsigned AddTimeValue(double time,double value);
  unsigned AddTimeValue(double time,double value,double value2){ const unsigned idx=AddTimeValue(time,value); SetValue(idx,1,value2); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3){ const unsigned idx=AddTimeValue(time,value,value2); SetValue(idx,2,value3); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3,double value4){ const unsigned idx=AddTimeValue(time,value,value2,value3); SetValue(idx,3,value4); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3,double value4,double value5){ const unsigned idx=AddTimeValue(time,value,value2,value3,value4); SetValue(idx,4,value5); return(idx); }
  unsigned AddTimeValue(double time,double value,double value2,double value3,double value4,double value5,double value6){ const unsigned idx=AddTimeValue(time,value,value2,value3,value4,value5); SetValue(idx,5,value6); return(idx); }
  void SetValue(unsigned idx,unsigned cvalue,double value);

  void SetTimeValue(unsigned idx,double time,double value);
  void SetTimeValue(unsigned idx,double time,double value,double value2){  SetTimeValue(idx,time,value); SetValue(idx,1,value2);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3){  SetTimeValue(idx,time,value,value2); SetValue(idx,2,value3);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3,double value4){  SetTimeValue(idx,time,value,value2,value3); SetValue(idx,3,value4);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3,double value4,double value5){  SetTimeValue(idx,time,value,value2,value3,value4); SetValue(idx,4,value5);  }
  void SetTimeValue(unsigned idx,double time,double value,double value2,double value3,double value4,double value5,double value6){  SetTimeValue(idx,time,value,value2,value3,value4,value5); SetValue(idx,5,value6);  }

  void LoadFile(std::string file);
  std::string GetFile()const{ return(File); };

  unsigned GetCount()const{ return(Count); }
  void VisuData();

  //-For simple use.
  double GetValue(double timestep,unsigned cvalue=0);
  bool GetNewInterval()const{ return(NewInterval); }

  //-For external use.
  void FindTime(double timestep);
  unsigned GetPos()const{ return(Position); };
  unsigned GetPosNext()const{ return(PositionNext); };
  double GetTimeByIdx(unsigned idx)const;
  double GetValueByIdx(unsigned idx,unsigned cvalue=0)const;
};

#endif


