//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Gestiona el calculo automatico de BlockSize para ejecucion de kernels GPU. (18-02-2016)
//:# - Uso de JSaveCsv2 para generar ficheros csv. (12-03-2018)
//:#############################################################################

/// \file JBlockSizeAuto.h \brief Declares the class \ref JBlockSizeAuto.

#ifndef _JBlockSizeAuto_
#define _JBlockSizeAuto_

#include "JObject.h"
#include "TypesDef.h"
#include "JMeanValues.h"
#include <string>
#include <vector>

class JLog2;

//##############################################################################
//# JBlockSizeAutoKer
//##############################################################################
/// \brief Manages the automatic computation of optimum BlockSize for each kernel interaction.

class JBlockSizeAutoKer : protected JObject
{
public:
  const std::string Name; ///< Kernel name.
  const int BsDef;  ///< Default BlockSize.
  const int BsMin;  ///< Minimum BlockSize.
  const int BsInc;  ///< Jump between two blocksize values.
  const int BsNum;  ///< Number of BlockSize to test.
protected:
  JLog2 *Log;

  unsigned Nrun;   ///< Number of executions.
  unsigned BsSel;  ///< Indicates optimum Blocksize.

  static const int REMOVESTART=5; //- Despues de cuantas ejecuciones empieza a descartar valores.
  static const int REMOVEPRC=20;  //- Porcentaje de descartes por ejecucion.
  static const int REMOVELIMIT=5; //- Numero minimo de valores sin descartar.

  int NumActive;       ///< Number of active values.
  bool *BsActive;      ///< Indicates if Blocksize is active [BsNum].
  float *Times;        ///< Saves times of test [BsNum].
  float *OverMean;     ///< Stores exponential mean overhead [BsNum].

  static const int MEANDEPTH=10;
  JMeanValue *MeanTot;
  JMeanMoving *MeanExp;

  static const int SAVEINFO=0; ///< Saves statistical data in CSV format or not.
  bool InfoDataSaved;          ///< Indicates if data was saved.
  unsigned InfoDataSizeLine;   ///< Number of floats per line.
  unsigned InfoDataLines;      ///< Number of lines for which memory was allocated.
  unsigned InfoDataCount;      ///< Number of used lines.
  float *InfoData;             ///< Buffer to store values [InfoDataSizeLine*InfoDataLines].

  void AllocateMemory(unsigned size);
  void AllocateInfoData(unsigned nlines);
  void SaveFileInfoData();
  void SaveInfoData(unsigned nstep,float timestep);


public:
  JBlockSizeAutoKer(JLog2 *log,std::string name,int bsmin,int bsnum,int bsinc,int bsdef);
  ~JBlockSizeAutoKer();
  void Reset();

  bool IsActive(unsigned ct)const{ return(int(ct)<BsNum && BsActive[ct]); }
  unsigned GetBs(unsigned ct)const{ return(unsigned(BsMin+BsInc*int(ct))); }

  void SetTime(unsigned ct,float t){ if(int(ct)<BsNum)Times[ct]=t; }
  void ProcessTimes(double timestep,unsigned nstep);
  unsigned GetOptimumBs()const{ return(BsSel); }
};


//##############################################################################
//# JBlockSizeAuto
//##############################################################################
/// \brief Manages the automatic computation of optimum Blocksize in kernel interactions.

class JBlockSizeAuto : protected JObject
{
protected:
  JLog2 *Log;
  unsigned StepsInterval;
  std::vector<JBlockSizeAutoKer*> Kernels;

  unsigned GetKernelByName(std::string name)const;

public:
  JBlockSizeAuto(JLog2 *log,unsigned steps);
  ~JBlockSizeAuto();
  void Reset();
  
  unsigned GetStepsInterval()const{ return(StepsInterval); }
  unsigned GetKernelsCount()const{ return(unsigned(Kernels.size())); }

  void AddKernel(std::string name,int bsmin,int bsnum,int bsinc,int bsdefault);

  JBlockSizeAutoKer* GetKernel(unsigned c){ return(c<GetKernelsCount()? Kernels[c]: NULL); }

  void ProcessTimes(double timestep,unsigned nstep);


};


#endif


