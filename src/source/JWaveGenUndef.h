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

/// \file JWaveGenUndef.h \brief Declares the empty class \ref JWaveGen.

#ifndef _JWaveGenUndef_
#define _JWaveGenUndef_

//##############################################################################
//# JWaveGen
//##############################################################################
/// \brief Implements wave generation for regular and irregular waves.

#ifdef DISABLE_WAVEGEN
class JWaveGen
{
public:
  StMotionData MotionData;
  JWaveGen(bool useomp,bool usegpu,JLog2 *log,std::string dirdata,const JXml *sxml
    ,const std::string &place,tdouble3 gravity3){}
  ~JWaveGen(){}
  static bool Available(){ return(false); }
  bool ConfigPaddleParts(word mkbound,word motionref,unsigned idbegin,unsigned np){ return(false); }
  void SetTimeMod(double timemod){}
  word GetPaddleMkbound(unsigned cp)const{ return(0); }
  void InitPaddle(unsigned cp,double timemax,double timepart
    ,const StWvgDimensions &wdims){}
  bool PaddleUseAwas(unsigned cp)const{ return(false); }
  void PaddleGetAwasInfo(unsigned cp,double coefmassdef,double massfluid
    ,double &masslimit,double &tstart,double &gdp,tdouble3 &point0,tdouble3 &point2)const{}
  void PaddleGaugeInit(unsigned cp,void *gswl){}
  void* PaddleGaugeGet(unsigned cp)const{ return(NULL); }
  void VisuConfig(std::string txhead,std::string txfoot){}

  void LoadLastGaugeResults(unsigned cp,double timestep,const tfloat3 &posswl,const tfloat3 &point0){}
  void GetUpdateGaugePoints(unsigned cp,bool updatepoints,tdouble3 &point0,tdouble3 &point2)const{}
  
  const StMotionData& GetMotionNull()const{ return(MotionData); };
  const StMotionData& GetMotion(bool svdata,unsigned cp,double timestep,double dt){ return(MotionData); }
  const StMotionData& GetMotionAce(bool svdata,unsigned cp,double timestep,double dt){ return(MotionData); }

  unsigned GetCount()const{ return(0); }
  bool UseAwas()const{ return(false); }
  bool WavesRegular() const{ return(false); }
  bool WavesSpectrum()const{ return(false); }
  bool WavesFocused() const{ return(false); }
  bool WavesFile()    const{ return(false); }
  bool WavesSolitary()const{ return(false); }
};
#endif

#endif


