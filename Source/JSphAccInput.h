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

/// \file JSphAccInput.h \brief Declares the class \ref JSphAccInput.

#ifndef _JSphAccInput_
#define _JSphAccInput_

#include "JObject.h"
#include "Types.h"
#include "JTimer.h"
#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>

class JLog2;
class JXml;
class TiXmlElement;

//##############################################################################
//# XML format in _FmtXML_AccInput.xml.
//##############################################################################

//##############################################################################
//# JSphAccInputMk
//##############################################################################
/// \brief Provides the force to be applied to different blocks of particles that is loaded from files.

class JSphAccInputMk : protected JObject
{
protected:
  JLog2* Log;

  static const unsigned SIZEMAX=104857600; ///<Maximum file size (100mb).
  static const unsigned SIZEINITIAL=100;

  word MkFluid;              ///<The MK values stored in the acceleration input file.
  bool GravityEnabled;       ///<Determines whether global gravity is enabled or disabled for this particle set SL
  tfloat3 AccCoG;            ///<The centre of gravity that will be used for angular acceleration calculations.
  std::string File;          ///<File of data.

  unsigned AccSize;          ///<Number of acceleration values that were allocated.
  unsigned AccCount;         ///<Number of acceleration values in each input file(s).
  float *AccTime;            ///<Variable acceleration time evolution as detailed in the input file.
  tfloat3 *AccLin;           ///<Linear acceleration variable to store values as they are read from the input files.
  tfloat3 *AccAng;           ///<Angular acceleration variable to store values as they are read from the input files.
  tfloat3 *VelAng;           ///<Angular velocity variable to store values as the angular acceleration values are read from the input files. SL
  tfloat3 *VelLin;           ///<Linear velocity variable to store values as the linear acceleration values are read from the input files. SL

  unsigned AccIndex;         ///<Current index for variable acceleration interpolation.

  tdouble3 CurrAccLin;        ///<The current interpolated values for linear acceleration.
  tdouble3 CurrAccAng;        ///<The current interpolated values for angular acceleration.
  tdouble3 CurrVelLin;        ///<The current interpolated values for linear velocity. SL
  tdouble3 CurrVelAng;        ///<The current interpolated values for angular velocity. SL


  void Reset();
  void Resize(unsigned size);
  void LoadFile(std::string file,double tmax);

public:
  JSphAccInputMk(JLog2* log,word mkfluid,bool genabled,tfloat3 acccentre,std::string file);
  ~JSphAccInputMk();
  long long GetAllocMemory()const;

  void Init(double tmax);
  void GetConfig(std::vector<std::string> &lines)const;

  word GetMkFluid()const{ return(MkFluid); }
  void GetAccValues(double timestep,unsigned &mkfluid,tdouble3 &acclin,tdouble3 &accang,tdouble3 &centre,tdouble3 &velang,tdouble3 &vellin,bool &setgravity); //SL: Added linear and angular velocity and set gravity flag
};

//##############################################################################
//# JSphAccInput
//##############################################################################
/// \brief Manages the application of external forces to different blocks of particles (with the same MK).

class JSphAccInput : protected JObject
{
protected:
  JLog2* Log;
  std::string DirData;
  std::vector<JSphAccInputMk*> Inputs;
  long long MemSize;

  void Reset();
  bool ExistMk(word mkfluid)const;
  void LoadXml(JXml *sxml,const std::string &place);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);

public:
  JSphAccInput(JLog2* log,const std::string &dirdata,JXml *sxml,const std::string &place);
  ~JSphAccInput();
  long long GetAllocMemory()const{ return(MemSize); }

  void Init(double tmax);
  void VisuConfig(std::string txhead,std::string txfoot)const;

  unsigned GetCount()const{ return(unsigned(Inputs.size())); };
  void GetAccValues(unsigned cfile,double timestep,unsigned &mkfluid,tdouble3 &acclin,tdouble3 &accang,tdouble3 &centre,tdouble3 &velang,tdouble3 &vellin,bool &setgravity); //SL: Added linear and angular velocity and set gravity flag
};

#endif


