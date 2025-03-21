//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JMotionData.h \brief Declares the classes \ref JMotionData.

#ifndef _JMotionData_
#define _JMotionData_

#include "TypesDef.h"
#include "JObject.h"

class JXml;
class TiXmlNode;

//##############################################################################
//# JMotionDataMov
//##############################################################################
/// \brief Stores data for rectilinear movement in time.
class JMotionDataMov : protected JObject
{
private:
  static const unsigned SIZEMAX=104857600; ///<Maximum file size (100mb).

  std::vector<double>   Times;   //-Times.
  std::vector<tdouble3> DataPos; //-Position in time.

  void Reset();
  void LoadFilePos(std::string dirdata,std::string file,const int fields
    ,const int fieldtime,const int fieldx,const int fieldy,const int fieldz);

public:
  JMotionDataMov(std::string dirdata,std::string file,int fields,int fieldtime
    ,int fieldx,int fieldy,int fieldz);
  ~JMotionDataMov();

  unsigned GetCount()const{ return(unsigned(Times.size())); }
  const double*   GetTimes()const{ return(Times.data()); }
  const tdouble3* GetValuesPos()const{ return(DataPos.data()); }
};

//##############################################################################
//# JMotionDataRotAxis
//##############################################################################
/// \brief Stores data for rotation according to an axis in time.
class JMotionDataRotAxis : protected JObject
{
private:
  static const unsigned SIZEMAX=104857600; ///<Maximum file size (100mb).

  std::vector<double> Times;   //-Times.
  std::vector<double> DataAng; //-Angle in time (degrees).

  void Reset();
  void LoadFileAng(std::string dirdata,std::string file,bool angdegrees);

public:
  JMotionDataRotAxis(std::string dirdata,std::string file,bool angdegrees);
  ~JMotionDataRotAxis();

  unsigned GetCount()const{ return(unsigned(Times.size())); }
  const double* GetTimes()const{ return(Times.data()); }
  const double* GetValuesAng()const{ return(DataAng.data()); }
};

//##############################################################################
//# JMotionDataRotEuler
//##############################################################################
/// \brief Stores data for Euler rotation in time (3 angles).
class JMotionDataRotEuler : protected JObject
{
private:
  static const unsigned SIZEMAX=104857600; ///<Maximum file size (100mb).

  std::vector<double>   Times;   //-Times.
  std::vector<tdouble3> DataAng; //-Angles in time (degrees).

  void Reset();
  void LoadFile3Ang(std::string dirdata,std::string file,const int fields
    ,const int fieldtime,const int fieldang1,const int fieldang2
    ,const int fieldang3,bool angdegrees);

public:
  JMotionDataRotEuler(std::string dirdata,std::string file,int fields
    ,int fieldtime,int fieldang1,int fieldang2,int fieldang3,bool angdegrees);
  ~JMotionDataRotEuler();

  unsigned GetCount()const{ return(unsigned(Times.size())); }
  const double* GetTimes()const{ return(Times.data()); }
  const tdouble3* GetValuesAng3()const{ return(DataAng.data()); }
};

//##############################################################################
//# JMotionDataPath
//##############################################################################
/// \brief Stores data for path movement in time (displacement and Euler rotation).
class JMotionDataPath : protected JObject
{
private:
  static const unsigned SIZEMAX=104857600; ///<Maximum file size (100mb).

  std::vector<double>   Times;   //-Times.
  std::vector<tdouble3> DataPos; //-Position in time.
  std::vector<tdouble3> DataAng; //-Angles in time (degrees).

  void Reset();
  void LoadPathFile(std::string dirdata,std::string file,const int fields
    ,const int fieldtime,const int fieldx,const int fieldy,const int fieldz
    ,const int fieldang1,const int fieldang2,const int fieldang3,bool angdegrees);

public:
  JMotionDataPath(std::string dirdata,std::string file,int fields
    ,int fieldtime,int fieldx,int fieldy,int fieldz,int fieldang1,int fieldang2
    ,int fieldang3,bool angdegrees);
  ~JMotionDataPath();

  unsigned GetCount()const{ return(unsigned(Times.size())); }
  const double* GetTimes()const{ return(Times.data()); }
  const tdouble3* GetValuesPos()const{ return(DataPos.data()); }
  const tdouble3* GetValuesAng3()const{ return(DataAng.data()); }
};

#endif

