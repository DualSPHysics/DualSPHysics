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

/// \file JChronoGeoObj.h \brief Declares the class \ref JChronoGeoObj.

#ifndef _JChronoGeoObj_
#define _JChronoGeoObj_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"

class JBinaryData;

//##############################################################################
//# JChronoGeoObj
//##############################################################################
/// \brief Creates OBJ files for Chrono collisions.

class JChronoGeoObj : protected JObject
{
private:
  static const unsigned FmtVersionDef=241130;   ///<Version de formato by default. Version of format by default.
  unsigned FmtVersion;    ///<Version de formato. Version of format.

  std::string CaseName;
  std::string DirData;
  std::string FileData;
  word MkBoundFirst;

  JBinaryData* Bdat;

private:
  void SavesFileObj(std::string fname,byte normalmode,unsigned npos
    ,const tfloat3* pos,unsigned ntri,const tuint3* vtri
    ,unsigned nqua,const tuint4* vqua)const;

private:
  JChronoGeoObj();
  ~JChronoGeoObj();
  void Reset();
  void LoadFile(std::string casename);

  unsigned CreateMkObj(std::string filein,std::string filesout
    ,const std::vector<unsigned>& mkbounds,byte normalmode);

public:
  static unsigned CreateOBJsByMk(std::string casename
    ,std::string filein,std::string filesout
    ,const std::vector<unsigned>& mkbounds,byte normalmode);

};

#endif

