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

/// \file JDsphConfig.cpp \brief Implements the class \ref JDsphConfig.

#include "JDsphConfig.h"
#include "Functions.h"
#include "JXml.h"
#include <algorithm>

using namespace std;

//##############################################################################
//# JDsphConfig
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsphConfig::JDsphConfig(){
  ClassName="JDsphConfig";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsphConfig::~JDsphConfig(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsphConfig::Reset(){
  FileCfg="";
  CreateDirs=-1;
  CsvSeparator=-1;
}

//==============================================================================
/// Initialisation of log file.
//==============================================================================
void JDsphConfig::Init(std::string path){
  Reset();
  const string file=fun::GetDirWithSlash(path)+"DsphConfig.xml";
  if(fun::FileExists(file)){
    //-Loads XML file.
    JXml sxml;
    sxml.LoadFile(file);
    const string place="dsphconfig.common";
    TiXmlNode* node=sxml.GetNode(place,false);
    if(!node)Run_ExceptioonFile(string("Cannot find the element \'")+place+"\'.",file);
    //-Reads configuration values in XML file.
    CreateDirs=sxml.ReadElementInt(node,"createdirs","v",true,-1);
    CsvSeparator=sxml.ReadElementInt(node,"csvseparator","v",true,-1);
    //printf("JDsphConfig::CreateDirs: %d\n",CreateDirs);
    //printf("JDsphConfig::csvseparator: %d\n",CsvSeparator);
    //-Stores name of configuration file.
    FileCfg=file;
  }
}


