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

/// \file JMotionData.cpp \brief Implements the classes \ref JMotionData.

#include "JMotionData.h"
#include "JReadDatafile.h"

//##############################################################################
//# JMotionDataMov
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMotionDataMov::JMotionDataMov(std::string dirdata,std::string file,int fields
  ,int fieldtime,int fieldx,int fieldy,int fieldz)
{
  ClassName="JMotionDataMov";
  Reset();
  LoadFilePos(dirdata,file,fields,fieldtime,fieldx,fieldy,fieldz);
}

//==============================================================================
/// Destructor.
//==============================================================================
JMotionDataMov::~JMotionDataMov(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JMotionDataMov::Reset(){
  Times.clear();
  DataPos.clear();
}

//==============================================================================
/// Loads times and motion data from file.
//==============================================================================
void JMotionDataMov::LoadFilePos(std::string dirdata,std::string file
  ,const int fields,const int fieldtime,const int fieldx,const int fieldy
  ,const int fieldz)
{
  //-Adds dirdata when only the filename is provided.
  if(int(file.find("/"))<0 && int(file.find("\\"))<0)file=dirdata+file;
  //-Loads file with data.
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  //printf("----->rows:%u\n",rows);
  //-Allocates memory.
  Times.reserve(rows);
  DataPos.reserve(rows);
  //-Loads data arrays.
  for(unsigned r=0;r<rows;r++){
    double time,posx=0,posy=0,posz=0;
    for(int f=0;f<fields;f++){
      const double v=rdat.ReadNextDouble(f!=0);
      //printf("[%u]------>  v[%u]:%f\n",r,f,v);
      if(f==fieldtime)time=v;
      if(f==fieldx)posx=v;
      if(f==fieldy)posy=v;
      if(f==fieldz)posz=v;
    }
    Times.push_back(time);
    DataPos.push_back(TDouble3(posx,posy,posz));
    //printf("[%u]>  t:%f  x:%f\n",r,time,posx);
  }
  //-Checks minimum data.
  if(GetCount()<2)Run_ExceptioonFile("Data is required for at least two times.",file);
}


//##############################################################################
//# JMotionDataRotAxis
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMotionDataRotAxis::JMotionDataRotAxis(std::string dirdata,std::string file
  ,bool angdegrees)
{
  ClassName="JMotionDataRotAxis";
  Reset();
  LoadFileAng(dirdata,file,angdegrees);
}

//==============================================================================
/// Destructor.
//==============================================================================
JMotionDataRotAxis::~JMotionDataRotAxis(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JMotionDataRotAxis::Reset(){
  Times.clear();
  DataAng.clear();
}

//==============================================================================
/// Loads times and motion data from file.
//==============================================================================
void JMotionDataRotAxis::LoadFileAng(std::string dirdata,std::string file
  ,bool angdegrees)
{
  //-Adds dirdata when only the filename is provided.
  if(int(file.find("/"))<0 && int(file.find("\\"))<0)file=dirdata+file;
  //-Loads file with data.
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  //-Allocates memory.
  Times.reserve(rows);
  DataAng.reserve(rows);
  //-Loads data arrays.
  for(unsigned r=0;r<rows;r++){
    Times.push_back(rdat.ReadNextDouble());
    DataAng.push_back(rdat.ReadNextDouble());
  }
  //-Changes angles from radians to degrees.
  if(!angdegrees)for(unsigned r=0;r<rows;r++)DataAng[r]=DataAng[r]*TODEG;
  //-Checks minimum data.
  if(GetCount()<2)Run_ExceptioonFile("Data is required for at least two times.",file);
}


//##############################################################################
//# JMotionDataRotEuler
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMotionDataRotEuler::JMotionDataRotEuler(std::string dirdata,std::string file
  ,int fields,int fieldtime,int fieldang1,int fieldang2,int fieldang3
  ,bool angdegrees)
{
  ClassName="JMotionDataRotEuler";
  Reset();
  LoadFile3Ang(dirdata,file,fields,fieldtime,fieldang1,fieldang2,fieldang3
    ,angdegrees);
}

//==============================================================================
/// Destructor.
//==============================================================================
JMotionDataRotEuler::~JMotionDataRotEuler(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JMotionDataRotEuler::Reset(){
  Times.clear();
  DataAng.clear();
}

//==============================================================================
/// Loads times and motion data from file.
//==============================================================================
void JMotionDataRotEuler::LoadFile3Ang(std::string dirdata,std::string file
  ,const int fields,const int fieldtime,const int fieldang1,const int fieldang2
  ,const int fieldang3,bool angdegrees)
{
  //-Adds dirdata when only the filename is provided.
  if(int(file.find("/"))<0 && int(file.find("\\"))<0)file=dirdata+file;
  //-Loads file with data.
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  //-Allocates memory.
  Times.reserve(rows);
  DataAng.reserve(rows);
  //-Loads data arrays.
  for(unsigned r=0;r<rows;r++){
    double time,ang1=0,ang2=0,ang3=0;
    for(int f=0;f<fields;f++){
      const double v=rdat.ReadNextDouble(f!=0);
      if(f==fieldtime)time=v;
      if(f==fieldang1)ang1=v;
      if(f==fieldang2)ang2=v;
      if(f==fieldang3)ang3=v;
    }
    //printf("[%u]>  t:%f  x:%f y:%f z:%f\n",r,time,ang1,ang2,ang3);
    Times.push_back(time);
    DataAng.push_back(TDouble3(ang1,ang2,ang3));
  }
  //-Changes angles from radians to degrees.
  if(!angdegrees)for(unsigned r=0;r<rows;r++)DataAng[r]=DataAng[r]*TODEG;
  //-Checks minimum data.
  if(GetCount()<2)Run_ExceptioonFile("Data is required for at least two times.",file);
}


//##############################################################################
//# JMotionDataRotEuler
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMotionDataPath::JMotionDataPath(std::string dirdata,std::string file
  ,int fields,int fieldtime,int fieldx,int fieldy,int fieldz,int fieldang1
  ,int fieldang2,int fieldang3,bool angdegrees)
{
  ClassName="JMotionDataPath";
  Reset();
  LoadPathFile(dirdata,file,fields,fieldtime,fieldx,fieldy,fieldz,fieldang1
    ,fieldang2,fieldang3,angdegrees);
}

//==============================================================================
/// Destructor.
//==============================================================================
JMotionDataPath::~JMotionDataPath(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JMotionDataPath::Reset(){
  Times.clear();
  DataPos.clear();
  DataAng.clear();
}

//==============================================================================
/// Loads times and motion data from file.
//==============================================================================
void JMotionDataPath::LoadPathFile(std::string dirdata,std::string file
  ,const int fields,const int fieldtime,const int fieldx,const int fieldy
  ,const int fieldz,const int fieldang1,const int fieldang2,const int fieldang3
  ,bool angdegrees)
{
  //-Adds dirdata when only the filename is provided.
  if(int(file.find("/"))<0 && int(file.find("\\"))<0)file=dirdata+file;
  //-Loads file with data.
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  //-Allocates memory.
  Times.reserve(rows);
  DataPos.reserve(rows);
  DataAng.reserve(rows);
  //-Loads data arrays.
  for(unsigned r=0;r<rows;r++){
    double time,posx=0,posy=0,posz=0,ang1=0,ang2=0,ang3=0;
    for(int f=0;f<fields;f++){
      const double v=rdat.ReadNextDouble(f!=0);
      if(f==fieldtime)time=v;
      if(f==fieldx)posx=v;
      if(f==fieldy)posy=v;
      if(f==fieldz)posz=v;
      if(f==fieldang1)ang1=v;
      if(f==fieldang2)ang2=v;
      if(f==fieldang3)ang3=v;
    }
    Times.push_back(time);
    DataPos.push_back(TDouble3(posx,posy,posz));
    DataAng.push_back(TDouble3(ang1,ang2,ang3));
  }
  //-Changes angles from radians to degrees.
  if(!angdegrees)for(unsigned r=0;r<rows;r++)DataAng[r]=DataAng[r]*TODEG;
  //-Checks minimum data.
  if(GetCount()<2)Run_ExceptioonFile("Data is required for at least two times.",file);
}

