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

/// \file JException.cpp \brief Implements the class \ref JException.

#include "JException.h"
#include "Functions.h"
#include <cstdio>

//==============================================================================
/// Constructor.
/// \param classname Name of the class that throws the exception.
/// \param method Name of the method that throws the exception.
/// \param text Text of the exception.
/// \param file Name of the file related to the exception.
//==============================================================================
JException::JException(const std::string &classname,const std::string &method
  ,const std::string &text,const std::string &file)
{
  ExName="JException";
  SrcFile="";
  SrcLine=0;
  ClassName=classname;
  Method=method;
  Text=text;
  File=file;
  Print();
}

//==============================================================================
/// Constructor.
/// \param srcfile Name of the file that throws the exception.
/// \param srcline Number of the line that throws the exception.
/// \param classname Name of the class that throws the exception.
/// \param method Name of the method that throws the exception.
/// \param text Text of the exception.
/// \param file Name of the file related to the exception.
//==============================================================================
JException::JException(const std::string &srcfile,int srcline
  ,const std::string &classname,const std::string &method
  ,const std::string &text,const std::string &file)
{
  ExName="JException";
  SrcFile=srcfile;
  SrcLine=srcline;
  ClassName=classname;
  Method=method;
  Text=text;
  File=file;
  Print();
}

//==============================================================================
/// Returns the complete text message with the information of the exception. 
//==============================================================================
std::string JException::ToStr()const{
  std::string tx;
  if(SrcFile.empty())tx=fun::PrintStr("\n*** Exception (%s::%s)\n",ClassName.c_str(),Method.c_str());
  else tx=fun::PrintStr("\n*** Exception (%s::%s) at %s:%d\n",ClassName.c_str(),Method.c_str(),SrcFile.c_str(),SrcLine);
  if(!Text.empty())tx=tx+fun::PrintStr("Text: %s\n",Text.c_str());
  if(!File.empty())tx=tx+fun::PrintStr("File: %s\n",File.c_str());
  return(tx);
}

//==============================================================================
/// Visualises the exception message in console.
//==============================================================================
void JException::Print()const{
  printf("%s\n",ToStr().c_str());
  fflush(stdout);
}


