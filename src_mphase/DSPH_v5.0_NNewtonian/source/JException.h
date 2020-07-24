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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - El constructor muestra mensaje por pantalla de forma automatica. (14-09-2019)
//:# - Gestion de excepciones mejorada.  (15-09-2019)
//:#############################################################################

/// \file JException.h \brief Declares the class \ref JException.

#ifndef _JException_
#define _JException_

#include <string>

//##############################################################################
//# JException
//##############################################################################
/// \brief Defines exceptions with the information of the class and method.

class JException : public std::exception{
 protected:
  std::string ExName;    ///<Name of the exception. 
  std::string SrcFile;   ///<Name of the file that generated an exception. 
  int SrcLine;           ///<Number of the line that generated an exception. 
  std::string ClassName; ///<Name of the class that generated an exception. 
  std::string Method;    ///<Name of the method that generated an exception. 
  std::string Text;      ///<Text of the exception.
  std::string File;      ///<File related to the exception.
 public:
  JException(const std::string &classname,const std::string &method,const std::string &text,const std::string &file);
  JException(const std::string &srcfile,int srcline,const std::string &classname,const std::string &method,const std::string &text,const std::string &file);
  ~JException() throw(){}  ///<Destructor of objects.
  std::string ToStr()const;
  void Print()const;
  virtual const char* what() const throw(){ 
    static char tx[2048];
    std::string tex=ToStr();
    for(unsigned c=0;c<=tex.size() && c<2047;c++)tx[c]=tex[c]; tx[2047]='\0';
    return(tx);
  } 
};

#endif


