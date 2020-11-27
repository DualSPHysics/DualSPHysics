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
//:# Clase para facilitar la lectura de ficheros de datos en texto. (03-11-2015)
//:# - Los valores pueden estar delimitados por ',' ';' o tabs-spaces.
//:# - El delimitador se identifica automaticamente.
//:# - Cuenta numero de lineas de datos y de comentarios.
//:# - Las lineas de comentarios empiezan por '#'.
//:# - Permite extraer valores int, unsigned, float, double y string.
//:# - Mantiene columna y fila del ultimo valor.
//:# - En caso de usar el separador \t los espacios y tabulaciones seguidas se
//:#   remplazan por una sola tabulacion. Tambien elimina espacions y tabulaciones
//:#   al principio y final de la linea. (16-12-2015)
//:# - Ignora lineas vacias al final del fichero. (16-12-2015)
//:# - Error corregido en ProcessSpaces(). (17-12-2015)
//:# - Nuevos metodos: ReadNextInt3() y ReadNextUnsigned3(). (27-04-2016)
//:# - Nuevos metodo: RemoveChar(). (20-06-2016)
//:# - Nuevos metodo: Find(),FindValueStr(),FindValueDbl(). (02-07-2018)
//:# - Muestra extension maxima de fichero permitida en caso de ser insuficiente. (21-10-2019)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Nuevo metodo: ReadNextBool(). (19-08-2020)
//:#############################################################################

/// \file JReadDatafile.h \brief Declares the class \ref JReadDatafile.

#ifndef _JReadDatafile_
#define _JReadDatafile_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"

//##############################################################################
//# JReadDatafile
//##############################################################################
/// \brief Allows reading data in ASCII files.
// Clase para facilitar la lectura de ficheros de datos en texto.

class JReadDatafile  : protected JObject
{
private:
  std::string File;      ///< Name of file.

  unsigned SizeFile;     ///< Size of file.
  unsigned Size;         ///< Size of data.
  char *Data;            ///< Data from file.

  int LineCount;         ///< Number of lines.
  unsigned *LineBegin;   ///< Inicio de cada linea [LineCount+1].
  int RemLineCount;      ///< Number of remark lines.

  std::string Sep;       ///< Value separator.

  int ReadLin;
  int ReadLinValue;
  std::string ReadLine;
  std::string ReadValue;

  void ProcessSpaces();
  void ProcessLines();
  void ResetReadLine();

public:
  JReadDatafile();
  ~JReadDatafile();
  void Reset();

  void LoadFile(const std::string &file,unsigned maxsize=1048576000);
  void RemoveChar(char let);

  unsigned Lines()const{ return(LineCount); }
  unsigned RemLines()const{ return(RemLineCount); }

  void SetReadLine(int line);

  std::string GetLine(int line)const;
  
  std::string ReadNextValue(bool in_line=false);

  bool ReadNextBool(bool in_line=false);

  double ReadNextDouble(bool in_line=false);
  tdouble3 ReadNextDouble3(bool in_line=false);
  
  float ReadNextFloat(bool in_line=false){    return(float(ReadNextDouble(in_line)));      }
  tfloat3 ReadNextFloat3(bool in_line=false){ return(ToTFloat3(ReadNextDouble3(in_line))); }

  int ReadNextInt(bool in_line=false);
  tint3 ReadNextInt3(bool in_line=false);

  unsigned ReadNextUnsigned(bool in_line=false){ return(unsigned(ReadNextInt(in_line))); }
  tuint3 ReadNextUnsigned3(bool in_line=false);

  int GetReadLin()const{           return(ReadLin);      }
  int GetReadLinValue()const{      return(ReadLinValue); }
  std::string GetReadValue()const{ return(ReadValue);    }


  tint2 Find(std::string key,int firstline=0)const;
  std::string FindValueStr(std::string key,bool optional=false,std::string valdef="")const;
  double FindValueDbl(std::string key,bool optional=false,double valdef=0)const;

};

#endif


