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

//#############################################################################
//# Cambios:
//# =========
//# - Nuevos metodos LoadFloat3() y LoadDouble3(). (27-06-2020)
//# - Actualiza codigo para DualSPHysics. (28-06-2020)
//# - Nuevas funciones LoadFloats() y LoadDoubles(). (12-08-2020)
//# - Ignora parametros vacios. (10-09-2020)
//:# - El uso de JDsphConfig o no se define en JCfgRunBaseDef.h. (20-04-2021)
//#############################################################################

#ifndef _JCfgRunBase_
#define _JCfgRunBase_

#pragma warning(disable : 4996) //Anula sprintf() deprecated.

#include "JCfgRunBaseDef.h"
#include "TypesDef.h"
#include "Functions.h"
#include "JObject.h"
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>

//==============================================================================
//##############################################################################
//==============================================================================
class JCfgRunBase : protected JObject
{
protected:
  int ParmDef;
  void LoadDsphConfig(std::string path);

  static unsigned LoadFloats (std::string txopt,float  def,unsigned nv,std::vector<float > &vv);
  static unsigned LoadDoubles(std::string txopt,double def,unsigned nv,std::vector<double> &vv);
  static unsigned LoadFloat3 (std::string txopt,float def,unsigned nv,tfloat3 *v);
  static unsigned LoadDouble3(std::string txopt,double def,unsigned nv,tdouble3 *v);
  static void LoadFloat3(std::string txopt,float def,tfloat3 &v1);
  static void LoadFloat6(std::string txopt,float def,tfloat3 &v1,tfloat3 &v2);
  static void LoadDouble3(std::string txopt,double def,tdouble3 &v1);
  static void LoadDouble6(std::string txopt,double def,tdouble3 &v1,tdouble3 &v2);

  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull,std::string &txopt1,std::string &txopt2,std::string &txopt3,std::string &txopt4)const;
  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull,std::string &txopt1,std::string &txopt2,std::string &txopt3)const{
    std::string tx4; SplitsOpts(opt,txword,txoptfull,txopt1,txopt2,txopt3,tx4);
  }
  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull,std::string &txopt1,std::string &txopt2)const{
    std::string tx3,tx4; SplitsOpts(opt,txword,txoptfull,txopt1,txopt2,tx3,tx4);
  }
  void SplitsOpts(const std::string &opt,std::string &txword,std::string &txoptfull)const{
    std::string tx1,tx2,tx3,tx4; SplitsOpts(opt,txword,txoptfull,tx1,tx2,tx3,tx4);
  }

public:
  const bool NoParms; ///<Allows zero parameters without showing help.
  bool PrintInfo;

  //-General configuration from DsphConfig.xml
  bool CreateDirs;   ///<Creates full path for output files (true by default).
  bool CsvSepComa;   ///<Separator character in CSV files (false=semicolon, true=coma).

public:
  JCfgRunBase(bool noparms=false);
  void Reset();

  void LoadArgv(int argc,char** argv);
  void LoadFile(std::string fname,int lv);
  void ErrorParm(const std::string &opt,int optc,int lv,const std::string &file)const;
  void ErrorParmText(const std::string &text,int optc,int lv,const std::string &file)const;

  virtual void VisuInfo()const=0;
  virtual void VisuConfig()const=0;
  virtual void LoadOpts(std::string *optlis,int optn,int lv,const std::string &file)=0;
  virtual void ValidaCfg()=0;
};

#endif

