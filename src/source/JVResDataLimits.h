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

/// \file JVResDataLimits.h \brief Declares the class \ref JVResDataLimits.

#ifndef _JVResDataLimits_
#define _JVResDataLimits_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"
#include "JBoxDef.h"

class JCaseVRes;
class JBinaryData;

//##############################################################################
//# JVResDataLimitsZone
//##############################################################################
/// \brief Limtis of zone and subzones.
class JVResDataLimitsZone : protected JObject
{
public:
  typedef struct{
    unsigned Id;
    JBoxDef InnLimit;
    JBoxDef InnLimit2h;
  }StSubZone;

public:
  const unsigned Id;
  const unsigned ParentId;   ///<Parent id or UINT_MAX.
  const bool OutLimits;      ///<Outter limits are defined (ParentId!=UINT_MAX).
  const int TrackingMk;      ///<Mk for tracking (-1 by default).
  const int TrackingMkBound; ///<MkBound for tracking (-1 by default).

protected:
  //-Outter limits when it is a subdomain.
  JBoxDef OutLimit;
  JBoxDef OutLimit2h;

  //-Inner limits when it has subdomains.
  std::vector<StSubZone> SubZones;

public:
  JVResDataLimitsZone(unsigned id,unsigned parentid,int trackmk,int trackmkb);
  ~JVResDataLimitsZone();
  void Reset();

  void SetOutLimits(const JBoxDef& outlimit,const JBoxDef& outlimit2h);
  void AddInnLimits(unsigned id,const JBoxDef& innlimit,const JBoxDef& innlimit2h);

  JBoxDef GetOutLimit()const{ return(OutLimit); }
  JBoxDef GetOutLimit2h()const{ return(OutLimit2h); }

  unsigned SubCount()const{ return(unsigned(SubZones.size())); }
  unsigned GetInnId(unsigned csub)const{ return(SubZones[csub].Id); }
  JBoxDef  GetInnLimit  (unsigned csub)const{ return(SubZones[csub].InnLimit  ); }
  JBoxDef  GetInnLimit2h(unsigned csub)const{ return(SubZones[csub].InnLimit2h); }
};

//##############################################################################
//# JVResDataLimits
//##############################################################################
/// \brief Manages the information of VRes zones.
class JVResDataLimits : protected JObject
{
private:
  static const unsigned FmtVersionDef=240604; ///<Version de formato by default. Version of format by default.
  unsigned FmtVersion;    ///<Version de formato. Version of format.

  bool   Is2D;            ///<2-D domain.
  double Posy2D;          ///<Y-position for 2-D domain.
  std::string DirData;
  std::vector<JVResDataLimitsZone*> Zones;

private:
  JVResDataLimitsZone* AddZone(unsigned id,unsigned parentid,int trackmk,int trackmkb);
  static void GetBoxPoints2d(JBoxDef box,tdouble3* vpt);
  static void GetBoxPoints3d(JBoxDef box,tdouble3* ptvec);

  void    SaveBoxDef(std::string name,std::string subname,const JBoxDef& box,JBinaryData* item)const;
  JBoxDef LoadBoxDef(std::string name,std::string subname,const JBinaryData* item)const;
  void LoadFileFmt240102(JBinaryData* bdat)const;
public:
  JVResDataLimits();
  ~JVResDataLimits();
  void Reset();

  static std::string GetFileName(std::string dir="");

  void LoadData(const JCaseVRes* vrdat);
  void SaveFile(std::string dir);
  void LoadFile(std::string dir);

  unsigned Count()const{ return(unsigned(Zones.size())); }
  const JVResDataLimitsZone* GetZone(unsigned id)const;

  void SaveVtkLimits(std::string fname,bool onefile)const;
};

#endif

