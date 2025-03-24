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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Class to loads basic Variable Resolution configuration. (16-06-2022)
//:# - Avoids VTK shapes of fixed limits. (16-05-2024)
//:# - New general-purpose section <extra>. (31-05-2024)
//:# - Use JSpVtkShape instead of JVtk Lib. (10-12-2024)
//:# - Removes SaveVtkPoints(). (10-12-2024)
//:#############################################################################

/// \file JCaseVRes.h \brief Declares the class \ref JCaseVRes.

#ifndef _JCaseVRes_
#define _JCaseVRes_

#include "JObject.h"
#include "JBoxDef.h"
#include "TypesDef.h"
#include <string>
#include <vector>
#include <climits>

class JXml;
class TiXmlElement;
class JCaseVResBase;
class JCaseVRes_Box;

///Shapes of buffer zones.
typedef enum{ 
  MRSH_Null=0
 ,MRSH_Box
}TpMRShape; 

//##############################################################################
//# JCaseVResExtra
//##############################################################################
/// \brief class for general-purpose section.
class JCaseVResExtra : protected JObject 
{
protected:
  std::vector<std::string> Keys;
  std::vector<std::string> Values;

public:
  JCaseVResExtra();
  JCaseVResExtra(const JCaseVResExtra& src);
  ~JCaseVResExtra();
  JCaseVResExtra& operator=(const JCaseVResExtra& src);
  void Reset();

  void AddValue(const std::string& key,const std::string& value);

  unsigned Count()const{ return(unsigned(Keys.size())); }
  std::string GetKey(unsigned idx)const;
  std::string GetValue(unsigned idx)const;

  unsigned FindKey(const std::string& key)const;
  bool ExistsValue(const std::string& key)const{
    return(FindKey(key)!=UINT_MAX);
  }

  std::string GetValueStr (const std::string& key,bool optional=false
    ,std::string valdef="")const;
  int         GetValueInt (const std::string& key,bool optional=false
    ,int valdef=0)const;
  unsigned    GetValueUint(const std::string& key,bool optional=false
    ,unsigned valdef=0)const;
  double      GetValueDbl (const std::string& key,bool optional=false
    ,double valdef=0)const;
  tdouble3    GetValueDbl3(const std::string& key,bool optional=false
    ,tdouble3 valdef=TDouble3(0))const{
      return(TDouble3(GetValueDbl(key+".x",optional,valdef.x),
                      GetValueDbl(key+".y",optional,valdef.y),
                      GetValueDbl(key+".z",optional,valdef.z)));
  }

  void Print(std::string caption="")const;
};

//##############################################################################
//# JCaseVResBase
//##############################################################################
/// \brief Base class for buffer zone configuration.
class JCaseVResBase : protected JObject 
{
public:
  const TpMRShape Shape;
  const int Id;
  JCaseVResBase* const Parent;
  const double Dp;
  const bool   Is2D;    ///< 2-D domain.
  const double Posy2D;  ///< Y-position for 2-D domain.

  std::vector<JCaseVResBase*> SubZones;

  //-Number of particles for summary.
  ullong NpFixed;
  ullong NpMoving;
  ullong NpFloating;
  ullong NpFluid;

//-Attributes for computing domain limits in definition step.
protected:
  const double Hdp;
  const double BuffSizeh;
  const double Overlaph;
  const std::string FileRow;

protected:
  int TrackingMk;       ///< Mk for tracking (-1 by default).
  int TrackingMkBound;  ///< MkBound for tracking (-1 by default).

  //-Variables from simulation domain configuration.
  bool SimUseParent;        ///<Activate the use of parent subdomain.
  bool SimPosXml;           ///<Is TRUE when SDomainPosmin and SDomainPosmax are defined.
  std::string SimPosmin[3]; ///<Values x,y,z of <simulationdomain><posmin>
  std::string SimPosmax[3]; ///<Values x,y,z of <simulationdomain><posmax>

  bool SimPartsPos;         ///<Is TRUE when SimPartsPosmin and SimPartsPosmax are configured.
  tdouble3 SimPartsPosmin;  ///<Minimum position of generated particles.
  tdouble3 SimPartsPosmax;  ///<Maximum position of generated particles.

  JCaseVResExtra Extra_Data;

protected:
  void WriteXmlExtra(JXml* sxml,TiXmlElement* ele)const;

public:
  JCaseVResBase(TpMRShape shape,int id,JCaseVResBase* parent,bool is2d
    ,double posy2d,double hdp,double dp,double bsizeh,double overlaph
    ,std::string filerow);
  virtual ~JCaseVResBase();
  void SetExtraData(const JCaseVResExtra& edata);

  void NpReset();
  void NpSet(ullong npfixed,ullong npmoving,ullong npfloating,ullong npfluid);

  void SimReset();
  void SimReadXmlDef(const JXml* sxml,const TiXmlElement* node);
  void SimConfigPartsPos(tdouble3 partpmin,tdouble3 partpmax);
  double SimComputeValue(int mode,double vdef,double vsize,double vmod)const;
  void SimWriteXmlRun(JXml* sxml,TiXmlElement* ele)const;

  void TrackingDisable();
  void TrackingConfig(int mkbound,int mk);
  bool TrackingIsActive()const{ return(TrackingMk>0); }
  int  GetTrackingMk()const{ return(TrackingMk); }
  int  GetTrackingMkBound()const{ return(TrackingMkBound); }

  std::string GetSubName()const{ return(Shape==MRSH_Box? "Box": "UNKNOWN"); }
  unsigned Count()const{ return(unsigned(SubZones.size())); }
  const JCaseVResBase* GetSubZone(unsigned idx)const;
  const JCaseVRes_Box* GetSubZoneBox(unsigned idx)const;
  std::string GetSubZonesStr()const;

  JCaseVResExtra ExtraData()const{ return(Extra_Data); }
};

//##############################################################################
//# JCaseVRes_Box
//##############################################################################
/// \brief Box buffer zone configuration.
class JCaseVRes_Box : public JCaseVResBase
{
protected:
  //-Configuration values.
  tdouble3 ConfigRef;  ///<Reference position according to the initial configuration.
  JBoxDef  ConfigBox;  ///<Limits according to the initial configuration.

  JBoxDef  PtBox;      ///<Inner limits of subdomain (close to ConfigMin). 
  JBoxDef  BuffBox;    ///<Limits of buffer around the subdomain (PtMin - buffsizeh x h). 
  JBoxDef  FixedBox;   ///<Outer limits of fixed zone around the subdomain (BuffMin - 2.0 x h). 

  JBoxDef  ParentBuffIniBox;   ///<Inner limits of parent buffer within the subdomain (close to PtMin + overlaph x h). 
  JBoxDef  ParentBuffEndBox;   ///<Outer limits of parent buffer within the subdomain (ParentBuffIniMin + parent.buffsizeh x parent.h). 
  JBoxDef  ParentFixedBox;     ///<Limits of parent fixed zone within the subdomain (ParentBuffEndMin + 2.0 x parent.h). 

protected:
  std::string CheckSubOutDomainBox(const std::string& tex
    ,const JBoxDef& subdomain,const JBoxDef& domain);
  tdouble3 FitVectorDp(const tdouble3& v,double dp)const;
  void FitDomain(tdouble3 configmin,tdouble3 configmax
    ,tdouble3& ptmin,tdouble3& ptmax)const;
  JBoxDef FitDomainBox(const JBoxDef& configbox)const;
  
  void ComputeExtendingLimits(const tdouble3& pmin0,const tdouble3& pmax0
    ,double dp,double resizemin,tdouble3& pmin,tdouble3& pmax)const;
  JBoxDef ComputeExtendingLimitsBox(const JBoxDef& box,double dp
    ,double resizemin)const;
  void ComputeReductionLimits(const tdouble3& pmin0,const tdouble3& pmax0
    ,double dp,double resizemin,tdouble3& pmin,tdouble3& pmax)const;

  void ComputeParentBufferLimits(const JBoxDef& ptbox,tdouble3& pinimin,tdouble3& pinimax
    ,tdouble3& pendmin,tdouble3& pendmax,tdouble3& pfixmin,tdouble3& pfixmax)const;

  void ComputeParentBufferLimitsBox(JBoxDef& pini,JBoxDef& pend,JBoxDef& pfix)const;

  void WriteXmlRunBox(JXml* sxml,TiXmlElement* ele,JBoxDef box,double resize)const;
  JBoxDef ReadXmlRunBox(const JXml* sxml,const TiXmlElement* ele,bool is2d)const;

public:
  JCaseVRes_Box(int id,JCaseVResBase* parent,double hdp
    ,double dp ,double bsizeh,double overlaph,std::string filerow
    ,tdouble3 configref,const JBoxDef& configbox);
  JCaseVRes_Box(int id,JCaseVResBase* parent,bool is2d,double posy2d
    ,double dp,const JXml* sxml,TiXmlElement* item);
  void Reset();

  void WriteXmlRun(JXml* sxml,TiXmlElement* ele,bool halfdp)const;
  void ReadXmlRun(const JXml* sxml,TiXmlElement* item);

  tdouble3 GetConfigRef()const{ return(ConfigRef); }

  tdouble3 GetPtRef()const{ return(ConfigRef); }
  JBoxDef  GetPtBox()const{ return(PtBox); }

  JBoxDef  GetBuffBox()const{ return(BuffBox); }

  JBoxDef  GetFixedBox()const{ return(FixedBox); }

  JBoxDef  GetParentBuffIniBox()const{ return(ParentBuffIniBox); }
  JBoxDef  GetParentBuffEndBox()const{ return(ParentBuffEndBox); }
  JBoxDef  GetParentFixedBox()const{ return(ParentFixedBox); }

  bool UseSimpleBoxes()const;
};

//##############################################################################
//# JCaseVRes
//##############################################################################
/// \brief Manages the variable resolution configuration.
class JCaseVRes : protected JObject
{
private:
  bool   Is2D;        ///<2-D domain.
  double Posy2D;      ///<Y-position for 2-D domain.

  //-Variables for speedsound and density gradient. 
  double HSwl;        ///<Maximum height of the volume of fluid computed by GenCase for 1st subdomain.
  double DepthDismin; ///<Depth value for density gradient.

  //-Variables for execution section. 
  bool RunData;       ///<Execution data instead of definition data is loaded.
  unsigned SelId;     ///<Id selected from execution data (UINT_MAX for RunData==false).

  std::vector<JCaseVResBase*> Zones;

private:
  static void GetBoxPoints2d(JBoxDef box,double resize,tdouble3* vpt);
  static void GetBoxPoints3d(JBoxDef box,double resize,tdouble3* ptvec);
  tdouble3 LoadXmlPtref(const JXml* sxml,double dp)const;
  void ReadXmlTransform(const JXml* sxml,const TiXmlElement* ele
    ,bool is2d,tdouble3& mov,tdouble3& rot,tdouble3& rotcen);

  void ReadXmlDef(JCaseVResBase* parent,const JXml* sxml
    ,TiXmlElement* lis,int mkboundfirst,double hdp);
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele
    ,JCaseVResExtra& edata,std::string kparent="");
  bool CheckOutDomain(const tdouble3& ps,const tdouble3& pmin,const tdouble3& pmax){
    return(ps.x<pmin.x || ps.y<pmin.y || ps.z<pmin.z 
      || ps.x>pmax.x || ps.y>pmax.y || ps.z>pmax.z);
  }
  void ReadXmlRun(JCaseVResBase* parent,const JXml* sxml,TiXmlElement* lis);

public:
  JCaseVRes();
  ~JCaseVRes();
  void Reset();

  void SetHSwl(double hswl){ HSwl=hswl; }
  double GetHSwl()const{ return(HSwl); }
  void SetDepthDismin(double v){ DepthDismin=v; }
  double GetDepthDismin()const{ return(DepthDismin); }

  bool   GetIs2D()const{ return(Is2D); }
  double GetPosy2D()const{ return(Posy2D); }

  bool GetRunData()const{ return(RunData); }
  unsigned GetSelId()const{ return(SelId); }
  unsigned Count()const{ return(unsigned(Zones.size())); }
  const JCaseVRes_Box* GetZoneBox(unsigned id)const;

  //-Methods for GenCase.
  void LoadXmlDef(const JXml* sxml,const std::string& place
    ,int mkboundfirst,double hdp,double dp,tdouble3 ptref
    ,tdouble3 ptmin,tdouble3 ptmax);
  void SaveXmlRun(JXml* sxml,const std::string& place,unsigned vresid,bool halfdp)const;
  void SaveXmlSimRun(JXml* sxml,const std::string& place,unsigned vresid)const;
  void NpSet(unsigned id,ullong npfixed,ullong npmoving,ullong npfloating,ullong npfluid);

  //-Methods for DualSPHysics.
  unsigned LoadFileXmlRun(const std::string& file,const std::string& place);
  unsigned LoadXmlRun(const JXml* sxml,const std::string& place);

  void SaveVtkDomains(std::string fname,bool onefile,bool halfdp)const;
  void SaveVtkLimits(std::string fname,bool halfdp)const;
};

#endif

