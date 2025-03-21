//HEAD_DSPH
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
//:# - Gestiona el filtro de particulas a guardar en los ficheros PART_XXXX. (03-04-2023)
//:#############################################################################

/// \file JDsOutputParts.h \brief Declares the class \ref JDsOutputParts.

#ifndef _JDsOutputParts_
#define _JDsOutputParts_

#include "JObject.h"
#include "DualSphDef.h"
#include <string>
#include <vector>
#include <climits>
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

class JXml;
class TiXmlElement;
class JSpVtkShape;
class JSphMk;

//##############################################################################
//# XML format in _FmtXML_OutputParts.xml.
//##############################################################################

//##############################################################################
//# JDsOutputPartsOp
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp : public JObject
{
public:
  ///Types of damping configurations.
  typedef enum{ 
    OPA_Init=1
   ,OPA_Pos=2
   ,OPA_Plane=3
   ,OPA_Sphere=4
   ,OPA_Cylinder=5
   ,OPA_Type=6
   ,OPA_Mk=7
   ,OPA_GroupFin=8
  }TpOutputParts; 

  ///Returns damping type as string.
  static std::string GetNameType(TpOutputParts type){
    switch(type){
      case OPA_Init:     return("Init");
      case OPA_Pos:      return("Pos");
      case OPA_Plane:    return("Plane");
      case OPA_Sphere:   return("Sphere");
      case OPA_Cylinder: return("Cylinder");
      case OPA_Type:     return("Type");
      case OPA_Mk:       return("Mk"); 
      case OPA_GroupFin: return("GroupFin");
    }
    return("???");
  }

public:
  const unsigned Id;
  const TpOutputParts Type;
  const bool Inverse;       ///<Inverse selection.
  const bool CombineAnd;    ///<Combine with result using AND (instead of OR).
  const byte BitResult;     ///<Bit to save result (0-7).
  const byte BitResultPrev; ///<Bit of previous result to combine (0-7).

  const unsigned FtFollowId;   ///<Id of floating to follow in FtObjs[] (disabled when UINT_MAX).
  const unsigned FtFollowMkb;  ///<Mkbound of floating to follow (disabled when UINT_MAX).
  const tdouble3 FtFollowCen0; ///<Initial center of floating to follow.

public:
  JDsOutputPartsOp(unsigned id,TpOutputParts type
    ,bool inverse,bool cmband,byte bitres,byte bitprev
    ,unsigned ftid,unsigned ftmkb=UINT_MAX,tdouble3 ftcen=TDouble3(0))
    :Id(id),Type(type),Inverse(inverse),CombineAnd(cmband)
    ,BitResult(bitres),BitResultPrev(bitprev)
    ,FtFollowId(ftid),FtFollowMkb(ftmkb),FtFollowCen0(ftcen)
  { 
    ClassName=std::string("JDsOutputPartsOp_")+GetNameType(type);
  } 
  virtual ~JDsOutputPartsOp(){ DestructorActive=true; }
  virtual void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele)=0;
  virtual void UpdateFtPos(const tdouble3& ftcenter)=0;
  virtual void GetConfig(std::vector<std::string>& lines)const=0;
  virtual void SaveVtkConfig(double size,JSpVtkShape* ss)const{ };

  virtual void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const=0;
#ifdef _WITHGPU
  virtual void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const=0;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_Init
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_Init : public JDsOutputPartsOp
{
protected:
  const bool SelAll;

public:
  JDsOutputPartsOp_Init(unsigned id,bool selall,byte bitres
    ,const JXml* sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Init,false,false,bitres,0,UINT_MAX)
    ,SelAll(selall)
  { ReadXmlExtra(sxml,ele); } 
  void Reset(){}
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){ Reset(); }
  void UpdateFtPos(const tdouble3& ftcenter){ };
  void GetConfig(std::vector<std::string>& lines)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_GroupFin
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_GroupFin : public JDsOutputPartsOp
{
protected:
  const bool SelAll;

public:
  JDsOutputPartsOp_GroupFin(unsigned id,bool selall
    ,bool inverse,bool cmband,byte bitres,byte bitprev
    ,const JXml* sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_GroupFin,inverse,cmband,bitres,bitprev,UINT_MAX)
    ,SelAll(selall)
  { ReadXmlExtra(sxml,ele); } 
  void Reset(){}
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){ Reset(); }
  void UpdateFtPos(const tdouble3& ftcenter){ };
  void GetConfig(std::vector<std::string>& lines)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_Pos
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_Pos : public JDsOutputPartsOp
{
protected:
  tdouble3 PosMinInit;
  tdouble3 PosMaxInit;
  tdouble3 PosMin;
  tdouble3 PosMax;

public:
  JDsOutputPartsOp_Pos(unsigned id,bool inverse,bool cmband
    ,byte bitres,unsigned ftid,unsigned ftmkb,tdouble3 ftcen
    ,const JXml* sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Pos,inverse,cmband,bitres,0
      ,ftid,ftmkb,ftcen)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void UpdateFtPos(const tdouble3& ftcenter);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JSpVtkShape* ss)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_Plane
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_Plane: public JDsOutputPartsOp
{
protected:
  tdouble3 PointInit;
  tdouble3 Point;
  tdouble3 Vector;
  double Distance;
  tplane3d Plane;

public:
  JDsOutputPartsOp_Plane(unsigned id,bool inverse,bool cmband
    ,byte bitres,unsigned ftid,unsigned ftmkb,tdouble3 ftcen
    ,const JXml* sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Plane,inverse,cmband,bitres,0
      ,ftid,ftmkb,ftcen)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void UpdateFtPos(const tdouble3& ftcenter);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JSpVtkShape* ss)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_Sphere
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_Sphere: public JDsOutputPartsOp
{
protected:
  tdouble3 CentreInit; ///<Initial centre.
  tdouble3 Centre;     ///<Updated centre.
  double Radius;

public:
  JDsOutputPartsOp_Sphere(unsigned id,bool inverse,bool cmband
    ,byte bitres,unsigned ftid,unsigned ftmkb,tdouble3 ftcen
    ,const JXml* sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Sphere,inverse,cmband,bitres,0
      ,ftid,ftmkb,ftcen)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void UpdateFtPos(const tdouble3& ftcenter);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JSpVtkShape* ss)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_Cylinder
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_Cylinder: public JDsOutputPartsOp
{
protected:
  tdouble3 Point1Init;
  tdouble3 Point2Init;
  tdouble3 Point1;
  tdouble3 Point2;
  double Radius;
  tplane3d Plane;
  double Distance;

public:
  JDsOutputPartsOp_Cylinder(unsigned id,bool inverse,bool cmband
    ,byte bitres,unsigned ftid,unsigned ftmkb,tdouble3 ftcen
    ,const JXml* sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Cylinder,inverse,cmband,bitres,0
      ,ftid,ftmkb,ftcen)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void UpdateFtPos(const tdouble3& ftcenter);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JSpVtkShape* ss)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_Type
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_Type : public JDsOutputPartsOp
{
protected:
  byte Types;

public:
  JDsOutputPartsOp_Type(unsigned id,bool inverse,bool cmband
    ,byte bitres,const JXml* sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Type,inverse,cmband,bitres,0,UINT_MAX)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void UpdateFtPos(const tdouble3& ftcenter){ };
  void GetConfig(std::vector<std::string>& lines)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputPartsOp_Mk
//##############################################################################
/// \brief Base clase for filter configurations.
class JDsOutputPartsOp_Mk : public JDsOutputPartsOp
{
protected:
  std::string MkDef;
  unsigned Mk1,Mk2;
  unsigned MkCode1,MkCode2;
  const JSphMk* MkInfo;

public:
  JDsOutputPartsOp_Mk(unsigned id,bool inverse,bool cmband
    ,byte bitres,const JXml* sxml,const TiXmlElement* ele
    ,const JSphMk* mkinfo)
    :JDsOutputPartsOp(id,OPA_Mk,inverse,cmband,bitres,0,UINT_MAX),MkInfo(mkinfo)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void UpdateFtPos(const tdouble3& ftcenter){ };
  void GetConfig(std::vector<std::string>& lines)const;

  void ComputeFilterCpu(unsigned np,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif
};


//##############################################################################
//# JDsOutputParts
//##############################################################################
/// \brief Manage the filter to save particles in files PART_XXXX.
class JDsOutputParts : protected JObject
{
private:
  const bool Cpu;
  const double VtkSize;
  int Nparts;  ///<Ignore filters every nparts output files (default=0).

  std::vector<JDsOutputPartsOp*> List;

  void ReadXml(const JXml* sxml,const TiXmlElement* ele,byte level
    ,const JSphMk* mkinfo,unsigned ftcount,const StFloatingData* ftobjs);
  void SaveVtkConfig(double size)const;

public:
  JDsOutputParts(bool cpu,double vtksize);

  ~JDsOutputParts();
  void Reset();

  unsigned Count()const{ return(unsigned(List.size())); }

  void LoadXml(const JXml* sxml,const std::string& place
    ,const JSphMk* mkinfo,unsigned ftcount,const StFloatingData* ftobjs);
  void GetConfig(std::string txhead,std::string txfoot
    ,std::vector<std::string>& lines)const;

  bool CheckFilters(int cpart)const;
  void UpdateFtPos(unsigned ftcount,const StFloatingData* ftobjs);

  void ComputeFilterCpu(unsigned n,const unsigned* idp
    ,const tdouble3* pos,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif

};

#endif


