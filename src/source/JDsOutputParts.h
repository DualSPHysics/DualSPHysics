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
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

class JXml;
class TiXmlElement;
class JVtkLib;

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
  }TpOutputParts; 

  ///Returns damping type as string.
  static std::string GetNameType(TpOutputParts type){
    switch(type){
      case OPA_Init:     return("Init");
      case OPA_Pos:      return("Pos");
      case OPA_Plane:    return("Plane");
      case OPA_Sphere:   return("Sphere");
      case OPA_Cylinder: return("Cylinder");
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

public:
  JDsOutputPartsOp(unsigned id,TpOutputParts type
    ,bool inverse,bool cmband,byte bitres,byte bitprev)
    :Id(id),Type(type),Inverse(inverse),CombineAnd(cmband)
    ,BitResult(bitres),BitResultPrev(bitprev)
  { 
    ClassName=std::string("JDsOutputPartsOp_")+GetNameType(type);
  } 
  virtual ~JDsOutputPartsOp(){ DestructorActive=true; }
  virtual void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele)=0;
  virtual void GetConfig(std::vector<std::string>& lines)const=0;
  virtual void SaveVtkConfig(double size,JVtkLib* sh)const{ };

  virtual void ComputeFilterCpu(unsigned np,const tdouble3* pos
    ,const typecode* code,byte* sel)const=0;
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
    ,const JXml *sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Init,false,false,bitres,0)
    ,SelAll(selall)
  { ReadXmlExtra(sxml,ele); } 
  void Reset(){}
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele){ Reset(); }
  void GetConfig(std::vector<std::string>& lines)const;

  void ComputeFilterCpu(unsigned np,const tdouble3* pos
    ,const typecode* code,byte* sel)const;
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
  tdouble3 PosMin;
  tdouble3 PosMax;

public:
  JDsOutputPartsOp_Pos(unsigned id,bool invert,bool cmband
    ,byte bitres,const JXml *sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Pos,invert,cmband,bitres,0)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JVtkLib* sh)const;

  void ComputeFilterCpu(unsigned np,const tdouble3* pos
    ,const typecode* code,byte* sel)const;
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
  tdouble3 Point;
  tdouble3 Vector;
  double Distance;
  tplane3d Plane;

public:
  JDsOutputPartsOp_Plane(unsigned id,bool invert,bool cmband
    ,byte bitres,const JXml *sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Plane,invert,cmband,bitres,0)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JVtkLib* sh)const;

  void ComputeFilterCpu(unsigned np,const tdouble3* pos
    ,const typecode* code,byte* sel)const;
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
  tdouble3 Centre;
  double Radius;

public:
  JDsOutputPartsOp_Sphere(unsigned id,bool invert,bool cmband
    ,byte bitres,const JXml *sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Sphere,invert,cmband,bitres,0)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JVtkLib* sh)const;

  void ComputeFilterCpu(unsigned np,const tdouble3* pos
    ,const typecode* code,byte* sel)const;
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
  tdouble3 Point1;
  tdouble3 Point2;
  double Radius;
  tplane3d Plane;
  double Distance;

public:
  JDsOutputPartsOp_Cylinder(unsigned id,bool invert,bool cmband
    ,byte bitres,const JXml *sxml,const TiXmlElement* ele)
    :JDsOutputPartsOp(id,OPA_Cylinder,invert,cmband,bitres,0)
  { ReadXmlExtra(sxml,ele); } 
  void Reset();
  void ReadXmlExtra(const JXml* sxml,const TiXmlElement* ele);
  void GetConfig(std::vector<std::string>& lines)const;
  void SaveVtkConfig(double size,JVtkLib* sh)const;

  void ComputeFilterCpu(unsigned np,const tdouble3* pos
    ,const typecode* code,byte* sel)const;
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

  std::vector<JDsOutputPartsOp*> List;

  void ReadXml(const JXml* sxml,const TiXmlElement* ele);
  void SaveVtkConfig(double size)const;

public:
  JDsOutputParts(bool cpu,double vtksize);
  ~JDsOutputParts();
  void Reset();

  unsigned Count()const{ return(unsigned(List.size())); }

  void LoadXml(const JXml* sxml,const std::string& place);
  void GetConfig(std::string txhead,std::string txfoot
    ,std::vector<std::string>& lines)const;

  void ComputeFilterCpu(unsigned n,const tdouble3* pos
    ,const typecode* code,byte* sel)const;
#ifdef _WITHGPU
  void ComputeFilterGpu(unsigned np,const double2* posxy
    ,const double*posz,const typecode* code,byte*sel)const;
#endif

};

#endif


