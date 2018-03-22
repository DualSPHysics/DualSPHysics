//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Se guarda el MkBoundFirst y el MkFluidFirst para poder establecer una
//:#   relacion con el mk de los poligonos creados. (27-11-2010)
//:# - Nuevos metodos LoadFileXml() y SaveFileXml() para cargar o generar un
//:#   fichero xml de forma directa. (27-11-2010)
//:# - Se guarda el Mk absoluto para cada bloque de particulas para facilitar
//:#   su utilizacion en herramientas de postprocesing. (21-01-2011)
//:# - Se añadieron dos nuevas variables para los floating bodies: Velini y 
//:#   Omegaini. (25-01-2011)
//:# - Cambio de nombre, de JParticles a JSpaceParts. (09-02-2012)
//:# - Traduccion de comentarios al ingles. (10-02-2012)
//:# - El MK pasa a ser de tipo word. (23-11-2013)
//:# - Las variables float pasan a ser double. (23-11-2013)
//:# - Clase JSpacePartsGetMk para calcular Mk a partir de Id. (23-11-2013)
//:# - Gestiona Properties para bloques de particulas. (14-12-2013)
//:# - Introduccion de valores por defecto en metodos GetSubValue(). (24-07-2014)
//:# - Permite cambiar el numero de particulas de cada bloque. (13-08-2014)
//:# - JSpacePartsGetMk devuelve el ultimo mk de fluido para las particulas 
//:#   creadas con splitting. (12-05-2015)
//:# - Se graba el valor masspart para floating bodies en el fichero XML. (05-01-2016)
//:# - Se graban los datos de floatingns con las unidades. (29-01-2016)
//:# - Genera resumen de particulas y bloques de MK. (03-08-2017)
//:# - Inertia de floatings se guarda como tmatrix3d. (29-11-2017)
//:#############################################################################

/// \file JSpaceParts.h \brief Declares the class \ref JSpaceParts.

#ifndef _JSpaceParts_
#define _JSpaceParts_

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "JObject.h"
#include "TypesDef.h"

class JXml;
class TiXmlElement;
class JSpaceProperties;

///<Types of particles.
typedef enum{ PT_Fixed=1,PT_Moving=2,PT_Floating=4,PT_Fluid=5 }TpParticles; 

//##############################################################################
//# JSpacePartBlock
//##############################################################################
/// \brief Manages the info of each block of particles from the input XML file.
class JSpacePartBlock : public JObject
{
private:
  const JSpaceProperties* Properties;   ///<Pointer to properties object.
  std::string Props;                    ///<Assigned properties.
  word Mk;                              ///<Absolute label.
  word MkType;                          ///<Label of block fluid or bound.
  unsigned Begin;                       ///<Id of the first particle of the block.
  unsigned Count;                       ///<Number of particles.

public:
  const TpParticles Type;    ///<Type of particle.
  const bool Bound;          ///<Indicates whether a particle is boundary or not.

  JSpacePartBlock(const JSpaceProperties* properties,TpParticles type,const char* name,bool bound,word mktype=0,unsigned begin=0,unsigned count=0):Properties(properties),Type(type),Bound(bound),MkType(mktype),Begin(begin),Count(count){ 
    ClassName=std::string("JSpacePartBlock_")+name;
  } 
  virtual ~JSpacePartBlock(){ DestructorActive=true; }
  void UpdateProperty();
  void ConfigMk(word mkfirst){ Mk=MkType+mkfirst; }
  std::string GetNameXml()const;
  unsigned GetBegin()const{ return(Begin); }
  unsigned GetCount()const{ return(Count); }
  word GetMkType()const{ return(MkType); }
  word GetMk()const{ return(Mk); }
  std::string GetProperty()const{ return(Props); }
  virtual void ReadXml(JXml *sxml,TiXmlElement* ele);
  virtual TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;

  void SetBegin(unsigned begin){ Begin=begin; }
  void SetCount(unsigned count){ Count=count; }

  //-Returns values of properties.
  unsigned GetValuesCount()const;
  std::string GetValueName(unsigned idx)const;
  std::string GetValueStr(unsigned idx)const;
  int GetValueInt(unsigned idx)const{        return(atoi(GetValueStr(idx).c_str())); }
  unsigned GetValueUint(unsigned idx)const{  return(unsigned(GetValueInt(idx)));     }
  double GetValueDouble(unsigned idx)const{  return(atof(GetValueStr(idx).c_str())); }
  float GetValueFloat(unsigned idx)const{    return(float(GetValueDouble(idx)));     }
  bool ExistsValue(std::string name)const;
  std::string GetValueStr(std::string name)const;
  int GetValueInt(std::string name)const{        return(atoi(GetValueStr(name).c_str())); }
  unsigned GetValueUint(std::string name)const{  return(unsigned(GetValueInt(name)));     }
  double GetValueDouble(std::string name)const{  return(atof(GetValueStr(name).c_str())); }
  float GetValueFloat(std::string name)const{    return(float(GetValueDouble(name)));     }
  
  //-Returns subvalues of properties.
  unsigned GetSubValuesCount(unsigned idx)const;
  std::string GetSubValueName(unsigned idx,unsigned subidx)const;
  std::string GetSubValueStr(unsigned idx,unsigned subidx)const;
  int GetSubValueInt(unsigned idx,unsigned subidx)const{        return(atoi(GetSubValueStr(idx,subidx).c_str())); }
  unsigned GetSubValueUint(unsigned idx,unsigned subidx)const{  return(unsigned(GetSubValueInt(idx,subidx)));     }
  double GetSubValueDouble(unsigned idx,unsigned subidx)const{  return(atof(GetSubValueStr(idx,subidx).c_str())); }
  float GetSubValueFloat(unsigned idx,unsigned subidx)const{    return(float(GetSubValueDouble(idx,subidx)));     }
  bool ExistsSubValue(std::string name,std::string subname)const;

  std::string GetSubValueStr(std::string name,std::string subname,bool optional=false,std::string valdef="")const;
  int GetSubValueInt(std::string name,std::string subname,bool optional=false,int valdef=0)const;
  unsigned GetSubValueUint(std::string name,std::string subname,bool optional=false,unsigned valdef=0)const{  return(unsigned(GetSubValueInt(name,subname,optional,int(valdef))));  }
  double GetSubValueDouble(std::string name,std::string subname,bool optional=false,double valdef=0)const;
  float GetSubValueFloat(std::string name,std::string subname,bool optional=false,float valdef=0)const{  return(float(GetSubValueDouble(name,subname,optional,valdef)));  }
  tdouble3 GetSubValueDouble3(std::string name)const{ return(TDouble3(GetSubValueDouble(name,"x"),GetSubValueDouble(name,"y"),GetSubValueDouble(name,"z"))); }
  tfloat3 GetSubValueFloat3(std::string name)const{ return(ToTFloat3(GetSubValueDouble3(name))); }
};

//##############################################################################
//# JSpacePartBlock_Fixed
//##############################################################################
/// Manages the info of boundary fixed particles.
class JSpacePartBlock_Fixed : public JSpacePartBlock
{
public:
  JSpacePartBlock_Fixed(const JSpaceProperties* properties,word mktype,unsigned begin,unsigned count):JSpacePartBlock(properties,PT_Fixed,"Fixed",true,mktype,begin,count){}
  JSpacePartBlock_Fixed(const JSpaceProperties* properties,JXml *sxml,TiXmlElement* ele):JSpacePartBlock(properties,PT_Fixed,"Fixed",true){ ReadXml(sxml,ele); }
};  

//##############################################################################
//# JSpacePartBlock_Moving
//##############################################################################
/// Manages the info of boundary moving particles.
class JSpacePartBlock_Moving : public JSpacePartBlock
{
private:
  unsigned RefMotion;
public:
  JSpacePartBlock_Moving(const JSpaceProperties* properties,word mktype,unsigned begin,unsigned count,unsigned refmotion):JSpacePartBlock(properties,PT_Moving,"Moving",true,mktype,begin,count),RefMotion(refmotion){}
  JSpacePartBlock_Moving(const JSpaceProperties* properties,JXml *sxml,TiXmlElement* ele):JSpacePartBlock(properties,PT_Moving,"Moving",true){ ReadXml(sxml,ele); }
  unsigned GetRefMotion()const{ return(RefMotion); }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;
};  

//##############################################################################
//# JSpacePartBlock_Floating
//##############################################################################
/// Manages the info of floating particles.
class JSpacePartBlock_Floating : public JSpacePartBlock
{
private:
  double Massbody;
  tdouble3 Center;
  tmatrix3d Inertia;
  tdouble3 Velini;
  tdouble3 Omegaini;
public:
  JSpacePartBlock_Floating(const JSpaceProperties* properties,word mktype,unsigned begin,unsigned count,double massbody,const tdouble3& center,const tmatrix3d& inertia,const tdouble3& velini,const tdouble3& omegaini):JSpacePartBlock(properties,PT_Floating,"Floating",true,mktype,begin,count),Massbody(massbody),Center(center),Inertia(inertia),Velini(velini),Omegaini(omegaini){}
  JSpacePartBlock_Floating(const JSpaceProperties* properties,JXml *sxml,TiXmlElement* ele):JSpacePartBlock(properties,PT_Floating,"Floating",true){ ReadXml(sxml,ele); }
  double GetMassbody()const{ return(Massbody); }
  tdouble3 GetCenter()const{ return(Center); }
  tmatrix3d GetInertia()const{ return(Inertia); }
  tdouble3 GetVelini()const{ return(Velini); }
  tdouble3 GetOmegaini()const{ return(Omegaini); }
  void ReadXml(JXml *sxml,TiXmlElement* ele);
  TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;
};  

//##############################################################################
//# JSpacePartBlock_Fluid
//##############################################################################
/// Manages the info of fluid particles.
class JSpacePartBlock_Fluid : public JSpacePartBlock
{
public:
  JSpacePartBlock_Fluid(const JSpaceProperties* properties,word mktype,unsigned begin,unsigned count):JSpacePartBlock(properties,PT_Fluid,"Fluid",false,mktype,begin,count){}
  JSpacePartBlock_Fluid(const JSpaceProperties* properties,JXml *sxml,TiXmlElement* ele):JSpacePartBlock(properties,PT_Fluid,"Fluid",false){ ReadXml(sxml,ele); }
};  

//##############################################################################
//# JSpaceParts
//##############################################################################
/// \brief Manages the info of particles from the input XML file.

class JSpaceParts  : protected JObject
{
public:
  /// Structure with summary of particle information.
  typedef struct {
    unsigned np[4];
    unsigned idini[4];
    unsigned idlast[4];
    unsigned nmk[4];
    std::string mklist[4];
  }StSummaryData;

private:
  std::vector<JSpacePartBlock*> Blocks;
  unsigned Begin;
  TpParticles LastType;
  word MkBoundFirst,MkFluidFirst;
  JSpaceProperties* Properties;
  
  unsigned GetBegin()const{ return(Begin); }
  JSpacePartBlock* GetByMkType(bool bound,word mktype)const;
  void Add(JSpacePartBlock* block);
  void ReadXml(JXml *sxml,TiXmlElement* lis);
  void WriteXml(JXml *sxml,TiXmlElement* lis)const;
  void WriteXmlSummary(JXml *sxml,TiXmlElement* ele)const;

public:
  JSpaceParts();
  ~JSpaceParts();
  void Reset();
  unsigned Count(TpParticles type)const;
  unsigned Count()const;
  unsigned CountBlocks()const{ return(unsigned(Blocks.size())); }
  unsigned CountBlocks(TpParticles type)const;
  const JSpacePartBlock& GetBlock(unsigned pos)const;

  void LoadFileXml(const std::string &file,const std::string &path);
  void SaveFileXml(const std::string &file,const std::string &path,bool newfile=true)const;
  void LoadXml(JXml *sxml,const std::string &place);
  void SaveXml(JXml *sxml,const std::string &place)const;

  void SetMkFirst(word boundfirst,word fluidfirst);
  word GetMkBoundFirst()const{ return(MkBoundFirst); }
  word GetMkFluidFirst()const{ return(MkFluidFirst); }

  void AddFixed(word mktype,unsigned count){ Add(new JSpacePartBlock_Fixed(Properties,mktype,GetBegin(),count)); }
  void AddMoving(word mktype,unsigned count,unsigned refmotion){ Add(new JSpacePartBlock_Moving(Properties,mktype,GetBegin(),count,refmotion)); }
  void AddFloating(word mktype,unsigned count,double massbody,const tdouble3& center,const tmatrix3d& inertia,const tdouble3& velini,const tdouble3& omegaini){ Add(new JSpacePartBlock_Floating(Properties,mktype,GetBegin(),count,massbody,center,inertia,velini,omegaini)); }
  void AddFluid(word mktype,unsigned count){ Add(new JSpacePartBlock_Fluid(Properties,mktype,GetBegin(),count)); }

  void SetBlockSize(unsigned pos,unsigned np);

  void LoadProperties(const JSpaceProperties *props);

  std::string GetMkList(TpParticles type)const;

  JSpaceParts::StSummaryData GetSummaryData()const;
  void GetParticleSummary(std::vector<std::string> &out)const;
};


//##############################################################################
//# JSpacePartsGetMk
//##############################################################################
/// \brief Compute Mk value from Id.

class JSpacePartsGetMk  : protected JObject
{
private:
  const bool Splitting;
  unsigned MkCount;
  unsigned *MkRange;
  word *MkValue;
  word MkSplitting;

  void Config(const JSpaceParts *sparts);
public:
  JSpacePartsGetMk(const JSpaceParts *sparts,bool splitting);
  ~JSpacePartsGetMk();
  void Reset();
  unsigned GetMkCount()const{ return(MkCount); }

  ///Returns MK from Id value.
  inline word IdToMk(unsigned id)const{
    unsigned c=0;
    for(;c<MkCount&&id>=MkRange[c];c++);
    return(c<MkCount? MkValue[c]: MkSplitting);
  }
};

#endif


