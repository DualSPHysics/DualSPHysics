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
//:# - Se guarda el MkBoundFirst y el MkFluidFirst para poder establecer una
//:#   relacion con el mk de los poligonos creados. (27-11-2010)
//:# - Nuevos metodos LoadFileXml() y SaveFileXml() para cargar o generar un
//:#   fichero xml de forma directa. (27-11-2010)
//:# - Se guarda el Mk absoluto para cada bloque de particulas para facilitar
//:#   su utilizacion en herramientas de postprocesing. (21-01-2011)
//:# - Se anhadieron dos nuevas variables para los floating bodies: Velini y 
//:#   Omegaini. (25-01-2011)
//:# - Cambio de nombre, de JParticles a J.SpaceParts. (09-02-2012)
//:# - Traduccion de comentarios al ingles. (10-02-2012)
//:# - El MK pasa a ser de tipo word. (23-11-2013)
//:# - Las variables float pasan a ser double. (23-11-2013)
//:# - Clase JCasePartsGetMk para calcular Mk a partir de Id. (23-11-2013)
//:# - Gestiona Properties para bloques de particulas. (14-12-2013)
//:# - Introduccion de valores por defecto en metodos GetSubValue(). (24-07-2014)
//:# - Permite cambiar el numero de particulas de cada bloque. (13-08-2014)
//:# - JCasePartsGetMk devuelve el ultimo mk de fluido para las particulas 
//:#   creadas con splitting. (12-05-2015)
//:# - Se graba el valor masspart para floating bodies en el fichero XML. (05-01-2016)
//:# - Se graban los datos de floatingns con las unidades. (29-01-2016)
//:# - Genera resumen de particulas y bloques de MK. (03-08-2017)
//:# - Inertia de floatings se guarda como tmatrix3d. (29-11-2017)
//:# - Define SPACEPARTS_USE_INERTIA_OLD para grabar el vector inertia como en  
//:#   la version v4.0.055. (07-03-2018)
//:# - Se pasa a usar los tipos de particulas definidos en JParticlesDef.h. (23-03-2018)
//:# - Metodos para obtener informacion de particulas. (25-04-2018)
//:# - Se cambian variables Velini y Omegaini por LinearVelini y AngularVelini. (13-02-2019)
//:# - Permite definir restricciones al calcular aceleracion lineal (TranslationFree)
//:#   y angular (RotationFree). (13-02-2019)
//:# - Permite definir una lista de velocidades peredefinidas (lineal y angular) en el XML
//:#   o en un fichero independiente. (21-02-2019)
//:# - Gestiona posicion minima y maxima de particulas. (21-03-2019)
//:# - Uso de JCasePartsDef.h para omitir codigo de JCaseProperties, 
//:#   JLinearValue y JReadDatafile. (16-05-2019)
//:# - Cambia 'translation' y 'rotation' por 'translationDOF' y 'rotationDOF' 
//:#   manteniendo compatibilidad. (14-12-2019)
//:# - Los floating incluyen el valor masspart. (10-03-2020)
//:# - Objeto JXml pasado como const para operaciones de lectura. (17-03-2020)  
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:# - Cambio de nombre de J.SpaceParts a J.CaseParts. (28-06-2020)
//:# - Nuevas opciones LinearForce y AngularForce para floatings. (30-06-2020)
//:# - Se definio BaseNameSize para obtener nombre de subclase. (02-07-2020)
//:#############################################################################

/// \file JCaseParts.h \brief Declares the class \ref JCaseParts.

#ifndef _JCaseParts_
#define _JCaseParts_

#include "JObject.h"
#include "TypesDef.h"
#include "JCasePartsDef.h"
#include "JParticlesDef.h"

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

//#define SPACEPARTS_USE_INERTIA_OLD  //-Saves the inertia vector as in GenCase v4.0.055.

class JXml;
class TiXmlElement;
#ifdef JCaseParts_UseProps
  class JCaseProperties;
#else
  #define JCaseProperties void
#endif

#ifdef JCaseParts_UseLinearv
  class JLinearValue;
#else
  #define JLinearValue void
#endif

//##############################################################################
//# JCasePartBlock
//##############################################################################
/// \brief Manages the info of each block of particles from the input XML file.
class JCasePartBlock : public JObject
{
private:
  const JCaseProperties* Properties;   ///<Pointer to properties object.
  std::string Props;                    ///<Assigned properties.
  word Mk;                              ///<Absolute label.
  word MkType;                          ///<Label of block fluid or bound.
  unsigned Begin;                       ///<Id of the first particle of the block.
  unsigned Count;                       ///<Number of particles.

public:
  const unsigned BaseNameSize;
  const TpParticles Type;    ///<Type of particle.
  const bool Bound;          ///<Indicates whether a particle is boundary or not.

  JCasePartBlock(const JCaseProperties* properties,TpParticles type,const char* name
    ,word mktype=0,unsigned begin=0,unsigned count=0)
    :BaseNameSize(unsigned(std::string("JCasePartBlock_").size()))
    ,Properties(properties),Type(type),Bound(IsBound(type)),MkType(mktype),Begin(begin),Count(count)
  { 
    ClassName=std::string("JCasePartBlock_")+name;
  } 
  virtual ~JCasePartBlock(){ DestructorActive=true; }
#ifdef JCaseParts_UseProps
  void UpdateProperty();
#endif
  void ConfigMk(word mkfirst){ Mk=MkType+mkfirst; }
  std::string GetNameXml()const;
  unsigned GetBegin()const{ return(Begin); }
  unsigned GetCount()const{ return(Count); }
  word GetMkType()const{ return(MkType); }
  word GetMk()const{ return(Mk); }
  std::string GetProperty()const{ return(Props); }
  virtual void ReadXml(const JXml *sxml,TiXmlElement* ele);
  virtual TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;

  void SetBegin(unsigned begin){ Begin=begin; }
  void SetCount(unsigned count){ Count=count; }

#ifdef JCaseParts_UseProps
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
#endif
};

//##############################################################################
//# JCasePartBlock_Fixed
//##############################################################################
/// Manages the info of boundary fixed particles.
class JCasePartBlock_Fixed : public JCasePartBlock
{
public:
  JCasePartBlock_Fixed(const JCaseProperties* properties,word mktype,unsigned begin,unsigned count)
    :JCasePartBlock(properties,TpPartFixed,"Fixed",mktype,begin,count){}
  JCasePartBlock_Fixed(const JCaseProperties* properties,const JXml *sxml,TiXmlElement* ele)
    :JCasePartBlock(properties,TpPartFixed,"Fixed"){ ReadXml(sxml,ele); }
};  

//##############################################################################
//# JCasePartBlock_Moving
//##############################################################################
/// Manages the info of boundary moving particles.
class JCasePartBlock_Moving : public JCasePartBlock
{
private:
  unsigned RefMotion;
public:
  JCasePartBlock_Moving(const JCaseProperties* properties,word mktype,unsigned begin
    ,unsigned count,unsigned refmotion)
    :JCasePartBlock(properties,TpPartMoving,"Moving",mktype,begin,count),RefMotion(refmotion){}
  JCasePartBlock_Moving(const JCaseProperties* properties,const JXml *sxml,TiXmlElement* ele)
    :JCasePartBlock(properties,TpPartMoving,"Moving"){ ReadXml(sxml,ele); }
  unsigned GetRefMotion()const{ return(RefMotion); }
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;
};  

//##############################################################################
//# JCasePartBlock_Floating
//##############################################################################
/// Manages the info of floating particles.
class JCasePartBlock_Floating : public JCasePartBlock
{
private:
  double Massbody;
  double Masspart;
  tdouble3 Center;
  tmatrix3d Inertia;
  tint3 TranslationFree;
  tint3 RotationFree;
  tdouble3 LinearVelini;
  tdouble3 AngularVelini;
  JLinearValue *LinearVel;
  JLinearValue *AngularVel;
  JLinearValue *LinearForce;
  JLinearValue *AngularForce;

public:
  JCasePartBlock_Floating(const JCaseProperties* properties
    ,word mktype,unsigned begin,unsigned count,double massbody,double masspart
    ,const tdouble3& center,const tmatrix3d& inertia
    ,const tint3 &translationfree,const tint3 &rotationfree
    ,const tdouble3 &linvelini,const tdouble3 &angvelini
    ,const JLinearValue *linvel  ,const JLinearValue *angvel
    ,const JLinearValue *linforce,const JLinearValue *angforce);
  JCasePartBlock_Floating(const JCaseProperties* properties,const JXml *sxml,TiXmlElement* ele);
  ~JCasePartBlock_Floating();
  double        GetMassbody()       const{ return(Massbody); }
  double        GetMasspart()       const{ return(Masspart); }
  tdouble3      GetCenter()         const{ return(Center); }
  tmatrix3d     GetInertia()        const{ return(Inertia); }
  tint3         GetTranslationFree()const{ return(TranslationFree); }
  tint3         GetRotationFree()   const{ return(RotationFree); }
  tdouble3      GetLinearVelini()   const{ return(LinearVelini); }
  tdouble3      GetAngularVelini()  const{ return(AngularVelini); }
  JLinearValue* GetLinearVel()      const{ return(LinearVel); }
  JLinearValue* GetAngularVel()     const{ return(AngularVel); }
  JLinearValue* GetLinearForce()    const{ return(LinearForce); }
  JLinearValue* GetAngularForce()   const{ return(AngularForce); }
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  TiXmlElement* WriteXml(JXml *sxml,TiXmlElement* ele)const;
};  

//##############################################################################
//# JCasePartBlock_Fluid
//##############################################################################
/// Manages the info of fluid particles.
class JCasePartBlock_Fluid : public JCasePartBlock
{
public:
  JCasePartBlock_Fluid(const JCaseProperties* properties,word mktype,unsigned begin,unsigned count)
    :JCasePartBlock(properties,TpPartFluid,"Fluid",mktype,begin,count){}
  JCasePartBlock_Fluid(const JCaseProperties* properties,const JXml *sxml,TiXmlElement* ele)
    :JCasePartBlock(properties,TpPartFluid,"Fluid"){ ReadXml(sxml,ele); }
};  

//##############################################################################
//# JCaseParts
//##############################################################################
/// \brief Manages the info of particles from the input XML file.

class JCaseParts  : protected JObject
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
  std::vector<JCasePartBlock*> Blocks;
  unsigned Begin;
  TpParticles LastType;
  word MkBoundFirst,MkFluidFirst;
  tdouble3 Posmin;  ///<Minimum position of particles.
  tdouble3 Posmax;  ///<Maximum position of particles.

  JCaseProperties* Properties;

  unsigned GetBegin()const{ return(Begin); }
  JCasePartBlock* GetByMkType(bool bound,word mktype)const;
  void Add(JCasePartBlock* block);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void WriteXml(JXml *sxml,TiXmlElement* lis)const;
  void WriteXmlSummary(JXml *sxml,TiXmlElement* ele)const;

public:
  JCaseParts();
  ~JCaseParts();
  void Reset();
  unsigned Count(TpParticles type)const;
  unsigned Count()const;
  unsigned CountBlocks()const{ return(unsigned(Blocks.size())); }
  unsigned CountBlocks(TpParticles type)const;
  const JCasePartBlock& GetBlock(unsigned pos)const;

  bool CheckNparticles(unsigned casenfixed,unsigned casenmoving,unsigned casenfloat,unsigned casenfluid)const;
  bool CheckNparticles(ullong casenfixed,ullong casenmoving,ullong casenfloat,ullong casenfluid)const;

  void LoadFileXml(const std::string &file,const std::string &path);
  void SaveFileXml(const std::string &file,const std::string &path,bool newfile=true)const;
  void LoadXml(const JXml *sxml,const std::string &place);
  void SaveXml(JXml *sxml,const std::string &place)const;

  void SetMkFirst(word boundfirst,word fluidfirst);
  word GetMkBoundFirst()const{ return(MkBoundFirst); }
  word GetMkFluidFirst()const{ return(MkFluidFirst); }

  void SetPosDomain(tdouble3 posmin,tdouble3 posmax){ Posmin=posmin; Posmax=posmax; }
  tdouble3 GetPosmin()const{ return(Posmin); }
  tdouble3 GetPosmax()const{ return(Posmax); }

  void AddFixed(word mktype,unsigned count){ 
    Add(new JCasePartBlock_Fixed(Properties,mktype,GetBegin(),count));
  }

  void AddMoving(word mktype,unsigned count,unsigned refmotion){ 
    Add(new JCasePartBlock_Moving(Properties,mktype,GetBegin(),count,refmotion));
  }

  void AddFloating(word mktype,unsigned count,double massbody,double masspart
    ,const tdouble3 &center,const tmatrix3d &inertia
    ,const tint3 &translationfree,const tint3 &rotationfree
    ,const tdouble3 &linvelini,const tdouble3 &angvelini
    ,const JLinearValue *linvel,const JLinearValue *angvel
    ,const JLinearValue *linforce,const JLinearValue *angforce)
  { 
    Add(new JCasePartBlock_Floating(Properties,mktype,GetBegin(),count
      ,massbody,masspart,center,inertia,translationfree,rotationfree
      ,linvelini,angvelini,linvel,angvel,linforce,angforce)); 
  }

  void AddFluid(word mktype,unsigned count){ 
    Add(new JCasePartBlock_Fluid(Properties,mktype,GetBegin(),count));
  }

  void SetBlockSize(unsigned pos,unsigned np);

#ifdef JCaseParts_UseProps
  void LoadProperties(const JCaseProperties *props);
#endif

  std::string GetMkList(TpParticles type)const;

  JCaseParts::StSummaryData GetSummaryData()const;
  void GetParticleSummary(std::vector<std::string> &out)const;

  void GetParticlesInfo(std::vector<std::string> &out)const;
  void VisuParticlesInfo()const;

  bool UseImposedFtVel()const;
  bool UseAddedFtForce()const;
};


//##############################################################################
//# JCasePartsGetMk
//##############################################################################
/// \brief Compute Mk value from Id.

class JCasePartsGetMk  : protected JObject
{
private:
  const bool Splitting;
  unsigned MkCount;
  unsigned *MkRange;
  word *MkValue;
  word MkSplitting;

  void Config(const JCaseParts *sparts);
public:
  JCasePartsGetMk(const JCaseParts *sparts,bool splitting);
  ~JCasePartsGetMk();
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


