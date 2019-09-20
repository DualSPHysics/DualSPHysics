/*
 <DUALSPHYSICS>  Copyright (c) 2019, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

//#############################################################################
//# Cambios:
//# =========
//# - Clase creada para almacenar la configuracion de los distintos cuerpos y
//#   sus uniones para usar con Chrono. (05-05-2016)
//#############################################################################

/// \file JChronoData.h \brief Declares the class \ref JChronoData.

#ifndef _JChronoData_
#define _JChronoData_

#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>
#include "TypesDef.h"

//-Defines for normal exceptions.
#ifndef Run_Exceptioon
#define Run_Exceptioon(msg) RunExceptioon(__FILE__,__LINE__,ClassName,__func__,msg)
#endif
#ifndef Run_ExceptioonFile
#define Run_ExceptioonFile(msg,file) RunExceptioon(__FILE__,__LINE__,ClassName,__func__,msg,file)
#endif

class JChBody;
class JChLink;

//##############################################################################
//# JChBase
//##############################################################################
/// \brief Defines objects with methods that throw exceptions.
class JChBase
{
protected:
  std::string ClassName;           ///<Name of the class.
  void RunExceptioon(const std::string &srcfile,int srcline
    ,const std::string &classname,const std::string &method
    ,const std::string &msg,const std::string &file="")const;
public:  
  JChBase():ClassName("JChBase"){} ///<Constructor of objects.
};

//##############################################################################
//# JChValues
//##############################################################################
/// \brief Manages the info stored in values section.
class JChValues : protected JChBase
{
public:
  ///<Types of value.
  typedef enum{ TP_Text,TP_Int,TP_Uint,TP_Double,TP_Int3,TP_Uint3,TP_Double3 }TpValue; 

  ///Structure that describes the information of a value.
  typedef struct{
    std::string name;
    TpValue type;
    std::string vtext;
    union{
      int vint;
      unsigned vuint;
      double vdouble;
      tint3 vint3;
      tuint3 vuint3;
      tfloat3 vfloat3;
      tdouble3 vdouble3;   //-The largest element to set zero.
    };
  }StValue;

  static std::string TypeToStr(TpValue type);

private:
  std::vector<StValue> LisValue;

  StValue GetValueVoid(TpValue type=TP_Int,const std::string name="")const;
  void AddValue(const StValue &va);
  unsigned CheckValueType(const std::string &name,TpValue type,bool optional)const;

public:
  JChValues(); 
  ~JChValues();
  void Reset();

  unsigned GetCount()const{ return(unsigned(LisValue.size())); }
  unsigned IndexByName(const std::string &name)const;

  void AddValueStr    (const std::string &name,std::string v);
  void AddValueInt    (const std::string &name,int         v);
  void AddValueUint   (const std::string &name,unsigned    v);
  void AddValueDouble (const std::string &name,double      v);
  void AddValueInt3   (const std::string &name,tint3       v);
  void AddValueUint3  (const std::string &name,tuint3      v);
  void AddValueDouble3(const std::string &name,tdouble3    v);

  bool ExistsValue(const std::string &name)const;
  bool ExistsValue(const std::string &name,TpValue type)const;

  std::string GetValueStr    (const std::string &name,bool optional=false,std::string v=""         )const;
  int         GetValueInt    (const std::string &name,bool optional=false,int         v=0          )const;
  unsigned    GetValueUint   (const std::string &name,bool optional=false,unsigned    v=0          )const;
  double      GetValueDouble (const std::string &name,bool optional=false,double      v=0          )const;
  tint3       GetValueInt3   (const std::string &name,bool optional=false,tint3       v=TInt3(0)   )const;
  tuint3      GetValueUint3  (const std::string &name,bool optional=false,tuint3      v=TUint3(0)  )const;
  tdouble3    GetValueDouble3(const std::string &name,bool optional=false,tdouble3    v=TDouble3(0))const;

  const StValue* GetValue(unsigned ipos)const;

};

//##############################################################################
//# JChObject
//##############################################################################
/// \brief Defines objects with methods that throw exceptions.
class JChObject : protected JChBase
{
protected:
  JChValues Values;

public:  
  JChObject(){ ClassName="JChBase"; } ///<Constructor of objects.

  JChValues* GetValuesPtr()const{ return((JChValues*)&Values); }

  bool ExistsValue(const std::string &name)const{ return(Values.ExistsValue(name)); }
  bool ExistsValue(const std::string &name,JChValues::TpValue type)const{ return(Values.ExistsValue(name,type)); }

  std::string GetValueStr    (const std::string &name,bool optional=false,std::string v=""         )const{ return(GetValueStr    (name,optional,v)); }
  int         GetValueInt    (const std::string &name,bool optional=false,int         v=0          )const{ return(GetValueInt    (name,optional,v)); }
  unsigned    GetValueUint   (const std::string &name,bool optional=false,unsigned    v=0          )const{ return(GetValueUint   (name,optional,v)); }
  double      GetValueDouble (const std::string &name,bool optional=false,double      v=0          )const{ return(GetValueDouble (name,optional,v)); }
  tint3       GetValueInt3   (const std::string &name,bool optional=false,tint3       v=TInt3(0)   )const{ return(GetValueInt3   (name,optional,v)); }
  tuint3      GetValueUint3  (const std::string &name,bool optional=false,tuint3      v=TUint3(0)  )const{ return(GetValueUint3  (name,optional,v)); }
  tdouble3    GetValueDouble3(const std::string &name,bool optional=false,tdouble3    v=TDouble3(0))const{ return(GetValueDouble3(name,optional,v)); }
};


//##############################################################################
//# JChBody
//##############################################################################
/// \brief Manages the info of each body.
class JChBody : public JChObject
{
public:
  ///<Types of body.
  typedef enum{ BD_Floating,BD_Moving,BD_Fixed }TpBody; 
  typedef enum{ NorNull,NorOriginal,NorInvert,NorTwoFace }TpModelNormal; 

  static std::string TypeToStr(TpBody type);
  static std::string NormalToStr(TpModelNormal tnor);

private:
  std::vector<JChLink*> LinkRefs;

protected:
  double Mass;
  tdouble3 Center;
  tmatrix3d Inertia;
  bool MotionFree;
  tint3 TranslationFree;
  tint3 RotationFree;
  tfloat3 LinearVelini;
  tfloat3 AngularVelini;
  std::string ModelFile;
  TpModelNormal ModelNormal;

  //DVI parameters
  float Kfric;
  float Restitu;
  //Extra DEM parameters
  //float Young;    //(DEM)
  //float Poisson;  //(DEM)

  void CopyFrom(const JChBody &src);

public:
  const unsigned Idb;
  const std::string IdName;
  const TpBody Type;     
  const word MkBound;

  JChBody(unsigned idb,std::string idname,TpBody type,word mkbound); 
  virtual ~JChBody();
  virtual void Reset();

  void ResetRefs(){ LinkRefs.clear(); }
  void AddLinkRef(const JChLink* link);
  unsigned GetLinkRefCount()const{ return(unsigned(LinkRefs.size())); };
  const JChLink* GetLinkRef(unsigned ipos)const;

  double      GetMass()       const{ return(Mass);        } 
  tdouble3    GetCenter()     const{ return(Center);      } 
  tmatrix3d   GetInertia()    const{ return(Inertia);     } 

  bool        GetMotionFree()     const{ return(MotionFree);      }
  tint3       GetTranslationFree()const{ return(TranslationFree); }
  tint3       GetRotationFree()   const{ return(RotationFree);    }

  tfloat3     GetLinearVelini()   const{ return(LinearVelini);    }
  tfloat3     GetAngularVelini()  const{ return(AngularVelini);   }

  std::string GetModelFile()  const{ return(ModelFile);   }
  TpModelNormal GetModelNormal()const{ return(ModelNormal); }
//float       GetYoung()      const{ return(Young);       }  //(DEM)
//float       GetPoisson()    const{ return(Poisson);     }  //(DEM)
  float       GetKfric()      const{ return(Kfric);       }
  float       GetRestitu()    const{ return(Restitu);     }

  void SetModel(const std::string &file,TpModelNormal normal){ ModelFile=file; ModelNormal=normal; } 
  void SetDVIData(float kfric,float restitu);
//void SetDEMData(float kfric,float restitu,float young,float poisson);  //(DEM)
};

//##############################################################################
//# JChBodyFloating
//##############################################################################
/// \brief Manages the info of a floating body.
class JChBodyFloating : public JChBody
{
protected:
  bool InputData;
  tfloat3 InputFace;
  tfloat3 InputFomegaAce;

  tdouble3 OutputCenter;
  tfloat3 OutputVel;
  tfloat3 OutputOmega;

public:
  JChBodyFloating(unsigned idb,std::string idname,word mkbound);
  JChBodyFloating(const JChBodyFloating &src);
  void Reset();
  void SetFloatingData(double mass,tdouble3 center,tmatrix3d inertia
    ,tint3 translationfree,tint3 rotationfree,tfloat3 linvelini,tfloat3 angvelini);

  void ResetInputData(){ InputData=false; }
  void SetInputData(const tfloat3 &face,const tfloat3 &fomegaace){ InputData=true; InputFace=face; InputFomegaAce=fomegaace; }
  void SetOutputData(const tdouble3 &center,const tfloat3 &vel,const tfloat3 &omega){ OutputCenter=center; OutputVel=vel; OutputOmega=omega; }

  bool     GetInputData()      const{ return(InputData);       }
  tfloat3  GetInputFace()      const{ return(InputFace);       }
  tfloat3  GetInputFomegaAce() const{ return(InputFomegaAce);  }
  tdouble3 GetOutputCenter()   const{ return(OutputCenter);    }
  tfloat3  GetOutputVel()      const{ return(OutputVel);       }
  tfloat3  GetOutputOmega()    const{ return(OutputOmega);     }
};

//##############################################################################
//# JChBodyMoving
//##############################################################################
/// \brief Manages the info of a moving body.
class JChBodyMoving : public JChBody
{
public:
  typedef enum{ MV_None,MV_Simple,MV_Matrix }TpMotion; 

protected:
  TpMotion MotionType;    ///<Type of motion: linear or rotational.
  tdouble3 MotionSimple;  ///<Linear displacement during dt.
  tmatrix4d MotionMatrix; ///<Trasformation matrix for displacement during dt.
  double MotionDt;        ///<Stores the total dt even a predictor-corrector is used.

public:
  JChBodyMoving(unsigned idb,std::string idname,word mkbound,double mass);
  JChBodyMoving(const JChBodyMoving &src);
  void Reset();
  void SetInitialCenter(const tdouble3  &center){ Center=center; }

  void ResetMotion(){ MotionType=MV_None; }
  void SetMotionSimple(double dt,const tdouble3  &msimple){ MotionType=MV_Simple; MotionDt=dt; MotionSimple=msimple; }
  void SetMotionMatrix(double dt,const tmatrix4d &mmatrix){ MotionType=MV_Matrix; MotionDt=dt; MotionMatrix=mmatrix; }

  TpMotion GetMotionType()const{ return(MotionType); }
  tdouble3 GetMotionSimple()const{ return(MotionSimple); }
  tmatrix4d GetMotionMatrix()const{ return(MotionMatrix); }
  double GetMotionDt()const{ return(MotionDt); }
};

//##############################################################################
//# JChBodyFixed
//##############################################################################
/// \brief Manages the info of a fixed body.
class JChBodyFixed : public JChBody
{
public:
  JChBodyFixed(unsigned idb,std::string idname,word mkbound);
  JChBodyFixed(const JChBodyFixed &src);
  void Reset();
};

//##############################################################################
//# JChLink
//##############################################################################
/// \brief Manages the info of each link between bodies.
class JChLink : public JChObject
{
public:
  ///<Types of link.
  typedef enum{ LK_Hinge, LK_Spheric, LK_PointLine, LK_LinearSpring }TpLink;


  /// Structure with parameters to create VTK of spring.
  typedef struct StrSaveSpring{
    float radius;    //-Spring radius.
    float length;    //-Length for each revolution.
    int nside;       //-Number of sections for each revolution.

    StrSaveSpring(){ reset(); }
    StrSaveSpring(float xradius,float xlength,int xnside)
    { 
      reset(); radius=xradius; length=xlength; nside=xnside;
    }
    void reset(){ 
      radius=3;
      length=1.f;
      nside=16;
    }
  }StSaveSpring;


  static std::string TypeToStr(TpLink type);

protected:
  std::vector<JChBody*> BodyRefs;

  double Stiffness;   //-Stiffness.
  double Damping;     //-Damping.

  tdouble3 RotPoint;  //-Point for rotation.
  tdouble3 RotVector; //-Vector for rotation.  

  tdouble3 Pointfb0;  //-Point from fb0
  tdouble3 Pointfb1;  //-Point from fb1

  tdouble3 SlidingVector; //-Vector direction for sliding axis (for pointline).  
  tdouble3 RotVector2;    //-Second vector for rotation (for pointline).  

  double RestLength;  //-Rest length for spring.
  StSaveSpring SvSpring; //-Configuration to save vtk spring.
  double CoulombDamping; //-Coulomb damping value. Not applied when it is zero. (default=0).

protected:
  void CopyFrom(const JChLink &src);

public:
  const unsigned IdBody1;
  const unsigned IdBody2;
  const std::string Name;
  const TpLink Type;     

  JChLink(std::string name,TpLink type,unsigned idbody1,unsigned idbody2); 
  virtual ~JChLink();
  void Reset();

  void ResetRefs(){ BodyRefs.clear(); }
  void AddBodyRef(const JChBody* body);
  unsigned GetBodyRefCount()const{ return(unsigned(BodyRefs.size())); };
  const JChBody* GetBodyRef(unsigned ipos)const;
  
  double GetStiffness()const{ return(Stiffness); } 
  double GetDamping()  const{ return(Damping);   }

  void SetStiffness(double v){ Stiffness=v; } 
  void SetDamping  (double v){ Damping  =v; } 
};

//##############################################################################
//# JChLinkHinge
//##############################################################################
/// \brief Manages the info of a link.
class JChLinkHinge : public JChLink
{
public:
  JChLinkHinge(std::string name,unsigned idbody1,unsigned idbody2);
  JChLinkHinge(const JChLinkHinge &src);

  tdouble3 GetRotPoint() const{ return(RotPoint); } 
  tdouble3 GetRotVector()const{ return(RotVector); }  

  void SetRotPoint (tdouble3 v){ RotPoint =v; }
  void SetRotVector(tdouble3 v){ RotVector=v; }  
};

//##############################################################################
//# JChLinkSpheric
//##############################################################################
/// \brief Manages the info of a link.
class JChLinkSpheric : public JChLink
{
public:
  JChLinkSpheric(std::string name, unsigned idbody1, unsigned idbody2);
  JChLinkSpheric(const JChLinkSpheric &src);

  tdouble3 GetRotPoint()const{ return(RotPoint); } 

  void SetRotPoint(tdouble3 v){ RotPoint=v; }
};

//##############################################################################
//# JChPointLine
//##############################################################################
/// \brief Manages the info of a link.
class JChLinkPointLine : public JChLink
{
public:
  JChLinkPointLine(std::string name, unsigned idbody1, unsigned idbody2);
  JChLinkPointLine(const JChLinkPointLine &src);

  tdouble3 GetSlidingVector()const{ return(SlidingVector); }  
  tdouble3 GetRotPoint() const{ return(RotPoint); } 
  tdouble3 GetRotVector()const{ return(RotVector); }  
  tdouble3 GetRotVector2()const{ return(RotVector2); }  

  void SetSlidingVector(tdouble3 v){ SlidingVector=v; }  
  void SetRotPoint (tdouble3 v){ RotPoint=v; }
  void SetRotVector(tdouble3 v){ RotVector=v; }  
  void SetRotVector2(tdouble3 v){ RotVector2=v; }  
};

//##############################################################################
//# JChLinearSpring
//##############################################################################
/// \brief Manages the info of a link.
class JChLinkLinearSpring : public JChLink
{
public:
  JChLinkLinearSpring(std::string name, unsigned idbody1, unsigned idbody2);
  JChLinkLinearSpring(const JChLinkLinearSpring &src);

  tdouble3 GetPointfb0()const{ return(Pointfb0); }
  tdouble3 GetPointfb1()const{ return(Pointfb1); }
  double GetRestLength()const{ return(RestLength); }
  StSaveSpring GetSvSpring()const{ return(SvSpring); }
  double GetCoulombDamping()const{ return(CoulombDamping); }

  void SetPointfb0  (tdouble3 v){ Pointfb0 =v; }
  void SetPointfb1  (tdouble3 v){ Pointfb1 =v; }
  void SetRestLength(double   v){ RestLength=v; }
  void SetSvSpring  (const StSaveSpring &v){ SvSpring=v; }
  void SetCoulombDamping(double v){ CoulombDamping=v; }
};

//##############################################################################
//# JChronoData
//##############################################################################
/// \brief Manages the info of Chrono from the input XML file.
class JChronoData : protected JChBase
{
private:
  std::string DataDir;
  bool UseNSCChrono;
  std::vector<JChBody*> LisBody;
  std::vector<JChLink*> LisLink;
  double Dp;
  float CollisionDp;

public:
  JChronoData();
  JChronoData(const JChronoData &src);
  ~JChronoData();
  JChronoData& operator=(const JChronoData &src);
  void Reset();
  void Prepare();
  void CheckData();

  void SetDataDir(const std::string &datadir){ DataDir=datadir; }
  std::string GetDataDir()const{ return(DataDir); }

  void SetUseNSCChrono(bool useNSCchrono){ UseNSCChrono=useNSCchrono; }
  bool GetUseNSCChrono()const{ return(UseNSCChrono); }

  void SetDp(double v){ Dp=v; }
  double GetDp()const{ return(Dp); }

  void SetCollisionDp(float v){ CollisionDp=v; }
  float GetCollisionDp()const{ return(CollisionDp); }

  std::string CheckAddBodyError(unsigned idb,std::string idname,word mkbound)const;
  
  JChBodyFloating* AddBodyFloating(unsigned idb,std::string idname,word mkbound,std::string fileinfo="");
  JChBodyMoving*   AddBodyMoving  (unsigned idb,std::string idname,word mkbound,double mass,std::string fileinfo="");
  JChBodyFixed*    AddBodyFixed   (unsigned idb,std::string idname,word mkbound,std::string fileinfo="");

  JChLinkHinge*        AddLinkHinge       (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkSpheric*      AddLinkSpheric     (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkPointLine*    AddLinkPointLine   (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkLinearSpring* AddLinkLinearSpring(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  
  unsigned GetBodyCount()const{ return(unsigned(LisBody.size())); }
  unsigned GetLinkCount()const{ return(unsigned(LisLink.size())); }

  unsigned BodyIndexById(unsigned idb)const;
  unsigned BodyIndexByMkBound(word mkbound)const;
  unsigned BodyIndexByName(const std::string &name)const;
  unsigned LinkIndexByName(const std::string &name)const;

  const JChBody* GetBody(unsigned ipos)const;
  const JChBody* GetBodyByMk(word mkbound)const;
  const JChLink* GetLink(unsigned ipos)const;

  
  const JChBodyFloating* GetBodyFloating(word mkbound)const;
  const JChBodyMoving*   GetBodyMoving  (word mkbound)const;
  const JChBodyFixed*    GetBodyFixed   (word mkbound)const;
};



#endif




