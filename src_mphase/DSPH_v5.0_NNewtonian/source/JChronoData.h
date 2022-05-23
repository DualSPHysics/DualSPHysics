/*
 <DUALSPHYSICS>  Copyright (c) 2020, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/).

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics.

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>.
*/

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase creada para almacenar la configuracion de los distintos cuerpos y
//:#   sus uniones para usar con Chrono. (05-05-2016).
//:# - Nuevo enumerado para elegir el tipo de Solver (23-12-2019).
//:# - Funcions para incorporar los parametros Young's Modulus y Poisson Ration 
//:#   a los objetos de Chrono (16-01-2020).
//:# - Nuevo objeto link_pulley (24-02-2020).
//:# - Nuevo objeto link_coulombdamping (27-02-2020).
//:# - Nuevos objetos para simular Finite Element Analysis (FEA) en Chrono 
//:#   [JChBodyFEA,JChCable,JChBeamANCF,JChBeamEuler,JChBeamIGA] (24-04-2020).
//:# - Modificacion del tipo de array que contienen informacion de los nodos
//:#   de JChBodyFEA (14-05-2020).
//:# - Se incorporta un array con los deplazamientos de los nodos para objetos 
//:#   FEA. (15-05-2020).
//:# - Se modifica el tipo de datos del array Forces de JChBodyFEA (11-06-2020).
//:# - Se crea una structura StPoints para la gestion de los extremos de
//:#   los objetos JChBodyFEA (15-07-2020).
//:# - Se incorpora la posibilidad de establecer extremos fijos para los
//:#   objetos JChBodyFEA (15-07-2020).
//:# - Permite la simulacion de ChLinks con coeficientes variables de stiffness 
//:#   y damping (04-10-2020).
//:# - Variable density para objetos JChBodyFEA (09-10-2020).
//:# - Se utilizan los parametros DVI y DEM para la definicion de superficies
//:#   en objetos JChBodyFEA (09-10-2020).
//:# - Variable booleana UseCollision para habilitar/deshabilitar colisiones 
//:#   en Chrono (13-10-2020).
//:# - Se renombra las variables Positions y Forces por Pos y Force respectivamente 
//:#   en los objetos JChBodyFEA (13-10-2020).
//:# - Variable Torque para aplicar fuerzas de torsion en los nodos de los
//:#   objetos JChBodyFEA (13-10-2020).
//:# - Variable ImposedFric para indicar que se usara el valor de Kfric
//:#   de un objeto en la colision entre dos objetos (30-10-2020).
//:# - Enumerado TpCollisionFEA en objetos JChronoFEA para indicar el modo
//:#   de colision en objetos deformables (06-11-2020).
//:# - Nueva clase para administrar links tipo JChLinkPointFrame 
//:#   para conectar FEA nodos a objetos rigidos (09-11-2020).
//:# - Funcion ResetForces() para restablecer a 0 las fuerzas del array de los
//:#   objetos FEA (06-12-2020).
//:# - Se incorpora un vector de coeficientes de fuerzas para escalar las fuerzas
//:#   para cada objeto (10-03-2021).
//:#############################################################################

/// \file JChronoData.h \brief Declares the class \ref JChronoData.

#ifndef _JChronoData_
#define _JChronoData_

#include "TypesDef.h"
#include "RunExceptionDef.h"
#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <cstdlib>
#include <cmath>

class JChBody;
class JChLink;

//##############################################################################
//# JChBase
//##############################################################################
/// \brief Defines objects with methods that throw exceptions.
class JChBase
{
protected:
  std::string ClassName;  ///<Name of the class.
  void RunExceptioon(const std::string &srcfile,int srcline
    ,const std::string &classname,const std::string &method
    ,const std::string &msg,const std::string &file="")const;
public:
  JChBase() : ClassName("JChBase"){} ///<Constructor of objects.
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
      tdouble3 vdouble3;   ///<The largest element to set zero.
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
  bool UseFEA;
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
  tfloat3 ScaleForce; ///<Scale the forces of each object in all directions.

  //-Parameters for collisions.
  bool ImposeFric; ///< Indicates that its Sfric and Kfric coefficients will be used in collisions between two bodies.
  float Kfric; ///< Kinetic friction coefficient
  float Sfric; ///< Static friction coefficient
  float Restitu; ///< Restitution coefficient
  //-Extra parameters for collisions using Smooth Contacts.
  float Young; ///< Young modulus
  float Poisson; ///< Poisson ratio

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

  std::string   GetModelFile()   const{ return(ModelFile);  }
  TpModelNormal GetModelNormal() const{ return(ModelNormal); }
  float         GetYoung()       const{ return(Young);       }
  float         GetPoisson()     const{ return(Poisson);     }
  float         GetKfric()       const{ return(Kfric);       }
  float         GetSfric()       const{ return(Sfric);       }
  float         GetRestitu()     const{ return(Restitu);     }
  bool          GetUseFEA()      const{ return(UseFEA);      } 

  bool          GetImposeFric()  const{ return(ImposeFric);  } 
  void          SetImposeFric(const bool i)  { ImposeFric=i; } 

  void SetModel(const std::string &file,TpModelNormal normal){ ModelFile=file; ModelNormal=normal; }
  void SetUseFEA(const bool u){ UseFEA=u; } 
  void SetCollisionData(float kfric,float sfric,float restitu,float young,float poisson);

  
  void    SetScaleForce(tfloat3 s){        ScaleForce=s;}
  tfloat3 GetScaleForce()    const{ return(ScaleForce); }
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
  tfloat3 InputLinearVel;
  tfloat3 InputAngularVel;

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
  void SetInputData(const tfloat3 &face,const tfloat3 &fomegaace)  { InputData=true; InputFace=face; InputFomegaAce=fomegaace; }
  void SetInputDataVel(const tfloat3 &vlin,const tfloat3 &vang)    { InputLinearVel=vlin; InputAngularVel=vang; }
  void SetOutputData(const tdouble3 &center,const tfloat3 &vel,const tfloat3 &omega){ OutputCenter=center; OutputVel=vel; OutputOmega=omega; }

  bool     GetInputData()        const{ return(InputData);        }
  tfloat3  GetInputFace()        const{ return(InputFace);        }
  tfloat3  GetInputFomegaAce()   const{ return(InputFomegaAce);   }
  tfloat3  GetInputLinearVel()   const{ return(InputLinearVel);   }
  tfloat3  GetInputAngularVel()  const{ return(InputAngularVel);  }
  tdouble3 GetOutputCenter()     const{ return(OutputCenter);     }
  tfloat3  GetOutputVel()        const{ return(OutputVel);        }
  tfloat3  GetOutputOmega()      const{ return(OutputOmega);      }
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

  TpMotion GetMotionType()    const{ return(MotionType);    } 
  tdouble3 GetMotionSimple()  const{ return(MotionSimple);  }
  tmatrix4d GetMotionMatrix() const{ return(MotionMatrix);  }
  double GetMotionDt()        const{ return(MotionDt);      }
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
  typedef enum{ LK_Hinge,LK_Spheric,LK_PointLine,LK_LinearSpring,LK_CoulombDamping,LK_Pulley,LK_PointFrame }TpLink;
  
  /// Structure with parameters to create VTK of spring.
  typedef struct StrSaveSpring{
    float radius;    ///<Spring radius.
    float length;    ///<Length for each revolution.
    int nside;       ///<Number of sections for each revolution.

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

  bool VariableK; ///<Indicates if it is using a variable stiffness coefficient along time.
  bool VariableC; ///<Indicates if it is using a variable damping coefficient along time.

  double Stiffness;	///<Stiffness.
  double Damping;	///<Damping.

  tdouble3 RotPoint;	///<Point for rotation.
  tdouble3 RotVector;	///<Vector for rotation.  

  tdouble3 Pointfb0;	///<Point from fb0
  tdouble3 Pointfb1;	///<Point from fb1

  tdouble3 SlidingVector; ///<Vector direction for sliding axis (for pointline).  
  tdouble3 RotVector2;    ///<Second vector for rotation (for pointline).  

  double RestLength;      ///<Rest length for spring.
  StSaveSpring SvSpring;  ///<Configuration to save vtk spring.
  double CoulombDamping;  ///<Coulomb damping value. Not applied when it is zero. (default=0).

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
  
  bool GetVariableK()  const{ return(VariableK); } 
  bool GetVariableC()  const{ return(VariableC);   }

  void SetStiffness(double v){ Stiffness=v; } 
  void SetDamping  (double v){ Damping=v;   } 
  
  void SetVariableK(bool v){ VariableK=v; } 
  void SetVariableC(bool v){ VariableC=v;   } 
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

  tdouble3 GetRotPoint() const{ return(RotPoint);   } 
  tdouble3 GetRotVector()const{ return(RotVector);  }  

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
  JChLinkSpheric(std::string name,unsigned idbody1,unsigned idbody2);
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
  JChLinkPointLine(std::string name,unsigned idbody1,unsigned idbody2);
  JChLinkPointLine(const JChLinkPointLine &src);

  tdouble3 GetSlidingVector() const{ return(SlidingVector); }  
  tdouble3 GetRotPoint()      const{ return(RotPoint);      } 
  tdouble3 GetRotVector()     const{ return(RotVector);     }  
  tdouble3 GetRotVector2()    const{ return(RotVector2);    }  

  void SetSlidingVector(tdouble3 v){  SlidingVector=v;  }  
  void SetRotPoint (tdouble3 v){      RotPoint=v;       }
  void SetRotVector(tdouble3 v){      RotVector=v;      }  
  void SetRotVector2(tdouble3 v){     RotVector2=v;     }  
};

//##############################################################################
//# JChLinearSpring
//##############################################################################
/// \brief Manages the info of a link.
class JChLinkLinearSpring : public JChLink
{
public:
  JChLinkLinearSpring(std::string name,unsigned idbody1,unsigned idbody2);
  JChLinkLinearSpring(const JChLinkLinearSpring &src);

  tdouble3 GetPointfb0()const{ return(Pointfb0); }
  tdouble3 GetPointfb1()const{ return(Pointfb1); }
  double GetRestLength()const{ return(RestLength); }
  StSaveSpring GetSvSpring()const{ return(SvSpring); }

  void SetPointfb0  (tdouble3 v){ Pointfb0 =v; }
  void SetPointfb1  (tdouble3 v){ Pointfb1 =v; }
  void SetRestLength(double   v){ RestLength=v; }
  void SetSvSpring  (const StSaveSpring &v){ SvSpring=v; }
};

//##############################################################################
//# JChLinkCoulombDamping
//##############################################################################
/// \brief Manages the info of a link.
class JChLinkCoulombDamping : public JChLink
{
public:
  JChLinkCoulombDamping(std::string name,unsigned idbody1,unsigned idbody2);
  JChLinkCoulombDamping(const JChLinkCoulombDamping &src);

  tdouble3 GetPointfb0()const{ return(Pointfb0); }
  tdouble3 GetPointfb1()const{ return(Pointfb1); }
  double GetRestLength()const{ return(RestLength); }
  double GetCoulombDamping()const{ return(CoulombDamping); }
  StSaveSpring GetSvSpring()const{ return(SvSpring); }

  void SetPointfb0  (tdouble3 v){ Pointfb0 =v; }
  void SetPointfb1  (tdouble3 v){ Pointfb1 =v; }
  void SetRestLength(double   v){ RestLength=v; }
  void SetCoulombDamping(double v){ CoulombDamping=v; }
  void SetSvSpring  (const StSaveSpring &v){ SvSpring=v; }

};

//##############################################################################
//# JChLinkPulley
//##############################################################################
/// \brief Manages the info of a link.
class JChLinkPulley : public JChLink
{
private:
  float Radius;
  float Radius2;
  void CopyFrom(const JChLinkPulley &src);
public:
  JChLinkPulley(std::string name,unsigned idbody1,unsigned idbody2);
  JChLinkPulley(const JChLinkPulley &src);

  tdouble3 GetRotPoint() const{ return(RotPoint);  } 
  tdouble3 GetRotVector()const{ return(RotVector); }  
  float    GetRadius()   const{ return(Radius);    } 
  float    GetRadius2()  const{ return(Radius2);   }  

  void SetRotPoint  (tdouble3 v){ RotPoint =v; }
  void SetRotVector (tdouble3 v){ RotVector=v; }  
  void SetRadius    (float r)   { Radius=r;    }
  void SetRadius2   (float r)   { Radius2=r;   }  
};

//##############################################################################
//# JChronoData
//##############################################################################
/// \brief Manages the info of Chrono from the input XML file.
class JChronoData : protected JChBase
{
public:
  typedef enum{ NSC,SMC }               TpContactMethod;
  typedef enum{ BB=0,APGD=1,APGDREF=2 } DSolverType;     //Allowed Solvers

private:
  std::string Mode;         ///<Chrono execution mode
  //unsigned MaxIter;         ///<Indicates the maximun number of iterations for the solver
  std::string DataDir;
  bool UseNSCChrono;
  std::vector<JChBody*> LisBody;
  std::vector<JChLink*> LisLink;
  double Dp;
  float CollisionDp;
  unsigned Solver;
  int OmpThreads;
  bool UseFEA;                  //<vs_chronoo_fea>
  TpContactMethod ContactMethod;
  bool UseVariableCoeff;
  bool UseCollision;

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

  void SetUseCollision(bool c){ UseCollision=c; }
  bool GetUseCollision()const{ return(UseCollision); }

  void SetUseNSCChrono(bool useNSCchrono){ UseNSCChrono=useNSCchrono; }
  bool GetUseNSCChrono()const{ return(UseNSCChrono); }

  void SetDp(double v){ Dp=v; }
  double GetDp()const{ return(Dp); }

  void SetCollisionDp(float v){ CollisionDp=v; }
  float GetCollisionDp()const{ return(CollisionDp); }

  void SetSolver(unsigned s){ Solver=s; }
  unsigned GetSolver()const{ return(Solver); }

  void SetOmpThreads(unsigned th){ OmpThreads=th; }
  int GetOmpThreads()const{ return(OmpThreads); }

  void SetMode(std::string m){ Mode=m; }
  std::string GetMode()const{ return(Mode); }

  //void SetMaxIter(unsigned i){ MaxIter=i; }
  //unsigned GetMaxIter()const{ return(MaxIter); }

  void SetUseFEA(bool e){ UseFEA=e; }          //<vs_chronoo_fea>
  bool GetUseFEA()const{ return(UseFEA); }     //<vs_chronoo_fea>

  void SetContactMethod(TpContactMethod c){ ContactMethod=c; }
  TpContactMethod GetContactMethod()const{ return(ContactMethod); }

  std::string SolverToStr()const;
  std::string ContactMethodToStr()const;

  std::string CheckAddBodyError(unsigned idb,std::string idname,word mkbound)const;

  JChBodyFloating* AddBodyFloating(unsigned idb,std::string idname,word mkbound,std::string fileinfo="");
  JChBodyMoving*   AddBodyMoving  (unsigned idb,std::string idname,word mkbound,double mass,std::string fileinfo="");
  JChBodyFixed*    AddBodyFixed   (unsigned idb,std::string idname,word mkbound,std::string fileinfo="");

  JChLinkHinge*          AddLinkHinge         (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkPulley*         AddLinkPulley        (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkSpheric*        AddLinkSpheric       (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkPointLine*      AddLinkPointLine     (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkLinearSpring*   AddLinkLinearSpring  (std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");
  JChLinkCoulombDamping* AddLinkCoulombDamping(std::string name,unsigned idbody1,unsigned idbody2,std::string fileinfo="");

  unsigned GetBodyCount()const{ return(unsigned(LisBody.size())); }
  unsigned GetLinkCount()const{ return(unsigned(LisLink.size())); }

  unsigned BodyIndexById(unsigned idb)const;
  unsigned BodyIndexByMkBound(word mkbound)const;
  unsigned BodyIndexByName(const std::string &name)const;
  unsigned LinkIndexByName(const std::string &name)const;

  const JChBody* GetBody(unsigned ipos)const;
  const JChBody* GetBodyByMk(word mkbound)const;
  JChLink* GetLink(unsigned ipos)const;
  unsigned GetPosLink(const JChLink *link)const;

  std::vector<const JChLink*> GetLinkByTp(JChLink::TpLink type)const;
  bool BodyBelongsLink(unsigned idbody1,JChLink::TpLink type)const;

  const JChBodyFloating* GetBodyFloating(word mkbound)const;
  const JChBodyMoving*   GetBodyMoving  (word mkbound)const;
  const JChBodyFixed*    GetBodyFixed   (word mkbound)const;

  bool GetUseVariableCoeff()const{return(UseVariableCoeff);}
  void SetUseVariableCoeff(bool v){UseVariableCoeff=v;}

};
#endif