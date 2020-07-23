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
//:# - Creacion de clase para gestionar informacion relativa a datos de simulacion. (23-03-2018)
//:# - Metodos para obtener informacion de particulas. (25-04-2018)
//:# - Establece PERI_Unknown por defecto. (27-04-2018)
//:# - Improved definition of the periodic conditions. (27-04-2018)
//:# - Incluye informacion de Symmetry. (13-05-2018)
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:#############################################################################

/// \file JPartDataHead.h \brief Declares the class \ref JPartDataHead.

#ifndef _JPartDataHead_
#define _JPartDataHead_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"
#include "JParticlesDef.h"
#include "JPeriodicDef.h"

//##############################################################################
//# JPartDataHeadMkBlock
//##############################################################################
/// \brief Stores the info of a block of particles.
class JPartDataHeadMkBlock
{
public:
  bool Bound;         ///<Indicates whether a particle is boundary or not.
  TpParticles Type;   ///<Type of particle.
  unsigned Mk;        ///<Absolute label.
  unsigned MkType;    ///<Label of block fluid or bound.
  unsigned Begin;     ///<Id of the first particle of the block.
  unsigned Count;     ///<Number of particles.

  JPartDataHeadMkBlock(TpParticles type,unsigned mk,unsigned mktype,unsigned begin,unsigned count)
    :Bound(IsBound(type)),Type(type),Mk(mk),MkType(mktype),Begin(begin),Count(count){ }
};

//##############################################################################
//# JPartDataHead
//##############################################################################
/// \brief Manages the information of simulation data.

class JPartDataHead : protected JObject
{
public:
  ///Types of viscosity treatment.
  typedef enum{ 
    VISCO_LaminarSPS=2,        ///<Laminar viscosity and Sub-Partice Scale Turbulence.
    VISCO_Artificial=1,        ///<Artificial viscosity.
    VISCO_None=0 
  }TpVisco;            

private:
  static const unsigned FmtVersionDef=180324;    ///<Version de formato by default. Version of format by default.
  unsigned FmtVersion;    ///<Version de formato. Version of format.

  std::string DirData;

  //-General variables.
  std::string AppName;
  std::string Date;
  std::string RunCode;

  std::string CaseName;
  bool Data2d;           ///<Toggles 2D simulation.
  double Data2dPosY;     ///<Y value in 2D simulations.
  unsigned Npiece;       ///<Total number of file pieces.
  unsigned FirstPart;    ///<First PART number.

  tdouble3 CasePosMin;   ///<Lower particle limit of the case in the initial instant.
  tdouble3 CasePosMax;   ///<Upper particle limit of the case in the initial instant.

  double Dp;
  double H;
  double B;
  double Gamma;
  double RhopZero;
  double MassBound;  
  double MassFluid;
  tfloat3 Gravity;

  bool NpDynamic;        ///<CaseNp can increase.
  bool ReuseIds;         ///<Id of particles excluded values ​​are reused.

  tdouble3 MapPosMin;    ///<Lower limit of simulation + edge (KernelSize) if periodic conditions.
  tdouble3 MapPosMax;    ///<Upper limit of simulation + edge (KernelSize) if periodic conditions.

  TpPeri PeriMode;
  tdouble3 PeriXinc;     ///<Value that is added at the outer limit to modify the position.
  tdouble3 PeriYinc;     ///<Value that is added at the outer limit to modify the position.
  tdouble3 PeriZinc;     ///<Value that is added at the outer limit to modify the position.

  TpVisco ViscoType;     ///<Viscosity type: Artificial,... 
  float ViscoValue;      ///<Viscosity value. 
  float ViscoBoundFactor;///<For boundary interaction use ViscoValue*ViscoBoundFactor.

  bool Symmetry;         ///<Use of symmetry according plane y=0.

  bool Splitting;        ///<Use of Splitting.

  //-Mk blocks data.
  std::vector<JPartDataHeadMkBlock> MkList;
  unsigned MkListSize;   ///<Total number of Mk blocks.
  unsigned MkListFixed;  ///<Number of Mk blocks of fixed type.
  unsigned MkListMoving; ///<Number of Mk blocks of moving type.
  unsigned MkListFloat;  ///<Number of Mk blocks of floating type.
  unsigned MkListBound;  ///<Number of Mk blocks of boundary types. MkListBound=MkListFixed+MkListMoving+MkListFloat
  unsigned MkListFluid;  ///<Number of Mk blocks of fluid type.

  unsigned MkBoundFirst; ///<First Mk for boundary blocks (Mk=MkBound+MkBoundFirst).
  unsigned MkFluidFirst; ///<First Mk for fluid blocks (Mk=MkFluid+MkFluidFirst).

  //-Number of particles.
  unsigned CaseNp;       ///<Number of total particles of initial PART.  
  unsigned CaseNfixed;   ///<Number of fixed boundary particles. 
  unsigned CaseNmoving;  ///<Number of moving boundary particles. 
  unsigned CaseNfloat;   ///<Number of floating boundary particles. 
  unsigned CaseNfluid;   ///<Number of fluid particles (including the excluded ones). 

  void UptadeMkNumbers();

public:
  JPartDataHead();
  ~JPartDataHead();
  void Reset();
  void ConfigBasic(std::string runcode,std::string appname
    ,std::string casename,tdouble3 caseposmin,tdouble3 caseposmax
    ,bool data2d,double data2dposy,unsigned npieces,unsigned firstpart);
  void ConfigParticles(TpParticles type,unsigned mk,unsigned mktype,unsigned begin,unsigned count);
  void ConfigCtes(double dp,double h,double b,double rhop0,double gamma
    ,double massbound,double massfluid,tfloat3 gravity);
  void ConfigSimNp(bool npdynamic=false,bool reuseids=false);
  void ConfigSimMap(tdouble3 mapposmin,tdouble3 mapposmax);
  void ConfigSimPeri(TpPeri tperi,tdouble3 perixinc,tdouble3 periyinc,tdouble3 perizinc);
  void ConfigVisco(JPartDataHead::TpVisco type,float value,float boundfactor);
  void ConfigSymmetry(bool symmetry);
  void ConfigSplitting(bool splitting);

  static std::string GetFileName(std::string dir="");

  void LoadFile(std::string dir);
  void SaveFile(std::string dir);

  void GetParticlesInfo(std::vector<std::string> &out)const;
  void VisuParticlesInfo()const;

//-Methods for Mk blocks.
//------------------------
  unsigned MkBlockCount()const{ return(MkListSize); };
  const JPartDataHeadMkBlock& Mkblock(unsigned c)const{ return(MkList[c]); }

  unsigned GetMkBlockByMk(unsigned mk)const;
  unsigned GetMkBlockByMkBound(unsigned mkbound)const;
  unsigned GetMkBlockByMkFluid(unsigned mkfluid)const;
  inline unsigned GetMkBlockById(unsigned id)const;

  unsigned GetMkBoundFirst()const{ return(MkBoundFirst); }
  unsigned GetMkFluidFirst()const{ return(MkFluidFirst); }

//-Returns general values.
//--------------------------
  std::string GetAppName()   const{ return(AppName);    };
  std::string GetDate()      const{ return(Date);       };
  std::string GetRunCode()   const{ return(RunCode);    };
  std::string GetCaseName()  const{ return(CaseName);   };
  bool        GetData2d()    const{ return(Data2d);     };
  double      GetData2dPosY()const{ return(Data2dPosY); };
  unsigned    GetNpiece()    const{ return(Npiece);     };
  unsigned    GetFirstPart() const{ return(FirstPart);  };


  unsigned GetCaseNp()     const{ return(CaseNp);      };
  unsigned GetCaseNfixed() const{ return(CaseNfixed);  };
  unsigned GetCaseNmoving()const{ return(CaseNmoving); };
  unsigned GetCaseNfloat() const{ return(CaseNfloat);  };
  unsigned GetCaseNfluid() const{ return(CaseNfluid);  };
  tdouble3 GetCasePosMin() const{ return(CasePosMin);  };
  tdouble3 GetCasePosMax() const{ return(CasePosMax);  };

  bool GetNpDynamic()const{ return(NpDynamic); };
  bool GetReuseIds() const{ return(ReuseIds);  };

  tdouble3 GetMapPosMin()const{ return(MapPosMin); };
  tdouble3 GetMapPosMax()const{ return(MapPosMax); };

  TpPeri GetPeriMode()const{ return(PeriMode); };
  tdouble3 GetPeriXinc()const{ return(PeriXinc);   };
  tdouble3 GetPeriYinc()const{ return(PeriYinc);   };
  tdouble3 GetPeriZinc()const{ return(PeriZinc);   };

  TpVisco GetViscoType()       const{ return(ViscoType);        };
  float   GetViscoValue()      const{ return(ViscoValue);       };
  float   GetViscoBoundFactor()const{ return(ViscoBoundFactor); };

  bool GetSymmetry()const{ return(Symmetry); };

  bool GetSplitting()const{ return(Splitting); };

  double GetDp()       const{ return(Dp);        };
  double GetH()        const{ return(H);         };
  double GetB()        const{ return(B);         };
  double GetGamma()    const{ return(Gamma);     };
  double GetRhopZero() const{ return(RhopZero);  };
  double GetMassBound()const{ return(MassBound); };
  double GetMassFluid()const{ return(MassFluid); };
  tfloat3 GetGravity() const{ return(Gravity);   };

};



#endif


