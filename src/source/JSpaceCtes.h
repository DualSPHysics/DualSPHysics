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
//:# - Se paso a usar double en lugar de float. (25-11-2013)
//:# - El valor de Eps pasa a ser opcional para mantener compatibilidad. (08-01-2015)
//:# - Se cambio Coefficient por CoefH pero manteniendo compatibilidad. (08-01-2015)
//:# - Se anhadio SpeedSound para asignar valor de forma explicita. (08-01-2015)
//:# - Se anhadieron comentarios al escribir el XML. (08-01-2015)
//:# - Se ampliaron los limites de CFLnumber de (0.1-0.5) a (0.001-1). (08-01-2015)
//:# - <speedsystem> y <speedsound> pasan a ser opcionales. (20-01-2015)
//:# - <eps> solo se pasa a <constants> cuando este definido en <constantsdef>. (20-01-2015)
//:# - Se muestran unidades de las constantes. (15-12-2015)
//:# - Nueva funcion estatica para calcular constantes. (19-01-2016)
//:# - Ahora se puede definir <coefh> o <hdp> pero no ambos a la vez. (29-01-2016)
//:# - Se cambio las unidades de la constante B a kg/(m*s^2). (08-06-2016)
//:# - Se cambio las unidades de la constante B a Pascal (Pa). (07-06-2017)
//:# - Save example of <hdp> in <constantsdef> of template. (07-04-2019)  
//:# - <lattice> pasa a ser opcional con el valor 1 por defecto. (03-03-2020)  
//:# - Mejora uso de excepciones. (03-03-2020)  
//:# - Objeto JXml pasado como const para operaciones de lectura. (17-03-2020)  
//:# - Mejora la gestion de excepciones. (06-05-2020)
//:#############################################################################

/// \file JSpaceCtes.h \brief Declares the class \ref JSpaceCtes.

#ifndef _JSpaceCtes_
#define _JSpaceCtes_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"

class JXml;
class TiXmlElement;

//##############################################################################
//# JSpaceCtes
//##############################################################################
/// \brief Manages the info of constants from the input XML file.

class JSpaceCtes : protected JObject 
{
public:
  ///// Defines structure to calculate constants. 
  //typedef struct StrConstants{
  //  bool data2d;
  //  double data2dposy;
  //  tdouble3 gravity;
  //  double dp,coefh,coefhdp;
  //  double hswl,speedsystem,coefsound,speedsound;
  //  double gamma,rhop0;
  //  double cteh,cteb;
  //  double massbound;
  //  double massfluid;

  //  StrConstants(){ Clear(); }
  //  StrConstants(bool vdata2d,double vdata2dposy,tdouble3 vgravity,double vdp,double vcoefh,double vcoefhdp,double vhswl
  //    ,double vspeedsystem,double vcoefsound,double vspeedsound,double vgamma,double vrhop0
  //    ,double vcteh,double vcteb,double vmassbound,double vmassfluid)
  //  {
  //    data2d=vdata2d; data2dposy=vdata2dposy; gravity=vgravity; dp=vdp; coefh=vcoefh; coefhdp=vcoefhdp;
  //    hswl=vhswl; speedsystem=vspeedsystem; coefsound=vcoefsound; speedsound=vspeedsound; gamma=vgamma; 
  //    rhop0=vrhop0; cteh=vcteh; cteb=vcteb; massbound=vmassbound; massfluid=vmassfluid;
  //  }
  //  void Clear(){ 
  //    data2d=false; data2dposy=0; gravity=TDouble3(0);
  //    dp=hswl=speedsystem=coefsound=speedsound=coefh=coefhdp=gamma=rhop0=cteh=cteb=massbound=massfluid=0;
  //  }
  //}StConstants;

private:
  bool Data2DDefined;     ///<Toggles 2D simulation (cancels forces in Y axis).
  bool Data2D;            ///<Data dimension (2D, 3D)) 2D simulation (cancels forces in Y axis).
  double Data2DPosY;      ///<Y value in 2D simulations.    
  int LatticeBound;       ///<Lattice to create boundary particles on its nodes.
  int LatticeFluid;       ///<Lattice to create fluid particles on its nodes.
  tdouble3 Gravity;       ///<Gravity acceleration.
  double CFLnumber;       ///<CFL number (0.001-1).
  bool HSwlAuto;          ///<Activates the automatic computation of H_Swl.
  double HSwl;            ///<Maximum height of the volume of fluid.
  bool SpeedSystemAuto;   ///<Activates the automatic computation of SpeedSystem.
  double SpeedSystem;     ///<Maximum system speed.
  double CoefSound;       ///<Coefficient to multiply speedsystem.
  bool SpeedSoundAuto;    ///<Activates the automatic computation of SpeedSound.
  double SpeedSound;      ///<Speed of sound to use in the simulation (by default speedofsound=coefsound*speedsystem).

  double CoefH;           ///<Coefficient to calculate the smoothing length H (H=coefficient*sqrt(3*dp^2) in 3D).
  double CoefHdp;         ///<Relationship between h and dp. (it is optional).
  double Gamma;           ///<Polytropic constant. (1-7).
  double Rhop0;           ///<Density of reference.

  double Eps;             ///<Epsilon constant for XSPH variant.
  bool EpsDefined;        ///<Epsilon was defined in constantsdef.

  bool HAuto;             ///<Activates the automatic computation of H.
  bool BAuto;             ///<Activates the automatic computation of B.
  bool MassBoundAuto;     ///<Activates the automatic computation of MassBound.
  bool MassFluidAuto;     ///<Activates the automatic computation of MassFluid.
  double H;               ///<Smoothing length.
  double B;               ///<Constant that sets a limit for the maximum change in density.
  double MassBound;       ///<Mass of a boundary particle.
  double MassFluid;       ///<Mass of a fluid particle.

  //-Computed values:
  double Dp;              ///<Inter-particle distance.

  void ReadXmlElementAuto(JXml *sxml,TiXmlElement* node,bool optional,std::string name,double &value,bool &valueauto);
  void WriteXmlElementAuto(JXml *sxml,TiXmlElement* node,std::string name,double value,bool valueauto,std::string comment="",std::string unitscomment="")const;

  void WriteXmlElementComment(TiXmlElement* ele,std::string comment="",std::string unitscomment="")const;

  void ReadXmlDef(JXml *sxml,TiXmlElement* ele);
  void WriteXmlDef(JXml *sxml,TiXmlElement* ele,bool svtemplate)const;
  void ReadXmlRun(const JXml *sxml,TiXmlElement* ele);
  void WriteXmlRun(JXml *sxml,TiXmlElement* ele)const;

public:
  //static StConstants CalcConstans(StConstants cte);
  JSpaceCtes();
  void Reset();
  void LoadDefault();
  void LoadXmlDef(JXml *sxml,const std::string &place);
  void SaveXmlDef(JXml *sxml,const std::string &place,bool svtemplate=false)const;
  void LoadXmlRun(const JXml *sxml,const std::string &place);
  void SaveXmlRun(JXml *sxml,const std::string &place)const;

  bool GetData2D()const{ return(Data2D); }
  double GetData2DPosY()const{ return(Data2DPosY); }
  int GetLatticeBound()const{ return(LatticeBound); }
  int GetLatticeFluid()const{ return(LatticeFluid); }
  tdouble3 GetGravity()const{ return(Gravity); }
  double GetCFLnumber()const{ return(CFLnumber); }
  bool GetHSwlAuto()const{ return(HSwlAuto); }
  double GetHSwl()const{ return(HSwl); }
  bool GetSpeedSystemAuto()const{ return(SpeedSystemAuto); }
  double GetSpeedSystem()const{ return(SpeedSystem); }
  double GetCoefSound()const{ return(CoefSound); }
  bool GetSpeedSoundAuto()const{ return(SpeedSoundAuto); }
  double GetSpeedSound()const{ return(SpeedSound); }
  double GetCoefH()const{ return(CoefH); }
  double GetCoefHdp()const{ return(CoefHdp); }
  double GetCoefficient()const{ return(GetCoefH()); }
  double GetGamma()const{ return(Gamma); }
  double GetRhop0()const{ return(Rhop0); }
  double GetEps()const{ return(Eps); }

  void SetData2D(bool data2d,double data2dposy=0){ Data2D=data2d; Data2DPosY=(data2d? data2dposy: 0); Data2DDefined=true; }
  void SetLatticeBound(bool simple){ LatticeBound=(simple? 1: 2); }
  void SetLatticeFluid(bool simple){ LatticeFluid=(simple? 1: 2); }
  void SetGravity(const tdouble3& g){ Gravity=g; }
  void SetCFLnumber(double v){ 
    if(!v)Run_Exceptioon("Value cannot be zero.");
    if(v>1)Run_Exceptioon("Value cannot be greater than 1.");
    CFLnumber=v;
  }
  void SetHSwlAuto(bool on){ HSwlAuto=on; }
  void SetHSwl(double v){ HSwl=v; }
  void SetSpeedSystemAuto(bool on){ SpeedSystemAuto=on; }
  void SetSpeedSystem(double v){ SpeedSystem=v; }
  void SetCoefSound(double v){ CoefSound=v; }
  void SetSpeedSoundAuto(bool on){ SpeedSoundAuto=on; }
  void SetSpeedSound(double v){ SpeedSound=v; }
  void SetCoefH(double v){ CoefH=v; CoefHdp=0; }
  void SetCoefHdp(double v){ if(v){ CoefHdp=v; CoefH=0; } }
  void SetCoefficient(double v){ SetCoefH(v); }
  void SetGamma(double v){ Gamma=v; }
  void SetRhop0(double v){ Rhop0=v; }
  void SetEps(double v){ Eps=v; }

  bool GetHAuto()const{ return(HAuto); }
  bool GetBAuto()const{ return(BAuto); }
  bool GetMassBoundAuto()const{ return(MassBoundAuto); }
  bool GetMassFluidAuto()const{ return(MassFluidAuto); }
  double GetH()const{ return(H); }
  double GetB()const{ return(B); }
  double GetMassBound()const{ return(MassBound); }
  double GetMassFluid()const{ return(MassFluid); }

  void SetHAuto(bool on){ HAuto=on; }
  void SetBAuto(bool on){ BAuto=on; }
  void SetMassBoundAuto(bool on){ MassBoundAuto=on; }
  void SetMassFluidAuto(bool on){ MassFluidAuto=on; }
  void SetH(double v){ H=v; }
  void SetB(double v){ B=v; }
  void SetMassBound(double v){ MassBound=v; }
  void SetMassFluid(double v){ MassFluid=v; }

  double GetDp()const{ return(Dp); }
  void SetDp(double v){ Dp=v; }

  double ComputeFinalH(bool data2d,double dp)const;
};

#endif


