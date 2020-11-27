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
//:# - GetAccValues() guarda entrada y salida para evitar calculos con llamadas 
//:#   consecutivas iguales. (13-09-2019)
//:# - Objeto JXml pasado como const para operaciones de lectura. (18-03-2020)  
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (18-03-2020)  
//:# - Permite definir un intervalo de tiempo para su activacion. (26-04-2020)  
//:# - Cambio de nombre de J.SphAccInput a J.DsAccInput. (28-06-2020)
//:#############################################################################

/// \file JDsAccInput.h \brief Declares the class \ref JDsAccInput.

#ifndef _JDsAccInput_
#define _JDsAccInput_

#include "JObject.h"
#include "DualSphDef.h"
#ifdef _WITHGPU
  #include <cuda_runtime_api.h>
#endif

#include <string>
#include <vector>

class JLog2;
class JXml;
class TiXmlElement;
class JLinearValue;
class JSphMk;

//##############################################################################
//# XML format in _FmtXML_AccInput.xml.
//##############################################################################

//##############################################################################
//# JDsAccInputMk
//##############################################################################
/// \brief Provides the force to be applied to different blocks of particles that is loaded from files.

class JDsAccInputMk : protected JObject
{
public:
  const unsigned Idx;     ///<Index of configuration.
  const bool Bound;       ///<Type of target particles (boundary floating or fluid).
  const word MkType1;     ///<The MK bound or fluid to select the target particles.
  const word MkType2;     ///<Final range of MK values to select the target particles.
  const double TimeIni;   ///<Initial time of activation (default=0).
  const double TimeEnd;   ///<Final time of activation (default=DBL_MAX).
  const bool GravityEnabled;  ///<Determines whether global gravity is enabled or disabled for this particle set SL
  const tfloat3 AccCoG;       ///<The centre of gravity that will be used for angular acceleration calculations.

protected:
  JLog2* Log;

  typecode CodeSel1; ///<First code of target particles (TypeAndValue).
  typecode CodeSel2; ///<Last code of target particles (TypeAndValue).

  JLinearValue *AceData; ///<Input acceleration data.
  JLinearValue *VelData; ///<Input velocity data.

  double LastTimestepInput; ///<Saves the last value used with GetAccValues().
  StAceInput LastOutput;    ///<Saves the last value returned by GetAccValues().

  void Reset();

public:
  JDsAccInputMk(unsigned idx,bool bound,word mktype1,word mktype2
    ,double tini,double tend,bool genabled,tfloat3 acccentre
    ,const JLinearValue &acedata,const JLinearValue &veldata);
  ~JDsAccInputMk();
  long long GetAllocMemory()const;

  void ConfigCodeSel(typecode codesel1,typecode codesel2){ 
    CodeSel1=codesel1; CodeSel2=codesel2; 
  }
  typecode GetCodeSel1()const{ return(CodeSel1); };
  typecode GetCodeSel2()const{ return(CodeSel2); };

  void GetConfig(std::vector<std::string> &lines)const;

  //word GetMkFluid()const{ return(MkFluid); }
  const StAceInput& GetAccValues(double timestep); //SL: Added linear and angular velocity and set gravity flag
};

//##############################################################################
//# JDsAccInput
//##############################################################################
/// \brief Manages the application of external forces to different blocks of particles (with the same MK).

class JDsAccInput : protected JObject
{
protected:
  JLog2* Log;
  std::string DirData;
  std::vector<JDsAccInputMk*> Inputs;
  long long MemSize;

  void Reset();
  bool ExistMk(bool bound,word mktype)const;
  void LoadXml(const JXml *sxml,const std::string &place);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void ComputeVelocity(const JLinearValue &acedata,JLinearValue &veldata)const;

public:
  JDsAccInput(const std::string &dirdata,const JXml *sxml,const std::string &place);
  ~JDsAccInput();
  long long GetAllocMemory()const{ return(MemSize); }

  void VisuConfig(std::string txhead,std::string txfoot)const;
  void Init(const JSphMk *mkinfo);

  unsigned GetCount()const{ return(unsigned(Inputs.size())); };
  const StAceInput& GetAccValues(unsigned cinput,double timestep); //SL: Added linear and angular velocity and set gravity flag

  void RunCpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
    ,const typecode *code,const tdouble3 *pos,const tfloat4 *velrhop,tfloat3 *ace);

#ifdef _WITHGPU
  void RunGpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
    ,const typecode *code,const double2 *posxy,const double *posz,const float4 *velrhop,float3 *ace);
#endif
};

#endif


