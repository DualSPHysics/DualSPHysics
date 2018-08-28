//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2018 by Dr Jose M. Dominguez
 All rights reserved.

 DualSPHysics is an international collaboration between:
 - EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 - School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
 the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the distribution.
 * Neither the name of the DualSPHysics nor the names of its contributors may be used to endorse or promote products derived 
   from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
 SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

//#############################################################################
//# Cambios:
//# =========
//# - Clase para gestionar uno o varios pistones de velocidad-X. (10-05-2014)
//# - Remplaza long long por llong. (01-10-2015)
//#############################################################################

#ifndef _JMLPistons_
#define _JMLPistons_

#include "TypesDef.h"
#include "JObject.h"
#include "Types.h"
#include <string>
#include <vector>

class JXml;
class TiXmlElement;
class JLog2;
class JMLPiston1D;
class JMLPiston2D;
class JMLPistonsGpu;

//##############################################################################
//# XML format in _FmtXML_MLPistons.xml.
//##############################################################################

//##############################################################################
//# JMLPistons
//##############################################################################
class JMLPistons : protected JObject
{
public:
/// Structure with motion information of JMLPiston2D.
  typedef struct {
    unsigned np;
    unsigned idbegin;
    double posymin;
    double poszmin;
    unsigned poszcount;
    const double *movyz;
    const double *velyz;
  }StMotionInfoPiston2D;

private:
  JLog2 *Log;
  const bool UseGpu;
  const bool Cpu;
  const std::string DirData;
  double TimeMod;       ///<Modifies the timestep for piston motion.

  std::vector<JMLPiston1D*> List1d;
  std::vector<JMLPiston2D*> List2d;

  bool PistonsReady;   //-Indica que el objeto ya esta preparado para el calculo de movimiento.

  //-Vars para gestion de pistones 1D de forma conjunta.
  unsigned Piston1dCount;
  double PoszMin;      //-Posicion Z minima de todos los pistones.
  unsigned PoszCount;  //-Numero de posiciones Z de los pistones (comun para todos);

  static const unsigned PistonIdSize=2048;
  byte *PistonId;      //-Contiene el id del piston en funcion del codigo de particula moving [PistonIdSize].
  double *MovVel;      //-Almacena Movx y Velx de todos los pistones 1D (only GPU) [PoszCount*Piston1dCount*2].

  JMLPistonsGpu* PistonsGpu; //-Object for GPU use.

  bool Use_AwasZsurf;  //-AWAS-Zsurf.

  void ReadXml(JXml *sxml,TiXmlElement* lis);
  void PreparePiston1d(double dp,unsigned np,const unsigned *idp,const tdouble3 *pos);

public:
  JMLPistons(bool usegpu,JLog2* log,std::string dirdata);
  ~JMLPistons();
  void Reset();
  llong GetAllocMemoryCpu()const;
  llong GetAllocMemoryGpu()const;

  void LoadFileXml(const std::string &filexml,const std::string &place);
  void LoadXml(JXml *sxml,const std::string &place);
  bool ConfigPiston(word mkbound,word pistonid,unsigned idbegin,unsigned np,double timemax);
  bool ExistsPistonByMk(word mkbound);
  void CheckPistons();
  void PreparePiston(double dp,unsigned np,const unsigned *idp,const tdouble3 *pos);

  void VisuConfig(std::string txhead,std::string txfoot);

  void SetTimeMod(double timemod);
  double GetTimeMod()const{ return(TimeMod); }

  void CalculateMotion1d(double time);

  unsigned GetPiston1dCount()const{ return(unsigned(List1d.size())); }
  JMLPiston1D* GetPiston1d(unsigned num)const;

  double GetPoszMin()const{ return(PoszMin); }
  unsigned GetPoszCount()const{ return(PoszCount); }
  const byte* GetPistonId()const{ return(PistonId); }
  const double* GetMovx()const{ return(MovVel); }
  const double* GetVelx()const{ return(MovVel? MovVel+(PoszCount*Piston1dCount): NULL); }

  const byte* GetPistonIdGpu()const;
  const double* GetMovxGpu()const;
  const double* GetVelxGpu()const;


  unsigned GetPiston2dCount()const{ return(unsigned(List2d.size())); }
  JMLPiston2D* GetPiston2d(unsigned num)const;

  StMotionInfoPiston2D CalculateMotion2d(unsigned num,double time);

//  bool UseAwasZsurf()const{ return(Use_AwasZsurf); }  ->awaspdte
//  void InitAwas(tfloat3 gravity,bool simulate2d,TpCellOrder cellorder,float massfluid,double dp,float dosh,float scell,int hdiv,tdouble3 domposmin,tdouble3 domrealposmin,tdouble3 domrealposmax);
//  void RunAwasCpu(double timestep,bool svdata,tuint3 ncells,tuint3 cellmin
//    ,const unsigned *begincell,const tdouble3 *pos,const typecode *code,const tfloat4 *velrhop);
//#ifdef _WITHGPU
//  void RunAwasGpu(double timestep,bool svdata,tuint3 ncells,tuint3 cellmin
//    ,const int2 *beginendcell,const double2 *posxy,const double *posz,const typecode *code,const float4 *velrhop,float3 *aux);
//#endif
};


#endif


