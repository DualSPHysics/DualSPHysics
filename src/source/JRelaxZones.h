//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2020 by Dr. Jose M. Dominguez
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
//# - Clase para la creacion de oleaje actuando sobre la velocidad de las 
//#   particulas. (28-10-2016)
//# - Comprueba opcion active en elementos de primer y segundo nivel. (19-03-2020)  
//#############################################################################

#ifndef _JRelaxZones_
#define _JRelaxZones_

#include "JObject.h"
#include "DualSphDef.h"
#include <string>
#include <vector>

//#define DISABLE_RZ     ///<It allows compile without RZ library.

class JXml;
class TiXmlElement;
class JLog2;
class JRelaxZone;
class JRelaxZoneExternal;
class JRelaxZoneUniform;

//##############################################################################
//# XML format in _FmtXML_SourceGeneration.xml.
//##############################################################################

//##############################################################################
//# JRelaxZones
//##############################################################################

#ifdef DISABLE_RZ
class JRelaxZones : protected JObject
{
public:
  JRelaxZones(bool useomp,bool usegpu,JLog2* log,std::string dirdata
    ,bool withfloatings,unsigned fluidbeginidp,tdouble3 gravity3){}
  ~JRelaxZones(){}
  void Reset(){}
  static bool Available(){ return(false); }

  void LoadFileXml(const std::string &filexml,const std::string &place){}
  void LoadXml(const JXml *sxml,const std::string &place){}

  void Init(std::string dircase,double timemax,double dp){}

  void VisuConfig(std::string txhead,std::string txfoot){}

  unsigned GetCount()const{ return(0); }
  unsigned GetCountExternal()const{ return(0); }
  unsigned GetCountUniform()const{ return(0); }

  void SetFluidVel(double timestep,double dt,unsigned n,unsigned pini
    ,const tdouble3 *pos,const unsigned *idp,tfloat4 *velrhop)const{}
  void SetFluidVelGpu(double timestep,double dt,unsigned n,unsigned pini
    ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop)const{}
};


#else
class JRelaxZones : protected JObject
{
private:
  const bool UseOmp;
  const bool UseGpu;
  JLog2 *Log;
  const std::string DirData;
  const bool WithFloatings;
  const unsigned FluidBeginIdp; ///<Idp for first fluid particle.
  const double Gravity;    ///<Gravity value (always positive).

  std::vector<JRelaxZone*> List;
  std::vector<JRelaxZoneExternal*> ListExternal; //-For external velocity (SWASH).
  std::vector<JRelaxZoneUniform*> ListUniform;   //-For uniform velocity.
  
  void ReadXml_RelaxZone(bool wavereg,const JXml *sxml,TiXmlElement* ele);
  void ReadXml_RelaxZoneExternal(bool onedim,const JXml *sxml,TiXmlElement* ele);
  void ReadXml_RelaxZoneUniform(const JXml *sxml,TiXmlElement* ele);
  void ReadXml(const JXml *sxml,TiXmlElement* lis);
  void SaveDomainVtk(std::string fname,tdouble3 center,float width,double swl,double depth);

public:
  JRelaxZones(bool useomp,bool usegpu,JLog2 *log,std::string dirdata
    ,bool withfloatings,unsigned fluidbeginidp,tdouble3 gravity3);
  ~JRelaxZones();
  void Reset();
  static bool Available(){ return(true); }

  void LoadFileXml(const std::string &filexml,const std::string &place);
  void LoadXml(const JXml *sxml,const std::string &place);

  void Init(std::string dircase,double timemax,double dp);

  void VisuConfig(std::string txhead,std::string txfoot);

  unsigned GetCount()const{ return(unsigned(List.size())); }
  unsigned GetCountExternal()const{ return(unsigned(ListExternal.size())); }
  unsigned GetCountUniform()const{ return(unsigned(ListUniform.size())); }

  void SetFluidVel(double timestep,double dt,unsigned n,unsigned pini
    ,const tdouble3 *pos,const unsigned *idp,tfloat4 *velrhop)const;
  void SetFluidVelGpu(double timestep,double dt,unsigned n,unsigned pini
    ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop)const;
};
#endif

#endif

