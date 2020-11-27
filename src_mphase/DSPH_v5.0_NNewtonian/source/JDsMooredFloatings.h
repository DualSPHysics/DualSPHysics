//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2020 by Dr Jose M. Dominguez
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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Gestion de moorings usando la libreria MoorDyn. (06-02-2017)
//:# - Inclusion en codigo DualSPHysics v4.2.073. (13-08-2018)
//:# - Actualiza coupling con MoorDyn. (23-12-2019)
//:# - Saves VTK files in MooringsVtk directory. (24-12-2019)
//:# - Muestra warning en el caso de que la gravedad no sea definida solo en Z
//:#   y error en caso de ser cero. (12-03-2020)
//:# - Objeto JXml pasado como const para operaciones de lectura. (19-03-2020)  
//:# - Comprueba opcion active en elementos de primer y segundo nivel. (19-03-2020)  
//:# - Cambio de nombre de J.MooredFloatings a J.DsMooredFloatings. (28-06-2020)
//:#############################################################################

/// \file JDsMooredFloatings.h \brief Declares the class \ref JDsMooredFloatings.

#ifndef _JDsMooredFloatings_
#define _JDsMooredFloatings_

#include "JObject.h"
#include "DualSphDef.h"
#include <string>
#include <vector>
#include <cmath>

class JLog2;
class JXml;
class TiXmlElement;
class JDsFtForcePoints;

//##############################################################################
//# XML format in _FmtXML_Moorings.xml.
//##############################################################################

//##############################################################################
//# JDsMooredFloating
//##############################################################################
/// \brief Manages a moored floating body.

class JDsMooredFloating : protected JObject
{
public:
  typedef struct StrLinkData{
    unsigned ftid;     ///<Number of floating body.
    unsigned fairnum;  ///<Number of fairlead.
    tdouble3 linkpos;  ///<Link point in the floating (initial position).
    word ptid;         ///<Ptid in JDsFtForcePoints object.
    StrLinkData(unsigned xftid,unsigned xfairnum,const tdouble3 &xlinkpos,word xptid){ 
      ftid=xftid; fairnum=xfairnum; linkpos=xlinkpos; ptid=xptid;
    }
  }StLinkData;

private:
  JLog2 *Log;
  unsigned FtIdx;  ///<Number of floating with moorings.
  unsigned FtId;   ///<Number of floating in the general list of floatings.

  std::vector<StLinkData> Fairleads;

public:
  const word FloatingMk;   ///<Mkbound of the Floating body the mooring is linked to.

  JDsMooredFloating(word fmk);
  ~JDsMooredFloating();
  void Reset();

  void ConfigIds(unsigned ftidx,unsigned ftid);
  void AddFairlead(unsigned fairnum,const tdouble3 &linkpos,word ptid);
  
  void VisuConfig()const;

  unsigned GetFtIdx()const{ return(FtIdx); }
  unsigned GetFtId()const{ return(FtId); }
  unsigned Count()const{ return(unsigned(Fairleads.size())); }
  StLinkData GetFairlead(unsigned fairnum)const;

};

//##############################################################################
//# JDsMooredFloating
//##############################################################################
/// \brief Manages set of moored floating bodies.

class JDsMooredFloatings : protected JObject
{
private:
  JLog2 *Log;
  const std::string DirCase;
  const std::string CaseName;
  const tfloat3 Gravity;
  std::string FileLines;
  std::string MoordynDir;   ///<Work directory for MoorDyn.

  bool SvVtkMoorings;  ///<Saves vtk with moorings (def=true).
  bool SvCsvPoints;    ///<Saves csv with link points (def=true). 
  bool SvVtkPoints;    ///<Saves vtk with link points (def=false).

  std::vector<JDsMooredFloating*> Floatings; ///<List of floatings with moorings.

  bool MoorDynReady;   ///<Indicate if MoorDyn was initializated.

  bool FairArrays;          ///<Indicate if the fairlead arrays was initializated.
  unsigned FairNftm;        ///<Number of moored floatings.
  unsigned* FairFtmNum;     ///<Number of fairleads per each moored floating [FairNftm].
  double*** FairleadPos;    ///<Stores link positions  [FairNftm][FairFtmNum[cf]][3].  
  double*** FairleadVel;    ///<Stores link velocities [FairNftm][FairFtmNum[cf]][3].  
  double*** FairleadForce;  ///<Stores link forces     [FairNftm][FairFtmNum[cf]][3].  

  //double SaveDataTime;   ///<Saves CSV with data exchange (0=all steps, <0:none).
  //double NextTime;
  //double LastTimeOk;

  unsigned GetFloatingByMk(word mkbound)const;
  void ReadXml(const JXml *sxml,TiXmlElement* ele);
  void ConfigFloatings(unsigned ftcount,const StFloatingData *ftdata);
  void AllocFairMemory();
  void FreeFairMemory();

public:
  JDsMooredFloatings(std::string dircase,std::string casename,tfloat3 gravity);
  ~JDsMooredFloatings();
  void Reset();
  void LoadXml(const JXml *sxml,const std::string &place);
  void Config(unsigned ftcount,const StFloatingData *ftdata,JDsFtForcePoints *forcepoints);

  void VisuConfig(std::string txhead,std::string txfoot)const;

  void ComputeForces(unsigned nstep,double timestep,double dt,JDsFtForcePoints *forcepoints);

  void SaveVtkMoorings(unsigned numfile)const;
  void SaveData(unsigned numfile)const{ if(SvVtkMoorings)SaveVtkMoorings(numfile); }

  unsigned Count()const{ return(unsigned(Floatings.size())); }
};


#endif


