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
//:# - Gestiona la aplicacion de fuerzas externas a floatings. (07-02-2017)
//:# - Inclusion en codigo DualSPHysics v4.2.073. (13-08-2018)
//:# - Cambia uso de JFormatFiles2 por JVtkLib. (12-12-2019)
//:# - Actualiza coupling con MoorDyn. (23-12-2019)
//:# - Saves VTK files in MooringsVtk directory. (24-12-2019)
//:# - Cambio de nombre de J.SphFtForcePoints a J.DsFtForcePoints. (28-06-2020)
//#############################################################################

/// \file JDsFtForcePoints.h \brief Declares the class \ref JDsFtForcePoints.

#ifndef _JDsFtForcePoints_
#define _JDsFtForcePoints_

#include "JObject.h"
#include "DualSphDef.h"
#include <string>
#include <vector>
#include <cmath>

class JLog2;
class JSphMk;

//##############################################################################
//# JDsFtForcePoints
//##############################################################################
/// \brief Allows external forces to be imposed on floating objects.

class JDsFtForcePoints : protected JObject
{
private:
  JLog2 *Log;
  const bool Cpu;
  const double Dp;     ///<Initial distance between particles [m].

  bool SvCsvPoints;    ///<Saves csv with link points (def=true). 
  bool SvVtkPoints;    ///<Saves vtk with link points (def=false).

  byte PeriActive;
  bool PeriX;          ///<Periodic conditions in X.
  bool PeriY;          ///<Periodic conditions in Y.
  bool PeriZ;          ///<Periodic conditions in Z.
  tdouble3 PeriXinc;   ///<Value that is added at the outer limit to modify the position.
  tdouble3 PeriYinc;   ///<Value that is added at the outer limit to modify the position.
  tdouble3 PeriZinc;   ///<Value that is added at the outer limit to modify the position.

  const word FtCount; ///<Number of floating objects.
  float *FtRadius;    ///<Maximum distance between particles and center [FtCount].
  float *FtMass;      ///<Mass of the floating object [FtCount].

  word *FtCountPt;    ///<Number of points for each floating [FtCount].
  word *FtBeginPt;    ///<Begin of points for each floating [FtCount].

  //-Variables to apply forces on floating bodies.
  word PtCount;       ///<Total number of points.
  word *PtId;         ///<Id for points (for initialization) [PtCount].
  word *PtFtid;       ///<Id of floating (for initialization) [PtCount].
  float *PartDist;    ///<Initial distance to nearest particle [PtCount].

  double TimeStep;    ///<Current instant of the simulation.
  tdouble3 *PtPos;    ///<Application position [PtCount].
  tfloat3 *PtVel;     ///<Velocity [PtCount].
  tfloat3 *PtForce;   ///<Force [PtCount].

  word SelFtCount;      ///<Number of floating objects with force points.  
  word *SelFtIndex;     ///<Index to selected floatings. Uses USHRT_MAX when index is invalid [FtCount].
  word *SelFtid;        ///<Id of floating [SelFtCount].
  tdouble3 *SelFtCenter;///<Center of the selected floating [SelFtCount].
  tfloat3 *SelFtAce;    ///<Total ace of each selected floating [SelFtCount].
  tfloat3 *SelFtOmega;  ///<Total omega of each selected floating [SelFtCount].

  void AllocMemoryFt(word ftcount);
  void AllocMemorySelFt(word selftcount);
  void ResizeMemoryPt(word ptcount);
  inline tfloat3 FtPeriodicDist(const tdouble3 &pos,const tdouble3 &center,float radius)const;
  void ConfigPeri(byte periactive,bool perix,bool periy,bool periz,tdouble3 perixinc,tdouble3 periyinc,tdouble3 perizinc);
  void SaveVtkPoints(unsigned numfile)const;
  void SaveCsvPoints(unsigned numfile)const;

public:
  JDsFtForcePoints(bool iscpu,double dp,unsigned ftcount);
  ~JDsFtForcePoints();
  void Reset();
  llong GetAllocMemory()const;

  void SetSaveData(bool savecsv,bool savevtk){ SvCsvPoints|=savecsv; SvVtkPoints|=savevtk; }
  word AddPoint(unsigned ftid,const tdouble3 &pos);
  word GetIdx(word ptid)const;

  void Config(unsigned ftcount,const StFloatingData *ftdata,byte periactive,bool perix,bool periy,bool periz,tdouble3 perixinc,tdouble3 periyinc,tdouble3 perizinc);
  void CheckPoints(const JSphMk *mkinfo,unsigned np,const unsigned *idp,const tdouble3 *pos);
  void VisuConfig(std::string txhead,std::string txfoot,unsigned ftcount,const StFloatingData *ftdata)const;

  void UpdatePoints(double timestep,double dt,const StFloatingData *ftdata);

  tdouble3 GetPos(word idx)const{ return(PtPos[idx]); }
  tfloat3 GetVel(word idx)const{ return(PtVel[idx]); }

  void SetForce(word idx,const tfloat3 &force){ PtForce[idx]=force; }
  void ComputeFtMotion();
  void GetFtMotionData(StFtoForces *ftoforces)const;

  void SaveData(unsigned numfile)const;
};


#endif


