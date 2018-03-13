//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# Descripcion:
//:# =============
//:# Clase para cargar ficheros en formato BI4.
//:#
//:# Cambios:
//:# =========
//:# - Carga datos iniciales de particulas. (27-07-2012)
//:# - Carga datos iniciales de particulas. (27-07-2012)
//:# - Carga particulas excluidas de BI2 y BI3. (05-03-2013)
//:# - Arranque de simulaciones iniciadas a partir de BI2 y BI3. (05-03-2013)
//:# - Las funcion GetAllocMemory() devuelve long long. (05-04-2013)
//:# - Remplaza long long por llong. (01-10-2015)
//:# - En el constructor se puede indicar si se usa OMP. (07-07-2016)
//:# - Documentacion del codigo en ingles. (08-08-2017)
//:#############################################################################

/// \file JPartsLoad4.h \brief Declares the class \ref JPartsLoad4.

#ifndef _JPartsLoad4_
#define _JPartsLoad4_

#include "TypesDef.h"
#include "JObject.h"
#include <cstring>

//##############################################################################
//# JPartsLoad4
//##############################################################################
/// \brief Manages the initial load of particle data.

class JPartsLoad4 : protected JObject
{
public:
  typedef enum{ PERI_None=0,PERI_X=1,PERI_Y=2,PERI_Z=4,PERI_XY=3,PERI_XZ=5,PERI_YZ=6,PERI_Unknown=96 }TpPeri; 

protected:
  const bool UseOmp;

  unsigned Npiece;
  bool Simulate2D;         ///<Indicates 2D simulation.
  double Simulate2DPosY;   ///<Y value in 2D simulations.
  bool NpDynamic;          ///<CaseNp can increase.

  ullong CaseNp;           ///<Number of total particles.  
  ullong CaseNfixed;       ///<Number of fixed boundary particles. 
  ullong CaseNmoving;      ///<Number of moving boundary particles. 
  ullong CaseNfloat;       ///<Number of floating boundary particles. 
  ullong CaseNfluid;       ///<Number of fluid particles (including the excluded ones). 

  TpPeri PeriMode;
  tdouble3 PeriXinc;
  tdouble3 PeriYinc;
  tdouble3 PeriZinc;

  bool MapSize;                  ///<Indicates whether MapPosMin and MapPosMax are valid. | Indica si MapPosMin y MapPosMax son validos.
  tdouble3 MapPosMin,MapPosMax;  ///<Domain limits that already include the border. | Limites del dominio que ya incluyen el borde.

  tdouble3 CasePosMin,CasePosMax;

  unsigned PartBegin;
  double PartBeginTimeStep;
  ullong PartBeginTotalNp;        ///<Total number of simulated particles.

  //-Variables for particles.
  unsigned Count;    //-Number of particles.
  unsigned *Idp;
  tdouble3 *Pos;
  tfloat4 *VelRhop;

  void AllocMemory(unsigned count);
  template<typename T> T* SortParticles(const unsigned *vsort,unsigned count,T *v)const;
  void SortParticles();
  void CalculateCasePos();

public:
  JPartsLoad4(bool useomp);
  ~JPartsLoad4();
  void Reset();

  void LoadParticles(const std::string &casedir,const std::string &casename,unsigned partbegin,const std::string &casedirbegin);
  void CheckConfig(ullong casenp,ullong casenfixed,ullong casenmoving,ullong casenfloat,ullong casenfluid,bool perix,bool periy,bool periz)const;
  void CheckConfig(ullong casenp,ullong casenfixed,ullong casenmoving,ullong casenfloat,ullong casenfluid)const;
  void RemoveBoundary();

  unsigned GetCount()const{ return(Count); }

  bool MapSizeLoaded()const{ return(MapSize); }
  void GetMapSize(tdouble3 &mapmin,tdouble3 &mapmax)const;
  void CalculeLimits(double border,double borderperi,bool perix,bool periy,bool periz,tdouble3 &mapmin,tdouble3 &mapmax);

  bool GetSimulate2D()const{ return(Simulate2D); }
  double GetSimulate2DPosY()const{ return(Simulate2DPosY); }
  double GetPartBeginTimeStep()const{ return(PartBeginTimeStep); }
  ullong GetPartBeginTotalNp()const{ return(PartBeginTotalNp); }

  const unsigned* GetIdp(){ return(Idp); }
  const tdouble3* GetPos(){ return(Pos); }
  const tfloat4* GetVelRhop(){ return(VelRhop); }

  tdouble3 GetCasePosMin()const{ return(CasePosMin); }
  tdouble3 GetCasePosMax()const{ return(CasePosMax); }

  llong GetAllocMemory()const;
};

#endif


