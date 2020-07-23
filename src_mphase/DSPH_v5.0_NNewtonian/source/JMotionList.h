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

/// \file JMotionList.h \brief Declares the classes \ref JMotionList and \ref JMotionListData.

#ifndef _JMotionList_
#define _JMotionList_

#include "TypesDef.h"
#include "JObject.h"
#include "JMatrix4.h"
#include <string>
#include <cstdlib>
#include <ctime>


//##############################################################################
//# JMotionListData
//##############################################################################
/// \brief Manages the state of different motions.

class JMotionListData : protected JObject
{
private:
  bool Active;
  //byte State;  //-0:Sin movimiento, 1:Movimiento finalizado, 2:Inicia movimiento, 3:Continua movimiento
  bool TypeSimple;
  tdouble3 MvSimple1;
  tmatrix4d MvMatrix1;
  tdouble3 MvSimple2;
  tmatrix4d MvMatrix2;

  tdouble3 VelSimple;
  tdouble3 AceSimple;

public:
  JMotionListData();
  void PreMotion(){  Active=false;  }

  //-Only for simple motion.
  void Sp_Movedt(const tdouble3  &mvsimple,double dt);
  void Sp_Movedt(const tmatrix4d &mvmatrix,double dt);

  //-Only for double dt acceleration motion.
  void Ace2_Move1dt(const tdouble3 &mvsimple);
  void Ace2_Move2dt(const tdouble3 &mvsimple);
  void Ace2_Move1dt(const tmatrix4d &mvmatrix);
  void Ace2_Move2dt(const tmatrix4d &mvmatrix);
  void Ace2_PosMotion(double dt);

  //bool GetActive()    const{ return(Active);     }
  //bool GetTypeSimple()const{ return(TypeSimple); }
  bool GetData(bool &typesimple,tdouble3 &simplemov,tdouble3 &simplevel
    ,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2)const;
  bool GetData(bool &typesimple,tdouble3 &simplemov,tmatrix4d &matmov)const;
};

//##############################################################################
//# JMotionList
//##############################################################################
/// \brief Manages the list of different motions.

class JMotionList : protected JObject
{
private:
  JMotionListData* MotionData;

public:
  const unsigned Nref;
  double TimeStep;

  JMotionList(int nref);
  ~JMotionList();

  void PreMotion(){ for(unsigned c=0;c<Nref;c++)MotionData[c].PreMotion(); }

  //-Only for simple motion.
  void Sp_Movedt(unsigned ref,const tdouble3  &mvsimple,double dt){  MotionData[ref].Sp_Movedt(mvsimple,dt);  }
  void Sp_Movedt(unsigned ref,const tmatrix4d &mvmatrix,double dt){  MotionData[ref].Sp_Movedt(mvmatrix,dt);  }

  //-Only for double dt acceleration motion.
  void Ace2_Move1dt(unsigned ref,const tdouble3  &mvsimple){  MotionData[ref].Ace2_Move1dt(mvsimple);  }
  void Ace2_Move2dt(unsigned ref,const tdouble3  &mvsimple){  MotionData[ref].Ace2_Move2dt(mvsimple);  }
  void Ace2_Move1dt(unsigned ref,const tmatrix4d &mvmatrix){  MotionData[ref].Ace2_Move1dt(mvmatrix);  }
  void Ace2_Move2dt(unsigned ref,const tmatrix4d &mvmatrix){  MotionData[ref].Ace2_Move2dt(mvmatrix);  }
  void Ace2_PosMotion(double dt){ for(unsigned c=0;c<Nref;c++)MotionData[c].Ace2_PosMotion(dt); }

  bool GetData(unsigned ref,bool &typesimple,tdouble3 &simplemov,tdouble3 &simplevel,tdouble3 &simpleace
    ,tmatrix4d &matmov,tmatrix4d &matmov2)const;
  bool GetData(unsigned ref,bool &typesimple,tdouble3 &simplemov,tmatrix4d &matmov)const;
};


#endif

