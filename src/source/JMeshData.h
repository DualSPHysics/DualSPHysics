//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - Almacena datos de un malla para un instante y ofrece metodos para su
//:#   gestion. (04-08-2020)
//:#############################################################################

/// \file JMeshData.h \brief Declares the class \ref JMeshData.

#ifndef _JMeshData_
#define _JMeshData_

#include <string>
#include "JObject.h"
#include "JMeshDataDef.h"

class JDataArrays;
namespace jcsv{
  class JSaveCsv2;
}

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshData
//##############################################################################
/// \brief Manages measurements on a grid of positions of one instant.

class JMeshData : protected JObject
{
public:
  static const int TAGDATA12=12;
private:
  bool UseMesh;     ///<Define use of mesh instead of list of points.
  //-Mesh definition (UseMesh==true).
  tdouble3 PtRef;   ///<Initial mesh position.
  tdouble3 Vdp1;    ///<First axis vector with distance between points.
  tdouble3 Vdp2;    ///<Second axis vector with distance between points.
  tdouble3 Vdp3;    ///<Third axis vector with distance between points (used for swl calculation).
  //-Points definition (UseMesh==false).
  tdouble3* PtPos;  ///<Positions [Npt12*Npt3 = Npt].
  //-Number of points per axis.
  unsigned Npt1;    ///<Number of points for first axis.
  unsigned Npt2;    ///<Number of points for second axis.
  unsigned Npt3;    ///<Number of points for third axis.
  unsigned Npt;     ///<Total number of points (Npt=Npt1*Npt2*Npt3).
  unsigned Npt12;   ///<Total number of points (Npt12=Npt1*Npt2).

  tfloat3 DirDat;  ///<Direction vector for computed linear velocity or other variables.

  JDataArrays* Arrays;

  double TimeStep;    ///<Time of of data [s].
private:
  static void RunExceptioonStatic(const std::string& srcfile,int srcline
    ,const std::string& method
    ,const std::string& msg,const std::string& file="");

public:
  //void DbgRotateXtoY(){
  //  unsigned aux=Npt1; Npt1=Npt2; Npt2=aux;
  //  tdouble3 vux=Vdp1; Vdp1=Vdp2; Vdp2=vux;
  //  PrintMeshPts("Changed::",GetMeshPt());
  //}

public:
  JMeshData();
  ~JMeshData();
  void Reset(bool createarrays);
  void CopyDataFrom(const JMeshData* obj);
  void ConfigMesh(const jmsh::StMeshPts& mesh,double timestep=0
    ,std::string varlist="");
  void ConfigPtos(unsigned npt,const tdouble3* pos
    ,tfloat3 dirdat,double timestep=0,std::string varlist="");
  void AddVarList(std::string varlist);

  std::string GetVarList(bool novar12,bool var12)const;

  bool EqualStructure(const JMeshData& mdat)const;
  bool EqualStructure(const JMeshData* mdat)const;

  void ClearData();

  void ReverseData(std::string varlist);
  void SetMulData(std::string varlist,double v2);
  void SetAddData(std::string varlist,double v2);

  tdouble3 GetPtRef()const{ return(PtRef); }
  void SetPtRef(const tdouble3& ps){ PtRef=ps; }

  double GetTimeStep()const{ return(TimeStep); }
  void SetTimeStep(double timestep){ TimeStep=timestep; }
  StMeshPts GetMeshPt()const;
  bool GetUseMesh()const{ return(UseMesh); }

  unsigned GetNpt1() const{ return(Npt1); }
  unsigned GetNpt2() const{ return(Npt2); }
  unsigned GetNpt3() const{ return(Npt3); }
  unsigned GetNpt12()const{ return(Npt12); }
  unsigned GetNpt()  const{ return(Npt); }

  tdouble3* GetPtPos()const{ return(PtPos); }

  float*    CreateArrayPtrFloat  (std::string fullname,int tag,unsigned size,float*    data=NULL);
  tfloat3*  CreateArrayPtrFloat3 (std::string fullname,int tag,unsigned size,tfloat3*  data=NULL);
  double*   CreateArrayPtrDouble (std::string fullname,int tag,unsigned size,double*   data=NULL);
  tdouble3* CreateArrayPtrDouble3(std::string fullname,int tag,unsigned size,tdouble3* data=NULL);

  const JDataArrays* GetArrays()const{ return(Arrays); } 
  unsigned GetVarCount()const;
  std::string GetVarName(unsigned idx)const;
  int GetVarTag(unsigned idx)const;

  tfloat3* GetVarFloat3(const std::string& name);
  float*   GetVarFloat (const std::string& name);

  void GetPosf(unsigned np,tfloat3*  pos)const;
  void GetPosd(unsigned np,tdouble3* pos)const;

  static void GetMeshPos(const StMeshPts& m,unsigned np,tfloat3*  pos);
  static void GetMeshPos(const StMeshPts& m,unsigned np,tdouble3* pos);

  static StMeshPts MakeMeshPt(const StMeshBasic& b);
  static void PrintMeshPts(std::string text,const StMeshPts& m);
};

}

#endif


