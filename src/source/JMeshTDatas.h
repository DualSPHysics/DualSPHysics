//HEAD_DSCODES
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
//:# Descripcion:
//:# =============
//:# Clase para cargar y procesar mesh de datos.
//:# - Implementacion. (13-08-2020)
//:# - Nuevos metodos IntpComputeFloat() y IntpComputeFloat3(). (12-09-2020)
//:# - Soporta definicion de loops desde la carga de datos por fichero. (03-08-2022)
//:#############################################################################

/// \file JMeshTDatas.h \brief Declares the class \ref JMeshTDatas.

#ifndef _JMeshTDatas_
#define _JMeshTDatas_

#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"
#include "JMeshDataDef.h"
#include <string>
#include <vector>
#include <fstream>
#include <cfloat>

class JDataArrays;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

class JMeshData;

//##############################################################################
//# JMeshTDatas
//##############################################################################
/// \brief Allows loading and processing mesh data.

class JMeshTDatas : protected JObject
{
 protected:
  const std::string AppName;  ///<Application Name.
  std::string FileIn;    ///<Name of input file.

  std::vector<JMeshData*> Datas;

  double LoopTmax; ///<Final time of loop (DBL_MAX=disabled).
  double LoopTsub; ///<Time interval to begin time of loop (0=disabled).
  double LoopTbeg; ///<Begin time of loop (LoopTmax-LoopTsub).

  //-Auxiliary array with times for interpolation.
  unsigned SizeTimes;
  double* Times;

  //-Variables with last search.
  double TimeStep;
  unsigned Position;
  unsigned PositionNext;
  double TimePre;
  double TimeNext;
  double TimeFactor;

  //-Variables for interpolation I (from Datas[] computed by IntpPrepare()).
  bool InPrepared;
  StMeshPts InMesh;
  bool Dir1,Dir2,Dir3;
  byte InMode;
  byte InMode12;
  unsigned Npt1,Npt2,Npt3,Npt12,Npt;
  tplane3d Pla1;
  tplane3d Pla2;
  tplane3d Pla3;
  tdouble3 Dp;
  tuint3 Cmax;

  /// Modes of interploation configuration.
  typedef enum{ 
    IntpNone=0  ///<Unconfigured.
   ,IntpMesh=1  ///<Configured for fixed mesh mode.
   ,IntpPos=2   ///<Configured for fixed list of positions mode.
   ,IntpVar=3   ///<Configured for variable list of positions mode.
  }TpIntpMode; 

  //-Variables for interpolation II (from fixed positions computed by IntpAllocMemory()).
  TpIntpMode InConfigured;  ///<Configure mode: IntpNone:no configured, IntpMesh:mesh mode, IntpPos:pos list mode.
  StMeshPts OutMesh;     ///<Mesh definition for interpolation points (InConfigured==IntpMesh).
  tdouble3* OutPos;      ///<Position of interpolation points (InConfigured==IntpPos). [NposTot]
  unsigned Npos12;
  unsigned NposTot;
  tuint2*   PosCell1;  ///<Index in data array for values 0 and 1 in axix 1. [NposTot]
  tuint2*   PosCell2;  ///<Index in data array for values 0 and 1 in axix 2. [NposTot]
  tuint2*   PosCell3;  ///<Index in data array for values 0 and 1 in axix 3. [NposTot]
  tfloat3*  PosFdis;   ///<Factor distance from corner of cell. [NposTot]

  unsigned Ctime1;       ///<Time of position-interpolation data in Arrays1. It is usually Times[Ctime1]<=timestep
  JDataArrays* Arrays1;  ///<Auxiliary memory for position-interpolation. [Npos12 or NposTot]

  unsigned Ctime2;       ///<Time of position-interpolation data in Arrays1. It is usually Ctime1+1
  JDataArrays* Arrays2;  ///<Auxiliary memory for position-interpolation. [Npos12 or NposTot]

  //-Variables for interpolation III (from free positions computed by IntpConfigFreePos()).
  double FrTimeStep;
  double FrTimeFactor;
  unsigned FrCtime;
  bool FrMode12;
  unsigned FrSize;
  TpTypeData FrType;
  void* FrData;

  byte FrMode;
  unsigned FrCmax1,FrCmax2,FrCmax3;
  tplane3d FrPla1,FrPla2,FrPla3;
  unsigned FrNum1,FrNum2;



 public:
  bool PrintTime;          ///<Prints estimated finished time for some time-consuming tasks.

 protected:
  void ClearDatas();

  void ResetTimes();
  void PrepareTimes();
  void IntpReset();

  tplane3d ComputeAxisPlane(const tdouble3& ps,const tdouble3& vdp)const;
  void IntpPrepare();
  void IntpFreeMemory();
  void IntpAllocMemory(TpIntpMode mode,unsigned npos12,unsigned npostot);
  void IntpPrepareArray(unsigned npos12,unsigned npostot,JDataArrays* arrays)const;
  inline float IntpComputePosCellDist(const tdouble3& ps,const tplane3d& pla,unsigned cmax,tuint2& cell)const;

  void FindTime(const double timestep);
  template<class T> void IntpComputePosVar(unsigned ctime,bool data12,const T* data,unsigned npos,T* res);
  void IntpComputePosArray(unsigned ctime,JDataArrays* ars);
  void IntpComputePosIntervals(double t);

  unsigned IntpCheckVarName(const std::string& varname,unsigned size)const;
  unsigned IntpCheckVarNum(unsigned varnum,unsigned size)const;
  template<class T,class T1> void IntpComputeVarT(double timestep,unsigned varnum,T* res);

  //-Protected methods for Free-Interpolation.
  void IntpAllocFrMemory(unsigned size,TpTypeData type);
  template<class T> void IntpComputeTime(double timestep);
  template<class T> inline T IntpComputeVarFr1(const tdouble3& ps)const;
  template<class T> inline T IntpComputeVarFr2(const tdouble3& ps)const;
  template<class T> inline T IntpComputeVarFr3(const tdouble3& ps)const;

  void  IntpComputeTime_f(double timestep);
  float IntpComputeVarFr1_f(const tdouble3& ps);
  float IntpComputeVarFr2_f(const tdouble3& ps);
  float IntpComputeVarFr3_f(const tdouble3& ps);

  void  IntpComputeTime_f3(double timestep);
  tfloat3 IntpComputeVarFr1_f3(const tdouble3& ps);
  tfloat3 IntpComputeVarFr2_f3(const tdouble3& ps);
  tfloat3 IntpComputeVarFr3_f3(const tdouble3& ps);

 public:
  JMeshTDatas(std::string appname="",bool printtime=false);
  ~JMeshTDatas();
  void Reset();

  //llong GetAllocMemory()const;

  void LoadFile(const std::string file,std::string varlist=""
    ,double tmin=DBL_MAX,double tmax=DBL_MAX
    ,double looptmax=DBL_MAX,double looptbegin=DBL_MAX
    ,tdouble3 pmin=TDouble3(DBL_MAX),tdouble3 pmax=TDouble3(DBL_MAX)
    ,double settime=0,tdouble3 setpos=TDouble3(0));

  void AddData(JMeshData* mdat);

  void ReverseData(std::string varlist);
  void SetMulData(std::string varlist,double v2);
  void SetAddData(std::string varlist,double v2);

  unsigned GetCount()const{ return(unsigned(Datas.size())); }
  const JMeshData* GetData(unsigned cdata)const{ return(cdata<GetCount()? Datas[cdata]: NULL); }

  double GetTime(unsigned ct)const;
  unsigned GetTimes(std::vector<double>& vtimes)const;
  StMeshPts GetMeshPt()const;

  std::string GetFileIn()const{ return(FileIn); }

  double GetLoopTmax()const{ return(LoopTmax); };
  double GetLoopTsub()const{ return(LoopTsub); };
  double GetLoopTbeg()const{ return(LoopTbeg); };

  //-Configures interpolation (for fixed positions).
  void IntpConfigFixedMesh(const StMeshPts& m);
  void IntpConfigFixedPos(unsigned npos,const tdouble3* pos);

  //-Computes interpolation (for fixed positions).
  void IntpComputeVar(double timestep,unsigned varnum,unsigned size,float*    result);
  void IntpComputeVar(double timestep,unsigned varnum,unsigned size,tfloat3*  result);
  void IntpComputeVar(double timestep,unsigned varnum,unsigned size,double*   result);
  void IntpComputeVar(double timestep,unsigned varnum,unsigned size,tdouble3* result);

  void IntpComputeVar(double timestep,const std::string& varname,unsigned size,float*    result);
  void IntpComputeVar(double timestep,const std::string& varname,unsigned size,tfloat3*  result);
  void IntpComputeVar(double timestep,const std::string& varname,unsigned size,double*   result);
  void IntpComputeVar(double timestep,const std::string& varname,unsigned size,tdouble3* result);
  
  //-Configures interpolation (for free positions).
  void IntpConfigFreePos();
  void IntpComputeInTime(double timestep);

  void IntpCompute(unsigned np,const tdouble3* pos,float*   result)const;
  void IntpCompute(unsigned np,const tdouble3* pos,tfloat3* result)const;

  float   IntpComputeFloat (tdouble3 pos)const;
  tfloat3 IntpComputeFloat3(tdouble3 pos)const;


  //-Saves output files.
  void SaveVtkScheme(std::string file)const;
  void SaveVtks(std::string file)const;
  void SaveCsv(std::string file,bool csvsepcoma=false)const;
  void SaveBin(std::string file)const;
  static void SaveTemplate(std::string filename="MeshData_Template"
    ,std::string appname="",bool csvsepcoma=false);

};



}

#endif


