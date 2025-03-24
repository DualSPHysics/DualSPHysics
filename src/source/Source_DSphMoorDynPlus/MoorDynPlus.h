/*===================================================================================
<MOORDYNPLUS> Copyright (c) 2025 by Ivan Martinez-Estevez

Ivan Martinez-Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
Matt Hall (github.com/mattEhall)

This file is part of MoorDynPlus. MoorDynPlus is free software: you can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

Linking the MoorDynPlus library statically or dynamically with other modules is
making a combined work based on this library. Thus, the terms and conditions
of the GNU General Public License cover the whole combination. As a special
exception, the copyright holders of MoorDynPlus give you permission to dynamically
link this library with the program DualSPHysics to produce a combined model
featuring the capabilities of both DualSPHysics and MoorDynPlus. This exception
is strictly limited to linking between the compiled MoorDynPlus library and
DualSPHysics. It does not extend to other programs or the use of the MoorDynPlus
source code beyond the stipulations of the GPL. When the exception is used,
this paragraph must be included in the copyright notice.

MoorDynPlus is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.

You should have received a copy of the GNU General Public License along with
MoorDynPlus. If not, see <http://www.gnu.org/licenses/>.
===================================================================================*/

/// \file MoorDynPlus.h \brief Defines the class \ref MoorDynPlus.

#ifndef _MoorDynPlus_
#define _MoorDynPlus_

#include "TypesDef.h"
#include "Functions.h"
#include "JXml.h"
#include "JObject.h"
#include "Misc.h"
#include "IEnvironment.h"
#include "IBody.h"
#include "IOutput.h"
#include "JLog2.h"
#include <stdio.h>
#include <string>

 class MoorDynPlus : protected JObject
{
private:
  JLog2* Log;

  std::string FileIN; ///<Name of input file.
  std::string Dirout; ///<Default IOutput directory (Is replaced when an external software calls to LinesInit() )

  bool Coupled; ///<Used to know if MoorDynPlus is coupled with external software
  
  std::string RootNodeXml; ///<Root element of Xml for read MoorDynPlus configuration
  IEnvironment* Env;       ///<IEnvironment Object
  ILineBase* LineDef;      ///<Pointer to default parameters of the lines
 
  std::vector<IBody*> Bodies;            ///<Array of pointers to Bodies
  std::vector<ILine*> Lines;             ///<Array of pointers to Lines
  std::vector<ILine*> OutLines;          ///<Array of pointers to Lines which will be stored into the csv file
  std::vector<IConnection*> Connections; ///<Array of pointers to Connections
  std::vector<unsigned> References;     ///<Stores all references of Bodies and connects
  
  unsigned BodyCount;       ///<Number of Bodies
  unsigned LineCount;       ///<Number of Lines
  unsigned FairCount;       ///<Number of fairlead connections
  unsigned AnchCount;       ///<Number of anchor connections
  unsigned ConnCount;       ///<Number of "connect" connections
  unsigned ConnectionCount; ///<Number of IConnection objects

  bool XmlReaded;           ///<Used for know if Xml already was read
  bool HeaderProg;          ///<Used to know if header progress was print
  bool FairArrays;          ///<Indicate if the fairlead arrays was initialised.
  bool ItArrays;            ///<Indicate if the integration arrays to store the stats were initialised.
  bool UseWaves;
  bool DsphLog;             ///<True if DualSPHysics sends a JLog2 object, else false
  
  double TimeMax;
  double DtOut;
  double NextTime;
  unsigned Part;
  IOutput* OutP; ///<Pointer to Properties of the outputs

  double*** FairleadPos; ///<Stores link positions.  
  double*** FairleadVel; ///<Stores link velocities. 
  double*** FairleadFor; ///<Stores link forces.     

  //-Static vectors for fairleads
  std::vector<double> FlinesS;          ///<Net line force vector (6-DOF)-retains last solution for use when inputted dt is zero (such as during FAST predictor steps) when the model does not time step
  std::vector<std::vector<double>> rFairRel; ///<Fairlead locations relative to platform ref point but in inertial orientation
  double** rFairtS;          ///<Fairlead locations in turbine/platform coordinates
  double** rFairi;           ///<Fairlead locations in inertial reference frame
  double** rdFairi;          ///<Fairlead velocities in inertial reference frame
  std::vector<int>  FairIs;  ///<Vector of fairlead connection indices in ConnectList vector
  std::vector<int> ConnIs;   ///<Vector of connect connection indices in ConnectList vector
  int NStateVariables;       ///<Used for add six state variables for each "connect"-type IConnection 
  double** Ffair;            ///<Pointer to 2-d array holding fairlead forces
  unsigned NX;               ///<Size of state vector array
  double* X0;                ///<Comprising global state vector

  //-State vectors things for Rk2/rk4 integration 
  double* Xt; ///<Temporal states of positions for RHSmaster
  double* F0; ///<Temporal state of velocities for RHSmaster
  double* F1; ///<Temporal state of velocities for RHSmaster
  //double* F2;
  //double* F3;  
  tmatrix3d TransMat; ///<Calculate TransMat 

  void Reset();
  void Setup();
  void AllocateMemoryIntegration(const unsigned size);
  void AllocateArrays();
  void LoadXml();
  void ReadXml(JXml* xml); 
  void CheckDepth(IEnvironment* env,IBody* body); 
  void CreateOutDirectory(); 
  void PrepareState();
  //void StoreConnPositions(); 
  void CheckReferences(const unsigned ref);
  void CheckConnReferences(const unsigned ref,std::vector<IConnection*> conns); 
  void InitializeSystem(double X[],double XD[]);
  void DoDynamicRelaxIC();
  void AllOutput(double t,double dtC);
  void UpdatePos(const double divisor);
  void PrintProgress(double time); 
  void PrintHeader();
  double PrintPercent(double time);
  void ConfigBodies();
  void ConfigConnections(std::vector<IConnection*> connects); 
  void CheckFtIds(const unsigned ids[],const unsigned numFt)const;
  void CheckBody(const unsigned ftid)const;
  void CheckLine(const unsigned lineid)const;
  void CheckBodyLine(const unsigned ftid,const unsigned lineid)const;
  void VisuConfig()const;

  std::vector<ILine*> GetLines(IBody* body)const; 
  std::vector<ILine*> GetLines(const unsigned ref)const;
  std::vector<ILine*> GetLinesNoBody()const;
  IBody* GetBody(const unsigned ref)const; 
  IBody* GetBody(ILine* line)const; 
  std::vector<IConnection*> GetConnections(const TpConnection type)const;
  void AddAttachedLines();
  void RHSmaster(const double X[],double Xd[],const double t,const double dt); 
  void Rk2(double* t0,const double dt);
  void CheckBodyRefs();
  void RemoveLine(ILine* line); 
  void SaveVtk(double timestep,double dt);
  void VisuIniTen()const;
  std::vector<IConnection*> GetConnections(){return Connections;} ///<Returns a pointer to pointer of all connections

public:
  MoorDynPlus(); 
  virtual ~MoorDynPlus();

  void Run(); 
  int LinesInit(double X[],double XD[],const std::string filexml
    ,const std::string nodexml,const char* dir,const tfloat3 gravity
    ,const double tmax,const double dtout,const unsigned ftmkbound[]=NULL
    ,const unsigned numFts=0);
  int FairleadsCalc(const unsigned numFts,double*** rFairIn,double*** rdFairIn,double*** fFairIn,double* t_in,double* dt_in);
  unsigned GetNFairs(const unsigned ftid)const;
  unsigned GetSegsCount(const unsigned lineid)const;
  unsigned GetSegsCount(const unsigned ftid,const unsigned lineid)const; 
  double GetFairTen(const unsigned lineid)const;
  int GetNodePos(const unsigned lineid,const int nodeid,double pos[3]);
  tdouble3 GetNodePos(const unsigned lineid,const int nodeid)const;
  int GetNodePosLink(const unsigned ftid,const unsigned lineid,double pos[3]);
  tdouble3 GetNodePosLink(const unsigned ftid,const unsigned lineid)const;
  unsigned GetBodyReference(const unsigned ftid=0)const;
  int LinesClose();
  void LogInit(JLog2* log);
  void CreateLog();

  unsigned GetNLines()   {return LineCount; } ///<Returns the number of lines created of a Mooring
  unsigned GetNBodies() { return BodyCount; } ///<Returns the number of Moorings
};
#endif //!MoorDynPlus