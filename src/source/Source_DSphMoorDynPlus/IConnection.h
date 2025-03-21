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

///<\file IConnection.h \brief Defines the class \ref IConnection.

#ifndef _Connection_
#define _Connection_

#include "Misc.h"
#include "JXml.h"
#include "JObject.h"
#include "TypesMoorDynPlus.h"
#include "JLog2.h"

class IEnvironment;
class ILine;

class IConnection : protected JObject
{
private:
  JLog2* Log;

  const TpConnection Type; ///<Defining whether fixed 0,vessel 1,or Connect 2   

  bool Disabled;
  unsigned Node;      ///<Indicates which is its position in the line
  bool UseMass;       ///<Indicates if the Mass was introduced
  bool UseVol;        ///<Indicates if the Volumn was introduced
  unsigned NumExec;   ///<Indicates the number of executions of CalculateMassVol()
  unsigned Ref;       ///<Reference 
  std::vector<ILine*> Attached;///<Pointers to lines attached to this IConnection Node
  std::vector<unsigned> Top;  ///<Which end of line are we attached to? 1=top/Fairlead,0=bottom/anchor
  unsigned Idn;          ///<Identifier
  unsigned AttachCount;  ///<Number of attached lines
  
  tdouble3 Pos;      ///<Initial positions user-defined
  tdouble3 ForceIni; ///<Initial forces user-defined

  double ConM;    ///<Initial Mass
  double ConV;    ///<Initial Velocity
  double ConCdA;  ///<Product of drag coefficient and frontal area
  double ConCa;   ///<Added mass coefficient

  IEnvironment*  Env; ///<Pointer to environmental settings

  //-Common properties with line internal Nodes   
  tdouble3 R;      ///<Node position [x/y/z] // double* to continuous state or input state or parameter/other-state
  tdouble3 Rd;     ///<Node velocity[x/y/z]  // double* to continuous state or input state or parameter/other-state
  tdouble3 Q;      ///<Unit tangent vector for end Node (unused)
  tdouble3 R_ves;  ///<Fairlead position for vessel Node types [x/y/z]
  tdouble3 Rd_ves; ///<Fairlead velocity for vessel Node  types [x/y/z]
  tdouble3 Fnet;   ///<Total force on Node
  tdouble3 Fnet_i; ///<
  tdouble3 RHS;    ///<RHS of state-space equation (Forces divided by mass matrix)
  tmatrix3d S;     ///<Inverse mass matrices (3x3) for each Node
  tmatrix3d M;     ///<Node mass+added mass matrices
  tmatrix3d M_i;   ///<New Mass calculated 
  
  double TimeStep; ///<Simulation time
  double T0;       ///<Simulation time current integration was started at (used for BC function)
  
  void ReadXml(JXml* sxml,TiXmlElement* lis,const int id); 
  void Reset(); 
  void CalculateMassVol();

public:
  IConnection(JLog2* log,TpConnection tpc);
  virtual ~IConnection();

  void Setup();
  void GetConnectState(tdouble3& r_out,tdouble3& rd_out){r_out=R; rd_out=Rd;} ///<Function to return connection position and velocity to ILine object
  tdouble3 GetFnet()const{return Fnet;} ///<Function to return net force on fairlead (just to allow public reading of Fnet)
  void InitializeFairlead(const tdouble3 pX,const tmatrix3d& transMat);
  void InitializeConnect(double* X);
  void GetNetForceAndMass();
  void DoRHS(const double* X,double* Xd,const double time);
  void InitiateStep(const tdouble3 rFairIn,const tdouble3 rdFairIn,const double time);
  void UpdateFairlead(const double time);
  void LoadXml(JXml* sxml,const std::string& place,TiXmlElement*  elec,const int id);
  void AddLineToConnect(ILine* line);
  tdouble3 GetPositions()const {return R;}  ///<Return a TDouble3 with th positions of connection
  void SetPositions(const tdouble3& pos) {  R=pos;} ///<Stores the new positions [X Y Z]
  void ResetTension(){Fnet=TDouble3(0);}    ///<Restores the tension for the IConnection when it is disabled; 
  
  unsigned     GetId()const       { return Idn;    } ///<Returns the number of IConnection
  unsigned     GetNode()const     { return Node;   } ///<Returns the node position
  unsigned     GetRef()const      { return Ref;    } ///<Returns the reference of connect
  TpConnection GetType()const     { return Type;   } ///<Returns the connection Type
  double       GetX()const        { return Pos.x;  } ///<Returns the initial position in X
  double       GetY()const        { return Pos.y;  } ///<Returns the initial position in Y
  double       GetZ()const        { return Pos.z;  } ///<Returns the initial position in Z
  bool         GetDisabled()const {return Disabled;} ///<Returns true if the line is broken,else returns false
  std::string  GetTypeStr()const; 

  void SetId      (const unsigned idn)    { Idn=idn;     } ///<Stores the number id
  void SetNode    (const unsigned node)   { Node=node;   } ///<Stores the node position
  void SetRef     (unsigned ref)          { Ref=ref;     } ///<Stores the reference of connect
  void SetEnv     (IEnvironment* env_in)   { Env=env_in;  } ///<Stores a pointer to IEnvironment settings
  void SetDisabled(const bool b)          { Disabled=b;  } ///<Stores the new value for Disabled variable

  bool IsVessel() const{return (Type==TpConnection::vessel) ;}
  bool IsFixed()  const{return (Type==TpConnection::fixed)  ;}
  bool IsConnect()const{return (Type==TpConnection::connect);}

  bool operator==(IConnection* conn)const{ return ((conn->GetId()==Idn)&&(conn->GetType()==Type)); } ///<Compare two connections and return true when the Number and Type are matching

  void VisuProperties()const;

#ifdef USEGL
  void drawGL(void);
#endif
};

#endif
