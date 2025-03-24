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

/// \file ILine.h \brief Defines the class \ref ILine.

#ifndef _Line_
#define _Line_

#include "Misc.h"
#include "JXml.h"
#include "JObject.h"
#include "JLog2.h"
#include "TypesMoorDynPlus.h"

#include <limits.h>
#include <float.h>

class IBody;
class IConnection;
class IEnvironment;
class IQSlines;

//********************************************************************************
/// Class ILineBase
//********************************************************************************
class ILineBase : protected JObject
{

private:
  void Reset(); /// Restores the attributes of the class
  void ReadXml(JXml* sxml,TiXmlElement* lis,const std::string& place); ///Reads the input file

public:
  double E;            ///<Stiffness (Young's modulus) [N/m^2]
  double A;            ///<Area of ILine [m^2]
  double Ea;           ///<Stiffness of the line; product of elasticity modulus and cross-sectional area [N]
  double D;            ///<Diameter of ILine
  double Can;          ///<Default=1.0
  double Cat;          ///<Default=0.0
  double Cdn;          ///<Default=1.6
  double Cdt;          ///<Default=0.05
  double MassDenInAir; ///<Linear weight in air
  double C;            ///<The line internal damping (Ns)
  double BkTension;    ///<Maximum value of tension allowed
  unsigned N;          ///<Number of line segments 
  double UnstrLen;     ///<Length of line

  std::string Channels;///<Node output properties
  
  ILineBase(); 
  virtual ~ILineBase();

  void LoadXml(JXml* sxml,const std::string& place);
  void CheckChannels(JXml* sxml,TiXmlElement*  elel);
};

//********************************************************************************
/// Class ILine
//********************************************************************************
class ILine : protected ILineBase
{
  // here is the numbering scheme (N segments per line):
  // [connect (node 0)]  --- segment 0 --- [ node 1 ] --- seg 1 --- [node2] --- ... --- seg n-2 --- [node n-1] --- seg n-1 ---  [connect (node N)]
private:
  double DtOut;        ///<Time step for write
  bool OutputLine;     ///<Indticates if this line will be stored in the csv file. default=true
  bool OutputNode;     ///<Indticates if the nodes of this line will be stored in the csv files. default=true
  std::string DirOut;

  bool Disabled;       ///<Indicates if the line was broken
  bool ConnectFlag;    ///<Indicates if the line has a connect connection
  bool FairFlag;       ///<Indicates if the line has a fair connection

  JLog2* Log;
  IQSlines* Qslines;

  double Depth;          ///<Depth of the water (m).
  double Zmin;           ///<Gets the minimum position at z for contacts with ground
  std::string FileName; ///<Pointer to file output name
  IEnvironment* Env;      ///<Pointer to environmental settings
  unsigned* BodyId;      ///<Id of body and Floating object (used when there is a coupling with external software)
  std::vector<tdouble3> R;  ///<Node positions [i][x/y/z]
  std::vector<tdouble3> Rd; ///<Node velocities [i][x/y/z]
  std::vector<tdouble3> Q;  ///<Unit tangent vectors for each segment

  std::vector<IConnection*> Connections; ///<Array of Connections
  IConnection* AnchConnect; ///<Pointer to anchor connection
  IConnection* FairConnect; ///<pointer to fairlead connection

  unsigned Idl;             ///<ILine "number" id
  unsigned ConnectionCount; ///<Total number of IConnection objects
  unsigned FairCount;       ///<Number of fairlead connections
  unsigned AnchCount;       ///<Number of anchor connections
  unsigned ConnCount;       ///<Number of "connect" connections
  unsigned StartIndex;      ///<Start index for each node for States array 
  bool UseWaves;            ///<Flag indicating whether wave kinematics will be considered for this line

  double TimeStep; ///<Simulation time
  double T0;       ///<Simulation time current integration was started at (used for BC function)

  // declaring these here as they caused problems if declared in the loop
  tdouble3 vi; ///<relative velocity
  tdouble3 vp; ///<transverse component of relative velocity
  tdouble3 vq; ///<axial component of relative velocity
  tdouble3 ap; ///<transverse component of absolute acceleration-HAD TO MOVE THIS UP FROM LOWER DOWN (USED TO CRASH AFTER I=3)
  tdouble3 aq; ///<axial component of absolute acceleration

  // forces 
  std::vector<tdouble3> T;  ///<
  std::vector<tdouble3> Td; ///<
  std::vector<double> Tmag; ///<Segment tension magnitude 
  std::vector<tdouble3> W;  ///<Node weight 

  std::vector<tdouble3> Dp; ///<Node drag (transverse)
  std::vector<tdouble3> Dq; ///<Node drag (axial)
  std::vector<tdouble3> Ap; ///<Node added mass forcing (transverse)
  std::vector<tdouble3> Aq; ///<Node added mass forcing (axial)
  std::vector<tdouble3> B; ///<Node bottom contact force

  std::vector<tdouble3> Fnet; ///<total force on node

  std::vector<tmatrix3d> S; ///<Inverse mass matrices (3x3) for each node
  std::vector<tmatrix3d> M; ///<Node mass+added mass matrix

  std::vector<double> F;    ///<VOF scalar for each segment (1=fully submerged,0=out of water)
  std::vector<double> L;    ///< ILine unstretched segment lengths
  std::vector<double> Lstr; ///<Stretched lengths (m)
  std::vector<double> Ldstr;///<Rate of stretch (m/s)

  double Rho;    ///<line density
  double ReFac;  ///<
  double MassDenInAir;   ///<Linear weight in air [ w in ILineBase ]

  std::vector<double> V; ///<ILine segment volume

  //-set up output arrays,at each node i:
  std::vector<tdouble3> U;  ///<Wave velocities  
  std::vector<tdouble3> Ud; ///<Wave accelerations
  std::vector<double> Zeta; ///<Free surface elevation

  //-file stuff
  //-new additions for handling waves in-object and precalculating them  (not necessarily used right now)
  int WaveMod;   ///< 
  int WaveStMod; ///<
  double Hs;     ///< 
  double Tp;     ///<
  double Gamma;  ///<
  float Beta;    ///<Wave heading

  int Nw;               ///<Number of wave frequency components considered    //OK AS INT???
  std::vector<float> w; ///<
  std::vector<float> k; ///<
  float Dw;             ///<FAST's dw (really small typically)

  std::vector<misc::floatC> ZetaC0;           ///<Fourier transform of wave elevation at origin
  std::vector<misc::floatC> ZetaC;            ///<Fourier transform of wave elevation
  std::vector<std::vector<misc::floatC>> UC;  ///<Fourier transform of wave velocities
  std::vector<std::vector<misc::floatC>> UdC; ///<Fourier transform of wave accelerations

  std::vector<misc::doubleC> WGNC;///<
  std::vector<double> WaveS2Sdd;  ///<
  misc::doubleC WGNC_Fact;        ///<This factor is needed by the discrete time inverse Fourier transform to ensure that the time series WGN process has unit variance // sqrt( PI/(dw*WaveDT) );   
  double S2Sd_Fact;               ///<This factor is also needed by the discrete time inverse Fourier transform // 1.0/WaveDT;
  double*  Ucurrent;              ///<Constant uniform current to add (three components)

  // new additions for precalculating wave quantities
  std::vector<std::vector<double>> zetaTS;             ///<Time series of wave elevations above each node
  std::vector<std::vector<double>> FTS;                ///<
  std::vector<std::vector< std::vector<double>> > UTS; ///<
  std::vector<std::vector< std::vector<double>> > UdTS;///<
  int Nt;          ///<Number of wave time steps
  double WaveDT;   ///<Wave time step size (s)
  int Ts0;         ///<Time step index used for interpolating wave kinematics time series data (put here so it's persistent)
  std::vector<double> tTS; ///<Time step vector

  void Reset(); 
  void ReadXml(JXml* sxml,const std::string& place,TiXmlElement* eleL
    ,const unsigned lineid,std::vector<IConnection*> connects,ILineBase* lineDef=NULL); 
  void AllocateMemoryVectors();
  void CalculateMassMatrix();
  void CheckInputParameters();
  void CheckDtMRelationship();
  void ConfigLine(std::vector<IConnection*> connections);
  void CopyFrom(const ILineBase* lineDef);
  double CalculateVOFSeg(double n1z,double n2z,double hswl)const;

  IConnection* SearchConnect(const unsigned ref,std::vector<IConnection*> connects); 

public:
  ILine(JLog2* log);
  virtual ~ILine();

  void LoadXml(JXml* sxml,const std::string& place,TiXmlElement* eleL
    ,const unsigned num_tag,std::vector<IConnection*> connects,ILineBase* lineDef=NULL);
  void Setup(std::string& dir,IEnvironment* env_in,IBody* body=NULL);
  void SetupWaves(double* Ucurrent_in,float dt_in);
  void Initialize(double* X);
  
  double   GetNodeTen       (unsigned i)const; 
  tdouble3 GetNodePos       (unsigned i)const;
  unsigned GetNConns        (const unsigned ref)const;
  double   GetTensionOutput (const unsigned n)const; 
  tdouble3 GetPositionOutput(const unsigned n)const;
  tdouble3 GetVelocityOutput(const unsigned n)const;
  tdouble3 GetForceOutput   (const unsigned n)const; 

  std::vector<IConnection*> GetConnections(const TpConnection type)const; 
  std::vector<IConnection*> GetFairs(const unsigned ref)const;
  std::vector<IConnection*> GetConnConnects()const;
  std::vector<IConnection*> GetConnConnects(const unsigned ref)const; 

  void ScaleDrag(double scaler);
  void DoRHS(const double* X,double* Xd,const double time,const double dt);
  void WriteOutput(double timestep,double dt,double start,double finish,double next);
  void GetFASTtens(float* FairHTen,float* FairVTen,float* AnchHTen,float* AnchVTen); 
  void GetAnchStuff(tdouble3 & fnet,tmatrix3d& mass); 
  void GetFairStuff(tdouble3 & fnet,tmatrix3d& mass);
  void BreakLine();
  void CheckTension();
  void ConfigDepth(const IEnvironment* env, const IBody* b);
  void VisuConfig()const;
  void VisuIniTen()const;

  IConnection* GetAnchConnect()const  {return AnchConnect;} ///<Returns a pointer to Anch connect
  IConnection* GetFairConnect()const  {return FairConnect;} ///<Returns a pointer to Fair connect

  double   GetRho()const                 { return Rho;             } ///<Returns the density
  unsigned GetNConnections()const        { return ConnectionCount; } ///<Returns the number of connections
  unsigned GetNFairs()const              { return FairCount;       } ///<Returns the number of Fairleads
  unsigned GetNAnchs()const              { return AnchCount;       } ///<Returns the number of Anchors
  unsigned GetNConns()const              { return ConnCount;       } ///<Returns the number of Connects
  unsigned GetId()const                  { return Idl;             } ///<Returns the number of ILine
  unsigned GetN()const                   { return N;               } ///<Returns N (number of segments)  
  unsigned GetIndex()const               { return StartIndex;      } ///<Returns the Start Index
  double   GetLength()const              { return UnstrLen;        } ///<Returns the length
  bool     GetDisabled()const            { return Disabled;        } ///<Returns the break tension value
  double   GetBreakTension()const        { return BkTension;       } ///<Returns true if the line is broken,else returns false
  bool     GetUseConnect()const          { return ConnectFlag;     } ///<Returns true if this line has a connect connection
  bool     GetUseFair()const             { return FairFlag;        } ///<Returns true if this line has a fair connection
  bool     GetUseOutput()const           { return OutputLine;      } ///<Returns true if will store data of this line
  double   GetM(const unsigned seg)const { return M[seg].a11;      } ///<Returns the Mass of one segment
  double   GetV(const unsigned seg)const { return V[seg];          } ///<Returns the Volumn of one segment

  void SetIndex(const unsigned index){ StartIndex=index; } ///<Stores the new Index
  void SetDisabled(const bool b)     { Disabled=b;       } ///<Stores the new value for Disabled variable
  void SetDtOut(const double d)      { DtOut=d;          } ///<Stores the new value for time step to write node properties
  void SetTime(double time)          { TimeStep=time;    } /// Function to reset time after IC generation

  std::vector<IConnection*> GetConnections()const{return Connections;} ///<Returns the array of Connections
 
  void VisuProperties(unsigned ini=0,unsigned end=UINT_MAX)const;

#ifdef USEGL  
  void drawGL(void);
  void drawGL2(void);
#endif
};

#endif


