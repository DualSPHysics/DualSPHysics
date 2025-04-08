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

/// \file ILine.cpp \brief Implements the class \ref ILine.

#include "ILine.h"
#include "IEnvironment.h"
#include "IConnection.h"
#include "IBody.h"
#include "IQSlines.h" // the c++ version of quasi-static model Catenary
#include "FunMoorDynPlus.h"
#include "Functions.h"
#include "JSaveCsv2.h"
#include "FunGeo3d.h"
#include "FunctionsMath.h"
#include <string>
#include <algorithm>


//*******************************************************************************
//************************** ILineBase Section ***********************************
//*******************************************************************************

//==============================================================================
/// Constructor
//==============================================================================
ILineBase::ILineBase() {
  Reset();
  ClassName="ILineBase";
}
//==============================================================================
/// Destructor
//==============================================================================
ILineBase::~ILineBase() {
  //Reset();
  A=0;
  E=0;
  Ea=0;
  D=0;
  Can=0;
  Cat=0;
  Cdn=0;
  Cdt=0;
  MassDenInAir=0;
  UnstrLen=0;
  N=0;
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void ILineBase::Reset() {
  A=DBL_MAX;
  E=DBL_MAX;
  Ea=DBL_MAX;
  D=DBL_MAX;
  Can=DBL_MAX;
  Cat=DBL_MAX;
  Cdn=DBL_MAX;
  Cdt=DBL_MAX;
  MassDenInAir=DBL_MAX;
  C=DBL_MAX;
  Channels="";
  BkTension=0.0;
  UnstrLen=DBL_MAX;
  N=0;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void ILineBase::LoadXml(JXml* sxml,const std::string& place) {
  std::string function="LoadXml";
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node )Run_Exceptioon("Cannot find the element \'"+place+"\'.");
  ReadXml(sxml,node->ToElement(),place);
}
//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void ILineBase::ReadXml(JXml* sxml,TiXmlElement* lis,const std::string& place) {
  //-Loads linedefault zone.
  TiXmlElement* elelDef=lis->FirstChildElement("linedefault");
    
  sxml->CheckElementNames(elelDef,false,"ea e diameter massDenInAir ba can cat cdn cdt breaktension length segments outputFlags");

  // If wants overwrite default parameters
  D=sxml->ReadElementDouble(elelDef,"diameter","value",false);
  Ea=sxml->ReadElementDouble(elelDef,"ea","value",true,DBL_MAX);
  E=sxml->ReadElementDouble(elelDef,"e","value",true,DBL_MAX);
  C=sxml->ReadElementDouble(elelDef,"ba","value",true,-0.8);
  Can=sxml->ReadElementDouble(elelDef,"can","value",true,1.0);
  Cat=sxml->ReadElementDouble(elelDef,"cat","value",true,0.0);
  Cdn=sxml->ReadElementDouble(elelDef,"cdn","value",true,1.6);
  Cdt=sxml->ReadElementDouble(elelDef,"cdt","value",true,0.05);
  MassDenInAir=sxml->ReadElementDouble(elelDef,"massDenInAir","value",false);
  BkTension=sxml->ReadElementDouble(elelDef,"breaktension","value",true,0.0); 
  UnstrLen=sxml->ReadElementDouble(elelDef,"length","value",true,DBL_MAX);
  N=sxml->ReadElementUnsigned(elelDef,"segments","value",true,0); 
  Channels=sxml->ReadElementStr(elelDef,"outputFlags","value",true,"-");

  if(Ea==DBL_MAX && E==DBL_MAX) Run_ExceptioonFile(fun::PrintStr("The parameter <ea> or <e> should be introduced. Just define one of them."),sxml->ErrGetFileRow(elelDef));
  if(Ea!=DBL_MAX && E!=DBL_MAX) Run_ExceptioonFile(fun::PrintStr("The parameters <ea> and <e> cannot be defined at the same time. Just define one of them."),sxml->ErrGetFileRow(elelDef));

  CheckChannels(sxml,elelDef);
}

//==============================================================================
/// Check if the output flag is in the list of possibilities.
//==============================================================================
void ILineBase::CheckChannels(JXml* sxml,TiXmlElement* eleLine) {
  std::string function="CheckChannels";

  bool equals=false;
  std::string chandef=OutChannels;
  std::string chaninp=Channels;
 
  for(unsigned i=0; i<chaninp.length(); i++) { //for each input channels insert by user
    for(unsigned j=0; j<chandef.length(); j++) { //for each char of nodes
      if(chandef[j]==chaninp[i]) { //If both of them are equals
        equals=true;
        break;
      }
    }
    if(!equals)Run_Exceptioon(fun::PrintStr("Value \'%c\' in <ouputsFlags> is not allowed.",chaninp[i]));
    equals=false;
  }
}

//*******************************************************************************
//****************************** ILine Section ***********************************
//*******************************************************************************

// here is the new numbering scheme (N segments per line)
//   [connect (node 0)]  --- segment 0 --- [ node 1 ] --- seg 1 --- [node2] --- ... --- seg n-2 --- [node n-1] --- seg n-1 ---  [connect (node N)]

//==============================================================================
/// Constructor
//==============================================================================
ILine::ILine(JLog2* log):Log(log) {
  Reset();
  Qslines=new IQSlines(log);
  ClassName="ILine";
}
//==============================================================================
/// Destructor.
//==============================================================================
ILine::~ILine() {
  Log=NULL;
  delete Qslines;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void ILine::Reset() {
  Connections.clear();
  ConnectionCount=0; // total Idl of IConnection objects
  FairCount=0;       // Number of fairlead connections
  AnchCount=0;       // Number of anchor connections
  ConnCount=0;       // Number of "connect" connections
  UseWaves=false;    // start off with wave kinematics disabled.  Can be enabled after initial conditions are found and wave kinematics are calculated  
  StartIndex=0;
  ConnectFlag=false;
  FairFlag=false;
  OutputLine=false;
  OutputNode=false;
  BodyId=NULL;
  FileName="";
  Env=NULL;
  AnchConnect=NULL;
  FairConnect=NULL;
  Qslines=NULL;
  BkTension=0.0;
  Disabled=false;
  DtOut=0.0;
  Zmin=0;
  Depth=DBL_MAX;
  DirOut="";
}

//==============================================================================
/// Allocates memory of vectors
//==============================================================================
void ILine::AllocateMemoryVectors() {
  // =============== size vectors =========================
  R.resize(N+1 ,TDouble3(0)); //< node positions [i][x/y/z]
  Rd.resize(N+1,TDouble3(0)); //< node velocities [i][x/y/z]
  Q.resize(N+1 ,TDouble3(0)); //< unit tangent vectors for each node

  // forces 
  T.resize(N,TDouble3(0));  // line tensions
  Td.resize(N,TDouble3(0));   // line damping forces
  //  Tmag.resize(N,0.0);        // segment tension magnitudes << hardly used
  W.resize(N+1,TDouble3(0));  // node weights

  Dp.resize(N+1,TDouble3(0));    // node drag (transverse)
  Dq.resize(N+1,TDouble3(0));    // node drag (axial)
  Ap.resize(N+1,TDouble3(0));    // node added mass forcing (transverse)
  Aq.resize(N+1,TDouble3(0));    // node added mass forcing (axial)
  B.resize(N+1,TDouble3(0));    // node bottom contact force
  Fnet.resize(N+1,TDouble3(0));  // total force on node

  S.resize(N+1,TMatrix3d(0));  // inverse mass matrices (3x3) for each node
  M.resize(N+1,TMatrix3d(0));  // mass matrices (3x3) for each node

  L.resize(N,0.0);       // line unstretched segment lengths
  Lstr.resize(N,0.0);     // stretched lengths
  Ldstr.resize(N,0.0);     // rate of stretch
  V.resize(N,0.0);      // volume?

  Zeta.resize(N+1,0.0);          // wave elevation above each node
  F.resize(N+1,0.0);   // fixed 2014-12-07  // VOF scalar for each NODE (mean of two half adjacent segments) (1=fully submerged,0=out of water)
  U.resize(N+1,TDouble3(0));     // wave velocities
  Ud.resize(N+1,TDouble3(0));;   // wave accelerations
}

//==============================================================================
/// Copy the parameters from ILineBase.
//==============================================================================
void ILine::CopyFrom(const ILineBase* lineDef){
  D=lineDef->D;
  E=lineDef->E;
  Ea=lineDef->Ea;
  C=lineDef->C;
  Can=lineDef->Can;
  Cat=lineDef->Cat;
  Cdn=lineDef->Cdn;
  Cdt=lineDef->Cdt;
  Channels=lineDef->Channels;
  MassDenInAir=lineDef->MassDenInAir;
  BkTension=lineDef->BkTension;
  UnstrLen=lineDef->UnstrLen;
  N=lineDef->N;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void ILine::LoadXml(JXml* sxml,const std::string& place,TiXmlElement* eleL,const unsigned lineid,std::vector<IConnection*> connects,ILineBase* lineDef) {
  std::string function="LoadXml";  
  //TiXmlNode* node=sxml->GetNode(place,false);
  //if(!node)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
  if(!eleL)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
  ReadXml(sxml,place,eleL,lineid,connects,lineDef);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void ILine::ReadXml(JXml* sxml,const std::string& place_in,TiXmlElement* eleLine,const unsigned lineid,std::vector<IConnection*> connects,ILineBase* lineDef) {
  std::string function="ReadXml";
  std::string place=place_in;
  
  Idl=lineid; // stores a secuencial Number of lines in the program for each body
  
  sxml->CheckElementNames(eleLine,false,"ea e diameter massDenInAir ba can cat cdn cdt breaktension outputFlags length segments depth *fix *fixconnection *vessel *vesselconnection *connect");

  //-If there is line tag
  if(eleLine) { 
    //-If there are default parameters
    if(lineDef) CopyFrom(lineDef);
  
    //-If exists,then overwrites the default parameters
    if(sxml->ExistsElement(eleLine,"ea"))           { Ea=sxml->ReadElementDouble(eleLine,"ea","value",false); }
    if(sxml->ExistsElement(eleLine,"e"))            { E=sxml->ReadElementDouble(eleLine,"e","value",false); }
    if(sxml->ExistsElement(eleLine,"massDenInAir")) { MassDenInAir=sxml->ReadElementDouble(eleLine,"massDenInAir","value",false); }
    if(sxml->ExistsElement(eleLine,"diameter"))     { D=sxml->ReadElementDouble(eleLine,"diameter","value",false); }
    if(sxml->ExistsElement(eleLine,"ba"))           { C=sxml->ReadElementDouble(eleLine,"ba","value",false); }
    if(sxml->ExistsElement(eleLine,"can"))          { Can=sxml->ReadElementDouble(eleLine,"can","value",false); }
    if(sxml->ExistsElement(eleLine,"cat"))          { Cat=sxml->ReadElementDouble(eleLine,"cat","value",false); }
    if(sxml->ExistsElement(eleLine,"cdn"))          { Cdn=sxml->ReadElementDouble(eleLine,"cdn","value",false); }
    if(sxml->ExistsElement(eleLine,"cdt"))          { Cdt=sxml->ReadElementDouble(eleLine,"cdt","value",false); }
    if(sxml->ExistsElement(eleLine,"outputFlags"))  { Channels=sxml->ReadElementStr(eleLine,"outputFlags","value",false); }
    if(sxml->ExistsElement(eleLine,"breaktension")) { BkTension=sxml->ReadElementDouble(eleLine,"breaktension","value",false); }
    if(sxml->ExistsElement(eleLine, "depth"))       { Depth=sxml->ReadElementDouble(eleLine,"depth","value",true,DBL_MAX); }
    if(sxml->ExistsElement(eleLine, "length"))      { UnstrLen=sxml->ReadElementDouble(eleLine,"length","value",true,DBL_MAX); }
    if(sxml->ExistsElement(eleLine, "segments"))    { N=sxml->ReadElementUnsigned(eleLine,"segments","value",true,0); }

    if(UnstrLen==DBL_MAX)Run_Exceptioon(fun::PrintStr("The unstretched length cannot be 0 (line_id=%d)",Idl));
    if(N==0)             Run_Exceptioon(fun::PrintStr("The number of segments cannot be 0 (line_id=%d)",Idl));

    //-Checks if the output Channels are right
    CheckChannels(sxml,eleLine);

    //-Stores the elements of the line
    OutputLine=sxml->GetAttributeBool(eleLine,"savedata",true,true);

    //PRINT: look num of line executed
    if(Debug) Log->PrintfDbg("\tReadLine()-num: %d\n",Idl);

    //-Reads the connections of the ILine
    std::vector<IConnection*> connections_temp;
    ConnectionCount=0;
    TiXmlElement* ele=eleLine->FirstChildElement(); 
    while(ele && ConnectionCount<=2){
      const std::string elename=ele->Value();
        if(elename.length()>2 && (fun::StrLower(elename.substr(0,3))=="fix") && sxml->CheckElementActive(ele)){
        IConnection* c=new IConnection(Log,TpConnection::fixed);
        c->LoadXml(sxml,place,ele,(unsigned)connections_temp.size());//Build Fix IConnection
        connections_temp.push_back(c); //Pointer to Fix connection
        AnchCount++; ConnectionCount++;
      }else if(elename.length()>5 && (fun::StrLower(elename.substr(0,6))=="vessel") && sxml->CheckElementActive(ele)){
        IConnection* c=new IConnection(Log,TpConnection::vessel);
        c->LoadXml(sxml,place,ele,(unsigned)connections_temp.size()); //Build VesselConnection
        connections_temp.push_back(c); //Pointer to Vessel connection
        FairFlag=true; FairCount++; ConnectionCount++;
      }else if(elename.length()>6 && (fun::StrLower(elename)=="connect") && sxml->CheckElementActive(ele)){
        if(connects.size()>0) {
          IConnection* c=SearchConnect(sxml->GetAttributeInt(ele,"conref",false),connects);
          connections_temp.push_back(c);
          ConnectFlag=true; ConnCount++;  ConnectionCount++;
        }else Run_Exceptioon("There are no \'connect\' connections declared. Define the section \'<connects>\' .");
      }
      ele=ele->NextSiblingElement();
    }
    
    //-Stores the number of connections read
    ConnectionCount=(unsigned)connections_temp.size();
    if(ConnectionCount!=2) Run_Exceptioon(fun::PrintStr("Only 2 connections are allowed. Found: %d (line_id=%d)",ConnectionCount,Idl));
    ConfigLine(connections_temp);
    connections_temp.clear();

    //PRINT: look Number and types of Connections
    if(Debug) {
      for(unsigned i=0; i<ConnectionCount; i++) {
        Log->PrintfDbg("\t\tconnection: %d Type: %d\n",connections_temp[i]->GetId(),connections_temp[i]->GetType());
        Log->PrintfDbg("\t\t\tx= %.4f ; y= %.4f ; z= %.4f\n",connections_temp[i]->GetX(),connections_temp[i]->GetY(),connections_temp[i]->GetZ());
      }
    }
    CheckInputParameters();
  }
}
//==============================================================================
/// Configures the start and end IConnection
//==============================================================================
void ILine::ConfigLine(vector<IConnection*> connections){
  std::string function="ConfigLine";
  IConnection* con0=connections[0];
  IConnection* con1=connections[1];

  if     (con0->IsFixed()) { AnchConnect=con0; FairConnect=con1; }
  else if(con1->IsFixed()) { AnchConnect=con1; FairConnect=con0; }

  if((con0->IsVessel() && con1->IsVessel()) || (con0->IsConnect() && con1->IsConnect())){ 
    if(con0->GetZ()<=con1->GetZ()) { AnchConnect=con0; FairConnect=con1; }
    else                           { AnchConnect=con1; FairConnect=con0; }
  }
  else if(con0->IsVessel()){ FairConnect=con0; AnchConnect=con1;}
  else if(con1->IsVessel()){ FairConnect=con1; AnchConnect=con0;}
  
  if (!AnchConnect || !FairConnect) Run_Exceptioon(fun::PrintStr("A connection is missing (line_id=%d)",Idl));
  //-Sets the id number of connection
  //AnchConnect->SetId(0);  
  //FairConnect->SetId(1);  

  //-Sets the node position
  AnchConnect->SetNode(0);
  FairConnect->SetNode(N);
  //-Stores the Connections
  Connections.clear();
  Connections.push_back(AnchConnect); //assign Start IConnection to the array
  Connections.push_back(FairConnect); //assign End IConnection to the array
}

//====================================================================================
/// Checks if the connect IConnection selected in Lines was configurated previously
//====================================================================================
IConnection* ILine::SearchConnect(const unsigned ref,std::vector<IConnection*> connects){
  std::string function="SearchConnect";
  IConnection* conn=NULL;

  for(unsigned c=0;c<connects.size()&&!conn;c++)
    if(connects[c]->GetRef()==ref)conn=connects[c];
  
  if(!conn)Run_Exceptioon("The connect with ref: "+fmdp::ToString(ref)+" is missing in <connects>." );
 
  return (conn);
}

//==============================================================================
/// Checks if the input parameters are into the range
//==============================================================================
void ILine::CheckInputParameters() {
  std::string function="CheckInputParameters";
  
  //-Checks distance between 2 points
  const double distance=sqrt((pow((abs(AnchConnect->GetX())-abs(FairConnect->GetX())),2))
   +(pow((abs(AnchConnect->GetY())-abs(FairConnect->GetY())),2))
   +(pow((abs(AnchConnect->GetZ())-abs(FairConnect->GetZ())),2))
  ); //distance between nodeAnch and nodeFair
  
  if(UnstrLen<distance) { 
    string tex="The unstretched length of the line "+fmdp::ToString(Idl) 
     +" should be at least of "+fmdp::ToString(distance)+"m.";
    //Run_Exceptioon(tex); 
    Log->PrintfWarning("%s",tex.c_str());
  }
}

//==============================================================================
/// Returns a pointer to Connections by type
//==============================================================================
std::vector<IConnection*> ILine::GetConnections(const TpConnection type)const{
  std::vector<IConnection*> conns;
  for(unsigned c=0; c<ConnectionCount; c++) {
    if(Connections[c]->GetType()==type) {
      conns.push_back(Connections[c]);
    }
  }
  return conns;
}

//==============================================================================
/// Returns the array of Connections by type
//==============================================================================
std::vector<IConnection*> ILine::GetFairs(const unsigned ref)const{
  std::vector<IConnection*> conns;
  for(unsigned c=0; c<ConnectionCount; c++) {
    if( (Connections[c]->IsVessel())
      &&(Connections[c]->GetRef()==ref) ) {
      conns.push_back(Connections[c]);
    }
  }
  return conns;
}

//==============================================================================
/// Returns a pointer to Connecs connect
//==============================================================================
std::vector<IConnection*> ILine::GetConnConnects()const{
  std::vector<IConnection*> conns;
  for(unsigned c=0; c<ConnectionCount; c++) {
    if(Connections[c]->IsConnect()) {
      conns.push_back(Connections[c]);
    }
  }
  return conns;
}

//==============================================================================
/// Returns the number of Connects by reference
//=============================================================================
unsigned ILine::GetNConns(const unsigned ref)const{
  std::vector<IConnection*> conns=GetConnConnects();
  int nc=0;
  for(unsigned c=0; c<conns.size(); c++) {
    if(conns[c]->GetRef()==ref) nc++;
  }
  return nc;
}

//==============================================================================
/// Returns the "Connects" by reference
//==============================================================================
std::vector<IConnection*> ILine::GetConnConnects(const unsigned ref) const{
  std::vector<IConnection*> conns;
  for(unsigned c=0; c<ConnectionCount; c++) {
    if(Connections[c]->GetRef()==ref) conns.push_back(Connections[c]);
  }
  return conns;
}

//==============================================================================
/// Check if the relationship between dtM, N and EA is right
//==============================================================================
void ILine::CheckDtMRelationship() {
  const double dtini=Env->GetDtM0();
  const double mass_w=abs(Rho-Env->GetRho_w())*A; //Mass in water per length
  const double fn_rads=((2*N)/UnstrLen)*(sqrt(Ea/mass_w));
  const double fn_eq=fmdp::RadToHertz(fn_rads)*10; //Frequency in Hertz
  const double dtnew=(1/fn_eq)/10; // Decrease dt to guarantee the working system
  Env->SetDtM0(min(dtini,dtnew));
}

//==============================================================================
///  Makes the setup for this ILine. Initializes this and creates output file
//==============================================================================
void ILine::Setup(std::string& dir,IEnvironment* env_in,IBody* body_in) {
  Env=env_in;
  BodyId=NULL;
  DirOut=dir;
  //-If this line is connected to a body
  if(body_in)BodyId=body_in->GetPtrRef();

  ConfigDepth(env_in,body_in);
  
  Zmin=(Env->GetFreeSurface()-(Depth)); ///<Gets the minimum position at z for contacts with ground
  A=(PI/4.*D*D); //<Area 

  //PRINT: look when ILine::Setup() is executed
  if(Debug) Log->PrintfDbg("\tLine::setup()-Number: %d\n",Idl);

  // ================== set up properties ===========  

  //-Configures ILine parameters
  Rho=MassDenInAir/A;
  if(E==DBL_MAX)      E=Ea/A; //<Stiffness per area (Young's modulus) [N/m^2]
  else if(Ea==DBL_MAX)Ea=E*A; //<Stiffness of the line; product of elasticity modulus and cross-sectional area [N]

  const double temp_c=C; //Stores BA/-Zeta – the line internal damping (Ns)
  C=C/A;

  if(Env->GetDtMauto())CheckDtMRelationship();
  AllocateMemoryVectors();

  //- Sets the wave elevation at the free surface
  std::fill(Zeta.begin(),Zeta.end(),Env->GetFreeSurface());

  //- Automatic internal damping option (if negative BA provided,as damping ratio)
  if(temp_c<0) {
    const double zEta=-temp_c; // desired damping ratio
    C=zEta*UnstrLen/N*sqrt(E*Rho);   // Rho=w/A
    if(misc::wordy>1) Log->Printf("Line %u damping set to %.f Ns.\n",Idl,C); //wordy is in Misc.h
  }

  for(unsigned i=0; i<N; i++) {
    L[i]=UnstrLen/double(N);  // distribute line length evenly over segments
    V[i]=L[i]*0.25*PI*D*D;
  }
}

//==============================================================================
/// Initialize wave parameters for no waves situation
//==============================================================================
void ILine::SetupWaves(double* Ucurrent_in,float dt_in) {
  Ucurrent=Ucurrent_in;
  WaveDT=dt_in; // new variable for wave time step (should be same as WaveDT I think...)  

  if(Env->GetWaveKin()>0) Log->PrintWarning("Dummy function called when waves are supposed to be enabled!");

  if (misc::wordy>1) {
    Log->Printf("Setting up wave variables for ILine %u",Idl);
    Log->Printf("Nt= %d, and WaveDT= %.f, Depth= %.f\n",Nt,WaveDT,Depth);
  }

  Nw=0; // use no components since not worrying about wave kinematics
  Nt=2; // this is a new variable containing the Number of wave time steps to be calculated

  WGNC_Fact=1.0;
  S2Sd_Fact=1.0;
  // resize the new time series vectors
  zetaTS.resize(N+1,vector<double>(Nt,0.));
  FTS.resize(N+1,vector<double>(Nt,0.));
  UTS.resize(N+1,vector<vector< double>>(Nt,vector<double>(3,0.)));
  UdTS.resize(N+1,vector<vector< double>>(Nt,vector<double>(3,0.)));
  tTS.resize(Nt,0.);
}

//==============================================================================
/// Get ICs for line using quasi-static approach
//==============================================================================
void ILine::Initialize(double* X) {
  std::string function="Initialize";

  // set end node positions and velocities from connect objects
  //TODO: check these functions
  AnchConnect->GetConnectState(R[0],Rd[0]);
  FairConnect->GetConnectState(R[N],Rd[N]);

  //  cout << "IT: " << Zmin << ";  X= "<< R[0].x<< " ; Y= "<< R[0].y << " ; Z= "<< R[0].z << "-I: " <<0 <<endl ;
  if(-(Depth)>R[0].z) { 
    string tex=fun::PrintStr("Error: water depth is shallower than the anchor height (line_id=%d)",Idl);
    Run_Exceptioon(tex); 
  }

  // try to calculate initial line profile using catenary routine (from FAST v.7)
  // note: much of this function is adapted from the FAST source code
  //*
  // input variables for the Catenary function
  const double XF=sqrt(std::pow(std::abs((R[N].x-R[0].x)),2.0)+std::pow(std::abs((R[N].y-R[0].y)),2.0)); // quasi-static body line coordinate system (vertical plane with corners at anchor and fairlead) 
  const double ZF=std::abs(R[N].z-R[0].z);
  const double W=((Rho-Env->GetRho_w())*A)*(Env->GetG());
  const double CB=0.;
  const double Tol=0.00001;

  vector<double> snodes(N+1,0.0);             // locations of line nodes along line length-evenly distributed here 
  for(unsigned i=1; i<=N; i++) { snodes[i]=snodes[i-1]+L[i-1]; }
  snodes[N]=UnstrLen;                 // double check to ensure the last node does not surpass the line length

  // output variables
  double HF,VF,HA,VA,COSPhi,SINPhi;
  vector<double> Xl(N+1,0.0); // x location of line nodes
  vector<double> Zl(N+1,0.0);
  vector<double> Te(N+1,0.0);

  if(XF==0.0) {  // if the current body line is exactly vertical; thus, the solution below is ill-conditioned because the orientation is undefined; so set it such that the tensions and nodal positions are only vertical
    COSPhi=0.0;   SINPhi=0.0;
  }
  else {   // The current body line must not be vertical; use simple trigonometry
    COSPhi=(R[N].x-R[0].x)/XF;
    SINPhi=(R[N].y-R[0].y)/XF;
  }

  const int success=Qslines->Catenary(XF,ZF,UnstrLen,Ea,W,CB,Tol,&HF,&VF,&HA,&VA,N,snodes,Xl,Zl,Te,Idl);

  if(success >= 0) {  // assign the resulting line positions to the model
    for(unsigned i=1; i<N; i++) {
      R[i].x=R[0].x+Xl[i]*COSPhi;
      R[i].y=R[0].y+Xl[i]*SINPhi;
      R[i].z=R[0].z+Zl[i];
    }
  }
  else {  // otherwise just stretch the nodes between the endpoints linearly and hope for the best
    if(misc::wordy>0)Log->Printf("Catenary IC gen failed for ILine %u, so using linear node spacing.",Idl);
    for(unsigned i=1; i<N; i++) {
      R[i].x=R[0].x+(R[N].x-R[0].x)*(float(i)/float(N));
      R[i].y=R[0].y+(R[N].y-R[0].y)*(float(i)/float(N));
      R[i].z=R[0].z+(R[N].z-R[0].z)*(float(i)/float(N));
    }
  }
  // also assign the resulting internal node positions to the integrator initial state vector! (velocities leave at 0)
  for(unsigned i=1; i<N; i++) {
    const unsigned vidx=3*i-3;
    const unsigned pidx=3*N-3+vidx;
    X[pidx+0]=R[i].x; X[pidx+1]=R[i].y; X[pidx+2]=R[i].z; // positions
    X[vidx+0]=0.0;    X[vidx+0]=0.0;    X[vidx+0]=0.0;       // velocities=0
  }
  // now we need to return to the integrator for the dynamic relaxation stuff
  return;
}

//==============================================================================
/// Store the deoth for this line
//==============================================================================
void ILine::ConfigDepth(const IEnvironment* env,const IBody* b) {
  if (Depth==DBL_MAX) {
    if(env)Depth=env->GetWaterDepth();
    if(b)  Depth=b->GetDepth()!=DBL_MAX?b->GetDepth():Depth;
  }
  if(Depth==DBL_MAX)Run_Exceptioon(fun::PrintStr("The depth for the line %d is missing.", Idl));
}


//==============================================================================
/// Smart (selective) function to get tension at any node including fairlead or anchor 
/// (accounting for weight in these latter cases) 
//==============================================================================
double ILine::GetNodeTen(unsigned i)const {
  double NodeTen=0.0;
  if (i==0 || i==N) NodeTen=sqrt(Fnet[i].x*Fnet[i].x+Fnet[i].y*Fnet[i].y+(Fnet[i].z+M[i].a11*(-Env->GetG()))*(Fnet[i].z+M[i].a11*(-Env->GetG())));
  else {
    tdouble3 Tmag_squared=(T[i]+T[i-1])*(T[i]+T[i-1])*0.25;  // take average of tension in adjacent segments 
    NodeTen=sqrt(Tmag_squared.x+Tmag_squared.y+Tmag_squared.z);    //     previously used: NodeTen=0.5*(Tmag[i-1]+Tmag[i]); // should add damping in here too <<<<<<<<<<<<<
  }
  return NodeTen;
}

//==============================================================================
/// Function to return the position of any node along the line 
//==============================================================================
tdouble3 ILine::GetNodePos(unsigned idnode)const {
  if((idnode<0) && (idnode>N)) Run_Exceptioon(fun::PrintStr("The node %d is out of range for the line %d.",idnode,Idl));
  return TDouble3(R[idnode].x,R[idnode].y,R[idnode].z); 
}

//==============================================================================
/// FASTv7 style line tension outputs
//==============================================================================
void ILine::GetFASTtens(float* FairHTen,float* FairVTen,float* AnchHTen,float* AnchVTen){
  *FairHTen=float(sqrt(Fnet[N].x*Fnet[N].x+Fnet[N].y*Fnet[N].y));
  *FairVTen=float((Fnet[N].z+M[N].a11*(-Env->GetG())));
  *AnchHTen=float(sqrt(Fnet[0].x*Fnet[0].x+Fnet[0].y*Fnet[0].y));
  *AnchVTen=float((Fnet[0].z+M[0].a11*(-Env->GetG())));
}

//==============================================================================
/// Writes Forces and Mass into the vectors, which are passed how parameters
//==============================================================================
void ILine::GetAnchStuff(tdouble3& fnet,tmatrix3d& mass) {
  fnet=Fnet[0];
  mass=M[0];
}

//==============================================================================
/// Writes Forces and Mass into the vectors, which are passed how parameters
//==============================================================================
void ILine::GetFairStuff(tdouble3& fnet,tmatrix3d& mass) {
  fnet=Fnet[N];
  mass=M[N];
}

//==============================================================================
/// Returns the position (X,Y,Z) of one IConnection of this ILine
//==============================================================================
tdouble3 ILine::GetPositionOutput(const unsigned n)const {
  return R[n];
}

//==============================================================================
/// Returns the velocity (X,Y,Z)
//==============================================================================
tdouble3 ILine::GetVelocityOutput(const unsigned n)const {
  return Rd[n];
}

//==============================================================================
/// Returns the force (X,Y,Z)
//==============================================================================
tdouble3 ILine::GetForceOutput(const unsigned n)const {
  return TDouble3(Fnet[n].x,Fnet[n].y,Fnet[n].z);
}

//==============================================================================
/// Returns the Tension
//==============================================================================
double ILine::GetTensionOutput(const unsigned n)const {
  return GetNodeTen(n);
}

//==============================================================================
/// Function for boosting drag coefficients during IC generation  
//==============================================================================
void ILine::ScaleDrag(double scaler) {
  Cdn=Cdn*scaler;
  Cdt=Cdt*scaler;
}

//==============================================================================
// This is the big function that updates the States
//==============================================================================
void ILine::DoRHS(const double* X,double* Xd,const double time,const double dt) {
  std::string function="DoRHS";

  TimeStep=time;

  AnchConnect->GetConnectState(R[0],Rd[0]);
  FairConnect->GetConnectState(R[N],Rd[N]);

  //-Set interior node positions and velocities
  for(unsigned i=1; i<N; i++) {
    const unsigned vidx=3*i-3;
    const unsigned pidx=3*N-3+vidx;

    R[i] =TDouble3(X[pidx+0],X[pidx+1],X[pidx+2]); // get positions
    Rd[i]=TDouble3(X[vidx+0],X[vidx+1],X[vidx+2]); // get velocities
  }
  //-Calculate current (Stretched) segment lengths
  for(unsigned i=0; i<N; i++) {
    const tdouble3 lstr=(R[i+1]-R[i])*(R[i+1]-R[i]);
    const double lstr_squared=lstr.x+lstr.y+lstr.z;
    Lstr[i]=sqrt(lstr_squared);   // stretched segment length
    const tdouble3 ldstr=(R[i+1]-R[i])*(Rd[i+1]-Rd[i]);
    const double ldstr_top=ldstr.x+ldstr.y+ldstr.z;
    if(fun::IsNAN(ldstr_top)) Run_Exceptioon(fun::PrintStr("A NaN value detected in the strain rate of segment %d in the line %d (t=%gs).",i,Idl,TimeStep)); 
    if(Lstr[i]==double(0))Run_Exceptioon(fun::PrintStr("The stretched length of segment %d in the line %d cannot be 0 (t=%gs).",i,Idl,TimeStep)); 
    Ldstr[i]=ldstr_top/Lstr[i];   // strain rate of segment
    V[i]=PI/4.* (D*D*L[i]);    // volume attributed to segment
  }
  //-Calculate unit tangent vectors (Q) for each node (including ends)  note: I think these are pointing toward 0 rather than N!
  for(unsigned i=0; i<=N; i++) {
    if     (i==0) misc::unitvector(Q[i],R[i+1],R[i  ]); // compute unit vector Q
    else if(i==N) misc::unitvector(Q[i],R[i  ],R[i-1]); // compute unit vector Q
    else          misc::unitvector(Q[i],R[i+1],R[i-1]); // compute unit vector Q ... using adjacent two nodes!
  }

  //============================================================================================
  // --------------------------------- apply wave kinematics ------------------------------------
  std::fill(U.begin()   ,U.end() ,TDouble3(0));
  std::fill(Ud.begin()  ,Ud.end(),TDouble3(0));
  std::fill(F.begin()   ,F.end() ,0);

  // loop through the segments
  for(unsigned i=0; i<N; i++) {
    double eta=0.5*(Zeta[i]+Zeta[i+1]);
    F[i]=CalculateVOFSeg(R[i].z,R[i+1].z,eta);
  }

  CalculateMassMatrix(); //Calls calculate Mass matrix

  // ============  CALCULATE FORCES ON EACH NODE ===============================
  // loop through the segments
  for(unsigned i=0; i<N; i++) {

    // line tension
    if(Lstr[i]/L[i]>1.0) {  T[i]=(R[i+1]-R[i])*Ea*(1./L[i]-1./Lstr[i]); }
    else                 {  T[i]=TDouble3(0); }// cable can't "push"
    // line internal damping force
    Td[i]=(R[i+1]-R[i])*C*A*(Ldstr[i]/L[i])/Lstr[i]; 
  }

  const double rhow=Env->GetRho_w();
  const double g   =Env->GetG();
  const double kb  =Env->GetKb();
  const double cb  =Env->GetCb();
  const double d2  =D*D;
  // loop through the nodes
  for(unsigned i=0; i<=N; i++) {
  // submerged weight (including buoyancy)
    if     (i==0) W[i].z=PI/8.*(d2*L[i]  *(Rho-F[i]*rhow))  *(-g);
    else if(i==N) W[i].z=PI/8.*(d2*L[i-1]*(Rho-F[i-1]*rhow))*(-g);// missing the "W[i][2] =" previously!
    else          W[i].z=PI/8.*(d2*L[i]  *(Rho-F[i]*rhow)+d2*L[i-1]*(Rho-F[i-1]*rhow))*(-g);
   
    vi=U[i]-Rd[i]; // relative flow velocity over node
    vq=Q[i]*fgeo::ProductScalar(vi,Q[i]);
    vp=vi-vq;          // transverse relative flow component

    // flow velocity calculations       
    const tdouble3 vq_squared=vq*vq;
    const tdouble3 vp_squared=vp*vp;

    const double vp_mag=sqrt(vp_squared.x+vp_squared.y+vp_squared.z);
    const double vq_mag=sqrt(vq_squared.x+vq_squared.y+vq_squared.z);

    // transverse drag    
    if     (i==0) Dp[i]=vp*0.5*rhow*Cdn*(F[i]  *D*L[i])                /2.*vp_mag;
    else if(i==N) Dp[i]=vp*0.5*rhow*Cdn*(F[i-1]*D*L[i-1])              /2.*vp_mag;
    else          Dp[i]=vp*0.5*rhow*Cdn*(F[i]  *D*L[i]+F[i-1]*D*L[i-1])/2.*vp_mag;

    // tangential drag    
    if     (i==0) Dq[i]=vq*0.5*rhow*Cdt*PI*(F[i]  *D*L[i])                /2.*vq_mag;
    else if(i==N) Dq[i]=vq*0.5*rhow*Cdt*PI*(F[i-1]*D*L[i-1])              /2.*vq_mag;
    else          Dq[i]=vq*0.5*rhow*Cdt*PI*(F[i]  *D*L[i]+F[i-1]*D*L[i-1])/2.*vq_mag;
    // acceleration calculations          
    aq=Q[i]*fgeo::ProductScalar(Ud[i],Q[i]); // tangential component of fluid acceleration
    ap=Ud[i]-aq; // normal component of fluid acceleration

    // transverse Froude-Krylov force
    if     (i==0) Ap[i]=ap*rhow*(1.+Can)*0.5*(V[i]);
    else if(i==N) Ap[i]=ap*rhow*(1.+Can)*0.5*(V[i-1]);
    else          Ap[i]=ap*rhow*(1.+Can)*0.5*(V[i]+V[i-1]);

    // tangential Froude-Krylov force          
    if     (i==0) Aq[i]=ap*rhow*(1.+Cat)*0.5*(V[i]);
    else if(i==N) Aq[i]=ap*rhow*(1.+Cat)*0.5*(V[i-1]);
    else          Aq[i]=ap*rhow*(1.+Cat)*0.5*(V[i]+V[i-1]);

    // bottom contact (stiffness and damping, vertical-only for now)-updated for general case of potentially anchor or fairlead end in contact
    if(R[i].z<Zmin) {
      if     (i==0) B[i].z=((Zmin-R[i].z)*kb-Rd[i].z*cb)*0.5*(L[i]);//B[i].z=((Zmin-R[i].z)*kb-Rd[i].z*cb)*0.5*(D*L[i-1]); //-ERROR if i==0 => L[i-1] -> segmentation fault
      else if(i==N) B[i].z=((Zmin-R[i].z)*kb-Rd[i].z*cb)*0.5*(D*L[i]);
      else          B[i].z=((Zmin-R[i].z)*kb-Rd[i].z*cb)*0.5*(D*L[i]+D*L[i-1]);

      // new rough-draft addition of seabed friction
      //double FrictionCoefficient=0.5;                  // just using one coefficient to start with          
      const double FrictionMax=abs(B[i].z)*(Env->GetFricCoeff()); // dynamic friction force saturation level based on bottom contact force
      // saturated damping approach to applying friction, for now
      double BottomVel=sqrt(Rd[i].x*Rd[i].x+Rd[i].y*Rd[i].y); // velocity of node along sea bed

      double FrictionForce=BottomVel*(Env->GetFricCoeff())*(Env->GetFricDamp()); // some arbitrary damping scaling thing at end
      if(FrictionForce>Env->GetStatDynFricScale()*FrictionMax) { FrictionForce=FrictionMax; }    // saturate (quickly) to static/dynamic friction force level 
      // apply force in correct directions -- opposing direction of motion
      // could add ifs in here to handle end nodes
      if(BottomVel==double(0)) BottomVel=DBL_MIN; 
      if(std::isinf(BottomVel)) Run_Exceptioon(fun::PrintStr("The velocity of node %d in the line %d along sea bed cannot be infinite (t=%gs). Set the dtM to \'auto\' or reduce it.",i,Idl,TimeStep)); 
      
      B[i].x=-FrictionForce*Rd[i].x/BottomVel;
      B[i].y=-FrictionForce*Rd[i].y/BottomVel;
    }
    else  B[i]=TDouble3(0);
    
    // total forces
    if     (i==0) Fnet[i]= T[i  ]        +Td[i  ]      +W[i]+(Dp[i]+Dq[i]+Ap[i]+Aq[i])+B[i];
    else if(i==N) Fnet[i]=(T[i-1]*(-1))  -Td[i-1]      +W[i]+(Dp[i]+Dq[i]+Ap[i]+Aq[i])+B[i];
    else          Fnet[i]= T[i  ]-T [i-1]+Td[i]-Td[i-1]+W[i]+(Dp[i]+Dq[i]+Ap[i]+Aq[i])+B[i];
    
    if(Fnet[i]==TDouble3(0)) Fnet[i]=TDouble3(-DBL_MIN);
    if(fun::IsNAN(Fnet[i].x) || fun::IsNAN(Fnet[i].y) || fun::IsNAN(Fnet[i].z)) Run_Exceptioon(fun::PrintStr("A NaN value detected in the tension of the node %d in the line %d (t=%gs).",i,Idl,TimeStep)); 
  }

  // loop through internal nodes and update their States
  for(unsigned i=1; i<N; i++) {
    const unsigned vidx=3*i-3;
    const unsigned offset=3*N-3;
    fmdp::UpdateStates(offset,X+vidx,Xd+vidx,S[i],Fnet[i]);

    if(Debug)VisuProperties(i);
  }
  //-Checks the tension in the connections
  if(BkTension&&!Disabled)CheckTension(); 
}

//============================================================================================
//Calculate the VOF of each segment
//============================================================================================
double ILine::CalculateVOFSeg(double n1z,double n2z,double hswl)const{
  const double dn1z=n1z-hswl;
  const double dn2z=n2z-hswl;
  double f=0;
  if      (n1z<=0 && n2z<0)f=1; //Sumerged
  else if (n1z >0 && n2z>0)f=0; //Above fluid
  else f=0.5; //Partially submerged
  return f;
}

//============================================================================================
/// Checks if some IConnection reached the maximum allowed value of tension
//============================================================================================
void ILine::CheckTension(){
  const double tenF=GetNodeTen(N);
  const double tenA=GetNodeTen(0);
  if(std::abs(tenF)>=BkTension||std::abs(tenA)>=BkTension) BreakLine();
}

//============================================================================================
/// Disables the ILine and its Connections when the maximum value of tension is reached. 
/// Throws a warning.
//============================================================================================
void ILine::BreakLine(){
  FairConnect->SetDisabled(true); FairConnect->ResetTension();
  AnchConnect->SetDisabled(true); AnchConnect->ResetTension();
  for(unsigned i=0; i<N; i++) Fnet[i]=TDouble3(0);
  if(!Disabled)Log->PrintfWarning("The line %d reached the maximum value of tension (%gN). It will be disabled.",Idl,BkTension);
  Disabled=true;
}

//============================================================================================
/// Store in M the mass of nodes
//============================================================================================
void ILine::CalculateMassMatrix() {
  // calculate mass matrix 
  const double d2=D*D;
  for(unsigned i=0; i<=N; i++) {
    double m_i; // node mass
    double v_i; // node submerged volume 

    if     (i==0) { m_i=PI/8.* d2*L[0]  *Rho;                v_i=0.5*F[i]   *V[i];              }
    else if(i==N) { m_i=PI/8.* d2*L[N-1]*Rho;                v_i=0.5*F[i-1] *V[i-1];            }
    else          { m_i=PI/8.*(d2       *Rho*(L[i]+L[i-1])); v_i=0.5*(F[i-1]*V[i-1]+F[i]*V[i]); }

    // make node mass matrix
    const tmatrix3d imat=TMatrix3d();//<identity 3x3
    const tmatrix3d qmat=fmath::ProductMatrix3x3(Q[i],Q[i]);
    M[i]=imat*m_i+((imat-qmat)*Can+qmat*Cat)*Env->GetRho_w()*v_i;
    S[i]=fmath::InverseMatrix3x3(M[i]);// invert node mass matrix (written to S[i][3x3])  
  }
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void ILine::VisuConfig()const{	
  Log->Printf("  Line %d:",Idl);	
  Log->Printf("    UnstrLength...: %g [m]",UnstrLen);
  Log->Printf("    Segments......: %d",N);
  Log->Printf("    Diameter......: %g [m]",D);
  Log->Printf("    E (Stiffness).: %g [N/m^2]",E);
  Log->Printf("    Area..........: %g [m^2]",A);
  Log->Printf("    EA (E*Area)...: %g [N]",Ea);
  Log->Printf("    ba (Damping)..: %g [Ns]",C);
  Log->Printf("    MassDenInAir..: %g [kg/m]",MassDenInAir);
  Log->Printf("    Density.......: %g [kg/m^3]",Rho);
  Log->Printf("    Break Tension.: %g [N]",BkTension);
  Log->Printf("    Can...........: %g",Can);
  Log->Printf("    Cat...........: %g",Cat);
  Log->Printf("    Cdn...........: %g",Cdn);
  Log->Printf("    Cdt...........: %g",Cdt);
  Log->Print ("    Connections...: ");
  Log->Printf("      Type:%s - Point(%s)",AnchConnect->GetTypeStr().c_str(),fun::Double3gStr(AnchConnect->GetPositions()).c_str());
  Log->Printf("      Type:%s - Point(%s)",FairConnect->GetTypeStr().c_str(),fun::Double3gStr(FairConnect->GetPositions()).c_str());
}

//==============================================================================
/// Shows initial tension in the nodes
//==============================================================================
void ILine::VisuIniTen()const{
  Log->Printf("  Line %d:",Idl);
  if(AnchConnect->GetType()!=TpConnection::connect)Log->Printf("    Point(%s): %g [N]",fun::Double3gStr(AnchConnect->GetPositions()).c_str(),GetTensionOutput(AnchConnect->GetNode()));
  if(FairConnect->GetType()!=TpConnection::connect)Log->Printf("    Point(%s): %g [N]",fun::Double3gStr(FairConnect->GetPositions()).c_str(),GetTensionOutput(FairConnect->GetNode()));
}

//============================================================================================
/// Write output file for line  (accepts time parameter since retained time value (t) 
/// will be behind by one line time step
//============================================================================================
void ILine::WriteOutput(double timestep,double dt,double start,double finish,double next) {
  // run through output flags
  // if channel is flagged for output, write to file.
  // Flags changed to just be one character (case sensitive) per output flag.  To match FASTv8 version.
  // create output file for writing output (and write channel header and units lines) if applicable
  std::string function="WriteOutput";
  const double t=timestep;
  bool savedata=false;
  if(start<=timestep && timestep<=finish) savedata=next<=timestep;

  if(savedata) { // check it's not null.  Null signals no individual line output files
    std::string chaninp=Channels;
    if((chaninp.length()>=1) && chaninp[0]!='-') {
      for(unsigned c=0; c<(unsigned)chaninp.length(); c++) {
        const char ch=chaninp[c]; //takes a character
        //-Save Positions
        if(ch=='p'){
          const string file=DirOut+fun::PrintStr("Line_%u_pos.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_pos.csv", "Saves position of the nodes.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<=N; i++) scsv << fun::PrintStr("pos_%d.x [m];pos_%d.y [m];pos_%d.z [m]",i,i,i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<=N; i++) scsv << R[i];
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
        //-Save Velocities
        if(ch=='v'){
          const string file=DirOut+fun::PrintStr("Line_%u_vel.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_vel.csv", "Saves velocities of the nodes.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<=N; i++) scsv << fun::PrintStr("vel_%d.x [m/s];vel_%d.y [m/s];vel_%d.z [m/s]",i,i,i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<=N; i++) scsv << Rd[i];
          scsv << jcsv::Endl();
          scsv.SaveData();
        }

        //-Wave velocities
        if(ch=='U'){
          const string file=DirOut+fun::PrintStr("Line_%u_wavevel.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_wavevel.csv", "Saves wave velocities of the nodes.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<=N; i++) scsv << fun::PrintStr("wvel_%d.x [m/s];wvel_%d.y [m/s];wvel_%d.z [m/s]",i,i,i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<=N; i++) scsv << U[i];
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
        //-Hydro force
        if(ch=='D'){
          const string file=DirOut+fun::PrintStr("Line_%u_hydroforce.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_hydroforce.csv", "Saves hydro forces of the nodes.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<=N; i++) scsv << fun::PrintStr("hforce_%d.x [N];hforce_%d.y [N];hforce_%d.z [N]",i,i,i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<=N; i++) scsv << Dp[i]+Dq[i]+Ap[i]+Aq[i];
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
        //-Segment internal damping force
        if(ch=='c'){
          const string file=DirOut+fun::PrintStr("Line_%u_damping.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_damping.csv", "Saves internal damping of the segments.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<N; i++) scsv << fun::PrintStr("c_%d.x [N];c_%d.y [N];c_%d.z [N]",i,i,i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<N; i++) scsv << Dp[i]+Dq[i]+Ap[i]+Aq[i];
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
        //-Segment tensions
        if(ch=='t'){
          const string file=DirOut+fun::PrintStr("Line_%u_tension.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_tension.csv", "Saves tension of the segments.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<N; i++) scsv << fun::PrintStr("tension_%d [N]",i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<N; i++){
            const tdouble3 Tmag_squared=T[i]*T[i];
            scsv << sqrt(Tmag_squared.x+Tmag_squared.y+Tmag_squared.z);
          }
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
        //-Segment strains
        if(ch=='s'){
          const string file=DirOut+fun::PrintStr("Line_%u_strain.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_strain.csv", "Saves strain of the segments.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<N; i++) scsv << fun::PrintStr("St_%d [-/s]",i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<N; i++) scsv << Lstr[i]/L[i]-1.0;
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
        //-Segment strain rates
        if(ch=='d'){
          const string file=DirOut+fun::PrintStr("Line_%u_strainrate.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_strainrate.csv", "Saves strain rates of the segments.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<N; i++) scsv << fun::PrintStr("dSt_%d [-/s]",i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<N; i++) scsv << Ldstr[i]/L[i];
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
        //-Bottom contact force
        if(ch=='b'){
          const string file=DirOut+fun::PrintStr("Line_%u_bottomcontact.csv",Idl);
          jcsv::JSaveCsv2 scsv(file,true,false);
          if(!scsv.GetAppendMode()){
            Log->AddFileInfo("Line_XX_bottomcontact.csv", "Saves bottom contact force of the nodes.");
            //-Saves head.
            scsv.SetHead();
            scsv << "time [s]";
            for(unsigned i=0; i<=N; i++) scsv << fun::PrintStr("B_%d.x [N];B_%d.y [N];B_%d.z [N]",i,i,i);
            scsv << jcsv::Endl();
          }
          //-Saves data.
          scsv.SetData();
          scsv << t ;
          for(unsigned i=0; i<=N; i++)  scsv << B[i];
          scsv << jcsv::Endl();
          scsv.SaveData();
        }
      }
    }
  }
}

//============================================================================================
/// Visualises the properties of the line
//============================================================================================
void ILine::VisuProperties(unsigned ini,unsigned end)const{
  const unsigned endn=end==UINT_MAX?ini+1:end+1;
  if(ini>=0 && ini<endn && endn>0&& endn<=N+1){
    printf("Line=%u:\n",Idl);
    printf("  Lstr:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g)\n"      ,i,Lstr[i]);
    printf("  Ldstr:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g)\n"      ,i,Ldstr[i]);
    printf("  L...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g)\n"      ,i,   L[i]);
    printf("  F...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g)\n"      ,i,   F[i]);
    printf("  R...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,   R[i].x,   R[i].y,   R[i].z);
    printf("  Rd..:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,  Rd[i].x,  Rd[i].y,  Rd[i].z);
    printf("  Q...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,   Q[i].x,   Q[i].y,   Q[i].z);
    printf("  W...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,   W[i].x,   W[i].y,   W[i].z);
    printf("  U...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,   U[i].x,   U[i].y,   U[i].z);
    printf("  Ud..:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,  Ud[i].x,  Ud[i].y,  Ud[i].z);
    printf("  Fnet:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,Fnet[i].x,Fnet[i].y,Fnet[i].z);
    printf("  B...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,   B[i].x,   B[i].y,   B[i].z);
    printf("  T...:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,   T[i].x,   T[i].y,   T[i].z);
    printf("  Td..:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,  Td[i].x,  Td[i].y,  Td[i].z);
    printf("  Ap..:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,  Ap[i].x,  Ap[i].y,  Ap[i].z);
    printf("  Aq..:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,  Aq[i].x,  Aq[i].y,  Aq[i].z);
    printf("  Dp..:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,  Dp[i].x,  Dp[i].y,  Dp[i].z);
    printf("  Dq..:\n");for(unsigned i=ini; i<endn; i++) printf("    [%u]:(%g,%g,%g)\n",i,  Dq[i].x,  Dq[i].y,  Dq[i].z);
  }
}


// new function to draw instantaneous line positions in openGL context
#ifdef USEGL
void ILine::drawGL(void) {
  double maxTen=0.0;
  double normTen;
  double rgb[3];
  for(unsigned i=0; i<=N; i++) {
    double newTen=GetNodeTen(i);
    if(newTen>maxTen) { maxTen=newTen; }
  }

  glColor3f(0.5,0.5,1.0);
  glBegin(GL_LINE_STRIP);
  for(unsigned i=0; i<=N; i++) {
    glVertex3d(R[i].x,R[i].y,R[i].z);
    if(i<N) {
      normTen=GetNodeTen(i)/maxTen;
      ColorMap(normTen,rgb);
      glColor3d(rgb[0],rgb[1],rgb[2]);
    }
  }
  glEnd();
}

void ILine::drawGL2(void) {
  double maxTen=0.0;
  double normTen;
  double rgb[3];
  for(unsigned i=0; i<=N; i++) {
    double newTen=GetNodeTen(i);
    if(newTen>maxTen) { maxTen=newTen; }
  }
  // line
  for(unsigned i=0; i<N; i++) {
    normTen=0.2+0.8*pow(GetNodeTen(i)/maxTen,4.0);
    ColorMap(normTen,rgb);
    glColor3d(rgb[0],rgb[1],rgb[2]);
    Cylinder(R[i].x,R[i].y,R[i].z,R[i+1][0],R[i+1][1],R[i+1][2],27,0.5);
  }
  // velocity vectors
  for(unsigned i=0; i<=N; i++) {
    glColor3d(0.0,0.2,0.8);
    double vscal=5.0;
    Arrow(R[i].x,R[i].y,R[i].z,vscal*Rd[i].x,vscal*Rd[i].y,vscal*Rd[i].z,0.1,0.7);
  }
}
#endif