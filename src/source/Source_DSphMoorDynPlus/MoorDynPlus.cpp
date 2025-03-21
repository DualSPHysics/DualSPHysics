/*===================================================================================
<MOORDYNPLUS> Copyright (c) 2025 by Ivan Martinez-Estevez

Ivan Martinez Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
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

/// \file MoorDynPlus.cpp \brief Implements the class \ref MoorDynPlus.

#include "MoorDynPlus.h"
#include "JLog2.h"
#include "FunMoorDynPlus.h"
#include "IApi.h"
#ifdef DDISABLE_DSPH
#include "JSpVtkShape.h"
#endif // DDISABLE_DSPH
#include "FunGeo3d.h"

#include <string>
#include <time.h>

//==============================================================================
/// Constructor.
//==============================================================================
MoorDynPlus::MoorDynPlus() {
  FairArrays=false; ItArrays=false;
  NStateVariables=6;
  FileIN="moordynplus.xml";
  RootNodeXml="moordynplus";
  Env=NULL; OutP=NULL;
  LineDef=NULL; Log=NULL;
  Coupled=true; DsphLog=false;
  Part=0;
#ifdef WIN32
  Dirout="\\MoorDynPlus_out"; ///<Default output directory (Is replaced when an external software calls to LinesInit() )
#else
  Dirout="/MoorDynPlus_out"; ///<Default output directory (Is replaced when an external software calls to LinesInit() )
#endif // WIN32
  Reset();
  ClassName="MoorDynPlus";
}

//==============================================================================
/// Destructor.
//==============================================================================
MoorDynPlus::~MoorDynPlus() {
  fmdp::FreeVector(Lines);
  fmdp::FreeVector(Connections); 
  fmdp::FreeVector(Bodies); 
  fmdp::FreeArray2D(Ffair,FairCount);
  fmdp::FreeArray2D(rFairi,FairCount);
  fmdp::FreeArray2D(rdFairi,FairCount);
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void MoorDynPlus::Reset() {
  delete Env;     Env=NULL;
  delete OutP;    OutP=NULL;
  delete LineDef; LineDef=NULL;
  if(!DsphLog) delete Log;
  if(FairArrays) {
    if(FairleadPos) fmdp::FreeArray3D(FairleadPos,BodyCount);  
    if(FairleadVel) fmdp::FreeArray3D(FairleadVel,BodyCount);  
    if(FairleadFor) fmdp::FreeArray3D(FairleadFor,BodyCount);  
  }
  Log=NULL;         FairleadPos=NULL; 
  FairleadVel=NULL; FairleadFor=NULL;
  References.clear();  OutLines.clear();
  XmlReaded=FairArrays=HeaderProg=UseWaves=false;
  BodyCount=LineCount=ConnectionCount=0;
  FairCount=AnchCount=ConnCount=0;
  TimeMax=DtOut=0;
  if(ItArrays){
    delete[] X0;    X0=NULL;
    delete[] F0;    F0=NULL;
    delete[] F1;    F1=NULL;
    //delete[] F2;  F2=NULL;
    //delete[] F3;  F3=NULL;
    delete[] Xt;    Xt=NULL;
  }
  ItArrays=false;
  Debug=false;
  NX=0;
  TransMat=TMatrix3d(0);
  #ifdef DISABLE_DSPH
  delete Log;  Log=NULL;
  #endif //DISABLE_DSPH
}

//==============================================================================
/// Creates the output directory for store LineX_BodyX.txt files.
//==============================================================================
void MoorDynPlus::CreateOutDirectory() {
  Dirout=fun::GetCurrentDir()+Dirout;
  fun::Mkdir(Dirout);

#ifdef WIN32
  Dirout=Dirout+"\\";
#else
  Dirout=Dirout+"/";
#endif
  //PRINT: look the output directory
  if(Debug) Log->PrintfDbg("Output directory: %s\n",Dirout.c_str());
}

//==============================================================================
/// Creates a JLog2 system if DualSPHysics doesn't send its JLog2 object
//==============================================================================
void MoorDynPlus::CreateLog(){
  Log=new JLog2;
  Log->Init(Dirout+"Run_MoorDynPlus.out");
  Log->AddFileInfo(Dirout+"Run_MoorDynPlus.out","Log file of the MoorDynPlus simulation.");
}

//==============================================================================
/// Allocates memory for arrays of integration
//==============================================================================
void MoorDynPlus::AllocateMemoryIntegration(const unsigned size) {
  std::string function="AllocateMemoryIntegration"; 

  if(size) {
    try {
      X0=new double[size];
      F0=new double[size];
      F1=new double[size];
      //F2=new double[size];  
      //F3=new double[size];                                        
      Xt=new double[size];
    }
    catch(std::bad_alloc const&) {
      Run_Exceptioon("Could not allocate the requested memory.");
    }
    ItArrays=true;
    std::memset(X0,0,sizeof(double)*size);
    std::memset(F0,0,sizeof(double)*size);
    std::memset(F1,0,sizeof(double)*size);
    //std::memset(F2,0,sizeof(double)*size);
    //std::memset(F3,0,sizeof(double)*size);
    std::memset(Xt,0,sizeof(double)*size);
  }
}

//==============================================================================
/// Allocates memory for persistent vectors
//==============================================================================
//void MoorDynPlus::AllocateMemoryPersistent(){
//FlinesS.resize(nStateVariables);
//rFairtS=GetPtrToPtr(FairCount,3);
//rFairRel=GetPtrToPtr(FairCount,3);
//
//
////rFairi.resize (FairCount);  // after applying platform DOF ICs,should eventually pass this rather than rFairt to ILine.setup()
////rdFairi.resize(FairCount);
//}

//==============================================================================
///Receives the log to store (in run.out) and print messages
//==============================================================================
void MoorDynPlus::LogInit(JLog2* log) {
  DsphLog=true;
  if(log)Log=log;
  else Run_Exceptioon("JLog2 object is not initialised.");
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void MoorDynPlus::LoadXml() {
  std::string function="LoadXml"; 
  XmlReaded=true;
  if(!fun::FileExists(FileIN))Run_ExceptioonFile("MoorDynPlus configuration was not found.",FileIN);
  JXml xml; xml.LoadFile(FileIN);
  const TiXmlNode* node=xml.GetNode(RootNodeXml,false);
  if(!node)Run_Exceptioon("Cannot find the element \'"+RootNodeXml+"\'.");
  ReadXml(&xml);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void MoorDynPlus::ReadXml(JXml* sxml) {
  std::string function="ReadXml";
  std::string place=RootNodeXml;
  TiXmlNode* root;
  if(sxml->GetNode(place,false)) root=(sxml->GetNode(RootNodeXml,false));
  else Run_Exceptioon("Cannot find the element \'"+RootNodeXml+"\'.");

  //TiXmlElement* rootele=(root->ToElement())->FirstChildElement(fun::PrintStr("%s",RootNodeXml).c_str());
  const TiXmlElement* rootele=root->ToElement();
  sxml->CheckElementNames(rootele,true,"solverOptions waves bodies connects lines savedata output");

  //-Get mode
  Debug=sxml->GetAttributeBool(rootele,"debug",true,false);

  //-Configuration of IEnvironment
  place=RootNodeXml+".solverOptions";
  if(sxml->GetNode(place,false)) {
    Env=new IEnvironment(Log);
    Env->LoadXml(sxml,RootNodeXml);
    //PRINT: look when reads solverOptions tag
    if(Debug) Log->PrintfDbg("(%s)-Environment\n",function.c_str());
  }else Run_Exceptioon("Section <solverOptions> is missing.");

  //-Configuration of Bodies
  place=RootNodeXml+".bodies";
  if(sxml->GetNode(place,false)) {
    const unsigned nBodiesTag=sxml->CountElements(root,"bodies");
    if(nBodiesTag>1)Run_Exceptioon("Error: There are more \"bodies\" section than one");
    TiXmlElement* eleBodies=(sxml->GetNode(place,false))->ToElement();
    if(Debug) Log->PrintfDbg("(%s)-Body-Num Bodies: %d\n",function.c_str(),BodyCount); 
    BodyCount=sxml->CountElements(sxml->GetNode(place,false),"body");  

    TiXmlElement* eleBody=eleBodies->FirstChildElement("body");
    for(unsigned b=0; b<BodyCount; b++){
      IBody* body=new IBody(Log);
      body->LoadXml(sxml,place,eleBody,unsigned(Bodies.size()));
      References.push_back(body->GetId());
      CheckReferences(body->GetRef()); //check if exist someone reference equals than this
      CheckDepth(Env,body); //Checks if the user inserts a depth for a body (if not, stores a IEnvironment detph)
      Bodies.push_back(body); //stores the IBody
      if(b+1<BodyCount)eleBody=sxml->GetNextElement(eleBody,"body");
    }
    ConfigBodies();
  }

  //-Configuration of connects
  place=RootNodeXml+".connects";
  std::vector<IConnection*> connects;
  if(sxml->GetNode(place,false)) {  
    //-Count "connects" blocks  
    const unsigned nConnectsTag=sxml->CountElements(root,"connects");
    if(nConnectsTag>1){Run_Exceptioon("Error: There are more \"connects\" section than one");}
    //-Count connect tags  
    TiXmlElement* eleConnects=(sxml->GetNode(place,false))->ToElement();
    const unsigned ConnectionCount=sxml->CountElements(sxml->GetNode(place,false),"connect");
    if(Debug) Log->PrintfDbg("(%s)-Num Connects: %d\n",function.c_str() ,ConnectionCount);
    TiXmlElement* eleConn=eleConnects->FirstChildElement("connect");
    for(unsigned c=0; c<ConnectionCount; c++){
      IConnection* conn=new IConnection(Log,TpConnection::connect);
      conn->LoadXml(sxml,place,eleConn,unsigned(connects.size()));
      CheckConnReferences(conn->GetRef(),connects); //check if exist someone reference equals than this
      connects.push_back(conn); //stores the connect
      if(c+1<ConnectionCount)eleConn=sxml->GetNextElement(eleConn,"connect");
    }
  }

  //-Configuration of lines
  place=RootNodeXml+".lines";
  if(!sxml->GetNode(place,false)) {Run_Exceptioon("Cannot find the element \'"+place+"\'.");  }
  else{
    const unsigned nLinesTag=sxml->CountElements(root,"lines");
    if(nLinesTag>1){Run_Exceptioon("Error: There are more \"lines\" section than one");}
    TiXmlElement* eleLines=(sxml->GetNode(place,false))->ToElement();
    //-Configuration of ILineBase
    std::string placeLinDef=place+".linedefault";
    LineDef=NULL;
    if(sxml->GetNode(placeLinDef,false)) {
      unsigned nLinesDefTag=sxml->CountElements(sxml->GetNode(place,false),"linedefault");
      if(nLinesDefTag>1){Run_Exceptioon("Error: There are more \"linedefault\" section than one");}
      LineDef=new ILineBase();
      LineDef->LoadXml(sxml,place);
    }

    //-Configuration of ILine
    std::string placeLine=place+".line";
    unsigned lcount=sxml->CountElements(sxml->GetNode(place,false),"line");
    if(Debug) Log->PrintfDbg("(%s)-Num Lines: %d\n",function.c_str(),lcount); 
    TiXmlElement* eleLine=eleLines->FirstChildElement("line");
    for(unsigned l=0; l<lcount; l++){
      ILine* line=new ILine(Log);
      line->LoadXml(sxml,placeLine,eleLine,(unsigned)Lines.size(),connects,LineDef);
      Lines.push_back(line); 
      if(line->GetUseOutput()){OutLines.push_back(line);} //Stores the lines which will store the data in csv file
      Connections.push_back(line->GetAnchConnect());
      Connections.push_back(line->GetFairConnect());
      if(l+1<lcount)eleLine=sxml->GetNextElement(eleLine,"line");
    }
    LineCount=(unsigned)Lines.size();
  }
  ConfigConnections(connects);
  CheckBodyRefs();
  //-Configuration of IOutput
  place=RootNodeXml+".savedata";
  std::string place2=RootNodeXml+".output"; //-For compatibility with old versions
  if(sxml->GetNode(place,false)||sxml->GetNode(place2,false)) {
    OutP=new IOutput(TimeMax,DtOut);
    OutP->LoadXml(sxml,RootNodeXml,OutLines);
  }else OutP=NULL;
}

//==================================================================================
// Checks if there are Fair connections with IBody reference which doesn't exist
//==================================================================================
void MoorDynPlus::CheckBodyRefs(){
  std::string function="CheckBodyRefs";  
  std::vector<IConnection*> fairs=GetConnections(TpConnection::vessel);
  for(unsigned c=0; c<fairs.size(); c++){
    const unsigned bRef=fairs[c]->GetRef(); //Gets the body reference for each fairlead 
    if(!GetBody(bRef)){ Run_Exceptioon("Error: Cannot find someone IBody with id "+fmdp::ToString<unsigned>(bRef)+"."); }
  }
}

//==================================================================================
/// Checks if all ids of floatings are matching with ids of Bodies
//==================================================================================
void MoorDynPlus::CheckFtIds(const unsigned ids[],const unsigned numFt)const {
  std::string function="CheckFtIds";
  for(unsigned i=0; i<numFt; i++) {
    const IBody* body=GetBody(ids[i]);
    if(!body) Run_Exceptioon("Error: Cannot find someone IBody with id "+fmdp::ToString<unsigned>(ids[i])+"."); 
  }
}

//==============================================================================
/// Searches if Bodies arrays contains some NULL 
/// ( When user doesn't write lines for any IBody )
/// If founds some NULL deletes it and updates number of Bodies
//==============================================================================
void MoorDynPlus::ConfigBodies() {
  //-Searches Null and deletes it
  unsigned numBodies=(unsigned)Bodies.size();
  for(unsigned m=0; m<numBodies; m++) {
    if(Bodies[m]==NULL) {
      for(unsigned mn=m; mn<numBodies-1; mn++) {Bodies[mn]=Bodies[mn+1];}
      BodyCount--;
    }
  }
  //-Sort by id
  numBodies=(unsigned)Bodies.size(); //Updates if was modify
  for(unsigned c=0;c<numBodies-1;c++) {
    for(unsigned c2=c+1;c2<numBodies;c2++) {
      if(Bodies[c]->GetRef()>Bodies[c2]->GetRef()) { //Sort by ftid and changes Number
        IBody* mo=Bodies[c];
        Bodies[c]=Bodies[c2];
        Bodies[c]->SetId(c2);
        Bodies[c2]=mo;
        Bodies[c2]->SetId(c);
      }
    }
  }
  BodyCount=(unsigned)Bodies.size();
}

//==============================================================================
/// Configures the Connections
//==============================================================================
void MoorDynPlus::ConfigConnections(std::vector<IConnection*> connects){
  std::string function="ConfigConnections";
  //-Configures connect Connections
  std::vector<IConnection*> conns_temp;
  for(unsigned c=0; c<connects.size(); c++){
    IConnection* conn=connects[c]; //Gets the current connect
    const unsigned ref=conn->GetRef();
    for(unsigned l=0; l<LineCount; l++){
      const ILine* line=Lines[l]; //Gets the current line
      if(line->GetUseConnect()){ 
        std::vector<IConnection*> connsLine=line->GetConnections(TpConnection::connect); //Gets connects
        for(unsigned c=0; c<line->GetNConns(); c++){
          if(connsLine[c]->GetRef()==ref){ 
            if(std::find(conns_temp.begin(),conns_temp.end(),connsLine[c])==conns_temp.end()) {
              conns_temp.push_back(conn); //-Attached does not contain line
              break; // Goes out of the loop
            }
          }
        }
      }
    }
  }
  //-Checks if all connects are used
  if(connects.size()!=connects.size()){Run_Exceptioon("Error: There are \"connects\" not used by the lines.");}
  //-Stores the number of connects used
  ConnCount=(unsigned)conns_temp.size();
  //-Configures anch and fair Connections
  std::vector<IConnection*> connections_temp; 
  for(unsigned l=0; l<LineCount; l++){
    const ILine* line=Lines[l]; //gets the current line
    std::vector<IConnection*> connections=line->GetConnections(); //Gets all connections
    for(unsigned c=0; c<line->GetNConnections(); c++){
      IConnection* connection=connections[c]; //Gets the current connect
      if(connection->GetType()!=TpConnection::connect)connections_temp.push_back(connection);//If the connection isn not a connect
    }
  }
  //-Stores the number of all connections
  ConnectionCount=(unsigned)connections_temp.size()+ConnCount;
  Connections.clear();
  //-Orders the connections
  //Stores connect Connections
  for(unsigned c=0; c<ConnCount; c++){ 
    Connections.push_back(conns_temp[c]);
    Connections[c]->SetId(c);
  }
  //Stores anch and fair Connections
  for(unsigned c=ConnCount; c<ConnectionCount; c++){
    Connections.push_back(connections_temp[c-ConnCount]);
    Connections[c]->SetId(c);
    if(Connections[c]->IsFixed()){AnchCount++;}
    else if(Connections[c]->IsVessel()){FairCount++;}
  }
  ConnectionCount=(unsigned)Connections.size();
}

//==============================================================================
/// Checks if exists equals references
//==============================================================================
void MoorDynPlus::CheckReferences(const unsigned ref) {
  std::string function="CheckReferences";
  for(unsigned m=0; m<Bodies.size(); m++){
    if(Bodies[m] && Bodies[m]->GetRef()==ref) Run_Exceptioon(fun::PrintStr("Body ref=%s is duplicated.",fmdp::ToString(ref).c_str()));
  }
}

//==============================================================================
/// Checks if exists equals references
//==============================================================================
void MoorDynPlus::CheckConnReferences(const unsigned ref,std::vector<IConnection*> conns) {
  std::string function="CheckConnReferences";
  for(unsigned m=0; m<(unsigned)conns.size(); m++){
    if(conns[m] && conns[m]->GetRef()==ref) Run_Exceptioon(fun::PrintStr("Connect ref=%s is duplicated.",fmdp::ToString(ref).c_str()));
  }
}

//==============================================================================
/// Checks depth value of body and changes it if doesn't initialized
//==============================================================================
void MoorDynPlus::CheckDepth(IEnvironment* env,IBody* body) {
  if(body->GetDepth()==DBL_MAX) {
    body->SetDepth(env->GetWaterDepth());
    //PRINT: look if depth is changed because this tag doesn't exists
    if(Debug) Log->PrintfDbg("\n**** The wather depth of IBody %d is %.2f( WaterDepth chosen in solverOptions ) \n",body->GetRef(),body->GetDepth());
  }
}

//==============================================================================
/// Return a pointer to array of IBody
//==============================================================================
IBody* MoorDynPlus::GetBody(const unsigned ref)const {
  IBody* body=NULL;
  for(unsigned m=0; m<BodyCount && !body; m++) if(Bodies[m]->GetRef()==ref) body=Bodies[m];
  return body;
}

//==============================================================================
/// Returns a pointer to IBody that belongs to the line sent
//==============================================================================
IBody* MoorDynPlus::GetBody(ILine* line)const {
  IBody* body=NULL;
  std::vector<IConnection*> fairs=line->GetConnections(TpConnection::vessel); //Gets Fairs
  for(unsigned c=0; c<fairs.size() && !body; c++){ //For each fair
    const unsigned ref=fairs[c]->GetRef(); //Gets ref of the current fairlead
    for(unsigned m=0; m<BodyCount; m++) if(Bodies[m]->GetRef()==ref) body=Bodies[m];
  }
  return body;
}

//==============================================================================
/// Return a pointer to array of lines of IBody 
//==============================================================================
std::vector<ILine*> MoorDynPlus::GetLines(IBody* body_in)const {
  std::vector<ILine*> lines;
  const unsigned nL=body_in->GetNLines(); // num lines

  const IBody* body=body_in; //Current body
  for(unsigned l=0; l<nL; l++) lines.push_back(body->GetLines()[l]);

  return lines;
}

//==============================================================================
/// Return a pointer to array of lines of IBody by id
//==============================================================================
std::vector<ILine*> MoorDynPlus::GetLines(const unsigned ref)const {
  std::vector<ILine*> lines;
  for(unsigned l=0; l<Lines.size(); l++) {
    ILine* line=Lines[l];
    if(line->GetUseFair()){
      std::vector<IConnection*> connsF=line->GetConnections(TpConnection::vessel);
      for(unsigned c=0; c<connsF.size(); c++){
        if(connsF[c]->GetRef()==ref){lines.push_back(line);break;}
      }
    }
  }
  return lines;
}

//==============================================================================
/// Return a pointer to array of lines without body
//==============================================================================
std::vector<ILine*> MoorDynPlus::GetLinesNoBody()const {
  std::vector<ILine*> lines;
  for(unsigned l=0; l<Lines.size(); l++) {
    ILine* line=Lines[l];
    if(!line->GetUseFair())lines.push_back(line);  
  }
  return lines;
}

//==============================================================================
/// Return a pointer to array of connections of a Type
//==============================================================================
std::vector<IConnection*> MoorDynPlus::GetConnections(const TpConnection type)const {
  std::vector<IConnection*> connections;//Stores the connections by type
  for(unsigned c=0; c<ConnectionCount; c++) {
    IConnection* conn=Connections[c];
    if(conn->GetType()==type){connections.push_back(conn);}
  }
  return connections;
}

//==============================================================================
/// Removes the broken line from the system
//==============================================================================
void MoorDynPlus::RemoveLine(ILine* line){
  if(LineCount){
    Lines.erase(std::remove(Lines.begin(),Lines.end(),line),Lines.end());
    LineCount--;
  }
}

//==============================================================================
/// Searches the connected lines with a connect IConnection
//==============================================================================
void MoorDynPlus::AddAttachedLines(){
  std::vector<IConnection*> conns=GetConnections(TpConnection::connect); //-Catches the connects connections
  for(unsigned cc=0; cc<ConnCount; cc++){
    IConnection* conn=conns[cc]; //For each connect connection
    const unsigned ref=conn->GetRef();
    for(unsigned l=0; l<LineCount; l++){ // Searches the connected lines to this connection
      ILine* line=Lines[l];
      //-If there are connects connections with the same reference,add the line to this connect
      if(line->GetNConns(ref)>0)conn->AddLineToConnect(line);
    }
  }
}

//==============================================================================
/// Calls to setups of all objects.
//==============================================================================
void MoorDynPlus::Setup() {
  std::string function="Setup";
  double* ucurr=fmdp::GetArray(3);

  //-Setup of Connections
  for(unsigned c=0;c<ConnectionCount;c++)if(Connections[c])Connections[c]->Setup();

  //-Setup of lines
  for(unsigned l=0; l<LineCount; l++) {
    ILine* line=Lines[l]; //current line    
    if(Debug) Log->PrintfDbg("\tLine: %d\n",l);
    IBody* body=NULL;
    if(line->GetUseFair())body=GetBody(line); //-Gets its body
    line->Setup(Dirout,Env,body);
    line->SetupWaves(ucurr,0.0);
    for(unsigned c=0; c<line->GetNConnections(); c++) {
      IConnection* conn=line->GetConnections()[c]; //current IConnection
      conn->SetEnv(Env);
      if((conn->GetType()<TpConnection::fixed)||(conn->GetType()>TpConnection::connect)){
        Run_Exceptioon("Error: Wrong IConnection type found: "+ fmdp::ToString(conn->GetType())+" in ILine: "+fmdp::ToString(line->GetId())+ ".");
      }
      conn->AddLineToConnect(line);
    }
  }
  //- Check if dtm was adjusted
  if(Env->GetDtMmod()) Log->PrintfWarning("Adjusted dtM to %g for fairlead calculations.",Env->GetDtM0());

  //-Setup of Bodies  
  for(unsigned m=0; m<BodyCount; m++) {
    IBody* body=Bodies[m]; //current body
    std::vector<ILine*> lines=GetLines(body->GetRef());
    body->Setup(Dirout,Env,lines);
    //PRINT: look num body which is computed
    if(Debug) Log->PrintfDbg("\nBody: %d\n",m);
  }

  //-Setup of Outputs
  if(OutP){OutP->Setup(Dirout);}

  //-Attaches the lines with connects
  AddAttachedLines();

  //Free array
  delete[] ucurr;
}

//==================================================================================
/// Prepares the States
//==================================================================================
void MoorDynPlus::PrepareState() {
  unsigned n=ConnCount*NStateVariables;
  for(unsigned l=0; l<LineCount; l++) {
    ILine* line=Lines[l];
    line->SetIndex(n);      // assign start index of each line
    n+=(line->GetN()-1)*6;  // add 6 state variables for each internal node
  }
  // make state vector  
  NX=n;  // size of state vector array
  if(misc::wordy>1) Log->Printf("Creating state vectors of size %d \n",NX);

  AllocateMemoryIntegration(NX);  //reserves memory for integration arrays arrays with size NX

  // make array used for passing fairlead kinematics and forces between fairlead- and platform-centric interface functions
  Ffair  =fmdp::GetArray2D(FairCount,3);
  rFairi =fmdp::GetArray2D(FairCount,3);
  rdFairi=fmdp::GetArray2D(FairCount,3);
}

//==================================================================================
/// Inserts in a vector the connections positions
//==================================================================================
// void MoorDynPlus::StoreConnPositions() {
//   const TpConnection type=TpConnection::vessel; // Fair connections
//   unsigned count=0; // Number of connections to store in array
//   //-Select the Number of connections of the type introduced
//   if     (type==TpConnection::fixed)   { count=AnchCount; } //Fixed
//   else if(type==TpConnection::vessel)  { count=FairCount; } //Vessel
//   else if(type==TpConnection::connect) { count=ConnCount; } //Connects

//   rFairtS=fmdp::GetArray2D(count,3); //returns a double** 
//   std::vector<IConnection*> connections=GetConnections(type); // returns a pointer to connections of type selected
  
//   for(unsigned c=0;c<count;c++) {
//     IConnection* connection=connections[c];
//     rFairtS[c][0]=connection->GetX();
//     rFairtS[c][1]=connection->GetY();
//     rFairtS[c][2]=connection->GetZ();
//   }
// }

//==================================================================================
/// Initialize system,including trying catenary IC gen of Lines
//==================================================================================
void MoorDynPlus::InitializeSystem(double X[],double XD[]) {
  double* x=X;
  // double* xd=XD;

  misc::RotMat(x[3],x[4],x[5],TransMat);
  std::vector<IConnection*> connAnch=GetConnections (TpConnection::fixed);
  std::vector<IConnection*> connFair=GetConnections (TpConnection::vessel);
  std::vector<IConnection*> connConns=GetConnections(TpConnection::connect);

  //-Set positions of fairleads based on inputted platform position
  for(unsigned f=0; f<FairCount; f++) connFair[f]->InitializeFairlead(TDouble3(x[0],x[1],x[2]),TransMat); //FairLeads

  //-For connect types,write the coordinates to the state vector
  for(unsigned c=0; c<ConnCount; c++)  connConns[c]->InitializeConnect(X0+6*c); //Connects. 6 for each connect

  //-Go through lines and Initialize internal node positions using quasi-static model
  for(unsigned l=0; l<LineCount; l++) {
    ILine* line=Lines[l];
    line->Initialize(X0+line->GetIndex());   // lines
  }
}

//==================================================================================
/// Finalizing ICs using dynamic relaxation
//==================================================================================
void MoorDynPlus::DoDynamicRelaxIC() {
  Log->Printf("\nFinalizing ICs using dynamic relaxation( %.2f X normal drag ).",Env->GetICDfac());

  for(unsigned l=0; l<LineCount; l++) { Lines[l]->ScaleDrag(Env->GetICDfac()); } // boost drag coefficient

  int niic=(int)round(Env->GetICTmax()/Env->GetICdt());  // max Number of IC gen time steps
  tdouble3 ffair; // array to temporarily store fairlead force components
  double* tensions=fmdp::GetArray(FairCount*  3*  niic);  // pointer to double for store tensions for analyzing convergence
  double* fairTens=fmdp::GetArray(FairCount); // pointer to double for store tensions for analyzing convergence
  double* fairTensLast=fmdp::GetArray(FairCount); // pointer to double for store tensions for analyzing convergence
  double* fairTensLast2=fmdp::GetArray(FairCount); // pointer to double for store tensions for analyzing convergence

  unsigned lf=0; // Fairlead index
  float percent=0.0; // Percent of Fairlead tensions convergence
  unsigned progress; // Used to store the current progress of the simulation
  double t=0.0; // Convergence current time
  unsigned exp=0; // Exponent to round the Percent
  // round to get appropriate mooring model time step
  const int NdtM=(int)(ceil(Env->GetICdt()/Env->GetDtM0())); // Number of mooring model time steps per outer time step(10%) to acelerate
  unsigned step=(unsigned)(NdtM >= 1e6 ? NdtM*1e-1: NdtM*25e-2); //Adjust the step to print progress during the loop for IC generation
  const double dtM=Env->GetICdt()/NdtM;  // mooring model time step size (s)
  int iic; //Main loop iterator
  std::vector< std::vector<double>> fairTension; //Store [time,FairTension] for each step
  // loop through IC generation time analysis time steps
  for(iic=0; iic<niic; iic++) {
    progress=iic*100/niic;  // Get current generation progress
    IApi::ShowProgressBar(progress);    // Show progress by console
    t=iic*Env->GetICdt();// IC gen time(s).  << is this a robust way to handle time progression?
    //#ifdef _OPENMP
    //#pragma omp parallel for schedule(static)
    //#endif
    // loop through line integration time steps
    for(int its=0; its<NdtM; its++) {
      if(its%step==0) IApi::ShowProgressBar((unsigned)(progress+((its*(100.f/niic)/float(NdtM)))));
      Rk2(&t,dtM);// call RK2 time integrator (which calls the model)  
    }   

    // store previous fairlead tensions for comparison
    for(lf=0; lf<FairCount; lf++) {
      fairTensLast2[lf]=fairTensLast[lf];
      fairTensLast[lf]=fairTens[lf];
    }

    std::vector<IConnection*> conFair=GetConnections(TpConnection::vessel); //Returns Fairlead Connections
    // go through connections to get fairlead forces
    for(lf=0; lf<FairCount; lf++) {
      ffair=conFair[lf]->GetFnet();
      fairTens[lf]=sqrt(fgeo::ProductScalar(ffair,ffair));
    }
    //Store fairlead tensions for print by console
    int size=static_cast<int>(fairTension.size());
    fairTension.resize(size+1,std::vector<double>(2,0.0));
    fairTension[size][0]=t;
    fairTension[size][1]=fairTens[0];

     // check for convergence (compare current tension at each fairlead with previous two values)
    if(iic>2) {
      for(lf=0; lf<FairCount; lf++)if((abs(fairTens[lf]/fairTensLast[lf]-1.0)>Env->GetICthresh()) ||(abs(fairTens[lf]/fairTensLast2[lf]-1.0)>Env->GetICthresh()))break;
      if(lf==FairCount) {  // if we made it with all cases satisfying the threshold
        percent=float((100.0f*Env->GetICthresh())); 
        exp=(unsigned)(round(log10(percent)));
        break; // break out of the time stepping loop
      }
    }
  }
  IApi::ShowProgressBar();
  if(lf==FairCount) Log->Printf("\nFairlead tensions converged to %.*f %% after %.f seconds.",(exp*(-1)),percent,t);  
  for(unsigned l=0; l<LineCount; l++) {
    Lines[l]->ScaleDrag(1.0/Env->GetICDfac()); // restore drag coefficients
    Lines[l]->SetTime(0.0);    // reset time to zero so first output line doesn't start at>0 s
  }
  for(int i=0; i<static_cast<int>(fairTension.size()); i++) Log->Printf("  t=%.f s,tension at first fairlead is %.4f N.",fairTension[i][0],fairTension[i][1]);   // write status update and send cursor back to start of line

  // @mth: new approach to be implemented
  //  // ------------------------- calculate wave time series if needed -------------------
  //  if(env.WaveKin==2)
  //  {
  //    for(unsigned l=0; l<LineCount; l++) 
  //      LineList[l].makeWaveKinematics( 0.0 );
  //  }
  AllOutput(0,dtM);   // writes outputs for t=0

  //Free and clear vectors
  fairTension.clear();
  delete[] tensions;      delete[] fairTens;
  delete[] fairTensLast;  delete[] fairTensLast2;
}

//==============================================================================
/// Runge-Kutta 2 integration routine (integrates X0 and time)
//==============================================================================
void MoorDynPlus::Rk2(double* t0,const double dt){
  std::string function="Rk2";
  std::string textExc="NaN value detected in MoorDynPlus state at dynamic relaxation time";
  unsigned step;
  //-Gets derivatives at t0. F0=f( t0,x0 );
  step=1;RHSmaster(X0,F0,*t0,dt);         
  //-Check for Not-A-Number
  if(fmdp::CheckForNan(X0,NX)) Run_Exceptioon(fun::PrintStr("%d %s %g s.",step,textExc.c_str(),*t0));
  //-Integrates to t0 +dt/2. x1=x0+dt*F0/2.0;
  for(unsigned i=0; i<NX; i++) Xt[i]=X0[i]+0.5*dt*F0[i];  
  //-Gets derivatives at t0 +dt/2.  F1=f( t1,x1 );
  step=2;RHSmaster(Xt,F1,*t0+0.5*dt,dt);
  //-Integrates X0 to t0+dt
  for(unsigned i=0; i<NX; i++) X0[i]=X0[i]+dt*F1[i]; 
  //-Check for Not-A-Number
  if(fmdp::CheckForNan(X0,NX)) Run_Exceptioon(fun::PrintStr("%d %s %g s.",step,textExc.c_str(),*t0));
  //-Update time   
  *t0=*t0+dt;    
}

//==============================================================================
/// Master function to handle time stepping (updated in v1.0.1 to follow MoorDynPlus F)
//==============================================================================
void MoorDynPlus::RHSmaster(const double X[],double Xd[],const double t,const double dt){
  //std::vector<IConnection*> connAnch=GetConnections(TpConnection::fixed  ); ///<Anchors connections //This kind of connection doesn't used here
  std::vector<IConnection*> connFairs=GetConnections (TpConnection::vessel ); ///<Fairs connections
  std::vector<IConnection*> connConns=GetConnections (TpConnection::connect); ///<Connects connections
  std::vector<IConnection*> connections=GetConnections(); ///<All connections
  // extrapolate instaneous fairlead positions
  for(unsigned l=0; l<FairCount; l++) { connFairs[l]->UpdateFairlead(t); }
  // calculate forces on all connection objects  
  for(unsigned l=0; l<ConnectionCount; l++) { connections[l]->GetNetForceAndMass(); }
  // calculate connect dynamics(including contributions from latest line dynamics,above,as well as hydrodynamic forces)
  for(unsigned l=0; l<ConnCount; l++) { connConns[l]->DoRHS((X+6*l),(Xd+6*l),t); }
  // calculate line dynamics
  for(unsigned l=0; l<LineCount; l++) {
    ILine* line=Lines[l];
    if(!line->GetDisabled())line->DoRHS((X+line->GetIndex()),(Xd+line->GetIndex()),t,dt);
  }
}

//==================================================================================
///Function to Configure the initial boundary conditions.
/// Used by external software coupling
//==================================================================================
int MoorDynPlus::LinesInit(double X[],double XD[],const std::string filexml
    ,const std::string nodexml,const char* dir
    ,const tfloat3 gravity,const double tmax
    ,const double dtout,const unsigned ftmkbound[]
    ,const unsigned numFts) {
  std::string function="LinesInit";
  std::string dir_in=dir;
  if((numFts==0) || (!ftmkbound)) Run_Exceptioon("The number of driven-object floatings cannot be 0."); 
  if(dir_in!="") Dirout=dir;
  else CreateOutDirectory();
  Log->Print(IApi::GetLicense());
  FileIN=filexml;
  RootNodeXml=nodexml;
  TimeMax=tmax;
  DtOut=dtout;
  if(!XmlReaded) LoadXml();  //Check if Xml was read before
  Env->SetG(gravity);
  Env->SetTimeMax(tmax);
  Setup();
  CheckFtIds(ftmkbound,numFts);
  PrepareState();
  //StoreConnPositions();
  InitializeSystem(X,XD);
  VisuConfig();
  DoDynamicRelaxIC();
  VisuIniTen();
  return 0;
}

//==================================================================================
/// This function now handles the assignment of fairlead boundary conditions, time stepping, and collection of resulting forces at fairleads
/// It can also be called externally for fairlead-centric coupling.
//==================================================================================
int MoorDynPlus::FairleadsCalc(const unsigned numFts,double*** rFairIn,double*** rdFairIn,double*** fFairIn,double* t_in,double* dt_in) {
  std::string function="FairleadsCalc";

  double t=*t_in;    // this is the current time
  const double dtC=*dt_in;  // this is the coupling time step

  if(dtC>0) { // if DT>0, do simulation, otherwise leave passed fFairs unadjusted.
    // send latest fairlead kinematics to fairlead objects
    for(unsigned b=0; b<numFts; b++){
      CheckBody(b);
      const IBody* body=Bodies[b];
      std::vector<IConnection*> connFairs=body->GetFairs(); //Get Fair connections of IBody passed
      const unsigned nFairs_m=body->GetNFairs(); //Get num of Fairs connections of this body
      for(unsigned l=0; l<nFairs_m; l++) {
        const tdouble3 rfairin=TDouble3(rFairIn[b][l][0],rFairIn[b][l][1],rFairIn[b][l][2]);
        const tdouble3 rdfairin=TDouble3(rdFairIn[b][l][0],rdFairIn[b][l][1],rdFairIn[b][l][2]);
        if(!connFairs[l]->GetDisabled()) connFairs[l]->InitiateStep(rfairin,rdfairin,t); 
      }
    }
    // round to get appropriate mooring model time step
    const int NdtM=(int)ceil(dtC/Env->GetDtM0());   // Number of mooring model time steps per outer time step
    if(NdtM<1) Run_Exceptioon(fun::PrintStr("The coupling time step is lower than dtM. (%g<%g).",dtC,Env->GetDtM0()));

    const double dtM=dtC/NdtM;// mooring model time step size (s)

    // loop through line integration time steps (integrate solution forward by dtC)
    for(int its=0; its<NdtM; its++) { Rk2(&t,dtM); } // call RK2 time integrator (which calls the model)
    //for(unsigned nl=0;nl<LineCount;nl++){ if(!Lines[nl]->GetDisabled()&&Lines[nl]->GetBreakTension())Lines[nl]->CheckTension();}                         
    for(unsigned b=0; b<numFts; b++){
      // go through connections to get fairlead forces  
      CheckBody(b);
      const IBody* body=Bodies[b];
      std::vector<IConnection*> connFairs=body->GetFairs(); //Get Fair connections of IBody passed
      const unsigned nFairs_m=body->GetNFairs(); //Get num of Fairs connections of this body  
      for(unsigned l=0; l<nFairs_m; l++) { 
        if(!connFairs[l]->GetDisabled()){
          tdouble3 ffairin=connFairs[l]->GetFnet();   
          fFairIn[b][l][0]=ffairin.x;
          fFairIn[b][l][1]=ffairin.y;
          fFairIn[b][l][2]=ffairin.z;
        }
      }
    }
    AllOutput(t,dtM);   // write outputs
  }
  return 0;
}

//==================================================================================
/// Checks if floating id and body exists. Launches exception in error case
//==================================================================================
void MoorDynPlus::CheckBody(const unsigned ftid)const {
  std::string function="CheckBody";
  if((ftid<0) ||(ftid>=BodyCount)) Run_Exceptioon(fun::PrintStr("Body id=%d is out of range. Expected [0,%d].",ftid,BodyCount-1));
  //-If the body does not exist
  if(!Bodies[ftid]) Run_Exceptioon(fun::PrintStr("Body id=%d does not exist.",ftid));;
}

//==================================================================================
/// Checks if floating and line id exists. Launches exception in error case
//==================================================================================
void MoorDynPlus::CheckBodyLine(const unsigned ftid,const unsigned lineid)const {
  std::string function="CheckBodyLine";
  CheckBody(ftid);
  CheckLine(lineid);
}

//==================================================================================
/// Checks if floating id and body exists. Returns true in error case
//==================================================================================
void MoorDynPlus::CheckLine(const unsigned lineid)const {
  std::string function="CheckLine";
  if((lineid<0) || (lineid>=LineCount)) Run_Exceptioon(fun::PrintStr("Line id=%d is out of range. Expected [0,%d].",lineid,LineCount-1));
  const ILine* line=Lines[lineid]; 
  //-If the line doesn't exist
  if(!line) Run_Exceptioon(fun::PrintStr("Line id=%d does not exist.",lineid));
}

//==================================================================================
///Returns the number of fairleads (Vessel connections) of a Mooring created
//==================================================================================
unsigned MoorDynPlus::GetNFairs(const unsigned ftid)const {
  std::string function="GetNFairs";
  CheckBody(ftid);
  unsigned numF=0;
  const IBody* body=Bodies[ftid];
  numF=body->GetNFairs();

  return numF;
}

//==================================================================================
/// Returns the Number of segments of the line of floating selected
//==================================================================================
unsigned MoorDynPlus::GetSegsCount(const unsigned lineid)const {
  std::string function="GetSegsCount";
  CheckLine(lineid);
  unsigned numSeg=0;
  const ILine* line=Lines[lineid]; 
  numSeg=line->GetN();
  return numSeg;
}

//==================================================================================
/// Returns the Number of segments of the line of floating selected
//==================================================================================
unsigned MoorDynPlus::GetSegsCount(const unsigned ftid,const unsigned lineid)const {
  std::string function="GetSegsCount";
  CheckBodyLine(ftid,lineid);
  unsigned numSeg=0;
  const IBody* body=Bodies[ftid];   
  const ILine* line=body->GetLines()[lineid]; 
  numSeg=line->GetN();
  return numSeg;
}

//==================================================================================
/// Returns the Tension of the Fair Selected
//==================================================================================
double MoorDynPlus::GetFairTen(const unsigned lineid)const {
  std::string function="GetFairTen";
  CheckLine(lineid);
  const ILine* line=Lines[lineid]; 
  const int nNodes=line->GetN();
  const double fairTen=line->GetNodeTen(nNodes);
  return fairTen;
}

//==================================================================================
/// Writes the position of Node selected in pos
//==================================================================================
int MoorDynPlus::GetNodePos(const unsigned lineid,const int nodeid,double pos[3]) {
  std::string function="GetNodePos";
  CheckLine(lineid);
  const tdouble3 nodepos=GetNodePos(lineid,nodeid);
  pos[0]=nodepos.x;
  pos[1]=nodepos.y;
  pos[2]=nodepos.z;
  return 0;
}

//==================================================================================
/// Returns the position of Node selected
//==================================================================================
tdouble3 MoorDynPlus::GetNodePos(const unsigned lineid,const int nodeid)const {
  std::string function="GetNodePos";
  CheckLine(lineid);
  const ILine* line=Lines[lineid]; 
  return line->GetNodePos(nodeid);  // call line member function to fill in coordinates of the node of interest 
}

//==================================================================================
/// Returns the position of Node selected
//==================================================================================
int MoorDynPlus::GetNodePosLink(const unsigned ftid,const unsigned lineid,double pos[3]) {
  std::string function="GetNodePosLink";
  CheckBodyLine(ftid,lineid);
  const IBody* body=Bodies[ftid];   
  ILine* line=body->GetLines()[lineid]; 
  std::vector<IConnection*> fairs=line->GetFairs(body->GetRef()); //Gets the fairleads connected to this body
  const unsigned nodeid=fairs[0]->GetNode();//Gets the node number of the fairlead in the line
  GetNodePos(lineid,nodeid,pos);
  return 0;
}

//==================================================================================
/// Returns the position of Node selected
//==================================================================================
tdouble3 MoorDynPlus::GetNodePosLink(const unsigned ftid,const unsigned lineid)const {
  std::string function="GetNodePosLink";
  CheckBodyLine(ftid,lineid);
  const IBody* body=Bodies[ftid];
  ILine* line=body->GetLines()[lineid]; 
  std::vector<IConnection*> fairs=line->GetFairs(body->GetRef()); //Gets the fairleads connected to this body
  const unsigned nodeid=fairs[0]->GetNode();//Gets the node number of the fairlead in the line
  return line->GetNodePos(nodeid);
}

//=================================================================================
///Returns the body reference of the ftid passed  
//=================================================================================
unsigned MoorDynPlus::GetBodyReference(const unsigned ftid)const {
  std::string function="GetBodyReference";
  CheckBody(ftid);
  const IBody* body=Bodies[ftid];   
  const unsigned ref=body?body->GetRef():UINT_MAX;
  return ref;
}
//==================================================================================
/// Calls to writer function of the IBody passed
//==================================================================================
void MoorDynPlus::AllOutput(double t,double dtC) {
  //-Writes IConnection output of lines
  if(OutP){
    const double start=OutP->GetStartTime();
    const double finish=OutP->GetEndTime();
    const double next=OutP->GetNextTime();
    OutP->SaveCsv(t,dtC);
    //-Writes individual line output files (of Nodes)
    for(unsigned l=0; l<LineCount; l++)
      Lines[l]->WriteOutput(t,dtC,start,finish,next);
  }
  if(!Coupled) SaveVtk(t,dtC);
}

//==============================================================================
/// Write the vtk for each mooring line.
//==============================================================================
void MoorDynPlus::SaveVtk(double timestep,double dt) {
#ifdef DDISABLE_DSPH
  const double start=0.;  //-Start time. | Tiempo de inicio.
  const double finish=TimeMax;   //-End time. | Tiempo de finalizacion.
  bool savedata=false;
  if(start<=timestep && timestep<=finish) savedata=NextTime<=timestep;
  if(savedata){
    JSpVtkShape ss;
    unsigned vpossize=0;
    tdouble3* vpos=NULL;
    const unsigned nlines=LineCount;
    for(unsigned cl=0;cl<nlines;cl++){
      const unsigned nodes=GetSegsCount(cl)+1;
      if(nodes>vpossize){
        vpossize=nodes;
        vpos=fun::ResizeAlloc(vpos,0,vpossize);
      }
      for(unsigned cn=0;cn<nodes;cn++)vpos[cn]=GetNodePos(cl,cn);
      ss.AddLines(nodes,vpos,word(cl));
      ss.AddPoints(nodes,vpos,word(cl));
    }
    delete[] vpos; vpos=NULL;
    //-Genera fichero VTK.
    Log->AddFileInfo(Dirout+"MooringsVtk/Moorings_????.vtk","Saves VTK file with moorings.");
    const std::string file=Dirout+fun::FileNameSec("MooringsVtk/Moorings.vtk",Part);
    ss.SaveVtk(file,"Line");
    Part++;
  }
#endif // DDISABLE_DSPH
}

//==================================================================================
/// Show the initial tensions in the nodes
//==================================================================================
void MoorDynPlus::VisuIniTen()const{
  Log->Print("");
  Log->Printf("Initial node tensions:");
  for(unsigned l=0;l<LineCount;l++)Lines[l]->VisuIniTen();
  Log->Print("");
}

//==================================================================================
/// Show the set up of the lines
//==================================================================================
void MoorDynPlus::VisuConfig()const{
  //-System
  Log->Printf("Creating mooring system:");
  Log->Printf("  Bodies......: %d",BodyCount);
  Log->Printf("  Lines.......: %d",LineCount);
  Log->Printf("  Connections.: %d",ConnectionCount);
  Log->Printf("    Fairleads.: %d",FairCount); 
  Log->Printf("    Anchors...: %d",AnchCount); 
  Log->Printf("    Connects..: %d",ConnCount);
  //-IEnvironment
  Env->VisuConfig();
  //-Lines
  Log->Printf("Line properties:");
  for(unsigned l=0;l<LineCount;l++)Lines[l]->VisuConfig();
}

//==================================================================================
/// Function for finish the program and free the memory
//==================================================================================
int MoorDynPlus::LinesClose() {
  //PRINT
  Log->Print("\nLines Close.");
  return 0;
}

//==================================================================================
/// Main function of MoorDynPlus. Order the necessaries executions. 
/// This function is executed when this library isn't coupling to external software.
//==================================================================================
void MoorDynPlus::Run() {
  //Only calls CreateOutDirectory() function if executes this library how module. 
  //When is coupling with DSPH,use the DSPH output directory
  Coupled=false;
  CreateOutDirectory();
  if(!Log)CreateLog();
  Log->Print(IApi::GetLicense());
  if(!XmlReaded) LoadXml(); 
  Setup();
  PrepareState();
  //StoreConnPositions();
  //****** For init the lines ******
  double* X=new double[6]; //Positions
  double* XD=new double[6]; //Velocities
  std::memset(X,0,sizeof(double)*6);
  std::memset(XD,0,sizeof(double)*6);
  //********************************
  InitializeSystem(X,XD);
  VisuConfig();
  DoDynamicRelaxIC();
  VisuIniTen();
  //*********************************
  AllocateArrays();
  //*********************************** FAIRLEADS CALC ***********************************
  double t=NextTime=0; ///<Initial time
  TimeMax=OutP->GetEndTime(); ///<Initial time
  DtOut=OutP->GetDtOut(); ///<Out time step
  double dtm=Env->GetDtM0(); ///<Time step
  const double divisor=0.00001; ///<Used for increment positions and to make the movement
  while(t<Env->GetTimeMax()+dtm) {//loop time
    for(unsigned c=0;c<BodyCount;c++) {
      const IBody* b=Bodies[c];
      const unsigned nFairs=b->GetNFairs();
      std::vector<IConnection*> fairs=b->GetFairs();
      for(unsigned cp=0;cp<nFairs;cp++) {
        // const unsigned ptid=fairs[cp]->GetId();
        const tdouble3 pos=fairs[cp]->GetPositions();
        const tfloat3  vel=TFloat3(0);
        FairleadPos[c][cp][0]=pos.x;  FairleadPos[c][cp][1]=pos.y;  FairleadPos[c][cp][2]=pos.z;
        FairleadVel[c][cp][0]=vel.x;  FairleadVel[c][cp][1]=vel.y;  FairleadVel[c][cp][2]=vel.z;
      }
    }
    FairleadsCalc(BodyCount,FairleadPos,FairleadVel,FairleadFor,&t,&dtm);
    UpdatePos(divisor); //Update positions

    if((NextTime<=t) || (t>=Env->GetTimeMax())) {
      PrintProgress(t);  // if output should print to console
      NextTime+=DtOut;
    }
    t+=dtm;
  }
  delete[] X; delete[] XD;

  LinesClose();
  Log->Print("\n Execution completed successfully. Press any key to continue..."); getchar();
  return;
}

//==================================================================================
/// Reserves memory for the arrays of positions, velocities and forces.
//==================================================================================
void MoorDynPlus::AllocateArrays(){
  FairArrays=true;
  try{
    FairleadPos=new double**[BodyCount];
    FairleadVel=new double**[BodyCount];
    FairleadFor=new double**[BodyCount];
    for(unsigned c=0;c<BodyCount;c++) {
      IBody* b=Bodies[c];
      unsigned nFairs=b->GetNFairs();
      FairleadPos[c]=new double*[nFairs];
      FairleadVel[c]=new double*[nFairs];
      FairleadFor[c]=new double*[nFairs];
      for(unsigned cp=0;cp<nFairs;cp++) {
        FairleadPos[c][cp]=new double[3];
        FairleadVel[c][cp]=new double[3];
        FairleadFor[c][cp]=new double[3];
      }
    }
  }catch(std::bad_alloc const&) {
      Run_Exceptioon("Could not allocate the requested memory.");
    }
}
//==================================================================================
/// Function for simulated the driven-object floating
//==================================================================================
void MoorDynPlus::UpdatePos(const double divisor) {
  for(unsigned b=0; b<BodyCount; b++) {
    const IBody* body=Bodies[b];
    for(unsigned l=0; l<body->GetNFairs(); l++) {
      IConnection* fairlead=body->GetFairs()[l];
      const tdouble3 currentPos=fairlead->GetPositions();
      fairlead->SetPositions(TDouble3(currentPos.x,currentPos.y,(currentPos.z+divisor)));
    }
  }
}

//==============================================================================
/// Prints the execution progress in the console
//==============================================================================
void MoorDynPlus::PrintProgress(double time) {
  if(!HeaderProg) {
    PrintHeader();
    HeaderProg=true;
  }
  Log->Printf("%.2f\t\t%.2f\t\t%.2f %%",time,Env->GetTimeMax(),PrintPercent(time));
}

//==============================================================================
/// Prints the header for the execution progress
//==============================================================================
void MoorDynPlus::PrintHeader() {
  Log->Print("\n==================================================");
  Log->Print("Time [s]\tTimeMax [s]\tProgress [%]");
  Log->Print("==================================================");
}

//==============================================================================
/// Prints the completed perecent of the simulation
//==============================================================================
double MoorDynPlus::PrintPercent(double time) {
  return(100*  time)/Env->GetTimeMax();
}