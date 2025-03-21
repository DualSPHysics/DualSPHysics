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

/// \file IConnection.cpp \brief Implements the class \ref IConnection.

#include "IConnection.h"
#include "IEnvironment.h"
#include "ILine.h"
#include "FunMoorDynPlus.h"
#include "Functions.h"
#include "FunctionsMath.h"

//==============================================================================
/// Constructor.
//==============================================================================
IConnection::IConnection(JLog2* log,TpConnection tpc):Log(log),Type(tpc){
	Reset();
	ClassName="IConnection";
}

//==============================================================================
/// Destructor.
//==============================================================================
IConnection::~IConnection(){
	Reset();
	Log=NULL;
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void IConnection::Reset() {
  UseMass=Disabled=UseVol=false;

  TimeStep=0; NumExec=0;
  ConM=0; ConV=0; ConCdA=0; ConCa=0;
  Ref=0; AttachCount=0; AttachCount=0;

  Pos=TDouble3(0);  ForceIni=TDouble3(0);
  R=TDouble3(0);    R_ves=TDouble3(0); 
  Rd=TDouble3(0);   Rd_ves=TDouble3(0);	
  Fnet=TDouble3(0); Fnet_i=TDouble3(0);
  Q=TDouble3(0);    RHS=TDouble3(0);    

  M=M_i=S=TMatrix3d(0);

  Attached.clear(); Top.clear();
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void IConnection::LoadXml(JXml* sxml,const std::string& place,TiXmlElement* elec,const int id) {
	std::string function="LoadXml";
	//TiXmlNode* node=sxml->GetNode(place,false);
	//if(!node)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
	if(!elec)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
	ReadXml(sxml,elec,id);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void IConnection::ReadXml(JXml* sxml,TiXmlElement* eleConn,const int id) {
	std::string function="ReadXml";
	//-Loads connections zones.
	std::string value=eleConn->Value(); // Text of connection tag

	//-If the current connection tag exists
	if(eleConn) {
		//-Stores the IConnection Type [Anch=0; Fair=1 ; Con=2]
		if     (Type==TpConnection::vessel)  Ref=sxml->GetAttributeUnsigned(eleConn,"bodyref",false);
		else if(Type==TpConnection::connect) Ref=sxml->GetAttributeUnsigned(eleConn,"conref",false); 

		Idn=id; //Stores a sequential number of connections for each connection
		//PRINT: Look number of connection stored
		if(Debug) Log->PrintfDbg("\t\tnumber of connection: %d\n",Idn);
		
		//-Gets positions
		Pos=sxml->GetAttributeDouble3(eleConn); 

	 	//-Gets Mass and volum attributes
		ConM=sxml->GetAttributeDouble(eleConn,"M",true,0);
		ConV=sxml->GetAttributeDouble(eleConn,"V",true,0);
		if(ConM!=0)UseMass=true;
		if(ConV!=0)UseVol=true;

		//-Gets forces
		//if(sxml->ExistsElement(eleConn,"forceini")){
		//	TiXmlElement* eleConnF=eleConn->FirstChildElement("forceini");
		//	ForceIni=sxml->GetAttributeDouble3(eleConnF);
		//}

		//-Stores mass and drag coefficients (Only type=2 [Connects])
		if(IsConnect()){
			//-Gets Coefficients
			if(sxml->ExistsElement(eleConn,"coeffini")){
				const TiXmlElement* eleCoeff=eleConn->FirstChildElement("coeffini");
				ConCdA=sxml->GetAttributeDouble(eleCoeff,"cda",false);
				ConCa=sxml->GetAttributeDouble(eleCoeff,"ca",false);
			}
		}
	}
}

//==============================================================================
///  Makes the setup for this IConnection. Initializes the position of this
//==============================================================================
void IConnection::Setup() {
	//PRINT: look when IConnection::Setup() is executed
	if(Debug) Log->PrintfDbg("\t\tConnection::setup() - Number: %d\n",Idn); 
  // start off position at that specified in input file 
	// (will be starting point for connect connections
	// and the permanent location of anchor connections.)	
	R=Pos;
}

//==============================================================================
// This function handles assigning a line to a connection node
//==============================================================================
void IConnection::AddLineToConnect(ILine* line) {
	if(misc::wordy>0) Log->Printf("L %u ->N %u \n",line->GetId(),Idn); 
	if(std::find(Attached.begin(),Attached.end(),line)==Attached.end()) {
		//-Attached doesn't contain line
		unsigned top=line->GetN(); // Gets the number of segments
		const IConnection* anch=line->GetAnchConnect();
		if(anch==this)top=0; //If the connect is an "anchor", stores the segment=0
		if(AttachCount<10) { // this is currently just a maximum imposed by a fixed array size.  could be improved.
			Attached.push_back(line);
			Top.push_back(top);
			AttachCount += 1;
		}
	}
}

////==============================================================================
///// Function to return net force on fairlead (just to allow public reading of Fnet)
////==============================================================================
//void IConnection::GetFnet(tdouble3& fnet_out) {
//	fnet_out=Fnet;
//	//fnet[2] += M[0][0]*(-env.g); // add weight  NO this is already in Fnet !!! (removed Oct 20)
//	// should it include inertial "force"?  i.e.	for(int J=0; J<3; J++)  fnet += M[I][J]*acceleration[J] 	
//}

//==============================================================================
/// Initialize Positions and velocities of fairlead. 
/// pX -> Positions sent by LinesInit()
//==============================================================================
void IConnection::InitializeFairlead(const tdouble3 pX,const tmatrix3d& transMat) {
	std::string function="InitializeFairlead";

	if(IsVessel()) { //If is vessel Type
		R.x=transMat.a11*Pos.x+transMat.a12*Pos.y+transMat.a13*Pos.z+pX.x;	// x
		R.y=transMat.a21*Pos.x+transMat.a22*Pos.y+transMat.a23*Pos.z+pX.y;	// y
		R.z=transMat.a31*Pos.x+transMat.a32*Pos.y+transMat.a33*Pos.z+pX.z;	// z
		// also specify the above in the "vessel" arrays that prescribe the kinematics over the following time steps, for IC gen
		R_ves=R; 
		Rd=TDouble3(0);Rd_ves=TDouble3(0);
	}else Run_Exceptioon(fun::PrintStr("Wrong connection Type given. Something's not right.")); 
}

//==============================================================================
/// Initialize the positions and velocities of this IConnection
//==============================================================================
void IConnection::InitializeConnect(double* X) {
	std::string function="InitializeConnect";
	if(IsConnect()) {  // If is connect Type
		// assign initial node kinematics to state vector
		X[3+0]=R.x; X[3+1]=R.y; X[3+2]=R.z;
		X[  0]=Rd.x;X[  1]=Rd.y;X[  2]=Rd.z;
	}else	Run_Exceptioon(fun::PrintStr("Wrong connection Type given. Something's not right.")); 
}

//==============================================================================
/// Helper function to sum forces and mass from attached lines 
/// Used for connect dynamics and fair/anch tensions
//==============================================================================
void IConnection::GetNetForceAndMass() {
	// loop through each connected line,summing to get the final result
	const unsigned Nl=AttachCount;				// Number of attached line segments
	// clear before re-summing	
	Fnet=TDouble3(0);	M=TMatrix3d(0);

	// loop through attached lines
	for(unsigned l=0;l<Nl;l++) {	// get quantities
		if(Top[l]==0) (Attached[l])->GetAnchStuff(Fnet_i,M_i);// if attached to bottom/anchor of a line...
		else 	        (Attached[l])->GetFairStuff(Fnet_i,M_i);// attached to top/fairlead
		// sum quantities
		Fnet=Fnet+Fnet_i;
		M=M+M_i; 
	}
	//-Checks if the Mass and the Volum were initialized and calculates it (only connects)
	if(IsConnect()) CalculateMassVol(); 

	// add Constant quantities for IConnection if applicable from input file
	Fnet.x+=ForceIni.x;
	Fnet.y+=ForceIni.y;
	Fnet.z+=ForceIni.z+ConV*(Env->GetRho_w())*(Env->GetG())-ConM*(Env->GetG());
	M=M+ConM;
}

//==============================================================================
/// This function calculates the Mass and Volumn of the Connect Connections if
/// the user didn't introduce value for them
//==============================================================================
void IConnection::CalculateMassVol(){
	//-Calculates Mass
	if(!UseMass){
		ConM=0;
		unsigned numF=0; //Number of connects that are "Fairs"
		for(unsigned lm=0; lm<AttachCount; lm++){
			const unsigned posm=Top[lm]; //Gets the node of this connection in the current line
			if(posm>0){
				ConM+=Attached[lm]->GetM(posm); 
				numF++;
			} //If the connect is an "fairlead"
		}
		ConM /= numF; //Divides by number of lines that have Connect<------>Anchor Connections
	}
	//-Calculates Volumn
	if(!UseVol){
		ConV=0;
		for(unsigned lv=0; lv<AttachCount; lv++){
			const unsigned posv=Top[lv];//Gets the node of this connection in the current line
			if(posv>0) ConV+=Attached[lv]->GetV(0)*posv*Attached[lv]->GetLength();
		}
	}
	if(NumExec>2){UseMass=UseVol=true; }
	else NumExec++; 
}

//==============================================================================
/// This is the function that updates the States. Also includes hydrodynamic forces
//==============================================================================
void IConnection::DoRHS(const double* X,double* Xd,const double time) {
	std::string function="DoRHS";

	TimeStep=time;

	// below now called by master RHS function, for all connection types
	//GetNetForceAndMass(); // in all cases, call to get net force and mass of connection (especially for returning fairlead forces to FAST)

	// ------ behavior dependent on connect Type -------

 //	if(IsFixed()){ // fixed Type
 //		R[0]=ConX;
 //		R[1]=ConY;
 //		R[2]=ConZ;
 //		for(int I=0; I<3; I++){ Rd[I]=0; }
 //	}else if(IsVessel()) {// vessel (moves with platform)					
 //		// set fairlead position and velocity based on BCs (linear model for now)
 //		for(int J=0; J<3; J++)  {
 //			R[J]=R_ves[J]+Rd_ves[J]*(TimeStep-T0);
 //			Rd[J]=Rd_ves[J];
 //		}		
 //		
 //		// assign States
 //		for(int I=0; I<3; I++)  {
 //			Xd[3+I]=Rd[I];  // velocities - these are unused in integration
 //			Xd[I]=0.;     // accelerations - these are unused in integration
 //		}
 //	}			
	if(IsConnect()) { // "connect" Type
		if(TimeStep==0) {   // with current IC gen approach, we skip the first call to the line objects, because they're set AFTER the call to the connects
			// above is no longer true!!! <<<
			for(int I=0; I<3; I++) {
				Xd[3+I]=X[I];  	// velocities - these are unused in integration
				Xd[I]=0.;     		// accelerations - these are unused in integration
				fmdp::UpdateStates(3,X,Xd,TMatrix3d(0),TDouble3(0));
			}
		}else {
			// from state values, get R and rdot values 
			R =TDouble3(X[3+0],X[3+1],X[3+2]); // get positions
			Rd=TDouble3(X[  0],X[  1],X[  2]); // get velocities

			//cout << "ConRHS: m: " << M[0][0] << ", f: " << Fnet[0] << " " << Fnet[1] << " " << Fnet[2] << endl;
			// add dynamic quantities for connection as specified in input file (feature added 2015/01/15)
			Fnet.x-=0.5*(Env->GetRho_w())*Rd.x*abs(Rd.x)*ConCdA;
			Fnet.y-=0.5*(Env->GetRho_w())*Rd.y*abs(Rd.y)*ConCdA;
			Fnet.z-=0.5*(Env->GetRho_w())*Rd.z*abs(Rd.z)*ConCdA;
			M=M+ConV*(Env->GetRho_w())*ConCa;
			// invert node mass matrix
			S=fmath::InverseMatrix3x3(M);// invert node mass matrix (written to S[3x3]) 
			fmdp::UpdateStates(3,X,Xd,S,Fnet);
		}
	}else	Run_Exceptioon(fun::PrintStr("Wrong connection Type given. Something's not right.")); 
}

//==============================================================================
/// Called at the beginning of each coupling step to update the boundary conditions 
/// (fairlead kinematics) for the proceeding line time steps
//==============================================================================
void IConnection::InitiateStep(const tdouble3 rFairIn,const tdouble3 rdFairIn
   ,const double time) {
	T0=time; // set start time for BC functions
	if(IsVessel()) {  // if is vessel Type 
		// update values to fairlead position and velocity functions (fn of time)
		R_ves=rFairIn;
		Rd_ves=rdFairIn;
	}
	// do I want to get precalculated values here at each FAST time step or at each line time step?
}

//==============================================================================
/// Updates the positions and velocities for the tine passed
//==============================================================================
void IConnection::UpdateFairlead(const double time) {
	std::string function="UpdateFairlead";
	TimeStep=time;
	if(IsVessel()) { // vessel (moves with platform)					
		// set fairlead position and velocity based on BCs (linear model for now)
		R=R_ves+Rd_ves*(time-T0);
		Rd=Rd_ves;
	}else	Run_Exceptioon(fun::PrintStr("Wrong connection Type given. Something's not right.")); 
}

//==============================================================================
/// Returns the connection Type
//==============================================================================
std::string IConnection::GetTypeStr()const {
	std::string tp="";
	switch (Type)
	{
		case fixed  : tp="fixed" ; break;
		case vessel : tp="vessel"; break;
		case connect: tp="connect";break;
		default: Run_Exceptioon(fun::PrintStr("Wrong connection Type given. Something's not right.")); 	break;
	}
	return (tp);
}

//==============================================================================
/// Visualises the properties of the connection
//==============================================================================
void IConnection::VisuProperties()const{
	printf("Connection=%u:\n",Idn);
	printf("  R...:"); printf(" (%g,%g,%g)\n",  R.x,  R.y,   R.z);
	printf("  Rd..:"); printf(" (%g,%g,%g)\n",  Rd.x,Rd.y,Rd.z);
	printf("  Q...:"); printf(" (%g,%g,%g)\n",   Q.x,   Q.y,   Q.z);
	printf("  Fnet:"); printf(" (%g,%g,%g)\n",Fnet.x,Fnet.y,Fnet.z);
}

// new function to draw instantaneous line positions in openGL context
#ifdef USEGL
void IConnection::drawGL(void)
{
	double radius=pow(ConV/(4/3*PI),0.33333);  //ConV
	Sphere(R.x,R.y,R.z,radius);
}
#endif