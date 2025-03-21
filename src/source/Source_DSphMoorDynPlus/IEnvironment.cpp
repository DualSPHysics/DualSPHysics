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

/// \file IEnvironment.cpp \brief Implements the class \ref IEnvironment.

#include "IEnvironment.h"
#include <cmath>        
#include <cstdlib>

//==============================================================================
/// Constructor.
//==============================================================================
IEnvironment::IEnvironment(JLog2* log):Log(log) {
	Reset();
	ClassName="IEnvironment";
}

//==============================================================================
/// Destructor.
//==============================================================================
IEnvironment::~IEnvironment() {
	Reset();
	Log=NULL;
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void IEnvironment::Reset() {
	G=0.0; // gravity
	TimeMax=0.0; // Time of simulation
	WtrDpth=0.0; // Water depth
	Rho_w=0.0; // density
	Kb=0.0; // bottom stiffness (Pa/m)
	Cb=0.0; // bottom damping   (Pa/m/s)
	WaveKin=0; // wave kinematics flag (0=off, >0=on)
	WriteUnits=0; // a global switch for whether to show the units line in the output files (1, default), or skip it (0)
	FricCoeff=0.0; // general bottom friction coefficient, as a start
	FricDamp=0.0; // a damping coefficient used to model the friction at speeds near zero
	StatDynFricScale=0.0; // a ratio of static to dynamic friction (=mu_static/mu_dynamic)
	ICDfac=0.0; // factor by which to boost drag coefficients during dynamic relaxation IC generation
	ICdt=0.0; // convergence analysis time step for IC generation
	ICTmax=0.0; // max time for IC generation
	ICthresh=0.0; // threshold for relative change in tensions to call it converged
	DtM0=0.0; // value for desired mooring model time step
	FreeSurface=0.0;
	DtMmod=false;
	DtMauto=true;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void IEnvironment::LoadXml(JXml* sxml,const std::string& place) {
	TiXmlNode* node=sxml->GetNode(place,false);
	if(!node)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
	ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void IEnvironment::ReadXml(JXml* sxml,TiXmlElement* lis) {
	const std::string function="ReadXml";
	//-Loads solverOptions zone.
	TiXmlElement* ele=lis->FirstChildElement("solverOptions");
	if(ele) { //Only one solverOptions
#ifdef DISABLE_DSPH
		sxml->CheckElementNames(ele,false,"gravity timeMax waterDepth kBot cBot waveKin writeUnits frictionCoefficient fricDamp statDynFricScale rho cdScaleIC dtIC tmaxIC threshIC dtM freesurface");
		G=std::abs(sxml->ReadElementDouble(ele,"gravity","value",true,9.81));
		TimeMax=sxml->ReadElementDouble(ele,"timeMax","value",true,UINT32_MAX);  // Time of simulation	
#else
  		sxml->CheckElementNames(ele,false,"waterDepth kBot cBot waveKin writeUnits frictionCoefficient fricDamp statDynFricScale rho cdScaleIC dtIC tmaxIC threshIC dtM freesurface");
#endif // DISABLE_DSPH

		WtrDpth=sxml->ReadElementDouble(ele,"waterDepth","value",false);
		Kb=sxml->ReadElementDouble(ele,"kBot","value",true,3.0e6);
		Cb=sxml->ReadElementDouble(ele,"cBot","value",true,3.0e5);
		std::string writeUnits_str=sxml->ReadElementStr(ele,"writeUnits","value",true,"yes");
		std::transform(writeUnits_str.begin(),writeUnits_str.end(),writeUnits_str.begin(),::tolower);
		const std::size_t yes=writeUnits_str.find("y");
		const std::size_t no=writeUnits_str.find("n");
		WriteUnits=(yes!=std::string::npos?1:no!=std::string::npos?0:1);
		FricCoeff=sxml->ReadElementDouble(ele,"frictionCoefficient","value",true,0.0); // general bottom friction coefficient,as a start
		FricDamp=sxml->ReadElementDouble(ele,"fricDamp","value",true,200.0); // a damping coefficient used to model the friction at speeds near zero
		StatDynFricScale=sxml->ReadElementDouble(ele,"statDynFricScale","value",true,1.0); // a ratio of static to dynamic friction (=mu_static/mu_dynamic)
		Rho_w=sxml->ReadElementDouble(ele,"rho","value",true,1000);
		ICDfac=sxml->ReadElementDouble(ele,"cdScaleIC","value",true,5); // factor by which to boost drag coefficients during dynamic relaxation IC generation
		ICdt=sxml->ReadElementDouble(ele,"dtIC","value",true,1.0); // convergence analysis time step for IC generation
		ICTmax=sxml->ReadElementDouble(ele,"tmaxIC","value",true,0); // max time for IC generation
		ICthresh=sxml->ReadElementDouble(ele,"threshIC","value",true,0.001); // threshold for relative change in tensions to call it converged
		FreeSurface=sxml->ReadElementDouble(ele,"freesurface","value",true,0.0);// Water free surface
		DtM0=sxml->ReadElementDouble(ele,"dtM","value",true,0.0001);  // default value for desired mooring model time step		
		if(sxml->ExistsElement(ele,"dtM"))DtMauto=sxml->GetAttributeBool(ele->FirstChildElement("dtM"),"auto",true,true);
		WaveKin=sxml->ReadElementInt(ele,"waveKin","value",true,0); // Default=0
		ele=ele->NextSiblingElement();
	}
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void IEnvironment::VisuConfig()const{	
  Log->Print ("Environment properties:");	
  Log->Printf("  DtM..........: %g [s] %s",DtM0,DtMauto?"- Auto":"");
  Log->Printf("  Water Depth..: %g [m]",WtrDpth);
  Log->Printf("  Free Surface.: %g [m]",FreeSurface);
  Log->Printf("  Z Gravity....: %g [m/s^2]",G);
}