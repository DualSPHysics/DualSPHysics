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

/// \file IBody.cpp \brief Implements the class \ref IBody.

#include "IBody.h"

//==============================================================================
/// Constructor.
//==============================================================================

IBody::IBody(JLog2* log):Log(log) {
	Reset();
	ClassName="IBody";
}

//==============================================================================
/// Destructor.
//==============================================================================
IBody::~IBody() {
	Reset();
	Lines.clear();
	Log=NULL;
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void IBody::Reset() {
	LineCount=0;
	FairCount=0;
	Env=NULL;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void IBody::LoadXml(JXml* sxml,const std::string& place,TiXmlElement* eleb,const unsigned numBody) {
	std::string function="LoadXml";
	TiXmlNode* node=sxml->GetNode(place,false);
	if(!node)Run_Exceptioon("Cannot find the element \'"+place+"\'.");
	if(!eleb)Run_Exceptioon("Cannot find the element \'"+place+"."+eleb->Value()+"\'.");	
	ReadXml(sxml,place,eleb,numBody);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void IBody::ReadXml(JXml* sxml,const std::string& place_in,TiXmlElement* eleBody,const unsigned numBody) {
	//-Loads mooring zones.
	//PRINT: look the Number of Mooting tag
	sxml->CheckElementNames(eleBody,false,"depth");

	if(Debug) Log->PrintfDbg("numBody - IBody: %d\n",numBody);

	if(eleBody) { 
		Idb=numBody; //ID of mooring
		Ref=sxml->GetAttributeInt(eleBody,"ref",true,Idb); //If "ref" doesn't exist,stores Number value of creation
		Depth=sxml->ReadElementDouble(eleBody,"depth","value",true,DBL_MAX); //If depth doesn't exist, stores more later the WaterDepth of SolverOptions
	}
}

//==============================================================================
/// Returns a pointer of connections of this IBody for type
//==============================================================================
ILine* IBody::GetLine(unsigned idl)const {
	ILine* line=NULL;
	for(unsigned l=0; l<LineCount && !line; l++)if(Lines[l]->GetId()==idl)line=Lines[l];
	return line;
}

//==========================================================================================
/// Makes the setup for this IBody. Initializes the Outptfile if exists this option in XML
//==========================================================================================
void IBody::Setup(std::string &dir,IEnvironment* env_in,std::vector<ILine*> lines_in) {
	std::string function="Setup";
	Env=env_in;
	Lines=lines_in;
	LineCount=unsigned(Lines.size());
	//-Stores the fairs connected to this body
	for(unsigned l=0; l<LineCount; l++) {
		ILine* line=Lines[l];
		if(line->GetUseFair()){
			std::vector<IConnection*> fairs=line->GetFairs(Ref);
			for(unsigned c=0; c<fairs.size(); c++){Fairs.push_back(fairs[c]);}
		}
	}
	//-Stores the number of fairs
	FairCount=unsigned(Fairs.size());
}

//============================================================================================
/// Visualises the properties of the line
//============================================================================================
void IBody::VisuProperties()const{
	printf("Body=%u:\n",Ref);
	printf("  Idb......:%u:\n",Idb);
	printf("  FairCount:%u:\n",FairCount);
	printf("  LineCount:%u:\n",LineCount);
	printf("  Depth....:%g:\n",Depth);
}


