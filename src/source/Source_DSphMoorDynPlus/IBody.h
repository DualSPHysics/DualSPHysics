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

/// \file IBody.h \brief Defines the class \ref IBody.

#ifndef _Body_
#define _Body_

#include "JXml.h"
#include "ILine.h"
#include "IEnvironment.h"
#include "IConnection.h"
#include "JObject.h"
#include "TypesMoorDynPlus.h"

class IBody : protected JObject
{
private:	
	JLog2* Log;
	unsigned Ref;                  ///<Identifier of body introduced in xml file
	unsigned Idb;                  ///<Identifier of body. Secuential by creation [0, 1, 2...]
	unsigned FairCount;            ///<Number of fairleads
	unsigned LineCount;            ///<Number of lines of this IBody
	double Depth;                  ///<Stores the Depth of the water (m)

	std::vector<ILine*>Lines;       ///<Vector of pointers to lines of this IBody
	std::vector<IConnection*>Fairs; ///<Vector of pointers to Fairs of its lines
	IEnvironment* Env;              ///<Pointer to environmental settings

	void Reset(); 
	void ReadXml(JXml* sxml
    ,const std::string& place,TiXmlElement* eleb
		,const unsigned num_tag); 

public:
	IBody(JLog2* log); /// Constructor
	virtual ~IBody(); /// Destructor
	
	void LoadXml(JXml* sxml,const std::string& place
	  ,TiXmlElement* eleb,const unsigned num_tag);
	void Setup(std::string& dir,IEnvironment* env_in
	  ,std::vector<ILine*> lines); 

	ILine* GetLine(unsigned idl)const;

	std::vector<ILine*>       GetLines()const { return Lines; } ///<Returns a vector of Lines
	std::vector<IConnection*> GetFairs()const { return Fairs; } ///<Returns a pointer of connections of this body

	unsigned  GetNFairs()const{ return FairCount; } ///<Returns the number of Fairleads
	unsigned  GetId()const    { return Idb;       } ///<Returns the number of IBody
	unsigned  GetNLines()const{ return LineCount; } ///<Returns the number of IBody
	unsigned  GetRef()const   { return Ref;       } ///<Returns the floating id of IBody
	unsigned  GetFtMk()const  { return GetRef();  } ///<Returns the mk value

	unsigned* GetPtrRef()     { return &Ref;      } ///<Returns a pointer to floating id of IBody
	double    GetDepth()const { return Depth;     } ///<Returns a pointer of the depth
	
	void SetId(const unsigned idb) { Idb=idb;     } ///<Stores the new number
	void SetDepth (double value) { Depth =value;} ///<Sets depth value

	void VisuProperties()const;
};

#endif