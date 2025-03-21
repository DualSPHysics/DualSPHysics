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

/// \file IQSlines.h \brief Defines the class \ref IQSlines.


#ifndef _QSlines_
#define _QSlines_

#include <vector>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <algorithm>

#include "JObject.h"
#include "JLog2.h"

using namespace std;

class IQSlines: protected JObject
{
public:
	IQSlines(JLog2* log):Log(log){} /// Constructor
	virtual ~IQSlines(){};
	int Catenary(double XF,double ZF,double L,double EA,double W,double CB,double Tol
		,double* HFout,double* VFout,double* HAout,double* VAout,int Nnodes,std::vector<double>& s
		,std::vector<double>& X,std::vector<double>& Z,std::vector<double>& Te,const unsigned lineid); /// This function return the positions and tensions of a single mooring line
private:
	static const int longwinded=0;	///<Switch to turn on excessive output for locating crashes
	JLog2* Log; ///<Shows and stores messages
};
#endif //!IQSlines
