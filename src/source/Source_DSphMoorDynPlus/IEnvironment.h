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

///<\file IEnvironment.h \brief Defines the class \ref IEnvironment.

#ifndef _Environment_
#define _Environment_

#include "JXml.h"
#include "TypesDef.h"
#include "JObject.h"
#include "JLog2.h"
#include <algorithm>
#include <cmath>
#include <iostream>

class IEnvironment : protected JObject
{
private:
	JLog2* Log;
	bool DtMmod;             ///<True if the dtM0 was adjusted.
	bool DtMauto;            ///<True if the dtM0 is adjusted automatically.
	double TimeMax;          ///<Time of simulation
	double G;                ///<Gravity (m/s^2)
	double WtrDpth;          ///<Depth of the water (m)
	double Rho_w;            ///<Density
	double Kb;               ///<Bottom stiffness (Pa/m)
	double Cb;               ///<Bottom damping   (Pa/m/s)
	int WaveKin;             ///<Wave kinematics flag (0=off,gt than 0=on)
	int WriteUnits;          ///<Global switch for whether to show the units line in the output files (1,default),or skip it (0)
	double FricCoeff;        ///<General bottom friction coefficient,as a start
	double FricDamp;         ///<Damping coefficient used to model the friction at speeds near zero
	double StatDynFricScale; ///<Ratio of static to dynamic friction (=mu_static/mu_dynamic)
	double ICDfac;           ///<Factor by which to boost drag coefficients during dynamic relaxation IC generation
	double ICdt;             ///<Convergence analysis time step for IC generation
	double ICTmax;           ///<Max time for IC generation
	double ICthresh;         ///<Threshold for relative change in tensions to call it converged
	double DtM0;             ///<Default value for desired mooring model time step
	double FreeSurface;      ///<Z position of water free surface	
	
	void ReadXml(JXml* sxml,TiXmlElement* lis); ///<Reads the Xml file

public:
	IEnvironment(JLog2* log); 
	virtual ~IEnvironment(); 
	
	void Reset();
	void LoadXml(JXml* sxml,const std::string& place);
  void VisuConfig()const; 

	bool   GetDtMmod()   const { return DtMmod;          }; ///<Returns if the dtM value was modified.
	bool   GetDtMauto()  const { return DtMauto;         }; ///<Returns if the dtM value is calculated automatically.
	double GetTimeMax()  const { return TimeMax;         }; ///<Time of simulation
	double GetG()        const { return G;               }; ///<Returns the Gravity		
	double GetWaterDepth()const { return WtrDpth;        }; ///<Returns the depth of the water
	double GetRho_w()    const { return Rho_w;           }; ///<Returns the density	
	double GetKb()       const { return Kb;              }; ///<Returns the bottom stiffness  		
	double GetCb()       const { return Cb;              }; ///<Returns the botom damping     			
	int GetWaveKin()     const { return WaveKin;         }; ///<Returns the waves kinematics flag	
	int GetWriteUnits()  const { return WriteUnits;      }; ///<Returns the Global switch for whether to show the units		
	double GetFricCoeff()const { return FricCoeff;       }; ///<Returns the General bottom friction coefficient
	double GetFricDamp() const { return FricDamp;        }; ///<Returns the Damping coefficient 
	double GetICDfac()   const { return ICDfac;          }; ///<Returns the factor by which to boost drag coefficients during dynamic relaxation IC generation
	double GetICdt()     const { return ICdt;            }; ///<Returns the Convergence analysis time step for IC generation
	double GetICTmax()   const { return ICTmax;          }; ///<Returns the Max time for IC generation
	double GetICthresh() const { return ICthresh;        }; ///<Returns the Threshold for relative change in tensions to call it converged
	double GetDtM0()     const { return DtM0;            }; ///<Returns the Default value for desired mooring model time step	
	double GetFreeSurface()const { return FreeSurface;   }; ///<Returns the Z position of water free surface	
	double GetStatDynFricScale()const { return StatDynFricScale; }; ///<Returns the Ratio of static to dynamic friction 

	void SetDtM0       (const double dtm){ DtM0=dtm; DtMmod=true;}; ///<Sets a new value of Dtm0
  void SetG          (const tfloat3 g) { G=std::abs(g.z);};       ///<Sets the Gravity		
  void SetTimeMax    (const double t)  { TimeMax=t;      };       ///<Sets Time of simulation	
	void SetFreeSurface(const double fs) { FreeSurface=fs; };       ///<Sets the Z position of water free surface
};


#endif