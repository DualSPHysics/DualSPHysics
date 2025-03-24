//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphVResZone.h \brief Declares the class \ref JSphVResZone.

#ifndef _JSphVResZoneh_
#define _JSphVResZoneh_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#include "JMatrix4.h"

class JXml;
class TiXmlElement;
class JLog2;

class JSphVResZone : protected JObject
{
private:
	const bool Cpu;
	const StCteSph CSP; ///< Structure with main SPH constants values and configurations.
	JLog2 *Log;

	tdouble3 BoxLimitMin;
	tdouble3 BoxLimitMax;
	JMatrix4d Mat;
	bool IsSimple;

	tuint3 NPoints;
	int Ntot;

	bool TrackingisActive=false;
	unsigned TrackingMk=0;

	std::vector<tdouble3> Points;
	std::vector<tdouble3> Normals;

  double    Width;
	tdouble3  BoxLimitMinInner;
	tdouble3  BoxLimitMaxInner;
	tdouble3  BoxLimitMinOuter;
	tdouble3  BoxLimitMaxOuter;
	tdouble3  BoxLimitMinMid;
	tdouble3  BoxLimitMaxMid;
	tdouble3  BoxSize;
	unsigned  Zone;
	tdouble3  Origin;
	bool      Inner;

	tdouble3 MapRealPosMin;
	tdouble3 MapRealPosMax;

	
	void Reset();
	void ComputeInterface();


public:
	JSphVResZone(bool cpu, const StCteSph &csp,bool inner,unsigned zone
    ,tdouble3 boxlimitmininner,tdouble3 boxlimitmaxinner
    ,tdouble3 boxlimitminouter,tdouble3 boxlimitmaxouter
    ,tdouble3 boxlimitminmid,tdouble3 boxlimitmaxmid
		,tdouble3 maprealposmin,tdouble3 maprealposmax
		,bool trackingisactive,unsigned trackingmk,JMatrix4d mat,bool issimple);
  ~JSphVResZone();


	void Config();
	bool InZoneBoxOuter(const tdouble3 &ps) const;
	bool InZoneBoxInner(const tdouble3 &ps) const;
	bool InZoneBox(const tdouble3 &ps) const;
	bool InZoneBoxMid(const tdouble3 &ps) const;
	bool is_Out(const tdouble3 &ps) const;
	bool is_Normal(const tdouble3 &ps) const;
	int getCount() const { return Ntot; };
	
	unsigned getZone()const {return Zone; };
	
	tdouble3 *getPoints(){return Points.data();};
	tdouble3 *getNormals(){return Normals.data();};

	
	tdouble3 getNormal(const tdouble3 &ps);
	double getmass(){return CSP.massfluid;};

  tdouble3    GetBoxLimitMinInner(){return BoxLimitMinInner;};
  tdouble3    GetBoxLimitMaxInner(){return BoxLimitMaxInner;};
  tdouble3    GetBoxLimitMinOuter(){return BoxLimitMinOuter;};
  tdouble3    GetBoxLimitMaxOuter(){return BoxLimitMaxOuter;};
  tdouble3    GetBoxLimitMinMid(){return BoxLimitMinMid;};
  tdouble3    GetBoxLimitMaxMid(){return BoxLimitMaxMid;};

	tdouble3    GetBoxLimitMin(){return BoxLimitMinInner;};
	tdouble3    GetBoxLimitMax(){return BoxLimitMaxInner;};
	tdouble3    GetBoxDomMin(){return (Inner? BoxLimitMinOuter: BoxLimitMinInner);};
	tdouble3    GetBoxDomMax(){return (Inner? BoxLimitMaxOuter: BoxLimitMaxInner);};

	float 		  GetWidth(){return float(Width);};
	bool 		    GetInner(){return Inner;};
	bool        TrackingIsActive()const{ return(TrackingisActive); }
  unsigned    GetTrackingMk()const{ return(TrackingMk); }
	JMatrix4d   GetMat(){return Mat;};

};

#endif