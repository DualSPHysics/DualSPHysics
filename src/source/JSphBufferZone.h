/*
 * JSphBufferZone.h
 *
 *  Created on: Nov 9, 2021
 *      Author: francesco
 */

#ifndef JSPHBUFFERZONE_H_
#define JSPHBUFFERZONE_H_

#include <string>
#include <vector>
#include "JObject.h"
#include "DualSphDef.h"
#include "JMatrix4.h"

typedef struct
{
	tdouble3 pointmin, planesize, normal;
} bufferplane;

inline bufferplane Bufferplane(tdouble3 pointmin, tdouble3 planesize, tdouble3 normal)
{
	bufferplane d = {pointmin,
					 planesize,
					 normal};
	return (d);
}

class JXml;
class TiXmlElement;
class JLog2;

class JSphBufferZone : protected JObject
{
private:
	const bool Cpu;
	const StCteSph CSP; ///< Structure with main SPH constants values and configurations.
	JLog2 *Log;

	tdouble3 BoxLimitMin;
	tdouble3 BoxLimitMax;
	JMatrix4d Mat;
	bool IsSimple;

	unsigned Nx = 0;
	unsigned Nz = 0;
	tuint3 NPoints;
	int Ntot;
	tdouble3 *Points = nullptr;
	tdouble3 *Normals = nullptr;

	bool TrackingisActive=false;
	unsigned TrackingMk=0;

	std::vector<tdouble3> PointsIn;
	std::vector<tdouble3> NormalIn;

	
	void CalculateNpoints();
	void AllocateMemory();
	void CalculatePointsNormals();
	void Reset();
	void ComputeInterface();


public:
	JSphBufferZone(bool cpu, const StCteSph &csp, bool inner, unsigned zone, tdouble3 boxlimitmininner, tdouble3 boxlimitmaxinner,
                tdouble3 boxlimitminouter, tdouble3 boxlimitmaxouter, tdouble3 boxlimitminmid, tdouble3 boxlimitmaxmid
								,bool trackingisactive,unsigned trackingmk,JMatrix4d mat,bool issimple);
	void Config();
	bool InZoneBoxOuter(const tdouble3 &ps) const;
	bool InZoneBoxInner(const tdouble3 &ps) const;
	bool InZoneBox(const tdouble3 &ps) const;
	bool InZoneBoxMid(const tdouble3 &ps) const;
	bool is_Out(const tdouble3 &ps) const;
	bool is_Normal(const tdouble3 &ps) const;
	int getCount() const { return Ntot; };
	unsigned getZone() const { return Zone; };
	tdouble3 *getPoints() { return Points; } ;
	tdouble3 *getNormals() { return Normals; };
	unsigned Partout = 0;
	unsigned PartNew = 0;
	unsigned PartFluidToBuffer = 0;
	unsigned PartBufferToFluid = 0;

	double Width;
	tdouble3 BoxLimitMinInner;
	tdouble3 BoxLimitMaxInner;
	tdouble3 BoxLimitMinOuter;
	tdouble3 BoxLimitMaxOuter;
	tdouble3 BoxLimitMinMid;
	tdouble3 BoxLimitMaxMid;
	tdouble3 BoxSize;
	unsigned Zone;
	tdouble3 Origin;
	bool Inner;
	tdouble3 getNormal(const tdouble3 &ps);
	double getmass(){return CSP.massfluid;};


	tdouble3 getBoxLimitMin(){return BoxLimitMinInner;};
	tdouble3 getBoxLimitMax(){return BoxLimitMaxInner;};
	tdouble3 getBoxDomMin(){return (Inner? BoxLimitMinOuter: BoxLimitMinInner);};
	tdouble3 getBoxDomMax(){return (Inner? BoxLimitMaxOuter: BoxLimitMaxInner);};
	float 	getWidth(){return float(Width);};
	bool 		getInner(){return Inner;};
	bool TrackingIsActive()const{ return(TrackingisActive); }
  unsigned  GetTrackingMk()const{ return(TrackingMk); }
	JMatrix4d GetMat(){return Mat;};

};

#endif /* JSPHBUFFERZONE_H_ */
