//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JComputeMotionRef.cpp \brief Implements the class \ref JComputeMotionRef.

#include "JComputeMotionRef.h"
#include "FunGeo3d.h"

#include <algorithm>
#include <climits>
#include <cfloat>
#include <cstring>

using namespace std;

//##############################################################################
//# JComputeMotionRef
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JComputeMotionRef::JComputeMotionRef(){
  ClassName="JComputeMotionRef";
  PartRidp=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JComputeMotionRef::~JComputeMotionRef(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JComputeMotionRef::Reset(){
  MkBoundFirst=0;
  Mks.clear();
  CaseNbound=CaseNfixed=CaseNmoving=CaseNfloat=0;
  delete[] PartRidp; PartRidp=NULL;
}

//==============================================================================
// Returns index of block or UINT_MAX when it is not found.
//==============================================================================
unsigned JComputeMotionRef::IdxMkBlock(word mkbound)const{
  const unsigned n=CountMkBlock();
  unsigned c=0;
  for(;c<n && Mks[c].mkbound!=mkbound;c++);
  return(c<n? c: UINT_MAX);
}

//==============================================================================
// Adds Mk data to Mks vector.
//==============================================================================
void JComputeMotionRef::AddMkBlock(word mk,word mkbound,unsigned begin,unsigned np){
  if(IdxMkBlock(mkbound)!=UINT_MAX)Run_Exceptioon("Mkbound is already defined.");
  MkBoundFirst=mk-mkbound;
  StMkMotionData v;
  v.mkbound=mkbound;
  v.begin=begin;
  v.np=np;
  v.nid=0;
  for(unsigned c=0;c<3;c++){
    v.id[c]=UINT_MAX;
    v.ps[c]=TDouble3(DBL_MAX);
    v.dis[c]=DBL_MAX;
  }
  Mks.push_back(v);
}

//==============================================================================
// Adds mks data to Mks vector.
//==============================================================================
void JComputeMotionRef::AddMkBlocks(word mkboundfirst
  ,const std::vector<StMkMotionData> &mks)
{
  const unsigned n=unsigned(mks.size());
  for(unsigned c=0;c<n;c++){
    const StMkMotionData &v=mks[c];
    AddMkBlock(v.mkbound+mkboundfirst,v.mkbound,v.begin,v.np);
  }
}

//==============================================================================
// Adds mks[] data to Mks vector.
//==============================================================================
void JComputeMotionRef::AddMkBlocks(word mkboundfirst,unsigned nmk
  ,const StMkMotionData *mks)
{
  for(unsigned c=0;c<nmk;c++){
    const StMkMotionData &v=mks[c];
    AddMkBlock(v.mkbound+mkboundfirst,v.mkbound,v.begin,v.np);
  }
}

//==============================================================================
// Computes RIdp to look for particles by Id.
//==============================================================================
unsigned* JComputeMotionRef::ComputeRidp(unsigned np,const unsigned *idp){
  delete[] PartRidp; PartRidp=NULL;
  PartRidp=new unsigned[CaseNbound];
  memset(PartRidp,255,sizeof(unsigned)*CaseNbound);
  for(unsigned p=0;p<np;p++)if(idp[p]<CaseNbound)PartRidp[idp[p]]=p;
  return(PartRidp);
}

//==============================================================================
// Find reference points to compute motion.
//==============================================================================
void JComputeMotionRef::ComputeRefPoints(unsigned casenfixed,unsigned casenmoving
  ,unsigned casenfloat,unsigned np,const unsigned *idp,const tdouble3 *posd
  ,const tfloat3 *posf,const unsigned *ridp)
{
  if(!CountMkBlock())Run_Exceptioon("No MK data defined.");
  const bool possimple=(posd==NULL);
  //-Computes ridp when it is NULL.
  CaseNfixed=casenfixed;
  CaseNmoving=casenmoving;
  CaseNfloat=casenfloat;
  CaseNbound=CaseNfixed+CaseNmoving+CaseNfloat;
  if(ridp==NULL)ridp=ComputeRidp(np,idp);

  //-Busca puntos de referencia para calcular movimiento
  //-----------------------------------------------------
  const unsigned nmk=unsigned(Mks.size());
  for(unsigned cmk=0;cmk<nmk;cmk++){
    unsigned pini=Mks[cmk].begin;
    unsigned pfin=pini+Mks[cmk].np;
    //-Computes domain size of the body.
    tdouble3 posmin=TDouble3(DBL_MAX),posmax=TDouble3(-DBL_MAX);
    for(unsigned p=pini;p<pfin;p++)if(ridp[p]<UINT_MAX){
      const unsigned rp=ridp[p];
      const tdouble3 ps=(possimple? ToTDouble3(posf[rp]): posd[rp]);
      if(posmin.x>ps.x)posmin.x=ps.x;
      if(posmin.y>ps.y)posmin.y=ps.y;
      if(posmin.z>ps.z)posmin.z=ps.z;
      if(posmax.x<ps.x)posmax.x=ps.x;
      if(posmax.y<ps.y)posmax.y=ps.y;
      if(posmax.z<ps.z)posmax.z=ps.z;
    }
    //printf("Mk[%d] (%g,%g,%g)-(%g,%g,%g)\n",cmk,posmin.x,posmin.y,posmin.z,posmax.x,posmax.y,posmax.z);
    //-Look for the particles closest to the 8 corners.
    tdouble3 poslim[8];
    unsigned idlim[8];
    unsigned nlim=0;
    for(unsigned clim=0;clim<8;clim++){
      tdouble3 pm;
      switch(clim){
      case 0: pm=TDouble3(posmin.x,posmin.y,posmin.z); break;
      case 1: pm=TDouble3(posmax.x,posmin.y,posmin.z); break;
      case 2: pm=TDouble3(posmax.x,posmax.y,posmin.z); break;
      case 3: pm=TDouble3(posmin.x,posmax.y,posmin.z); break;
      case 4: pm=TDouble3(posmin.x,posmin.y,posmax.z); break;
      case 5: pm=TDouble3(posmax.x,posmin.y,posmax.z); break;
      case 6: pm=TDouble3(posmax.x,posmax.y,posmax.z); break;
      case 7: pm=TDouble3(posmin.x,posmax.y,posmax.z); break;
      }
      unsigned id=pini;
      tdouble3 pslim;
      double dis=DBL_MAX;
      for(unsigned p=pini;p<pfin;p++)if(ridp[p]<UINT_MAX){
        const unsigned rp=ridp[p];
        const tdouble3 ps=(possimple? ToTDouble3(posf[rp]): posd[rp]);
        double dis2=fgeo::PointsDist(ps,pm);
        if(dis>dis2){ dis=dis2; id=p; pslim=ps; }
      }
      //-Discards repeated particles.
      bool rep=false;
      for(unsigned cr=0;cr<nlim && !rep;cr++)rep=(idlim[cr]==id);
      if(!rep){
        //printf("idlim[%u]=%u (%g,%g,%g)\n",nlim,id,pos[id].x,pos[id].y,pos[id].z);
        idlim[nlim]=id; poslim[nlim]=pslim; nlim++;
      }
    }
    //-Among the selected ones, look for the one furthest away from the first one.
    //printf("p1[%u]=(%g,%g,%g)\n",idlim[0],pos[idlim[0]].x,pos[idlim[0]].y,pos[idlim[0]].z);
    if(nlim>1){
      tdouble3 ps1=poslim[0];
      unsigned id=1;
      tdouble3 ps2=poslim[1];
      double dis=fgeo::PointsDist(ps1,ps2);
      for(unsigned clim=2;clim<nlim;clim++){
        ps2=poslim[clim];
        double dis2=fgeo::PointsDist(ps1,ps2);
        if(dis<dis2){ dis=dis2; id=clim; }
      }
      unsigned aux=idlim[1];   idlim[1]=idlim[id];   idlim[id]=aux;
      tdouble3 aux2=poslim[1]; poslim[1]=poslim[id]; poslim[id]=aux2;
      //printf("p2[%u]=(%g,%g,%g)\n",idlim[1],pos[idlim[1]].x,pos[idlim[1]].y,pos[idlim[1]].z);
    }
    //-Among the selected ones, look for the one furthest away from the straight line between the two selected points.
    if(nlim>2){
      tdouble3 ps1=poslim[0];
      tdouble3 ps2=poslim[1];
      unsigned id=2;
      tdouble3 ps3=poslim[id];
      double dis=fgeo::LinePointDist(ps1,ps2,ps3);
      for(unsigned clim=3;clim<nlim;clim++){
        ps3=poslim[clim];
        double dis2=fgeo::LinePointDist(ps1,ps2,ps3);
        if(dis<dis2){ dis=dis2; id=clim; }
      }
      unsigned aux=idlim[2];   idlim[2]=idlim[id];   idlim[id]=aux;
      tdouble3 aux2=poslim[2]; poslim[2]=poslim[id]; poslim[id]=aux2;
      //printf("p3[%u]=(%g,%g,%g)\n",idlim[2],pos[idlim[2]].x,pos[idlim[2]].y,pos[idlim[2]].z);
    }
    //-Saves the selected particles and their distance from the first one for periodic conditions.
    Mks[cmk].nid=(nlim>3? 3: nlim);
    for(unsigned clim=0;clim<Mks[cmk].nid;clim++){
      Mks[cmk].id[clim]=idlim[clim];
      Mks[cmk].ps[clim]=poslim[clim];
      Mks[cmk].dis[clim]=fgeo::PointsDist(poslim[0],poslim[clim])*1.1;
      //printf("p(%d)[%u]=(%g,%g,%g)\n",clim,mkparts[cmk].id[clim],mkparts[cmk].ps[clim].x,mkparts[cmk].ps[clim].y,mkparts[cmk].ps[clim].z);
    }
  }
}

//==============================================================================
// Copy motion reference data to mks.
//==============================================================================
void JComputeMotionRef::GetMotionRefData(StMkMotionData &v)const{
  const unsigned idx=IdxMkBlock(v.mkbound);
  if(idx==UINT_MAX)Run_Exceptioon("Mkbound value is was not found.");
  const StMkMotionData &v0=Mks[idx];
  //-Copy data.
  v.nid=v0.nid;
  for(unsigned cv=0;cv<3;cv++){
    v.id [cv]=v0.id [cv];
    v.ps [cv]=v0.ps [cv];
    v.dis[cv]=v0.dis[cv];
  }
}

//==============================================================================
// Copy motion reference data to mks.
//==============================================================================
void JComputeMotionRef::GetMotionRef(std::vector<StMkMotionData> &mks)const{
  const unsigned n=unsigned(mks.size());
  for(unsigned c=0;c<n;c++)GetMotionRefData(mks[c]);
}

//==============================================================================
// Copy motion reference data to mks.
//==============================================================================
void JComputeMotionRef::GetMotionRef(unsigned nmk,StMkMotionData *mks)const{
  for(unsigned c=0;c<nmk;c++)GetMotionRefData(mks[c]);
}
