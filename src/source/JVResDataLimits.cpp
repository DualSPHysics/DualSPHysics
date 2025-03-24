//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2025 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JVResDataLimits.cpp \brief Implements the class \ref JVResDataLimits.

#include "JVResDataLimits.h"
#include "JBinaryData.h"
#include "JCaseVRes.h"
#include "JSpVtkShape.h"
#include "Functions.h"

#include <algorithm>
#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JVResDataLimitsZone
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JVResDataLimitsZone::JVResDataLimitsZone(unsigned id,unsigned parentid
  ,int trackmk,int trackmkb)
  :Id(id),ParentId(parentid),OutLimits(ParentId!=UINT_MAX)
  ,TrackingMk(trackmk),TrackingMkBound(trackmkb)
{
  ClassName="JVResDataLimitsZone";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JVResDataLimitsZone::~JVResDataLimitsZone(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JVResDataLimitsZone::Reset(){
  OutLimit.Reset();
  OutLimit2h.Reset();
  SubZones.clear();
}

//==============================================================================
/// Configures outter limits.
//==============================================================================
void JVResDataLimitsZone::SetOutLimits(const JBoxDef& outlimit
  ,const JBoxDef& outlimit2h)
{
  if(!OutLimits)Run_Exceptioon("Outter limits are invalid for current zone.");
  OutLimit=outlimit;
  OutLimit2h=outlimit2h;
}

//==============================================================================
/// Configures inner limits.
//==============================================================================
void JVResDataLimitsZone::AddInnLimits(unsigned id,const JBoxDef& innlimit
  ,const JBoxDef& innlimit2h)
{
  const StSubZone subzo={id,innlimit,innlimit2h};
  SubZones.push_back(subzo);
}

//##############################################################################
//# JVResDataLimits
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JVResDataLimits::JVResDataLimits(){
  ClassName="JVResDataLimits";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JVResDataLimits::~JVResDataLimits(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JVResDataLimits::Reset(){
  FmtVersion=FmtVersionDef;
  Is2D=false;
  Posy2D=0;
  DirData="";
  for(unsigned c=0;c<Count();c++)delete Zones[c];
  Zones.clear();
}

//==============================================================================
/// Adds zone data.
//==============================================================================
JVResDataLimitsZone* JVResDataLimits::AddZone(unsigned id,unsigned parentid
  ,int trackmk,int trackmkb)
{
  JVResDataLimitsZone* pzo=new JVResDataLimitsZone(id,parentid,trackmk,trackmkb);
  Zones.push_back(pzo);
  return(pzo);
}

//==============================================================================
/// Returns file name and path.
//==============================================================================
std::string JVResDataLimits::GetFileName(std::string dir){
  return(fun::GetDirWithSlash(dir)+"VResData.ibi4");
}

//==============================================================================
/// Loads data from object JCaseVRes.
//==============================================================================
void JVResDataLimits::LoadData(const JCaseVRes* vrdat){
  Reset();
  Is2D  =vrdat->GetIs2D();
  Posy2D=vrdat->GetPosy2D();
  const unsigned nz=vrdat->Count();
  for(unsigned cz=0;cz<nz;cz++){
    const JCaseVRes_Box* zbox=vrdat->GetZoneBox(cz);
    const unsigned parentid=(zbox->Parent? zbox->Parent->Id: UINT_MAX);
    const int trkmk=zbox->GetTrackingMk();
    const int trkmkb=zbox->GetTrackingMkBound();
    JVResDataLimitsZone* zo=new JVResDataLimitsZone(cz,parentid,trkmk,trkmkb);
    //-Defines outter limits.
    if(zo->OutLimits){
      zo->SetOutLimits(zbox->GetPtBox(),zbox->GetBuffBox());
    }
    //-Defines subzones with inner limits.
    const unsigned nsub=zbox->Count();
    for(unsigned csub=0;csub<nsub;csub++){
      const JCaseVRes_Box* subbox=zbox->GetSubZoneBox(csub);
      const unsigned idsub=unsigned(subbox->Id);
      zo->AddInnLimits(idsub,subbox->GetPtBox(),subbox->GetParentBuffEndBox());
    }
    //-Adds new zone.
    Zones.push_back(zo);
  }
}

//==============================================================================
/// Saves JBoxDef data in binary item.
//==============================================================================
void JVResDataLimits::SaveBoxDef(std::string name,std::string subname
  ,const JBoxDef& box,JBinaryData* item)const
{
  item->SetvDouble3(name+"Min"+subname,box.GetPosMin());
  if(box.IsSimple()){
    item->SetvDouble3(name+"Max"+subname,box.GetPosMax());
  }
  else{
    item->SetvDouble3(name+"Vx"+subname,box.GetVx());
    if(!box.GetIs2D())item->SetvDouble3(name+"Vy"+subname,box.GetVy());
    item->SetvDouble3(name+"Vz"+subname,box.GetVz());
  }
}

//==============================================================================
/// Saves binary file VResData.ibi4.
//==============================================================================
void JVResDataLimits::SaveFile(std::string dir){
  DirData=fun::GetDirWithSlash(dir);
  JBinaryData bdat(ClassName);
  bdat.SetvUint("FmtVersion",FmtVersionDef);
  //-Saves general variables.
  bdat.SetvBool("Is2D",Is2D);
  bdat.SetvDouble("Posy2D",Posy2D);
  const unsigned nz=Count();
  bdat.SetvUint("ZonesCount",nz);
  for(unsigned cz=0;cz<nz;cz++){
    const JVResDataLimitsZone* zo=Zones[cz];
    JBinaryData* item=bdat.CreateItem(fun::PrintStr("Zone_%u",cz));
    item->SetvUint("Id",zo->Id);
    //-Saves outter limits.
    if(zo->OutLimits){
      item->SetvUint("ParentId",zo->ParentId);
      if(zo->TrackingMk>0){
        item->SetvInt("TrackingMk",zo->TrackingMk);
        item->SetvInt("TrackingMkBound",zo->TrackingMkBound);
      }
      SaveBoxDef("OutLimit"  ,"",zo->GetOutLimit(),item);
      SaveBoxDef("OutLimit2h","",zo->GetOutLimit2h(),item);
    }
    //-Saves subzones with inner limits.
    const unsigned nsub=zo->SubCount();
    item->SetvUint("SubzonesCount",nsub);
    for(unsigned csub=0;csub<nsub;csub++){
      const string subname=fun::PrintStr("_z%u",csub);
      item->SetvUint(string("Id")+subname,zo->GetInnId(csub));
      SaveBoxDef("InnLimit"  ,subname,zo->GetInnLimit  (csub),item);
      SaveBoxDef("InnLimit2h",subname,zo->GetInnLimit2h(csub),item);
    }
  }
  //-Saves information of particle blocks.
  bdat.SaveFile(DirData+GetFileName(),false,true);
}

//==============================================================================
/// Loads JBoxDef data in binary item.
//==============================================================================
JBoxDef JVResDataLimits::LoadBoxDef(std::string name,std::string subname
  ,const JBinaryData* item)const
{
  JBoxDef box;

  const tdouble3 vmin=item->GetvDouble3(name+"Min"+subname);
  if(item->ExistsValue(name+"Max"+subname)){
    const tdouble3 vmax=item->GetvDouble3(name+"Max"+subname);
    box=JBoxDef(vmin,vmax,Is2D);
  }
  else{
    const tdouble3 vx=item->GetvDouble3(name+"Vx"+subname);
    const tdouble3 vy=item->GetvDouble3(name+"Vy"+subname,Is2D);
    const tdouble3 vz=item->GetvDouble3(name+"Vz"+subname);
    box=JBoxDef(vmin,vx,vy,vz,Is2D);
  }
  return(box);
}

//==============================================================================
/// Updates the bdat object to support the expected format version.
//==============================================================================
void JVResDataLimits::LoadFileFmt240102(JBinaryData* bdat)const{
  const unsigned nz=bdat->GetvUint("ZonesCount");
  bool chk2d=false;
  bool is2d=false;
  double posy2d=0;
  for(unsigned cz=0;cz<nz && !chk2d;cz++){
    const JBinaryData* item=bdat->GetItem(fun::PrintStr("Zone_%u",cz));
    const unsigned parentid=item->GetvUint("ParentId",true,UINT_MAX);
    //-Defines outter limits.
    if(parentid!=UINT_MAX){
      const tdouble3 vmin  =item->GetvDouble3("OutLimitMin");
      const tdouble3 vmax  =item->GetvDouble3("OutLimitMax");
      is2d=(vmin.y==vmax.y);
      if(is2d)posy2d=vmin.y;
      chk2d=true;
    }
  }
  //-Updates bdat object.
  bdat->SetvBool("Is2D",is2d);
  bdat->SetvDouble("Posy2D",posy2d);
}

//==============================================================================
/// Loads binary file VResData.ibi4.
//==============================================================================
void JVResDataLimits::LoadFile(std::string dir){
  Reset();
  DirData=fun::GetDirWithSlash(dir);
  const string file=DirData+GetFileName();
  JBinaryData bdat;
  bdat.LoadFile(file,ClassName);
  FmtVersion=bdat.GetvUint("FmtVersion");
  //-Check file format version.
  const int fmt240102=240102;
  if(FmtVersion!=FmtVersionDef && FmtVersion!=fmt240102)Run_ExceptioonFile(
    fun::PrintStr("Expected format file is %d. Format file %d is unknown."
    ,FmtVersionDef,FmtVersion),file);
  //-Updates the bdat object to support the expected format version.
  if(FmtVersion==fmt240102)LoadFileFmt240102(&bdat);

  //-Loads general variables.
  Is2D=bdat.GetvBool("Is2D");
  Posy2D=bdat.GetvDouble("Posy2D");
  const unsigned nz=bdat.GetvUint("ZonesCount");
  for(unsigned cz=0;cz<nz;cz++){
    const JBinaryData* item=bdat.GetItem(fun::PrintStr("Zone_%u",cz));
    const unsigned id=item->GetvUint("Id");
    const unsigned parentid=item->GetvUint("ParentId",true,UINT_MAX);
    const int trkmk=item->GetvInt("TrackingMk",true,-1);
    const int trkmkb=(trkmk<=0? -1: item->GetvInt("TrackingMkBound",true,-1));
    JVResDataLimitsZone* zo=new JVResDataLimitsZone(id,parentid,trkmk,trkmkb);
    //-Defines outter limits.
    if(zo->OutLimits){
      zo->SetOutLimits(LoadBoxDef("OutLimit","",item),
                       LoadBoxDef("OutLimit2h","",item));
    }
    //-Defines subzones with inner limits.
    const unsigned nsub=item->GetvUint("SubzonesCount");
    for(unsigned csub=0;csub<nsub;csub++){
      const string subname=fun::PrintStr("_z%u",csub);
      const unsigned idsub =item->GetvUint(string("Id")+subname);
      zo->AddInnLimits(idsub,LoadBoxDef("InnLimit",subname,item),
                             LoadBoxDef("InnLimit2h",subname,item));
    }
    //-Adds new zone.
    Zones.push_back(zo);
  }
}

//==============================================================================
// Returns the requested zone.
//==============================================================================
const JVResDataLimitsZone* JVResDataLimits::GetZone(unsigned id)const{
  if(id>=Count())Run_Exceptioon(fun::PrintStr("The requested zone (%u) does not exist.",id));
  return((const JVResDataLimitsZone*)Zones[id]);
}

//==============================================================================
/// Load vertices of 2-D quad from JBoxDef object.
//==============================================================================
void JVResDataLimits::GetBoxPoints2d(JBoxDef box,tdouble3* vpt){
  tdouble3 points[8];
  box.GetBoxPoints(points);
  vpt[0]=points[0];
  vpt[1]=points[1];
  vpt[2]=points[5];
  vpt[3]=points[4];
}

//==============================================================================
/// Load vertice and vectors of 3-D box from JBoxDef object.
//==============================================================================
void JVResDataLimits::GetBoxPoints3d(JBoxDef box,tdouble3* ptvec){
  ptvec[0]=box.GetPosMin();
  ptvec[1]=box.GetVx();
  ptvec[2]=box.GetVy();
  ptvec[3]=box.GetVz();
}

//==============================================================================
// Saves VTK file with limits of zones.
//==============================================================================
void JVResDataLimits::SaveVtkLimits(std::string fname,bool onefile)const{
  string file=fun::GetWithoutExtension(fname);
  const unsigned nz=Count();
  if(nz<2)Run_ExceptioonFile("No subdomains are available.",file+".vtk");
  JSpVtkShape ss_1;
  if(Is2D){
    for(unsigned id=0;id<nz;id++){
      const word wid=word(id);
      JSpVtkShape ss_2;
      JSpVtkShape& ss=(onefile? ss_1: ss_2);
      const JVResDataLimitsZone* pzone=GetZone(id);
      //-Draws outter limits of zone.
      if(pzone->OutLimits){
        tdouble3 pt[4],pt2[4];
        GetBoxPoints2d(pzone->GetOutLimit(),pt);
        GetBoxPoints2d(pzone->GetOutLimit2h(),pt2);
        ss.AddQuadWire(pt[0],pt[1],pt[2],pt[3],wid);
        ss.AddQuadWire(pt2[0],pt2[1],pt2[2],pt2[3],wid);
        for(unsigned c=0;c<4;c++){
          ss.AddLine(pt[c],pt2[c],wid);
        }
      }
      //-Draws inner limits of zone.
      const unsigned nsub=pzone->SubCount();
      for(unsigned csub=0;csub<nsub;csub++){
        tdouble3 pt[4],pt2[4];
        GetBoxPoints2d(pzone->GetInnLimit(csub),pt);
        GetBoxPoints2d(pzone->GetInnLimit2h(csub),pt2);
        ss.AddQuadWire(pt[0],pt[1],pt[2],pt[3],wid);
        ss.AddQuadWire(pt2[0],pt2[1],pt2[2],pt2[3],wid);
        for(unsigned c=0;c<4;c++){
          ss.AddLine(pt[c],pt2[c],wid);
        }
      }
      if(!onefile)ss_2.SaveVtk(file+fun::PrintStr("%02d.vtk",id),"Zone");
    }
    if(onefile)ss_1.SaveVtk(file+".vtk","Zone");
  }
  else{
    for(unsigned id=0;id<nz;id++){
      const word wid=word(id);
      JSpVtkShape ss_2;
      JSpVtkShape& ss=(onefile? ss_1: ss_2);
      const JVResDataLimitsZone* pzone=GetZone(id);
      //-Draws outter limits of zone.
      if(pzone->OutLimits){
        tdouble3 ptvec[4];
        GetBoxPoints3d(pzone->GetOutLimit(),ptvec);
        ss.AddBoxSizeVecWire(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
        tdouble3 ptvec2[4];
        GetBoxPoints3d(pzone->GetOutLimit2h(),ptvec2);
        ss.AddBoxSizeVecWire(ptvec2[0],ptvec2[1],ptvec2[2],ptvec2[3],wid);
        //for(unsigned c=0;c<4;c++){
        //  const tdouble3 sx=TDouble3(ptb.x-pta.x,0,0);
        //  const tdouble3 sy=TDouble3(0,ptb.y-pta.y,0);
        //  const tdouble3 sz=TDouble3(0,0,ptb.z-pta.z);
        //  const tdouble3 sx2=TDouble3(pt2b.x-pt2a.x,0,0);
        //  const tdouble3 sy2=TDouble3(0,pt2b.y-pt2a.y,0);
        //  const tdouble3 sz2=TDouble3(0,0,pt2b.z-pt2a.z);
        //  ss.AddLine(pta         ,pt2a            ,wid);
        //  ss.AddLine(pta+sx      ,pt2a+sx2        ,wid);
        //  ss.AddLine(pta+sy      ,pt2a+sy2        ,wid);
        //  ss.AddLine(pta+sx+sy   ,pt2a+sx2+sy2    ,wid);
        //  ss.AddLine(pta+sz      ,pt2a+sz2        ,wid);
        //  ss.AddLine(pta+sz+sx   ,pt2a+sz2+sx2    ,wid);
        //  ss.AddLine(pta+sz+sy   ,pt2a+sz2+sy2    ,wid);
        //  ss.AddLine(pta+sz+sx+sy,pt2a+sz2+sx2+sy2,wid);
        //}
      }
      //-Draws inner limits of zone.
      const unsigned nsub=pzone->SubCount();
      for(unsigned csub=0;csub<nsub;csub++){
        tdouble3 ptvec[4];
        GetBoxPoints3d(pzone->GetInnLimit(csub),ptvec);
        ss.AddBoxSizeVecWire(ptvec[0],ptvec[1],ptvec[2],ptvec[3],wid);
        tdouble3 ptvec2[4];
        GetBoxPoints3d(pzone->GetInnLimit2h(csub),ptvec2);
        ss.AddBoxSizeVecWire(ptvec2[0],ptvec2[1],ptvec2[2],ptvec2[3],wid);
        //for(unsigned c=0;c<4;c++){
        //  const tdouble3 sx=TDouble3(ptb.x-pta.x,0,0);
        //  const tdouble3 sy=TDouble3(0,ptb.y-pta.y,0);
        //  const tdouble3 sz=TDouble3(0,0,ptb.z-pta.z);
        //  const tdouble3 sx2=TDouble3(pt2b.x-pt2a.x,0,0);
        //  const tdouble3 sy2=TDouble3(0,pt2b.y-pt2a.y,0);
        //  const tdouble3 sz2=TDouble3(0,0,pt2b.z-pt2a.z);
        //  ss.AddLine(pta         ,pt2a            ,wid);
        //  ss.AddLine(pta+sx      ,pt2a+sx2        ,wid);
        //  ss.AddLine(pta+sy      ,pt2a+sy2        ,wid);
        //  ss.AddLine(pta+sx+sy   ,pt2a+sx2+sy2    ,wid);
        //  ss.AddLine(pta+sz      ,pt2a+sz2        ,wid);
        //  ss.AddLine(pta+sz+sx   ,pt2a+sz2+sx2    ,wid);
        //  ss.AddLine(pta+sz+sy   ,pt2a+sz2+sy2    ,wid);
        //  ss.AddLine(pta+sz+sx+sy,pt2a+sz2+sx2+sy2,wid);
        //}
      }
      if(!onefile)ss_2.SaveVtk(file+fun::PrintStr("%02d.vtk",id),"Zone");
    }
    ss_1.SaveVtk(file+".vtk","Zone");
  }
}

