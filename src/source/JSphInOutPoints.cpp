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

/// \file JSphInOutPoints.cpp \brief Implements the class \ref JSphInOutPoints.

#include "JSphInOutPoints.h"
#include "JSphMk.h"
#include "JDsPartsInit.h"
#include "JXml.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JMatrix4.h"
#include "JDataArrays.h"
#include "JSpVtkData.h"
#include "JSpVtkShape.h"
 
#include <cfloat>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSphInOutPoints
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutPoints::JSphInOutPoints(bool simulate2d,double simulate2dposy,byte layers
  ,double dp,double initialmove,tdouble3 posmin,tdouble3 posmax)
  :Log(AppInfo.LogPtr()),Simulate2D(simulate2d),Simulate2DPosY(simulate2dposy),Layers(layers)
  ,Dp(dp),InitialMove(initialmove),MapRealPosMin(posmin),MapRealPosMax(posmax)
{
  ClassName="JSphInOutPoints";
  Points=NULL;
  PointsInit=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOutPoints::~JSphInOutPoints(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInOutPoints::Reset(){
  ConfigInfo.clear();
  XmlShape="";
  CircleRadius=0;
  Direction=TDouble3(0);
  ResetPoints();
  for(unsigned c=0;c<10;c++)PtDom[c]=TDouble3(DBL_MAX);
  ZonePosMin=ZonePosMax=TDouble3(0);
}

//==============================================================================
/// Frees memory of Points.
//==============================================================================
void JSphInOutPoints::ResetPoints(){
  delete[] Points; Points=NULL;
  delete[] PointsInit; PointsInit=NULL;
  Size=Count=0;
}

//==============================================================================
/// Increases the allocated memory for points.
//==============================================================================
void JSphInOutPoints::ResizeMemory(unsigned newnpt){
  if(Size-Count<newnpt){
    try{
      Points    =fun::ResizeAlloc(Points    ,Count,Count+newnpt);
      PointsInit=fun::ResizeAlloc(PointsInit,0    ,Count+newnpt);
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
    Size=Count+newnpt;
  }
}

//==============================================================================
/// Returns matrix for Y-rotation in 2D.
//==============================================================================
JMatrix4d JSphInOutPoints::ReadRotate2D(const JXml* sxml,TiXmlElement* ele
  ,const tdouble3& pt,std::string& rotationinfo)
{
  sxml->CheckAttributeNames(ele,"rotate","angle centerx centerz anglesunits");
  JMatrix4d m;
  double rotate=sxml->ReadElementDouble(ele,"rotate","angle",true);
  if(rotate){
    double centerx=pt.x;
    double centerz=pt.z;
    const TiXmlElement* erot=sxml->GetFirstElement(ele,"rotate");
    if(sxml->ExistsAttribute(erot,"centerx") || sxml->ExistsAttribute(erot,"centerz")){
      centerx=sxml->ReadElementDouble(ele,"rotate","centerx");
      centerz=sxml->ReadElementDouble(ele,"rotate","centerz");
    }
    string angunits=fun::StrLower(sxml->ReadElementStr(ele,"rotate","anglesunits"));
    if(angunits=="radians")rotate=rotate*TODEG;
    else if(angunits!="degrees")sxml->ErrReadElement(ele,"rotate",false,"The value anglesunits must be \"degrees\" or \"radians\"."); 
    m=JMatrix4d::MatrixRot(rotate,TDouble3(centerx,0,centerz),TDouble3(centerx,1,centerz));
    //-Return config information about rotation.
    rotationinfo=fun::PrintStr("  Rotated: angle:%g centerx:%g centerz:%g",rotate,centerx,centerz);
  }
  return(m);
}

//==============================================================================
/// Returns direction vector starting from direction in text.
/// Returns (0,0,0) for custom direction and (DBL_MAX,...) for invalid direction.
//==============================================================================
tdouble3 JSphInOutPoints::DirectionFromStr(const std::string& strdir)const{
  bool custom=false;
  tdouble3 dir=TDouble3(0);
  if     (strdir=="top"   )dir.z= 1;
  else if(strdir=="bottom")dir.z=-1;
  else if(strdir=="right" )dir.x= 1;
  else if(strdir=="left"  )dir.x=-1;
  else if(strdir=="back"  )dir.y= 1;
  else if(strdir=="front" )dir.y=-1;
  else if(strdir=="custom")custom=true;
  if(dir==TDouble3(0) && !custom)dir=TDouble3(DBL_MAX);    
  return(dir);
}

//==============================================================================
/// Returns direction vector starting from direction in text.
//==============================================================================
std::string JSphInOutPoints::CheckParticlesDirection(const JSphMkBlock* pmk
  ,const tdouble3& dir)const
{
  const tdouble3 size=pmk->GetPosMax()-pmk->GetPosMin();
  std::string error;
  const unsigned ndir=(dir.z!=0? 1: 0)+(dir.y!=0? 1: 0)+(dir.x!=0? 1: 0);
  if(!ndir)error="The direction is null.";
  else if(ndir>1)error="The direction is not parallel to the axes.";
  else if(dir.y!=0 && Simulate2D)error="The direction in Y is invalid for simulations 2D.";
  else if(dir.z!=0 && size.z!=0)error="Domain size of particles in direction Z is not zero.";
  else if(dir.y!=0 && size.y!=0)error="Domain size of particles in direction Y is not zero.";
  else if(dir.x!=0 && size.x!=0)error="Domain size of particles in direction X is not zero.";
  return(error);
}

//==============================================================================
/// Returns direction according to the particles for 3D inlet.
//==============================================================================
tdouble3 JSphInOutPoints::ComputeDir3dFromParts(unsigned mkfluid,
  const JDsPartsInit* partsdata)const
{
  tdouble3 dir=TDouble3(0);
  if(partsdata && partsdata->GetNp()>2){
    const JSphMk* mkinfo=partsdata->GetMkInfo();
    const unsigned cmk=mkinfo->GetMkBlockByMkFluid(mkfluid);
    const JSphMkBlock* pmk=(cmk<mkinfo->Size()? mkinfo->Mkblock(cmk): NULL);
    if(pmk && pmk->Count>2){
      vector<tdouble3> vpoints;
      tdouble3 pm=TDouble3(0);
      {
        const typecode rcode=pmk->Code;
        const unsigned np=partsdata->GetNp();
        const typecode* code=partsdata->GetCode();
        const tdouble3* pos=partsdata->GetPos();
        for(unsigned p=0;p<np;p++)if(code[p]==rcode){
          vpoints.push_back(pos[p]);
          pm=pm+pos[p];
        }
      }
      const unsigned np=unsigned(vpoints.size());
      //-Looks for the furthest point from the midpoint.
      pm=pm/np;
      double dist2=0;
      unsigned cp=0;
      for(unsigned p=0;p<np;p++){
        const double d2=fgeo::PointsDist2(pm,vpoints[p]);
        if(dist2<d2){
          dist2=d2;
          cp=p;
        }
      }
      //-Looks for the furthest point from pt1.
      const tdouble3 pt1=vpoints[cp];
      dist2=0;
      cp=0;
      for(unsigned p=0;p<np;p++){
        const double d2=fgeo::PointsDist2(pt1,vpoints[p]);
        if(dist2<d2){
          dist2=d2;
          cp=p;
        }
      }
      //-Looks for the furthest point from pt1 and pt2.
      const tdouble3 pt2=vpoints[cp];
      dist2=0;
      cp=0;
      const double PQx=pt2.x-pt1.x;
      const double PQy=pt2.y-pt1.y;
      const double PQz=pt2.z-pt1.z;
      for(unsigned p=0;p<np;p++){
        const tdouble3 p3=vpoints[p];
        const double PRx=p3.x-pt1.x;
        const double PRy=p3.y-pt1.y;
        const double PRz=p3.z-pt1.z;
        const double Vi=PQy*PRz-PRy*PQz;
        const double Vj=-(PQx*PRz-PRx*PQz);
        const double Vk=PQx*PRy-PRx*PQy;
        const double area2=(Vi*Vi+Vj*Vj+Vk*Vk);
        if(dist2<area2){
          dist2=area2;
          cp=p;
        }
      }
      const tdouble3 pt3=vpoints[cp];
      //-Computes normal from selected points.
      dir=fgeo::VecUnitary(fgeo::ProductVec(pt1-pt2,pt2-pt3));
    }
  }
  return(dir);
}

//==============================================================================
/// Creates points starting from special fluid particles for 2D and 3D inlet.
//==============================================================================
void JSphInOutPoints::Create2d3d_Particles(const JXml* sxml,TiXmlElement* ele
  ,const JDsPartsInit* partsdata)
{
  const std::string xmlrow=sxml->ErrGetFileRow(ele);
  if(Count)Run_ExceptioonFile("Only one description zone is allowed for inlet/outlet points",xmlrow);
  //-Checks XML elements.
  if(Simulate2D)sxml->CheckElementNames(ele,true,"direction rotate");
  else sxml->CheckElementNames(ele,true,"direction rotateaxis rotateadv");
  //-Load basic data.
  const unsigned mkfluid=sxml->GetAttributeUint(ele,"mkfluid");
  const string strdir=fun::StrLower(sxml->GetAttributeStr(ele,"direction"));
  //-Load direction.
  tdouble3 dir=DirectionFromStr(strdir);
  if(dir.x==DBL_MAX)sxml->ErrReadAtrib(ele,"direction",false);
  //-Loads custom direction.
  const bool dircustom=(dir==TDouble3(0));
  if(dircustom){
    dir=TDouble3(0);
    if(sxml->ExistsElement(ele,"direction")){
      dir.x=sxml->ReadElementFloat(ele,"direction","x");
      dir.z=sxml->ReadElementFloat(ele,"direction","z");
      if(!Simulate2D)dir.y=sxml->ReadElementFloat(ele,"direction","y",true);
    }
    else if(!Simulate2D){
      dir=ComputeDir3dFromParts(mkfluid,partsdata);
      string msg="Custom direction vector definition is missing. ";
      if(dir==TDouble3(0))msg=msg+"A valid vector could not be calculated automatically from the selected particles.";
      else msg=msg+fun::PrintStr("The vector (%g,%g,%g) calculated from the selected particles is suggested."
        ,dir.x,dir.y,dir.z);
      sxml->ErrReadElement(ele,"direction",false,msg);
    }
    if(dir==TDouble3(0))sxml->ErrReadAtrib(ele,"direction",false,"Custom direction vector is null.");
  }
  else if(sxml->ExistsElement(ele,"direction"))sxml->ErrReadElement(ele,"direction"
    ,false,"Direction definition is valid for custom option only.");
  if(dir==TDouble3(0))sxml->ErrReadAtrib(ele,"direction",false,"Direction vector is null.");

  //-Check mkfluid information.
  if(!partsdata || partsdata->GetNp()==0)
    Run_ExceptioonFile("No particles data to define inout points.",xmlrow);
  const JSphMk* mkinfo=partsdata->GetMkInfo();
  const unsigned cmk=mkinfo->GetMkBlockByMkFluid(mkfluid);
  const JSphMkBlock* pmk=(cmk<mkinfo->Size()? mkinfo->Mkblock(cmk): NULL);
  if(!pmk || !pmk->Count)sxml->ErrReadAtrib(ele,"mkfluid",false
    ,fun::PrintStr("No particles data with mkfluid=%u to define inout points.",mkfluid));

  //-Adds config information.
  const tdouble3 pmin=pmk->GetPosMin(),pmax=pmk->GetPosMax();
  if(Simulate2D)ConfigInfo.push_back(fun::PrintStr("Initial particles mkfluid=%u (x,z): (%g,%g)-(%g,%g)"
    ,mkfluid,pmin.x,pmin.z,pmax.x,pmax.z));
  else          ConfigInfo.push_back(fun::PrintStr("Initial particles mkfluid=%u (x,y,z): (%g,%g,%g)-(%g,%g,%g)"
    ,mkfluid,pmin.x,pmin.y,pmin.z,pmax.x,pmax.y,pmax.z));

  //-Creates points.
  {
    //-Calculates nunmber of inlet points.
    unsigned npt=pmk->Count;
    //-Resize memory when its size is not enough.
    ResizeMemory(npt);
    //-Get position of points.
    const typecode rcode=pmk->Code;
    const unsigned np=partsdata->GetNp();
    const typecode* code=partsdata->GetCode();
    const tdouble3* pos=partsdata->GetPos();
    unsigned cp=0;
    for(unsigned p=0;p<np;p++)if(code[p]==rcode){
      if(cp<npt)Points[cp]=pos[p];
      cp++;
    }
    if(cp!=npt)Run_ExceptioonFile("Error in number of particles.",xmlrow);
    Count=npt;
  }

  //-Loads 2-D or 3-D rotation configuration from XML.
  string rotinfo;
  const JMatrix4d m=(Simulate2D? ReadRotate2D(sxml,ele,TDouble3(0),rotinfo): 
                                 ReadRotate3D(sxml,ele,rotinfo));
  //-Applies rotation to points and dir vector.
  if(!rotinfo.empty()){
    ConfigInfo.push_back(rotinfo);
    const tdouble3 p0=Points[0];
    for(unsigned p=0;p<Count;p++)Points[p]=m.MulPoint(Points[p]);
    dir=fgeo::VecUnitary(m.MulPoint(p0+dir)-m.MulPoint(p0));
  }

  //-Updates Direction.
  if(dir==TDouble3(0))sxml->ErrReadElement(ele,"direction",false,"Direction is not a valid vector.");
  Direction=fgeo::VecUnitary(dir);

  //-Check points in line for 2D.
  if(Simulate2D)CheckPointsInLine(Direction,Count,Points,xmlrow);
  //-Check points in plane for 3D.
  CheckPointsInPlane(Direction,Count,Points,xmlrow);
  //-Compute domain limits of zone starting from inout points.
  ComputeDomainFromPoints();
}

//==============================================================================
/// Creates points in a line for 2D inlet.
//==============================================================================
void JSphInOutPoints::Create2d_Line(const JXml* sxml,TiXmlElement* ele){
  const std::string xmlrow=sxml->ErrGetFileRow(ele);
  if(Count)Run_ExceptioonFile("Only one description zone is allowed for inlet/outlet points",xmlrow);
  //-Checks XML elements.
  sxml->CheckElementNames(ele,true,"point point2 direction rotate");
  //-Load basic data.
  double px1=sxml->ReadElementFloat(ele,"point","x");
  double pz1=sxml->ReadElementFloat(ele,"point","z");
  double px2=sxml->ReadElementFloat(ele,"point2","x");
  double pz2=sxml->ReadElementFloat(ele,"point2","z");
  //-Adds config information.
  ConfigInfo.push_back(fun::PrintStr("Line(x,z): (%g,%g)-(%g,%g)",px1,pz1,px2,pz2));
  //-Load or compute direction.
  tdouble3 dir=TDouble3(0);
  if(sxml->ExistsElement(ele,"direction")){
    dir.x=sxml->ReadElementFloat(ele,"direction","x");
    dir.z=sxml->ReadElementFloat(ele,"direction","z");
    dir=fgeo::VecUnitary(dir);
  }
  else{
    const tdouble3 v=TDouble3(px2-px1,0,pz2-pz1);
    dir=fgeo::VecUnitary(TDouble3(-v.z,0,v.x));
  }
  if(dir==TDouble3(0))sxml->ErrReadElement(ele,"direction",false
    ,"Loaded or computed direction is not a valid vector.");

  //-Creates points.
  {
    //-Calculates nunmber of inlet points.
    const tdouble3 pt1=TDouble3(px1,Simulate2DPosY,pz1);
    const tdouble3 pt2=TDouble3(px2,Simulate2DPosY,pz2);
    const double dist=fgeo::PointsDist(pt2,pt1);
    double dppoints=Dp;
    unsigned npt=unsigned(dist/dppoints);
    if(dist-Dp*npt>Dp*0.95){
      npt++;
      dppoints=dist/double(npt);
    }
    npt++;
    //-Resize memory when its size is not enough.
    ResizeMemory(npt);
    //-Calculates position of points.
    const tdouble3 vec=fgeo::VecUnitary(pt2-pt1)*dppoints;
    if(InitialMove)Run_Exceptioon("InitialMove is non-zero...");
    //const tdouble3 inimove=(Direction*InitialMove);
    for(unsigned p=0;p<npt;p++)Points[p]=pt1+(vec*double(p));
    Count=npt;
  }

  //-Loads 2-D rotation configuration from XML.
  string rotinfo;
  const JMatrix4d m=ReadRotate2D(sxml,ele,TDouble3(px1,0,pz1),rotinfo);
  //-Applies rotation to points and dir vector.
  if(!rotinfo.empty()){
    ConfigInfo.push_back(rotinfo);
    const tdouble3 p0=Points[0];
    for(unsigned p=0;p<Count;p++)Points[p]=m.MulPoint(Points[p]);
    dir=fgeo::VecUnitary(m.MulPoint(p0+dir)-m.MulPoint(p0));
  }

  //-Updates Direction.
  if(dir==TDouble3(0))sxml->ErrReadElement(ele,"direction",false,"Direction is not a valid vector.");
  Direction=fgeo::VecUnitary(dir);

  //-Check points in line 2D.
  CheckPointsInLine(Direction,Count,Points,xmlrow);
  //-Compute domain limits of zone starting from inout points.
  ComputeDomainFromPoints();
}

//==============================================================================
/// Returns matrix for rotation in 3D.
//==============================================================================
JMatrix4d JSphInOutPoints::ReadRotate3D(const JXml* sxml,TiXmlElement* ele
  ,std::string& rotationinfo)
{
  JMatrix4d m;
  const bool rotaxis=sxml->ExistsElement(ele,"rotateaxis");
  const bool rotadv =sxml->ExistsElement(ele,"rotateadv");
  if(rotaxis && rotadv)sxml->ErrReadElement(ele,"rotateaxis",false
    ,"Only one rotation configuration is valid (<rotateaxis> or <rotateadv>)."); 
  if(rotaxis){
    double rotate=sxml->ReadElementDouble(ele,"rotateaxis","angle",true);
    if(rotate){
      const string angunits=fun::StrLower(sxml->ReadElementStr(ele,"rotateaxis","anglesunits"));
      if(angunits=="radians")rotate=rotate*TODEG;
      else if(angunits!="degrees")sxml->ErrReadElement(ele,"rotate",false
        ,"The value anglesunits must be \"degrees\" or \"radians\"."); 
      TiXmlElement* rot=ele->FirstChildElement("rotateaxis");
      tdouble3 pt1=sxml->ReadElementDouble3(rot,"point1");
      tdouble3 pt2=sxml->ReadElementDouble3(rot,"point2");
      m=JMatrix4d::MatrixRot(rotate,pt1,pt2);
      //-Return config information about rotation.
      rotationinfo=fun::PrintStr("  Rotated: angle:%g axis:%s",rotate,fun::Double3gRangeStr(pt1,pt2).c_str());
    }
  }
  if(rotadv){
    //-Loads angles.
    double angle1=sxml->ReadElementDouble(ele,"rotateadv","angle1",false);
    double angle2=sxml->ReadElementDouble(ele,"rotateadv","angle2",true);
    double angle3=sxml->ReadElementDouble(ele,"rotateadv","angle3",true);
    if(angle1 || angle2 || angle3){
      //-Loads angles units.
      const string angunits=fun::StrLower(sxml->ReadElementStr(ele,"rotateadv","anglesunits"));
      if(angunits=="radians"){
        angle1=angle1*TODEG;
        angle2=angle2*TODEG;
        angle3=angle3*TODEG;
      }
      else if(angunits!="degrees")sxml->ErrReadElement(ele,"rotateadv",false
        ,"The value anglesunits must be \"degrees\" or \"radians\"."); 
      //-Loads axes and intrinsic.
      const string axes=fun::StrUpper(sxml->ReadElementStr(ele,"rotateadv","axes"));
      if(axes.size()!=3 || !JMatrix4d::CheckRotateAxes(axes.c_str()))
        sxml->ErrReadElement(ele,"rotateadv",false,"The axes value is invalid."); 
      const bool intrinsic=sxml->ReadElementBool(ele,"rotateadv","intrinsic");
      //-Loads center.
      TiXmlElement* rot=ele->FirstChildElement("rotateadv");
      const tdouble3 center=sxml->ReadElementDouble3(rot,"center",true,TDouble3(0));
      //-Compute rotation matrix.
      if(center==TDouble3(0))m=JMatrix4d::MatrixRotate(angle1,angle2,angle3,axes.c_str(),intrinsic);
      else m=JMatrix4d::MatrixRotateCen(center,angle1,angle2,angle3,axes.c_str(),intrinsic);
      //-Return config information about rotation.
      rotationinfo=fun::PrintStr("  Rotated: angles:(%g,%g,%g) axis:%s intrinsic:%s"
        ,angle1,angle2,angle3,axes.c_str(),(intrinsic? "True": "False"));
      if(center!=TDouble3(0))rotationinfo=rotationinfo+fun::PrintStr(" center:(%s)",fun::Double3gStr(center).c_str());
    }
  }
  return(m);
}

//==============================================================================
/// Creates points in a box for 3D inlet.
//==============================================================================
void JSphInOutPoints::Create3d_Box(const JXml* sxml,TiXmlElement* ele){
  const std::string xmlrow=sxml->ErrGetFileRow(ele);
  if(Count)Run_ExceptioonFile("Only one description zone is allowed for inlet/outlet points",xmlrow);
  //-Checks XML elements.
  sxml->CheckElementNames(ele,true,"point size direction rotateaxis rotateadv");
  //-Load basic data.
  tdouble3 pt0=sxml->ReadElementDouble3(ele,"point");
  tdouble3 spt=sxml->ReadElementDouble3(ele,"size");
  if(spt.x!=0 && spt.y!=0 && spt.z!=0)sxml->ErrReadElement(ele,"size",false,"One size of axis must be zero."); 
  tdouble3 ptx=pt0+TDouble3(spt.x,0,0);
  tdouble3 pty=pt0+TDouble3(0,spt.y,0);
  tdouble3 ptz=pt0+TDouble3(0,0,spt.z);
  //-Adds config information.
  ConfigInfo.push_back(fun::PrintStr("Box: p0:(%s) size:(%s)",fun::Double3gStr(pt0).c_str(),fun::Double3gStr(spt).c_str()));
  //-Load direction.
  tdouble3 dir=fgeo::VecUnitary(sxml->ReadElementDouble3(ele,"direction"));
  if(dir==TDouble3(0))sxml->ErrReadElement(ele,"direction",false,"Loaded direction is not a valid vector.");

  //-Creates points.
  {
    //-Calculates nunmber of inlet points.
    double dpx=Dp,dpy=Dp,dpz=Dp;
    unsigned npx=unsigned(spt.x/dpx);
    unsigned npy=unsigned(spt.y/dpy);
    unsigned npz=unsigned(spt.z/dpz);
    if(spt.x-dpx*npx>dpx*0.95){ npx++; dpx=spt.x/double(npx); }
    if(spt.y-dpy*npy>dpy*0.95){ npy++; dpy=spt.y/double(npy); }
    if(spt.z-dpz*npz>dpz*0.95){ npz++; dpz=spt.z/double(npz); }
    npx++; npy++; npz++;
    unsigned npt=npx*npy*npz;
    //-Resize memory when its size is not enough.
    ResizeMemory(npt);
    //-Calculates position of points.
    const tdouble3 vx=(spt.x? fgeo::VecUnitary(ptx-pt0)*dpx: TDouble3(0));
    const tdouble3 vy=(spt.y? fgeo::VecUnitary(pty-pt0)*dpy: TDouble3(0));
    const tdouble3 vz=(spt.z? fgeo::VecUnitary(ptz-pt0)*dpz: TDouble3(0));
    unsigned p=0;
    for(unsigned pz=0;pz<npz;pz++)for(unsigned py=0;py<npy;py++)for(unsigned px=0;px<npx;px++){
      Points[p]=pt0+(vx*px)+(vy*py)+(vz*pz);
      p++;
    }
    if(p!=npt)Run_Exceptioon("Error calculating number of points.");
    Count=npt;
  }

  //-Loads 3-D rotation configuration from XML.
  string rotinfo;
  const JMatrix4d m=ReadRotate3D(sxml,ele,rotinfo);
  //-Applies rotation to points and dir vector.
  if(!rotinfo.empty()){
    ConfigInfo.push_back(rotinfo);
    const tdouble3 p0=Points[0];
    for(unsigned p=0;p<Count;p++)Points[p]=m.MulPoint(Points[p]);
    dir=fgeo::VecUnitary(m.MulPoint(p0+dir)-m.MulPoint(p0));
  }

  //-Updates Direction.
  if(dir==TDouble3(0))sxml->ErrReadElement(ele,"direction",false,"Direction is not a valid vector.");
  Direction=fgeo::VecUnitary(dir);

  //-Check points in plane 3D.
  CheckPointsInPlane(Direction,Count,Points,xmlrow);
  //-Compute domain limits of zone starting from inout points.
  ComputeDomainFromPoints();
}

//==============================================================================
/// Creates points in a circle.
//==============================================================================
void JSphInOutPoints::Create3d_Circle(const JXml* sxml,TiXmlElement* ele){
  const std::string xmlrow=sxml->ErrGetFileRow(ele);
  if(Count)Run_ExceptioonFile("Only one description zone is allowed for inlet/outlet points",xmlrow);
  //-Checks XML elements.
  sxml->CheckElementNames(ele,true,"point radius direction rotateaxis rotateadv");
  //-Load basic data.
  const tdouble3 pt0=sxml->ReadElementDouble3(ele,"point");
  const double radius=sxml->ReadElementDouble(ele,"radius","v");
  CircleRadius=radius;
  tdouble3 pcen=pt0;
  //-Adds config information.
  ConfigInfo.push_back(fun::PrintStr("Circle: center:(%s) radius:%g",fun::Double3gStr(pcen).c_str(),radius));
  //-Load direction.
  tdouble3 dir=fgeo::VecUnitary(sxml->ReadElementDouble3(ele,"direction"));

  //-Creates points.
  {
    //-Calculates nunmber of inlet points.
    unsigned npt=0;
    const unsigned nra=max(1u,unsigned((radius+Dp*0.25)/Dp));
    const double rasum=radius/nra;
    for(unsigned cr=0;cr<=nra;cr++){
      const double ra=rasum*cr;
      const unsigned nang=max(1u,unsigned((TWOPI*ra+Dp/2)/Dp));
      npt+=nang;
    }
    //-Resize memory when its size is not enough.
    ResizeMemory(npt);
    //-Calculates position of points.
    unsigned p=0;
    for(unsigned cr=0;cr<=nra;cr++){
      const double ra=rasum*cr;
      const unsigned nang=max(1u,unsigned((TWOPI*ra+Dp/2)/Dp));
      const double angsum=TWOPI/nang; //-In radians.
      const tdouble3 v2=fgeo::VecOrthogonal2(dir,1,true);
      const tdouble3 v3=fgeo::ProductVec(dir,v2);
      for(unsigned c=0;c<nang;c++){
        const double angle=angsum*c;
        Points[p]=pt0 + (v2*(ra*cos(angle))) + (v3*(ra*sin(angle)));
        p++;
      }
    }
    if(p!=npt)Run_Exceptioon("Error calculating number of points.");
    Count=npt;
  }

  //-Loads 3-D rotation configuration from XML.
  string rotinfo;
  const JMatrix4d m=ReadRotate3D(sxml,ele,rotinfo);
  //-Applies rotation to points and dir vector.
  if(!rotinfo.empty()){
    ConfigInfo.push_back(rotinfo);
    const tdouble3 p0=Points[0];
    for(unsigned p=0;p<Count;p++)Points[p]=m.MulPoint(Points[p]);
    dir=fgeo::VecUnitary(m.MulPoint(p0+dir)-m.MulPoint(p0));
  }

  //-Updates Direction.
  if(dir==TDouble3(0))sxml->ErrReadElement(ele,"direction",false,"Direction is not a valid vector.");
  Direction=fgeo::VecUnitary(dir);

  //-Check points in plane 3D.
  CheckPointsInPlane(Direction,Count,Points,xmlrow);
  //-Compute domain limits of zone starting from inout points.
  ComputeDomainFromPoints();
}

//==============================================================================
/// Reads definition of inlet points in the XML node and creates points.
//==============================================================================
void JSphInOutPoints::CreatePoints(const JXml* sxml,TiXmlElement* lis
  ,const JDsPartsInit* partsdata)
{
  string xmlrow=sxml->ErrGetFileRow(lis);
  XmlShape="";
  TiXmlElement* ele=lis->FirstChildElement();
  while(ele){
    string cmd=ele->Value();
    if(cmd.length() && cmd[0]!='_'){
      if(!XmlShape.empty())Run_ExceptioonFile("Only one inlet shape is allowed per zone.",xmlrow);
      if(Simulate2D){//-Loads inflow definition for 2D.
        if(cmd=="particles")Create2d3d_Particles(sxml,ele,partsdata);
        else if(cmd=="line")Create2d_Line(sxml,ele);
        else sxml->ErrReadElement(ele,cmd,false);
      }
      else{//-Loads inflow definition for 3D.
        if(cmd=="particles")Create2d3d_Particles(sxml,ele,partsdata);
        else if(cmd=="box")Create3d_Box(sxml,ele);
        else if(cmd=="circle")Create3d_Circle(sxml,ele);
        else sxml->ErrReadElement(ele,cmd,false);
      }
      XmlShape=fun::StrUpper(cmd);
    }
    ele=ele->NextSiblingElement();
  }
  //-Checks direction and position of points in simulation domain.
  CheckPoints(xmlrow);
  //-Computes domain limits of inlet points.
  ComputeDomainLimits(ZonePosMin,ZonePosMax);
  //-Set PointsInit[] to zero.
  memset(PointsInit,0,sizeof(byte)*Size);
}

//==============================================================================
/// Compute domain limits from inout points.
//==============================================================================
void JSphInOutPoints::ComputeDomainLimits(tdouble3& posmin,tdouble3& posmax)const{
  if(Count==0)Run_Exceptioon("There are not defined points.");
  tdouble3 pmin=Points[0],pmax=Points[0];
  //-Calculates minimum and maximum position of inout points. 
  for(unsigned p=1;p<Count;p++){
    const tdouble3 ps=Points[p];
    if(pmin.x>ps.x)pmin.x=ps.x;
    if(pmin.y>ps.y)pmin.y=ps.y;
    if(pmin.z>ps.z)pmin.z=ps.z;
    if(pmax.x<ps.x)pmax.x=ps.x;
    if(pmax.y<ps.y)pmax.y=ps.y;
    if(pmax.z<ps.z)pmax.z=ps.z;
  }
  posmin=pmin; posmax=pmax;
}

//==============================================================================
/// Find end points of the inlet line in list of points.
//==============================================================================
void JSphInOutPoints::FindLineExtremes(tdouble3 dirline,unsigned np
  ,const tdouble3* points,tdouble3& ptmin,tdouble3& ptmax)const
{
  const tplane3d pla=fgeo::PlaneNormalized(fgeo::PlanePtVec(points[0],dirline));
  unsigned pmin=0,pmax=0;
  double dmin=(pla.a*points[0].x+pla.b*points[0].y+pla.c*points[0].z+pla.d);
  double dmax=dmin;
  for(unsigned p=1;p<np;p++){
    const tdouble3 pt=points[p];
    const double d=(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d);
    if(dmax<d){ dmax=d; pmax=p; }
    if(dmin>d){ dmin=d; pmin=p; }
  }
  ptmin=points[pmin];
  ptmax=points[pmax];
}

//==============================================================================
/// Compute packing box adjusted to the inlet points.
//==============================================================================
void JSphInOutPoints::ComputePackingBox(tdouble3 direction,unsigned np
  ,const tdouble3* points,tdouble3& boxpt,tdouble3& boxv1,tdouble3& boxv2
  ,tdouble3& boxv3)const
{
  //-Computes midpoint. 
  tdouble3 midp=TDouble3(0);
  for(unsigned p=0;p<np;p++)midp=midp+points[p];
  midp=midp/np;
  boxpt=midp;
  //-Computes axes orthogonal to the direction.
  boxv1=fgeo::VecUnitary(direction);
  boxv2=fgeo::VecOrthogonal2(boxv1,1,true);
  boxv3=fgeo::ProductVec(boxv1,boxv2);
  //-Computes minimum and maximum position for each axis.
  double angle=0;
  tdouble3 vec2=boxv2,vec3=boxv3;
  double surfmin=DBL_MAX;
  double dist1min=0,dist1max=0;
  double dist2min=0,dist2max=0;
  double dist3min=0,dist3max=0;
  const tdouble3 p0=boxpt; 
  const tdouble3 p1=boxpt+boxv1; 
  for(unsigned cang=0;cang<180;cang++){
    const double ang=double(cang)/2;
    tdouble3 v2=boxv2;
    tdouble3 v3=boxv3;
    if(ang){
      const JMatrix4d m=JMatrix4d::MatrixRot(ang,p0,p1);
      v2=fgeo::VecUnitary(m.MulPoint(p0+v2)-p0);
      v3=fgeo::VecUnitary(m.MulPoint(p0+v3)-p0);
    }
    double d2min,d2max;
    double d3min,d3max;
    ComputeVectorExtremes(np,points,p0,v2,d2min,d2max);
    ComputeVectorExtremes(np,points,p0,v3,d3min,d3max);
    const double d2=d2max-d2min;
    const double d3=d3max-d3min;
    const double surf=d2*d3;
    if(surf<surfmin || cang==0){
      angle=ang;  vec2=v2;  vec3=v3;
      surfmin=surf;
      dist2min=d2min; dist2max=d2max;
      dist3min=d3min; dist3max=d3max;
      ComputeVectorExtremes(np,points,p0,boxv1,dist1min,dist1max);
    }
  }
  //-Computes packing box to return.
  boxpt=boxpt+(boxv1*dist1min)+(vec2*dist2min)+(vec3*dist3min);
  boxv1=boxv1*(dist1max-dist1min);
  boxv2=vec2*(dist2max-dist2min);
  boxv3=vec3*(dist3max-dist3min);
  //printf("====--> ANG:%g  d2=(%g,%g)  d3=(%g,%g)  surf:%g\n\n",
  //  angle,dist2min,dist2max,dist3min,dist3max,surfmin);
}

//==============================================================================
/// Compute domain limits of zone from inout points (for 2D and particles-3D).
//==============================================================================
void JSphInOutPoints::ComputeDomainFromPoints(){
  //-Calculates minimum and maximum position of inout points. 
  tdouble3 pmin=TDouble3(DBL_MAX),pmax=TDouble3(-DBL_MAX);
  ComputeDomainLimits(pmin,pmax);
  const tdouble3 dir1=Direction*(Dp/2); //-Vector direction.
  const unsigned nfar=2*(Layers)+1;
  if(Simulate2D){
    const tdouble3 dirline=fgeo::VecUnitary(TDouble3(Direction.z,0,-Direction.x));
    //-Find end points of the inlet line in list of points.
    tdouble3 pmin,pmax;
    FindLineExtremes(dirline,Count,Points,pmin,pmax);
    const tdouble3 dir2=dirline*(Dp/2);
    const tdouble3 pnearmin=pmin+dir1-dir2;
    const tdouble3 pnearmax=pmax+dir1+dir2;
    const tdouble3 pfarmin=pnearmin-(dir1*nfar);//-Increases one layer and adds other dp for borders.
    const tdouble3 pfarmax=pnearmax-(dir1*nfar);//-Increases one layer and adds other dp for borders.
    const tdouble3 ptzero=pnearmin;
    const tdouble3 dirup=pnearmax-ptzero;
    const tdouble3 dirfar=pfarmin-ptzero;
    const tdouble3 dirlat=TDouble3(0);
    //-Stores domain points.
    PtDom[0]=ptzero;
    PtDom[1]=PtDom[0]+dirup;
    PtDom[2]=PtDom[0]+dirup+dirfar;
    PtDom[3]=PtDom[0]+dirfar;
    PtDom[4]=PtDom[5]=PtDom[6]=PtDom[7]=TDouble3(0);
    //-Stores normal points.
    PtDom[8]=(PtDom[1]+PtDom[0])/2;
    PtDom[9]=PtDom[8]+(Direction*(Dp*Layers));
  }
  else{
    tdouble3 boxpt,boxv1,boxv2,boxv3;
    ComputePackingBox(Direction,Count,Points,boxpt,boxv1,boxv2,boxv3);
    const tdouble3 v2=fgeo::VecUnitary(boxv2);
    const tdouble3 v3=fgeo::VecUnitary(boxv3);
    boxv2=boxv2+v2*Dp;
    boxv3=boxv3+v3*Dp;
    //-Stores domain points.
    PtDom[0]=boxpt+dir1-v2*(Dp/2)-v3*(Dp/2);
    PtDom[1]=PtDom[0]+boxv2;
    PtDom[2]=PtDom[0]+boxv2+boxv3;
    PtDom[3]=PtDom[0]+boxv3;
    PtDom[4]=PtDom[0]-(dir1*nfar);
    PtDom[5]=PtDom[1]-(dir1*nfar);
    PtDom[6]=PtDom[2]-(dir1*nfar);
    PtDom[7]=PtDom[3]-(dir1*nfar);
    //-Stores normal points.
    PtDom[8]=(PtDom[0]+PtDom[2])/2;
    PtDom[9]=PtDom[8]+(Direction*(Dp*Layers));
  }
}

//==============================================================================
/// Checks direction and position of points in simulation domain.
//==============================================================================
void JSphInOutPoints::CheckPointsInLine(tdouble3 dir,unsigned np
  ,const tdouble3* points,string xmlrow)const
{
  if(!np)Run_ExceptioonFile("There are no points",xmlrow);
  //-Computes maximum distance to plane.
  const tplane3d pla=fgeo::PlaneNormalized(fgeo::PlanePtVec(points[0],dir));
  double dmax=fabs(pla.a*points[0].x+pla.b*points[0].y+pla.c*points[0].z+pla.d);
  for(unsigned p=1;p<np;p++){
    const tdouble3 pt=points[p];
    const double d=fabs(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d);
    if(dmax<d)dmax=d;
  }
  //-Show error.
  if(dmax>Dp*0.01){
    //-Find end points of the inlet line in list of points.
    const tdouble3 dirline=fgeo::VecUnitary(TDouble3(dir.z,0,-dir.x));
    tdouble3 pmin,pmax;
    FindLineExtremes(dirline,Count,Points,pmin,pmax);
    const tplane3d pla2=fgeo::PlaneNormalized(fgeo::PlanePtVec(pmin,dir));
    //-Creates VTK file with inlet line.
    {
      JSpVtkShape ss;
      ss.AddLine(pmin,pmin+dirline*fgeo::PointsDist(pmin,pmax));
      ss.SaveVtk(AppInfo.GetDirOut()+"ErrorInOut_InoutPointsInLine_Line.vtk");
    }
    //-Creates VTK file with inlet points.
    JDataArrays arrays;
    arrays.AddArray("Pos",np,points,false);
    double* dist=arrays.CreateArrayPtrDouble("Distance",np,false);
    for(unsigned p=0;p<np;p++){
      const tdouble3 pt=points[p];
      dist[p]=fabs(pla2.a*pt.x+pla2.b*pt.y+pla2.c*pt.z+pla2.d);
    } 
    const string file=AppInfo.GetDirOut()+"ErrorInOut_InoutPointsInLine_Points.vtk";
    JSpVtkData::Save(file,arrays,"Pos");
    Log->AddFileInfo(file,"Saves invalid InOut points (too far from the inlet line).");
    Run_ExceptioonFile(fun::PrintStr("All inlet points must be on the same line orthogonal to the direction vector. Check distance to the inlet line in the VTK file \'%s\'.",file.c_str()),xmlrow);
  }
}

//==============================================================================
// Crea VTK simple para pruebas.
//==============================================================================
void JSphInOutPoints::ComputeVectorExtremes(unsigned np,const tdouble3* points
  ,tdouble3 pt0,tdouble3 vec,double& dmin,double& dmax)const
{
  if(!np)Run_Exceptioon("There are no points");
  const tplane3d pla=fgeo::PlaneNormalized(fgeo::PlanePtVec(pt0,vec));
  dmin=(pla.a*pt0.x+pla.b*pt0.y+pla.c*pt0.z+pla.d);
  dmax=dmin;
  for(unsigned p=1;p<np;p++){
    const tdouble3 pt=points[p];
    const double d=(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d);
    if(dmax<d)dmax=d;
    if(dmin>d)dmin=d;
  }
}

//==============================================================================
/// Checks direction and position of points in simulation domain.
//==============================================================================
void JSphInOutPoints::CheckPointsInPlane(tdouble3 dir,unsigned np
  ,const tdouble3* points,string xmlrow)const
{
  if(!np)Run_ExceptioonFile("There are no points",xmlrow);
  //-Computes maximum distance to plane.
  const tplane3d pla=fgeo::PlaneNormalized(fgeo::PlanePtVec(points[0],dir));
  double dmax=fabs(pla.a*points[0].x+pla.b*points[0].y+pla.c*points[0].z+pla.d);
  for(unsigned p=1;p<np;p++){
    const tdouble3 pt=points[p];
    const double d=fabs(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d);
    if(dmax<d)dmax=d;
  }
  //-Show error.
  if(dmax>Dp*0.01){
    const tdouble3 p0=points[0];
    const tdouble3 v2=fgeo::VecOrthogonal2(dir,1,true);
    const tdouble3 v3=fgeo::ProductVec(dir,v2);
    double d2min,d2max;
    double d3min,d3max;
    ComputeVectorExtremes(np,points,p0,v2,d2min,d2max);
    ComputeVectorExtremes(np,points,p0,v3,d3min,d3max);
    //-Creates VTK file with inlet plane.
    {
      JSpVtkShape ss;
      const tdouble3 p1=p0+v2*d2min+v3*d3min;
      const tdouble3 p2=p0+v2*d2max+v3*d3min;
      const tdouble3 p3=p0+v2*d2max+v3*d3max;
      const tdouble3 p4=p0+v2*d2min+v3*d3max;
      ss.AddQuad(p1,p2,p3,p4);
      const tdouble3 pm=(p1+p2+p3+p4)/4.;
      ss.AddLine(pm,pm+dir*(Dp*5));
      ss.SaveVtk(AppInfo.GetDirOut()+"ErrorInOut_InoutPointsInPlane_Plane.vtk");
    }
    //-Creates VTK file with inlet points.
    JDataArrays arrays;
    arrays.AddArray("Pos",np,points,false);
    double* dist=arrays.CreateArrayPtrDouble("Distance",np,false);
    for(unsigned p=0;p<np;p++){
      const tdouble3 pt=points[p];
      dist[p]=fabs(pla.a*pt.x+pla.b*pt.y+pla.c*pt.z+pla.d);
    } 
    const string file=AppInfo.GetDirOut()+"ErrorInOut_InoutPointsInPlane_Points.vtk";
    JSpVtkData::Save(file,arrays,"Pos");
    Log->AddFileInfo(file,"Saves invalid InOut points (too far from the inlet plane).");
    Run_ExceptioonFile(fun::PrintStr("All inlet points must be on the same plane orthogonal to the direction vector. Check distance to the inlet plane in the VTK file \'%s\'.",file.c_str()),xmlrow);
  }
}

//==============================================================================
/// Checks direction and position of points in simulation domain.
//==============================================================================
void JSphInOutPoints::CheckPoints(const std::string& xmlrow){
  //-Checks direction.
  if(Simulate2D && Direction.y!=0)Run_ExceptioonFile("Direction.y is not zero.",xmlrow);
  if(Direction==TDouble3(0))Run_ExceptioonFile("Direction vector is zero.",xmlrow);
  if(Count==0)Run_ExceptioonFile("There are not defined points.",xmlrow);
  //-Checks domain using one layer more because it is the real limit.
  bool error=false;
  tdouble3 errpt;
  for(unsigned c=0;c<=Layers && !error;c++){
    const tdouble3 sub=(Direction*double(Dp*c+InitialMove));
    for(unsigned p=0;p<Count && !error;p++){
      const tdouble3 ps=Points[p]-sub;
      if(!(MapRealPosMin<ps && ps<MapRealPosMax)){ error=true; errpt=ps; }
    }
  }
  //-If it fails then saves VTK with points.
  if(error){
    //-Allocates memory.
    const unsigned np=Count*(Layers+1);
    tfloat3* pos=new tfloat3[np];
    byte*    layer=new byte[np];
    byte*    outside=new byte[np];
    //-Loads point data.
    for(unsigned c=0;c<=Layers;c++){
      const tdouble3 sub=(Direction*double(Dp*c+InitialMove));
      for(unsigned p=0;p<Count;p++){
        const tdouble3 ps=Points[p]-sub;
        const unsigned pp=Count*c+p;
        pos[pp]=ToTFloat3(ps);
        layer[pp]=byte(c);
        outside[pp]=byte(!(MapRealPosMin<ps && ps<MapRealPosMax)? 1: 0);
      }
    }
    //-Creates VTK file.
    JDataArrays arrays;
    arrays.AddArray("Pos",np,pos,false);
    if(outside)arrays.AddArray("Outside",np,outside,false);
    if(layer)  arrays.AddArray("Layer"  ,np,layer  ,false);
    const string file=AppInfo.GetDirOut()+"ErrorInOut_InoutPoints.vtk";
    JSpVtkData::Save(file,arrays,"Pos");
    Log->AddFileInfo(file,"Saves invalid InOut points (outside the domain).");
    //-Frees memory.
    delete[] pos;     pos=NULL;
    delete[] layer;   layer=NULL;
    delete[] outside; outside=NULL;
    Run_ExceptioonFile(fun::PrintStr("Point for inlet conditions with position (%g,%g,%g) is outside the domain. Checks the VTK file \'%s\'.",errpt.x,errpt.y,errpt.z,file.c_str()),xmlrow);
  }
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphInOutPoints::GetConfig(std::vector<std::string>& lines)const{
  lines.push_back("InOut point definition:");
  for(unsigned i=0;i<unsigned(ConfigInfo.size());i++)lines.push_back(string("  ")+ConfigInfo[i]);
}

//==============================================================================
/// Activates or deactivates initial points.
//==============================================================================
void JSphInOutPoints::SetPointsInit(bool active){
  memset(PointsInit,(active? 1: 0),sizeof(byte)*Count);
}

//==============================================================================
/// Counts and returns number of valid initial points.
//==============================================================================
unsigned JSphInOutPoints::CountPointsInit()const{
  unsigned n=0;
  for(unsigned p=0;p<Count;p++)if(PointsInit[p])n++;
  return(n);
}

//==============================================================================
/// Returns border points of the domain of inlet points.
//==============================================================================
void JSphInOutPoints::GetPtDomain(std::vector<tdouble3>& ptdom)const{
  ptdom.clear();
  for(unsigned p=0;p<10;p++)ptdom.push_back(PtDom[p]);
}


