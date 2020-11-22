//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2020 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
#include "JVtkLib.h"
 
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
/// Returns matrix for rotation in 2D.
//==============================================================================
JMatrix4d JSphInOutPoints::ReadRotate2D(const JXml *sxml,TiXmlElement* ele,const tdouble3 &pt){
  double rotate=sxml->ReadElementDouble(ele,"rotate","angle",true);
  string angunits=fun::StrLower(sxml->ReadElementStr(ele,"rotate","anglesunits"));
  if(angunits=="radians")rotate=rotate*TODEG;
  else if(angunits!="degrees")sxml->ErrReadElement(ele,"rotate",false,"The value anglesunits must be \"degrees\" or \"radians\"."); 
  const JMatrix4d m=JMatrix4d::MatrixRot(rotate,TDouble3(pt.x,0,pt.z),TDouble3(pt.x,1,pt.z));
  //-Adds config information about rotation.
  ConfigInfo.push_back(fun::PrintStr("  rotate: angle:%g",rotate));
  return(m);
}

//==============================================================================
/// Returns direction vector starting from direction in text.
//==============================================================================
tdouble3 JSphInOutPoints::DirectionFromStr(const std::string &strdir)const{
  tdouble3 dir=TDouble3(0);
  if     (strdir=="top"   )dir.z= 1;
  else if(strdir=="bottom")dir.z=-1;
  else if(strdir=="right" )dir.x= 1;
  else if(strdir=="left"  )dir.x=-1;
  else if(strdir=="back"  )dir.y= 1;
  else if(strdir=="front" )dir.y=-1;
  return(dir);
}

//==============================================================================
/// Returns direction vector starting from direction in text.
//==============================================================================
std::string JSphInOutPoints::CheckParticlesDirection(const JSphMkBlock *pmk,const tdouble3 &dir)const{
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
/// Creates points starting from special fluid particles.
//==============================================================================
void JSphInOutPoints::Create2d3d_Particles(const JXml *sxml,TiXmlElement* ele
  ,const JDsPartsInit *partsdata)
{
  if(Count)Run_ExceptioonFile("There are previous definitions of inout points.",sxml->ErrGetFileRow(ele));
  unsigned mkfluid=sxml->GetAttributeUint(ele,"mkfluid");
  string strdir=sxml->GetAttributeStr(ele,"direction");

  //-Load direction.
  const tdouble3 dir=DirectionFromStr(strdir);
  if(dir==TDouble3(0))sxml->ErrReadAtrib(ele,"direction",false);
  Direction=fgeo::VecUnitary(dir);

  //-Check mkfluid information.
  if(!partsdata || partsdata->GetNp()==0)Run_Exceptioon("No particles data to define inout points.");
  const JSphMk* mkinfo=partsdata->GetMkInfo();
  const unsigned cmk=mkinfo->GetMkBlockByMkFluid(mkfluid);
  const JSphMkBlock* pmk=(cmk<mkinfo->Size()? mkinfo->Mkblock(cmk): NULL);
  if(!pmk || !pmk->Count)sxml->ErrReadAtrib(ele,"mkfluid",false,fun::PrintStr("No particles data with mkfluid=%u to define inout points.",mkfluid));
  std::string error=CheckParticlesDirection(pmk,dir);
  if(!error.empty())Run_ExceptioonFile(error,sxml->ErrGetFileRow(ele));

  //-Adds config information.
  const tdouble3 pmin=pmk->GetPosMin(),pmax=pmk->GetPosMax();
  if(Simulate2D)ConfigInfo.push_back(fun::PrintStr("Initial particles mkfluid=%u (x,z): (%g,%g)-(%g,%g)",mkfluid,pmin.x,pmin.z,pmax.x,pmax.z));
  else          ConfigInfo.push_back(fun::PrintStr("Initial particles mkfluid=%u (x,y,z): (%g,%g,%g)-(%g,%g,%g)",mkfluid,pmin.x,pmin.y,pmin.z,pmax.x,pmax.y,pmax.z));

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
      if(cp<npt)Points[Count+cp]=pos[p];
      cp++;
    }
    if(cp!=npt)Run_ExceptioonFile("Error in number of particles.",sxml->ErrGetFileRow(ele));
    Count+=npt;
  }
  //-Compute domain limits of zone starting from inout points.
  ComputeDomainFromPoints();
}

//==============================================================================
/// Creates points in a line.
//==============================================================================
void JSphInOutPoints::Create2d_Line(const JXml *sxml,TiXmlElement* ele){
  if(Count)Run_Exceptioon("Only one description zone is allowed for inlet/outlet points.");
  //-Load basic data.
  double px1=sxml->ReadElementFloat(ele,"point","x");
  double pz1=sxml->ReadElementFloat(ele,"point","z");
  double px2=sxml->ReadElementFloat(ele,"point2","x");
  double pz2=sxml->ReadElementFloat(ele,"point2","z");
  //-Adds config information.
  ConfigInfo.push_back(fun::PrintStr("Line(x,z): (%g,%g)-(%g,%g)",px1,pz1,px2,pz2));
  //-Load direction.
  tdouble3 dir=TDouble3(0);
  if(sxml->ExistsElement(ele,"direction")){
    dir.x=sxml->ReadElementFloat(ele,"direction","x");
    dir.z=sxml->ReadElementFloat(ele,"direction","z");
  }
  //-Applies rotation to inlet definition.
  double rotate=sxml->ReadElementDouble(ele,"rotate","angle",true);
  if(rotate){
    const JMatrix4d m=ReadRotate2D(sxml,ele,TDouble3(px1,0,pz1));
    tdouble3 pt1=TDouble3(px1,Simulate2DPosY,pz1);
    tdouble3 pt2=TDouble3(px2,Simulate2DPosY,pz2);
    tdouble3 pdir=pt1+dir;
    pt1=m.MulPoint(pt1);
    pt2=m.MulPoint(pt2);
    pdir=m.MulPoint(pdir);
    if(dir!=TDouble3(0))dir=pdir-pt1;
    px1=pt1.x; pz1=pt1.z;
    px2=pt2.x; pz2=pt2.z;
  }
  //-Updates Direction.
  if(dir!=TDouble3(0)){
    dir=fgeo::VecUnitary(dir);
    if(Direction!=TDouble3(0) && Direction!=dir)sxml->ErrReadElement(ele,"direction",false,"Direction is already defined.");
    Direction=dir;
  }
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
    const tdouble3 inimove=(Direction*InitialMove);
    for(unsigned p=0;p<npt;p++){
      const tdouble3 ps=pt1+(vec*TDouble3(p));
      Points[Count+p]=ps+inimove;
    }
    Count+=npt;
  }
  //-Compute domain limits of zone starting from inout points.
  ComputeDomainFromPoints();
}

//==============================================================================
/// Returns matrix for rotation in 3D.
//==============================================================================
JMatrix4d JSphInOutPoints::ReadRotate3D(const JXml *sxml,TiXmlElement* ele){
  double rotate=sxml->ReadElementDouble(ele,"rotateaxis","angle",true);
  string angunits=fun::StrLower(sxml->ReadElementStr(ele,"rotateaxis","anglesunits"));
  if(angunits=="radians")rotate=rotate*TODEG;
  else if(angunits!="degrees")sxml->ErrReadElement(ele,"rotate",false,"The value anglesunits must be \"degrees\" or \"radians\"."); 
  TiXmlElement* rot=ele->FirstChildElement("rotateaxis");
  tdouble3 pt1=sxml->ReadElementDouble3(rot,"point1");
  tdouble3 pt2=sxml->ReadElementDouble3(rot,"point2");
  const JMatrix4d m=JMatrix4d::MatrixRot(rotate,pt1,pt2);
  //-Adds config information about rotation.
  ConfigInfo.push_back(fun::PrintStr("  rotate: angle:%g axis:%s",rotate,fun::Double3gRangeStr(pt1,pt2).c_str()));
  return(m);
}

//==============================================================================
/// Creates points in a box.
//==============================================================================
void JSphInOutPoints::Create3d_Box(const JXml *sxml,TiXmlElement* ele){
  if(Count)Run_Exceptioon("Only one description zone is allowed for inlet/outlet points.");
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

  //-Applies rotation to inlet definition.
  double rotate=sxml->ReadElementDouble(ele,"rotateaxis","angle",true);
  if(rotate){
    const JMatrix4d m=ReadRotate3D(sxml,ele);
    tdouble3 pdir=pt0+dir;
    pt0=m.MulPoint(pt0);
    ptx=m.MulPoint(ptx);
    pty=m.MulPoint(pty);
    ptz=m.MulPoint(ptz);
    pdir=m.MulPoint(pdir);
    if(dir!=TDouble3(0))dir=pdir-pt0;
  }
  //-Updates Direction.
  if(dir!=TDouble3(0)){
    dir=fgeo::VecUnitary(dir);
    if(Direction!=TDouble3(0) && Direction!=dir)sxml->ErrReadElement(ele,"direction",false,"Direction is already defined.");
    Direction=dir;
  }
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
    const tdouble3 inimove=(Direction*InitialMove);
    for(unsigned pz=0;pz<npz;pz++)for(unsigned py=0;py<npy;py++)for(unsigned px=0;px<npx;px++){
      const tdouble3 ps=pt0+(vx*px)+(vy*py)+(vz*pz);
      Points[Count+p]=ps+inimove;
      p++;
    }
    if(p!=npt || Count+p>Size)Run_Exceptioon("Error calculating number of points.");
    Count+=npt;

    //-Compute domain limits of zone.
    {
      tdouble3 vxx,vyy;
      if(npz==1){      vxx=vx*(npx); vyy=vy*(npy); } //-Adds one dp.
      else if(npy==1){ vxx=vx*(npx); vyy=vz*(npz); } //-Adds one dp.
      else if(npx==1){ vxx=vy*(npy); vyy=vz*(npz); } //-Adds one dp.
      else Run_Exceptioon("Configuration is invalid to calculate domain.");
      const tdouble3 dir1=Direction*(Dp/2);
      const unsigned nfar=2*(Layers)+1;
      PtDom[0]=pt0+dir1-(fgeo::VecUnitary(vxx)*(Dp/2))-(fgeo::VecUnitary(vyy)*(Dp/2));
      PtDom[1]=PtDom[0]+vxx;
      PtDom[2]=PtDom[1]+vyy;
      PtDom[3]=PtDom[0]+vyy;
      PtDom[4]=PtDom[0]-(dir1*nfar);
      PtDom[5]=PtDom[1]-(dir1*nfar);
      PtDom[6]=PtDom[2]-(dir1*nfar);
      PtDom[7]=PtDom[3]-(dir1*nfar);
      PtDom[8]=(PtDom[0]+PtDom[2])/2;
      PtDom[9]=PtDom[8]+(Direction*(Dp*Layers));
    }
  }
}

//==============================================================================
/// Creates points in a circle.
//==============================================================================
void JSphInOutPoints::Create3d_Circle(const JXml *sxml,TiXmlElement* ele){
  if(Count)Run_Exceptioon("Only one description zone is allowed for inlet/outlet points.");
  //-Load basic data.
  const tdouble3 pt0=sxml->ReadElementDouble3(ele,"point");
  const double radius=sxml->ReadElementDouble(ele,"radius","v");
  tdouble3 pcen=pt0;
  //-Load direction.
  tdouble3 dir=TDouble3(0);
  if(sxml->ExistsElement(ele,"direction"))dir=sxml->ReadElementDouble3(ele,"direction");
  //-Applies rotation to inlet direction.
  JMatrix4d m;
  double rotate=sxml->ReadElementDouble(ele,"rotateaxis","angle",true);
  if(rotate){
    m=ReadRotate3D(sxml,ele);
    tdouble3 pdir=pt0+dir;
    pdir=m.MulPoint(pdir);
    pcen=m.MulPoint(pt0);
    if(dir!=TDouble3(0))dir=pdir-pcen;
  }
  //-Adds config information.
  ConfigInfo.push_back(fun::PrintStr("Circle: center:(%s) radius:%g",fun::Double3gStr(pcen).c_str(),radius));
  //-Updates Direction.
  if(dir!=TDouble3(0)){
    dir=fgeo::VecUnitary(dir);
    if(Direction!=TDouble3(0) && Direction!=dir)sxml->ErrReadElement(ele,"direction",false,"Direction is already defined.");
    Direction=dir;
  }
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
    const tdouble3 inimove=(Direction*InitialMove);
    unsigned p=0;
    for(unsigned cr=0;cr<=nra;cr++){
      const double ra=rasum*cr;
      const unsigned nang=max(1u,unsigned((TWOPI*ra+Dp/2)/Dp));
      const double angsum=TWOPI/nang; //-In radians.
      for(unsigned c=0;c<nang;c++){
        const tdouble3 ps=pt0+TDouble3(ra*cos(angsum*c),0,ra*sin(angsum*c));
        Points[Count+p]=(rotate? m.MulPoint(ps): ps)+inimove;
        p++;
      }
    }
    if(p!=npt || Count+p>Size)Run_Exceptioon("Error calculating number of points.");
    Count+=npt;

    //-Compute domain limits of zone.
    {
      const double ra=rasum*nra;
      tdouble3 vxx=TDouble3(ra*cos(0),0,ra*sin(0))*2;
      tdouble3 vyy=TDouble3(ra*cos(PIHALF),0,ra*sin(PIHALF))*2;
      PtDom[0]=pt0-(vxx/2)-(vyy/2)-(fgeo::VecUnitary(vxx)*(Dp/2))-(fgeo::VecUnitary(vyy)*(Dp/2));
      vxx=vxx+(fgeo::VecUnitary(vxx)*Dp);//-Adds one dp.
      vyy=vyy+(fgeo::VecUnitary(vyy)*Dp);//-Adds one dp.
      PtDom[1]=PtDom[0]+vxx;
      PtDom[2]=PtDom[1]+vyy;
      PtDom[3]=PtDom[0]+vyy;
      if(rotate)for(unsigned c=0;c<4;c++)PtDom[c]=m.MulPoint(PtDom[c]);
      const tdouble3 dir1=Direction*(Dp/2);
      const unsigned nfar=2*(Layers)+1;
      for(unsigned c=0;c<4;c++)PtDom[c]=PtDom[c]+dir1;
      PtDom[4]=PtDom[0]-(dir1*nfar);
      PtDom[5]=PtDom[1]-(dir1*nfar);
      PtDom[6]=PtDom[2]-(dir1*nfar);
      PtDom[7]=PtDom[3]-(dir1*nfar);
      PtDom[8]=(PtDom[0]+PtDom[2])/2;
      PtDom[9]=PtDom[8]+(Direction*(Dp*Layers));
    }
  }
}

//==============================================================================
/// Reads definition of inlet points in the XML node and creates points.
//==============================================================================
void JSphInOutPoints::CreatePoints(const JXml *sxml,TiXmlElement* lis
  ,const JDsPartsInit *partsdata)
{
  string xmlrow=sxml->ErrGetFileRow(lis);
  TiXmlElement* ele=lis->FirstChildElement();
  while(ele){
    string cmd=ele->Value();
    if(cmd.length()&&cmd[0]!='_'){
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
void JSphInOutPoints::ComputeDomainLimits(tdouble3 &posmin,tdouble3 &posmax)const{
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
/// Compute domain limits of zone from inout points (for 2D and particles-3D).
//==============================================================================
void JSphInOutPoints::ComputeDomainFromPoints(){
  //-Calculates minimum and maximum position of inout points. 
  tdouble3 pmin=TDouble3(DBL_MAX),pmax=TDouble3(-DBL_MAX);
  ComputeDomainLimits(pmin,pmax);
  const tdouble3 dir1=Direction*(Dp/2); //-Vector direction.
  const unsigned nfar=2*(Layers)+1;
  if(Simulate2D){
    const tdouble3 dir2=fgeo::VecUnitary(pmax-pmin)*(Dp/2); //-Vector normal direction.
    tdouble3 pnearmin=pmin+dir1-dir2;
    tdouble3 pnearmax=pmax+dir1+dir2;
    tdouble3 pfarmin=pnearmin-(dir1*nfar);//-Increases one layer and adds other dp for borders.
    tdouble3 pfarmax=pnearmax-(dir1*nfar);//-Increases one layer and adds other dp for borders.
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
    byte dir=0;
         if(Direction.z<Direction.x && Direction.z<Direction.y)dir=1;//-Bottom.
    else if(Direction.z>Direction.x && Direction.z>Direction.y)dir=2;//-Top.
    else if(Direction.x<Direction.y && Direction.x<Direction.z)dir=3;//-Left.
    else if(Direction.x>Direction.y && Direction.x>Direction.z)dir=4;//-Right.
    else if(Direction.y<Direction.x && Direction.y<Direction.z)dir=5;//-Front.
    else if(Direction.y>Direction.x && Direction.y>Direction.z)dir=6;//-Back.
    else Run_Exceptioon("The direction is invalid.");
    const double dph=Dp/2;
    tdouble3 p0=TDouble3(pmin.x-dph,pmin.y-dph,pmin.z-dph);
    double sx=pmax.x-pmin.x+Dp;
    double sy=pmax.y-pmin.y+Dp;
    double sz=pmax.z-pmin.z+Dp;
    if(dir==3){//-Left.
      sx=dph*nfar;
    }
    if(dir==4){//-Right.
      sx=dph*nfar;
      p0.x=pmin.x+dph-sx;
    }
    if(dir==5){//-Front.
      sy=dph*nfar;
    } 
    if(dir==6){//-Back.
      sy=dph*nfar;
      p0.y=pmin.y+dph-sy;
    } 
    if(dir==1){//-Bottom.
      sz=dph*nfar;
    } 
    if(dir==2){//-Top.
      sz=dph*nfar;
      p0.z=pmin.z+dph-sz;
    } 
    //-Stores domain points.
    PtDom[0]=p0;
    PtDom[1]=PtDom[0]+TDouble3(sx,0 ,0 );
    PtDom[2]=PtDom[1]+TDouble3(0 ,0 ,sz);
    PtDom[3]=PtDom[0]+TDouble3(0 ,0 ,sz);
    PtDom[4]=PtDom[0]+TDouble3(0 ,sy,0 );
    PtDom[5]=PtDom[1]+TDouble3(0 ,sy,0 );
    PtDom[6]=PtDom[2]+TDouble3(0 ,sy,0 );
    PtDom[7]=PtDom[3]+TDouble3(0 ,sy,0 );
    //-Stores normal points.
    if(dir==3)PtDom[8]=(PtDom[7]+PtDom[0])/2;//-Left.
    if(dir==4)PtDom[8]=(PtDom[6]+PtDom[1])/2;//-Right.
    if(dir==5)PtDom[8]=(PtDom[2]+PtDom[0])/2;//-Front.
    if(dir==6)PtDom[8]=(PtDom[6]+PtDom[4])/2;//-Back.
    if(dir==1)PtDom[8]=(PtDom[5]+PtDom[0])/2;//-Bottom.
    if(dir==2)PtDom[8]=(PtDom[6]+PtDom[3])/2;//-Top.
    PtDom[9]=PtDom[8]+(Direction*(Dp*Layers));
  }
}

//==============================================================================
/// Checks direction and position of points in simulation domain.
//==============================================================================
void JSphInOutPoints::CheckPoints(const std::string &xmlrow){
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
    tfloat3 *pos=new tfloat3[np];
    byte *layer=new byte[np];
    byte *outside=new byte[np];
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
    JVtkLib::SaveVtkData(file,arrays,"Pos");
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
void JSphInOutPoints::GetConfig(std::vector<std::string> &lines)const{
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
void JSphInOutPoints::GetPtDomain(std::vector<tdouble3> &ptdom)const{
  ptdom.clear();
  for(unsigned p=0;p<10;p++)ptdom.push_back(PtDom[p]);
}


