//HEAD_DSPH
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

/// \file JSphShifting.cpp \brief Implements the class \ref JSphShifting.

#include "JSphShifting.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "JXml.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JDataArrays.h"
#include "JVtkLib.h"
#ifdef _WITHGPU
  #include "JSphShifting_ker.h"
  //#include "FunctionsCuda.h"
#endif
#include <cfloat>
#include <algorithm>

using std::string;

//##############################################################################
//# JSphShiftingZone
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphShiftingZone::JSphShiftingZone(unsigned id,const tdouble3 &posref,const tdouble3 &vx
  ,const tdouble3 &vy,const tdouble3 &vz):Id(id)
{
  ClassName="JSphShifting";
  Reset();
  PosRef=posref;
  Vecx=vx;
  Vecy=vy;
  Vecz=vz;
  PrepareZone();
}
//==============================================================================
/// Destructor.
//==============================================================================
JSphShiftingZone::~JSphShiftingZone(){
  DestructorActive=true;
  Reset();
}
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphShiftingZone::Reset(){
  PosRef=Vecx=Vecy=Vecz=TDouble3(0);
  UsePosMax=true;
  PosMax=TDouble3(0);
  DomPlax=DomPlay=DomPlaz=TPlane3d(0);
  DomPladis=TDouble3(0);
}
//==============================================================================
/// Prepare zone for calculation.
//==============================================================================
void JSphShiftingZone::PrepareZone(){
  UsePosMax=(!Vecx.y && !Vecx.z && !Vecy.x && !Vecy.z && !Vecz.x && !Vecz.y);
  if(UsePosMax)PosMax=PosRef+TDouble3(Vecx.x,Vecy.y,Vecz.z);
  else fgeo::PlanesDomain(PosRef,Vecx,Vecy,Vecz,DomPlax,DomPlay,DomPlaz,DomPladis);
}


//##############################################################################
//# JSphShifting
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphShifting::JSphShifting(bool simulate2d,double dp,float kernelh)
  :Log(AppInfo.LogPtr()),Simulate2D(simulate2d),Dp(dp),KernelH(kernelh)
{
  ClassName="JSphShifting";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphShifting::~JSphShifting(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphShifting::Reset(){
  ConfigBasic(SHIFT_None,0,0);
  ZonesXml=ZonesPosmax=0;
  for(unsigned c=0;c<GetCount();c++)delete Zones[c];
  Zones.clear();
}

//==============================================================================
/// Returns value of Shifting mode in text format.
/// Devuelve el valor de Shifting mode en texto.
//==============================================================================
std::string JSphShifting::GetShiftingModeStr()const{
  string tx;
  if(ShiftMode==SHIFT_None)tx="None";
  else if(ShiftMode==SHIFT_NoBound)tx="NoBound";
  else if(ShiftMode==SHIFT_NoFixed)tx="NoFixed";
  else if(ShiftMode==SHIFT_Full)tx="Full";
  else tx="???";
  return(tx);
}

//==============================================================================
/// Explicit configuration from execution parameters.
//==============================================================================
void JSphShifting::ConfigBasic(TpShifting shiftmode,float shiftcoef,float shifttfs){
  ShiftMode=shiftmode;
  ShiftCoef=shiftcoef;
  ShiftTFS=shifttfs;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSphShifting::LoadXml(const JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Returns matrix for rotation in 3D.
//==============================================================================
JMatrix4d JSphShifting::ReadRotate3D(const JXml *sxml,TiXmlElement* ele){
  double rotate=sxml->ReadElementDouble(ele,"rotateaxis","angle",true);
  string angunits=fun::StrLower(sxml->ReadElementStr(ele,"rotateaxis","anglesunits"));
  if(angunits=="radians")rotate=rotate*TODEG;
  else if(angunits!="degrees")sxml->ErrReadElement(ele,"rotate",false,"The value anglesunits must be \"degrees\" or \"radians\"."); 
  TiXmlElement* rot=ele->FirstChildElement("rotateaxis");
  sxml->CheckElementNames(rot,true,"point1 point2");
  tdouble3 pt1=sxml->ReadElementDouble3(rot,"point1");
  tdouble3 pt2=sxml->ReadElementDouble3(rot,"point2");
  const JMatrix4d m=JMatrix4d::MatrixRot(rotate,pt1,pt2);
  return(m);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JSphShifting::ReadXml(const JXml *sxml,TiXmlElement* lis){
  const int smode=sxml->ReadElementInt(lis,"mode","value");
  switch(smode){
    case 0:  ShiftMode=SHIFT_None;     break;
    case 1:  ShiftMode=SHIFT_NoBound;  break;
    case 2:  ShiftMode=SHIFT_NoFixed;  break;
    case 3:  ShiftMode=SHIFT_Full;     break;
    default: Run_ExceptioonFile("Shifting mode is not valid.",sxml->ErrGetFileRow(lis,"mode"));
  }
  //-Loads other configurations when shifting is activated.
  if(ShiftMode!=SHIFT_None){
    sxml->CheckElementNames(lis,true,"mode coefficient fsthreshold *shiftingzone");
    ShiftCoef=sxml->ReadElementFloat(lis,"coefficient","value",true,-2.f);  
    ShiftTFS =sxml->ReadElementFloat(lis,"fsthreshold","value",true,0);  
    //-Loads shifting zones.
    TiXmlElement* ele=lis->FirstChildElement("shiftingzone"); 
    while(ele){
      if(sxml->CheckElementActive(ele)){
        JSphShiftingZone* pzo=NULL;
        //StShifting zo;
        //memset(&zo,0,sizeof(StShifting));
        //zo.useposmax=false;
        tdouble3 posref=sxml->ReadElementDouble3(ele,"posref");
        tdouble3 vecx=TDouble3(0),vecy=TDouble3(0),vecz=TDouble3(0);
        if(sxml->ExistsElement(ele,"size")){
          sxml->CheckElementNames(ele,true,"posref size rotateaxis");
          const tdouble3 size=sxml->ReadElementDouble3(ele,"size");
          vecx=TDouble3(size.x,0,0);
          vecy=TDouble3(0,size.y,0);
          vecz=TDouble3(0,0,size.z);
        }
        else{
          sxml->CheckElementNames(ele,true,"posref vecx vecy vecz rotateaxis");
          vecx=sxml->ReadElementDouble3(ele,"vecx");
          vecy=sxml->ReadElementDouble3(ele,"vecy");
          vecz=sxml->ReadElementDouble3(ele,"vecz");
        }
        //-Applies rotation to domain definition.
        if(sxml->ReadElementDouble(ele,"rotateaxis","angle",true)){
          const JMatrix4d m=ReadRotate3D(sxml,ele);
          const tdouble3 px=m.MulPoint(posref+vecx);
          const tdouble3 py=m.MulPoint(posref+vecy);
          const tdouble3 pz=m.MulPoint(posref+vecz);
          posref=m.MulPoint(posref);
          vecx=px-posref;
          vecy=py-posref;
          vecz=pz-posref;
        }
        AddZone(true,posref,vecx,vecy,vecz);
      }
      ele=ele->NextSiblingElement("shiftingzone");
    }
  }
}

//==============================================================================
/// Add configurantion for a new shifting zone.
//==============================================================================
void JSphShifting::AddZone(bool fromxml,const tdouble3 &posref
  ,const tdouble3 &vx,const tdouble3 &vy,const tdouble3 &vz)
{
  const unsigned id=GetCount();
  JSphShiftingZone* zo=new JSphShiftingZone(id,posref,vx,vy,vz);
  if(zo->GetUsePosMax()){
    Zones.insert(Zones.begin()+ZonesPosmax,zo);
    ZonesPosmax++;
  }
  else Zones.push_back(zo);
  if(fromxml)ZonesXml++;
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSphShifting::VisuConfig(std::string txhead,std::string txfoot){
  if(!txhead.empty())Log->Print(txhead);
  Log->Print(fun::VarStr("Shifting",GetShiftingModeStr()));
  Log->Print(fun::VarStr("  ShiftCoef",ShiftCoef));
  if(ShiftTFS)Log->Print(fun::VarStr("  ShiftTFS",ShiftTFS));
  else Log->Printf("  ShiftTFS=%g (disabled)",ShiftTFS);
  if(!GetCount())Log->Print(fun::VarStr("  ShiftDomain","Full"));
  else{
    string tx;
    if(ZonesXml)tx=fun::PrintStr("%u from XML file",ZonesXml);
    if(GetCount()>ZonesXml){
      if(!tx.empty())tx=tx+"and ";
      tx=tx+fun::PrintStr("%u from in/out configuration",GetCount()-ZonesXml);
    }
    Log->Print(fun::VarStr("  ShiftDomain",string("Zones (")+tx+")"));
    for(unsigned c=0;c<GetCount();c++){
      const JSphShiftingZone* zo=Zones[c];
      if(zo->GetUsePosMax())Log->Printf("  Zone_%d  %s",c,fun::Double3gRangeStr(zo->GetPosMin(),zo->GetPosMax()).c_str());
      else{
        tdouble3 p=zo->GetPosMin();
        Log->Printf("  Zone_%d  posref:(%g,%g,%g)",zo->Id,p.x,p.y,p.z);
        p=zo->GetVecx();  Log->Printf("    vx:(%g,%g,%g)",p.x,p.y,p.z);
        p=zo->GetVecy();  Log->Printf("    vy:(%g,%g,%g)",p.x,p.y,p.z);
        p=zo->GetVecz();  Log->Printf("    vz:(%g,%g,%g)",p.x,p.y,p.z);
      }
    }
  }
  if(!txfoot.empty())Log->Print(txfoot);
  SaveVtkConfig();
}

//==============================================================================
/// Returns configuration information in one string line.
//==============================================================================
std::string JSphShifting::GetConfigInfo()const{
  string ret;
  if(ShiftMode!=SHIFT_None){
    ret=string("Shifting(")+GetShiftingModeStr()
      +fun::PrintStr(",%g,%g,%s)",ShiftCoef,ShiftTFS,(GetCount()? "Zones": "Full"));
  }
  return(ret);
}

//==============================================================================
/// Saves VTK file with zones configuration.
//==============================================================================
void JSphShifting::SaveVtkConfig()const{
  const unsigned nz=GetCount();
  if(nz){
    JVtkLib sh;
    const unsigned nz=GetCount();
    for(unsigned c=0;c<nz;c++){
      const JSphShiftingZone* zo=Zones[c];
      sh.AddShapeBox(zo->GetPosMin(),zo->GetVecx(),zo->GetVecy(),zo->GetVecz(),c);
    }
    const string filevtk=AppInfo.GetDirOut()+"CfgShifting_Zones.vtk";
    sh.SaveShapeVtk(filevtk,"ZoneId");
    if(nz==ZonesXml || !ZonesXml)Log->AddFileInfo(filevtk,"Saves VTK file with Shifting zones.");
  }
}

//==============================================================================
/// Select particles for shifting according min-max configuration
/// Selecciona particulas para shifting segun configuracion min-max.
//==============================================================================
template<bool first,bool dbl> void JSphShifting::InitCpuPosMax(unsigned n,unsigned pini
  ,const tdouble3& pmin1,const tdouble3& pmax1,const tdouble3& pmin2,const tdouble3& pmax2
  ,const tdouble3* pos,tfloat4* shiftposfs)const
{
  const int ppini=int(pini),ppfin=ppini+int(n),npf=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ppini;p<ppfin;p++){
    const tdouble3 ps=pos[p];
    if(fgeo::PointInMinMax(ps,pmin1,pmax1) || (dbl && fgeo::PointInMinMax(ps,pmin2,pmax2))){
      shiftposfs[p]=TFloat4(0);
    }
    else if(first)shiftposfs[p]=TFloat4(FLT_MAX);
  }
}

//==============================================================================
/// Select particles for shifting according planes configuration
/// Selecciona particulas para shifting segun configuracion de planos.
//==============================================================================
template<bool first,bool dbl> void JSphShifting::InitCpuPlanes(unsigned n,unsigned pini
  ,const tplane3d& plax1,const tplane3d& play1,const tplane3d& plaz1,const tdouble3& pladis1
  ,const tplane3d& plax2,const tplane3d& play2,const tplane3d& plaz2,const tdouble3& pladis2
  ,const tdouble3* pos,tfloat4* shiftposfs)const
{
  const int ppini=int(pini),ppfin=ppini+int(n),npf=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ppini;p<ppfin;p++){
    const tdouble3 ps=pos[p];
    if(fgeo::PlanesDomainCheck(ps,plax1,play1,plaz1,pladis1) || (dbl && fgeo::PlanesDomainCheck(ps,plax2,play2,plaz2,pladis2))){
      shiftposfs[p]=TFloat4(0);
    }
    else if(first)shiftposfs[p]=TFloat4(FLT_MAX);
  }
}

//==============================================================================
/// Select particles for shifting and initialize shiftposfs[].
/// Selecciona particulas para shifting e inicializa shiftposfs[].
//==============================================================================
void JSphShifting::InitCpu(unsigned n,unsigned pini,const tdouble3* pos,tfloat4* shiftposfs)const{
  const unsigned nz=GetCount();
  if(!nz)memset(shiftposfs+pini,0,sizeof(tfloat4)*n);   //shiftposfs[]=0
  else{
    //-Zones defined by position min-max.
    unsigned cz=0;
    unsigned nz=ZonesPosmax;
    while(nz){
      const JSphShiftingZone* zo1=Zones[cz];
      const JSphShiftingZone* zo2=(nz>1? Zones[cz+1]: NULL);
      if(!cz){ const bool first=true; 
        if(nz>1)InitCpuPosMax<first,true >(n,pini,zo1->GetPosMin(),zo1->GetPosMax(),zo2->GetPosMin(),zo2->GetPosMax(),pos,shiftposfs);
        else    InitCpuPosMax<first,false>(n,pini,zo1->GetPosMin(),zo1->GetPosMax(),TDouble3(0)     ,TDouble3(0)     ,pos,shiftposfs);
      }else{   const bool first=false; 
        if(nz>1)InitCpuPosMax<first,true >(n,pini,zo1->GetPosMin(),zo1->GetPosMax(),zo2->GetPosMin(),zo2->GetPosMax(),pos,shiftposfs);
        else    InitCpuPosMax<first,false>(n,pini,zo1->GetPosMin(),zo1->GetPosMax(),TDouble3(0)     ,TDouble3(0)     ,pos,shiftposfs);
      }
      if(nz>1){ nz-=2; cz+=2; }
      else    { nz-=1; cz+=1; }
    }
    //-Zones defined by planes.
    nz=GetCount()-ZonesPosmax;
    while(nz){
      const JSphShiftingZone* zo1=Zones[cz];
      const JSphShiftingZone* zo2=(nz>1? Zones[cz+1]: NULL);
      if(!cz){ const bool first=true; 
        if(nz>1)InitCpuPlanes<first,true >(n,pini,zo1->GetDomPlax(),zo1->GetDomPlay(),zo1->GetDomPlaz(),zo1->GetDomPladis(),zo2->GetDomPlax(),zo2->GetDomPlay(),zo2->GetDomPlaz(),zo2->GetDomPladis(),pos,shiftposfs);
        else    InitCpuPlanes<first,false>(n,pini,zo1->GetDomPlax(),zo1->GetDomPlay(),zo1->GetDomPlaz(),zo1->GetDomPladis(),TPlane3d(0)      ,TPlane3d(0)      ,TPlane3d(0)      ,TDouble3(0)        ,pos,shiftposfs);
      }else{   const bool first=false; 
        if(nz>1)InitCpuPlanes<first,true >(n,pini,zo1->GetDomPlax(),zo1->GetDomPlay(),zo1->GetDomPlaz(),zo1->GetDomPladis(),zo2->GetDomPlax(),zo2->GetDomPlay(),zo2->GetDomPlaz(),zo2->GetDomPladis(),pos,shiftposfs);
        else    InitCpuPlanes<first,false>(n,pini,zo1->GetDomPlax(),zo1->GetDomPlay(),zo1->GetDomPlaz(),zo1->GetDomPladis(),TPlane3d(0)      ,TPlane3d(0)      ,TPlane3d(0)      ,TDouble3(0)        ,pos,shiftposfs);
      }
      if(nz>1){ nz-=2; cz+=2; }
      else    { nz-=1; cz+=1; }
    }
  }
}

//==============================================================================
/// Calculate final Shifting for particles' position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphShifting::RunCpu(unsigned n,unsigned pini,double dt,const tfloat4* velrhop
  ,tfloat4* shiftposfs)const
{
  const double coefumagn=dt*ShiftCoef*KernelH;
  const double coeftfs=(Simulate2D? 2.0: 3.0)-ShiftTFS;
  const float maxdist=float(Dp*0.1); //-Max shifting distance permitted (recommended).
  const int ppini=int(pini),ppfin=pini+int(n),npf=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(npf>OMP_LIMIT_COMPUTELIGHT)
  #endif
  for(int p=ppini;p<ppfin;p++){
    const tfloat4 rs=shiftposfs[p];
    if(rs.x!=FLT_MAX){
      const double vx=double(velrhop[p].x);
      const double vy=double(velrhop[p].y);
      const double vz=double(velrhop[p].z);
      double umagn=coefumagn*sqrt(vx*vx+vy*vy+vz*vz);
      if(ShiftTFS){
        if(rs.w<ShiftTFS)umagn=0;
        else umagn*=(double(rs.w)-ShiftTFS)/coeftfs;
      }
      const float shiftdistx=float(double(rs.x)*umagn);
      const float shiftdisty=float(double(rs.y)*umagn);
      const float shiftdistz=float(double(rs.z)*umagn);
      shiftposfs[p].x=(shiftdistx<maxdist? shiftdistx: maxdist);
      shiftposfs[p].y=(shiftdisty<maxdist? shiftdisty: maxdist);
      shiftposfs[p].z=(shiftdistz<maxdist? shiftdistz: maxdist);
    }
    else shiftposfs[p]=TFloat4(0); //-Cancels shifting close to the boundaries. | Anula shifting por proximidad del contorno. 
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Select particles for shifting and initialize shiftposfs[].
/// Selecciona particulas para shifting e inicializa shiftposfs[].
//==============================================================================
void JSphShifting::InitGpu(unsigned n,unsigned pini,const double2* posxy,const double* posz
  ,float4* shiftposfs,cudaStream_t stm)const
{
  const unsigned nz=GetCount();
  if(!nz)cudaMemsetAsync(shiftposfs+pini,0,sizeof(float4)*n,stm);  //ShiftPosfsg[]=0
  else{
    //-Zones defined by position min-max.
    unsigned cz=0;
    unsigned nz=ZonesPosmax;
    while(nz){
      const JSphShiftingZone* zo1=Zones[cz];
      const JSphShiftingZone* zo2=(nz>1? Zones[cz+1]: NULL);
      const bool tfirst=(!cz);
      if(nz>1)cushift::InitGpuPosMax(tfirst,true ,n,pini,zo1->GetPosMin(),zo1->GetPosMax(),zo2->GetPosMin(),zo2->GetPosMax(),posxy,posz,shiftposfs,stm);
      else    cushift::InitGpuPosMax(tfirst,false,n,pini,zo1->GetPosMin(),zo1->GetPosMax(),TDouble3(0)     ,TDouble3(0)     ,posxy,posz,shiftposfs,stm);
      if(nz>1){ nz-=2; cz+=2; }
      else    { nz-=1; cz+=1; }
    }
    //-Zones defined by planes.
    nz=GetCount()-ZonesPosmax;
    while(nz){
      const JSphShiftingZone* zo1=Zones[cz];
      const JSphShiftingZone* zo2=(nz>1? Zones[cz+1]: NULL);
      const bool tfirst=(!cz);
      if(nz>1)cushift::InitGpuPlanes(tfirst,true ,n,pini,zo1->GetDomPlax(),zo1->GetDomPlay(),zo1->GetDomPlaz(),zo1->GetDomPladis(),zo2->GetDomPlax(),zo2->GetDomPlay(),zo2->GetDomPlaz(),zo2->GetDomPladis(),posxy,posz,shiftposfs,stm);
      else    cushift::InitGpuPlanes(tfirst,false,n,pini,zo1->GetDomPlax(),zo1->GetDomPlay(),zo1->GetDomPlaz(),zo1->GetDomPladis(),TPlane3d(0)      ,TPlane3d(0)      ,TPlane3d(0)      ,TDouble3(0)        ,posxy,posz,shiftposfs,stm);
      if(nz>1){ nz-=2; cz+=2; }
      else    { nz-=1; cz+=1; }
    }
  }
}

//==============================================================================
/// Calculate final Shifting for particles' position.
/// Calcula Shifting final para posicion de particulas.
//==============================================================================
void JSphShifting::RunGpu(unsigned n,unsigned pini,double dt,const float4* velrhop
  ,float4* shiftposfs,cudaStream_t stm)const
{
  const double coefumagn=dt*ShiftCoef*KernelH;
  const float coeftfs=(Simulate2D? 2.0f: 3.0f)-ShiftTFS;
  const float maxdist=float(Dp*0.1); //-Max shifting distance permitted (recommended).
  cushift::RunShifting(n,pini,dt,coefumagn,ShiftTFS,coeftfs,maxdist,velrhop,shiftposfs,stm);
}
#endif



