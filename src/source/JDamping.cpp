//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JDamping.cpp \brief Implements the class \ref JDamping.

#include "JDamping.h"
#include "JLog2.h"
#include "JXml.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include <cfloat>
#include <algorithm>

using std::string;

//##############################################################################
//# JDamping
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDamping::JDamping(JLog2* log):Log(log){
  ClassName="JDamping";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDamping::~JDamping(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDamping::Reset(){
  CellOrder=ORDER_None;
  List.clear();
}

//==============================================================================
/// Configures object.
//==============================================================================
void JDamping::Config(TpCellOrder cellorder){
  CellOrder=cellorder;
  if(CellOrder!=ORDER_XYZ)RunException("ConfigOrder","Only order XYZ is valid for now...");
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JDamping::LoadXml(JXml *sxml,const std::string &place){
  List.clear();
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)RunException("LoadXml",std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JDamping::ReadXml(JXml *sxml,TiXmlElement* lis){
  const char met[]="ReadXml";
  //-Loads damping zones.
  TiXmlElement* ele=lis->FirstChildElement("dampingzone"); 
  while(ele){
    StDamping da;
    memset(&da,0,sizeof(StDamping));
    da.limitmin=sxml->ReadElementDouble3(ele,"limitmin");
    da.limitmax=sxml->ReadElementDouble3(ele,"limitmax");
    da.overlimit=sxml->ReadElementFloat(ele,"overlimit","value");
    da.redumax=sxml->ReadElementFloat(ele,"redumax","value",true,10);
    da.factorxyz=TFloat3(1);
    if(sxml->ExistsElement(ele,"factorxyz")){
      TiXmlElement* elec=sxml->GetFirstElement(ele,"factorxyz",true);
      da.factorxyz.x=sxml->GetAttributeFloat(elec,"x",true,1.f);
      da.factorxyz.y=sxml->GetAttributeFloat(elec,"y",true,1.f);
      da.factorxyz.z=sxml->GetAttributeFloat(elec,"z",true,1.f);
      da.factorxyz=MinValues(da.factorxyz,TFloat3(1.f));
      da.factorxyz=MaxValues(da.factorxyz,TFloat3(0));
    }
    //-Loads domain limits.
    da.usedomain=false;
    TiXmlElement* dom=ele->FirstChildElement("domain");
    if(dom){
      da.usedomain=true;
      tdouble3 vp[4];
      double h=sxml->ReadElementFloat(ele,"domain","height");
      sxml->ReadArrayDouble3(dom,"point",vp,4);
      //-Obtains minimum and maximum Z.
      double zmin=vp[0].z;
      double zmax=vp[0].z;
      for(unsigned c=1;c<4;c++){
        if(zmin>vp[c].z)zmin=vp[c].z;
        if(zmax<vp[c].z)zmax=vp[c].z;
      }
      da.domzmin=zmin;
      da.domzmax=zmax+h;
      //-Sorts points 1-3 depending on the angle with the point 0.
      //-Ordena puntos 1-3 en funcion del angulo con el punto 0.
      double angles[4];
      for(unsigned c=1;c<4;c++){
        angles[c]=atan2((vp[0].y-vp[c].y),(vp[0].x-vp[c].x))*TODEG; 
        if(angles[c]<0)angles[c]+=360.;
        //:printf(" pt=(%s) - pt_%d=(%s) ang:%g\n",fun::Double3Str(vp[0]).c_str(),c,fun::Double3Str(vp[c]).c_str(),angles[c]);
      }
      for(unsigned c=1;c<3;c++)for(unsigned c2=c+1;c2<4;c2++)if(angles[c]>angles[c2]){
        double aux=angles[c]; angles[c]=angles[c2]; angles[c2]=aux;
        tdouble3 pt=vp[c]; vp[c]=vp[c2]; vp[c2]=pt;
      }
      //:for(unsigned c=1;c<4;c++)printf("++ pt=(%s) - pt=(%s) ang:%g\n",fun::Double3Str(vp[0]).c_str(),fun::Double3Str(vp[c]).c_str(),angles[c]);

      //-Calculates lateral planes.
      //-Calcula planos laterales.
      for(unsigned c=0;c<4;c++)vp[c].z=0;
      da.dompla0=fmath::Plane3Pt(vp[0],vp[1],TDouble3(vp[1].x,vp[1].y,1));
      da.dompla1=fmath::Plane3Pt(vp[1],vp[2],TDouble3(vp[2].x,vp[2].y,1));
      da.dompla2=fmath::Plane3Pt(vp[2],vp[3],TDouble3(vp[3].x,vp[3].y,1));
      da.dompla3=fmath::Plane3Pt(vp[3],vp[0],TDouble3(vp[0].x,vp[0].y,1));
      //:tdouble3 pt=sxml->ReadElementDouble3(ele,"pt");
      //:printf("++> pt=(%s)\n",fun::Double3Str(pt).c_str());
      //:printf("++> pla0:%f\n",fmath::PointPlane(da.dompla0,pt));
      //:printf("++> pla1:%f\n",fmath::PointPlane(da.dompla1,pt));
      //:printf("++> pla2:%f\n",fmath::PointPlane(da.dompla2,pt));
      //:printf("++> pla3:%f\n",fmath::PointPlane(da.dompla3,pt));
      //:bool inside=(fmath::PointPlane(da.dompla0,pt)<=0 && fmath::PointPlane(da.dompla1,pt)<=0 && fmath::PointPlane(da.dompla2,pt)<=0 && fmath::PointPlane(da.dompla3,pt)<=0);
      //:if(inside)printf("++> DENTRO\n"); else printf("++> fuera\n");
      //:exit(1);
    }
    //-Processes data entry.
    //-Procesa entrada de datos.
    {
      tdouble3 pt=da.limitmin;
      tdouble3 vec=da.limitmax-da.limitmin;
      da.dist=(float)fmath::DistPoint(vec);
      da.plane=fmath::PlanePtVec(pt,vec);
    }
    List.push_back(da);
    ele=ele->NextSiblingElement("dampingzone");
  }
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JDamping::VisuConfig(std::string txhead,std::string txfoot){
  if(!txhead.empty())Log->Print(txhead);
  for(unsigned c=0;c<GetCount();c++){
    const StDamping* da=GetDampingZone(c);
    Log->Printf("Damping zone_%u",c);
    Log->Printf("  LimitPoints: %s overlimit:%f",fun::Double3gRangeStr(da->limitmin,da->limitmax).c_str(),da->overlimit);
    Log->Printf("  LimitDist..: %g",da->dist);
    Log->Printf("  ReduMax....: %g",da->redumax);
    Log->Printf("  Factorxyz..: (%s)",fun::Float3gStr(da->factorxyz).c_str());

  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Returns the information of a block of particles.
/// Devuelve la informacion de un bloque de particulas.
//==============================================================================
const JDamping::StDamping* JDamping::GetDampingZone(unsigned c)const{
  if(c>=GetCount())RunException("GetDampingZone","The requested damping zone is not valid.");
  return(&(List[c]));
}

//==============================================================================
/// Applies Damping to indicated particles without domain limits.
/// Aplica Damping a las particulas indicadas sin limites de dominio.
//==============================================================================
void JDamping::ComputeDamping(const JDamping::StDamping &da,double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  const tdouble4 plane=da.plane;
  const float dist=da.dist;
  const float over=da.overlimit;
  const tfloat3 factorxyz=da.factorxyz;
  const float redumax=da.redumax;
  const int inp=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(inp>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<inp;p++){
    const unsigned p1=pini+unsigned(p);
    bool ok=true;
    if(code){//-Ignores floating and periodic particles. | Descarta particulas floating o periodicas.
      const typecode cod=code[p1];
      ok=(CODE_IsNormal(cod) && CODE_IsFluid(cod));
    }
    if(ok){
      const tdouble3 ps=pos[p1];
      double vdis=fmath::PointPlane(plane,ps);
      if(0<vdis && vdis<=dist+over){
        const double fdis=(vdis>=dist? 1.: vdis/dist);
        const double redudt=dt*(fdis*fdis)*redumax;
        double redudtx=(1.-redudt*factorxyz.x);
        double redudty=(1.-redudt*factorxyz.y);
        double redudtz=(1.-redudt*factorxyz.z);
        redudtx=(redudtx<0? 0.: redudtx);
        redudty=(redudty<0? 0.: redudty);
        redudtz=(redudtz<0? 0.: redudtz);
        velrhop[p1].x=float(redudtx*velrhop[p1].x);
        velrhop[p1].y=float(redudty*velrhop[p1].y);
        velrhop[p1].z=float(redudtz*velrhop[p1].z);
      }
    }
  }
}

//==============================================================================
/// Applies Damping to the indicated particles within the domain delimited by 
/// domplaX planes.
///
/// Aplica Damping a las particulas indicadas dentro del dominio delimitado por 
/// los planos domplaX.
//==============================================================================
void JDamping::ComputeDampingPla(const JDamping::StDamping &da,double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  const double zmin=da.domzmin;
  const double zmax=da.domzmax;
  const tdouble4 pla0=da.dompla0;
  const tdouble4 pla1=da.dompla1;
  const tdouble4 pla2=da.dompla2;
  const tdouble4 pla3=da.dompla3;
  const tdouble4 plane=da.plane;
  const float dist=da.dist;
  const float over=da.overlimit;
  const tfloat3 factorxyz=da.factorxyz;
  const float redumax=da.redumax;
  const int inp=int(n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (static) if(inp>OMP_LIMIT_COMPUTEMEDIUM)
  #endif
  for(int p=0;p<inp;p++){
    const unsigned p1=pini+unsigned(p);
    bool ok=true;
    if(code){//-Ignores floating and periodic particles. | Descarta particulas floating o periodicas.
      const typecode cod=code[p1];
      ok=(CODE_IsNormal(cod) && CODE_IsFluid(cod));
    }
    if(ok){
      const tdouble3 ps=pos[p1];
      //-Check if it is within the domain. | Comprueba si esta dentro del dominio.
      double vdis=fmath::PointPlane(plane,ps);
      if(0<vdis && vdis<=dist+over){
        if(ps.z>=zmin && ps.z<=zmax && fmath::PointPlane(pla0,ps)<=0 && fmath::PointPlane(pla1,ps)<=0 && fmath::PointPlane(pla2,ps)<=0 && fmath::PointPlane(pla3,ps)<=0){
          const double fdis=(vdis>=dist? 1.: vdis/dist);
          const double redudt=dt*(fdis*fdis)*redumax;
          double redudtx=(1.-redudt*factorxyz.x);
          double redudty=(1.-redudt*factorxyz.y);
          double redudtz=(1.-redudt*factorxyz.z);
          redudtx=(redudtx<0? 0.: redudtx);
          redudty=(redudty<0? 0.: redudty);
          redudtz=(redudtz<0? 0.: redudtz);
          velrhop[p1].x=float(redudtx*velrhop[p1].x);
          velrhop[p1].y=float(redudty*velrhop[p1].y);
          velrhop[p1].z=float(redudtz*velrhop[p1].z);
        }
      }
    }
  }
}

//==============================================================================
/// Applies Damping to the indicated particles.
/// Aplica Damping a las particulas indicadas.
//==============================================================================
//:#include "JSaveCsv.h"
//:#include "JFormatFiles2.h"
void JDamping::ComputeDamping(double timestep,double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  for(unsigned c=0;c<GetCount();c++){
    if(!List[c].usedomain)ComputeDamping(List[c],dt,n,pini,pos,code,velrhop);
    else ComputeDampingPla(List[c],dt,n,pini,pos,code,velrhop);
  }
  /*:
  //if(0){//dg
  //  tdouble3 pmin=TDouble3(0);
  //  tdouble3 pmax=TDouble3(10,10,0);
  //  double dp=0.01;
  //  unsigned npx=unsigned((pmax.x-pmin.x)/dp)+1;
  //  unsigned npy=unsigned((pmax.y-pmin.y)/dp)+1;
  //  unsigned np=npx*npy;
  //  tdouble3 *pos=new tdouble3[np];
  //  tfloat4 *velrhop=new tfloat4[np];
  //  unsigned cp=0;
  //  for(unsigned cy=0;cy<npy;cy++)for(unsigned cx=0;cx<npx;cx++){
  //    pos[cp]=TDouble3(dp*cx,dp*cy,0);
  //    velrhop[cp]=TFloat4(1);
  //    cp++;
  //  }
  //  for(unsigned cd=0;cd<GetCount();cd++){
  //    if(!List[cd].usedomain)ComputeDamping(List[cd],1.f/List[cd].redumax,np,0,pos,NULL,velrhop);
  //    else ComputeDampingPla(List[cd],1.f/List[cd].redumax,np,0,pos,NULL,velrhop);
  //  }
  //  //-Define campos a grabar.
  //  tfloat3 *posf=new tfloat3[np];
  //  tfloat3 *vel=new tfloat3[np];
  //  tfloat3 *velx=new tfloat3[np];
  //  tfloat3 *vely=new tfloat3[np];
  //  tfloat3 *velz=new tfloat3[np];
  //  for(unsigned c=0;c<np;c++){
  //    tfloat3 ps=ToTFloat3(pos[c]);
  //    ps.z=velrhop[c].x;
  //    posf[c]=ps;
  //    const tfloat3 v=TFloat3(velrhop[c].x,velrhop[c].y,velrhop[c].z);
  //    vel[c]=v;
  //    velx[c]=TFloat3(0,0,v.x);
  //    vely[c]=TFloat3(0,0,v.y);
  //    velz[c]=TFloat3(0,0,v.z);
  //  }
  //  JFormatFiles2::StScalarData fields[8];
  //  unsigned nfields=0;
  //  if(vel){   fields[nfields]=JFormatFiles2::DefineField("Vel" ,JFormatFiles2::Float32,3,vel);   nfields++; }
  //  if(velx){  fields[nfields]=JFormatFiles2::DefineField("Velx",JFormatFiles2::Float32,3,velx);  nfields++; }
  //  if(vely){  fields[nfields]=JFormatFiles2::DefineField("Vely",JFormatFiles2::Float32,3,vely);  nfields++; }
  //  if(velz){  fields[nfields]=JFormatFiles2::DefineField("Velz",JFormatFiles2::Float32,3,velz);  nfields++; }
  //  JFormatFiles2::SaveVtk("__Damping_00.vtk",np,posf,nfields,fields);
  //  exit(0);
  //}

  //if(1){//dg
  //  double x1=3,x2=7,ix=0.1;
  //  unsigned np=unsigned((x2-x1)/ix+1);
  //  tdouble3 *vpos=new tdouble3[np];
  //  tfloat4 *vvel=new tfloat4[np];
  //  for(unsigned p=0;p<np;p++){
  //    vpos[p]=TDouble3(x1+ix*p,0,0);
  //    vvel[p]=TFloat4(10.f,1.f,0.1f,1000.f);
  //  }
  //  if(1){
  //    JSaveCsv sv(Log->GetDirOut()+"_Damping_0.csv",false);
  //    sv.AddHead("X;VelX;VelY;VelZ");
  //    ComputeDamping(List[c],1,np,0,vpos,NULL,vvel);
  //    for(unsigned p=0;p<np;p++){
  //      sv.AddValuesf("%g;%f;%f;%f",vpos[p].x,vvel[p].x,vvel[p].y,vvel[p].z);
  //      sv.AddEndl();
  //    }
  //  }
  //  JSaveCsv sv(Log->GetDirOut()+"_Damping_1.csv",false);
  //  string txhead="time";
  //  for(unsigned p=0;p<np;p++)txhead=txhead+fun::PrintStr(";x_%.1f",vpos[p].x);
  //  sv.AddHead(txhead);
  //  for(unsigned p=0;p<np;p++){
  //    vpos[p]=TDouble3(x1+ix*p,0,0);
  //    vvel[p]=TFloat4(10.f,1.f,0.1f,1000.f);
  //  }
  //  double tmax=10,dt=0.01;
  //  for(double t=0;t<tmax+dt/2;t+=dt){
  //    sv.AddValue(t);
  //    for(unsigned p=0;p<np;p++)sv.AddValue(vvel[p].y);
  //    sv.AddEndl();
  //    ComputeDamping(List[c],dt,np,0,vpos,NULL,vvel);
  //  }
  //  exit(0);
  //}
  :*/
}


