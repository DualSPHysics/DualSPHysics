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

/// \file JDsDamping.cpp \brief Implements the class \ref JDsDamping.

#include "JDsDamping.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "JXml.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JVtkLib.h"
#include <cfloat>
#include <algorithm>

#ifdef _WITHGPU
#include "JSphGpu_ker.h"
#include "FunctionsBasic_iker.h"
#endif

using std::string;

//##############################################################################
//# JDsDampingOp_Plane
//##############################################################################
//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsDampingOp_Plane::Reset(){
  JDsDampingOp::Reset();
  LimitMin=LimitMax=TDouble3(0);
  UseDomain=false;
  DomzMin=DomzMax=0;
  DomPt0=DomPt1=DomPt2=DomPt3=TDouble2(0);
  DomPla0=DomPla1=DomPla2=DomPla3=TPlane3d(0);
  Plane=TPlane3d(0);
  Dist=0;
}

//==============================================================================
/// Calculates lateral planes.
//==============================================================================
void JDsDampingOp_Plane::ComputeDomPlanes(){
  tdouble2 vp[4]{DomPt0,DomPt1,DomPt2,DomPt3};
  //-Sorts points 1-3 depending on the angle with the point 0.
  double angles[4];
  for(unsigned c=1;c<4;c++){
    angles[c]=atan2((vp[0].y-vp[c].y),(vp[0].x-vp[c].x))*TODEG; 
    if(angles[c]<0)angles[c]+=360.;
    //:printf(" pt=(%s) - pt_%d=(%s) ang:%g\n",fun::Double3Str(vp[0]).c_str(),c,fun::Double3Str(vp[c]).c_str(),angles[c]);
  }
  for(unsigned c=1;c<3;c++)for(unsigned c2=c+1;c2<4;c2++)if(angles[c]>angles[c2]){
    double aux=angles[c]; angles[c]=angles[c2]; angles[c2]=aux;
    tdouble2 pt=vp[c]; vp[c]=vp[c2]; vp[c2]=pt;
  }
  //:for(unsigned c=1;c<4;c++)printf("++ pt=(%s) - pt=(%s) ang:%g\n",fun::Double3Str(vp[0]).c_str(),fun::Double3Str(vp[c]).c_str(),angles[c]);

  //-Calculates lateral planes.
  const tdouble3 vp0=TDouble3(vp[0].x,vp[0].y,0);
  const tdouble3 vp1=TDouble3(vp[1].x,vp[1].y,0);
  const tdouble3 vp2=TDouble3(vp[2].x,vp[2].y,0);
  const tdouble3 vp3=TDouble3(vp[3].x,vp[3].y,0);
  DomPla0=fgeo::Plane3Pt(vp0,vp1,TDouble3(vp1.x,vp1.y,1));
  DomPla1=fgeo::Plane3Pt(vp1,vp2,TDouble3(vp2.x,vp2.y,1));
  DomPla2=fgeo::Plane3Pt(vp2,vp3,TDouble3(vp3.x,vp3.y,1));
  DomPla3=fgeo::Plane3Pt(vp3,vp0,TDouble3(vp0.x,vp0.y,1));
  //:tdouble3 pt=sxml->ReadElementDouble3(ele,"pt");
  //:printf("++> pt=(%s)\n",fun::Double3Str(pt).c_str());
  //:printf("++> pla0:%f\n",fgeo::PlanePoint(da.dompla0,pt));
  //:printf("++> pla1:%f\n",fgeo::PlanePoint(da.dompla1,pt));
  //:printf("++> pla2:%f\n",fgeo::PlanePoint(da.dompla2,pt));
  //:printf("++> pla3:%f\n",fgeo::PlanePoint(da.dompla3,pt));
  //:bool inside=(fgeo::PlanePoint(da.dompla0,pt)<=0 && fgeo::PlanePoint(da.dompla1,pt)<=0 && fgeo::PlanePoint(da.dompla2,pt)<=0 && fgeo::PlanePoint(da.dompla3,pt)<=0);
  //:if(inside)printf("++> DENTRO\n"); else printf("++> fuera\n");
  //:exit(1);
}

//==============================================================================
/// Reads particles information in xml format.
//==============================================================================
void JDsDampingOp_Plane::ReadXml(const JXml *sxml,TiXmlElement* ele){
  sxml->CheckElementNames(ele,true,"overlimit redumax factorxyz limitmin limitmax domain");
  //-General options.
  OverLimit=sxml->ReadElementFloat(ele,"overlimit","value");
  ReduMax=sxml->ReadElementFloat(ele,"redumax","value",true,10);
  Factorxyz.x=sxml->ReadElementFloat(ele,"factorxyz","x",true,1);
  Factorxyz.y=sxml->ReadElementFloat(ele,"factorxyz","y",true,1);
  Factorxyz.z=sxml->ReadElementFloat(ele,"factorxyz","z",true,1);
  //-Specific options.
  LimitMin=sxml->ReadElementDouble3(ele,"limitmin");
  LimitMax=sxml->ReadElementDouble3(ele,"limitmax");
  //-Loads domain limits.
  UseDomain=false;
  TiXmlElement* dom=ele->FirstChildElement("domain");
  if(dom){
    UseDomain=true;
    //-Obtains minimum and maximum Z.
    DomzMin=sxml->ReadElementDouble(ele,"domain","zmin");
    DomzMax=sxml->ReadElementDouble(ele,"domain","zmax");
    //-Obtains limit points.
    DomPt0=TDouble2(sxml->ReadElementDouble(dom,"point1","x"),sxml->ReadElementDouble(dom,"point1","y"));
    DomPt1=TDouble2(sxml->ReadElementDouble(dom,"point2","x"),sxml->ReadElementDouble(dom,"point2","y"));
    DomPt2=TDouble2(sxml->ReadElementDouble(dom,"point3","x"),sxml->ReadElementDouble(dom,"point3","y"));
    DomPt3=TDouble2(sxml->ReadElementDouble(dom,"point4","x"),sxml->ReadElementDouble(dom,"point4","y"));
    //-Calculates lateral planes.
    ComputeDomPlanes();
  }
  //-Processes input configuration.
  {
    const tdouble3 pt=LimitMin;
    const tdouble3 vec=LimitMax-LimitMin;
    Dist=float(fgeo::PointDist(vec));
    Plane=fgeo::PlanePtVec(pt,vec);
  }
}

//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsDampingOp_Plane::SaveVtkConfig(double dp,JVtkLib *sh)const{
  const double sizequad=dp*16;
  const double sizedir=dp*4;
  const int cv=int(Id);
  const tdouble3 ps=LimitMin;
  const tdouble3 ve=fgeo::VecUnitary(LimitMax-ps);
  //-Adds limit quad.
  sh->AddShapeQuad(LimitMin,ve,sizequad,cv);
  sh->AddShapeQuadWire(LimitMin,ve,sizequad,cv);
  //-Adds limitmax quad.
  sh->AddShapeQuad(LimitMax,ve,sizequad/3,cv);
  sh->AddShapeQuadWire(LimitMax,ve,sizequad/3,cv);
  //-Adds overlimit quad.
  const tdouble3 pt=LimitMax+(ve*double(OverLimit));
  sh->AddShapeQuad(pt,ve,sizequad/3,cv);
  sh->AddShapeQuadWire(pt,ve,sizequad/3,cv);
  //-Adds direction line.
  sh->AddShapeLine(LimitMin,pt,cv);
  if(UseDomain){
    const tdouble3 p0=TDouble3(DomPt0.x,DomPt0.y,DomzMin);
    const tdouble3 p1=TDouble3(DomPt1.x,DomPt1.y,DomzMin);
    const tdouble3 p2=TDouble3(DomPt2.x,DomPt2.y,DomzMin);
    const tdouble3 p3=TDouble3(DomPt3.x,DomPt3.y,DomzMin);
    const tdouble3 q0=TDouble3(DomPt0.x,DomPt0.y,DomzMax);
    const tdouble3 q1=TDouble3(DomPt1.x,DomPt1.y,DomzMax);
    const tdouble3 q2=TDouble3(DomPt2.x,DomPt2.y,DomzMax);
    const tdouble3 q3=TDouble3(DomPt3.x,DomPt3.y,DomzMax);
    sh->AddShapeQuad(p0,p1,p2,p3,cv); //-Bottom.
    sh->AddShapeQuad(q3,q2,q1,q0,cv); //-Top.
    sh->AddShapeQuad(p0,p1,q1,q0,cv);
    sh->AddShapeQuad(p1,p2,q2,q1,cv);
    sh->AddShapeQuad(p2,p3,q3,q2,cv);
    sh->AddShapeQuad(p3,p0,q0,q3,cv);
  }
}

//==============================================================================
/// Returns strings with configuration.
//==============================================================================
void JDsDampingOp_Plane::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back(fun::PrintStr("Damping zone_%u (%s): ",Id,GetNameType(Type).c_str()));
  lines.push_back(fun::PrintStr("  LimitPoints: %s overlimit:%f",fun::Double3gRangeStr(LimitMin,LimitMax).c_str(),OverLimit));
  lines.push_back(fun::PrintStr("  LimitDist..: %g",Dist));
  lines.push_back(fun::PrintStr("  ReduMax....: %g",ReduMax));
  lines.push_back(fun::PrintStr("  Factorxyz..: (%g,%g,%g)",Factorxyz.x,Factorxyz.y,Factorxyz.z));
}

//==============================================================================
/// Applies Damping to the indicated particles within the domain delimited by 
/// domplaX planes.
//==============================================================================
void JDsDampingOp_Plane::ComputeDampingCpu(double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  const bool usedomain=UseDomain;
  const double zmin=DomzMin;
  const double zmax=DomzMax;
  const tplane3d pla0=DomPla0;
  const tplane3d pla1=DomPla1;
  const tplane3d pla2=DomPla2;
  const tplane3d pla3=DomPla3;
  const tplane3d plane=Plane;
  const float dist=Dist;
  const float over=OverLimit;
  const tfloat3 factorxyz=Factorxyz;
  const float redumax=ReduMax;
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
      double vdis=fgeo::PlanePoint(plane,ps);
      if(0<vdis && vdis<=dist+over){
        if(!usedomain || (ps.z>=zmin && ps.z<=zmax && fgeo::PlanePoint(pla0,ps)<=0 && fgeo::PlanePoint(pla1,ps)<=0 && fgeo::PlanePoint(pla2,ps)<=0 && fgeo::PlanePoint(pla3,ps)<=0)){
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

#ifdef _WITHGPU
//==============================================================================
/// Applies Damping to the indicated particles within the domain delimited by 
/// domplaX planes.
//==============================================================================
void JDsDampingOp_Plane::ComputeDampingGpu(double dt,unsigned n,unsigned pini
  ,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop)const
{
  const float dist=Dist;
  const float over=OverLimit;
  const float3 factorxyz=Float3(Factorxyz);
  const float redumax=ReduMax;
  const double4 plane=Double4(Plane);
  if(!UseDomain){
    cusph::ComputeDampingPlane(dt,plane,dist,over,factorxyz,redumax,n,pini,posxy,posz,code,velrhop);
  }
  else{
    const double zmin=DomzMin;
    const double zmax=DomzMax;
    const double4 pla0=Double4(DomPla0);
    const double4 pla1=Double4(DomPla1);
    const double4 pla2=Double4(DomPla2);
    const double4 pla3=Double4(DomPla3);
    cusph::ComputeDampingPlaneDom(dt,plane,dist,over,factorxyz,redumax,zmin,zmax
      ,pla0,pla1,pla2,pla3,n,pini,posxy,posz,code,velrhop);
  }
}
#endif


//##############################################################################
//# JDsDamping
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsDamping::JDsDamping(double dp):Log(AppInfo.LogPtr()),Dp(dp){
  ClassName="JDsDamping";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsDamping::~JDsDamping(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsDamping::Reset(){
  for(unsigned c=0;c<Count();c++)delete List[c];
  List.clear();
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JDsDamping::LoadXml(const JXml *sxml,const std::string &place){
  Reset();
  TiXmlNode* node=sxml->GetNodeSimple(place,false);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JDsDamping::ReadXml(const JXml *sxml,TiXmlElement* lis){
  //-Loads damping zones.
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    string cmd=ele->Value();
    if(cmd.length() && cmd[0]!='_' && sxml->CheckElementActive(ele)){
      //printf("-----------> [%s]\n",cmd.c_str());
      if(cmd=="dampingzone"){  JDsDampingOp_Plane *dmp=new JDsDampingOp_Plane(Count(),sxml,ele); List.push_back(dmp);  }
      else sxml->ErrReadElement(ele,cmd,false);
    }
    ele=ele->NextSiblingElement();
  }
  //-Saves VTK file with scheme of configuration.
  SaveVtkConfig(Dp);
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JDsDamping::VisuConfig(std::string txhead,std::string txfoot){
  if(!txhead.empty())Log->Print(txhead);
  for(unsigned c=0;c<Count();c++){
    std::vector<std::string> lines;
    List[c]->GetConfig(lines);
    Log->Print(lines);
    //for(unsigned i=0;i<unsigned(lines.size());i++)Log->Print(string("  ")+lines[i]);
  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsDamping::SaveVtkConfig(double dp)const{
  if(Count()){
    JVtkLib sh;
    for(unsigned c=0;c<Count();c++)List[c]->SaveVtkConfig(dp,&sh);
    string filevtk=AppInfo.GetDirOut()+"CfgDamping_Scheme.vtk";
    sh.SaveShapeVtk(filevtk,(Count()>1? "num": ""));
    Log->AddFileInfo(filevtk,"Saves VTK file with Damping configurations.");
  }
}

//==============================================================================
/// Applies Damping to the indicated particles.
//==============================================================================
//:#include "JSaveCsv.h"
void JDsDamping::ComputeDampingCpu(double timestep,double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  for(unsigned c=0;c<Count();c++){
    List[c]->ComputeDampingCpu(dt,n,pini,pos,code,velrhop);
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Applies Damping to the indicated particles.
//==============================================================================
void JDsDamping::ComputeDampingGpu(double timestep,double dt,unsigned n,unsigned pini
  ,const double2 *posxy,const double *posz,const typecode *code,float4 *velrhop)const
{
  for(unsigned c=0;c<Count();c++){
    List[c]->ComputeDampingGpu(dt,n,pini,posxy,posz,code,velrhop);
  }
}
#endif

