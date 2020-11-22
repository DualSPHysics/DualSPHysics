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

using std::string;

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
  List.clear();
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JDsDamping::LoadXml(const JXml *sxml,const std::string &place){
  List.clear();
  TiXmlNode* node=sxml->GetNodeSimple(place,false);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JDsDamping::ReadXml(const JXml *sxml,TiXmlElement* lis){
  //-Loads damping zones.
  TiXmlElement* ele=lis->FirstChildElement("dampingzone"); 
  while(ele){
    if(sxml->CheckElementActive(ele)){
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
        //-Obtains minimum and maximum Z.
        da.domzmin=sxml->ReadElementDouble(ele,"domain","zmin");
        da.domzmax=sxml->ReadElementDouble(ele,"domain","zmax");
        //-Obtains limit points.
        tdouble2 vp[4];
        vp[0].x=sxml->ReadElementDouble(dom,"point1","x");
        vp[0].y=sxml->ReadElementDouble(dom,"point1","y");
        vp[1].x=sxml->ReadElementDouble(dom,"point2","x");
        vp[1].y=sxml->ReadElementDouble(dom,"point2","y");
        vp[2].x=sxml->ReadElementDouble(dom,"point3","x");
        vp[2].y=sxml->ReadElementDouble(dom,"point3","y");
        vp[3].x=sxml->ReadElementDouble(dom,"point4","x");
        vp[3].y=sxml->ReadElementDouble(dom,"point4","y");
        da.dompt0=vp[0];
        da.dompt1=vp[1];
        da.dompt2=vp[2];
        da.dompt3=vp[3];
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
          tdouble2 pt=vp[c]; vp[c]=vp[c2]; vp[c2]=pt;
        }
        //:for(unsigned c=1;c<4;c++)printf("++ pt=(%s) - pt=(%s) ang:%g\n",fun::Double3Str(vp[0]).c_str(),fun::Double3Str(vp[c]).c_str(),angles[c]);

        //-Calculates lateral planes.
        //-Calcula planos laterales.
        const tdouble3 vp0=TDouble3(vp[0].x,vp[0].y,0);
        const tdouble3 vp1=TDouble3(vp[1].x,vp[1].y,0);
        const tdouble3 vp2=TDouble3(vp[2].x,vp[2].y,0);
        const tdouble3 vp3=TDouble3(vp[3].x,vp[3].y,0);
        da.dompla0=fgeo::Plane3Pt(vp0,vp1,TDouble3(vp1.x,vp1.y,1));
        da.dompla1=fgeo::Plane3Pt(vp1,vp2,TDouble3(vp2.x,vp2.y,1));
        da.dompla2=fgeo::Plane3Pt(vp2,vp3,TDouble3(vp3.x,vp3.y,1));
        da.dompla3=fgeo::Plane3Pt(vp3,vp0,TDouble3(vp0.x,vp0.y,1));
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
      //-Processes data entry.
      //-Procesa entrada de datos.
      {
        tdouble3 pt=da.limitmin;
        tdouble3 vec=da.limitmax-da.limitmin;
        da.dist=(float)fgeo::PointDist(vec);
        da.plane=fgeo::PlanePtVec(pt,vec);
      }
      List.push_back(da);
    }
    ele=ele->NextSiblingElement("dampingzone");
  }
  SaveVtkConfig(Dp);
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JDsDamping::VisuConfig(std::string txhead,std::string txfoot){
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
/// Saves VTK file with scheme of configuration.
//==============================================================================
void JDsDamping::SaveVtkConfig(double dp)const{
  const double sizequad=dp*16;
  const double sizedir=dp*4;
  JVtkLib sh;
  for(unsigned c=0;c<GetCount();c++){
    const int cv=int(c);
    const StDamping &zo=List[c];
    const tdouble3 ps=zo.limitmin;
    const tdouble3 ve=fgeo::VecUnitary(zo.limitmax-ps);
    //-Adds limit quad.
    sh.AddShapeQuad(zo.limitmin,ve,sizequad,cv);
    sh.AddShapeQuadWire(zo.limitmin,ve,sizequad,cv);
    //-Adds limitmax quad.
    sh.AddShapeQuad(zo.limitmax,ve,sizequad/3,cv);
    sh.AddShapeQuadWire(zo.limitmax,ve,sizequad/3,cv);
    //-Adds overlimit quad.
    const tdouble3 pt=zo.limitmax+(ve*double(zo.overlimit));
    sh.AddShapeQuad(pt,ve,sizequad/3,cv);
    sh.AddShapeQuadWire(pt,ve,sizequad/3,cv);
    //-Adds direction line.
    sh.AddShapeLine(zo.limitmin,pt,cv);
    if(zo.usedomain){
      const tdouble3 p0=TDouble3(zo.dompt0.x,zo.dompt0.y,zo.domzmin);
      const tdouble3 p1=TDouble3(zo.dompt1.x,zo.dompt1.y,zo.domzmin);
      const tdouble3 p2=TDouble3(zo.dompt2.x,zo.dompt2.y,zo.domzmin);
      const tdouble3 p3=TDouble3(zo.dompt3.x,zo.dompt3.y,zo.domzmin);
      const tdouble3 q0=TDouble3(zo.dompt0.x,zo.dompt0.y,zo.domzmax);
      const tdouble3 q1=TDouble3(zo.dompt1.x,zo.dompt1.y,zo.domzmax);
      const tdouble3 q2=TDouble3(zo.dompt2.x,zo.dompt2.y,zo.domzmax);
      const tdouble3 q3=TDouble3(zo.dompt3.x,zo.dompt3.y,zo.domzmax);
      sh.AddShapeQuad(p0,p1,p2,p3,cv); //-Bottom.
      sh.AddShapeQuad(q3,q2,q1,q0,cv); //-Top.
      sh.AddShapeQuad(p0,p1,q1,q0,cv);
      sh.AddShapeQuad(p1,p2,q2,q1,cv);
      sh.AddShapeQuad(p2,p3,q3,q2,cv);
      sh.AddShapeQuad(p3,p0,q0,q3,cv);
    }
  }
  if(GetCount()){
    string filevtk=AppInfo.GetDirOut()+"CfgDamping_Scheme.vtk";
    sh.SaveShapeVtk(filevtk,(GetCount()>1? "num": ""));
    Log->AddFileInfo(filevtk,"Saves VTK file with Damping configurations.");
  }
}

//==============================================================================
/// Returns the information of a block of particles.
/// Devuelve la informacion de un bloque de particulas.
//==============================================================================
const JDsDamping::StDamping* JDsDamping::GetDampingZone(unsigned c)const{
  if(c>=GetCount())Run_Exceptioon("The requested damping zone is not valid.");
  return(&(List[c]));
}

//==============================================================================
/// Applies Damping to indicated particles without domain limits.
/// Aplica Damping a las particulas indicadas sin limites de dominio.
//==============================================================================
void JDsDamping::ComputeDamping(const JDsDamping::StDamping &da,double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  const tplane3d plane=da.plane;
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
      double vdis=fgeo::PlanePoint(plane,ps);
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
void JDsDamping::ComputeDampingPla(const JDsDamping::StDamping &da,double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  const double zmin=da.domzmin;
  const double zmax=da.domzmax;
  const tplane3d pla0=da.dompla0;
  const tplane3d pla1=da.dompla1;
  const tplane3d pla2=da.dompla2;
  const tplane3d pla3=da.dompla3;
  const tplane3d plane=da.plane;
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
      double vdis=fgeo::PlanePoint(plane,ps);
      if(0<vdis && vdis<=dist+over){
        if(ps.z>=zmin && ps.z<=zmax && fgeo::PlanePoint(pla0,ps)<=0 && fgeo::PlanePoint(pla1,ps)<=0 && fgeo::PlanePoint(pla2,ps)<=0 && fgeo::PlanePoint(pla3,ps)<=0){
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
void JDsDamping::ComputeDamping(double timestep,double dt,unsigned n,unsigned pini
  ,const tdouble3 *pos,const typecode *code,tfloat4 *velrhop)const
{
  for(unsigned c=0;c<GetCount();c++){
    if(!List[c].usedomain)ComputeDamping(List[c],dt,n,pini,pos,code,velrhop);
    else ComputeDampingPla(List[c],dt,n,pini,pos,code,velrhop);
  }
}


