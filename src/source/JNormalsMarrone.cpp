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

/// \file JNormalsMarrone.cpp \brief Implements the class \ref JNormalsMarrone.

#include "JNormalsMarrone.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JPartDataBi4.h"
#include "JPartNormalData.h"
#include "JDataArrays.h"
#include "JVtkLib.h"

#include <cstring>
#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JNormalsMarrone
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JNormalsMarrone::JNormalsMarrone(){
  ClassName="JNormalsMarrone";
  PartPos=NULL;  PartNor=NULL;
  NormalBegin=NULL;  
  Normals=NULL;  NormalsDist=NULL;
  OutVecs=NULL;  OutVecsDist=NULL;
  ExternalData=false;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JNormalsMarrone::~JNormalsMarrone(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JNormalsMarrone::Reset(){
  CaseDir="";
  CaseName="";
  DirOut="";
  ResetParts();
  ResetNormals();
  ExternalData=false;
}

//==============================================================================
/// Initialization of variables about particles.
//==============================================================================
void JNormalsMarrone::ResetParts(){
  Data2D=false;
  Data2DPosY=0;
  H=0;
  Dp=0;
  AllocParts(0);
}

//==============================================================================
/// Initialization of variables about normals.
//==============================================================================
void JNormalsMarrone::ResetNormals(){
  Dist=0;
  AllocNormals(0);
}

//==============================================================================
/// Allocates or frees memory for particles.
//==============================================================================
void JNormalsMarrone::AllocParts(unsigned size){
  if(ExternalData){
    PartPos=NULL;
    PartNor=NULL;
  }
  else{
    delete[] PartPos; PartPos=NULL;
    delete[] PartNor; PartNor=NULL;
    SizePart=size;
    if(SizePart){
      try{
        PartPos=new tdouble3[SizePart];
        PartNor=new tdouble3[SizePart];
      }
      catch(const std::bad_alloc){
        Run_Exceptioon("Could not allocate the requested memory.");
      }
    }
  } 
}

//==============================================================================
/// Allocates or frees memory for normals data.
//==============================================================================
void JNormalsMarrone::AllocNormals(unsigned size){
  if(ExternalData){
    NormalBegin=NULL;
    Normals=NULL;
    NormalsDist=NULL;
    OutVecs=NULL;
    OutVecsDist=NULL;
  }
  else{
    delete[] NormalBegin;  NormalBegin=NULL;
    delete[] Normals;      Normals=NULL;
    delete[] NormalsDist;  NormalsDist=NULL;
    delete[] OutVecs;      OutVecs=NULL;
    delete[] OutVecsDist;  OutVecsDist=NULL;
    SizeNor=size;
    if(SizePart){
      try{
        NormalBegin=new unsigned[SizePart+1];
        Normals    =new tdouble3[SizeNor];
        NormalsDist=new double  [SizeNor];
        OutVecs    =new tdouble3[SizeNor];
        OutVecsDist=new double  [SizeNor];
      }
      catch(const std::bad_alloc){
        Run_Exceptioon("Could not allocate the requested memory.");
      }
    } 
  }
}

//==============================================================================
/// Load boundary particles.
//==============================================================================
void JNormalsMarrone::LoadBoundParticles(){
  ResetParts();
  JPartDataBi4 pd;
  //-Loads file piece_0 and obtains configuration.
  //-Carga fichero piece_0 y obtiene configuracion.
  {
    const string file1=CaseDir+JPartDataBi4::GetFileNameCase(CaseName,0,1);
    if(fun::FileExists(file1))pd.LoadFileCase(CaseDir,CaseName,0,1);
    else if(fun::FileExists(CaseDir+JPartDataBi4::GetFileNameCase(CaseName,0,2)))pd.LoadFileCase(CaseDir,CaseName,0,2);
    else Run_ExceptioonFile("File of the particles was not found.",file1);
  }
  //-Obtains configuration. | Obtiene configuracion.
  const unsigned npiece=pd.GetNpiece();
  Data2D=pd.Get_Data2d();
  Data2DPosY=(Data2D? pd.Get_Data2dPosY(): 0);
  H=pd.Get_H();
  Dp=pd.Get_Dp();
  unsigned casenbound=unsigned(pd.Get_CaseNfixed()+pd.Get_CaseNmoving()+pd.Get_CaseNfloat());
  const bool possingle=pd.Get_PosSimple();
  if(!pd.Get_IdpSimple())Run_Exceptioon("Only Idp (32 bits) is valid at the moment.");
  //-Allocates memory.
  AllocParts(casenbound);
  for(unsigned p=0;p<SizePart;p++){
    PartPos[p]=TDouble3(DBL_MAX);
    PartNor[p]=TDouble3(0);
  }
  //-Loads particles.
  for(unsigned piece=0;piece<npiece;piece++){
    if(piece)pd.LoadFileCase(CaseDir,CaseName,piece,npiece);
    const unsigned npok=pd.Get_Npok();
    if(npok){
      unsigned *vidp=new unsigned[npok];
      pd.Get_Idp(npok,vidp);
      if(possingle){
        tfloat3 *vposf=new tfloat3[npok];
        pd.Get_Pos(npok,vposf);
        for(unsigned p=0;p<npok;p++){
          const unsigned idp=vidp[p];
          if(idp<casenbound)PartPos[idp]=ToTDouble3(vposf[p]);
        }
        delete[] vposf; vposf=NULL;
      }
      else{
        tdouble3 *vposd=new tdouble3[npok];
        pd.Get_Posd(npok,vposd);
        for(unsigned p=0;p<npok;p++){
          const unsigned idp=vidp[p];
          if(idp<casenbound)PartPos[idp]=vposd[p];
        }
        delete[] vposd; vposd=NULL;
      }
      delete[] vidp; vidp=NULL;
    }
  }
  //-In simulations 2D, if PosY is invalid then calculates starting from position of particles.
  if(Data2DPosY==DBL_MAX){
    if(!SizePart)Run_Exceptioon("Number of particles is invalid to calculates Y in 2D simulations.");
    Data2DPosY=PartPos[0].y;
  }
  //-Checks data of boundary particles.
  for(unsigned p=0;p<SizePart;p++)if(PartPos[p].x==DBL_MAX)Run_Exceptioon("Some postion of bound particles is invalid.");
}

//==============================================================================
/// Load normal data.
//==============================================================================
void JNormalsMarrone::LoadNormalData(){
  ResetNormals();
  JPartNormalData nd;
  nd.LoadFile(CaseDir+CaseName);
  if(nd.GetCountNormals()==0)Run_ExceptioonFile("Normal data is missing in NBI4 file.",GetNormalDataFile(CaseDir+CaseName));
  const unsigned nb=nd.GetNbound(); //-Number of particles with normal data.
  if(nb>SizePart)Run_ExceptioonFile("Normal data for more particles than expected in NBI4 file.",GetNormalDataFile(CaseDir+CaseName));
  if(nb<SizePart)SizePart=nb;
  Dist=nd.GetDist();
  AllocNormals(nd.GetCountNormals());
  memcpy(NormalBegin,nd.GetNormalBegin(),sizeof(unsigned)*(SizePart+1));
  memcpy(Normals    ,nd.GetNormals()    ,sizeof(tdouble3)*SizeNor);
  memcpy(NormalsDist,nd.GetNormalsDist(),sizeof(double)  *SizeNor);
  memcpy(OutVecs    ,nd.GetOutVecs()    ,sizeof(tdouble3)*SizeNor);
  memcpy(OutVecsDist,nd.GetOutVecsDist(),sizeof(double)  *SizeNor);
}

//==============================================================================
/// Saves particles with normals (for debug).
//==============================================================================
void JNormalsMarrone::SaveVtkNormalData(){
  const string file=DirOut+CaseName+"_NorData.vtk";
  //-Allocates memory.
  tfloat3  *vpos=new tfloat3[SizeNor];
  unsigned *vidp=new unsigned[SizeNor];
  byte    *vnum=new byte[SizeNor];
  tfloat3 *vnor=new tfloat3[SizeNor];
  tfloat3 *vout=new tfloat3[SizeNor];
  float   *vnordist=new float[SizeNor];
  float   *voutdist=new float[SizeNor];
  byte    *vnorsign=new byte[SizeNor];
  byte    *voutsign=new byte[SizeNor];
  //-Compute data.
  unsigned np=0;
  for(unsigned p=0;p<SizePart;p++){
    const tdouble3 ps=PartPos[p];
    const unsigned cini=NormalBegin[p];
    const unsigned cfin=NormalBegin[p+1];
    for(unsigned c=cini;c<cfin;c++){
      const tdouble3 n=Normals[c];
      const tdouble3 t=OutVecs[c];
      const double  nd=NormalsDist[c];
      const double  td=OutVecsDist[c];
      vpos[np]=ToTFloat3(ps);
      vidp[np]=p;
      vnum[np]=byte(c-cini);
      vnor[np]=ToTFloat3(TDouble3(n.x,n.y,n.z)*nd);
      vout[np]=ToTFloat3(TDouble3(t.x,t.y,t.z)*td);
      vnordist[np]=float(nd);
      voutdist[np]=float(td);
      vnorsign[np]=(nd>=0? 0: 1);
      voutsign[np]=(td>=0? 0: 1);
      np++;
    }
  }
  //-Saves VTK file and frees memory.
  JDataArrays arrays;
  arrays.AddArray("Pos",np,vpos,true);
  if(vidp)    arrays.AddArray("Idp"       ,np,vidp    ,true);
  if(vnum)    arrays.AddArray("NormalNum" ,np,vnum    ,true);
  if(vnor)    arrays.AddArray("Normal"    ,np,vnor    ,true);
  if(vout)    arrays.AddArray("OutVec"    ,np,vout    ,true);
  if(vnordist)arrays.AddArray("NormalDist",np,vnordist,true);
  if(voutdist)arrays.AddArray("OutVecDist",np,voutdist,true);
  if(vnorsign)arrays.AddArray("NormalSign",np,vnorsign,true);
  if(voutsign)arrays.AddArray("OutVecSign",np,voutsign,true);
  JVtkLib::SaveVtkData(file,arrays,"Pos");
  arrays.Reset();
}

//==============================================================================
/// Returns the file name with normal data.
//==============================================================================
std::string JNormalsMarrone::GetNormalDataFile(std::string casename){
  const string casedir=fun::GetDirWithSlash(fun::GetDirParent(casename));
  const string casenam=fun::GetWithoutExtension(fun::GetFile(casename));
  return(JPartNormalData::GetNormalDataFile(casedir+casenam));
}

//==============================================================================
/// Computes normals of case according configuration in XML file.
//==============================================================================
void JNormalsMarrone::RunCase(std::string casename,std::string dirout,bool savevtk){
  Reset();
  CaseDir=fun::GetDirWithSlash(fun::GetDirParent(casename));
  CaseName=fun::GetWithoutExtension(fun::GetFile(casename));
  DirOut=fun::GetDirWithSlash(!dirout.empty()? dirout: CaseDir);
  LoadBoundParticles();
  LoadNormalData();
  //SaveVtkNormalData(); //-For debug.
  ComputeNormalsMarrone();
  if(savevtk)SaveVtkNormalFinal(DirOut+CaseName+"_NorMarrone.vtk");
}

//==============================================================================
/// Computes normals of case according parameters.
//==============================================================================
void JNormalsMarrone::RunData(std::string casename,std::string dirout,bool savevtk
  ,bool data2d,double data2dposy,double h,double dp,unsigned sizepart,tdouble3 *partpos
  ,double dist,unsigned sizenor,unsigned *norbegin,tdouble3 *normals,double *normalsdist
  ,tdouble3 *outvecs,double *outvecsdist,tdouble3 *partnor)
{
  Reset();
  CaseDir=fun::GetDirWithSlash(fun::GetDirParent(casename));
  CaseName=fun::GetWithoutExtension(fun::GetFile(casename));
  DirOut=fun::GetDirWithSlash(!dirout.empty()? dirout: CaseDir);
  ExternalData=true;
  //-Particle data.
  Data2D=data2d;
  Data2DPosY=data2dposy;
  H=h;
  Dp=dp;
  SizePart=sizepart;
  PartPos=partpos;
  PartNor=partnor;
  //-Normal data.
  Dist=dist;
  SizeNor=sizenor;
  NormalBegin=norbegin;
  Normals=normals;
  NormalsDist=normalsdist;
  OutVecs=outvecs;
  OutVecsDist=outvecsdist;
  //SaveVtkNormalData(); //-For debug.
  ComputeNormalsMarrone();
  if(savevtk)SaveVtkNormalFinal(DirOut+CaseName+"_NorMarrone.vtk");
}

//==============================================================================
/// Computes normals for Marrone BC.
//==============================================================================
unsigned JNormalsMarrone::MaxNormalsByPart()const{
  unsigned nmax=0;
  for(unsigned p=0;p<SizePart;p++){
    const tdouble3 ps=PartPos[p];
    const unsigned cini=NormalBegin[p];
    const unsigned cfin=NormalBegin[p+1];
    const unsigned n=cfin-cini;
    nmax=max(nmax,n);
  }
  return(nmax);
}

//==============================================================================
/// Computes normals for Marrone BC.
//==============================================================================
void JNormalsMarrone::ComputeNormalsMarrone(){
  const double tolerance=Dp*0.0001;
  const unsigned nmax=MaxNormalsByPart();
  tdouble3* vnor=new tdouble3[nmax];
  double*   dnor=new double  [nmax];
  tdouble3* vout=new tdouble3[nmax];
  double*   dout=new double  [nmax];
  tdouble3* plim=new tdouble3[nmax];
  unsigned* cmin=new unsigned[nmax];
  //printf("AAA_000 nmax:%d\n",nmax);
  for(unsigned p=0;p<SizePart;p++){
    const tdouble3 ps=PartPos[p];
    const unsigned cini=NormalBegin[p];
    const unsigned cfin=NormalBegin[p+1];
    tdouble3 nor=TDouble3(0);
    bool nordone=false;
    //-For debug.
    const bool DG=false;//(p==0); (p==56)
    if(DG)for(unsigned c=cini;c<cfin;c++){
      const tdouble3 n=Normals[c];
      const tdouble3 t=OutVecs[c];
      const double  nd=NormalsDist[c];
      const double  td=OutVecsDist[c];
      printf("pos[%u]:(%g,%g,%g) n:(%g,%g,%g):%g  t:(%g,%g,%g):%g\n",p,ps.x,ps.y,ps.z,n.x,n.y,n.z,nd,t.x,t.y,t.z,td);
    }
    //-Busca distancia minima cuando la proyeccion de la posicion esta dentro del shape.
    double norindmin=DBL_MAX;
    double noroutdmin=DBL_MAX;
    for(unsigned c=cini;c<cfin;c++)if(NormalsDist[c]>0){
      if(fun::IsGtEqual(OutVecsDist[c],0,tolerance)){//-The projection position is in the shape.
        if(norindmin>NormalsDist[c])norindmin=NormalsDist[c];
      }
      else if(noroutdmin>NormalsDist[c])noroutdmin=NormalsDist[c];
    }
    const bool norin=(norindmin<=Dist);//-Maximum distance allows is Dist (KernelSize).
    if(DG)printf("-> norindmin:%g  noroutdmin:%g\n",norindmin,noroutdmin);
    //-Almacena informacion de normales distintas.
    unsigned nsel=0;
    for(unsigned c=cini;c<cfin;c++)if(NormalsDist[c]>0 && (!norin || fun::IsGtEqual(OutVecsDist[c],0,tolerance))){
      unsigned cn=0;
      for(;cn<nsel && !fun::IsEqual(vnor[cn],Normals[c],tolerance);cn++);
      if(cn<nsel && !norin){//-Si no hay normales in shape se queda con la que proporcione un punto de limite mas cercano.
        tdouble3 ptlim=ps+(Normals[c]*NormalsDist[c]);   //-Projection position in shape.
        ptlim=ptlim+(OutVecs[c]*OutVecsDist[c]);         //-Corner position.
        double d=fgeo::PointDist(ptlim-ps);              //-Distance to corner position.
        if(d<fgeo::PointDist(plim[cn]-ps)){
          vnor[cn]=Normals[c];
          dnor[cn]=NormalsDist[c];
          vout[cn]=OutVecs[c];
          dout[cn]=OutVecsDist[c];
          plim[cn]=ptlim;
        }
      }
      if(cn>=nsel){
        vnor[nsel]=Normals[c];
        dnor[nsel]=NormalsDist[c];
        vout[nsel]=OutVecs[c];
        dout[nsel]=OutVecsDist[c];
        plim[nsel]=ps+(Normals[c]*NormalsDist[c]);                //-Projection position in shape.
        if(!norin)plim[nsel]=plim[nsel]+(OutVecs[c]*OutVecsDist[c]); //-Corner position.
        nsel++;
      }
    }
    //-Shows information for debug.
    if(DG)for(unsigned c=0;c<nsel;c++){
      const tdouble3 n=vnor[c];
      const tdouble3 t=vout[c];
      const double  nd=dnor[c];
      const double  td=dout[c];
      const tdouble3 pp=plim[c];
      printf("Sel-> pos[%u]:(%g,%g,%g) n:(%g,%g,%g):%g  t:(%g,%g,%g):%g  plim:(%g,%g,%g)\n",p,ps.x,ps.y,ps.z,n.x,n.y,n.z,nd,t.x,t.y,t.z,td,pp.x,pp.y,pp.z);
    }
    //-Almacena en cmin[] las normales a distancia minima.
    unsigned ndismin=0;
    //for(unsigned cn=0;cn<nsel;cn++)if((norin && fun::IsEqual(norindmin,dnor[cn],tolerance)) || (!norin && fun::IsEqual(noroutdmin,dnor[cn],tolerance))){
    for(unsigned cn=0;cn<nsel;cn++)if(norin && fun::IsEqual(norindmin,dnor[cn],tolerance)){
      cmin[ndismin]=cn;
      ndismin++;
      if(DG){
        const tdouble3 n=vnor[cn];
        const tdouble3 t=vout[cn];
        const double  nd=dnor[cn];
        const double  td=dout[cn];
        const tdouble3 pp=plim[cn];
        printf("SelNear-> %u> n:(%g,%g,%g):%g  t:(%g,%g,%g):%g  plim:(%g,%g,%g)\n",cn,n.x,n.y,n.z,nd,t.x,t.y,t.z,td,pp.x,pp.y,pp.z);
      }
    }

    //-Compute final normal.
    if(norin){//-The projection position is in the shape.
      //-Calcula normal cuando hay dos o mas posiciones de interseccion a la misma distancia.
      if(ndismin>1){
        nor=TDouble3(0);
        for(unsigned cn=0;cn<ndismin;cn++)nor=nor+(vnor[cmin[cn]]*dnor[cmin[cn]]);
        if(ndismin){
          nor=nor/double(ndismin);
          const double snor=fgeo::PointDist(nor);
          //printf("==> pos[%u]:(%g,%g,%g)\n",p,ps.x,ps.y,ps.z);
          double snormin=DBL_MAX;
          for(unsigned cn=0;cn<ndismin;cn++){
            const unsigned cc=cmin[cn];
            const double d=dnor[cc]*dnor[cc]/fgeo::ProductScalar(nor,vnor[cc]*dnor[cc])*snor;//-Distancia entre pos y interseccion de nor con plano de vnor[cn].
            if(d>=0 && d<snormin)snormin=d;
          }
          if(snormin==DBL_MAX)snormin=0;
          nor=fgeo::VecUnitary(nor)*snormin;
        }
        nordone=true;
        if(fgeo::PointDist(nor)>Dist*1.5)nordone=false; //-Disables normals by projection greater than 1.5* maximum distance. 
        if(DG)printf("FinalNear-> [%u]> ndismin:%d  n:(%g,%g,%g):%g\n",p,ndismin,nor.x,nor.y,nor.z,fgeo::PointDist(nor));
      }
      //-Calcula normal cuando solo hay una posicion de interseccion o habiendo varias no dieron un resultado valido.
      if(!nordone && ndismin>=1){
        unsigned cn=cmin[0];
        nor=plim[cn]-ps;
        nordone=true;
      }
    }
    else if(nsel){//-The projection position is not in any shape.
      //-Compute minimum distance to border point.
      double dmin=DBL_MAX;
      unsigned nmin=0;
      tdouble3 pmin=TDouble3(0);
      bool diff=false;
      for(unsigned c=0;c<nsel;c++){
        const double d=fgeo::PointDist(plim[c]-ps);
        if(fun::IsEqual(d,dmin,tolerance)){
          nmin++;
          if(!fun::IsEqual(pmin,plim[c],tolerance))diff=true;
        }
        else if(dmin>d){ dmin=d; nmin=1; pmin=plim[c]; }
      }
      if(DG)printf("==> dmin:%g  nmin:%u  diff:%s\n",dmin,nmin,(diff? "true": "false"));
      if(nmin==3 && diff){//-Busca interseccion de planos.
        tplane3d pla[3];
        unsigned npla=0;
        for(unsigned c=0;c<nsel && npla<3;c++){
          const double d=fgeo::PointDist(plim[c]-ps);
          if(fun::IsEqual(d,dmin,tolerance)){
            pla[npla]=fgeo::PlanePtVec(ps+(vnor[c]*dnor[c]),vnor[c]);
            if(DG){
              const tdouble3 pt=ps+(vnor[c]*dnor[c]);
              const tdouble3 n=vnor[c];
              //printf("==> c:%u  pla:%u  n:(%g,%g,%g)  pt:(%g,%g,%g)\n",c,npla,n.x,n.y,n.z,pt.x,pt.y,pt.z);
            }
            npla++;
          }
        }
        if(npla==3){
          const tdouble3 pt=fgeo::PlanesIntersec(pla[0],pla[1],pla[2]);
          //printf("==> pla-intersect  pt:(%g,%g,%g)\n",pt.x,pt.y,pt.z);
          nor=pt-ps;
          nordone=true;
        }
      }
      else{//-Busca el plim mas cercano.
        double dsel=DBL_MAX;
        unsigned csel=0;
        for(unsigned c=0;c<nsel;c++){
          const double d=fgeo::PointDist(plim[c]-ps);
          if(dsel>d){
            dsel=d;
            csel=c;
          }
        }
        nor=plim[csel]-ps;
        nordone=true;
      }
    }

    //-Saves final normal.
    PartNor[p]=(nordone? nor: TDouble3(0));
    if(DG)if(nordone && fgeo::PointDist(nor)>Dist)printf("==> NorDist[%d]:%g    Dist:%f\n",p,fgeo::PointDist(nor),Dist);
  }
  //-Libera memoria dinamica.
  delete[] vnor; vnor=NULL;
  delete[] dnor; dnor=NULL;
  delete[] vout; vout=NULL;
  delete[] dout; dout=NULL;
  delete[] plim; plim=NULL;
  delete[] cmin; cmin=NULL;
}

//==============================================================================
/// Saves particles with normals (for debug).
//==============================================================================
void JNormalsMarrone::SaveVtkNormalFinal(std::string file){
  //-Allocates memory.
  unsigned *vidp=new unsigned[SizePart];
  tfloat3  *vpos=new tfloat3 [SizePart];
  tfloat3  *vnor=new tfloat3 [SizePart];
  float    *vtam=new float   [SizePart];
  //-Compute data.
  for(unsigned p=0;p<SizePart;p++){
    vidp[p]=p;
    vpos[p]=ToTFloat3(PartPos[p]);
    vnor[p]=ToTFloat3(PartNor[p]);
    vtam[p]=fgeo::PointDist(vnor[p]);
  }
  //-Saves VTK file and free memory.
  JDataArrays arrays;
  arrays.AddArray("Pos",SizePart,vpos,true);
  if(vidp)arrays.AddArray("Idp"       ,SizePart,vidp,true);
  if(vnor)arrays.AddArray("Normal"    ,SizePart,vnor,true);
  if(vtam)arrays.AddArray("NormalSize",SizePart,vtam,true);
  JVtkLib::SaveVtkData(file,arrays,"Pos");
  arrays.Reset();
}
