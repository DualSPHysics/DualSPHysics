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

/// \file JPartsLoad4.cpp \brief Implements the class \ref JPartsLoad4.

#include "JPartsLoad4.h"
#include "Functions.h"
#include "JPartDataBi4.h"
#include "JRadixSort.h"
#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JPartsLoad4
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JPartsLoad4::JPartsLoad4(bool useomp):UseOmp(useomp){
  ClassName="JPartsLoad4";
  Idp=NULL; Pos=NULL; VelRhop=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JPartsLoad4::~JPartsLoad4(){
  DestructorActive=true;
  AllocMemory(0);
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JPartsLoad4::Reset(){
  Simulate2D=false;
  Simulate2DPosY=0;
  CaseNp=CaseNfixed=CaseNmoving=CaseNfloat=CaseNfluid=0;
  PeriMode=PERI_Unknown;
  PeriXinc=PeriYinc=PeriZinc=TDouble3(0);
  MapSize=false;
  MapPosMin=MapPosMax=TDouble3(0);
  PartBegin=0;
  PartBeginTimeStep=0;
  SymplecticDtPre=0;
  DemDtForce=0;
  AllocMemory(0);
}

//==============================================================================
/// Resizes space for particle data.
/// Redimensiona espacio para datos de las particulas.
//==============================================================================
void JPartsLoad4::AllocMemory(unsigned count){
  Count=count;
  delete[] Idp;      Idp=NULL; 
  delete[] Pos;      Pos=NULL; 
  delete[] VelRhop;  VelRhop=NULL; 
  if(Count){
    try{
      Idp=new unsigned[Count];
      Pos=new tdouble3[Count];
      VelRhop=new tfloat4[Count];
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Could not allocate the requested memory.");
    }
  } 
}

//==============================================================================
/// Returns the reserved memory in CPU.
/// Devuelve la memoria reservada en CPU.
//==============================================================================
llong JPartsLoad4::GetAllocMemory()const{  
  llong s=0;
  //-Allocated in AllocMemory().
  if(Idp)s+=sizeof(unsigned)*Count;
  if(Pos)s+=sizeof(tdouble3)*Count;
  if(VelRhop)s+=sizeof(tfloat4)*Count;
  return(s);
}

//==============================================================================
/// Sorts values according to vsort[].
/// Ordena valores segun vsort[].
//==============================================================================
template<typename T> T* JPartsLoad4::SortParticles(const unsigned *vsort,unsigned count,T *v)const{
  T* v2=new T[count];
  for(unsigned p=0;p<count;p++)v2[p]=v[vsort[p]];
  delete[] v;
  return(v2);
}

//==============================================================================
/// Checks order of boundary particles.
/// Comprubeba orden de particulas de contorno.
//==============================================================================
void JPartsLoad4::CheckSortParticles(){
  const unsigned nbound=unsigned(CaseNfixed+CaseNmoving);
  if(nbound){
    unsigned lastbound=0;
    //-Computes position of last boundary particle.
    for(unsigned p=0;p<Count;p++)if(Idp[p]<nbound && p>lastbound)lastbound=p;
    if(lastbound+1!=nbound)Run_Exceptioon("Order of boundary (fixed and moving) particles is invalid.");
  }
}

//==============================================================================
/// Sorts values according to Idp[].
/// Ordena particulas por Idp[].
//==============================================================================
void JPartsLoad4::SortParticles(){
  //-Checks order. | Comprueba orden.
  bool sorted=true;
  for(unsigned p=1;p<Count && sorted;p++)sorted=(Idp[p-1]<Idp[p]);
  if(!sorted){
    //-Sorts points according to id. | Ordena puntos segun id.
    JRadixSort rs(UseOmp);
    rs.Sort(true,Count,Idp);
    rs.SortData(Count,Pos,Pos);
    rs.SortData(Count,VelRhop,VelRhop);
  }
}

//==============================================================================
/// It loads particles of bi4 file and it orders them by Id.
/// Carga particulas de fichero bi4 y las ordena por Id.
//==============================================================================
void JPartsLoad4::LoadParticles(const std::string &casedir,const std::string &casename
  ,unsigned partbegin,const std::string &casedirbegin)
{
  Reset();
  PartBegin=partbegin;
  JPartDataBi4 pd;
  //-Loads file piece_0 and obtains configuration.
  //-Carga fichero piece_0 y obtiene configuracion.
  const string dir=fun::GetDirWithSlash(!PartBegin? casedir: casedirbegin);
  if(!PartBegin){
    const string file1=dir+JPartDataBi4::GetFileNameCase(casename,0,1);
    if(fun::FileExists(file1))pd.LoadFileCase(dir,casename,0,1);
    else if(fun::FileExists(dir+JPartDataBi4::GetFileNameCase(casename,0,2)))pd.LoadFileCase(dir,casename,0,2);
    else Run_ExceptioonFile("File of the particles was not found.",file1);
  }
  else{
    const string file1=dir+JPartDataBi4::GetFileNamePart(PartBegin,0,1);
    if(fun::FileExists(file1))pd.LoadFilePart(dir,PartBegin,0,1);
    else if(fun::FileExists(dir+JPartDataBi4::GetFileNamePart(PartBegin,0,2)))pd.LoadFilePart(dir,PartBegin,0,2);
    else Run_ExceptioonFile("File of the particles was not found.",file1);
  }
  //-Obtains configuration. | Obtiene configuracion.
  PartBeginTimeStep=(!PartBegin? 0: pd.Get_TimeStep());
  Npiece=pd.GetNpiece();
  Simulate2D=pd.Get_Data2d();
  Simulate2DPosY=(Simulate2D? pd.Get_Data2dPosY(): 0);
  NpDynamic=pd.Get_NpDynamic();
  PartBeginTotalNp=(NpDynamic? pd.Get_NpTotal(): 0);
  CaseNp=pd.Get_CaseNp();
  CaseNfixed=pd.Get_CaseNfixed();
  CaseNmoving=pd.Get_CaseNmoving();
  CaseNfloat=pd.Get_CaseNfloat();
  CaseNfluid=pd.Get_CaseNfluid();
  PeriMode=(!PartBegin? PERI_Unknown: pd.Get_PeriMode());
  PeriXinc=pd.Get_PeriXinc();
  PeriYinc=pd.Get_PeriYinc();
  PeriZinc=pd.Get_PeriZinc();
  MapPosMin=pd.Get_MapPosMin();
  MapPosMax=pd.Get_MapPosMax();
  MapSize=(MapPosMin!=MapPosMax);
  CasePosMin=pd.Get_CasePosMin();
  CasePosMax=pd.Get_CasePosMax();
  const bool possingle=pd.Get_PosSimple();
  if(!pd.Get_IdpSimple())Run_Exceptioon("Only Idp (32 bits) is valid at the moment.");
  //-Loads data for restarting.
  if(PartBegin){
    SymplecticDtPre=pd.GetPart()->GetvDouble("SymplecticDtPre",true,0);
    DemDtForce=pd.GetPart()->GetvDouble("DemDtForce",true,0);
  }
  //-Calculates number of particles. | Calcula numero de particulas.
  unsigned sizetot=pd.Get_Npok();
  for(unsigned piece=1;piece<Npiece;piece++){
    JPartDataBi4 pd2;
    if(!PartBegin)pd2.LoadFileCase(dir,casename,piece,Npiece);
    else pd2.LoadFilePart(dir,PartBegin,piece,Npiece);
    sizetot+=pd.Get_Npok();
  }
  //-Allocates memory.
  AllocMemory(sizetot);
  //-Loads particles.
  {
    unsigned ntot=0;
    unsigned auxsize=0;
    tfloat3 *auxf3=NULL;
    float *auxf=NULL;
    for(unsigned piece=0;piece<Npiece;piece++){
      if(piece){
        if(!PartBegin)pd.LoadFileCase(dir,casename,piece,Npiece);
        else pd.LoadFilePart(dir,PartBegin,piece,Npiece);
      }
      const unsigned npok=pd.Get_Npok();
      if(npok){
        if(auxsize<npok){
          auxsize=npok;
          delete[] auxf3; auxf3=NULL;
          delete[] auxf;  auxf=NULL;
          auxf3=new tfloat3[auxsize];
          auxf=new float[auxsize];
        }
        if(possingle){
          pd.Get_Pos(npok,auxf3);
          for(unsigned p=0;p<npok;p++)Pos[ntot+p]=ToTDouble3(auxf3[p]);
        }
        else pd.Get_Posd(npok,Pos+ntot);
        pd.Get_Idp(npok,Idp+ntot);  
        pd.Get_Vel(npok,auxf3);  
        pd.Get_Rhop(npok,auxf);  
        for(unsigned p=0;p<npok;p++)VelRhop[ntot+p]=TFloat4(auxf3[p].x,auxf3[p].y,auxf3[p].z,auxf[p]);
      }
      ntot+=npok;
    }
    delete[] auxf3; auxf3=NULL;
    delete[] auxf;  auxf=NULL;
  }
  //-In simulations 2D, if PosY is invalid then calculates starting from position of particles.
  if(Simulate2DPosY==DBL_MAX){
    if(!sizetot)Run_Exceptioon("Number of particles is invalid to calculates Y in 2D simulations.");
    Simulate2DPosY=Pos[0].y;
  }
  //-Checks order of boundary particles.
  CheckSortParticles();
  //-Sorts particles according to Id. | Ordena particulas por Id.
  //SortParticles();
}

//==============================================================================
/// Check validity of loaded configuration or throw exception.
/// Comprueba validez de la configuracion cargada o lanza excepcion.
//==============================================================================
void JPartsLoad4::CheckConfig(ullong casenp,ullong casenfixed,ullong casenmoving
  ,ullong casenfloat,ullong casenfluid,bool simulate2d,double simulate2dposy,TpPeri tperi)const
{
  CheckConfig(casenp,casenfixed,casenmoving,casenfloat,casenfluid);
  if(simulate2d!=Simulate2D || Simulate2DPosY!=simulate2dposy)Run_Exceptioon("Data file does not match the dimension of the case (2D/3D).");
  //-Obtains periodic mode and compares with loaded file.
  if(tperi!=PeriMode && PeriMode!=PERI_Unknown)Run_Exceptioon("Data file does not match the periodic configuration of the case.");
}

//==============================================================================
/// Check validity of loaded configuration or throw exception.
/// Comprueba validez de la configuracion cargada o lanza excepcion.
//==============================================================================
void JPartsLoad4::CheckConfig(ullong casenp,ullong casenfixed,ullong casenmoving
  ,ullong casenfloat,ullong casenfluid)const
{
  if(casenp!=CaseNp || casenfixed!=CaseNfixed || casenmoving!=CaseNmoving || casenfloat!=CaseNfloat || casenfluid!=CaseNfluid)
    Run_Exceptioon("Particle number does not match the configuration of the case.");
}

//==============================================================================
/// Removes boundary particles.
/// Elimina particulas de contorno.
//==============================================================================
void JPartsLoad4::RemoveBoundary(){
  const unsigned casenbound=unsigned(CaseNp-CaseNfluid);
  //-Calculate number of boundary particles. 
  unsigned nbound=0;
  for(;nbound<Count && Idp[nbound]<casenbound;nbound++);
  //-Saves old pointers and allocates new memory.
  unsigned count0=Count;
  unsigned *idp0=Idp;        Idp=NULL;
  tdouble3 *pos0=Pos;        Pos=NULL;
  tfloat4 *velrhop0=VelRhop; VelRhop=NULL;
  AllocMemory(count0-nbound);
  //-Copies data in new pointers.
  memcpy(Idp,idp0+nbound,sizeof(unsigned)*Count);
  memcpy(Pos,pos0+nbound,sizeof(tdouble3)*Count);
  memcpy(VelRhop,velrhop0+nbound,sizeof(tfloat4)*Count);
  //-Frees old pointers.
  delete[] idp0;      idp0=NULL; 
  delete[] pos0;      pos0=NULL; 
  delete[] velrhop0;  velrhop0=NULL; 
}

//==============================================================================
/// Returns the limits of the map and if they are not valid throws exception.
/// Devuelve los limites del mapa y si no son validos genera excepcion.
//==============================================================================
void JPartsLoad4::GetMapSize(tdouble3 &mapmin,tdouble3 &mapmax)const{
  if(!MapSizeLoaded())Run_Exceptioon("The MapSize information is invalid.");
  mapmin=MapPosMin; mapmax=MapPosMax;
}

//==============================================================================
/// Calculates limits of the loaded particles.
/// Calcula limites de las particulas cargadas.
//==============================================================================
void JPartsLoad4::CalculateCasePos(){
  if(!PartBegin)Run_Exceptioon("The limits of the initial case cannot be calculated from a file PART.");
  tdouble3 pmin=TDouble3(DBL_MAX),pmax=TDouble3(-DBL_MAX);
  //-Calculates minimum and maximum position. 
  //-Calcula posicion minima y maxima. 
  for(unsigned p=0;p<Count;p++){
    const tdouble3 ps=Pos[p];
    if(pmin.x>ps.x)pmin.x=ps.x;
    if(pmin.y>ps.y)pmin.y=ps.y;
    if(pmin.z>ps.z)pmin.z=ps.z;
    if(pmax.x<ps.x)pmax.x=ps.x;
    if(pmax.y<ps.y)pmax.y=ps.y;
    if(pmax.z<ps.z)pmax.z=ps.z;
  }
  CasePosMin=pmin; CasePosMax=pmax;
}

//==============================================================================
/// Calculates and returns particle limits with the indicated border.
/// Calcula y devuelve limites de particulas con el borde indicado.
//==============================================================================
void JPartsLoad4::CalculeLimits(double border,double borderperi,bool perix,bool periy,bool periz,tdouble3 &mapmin,tdouble3 &mapmax){
  if(CasePosMin==CasePosMax)CalculateCasePos();
  tdouble3 bor=TDouble3(border);
  if(perix)bor.x=borderperi;
  if(periy)bor.y=borderperi;
  if(periz)bor.z=borderperi;
  mapmin=CasePosMin-bor;
  mapmax=CasePosMax+bor;
}


