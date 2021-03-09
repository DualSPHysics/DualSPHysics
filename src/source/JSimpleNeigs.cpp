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

/// \file JSimpleNeigs.cpp \brief Implements the class \ref JSimpleNeigs.

#include "JSimpleNeigs.h"
#include "Functions.h"
#include <cstring>
#include <cmath>
#include <cfloat>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSimpleNeigs
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSimpleNeigs::JSimpleNeigs(unsigned np,const tdouble3* pos,double scell):Np(np),Pos(pos),Scell(scell){
  ClassName="JSimpleNeigs";
  PosInCell=NULL;
  BeginCell=NULL;
  SelectPos=NULL;
  Reset();
  CreateMapCells();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSimpleNeigs::~JSimpleNeigs(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JSimpleNeigs::Reset(){
  PosMin=PosMax=TDouble3(0);
  Ncx=Ncy=Ncz=Nsheet=Nct=0;
  CountSelect=SizeSelect=0;
  delete[] PosInCell; PosInCell=NULL;
  delete[] BeginCell; BeginCell=NULL;
  delete[] SelectPos; SelectPos=NULL;
}

//==============================================================================
/// Returns the allocated memory.
/// Devuelve la memoria reservada.
//==============================================================================
unsigned JSimpleNeigs::GetAllocMemory()const{
  unsigned s=0;
  if(PosInCell)s+=sizeof(unsigned)*Np;
  if(BeginCell)s+=sizeof(unsigned)*(Nct+1);
  if(SelectPos)s+=sizeof(unsigned)*SizeSelect;
  return(s);
}

//==============================================================================
/// Defines division of domain in cells.
/// Definie la division del dominio en celdas.
//==============================================================================
void JSimpleNeigs::DefineMapCells(){
  //-Calculates minimum and maximum position. 
  tdouble3 pmin=TDouble3(DBL_MAX),pmax=TDouble3(-DBL_MAX);
  for(unsigned p=0;p<Np;p++){
    const tdouble3 ps=Pos[p];
    if(pmin.x>ps.x)pmin.x=ps.x;
    if(pmin.y>ps.y)pmin.y=ps.y;
    if(pmin.z>ps.z)pmin.z=ps.z;
    if(pmax.x<ps.x)pmax.x=ps.x;
    if(pmax.y<ps.y)pmax.y=ps.y;
    if(pmax.z<ps.z)pmax.z=ps.z;
  }
  //-Adds border to minimum and maximum.
  const double border=Scell*0.1;
  PosMin=pmin-TDouble3(border); 
  PosMax=pmax+TDouble3(border);
  //-Calculates number o cells.
  tdouble3 mapsize=(PosMax-PosMin)/TDouble3(Scell);
  Ncx=int(ceil(mapsize.x));
  Ncy=int(ceil(mapsize.y));
  Ncz=int(ceil(mapsize.z));
  //printf("==>  ncx:%u ncy:%u ncz:%u\n",Ncx,Ncy,Ncz);
  Nsheet=Ncx*Ncy; Nct=Nsheet*Ncz;
  ullong nct0=ullong(Ncx)*ullong(Ncy)*ullong(Ncz);
  const ullong ncmax=ullong(1024*1024*512);
  if(nct0>ncmax)Run_Exceptioon(fun::PrintStr("Number of cells (%d x %d x %d) is too high.",Ncx,Ncy,Ncz));
  if(ullong(Nct)!=nct0)Run_Exceptioon(fun::PrintStr("Number of cells (%d x %d x %d) is invalid.",Ncx,Ncy,Ncz));
}

//==============================================================================
/// Creates map of cells.
/// Crea mapa de celdas.
//==============================================================================
void JSimpleNeigs::CreateMapCells(){
  if(Np==0)Run_Exceptioon("Nummber of positions is zero.");
  if(Scell<=0)Run_Exceptioon("Size of cells is invalid.");
  DefineMapCells();
  //-Allocate memory.
  unsigned *poscell=NULL;
  unsigned *npcell=NULL;
  try{
    poscell=new unsigned[Np];
    npcell=new unsigned[Nct];
    PosInCell=new unsigned[Np];
    BeginCell=new unsigned[Nct+1];
    SizeSelect=100;
    SelectPos=new unsigned[SizeSelect];
  }
  catch(const std::bad_alloc){
    Run_Exceptioon("Could not allocate the requested memory.");
  }
  //-Computes cell of each particle and number of postions for each cell.
  memset(npcell,0,sizeof(unsigned)*Nct);
  bool error=false;
  for(unsigned p=0;p<Np;p++){
    const unsigned cel=GetCell(GetCell3(Pos[p]));
    if(int(cel)<=Nct){
      poscell[p]=cel;
      npcell[cel]++;
    }
    else error=true;
  }
  if(error)Run_Exceptioon("Some position is outside the defined domain.");
  //-Computes BeginCell[].
  BeginCell[0]=0;
  for(int c=0;c<Nct;c++)BeginCell[c+1]=BeginCell[c]+npcell[c];
  //-Computes PosInCell[].
  memset(npcell,0,sizeof(unsigned)*Nct);
  for(unsigned p=0;p<Np;p++){
    const unsigned cel=poscell[p];
    PosInCell[BeginCell[cel]+npcell[cel]]=p;
    npcell[cel]++;
  }
  //-Free auxiliary memory.
  delete[] poscell; poscell=NULL;
  delete[] npcell;  npcell=NULL;
}

//==============================================================================
/// Return cell limits for interaction starting from position.
/// Devuelve limites de celdas para interaccion a partir de posicion.
//==============================================================================
void JSimpleNeigs::GetNearbyCells(const tdouble3 &ps,double dist,tint3 &celmin,tint3 &celmax)const{
  //-Compute distance in cells.
  const int celdist=int(dist/Scell)+1;
  //-Get cell coordinates of position.
  const tint3 cel=GetCell3(ps);
  //-Compute cell domain to look for.
  celmin=MaxValues(TInt3(0),cel-TInt3(celdist));
  celmax=MinValues(TInt3(Ncx-1,Ncy-1,Ncz-1),cel+TInt3(celdist));
}

//==============================================================================
/// Add position to SelectPos[].
/// Incluye posicion en SelectPos[].
//==============================================================================
void JSimpleNeigs::SelectAdd(unsigned p){
  if(CountSelect>=SizeSelect){
    SizeSelect+=200;
    SelectPos=fun::ResizeAlloc(SelectPos,CountSelect,SizeSelect);
  }
  SelectPos[CountSelect]=p;
  CountSelect++;
}

//==============================================================================
/// Store nearby positions in SelectPos[] and returns number of selected positions.
/// Guarda las posiciones cercanas en SelectPos[] y devuelve el numero de posiciones seleccionadas.
//==============================================================================
unsigned JSimpleNeigs::NearbyPositions(const tdouble3 &ps,unsigned pignore,double dist){
  const double dist2=dist*dist;
  CountSelect=0;
  //printf("==> pos:(%f,%f,%f)\n",ps.x,ps.y,ps.z);
  tint3 celmin,celmax;
  GetNearbyCells(ps,dist,celmin,celmax);
  //printf("==> NearbyCells:%s\n",fun::Int3RangeStr(celmin,celmax).c_str());
  for(int cz=celmin.z;cz<=celmax.z;cz++)for(int cy=celmin.y;cy<=celmax.y;cy++){
    const unsigned cmin=GetCell(TInt3(celmin.x,cy,cz));
    const unsigned cmax=GetCell(TInt3(celmax.x,cy,cz));
    const unsigned pini=BeginCell[cmin];
    const unsigned pfin=BeginCell[cmax+1];
    for(unsigned cp=pini;cp<pfin;cp++){
      const unsigned p=PosInCell[cp];
      const tdouble3 ds=ps-Pos[p];
      if(ds.x*ds.x+ds.y*ds.y+ds.z*ds.z<=dist2 && p!=pignore)SelectAdd(p);
    }
  }
  return(CountSelect);
}

//==============================================================================
/// Store nearby positions in vector vsel and returns number of selected positions.
/// Guarda las posiciones cercanas en vsel y devuelve el numero de posiciones seleccionadas.
//==============================================================================
unsigned JSimpleNeigs::NearbyPositionsLt(const tdouble3 &ps,unsigned pignore
  ,double dist,std::vector<unsigned> &vsel)const
{
  vsel.clear();
  const double dist2=dist*dist;
  //printf("==> pos:(%f,%f,%f)\n",ps.x,ps.y,ps.z);
  tint3 celmin,celmax;
  GetNearbyCells(ps,dist,celmin,celmax);
  //printf("==> NearbyCells:%s\n",fun::Int3RangeStr(celmin,celmax).c_str());
  for(int cz=celmin.z;cz<=celmax.z;cz++)for(int cy=celmin.y;cy<=celmax.y;cy++){
    const unsigned cmin=GetCell(TInt3(celmin.x,cy,cz));
    const unsigned cmax=GetCell(TInt3(celmax.x,cy,cz));
    const unsigned pini=BeginCell[cmin];
    const unsigned pfin=BeginCell[cmax+1];
    for(unsigned cp=pini;cp<pfin;cp++){
      const unsigned p=PosInCell[cp];
      const tdouble3 ds=ps-Pos[p];
      if(ds.x*ds.x+ds.y*ds.y+ds.z*ds.z<dist2 && p!=pignore)vsel.push_back(p);
    }
  }
  return(unsigned(vsel.size()));
}

