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

/// \file JSpVtkShape.cpp \brief Implements the class \ref JSpVtkShape.

#include "JSpVtkShape.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JMatrix4.h"

#include <fstream>
#include <cstring>

using namespace std;

//##############################################################################
//# JSpVtkShape
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSpVtkShape::JSpVtkShape(){
  ClassName="JSpVtkShape";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSpVtkShape::~JSpVtkShape(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSpVtkShape::Reset(){
  NpTp[0]=NpTp[1]=NpTp[2]=0;
  NcTp[0]=NcTp[1]=NcTp[2]=0;
  Ptos.clear();
  Cells.clear();
}

//==============================================================================
/// Saves Ptos data.
//==============================================================================
void JSpVtkShape::SavePtos(std::ofstream& pf)const{
  const unsigned np=PtCount();
  const unsigned nc=CellCount();
  const tfloat3* vpto=(const tfloat3*)Ptos.data();
  tfloat3* vp=new tfloat3[np];
  unsigned cpw[3];
  cpw[0]=0;
  cpw[1]=NpTp[0];
  cpw[2]=NpTp[0]+NpTp[1];
  unsigned cpr=0;
  for(unsigned cc=0;cc<nc;cc++){
     const StCell cell=Cells[cc];
     memcpy(vp+cpw[cell.tp],vpto+cpr,sizeof(tfloat3)*cell.np);
     cpr+=cell.np;
     cpw[cell.tp]+=cell.np;
  }
  pf << "POINTS " << np << " float" << endl;
  if(fun::GetByteOrder()==fun::LittleEndian)
    fun::SetByteOrder32((int*)vp,np*3);
  pf.write((const char*)vp,sizeof(tfloat3)*np);
  pf << endl;
  delete[] vp;
}

//==============================================================================
/// Saves points data.
//==============================================================================
void JSpVtkShape::SaveVerts(std::ofstream& pf)const{
  if(NcTp[0]){
    pf << "VERTICES " << NcTp[0] << " " << NcTp[0]*2 << endl;
    const unsigned size=NcTp[0]*2;
    unsigned* data=new unsigned[size];
    const unsigned nc=CellCount();
    unsigned cp=0,n=0;
    for(unsigned cc=0;cc<nc;cc++)if(Cells[cc].tp==0){
      data[cp++]=1;
      data[cp++]=n++;
    }
    if(fun::GetByteOrder()==fun::LittleEndian)
      fun::SetByteOrder32((int*)data,size);
    pf.write((char*)data,sizeof(int)*size);
    pf << endl;
    delete[] data;
  }
}

//==============================================================================
/// Saves lines data.
//==============================================================================
void JSpVtkShape::SaveLines(std::ofstream& pf)const{
  if(NcTp[1]){
    pf << "LINES " << NcTp[1] << " " << NcTp[1]*3 << endl;
    const unsigned size=NcTp[1]*3;
    unsigned* data=new unsigned[size];
    const unsigned nc=CellCount();
    unsigned cp=0,n=0,n0=NpTp[0];
    for(unsigned cc=0;cc<nc;cc++)if(Cells[cc].tp==1){
      data[cp++]=2;        
      data[cp++]=n0+n*2;   
      data[cp++]=n0+n*2+1; 
      n++;
    }
    if(fun::GetByteOrder()==fun::LittleEndian)
      fun::SetByteOrder32((int*)data,size);
    pf.write((char*)data,sizeof(int)*size);
    pf << endl;
    delete[] data;
  }
}

//==============================================================================
/// Saves polygon data.
//==============================================================================
void JSpVtkShape::SavePoly(std::ofstream& pf)const{
  if(NcTp[2]){
    const unsigned nc=CellCount();
    unsigned size=NcTp[2];
    for(unsigned cc=0;cc<nc;cc++)if(Cells[cc].tp==2)size+=Cells[cc].np;
    pf << "POLYGONS " << NcTp[2] << " " << size << endl;
    unsigned* data=new unsigned[size];
    unsigned cp=0,n0=NpTp[0]+NpTp[1];
    for(unsigned cc=0;cc<nc;cc++)if(Cells[cc].tp==2){
      const StCell cell=Cells[cc];
      data[cp++]=cell.np;
      for(unsigned cv=0;cv<cell.np;cv++)data[cp++]=n0++;
    }
    if(fun::GetByteOrder()==fun::LittleEndian)
      fun::SetByteOrder32((int*)data,size);
    pf.write((char*)data,sizeof(int)*size);
    pf << endl;
    delete[] data;
  }
}

//==============================================================================
/// Saves lines data.
//==============================================================================
void JSpVtkShape::SaveData(const std::string& vname,std::ofstream& pf)const{
  if(!vname.empty()){
    const unsigned nc=CellCount();
    pf << "CELL_DATA " << nc << endl;
    pf << "SCALARS " << vname << " int" << endl;
    pf << "LOOKUP_TABLE default" << endl;
    int* data=new int[nc];
    unsigned cpw[3];
    cpw[0]=0;
    cpw[1]=NcTp[0];
    cpw[2]=NcTp[0]+NcTp[1];
    unsigned cpr=0;
    for(unsigned cc=0;cc<nc;cc++){
       const StCell cell=Cells[cc];
       data[cpw[cell.tp]++]=cell.v;
    }
    if(fun::GetByteOrder()==fun::LittleEndian)
      fun::SetByteOrder32((int*)data,nc);
    pf.write((char*)data,sizeof(int)*nc);
    delete[] data;
    pf << endl;
  }
}

//==============================================================================
/// Saves VTK file.
//==============================================================================
void JSpVtkShape::SaveVtk(std::string file,std::string vname)const{
  if(fun::GetByteOrder()!=fun::LittleEndian)
    Run_ExceptioonFile("Big-Endian mode is not supported for now.",file);
  const unsigned np=PtCount();
  const bool empty=(np==0);
  fun::MkdirPath(fun::GetDirParent(file));
  std::ofstream pf;
  pf.open(file.c_str(),ios::binary|ios::out);
  if(pf){
    pf << "# vtk DataFile Version 3.0" << endl;
    pf << "vtk output" << endl;
    pf << "BINARY" << endl;
    pf << "DATASET POLYDATA" << endl;
    if(empty)pf << "POINTS 0 float" << endl;
    else{
      SavePtos(pf);
      SaveVerts(pf);
      SaveLines(pf);
      SavePoly(pf);
      SaveData(vname,pf);
    }
    if(pf.fail())Run_ExceptioonFile("File writing failure.",file);
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",file);
}

//==============================================================================
/// Add point.
//==============================================================================
void JSpVtkShape::AddPoint(const tfloat3& p1,word v){
  Ptos.push_back(p1);
  const StCell cell={0,1,v};
  Cells.push_back(cell);
  NpTp[cell.tp]+=cell.np;
  NcTp[cell.tp]++;
}
//==============================================================================
/// Add point.
//==============================================================================
void JSpVtkShape::AddPoint(const tdouble3& p1,word v){
  AddPoint(ToTFloat3(p1),v);
}

//==============================================================================
/// Add points.
//==============================================================================
void JSpVtkShape::AddPoints(unsigned npoints,const tfloat3* points,word v){
  for(unsigned c=0;c<npoints;c++)AddPoint(points[c],v);
}
//==============================================================================
/// Add lines.
//==============================================================================
void JSpVtkShape::AddPoints(unsigned npoints,const tdouble3* points,word v){
  for(unsigned c=0;c<npoints;c++)AddPoint(points[c],v);
}

//==============================================================================
/// Add line.
//==============================================================================
void JSpVtkShape::AddLine(const tfloat3& p1,const tfloat3& p2,word v){
  Ptos.push_back(p1);
  Ptos.push_back(p2);
  const StCell cell={1,2,v};
  Cells.push_back(cell);
  NpTp[cell.tp]+=cell.np;
  NcTp[cell.tp]++;
}
//==============================================================================
/// Add line.
//==============================================================================
void JSpVtkShape::AddLine(const tdouble3& p1,const tdouble3& p2,word v){
  AddLine(ToTFloat3(p1),ToTFloat3(p2),v);
}

//==============================================================================
/// Add lines.
//==============================================================================
void JSpVtkShape::AddLines(unsigned npoints,const tfloat3* points,word v){
  for(unsigned c=0;c<npoints-1;c++)AddLine(points[c],points[c+1],v);
}
//==============================================================================
/// Add lines.
//==============================================================================
void JSpVtkShape::AddLines(unsigned npoints,const tdouble3* points,word v){
  for(unsigned c=0;c<npoints-1;c++)AddLine(points[c],points[c+1],v);
}

//==============================================================================
/// Add triangle.
//==============================================================================
void JSpVtkShape::AddTriangle(const tfloat3& p1,const tfloat3& p2
  ,const tfloat3& p3,word v)
{
  Ptos.push_back(p1);
  Ptos.push_back(p2);
  Ptos.push_back(p3);
  const StCell cell={2,3,v};
  Cells.push_back(cell);
  NpTp[cell.tp]+=cell.np;
  NcTp[cell.tp]++;
}
//==============================================================================
/// Add triangle.
//==============================================================================
void JSpVtkShape::AddTriangle(const tdouble3& p1,const tdouble3& p2
  ,const tdouble3& p3,word v)
{
  AddTriangle(ToTFloat3(p1),ToTFloat3(p2),ToTFloat3(p3),v);
}

//==============================================================================
/// Add triangle (wire mode).
//==============================================================================
void JSpVtkShape::AddTriangleWire(const tfloat3& p1,const tfloat3& p2
  ,const tfloat3& p3,word v)
{
  AddLine(p1,p2,v);
  AddLine(p2,p3,v);
  AddLine(p3,p1,v);
}
//==============================================================================
/// Add triangle (wire mode).
//==============================================================================
void JSpVtkShape::AddTriangleWire(const tdouble3& p1,const tdouble3& p2
  ,const tdouble3& p3,word v)
{
  AddTriangle(ToTFloat3(p1),ToTFloat3(p2),ToTFloat3(p3),v);
}

//==============================================================================
/// Add quad.
//==============================================================================
void JSpVtkShape::AddQuad(const tfloat3& p1,const tfloat3& p2
  ,const tfloat3& p3,const tfloat3& p4,word v)
{
  Ptos.push_back(p1);
  Ptos.push_back(p2);
  Ptos.push_back(p3);
  Ptos.push_back(p4);
  const StCell cell={2,4,v};
  Cells.push_back(cell);
  NpTp[cell.tp]+=cell.np;
  NcTp[cell.tp]++;
}
//==============================================================================
/// Add quad.
//==============================================================================
void JSpVtkShape::AddQuad(const tdouble3& p1,const tdouble3& p2
  ,const tdouble3& p3,const tdouble3& p4,word v)
{
  AddQuad(ToTFloat3(p1),ToTFloat3(p2),ToTFloat3(p3),ToTFloat3(p4),v);
}

//==============================================================================
/// Add quad.
//==============================================================================
void JSpVtkShape::AddQuadWire(const tfloat3& p1,const tfloat3& p2
  ,const tfloat3& p3,const tfloat3& p4,word v)
{
  AddLine(p1,p2,v);
  AddLine(p2,p3,v);
  AddLine(p3,p4,v);
  AddLine(p4,p1,v);
}
//==============================================================================
/// Add quad.
//==============================================================================
void JSpVtkShape::AddQuadWire(const tdouble3& p1,const tdouble3& p2
  ,const tdouble3& p3,const tdouble3& p4,word v)
{
  AddQuadWire(ToTFloat3(p1),ToTFloat3(p2),ToTFloat3(p3),ToTFloat3(p4),v);
}

//==============================================================================
/// Add quad with orthogonal vector.
//==============================================================================
void JSpVtkShape::AddQuadOrtho(const tfloat3& p1,const tfloat3& v1
  ,float sidesize,word v)
{
  const tfloat3 vdir=fgeo::VecUnitary(v1);
  const tfloat3 vort1=fgeo::VecOrthogonal2(vdir,sidesize,true);
  const tfloat3 vort2=fgeo::VecUnitary(fgeo::ProductVec(vdir,vort1))*sidesize;
  const tfloat3 pos0=p1-(vort1/2)-(vort2/2);
  AddQuad(pos0,pos0+vort1,pos0+vort1+vort2,pos0+vort2,v);
}
//==============================================================================
/// Add quad with orthogonal vector.
//==============================================================================
void JSpVtkShape::AddQuadOrtho(const tdouble3& p1,const tdouble3& v1
  ,double sidesize,word v)
{
  AddQuadOrtho(ToTFloat3(p1),ToTFloat3(v1),float(sidesize),v);
}

//==============================================================================
/// Add quad with orthogonal vector (wire mode).
//==============================================================================
void JSpVtkShape::AddQuadOrthoWire(const tfloat3& p1,const tfloat3& v1
  ,float sidesize,word v)
{
  const tfloat3 vdir=fgeo::VecUnitary(v1);
  const tfloat3 vort1=fgeo::VecOrthogonal2(vdir,sidesize,true);
  const tfloat3 vort2=fgeo::VecUnitary(fgeo::ProductVec(vdir,vort1))*sidesize;
  const tfloat3 pos0=p1-(vort1/2)-(vort2/2);
  AddQuadWire(pos0,pos0+vort1,pos0+vort1+vort2,pos0+vort2,v);
}
//==============================================================================
/// Add quad with orthogonal vector (wire mode).
//==============================================================================
void JSpVtkShape::AddQuadOrthoWire(const tdouble3& p1,const tdouble3& v1
  ,double sidesize,word v)
{
  AddQuadOrtho(ToTFloat3(p1),ToTFloat3(v1),float(sidesize),v);
}

//==============================================================================
/// Add polygon.
//==============================================================================
void JSpVtkShape::AddPolygon(unsigned np,const tfloat3* vp,word v){
  for(unsigned p=0;p<np;p++)Ptos.push_back(vp[p]);
  const StCell cell={2,word(np),v};
  Cells.push_back(cell);
  NpTp[cell.tp]+=cell.np;
  NcTp[cell.tp]++;
}

//==============================================================================
/// Add box with orthogonal size vectors.
//==============================================================================
void JSpVtkShape::AddBoxSizeVec(const tfloat3& p0,const tfloat3& sizex
  ,const tfloat3& sizey,const tfloat3& sizez,word v)
{
  tfloat3 pos[8];
  pos[0]=p0;
  pos[1]=p0+sizex;
  pos[2]=pos[1]+sizey;
  pos[3]=p0+sizey;
  pos[4]=p0+sizez;
  pos[5]=pos[1]+sizez;
  pos[6]=pos[2]+sizez;
  pos[7]=pos[3]+sizez;
  AddQuad(pos[0],pos[3],pos[2],pos[1],v);
  AddQuad(pos[4],pos[5],pos[6],pos[7],v);
  AddQuad(pos[0],pos[1],pos[5],pos[4],v);
  AddQuad(pos[2],pos[3],pos[7],pos[6],v);
  AddQuad(pos[1],pos[2],pos[6],pos[5],v);
  AddQuad(pos[0],pos[4],pos[7],pos[3],v);
}
//==============================================================================
/// Add box with size vectors.
//==============================================================================
void JSpVtkShape::AddBoxSizeVec(const tdouble3& p0,const tdouble3& sizex
  ,const tdouble3& sizey,const tdouble3& sizez,word v)
{
  AddBoxSizeVec(ToTFloat3(p0),ToTFloat3(sizex),ToTFloat3(sizey)
    ,ToTFloat3(sizez),v);
}

//==============================================================================
/// Add box with orthogonal size vectors (wire mode).
//==============================================================================
void JSpVtkShape::AddBoxSizeVecWire(const tfloat3& p0,const tfloat3& sizex
  ,const tfloat3& sizey,const tfloat3& sizez,word v)
{
  tfloat3 pos[8];
  pos[0]=p0;
  pos[1]=p0+sizex;
  pos[2]=pos[1]+sizey;
  pos[3]=p0+sizey;
  pos[4]=p0+sizez;
  pos[5]=pos[1]+sizez;
  pos[6]=pos[2]+sizez;
  pos[7]=pos[3]+sizez;
  AddQuadWire(pos[0],pos[3],pos[2],pos[1],v);
  AddQuadWire(pos[4],pos[5],pos[6],pos[7],v);
  AddQuadWire(pos[0],pos[1],pos[5],pos[4],v);
  AddQuadWire(pos[2],pos[3],pos[7],pos[6],v);
  AddQuadWire(pos[1],pos[2],pos[6],pos[5],v);
  AddQuadWire(pos[0],pos[4],pos[7],pos[3],v);
}
//==============================================================================
/// Add box with size vectors (wire mode).
//==============================================================================
void JSpVtkShape::AddBoxSizeVecWire(const tdouble3& p0,const tdouble3& sizex
  ,const tdouble3& sizey,const tdouble3& sizez,word v)
{
  AddBoxSizeVecWire(ToTFloat3(p0),ToTFloat3(sizex),ToTFloat3(sizey)
    ,ToTFloat3(sizez),v);
}

//==============================================================================
/// Add orthogonal box.
//==============================================================================
void JSpVtkShape::AddBoxSize(const tfloat3& p0,const tfloat3& sizexyz,word v){
  AddBoxSizeVec(p0,TFloat3(sizexyz.x,0,0),TFloat3(0,sizexyz.y,0),TFloat3(0,0,sizexyz.z),v);
}
//==============================================================================
/// Add orthogonal box.
//==============================================================================
void JSpVtkShape::AddBoxSize(const tdouble3& p0,const tdouble3& sizexyz,word v)
{
  AddBoxSize(ToTFloat3(p0),ToTFloat3(sizexyz),v);
}

//==============================================================================
/// Add orthogonal box (wire mode).
//==============================================================================
void JSpVtkShape::AddBoxSizeWire(const tfloat3& p0,const tfloat3& sizexyz,word v){
  AddBoxSizeVecWire(p0,TFloat3(sizexyz.x,0,0),TFloat3(0,sizexyz.y,0),TFloat3(0,0,sizexyz.z));
}
//==============================================================================
/// Add orthogonal box (wire mode).
//==============================================================================
void JSpVtkShape::AddBoxSizeWire(const tdouble3& p0,const tdouble3& sizexyz,word v)
{
  AddBoxSizeWire(ToTFloat3(p0),ToTFloat3(sizexyz),v);
}

//==============================================================================
/// Add box (4pt front + 4pt back).
//==============================================================================
void JSpVtkShape::AddBoxFront(const tfloat3& p0,const tfloat3& px
  ,const tfloat3& pxz,const tfloat3& pz,const tfloat3& py,const tfloat3& pyx
  ,const tfloat3& pyxz,const tfloat3& pyz,word v)
{
  AddBoxSizeVec(p0,px-p0,py-p0,pz-p0,v);
}
//==============================================================================
/// Add box (4pt front + 4pt back).
//==============================================================================
void JSpVtkShape::AddBoxFront(const tdouble3& p0,const tdouble3& px
  ,const tdouble3& pxz,const tdouble3& pz,const tdouble3& py,const tdouble3& pyx
  ,const tdouble3& pyxz,const tdouble3& pyz,word v)
{
  AddBoxFront(ToTFloat3(p0),ToTFloat3(px),ToTFloat3(pxz),ToTFloat3(pz)
    ,ToTFloat3(py),ToTFloat3(pyx),ToTFloat3(pyxz),ToTFloat3(pyz),v);
}

//==============================================================================
/// Add boxes.
//==============================================================================
void JSpVtkShape::AddBoxes(unsigned nbox,const tfloat3* vbox,float sizemin){
  for(unsigned d=0;d<nbox;d++){
    const tfloat3 pp=vbox[d*2];
    tfloat3 size=vbox[d*2+1]-pp;
    if(!size.x)size.x=sizemin;
    if(!size.y)size.y=sizemin;
    if(!size.z)size.z=sizemin;
    AddBoxSize(pp,size,word(d));
  }
}
//==============================================================================
/// Add boxes.
//==============================================================================
void JSpVtkShape::AddBoxes(unsigned nbox,const tdouble3* vbox,double sizemin){
  for(unsigned d=0;d<nbox;d++){
    const tfloat3 pp=ToTFloat3(vbox[d*2]);
    tfloat3 size=ToTFloat3(vbox[d*2+1])-pp;
    if(!size.x)size.x=float(sizemin);
    if(!size.y)size.y=float(sizemin);
    if(!size.z)size.z=float(sizemin);
    AddBoxSize(pp,size,word(d));
  }
}

//==============================================================================
/// Add cylinder.
//==============================================================================
void JSpVtkShape::AddCylinder(const tfloat3& cen1,const tfloat3& cen2
  ,float radius,word v)
{
  const int nside=28;
  tfloat3 vcen=TFloat3(cen2.x-cen1.x,cen2.y-cen1.y,cen2.z-cen1.z);
  tfloat3 vp=fgeo::VecOrthogonal(vcen,radius);
  tfloat3 p1b=TFloat3(cen1.x+vp.x,cen1.y+vp.y,cen1.z+vp.z);
  tfloat3 p2b=TFloat3(cen2.x+vp.x,cen2.y+vp.y,cen2.z+vp.z);
  tfloat3* pos=new tfloat3[56];
  pos[0]=p1b; pos[nside]=p2b;
  JMatrix4f mt=JMatrix4f::MatrixRot(360.f/28,cen1,cen2);
  for(int c=1;c<nside;c++){
    p1b=mt.MulPoint(p1b);
    p2b=mt.MulPoint(p2b);
    pos[c]=p1b; pos[nside+c]=p2b;
  }
  for(int c=0;c<nside;c++){
    AddTriangle(pos[c],pos[(c+1)%nside],cen1,v);
    AddTriangle(pos[nside+c],cen2,pos[nside+((c+1)%nside)],v);
    AddQuad(pos[(c+1)%nside],pos[c],pos[nside+c],pos[nside+((c+1)%nside)],v);
  }
  delete[] pos;
}
//==============================================================================
/// Add cylinder.
//==============================================================================
void JSpVtkShape::AddCylinder(const tdouble3& cen1,const tdouble3& cen2
  ,double radius,word v)
{
  AddCylinder(ToTFloat3(cen1),ToTFloat3(cen2),float(radius),v);
}

//==============================================================================
/// Add cylinder (wire mode).
//==============================================================================
void JSpVtkShape::AddCylinderWire(const tfloat3& cen1,const tfloat3& cen2
  ,float radius,word v)
{
  const int nside=28;
  tfloat3 vcen=TFloat3(cen2.x-cen1.x,cen2.y-cen1.y,cen2.z-cen1.z);
  tfloat3 vp=fgeo::VecOrthogonal(vcen,radius);
  tfloat3 p1b=TFloat3(cen1.x+vp.x,cen1.y+vp.y,cen1.z+vp.z);
  tfloat3 p2b=TFloat3(cen2.x+vp.x,cen2.y+vp.y,cen2.z+vp.z);
  tfloat3* pos=new tfloat3[56];
  pos[0]=p1b; pos[nside]=p2b;
  JMatrix4f mt=JMatrix4f::MatrixRot(360.f/28,cen1,cen2);
  for(int c=1;c<nside;c++){
    p1b=mt.MulPoint(p1b);
    p2b=mt.MulPoint(p2b);
    pos[c]=p1b; pos[nside+c]=p2b;
  }
  for(int c=0;c<nside;c++){
    AddTriangleWire(pos[c],pos[(c+1)%nside],cen1,v);
    AddTriangleWire(pos[nside+c],cen2,pos[nside+((c+1)%nside)],v);
    AddQuadWire(pos[(c+1)%nside],pos[c],pos[nside+c],pos[nside+((c+1)%nside)],v);
  }
  delete[] pos;
}
//==============================================================================
/// Add cylinder (wire mode).
//==============================================================================
void JSpVtkShape::AddCylinderWire(const tdouble3& cen1,const tdouble3& cen2
  ,double radius,word v)
{
  AddCylinderWire(ToTFloat3(cen1),ToTFloat3(cen2),float(radius),v);
}

//==============================================================================
/// Add sphere.
//==============================================================================
void JSpVtkShape::AddSphere(const tfloat3& pcen,float radius,word v){
  const int nside=28,nsidez=14;
  JMatrix4f mtz=JMatrix4f::MatrixRot(360.f/nside,pcen,TFloat3(pcen.x,pcen.y,pcen.z+1));
  JMatrix4f mtx=JMatrix4f::MatrixRot(180.f/(nsidez*2),pcen,TFloat3(pcen.x+1,pcen.y,pcen.z));
  tfloat3* pos=new tfloat3[nside*2];
  tfloat3* vpt=pos;
  tfloat3* vpt2=pos+nside;
  tfloat3 ptz1=TFloat3(pcen.x,pcen.y,pcen.z-radius);
  tfloat3 ptz2=TFloat3(pcen.x,pcen.y,pcen.z+radius);
  for(int c=0;c<nside;c++)vpt[c]=ptz1;
  for(int cz=1-nsidez;cz<=nsidez;cz++){
    if(cz==nsidez)for(int c=0;c<nside;c++)vpt2[c]=ptz2;
    else{
      vpt2[0]=mtx.MulPoint(vpt[0]);
      for(int c=1;c<nside;c++)vpt2[c]=mtz.MulPoint(vpt2[c-1]);
    }
    if(cz==1-nsidez)for(int c=0;c<nside;c++)AddTriangle(vpt[c],vpt2[c],vpt2[(c+1)%nside],v);
    else if(cz==nsidez)for(int c=0;c<nside;c++)AddTriangle(vpt[(c+1)%nside],vpt[c],vpt2[c],v);
    else for(int c=0;c<nside;c++)AddQuad(vpt[(c+1)%nside],vpt[c],vpt2[c],vpt2[(c+1)%nside],v);
    std::swap(vpt,vpt2);
  }
  delete[] pos;
}

//==============================================================================
/// Add sphere.
//==============================================================================
void JSpVtkShape::AddSphere(const tdouble3& pcen,double radius,word v){
  AddSphere(ToTFloat3(pcen),float(radius),v);
}

//==============================================================================
/// Add sphere (wire mode).
//==============================================================================
void JSpVtkShape::AddSphereWire(const tfloat3& pcen,float radius,word v){
  const int nside=28,nsidez=14;
  JMatrix4f mtz=JMatrix4f::MatrixRot(360.f/nside,pcen,TFloat3(pcen.x,pcen.y,pcen.z+1));
  JMatrix4f mtx=JMatrix4f::MatrixRot(180.f/(nsidez*2),pcen,TFloat3(pcen.x+1,pcen.y,pcen.z));
  tfloat3* pos=new tfloat3[nside*2];
  tfloat3* vpt=pos;
  tfloat3* vpt2=pos+nside;
  tfloat3 ptz1=TFloat3(pcen.x,pcen.y,pcen.z-radius);
  tfloat3 ptz2=TFloat3(pcen.x,pcen.y,pcen.z+radius);
  for(int c=0;c<nside;c++)vpt[c]=ptz1;
  for(int cz=1-nsidez;cz<=nsidez;cz++){
    if(cz==nsidez)for(int c=0;c<nside;c++)vpt2[c]=ptz2;
    else{
      vpt2[0]=mtx.MulPoint(vpt[0]);
      for(int c=1;c<nside;c++)vpt2[c]=mtz.MulPoint(vpt2[c-1]);
    }
    if(cz==1-nsidez)for(int c=0;c<nside;c++)AddTriangleWire(vpt[c],vpt2[c],vpt2[(c+1)%nside],v);
    else if(cz==nsidez)for(int c=0;c<nside;c++)AddTriangleWire(vpt[(c+1)%nside],vpt[c],vpt2[c],v);
    else for(int c=0;c<nside;c++)AddQuadWire(vpt[(c+1)%nside],vpt[c],vpt2[c],vpt2[(c+1)%nside],v);
    std::swap(vpt,vpt2);
  }
  delete[] pos;
}

//==============================================================================
/// Add sphere (wire mode).
//==============================================================================
void JSpVtkShape::AddSphereWire(const tdouble3& pcen,double radius,word v){
  AddSphereWire(ToTFloat3(pcen),float(radius),v);
}

//==============================================================================
/// Add spring.
//==============================================================================
void JSpVtkShape::AddSpring(const tfloat3& point1,const tfloat3& point2
  ,float restlength,float scalesize,float radius,float revlength,word v)
{
  const float cornersout=radius/2.f,cornersin=radius/4.f;
  const double ss=scalesize;
  const tdouble3 pt0=ToTDouble3(point1);
  const tdouble3 ptf=ToTDouble3(point2);
  const double scornerout=ss*cornersout;
  if(fgeo::PointsDist(pt0,ptf)>2.*ss*(cornersin+scornerout)){
    const tdouble3 vu=fgeo::VecUnitary(ptf-pt0);
    const tdouble3 pt1=pt0+(vu*(ss*cornersout));
    const tdouble3 pt2=ptf-(vu*(ss*cornersout));
    tdouble3 vr=fgeo::VecOrthogonal(vu,ss*radius);
    const unsigned nc=100;
    const JMatrix4d mr=JMatrix4d::MatrixRot(360./nc,pt1-pt1,pt2-pt1);
    tdouble3 vr2=vr;
    for(unsigned c=1;c<nc;c++){
      vr2=mr.MulPoint(vr2);
      if(vr.x>vr2.x || (vr.x==vr2.x && vr.y>vr2.y) || (vr.x==vr2.x && vr.y==vr2.y && vr.z>vr2.z))vr=vr2;
    }
    const JMatrix4d mt=JMatrix4d::MatrixRot(360./28,pt1,pt2);
    const double rlen=fgeo::PointsDist(pt1,pt2);
    const double restlen=restlength-scornerout*2.;
    const double rl=(rlen/restlen);
    const double rss=ss*rl;
    const double restlenin=double(restlen)-(ss*2.*cornersin);
    unsigned nrev=unsigned(restlenin/(ss*revlength))+1;
    unsigned nstep=nrev*28;
    const double sstep=rl*(restlenin/nstep);
    tdouble3 ptlast=pt0;
    tdouble3 pt=pt1;
    AddLine(ptlast,pt,v); ptlast=pt;
    pt=ptlast+(vu*(rss*cornersin))+vr;
    AddLine(ptlast,pt,v); ptlast=pt;
    for(unsigned c=1;c<=nstep;c++){
      pt=mt.MulPoint(pt)+(vu*sstep);
      AddLine(ptlast,pt,v); ptlast=pt;
    }
    pt=pt2;
    AddLine(ptlast,pt,v); ptlast=pt;
    pt=ptf;
    AddLine(ptlast,pt,v); ptlast=pt;
  }
  else AddLine(point1,point2,v);
}
//==============================================================================
/// Add spring.
//==============================================================================
void JSpVtkShape::AddSpring(const tdouble3& point1,const tdouble3& point2
  ,double restlength,double scalesize,double radius,double revlength,word v)
{
  AddSpring(ToTFloat3(point1),ToTFloat3(point2),float(restlength),float(scalesize)
    ,float(radius),float(revlength),v);
}

//==============================================================================
// Add cross.
//==============================================================================
void JSpVtkShape::AddCross(const tfloat3& pcen,float size,word v){
  AddLine(TFloat3(pcen.x-size,pcen.y,pcen.z),TFloat3(pcen.x+size,pcen.y,pcen.z),v);
  AddLine(TFloat3(pcen.x,pcen.y-size,pcen.z),TFloat3(pcen.x,pcen.y+size,pcen.z),v);
  AddLine(TFloat3(pcen.x,pcen.y,pcen.z-size),TFloat3(pcen.x,pcen.y,pcen.z+size),v);
}
//==============================================================================
// Add cross.
//==============================================================================
void JSpVtkShape::AddCross(const tdouble3& pcen,double size,word v){
  AddCross(ToTFloat3(pcen),float(size),v);
}

//==============================================================================
// Add circle.
//==============================================================================
void JSpVtkShape::AddCircle(const tfloat3& pcen,const tfloat3& vec,float radius
  ,word v)
{
  tfloat3 vp=fgeo::VecOrthogonal(vec,radius);
  tfloat3 pr=pcen+vp;
  tfloat3 pr2=pr;
  const JMatrix4f mt=JMatrix4f::MatrixRot(360.f/28,pcen,pcen+vec);
  for(int c=1;c<28;c++){
    tfloat3 pr3=mt.MulPoint(pr2);
    AddLine(pr2,pr3);
    pr2=pr3;
  }
  AddLine(pr2,pr);
}
//==============================================================================
// Add circle.
//==============================================================================
void JSpVtkShape::AddCircle(const tdouble3& pcen,const tdouble3& vec
  ,double radius,word v)
{
  AddCircle(ToTFloat3(pcen),ToTFloat3(vec),float(radius),v);
}




