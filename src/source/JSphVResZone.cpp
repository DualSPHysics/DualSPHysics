//HEAD_DSPH
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

/// \file JSphVResZone.cpp \brief Implements the class \ref JSphVResZone.

#include "JSphVResZone.h"

#include <cmath>

#include "Functions.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JXml.h"

JSphVResZone::JSphVResZone(bool cpu,const StCteSph &csp,bool inner,unsigned zone
  ,tdouble3 boxlimitmininner,tdouble3 boxlimitmaxinner
  ,tdouble3 boxlimitminouter,tdouble3 boxlimitmaxouter
  ,tdouble3 boxlimitminmid,tdouble3 boxlimitmaxmid
  ,tdouble3 maprealposmin,tdouble3 maprealposmax
  ,bool trackingisactive,unsigned trackingmk,JMatrix4d mat,bool issimple)
  : Log(AppInfo.LogPtr()),Cpu(cpu),CSP(csp),Inner(inner),Zone(zone+1)
  ,BoxLimitMinInner(boxlimitmininner),BoxLimitMaxInner(boxlimitmaxinner)
  ,BoxLimitMinOuter(boxlimitminouter),BoxLimitMaxOuter(boxlimitmaxouter)
  ,BoxLimitMinMid(boxlimitminmid),BoxLimitMaxMid(boxlimitmaxmid)
  ,MapRealPosMin(maprealposmin),MapRealPosMax(maprealposmax)
  ,TrackingisActive(trackingisactive),Mat(mat),IsSimple(issimple)
{
  ClassName = "JSphVResZone";

  if(!IsSimple){
    tmatrix4d mat=Mat.GetMatrix4d();
    tmatrix4d mat_inv= TMatrix4d(mat.a11,mat.a21,mat.a31,0,mat.a12,mat.a22,mat.a32,0,mat.a13,mat.a23,mat.a33,0,0,0,0,1);
    mat_inv.a14=-(mat.a14*mat.a11+mat.a24*mat.a21+mat.a34*mat.a31);
    mat_inv.a24=-(mat.a14*mat.a12+mat.a24*mat.a22+mat.a34*mat.a32);
    mat_inv.a34=-(mat.a14*mat.a13+mat.a24*mat.a23+mat.a34*mat.a33);
    BoxLimitMinInner= MatrixMulPoint(mat_inv,BoxLimitMinInner);       
    BoxLimitMaxInner= MatrixMulPoint(mat_inv,BoxLimitMaxInner);
    BoxLimitMinOuter= MatrixMulPoint(mat_inv,BoxLimitMinOuter);       
    BoxLimitMaxOuter= MatrixMulPoint(mat_inv,BoxLimitMaxOuter);      
    }

    BoxLimitMin=(Inner? BoxLimitMinInner : BoxLimitMinOuter);
    BoxLimitMax=(Inner? BoxLimitMaxInner : BoxLimitMaxOuter);
    BoxSize=(Inner? (BoxLimitMaxInner-BoxLimitMinInner) : (BoxLimitMaxOuter-BoxLimitMinOuter));
    Width = BoxLimitMaxOuter.x-BoxLimitMaxInner.x;
    Origin= (BoxLimitMaxInner+BoxLimitMinInner)/2.0;
    if(TrackingisActive) TrackingMk=trackingmk;

    Config();    
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphVResZone::~JSphVResZone()
{
  DestructorActive = true;
  Reset();
}

void JSphVResZone::Config(){    
    
  ComputeInterface();
  Ntot = static_cast<int>(Points.size());
  if(!IsSimple){
    tmatrix4d mat=Mat.GetMatrix4d();
    tmatrix4d mat_inv= TMatrix4d(mat.a11,mat.a21,mat.a31,0
      ,mat.a12,mat.a22,mat.a32,0,mat.a13,mat.a23,mat.a33,0,0,0,0,1);
    JMatrix4d matrot= JMatrix4d(mat_inv);
    Mat.MulArray(Ntot,Points.data());
    Mat.MulArray(Ntot,Normals.data());        
  }    
}

bool JSphVResZone::InZoneBoxOuter(const tdouble3 &ps) const {
    return (BoxLimitMinOuter.x <= ps.x && ps.x <= BoxLimitMaxOuter.x && BoxLimitMinOuter.y <= ps.y && ps.y <= BoxLimitMaxOuter.y && BoxLimitMinOuter.z <= ps.z && ps.z <= BoxLimitMaxOuter.z);
}
bool JSphVResZone::InZoneBoxMid(const tdouble3 &ps) const {
    return ((BoxLimitMinMid.x <= ps.x && ps.x <= BoxLimitMaxMid.x && BoxLimitMinMid.y <= ps.y && ps.y <= BoxLimitMaxMid.y && BoxLimitMinMid.z <= ps.z && ps.z <= BoxLimitMaxMid.z) != Inner);
}
bool JSphVResZone::InZoneBoxInner(const tdouble3 &ps) const {
    return (BoxLimitMinInner.x <= ps.x && ps.x <= BoxLimitMaxInner.x && BoxLimitMinInner.y <= ps.y && ps.y <= BoxLimitMaxInner.y && BoxLimitMinInner.z <= ps.z && ps.z <= BoxLimitMaxInner.z);
}
bool JSphVResZone::InZoneBox(const tdouble3 &ps) const {
    return (InZoneBoxOuter(ps) && !InZoneBoxInner(ps));
}

bool JSphVResZone::is_Out(const tdouble3 &ps) const {
    return ((!InZoneBoxOuter(ps) && Inner) || (InZoneBoxInner(ps) && !Inner));
}

bool JSphVResZone::is_Normal(const tdouble3 &ps) const {
    return ((InZoneBoxInner(ps) && Inner) || (!InZoneBoxOuter(ps) && !Inner));
}


void JSphVResZone::Reset() {
    Points.clear();
    Normals.clear();
}



void JSphVResZone::ComputeInterface() {
    
  tdouble3 length = (Inner ? (BoxLimitMaxOuter - BoxLimitMinOuter) : (BoxLimitMaxInner - BoxLimitMinInner));
  NPoints = TUint3(static_cast<unsigned>(std::round(length.x / CSP.dp)),static_cast<unsigned>(std::round(length.y / CSP.dp)), static_cast<unsigned>(std::round(length.z / CSP.dp)));
  tdouble3 shift = TDouble3(length.x - CSP.dp * (NPoints.x - 1), length.y - CSP.dp * (NPoints.y-1), length.z - CSP.dp * (NPoints.z - 1)) / 2.0;

  if(CSP.simulate2d){
    NPoints.y=0;
    shift.y=0;
  }

    // I don't remember what is does, I just know it is looping through six cube's faces!
  for (unsigned m=0;m<6;m++){
    int col[3]=     {(m%3)!=0,(m%3)!=1,(m%3)!=2};
    int colinv[3]=  {(m%3)==0,(m%3)==1,(m%3)==2};
    
    double shiftS[3] = {shift.x*col[0],shift.y*col[1],shift.z*col[2]};
    tdouble3 pointmin = (Inner ? BoxLimitMinOuter : BoxLimitMinInner);
    tdouble3 startpoint=pointmin;
    if(m>2)  startpoint=pointmin+TDouble3(double(colinv[0])*length.x,double(colinv[1])*length.y,double(colinv[2])*length.z);

    unsigned np[3]={NPoints.x*col[0],NPoints.y*col[1],NPoints.z*col[2]};

    for(int i=0;i<3;i++) np[i]=(np[i]!=0 ? np[i] : 1);
    
    if(CSP.simulate2d&&(m%3)==1)
            np[1] = 0;

    for(unsigned i=0;i<np[0];i++)
      for(unsigned j=0;j<np[1];j++)
        for(unsigned k=0;k<np[2];k++){
          tdouble3 point = TDouble3(startpoint.x + i * CSP.dp + shiftS[0], startpoint.y + j * CSP.dp + shiftS[1], startpoint.z + k * CSP.dp + shiftS[2]);
          tdouble3 normal = TDouble3(colinv[0], colinv[1], colinv[2]);

          normal = normal + (Inner==(m<3) ? normal*(-2.0) : TDouble3(0.0));

          // if(point>=MapRealPosMin && point<=MapRealPosMax){
            Points.push_back(point);
            Normals.push_back(normal);
          // }
        }
  }
}



