/*
 * JSphBufferZone.cpp
 *
 *  Created on: Nov 9, 2021
 *      Author: francesco
 */
#include "JSphBufferZone.h"

#include <cmath>

#include "Functions.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JXml.h"

JSphBufferZone::JSphBufferZone(bool cpu,const StCteSph &csp,bool inner,unsigned zone,tdouble3 boxlimitmininner,tdouble3 boxlimitmaxinner
  ,tdouble3 boxlimitminouter,tdouble3 boxlimitmaxouter,tdouble3 boxlimitminmid
  ,tdouble3 boxlimitmaxmid,bool trackingisactive,unsigned trackingmk,JMatrix4d mat,bool issimple)
  : Log(AppInfo.LogPtr()),Cpu(cpu),CSP(csp),Inner(inner),Zone(zone+1),BoxLimitMinInner(boxlimitmininner),BoxLimitMaxInner(boxlimitmaxinner)
  ,BoxLimitMinOuter(boxlimitminouter),BoxLimitMaxOuter(boxlimitmaxouter),BoxLimitMinMid(boxlimitminmid),BoxLimitMaxMid(boxlimitmaxmid)
  ,TrackingisActive(trackingisactive),Mat(mat),IsSimple(issimple)
  {
  ClassName = "BufferZone";

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

void JSphBufferZone::Config(){    
    ComputeInterface();
    Ntot = PointsIn.size();
    // CalculateNpoints();
    AllocateMemory();
    CalculatePointsNormals();
    if(!IsSimple){
        tmatrix4d mat=Mat.GetMatrix4d();
        tmatrix4d mat_inv= TMatrix4d(mat.a11,mat.a21,mat.a31,0,mat.a12,mat.a22,mat.a32,0,mat.a13,mat.a23,mat.a33,0,0,0,0,1);
        JMatrix4d matrot= JMatrix4d(mat_inv);
        Mat.MulArray(Ntot,Points);
        Mat.MulArray(Ntot,Normals);        
    }    
}

bool JSphBufferZone::InZoneBoxOuter(const tdouble3 &ps) const {
    return (BoxLimitMinOuter.x <= ps.x && ps.x <= BoxLimitMaxOuter.x && BoxLimitMinOuter.y <= ps.y && ps.y <= BoxLimitMaxOuter.y && BoxLimitMinOuter.z <= ps.z && ps.z <= BoxLimitMaxOuter.z);
}
bool JSphBufferZone::InZoneBoxMid(const tdouble3 &ps) const {
    return ((BoxLimitMinMid.x <= ps.x && ps.x <= BoxLimitMaxMid.x && BoxLimitMinMid.y <= ps.y && ps.y <= BoxLimitMaxMid.y && BoxLimitMinMid.z <= ps.z && ps.z <= BoxLimitMaxMid.z) != Inner);
}
bool JSphBufferZone::InZoneBoxInner(const tdouble3 &ps) const {
    return (BoxLimitMinInner.x <= ps.x && ps.x <= BoxLimitMaxInner.x && BoxLimitMinInner.y <= ps.y && ps.y <= BoxLimitMaxInner.y && BoxLimitMinInner.z <= ps.z && ps.z <= BoxLimitMaxInner.z);
}
bool JSphBufferZone::InZoneBox(const tdouble3 &ps) const {
    return (InZoneBoxOuter(ps) && !InZoneBoxInner(ps));
}

void JSphBufferZone::CalculateNpoints() {
    tdouble3 length = (Inner ? (BoxLimitMaxOuter - BoxLimitMinOuter) : (BoxLimitMaxInner - BoxLimitMinInner));
    NPoints = TUint3(std::round(length.x / CSP.dp), std::round(length.y / CSP.dp), std::round(length.z / CSP.dp));
    Ntot = NPoints.x + NPoints.x + NPoints.y + NPoints.y + NPoints.z + NPoints.z;
}

void JSphBufferZone::AllocateMemory() {
    Points = new tdouble3[Ntot];
    Normals = new tdouble3[Ntot];
}

void JSphBufferZone::CalculatePointsNormals() {
  for (unsigned i = 0; i < Ntot; i++){
    Points[i] = PointsIn[i];
    Normals[i] = NormalIn[i];
  }
}

tdouble3 JSphBufferZone::getNormal(const tdouble3 &ps) {
    tdouble3 dis = ps - Origin;
    tdouble3 boxsize = (Inner ? BoxSize : BoxSize - Width);
    if (abs(dis.x) > boxsize.x / 2.0 && abs(dis.z) > boxsize.z / 2.0)
        return (TDouble3(1, 1, 1));
    else if (abs(dis.x) > boxsize.x / 2.0)
        return TDouble3(1.0, 0.0, 0.0);
    else if (boxsize.z / 2.0 < abs(dis.z))
        return TDouble3(0.0, 0.0, 1.0);
    else
        return TDouble3(0.0, 0.0, 0.0);
}

void JSphBufferZone::Reset() {
    delete[] Points;
    delete[] Normals;
    Normals = nullptr;
    Points = nullptr;
}

bool JSphBufferZone::is_Out(const tdouble3 &ps) const {
    return ((!InZoneBoxOuter(ps) && Inner) || (InZoneBoxInner(ps) && !Inner));
}
//
bool JSphBufferZone::is_Normal(const tdouble3 &ps) const {
    return ((InZoneBoxInner(ps) && Inner) || (!InZoneBoxOuter(ps) && !Inner));
}

void JSphBufferZone::ComputeInterface() {
    tdouble3 length = (Inner ? (BoxLimitMaxOuter - BoxLimitMinOuter) : (BoxLimitMaxInner - BoxLimitMinInner));
    NPoints = TUint3(std::round(length.x / CSP.dp), std::round(length.y / CSP.dp), std::round(length.z / CSP.dp));
    tdouble3 shift = TDouble3(length.x - CSP.dp * (NPoints.x - 1), length.y - CSP.dp * (NPoints.y-1), length.z - CSP.dp * (NPoints.z - 1)) / 2.0;

    if(CSP.simulate2d){
        NPoints.y=0;
        shift.y=0;
    }

    // outer loop
    for (unsigned m = 0; m < 6; m++) {
        // permutation column
        int col[3] = {(m % 3) != 0, (m % 3) != 1, (m % 3) != 2};
        int colinv[3] = {(m % 3) == 0, (m % 3) == 1, (m % 3) == 2};
        double shiftS[3] = {shift.x * col[0], shift.y * col[1], shift.z * col[2]};
        tdouble3 pointmin = (Inner ? BoxLimitMinOuter : BoxLimitMinInner);
        tdouble3 startpoint = pointmin;
        if (m > 2) startpoint = pointmin + TDouble3(double(colinv[0]) * length.x, double(colinv[1]) * length.y, double(colinv[2]) * length.z);

        unsigned np[3] = {NPoints.x * col[0], NPoints.y * col[1], NPoints.z * col[2]};

        for (int i = 0; i < 3; i++)
            np[i] = (np[i] != 0 ? np[i] : 1);
        if (CSP.simulate2d && (m % 3) == 1)
            np[1] = 0;

        for (int i = 0; i < np[0]; i++)
            for (int j = 0; j < np[1]; j++)
                for (int k = 0; k < np[2]; k++) {
                    tdouble3 point = TDouble3(startpoint.x + i * CSP.dp + shiftS[0], startpoint.y + j * CSP.dp + shiftS[1], startpoint.z + k * CSP.dp + shiftS[2]);
                    tdouble3 normal = TDouble3(colinv[0], colinv[1], colinv[2]);
                    // if (CSP.simulate2d) point.y = 0.0;

                    normal = normal + (Inner == (m < 3) ? normal * (-2.0) : TDouble3(0.0));
                    // if (CheckPointInterface(point))
                    // {
                    PointsIn.push_back(point);
                    NormalIn.push_back(normal);
                    // }
                }
    }
}



