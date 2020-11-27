//HEAD_DSCODES
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

/// \file JOutputCsv.cpp \brief Implements the class \ref JOutputCsv.

#include "JOutputCsv.h"
#include "JDataArrays.h"
#include "Functions.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using namespace std;

//##############################################################################
//# JOutputCsv
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JOutputCsv::JOutputCsv(bool csvsepcoma,bool createpath)
  :CsvSepComa(csvsepcoma),CreatPath(createpath)
{
  ClassName="JOutputCsv";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JOutputCsv::~JOutputCsv(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialization of variables.
//==============================================================================
void JOutputCsv::Reset(){
  FileName="";
}

//==============================================================================
/// Stores data in CSV format.
//==============================================================================
void JOutputCsv::SaveCsv(std::string fname,const JDataArrays &arrays,std::string head){
  if(fun::GetExtension(fname).empty())FileName=fname=fun::AddExtension(fname,".csv");
  const char csvsep=(CsvSepComa? ',': ';');
  const unsigned nf=arrays.Count();
  const unsigned nv=arrays.GetDataCount(true);
  if(nv!=arrays.GetDataCount(false))Run_ExceptioonFile("The number of values in arrays is not the same.",fname);
  if(CreatPath)fun::MkdirPath(fun::GetDirParent(fname));
  ofstream pf;
  pf.open(fname.c_str());
  if(pf){
    //-Saves lines of head.
    while(!head.empty())pf << fun::StrSplit("\n",head) << endl;
    //-Saves head.
    for(unsigned cf=0;cf<nf;cf++){
      const string keyname=arrays.GetArrayCte(cf).keyname;
      const string units=arrays.GetArrayUnits(cf);
      const int dim=arrays.GetArrayDim(cf);
      if(dim==1)pf << keyname << units << csvsep;
      else if(dim==3)for(unsigned c=0;c<3;c++)pf << keyname << (c? (c==2? ".z": ".y"): ".x") << units << csvsep;
      else Run_ExceptioonFile(fun::PrintStr("Dimension %d of array \'%s\' is invalid.",dim,keyname.c_str()),fname);
    }
    pf << endl;
    //-Creates vector with output format of arrays.
    std::vector<std::string> outfmt;
    for(unsigned cf=0;cf<nf;cf++){
      string fmt=arrays.GetArrayFmt(cf);
      const int dim=arrays.GetArrayDim(cf);
      if(dim==3)fmt=fmt+csvsep+fmt+csvsep+fmt;
      else if(dim!=1)Run_ExceptioonFile(fun::PrintStr("Dimension %d of array \'%s\' is invalid.",dim,arrays.GetArrayCte(cf).keyname.c_str()),fname);
      outfmt.push_back(fmt);
    }
    //-Saves data.
    for(unsigned cv=0;cv<nv;cv++){
      for(unsigned cf=0;cf<nf;cf++){
        const JDataArrays::StDataArray& ar=arrays.GetArrayCte(cf);
        string fmt=outfmt[cf];
        string values;
        switch(ar.type){
          case TypeUchar:  { const byte     *v=(byte    *)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv])+csvsep; }break;
          case TypeUshort: { const word     *v=(word    *)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv])+csvsep; }break;
          case TypeUint:   { const unsigned *v=(unsigned*)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv])+csvsep; }break;
          case TypeFloat:  { const float    *v=(float   *)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv])+csvsep; }break;
          case TypeDouble: { const double   *v=(double  *)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv])+csvsep; }break;
          case TypeUint3:  { const tuint3   *v=(tuint3  *)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv].x,v[cv].y,v[cv].z)+csvsep; }break;
          case TypeFloat3: { const tfloat3  *v=(tfloat3 *)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv].x,v[cv].y,v[cv].z)+csvsep; }break;
          case TypeDouble3:{ const tdouble3 *v=(tdouble3*)ar.ptr;  values=fun::PrintStr(fmt.c_str(),v[cv].x,v[cv].y,v[cv].z)+csvsep; }break;
          default: Run_ExceptioonFile(fun::PrintStr("Type of array \'%s\' is invalid.",TypeToStr(ar.type)),fname);
        }
        pf << values;
      }
      pf << endl;
    }
    if(pf.fail())Run_ExceptioonFile("File writing failure.",fname);
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",fname);
}

//==============================================================================
/// Calculates statistic information of unitary arrays.
//==============================================================================
template<typename T> void JOutputCsv::CalculateStatsArray1(unsigned ndata,T *data
    ,double &valmin,double &valmax,double &valmean)const
{
  double vmin=0,vmax=0,vmea=0;
  if(ndata)vmin=vmax=double(data[0]);
  for(unsigned p=0;p<ndata;p++){
    const double v=double(data[p]);
    vmin=(vmin<v? vmin: v);
    vmax=(vmax>v? vmax: v);
    vmea=vmea+v;
  }
  vmea=vmea/double(ndata);
  //-Saves results.
  valmin=vmin; valmax=vmax; valmean=vmea;
}

//==============================================================================
/// Calculates statistic information of triple arrays.
//==============================================================================
template<typename T> void JOutputCsv::CalculateStatsArray3(unsigned ndata,T *data
  ,tdouble4 &valmin,tdouble4 &valmax,tdouble4 &valmean)const
{
  tdouble4 vmin=TDouble4(0);
  tdouble4 vmax=TDouble4(0);
  tdouble4 vmea=TDouble4(0);
  if(ndata){
    const T vv=data[0];
    const double vx=double(vv.x);
    const double vy=double(vv.y);
    const double vz=double(vv.z);
    const tdouble4 v=TDouble4(vx,vy,vz,sqrt(vx*vx+vy*vy+vz*vz));
    vmin=vmax=v;
  }
  for(unsigned p=0;p<ndata;p++){
    const T vv=data[p];
    const double vx=double(vv.x);
    const double vy=double(vv.y);
    const double vz=double(vv.z);
    const tdouble4 v=TDouble4(vx,vy,vz,sqrt(vx*vx+vy*vy+vz*vz));
    vmin=MinValues(vmin,v);
    vmax=MaxValues(vmax,v);
    vmea=vmea+v;
  }
  vmea=vmea/TDouble4(ndata);
  //-Saves results.
  valmin=vmin; valmax=vmax; valmean=vmea;
}

////==============================================================================
///// Compute and stores statistic information in CSV format.
////==============================================================================
//void JOutputCsv::SaveStatsCsv(std::string fname,bool create,int part,double timestep,const JDataArrays &arrays,std::string head){
//  if(fun::GetExtension(fname).empty())FileName=fname=fun::AddExtension(fname,".csv");
//  const char csvsep=(CsvSepComa? ',': ';');
//  const unsigned nf=arrays.GetCount();
//  const unsigned nv=arrays.GetDataCount();
//  if(CreateDirs)fun::MkdirPath(fun::GetDirParent(fname));
//  ofstream pf;
//  if(create)pf.open(fname.c_str());
//  else pf.open(fname.c_str(),ios_base::app);
//  if(pf){
//    if(create){
//      //-Saves lines of head.
//      while(!head.empty())pf << fun::StrSplit("\n",head) << endl;
//      //-Saves head.
//      pf << fun::StrCsvSep(CsvSepComa,"Part;Time [s];Count;");
//      for(unsigned cf=0;cf<nf;cf++){
//        const JDataArray* ar=arrays.GetArray(cf);
//        if(ar->GetPointer() && !ar->GetHidden()){
//          const string units=ar->GetUnits();
//          const string name=ar->GetName();
//          if(fun::StrLower(name)!="pos"){
//            pf << name << "_min"  << units << csvsep;
//            pf << name << "_max"  << units << csvsep;
//            pf << name << "_mean" << units << csvsep;
//          }
//          if(ar->IsTriple()){
//            pf << name << "_Xmin"  << units << csvsep;
//            pf << name << "_Xmax"  << units << csvsep;
//            pf << name << "_Xmean" << units << csvsep;
//            pf << name << "_Ymin"  << units << csvsep;
//            pf << name << "_Ymax"  << units << csvsep;
//            pf << name << "_Ymean" << units << csvsep;
//            pf << name << "_Zmin"  << units << csvsep;
//            pf << name << "_Zmax"  << units << csvsep;
//            pf << name << "_Zmean" << units << csvsep;
//          }
//        }
//      }
//      pf << endl;
//    }
//    //-Saves data.
//    pf << fun::PrintStr("%u",part)+csvsep;
//    pf << fun::PrintStr("%20.12E",timestep)+csvsep;
//    pf << fun::PrintStr("%u",nv)+csvsep;
//    for(unsigned cf=0;cf<nf;cf++){
//      const JDataArray* ar=arrays.GetArray(cf);
//      if(ar->GetPointer() && !ar->GetHidden()){
//        string fmt=JDataArraysDef::GetFmt(JDataArraysDef::TpDouble);
//        fmt=fmt+csvsep+fmt+csvsep+fmt;
//        double vmin,vmax,vmea;
//        tdouble4 vmin3,vmax3,vmea3;
//        switch(ar->GetType()){
//          case JDataArraysDef::TpChar:   { const char     *v=ar->GetPointerChar();     CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpUchar:  { const byte     *v=ar->GetPointerUchar();    CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpShort:  { const short    *v=ar->GetPointerShort();    CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpUshort: { const word     *v=ar->GetPointerUshort();   CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpInt:    { const int      *v=ar->GetPointerInt();      CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpUint:   { const unsigned *v=ar->GetPointerUint();     CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpLlong:  { const llong    *v=ar->GetPointerLlong();    CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpUllong: { const ullong   *v=ar->GetPointerUllong();   CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpFloat:  { const float    *v=ar->GetPointerFloat();    CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpDouble: { const double   *v=ar->GetPointerDouble();   CalculateStatsArray1(nv,v,vmin,vmax,vmea); }break;
//          case JDataArraysDef::TpInt3:   { const tint3    *v=ar->GetPointerInt3();     CalculateStatsArray3(nv,v,vmin3,vmax3,vmea3); }break;
//          case JDataArraysDef::TpUint3:  { const tuint3   *v=ar->GetPointerUint3();    CalculateStatsArray3(nv,v,vmin3,vmax3,vmea3); }break;
//          case JDataArraysDef::TpFloat3: { const tfloat3  *v=ar->GetPointerFloat3();   CalculateStatsArray3(nv,v,vmin3,vmax3,vmea3); }break;
//          case JDataArraysDef::TpDouble3:{ const tdouble3 *v=ar->GetPointerDouble3();  CalculateStatsArray3(nv,v,vmin3,vmax3,vmea3); }break;
//          default: Run_Exceptioon("Type of array is invalid.");
//        }
//        if(!ar->IsTriple())pf << fun::PrintStr(fmt.c_str(),vmin,vmax,vmea)+csvsep;
//        else{
//          if(fun::StrLower(ar->GetName())!="pos"){
//            pf << fun::PrintStr(fmt.c_str(),vmin3.w,vmax3.w,vmea3.w)+csvsep;
//          }
//          pf << fun::PrintStr(fmt.c_str(),vmin3.x,vmax3.x,vmea3.x)+csvsep;
//          pf << fun::PrintStr(fmt.c_str(),vmin3.y,vmax3.y,vmea3.y)+csvsep;
//          pf << fun::PrintStr(fmt.c_str(),vmin3.z,vmax3.z,vmea3.z)+csvsep;
//        }
//      }
//    }
//    pf << endl;
//    if(pf.fail())Run_ExceptioonFile("File writing failure.",fname);
//    pf.close();
//  }
//  else Run_ExceptioonFile("Cannot open the file.",fname);
//}




