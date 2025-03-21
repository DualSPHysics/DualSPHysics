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

/// \file JSpVtkData.cpp \brief Implements the class \ref JSpVtkData.

#include "JSpVtkData.h"
#include "Functions.h"
#include "JDataArrays.h"

#include <fstream>
#include <cstring>

using namespace std;

//##############################################################################
//# JSpVtkData
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSpVtkData::JSpVtkData(){
  ClassName="JSpVtkData";
}

//==============================================================================
/// Destructor.
//==============================================================================
JSpVtkData::~JSpVtkData(){
  DestructorActive=true;
}

//==============================================================================
/// Saves array data in VTK file.
//==============================================================================
void JSpVtkData::SaveVtk(std::string file,const JDataArrays& arrays
  ,std::string posfield)
{
  if(fun::GetByteOrder()!=fun::LittleEndian)
    Run_ExceptioonFile("Big-Endian mode is not supported for now.",file);
  const unsigned np=arrays.GetDataCount(true);
  if(np!=arrays.GetDataCount(false))
    Run_ExceptioonFile("The number of values in arrays does not match.",file);
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
      tfloat3* vp=new tfloat3[np];
      const JDataArrays::StDataArray& arpos=arrays.GetArrayCte(posfield);
      if(arpos.type==TypeFloat3)memcpy(vp,arpos.ptr,sizeof(tfloat3)*np);
      else if(arpos.type==TypeDouble3){
        const tdouble3* pos=(tdouble3*)arpos.ptr;
        for(unsigned p=0;p<np;p++)vp[p]=ToTFloat3(pos[p]);
      }
      pf << "POINTS " << np << " float" << endl;
      const bool lite=(fun::GetByteOrder()==fun::LittleEndian);
      if(lite)fun::SetByteOrder32((int*)vp,np*3);
      pf.write((const char*)vp,sizeof(tfloat3)*np);
      pf << endl;
      pf << "VERTICES " << np << " " << np*2 << endl;
      unsigned* data=(unsigned*)vp;
      for(unsigned p=0,pp=0;p<np;p++){
        data[pp++]=1;
        data[pp++]=p;
      }
      if(lite)fun::SetByteOrder32((int*)data,np*2);
      pf.write((const char*)data,sizeof(int)*np*2);
      pf << endl;
      const unsigned nd=arrays.Count();
      unsigned nf=0;
      for(unsigned cd=0;cd<nd;cd++){
        const JDataArrays::StDataArray& ar=arrays.GetArrayCte(cd);
        if(ar.ptr && ar.keyname!=posfield){
          if(!nf)pf << "POINT_DATA " << np << endl;
          const TpTypeData tp=ar.type,tpa=TypeParent(ar.type);
          const unsigned stp=TypeSizeU(tp),stpa=TypeSizeU(tpa);
          const int dim=TypeDim(tp);
          if(dim!=1 && dim!=3)Run_ExceptioonFile(fun::PrintStr(
            "Dimension of \'%s\' is invalid.",ar.keyname.c_str()),file);
          const unsigned npa=np*dim;
          string tpstr=(stp==1 || stp==2? " short": (tpa==TypeUint || tpa==TypeInt? " int": " float"));
          if(stp==1){
            short* dats=(short*)vp;
            const byte* ptr=(byte*)ar.ptr;
            for(unsigned p=0;p<np;p++)dats[p]=short(ptr[p]);
          }
          else if(stp==2)memcpy(vp,ar.ptr,sizeof(short)*np);
          else if(tpa==TypeUint || tpa==TypeInt)memcpy(vp,ar.ptr,sizeof(int)*npa);
          else if(tpa==TypeDouble){
            float* datf=(float*)vp;
            const double* ptr=(double*)ar.ptr;
            for(unsigned p=0;p<npa;p++)datf[p]=float(ptr[p]);
          }
          else if(tpa==TypeFloat)memcpy(vp,ar.ptr,sizeof(float)*npa);
          else Run_ExceptioonFile(fun::PrintStr(
            "Array type %s of \'%s\' is invalid.",TypeToStr(tp),ar.keyname.c_str()),file);
          if(lite){
            if(stpa<=2)fun::SetByteOrder16((short*)vp,npa);
            else       fun::SetByteOrder32((int*)vp,npa);
          } 

          if(!nf)  pf << "SCALARS " << ar.keyname << tpstr << (dim>1? " 3": "") << endl;
          if(nf==1)pf << "FIELD FieldData " << nd-2 << endl;
          if(!nf)pf << "LOOKUP_TABLE default" << endl;
          else pf << ar.keyname << " " << dim << " " << np << tpstr << endl;
          pf.write((const char*)vp,size_t(stpa<=2? 2: 4)*npa);
          pf << endl;
          nf++;
        }
      }
      delete[] vp;
    }
    if(pf.fail())Run_ExceptioonFile("File writing failure.",file);
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",file);
}

//==============================================================================
/// Saves array data in VTK file.
//==============================================================================
void JSpVtkData::Save(std::string file,const JDataArrays& arrays
  ,std::string posfield)
{
  JSpVtkData svtk;
  svtk.SaveVtk(file,arrays,posfield);
}
