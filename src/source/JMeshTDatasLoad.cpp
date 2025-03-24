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

/// \file JMeshTDatasLoad.cpp \brief Implements the classes JMeshTDatasLoad.

#include "JMeshTDatasLoad.h"
#include "JMeshData.h"
#include "JMeshTDatasSave.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JBinaryData.h"
#include "JReadDatafile.h"
#include "JDataArrays.h"
#include "JTimeControl.h"

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>

using namespace std;

/// Set of definitions, functions and clases for JMeshData codes.
namespace jmsh{

//##############################################################################
//# JMeshTDatasLoad
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JMeshTDatasLoad::JMeshTDatasLoad(double tmin,double tmax
  ,double looptmax,double looptbegin
  ,tdouble3 pmin,tdouble3 pmax
  ,double settime,tdouble3 setpos)
  :FilterTime(tmin!=DBL_MAX || tmax!=DBL_MAX || looptmax!=DBL_MAX)
  ,TimeMin(tmin!=DBL_MAX? tmin: (looptmax!=DBL_MAX? -DBL_MAX: DBL_MAX))
  ,TimeMax(min(tmax,looptmax))
  ,LoopTmaxRq(looptmax),LoopTbeginRq(looptbegin)
  ,FilterPos(pmin.x!=DBL_MAX),PosMin(pmin),PosMax(pmax)
  ,SetTime(settime),SetPos(setpos)
{
  ClassName="JMeshTDatasLoad";
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JMeshTDatasLoad::~JMeshTDatasLoad(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JMeshTDatasLoad::Reset(){
  FormatVer=FormatVerDef;
  FileIn="";
  LoopTmax=DBL_MAX; 
  LoopTsub=0;
}

//==============================================================================
/// Loads data definitions from binary file (.mbi4) and check selected data.
//==============================================================================
void JMeshTDatasLoad::LoadDataDefs(const JBinaryData* bdat,std::string varlist
  ,std::vector<StMeshData>& datadef)
{
  //-Prepares data selection.
  varlist=fun::StrWithoutChar(varlist,' ');
  vector<string> vnames;
  unsigned nnames=fun::VectorSplitStr(",",varlist,vnames);
  if(!varlist.empty())varlist=string(",")+fun::StrLower(varlist)+",";
  //-Loads data definitions and filter selected data in varlist.
  string varlist2;
  datadef.clear();
  const unsigned nd=bdat->GetvUint("Data_Count");
  for(unsigned c=0;c<nd;c++){
    const string name=bdat->GetvText(fun::PrintStr("DataName_%02u",c));
    if(varlist.empty() || int(varlist.find(string(",")+fun::StrLower(name)+","))>=0){
      varlist2=varlist2+string(",")+fun::StrLower(name);
      StMeshData dd;
      dd.name=name;
      dd.type=(TpTypeData)bdat->GetvUint(fun::PrintStr("DataType_%02u",c));
      dd.size=bdat->GetvUint(fun::PrintStr("DataSize_%02u",c));
      dd.tag =bdat->GetvInt (fun::PrintStr("DataTag_%02u" ,c));
      datadef.push_back(dd);
    }
  }
  //-Check requested data.
  for(unsigned c=0;c<nnames;c++)if(int(varlist2.find(fun::StrLower(vnames[c])))<0){
    Run_ExceptioonFile(fun::PrintStr("The requested data \'%s\' is missing.",vnames[c].c_str()),FileIn);
  }
  //-Check no data.
  if(varlist2.empty())Run_ExceptioonFile("No data in file.",FileIn);
}

//==============================================================================
/// Loads data from binary file (.mbi4).
//==============================================================================
void JMeshTDatasLoad::LoadFileBin(std::vector<JMeshData*>& vdata
  ,std::string varlist)
{
  JBinaryData bdat;
  bdat.LoadFile(FileIn,"JMeshTDatas");
  FormatVer=bdat.GetvUint("FormatVer");
  if(FormatVer<FormatVerDef)Run_ExceptioonFile(fun::PrintStr("Format version \'%u\' of file is not supported.",FormatVer),FileIn);
  //-Loads grid data definition.
  jmsh::StMeshPts m;
  memset(&m,0,sizeof(jmsh::StMeshPts));
  m.ptref =bdat.GetvDouble3("MeshPtRef");
  m.vdp1  =bdat.GetvDouble3("MeshVdp1");
  m.vdp2  =bdat.GetvDouble3("MeshVdp2");
  m.vdp3  =bdat.GetvDouble3("MeshVdp3");
  m.npt1  =bdat.GetvUint   ("MeshNpt1");
  m.npt2  =bdat.GetvUint   ("MeshNpt2");
  m.npt3  =bdat.GetvUint   ("MeshNpt3");
  m.npt   =m.npt1*m.npt2*m.npt3; 
  m.dirdat=bdat.GetvFloat3 ("MeshDirDat");
  m.dirdat=fgeo::VecUnitary(m.dirdat);
  //-Configuration for loop mode.
  const bool loopmode=(LoopTmaxRq!=DBL_MAX);
  const double looptbeg=LoopTbeginRq;
  unsigned loopct=0;
  double loopct_time=0;
  double loopct_tdif=DBL_MAX;
  //-Loads data definitions and filter selected data in varlist.
  std::vector<StMeshData> datadefs;
  LoadDataDefs(&bdat,varlist,datadefs);
  const unsigned ndatas=unsigned(datadefs.size());
  //-Loads datatimes.
  bdat.LoadFileListApp(FileIn,"JMeshTDatas",false);
  const size_t ntimes=bdat.GetItemsCount()-1;
  //printf("==> times:%u\n",ntimes);
  bool nextct=true;
  for(size_t ct=0;ct<ntimes && nextct;ct++){
    JBinaryData* ite=bdat.GetItem(ct+1);
    if(!ite->ExistsValue("Ntimes")){//-SigleData version.
      const double t=ite->GetvDouble("TimeStep");
      if(!FilterTime || (TimeMin<=t && t<=TimeMax)){
        ///printf("==>    [%d]:%f\n",ct,t);
        JMeshData* mdat=new JMeshData();
        mdat->ConfigMesh(m,t,"");
        for(unsigned cd=0;cd<ndatas;cd++){
          const StMeshData& d=datadefs[cd];
          switch(d.type){
            case TypeFloat:    ite->CopyArrayData(d.name,d.size,mdat->CreateArrayPtrFloat  (d.name,d.tag,d.size));   break;
            case TypeDouble:   ite->CopyArrayData(d.name,d.size,mdat->CreateArrayPtrDouble (d.name,d.tag,d.size));   break;
            case TypeFloat3:   ite->CopyArrayData(d.name,d.size,mdat->CreateArrayPtrFloat3 (d.name,d.tag,d.size));   break;
            case TypeDouble3:  ite->CopyArrayData(d.name,d.size,mdat->CreateArrayPtrDouble3(d.name,d.tag,d.size));   break;
            default: Run_ExceptioonFile(fun::PrintStr("Type of data \'%s\' is invalid.",TypeToStr(d.type)),FileIn);
          }
        }
        vdata.push_back(mdat);
        //-Configuration for loop mode.
        if(loopmode){
          const double tdif=fabs(t-looptbeg);
          if(tdif<loopct_tdif){
            loopct=unsigned(vdata.size())-1;
            loopct_time=t;
            loopct_tdif=tdif;
          }
        }
      }
      else if(t>TimeMax)nextct=false;
    }
    else{//-MultiData version.
      const unsigned nt2=ite->GetvUint("Ntimes");
      //printf("==> ct:%u/%u nt2:%u\n",ct,ntimes,nt2);
      double* vtimes=new double[nt2];
      ite->CopyArrayData("TimeSteps",nt2,vtimes);
      for(unsigned ct2=0;ct2<nt2;ct2++){
        const double t=vtimes[ct2];
        //printf("====> ct2:%u/%u time:%f\n",ct2,nt2,t);
        if(!FilterTime || (TimeMin<=t && t<=TimeMax)){
          ///printf("==>    [%d]:%f\n",ct,t);
          JMeshData* mdat=new JMeshData();
          mdat->ConfigMesh(m,t,"");
          for(unsigned cd=0;cd<ndatas;cd++){
            const StMeshData& d=datadefs[cd];
            switch(d.type){
              case TypeFloat:{
                JBinaryDataArray* ar=ite->CheckCopyArrayData(d.name,d.size*nt2,JBinaryDataDef::DatFloat);
                const float* vsrc=(const float*)ar->GetDataPointer();
                float* vdst=mdat->CreateArrayPtrFloat(d.name,d.tag,d.size);
                memcpy(vdst,vsrc+(d.size*ct2),sizeof(float)*d.size);
              }break;
              case TypeDouble:{
                JBinaryDataArray* ar=ite->CheckCopyArrayData(d.name,d.size*nt2,JBinaryDataDef::DatDouble);
                const double* vsrc=(const double*)ar->GetDataPointer();
                double* vdst=mdat->CreateArrayPtrDouble(d.name,d.tag,d.size);
                memcpy(vdst,vsrc+(d.size*ct2),sizeof(double)*d.size);
              }break;
              case TypeFloat3:{
                JBinaryDataArray* ar=ite->CheckCopyArrayData(d.name,d.size*nt2,JBinaryDataDef::DatFloat3);
                const tfloat3* vsrc=(const tfloat3*)ar->GetDataPointer();
                tfloat3* vdst=mdat->CreateArrayPtrFloat3(d.name,d.tag,d.size);
                memcpy(vdst,vsrc+(d.size*ct2),sizeof(tfloat3)*d.size);
              }break;
              case TypeDouble3:{
                JBinaryDataArray* ar=ite->CheckCopyArrayData(d.name,d.size*nt2,JBinaryDataDef::DatDouble3);
                const tdouble3* vsrc=(const tdouble3*)ar->GetDataPointer();
                tdouble3* vdst=mdat->CreateArrayPtrDouble3(d.name,d.tag,d.size);
                memcpy(vdst,vsrc+(d.size*ct2),sizeof(tdouble3)*d.size);
              }break;
              default: Run_ExceptioonFile(fun::PrintStr("Type of data \'%s\' is invalid.",TypeToStr(d.type)),FileIn);
            }
          }
          vdata.push_back(mdat);
          //-Configuration for loop mode.
          if(loopmode){
            const double tdif=fabs(t-looptbeg);
            if(tdif<loopct_tdif){
              loopct=unsigned(vdata.size())-1;
              loopct_time=t;
              loopct_tdif=tdif;
            }
          }
        }
      }
      delete[] vtimes; vtimes=NULL;
    }
  }
  //-Final configuration for loop mode.
  if(loopmode){
    const unsigned ndata=unsigned(vdata.size());
    if(loopct_tdif==DBL_MAX)Run_ExceptioonFile("The beginning time for loop configuration not found.",FileIn);
    if(loopct+1>=ndata)Run_ExceptioonFile("The final time is not higher than beginning time for loop configuration.",FileIn);
    if(ndata){
      const JMeshData* mdat0=vdata[loopct];
      JMeshData* mdat=vdata[ndata-1];
      const double tbeg=mdat0->GetTimeStep();
      LoopTmax=mdat->GetTimeStep();
      LoopTsub=LoopTmax-tbeg;
      mdat->CopyDataFrom(mdat0);
      mdat->SetTimeStep(LoopTmax);
    }
  }
}

//==============================================================================
/// Loads data from CSV file.
//==============================================================================
void JMeshTDatasLoad::LoadFileCsv(std::vector<JMeshData*>& vdata){
  JReadDatafile rdat;
  rdat.LoadFile(FileIn);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  //printf("=> Rows:%u\n",rows);
  if(rows<4)Run_ExceptioonFile("The file does not have enough lines.",FileIn);
  //-Checks format file.
  {
    string fmth=rdat.ReadNextValue(true);
    rdat.SetReadLine(1);
    string fmtv=rdat.ReadNextValue(true);
    //printf("=> [%s]:[%s]\n",fmth.c_str(),fmtv.c_str());
    const string fmt=fun::PrintStr("JMeshTDatas-%u",FormatVerDef);
    if(fmth!="Format" || fmtv!=fmt)Run_ExceptioonFile(fun::PrintStr("Format file \'%s=%s\' does not match \'Format=%s\'.",fmth.c_str(),fmtv.c_str(),fmt.c_str()),FileIn);
  }
  //-Loads data definition.
  const string datname =rdat.ReadNextValue(true);
  const string datunits=rdat.ReadNextValue(true);
  const string dattypex=fun::StrTrim(fun::StrLower(rdat.ReadNextValue(true)));
  const bool   dat12   =rdat.ReadNextBool(true);
  TpTypeData dattype=TypeNull;
       if(dattypex==TypeToStr(TypeFloat  ))dattype=TypeFloat;
  else if(dattypex==TypeToStr(TypeFloat3 ))dattype=TypeFloat3;
  else if(dattypex==TypeToStr(TypeDouble ))dattype=TypeDouble;
  else if(dattypex==TypeToStr(TypeDouble3))dattype=TypeDouble3;
  else Run_ExceptioonFile(fun::PrintStr("Type of data \'%s\' is invalid.",dattypex.c_str()),FileIn);
  //-Loads mesh definition.
  StMeshPts m;
  memset(&m,0,sizeof(StMeshPts));
  m.npt1  =rdat.ReadNextUnsigned(true);
  m.npt2  =rdat.ReadNextUnsigned(true);
  m.npt3  =rdat.ReadNextUnsigned(true);
  m.ptref =rdat.ReadNextDouble3(true);
  m.vdp1  =rdat.ReadNextDouble3(true);
  m.vdp2  =rdat.ReadNextDouble3(true);
  m.vdp3  =rdat.ReadNextDouble3(true);
  m.dirdat=rdat.ReadNextFloat3 (true);
  m.dirdat=fgeo::VecUnitary(m.dirdat);
  m.npt=m.npt1*m.npt2*m.npt3;
  //-Configuration for loop mode.
  const bool loopmode=(LoopTmaxRq!=DBL_MAX);
  const double looptbeg=LoopTbeginRq;
  unsigned loopct=0;
  double loopct_time=0;
  double loopct_tdif=DBL_MAX;
  //printf("=> [%s] npt:%u,%u,%u  ptref:(%g,%g,%g)  vdir:(%g,%g,%g)\n",datname.c_str(),m.npt1,m.npt2,m.npt3,m.ptref.x,m.ptref.y,m.ptref.z,m.dirdat.x,m.dirdat.y,m.dirdat.z);
  //-Loads data times.
  const int dattag=(dat12? JMeshData::TAGDATA12: 0);
  const unsigned datsize=(dat12? m.npt1*m.npt2: m.npt);
  const bool printtime=true;
  const bool tc_rows=(!FilterTime || TimeMax==DBL_MAX);
  JTimeControl tc(5,60);
  for(unsigned cr=3;cr<rows;cr++){
    rdat.SetReadLine(cr);
    const double t=rdat.ReadNextDouble(true);
    if(!FilterTime || (TimeMin<=t && t<=TimeMax)){
      //printf("==>    [%d]:%f\n",cr,t);
      JMeshData* mdat=new JMeshData();
      mdat->ConfigMesh(m,t);
      switch(dattype){
        case TypeFloat:{
          float* ptr=mdat->CreateArrayPtrFloat(datname,dattag,datsize);
          for(unsigned c=0;c<datsize;c++)*(ptr++)=rdat.ReadNextFloat(true);
        }break;
        case TypeFloat3:{
          tfloat3* ptr=mdat->CreateArrayPtrFloat3(datname,dattag,datsize);
          for(unsigned c=0;c<datsize;c++)*(ptr++)=rdat.ReadNextFloat3(true);
        }break;
        case TypeDouble:{
          double* ptr=mdat->CreateArrayPtrDouble(datname,dattag,datsize);
          for(unsigned c=0;c<datsize;c++)*(ptr++)=rdat.ReadNextDouble(true);
        }break;
        case TypeDouble3:{
          tdouble3* ptr=mdat->CreateArrayPtrDouble3(datname,dattag,datsize);
          for(unsigned c=0;c<datsize;c++)*(ptr++)=rdat.ReadNextDouble3(true);
        }break;
        default: Run_ExceptioonFile(fun::PrintStr("Type of data \'%s\' is invalid.",dattypex.c_str()),FileIn);
      }
      vdata.push_back(mdat);
      //-Configuration for loop mode.
      if(loopmode){
        const double tdif=fabs(t-looptbeg);
        if(tdif<loopct_tdif){
          loopct=unsigned(vdata.size())-1;
          loopct_time=t;
          loopct_tdif=tdif;
        }
      }
      if(printtime && tc.CheckTime()){
        const double tt=(tc_rows? double(cr)/double(rows+1): t/TimeMax);
        printf("  %s\n",tc.GetInfoFinish(tt).c_str());
      }
    }
  }
  //-Final configuration for loop mode.
  if(loopmode){
    const unsigned ndata=unsigned(vdata.size());
    if(loopct_tdif==DBL_MAX)Run_ExceptioonFile("The beginning time for loop configuration not found.",FileIn);
    if(loopct+1>=ndata)Run_ExceptioonFile("The final time is not higher than beginning time for loop configuration.",FileIn);
    if(ndata){
      const JMeshData* mdat0=vdata[loopct];
      JMeshData* mdat=vdata[ndata-1];
      const double tbeg=mdat0->GetTimeStep();
      LoopTmax=mdat->GetTimeStep();
      LoopTsub=LoopTmax-tbeg;
      mdat->CopyDataFrom(mdat0);
      mdat->SetTimeStep(LoopTmax);
    }
  }
}

//==============================================================================
/// Computes range of points according to PosMin and PosMax.
//==============================================================================
StMeshPts JMeshTDatasLoad::ComputeFilterPosRanges(StMeshPts mp,unsigned& c1ini
  ,unsigned& c1fin,unsigned& c2ini,unsigned& c2fin,unsigned& c3ini,unsigned& c3fin)const
{
  const unsigned npt1=mp.npt1;
  const unsigned npt2=mp.npt2;
  const unsigned npt3=mp.npt3;
  unsigned* vnum1=new unsigned[npt1];
  unsigned* vnum2=new unsigned[npt2];
  unsigned* vnum3=new unsigned[npt3];
  memset(vnum1,0,sizeof(unsigned)*npt1);
  memset(vnum2,0,sizeof(unsigned)*npt2);
  memset(vnum3,0,sizeof(unsigned)*npt3);
  for(unsigned c3=0;c3<npt3;c3++)for(unsigned c2=0;c2<npt2;c2++)for(unsigned c1=0;c1<npt1;c1++){
    const tdouble3 pt=mp.ptref+(mp.vdp1*c1)+(mp.vdp2*c2)+(mp.vdp3*c3);
    if(PosMin<=pt && pt<=PosMax){//-Position in selected domain.
      vnum1[c1]++;
      vnum2[c2]++;
      vnum3[c3]++;
    }
  }
  //printf("%s\n",fun::VarStr("vnum1",npt1,vnum1).c_str());
  //printf("%s\n",fun::VarStr("vnum2",npt2,vnum2).c_str());
  //printf("%s\n",fun::VarStr("vnum3",npt3,vnum3).c_str());
  //-Compute range of selected points in all directions.  
  auto frange=[](unsigned& cini,unsigned& cfin,unsigned snum,const unsigned* vnum){
    unsigned np=0;
    for(unsigned c=0;c<snum;c++){
      if(vnum[c]){
        if(!np)cini=c;
        np++;
      }
      else if(np)c=snum;
    }
    cfin=cini+np;
  };
  c1ini=c2ini=c3ini=c1fin=c2fin=c3fin=0;
  frange(c1ini,c1fin,npt1,vnum1);
  frange(c2ini,c2fin,npt2,vnum2);
  frange(c3ini,c3fin,npt3,vnum3);
  //printf("1. cini:%u-%u  np:%u /%u \n",c1ini,c1fin,c1fin-c1ini,npt1);
  //printf("2. cini:%u-%u  np:%u /%u \n",c2ini,c2fin,c2fin-c2ini,npt2);
  //printf("3. cini:%u-%u  np:%u /%u \n",c3ini,c3fin,c3fin-c3ini,npt3);
  //-Config new grid of points according to selected ranges.  
  StMeshPts m2=mp;
  m2.ptref=mp.ptref+(mp.vdp1*c1ini)+(mp.vdp2*c2ini)+(mp.vdp3*c3ini);
  m2.npt1=c1fin-c1ini;
  m2.npt2=c2fin-c2ini;
  m2.npt3=c3fin-c3ini;
  m2.npt=m2.npt1*m2.npt2*m2.npt3;
  //JMeshTDatasSave::SaveVtkScheme("_DG_FilterPos_Scheme.vtk",g2);
  return(m2);
}

//==============================================================================
/// Applies OnlyPos filter.
//==============================================================================
void JMeshTDatasLoad::RunFilterPos(std::vector<JMeshData*>& vdata)const{
  const unsigned ndata=unsigned(vdata.size());
  if(ndata){
    const StMeshPts m1=vdata[0]->GetMeshPt();
    //-Compute new mesh definition.
    unsigned c1ini,c1fin,c2ini,c2fin,c3ini,c3fin;
    const StMeshPts m2=ComputeFilterPosRanges(m1,c1ini,c1fin,c2ini,c2fin,c3ini,c3fin);
    //-Filter positions.
    if(m2.npt1!=m1.npt1 || m2.npt2!=m1.npt2 || m2.npt3!=m1.npt3){
      const unsigned size  =m2.npt;
      const unsigned size12=m2.npt1*m2.npt2;
      const unsigned sdata=unsigned(vdata.size());
      for(unsigned cd=0;cd<sdata;cd++){
        JMeshData* mdat1=vdata[cd];
        JMeshData* mdat2=new JMeshData();
        //JMeshData::PrintMeshPts(m2);
        mdat2->ConfigMesh(m2,mdat1->GetTimeStep());
        const unsigned na=mdat1->GetVarCount();
        for(unsigned ca=0;ca<na;ca++){
          const JDataArrays::StDataArray& ar1=mdat1->GetArrays()->GetArrayCte(ca);
          const bool dat12=(ar1.tag==JMeshData::TAGDATA12);
          const unsigned ss=(dat12? size12: size);
          const unsigned n1=m1.npt1;
          const unsigned n12=(dat12? 0: m1.npt1*m1.npt2);
          const unsigned np3=(dat12? 1: m2.npt3);
          switch(ar1.type){
            case TypeFloat:{
              //printf("=> float m2.npt:(%u,%u,%u)   n1:%u n12:%u\n",m2.npt1,m2.npt2,np3,n1,n12);
              const float* pdat1=(float*)ar1.ptr;
              float* pdat2=mdat2->CreateArrayPtrFloat(ar1.fullname,ar1.tag,ss);
              for(unsigned c3=0;c3<np3;c3++)for(unsigned c2=0;c2<m2.npt2;c2++)for(unsigned c1=0;c1<m2.npt1;c1++){
                *(pdat2++)=pdat1[(c1+c1ini) + (c2+c2ini)*n1 + (c3+c3ini)*n12];
              }
            }break;
            case TypeFloat3:{
              const tfloat3* pdat1=(tfloat3*)ar1.ptr;
              tfloat3* pdat2=mdat2->CreateArrayPtrFloat3(ar1.fullname,ar1.tag,ss);
              for(unsigned c3=0;c3<np3;c3++)for(unsigned c2=0;c2<m2.npt2;c2++)for(unsigned c1=0;c1<m2.npt1;c1++){
                *(pdat2++)=pdat1[(c1+c1ini) + (c2+c2ini)*n1 + (c3+c3ini)*n12];
              }
            }break;
            case TypeDouble:{
              const double* pdat1=(double*)ar1.ptr;
              double* pdat2=mdat2->CreateArrayPtrDouble(ar1.fullname,ar1.tag,ss);
              for(unsigned c3=0;c3<np3;c3++)for(unsigned c2=0;c2<m2.npt2;c2++)for(unsigned c1=0;c1<m2.npt1;c1++){
                *(pdat2++)=pdat1[(c1+c1ini) + (c2+c2ini)*n1 + (c3+c3ini)*n12];
              }
            }break;
            case TypeDouble3:{
              const tdouble3* pdat1=(tdouble3*)ar1.ptr;
              tdouble3* pdat2=mdat2->CreateArrayPtrDouble3(ar1.fullname,ar1.tag,ss);
              for(unsigned c3=0;c3<np3;c3++)for(unsigned c2=0;c2<m2.npt2;c2++)for(unsigned c1=0;c1<m2.npt1;c1++){
                *(pdat2++)=pdat1[(c1+c1ini) + (c2+c2ini)*n1 + (c3+c3ini)*n12];
              }
            }break;
            default: Run_ExceptioonFile(fun::PrintStr("Type of data \'%s\' is invalid.",ar1.keyname.c_str()),FileIn);
          }
        }
        vdata[cd]=mdat2;
        delete mdat1;
      }    
    }
  }
}

//==============================================================================
/// Applies time and position modifications.
//==============================================================================
void JMeshTDatasLoad::SetTimePos(std::vector<JMeshData*>& vdata)const{
  //printf("SetTimePos> t:%f  ps:(%f,%f,%f)\n",SetTime,SetPos.x,SetPos.y,SetPos.z);
  const bool stim=(SetTime!=0);
  const bool spos=(SetPos!=TDouble3(0));
  if(stim || spos){
    const unsigned ndata=unsigned(vdata.size());
    for(unsigned c=0;c<ndata;c++){
      if(stim)vdata[c]->SetTimeStep(vdata[c]->GetTimeStep()+SetTime);
      if(spos)vdata[c]->SetPtRef(vdata[c]->GetPtRef()+SetPos);
    }
  }
}

//==============================================================================
/// Loads data from a file (.mbi4 or .csv).
//==============================================================================
void JMeshTDatasLoad::LoadFile(const std::string file
  ,std::vector<JMeshData*>& vdata,std::string varlist)
{
  Reset();
  FileIn=file;
  if(!vdata.empty())Run_ExceptioonFile("Vector for data is not empty.",FileIn);
  //printf("=> filein_:[%s]\n",FileIn.c_str());
  const string fext=fun::StrLower(fun::GetExtension(FileIn));
  if(fext=="mbi4")LoadFileBin(vdata,varlist);
  else if(fext=="csv")LoadFileCsv(vdata);
  else Run_ExceptioonFile("File extension is invalid. Only CSV and MBI4 files are supported.",FileIn);
  //-Applies OnlyPos filter.
  if(FilterPos)RunFilterPos(vdata);
  //-Applies time and position modifications.
  SetTimePos(vdata);
}


}
