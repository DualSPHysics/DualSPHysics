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

/// \file JDsPips.cpp \brief Implements the class \ref JDsPips.

#include "JDsPips.h"
#include "JCellSearch_inline.h"
#include "JLog2.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "JSaveCsv2.h"
#ifdef _WITHGPU
  #include "JDsPips_ker.h"
#endif

#include <climits>
#include <cfloat>
#include <cstring>

using namespace std;

//##############################################################################
//# JDsPips
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsPips::JDsPips(bool cpu,unsigned stepsnum,bool svdata,unsigned ntimes)
  :Log(AppInfo.LogPtr()),Cpu(cpu),StepsNum(stepsnum),SvData(svdata),Ntimes(ntimes)
{
  ClassName="JDsPips";
  NextNstep=0;
  NewData=0;
  SizeResultAux=0;
  ResultAux=NULL;
  //Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsPips::~JDsPips(){
  DestructorActive=true;
  //Reset();
  SizeResultAux=0;
  delete[] ResultAux; ResultAux=NULL;
}

//==============================================================================
/// Returns the allocated memory.
//==============================================================================
long long JDsPips::GetAllocMemory()const{
  return(sizeof(StPipsInfo)*int(Data.capacity()));
}

//==============================================================================
/// Returns GigaPIs of requested interval data (Ntimes is applied).
//==============================================================================
double JDsPips::GetGPIs(unsigned cdata)const{
  double gpis=0;
  if(cdata>0 && cdata<unsigned(Data.size())){
    const ullong pi0=Data[cdata-1].pirb + Data[cdata-1].pirf;
    const ullong pi1=Data[cdata  ].pirb + Data[cdata  ].pirf;
    const unsigned nsteps=Data[cdata].nstep-Data[cdata-1].nstep+1;
    gpis=(((double(pi0)/2 + double(pi1)/2)*nsteps - double(pi1))/1e9)*Ntimes;
  }
  return(gpis);
}

//==============================================================================
/// Returns fluid or bound GigaPIs of requested interval data (Ntimes is applied).
//==============================================================================
double JDsPips::GetGPIsType(unsigned cdata,bool fluid)const{
  double gpis=0;
  if(cdata>0 && cdata<unsigned(Data.size())){
    const ullong pi0=(fluid? Data[cdata-1].pirf: Data[cdata-1].pirb);
    const ullong pi1=(fluid? Data[cdata  ].pirf: Data[cdata  ].pirb);
    const unsigned nsteps=Data[cdata].nstep-Data[cdata-1].nstep+1;
    gpis=(((double(pi0)/2 + double(pi1)/2)*nsteps - double(pi1))/1e9)*Ntimes;
  }
  return(gpis);
}

//==============================================================================
/// Returns fluid or bound check GigaPIs of requested interval data (Ntimes is applied).
//==============================================================================
double JDsPips::GetCheckGPIsType(unsigned cdata,bool fluid)const{
  double gpis=0;
  if(cdata>0 && cdata<unsigned(Data.size())){
    const ullong pi0=(fluid? Data[cdata-1].picf: Data[cdata-1].picb);
    const ullong pi1=(fluid? Data[cdata  ].picf: Data[cdata  ].picb);
    const unsigned nsteps=Data[cdata].nstep-Data[cdata-1].nstep+1;
    gpis=(((double(pi0)/2 + double(pi1)/2)*nsteps - double(pi1))/1e9)*Ntimes;
  }
  return(gpis);
}

//==============================================================================
/// Saves stored data in CSV file "PIPS.csv"
//==============================================================================
void JDsPips::SaveData(){
  const unsigned ndata=unsigned(Data.size());
  if(ndata>NewData){
    const string file=AppInfo.GetDirOut()+"PIPS.csv";
    jcsv::JSaveCsv2 scsv(file,true,AppInfo.GetCsvSepComa());
    if(scsv.GetAppendMode())Log->AddFileInfo(file,"Saves CSV with PIPS data in detail.");
    scsv.SetHead();
    scsv << "Time [s];RunTime [s];Nstep";
    scsv << "RealFluid [PIs];RealBound [PIs];ChkFluid [PIs];ChkBound [PIs]";
    scsv << "GPIPS;SumRealFluid [GPIs];SumRealBound [GPIs];SumChkFluid [GPIs];SumChkBound [GPIs]" << jcsv::Endl();
    scsv.SetData();
    for(unsigned c=NewData;c<ndata;c++){
      const StPipsInfo& v=Data[c];
      const double rf=GetGPIsType(c,true );
      const double rb=GetGPIsType(c,false);
      const double cf=GetCheckGPIsType(c,true );
      const double cb=GetCheckGPIsType(c,false);
      const double gpis=GetGPIs(c);
      const double t=(c>0? Data[c].tsim - Data[c-1].tsim: 1);
      const double gpips=gpis/t;
      scsv << v.tstep << v.tsim << v.nstep << v.pirf << v.pirb << v.picf << v.picb;
      scsv << gpips << rf << rb << cf << cb << jcsv::Endl();
    }
    NewData=ndata;
    scsv.SaveData(true);
  }
}

//==============================================================================
/// Returns total PIs of fluid and boundary.
//==============================================================================
tdouble2 JDsPips::GetTotalPIs()const{
  const unsigned ndata=unsigned(Data.size());
  double totgpisf=0,totgpisb=0;
  for(unsigned c=1;c<ndata;c++){
    totgpisf+=GetGPIsType(c,true );
    totgpisb+=GetGPIsType(c,false);
  }
  return(TDouble2(totgpisf,totgpisb));
}

//==============================================================================
/// Returns GigaPIPS of simulation.
//==============================================================================
double JDsPips::GetGPIPS(double tsim)const{
  const tdouble2 totgpis=GetTotalPIs();
  return(tsim>0? (totgpis.x+totgpis.y)/tsim: 0);
}

//==============================================================================
/// Returns total PIs information as string.
/// Returns: 0.123456 HPIs (0.1234e12 + 0.1234e12)
//==============================================================================
std::string JDsPips::GetTotalPIsInfo()const{
  const tdouble2 totgpis=GetTotalPIs();
  double v=(totgpis.x+totgpis.y);
  char unit='G';
  if(v<0.1 ){ v*=1.e3; unit='M'; }
  if(v<0.1 ){ v*=1.e3; unit='k'; }
  if(v>1.e3){ v/=1.e3; unit='T'; }
  if(v>1.e3){ v/=1.e3; unit='P'; }
  if(v>1.e3){ v/=1.e3; unit='E'; }
  if(v>1.e3){ v/=1.e3; unit='Z'; }
  if(v>1.e3){ v/=1.e3; unit='Y'; }
  string vn=fun::DoublexStr(v,fun::NaturalFmt(v,12,true).c_str());
  return(fun::PrintStr("%s %cPIs (%.4e + %.4e)",vn.c_str(),unit,totgpis.x*1e9,totgpis.y*1e9));
}

//==============================================================================
/// Compute number of particle interactions on CPU.
//==============================================================================
void JDsPips::ComputeCpu(unsigned nstep,double tstep,double tsim
  ,const StCteSph &csp,int ompthreads
  ,unsigned np,unsigned npb,unsigned npbok
  ,const StDivDataCpu &dvd,const unsigned *dcell,const tdouble3 *pos)
{
  //-Compute bound & fluid PIs.
  ullong npith[OMP_MAXTHREADS*OMP_STRIDE];
  memset(npith,0,sizeof(ullong)*OMP_MAXTHREADS*OMP_STRIDE);
  //-Counts bound PIs;
  const int inpbok=int(npbok);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=0;p1<inpbok;p1++){
    unsigned picb=0,pirb=0;
    const tdouble3 posp1=pos[p1];
    //-Search for fluid neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(dcell[p1],false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float rr2=nsearch::Distance2(posp1,pos[p2]);
        if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO)pirb++;
        picb++;
      }
    }
    //-Sum results.
    const int cth=omp_get_thread_num()*OMP_STRIDE;
    npith[cth  ]+=picb;
    npith[cth+1]+=pirb;
  }
  //-Counts fluid PIs;
  const int inpb=int(npb);
  const int inp =int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=inpb;p1<inp;p1++){
    unsigned picf=0,pirf=0;
    const tdouble3 posp1=pos[p1];
    //-Search for bound & fluid neighbours in adjacent cells.
    for(byte tpfluid=0;tpfluid<=1;tpfluid++){
      const StNgSearch ngs=nsearch::Init(dcell[p1],!tpfluid,dvd);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float rr2=nsearch::Distance2(posp1,pos[p2]);
          if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO)pirf++;
          picf++;
        }
      }
    }
    //-Sum results.
    const int cth=omp_get_thread_num()*OMP_STRIDE;
    npith[cth+2]+=picf;
    npith[cth+3]+=pirf;
  }
  //-Sum total results.
  StPipsInfo v;
  v.nstep=nstep;
  v.tstep=tstep;
  v.tsim= tsim;
  v.pirf=v.pirb=v.picf=v.picb=0;
  for(int th=0;th<ompthreads;th++){
    int cth=th*OMP_STRIDE;
    v.picb+=npith[cth++];
    v.pirb+=npith[cth++];
    v.picf+=npith[cth++];
    v.pirf+=npith[cth  ];
  }
  //-Stores results.
  Data.push_back(v);
  NextNstep+=StepsNum;
  
  //Log->Printf("%u> PIf: %llu/%llu   PIb: %llu/%llu",v.nstep,v.pirf,v.picf,v.pirb,v.picb);
}

#ifdef _WITHGPU
//==============================================================================
/// Compute number of particle interactions on GPU.
//==============================================================================
void JDsPips::ComputeGpu(unsigned nstep,double tstep,double tsim
  ,unsigned np,unsigned npb,unsigned npbok
  ,const StDivDataGpu &dvd,const unsigned *dcell,const float4 *poscell
  ,unsigned sauxmem,unsigned *auxmem)
{
  //-First calculation with reduction and the result is saved as uint4.
  const unsigned n1=npbok+(np-npb);
  const unsigned n1size=cupips::InteractionNgSize_1st(n1);
  if(n1size*4>sauxmem)Run_Exceptioon("Auxiliary memory is not enough for 1st calculation level.");
  //Log->Printf("\n==> np:%u  npb:%u  npbok:%u  n1:%u  n1size:%u",np,npb,npbok,n1,n1size);
  cupips::InteractionNg_1st(npbok,0,np-npb,npb,dvd,dcell,poscell,(uint4*)auxmem);

  //-Second reduction and the result is saved as uint4.
  const unsigned n2=n1size;
  const unsigned n2size=cupips::InteractionNgSize_2nd(n2);
  if((n1size+n2size)*4>sauxmem)Run_Exceptioon("Auxiliary memory is not enough for 2nd calculation level.");
  unsigned *auxmem2=auxmem+(n1size*4);
  cupips::InteractionNg_2nd(n2,(const uint4*)auxmem,(uint4*)auxmem2);

  //-Last reduction and the result is saved as ullong*4.
  const unsigned n3=n2size;
  const unsigned n3size=cupips::InteractionNgSize_3th(n3);
  if((n1size+n2size)*4+n3size*8>sauxmem)Run_Exceptioon("Auxiliary memory is not enough for 3th calculation level.");
  ullong *auxmem3=(ullong*)(auxmem2+(n2size*4));
  cupips::InteractionNg_3th(n3,(const uint4*)auxmem2,auxmem3);
  //-Allocates memory to final results.
  //Log->Printf("\n==> n3:%u  n3size:%u",n3,n3size);
  if(SizeResultAux<n3size*4){
    delete[] ResultAux; ResultAux=NULL;
    SizeResultAux=n3size*4;
    ResultAux=new ullong[SizeResultAux];
  }
  //-Copy final results from GPU memory.
  if(n3size)cudaMemcpy(ResultAux,auxmem3,sizeof(ullong)*n3size*4,cudaMemcpyDeviceToHost);
  //-Sum total results.
  ullong totnrf=0,totnrb=0,totncf=0,totncb=0;
  for(unsigned c=0,c4=0;c<n3size;c++){
    totnrf+=ResultAux[c4++];
    totnrb+=ResultAux[c4++];
    totncf+=ResultAux[c4++];
    totncb+=ResultAux[c4++];
  }
  StPipsInfo v;
  v.nstep=nstep;
  v.tstep=tstep;
  v.tsim= tsim;
  v.pirf=totnrf;
  v.pirb=totnrb;
  v.picf=totncf;
  v.picb=totncb;
  //-Stores results.
  Data.push_back(v);
  NextNstep+=StepsNum;
  
  //Log->Printf("%u> PIf: %llu/%llu   PIb: %llu/%llu",v.nstep,v.pirf,v.picf,v.pirb,v.picb);
}
#endif
