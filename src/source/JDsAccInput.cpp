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

/// \file JDsAccInput.cpp \brief Implements the class \ref JDsAccInput.

#include "JDsAccInput.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JXml.h"
#include "Functions.h"
#include "JLinearValue.h"
#include "JRangeFilter.h"
#include "JSphMk.h"
#ifdef _WITHGPU
  #include "JDsAccInput_ker.h"
#endif

#include <climits>
#include <cfloat>

using namespace std;

//##############################################################################
//# JDsAccInputMk
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsAccInputMk::JDsAccInputMk(unsigned idx,bool bound,word mktype1,word mktype2
  ,double tini,double tend,bool genabled,tfloat3 acccentre,const JLinearValue &acedata
  ,const JLinearValue &veldata)
  :Log(AppInfo.LogPtr()),Idx(idx),Bound(bound),MkType1(mktype1),MkType2(mktype2)
  ,TimeIni(tini),TimeEnd(tend),GravityEnabled(genabled),AccCoG(acccentre)
{
  ClassName="JDsAccInputMk";
  AceData=NULL;
  VelData=NULL;
  Reset();
  AceData=new JLinearValue(acedata);
  VelData=new JLinearValue(veldata);
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsAccInputMk::~JDsAccInputMk(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsAccInputMk::Reset(){
  CodeSel1=CodeSel2=0;
  delete AceData; AceData=NULL;
  delete VelData; VelData=NULL;
  LastTimestepInput=-1;
  memset(&LastOutput,0,sizeof(StAceInput));
}

//==============================================================================
/// Returns the allocated memory.
//==============================================================================
long long JDsAccInputMk::GetAllocMemory()const{
  return(AceData->GetAllocMemory()+VelData->GetAllocMemory());
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JDsAccInputMk::GetConfig(std::vector<std::string> &lines)const{
  if(MkType1==MkType2)lines.push_back(fun::PrintStr("Input_%u (%s:%u)",Idx,(Bound? "mkbound": "mkfluid"),MkType1));
  else                lines.push_back(fun::PrintStr("Input_%u (%s:%u - %u)",Idx,(Bound? "mkbound": "mkfluid"),MkType1,MkType2));
  if(!AceData->GetFile().empty())lines.push_back(fun::PrintStr("  Data file.....: %s",AceData->GetFile().c_str()));
  lines.push_back(fun::PrintStr("  Time interval.: %g - %s [s]",TimeIni,fun::DoublexStr(TimeEnd).c_str()));
  lines.push_back(fun::PrintStr("  Global gravity: %s",(GravityEnabled? "True": "False")));
  lines.push_back(fun::PrintStr("  Acc center....: (%g,%g,%g) [m]",AccCoG.x,AccCoG.y,AccCoG.z));
}

//=================================================================================================================
/// Returns interpolation variable acceleration values. SL: Added angular and linear velocity and set gravity flag
//=================================================================================================================
const StAceInput& JDsAccInputMk::GetAccValues(double timestep){
  if(LastTimestepInput>=0 && timestep==LastTimestepInput)return(LastOutput);
  LastTimestepInput=timestep;
  //Return values.
  if(TimeIni<=timestep && timestep<=TimeEnd){
    LastOutput.codesel1=CodeSel1;
    LastOutput.codesel2=CodeSel2;
    LastOutput.centre=ToTDouble3(AccCoG);
    LastOutput.setgravity=GravityEnabled;
    AceData->GetValue3d3d(timestep,LastOutput.acclin,LastOutput.accang);
    VelData->GetValue3d3d(timestep,LastOutput.vellin,LastOutput.velang);
  }
  else LastOutput.codesel1=UINT_MAX;
  return(LastOutput);
}

//##############################################################################
//# JDsAccInput
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JDsAccInput::JDsAccInput(const std::string &dirdata,const JXml *sxml
  ,const std::string &place):Log(AppInfo.LogPtr()),DirData(dirdata)
{
  ClassName="JDsAccInput";
  Reset();
  LoadXml(sxml,place);
}

//==============================================================================
/// Destructor.
//==============================================================================
JDsAccInput::~JDsAccInput(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JDsAccInput::Reset(){
  for(unsigned c=0;c<Inputs.size();c++)delete Inputs[c];
  Inputs.clear();
  MemSize=0;
}

//==============================================================================
/// Returns true if mktype value is already configured.
//==============================================================================
bool JDsAccInput::ExistMk(bool bound,word mktype)const{
  bool ret=false;
  for(unsigned c=0;c<Inputs.size() && !ret;c++)
    ret=(Inputs[c]->Bound==bound && Inputs[c]->MkType1<=mktype && mktype<=Inputs[c]->MkType2);
  return(ret);
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JDsAccInput::LoadXml(const JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JDsAccInput::ReadXml(const JXml *sxml,TiXmlElement* lis){
  //-Loads list of inputs.
  TiXmlElement* ele=lis->FirstChildElement("accinput"); 
  while(ele){
    if(sxml->CheckElementActive(ele)){
      //-Check XML configuration.
      if(sxml->ExistsElement(ele,"mkfluid" ))Run_ExceptioonFile("Element <mkfluid> is invalid for current version. Update the XML file.",sxml->ErrGetFileRow(ele));
      if(sxml->ExistsElement(ele,"datafile"))Run_ExceptioonFile("Element <datafile> is invalid for current version. Update the XML file.",sxml->ErrGetFileRow(ele));
      sxml->CheckElementNames(ele,true,"time acccentre globalgravity acctimes acctimesfile");
      if(sxml->ExistsElement(ele,"acctimes") && sxml->ExistsElement(ele,"acctimesfile"))
        Run_ExceptioonFile("Only <acctimes> or <acctimesfile> is valid but not both.",sxml->ErrGetFileRow(ele));
      if(sxml->ExistsAttribute(ele,"mkbound") && sxml->ExistsAttribute(ele,"mkfluid"))
        Run_ExceptioonFile("Only mkbound or mkfluid values is valid but not both.",sxml->ErrGetFileRow(ele));
      //-Load general configuration.
      const double tini=sxml->ReadElementDouble(ele,"time","start",true,0);
      const double tend=sxml->ReadElementDouble(ele,"time","end",true,DBL_MAX);
      const tfloat3 acccentre=sxml->ReadElementFloat3(ele,"acccentre");
      const bool genabled=sxml->ReadElementBool(ele,"globalgravity","value");
      const bool bound=sxml->ExistsAttribute(ele,"mkbound");
      word mktype1=USHRT_MAX;
      word mktype2=USHRT_MAX;
      string strmktype=sxml->GetAttributeStrSimple(ele,(bound? "mkbound": "mkfluid"));
      if(strmktype.empty() || fun::StrIsIntegerNumber(strmktype) || strmktype[0]=='#'){
        mktype1=mktype2=sxml->GetAttributeWord(ele,(bound? "mkbound": "mkfluid"));
      }
      //-Loads acceleration values.
      JLinearValue acedata(6,false,true);
      if(sxml->ExistsElement(ele,"acctimes"))acedata.ReadXmlValues(sxml,ele,"acctimes","timevalue","time:linx:liny:linz:angx:angy:angz");
      else acedata.LoadFile(DirData+sxml->ReadElementStr(ele,"acctimesfile","value"));
      //-Computes velocity values.
      JLinearValue veldata(6,false);
      ComputeVelocity(acedata,veldata);
      //-Create input configurations.
      if(mktype1!=USHRT_MAX){
        if(ExistMk(bound,mktype1))Run_ExceptioonFile(fun::PrintStr("An input already exists for the same %s=%u.",(bound? "mkbound": "mkfluid"),mktype1),sxml->ErrGetFileRow(ele));
        JDsAccInputMk *input=new JDsAccInputMk(GetCount(),bound,mktype1,mktype2,tini,tend,genabled,acccentre,acedata,veldata);
        Inputs.push_back(input);
      }
      else{//-Check range of mkvalues.
        std::vector<unsigned> vmk;
        JRangeFilter rg(strmktype);
        rg.GetValues(vmk);
        const unsigned nmk=unsigned(vmk.size());
        for(unsigned c=0;c<nmk;c++)if(ExistMk(bound,word(vmk[c])))Run_ExceptioonFile(fun::PrintStr("An input already exists for the same %s=%u.",(bound? "mkbound": "mkfluid"),vmk[c]),sxml->ErrGetFileRow(ele));
        if(!nmk)Run_ExceptioonFile((bound? "The mkbound is invalid.": "The mkfluid is invalid."),sxml->ErrGetFileRow(ele));
        mktype1=mktype2=word(vmk[0]);
        for(unsigned c=1;c<nmk;c++){
          word v=word(vmk[c]);
          if(mktype2+1==v)mktype2=v;
          else{
            JDsAccInputMk *input=new JDsAccInputMk(GetCount(),bound,mktype1,mktype2,tini,tend,genabled,acccentre,acedata,veldata);
            Inputs.push_back(input);
            mktype1=mktype2=v;
          }
        }
        JDsAccInputMk *input=new JDsAccInputMk(GetCount(),bound,mktype1,mktype2,tini,tend,genabled,acccentre,acedata,veldata);
        Inputs.push_back(input);
      }
    }
    ele=ele->NextSiblingElement("accinput");
  }
  //-Calculate allocated memory.
  for(unsigned c=0;c<Inputs.size();c++){
    MemSize+=Inputs[c]->GetAllocMemory();
  }
}

//==============================================================================
/// Compute velocity starting from acceleration.
//==============================================================================
void JDsAccInput::ComputeVelocity(const JLinearValue &acedata,JLinearValue &veldata)const{
  const unsigned nt=acedata.GetCount();
  double atime0=0;
  tdouble3 vellin0=TDouble3(0);
  tdouble3 velang0=TDouble3(0);
  for(unsigned ct=0;ct<nt;ct++){
    const double atime=acedata.GetTimeByIdx(ct);        //Time.
    const tdouble3 acclin=acedata.GetValue3ByIdx(ct,0); //Acc Linear.
    const tdouble3 accang=acedata.GetValue3ByIdx(ct,1); //Acc Angular.
    //Log->Printf("t:%g  lin:(%g  %g  %g)  ang:(%g  %g  %g)",atime,acclin.x,acclin.y,acclin.z,accang.x,accang.y,accang.z);
    //SL: Calculate linear velocity vector based on acceleration and time data loaded
    tdouble3 currvellin=TDouble3(0); //SL: New linear velocity variable
    if(ct>0){ //SL: Angular velocity is always zero at time zero
      const double dt=atime-atime0;
      currvellin.x=vellin0.x+(acclin.x*dt);
      currvellin.y=vellin0.y+(acclin.y*dt);
      currvellin.z=vellin0.z+(acclin.z*dt);
    }
    //SL: Calculate angular velocity vector based on acceleration and time data loaded
    tdouble3 currvelang=TDouble3(0); //SL: New angular velocity variable
    if(ct>0){ //SL: Angular velocity is always zero at time zero
      const double dt=atime-atime0;
      currvelang.x=velang0.x+(accang.x*dt);
      currvelang.y=velang0.y+(accang.y*dt);
      currvelang.z=velang0.z+(accang.z*dt);
    }
    //SL: Save the calculated linar and angular velocity.
    veldata.AddTimeValue(atime,currvellin.x,currvellin.y,currvellin.z,currvelang.x,currvelang.y,currvelang.z);
    atime0=atime;
    vellin0=currvellin;
    velang0=currvelang;
  }
  //Check that at least 2 values were given or interpolation will be impossible.
  if(nt<2)Run_Exceptioon("Cannot be less than two registers in variable acceleration data.");
  //Check that the final value for time is not smaller than the final simulation time.
  //if(double(AccTime[AccCount-1])<tmax)Run_ExceptioonFile(fun::PrintStr("Final time (%g) is less than total simulation time in variable acceleration file.",AccTime[AccCount-1]),file);
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JDsAccInput::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  for(unsigned c=0;c<GetCount();c++){
    const JDsAccInputMk* ip=Inputs[c];
    std::vector<std::string> lines;
    ip->GetConfig(lines);
    for(unsigned i=0;i<unsigned(lines.size());i++)Log->Print(lines[i]);
  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Checks config according to Mk of loaded particles.
//==============================================================================
void JDsAccInput::Init(const JSphMk *mkinfo){
  for(unsigned c=0;c<GetCount();c++){
    JDsAccInputMk* ip=Inputs[c];
    //-Check configuration.
    if(ip->Bound){//-For floating bodies.
      for(word mktype=ip->MkType1;mktype<=ip->MkType2;mktype++){
        unsigned cmk=mkinfo->GetMkBlockByMkBound(mktype);
        if(cmk>=mkinfo->Size())Run_Exceptioon(fun::PrintStr("The MkBound %u with imposed acceleration is not present in the simulation.",mktype));
        if(mkinfo->Mkblock(cmk)->Type!=TpPartFloating)Run_Exceptioon(fun::PrintStr("The MkBound %u with imposed acceleration is not floating body.",mktype));
      }
    }
    else{//-For fluid.
      for(word mktype=ip->MkType1;mktype<=ip->MkType2;mktype++){
        unsigned cmk=mkinfo->GetMkBlockByMkFluid(mktype);
        if(cmk>=mkinfo->Size())Run_Exceptioon(fun::PrintStr("The MkFluid %u with imposed acceleration is not present in the simulation.",mktype));
      }
    }
    //-Config code TypeAndValue.
    typecode cod1,cod2;
    if(ip->Bound){//-For floating bodies.
      const unsigned cmk1=mkinfo->GetMkBlockByMkBound(ip->MkType1);
      const unsigned cmk2=mkinfo->GetMkBlockByMkBound(ip->MkType2);
      cod1=mkinfo->Mkblock(cmk1)->Code;
      cod2=mkinfo->Mkblock(cmk2)->Code;
      //Log->Printf("==> MkBound:%u-%u  code:%u-%u",ip->MkType1,ip->MkType2,cod1,cod2);
    }
    else{//-For fluid.
      const unsigned cmk1=mkinfo->GetMkBlockByMkFluid(ip->MkType1);
      const unsigned cmk2=mkinfo->GetMkBlockByMkFluid(ip->MkType2);
      cod1=mkinfo->Mkblock(cmk1)->Code;
      cod2=mkinfo->Mkblock(cmk2)->Code;
      //Log->Printf("==> MkFluid:%u-%u  code:%u-%u",ip->MkType1,ip->MkType2,cod1,cod2);
    }
    ip->ConfigCodeSel(cod1,cod2);
  }
}

//=====================================================================================================================================================
/// Returns interpolation variable acceleration values. SL: Corrected spelling mistake in exception and added angular velocity and global gravity flag
//=====================================================================================================================================================
const StAceInput& JDsAccInput::GetAccValues(unsigned cinput,double timestep){
  if(cinput>=GetCount())Run_Exceptioon("The number of input data for variable acceleration is invalid.");
  return(Inputs[cinput]->GetAccValues(timestep));
}

//==============================================================================
/// Adds variable acceleration from input configurations.
//==============================================================================
void JDsAccInput::RunCpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
  ,const typecode *code,const tdouble3 *pos,const tfloat4 *velrhop,tfloat3 *ace)
{
  for(unsigned c=0;c<GetCount();c++){
    const StAceInput v=GetAccValues(c,timestep);
    if(v.codesel1!=UINT_MAX){
      const bool withaccang=(v.accang.x!=0 || v.accang.y!=0 || v.accang.z!=0);
      //const typecode codesel=typecode(v.mkfluid);
      const typecode codesel1=typecode(v.codesel1);
      const typecode codesel2=typecode(v.codesel2);
      const int ppini=int(pini),ppfin=pini+int(n);
      #ifdef OMP_USE
        #pragma omp parallel for schedule (static)
      #endif
      for(int p=ppini;p<ppfin;p++){//-Iterates through the fluid particles.
        //-Checks if the current particle is part of the particle set by its MK.
        const typecode tav=CODE_GetTypeAndValue(code[p]);
        if(codesel1<=tav && tav<=codesel2){
          tdouble3 acc=ToTDouble3(ace[p]);
          acc=acc+v.acclin;                             //-Adds linear acceleration.
          if(!v.setgravity)acc=acc-ToTDouble3(gravity); //-Subtract global gravity from the acceleration if it is set in the input file
          if(withaccang){                               //-Adds angular acceleration.
            const tdouble3 dc=pos[p]-v.centre;
            const tdouble3 vel=TDouble3(velrhop[p].x,velrhop[p].y,velrhop[p].z);//-Get the current particle's velocity

            //-Calculate angular acceleration ((Dw/Dt) x (r_i - r)) + (w x (w x (r_i - r))) + (2w x (v_i - v))
            //(Dw/Dt) x (r_i - r) (term1)
            acc.x+=(v.accang.y*dc.z)-(v.accang.z*dc.y);
            acc.y+=(v.accang.z*dc.x)-(v.accang.x*dc.z);
            acc.z+=(v.accang.x*dc.y)-(v.accang.y*dc.x);

            //-Centripetal acceleration (term2)
            //-First find w x (r_i - r))
            const double innerx=(v.velang.y*dc.z)-(v.velang.z*dc.y);
            const double innery=(v.velang.z*dc.x)-(v.velang.x*dc.z);
            const double innerz=(v.velang.x*dc.y)-(v.velang.y*dc.x);
            //-Find w x inner.
            acc.x+=(v.velang.y*innerz)-(v.velang.z*innery);
            acc.y+=(v.velang.z*innerx)-(v.velang.x*innerz);
            acc.z+=(v.velang.x*innery)-(v.velang.y*innerx);

            //-Coriolis acceleration 2w x (v_i - v) (term3)
            acc.x+=((2.0*v.velang.y)*vel.z)-((2.0*v.velang.z)*(vel.y-v.vellin.y));
            acc.y+=((2.0*v.velang.z)*vel.x)-((2.0*v.velang.x)*(vel.z-v.vellin.z));
            acc.z+=((2.0*v.velang.x)*vel.y)-((2.0*v.velang.y)*(vel.x-v.vellin.x));
          }
          //-Stores the new acceleration value.
          ace[p]=ToTFloat3(acc);
        }
      }
    }
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Adds variable acceleration from input configurations.
//==============================================================================
void JDsAccInput::RunGpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
    ,const typecode *code,const double2 *posxy,const double *posz,const float4 *velrhop,float3 *ace)
{
  for(unsigned c=0;c<GetCount();c++){
    const StAceInput v=GetAccValues(c,timestep);
    if(v.codesel1!=UINT_MAX){
      const typecode codesel1=typecode(v.codesel1);
      const typecode codesel2=typecode(v.codesel2);
      cuaccin::AddAccInput(n,pini,codesel1,codesel2,v.acclin,v.accang,v.centre
        ,v.velang,v.vellin,v.setgravity,gravity,code,posxy,posz,velrhop,ace,NULL);
    }
  }
}
#endif

