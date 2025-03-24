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

/// \file JDsGaugeSystem.cpp \brief Implements the class \ref JGaugeSystem.

#include "JDsGaugeSystem.h"
#include "FunSphKernel.h"
#include "JLog2.h"
#include "JXml.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JSphMk.h"
#include "JDataArrays.h"
#include "JSpVtkData.h"
#include <cfloat>
#include <climits>
#include <algorithm>
#ifdef _WITHGPU
 #include "FunctionsCuda.h"
#endif

using namespace std;

//##############################################################################
//# JGaugeSystem
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JGaugeSystem::JGaugeSystem(int gpucount):Log(AppInfo.LogPtr())
  ,Cpu(gpucount==0),GpuCount(gpucount)
{
  ClassName="JGaugeSystem";
  if(unsigned(GpuCount)>MAXGPUS)Run_Exceptioon("Number of GPUs is invalid.");
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JGaugeSystem::~JGaugeSystem(){
  DestructorActive=true;
  #ifdef _WITHGPU
    if(GpuCount==1)FreeMemoryGpu(0);
  #endif
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JGaugeSystem::Reset(){
  ConfiguredCtes=false;
  CSP=CteSphNull();
  TimeMax=TimePart=0;
  Scell=0;
  ScellDiv=0;
  MapPosMin=TDouble3(0);
  DomPosMin=DomPosMax=TDouble3(0);
  ResetCfgDefault();
  //-Free data CPU.
  DataCpu=JGaugeItem::StrDataCpu();
  //-Free data GPU.
  #ifdef _WITHGPU
    if(GpuCount==1){
      FreeMemoryGpu(0);
      DataGpu[0]=JGaugeItem::StrDataGpu();
    }
  #endif
  //-Removes gauges.
  for(unsigned cg=0;cg<GetCount();cg++)delete Gauges[cg];
  Gauges.clear();
}

#ifdef _WITHGPU
//==============================================================================
/// Check if GPU memory was allocated for some gauge and gpu.
//==============================================================================
bool JGaugeSystem::AllocatedMemoryGpu()const{
  bool allocated=false;
  for(unsigned cg=0;cg<GetCount();cg++)for(int id=0;id<GpuCount;id++){
    allocated=(allocated || Gauges[cg]->AllocatedGpuMemory(id));
  }
  return(allocated);
}
//==============================================================================
/// Free GPU memory on gauges.
//==============================================================================
void JGaugeSystem::FreeMemoryGpu(int id){
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  //-Free gauges auxiliary memory on GPU.
  for(unsigned cg=0;cg<GetCount();cg++){
    Gauges[cg]->FreeGpuMemory(id);
  }
}
#endif

//==============================================================================
/// Initialisation of CfgDefault.
//==============================================================================
void JGaugeSystem::ResetCfgDefault(){
  CfgDefault.savevtkpart=false;
  CfgDefault.computedt=TimePart;
  CfgDefault.computestart=0;
  CfgDefault.computeend=TimeMax;
  CfgDefault.output=false;
  CfgDefault.outputdt=TimePart;
  CfgDefault.outputstart=0;
  CfgDefault.outputend=TimeMax;
}

//==============================================================================
/// Configures constants.
//==============================================================================
void JGaugeSystem::ConfigCtes(const StCteSph& csp,double timemax,double timepart
  ,float scell,int scelldiv,tdouble3 mapposmin,tdouble3 domposmin
  ,tdouble3 domposmax)
{
  CSP=csp;
  //-Wendland kernel is used when Cubic is selected.
  if(CSP.tkernel==KERNEL_Cubic){
    Log->PrintfWarning("The kernel Cubic Spline is not available in GaugeSystem, so kernel Wendland is used.");
    if(!CSP.kwend.awen || !CSP.kwend.bwenh)Run_Exceptioon("Constants of kernel Wendland are not defined.");
  }
  TimeMax=timemax;
  TimePart=timepart;
  Scell=scell;
  ScellDiv=scelldiv;
  MapPosMin=mapposmin;
  DomPosMin=domposmin;
  DomPosMax=domposmax;
  ConfiguredCtes=true;
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JGaugeSystem::LoadXml(const JXml* sxml,const std::string& place
  ,const JSphMk* mkinfo)
{
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement(),mkinfo);
}

//==============================================================================
/// Loads points of line point1-point2.
//==============================================================================
void JGaugeSystem::LoadLinePoints(double coefdp,const tdouble3& point1
  ,const tdouble3& point2,std::vector<tdouble3>& points,const std::string& ref)const
{
  const double dis=fgeo::PointsDist(point1,point2);
  const double dp=CSP.dp*coefdp;
  if(dis<dp)Run_ExceptioonFile("The line length is less than the indicated dp.",ref);
  unsigned count=unsigned(dis/dp);
  if(dis-(dp*count)>=dp*0.1)count++;
  count++;
  if(count<2)count++;
  LoadLinePoints(count,point1,point2,points,ref);
}

//==============================================================================
/// Loads points of line point1-point2.
//==============================================================================
void JGaugeSystem::LoadLinePoints(unsigned count,const tdouble3& point1
  ,const tdouble3& point2,std::vector<tdouble3>& points,const std::string& ref)const
{
  const double dis=fgeo::PointsDist(point1,point2);
  const double dp=dis/(count-1);
  const tdouble3 dir=fgeo::VecUnitary(point2-point1);
  for(unsigned c=0;c<count;c++){
    tdouble3 pt=point1+(dir*(dp*c));
    points.push_back(pt);
  }
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JGaugeSystem::LoadPoints(JXml* sxml,TiXmlElement* lis
  ,std::vector<tdouble3>& points)const
{
  //-Loads points.
  TiXmlElement* ele=lis->FirstChildElement("point"); 
  while(ele){
    points.push_back(sxml->GetAttributeDouble3(ele));
    ele=ele->NextSiblingElement("point");
  }
  //-Loads lines.
  ele=lis->FirstChildElement("line"); 
  while(ele){
    double coefdp=sxml->GetAttributeDouble(ele,"coefdp",true,DBL_MAX);
    unsigned count=sxml->GetAttributeUint(ele,"count",true,0);
    const tdouble3 pt1=sxml->ReadElementDouble3(ele,"point1");
    const tdouble3 pt2=sxml->ReadElementDouble3(ele,"point2");
    if(count!=0 && coefdp!=DBL_MAX)Run_ExceptioonFile("Only \'coefdp\' or \'count\' definition are valid, but not both.",sxml->ErrGetFileRow(ele));

    if(!count)LoadLinePoints(coefdp,pt1,pt2,points,sxml->ErrGetFileRow(ele));
    else LoadLinePoints(count,pt1,pt2,points,sxml->ErrGetFileRow(ele));
    ele=ele->NextSiblingElement("line");
  }

}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
JGaugeItem::StDefault JGaugeSystem::ReadXmlCommon(const JXml* sxml
  ,TiXmlElement* ele)const
{
  JGaugeItem::StDefault cfg=CfgDefault;
  if(ele){
    cfg.savevtkpart =sxml->ReadElementBool  (ele,"savevtkpart","value",true,CfgDefault.savevtkpart);
    cfg.computedt   =sxml->ReadElementDouble(ele,"computedt"  ,"value",true,CfgDefault.computedt);
    cfg.computestart=sxml->ReadElementDouble(ele,"computetime","start",true,CfgDefault.computestart);
    cfg.computeend  =sxml->ReadElementDouble(ele,"computetime","end"  ,true,CfgDefault.computeend);
    cfg.output      =sxml->ReadElementBool  (ele,"output"     ,"value",true,CfgDefault.output);
    cfg.outputdt    =sxml->ReadElementDouble(ele,"outputdt"   ,"value",true,CfgDefault.outputdt);
    cfg.outputstart =sxml->ReadElementDouble(ele,"outputtime" ,"start",true,CfgDefault.outputstart);
    cfg.outputend   =sxml->ReadElementDouble(ele,"outputtime" ,"end"  ,true,CfgDefault.outputend);
  }
  return(cfg);
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JGaugeSystem::ReadXml(const JXml* sxml,TiXmlElement* lis,const JSphMk* mkinfo){
  if(!ConfiguredCtes)Run_Exceptioon("The object is not yet configured.");
  //-Loads default configuration.
  ResetCfgDefault();
  CfgDefault=ReadXmlCommon(sxml,lis->FirstChildElement("default"));
  //-Loads gauge definitions.
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    string cmd=ele->Value();
    if(cmd.length() && cmd[0]!='_' && cmd!="default"){
      if(sxml->CheckElementActive(ele)){
        const string name=sxml->GetAttributeStr(ele,"name");
        if(GetGaugeIdx(name)!=UINT_MAX)Run_ExceptioonFile(fun::PrintStr("The name \'%s\' already exists.",name.c_str()),sxml->ErrGetFileRow(ele));
        const JGaugeItem::StDefault cfg=ReadXmlCommon(sxml,ele);
        //-Loads points
        //std::vector<tdouble3> points;
        //LoadPoints(sxml,ele,points);
        JGaugeItem* gau=NULL;
        if(cmd=="velocity"){
          const tdouble3 point=sxml->ReadElementDouble3(ele,"point");
          gau=AddGaugeVel(name,cfg.computestart,cfg.computeend,cfg.computedt
            ,true,point);
        }
        else if(cmd=="swl"){
          //-Reads masslimit.
          float masslimit=0;
          switch(sxml->CheckElementAttributes(ele,"masslimit","value coef",true,true)){
            case 1:  masslimit=sxml->ReadElementFloat(ele,"masslimit","value");           break;
            case 2:  masslimit=CSP.massfluid*sxml->ReadElementFloat(ele,"masslimit","coef");  break;
            case 0:  masslimit=CSP.massfluid*(CSP.simulate2d? 0.4f: 0.5f);                        break;
          }
          if(masslimit<=0)Run_ExceptioonFile(fun::PrintStr("The masslimit (%f) is invalid.",masslimit),sxml->ErrGetFileRow(ele));
          //-Reads pointdp.
          double pointdp=0;
          switch(sxml->CheckElementAttributes(ele,"pointdp","value coefdp",true,true)){
            case 0:
            case 1:  pointdp=sxml->ReadElementFloat(ele,"pointdp","value");      break;
            case 2:  pointdp=CSP.dp*sxml->ReadElementFloat(ele,"pointdp","coefdp");  break;
          }
          if(pointdp<=0)Run_ExceptioonFile(fun::PrintStr("The pointdp (%f) is invalid.",pointdp),sxml->ErrGetFileRow(ele));
          //-Reads point0 and point2.
          const tdouble3 pt0=sxml->ReadElementDouble3(ele,"point0");
          const tdouble3 pt2=sxml->ReadElementDouble3(ele,"point2");
          gau=AddGaugeSwl(name,cfg.computestart,cfg.computeend,cfg.computedt
            ,true,pt0,pt2,pointdp,masslimit);
        }
        else if(cmd=="maxz"){
          const tdouble3 pt0=sxml->ReadElementDouble3(ele,"point0");
          const float height=sxml->ReadElementFloat(ele,"height","value");
          //-Reads distlimit.
          float distlimit=0;
          switch(sxml->CheckElementAttributes(ele,"distlimit","value coefdp coefh",true,true)){
            case 0:
            case 1:  distlimit=sxml->ReadElementFloat(ele,"distlimit","value");             break;
            case 2:  distlimit=float(CSP.dp*sxml->ReadElementFloat(ele,"distlimit","coefdp"));  break;
            case 3:  distlimit=float(CSP.kernelh*sxml->ReadElementFloat(ele,"distlimit","coefh"));   break;
          }
          if(distlimit<=0)Run_ExceptioonFile(fun::PrintStr("The distlimit (%f) is invalid.",distlimit),sxml->ErrGetFileRow(ele));
          gau=AddGaugeMaxZ(name,cfg.computestart,cfg.computeend,cfg.computedt
            ,true,pt0,height,distlimit);
        }
        else if(cmd=="force"){
          const word mkbound=(word)sxml->ReadElementUnsigned(ele,"target","mkbound");
          gau=AddGaugeForce(name,cfg.computestart,cfg.computeend,cfg.computedt
            ,true,mkinfo,mkbound);
        }
        else if(cmd=="mesh"){ //<vs_meeshdat_ini>
          sxml->CheckElementNames(ele,true,"outputdata outputfmt buffersize dirdat point vec1 vec2 vec3 size1 size2 size3 kclimit kcdummy masslimit savevtkpart computedt computetime output outputdt outputtime");
          jmsh::StMeshBasic mesh={TDouble3(0),TDouble3(0),TDouble3(0),TDouble3(0),0,0,0,0,0,0,TFloat3(0)};
          mesh.ptref=sxml->ReadElementDouble3(ele,"point");
          mesh.vec1=sxml->ReadElementDouble3(ele,"vec1",true);
          mesh.vec2=sxml->ReadElementDouble3(ele,"vec2",true);
          mesh.vec3=sxml->ReadElementDouble3(ele,"vec3",true);
          if(mesh.vec1!=TDouble3(0)){
            mesh.dis1  =sxml->ReadElementDouble(ele,"size1","length");
            mesh.dispt1=sxml->ReadElementDouble(ele,"size1","distpt");
          }
          if(mesh.vec2!=TDouble3(0)){
            mesh.dis2  =sxml->ReadElementDouble(ele,"size2","length");
            mesh.dispt2=sxml->ReadElementDouble(ele,"size2","distpt");
          }
          if(mesh.vec3!=TDouble3(0)){
            mesh.dis3  =sxml->ReadElementDouble(ele,"size3","length");
            mesh.dispt3=sxml->ReadElementDouble(ele,"size3","distpt");
          }
          mesh.dirdat=sxml->ReadElementFloat3(ele,"dirdat",true,TDouble3(1,0,0));
          mesh.dirdat=fgeo::VecUnitary(mesh.dirdat);
          string outdata=sxml->ReadElementStr(ele,"outputdata","value");
          outdata=JGaugeMesh::CorrectDataList(outdata);
          if(outdata.empty())Run_ExceptioonFile("No output data requested.",sxml->ErrGetFileRow(ele,"outputdata"));
          if(outdata=="error")Run_ExceptioonFile("Output data requested is invalid.",sxml->ErrGetFileRow(ele,"outputdata"));
          unsigned tfmt=0;
          {
            string txdata=sxml->ReadElementStr(ele,"outputfmt","value",true,"csv");
            txdata=fun::StrLower(fun::StrWithoutChar(txdata,' '));
            std::vector<std::string> vstr;
            unsigned nstr=fun::VectorSplitStr(",",txdata,vstr);
            for(unsigned c=0;c<nstr;c++)if(!vstr[c].empty()){
                   if(vstr[c]==jmsh::GetTpFormatKey(jmsh::TpFmtBin))tfmt|=jmsh::TpFmtBin;
              else if(vstr[c]==jmsh::GetTpFormatKey(jmsh::TpFmtCsv))tfmt|=jmsh::TpFmtCsv;
              else Run_ExceptioonFile("The output format is invalid.",sxml->ErrGetFileRow(ele,"outputfmt"));
            }
          }
          //-Reads buffersize.
          unsigned buffersize=sxml->ReadElementUnsigned(ele,"buffersize","value",true,0);
          //-Reads kclimit.
          float kclimit=0.5f;
          if(fun::StrLower(sxml->ReadElementStr(ele,"kclimit","value",true))=="none")kclimit=FLT_MAX;
          else kclimit=sxml->ReadElementFloat(ele,"kclimit","value",true,0.5f);
          //-Reads kcdummy.
          float kcdummy=0;
          if(fun::StrLower(sxml->ReadElementStr(ele,"kcdummy","value",true))=="none")kcdummy=FLT_MAX;
          else kcdummy=sxml->ReadElementFloat(ele,"kcdummy","value",true,0);
          //-Reads masslimit.
          float masslimit=0;
          switch(sxml->CheckElementAttributes(ele,"masslimit","value coef",true,true)){
            case 1:  masslimit=sxml->ReadElementFloat(ele,"masslimit","value");           break;
            case 2:  masslimit=CSP.massfluid*sxml->ReadElementFloat(ele,"masslimit","coef");  break;
            case 0:  masslimit=CSP.massfluid*(CSP.simulate2d? 0.4f: 0.5f);                        break;
          }
          if(masslimit<=0)Run_ExceptioonFile(fun::PrintStr("The masslimit (%f) is invalid.",masslimit),sxml->ErrGetFileRow(ele));
          //-Creates gauge object.
          gau=AddGaugeMesh(name,cfg.computestart,cfg.computeend,cfg.computedt
            ,true,mesh,outdata,tfmt,buffersize,kclimit,kcdummy,masslimit);
        }  //<vs_meeshdat_end>
        else Run_ExceptioonFile(fun::PrintStr("Gauge type \'%s\' is invalid.",cmd.c_str()),sxml->ErrGetFileRow(ele));
        gau->SetSaveVtkPart(cfg.savevtkpart);
        //gau->ConfigComputeTiming(cfg.computestart,cfg.computeend,cfg.computedt);
        gau->ConfigOutputTiming (cfg.output,cfg.outputstart,cfg.outputend,cfg.outputdt);
      }
    }
    ele=ele->NextSiblingElement();
  }
  SaveVtkInitPoints();
}

//==============================================================================
/// Creates new gauge-Velocity and returns pointer.
//==============================================================================
JGaugeVelocity* JGaugeSystem::AddGaugeVel(std::string name,double computestart
  ,double computeend,double computedt,bool fixed,const tdouble3& point)
{
  if(GetGaugeIdx(name)!=UINT_MAX)Run_Exceptioon(fun::PrintStr("The name \'%s\' already exists.",name.c_str()));
  //-Creates object.
  JGaugeVelocity* gau=new JGaugeVelocity(GetCount(),name,point,GpuCount);
  gau->Config(CSP,Scell,ScellDiv,MapPosMin,DomPosMin,DomPosMax);
  gau->ConfigDomMCel(fixed);
  gau->ConfigComputeTiming(computestart,computeend,computedt);
  //-Uses common configuration.
  gau->SetSaveVtkPart(CfgDefault.savevtkpart);
  gau->ConfigOutputTiming(CfgDefault.output,CfgDefault.outputstart
    ,CfgDefault.outputend,CfgDefault.outputdt);
  Gauges.push_back(gau);
  return(gau);
}

//==============================================================================
/// Creates new gauge-SWL and returns pointer.
//==============================================================================
JGaugeSwl* JGaugeSystem::AddGaugeSwl(std::string name,double computestart
  ,double computeend,double computedt,bool fixed
  ,tdouble3 point0,tdouble3 point2,double pointdp,float masslimit)
{
  if(GetGaugeIdx(name)!=UINT_MAX)Run_Exceptioon(fun::PrintStr("The name \'%s\' already exists.",name.c_str()));
  if(masslimit<=0)masslimit=CSP.massfluid*(CSP.simulate2d? 0.4f: 0.5f);
  //-Creates object.
  JGaugeSwl* gau=new JGaugeSwl(GetCount(),name,point0,point2,pointdp,masslimit,GpuCount);
  gau->Config(CSP,Scell,ScellDiv,MapPosMin,DomPosMin,DomPosMax);
  gau->ConfigDomMCel(fixed);
  gau->ConfigComputeTiming(computestart,computeend,computedt);
  //-Uses common configuration.
  gau->SetSaveVtkPart(CfgDefault.savevtkpart);
  gau->ConfigOutputTiming(CfgDefault.output,CfgDefault.outputstart
    ,CfgDefault.outputend,CfgDefault.outputdt);
  Gauges.push_back(gau);
  return(gau);
}

//==============================================================================
/// Creates new gauge-MaxZ and returns pointer.
//==============================================================================
JGaugeMaxZ* JGaugeSystem::AddGaugeMaxZ(std::string name,double computestart
  ,double computeend,double computedt,bool fixed
  ,tdouble3 point0,double height,float distlimit)
{
  if(GetGaugeIdx(name)!=UINT_MAX)Run_Exceptioon(fun::PrintStr("The name \'%s\' already exists.",name.c_str()));
  //-Creates object.
  JGaugeMaxZ* gau=new JGaugeMaxZ(GetCount(),name,point0,height,distlimit,GpuCount);
  gau->Config(CSP,Scell,ScellDiv,MapPosMin,DomPosMin,DomPosMax);
  gau->ConfigDomMCel(fixed);
  gau->ConfigComputeTiming(computestart,computeend,computedt);
  //-Uses common configuration.
  gau->SetSaveVtkPart(CfgDefault.savevtkpart);
  gau->ConfigOutputTiming(CfgDefault.output,CfgDefault.outputstart
    ,CfgDefault.outputend,CfgDefault.outputdt);
  Gauges.push_back(gau);
  return(gau);
}

//<vs_meeshdat_ini>
//==============================================================================
/// Creates new gauge-Mesh and returns pointer.
//==============================================================================
JGaugeMesh* JGaugeSystem::AddGaugeMesh(std::string name,double computestart
  ,double computeend,double computedt,bool fixed
  ,const jmsh::StMeshBasic& meshbas,std::string outdata,unsigned tfmt
  ,unsigned buffersize,float kclimit,float kcdummy,float masslimit)
{
  if(GetGaugeIdx(name)!=UINT_MAX)Run_Exceptioon(fun::PrintStr("The name \'%s\' already exists.",name.c_str()));
  if(masslimit<=0)masslimit=CSP.massfluid*(CSP.simulate2d? 0.4f: 0.5f);
  //-Creates object.
  JGaugeMesh* gau=new JGaugeMesh(GetCount(),name,meshbas,outdata,tfmt
    ,buffersize,kclimit,kcdummy,masslimit,GpuCount);
  gau->Config(CSP,Scell,ScellDiv,MapPosMin,DomPosMin,DomPosMax);
  gau->ConfigDomMCel(fixed);
  gau->ConfigComputeTiming(computestart,computeend,computedt);
  //-Uses common configuration.
  gau->SetSaveVtkPart(CfgDefault.savevtkpart);
  gau->ConfigOutputTiming(CfgDefault.output,CfgDefault.outputstart
    ,CfgDefault.outputend,CfgDefault.outputdt);
  Gauges.push_back(gau);
  return(gau);
}
//<vs_meeshdat_end>

//==============================================================================
/// Creates new gauge-Force and returns pointer.
//==============================================================================
JGaugeForce* JGaugeSystem::AddGaugeForce(std::string name,double computestart
  ,double computeend,double computedt,bool fixed
  ,const JSphMk* mkinfo,word mkbound)
{
  if(GetGaugeIdx(name)!=UINT_MAX)Run_Exceptioon(fun::PrintStr("The name \'%s\' already exists.",name.c_str()));
  //-Obtains data from mkbound particles.
  const unsigned cmk=mkinfo->GetMkBlockByMkBound(mkbound);
  if(cmk>=mkinfo->Size())Run_Exceptioon(fun::PrintStr("Error loading boundary objects. Mkbound=%u is unknown.",mkbound));
  const JSphMkBlock* mkb=mkinfo->Mkblock(cmk);
  const TpParticles typeparts=mkb->Type;
  if(typeparts!=TpPartFixed && typeparts!=TpPartMoving)
    Run_Exceptioon(fun::PrintStr("Type of boundary particles (Mkbound=%u) is invalid. Only fixed or moving particles are allowed.",mkbound));
  const unsigned idbegin=mkb->Begin;
  const unsigned count=mkb->Count;
  const typecode code=mkb->Code;
  const tfloat3 center=ToTFloat3((mkb->GetPosMin()+mkb->GetPosMax())/TDouble3(2));
  //-Creates object.
  JGaugeForce* gau=new JGaugeForce(GetCount(),name,mkbound,typeparts,idbegin
    ,count,code,center,GpuCount);
  gau->Config(CSP,Scell,ScellDiv,MapPosMin,DomPosMin,DomPosMax);
  gau->ConfigDomMCel(fixed);
  gau->ConfigComputeTiming(computestart,computeend,computedt);
  //-Uses common configuration.
  gau->SetSaveVtkPart(CfgDefault.savevtkpart);
  gau->ConfigOutputTiming(CfgDefault.output,CfgDefault.outputstart
    ,CfgDefault.outputend,CfgDefault.outputdt);
  Gauges.push_back(gau);
  return(gau);
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JGaugeSystem::VisuConfig(std::string txhead,std::string txfoot){
  SaveVtkInitPoints(); //-Includes gauges defined by coding.
  if(!txhead.empty())Log->Print(txhead);
  for(unsigned cg=0;cg<GetCount();cg++){
    const JGaugeItem* gau=Gauges[cg];
    Log->Printf("Guage_%u: \'%s\'",gau->Idx,gau->Name.c_str());
    std::vector<std::string> lines;
    gau->GetConfig(lines);
    for(unsigned i=0;i<unsigned(lines.size());i++)Log->Print(string("  ")+lines[i]);
  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Saves VTK file with points.
//==============================================================================
void JGaugeSystem::SaveVtkInitPoints()const{
  //-Save individual schemes for complex gauges.
  for(unsigned cg=0;cg<GetCount();cg++)Gauges[cg]->SaveVtkScheme();
  //-Save VKT file with all initial gauge points.
  unsigned* vidx=NULL;
  unsigned* vtype=NULL;
  byte* vout=NULL;
  std::vector<tfloat3> points;
  unsigned ndata=0;
  for(unsigned cg=0;cg<GetCount();cg++){
    const unsigned idx=Gauges[cg]->Idx;
    const unsigned type=Gauges[cg]->Type;
    const unsigned np=Gauges[cg]->GetPointDef(points);
    //-Resizes allocated memory.
    vidx =fun::ResizeAlloc(vidx ,ndata,ndata+np);
    vtype=fun::ResizeAlloc(vtype,ndata,ndata+np);
    vout =fun::ResizeAlloc(vout ,ndata,ndata+np);
    for(unsigned p=0;p<np;p++){
      const unsigned pp=ndata+p;
      vidx[pp]=idx;
      vtype[pp]=type;
      const tdouble3 ps=ToTDouble3(points[pp]);
      vout[pp]=(DomPosMin<=ps && ps<DomPosMax? 0: 1);
    }
    ndata+=np;
  }
  //-Prepares data.
  JDataArrays arrays;
  arrays.AddArray("Pos",ndata,points.data(),false);
  arrays.AddArray("Idx",ndata,vidx,false);
  arrays.AddArray("Type",ndata,vtype,false);
  arrays.AddArray("Out",ndata,vout,false);
  const string filevtk=AppInfo.GetDirOut()+"CfgGauge_InitPoints.vtk";
  Log->AddFileInfo(filevtk,"Saves points used for gauge calculations (by JGaugeSystem).");
  JSpVtkData::Save(filevtk,arrays,"Pos");
  arrays.Reset();
  //-Frees memory.
  delete[] vidx;
  delete[] vtype;
  delete[] vout;
}

//==============================================================================
/// Returns idx of gauge with indicated name (UINT_MAX: name did not exist).
/// Devuelve idx del gauge con el nombre indicado (UINT_MAX: no existe).
//==============================================================================
unsigned JGaugeSystem::GetGaugeIdx(const std::string& name)const{
  unsigned idx=0;
  const unsigned n=GetCount();
  for(;idx<n && Gauges[idx]->Name!=name;idx++);
  return(idx<n? idx: UINT_MAX);
}

//==============================================================================
/// Returns the requested gauge object.
/// Devuelve el objeto gauge solicitado.
//==============================================================================
JGaugeItem* JGaugeSystem::GetGauge(unsigned c)const{
  if(c>=GetCount())Run_Exceptioon("The requested gauge is not valid.");
  return(Gauges[c]);
}

//==============================================================================
/// Configures particle arrays (on CPU).
//==============================================================================
void JGaugeSystem::ConfigArraysCpu(const acdouble3* pos,const actypecode* code
  ,const acuint* idp,const acfloat4* velrho)
{
  DataCpu.ConfigArrays(pos,code,idp,velrho);
}

//==============================================================================
/// Updates results on gauges (on CPU).
//==============================================================================
void JGaugeSystem::CalculeCpu(double timestep,const StDivDataCpu& dvd
  ,unsigned npbok,unsigned npb,unsigned np,bool savedivstate)
{
  const unsigned ng=GetCount();
  DataCpu.SetDivState(timestep,dvd,npbok,npb,np);
  for(unsigned cg=0;cg<ng;cg++){
    JGaugeItem* gau=Gauges[cg];
    if(gau->Update(timestep))gau->CalculeCpu(DataCpu);
  }
  //-Clear divide state.
  if(!savedivstate)DataCpu.ClearDivState();
}

//==============================================================================
/// Updates results on requested gauge using last input data. (on CPU).
//==============================================================================
void JGaugeSystem::CalculeLastInputCpu(std::string gaugename){
  const unsigned idx=GetGaugeIdx(gaugename);
  if(idx==UINT_MAX)Run_Exceptioon(fun::PrintStr("Requested gauge \'%s\' is missing.",gaugename.c_str()));
  if(!DataCpu.divstate)Run_Exceptioon(fun::PrintStr("Divide state to compute gauge \'%s\' is not available.",gaugename.c_str()));
  Gauges[idx]->CalculeCpu(DataCpu);
}

#ifdef _WITHGPU
//==============================================================================
/// Configures particle arrays (on GPU).
//==============================================================================
void JGaugeSystem::ConfigArraysGpu(int id,const agdouble2* posxy
  ,const agdouble* posz,const agtypecode* code,const aguint* idp
  ,const agfloat4* velrho)
{
  if(id>=GpuCount)Run_Exceptioon("Id is invalid.");
  DataGpu[id].ConfigArrays(posxy,posz,code,idp,velrho);
}

//==============================================================================
/// Updates results on gauges (on GPU).
//==============================================================================
void JGaugeSystem::CalculeGpu(int id,double timestep,const StDivDataGpu& dvd
  ,unsigned npbok,unsigned npb,unsigned np,bool savedivstate)
{
  //-Compute measures.
  const unsigned ng=GetCount();
  DataGpu[id].SetDivState(timestep,dvd,npbok,npb,np);
  for(unsigned cg=0;cg<ng;cg++){
    JGaugeItem* gau=Gauges[cg];
    if(gau->Update(timestep))gau->CalculeGpu(DataGpu[id]);
  }
  //-Clear divide state.
  if(!savedivstate)DataGpu[id].ClearDivState();
}

//==============================================================================
/// Updates results on requested gauge using last input data. (on GPU).
//==============================================================================
void JGaugeSystem::CalculeLastInputGpu(int id,std::string gaugename){
  const unsigned idx=GetGaugeIdx(gaugename);
  if(idx==UINT_MAX)Run_Exceptioon(fun::PrintStr("Requested gauge \'%s\' is missing.",gaugename.c_str()));
  if(!DataCpu.divstate)Run_Exceptioon(fun::PrintStr("Input state to compute gauge \'%s\' is not available.",gaugename.c_str()));
  Gauges[idx]->CalculeGpu(DataGpu[id]);
}
#endif

//==============================================================================
/// Updates results on requested gauge using last input data. (on CPU or GPU).
//==============================================================================
void JGaugeSystem::CalculeLastInput(int id,std::string gaugename){
  if(Cpu)CalculeLastInputCpu(gaugename);
  #ifdef _WITHGPU
  else{
    if(GpuCount>1)Run_Exceptioon("Error id look for \'gaugesystem->CalculeLastInput(0,\'");
    CalculeLastInputGpu(id,gaugename);
  }
  #endif
}

//==============================================================================
/// Saves results in VTK and/or CSV file.
//==============================================================================
void JGaugeSystem::SaveResults(unsigned cpart){
  const unsigned ng=GetCount();
  for(unsigned cg=0;cg<ng;cg++)Gauges[cg]->SaveResults(cpart);
}

