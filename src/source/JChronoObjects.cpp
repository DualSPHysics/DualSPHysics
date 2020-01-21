/*
 <DUALSPHYSICS>  Copyright (c) 2016, Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or (at your option) any later version. 

 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. 

 You should have received a copy of the GNU General Public License, along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JChronoObjects.cpp \brief Implements the class \ref JChronoObjects.

#include "JChronoObjects.h"
#include "DSPHChronoLib.h"
#include "JChronoData.h"
#include "Functions.h"
#include "FunctionsGeo3d.h"
#include "JLog2.h"
#include "JXml.h"
#include "JSpaceParts.h"
#include "JAppInfo.h"
#include "JSaveCsv2.h"
#include "JVtkLib.h"
#include "JRangeFilter.h"
#include "JSpaceVtkOut.h"
#include "JSphMk.h"
#include <cstring>
#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

#ifndef DISABLE_CHRONO

//##############################################################################
//# JChronoObjects
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JChronoObjects::JChronoObjects(JLog2* log,const std::string &dirdata,const std::string &casename
 ,JXml *sxml,const std::string &place,double dp,word mkboundfirst)
 :Log(log),DirData(dirdata),CaseName(casename),Dp(dp),MkBoundFirst(mkboundfirst),UseDVI(true)
{
  ClassName="JChronoObjects";
  ChronoDataXml=NULL;
  ChronoLib=NULL;
  Ptr_VtkSimple_AutoActual=NULL;
  Ptr_VtkSimple_AutoDp=NULL;
  Reset();
  LoadXml(sxml,place);
}

//==============================================================================
/// Destructor.
//==============================================================================
JChronoObjects::~JChronoObjects(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JChronoObjects::Reset(){
  WithMotion=false;
  delete ChronoDataXml; ChronoDataXml=NULL;
  delete ChronoLib; ChronoLib=NULL;
  CollisionDp=0;
  SchemeScale=1;
  SaveDataTime=NextTime=0;
  LastTimeOk=-1;
  JVtkLib::DeleteMkShapes(Ptr_VtkSimple_AutoActual); Ptr_VtkSimple_AutoActual=NULL;
  JVtkLib::DeleteMkShapes(Ptr_VtkSimple_AutoDp);     Ptr_VtkSimple_AutoDp=NULL;
}

//==============================================================================
/// Returns TRUE when a chrono object with geometry is defined for indicated mkbound.
//==============================================================================
bool JChronoObjects::UseDataDVI(word mkbound)const{
  const unsigned ipos=ChronoDataXml->BodyIndexByMkBound(mkbound);
  return(ipos<ChronoDataXml->GetBodyCount() && !ChronoDataXml->GetBody(ipos)->GetModelFile().empty());
}

//==============================================================================
/// Configures body with data from floating body. 
/// Returns TRUE when it is a floating with Chrono.
//==============================================================================
bool JChronoObjects::ConfigBodyFloating(word mkbound,double mass
  ,const tdouble3 &center,const tmatrix3d &inertia
  ,const tint3 &translationfree,const tint3 &rotationfree
  ,const tfloat3 &linvelini,const tfloat3 &angvelini)
{
  JChBodyFloating* body=(JChBodyFloating*)ChronoDataXml->GetBodyFloating(mkbound);
  if(body)body->SetFloatingData(mass,center,inertia,translationfree,rotationfree,linvelini,angvelini);
  return(body!=NULL);
}

//==============================================================================
/// Configures Floating body with DVI parameters. 
//==============================================================================
void JChronoObjects::ConfigDataDVIBodyFloating(word mkbound,float kfric,float restitu){
  JChBodyFloating* body=(JChBodyFloating*)ChronoDataXml->GetBodyFloating(mkbound);
  if(body)body->SetDVIData(kfric,restitu); 
}

//==============================================================================
/// Configures Moving body with DVI parameters. 
//==============================================================================
void JChronoObjects::ConfigDataDVIBodyMoving(word mkbound,float kfric,float restitu){
  JChBodyMoving* body=(JChBodyMoving*)ChronoDataXml->GetBodyMoving(mkbound);
  if(body)body->SetDVIData(kfric,restitu); 
}

//==============================================================================
/// Configures Fixed body with DVI parameters. 
//==============================================================================
void JChronoObjects::ConfigDataDVIBodyFixed(word mkbound,float kfric,float restitu){
  JChBodyFixed* body=(JChBodyFixed*)ChronoDataXml->GetBodyFixed(mkbound);
  if(body)body->SetDVIData(kfric,restitu); 
}

//==============================================================================
/// Loads data from XML file.
//==============================================================================
void JChronoObjects::LoadXml(JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of bodies and links in the XML node.
//==============================================================================
JChBody::TpModelNormal ReadXmlModelNormal(const JXml *sxml,TiXmlElement* ele){
  string mnormal=sxml->GetAttributeStr(ele,"modelnormal",true,"original");
  JChBody::TpModelNormal tnormal;
  if(fun::StrLower(mnormal)=="original")tnormal=JChBody::NorOriginal;
  else if(fun::StrLower(mnormal)=="invert")tnormal=JChBody::NorInvert;
  else if(fun::StrLower(mnormal)=="twoface")tnormal=JChBody::NorTwoFace;
  else sxml->ErrReadAtrib(ele,"modelnormal",false,"The value is invalid");
  return(tnormal);
}

//==============================================================================
/// Reads list of bodies and links in the XML node.
//==============================================================================
std::string JChronoObjects::ReadXmlModelFile(const JXml *sxml,TiXmlElement* ele)const{
  string mfile=sxml->GetAttributeStr(ele,"modelfile",true);
  if(!mfile.empty()){
    mfile=fun::StrReplace(mfile,"[CaseName]",CaseName);
    mfile=fun::StrReplace(mfile,"[casename]",CaseName);
  }
  return(mfile);
}

//==============================================================================
/// Loads pointer model data for AutoActual configuration.
//==============================================================================
void JChronoObjects::LoadPtrAutoActual(const JXml *sxml,std::string xmlrow){
  if(!JVtkLib::Available())Run_Exceptioon("Code for VTK format files is not included in the current compilation, so CHRONO collisions according to VTK geometry are not supported.");
  if(Ptr_VtkSimple_AutoActual==NULL){
    JSpaceVtkOut vtkout;
    vtkout.LoadXml(sxml,"case.execution.vtkout",false);
    std::vector<std::string> vtkfiles;
    vtkout.GetFiles("_Actual.vtk",vtkfiles);
    if(vtkfiles.size()==0)Run_ExceptioonFile("Actual VTK files not found for configuration AutoActual",xmlrow);
    for(unsigned c=0;c<unsigned(vtkfiles.size());c++)vtkfiles[c]=DirData+vtkfiles[c];
    Ptr_VtkSimple_AutoActual=JVtkLib::CreateMkShapes(vtkfiles);
  }
}

//==============================================================================
/// Loads pointer model data for AutoDp configuration.
//==============================================================================
void JChronoObjects::LoadPtrAutoDp(const JXml *sxml,std::string xmlrow){
  if(!JVtkLib::Available())Run_Exceptioon("Code for VTK format files is not included in the current compilation, so CHRONO collisions according to VTK geometry are not supported.");
  if(Ptr_VtkSimple_AutoDp==NULL){
    JSpaceVtkOut vtkout;
    vtkout.LoadXml(sxml,"case.execution.vtkout",false);
    std::vector<std::string> vtkfiles;
    vtkout.GetFiles("_Dp.vtk",vtkfiles);
    if(vtkfiles.size()==0)Run_ExceptioonFile("Dp VTK files not found for configuration AutoDp",xmlrow);
    for(unsigned c=0;c<unsigned(vtkfiles.size());c++)vtkfiles[c]=DirData+vtkfiles[c];
    Ptr_VtkSimple_AutoDp=JVtkLib::CreateMkShapes(vtkfiles);
  }
}

//==============================================================================
/// Creates and return obj file with body geometry.
//==============================================================================
void JChronoObjects::CreateObjFiles(std::string idname,const std::vector<unsigned> &mkbounds
  ,std::string datadir,std::string mfile,byte normalmode,std::string diroutobj,std::string xmlrow)
{
  const unsigned nmkbounds=unsigned(mkbounds.size());
  JVtkLib::TpModeNormal tnormal=JVtkLib::NorNULL;
  switch((JChBody::TpModelNormal)normalmode){
    case JChBody::NorOriginal: tnormal=JVtkLib::NorOriginal; break;
    case JChBody::NorInvert:   tnormal=JVtkLib::NorInvert;   break;
    case JChBody::NorTwoFace:  tnormal=JVtkLib::NorTwoFace;  break;
    default: Run_ExceptioonFile("Normal mode is unknown.",xmlrow);
  }
  const string mfext=fun::StrLower(fun::GetExtension(mfile));
  const bool auto_actual=fun::StrUpper(mfile)=="AUTOACTUAL";
  const bool auto_dp=fun::StrUpper(mfile)=="AUTODP";
  if(mfext=="obj"){
    if(nmkbounds>1)Run_ExceptioonFile("Model file OBJ is invalid for several MK values.",xmlrow);
    if(tnormal!=JVtkLib::NorOriginal)Run_ExceptioonFile("Only ModelNormal=Original is valid for Model file OBJ.",xmlrow);
    const string filein=datadir+mfile;
    const string fileout=diroutobj+idname+fun::PrintStr("_mkb%04u.obj",word(mkbounds[0]));
    //Log->Printf("-----> filein:[%s] -> fout[%s]",filein.c_str(),fileout.c_str());
    if(!fun::FileExists(filein))Run_ExceptioonFile("Error: File was not found.",filein);
    if(fun::CpyFile(filein,fileout))Run_ExceptioonFile("Error: File could not be created.",fileout);
  }
  else if(mfext=="vtk" || auto_actual || auto_dp){
    void* ptr=NULL;
    if(auto_actual){
      if(Ptr_VtkSimple_AutoActual==NULL)Run_Exceptioon("Error: Ptr_VtkSimple_AutoActual is not defined.");
      ptr=Ptr_VtkSimple_AutoActual;
    }
    if(auto_dp){
      if(Ptr_VtkSimple_AutoDp==NULL)Run_Exceptioon("Error: Ptr_VtkSimple_AutoDp is not defined.");
      ptr=Ptr_VtkSimple_AutoDp;
    }
    const string filein=(ptr==NULL? datadir+mfile: "");
    JVtkLib::CreateOBJsByMk(ptr,filein,diroutobj+idname,mkbounds,MkBoundFirst,tnormal);
  }
  else Run_ExceptioonFile("Model file format is invalid.",xmlrow);
}

//==============================================================================
/// Reads list of bodies and links in the XML node.
//==============================================================================
void JChronoObjects::ReadXml(const JXml *sxml,TiXmlElement* lis){
  Reset();
  //-Checks XML elements.
  sxml->CheckElementNames(lis,false,"savedata schemescale collisiondp bodyfloating bodymoving bodyfixed link_hinge link_spheric link_pointline link_linearspring");

  ChronoDataXml=new JChronoData;
  ChronoDataXml->SetUseNSCChrono(UseDVI);
  ChronoDataXml->SetDp(Dp);
  //-Create dir for OBJ files with geometry.
  const string diroutobj=AppInfo.GetDirOut()+"chrono_objs/";
  fun::MkdirPath(diroutobj);
  Log->AddFileInfo(diroutobj+"XXXX_mkXXXX.obj","Geometry files used for Chrono interaction.");
  ChronoDataXml->SetDataDir(diroutobj);
  //-Loads configuration to save CSV file for debug.
  SaveDataTime=sxml->ReadElementFloat(lis,"savedata","value",true,-1.f);
  //-Loads allowed collision overlap according Dp.
  CollisionDp=sxml->ReadElementFloat(lis,"collisiondp","value",true,0.5f);
  ChronoDataXml->SetCollisionDp(CollisionDp);
  //-Loads scale value to create initial scheme of configuration.
  SchemeScale=sxml->ReadElementFloat(lis,"schemescale","value",true,1);

  //-Loads body elements.
  unsigned nextidbody=0;
  TiXmlElement* ele=lis->FirstChildElement(); 
  while(ele){
    const std::string elename=ele->Value();
    if(elename.length()>4 && elename.substr(0,4)=="body"){
      const string xmlrow=sxml->ErrGetFileRow(ele);
      string idnamebase=sxml->GetAttributeStr(ele,"id");
      //word mkbound=sxml->GetAttributeWord(ele,"mkbound");
      std::vector<unsigned> mkbounds;
      JRangeFilter rg(sxml->GetAttributeStr(ele,"mkbound"));
      //Log->Printf("-----> mkbounds:[%s]",rg.ToString().c_str());
      rg.GetValues(mkbounds);
      const unsigned nmkbounds=unsigned(mkbounds.size());
      //-Creates obj files with geometry.
      const string mfilebase=ReadXmlModelFile(sxml,ele);
      if(fun::StrUpper(mfilebase)=="AUTOACTUAL" && Ptr_VtkSimple_AutoActual==NULL)LoadPtrAutoActual(sxml,xmlrow);
      if(fun::StrUpper(mfilebase)=="AUTODP"     && Ptr_VtkSimple_AutoDp    ==NULL)LoadPtrAutoDp    (sxml,xmlrow);
      const JChBody::TpModelNormal tnormal=ReadXmlModelNormal(sxml,ele);
      if(!mfilebase.empty())CreateObjFiles(idnamebase,mkbounds,DirData,mfilebase,byte(tnormal),diroutobj,xmlrow);
      //-Creates a body object for each MK value in mkbounds[].
      for(unsigned cmk=0;cmk<nmkbounds;cmk++){
        const word mkbound=word(mkbounds[cmk]);
        //-Creates body object.
        unsigned idb=nextidbody; nextidbody++;
        const string idname=(nmkbounds>1? idnamebase+fun::UintStr(mkbound): idnamebase);
        const string mfile=(!mfilebase.empty()? idnamebase+fun::PrintStr("_mkb%04u.obj",mkbound): ""); 
        if(elename=="bodyfloating"){
          //Log->Printf("----> AddBodyFloating>> \'%s\' mkb:%u mf:[%s]",idname.c_str(),mkbound,mfile.c_str());
          JChBodyFloating *body=ChronoDataXml->AddBodyFloating(idb,idname,mkbound,xmlrow);
          body->SetModel(mfile,tnormal);
          ReadXmlValues(sxml,ele->FirstChildElement("values"),body->GetValuesPtr());
        }
        else if(elename=="bodymoving"){
          Run_Exceptioon("The use of predefined moving objects (<bodymoving>) is disabled since it does not work properly.");
          //Log->Printf("----> AddBodyMoving>> \'%s\' mkb:%u mf:[%s]",idname.c_str(),mkbound,mfile.c_str());
          const double mass=sxml->GetAttributeDouble(ele,"massbody");
          JChBodyMoving *body=ChronoDataXml->AddBodyMoving(idb,idname,mkbound,mass,xmlrow);
          body->SetModel(mfile,tnormal);
          ReadXmlValues(sxml,ele->FirstChildElement("values"),body->GetValuesPtr());
          if(!mfile.empty())WithMotion=true;
        }
        else if(elename=="bodyfixed"){
          //Log->Printf("----> AddBodyFixed>> \'%s\' mkb:%u mf:[%s]",idname.c_str(),mkbound,mfile.c_str());
          JChBodyFixed *body=ChronoDataXml->AddBodyFixed(idb,idname,mkbound,xmlrow);
          body->SetModel(mfile,tnormal);
          ReadXmlValues(sxml,ele->FirstChildElement("values"),body->GetValuesPtr());
        }
        else sxml->ErrReadElement(ele,elename,false);
      }
    }
    ele=ele->NextSiblingElement();
  }

  //-Loads link elements.
  ele=lis->FirstChildElement(); 
  while(ele){
    const std::string elename=ele->Value();
    if(elename.length()>5 && elename.substr(0,5)=="link_"){
      const string xmlrow=sxml->ErrGetFileRow(ele);
      //-Identify body1.
      const string idnamebody1=sxml->GetAttributeStr(ele,"idbody1");
      unsigned idx1=ChronoDataXml->BodyIndexByName(idnamebody1);
      if(idx1==UINT_MAX)Run_ExceptioonFile(fun::PrintStr("idbody1 \'%s\' is not found.",idnamebody1.c_str()),xmlrow);
      const unsigned idbody1=ChronoDataXml->GetBody(idx1)->Idb;
      //-Identify body2.
      const string idnamebody2=sxml->GetAttributeStr(ele,"idbody2",true,"NULL");
      unsigned idbody2=UINT_MAX;
      if(idnamebody2!="NULL"){
        unsigned idx2=ChronoDataXml->BodyIndexByName(idnamebody2);
        if(idx2==UINT_MAX)Run_ExceptioonFile(fun::PrintStr("idbody2 \'%s\' is not found.",idnamebody2.c_str()),xmlrow);
        idbody2=ChronoDataXml->GetBody(idx2)->Idb;
      }
      //-Defines link name.
      string name=sxml->GetAttributeStr(ele,"name",true);
      if(name.empty()){
        if(idbody2!=UINT_MAX)name=fun::PrintStr("Link_%s_%s",idnamebody1.c_str(),idnamebody2.c_str());
        else name=fun::PrintStr("Link_%s",idnamebody1.c_str());
        if(ChronoDataXml->LinkIndexByName(name)!=UINT_MAX){//-Adds number to name when the name already exists.
          int num=2;
          string name2=name+"_"+fun::IntStr(num);
          while(ChronoDataXml->LinkIndexByName(name2)!=UINT_MAX){
            num++;
            name2=name+"_"+fun::IntStr(num);
          }
          name=name2;
        }
      }
      //-Creates link object.
      if(elename=="link_hinge"){
        JChLinkHinge *link=ChronoDataXml->AddLinkHinge(name,idbody1,idbody2,xmlrow);
        link->SetRotPoint (sxml->ReadElementDouble3(ele,"rotpoint" ));
        link->SetRotVector(sxml->ReadElementDouble3(ele,"rotvector"));
        link->SetStiffness(sxml->ReadElementDouble (ele,"stiffness","value"));
        link->SetDamping  (sxml->ReadElementDouble (ele,"damping"  ,"value"));
        ReadXmlValues(sxml,ele->FirstChildElement("values"),link->GetValuesPtr());
      }
      else if(elename=="link_spheric"){
        JChLinkSpheric *link=ChronoDataXml->AddLinkSpheric(name,idbody1,idbody2,xmlrow);
        link->SetRotPoint (sxml->ReadElementDouble3(ele,"rotpoint"));
        link->SetStiffness(sxml->ReadElementDouble (ele,"stiffness","value"));
        link->SetDamping  (sxml->ReadElementDouble (ele,"damping"  ,"value"));
        ReadXmlValues(sxml,ele->FirstChildElement("values"),link->GetValuesPtr());
      }
      else if(elename=="link_pointline"){
        if(idnamebody2!="NULL")Run_ExceptioonFile("Link-PointLine only uses one body.",xmlrow);
        JChLinkPointLine *link=ChronoDataXml->AddLinkPointLine(name,idbody1,idbody2,xmlrow);
        link->SetSlidingVector(sxml->ReadElementDouble3(ele,"slidingvector"));
        link->SetRotPoint (sxml->ReadElementDouble3(ele,"rotpoint"));
        link->SetRotVector(sxml->ExistsElement(ele,"rotvector")? sxml->ReadElementDouble3(ele,"rotvector"): TDouble3(0));
        link->SetRotVector2(sxml->ExistsElement(ele,"rotvector") && sxml->ExistsElement(ele,"rotvector2")? sxml->ReadElementDouble3(ele,"rotvector2"): TDouble3(0));
        link->SetStiffness(sxml->ReadElementDouble (ele,"stiffness","value"));
        link->SetDamping  (sxml->ReadElementDouble (ele,"damping"  ,"value"));
        ReadXmlValues(sxml,ele->FirstChildElement("values"),link->GetValuesPtr());
      }
      else if(elename=="link_linearspring"){
        JChLinkLinearSpring *link=ChronoDataXml->AddLinkLinearSpring(name,idbody1,idbody2,xmlrow);
        link->SetPointfb0  (sxml->ReadElementDouble3(ele,"point_fb1"));
        link->SetPointfb1  (sxml->ReadElementDouble3(ele,"point_fb2"));
        link->SetStiffness (sxml->ReadElementDouble (ele,"stiffness"  ,"value"));
        link->SetDamping   (sxml->ReadElementDouble (ele,"damping"    ,"value"));
        link->SetRestLength(sxml->ReadElementDouble (ele,"rest_length","value"));
        link->SetCoulombDamping(sxml->ReadElementDouble(ele,"coulombdamping","value",true,0));
        ReadXmlValues(sxml,ele->FirstChildElement("values"),link->GetValuesPtr());
        TiXmlElement* ele2=ele->FirstChildElement("savevtk");
        if(ele2){//-Configuration to save vtk spring.
          JChLink::StSaveSpring cfg;
          cfg.radius=sxml->ReadElementFloat(ele2,"radius","value",true,cfg.radius);
          cfg.length=sxml->ReadElementFloat(ele2,"length","value",true,cfg.length);
          cfg.nside =sxml->ReadElementInt  (ele2,"nside" ,"value",true,cfg.nside );
          link->SetSvSpring(cfg);
        }
      }
      else sxml->ErrReadElement(ele,elename,false);
    }
    ele=ele->NextSiblingElement();
  }
  //-Prepares object ChronoDataXml.
  ChronoDataXml->Prepare();
  NextTime=(SaveDataTime>=0? 0: DBL_MAX);
  LastTimeOk=-1;
}

//==============================================================================
/// Reads list of values in the XML node.
//==============================================================================
void JChronoObjects::ReadXmlValues(const JXml *sxml,TiXmlElement* lis,JChValues* values){
  if(lis){
    //-Loads elements bodyfloating.
    TiXmlElement* ele=lis->FirstChildElement(); 
    while(ele){
      const string cmd=ele->Value();
      if(cmd.length()&&cmd[0]!='_'){
        const bool triple=(sxml->ExistsAttribute(ele,"x") && sxml->ExistsAttribute(ele,"y") && sxml->ExistsAttribute(ele,"z"));
        if(cmd=="vstr")values->AddValueStr(sxml->GetAttributeStr(ele,"name"),sxml->GetAttributeStr(ele,"v"));
        else if(cmd=="vint"){
          const string name=sxml->GetAttributeStr(ele,"name");
          if(triple)values->AddValueInt3(name,sxml->GetAttributeInt3(ele));
          else      values->AddValueInt (name,sxml->GetAttributeInt (ele,"v"));
        }
        else if(cmd=="vuint"){
          const string name=sxml->GetAttributeStr(ele,"name");
          if(triple)values->AddValueUint3(name,sxml->GetAttributeUint3(ele));
          else      values->AddValueUint (name,sxml->GetAttributeUint (ele,"v"));
        }
        else if(cmd=="vreal"){
          const string name=sxml->GetAttributeStr(ele,"name");
          if(triple)values->AddValueDouble3(name,sxml->GetAttributeDouble3(ele));
          else      values->AddValueDouble (name,sxml->GetAttributeDouble (ele,"v"));
        }
        else sxml->ErrReadElement(ele,cmd,false);
      }
      ele=ele->NextSiblingElement();
    }
  }
}

//==============================================================================
/// Creates VTK file with the scheme of Chrono objects.
//==============================================================================
void JChronoObjects::SaveVtkScheme()const{
  const double ds=Dp*SchemeScale;
  JVtkLib sh;
  //-Represents floating objects.
  for(unsigned c=0;c<ChronoDataXml->GetBodyCount();c++){
    const JChBody* body=ChronoDataXml->GetBody(c);
    if(body->Type==JChBody::BD_Floating){
      const word mk=body->MkBound+MkBoundFirst;
      const tdouble3 center=body->GetCenter();
      sh.AddShapeCross(center,ds*2,mk);
    }
  }
  //-Represents links.
  for(unsigned c=0;c<ChronoDataXml->GetLinkCount();c++){
    const JChLink *link=ChronoDataXml->GetLink(c);
    if(link->GetBodyRefCount()<1)Run_Exceptioon("Link without body reference.");
    const word mk=link->GetBodyRef(0)->MkBound+MkBoundFirst;
    switch(link->Type){
      case JChLink::LK_Hinge:{
        const JChLinkHinge* linktype=(const JChLinkHinge*)link;
        const tdouble3 pt=linktype->GetRotPoint();
        const tdouble3 v=fgeo::VecUnitary(linktype->GetRotVector())*(ds*2);
        sh.AddShapeCylinder(pt-v,pt+v,ds*1.5,16,mk);
      }break;
      case JChLink::LK_Spheric:{
        const JChLinkSpheric* linktype=(const JChLinkSpheric*)link;
        sh.AddShapeSphere(linktype->GetRotPoint(),ds*2,16,mk);
      }break;
      case JChLink::LK_PointLine:{
        const JChLinkPointLine* linktype=(const JChLinkPointLine*)link;
        const tdouble3 pt=linktype->GetRotPoint();
        const tdouble3 v=fgeo::VecUnitary(linktype->GetSlidingVector())*(ds*4);
        const tdouble3 vrot=fgeo::VecUnitary(linktype->GetRotVector())*(ds*4);
        const tdouble3 vrot2=fgeo::VecUnitary(linktype->GetRotVector2())*(ds*4);
        if(vrot==TDouble3(0))sh.AddShapeSphere(linktype->GetRotPoint(),ds*2,16,mk);
        else{ 
          sh.AddShapeCylinder(pt-vrot,pt+vrot,ds*1.5,16,mk);
          if(vrot2!=TDouble3(0))sh.AddShapeCylinder(pt-vrot2,pt+vrot2,ds*1.5,16,mk);
        }
        sh.AddShapeLine(pt-(v*4.),pt+(v*4.),mk);
      }break;
      case JChLink::LK_LinearSpring:{
        const JChLinkLinearSpring* linktype=(const JChLinkLinearSpring*)link;
        if(link->GetBodyRefCount()<2)Run_Exceptioon("Link without two bodies reference.");
        const word mk1=link->GetBodyRef(1)->MkBound+MkBoundFirst;
        const tdouble3 pt0=linktype->GetPointfb0();
        const tdouble3 pt1=linktype->GetPointfb1();
        const JChLink::StSaveSpring cfg=linktype->GetSvSpring();
        if(cfg.nside>1){
          const double radius=cfg.radius,revlength=cfg.length;
          const int nsides=cfg.nside;
          const double cornersout=cfg.radius/2.f,cornersin=cfg.radius/4.f;
          sh.AddShapeSpring(pt0,pt1,linktype->GetRestLength(),ds,cornersout,cornersin,radius,revlength,nsides,mk);
        }
        else sh.AddShapeLine(pt0,pt1,mk);
        sh.AddShapeBoxSize(pt0-TDouble3(ds),TDouble3(ds*2),mk);
        sh.AddShapeBoxSize(pt1-TDouble3(ds),TDouble3(ds*2),mk1);
      }break;
      default: Run_Exceptioon("Type of link is not supported.");
    }
  }
  const string filevtk=AppInfo.GetDirOut()+"CfgChrono_Scheme.vtk";
  Log->AddFileInfo(filevtk,"Saves VTK file with scheme of Chrono objects and links between objects.");
  sh.SaveShapeVtk(filevtk,"Mk");
}

//==============================================================================
/// Configures center of moving bodies starting from particles domains.
//==============================================================================
void JChronoObjects::ConfigMovingBodies(const JSphMk* mkinfo){
  for(unsigned c=0;c<ChronoDataXml->GetBodyCount();c++){
    const JChBody* body=ChronoDataXml->GetBody(c);
    if(body->Type==JChBody::BD_Moving){
      unsigned cb=mkinfo->GetMkBlockByMkBound(body->MkBound);
      if(cb>=mkinfo->Size())Run_Exceptioon(fun::PrintStr("Center of body \'%s\' (mkbound=%u) is not available.",body->IdName.c_str(),body->MkBound));
      const tdouble3 pcen=(mkinfo->Mkblock(cb)->GetPosMin()+mkinfo->Mkblock(cb)->GetPosMax())/2.;
      ((JChBodyMoving*)body)->SetInitialCenter(pcen);
    }
  }
}

//==============================================================================
/// Configures and reads floating data from XML file.
//==============================================================================
void JChronoObjects::Init(bool simulate2d,const JSphMk* mkinfo){
  //-Updates center of moving objects.
  ConfigMovingBodies(mkinfo);
  //-Checks data in ChronoData.
  ChronoDataXml->CheckData();
  //-Creates VTK file with the scheme of Chrono objects.
  SaveVtkScheme();
  //-Creates and configures object ChronoLib.
  ChronoLib=new DSPHChronoLib();
  const bool svforces=SaveDataTime>=0;
  ChronoLib->Config(AppInfo.GetDirOut(),svforces,simulate2d,*ChronoDataXml);
  if(svforces){
    Log->AddFileInfo("ChronoBody_forces.csv","Saves forces for each body.");
    Log->AddFileInfo("ChronoLink_forces.csv","Saves forces for each link.");
  }
  delete ChronoDataXml; ChronoDataXml=NULL;
  ChronoLib->Config_Inertia();
}

//==============================================================================
/// Shows values configuration using Log.
//==============================================================================
void JChronoObjects::VisuValues(const JChValues *values)const{
  if(values->GetCount()){
    Log->Printf("    Values...: %u",values->GetCount());
    int lenmax=0;
    for(byte mode=0;mode<2;mode++){
      for(unsigned c=0;c<values->GetCount();c++){
        const JChValues::StValue* v=values->GetValue(c);
        string vtx;
        switch(v->type){
          case JChValues::TP_Text:     vtx=v->vtext;                      break;
          case JChValues::TP_Int:      vtx=fun::IntStr(v->vint);          break;
          case JChValues::TP_Uint:     vtx=fun::UintStr(v->vuint);        break;
          case JChValues::TP_Double:   vtx=fun::DoubleStr(v->vdouble);    break;
          case JChValues::TP_Int3:     vtx=fun::Int3Str(v->vint3);        break;
          case JChValues::TP_Uint3:    vtx=fun::Uint3Str(v->vuint3);      break;
          case JChValues::TP_Double3:  vtx=fun::Double3Str(v->vdouble3);  break;
          default: Run_Exceptioon("Type of value is invalid.");
        }
        if(mode==0){
          int len=(int)fun::PrintStr("      %s <%s>",v->name.c_str(),JChValues::TypeToStr(v->type).c_str()).size();
          lenmax=max(lenmax,len);
        }
        else{
          string tx=fun::PrintStr("      %s <%s>",v->name.c_str(),JChValues::TypeToStr(v->type).c_str());
          while(tx.size()<lenmax)tx=tx+".";
          Log->Printf("%s: %s",tx.c_str(),vtx.c_str());
        }
      }
    }
  }
}

//==============================================================================
/// Shows body configuration using Log.
//==============================================================================
void JChronoObjects::VisuBody(const JChBody *body)const{
  Log->Printf("  Body_%04u \"%s\" -  type: %s",body->Idb,body->IdName.c_str(),body->TypeToStr(body->Type).c_str());
  if(body->Type == JChBody::BD_Floating){
    Log->Printf("    MkBound....: %u",((const JChBodyFloating *)body)->MkBound);
    Log->Printf("    Mass.......: %g",body->GetMass());
    Log->Printf("    Center.....: (%s)",fun::Double3gStr(body->GetCenter()).c_str());
    const tmatrix3f inert=ToTMatrix3f(body->GetInertia());
    Log->Printf("    Inertia....: (%g,%g,%g) (xx,yy,zz)",inert.a11,inert.a22,inert.a33);
    Log->Printf("    Inertia....: (%g,%g,%g) (xy,yz,xz)",inert.a12,inert.a23,inert.a13);
    if(!body->GetMotionFree()){
      const tint3 m=body->GetTranslationFree();
      const tint3 r=body->GetRotationFree();
      Log->Printf("    MotionFree.: Translation:(%d,%d,%d) Rotation:(%d,%d,%d)",m.x,m.y,m.z,r.x,r.y,r.z);
    }
    if(body->GetLinearVelini()!=TFloat3(0) || body->GetAngularVelini()!=TFloat3(0)){
      const tfloat3 v=body->GetLinearVelini();
      const tfloat3 a=body->GetAngularVelini();
      if(v!=TFloat3(0))Log->Printf("    LinearVel0.: (%g,%g,%g) [m/s]",v.x,v.y,v.z);
      if(a!=TFloat3(0))Log->Printf("    AngularVel0: (%g,%g,%g) [rad/s]",a.x,a.y,a.z);
    }
  }
  if(body->Type==JChBody::BD_Moving){
    Log->Printf("    MkBound....: %u",((const JChBodyFixed *)body)->MkBound);
    Log->Printf("    Mass.......: %g",body->GetMass());
  }
  if(body->Type==JChBody::BD_Fixed){
    Log->Printf("    MkBound....: %u",((const JChBodyFixed *)body)->MkBound);
  }
  if(!body->GetModelFile().empty()){
    Log->Printf("    Kfric......: %g",body->GetKfric());
    Log->Printf("    Restitution: %g",body->GetRestitu());
    Log->Printf("    ModelFile..: %s",body->GetModelFile().c_str());
    Log->Printf("    ModelNormal: %s",JChBody::NormalToStr(body->GetModelNormal()).c_str());
  }
  VisuValues(body->GetValuesPtr());
  if(body->GetLinkRefCount()){
    Log->Printf("    Links......: %u",body->GetLinkRefCount());
    for(unsigned c=0;c<body->GetLinkRefCount();c++){
      Log->Printf("      %s",body->GetLinkRef(c)->Name.c_str());
    }
  }
}

//==============================================================================
/// Shows link configuration using Log.
//==============================================================================
void JChronoObjects::VisuLink(const JChLink *link)const{
  Log->Printf("  Link \"%s\" -  type: %s",link->Name.c_str(),link->TypeToStr(link->Type).c_str());
  switch(link->Type){
    case JChLink::LK_Hinge:{
      const JChLinkHinge* linktype=(const JChLinkHinge*)link;
      Log->Printf("    Rotation point: (%s)",fun::Double3gStr(linktype->GetRotPoint()).c_str());
      Log->Printf("    Rotation axis.: (%s)",fun::Double3gStr(linktype->GetRotVector()).c_str());
    }break;
    case JChLink::LK_Spheric:{
      const JChLinkSpheric* linktype=(const JChLinkSpheric*)link;
      Log->Printf("    Rotation point: (%s)",fun::Double3gStr(linktype->GetRotPoint()).c_str());
    }break;
    case JChLink::LK_PointLine:{
      const JChLinkPointLine* linktype=(const JChLinkPointLine*)link;
      Log->Printf("    Sliding vector: (%s)",fun::Double3gStr(linktype->GetSlidingVector()).c_str());
      Log->Printf("    Rotation point: (%s)",fun::Double3gStr(linktype->GetRotPoint()).c_str());
      if(linktype->GetRotVector() !=TDouble3(0))Log->Printf("    Rotation axis.: (%s)",fun::Double3gStr(linktype->GetRotVector()).c_str());
      if(linktype->GetRotVector2()!=TDouble3(0))Log->Printf("    Rotation axis2: (%s)",fun::Double3gStr(linktype->GetRotVector2()).c_str());
    }break;
    case JChLink::LK_LinearSpring:{
      const JChLinkLinearSpring* linktype=(const JChLinkLinearSpring*)link;
      Log->Printf("    Point Body 1..: (%s)",fun::Double3gStr(linktype->GetPointfb0()).c_str());
      Log->Printf("    Point Body 2..: (%s)",fun::Double3gStr(linktype->GetPointfb1()).c_str());
      Log->Printf("    Rest length...: %g", linktype->GetRestLength());
      if(linktype->GetCoulombDamping())Log->Printf("    CoulombDamping: %g", linktype->GetCoulombDamping());
    }break;
    default: Run_Exceptioon("Type of link is not supported.");
  }
  Log->Printf("    Stiffness.....: %g",link->GetStiffness());
  Log->Printf("    Damping.......: %g",link->GetDamping());
  VisuValues(link->GetValuesPtr());
  if(link->GetBodyRefCount()){
    Log->Printf("    Bodies........: %u",link->GetBodyRefCount());
    for(unsigned c=0;c<link->GetBodyRefCount();c++){
      Log->Printf("      Body_%04u \"%s\"",link->GetBodyRef(c)->Idb,link->GetBodyRef(c)->IdName.c_str());
    }
  }
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JChronoObjects::VisuConfig(std::string txhead, std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  const JChronoData* chdata=ChronoLib->GetChronoData();
  Log->Printf("  DSPHChrono version: %s",ChronoLib->version.c_str());
  Log->Printf("  Data directory.....: [%s]",chdata->GetDataDir().c_str());
  Log->Printf("  Bodies.............: %d",chdata->GetBodyCount());
  Log->Printf("  Links..............: %u",chdata->GetLinkCount());
  for(unsigned c=0;c<chdata->GetBodyCount();c++)VisuBody(chdata->GetBody(c));
  for(unsigned c=0;c<chdata->GetLinkCount();c++)VisuLink(chdata->GetLink(c));
  if(!txfoot.empty())Log->Print(txfoot);
}

//==============================================================================
/// Loads floating data to calculate coupling with Chrono.
//==============================================================================
void JChronoObjects::SetFtData(word mkbound,const tfloat3 &face,const tfloat3 &fomegaace){
  if(!ChronoLib->SetFtData(mkbound,face,fomegaace))Run_Exceptioon("Error running Chrono library.");
}

//==============================================================================
/// Obtains floating data from coupling with Chrono.
//==============================================================================
void JChronoObjects::GetFtData(word mkbound,tdouble3 &fcenter,tfloat3 &fvel,tfloat3 &fomega)const{
  if(!ChronoLib->GetFtData(mkbound,fcenter,fvel,fomega))Run_Exceptioon("Error running Chrono library.");
}

//==============================================================================
/// Loads motion data to calculate coupling with Chrono.
//==============================================================================
void JChronoObjects::SetMovingData(word mkbound,bool simple,const tdouble3 &msimple,const tmatrix4d &mmatrix,double stepdt){
  if(!ChronoLib->SetMovingData(mkbound,simple,msimple,mmatrix,stepdt))Run_Exceptioon("Error running Chrono library.");
}

//==============================================================================
/// Computes a single timestep with Chrono for the system
//==============================================================================
void JChronoObjects::RunChrono(unsigned nstep,double timestep,double dt,bool predictor){
  if(!ChronoLib->RunChrono(timestep,dt,predictor))Run_Exceptioon("Error running Chrono library.");
  //-Saves floating body data in CSV files.
  if((LastTimeOk==timestep || NextTime<=timestep) && (SaveDataTime==0 || !predictor)){
    const JChronoData* chdata=ChronoLib->GetChronoData();
    for(unsigned cb=0;cb<chdata->GetBodyCount();cb++)if(chdata->GetBody(cb)->Type==JChBody::BD_Floating){
      const JChBodyFloating *body=(JChBodyFloating*)chdata->GetBody(cb);
      const string file=AppInfo.GetDirOut()+fun::PrintStr("ChronoExchange_mkbound_%u.csv",body->MkBound);
      jcsv::JSaveCsv2 scsv(file,true,AppInfo.GetCsvSepComa());
      if(!scsv.GetAppendMode()){
        Log->AddFileInfo("ChronoExchange_mkbound_XX.csv","Saves information of data exchange between DualSPHysics and Chrono library for each body.");
        //-Saves head.
        scsv.SetHead();
        scsv << "nstep;time [s];dt [s];predictor;face.x [m/s^2];face.y [m/s^2];face.z [m/s^2]";
        scsv << "fomegaace.x [rad/s^2];fomegaace.y [rad/s^2];fomegaace.z [rad/s^2];fvel.x [m/s];fvel.y [m/s];fvel.z [m/s]";
        scsv << "fcenter.x [m];fcenter.y [m];fcenter.z [m];fomega.x [rad/s];fomega.y [rad/s];fomega.z [rad/s]" << jcsv::Endl();
      }
      //-Saves data.
      scsv.SetData();
      scsv << nstep << timestep << dt << (predictor? "True": "False");
      scsv << body->GetInputFace();
      scsv << body->GetInputFomegaAce();
      scsv << body->GetOutputVel();
      scsv << body->GetOutputCenter();
      scsv << body->GetOutputOmega();
      scsv << jcsv::Endl();
      scsv.SaveData();
      //-Recalculates NextTime.
      if(LastTimeOk!=timestep){
        if(SaveDataTime>0)while(NextTime<=timestep)NextTime+=SaveDataTime;
        LastTimeOk=timestep;
      }
    }
    //-Saves forces for each body and link (link_forces.csv, body_forces.csv).
    ChronoLib->SaveForces();
  }
  //if(1){
  //  tdouble3 pcen;
  //  ChronoLib->GetBodyCenter("ball",pcen);
  //  Log->Printf("RunChrono----> timestep:%f  dt:%f  ball.center:(%f,%f,%f)",timestep,dt,pcen.x,pcen.y,pcen.z);
  //}
}

//==============================================================================
/// Saves special data for each PART.
//==============================================================================
void JChronoObjects::SavePart(int part){
  const double ds=Dp*SchemeScale;
  const JChronoData* chdata=ChronoLib->GetChronoData();
  //-Saves VTK of LinearSpring links.
  if(1){
    bool save=false;
    JVtkLib sh;
    for(unsigned c=0;c<chdata->GetLinkCount();c++)if(chdata->GetLink(c)->Type==JChLink::LK_LinearSpring){
      const JChLinkLinearSpring* linktype=(const JChLinkLinearSpring*)chdata->GetLink(c);
      if(linktype->GetBodyRefCount()<1)Run_Exceptioon("Link without body reference.");
      const word mk=linktype->GetBodyRef(0)->MkBound+MkBoundFirst;
      tdouble3 p1,p2;
      if(ChronoLib->GetSpringLinkPositions(linktype->Name,p1,p2))Run_Exceptioon("Error running GetSpringLinkPositions() of Chrono library.");
      //Log->Printf("---> SpringPos: (%f,%f,%f) (%f,%f,%f)\n",p1.x,p1.y,p1.z,p2.x,p2.y,p2.z);
      const JChLink::StSaveSpring cfg=linktype->GetSvSpring();
      if(cfg.nside>1){
        const double radius=cfg.radius,revlength=cfg.length;
        const int nsides=cfg.nside;
        const double cornersout=cfg.radius/2.f,cornersin=cfg.radius/4.f;
        //const double restlen=linktype->GetRestLength();
        const double restlen=ChronoLib->GetSpringLinkRestLength(linktype->Name); //-Is necessary when it is variable.
        sh.AddShapeSpring(p1,p2,restlen,ds,cornersout,cornersin,radius,revlength,nsides,mk);
        save=true;
      }
      else if(cfg.nside==1){
        sh.AddShapeLine(p1,p2,mk);
        save=true;
      }
    }
    if(save){
      Log->AddFileInfo("data/Chrono_Springs_????.vtk","Saves VTK file with representation of Chrono springs.");
      const string filevtk=AppInfo.GetDirDataOut()+fun::FileNameSec("Chrono_Springs.vtk",part);
      sh.SaveShapeVtk(filevtk,"Mk");
    }
  }
}

#endif




