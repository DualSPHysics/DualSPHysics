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

/// \file JSphInOutVel.cpp \brief Implements the class \ref JSphInOutVel.

#include "JSphInOutVel.h"
#include "Functions.h"
#include "FunctionsMath.h"
#include "JAppInfo.h"
#include "JLog2.h"
#include "JXml.h"
#include "JLinearValue.h"
#include "JSphInOutVelAwas.h"
#include "JSphInOutGridData.h"

#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JSphGpu_InOut_iker.h"
#endif

#include <cstring>
#include <cfloat>
#include <climits>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSphInOutVel
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutVel::JSphInOutVel(bool cpu,unsigned idzone,const StCteSph &csp
  ,tdouble3 direction,tdouble3 ptplane,tdouble3 zoneposmin,tdouble3 zoneposmax)
  :Log(AppInfo.LogPtr()),Cpu(cpu),IdZone(idzone),CSP(csp),Direction(direction)
  ,PtPlane(ptplane),ZonePosMin(zoneposmin),ZonePosMax(zoneposmax)
{
  ClassName="JSphInOutVel";

  InputTimeVel=NULL;
  AwasVel=NULL;
  InputVelGrid=NULL;

  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOutVel::~JSphInOutVel(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInOutVel::Reset(){
  VelMode=InVelM_Fixed;
  VelProfile=InVelP_Uniform;
  VelBehaviour=InVelB_Unknown;

  VelMin=VelMax=0;

  CurrentTime=DBL_MAX;
  CurrentCoefs0=CurrentCoefs1=TFloat4(0);

  //-Velocity definition for Fixed and Variable velocity.
  InputVel=InputVel2=InputVel3=0;
  InputVelPosz=InputVelPosz2=InputVelPosz3=0;

  delete InputTimeVel; InputTimeVel=NULL;
  TimeVelIdx0=TimeVelIdx1=UINT_MAX; 

  FlowActive=false;
  FlowRatio=0;
  FlowPointsOk=0;
  FlowToVel=1.f;

  SaveVelProfile=true;

  delete AwasVel; AwasVel=NULL;

  delete InputVelGrid; InputVelGrid=NULL;

}

//==============================================================================
/// Reads initial configuration in the XML node.
//==============================================================================
TpInVelMode JSphInOutVel::ReadXml(const JXml *sxml,TiXmlElement *ele
  ,const std::string &dirdatafile,JGaugeSystem *gaugesystem,double maprealposminy)
{
  TiXmlElement *xele=ele->FirstChildElement("imposevelocity");
  if(xele){
    const unsigned mode=sxml->GetAttributeUint(xele,"mode",true);
    switch(mode){
      case 0:  VelMode=InVelM_Fixed;         break;
      case 1:  VelMode=InVelM_Variable;      break;
      case 2:  VelMode=InVelM_Extrapolated;  break;
      case 3:  VelMode=InVelM_Interpolated;  break;
      default: sxml->ErrReadAtrib(xele,"mode",false,"Value is not valid.");
    }
//-Fixed velocity.
    if(VelMode==InVelM_Fixed){
      const byte vel1=(sxml->ExistsElement(xele,"velocity" )? 1: 0);
      const byte vel2=(sxml->ExistsElement(xele,"velocity2")? 1: 0);
      const byte vel3=(sxml->ExistsElement(xele,"velocity3")? 1: 0);
      if(vel1+vel2+vel3>1)sxml->ErrReadElement(xele,"velocity",false,"Several definitions for velocity were found.");
      if(vel1 || vel1+vel2+vel3==0){
        VelProfile=InVelP_Uniform;
        sxml->CheckAttributeNames(xele,"velocity","time v comment units_comment");
        InputVel=sxml->ReadElementFloat(xele,"velocity","v");
      }
      if(vel2){
        VelProfile=InVelP_Linear;
        sxml->CheckAttributeNames(xele,"velocity2","time v v2 z z2 comment units_comment");
        InputVel=sxml->ReadElementFloat(xele,"velocity2","v");
        InputVel2=sxml->ReadElementFloat(xele,"velocity2","v2");
        InputVelPosz=sxml->ReadElementFloat(xele,"velocity2","z");
        InputVelPosz2=sxml->ReadElementFloat(xele,"velocity2","z2");
      }
      if(vel3){
        VelProfile=InVelP_Parabolic;
        sxml->CheckAttributeNames(xele,"velocity3","time v v2 v3 z z2 z3 comment units_comment");
        InputVel =sxml->ReadElementFloat(xele,"velocity3","v");
        InputVel2=sxml->ReadElementFloat(xele,"velocity3","v2");
        InputVel3=sxml->ReadElementFloat(xele,"velocity3","v3");
        InputVelPosz =sxml->ReadElementFloat(xele,"velocity3","z");
        InputVelPosz2=sxml->ReadElementFloat(xele,"velocity3","z2");
        InputVelPosz3=sxml->ReadElementFloat(xele,"velocity3","z3");
      }
      sxml->CheckElementNames(xele,true,"*velocity *velocity2 *velocity3 *flowvelocity");
    }
//-Variable velocity.
    else if(VelMode==InVelM_Variable){
      const unsigned NTMIN=50;
      const byte vel1= (sxml->ExistsElement(xele,"velocitytimes" )? 1: 0);
      const byte vel2= (sxml->ExistsElement(xele,"velocitytimes2")? 1: 0);
      const byte vel3= (sxml->ExistsElement(xele,"velocitytimes3")? 1: 0);
      const byte vel1f=(sxml->ExistsElement(xele,"velocityfile"  )? 1: 0);
      const byte vel2f=(sxml->ExistsElement(xele,"velocityfile2" )? 1: 0);
      const byte vel3f=(sxml->ExistsElement(xele,"velocityfile3" )? 1: 0);
      if(vel1+vel2+vel3+vel1f+vel2f+vel3f>1)sxml->ErrReadElement(xele,"velocitytimes/velocityfile",false,"Several definitions for velocity were found.");
      if(vel1+vel2+vel3+vel1f+vel2f+vel3f==0)sxml->ErrReadElement(xele,"velocitytimes/velocityfile",false,"No definitions for variable velocity were found.");
      if(vel1 || vel1f){
        VelProfile=InVelP_Uniform;
        InputTimeVel=new JLinearValue(1);
        TiXmlElement* xlis=xele->FirstChildElement("velocitytimes");
        if(xlis){
          InputTimeVel->SetSize(NTMIN);
          TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
          while(elet){
            sxml->CheckAttributeNames(elet,"time v");
            double t=sxml->GetAttributeDouble(elet,"time");
            double v=sxml->GetAttributeDouble(elet,"v");
            InputTimeVel->AddTimeValue(t,v);
            elet=elet->NextSiblingElement("timevalue");
          }
        }
        else{
          InputTimeVel->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"velocityfile","file"));
        }
      }
      if(vel2 || vel2f){
        VelProfile=InVelP_Linear;
        InputTimeVel=new JLinearValue(4);
        TiXmlElement* xlis=xele->FirstChildElement("velocitytimes2");
        if(xlis){
          InputTimeVel->SetSize(NTMIN);
          TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
          while(elet){
            sxml->CheckAttributeNames(elet,"time v v2 z z2");
            double t=sxml->GetAttributeDouble(elet,"time");
            double v=sxml->GetAttributeDouble(elet,"v");
            double v2=sxml->GetAttributeDouble(elet,"v2");
            double z=sxml->GetAttributeDouble(elet,"z");
            double z2=sxml->GetAttributeDouble(elet,"z2");
            InputTimeVel->AddTimeValue(t,v,v2,z,z2);
            elet=elet->NextSiblingElement("timevalue");
          }
        }
        else{
          InputTimeVel->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"velocityfile2","file"));
        }
      }
      if(vel3 || vel3f){
        VelProfile=InVelP_Parabolic;
        InputTimeVel=new JLinearValue(6);
        TiXmlElement* xlis=xele->FirstChildElement("velocitytimes3");
        if(xlis){
          InputTimeVel->SetSize(NTMIN);
          TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
          while(elet){
            sxml->CheckAttributeNames(elet,"time v v2 v3 z z2 z3");
            double t=sxml->GetAttributeDouble(elet,"time");
            double v=sxml->GetAttributeDouble(elet,"v");
            double v2=sxml->GetAttributeDouble(elet,"v2");
            double v3=sxml->GetAttributeDouble(elet,"v3");
            double z=sxml->GetAttributeDouble(elet,"z");
            double z2=sxml->GetAttributeDouble(elet,"z2");
            double z3=sxml->GetAttributeDouble(elet,"z3");
            InputTimeVel->AddTimeValue(t,v,v2,v3,z,z2,z3);
            elet=elet->NextSiblingElement("timevalue");
          }
        }
        else{
          InputTimeVel->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"velocityfile3","file"));
        }
      }
      if(!InputTimeVel || InputTimeVel->GetCount()<1)Run_Exceptioon("There are not inlet/outlet velocity values.");
      sxml->CheckElementNames(xele,true,"*velocitytimes *velocitytimes2 *velocitytimes3 velocityfile velocityfile2 velocityfile3 *flowvelocity");
    }
//-Interpolated velocity with grid-data (old version).
    else if(VelMode==InVelM_Interpolated){
      const string checklist="gridveldata gridposzero awas";
      sxml->CheckElementNames(xele,true,checklist);
      InputVelGrid=new JSphInOutGridData();
      InputVelGrid->ConfigFromFile(dirdatafile+sxml->ReadElementStr(xele,"gridveldata","file"));
      double xmin=sxml->ReadElementDouble(xele,"gridposzero","x",true,0);
      double zmin=sxml->ReadElementDouble(xele,"gridposzero","z",true,0);
      InputVelGrid->SetPosMin(TDouble3(xmin,(CSP.simulate2d? CSP.simulate2dposy: maprealposminy),zmin));
      //-Load AWAS configuration.
      if(sxml->CheckElementActive(xele,"awas")){
        AwasVel=new JSphInOutVelAwas(IdZone,PtPlane.x,Direction,CSP.gravity.z
          ,dirdatafile,gaugesystem,sxml,xele->FirstChildElement("awas"));
      }
    }
    else if(VelMode!=InVelM_Extrapolated)Run_Exceptioon("Inlet/outlet velocity profile is unknown.");
    //-Loads flow configuration.
    {
      sxml->CheckAttributeNames(xele,"flowvelocity","active ratio comment");
      FlowActive=sxml->ReadElementBool(xele,"flowvelocity","active",true,false);
      if(FlowActive){
        if(VelMode!=InVelM_Fixed && VelMode!=InVelM_Variable)sxml->ErrReadElement(xele,"flowvelocity",false,"The use of flow velocity is only supported by fixed or variable velocity.");
        if(VelProfile!=InVelP_Uniform)sxml->ErrReadElement(xele,"flowvelocity",false,"The use of flow velocity is only supported by uniform velocity profile.");
        FlowRatio=sxml->ReadElementFloat(xele,"flowvelocity","ratio",true,1.f);
      }
    }
  }
  else sxml->ErrReadElement(ele,"imposevelocity",true);

  //-Defines velocity limits and behaviour according to the velocity limits.
  DefineBehaviour();
  //-Computes initial velocity data.
  ComputeInitialVel();
  return(VelMode);
}

//==============================================================================
/// Calculates minimum and maximum velocity according velocity configuration.
/// Returns -FLT_MAX and FLT_MAX when velocity is unknown.
//==============================================================================
void JSphInOutVel::CalculateVelMinMax(float &velmin,float &velmax)const{
  velmin=FLT_MAX;
  velmax=-FLT_MAX;
  if(VelMode==InVelM_Fixed){
    switch(VelProfile){
      case InVelP_Uniform:
        velmin=velmax=InputVel;  
      break;
      case InVelP_Linear:    
        velmin=min(InputVel,InputVel2);
        velmax=max(InputVel,InputVel2);  
      break;
      case InVelP_Parabolic:    
        velmin=min(InputVel,min(InputVel2,InputVel3));
        velmax=max(InputVel,max(InputVel2,InputVel3));
      break;
      default: Run_Exceptioon("Velocity profile is unknown.");
    }
  }
  else if(VelMode==InVelM_Variable){
    const unsigned count=InputTimeVel->GetCount();
    switch(VelProfile){
      case InVelP_Uniform:
        for(unsigned c=0;c<count;c++){
          const float v=float(InputTimeVel->GetValueByIdx(c));
          velmin=min(velmin,v);
          velmax=max(velmax,v);
        }
      break;
      case InVelP_Linear:    
        for(unsigned c=0;c<count;c++){
          const float v0=float(InputTimeVel->GetValueByIdx(c,0));
          const float v1=float(InputTimeVel->GetValueByIdx(c,1));
          velmin=min(velmin,min(v0,v1));
          velmax=max(velmax,max(v0,v1));
        }
      break;
      case InVelP_Parabolic:    
        for(unsigned c=0;c<count;c++){
          const float v0=float(InputTimeVel->GetValueByIdx(c,0));
          const float v1=float(InputTimeVel->GetValueByIdx(c,1));
          const float v2=float(InputTimeVel->GetValueByIdx(c,2));
          velmin=min(velmin,min(v0,min(v1,v2)));
          velmax=max(velmax,max(v0,max(v1,v2)));
        }
      break;
      default: Run_Exceptioon("Velocity profile is unknown.");
    }
  }
  else if(VelMode==InVelM_Extrapolated || VelMode==InVelM_Interpolated){
    velmin=-FLT_MAX;
    velmax=FLT_MAX;
  }
  else Run_Exceptioon("Velocity mode is unknown.");
}

//==============================================================================
/// Defines behaviour according to the velocity limits.
//==============================================================================
void JSphInOutVel::DefineBehaviour(){
  //-Calculates velocity limits.
  CalculateVelMinMax(VelMin,VelMax);
  //-Defines behaviour according to the velocity limits.
       if(VelMin>=0 && VelMax>=0)VelBehaviour=InVelB_Inlet;
  else if(VelMin<=0 && VelMax<=0)VelBehaviour=InVelB_Outlet;
  else if(VelMin<0 && VelMax>=0 && VelMin!=-FLT_MAX && VelMax!=FLT_MAX)VelBehaviour=InVelB_Reverse;
  else VelBehaviour=InVelB_Unknown;
}

//==============================================================================
/// Returns behaviour as string.
//==============================================================================
std::string JSphInOutVel::GetInletBehaviourName()const{
  string tx="Unknown";
  switch(VelBehaviour){
    case InVelB_Inlet:   tx="Inlet (velocity>=0)";  break;
    case InVelB_Outlet:  tx="Outlet behaviour (velocity<0)";  break;
    case InVelB_Reverse: tx="Inlet & Outlet behaviour (velocity>0 and velocity<0)";  break;
  }
  return(tx);
}

//==============================================================================
/// Configures the conversion from velocity m/s to flow l/s.
//==============================================================================
void JSphInOutVel::ConfigFlowToVel(unsigned initnptok){
  FlowPointsOk=initnptok;
  const double volpart=(CSP.simulate2d? CSP.dp*CSP.dp: CSP.dp*CSP.dp*CSP.dp)*FlowRatio; //-Volume of particle.
  const double npartlitre=(CSP.simulate2d? 0.01: 0.001)/volpart; //-Particles per litre.
  FlowToVel=float(((npartlitre/FlowPointsOk)*CSP.dp));
  //-Updates velocity using modified FlowToVel.
  if(VelMode==InVelM_Variable)TimeVelIdx0=TimeVelIdx1=UINT_MAX; 
  ComputeInitialVel();
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphInOutVel::GetConfig(std::vector<std::string> &lines)const{
  const bool simulate2d=CSP.simulate2d;
  lines.push_back(fun::PrintStr("Velocity mode: %s",TpInVelModeText(VelMode)));
  if(VelMode==InVelM_Fixed || VelMode==InVelM_Variable){
    if(FlowActive){
      lines.push_back(fun::PrintStr("  Flow velocity configuration: True  (flow to velocity: %g m/l  volratio:%g)",FlowToVel,FlowRatio));
      lines.push_back(fun::PrintStr("  Active InOut points: %u",FlowPointsOk));
    }
    if(VelMode==InVelM_Variable && !InputTimeVel->GetFile().empty())lines.push_back(fun::PrintStr("  Velocity file: %s",InputTimeVel->GetFile().c_str()));
    if(FlowActive){
      if(VelProfile==InVelP_Uniform && VelMode==InVelM_Fixed)
        lines.push_back(fun::PrintStr("  Velocity profile: Uniform %g m/s  (%g l/s)",InputVel*FlowToVel,InputVel));
    }
    else{
      if(VelProfile==InVelP_Uniform && VelMode==InVelM_Fixed)
        lines.push_back(fun::PrintStr("  Velocity profile: Uniform %g m/s",InputVel));
    }
    if(VelProfile==InVelP_Linear   )lines.push_back(fun::PrintStr("  Velocity profile: Linear %g(z=%g), %g(z=%g)",InputVel,InputVelPosz,InputVel2,InputVelPosz2));
    if(VelProfile==InVelP_Parabolic)lines.push_back(fun::PrintStr("  Velocity profile: Parabolic %g(z=%g), %g(z=%g), %g(z=%g)",InputVel,InputVelPosz,InputVel2,InputVelPosz2,InputVel3,InputVelPosz3));
  }
  else if(VelMode==InVelM_Interpolated){
    lines.push_back(fun::PrintStr("  Velocity file: %s",InputVelGrid->GetFile().c_str()));
    lines.push_back(fun::PrintStr("  Reset Z velocity: %s","True"));
    if(AwasVel)AwasVel->GetConfig(lines);
  }
  else if(VelMode!=InVelM_Extrapolated)Run_Exceptioon("Velocity mode is unknown.");
}

//==============================================================================
/// Computes the initial velocity data.
//==============================================================================
void JSphInOutVel::ComputeInitialVel(){
  CurrentCoefs0=CurrentCoefs1=TFloat4(0);
  if(VelMode==InVelM_Fixed){
    if(VelProfile==InVelP_Uniform){
      CurrentCoefs0=TFloat4(InputVel*FlowToVel,0,0,0);
    }
    else if(VelProfile==InVelP_Linear){
      const float m=(InputVel2-InputVel)/(InputVelPosz2-InputVelPosz);
      const float b=InputVel-m*InputVelPosz;
      CurrentCoefs0=TFloat4(m,b,0,0);
    }
    else if(VelProfile==InVelP_Parabolic){
      const tmatrix3f inv=fmath::InverseMatrix3x3(TMatrix3f(InputVelPosz*InputVelPosz,InputVelPosz,1,InputVelPosz2*InputVelPosz2,InputVelPosz2,1,InputVelPosz3*InputVelPosz3,InputVelPosz3,1));
      const float a=inv.a11*InputVel+inv.a12*InputVel2+inv.a13*InputVel3;
      const float b=inv.a21*InputVel+inv.a22*InputVel2+inv.a23*InputVel3;
      const float c=inv.a31*InputVel+inv.a32*InputVel2+inv.a33*InputVel3;
      CurrentCoefs0=TFloat4(a,b,c,0);
    }
    else Run_Exceptioon("Inlet/outlet velocity profile is unknown.");
  }
  UpdateVel(0);
}

//==============================================================================
/// Updates velocity coefficients for variable velocity according timestep. 
//==============================================================================
void JSphInOutVel::UpdateVelVariable(double timestep){
  InputTimeVel->FindTime(timestep);
  if(TimeVelIdx0!=InputTimeVel->GetPos() || TimeVelIdx1!=InputTimeVel->GetPosNext()){
    TimeVelIdx0=InputTimeVel->GetPos();
    TimeVelIdx1=InputTimeVel->GetPosNext();
    const float t =(float)InputTimeVel->GetTimeByIdx(TimeVelIdx0);
    const float t2=(float)InputTimeVel->GetTimeByIdx(TimeVelIdx1);
    if(VelProfile==InVelP_Uniform){
      CurrentCoefs0=TFloat4((float)InputTimeVel->GetValueByIdx(TimeVelIdx0)*FlowToVel,0,0,t);
      CurrentCoefs1=TFloat4((float)InputTimeVel->GetValueByIdx(TimeVelIdx1)*FlowToVel,0,0,t2);
    }
    else if(VelProfile==InVelP_Linear){
      for(unsigned cc=0;cc<=1;cc++){
        const unsigned idx=(!cc? TimeVelIdx0: TimeVelIdx1);
        const float v =(float)InputTimeVel->GetValueByIdx(idx,0);
        const float v2=(float)InputTimeVel->GetValueByIdx(idx,1);
        const float z =(float)InputTimeVel->GetValueByIdx(idx,2);
        const float z2=(float)InputTimeVel->GetValueByIdx(idx,3);
        const float m=(v2-v)/(z2-z);
        const float b=v-m*z;
        if(!cc)CurrentCoefs0=TFloat4(m,b,0,t);
        else   CurrentCoefs1=TFloat4(m,b,0,t2);
      }
    }
    else if(VelProfile==InVelP_Parabolic){
      for(unsigned cc=0;cc<=1;cc++){
        const unsigned idx=(!cc? TimeVelIdx0: TimeVelIdx1);
        const float v =(float)InputTimeVel->GetValueByIdx(idx,0);
        const float v2=(float)InputTimeVel->GetValueByIdx(idx,1);
        const float v3=(float)InputTimeVel->GetValueByIdx(idx,2);
        const float z =(float)InputTimeVel->GetValueByIdx(idx,3);
        const float z2=(float)InputTimeVel->GetValueByIdx(idx,4);
        const float z3=(float)InputTimeVel->GetValueByIdx(idx,5);
        const tmatrix3f inv=fmath::InverseMatrix3x3(TMatrix3f(z*z,z,1,z2*z2,z2,1,z3*z3,z3,1));
        const float a=inv.a11*v+inv.a12*v2+inv.a13*v3;
        const float b=inv.a21*v+inv.a22*v2+inv.a23*v3;
        const float c=inv.a31*v+inv.a32*v2+inv.a33*v3;
        if(!cc)CurrentCoefs0=TFloat4(a,b,c,t);
        else   CurrentCoefs1=TFloat4(a,b,c,t2);
      }
    }
    else Run_Exceptioon("Inlet/outlet velocity profile is unknown.");
  }
}

//==============================================================================
/// Updates velocity according to the timestep. 
//==============================================================================
void JSphInOutVel::UpdateVel(double timestep){
  if(InputTimeVel)UpdateVelVariable(timestep);
  //else if(VelMode==InVelM_Interpolated && InputVelGrid) Nothing to do.
  CurrentTime=timestep;
}

//==============================================================================
/// Applies interpolated velocity to inout fluid according to the timestep on CPU. 
//==============================================================================
void JSphInOutVel::UpdateVelInterpolateCpu(double timestep,unsigned nplist
  ,const int *plist,const tdouble3 *pos,const typecode *code
  ,const unsigned *idp,tfloat4 *velrhop)
{
  if(InputVelGrid){
    const float velcorr=(AwasVel? AwasVel->GetVelCorr(timestep): 0);
    if(InputVelGrid->GetNx()==1)
         InputVelGrid->InterpolateZVelCpu(timestep,byte(IdZone),nplist,plist,pos,code,idp,velrhop,velcorr);
    else InputVelGrid->InterpolateVelCpu (timestep,byte(IdZone),nplist,plist,pos,code,idp,velrhop,velcorr);
    //-Removes interpolated Z velocity of inlet/outlet particles.
    if(InputVelGrid->GetUseVelz()){
      const int n=int(nplist);
      #ifdef OMP_USE
        #pragma omp parallel for schedule (static) if(n>OMP_LIMIT_COMPUTELIGHT)
      #endif
      for(int cp=0;cp<n;cp++){
        const unsigned p=plist[cp];
        if(IdZone==CODE_GetIzoneFluidInout(code[p]))velrhop[p].z=0;
      }
    }
  }
}

#ifdef _WITHGPU
//==============================================================================
/// Applies interpolated velocity to inout fluid according to the timestep on GPU. 
//==============================================================================
void JSphInOutVel::UpdateVelInterpolateGpu(double timestep,unsigned nplist,const int *plist
  ,const double2 *posxyg,const double *poszg,const typecode *codeg
  ,const unsigned *idpg,float4 *velrhopg)
{
  if(InputVelGrid){
    const float velcorr=(AwasVel? AwasVel->GetVelCorr(timestep): 0);
    if(InputVelGrid->GetNx()==1)
         InputVelGrid->InterpolateZVelGpu(timestep,byte(IdZone),nplist,plist,posxyg,poszg,codeg,idpg,velrhopg,velcorr);
    else Run_Exceptioon("GPU code was not implemented for nx>1.");
    //-Removes interpolated Z velocity of inlet/outlet particles.
    if(InputVelGrid->GetUseVelz()){
      cusphinout::InOutInterpolateResetZVel(IdZone,nplist,plist,codeg,velrhopg);
    }
  }
}
#endif

//==============================================================================
/// Saves CSV file with AWAS information.
//==============================================================================
void JSphInOutVel::SaveAwasVelCsv(){
  if(AwasVel)AwasVel->SaveCsvData();
}

//==============================================================================
/// Saves grid nodes in VTK file.
//==============================================================================
void JSphInOutVel::SaveVtkVelGrid(){
  if(InputVelGrid){
    const string filevtk=AppInfo.GetDirOut()+fun::PrintStr("CfgInOut_VelGrid_i%d.vtk",IdZone);
    Log->AddFileInfo(filevtk,"Saves VTK file with InputVelGrid nodes (by JSphInOut).");
    InputVelGrid->SaveDataVtk(filevtk,0);
  }
}

