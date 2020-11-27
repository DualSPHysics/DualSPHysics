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

/// \file JSphInOutZsurf.cpp \brief Implements the class \ref JSphInOutZsurf.

#include "JSphInOutZsurf.h"
#include "Functions.h"
#include "FunGeo3d.h"
#include "JXml.h"
#include "JLog2.h"
#include "JLinearValue.h"
#include "JDsGaugeSystem.h"
#include "JVtkLib.h"
#include "JAppInfo.h"

#ifdef _WITHGPU
  #include "FunctionsCuda.h"
  #include "JSphGpu_InOut_iker.h"
#endif

#include <cstring>
#include <cfloat>
#include <algorithm>

using namespace std;

//##############################################################################
//# JSphInOutZsurf
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphInOutZsurf::JSphInOutZsurf(bool cpu,unsigned idzone,const StCteSph &csp
  ,tdouble3 direction,tdouble3 zoneposmin,tdouble3 zoneposmax)
  :Log(AppInfo.LogPtr()),Cpu(cpu),IdZone(idzone),CSP(csp),Direction(direction)
  ,ZonePosMin(zoneposmin),ZonePosMax(zoneposmax)
{
  ClassName="JSphInOutZsurf";
  InputTime=NULL;
  Times=NULL;  TimesZsurf=NULL;  TimesZsurfg=NULL;
  CurrentExternal=false;
  CurrentZsurf=NULL;  CurrentZsurfg=NULL;
  GaugeSwl=NULL;
  OldCode=true;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphInOutZsurf::~JSphInOutZsurf(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphInOutZsurf::Reset(){
  ZsurfMode=InZsurf_Undefined;
  UniformZsurf=true;

  InputZbottom=InputZsurf=FLT_MAX;
  RemoveZsurf=false;
  SvVtkZsurf=false;
  ZsurfMin=0;
  ZsurfFit=0;
  delete InputTime; InputTime=NULL;

  memset(&ZsurfResults,0,sizeof(StZsurfResult));

  GaugePt0=GaugePtz=GaugePtx=TDouble3(0);
  Dptz=Dptx=0;
  Nptz=Nptx=0;
  PlaDisx=TPlane3d(0);

  GaugeSwl=NULL;

  ResetTimes();
}

//==============================================================================
/// Reads initial configuration in the XML node.
//==============================================================================
TpInZsurfMode JSphInOutZsurf::ReadXml(const JXml *sxml,TiXmlElement* ele
  ,const std::string &dirdatafile,JGaugeSystem *gaugesystem)
{
  ZsurfMode=InZsurf_Undefined;
  UniformZsurf=true;
  InputZsurf=FLT_MAX;
  TiXmlElement* xele=ele->FirstChildElement("imposezsurf");
  if(xele){
    const unsigned mode=sxml->GetAttributeUint(xele,"mode",true);
    switch(mode){
      case 0:  ZsurfMode=InZsurf_Fixed;       break;
      case 1:  ZsurfMode=InZsurf_Variable;    break;
      case 2:  ZsurfMode=InZsurf_Calculated;  break;
      default: sxml->ErrReadAtrib(xele,"mode",false,"Value is not valid.");
    }
    //-Load general parameters.
    //InputZbottom=sxml->ReadElementFloat(xele,"zbottom","value"); 
    InputZbottom=float(ZonePosMin.z); //-It is automatically calculated.
    RemoveZsurf=sxml->ReadElementBool(xele,"remove","value",true,false);
    if(ZsurfMode==InZsurf_Fixed){
      sxml->CheckElementNames(xele,true,"remove zsurf zbottom");
      if(sxml->ExistsElement(xele,"zbottom"))Log->PrintfWarning("Option \'zbottom\' is ignored at %s",sxml->ErrGetFileRow(xele,"zbottom").c_str());
      InputZsurf=sxml->ReadElementFloat(xele,"zsurf","value");
      UniformZsurf=true;
    }
    else if(ZsurfMode==InZsurf_Variable){
      sxml->CheckElementNames(xele,true,"remove savevtk zsurftimes zsurffile zbottom");
      if(sxml->ExistsElement(xele,"zbottom"))Log->PrintfWarning("Option \'zbottom\' is ignored at %s",sxml->ErrGetFileRow(xele,"zbottom").c_str());
      SvVtkZsurf=sxml->ReadElementBool(xele,"savevtk","value",true,false);
      const byte zlist=(sxml->ExistsElement(xele,"zsurftimes")? 1: 0);
      const byte zfile=(sxml->ExistsElement(xele,"zsurffile" )? 1: 0);
      if(zlist+zfile>1)sxml->ErrReadElement(xele,"zsurftimes/zsurffile",false,"Several definitions for zsurf were found.");
      if(zlist+zfile==0)sxml->ErrReadElement(xele,"zsurftimes/zsurffile",false,"No definitions for variable zsurf were found.");
      //-Loads list of time-values.
      InputTime=new JLinearValue(1);
      TiXmlElement* xlis=xele->FirstChildElement("zsurftimes");
      if(xlis){
        const unsigned NTMIN=50;
        InputTime->SetSize(NTMIN);
        TiXmlElement* elet=xlis->FirstChildElement("timevalue"); 
        while(elet){
          double t=sxml->GetAttributeDouble(elet,"time");
          double v=sxml->GetAttributeDouble(elet,"zsurf");
          InputTime->AddTimeValue(t,v);
          elet=elet->NextSiblingElement("timevalue");
        }
      }
      else{
        InputTime->LoadFile(dirdatafile+sxml->ReadElementStr(xele,"zsurffile","file"));
      }
      InputZsurf=(float)InputTime->GetValue(0);
      UniformZsurf=true;
    }
    else if(ZsurfMode==InZsurf_Calculated){
      sxml->CheckElementNames(xele,true,"remove savevtk zsurfmin zsurffit zsurf zbottom");
      if(sxml->ExistsElement(xele,"zsurf"  ))Log->PrintfWarning("Option \'zsurf\' is ignored at %s"  ,sxml->ErrGetFileRow(xele,"zsurf").c_str());
      if(sxml->ExistsElement(xele,"zbottom"))Log->PrintfWarning("Option \'zbottom\' is ignored at %s",sxml->ErrGetFileRow(xele,"zbottom").c_str());
      //InputZsurf=sxml->ReadElementFloat(xele,"zsurf","value");
      SvVtkZsurf=sxml->ReadElementBool(xele,"savevtk","value",true,false);
      ZsurfMin=sxml->ReadElementDouble(xele,"zsurfmin","value",true,ZonePosMin.z);
      ZsurfFit=-float(CSP.dp/2);
      switch(sxml->CheckAttributes(xele,"zsurffit","value valuedp",true)){
        case 1:  ZsurfFit=-sxml->ReadElementFloat(xele,"zsurffit","value");                  break;
        case 2:  ZsurfFit=-sxml->ReadElementFloat(xele,"zsurffit","valuedp")*float(CSP.dp);  break;
      }
      if(ZsurfMin<ZonePosMin.z)ZsurfMin=ZonePosMin.z;
      UniformZsurf=ConfigGaugeZsurf(gaugesystem);
      if(GaugeSwl){
        gaugesystem->CalculeLastInput(GaugeSwl->Name);
        UpdateZsurf(-1);
      }
      else Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: Initial Zsurf is invalid.",IdZone));
    }
  }
  ConfigZsurfResults();
  return(ZsurfMode);
}

//==============================================================================
/// Interpolate zsurf in timestep and stores results in CurrentZsurf[] 
/// or CurrentZsurfg[].
//==============================================================================
void JSphInOutZsurf::InterpolateZsurfTime(double timestep,bool full){
  const unsigned tpos=TimePosition;
  const double tfactor=TimeFactor;
  FindTime(timestep);
  if(tpos!=TimePosition || tfactor!=TimeFactor){
    const float tf=float(TimeFactor);
    //-Computes zsurf values for timestep on CPU.
    if(full || Cpu){
      const float *ptr0=TimesZsurf+(Nptx*TimePosition);
      const float *ptr1=TimesZsurf+(Nptx*TimePositionNext);
      for(unsigned p=0;p<Nptx;p++){
        const float v0=ptr0[p];
        CurrentZsurf[p]=(ptr1[p]-v0)*tf+v0;
      }
    }
  }
  CurrentTime=timestep;
}

//==============================================================================
/// Initialise data about times.
//==============================================================================
void JSphInOutZsurf::ResetTimes(){
  if(CurrentExternal){
    CurrentZsurf=NULL;
    CurrentZsurfg=NULL;
    CurrentExternal=false;
  }
  TimeCount=0;
  delete[] Times; Times=NULL;
  delete[] TimesZsurf; TimesZsurf=NULL;
  TimeStep=TimePre=TimeNext=TimeFactor=0;
  TimePosition=TimePositionNext=UINT_MAX;
  CurrentTime=DBL_MAX;
  delete[] CurrentZsurf; CurrentZsurf=NULL;
#ifdef _WITHGPU
  if(!Cpu){
    cudaFree(TimesZsurfg);    TimesZsurfg=NULL;
    cudaFree(CurrentZsurfg);  CurrentZsurfg=NULL;
  }
#endif
}

//==============================================================================
/// Select previous and following times with data.
//==============================================================================
void JSphInOutZsurf::FindTime(double timestep){
  if(timestep!=TimeStep || TimePosition==UINT_MAX){
    unsigned pos=(TimePosition==UINT_MAX? 0: TimePosition);
    unsigned posnext=pos;
    double tpre=Times[pos];
    double tnext=tpre;
    if(TimeCount>1){
      while(tpre>=timestep && pos>0){//-Retrocede.
        pos--;
        tpre=Times[pos];
      }
      posnext=(pos+1<TimeCount? pos+1: pos);
      tnext=Times[posnext];
      while(tnext<timestep && posnext+1<TimeCount){//-Avanza.
        posnext++;
        tnext=Times[posnext];
      }
      if(posnext-pos>1){
        pos=posnext-1;
        tpre=Times[pos];
      }
    }
    TimeStep=timestep;
    TimePosition=pos;
    TimePositionNext=posnext;
    TimePre=tpre;
    TimeNext=tnext;
    const double tdif=TimeNext-TimePre;
    TimeFactor=(tdif? (timestep-TimePre)/tdif: 0);
    if(TimeFactor<0.)TimeFactor=0.;
    if(TimeFactor>1.)TimeFactor=1.;
  }
}

//==============================================================================
/// Configures zsurf line definition for InZsurf_Calculated or InZsurf_Variable
/// with mesh-data.
//==============================================================================
void JSphInOutZsurf::ComputeZsurfLine(bool forgauge,bool forceuniform){
  //-Computes reference points (only for horizontal direction).
  if(fabs(Direction.z)>0.0001)Run_Exceptioon(fun::PrintStr("Inlet/outlet %d: Direction is invalid for Zsurf calulation. Only horizontal direction is valid.",IdZone));
  if(CSP.simulate2d){
    GaugePt0=ZonePosMin;
    GaugePtz=ZonePosMax;
    GaugePtx=GaugePt0;
  }
  else{
    if(Direction.x>=0){
      if(Direction.y<=0){
        GaugePt0=ZonePosMin;
        GaugePtx=TDouble3(ZonePosMax.x,ZonePosMax.y,ZonePosMin.z);
      }
      else{
        GaugePt0=TDouble3(ZonePosMin.x,ZonePosMax.y,ZonePosMin.z);
        GaugePtx=TDouble3(ZonePosMax.x,ZonePosMin.y,ZonePosMin.z);
      }
    }
    else{
      if(Direction.y<=0){
        GaugePt0=TDouble3(ZonePosMax.x,ZonePosMin.y,ZonePosMin.z);
        GaugePtx=TDouble3(ZonePosMin.x,ZonePosMax.y,ZonePosMin.z);
      }
      else{
        GaugePt0=TDouble3(ZonePosMax.x,ZonePosMax.y,ZonePosMin.z);
        GaugePtx=TDouble3(ZonePosMin.x,ZonePosMin.y,ZonePosMin.z);
      }
    }
    GaugePtz=TDouble3(GaugePt0.x,GaugePt0.y,ZonePosMax.z);
  }
  //-Compute centered position.
  if(forceuniform && !CSP.simulate2d){
    GaugePtz=GaugePtz-GaugePt0;
    GaugePt0=GaugePtx=(GaugePt0+GaugePtx)/2.;
    GaugePtz=GaugePtz+GaugePt0;
  }
  //-Move points to fluid.
  if(forgauge){
    const double dist=double(CSP.kernelsize);
    const tdouble3 ddir=Direction*dist; 
    GaugePt0.z=GaugePtx.z=ZsurfMin;
    GaugePt0=GaugePt0+ddir;
    GaugePtz=GaugePtz+ddir;
    GaugePtx=GaugePtx+ddir;
  }
  if(CSP.simulate2d)GaugePt0.y=GaugePtz.y=GaugePtx.y=CSP.simulate2dposy;
  //-Computes other zsurf line data.
  Dptz=CSP.dp;
  Dptx=CSP.dp;
  Nptz=unsigned((fgeo::PointsDist(GaugePtz,GaugePt0)+Dptz*0.01)/Dptz)+1;
  Nptx=unsigned((fgeo::PointsDist(GaugePtx,GaugePt0)+Dptx*0.01)/Dptx)+1;
  if(Nptx>1){
    //Log->Printf("----> i:%u  nptx:%u\n",IdZone,Nptx);
    const tdouble3 pt=GaugePt0-(fgeo::VecUnitary(GaugePtx-GaugePt0)*(Dptx/2.));
    PlaDisx=fgeo::PlaneAxisDist(pt,GaugePtx-pt,Dptx);
  }

}

//==============================================================================
/// Configures gauge to measure the zsurf. Returns true when zsurf is uniform.
//==============================================================================
bool JSphInOutZsurf::ConfigGaugeZsurf(JGaugeSystem *gaugesystem){
  ComputeZsurfLine(true,OldCode);
  const string gname=fun::PrintStr("Inlet_i%d_zsurf",IdZone);
  if(OldCode){
    GaugeSwl=gaugesystem->AddGaugeSwl(gname,0,DBL_MAX,0,GaugePt0,GaugePtz,Dptz,0);
  }
  gaugesystem->SaveVtkInitPoints();
  UniformZsurf=(Nptx==1);
  return(UniformZsurf);
}

//==============================================================================
/// Configures ZsurfResults to send.
//==============================================================================
void JSphInOutZsurf::ConfigZsurfResults(){
  ZsurfResults.npt=1;
  ZsurfResults.zsurf=&InputZsurf;
  ZsurfResults.pt=(ZonePosMin+ZonePosMax)/2.;
  ZsurfResults.vdp=TDouble3(0);
  ZsurfResults.direction=Direction;
  if(!UniformZsurf){
    ZsurfResults.npt=Nptx;
    ZsurfResults.pt=GaugePt0;
    ZsurfResults.vdp=fgeo::VecUnitary(GaugePtx-GaugePt0)*Dptx;
    if(ZsurfMode==InZsurf_Calculated){
      const tdouble3 ddir=Direction*double(CSP.kernelsize); 
      ZsurfResults.pt=GaugePt0-ddir;
      //-Nothing to do since CurrentZsurf (for CPU and GPU) and CurrentZsurfg (only for GPU) contains the zsurf results.
    }
    #ifdef _WITHGPU
      if(!Cpu && CurrentZsurfg!=NULL)cudaMemcpy(CurrentZsurf,CurrentZsurfg,sizeof(float)*Nptx,cudaMemcpyDeviceToHost);
    #endif
    ZsurfResults.zsurf=CurrentZsurf;
  }
  else if(ZsurfMode==InZsurf_Calculated){//Uniform
    ZsurfResults.npt=1;
    ZsurfResults.zsurf=&InputZsurf;
    ZsurfResults.pt=GaugePt0;
    ZsurfResults.vdp=TDouble3(CSP.kernelsize);
  }
}

//==============================================================================
/// Returns ZsurfResults with current data about zsurf.
//==============================================================================
const StZsurfResult& JSphInOutZsurf::GetZsurfResults()const{
  #ifdef _WITHGPU
    if(!Cpu && CurrentZsurfg!=NULL)cudaMemcpy(CurrentZsurf,CurrentZsurfg,sizeof(float)*Nptx,cudaMemcpyDeviceToHost);
  #endif
  return(ZsurfResults);
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphInOutZsurf::GetConfig(std::vector<std::string> &lines)const{
  const bool simulate2d=CSP.simulate2d;
  lines.push_back(fun::PrintStr("Z-Surface mode: %s (%s)",TpInZsurfModeText(ZsurfMode),(UniformZsurf? "uniform": "non-uniform")));
  if(ZsurfMode!=InZsurf_Undefined){
    if(ZsurfMode==InZsurf_Variable){
      if(InputTime && !InputTime->GetFile().empty())lines.push_back(fun::PrintStr("  Zsurf file: %s",InputTime->GetFile().c_str()));
      if(TimeCount){
        double size=(double(sizeof(double)+sizeof(float)*Nptx)*TimeCount)/(1024*1024);
        lines.push_back(fun::PrintStr("  Z-Surface times: %u  (%.1f MB)",TimeCount,size));
      }
    }
    if(InputZsurf!=FLT_MAX)  lines.push_back(fun::PrintStr("  Z-Surface value (initial): %g",InputZsurf));
    else                     lines.push_back(fun::PrintStr("  Z-Surface value (initial): Undefined"));
    if(InputZbottom!=FLT_MAX)lines.push_back(fun::PrintStr("  Z-Bottom value: %g",InputZbottom));
    else                     lines.push_back(fun::PrintStr("  Z-Bottom value (initial): Undefined"));
    if(ZsurfMode==InZsurf_Variable || ZsurfMode==InZsurf_Calculated){
      lines.push_back(fun::PrintStr("  Z-Surface Remove: %s",(RemoveZsurf? "True": "False")));
      lines.push_back(fun::PrintStr("  Z-Surface SaveVTK: %s",(SvVtkZsurf? "True": "False")));
    }
    if(ZsurfMode==InZsurf_Calculated){
      lines.push_back(fun::PrintStr("  Z-Surface Minimum: %g",ZsurfMin));
      lines.push_back(fun::PrintStr("  Z-Surface final Fit: %g (%g x Dp)",ZsurfFit,ZsurfFit/float(CSP.dp)));
    }
  }
}

//==============================================================================
/// Activates or deactivates points according to zsurf.
//==============================================================================
void JSphInOutZsurf::SetInitialPoints(unsigned npt,const tdouble3 *ptpos,byte *ptok)const{
  if(UniformZsurf){
    if(InputZsurf==FLT_MAX)memset(ptok,1,sizeof(byte)*npt);
    else for(unsigned p=0;p<npt;p++)ptok[p]=(float(ptpos[p].z)<=InputZsurf? 1: 0);
  }
  else Run_Exceptioon("Valid non-uniform zsurf is missing.");
}

//==============================================================================
/// Updates zsurf according timestep. 
/// Returns the uniform zsurf or FLT_MAX when it depends on position.
//==============================================================================
float JSphInOutZsurf::UpdateZsurf(double timestep){
  if(InputTime)InputZsurf=(float)InputTime->GetValue(timestep);  
  else if(ZsurfMode==InZsurf_Calculated){
    if(GaugeSwl){
      const bool active=GaugeSwl->GetResult().modified;
      if(active)InputZsurf=GaugeSwl->GetResult().posswl.z+ZsurfFit;
    }
  }
  return(InputZsurf);
}

//==============================================================================
/// Returns current zsurf non-uniform (according to position).
//==============================================================================
float JSphInOutZsurf::GetCurrentZsurfNonUniform(const tdouble3 &ps)const{
  const float dx=float(fgeo::PlanePoint(PlaDisx,ps));
  const unsigned cx=(dx<=0? 0: unsigned(dx));
  return(CurrentZsurf[(cx<Nptx? cx: Nptx-1)]);
}


