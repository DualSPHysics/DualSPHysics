//HEAD_DSTOOLS
/* 
 <DualSPHysics codes>  Copyright (c) 2020 by Dr Jose M. Dominguez
 All rights reserved.

 DualSPHysics is an international collaboration between:
 - EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 - School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 Redistribution and use in source and binary forms, with or without modification, are permitted provided that 
 the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the distribution.
 * Neither the name of the DualSPHysics nor the names of its contributors may be used to endorse or promote products derived 
   from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, 
 BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT 
 SHALL THE COPYRIGHT OWNER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "JWaveAwasZsurf.h"
#include "JDsGaugeSystem.h"
#include "JSphMk.h"
#include "JLog2.h"
#include "JXml.h"
#include "JAppInfo.h"
#include "Functions.h"
#include "JSaveCsv2.h"
//:#include "JShapeVtk.h"

#include <cstring>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <climits>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

//##############################################################################
//# JWaveAwasZsurf
//##############################################################################
//==============================================================================
// Constructor de objeto.
//==============================================================================
JWaveAwasZsurf::JWaveAwasZsurf(word mkbound,const JXml *sxml,TiXmlElement* lis)
  :Log(AppInfo.LogPtr()),MkBound(mkbound)
{
  ClassName="JWaveAwasZsurf";
  GaugeSwl=NULL;
  FileSurf=NULL; FileData=NULL; 
  FileSurfStep=NULL; 
  Reset();
  ReadXml(sxml,lis);
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JWaveAwasZsurf::Reset(){
  Xpos0=DBL_MAX;
  Scell=0;
  StartAwas=Swl=GaugeX=GaugeXh=GaugeXdp=GaugeY=GaugeZmin=GaugeZmax=GaugeDpXml=CoefMassLimit=LimitAceFactor=0;
  CorrCoefStroke=CorrCoefPeriod=CorrPowerFunc=0;
  CorrLimit=CorrTime=CorrStart=0; CorrPositive=false;
  Ele2ndOrder=false;
  Xpos0=DBL_MAX;
  GravityZ=Depth=0;
  Zsurf=Xpiston=Upiston=0;
  NumSave=0;
  SaveData=0;
  if(FileSurf)delete FileSurf; FileSurf=NULL;
  if(FileData)delete FileData; FileData=NULL;
  if(FileSurfStep)SaveFileSurfStep();
  if(FileSurfStep)delete FileSurfStep; FileSurfStep=NULL;
  GaugeSwl=NULL;
  LimMotionAceMax=LimAceMax=0;
  LimVel=LimAcePre=LimAce=DBL_MAX;
  LimAceStats.Reset();
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JWaveAwasZsurf::ReadXml(const JXml *sxml,TiXmlElement* ele){
  StartAwas=sxml->ReadElementDouble(ele,"startawas","value",true,DBL_MAX);
  CorrCoefStroke=sxml->ReadElementDouble(ele,"correction","coefstroke",true,0);
  CorrCoefPeriod=sxml->ReadElementDouble(ele,"correction","coefperiod",true,0);
  CorrPowerFunc=sxml->ReadElementDouble(ele,"correction","powerfunc",true,1);
  if(CorrPowerFunc<1)Run_Exceptioon("Value powerfunc must not be lower than 1.");
  if(CorrCoefStroke || CorrCoefPeriod){
    if(CorrCoefStroke<1)Run_Exceptioon("CoefStroke must be higher than one.");
    if(CorrCoefPeriod<=0)Run_Exceptioon("CoefPeriod is invalid.");
  }
  Swl=sxml->ReadElementDouble(ele,"swl","value");
  word waveorder=(word)sxml->ReadElementUnsigned(ele,"elevation","value",true,2);
  Ele2ndOrder=false;
  if(waveorder==1)Ele2ndOrder=false;
  else if(waveorder==2)Ele2ndOrder=true;
  else Run_Exceptioon("Order wave to calculate elevation is invalid.");
  GaugeX=sxml->ReadElementDouble(ele,"gaugex","value",true,DBL_MAX);
  GaugeXh=sxml->ReadElementDouble(ele,"gaugex","valueh",true,DBL_MAX);
  GaugeXdp=sxml->ReadElementDouble(ele,"gaugex","valuedp",true,DBL_MAX);
  if((GaugeX!=DBL_MAX? 1: 0)+(GaugeXh!=DBL_MAX? 1: 0)+(GaugeXdp!=DBL_MAX? 1: 0)>1){
    sxml->ErrReadElement(sxml->GetFirstElement(ele,"gaugex",false),"gaugex",false,"GaugeX value must be defined only one time (using value, valueh or valuedp).");
  }
  GaugeY=sxml->ReadElementDouble(ele,"gaugey","value");
  GaugeZmin=sxml->ReadElementDouble(ele,"gaugezmin","value",true,DBL_MAX);
  GaugeZmax=sxml->ReadElementDouble(ele,"gaugezmax","value",true,-DBL_MAX);
  GaugeDpXml=sxml->ReadElementDouble(ele,"gaugedp","value",true,0.1);
  CoefMassLimit=sxml->ReadElementDouble(ele,"coefmasslimit","value",true,DBL_MAX);
  SaveData=(byte)sxml->ReadElementUnsigned(ele,"savedata","value",true,0);
  if(SaveData>3)sxml->ErrReadElement(ele,"savedata",false,"Value out of range.");
  LimitAceFactor=sxml->ReadElementDouble(ele,"limitace","value",true,2);
}

//==============================================================================
// Prepara inicio de ejecucion.
//==============================================================================
void JWaveAwasZsurf::Init(JGaugeSystem *gaugesystem,const JSphMk *mkinfo,tdouble3 pistondir
  ,double gravityz,double depth,double wavestart,double waveheight,double waveperiod
  ,double initphase,double timeramp,double wavelength,double amplitude,double maxacepiston)
{
  const unsigned cmk=mkinfo->GetMkBlockByMkBound(MkBound);
  if(cmk<mkinfo->Size()){
    //if(pistondir.z!=0 || (pistondir.x!=0 && pistondir.y!=0))Run_Exceptioon("AWAS can only be used when the piston only moves in direction X or Y.");
    if(pistondir.z!=0 || pistondir.y!=0 || pistondir.x<=0)Run_Exceptioon("AWAS can only be used when the piston only moves in direction X+.");
    Xpos0=mkinfo->Mkblock(cmk)->GetPosMax().x;
  }
  if(Xpos0==DBL_MAX)Run_Exceptioon("Initial position of piston for AWAS is unknown.");
  //-Constants for waves.
  GravityZ=gravityz;
  Depth=depth;
  //-Other constants.
  Scell=gaugesystem->GetScell();
  if(StartAwas==DBL_MAX)StartAwas=wavestart+timeramp;
  if(StartAwas<wavestart+timeramp)Run_Exceptioon("The use of AWAS must begin after initial ramp for waves generation.");
  GaugeDp=GaugeDpXml*gaugesystem->GetDp();
  if(GaugeZmin<gaugesystem->GetDomPosMin().z)GaugeZmin=gaugesystem->GetDomPosMin().z;
  if(GaugeZmax>gaugesystem->GetDomPosMax().z)GaugeZmax=gaugesystem->GetDomPosMax().z;
  unsigned gaugenmax=unsigned((GaugeZmax-GaugeZmin)/GaugeDp);
  GaugeZmax=GaugeZmin+GaugeDp*gaugenmax;
  if(CoefMassLimit==DBL_MAX)CoefMassLimit=(gaugesystem->GetSimulate2D()? 0.4: 0.5);
  MassLimit=float(CoefMassLimit*gaugesystem->GetMassFluid());
  if(GaugeX==DBL_MAX)GaugeX=gaugesystem->GetDp()*5;
  if(GaugeXdp!=DBL_MAX)GaugeX=GaugeXdp*gaugesystem->GetDp();
  if(GaugeXh!=DBL_MAX)GaugeX=GaugeXh*gaugesystem->GetKernelH();
  //-Configura drift correction.
  if(CorrCoefStroke){
    CorrLimit=CorrCoefStroke*amplitude;
    CorrTime=CorrCoefPeriod*waveperiod;
  }
  //-Initialization of variables.
  SvDataStep=false;
  Ztarget=0;
  Zsurf=GaugeZmin-Swl;
  Xpiston=Xpos0;
  Upiston=0;
  TimeStepM1=ZsurfM1=0;
  VarNum=0;
  VarZsurf=VarDt=VarZsurfDt=TDouble3(0);
  SurfStepCount=0;
  if(SaveData)CreateFileSurf(initphase==DBL_MAX,waveperiod,waveheight,wavelength,amplitude,initphase);
  if(SaveData==3)CreateFileSurfStep(0);
  //-Correction for extreme movements.
  LimMotionAceMax=maxacepiston;
//  const double pulso=TWOPI/waveperiod;
//  LimMotionAceMax=fabs(-amplitude*pulso*pulso*sin(pulso*waveperiod/4));
  if(!LimMotionAceMax)LimitAceFactor=0;
  LimAceMax=(LimitAceFactor? LimMotionAceMax*LimitAceFactor: DBL_MAX);
  //-Configura medidor de Eta.
  const tdouble3 point0=TDouble3(Xpiston+GaugeX,GaugeY,GaugeZmin);
  const tdouble3 point2=TDouble3(Xpiston+GaugeX,GaugeY,GaugeZmax);
  GaugeSwl=gaugesystem->AddGaugeSwl(fun::PrintStr("AwasMkb%02u",MkBound),StartAwas,DBL_MAX,0,point0,point2,GaugeDp,MassLimit);
}

//==============================================================================
// Introduce informacion de configuracion en lines.
//==============================================================================
void JWaveAwasZsurf::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back("AWAS-Zsurf configuration:");
  lines.push_back(fun::PrintStr("  StartAWAS: %f [s]",StartAwas));
  lines.push_back(fun::PrintStr("  Order wave theory: %s",(Ele2ndOrder? "2nd": "1st")));
  lines.push_back(fun::PrintStr("  SWL: %f [m]",Swl));
  lines.push_back(fun::PrintStr("  Piston position: %f [m]",Xpos0));
  lines.push_back(fun::PrintStr("  GaugePos: (piston + %g, %g, %g -- %g) GaugeDp:%f [m]",GaugeX,GaugeY,GaugeZmin,GaugeZmax,GaugeDp));
  lines.push_back(fun::PrintStr("  MassLimit: %f [kg] (%f*MassParticle)",MassLimit,CoefMassLimit));
  lines.push_back(fun::PrintStr("  LimitAce: %f [m/s^2] (peak theoretical: %f [m/s^2])",(LimAceMax!=DBL_MAX? LimAceMax: 0),LimMotionAceMax));
  if(!CorrCoefStroke)lines.push_back("  Correction: disabled");
  else               lines.push_back(fun::PrintStr("  Correction: enabled  (coefstroke:%g (%g [m])  coefperiod:%g (%g [s])  powerfunc:%g)",CorrCoefStroke,CorrLimit,CorrCoefPeriod,CorrTime,CorrPowerFunc));
}

//using namespace jcsv;
//==============================================================================
/// Crea fichero con informacion de free-surface.
//==============================================================================
void JWaveAwasZsurf::CreateFileSurf(bool spectrum,double waveperiod,double waveheight,double wavelength,double amplitude,double initphase){
  const string file=AppInfo.GetDirOut()+fun::PrintStr("AWASFreeSurf_mkbound_%u.csv",MkBound);
  Log->AddFileInfo(file,"Saves free surface and acceleration information about AWAS operation.");
  FileSurf=new jcsv::JSaveCsv2(file,false,AppInfo.GetCsvSepComa());
  //-Saves head.
  *FileSurf << "Order wave theory:" << (Ele2ndOrder? "2nd": "1st") << jcsv::Endl();
  *FileSurf << "SWL [m]:" << Swl << jcsv::Endl();
  *FileSurf << "Depth [m]:" << Depth << jcsv::Endl();
  if(spectrum){
    *FileSurf << "Wave period peak [s]:" << waveperiod << jcsv::Endl();
    *FileSurf << "Wave height significant [m]:" << waveheight << jcsv::Endl();
    *FileSurf << "Wave length [m]:"   << "?" << jcsv::Endl();
    *FileSurf << "Amplitude [m]:"     << "?" << jcsv::Endl();
    *FileSurf << "Initial phase:" << "?" << jcsv::Endl();
  }
  else{
    *FileSurf << "Wave period [s]:"   << waveperiod << jcsv::Endl();
    *FileSurf << "Wave height [m]:"   << waveheight << jcsv::Endl();
    *FileSurf << "Wave length [m]:"   << wavelength << jcsv::Endl();
    *FileSurf << "Amplitude [m]:"     << amplitude  << jcsv::Endl();
    *FileSurf << "Initial phase:" << initphase  << jcsv::Endl();

  }
  *FileSurf << "X0 [m]:" << Xpos0 << jcsv::Endl();
  *FileSurf << " "   << jcsv::Endl();
  *FileSurf << "Timestep [s];Xsurf [m];Zsurf [m];ZsurfTarget [m];ZminCheck [m];AceMean [m/s^2];AceMin [m/s^2];AceMax [m/s^2]" << jcsv::Endl();
}

//==============================================================================
/// Crea fichero con informacion de free-surface por step.
//==============================================================================
void JWaveAwasZsurf::CreateFileSurfStep(unsigned numfile){
  SurfStepFile=numfile;
  delete FileSurfStep;
  Log->AddFileInfo(AppInfo.GetDirOut()+"AWASFreeSurfStep_????.csv","Saves free surface and acceleration information about AWAS operation for each simulation step.");
  FileSurfStep=new jcsv::JSaveCsv2(AppInfo.GetDirOut()+fun::FileNameSec("AWASFreeSurfStep.csv",SurfStepFile),false,AppInfo.GetCsvSepComa());
  //-Saves head.
  FileSurfStep->SetHead();
  *FileSurfStep << "Timestep [s];Xsurf [m];Zsurf [m];ZsurfTarget [m];AcePre [m/s^2];Ace [m/s^2]" << jcsv::Endl();
}

//==============================================================================
/// Graba informacion de free-surface por step.
//==============================================================================
void JWaveAwasZsurf::SaveFileSurfStep(){
  FileSurfStep->SetData();
  *FileSurfStep << jcsv::Fmt(jcsv::TpDouble1,"%20.12E");
  for(unsigned c=0;c<SurfStepCount;c++){
    const tdouble4 v=SurfStepData[c];
    *FileSurfStep << v.x << v.y << v.z << v.w;
    {
      const tdouble2 w=SurfStepData2[c];
      *FileSurfStep << w.x << w.y;
    }
    *FileSurfStep << jcsv::Endl();
  }
  FileSurfStep->SaveData();
  SurfStepCount=0;
}

//==============================================================================
/// Graba informacion de free-surface por step.
//==============================================================================
void JWaveAwasZsurf::SaveFileSurfStep(double timestep,double xpos,double zsurf,double zsurftarget){
  const double tmax=(SurfStepFile+1)*5;
  if(timestep>=tmax || SvDataStep || SurfStepCount==SurfStepSize)SaveFileSurfStep();
  if(timestep>=tmax)CreateFileSurfStep(SurfStepFile+1);
  SurfStepData[SurfStepCount]=TDouble4(timestep,xpos,zsurf,zsurftarget);
  SurfStepData2[SurfStepCount]=TDouble2(LimAcePre,LimAce);
  SurfStepCount++;
}

//==============================================================================
/// Graba informacion de free-surface.
//==============================================================================
void JWaveAwasZsurf::SaveInfoSurf(){
  const tfloat3 point0=GaugeSwl->GetResult().point0;
  FileSurf->SetData();
  *FileSurf << GaugeSwl->GetResult().timestep << point0.x << Zsurf << Ztarget << point0.z;
  if(LimAceStats.GetValues())*FileSurf << LimAceStats.GetMean() << LimAceStats.GetMin() << LimAceStats.GetMax();
  else *FileSurf << "" << "";
  //:if(SaveVarInfo){
  //:  FileSurf->AddValuesf("%u;%20.12E;%20.12E;%20.12E;%20.12E;%20.12E;%20.12E;%20.12E;%20.12E;%20.12E",VarNum,VarZsurf.y,VarZsurf.z,VarZsurf.x,VarDt.y,VarDt.z,VarDt.x,VarZsurfDt.y,VarZsurfDt.z,VarZsurfDt.x);
  //:  VarNum=0;
  //:}
  *FileSurf << jcsv::Endl();
  FileSurf->SaveData();
  LimAceStats.Reset();
}

//==============================================================================
// Anota la superficie libre en fluido y la esperada para correccion de AWAS.
//==============================================================================
void JWaveAwasZsurf::RunAwas(double timestep,double ztarget,bool svdata){
  Ztarget=ztarget; Zsurf=GaugeSwl->GetResult().posswl.z-Swl; SvDataStep=svdata;
  if(SaveData){
    TimeStepM1=timestep; ZsurfM1=Zsurf;
    if(SvDataStep)SaveInfoSurf();
    if(FileSurfStep)SaveFileSurfStep(GaugeSwl->GetResult().timestep,GaugeSwl->GetResult().point0.x,Zsurf,Ztarget);
  }
}

//==============================================================================
// Devuelve la correccion del desplazamiento para el intervalo de tiempo indicado.
//==============================================================================
double JWaveAwasZsurf::GetPosCorrection(double timestep,double ztarget,bool svdata,double dt,double xmov,bool dgvoid){
  RunAwas(timestep,ztarget,svdata);
  const bool DG=false;//(timestep>0.499);
  if(SvDataStep && SaveData>=2 && !FileData){
    const string file=AppInfo.GetDirOut()+"AWASData.csv";
    Log->AddFileInfo(file,"Saves information about the operation of the AWAS.");
    FileData=new jcsv::JSaveCsv2(file,false,AppInfo.GetCsvSepComa());
    FileData->SetHead();
    *FileData << "Timestep [s];dt [s];x(t) [m];u(t) [m/s];z(t) [m];x(t+dt) [m];u(t+dt) [m/s];ztarget(t) [m];ur(t) [m/s];utarget(t) [m/s];movawas [m]" << jcsv::Endl();
    FileData->SetData();
  }
  const double utarget=xmov/dt;
  double movawas=0;
  double up=utarget;
  const double t=timestep+dt;
  if(DG){
    Log->Printf("\naz --> utarget:%f ",utarget);
  }

  if(timestep>StartAwas){//-Con correccion AWAS.
    double zr=Ztarget-Zsurf;
    double ur=zr*sqrt(GravityZ/Depth);
    up=utarget+ur;
    double xp=Xpiston+(up+Upiston)*dt/2;
    movawas=(xp-Xpiston)-xmov;
    if(DG){
      Log->Printf("az --> zr:%f ur:%f up:%f xp:%f ma:%f",zr,ur,up,xp,movawas);
    }

    if(SvDataStep && FileData){
      //double utarget2=DosPiT*Amplitude*cos(DosPiT*timestep+InitPhase);
      //FileData->AddValuesf("%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f",timestep,dt,Xpiston,Upiston,Zsurf,xp,up,Ztarget,ur,utarget,movawas,utarget2);
      *FileData << timestep << dt << Xpiston << Upiston << Zsurf << xp << up << Ztarget << ur << utarget << movawas << jcsv::Endl();
      FileData->SaveData();
    }
  }
  else if(SvDataStep && FileData){
    *FileData << timestep << dt << Xpiston << Upiston << Zsurf << (Xpiston+xmov) << utarget << jcsv::Endl();
    FileData->SaveData();
  }

  //-Applies drift correction.
  if(CorrLimit){//-Drift correction enabled.
    const double xp=(Xpiston+(xmov+movawas))-Xpos0;
    if(!CorrStart && fabs(xp)>CorrLimit){
      CorrStart=timestep;
      CorrPositive=(xp<0);
      Log->Printf("**Drift Correction> time:%f -> %f xp:%f (limit:%f)",timestep,timestep+CorrTime,xp,CorrLimit);
    }
    if(CorrStart){//-Correccion en curso.
      if(timestep-CorrStart<CorrTime){//-Aplica correccion para dt indicado.
        double corrdt=GetDriftCorrection(timestep-CorrStart+dt)-GetDriftCorrection(timestep-CorrStart);
        movawas=(CorrPositive? movawas+corrdt:  movawas-corrdt);
        //Log->Printf("==> Drift correction. t:%f corrdt:%f",timestep,(CorrPositive? corrdt:  -corrdt));
      }
      else CorrStart=0;
    }
  }

  //-Limits exteme movements.
  {
    const double dx=xmov+movawas;
    const double limvel0=LimVel,limace0=LimAce;
    LimVel=dx/dt;
    LimAcePre=LimAce=(limvel0!=DBL_MAX? (LimVel-limvel0)/dt: DBL_MAX);
    if(timestep>StartAwas && LimAce!=DBL_MAX && fabs(LimAce)>LimAceMax){//-AWAS activado y ace fuera de limites.
      const double dx2=(LimAce>=0? LimAceMax: -LimAceMax)*dt*dt + limvel0*dt;
      //Log->Printf("\n---> limace0:%f limvel0:%f dt:%f dx:%f",limace0,limvel0,dt,dx);
      movawas=dx2-xmov;
      LimVel =dx2/dt;
      LimAce =(LimVel-limvel0)/dt;
      //Log->Printf("     LimAcePre:%f LimAceMax:%f LimAce:%f dx2:%f",LimAcePre,LimAceMax,LimAce,dx2);
    }
    if(limace0!=DBL_MAX && LimAce!=DBL_MAX)LimAceStats.AddValue(LimAce);
  }

  //-Anula correccion.
  if(dgvoid)movawas=0;
  //-Actualiza posicion y velocidad del piston.
  Xpiston+=(xmov+movawas);
  Upiston=up;

  //-Actualiza posicion de gauge.
  UpdateGaugePoints();

  if(SvDataStep && SaveData)NumSave++;
  return(movawas);
}

//==============================================================================
/// Calculates measurement limits according cell size and last measurement.
/// Calcula limites de medicion segun tamanho de celdas y ultima medicion.
//==============================================================================
void JWaveAwasZsurf::UpdateGaugePoints(){
  double zmin=GaugeZmin,zmax=GaugeZmax;
  if(GaugeSwl->GetResult().modified){
    const double zsurf0=Zsurf+Swl;
    zmin=max(zsurf0-Scell,zmin);
    zmax=min(zsurf0+Scell,zmax);
  }
  //-Actualiza posicion de gauge.
  const tdouble3 point0=TDouble3(Xpiston+GaugeX,GaugeY,zmin);
  const tdouble3 point2=TDouble3(Xpiston+GaugeX,GaugeY,zmax);
  GaugeSwl->SetPoints(point0,point2,GaugeDp);
}

//==============================================================================
// Devuelve la correccion acumulada para el tiempo de correccion indicada.
// Cuando timecorr>=CorrTime entonces returns CorrLimit.
//==============================================================================
double JWaveAwasZsurf::GetDriftCorrection(double timecorr)const{
  const double ft=min(1.,timecorr/CorrTime);
  const double fc=1.-pow(1.-ft,CorrPowerFunc);
  return(CorrLimit*fc);
}
