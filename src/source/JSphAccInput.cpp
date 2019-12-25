//HEAD_DSPH
/*
 <DUALSPHYSICS>  Copyright (c) 2019 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

/// \file JSphAccInput.cpp \brief Implements the class \ref JSphAccInput.

#include "JSphAccInput.h"
#include "JLog2.h"
#include "JXml.h"
#include "Functions.h"
#include "JReadDatafile.h"
#ifdef _WITHGPU
  #include "JSphAccInput_ker.h"
#endif

using namespace std;

//##############################################################################
//# JSphAccInputMk
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphAccInputMk::JSphAccInputMk(JLog2* log,word mkfluid,bool genabled,tfloat3 acccentre,std::string file):Log(log){
  ClassName="JSphAccInputMk";
  MkFluid=mkfluid;
  GravityEnabled=genabled;
  AccCoG=acccentre;
  File=file;
  AccTime=NULL;
  AccLin=NULL;
  AccAng=NULL;
  VelLin=NULL;
  VelAng=NULL;
  Reset();
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphAccInputMk::~JSphAccInputMk(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphAccInputMk::Reset(){
  AccSize=AccCount=0;
  delete[] AccTime;  AccTime=NULL;
  delete[] AccLin;   AccLin=NULL;
  delete[] AccAng;   AccAng=NULL;
  delete[] VelLin;   VelLin=NULL;
  delete[] VelAng;   VelAng=NULL;
  AccIndex=0;
  CurrAccLin=CurrAccAng=CurrVelLin=CurrVelAng=TDouble3(0);
  LastTimestepInput=-1;
  memset(&LastOutput,0,sizeof(StAceInput));
}

//==============================================================================
/// Resizes memory space for values.
//==============================================================================
void JSphAccInputMk::Resize(unsigned size){
  AccTime=fun::ResizeAlloc(AccTime,AccCount,size);
  AccLin=fun::ResizeAlloc(AccLin,AccCount,size);
  AccAng=fun::ResizeAlloc(AccAng,AccCount,size);
  VelLin=fun::ResizeAlloc(VelLin,AccCount,size);
  VelAng=fun::ResizeAlloc(VelAng,AccCount,size);
  AccSize=size;
}

//==============================================================================
/// Returns the allocated memory.
//==============================================================================
long long JSphAccInputMk::GetAllocMemory()const{
  long long s=0;
  if(AccTime)s+=sizeof(float)*AccSize;
  if(AccLin)s+=sizeof(tfloat3)*AccSize;
  if(AccAng)s+=sizeof(tfloat3)*AccSize;
  if(VelLin)s+=sizeof(tfloat3)*AccSize;
  if(VelAng)s+=sizeof(tfloat3)*AccSize;
  return(s);
}

//==============================================================================
/// Reads data from an external file to enable time-dependent acceleration of specific particles.
//==============================================================================
void JSphAccInputMk::LoadFile(std::string file,double tmax){
  Reset();
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  Resize(rows);
  for(unsigned r=0;r<rows;r++){
    const float atime=rdat.ReadNextFloat();     //Time.
    const tfloat3 acclin=rdat.ReadNextFloat3(); //Acc Linear.
    const tfloat3 accang=rdat.ReadNextFloat3(); //Acc Angular.
    //Log->Printf("t:%g  lin:(%g  %g  %g)  ang:(%g  %g  %g)",atime,acclin.x,acclin.y,acclin.z,accang.x,accang.y,accang.z);
    //Save the loaded time value.
    AccTime[AccCount]=atime;
    //Save the loaded gravity vector.
    AccLin[AccCount]=acclin;
    //SL: Calculate angular velocity vector based on acceleration and time data loaded
    tfloat3 currvellin=TFloat3(0.0f); //SL: New linear velocity variable
    if(AccCount>0){ //SL: Angular velocity is always zero at time zero
      const float dt=AccTime[AccCount]-AccTime[AccCount-1];
      currvellin.x=VelLin[AccCount-1].x+(AccLin[AccCount].x*dt);
      currvellin.y=VelLin[AccCount-1].y+(AccLin[AccCount].y*dt);
      currvellin.z=VelLin[AccCount-1].z+(AccLin[AccCount].z*dt);
    }
    //SL: Save the calculated angular velocity vector
    VelLin[AccCount]=currvellin;
    //Save the loaded angular velocity vector (may be zero).
    AccAng[AccCount]=accang;
    //SL: Calculate angular velocity vector based on acceleration and time data loaded
    tfloat3 currvelang=TFloat3(0.0f); //SL: New angular velocity variable
    if(AccCount>0){ //SL: Angular velocity is always zero at time zero
      const float dt=AccTime[AccCount]-AccTime[AccCount-1];
      currvelang.x=VelAng[AccCount-1].x+(AccAng[AccCount].x*dt);
      currvelang.y=VelAng[AccCount-1].y+(AccAng[AccCount].y*dt);
      currvelang.z=VelAng[AccCount-1].z+(AccAng[AccCount].z*dt);
    }
    //SL: Save the calculated angular velocity vector
    VelAng[AccCount]=currvelang;
    AccCount++;      //Increment the global line counter.
  }
  //Check that at least 2 values were given or interpolation will be impossible.
  if(AccCount<2)Run_ExceptioonFile("Cannot be less than two positions in variable acceleration file.",file);
  //Check that the final value for time is not smaller than the final simulation time.
  if(double(AccTime[AccCount-1])<tmax)Run_ExceptioonFile(fun::PrintStr("Final time (%g) is less than total simulation time in variable acceleration file.",AccTime[AccCount-1]),file);
}

//==============================================================================
/// Initialisation of object for execution.
//==============================================================================
void JSphAccInputMk::Init(double tmax){
  LoadFile(File,tmax);
}

//==============================================================================
/// Loads lines with configuration information.
//==============================================================================
void JSphAccInputMk::GetConfig(std::vector<std::string> &lines)const{
  lines.push_back(fun::PrintStr("Datafile: %s",File.c_str()));
  lines.push_back(fun::PrintStr("Global gravity: %s",(GravityEnabled? "True": "False")));
  lines.push_back(fun::PrintStr("Acc center: %g,%g,%g",AccCoG.x,AccCoG.y,AccCoG.z));
}

//=================================================================================================================
/// Returns interpolation variable acceleration values. SL: Added angular and linear velocity and set gravity flag
//=================================================================================================================
const StAceInput& JSphAccInputMk::GetAccValues(double timestep){
  if(LastTimestepInput>=0 && timestep==LastTimestepInput)return(LastOutput);
  LastTimestepInput=timestep;
  double currtime=AccTime[AccIndex];
  //Find the next nearest time value compared to the current simulation time (current value used if still appropriate).
  while((AccIndex<(AccCount-1))&&(timestep>=currtime)){
    AccIndex++;                   //Increment the index.
    currtime=AccTime[AccIndex];   //Get the next value for time.
  }
  //Not yet reached the final time value, so interpolate new values.
  if(AccIndex>0 && AccIndex<AccCount){
    const double prevtime=AccTime[AccIndex-1];    //Get the previous value for time.
    //Calculate a scaling factor for time.
    const double tfactor=(timestep-prevtime)/(currtime-prevtime);
    //Interpolate and store new value for linear accelerations. (SL: changed variable names to make sense for angular velocity use)
    tdouble3 currval=ToTDouble3(AccLin[AccIndex]);
    tdouble3 prevval=ToTDouble3(AccLin[AccIndex-1]);
    CurrAccLin.x=prevval.x+(tfactor*(currval.x-prevval.x));
    CurrAccLin.y=prevval.y+(tfactor*(currval.y-prevval.y));
    CurrAccLin.z=prevval.z+(tfactor*(currval.z-prevval.z));
    //Interpolate and store new value for angular accelerations.
    currval=ToTDouble3(AccAng[AccIndex]);
    prevval=ToTDouble3(AccAng[AccIndex-1]);
    CurrAccAng.x=prevval.x+(tfactor*(currval.x-prevval.x));
    CurrAccAng.y=prevval.y+(tfactor*(currval.y-prevval.y));
    CurrAccAng.z=prevval.z+(tfactor*(currval.z-prevval.z));
    //SL: Interpolate and store new value for linear velocity.
    currval=ToTDouble3(VelLin[AccIndex]);
    prevval=ToTDouble3(VelLin[AccIndex-1]);
    CurrVelLin.x=prevval.x+(tfactor*(currval.x-prevval.x));
    CurrVelLin.y=prevval.y+(tfactor*(currval.y-prevval.y));
    CurrVelLin.z=prevval.z+(tfactor*(currval.z-prevval.z));
    //SL: Interpolate and store new value for angular velocity.
    currval=ToTDouble3(VelAng[AccIndex]);
    prevval=ToTDouble3(VelAng[AccIndex-1]);
    CurrVelAng.x=prevval.x+(tfactor*(currval.x-prevval.x));
    CurrVelAng.y=prevval.y+(tfactor*(currval.y-prevval.y));
    CurrVelAng.z=prevval.z+(tfactor*(currval.z-prevval.z));
  }
  else{ //Reached the final time value, truncate to that value.
    const unsigned index=(AccIndex>0? AccIndex-1: 0);
    CurrAccLin=ToTDouble3(AccLin[index]);     //Get the last position for linear acceleration.
    CurrAccAng=ToTDouble3(AccAng[index]);     //Get the last position for angular acceleration.
    CurrVelLin=ToTDouble3(VelLin[index]);     //SL: Get the last position for angular velocity.
    CurrVelAng=ToTDouble3(VelAng[index]);     //SL: Get the last position for angular velocity.
  }
  //Return values.
  LastOutput.mkfluid=MkFluid;
  LastOutput.acclin=CurrAccLin;
  LastOutput.accang=CurrAccAng;
  LastOutput.centre=ToTDouble3(AccCoG);
  LastOutput.vellin=CurrVelLin; //SL: Added linear velocity
  LastOutput.velang=CurrVelAng; //SL: Added angular velocity
  LastOutput.setgravity=GravityEnabled; //SL: Added set gravity flag
  return(LastOutput);
}

//##############################################################################
//# JSphAccInput
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphAccInput::JSphAccInput(JLog2* log,const std::string &dirdata,JXml *sxml,const std::string &place):Log(log),DirData(dirdata){
  ClassName="JSphAccInput";
  Reset();
  LoadXml(sxml,place);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphAccInput::~JSphAccInput(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphAccInput::Reset(){
  for(unsigned c=0;c<Inputs.size();c++)delete Inputs[c];
  Inputs.clear();
  MemSize=0;
}

//==============================================================================
/// Returns true if mkfluid value is already configured.
//==============================================================================
bool JSphAccInput::ExistMk(word mkfluid)const{
  bool ret=false;
  for(unsigned c=0;c<Inputs.size() && !ret;c++)ret=(Inputs[c]->GetMkFluid()==mkfluid);
  return(ret);
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSphAccInput::LoadXml(JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNode(place,false);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JSphAccInput::ReadXml(const JXml *sxml,TiXmlElement* lis){
  //-Loads list of inputs.
  TiXmlElement* ele=lis->FirstChildElement("accinput"); 
  while(ele){
    word mkfluid=(word)sxml->ReadElementUnsigned(ele,"mkfluid","value");
    tfloat3 acccentre=sxml->ReadElementFloat3(ele,"acccentre");
    bool genabled=sxml->ReadElementBool(ele,"globalgravity","value");
    std::string file=DirData+sxml->ReadElementStr(ele,"datafile","value");
    if(ExistMk(mkfluid))Run_Exceptioon("An input already exists for the same mkfluid.");
    JSphAccInputMk *input=new JSphAccInputMk(Log,mkfluid,genabled,acccentre,file);
    Inputs.push_back(input);
    ele=ele->NextSiblingElement("accinput");
  }
}

//==============================================================================
/// Method to load data from input files for variable acceleration.
//==============================================================================
void JSphAccInput::Init(double tmax){
  MemSize=0;
  for(unsigned c=0;c<Inputs.size();c++){
    Inputs[c]->Init(tmax);
    MemSize+=Inputs[c]->GetAllocMemory();
  }
}

//==============================================================================
/// Shows object configuration using Log.
//==============================================================================
void JSphAccInput::VisuConfig(std::string txhead,std::string txfoot)const{
  if(!txhead.empty())Log->Print(txhead);
  for(unsigned c=0;c<GetCount();c++){
    const JSphAccInputMk* ip=Inputs[c];
    Log->Printf("Input_%u (mkfluid:%u)",c,ip->GetMkFluid());
    std::vector<std::string> lines;
    ip->GetConfig(lines);
    for(unsigned i=0;i<unsigned(lines.size());i++)Log->Print(string("  ")+lines[i]);
  }
  if(!txfoot.empty())Log->Print(txfoot);
}

//=====================================================================================================================================================
/// Returns interpolation variable acceleration values. SL: Corrected spelling mistake in exception and added angular velocity and global gravity flag
//=====================================================================================================================================================
const StAceInput& JSphAccInput::GetAccValues(unsigned cfile,double timestep){
  if(cfile>=GetCount())Run_Exceptioon("The number of input file for variable acceleration is invalid.");
  return(Inputs[cfile]->GetAccValues(timestep));
}

//==============================================================================
/// Adds variable acceleration from input configurations.
//==============================================================================
void JSphAccInput::RunCpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
  ,const typecode *code,const tdouble3 *pos,const tfloat4 *velrhop,tfloat3 *ace)
{
  for(unsigned c=0;c<GetCount();c++){
    const StAceInput v=GetAccValues(c,timestep);
    const bool withaccang=(v.accang.x!=0 || v.accang.y!=0 || v.accang.z!=0);
    const typecode codesel=typecode(v.mkfluid);
    const int ppini=int(pini),ppfin=pini+int(n);
    #ifdef OMP_USE
      #pragma omp parallel for schedule (static)
    #endif
    for(int p=ppini;p<ppfin;p++){//-Iterates through the fluid particles.
      //-Checks if the current particle is part of the particle set by its MK.
      if(CODE_GetTypeValue(code[p])==codesel){
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

#ifdef _WITHGPU
//==============================================================================
/// Adds variable acceleration from input configurations.
//==============================================================================
void JSphAccInput::RunGpu(double timestep,tfloat3 gravity,unsigned n,unsigned pini
    ,const typecode *code,const double2 *posxy,const double *posz,const float4 *velrhop,float3 *ace)
{
  for(unsigned c=0;c<GetCount();c++){
    const StAceInput v=GetAccValues(c,timestep);
    const typecode codesel=typecode(v.mkfluid);
    cuaccin::AddAccInput(n,pini,codesel,v.acclin,v.accang,v.centre,v.velang,v.vellin,v.setgravity,gravity
      ,code,posxy,posz,velrhop,ace,NULL);
  }
}
#endif

