//HEAD_DSCODES
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

/// \file JMotionObj.cpp \brief Implements the classes \ref JMotionObj and \ref JMotionMovActive.

#include "JMotionObj.h"
#include "JMotionMov.h"
#include "JMotionEvent.h"
#include "JMotion.h"
#include "Functions.h"
#include "JXml.h"
#include "JReadDatafile.h"
#include <fstream>
#include <cstring>

using namespace std;
//using std::ios;

//##############################################################################
//# JMotionMovActive
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JMotionMovActive::JMotionMovActive(double start,double eventfinish,JMotionMov* mov):EventFinish(eventfinish){
  ClassName="JMotionMovActive";
  DfTimes=NULL; DfPos=NULL; DfAng=NULL;
  Mov=mov;
  Start=start;
  Flash=(Mov->Time<0);
  Finish=(Flash? Start: Start+Mov->Time);
  if(eventfinish>=0&&eventfinish<Finish)Finish=eventfinish;
  ConfigData();
  Del=false;
}
//==============================================================================
// Destructor.
//==============================================================================
JMotionMovActive::~JMotionMovActive(){
  DestructorActive=true;
  DfReset();
}

//==============================================================================
// Resetea memoria asignada y datos de movimientos RectFile y RotFile.
//==============================================================================
void JMotionMovActive::DfReset(){
  DfTimes=NULL;
  DfPos=NULL;
  DfAng=NULL;
  DfPosType=false;
  DfCount=0; DfIndex=0;
  DfLastPos=TDouble3(0);
  DfLastAng=0;
}

//==============================================================================
// Configura datos del mov activo
//==============================================================================
void JMotionMovActive::ConfigData(){
  DfReset();
  Vel=TDouble3(0);  VelAng=0;  Phase=TDouble3(0);  PhaseUni=0;
  switch(Mov->Type){
    case JMotionMov::Rectilinear:           Vel=((JMotionMovRect*)Mov)->Vel;            break;
    case JMotionMov::RectilinearAce:        Vel=((JMotionMovRectAce*)Mov)->Vel;         break;
    case JMotionMov::Rotation:              VelAng=((JMotionMovRot*)Mov)->VelAng;       break;
    case JMotionMov::RotationAce:           VelAng=((JMotionMovRotAce*)Mov)->VelAng;    break;
    case JMotionMov::Circular:              VelAng=((JMotionMovCir*)Mov)->VelAng;       break;
    case JMotionMov::CircularAce:           VelAng=((JMotionMovCirAce*)Mov)->VelAng;    break;
    case JMotionMov::RectilinearSinusoidal: Phase=((JMotionMovRectSinu*)Mov)->Phase;    break;
    case JMotionMov::RotationSinusoidal:    PhaseUni=((JMotionMovRotSinu*)Mov)->Phase;  break;
    case JMotionMov::CircularSinusoidal:    PhaseUni=((JMotionMovCirSinu*)Mov)->Phase;  break;
    case JMotionMov::RectilinearFile:       DfConfig(true);                             break;
    case JMotionMov::RotationFile:          DfConfig(false);                            break;
  }
}

//==============================================================================
// Carga y configura datos del movimiento a partir de fichero de datos.
//==============================================================================
void JMotionMovActive::DfConfig(bool postype){
  DfReset();
  DfPosType=postype;
  if(DfPosType){
    JMotionMovRectFile* mv=(JMotionMovRectFile*)Mov;
    mv->PrepareData();
    DfCount=mv->GetCount();
    DfTimes=mv->GetTimes();
    DfPos=mv->GetValuesPos();
    DfLastPos=DfPos[0];
  }
  else{
    JMotionMovRotFile* mv=(JMotionMovRotFile*)Mov;
    mv->PrepareData();
    DfCount=mv->GetCount();
    DfTimes=mv->GetTimes();
    DfAng=mv->GetValuesAng();
    DfLastAng=DfAng[0];
  }
}

//==============================================================================
// Devuelve posicion de times[] mas proxima por debajo de t o la ultima cuando
// t es mayor que el ultimo valor.
//==============================================================================
unsigned JMotionMovActive::BinarySearch(unsigned size,const double *times,double t){
  //unsigned n=0;
  unsigned ret=0;
  if(size>1){
    int ccen,cmin=0,cmax=int(size)-1;
    //printf("bs> cmin:%d cmax:%d\n",cmin,cmax);
    while(cmin<=cmax){
      ccen=((cmax-cmin)/2)+cmin;
      if(times[ccen]==t)cmin=cmax+1;
      else if(t<times[ccen])cmax=ccen-1;
      else cmin=ccen+1;
      //printf("bs> ccen:%d time:%f cmin:%d cmax:%d\n",ccen,times[ccen],cmin,cmax);
      //n++;
    }
    //printf("bs> ccen_ini:%d t:%f\n",ccen,t);
    if(ccen && ccen<int(size) && t<times[ccen])ccen--;
    //if(t>times[ccen])ccen--;
    //printf("bs> ccen_fin:%d t:%f\n",ccen,t);
    ret=unsigned(ccen);
    while(ret && t==times[ret-1])ret--;
  }
  //printf("bs> ret:%u\n",ret);
  //printf("bs> t:%f n:%u \n",t,n);
  return(ret);
}

//==============================================================================
// Devuelve la siguiente posicion en funcion del instante indicado
//==============================================================================
tdouble3 JMotionMovActive::DfGetNewPos(double t){
  tdouble3 newpos;
  if(DfIndex==0)DfIndex=BinarySearch(DfCount,DfTimes,t);
  while(DfIndex<DfCount&&t>DfTimes[DfIndex])DfIndex++;
  if(DfIndex>=DfCount)newpos=DfPos[DfCount-1];//-Mayor al instante final se queda con la ultima posicion.
  else{
    const double tfactor=(t-DfTimes[DfIndex-1])/(DfTimes[DfIndex]-DfTimes[DfIndex-1]);
    tdouble3 pos0=DfPos[DfIndex-1],pos=DfPos[DfIndex];
    double x=(((JMotionMovRectFile*)Mov)->FieldX>=0? pos0.x+tfactor*(pos.x-pos0.x): 0);
    double y=(((JMotionMovRectFile*)Mov)->FieldY>=0? pos0.y+tfactor*(pos.y-pos0.y): 0);
    double z=(((JMotionMovRectFile*)Mov)->FieldZ>=0? pos0.z+tfactor*(pos.z-pos0.z): 0);
    newpos=TDouble3(x,y,z);
    //newpos=DfPosX[DfIndex-1]+tfactor*(DfPosX[DfIndex]-DfPosX[DfIndex-1]);
  }
  //printf("index:%u  t:%g  newpos:%f\n",DfIndex,t,newpos);
  return(newpos);
}

//==============================================================================
// Devuelve la siguiente angulo en funcion del instante indicado
//==============================================================================
double JMotionMovActive::DfGetNewAng(double t){
  double newang;
  if(DfIndex==0)DfIndex=BinarySearch(DfCount,DfTimes,t);
  while(DfIndex<DfCount&&t>DfTimes[DfIndex])DfIndex++;
  if(DfIndex>=DfCount)newang=DfAng[DfCount-1];//-Mayor al instante final se queda con el ultimo angulo.
  else{
    const double tfactor=(t-DfTimes[DfIndex-1])/(DfTimes[DfIndex]-DfTimes[DfIndex-1]);
    double ang0=DfAng[DfIndex-1],ang=DfAng[DfIndex];
    newang=ang0+tfactor*(ang-ang0);
  }
  //printf("index:%u  t:%g  newpos:%f\n",DfIndex,t,newpos);
  return(newang);
}

//==============================================================================
// Pasa al siguiente movimiento enlazado
//==============================================================================
void JMotionMovActive::NextMov(){
  if(Mov->NextMov){
    if(!Flash)Start+=Mov->Time;
    Mov=Mov->NextMov;
    Flash=(Mov->Time<0);
    Finish=(Flash? Start: Start+Mov->Time);
    if(EventFinish>=0&&EventFinish<Finish)Finish=EventFinish;
    //-Uso de velocidad previa...
    tdouble3 velp=Vel,phasep=Phase;
    double velangp=VelAng,phaseunip=PhaseUni;
    ConfigData();
    switch(Mov->Type){
      case JMotionMov::RectilinearAce: if(((JMotionMovRectAce*)Mov)->VelPrev)Vel=velp;       break;
      case JMotionMov::RotationAce:    if(((JMotionMovRotAce*)Mov)->VelPrev)VelAng=velangp;  break;
      case JMotionMov::CircularAce:    if(((JMotionMovCirAce*)Mov)->VelPrev)VelAng=velangp;  break;
      case JMotionMov::RectilinearSinusoidal: if(((JMotionMovRectSinu*)Mov)->PhasePrev)Phase=phasep;       break;
      case JMotionMov::RotationSinusoidal:    if(((JMotionMovRotSinu*)Mov)->PhasePrev)PhaseUni=phaseunip;  break;
      case JMotionMov::CircularSinusoidal:    if(((JMotionMovCirSinu*)Mov)->PhasePrev)PhaseUni=phaseunip;  break;
    }
  }
}


//##############################################################################
//# JMotionObj
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JMotionObj::JMotionObj(unsigned id,JMotionObj* parent,int ref):Id(id),Parent(parent),Ref(ref){
  ClassName="JMotionObj";
  //printf("<New-Obj:%d>\n",ref);
  Reset();
}

//==============================================================================
// Destructor.
//==============================================================================
JMotionObj::~JMotionObj(){
  DestructorActive=true;
  //printf("<DEL-obj:%d>\n",Ref);
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JMotionObj::Reset(){
  //-Libera objetos hijos
  for(unsigned c=0;c<Children.size();c++)delete Children[c];
  Children.clear();
  //-Libera movimientos
  for(unsigned c=0;c<Movs.size();c++)delete Movs[c];
  Movs.clear();
  //-Libera ejes
  for(unsigned c=0;c<Axis.size();c++)delete Axis[c];
  Axis.clear();
  //-Libera eventos
  for(unsigned c=0;c<Events.size();c++)delete Events[c];
  Events.clear();
  //-Libera movimientos activos
  for(unsigned c=0;c<ActiveMovs.size();c++)delete ActiveMovs[c];
  ActiveMovs.clear();
  Active=false;
  Moving=false;
  //Pos.Reset();
  ModPos.Reset();
}

//==============================================================================
// Devuelve puntero al objeto con el id indicado
//==============================================================================
JMotionObj* JMotionObj::ObjGetPointer(unsigned id){
  JMotionObj *obj=(Id==id? this: NULL);
  for(unsigned c=0;c<Children.size()&&!obj;c++)obj=Children[c]->ObjGetPointer(id);
  return(obj);
}

//==============================================================================
// Devuelve puntero al objeto con la referencia indicada
//==============================================================================
JMotionObj* JMotionObj::ObjGetPointerByRef(int ref){
  JMotionObj *obj=(Ref==ref? this: NULL);
  for(unsigned c=0;c<Children.size()&&!obj;c++)obj=Children[c]->ObjGetPointerByRef(ref);
  return(obj);
}

//==============================================================================
// Devuelve puntero al objeto solicitado
//==============================================================================
JMotionMov* JMotionObj::MovGetPointer(unsigned id)const{
  JMotionMov *mov=NULL;
  for(unsigned c=0;c<Movs.size()&&!mov;c++)if(Movs[c]->Id==id)mov=Movs[c];
  return(mov);
}

//==============================================================================
// Devuelve puntero al eje solicitado o null sino existe
//==============================================================================
JMotionAxis* JMotionObj::AxisGetPointer(const tdouble3 &p1,const tdouble3 &p2)const{
  JMotionAxis *axis=NULL;
  for(unsigned c=0;c<Axis.size()&&!axis;c++)if(Axis[c]->Equals(p1,p2))axis=Axis[c];
  return(axis);
}

//==============================================================================
// Devuelve la posicion del movimiento solicitado.
//==============================================================================
int JMotionObj::GetPosMov(JMotionMov* mv)const{
  int pos=-1;
  for(int c=0;c<int(Movs.size())&&pos<0;c++)if(Movs[c]==mv)pos=c;
  if(pos<0)Run_Exceptioon("Cannot find the requested movement.");
  return(pos);
}

//==============================================================================
// Anhade un objeto hijo, movimiento, evento o eje.
//==============================================================================
void JMotionObj::AddChild(JMotionObj* obj){   Children.push_back(obj); }
void JMotionObj::AddMov(JMotionMov* mov){     Movs.push_back(mov);     }
void JMotionObj::AddAxis(JMotionAxis* axis){  Axis.push_back(axis);    }
void JMotionObj::AddEvent(JMotionEvent* evt){ Events.push_back(evt);   }

//==============================================================================
// Comprueba y establece enlaces entre movimientos.
//==============================================================================
void JMotionObj::LinkMovs(){
  for(unsigned c=0;c<Movs.size();c++){
    JMotionMov* mov=Movs[c];
    if(mov->NextId){
      mov->SetNextMov(MovGetPointer(mov->NextId));
      if(!mov->NextMov)Run_Exceptioon(fun::PrintStr("The movement with id=%u refers to another non-existent movement (id=%u) within the object.",mov->Id,mov->NextId)); 
    }
  }
  for(unsigned c=0;c<Children.size();c++)Children[c]->LinkMovs();
}

//==============================================================================
// Devuelve el numero total de objetos hijos (todos los niveles)
//==============================================================================
unsigned JMotionObj::ChildrenCount(){
  unsigned n=0;
  for(unsigned c=0;c<Children.size();c++)n+=1+Children[c]->ChildrenCount();
  return(n);
}

//==============================================================================
// Coloca movimiento de evento en lista de movs activos. 
// Devuelve true si paso ahora a estar activo.
//==============================================================================
void JMotionObj::BeginEvent(double start,double eventfinish,JMotionMov* mov){
  JMotionMovActive* amov=new JMotionMovActive(start,eventfinish,mov);
  if(!amov)Run_Exceptioon("Cannot allocate the requested memory.");
  ActiveMovs.push_back(amov);
  //-Pasa a ser un objeto activo el y sus antecesores
  JMotionObj* obj=this;
  while(obj&&!obj->Active){ obj->Active=true; obj=(JMotionObj*)obj->Parent; }
}

//==============================================================================
// Reinicia ejecucion al timestep 0.
//==============================================================================
void JMotionObj::ResetTime(){
  //-Libera movimientos activos
  for(unsigned c=0;c<ActiveMovs.size();c++)delete ActiveMovs[c];
  ActiveMovs.clear();
  Active=false;
  Moving=false;
  //Pos.Reset();
  ModPos.Reset();
  //-Aplica reinicio a los hijos.
  for(unsigned c=0;c<Children.size();c++)Children[c]->ResetTime();
  //-Aplica reinicio a Axis[].
  for(unsigned c=0;c<Axis.size();c++)Axis[c]->ResetTime();
}

//==============================================================================
// Calcula desplazamiento de objeto
// Devuelve true si el objeto o alguno de sus hijos esta activo.
//==============================================================================
bool JMotionObj::ProcesTime(double timestep,double dt,JMotionObj** lismov,unsigned &lismovcount,JMotionObj** lisstop,unsigned &lisstopcount){
  //printf("ProcesTime-> %f %f  \n",timestep,dt);

//printf("\nProcesTime>\n");  //printf("timestep: %G    dt: %G\n",timestep,dt);
  Active=true;
  bool modif=false;

  if(Parent&&Parent->Moving){//-Aplica movimiento del padre a todos los ejes
    //printf("ProcesTime-> Parent->Moving\n");
    int n=int(Axis.size());
    for(int c=0;c<n;c++)Parent->ModPos.PointsMove(Axis[c]->P1,Axis[c]->P2);
    /*for(int c=0;c<n;c++){
      //printf("Ref[%d]>> Axisp1:(%g,%g,%g)-",Ref,Axis[c]->P1.x,Axis[c]->P1.y,Axis[c]->P1.z);
      Parent->ModPos.PointsMove(Axis[c]->P1,Axis[c]->P2);
      //Axis[c]->P1=Parent->ModPos.PointMove(Axis[c]->P1);
      //Axis[c]->P2=Parent->ModPos.PointMove(Axis[c]->P2);
      //printf("(%g,%g,%g)\n",Axis[c]->P1.x,Axis[c]->P1.y,Axis[c]->P1.z);
      //Parent->ModPos.GetMatrix().Print("ModAxis");
    }*/
  }

  int na=int(ActiveMovs.size());
  //printf("ProcesTime-> ActiveMovs.size():%d\n",ActiveMovs.size());
  if(na){
    double tstepfin=timestep+dt;
    ModPos.Reset();
//ModPos.ToMatrix();
    for(int ca=0;ca<na;ca++){
      JMotionMovActive* amov=ActiveMovs[ca];
      if(amov->Del){//-Elimina mov-activo marcado para borrar.
//printf("*** Eliminacion de movimiento activo timestep: %G   dt: %G\n",timestep,dt);
        ActiveMovs.erase(ActiveMovs.begin()+ca);
        delete amov;
        ca--; na--;
      }
      else{//-Procesa mov-activo.
        bool rep;
        double dt2=dt;
        double timestep2=timestep;
        do{
//printf("timestep2: %G    dt2: %G\n",timestep2,dt2);
          rep=false;
          JMotionMov* mov=amov->Mov;
//printf("mstart: %G    mfinish: %G\n",amov->Start,amov->Finish);
          //-Calcula dt ajustado al movimiento y el dt sobrante para el siguiente movimiento si lo hubiera.
          double dtmov=(tstepfin>amov->Finish? amov->Finish-timestep2: dt2);
          double dtover=dt2-dtmov;
//printf("dtmov: %G    dtover: %G\n",dtmov,dtover);
          //-Ajuste del dtmov por inicio de movimiento.
          if(timestep2<amov->Start){
            dtmov-=(amov->Start-timestep2);
//printf("ReajusteDeInicio dtmov: %G     dtmov2: %G\n",-(amov->Start-timestep2),dtmov);
          }
          //-Calcula movimiento.
          if(dtmov>0||amov->Flash)switch(mov->Type){
            case JMotionMov::Rectilinear:{
              const JMotionMovRect *mv=(JMotionMovRect*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              ModPos.Move(TDouble3(mv->Vel.x*t,mv->Vel.y*t,mv->Vel.z*t));
              modif=true;
//            printf("***JMotionMov::Rectilinear\n");
            }break;
            case JMotionMov::RectilinearAce:{
              const JMotionMovRectAce *mv=(JMotionMovRectAce*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              tdouble3 at=TDouble3(mv->Ace.x*t,mv->Ace.y*t,mv->Ace.z*t);
              ModPos.Move(TDouble3(amov->Vel.x*t+0.5f*at.x*t,amov->Vel.y*t+0.5f*at.y*t,amov->Vel.z*t+0.5f*at.z*t));
              amov->Vel=TDouble3(amov->Vel.x+at.x,amov->Vel.y+at.y,amov->Vel.z+at.z);
              modif=true;
            }break;
            case JMotionMov::Rotation:{
              const JMotionMovRot *mv=(JMotionMovRot*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              //printf("ProcesTime-> Rotation> t:%f VelAng:%f \n",t,mv->VelAng);
              //printf("ProcesTime-> AxisP1:(%f,%f,%f) AxisP2:(%f,%f,%f) \n",mv->Axis->P1.x,mv->Axis->P1.y,mv->Axis->P1.z,mv->Axis->P2.x,mv->Axis->P2.y,mv->Axis->P2.z);
              ModPos.Rotate(mv->VelAng*t,mv->Axis->P1,mv->Axis->P2);
              modif=true;
            }break;
            case JMotionMov::RotationAce:{
              const JMotionMovRotAce *mv=(JMotionMovRotAce*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              double at=mv->AceAng*t;
              ModPos.Rotate(amov->VelAng*t+0.5f*at*t,mv->Axis->P1,mv->Axis->P2);
              amov->VelAng+=at;
              modif=true;
            }break;
            case JMotionMov::Circular:{
              const JMotionMovCir *mv=(JMotionMovCir*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              JMatrix4d m=JMatrix4d::MatrixRot(mv->VelAng*t,mv->Axis->P1,mv->Axis->P2);
              tdouble3 &p1=*(const_cast<tdouble3*>(&mv->Ref->P1));
              tdouble3 &p2=*(const_cast<tdouble3*>(&mv->Ref->P2));
              p2=m.MulPoint(p1);
              ModPos.Move(TDouble3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z));
              p1=p2;
              modif=true;
            }break;
            case JMotionMov::CircularAce:{
              const JMotionMovCirAce *mv=(JMotionMovCirAce*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              double at=mv->AceAng*t;
              JMatrix4d m=JMatrix4d::MatrixRot(amov->VelAng*t+0.5f*at*t,mv->Axis->P1,mv->Axis->P2);
              tdouble3 &p1=*(const_cast<tdouble3*>(&mv->Ref->P1));
              tdouble3 &p2=*(const_cast<tdouble3*>(&mv->Ref->P2));
              p2=m.MulPoint(p1);
              ModPos.Move(TDouble3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z));
              p1=p2;
              amov->VelAng+=at;
              modif=true;
            }break;
            case JMotionMov::RectilinearSinusoidal:{
              const JMotionMovRectSinu *mv=(JMotionMovRectSinu*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              tdouble3 ph=amov->Phase;
              tdouble3 p1=TDouble3(0),p2=TDouble3(0);
              if(mv->Ampl.x){
                p1.x=mv->Ampl.x*sin(ph.x); ph.x+=double(mv->Freq.x*TWOPI*t); p2.x=mv->Ampl.x*sin(ph.x);
              }
              if(mv->Ampl.y){
                p1.y=mv->Ampl.y*sin(ph.y); ph.y+=double(mv->Freq.y*TWOPI*t); p2.y=mv->Ampl.y*sin(ph.y);
              }
              if(mv->Ampl.z){
                p1.z=mv->Ampl.z*sin(ph.z); ph.z+=double(mv->Freq.z*TWOPI*t); p2.z=mv->Ampl.z*sin(ph.z);
              }
              ModPos.Move(TDouble3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z));
              amov->Phase=ph;
              //ModPos.Move(TDouble3(0,0.1f,0));
              modif=true;
            }break;
            case JMotionMov::RotationSinusoidal:{
              const JMotionMovRotSinu *mv=(JMotionMovRotSinu*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              double ph=amov->PhaseUni;
              double ang=mv->Ampl*sin(ph); 
              ph+=double(mv->Freq*(PI+PI)*t); 
              ang=mv->Ampl*sin(ph)-ang;
              ModPos.Rotate(ang,mv->Axis->P1,mv->Axis->P2);
              amov->PhaseUni=ph;
              modif=true;
            }break;
            case JMotionMov::CircularSinusoidal:{
              const JMotionMovCirSinu *mv=(JMotionMovCirSinu*)mov;
              double t=(amov->Flash? -mv->Time: dtmov);
              double ph=amov->PhaseUni;
              double ang=mv->Ampl*sin(ph); 
              ph+=double(mv->Freq*(PI+PI)*t); 
              ang=mv->Ampl*sin(ph)-ang;
              JMatrix4d m=JMatrix4d::MatrixRot(ang,mv->Axis->P1,mv->Axis->P2);
              tdouble3 &p1=*(const_cast<tdouble3*>(&mv->Ref->P1));
              tdouble3 &p2=*(const_cast<tdouble3*>(&mv->Ref->P2));
              p2=m.MulPoint(p1);
              ModPos.Move(TDouble3(p2.x-p1.x,p2.y-p1.y,p2.z-p1.z));
              p1=p2;
              amov->PhaseUni=ph;
              modif=true;
            }break;
            case JMotionMov::RectilinearFile:{
              const JMotionMovRectFile *mv=(JMotionMovRectFile*)mov;
              double t=timestep-amov->Start;
              if(t<0)t=0;
              t+=amov->DfTimes[0];
              t+=(amov->Flash? -mv->Time: dtmov);
              tdouble3 newpos=amov->DfGetNewPos(t);
              ModPos.Move(TDouble3(newpos.x-amov->DfLastPos.x,newpos.y-amov->DfLastPos.y,newpos.z-amov->DfLastPos.z));
              amov->DfLastPos=newpos;
              modif=true;
            }break;
            case JMotionMov::RotationFile:{
              const JMotionMovRotFile *mv=(JMotionMovRotFile*)mov;
              double t=timestep-amov->Start;
              if(t<0)t=0;
              t+=amov->DfTimes[0];
              t+=(amov->Flash? -mv->Time: dtmov);
              double newang=amov->DfGetNewAng(t);
              double ang=newang-amov->DfLastAng;
              ModPos.Rotate(newang-amov->DfLastAng,mv->Axis->P1,mv->Axis->P2);
              amov->DfLastAng=newang;
              modif=true;
              //printf(" PT>> t:%f ang:%f newang:%f\n",t,ang,newang);
            }break;
          }
          //-Cambia al movimiento enlazado con el actual para terminar de consumir el dt.
          if((dtover>0||dtmov==0)&&mov->NextMov!=NULL){
            amov->NextMov();
            dt2=dtover;
            timestep2=amov->Start;
            if(timestep2<=amov->Finish||amov->Flash)rep=true;
          }
        }while(rep);
        if(tstepfin>amov->Finish)amov->Del=true;//-Lo marca para que sea borrado en el siguiente ProcesTime().
      }
    }
  }
  if(Parent&&Parent->Moving){//-Aplica movimiento del padre
    if(modif)ModPos.MoveMix(Parent->ModPos); else ModPos=Parent->ModPos;
    modif=true;
  }
  int nc=int(Children.size());
  if(modif){
    if(Ref>=0){ lismov[lismovcount]=this; lismovcount++; }
    //if(nc)Pos.MoveMix(ModPos); //-Solo mantiene el desplazamiento acumulado cuando tiene hijos q puedan usarlo.
    Moving=true;
  }
  else if(Moving){
    if(Ref>=0){ lisstop[lisstopcount]=this; lisstopcount++; }
    Moving=false;
    //printf("Motion> Obj[%d] goes to inactive...t_tmax:%f %f\n",Id,timestep,MaxTimestep);
  }
  else if(!na)Active=false; //-No tiene ningun movimiento activo y moving ya es false.
  for(int c=0;c<nc;c++)Active|=Children[c]->ProcesTime(timestep,dt,lismov,lismovcount,lisstop,lisstopcount);
  return(Active);
}

//==============================================================================
// Devuelve datos de movimiento.
//==============================================================================
bool JMotionObj::GetMov(unsigned &ref,tdouble3 &mvsimple,JMatrix4d &mvmatrix)const{
  ref=Ref;
  bool simple=ModPos.IsSimple();
  if(simple)mvsimple=ModPos.GetSimple();
  else mvmatrix=ModPos.GetMatrix();
  return(simple);
}

//==============================================================================
// Devuelve la mayor referencia de objeto definida.
//==============================================================================
int JMotionObj::GetMaxRef()const{
  int ref=Ref;
  for(unsigned c=0;c<Children.size();c++){
    int rf=Children[c]->GetMaxRef();
    if(rf>ref)ref=rf;
  }
  return(ref);
}

//==============================================================================
// Loads list of references.
//==============================================================================
void JMotionObj::GetRefs(std::vector<int> &refs)const{
  if(Ref>=0){
    unsigned c=0;
    for(;c<unsigned(refs.size()) && refs[c]!=Ref;c++);
    if(c>=unsigned(refs.size()))refs.push_back(Ref);
  }
  for(unsigned c=0;c<Children.size();c++)Children[c]->GetRefs(refs);
}

//==============================================================================
// Copia la configuracion de los movimientos a mot.
//==============================================================================
void JMotionObj::CopyConfigMovs(JMotion &mot)const{
  for(int c=0;c<int(Movs.size());c++){
    switch(Movs[c]->Type){
      case JMotionMov::Wait:{
        JMotionMovWait *mv=(JMotionMovWait*)Movs[c];
        mot.MovAddWait(Id,mv->Id,mv->NextId,mv->Time);
      }break; 
      case JMotionMov::Rectilinear:{
        JMotionMovRect *mv=(JMotionMovRect*)Movs[c];
        mot.MovAddRectilinear(Id,mv->Id,mv->NextId,mv->Time,mv->Vel);
      }break; 
      case JMotionMov::RectilinearAce:{
        JMotionMovRectAce *mv=(JMotionMovRectAce*)Movs[c];
        mot.MovAddRectilinearAce(Id,mv->Id,mv->NextId,mv->Time,mv->Ace,mv->Vel,mv->VelPrev);
      }break; 
      case JMotionMov::Rotation:{
        JMotionMovRot *mv=(JMotionMovRot*)Movs[c];
        mot.MovAddRotation(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Axis->P1,mv->Axis->P2,mv->VelAng,false);
      }break;
      case JMotionMov::RotationAce:{
        JMotionMovRotAce *mv=(JMotionMovRotAce*)Movs[c];
        mot.MovAddRotationAce(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Axis->P1,mv->Axis->P2,mv->AceAng,mv->VelAng,mv->VelPrev,false);
      }break;
      case JMotionMov::Circular:{
        JMotionMovCir *mv=(JMotionMovCir*)Movs[c];
        mot.MovAddCircular(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Axis->P1,mv->Axis->P2,mv->Ref->P1,mv->VelAng,false);
      }break;
      case JMotionMov::CircularAce:{
        JMotionMovCirAce *mv=(JMotionMovCirAce*)Movs[c];
        mot.MovAddCircularAce(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Axis->P1,mv->Axis->P2,mv->Ref->P1,mv->AceAng,mv->VelAng,mv->VelPrev,false);
      }break;
      case JMotionMov::RectilinearSinusoidal:{
        JMotionMovRectSinu *mv=(JMotionMovRectSinu*)Movs[c];
        mot.MovAddRecSinu(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Freq,mv->Ampl,mv->Phase,mv->PhasePrev,false);
      }break;
      case JMotionMov::RotationSinusoidal:{
        JMotionMovRotSinu *mv=(JMotionMovRotSinu*)Movs[c];
        mot.MovAddRotSinu(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Axis->P1,mv->Axis->P2,mv->Freq,mv->Ampl,mv->Phase,mv->PhasePrev,false);
      }break;
      case JMotionMov::CircularSinusoidal:{
        JMotionMovCirSinu *mv=(JMotionMovCirSinu*)Movs[c];
        mot.MovAddCirSinu(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Axis->P1,mv->Axis->P2,mv->Ref->P1,mv->Freq,mv->Ampl,mv->Phase,mv->PhasePrev,false);
      }break;
      case JMotionMov::RectilinearFile:{
        JMotionMovRectFile *mv=(JMotionMovRectFile*)Movs[c];
        mot.MovAddRectilinearFile(Id,mv->Id,mv->NextId,mv->Time,mv->File,mv->Fields,mv->FieldTime,mv->FieldX,mv->FieldY,mv->FieldZ);
      }break;
      case JMotionMov::RotationFile:{
        JMotionMovRotFile *mv=(JMotionMovRotFile*)Movs[c];
        mot.MovAddRotationFile(Id,mv->Id,mv->NextId,mv->Time,mv->AngDegrees,mv->Axis->P1,mv->Axis->P2,mv->File);
      }break;
      case JMotionMov::Nulo:{
        JMotionMovNull *mv=(JMotionMovNull*)Movs[c];
        mot.MovAddNull(Id,mv->Id);
      }break;
      default: Run_Exceptioon("Unrecognised movement type.");
    }
  }   
}

//==============================================================================
// Copia los datos de configuracion (Objs,Movs,Events) a mot.
//==============================================================================
void JMotionObj::CopyConfig(JMotion &mot)const{
  mot.ObjAdd(Id,(Parent? Parent->Id: 0),Ref);
  CopyConfigMovs(mot);
  for(int c=0;c<int(Children.size());c++)Children[c]->CopyConfig(mot);
}

//==============================================================================
// Copia configuracion a mot intercambiando las referencias por su nuevo valor.
// Solo copia los datos de configuracion (Objs,Movs,Events).
// Las referencias que no aparezcan pasan a ser -1.
//==============================================================================
void JMotionObj::CopyChangeRef(JMotion &mot,const int* ref,const int* refnew,unsigned refcount)const{
  int ref2=-1;
  for(unsigned c=0;c<refcount&&ref2<0;c++)if(Ref==ref[c])ref2=refnew[c];
  mot.ObjAdd(Id,(Parent? Parent->Id: 0),ref2);
  CopyConfigMovs(mot);
  for(int c=0;c<int(Children.size());c++)Children[c]->CopyChangeRef(mot,ref,refnew,refcount);
}

//==============================================================================
// Devuelve true si existe el objeto solicitado
//==============================================================================
bool JMotionObj::ExistsObj(JMotionObj* obj)const{
  bool ret=(this==obj);
  for(unsigned c=0;c<Children.size()&&!ret;c++)ret=Children[c]->ExistsObj(obj);
  return(ret);
}

//==============================================================================
// Optimiza configuracion eliminando objetos, movimientos y eventos inutiles.
// Devuelve true cuando el objeto es inutil y debe borrarse.
//==============================================================================
bool JMotionObj::Optimize(){
  //-Elimina movimientos inutiles
  if(Movs.size()){
    int nm=int(Movs.size());
    byte *use=new byte[nm];
    for(int c=0;c<nm;c++)use[c]=0;
    //-Marca movimientos referenciados
    for(int c=0;c<int(Events.size());c++)use[GetPosMov(Events[c]->Mov)]=1;
    bool run;
    do{
      run=false;
      for(int c=0;c<nm;c++)if(use[c]==1){
        use[c]=2;
        if(Movs[c]->NextMov){
          int cm=GetPosMov(Movs[c]->NextMov);
          if(!use[cm]){ use[cm]=1; run=true; }
        }
      }
    }while(run);
    //-Elimina movimientos no referenciados
    for(int c=nm-1;c>=0;c--)if(!use[c]){
      //printf("Borrando mov(o_%d,m_%d)\n",Id,Movs[c]->Id);
      delete Movs[c];
      Movs.erase(Movs.begin()+c);
    }
    delete[] use;
  }
  //-Elimina objetos-hijo inutiles
  for(int c=int(Children.size())-1;c>=0;c--){
    if(Children[c]->Optimize()){
      //printf("***Borrando obj(o_%d,r_%d)\n",Children[c]->Id,Children[c]->Ref);
      delete Children[c];
      Children.erase(Children.begin()+c);
    }
  }
  //-Comprueba si los antecesores tienen movimientos
  bool parentmovs=false;
  JMotionObj* obj=(JMotionObj*)Parent; 
  while(obj&&!parentmovs){
    if(obj->Movs.size())parentmovs=true;
    obj=(JMotionObj*)obj->Parent;
  }
  return(Children.size()==0&&(Ref<0||(Movs.size()==0&&!parentmovs)));
}

//==============================================================================
// Guarda configuracion de evento en formato xml
//==============================================================================
void JMotionObj::WriteXml(TiXmlNode* node,const JMotionEvent &evt)const{
  TiXmlElement item("begin");
  JXml::AddAttribute(&item,"mov",int(evt.Mov->Id));
  JXml::AddAttribute(&item,"start",evt.TimeStart);
  if(evt.TimeFinish>=0)JXml::AddAttribute(&item,"finish",evt.TimeFinish);
  node->InsertEndChild(item);
}

//==============================================================================
// Guarda configuracion de objeto en formato xml
//==============================================================================
void JMotionObj::WriteXml(TiXmlNode* node)const{
  TiXmlNode* node2=node->InsertEndChild(Ref>=0? JXml::MakeElementAttrib("objreal","ref",Ref): TiXmlElement("obj"));
  for(unsigned c=0;c<Events.size();c++)WriteXml(node2,*Events[c]);
  for(unsigned c=0;c<Movs.size();c++)Movs[c]->WriteXml(node2);
  for(unsigned c=0;c<Children.size();c++)Children[c]->WriteXml(node2);
}

