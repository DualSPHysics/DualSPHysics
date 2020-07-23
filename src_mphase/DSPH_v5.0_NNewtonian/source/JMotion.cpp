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

/// \file JMotion.cpp \brief Implements the class \ref JMotion.

#include "JMotion.h"
#include "JMotionMov.h"
#include "JMotionObj.h"
#include "JMotionEvent.h"
#include "JMotionList.h"
#include "Functions.h"
#include "JXml.h"
#include <algorithm>

using namespace std;

//##############################################################################
//# JMotion
//##############################################################################
//==============================================================================
// Constructor.
//==============================================================================
JMotion::JMotion(){
  ClassName="JMotion";
  LisMov=NULL;
  LisStop=NULL;
  MotList=NULL;
  Reset();
}

//==============================================================================
// Destructor.
//==============================================================================
JMotion::~JMotion(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JMotion::Reset(){
  //-Libera eventos
  Events.clear();
  //-Libera objetos
  for(unsigned c=0;c<Objs.size();c++)delete Objs[c];
  Objs.clear();
  //-Arrays de objetos en movimiento o parados recientes
  if(LisMov){ delete[] LisMov; LisMov=NULL; }
  LisMovCount=0;
  if(LisStop){ delete[] LisStop; LisStop=NULL; }
  LisStopCount=0;
  ObjCount=0;
  Prepared=false;
  EventNext=0;
  ObjsActive=false;
  DirData="";
  delete MotList; MotList=NULL;
}

//==============================================================================
// Configura el directorio para datos de entrada.
//==============================================================================
void JMotion::SetDirData(const std::string &dirdata){
  DirData=fun::GetDirWithSlash(dirdata);
}

//==============================================================================
// Devuelve puntero al objeto con el id indicado
//==============================================================================
JMotionObj* JMotion::ObjGetPointer(unsigned id)const{
  JMotionObj *obj=NULL;
  for(unsigned c=0;c<Objs.size()&&!obj;c++)obj=Objs[c]->ObjGetPointer(id);
  return(obj);
}

//==============================================================================
// Devuelve puntero al objeto con la referencia indicada
//==============================================================================
JMotionObj* JMotion::ObjGetPointerByRef(int ref)const{
  JMotionObj *obj=NULL;
  for(unsigned c=0;c<Objs.size()&&!obj;c++)obj=Objs[c]->ObjGetPointerByRef(ref);
  return(obj);
}

//==============================================================================
// Devuelve true si existe el objeto solicitado
//==============================================================================
bool JMotion::ExistsObj(JMotionObj* obj)const{
  bool ret=false;
  for(unsigned c=0;c<Objs.size()&&!ret;c++)ret=Objs[c]->ExistsObj(obj);
  return(ret);
}

//==============================================================================
// Anhade un objeto
// - id: Es el identificador del objeto tiene que ser mayor que cero.
// - idparent: Id del objeto padre, es 0 cuando no tiene padre.
// - ref: Referencia (Mk) al objeto real que se mueve. Menor que cero cuando se 
//   trata de objetos virtuales.
//==============================================================================
void JMotion::ObjAdd(unsigned id,unsigned idparent,int ref){
  if(Prepared)Run_Exceptioon("Invalid method in execution mode.");
  if(!id)Run_Exceptioon("Cannot add an object with zero id.");
  if(ObjGetPointer(id))Run_Exceptioon(fun::PrintStr("Cannot add an object with an existing id=%u.",id));
  if(ref>=0&&ObjGetPointerByRef(ref))Run_Exceptioon(fun::PrintStr("Cannot add a new object with an existing real reference (ref=%d) (higher than or equal to zero).",ref));
  JMotionObj* parent=(idparent? ObjGetPointer(idparent): NULL);
  if(idparent&&!parent)Run_Exceptioon(fun::PrintStr("Parent object with id=%u is missing.",idparent));
  JMotionObj* obj=new JMotionObj(id,parent,ref);
  if(!obj)Run_Exceptioon("Could not allocate the requested memory.");
  if(!parent)Objs.push_back(obj);
  else parent->AddChild(obj);
}

//==============================================================================
// Anhade un evento
// - timefinish: Indica el tiempo maximo del movimiento/s iniciados por este 
//   evento, si es menor que cero se ignora.
//==============================================================================
void JMotion::EventAdd(unsigned objid,unsigned movid,double timestart,double timefinish){
  if(Prepared)Run_Exceptioon("Invalid method in execution mode.");
  JMotionObj* obj=ObjGetPointer(objid);
  if(!obj)Run_Exceptioon(fun::PrintStr("Missing object with id=%u.",objid));
  JMotionMov* mov=obj->MovGetPointer(movid);
  if(!mov)Run_Exceptioon(fun::PrintStr("Missing movement (id=%u) inside the object with id=%u.",movid,obj->Id));
  JMotionEvent* evt=new JMotionEvent(obj,mov,timestart,timefinish);
  if(!evt)Run_Exceptioon("Cannot allocate the requested memory.");
  obj->AddEvent(evt);
  Events.push_back(evt);
}

//==============================================================================
// Anhade (si es necesario) un nuevo eje para un objeto y lo devuelve.
//==============================================================================
JMotionAxis* JMotion::AxisAdd(unsigned objid,const tdouble3 &p1,const tdouble3 &p2){
  if(Prepared)Run_Exceptioon("Invalid method in execution mode.");
  JMotionAxis* axis=NULL;
  JMotionObj* obj=ObjGetPointer(objid);
  if(!obj)Run_Exceptioon(fun::PrintStr("Missing object with id=%u.",objid));
  if(p1!=p2)axis=obj->AxisGetPointer(p1,p2);//-No reutiliza ejes usados como referencia (p1==p2) pq esto se modifican.
  if(!axis){
    axis=new JMotionAxis(p1,p2);
    if(!axis)Run_Exceptioon("Cannot allocate the requested memory.");
    obj->AddAxis(axis);
  }   
  return(axis);
}

//==============================================================================
// Anhade un nuevo movimiento
//==============================================================================
void JMotion::MovAdd(unsigned objid,JMotionMov* mov){
  if(Prepared)Run_Exceptioon("Invalid method in execution mode.");
  if(!mov)Run_Exceptioon("Cannot allocated the requested memory.");
  if(!mov->Id)Run_Exceptioon("The movement id must be greater than zero.");
  JMotionObj* obj=ObjGetPointer(objid);
  if(!obj)Run_Exceptioon("Missing object.");
  if(obj->MovGetPointer(mov->Id))Run_Exceptioon("Cannot add a movement with a existing id inside the object.");
  obj->AddMov(mov);
}

//==============================================================================
// Anhade un tiempo de espera
//==============================================================================
void JMotion::MovAddWait(unsigned objid,unsigned id,unsigned nextid,double time){
  if(time<0)Run_Exceptioon("Wating times lenght lower than zero are not allowed.");
  MovAdd(objid,new JMotionMovWait(id,nextid,time));
}
//==============================================================================
// Anhade un desplazamiento instantaneo
//==============================================================================
void JMotion::MovAddTeleport(unsigned objid,unsigned id,unsigned nextid,const tdouble3 &mpos){
  MovAdd(objid,new JMotionMovRect(id,nextid,-1,mpos));
}
//==============================================================================
// Anhade un movimiento rectilineo uniforme
//==============================================================================
void JMotion::MovAddRectilinear(unsigned objid,unsigned id,unsigned nextid,double time,const tdouble3 &vel){
  MovAdd(objid,new JMotionMovRect(id,nextid,time,vel));
}
//==============================================================================
// Anhade un movimiento rectilineo uniformemente acelerado
//==============================================================================
void JMotion::MovAddRectilinearAce(unsigned objid,unsigned id,unsigned nextid,double time,const tdouble3 &ace,const tdouble3 &vel,bool velpre){
  MovAdd(objid,new JMotionMovRectAce(id,nextid,time,ace,vel,velpre));
}
//==============================================================================
// Anhade un movimiento rotacion
//==============================================================================
void JMotion::MovAddRotation(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees,const tdouble3 &axisp1,const tdouble3 &axisp2,double velang,bool useangdegrees){
  if(useangdegrees && !angdegrees)velang=velang*TODEG;
  MovAdd(objid,new JMotionMovRot(id,nextid,time,angdegrees,AxisAdd(objid,axisp1,axisp2),velang));
}
//==============================================================================
// Anhade un movimiento rotacion uniformemente acelerado
//==============================================================================
void JMotion::MovAddRotationAce(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees,const tdouble3 &axisp1,const tdouble3 &axisp2,double aceang,double velang,bool velpre,bool useangdegrees){
  if(useangdegrees && !angdegrees){
    aceang=aceang*TODEG;
    velang=velang*TODEG;
  }
  MovAdd(objid,new JMotionMovRotAce(id,nextid,time,angdegrees,AxisAdd(objid,axisp1,axisp2),aceang,velang,velpre));
}
//==============================================================================
// Anhade un movimiento circular sin rotacion
//==============================================================================
void JMotion::MovAddCircular(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees,const tdouble3 &axisp1,const tdouble3 &axisp2,const tdouble3 &ref,double velang,bool useangdegrees){
  if(useangdegrees && !angdegrees)velang=velang*TODEG;
  MovAdd(objid,new JMotionMovCir(id,nextid,time,angdegrees,AxisAdd(objid,axisp1,axisp2),AxisAdd(objid,ref,ref),velang));
}
//==============================================================================
// Anhade un movimiento circular sin rotacion uniformemente acelerado
//==============================================================================
void JMotion::MovAddCircularAce(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees,const tdouble3 &axisp1,const tdouble3 &axisp2,const tdouble3 &ref,double aceang,double velang,bool velpre,bool useangdegrees){
  if(useangdegrees && !angdegrees){
    aceang=aceang*TODEG;
    velang=velang*TODEG;
  }
  MovAdd(objid,new JMotionMovCirAce(id,nextid,time,angdegrees,AxisAdd(objid,axisp1,axisp2),AxisAdd(objid,ref,ref),aceang,velang,velpre));
}
//==============================================================================
// Anhade un movimiento rectilineo sinusoidal
//==============================================================================
void JMotion::MovAddRecSinu(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees,const tdouble3 &freq,const tdouble3 &ampl,tdouble3 phase,bool phaseprev,bool useangdegrees){
  if(useangdegrees && angdegrees)phase=phase*TDouble3(TORAD); //-Convierte a radianes.
  MovAdd(objid,new JMotionMovRectSinu(id,nextid,time,angdegrees,freq,ampl,phase,phaseprev));
}
//==============================================================================
// Anhade un movimiento de rotacion sinusoidal
//==============================================================================
void JMotion::MovAddRotSinu(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees,const tdouble3 &axisp1,const tdouble3 &axisp2,double freq,double ampl,double phase,bool phaseprev,bool useangdegrees){
  if(useangdegrees && !angdegrees)ampl=ampl*TODEG;
  if(useangdegrees && angdegrees)phase=phase*TORAD; //-Convierte a radianes.
  MovAdd(objid,new JMotionMovRotSinu(id,nextid,time,angdegrees,AxisAdd(objid,axisp1,axisp2),freq,ampl,phase,phaseprev));
}
//==============================================================================
// Anhade un movimiento circular sinusoidal
//==============================================================================
void JMotion::MovAddCirSinu(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees,const tdouble3 &axisp1,const tdouble3 &axisp2,const tdouble3 &ref,double freq,double ampl,double phase,bool phaseprev,bool useangdegrees){
  if(useangdegrees && !angdegrees)ampl=ampl*TODEG;
  if(useangdegrees && angdegrees)phase=phase*TORAD; //-Convierte a radianes.
  MovAdd(objid,new JMotionMovCirSinu(id,nextid,time,angdegrees,AxisAdd(objid,axisp1,axisp2),AxisAdd(objid,ref,ref),freq,ampl,phase,phaseprev));
}
//==============================================================================
// Anhade un movimiento rectilineo a partir de datos de un fichero.
// - fields: Numero total de campos en el fichero.
// - fieldtime: Posicion del campo time dentro de fields.
// - fieldx: Posicion del campo x dentro de fields (menor que 0 se ignora).
// - fieldy: Posicion del campo y dentro de fields (menor que 0 se ignora).
// - fieldz: Posicion del campo z dentro de fields (menor que 0 se ignora).
//==============================================================================
void JMotion::MovAddRectilinearFile(unsigned objid,unsigned id,unsigned nextid,double time,const std::string &file,int fields,int fieldtime,int fieldx,int fieldy,int fieldz){
  if(fieldtime<0)Run_Exceptioon("The \'time\' is not defined.");
  if(fieldtime>=0 && fieldtime>=fields)Run_Exceptioon("the position of field \'time\' is invalid.");
  if(fieldx<0 && fieldy<0 && fieldz<0)Run_Exceptioon("You need at least one position field.");
  if(fieldx>=0 && fieldx>=fields)Run_Exceptioon("the position of field \'x\' is invalid.");
  if(fieldy>=0 && fieldy>=fields)Run_Exceptioon("the position of field \'y\' is invalid.");
  if(fieldz>=0 && fieldz>=fields)Run_Exceptioon("the position of field \'z\' is invalid.");
  MovAdd(objid,new JMotionMovRectFile(id,nextid,time,&DirData,file,fields,fieldtime,fieldx,fieldy,fieldz));
}
//==============================================================================
// Anhade un movimiento de rotacion a partir de datos de un fichero.
//==============================================================================
void JMotion::MovAddRotationFile(unsigned objid,unsigned id,unsigned nextid,double time,bool angdegrees
  ,const tdouble3 &axisp1,const tdouble3 &axisp2,const std::string &file)
{
  MovAdd(objid,new JMotionMovRotFile(id,nextid,time,angdegrees,AxisAdd(objid,axisp1,axisp2),&DirData,file));
}

//==============================================================================
// Anhade un movimiento nulo
//==============================================================================
void JMotion::MovAddNull(unsigned objid,unsigned id){
  MovAdd(objid,new JMotionMovNull(id));
}


//==============================================================================
// Revisa y actualiza vinculos entre movimientos de los objetos.
//==============================================================================
void JMotion::CheckLinkMovs()const{
  for(unsigned c=0;c<Objs.size();c++)Objs[c]->LinkMovs(); 
}

//==============================================================================
// Prepara objeto revisando y actualizando vinculos.
//==============================================================================
void JMotion::Prepare(){
  if(Prepared)Run_Exceptioon("Invalid method in execution mode.");
  CheckLinkMovs();
  //-Ordena eventos de ultimo a primero.
  if(Events.size())for(unsigned c=0;c<Events.size()-1;c++)for(unsigned c2=c+1;c2<Events.size();c2++)if(Events[c]->TimeStart<Events[c2]->TimeStart)swap(Events[c],Events[c2]);
//for(unsigned c=0;c<Events.size();c++)printf("Evt[%d].start: %G\n",c,Events[c]->TimeStart);
  EventNext=int(Events.size())-1;
  //-Contabiliza el numero de objetos totales
  ObjCount=0;
  for(unsigned c=0;c<Objs.size();c++)ObjCount+=1+Objs[c]->ChildrenCount(); 
  LisMov=new JMotionObj*[ObjCount];
  LisStop=new JMotionObj*[ObjCount];
  Prepared=true;
  ObjsActive=false;
}

//==============================================================================
// Crea y prepara objeto MotList.
//==============================================================================
void JMotion::CreateMotList(){
  //-Carga lista de referencias.
  std::vector<int> refs;
  for(unsigned c=0;c<Objs.size();c++)Objs[c]->GetRefs(refs);
  //-Ordena lista de referencias.
  const int nref=int(refs.size());
  for(int c=0;c<nref-1;c++)for(int c2=c+1;c2<nref;c2++)if(refs[c]>refs[c2]){ 
    int aux=refs[c]; refs[c]=refs[c2]; refs[c2]=aux; 
  }
  //-Calcula y comprueba referencia maxima.
  int maxref=0;
  for(int c=0;c<nref;c++){
    maxref=max(maxref,refs[c]);
    //printf("---> %d \n",refs[c]);
  }
  if(maxref+1!=nref)Run_Exceptioon("Motion references are no consecutives.");
  //-Creates object MotList.
  MotList=new JMotionList(nref);
}

//==============================================================================
// Sistema original de calculo pero guardando resultados en MotList igual que ProcesTimeAce().
// Revisa lista de eventos para crear nuevos movimientos activos.
// Devuelve true si hay movimientos activos.
//==============================================================================
bool JMotion::ProcesTimeSimple(double timestep,double dt){
  if(MotList==NULL)CreateMotList();
  MotList->PreMotion();
  bool movs=false;
  if(MotList->TimeStep>timestep)Run_Exceptioon("The previous timestep+dt is higher than the requested timestep. It is invalid for simple mode.");
  //-Procesa primer movimiento de dt.
  if(ProcesTime(timestep,dt)){
    movs=true;
    const unsigned nmove=GetMovCount();
    for(unsigned c=0;c<nmove;c++){
      unsigned ref;
      tdouble3 mvsimple;
      JMatrix4d mvmatrix;
      if(GetMov(c,ref,mvsimple,mvmatrix))MotList->Sp_Movedt(ref,mvsimple,dt);//-Simple movement. | Movimiento simple.
      else MotList->Sp_Movedt(ref,mvmatrix.GetMatrix4d(),dt); //-Movement using a matrix. | Movimiento con matriz.
    }
  }
  MotList->TimeStep=timestep+dt;
  return(movs);
}

//==============================================================================
// Sistema nuevo para calcular velocidad y aceleracion, reiniciando calculo y 
// usando posiciones anteriores y futuras.
// Revisa lista de eventos para crear nuevos movimientos activos.
// Devuelve true si hay movimientos activos.
//==============================================================================
bool JMotion::ProcesTimeAce(double timestep,double dt){
  if(MotList==NULL)CreateMotList();
  MotList->PreMotion();
  bool movs=false;
  //-Procesa primer movimiento de dt.
  ResetTime(timestep);
  if(ProcesTime(timestep,dt)){
    movs=true;
    const unsigned nmove=GetMovCount();
    for(unsigned c=0;c<nmove;c++){
      unsigned ref;
      tdouble3 mvsimple;
      JMatrix4d mvmatrix;
      if(GetMov(c,ref,mvsimple,mvmatrix))MotList->Ace2_Move1dt(ref,mvsimple);//-Simple movement. | Movimiento simple.
      else MotList->Ace2_Move1dt(ref,mvmatrix.GetMatrix4d()); //-Movement using a matrix. | Movimiento con matriz.
    }
  }
  //-Procesa segundo movimiento de 2*dt.
  ResetTime(timestep);
  if(ProcesTime(timestep,dt+dt)){
    movs=true;
    const unsigned nmove=GetMovCount();
    for(unsigned c=0;c<nmove;c++){
      unsigned ref;
      tdouble3 mvsimple;
      JMatrix4d mvmatrix;
      if(GetMov(c,ref,mvsimple,mvmatrix))MotList->Ace2_Move2dt(ref,mvsimple);//-Simple movement. | Movimiento simple.
      else MotList->Ace2_Move2dt(ref,mvmatrix.GetMatrix4d()); //-Movement using a matrix. | Movimiento con matriz.
    }
  }
  MotList->TimeStep=timestep+dt+dt;
  //-Procesa movimientos.
  MotList->Ace2_PosMotion(dt);
  return(movs);
}

//==============================================================================
// Nuevo metodo para devolver resultados de ProcesTimeSimple() o ProcesTimeAce().
// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JMotion::ProcesTimeGetData(unsigned ref,bool &typesimple,tdouble3 &simplemov
  ,tdouble3 &simplevel,tdouble3 &simpleace,tmatrix4d &matmov,tmatrix4d &matmov2)const
{
  return(MotList->GetData(ref,typesimple,simplemov,simplevel,simpleace,matmov,matmov2));
}

//==============================================================================
// Nuevo metodo para devolver resultados de ProcesTimeSimple() o ProcesTimeAce().
// Returns data of one moving object. Returns true when the motion is active.
//==============================================================================
bool JMotion::ProcesTimeGetData(unsigned ref,bool &typesimple,tdouble3 &simplemov,tmatrix4d &matmov)const
{
  return(MotList->GetData(ref,typesimple,simplemov,matmov));
}

//==============================================================================
// Reinicia ejecucion al timestep indicado.
//==============================================================================
void JMotion::ResetTime(double timestep){
  EventNext=int(Events.size())-1;
  LisMovCount=0;
  LisStopCount=0;
  for(unsigned c=0;c<Objs.size();c++)Objs[c]->ResetTime(); 
  ObjsActive=false;
  //-Procesa movimiento hasta timestep.
  if(timestep>0)ProcesTime(0,timestep);
}

//==============================================================================
// Revisa lista de eventos para crear nuevos movimientos activos.
// Devuelve true si hay movimientos activos.
//==============================================================================
bool JMotion::ProcesTime(double timestep,double dt){
  //printf("ProcesTime> timestep:%f dt:%f  \n",timestep,dt);
  if(!Prepared)Run_Exceptioon("Invalid method in initialization mode.");
  //-Comprueba eventos para activar nuevos movimientos.
  bool looking=true;
  if(EventNext>=0)for(int c=EventNext;c>=0&&looking;c--){
    if(Events[c]->TimeStart<(timestep+dt)){
      JMotionEvent* evt=Events[c];
      //printf("Motion> Time %g  Empieza mov -> start:%g  obj:%d  mov:%d\n",timestep,evt->TimeStart,evt->Obj->Id,evt->Mov->Id);
      evt->Obj->BeginEvent(evt->TimeStart,evt->TimeFinish,evt->Mov);
      EventNext--;
      ObjsActive=true;
      //printf("EventNext-----> %d\n",EventNext);
    }
    else looking=false;
  }
  //-Calcula movimiento de cada objeto y registra objetos en movimiento o recien parados.
  if(ObjsActive){
    ObjsActive=false;
    LisMovCount=0; LisStopCount=0;
    for(unsigned c=0;c<Objs.size();c++)if(Objs[c]->Active)ObjsActive|=Objs[c]->ProcesTime(timestep,dt,LisMov,LisMovCount,LisStop,LisStopCount);
  }
  return(ObjsActive);//return((LisMovCount+LisStopCount)!=0);
}

//==============================================================================
// Devuelve datos de movimiento de un objeto.
//==============================================================================
bool JMotion::GetMov(unsigned pos,unsigned &ref,tdouble3 &mvsimple,JMatrix4d &mvmatrix)const{
  if(pos>=LisMovCount)Run_Exceptioon("The requested movement is invalid.");
  return(LisMov[pos]->GetMov(ref,mvsimple,mvmatrix));
}

//==============================================================================
// Devuelve la referencia del objeto parado.
//==============================================================================
unsigned JMotion::GetStopRef(unsigned pos)const{
  if(pos>=LisStopCount)Run_Exceptioon("End of invalid requested movement.");
  return(LisStop[pos]->Ref);
}

//==============================================================================
// Devuelve la mayor referencia de objeto definida.
//==============================================================================
int JMotion::GetMaxRef()const{
  int ref=-1;
  for(unsigned c=0;c<Objs.size();c++){
    int rf=Objs[c]->GetMaxRef();
    if(rf>ref)ref=rf;
  }
  return(ref);
}

//==============================================================================
// Copia los datos de configuracion (Objs,Movs,Events) a mot.
//==============================================================================
void JMotion::CopyConfig(JMotion &mot)const{
  mot.Reset();
  mot.DirData=DirData;
  for(int c=0;c<int(Objs.size());c++)Objs[c]->CopyConfig(mot);
  for(int c=0;c<int(Events.size());c++)mot.EventAdd(Events[c]->Obj->Id,Events[c]->Mov->Id,Events[c]->TimeStart,Events[c]->TimeFinish);
  mot.CheckLinkMovs();
}

//==============================================================================
// Copia configuracion a mot intercambiando las referencias por su nuevo valor.
// Solo copia los datos de configuracion (Objs,Movs,Events).
// Las referencias que no aparezcan pasan a ser -1.
//==============================================================================
void JMotion::CopyChangeRef(JMotion &mot,const int* ref,const int* refnew,unsigned refcount)const{
  mot.Reset();
  for(int c=0;c<int(Objs.size());c++)Objs[c]->CopyChangeRef(mot,ref,refnew,refcount);
  for(int c=0;c<int(Events.size());c++)mot.EventAdd(Events[c]->Obj->Id,Events[c]->Mov->Id,Events[c]->TimeStart,Events[c]->TimeFinish);
  mot.CheckLinkMovs();
}

//==============================================================================
// Optimiza configuracion eliminando objetos, movimientos y eventos inutiles.
//==============================================================================
void JMotion::Optimize(){
  //-Elimina objetos inutiles
  for(int c=int(Objs.size())-1;c>=0;c--){
    if(Objs[c]->Optimize()){
      //printf("+++Borrando obj(o_%d,r_%d)\n",Objs[c]->Id,Objs[c]->Ref);
      delete Objs[c];
      Objs.erase(Objs.begin()+c);
    }
  }
  //-Elimina eventos inutiles
  for(int c=int(Events.size())-1;c>=0;c--){
    if(!ExistsObj(Events[c]->Obj)){
      //printf("Borrando evt[%d]\n",c);
      //delete Events[c];  //-Los eventos se liberan dentro de cada JMotionObj propietario.
      Events.erase(Events.begin()+c);
    }
  } 
}

//==============================================================================
// Guarda configuracion de motion en formato xml
//==============================================================================
void JMotion::WriteXml(JXml *jxml,const std::string &path)const{
  jxml->RemoveNode(path);
  CheckLinkMovs();
  TiXmlNode* node=jxml->GetNode(path,true);
  for(unsigned c=0;c<Objs.size();c++)Objs[c]->WriteXml(node);
}

//==============================================================================
// Carga configuracion de motion en formato xml
//==============================================================================
void JMotion::ReadXml(const std::string &dirdata,JXml *jxml,TiXmlNode* node,unsigned &id,unsigned idp){
  TiXmlElement* ele=node->FirstChildElement(); 
  while(ele){
    string name=ele->Value();
    if(name.length()&&name[0]!='_'){
      if(name=="obj"||name=="objreal"){
        id++; ObjAdd(id,idp,(name=="objreal"? jxml->GetAttributeInt(ele,"ref"): -1));
        //printf("ObjAdd(%d,%d,%d)\n",id,idp,(name=="objreal"? jxml->GetAttributeInt(ele,"ref"): -1));
        ReadXml(dirdata,jxml,ele,id,id);
      }
      else if(name=="wait"||name=="mvrect"||name=="mvrectace"||name=="mvrot"||name=="mvrotace"||name=="mvcir"||name=="mvcirace"||name=="mvrectsinu"||name=="mvrotsinu"||name=="mvcirsinu"||name=="mvpredef"||name=="mvfile"||name=="mvrectfile"||name=="mvrotfile"||name=="mvnull"){
        int mvid=jxml->GetAttributeInt(ele,"id");
        double time=0;
        int nextid=0;
        bool angdegrees=true;
        if(name!="mvnull"){
          time=jxml->GetAttributeFloat(ele,"duration");
          nextid=jxml->GetAttributeInt(ele,"next",true,0);
          string txunits=jxml->GetAttributeStr(ele,"anglesunits",true,"degrees");
          if(txunits!="degrees" && txunits!="radians")jxml->ErrReadAtrib(ele,"anglesunits",false);
          angdegrees=(txunits=="degrees");
        }
        if(name=="wait"){
          MovAddWait(idp,mvid,nextid,time);
          //printf("MovAddWait(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvrect"){
          tdouble3 vel=jxml->ReadElementDouble3(ele,"vel");
          MovAddRectilinear(idp,mvid,nextid,time,vel);
          //printf("MovAddRectilinear(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvrectace"){
          tdouble3 ace=jxml->ReadElementDouble3(ele,"ace");
          bool velpre=(ele->FirstChildElement("velini")==NULL);
          tdouble3 velini=(!velpre? jxml->ReadElementDouble3(ele,"velini"): TDouble3(0));
          MovAddRectilinearAce(idp,mvid,nextid,time,ace,velini,velpre);
          //printf("MovAddRectilinearAce(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvrot"){
          tdouble3 axisp1=jxml->ReadElementDouble3(ele,"axisp1");
          tdouble3 axisp2=jxml->ReadElementDouble3(ele,"axisp2");
          double vel=jxml->ReadElementDouble(ele,"vel","ang");
          MovAddRotation(idp,mvid,nextid,time,angdegrees,axisp1,axisp2,vel,true);
          //printf("MovAddRotation(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvrotace"){
          tdouble3 axisp1=jxml->ReadElementDouble3(ele,"axisp1");
          tdouble3 axisp2=jxml->ReadElementDouble3(ele,"axisp2");
          double ace=jxml->ReadElementDouble(ele,"ace","ang");
          bool velpre=(ele->FirstChildElement("velini")==NULL);
          double velini=(!velpre? jxml->ReadElementDouble(ele,"velini","ang"): 0);
          MovAddRotationAce(idp,mvid,nextid,time,angdegrees,axisp1,axisp2,ace,velini,velpre,true);      
          //printf("MovAddRotationAce(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvcir"){
          tdouble3 axisp1=jxml->ReadElementDouble3(ele,"axisp1");
          tdouble3 axisp2=jxml->ReadElementDouble3(ele,"axisp2");
          tdouble3 ref=jxml->ReadElementDouble3(ele,"ref");
          double vel=jxml->ReadElementDouble(ele,"vel","ang");
          MovAddCircular(idp,mvid,nextid,time,angdegrees,axisp1,axisp2,ref,vel,true);       
          //printf("MovAddCircular(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvcirace"){
          tdouble3 axisp1=jxml->ReadElementDouble3(ele,"axisp1");
          tdouble3 axisp2=jxml->ReadElementDouble3(ele,"axisp2");
          tdouble3 ref=jxml->ReadElementDouble3(ele,"ref");
          double ace=jxml->ReadElementDouble(ele,"ace","ang");
          bool velpre=(ele->FirstChildElement("velini")==NULL);
          double velini=(!velpre? jxml->ReadElementDouble(ele,"velini","ang"): 0);
          MovAddCircularAce(idp,mvid,nextid,time,angdegrees,axisp1,axisp2,ref,ace,velini,velpre,true);      
          //printf("MovAddCircularAce(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvrectsinu"){
          tdouble3 freq=jxml->ReadElementDouble3(ele,"freq");
          tdouble3 ampl=jxml->ReadElementDouble3(ele,"ampl");
          bool phaseprev=(ele->FirstChildElement("phase")==NULL);
          tdouble3 phase=(!phaseprev? jxml->ReadElementDouble3(ele,"phase"): TDouble3(0));
          MovAddRecSinu(idp,mvid,nextid,time,angdegrees,freq,ampl,phase,phaseprev,true);        
          //printf("MovAddRecSinu(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvrotsinu"){
          tdouble3 axisp1=jxml->ReadElementDouble3(ele,"axisp1");
          tdouble3 axisp2=jxml->ReadElementDouble3(ele,"axisp2");
          double freq=jxml->ReadElementDouble(ele,"freq","v");
          double ampl=jxml->ReadElementDouble(ele,"ampl","v");
          bool phaseprev=(ele->FirstChildElement("phase")==NULL);
          double phase=(!phaseprev? jxml->ReadElementDouble(ele,"phase","v"): 0);
          MovAddRotSinu(idp,mvid,nextid,time,angdegrees,axisp1,axisp2,freq,ampl,phase,phaseprev,true);      
          //printf("MovAddRotSinu(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvcirsinu"){
          tdouble3 axisp1=jxml->ReadElementDouble3(ele,"axisp1");
          tdouble3 axisp2=jxml->ReadElementDouble3(ele,"axisp2");
          tdouble3 ref=jxml->ReadElementDouble3(ele,"ref");
          double freq=jxml->ReadElementDouble(ele,"freq","v");
          double ampl=jxml->ReadElementDouble(ele,"ampl","v");
          bool phaseprev=(ele->FirstChildElement("phase")==NULL);
          double phase=(!phaseprev? jxml->ReadElementDouble(ele,"phase","v"): 0);
          MovAddCirSinu(idp,mvid,nextid,time,angdegrees,axisp1,axisp2,ref,freq,ampl,phase,phaseprev,true);      
          //printf("MovAddCirSinu(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvpredef" || name=="mvfile" || name=="mvrectfile"){
          TiXmlElement* efile=jxml->GetFirstElement(ele,"file");
          string file=jxml->GetAttributeStr(efile,"name");
          int fields=jxml->GetAttributeInt(efile,"fields");
          int fieldtime=jxml->GetAttributeInt(efile,"fieldtime");
          int fieldx=jxml->GetAttributeInt(efile,"fieldx",true,-1);
          int fieldy=jxml->GetAttributeInt(efile,"fieldy",true,-1);
          int fieldz=jxml->GetAttributeInt(efile,"fieldz",true,-1);
          MovAddRectilinearFile(idp,mvid,nextid,time,file,fields,fieldtime,fieldx,fieldy,fieldz);
          //printf("MovAddRectilinearFile(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvrotfile"){
          tdouble3 axisp1=jxml->ReadElementDouble3(ele,"axisp1");
          tdouble3 axisp2=jxml->ReadElementDouble3(ele,"axisp2");
          TiXmlElement* efile=jxml->GetFirstElement(ele,"file");
          string file=jxml->GetAttributeStr(efile,"name");
          MovAddRotationFile(idp,mvid,nextid,time,angdegrees,axisp1,axisp2,file);      
          //printf("MovAddRotationFile(%d,%d,%d,%g)\n",id,mvid,nextid,time);
        }
        else if(name=="mvnull"){
          MovAddNull(idp,mvid);
          //printf("MovAddNull(%d,%d)\n",id,mvid);
        }
      }
      else if(!idp||name!="begin")jxml->ErrReadElement(ele,name,false);
    }
    ele=ele->NextSiblingElement();
  }
  //-Carga de eventos
  if(idp){
    TiXmlElement* ele=node->FirstChildElement("begin"); 
    while(ele){
      int mvid=jxml->GetAttributeInt(ele,"mov");
      double start=jxml->GetAttributeFloat(ele,"start");
      double finish=jxml->GetAttributeFloat(ele,"finish",true,-1);
      EventAdd(idp,mvid,start,finish);
      //printf("EventAdd(%d,%d,%g,%g)\n",idp,mvid,start,finish);
      ele=ele->NextSiblingElement("begin");
    }
  }
}

//==============================================================================
// Carga configuracion de motion en formato xml
//==============================================================================
void JMotion::ReadXml(const std::string &dirdata,JXml *jxml,const std::string &path,bool checkexists){
  Reset();
  DirData=fun::GetDirWithSlash(dirdata);
  unsigned id=0;
  TiXmlNode* node=jxml->GetNode(path,false);
  if(node){
    ReadXml(DirData,jxml,node,id,0);
    CheckLinkMovs();
  }
  else if(checkexists)Run_Exceptioon(string("Cannot fin the element \'")+path+"\'.");
}

//==============================================================================
// Carga configuracion de motion en formato xml de un fichero
//==============================================================================
void JMotion::LoadFileXml(const std::string &dirdata,const std::string &file,const string &path){
  JXml jxml;
  jxml.LoadFile(file);
  ReadXml(dirdata,&jxml,path);
}

//==============================================================================
// Graba configuracion de motion en formato xml en un fichero
//==============================================================================
void JMotion::SaveFileXml(const std::string &file,const std::string &path,bool newfile)const{
  JXml jxml;
  if(!newfile)jxml.LoadFile(file);
  WriteXml(&jxml,path);
  jxml.SaveFile(file);
}

