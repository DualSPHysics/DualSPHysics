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

/// \file JMotionMov.cpp \brief Implements the classes \ref JMotionMov.

#include "JMotionMov.h"
#include "JXml.h"
#include "JReadDatafile.h"

//##############################################################################
//# JMotionMovXXX
//##############################################################################
//==============================================================================
// Guarda configuracion de movimiento en formato xml
//==============================================================================
void JMotionMovWait::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("wait"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  node->InsertEndChild(item);
}
//==============================================================================
void JMotionMovRect::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrect"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  JXml::AddAttribute(JXml::AddElementDouble3(node2,"vel",Vel),"units_comment","m/s");
}
//==============================================================================
void JMotionMovRectAce::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrectace"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  if(!VelPrev)JXml::AddAttribute(JXml::AddElementDouble3(node2,"velini",Vel),"units_comment","m/s");
  JXml::AddAttribute(JXml::AddElementDouble3(node2,"ace",Ace),"units_comment","m/s^2");
}
//==============================================================================
void JMotionMovRot::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrot"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  const double velang=(AngDegrees? VelAng: VelAng*TORAD);
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"vel","ang",velang),"units_comment",(AngDegrees? "degrees/s": "radians/s"));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}
//==============================================================================
void JMotionMovRotAce::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrotace"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  const double aceang=(AngDegrees? AceAng: AceAng*TORAD);
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"ace","ang",aceang),"units_comment",(AngDegrees? "degrees/s^2": "radians/s^2"));
  const double velang=(AngDegrees? VelAng: VelAng*TORAD);
  if(!VelPrev)JXml::AddAttribute(JXml::AddElementAttrib(node2,"velini","ang",velang),"units_comment",(AngDegrees? "degrees/s": "radians/s"));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}

//==============================================================================
void JMotionMovCir::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvcir"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  const double velang=(AngDegrees? VelAng: VelAng*TORAD);
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"vel","ang",velang),"units_comment",(AngDegrees? "degrees/s": "radians/s"));
  //JXml::AddAttribute(JXml::AddElementDouble3(node2,"ref",Ref->P1),"units_comment","metres (m)");
  node2->InsertEndChild(JXml::MakeElementDouble3("ref",Ref->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}
//==============================================================================
void JMotionMovCirAce::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvcirace"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  const double aceang=(AngDegrees? AceAng: AceAng*TORAD);
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"ace","ang",aceang),"units_comment",(AngDegrees? "degrees/s^2": "radians/s^2"));
  const double velang=(AngDegrees? VelAng: VelAng*TORAD);
  if(!VelPrev)JXml::AddAttribute(JXml::AddElementAttrib(node2,"velini","ang",velang),"units_comment",(AngDegrees? "degrees/s": "radians/s"));
  //JXml::AddAttribute(JXml::AddElementDouble3(node2,"ref",Ref->P1),"units_comment","metres (m)");
  node2->InsertEndChild(JXml::MakeElementDouble3("ref",Ref->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}

//==============================================================================
void JMotionMovRectSinu::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrectsinu"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  JXml::AddAttribute(JXml::AddElementDouble3(node2,"freq",Freq),"units_comment","1/s");
  JXml::AddAttribute(JXml::AddElementDouble3(node2,"ampl",Ampl),"units_comment","metres (m)");
  const tdouble3 phase=(!AngDegrees? Phase: Phase*TODEG);
  if(!PhasePrev)JXml::AddAttribute(JXml::AddElementDouble3(node2,"phase",phase),"units_comment",(AngDegrees? "degrees": "radians"));
}

//==============================================================================
void JMotionMovRotSinu::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrotsinu"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"freq","v",Freq),"units_comment","1/s");
  const double ampl=(AngDegrees? Ampl: Ampl*TORAD);
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"ampl","v",ampl),"units_comment",(AngDegrees? "degrees": "radians"));
  const double phase=(!AngDegrees? Phase: Phase*TODEG);
  if(!PhasePrev)JXml::AddAttribute(JXml::AddElementAttrib(node2,"phase","v",phase),"units_comment",(AngDegrees? "degrees": "radians"));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}

//==============================================================================
void JMotionMovCirSinu::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvcirsinu"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"freq","v",Freq),"units_comment","1/s");
  const double ampl=(AngDegrees? Ampl: Ampl*TORAD);
  JXml::AddAttribute(JXml::AddElementAttrib(node2,"ampl","v",ampl),"units_comment",(AngDegrees? "degrees": "radians"));
  const double phase=(!AngDegrees? Phase: Phase*TODEG);
  if(!PhasePrev)JXml::AddAttribute(JXml::AddElementAttrib(node2,"phase","v",phase),"units_comment",(AngDegrees? "degrees": "radians"));
  node2->InsertEndChild(JXml::MakeElementDouble3("ref",Ref->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}

//==============================================================================
void JMotionMovRectFile::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrectfile"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  //-Define contenido de movimiento.
  TiXmlElement item2("file");
  JXml::AddAttribute(&item2,"name",File);
  JXml::AddAttribute(&item2,"fields",Fields);
  JXml::AddAttribute(&item2,"fieldtime",FieldTime);
  if(FieldX>=0)JXml::AddAttribute(&item2,"fieldx",FieldX);
  if(FieldY>=0)JXml::AddAttribute(&item2,"fieldy",FieldY);
  if(FieldZ>=0)JXml::AddAttribute(&item2,"fieldz",FieldZ);
  TiXmlNode* node2=node->InsertEndChild(item);
  node2->InsertEndChild(item2);
}

//==============================================================================
void JMotionMovRotFile::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrotfile"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Define contenido de movimiento.
  node2->InsertEndChild(JXml::MakeElementAttrib("file","name",File));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}

//==============================================================================
void JMotionMovNull::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvnull"); JXml::AddAttribute(&item,"id",int(Id));
  node->InsertEndChild(item);
}


//##############################################################################
//# JMotionDataFile
//##############################################################################
//==============================================================================
// Constructor for positions.
//==============================================================================
JMotionDataFile::JMotionDataFile(std::string dirdata,std::string file,const int fields,const int fieldtime,const int fieldx,const int fieldy,const int fieldz):PosType(true){
  ClassName="JMotionDataFile";
  //printf("--------> NEW  JMotionDataFile\n");
  Times=NULL; ValuesPos=NULL; ValuesAng=NULL;
  Reset();
  LoadFilePos(dirdata,file,fields,fieldtime,fieldx,fieldy,fieldz);
}

//==============================================================================
// Constructor for angles.
//==============================================================================
JMotionDataFile::JMotionDataFile(std::string dirdata,std::string file,bool angdegrees):PosType(false){
  ClassName="JMotionDataFile";
  //printf("--------> NEW  JMotionDataFile\n");
  Times=NULL; ValuesPos=NULL; ValuesAng=NULL;
  Reset();
  LoadFileAng(dirdata,file,angdegrees);
}

//==============================================================================
// Destructor.
//==============================================================================
JMotionDataFile::~JMotionDataFile(){
  DestructorActive=true;
  //printf("--------> DELETE  JMotionDataFile\n");
  Reset();
}

//==============================================================================
// Initialization of variables.
//==============================================================================
void JMotionDataFile::Reset(){
  Size=Count=0;
  delete[] Times;     Times=NULL;
  delete[] ValuesPos; ValuesPos=NULL;
  delete[] ValuesAng; ValuesAng=NULL;
}

//==============================================================================
// Modifica la asignacion de memoria para los vectores de MovRectFile y MovRotFile
//==============================================================================
void JMotionDataFile::Resize(unsigned size){
  //printf("JMotionMovActive::DfResize> %u --> %u\n",DfSize,size);
  if(size>Size){
    double* times=new double[size];
    tdouble3* vpos=NULL;
    double* vang=NULL;
    if(PosType)vpos=new tdouble3[size];
    else         vang=new double[size];
    if(Size){
      memcpy(times,Times,sizeof(double)*Size);
      if(PosType)memcpy(vpos,ValuesPos,sizeof(tdouble3)*Size);
      else       memcpy(vang,ValuesAng,sizeof(double)  *Size);
      delete[] Times;
      delete[] ValuesPos;
      delete[] ValuesAng;
    }
    Times=times;
    ValuesPos=vpos;
    ValuesAng=vang;
    Size=size;
  }
}

//==============================================================================
// Carga y configura datos del movimiento a partir de fichero de datos.
//==============================================================================
void JMotionDataFile::LoadFilePos(std::string dirdata,std::string file,const int fields,const int fieldtime,const int fieldx,const int fieldy,const int fieldz){
  if(int(file.find("/"))<0 && int(file.find("\\"))<0)file=dirdata+file; //-If only name of file then add Dirin.
  //printf("LoadFilePos>  file:[%s]\n",file.c_str());

  //-Carga datos del fichero.
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  Resize(rows);
  //printf("----->rows:%u\n",rows);
  for(unsigned r=0;r<rows;r++){
    double time,posx=0,posy=0,posz=0;
    for(int f=0;f<fields;f++){
      const double v=rdat.ReadNextDouble(f!=0);
      //printf("[%u]------>  v[%u]:%f\n",r,f,v);
      if(f==fieldtime)time=v;
      if(f==fieldx)posx=v;
      if(f==fieldy)posy=v;
      if(f==fieldz)posz=v;
    }
    Times[r]=time; ValuesPos[r]=TDouble3(posx,posy,posz);
    //printf("[%u]>  t:%f  x:%f\n",r,time,posx);
  }
  Count=rows;
  if(Count<2)Run_ExceptioonFile("Cannot be less than two positions.",file);
}

//==============================================================================
// Carga y configura datos del movimiento a partir de fichero de datos.
//==============================================================================
void JMotionDataFile::LoadFileAng(std::string dirdata,std::string file,bool angdegrees){
  if(int(file.find("/"))<0 && int(file.find("\\"))<0)file=dirdata+file; //-If only name of file then add Dirin.

  //-Carga datos del fichero.
  JReadDatafile rdat;
  rdat.LoadFile(file,SIZEMAX);
  const unsigned rows=rdat.Lines()-rdat.RemLines();
  Resize(rows);
  for(unsigned r=0;r<rows;r++){
    Times[r]=rdat.ReadNextDouble();
    ValuesAng[r]=rdat.ReadNextDouble();
    //printf("FileData[%u]>  t:%f  ang:%f\n",r,time,ang);
  }
  if(!angdegrees)for(unsigned r=0;r<rows;r++)ValuesAng[r]=ValuesAng[r]*TODEG; //-Convierte de radianes a grados.
  Count=rows;
  if(Count<2)Run_ExceptioonFile("Cannot be less than two angles.",file);
}
