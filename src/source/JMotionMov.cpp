//HEAD_DSCODES
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

/// \file JMotionMov.cpp \brief Implements the classes \ref JMotionMov.

#include "JMotionMov.h"
#include "JXml.h"

//##############################################################################
//# JMotionMovWait
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
//==============================================================================
void JMotionMovWait::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("wait"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  node->InsertEndChild(item);
}

//##############################################################################
//# JMotionMovRect
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovRectAce
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovRot
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovRotAce
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovCir
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovCirAce
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovRectSinu
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovRotSinu
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovCirSinu
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
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

//##############################################################################
//# JMotionMovRectFile
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
//==============================================================================
void JMotionMovRectFile::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrectfile"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Saves file configuration.
  TiXmlElement item2("file");
  JXml::AddAttribute(&item2,"name",File);
  JXml::AddAttribute(&item2,"fields",Fields);
  JXml::AddAttribute(&item2,"fieldtime",FieldTime);
  if(FieldX>=0)JXml::AddAttribute(&item2,"fieldx",FieldX);
  if(FieldY>=0)JXml::AddAttribute(&item2,"fieldy",FieldY);
  if(FieldZ>=0)JXml::AddAttribute(&item2,"fieldz",FieldZ);
  node2->InsertEndChild(item2);
}

//##############################################################################
//# JMotionMovRotFile
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
//==============================================================================
void JMotionMovRotFile::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrotfile"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Saves file configuration.
  node2->InsertEndChild(JXml::MakeElementAttrib("file","name",File));
  //-Saves other configuration.
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp1",Axis->P1));
  node2->InsertEndChild(JXml::MakeElementDouble3("axisp2",Axis->P2));
}

//##############################################################################
//# JMotionMovRotAdvFile
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
//==============================================================================
void JMotionMovRotAdvFile::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvrotadvfile"); 
  JXml::AddAttribute(&item,"id",int(Id)); 
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Saves file configuration.
  TiXmlElement item2("file");
  JXml::AddAttribute(&item2,"name",File);
  JXml::AddAttribute(&item2,"fields",Fields);
  JXml::AddAttribute(&item2,"fieldtime",FieldTime);
  if(FieldAng1>=0)JXml::AddAttribute(&item2,"fieldang1",FieldAng1);
  if(FieldAng2>=0)JXml::AddAttribute(&item2,"fieldang2",FieldAng2);
  if(FieldAng3>=0)JXml::AddAttribute(&item2,"fieldang3",FieldAng3);
  node2->InsertEndChild(item2);
  //-Saves other configuration.
  node2->InsertEndChild(JXml::MakeElementDouble3("center",Center));
  node2->InsertEndChild(JXml::MakeElementAttrib("intrinsic","value",Intrinsic));
  node2->InsertEndChild(JXml::MakeElementAttrib("axes","value",Axes));
}

//##############################################################################
//# JMotionMovPathFile
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
//==============================================================================
void JMotionMovPathFile::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvpathfile");
  JXml::AddAttribute(&item,"id",int(Id));
  JXml::AddAttribute(&item,"duration",Time);
  if(NextId)JXml::AddAttribute(&item,"next",int(NextId));
  JXml::AddAttribute(&item,"anglesunits",(AngDegrees? "degrees": "radians"));
  TiXmlNode* node2=node->InsertEndChild(item);
  //-Saves file configuration.
  TiXmlElement item2("file");
  JXml::AddAttribute(&item2,"name",File);
  JXml::AddAttribute(&item2,"fields",Fields);
  JXml::AddAttribute(&item2,"fieldtime",FieldTime);
  if(FieldX>=0)JXml::AddAttribute(&item2,"fieldx",FieldX);
  if(FieldY>=0)JXml::AddAttribute(&item2,"fieldy",FieldY);
  if(FieldZ>=0)JXml::AddAttribute(&item2,"fieldz",FieldZ);
  if(FieldAng1>=0)JXml::AddAttribute(&item2,"fieldang1",FieldAng1);
  if(FieldAng2>=0)JXml::AddAttribute(&item2,"fieldang2",FieldAng2);
  if(FieldAng3>=0)JXml::AddAttribute(&item2,"fieldang3",FieldAng3);
  node2->InsertEndChild(item2);
  //-Saves other configuration.
  node2->InsertEndChild(JXml::MakeElementDouble3("center",Center));
  node2->InsertEndChild(JXml::MakeElementAttrib("movecenter","value",MoveCenter));
  node2->InsertEndChild(JXml::MakeElementAttrib("intrinsic","value",Intrinsic));
  node2->InsertEndChild(JXml::MakeElementAttrib("axes","value",Axes));
}

//##############################################################################
//# JMotionMovNull
//##############################################################################
//==============================================================================
/// Saves configuration in XML.
//==============================================================================
void JMotionMovNull::WriteXml(TiXmlNode* node)const{
  TiXmlElement item("mvnull"); JXml::AddAttribute(&item,"id",int(Id));
  node->InsertEndChild(item);
}

