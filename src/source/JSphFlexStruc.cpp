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

/// \file JSphFlexStruc.cpp \brief Implements the class \ref JSphFlexStruc.

#include "JSphFlexStruc.h"
#include "JAppInfo.h"
#include "JXml.h"
#include "Functions.h"
#include "JSphMk.h"

using namespace std;

//##############################################################################
//# JSphFlexStrucBody
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphFlexStrucBody::JSphFlexStrucBody(
    unsigned idbody,word mkbound,word mkclamp,float particlevolume,float density,double youngmod,double poissratio,TpConstitModel constitmodel,float hgfactor)
    :Log(AppInfo.LogPtr()),IdBody(idbody),MkBound(mkbound),MkClamp(mkclamp)
{
  ClassName="JSphFlexStrucBody";
  Reset();
  ParticleVolume=particlevolume;
  Density=density;
  YoungMod=static_cast<float>(youngmod);
  PoissRatio=static_cast<float>(poissratio);
  ConstitModel=constitmodel;
  HgFactor=hgfactor;
  //-Compute the constitutive model matrix.
  if(ConstitModel==CONSTITMODEL_PlaneStrain || ConstitModel==CONSTITMODEL_SVK){
    //-Calculate coefficients for plane strain and St Venant Kirchoff models
    double mulFactor=youngmod/((1.0+poissratio)*(1.0-2.0*poissratio));
    double c1=mulFactor*(1.0-poissratio);
    double c2=mulFactor*poissratio;
    double c3=mulFactor*(1.0-2.0*poissratio);
    //-Fill rows and columns for 2D.
    ConstitMatrix.a11=ConstitMatrix.a33=static_cast<float>(c1);
    ConstitMatrix.a13=ConstitMatrix.a31=static_cast<float>(c2);
    ConstitMatrix.a55=static_cast<float>(c3);
    //-If St Venant-Kirchoff then fill the 3D rows too.
    if(ConstitModel==CONSTITMODEL_SVK){
      ConstitMatrix.a22=static_cast<float>(c1);
      ConstitMatrix.a12=ConstitMatrix.a21=ConstitMatrix.a23=ConstitMatrix.a32=static_cast<float>(c2);
      ConstitMatrix.a44=ConstitMatrix.a66=static_cast<float>(c3);
    }
  }
  else if(ConstitModel==CONSTITMODEL_PlaneStress){
    //-Calculate coefficients for plane stress.
    double mulFactor=youngmod/(1.0-poissratio*poissratio);
    double c1=mulFactor*1.0;
    double c2=mulFactor*poissratio;
    double c3=mulFactor*(1.0-poissratio);
    ConstitMatrix.a11=ConstitMatrix.a33=static_cast<float>(c1);
    ConstitMatrix.a13=ConstitMatrix.a31=static_cast<float>(c2);
    ConstitMatrix.a55=static_cast<float>(c3);
  }
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphFlexStrucBody::~JSphFlexStrucBody(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphFlexStrucBody::Reset(){
  BoundCode=0;
  ClampCode=0;
  ParticleVolume=0;
  Density=0;
  YoungMod=0;
  PoissRatio=0;
  ConstitModel=CONSTITMODEL_None;
  HgFactor=0;
  ConstitMatrix=TMatrix6f(0);
}

//==============================================================================
/// Configures BoundCode.
//==============================================================================
void JSphFlexStrucBody::ConfigBoundCode(typecode boundcode){
  if(BoundCode)Run_Exceptioon(fun::PrintStr("BoundCode was already configured for mkbound=%u.",MkBound));
  BoundCode=boundcode;
}

//==============================================================================
/// Configures ClampCode.
//==============================================================================
void JSphFlexStrucBody::ConfigClampCode(typecode clampcode){
  if(ClampCode)Run_Exceptioon(fun::PrintStr("ClampCode was already configured for mkbound=%u.",MkBound));
  ClampCode=clampcode;
}


//##############################################################################
//# JSphFlexStruc
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JSphFlexStruc::JSphFlexStruc(bool simulate2d,double dp,JXml *sxml,const std::string &place,const JSphMk *mkinfo)
    :Log(AppInfo.LogPtr()),Simulate2D(simulate2d),Dp(dp)
{
  ClassName="JSphFlexStruc";
  Reset();
  LoadXml(sxml,place);
  ConfigBoundCode(mkinfo);
  ConfigClampCode(mkinfo);
}

//==============================================================================
/// Destructor.
//==============================================================================
JSphFlexStruc::~JSphFlexStruc(){
  DestructorActive=true;
  Reset();
}

//==============================================================================
/// Initialisation of variables.
//==============================================================================
void JSphFlexStruc::Reset(){
  for(unsigned c=0;c<GetCount();c++)delete List[c];
  List.clear();
}

//==============================================================================
/// Returns true if mkbound value is already configured.
//==============================================================================
bool JSphFlexStruc::ExistMk(word mkbound)const{
  bool ret=false;
  for(unsigned c=0;c<List.size() && !ret;c++)ret=(List[c]->MkBound==mkbound);
  return(ret);
}

//==============================================================================
/// Loads initial conditions of XML object.
//==============================================================================
void JSphFlexStruc::LoadXml(const JXml *sxml,const std::string &place){
  TiXmlNode* node=sxml->GetNodeSimple(place);
  if(!node)Run_Exceptioon(std::string("Cannot find the element \'")+place+"\'.");
  if(sxml->CheckNodeActive(node))ReadXml(sxml,node->ToElement());
}

//==============================================================================
/// Reads list of initial conditions in the XML node.
//==============================================================================
void JSphFlexStruc::ReadXml(const JXml *sxml,TiXmlElement* lis){
  //-Loads flexible structure body elements.
  const unsigned idmax=CODE_TYPE_FLEXSTRUCCLAMP_MASK-1;
  TiXmlElement* ele=lis->FirstChildElement("flexstrucbody");
  while(ele){
    if(sxml->CheckElementActive(ele)){
      const unsigned id=GetCount();
      if(id>idmax)Run_Exceptioon("Maximum number of flexible structure bodies has been reached.");
      word mkbound=sxml->GetAttributeWord(ele,"mkbound");
      if(ExistMk(mkbound))Run_Exceptioon(fun::PrintStr("An input already exists for the same mkbound=%u.",mkbound));
      word mkclamp=sxml->GetAttributeWord(ele,"mkclamp");
      float particlevolume=static_cast<float>(Simulate2D? Dp*Dp: Dp*Dp*Dp);
      float density=sxml->ReadElementFloat(ele,"density","value");
      double youngmod=sxml->ReadElementDouble(ele,"youngmod","value");
      double poissratio=sxml->ReadElementDouble(ele,"poissratio","value");
      TpConstitModel constitmodel;
      switch(sxml->ReadElementUnsigned(ele,"constitmodel","value")){
        case 1: constitmodel=CONSTITMODEL_PlaneStrain;  break;
        case 2: constitmodel=CONSTITMODEL_PlaneStress;  break;
        case 3: constitmodel=CONSTITMODEL_SVK;          break;
        default: Run_Exceptioon(fun::PrintStr("Constitutive model for body %u is not valid.",id));
      }
      if(Simulate2D && (constitmodel!=CONSTITMODEL_PlaneStrain && constitmodel!=CONSTITMODEL_PlaneStress))
        Run_Exceptioon(fun::PrintStr("Constitutive model for body %u is not valid for 2D simulations.",id));
      else if(!Simulate2D && constitmodel!=CONSTITMODEL_SVK)
        Run_Exceptioon(fun::PrintStr("Constitutive model for body %u is not valid for 3D simulations.",id));
      float hgfactor=sxml->ReadElementFloat(ele,"hgfactor","value",true,0);
      JSphFlexStrucBody* body=new JSphFlexStrucBody(id,mkbound,mkclamp,particlevolume,density,youngmod,poissratio,constitmodel,hgfactor);
      List.push_back(body);
    }
    ele=ele->NextSiblingElement("flexstrucbody");
  }
}

//==============================================================================
/// Configures BoundCode for each body.
//==============================================================================
void JSphFlexStruc::ConfigBoundCode(const JSphMk *mkinfo){
  for(unsigned c=0;c<GetCount();c++){
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(List[c]->MkBound);
    if(cmk<mkinfo->Size() && (CODE_IsMoving(mkinfo->Mkblock(cmk)->Code))){
      List[c]->ConfigBoundCode(mkinfo->Mkblock(cmk)->Code);
    }
    else Run_Exceptioon(fun::PrintStr("MkBound value for mkbound=%u is not a valid Mk moving boundary.",List[c]->MkBound));
  }
}

//==============================================================================
/// Configures ClampCode for each body.
//==============================================================================
void JSphFlexStruc::ConfigClampCode(const JSphMk *mkinfo){
  for(unsigned c=0;c<GetCount();c++){
    if(ExistMk(List[c]->MkClamp))Run_Exceptioon(fun::PrintStr("MkClamp for mkbound=%u cannot be a flexible structure.",List[c]->MkBound));
    const unsigned cmk=mkinfo->GetMkBlockByMkBound(List[c]->MkClamp);
    if(cmk<mkinfo->Size()){
      List[c]->ConfigClampCode(mkinfo->Mkblock(cmk)->Code);
    }
    else Run_Exceptioon(fun::PrintStr("MkClamp value for mkbound=%u is not a valid Mk boundary.",List[c]->MkBound));
  }
}

//==============================================================================
/// Configures particle codings for flexible structures.
/// Configura codificaciones de partículas para estructuras flexibles.
//==============================================================================
void JSphFlexStruc::ConfigCode(unsigned npb,typecode *code){
  for(unsigned c=0;c<GetCount();c++){
    typecode bcode=List[c]->GetBoundCode();
    for(unsigned p=0;p<npb;p++){
      if(code[p]==bcode)code[p]=typecode(CODE_ToFlexStrucFlex(code[p],List[c]->IdBody));
    }
  }
}

//==============================================================================
/// Get maximum initial sound speed across all flexible structures.
/// Obtenga la máxima velocidad de sonido inicial en todas las estructuras flexibles.
//==============================================================================
double JSphFlexStruc::GetInitialSoundSpeed(){
  double cs0=0.0;
  for(unsigned c=0;c<GetCount();c++)cs0=max(cs0,double(List[c]->GetSoundSpeed()));
  return cs0;
}

//==============================================================================
/// Sets mass density for flexible structures.
/// Establece la densidad de masa para estructuras flexibles.
//==============================================================================
void JSphFlexStruc::SetDensity(unsigned npb,const typecode *code,float *rhos){
  for(unsigned p=0;p<npb;p++)rhos[p]=(CODE_IsFlexStrucFlex(code[p])? GetBody(CODE_GetIbodyFlexStruc(code[p]))->GetDensity(): 0);
}
