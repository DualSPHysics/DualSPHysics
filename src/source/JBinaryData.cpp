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

/// \file JBinaryData.cpp \brief Implements the class \ref JBinaryData.

#include "JBinaryData.h"
#include "Functions.h"

#include <fstream>
#include <cmath>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <algorithm>

using namespace std;

const std::string JBinaryData::CodeItemDef="\nITEM\n";
const std::string JBinaryData::CodeValuesDef="\nVALUES";
const std::string JBinaryData::CodeArrayDef="\nARRAY";

//##############################################################################
//# JBinaryDataDef
//##############################################################################
//==============================================================================
/// Devuelve tipo de datos en texto.
/// Returns data type text.
//==============================================================================
std::string JBinaryDataDef::TypeToStr(TpData type){
  string tx="";
  switch(type){
    case JBinaryDataDef::DatText:     tx="text";     break;
    case JBinaryDataDef::DatBool:     tx="bool";     break;
    case JBinaryDataDef::DatChar:     tx="char";     break;
    case JBinaryDataDef::DatUchar:    tx="uchar";    break;
    case JBinaryDataDef::DatShort:    tx="short";    break;
    case JBinaryDataDef::DatUshort:   tx="ushort";   break;
    case JBinaryDataDef::DatInt:      tx="int";      break;
    case JBinaryDataDef::DatUint:     tx="uint";     break;
    case JBinaryDataDef::DatLlong:    tx="llong";    break;
    case JBinaryDataDef::DatUllong:   tx="ullong";   break;
    case JBinaryDataDef::DatFloat:    tx="float";    break;
    case JBinaryDataDef::DatDouble:   tx="double";   break;
    case JBinaryDataDef::DatInt3:     tx="int3";     break;
    case JBinaryDataDef::DatUint3:    tx="uint3";    break;
    case JBinaryDataDef::DatFloat3:   tx="float3";   break;
    case JBinaryDataDef::DatDouble3:  tx="double3";  break;
  }
  return(tx);
}

//==============================================================================
/// Devuelve tamanho del tipo de datos.
/// Returns size of the data type.
//==============================================================================
size_t JBinaryDataDef::SizeOfType(TpData type){
  size_t ret=0;
  switch(type){
    //case JBinaryDataDef::DatText:
    case JBinaryDataDef::DatBool:     ret=sizeof(int);             break;
    case JBinaryDataDef::DatChar:     ret=sizeof(char);            break;
    case JBinaryDataDef::DatUchar:    ret=sizeof(unsigned char);   break;
    case JBinaryDataDef::DatShort:    ret=sizeof(short);           break;
    case JBinaryDataDef::DatUshort:   ret=sizeof(unsigned short);  break;
    case JBinaryDataDef::DatInt:      ret=sizeof(int);             break;
    case JBinaryDataDef::DatUint:     ret=sizeof(unsigned);        break;
    case JBinaryDataDef::DatLlong:    ret=sizeof(llong);           break;
    case JBinaryDataDef::DatUllong:   ret=sizeof(ullong);          break;
    case JBinaryDataDef::DatFloat:    ret=sizeof(float);           break;
    case JBinaryDataDef::DatDouble:   ret=sizeof(double);          break;
    case JBinaryDataDef::DatInt3:     ret=sizeof(tint3);           break;
    case JBinaryDataDef::DatUint3:    ret=sizeof(tuint3);          break;
    case JBinaryDataDef::DatFloat3:   ret=sizeof(tfloat3);         break;
    case JBinaryDataDef::DatDouble3:  ret=sizeof(tdouble3);        break;
  }
  return(ret);
}

//==============================================================================
/// Devuelve true cuando el tipo es triple.
/// Returns true when the type is triple.
//==============================================================================
bool JBinaryDataDef::TypeIsTriple(TpData type){
  bool ret=false;
  switch(type){
    //case JBinaryDataDef::DatText:
    case JBinaryDataDef::DatBool:
    case JBinaryDataDef::DatChar:
    case JBinaryDataDef::DatUchar:
    case JBinaryDataDef::DatShort:
    case JBinaryDataDef::DatUshort:
    case JBinaryDataDef::DatInt:
    case JBinaryDataDef::DatUint:
    case JBinaryDataDef::DatLlong:
    case JBinaryDataDef::DatUllong:
    case JBinaryDataDef::DatFloat:
    case JBinaryDataDef::DatDouble:   
      ret=false;
    break;
    case JBinaryDataDef::DatInt3:
    case JBinaryDataDef::DatUint3:
    case JBinaryDataDef::DatFloat3:
    case JBinaryDataDef::DatDouble3:
      ret=true;
    break;
  }
  return(ret);
}

//##############################################################################
//# JBinaryDataArray
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JBinaryDataArray::JBinaryDataArray(JBinaryData* parent,const std::string &name,JBinaryDataDef::TpData type):Type(type){ 
  ClassName="JBinaryDataArray";
  Parent=parent;
  Name=name;
  Hide=false;
  Pointer=NULL;
  ExternalPointer=false;
  Count=Size=0;
  ClearFileData();
}

//==============================================================================
/// Destructor.
//==============================================================================
JBinaryDataArray::~JBinaryDataArray(){
  DestructorActive=true;
  FreeMemory();
}

//==============================================================================
/// Devuelve la cantidad de memoria reservada.
/// Returns the amount of memory reserved.
//==============================================================================
llong JBinaryDataArray::GetAllocMemory()const{
  return(Pointer&&!ExternalPointer? Size*JBinaryDataDef::SizeOfType(Type): 0);
}

//==============================================================================
/// Cambia nombre de array comprobando que no existe otro array o value con el mismo nombre.
/// Exception to ensure that there is no duplictaes arrays or values.
//==============================================================================
void JBinaryDataArray::SetName(const std::string &name){
  if(Parent->ExistsValue(name))Run_Exceptioon("There is already a value with the name given.");
  if(Parent->GetArray(name)!=NULL)Run_Exceptioon("There is already an array with the name given.");
  if(Parent->GetItem(name)!=NULL)Run_Exceptioon("There is already an item with the name given.");
  Name=name;
}

//==============================================================================
/// Libera memoria asignada al puntero indicado.
/// Frees memory allocated to the specified pointer.
//==============================================================================
void JBinaryDataArray::FreePointer(void* ptr)const{
  if(ptr)switch(Type){
    case JBinaryDataDef::DatText:     delete[] (string*)ptr;          break;
    case JBinaryDataDef::DatBool:     delete[] (bool*)ptr;            break;
    case JBinaryDataDef::DatInt:      delete[] (int*)ptr;             break;
    case JBinaryDataDef::DatUint:     delete[] (unsigned int*)ptr;    break;
    case JBinaryDataDef::DatChar:     delete[] (char*)ptr;            break;
    case JBinaryDataDef::DatUchar:    delete[] (unsigned char*)ptr;   break;
    case JBinaryDataDef::DatShort:    delete[] (short*)ptr;           break;
    case JBinaryDataDef::DatUshort:   delete[] (unsigned short*)ptr;  break;
    case JBinaryDataDef::DatLlong:    delete[] (llong*)ptr;           break;
    case JBinaryDataDef::DatUllong:   delete[] (ullong*)ptr;          break;
    case JBinaryDataDef::DatFloat:    delete[] (float*)ptr;           break;
    case JBinaryDataDef::DatDouble:   delete[] (double*)ptr;          break;
    case JBinaryDataDef::DatInt3:     delete[] (tint3*)ptr;           break;  
    case JBinaryDataDef::DatUint3:    delete[] (tuint3*)ptr;          break;
    case JBinaryDataDef::DatFloat3:   delete[] (tfloat3*)ptr;         break;
    case JBinaryDataDef::DatDouble3:  delete[] (tdouble3*)ptr;        break;
    default: Run_Exceptioon("Type of array invalid.");
  }
}

//==============================================================================
/// Devuelve puntero con la memoria asiganda.
/// Returns pointer to the allocated memory.
//==============================================================================
void* JBinaryDataArray::AllocPointer(unsigned size)const{
  void* ptr=NULL;
  if(size){
    try{
      switch(Type){
        case JBinaryDataDef::DatText:     ptr=new string[size];          break;
        case JBinaryDataDef::DatBool:     ptr=new bool[size];            break;
        case JBinaryDataDef::DatInt:      ptr=new int[size];             break;
        case JBinaryDataDef::DatUint:     ptr=new unsigned int[size];    break;
        case JBinaryDataDef::DatChar:     ptr=new char[size];            break;
        case JBinaryDataDef::DatUchar:    ptr=new unsigned char[size];   break;
        case JBinaryDataDef::DatShort:    ptr=new short[size];           break;
        case JBinaryDataDef::DatUshort:   ptr=new unsigned short[size];  break;
        case JBinaryDataDef::DatLlong:    ptr=new llong[size];           break;
        case JBinaryDataDef::DatUllong:   ptr=new ullong[size];          break;
        case JBinaryDataDef::DatFloat:    ptr=new float[size];           break;
        case JBinaryDataDef::DatDouble:   ptr=new double[size];          break;
        case JBinaryDataDef::DatInt3:     ptr=new tint3[size];           break;
        case JBinaryDataDef::DatUint3:    ptr=new tuint3[size];          break;
        case JBinaryDataDef::DatFloat3:   ptr=new tfloat3[size];         break;
        case JBinaryDataDef::DatDouble3:  ptr=new tdouble3[size];        break;
        default: Run_Exceptioon("Type of array invalid.");
      }
    }
    catch(const std::bad_alloc){
      Run_Exceptioon("Cannot allocate the requested memory.");
    }
  }
  return(ptr);
}

//==============================================================================
/// Comprueba memoria disponible y redimensiona array si hace falta.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
/// Check available memory array and resize if necessary.
/// If ExternalPointer will not allow to resize the allocated memory.
//==============================================================================
void JBinaryDataArray::CheckMemory(unsigned count,bool resize){
  if(count){
    //-Reserva memoria si fuese necesario.
    //-Allocates memory if necessary.
    if(!Pointer){
      if(!resize)Run_Exceptioon("Memory no allocated.");
      AllocMemory(count);
    }
    if(Count+count>Size){
      if(ExternalPointer)Run_Exceptioon("Allocated memory in external pointer is not enough.");
      if(!resize)Run_Exceptioon("Allocated memory is not enough.");
      AllocMemory(Count+count,true);
    }
  }
}

//==============================================================================
/// Extrae datos del ptr indicado.
/// Extract data from ptr indicated.
//==============================================================================
void JBinaryDataArray::OutData(unsigned &count,unsigned size,const byte *ptr,byte *dat,unsigned sdat)const{
  const unsigned count2=count+sdat;
  if(count2>size)Run_Exceptioon("Overflow in reading data.");
  memcpy(dat,ptr+count,sdat);
  count=count2;
}

//==============================================================================
/// Extrae string de ptr.
/// Extract string ptr.
//==============================================================================
std::string JBinaryDataArray::OutStr(unsigned &count,unsigned size,const byte *ptr)const{
  unsigned len=OutUint(count,size,ptr);
  string tex;
  tex.resize(len);
  const unsigned count2=count+len;
  if(count2>size)Run_Exceptioon("Overflow in reading data.");
  memcpy((char*)tex.c_str(),ptr+count,len);
  count=count2;
  return(tex);
}


//==============================================================================
/// Libera memoria asignada.
/// Frees allocated memory.
//==============================================================================
void JBinaryDataArray::FreeMemory(){
  if(Pointer&&!ExternalPointer)FreePointer(Pointer);
  Pointer=NULL;
  ExternalPointer=false;
  Count=0;
  Size=0;
}

//==============================================================================
/// Asigna memoria para los elementos indicados.
/// Allocate memory for the elements indicated.
//==============================================================================
void JBinaryDataArray::AllocMemory(unsigned size,bool savedata){
  if(Count&&savedata&&size){
    if(ExternalPointer)Run_Exceptioon("External pointer can not be resized.");
    const unsigned count2=min(Count,size);
    void *ptr=AllocPointer(size);
    if(Type==JBinaryDataDef::DatText){//-String array.
      string *strings1=(string*)Pointer;
      string *strings2=(string*)ptr;
      for(unsigned c=0;c<count2;c++)strings2[c]=strings1[c];
    }
    else memcpy((byte*)ptr,(byte*)Pointer,JBinaryDataDef::SizeOfType(Type)*count2);
    FreeMemory();
    Pointer=ptr;
    Count=count2;
    Size=size;
  }
  else{
    FreeMemory();
    Size=size;
    if(Size)Pointer=AllocPointer(Size);
  }
}

//==============================================================================
/// Asigna memoria para los elementos indicados.
/// Allocate memory for the elements indicated.
//==============================================================================
void JBinaryDataArray::ConfigExternalMemory(unsigned size,void* pointer){
  FreeMemory();
  ExternalPointer=true;
  Pointer=pointer;
  Size=size;
  Count=0;
}

//==============================================================================
/// Configura acceso a datos en fichero.
/// Set file data access.
//==============================================================================
void JBinaryDataArray::ConfigFileData(llong filepos,unsigned datacount,unsigned datasize){
  FreeMemory();
  FileDataPos=filepos; FileDataCount=datacount; FileDataSize=datasize;
}

//==============================================================================
/// Borra datos de acceso a datos en fichero.
/// Delete data file data access.
//==============================================================================
void JBinaryDataArray::ClearFileData(){
  FileDataPos=-1; FileDataCount=FileDataSize=0;
}

//==============================================================================
/// Carga contenido de fichero abierto con OpenFileStructure().
/// Load open file content with OpenFileStructure (). 
//==============================================================================
void JBinaryDataArray::ReadFileData(bool resize){
  //printf("ReadFileData Parent_name:[%s] p:%p\n",Parent->GetName().c_str(),Parent);
  //printf("ReadFileData Parent2_name:[%s] p:%p\n",(Parent->GetParent()? Parent->GetParent()->GetName().c_str(): "none"),Parent->GetParent());
  //printf("ReadFileData root_name:[%s] p:%p\n",Parent->GetItemRoot()->GetName().c_str(),Parent->GetItemRoot());
  ifstream *pf=Parent->GetItemRoot()->GetFileStructure();
  if(!pf||!pf->is_open())Run_Exceptioon("The file with data is not available.");
  //printf("ReadFileData[%s]> fpos:%llu count:%u size:%u\n",Name.c_str(),FileDataPos,FileDataCount,FileDataSize);
  if(FileDataPos<0)Run_Exceptioon("The access information to data file is not available.");
  pf->seekg(FileDataPos,ios::beg);
  ReadData(FileDataCount,FileDataSize,pf,resize);
}

//==============================================================================
/// Anhade elementos al array de un fichero.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
/// Add elements to the array of a file. 
/// If ExternalPointer will not allow to resize the allocated memory.
//==============================================================================
void JBinaryDataArray::ReadData(unsigned count,unsigned size,std::ifstream *pf,bool resize){
  if(count){
    //-Reserva memoria si fuese necesario.
    CheckMemory(count,resize);
    //-Carga datos de fichero.
    if(GetType()==JBinaryDataDef::DatText){//-String Array.
      byte *buf=new byte[size];
      pf->read((char*)buf,size);
      unsigned cbuf=0;
      for(unsigned c=0;c<count;c++)AddText(OutStr(cbuf,size,buf),false);
      delete[] buf;
    }
    else{
      const unsigned stype=(unsigned)JBinaryDataDef::SizeOfType(Type);
      const unsigned sdat=stype*count;
      const unsigned cdat=stype*Count;
      pf->read(((char*)Pointer)+cdat,sdat);
      Count+=count;
    }
  }
}

//==============================================================================
/// Anhade elementos al array.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
/// Add elements to the array.
/// If ExternalPointer will not allow to resize the allocated memory.
//==============================================================================
void JBinaryDataArray::AddData(unsigned count,const void* data,bool resize){
  if(count){
    //-Reserva memoria si fuese necesario.
    //-Allocates memory if necessary.
    CheckMemory(count,resize);
    //-Anhade datos al puntero.
    //-Add data to the pointer.
    if(Type==JBinaryDataDef::DatText){
      string *strings=(string*)Pointer;
      string *strings2=(string*)data;
      for(unsigned c=0;c<count;c++)strings[Count+c]=strings2[c];
    }
    else{
      const unsigned stype=(unsigned)JBinaryDataDef::SizeOfType(Type);
      const unsigned sdat=stype*count;
      const unsigned cdat=stype*Count;
      memcpy(((byte*)Pointer)+cdat,(byte*)data,sdat);
    }
    Count+=count;
  }
}

//==============================================================================
/// Guarda datos como contenido del array.
/// Save data as contents of the array.
//==============================================================================
void JBinaryDataArray::SetData(unsigned count,const void* data,bool externalpointer){
  FreeMemory();
  if(externalpointer){
    Pointer=(void*)data;
    ExternalPointer=true;
    Size=count;
    Count=count;
  }
  else AddData(count,data,true);
}

//==============================================================================
/// Anhade un string al array.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
/// Add a string to array.
/// If ExternalPointer will not allow to resize the allocated memory.
//==============================================================================
void JBinaryDataArray::AddText(const std::string &str,bool resize){
  if(Type!=JBinaryDataDef::DatText)Run_Exceptioon("Type of array is not Text.");
  unsigned count=1;
  if(count){
    //-Reserva memoria si fuese necesario.
    //-Allocates memory if necessary.
    CheckMemory(count,resize);
    //-Anhade string al array.
    //-Add string to array.
    string *strings=(string*)Pointer;
    strings[Count]=str;
    Count+=count;
  }
}

//==============================================================================
/// Anhade array de strings al array.
/// Si es ExternalPointer no permite redimensionar la memoria asignada.
/// Add array of strings to array.
/// If ExternalPointer will not allow to resize the allocated memory.
//==============================================================================
void JBinaryDataArray::AddTexts(unsigned count,const std::string *strs,bool resize){
  if(Type!=JBinaryDataDef::DatText)Run_Exceptioon("Type of array is not Text.");
  if(count)AddData(count,(const void*)strs,resize);
}

//==============================================================================
/// Devuelve puntero de datos comprobando que tenga datos.
/// Returns pointer and check is pointer is populated with data.
//==============================================================================
const void* JBinaryDataArray::GetDataPointer()const{
  if(!DataInPointer())Run_Exceptioon("There are not available data in pointer.");
  return(Pointer);
}

//==============================================================================
/// Copia datos de Pointer o FileData al puntero indicado y devuelve el numero
/// de elementos.
/// Copies data of the pointer or FileData indicated pointer and returns
/// the number of elements.
//==============================================================================
unsigned JBinaryDataArray::GetDataCopy(unsigned size,void* pointer)const{
  if(!DataInPointer()&&!DataInFile())Run_Exceptioon("There are not available data in Pointer or FileData.");
  const size_t stype=JBinaryDataDef::SizeOfType(GetType());
  if(!stype)Run_Exceptioon("Type of array is invalid for this function.");
  unsigned count=0;
  if(DataInPointer()){
    count=GetCount();
    if(size>=count)memcpy(pointer,Pointer,stype*count);
  }
  else{
    count=FileDataCount;
    if(size>=count){
      ifstream *pf=Parent->GetItemRoot()->GetFileStructure();
      if(!pf||!pf->is_open())Run_Exceptioon("The file with data is not available.");
      pf->seekg(FileDataPos,ios::beg);
      count=FileDataCount;
      pf->read((char*)pointer,stype*count);
    }
  }
  if(size<count)Run_Exceptioon("Size of array is not enough to store all data.");
  return(count);
}


//##############################################################################
//# JBinaryData
//##############################################################################
//==============================================================================
/// Constructor.
//==============================================================================
JBinaryData::JBinaryData(std::string name):Name(name){ 
  ClassName="JBinaryData";
  Parent=NULL;
  FileStructure=NULL;
  ValuesData=NULL;
  ValuesCacheReset();
  HideAll=HideValues=false;
  FmtFloat="%.7E";
  FmtDouble="%.15E";
}

//==============================================================================
/// Constructor de copias.
/// Copu of constructor
//==============================================================================
JBinaryData::JBinaryData(const JBinaryData &src){
  ClassName="JBinaryData";
  Parent=NULL;
  FileStructure=NULL;
  ValuesData=NULL;
  ValuesCacheReset();
  *this=src;
}

//==============================================================================
/// Destructor.
//==============================================================================
JBinaryData::~JBinaryData(){
  DestructorActive=true;
  Clear();
}

//==============================================================================
/// Sobrecarga del operador de asignacion.
/// Overload assignment operator.
//==============================================================================
JBinaryData& JBinaryData::operator=(const JBinaryData &src){
  if(this!=&src){
    unsigned size=src.GetSizeDataConst(true);
    byte *dat=new byte[size];
    src.SaveDataConst(size,dat,true);
    LoadData(size,dat);
    delete[] dat;
  }
  return(*this);
}

//==============================================================================
/// Elimina todo el contenido (values, arrays e items).
/// Deletes all contents (values, arrays and items).
//==============================================================================
void JBinaryData::Clear(){
  RemoveValues();
  RemoveArrays();
  RemoveItems();
  ValuesCacheReset();
  CloseFileStructure();
}

//==============================================================================
/// Devuelve la cantidad de memoria reservada por el objeto y subobjetos.
/// Returns the amount of memory allocated by the object and sub-objects.
//==============================================================================
llong JBinaryData::GetAllocMemory()const{
  llong s=0; 
  for(unsigned c=0;c<Arrays.size();c++)s+=Arrays[c]->GetAllocMemory(); //-Memoria de arrays no externos.
  s+=ValuesSize;                                                       //-Memoria para cache de values.
  for(unsigned c=0;c<Items.size();c++)s+=Items[c]->GetAllocMemory();   //-Memoria de Items descencientes.
  return(s);
}

//==============================================================================
/// Elimina cache de values.
/// Removes cache values.
//==============================================================================
void JBinaryData::ValuesCacheReset(){
  delete[] ValuesData; ValuesData=NULL;
  ValuesSize=0;
  ValuesModif=true;
}

//==============================================================================
/// Prepara cache de values.
/// Prepare cache values.
//==============================================================================
void JBinaryData::ValuesCachePrepare(bool down){
  if(ValuesModif){
    ValuesCacheReset();
    ValuesSize=GetSizeValues();
    ValuesData=new byte[ValuesSize];
    unsigned count=0;
    SaveValues(count,ValuesSize,ValuesData);
    ValuesModif=false;
  }
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->ValuesCachePrepare(true);
}

//==============================================================================
/// Comprueba existencia y tipo de valor. Devuelve posicion de valor (-1 no exite).
/// Check existence of value and type. Returns position (index) value (-1 if not exist).
//==============================================================================
int JBinaryData::CheckGetValue(const std::string &name,bool optional,JBinaryDataDef::TpData type)const{
  int idx=GetValueIndex(name);
  //if(idx<0&&GetArrayIndex(name)>=0)Run_Exceptioon(string("The value ")+name+" is an array.");
  if(!optional&&idx<0)Run_Exceptioon(string("Value ")+name+" not found.");
  if(idx>=0&&Values[idx].type!=type)Run_Exceptioon(string("Type of value ")+name+" invalid.");
  return(idx);
}

//==============================================================================
/// Comprueba tipo de valor, y sino existe lo crea. Devuelve posicion de valor.
/// Check value type, create if it does not exists, Returns position value.
//==============================================================================
int JBinaryData::CheckSetValue(const std::string &name,JBinaryDataDef::TpData type){
  int idx=GetValueIndex(name);
  if(idx<0&&GetArray(name)!=NULL)Run_Exceptioon(string("The value ")+name+" is an array.");
  if(idx<0&&GetItem(name)!=NULL)Run_Exceptioon(string("The value ")+name+" is an item.");
  if(idx>=0&&Values[idx].type!=type)Run_Exceptioon(string("Type of value ")+name+" invalid.");
  if(idx<0){
    StValue v;
    if(name.length()>120)Run_Exceptioon(string("The name of value ")+name+"  is too large.");
    if(name.empty())Run_Exceptioon(string("The name of value ")+name+"  is empty.");
    ResetValue(name,type,v);
    Values.push_back(v);
    idx=int(Values.size())-1;
  }
  ValuesModif=true;
  return(idx);
}

//==============================================================================
/// Reset del valor pasado como referencia.
/// Reset the value passed by reference.
//==============================================================================
void JBinaryData::ResetValue(const std::string &name,JBinaryDataDef::TpData type,JBinaryData::StValue &v){
  v.name=name; v.type=type; v.vdouble3=TDouble3(0);
}

//==============================================================================
/// Devuelve Value en formato XML.
/// Value returned in XML format.
//==============================================================================
std::string JBinaryData::ValueToXml(const StValue &v)const{
  string tx=JBinaryDataDef::TypeToStr(v.type);
  if(tx.empty())Run_Exceptioon("Name of type invalid.");
  tx=string("<")+tx+" name=\""+v.name+"\" ";
  switch(v.type){
    case JBinaryDataDef::DatText:     tx=tx+"v=\""+ v.vtext +"\" />";                         break;
    case JBinaryDataDef::DatBool:     tx=tx+"v=\""+ (v.vint? "1": "0") +"\" />";              break;    
    case JBinaryDataDef::DatChar:     tx=tx+"v=\""+ fun::IntStr(v.vchar)  +"\" />";           break;    
    case JBinaryDataDef::DatUchar:    tx=tx+"v=\""+ fun::UintStr(v.vuchar)  +"\" />";         break;    
    case JBinaryDataDef::DatShort:    tx=tx+"v=\""+ fun::IntStr(v.vshort) +"\" />";           break;
    case JBinaryDataDef::DatUshort:   tx=tx+"v=\""+ fun::UintStr(v.vushort) +"\" />";         break;
    case JBinaryDataDef::DatInt:      tx=tx+"v=\""+ fun::IntStr(v.vint) +"\" />";             break;
    case JBinaryDataDef::DatUint:     tx=tx+"v=\""+ fun::UintStr(v.vuint) +"\" />";           break;
    case JBinaryDataDef::DatLlong:    tx=tx+"v=\""+ fun::LongStr(v.vllong) +"\" />";          break;
    case JBinaryDataDef::DatUllong:   tx=tx+"v=\""+ fun::UlongStr(v.vullong) +"\" />";        break;
    case JBinaryDataDef::DatFloat:    tx=tx+"v=\""+ fun::FloatStr(v.vfloat,FmtFloat.c_str()) +"\" />";     break;
    case JBinaryDataDef::DatDouble:   tx=tx+"v=\""+ fun::DoubleStr(v.vdouble,FmtDouble.c_str()) +"\" />";  break;
    case JBinaryDataDef::DatInt3:     tx=tx+"x=\""+ fun::IntStr(v.vint3.x) +"\" y=\""+ fun::IntStr(v.vint3.y) +"\" z=\""+ fun::IntStr(v.vint3.z) +"\" />";                                   break;
    case JBinaryDataDef::DatUint3:    tx=tx+"x=\""+ fun::UintStr(v.vuint3.x) +"\" y=\""+ fun::UintStr(v.vuint3.y) +"\" z=\""+ fun::UintStr(v.vuint3.z) +"\" />";                             break;
    case JBinaryDataDef::DatFloat3:   tx=tx+"x=\""+ fun::FloatStr(v.vfloat3.x,FmtFloat.c_str()) +"\" y=\""+ fun::FloatStr(v.vfloat3.y,FmtFloat.c_str()) +"\" z=\""+ fun::FloatStr(v.vfloat3.z,FmtFloat.c_str()) +"\" />";           break;  //"%.7E"
    case JBinaryDataDef::DatDouble3:  tx=tx+"x=\""+ fun::DoubleStr(v.vdouble3.x,FmtDouble.c_str()) +"\" y=\""+ fun::DoubleStr(v.vdouble3.y,FmtDouble.c_str()) +"\" z=\""+ fun::DoubleStr(v.vdouble3.z,FmtDouble.c_str()) +"\" />";  break;  //"%.15E"
    default: Run_Exceptioon("Type of value invalid.");
  }
  return(tx);
}



//==============================================================================
/// Extrae datos del ptr indicado.
/// Extract data from ptr indicated.
//==============================================================================
void JBinaryData::OutData(unsigned &count,unsigned size,const byte *ptr,byte *dat,unsigned sdat)const{
  const unsigned count2=count+sdat;
  if(count2>size)Run_Exceptioon("Overflow in reading data.");
  memcpy(dat,ptr+count,sdat);
  count=count2;
}

//==============================================================================
/// Extrae string de ptr.
/// Remove string ptr.
//==============================================================================
std::string JBinaryData::OutStr(unsigned &count,unsigned size,const byte *ptr)const{
  unsigned len=OutUint(count,size,ptr);
  string tex;
  tex.resize(len);
  const unsigned count2=count+len;
  if(count2>size)Run_Exceptioon("Overflow in reading data.");
  memcpy((char*)tex.c_str(),ptr+count,len);
  count=count2;
  return(tex);
}

//==============================================================================
/// Introduce datos en ptr.
/// Put data in ptr.
//==============================================================================
void JBinaryData::InData(unsigned &count,unsigned size,byte *ptr,const byte *dat,unsigned sdat)const{
  if(ptr){
    if(count+sdat>size)Run_Exceptioon("Insufficient memory for data.");
    memcpy(ptr+count,dat,sdat);
  }
  if(count+sdat<count)Run_Exceptioon("Size of data is too huge.");
  count+=sdat;
}

//==============================================================================
/// Introduce string en ptr.
/// Put string in ptr.
//==============================================================================
void JBinaryData::InStr(unsigned &count,unsigned size,byte *ptr,const std::string &cad)const{
  InUint(count,size,ptr,(unsigned)cad.length());
  InData(count,size,ptr,(byte*)cad.c_str(),unsigned(cad.length()));
}

//==============================================================================
/// Introduce Value en ptr.
/// Put Value in ptr.
//==============================================================================
void JBinaryData::InValue(unsigned &count,unsigned size,byte *ptr,const StValue &v)const{
  InStr(count,size,ptr,v.name);
  InInt(count,size,ptr,int(v.type));
  switch(v.type){
    case JBinaryDataDef::DatText:     InStr    (count,size,ptr,v.vtext);     break;
    case JBinaryDataDef::DatBool:     InBool   (count,size,ptr,v.vint!=0);   break;    
    case JBinaryDataDef::DatChar:     InChar   (count,size,ptr,v.vchar);     break;    
    case JBinaryDataDef::DatUchar:    InUchar  (count,size,ptr,v.vuchar);    break;    
    case JBinaryDataDef::DatShort:    InShort  (count,size,ptr,v.vshort);    break;    
    case JBinaryDataDef::DatUshort:   InUshort (count,size,ptr,v.vushort);   break;    
    case JBinaryDataDef::DatInt:      InInt    (count,size,ptr,v.vint);      break;    
    case JBinaryDataDef::DatUint:     InUint   (count,size,ptr,v.vuint);     break;    
    case JBinaryDataDef::DatLlong:    InLlong  (count,size,ptr,v.vllong);    break;     
    case JBinaryDataDef::DatUllong:   InUllong (count,size,ptr,v.vullong);   break;    
    case JBinaryDataDef::DatFloat:    InFloat  (count,size,ptr,v.vfloat);    break;    
    case JBinaryDataDef::DatDouble:   InDouble (count,size,ptr,v.vdouble);   break;    
    case JBinaryDataDef::DatInt3:     InInt3   (count,size,ptr,v.vint3);     break;    
    case JBinaryDataDef::DatUint3:    InUint3  (count,size,ptr,v.vuint3);    break;    
    case JBinaryDataDef::DatFloat3:   InFloat3 (count,size,ptr,v.vfloat3);   break;    
    case JBinaryDataDef::DatDouble3:  InDouble3(count,size,ptr,v.vdouble3);  break;    
    default: Run_Exceptioon("Type of value invalid.");
  }
}

//==============================================================================
/// Extrae Value en ptr.
/// Extract Value in ptr.
//==============================================================================
void JBinaryData::OutValue(unsigned &count,unsigned size,const byte *ptr){
  string name=OutStr(count,size,ptr);
  JBinaryDataDef::TpData type=(JBinaryDataDef::TpData)OutInt(count,size,ptr);
  switch(type){
    case JBinaryDataDef::DatText:     SetvText   (name,OutStr    (count,size,ptr));   break;
    case JBinaryDataDef::DatBool:     SetvBool   (name,OutBool   (count,size,ptr));   break;
    case JBinaryDataDef::DatChar:     SetvChar   (name,OutChar   (count,size,ptr));   break;
    case JBinaryDataDef::DatUchar:    SetvUchar  (name,OutUchar  (count,size,ptr));   break;
    case JBinaryDataDef::DatShort:    SetvShort  (name,OutShort  (count,size,ptr));   break;
    case JBinaryDataDef::DatUshort:   SetvUshort (name,OutUshort (count,size,ptr));   break;
    case JBinaryDataDef::DatInt:      SetvInt    (name,OutInt    (count,size,ptr));   break;
    case JBinaryDataDef::DatUint:     SetvUint   (name,OutUint   (count,size,ptr));   break;
    case JBinaryDataDef::DatLlong:    SetvLlong  (name,OutLlong  (count,size,ptr));   break;
    case JBinaryDataDef::DatUllong:   SetvUllong (name,OutUllong (count,size,ptr));   break;
    case JBinaryDataDef::DatFloat:    SetvFloat  (name,OutFloat  (count,size,ptr));   break;
    case JBinaryDataDef::DatDouble:   SetvDouble (name,OutDouble (count,size,ptr));   break;
    case JBinaryDataDef::DatInt3:     SetvInt3   (name,OutInt3   (count,size,ptr));   break;
    case JBinaryDataDef::DatUint3:    SetvUint3  (name,OutUint3  (count,size,ptr));   break;
    case JBinaryDataDef::DatFloat3:   SetvFloat3 (name,OutFloat3 (count,size,ptr));   break;
    case JBinaryDataDef::DatDouble3:  SetvDouble3(name,OutDouble3(count,size,ptr));   break;
    default: Run_Exceptioon("Type of value invalid.");
  }
}

//==============================================================================
/// Introduce datos basicos de Array en ptr.
/// Put basic data Array in ptr.
//==============================================================================
void JBinaryData::InArrayBase(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const{
  InStr(count,size,ptr,CodeArrayDef);
  InStr(count,size,ptr,ar->GetName());
  InBool(count,size,ptr,ar->GetHide());
  InInt(count,size,ptr,int(ar->GetType()));
  InUint(count,size,ptr,ar->GetCount());
  //-Calcula e introduce size de los datos del array.
  unsigned sizearraydata=0;
  InArrayData(sizearraydata,0,NULL,ar);
  InUint(count,size,ptr,sizearraydata);
}
//==============================================================================
/// Introduce contendido de Array en ptr.
/// Put ptr Array content.
//==============================================================================
void JBinaryData::InArrayData(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const{
  const JBinaryDataDef::TpData type=ar->GetType();
  const unsigned num=ar->GetCount();
  const void* pointer=ar->GetPointer();
  if(num&&!pointer)Run_Exceptioon("Pointer of array with data is invalid.");
  if(type==JBinaryDataDef::DatText){//-Array de strings.
    const string *list=(string*)pointer;
    for(unsigned c=0;c<num;c++)InStr(count,size,ptr,list[c]);
  }
  else{//-Array de tipos basicos.
    unsigned sizetype=(unsigned)JBinaryDataDef::SizeOfType(ar->GetType());
    InData(count,size,ptr,(byte*)pointer,sizetype*num);
  }
}

//==============================================================================
/// Introduce Array en ptr.
/// Put Array in ptr.
//==============================================================================
void JBinaryData::InArray(unsigned &count,unsigned size,byte *ptr,const JBinaryDataArray *ar)const{
  //-Calcula size de la definicion del array.
  unsigned sizearray=0;
  InArrayBase(sizearray,0,NULL,ar);
  //-Introduce propiedades de array en ptr.
  InUint(count,size,ptr,sizearray);
  InArrayBase(count,size,ptr,ar);
  //-Introduce contenido del array.
  InArrayData(count,size,ptr,ar);
}

//==============================================================================
/// Introduce datos basicos de Item en ptr.
/// Put basic data Item in ptr.
//==============================================================================
void JBinaryData::InItemBase(unsigned &count,unsigned size,byte *ptr,bool all)const{
  InStr(count,size,ptr,CodeItemDef);
  InStr(count,size,ptr,GetName());
  InBool(count,size,ptr,GetHide());     
  InBool(count,size,ptr,GetHideValues());
  InStr(count,size,ptr,GetFmtFloat());
  InStr(count,size,ptr,GetFmtDouble());
  InUint(count,size,ptr,(all? GetArraysCount(): GetVisibleArraysCount()));
  InUint(count,size,ptr,(all? GetItemsCount(): GetVisibleItemsCount()));
  if(all||!HideValues)InUint(count,size,ptr,(ValuesModif? GetSizeValues(): ValuesSize));
  else InUint(count,size,ptr,0);
}

//==============================================================================
/// Introduce Item en ptr.
/// Put Item in ptr.
//==============================================================================
void JBinaryData::InItem(unsigned &count,unsigned size,byte *ptr,bool all)const{
  //-Calcula size de la definicion del item.
  unsigned sizeitem=0;
  InItemBase(sizeitem,0,NULL,all);
  //-Introduce propiedades de item en ptr.
  InUint(count,size,ptr,sizeitem);
  InItemBase(count,size,ptr,all);
  //-Introduce values en ptr.
  if(all||!HideValues){
    if(ValuesModif)SaveValues(count,size,ptr);//-Cache no valida.
    else InData(count,size,ptr,ValuesData,ValuesSize);//-Cache actualizada.
  }
  //-Introduce arrays en ptr.
  for(unsigned c=0;c<Arrays.size();c++)if(all||!Arrays[c]->GetHide())InArray(count,size,ptr,Arrays[c]);
  //-Introduce items en ptr.
  for(unsigned c=0;c<Items.size();c++)if(all||!Items[c]->GetHide())Items[c]->InItem(count,size,ptr,all);
}

//==============================================================================
/// Extrae datos basicos del Array de ptr.
/// Extract basic data from ptr Array 
//==============================================================================
JBinaryDataArray* JBinaryData::OutArrayBase(unsigned &count,unsigned size,const byte *ptr,unsigned &countdata,unsigned &sizedata){
  if(OutStr(count,size,ptr)!=CodeArrayDef)Run_Exceptioon("Validation code is invalid.");
  string name=OutStr(count,size,ptr);
  bool hide=OutBool(count,size,ptr);
  JBinaryDataDef::TpData type=(JBinaryDataDef::TpData)OutInt(count,size,ptr);
  countdata=OutUint(count,size,ptr);
  sizedata=OutUint(count,size,ptr);
  if(type!=JBinaryDataDef::DatText&&sizedata!=JBinaryDataDef::SizeOfType(type)*countdata)Run_Exceptioon("Size of data is invalid.");
  //-Crea array.
  JBinaryDataArray *ar=CreateArray(name,type);
  ar->SetHide(hide);
  return(ar);
}

//==============================================================================
/// Extrae contenido de Array de ptr.
/// Extract the contents of the ptr Array
//==============================================================================
void JBinaryData::OutArrayData(unsigned &count,unsigned size,const byte *ptr,JBinaryDataArray *ar,unsigned countdata,unsigned sizedata){
  if(ar->GetType()==JBinaryDataDef::DatText){//-Array de strings.
    ar->AllocMemory(countdata);
    for(unsigned c=0;c<countdata;c++)ar->AddText(OutStr(count,size,ptr),false);
  }
  else{
    //-Comprueba que los datos del array estan disponibles.
    //-Checks that the data array is available.
    unsigned count2=count+sizedata;
    if(count2>size)Run_Exceptioon("Overflow in reading data.");
    //-Extrae datos para el array.
    //-Extracts the data for the array.
    ar->AddData(countdata,ptr+count,true);
    count=count2;
  }
}

//==============================================================================
/// Extrae Array de ptr.
/// Extracts the Array of the ptr.
//==============================================================================
void JBinaryData::OutArray(unsigned &count,unsigned size,const byte *ptr){
  //-Crea y configura array a partir de ptr.
  //-Creates and configures array from ptr 
  const unsigned sizearraydef=OutUint(count,size,ptr);
  unsigned countdata,sizedata;
  JBinaryDataArray *ar=OutArrayBase(count,size,ptr,countdata,sizedata);
  //-Extrae contenido del array.
  //-Extract contents of the array.
  OutArrayData(count,size,ptr,ar,countdata,sizedata);
}

//==============================================================================
/// Extrae propiedades basicas de Item de ptr.
/// Extracts basic properties from the ptr Item
//==============================================================================
JBinaryData* JBinaryData::OutItemBase(unsigned &count,unsigned size,const byte *ptr,bool create,unsigned &narrays,unsigned &nitems,unsigned &sizevalues){
  if(OutStr(count,size,ptr)!=CodeItemDef)Run_Exceptioon("Validation code is invalid.");
  JBinaryData* item=this;
  if(create)item=CreateItem(OutStr(count,size,ptr));
  else item->SetName(OutStr(count,size,ptr));
  item->SetHide(OutBool(count,size,ptr));
  item->SetHideValues(OutBool(count,size,ptr),false);
  item->SetFmtFloat(OutStr(count,size,ptr),false);
  item->SetFmtDouble(OutStr(count,size,ptr),false);
  narrays=OutUint(count,size,ptr);
  nitems=OutUint(count,size,ptr);
  sizevalues=OutUint(count,size,ptr);
  return(item);
}

//==============================================================================
/// Extrae Item de ptr.
/// Extracts the Item ptr.
//==============================================================================
void JBinaryData::OutItem(unsigned &count,unsigned size,const byte *ptr,bool create){
  //-Extrae propiedades del item.
  //-Extract item properties
  const unsigned sizeitemdef=OutUint(count,size,ptr);
  unsigned narrays,nitems,sizevalues;
  JBinaryData* item=OutItemBase(count,size,ptr,create,narrays,nitems,sizevalues);
  //-Extrae values del item.
  //-Extract values of the item.
  if(sizevalues){
    if(OutStr(count,size,ptr)!=CodeValuesDef)Run_Exceptioon("Validation code is invalid.");
    unsigned num=OutUint(count,size,ptr);
    for(unsigned c=0;c<num;c++)item->OutValue(count,size,ptr);
  }
  //-Extrae arrays del item.
  //-Extract arrays from item.
  for(unsigned c=0;c<narrays;c++)item->OutArray(count,size,ptr);
  //-Extrae items del item.
  //-Extract items from item.
  for(unsigned c=0;c<nitems;c++)item->OutItem(count,size,ptr,true);
}

//==============================================================================
/// Devuelve el tamanho necesario para almacenar todos los values del item.
/// Returns the size necessary to store all the values in the item 
//==============================================================================
unsigned JBinaryData::GetSizeValues()const{
  unsigned count=0;
  SaveValues(count,0,NULL);
  return(count);
}

//==============================================================================
/// Almacena datos de values en ptr.
/// Data values stored in ptr.
//==============================================================================
void JBinaryData::SaveValues(unsigned &count,unsigned size,byte *ptr)const{
  unsigned num=unsigned(Values.size());
  InStr(count,size,ptr,CodeValuesDef);
  InUint(count,size,ptr,num);
  for(unsigned c=0;c<num;c++)InValue(count,size,ptr,Values[c]);
}



//==============================================================================
/// Graba contenido de array en fichero.
/// Saves the contents of array into the file.
//==============================================================================
void JBinaryData::WriteArrayData(std::fstream *pf,const JBinaryDataArray *ar)const{
  const JBinaryDataDef::TpData type=ar->GetType();
  const unsigned countdata=ar->GetCount();
  const void* pointer=ar->GetPointer();
  if(countdata&&!pointer)Run_Exceptioon("Pointer of array with data is invalid.");
  if(type==JBinaryDataDef::DatText){//-Array de strings. Stings Array
    const string *list=(string*)pointer;
    unsigned sbuf=0;
    for(unsigned c=0;c<countdata;c++)InStr(sbuf,0,NULL,list[c]);//-Calcula size de buffer. Calculate buffer size.
    byte *buf=new byte[sbuf];
    unsigned cbuf=0;
    for(unsigned c=0;c<countdata;c++)InStr(cbuf,sbuf,buf,list[c]);//-Copia en buffer. Copy buffer.
    pf->write((char*)buf,cbuf);
    delete[] buf;
  }
  else{//-Array de tipos basicos. Array of basic types.
    unsigned sizetype=(unsigned)JBinaryDataDef::SizeOfType(ar->GetType());
    pf->write((char*)pointer,sizetype*countdata);
  }
}

//==============================================================================
/// Graba Array en fichero.
/// Saves the Array in the file. 
//==============================================================================
void JBinaryData::WriteArray(std::fstream *pf,unsigned sbuf,byte *buf,const JBinaryDataArray *ar)const{
  //-Calcula size de la definicion del array.
  unsigned sizearray=0;
  InArrayBase(sizearray,0,NULL,ar);
  //-Graba propiedades de array. Saves properties of array.
  unsigned cbuf=0;
  InUint(cbuf,sbuf,buf,sizearray);
  InArrayBase(cbuf,sbuf,buf,ar);
  pf->write((char*)buf,cbuf);
  //-Graba contenido del array. Saves contents of array.
  WriteArrayData(pf,ar);
}

//==============================================================================
/// Graba Item en fichero.
/// Saves items to file.
//==============================================================================
void JBinaryData::WriteItem(std::fstream *pf,unsigned sbuf,byte *buf,bool all)const{
  //-Calcula size de la definicion del item.
  //-Calculates the size of the item's definition.
  unsigned sizeitem=0;
  InItemBase(sizeitem,0,NULL,all);
  //-Graba propiedades de item en fichero.
  //-Saves item properties in the file.
  {
    unsigned cbuf=0;
    InUint(cbuf,sbuf,buf,sizeitem);
    InItemBase(cbuf,sbuf,buf,all);
    pf->write((char*)buf,cbuf);
  }
  //-Graba values. 
  //-Save values
  if(all||!GetHideValues())pf->write((char*)ValuesData,ValuesSize);
  //-Graba arrays.
  //-Save arrays.
  for(unsigned c=0;c<Arrays.size();c++)if(all||!Arrays[c]->GetHide())WriteArray(pf,sbuf,buf,Arrays[c]);
  //-Graba items.
  //-Save items.
  for(unsigned c=0;c<Items.size();c++)if(all||!Items[c]->GetHide())Items[c]->WriteItem(pf,sbuf,buf,all);
}


//==============================================================================
/// Devuelve unsigned leido de fichero.
/// Returns the unsigned value read from the file.
//==============================================================================
unsigned JBinaryData::ReadUint(std::ifstream *pf)const{
  unsigned v=0;
  pf->read((char*)&v,sizeof(unsigned));
  return(v);
}

//==============================================================================
/// Carga datos de array de fichero.
/// Loads data to array from the file.
//==============================================================================
void JBinaryData::ReadArrayData(std::ifstream *pf,JBinaryDataArray *ar,unsigned countdata,unsigned sizedata,bool loadarraysdata){
  const JBinaryDataDef::TpData type=ar->GetType();
  if(loadarraysdata)ar->ReadData(countdata,sizedata,pf,true);
  else{
    ar->ConfigFileData((llong)pf->tellg(),countdata,sizedata);  
    pf->seekg(sizedata,ios::cur);
  }
}

//==============================================================================
/// Carga array de fichero.
/// Loads array from the file.
//==============================================================================
void JBinaryData::ReadArray(std::ifstream *pf,unsigned sbuf,byte *buf,bool loadarraysdata){
  //-Carga propiedades del array. 
  //-Load array properties.
  const unsigned sizearraydef=ReadUint(pf);
  pf->read((char*)buf,sizearraydef);
  unsigned countdata,sizedata;
  unsigned cbuf=0;
  JBinaryDataArray *ar=OutArrayBase(cbuf,sizearraydef,buf,countdata,sizedata);
  //-Extrae contenido del array.
  //-Extract contents of the array.
  ReadArrayData(pf,ar,countdata,sizedata,loadarraysdata);
}

//==============================================================================
/// Carga Item de fichero.
/// Loads item from the file.
//==============================================================================
void JBinaryData::ReadItem(std::ifstream *pf,unsigned sbuf,byte *buf,bool create,bool loadarraysdata){
  //-Carga propiedades del item.
  //-Load item properties.
  const unsigned sizeitemdef=ReadUint(pf);
  pf->read((char*)buf,sizeitemdef);
  unsigned narrays,nitems,sizevalues;
  unsigned cbuf=0;
  JBinaryData* item=OutItemBase(cbuf,sizeitemdef,buf,create,narrays,nitems,sizevalues);
  //-Carga values del item.
  //-Loading item values.
  if(sizevalues){
    byte* buf2=(sbuf>=sizevalues? buf: NULL);
    if(!buf2)buf2=new byte[sizevalues];
    pf->read((char*)buf2,sizevalues);
    unsigned cbuf2=0;
    if(OutStr(cbuf2,sizevalues,buf2)!=CodeValuesDef)Run_Exceptioon("Validation code is invalid.");
    unsigned num=OutUint(cbuf2,sizevalues,buf2);
    for(unsigned c=0;c<num;c++)item->OutValue(cbuf2,sizevalues,buf2);
    if(buf!=buf2)delete[] buf2;
  }
  //-Carga arrays del item.
  // Load arrays from the item
  for(unsigned c=0;c<narrays;c++)item->ReadArray(pf,sbuf,buf,loadarraysdata);
  //-Carga items del item.
  //-Loads items from the item.
  for(unsigned c=0;c<nitems;c++)item->ReadItem(pf,sbuf,buf,true,loadarraysdata);
}

//==============================================================================
/// Genera cabecera para fichero.
/// Generate header for file.
//==============================================================================
JBinaryData::StHeadFmtBin JBinaryData::MakeFileHead(const std::string &filecode)const{
  StHeadFmtBin hfmt; 
  memset(&hfmt,0,sizeof(StHeadFmtBin));
  string titu=string("#FileJBD ")+filecode;
  unsigned stitu=min(58u,unsigned(titu.size()));
  for(unsigned c=0;c<stitu;c++)hfmt.titu[c]=titu[c];
  for(unsigned c=stitu;c<58;c++)hfmt.titu[c]=' ';
  hfmt.titu[58]='\n';
  hfmt.byteorder=byte(fun::GetByteOrder());
  return(hfmt);
}

//==============================================================================
/// Devuelve tamanho de fichero y su cabecera.
/// Si el fichero no contiene una cabecera devuelve 0.
/// Returns file size and its header.
/// If the file does not contain a header returns 0.
//==============================================================================
unsigned JBinaryData::GetFileHead(std::ifstream *pf,JBinaryData::StHeadFmtBin &head)const{
  //-Obtiene size del fichero.
  //-Gets file size.
  pf->seekg(0,ios::end);
  const unsigned fsize=(unsigned)pf->tellg();
  pf->seekg(0,ios::beg);
  //-Lee cabecera basica.
  //-Reads basic header.
  if(fsize>=sizeof(StHeadFmtBin))pf->read((char*)&head,sizeof(StHeadFmtBin));
  else memset(&head,0,sizeof(StHeadFmtBin));
  return(fsize);
}

//==============================================================================
/// Comprueba formato de cabecera con filecode y bitorder.
/// Genera excepcion en caso de error.
/// Check the format of the file header and bitorder.
/// Generates exception on error.
//==============================================================================
void JBinaryData::CheckHead(const std::string &file,const StHeadFmtBin &head,const std::string &filecode)const{
  int err=0;
  //-Coprueba formato de cabecera y filecode.
  //-Check header format and filecode.
  if(!err){
    StHeadFmtBin head2=MakeFileHead(filecode);
    unsigned c=0;
    for(;head.titu[c]==head2.titu[c]&&c<60;c++);
    if(c<9)err=2;
    //else if(!filecode.empty()&&c<60)err=3;
    else if(!filecode.empty()&&c<50)err=3;    //<----- Preliminary solution for strange error on some Linux...
  }
  //-Coprueba orden de bytes.
  //-Check byte order.
  if(!err){
    byte byteorder=byte(fun::GetByteOrder());
    if(head.byteorder!=byte(fun::GetByteOrder()))err=1;
  }
  if(err==1)Run_ExceptioonFile("The byte-order in file is invalid.",file);
  else if(err==2)Run_ExceptioonFile("The format file JBinaryData is invalid.",file);
  else if(err==3)Run_ExceptioonFile("The file code is invalid.",file);
}

//==============================================================================
/// Comprueba formato de cabecera con filecode y bitorder.
/// Si el fichero este vacio tambien genera excepcion.
/// Devuelve size de fichero.
/// Check with filecode header format and bitorder.
/// If the file is empty, it generates exception.
/// Returns file size.
//==============================================================================
unsigned JBinaryData::CheckFileHead(const std::string &file,std::ifstream *pf,const std::string &filecode)const{
  JBinaryData::StHeadFmtBin head;
  //-Obtiene size y cabecera del fichero.
  //-Get size and file header.
  const unsigned fsize=GetFileHead(pf,head);
  //-Comprueba validez de cabecera.
  //-Check for valid header.
  CheckHead(file,head,filecode);
  return(fsize);
}

//==============================================================================
/// Comprueba formato de cabecera con filecode y bitorder.
/// En caso de que el fichero este vacio no genera excepcion.
/// Devuelve size de fichero.
/// Check with filecode header format and bitorder.
/// If the file is empty generates no exception.
/// Returns file size.
//==============================================================================
unsigned JBinaryData::CheckFileListHead(const std::string &file,std::fstream *pf,const std::string &filecode)const{
  //-Obtiene size del fichero.
  //-Gets file size.
  pf->seekg(0,ios::end);
  const unsigned fsize=(unsigned)pf->tellg();   //printf("CheckFileHead> FileSize:%u\n",fsize);
  pf->seekg(0,ios::beg);
  //-Lee cabecera basica y comprueba validez.
  //-Reads basic header and checks validity.
  StHeadFmtBin head;
  if(fsize>=sizeof(StHeadFmtBin)){
    pf->read((char*)&head,sizeof(StHeadFmtBin));
    //-Check for valid header.
    CheckHead(file,head,filecode);
  }
  return(fsize);
}

//==============================================================================
/// Graba contenido en fichero XML.
/// Saves contents to XML file.
//==============================================================================
void JBinaryData::WriteFileXmlArray(const std::string &tabs,std::ofstream* pf,bool svarrays,const JBinaryDataArray* ar)const{
  const JBinaryDataDef::TpData type=ar->GetType();
  const unsigned size=ar->GetSize();
  const unsigned count=ar->GetCount();
  string tx=JBinaryDataDef::TypeToStr(type);
  if(tx.empty())Run_Exceptioon("Name of type invalid.");
  string res=string("<array_")+tx+" name=\""+ar->GetName()+"\" size=\""+fun::UintStr(size)+"\" count=\""+fun::UintStr(count)+"\" hide=\""+ (ar->GetHide()? '1': '0') +"\"";
  if(!svarrays){
    res=res+"/>";
    (*pf) << tabs << res << endl;
  }
  else{
    res=res+">";
    (*pf) << tabs << res << endl;
    const void *data=ar->GetDataPointer();
    JBinaryData::StValue v;
    //ResetValue("",type,v);
    v.type=type;
    for(unsigned c=0;c<count;c++){
      v.name=fun::UintStr(c);
      switch(v.type){
        case JBinaryDataDef::DatText:     v.vtext=   ((const string *)       data)[c];   break;
        case JBinaryDataDef::DatBool:     v.vint=   (((const bool *)data)[c]? 1: 0);     break;
        case JBinaryDataDef::DatChar:     v.vchar=   ((const char *)         data)[c];   break;
        case JBinaryDataDef::DatUchar:    v.vuchar=  ((const unsigned char *)data)[c];   break;
        case JBinaryDataDef::DatShort:    v.vshort=  ((const short *)        data)[c];   break;
        case JBinaryDataDef::DatUshort:   v.vushort= ((const word *)         data)[c];   break;
        case JBinaryDataDef::DatInt:      v.vint=    ((const int *)          data)[c];   break;   
        case JBinaryDataDef::DatUint:     v.vuint=   ((const unsigned *)     data)[c];   break;   
        case JBinaryDataDef::DatLlong:    v.vllong=  ((const llong *)        data)[c];   break;   
        case JBinaryDataDef::DatUllong:   v.vullong= ((const ullong *)       data)[c];   break;   
        case JBinaryDataDef::DatFloat:    v.vfloat=  ((const float *)        data)[c];   break;   
        case JBinaryDataDef::DatDouble:   v.vdouble= ((const double *)       data)[c];   break;   
        case JBinaryDataDef::DatInt3:     v.vint3=   ((const tint3 *)        data)[c];   break;   
        case JBinaryDataDef::DatUint3:    v.vuint3=  ((const tuint3 *)       data)[c];   break;   
        case JBinaryDataDef::DatFloat3:   v.vfloat3= ((const tfloat3 *)      data)[c];   break;   
        case JBinaryDataDef::DatDouble3:  v.vdouble3=((const tdouble3 *)     data)[c];   break;   
        default: Run_Exceptioon("Type of value invalid.");
      }
      (*pf) << tabs << "\t" << ValueToXml(v) << endl;
    }
    (*pf) << tabs << string("</array_")+tx+">" << endl;
  }
}

//==============================================================================
/// Graba contenido en fichero XML.
/// Saves contents to XML file.
//==============================================================================
void JBinaryData::WriteFileXml(const std::string &tabs,std::ofstream* pf,bool svarrays)const{
  (*pf) << tabs << "<item name=\"" << GetName() << "\" hide=\"" << (GetHide()? '1': '0') << "\" hidevalues=\"" << (GetHideValues()? '1': '0') << "\">" << endl;
  for(unsigned c=0;c<Values.size();c++)(*pf) << tabs << "\t" << ValueToXml(Values[c]) << endl;
  for(unsigned c=0;c<Arrays.size();c++)WriteFileXmlArray(tabs+"\t",pf,svarrays,Arrays[c]);
  for(unsigned c=0;c<Items.size();c++)Items[c]->WriteFileXml(tabs+"\t",pf,svarrays);
  (*pf) << tabs << "</item>" << endl;
}

//==============================================================================
/// Cambia nombre de objeto comprobando que no tenga hermanos con el mismo nombre.
/// Changes the name of the object checking that there is no other object with 
/// the same name 
//==============================================================================
void JBinaryData::SetName(const std::string &name){
  if(Parent){
    if(Parent->ExistsValue(name))Run_Exceptioon("There is already a value with the name given.");
    if(Parent->GetArray(name)!=NULL)Run_Exceptioon("There is already an array with the name given.");
    if(Parent->GetItem(name)!=NULL)Run_Exceptioon("There is already an item with the name given.");
  }
  Name=name;
}

//==============================================================================
/// Cambia oculatacion de values.
/// Change the SetHide of the values
//==============================================================================
void JBinaryData::SetHideValues(bool hide,bool down){
  HideValues=hide;
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetHideValues(hide,true);
}

//==============================================================================
/// Cambia oculatacion de arrays.
/// Change the SetHide of the arrays
//==============================================================================
void JBinaryData::SetHideArrays(bool hide,bool down){
  for(unsigned c=0;c<Arrays.size();c++){
    Arrays[c]->SetHide(hide);
    if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetHideArrays(hide,true);
  }
}

//==============================================================================
/// Cambia oculatacion de items.
/// Change the SetHide of the items
//==============================================================================
void JBinaryData::SetHideItems(bool hide,bool down){
  for(unsigned c=0;c<Items.size();c++){
    Items[c]->SetHide(hide);
    if(down)Items[c]->SetHideItems(hide,true);
  }
}

//==============================================================================
/// Cambia formato de texto para valores float.
/// Change text format for floats.
//==============================================================================
void JBinaryData::SetFmtFloat(const std::string &fmt,bool down){
  FmtFloat=fmt;
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetFmtFloat(fmt,true);
}

//==============================================================================
/// Cambia formato de texto para valores double.
/// Change text format for doubles.
//==============================================================================
void JBinaryData::SetFmtDouble(const std::string &fmt,bool down){
  FmtDouble=fmt;
  if(down)for(unsigned c=0;c<Items.size();c++)Items[c]->SetFmtDouble(fmt,true);
}

//==============================================================================
/// Devuelve el tamanho necesario para almacenar todos los datos del item y descendientes.
/// Con all activado se incluyen tambien los elementos ocultos.
/// Returns the size needed to store the data from the item and "descendants"
/// With bool "all" true the hidden elements are also included. 
//==============================================================================
unsigned JBinaryData::GetSizeDataConst(bool all)const{
  unsigned count=0;
  InItem(count,0,NULL,all);
  return(count);
}

//==============================================================================
/// Almacena datos de item en ptr y devuelve los bytes almacenados.
/// Con all activado se incluyen tambien los elementos ocultos.
/// Stores data of item in ptr and and returns the bytes stored.
// With bool "all" true the hidden elements are also included.
//==============================================================================
unsigned JBinaryData::SaveDataConst(unsigned size,byte* ptr,bool all)const{
  if(!ptr)Run_Exceptioon("The pointer is invalid.");
  unsigned count=0;
  InItem(count,size,ptr,all);
  return(count);
}

//==============================================================================
/// Devuelve el tamanho necesario para almacenar todos los del item y 
/// descendientes. Actualiza la cache de values  para mejorar el rendimiento en 
/// operaciones posteriores.
/// Con all activado se incluyen tambien los elementos ocultos.
/// Returns the size needed to store all the item and descendants.
/// Updates the cache of values to improve performance for later operations.
/// With bool "all" true the hidden items are also included.
//==============================================================================
unsigned JBinaryData::GetSizeData(bool all){
  ValuesCachePrepare(true);
  unsigned count=0;
  InItem(count,0,NULL,all);
  return(count);
}

//==============================================================================
/// Almacena datos de item en ptr y devuelve los bytes almacenados. Actualiza la
/// cache de values para mejorar el rendimiento en operaciones posteriores.
/// Con all activado se incluyen tambien los elementos ocultos.
/// Stores data of item in ptr and returns the bytes stored. 
/// Updates cache values to improve performance in subsequent operations.
/// Activated Also included all the hidden items.
/// With bool "all" true  the hidden items are also included.
//==============================================================================
unsigned JBinaryData::SaveData(unsigned size,byte* ptr,bool all){
  if(!ptr)Run_Exceptioon("The pointer is invalid.");
  ValuesCachePrepare(true);
  unsigned count=0;
  InItem(count,size,ptr,all);
  return(count);
}

//==============================================================================
/// Carga datos desde ptr.
/// Loads data from ptr.
//==============================================================================
void JBinaryData::LoadData(unsigned size,const byte* ptr){
  if(!ptr)Run_Exceptioon("The pointer is invalid.");
  Clear(); //-Limpia contenido de objeto. Clean object content.
  unsigned count=0;
  OutItem(count,size,ptr,false);
}


//==============================================================================
/// Graba datos en fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido si no hay arrays grandes.
/// Con all activado se incluyen tambien los elementos ocultos.
/// Save data file.
/// With a buffer memory used for all data. It consumes more memory but
/// it is faster if there are large arrays.
/// With bool "all" true  the hidden items are also included.
//==============================================================================
void JBinaryData::SaveFileData(std::fstream *pf,bool head,const std::string &filecode,bool memory,bool all)const{
  if(head){//-Graba cabecera basica.
    StHeadFmtBin head=MakeFileHead(filecode); 
    pf->write((char*)&head,sizeof(StHeadFmtBin));
  }
  //-Graba datos. Save data.
  if(memory){//-Graba datos desde memoria. Write data from memory.
    const unsigned sbuf=GetSizeDataConst(all);
    //printf("SaveFile> sbuf:%u\n",sbuf);
    byte *buf=new byte[sbuf];
    unsigned sbuf2=SaveDataConst(sbuf,buf,all);
    //printf("SaveFile> sbuf2:%u\n",sbuf2);
    pf->write((char*)buf,sbuf);
    delete[] buf;
  }
  else{//-Graba datos directamente. Save data directly.
    const unsigned sbuf=1024;
    byte buf[sbuf];
    WriteItem(pf,sbuf,buf,all);
  }
}

//==============================================================================
/// Graba datos en fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido si no hay arrays grandes.
/// Con all activado se incluyen tambien los elementos ocultos.
/// Save data file.
/// With a buffer memory used for all data. It consumes more memory but
/// It is faster if there are large arrays.
/// With bool "all" true  the hidden items are also included.
//==============================================================================
void JBinaryData::SaveFile(const std::string &file,bool memory,bool all){
  ValuesCachePrepare(true);
  fstream pf;
  pf.open(file.c_str(),ios::binary|ios::out);
  if(pf){
    SaveFileData(&pf,true,Name,memory,all);
    if(pf.fail())Run_ExceptioonFile("File writing failure.",file);
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",file);
}

//==============================================================================
/// Carga datos de un fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido cuando no hay arrays grandes.
/// Load data from a file.
/// With a buffer memory used for all data. It consumes more memory but
/// it is faster if there are large arrays.
//==============================================================================
void JBinaryData::LoadFile(const std::string &file,const std::string &filecode,bool memory){
  Clear(); //-Limpia contenido de objeto. Clean object content.
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    const unsigned fsize=CheckFileHead(file,&pf,filecode);
    //-Carga datos.
    if(memory){//-Carga datos desde memoria. Write data from memory.
      const unsigned sbuf=fsize-sizeof(StHeadFmtBin);
      //printf("LoadFile> sbuf:%u\n",sbuf);
      byte *buf=new byte[sbuf];
      pf.read((char*)buf,sbuf);
      LoadData(sbuf,buf);
      delete[] buf;
    }
    else{//-Carga datos directamente. Save data directly.
      const unsigned sbuf=1024;
      byte buf[sbuf];
      ReadItem(&pf,sbuf,buf,false,true);
    }
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",file);
}

//==============================================================================
/// Graba Item en fichero a continuacion de los items que ya exitan.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido si no hay arrays grandes.
/// Con all activado se incluyen tambien los elementos ocultos.

/// Save Item in file below from the items that already existed.
/// With a buffer memory used for all data. It consumes more memory but
/// it is faster if there are large arrays.
/// With bool "all" true the hidden items are also included.
//==============================================================================
void JBinaryData::SaveFileListApp(const std::string &file,const std::string &filecode,bool memory,bool all){
  ValuesCachePrepare(true);
  fstream pf;
  if(fun::FileExists(file))pf.open(file.c_str(),ios::binary|ios::out|ios::in|ios::app);
  else pf.open(file.c_str(),ios::binary|ios::out);
  if(pf){
    const unsigned fsize=CheckFileListHead(file,&pf,filecode);
    pf.seekp(0,pf.end);
    //-Graba datos de parent. Save parent data.
    if(!fsize)Parent->SaveFileData(&pf,true,filecode,memory,all);
    //-Graba datos de item. Save item data.
    SaveFileData(&pf,false,filecode,memory,all);
    if(pf.fail())Run_ExceptioonFile("File writing failure.",file);
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",file);
}

//==============================================================================
/// Carga datos de un fichero.
/// Con memory utiliza un buffer para todos los datos. Consume mas memoria pero
/// es mas rapido cuando no hay arrays grandes.
/// Load data from a file.
/// With a buffer memory used for all data. It consumes more memory but
/// it is faster if there are large arrays.
//==============================================================================
void JBinaryData::LoadFileListApp(const std::string &file,const std::string &filecode,bool memory){
  Clear(); //-Limpia contenido de objeto. Clean object content.
  ifstream pf;
  pf.open(file.c_str(),ios::binary|ios::in);
  if(pf){
    SetName(filecode);
    const unsigned fsize=CheckFileHead(file,&pf,filecode);
    unsigned pfile=(unsigned)pf.tellg();
    while(pfile<fsize){
      //-Carga datos.
      if(memory){//-Carga datos desde memoria. Loads data from memory.
        const unsigned sbuf=fsize-sizeof(StHeadFmtBin);
        //printf("LoadFile> sbuf:%u\n",sbuf);
        byte *buf=new byte[sbuf];
        pf.read((char*)buf,sbuf);
        unsigned cbuf=0;
        while(cbuf<sbuf){
          OutItem(cbuf,sbuf,buf,true);
          //-Renombra ultimo item leido. Renames last read item.
          unsigned lastitem=GetItemsCount()-1;
          JBinaryData* ite=GetItem(lastitem);
          ite->SetName(fun::PrintStr("LS%04u_",lastitem)+ite->GetName());
        }
        delete[] buf;
      }
      else{//-Carga datos directamente. Load data directly.
        const unsigned sbuf=1024;
        byte buf[sbuf];
        ReadItem(&pf,sbuf,buf,true,true);
        //-Renombra ultimo item leido. Renames last read item.
        unsigned lastitem=GetItemsCount()-1;
        JBinaryData* ite=GetItem(lastitem);
        ite->SetName(fun::PrintStr("LS%04u_",lastitem)+ite->GetName());
      }
      pfile=(unsigned)pf.tellg();
    }
    pf.close();
  }
  else Run_ExceptioonFile("Cannot open the file.",file);
}

//==============================================================================
/// Abre fichero y carga estructura de datos pero sin cargar el contenido de los
/// arrays.
/// Open file and load data structure but without loading the contents of the
/// arrays.
//==============================================================================
void JBinaryData::OpenFileStructure(const std::string &file,const std::string &filecode){
  if(Parent)Run_Exceptioon("Item is not root.");
  Clear(); //-Limpia contenido de objeto. Clean object content.
  FileStructure=new ifstream;
  FileStructure->open(file.c_str(),ios::binary|ios::in);
  if(*FileStructure){
    const unsigned fsize=CheckFileHead(file,FileStructure,filecode);
    const unsigned sbuf=1024;
    byte buf[sbuf];
    ReadItem(FileStructure,sbuf,buf,false,false);
  }
  else{
    CloseFileStructure();
    Run_ExceptioonFile("Cannot open the file.",file);
  }
}

//==============================================================================
/// Cierra fichero abierto con OpenFileStructure().
/// Close open file using OpenFileStructure ().
//==============================================================================
void JBinaryData::CloseFileStructure(){
  if(FileStructure&&FileStructure->is_open())FileStructure->close();
  delete FileStructure; FileStructure=NULL;
}

//==============================================================================
/// Devuelve puntero al fichero abierto con OpenFileStructure().
/// Returns pointer to open file with OpenFileStructure ().
//==============================================================================
std::ifstream* JBinaryData::GetFileStructure()const{
  if(Parent)Run_Exceptioon("Item is not root.");
  return(FileStructure);
}

//==============================================================================
/// Graba contenido en fichero XML.
/// Record XML file content.
//==============================================================================
void JBinaryData::SaveFileXml(std::string file,bool svarrays,const std::string &head)const{
  file=fun::GetWithoutExtension(file)+".xml";
  ofstream pf;
  pf.open(file.c_str());
  if(pf){
    pf << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
    pf << "<data" << head <<" date=\"" << fun::GetDateTime() << "\">" << endl;
    WriteFileXml("\t",&pf,svarrays);
    pf << "</data>" << endl;
    if(pf.fail())Run_ExceptioonFile("Failed writing to file.",file);
    pf.close();
  }
  else Run_ExceptioonFile("File could not be opened.",file);
}

//==============================================================================
/// Devuelve item principal.
/// Returns main item.
//==============================================================================
JBinaryData* JBinaryData::GetItemRoot(){
  return(Parent? Parent->GetItemRoot(): this);
}

//==============================================================================
/// Devuelve el numero de items no marcados como ocultos.
/// Returns the number of items not marked as hidden.
//==============================================================================
unsigned JBinaryData::GetVisibleItemsCount()const{
  unsigned num=0;
  for(unsigned c=0;c<Items.size();c++)if(!Items[c]->GetHide())num++;
  return(num);
}

//==============================================================================
/// Devuelve indice del item con el nombre indicado o -1 si no existe.
//==============================================================================
int JBinaryData::GetItemIndex(const std::string &name){
  int idx=-1;
  for(unsigned c=0;c<Items.size()&&idx<0;c++)if(Items[c]->Name==name)idx=c;
  return(idx);
}

//==============================================================================
/// Devuelve item con el nombre indicado o NULL si no existe.
/// Returns item with the specified name or NULL if not present.
//==============================================================================
JBinaryData* JBinaryData::GetItem(const std::string &name){
  JBinaryData* ret=NULL;
  for(unsigned c=0;c<Items.size()&&!ret;c++)if(Items[c]->Name==name)ret=Items[c];
  return(ret);
}

//==============================================================================
/// Devuelve item segun el indice indicado o NULL si no existe.
/// Returns index according to the indicated item or NULL if no index exists.
//==============================================================================
JBinaryData* JBinaryData::GetItem(unsigned index){
  return(index>=GetItemsCount()? NULL: Items[index]);
}

//==============================================================================
/// Crea y devuelve item con el nombre. Genera excepcion si ya existe.
/// Creates and returns item with the name. It generates exception if it already exists.
//==============================================================================
JBinaryData* JBinaryData::CreateItem(const std::string &name){
  if(GetItem(name)!=NULL)Run_Exceptioon("There is already an item with the name given.");
  JBinaryData *item=new JBinaryData(name);
  item->Parent=this;
  Items.push_back(item);
  return(item);
}

//==============================================================================
/// Elimina el item indicado.
/// Deletes the item indicated.
//==============================================================================
void JBinaryData::RemoveItem(const std::string &name){
  int idx=GetItemIndex(name);
  if(idx>=0){
    JBinaryData* item=Items[idx];
    Items.erase(Items.begin()+idx);
    delete item;
  }
}

//==============================================================================
/// Elimina todos los items almacenados.
/// Remove all stored items.
//==============================================================================
void JBinaryData::RemoveItems(){
  for(unsigned c=0;c<Items.size();c++)delete Items[c];
  Items.clear();
}

//==============================================================================
/// Devuelve el numero de arrays no marcados como ocultos.
/// Returns the number of arrays not marked as hidden.
//==============================================================================
unsigned JBinaryData::GetVisibleArraysCount()const{
  unsigned num=0;
  for(unsigned c=0;c<Arrays.size();c++)if(!Arrays[c]->GetHide())num++;
  return(num);
}

//==============================================================================
/// Devuelve posicion de la variable solicitada, -1 en caso de no existir.
/// Returns the requested position variable, -1 if does not exist.
//==============================================================================
int JBinaryData::GetArrayIndex(const std::string &name)const{
  int idx=-1; 
  for(unsigned c=0;c<Arrays.size()&&idx<0;c++)if(Arrays[c]->GetName()==name)idx=c;
  return(idx);
}

//==============================================================================
/// Devuelve array con el nombre indicado o NULL si no existe.
/// Returns array with the specified name or NULL if does not exist.
//==============================================================================
JBinaryDataArray* JBinaryData::GetArray(const std::string &name){
  JBinaryDataArray* ret=NULL;
  for(unsigned c=0;c<Arrays.size()&&!ret;c++)if(Arrays[c]->GetName()==name)ret=Arrays[c];
  return(ret);
}

//==============================================================================
/// Devuelve array segun el indice indicado o NULL si no existe.
/// Returns array according to the specified index or NULL if not present.
//==============================================================================
JBinaryDataArray* JBinaryData::GetArray(unsigned index){
  return(index>=GetArraysCount()? NULL: Arrays[index]);
}

//==============================================================================
/// Crea y devuelve array con el nombre. Genera excepcion si ya existe.
/// Creates and returns array with the name. It generates exception if it already exists.
//==============================================================================
JBinaryDataArray* JBinaryData::CreateArray(const std::string &name,JBinaryDataDef::TpData type){
  if(GetItem(name)!=NULL)Run_Exceptioon("There is already an array with the name given.");
  JBinaryDataArray *ar=new JBinaryDataArray(this,name,type);
  Arrays.push_back(ar);
  return(ar);
}

//==============================================================================
/// Crea y devuelve array con datos.
/// Creates and returns array data.
//==============================================================================
JBinaryDataArray* JBinaryData::CreateArray(const std::string &name,JBinaryDataDef::TpData type,unsigned count,const void *data,bool externalpointer){
  JBinaryDataArray *ar=CreateArray(name,type);
  ar->SetData(count,data,externalpointer);
  return(ar);
}

//==============================================================================
/// Elimina el array indicado.
/// Removes the specified array.
//==============================================================================
void JBinaryData::RemoveArray(const std::string &name){
  int idx=GetArrayIndex(name);
  if(idx>=0){
    JBinaryDataArray* ar=Arrays[idx];
    Arrays.erase(Arrays.begin()+idx);
    delete ar;
  }
}

//==============================================================================
/// Elimina todos los arrays almacenados.
/// Remove all storage arrays.
//==============================================================================
void JBinaryData::RemoveArrays(){
  for(unsigned c=0;c<Arrays.size();c++)delete Arrays[c];
  Arrays.clear();
}

//==============================================================================
/// Comprueba copia datos de array despues de comprobar tipo y tamanho.
/// Check copy data from array after checking type and size.
//==============================================================================
JBinaryDataArray* JBinaryData::CheckCopyArrayData(const std::string &name,unsigned size,JBinaryDataDef::TpData type){
  JBinaryDataArray* ar=GetArray(name);
  if(!ar)Run_Exceptioon(fun::PrintStr("Array \'%s\' is not available.",name.c_str()));
  if(ar->GetType()!=type)Run_Exceptioon(fun::PrintStr("Type of array \'%s\' is not %s.",name.c_str(),JBinaryDataDef::TypeToStr(type).c_str()));
  if(ar->GetCount()!=size)Run_Exceptioon(fun::PrintStr("Size of array \'%s\' does not match.",name.c_str()));
  return(ar);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,char           *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatChar)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,unsigned char  *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatUchar)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,short          *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatShort)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,unsigned short *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatUshort)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,int            *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatInt)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,unsigned       *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatUint)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,llong          *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatLlong)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,ullong         *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatUllong)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,float          *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatFloat)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,double         *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatDouble)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,tint3          *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatInt3)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,tuint3         *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatUint3)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,tfloat3        *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatFloat3)->GetDataCopy(size,ptr);
}
//==============================================================================
/// Copia datos de array despues de comprobar tipo y tamanho.
/// Copy data from array after checking type and size.
//==============================================================================
void JBinaryData::CopyArrayData(const std::string &name,unsigned size,tdouble3       *ptr){
  CheckCopyArrayData(name,size,JBinaryDataDef::DatDouble3)->GetDataCopy(size,ptr);
}

//==============================================================================
/// Devuelve posicion de la variable solicitada, -1 en caso de no existir.
/// Returns the requested position variable, -1 if does not exist.
//==============================================================================
int JBinaryData::GetValueIndex(const std::string &name)const{
  int pos=-1; 
  for(unsigned c=0;c<Values.size()&&pos<0;c++)if(Values[c].name==name)pos=c;
  return(pos);
}

//==============================================================================
/// Devuelve el nombre del value solicitado (vacio si no existe).
/// Returns the name of the requested value (empty if does not exist).
//==============================================================================
std::string JBinaryData::NameOfValue(unsigned index)const{
  return(index>=GetValuesCount()? "": Values[index].name);
}

//==============================================================================
/// Devuelve el tipo del value solicitado (DatNull si no existe).
/// Returns the type of value requested (DatNull if any).
//==============================================================================
JBinaryDataDef::TpData JBinaryData::TypeOfValue(const std::string &name)const{
  int idx=GetValueIndex(name);
  return(idx<0? JBinaryDataDef::DatNull: Values[idx].type);
}

//==============================================================================
/// Devuelve el tipo del value solicitado (DatNull si no existe).
/// Returns the type of value requested (DatNull if any).
//==============================================================================
JBinaryDataDef::TpData JBinaryData::TypeOfValue(unsigned index)const{
  return(index>=GetValuesCount()? JBinaryDataDef::DatNull: Values[index].type);
}

//==============================================================================
/// Indica si existe el value solicitado.
/// Indicates whether the requested value exists.
//==============================================================================
bool JBinaryData::ExistsValue(const std::string &name)const{
  return(GetValueIndex(name)>=0);
}

//==============================================================================
/// Indica si existe el value solicitado del tipo indicado.
/// It indicates the existance of the requested value of the type.
//==============================================================================
bool JBinaryData::ExistsValue(const std::string &name,JBinaryDataDef::TpData type)const{
  int idx=GetValueIndex(name);
  return(idx>=0&&Values[idx].type==type);
}

//==============================================================================
/// Elimina el value indicado.
/// Deletes the value indicated.
//==============================================================================
void JBinaryData::RemoveValue(const std::string &name){
  int idx=GetValueIndex(name);
  if(idx>=0)Values.erase(Values.begin()+idx);
  ValuesModif=true;
}

//==============================================================================
/// Elimina todos los values almacenados.
/// Removes all stored values.
//==============================================================================
void JBinaryData::RemoveValues(){
  Values.clear();
  ValuesCacheReset();
}

//==============================================================================
/// Devuelve el valor solicitado de tipo texto.
/// Returns the requested type value text.
//==============================================================================
std::string JBinaryData::GetvText(const std::string &name,bool optional,std::string valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatText);
  return(pos<0? valdef: Values[pos].vtext);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo bool.
/// Returns the requested value of bool type.
//==============================================================================
bool JBinaryData::GetvBool(const std::string &name,bool optional,bool valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatBool);
  return(pos<0? valdef: Values[pos].vint!=0);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo char.
/// Returns the requested value of char type .
//==============================================================================
char JBinaryData::GetvChar(const std::string &name,bool optional,char valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatChar);
  return(pos<0? valdef: Values[pos].vchar);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned char.
/// Returns the requested value of unsigned char type.
//==============================================================================
unsigned char JBinaryData::GetvUchar(const std::string &name,bool optional,unsigned char valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUchar);
  return(pos<0? valdef: Values[pos].vuchar);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo short.
/// Returns the requested value of type short.
//==============================================================================
short JBinaryData::GetvShort(const std::string &name,bool optional,short valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatShort);
  return(pos<0? valdef: Values[pos].vshort);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned short.
/// Returns the requested value of type unsigned short.
//==============================================================================
unsigned short JBinaryData::GetvUshort(const std::string &name,bool optional,unsigned short valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUshort);
  return(pos<0? valdef: Values[pos].vushort);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo int.
/// Returns the requested value of type int.
//==============================================================================
int JBinaryData::GetvInt(const std::string &name,bool optional,int valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatInt);
  return(pos<0? valdef: Values[pos].vint);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned.
/// Returns the requested value of unsigned type.
//==============================================================================
unsigned JBinaryData::GetvUint(const std::string &name,bool optional,unsigned valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUint);
  return(pos<0? valdef: Values[pos].vuint);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo long long.
/// Returns the requested value long long.
//==============================================================================
llong JBinaryData::GetvLlong(const std::string &name,bool optional,llong valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatLlong);
  return(pos<0? valdef: Values[pos].vllong);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo unsigned long long.
/// Returns the requested value of type unsigned long long.
//==============================================================================
ullong JBinaryData::GetvUllong(const std::string &name,bool optional,ullong valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUllong);
  return(pos<0? valdef: Values[pos].vullong);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo float.
/// Returns the requested value of type float.
//==============================================================================
float JBinaryData::GetvFloat(const std::string &name,bool optional,float valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatFloat);
  return(pos<0? valdef: Values[pos].vfloat);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo double.
/// Returns the requested value of type double.
//==============================================================================
double JBinaryData::GetvDouble(const std::string &name,bool optional,double valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatDouble);
  return(pos<0? valdef: Values[pos].vdouble);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo tint3.
/// Returns the value of tint3 type requested.
//==============================================================================
tint3 JBinaryData::GetvInt3(const std::string &name,bool optional,tint3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatInt3);
  return(pos<0? valdef: Values[pos].vint3);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo tuint3.
/// Returns the value of tint3 type requested.
//==============================================================================
tuint3 JBinaryData::GetvUint3(const std::string &name,bool optional,tuint3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatUint3);
  return(pos<0? valdef: Values[pos].vuint3);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo tfloat3.
/// Returns the requested value of type tfloat3.
//==============================================================================
tfloat3 JBinaryData::GetvFloat3(const std::string &name,bool optional,tfloat3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatFloat3);
  return(pos<0? valdef: Values[pos].vfloat3);
}

//==============================================================================
/// Devuelve el valor solicitado de tipo tdouble3.
/// Returns the requested value of type tdouble3.
//==============================================================================
tdouble3 JBinaryData::GetvDouble3(const std::string &name,bool optional,tdouble3 valdef)const{
  int pos=CheckGetValue(name,optional,JBinaryDataDef::DatDouble3);
  return(pos<0? valdef: Values[pos].vdouble3);
}

//==============================================================================
/// Crea o modifica un valor de tipo texto.
/// Creates or modifies a text value.
//==============================================================================
void JBinaryData::SetvText(const std::string &name,const std::string &v){
  Values[CheckSetValue(name,JBinaryDataDef::DatText)].vtext=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo bool.
/// Creates or modifies a value of type bool.
//==============================================================================
void JBinaryData::SetvBool(const std::string &name,bool v){
  Values[CheckSetValue(name,JBinaryDataDef::DatBool)].vint=(v? 1: 0);
}

//==============================================================================
/// Crea o modifica un valor de tipo char.
/// Creates or modifies a value of type char.
//==============================================================================
void JBinaryData::SetvChar(const std::string &name,char v){
  Values[CheckSetValue(name,JBinaryDataDef::DatChar)].vchar=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo unsigned char.
/// Creates or modifies a value of type unsigned char.
//==============================================================================
void JBinaryData::SetvUchar(const std::string &name,unsigned char v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUchar)].vuchar=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo short.
/// Creates or modifies a value of type short.
//==============================================================================
void JBinaryData::SetvShort(const std::string &name,short v){
  Values[CheckSetValue(name,JBinaryDataDef::DatShort)].vshort=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo unsigned short.
/// Creates or modifies a value of type unsigned short
//==============================================================================
void JBinaryData::SetvUshort(const std::string &name,unsigned short v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUshort)].vushort=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo int.
/// Creates or modifies a value of type int.
//==============================================================================
void JBinaryData::SetvInt(const std::string &name,int v){
  Values[CheckSetValue(name,JBinaryDataDef::DatInt)].vint=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo unsigned.
/// Creates or modifies a value of type unsigned.
//==============================================================================
void JBinaryData::SetvUint(const std::string &name,unsigned v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUint)].vuint=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo long long.
/// Creates or modifies a value of type long long.
//==============================================================================
void JBinaryData::SetvLlong(const std::string &name,llong v){
  Values[CheckSetValue(name,JBinaryDataDef::DatLlong)].vllong=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo unsigned long long.
/// Creates or modifies a value of type unsigned long long.
//==============================================================================
void JBinaryData::SetvUllong(const std::string &name,ullong v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUllong)].vullong=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo float.
/// Creates or modifies a value of type float.
//==============================================================================
void JBinaryData::SetvFloat(const std::string &name,float v){
  Values[CheckSetValue(name,JBinaryDataDef::DatFloat)].vfloat=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo double.
/// Creates or modifies a value of type double.
//==============================================================================
void JBinaryData::SetvDouble(const std::string &name,double v){
  Values[CheckSetValue(name,JBinaryDataDef::DatDouble)].vdouble=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo tint3.
/// Creates or modifies a value of type tint3.
//==============================================================================
void JBinaryData::SetvInt3(const std::string &name,tint3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatInt3)].vint3=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo tuint3.
/// Creates or modifies a value of type tuint3.
//==============================================================================
void JBinaryData::SetvUint3(const std::string &name,tuint3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatUint3)].vuint3=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo tfloat3.
/// Creates or modifies a value of type tfloat3.
//==============================================================================
void JBinaryData::SetvFloat3(const std::string &name,tfloat3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatFloat3)].vfloat3=v;
}

//==============================================================================
/// Crea o modifica un valor de tipo tdouble3.
/// Creates or modifies a value of type tdouble3.
//==============================================================================
void JBinaryData::SetvDouble3(const std::string &name,tdouble3 v){
  Values[CheckSetValue(name,JBinaryDataDef::DatDouble3)].vdouble3=v;
}


