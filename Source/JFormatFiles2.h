//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2017 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

 EPHYSLAB Environmental Physics Laboratory, Universidade de Vigo, Ourense, Spain.
 School of Mechanical, Aerospace and Civil Engineering, University of Manchester, Manchester, U.K.

 This file is part of DualSPHysics. 

 DualSPHysics is free software: you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License 
 as published by the Free Software Foundation; either version 2.1 of the License, or (at your option) any later version.
 
 DualSPHysics is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details. 

 You should have received a copy of the GNU Lesser General Public License along with DualSPHysics. If not, see <http://www.gnu.org/licenses/>. 
*/

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Clase creada para sustituir a JFormatFiles. (15-12-2010)
//:# - Nuevas funciones para generar ficheros con datos de puntos en formato
//:#   CSV y ASCII:  SaveCsvPointsVar(), SaveCsvPointsVar3(), SaveAscPointsVar(), 
//:#   SaveAscPointsVar3() (22-12-2010)
//:# - Nueva variable ACE en formatos VTK, CSV y ASCII. (16-06-2011)
//:# - Nueva variable Vorticity en formatos VTK, CSV y ASCII. (13-07-2011)
//:# - Se separa la variable ACE en parte positiva y negativa en formatos CSV y 
//:#   ASCII (ParticlesToCsv2() y ParticlesToAscii()). (28-09-2011)
//:# - En los VTK de puntos, cada punto se guarda como una celda. Esto permite
//:#   el uso de algunas operaciones como Threshold en Paraview. Sin embargo,
//:#   aumenta el tamaño en 3 bytes por partícula (entorno al 11%). (03-12-2011)
//:# - Nueva funcion SaveVtkCells() para generar ficheros VTK con una malla de 
//:#   celdas segun los parametros especificados. (13-12-2011)
//:# - Nueva funcion PointsToVtk() para generar ficheros VTK de puntos a partir 
//:#   de datos codificados en un JBuffer. (30-12-2011)
//:# - Nueva funcion CellsToVtk() para generar ficheros VTK de celdas sin datos
//:#   a partir de datos codificados en un JBuffer. (09-01-2012)
//:# - Nueva funcion SaveVtkDomain() para generar ficheros VTK con una serie
//:#   de dominios. (09-01-2012)
//:# - En PointsToVtk() si el numero de puntos es cero, se crea uno para evitar
//:#   un error al generar el fichero VTK. (18-01-2012)
//:# - Traduccion de comentarios al ingles. (10-02-2012)
//:# - Nuevas funciones (SaveVtk,SaveCsv) para crear ficheros VTK mas 
//:#   generales. (17-05-2012)
//:# - Nueva funcion SaveVtkBasic como ejemplo de uso de SaveVtk. (06-06-2012)
//:# - Funcion SaveVtkBox() para generar ficheros VTK con una serie de 
//:#   cajas. (06-06-2012)
//:# - Los metodos PointsToVtk() y CellsToVtk() se eleminaron de esta clase y
//:#   pasaron a la clase JBufferToVtk. (06-06-2012)
//:# - Error corregido en SaveVtkBasic() cuando se usaban arrays nulos. (06-10-2012)
//:# - Se elimino un posible problema de precision en el calculo de los vertices
//:#   de una caja. (15-05-2013)
//:# - Permite añadir una cabecera de fichero usando SaveCsv(...,head). (20-11-2013)
//:# - Nuevos metodos SaveCsvPointsVar(), SaveCsvPointsVar3(), SaveAscPointsVar()
//:#   y SaveAscPointsVar3() para grabar la posicion como tdouble3. (30-12-2013)
//:# - Nuevo metodo SaveAscii() mediente StScalarData* fields. (16-01-2014)
//:# - Nuevos metodos XXXStats() que permiten calcular y grabar ficheros CSV
//:#   con valores de minimo, maximo y media. (16-01-2014)
//:# - Error corregido en funcion DefineStatsField(). (18-02-2015)
//:# - Funciones SaveCsv() y SaveAsc() para POS en doble precision. (24-03-2015)
//:# - En funciones SaveCsv() se indica el nombre del campo. (03-03-2016)
//:# - Nueva clase JFormatFiles2Data para simplificar la generacion de ficheros
//:#   VTK y CSV. (07-04-2017)
//:# - New parameter csvsep to configure separator in CSV files. (23-10-2017)
//:#############################################################################

/// \file JFormatFiles2.h \brief Declares the class \ref JFormatFiles2.

#ifndef _JFormatFiles2_
#define _JFormatFiles2_

#include "TypesDef.h"
#include "JObject.h"
#include <string>
#include <cstring>
#include <string>

//##############################################################################
//# JFormatFiles2
//##############################################################################
/// \brief Provides functions to store particle data in formats VTK, CSV, ASCII.

class JFormatFiles2
{
public:

  typedef enum{ UChar8,Char8,UShort16,Short16,UInt32,Int32,Float32,Double64,ULlong64,Llong64,TpDataNull }TpData;

  /// Structure with the information of an array of particle data to be stored in CSV or VTK format.
  typedef struct {
    std::string name;
    TpData type;
    unsigned comp;
    const void *pointer;
  }StScalarData;

  /// Strucutre with the information of an array to calculate and save statistic information.
  typedef struct {
    //-Data of arrays.
    std::string name;
    TpData type;
    unsigned comp;
    const void *pointer;
    //-Selection of results.
    bool selmin,selmax,selmean;
    bool seltotal,selcomponent;
    //-Results.
    ullong num;
    double min,max,mean;
    double min1,max1,mean1;
    double min2,max2,mean2;
    double min3,max3,mean3;
  }StStatistics;

  //==============================================================================
  /// Throws a simple exception.
  //==============================================================================
  static void RunException(std::string method,std::string msg);
  
  //==============================================================================
  /// Throws an exception related to a file.
  //==============================================================================  
  static void RunException(std::string method,std::string msg,std::string file);

  //==============================================================================
  /// Returns the definition of fields.
  //==============================================================================
  static StScalarData DefineField(const std::string &name,TpData type,unsigned comp,const void *pointer=NULL){
    StScalarData f; f.name=name; f.type=type; f.comp=comp; f.pointer=pointer;
    return(f);
  }

  //==============================================================================
  /// Checks the definition of fields.
  //==============================================================================
  static void CheckFields(unsigned nfields,const StScalarData* fields);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void SaveVtk(std::string fname,unsigned np
    ,const tfloat3* pos,unsigned nfields,const StScalarData* fields);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void SaveVtk(std::string fname,unsigned np
    ,const tdouble3* pos,unsigned nfields,const StScalarData* fields);
  
  //==============================================================================
  /// Stores data in CSV format.
  //============================================================================== 
  static void SaveCsv(std::string fname,bool csvsepcoma,unsigned np
    ,const tfloat3* pos,unsigned nfields,const StScalarData* fields,std::string head="");
  
  //==============================================================================
  /// Stores data in CSV format.
  //============================================================================== 
  static void SaveCsv(std::string fname,bool csvsepcoma,unsigned np
    ,const tdouble3* pos,unsigned nfields,const StScalarData* fields,std::string head="");
  
  //==============================================================================
  /// Stores data in ASCII format.
  //============================================================================== 
  static void SaveAscii(std::string fname,unsigned np,const tfloat3* pos,const tdouble3* posd
    ,unsigned nfields,const StScalarData* fields,std::string head="");
  
  //==============================================================================
  /// Returns the definition of statistics fields.
  //==============================================================================
  static StStatistics DefineStatsField(const std::string &name
    ,TpData type,unsigned comp,const void *pointer
    ,bool selmin=true,bool selmax=true,bool selmean=true,bool seltotal=true,bool selcomponent=true);

  //==============================================================================
  /// Checks the definition of fields.
  //==============================================================================
  static void CheckStats(unsigned nfields,const StStatistics* fields);

  //==============================================================================
  /// Calculates statistic information of arrays.
  //============================================================================== 
  static void CalculateStats(unsigned np,unsigned nfields,StStatistics* fields);

  //==============================================================================
  /// Calculates and save statistic information of arrays.
  //============================================================================== 
  static void SaveStats(std::string fname,bool csvsepcoma,bool firstdata,unsigned part
    ,double timestep,unsigned np,unsigned nfields,StStatistics* fields,std::string head);

  //==============================================================================
  /// Stores data of basic particles in VTK format (example code).
  //============================================================================== 
  static void SaveVtkBasic(std::string fname,unsigned np
    ,const tfloat3* pos,const unsigned* idp,const tfloat3* vel,const float* rhop)
  {
    StScalarData fields[3];
    unsigned nfields=0;
    if(idp){  fields[nfields]=DefineField("Idp" ,UInt32 ,1,idp);  nfields++; }
    if(vel){  fields[nfields]=DefineField("Vel" ,Float32,3,vel);  nfields++; }
    if(rhop){ fields[nfields]=DefineField("Rhop",Float32,1,rhop); nfields++; }
    SaveVtk(fname,np,pos,nfields,fields);
  }

  //==============================================================================
  /// Stores data in CSV format.
  //==============================================================================
  static void ParticlesToCsv(std::string fname,bool csvsepcoma,unsigned np
    ,unsigned nfixed,unsigned nmoving,unsigned nfloat,unsigned nfluid,unsigned nfluidout,float timestep
    ,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass
    ,const unsigned *id,const byte *type,const byte *mk,const tfloat3 *ace,const tfloat3 *vor);
  
  //==============================================================================
  /// Stores data in ASCII format.
  //==============================================================================  
  static void ParticlesToAscii(std::string fname,unsigned np
    ,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press
    ,const float *mass,const unsigned *id,const byte *type,const byte *mk
    ,const tfloat3 *ace,const tfloat3 *vor);
  
  //==============================================================================
  /// Stores data in CSV format (splits positive and negative part of Ace).
  //==============================================================================  
  static void ParticlesToCsv2(std::string fname,bool csvsepcoma,unsigned np
    ,unsigned nfixed,unsigned nmoving,unsigned nfloat,unsigned nfluid,unsigned nfluidout,float timestep
    ,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press,const float *mass
    ,const unsigned *id,const byte *type,const byte *mk,const tfloat3 *acepos,const tfloat3 *aceneg,const tfloat3 *vor);

  //==============================================================================
  /// Stores data in ASCII format (splits positive and negative part of Ace).
  //==============================================================================
  static void ParticlesToAscii2(std::string fname,unsigned np
    ,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press
    ,const float *mass,const unsigned *id,const byte *type,const byte *mk
    ,const tfloat3 *acepos,const tfloat3 *aceneg,const tfloat3 *vor);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void ParticlesToVtk(std::string fname,unsigned np
    ,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press
    ,const float *mass,const unsigned *id,const byte *type,const byte *mk
    ,const tfloat3 *ace,const tfloat3 *vor,int domain=0);

  //==============================================================================
  /// Stores data in VTK format (position is double).
  //============================================================================== 
  static void ParticlesToVtk(std::string fname,unsigned np
    ,const tdouble3 *pos,const tfloat3 *vel,const float *rhop,const float *press
    ,const float *mass,const unsigned *id,const byte *type,const byte *mk
    ,const tfloat3 *ace,const tfloat3 *vor,int domain=0);

  //==============================================================================
  /// Stores data in VTK format for variables of type float.
  //==============================================================================
  static void ParticlesToVtkFloat(std::string fname,unsigned np
    ,const tfloat3 *pos,const tfloat3 *vel,const float *rhop,const float *press
    ,const float *mass,const float *id,const float *type,const float *mk
    ,const tfloat3 *ace,const tfloat3 *vor);

  //==============================================================================
  /// Stores information of points in  VTK format.
  //==============================================================================
  static void SaveVtkPointsVar(std::string fname,unsigned np
    ,const tfloat3 *pos,const std::string &varname,const float *var);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type float.
  //==============================================================================
  static void SaveCsvPointsVar(const std::string &fname,bool csvsepcoma
    ,const std::string &dataname,int part,double timestep,unsigned np
    ,const tfloat3* pos,const float* data,bool first=false);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type float3.
  //==============================================================================
  static void SaveCsvPointsVar3(const std::string &fname,bool csvsepcoma
    ,const std::string &dataname,int part,double timestep,unsigned np
    ,const tfloat3* pos,const tfloat3* data,bool first=false);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type float.
  //==============================================================================
  static void SaveCsvPointsVar(const std::string &fname,bool csvsepcoma
    ,const std::string &dataname,int part,double timestep,unsigned np
    ,const tdouble3* pos,const float* data,bool first=false);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type double.
  //==============================================================================
  static void SaveCsvPointsVar(const std::string &fname,bool csvsepcoma
    ,const std::string &dataname,int part,double timestep,unsigned np
    ,const tdouble3* pos,const double* data,bool first=false);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type float3.
  //==============================================================================
  static void SaveCsvPointsVar3(const std::string &fname,bool csvsepcoma
    ,const std::string &dataname,int part,double timestep,unsigned np
    ,const tdouble3* pos,const tfloat3* data,bool first=false);

  //==============================================================================
  /// Stores information of points in CSV file for variables of type double3.
  //==============================================================================
  static void SaveCsvPointsVar3(const std::string &fname,bool csvsepcoma
    ,const std::string &dataname,int part,double timestep,unsigned np
    ,const tdouble3* pos,const tdouble3* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type float.
  //==============================================================================  
  static void SaveAscPointsVar(const std::string &fname,double timestep,unsigned np
    ,const tfloat3* pos,const float* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type float3.
  //============================================================================== 
  static void SaveAscPointsVar3(const std::string &fname,double timestep,unsigned np
    ,const tfloat3* pos,const tfloat3* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type float.
  //==============================================================================  
  static void SaveAscPointsVar(const std::string &fname,double timestep,unsigned np
    ,const tdouble3* pos,const float* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type double.
  //==============================================================================  
  static void SaveAscPointsVar(const std::string &fname,double timestep,unsigned np
    ,const tdouble3* pos,const double* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type float3.
  //============================================================================== 
  static void SaveAscPointsVar3(const std::string &fname,double timestep,unsigned np
    ,const tdouble3* pos,const tfloat3* data,bool first=false);

  //==============================================================================
  /// Stores information of points in ASCII file for variables of type double3.
  //============================================================================== 
  static void SaveAscPointsVar3(const std::string &fname,double timestep,unsigned np
    ,const tdouble3* pos,const tdouble3* data,bool first=false);

  //==============================================================================
  /// Stores time and position for predefined motion.
  //==============================================================================
  static void SaveMotionPredef(const std::string &fname,unsigned np,const float *time,const tfloat3 *pos);

  //==============================================================================
  /// Stores time and position for predefined motion.
  //==============================================================================
  static void SaveMotionPredef(const std::string &fname,unsigned np,const float *time,const float *pos);

  //==============================================================================
  /// Generates a VTK file with map cells.
  //==============================================================================
  static void SaveVtkCells(const std::string &fname,const tfloat3 &posmin,const tuint3 &cells,float scell);

  //==============================================================================
  /// Generates a VTK file with boxes.
  //==============================================================================
  static void SaveVtkBoxes(const std::string &fname,unsigned nbox,const tfloat3 *vbox,float sizemin=0);
};


//##############################################################################
//# JFormatFiles2Data
//##############################################################################
/// \brief Stores the definition of fields and implements some basic functionalities.

class JFormatFiles2Data : protected JObject
{
private:
  const bool CsvSepComa;
  unsigned Np;
  tfloat3  *Posf;
  tdouble3 *Posd;

  unsigned FieldsSize;
  unsigned FieldsCount;
  JFormatFiles2::StScalarData *Fields;

public:
  ///Constructor.
  JFormatFiles2Data(bool csvsepcoma):CsvSepComa(csvsepcoma){
    ClassName="JFormatFiles2Data";
    Posf=NULL;  Posd=NULL;  Fields=NULL;
    Reset();
  }
  ///Destructor.
  ~JFormatFiles2Data(){ Reset(); }
  ///Initialization of variables.
  void Reset(){
    Np=0;  Posf=NULL;  Posd=NULL;
    ResizeFields(0);
  }
  ///Adds position of particles (tfloat3).
  void SetPos(unsigned np,tfloat3 *pos){
    Np=np;  Posf=pos;  Posd=NULL;
  }
  ///Adds position of particles (tdouble3).
  void SetPos(unsigned np,tdouble3 *pos){
    Np=np;  Posf=NULL;  Posd=pos;
  }
  ///Initialization of variables.
  void ResizeFields(unsigned size){
    JFormatFiles2::StScalarData *fieldsold=Fields;
    Fields=NULL;
    FieldsSize=size;
    FieldsCount=(FieldsCount<FieldsSize? FieldsCount: FieldsSize);
    if(FieldsSize){
      Fields=new JFormatFiles2::StScalarData[FieldsSize];
      for(unsigned c=0;c<FieldsCount;c++)Fields[c]=fieldsold[c];
    }
    delete[] fieldsold;
  }
  ///Adds generic scalar defined by the user.
  void AddScalar(const std::string &name,JFormatFiles2::TpData type,unsigned comp,void *ptr){ 
    if(FieldsCount==FieldsSize)ResizeFields(FieldsSize+10);
    Fields[FieldsCount]=JFormatFiles2::DefineField(name,type,comp,ptr); 
    FieldsCount++;
  }
  void AddScalar(const std::string &name,byte     *ptr){ AddScalar(name,JFormatFiles2::UChar8  ,1,ptr); }
  void AddScalar(const std::string &name,word     *ptr){ AddScalar(name,JFormatFiles2::UShort16,1,ptr); }
  void AddScalar(const std::string &name,unsigned *ptr){ AddScalar(name,JFormatFiles2::UInt32  ,1,ptr); }
  void AddScalar(const std::string &name,float    *ptr){ AddScalar(name,JFormatFiles2::Float32 ,1,ptr); }
  void AddScalar(const std::string &name,double   *ptr){ AddScalar(name,JFormatFiles2::Double64,1,ptr); }
  void AddScalar(const std::string &name,tuint3   *ptr){ AddScalar(name,JFormatFiles2::UInt32  ,3,ptr); }
  void AddScalar(const std::string &name,tfloat3  *ptr){ AddScalar(name,JFormatFiles2::Float32 ,3,ptr); }
  void AddScalar(const std::string &name,tdouble3 *ptr){ AddScalar(name,JFormatFiles2::Double64,3,ptr); }
  ///Frees all pointers using delete[].
  void DeleteMemory(){
    delete[] Posf;  Posf=NULL;
    delete[] Posd;  Posd=NULL;
    for(unsigned c=0;c<FieldsCount;c++){
      JFormatFiles2::StScalarData dat=Fields[c];
           if(dat.type==JFormatFiles2::UChar8   && dat.comp==1)delete[] (byte    *)dat.pointer;
      else if(dat.type==JFormatFiles2::UShort16 && dat.comp==1)delete[] (word    *)dat.pointer;
      else if(dat.type==JFormatFiles2::UInt32   && dat.comp==1)delete[] (unsigned*)dat.pointer;
      else if(dat.type==JFormatFiles2::Float32  && dat.comp==1)delete[] (float   *)dat.pointer;
      else if(dat.type==JFormatFiles2::Double64 && dat.comp==1)delete[] (double  *)dat.pointer;
      else if(dat.type==JFormatFiles2::UInt32   && dat.comp==3)delete[] (tuint3  *)dat.pointer;
      else if(dat.type==JFormatFiles2::Float32  && dat.comp==3)delete[] (tfloat3 *)dat.pointer;
      else if(dat.type==JFormatFiles2::Double64 && dat.comp==3)delete[] (tdouble3*)dat.pointer;
      else RunException("FreeMemory","Type of pointer is unknown.");
    }
    Reset();
  }
  ///Saves VTK file.
  void SaveVtk(const std::string &file){
    if(Posf)JFormatFiles2::SaveVtk(file,Np,Posf,FieldsCount,Fields);
    if(Posd)JFormatFiles2::SaveVtk(file,Np,Posd,FieldsCount,Fields);
  }
  ///Saves CSV file.
  void SaveCsv(const std::string &file){
    if(Posf)JFormatFiles2::SaveCsv(file,CsvSepComa,Np,Posf,FieldsCount,Fields);
    if(Posd)JFormatFiles2::SaveCsv(file,CsvSepComa,Np,Posd,FieldsCount,Fields);
  }
};


#endif


