//HEAD_DSCODES
/*
 <DUALSPHYSICS>  Copyright (c) 2018 by Dr Jose M. Dominguez et al. (see http://dual.sphysics.org/index.php/developers/). 

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
//:# - New functions to create shapes in VTK files. (25-01-2018)
//:# - Functions to create VTK/CSV files starting from vector<StScalarData>. (25-01-2018)
//:# - Se escriben las unidades en las cabeceras de los ficheros CSV. (26-04-2018)
//:# - Nuevo metodo AddShape_Sphere() para generar esferas a partir de quads. (06-07-2018)
//:# - Nuevos metodos CreateShapesMk(), DeleteShapesMk() y CreateOBJsByMk(). (06-08-2018)
//:# - Nuevos metodos AddShape_Cylinder(), AddShape_Cross(). (10-08-2018)
//:# - Nuevos metodos AddShape_Circle(), AddShape_Spring(). (13-08-2018)
//:# - El metodo SaveVtkShapes() combina lineas consecutivas con mismos valores. (13-08-2018)
//:# - Se permite definir formato de salida y unidades en StScalarData. (12-09-2018)
//:# - Nuevo metodo DeleteFields() para liberar memoria dinamica. (12-09-2018)
//:# - Se elimina codigo JFormatFiles2Data por falta de uso. (12-09-2018)
//:# - Nuevos metodos AddShape_TriangleLines(), AddShape_QuadLines(). (21-12-2018)
//:#############################################################################
  

/// \file JFormatFiles2.h \brief Declares the class \ref JFormatFiles2.

#ifndef _JFormatFiles2_
#define _JFormatFiles2_

#include "TypesDef.h"
#include "JObject.h"
#include <string>
#include <cstring>
#include <string>
#include <vector>

//##############################################################################
//# JFormatFiles2
//##############################################################################
/// \brief Provides functions to store particle data in formats VTK, CSV, ASCII and geometric elements in VTK format.

class JFormatFiles2
{
public:

  /// Modes to define the normals.
  typedef enum{ NorNULL,NorOriginal,NorInvert,NorTwoFace }TpModeNormal; 

  /// Data types.
  typedef enum{ UChar8,Char8,UShort16,Short16,UInt32,Int32,Float32,Double64,ULlong64,Llong64,TpDataNull }TpData;

  /// Structure with the information of an array of particle data to be stored in CSV or VTK format.
  typedef struct {
    std::string name; //"name:outputformat:units", e.g. "velocity:%f:m/s"
    std::string fmt;
    std::string units;
    TpData type;
    unsigned comp;
    const void *pointer;
    bool delptr;
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

  /// Types of shape.
  typedef enum{ 
     ShLine
    ,ShTriangle
    ,ShQuad
    ,ShBox
    ,ShSphere
    ,ShNull 
  }TpShape;

  /// Structure with data of one shape to be stored in VTK format.
  typedef struct StrShapeData{
    TpShape tshape;
    tfloat3 vpt[8];
    int value;
    float valuef;
    StrShapeData(){ reset(); }
    StrShapeData(TpShape tsh,int vi,float vf,tfloat3 p0=TFloat3(0)
      ,tfloat3 p1=TFloat3(0),tfloat3 p2=TFloat3(0),tfloat3 p3=TFloat3(0)
      ,tfloat3 p4=TFloat3(0),tfloat3 p5=TFloat3(0),tfloat3 p6=TFloat3(0)
      ,tfloat3 p7=TFloat3(0))
    { 
      reset(); tshape=tsh; value=vi; valuef=vf; vpt[0]=p0; vpt[1]=p1; 
      vpt[2]=p2; vpt[3]=p3; vpt[4]=p4; vpt[5]=p5; vpt[6]=p6; vpt[7]=p7;
    }
    void reset(){ 
      tshape=ShNull; for(unsigned c=0;c<8;c++)vpt[c]=TFloat3(0);
      value=0; valuef=0; 
    }
  }StShapeData;

  /// Structure with parameters to create spring.
  typedef struct StrShapeSpring{
    float cornersout; //-Size of corner.
    float cornersin;  //-Size of corner (inside).
    float radius;     //-Spring radius.
    float length;     //-Length for each revolution.
    int nside;        //-Number of sections for each revolution.

    StrShapeSpring(){ reset(); }
    StrShapeSpring(float xcornersout,float xcornersin,float xradius,float xlength,int xnside)
    { 
      reset(); cornersout=xcornersout; cornersin=xcornersin; 
      radius=xradius; length=xlength; nside=xnside;
    }
    void reset(){ 
      cornersout=0.75f;
      cornersin=0.25f;
      radius=3;
      length=1.f;
      nside=16;
    }
  }StShapeSpring;


  //==============================================================================
  /// Throws a simple exception.
  //==============================================================================
  static void RunException(std::string method,std::string msg);
  
  //==============================================================================
  /// Throws an exception related to a file.
  //==============================================================================  
  static void RunException(std::string method,std::string msg,std::string file);


  //==============================================================================
  /// Returns units according variable name. E.g.: Vel -> " [m/s]"
  //==============================================================================
  static std::string GetUnits(const std::string &varname);

  //==============================================================================
  /// Defines automatically the output format and units of field.
  //==============================================================================
  static void DefineFieldFormat(StScalarData &field);

  //==============================================================================
  /// Returns the definition of fields.
  //==============================================================================
  static StScalarData DefineField(const std::string &name,TpData type,unsigned comp,const void *pointer=NULL){
    StScalarData f; f.name=name; f.type=type; f.comp=comp; f.pointer=pointer; f.delptr=false;
    f.fmt=f.units="";
    DefineFieldFormat(f);
    return(f);
  }

  //==============================================================================
  /// Returns the definition of fields.
  //==============================================================================
  static StScalarData DefineFieldDel(const std::string &name,TpData type,unsigned comp,const void *pointer=NULL){
    StScalarData f; f.name=name; f.type=type; f.comp=comp; f.pointer=pointer; f.delptr=true;
    f.fmt=f.units="";
    DefineFieldFormat(f);
    return(f);
  }

  //==============================================================================
  /// Checks the definition of fields.
  //==============================================================================
  static void CheckFields(unsigned nfields,const StScalarData* fields);

  //==============================================================================
  /// Checks the definition of fields.
  //==============================================================================
  static void CheckFields(const std::vector<StScalarData> &fields);

  //==============================================================================
  /// Delete dynamic memory of pointers with delptr=true.
  //==============================================================================
  static void DeleteFields(std::vector<StScalarData> &fields);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void SaveVtk(std::string fname,unsigned np
    ,const tfloat3* pos,unsigned nfields,const StScalarData* fields);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void SaveVtk(std::string fname,unsigned np
    ,const tfloat3* pos,const std::vector<StScalarData> &fields);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void SaveVtk(std::string fname,unsigned np
    ,const tdouble3* pos,unsigned nfields,const StScalarData* fields);

  //==============================================================================
  /// Stores data in VTK format.
  //============================================================================== 
  static void SaveVtk(std::string fname,unsigned np
    ,const tdouble3* pos,const std::vector<StScalarData> &fields);
  
  //==============================================================================
  /// Stores data in CSV format.
  //============================================================================== 
  static void SaveCsv(std::string fname,bool csvsepcoma,unsigned np
    ,const tfloat3* pos,unsigned nfields,const StScalarData* fields,std::string head="");
  
  //==============================================================================
  /// Stores data in CSV format.
  //============================================================================== 
  static void SaveCsv(std::string fname,bool csvsepcoma,unsigned np
    ,const tfloat3* pos,const std::vector<StScalarData> &fields,std::string head="");
  
  //==============================================================================
  /// Stores data in CSV format.
  //============================================================================== 
  static void SaveCsv(std::string fname,bool csvsepcoma,unsigned np
    ,const tdouble3* pos,unsigned nfields,const StScalarData* fields,std::string head="");
  
  //==============================================================================
  /// Stores data in CSV format.
  //============================================================================== 
  static void SaveCsv(std::string fname,bool csvsepcoma,unsigned np
    ,const tdouble3* pos,const std::vector<StScalarData> &fields,std::string head="");

  //==============================================================================
  /// Stores data in ASCII format.
  //============================================================================== 
  static void SaveAscii(std::string fname,unsigned np,const tfloat3* pos,const tdouble3* posd
    ,unsigned nfields,const StScalarData* fields,std::string head="");

  //==============================================================================
  /// Stores data in ASCII format.
  //============================================================================== 
  static void SaveAscii(std::string fname,unsigned np,const tfloat3* pos,const tdouble3* posd
    ,const std::vector<StScalarData> &fields,std::string head="");
  
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
    std::vector<StScalarData> fields;
    if(idp) fields.push_back(DefineField("Idp" ,UInt32 ,1,idp));
    if(vel) fields.push_back(DefineField("Vel" ,Float32,3,vel));
    if(rhop)fields.push_back(DefineField("Rhop",Float32,1,rhop));
    SaveVtk(fname,np,pos,fields);
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


  //##############################################################################
  //# Functions to create shapes in VTK files.
  //##############################################################################
  //==============================================================================
  /// Returns the definition of shape line.
  //==============================================================================
  static StShapeData DefineShape_Line(const tfloat3 &pt1,const tfloat3 &pt2,int value,float valuef){
    return(StrShapeData(ShLine,value,valuef,pt1,pt2));
  }

  //==============================================================================
  /// Returns the definition of shape line.
  //==============================================================================
  static StShapeData DefineShape_Line(const tdouble3 &pt1,const tdouble3 &pt2,int value,float valuef){
    return(DefineShape_Line(ToTFloat3(pt1),ToTFloat3(pt2),value,valuef));
  }

  //==============================================================================
  /// Returns the definition of shape triangle.
  //==============================================================================
  static StShapeData DefineShape_Triangle(const tfloat3 &pt1,const tfloat3 &pt2
    ,const tfloat3 &pt3,int value,float valuef)
  {
    return(StrShapeData(ShTriangle,value,valuef,pt1,pt2,pt3));
  }

  //==============================================================================
  /// Returns the definition of shape triangle.
  //==============================================================================
  static StShapeData DefineShape_Triangle(const tdouble3 &pt1,const tdouble3 &pt2
    ,const tdouble3 &pt3,int value,float valuef)
  {
    return(StrShapeData(ShTriangle,value,valuef,ToTFloat3(pt1),ToTFloat3(pt2),ToTFloat3(pt3)));
  }

  //==============================================================================
  /// Returns the definition of triangle using lines.
  //==============================================================================
  static void AddShape_TriangleLines(std::vector<StShapeData> &shapes
    ,const tfloat3 &pt1,const tfloat3 &pt2,const tfloat3 &pt3,int value,float valuef);
  
  //==============================================================================
  /// Returns the definition of triangle using lines.
  //==============================================================================
  static void AddShape_TriangleLines(std::vector<StShapeData> &shapes
    ,const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,int value,float valuef)
  {
    AddShape_TriangleLines(shapes,ToTFloat3(pt1),ToTFloat3(pt2),ToTFloat3(pt3),value,valuef);
  }

  //==============================================================================
  /// Returns the definition of shape quad.
  //==============================================================================
  static StShapeData DefineShape_Quad(const tfloat3 &pt1,const tfloat3 &pt2
    ,const tfloat3 &pt3,const tfloat3 &pt4,int value,float valuef)
  {
    return(StrShapeData(ShQuad,value,valuef,pt1,pt2,pt3,pt4));
  }

  //==============================================================================
  /// Returns the definition of shape quad.
  //==============================================================================
  static StShapeData DefineShape_Quad(const tdouble3 &pt1,const tdouble3 &pt2
    ,const tdouble3 &pt3,const tdouble3 &pt4,int value,float valuef)
  {
    return(StrShapeData(ShQuad,value,valuef,ToTFloat3(pt1),ToTFloat3(pt2),ToTFloat3(pt3),ToTFloat3(pt4)));
  }

  //==============================================================================
  /// Returns the definition of shape quad.
  //==============================================================================
  static StShapeData DefineShape_Quad(const tfloat3 &pt,const tfloat3 &vec,float size,int value,float valuef);
  
  //==============================================================================
  /// Returns the definition of quad using lines.
  //==============================================================================
  static StShapeData DefineShape_Quad(const tdouble3 &pt,const tdouble3 &vec,double size,int value,float valuef)
  {
    return(DefineShape_Quad(ToTFloat3(pt),ToTFloat3(vec),float(size),value,valuef));
  }

  //==============================================================================
  /// Returns the definition of quad using lines.
  //==============================================================================
  static void AddShape_QuadLines(std::vector<StShapeData> &shapes
    ,const tfloat3 &pt1,const tfloat3 &pt2,const tfloat3 &pt3,const tfloat3 &pt4,int value,float valuef);
  
  //==============================================================================
  /// Returns the definition of quad using lines.
  //==============================================================================
  static void AddShape_QuadLines(std::vector<StShapeData> &shapes
    ,const tdouble3 &pt1,const tdouble3 &pt2,const tdouble3 &pt3,const tdouble3 &pt4,int value,float valuef)
  {
    AddShape_QuadLines(shapes,ToTFloat3(pt1),ToTFloat3(pt2),ToTFloat3(pt3),ToTFloat3(pt4),value,valuef);
  }

  //==============================================================================
  /// Returns the definition of quad using lines.
  //==============================================================================
  static void AddShape_QuadLines(std::vector<StShapeData> &shapes
    ,const tfloat3 &pt,const tfloat3 &vec,float size,int value,float valuef);
  
  //==============================================================================
  /// Returns the definition of quad using lines.
  //==============================================================================
  static void AddShape_QuadLines(std::vector<StShapeData> &shapes
    ,const tdouble3 &pt,const tdouble3 &vec,double size,int value,float valuef)
  {
    AddShape_QuadLines(shapes,ToTFloat3(pt),ToTFloat3(vec),float(size),value,valuef);
  }

  //==============================================================================
  /// Returns the definition of shape box (8pt).
  /// pt1: left  - front - bottom.
  /// pt2: right - front - bottom.
  /// pt3: left  - front - top.
  /// pt4: right - front - top.
  /// pt5: left  - back  - bottom.
  /// pt6: right - back  - bottom.
  /// pt7: left  - back  - top.
  /// pt8: right - back  - top.
  //==============================================================================
  static StShapeData DefineShape_Box(const tfloat3 &pt1,const tfloat3 &pt2
    ,const tfloat3 &pt3,const tfloat3 &pt4,const tfloat3 &pt5,const tfloat3 &pt6
    ,const tfloat3 &pt7,const tfloat3 &pt8,int value,float valuef)
  {
    return(StrShapeData(ShBox,value,valuef,pt1,pt2,pt3,pt4,pt5,pt6,pt7,pt8));
  }
  //==============================================================================
  /// Returns the definition of shape box (8pt).
  //==============================================================================
  static StShapeData DefineShape_Box(const tdouble3 &pt1,const tdouble3 &pt2
    ,const tdouble3 &pt3,const tdouble3 &pt4,const tdouble3 &pt5,const tdouble3 &pt6
    ,const tdouble3 &pt7,const tdouble3 &pt8,int value,float valuef)
  {
    return(DefineShape_Box(ToTFloat3(pt1),ToTFloat3(pt2),ToTFloat3(pt3),ToTFloat3(pt4),ToTFloat3(pt5),ToTFloat3(pt6),ToTFloat3(pt7),ToTFloat3(pt8),value,valuef));
  }

  //==============================================================================
  /// Returns the definition of shape box (1pt + size3).
  //==============================================================================
  static StShapeData DefineShape_Box(const tfloat3 &pt1,const tfloat3 &size
    ,int value,float valuef)
  {
    const tfloat3 p0=pt1;
    const tfloat3 p1=p0+TFloat3(size.x,0,0);
    const tfloat3 p2=p1+TFloat3(0,0,size.z);
    const tfloat3 p3=p0+TFloat3(0,0,size.z);
    const tfloat3 p4=p0+TFloat3(0,size.y,0);
    const tfloat3 p5=p1+TFloat3(0,size.y,0);
    const tfloat3 p6=p2+TFloat3(0,size.y,0);
    const tfloat3 p7=p3+TFloat3(0,size.y,0);
    return(StrShapeData(ShBox,value,valuef,p0,p1,p2,p3,p4,p5,p6,p7));
  }
    //==============================================================================
  /// Returns the definition of shape box (1pt + size3).
  //==============================================================================
  static StShapeData DefineShape_Box(const tdouble3 &pt1,const tdouble3 &size
    ,int value,float valuef)
  {
    return(DefineShape_Box(ToTFloat3(pt1),ToTFloat3(size),value,valuef));
  }

  //==============================================================================
  /// Returns the definition of shape box (1pt + 3x vec3).
  //==============================================================================
  static StShapeData DefineShape_Box(const tfloat3 &pt1,const tfloat3 &vx
    ,const tfloat3 &vy,const tfloat3 &vz,int value,float valuef)
  {
    const tfloat3 p0=pt1;
    const tfloat3 p1=p0+vx;
    const tfloat3 p2=p1+vz;
    const tfloat3 p3=p0+vz;
    const tfloat3 p4=p0+vy;
    const tfloat3 p5=p1+vy;
    const tfloat3 p6=p2+vy;
    const tfloat3 p7=p3+vy;
    return(StrShapeData(ShBox,value,valuef,p0,p1,p2,p3,p4,p5,p6,p7));
  }
  //==============================================================================
  /// Returns the definition of shape box (1pt + 3x vec3).
  //==============================================================================
  static StShapeData DefineShape_Box(const tdouble3 &pt1,const tdouble3 &vx
    ,const tdouble3 &vy,const tdouble3 &vz,int value,float valuef)
  {
    return(DefineShape_Box(ToTFloat3(pt1),ToTFloat3(vx),ToTFloat3(vy),ToTFloat3(vz),value,valuef));
  }

  //==============================================================================
  /// Adds triangles/lines for the definition of a circle/circumference.
  //==============================================================================
  static void AddShape_Circle(std::vector<StShapeData> &shapes,bool circle
    ,const tfloat3 &pt,float radius,const tfloat3 &vec,int nside,int value,float valuef);
  //==============================================================================
  /// Adds triangles/lines for the definition of a circle/circumference.
  //==============================================================================
  static void AddShape_Circle(std::vector<StShapeData> &shapes,bool circle
    ,const tdouble3 &pt,double radius,const tdouble3 &vec,int nside,int value,float valuef)
  {
    AddShape_Circle(shapes,circle,ToTFloat3(pt),float(radius),ToTFloat3(vec),nside,value,valuef);
  }

  //==============================================================================
  /// Adds quads for the definition of a sphere.
  //==============================================================================
  static void AddShape_Sphere(std::vector<StShapeData> &shapes
    ,const tfloat3 &pt,float radius,int nside,int value,float valuef);
  //==============================================================================
  /// Adds quads for the definition of a sphere.
  //==============================================================================
  static void AddShape_Sphere(std::vector<StShapeData> &shapes
    ,const tdouble3 &pt,double radius,int nside,int value,float valuef)
  {
    AddShape_Sphere(shapes,ToTFloat3(pt),float(radius),nside,value,valuef);
  }

  //==============================================================================
  /// Adds triangles and quads for the definition of a cylinder.
  //==============================================================================
  static void AddShape_Cylinder(std::vector<JFormatFiles2::StShapeData> &shapes
    ,const tfloat3 &p1,const tfloat3 &p2,float radius,int nside,unsigned maskfaceshide,int value,float valuef);
  //==============================================================================
  /// Adds triangles and quads for the definition of a cylinder.
  //==============================================================================
  static void AddShape_Cylinder(std::vector<StShapeData> &shapes
    ,const tdouble3 &p1,const tdouble3 &p2,double radius,int nside,unsigned maskfaceshide,int value,float valuef)
  {
    AddShape_Cylinder(shapes,ToTFloat3(p1),ToTFloat3(p2),float(radius),nside,maskfaceshide,value,valuef);
  }

  //==============================================================================
  /// Adds lines for the definition of a cross.
  //==============================================================================
  static void AddShape_Cross(std::vector<StShapeData> &shapes
    ,const tfloat3 &pt,float radius,int value,float valuef)
  {
    shapes.push_back(DefineShape_Line(TFloat3(pt.x-radius,pt.y,pt.z),TFloat3(pt.x+radius,pt.y,pt.z),value,valuef));
    shapes.push_back(DefineShape_Line(TFloat3(pt.x,pt.y-radius,pt.z),TFloat3(pt.x,pt.y+radius,pt.z),value,valuef));
    shapes.push_back(DefineShape_Line(TFloat3(pt.x,pt.y,pt.z-radius),TFloat3(pt.x,pt.y,pt.z+radius),value,valuef));
  }
  //==============================================================================
  /// Adds lines for the definition of a cross.
  //==============================================================================
  static void AddShape_Cross(std::vector<StShapeData> &shapes
    ,const tdouble3 &pt,double radius,int value,float valuef)
  {
    AddShape_Cross(shapes,ToTFloat3(pt),float(radius),value,valuef);
  }
  
  //==============================================================================
  /// Adds lines for the definition of a spring.
  //==============================================================================
  static void AddShape_Spring(std::vector<StShapeData> &shapes
    ,const tfloat3 &pt1,const tfloat3 &p2,float restlength,float scalesize
    ,StShapeSpring config,int value,float valuef);
  //==============================================================================
  /// Adds lines for the definition of a spring.
  //==============================================================================
  static void AddShape_Spring(std::vector<StShapeData> &shapes
    ,const tdouble3 &pt1,const tdouble3 &pt2,double restlength,double scalesize
    ,StShapeSpring config,int value,float valuef)
  {
    AddShape_Spring(shapes,ToTFloat3(pt1),ToTFloat3(pt2),float(restlength),float(scalesize),config,value,valuef);
  }

  //==============================================================================
  /// Generates a VTK file with shapes.
  //============================================================================== 
  static void SaveVtkShapes(std::string fname,const std::string &valuename
    ,const std::string &valuefname,const std::vector<StShapeData> &shapes);


  //##############################################################################
  //# Functions to create OBJ files starting from VTK files.
  //##############################################################################
  //==============================================================================
  /// Creates object with geometry (triangles and quads) and mk data from VTK files.
  //==============================================================================
  static void* CreateShapesMk(const std::vector<std::string> &vtkfiles);

  //==============================================================================
  /// Frees object with geometry and mk data from VTK files.
  //==============================================================================
  static void DeleteShapesMk(void* ptr_vtksimple);


  //==============================================================================
  /// Creates OBJ file with MK geometry in VTK file. Returns not zero in case of error.
  //==============================================================================
  static void CreateOBJsByMk(void* ptr_vtksimple,std::string filein,std::string filesout
    ,const std::vector<unsigned> &mkbounds,unsigned mkboundfirst,TpModeNormal normalmode);

};

/*
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
  ~JFormatFiles2Data(){ DestructorActive=true; Reset(); }
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
*/

#endif


