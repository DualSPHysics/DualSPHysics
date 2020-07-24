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

//:#############################################################################
//:# Cambios:
//:# =========
//:# - Creacion de clase para gestionar almacenar normales de particulas contorno. (07-04-2019)
//:# - Graba el tipo normales almacenadas en PartNormals[]. (05-06-2019)
//:# - Cambia nombre de fichero con datos de normales de x_Normals.nbi4 a 
//:#   x_NormalData.nbi4. (30-07-2019)
//:# - Error corregido al cargar PartNormals en LoadFile(). (06-08-2019)
//:# - Error corregido en AddNormalData(). (06-08-2019)
//:# - Permite almacenar solo las normales finales. (07-08-2019)
//:# - Usa macro para lanzar excepciones. (12-12-2019)
//:#############################################################################

/// \file JPartNormalData.h \brief Declares the class \ref JPartNormalData.

#ifndef _JPartNormalData_
#define _JPartNormalData_

#include <string>
#include <vector>
#include "JObject.h"
#include "TypesDef.h"
#include "JBinaryData.h"

//##############################################################################
//# JPartNormalData
//##############################################################################
/// \brief Manages the information of normals in boundary particles.

class JPartNormalData : protected JObject
{
private:
  static const unsigned FmtVersionDef=190407;    ///<Version de formato by default. Version of format by default.
  unsigned FmtVersion;    ///<Version de formato. Version of format.

  std::string DirData;

  //-General variables.
  std::string AppName;
  std::string Date;

  std::string CaseName;
  bool Data2d;           ///<Toggles 2D simulation.
  double Data2dPosY;     ///<Y value in 2D simulations.
  double Dp;
  double H;
  double Dist;           ///<Distance used for calculating normal data (tipically KernelSize).
  bool FtSupport;        ///<Enables support for floating bodies.

  std::string PartNormalsName; ///<Name of normals approach used for PartNormals (Mean or Marrone).
  unsigned Nbound;       ///<Number boudary particles.
  tdouble3 *PartNormals; ///<Final normals of particles [Nbound]

  unsigned CountNormals; ///<Total number of normals.
  unsigned *NormalBegin; ///<Normals for each particle [Nbound+1]
  tdouble3 *Normals;     ///<Unitary normal to shape [CountNormals]
  double   *NormalsDist;  ///<Distance to shape [CountNormals]
  tdouble3 *OutVecs;     ///<Unitary vector to limit of shape [CountNormals]
  double   *OutVecsDist;  ///<Distance to limit of shape [CountNormals]

  void AllocNormals(unsigned nbound,unsigned countnor,bool usepartnormals);
  bool ArrayExists(JBinaryData *bd,std::string name)const;
  JBinaryDataArray* GetArray(JBinaryData *bd,std::string name)const;
  JBinaryDataArray* GetArray(JBinaryData *bd,std::string name,JBinaryDataDef::TpData type)const;

public:
  JPartNormalData();
  ~JPartNormalData();
  void Reset();
  void ConfigBasic(std::string appname,std::string casename
    ,bool data2d,double data2dposy,double dp,double h,double dist,bool ftsupport);
  void AddNormalData(
    std::string partnorname,unsigned nbound,const tdouble3 *partnor
    ,const unsigned *norbegin,unsigned countnor
    ,const tdouble3 *nordata,const double *nordist
    ,const tdouble3 *outdata,const double *outdist);
  void AddNormalData(
    std::string partnorname,unsigned nbound,const tdouble3 *partnor);

  static std::string GetFileName(std::string casename,std::string dir="");
  static std::string GetNormalDataFile(std::string casename);

  void LoadFile(std::string casename);
  void SaveFile(std::string dir);

//-Returns general values.
//--------------------------
  std::string GetAppName()   const{ return(AppName);    }
  std::string GetDate()      const{ return(Date);       }
  std::string GetCaseName()  const{ return(CaseName);   }
  bool        GetData2d()    const{ return(Data2d);     }
  double      GetData2dPosY()const{ return(Data2dPosY); }
  double      GetDp()        const{ return(Dp);         }
  double      GetH()         const{ return(H);          }

//-Returns normal data.
//-----------------------
  double      GetDist()      const{ return(Dist);       }
  bool        GetFtSupport() const{ return(FtSupport);  }
  unsigned    GetNbound()    const{ return(Nbound);     }
  std::string GetPartNormalsName()const{ return(PartNormalsName); }
  const tdouble3* GetPartNormals()const{ return(PartNormals);     }

  const unsigned* GetNormalBegin()const{ return(NormalBegin); }
  unsigned   GetCountNormals()const{ return(CountNormals); }
  const tdouble3* GetNormals()    const{ return(Normals);     }
  const double*   GetNormalsDist()const{ return(NormalsDist); }
  const tdouble3* GetOutVecs()    const{ return(OutVecs);     }
  const double*   GetOutVecsDist()const{ return(OutVecsDist); }


};

#endif


