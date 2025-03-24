/*===================================================================================
<MOORDYNPLUS> Copyright (c) 2025 by Ivan Martinez-Estevez

Ivan Martinez-Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
Matt Hall (github.com/mattEhall)

This file is part of MoorDynPlus. MoorDynPlus is free software: you can 
redistribute it and/or modify it under the terms of the GNU General Public 
License as published by the Free Software Foundation, either version 3 of 
the License, or (at your option) any later version.

Linking the MoorDynPlus library statically or dynamically with other modules is
making a combined work based on this library. Thus, the terms and conditions
of the GNU General Public License cover the whole combination. As a special
exception, the copyright holders of MoorDynPlus give you permission to dynamically
link this library with the program DualSPHysics to produce a combined model
featuring the capabilities of both DualSPHysics and MoorDynPlus. This exception
is strictly limited to linking between the compiled MoorDynPlus library and
DualSPHysics. It does not extend to other programs or the use of the MoorDynPlus
source code beyond the stipulations of the GPL. When the exception is used,
this paragraph must be included in the copyright notice.

MoorDynPlus is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for details.

You should have received a copy of the GNU General Public License along with
MoorDynPlus. If not, see <http://www.gnu.org/licenses/>.
===================================================================================*/

/// \file IOutput.h \brief Defines the class \ref IOutput.

#ifndef _Output_
#define _Output_

#include "JObject.h"
#include "JXml.h"
#include "ILine.h"
#include "IConnection.h"
#include "TypesMoorDynPlus.h"


class IOutProperties : protected JObject
{
private:
  TpMagnitude Mag;            ///<IOutput magnitude [tension,force,velocity,position]
  std::vector<ILine*> Lines; ///<Array of pointers of lines selected
  
  unsigned LineCount;       ///<Num of lines selected
  std::string Units;        ///<Pointer to units characters
  unsigned NodeID;          ///<Node identifier [0,1,2,...,N]
 
  void Reset(); ///<Restores attributes

public:
  IOutProperties();
  virtual ~IOutProperties();
 
  unsigned           GetNodeID()const   { return NodeID;    }; ///<Returns the Id of Node
  TpMagnitude        GetMagnitude()const{ return Mag;       }; ///<Returns the TpMagnitude 
  std::vector<ILine*>GetLines()const    { return Lines;     }; ///<Returns an array of Lines
  std::string        GetUnits()const    { return Units;     }; ///<Returns a pointer to Units
  unsigned           GetNLines()const   { return LineCount; }; ///<Returns the number of lines
 
  void SetNodeID   (const unsigned id)        { NodeID=id;     }; ///<Sets the Id of Node
  void SetMagnitude(const TpMagnitude mag)    { Mag=mag;       }; ///<Sets the Magnitude  
  void SetLines    (std::vector<ILine*> lines){ Lines=lines;   }; ///<Sets the Lines
  void SetNLines   (const unsigned num)       { LineCount=num; }; ///<Sets the number of lines
  void SetUnits    (const std::string units)  { Units=units;   }; ///<Sets the units
};

class IOutput : protected JObject
{
private:
  unsigned LineCount;   ///<Number of lines
  unsigned TypeCount;   ///<Number of types for write
  double DtOut;         ///<Time step for write
  double TimeStart;     ///<Start time for begin to write
  double TimeMax;       ///<Time of simulation
  bool IsFirst;         ///<IsFirst time to save data
  double NextTime;      ///<Next time to save data
  std::string FileName; ///<Name of the csv file to save data
  std::string DataDir;  ///<Path of the folder to save data

  std::vector<IOutProperties*> OutProps; ///<Vector of output properties

  void Reset();
  void ReadXml(JXml* sxml,TiXmlElement* ele,std::vector<ILine*> lines); 
  void SetLinesSelected(IOutProperties* oProps,std::vector<ILine*> lines); 

public:
  IOutput(double timemax,double dtout);
  virtual ~IOutput();
  void LoadXml(JXml* sxml,const std::string& place,std::vector<ILine*> lines);
  void Setup(std::string dir);
  void SaveCsv(double timestep,double dt);
  void SaveCsvTen(const std::string path,const IOutProperties* out_p,double timestep,double dt);
  void SaveCsvFor(const std::string path,const IOutProperties* out_p,double timestep,double dt);
  void SaveCsvVel(const std::string path,const IOutProperties* out_p,double timestep,double dt);
  void SaveCsvPos(const std::string path,const IOutProperties* out_p,double timestep,double dt);

  unsigned GetLineCount()const { return LineCount; }; ///<Returns the number of lines                    
  unsigned GetTypeCount()const { return TypeCount; }; ///<Returns the number of types    
  double   GetDtOut()const     { return DtOut;     }; ///<Returns the DtOut
  double   GetStartTime()const { return TimeStart; }; ///<Returns the StartTime
  double   GetEndTime()const   { return TimeMax;   }; ///<Returns the EndTime
  double   GetNextTime()const  { return NextTime;  }; ///<Returns the next time
  std::vector<IOutProperties*> GetOutProps()const{ return OutProps; }; ///<Returns an array of IOutProperties  

  void SetNLines(const unsigned numLines) { LineCount=numLines; }; ///<Sets the number of Lines    
};
#endif //!IOutput
