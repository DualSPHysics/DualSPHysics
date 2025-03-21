//HEAD_DSPH
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

/// \file DualSphDef.h \brief Includes definitions and specific types for DualSPHysics program.

#ifndef _DualSphDef_
#define _DualSphDef_

#include "TypesDef.h"
#include "JAppInfoDef.h"
#include "JParticlesDef.h"
#include "JPeriodicDef.h"
#include "JViscosityDef.h"
#include "FunSphKernelDef.h"
#include "JMeshDataDef.h" //<vs_meeshdat>
#include "JDsDcellDef.h"  //-Includes definitions for cell codificantion as unsigned (32 bits) value.

#include "OmpDefs.h"
#include <algorithm>



// Developing features to disable by compilation directives.
//------------------------------------------------------------

//------------------------------------------------------------


//#define DISABLE_KERNELS_EXTRA  ///<Compiles Wendland kernel and ignores the rest of the SPH kernels.

//#define DISABLE_TIMERS     ///<Compiles without timers. | Compilado sin timers.
//#define DISABLE_BSMODES    ///<compiles without advanced BlockSize modes.


//-Removes dependencies from precompiled libraries.
//#define DISABLE_CHRONO     ///<It allows compile without ChronoLib library (dsphchrono.dll, ChronoEngine.dll and ChronoEngine_parallel.dll).
#define DISABLE_CHRONO_OMP   ///<It allows compile without parallel module of Chrono (ignores ChronoEngine_parallel.dll).
//#define DISABLE_WAVEGEN    ///<It allows compile without Wave-Paddles, Multi-Layer Pistons and Relaxation Zones libraries.
//#define DISABLE_MOORDYNPLUS    ///<It allows compile without LibDSphMoorDynPlus library.


//-Defines AVAILABLE_CHRONO and AVAILABLE_CHRONO_OMP when these features are compiled.
#ifdef DISABLE_CHRONO
  #ifndef DISABLE_CHRONO_OMP
	  #define DISABLE_CHRONO_OMP
  #endif
  #define AVAILABLE_CHRONO false
  #define AVAILABLE_CHRONO_OMP false
#else
  #define AVAILABLE_CHRONO true
  #ifdef DISABLE_CHRONO_OMP
    #define AVAILABLE_CHRONO_OMP false
  #else
    #define AVAILABLE_CHRONO_OMP true
  #endif  
#endif

//-Defines AVAILABLE_WAVEGEN when this feature is compiled.
#ifdef DISABLE_WAVEGEN
  #define AVAILABLE_WAVEGEN false
  #define DISABLE_MLPISTONS
  #define DISABLE_RZ
#else
  #define AVAILABLE_WAVEGEN true
#endif

//-Defines AVAILABLE_MOORDYNPLUS when this feature is compiled.
#ifdef DISABLE_MOORDYNPLUS
  #define AVAILABLE_MOORDYNPLUS false
#else
  #define AVAILABLE_MOORDYNPLUS true
#endif

//-Defines AVAILABLE_GPU when this feature is compiled.
#ifdef _WITHGPU
  #define AVAILABLE_GPU true
#else
  #define AVAILABLE_GPU false
#endif

//-Defines AVAILABLE_MGPU when this feature is compiled.
#ifdef _WITHMGPU
  #define AVAILABLE_MGPU true
#else
  #define AVAILABLE_MGPU false
#endif


#define DELTA_HEAVYFLOATING  ///<Applies DDT to fluid particles interacting with floatings with higher density (massp>MassFluid*1.2). | Aplica DDT a fluido que interaccionan con floatings pesados (massp>MassFluid*1.2). NO_COMENTARIO

#define CELLDIV_OVERMEMORYNP 0.05f   ///<Memory that is reserved for the particle management in JCellDivGpu. | Memoria que se reserva de mas para la gestion de particulas en JCellDivGpu.
#define CELLDIV_OVERMEMORYCELLS 2    ///<Number of cells in each dimension is increased to allocate memory for JCellDivGpu cells. | Numero celdas que se incrementa en cada dimension al reservar memoria para celdas en JCellDivGpu.
#define CELLDIV_OVERMEMORYNCELLS 6553598 ///<Number of cells (6.5M for 100 MiB on GPU) to increase memory allocation for JCellDivGpu cells.
#define PERIODIC_OVERMEMORYNP 0.05f  ///<Memory reserved for the creation of periodic particles in JSphGpuSingle::RunPeriodic(). | Mermoria que se reserva de mas para la creacion de particulas periodicas en JSphGpuSingle::RunPeriodic().
#define PARTICLES_OVERMEMORY_MIN 128 ///<Minimum over memory allocated on CPU or GPU according number of particles.

#define BORDER_MAP 0.05

#define ALMOSTZERO 1e-18f

#define BSIZE_FORCES 128  ///<Blocksize for particle interaction (default=128).

#define AVAILABLE_FLEXSTRUC true //<vs_flexstruc>
#define MAX_NUM_MKCLAMP 8 ///<Maximum number of mkclamps for each flexible structure. <vs_flexstruc>

// #define AVAILABLE_DIVCLEAN  //<vs_divclean>
//#define CODE_SIZE4  //-Enables or disables the use of unsigned type (32 bits) for code (allows valid 65530 MKs). | Activa o desactiva el uso de unsigned (32 bits) para code (permite 65530 MKs validos).
#ifdef CODE_SIZE4
  #define CODE_MKRANGEMAX 65530        //-Maximum valid MK value. | Valor maximo de MK valido.
  typedef unsigned typecode;           //-Type of the variable code using 4 bytes.
  //-Code of the particles:
  #define CODE_MASKSPECIAL 0xe0000000  //-Bits de special: 11100000 00000000 00000000 00000000              | Special bits: 11100000 00000000 00000000 00000000
  #define CODE_MASKSPECIALMV 29        //-Displacemente in bits of special bits.
  #define CODE_NORMAL   0x0            //-0  Particulas normales no excluidas.                              | 0 Normal particles (not excluded)
  #define CODE_PERIODIC 0x20000000     //-1  Particulas duplicadas por periodicas.                          | 1 Duplicate particles for periodic
  #define CODE_OUTIGNORE 0x40000000    //-2  Marca particulas que se van a ignorar en el siguiente divide.  | 2 Brands particles to be ignored in the next division.
  #define CODE_OUTMOV 0x60000000       //-3  Particulas normales excluidas por movimiento.                  | 3 Normal particles excluded for motion
  #define CODE_OUTPOS 0x80000000       //-4  Particulas normales excluidas por posicion.                    | 4 Normal particles excluded for position
  #define CODE_OUTRHO 0xA0000000       //-5  Particulas normales excluidas por densidad.                    | 5 Normal particles excluded for density

  #define CODE_MASKTYPEVALUE 0x0003ffff //-Bits for type: 00000000 00000011 11111111 11111111
  #define CODE_MASKTYPE 0x00030000      //-Bits for type: 00000000 00000011 00000000 00000000
  #define CODE_TYPE_FIXED 0x0           //---Particles fixed:       0- 65535
  #define CODE_TYPE_MOVING 0x00010000   //---Particles moving:  65536-131071
  #define CODE_TYPE_FLOATING 0x00020000 //---Particles float:  131072-196607
  #define CODE_TYPE_FLUID 0x00030000    //---Particles fluid:  196608-262143
  #define CODE_MASKVALUE 0x00000ffff    //-Bits type-value: 0000 0111 1111 1111  Range:0-65535

  #define CODE_TYPE_FLUID_LIMITFREE 0x0003ffdf //---Last normal fluid code: 262111
  #define CODE_TYPE_FLUID_INOUT     0x0003ffe0 //---First inlet/outlet code: 262112 (16 different codes for InOut zones + 16 to select input particles).
  #define CODE_TYPE_FLUID_INOUTNUM  16      //---Maximum number of valid inlet/outlet zones.
  #define CODE_TYPE_FLUID_INOUTMASK 31      //---Mask to obtain zone value.
  #define CODE_TYPE_FLUID_INOUT015MASK 15   //---Mask to obtain zone value (ignore extra bit).
#else
  #define CODE_MKRANGEMAX 2047      //-Maximum valid MK value. | Valor maximo de MK valido.
  typedef word typecode;            //-Type of the variable code using 2 bytes.
  //-Code of the particles:
  #define CODE_MASKSPECIAL 0xe000   //-Bits de special:     1110 0000 0000 0000                          | Special bits: 1110 0000 0000 0000
  #define CODE_MASKSPECIALMV 13     //-Displacemente in bits of special bits.
  #define CODE_NORMAL   0x0         //-0  Particulas normales no excluidas.                              | 0 Normal particles (not excluded)
  #define CODE_PERIODIC 0x2000      //-1  Particulas duplicadas por periodicas.                          | 1 Duplicate particles for periodic
  #define CODE_OUTIGNORE 0x4000     //-2  Marca particulas que se van a ignorar en el siguiente divide.  | 2 Brands particles to be ignored in the next division.
  #define CODE_OUTMOV 0x6000        //-3  Particulas normales excluidas por movimiento.                  | 3 Normal particles excluded for motion
  #define CODE_OUTPOS 0x8000        //-4  Particulas normales excluidas por posicion.                    | 4 Normal particles excluded for position
  #define CODE_OUTRHO 0xA000        //-5  Particulas normales excluidas por densidad.                    | 5 Normal particles excluded for density
  //#define CODE_SPECIAL1 0xC000    //-6  Por ejemplo, CODE_DOMAINPREV pertenece a Proceso-1             | 6 For example, CODE_DOMAINPREV belongs to Process-1
  //#define CODE_SPECIAL2 0xE000    //-7  Por ejemplo, CODE_DOMAINNEXT pertenece a Proceso+1             | 7 For example, CODE_DOMAINNEXT belongs to Process+1

  #define CODE_MASKTYPEVALUE 0x1fff //-Bits for type: 0001 1111 1111 1111
  #define CODE_MASKTYPE 0x1800      //-Bits for type: 0001 1000 0000 0000
  #define CODE_TYPE_FIXED 0x0       //---Particles fixed:  0-2047                                         
  #define CODE_TYPE_MOVING 0x800    //---Particles moving: 2048-4095                                      
  #define CODE_TYPE_FLOATING 0x1000 //---Particles float:  4096-6143                                      
  #define CODE_TYPE_FLUID 0x1800    //---Particles fluid:  6144-8191                                      
  #define CODE_MASKVALUE 0x7ff      //-Bits type-value: 0000 0111 1111 1111  Range:0-2047

  //<vs_flexstruc_ini>
  #define CODE_TYPE_FLEXSTRUC_MASK 0x07e0         //-Bits for flexible structure: 0000 0111 1110 0000
  #define CODE_TYPE_FLEXSTRUCCLAMP_MASK 0x0010    //-Bits for flexible structure clamp: 0000 0000 0001 0000
  //<vs_flexstruc_end>

  #define CODE_TYPE_FLUID_LIMITFREE 0x1fdf  //---Last normal fluid code: 8159
  #define CODE_TYPE_FLUID_INOUT     0x1fe0  //---First inlet/outlet code: 8160 (16 different codes for InOut zones + 16 to select input particles).
  #define CODE_TYPE_FLUID_INOUTNUM  16      //---Maximum number of valid inlet/outlet zones.
  #define CODE_TYPE_FLUID_INOUTMASK 31      //---Mask to obtain zone value.
  #define CODE_TYPE_FLUID_INOUT015MASK 15   //---Mask to obtain zone value (ignore extra bit).
#endif

#define CODE_SetNormal(code)    (code&(~CODE_MASKSPECIAL))
#define CODE_SetPeriodic(code)  (CODE_SetNormal(code)|CODE_PERIODIC)
#define CODE_SetOutIgnore(code) (CODE_SetNormal(code)|CODE_OUTIGNORE)
#define CODE_SetOutPos(code)    (CODE_SetNormal(code)|CODE_OUTPOS)
#define CODE_SetOutMov(code)    (CODE_SetNormal(code)|CODE_OUTMOV)
#define CODE_SetOutRho(code)    (CODE_SetNormal(code)|CODE_OUTRHO)
#define CODE_GetSpecialValue(code) (code&CODE_MASKSPECIAL)
#define CODE_GetSpecialByte(code) (CODE_GetSpecialValue(code)>>CODE_MASKSPECIALMV)

#define CODE_GetType(code) (code&CODE_MASKTYPE)
#define CODE_GetTypeValue(code) (code&CODE_MASKVALUE)
#define CODE_GetTypeAndValue(code) (code&CODE_MASKTYPEVALUE)

#define CODE_IsNormal(code)    (CODE_GetSpecialValue(code)==CODE_NORMAL)
#define CODE_IsPeriodic(code)  (CODE_GetSpecialValue(code)==CODE_PERIODIC)
#define CODE_IsNotOut(code)    (CODE_GetSpecialValue(code)<=CODE_PERIODIC)
#define CODE_IsOutRho(code)   (CODE_GetSpecialValue(code)==CODE_OUTRHO)
#define CODE_IsOutIgnore(code) (CODE_GetSpecialValue(code)==CODE_OUTIGNORE)

#define CODE_IsFixed(code)    (CODE_GetType(code)==CODE_TYPE_FIXED)
#define CODE_IsMoving(code)   (CODE_GetType(code)==CODE_TYPE_MOVING)
#define CODE_IsFloating(code) (CODE_GetType(code)==CODE_TYPE_FLOATING)
#define CODE_IsFluid(code)    (CODE_GetType(code)==CODE_TYPE_FLUID)
#define CODE_IsNotFluid(code) (CODE_GetType(code)!=CODE_TYPE_FLUID)

//<vs_flexstruc_ini>
#define CODE_IsFlexStrucAny(code)   ((CODE_IsFixed(code)||CODE_IsMoving(code)) && (code&CODE_TYPE_FLEXSTRUC_MASK)==CODE_TYPE_FLEXSTRUC_MASK)
#define CODE_IsFlexStrucFlex(code)  (CODE_IsMoving(code) && (code&(CODE_TYPE_FLEXSTRUC_MASK|CODE_TYPE_FLEXSTRUCCLAMP_MASK))==CODE_TYPE_FLEXSTRUC_MASK)
#define CODE_IsFlexStrucClamp(code) ((CODE_IsFixed(code)||CODE_IsMoving(code)) && (code&(CODE_TYPE_FLEXSTRUC_MASK|CODE_TYPE_FLEXSTRUCCLAMP_MASK))==(CODE_TYPE_FLEXSTRUC_MASK|CODE_TYPE_FLEXSTRUCCLAMP_MASK))

#define CODE_ToFlexStrucFlex(code,ibody)  ((code&(CODE_MASKSPECIAL|CODE_MASKTYPE))|CODE_TYPE_FLEXSTRUC_MASK|ibody)
#define CODE_ToFlexStrucClamp(code,ibody) ((code&(CODE_MASKSPECIAL|CODE_MASKTYPE))|CODE_TYPE_FLEXSTRUC_MASK|CODE_TYPE_FLEXSTRUCCLAMP_MASK|ibody)
#define CODE_GetIbodyFlexStruc(code)      (code&(~(CODE_MASKSPECIAL|CODE_MASKTYPE|CODE_TYPE_FLEXSTRUC_MASK|CODE_TYPE_FLEXSTRUCCLAMP_MASK)))
//<vs_flexstruc_end>

//#define CODE_IsFluidInout(code)    (CODE_IsFluid(code) && CODE_GetTypeAndValue(code)>=CODE_TYPE_FLUID_INOUT)
#define CODE_IsFluidInout(code)    (CODE_GetTypeAndValue(code)>=CODE_TYPE_FLUID_INOUT)
#define CODE_IsFluidNotInout(code) (CODE_IsFluid(code) && CODE_GetTypeAndValue(code)< CODE_TYPE_FLUID_INOUT)

#define CODE_ToFluidInout(code,izone) (code&(~CODE_MASKTYPEVALUE))|(CODE_TYPE_FLUID_INOUT|izone)
#define CODE_GetIzoneFluidInout(code) (code&CODE_TYPE_FLUID_INOUTMASK)

//<vs_vrres_ini>
  #define CODE_TYPE_FLUID_BUFFER 	0x1fc0    //First buffer code:8128
  #define CODE_TYPE_FLUID_BUFFERNUM 	16
  #define CODE_TYPE_FLUID_BUFFERMASK 	31
  #define CODE_TYPE_FLUID_BUFFER015MASK 15
    #define CODE_TYPE_FLUID_FIXED 0x1f80


  #define CODE_IsFluidBuffer(code)    (CODE_GetTypeAndValue(code)>=CODE_TYPE_FLUID_BUFFER && CODE_GetTypeAndValue(code)<CODE_TYPE_FLUID_INOUT)
  #define CODE_IsFluidNotBuffer(code) (CODE_IsFluid(code) && CODE_GetTypeAndValue(code)< CODE_TYPE_FLUID_BUFFER)
  #define CODE_ToFluidBuffer(code,izone) (code&(~CODE_MASKTYPEVALUE))|(CODE_TYPE_FLUID_BUFFER|izone)
  #define CODE_GetIzoneFluidBuffer(code) (code&CODE_TYPE_FLUID_BUFFERMASK)
  #define CODE_IsFluidFixed(code)    (CODE_GetTypeAndValue(code)>=CODE_TYPE_FLUID_FIXED && CODE_GetTypeAndValue(code)<CODE_TYPE_FLUID_BUFFER)
  #define CODE_ToFluidFixed(code,izone) (code&(~CODE_MASKTYPEVALUE))|(CODE_TYPE_FLUID_FIXED|izone)
  #define CODE_GetIzoneFluidFixed(code) (code&CODE_TYPE_FLUID_BUFFERMASK)
//<vs_vrres_end>

//<vs_flexstruc_ini>
///Structure with the information of the flexible structure.
typedef struct{
  unsigned nc;                          ///<Number of clamping objects.
  typecode clampcode[MAX_NUM_MKCLAMP];  ///<Code for clamping particles.
  float vol0;                           ///<Initial particle volume.
  float rho0;                           ///<Initial particle density.
  float youngmod;                       ///<Young's modulus.
  float poissratio;                     ///<Poisson ratio.
  float hgfactor;                       ///<Hourglass correction factor.
  tmatrix6f cmat;                       ///<Constitutive matrix.
}StFlexStrucData;
//<vs_flexstruc_end>

///Structure with the information of the floating object.
typedef struct{
  word mkbound;     ///<MkBound of the floating object.
  unsigned begin;   ///<First particle of the floating object.
  unsigned count;   ///<Number of floating objects.
  float mass;       ///<Mass of the floating object (units:Kg).
  float massp;      ///<Mass of the particle of the floating object (units:Kg).
  float radius;     ///<Maximum distance between particles and center (units:m).
  byte constraints; ///<Translation and rotation restrictions (combination of TpFtConstrains values).
  tdouble3 center;  ///<Center of the floating object (units:m).
  tfloat3 angles;   ///<Rotation angles from center (angle xz, angle yz, angle xy) (units:Rad).
  tfloat3 fvel;     ///<Linear velocity of the floating object (units:m/s).
  tfloat3 fomega;   ///<Angular velocity of the floating object (units:rad/s).
  tfloat3 facelin;  ///<Linear acceleration of the floating object computed from velocity difference (units:m/s^2).
  tfloat3 faceang;  ///<Angular acceleration of the floating object computed from velocity difference (units:rad/s^2).
  tmatrix3f inertiaini; ///<Initial state inertia tensor in world coordinates (computed or user-given).
  bool usechrono;   ///<Activates the use of Chrono library.
  //-Additional data stored for output results only.
  tfloat3 extforcelin;  ///<External linear forces (moorings and imposed forces) (units:N).
  tfloat3 extforceang;  ///<External angular forces (moorings and imposed forces) (units:N*m*rad).
  tfloat3 fluforcelin;  ///<Linear forces from fluid (sum in eq.48 at Dominguez et al 2022) (units:N).
  tfloat3 fluforceang;  ///<Angular forces from fluid (sum in eq.49 at Dominguez et al 2022) (units:N*m*rad).
  tfloat3 preacelin;    ///<Linear acceleration before constraints (includes external forces and gravity) (units:m/s^2).
  tfloat3 preaceang;    ///<Angular acceleration before constraints (multiplied by rotated inertia tensor) (units:rad/s^2).
}StFloatingData;

///Structure with the information of the solid object for DEM interaction (Discrete Element Method).
typedef struct{ //(DEM)
  float mass;         ///<Mass of the object (units:Kg).
  float massp;        ///<Mass of the particle of the floating object (units:Kg).
  float young;        ///<Young Modulus of the floating object (units:N/m2).
  float poisson;      ///<Poisson coefficient of the floating object (units:-).
  float kfric;        ///<Kinetic friction coefficient of the floating object (units:-).
  float sfric;        ///<Static friction coefficient of the floating object (units:-).
  float tau;          ///<Value of (1-poisson^2)/young (units:-).
  float restitu;      ///<Restitution Coefficient (units:-).
}StDemData;

///Structure that saves extra information about the execution.
typedef struct StrInfoPartPlus{
  unsigned nct;        ///<Number of cells used in the divide.
  unsigned nctsize;    ///<Number of supported cells with memory allocated.
  unsigned npsim;      ///<Number of particles used in simulation (normal + periodic particles).
  unsigned npsize;     ///<Number of supported particles with memory allocated.
  unsigned npnormal;   ///<Number of normal particles used in simulation (without periodic particles).
  unsigned npsave;     ///<Number of selected particles to save.
  unsigned npnew;      ///<Number of new fluid particles (inlet conditions).
  unsigned noutpos;    ///<Number of excluded particles due to invalid position.
  unsigned noutrho;    ///<Number of excluded particles due to invalid density.
  unsigned noutmov;    ///<Number of excluded particles due to invalid movement.
  unsigned npbin;      ///<Number of boundary particles within the area of the divide (includes periodic particles).
  unsigned npbout;     ///<Number of boundary particles outside of the area of the divide (includes periodic particles).
  unsigned npf;        ///<Number of floating+fluid particles (includes periodic particles).
  unsigned npbper;     ///<Number of periodic boundary particles (inside and outside the area of the split).
  unsigned npfper;     ///<Number of periodic floating+fluid particles.
  llong memorycpualloc;
  llong memorynpalloc;
  llong memorynpused;
  llong memorynctalloc;
  llong memorynctused;
  double timesim;      ///<Seconds from the start of the simulation (after loading the initial data).                    | Segundos desde el inicio de la simulacion (despues de cargar los datos iniciales).
  bool gpudata;

  StrInfoPartPlus(){ 
    nct=nctsize=0;
    npsim=npsize=npnormal=npsave=0;
    npnew=0;
    noutpos=noutrho=noutmov=0;
    npbin=npbout=npf=npbper=npfper=0;
    memorycpualloc=memorynpalloc=memorynpused=0;
    memorynctalloc=memorynctused=0;
    timesim=0;
    gpudata=false;
  }
  /// Stores basic values displayed in JSph::SaveData().
  void SetBasic(unsigned npsim,unsigned npnormal,unsigned npf,unsigned nct){
    this->npsim=npsim;
    this->npnormal=npnormal;
    this->npf=npf;
    this->nct=nct;
  }
  void SetNct(unsigned nct,unsigned nctsize){
    this->nct=nct; this->nctsize=nctsize;
  }
  void SetNp(unsigned npsim,unsigned npsize,unsigned npnormal,unsigned npsave){
    this->npsim=npsim;       this->npsize=npsize; 
    this->npnormal=npnormal; this->npsave=npsave;
  }
  void SetNout(unsigned noutpos,unsigned noutrho,unsigned noutmov){
    this->noutpos=noutpos; this->noutrho=noutrho; this->noutmov=noutmov;
  }
  void SetNpExtra(unsigned npbin,unsigned npbout,unsigned npf
    ,unsigned npbper,unsigned npfper)
  {
    this->npbin=npbin;   this->npbout=npbout;  this->npf=npf;
    this->npbper=npbper; this->npfper=npfper;
  }
  //-Adding values functions.
  void AddBasic(const StrInfoPartPlus& v){
    this->npsim+=v.npsim;
    this->npnormal+=v.npnormal;
    this->npf+=v.npf;
    this->nct+=v.nct;
  }
  void AddNct(const StrInfoPartPlus& v){
    this->nct+=v.nct; this->nctsize+=v.nctsize;
  }
  void AddNp(const StrInfoPartPlus& v){
    this->npsim+=v.npsim;       this->npsize+=v.npsize; 
    this->npnormal+=v.npnormal; this->npsave+=v.npsave;
    this->npnew+=v.npnew;
  }
  void AddNpExtra(const StrInfoPartPlus& v)
  {
    this->npbin+=v.npbin;   this->npbout+=v.npbout;  this->npf+=v.npf;
    this->npbper+=v.npbper; this->npfper+=v.npfper;
  }
  void AddMemory(const StrInfoPartPlus& v)
  {
    memorycpualloc+=v.memorycpualloc;
    memorynpalloc +=v.memorynpalloc;   memorynpused +=v.memorynpused;
    memorynctalloc+=v.memorynctalloc;  memorynctused+=v.memorynctused;
  }
}StInfoPartPlus;

///Structure that stores the maximum values (or almost) achieved during the simulation.
typedef struct StrMaxNumbers{
  llong memcpu;       ///<Amount of reserved CPU memory.
  llong memgpu;       ///<Amount of reserved GPU memory.
  llong memgpunct;    ///<Amount of reserved GPU memory for cells in memgpu.
  unsigned particles; ///<Maximum number of particles in simulation.
  unsigned cells;     ///<Maximum number of cells in simulation.
  StrMaxNumbers(){ Clear(); }
  StrMaxNumbers(llong memcpu,llong memgpu,llong memgpunct
    ,unsigned particles,unsigned cells)
  {
    this->memcpu=memcpu; 
    this->memgpu=memgpu; 
    this->memgpunct=memgpunct;
    this->particles=particles; 
    this->cells=cells;
  }
  void Clear(){ 
    memcpu=memgpu=memgpunct=0;
    particles=cells=0;
  }
}StMaxNumbers;

///Structure with values to apply external acceleration to fluid.
typedef struct{
  unsigned codesel1;  ///<First code for application (It is disabled with UINT_MAX).
  unsigned codesel2;  ///<Last code for application.
  tdouble3 acclin;
  tdouble3 accang;
  tdouble3 centre;
  tdouble3 velang;
  tdouble3 vellin;
  bool setgravity;
}StAceInput;

///Controls the output of information on the screen and/or log.
typedef enum{ 
  MOUT_ScrFile=3,  ///<Output on the screen and log.
  MOUT_File=2,     ///<Output in the log.
  MOUT_Screen=1,   ///<Output on the screen.
  MOUT_None=0      ///<No output.
}TpModeOut;   

///Data output options.
typedef enum{ 
  SDAT_Binx=1,       ///<BYNARY format .bi2
  SDAT_Vtk=2,        ///<VTK format .vtk
  SDAT_Csv=4,        ///<CSV format .csv
  SDAT_Info=8,
  SDAT_None=0 
}TpSaveDat; 

///Types of step algorithm.
typedef enum{ 
  STEP_Symplectic=2,  ///<Symplectic algorithm.
  STEP_Verlet=1,      ///<Verlet algorithm.
  STEP_None=0 
}TpStep;                    

///Types of kernel function.
typedef enum{ 
  KERNEL_Wendland=2,   ///<Wendland kernel.
  KERNEL_Cubic=1,      ///<Cubic Spline kernel.
  KERNEL_None=0 
}TpKernel;                  

///Types of viscosity treatment.
//-Defined in JViscosityDef.h

///Types of boundary conditions.
typedef enum{ 
  BC_DBC=1    ///<Dynamic Boundary Condition (DBC).
 ,BC_MDBC=2   ///<mDBC.
}TpBoundary;

///Slip modes for mDBC. 
typedef enum{ 
  SLIP_None=0      ///<mDBC is not used.
 ,SLIP_Vel0=1      ///<mDBC original: DBC vel=0
 ,SLIP_NoSlip=2    ///<mDBC2 slip mode: No-slip
 ,SLIP_FreeSlip=3  ///<mDBC2 slip mode: Free slip (in development).
}TpSlipMode;

typedef enum{
  MDBC2_None=0    ///<mDBC original
 ,MDBC2_Std=1     ///mDBC2 
 ,MDBC2_NoPen=2   ///mDBC2 No Penetration.
}TpMdbc2Mode;

#define MDBC2_KEEPVEL  ///<En MdbcBoundCorrection() no modifica velrho para usarlo en interaccion en lugar de velmotion.

#define BMODE_DBC 0       ///<Boundary particle use DBC.
#define BMODE_MDBC2 1     ///<Boundary particle use mDBC2.
#define BMODE_MDBC2OFF 2  ///<Boundary particle use mDBC2 but mass is disabled.

///Types of interaction step.
typedef enum{ 
  INTERSTEP_None=0,         
  INTERSTEP_Verlet=1,       ///<Interaction to compute forces using the Verlet algorithm.
  INTERSTEP_SymPredictor=2, ///<Interaction to compute forces using the Symplectic algorithm (predictor step). 
  INTERSTEP_SymCorrector=3  ///<Interaction to compute forces using the Symplectic algorithm (corrector step). 
}TpInterStep;

///Types of density diffussion term.
typedef enum{ 
  DDT_None=0 
 ,DDT_DDT=1      ///<Density Diffussion Term 1 (Molteni and Colagrossi 2009). It is only applied to inner fluid particles.
 ,DDT_DDT2=2     ///<Density Diffussion Term 2 (Fourtakas et al 2019). It is only applied to inner fluid particles.
 ,DDT_DDT2Full=3 ///<Density Diffussion Term 2 (Fourtakas et al 2019). It is applied to all fluid particles.
}TpDensity;

///Types of Shifting applied to fluid particles. 
typedef enum{
  SHIFT_None=0     ///<Shifting is not applied.
 ,SHIFT_NoBound=1  ///<Shifting is applied to fluid particles except those that interact with all boundaries.
 ,SHIFT_NoFixed=2  ///<Shifting is applied to fluid particles except those that interact with fixed boundaries.
 ,SHIFT_Full=3     ///<Shifting is applied to all fluid particles.
 ,SHIFT_FS=4       ///<Advanced shifting for free-surface detection. //<vs_advshift>
}TpShifting; 

///Structure with main SPH constants and configurations.
typedef struct{
  bool simulate2d;          ///<Toggles 2D simulation (cancels forces in Y axis).
  double simulate2dposy;    ///<Y value in 2D simulations.

  TpKernel tkernel;               ///<Kernel type: Cubic or Wendland.
  fsph::StKCubicCte      kcubic;  ///<Constants for the Cubic Spline kernel.
  fsph::StKWendlandCte   kwend;   ///<Constants for the Wendland kernel.

  float kernelh;            ///<The smoothing length of SPH kernel (h) [m].
  float cteb;               ///<Constant used in the state equation [Pa].
  float gamma;              ///<Politropic constant for water used in the state equation.
  float rhopzero;           ///<Reference density of the fluid [kg/m3].
  double dp;                ///<Initial distance between particles [m].
  float massfluid;          ///<Reference mass of the fluid particle [kg].
  float massbound;          ///<Reference mass of the general boundary particle [kg].
  tfloat3 gravity;          ///<Gravitational acceleration [m/s^2].

  //-Constants for computation (computed starting from previous constants).
  float kernelsize;         ///<Maximum interaction distance between particles (2h) (KernelK*KernelH).
  float kernelsize2;        ///<Maximum interaction distance squared (2h*2h) (KernelSize^2).
  double cs0;               ///<Speed of sound at the reference density.
  float eta2;               ///<Constant related to H (Eta2=(h*0.1)*(h*0.1)).
  //-Other constants for computation of special features.
  float spssmag;            ///<Smagorinsky constant used in SPS turbulence model.
  float spsblin;            ///<Blin constant used in the SPS turbulence model.
  float ddtkhcte;           ///<Store fixed constant DDTkh.
  float ddtkh;              ///<Constant for DDT1 & DDT2. DDTkh=DDTValue*KernelSize
  float ddtgz;              ///<Constant for DDT2.        DDTgz=RhopZero*Gravity.z/CteB
}StCteSph;

///Returns empty StCteSph structure.
inline StCteSph CteSphNull(){
  StCteSph c={false,0,KERNEL_None
    ,{0,0,0,0,0,0,0,0}
    ,{0,0}
    ,0,0,0,0,0,0,0,{0,0,0},0,0,0,0,0,0,0,0,0};
  return(c);
}

///Rigid Algorithm for floating collisions.
typedef enum{ 
   FTRIGID_Free=0    ///<Collision-free               (FtMode=FTMODE_Ext, UseDEM=false ). 
  ,FTRIGID_Sph=1     ///<Collision in terms of SPH    (FtMode=FTMODE_Sph, UseDEM=false ). 
  ,FTRIGID_Dem=2     ///<Collision in terms of DEM    (FtMode=FTMODE_Ext, UseDEM=true  ). 
  ,FTRIGID_Chrono=3  ///<Collision in terms of Chrono (FtMode=FTMODE_Ext, UseDEM=true, UseChrono=true). 
}TpRigidMode;  

///Returns the name of the RigidMode in text format.
inline const char* GetNameRigidMode(TpRigidMode rigidmode){
  switch(rigidmode){
    case FTRIGID_Free:   return("Collision-Free");
    case FTRIGID_Sph:    return("SPH");
    case FTRIGID_Dem:    return("DCDEM");
    case FTRIGID_Chrono: return("CHRONO");
  }
  return("???");
}

///Interaction mode for floatings.
typedef enum{ 
   FTMODE_None=0            ///<There are not floatings.
  ,FTMODE_Sph=1             ///<Interaction between floatings and boundaries in terms of SPH.
  ,FTMODE_Ext=2             ///<Interaction between floatings and boundaries in terms of DEM or CHRONO.
}TpFtMode;  

#define USE_FLOATING (ftmode!=FTMODE_None)
#define USE_NOFLOATING (ftmode==FTMODE_None)
#define USE_FTEXTERNAL (ftmode==FTMODE_Ext)


///Mask values for translation or rotation constraints applied to floating bodies.
typedef enum{ 
  FTCON_Free=0,     ///<No translation or rotation constraints.
  FTCON_MoveX=1,    ///<Translation in X is avoided.
  FTCON_MoveY=2,    ///<Translation in Y is avoided.
  FTCON_MoveZ=4,    ///<Translation in Z is avoided.
  FTCON_RotateX=8,  ///<Rotation in X is avoided.
  FTCON_RotateY=16, ///<Rotation in Y is avoided.
  FTCON_RotateZ=32  ///<Rotation in Z is avoided.
}TpFtConstrains;

///Returns combination of TpFtConstrains values to define the constraints.
inline byte ComputeConstraintsValue(const tint3& translationfree,const tint3& rotationfree){
  return((translationfree.x? 0: FTCON_MoveX)
        +(translationfree.y? 0: FTCON_MoveY)
        +(translationfree.z? 0: FTCON_MoveZ)
        +(rotationfree.x   ? 0: FTCON_RotateX)
        +(rotationfree.y   ? 0: FTCON_RotateY)
        +(rotationfree.z   ? 0: FTCON_RotateZ));
}

///Applies constraints.
inline void ApplyConstraints(byte constraints,tfloat3& linear,tfloat3& angular){
  if(constraints&FTCON_MoveX  )linear.x=0;
  if(constraints&FTCON_MoveY  )linear.y=0;
  if(constraints&FTCON_MoveZ  )linear.z=0;
  if(constraints&FTCON_RotateX)angular.x=0;
  if(constraints&FTCON_RotateY)angular.y=0;
  if(constraints&FTCON_RotateZ)angular.z=0;
}

//<vs_flexstruc_ini>
///Constitutive model for flexible structures.
typedef enum{
  CONSTITMODEL_None=1,          ///<None.
  CONSTITMODEL_PlaneStrain=1,   ///<Plane strain.
  CONSTITMODEL_PlaneStress=2,   ///<Plane stress.
  CONSTITMODEL_SVK=3            ///<St. Venant Kirchhoff.
}TpConstitModel;
//<vs_flexstruc_end>


///Modes of cells division.
typedef enum{ 
   CELLMODE_None=0
  ,CELLMODE_Full=1    ///<Cells of size KernelSize (maximum interaction distance).
  ,CELLMODE_Half=2    ///<Cells of size KernelSize/2.
}TpCellMode; 

///Returns the name of the CellMode in text format.
inline const char* GetNameCellMode(TpCellMode cellmode){
  switch(cellmode){
    case CELLMODE_Full:   return("Full");
    case CELLMODE_Half:   return("Half");
  }
  return("???");
}


///Domain division mode.
typedef enum{ 
  MGDIV_None=0      ///<Not specified. 
 ,MGDIV_X=1         ///<Main division in X direction.
 ,MGDIV_Y=2         ///<Main division in Y direction.
 ,MGDIV_Z=3         ///<Main division in Z direction (used for Single-GPU).
}TpMgDivMode;  

///Returns the name of division mode in text.
inline const char* GetNameDivision(TpMgDivMode axis){
  switch(axis){
    case MGDIV_None:  return("None");
    case MGDIV_X:     return("X");
    case MGDIV_Y:     return("Y");
    case MGDIV_Z:     return("Z");
  }
  return("???");
}


//##############################################################################
// Codification of global cells (of size PosCellSize) according to the position for particle interaction using pos-cell method on GPU.
//##############################################################################
//#define PSCEL_CONFIG_USER //-Configuration defined by user for special simulation cases.
#ifdef PSCEL_CONFIG_USER
  // Place reserved for configuration defined by user for special simulation cases.
  // Place reserved for configuration defined by user for special simulation cases.
  // Place reserved for configuration defined by user for special simulation cases.
#else
  #define PSCEL_CONFIG_13_10_9 //-Typical configuration of cells in X, Y and Z for 3-D simulations.
  //#define PSCEL_CONFIG_17_2_13 //-Configuration to maximize cells in X and Z for 2-D simulations.

  //-Cell configuration 13_10_9:
  #ifdef PSCEL_CONFIG_13_10_9
    #define PSCEL1_X 0xfff80000  //-Mask of bits for cell X: 13 bits for 8,192 ->  1111 1111 1111 1000   0000 0000 0000 0000
    #define PSCEL1_Y 0x0007fe00  //-Mask of bits for cell Y: 10 bits for 1,024 ->  0000 0000 0000 0111   1111 1110 0000 0000
    #define PSCEL1_Z 0x000001ff  //-Mask of bits for cell Z:  9 bits for   512 ->  0000 0000 0000 0000   0000 0001 1111 1111
    #define PSCEL1_MOVX 19       //-Displacement to obaint X cell.
    #define PSCEL1_MOVY 9        //-Displacement to obaint Y cell.
  #endif
  //-Cell configuration 17_2_13:
  #ifdef PSCEL_CONFIG_17_2_13
    #define PSCEL1_X 0xffff8000  //-Mask of bits for cell X: 17 bits for 131,072 ->  1111 1111 1111 1111   1000 0000 0000 0000
    #define PSCEL1_Y 0x00006000  //-Mask of bits for cell Y:  2 bits for       4 ->  0000 0000 0000 0000   0110 0000 0000 0000
    #define PSCEL1_Z 0x00001fff  //-Mask of bits for cell Z: 13 bits for   8,192 ->  0000 0000 0000 0000   0001 1111 1111 1111
    #define PSCEL1_MOVX 15       //-Displacement to obaint X cell.
    #define PSCEL1_MOVY 13       //-Displacement to obaint Y cell.
  #endif
#endif

//-Cell configuration in use:
#define PSCEL_X PSCEL1_X        //-Selected mask of bits for cell X
#define PSCEL_Y PSCEL1_Y        //-Selected mask of bits for cell Y
#define PSCEL_Z PSCEL1_Z        //-Selected mask of bits for cell Z
#define PSCEL_MOVX PSCEL1_MOVX  //-Selected displacement to obaint X cell.
#define PSCEL_MOVY PSCEL1_MOVY  //-Selected displacement to obaint Y cell.
//-Cell methods:
#define PSCEL_GetPartX(cel)  (cel&PSCEL_X)                          ///<Returns bits of X cell.
#define PSCEL_GetPartY(cel)  (cel&PSCEL_Y)                          ///<Returns bits of Y cell.
#define PSCEL_GetPartZ(cel)  (cel&PSCEL_Z)                          ///<Returns bits of Z cell.
#define PSCEL_GetX(cel)      (cel>>PSCEL_MOVX)                      ///<Returns X coordinate of the cell.
#define PSCEL_GetY(cel)      (PSCEL_GetPartY(cel)>>PSCEL_MOVY)      ///<Returns Y coordinate of the cell.
#define PSCEL_GetZ(cel)      (cel&PSCEL_Z)                          ///<Returns Z coordinate of the cell.
#define PSCEL_Code(cx,cy,cz) ((cx<<PSCEL_MOVX)|(cy<<PSCEL_MOVY)|cz) ///<Returns code for cell (cx,cy,cz).

//-Returns cell coordinates from float number (for GPU).
#define PSCEL_GetfX(celf) (int(PSCEL_GetX(__float_as_uint(celf)))) ///<Returns X coordinate of the cell (for GPU).
#define PSCEL_GetfY(celf) (int(PSCEL_GetY(__float_as_uint(celf)))) ///<Returns Y coordinate of the cell (for GPU).
#define PSCEL_GetfZ(celf) (int(PSCEL_GetZ(__float_as_uint(celf)))) ///<Returns Z coordinate of the cell (for GPU).


#endif


