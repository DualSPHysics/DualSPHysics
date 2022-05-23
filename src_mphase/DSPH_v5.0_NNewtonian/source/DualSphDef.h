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

/// \file DualSphDef.h \brief Includes definitions and specific types for DualSPHysics program.

#ifndef _DualSphDef_
#define _DualSphDef_

#include "TypesDef.h"
#include "JParticlesDef.h"
#include "JPeriodicDef.h"
#include "FunSphKernelDef.h"
#include "OmpDefs.h"
#include <algorithm>

//#define DISABLE_KERNELS_EXTRA  ///<Compiles Wendland kernel and ignores the rest of the SPH kernels.

#define DISABLE_MDBC_EXTRAMODES  ///<Compiles only slipmode=SLIP_Vel0 and ignores the rest since they are not ready.

//#define DISABLE_TIMERS     ///<Compiles without timers. | Compilado sin timers.
//#define DISABLE_BSMODES    ///<compiles without advanced BlockSize modes.


//-Removes dependencies from precompiled libraries.
#include "JNumexLibDef.h"    //Defines DISABLE_NUMEXLIB to compile without Numex library.
#include "JVtkLibDef.h"      //Defines DISABLE_VTKLIB to compile without VTK library.
//#define DISABLE_CHRONO     ///<It allows compile without ChronoLib library (dsphchrono.dll, ChronoEngine.dll and ChronoEngine_parallel.dll).
#define DISABLE_CHRONO_OMP   ///<It allows compile without parallel module of Chrono (ignores ChronoEngine_parallel.dll).
//#define DISABLE_WAVEGEN    ///<It allows compile without Wave-Paddles, Multi-Layer Pistons and Relaxation Zones libraries.
//#define DISABLE_MOORDYN    ///<It allows compile without LibDSphMoorDyn library.


//-Defines AVAILABLE_VTKLIB when this feature is compiled.
#ifdef DISABLE_VTKLIB
  #define AVAILABLE_VTKLIB false
#else
  #define AVAILABLE_VTKLIB true
#endif

//-Defines AVAILABLE_NUMEXLIB when this feature is compiled.
#ifdef DISABLE_NUMEXLIB
  #define AVAILABLE_NUMEXLIB false
#else
  #define AVAILABLE_NUMEXLIB true
#endif

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

//-Defines AVAILABLE_MOORDYN when this feature is compiled.
#ifdef DISABLE_MOORDYN
  #define AVAILABLE_MOORDYN false
#else
  #define AVAILABLE_MOORDYN true
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
#define CELLDIV_OVERMEMORYCELLS 1    ///<Number of cells in each dimension is increased to allocate memory for JCellDivGpu cells. | Numero celdas que se incrementa en cada dimension al reservar memoria para celdas en JCellDivGpu.
#define PERIODIC_OVERMEMORYNP 0.05f  ///<Memory reserved for the creation of periodic particles in JSphGpuSingle::RunPeriodic(). | Mermoria que se reserva de mas para la creacion de particulas periodicas en JSphGpuSingle::RunPeriodic().
#define PARTICLES_OVERMEMORY_MIN 128 ///<Minimum over memory allocated on CPU or GPU according number of particles.

#define BORDER_MAP 0.05

#define ALMOSTZERO 1e-18f

#define BSIZE_FORCES 128  ///<Blocksize for particle interaction (default=128).

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
  #define CODE_OUTMOVE 0x60000000      //-3  Particulas normales excluidas por movimiento.                  | 3 Normal particles excluded for motion
  #define CODE_OUTPOS 0x80000000       //-4  Particulas normales excluidas por posicion.                    | 4 Normal particles excluded for position
  #define CODE_OUTRHOP 0xA0000000      //-5  Particulas normales excluidas por densidad.                    | 5 Normal particles excluded for density

  #define CODE_MASKTYPEVALUE 0x0003ffff //-Bits for type: 00000000 00000011 11111111 11111111
  #define CODE_MASKTYPE 0x00030000      //-Bits for type: 00000000 00000011 00000000 00000000
  #define CODE_TYPE_FIXED 0x0           //---Particles fixed:       0- 65535
  #define CODE_TYPE_MOVING 0x00010000   //---Particles moving:  65536-131071
  #define CODE_TYPE_FLOATING 0x00020000 //---Particles float:  131072-196607
  #define CODE_TYPE_FLUID 0x00030000    //---Particles fluid:  196608-262143
  #define CODE_MASKVALUE 0x00000ffff    //-Bits type-value: 0000 0111 1111 1111  Range:0-65535

  #define CODE_TYPE_FLUID_LIMITFREE 0x0003ffdf //---Last normal fluid code: 262111
  #define CODE_TYPE_FLUID_INOUT     0x0003ffe0 //---First inlet/outlet code: 262112 (16 different codes for InOut zones + 16 to select input particles).
  #define CODE_TYPE_FLUID_INOUTNUM  16      //---Maximum number of vali d inlet/outlet zones.
  #define CODE_TYPE_FLUID_INOUTMASK 31      //---Mask to obtain zone value.
#else
  #define CODE_MKRANGEMAX 2047      //-Maximum valid MK value. | Valor maximo de MK valido.
  typedef word typecode;            //-Type of the variable code using 2 bytes.
  //-Code of the particles:
  #define CODE_MASKSPECIAL 0xe000   //-Bits de special:     1110 0000 0000 0000                          | Special bits: 1110 0000 0000 0000
  #define CODE_MASKSPECIALMV 13     //-Displacemente in bits of special bits.
  #define CODE_NORMAL   0x0         //-0  Particulas normales no excluidas.                              | 0 Normal particles (not excluded)
  #define CODE_PERIODIC 0x2000      //-1  Particulas duplicadas por periodicas.                          | 1 Duplicate particles for periodic
  #define CODE_OUTIGNORE 0x4000     //-2  Marca particulas que se van a ignorar en el siguiente divide.  | 2 Brands particles to be ignored in the next division.
  #define CODE_OUTMOVE 0x6000       //-3  Particulas normales excluidas por movimiento.                  | 3 Normal particles excluded for motion
  #define CODE_OUTPOS 0x8000        //-4  Particulas normales excluidas por posicion.                    | 4 Normal particles excluded for position
  #define CODE_OUTRHOP 0xA000       //-5  Particulas normales excluidas por densidad.                    | 5 Normal particles excluded for density
  //#define CODE_SPECIAL1 0xC000    //-6  Por ejemplo, CODE_DOMAINPREV pertenece a Proceso-1             | 6 For example, CODE_DOMAINPREV belongs to Process-1
  //#define CODE_SPECIAL2 0xE000    //-7  Por ejemplo, CODE_DOMAINNEXT pertenece a Proceso+1             | 7 For example, CODE_DOMAINNEXT belongs to Process+1

  #define CODE_MASKTYPEVALUE 0x1fff //-Bits for type: 0001 1111 1111 1111
  #define CODE_MASKTYPE 0x1800      //-Bits for type: 0001 1000 0000 0000
  #define CODE_TYPE_FIXED 0x0       //---Particles fixed:  0-2047                                         
  #define CODE_TYPE_MOVING 0x800    //---Particles moving: 2048-4095                                      
  #define CODE_TYPE_FLOATING 0x1000 //---Particles float:  4096-6143                                      
  #define CODE_TYPE_FLUID 0x1800    //---Particles fluid:  6144-8191                                      
  #define CODE_MASKVALUE 0x7ff      //-Bits type-value: 0000 0111 1111 1111  Range:0-2047

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
#define CODE_SetOutMove(code)   (CODE_SetNormal(code)|CODE_OUTMOVE)
#define CODE_SetOutRhop(code)   (CODE_SetNormal(code)|CODE_OUTRHOP)
#define CODE_GetSpecialValue(code) (code&CODE_MASKSPECIAL)
#define CODE_GetSpecialByte(code) (CODE_GetSpecialValue(code)>>CODE_MASKSPECIALMV)

#define CODE_GetType(code) (code&CODE_MASKTYPE)
#define CODE_GetTypeValue(code) (code&CODE_MASKVALUE)
#define CODE_GetTypeAndValue(code) (code&CODE_MASKTYPEVALUE)

#define CODE_IsNormal(code)    (CODE_GetSpecialValue(code)==CODE_NORMAL)
#define CODE_IsPeriodic(code)  (CODE_GetSpecialValue(code)==CODE_PERIODIC)
#define CODE_IsNotOut(code)    (CODE_GetSpecialValue(code)<=CODE_PERIODIC)
#define CODE_IsOutRhop(code)   (CODE_GetSpecialValue(code)==CODE_OUTRHOP)
#define CODE_IsOutIgnore(code) (CODE_GetSpecialValue(code)==CODE_OUTIGNORE)

#define CODE_IsFixed(code)    (CODE_GetType(code)==CODE_TYPE_FIXED)
#define CODE_IsMoving(code)   (CODE_GetType(code)==CODE_TYPE_MOVING)
#define CODE_IsFloating(code) (CODE_GetType(code)==CODE_TYPE_FLOATING)
#define CODE_IsFluid(code)    (CODE_GetType(code)==CODE_TYPE_FLUID)
#define CODE_IsNotFluid(code) (CODE_GetType(code)!=CODE_TYPE_FLUID)

//#define CODE_IsFluidInout(code)    (CODE_IsFluid(code) && CODE_GetTypeAndValue(code)>=CODE_TYPE_FLUID_INOUT)
#define CODE_IsFluidInout(code)    (CODE_GetTypeAndValue(code)>=CODE_TYPE_FLUID_INOUT)
#define CODE_IsFluidNotInout(code) (CODE_IsFluid(code) && CODE_GetTypeAndValue(code)< CODE_TYPE_FLUID_INOUT)

#define CODE_ToFluidInout(code,izone) (code&(~CODE_MASKTYPEVALUE))|(CODE_TYPE_FLUID_INOUT|izone)
#define CODE_GetIzoneFluidInout(code) (code&CODE_TYPE_FLUID_INOUTMASK)


///Defines type of movement.
typedef enum{ 
  MOTT_None=0,    ///<No movement.
  MOTT_Linear=1,  ///<Linear movement.
  MOTT_Matrix=2   ///<Matrix movement (for rotations).
}TpMotionType;   

///Structure with the information for moving particles (lineal and matrix movement).
typedef struct{
  word ref;            ///<Idx of moving object.
  word mkbound;        ///<MkBound of moving particles.
  unsigned idbegin;    ///<First id of moving particles.
  unsigned count;      ///<Number of moving particles.
  TpMotionType type;   ///<Type of motion (none, linear, matrix).
  tdouble3 linmov;     ///<Linear displacement to apply to the particles position.
  tdouble3 linvel;     ///<Linear velocity for particles.
  tdouble3 linace;     ///<Linear acceleration for particles (when acceleration movement is computed).
  tmatrix4d matmov;    ///<Matrix transformation to apply to the particles position.
  tmatrix4d matmov2;   ///Matrix transformation to compute acceleration of particles (when acceleration movement is computed).
}StMotionData;

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
}StFloatingData;

///Structure with the information of the floating object in forces calculation.
typedef struct{
  tfloat3 face;       ///<Sum of particle acceleration (units:m/s2). | Sumatorio de ace de particulas.
  tfloat3 fomegaace;  ///<Angular acceleration of the floating object (units:rad/s2). | Aceleracion angular del objecto floating.
}StFtoForces;

///Structure with the information of the floating object in forces calculation.
typedef struct{
  tfloat3 fomegares;   ///<Calculated angular velocity to upadte floating body (units:rad/s).
  tfloat3 fvelres;     ///<Calculated linear velocity to upadte floating body (units:m/s).
  tdouble3 fcenterres; ///<Calculated center to upadte floating body (units:m).
}StFtoForcesRes;

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

///Structure that stores the maximum values (or almost) achieved during the simulation.
typedef struct StrMaxNumbers{
  llong memcpu;       ///<Amount of reserved CPU memory. | Cantidad de memoria Cpu reservada.            
  llong memgpu;       ///<Amount of reserved GPU memory. | Cantidad de memoria Gpu reservada.
  unsigned particles; ///<Maximum number of particles.   | Numero maximo de particulas.
  unsigned cells;     ///<Maximum number of cells.       | Numero maximo de celdas.                   
  StrMaxNumbers(){ Clear(); }
  StrMaxNumbers(llong vmemcpu,llong vmemgpu,unsigned vparticles,unsigned vcells){
    memcpu=vmemcpu; memgpu=vmemgpu; particles=vparticles; cells=vcells;
  }
  void Clear(){ 
    memcpu=memgpu=0; particles=cells=0;
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


//<vs_non-Newtonian_ini>
///Structure that holds const. eq. multiphase constants
typedef struct {  
    word mkfluid;				///<Mk of the phase
    unsigned idbegin;		///<First id of phase.
    unsigned count;			///<Number of particles in phase.
    float visco;				///<viscosity of the phase.    
    float tau_max;			///<maximum tau of the phase.
	  float Bi_multi;			///<viscosity multiplier for bi-visocity model of the phase.
    float m_NN;					///<HBP model n parameter
    float n_NN;					///<HBP model m parameter
    float coh;					///<Not use in this current version
    float phi;					///<Not use in this current version
    float DP_alpha;			///<Not use in this current version
    float DP_kappa;			///<Not use in this current version
    float tau_yield;		///<Yield strength of phase
    unsigned phasetype; ///<Typer of phase
}StPhaseCte;

///Structure PhaseArray holds physical proprerties of phases
typedef struct {
    int phaseid;		///<ID of the phase
    float rho;			///<Density of the phase
    float mass;			///<Mass of the phase
    float Cs0;			///<Speed of sound of phase
    float CteB;			///<B coefficient in EOS of the phase
    float Gamma;		///<Polytropic index of the phase
}StPhaseArray;
//<vs_non-Newtonian_end>


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

//<vs_non-Newtonian_ini>
///Types of velocity gradient calculation.
typedef enum {
  VELGRAD_None = 0,           
  VELGRAD_FDA = 1,      ///<Velocity gradient by Shao etal 2003, eq 20.
  VELGRAD_SPH = 2       ///<Velocity gradient by Fourtakas etal 2016.
}TpVelGrad;
//<vs_non-Newtonian_end>

///Types of viscosity treatment.
typedef enum{ 
  VISCO_ConstEq = 3,				 ///<Constitutive equation for viscous forces.  //<vs_non-Newtonian
  VISCO_LaminarSPS=2,        ///<Laminar viscosity and Sub-Partice Scale Turbulence.
  VISCO_Artificial=1,        ///<Artificial viscosity.
  VISCO_None=0 
}TpVisco;            

///Types of boundary conditions.
typedef enum{ 
  BC_MDBC=2,   ///<M-DBC.
  BC_DBC=1     ///<Dynamic Boundary Condition (DBC).
}TpBoundary;

///Types of boundary conditions. 
typedef enum{ 
  SLIP_FreeSlip=3,  ///<Free slip
  SLIP_NoSlip=2,    ///<No-slip
  SLIP_Vel0=1       ///<DBC vel=0
}TpSlipMode;

///Types of interaction step.
typedef enum{ 
  INTERSTEP_None=0,         
  INTERSTEP_Verlet=1,       ///<Interaction to compute forces using the Verlet algorithm.
  INTERSTEP_SymPredictor=2, ///<Interaction to compute forces using the Symplectic algorithm (predictor step). 
  INTERSTEP_SymCorrector=3  ///<Interaction to compute forces using the Symplectic algorithm (corrector step). 
}TpInterStep;

///Types of density diffussion term.
typedef enum{ 
  DDT_DDT2Full=3, ///<Density Diffussion Term 2 (Fourtakas et al 2019). It is applied to all fluid particles.
  DDT_DDT2=2,     ///<Density Diffussion Term 2 (Fourtakas et al 2019). It is only applied to inner fluid particles.
  DDT_DDT=1,      ///<Density Diffussion Term. It is only applied to inner fluid particles.
  DDT_None=0 
}TpDensity;

///Types of Shifting applied to fluid particles. 
typedef enum{
  SHIFT_Full=3,             ///<Shifting is applied to all fluid particles.
  SHIFT_NoFixed=2,          ///<Shifting is applied to fluid particles except those that interact with fixed boundaries.
  SHIFT_NoBound=1,          ///<Shifting is applied to fluid particles except those that interact with all boundaries.
  SHIFT_None=0              ///<Shifting is not applied.
}TpShifting; 



///Structure with main SPH constants and configurations.
typedef struct{
  bool simulate2d;          ///<Toggles 2D simulation (cancels forces in Y axis).
  double simulate2dposy;    ///<Y value in 2D simulations.

  TpKernel tkernel;               ///<Kernel type: Cubic or Wendland.
  fsph::StKCubicCte      kcubic;  ///<Constants for the Cubic Spline kernel.
  fsph::StKWendlandCte   kwend;   ///<Constants for the Wendland kernel.

  float kernelh;            ///<The smoothing length of SPH kernel [m].
  float cteb;               ///<Constant used in the state equation [Pa].
  float gamma;              ///<Politropic constant for water used in the state equation.
  float rhopzero;           ///<Reference density of the fluid [kg/m3].
  double dp;                ///<Initial distance between particles [m].
  float massfluid;          ///<Reference mass of the fluid particle [kg].
  float massbound;          ///<Reference mass of the general boundary particle [kg].
  tfloat3 gravity;          ///<Gravitational acceleration [m/s^2].

  //-Constants for computation (computed starting from previous constants).
  float kernelsize;         ///<Maximum interaction distance between particles (KernelK*KernelH).
  float kernelsize2;        ///<Maximum interaction distance squared (KernelSize^2).
  double cs0;               ///<Speed of sound at the reference density.
  float eta2;               ///<Constant related to H (Eta2=(h*0.1)*(h*0.1)).
  float spssmag;            ///<Smagorinsky constant used in SPS turbulence model.
  float spsblin;            ///<Blin constant used in the SPS turbulence model.
  float ddtkh;              ///<Constant for DDT1 & DDT2. DDTkh=DDTValue*KernelSize
  float ddtgz;              ///<Constant for DDT2.        DDTgz=RhopZero*Gravity.z/CteB
}StCteSph;

///Returns empty StCteSph structure.
inline StCteSph CteSphNull(){
  StCteSph c={false,0,KERNEL_None
    ,{0,0,0,0,0,0,0,0}
    ,{0,0}
    ,0,0,0,0,0,0,0,{0,0,0},0,0,0,0,0,0,0,0};
  return(c);
}

///Interaction mode for floatings.
typedef enum{ 
  FTMODE_None=0,            ///<There are not floatings.
  FTMODE_Sph=1,             ///<Interaction between floatings and boundaries in terms of SPH.
  FTMODE_Ext=2              ///<Interaction between floatings and boundaries in terms of DEM or CHRONO.
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
inline byte ComputeConstraintsValue(const tint3 &translationfree,const tint3 &rotationfree){
  return((translationfree.x? 0: FTCON_MoveX)
        +(translationfree.y? 0: FTCON_MoveY)
        +(translationfree.z? 0: FTCON_MoveZ)
        +(rotationfree.x   ? 0: FTCON_RotateX)
        +(rotationfree.y   ? 0: FTCON_RotateY)
        +(rotationfree.z   ? 0: FTCON_RotateZ));
}

///Applies constraints.
inline void ApplyConstraints(byte constraints,tfloat3 &linear,tfloat3 &angular){
  if(constraints&FTCON_MoveX  )linear.x=0;
  if(constraints&FTCON_MoveY  )linear.y=0;
  if(constraints&FTCON_MoveZ  )linear.z=0;
  if(constraints&FTCON_RotateX)angular.x=0;
  if(constraints&FTCON_RotateY)angular.y=0;
  if(constraints&FTCON_RotateZ)angular.z=0;
}


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
  MGDIV_None=0,      ///<Not specified. 
  MGDIV_X=1,         ///<Main division in X direction.
  MGDIV_Y=2,         ///<Main division in Y direction.
  MGDIV_Z=3          ///<Main division in Z direction.
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


///Codification of local cells according to the position.
#define PC__CodeMapOut   0xffffffff
#define PC__CodeSpecial  0x80000000
#define PC__CodeDomLeft  0x80000000
#define PC__CodeDomRight 0x80000001
#define PC__GetCode(sx,sy,sz) (((sx+1)<<25)|(sy<<20)|(sz<<15)|((sy+sz)<<10)|((sx+1+sz)<<5)|(sx+1+sy))  //-Clave de codificacion (orden de valores: sx,sy,sz,sy+sz,sx+sz,sx+sy). | Encryption key (order of values: sx,sy,sz,sy+sz,sz+sx,sx+sy).
#define PC__GetSx(cc) (cc>>25)       //-Numero de bits para coordenada X de celda. | Number of bits for X coordinate cell.
#define PC__GetSy(cc) ((cc>>20)&31)  //-Numero de bits para coordenada Y de celda. | Number of bits for Y coordinate cell.
#define PC__GetSz(cc) ((cc>>15)&31)  //-Numero de bits para coordenada Z de celda. | Number of bits for Z coordinate cell.
#define PC__Cellx(cc,cel) ((*((unsigned*)&cel))>>((cc>>10)&31))             //-Coordenada X de celda. | X coordinate of the cell.
#define PC__Celly(cc,cel) (((*((unsigned*)&cel))<<(cc>>25))>>((cc>>5)&31))  //-Coordenada Y de celda. | Y coordinate of the cell.
#define PC__Cellz(cc,cel) (((*((unsigned*)&cel))<<(cc&31))>>(cc&31))        //-Coordenada Z de celda. | Z coordinate of the cell.
#define PC__Cell(cc,cx,cy,cz) ((cx<<((cc>>10)&31))|(cy<<((cc>>15)&31))|cz)  //-Valor de celda para cx, cy y cz. | Cell value for cx,cy and cz.
#define PC__MaxCellx(cc) ((0xffffffff>>((cc>>10)&31))>>1)           //-Coordenada X de celda maxima. | Maximum X coordinate of the cell.
#define PC__MaxCelly(cc) ((0xffffffff<<(cc>>25))>>((cc>>5)&31))     //-Coordenada Y de celda maxima. | Maximum Y coordinate of the cell.
#define PC__MaxCellz(cc) ((0xffffffff<<(cc&31))>>(cc&31))           //-Coordenada Z de celda maxima. | Maximum Z coordinate of the cell.


///Codification of global cells (of size PosCellSize) according to the position for particle interaction using pos-cell method.
//#define CEL_CONFIG_USER //-Configuration defined by user for special simulation cases.
#ifdef CEL_CONFIG_USER
  // Place reserved for configuration defined by user for special simulation cases.
  // Place reserved for configuration defined by user for special simulation cases.
  // Place reserved for configuration defined by user for special simulation cases.
#else
  #define CEL_CONFIG_13_10_9 //-Typical configuration of cells in X, Y and Z for 3-D simulations.
  //#define CEL_CONFIG_17_2_13 //-Configuration to maximize cells in X and Z for 2-D simulations.

  //-Cell configuration 13_10_9:
  #ifdef CEL_CONFIG_13_10_9
    #define CEL1_X 0xfff80000  //-Mask of bits for cell X: 13 bits for 8,192 ->  1111 1111 1111 1000   0000 0000 0000 0000
    #define CEL1_Y 0x0007fe00  //-Mask of bits for cell Y: 10 bits for 1,024 ->  0000 0000 0000 0111   1111 1110 0000 0000
    #define CEL1_Z 0x000001ff  //-Mask of bits for cell Z:  9 bits for   512 ->  0000 0000 0000 0000   0000 0001 1111 1111
    #define CEL1_MOVX 19       //-Displacement to obaint X cell.
    #define CEL1_MOVY 9        //-Displacement to obaint Y cell.
  #endif
  //-Cell configuration 17_2_13:
  #ifdef CEL_CONFIG_17_2_13
    #define CEL1_X 0xffff8000  //-Mask of bits for cell X: 17 bits for 131,072 ->  1111 1111 1111 1111   1000 0000 0000 0000
    #define CEL1_Y 0x00006000  //-Mask of bits for cell Y:  2 bits for       4 ->  0000 0000 0000 0000   0110 0000 0000 0000
    #define CEL1_Z 0x00001fff  //-Mask of bits for cell Z: 13 bits for   8,192 ->  0000 0000 0000 0000   0001 1111 1111 1111
    #define CEL1_MOVX 15       //-Displacement to obaint X cell.
    #define CEL1_MOVY 13       //-Displacement to obaint Y cell.
  #endif
#endif

//-Cell configuration in use:
#define CEL_X CEL1_X        //-Selected mask of bits for cell X
#define CEL_Y CEL1_Y        //-Selected mask of bits for cell Y
#define CEL_Z CEL1_Z        //-Selected mask of bits for cell Z
#define CEL_MOVX CEL1_MOVX  //-Selected displacement to obaint X cell.
#define CEL_MOVY CEL1_MOVY  //-Selected displacement to obaint Y cell.
//-Cell methods:
#define CEL_GetPartX(cel)  (cel&CEL_X)                        //-Returns bits of X cell.
#define CEL_GetPartY(cel)  (cel&CEL_Y)                        //-Returns bits of Y cell.
#define CEL_GetPartZ(cel)  (cel&CEL_Z)                        //-Returns bits of Z cell.
#define CEL_GetX(cel)      (cel>>CEL_MOVX)                    //-Returns X coordinate of the cell.
#define CEL_GetY(cel)      (CEL_GetPartY(cel)>>CEL_MOVY)      //-Returns Y coordinate of the cell.
#define CEL_GetZ(cel)      (cel&CEL_Z)                        //-Returns Z coordinate of the cell.
#define CEL_Code(cx,cy,cz) ((cx<<CEL_MOVX)|(cy<<CEL_MOVY)|cz) //-Returns code for cell (cx,cy,cz).


#endif


