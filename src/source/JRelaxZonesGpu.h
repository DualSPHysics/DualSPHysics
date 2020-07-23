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

//#############################################################################
//# Cambios:
//# =========
//# - Clase para uso de GPU por parte de Relaxation Zone. (28-05-2018)
//#############################################################################

/// \file JRelaxZonesGpu.h \brief Declares the class \ref JRelaxZonesGpu.

#ifndef _JRelaxZonesGpu_
#define _JRelaxZonesGpu_

#include "TypesDef.h"
#include "JObject.h"


//##############################################################################
//# JRelaxZoneRegularGpu
//##############################################################################
/// \brief Manages the GPU computations for regular waves using RZ (JRelaxZoneRegular class).

class JRelaxZoneRegularGpu : protected JObject
{
public:
  JRelaxZoneRegularGpu();
  inline llong GetAllocMemoryGpu()const{ return(0); }

  void SetFluidVel(unsigned n,unsigned pini,bool order2,bool subdrift
    ,double centerx,float widthhalf,float coeffx,float coeffz
    ,double falpha,double fbeta,double fsub,double fdiv
    ,double timewave,double swl,double kl,double sinhkld
    ,double wpf,double cta,double depth,double framp
    ,double ct2,double sinhkld4
    ,double ctd,double ctd2,unsigned fluidbeginidp
    ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop);

};


//##############################################################################
//# JRelaxZonesSpectrumGpu
//##############################################################################
/// \brief Manages the GPU computations for irregular waves using RZ (JRelaxZoneSpectrum class).

class JRelaxZoneSpectrumGpu : protected JObject
{
private:
  bool WavesOnGpu;
  double *WaveKlg;     //-TWOPI / Wave length [WavesCount].
  double *WaveAmpg;    //-Amplitud [WavesCount].
  double *WaveFangg;   //-Frecuencia angular [WavesCount]. (WaveFang=Freq*TWOPI)
  double *WavePhaseg;  //-Initial phase [WavesCount].

  llong MemGpuFixed;  
  void AllocMemoryGpu(unsigned wavecount);

public:
  JRelaxZoneSpectrumGpu();
  ~JRelaxZoneSpectrumGpu();
  void FreeMemoryGpu();
  inline bool GetWavesOnGpu()const{ return(WavesOnGpu); }

  void PrepareWaveDataGpu(unsigned wavecount
    ,const double *kl,const double *amp,const double *fang,const double *phase);

  inline llong GetAllocMemoryGpu()const{ return(MemGpuFixed); }

  void SetFluidVelSpectrumSub(unsigned n,unsigned pini
    ,double centerx,float widthhalf,float coeffx,float coeffz
    ,double falpha,double fbeta,double fsub,double fdiv
    ,double timewave,double swl,double depth,double framp,unsigned wavecount
    ,unsigned fluidbeginidp
    ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop
    ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2);

};


//##############################################################################
//# JRelaxZonesExternalGpu
//##############################################################################
/// \brief Manages the GPU computations for external velocity using RZ (JRelaxZoneExternal class).

class JRelaxZonesExternalGpu : protected JObject
{
private:
  bool GpuReady;
  double *VelXg;
  double *VelZg;

  llong MemGpuFixed;  
  void AllocMemoryGpu(unsigned size,bool loadvelz);

public:
  JRelaxZonesExternalGpu();
  ~JRelaxZonesExternalGpu();
  void FreeMemoryGpu();

  void PrepareDataGpu(unsigned size,bool loadvelz,const double *velx,const double *velz);

  inline llong GetAllocMemoryGpu()const{ return(MemGpuFixed); }

  void SetFluidVelExternal(unsigned n,unsigned pini
    ,double centerx,float widthhalf,float coeffx,float coeffz
    ,double falpha,double fbeta,double fsub,double fdiv
    ,double pxmin,double pymin,double pzmin
    ,double dpx,double dpy,double dpz
    ,unsigned npx1,unsigned npy1,unsigned npz1
    ,unsigned fluidbeginidp
    ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop
    ,bool subdrift,double fun,double ctd,double ctd2,double ctd_2,double ctd2_2,double bottom);

};


//##############################################################################
//# JRelaxZoneUniformGpu
//##############################################################################
/// \brief Manages the GPU computations for RZ uniform application (JRelaxZoneUniform class).

class JRelaxZoneUniformGpu : protected JObject
{
public:
  JRelaxZoneUniformGpu();
  inline llong GetAllocMemoryGpu()const{ return(0); }

  void SetFluidVelUniform(unsigned n,unsigned pini
    ,const tfloat3 &vt,const tfloat4 &cenpla
    ,const tfloat4 &dompla1,const tfloat4 &dompla2,const tfloat4 &dompla3
    ,const float domsize1,const float domsize2,const float domsize3,float widthhalf
    ,float coeff,double falpha,double fbeta,double fsub,double fdiv,unsigned fluidbeginidp
    ,const tdouble2 *posxy,const double *posz,const unsigned *idp,tfloat4 *velrhop);

};


#endif

