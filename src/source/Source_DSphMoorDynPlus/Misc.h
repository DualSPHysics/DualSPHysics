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

/// \file Misc.h \brief Defines the class \ref Misc.

#ifndef _Misc_
#define _Misc_

#ifndef LINUX
#define LINUX
#endif // !LINUX

#include "TypesDef.h"

#include <iostream>
#include <vector>
#include <string>
#include <stdlib.h>
#include <cmath>
#include <complex>

#include <fstream>
#include <sstream>
#include <cstring>

#include <memory>

 //#ifdef USEGL
 // #include <GL/gl.h>  // for openGL drawing option
 // #include <GL/glu.h> // used in arrow function
 //#endif

#include "kiss_fft.h"  // used for any wave kinematics functions

#ifdef OSX
#include <sys/uio.h>
#elif defined LINUX

#else
#include <windows.h>  // these are for guicon function RedirectIOToConsole
#include <io.h>
#endif

#include <stdio.h>
#include <fcntl.h>
#include <iostream>
#include <fstream>

// note: this file contains the struct definitions for environmental and line/connect properties


// from Iñaki Zabala
namespace misc{

#ifdef _MSC_VER
template<typename T> static inline T round(T val) { return floor(val + 0.5); }
#endif

typedef std::complex<double> doubleC; ///<Make shorthand for complex double type
typedef std::complex<float>  floatC;  ///<Make shorthand for complex float type


const doubleC i1(0.,1.); ///<Set imaginary number 1
const floatC i1f(0.,1.); ///<Set imaginary number 1
constexpr int wordy=0; ///<Flag to enable excessive output if is greater than 0 for troubleshooting

// below are function prototypes for misc functions

double eye(int I,int J); /// simple convenience function for identity matrix
void unitvector(std::vector<double> & u,std::vector<double> & r1,std::vector<double> & r2); /// Returns unit vector (u) in direction from r1 to r2
void unitvector(tdouble3& u,tdouble3& r1,tdouble3& r2); /// Returns unit vector (u) in direction from r1 to r2
void inverse3by3(std::vector<std::vector<double>> & minv,std::vector<std::vector<double>> & m); /// Computes the inverse of a matrix m
void RotMat(double x1,double x2,double x3,tmatrix3d& transMat); /// Creates rotation matrix  (row major order?)
double dotprod(std::vector<double>& A,std::vector<double>& B); /// Computes dot product and returns it
double dotprod(double A[],std::vector<double>& B); /// Computes dot product and returns it
// calculate wave number from frequency,g,and depth (credit: FAST source)
float WaveNumber(float Omega,float g,float h); /// Calculates wave number from frequency,g,and depth (credit: FAST source)
float JONSWAP(float Omega,float Hs,float Tp,float Gamma); /// Compute the JONSWAP wave spectrum
float SINHNumOvrSIHNDen(float k,float h,float z); /// Computes the hyperbolic numerator over denominator
float COSHNumOvrSIHNDen(float k,float h,float z); /// Computes the hyperbolic numerator over denominator
void reverse(double* data,int datasize); /// Flip or reverse function
void doIIR(double* in,double* out,int dataSize,double* a,double* b,int kernelSize); /// IIR filter function
void doSSfilter(double* in,double* out,int dataSize,double* a,double* beta,double b0,int kernelSize); /// State Space filter function - a good reference is https://ccrma.stanford.edu/~jos/fp/Converting_State_Space_Form_Hand.html
}
#endif //!Misc
