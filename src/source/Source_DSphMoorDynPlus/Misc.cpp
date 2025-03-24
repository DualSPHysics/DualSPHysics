/*===================================================================================
<MOORDYNPLUS> Copyright (c) 2025 by Ivan Martinez-Estevez

Ivan Martinez Estevez and Jose M. Dominguez (Universidade de Vigo, Spain)
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

#include "Misc.h"
#include "TypesMoorDynPlus.h"

namespace misc{
//======================================================================
/// simple convenience function for identity matrix
//======================================================================
double eye(int I,int J) {
  if(I == J) { return 1.0; }
  else { return 0.0; }
}

//======================================================================
/// Returns unit vector (u) in direction from r1 to r2
//======================================================================
void unitvector(std::vector<double> &u,std::vector<double> & r1,std::vector<double> & r2) {
  double length_squared=0.0;

  for(int J=0; J<3; J++) { length_squared+=(r2[J]-r1[J])*(r2[J]-r1[J]); }

  double length=sqrt(length_squared);

  for(int J=0; J<3; J++) { u[J]=(r2[J]-r1[J])/length; }// write to unit vector

  return;
}

//======================================================================
/// Returns unit vector (u) in direction from r1 to r2
//======================================================================
void unitvector(tdouble3& u,tdouble3& r1,tdouble3& r2) {
  const tdouble3 length_squared=(r2-r1)*(r2-r1);
  double length=sqrt(length_squared.x+length_squared.y+length_squared.z);
  u=(r2-r1)/length; // write to unit vector
  return;
}

//======================================================================
/// Computes the inverse of a matrix m
//======================================================================
void inverse3by3(std::vector<std::vector<double>> & minv,std::vector<std::vector<double>> & m) {
  // computes the inverse of a matrix m
  double det=m[0][0]*(m[1][1]*m[2][2]-m[2][1]*m[1][2]) -
    m[0][1]*(m[1][0]*m[2][2]-m[1][2]*m[2][0]) +
    m[0][2]*(m[1][0]*m[2][1]-m[1][1]*m[2][0]);

  double invdet=1/det;
  //cout << endl << "inverted matrix: ";
  minv[0][0]=(m[1][1]*m[2][2]-m[2][1]*m[1][2])*invdet; //cout << minv[0][0];
  minv[0][1]=(m[0][2]*m[2][1]-m[0][1]*m[2][2])*invdet; //cout << minv[0][1];
  minv[0][2]=(m[0][1]*m[1][2]-m[0][2]*m[1][1])*invdet; //cout << minv[0][2];
  minv[1][0]=(m[1][2]*m[2][0]-m[1][0]*m[2][2])*invdet; //cout << minv[1][0];
  minv[1][1]=(m[0][0]*m[2][2]-m[0][2]*m[2][0])*invdet; //cout << minv[1][1];
  minv[1][2]=(m[1][0]*m[0][2]-m[0][0]*m[1][2])*invdet; //cout << minv[1][2];
  minv[2][0]=(m[1][0]*m[2][1]-m[2][0]*m[1][1])*invdet; //cout << minv[2][0];
  minv[2][1]=(m[2][0]*m[0][1]-m[0][0]*m[2][1])*invdet; //cout << minv[2][1];
  minv[2][2]=(m[0][0]*m[1][1]-m[1][0]*m[0][1])*invdet; //cout << minv[2][2] << endl;

}

//======================================================================
/// Creates rotation matrix  (row major order?)
//======================================================================
void RotMat(double x2,double x1,double x3,tmatrix3d& transMat) {
  // note above swapping of x1 and x2 to deal with weird coordinate system from FAST convention
  // ( x2 is roll,x1 is pitch,x3 is yaw )

  const float s1=float(sin(x1));
  const float c1=float(cos(x1));
  const float s2=float(sin(x2));
  const float c2=float(cos(x2));
  const float s3=float(sin(x3));
  const float c3=float(cos(x3));

  //rmat=transpose([ 1 0 0; 0 c1 -s1; 0 s1 c1]*[c2 0 s2; 0 1 0; -s2 0 c2]*[c3 s3 0; -s3 c3 0; 0 0 1]);
  //rmat=[1 rz -ry; -rz 1 rx; ry -rx 1];   % rotation matrix

  transMat.a11=c1*c3 + s1*s2*s3;
  transMat.a12=c3*s1*s2-c1*s3;
  transMat.a13=c2*s1;
  transMat.a21=c2*s3;
  transMat.a22=c2*c3;
  transMat.a23=-s2;
  transMat.a31=c1*s2*s3-c3*s1;
  transMat.a32=s1*s3 + c1*c3*s2;
  transMat.a33=c1*c2;
}

//======================================================================
/// Computes dot product and returns it
//======================================================================
double dotprod(std::vector<double>& A,std::vector<double>& B) {
  double ans=0.;
  for(int i=0; i<static_cast<int>(A.size()); i++) ans+=A[i]*B[i];
  return ans;
}

//======================================================================
/// Computes dot product and returns it
//======================================================================
double dotprod(double A[],std::vector<double>& B) {
  double ans=0.;
  for(int i=0; i<static_cast<int>(B.size()); i++) ans+=A[i]*B[i];
  return ans;
}

//======================================================================
/// Calculate wave number from frequency,g,and depth (credit: FAST source)
//======================================================================
float WaveNumber(float Omega,float g,float h) {
  // 
  // This FUNCTION solves the finite depth dispersion relationship:
  // 
  //                   k*tanh(k*h)=(Omega^2)/g
  // 
  // for k, the wavenumber (WaveNumber) given the frequency, Omega,
  // gravitational constant, g, and water depth, h, as inputs.  A
  // high order initial guess is used in conjunction with a quadratic
  // Newton's method for the solution with seven significant digits
  // accuracy using only one iteration pass.  The method is due to
  // Professor J.N. Newman of M.I.T. as found in routine EIGVAL of
  // the SWIM-MOTION-LINES (SML) software package in source file
  // Solve.f of the SWIM module.
  // 
  // Compute the wavenumber, unless Omega is zero, in which case, return
  //   zero:
  // 
  float k,X0;

  if(Omega == 0.0) {   // When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.
    k=0.0;
    return k;
  }
  else {     // Omega > 0.0 solve for the wavenumber as usual.
    float C=Omega*Omega*h/g;
    float CC=C*C;

    // Find X0:
    if(C <= 2.0) { X0=sqrt(C)*(1.0f + C*(0.169f + (0.031f*C))); }
    else {
      float E2=float(exp(-2.0*C));
      X0=C*(1.0f + (E2*(2.0f-(12.0f*E2))));
    }

    // Find the WaveNumber:

    if(C <= 4.8) {
      float C2=CC-X0*X0;
      float A=1.0f/(C-C2);
      float B=A*((0.5f *log((X0 + C)/(X0-C)))-X0);

      k=(X0-(B*C2*(1.0f + (A*B*C*X0))))/h;
    }
    else { k=X0/h; }

    if(Omega<0) { k=-k; }  // @mth: modified to return negative k for negative Omega
    return k;
  }
}

//======================================================================
/// Compute the JONSWAP wave spectrum
//======================================================================
float JONSWAP(float Omega,float Hs,float Tp,float Gamma) {

  // This FUNCTION computes the JOint North Sea WAve Project
  // (JONSWAP) representation of the one-sided power spectral density
  // or wave spectrum given the frequency, Omega, peak shape
  // parameter, Gamma, significant wave height, Hs, and peak spectral
  // period, Tp, as inputs.  If the value of Gamma is 1.0, the
  // Pierson-Moskowitz wave spectrum is returned.
  //
  // There are several different versions of the JONSWAP spectrum
  // formula.  This version is based on the one documented in the
  // IEC61400-3 wind turbine design standard for offshore wind
  // turbines.

  // Local Variables:

  float Alpha;                          // Exponent on Gamma used in the spectral formulation (-)
  float C;                              // Normalising factor used in the spectral formulation (-)
  float f;                              // Wave frequency (Hz)
  float fp;                             // Peak spectral frequency (Hz)
  float fpOvrf4;                        // (fp/f)^4
  float Sigma;                          // Scaling factor used in the spectral formulation (-)

  // Compute the JONSWAP wave spectrum, unless Omega is zero, in which case,
  //   return zero:

  if(Omega == 0.0) { return 0.0; }  // When .TRUE., the formulation below is ill-conditioned; thus, the known value of zero is returned.
  else {  // Omega > 0.0; forumulate the JONSWAP spectrum.
    // Compute the wave frequency and peak spectral frequency in Hz:
    f=float(1.f/2.f/float(PI*Omega));
    fp=1.f/Tp;
    fpOvrf4=pow((fp/f),4.0f);

    // Compute the normalising factor:
    C=1.0f-(0.287f*log(Gamma));

    // Compute Alpha:
    if(f <= fp) { Sigma=float(0.07); }
    else { Sigma=float(0.09); }
    Alpha=exp((-0.5f *(((f/fp)-1.0f)/Sigma)*(((f/fp)-1.0f)/Sigma)));
    // Compute the wave spectrum:
    return (1.f/2.f/float(PI)*C*(0.3125f*Hs*Hs*fpOvrf4/f)*exp((-1.25f*fpOvrf4))*(pow(Gamma,Alpha)));
  }
}

//======================================================================
/// Compute the hyperbolic numerator over denominator
//======================================================================
float SINHNumOvrSIHNDen(float k,float h,float z) {
  //bjj: note, MLB had issues with the IVF 12 compiler inline function expansion
  // and causing "unexpected results" (NaN) here. possible error in optimazations.
  //MLB: I turned off in-line function expansion to eliminate the NaN issue with the COSH/SUNH

   // This FUNCTION computes the shallow water hyperbolic numerator
   // over denominator term in the wave kinematics expressions:
   //
   //                    SINH( k*( z + h ) )/SINH( k*h )
   //
   // given the wave number, k, water depth, h, and elevation z, as
   // inputs.

   // Compute the hyperbolic numerator over denominator:
  float SINHNumOvrSIHNDen;

  if(k == 0.0) { SINHNumOvrSIHNDen=1.0; }  // When .TRUE., the shallow water formulation is ill-conditioned; thus, the known value of unity is returned.
  else if(k*h > 89.4) { SINHNumOvrSIHNDen=exp(k*z); }   // When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, with h > 14.23*wavelength (since k=2*Pi/wavelength) we can use the numerically-stable deep water formulation instead.
  else if(-k*h > 89.4) { SINHNumOvrSIHNDen=-exp(-k*z); }  // @mth: added negative k case
  else { SINHNumOvrSIHNDen=sinh(k*(z + h))/sinh(k*h); } // 0<k*h <= 89.4; use the shallow water formulation.

  return SINHNumOvrSIHNDen;
}
//======================================================================
/// Compute the hyperbolic numerator over denominator
//======================================================================
float COSHNumOvrSIHNDen(float k,float h,float z) {
  // This FUNCTION computes the shallow water hyperbolic numerator
  // over denominator term in the wave kinematics expressions:
  //
  //                    COSH( k*( z + h ) )/SINH( k*h )
  //
  // given the wave number, k, water depth, h, and elevation z, as
  // inputs.

  // Compute the hyperbolic numerator over denominator:
  float COSHNumOvrSIHNDen;
  //k=abs(k);  //@mth: added 

  if(k == 0.0) { COSHNumOvrSIHNDen=99999; } // When .TRUE., the shallow water formulation is ill-conditioned; thus, HUGE(k) is returned to approximate the known value of infinity.
  else if(k*h > 89.4) { COSHNumOvrSIHNDen=exp(k*z); }// When .TRUE., the shallow water formulation will trigger a floating point overflow error; however, with h > 14.23*wavelength (since k=2*Pi/wavelength) we can use the numerically-stable deep water formulation instead.  
  else if(-k*h > 89.4) { COSHNumOvrSIHNDen=-exp(-k*z); } // @mth: added negative k case
  else { COSHNumOvrSIHNDen=cosh(k*(z + h))/sinh(k*h); }  // 0<k*h <= 89.4; use the shallow water formulation.

  return COSHNumOvrSIHNDen;
}

// ==================== some filtering functions used in OC4 study =============================

//======================================================================
/// Flip or reverse function
//======================================================================
void reverse(double* data,int datasize) {
  int i;
  double temp;
  for(i=0; i<floor(datasize/2); i++) {
    temp=data[i];
    data[i]=data[datasize-1-i];
    data[datasize-1-i]=temp;
  }
  return;
}

//======================================================================
/// IIR filter function
//======================================================================
void doIIR(double* in,double* out,int dataSize,double* a,double* b,int kernelSize) {

  // kernelSize is size of vectors a and b, which contain the coefficients.  Normally a[0]=1.
  int i,k;
  // filtering from out[0] to out[kernelSize-2]
  for(i=0; i<kernelSize-1; ++i) {
    out[i]=in[kernelSize-1];          // maybe this hack will work...
//    out[i]=in[i]*b[0];                             // using just part of a filter kernel doesn't really work
//    for(k=1; k<i+1; k++)
//      out[i]+=in[i-k]*b[k]-out[i-k]*a[k];
//    out[i]=out[i]/a[0];
  }

  // filtering from out[kernelSize-1] to out[dataSize-1] (last)
  for(i=kernelSize-1; i<dataSize; i++) {
    out[i]=in[i]*b[0];
    for(k=1; k<kernelSize; k++) { out[i]+=in[i-k]*b[k]-out[i-k]*a[k]; }
    out[i]=out[i]/a[0];
  }
  return;
}

//======================================================================
/// State Space filter function-a good reference is https://ccrma.stanford.edu/~jos/fp/Converting_State_Space_Form_Hand.html
//=====================================================================
void doSSfilter(double* in,double* out,int dataSize,double* a,double* beta,double b0,int kernelSize) {
  // kernelSize is size of vectors a and beta (excluding 0th entry!),which contain the coefficients.

  int i,k;
  double* fstates=(double*)malloc(dataSize*sizeof(double)); // create state vector

  // filtering from out[0] to out[kernelSize-2]
  for(i=0; i<kernelSize; ++i) {
    out[i]=in[kernelSize];      // maybe this hack will work...
    fstates[i]=0.0;
  }

  // filtering from out[kernelSize] to out[dataSize-1] (last)
  for(i=kernelSize; i<dataSize; i++) {
    fstates[i]=in[i-1];
    out[i]=b0*in[i];
    for(k=0; k<kernelSize; k++) { fstates[i]+=-1.0*a[k]*fstates[i-1-k]; }
    for(k=0; k<kernelSize; k++) { out[i]+=beta[k]*fstates[i-k]; }
  }

  free(fstates);
  return;
}
} //end misc