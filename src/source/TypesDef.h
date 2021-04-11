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

/// \file TypesDef.h \brief Declares general types and functions for the entire application.

#ifndef _TypesDef_
#define _TypesDef_

#define PI 3.14159265358979323846      ///<Value of cte PI. 
#define TWOPI 6.28318530717958647692   ///<Value of cte PI*2. 
#define PIHALF 1.57079632679489661923  ///<Value of cte PI/2. 
#define TORAD 0.017453292519943295769  ///<Constant for conversion to radians. rad=degrees*TORAD (TORAD=PI/180)
#define TODEG 57.29577951308232087684  ///<Constant for conversion to degrees. degrees=rad*TODEG (TODEG=180/PI)
#define EULER 2.71828182845904523536   ///<Value of cte E (Euler's number or Napier's constant), E=std::exp(1.0);

typedef unsigned char byte;
typedef unsigned short word;
typedef long long llong;
typedef unsigned long long ullong;

//##############################################################################
//# Basic data type definitions.
//##############################################################################

///Structure of 2 variables of type unsigned.
typedef struct{
  int x,y;
}tint2;

inline tint2 TInt2(int v){ tint2 p={v,v}; return(p); }
inline tint2 TInt2(int x,int y){ tint2 p={x,y}; return(p); }
inline bool  operator ==(const tint2& a,const tint2& b){ return(a.x==b.x && a.y==b.y); }
inline bool  operator !=(const tint2& a,const tint2& b){ return(a.x!=b.x || a.y!=b.y); }
inline bool  operator  <(const tint2& a,const tint2& b){ return(a.x<b.x && a.y<b.y); }
inline bool  operator  >(const tint2& a,const tint2& b){ return(a.x>b.x && a.y>b.y); }
inline bool  operator <=(const tint2& a,const tint2& b){ return(a.x<=b.x && a.y<=b.y); }
inline bool  operator >=(const tint2& a,const tint2& b){ return(a.x>=b.x && a.y>=b.y); }
inline tint2 operator  +(const tint2& a,const tint2& b){ return(TInt2(a.x+b.x,a.y+b.y)); }
inline tint2 operator  -(const tint2& a,const tint2& b){ return(TInt2(a.x-b.x,a.y-b.y)); }
inline tint2 operator  *(const tint2& a,const tint2& b){ return(TInt2(a.x*b.x,a.y*b.y)); }
inline tint2 operator  /(const tint2& a,const tint2& b){ return(TInt2(a.x/b.x,a.y/b.y)); }
inline tint2 operator  +(const tint2& a,const int& b){ return(TInt2(a.x+b,a.y+b)); }
inline tint2 operator  -(const tint2& a,const int& b){ return(TInt2(a.x-b,a.y-b)); }
inline tint2 operator  *(const tint2& a,const int& b){ return(TInt2(a.x*b,a.y*b)); }
inline tint2 operator  /(const tint2& a,const int& b){ return(TInt2(a.x/b,a.y/b)); }
inline tint2 MinValues(const tint2& a,const tint2& b){ return(TInt2((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y))); }
inline tint2 MaxValues(const tint2& a,const tint2& b){ return(TInt2((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y))); }
inline int TInt2Get(const tint2& a,int c){ return(!c? a.x: a.y); }
inline tint2 TInt2Set(const tint2& a,int c,int v){ return(TInt2((c? a.x: v),(c!=1? a.y: v))); }


///Structure of 2 variables of type unsigned.
typedef struct{
  unsigned x,y;
}tuint2;

inline tuint2 TUint2(unsigned v){ tuint2 p={v,v}; return(p); }
inline tuint2 TUint2(unsigned x,unsigned y){ tuint2 p={x,y}; return(p); }
inline bool   operator ==(const tuint2& a,const tuint2& b){ return(a.x==b.x && a.y==b.y); }
inline bool   operator !=(const tuint2& a,const tuint2& b){ return(a.x!=b.x || a.y!=b.y); }
inline bool   operator  <(const tuint2& a,const tuint2& b){ return(a.x<b.x && a.y<b.y); }
inline bool   operator  >(const tuint2& a,const tuint2& b){ return(a.x>b.x && a.y>b.y); }
inline bool   operator <=(const tuint2& a,const tuint2& b){ return(a.x<=b.x && a.y<=b.y); }
inline bool   operator >=(const tuint2& a,const tuint2& b){ return(a.x>=b.x && a.y>=b.y); }
inline tuint2 operator  +(const tuint2& a,const tuint2& b){ return(TUint2(a.x+b.x,a.y+b.y)); }
inline tuint2 operator  -(const tuint2& a,const tuint2& b){ return(TUint2(a.x-b.x,a.y-b.y)); }
inline tuint2 operator  *(const tuint2& a,const tuint2& b){ return(TUint2(a.x*b.x,a.y*b.y)); }
inline tuint2 operator  /(const tuint2& a,const tuint2& b){ return(TUint2(a.x/b.x,a.y/b.y)); }
inline tuint2 operator  +(const tuint2& a,const unsigned& b){ return(TUint2(a.x+b,a.y+b)); }
inline tuint2 operator  -(const tuint2& a,const unsigned& b){ return(TUint2(a.x-b,a.y-b)); }
inline tuint2 operator  *(const tuint2& a,const unsigned& b){ return(TUint2(a.x*b,a.y*b)); }
inline tuint2 operator  /(const tuint2& a,const unsigned& b){ return(TUint2(a.x/b,a.y/b)); }
inline tuint2 MinValues(const tuint2& a,const tuint2& b){ return(TUint2((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y))); }
inline tuint2 MaxValues(const tuint2& a,const tuint2& b){ return(TUint2((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y))); }
inline unsigned TUint2Get(const tuint2& a,unsigned c){ return(!c? a.x: a.y); }
inline tuint2 TUint2Set(const tuint2& a,unsigned c,unsigned v){ return(TUint2((c? a.x: v),(c!=1? a.y: v))); }


///Structure of 3 variables of type int.
typedef struct{
  int x,y,z;
}tint3;

inline tint3 TInt3(int v){ tint3 p={v,v,v}; return(p); }
inline tint3 TInt3(int x,int y,int z){ tint3 p={x,y,z}; return(p); }
inline bool   operator ==(const tint3& a, const tint3& b){ return(a.x==b.x && a.y==b.y && a.z==b.z); }
inline bool   operator !=(const tint3& a, const tint3& b){ return(a.x!=b.x || a.y!=b.y || a.z!=b.z); }
inline bool   operator  <(const tint3& a, const tint3& b){ return(a.x <b.x && a.y <b.y && a.z <b.z); }
inline bool   operator  >(const tint3& a, const tint3& b){ return(a.x >b.x && a.y >b.y && a.z >b.z); }
inline bool   operator <=(const tint3& a, const tint3& b){ return(a.x<=b.x && a.y<=b.y && a.z<=b.z); }
inline bool   operator >=(const tint3& a, const tint3& b){ return(a.x>=b.x && a.y>=b.y && a.z>=b.z); }
inline tint3 operator  +(const tint3& a, const tint3& b){ return(TInt3(a.x+b.x,a.y+b.y,a.z+b.z)); }
inline tint3 operator  -(const tint3& a, const tint3& b){ return(TInt3(a.x-b.x,a.y-b.y,a.z-b.z)); }
inline tint3 operator  *(const tint3& a, const tint3& b){ return(TInt3(a.x*b.x,a.y*b.y,a.z*b.z)); }
inline tint3 operator  /(const tint3& a, const tint3& b){ return(TInt3(a.x/b.x,a.y/b.y,a.z/b.z)); }
inline tint3 operator  +(const tint3& a, const int& b){ return(TInt3(a.x+b,a.y+b,a.z+b)); }
inline tint3 operator  -(const tint3& a, const int& b){ return(TInt3(a.x-b,a.y-b,a.z-b)); }
inline tint3 operator  *(const tint3& a, const int& b){ return(TInt3(a.x*b,a.y*b,a.z*b)); }
inline tint3 operator  /(const tint3& a, const int& b){ return(TInt3(a.x/b,a.y/b,a.z/b)); }
inline tint3 MinValues(const tint3& a, const tint3& b){ return(TInt3((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y),(a.z<=b.z? a.z: b.z))); }
inline tint3 MaxValues(const tint3& a, const tint3& b){ return(TInt3((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y),(a.z>=b.z? a.z: b.z))); }
inline int TInt3Get(const tint3& a,int c){ return(!c? a.x: (c==1? a.y: a.z)); }
inline tint3 TInt3Set(const tint3& a,int c,int v){ return(TInt3((c? a.x: v),(c!=1? a.y: v),(c!=2? a.z: v))); }


///Structure of 3 variables of type unsigned.
typedef struct{
  unsigned x,y,z;
}tuint3;

inline tuint3 TUint3(unsigned v){ tuint3 p={v,v,v}; return(p); }
inline tuint3 TUint3(unsigned x,unsigned y,unsigned z){ tuint3 p={x,y,z}; return(p); }
inline bool   operator ==(const tuint3& a, const tuint3& b){ return(a.x==b.x && a.y==b.y && a.z==b.z); }
inline bool   operator !=(const tuint3& a, const tuint3& b){ return(a.x!=b.x || a.y!=b.y || a.z!=b.z); }
inline bool   operator  <(const tuint3& a, const tuint3& b){ return(a.x <b.x && a.y <b.y && a.z <b.z); }
inline bool   operator  >(const tuint3& a, const tuint3& b){ return(a.x >b.x && a.y >b.y && a.z >b.z); }
inline bool   operator <=(const tuint3& a, const tuint3& b){ return(a.x<=b.x && a.y<=b.y && a.z<=b.z); }
inline bool   operator >=(const tuint3& a, const tuint3& b){ return(a.x>=b.x && a.y>=b.y && a.z>=b.z); }
inline tuint3 operator  +(const tuint3& a, const tuint3& b){ return(TUint3(a.x+b.x,a.y+b.y,a.z+b.z)); }
inline tuint3 operator  -(const tuint3& a, const tuint3& b){ return(TUint3(a.x-b.x,a.y-b.y,a.z-b.z)); }
inline tuint3 operator  *(const tuint3& a, const tuint3& b){ return(TUint3(a.x*b.x,a.y*b.y,a.z*b.z)); }
inline tuint3 operator  /(const tuint3& a, const tuint3& b){ return(TUint3(a.x/b.x,a.y/b.y,a.z/b.z)); }
inline tuint3 operator  +(const tuint3& a, const unsigned& b){ return(TUint3(a.x+b,a.y+b,a.z+b)); }
inline tuint3 operator  -(const tuint3& a, const unsigned& b){ return(TUint3(a.x-b,a.y-b,a.z-b)); }
inline tuint3 operator  *(const tuint3& a, const unsigned& b){ return(TUint3(a.x*b,a.y*b,a.z*b)); }
inline tuint3 operator  /(const tuint3& a, const unsigned& b){ return(TUint3(a.x/b,a.y/b,a.z/b)); }
inline tuint3 MinValues(const tuint3& a, const tuint3& b){ return(TUint3((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y),(a.z<=b.z? a.z: b.z))); }
inline tuint3 MaxValues(const tuint3& a, const tuint3& b){ return(TUint3((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y),(a.z>=b.z? a.z: b.z))); }
inline unsigned TUint3Get(const tuint3& a,unsigned c){ return(!c? a.x: (c==1? a.y: a.z)); }
inline tuint3 TUint3Set(const tuint3& a,unsigned c,unsigned v){ return(TUint3((c? a.x: v),(c!=1? a.y: v),(c!=2? a.z: v))); }


///Structure of 2 variables of type float.
typedef struct{
  float x,y;
}tfloat2;

inline tfloat2 TFloat2(float v){ tfloat2 p={v,v}; return(p); }
inline tfloat2 TFloat2(float x,float y){ tfloat2 p={x,y}; return(p); }
inline bool operator ==(const tfloat2& a, const tfloat2& b){ return(a.x==b.x&&a.y==b.y); }
inline bool operator !=(const tfloat2& a, const tfloat2& b){ return(a.x!=b.x||a.y!=b.y); }
inline bool operator <(const tfloat2& a, const tfloat2& b){ return(a.x<b.x&&a.y<b.y); }
inline bool operator >(const tfloat2& a, const tfloat2& b){ return(a.x>b.x&&a.y>b.y); }
inline bool operator <=(const tfloat2& a, const tfloat2& b){ return(a.x<=b.x&&a.y<=b.y); }
inline bool operator >=(const tfloat2& a, const tfloat2& b){ return(a.x>=b.x&&a.y>=b.y); }
inline tfloat2 operator +(const tfloat2& a, const tfloat2& b){ return(TFloat2(a.x+b.x,a.y+b.y)); }
inline tfloat2 operator -(const tfloat2& a, const tfloat2& b){ return(TFloat2(a.x-b.x,a.y-b.y)); }
inline tfloat2 operator *(const tfloat2& a, const tfloat2& b){ return(TFloat2(a.x*b.x,a.y*b.y)); }
inline tfloat2 operator /(const tfloat2& a, const tfloat2& b){ return(TFloat2(a.x/b.x,a.y/b.y)); }
inline tfloat2 operator +(const tfloat2& a, const float& b){ return(TFloat2(a.x+b,a.y+b)); }
inline tfloat2 operator -(const tfloat2& a, const float& b){ return(TFloat2(a.x-b,a.y-b)); }
inline tfloat2 operator *(const tfloat2& a, const float& b){ return(TFloat2(a.x*b,a.y*b)); }
inline tfloat2 operator /(const tfloat2& a, const float& b){ return(TFloat2(a.x/b,a.y/b)); }
inline tfloat2 MinValues(const tfloat2& a, const tfloat2& b){ return(TFloat2((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y))); }
inline tfloat2 MaxValues(const tfloat2& a, const tfloat2& b){ return(TFloat2((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y))); }


///Structure of 3 variables of type float.
typedef struct{
  float x,y,z;
}tfloat3;

inline tfloat3 TFloat3(float v){ tfloat3 p={v,v,v}; return(p); }
inline tfloat3 TFloat3(float x,float y,float z){ tfloat3 p={x,y,z}; return(p); }
inline bool operator ==(const tfloat3& a, const tfloat3& b){ return(a.x==b.x&&a.y==b.y&&a.z==b.z); }
inline bool operator !=(const tfloat3& a, const tfloat3& b){ return(a.x!=b.x||a.y!=b.y||a.z!=b.z); }
inline bool operator <(const tfloat3& a, const tfloat3& b){ return(a.x<b.x&&a.y<b.y&&a.z<b.z); }
inline bool operator >(const tfloat3& a, const tfloat3& b){ return(a.x>b.x&&a.y>b.y&&a.z>b.z); }
inline bool operator <=(const tfloat3& a, const tfloat3& b){ return(a.x<=b.x&&a.y<=b.y&&a.z<=b.z); }
inline bool operator >=(const tfloat3& a, const tfloat3& b){ return(a.x>=b.x&&a.y>=b.y&&a.z>=b.z); }
inline tfloat3 operator +(const tfloat3& a, const tfloat3& b){ return(TFloat3(a.x+b.x,a.y+b.y,a.z+b.z)); }
inline tfloat3 operator -(const tfloat3& a, const tfloat3& b){ return(TFloat3(a.x-b.x,a.y-b.y,a.z-b.z)); }
inline tfloat3 operator *(const tfloat3& a, const tfloat3& b){ return(TFloat3(a.x*b.x,a.y*b.y,a.z*b.z)); }
inline tfloat3 operator /(const tfloat3& a, const tfloat3& b){ return(TFloat3(a.x/b.x,a.y/b.y,a.z/b.z)); }
inline tfloat3 operator +(const tfloat3& a, const float& b){ return(TFloat3(a.x+b,a.y+b,a.z+b)); }
inline tfloat3 operator -(const tfloat3& a, const float& b){ return(TFloat3(a.x-b,a.y-b,a.z-b)); }
inline tfloat3 operator *(const tfloat3& a, const float& b){ return(TFloat3(a.x*b,a.y*b,a.z*b)); }
inline tfloat3 operator /(const tfloat3& a, const float& b){ return(TFloat3(a.x/b,a.y/b,a.z/b)); }
inline tfloat3 MinValues(const tfloat3& a, const tfloat3& b){ return(TFloat3((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y),(a.z<=b.z? a.z: b.z))); }
inline tfloat3 MaxValues(const tfloat3& a, const tfloat3& b){ return(TFloat3((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y),(a.z>=b.z? a.z: b.z))); }
inline float TFloat3Get(const tfloat3& a,unsigned c){ return(!c? a.x: (c==1? a.y: a.z)); }
inline tfloat3 TFloat3Set(const tfloat3& a,unsigned c,float v){ return(TFloat3((c? a.x: v),(c!=1? a.y: v),(c!=2? a.z: v))); }


///Structure of 2 variables of type double.
typedef struct{
  double x,y;
}tdouble2;

inline tdouble2 TDouble2(double v){ tdouble2 p={v,v}; return(p); }
inline tdouble2 TDouble2(double x,double y){ tdouble2 p={x,y}; return(p); }
inline bool operator ==(const tdouble2& a, const tdouble2& b){ return(a.x==b.x&&a.y==b.y); }
inline bool operator !=(const tdouble2& a, const tdouble2& b){ return(a.x!=b.x||a.y!=b.y); }
inline tdouble2 operator +(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x+b.x,a.y+b.y)); }
inline tdouble2 operator -(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x-b.x,a.y-b.y)); }
inline tdouble2 operator *(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x*b.x,a.y*b.y)); }
inline tdouble2 operator /(const tdouble2& a, const tdouble2& b){ return(TDouble2(a.x/b.x,a.y/b.y)); }
inline tdouble2 operator +(const tdouble2& a, const double& b){ return(TDouble2(a.x+b,a.y+b)); }
inline tdouble2 operator -(const tdouble2& a, const double& b){ return(TDouble2(a.x-b,a.y-b)); }
inline tdouble2 operator *(const tdouble2& a, const double& b){ return(TDouble2(a.x*b,a.y*b)); }
inline tdouble2 operator /(const tdouble2& a, const double& b){ return(TDouble2(a.x/b,a.y/b)); }
inline tdouble2 MinValues (const tdouble2& a, const tdouble2& b){ return(TDouble2((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y))); }
inline tdouble2 MaxValues (const tdouble2& a, const tdouble2& b){ return(TDouble2((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y))); }


///Structure of 3 variables of type double.
typedef struct{
  double x,y,z;
}tdouble3;

inline tdouble3 TDouble3(double v){ tdouble3 p={v,v,v}; return(p); }
inline tdouble3 TDouble3(double x,double y,double z){ tdouble3 p={x,y,z}; return(p); }
inline bool operator ==(const tdouble3& a, const tdouble3& b){ return(a.x==b.x&&a.y==b.y&&a.z==b.z); }
inline bool operator !=(const tdouble3& a, const tdouble3& b){ return(a.x!=b.x||a.y!=b.y||a.z!=b.z); }
inline bool operator <(const tdouble3& a, const tdouble3& b){ return(a.x<b.x&&a.y<b.y&&a.z<b.z); }
inline bool operator >(const tdouble3& a, const tdouble3& b){ return(a.x>b.x&&a.y>b.y&&a.z>b.z); }
inline bool operator <=(const tdouble3& a, const tdouble3& b){ return(a.x<=b.x&&a.y<=b.y&&a.z<=b.z); }
inline bool operator >=(const tdouble3& a, const tdouble3& b){ return(a.x>=b.x&&a.y>=b.y&&a.z>=b.z); }
inline tdouble3 operator +(const tdouble3& a, const tdouble3& b){ return(TDouble3(a.x+b.x,a.y+b.y,a.z+b.z)); }
inline tdouble3 operator -(const tdouble3& a, const tdouble3& b){ return(TDouble3(a.x-b.x,a.y-b.y,a.z-b.z)); }
inline tdouble3 operator *(const tdouble3& a, const tdouble3& b){ return(TDouble3(a.x*b.x,a.y*b.y,a.z*b.z)); }
inline tdouble3 operator /(const tdouble3& a, const tdouble3& b){ return(TDouble3(a.x/b.x,a.y/b.y,a.z/b.z)); }
inline tdouble3 operator +(const tdouble3& a, const double& b){ return(TDouble3(a.x+b,a.y+b,a.z+b)); }
inline tdouble3 operator -(const tdouble3& a, const double& b){ return(TDouble3(a.x-b,a.y-b,a.z-b)); }
inline tdouble3 operator *(const tdouble3& a, const double& b){ return(TDouble3(a.x*b,a.y*b,a.z*b)); }
inline tdouble3 operator /(const tdouble3& a, const double& b){ return(TDouble3(a.x/b,a.y/b,a.z/b)); }
inline tdouble3 MinValues(const tdouble3& a, const tdouble3& b){ return(TDouble3((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y),(a.z<=b.z? a.z: b.z))); }
inline tdouble3 MaxValues(const tdouble3& a, const tdouble3& b){ return(TDouble3((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y),(a.z>=b.z? a.z: b.z))); }

///Converts \ref tuint3 to \ref tint3.
inline tint3 ToTInt3(const tuint3& v){ return(TInt3(int(v.x),int(v.y),int(v.z))); }
///Converts \ref tint3 to \ref tuint3.
inline tuint3 ToTUint3(const tint3& v){ return(TUint3(unsigned(v.x),unsigned(v.y),unsigned(v.z))); }

///Converts \ref tdouble2 to \ref tfloat2.
inline tfloat2 ToTFloat2(const tdouble2& v){ return(TFloat2(float(v.x),float(v.y))); }
///Converts \ref tfloat2 to \ref tdouble2.
inline tdouble2 ToTDouble2(const tfloat2& v){ return(TDouble2(v.x,v.y)); }

///Converts \ref tdouble3 to \ref tfloat3.
inline tfloat3 ToTFloat3(const tdouble3& v){ return(TFloat3(float(v.x),float(v.y),float(v.z))); }
///Converts \ref tfloat3 to \ref tdouble3.
inline tdouble3 ToTDouble3(const tfloat3& v){ return(TDouble3(v.x,v.y,v.z)); }


///Structure of 4 variables of type int.
typedef struct{
  int x,y,z,w;
}tint4;

inline tint4 TInt4(int v){ tint4 p={v,v,v,v}; return(p); }
inline tint4 TInt4(int x,int y,int z,int w){ tint4 p={x,y,z,w}; return(p); }
inline bool operator ==(const tint4& a, const tint4& b){ return(a.x==b.x&&a.y==b.y&&a.z==b.z&&a.w==b.w); }
inline bool operator !=(const tint4& a, const tint4& b){ return(a.x!=b.x||a.y!=b.y||a.z!=b.z||a.w!=b.w); }
inline tint4 operator +(const tint4& a, const tint4& b){ return(TInt4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w)); }
inline tint4 operator -(const tint4& a, const tint4& b){ return(TInt4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w)); }
inline tint4 operator *(const tint4& a, const tint4& b){ return(TInt4(a.x*b.x,a.y*b.y,a.z*b.z,a.w*b.w)); }
inline tint4 operator /(const tint4& a, const tint4& b){ return(TInt4(a.x/b.x,a.y/b.y,a.z/b.z,a.w/b.w)); }


///Structure of 4 variables of type unsigned.
typedef struct{
  unsigned x,y,z,w;
}tuint4;

inline tuint4 TUint4(unsigned v){ tuint4 p={v,v,v,v}; return(p); }
inline tuint4 TUint4(unsigned x,unsigned y,unsigned z,unsigned w){ tuint4 p={x,y,z,w}; return(p); }
inline bool operator ==(const tuint4& a, const tuint4& b){ return(a.x==b.x&&a.y==b.y&&a.z==b.z&&a.w==b.w); }
inline bool operator !=(const tuint4& a, const tuint4& b){ return(a.x!=b.x||a.y!=b.y||a.z!=b.z||a.w!=b.w); }
inline tuint4 operator +(const tuint4& a, const tuint4& b){ return(TUint4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w)); }
inline tuint4 operator -(const tuint4& a, const tuint4& b){ return(TUint4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w)); }
inline tuint4 operator *(const tuint4& a, const tuint4& b){ return(TUint4(a.x*b.x,a.y*b.y,a.z*b.z,a.w*b.w)); }
inline tuint4 operator /(const tuint4& a, const tuint4& b){ return(TUint4(a.x/b.x,a.y/b.y,a.z/b.z,a.w/b.w)); }


///Structure of 4 variables of type float.
typedef struct{
  float x,y,z,w;
}tfloat4;

inline tfloat4 TFloat4(float v){ tfloat4 p={v,v,v,v}; return(p); }
inline tfloat4 TFloat4(float x,float y,float z,float w){ tfloat4 p={x,y,z,w}; return(p); }
inline bool operator ==(const tfloat4& a, const tfloat4& b){ return(a.x==b.x&&a.y==b.y&&a.z==b.z&&a.w==b.w); }
inline bool operator !=(const tfloat4& a, const tfloat4& b){ return(a.x!=b.x||a.y!=b.y||a.z!=b.z||a.w!=b.w); }
inline tfloat4 operator +(const tfloat4& a, const tfloat4& b){ return(TFloat4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w)); }
inline tfloat4 operator -(const tfloat4& a, const tfloat4& b){ return(TFloat4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w)); }
inline tfloat4 operator *(const tfloat4& a, const tfloat4& b){ return(TFloat4(a.x*b.x,a.y*b.y,a.z*b.z,a.w*b.w)); }
inline tfloat4 operator /(const tfloat4& a, const tfloat4& b){ return(TFloat4(a.x/b.x,a.y/b.y,a.z/b.z,a.w/b.w)); }
inline tfloat4 MinValues (const tfloat4& a, const tfloat4& b){ return(TFloat4((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y),(a.z<=b.z? a.z: b.z),(a.w<=b.w? a.w: b.w))); }
inline tfloat4 MaxValues (const tfloat4& a, const tfloat4& b){ return(TFloat4((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y),(a.z>=b.z? a.z: b.z),(a.w>=b.w? a.w: b.w))); }


///Structure of 4 variables of type float.
typedef struct{
  double x,y,z,w;
}tdouble4;

inline tdouble4 TDouble4(double v){ tdouble4 p={v,v,v,v}; return(p); }
inline tdouble4 TDouble4(double x,double y,double z,double w){ tdouble4 p={x,y,z,w}; return(p); }
inline bool operator ==(const tdouble4& a, const tdouble4& b){ return(a.x==b.x&&a.y==b.y&&a.z==b.z&&a.w==b.w); }
inline bool operator !=(const tdouble4& a, const tdouble4& b){ return(a.x!=b.x||a.y!=b.y||a.z!=b.z||a.w!=b.w); }
inline tdouble4 operator +(const tdouble4& a, const tdouble4& b){ return(TDouble4(a.x+b.x,a.y+b.y,a.z+b.z,a.w+b.w)); }
inline tdouble4 operator -(const tdouble4& a, const tdouble4& b){ return(TDouble4(a.x-b.x,a.y-b.y,a.z-b.z,a.w-b.w)); }
inline tdouble4 operator *(const tdouble4& a, const tdouble4& b){ return(TDouble4(a.x*b.x,a.y*b.y,a.z*b.z,a.w*b.w)); }
inline tdouble4 operator /(const tdouble4& a, const tdouble4& b){ return(TDouble4(a.x/b.x,a.y/b.y,a.z/b.z,a.w/b.w)); }
inline tdouble4 MinValues (const tdouble4& a, const tdouble4& b){ return(TDouble4((a.x<=b.x? a.x: b.x),(a.y<=b.y? a.y: b.y),(a.z<=b.z? a.z: b.z),(a.w<=b.w? a.w: b.w))); }
inline tdouble4 MaxValues (const tdouble4& a, const tdouble4& b){ return(TDouble4((a.x>=b.x? a.x: b.x),(a.y>=b.y? a.y: b.y),(a.z>=b.z? a.z: b.z),(a.w>=b.w? a.w: b.w))); }

///Converts \ref tdouble4 to \ref tfloat4.
inline tfloat4 ToTFloat4(const tdouble4& v){ return(TFloat4(float(v.x),float(v.y),float(v.z),float(v.w))); }
///Converts \ref tfloat4 to \ref tdouble4.
inline tdouble4 ToTDouble4(const tfloat4& v){ return(TDouble4(v.x,v.y,v.z,v.w)); }


///Matrix of 2x2 values of type float.
typedef struct{
  float a11,a12;
  float a21,a22;
}tmatrix2f;

///Constructor of type \ref matrix2f.
inline tmatrix2f TMatrix2f(){ tmatrix2f m={1,0,0,1}; return(m); }
inline tmatrix2f TMatrix2f(float a11,float a12,float a21,float a22){ tmatrix2f m={a11,a12,a21,a22}; return(m); }
inline tmatrix2f TMatrix2f(float v){ tmatrix2f m={v,v,v,v}; return(m); }
inline bool operator ==(const tmatrix2f& a, const tmatrix2f& b){ return(a.a11==b.a11 && a.a12==b.a12 && a.a21==b.a21 && a.a22==b.a22); }
inline bool operator !=(const tmatrix2f& a, const tmatrix2f& b){ return(a.a11!=b.a11 || a.a12!=b.a12 || a.a21!=b.a21 || a.a22!=b.a22); }


///Matrix of 2x2 values of type double.
typedef struct{
  double a11,a12;
  double a21,a22;
}tmatrix2d;

///Constructor of type \ref matrix2d.
inline tmatrix2d TMatrix2d(){ tmatrix2d m={1,0,0,1}; return(m); }
inline tmatrix2d TMatrix2d(double a11,double a12,double a21,double a22){ tmatrix2d m={a11,a12,a21,a22}; return(m); }
inline tmatrix2d TMatrix2d(double v){ tmatrix2d m={v,v,v,v}; return(m); }
inline bool operator ==(const tmatrix2d& a, const tmatrix2d& b){ return(a.a11==b.a11 && a.a12==b.a12 && a.a21==b.a21 && a.a22==b.a22); }
inline bool operator !=(const tmatrix2d& a, const tmatrix2d& b){ return(a.a11!=b.a11 || a.a12!=b.a12 || a.a21!=b.a21 || a.a22!=b.a22); }

///Converts \ref matrix2d to \ref matrix2f.
inline tmatrix2f ToTMatrix2f(const tmatrix2d& v){ return(TMatrix2f(float(v.a11),float(v.a12),float(v.a21),float(v.a22))); }


///Matrix of 3x3 values of type float.
typedef struct{
  float a11,a12,a13;
  float a21,a22,a23;
  float a31,a32,a33;
}tmatrix3f;

///Constructor of type \ref matrix3f.
inline tmatrix3f TMatrix3f(){ tmatrix3f m={1,0,0,0,1,0,0,0,1}; return(m); }
inline tmatrix3f TMatrix3f(float a11,float a12,float a13,float a21,float a22,float a23,float a31,float a32,float a33){ tmatrix3f m={a11,a12,a13,a21,a22,a23,a31,a32,a33}; return(m); }
inline tmatrix3f TMatrix3f(float v){ tmatrix3f m={v,v,v,v,v,v,v,v,v}; return(m); }
inline bool operator ==(const tmatrix3f& a, const tmatrix3f& b){ return(a.a11==b.a11 && a.a12==b.a12 && a.a13==b.a13 && a.a21==b.a21 && a.a22==b.a22 && a.a23==b.a23 && a.a31==b.a31 && a.a32==b.a32 && a.a33==b.a33); }
inline bool operator !=(const tmatrix3f& a, const tmatrix3f& b){ return(a.a11!=b.a11 || a.a12!=b.a12 || a.a13!=b.a13 || a.a21!=b.a21 || a.a22!=b.a22 || a.a23!=b.a23 || a.a31!=b.a31 || a.a32!=b.a32 || a.a33!=b.a33); }


///Matrix of 3x3 values of type double.
typedef struct{
  double a11,a12,a13;
  double a21,a22,a23;
  double a31,a32,a33;
}tmatrix3d;

///Constructor of type \ref matrix3d.
inline tmatrix3d TMatrix3d(){ tmatrix3d m={1,0,0,0,1,0,0,0,1}; return(m); }
inline tmatrix3d TMatrix3d(double a11,double a12,double a13,double a21,double a22,double a23,double a31,double a32,double a33){ tmatrix3d m={a11,a12,a13,a21,a22,a23,a31,a32,a33}; return(m); }
inline tmatrix3d TMatrix3d(double v){ tmatrix3d m={v,v,v,v,v,v,v,v,v}; return(m); }
inline bool operator ==(const tmatrix3d& a, const tmatrix3d& b){ return(a.a11==b.a11 && a.a12==b.a12 && a.a13==b.a13 && a.a21==b.a21 && a.a22==b.a22 && a.a23==b.a23 && a.a31==b.a31 && a.a32==b.a32 && a.a33==b.a33); }
inline bool operator !=(const tmatrix3d& a, const tmatrix3d& b){ return(a.a11!=b.a11 || a.a12!=b.a12 || a.a13!=b.a13 || a.a21!=b.a21 || a.a22!=b.a22 || a.a23!=b.a23 || a.a31!=b.a31 || a.a32!=b.a32 || a.a33!=b.a33); }

///Converts \ref matrix3d to \ref matrix3f.
inline tmatrix3f ToTMatrix3f(const tmatrix3d& v){ return(TMatrix3f(float(v.a11),float(v.a12),float(v.a13),float(v.a21),float(v.a22),float(v.a23),float(v.a31),float(v.a32),float(v.a33))); }


///Matrix of 4x4 values of type float.
typedef struct{
  float a11,a12,a13,a14;
  float a21,a22,a23,a24;
  float a31,a32,a33,a34;
  float a41,a42,a43,a44;
}tmatrix4f;

///Constructor of type \ref matrix4f.
inline tmatrix4f TMatrix4f(){ tmatrix4f m={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}; return(m); }
inline tmatrix4f TMatrix4f(float a11,float a12,float a13,float a14,float a21,float a22,float a23,float a24,float a31,float a32,float a33,float a34,float a41,float a42,float a43,float a44){ tmatrix4f m={a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44}; return(m); }
inline tmatrix4f TMatrix4f(float v){ tmatrix4f m={v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v}; return(m); }
inline bool operator ==(const tmatrix4f& a, const tmatrix4f& b){ return(a.a11==b.a11 && a.a12==b.a12 && a.a13==b.a13 && a.a14==b.a14 && a.a21==b.a21 && a.a22==b.a22 && a.a23==b.a23 && a.a24==b.a24 && a.a31==b.a31 && a.a32==b.a32 && a.a33==b.a33 && a.a34==b.a34 && a.a41==b.a41 && a.a42==b.a42 && a.a43==b.a43 && a.a44==b.a44); }
inline bool operator !=(const tmatrix4f& a, const tmatrix4f& b){ return(a.a11!=b.a11 || a.a12!=b.a12 || a.a13!=b.a13 || a.a14!=b.a14 || a.a21!=b.a21 || a.a22!=b.a22 || a.a23!=b.a23 || a.a24!=b.a24 || a.a31!=b.a31 || a.a32!=b.a32 || a.a33!=b.a33 || a.a34!=b.a34 || a.a41!=b.a41 || a.a42!=b.a42 || a.a43!=b.a43 || a.a44!=b.a44); }
inline tfloat3 MatrixMulPoint(const tmatrix4f &m,const tfloat3 &p){ return(TFloat3(m.a11*p.x + m.a12*p.y + m.a13*p.z + m.a14, m.a21*p.x + m.a22*p.y + m.a23*p.z + m.a24, m.a31*p.x + m.a32*p.y + m.a33*p.z + m.a34)); }


///Matrix of 4x4 values of type double.
typedef struct{
  double a11,a12,a13,a14;
  double a21,a22,a23,a24;
  double a31,a32,a33,a34;
  double a41,a42,a43,a44;
}tmatrix4d;

///Constructor of type \ref matrix4d.
inline tmatrix4d TMatrix4d(){ tmatrix4d m={1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1}; return(m); }
inline tmatrix4d TMatrix4d(double a11,double a12,double a13,double a14,double a21,double a22,double a23,double a24,double a31,double a32,double a33,double a34,double a41,double a42,double a43,double a44){ tmatrix4d m={a11,a12,a13,a14,a21,a22,a23,a24,a31,a32,a33,a34,a41,a42,a43,a44}; return(m); }
inline tmatrix4d TMatrix4d(double v){ tmatrix4d m={v,v,v,v,v,v,v,v,v,v,v,v,v,v,v,v}; return(m); }
inline bool operator ==(const tmatrix4d& a, const tmatrix4d& b){ return(a.a11==b.a11 && a.a12==b.a12 && a.a13==b.a13 && a.a14==b.a14 && a.a21==b.a21 && a.a22==b.a22 && a.a23==b.a23 && a.a24==b.a24 && a.a31==b.a31 && a.a32==b.a32 && a.a33==b.a33 && a.a34==b.a34 && a.a41==b.a41 && a.a42==b.a42 && a.a43==b.a43 && a.a44==b.a44); }
inline bool operator !=(const tmatrix4d& a, const tmatrix4d& b){ return(a.a11!=b.a11 || a.a12!=b.a12 || a.a13!=b.a13 || a.a14!=b.a14 || a.a21!=b.a21 || a.a22!=b.a22 || a.a23!=b.a23 || a.a24!=b.a24 || a.a31!=b.a31 || a.a32!=b.a32 || a.a33!=b.a33 || a.a34!=b.a34 || a.a41!=b.a41 || a.a42!=b.a42 || a.a43!=b.a43 || a.a44!=b.a44); }
inline tdouble3 MatrixMulPoint(const tmatrix4d &m,const tdouble3 &p){ return(TDouble3(m.a11*p.x + m.a12*p.y + m.a13*p.z + m.a14, m.a21*p.x + m.a22*p.y + m.a23*p.z + m.a24, m.a31*p.x + m.a32*p.y + m.a33*p.z + m.a34)); }
inline tfloat3 MatrixMulPointNormal(const tmatrix4d &m,const tfloat3 &p){ return(ToTFloat3(TDouble3(m.a11*p.x + m.a12*p.y + m.a13*p.z, m.a21*p.x + m.a22*p.y + m.a23*p.z, m.a31*p.x + m.a32*p.y + m.a33*p.z))); }


///Symmetric matrix 3x3 of 6 values of type float.
typedef struct{
  float xx,xy,xz,yy,yz,zz;
}tsymatrix3f;
inline tsymatrix3f TSymMatrix3f(){ tsymatrix3f m={0,0,0,0,0,0}; return(m); }

///Symmetric matrix 4x4 of 10 values of type float.
typedef struct{
  float a11,a12,a13,a14 ,a22,a23,a24 ,a33,a34 ,a44;
}tsymatrix4f;
inline tsymatrix4f TSymMatrix4f(){ tsymatrix4f m={0,0,0,0,0,0,0,0,0,0}; return(m); }



//##############################################################################
//# Geometry type and functions
//##############################################################################
inline tdouble3 Point3dxy(const tdouble2 &p){ return(TDouble3(p.x,p.y,0)); }
inline tdouble3 Point3dxz(const tdouble2 &p){ return(TDouble3(p.x,0,p.y)); }
inline tdouble3 Point3dxy(const tfloat2  &p){ return(TDouble3(p.x,p.y,0)); }
inline tdouble3 Point3dxz(const tfloat2  &p){ return(TDouble3(p.x,0,p.y)); }
inline tfloat3  Point3fxy(const tfloat2  &p){ return(TFloat3 (p.x,p.y,0)); }
inline tfloat3  Point3fxz(const tfloat2  &p){ return(TFloat3 (p.x,0,p.y)); }
inline tfloat3  Point3fxy(const tdouble2 &p){ return(TFloat3 (float(p.x),float(p.y),0)); }
inline tfloat3  Point3fxz(const tdouble2 &p){ return(TFloat3 (float(p.x),0,float(p.y))); }

///Plane definition on 3D using double values.
typedef struct{
  double a,b,c,d;
}tplane3d;

///Plane definition on 3D using float values.
typedef struct{
  float a,b,c,d;
}tplane3f;


inline tplane3d TPlane3d(double v){ tplane3d p={v,v,v,v}; return(p); }
inline tplane3d TPlane3d(double a,double b,double c,double d){ tplane3d p={a,b,c,d}; return(p); }
inline tplane3d TPlane3d(const tdouble4 &v){ return(TPlane3d(v.x,v.y,v.z,v.w)); }
inline tplane3d TPlane3d(const tplane3f &v){ return(TPlane3d(v.a,v.b,v.c,v.d)); }

inline tplane3f TPlane3f(float v){ tplane3f p={v,v,v,v}; return(p); }
inline tplane3f TPlane3f(float a,float b,float c,float d){ tplane3f p={a,b,c,d}; return(p); }
inline tplane3f TPlane3f(const tfloat4 &v){ return(TPlane3f(v.x,v.y,v.z,v.w)); }
inline tplane3f TPlane3f(const tplane3d &v){ return(TPlane3f(float(v.a),float(v.b),float(v.c),float(v.d))); }

inline tfloat4  TPlane3fToTFloat4 (const tplane3f &v){ return(TFloat4(v.a,v.b,v.c,v.d)); }
inline tfloat4  TPlane3dToTFloat4 (const tplane3d &v){ return(TPlane3fToTFloat4(TPlane3f(v))); }
inline tdouble4 TPlane3fToTDouble4(const tplane3f &v){ return(TDouble4(v.a,v.b,v.c,v.d)); }
inline tdouble4 TPlane3dToTDouble4(const tplane3d &v){ return(TDouble4(v.a,v.b,v.c,v.d)); }

///Line definition on 3D using double values.
typedef struct{
  tdouble3 p; ///<Point of rect.
  tdouble3 v; ///<Vector of rect.
}tline3d;

inline tline3d TLine3d(tdouble3 pp,tdouble3 vv){ tline3d r={pp,vv}; return(r); }

///Line definition on 2D using double values.
typedef struct{
  double a,b,c;
}tline2d;

inline tline2d TLine2d(double a,double b,double c){ tline2d r={a,b,c}; return(r); }


//##############################################################################
//# Basic data types.
//##############################################################################

///Basic data types.
typedef enum{ 
  //-DimOfType() -> 0
  TypeNull=0,TypeText=1,TypeBool=2
  //-DimOfType() -> 1
  ,TypeChar  =10, TypeUchar  =11
  ,TypeShort =12, TypeUshort =13
  ,TypeInt   =14, TypeUint   =15
  ,TypeLlong =16, TypeUllong =17
  ,TypeFloat =18, TypeDouble =19
  //-DimOfType() -> 2
  ,TypeInt2  =40, TypeUint2  =41
  ,TypeFloat2=42, TypeDouble2=43 
  //-DimOfType() -> 3
  ,TypeInt3  =60, TypeUint3  =61
  ,TypeFloat3=62, TypeDouble3=63 
  //-DimOfType() -> 4
  ,TypeInt4  =80, TypeUint4  =81
  ,TypeFloat4=82, TypeDouble4=83 
  //-DimOfType() -> 6
  ,TypeSyMatrix3f=100
}TpTypeData; 

/// Returns data type in text.
inline const char* TypeToStr(TpTypeData type){
  switch(type){
    case TypeNull:       return("null");
    case TypeText:       return("text");
    case TypeBool:       return("bool");
    case TypeChar:       return("char");
    case TypeUchar:      return("uchar");
    case TypeShort:      return("short");
    case TypeUshort:     return("ushort");
    case TypeInt:        return("int");
    case TypeUint:       return("uint");
    case TypeLlong:      return("llong");
    case TypeUllong:     return("ullong");
    case TypeFloat:      return("float");
    case TypeDouble:     return("double");
    case TypeInt2:       return("int2");
    case TypeUint2:      return("uint2");
    case TypeFloat2:     return("float2");
    case TypeDouble2:    return("double2");
    case TypeInt3:       return("int3");
    case TypeUint3:      return("uint3");
    case TypeFloat3:     return("float3");
    case TypeDouble3:    return("double3");
    case TypeInt4:       return("int4");
    case TypeUint4:      return("uint4");
    case TypeFloat4:     return("float4");
    case TypeDouble4:    return("double4");
    case TypeSyMatrix3f: return("tsymatrix3f");
  }
  return("???");
}

///Returns size of the data type.
inline unsigned SizeOfType(TpTypeData type){
  switch(type){
    case TypeBool:       return(sizeof(int));
    case TypeChar:       return(sizeof(char));
    case TypeUchar:      return(sizeof(unsigned char));
    case TypeShort:      return(sizeof(short));
    case TypeUshort:     return(sizeof(unsigned short));
    case TypeInt:        return(sizeof(int));
    case TypeUint:       return(sizeof(unsigned));
    case TypeLlong:      return(sizeof(llong));
    case TypeUllong:     return(sizeof(ullong));
    case TypeFloat:      return(sizeof(float));
    case TypeDouble:     return(sizeof(double));
    case TypeInt2:       return(sizeof(tint2));
    case TypeUint2:      return(sizeof(tuint2));
    case TypeFloat2:     return(sizeof(tfloat2));
    case TypeDouble2:    return(sizeof(tdouble2));
    case TypeInt3:       return(sizeof(tint3));
    case TypeUint3:      return(sizeof(tuint3));
    case TypeFloat3:     return(sizeof(tfloat3));
    case TypeDouble3:    return(sizeof(tdouble3));
    case TypeInt4:       return(sizeof(tint4));
    case TypeUint4:      return(sizeof(tuint4));
    case TypeFloat4:     return(sizeof(tfloat4));
    case TypeDouble4:    return(sizeof(tdouble4));
    case TypeSyMatrix3f: return(sizeof(tsymatrix3f));
  }
  return(0);
}

///Returns number of components.
inline int DimOfType(TpTypeData type){
  return(type<TypeChar? 0: (type<TypeInt2? 1: (type<TypeInt3? 2: (type<TypeInt4? 3: (type<TypeSyMatrix3f? 4: 6)))));
}

//-Standard vector types.
#include <vector>
typedef std::vector<tfloat2>  vfloat2;
typedef std::vector<tfloat3>  vfloat3;
typedef std::vector<tdouble2> vdouble2;
typedef std::vector<tdouble3> vdouble3;

#endif


