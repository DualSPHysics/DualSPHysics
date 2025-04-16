/*
 * JSphCpuBuffer.cpp

 *
 *  Created on: Nov 9, 2021
 *      Author: francesco
 */



#include "JSphCpu.h"
#include "JCellSearch_inline.h"
#include "FunSphKernel.h"
#include "FunctionsMath.h"
#include "FunGeo3d.h"
#include "JSphCpuSingle_VRes.h"
#include "JSphCpu_VRes.h"
#include "RunExceptionDef.h"
#include "JMatrix4.h"
#include <climits>
#include <cstring>
#include <algorithm>



using namespace std;

namespace fvres{

//------------------------------------------------------------------------------
/// Solve a linear system of m_dim unknown and m_dim2 right hand sides with a LU Decomposition.
/// Return the reciprocal of the condition number.
//------------------------------------------------------------------------------
template<const int m_dim, const int m_dim2>
void LUdecomp_Single(double (&A)[m_dim*m_dim],int (&p)[m_dim]
  ,double (&b)[m_dim*m_dim2],double (&sol)[m_dim2],double &treshold){
  //-Compute the norm of matrix A
  double maxs = 0;
  for(int i=0;i<m_dim;i++){
    double sum = 0;
    for(int j=0;j<m_dim;j++) sum+=fabs(A[i*m_dim+j]);
    maxs =max(maxs, sum);
  }

  //- Initialize permutation array
  for(int i=0;i<m_dim;i++) p[i]=i;
  

  double maxv = 0;
  for(int i=0;i<m_dim;i++){
    int imax=i;
      for(int k=i;k<m_dim;k++){
        if(fabs(A[i+k*m_dim])>maxv){
          maxv=fabs(A[i+k*m_dim]);
          imax=k;
        }
      }

      if(imax!=i){
        int tmp =p[i];
        p[i]    =p[imax];
        p[imax] =tmp;

        for(int j=0;j<m_dim;j++) {
          double temp          =A[i*m_dim+j];
          A[i*m_dim+j]    =A[imax*m_dim+j];
          A[imax*m_dim+j] =temp;
        }
      }

        // LU decomposition
      for (int j=i+1;j<m_dim;j++){
        A[j*m_dim+i] /= A[i*m_dim+i];
        for (int k=i+1;k<m_dim;k++) A[j*m_dim+k]-=A[j*m_dim+i]*A[i*m_dim+k];
      }
  }

  double ia[m_dim*m_dim] = {0};

  //-Compute inverse matrix
  for (int j=0;j<m_dim;j++) {
    for (int i=0;i<m_dim;i++) {
      ia[i*m_dim+j]= p[i]==j ? 1.0 : 0.0;
      for (int k = 0;k<i;k++)ia[i*m_dim+j]-=A[i*m_dim+k]*ia[k*m_dim+j];
  }
    for (int i=m_dim-1;i>=0;i--){
      for (int k=i+1;k<m_dim;k++)
        ia[i*m_dim+j] -= A[i*m_dim+k]*ia[k*m_dim+j];
        ia[i*m_dim+j] /= A[i*m_dim+i];
    }
  }

  //-Compute norm of inverse matrix
  double maxs1 = 0;
  for (int i=0;i<m_dim;i++) {
    double sum=0;
    for (int j=0;j<m_dim;j++) sum += fabs(ia[i*m_dim+j]);
    maxs1 =max(maxs1,sum);
  }

  treshold =1.0f/(maxs*maxs1);

  //-Compute Solution array
  for (int k=0;k<m_dim2;k++)
    for (int i=0;i<m_dim;i++) sol[k]+=ia[i]*b[i+k*m_dim];
        
}

template<TpKernel tker,bool sim2d,TpVresOrder vrorder,TpVresMethod vrmethod> 
 void InteractionBufferExtrap(unsigned bufferpartcount,const int *bufferpart,
		StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,
		const typecode *code,const unsigned *idp,const tfloat4 *velrhop,const StCteSph csp
    ,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb,const float mrthreshold)
{
  const int n=int(bufferpartcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){

    const int p1=bufferpart[p];
    tdouble3 pos_p1=posb[p1];

    // Compute size of the reconstruction matrix and right hand sime.
    constexpr unsigned m_dim  = vrorder==VrOrder_2nd ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
    constexpr unsigned m_dim2 = sim2d ? 3: 4;

    double C[m_dim]{0};
    double C1[m_dim]{0};
    double D[m_dim2]{0};

    double A[m_dim*m_dim]{0};
    double B[m_dim*m_dim2]{0};

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(pos_p1,false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
    	const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
    	//-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
    	//---------------------------------------------------------------------------------------------
    	for(unsigned p2=pif.x;p2<pif.y;p2++){
    		double drx=double(pos_p1.x-pos[p2].x);
    		double dry=double(pos_p1.y-pos[p2].y);
    		double drz=double(pos_p1.z-pos[p2].z);
    		const double rr2=drx*drx+dry*dry+drz*drz;

    		if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) &&!CODE_IsFluidFixed(code[p2])){//-Only with fluid particles but not inout particles.
    			//-Computes kernel.
    			float fac;
    			float facc;
    			const double wab=fsph::GetKernel_WabFacFacc<tker>(csp,float(rr2),fac,facc);
    			double frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.
    			double frxx=facc*(drx*drx)/(rr2)+fac*(dry*dry+drz*drz)/rr2;
          double frzz=facc*(drz*drz)/(rr2)+fac*(dry*dry+drx*drx)/rr2;
          double fryy=facc*(dry*dry)/(rr2)+fac*(drz*drz+drx*drx)/rr2;
          double frxz=facc*(drx*drz)/(rr2)-fac*(drx*drz)/rr2;
          double frxy=facc*(drx*dry)/(rr2)-fac*(drx*dry)/rr2;
          double fryz=facc*(dry*drz)/(rr2)-fac*(dry*drz)/rr2;

    			
    			const tdouble4 velrhopp2=TDouble4(velrhop[p2].x,velrhop[p2].y,velrhop[p2].z,velrhop[p2].w);
    			//===== Get mass and volume of particle p2 =====
    			double massp2=csp.massfluid;
    			double volp2=massp2/velrhopp2.w;
           // Set up C and C1 based on order and sim2d
          if(vrmethod==VrMethod_Mls){
            drx=drx/csp.kernelh;
            dry=dry/csp.kernelh;
            drz=drz/csp.kernelh;
          }

          // Set up C and C1 based on order and sim2d
          if (vrorder==VrOrder_2nd){
            if (sim2d){
              double tempC[] ={(1.0),drx,drz,drx*drx*(0.5),drx*drz,drz*drz*(0.5)};
              double tempC1[]={wab,frx,frz,frxx,frxz,frzz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }else{
              double tempC[]   ={(1.0),drx,dry,drz,drx*drx*(0.5),drx*dry,dry*dry*(0.5),drx*drz,drz*drz*(0.5),dry*drz};
              double tempC1[]  ={wab,frx,fry,frz,frxx,frxy,fryy,frxz,frzz,fryz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }            
          }else{
            if (sim2d){
              double tempC[]   = {(1.0),-drx,-drz};
              double tempC1[]  = {wab,frx,frz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }else{
              double tempC[] = {(1.0),drx,dry,drz};
              double tempC1[] = {wab,frx,fry,frz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }
          }

          // Set up D
          if (sim2d){
            double tempD[] = {velrhop[p2].w,velrhop[p2].x,velrhop[p2].z};
            for (unsigned i =0;i<m_dim2;++i)D[i]=tempD[i];
          }else{
            double tempD[] = {velrhop[p2].w,velrhop[p2].x,velrhop[p2].y,velrhop[p2].z};
            for (unsigned i =0;i<m_dim2;++i)D[i]=tempD[i];
          }
          // Accumulate A and B matrices
          for (unsigned i=0;i<m_dim;i++)
            for (unsigned j=0;j<m_dim;j++){
              A[i*m_dim+j]+=C1[i]*C[j];
            }
              
          for (unsigned i=0; i<m_dim2;i++)
            for (unsigned j=0;j<m_dim;j++){
              B[i*m_dim+j]+=D[i]*C1[j];
          }
        }
      }
    }

    double shep = A[0];                                   ///<Shepard summation.
    double sol[m_dim2]{0};                                ///<Solution array;
    double treshold = 0.0;                                ///<condition number;
    double cond = 0.0;                                    ///<Scaled reciprocal condition number;
    double kernelh2=csp.kernelh*csp.kernelh;              ///<Scaling factor;
    double scaleh= 0.0;
      
    if(shep>0.05){
      if (vrorder==VrOrder_2nd){
        if(vrmethod==VrMethod_Liu) scaleh = kernelh2*kernelh2;
        else scaleh=1.0;

        int P[m_dim]{0};

        LUdecomp_Single<m_dim,m_dim2>(A,P,B,sol,treshold);

        cond=(1.0/treshold)*scaleh;
      }else if (vrorder==VrOrder_1st){
        if(vrmethod==VrMethod_Liu) scaleh = kernelh2;
        else scaleh=1.0;

        int P[m_dim]{0};

        LUdecomp_Single<m_dim,m_dim2>(A,P,B,sol,treshold);
        cond=(1.0/treshold)*scaleh;
      }
      if (cond>mrthreshold || vrorder==VrOrder_0th)
      {
        for (unsigned i=0;i<m_dim2;i++)
          sol[i]=B[i*m_dim]/shep;
      }

      if (sim2d)
      {      
        velrhopb[p1].w = static_cast<float>(sol[0]);
        velrhopb[p1].x = static_cast<float>(sol[1]);
        velrhopb[p1].z = static_cast<float>(sol[2]);
      } else {
        velrhopb[p1].w = static_cast<float>(sol[0]);
        velrhopb[p1].x = static_cast<float>(sol[1]);
        velrhopb[p1].y = static_cast<float>(sol[2]);
        velrhopb[p1].z = static_cast<float>(sol[3]);
      }
         
    }else{
      codeb[p1] = CODE_SetOutIgnore(codeb[p1]);
    }
    
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker,bool sim2d,TpVresOrder vrorder> 
void Interaction_BufferExtrapT(unsigned bufferpartcount,const int *bufferpart
  ,const stinterparmscb &t,tdouble3 *posb,tfloat4 *velrhopb,typecode *codeb
  ,const TpVresMethod vrmethod,float mrthreshold)
{	
  if(vrmethod==VrMethod_Liu) InteractionBufferExtrap<tker,sim2d,vrorder,VrMethod_Liu> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb,mrthreshold);
  else InteractionBufferExtrap<tker,sim2d,vrorder,VrMethod_Mls> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb,mrthreshold);
}

template<TpKernel tker> 
void Interaction_BufferExtrap_ct0(unsigned bufferpartcount,const int *bufferpart
  ,const stinterparmscb &t,tdouble3 *posb,tfloat4 *velrhopb,typecode *codeb
  ,const TpVresOrder vrorder,const TpVresMethod vrmethod,float mrthreshold)
{
	if(t.csp.simulate2d){
    switch(vrorder){
      case VrOrder_0th: Interaction_BufferExtrapT<tker,true,VrOrder_0th>(bufferpartcount,bufferpart,t,posb,velrhopb,codeb,vrmethod,mrthreshold); break;
      case VrOrder_1st: Interaction_BufferExtrapT<tker,true,VrOrder_1st>(bufferpartcount,bufferpart,t,posb,velrhopb,codeb,vrmethod,mrthreshold); break;
      case VrOrder_2nd: Interaction_BufferExtrapT<tker,true,VrOrder_2nd>(bufferpartcount,bufferpart,t,posb,velrhopb,codeb,vrmethod,mrthreshold); break;
    }
  }else{
    switch(vrorder){
      case VrOrder_0th: Interaction_BufferExtrapT<tker,false,VrOrder_0th>(bufferpartcount,bufferpart,t,posb,velrhopb,codeb,vrmethod,mrthreshold); break;
      case VrOrder_1st: Interaction_BufferExtrapT<tker,false,VrOrder_1st>(bufferpartcount,bufferpart,t,posb,velrhopb,codeb,vrmethod,mrthreshold); break;
      case VrOrder_2nd: Interaction_BufferExtrapT<tker,false,VrOrder_2nd>(bufferpartcount,bufferpart,t,posb,velrhopb,codeb,vrmethod,mrthreshold); break;
    }
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void Interaction_BufferExtrap(unsigned bufferpartcount,const int *bufferpart
  ,const stinterparmscb &t,tdouble3 *posb,tfloat4 *velrhopb,typecode *codeb
  ,const TpVresOrder order,const TpVresMethod vrmethod,float mrthreshold)
{
  switch(t.csp.tkernel){
    case KERNEL_Wendland:
      Interaction_BufferExtrap_ct0<KERNEL_Wendland> (bufferpartcount,bufferpart,t,posb,velrhopb,codeb,order,vrmethod,mrthreshold);  break;
    break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:
    	throw "KERNEL_Cubic is not available when variable resolution is active.";
    break;
#endif
    default: throw "Kernel unknown at Interaction_InOutExtrap().";
  }
}



template<TpKernel tker,bool sim2d,TpVresOrder vrorder,TpVresMethod vrmethod> 
void InteractionBufferExtrapFlux(const unsigned n,const int pini
		,StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos
		,const typecode *code,const unsigned *idp,const tfloat4 *velrhop
    ,const StCteSph csp,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velflux
    ,float *fluxes,double dp,double dt,float mrthreshold)
{
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<int(n);p++){

    const int p1=p+pini;
    tdouble3 pos_p1=ptpoints[p1];

    // Compute size of the reconstruction matrix and right hand sime.
    constexpr unsigned m_dim  = vrorder==VrOrder_2nd ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
    constexpr unsigned m_dim2 = sim2d ? 3: 4;

    double C[m_dim]{0};
    double C1[m_dim]{0};
    double D[m_dim2]{0};

    double A[m_dim*m_dim]{0};
    double B[m_dim*m_dim2]{0};


    double ShiftTFS=0;
    double mindist=1000000.0;
    double mindp=min(dp,csp.dp);

    const StNgSearch ngsb=nsearch::Init(pos_p1,true,dvd);
    for(int z=ngsb.zini;z<ngsb.zfin;z++)for(int y=ngsb.yini;y<ngsb.yfin;y++){
    	const tuint2 pif=nsearch::ParticleRange(y,z,ngsb,dvd);
    	//-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
    	//---------------------------------------------------------------------------------------------
    	for(unsigned p2=pif.x;p2<pif.y;p2++){
    		double drx=double(pos_p1.x-pos[p2].x);
    		double dry=double(pos_p1.y-pos[p2].y);
    		double drz=double(pos_p1.z-pos[p2].z);
    		const double rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO){//-Only with fluid particles but not inout particles.
    			//-Computes kernel.
    			float fac;
    			const double wab=fsph::GetKernel_WabFac<tker>(csp,float(rr2),fac);
    			double frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.

          double massp2=csp.massbound;
    			double volp2=massp2/csp.rhopzero;
          ShiftTFS-=volp2*(drx*frx+dry*fry+drz*frz);
        }
      }
    }

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(pos_p1,false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
    	const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
    	//-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
    	//---------------------------------------------------------------------------------------------
    	for(unsigned p2=pif.x;p2<pif.y;p2++){
    		double drx=double(pos_p1.x-pos[p2].x);
    		double dry=double(pos_p1.y-pos[p2].y);
    		double drz=double(pos_p1.z-pos[p2].z);
    		const double rr2=drx*drx+dry*dry+drz*drz;

    		if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) &&!CODE_IsFluidFixed(code[p2])){//-Only with fluid particles but not inout particles.
    			//-Computes kernel.
    			float fac=0.f;
    			float facc=0.f;
    			const double wab=fsph::GetKernel_WabFacFacc<tker>(csp,float(rr2),fac,facc);
    			double frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.
    			double frxx=facc*(drx*drx)/(rr2)+fac*(dry*dry+drz*drz)/rr2;
          double frzz=facc*(drz*drz)/(rr2)+fac*(dry*dry+drx*drx)/rr2;
          double fryy=facc*(dry*dry)/(rr2)+fac*(drz*drz+drx*drx)/rr2;
          double frxz=facc*(drx*drz)/(rr2)-fac*(drx*drz)/rr2;
          double frxy=facc*(drx*dry)/(rr2)-fac*(drx*dry)/rr2;
          double fryz=facc*(dry*drz)/(rr2)-fac*(dry*drz)/rr2;

    			
    			const tdouble4 velrhopp2=TDouble4(velrhop[p2].x,velrhop[p2].y,velrhop[p2].z,velrhop[p2].w);
    			//===== Get mass and volume of particle p2 =====
    			double massp2=csp.massfluid;
    			double volp2=massp2/velrhopp2.w;
          ShiftTFS-=volp2*(drx*frx+dry*fry+drz*frz);
          mindist=min(mindist,rr2);
          
          // Set up C and C1 based on order and sim2d
          if(vrmethod==VrMethod_Mls){
            drx=drx/csp.kernelh;
            dry=dry/csp.kernelh;
            drz=drz/csp.kernelh;
          }

          // Set up C and C1 based on order and sim2d
          if (vrorder==VrOrder_2nd){
            if (sim2d){
              double tempC[] ={(1.0),drx,drz,drx*drx*(0.5),drx*drz,drz*drz*(0.5)};
              double tempC1[]={wab,frx,frz,frxx,frxz,frzz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }else{
              double tempC[]   ={(1.0),drx,dry,drz,drx*drx*(0.5),drx*dry,dry*dry*(0.5),drx*drz,drz*drz*(0.5),dry*drz};
              double tempC1[]  ={wab,frx,fry,frz,frxx,frxy,fryy,frxz,frzz,fryz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }            
          }else{
            if (sim2d){
              double tempC[]   = {(1.0),-drx,-drz};
              double tempC1[]  = {wab,frx,frz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }else{
              double tempC[] = {(1.0),drx,dry,drz};
              double tempC1[] = {wab,frx,fry,frz};
              for (unsigned i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                if(vrmethod==VrMethod_Liu) C1[i] = volp2*tempC1[i];
                else C1[i] = volp2*wab*tempC[i];
              }
            }
          }

              // Set up D
          if (sim2d){
            double tempD[] = {velrhop[p2].w,velrhop[p2].x,velrhop[p2].z};
            for (unsigned i =0;i<m_dim2;++i)D[i]=tempD[i];
          }else{
            double tempD[] = {velrhop[p2].w,velrhop[p2].x,velrhop[p2].y,velrhop[p2].z};
            for (unsigned i =0;i<m_dim2;++i)D[i]=tempD[i];
          }
          // Accumulate A and B matrices
          for (unsigned i=0;i<m_dim;i++)
            for (unsigned j=0;j<m_dim;j++){
              A[i*m_dim+j]+=C1[i]*C[j];
            }
              
          for (unsigned i=0; i<m_dim2;i++)
            for (unsigned j=0;j<m_dim;j++){
              B[i*m_dim+j]+=D[i]*C1[j];
          }
        }
      }
    }

  
    double shep = A[0];                                   ///<Shepard summation.
    double sol[m_dim2]{0};                                ///<Solution array;
    double treshold = 0.0;                                ///<condition number;
    double cond = 0.0;                                    ///<Scaled reciprocal condition number;
    double kernelh2=csp.kernelh*csp.kernelh;              ///<Scaling factor;
    double scaleh= 0.0;
      
    if(shep>0.05){
      if (vrorder==VrOrder_2nd){
        if(vrmethod==VrMethod_Liu) scaleh = kernelh2*kernelh2;
        else scaleh=1.0;

        int P[m_dim]{0};

        LUdecomp_Single<m_dim,m_dim2>(A,P,B,sol,treshold);

        cond=(1.0/treshold)*scaleh;
      }else if (vrorder==VrOrder_1st){
        if(vrmethod==VrMethod_Liu) scaleh = kernelh2;
        else scaleh=1.0;

        int P[m_dim]{0};

        LUdecomp_Single<m_dim,m_dim2>(A,P,B,sol,treshold);
        cond=(1.0/treshold)*scaleh;
      }
      if (cond>mrthreshold || vrorder==VrOrder_0th)
      {
        for (unsigned i=0;i<m_dim2;i++)
          sol[i]=B[i*m_dim]/shep;
      }

      if (sim2d){
        if((ShiftTFS>1.5 || sqrt(mindist)<mindp) && fluxes[p1]<0.0)
        fluxes[p1]+=static_cast<float>(max(0.0,-sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt));
        else if((ShiftTFS>1.5 || sqrt(mindist)<mindp))
        fluxes[p1]+= static_cast<float>(-sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt);
      }else{
        if      ((ShiftTFS>2.75 || sqrt(mindist)<mindp)  &&fluxes[p1]<0.0)  fluxes[p1]+= static_cast<float>(max(0.0,-sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt));
        else if ((ShiftTFS>2.75 || sqrt(mindist)<mindp))                    fluxes[p1]+= static_cast<float>(-sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt);
          
      }
    } 
    
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker,bool sim2d,TpVresOrder vrorder> 
void Interaction_BufferExtrapFluxT(const stinterparmscb &t,StrDataVresCpu &vres
  ,double dp,double dt,const TpVresMethod vrmethod,float mrthreshold)
{
	if(vrmethod==VrMethod_Liu) InteractionBufferExtrapFlux<tker,sim2d,vrorder,VrMethod_Liu> (vres.ntot,vres.nini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,vres.points,vres.normals,vres.velmot,vres.mass,dp,dt,mrthreshold);
	else InteractionBufferExtrapFlux<tker,sim2d,vrorder,VrMethod_Mls>  (vres.ntot,vres.nini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,vres.points,vres.normals,vres.velmot,vres.mass,dp,dt,mrthreshold);
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker> 
void Interaction_BufferExtrapFlux_ct0(const stinterparmscb &t,StrDataVresCpu &vres
  ,double dp,double dt,const TpVresOrder vrorder,const TpVresMethod vrmethod,float mrthreshold)
{
	if(t.csp.simulate2d){
    switch(vrorder){
      case VrOrder_0th: Interaction_BufferExtrapFluxT<tker,true,VrOrder_0th>(t,vres,dp,dt,vrmethod,mrthreshold); break;
      case VrOrder_1st: Interaction_BufferExtrapFluxT<tker,true,VrOrder_1st>(t,vres,dp,dt,vrmethod,mrthreshold); break;
      case VrOrder_2nd: Interaction_BufferExtrapFluxT<tker,true,VrOrder_2nd>(t,vres,dp,dt,vrmethod,mrthreshold); break;
    }
  }else{
      switch(vrorder){
      case VrOrder_0th: Interaction_BufferExtrapFluxT<tker,false,VrOrder_0th>(t,vres,dp,dt,vrmethod,mrthreshold); break;
      case VrOrder_1st: Interaction_BufferExtrapFluxT<tker,false,VrOrder_1st>(t,vres,dp,dt,vrmethod,mrthreshold); break;
      case VrOrder_2nd: Interaction_BufferExtrapFluxT<tker,false,VrOrder_2nd>(t,vres,dp,dt,vrmethod,mrthreshold); break;
    }
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void Interaction_BufferExtrapFlux(const stinterparmscb &t,StrDataVresCpu &vres
  ,double dp,double dt,const TpVresOrder vrorder,const TpVresMethod vrmethod,float mrthreshold)
{
  switch(t.csp.tkernel){
    case KERNEL_Wendland:
      Interaction_BufferExtrapFlux_ct0<KERNEL_Wendland>(t,vres,dp,dt,vrorder,vrmethod,mrthreshold);
    break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:
    	throw "KERNEL_Cubic is not available when variable resolution is active.";
    break;
#endif
    default: throw "Kernel unknown at Interaction_InOutExtrap().";
  }
}

tdouble3 MovePoint(tdouble3 oldpos,const tmatrix4d& mat){
  tdouble3 newpos=TDouble3(0);
  oldpos.x-=mat.a14;  oldpos.y-=mat.a24;  oldpos.z-=mat.a34;
  newpos.x=oldpos.x*mat.a11+oldpos.y*mat.a21+oldpos.z*mat.a31;
  newpos.y=oldpos.x*mat.a12+oldpos.y*mat.a22+oldpos.z*mat.a32;
  newpos.z=oldpos.x*mat.a13+oldpos.y*mat.a23+oldpos.z*mat.a33;
  return(newpos);
}

void BufferShiftingCpu(unsigned n,unsigned pinit,const tdouble3 *pos
  ,tfloat4 *shiftpos,const typecode *code,StrGeomVresCpu& vresdata)
{
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<int(n);p1++){
    typecode rcode=code[p1];
	  if(CODE_IsFluidBuffer(rcode)){

		  const byte izone0=byte(CODE_GetIzoneFluidBuffer(rcode));
		  const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
      tdouble3 boxmin=vresdata.boxlimitmin[izone];
      tdouble3 boxmax=vresdata.boxlimitmax[izone];
      tdouble3 origin=(boxmax+boxmin)*0.5;
      tdouble3 boxsize=(boxmax-boxmin);

      bool tracking=vresdata.tracking[izone];
      tmatrix4d mat=vresdata.matmov[izone].GetMatrix4d();
      tdouble3 ps=pos[p1];
      if(tracking)ps=MovePoint(ps,mat);
      tdouble3 dis=ps-origin;
      tfloat4 shiftp1=shiftpos[p1];     

      if((fabs(dis.x)>boxsize.x/2.0 )){
        tfloat3 normal=TFloat3(static_cast<float>(mat.a11), static_cast<float>(mat.a21), static_cast<float>(mat.a31));
        shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                // shiftpos[p1].x=0.0f;

      } 
      if((fabs(dis.y)>boxsize.y/2.0)){
        tfloat3 normal=TFloat3(static_cast<float>(mat.a12), static_cast<float>(mat.a22), static_cast<float>(mat.a32));
        shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                // shiftpos[p1].y=0.0f;

      }
      if((fabs(dis.z)>boxsize.z/2.0)){
        tfloat3 normal=TFloat3(static_cast<float>((mat.a13)), static_cast<float>(mat.a23), static_cast<float>(mat.a33));
        shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                        // shiftpos[p1].z=0.0f;
      }
    }
  }
}


//==============================================================================
/// Creates list with free-surface particle (normal and periodic).
//==============================================================================
unsigned CountFreeSurfaceParticles(unsigned npf,unsigned pini
  ,const unsigned* fstype,unsigned* listp)
{
  unsigned count=0;
  const unsigned pfin=pini+npf;
  for(unsigned p=pini;p<pfin;p++){
    const unsigned fstypep=fstype[p];
    if(fstypep){//-It includes normal and periodic particles.
      listp[count]=p; count++;
    }
  }
  return(count);
}

//==============================================================================
/// Correct mass accumulated and avoid generation inside boundaries.
//==============================================================================
void CheckMassFlux(unsigned n,unsigned pinit,const StCteSph csp
  ,const StDivDataCpu& divdata,const unsigned* dcell,const tdouble3* pos
  ,const typecode* code,tdouble3* posb,tfloat3 *normals,float *fluxes)
{
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<int(n);p1++){

    tdouble3 posp1=posb[p1];

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(posp1,true,divdata);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
      const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);
      for(unsigned p2=pif.x;p2<pif.y;p2++){
        const float drx=float(posp1.x-pos[p2].x);
        const float dry=float(posp1.y-pos[p2].y);
        const float drz=float(posp1.z-pos[p2].z);
        const float rr2=drx*drx+dry*dry+drz*drz;
        if(rr2<csp.dp*csp.dp){
          const tfloat3 normalp1=normals[p1];
          const float norm =((-normalp1.x*drx-normalp1.y*dry-normalp1.z*drz)/sqrt(rr2));
          if(acos(norm)>0.785398) fluxes[p1]=0; 
        }      
      }
    }
  }
}

//==============================================================================
/// Perform interaction between particles: Fluid/Float-Fluid/Float or Fluid/Float-Bound
/// Realiza interaccion entre particulas: Fluid/Float-Fluid/Float or Fluid/Float-Bound
//==============================================================================
template<TpKernel tker,bool sim2d> void InteractionComputeFSNormals
  (unsigned np,unsigned pinit,const StCteSph csp
  ,StDivDataCpu divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,const unsigned* listp,unsigned* fstype,tfloat3* fsnormal,StrGeomVresCpu& vresdata)
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const unsigned p1=listp[p];

    float fs_treshold=0;                                //-Divergence of the position.
    tfloat3 gradc=TFloat3(0);                           //-Gradient of the concentration
    unsigned neigh=0;                                   //-Number of neighbours.
      
    tmatrix3f lcorr;        lcorr=TMatrix3f(0);         //-Correction matrix.
    tmatrix3f lcorr_inv;    lcorr_inv=TMatrix3f(0);     //-Inverse of the correction matrix.

    float Nzero=0;
    float ks2=csp.kernelsize2;
    double dp=csp.dp;
    if(sim2d)Nzero=float((3.141592)*ks2/(dp*dp));
    else     Nzero=float((4.f/3.f)*(3.141592)*ks2*ks2/(dp*dp*dp));
    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];

    bool boundp2=false;

    for(int b2=0;b2<2;b2++){
      if(b2==1)boundp2=true;

      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

        //-Interaction of Fluid with type Fluid or Bound.
        //-----------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(posp1.x-pos[p2].x);
          const float dry=float(posp1.y-pos[p2].y);
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO){
          //-Computes kernel.
            const float fac=fsph::GetKernel_Fac<tker>(csp,rr2);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
            const float rhopp2= float(velrho[p2].w);

            //===== Get mass of particle p2 ===== 
            float massp2=(boundp2? csp.massbound: csp.massfluid); //-Contiene masa de particula segun sea bound o fluid.

            // bool ftp2;
            // float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
            // const typecode cod=code[p2];
            // ftp2=CODE_IsFloating(cod);
            // ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);

            const float vol2=(float(massp2/rhopp2));
            neigh++;

            const float dot3=drx*frx+dry*fry+drz*frz;
            gradc.x+=vol2*frx;
            gradc.y+=vol2*fry;
            gradc.z+=vol2*frz;

            fs_treshold-=vol2*dot3;
            lcorr.a11+=-drx*frx*vol2; lcorr.a12+=-drx*fry*vol2; lcorr.a13+=-drx*frz*vol2;
            lcorr.a21+=-dry*frx*vol2; lcorr.a22+=-dry*fry*vol2; lcorr.a23+=-dry*frz*vol2;
            lcorr.a31+=-drz*frx*vol2; lcorr.a32+=-drz*fry*vol2; lcorr.a33+=-drz*frz*vol2;

            // const float wab=cufsph::GetKernel_Wab<KERNEL_Wendland>(rr2);
            // pou+=wab*vol2;  
          }
        }
      }
    }

    if(CODE_IsFluidBuffer(code[p1])){
    //-Loop through the local virtual stencil.
      double ks=csp.kernelsize;
      double dp=csp.dp;
      const byte izone0=byte(CODE_GetIzoneFluidBuffer(code[p1]));
      const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
      tdouble3 minpos=posp1-TDouble3(ks,sim2d? 0: ks,ks);
      tdouble3 maxpos=posp1+TDouble3(ks,sim2d? 0: ks,ks);
      for(double rx=minpos.x;rx<=maxpos.x; rx+=dp) for(double ry=minpos.y; ry<=maxpos.y; ry+=dp)
      for(double rz=minpos.z;rz<=maxpos.z; rz+=dp){
        const float drx=float(posp1.x-rx);
        const float dry=float(posp1.y-ry);
        const float drz=float(posp1.z-rz);
        const float rr2=drx*drx+dry*dry+drz*drz;          
        float rx1 = static_cast<float>(rx); float ry1 = static_cast<float>(ry); float rz1 = static_cast<float>(rz);
        if(vresdata.tracking[izone]){
          tmatrix4d mat=vresdata.matmov[izone].GetMatrix4d();
          rx1=static_cast<float>((rx-mat.a14)*mat.a11+(ry-mat.a24)*mat.a21+(rz-mat.a34)*mat.a31);
          ry1=static_cast<float>((rx-mat.a14)*mat.a12+(ry-mat.a24)*mat.a22+(rz-mat.a34)*mat.a32);
          rz1=static_cast<float>((rx-mat.a14)*mat.a13+(ry-mat.a24)*mat.a23+(rz-mat.a34)*mat.a33);
        }
        tdouble3 posp2=TDouble3(rx1,ry1,rz1);
        tdouble3 boxmin=vresdata.boxdommin[izone];
        tdouble3 boxmax=vresdata.boxdommax[izone];
        bool outside=(vresdata.inner[izone] ?!InZone(posp2,boxmin,boxmax): InZone(posp2,boxmin,boxmax));
        if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO && outside){
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(csp,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass of particle p2 ===== 
          float massp2=(boundp2? csp.massbound: csp.massfluid); //-Contiene masa de particula segun sea bound o fluid.
          const float vol2=(float(massp2/csp.rhopzero));
          neigh++;

          const float dot3=drx*frx+dry*fry+drz*frz;
          gradc.x+=vol2*frx;
          gradc.y+=vol2*fry;
          gradc.z+=vol2*frz;

          fs_treshold-=vol2*dot3;
          lcorr.a11+=-drx*frx*vol2; lcorr.a12+=-drx*fry*vol2; lcorr.a13+=-drx*frz*vol2;
          lcorr.a21+=-dry*frx*vol2; lcorr.a22+=-dry*fry*vol2; lcorr.a23+=-dry*frz*vol2;
          lcorr.a31+=-drz*frx*vol2; lcorr.a32+=-drz*fry*vol2; lcorr.a33+=-drz*frz*vol2;    
        }
      }
    }
    
    //-Find particles that are probably on the free-surface.
    unsigned fstypep1=0;
    if(sim2d){
      if(fs_treshold<1.7)fstypep1=2;
      if(fs_treshold<1.1 && Nzero/float(neigh)<0.4f)fstypep1=3;
    }
    else{
      if(fs_treshold<2.75)fstypep1=2;
      if(fs_treshold<1.8 && Nzero/float(neigh)<0.4f)fstypep1=3;
    }
    fstype[p1]=fstypep1;
    
    //-Calculation of the inverse of the correction matrix (Don't think there is a better way, create function for Determinant2x2 for clarity?).
    if(sim2d){
      tmatrix2f lcorr2d;
      tmatrix2f lcorr2d_inv;
      lcorr2d.a11=lcorr.a11; lcorr2d.a12=lcorr.a13;
      lcorr2d.a21=lcorr.a31; lcorr2d.a22=lcorr.a33;
      float lcorr_det=(lcorr2d.a11*lcorr2d.a22-lcorr2d.a12*lcorr2d.a21);
      lcorr2d_inv.a11=lcorr2d.a22/lcorr_det; lcorr2d_inv.a12=-lcorr2d.a12/lcorr_det; lcorr2d_inv.a22=lcorr2d.a11/lcorr_det; lcorr2d_inv.a21=-lcorr2d.a21/lcorr_det;
      lcorr_inv.a11=lcorr2d_inv.a11;  lcorr_inv.a13=lcorr2d_inv.a12;
      lcorr_inv.a31=lcorr2d_inv.a21;  lcorr_inv.a33=lcorr2d_inv.a22;
    }
    else{
      const float determ=fmath::Determinant3x3(lcorr);
      lcorr_inv=fmath::InverseMatrix3x3(lcorr,determ);
    }

    //-Correction of the gradient of concentration.
    tfloat3 gradc1=TFloat3(0);    
    gradc1.x=gradc.x*lcorr_inv.a11+gradc.y*lcorr_inv.a12+gradc.z*lcorr_inv.a13;
    gradc1.y=gradc.x*lcorr_inv.a21+gradc.y*lcorr_inv.a22+gradc.z*lcorr_inv.a23;
    gradc1.z=gradc.x*lcorr_inv.a31+gradc.y*lcorr_inv.a32+gradc.z*lcorr_inv.a33;    
    float gradc_norm=sqrt(gradc1.x*gradc1.x+gradc1.y*gradc1.y+gradc1.z*gradc1.z);
    //-Set normal.
    fsnormal[p1].x=-gradc1.x/gradc_norm;
    fsnormal[p1].y=-gradc1.y/gradc_norm;
    fsnormal[p1].z=-gradc1.z/gradc_norm;
  }
}

//==============================================================================
/// Computes free-surface particles and their normals.
//==============================================================================
void CallComputeFSNormals(const unsigned np,const unsigned npb
  ,const StCteSph csp,const StDivDataCpu& divdata,const unsigned* dcell
  ,const tdouble3* pos,const typecode* code,const tfloat4* velrho
  ,unsigned* fstype,tfloat3* fsnormal,unsigned* listp,StrGeomVresCpu& vresdata)
{
  const unsigned npf=np-npb;
  if(npf){
    //-Creates list with free-surface particle (normal and periodic).
    const unsigned count=CountFreeSurfaceParticles(npf,npb,fstype,listp);
    //-Computes normals on selected free-surface particles.
    if(count){
      if(csp.simulate2d){ const bool sim2d=true;
        switch(csp.tkernel){
          case KERNEL_Wendland:  InteractionComputeFSNormals<KERNEL_Wendland,sim2d>(count,npb,csp,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal,vresdata);  break;
          case KERNEL_Cubic:     InteractionComputeFSNormals<KERNEL_Cubic   ,sim2d>(count,npb,csp,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal,vresdata);  break;
          default: fun::Run_ExceptioonFun("Kernel unknown.");
        }
      }
      else{ const bool sim2d=false;
        switch(csp.tkernel){
          case KERNEL_Wendland:  InteractionComputeFSNormals<KERNEL_Wendland,sim2d>(count,npb,csp,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal,vresdata);  break;
          case KERNEL_Cubic:     InteractionComputeFSNormals<KERNEL_Cubic   ,sim2d>(count,npb,csp,divdata,dcell,pos,code,velrho,listp,fstype,fsnormal,vresdata);  break;
          default: fun::Run_ExceptioonFun("Kernel unknown.");
        }
      }
    }
  }
}

//==============================================================================
/// Interaction of Fluid-Fluid/Bound & Bound-Fluid.
/// Interaccion Fluid-Fluid/Bound & Bound-Fluid.
//==============================================================================
void InteractionCallScanUmbrellaRegion(unsigned np,unsigned pinit,const StCteSph csp
  ,StDivDataCpu divdata,const unsigned* dcell,const tdouble3* pos
  ,const typecode* code,const tfloat3* fsnormal,const unsigned* listp
  ,unsigned* fstype,StrGeomVresCpu& vresdata)
{
  //-Initialise execution with OpenMP. | Inicia ejecucion con OpenMP.
  const int n=int(np);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p=0;p<n;p++){
    const unsigned p1=listp[p];
    bool fs_flag=false;
    bool boundp2=false;

    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];

    for(int b2=0; b2<2;b2++){
      if(b2==1) boundp2=true;
      //-Search for neighbours in adjacent cells.
      const StNgSearch ngs=nsearch::Init(dcell[p1],boundp2,divdata);
      for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
        const tuint2 pif=nsearch::ParticleRange(y,z,ngs,divdata);

        //-Interaction of Fluid with type Fluid or Bound. | Interaccion de Fluid con varias Fluid o Bound.
        //------------------------------------------------------------------------------------------------
        for(unsigned p2=pif.x;p2<pif.y;p2++){
          const float drx=float(posp1.x-pos[p2].x);
          const float dry=float(posp1.y-pos[p2].y);
          const float drz=float(posp1.z-pos[p2].z);
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO){
            //-Computes kernel.
            const tfloat3 posq=fsnormal[p1]*csp.kernelh;
            if(rr2>2.f*csp.kernelh*csp.kernelh){
              const float drxq=-drx-posq.x;
              const float dryq=-dry-posq.y;
              const float drzq=-drz-posq.z;
              const float rrq=sqrt(drxq*drxq+dryq*dryq+drzq*drzq);
              if(rrq<csp.kernelh)fs_flag=true;
            }
            else{
              if(csp.simulate2d){
                const float drxq=-drx-posq.x;
                const float drzq=-drz-posq.z;
                const tfloat3 normalq=TFloat3(drxq*fsnormal[p1].x,0,drzq*fsnormal[p1].z);
                const tfloat3 tangq=TFloat3(-drxq*fsnormal[p1].z,0,drzq*fsnormal[p1].x);
                const float normalqnorm=sqrt(normalq.x*normalq.x+normalq.z*normalq.z);
                const float tangqnorm=sqrt(tangq.x*tangq.x+tangq.z*tangq.z);
                if(normalqnorm+tangqnorm<csp.kernelh)fs_flag=true;
              }
              else{
                float rrr=1.f/sqrt(rr2);
                const float arccosin=acos((-drx*fsnormal[p1].x*rrr-dry*fsnormal[p1].y*rrr-drz*fsnormal[p1].z*rrr));
                if(arccosin<0.785398f)fs_flag=true;
              }
            }
          }
        }
      }
    }


    if(CODE_IsFluidBuffer(code[p1])){
    //-Loop trough the local virtual stencil.
      double ks=csp.kernelsize;
      double dp=csp.dp;
      const byte izone0=byte(CODE_GetIzoneFluidBuffer(code[p1]));
      const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
      tdouble3 minpos=posp1-TDouble3(ks,csp.simulate2d? 0: ks,ks);
      tdouble3 maxpos=posp1+TDouble3(ks,csp.simulate2d? 0: ks,ks);
      for(double rx=minpos.x;rx<=maxpos.x; rx+=dp) for(double ry=minpos.y; ry<=maxpos.y; ry+=dp)
      for(double rz=minpos.z;rz<=maxpos.z; rz+=dp){
        const float drx=float(posp1.x-rx);
        const float dry=float(posp1.y-ry);
        const float drz=float(posp1.z-rz);
        const float rr2=drx*drx+dry*dry+drz*drz;          
        float rx1= static_cast<float>(rx); float ry1= static_cast<float>(ry); float rz1= static_cast<float>(rz);
        if(vresdata.tracking[izone]){
          tmatrix4d mat=vresdata.matmov[izone].GetMatrix4d();
          rx1=static_cast<float>((rx-mat.a14)*mat.a11+(ry-mat.a24)*mat.a21+(rz-mat.a34)*mat.a31);
          ry1=static_cast<float>((rx-mat.a14)*mat.a12+(ry-mat.a24)*mat.a22+(rz-mat.a34)*mat.a32);
          rz1=static_cast<float>((rx-mat.a14)*mat.a13+(ry-mat.a24)*mat.a23+(rz-mat.a34)*mat.a33);
        }
        tdouble3 posp2=TDouble3(rx1,ry1,rz1);
        tdouble3 boxmin=vresdata.boxdommin[izone];
        tdouble3 boxmax=vresdata.boxdommax[izone];
        bool outside=(vresdata.inner[izone] ?!InZone(posp2,boxmin,boxmax): InZone(posp2,boxmin,boxmax));
        if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO && outside){
          const tfloat3 posq=fsnormal[p1]*csp.kernelh;
          if(rr2>2.f*csp.kernelh*csp.kernelh){
            const float drxq=-drx-posq.x;
            const float dryq=-dry-posq.y;
            const float drzq=-drz-posq.z;
            const float rrq=sqrt(drxq*drxq+dryq*dryq+drzq*drzq);
            if(rrq<csp.kernelh)fs_flag=true;
          }
          else{
            if(csp.simulate2d){
              const float drxq=-drx-posq.x;
              const float drzq=-drz-posq.z;
              const tfloat3 normalq=TFloat3(drxq*fsnormal[p1].x,0,drzq*fsnormal[p1].z);
              const tfloat3 tangq=TFloat3(-drxq*fsnormal[p1].z,0,drzq*fsnormal[p1].x);
              const float normalqnorm=sqrt(normalq.x*normalq.x+normalq.z*normalq.z);
              const float tangqnorm=sqrt(tangq.x*tangq.x+tangq.z*tangq.z);
              if(normalqnorm+tangqnorm<csp.kernelh)fs_flag=true;
            }else{
              float rrr=1.f/sqrt(rr2);
              const float arccosin=acos((-drx*fsnormal[p1].x*rrr-dry*fsnormal[p1].y*rrr-drz*fsnormal[p1].z*rrr));
              if(arccosin<0.785398f)fs_flag=true;
            }
          }    
        }
      }
    }
    
    //-If particle was present in umbrella region, change the code of the particle.
    if(fs_flag && fstype[p1]==2)fstype[p1]=0;
    //-Periodic particle are internal by default.
    if(CODE_IsPeriodic(code[p1]))fstype[p1]=0;
  }
}

//==============================================================================
/// Scan Umbrella region to identify free-surface particle.
//==============================================================================
void CallScanUmbrellaRegion(const unsigned np,const unsigned npb
  ,const StCteSph csp,const StDivDataCpu& divdata
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,const tfloat3* fsnormal,unsigned* listp,unsigned* fstype,StrGeomVresCpu& vresdata)
{
  const unsigned npf=np-npb;
  if(npf){
    //-Obtain the list of particle that are probably on the free-surface (in ComputeUmbrellaRegion maybe is unnecessary).
    const unsigned count=CountFreeSurfaceParticles(npf,npb,fstype,listp);
    //-Scan Umbrella region on selected free-surface particles.
    if(count)InteractionCallScanUmbrellaRegion(count,npb,csp,divdata,dcell,pos,code,fsnormal,listp,fstype,vresdata);
  }
}


template<TpKernel tker,bool sim2d> void CorrectShiftBuff
  (const unsigned n,const unsigned pinit  ,const StCteSph csp
  ,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,tfloat4* shiftvel,unsigned* fstype,StrGeomVresCpu& vresdata){
  
  const int pfin=int(pinit+n);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif
  for(int p1=int(pinit);p1<pfin;p1++){
    
    //-Obtain data of particle p1.
    const tdouble3 posp1=pos[p1];
    tfloat4 shiftvelp1=shiftvel[p1];

    if(CODE_IsFluidBuffer(code[p1])){
    //-Loop trough the local virtual stencil.
      double ks=csp.kernelsize;
      double dp=csp.dp;
      const byte izone0=byte(CODE_GetIzoneFluidBuffer(code[p1]));
      const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
      tdouble3 minpos=posp1-TDouble3(ks,sim2d? 0: ks,ks)-TDouble3(dp/2,sim2d? 0:dp/2,dp/2);
      tdouble3 maxpos=posp1+TDouble3(ks,sim2d? 0: ks,ks)-TDouble3(dp/2,sim2d? 0:dp/2,dp/2);
      for(double rx=minpos.x;rx<=maxpos.x; rx+=dp) for(double ry=minpos.y; ry<=maxpos.y; ry+=dp)
      for(double rz=minpos.z;rz<=maxpos.z; rz+=dp){
        const float drx=float(posp1.x-rx);
        const float dry=float(posp1.y-ry);
        const float drz=float(posp1.z-rz);
        const float rr2=drx*drx+dry*dry+drz*drz;          
        float rx1 = static_cast<float>(rx); float ry1 = static_cast<float>(ry); float rz1 = static_cast<float>(rz);
        if(vresdata.tracking[izone]){
          tmatrix4d mat=vresdata.matmov[izone].GetMatrix4d();
          rx1=static_cast<float>((rx-mat.a14)*mat.a11+(ry-mat.a24)*mat.a21+(rz-mat.a34)*mat.a31);
          ry1=static_cast<float>((rx-mat.a14)*mat.a12+(ry-mat.a24)*mat.a22+(rz-mat.a34)*mat.a32);
          rz1=static_cast<float>((rx-mat.a14)*mat.a13+(ry-mat.a24)*mat.a23+(rz-mat.a34)*mat.a33);     
        }
        tdouble3 posp2=TDouble3(rx1,ry1,rz1);
        tdouble3 boxmin=vresdata.boxdommin[izone];
        tdouble3 boxmax=vresdata.boxdommax[izone];
        bool outside=(vresdata.inner[izone] ?!InZone(posp2,boxmin,boxmax): InZone(posp2,boxmin,boxmax));
        if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO && outside){
          //-Computes kernel.
          const float fac=fsph::GetKernel_Fac<tker>(csp,rr2);
          const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.

          //===== Get mass of particle p2 ===== 
          float massp2=csp.massfluid; //-Contiene masa de particula segun sea bound o fluid.
          const float massrho=(float(massp2/csp.rhopzero));

          //-Compute gradient of concentration and partition of unity.                  
          shiftvelp1.x+=massrho*frx;    
          shiftvelp1.y+=massrho*fry;
          shiftvelp1.z+=massrho*frz;          
          const float wab=fsph::GetKernel_Wab<tker>(csp,rr2);
          shiftvelp1.w+=wab*massrho;    
        }
      }
      shiftvel[p1]=shiftvelp1;
    }
  }
}


//==============================================================================
/// Computes free-surface particles and their normals.
//==============================================================================
void CallCorrectShiftBuff(const unsigned np,const unsigned npb
  ,const StCteSph csp,const unsigned* dcell,const tdouble3* pos,const typecode* code
  ,tfloat4* shiftvel,unsigned* fstype,StrGeomVresCpu& vresdata)
{
  const unsigned npf=np-npb;
  if(npf){
      if(csp.simulate2d){ const bool sim2d=true;
        switch(csp.tkernel){
          case KERNEL_Wendland:  CorrectShiftBuff<KERNEL_Wendland,sim2d>(npf,npb,csp,dcell,pos,code,shiftvel,fstype,vresdata);  break;
          case KERNEL_Cubic:     CorrectShiftBuff<KERNEL_Cubic   ,sim2d>(npf,npb,csp,dcell,pos,code,shiftvel,fstype,vresdata);  break;
          default: fun::Run_ExceptioonFun("Kernel unknown.");
        }
      }
      else{ const bool sim2d=false;
        switch(csp.tkernel){
          case KERNEL_Wendland:  CorrectShiftBuff<KERNEL_Wendland,sim2d>(npf,npb,csp,dcell,pos,code,shiftvel,fstype,vresdata);  break;
          case KERNEL_Cubic:     CorrectShiftBuff<KERNEL_Cubic   ,sim2d>(npf,npb,csp,dcell,pos,code,shiftvel,fstype,vresdata);  break;
          default: fun::Run_ExceptioonFun("Kernel unknown.");
        }
      }
  }
}
}
