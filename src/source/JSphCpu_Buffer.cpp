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
#include <climits>
#include <cstring>
#include <algorithm>



// using namespace std;

template<const int n,const int n1> bool LUdecomp(double *a,int *p,double *b,double *sol,double& treshold){


//norm matrix a
	double maxs=0;
#pragma unroll
	for (int i=0;i<n;i++){
		double sum=0;
#pragma unroll
		for(int j=0;j<n;j++) sum+=fabs(a[n*j+i]);
		maxs=std::max(maxs,sum);
	}





	double maxv=0;
#pragma unroll
	for(int i=0;i<n;i++){
		p[i]=i;
	}

#pragma unroll
	for(int i=0;i<n;i++){
		//pivoting
		maxv=0;
		int imax=i;
#pragma unroll
		for(int k=i;k<n;k++){
			if(std::abs(a[i+k*n])>maxv){
				maxv=std::abs(a[i+k*n]);
				imax=k;
			}
		}

		if(imax!=i){
			int j=p[i];
			p[i]=p[imax];
			p[imax]=j;
		}
#pragma unroll
		for(int j=0;j<n;j++){
			double temp=a[i*n+j];
			a[i*n+j]=a[imax*n+j];
			a[imax*n+j]=temp;
		}

		//LU decomp
#pragma unroll
		for (int j = i + 1; j < n; j++) {
			a[j*n+i] /= a[i*n+i];
#pragma unroll
			for (int k = i + 1; k < n; k++)
				a[j*n+k] -= a[j*n+i] * a[i*n+k];
		}

	}
	double ia[n*n]{0};
	//matrix inversion
#pragma unroll
	for (int j = 0; j < n; j++) {
#pragma unroll
		for (int i = 0; i < n; i++) {
			ia[i*n+j] = p[i] == j ? 1.0 : 0.0;
#pragma unroll
			for (int k = 0; k < i; k++)
				ia[i*n+j] -= a[i*n+k] * ia[k*n+j];
		}
#pragma unroll
		for (int i = n - 1; i >= 0; i--) {
#pragma unroll
			for (int k = i + 1; k < n; k++)
				ia[i*n+j] -= a[i*n+k] * ia[k*n+j];

			ia[i*n+j] /= a[i*n+i];
		}
	}

//norm of inv matrix
	double maxs1=0;
#pragma unroll
	for (int i=0;i<n;i++){
		double sum=0;
#pragma unroll
		for(int j=0;j<n;j++) sum+=std::abs(ia[n*j+i]);
		maxs1=std::max(maxs1,sum);
	}
	
  treshold=1.0/(maxs*maxs1);



	//solution
#pragma unroll
	for (int k = 0; k < n1; k++)
#pragma unroll
	for (int i = 0; i < n; i++) {
		sol[k] += ia[i]*b[i+k*n];
	}

	return true;

}

template<bool sim2d,TpKernel tker,unsigned order> void JSphCpu::InteractionBufferExtrap(unsigned bufferpartcount,const int *bufferpart,
		StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,
		const typecode *code,const unsigned *idp,const tfloat4 *velrhop,const StCteSph csp,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb)
{
    float scaleh=csp.kernelh*csp.kernelh*csp.kernelh*csp.kernelh;

  //Log->Printf("%u>++> InteractionInOutGhost_Double",Nstep);
  //-Inicia ejecucion con OpenMP.
  const int n=int(bufferpartcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif

  for(int p=0;p<n;p++){

    const int p1=bufferpart[p];
    tdouble3 pos_p1=posb[p1];

    constexpr unsigned nmatrix = (order == 0) ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
    constexpr unsigned nrhs = sim2d ? 3 : 4;

    double C[nmatrix]{0};
    double C1[nmatrix]{0};
    double D[nrhs]{0};


    double A[nmatrix*nmatrix]{0};
    double B[nmatrix*nrhs]{0};

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(pos_p1,false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
    	const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
    	//-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
    	//---------------------------------------------------------------------------------------------
    	for(unsigned p2=pif.x;p2<pif.y;p2++){
    		const double drx=-double(pos_p1.x-pos[p2].x);
    		const double dry=-double(pos_p1.y-pos[p2].y);
    		const double drz=-double(pos_p1.z-pos[p2].z);
    		const double rr2=drx*drx+dry*dry+drz*drz;

    		if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) &&!CODE_IsFluidFixed(code[p2])){//-Only with fluid particles but not inout particles.
    			//-Computes kernel.
    			float fac;
    			float facc;
    			const double wab=fsph::GetKernel_WabFac<tker>(csp,(rr2),fac);
    			double frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.
    			double frxx=fac*(dry*dry+drz*drz)/rr2;
        	double frzz=fac*(dry*dry+drx*drx)/rr2;
        	double fryy=fac*(drz*drz+drx*drx)/rr2;
        	double frxz=-fac*(drx*drz)/rr2;
        	double frxy=-fac*(drx*dry)/rr2;
        	double fryz=-fac*(dry*drz)/rr2;

    			
    			const tdouble4 velrhopp2=TDouble4(velrhop[p2].x,velrhop[p2].y,velrhop[p2].z,velrhop[p2].w);
    			//===== Get mass and volume of particle p2 =====
    			double massp2=csp.massfluid;
    			double volp2=massp2/velrhopp2.w;
          if constexpr (order == 0) {
    			  if constexpr (sim2d) {
              double tempC[] = {1.0f, drx, drz, drx * drx * 0.5f, drx * drz, drz * drz * 0.5f};
              double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz, volp2 * frxx, volp2 * frxz, volp2 * frzz};

              #pragma unroll
              for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
              }
            } else {
                double tempC[] = {1.0f, drx, dry, drz, drx * drx * 0.5f, drx * dry, dry * dry * 0.5f, drx * drz, drz * drz * 0.5f, dry * drz};
                double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz, volp2 * frxx, volp2 * frxy, volp2 * fryy, volp2 * frxz, volp2 * frzz, volp2 * fryz};

                #pragma unroll
                for (int i = 0; i < nmatrix; ++i) {
                    C[i] = tempC[i];
                    C1[i] = tempC1[i];
                }
            }
          } else {
              if constexpr (sim2d) {
                  double tempC[] = {1.0f, drx, drz};
                  double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz};

                  #pragma unroll
                  for (int i = 0; i < nmatrix; ++i) {
                      C[i] = tempC[i];
                      C1[i] = tempC1[i];
                  }
              } else {
                  double tempC[] = {1.0f, drx, dry, drz};
                  double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz};

                  #pragma unroll
                  for (int i = 0; i < nmatrix; ++i) {
                      C[i] = tempC[i];
                      C1[i] = tempC1[i];
                  }
              }
          }

          if constexpr (sim2d) {
            double tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].z};
            #pragma unroll
            for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
          } else {
            double tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].y, velrhop[p2].z};
            #pragma unroll
            for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
          }

              // Your computation logic
          #pragma unroll
          for (int i = 0; i < nmatrix; i++) {
            #pragma unroll
            for (int j = 0; j < nmatrix; j++) {
              A[i * nmatrix + j] += C1[i] * C[j];
            }
          }
            
          #pragma unroll
          for (int i = 0; i < nrhs; i++) {
            #pragma unroll
            for (int j = 0; j < nmatrix; j++) {
              B[i * nmatrix + j] += D[i] * C1[j];
            }
          }
        }
      }
    }

  double mrthreshold=100;
    double shep = A[0];
    double sol[nrhs]{0};
    double treshold=0;
    double cond=0;
      
      if (shep > 0.1){
        if (order==0){

          int P[nmatrix]{0};
          
          LUdecomp<nmatrix, nrhs>(A, P, B, sol, treshold);
    
    
          cond=(1.0f/treshold)*scaleh;

        } else if (order==1){

          int P[nmatrix]{0};


        LUdecomp<nmatrix, nrhs>(A, P, B, sol, treshold);
        cond=(1.0f/treshold)*CSP.kernelh*CSP.kernelh;


      }
      if (cond>mrthreshold || order==2){
         for (unsigned i = 0; i < nrhs; i++)
           sol[i] = B[i * nmatrix] / shep;
         }

      if (sim2d)
      {
         
        velrhopb[p1].w = sol[0];
        velrhopb[p1].x = sol[1];
        velrhopb[p1].z = sol[2];
      } else {
        velrhopb[p1].w = sol[0];
        velrhopb[p1].x = sol[1];
        velrhopb[p1].y = sol[2];
        velrhopb[p1].z = sol[3];
      }
         
      } else{
        codeb[p1] = CODE_SetOutIgnore(codeb[p1]);

      }
    
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker> void JSphCpu::Interaction_BufferExtrapT(unsigned bufferpartcount,const int *bufferpart
        ,const stinterparmscb &t,tdouble3*posb,tfloat4 *velrhopb,typecode *codeb,unsigned order)
{
	if(Simulate2D){const bool sim2d=true;
    if(order==0) InteractionBufferExtrap<sim2d ,tker,0> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb);
    if(order==1) InteractionBufferExtrap<sim2d ,tker,1> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb);
    if(order==2) InteractionBufferExtrap<sim2d ,tker,2> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb);

	}
	else{ const bool sim2d=false;
    if(order==0) InteractionBufferExtrap<sim2d ,tker,0> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb);
    if(order==1) InteractionBufferExtrap<sim2d ,tker,1> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb);
    if(order==2) InteractionBufferExtrap<sim2d ,tker,2> (bufferpartcount,bufferpart,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,posb,velrhopb,codeb);

	}

}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void JSphCpu::Interaction_BufferExtrap(unsigned bufferpartcount,const int *bufferpart
                                       ,const stinterparmscb &t,tdouble3 *posb,tfloat4 *velrhopb,typecode *codeb,unsigned order)
{
  switch(TKernel){
    // case KERNEL_Cubic:       Interaction_BufferExtrapT<KERNEL_Cubic>     (bufferpartcount,bufferpart,t,posb,velrhopb,order);  break;
    case KERNEL_Wendland:    Interaction_BufferExtrapT<KERNEL_Wendland>  (bufferpartcount,bufferpart,t,posb,velrhopb,codeb,order);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}



template<bool sim2d,TpKernel tker,unsigned order> void JSphCpu::InteractionBufferExtrapFlux(const unsigned n,const int pini,
		StDivDataCpu dvd,const unsigned *dcell,const tdouble3 *pos,
		const typecode *code,const unsigned *idp,const tfloat4 *velrhop,const StCteSph csp,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velflux,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold)
{
    float scaleh=csp.kernelh*csp.kernelh*csp.kernelh*csp.kernelh;

  //Log->Printf("%u>++> InteractionInOutGhost_Double",Nstep);
  //-Inicia ejecucion con OpenMP.
  // const int n=int(bufferpartcount);
  #ifdef OMP_USE
    #pragma omp parallel for schedule (guided)
  #endif

  for(int p=0;p<n;p++){

    const int p1=p+pini;
    tdouble3 pos_p1=ptpoints[p1];

    constexpr unsigned nmatrix = (order == 0) ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
    constexpr unsigned nrhs = sim2d ? 3 : 4;

    double C[nmatrix]{0};
    double C1[nmatrix]{0};
    double D[nrhs]{0};


    double A[nmatrix*nmatrix]{0};
    double B[nmatrix*nrhs]{0};


    float ShiftTFS=0;
    double mindist=1000000.0;
    float mindp=std::min(dp,csp.dp);

    //-Search for neighbours in adjacent cells.
    const StNgSearch ngs=nsearch::Init(pos_p1,false,dvd);
    for(int z=ngs.zini;z<ngs.zfin;z++)for(int y=ngs.yini;y<ngs.yfin;y++){
    	const tuint2 pif=nsearch::ParticleRange(y,z,ngs,dvd);
    	//-Interaction of boundary with type Fluid/Float | Interaccion de Bound con varias Fluid/Float.
    	//---------------------------------------------------------------------------------------------
    	for(unsigned p2=pif.x;p2<pif.y;p2++){
    		const double drx=-double(pos_p1.x-pos[p2].x);
    		const double dry=-double(pos_p1.y-pos[p2].y);
    		const double drz=-double(pos_p1.z-pos[p2].z);
    		const double rr2=drx*drx+dry*dry+drz*drz;

    		if(rr2<=csp.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) &&!CODE_IsFluidFixed(code[p2])){//-Only with fluid particles but not inout particles.
    			//-Computes kernel.
    			float fac;
    			float facc;
    			const double wab=fsph::GetKernel_WabFac<tker>(csp,(rr2),fac);
    			double frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.
    			double frxx=fac*(dry*dry+drz*drz)/rr2;
        	double frzz=fac*(dry*dry+drx*drx)/rr2;
        	double fryy=fac*(drz*drz+drx*drx)/rr2;
        	double frxz=-fac*(drx*drz)/rr2;
        	double frxy=-fac*(drx*dry)/rr2;
        	double fryz=-fac*(dry*drz)/rr2;

    			
    			const tdouble4 velrhopp2=TDouble4(velrhop[p2].x,velrhop[p2].y,velrhop[p2].z,velrhop[p2].w);
    			//===== Get mass and volume of particle p2 =====
    			double massp2=csp.massfluid;
    			double volp2=massp2/velrhopp2.w;
          ShiftTFS-=volp2*(drx*frx+dry*fry+drz*frz);
          mindist=std::min(mindist,rr2);
          if constexpr (order == 0) {
    			  if constexpr (sim2d) {
              double tempC[] = {1.0f, drx, drz, drx * drx * 0.5f, drx * drz, drz * drz * 0.5f};
              double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz, volp2 * frxx, volp2 * frxz, volp2 * frzz};

              #pragma unroll
              for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
              }
            } else {
                double tempC[] = {1.0f, drx, dry, drz, drx * drx * 0.5f, drx * dry, dry * dry * 0.5f, drx * drz, drz * drz * 0.5f, dry * drz};
                double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz, volp2 * frxx, volp2 * frxy, volp2 * fryy, volp2 * frxz, volp2 * frzz, volp2 * fryz};

                #pragma unroll
                for (int i = 0; i < nmatrix; ++i) {
                    C[i] = tempC[i];
                    C1[i] = tempC1[i];
                }
            }
          } else {
              if constexpr (sim2d) {
                  double tempC[] = {1.0f, drx, drz};
                  double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz};

                  #pragma unroll
                  for (int i = 0; i < nmatrix; ++i) {
                      C[i] = tempC[i];
                      C1[i] = tempC1[i];
                  }
              } else {
                  double tempC[] = {1.0f, drx, dry, drz};
                  double tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz};

                  #pragma unroll
                  for (int i = 0; i < nmatrix; ++i) {
                      C[i] = tempC[i];
                      C1[i] = tempC1[i];
                  }
              }
          }

          if constexpr (sim2d) {
            double tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].z};
            #pragma unroll
            for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
          } else {
            double tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].y, velrhop[p2].z};
            #pragma unroll
            for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
          }

              // Your computation logic
          #pragma unroll
          for (int i = 0; i < nmatrix; i++) {
            #pragma unroll
            for (int j = 0; j < nmatrix; j++) {
              A[i * nmatrix + j] += C1[i] * C[j];
            }
          }
            
          #pragma unroll
          for (int i = 0; i < nrhs; i++) {
            #pragma unroll
            for (int j = 0; j < nmatrix; j++) {
              B[i * nmatrix + j] += D[i] * C1[j];
            }
          }
        }
      }
    }

  
    double shep = A[0];
    double sol[nrhs]{0};
    double treshold=0;
    double cond=0;
      
      if (shep > 0.1){
        if (order==0){

          int P[nmatrix]{0};
          
          LUdecomp<nmatrix, nrhs>(A, P, B, sol, treshold);
    
    
          cond=(1.0f/treshold)*scaleh;

        } else if (order==1){

          int P[nmatrix]{0};


        LUdecomp<nmatrix, nrhs>(A, P, B, sol, treshold);
        cond=(1.0f/treshold)*CSP.kernelh*CSP.kernelh;


      }
      if (cond>mrthreshold || order==2){
         for (unsigned i = 0; i < nrhs; i++)
           sol[i] = B[i * nmatrix] / shep;
         }

      if (sim2d)
      {
        if ((ShiftTFS > 1.5 || sqrt(mindist) < mindp) && fluxes[p1] < 0.0 )
        fluxes[p1] += std::max(0.0, -sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt);
        else if ((ShiftTFS > 1.5 || sqrt(mindist) < mindp) )
        fluxes[p1] += -sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt;
      } else {
        if((ShiftTFS>2.75 || sqrt(mindist)<mindp)  && fluxes[p1]<0.0 ) fluxes[p1]+=std::max(0.0,-sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt);
          else if((ShiftTFS>2.75 || sqrt(mindist)<mindp) )               fluxes[p1]+= -sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt;
          
      }
         
      } 
    
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker> void JSphCpu::Interaction_BufferExtrapFluxT(const unsigned n,const int pini
  ,const stinterparmscb &t,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold)
{
	if(Simulate2D){const bool sim2d=true;
    if(mrorder==0) InteractionBufferExtrapFlux<sim2d ,tker,0> (n,pini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);
    if(mrorder==1) InteractionBufferExtrapFlux<sim2d ,tker,1> (n,pini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);
    if(mrorder==2) InteractionBufferExtrapFlux<sim2d ,tker,2> (n,pini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);

	}
	else{ const bool sim2d=false;
    if(mrorder==0) InteractionBufferExtrapFlux<sim2d ,tker,0> (n,pini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);
    if(mrorder==1) InteractionBufferExtrapFlux<sim2d ,tker,1> (n,pini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);
    if(mrorder==2) InteractionBufferExtrapFlux<sim2d ,tker,2> (n,pini,t.divdata,t.dcell,t.pos,t.code,t.idp,t.velrhop,t.csp,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);

	}

}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void JSphCpu::Interaction_BufferExtrapFlux(const unsigned n,const int pini
  ,const stinterparmscb &t,tdouble3 *ptpoints,tfloat3 *normals,tfloat3* velmot,float *fluxes
  ,unsigned mrorder,double dp,double dt,float mrthreshold)
{
  switch(TKernel){
    // case KERNEL_Cubic:       Interaction_BufferExtrapT<KERNEL_Cubic>  (n,pini,t,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);  break;
    case KERNEL_Wendland:    Interaction_BufferExtrapFluxT<KERNEL_Wendland>  (n,pini,t,ptpoints,normals,velmot,fluxes,mrorder,dp,dt,mrthreshold);  break;
    default: Run_Exceptioon("Kernel unknown.");
  }
}
