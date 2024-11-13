/*
 * JSphGpu_Buffer_iker.cu
 *
 *  Created on: Apr 25, 2022
 *      Author: francesco
 */

#include "JSphGpu_Buffer_iker.h"
#include "Functions.h"
#include "FunctionsCuda.h"
#include <cfloat>
#include <math_constants.h>
#include "DualSphDef.h"
#include "JSphGpu_ker.h"
#include "JSphBuffer.h"
// #define T(id) (threadIdx.x +blockDim.x*(id))
#define T(id) (id)
namespace cusphbuffer{

#include "FunctionsBasic_iker.h"
#include "FunctionsMath_iker.h"
#include "FunSphKernel_iker.h"
#include "FunctionsGeo3d_iker.h"


#undef _JCellSearch_iker_
#include "JCellSearch_iker.h"

__device__ void MovePoint(double2 rxy,double rz,double2& rxy_t,double& rz_t,tmatrix4f mat){
  rxy_t.x=rxy.x*mat.a11+rxy.y*mat.a21+rz*mat.a31+mat.a14;
  rxy_t.y=rxy.x*mat.a12+rxy.y*mat.a22+rz*mat.a32+mat.a24;
  rz_t   =rxy.x*mat.a13+rxy.y*mat.a23+rz*mat.a33+mat.a34;
}

__device__ bool KerBufferInZone(double2 rxy,double rz,double3 boxlimitmin,double3 boxlimitmax)
{
  return(boxlimitmin.x<=rxy.x && rxy.x<=boxlimitmax.x && boxlimitmin.y<=rxy.y && rxy.y<=boxlimitmax.y && boxlimitmin.z<=rz && rz<=boxlimitmax.z);
}

//------------------------------------------------------------------------------
/// Creates list with current inout particles and normal (no periodic) fluid in
/// inlet/outlet zones (update its code).
//------------------------------------------------------------------------------
__global__ void KerBufferCreateList(unsigned n,unsigned pini,const double3 boxlimitmininner,const double3 boxlimitmaxinner
  ,const double3 boxlimitminouter,const double3 boxlimitmaxouter,const bool inner,const double2 *posxy,const double *posz
  ,typecode *code,unsigned *listp,tmatrix4f mat,bool tracking,unsigned nzone)
{
  extern __shared__ unsigned slist[];
  //float *splanes=(float*)(slist+(n+1));
  if(!threadIdx.x)slist[0]=0;
  __syncthreads();
  const unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    const typecode rcode=code[p];
    if(CODE_IsNormal(rcode) && CODE_IsFluid(rcode)){//-It includes only normal fluid particles (no periodic).
      bool select=CODE_IsFluidBuffer(rcode);//-Particles already selected as InOut.
      if(!select){//-Particulas no periodicas y no marcadas como in/out.
        double2 rxy=posxy[p];
        double rz=posz[p];
        if(tracking){
          double2 rxy_t=make_double2(0,0);
          double rz_t=0.0;
          MovePoint(rxy,rz,rxy_t,rz_t,mat);
          rxy=rxy_t; rz=rz_t;
        }
        byte zone=255;
        if(KerBufferInZone(rxy,rz,boxlimitminouter,boxlimitmaxouter)&&!KerBufferInZone(rxy,rz,boxlimitmininner,boxlimitmaxinner))zone=byte(nzone);
          if(zone!=255){
            code[p]=CODE_ToFluidBuffer(rcode,nzone)|CODE_TYPE_FLUID_BUFFERNUM; //-Adds 16 to indicate new particle in zone.
            select=true;
          }
      } else{
    	  const byte izone0=byte(CODE_GetIzoneFluidBuffer(rcode));
    	  		const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
    	  		if(izone !=byte(nzone)) select=false;
      }
      if(select)slist[atomicAdd(slist,1)+1]=p; //-Add particle in the list.
    }
  }
  __syncthreads();
  const unsigned ns=slist[0];
  __syncthreads();
  if(!threadIdx.x && ns)slist[0]=atomicAdd((listp+n),ns);
  __syncthreads();
  if(threadIdx.x<ns){
    const unsigned cp=slist[0]+threadIdx.x;
    listp[cp]=slist[threadIdx.x+1];
  }
}


//==============================================================================
/// Creates list with current inout particles and normal (no periodic) fluid in
/// inlet/outlet zones (update its code).
//==============================================================================
unsigned BufferCreateList(bool stable,unsigned n,unsigned pini,const tdouble3 boxlimitmininner,const tdouble3 boxlimitmaxinner,
          const tdouble3 boxlimitminouter,const tdouble3 boxlimitmaxouter,const bool inner,const double2 *posxy,const double *posz
	        ,typecode *code,unsigned *listp,tmatrix4f mat,bool tracking,unsigned nzone)
{
  unsigned count=0;
  if(n){
    //-listp size list initialized to zero.
    //-Inicializa tamanho de lista listp a cero.
    cudaMemset(listp+n,0,sizeof(unsigned));
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    const unsigned smem=(SPHBSIZE+1)*sizeof(unsigned); //-All fluid particles can be in in/out area and one position for counter.
    KerBufferCreateList <<<sgrid,SPHBSIZE,smem>>> (n,pini,Double3(boxlimitmininner),Double3(boxlimitmaxinner),Double3(boxlimitminouter)
    	      ,Double3(boxlimitmaxouter),inner,posxy,posz,code,listp,mat,tracking,nzone);
    cudaMemcpy(&count,listp+n,sizeof(unsigned),cudaMemcpyDeviceToHost);
   // KerTestCTE <<<sgrid,SPHBSIZE,smem>>>();
    //-Reorders list when stable has been activated.
    //-Reordena lista cuando stable esta activado.
//    if(stable && count){ //-Does not affect results.
//      thrust::device_ptr<unsigned> dev_list(listp);
//      thrust::sort(dev_list,dev_list+count);
//    }
  //  cudaDeviceSynchronize();
  }
  return(count);
}



////------------------------------------------------------------------------------
///// Creates list with current inout particles and normal (no periodic) fluid in
///// inlet/outlet zones (update its code).
////------------------------------------------------------------------------------
//__global__ void KerBufferCreateList1(unsigned n,unsigned pini
//  ,const double3 boxlimitmininner,const double3 boxlimitmaxinner,const double3 boxlimitminouter,const double3 boxlimitmaxouter,const bool inner
//  ,const double2 *posxy,const double *posz
//  ,typecode *code,unsigned *listp,unsigned nzone)
//{
//  extern __shared__ unsigned slist[];
//  //float *splanes=(float*)(slist+(n+1));
//  if(!threadIdx.x)slist[0]=0;
//  __syncthreads();
//  const unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
//  if(pp<n){
//    const unsigned p=pp+pini;
//    const typecode rcode=code[p];
//    if(CODE_IsNormal(rcode) && CODE_IsFluid(rcode)){//-It includes only normal fluid particles (no periodic).
//      bool select=CODE_IsFluidBuffer(rcode);//-Particles already selected as InOut.
//      if(select){
//        const byte izone0=byte(CODE_GetIzoneFluidBuffer(rcode));
//    	  		const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
//    	  		if(izone !=byte(nzone)) select=false;
//      }
//    	  
//      }
//      if(select)slist[atomicAdd(slist,1)+1]=p; //-Add particle in the list.
//    }
//  
//  __syncthreads();
//  const unsigned ns=slist[0];
//  __syncthreads();
//  if(!threadIdx.x && ns)slist[0]=atomicAdd((listp+n),ns);
//  __syncthreads();
//  if(threadIdx.x<ns){
//    const unsigned cp=slist[0]+threadIdx.x;
//    listp[cp]=slist[threadIdx.x+1];
//  }
//}


//------------------------------------------------------------------------------
/// Creates list with current inout particles and normal (no periodic) fluid in
/// inlet/outlet zones (update its code).
//------------------------------------------------------------------------------
__global__ void KerBufferCreateListInit(unsigned n,unsigned pini
  ,const double3 boxlimitmininner,const double3 boxlimitmaxinner,const double3 boxlimitminouter,const double3 boxlimitmaxouter,const double3 boxlimitminmid,const double3 boxlimitmaxmid,const bool inner
  ,const double2 *posxy,const double *posz
  ,typecode *code,unsigned *listp,tmatrix4f mat,bool tracking,unsigned nzone)
{
  extern __shared__ unsigned slist[];
  //float *splanes=(float*)(slist+(n+1));
  if(!threadIdx.x)slist[0]=0;
  __syncthreads();
  bool inner1=inner;
  const unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(pp<n){
    const unsigned p=pp+pini;
    const typecode rcode=code[p];
    if(CODE_IsNormal(rcode) ){//-It includes only normal fluid particles (no periodic).
      bool select=CODE_IsFluidBuffer(rcode);//-Particles already selected as InOut.
      if(!select){//-Particulas no periodicas y no marcadas como in/out.
        double2 rxy=posxy[p];
        double rz=posz[p];
        if(tracking){
          double2 rxy_t=make_double2(0,0);
          double rz_t=0.0;
          MovePoint(rxy,rz,rxy_t,rz_t,mat);
          rxy=rxy_t; rz=rz_t;
        }
        byte zone=255;
        if(KerBufferInZone(rxy,rz,boxlimitminouter,boxlimitmaxouter)&&!KerBufferInZone(rxy,rz,boxlimitmininner,boxlimitmaxinner)&& CODE_IsFluid(rcode))zone=byte(nzone);
          if(zone!=255){
            code[p]=CODE_ToFluidBuffer(rcode,nzone)|CODE_TYPE_FLUID_BUFFERNUM; //-Adds 16 to indicate new particle in zone.
            select=true;
          }
        if((!KerBufferInZone(rxy,rz,boxlimitminouter,boxlimitmaxouter)&& inner) ||(KerBufferInZone(rxy,rz,boxlimitmininner,boxlimitmaxinner)&&!inner)){
        select=false;
        if(!(KerBufferInZone(rxy,rz,boxlimitminmid,boxlimitmaxmid)!=inner1)&& CODE_IsFluid(rcode)) {
        	code[p]=CODE_SetOutIgnore(rcode);//CODE_ToFluidFixed(rcode,nzone);
        }
        else code[p]=CODE_SetOutIgnore(rcode);
        }

      }
      if(select)slist[atomicAdd(slist,1)+1]=p; //-Add particle in the list.
    }
  }
  __syncthreads();
  const unsigned ns=slist[0];
  __syncthreads();
  if(!threadIdx.x && ns)slist[0]=atomicAdd((listp+n),ns);
  __syncthreads();
  if(threadIdx.x<ns){
    const unsigned cp=slist[0]+threadIdx.x;
    listp[cp]=slist[threadIdx.x+1];
  }
}


//==============================================================================
/// Creates list with current inout particles and normal (no periodic) fluid in
/// inlet/outlet zones (update its code).
//==============================================================================
unsigned BufferCreateListInit(bool stable,unsigned n,unsigned pini
  ,const tdouble3 boxlimitmininner,const tdouble3 boxlimitmaxinner,const tdouble3 boxlimitminouter,const tdouble3 boxlimitmaxouter,const tdouble3 boxlimitminmid,const tdouble3 boxlimitmaxmid,const bool inner,const double2 *posxy,const double *posz
  ,typecode *code,unsigned *listp,tmatrix4f mat,bool tracking,unsigned nzone)
{
  unsigned count=0;
  if(n){
    //-listp size list initialized to zero.
    //-Inicializa tamanho de lista listp a cero.
    cudaMemset(listp+n,0,sizeof(unsigned));
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    const unsigned smem=(SPHBSIZE+1)*sizeof(unsigned); //-All fluid particles can be in in/out area and one position for counter.
    KerBufferCreateListInit <<<sgrid,SPHBSIZE,smem>>> (n,pini,Double3(boxlimitmininner),Double3(boxlimitmaxinner),Double3(boxlimitminouter)
      ,Double3(boxlimitmaxouter),Double3(boxlimitminmid),Double3(boxlimitmaxmid),inner,posxy,posz,code,listp,mat,tracking,nzone);
    cudaMemcpy(&count,listp+n,sizeof(unsigned),cudaMemcpyDeviceToHost);
   // KerTestCTE <<<sgrid,SPHBSIZE,smem>>>();
    //-Reorders list when stable has been activated.
    //-Reordena lista cuando stable esta activado.
//    if(stable && count){ //-Does not affect results.
//      thrust::device_ptr<unsigned> dev_list(listp);
//      thrust::sort(dev_list,dev_list+count);
//    }
  //  cudaDeviceSynchronize();
  }
  return(count);
}




__global__ void KerBufferComputeStep(unsigned n,int *inoutpart,const double2 *posxy,const double *posz
  ,typecode *code,const double3 boxlimitmininner,const double3 boxlimitmaxinner,const double3 boxlimitminouter
  ,const double3 boxlimitmaxouter,const bool inner,tmatrix4f mat,bool tracking,unsigned nzone)
{
	const unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
	if(cp<n){
		typecode cod=0;
		const int p=inoutpart[cp];
		const typecode rcode=code[p];
		const byte izone0=byte(CODE_GetIzoneFluidBuffer(rcode));
		const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK); //-Substract 16 to obtain the actual zone (0-15).
		double2 rxy=posxy[p];
		double rz=(posz[p]);
    if(tracking){
          double2 rxy_t=make_double2(0,0);
          double rz_t=0.0;
          MovePoint(rxy,rz,rxy_t,rz_t,mat);
          rxy=rxy_t; rz=rz_t;
        }
		if(izone==byte(nzone)){
			if(izone0>=16){//-Normal fluid particle in zone buffer
				cod= rcode^0x10 ; //-Converts to inout particle or not.
				code[p]=cod;
			}
			else{//-Previous buffr fluid particle.
				if((!KerBufferInZone(rxy,rz,boxlimitminouter,boxlimitmaxouter)&& inner) ||(KerBufferInZone(rxy,rz,boxlimitmininner,boxlimitmaxinner)&&!inner)){
					cod=CODE_SetOutIgnore(rcode); //-Particle is moved out domain.
					code[p]=cod;
				}
				if((!KerBufferInZone(rxy,rz,boxlimitminouter,boxlimitmaxouter)&& !inner) ||(KerBufferInZone(rxy,rz,boxlimitmininner,boxlimitmaxinner)&&inner)){
					cod=CODE_TYPE_FLUID; //-Particle become normal;
					code[p]=cod;
				}
			}
		}
	}
}

//==============================================================================
/// Checks particle position.
/// If particle is moved to fluid zone then it changes to fluid particle and
/// it creates a new in/out particle.
/// If particle is moved out the domain then it changes to ignore particle.
//==============================================================================
void BufferComputeStep(unsigned n,int *inoutpart,const double2 *posxy,const double *posz
  ,typecode *code,const tdouble3 boxlimitmininner,const tdouble3 boxlimitmaxinner,const tdouble3 boxlimitminouter,const tdouble3 boxlimitmaxouter,const bool inner,tmatrix4f mat,bool tracking,unsigned nzone)
{
  if(n){
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    KerBufferComputeStep <<<sgrid,SPHBSIZE>>> (n,inoutpart,posxy,posz,code,Double3(boxlimitmininner),Double3(boxlimitmaxinner),Double3(boxlimitminouter)
    	      ,Double3(boxlimitmaxouter),inner,mat,tracking,nzone);
  }
}


template<const int n, const int n1>
__device__ void LUdecomp_Single(float *a, int *p, float *b, float *sol, float &treshold){
    // // Norm matrix a
    float maxs = 0;
    #pragma unroll
    for (int i = 0; i < n; i++) {
        float sum = 0;
        #pragma unroll
        for (int j = 0; j < n; j++) 
            sum += fabs(a[T(i * n + j)]);
        maxs = max(maxs, sum);
    }

    #pragma unroll
    for (int i = 0; i < n; i++) {
        p[i] = i;
    }

    float maxv = 0;
    #pragma unroll
    for (int i = 0; i < n; i++) {
        // Pivoting
        int imax = i;
        #pragma unroll
        for (int k = i; k < n; k++) {
            if (fabs(a[T(i + k * n)]) > maxv) {
                maxv = fabs(a[T(i + k * n)]);
                imax = k;
            }
        }

        if (imax != i) {
            int tmp = p[i];
            p[i] = p[imax];
            p[imax] = tmp;

            #pragma unroll
            for (int j = 0; j < n; j++) {
                double temp = a[T(i * n + j)];
                a[T(i * n + j)] = a[T(imax * n + j)];
                a[T(imax * n + j)] = temp;
            }
        }

        // LU decomposition
        #pragma unroll
        for (int j = i + 1; j < n; j++) {
            a[T(j * n + i)] /= a[T(i * n + i)];
            #pragma unroll
            for (int k = i + 1; k < n; k++)
                a[T(j * n + k)] -= a[T(j * n + i)] * a[T(i * n + k)];
        }
    }

    float ia[n * n] = {0};

    // Matrix inversion
    #pragma unroll
    for (int j = 0; j < n; j++) {
        #pragma unroll
        for (int i = 0; i < n; i++) {
            ia[i * n + j] = (p[i] == j) ? 1.0 : 0.0;
            #pragma unroll
            for (int k = 0; k < i; k++)
                ia[i * n + j] -= a[T(i * n + k)] * ia[k * n + j];
        }
        #pragma unroll
        for (int i = n - 1; i >= 0; i--) {
            #pragma unroll
            for (int k = i + 1; k < n; k++)
                ia[i * n + j] -= a[T(i * n + k)] * ia[k * n + j];
            ia[i * n + j] /= a[T(i * n + i)];
        }
    }

    // Norm of inv matrix
    float maxs1 = 0;
    #pragma unroll
    for (int i = 0; i < n; i++) {
        double sum = 0;
        #pragma unroll
        for (int j = 0; j < n; j++) 
            sum += fabs(ia[i * n + j]);
        maxs1 = max(maxs1, sum);
    }

    treshold = 1.0f / (maxs * maxs1);

    // Solution
    #pragma unroll
    for (int k = 0; k < n1; k++)
        #pragma unroll
        for (int i = 0; i < n; i++) {
            sol[k] += ia[i] * b[i + k * n];
        }
}
//------------------------------------------------------------------------------
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//------------------------------------------------------------------------------
template<bool sim2d,TpKernel tker,unsigned order> __global__ void KerInteractionBufferExtrap_Single
  (unsigned bufferpartcount,const int *bufferpart,int scelldiv,int4 nc,int3 cellzero,const int2 *beginendcellfluid,const float4* poscell, double3 mapposmin
  ,const double2 *posxy,const double *posz,const typecode *code,const unsigned *idp
  ,const float4 *velrhop,const double2 *posxyb,const double *poszb,float4 *velrhopg,double *rcond,typecode *code1,float mrthreshold)
{
  // extern __shared__ float A[];
  float scaleh=CTE.kernelh*CTE.kernelh*CTE.kernelh*CTE.kernelh;
  // if(!threadIdx.x)slist[0]=0;
  const unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(cp<bufferpartcount){
    const unsigned p1=bufferpart[cp];

      //-Calculates ghost node position.
      double3 posp1=make_double3(posxyb[p1].x,posxyb[p1].y,poszb[p1]);
      const float4 gpscellp1=cusph::KerComputePosCell(posp1,mapposmin,CTE.poscellsize);

   //   if(CODE_IsPeriodic(code[p1]))pos_p1=KerInteraction_PosNoPeriodic(pos_p1);

      //-Initializes variables for calculation.
      constexpr unsigned nmatrix = (order == 0) ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
      constexpr unsigned nrhs = sim2d ? 3 : 4;

      float C[nmatrix]{0};
      float C1[nmatrix]{0};
      float D[nrhs]{0};


      float A[nmatrix*nmatrix]{0};
      float B[nmatrix*nrhs]{0};
      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::InitCte(posp1.x,posp1.y,posp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      //-Interaction with fluids.
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
        if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
          const float4 pscellp2=poscell[p2];
          float drx=gpscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(gpscellp1.w)-PSCEL_GetfX(pscellp2.w));
          float dry=gpscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(gpscellp1.w)-PSCEL_GetfY(pscellp2.w));
          float drz=gpscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(gpscellp1.w)-PSCEL_GetfZ(pscellp2.w));
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) &&!CODE_IsFluidFixed(code[p2])){//-Only with fluid particles but not inout particles.
            //-Only Wendland or Cubic Spline kernel.
            //-Computes kernel.
        	  float fac;
        	  float facc;
        	  const float wab=cufsph::GetKernel_WabFacFacc<tker>(rr2,fac,facc);
        	  const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
        	  float frxx=facc*(drx*drx)/sqrt(rr2)+fac*(dry*dry+drz*drz)/rr2;
        	  float frzz=facc*(drz*drz)/sqrt(rr2)+fac*(dry*dry+drx*drx)/rr2;
        	  float fryy=facc*(dry*dry)/sqrt(rr2)+fac*(drz*drz+drx*drx)/rr2;
        	  float frxz=facc*(drx*drz)/sqrt(rr2)-fac*(drx*drz)/rr2;
        	  float frxy=facc*(drx*dry)/sqrt(rr2)-fac*(drx*dry)/rr2;
        	  float fryz=facc*(dry*drz)/sqrt(rr2)-fac*(dry*drz)/rr2;


//        	double A[3]{wab,frx,frz};
//        	double A2[3]{velrhop[p2].w,velrhop[p2].x,velrhop[p2].z};

            // float C[6]{1.0f,drx,drz,drx*drx*0.5,drx*drz,drz*drz*0.5};
            // float C1[6]{volp2*wab,volp2*frx,volp2*frz,volp2*frxx,volp2*frxz,volp2*frzz};

//        	  double	 A3D[4]{wab,frx,fry,frz};
//        	double   A23D[4]{velrhop[p2].w,velrhop[p2].x,velrhop[p2].y,velrhop[p2].z};

             float4 velrhopp2=velrhop[p2];
            //===== Get mass and volume of particle p2 =====
            float massp2=CTE.massf;
            float volp2=massp2/velrhopp2.w;
            if constexpr (order == 0) {
        if constexpr (sim2d) {
            float tempC[] = {1.0f, drx, drz, drx * drx * 0.5f, drx * drz, drz * drz * 0.5f};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz, volp2 * frxx, volp2 * frxz, volp2 * frzz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        } else {
            float tempC[] = {1.0f, drx, dry, drz, drx * drx * 0.5f, drx * dry, dry * dry * 0.5f, drx * drz, drz * drz * 0.5f, dry * drz};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz, volp2 * frxx, volp2 * frxy, volp2 * fryy, volp2 * frxz, volp2 * frzz, volp2 * fryz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        }
    } else {
        if constexpr (sim2d) {
            float tempC[] = {1.0f, drx, drz};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        } else {
            float tempC[] = {1.0f, drx, dry, drz};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        }
    }

    if constexpr (sim2d) {
        float tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].z};
        #pragma unroll
        for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
    } else {
        float tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].y, velrhop[p2].z};
        #pragma unroll
        for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
    }

    // Your computation logic
    #pragma unroll
    for (int i = 0; i < nmatrix; i++) {
        #pragma unroll
        for (int j = 0; j < nmatrix; j++) {
            A[T(i * nmatrix + j)] += C1[i] * C[j];
        }
    }
    // #pragma unroll
    // for (int i = 0; i < nmatrix; i++) {
    //     #pragma unroll
    //     for (int j = 0; j < nmatrix; j++) {
    //         A1[threadIdx.x +32*(i * nmatrix + j)] += C[i] * C1[j];
    //     }
    // }
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


      float shep = A[0];
      float sol[nrhs]{0};
      float treshold=0;
      float cond=0;
      
      if (shep > 0.05){
        if (order==0){

          int P[nmatrix]{0};
          
          LUdecomp_Single<nmatrix, nrhs>(A, P, B, sol, treshold);
    
    
          cond=(1.0f/treshold)*scaleh;

        } else if (order==1){

          int P[nmatrix]{0};

        LUdecomp_Single<nmatrix, nrhs>(A, P, B, sol, treshold);
        cond=(1.0f/treshold)*CTE.kernelh*CTE.kernelh;


      }
      if (cond>mrthreshold || order==2){
         for (unsigned i = 0; i < nrhs; i++)
           sol[i] = B[i * nmatrix] / shep;
        //  rcond[p1] = 1.0;
         }
 

      if (sim2d)
      {
         
        velrhopg[p1].w = sol[0];
        velrhopg[p1].x = sol[1];
        velrhopg[p1].z = sol[2];
      } else {
        velrhopg[p1].w = sol[0];
        velrhopg[p1].x = sol[1];
        velrhopg[p1].y = sol[2];
        velrhopg[p1].z = sol[3];
      }
         
      } else{
        code1[p1] = CODE_SetOutIgnore(code1[p1]);

      }

  }
}


//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker> void Interaction_BufferExtrapT(unsigned bufferpartcount,const int *bufferpart,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float4 *velrhop,double *rcond,typecode *code1,bool fastsingle, const unsigned order,float mrthreshold)
{
const StDivDataGpu &dvd=t.divdatag;
  const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(bufferpartcount){
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(bufferpartcount,bsize);
    if(t.simulate2d){ const bool sim2d=true;

      if(fastsingle){         
        if (order==0) KerInteractionBufferExtrap_Single    <sim2d,tker,0> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,rcond,code1,mrthreshold);
        if (order==1) KerInteractionBufferExtrap_Single    <sim2d,tker,1> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,rcond,code1,mrthreshold);
        if (order==2) KerInteractionBufferExtrap_Single    <sim2d,tker,2> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,rcond,code1,mrthreshold);
}
   

    }    else{           const bool sim2d=false;
      if(fastsingle){   
        if (order==0) KerInteractionBufferExtrap_Single    <sim2d,tker,0> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,rcond,code1,mrthreshold);
        if (order==1) KerInteractionBufferExtrap_Single    <sim2d,tker,1> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,rcond,code1,mrthreshold);
        if (order==2) KerInteractionBufferExtrap_Single    <sim2d,tker,2> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,rcond,code1,mrthreshold);
         }
  }
}
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void Interaction_BufferExtrap(unsigned bufferpartcount,const int *bufferpart,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float4 *velrhop,double *rcond,typecode *code1,bool fastsingle, const unsigned order,float mrthreshold)
{
  switch(t.tkernel){
    case KERNEL_Wendland:
      Interaction_BufferExtrapT<KERNEL_Wendland>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,rcond,code1,fastsingle,order,mrthreshold);
    break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:
    	Interaction_BufferExtrapT<KERNEL_Wendland>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,rcond,code1,fastsingle,order,mrthreshold);
    break;
#endif
    default: throw "Kernel unknown at Interaction_InOutExtrap().";
  }
}

//------------------------------------------------------------------------------
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//------------------------------------------------------------------------------
template<bool sim2d,TpKernel tker,unsigned order> __global__ void KerInteractionBufferExtrap_SingleFlux
  (unsigned bufferpartcount,unsigned pini,int scelldiv,int4 nc,int3 cellzero,const int2 *beginendcellfluid,const float4* poscell, double3 mapposmin
  ,const double2 *posxy,const double *posz,const typecode *code,const unsigned *idp
  ,const float4 *velrhop,const double2 *posxyb,const double *poszb,float3 *normals,float *fluxes,
  double dp,double dt,float3* velflux,float mrthreshold,const unsigned cellfluid)
{
    float scaleh=CTE.kernelh*CTE.kernelh*CTE.kernelh*CTE.kernelh;

  const unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(cp<bufferpartcount){
    const unsigned p1=cp+pini;

      //-Calculates ghost node position.
      double3 posp1=make_double3(posxyb[p1].x,posxyb[p1].y,poszb[p1]);
            const float4 gpscellp1=cusph::KerComputePosCell(posp1,mapposmin,CTE.poscellsize);

  //-Initializes variables for calculation.
      constexpr unsigned nmatrix = (order == 0) ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
      constexpr unsigned nrhs = sim2d ? 3 : 4;

      float C[nmatrix];
      float C1[nmatrix];
      float D[nrhs];



      float A[nmatrix*nmatrix]{0};
      float B[nmatrix*nrhs]{0};
         float ShiftTFS=0;
         float mindist=1000000.0;
         float mindp=min(dp,CTE.dp);
      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::InitCte(posp1.x,posp1.y,posp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      // ini3-=cellfluid; fin3-=cellfluid;
      //-Interaction with fluids.
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
        if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
          const float4 pscellp2=poscell[p2];
          float drx=gpscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(gpscellp1.w)-PSCEL_GetfX(pscellp2.w));
          float dry=gpscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(gpscellp1.w)-PSCEL_GetfY(pscellp2.w));
          float drz=gpscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(gpscellp1.w)-PSCEL_GetfZ(pscellp2.w));
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){//-Only with fluid particles but not inout particles.
            //-Only Wendland or Cubic Spline kernel.
            //-Computes kernel.
        	  float fac;
        	  float facc;
        	  const float wab=cufsph::GetKernel_WabFacFacc<tker>(rr2,fac,facc);
        	  const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
             float4 velrhopp2=velrhop[p2];
            //===== Get mass and volume of particle p2 =====
            float massp2=CTE.massf;
            float volp2=massp2/velrhopp2.w;
            ShiftTFS-=volp2*(drx*frx+dry*fry+drz*frz);
          }
        }
      }
      cunsearch::InitCte(posp1.x,posp1.y,posp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);
      ini3+=cellfluid; fin3+=cellfluid;

      //-Interaction with fluids.
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
        if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
          const float4 pscellp2=poscell[p2];
          float drx=gpscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(gpscellp1.w)-PSCEL_GetfX(pscellp2.w));
          float dry=gpscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(gpscellp1.w)-PSCEL_GetfY(pscellp2.w));
          float drz=gpscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(gpscellp1.w)-PSCEL_GetfZ(pscellp2.w));
          const float rr2=drx*drx+dry*dry+drz*drz;
          if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) && !CODE_IsFluidFixed(code[p2])){//-Only with fluid particles but not inout particles.
            //-Only Wendland or Cubic Spline kernel.
            //-Computes kernel.
        	  float fac;
        	  float facc;
        	  const float wab=cufsph::GetKernel_WabFacFacc<tker>(rr2,fac,facc);
        	  const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
        	  float frxx=facc*(drx*drx)/sqrt(rr2)+fac*(dry*dry+drz*drz)/rr2;
        	  float frzz=facc*(drz*drz)/sqrt(rr2)+fac*(dry*dry+drx*drx)/rr2;
        	  float fryy=facc*(dry*dry)/sqrt(rr2)+fac*(drz*drz+drx*drx)/rr2;
        	  float frxz=facc*(drx*drz)/sqrt(rr2)-fac*(drx*drz)/rr2;
        	  float frxy=facc*(drx*dry)/sqrt(rr2)-fac*(drx*dry)/rr2;
        	  float fryz=facc*(dry*drz)/sqrt(rr2)-fac*(dry*drz)/rr2;






        	

             float4 velrhopp2=velrhop[p2];
            //===== Get mass and volume of particle p2 =====
            float massp2=CTE.massf;
            float volp2=massp2/velrhopp2.w;
            ShiftTFS-=volp2*(drx*frx+dry*fry+drz*frz);
            mindist=min(mindist,rr2);


            if constexpr (order == 0) {
        if constexpr (sim2d) {
            float tempC[] = {1.0f, drx, drz, drx * drx * 0.5f, drx * drz, drz * drz * 0.5f};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz, volp2 * frxx, volp2 * frxz, volp2 * frzz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        } else {
            float tempC[] = {1.0f, drx, dry, drz, drx * drx * 0.5f, drx * dry, dry * dry * 0.5f, drx * drz, drz * drz * 0.5f, dry * drz};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz, volp2 * frxx, volp2 * frxy, volp2 * fryy, volp2 * frxz, volp2 * frzz, volp2 * fryz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        }
    } else {
        if constexpr (sim2d) {
            float tempC[] = {1.0f, drx, drz};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * frz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        } else {
            float tempC[] = {1.0f, drx, dry, drz};
            float tempC1[] = {volp2 * wab, volp2 * frx, volp2 * fry, volp2 * frz};

            #pragma unroll
            for (int i = 0; i < nmatrix; ++i) {
                C[i] = tempC[i];
                C1[i] = tempC1[i];
            }
        }
    }

    if constexpr (sim2d) {
        float tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].z};
        #pragma unroll
        for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
    } else {
        float tempD[] = {velrhop[p2].w, velrhop[p2].x, velrhop[p2].y, velrhop[p2].z};
        #pragma unroll
        for (int i = 0; i < nrhs; ++i) D[i] = tempD[i];
    }

    #pragma unroll
    for (int i = 0; i < nmatrix; i++) {
        #pragma unroll
        for (int j = 0; j < nmatrix; j++) {
            A[T(i * nmatrix + j)] += C1[i] * C[j];
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

    float shep = A[0];
      float sol[nrhs]{0};
      float treshold=0;
      float cond=0;
      if (shep > 0.05){
        if (order==0){

          int P[nmatrix]{0};
          
          LUdecomp_Single<nmatrix, nrhs>(A, P, B, sol, treshold);
    
    
          cond=(1.0f/treshold)*scaleh;

        } else if (order==1){

        //   float A[9]{0};
          int P[nmatrix]{0};


        LUdecomp_Single<nmatrix, nrhs>(A, P, B, sol, treshold);
        cond=(1.0f/treshold)*CTE.kernelh*CTE.kernelh;


      }
       if (cond>mrthreshold|| order==2){
         for (unsigned i = 0; i < nrhs; i++)
           sol[i] = B[i * nmatrix] / shep;
        //  rcond[p1] = 1.0;
         }
      if (sim2d)
      {
        if ((ShiftTFS > 1.5 || sqrt(mindist) < mindp) && fluxes[p1] < 0.0 )
        fluxes[p1] += max(0.0, -sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt);
        else if ((ShiftTFS > 1.5 || sqrt(mindist) < mindp) )
        fluxes[p1] += -sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt;
      } else {
        if((ShiftTFS>2.75 || sqrt(mindist)<mindp)  && fluxes[p1]<0.0 && posp1.z>0.1*CTE.dp ) fluxes[p1]+=max(0.0,-sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt);
          else if((ShiftTFS>2.75 || sqrt(mindist)<mindp)  && posp1.z>0.1*CTE.dp)               fluxes[p1]+= -sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt;
          
      }
         
      } 

  }
}


//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker> void Interaction_BufferExtrapFluxT(unsigned bufferpartcount,unsigned pini,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float3 *normals,float *fluxes,double dp,double dt,float3* velflux,bool fastsingle,const unsigned order,float mrthreshold)
{
const StDivDataGpu &dvd=t.divdatag;
  const int2* beginendcellfluid=dvd.beginendcell/* +dvd.cellfluid */;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(bufferpartcount){
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(bufferpartcount,bsize);
    if(t.simulate2d){ const bool sim2d=true;
      if (fastsingle) {
        if (order==0) KerInteractionBufferExtrap_SingleFlux   <sim2d,tker,0> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);
        if (order==1) KerInteractionBufferExtrap_SingleFlux   <sim2d,tker,1> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);
        if (order==2) KerInteractionBufferExtrap_SingleFlux   <sim2d,tker,2> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);
       }
    }
    else{           const bool sim2d=false;
      if (fastsingle) {
        if (order==0) KerInteractionBufferExtrap_SingleFlux   <sim2d,tker,0> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);
        if (order==1) KerInteractionBufferExtrap_SingleFlux   <sim2d,tker,1> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);
        if (order==2) KerInteractionBufferExtrap_SingleFlux   <sim2d,tker,2> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);
     }
    }
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void Interaction_BufferExtrapFlux(unsigned bufferpartcount,unsigned pini,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float3 *normals,float *fluxes,double dp,double dt,float3* velflux,bool fastsingle,const unsigned order,float mrthreshold)
{
  switch(t.tkernel){
    case KERNEL_Wendland:
      Interaction_BufferExtrapFluxT<KERNEL_Wendland>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,order,mrthreshold);
    break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:
    	Interaction_BufferExtrapFluxT<KERNEL_Wendland>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,order,mrthreshold);
    break;
#endif
    default: throw "Kernel unknown at Interaction_InOutExtrap().";
  }
}

//------------------------------------------------------------------------------
/// Create list for new inlet particles to create.
/// Crea lista de nuevas particulas inlet a crear.
//------------------------------------------------------------------------------
__global__ void KerBufferListCreate(unsigned n,unsigned pini,unsigned nmax,float *fluxes,int *bufferpart,double massf)
{
  extern __shared__ unsigned slist[];
  if(!threadIdx.x)slist[0]=0;
  __syncthreads();
  const unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(cp<n && fluxes[cp+pini]>massf){
	  fluxes[cp+pini]-=massf;
    slist[atomicAdd(slist,1)+1]=cp;
  }
  __syncthreads();
  const unsigned ns=slist[0];
  __syncthreads();
  if(!threadIdx.x && ns)slist[0]= atomicAdd((bufferpart+nmax),ns);
  __syncthreads();
  if(threadIdx.x<ns){
    const unsigned cp2=slist[0]+threadIdx.x;
    if(cp2<nmax)bufferpart[cp2]=slist[threadIdx.x+1];
  }
}

//==============================================================================
/// Create list for new inlet particles to create at end of inoutpart[].
/// Returns number of new particles to create.
///
/// Crea lista de nuevas particulas inlet a crear al final de inoutpart[].
/// Devuelve el numero de las nuevas particulas para crear.
//==============================================================================
unsigned BufferListCreate(bool stable,unsigned n,unsigned pini,unsigned nmax, float *fluxes,int *bufferpart,double massf)
{
  unsigned count=0;
  if(n){
    //-inoutpart size list initialized to zero.
    //-Inicializa tamanho de lista inoutpart a cero.
    cudaMemset(bufferpart+nmax,0,sizeof(unsigned));
    dim3 sgrid=GetSimpleGridSize(n,SPHBSIZE);
    const unsigned smem=(SPHBSIZE+1)*sizeof(unsigned); //-All fluid particles can be in in/out area and one position for counter.
    KerBufferListCreate <<<sgrid,SPHBSIZE,smem>>> (n,pini,nmax,fluxes,bufferpart,massf);
    cudaMemcpy(&count,bufferpart+nmax,sizeof(unsigned),cudaMemcpyDeviceToHost);
    //-Reorders list if it is valid and stable has been activated.
    //-Reordena lista si es valida y stable esta activado.
//    if(stable && count && count<=nmax){
//      thrust::device_ptr<unsigned> dev_list((unsigned*)inoutpart);
//      thrust::sort(dev_list+n,dev_list+n+count);
    }

  return(count);
}


//------------------------------------------------------------------------------
/// Creates new inlet particles to replace the particles moved to fluid domain.
//------------------------------------------------------------------------------
template<bool periactive> __global__ void KerBufferCreateNewPart(unsigned newn,unsigned pini
  ,int *bufferpart
  ,unsigned np,unsigned idnext,const float3 *normals,const float dp
  ,double2 *posxy,double *posz,double2 *posxyb,double *poszb,unsigned *dcell,typecode *code,unsigned *idp,float4 *velrhop,unsigned nzone)
{
  const unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(cp<newn){
    const int p=bufferpart[cp];
    const double dis=0.5*dp;
    const float3 normal=normals[p+pini];
    double2 rposxy=posxyb[p+pini];
    double rposz=poszb[p+pini];
    rposxy.x-=dis*normal.x;
    rposxy.y-=dis*normal.y;
    rposz   -=dis*normal.z;
    const unsigned p2=np+cp;
    code[p2]=CODE_ToFluidBuffer(CODE_TYPE_FLUID,byte(nzone));
    cusph::KerUpdatePos<false>(rposxy,rposz,0,0,0,false,p2,posxy,posz,dcell,code);
    idp[p2]=idnext+cp;
    velrhop[p2]=make_float4(0,0,0,1000);
  }
}

//==============================================================================
/// Creates new inlet particles to replace the particles moved to fluid domain.
//==============================================================================
void BufferCreateNewPart(byte periactive,unsigned newn,unsigned pini
		  ,int *bufferpart
		  ,unsigned np,unsigned idnext,const float3 *normals,const float dp
		  ,double2 *posxy,double *posz,double2 *posxyb,double *poszb,unsigned *dcell,typecode *code,unsigned *idp,float4 *velrhop,unsigned nzone)
{
  if(newn){
    dim3 sgrid=GetSimpleGridSize(newn,SPHBSIZE);
    if(periactive) KerBufferCreateNewPart<true> <<<sgrid,SPHBSIZE>>> (newn,pini,bufferpart,np,idnext,normals,dp,posxy,posz,posxyb,poszb,dcell,code,idp,velrhop,nzone);
    else KerBufferCreateNewPart<false> <<<sgrid,SPHBSIZE>>> (newn,pini,bufferpart,np,idnext,normals,dp,posxy,posz,posxyb,poszb,dcell,code,idp,velrhop,nzone);

  }
}

__global__ void KerMoveBufferZone(unsigned pini,unsigned ntot, double2 *posxy,double *posz,float3* normals,float3* velflux,double dt,tmatrix4d mat)
{
  const unsigned cp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
  if(cp<ntot){
    const unsigned p1=cp+pini;
    const double2 rxy=posxy[p1];
    const double rz=posz[p1];
    double2 rxy_t=make_double2(0,0);
    double rz_t=0.0;
    rxy_t.x=rxy.x*mat.a11+rxy.y*mat.a12+rz*mat.a13+mat.a14;
    rxy_t.y=rxy.x*mat.a21+rxy.y*mat.a22+rz*mat.a23+mat.a24;
    rz_t   =rxy.x*mat.a31+rxy.y*mat.a32+rz*mat.a33+mat.a34;
    if(dt!=0)velflux[p1]=make_float3((rxy_t.x-rxy.x)/dt,(rxy_t.y-rxy.y)/dt,(rz_t-rz)/dt);
    posxy[p1]=rxy_t;
    posz[p1]=rz_t;
    float3 normalp1=normals[p1];
    normals[p1].x=normalp1.x*mat.a11+normalp1.y*mat.a12+normalp1.z*mat.a13;
    normals[p1].y=normalp1.x*mat.a21+normalp1.y*mat.a22+normalp1.z*mat.a23;
    normals[p1].z=normalp1.x*mat.a31+normalp1.y*mat.a32+normalp1.z*mat.a33;
  }
}

void MoveBufferZone(unsigned pini,unsigned ntot, double2 *posxy,double *posz,float3* normals,float3* velflux,double dt,tmatrix4d mat,int zone)
{
	dim3 sgrid=GetSimpleGridSize(ntot,SPHBSIZE);
	KerMoveBufferZone<<<sgrid,SPHBSIZE>>> (pini,ntot,posxy,posz,normals,velflux,dt,mat);
}


//------------------------------------------------------------------------------
__global__ void KerBufferShiftingGpu(unsigned n,unsigned pini,const double2 *posxy,const double *posz
		  ,float4 *shiftpos,typecode *code,const double3* boxlimitmin,const double3* boxlimitmax
      ,const bool* tracking,const tmatrix4f* mat,const bool* inner)
{
  unsigned p=blockIdx.x*blockDim.x + threadIdx.x;
  if(p<n){
	  const unsigned p1=p+pini;
	  typecode rcode=code[p1];
	  if(CODE_IsFluidBuffer(rcode)){

		  const byte izone0=byte(CODE_GetIzoneFluidBuffer(rcode));
		  const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
		  double2 rxy=posxy[p1];
      double rz=posz[p1];
      double3 boxmin=boxlimitmin[izone];
      double3 boxmax=boxlimitmax[izone];
      double3 origin  =make_double3((boxmax.x+boxmin.x)*0.5f,(boxmax.y+boxmin.y)*0.5f,(boxmax.z+boxmin.z)*0.5f);
      double3 boxsize =make_double3(boxmax.x-boxmin.x,boxmax.y-boxmin.y,boxmax.z-boxmin.z);
      if(tracking[izone]){
          double2 rxy_t=make_double2(0,0);
          double rz_t=0.0;
          MovePoint(rxy,rz,rxy_t,rz_t,mat[izone]);
          rxy=rxy_t; rz=rz_t;
        }
		  double disx =rxy.x  -origin.x;
      double disy =rxy.y  -origin.y;
		  double disz =rz     -origin.z;
      float3 shiftp1=make_float3(shiftpos[p1].x,shiftpos[p1].y,shiftpos[p1].z);
      // if(fabs(disx)>boxsize.x/2.0 && fabs(disz)>boxsize.z/2.0 && fabs(disy)>boxsize.y/2.0){
      //   shiftpos[p1].x=0.0f;
      //   shiftpos[p1].y=0.0f;
      //   shiftpos[p1].z=0.0f;
      // }
      bool inn=inner[izone];
      if((fabs(disx)>boxsize.x/2.0 )){
        // float3 normal=make_float3(mat[izone].a11,mat[izone].a21,mat[izone].a31);
        // shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        // shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        // shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                shiftpos[p1].x=0.0f;

      } 
      if((fabs(disy)>boxsize.y/2.0)){
        // float3 normal=make_float3(mat[izone].a12,mat[izone].a22,mat[izone].a32);
        // shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        // shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        // shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                shiftpos[p1].y=0.0f;

      }
      if((fabs(disz)>boxsize.z/2.0)){
        // float3 normal=make_float3(mat[izone].a13,mat[izone].a23,mat[izone].a33);
        // shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        // shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        // shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                        shiftpos[p1].z=0.0f;

      }
      if((fabs(disx)>boxsize.x/2.0) && (fabs(disz)>boxsize.z/2.0)){
        shiftpos[p1].x=0.0f;
        shiftpos[p1].z=0.0f;

      }

		  

  }

}
}


//==============================================================================
/// Computes sub-particle stress tensor (Tau) for SPS turbulence model.
//==============================================================================
void BufferShiftingGpu(unsigned np,unsigned npb,const double2 *posxy,const double *posz
  ,float4 *shiftpos,typecode *code,StrGeomVresGpu* data,cudaStream_t stm)
{
  const unsigned npf=np-npb;
  if(npf){
    dim3 sgridf=GetSimpleGridSize(npf,SPHBSIZE);
    KerBufferShiftingGpu <<<sgridf,SPHBSIZE>>> (npf,npb,posxy,posz,shiftpos,code,data->boxlimitmin,data->boxlimitmax,data->tracking,data->matmov,data->inner);
  }
}

}
