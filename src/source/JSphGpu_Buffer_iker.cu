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
namespace cusphbuffer{

#include "FunctionsBasic_iker.h"
#include "FunctionsMath_iker.h"
#include "FunSphKernel_iker.h"
#include "FunctionsGeo3d_iker.h"
#include "FunSphEos_iker.h"


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
			if(izone0>=16){     //-Normal fluid particle in zone buffer
				cod= rcode^0x10 ; //-Converts to buffer particle or not.
				code[p]=cod;
			}
			else{//-Previous buffer fluid particle.
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


template<typename T = float,const int m_dim, const int m_dim2>
__device__ void LUdecomp_Single(T (&A)[m_dim*m_dim], int (&p)[m_dim], T (&b)[m_dim*m_dim2], T (&sol)[m_dim2], T &treshold){
  //-Compute the norm of matrix A
  T maxs = 0;
  for(int i=0;i<m_dim;i++){
    T sum = 0;
    for(int j=0;j<m_dim;j++) sum+=fabs(A[i*m_dim+j]);
    maxs =max(maxs, sum);
  }

  //- Initialize permutation array
  for(int i=0;i<m_dim;i++) p[i]=i;
  

  T maxv = 0;
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
          T temp          =A[i*m_dim+j];
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

  T ia[m_dim*m_dim] = {0};

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
  T maxs1 = 0;
  for (int i=0;i<m_dim;i++) {
    T sum=0;
    for (int j=0;j<m_dim;j++) sum += fabs(ia[i*m_dim+j]);
    maxs1 =max(maxs1,sum);
  }

  treshold =1.0f/(maxs*maxs1);

  //-Compute Solution array
  for (int k=0;k<m_dim2;k++)
    for (int i=0;i<m_dim;i++) sol[k]+=ia[i]*b[i+k*m_dim];
        
}
//------------------------------------------------------------------------------
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//------------------------------------------------------------------------------
template <TpKernel tker,bool sim2d,TpVresOrder vrorder, typename T = float>
__global__ void KerInteractionBufferExtrap_Single(
    unsigned bufferpartcount, const int *bufferpart, int scelldiv, int4 nc, int3 cellzero, const int2 *beginendcellfluid, const float4 *poscell, double3 mapposmin,
    const double2 *posxy, const double *posz, const typecode *code, const unsigned *idp,
    const float4 *velrhop, const double2 *posxyb, const double *poszb, float4 *velrhopg, typecode *code1, float mrthreshold)
{
  const unsigned cp=blockIdx.x*blockDim.x+threadIdx.x; // Number of particle.
  if(cp<bufferpartcount){
    const unsigned p1=bufferpart[cp];

    // Calculates ghost node position.
    double3 posp1 = make_double3(posxyb[p1].x, posxyb[p1].y, poszb[p1]);
    const float4 gpscellp1 = cusph::KerComputePosCell(posp1, mapposmin, CTE.poscellsize);

    // Compute size of the reconstruction matrix and right hand sime.
    constexpr unsigned m_dim  = vrorder==VrOrder_2nd ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
    constexpr unsigned m_dim2 = sim2d ? 3: 4;

    T C[m_dim]{0};
    T C1[m_dim]{0};
    T D[m_dim2]{0};

    T A[m_dim*m_dim]{0};
    T B[m_dim*m_dim2]{0};

    // Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    cunsearch::InitCte(posp1.x,posp1.y,posp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    //-Interaction with fluids.
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
      if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
        T drx,dry,drz;
        if (sizeof(T)==sizeof(float)){
          const float4 pscellp2 = poscell[p2];
          drx=gpscellp1.x-pscellp2.x+CTE.poscellsize*(PSCEL_GetfX(gpscellp1.w)-PSCEL_GetfX(pscellp2.w));
          dry=gpscellp1.y-pscellp2.y+CTE.poscellsize*(PSCEL_GetfY(gpscellp1.w)-PSCEL_GetfY(pscellp2.w));
          drz=gpscellp1.z-pscellp2.z+CTE.poscellsize*(PSCEL_GetfZ(gpscellp1.w)-PSCEL_GetfZ(pscellp2.w));
        }else{
          const double2 p2xy=posxy[p2];
          drx=T(posp1.x-p2xy.x);
          dry=T(posp1.y-p2xy.y);
          drz=T(posp1.z-posz[p2]);
        }
        const T rr2=drx*drx+dry*dry+drz*drz;
        if (rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) && !CODE_IsFluidFixed(code[p2])){
        // Computes kernel.
          float fac;
          float facc;
          const T wab=cufsph::GetKernel_WabFacFacc<tker>(rr2,fac,facc);
          const T frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.
          T frxx=facc*(drx*drx)/sqrt(rr2)+fac*(dry*dry+drz*drz)/rr2;
          T frzz=facc*(drz*drz)/sqrt(rr2)+fac*(dry*dry+drx*drx)/rr2;
          T fryy=facc*(dry*dry)/sqrt(rr2)+fac*(drz*drz+drx*drx)/rr2;
          T frxz=facc*(drx*drz)/sqrt(rr2)-fac*(drx*drz)/rr2;
          T frxy=facc*(drx*dry)/sqrt(rr2)-fac*(drx*dry)/rr2;
          T fryz=facc*(dry*drz)/sqrt(rr2)-fac*(dry*drz)/rr2;

          float4 velrhopp2 = velrhop[p2];
          // Get mass and volume of particle p2
          T massp2=T(CTE.massf);
          T volp2=massp2/T(velrhopp2.w);

          // Set up C and C1 based on order and sim2d
          if (vrorder==VrOrder_2nd){
            if (sim2d){
              T tempC[] ={T(1.0),drx,drz,drx*drx*T(0.5),drx*drz,drz*drz*T(0.5)};
              T tempC1[]={volp2*wab,volp2*frx,volp2*frz,volp2*frxx,volp2*frxz,volp2*frzz};

              for (int i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                C1[i] = tempC1[i];
              }
            }else{
              T tempC[]   ={T(1.0),drx,dry,drz,drx*drx*T(0.5),drx*dry,dry*dry*T(0.5),drx*drz,drz*drz*T(0.5),dry*drz};
              T tempC1[]  ={volp2*wab,volp2*frx,volp2*fry,volp2*frz,volp2*frxx,volp2*frxy,volp2*fryy,volp2*frxz,volp2*frzz,volp2*fryz};

              for (int i=0;i<m_dim;++i){
                C[i] = tempC[i];
                C1[i] = tempC1[i];
              }
            }
          }else{
            if (sim2d){
              T tempC[]   = {T(1.0),drx,drz};
              T tempC1[]  = {volp2*wab,volp2*frx,volp2*frz};

              for (int i=0;i<m_dim;++i){
                    C[i] = tempC[i];
                    C1[i] = tempC1[i];
              }
            }else{
              T tempC[] = {T(1.0),drx,dry,drz};
              T tempC1[] = {volp2*wab,volp2*frx,volp2*fry,volp2*frz};

              for (int i =0;i<m_dim;++i){
                C[i] = tempC[i];
                C1[i] = tempC1[i];
              }
            }
          }

              // Set up D
          if (sim2d){
            T tempD[] = {T(velrhop[p2].w),T(velrhop[p2].x),T(velrhop[p2].z)};
            for (int i =0;i<m_dim2;++i)D[i]=tempD[i];
          }else{
            T tempD[] = {T(velrhop[p2].w), T(velrhop[p2].x), T(velrhop[p2].y), T(velrhop[p2].z)};
            for (int i =0;i<m_dim2;++i)D[i]=tempD[i];
          }
          // Accumulate A and B matrices
          for (int i=0;i<m_dim;i++)
            for (int j=0;j<m_dim;j++){
              A[i*m_dim+j]+=C1[i]*C[j];
            }
              
          for (int i=0; i<m_dim2;i++)
            for (int j=0;j<m_dim;j++){
              B[i*m_dim+j]+=D[i]*C1[j];
          }
        }
      }        
    }

    T shep = A[0];                                  ///<Shepard summation.
    T sol[m_dim2]{0};                               ///<Solution array;
    T treshold = T(0);                              ///<condition number;
    T cond = T(0);                                  ///<Scaled reciprocal condition number;
    T kernelh2=CTE.kernelh*CTE.kernelh;             ///<Scaling factor;

    if(shep>T(0.05)){
      if (vrorder==VrOrder_2nd){
        T scaleh = kernelh2*kernelh2;

        int P[m_dim]{0};

        LUdecomp_Single<T,m_dim,m_dim2>(A,P,B,sol,treshold);

        cond=(T(1.0)/treshold)*scaleh;
      }else if (vrorder==VrOrder_1st){

        int P[m_dim]{0};

        LUdecomp_Single<T,m_dim,m_dim2>(A,P,B,sol,treshold);
        cond=(T(1.0)/treshold)*kernelh2;
      }
      if (cond>mrthreshold || vrorder==VrOrder_0th)
      {
        for (unsigned i=0;i<m_dim2;i++)
          sol[i]=B[i*m_dim]/shep;
      }

      if (sim2d){
        velrhopg[p1].w = float(sol[0]);
        velrhopg[p1].x = float(sol[1]);
        velrhopg[p1].z = float(sol[2]);
      }else{
        velrhopg[p1].w = float(sol[0]);
        velrhopg[p1].x = float(sol[1]);
        velrhopg[p1].y = float(sol[2]);
        velrhopg[p1].z = float(sol[3]);
      }
    }
    else{
      code1[p1] = CODE_SetOutIgnore(code1[p1]);
    }
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker,bool sim2d,TpVresOrder vrorder>
 void Interaction_BufferExtrapT(unsigned bufferpartcount,const int *bufferpart,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float4 *velrhop,typecode *code1,bool fastsingle,float mrthreshold)
{
  const StDivDataGpu &dvd=t.divdatag;
  const int2* beginendcellfluid=dvd.beginendcell+dvd.cellfluid;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(bufferpartcount){
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(bufferpartcount,bsize);
    if(fastsingle)  KerInteractionBufferExtrap_Single<tker,sim2d,vrorder,float> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero
      ,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,code1,mrthreshold);
    // else            KerInteractionBufferExtrap_Single<tker,sim2d,vrorder,double> <<<sgrid,bsize>>> (bufferpartcount,bufferpart,dvd.scelldiv,dvd.nc,dvd.cellzero
    //   ,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,velrhop,code1,mrthreshold);     
  }
}


//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker>
void Interaction_BufferExtrap_gt0(unsigned bufferpartcount,const int *bufferpart,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float4* velrhop,typecode *code1,bool fastsingle,const TpVresOrder vrorder,float mrthreshold)
{
  if(t.simulate2d){
    switch(vrorder){
      // case VrOrder_0th: Interaction_BufferExtrapT<tker,true,VrOrder_0th>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,mrthreshold); break;
      // case VrOrder_1st: Interaction_BufferExtrapT<tker,true,VrOrder_1st>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,mrthreshold); break;
      case VrOrder_2nd: Interaction_BufferExtrapT<tker,true,VrOrder_2nd>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,mrthreshold); break;
    }
  }else{
      switch(vrorder){
      // case VrOrder_0th: Interaction_BufferExtrapT<tker,false,VrOrder_0th>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,mrthreshold); break;
      // case VrOrder_1st: Interaction_BufferExtrapT<tker,false,VrOrder_1st>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,mrthreshold); break;
      // case VrOrder_2nd: Interaction_BufferExtrapT<tker,false,VrOrder_2nd>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,mrthreshold); break;
    }
  } 
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void Interaction_BufferExtrap(unsigned bufferpartcount,const int *bufferpart,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float4* velrhop,typecode *code1,bool fastsingle,const TpVresOrder order,float mrthreshold)
{
  switch(t.tkernel){
    case KERNEL_Wendland:
      Interaction_BufferExtrap_gt0<KERNEL_Wendland>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,order,mrthreshold);
    break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:
    	// Interaction_BufferExtrap_gt0<KERNEL_Wendland>(bufferpartcount,bufferpart,t,posxyb,poszb,velrhop,code1,fastsingle,order,mrthreshold);
    break;
#endif
    default: throw "Kernel unknown at Interaction_InOutExtrap().";
  }
}

//------------------------------------------------------------------------------
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//------------------------------------------------------------------------------
template <TpKernel tker,bool sim2d,TpVresOrder vrorder, typename T = float>
__global__ void KerInteractionBufferExtrap_SingleFlux
  (unsigned bufferpartcount,unsigned pini,int scelldiv,int4 nc,int3 cellzero,const int2 *beginendcellfluid,const float4* poscell, double3 mapposmin
  ,const double2 *posxy,const double *posz,const typecode *code,const unsigned *idp
  ,const float4 *velrhop,const double2 *posxyb,const double *poszb,float3 *normals,float *fluxes,
  double dp,double dt,float3* velflux,float mrthreshold,const unsigned cellfluid)
{
  const unsigned cp=blockIdx.x*blockDim.x+threadIdx.x; //-Number of particle.
  if(cp<bufferpartcount){
    const unsigned p1=cp+pini;

      //-Calculates ghost node position.
    double3 posp1 = make_double3(posxyb[p1].x, posxyb[p1].y, poszb[p1]);
    const float4 gpscellp1 = cusph::KerComputePosCell(posp1,mapposmin,CTE.poscellsize);

    // Compute size of the reconstruction matrix and right hand sime.
    constexpr unsigned m_dim  = vrorder==VrOrder_2nd ? (sim2d ? 6 : 10) : (sim2d ? 3 : 4);
    constexpr unsigned m_dim2 = sim2d ? 3: 4;

    T C[m_dim]{0};
    T C1[m_dim]{0};
    T D[m_dim2]{0};

    T A[m_dim*m_dim]{0};
    T B[m_dim*m_dim2]{0};



    
    T ShiftTFS=0;
    T mindist=1000000.0;
    T mindp=min(dp,CTE.dp);
    //-Obtains neighborhood search limits.
    int ini1,fin1,ini2,fin2,ini3,fin3;
    // cunsearch::InitCte(posp1.x,posp1.y,posp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

    // //-Interaction with fluids.
    // for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
    //   unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
    //   if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
    //     T drx,dry,drz;
    //     if (sizeof(T)==sizeof(float)){
    //       const float4 pscellp2 = poscell[p2];
    //       drx=gpscellp1.x-pscellp2.x+CTE.poscellsize*(PSCEL_GetfX(gpscellp1.w)-PSCEL_GetfX(pscellp2.w));
    //       dry=gpscellp1.y-pscellp2.y+CTE.poscellsize*(PSCEL_GetfY(gpscellp1.w)-PSCEL_GetfY(pscellp2.w));
    //       drz=gpscellp1.z-pscellp2.z+CTE.poscellsize*(PSCEL_GetfZ(gpscellp1.w)-PSCEL_GetfZ(pscellp2.w));
    //     }else{
    //       const double2 p2xy=posxy[p2];
    //       drx=T(posp1.x-p2xy.x);
    //       dry=T(posp1.y-p2xy.y);
    //       drz=T(posp1.z-posz[p2]);
    //     }
    //     const T rr2=drx*drx+dry*dry+drz*drz;
    //     if (rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
    //     // Computes kernel.
    //       float fac;
    //       float facc;
    //       const T wab=cufsph::GetKernel_WabFacFacc<tker>(rr2,fac,facc);
    //       const T frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.

    //       float4 velrhopp2 = velrhop[p2];
    //       // Get mass and volume of particle p2
    //       T massp2=T(CTE.massf);
    //       T volp2=massp2/T(velrhopp2.w);
    //       ShiftTFS-=volp2*(drx*frx+dry*fry+drz*frz);

    //     }
    //   }
    // }

    cunsearch::InitCte(posp1.x,posp1.y,posp1.z,scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);
    ini3+=cellfluid; fin3+=cellfluid;

    //-Interaction with fluids.
    for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
      unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,beginendcellfluid,pini,pfin);
      if(pfin)for(unsigned p2=pini;p2<pfin;p2++){
        T drx,dry,drz;
        if (sizeof(T)==sizeof(float)){
          const float4 pscellp2 = poscell[p2];
          drx=gpscellp1.x-pscellp2.x+CTE.poscellsize*(PSCEL_GetfX(gpscellp1.w)-PSCEL_GetfX(pscellp2.w));
          dry=gpscellp1.y-pscellp2.y+CTE.poscellsize*(PSCEL_GetfY(gpscellp1.w)-PSCEL_GetfY(pscellp2.w));
          drz=gpscellp1.z-pscellp2.z+CTE.poscellsize*(PSCEL_GetfZ(gpscellp1.w)-PSCEL_GetfZ(pscellp2.w));
        }else{
          const double2 p2xy=posxy[p2];
          drx=T(posp1.x-p2xy.x);
          dry=T(posp1.y-p2xy.y);
          drz=T(posp1.z-posz[p2]);
        }
        const T rr2=drx*drx+dry*dry+drz*drz;
        if (rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO && CODE_IsFluid(code[p2]) && !CODE_IsFluidBuffer(code[p2]) && !CODE_IsFluidFixed(code[p2])){
        // Computes kernel.
          float fac;
          float facc;
          const T wab=cufsph::GetKernel_WabFacFacc<tker>(rr2,fac,facc);
          const T frx=drx*fac,fry=dry*fac,frz=drz*fac; //-Gradients.
          T frxx=facc*(drx*drx)/sqrt(rr2)+fac*(dry*dry+drz*drz)/rr2;
          T frzz=facc*(drz*drz)/sqrt(rr2)+fac*(dry*dry+drx*drx)/rr2;
          T fryy=facc*(dry*dry)/sqrt(rr2)+fac*(drz*drz+drx*drx)/rr2;
          T frxz=facc*(drx*drz)/sqrt(rr2)-fac*(drx*drz)/rr2;
          T frxy=facc*(drx*dry)/sqrt(rr2)-fac*(drx*dry)/rr2;
          T fryz=facc*(dry*drz)/sqrt(rr2)-fac*(dry*drz)/rr2;

          float4 velrhopp2 = velrhop[p2];
          // Get mass and volume of particle p2
          T massp2=T(CTE.massf);
          T volp2=massp2/T(velrhopp2.w);
          ShiftTFS-=volp2*(drx*frx+dry*fry+drz*frz);
          mindist=min(mindist,rr2);


          // Set up C and C1 based on order and sim2d
          if (vrorder==VrOrder_2nd){
            if (sim2d){
              T tempC[] ={T(1.0),drx,drz,drx*drx*T(0.5),drx*drz,drz*drz*T(0.5)};
              T tempC1[]={volp2*wab,volp2*frx,volp2*frz,volp2*frxx,volp2*frxz,volp2*frzz};

              for (int i=0;i<m_dim;++i){
                C[i]  = tempC[i];
                C1[i] = tempC1[i];
              }
            }else{
              T tempC[]   ={T(1.0),drx,dry,drz,drx*drx*T(0.5),drx*dry,dry*dry*T(0.5),drx*drz,drz*drz*T(0.5),dry*drz};
              T tempC1[]  ={volp2*wab,volp2*frx,volp2*fry,volp2*frz,volp2*frxx,volp2*frxy,volp2*fryy,volp2*frxz,volp2*frzz,volp2*fryz};

              for (int i=0;i<m_dim;++i){
                C[i] = tempC[i];
                C1[i] = tempC1[i];
              }
            }
          }else{
            if (sim2d){
              T tempC[]   = {T(1.0),-drx,-drz};
              T tempC1[]  = {volp2*wab,volp2*frx,volp2*frz};

              for (int i=0;i<m_dim;++i){
                    C[i] = tempC[i];
                    C1[i] = tempC1[i];
              }
            }else{
              T tempC[] = {T(1.0),drx,dry,drz};
              T tempC1[] = {volp2*wab,volp2*frx,volp2*fry,volp2*frz};

              for (int i =0;i<m_dim;++i){
                C[i] = tempC[i];
                C1[i] = tempC1[i];
              }
            }
          }

              // Set up D
          if (sim2d){
            T tempD[] = {T(velrhop[p2].w),T(velrhop[p2].x),T(velrhop[p2].z)};
            for (int i =0;i<m_dim2;++i)D[i]=tempD[i];
          }else{
            T tempD[] = {T(velrhop[p2].w), T(velrhop[p2].x), T(velrhop[p2].y), T(velrhop[p2].z)};
            for (int i =0;i<m_dim2;++i)D[i]=tempD[i];
          }
          // Accumulate A and B matrices
          for (int i=0;i<m_dim;i++)
            for (int j=0;j<m_dim;j++){
              A[i*m_dim+j]+=C1[i]*C[j];
            }
              
          for (int i=0; i<m_dim2;i++)
            for (int j=0;j<m_dim;j++){
              B[i*m_dim+j]+=D[i]*C1[j];
          }
        }
      }        
    }

    T shep = A[0];                                  ///<Shepard summation.
    T sol[m_dim2]{0};                               ///<Solution array;
    T treshold = T(0);                              ///<condition number;
    T cond = T(0);                                  ///<Scaled reciprocal condition number;
    T kernelh2=CTE.kernelh*CTE.kernelh;             ///<Scaling factor;

    if(shep>T(0.05)){
      if (vrorder==VrOrder_2nd){
        T scaleh = kernelh2*kernelh2;

        int P[m_dim]{0};

        LUdecomp_Single<T,m_dim,m_dim2>(A,P,B,sol,treshold);

        cond=(T(1.0)/treshold)*scaleh;
      }else if (vrorder==VrOrder_1st){

        int P[m_dim]{0};

        LUdecomp_Single<T,m_dim,m_dim2>(A,P,B,sol,treshold);
        cond=(T(1.0)/treshold)*kernelh2;
      }
      if (cond>mrthreshold || vrorder==VrOrder_0th)
      {
        for (unsigned i=0;i<m_dim2;i++)
          sol[i]=B[i*m_dim]/shep;
      }
          // printf("%f %f %f %f\n",sol[0],sol[1],sol[2],dp);

      if (sim2d){
        if((ShiftTFS>1.5 || sqrt(mindist)<mindp) && fluxes[p1]<0.0)
        fluxes[p1]+=max(0.0,-sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt);
        else if((ShiftTFS>1.5 || sqrt(mindist)<mindp))
        fluxes[p1]+=-sol[0]*((-float(velflux[p1].x)+sol[1])*normals[p1].x+(-float(velflux[p1].z)+sol[2])*normals[p1].z)*dp*dt;
      }else{
        if      ((ShiftTFS>2.75 || sqrt(mindist)<mindp)  &&fluxes[p1]<0.0)  fluxes[p1]+=max(0.0,-sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt);
        else if ((ShiftTFS>2.75 || sqrt(mindist)<mindp))                    fluxes[p1]+= -sol[0]*((-velflux[p1].x+sol[1])*normals[p1].x+(-velflux[p1].y+sol[2])*normals[p1].y+(-velflux[p1].z+sol[3])*normals[p1].z)*dp*dp*dt;
          
      }
    }
  }
}


//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker,bool sim2d,TpVresOrder vrorder>
void Interaction_BufferExtrapFluxT(unsigned bufferpartcount,unsigned pini,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float3 *normals,float *fluxes,double dp,double dt,float3* velflux,bool fastsingle,float mrthreshold)
{
  const StDivDataGpu &dvd=t.divdatag;
  const int2* beginendcellfluid=dvd.beginendcell/* +dvd.cellfluid */;
  //-Interaction GhostBoundaryNodes-Fluid.
  if(bufferpartcount){
    const unsigned bsize=128;
    dim3 sgrid=GetSimpleGridSize(bufferpartcount,bsize);
    if(fastsingle)  KerInteractionBufferExtrap_SingleFlux<tker,sim2d,vrorder,float> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero
      ,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);
    // else KerInteractionBufferExtrap_SingleFlux<tker,sim2d,vrorder,double> <<<sgrid,bsize>>> (bufferpartcount,pini,dvd.scelldiv,dvd.nc,dvd.cellzero
    //   ,beginendcellfluid,t.poscell,Double3(t.mapposmin),t.posxy,t.posz,t.code,t.idp,t.velrho,posxyb,poszb,normals,fluxes,dp,dt,velflux,mrthreshold,dvd.cellfluid);    
  }
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
template<TpKernel tker>
void Interaction_BufferExtrapFlux_gt0(unsigned bufferpartcount,unsigned pini,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float3 *normals,float *fluxes,double dp,double dt,float3* velflux,bool fastsingle,const TpVresOrder vrorder,float mrthreshold)
{
  if(t.simulate2d){
    switch(vrorder){
      case VrOrder_0th: Interaction_BufferExtrapFluxT<tker,true,VrOrder_0th>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,mrthreshold); break;
      case VrOrder_1st: Interaction_BufferExtrapFluxT<tker,true,VrOrder_1st>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,mrthreshold); break;
      case VrOrder_2nd: Interaction_BufferExtrapFluxT<tker,true,VrOrder_2nd>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,mrthreshold); break;
    }
  }else{
      switch(vrorder){
      // case VrOrder_0th: Interaction_BufferExtrapFluxT<tker,false,VrOrder_0th>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,mrthreshold); break;
      // case VrOrder_1st: Interaction_BufferExtrapFluxT<tker,false,VrOrder_1st>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,mrthreshold); break;
      // case VrOrder_2nd: Interaction_BufferExtrapFluxT<tker,false,VrOrder_2nd>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,mrthreshold); break;
    }
  } 
}

//==============================================================================
/// Perform interaction between ghost inlet/outlet nodes and fluid particles. GhostNodes-Fluid
/// Realiza interaccion entre ghost inlet/outlet nodes y particulas de fluido. GhostNodes-Fluid
//==============================================================================
void Interaction_BufferExtrapFlux(unsigned bufferpartcount,unsigned pini,const StInterParmsbg &t,
		const double2 *posxyb,const double *poszb,float3 *normals,float *fluxes,double dp,double dt,float3* velflux,bool fastsingle,const TpVresOrder vrorder,float mrthreshold)
{
  switch(t.tkernel){
    case KERNEL_Wendland:
      Interaction_BufferExtrapFlux_gt0<KERNEL_Wendland>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,vrorder,mrthreshold);
    break;
#ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:
    	// Interaction_BufferExtrapFlux_gt0<KERNEL_Wendland>(bufferpartcount,pini,t,posxyb,poszb,normals,fluxes,dp,dt,velflux,fastsingle,vrorder,mrthreshold);
    break;
#endif
    default: throw "Kernel unknown at Interaction_InOutExtrap().";
  }
  cudaDeviceSynchronize();
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
        float3 normal=make_float3(mat[izone].a11,mat[izone].a21,mat[izone].a31);
        shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                // shiftpos[p1].x=0.0f;

      } 
      if((fabs(disy)>boxsize.y/2.0)){
        float3 normal=make_float3(mat[izone].a12,mat[izone].a22,mat[izone].a32);
        shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                // shiftpos[p1].y=0.0f;

      }
      if((fabs(disz)>boxsize.z/2.0)){
        float3 normal=make_float3(mat[izone].a13,mat[izone].a23,mat[izone].a33);
        shiftpos[p1].x=shiftp1.x-normal.x*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].y=shiftp1.y-normal.y*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
        shiftpos[p1].z=shiftp1.z-normal.z*(shiftp1.x*normal.x+shiftp1.y*normal.y+shiftp1.z*normal.z);
                        // shiftpos[p1].z=0.0f;

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

template<bool symm>  __device__ void KerComputeNormalsBufferBox(unsigned p1,const double3 posp1
  ,float massp2,float& fs_treshold,float3& gradc,tmatrix3f& lcorr,unsigned& neigh
  ,float& pou,const double3 boxlimitmin,const double3 boxlimitmax,const bool inner
  ,const tmatrix4f mat,const bool tracking,const bool sim2d)
{    
      
  float3 minpos=make_float3(0,0,0);
    minpos.x=posp1.x-CTE.kernelsize;
    minpos.y=(sim2d? 0.0f: posp1.y-CTE.kernelsize);
    minpos.z=posp1.z-CTE.kernelsize;
    float3 maxpos=make_float3(0,0,0);
    maxpos.x=posp1.x+CTE.kernelsize;
    maxpos.y=(sim2d? 0.0f: posp1.y+CTE.kernelsize);
    maxpos.z=posp1.z+CTE.kernelsize;

  float dp=CTE.dp;
  for (float rx=minpos.x; rx<=maxpos.x; rx+=dp) for (float ry=minpos.y; ry<=maxpos.y; ry+=dp)
    for (float rz=minpos.z; rz<=maxpos.z; rz+=dp){
      const float drx=float(posp1.x-rx);
      const float dry=float(posp1.y-ry);
      const float drz=float(posp1.z-rz);
      const float rr2=drx*drx+dry*dry+drz*drz;
      float rx1=rx; float ry1=ry; float rz1=rz;
      if(tracking){          
        rx1=rx*mat.a11+ry*mat.a21+rz*mat.a31+mat.a14;
        ry1=rx*mat.a12+ry*mat.a22+rz*mat.a32+mat.a24;
        rz1=rx*mat.a13+ry*mat.a23+rz*mat.a33+mat.a34;
      }
      bool outside=(inner ? !KerBufferInZone(make_double2(rx1,ry1),rz1,boxlimitmin,boxlimitmax) : KerBufferInZone(make_double2(rx1,ry1),rz1,boxlimitmin,boxlimitmax));
      if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO && outside){
        //-Computes kernel.
        const float fac=cufsph::GetKernel_Fac<KERNEL_Wendland>(rr2);
        const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
            
        const float vol2=massp2/CTE.rhopzero;
        neigh++;

        const float dot3=drx*frx+dry*fry+drz*frz;
        gradc.x+=vol2*frx;
        gradc.y+=vol2*fry;
        gradc.z+=vol2*frz;

        fs_treshold-=vol2*dot3;
        lcorr.a11+=-drx*frx*vol2; lcorr.a12+=-drx*fry*vol2; lcorr.a13+=-drx*frz*vol2;
        lcorr.a21+=-dry*frx*vol2; lcorr.a22+=-dry*fry*vol2; lcorr.a23+=-dry*frz*vol2;
        lcorr.a31+=-drz*frx*vol2; lcorr.a32+=-drz*fry*vol2; lcorr.a33+=-drz*frz*vol2;

        const float wab=cufsph::GetKernel_Wab<KERNEL_Wendland>(rr2);
        pou+=wab*vol2;       
      }
    }
  }



  template<bool symm>  __device__ void KerComputeNormalsBox(bool boundp2,unsigned p1
    ,const unsigned &pini,const unsigned &pfin,const float4 *poscell
    ,const float4* velrhop,const typecode* code,float massp2,const float4 &pscellp1
    ,float& fs_treshold,float3& gradc,tmatrix3f& lcorr,unsigned& neigh,float& pou,const float* ftomassp)
  {
    const float w0=cufsph::GetKernel_Wab<KERNEL_Wendland>(CTE.dp*CTE.dp);
    for(int p2=pini;p2<pfin;p2++){
    const float4 pscellp2=poscell[p2];
      float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
      float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
      float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
      if(symm)dry=pscellp1.y+pscellp2.y + CTE.poscellsize*PSCEL_GetfY(pscellp2.w); //<vs_syymmetry>
      const double rr2=drx*drx+dry*dry+drz*drz;
      if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){
        //-Computes kernel.
        const float fac=cufsph::GetKernel_Fac<KERNEL_Wendland>(rr2);
        const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
        // float4 velrhop2=velrhop[p2];
        // if(symm)velrhop2.y=-velrhop2.y; //<vs_syymmetry>
        const float rhopp2= float(velrhop[p2].w);
        //-Velocity derivative (Momentum equation).

        bool ftp2;
        float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
          const typecode cod=code[p2];
          ftp2=CODE_IsFloating(cod);
          ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);

        const float vol2=(ftp2 ? float(ftmassp2/rhopp2) : float(massp2/rhopp2));
        neigh++;

        const float dot3=drx*frx+dry*fry+drz*frz;
        gradc.x+=vol2*frx;
        gradc.y+=vol2*fry;
        gradc.z+=vol2*frz;

        fs_treshold-=vol2*dot3;
        lcorr.a11+=-drx*frx*vol2; lcorr.a12+=-drx*fry*vol2; lcorr.a13+=-drx*frz*vol2;
        lcorr.a21+=-dry*frx*vol2; lcorr.a22+=-dry*fry*vol2; lcorr.a23+=-dry*frz*vol2;
        lcorr.a31+=-drz*frx*vol2; lcorr.a32+=-drz*fry*vol2; lcorr.a33+=-drz*frz*vol2;

        const float wab=cufsph::GetKernel_Wab<KERNEL_Wendland>(rr2);
        pou+=wab*vol2;       

      }
    }
  }



  template<bool symm>
    __global__ void KerComputeNormals(unsigned n,unsigned pinit
    ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
    ,const float4 *poscell,const float4 *velrhop,const typecode *code
    ,unsigned* fstype,float3* fsnormal,bool simulate2d,float4* shiftposfs,const float* ftomassp,const unsigned* listp
    ,const double2* posxy,const double* posz,const double3* boxlimitmin,const double3* boxlimitmax,const bool* inner,const tmatrix4f* mat,const bool* tracking)
  {
    const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
    if(p<n){
      const unsigned p1=listp[p];   
      //-Obtains basic data of particle p1.
      const float4 pscellp1=poscell[p1];
      const float4 velrhop1=velrhop[p1];
      const bool rsymp1=(symm && PSCEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>
      
      float fs_treshold=0.f;                            //-Divergence of the position.
      float3 gradc=make_float3(0,0,0);                  //-Gradient of the concentration
      float pou=0;                                      //-Partition of unity.
      unsigned neigh=0;                                 //-Number of neighbours.
      
      tmatrix3f lcorr; cumath::Tmatrix3fReset(lcorr);             //-Correction matrix.
      tmatrix3f lcorr_inv; cumath::Tmatrix3fReset(lcorr_inv);     //-Inverse of the correction matrix.

      //-Calculate approx. number of neighbours when uniform distribution (in the future on the constant memory?)
      float Nzero=0.f;
      if(simulate2d){
        Nzero=(3.141592)*CTE.kernelsize2/(CTE.dp*CTE.dp);
      } else{
        Nzero=(4.f/3.f)*(3.141592)*CTE.kernelsize2*CTE.kernelsize/(CTE.dp*CTE.dp*CTE.dp);
      }

      //-Copy the value of shift to gradc. For single resolution is zero, but in Vres take in account virtual stencil.
      gradc=make_float3(shiftposfs[p1].x,shiftposfs[p1].y,shiftposfs[p1].z); 
      

    
      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      //-Interaction with fluids.
      ini3+=cellfluid; fin3+=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin){
          KerComputeNormalsBox<false> (false,p1,pini,pfin,poscell,velrhop,code,CTE.massf,pscellp1,fs_treshold,gradc,lcorr,neigh,pou,ftomassp);
          if(symm && rsymp1)KerComputeNormalsBox<true > (false,p1,pini,pfin,poscell,velrhop,code,CTE.massf,pscellp1,
                                                          fs_treshold,gradc,lcorr,neigh,pou,ftomassp); //<vs_syymmetry>
        }
      }

      // -Interaction with boundaries.
      ini3-=cellfluid; fin3-=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin){
          KerComputeNormalsBox<false> (false,p1,pini,pfin,poscell,velrhop,code,CTE.massf,pscellp1,fs_treshold,gradc,lcorr,neigh,pou,ftomassp);
          if(symm && rsymp1)KerComputeNormalsBox<true > (false,p1,pini,pfin,poscell,velrhop,code,CTE.massf,pscellp1,
                                                          fs_treshold,gradc,lcorr,neigh,pou,ftomassp); //<vs_syymmetry>
        }
      }

      // -Interaction with virtual stencil.
      if(CODE_IsFluidBuffer(code[p1])){
        const byte izone0=byte(CODE_GetIzoneFluidBuffer(code[p1]));
        const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
        const double3   boxlimmin =boxlimitmin[izone];
        const double3   boxlimmax =boxlimitmax[izone];
        const bool      inn       =inner      [izone];
        const tmatrix4f rmat      =mat        [izone];
        const bool      track     =tracking   [izone]; 
        const double3   posp1     =make_double3(posxy[p1].x,posxy[p1].y,posz[p1]); 
        KerComputeNormalsBufferBox<false> (p1,posp1,CTE.massf,fs_treshold,gradc,lcorr,neigh,pou,boxlimmin,boxlimmax,inn,rmat,track,simulate2d);
        if(symm && rsymp1)KerComputeNormalsBufferBox<false> (p1,posp1,CTE.massf,fs_treshold,gradc,lcorr,neigh,pou,boxlimmin,boxlimmax,inn,rmat,track,simulate2d); //<vs_syymmetry> 
      }

      unsigned fstypep1=0;
      if(simulate2d){
        if(fs_treshold<1.7) fstypep1=2;
        if(fs_treshold<1.1 && Nzero/float(neigh)<0.4f) fstypep1=3;
      } else {
        if(fs_treshold<2.75) fstypep1=2;
        if(fs_treshold<1.8 && Nzero/float(neigh)<0.4f) fstypep1=3;
      }
      fstype[p1]=fstypep1;

      //-Add the contribution of the particle itself
      pou+=cufsph::GetKernel_Wab<KERNEL_Wendland>(0.f)*CTE.massf/velrhop1.w;

      //-Calculation of the inverse of the correction matrix (Don't think there is a better way, create function for Determinant2x2 for clarity?).
      if(simulate2d){
        tmatrix2f lcorr2d;
        tmatrix2f lcorr2d_inv;
        lcorr2d.a11=lcorr.a11; lcorr2d.a12=lcorr.a13;
        lcorr2d.a21=lcorr.a31; lcorr2d.a22=lcorr.a33;
        float lcorr_det=(lcorr2d.a11*lcorr2d.a22-lcorr2d.a12*lcorr2d.a21);
        lcorr2d_inv.a11=lcorr2d.a22/lcorr_det; lcorr2d_inv.a12=-lcorr2d.a12/lcorr_det; lcorr2d_inv.a22=lcorr2d.a11/lcorr_det; lcorr2d_inv.a21=-lcorr2d.a21/lcorr_det;
        lcorr_inv.a11=lcorr2d_inv.a11;  lcorr_inv.a13=lcorr2d_inv.a12;
        lcorr_inv.a31=lcorr2d_inv.a21;  lcorr_inv.a33=lcorr2d_inv.a22;
      } else {
        const float determ = cumath::Determinant3x3(lcorr);
        lcorr_inv = cumath::InverseMatrix3x3(lcorr, determ);
      }

      //-Correction of the gradient of concentration and definition of the normal.
      float3 gradc1=make_float3(0,0,0);    
      gradc1.x=gradc.x*lcorr_inv.a11+gradc.y*lcorr_inv.a12+gradc.z*lcorr_inv.a13;
      gradc1.y=gradc.x*lcorr_inv.a21+gradc.y*lcorr_inv.a22+gradc.z*lcorr_inv.a23;
      gradc1.z=gradc.x*lcorr_inv.a31+gradc.y*lcorr_inv.a32+gradc.z*lcorr_inv.a33;    
      float gradc_norm=sqrt(gradc1.x*gradc1.x+gradc1.y*gradc1.y+gradc1.z*gradc1.z);
      fsnormal[p1].x=-gradc1.x/gradc_norm;
      fsnormal[p1].y=-gradc1.y/gradc_norm;
      fsnormal[p1].z=-gradc1.z/gradc_norm;
    }
  }

 

  __global__ void KerCountFreeSurface(unsigned n,unsigned pini
    ,unsigned* fs,unsigned* listp)
  {
    extern __shared__ unsigned slist[];
    //float* splanes=(float*)(slist+(n+1));
    if(!threadIdx.x)slist[0]=0;
    __syncthreads();
    const unsigned pp=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
    if(pp<n){
      const unsigned p=pp+pini;
      if(fs[p]>1) slist[atomicAdd(slist,1)+1]=p;
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

  void ComputeFSNormals(TpKernel tkernel,bool simulate2d,bool symmetry,unsigned bsfluid,unsigned fluidini,unsigned fluidnum
    ,StDivDataGpu& dvd,const unsigned* dcell,const double2* posxy,const double* posz
    ,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp,float4* shiftposfs
    ,unsigned* fstype,float3* fsnormal,unsigned* listp,StrGeomVresGpu* vresgdata,cudaStream_t stm)
  {
    unsigned count=0;

    //-Obtain the list of particle that are probably on the free-surface
    if(fluidnum){
      cudaMemset(listp+fluidnum,0,sizeof(unsigned));
      dim3 sgridf=GetSimpleGridSize(fluidnum,bsfluid);
      const unsigned smem=(bsfluid+1)*sizeof(unsigned); 
      KerCountFreeSurface <<<sgridf,bsfluid,smem,stm>>> (fluidnum,fluidini,fstype,listp);
    }
    cudaMemcpy(&count,listp+fluidnum,sizeof(unsigned),cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    if(count){
      dim3 sgridf=GetSimpleGridSize(count,bsfluid);
      if(symmetry) //<vs_syymmetry_ini>
        KerComputeNormals<true> <<<sgridf,bsfluid,0,stm>>> 
        (count,fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,dcell
            ,poscell,velrho,code,fstype,fsnormal,simulate2d,shiftposfs
            ,ftomassp,listp
            ,posxy,posz,vresgdata->boxdommin,vresgdata->boxdommax,vresgdata->inner,vresgdata->matmov,vresgdata->tracking);
      else //<vs_syymmetry_end>
        KerComputeNormals<false> <<<sgridf,bsfluid,0,stm>>> 
        (count,fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,dcell
            ,poscell,velrho,code,fstype,fsnormal,simulate2d,shiftposfs
            ,ftomassp,listp
            ,posxy,posz,vresgdata->boxdommin,vresgdata->boxdommax,vresgdata->inner,vresgdata->matmov,vresgdata->tracking);
    }
    
    cudaDeviceSynchronize();
    
  }

    template<bool symm>  __device__ void KerScanUmbrellaRegionBufferBox(unsigned p1
    ,const double3 posp1,bool& fs_flag,const float3* fsnormal
    ,const double3 boxlimitmin,const double3 boxlimitmax,const bool inner
    ,const tmatrix4f mat,const bool tracking,const bool sim2d)
  {    
    
    float3 minpos=make_float3(0,0,0);
    minpos.x=posp1.x-CTE.kernelsize;
    minpos.y=(sim2d? 0.0f: posp1.y-CTE.kernelsize);
    minpos.z=posp1.z-CTE.kernelsize;
    float3 maxpos=make_float3(0,0,0);
    maxpos.x=posp1.x+CTE.kernelsize;
    maxpos.y=(sim2d? 0.0f: posp1.y+CTE.kernelsize);
    maxpos.z=posp1.z+CTE.kernelsize;

    float dp=CTE.dp;
    for (float rx=minpos.x; rx<=maxpos.x; rx+=dp) for (float ry=minpos.y; ry<=maxpos.y; ry+=dp)
      for (float rz=minpos.z; rz<=maxpos.z; rz+=dp){
      const float drx=float(posp1.x-rx);
      const float dry=float(posp1.y-ry);
      const float drz=float(posp1.z-rz);
      const float rr2=drx*drx+dry*dry+drz*drz;
      float rx1=rx; float ry1=ry; float rz1=rz;
      if(tracking){          
        rx1=rx*mat.a11+ry*mat.a21+rz*mat.a31+mat.a14;
        ry1=rx*mat.a12+ry*mat.a22+rz*mat.a32+mat.a24;
        rz1=rx*mat.a13+ry*mat.a23+rz*mat.a33+mat.a34;
      }
      bool outside=(inner ? !KerBufferInZone(make_double2(rx1,ry1),rz1,boxlimitmin,boxlimitmax) : KerBufferInZone(make_double2(rx1,ry1),rz1,boxlimitmin,boxlimitmax));
      if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO && outside){

      const float3 posq=make_float3(fsnormal[p1].x*CTE.kernelh/* +posxy[p1].x */,fsnormal[p1].y*CTE.kernelh/* +posxy[p1].y */,fsnormal[p1].z*CTE.kernelh/* +posz[p1] */);

      if (rr2>2.f*CTE.kernelh*CTE.kernelh){
        // const float3 posp2=make_float3(posxy[p2].x,posxy[p2].y,posz[p2]);
        const float drxq=-drx-posq.x;
        const float dryq=-dry-posq.y;
        const float drzq=-drz-posq.z;
        const float rrq=sqrt(drxq*drxq+dryq*dryq+drzq*drzq);
        if(rrq<CTE.kernelh) fs_flag=true;
      } else {
        if(sim2d){
        const float drxq=-drx-posq.x;
        const float drzq=-drz-posq.z;
        const float3 normalq=make_float3(drxq*fsnormal[p1].x,0,drzq*fsnormal[p1].z);
        const float3 tangq=make_float3(-drxq*fsnormal[p1].z,0,drzq*fsnormal[p1].x);
        const float normalqnorm=sqrt(normalq.x*normalq.x+normalq.z*normalq.z);
        const float tangqnorm=sqrt(tangq.x*tangq.x+tangq.z*tangq.z);
        if (normalqnorm+tangqnorm<CTE.kernelh) fs_flag=true;
        } else{
        float rrr=1.f/sqrt(rr2);
        const float arccosin=acos((-drx*fsnormal[p1].x*rrr-dry*fsnormal[p1].y*rrr-drz*fsnormal[p1].z*rrr));
        if (arccosin < 0.785398) fs_flag=true;
        }
      }

    }
  }
  if(fs_flag) return;
  }

  template<bool symm>  __device__ void KerScanUmbrellaRegionBox(bool boundp2,unsigned p1
    ,const unsigned &pini,const unsigned &pfin,const float4 *poscell,const float4 &pscellp1
    ,bool& fs_flag,const float3* fsnormal,bool simulate2d)
  {
    for(int p2=pini;p2<pfin;p2++){
    const float4 pscellp2=poscell[p2];
      float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
      float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
      float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
      if(symm)dry=pscellp1.y+pscellp2.y + CTE.poscellsize*PSCEL_GetfY(pscellp2.w); //<vs_syymmetry>
      const double rr2=drx*drx+dry*dry+drz*drz;
      if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){

        const float3 posq=make_float3(fsnormal[p1].x*CTE.kernelh,fsnormal[p1].y*CTE.kernelh,fsnormal[p1].z*CTE.kernelh);

        if (rr2>2.f*CTE.kernelh*CTE.kernelh){
          const float drxq=-drx-posq.x;
          const float dryq=-dry-posq.y;
          const float drzq=-drz-posq.z;
          const float rrq=sqrt(drxq*drxq+dryq*dryq+drzq*drzq);
          if(rrq<CTE.kernelh) fs_flag=true;
        } else {
          if(simulate2d){
          const float drxq=-drx-posq.x;
          const float drzq=-drz-posq.z;
          const float3 normalq=make_float3(drxq*fsnormal[p1].x,0,drzq*fsnormal[p1].z);
          const float3 tangq=make_float3(-drxq*fsnormal[p1].z,0,drzq*fsnormal[p1].x);
          const float normalqnorm=sqrt(normalq.x*normalq.x+normalq.z*normalq.z);
          const float tangqnorm=sqrt(tangq.x*tangq.x+tangq.z*tangq.z);
          if (normalqnorm+tangqnorm<CTE.kernelh) fs_flag=true;
          } else{
            float rrr=1.f/sqrt(rr2);
          const float arccosin=acos((-drx*fsnormal[p1].x*rrr-dry*fsnormal[p1].y*rrr-drz*fsnormal[p1].z*rrr));
          if (arccosin < 0.785398) fs_flag=true;
          }
        }

      }
    }
    if(fs_flag) return;
  }

  template<bool symm>
  __global__ void KerScanUmbrellaRegion(unsigned n,unsigned pinit
    ,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid,const unsigned *dcell
    ,const float4 *poscell,const typecode* code,unsigned* fstype,float3* fsnormal,bool simulate2d,const unsigned* listp
    ,const double2* posxy,const double* posz,const double3* boxlimitmin,const double3* boxlimitmax,const bool* inner,const tmatrix4f* mat,const bool* tracking)
  {
    const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
    if(p<n){
      const unsigned p1=listp[p];   
      //-Obtains basic data of particle p1.
      const float4 pscellp1=poscell[p1];
      
      const bool rsymp1=(symm && PSCEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>

      bool fs_flag=false;
      int ini1,fin1,ini2,fin2,ini3,fin3;

      cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      //-Interaction with fluids.
      ini3+=cellfluid; fin3+=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin){
          KerScanUmbrellaRegionBox<false> (false,p1,pini,pfin,poscell,pscellp1,fs_flag,fsnormal,simulate2d);
          if(symm && rsymp1)KerScanUmbrellaRegionBox<true > (false,p1,pini,pfin,poscell,pscellp1,fs_flag,fsnormal,simulate2d); //<vs_syymmetry>
        }
      }

      ini3-=cellfluid; fin3-=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin){
          KerScanUmbrellaRegionBox<false> (true,p1,pini,pfin,poscell,pscellp1,fs_flag,fsnormal,simulate2d);
          if(symm && rsymp1)KerScanUmbrellaRegionBox<true > (true,p1,pini,pfin,poscell,pscellp1,fs_flag,fsnormal,simulate2d); //<vs_syymmetry>
        }
      }

      // -Interaction with virtual stencil.
      if(CODE_IsFluidBuffer(code[p1])){
        const byte izone0=byte(CODE_GetIzoneFluidBuffer(code[p1]));
        const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
        const double3   boxlimmin =boxlimitmin[izone];
        const double3   boxlimmax =boxlimitmax[izone];
        const bool      inn       =inner      [izone];
        const tmatrix4f rmat      =mat        [izone];
        const bool      track     =tracking   [izone]; 
        const double3   posp1     =make_double3(posxy[p1].x,posxy[p1].y,posz[p1]); 
        KerScanUmbrellaRegionBufferBox<false> (p1,posp1,fs_flag,fsnormal,boxlimmin,boxlimmax,inn,rmat,track,simulate2d);
        if(symm && rsymp1)KerScanUmbrellaRegionBufferBox<false> (p1,posp1,fs_flag,fsnormal,boxlimmin,boxlimmax,inn,rmat,track,simulate2d); //<vs_syymmetry> 
      }

      
      //-If particle was present in umbrella region, change the code of the particle.
      if(fs_flag && fstype[p1]==2) fstype[p1]=0;
      //-Periodic particle are internal by default.
      if(CODE_IsPeriodic(code[p1])) fstype[p1]=0;

    
    }
  }






  void ComputeUmbrellaRegion(TpKernel tkernel,bool simulate2d,bool symmetry,unsigned bsfluid,unsigned fluidini,unsigned fluidnum
    ,StDivDataGpu& dvd,const unsigned* dcell,const double2* posxy,const double* posz
    ,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp,float4* shiftposfs
    ,unsigned* fstype,float3* fsnormal,unsigned* listp,StrGeomVresGpu* vresgdata,cudaStream_t stm)
  {
    unsigned count=0;

    //-Obtain the list of particle that are probably on the free-surface (in ComputeUmbrellaRegion maybe is unnecessary).
    if(fluidnum){
      cudaMemset(listp+fluidnum,0,sizeof(unsigned));
      dim3 sgridf=GetSimpleGridSize(fluidnum,bsfluid);
      const unsigned smem=(bsfluid+1)*sizeof(unsigned); 
      KerCountFreeSurface <<<sgridf,bsfluid,smem,stm>>> (fluidnum,fluidini,fstype,listp);
    }
    cudaMemcpy(&count,listp+fluidnum,sizeof(unsigned),cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();

    if(count){
      dim3 sgridf=GetSimpleGridSize(count,bsfluid);
      if(symmetry) //<vs_syymmetry_ini>
        KerScanUmbrellaRegion<true> <<<sgridf,bsfluid,0,stm>>> 
        (count,fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,dcell
            ,poscell,code,fstype,fsnormal,simulate2d,listp
            ,posxy,posz,vresgdata->boxdommin,vresgdata->boxdommax,vresgdata->inner,vresgdata->matmov,vresgdata->tracking);
      else //<vs_syymmetry_end>
        KerScanUmbrellaRegion<false> <<<sgridf,bsfluid,0,stm>>> 
        (count,fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid,dcell
            ,poscell,code,fstype,fsnormal,simulate2d,listp
            ,posxy,posz,vresgdata->boxdommin,vresgdata->boxdommax,vresgdata->inner,vresgdata->matmov,vresgdata->tracking);
    }
    
    cudaDeviceSynchronize();
    
  }

    template<bool symm>  __device__ void KerComputeShiftingVelBufferBox(unsigned p1,const double3 posp1
    ,float massp2,float4& shiftvel,const double3 boxlimitmin,const double3 boxlimitmax
    ,const bool inner,const tmatrix4f mat,const bool tracking,const bool sim2d)
  {    
      
      float3 minpos=make_float3(0,0,0);
      minpos.x=posp1.x-CTE.kernelsize;
      minpos.y=(sim2d? 0.0f: posp1.y-CTE.kernelsize);
      minpos.z=posp1.z-CTE.kernelsize;
      float3 maxpos=make_float3(0,0,0);
      maxpos.x=posp1.x+CTE.kernelsize;
      maxpos.y=(sim2d? 0.0f: posp1.y+CTE.kernelsize);
      maxpos.z=posp1.z+CTE.kernelsize;

      float dp=CTE.dp;
      for (float rx=minpos.x; rx<=maxpos.x; rx+=dp) for (float ry=minpos.y; ry<=maxpos.y; ry+=dp)
        for (float rz=minpos.z; rz<=maxpos.z; rz+=dp){
          const float drx=float(posp1.x-rx);
          const float dry=float(posp1.y-ry);
          const float drz=float(posp1.z-rz);
          const float rr2=drx*drx+dry*dry+drz*drz;
          float rx1=rx; float ry1=ry; float rz1=rz;
          if(tracking){          
            rx1=rx*mat.a11+ry*mat.a21+rz*mat.a31+mat.a14;
            ry1=rx*mat.a12+ry*mat.a22+rz*mat.a32+mat.a24;
            rz1=rx*mat.a13+ry*mat.a23+rz*mat.a33+mat.a34;
          }
          bool outside=(inner ? !KerBufferInZone(make_double2(rx1,ry1),rz1,boxlimitmin,boxlimitmax) : KerBufferInZone(make_double2(rx1,ry1),rz1,boxlimitmin,boxlimitmax));
          if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO && outside){
            //-Computes kernel.
            const float fac=cufsph::GetKernel_Fac<KERNEL_Wendland>(rr2);
            const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
            
            const float vol2=massp2/CTE.rhopzero;

            shiftvel.x+=vol2*frx;    
            shiftvel.y+=vol2*fry;
            shiftvel.z+=vol2*frz;          
            const float wab=cufsph::GetKernel_Wab<KERNEL_Wendland>(rr2);
            shiftvel.w+=wab*vol2;      
      }
    }
  }



  template<TpKernel tker,bool simulate2d,bool shiftadv,bool symm>
  __device__ void KerPreLoopInteractionBox(bool boundp2,unsigned p1
    ,const unsigned &pini,const unsigned &pfin,const float4 *poscell,const float4 *velrhop
    ,const typecode *code,float massp2,const float4 &pscellp1,const float4 &velrhop1,const float* ftomassp
    ,float4 &shiftposf1,unsigned* fs,float3* fsnormal,bool& nearfs,float &mindist,float& maxarccos
    ,bool& bound_inter,float3& fsnormalp1,float& pou)
  {
    for(int p2=pini;p2<pfin;p2++){
      const float4 pscellp2=poscell[p2];
      float drx=pscellp1.x-pscellp2.x + CTE.poscellsize*(PSCEL_GetfX(pscellp1.w)-PSCEL_GetfX(pscellp2.w));
      float dry=pscellp1.y-pscellp2.y + CTE.poscellsize*(PSCEL_GetfY(pscellp1.w)-PSCEL_GetfY(pscellp2.w));
      float drz=pscellp1.z-pscellp2.z + CTE.poscellsize*(PSCEL_GetfZ(pscellp1.w)-PSCEL_GetfZ(pscellp2.w));
      if(symm)dry=pscellp1.y+pscellp2.y + CTE.poscellsize*PSCEL_GetfY(pscellp2.w); //<vs_syymmetry>
      const double rr2=drx*drx+dry*dry+drz*drz;
      if(rr2<=CTE.kernelsize2 && rr2>=ALMOSTZERO){

        const float fac=cufsph::GetKernel_Fac<KERNEL_Wendland>(rr2);
        const float frx=fac*drx,fry=fac*dry,frz=fac*drz; //-Gradients.
        float4 velrhop2=velrhop[p2];

        bool ftp2;
        float ftmassp2;    //-Contains mass of floating body or massf if fluid. | Contiene masa de particula floating o massp2 si es bound o fluid.
        const typecode cod=code[p2];
        ftp2=CODE_IsFloating(cod);
        ftmassp2=(ftp2? ftomassp[CODE_GetTypeValue(cod)]: massp2);
        
        if(shiftadv){
          const float massrho=(boundp2 ? CTE.massb/velrhop2.w : (ftmassp2)/velrhop2.w);

          //-Compute gradient of concentration and partition of unity.        
          shiftposf1.x+=massrho*frx;    
          shiftposf1.y+=massrho*fry;
          shiftposf1.z+=massrho*frz;          
          const float wab=cufsph::GetKernel_Wab<KERNEL_Wendland>(rr2);
          shiftposf1.w+=wab*massrho;


          //-Check if the particle is too close to solid or floating object
          if((boundp2 || ftp2)) bound_inter=true;

          //-Check if it close to the free-surface, calculate distance from free-surface and smoothing of free-surface normals.
          if(fs[p2]>1 && fs[p2]<3 && !boundp2 ) {
            nearfs=true;
            mindist=min(sqrt(rr2),mindist);
            pou+=wab*massrho;

            fsnormalp1.x+=fsnormal[p2].x*wab*massrho;
            fsnormalp1.y+=fsnormal[p2].y*wab*massrho;
            fsnormalp1.z+=fsnormal[p2].z*wab*massrho;
                        
          }

          //-Check maximum curvature.
          if(fs[p1]>1 && fs[p2]>1){
          const float norm1=sqrt(fsnormal[p1].x*fsnormal[p1].x+fsnormal[p1].y*fsnormal[p1].y+fsnormal[p1].z*fsnormal[p1].z);
          const float norm2=sqrt(fsnormal[p2].x*fsnormal[p2].x+fsnormal[p2].y*fsnormal[p2].y+fsnormal[p2].z*fsnormal[p2].z);
          maxarccos=max(maxarccos,(acos((fsnormal[p1].x*fsnormal[p2].x+fsnormal[p2].y*fsnormal[p1].y+fsnormal[p2].z*fsnormal[p1].z))));
          }
        }
      }
    }
  }

  template<TpKernel tker,bool simulate2d,bool shiftadv,bool symm>
  __global__ void KerPreLoopInteraction(unsigned n
    ,unsigned pinit,int scelldiv,int4 nc,int3 cellzero,const int2 *begincell,unsigned cellfluid
    ,const unsigned *dcell,const float4 *poscell,const float4 *velrhop,const typecode *code
    ,const float* ftomassp,float4* shiftvel,unsigned* fstype,float3* fsnormal,float* fsmindist
    ,const double2* posxy,const double* posz,const double3* boxlimitmin,const double3* boxlimitmax
    ,const bool* inner,const tmatrix4f* mat,const bool* tracking)
  {
    const unsigned p=blockIdx.x*blockDim.x + threadIdx.x; //-Number of particle.
    if(p<n){
      const unsigned p1=p+pinit;      //-Number of particle.
      //-Obtains basic data of particle p1.
      const float4 pscellp1=poscell[p1];
      const float4 velrhop1=velrhop[p1];
      const float pressp1=cufsph::ComputePressCte(velrhop1.w);
      const bool rsymp1=(symm && PSCEL_GetPartY(__float_as_uint(pscellp1.w))==0); //<vs_syymmetry>
      bool    nearfs=false;                     //-Bool for detecting near free-surface particles. <shiftImproved>
      float4  shiftposp1=make_float4(0,0,0,0);

      
      float mindist=CTE.kernelh;                //-Set Min Distance from free-surface to kernel radius. <shiftImproved>
      float maxarccos=0.0;                      //-Variable for identify high-curvature free-surface particle <shiftImproved>
      bool bound_inter=false;                   //-Variable for identify free-surface that interact with boundary <shiftImproved>
      float3 fsnormalp1=make_float3(0,0,0);     //-Normals for near free-surface particles <shiftImproved>
      unsigned fsp1=fstype[p1];                 //-Free-surface identification code: 0-internal, 1-close to free-surface, 2 free-surface, 3-isolated.
      float   pou=false;                        //-Partition of unity for normal correction.                      <ShiftingAdvanced>
    
      //-Obtains neighborhood search limits.
      int ini1,fin1,ini2,fin2,ini3,fin3;
      cunsearch::InitCte(dcell[p1],scelldiv,nc,cellzero,ini1,fin1,ini2,fin2,ini3,fin3);

      //-Interaction with fluids.
      ini3+=cellfluid; fin3+=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin){
          KerPreLoopInteractionBox<tker,simulate2d,shiftadv,false> (false,p1,pini,pfin,poscell,velrhop
            ,code,CTE.massf,pscellp1,velrhop1,ftomassp,shiftposp1
            ,fstype,fsnormal,nearfs,mindist,maxarccos,bound_inter,fsnormalp1,pou);
          if(symm && rsymp1)KerPreLoopInteractionBox<tker,simulate2d,shiftadv,true> (false,p1,pini,pfin,poscell,velrhop
            ,code,CTE.massf,pscellp1,velrhop1,ftomassp,shiftposp1
            ,fstype,fsnormal,nearfs,mindist,maxarccos,bound_inter,fsnormalp1,pou); //<vs_syymmetry>
        } 
      }


      //-Interaction with bound.
      ini3-=cellfluid; fin3-=cellfluid;
      for(int c3=ini3;c3<fin3;c3+=nc.w)for(int c2=ini2;c2<fin2;c2+=nc.x){
        unsigned pini,pfin=0;  cunsearch::ParticleRange(c2,c3,ini1,fin1,begincell,pini,pfin);
        if(pfin){
        KerPreLoopInteractionBox<tker,simulate2d,shiftadv,false> (true,p1,pini,pfin,poscell,velrhop
            ,code,CTE.massf,pscellp1,velrhop1,ftomassp,shiftposp1
            ,fstype,fsnormal,nearfs,mindist,maxarccos,bound_inter,fsnormalp1,pou);
        if(symm && rsymp1)KerPreLoopInteractionBox<tker,simulate2d,shiftadv,true> (true,p1,pini,pfin,poscell,velrhop
            ,code,CTE.massf,pscellp1,velrhop1,ftomassp,shiftposp1
            ,fstype,fsnormal,nearfs,mindist,maxarccos,bound_inter,fsnormalp1,pou); //<vs_syymmetry>
      }
      }



      if(shiftadv){

        // -Interaction with virtual stencil.
        if(CODE_IsFluidBuffer(code[p1])){
          const byte izone0=byte(CODE_GetIzoneFluidBuffer(code[p1]));
          const byte izone=(izone0&CODE_TYPE_FLUID_INOUT015MASK);
          const double3   boxlimmin =boxlimitmin[izone];
          const double3   boxlimmax =boxlimitmax[izone];
          const bool      inn       =inner      [izone];
          const tmatrix4f rmat      =mat        [izone];
          const bool      track     =tracking   [izone]; 
          const double3   posp1     =make_double3(posxy[p1].x,posxy[p1].y,posz[p1]); 
          KerComputeShiftingVelBufferBox<false> (p1,posp1,CTE.massf,shiftposp1,boxlimmin,boxlimmax,inn,rmat,track,simulate2d);
          if(symm && rsymp1)KerComputeShiftingVelBufferBox<false> (p1,posp1,CTE.massf,shiftposp1,boxlimmin,boxlimmax,inn,rmat,track,simulate2d); //<vs_syymmetry> 
        }

        shiftposp1.w+=cufsph::GetKernel_Wab<KERNEL_Wendland>(0.0)*CTE.massf/velrhop1.w;


      fsmindist[p1]=mindist;
      //-Assign correct code to near free-surface particle and correct their normals by Shepard's Correction.
      if(fsp1==0 && nearfs){
      if(pou>1e-6){
        fsnormalp1=make_float3(fsnormalp1.x,fsnormalp1.y,fsnormalp1.z);
        float norm=sqrt(fsnormalp1.x*fsnormalp1.x+fsnormalp1.y*fsnormalp1.y+fsnormalp1.z*fsnormalp1.z);
        fsnormal[p1]=make_float3(fsnormalp1.x/norm,fsnormalp1.y/norm,fsnormalp1.z/norm);
      }      
      fstype[p1]=1;
      if(bound_inter) fstype[p1]=3;
      }
      
      //-Check if free-surface particle interact with bound or has high-curvature.
      if(fsp1==2 && (bound_inter||maxarccos>0.52)) fstype[p1]=3;


      //-Compute shifting when <shiftImproved> true
      shiftvel[p1]=shiftposp1;
      }
      
    }
  }

  
  
  
  
  template<TpKernel tker,bool simulate2d,bool shiftadv> void PreLoopInteractionT3(bool symmetry
    ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd,const double2* posxy,const double* posz
    ,const unsigned* dcell,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp
    ,float4* shiftvel,unsigned* fstype,float3* fsnormal,float* fsmindist,StrGeomVresGpu* vresgdata,cudaStream_t stm)
{

  if(fluidnum){
    dim3 sgridf=GetSimpleGridSize(fluidnum,bsfluid);
    if(symmetry) //<vs_syymmetry_ini>
    KerPreLoopInteraction <tker,simulate2d,shiftadv,true> <<<sgridf,bsfluid,0,stm>>> 
      (fluidnum,fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid
      ,dcell,poscell,velrho,code,ftomassp,shiftvel,fstype,fsnormal,fsmindist
      ,posxy,posz,vresgdata->boxdommin,vresgdata->boxdommax,vresgdata->inner,vresgdata->matmov,vresgdata->tracking);
    else //<vs_syymmetry_end>
    KerPreLoopInteraction <tker,simulate2d,shiftadv,false> <<<sgridf,bsfluid,0,stm>>> 
      (fluidnum,fluidini,dvd.scelldiv,dvd.nc,dvd.cellzero,dvd.beginendcell,dvd.cellfluid
      ,dcell,poscell,velrho,code,ftomassp,shiftvel,fstype,fsnormal,fsmindist
      ,posxy,posz,vresgdata->boxdommin,vresgdata->boxdommax,vresgdata->inner,vresgdata->matmov,vresgdata->tracking);

  }
}

  template<TpKernel tker,bool simulate2d> void PreLoopInteractionT2(bool shiftadv,bool symmetry
    ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd,const double2* posxy,const double* posz
    ,const unsigned* dcell,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp
    ,float4* shiftvel,unsigned* fstype,float3* fsnormal,float* fsmindist,StrGeomVresGpu* vresgdata,cudaStream_t stm)
{
  if(shiftadv){
    PreLoopInteractionT3 <tker,simulate2d,true > (symmetry,bsfluid
        ,fluidnum,fluidini,dvd,posxy,posz,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,vresgdata,stm);
  }
  else{
    PreLoopInteractionT3 <tker,simulate2d,false> (symmetry,bsfluid
        ,fluidnum,fluidini,dvd,posxy,posz,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,vresgdata,stm);
  }
}

  template<TpKernel tker> void PreLoopInteractionT(bool simulate2d,bool shiftadv,bool symmetry
    ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd,const double2* posxy,const double* posz
    ,const unsigned* dcell,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp
    ,float4* shiftvel,unsigned* fstype,float3* fsnormal,float* fsmindist,StrGeomVresGpu* vresgdata,cudaStream_t stm)
{
  if(simulate2d){
    PreLoopInteractionT2 <tker,true > (shiftadv,symmetry,bsfluid
        ,fluidnum,fluidini,dvd,posxy,posz,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,vresgdata,stm);
  }
  else{
    PreLoopInteractionT2 <tker,false> (shiftadv,symmetry,bsfluid
        ,fluidnum,fluidini,dvd,posxy,posz,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,vresgdata,stm);
  }
}


   void PreLoopInteraction(TpKernel tkernel,bool simulate2d,bool shiftadv,bool symmetry
    ,unsigned bsfluid,unsigned fluidnum,unsigned fluidini,StDivDataGpu& dvd,const double2* posxy,const double* posz
    ,const unsigned* dcell,const float4* poscell,const float4* velrho,const typecode* code,const float* ftomassp
    ,float4* shiftvel,unsigned* fstype,float3* fsnormal,float* fsmindist,StrGeomVresGpu* vresgdata,cudaStream_t stm)
  {
    switch(tkernel){
    case KERNEL_Wendland:{ const TpKernel tker=KERNEL_Wendland;
      PreLoopInteractionT <tker> (simulate2d,shiftadv,symmetry,bsfluid
        ,fluidnum,fluidini,dvd,posxy,posz,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,vresgdata,stm);
    }break;
  #ifndef DISABLE_KERNELS_EXTRA
    case KERNEL_Cubic:{ const TpKernel tker=KERNEL_Cubic;
      PreLoopInteractionT <tker> (simulate2d,shiftadv,symmetry,bsfluid
        ,fluidnum,fluidini,dvd,posxy,posz,dcell,poscell,velrho,code,ftomassp
        ,shiftvel,fstype,fsnormal,fsmindist,vresgdata,stm);
    }break;
  #endif
    default: throw "Kernel unknown at Interaction_MdbcCorrection().";
  }
  
    cudaDeviceSynchronize();
    
  }

}
