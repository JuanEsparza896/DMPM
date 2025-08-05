#ifndef FUNCCOMP
#define FUNCCOMP

#include <iostream>
#include <fstream>

//para calculo de la energia
template<int Numthreads>
__global__ void Reduccionconwarps(double *arr,int nelementos,bool energia)
{
    __shared__ double s_mem[Numthreads];
    //Solo usamos un bloque
    int gid=threadIdx.x;
    int lane =threadIdx.x&31;
    int warp=threadIdx.x/32;
    int faltan(nelementos);
    double suma(0.0);
    int n=nelementos/Numthreads;
    n+=(nelementos%Numthreads)?1:0;
    for(int i=0;i<n;i++)suma+=(n*gid+i<nelementos)?(energia?arr[n*gid+i]:(arr[n*gid+i]*arr[n*gid+i])):0;
    __syncthreads();
    while(faltan>0){
        for(int i=16;i>=1;i/=2)suma+=__shfl_down_sync(0xffffffff,suma,i);
        faltan/=32;
        if(!lane)s_mem[warp]=suma;
        __syncthreads();
        suma=s_mem[gid];
    }
    suma/=nelementos;
    if(!gid)arr[0]=suma;
}

/*******************************************************************************************************/

//ARREGLO DE POSICIONES

/*******************************************************************************************************/

inline __device__ __host__ double3 CondPeriodicas(int3 condper,double3 caja,double3 dx,double3 cajai)
{
    //rint y signo hacen lo mismo pero rint es mas rapido
    //https://docs.nvidia.com/cuda/cuda-math-api/cuda_math_api/group__CUDA__MATH__DOUBLE.html#group__cuda__math__double_1ga3b8026edb2f2e441669845f0f3fa3bf7
    
    if(condper.x)dx.x-=rint(dx.x*cajai.x)*caja.x;
    if(condper.y)dx.y-=rint(dx.y*cajai.y)*caja.y;
    if(condper.z)dx.z-=rint(dx.z*cajai.z)*caja.z;

    return dx;
}

__host__ __device__ inline double Discuad(double3 var)
{
    return (var.x*var.x+var.y*var.y+var.z*var.z);
}

/*******************************************************************************************************/

//ARREGLO DE VELOCIDADES

/*******************************************************************************************************/

double CalculoEnergiaCinetica(uint np,double *v)
{
    double ec=0.0;
    for(int ip=0; ip<np; ip++)
        for(int id=0; id<nd; id++)
            ec += v[id+ip*nd] * v[id+ip*nd];    

    return(ec/np); 
}

void Velocidades(int np,double *v,double *a,double dt)
{
    for(int ip=0; ip<np;ip++)
        for(int id=0;id<nd;id++)
            v[id+ip*nd] += a[id+ip*nd] * 0.5 * dt;
}

/*******************************************************************************************************/

#endif
