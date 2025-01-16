#ifndef FUNCCOMP
#define FUNCCOMP

#include <iostream>
#include <fstream>

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

double CalculoEnergiaCinetica(int nd,uint np,double *v)
{
    double ec=0.0;
    for(int ip=0; ip<np; ip++)
        for(int id=0; id<nd; id++)
            ec += v[id+ip*nd] * v[id+ip*nd];    

    return(ec/np); 
}

void Velocidades(int np,int nd,double *v,double *a,double dt)
{
    for(int ip=0; ip<np;ip++)
        for(int id=0;id<nd;id++)
            v[id+ip*nd] += a[id+ip*nd] * 0.5 * dt;
}

void IAaD(int np,std::ofstream &ofasat,double *p,double *v,double *a,int nd)
{    
    int ip;
    ofasat << std::endl;
    for(ip=0; ip<np; ip++)
    {
        ofasat << ip << " " <<
        p[ip*nd] << " " << p[1+ip*nd] << " " << p[2+ip*nd] << " " <<
        v[ip*nd] << " " << v[1+ip*nd] << " " << v[2+ip*nd] << " " <<
        a[ip*nd] << " " << a[1+ip*nd] << " " << a[2+ip*nd] << std::endl;
    }
}

void ArchivosDeResultados(str dpsco,std::ofstream &ofasres,std::ofstream &ofasat,str op)
{
    str asat = dpsco + "/Posiciones"+op+".txt";
    str asres = dpsco + "/Res"+op+".txt";

    ofasat.open(asat.c_str());
    ofasres.open(asres.c_str());
}
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

__global__ void KernelPosiciones(uint np,double *p,double *v,double *a, double dt,double3 caja,double3 cajai,int3 condper)
{
    int ip=threadIdx.x+blockDim.x*blockIdx.x;
    double pix(0.0),piy(0.0),piz(0.0),vix(0.0),viy(0.0),viz(0.0),aix(0.0),aiy(0.0),aiz(0.0);
    pix=__ldg(p+3*ip);
    piy=__ldg(p+3*ip+1);
    piz=__ldg(p+3*ip+2);
    vix=__ldg(v+3*ip);
    viy=__ldg(v+3*ip+1);
    viz=__ldg(v+3*ip+2);
    aix=__ldg(a+3*ip);
    aiy=__ldg(a+3*ip+1);
    aiz=__ldg(a+3*ip+2);
    pix+=vix*dt+aix*dt*dt*0.5;
    piy+=viy*dt+aiy*dt*dt*0.5;
    piz+=viz*dt+aiz*dt*dt*0.5;
    
    p[ip*3]-=(condper.x)?rint(p[ip*3]*cajai.x)*caja.x:0;
    p[ip*3+1]-=(condper.y)?rint(p[ip*3+1]*cajai.y)*caja.y:0;
    p[ip*3+2]-=(condper.z)?rint(p[ip*3+2]*cajai.z)*caja.z:0;
}

__global__ void KernelVelocidades(uint np,double *v,double *a, double dt,double3 caja,double3 cajai,int3 condper)
{
    int ip=threadIdx.x+blockDim.x*blockIdx.x;
    double vix(0.0),viy(0.0),viz(0.0),aix(0.0),aiy(0.0),aiz(0.0);
    vix=__ldg(v+3*ip);
    viy=__ldg(v+3*ip+1);
    viz=__ldg(v+3*ip+2);
    aix=__ldg(a+3*ip);
    aiy=__ldg(a+3*ip+1);
    aiz=__ldg(a+3*ip+2);
    vix+=aix*dt;
    viy+=aiy*dt;
    viz+=aiz*dt;
    v[ip*3]=vix;
    v[ip*3+1]=viy;
    v[ip*3+2]=viz;
}


#endif