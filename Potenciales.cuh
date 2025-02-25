#ifndef HEADER_POT
#define HEADER_POT
#include "Definiciones.cuh"

void CalculoParamLJ(uint n_esp_p,uint nparam,double *param,double *M_param)
{
    //u(r)=4eps*((sig/r)12-(sig/r)6)
    double sig,d2,d6,eps;
    for(int i=0;i<n_esp_p;i++)
        for(int j=0;j<n_esp_p;j++){
            sig=0.5*(param[i*nparam]+param[j*nparam]);
            d2=sig*sig;
            d6=d2*d2*d2;
            eps=sqrt(param[i*nparam+1]*param[j*nparam+1]);
            M_param[n_esp_p*i+j]=4.0*eps*d6;
            M_param[n_esp_p*n_esp_p+n_esp_p*i+j]=M_param[n_esp_p*i+j]*d6;
        }
}

__device__ __forceinline__ double2 InteraccionLJ(int i, int j,uint n_esp_p,uint elem_M,double dis, double *M_param,bool nconf,bool eshift,double r2c=0.0)
{
    double apot=0.;
    double2 fuepotshift,val;
    bool gidj=i-j;

    double d2=gidj?(1./dis):0.;
    double d6=d2*d2*d2;
    double d12=d6*d6;
    double fue=6.0*(M_param[n_esp_p*n_esp_p + elem_M]*2.0*d12-M_param[elem_M]*d6)*d2;
    if(nconf){
        apot=M_param[n_esp_p*n_esp_p + elem_M]*d12-M_param[elem_M]*d6;
    }
    
    val.x=fue;
    val.y=apot;
    if(eshift&&gidj){
        d2=1/r2c;
        d6=d2*d2*d2;
        d12=d6*d6;
        fuepotshift.x=6.0*(M_param[n_esp_p*n_esp_p + elem_M]*2.0*d12-M_param[elem_M]*d6)*d2;
        fuepotshift.y=M_param[n_esp_p*n_esp_p + elem_M]*d12-M_param[elem_M]*d6;
        val.x+=fuepotshift.x;
        val.y+=fuepotshift.y;
    }
    return val;
}

void InicializarMatrizDeParametros(uint pot,uint n_esp_p,uint nparam,double *param,double *M_param)
{
    switch (pot)
    {
    case LennardJones:
        CalculoParamLJ(n_esp_p,nparam,param,M_param);
        break;
    
    case Yukawa:
        break;
    }
}

__device__ __forceinline__ double2 Interaccion(uint pot,uint n_esp_p,uint elem_M,int i, int j,double dis,double *M_param,bool nconf,bool eshift, double r2c=0.)
{
    double2 fuepot;
    fuepot.x=0.;
    fuepot.y=0.;
    switch (pot)
    {
    case LennardJones:
        fuepot = InteraccionLJ(i,j,n_esp_p,elem_M,dis,M_param,nconf,eshift,r2c);
        break;
    
    case Yukawa:
        break;
    }
    return fuepot;
}

#endif
