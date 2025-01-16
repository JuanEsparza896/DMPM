#ifndef SIMNO_HEADER
#define SIMNO_HEADER

#include <stdio.h>

#include "OperacionesTDatosCuda.cuh"
#include "OperacionesDeHilosyBloques.hpp"
#include "Definiciones.cuh"
#include "Potenciales.h"
#include "Funcionescompartidas.h"

__global__ void AceleracionesfFuerzasLJ(uint np,const double *p,int chp,double *param,double3 caja,double3 cajai,int3 condper,double *epot,double *a,bool nconf)
{
    int gid=threadIdx.x+blockDim.x*blockIdx.x;
    gid/=chp;
    int lane=threadIdx.x &(chp-1);
    //los param van a ir en s_mem
    extern __shared__ double s_mem[];

    s_mem[0]=param[0];
    s_mem[1]=param[1];
    __syncthreads();
    double dis,pot=0.0;
    double2 fuepot;
    double3 fuerza,dif;
    fuerza=InitDataType3<double3>(0.0,0.0,0.0);    
    if(gid>=np)return;
    double pix,piy,piz,pjx,pjy,pjz;
    double3 pi,pj;
    pix=__ldg(p+3*gid);
    piy=__ldg(p+3*gid+1);
    piz=__ldg(p+3*gid+2);
    pi=InitDataType3<double3>(pix,piy,piz);
    #pragma unroll
    for(int j=lane;j<np;j+=chp){

            pjx=__ldg(p+3*j);
            pjy=__ldg(p+3*j+1);
            pjz=__ldg(p+3*j+2);
            pj=InitDataType3<double3>(pjx,pjy,pjz);
            dif.x=pi.x-pj.x;
            dif.y=pi.y-pj.y;
            dif.z=pi.z-pj.z;
            dif=CondPeriodicas(condper,caja,dif,cajai);
            dis=Discuad(dif);
            fuepot=InteraccionLJ(gid,j,dis,s_mem,nconf,false);
            fuerza.x+=fuepot.x*dif.x;
            fuerza.y+=fuepot.x*dif.y;
            fuerza.z+=fuepot.x*dif.z;
            if(nconf){
                pot+=fuepot.y;
            } 
    }
    
    //reduccion de warps aplicado a fuerza y al potencial
    __syncwarp();
    for(int i=chp/2.0;i>=1;i/=2.0){
        fuerza.x+=__shfl_down_sync(FULL_MASK,fuerza.x,i,chp);
        fuerza.y+=__shfl_down_sync(FULL_MASK,fuerza.y,i,chp);
        fuerza.z+=__shfl_down_sync(FULL_MASK,fuerza.z,i,chp);
        if(nconf)pot+=__shfl_down_sync(FULL_MASK,pot,i,chp);
    }
    __syncwarp();
    
    
    if(!lane){
        if(nconf)epot[gid]=pot;
        a[3*gid]=fuerza.x;
        a[3*gid+1]=fuerza.y;
        a[3*gid+2]=fuerza.z;
    }
    
}

void Simulacion(uint np,int nd,double *p,double *v,double *a,double sig,double eps,
                double3 caja,double3 cajai,int3 condper,double temp,std::ofstream &ofasres,
                std::ofstream &ofasat,uint nc,double dt,double dens,uint ncp,int nhilos,int pot,int maxhilos,size_t memoria_global)
{
    int ncc=nc*ncp/100.0;
    clock_t ti, tf;
    double dtt;

    double ets=0.0,etp=0.0,ett=0.0;
    double ecs=0.0,ecp=0.0,ect=0.0;
    double eis=0.0,eip=0.0,eit=0.0;
    double temps;

    int ns,id,ip,ic=0;
    bool nconf=true;
    
    /***************************************************************************************************************************/
    
    int chp=HilosPorParticula(nhilos);

    int threadsperblock,blockspergrid;
    printf("\n\nPara la rutina de Aceleracion:\n");
    HilosenBloqueMultiplodeWarp(np,chp,blockspergrid,threadsperblock,MINBLOCKPERGRID,maxhilos);

    /***************************************************************************************************************************/
    //para el potencial de Lennard-Jones tenemos 2 parametros pl6 y plj12
    //para otros potenciales podemos llegar a tener mas parametros
    
    int nparam;
    switch(pot){
        case 1:
                nparam=2;
                break;
    }

    /***************************************************************************************************************************/
    //Memoria para arreglos en la CPU
    
    double *h_param=new double [nparam];
    double *h_eit=new double[1];

    /***************************************************************************************************************************/
    //Memoria para arreglos en la GPU

    double *d_eit,*d_p,*d_a,*param;
    //double *d_v;
    Errorcuda(cudaMalloc(&param,nparam*sizeof(double)),"param",0);
    Errorcuda(cudaMalloc(&d_eit,sizeof(double)*np),"d_eit",0);
    Errorcuda(cudaMalloc(&d_p,sizeof(double)*nd*np),"d_p",0);
    //Errorcuda(cudaMalloc(&d_v,sizeof(double)*nd*np),"d_p",0);
    Errorcuda(cudaMalloc(&d_a,sizeof(double)*nd*np),"d_a",0);
    //vemos cuantos byte ocupan nuestros arreglos en la gpu y el porcentaje

    size_t memoria_global_utilizada=(nparam+np*(1+2*nd))*sizeof(double);
    OcupacionDeMemoriaGlobal(memoria_global_utilizada,memoria_global);

    /***************************************************************************************************************************/
    //Mas optimizaciones de cuda

    /***************************************************************************************************************************/
    
    CalculoParamLJ(sig,eps,h_param);
    
    Errorcuda(cudaMemcpy(param,h_param,sizeof(double)*nparam,cudaMemcpyHostToDevice),"param",1);
    Errorcuda(cudaMemcpy(d_p,p,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"p",1);

    AceleracionesfFuerzasLJ<<<blockspergrid,threadsperblock>>>(np,d_p,chp,param,caja,cajai,condper,d_eit,d_a,nconf);
    Errorcuda(cudaGetLastError(),"PLJ",3);
    Errorcuda(cudaMemcpy(a,d_a,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);
    Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1);//Cuando es energia es true/1 calcula energia potencial, cuando es false/0 calcula la cinetica
    Errorcuda(cudaMemcpy(h_eit,d_eit,sizeof(double),cudaMemcpyDeviceToHost),"eit",2);
    
    eit=h_eit[0]/2.0;
    ect=CalculoEnergiaCinetica(nd,np,v);
    temp= ect/nd;
    ect=ect/2.0;
    ett = ect + eit;

    std::cout << "ett,ect,eit,temp " << ett << " " << ect << " " <<eit << " " << temp << std::endl;
    
    ti=clock();
    std::cout << " " << std::endl;
    std::cout << "Resultados parciales" << std::endl;
    std::cout << "ic,temp,dens,ett,ect,eit,dtt " << std::endl;
    ofasres << "Resultados parciales" << std::endl;
    ofasres << "ic,temp,dens,ett,ect,eit,dtt " << std::endl;
    eis=ns=0;
    ofasat << "ip,p[],v[],a[]"<< std::endl;
    

    for(ic=0;ic<nc;ic++){
        
        nconf=false;
        if(ic%ncc==0)nconf=true;

        for(ip=0; ip<np; ip++){
            for(id=0; id<nd; id++)p[id+ip*nd] +=(v[id+ip*nd] * dt + a[id+nd*ip]*dt*dt*0.5);
            if(p[ip*nd] > caja.x){p[ip*nd] -= caja.x;}
            if(p[ip*nd] < 0){p[ip*nd] += caja.x;}
            if(p[1+ip*nd] > caja.y){p[1+ip*nd] -= caja.y;}
            if(p[1+ip*nd] < 0){p[1+ip*nd] += caja.y;}
            if(p[2+ip*nd] > caja.z){p[2+ip*nd] -= caja.z;}
            if(p[2+ip*nd] < 0){p[2+ip*nd] += caja.z;}
        }

        Velocidades(np,nd,v,a,dt);

        Errorcuda(cudaMemcpy(d_p,p,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"p",1);
        AceleracionesfFuerzasLJ<<<blockspergrid,threadsperblock>>>(np,d_p,chp,param,caja,cajai,condper,d_eit,d_a,nconf);
        Errorcuda(cudaGetLastError(),"PLJ",3);
        Errorcuda(cudaMemcpy(a,d_a,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);    
        
        Velocidades(np,nd,v,a,dt);
        
        if(ic>0 && ic % ncc == 0){
            tf = clock();
            dtt =((double)(tf - ti))/CLOCKS_PER_SEC;
            Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1);
            Errorcuda(cudaMemcpy(h_eit,d_eit,sizeof(double),cudaMemcpyDeviceToHost),"eit",2);
            eit=h_eit[0]/2.0;
            ect = CalculoEnergiaCinetica(nd,np,v);
            temp = ect/nd;
            ect=ect/2.0;
            ett = ect + eit;

            ets = ets + ett;
            ecs = ecs + ect;
            eis = eis + eit;
            temps = temps + temp;
            ns++;

            std::cout << ic<< " "<< temp<< " "<< dens<< " "<< ett<<
            " "<< ect << " "<< eit << " " << dtt <<std::endl;

            ofasres << ic<< " "<< temp<< " "<< dens<< " "<< ett<<
            " "<< ect << " "<< eit << " " << dtt <<std::endl;

            dens=dens;
            IAaD(np,ofasat,p,v,a,nd);
        }
    }
    Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1);
    Errorcuda(cudaMemcpy(h_eit,d_eit,sizeof(double),cudaMemcpyDeviceToHost),"eit",2);
    ets += ett;
    ecs += ect;
    eis += h_eit[0]/2.0;
    temps += temp;
    ns++;

    etp =ets / ns;
    eip =eis / ns;
    ecp =ecs / ns;

    std::cout << " " << std::endl;
    std::cout << "Resultados finales" << std::endl;
    tf = clock();
    dtt = ((double)(tf - ti))/CLOCKS_PER_SEC;
    ect=CalculoEnergiaCinetica(np,nd,v);

    std::cout << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    std::cout << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;

    ofasres << std::endl;
    ofasres <<"Resultados finales"<<std::endl;
    ofasres << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    ofasres << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;
    

    IAaD(np,ofasat,p,v,a,nd);
    

    Errorcuda(cudaFree(d_eit),"d_eit",4);
    Errorcuda(cudaFree(d_p),"d_p",4);
    //Errorcuda(cudaFree(d_v),"d_v",4);
    Errorcuda(cudaFree(d_a),"d_a",4);
    Errorcuda(cudaFree(param),"param",4);
    cudaDeviceReset();
}
#endif