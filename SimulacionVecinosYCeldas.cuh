#ifndef SIMVECYCEL_HEADER
#define SIMVECYCEL_HEADER

#include <stdio.h>

#include "OperacionesDeHilosyBloques.hpp"
#include "Potenciales.h"
#include "Optimizaciones.cuh"

__global__ void AceleracionesfFuerzasLJVYC(uint np,const double *p,int chp,double *param,double3 caja,double3 cajai,int3 condper,double *epot,double *a,int *vec,unsigned int *nvec,int nmaxvec,double rc,bool nconf)
{
    int gid=threadIdx.x+blockDim.x*blockIdx.x;
    gid/=chp;
    int lane=threadIdx.x&(chp-1);
    int nvpt=nvec[gid];
    int j;
    extern __shared__ double s_mem[];

    s_mem[0]=param[0];
    s_mem[1]=param[1];
    __syncthreads();
    
    double dis,pot=0.0;
    double2 fuepot;
    double3 fuerza,dif;
    fuerza=InitDataType3<double3>(0.0,0.0,0.0);
    if(gid>=np)return;
    double r2c=rc*rc;
    double pix,piy,piz,pjx,pjy,pjz;
    double3 pi,pj;
    pix=__ldg(p+3*gid);
    piy=__ldg(p+3*gid+1);
    piz=__ldg(p+3*gid+2);
    pi=InitDataType3<double3>(pix,piy,piz);
    #pragma unroll
    for(int jp=lane;jp<nvpt;jp+=chp){
        j=vec[gid*nmaxvec+jp]; 
            pjx=__ldg(p+3*j);
            pjy=__ldg(p+3*j+1);
            pjz=__ldg(p+3*j+2);
            pj=InitDataType3<double3>(pjx,pjy,pjz);
            dif.x=pi.x-pj.x;
            dif.y=pi.y-pj.y;
            dif.z=pi.z-pj.z;
            dif=CondPeriodicas(condper,caja,dif,cajai);
            dis=Discuad(dif);
            if(dis<=r2c){
                fuepot=InteraccionLJ(gid,j,dis,s_mem,nconf,true,r2c);
                fuerza.x+=fuepot.x*dif.x;
                fuerza.y+=fuepot.x*dif.y;
                fuerza.z+=fuepot.x*dif.z;
                if(nconf){
                    pot+=fuepot.y;
                } 
            }    
    }
   
    __syncwarp();
    for(int i=chp/2.0;i>=1;i/=2.0){
        fuerza.x+=__shfl_down_sync(FULL_MASK,fuerza.x,i,chp);
        fuerza.y+=__shfl_down_sync(FULL_MASK,fuerza.y,i,chp);
        fuerza.z+=__shfl_down_sync(FULL_MASK,fuerza.z,i,chp);
        if(nconf)pot+=__shfl_down_sync(FULL_MASK,pot,i);
    }
    __syncwarp();
    
    
    if(!lane){
        if(nconf)epot[gid]=pot;
        a[3*gid]=fuerza.x;
        a[3*gid+1]=fuerza.y;
        a[3*gid+2]=fuerza.z;
    }
    
}

void SimulacionVYC(uint np,int nd,double *p,double *v,double *a,double sig,double eps,double3 caja,
                 double3 cajai,int3 condper,double temp,std::ofstream &ofasres,std::ofstream &ofasat,
                 uint nc,double dt,double dens,uint ncp,double rc,double rbuf,int nhilos,int pot,int maxhilos,size_t memoria_global)
{
    double rcc=rc+rbuf;
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
    //Cuantos hilos calculan la fuerza aplicada a una particula, tiene que cumplir warpsize%chp=0
    int chp=HilosPorParticula(nhilos);
    int chpc=1;

    //calculamos cuantos hilos requerimos para cada uno de los kernels
    int threadsperblock,blockspergrid,threadsperblockvec,blockspergridvec,blockspergridcel,threadsperblockcel;
    
    printf("\nPara la rutina de aceleraciones\n");
    HilosenBloqueMultiplodeWarp(np,chp,blockspergrid,threadsperblock,MINBLOCKPERGRID,maxhilos);

    printf("\nPara la rutina de celdas\n");
    HilosenBloqueMultiplodeWarp(np,chpc,blockspergridcel,threadsperblockcel,MINBLOCKPERGRID,maxhilos);

    printf("\nPara la rutina de vecinos\n");
    HilosenBloqueMultiplodeWarp(np,chp,blockspergridvec,threadsperblockvec,MINBLOCKPERGRID,maxhilos);
    /***************************************************************************************************************************/
    //para el potencial de Lennard-Jones tenemos 2 parametros pl6 y plj12
    //para otros potenciales podemos llegar a tener mas parametros
    int nparam;
    switch(pot){
        case 1:
        nparam=2;
        break;
    }
    printf("nparam: %d\n",nparam);
    
    int nmaxvec=CuantosVecCaben(rcc,sig,dens,np);
    double3 tamcel,invtamcel;
    int nceldas=CrearCeldas(caja,tamcel,invtamcel);
    int nmax_particulas_en_celda=CuantasPartEnCel(sig,tamcel);
    int dc=rc+0.5;
    int n_de_celdas_vecinas=pow(1+2*dc,3);

    /***************************************************************************************************************************/
    //Memoria para arreglos en la CPU
    double *h_param=new double [nparam];
    double *h_eit=new double[1];
    /***************************************************************************************************************************/
    //Memoria para arreglos de GPU
    cudaStream_t stream[2];
    for(int i=0;i<2;i++)cudaStreamCreate(&stream[i]);
    double *param,*d_eit,*d_p,*d_a;
    int *vec,*p_cel,*celdas_vecinas;
    unsigned int *nvec,*np_cel;
    int *h_celdas_vecinas=new int[n_de_celdas_vecinas*nceldas];
    
    Errorcuda(cudaMalloc(&param,nparam*sizeof(double)),"param",0);
    Errorcuda(cudaMalloc(&d_eit,sizeof(double)*np),"d_eit",0);
    Errorcuda(cudaMalloc(&d_p,sizeof(double)*nd*np),"d_p",0);
    Errorcuda(cudaMalloc(&d_a,sizeof(double)*nd*np),"d_a",0);
    Errorcuda(cudaMalloc(&nvec,sizeof(int)*np),"nvec",0);
    Errorcuda(cudaMalloc(&vec,sizeof(int)*nmaxvec*np),"vec",0);
    Errorcuda(cudaMalloc(&p_cel,sizeof(int)*nmax_particulas_en_celda*nceldas),"p_cel",0);
    Errorcuda(cudaMalloc(&np_cel,sizeof(int)*nceldas),"np_cel",0);
    Errorcuda(cudaMalloc(&celdas_vecinas,sizeof(int)*n_de_celdas_vecinas*nceldas),"celdas_vecinas",0);

    size_t memoria_global_utilizada=(nparam+np*(1+2*nd))*sizeof(double)+(np*(1+nmaxvec)+nceldas*(1+nmax_particulas_en_celda+n_de_celdas_vecinas))*sizeof(int);    
    OcupacionDeMemoriaGlobal(memoria_global_utilizada,memoria_global);
    /***************************************************************************************************************************/

    CalculoDeCeldasVecinas(caja,h_celdas_vecinas,dc,n_de_celdas_vecinas);
    Errorcuda(cudaMemcpy(celdas_vecinas,h_celdas_vecinas,sizeof(int)*n_de_celdas_vecinas*nceldas,cudaMemcpyHostToDevice),"celdas_vecinas",1);
    CalculoParamLJ(sig,eps,h_param);
    
    Errorcuda(cudaMemcpy(param,h_param,sizeof(double)*nparam,cudaMemcpyHostToDevice),"param",1);
    Errorcuda(cudaMemset(np_cel,0,sizeof(int)*nceldas),"memset",0);
    Errorcuda(cudaMemset(nvec,0,sizeof(int)*np),"memset",0);
    Errorcuda(cudaMemcpy(d_p,p,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"p",1);

    CalculoCeldas<<<blockspergridcel,threadsperblockcel>>>(np,nmax_particulas_en_celda,np_cel,p_cel,d_p,invtamcel,caja,nceldas);
    Errorcuda(cudaGetLastError(),"Calculo de celdas",3);
    cudaDeviceSynchronize();
    CalculoDeVecinosConCeldas<<<blockspergridvec,threadsperblockvec>>>(np,chp,nmaxvec,nmax_particulas_en_celda,vec,p_cel,np_cel,nvec,condper,caja,cajai,d_p,rc,n_de_celdas_vecinas,celdas_vecinas,invtamcel);
    Errorcuda(cudaGetLastError(),"Calculo de vecinos",3);
    cudaDeviceSynchronize();

    AceleracionesfFuerzasLJVYC<<<blockspergrid,threadsperblock>>>(np,d_p,chp,param,caja,cajai,condper,d_eit,d_a,vec,nvec,nmaxvec,rc,nconf);
    Errorcuda(cudaGetLastError(),"PLJ",3);
    Errorcuda(cudaMemcpy(a,d_a,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);

    Reduccionconwarps<1024><<<1,1024,0,stream[1]>>>(d_eit,np,1);//Cuando es energia es true/1 calcula energia potencial, cuando es false/0 calcula la cinetica
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

            //Ciclo de dimensiones
            for(id=0; id<nd; id++){
                p[id+ip*nd] +=(v[id+ip*nd] * dt + a[id+nd*ip]*dt*dt*0.5);              
            }
            if(p[ip*nd] > caja.x){p[ip*nd] -= caja.x;}
            if(p[ip*nd] < 0){p[ip*nd] += caja.x;}
            if(p[1+ip*nd] > caja.y){p[1+ip*nd] -= caja.y;}
            if(p[1+ip*nd] < 0){p[1+ip*nd] += caja.y;}
            if(p[2+ip*nd] > caja.z){p[2+ip*nd] -= caja.z;}
            if(p[2+ip*nd] < 0){p[2+ip*nd] += caja.z;}
        }

        Velocidades(np,nd,v,a,dt);
       
        Errorcuda(cudaMemcpy(d_p,p,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"p",1);
        if(nc%10==0){
            cudaMemset(np_cel,0,sizeof(int)*nceldas);
            cudaMemset(nvec,0,sizeof(int)*np);
            cudaDeviceSynchronize();
            CalculoCeldas<<<blockspergridcel,threadsperblockcel>>>(np,nmax_particulas_en_celda,np_cel,p_cel,d_p,invtamcel,caja,nceldas);
            Errorcuda(cudaGetLastError(),"Calculo de celdas",3);
            cudaDeviceSynchronize();
            CalculoDeVecinosConCeldas<<<blockspergridvec,threadsperblockvec>>>(np,chp,nmaxvec,nmax_particulas_en_celda,vec,p_cel,np_cel,nvec,condper,caja,cajai,d_p,rc,n_de_celdas_vecinas,celdas_vecinas,invtamcel);
            Errorcuda(cudaGetLastError(),"Calculo de vecinos",3);
            cudaDeviceSynchronize();
        }
        AceleracionesfFuerzasLJVYC<<<blockspergrid,threadsperblock>>>(np,d_p,chp,param,caja,cajai,condper,d_eit,d_a,vec,nvec,nmaxvec,rc,nconf);
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
    Reduccionconwarps<1024><<<1,1024,0,stream[1]>>>(d_eit,np,1);
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
    
      
    for(int i=0;i<2;i++)cudaStreamDestroy(stream[i]);
    Errorcuda(cudaFree(d_eit),"d_eit",4);
    Errorcuda(cudaFree(d_p),"d_eit",4);
    Errorcuda(cudaFree(d_a),"d_eit",4);
    Errorcuda(cudaFree(nvec),"d_eit",4);
    Errorcuda(cudaFree(vec),"d_eit",4);
    Errorcuda(cudaFree(p_cel),"p_cel",4);
    Errorcuda(cudaFree(np_cel),"np_cel",4);
    Errorcuda(cudaFree(celdas_vecinas),"celdas_vecinas",4);
    Errorcuda(cudaFree(param),"param",4);
    delete[] h_celdas_vecinas;
    delete[] h_eit;
    cudaDeviceReset();
}
#endif