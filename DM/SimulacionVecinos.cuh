#ifndef SIMVEC_HEADER
#define SIMVEC_HEADER

#include <stdio.h>

#include "../MISC/OperacionesDeHilosyBloques.cuh"
#include "Potenciales.cuh"
#include "Optimizaciones.cuh"
#include "PotencialesRestriccion.cuh"
#include "Termostatos.hpp"
#include "rattle.cuh"
#include "../MISC/PropGPU.cuh"

__global__ void AceleracionesfFuerzasLJV(uint np,uint n_esp_p,uint pot_int,uint *esp_de_p,const double *p,int chp,double *M_param,double3 caja,double3 cajai,int3 condper,double *epot,double *a,int *vec,unsigned int *nvec,int nmaxvec,double rc,bool nconf)
{
    int gid=threadIdx.x+blockDim.x*blockIdx.x;
    gid/=chp;
    int lane=threadIdx.x&(chp-1);
    int nvpt=nvec[gid];
    int j;
    uint esp1=0,esp2=0,elem_M = 0;
    
    double dis,pot=0.0;
    double2 fuepot;
    double3 fuerza,dif;
    fuerza=InitDataType3<double3,double>(0.0,0.0,0.0);
    if(gid>=np)return;
    double r2c=rc*rc;
    double pix,piy,piz,pjx,pjy,pjz;
    double3 pi,pj;
    pix=__ldg(p+3*gid);
    piy=__ldg(p+3*gid+1);
    piz=__ldg(p+3*gid+2);
    pi=InitDataType3<double3,double>(pix,piy,piz);
    esp1=esp_de_p[gid];
    #pragma unroll
    for(int jp=lane;jp<nvpt;jp+=chp){
        j=vec[gid*nmaxvec+jp];
        esp2=esp_de_p[j];
        elem_M = n_esp_p*esp1+esp2; 
        pjx=__ldg(p+3*j);
        pjy=__ldg(p+3*j+1);
        pjz=__ldg(p+3*j+2);
        pj=InitDataType3<double3,double>(pjx,pjy,pjz);
        dif.x=pi.x-pj.x;
        dif.y=pi.y-pj.y;
        dif.z=pi.z-pj.z;
        dif=CondPeriodicas(condper,caja,dif,cajai);
        dis=Discuad(dif);
        if(dis<=r2c){
            fuepot=Interaccion(pot_int,n_esp_p,elem_M,gid,j,dis,M_param,nconf,true,r2c);
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

void SimulacionV(uint nc,uint ncp,uint np,uint n_esp_p,uint n_esp_m,uint nparam,uint pot,uint max_p_en_esp_mr,
                 uint ensamble,uint termos,uint *esp_de_p,uint *M_int,uint *p_en_m,uint *n_m_esp_mr,uint *n_p_esp_m,
                 int nhilos,int maxhilos,float *constr,bool vibrante,bool p_o_m,uint3 *mad_de_p,int3 condper,double rc,double rbuf,double dens,
                 double dt,double temp_d,double param_termo,double *param,double *pos,double *vel,
                 double *acel,double *q_rat,double *dis_p_esp_mr_rep,double3 caja,double3 cajai,
                 size_t memoria_global,std::ofstream &ofasres,std::ofstream &ofasat)
{
    double rcc=rc+rbuf;
    int ncc=nc*ncp/100.0;
    clock_t ti, tf;
    double dtt;

    double ets=0.0,etp=0.0,ett=0.0;
    double ecs=0.0,ecp=0.0,ect=0.0;
    double eis=0.0,eip=0.0,eit=0.0;
    double temps,temp;

    int ns,id,ip,ic=0;

    bool nconf=true;
    /***************************************************************************************************************************/
    //Cuantos hilos calculan la fuerza aplicada a una particula, tiene que cumplir warpsize%chp=0
    int chp=HilosPorParticula(nhilos);

    //calculamos cuantos hilos requerimos para cada uno de los kernels
    int threadsperblock,blockspergrid,threadsperblockvec,blockspergridvec;
    
    printf("\nPara la rutina de aceleraciones\n");
    HilosenBloqueMultiplodeWarp(np,chp,blockspergrid,threadsperblock,MINBLOCKPERGRID,maxhilos);

    printf("\nPara la rutina de vecinos\n");
    HilosenBloqueMultiplodeWarp(np,chp,blockspergridvec,threadsperblockvec,MINBLOCKPERGRID,maxhilos);
    /***************************************************************************************************************************/
    //para el potencial de Lennard-Jones tenemos 2 parametros pl6 y plj12
    //para otros potenciales podemos llegar a tener mas parametros
    /***************************************************************************************************************************/
    //Memoria para arreglos en la CPU
    
    double *h_M_param=new double [nparam*n_esp_p*n_esp_p];
    double *h_eit=new double[1];
    double s_nh_a = 0.,s_nh_v=0.,s_nh_p=1.,g1=0.;

    InicializarMatrizDeParametros(pot,n_esp_p,nparam,param,h_M_param);
    int nmaxvec=CuantosVecCaben(rcc,param,dens,np,n_esp_p,nparam);
    /***************************************************************************************************************************/
    //Memoria para arreglos de GPU
    double *d_eit,*d_pos,*d_acel,*d_M_param;
    uint3 *d_mad_de_p;
    int *d_vec;
    uint *d_nvec,*d_esp_de_p,*d_M_int;

    Errorcuda(cudaMalloc(&d_vec,sizeof(int)*nmaxvec*np),"vec",0);

    Errorcuda(cudaMalloc(&d_mad_de_p,np*sizeof(uint3)),"mad_de_p",0);

    Errorcuda(cudaMalloc(&d_esp_de_p,np*sizeof(uint)),"esp_de_p",0);
    Errorcuda(cudaMalloc(&d_M_int,n_esp_p*n_esp_p*sizeof(uint)),"d_M_int",0);
    Errorcuda(cudaMalloc(&d_nvec,sizeof(uint)*np),"nvec",0);

    Errorcuda(cudaMalloc(&d_M_param,nparam*n_esp_p*n_esp_p*sizeof(double)),"param",0);
    Errorcuda(cudaMalloc(&d_eit,sizeof(double)*np),"d_eit",0);
    Errorcuda(cudaMalloc(&d_pos,sizeof(double)*nd*np),"d_p",0);
    Errorcuda(cudaMalloc(&d_acel,sizeof(double)*nd*np),"d_a",0);
    
    
    Errorcuda(cudaMemcpy(d_mad_de_p,mad_de_p,np*sizeof(uint3),cudaMemcpyHostToDevice),"mad_de_p",1);
    Errorcuda(cudaMemcpy(d_esp_de_p,esp_de_p,np*sizeof(uint),cudaMemcpyHostToDevice),"esp_de_p",1);
    Errorcuda(cudaMemcpy(d_M_int,M_int,n_esp_p*n_esp_p*sizeof(uint),cudaMemcpyHostToDevice),"M_int",1);

    Errorcuda(cudaMemcpy(d_M_param,h_M_param,sizeof(double)*nparam*n_esp_p*n_esp_p,cudaMemcpyHostToDevice),"M_param",1);
    Errorcuda(cudaMemcpy(d_pos,pos,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"pos",1);    
    

    size_t memoria_global_utilizada=(nparam+np*(1+2*nd))*sizeof(double)+(np*(1+nmaxvec))*sizeof(int);
    
    OcupacionDeMemoriaGlobal(memoria_global_utilizada,memoria_global);
    /***************************************************************************************************************************/

    cudaMemset(d_nvec,0,sizeof(int)*np);
    CalculoDeVecinos<<<blockspergridvec,threadsperblockvec>>>(np,n_esp_p,d_M_int,d_esp_de_p,chp,d_pos,d_vec,d_nvec,nmaxvec,condper,d_mad_de_p,caja,cajai,rcc);
    Errorcuda(cudaGetLastError(),"Calculo de vecinos",3);
    cudaDeviceSynchronize();

    AceleracionesfFuerzasLJV<<<blockspergrid,threadsperblock>>>(np,n_esp_p,pot,d_esp_de_p,d_pos,chp,d_M_param,caja,cajai,condper,d_eit,d_acel,d_vec,d_nvec,nmaxvec,rc,nconf);
    Errorcuda(cudaGetLastError(),"PLJ",3); 
    
    //Aqui podria considerar ponerlo en un stream para que se haga est la copia de la aceleraciones y el calculo de energia al mismo tiempo, ya que no dependen los porcesos uno de los otros
    Errorcuda(cudaMemcpy(acel,d_acel,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);
    Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1);//Cuando es energia es true/1 calcula energia potencial, cuando es false/0 calcula la cinetica
    Errorcuda(cudaMemcpy(h_eit,d_eit,sizeof(double),cudaMemcpyDeviceToHost),"eit",2);
    eit=0.;
    if(p_o_m&&vibrante)
        eit+=PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,condper,mad_de_p,constr[0],pos,acel,dis_p_esp_mr_rep,caja,cajai);
    
    eit+=h_eit[0]/2.0;
    ect=CalculoEnergiaCinetica(np,vel);
    temp= ect/nd;
    ect=ect/2.0;
    ett = ect + eit;

    ti=clock();
    std::cout << " " << std::endl;
    std::cout << "Resultados parciales" << std::endl;
    std::cout << "ic,temp,dens,ett,ect,eit,dtt " << std::endl;
    std::cout << ic<< " "<< temp<< " "<< dens<< " "<< ett<<" "<< ect << " "<< eit << " 0.0" <<std::endl;
    ofasres << "Resultados parciales" << std::endl;
    ofasres << "ic,temp,dens,ett,ect,eit,dtt " << std::endl;
    ofasres << ic<< " "<< temp<< " "<< dens<< " "<< ett<<" "<< ect << " "<< eit << " 0.0" <<std::endl;
    eis=ns=0;
    ofasat << "ip,p[],v[],a[]"<< std::endl;
    
    
    for(ic=0;ic<nc;ic++){
        nconf=false;
        if(ic%ncc==0)nconf=true;

        for(ip=0; ip<np; ip++)
            for(id=0; id<nd; id++)
                q_rat[id+ip*nd] =(vel[id+ip*nd] + acel[id+nd*ip]*dt*0.5);

        //aquí se hace el algoritmo de RATTLE
        if(!vibrante&&p_o_m)
            RattlePos(constr[1],n_esp_m,np,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,mad_de_p,condper,constr[0],dt,pos,q_rat,dis_p_esp_mr_rep,caja,cajai);

        for(ip=0; ip<np; ip++){   
            
            for(id=0; id<nd; id++)
                pos[id+ip*nd] += dt*q_rat[id+ip*nd];
            
            if(pos[ip*nd] > caja.x){pos[ip*nd] -= caja.x;}
            if(pos[ip*nd] < 0){pos[ip*nd] += caja.x;}
            if(pos[1+ip*nd] > caja.y){pos[1+ip*nd] -= caja.y;}
            if(pos[1+ip*nd] < 0){pos[1+ip*nd] += caja.y;}
            if(pos[2+ip*nd] > caja.z){pos[2+ip*nd] -= caja.z;}
            if(pos[2+ip*nd] < 0){pos[2+ip*nd] += caja.z;}
        }
        if(!vibrante&&(max_p_en_esp_mr-1)){
            for(ip=0; ip<np; ip++)
                for(id=0; id<nd; id++)
                     vel[id+ip*nd] = q_rat[id+ip*nd];

        }else{
            Velocidades(np,vel,acel,dt);
        }
       
        Errorcuda(cudaMemcpy(d_pos,pos,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"p",1);
        if(nc%10==0){
            cudaMemset(d_nvec,0,sizeof(int)*np);
            CalculoDeVecinos<<<blockspergridvec,threadsperblockvec>>>(np,n_esp_p,d_M_int,d_esp_de_p,chp,d_pos,d_vec,d_nvec,nmaxvec,condper,d_mad_de_p,caja,cajai,rcc);
            Errorcuda(cudaGetLastError(),"Calculo de vecinos",3);
            cudaDeviceSynchronize();
        }

        AceleracionesfFuerzasLJV<<<blockspergrid,threadsperblock>>>(np,n_esp_p,pot,d_esp_de_p,d_pos,chp,d_M_param,caja,cajai,condper,d_eit,d_acel,d_vec,d_nvec,nmaxvec,rc,nconf);
        Errorcuda(cudaGetLastError(),"PLJ",3);
        Errorcuda(cudaMemcpy(acel,d_acel,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);eit=0.;

        if(p_o_m&&vibrante)
            eit+=PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,condper,mad_de_p,constr[0],pos,acel,dis_p_esp_mr_rep,caja,cajai);    
        Velocidades(np,vel,acel,dt);
        //debido a que se necesitan F(t+dt) entonces RATTLEvel se realiza aquí
        if(!vibrante&&p_o_m)
            RattleVel(constr[1],n_esp_m,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,mad_de_p,constr[0],pos,vel,dis_p_esp_mr_rep);
        

        if(ensamble==1){
            ect = CalculoEnergiaCinetica(np,vel);
            g1=(ect-3.*temp_d)*np;
            temp = ect/nd;
            ect=ect/2.0;
            if(termos==4){
                for(int i=0;i<np;i++){
                    acel[nd*i]-=s_nh_v*vel[nd*i]/s_nh_p;
                    acel[nd*i+1]-=s_nh_v*vel[nd*i+1]/s_nh_p;
                    acel[nd*i+2]-=s_nh_v*vel[nd*i+2]/s_nh_p;
                }
                s_nh_p+=s_nh_v*dt+s_nh_a*dt*dt*0.5;
                s_nh_v+=s_nh_a*0.5*dt;
                s_nh_a = s_nh_v*s_nh_v/s_nh_p+g1*s_nh_p/param_termo;
                s_nh_v+=s_nh_a*0.5*dt;
            }else{
                Termostato(termos,np,temp,temp_d,dt,param_termo,vel);
            }
        }

        if(ic>0 && ic % ncc == 0){
            tf = clock();
            dtt =((double)(tf - ti))/CLOCKS_PER_SEC;
            Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1);
            Errorcuda(cudaMemcpy(h_eit,d_eit,sizeof(double),cudaMemcpyDeviceToHost),"eit",2);
            eit+=h_eit[0]/2.0;
            ett = ect + eit;
            ect = CalculoEnergiaCinetica(np,vel);
            temp = ect/nd;
            ect=ect/2.0;

            ets += ett;
            ecs += ect;
            eis += eit;
            temps += temp;
            ns++;

            std::cout << ic<< " "<< temp<< " "<< dens<< " "<< ett<<
            " "<< ect << " "<< eit << " " << dtt <<std::endl;

            ofasres << ic<< " "<< temp<< " "<< dens<< " "<< ett<<
            " "<< ect << " "<< eit << " " << dtt <<std::endl;

            dens=dens;
            IAaD(np,ofasat,pos,vel,acel);
        }
    }
    
    Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1);
    Errorcuda(cudaMemcpy(h_eit,d_eit,sizeof(double),cudaMemcpyDeviceToHost),"eit",2);
    eit=h_eit[0]/2.0;
    ett = ect + eit;
    ect = CalculoEnergiaCinetica(np,vel);
    temp = ect/nd;
    ect=ect/2.0;

    std::cout << ic<< " "<< temp<< " "<< dens<< " "<< ett<<
    " "<< ect << " "<< eit << " " << dtt <<std::endl;

    ofasres << ic<< " "<< temp<< " "<< dens<< " "<< ett<<
    " "<< ect << " "<< eit << " " << dtt <<std::endl;

    ets += ett;
    ecs += ect;
    eis += eit;
    temps += temp;
    ns++;
    
    etp =ets / ns;
    eip =eis / ns;
    ecp =ecs / ns;

    std::cout << " " << std::endl;
    std::cout << "Resultados finales" << std::endl;
    tf = clock();
    dtt = ((double)(tf - ti))/CLOCKS_PER_SEC;
    ect=CalculoEnergiaCinetica(np,vel);

    std::cout << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    std::cout << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;

    ofasres << std::endl;
    ofasres <<"Resultados finales"<<std::endl;
    ofasres << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    ofasres << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;
    

    IAaD(np,ofasat,pos,vel,acel);
    
      
    Errorcuda(cudaFree(d_vec),"vec",4);

    Errorcuda(cudaFree(d_mad_de_p),"mad_de_p",4);

    Errorcuda(cudaFree(d_esp_de_p),"esp_de_p",4);
    Errorcuda(cudaFree(d_M_int),"d_M_int",4);
    Errorcuda(cudaFree(d_nvec),"nvec",4);

    Errorcuda(cudaFree(d_M_param),"param",4);
    Errorcuda(cudaFree(d_eit),"d_eit",4);
    Errorcuda(cudaFree(d_pos),"d_p",4);
    Errorcuda(cudaFree(d_acel),"d_a",4);
    cudaDeviceReset();
}
#endif
