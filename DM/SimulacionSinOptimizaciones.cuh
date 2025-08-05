#ifndef SIMNO_HEADER
#define SIMNO_HEADER

#include <stdio.h>

#include "../MISC/OperacionesTDatosCuda.cuh"
#include "../MISC/OperacionesDeHilosyBloques.cuh"
#include "Potenciales.cuh"
#include "Termostatos.hpp"
#include "PotencialesRestriccion.cuh"
#include "rattle.cuh"
#include "../MISC/PropGPU.cuh"

//Nota, por ahora estoy considerando 2 loops para evitar hacer el calculo de fuerzas y potencial entre particulas que estan en misma molecula, esto evita tambien auto interaccion
//si la quito parece ser que si hay autointeraccion y eso no me gusta, lo arreglare pronto, segun yo por eso agregue la condicion en la rutina interaccionlj pero no funciona

__global__ void AceleracionesFuerzas(uint n_esp_p,uint pot,uint np,uint3 *m_de_p,uint *esp_de_p,uint *M_int,const double *pos,int chp,double *M_param,double3 caja,double3 cajai,int3 condper,double *epot,double *a,bool nconf)
{
    int gid=threadIdx.x+blockDim.x*blockIdx.x;
    gid/=chp;
    int lane=threadIdx.x &(chp-1);
    uint esp1=0,esp2=0,elem_M = 0;
    double dis,spot=0.0;
    double2 fuepot;
    double3 fuerza,dif;
    fuerza=InitDataType3<double3,double>(0.0,0.0,0.0);    
    if(gid>=np)return;
    double pix,piy,piz,pjx,pjy,pjz;
    double3 pi,pj;
    pix=__ldg(pos+3*gid);
    piy=__ldg(pos+3*gid+1);
    piz=__ldg(pos+3*gid+2);
    pi=InitDataType3<double3,double>(pix,piy,piz);
    esp1=esp_de_p[gid];

    #pragma unroll
    for(int j=lane;j<np;j+=chp){
        esp2=esp_de_p[j];
        elem_M = n_esp_p*esp1+esp2;
        if(M_int[n_esp_p*esp1+esp2]){
            if(m_de_p[gid].x-m_de_p[j].x){
                pjx=__ldg(pos+3*j);
                pjy=__ldg(pos+3*j+1);
                pjz=__ldg(pos+3*j+2);
                pj=InitDataType3<double3,double>(pjx,pjy,pjz);
                dif.x=pi.x-pj.x;
                dif.y=pi.y-pj.y;
                dif.z=pi.z-pj.z;
                dif=CondPeriodicas(condper,caja,dif,cajai);
                dis=Discuad(dif);
                fuepot=Interaccion(pot,n_esp_p,elem_M,gid,j,dis,M_param,nconf,false);
                
                fuerza.x+=fuepot.x*dif.x;
                fuerza.y+=fuepot.x*dif.y;
                fuerza.z+=fuepot.x*dif.z;
                
                if(nconf){
                    spot+=fuepot.y;
                }
            }
        }   
    }
    
    
    //reduccion de warps aplicado a fuerza y al potencial
    __syncwarp();
    for(int i=chp/2.0;i>=1;i/=2.0){
        fuerza.x+=__shfl_down_sync(FULL_MASK,fuerza.x,i,chp);
        fuerza.y+=__shfl_down_sync(FULL_MASK,fuerza.y,i,chp);
        fuerza.z+=__shfl_down_sync(FULL_MASK,fuerza.z,i,chp);
        if(nconf)spot+=__shfl_down_sync(FULL_MASK,spot,i,chp);
    }
    __syncwarp();
    
    
    if(!lane){
        if(nconf)epot[gid]=spot;
        a[3*gid]=fuerza.x;
        a[3*gid+1]=fuerza.y;
        a[3*gid+2]=fuerza.z;
    }
    
}

void Simulacion(uint nc,uint ncp,uint nhilos,uint np,uint nparam,uint pot,uint n_esp_p,uint n_esp_m,uint max_p_en_esp_mr,uint ensamble,uint termos,uint max_it,
                uint3 *mad_de_p,uint *esp_de_p,uint *n_m_esp_mr,uint *n_p_esp_m,uint *M_int,uint *p_en_m,int maxhilos,bool vibrante,
                size_t memoria_global,double dt,double dens,double kres,double temp_d,double param_termo,double tol,double *param,double *pos,double *q_rat,double *vel,
                double *acel,double *dis_p_esp_mr_rep,int3 condper,double3 caja, double3 cajai,std::ofstream &ofasres,std::ofstream &ofasat)
{
    int ncc=nc*ncp/100.0;
    clock_t ti, tf;
    double dtt;
    double temp;

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
    
    /***************************************************************************************************************************/
    //Memoria para arreglos en la CPU
    
    double *h_M_param=new double [nparam*n_esp_p*n_esp_p];
    double *h_eit=new double[1];
    double s_nh_a = 0.,s_nh_v=0.,s_nh_p=1.,g1=0.;

    InicializarMatrizDeParametros(pot,n_esp_p,nparam,param,h_M_param);

    /***************************************************************************************************************************/
    //Memoria para arreglos en la GPU
    uint *d_esp_de_p,*d_M_int;
    uint3 *d_mad_de_p;
    double *d_eit,*d_pos,*d_acel,*d_M_param;

    Errorcuda(cudaMalloc(&d_mad_de_p,np*sizeof(uint3)),"mad_de_p",0);

    Errorcuda(cudaMalloc(&d_esp_de_p,np*sizeof(uint)),"esp_de_p",0);
    Errorcuda(cudaMalloc(&d_M_int,n_esp_p*n_esp_p*sizeof(uint)),"d_M_int",0);

    Errorcuda(cudaMalloc(&d_M_param,nparam*n_esp_p*n_esp_p*sizeof(double)),"param",0);
    Errorcuda(cudaMalloc(&d_eit,sizeof(double)*np),"d_eit",0);
    Errorcuda(cudaMalloc(&d_pos,sizeof(double)*nd*np),"d_p",0);
    Errorcuda(cudaMalloc(&d_acel,sizeof(double)*nd*np),"d_a",0);
    //vemos cuantos byte ocupan nuestros arreglos en la gpu y el porcentaje

    size_t memoria_global_utilizada =   (nparam*n_esp_p*n_esp_p+np*(1+2*nd))*sizeof(double)
                                    +   (n_esp_p*n_esp_p+np)*sizeof(uint)
                                    +   np*sizeof(uint3);
    OcupacionDeMemoriaGlobal(memoria_global_utilizada,memoria_global);

    Errorcuda(cudaMemcpy(d_mad_de_p,mad_de_p,np*sizeof(uint3),cudaMemcpyHostToDevice),"mad_de_p",1);
    Errorcuda(cudaMemcpy(d_esp_de_p,esp_de_p,np*sizeof(uint),cudaMemcpyHostToDevice),"esp_de_p",1);
    Errorcuda(cudaMemcpy(d_M_int,M_int,n_esp_p*n_esp_p*sizeof(uint),cudaMemcpyHostToDevice),"M_int",1);

    Errorcuda(cudaMemcpy(d_M_param,h_M_param,sizeof(double)*nparam*n_esp_p*n_esp_p,cudaMemcpyHostToDevice),"M_param",1);
    Errorcuda(cudaMemcpy(d_pos,pos,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"pos",1);

    /***************************************************************************************************************************/
    //Mas optimizaciones de cuda

    /***************************************************************************************************************************/
    
    
    
    

    AceleracionesFuerzas<<<blockspergrid,threadsperblock>>>(n_esp_p,pot,np,d_mad_de_p,d_esp_de_p,d_M_int,d_pos,chp,d_M_param,caja,cajai,condper,d_eit,d_acel,nconf);
    Errorcuda(cudaGetLastError(),"PLJ",3);

    //Aqui podria considerar ponerlo en un stream para que se haga est la copia de la aceleraciones y el calculo de energia al mismo tiempo, ya que no dependen los porcesos uno de los otros
    Errorcuda(cudaMemcpy(acel,d_acel,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);
    if(max_p_en_esp_mr-1)PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,condper,mad_de_p,kres,pos,acel,dis_p_esp_mr_rep,caja,cajai);
    Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1);//Cuando es energia es true/1 calcula energia potencial, cuando es false/0 calcula la cinetica
    Errorcuda(cudaMemcpy(h_eit,d_eit,sizeof(double),cudaMemcpyDeviceToHost),"eit",2);
    
    eit=h_eit[0]/2.0;
    ect=CalculoEnergiaCinetica(np,vel);
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
    
    // el calculo de posiciones y la primeras velocidades son independientes entre si, ambos dependen de la mitad de las velocidades anterior y la aceleracion
    // talvez en el futuro se podrian hacer de forma asincrona
    for(ic=0;ic<nc;ic++){
        
        nconf=false;
        if(ic%ncc==0)nconf=true;

        for(ip=0; ip<np; ip++)
            for(id=0; id<nd; id++)
                q_rat[id+ip*nd] =(vel[id+ip*nd] + acel[id+nd*ip]*dt*0.5);

        //aquí se hace el algoritmo de RATTLE
        if(!vibrante&&(max_p_en_esp_mr-1))RattlePos(max_it,n_esp_m,np,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,mad_de_p,condper,tol,dt,pos,q_rat,dis_p_esp_mr_rep,caja,cajai);

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
        AceleracionesFuerzas<<<blockspergrid,threadsperblock>>>(n_esp_p,pot,np,d_mad_de_p,d_esp_de_p,d_M_int,d_pos,chp,d_M_param,caja,cajai,condper,d_eit,d_acel,nconf);
        Errorcuda(cudaGetLastError(),"PLJ",3);
        Errorcuda(cudaMemcpy(acel,d_acel,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);
        if((max_p_en_esp_mr-1)&&vibrante)PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,condper,mad_de_p,kres,pos,acel,dis_p_esp_mr_rep,caja,cajai);    
        Velocidades(np,vel,acel,dt);
        //debido a que se necesitan F(t+dt) entonces RATTLEvel se realiza aquí
        if(!vibrante&&(max_p_en_esp_mr-1)){
            RattleVel(max_it,n_esp_m,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,mad_de_p,tol,pos,vel,dis_p_esp_mr_rep);
        }
        
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
            eit=h_eit[0]/2.0;
            ett = ect + eit;
            ect = CalculoEnergiaCinetica(np,vel);
            temp = ect/nd;
            ect=ect/2.0;

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
            IAaD(np,ofasat,pos,vel,acel);
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
    ect=CalculoEnergiaCinetica(np,vel);

    std::cout << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    std::cout << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;

    ofasres << std::endl;
    ofasres <<"Resultados finales"<<std::endl;
    ofasres << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    ofasres << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;
    

    IAaD(np,ofasat,pos,vel,acel);
    
    Errorcuda(cudaFree(d_mad_de_p),"mad_de_p",4);
    Errorcuda(cudaFree(d_esp_de_p),"esp_de_p",4);
    Errorcuda(cudaFree(d_M_int),"M_int",4);

    Errorcuda(cudaFree(d_M_param),"param",4);
    Errorcuda(cudaFree(d_eit),"d_eit",4);
    Errorcuda(cudaFree(d_pos),"d_p",4);
    Errorcuda(cudaFree(d_acel),"d_a",4);
    cudaDeviceReset();
}
#endif
