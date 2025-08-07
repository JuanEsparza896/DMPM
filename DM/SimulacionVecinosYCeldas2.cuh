#ifndef SIMVECYCEL_HEADER
#define SIMVECYCEL_HEADER

#include <stdio.h>

#include "../MISC/OperacionesDeHilosyBloques.cuh"
#include "Potenciales.cuh"
#include "Optimizaciones.cuh"
#include "Termostatos.hpp"
#include "PotencialesRestriccion.cuh"
#include "rattle.cuh"
#include "../MISC/PropGPU.cuh"

__global__ void CalculoDeVecinosConCeldas2(uint np,uint n_de_cel_vec,uint nmax_p_en_cel,uint n_esp_p,uint *esp_de_p,uint *cel_vec,uint *np_cel,uint *p_en_cel,
    uint *M_int,uint *nvec,int chp,int nmaxvec,int *vecinos,double rc,double *p,int3 condper,uint3 *mad_de_p,double3 caja,double3 cajai,double3 invtamcel)
{
int particula,nv,celda,celv,j;
uint eLx,eLy,fa,esp1,esp2,nv;
double pix,piy,piz,pjx,pjy,pjz,dis;
const int lane= threadIdx.x & (chp-1);
double3 pi,pj,dif;
int3 pos;

particula = threadIdx.x+blockIdx.x*blockDim.x;
particula/=chp;

if(particula>=np)return;

int cuantas_celdas_por_hilo=n_de_cel_vec/chp;
int ccll=0;
if(n_de_cel_vec%chp)cuantas_celdas_por_hilo++;

eLx=caja.x,eLy=caja.y;
fa=eLx*eLy;

pix=__ldg(p+3*particula);
piy=__ldg(p+3*particula+1);
piz=__ldg(p+3*particula+2);
esp1=esp_de_p[particula];
pi=InitDataType3<double3,double>(pix,piy,piz);

if(pix==caja.x&&piy==caja.y&&piz==caja.z)celda=0;
pos.x=pix*invtamcel.x;
pos.y=piy*invtamcel.y;
pos.z=piz*invtamcel.z;
celda=pos.x+eLx*pos.y+fa*pos.z;
//if(!lane)printf("particula: %d,celda: %d\n",particula,celda);
bool done = false;
nv=np_cel[celda];

while(!done){
    unsigned int neighbor;
    unsigned char has_neighbor = 0;
    while(lane >=nv&& !done){
        
    }
}

#pragma unroll
for(int k=lane;k<n_de_cel_vec;k+=chp){
if(ccll>cuantas_celdas_por_hilo)break;
celv=cel_vec[celda*n_de_cel_vec+k];
for(int jp=0;jp<np_cel[celv];jp++){
j = p_en_cel[celv*nmax_p_en_cel+jp];
pjx = __ldg(p+3*j);
pjy = __ldg(p+3*j+1);
pjz = __ldg(p+3*j+2);
esp2=esp_de_p[j];
pj = InitDataType3<double3,double>(pjx,pjy,pjz);
dif.x=pi.x-pj.x;
dif.y=pi.y-pj.y;
dif.z=pi.z-pj.z;
dif=CondPeriodicas(condper,caja,dif,cajai);
dis=Discuad(dif);
if(dis<=rc*rc&&mad_de_p[particula].x!=mad_de_p[j].x&&M_int[n_esp_p*esp1+esp2]){
nv=atomicInc(&nvec[particula],FULL_MASK);
if(nv<nmaxvec)vecinos[particula*nmaxvec+nv]=j;
}
}
ccll++;
}
}

__global__ void AceleracionesfFuerzasLJVYC(uint np,uint n_esp_p,uint pot_int,uint *esp_de_p,const double *p,int chp,double *M_param,double3 caja,double3 cajai,int3 condper,double *epot,double *a,int *vec,unsigned int *nvec,int nmaxvec,double rc,bool nconf)
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

void SimulacionVYC(uint nc,uint ncp,uint np,uint n_esp_p,uint n_esp_m,uint nparam,uint pot_int,uint max_p_en_esp_mr,
                   uint ensamble,uint termos,uint max_it,uint *esp_de_p,uint *M_int,uint *p_en_m,uint *n_m_esp_mr,uint *n_p_esp_m,
                   int nhilos,int maxhilos,bool vibrante,double rc,double rbuf,double dens,double dt,double kres,
                   double temp_d,double param_termo,double tol,double *param,double *pos,double *vel,double *acel,double *q_rat,
                   double *dis_p_esp_mr_rep,uint3 *mad_de_p,int3 condper,double3 caja,double3 cajai,
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
    double3 tamcel,invtamcel;
    uint nceldas=CrearCeldas(caja,tamcel,invtamcel);
    uint nmax_p_en_cel=CuantasPartEnCel(n_esp_p,nparam,param,tamcel);
    int dc=rc+0.5;
    uint n_de_cel_vec=pow(1+2*dc,3);


    /***************************************************************************************************************************/
    //Memoria para arreglos en la CPU
    uint *h_cel_vec=new uint[n_de_cel_vec*nceldas];

    double *h_M_param=new double [nparam*n_esp_p*n_esp_p];
    double *h_eit=new double[1];
    double s_nh_a = 0.,s_nh_v=0.,s_nh_p=1.,g1=0.;
    
    InicializarMatrizDeParametros(pot_int,n_esp_p,nparam,param,h_M_param);
    CalculoDeCeldasVecinas(caja,h_cel_vec,dc,n_de_cel_vec);
    int nmaxvec=CuantosVecCaben(rcc,param,dens,np,n_esp_p,nparam);
    /***************************************************************************************************************************/
    
    cudaStream_t stream[2];
    for(int i=0;i<2;i++)cudaStreamCreate(&stream[i]);

    //Memoria para arreglos de GPU
    uint *d_esp_de_p,*d_M_int,*d_p_cel,*d_cel_vec,*d_nvec,*d_np_cel;
    int *d_vec;
    uint3 *d_mad_de_p;
    double *d_eit,*d_pos,*d_acel,*d_M_param;

    Errorcuda(cudaMalloc(&d_vec,sizeof(int)*nmaxvec*np),"vec",0);
    
    Errorcuda(cudaMalloc(&d_mad_de_p,np*sizeof(uint3)),"mad_de_p",0);

    Errorcuda(cudaMalloc(&d_esp_de_p,np*sizeof(uint)),"esp_de_p",0);
    Errorcuda(cudaMalloc(&d_M_int,n_esp_p*n_esp_p*sizeof(uint)),"d_M_int",0);
    Errorcuda(cudaMalloc(&d_p_cel,sizeof(uint)*nmax_p_en_cel*nceldas),"p_cel",0);
    Errorcuda(cudaMalloc(&d_cel_vec,sizeof(uint)*n_de_cel_vec*nceldas),"celdas_vecinas",0);
    Errorcuda(cudaMalloc(&d_np_cel,sizeof(uint)*nceldas),"np_cel",0);
    Errorcuda(cudaMalloc(&d_nvec,sizeof(uint)*np),"nvec",0);

    Errorcuda(cudaMalloc(&d_M_param,nparam*n_esp_p*n_esp_p*sizeof(double)),"param",0);
    Errorcuda(cudaMalloc(&d_eit,sizeof(double)*np),"d_eit",0);
    Errorcuda(cudaMalloc(&d_pos,sizeof(double)*nd*np),"d_p",0);
    Errorcuda(cudaMalloc(&d_acel,sizeof(double)*nd*np),"d_a",0);

    //vemos cuantos byte ocupan nuestros arreglos en la gpu y el porcentaje

    size_t memoria_global_utilizada =   (nparam*n_esp_p*n_esp_p+np*(2+2*nd))*sizeof(double)
                                    +   (n_esp_p*n_esp_p+np+nceldas*(1+nmax_p_en_cel+n_de_cel_vec))*sizeof(uint)
                                    +   np*sizeof(uint3)
                                    +   nmaxvec*np*sizeof(int);
    OcupacionDeMemoriaGlobal(memoria_global_utilizada,memoria_global);

    Errorcuda(cudaMemcpy(d_mad_de_p,mad_de_p,np*sizeof(uint3),cudaMemcpyHostToDevice),"mad_de_p",1);
    Errorcuda(cudaMemcpy(d_esp_de_p,esp_de_p,np*sizeof(uint),cudaMemcpyHostToDevice),"esp_de_p",1);
    Errorcuda(cudaMemcpy(d_M_int,M_int,n_esp_p*n_esp_p*sizeof(uint),cudaMemcpyHostToDevice),"M_int",1);

    Errorcuda(cudaMemcpy(d_M_param,h_M_param,sizeof(double)*nparam*n_esp_p*n_esp_p,cudaMemcpyHostToDevice),"M_param",1);
    Errorcuda(cudaMemcpy(d_pos,pos,sizeof(double)*nd*np,cudaMemcpyHostToDevice),"pos",1);
    Errorcuda(cudaMemcpy(d_cel_vec,h_cel_vec,sizeof(uint)*n_de_cel_vec*nceldas,cudaMemcpyHostToDevice),"celdas_vecinas",1);

    /***************************************************************************************************************************/

    Errorcuda(cudaMemset(d_nvec,0,sizeof(int)*np),"memset",0);
    Errorcuda(cudaMemset(d_np_cel,0,sizeof(int)*nceldas),"memset",0);
    CalculoCeldas<<<blockspergridcel,threadsperblockcel>>>(np,nmax_p_en_cel,d_np_cel,d_p_cel,d_pos,invtamcel,caja,nceldas);
    Errorcuda(cudaGetLastError(),"Calculo de celdas",3);
    cudaDeviceSynchronize();
    CalculoDeVecinosConCeldas<<<blockspergridvec,threadsperblockvec>>>(np,n_de_cel_vec,nmax_p_en_cel,n_esp_p,d_esp_de_p,d_cel_vec,d_np_cel,d_p_cel,d_M_int,d_nvec,chp,nmaxvec,d_vec,rcc,d_pos,condper,d_mad_de_p,caja,cajai,invtamcel);
    Errorcuda(cudaGetLastError(),"Calculo de vecinos",3);
    cudaDeviceSynchronize();

    AceleracionesfFuerzasLJVYC<<<blockspergrid,threadsperblock>>>(np,n_esp_p,pot_int,d_esp_de_p,d_pos,chp,d_M_param,caja,cajai,condper,d_eit,d_acel,d_vec,d_nvec,nmaxvec,rc,nconf);
    Errorcuda(cudaGetLastError(),"PLJ",3);
    Errorcuda(cudaMemcpy(acel,d_acel,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);
    if(max_p_en_esp_mr-1)PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,condper,mad_de_p,kres,pos,acel,dis_p_esp_mr_rep,caja,cajai);

    Reduccionconwarps<1024><<<1,1024,0,stream[1]>>>(d_eit,np,1);//Cuando es energia es true/1 calcula energia potencial, cuando es false/0 calcula la cinetica
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
        if(nc%10==0){
            cudaMemset(d_np_cel,0,sizeof(int)*nceldas);
            cudaMemset(d_nvec,0,sizeof(int)*np);
            CalculoCeldas<<<blockspergridcel,threadsperblockcel>>>(np,nmax_p_en_cel,d_np_cel,d_p_cel,d_pos,invtamcel,caja,nceldas);
            Errorcuda(cudaGetLastError(),"Calculo de celdas",3);
            cudaDeviceSynchronize();
            CalculoDeVecinosConCeldas<<<blockspergridvec,threadsperblockvec>>>(np,n_de_cel_vec,nmax_p_en_cel,n_esp_p,d_esp_de_p,d_cel_vec,d_np_cel,d_p_cel,d_M_int,d_nvec,chp,nmaxvec,d_vec,rcc,d_pos,condper,d_mad_de_p,caja,cajai,invtamcel);
            Errorcuda(cudaGetLastError(),"Calculo de vecinos",3);
            cudaDeviceSynchronize();
        }
         AceleracionesfFuerzasLJVYC<<<blockspergrid,threadsperblock>>>(np,n_esp_p,pot_int,d_esp_de_p,d_pos,chp,d_M_param,caja,cajai,condper,d_eit,d_acel,d_vec,d_nvec,nmaxvec,rc,nconf);
        Errorcuda(cudaGetLastError(),"PLJ",3);
        Errorcuda(cudaMemcpy(acel,d_acel,sizeof(double)*nd*np,cudaMemcpyDeviceToHost),"a",2);
        if(max_p_en_esp_mr-1)PotencialesDeRestriccion(n_esp_m,max_p_en_esp_mr,n_m_esp_mr,n_p_esp_m,p_en_m,condper,mad_de_p,kres,pos,acel,dis_p_esp_mr_rep,caja,cajai);

        Velocidades(np,vel,acel,dt);
        //debido a que se necesitan F(t+dt) entonces RATTLEvel se realiza aquí
        if(!vibrante&&(max_p_en_esp_mr-1)){
            RattleVel(max_it,n_esp_m,max_p_en_esp_mr,n_p_esp_m,n_m_esp_mr,p_en_m,mad_de_p,tol,pos,vel,dis_p_esp_mr_rep);
        }

        if(ensamble==1){
            ect = CalculoEnergiaCinetica(np,vel);
            g1=(ect*np-3*np*temp_d)/param_termo;
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
            ect = CalculoEnergiaCinetica(np,vel);
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
            IAaD(np,ofasat,pos,vel,acel);
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
    ect=CalculoEnergiaCinetica(np,vel);

    std::cout << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    std::cout << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;

    ofasres << std::endl;
    ofasres <<"Resultados finales"<<std::endl;
    ofasres << "nc,temp,dens,etp,ecp,eip,dtt"<<std::endl;
    ofasres << nc << " " << temp << " " << dens << " " << etp << " " << ecp << " " << eip << " " << dtt << std::endl;
    

    IAaD(np,ofasat,pos,vel,acel);
    
      
    for(int i=0;i<2;i++)cudaStreamDestroy(stream[i]);
    
    Errorcuda(cudaFree(d_vec),"vec",4);
    
    Errorcuda(cudaFree(d_mad_de_p),"mad_de_p",4);

    Errorcuda(cudaFree(d_esp_de_p),"esp_de_p",4);
    Errorcuda(cudaFree(d_M_int),"d_M_int",4);
    Errorcuda(cudaFree(d_p_cel),"p_cel",4);
    Errorcuda(cudaFree(d_cel_vec),"celdas_vecinas",4);
    Errorcuda(cudaFree(d_np_cel),"np_cel",4);
    Errorcuda(cudaFree(d_nvec),"nvec",4);

    Errorcuda(cudaFree(d_M_param),"param",4);
    Errorcuda(cudaFree(d_eit),"d_eit",4);
    Errorcuda(cudaFree(d_pos),"d_p",4);
    Errorcuda(cudaFree(d_acel),"d_a",4);
    delete[] h_cel_vec;
    delete[] h_eit;
    cudaDeviceReset();
}
#endif
