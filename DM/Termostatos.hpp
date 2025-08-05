#ifndef TERMOSTATOS_HEADER
#define TERMOSTATOS_HEADER
#include "MISC/Definiciones.cuh"
#include "FuncCompSim.cuh"
#include "MISC/FuncCompDeSimulacionMisc.h"
#include <stdlib.h>
#include <time.h>
#include <random>

void RescVel(uint np,double temp_ac,double temp_des,double *vel)
{
    double lambda=0.;
    lambda=sqrt(temp_des/temp_ac);
    for(uint i=0;i<np;i++){
        vel[i*nd]*=lambda;
        vel[i*nd+1]*=lambda;
        vel[i*nd+2]*=lambda;
    }
}

void BerendsenTermo(uint np,double temp_ac,double temp_des,double dt,double tau,double *vel)
{
    double lambda=0.;
    lambda=sqrt(1.+(dt/tau)*((temp_des/temp_ac)-1.));
    for(uint i=0;i<np;i++){
        vel[i*nd]*=lambda;
        vel[i*nd+1]*=lambda;
        vel[i*nd+2]*=lambda;
    }
}

void AndersenTermo(uint np,double nu,double dt,double temp_d,double *vel)
{
    double pro=nu*dt;
    double ran;
    uint num=time(NULL);
    srand(num);
    //fac es el factor que multiplica a la exponencial de la distribución de boltzmann
    //como kb=1.0 y la masa por ahora es 1.0 solo aparece el factor de la temperatura
    const double fac=sqrt(temp_d);
    for (uint i=0;i<np;i++)
    {
        ran=(rand() + 1.0) / (RAND_MAX + 1.0);
        if(ran<pro){
            vel[i*nd]=fac*rand_BM();
            vel[i*nd+1]=fac*rand_BM();
            vel[i*nd+2]=fac*rand_BM();
        }
    }
}

void BDPTermo(uint np,double ec_ac,double ec_d,double dt,double tau,double *vel)
{
    const double gdl = nd*np+1; //Grados de libertad
    const double c1 = exp(-dt/tau);
    const double c2 = ec_d*(1.-c1)/(ec_ac*gdl);
    uint num=time(NULL);
    srand(num);  

    double n_al=rand_BM();
    double n_all=n_al*n_al;
    double var=np-1;
    std::default_random_engine generator;
    std::gamma_distribution<double> distribution(var,1.);
    
    double n_al2=distribution(generator);
    double alpha = c1+(n_all+2*n_al2)*c2+2*sqrt(c1*c2*fabs(n_al));
    alpha = sqrt(alpha);
    for(uint i=0;i<np;i++){
        vel[i*nd]*=alpha;
        vel[i*nd+1]*=alpha;
        vel[i*nd+2]*=alpha;
    }
}

void Termostato(uint termos,uint np,double temp_ac,double temp_des,double dt, double param,double* vel)
{
    //el termostato de Nose-Hoover altera las ecuaciones de movimiento por lo que no se incluye aqui
    /*
    Cada termostato tiene un parámetro que le corresponde, excepto reescalamiento de velocidades
    Andersen:   nu
    Berendsen;  tau
    BDP:        tau
    NH:         M_s
    */
    switch (termos)
    {
    case 0:
        RescVel(np,temp_ac,temp_des,vel);
        break;
    case 1:
    AndersenTermo(np,param,dt,temp_des,vel);
        break;
    case 2:
    BerendsenTermo(np,temp_ac,temp_des,dt,param,vel);
        break;
    case 3:
    BDPTermo(np,temp_ac,temp_des,dt,param,vel);
        break;
    }
}

#endif
