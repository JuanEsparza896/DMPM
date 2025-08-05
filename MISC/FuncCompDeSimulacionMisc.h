#ifndef FUNC_COMP_SIM_MISC
#define FUNC_COMP_SIM_MISC

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>
#include <time.h> 

void IAaD(int np,std::ofstream &ofasat,double *p,double *v,double *a)
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

int signoF(double val){
    int val1;
    if(val>=0.)val1=val+0.5;
    else{val1=val-0.5;}
    return val1;
}

double rand_BM(uint randSeed)
{
    uint num=randSeed;
    if(randSeed==0)num=time(NULL);
    srand(num);
    const double pi=2*acos(0.);
    double u1 = (rand() + 1.0) / (RAND_MAX + 1.0);
    double u2 = (rand() + 1.0) / (RAND_MAX + 1.0);
    return sqrt(-2.*log(u1))*cos(2*pi*u2);
}

#endif
