#ifndef TERMOSTATOS_HEADER
#define TERMOSTATOS_HEADER
#include "Definiciones.cuh"

void NoseHoverTermo(unsigned int np,double &s_nh_a,double &s_nh_v,double &s_nh_p,double temp,double m_s_nh,double *vel)
{
    double g1;
    for(int i=0;i<np;i++)for(int id=0;i<nd;i++)g1 += vel[i*nd+id] * vel[i*nd+id];
    g1 -= (3*np + 1)*temp;
    s_nh_a = s_nh_v*s_nh_v/s_nh_p+g1*s_nh_p/m_s_nh;
}

#endif