#ifndef POT_REST_HEADER
#define POT_REST_HEADER

#include "MISC/Definiciones.cuh"
#include "FuncCompSim.cuh"

void RestriccionLongitudEnlace1(uint part,uint j,int3 condper,double kres,double *pos,double *acel,double dist_rep,double3 caja,double3 cajai)
{
    double3 dx;
    double dis,pot=0.,fac;

    dx.x=pos[part*nd]-pos[j*nd];
    dx.y=pos[part*nd+1]-pos[j*nd+1];
    dx.z=pos[part*nd+2]-pos[j*nd+2];

    if(condper.x)dx.x-=signoF(dx.x*cajai.x)*caja.x;
    if(condper.y)dx.y-=signoF(dx.y*cajai.y)*caja.y;
    if(condper.z)dx.z-=signoF(dx.z*cajai.z)*caja.z;

    dis=dx.x*dx.x+dx.y*dx.y+dx.z*dx.z;
    fac=(1.-(dist_rep/dis));
    fac*=kres;
    pot+=0.5*kres*(dis-dist_rep)*(dis-dist_rep);//no es 1/2 pq la contribución al potencial se esta contando 2 veces
    acel[part*nd]-=fac*dx.x;
    acel[part*nd + 1]-=fac*dx.y;
    acel[part*nd + 2]-=fac*dx.z;
    acel[j*nd]+=fac*dx.x;
    acel[j*nd + 1]+=fac*dx.y;
    acel[j*nd + 2]+=fac*dx.z;
}

void PotencialesDeRestriccion(uint n_esp_m,uint max_p_en_esp_mr,uint *n_m_esp_mr,uint *n_p_esp_m,uint *p_en_m,int3 condper,uint3 *mad_de_p,double kres,double *pos,double *acel,double *dis_p_esp_mr_rep,double3 caja,double3 cajai)
{
    /*
    En esta rutina se asume que todas las partículas en una molécula están conectadas mediante resortes
    !Revisar casos límite para el futuro de esta rutina.
    */

    uint mol=0,part=0;
    for(uint e_mol=0;e_mol<n_esp_m;e_mol++){
        if(n_p_esp_m[e_mol]>1) //Esto asegura que las moléculas de una partícula no ejecuten esto
        for(uint n_emol=0;n_emol<n_m_esp_mr[e_mol];n_emol++){
            for(uint i=0;i<n_p_esp_m[e_mol];i++){
                part=p_en_m[mol]+i;
                for(uint j=part+1;j<part+mad_de_p[i].z;j++){
                    if(part!=j)
                    RestriccionLongitudEnlace1(part,j,condper,kres,pos,acel,dis_p_esp_mr_rep[e_mol*max_p_en_esp_mr*max_p_en_esp_mr+i*max_p_en_esp_mr+j-p_en_m[i]],caja,cajai);
                }
            }
            mol++;
        }
    }

}
#endif
