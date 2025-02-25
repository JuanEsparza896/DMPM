#ifndef POT_REST_HEADER
#define POT_REST_HEADER

#include "Definiciones.cuh"
#include "Funcionescompartidas.h"
/*void RestriccionLongitudEnlace1(uint np,uint n_esp_m,uint max_p_en_esp_mr,uint kres,uint *m_de_esp_mr,uint3 condper,uint3 *mad_de_p,double *pos,double *acel,double3 caja,double3 cajai,double3 *pos_respecto_p_central)
{
    double3 dx;
    double dis,pot=0.;
    uint esp;
    uint *moleculas_de_especie_acumuladas=new uint[n_esp_m];
    for(int i=0;i<n_esp_m;i++)moleculas_de_especie_acumuladas[i]=m_de_esp_mr[i];

    for(int i=0;i<np;i++){
        esp=0;
        while(mad_de_p[i].x>=moleculas_de_especie_acumuladas[esp])esp++;
        //por como esta hecho el loop siempre se salta el caso donde ip no pertenezca a una molecula ya que i-mad[i].y=i+mad[i].x=i
        for(int j=i-mad_de_p[i].y;j<i+mad_de_p[i].z;j++){
            if(i!=j){
                dx.x=pos[i*nd]-pos[j*nd];
                dx.y=pos[i*nd+1]-pos[j*nd+1];
                dx.z=pos[i*nd+2]-pos[j*nd+2];
                if(condper.x)dx.x-=signoF(dx.x*cajai.x)*caja.x;
                if(condper.y)dx.y-=signoF(dx.y*cajai.y)*caja.y;
                if(condper.z)dx.z-=signoF(dx.z*cajai.z)*caja.z;
                dx.x-=pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[i].y].x-pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[j].y].x;
                dx.y-=pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[i].y].y-pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[j].y].y;
                dx.z-=pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[i].y].z-pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[j].y].z;
                dis=dx.x*dx.x+dx.y*dx.y+dx.z*dx.z;
                pot+=0.5*kres*dis;
                acel[i*nd]-=kres*dx.x;
                acel[i*nd + 1]-=kres*dx.y;
                acel[i*nd + 2]-=kres*dx.z;
            }
        }
    }
}*/

void RestriccionLongitudEnlace1(uint inicio,uint final,uint max_p_en_esp_mr,uint esp,uint *p_en_m,uint3 *mad_de_p,int3 condper,double kres,double *pos,double *acel,double3 *pos_respecto_p_central,double3 caja,double3 cajai)
{
    uint part;
    double3 dx;
    double dis,pot=0.;
    for(uint i=inicio;i<inicio+final;i++){
        part=p_en_m[inicio];
        for(int j=part-mad_de_p[i].y;j<part+mad_de_p[i].z;j++){
            if(i!=j){
                dx.x=pos[part*nd]-pos[j*nd];
                dx.y=pos[part*nd+1]-pos[j*nd+1];
                dx.z=pos[part*nd+2]-pos[j*nd+2];
                if(condper.x)dx.x-=signoF(dx.x*cajai.x)*caja.x;
                if(condper.y)dx.y-=signoF(dx.y*cajai.y)*caja.y;
                if(condper.z)dx.z-=signoF(dx.z*cajai.z)*caja.z;
                dx.x-=pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[part].y].x-pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[j].y].x;
                dx.y-=pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[part].y].y-pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[j].y].y;
                dx.z-=pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[part].y].z-pos_respecto_p_central[max_p_en_esp_mr*esp+mad_de_p[j].y].z;
                dis=dx.x*dx.x+dx.y*dx.y+dx.z*dx.z;
                pot+=0.5*kres*dis;
                acel[part*nd]-=kres*dx.x;
                acel[part*nd + 1]-=kres*dx.y;
                acel[part*nd + 2]-=kres*dx.z;
            }
        }
    }
}

void RestriccionAngulo1(uint inicio,uint final,uint *p_en_m)
{
    uint part1,part2,part3;
    for(uint i=inicio;i<inicio+final;i++){
        part1=p_en_m[inicio];

    }
}

void PotencialesDeRestriccion(uint n_esp_m,uint *m_de_esp_mr,uint *M_int_int,uint max_p_en_esp_mr, uint esp, uint *p_en_m, uint3 *mad_de_p, int3 condper, double kres, double *pos, double *acel, double3 *pos_respecto_p_central, double3 caja, double3 cajai)
{
    //Pronto se modificará esta rutina para tener un archivo donde se debe especificar que particulas en la especie molecular deben tener que restriccion, por ejemplo, si tenemos 4 particulas
    //pero solo particula 1-2 2-3 1-4 mantienen distancia de enlace fijo y 1-2-3 mantienen angulo, eso debe ser especificado.

    //Por ahora esto funciona donde todas las particulas en la molecula tienen un potencial restrictivo de longitud de enlace (Similar al algoritmo SHAKE/RATTLE)
    uint *primer_indice_de_mol=new uint[n_esp_m];
    primer_indice_de_mol[0]=0;
    for(int i=1;i<n_esp_m;i++)primer_indice_de_mol[i]=primer_indice_de_mol[i-1]+m_de_esp_mr[i];
    for(uint i=0;i<n_esp_m;i++){
        switch (M_int_int[i*n_esp_m])
        {
        case 0:
            break;
        
        case 1:
            RestriccionLongitudEnlace1(primer_indice_de_mol[i],m_de_esp_mr[i],max_p_en_esp_mr,i,p_en_m,mad_de_p,condper,kres,pos,acel,pos_respecto_p_central,caja,cajai);
            break;
        }
        switch (M_int_int[i*n_esp_m+1])
        {
        case 0:
            break;
        
        case 1:
            RestriccionAngulo1(primer_indice_de_mol[i],m_de_esp_mr[i],p_en_m);
            break;
        }
    }
}
#endif