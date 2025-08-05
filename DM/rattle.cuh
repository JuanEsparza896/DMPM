#ifndef RATTLE_HEADER
#define RATTLE_HEADER

#include "../MISC/Definiciones.cuh"
#include "FuncCompSim.cuh"
#include "../MISC/FuncCompDeSimulacionMisc.h"

void RattlePos(uint max_it,uint n_esp_m,uint np,uint max_p_en_esp_mr,uint *n_p_esp_m,uint *n_m_esp_mr,uint *p_en_m,uint3 *mad_de_p,int3 condper,double tol,double dt,double *pos,double *q_rat,double *dis_p_esp_mr_rep,double3 caja,double3 cajai)
{
    double dis,dif=0.,d_max=0.,grr=0.,srij=0.;
    uint mol=0,part=0;
    double3 dx;
    for(uint it=0;it<max_it;it++){
        mol=0;
        //Desde aquí
        for(uint e_mol=0;e_mol<n_esp_m;e_mol++){
            if(n_p_esp_m[e_mol]>1)
            for(uint n_emol=0;n_emol<n_m_esp_mr[e_mol];n_emol++){
                //hasta aquí,
                // lo único que simboliza es que se recorre cada molécula del sistema
                for(uint i=0;i<n_p_esp_m[e_mol];i++){
                    part=p_en_m[mol]+i;
                    for(uint j=part+1;j<part+mad_de_p[i].z;j++){
                        if(part!=j){
                            dx.x=pos[i*nd]+dt*q_rat[i*nd]-pos[j*nd]-dt*q_rat[j*nd];
                            dx.y=pos[i*nd+1]+dt*q_rat[i*nd+1]-pos[j*nd+1]-dt*q_rat[j*nd+1];
                            dx.z=pos[i*nd+2]+dt*q_rat[i*nd+2]-pos[j*nd+2]-dt*q_rat[j*nd+2];
                            dis=Discuad(dx);
                            if(condper.x)dx.x-=signoF(dx.x*cajai.x)*caja.x;
                            if(condper.y)dx.y-=signoF(dx.y*cajai.y)*caja.y;
                            if(condper.z)dx.z-=signoF(dx.z*cajai.z)*caja.z;

                            dif=dis - dis_p_esp_mr_rep[e_mol*max_p_en_esp_mr*max_p_en_esp_mr+i*max_p_en_esp_mr+j-p_en_m[i]]*dis_p_esp_mr_rep[e_mol*max_p_en_esp_mr*max_p_en_esp_mr+i*max_p_en_esp_mr+j-p_en_m[i]];  
                            d_max=fmax(d_max,fabs(dif));
                            if(dif>tol){
                                /*
                                  Si no se cumple la constriccion se busca una corrección g para ri y rj que haga que si se cumpla
                                  la forma de ri_t y rj_t son:
                                  ri_t=ri+dt*(q_rati-g*rij)
                                  rj_t=rj+dt*(q_ratj+g*rij)

                                  para la constricción de distancia se tiene la siguiente g
                                */
                               srij = dx.x*(pos[i*nd]-pos[j*nd])
                                    + dx.y*(pos[i*nd+1]-pos[j*nd+1])
                                    + dx.z*(pos[i*nd+2]-pos[j*nd+2]);
                                //si fueran distintos tipos de constricción aquí el valor de grr cambia, esto está en progreso
                                grr=dif/(2.*dt*(srij));
                                for(uint id=0;id<nd;id++){
                                    q_rat[i*nd+id]-=grr*(pos[i*nd+id]-pos[j*nd+id]);
                                    q_rat[j*nd+id]+=grr*(pos[i*nd+id]-pos[j*nd+id]);
                                }                         
                            }             
                        }
                    }
                }
                mol++;
            }
        }
        if(d_max<tol)break;
    }
    
}

void RattleVel(uint max_it,uint n_esp_m,uint max_p_en_esp_mr,uint *n_p_esp_m,uint *n_m_esp_mr,uint *p_en_m,uint3 *mad_de_p,double tol,double *pos,double *vel,double *dis_p_esp_mr_rep)
{
    uint mol=0,part;
    double d_max=0.,dot,ka=0.;
    double3 dx;
    for(uint it=0;it<max_it;it++){
        mol=0;
        for(uint emol=0;emol<n_esp_m;emol++){
            if(n_p_esp_m[emol]>1)
            for(uint n_emol=0;n_emol<n_m_esp_mr[emol];n_emol++){
                for(uint i=0;i<n_p_esp_m[emol];i++){
                    part=p_en_m[mol]+i;
                    for(uint j=part+1;j<part+mad_de_p[i].z;j++){
                        if(part!=j){
                            dx.x = (pos[i*nd]-pos[j*nd])*(vel[i*nd]-vel[j*nd]);
                            dx.y = (pos[i*nd+1]-pos[j*nd+1])*(vel[i*nd+1]-vel[j*nd+1]);
                            dx.z = (pos[i*nd+2]-pos[j*nd+2])*(vel[i*nd+1]-vel[j*nd+2]);
                            dot=Discuad(dx);
                            d_max=fmax(d_max,fabs(dot));
                            if(dot>tol){
                                ka=dot/(dis_p_esp_mr_rep[emol*max_p_en_esp_mr*max_p_en_esp_mr+i*max_p_en_esp_mr+j-p_en_m[i]]*dis_p_esp_mr_rep[emol*max_p_en_esp_mr*max_p_en_esp_mr+i*max_p_en_esp_mr+j-p_en_m[i]]);
                                for(uint id=0;id<nd;id++){
                                    vel[i*nd+id]-=ka*(pos[i*nd+id]-pos[j*nd+id]);
                                    vel[j*nd+id]+=ka*(pos[i*nd+id]-pos[j*nd+id]);
                                }
                            }
                        }
                    }
                }
                mol++;
            }
        }
        if(d_max<tol)break;
    }

}


#endif
