#ifndef OPTIM_HEADER
#define OPTIM_HEADER

#include <stdio.h>
#include "MISC/OperacionesTDatosCuda.cuh"
#include "MISC/Definiciones.cuh"
#include "FuncCompSim.cuh"


/**************************************************************************************************************************************************/
//Vecinos

int CuantosVecCaben(double rc,double *param,double dens,uint np,uint n_esp_p,uint nparam)
{
    double sig=1000.;
    for(int i=0;i<n_esp_p;i++)if(param[i*nparam]<=sig)sig=param[i*nparam];
    double rcc=rc+0.5*sig;
    /*
    La densidad maxima de esferas es:
    https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres
    consideramos cuantas esferas de radio sigma caben en un circulo de radio rcc
    0.75*((rcc*sig)/(0.5*sig))³=6*rcc³
    La expresion anterior es para el numero de esferas completas que caben en ese circulo/esfera
    Para los vecinos solo se considera la distancia entre centros de las particulas por lo que puede
    en el contenedor pueda incluir los centros de particulas que no caben completamente
    por lo que consideraremos mejor un radio del contenedor de rcc+(radio de la particula)
    Ademas debemos multiplicar por la densidad para considerar haya mas separacion inicial entre
    particulas por lo que caben menos
    */
    double calcvec = dens*6*rcc*rcc*rcc;
    int nmaxvec = calcvec;
    if(nmaxvec>=np)nmaxvec=np;
    printf("Numero de particulas %d\nproposicion del numero maximo de vecinos: %d\n",np,nmaxvec);
    printf("-----------------------------------------------\n");
    return nmaxvec;
}

__global__ void CalculoDeVecinos(uint np,uint n_esp_p,uint *M_int,uint *esp_de_p,int chp,const double *p, int *vecinos, unsigned int *n_vecinos, int nmaxvec,int3 condper,uint3 *mad_de_p, double3 caja,double3 cajai,double rc)
{
    int particula = threadIdx.x+blockIdx.x*blockDim.x;
    particula/=chp;
    uint esp1,esp2;
    if(particula>=np)return;
    const int lane= threadIdx.x & (chp-1);
    int nv;
    double pix,piy,piz,pjx,pjy,pjz,dis;
    pix=__ldg(p+3*particula);
    piy=__ldg(p+3*particula+1);
    piz=__ldg(p+3*particula+2);
    esp1=esp_de_p[particula];
    double3 pi=InitDataType3<double3,double>(pix,piy,piz);
    double3 pj,dif;
    for(int j=lane;j<np;j+=chp){
        pjx=__ldg(p+3*j);
        pjy=__ldg(p+3*j+1);
        pjz=__ldg(p+3*j+2);
        esp2=esp_de_p[j];
        pj=InitDataType3<double3,double>(pjx,pjy,pjz);
        dif.x=pi.x-pj.x;
        dif.y=pi.y-pj.y;
        dif.z=pi.z-pj.z;
        dif=CondPeriodicas(condper,caja,dif,cajai);
        dis=Discuad(dif);
        if(dis<=rc*rc&&mad_de_p[particula].x!=mad_de_p[j].x&&M_int[n_esp_p*esp1+esp2]){
            nv=atomicInc(&n_vecinos[particula],FULL_MASK);
            if(nv<nmaxvec)vecinos[particula*nmaxvec+nv]=j;
        }
    }
}

/**************************************************************************************************************************************************/
//Celdas


uint CuantasPartEnCel(uint n_esp_p,uint nparam,double *param, double3 tamcel)
{
    double sig=1000.;
    for(int i=0;i<n_esp_p;i++)if(param[i*nparam]<=sig)sig=param[i*nparam];
    uint num_part=1;
    /*
    para ver cuantas particulas caben en una celda consideraremos celdas cubicas cuyo 
    tamaño esta determinado por la longitud maxima de los distintos lados de la celda
    */
    double sigmax=0.0;
    sigmax=tamcel.x>=sigmax?tamcel.x:sigmax;
    sigmax=tamcel.y>=sigmax?tamcel.y:sigmax;
    sigmax=tamcel.z>=sigmax?tamcel.z:sigmax;
    /*
    https://doi.org/10.37236/1786 contiene un algoritmo para empacar esferas de un 
    mismo tamaño dentro de una caja del tamaño de la unidad en la tabla 2 d_n nos 
    indica la distancia maxima entre las 2 esferas, interpretando esto como el 
    diametro maximo tenemos cuantas particulas de diametro sigma podrían caber en 
    una celda unitaria.
    
    Por como generamos las celdas en este algoritmo la longitud de la celda es muy
    cercana a la unidad (longitud de la caja/entero mas cercano por la izquierda), 
    de todos modos podemos proponer una regla de 3:
        nuevo diametro = diametro reportado*longitud de la celda.
    A partir de la tabla 2 generamos un arreglo para de ahi sacar los valores:
    */

    double a[]=
        {
        2.0,sqrt(3.0),sqrt(2.0),sqrt(2.0),
        sqrt(5.0)/2.0,3*sqrt(2.0)/4.0,
        (-4.0+sqrt(10.0+4.0*sqrt(3.0)))/sqrt(3.0),
        1.0,sqrt(3.0)/2.0,3.0/4.0,0.710116382462,
        0.707106806467,sqrt(2.0)/2.0,sqrt(2.0)/2.0,
        5.0/8.0,0.606667120726,3.0*sqrt(2.0)/7.0,
        sqrt(13.0)/6.0,0.578209612716,0.554761174904,
        3.0/(2.0+2.0*sqrt(3.0)),3*sqrt(2.0)/8.0,
        0.523539214257,0.517638090205,0.505135865094,
        0.501074021252,1.0/2.0,0.471410634842,
        sqrt(2.0)/3.0,sqrt(2.0)/3.0,sqrt(2.0)/3.0,
        };
    int tam_a=sizeof(a)/sizeof(a[0]);
    for(int i=0;i<tam_a;i++)a[i]*=sigmax;

    //aqui comparamos la sigma actual con los valores
    if(sig>a[0]){
        printf("En las celdas cabe 1 particula\n");
        printf("-----------------------------------------------\n");
        return 1;
    }
    if(sig<a[tam_a-1]){
        printf("Las particulas son muy pequeñas\n");
        printf("-----------------------------------------------\n");
        exit(EXIT_FAILURE);
    }
    for(int i=0;i<tam_a;i++){
        if(sig<=a[i])num_part++;
    }
    printf("En la celda caben %d particulas\n",num_part);
    printf("-----------------------------------------------\n");
    return num_part;
}

uint CrearCeldas(double3 caja,double3 &tamcel,double3 &invtamcel)
{
    uint nceldas=0;
    int3 celdas;
    celdas.x=caja.x; printf("Celdas en x: %d\n",celdas.x);
    celdas.y=caja.y; printf("Celdas en y: %d\n",celdas.y);
    celdas.z=caja.z; printf("Celdas en z: %d\n",celdas.z);
    tamcel.x=caja.x/celdas.x;
    tamcel.y=caja.y/celdas.y;
    tamcel.z=caja.z/celdas.z;
    printf("Tamaño de celda minima: (%.4lf,%.4lf,%.4lf)\n",tamcel.x,tamcel.y,tamcel.z);
    printf("-----------------------------------------------\n");
    invtamcel=InvDataType3<double3>(tamcel);
    nceldas=celdas.x*celdas.y*celdas.z;
    return nceldas;
}

//TODO mejorar este algoritmo para la lista vaya desde la celda mas interna a la externa (tratando de redicir error computacional).
void CalculoDeCeldasVecinas(double3 caja,uint *veccel,short dc,uint nveccel)
{
    int ncx=0,ncy=0,ncz=0;
    int celv;
    int cont=0;
    uint eLx=caja.x,eLy=caja.y,eLz=caja.z;
    uint fa=eLx*eLy;
    int cel;

    for(int celz=0;celz<eLz;celz++){
        for(int cely=0;cely<eLy;cely++){
            for(int celx=0;celx<eLx;celx++){
                /********************************************/
                cont=0;
                cel=celx + eLx*cely + fa*celz;
                for(int cz=-dc;cz<=dc;cz++){
                    ncz=celz+cz;
                    if(ncz<0) ncz+=eLz;
                    if(ncz>=eLz) ncz-=eLz;
                    for(int cy=-dc;cy<=dc;cy++){
                        ncy=cely+cy;
                        if(ncy<0) ncy+=eLy;
                        if(ncy>=eLy) ncy-=eLy;
                        for(int cx=-dc;cx<=dc;cx++){
                            ncx=celx+cx;
                            if(ncx<0) ncx+=eLx;
                            if(ncx>=eLx) ncx-=eLx;
                            if(ncx>=0&&ncy>=0&&ncz>=0){
                                celv=ncx + eLx*ncy + fa*ncz;
                                veccel[nveccel*cel+cont]=celv;
                                cont++;
                            }
                        }
                    }
                }
                /********************************************/
            }
        }
    }
}

__global__ void CalculoCeldas(uint np,int nmax_particulas_en_celda,uint *num_particulas_en_celda,uint *particulas_en_celda,double *p,double3 invtamcel,double3 caja,uint nceldas)
{
    int gid=threadIdx.x+blockDim.x*blockIdx.x;
    if(gid>=np)return;
    int fac1=caja.x;
    int fac2=caja.y; fac2*=fac1;
    double px,py,pz;
    int3 pos;
    px=__ldg(p+gid*3);
    py=__ldg(p+gid*3+1);
    pz=__ldg(p+gid*3+2);
    if(px==caja.x&&py==caja.y&&pz==caja.z){px=0.0;py=0.0;pz=0.0;}
    pos.x=px*invtamcel.x;
    pos.y=py*invtamcel.y;
    pos.z=pz*invtamcel.z;
    int celda=pos.x+fac1*pos.y+fac2*pos.z;
    if(celda<0||celda>nceldas)printf("error celda: %d,gid: %d posiciones(%lf,%lf,%lf)\n",celda,gid,px,py,pz);
    int n=atomicInc(&num_particulas_en_celda[celda],FULL_MASK);
    if(n<nmax_particulas_en_celda)particulas_en_celda[celda*nmax_particulas_en_celda+n]=gid;
    if(num_particulas_en_celda[celda]>nmax_particulas_en_celda)printf("exceso de particulas %d en celda %d, n %d, particulas %d\n",gid,celda,n,num_particulas_en_celda[celda]);
}


/**************************************************************************************************************************************************/
//Mezcla

__global__ void CalculoDeVecinosConCeldas(uint np,uint n_de_cel_vec,uint nmax_p_en_cel,uint n_esp_p,uint *esp_de_p,uint *cel_vec,uint *np_cel,uint *p_en_cel,
                                          uint *M_int,uint *nvec,int chp,int nmaxvec,int *vecinos,double rc,double *p,int3 condper,uint3 *mad_de_p,double3 caja,double3 cajai,double3 invtamcel)
{
    int particula,nv,celda,celv,j;
    uint eLx,eLy,fa,esp1,esp2;
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
/**************************************************************************************************************************************************/
#endif
