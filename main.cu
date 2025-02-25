#include <iostream>
#include <string>
#include "SistemaInicial.hpp"
#include "OperacionesTDatosCuda.cuh"
#include "OperacionesDeHilosyBloques.cuh"
#include "SimulacionSinOptimizaciones.cuh"
#include "SimulacionVecinos.cuh"
#include "SimulacionCeldas.cuh"
#include "SimulacionVecinosYCeldas.cuh"


int main()
{
    /***********************************************************************/
    //Datos del programa
    std::string dir;        //directorio base
    str dpsco;              //directorio de resultados
    std::ofstream ofaedi;   //archivo de salida de datos iniciales
    std::ofstream ofapin;   //archivo de salida de posiciones Iniciales
    std::ofstream ofasres;  //archivo de salida de resultados
    std::ofstream ofasat;   //archivo de salida de posiciones de los atomos
    int maxhilos;           //numero maximo de hilos
    size_t memoria_global;  //memoria global disponible
    /***********************************************************************/
    //Datos de la corrida

    int3 condper;           //Condiciones de periodicidad
    int opt;                //optimizaciones vecinos, celdas o ambos
    int nhilos;             //cuantos hilos se encargan de realizar los 
                            //calculos para una particula
    uint nc;                //numero de configuraciones
    uint ncp;               //porcentaje de nc para el cual se calculan props
    double rc;              //radio de corte
    double dt;              //tamaño del paso de integracion
    /***********************************************************************/
    //Datos del sistema
    
    double dens;            //densidad
    double temp;            //caso NVT temperatura del baño
    double v0;              //rapidez maxima inicial de las particulas
    double3 caja;           //tamano de la caja de simulacion
    uint n_esp_m;           //numero de especies moleculares
    uint n_esp_p;           //numero de especies de particulas
    bool vibrante;          //determina que algoritmo de constricciones se usa, potenciales de restriccion o RATTLE
    uint *M_int;            //matriz de interaccion
    double3 celda_min;      //tamano de la celda minima
    /***********************************************************************/
    //Datos de moleculas
    uint nm;                //numero de moleculas
    uint *n_m_esp_mr;       //numero de moleculas de cierta especie molecular
    uint *n_p_esp_m;        //numero de particulas en cierta especie molecular
    uint max_p_en_esp_mr;   //maximo de los valores de n_p_esp_m
    uint *esp_p_en_esp_mr;  //especie de las particulas
    uint *M_int_int;        //matriz de interacciones internas
    uint *p_en_m;           //nos indica la primera particula dentro de cierta molecula 
    /***********************************************************************/
    //Datos de atomos

    uint np;
    uint nparam;
    uint *esp_p;            //arreglo con especie de cada particula
    uint3 *mad_de_p;           //particulas antes y despues de cierta particula en su molecula(Ver nota 1 antes de la inicializacion del sistema)
    uint pot;
    double *pos,*vel,*acel;
    double *param;
    double3 *pos_respecto_p_central;
    
    
    
    //temporales
    uint arr_temp=0;

    /*******************************************************************************************************************
    Notas
    
    Nota 1:
        Este arreglo funciona ya que al generar las moleculas sus particulas se crean consecutivamente, por ejemplo
        si la molecula m tiene 5 particulas, quiere decir que contiene las ip,ip+1,ip+2,ip+3,ip+4
        El arreglo pad_de_p funciona asi, 
            pad_de_p[ip].x a que molecula pertenece la particula
            pad_de_p[ip].y cuantas particulas hay antes de ip en la molecula
            pad_de_p[ip].z cuantas particulas hay despues de ip en la molecula
            
    Inicializacion del sistema
    *******************************************************************************************************************/
    dir="/home/gach/DMPM";
    condper=InitDataType3<int3>(1,1,1);
    PropiedadesGPU(maxhilos,memoria_global);
    LeerDatosSistema1(dir,n_esp_m,n_esp_p);
    n_m_esp_mr = new uint[n_esp_m];
    n_p_esp_m = new uint[n_esp_m];
    LeerDatosSistema2(dir,n_esp_m,n_esp_p,n_m_esp_mr,n_p_esp_m,np,nm);
    pos=new double[np*nd];
    vel=new double[np*nd];
    acel=new double[np*nd];
    p_en_m=new uint[nm];
    int cvec=0,ccel=0;
    LeerDatosCorrida(dir,nc,ncp,dt,temp,v0,rc,dens,pot,cvec,ccel,nhilos,vibrante,nparam);
    opt=cvec+2*ccel;
    param = new double[nparam*n_esp_p];
    LeerDatosAtomos(dir,param,nparam,n_esp_p);
    for(int i=0;i<n_esp_m;i++)
    if(n_p_esp_m[i]>=arr_temp)arr_temp=n_p_esp_m[i];
    max_p_en_esp_mr = arr_temp;arr_temp=0;
    esp_p_en_esp_mr = new uint[max_p_en_esp_mr*n_esp_m];
    pos_respecto_p_central = new double3[max_p_en_esp_mr*n_esp_m];
    LeerDatosMoleculas(dir,n_esp_m,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,pos_respecto_p_central);
    M_int = new uint[n_esp_p*n_esp_p];
    // el 2 es por que ahora solo hay 2 tipos de potenciales de restriccion: de enlace y de angulo, el numero sera 3 cuando se agregue torsion o 4 cuando agregue diedro, etc
    M_int_int = new uint[2*n_esp_m];
    LeerDatosInteraccion(dir,n_esp_p,M_int);
    LeerDatosInteraccionInterna(dir,n_esp_m,M_int_int);
    AbrirArchivos(dir,dens,n_esp_m,n_esp_p,n_m_esp_mr,nparam,param,ofaedi,ofapin,dpsco);
    ImpresionDeDatos(nc,ncp,dt,temp,v0,rc,cvec,ccel,pot,dens,n_esp_m,n_esp_p,n_m_esp_mr,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,pos_respecto_p_central);
    ImpresionDeDatosADisco(nc,ncp,dt,temp,v0,rc,cvec,ccel,pot,dens,n_esp_m,n_esp_p,n_m_esp_mr,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,pos_respecto_p_central,ofaedi);
    celda_min=CreandoCeldaMinima(n_esp_m,pos_respecto_p_central,max_p_en_esp_mr,param,nparam,esp_p_en_esp_mr,n_p_esp_m);
    double *centrar_m = new double[nd*n_esp_m];
    CentrarMoleculas(centrar_m,n_esp_m,n_p_esp_m,esp_p_en_esp_mr,max_p_en_esp_mr,nparam,param,pos_respecto_p_central);
    esp_p = new uint[np];
    mad_de_p = new uint3[np];
    ConfiguracionCubica(n_esp_m,n_m_esp_mr,n_p_esp_m,p_en_m,pos,pos_respecto_p_central,max_p_en_esp_mr,caja,centrar_m,celda_min,dens,ofapin,esp_p_en_esp_mr,nm,esp_p,mad_de_p);
    double3 cajai=InvDataType3<double3>(caja);
    InicializarVelocidades(v0,vel,np);
    /*******************************************************************************************************************/
    //Lo usamos para optimizaciones
    double rbuf=0.5;
    ArchivosDeResultados(dpsco,ofasres,ofasat,opt);
    
    switch(opt)
    {
        case 0:
        Simulacion(nc,ncp,nhilos,np,nparam,pot,n_esp_p,mad_de_p,esp_p,M_int,maxhilos,memoria_global,dt,dens,param,pos,vel,acel,condper,caja,cajai,ofasres,ofasat);
        break;
        case 1:
        SimulacionV(nc,ncp,np,n_esp_p,nparam,pot,esp_p,M_int,nhilos,maxhilos,mad_de_p,condper,rc,rbuf,dens,temp,dt,param,pos,vel,acel,caja,cajai,memoria_global,ofasres,ofasat);
        break;
        case 2:
        SimulacionC(nc,ncp,np,n_esp_p,nparam,pot,esp_p,M_int,nhilos,maxhilos,rc,temp,dt,dens,param,pos,vel,acel,mad_de_p,condper,caja,cajai,memoria_global,ofasat,ofasres);
        break;
        case 3:
        SimulacionVYC(nc,ncp,np,n_esp_p,nparam,pot,esp_p,M_int,nhilos,maxhilos,rc,rbuf,dens,temp,dt,param,pos,vel,acel,mad_de_p,condper,caja,cajai,memoria_global,ofasres,ofasat);
        break;
    }
    
    delete[] pos;
    delete[] vel;
    delete[] acel;
    return 0;
}
