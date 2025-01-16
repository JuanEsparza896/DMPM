#include <iostream>
#include <string>
#include "SistemaInicial.hpp"
#include "OperacionesTDatosCuda.cuh"
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
    int nd;                 //dimension en la que se trabaja
    /***********************************************************************/
    //Datos de atomos

    uint np;
    //pot ya no sirve, ahora necesitamos la matriz de interacciones
    //np necesita no ser dato en el documento de texto, viene de 
    //eps y sig varian dependiendo de la especie atomica 
    int pot;
    double *p,*v,*a;
    double eps,sig;
    
    
    



    //Inicializacion del sistema
    /*******************************************************************************************************************/
    dir="/home/gach/DMPM";
    condper=InitDataType3<int3>(1,1,1);
    PropiedadesGPU(maxhilos,memoria_global);
    LeerDatos(dir,dens,nd,np,nc,ncp,dt,temp,v0,rc,pot,eps,sig,dpsco,ofaedi,ofapin,opt,nhilos);
    p=new double[np*nd];
    v=new double[np*nd];
    a=new double[np*nd];
    Cuadrada(np,nd,sig,caja,dens,ofapin,p);
    double3 dcajai=InvDataType3<double3>(caja);
    ImprimirDatos(dens,nd,np,pot,sig,eps,caja,nc,ncp,dt,rc,ofaedi);
    VelocidadesInicialesalAzar(v0,v,np,nd);
    /*******************************************************************************************************************/
    double rbuf=0.5;
    ArchivosDeResultados(dpsco,ofasres,ofasat,opt);
    switch(opt)
    {
        case 0:
        
        Simulacion(np,nd,p,v,a,sig,eps,caja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,nhilos,pot,maxhilos,memoria_global);
        break;
        case 1:
        SimulacionV(np,nd,p,v,a,sig,eps,caja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,rc,rbuf,nhilos,pot,maxhilos,memoria_global);
        break;
        case 2:
        SimulacionC(np,nd,p,v,a,sig,eps,caja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,rc,rbuf,nhilos,pot,maxhilos,memoria_global);
        break;
        case 3:
        SimulacionVYC(np,nd,p,v,a,sig,eps,caja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,rc,rbuf,nhilos,pot,maxhilos,memoria_global);
        break;
    }
    delete[] p;
    delete[] v;
    delete[] a;
    return 0;
}
