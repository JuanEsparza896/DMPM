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
    std::string dir="/home/gach/DMP2";
    str dpsco;
    int3 condper=InitDataType3<int3>(1,1,1);
    double dens,dt,temp,v0,rc,eps,sig,cajax,cajay,cajaz;
    int nd,np,nc,ncp,pot,opt,nhilos;
    double *p,*v,*a;
    std::ofstream ofaedi,ofapin,ofasres,ofasat;
    int maxhilos;
    size_t memoria_global;
    PropiedadesGPU(maxhilos,memoria_global);

    //Inicializacion del sistema
    /*******************************************************************************************************************/
    LeerDatos(dir,dens,nd,np,nc,ncp,dt,temp,v0,rc,pot,eps,sig,dpsco,ofaedi,ofapin,opt,nhilos);
    p=new double[np*nd];
    v=new double[np*nd];
    a=new double[np*nd];
    Cuadrada(np,nd,sig,cajax,cajay,cajaz,dens,ofapin,p);
    double3 dcaja=InitDataType3<double3>(cajax,cajay,cajaz);
    double3 dcajai=InvDataType3<double3>(dcaja);
    ImprimirDatos(dens,nd,np,pot,sig,eps,cajax,cajay,cajaz,nc,ncp,dt,rc,ofaedi);
    VelocidadesInicialesalAzar(v0,v,np,nd);
    /*******************************************************************************************************************/
    double rbuf=0.5;
    switch(opt)
    {
        case 0:
        printf("\nNo se usan optimizaciones\n");
        ArchivosDeResultados(dpsco,ofasres,ofasat,"SinOptimizaciones");
        Simulacion(np,nd,p,v,a,sig,eps,dcaja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,nhilos,pot,maxhilos,memoria_global);
        break;
        case 1:
        printf("\nOptimizaciones: Vecinos\n");
        ArchivosDeResultados(dpsco,ofasres,ofasat,"Vecinos");
        SimulacionV(np,nd,p,v,a,sig,eps,dcaja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,rc,rbuf,nhilos,pot,maxhilos,memoria_global);
        break;
        case 2:
        printf("\nOptimizaciones: Celdas\n");
        ArchivosDeResultados(dpsco,ofasres,ofasat,"Celdas");
        SimulacionC(np,nd,p,v,a,sig,eps,dcaja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,rc,rbuf,nhilos,pot,maxhilos,memoria_global);
        break;
        case 3:
        printf("\nOptimizaciones: Vecinosy Celdas\n");
        ArchivosDeResultados(dpsco,ofasres,ofasat,"Vecinos_Celdas");
        SimulacionVYC(np,nd,p,v,a,sig,eps,dcaja,dcajai,condper,temp,ofasres,ofasat,nc,dt,dens,ncp,rc,rbuf,nhilos,pot,maxhilos,memoria_global);
        break;
    }
    delete[] p;
    delete[] v;
    delete[] a;
    return 0;
}
