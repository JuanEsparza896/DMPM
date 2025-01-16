#ifndef SIS_INICIAL_HEADER
#define SIS_INICIAL_HEADER
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
using str=std::string;
#define PAbrioArchivo(arch,farch) if(!farch){std::cout<<"error al abrir archivo "<<arch<<std::endl;exit(EXIT_FAILURE);}
#define POTLJ 1
#define POTYK 2

void LeerDatosSistema1(str dir,uint &n_esp_m,uint &n_esp_p)
{
    str f=dir+"/Datos/Datos_Sistema.txt";
    str temp;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> n_esp_m; iff >> temp;
    iff >> n_esp_p; iff >> temp;
    
}
void LeerDatosSistema2(str dir,uint n_esp_m,uint n_esp_p,uint *n_m_esp_mr,uint *n_p_esp_m,uint &np)
{
    str f=dir+"/Datos/Datos_Sistema.txt";
    str temp;
    uint imp;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> imp; iff >> temp;
    iff >> imp; iff >> temp;
    np=0;
    for(int i=0;i< n_esp_m;i++){
        iff >> imp; iff >> temp;
        n_m_esp_mr[i]=imp;
    }
    for(int i=0;i< n_esp_m;i++){
        iff >> imp; iff >> temp;
        n_p_esp_m[i]=imp;
        np+=n_m_esp_mr[i]*imp;
    }
}
void LeerArchivosIniciales(str dir,double &dens,int &nd,uint &np)
{
    str f=dir+"/Datos/DatosIniciales.txt";
    str temp;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> dens; iff >> temp;
    iff >> nd; iff >> temp;
    iff >> np; iff >> temp;
    iff.close();
}
/*
Definicion de las variables
nc numero de configuraciones que se exploran
ncp porcentaje de nc para el cual se imprimen  propiedades macroscopicas
dt tamaño del paso de integracion
temp temperatura del sistema en caso de usar NVT
v0 maximo de rapidez inicial de las particulas
rc a partir de que distancia se hacen las interacciones no se consideran
*/
void LeerDatosCorrida(str dir,uint &nc,uint &ncp,double &dt,double &temp,double &v0,double &rc,double &dens,int &pot,int &cvec,int &ccel,int &nhilos,bool &vibrante)
{
    str f=dir+"/Datos/Datos_Corrida.txt";
    str tp;
    std::ifstream iff(f);
    uint tmp;
    PAbrioArchivo(f,iff);
    iff >> nc; iff >> tp;
    iff >> ncp; iff >> tp;
    iff >> dt; iff >> tp;
    iff >> temp; iff >> tp;
    iff >> v0; iff >> tp;
    iff >> rc; iff >> tp;
    iff >> dens; iff >> tp;
    iff >> pot; iff >> tp;
    iff >> nhilos; iff >> tp;
    iff >> cvec; iff >> tp;
    iff >> ccel; iff >> tp;
    iff >> nhilos; iff >> tp;
    iff >> tmp; iff >> tp;
    if(tmp)vibrante=true; 
    iff.close();
}
/*
Definicion de las variables
sig,eps son parametros para el potencial de Lennard-Jones
*/
void LeerDatosAtomos(str dir,double &eps,double &sig)
{
    str f=dir+"/Datos/Datos_Interaccion.txt";
    str tp;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> eps; iff >> tp;
    iff >> sig; iff >> tp;
    iff.close();
}

void AbrirArchivos(str dir,str &dpsco,double dens,int nd,uint np,int pot,double param1,double param2,std::ofstream &ofaedi,std::ofstream &ofapin)
{
    //Creando Carpetas y abriendo archivos para escribir en ellos 
    str dpsc1 = dir + "/Corridas";
    //creando el directorio de corridas
    mkdir(dpsc1.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    str dpsc = dpsc1 +"/dens_"+ std::to_string(dens) + "_nd_" + std::to_string(nd) + "_np_" + std::to_string(np);
    
    switch (pot)
    {
    case POTLJ:
        dpsc+="_sig_"+ std::to_string(param1)+"_eps_"+ std::to_string(param2);
        break;
    case POTYK:
        /* code */
        break;

    }
    std::cout << "LOS ARCHIVOS DE ESTA CORRIDA SE GUARDAN EN: " << dpsc << std::endl;
    //creando el directorio de psc
    mkdir(dpsc.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    dpsco = dpsc+"/Resultados";
    mkdir(dpsco.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    str aedi = dpsc + "/DatosIniciales.txt";
    str apin = dpsc + "/PosIni.txt";

    ofaedi.open(aedi.c_str());
    ofapin.open(apin.c_str());
}

void Cuadrada(uint np,int nd,double sig,double3 &caja,double dens,std::ofstream &ofapin,double *p)
{
    //configuracion inicial cuadrada
    int ndivx=0,ndivy=0,ndivz=0;
    int ndiv=pow(double(np),1.0/double(nd))+0.5;

    double celx=sig,cely=sig,celz=sig;
    double fcx=celx/2,fcy=cely/2,fcz=celz/2;
    ndivx=ndiv;
    int prod=0;
    prod=np/ndivx;
    int ndiv2=0;
    ndiv2=pow(prod,1.0/(nd-1))+0.5;
    ndivy=ndiv2;
    if(nd>2){
        prod=ndivx*ndivy;
        ndivz=np/(ndivx*ndivy);
        if(np%prod!=0){
            ndivz++;
        }
    }
    else{
        ndivz=1;
    }
    caja.x = ndivx*celx;
    caja.y = ndivy*cely;
    caja.z = ndivz*celz;
    double di=0.0;
    di=pow(dens,1.0/nd);
    caja.x/=di;
    caja.y/=di;
    caja.z/=di;
    celx=caja.x/ndivx;
    cely=caja.y/ndivy;
    celz=caja.z/ndivz;

    std::cout <<"\ndimensiones de la caja: (" << caja.x <<","<< caja.y<<","<< caja.z <<")"<< std::endl;
    switch(nd){
        case 2: std::cout << "tamano de la caja: " << caja.x*caja.y << std::endl;
        break;
        case 3: std::cout << "tamano de la caja: " << caja.x*caja.y*caja.z << std::endl;
        break;
    }
    ofapin <<"ip"<< "\t\t" << "px:" << "\t\t" << "py:" << "\t\t" << "pz:" << std::endl;
    int ip=0;
    int x=0,y=0,z=0;
    for(ip=0;ip<np;ip++){
        if(x>=ndivx){
            x=0;
            y++;
        }
        if(y>=ndivy){
            y=0;
            z++;
        }
        p[nd*ip]=x*celx+fcx;
        p[nd*ip+1]=y*cely+fcy;
        p[nd*ip+2]=z*cely+fcz;
        x++;
        ofapin << ip << "\t\t" << p[nd*ip]<< "\t\t" <<p[1+nd*ip]<< "\t\t"<<p[2+nd*ip] << std::endl;
    }
}

void ImprimirDatos(double dens,int nd, uint np, int pot,double param1, double param2,double3 caja,int nc,uint ncp,double dt,double rc,std::ofstream &ofaedi)
{
    std::cout << "\nPROPIEDADES DEL SISTEMA\n";
    ofaedi << "\nPROPIEDADES DEL SISTEMA\n";
    std::cout << "-----------------------------------------------\n";
    ofaedi << "-----------------------------------------------\n";
    std::cout << "Densidad: " << dens << std::endl;
    ofaedi << "Densidad: " << dens << std::endl;
    std::cout << "Dimensiones: " << nd << std::endl;
    ofaedi << "Dimensiones: " << nd << std::endl;
    std::cout << "Numero de atomos: " << np << std::endl;
    ofaedi << "Numero de atomos: " << np << std::endl;
    std::cout << "Numero de configuraciones: " <<nc<< std::endl;
    ofaedi << "Numero de configuraciones: " <<nc<< std::endl;
    std::cout << "Porcentaje de nc para el cual imprimimos propiedades macroscopicas: " << ncp << std::endl;
    ofaedi << "Porcentaje de nc para el cual imprimimos propiedades macroscopicas: " << ncp << std::endl;
    std::cout << "Paso de integracion: " <<dt<< std::endl;
    ofaedi << "Paso de integracion: " <<dt<< std::endl;
    std::cout << "radio de corte: " << rc << std::endl;
    ofaedi << "radio de corte: " << rc << std::endl;
    switch (pot)
    {
    case POTLJ:
        std::cout << "Parametros del potencial de LJ" <<std::endl;
        ofaedi << "Parametros del potencial de LJ" <<std::endl;
        std::cout << "diametro: " << param1 << std::endl;
        ofaedi << "diametro: " << param1 << std::endl;
        std::cout << "profundidad del pozo de potencial: " << param2 << std::endl;
        ofaedi << "profundidad del pozo de potencial: " << param2 << std::endl;
        break;
    
    case POTYK:
        break;
    }
    ofaedi << "Dimensiones de la caja de simulacion: ("<< caja.x<<","<< caja.y<<","<< caja.z<<")\n";
    std::cout << "-----------------------------------------------\n";
    ofaedi << "-----------------------------------------------\n";
}
void LeerDatos(str dir,double &dens,int &nd,uint &np,uint &nc,uint &ncp,double &dt,double &temp,double &v0,double &rc,int &pot,double &eps,double &sig,str &dpsco,std::ofstream &ofaedi,std::ofstream &ofapin,int &opt,int &nhilos)
{
    int cvec=0,ccel=0;
    LeerArchivosIniciales(dir,dens,nd,np);
    LeerDatosCorrida(dir,nc,ncp,dt,temp,v0,rc,dens,pot,cvec,ccel,nhilos,vibrante);
    opt=cvec+2*ccel;
    switch (pot)
    {
    case POTLJ:
        LeerDatosLJ(dir,eps,sig);
        break;
    
    case POTYK:
        break;
    }
    AbrirArchivos(dir,dpsco,dens,nd,np,pot,sig,eps,ofaedi,ofapin);
}
void VelocidadesInicialesalAzar(double v0,double *v,int np,int nd)
{
    //velocidades Iniciales al azar
    double r;
    srand(1);
      for(int ip=0; ip<np; ip++)
      {
          for(int id=0;id<nd;id++){
            r=(rand() % 10000 /10000.)*2.0-1.0;
            v[id+ip*nd] = v0*r;
          }
      }

}

#endif
