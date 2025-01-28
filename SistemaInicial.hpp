#ifndef SIS_INICIAL_HEADER
#define SIS_INICIAL_HEADER
#include <string>
#include <iostream>
#include <fstream>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <iomanip>
#include <sstream>
#include "Definiciones.cuh"
#include "OperacionesTDatosCuda.cuh"
using str=std::string;
#define PAbrioArchivo(arch,farch) if(!farch){std::cout<<"error al abrir archivo "<<arch<<std::endl;exit(EXIT_FAILURE);}

void LeerDatosSistema1(str dir,uint &n_esp_m,uint &n_esp_p)
{
    str f=dir+"/Datos/Datos_Sistema.txt";
    str temp;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> n_esp_m; iff >> temp;
    iff >> n_esp_p; iff >> temp;
    iff.close();
}
void LeerDatosSistema2(str dir,uint n_esp_m,uint n_esp_p,uint *n_m_esp_mr,uint *n_p_esp_m,uint &np,uint &nm)
{
    str f=dir+"/Datos/Datos_Sistema.txt";
    str temp;
    uint imp;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> imp; iff >> temp;
    iff >> imp; iff >> temp;
    np=nm=0;
    for(int i=0;i< n_esp_m;i++){
        iff >> imp; iff >> temp;
        n_m_esp_mr[i]=imp;
        nm+=imp;
    }
    for(int i=0;i< n_esp_m;i++){
        iff >> imp; iff >> temp;
        n_p_esp_m[i]=imp;
        np+=n_m_esp_mr[i]*imp;
    }
    iff.close();
}
void LeerDatosCorrida(str dir,uint &nc,uint &ncp,double &dt,double &temp,double &v0,double &rc,double &dens,int &pot,int &cvec,int &ccel,int &nhilos,bool &vibrante,uint &nparam)
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
    iff >> tmp; iff >> tp;
    if(tmp)vibrante=true;
    iff >> nparam; iff >> tp;
    printf("nparam %d\n",nparam);
    iff.close();
}
void LeerDatosAtomos(str dir,double *param,uint n_param,uint n_esp_p)
{
    /*esta rutina asume que siempre el dato que esta en la primera columna de este archivo es el diametro de la especie atomica */
    str f=dir+"/Datos/Datos_Atomos.txt";
    str tp;
    double temporal;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> tp;iff >> tp;
    for(int j=0;j<n_esp_p;j++)
        for(int i=0;i<n_param;i++){
            iff >> temporal;
            param[n_param*j+i]=temporal;    
        }
    iff.close();
}
void LeerDatosMoleculas(str dir,uint n_esp_m,uint *n_p_esp_m,uint *esp_p_en_esp_mr,uint max_p_en_esp_mr,double3 *pos_respecto_p_central)
{
    //Nota este arreglo esp_p_en_esp_mr se debe borrar y reemplazar por otro ya que si el maximo y minimo de particulas en especies moleculares es grande se desperdicia mucho espacio
    //Nota2 no es cierto, este arreglo se borra ya que es remplazado por un arreglo de tamaño np donde viene la especie de cada particula
    uint tui;
    str f=dir+"/Datos/Datos_Moleculas.txt";
    str tp;
    double temporal;
    std::ifstream iff(f);
    PAbrioArchivo(f,iff);
    iff >> tp;iff >> tp;iff >> tp;
    for(int i=0;i<n_esp_m;i++){
        iff >> tui;
        for(int j=0;j<n_p_esp_m[i];j++){
            iff >> tui;
            esp_p_en_esp_mr[max_p_en_esp_mr*i+j]=tui;
            iff >> temporal;pos_respecto_p_central[max_p_en_esp_mr*i+j].x=temporal;
            iff >> temporal;pos_respecto_p_central[max_p_en_esp_mr*i+j].y=temporal;
            iff >> temporal;pos_respecto_p_central[max_p_en_esp_mr*i+j].z=temporal;
        }
    }
    iff.close();
}
void LeerDatosInteraccion(str dir,uint n_esp_p,uint *M_int)
{
    str f=dir+"/Datos/Datos_Interaccion.txt";
    str tp;
    std::ifstream iff(f);
    uint tmp;
    PAbrioArchivo(f,iff);
    for(int i=0;i<n_esp_p;i++)
    {
        for(int j=0;j<n_esp_p;j++)
        {
            iff >> tmp;
            M_int[n_esp_p*i+j]=tmp;
            printf("Interaccion de especies %d-%d: ",i,j);
            if(tmp)printf("Activa\n"); else printf("Inactiva\n"); 
        }
    }
    iff.close();
}
void AbrirArchivos(str directorio,double dens, uint nem, uint nea, uint *nme, uint nparam, double *param,std::ofstream &ofaedi,std::ofstream &ofapin,str &dpsco)
{
    /*****************************************/
    std::stringstream stream;
    std::stringstream *stream1;
    stream1=new std::stringstream[nparam*nea];
    str dpsc1,dpsc,s,aedi,apin;
    /*****************************************/
    dpsc1 = directorio + "/Corridas";
    mkdir(dpsc1.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    
    dpsc = dpsc1;
    printf("nparam: %d\n",nparam);
    
    stream << std::fixed << std::setprecision(2) << dens;
    s = stream.str();
    dpsc += "/dens_"+ s;
    dpsc += "_nd_" + std::to_string(nd) + "_nem_" + std::to_string(nem) + "_nea_" + std::to_string(nea);

    for(int i=0;i<nem;i++){
        dpsc+="_nme"+std::to_string(i)+"_"+std::to_string(nme[i]);
    }
    for(int j=0;j<nea;j++){
        dpsc += "_diametro"+std::to_string(j)+"_";
        std::cout << "s " << s<< std::endl;
        stream1[j*nparam] << std::fixed << std::setprecision(2) << param[j*nparam];
        s = stream1[j*nparam].str();
        std::cout << "s " << s<< std::endl;
        dpsc += s;
    }
    for(int j=0;j<nea;j++)
        for(int i=1;i<nparam;i++){
            dpsc += "_param_"+std::to_string(i)+"_"+std::to_string(j)+"_";
            stream1[nparam*j+i] << std::fixed << std::setprecision(2) << param[nparam*j+i];
            std::cout << "s " << s<< std::endl;
            s = stream1[nparam*j+i].str();
            std::cout << "s " << s<< std::endl;
            dpsc += s;    
        }
        printf("bandera interna\n");

    std::cout << "LOS ARCHIVOS DE ESTA CORRIDA SE GUARDAN EN: " << dpsc << std::endl;
    
    mkdir(dpsc.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    dpsco = dpsc+"/Resultados";
    mkdir(dpsco.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);

    aedi = dpsc + "/DatosIniciales.txt";
    apin = dpsc + "/Posiciones_Iniciales.txt";

    ofaedi.open(aedi.c_str());
    ofapin.open(apin.c_str());
}

void ImpresionDeDatos(int nc,int ncp,double dt,double temp,double v0,double rc,bool optvec,
                      bool optcel,int pot,double dens,uint n_esp_m,uint especies_atomicas,
                      uint *m_de_esp_mr,uint *p_en_esp_mr, uint *esp_de_p_en_m,
                      uint max_p_en_esp_mr,double3 *pos_respecto_p_central)
{
    int ncc=nc/ncp;
    printf("\n---------------------------------------------------------------------------\n\n");
    printf("Datos de Corrida:\n\n");
    printf("Configuraciones: %d\nCada cuantas configuraciones imprimimos propiedades:%d\n",nc,ncc);
    printf("Tamaño del paso de integracion: %.3lf\n",dt);
    printf("Temperatura a la que se debe mantener el sistema (CASO NVT): %.1lf\n",temp);
    printf("Rapidez inicial maxima de las particulas: %.1lf\n",v0);
    printf("Radio de corte (en caso de optimizaciones): %.1lf\n",rc);
    printf("Optimizacion de vecinos: ");
    if(optvec)printf("Activa\n");else{printf("Inactiva\n");}
    printf("Optimizacion de celdas: ");
    if(optcel)printf("Activa\n");else{printf("Inactiva\n");}
    printf("Potencial: ");
    switch(pot)
    {
        case LennardJones:printf("Lennard-Jones\n");break;
        case Yukawa:printf("Yukawa\n");break;
    }
    printf("\nDatos del Sistema:\n\n");
    printf("Densidad: %.2lf\n",dens);
    printf("Dimensiones: %d\n",nd);
    printf("\nDatos de Atomos y Moleculas\n\n");
    printf("Especies moleculares: %d\n",n_esp_m);
    printf("Especies atomicas: %d\n\n",especies_atomicas);
    for(int i=0;i<n_esp_m;i++)printf("Moleculas de la especie %d: %d\n",i,m_de_esp_mr[i]);
    for(int i=0;i<n_esp_m;i++)printf("Atomos en la especie molecular %d: %d\n",i,p_en_esp_mr[i]);
    printf("\n");
    for(int i=0;i<n_esp_m;i++)for(int j=0;j<p_en_esp_mr[i];j++)printf("Especie de los atomos en la especie molecular %d: %d\n",i,esp_de_p_en_m[i*max_p_en_esp_mr+j]);
    printf("Respecto al atomo central de la especie molecular correspondiente: \n\n");
    for(int i=0;i<n_esp_m;i++)
        for(int j=0;j<p_en_esp_mr[i];j++)
            printf("Posiciones de los atomos en la especie molecular %d: (%.4lf,%.4lf,%.4lf)\n",i,pos_respecto_p_central[max_p_en_esp_mr*i+j].x,pos_respecto_p_central[max_p_en_esp_mr*i+j].y,pos_respecto_p_central[max_p_en_esp_mr*i+j].z);
    printf("\n---------------------------------------------------------------------------\n\n");
}

void ImpresionDeDatosADisco(int nc,int ncp,double dt,double temp,double v0,double rc,bool optvec,
                            bool optcel,int pot,double dens,uint n_esp_m,uint n_esp_p,
                            uint *m_de_esp_mr,uint *p_en_esp_mr,uint *esp_de_p_en_m,
                            uint max_p_en_esp_mr,double3 *pos_respecto_p_central,std::ofstream &ofaedi)
{
    /**************/
    int ncc=nc/ncp;
    /**************/
    ofaedi << "Datos de Corrida:\n\n";
    ofaedi << "Configuraciones: " << nc << "\nCada cuantas configuraciones imprimimos propiedades:" << ncc << "\n";
    ofaedi << "Tamaño del paso de integracion:" << std::setprecision(3) << dt << "\n";
    ofaedi << "Temperatura a la que se debe mantener el sistema (CASO NVT):" << std::setprecision(2) << temp << "\n";
    ofaedi << "Rapidez inicial maxima de las particulas:" << std::setprecision(2) << v0 << "\n";
    ofaedi << "Radio de corte (en caso de optimizaciones):" << std::setprecision(2) << rc << "\n";
    ofaedi << "Optimizacion de vecinos: ";
    if(optvec)ofaedi << "Activa\n"; else{ofaedi << "Inactiva\n";}
    ofaedi << "Optimizacion de celdas: ";
    if(optcel)ofaedi << "Activa\n";else{ofaedi << "Inactiva\n";}
    ofaedi << "Potencial: ";
    switch(pot)
    {
        case LennardJones:ofaedi << "Lennard-Jones\n";break;
        case Yukawa:ofaedi << "Yukawa\n";break;
    }
    ofaedi << "\nDatos del Sistema:\n\n";
    ofaedi << "Densidad:" << std::setprecision(3) << dens << "\n";
    ofaedi << "Dimensiones:" << nd << "\n";
    ofaedi << "\nDatos de Atomos y Moleculas\n\n";
    ofaedi << "Especies moleculares:" << n_esp_m << "\n";
    ofaedi << "Especies atomicas: " << n_esp_p << "\n\n";
    for(int i=0;i<n_esp_m;i++)ofaedi << "Moleculas de la especie " << i << ": " << m_de_esp_mr[i] << "\n";
    for(int i=0;i<n_esp_m;i++)ofaedi << "Atomos en la especie molecular " << i << ": " << p_en_esp_mr[i] << "\n";
    ofaedi << "\n";
    for(int i=0;i<n_esp_m;i++)for(int j=0;j<p_en_esp_mr[i];j++)ofaedi << "Especie de los atomos en la especie molecular " << i << ": " << esp_de_p_en_m[i*max_p_en_esp_mr+j] << "\n";
    ofaedi << "Respecto al atomo central de la especie molecular correspondiente: \n\n";
    for(int i=0;i<n_esp_m;i++)for(int j=0;j<p_en_esp_mr[i];j++)ofaedi << "Posiciones de los atomos en la especie molecular " << i << ": (" << std::setprecision(4) << pos_respecto_p_central[max_p_en_esp_mr*i+j].x << "," << std::setprecision(4) << pos_respecto_p_central[max_p_en_esp_mr*i+j].y << "," << std::setprecision(4) << pos_respecto_p_central[max_p_en_esp_mr*i+j].z << ")\n";
    
}

void ArchivosDeResultados(str dpsco,std::ofstream &ofasres,std::ofstream &ofasat,int op)
{
    str k;
    switch(op)
    {
        case 0: k="SinOptimizaciones";break;
        case 1: k="Vecinos";break;
        case 2: k="Celdas";break;
        case 3: k="AmbasOptimizaciones";break;
    }
    std::cout << "Optimizacion: " << k << std::endl;
    str asat = dpsco + "/Posiciones_"+k+".txt";
    str asres = dpsco + "/Resultados_"+k+".txt";

    ofasat.open(asat.c_str());
    ofasres.open(asres.c_str());
}

double3 CreandoCeldaMinima(uint n_esp_m,double3 *pos_respecto_p_central,uint max_p_en_esp_mr,
                           double *param, uint nparam,uint *especies_de_atomos_en_molecula,uint *atomos_en_especie_molecular)
{
    
    /*
        [1]

        Por ahora la celda minima se hace de la siguiente manera
        tenemos nd-longitudes de la celda, y esas las separamos en direccion positiva y negativa,
        ahora tenemos 2*nd-longitudes

        Ejemplo de como se asignan las longitudes:
        Queremos calcular la longitud en x positiva (lxp),

        Tomamos la primera especie molecular y comparamos el valor actual de lxp con la posicion 
        de las particulas en x + su radio y si es mayor tenemos una nueva lxp.
        
        Para la longitud negativa en x (lxn),

        En vez de tomar la posicion en x + el radio consideramos posicion en x - radio esto es 
        menor que lxn tenemos nueva lxn.

        !Estas celdas unitarias son mejores cuando las moleculas son monodispersas.
    */

    /***************************************************************************************/
    double3 lon_p,lon_n;
    double3 comparador_pos,comparador_neg;
    int k=0;
    /***************************************************************************************/
    lon_p=InitDataType3<double3,double>(0.0,0.0,0.0);
    lon_n=InitDataType3<double3,double>(0.0,0.0,0.0);
    
    
    for(int i=0;i<n_esp_m;i++){
        for(int j=0;j<atomos_en_especie_molecular[i];j++){
            k=especies_de_atomos_en_molecula[i*max_p_en_esp_mr+j];
            //estas expresiones asumen que el primer parametro para cada especie atomica siempre es su diametro
            comparador_pos=InitDataType3<double3,double>(pos_respecto_p_central[max_p_en_esp_mr*i+j].x+0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].y+0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].z+0.5*param[k*nparam]);

            comparador_neg=InitDataType3<double3,double>(pos_respecto_p_central[max_p_en_esp_mr*i+j].x-0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].y-0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].z-0.5*param[k*nparam]);
            if(comparador_pos.x>=lon_p.x)lon_p.x=comparador_pos.x;
            if(comparador_pos.y>=lon_p.y)lon_p.y=comparador_pos.y;
            if(comparador_pos.z>=lon_p.z)lon_p.z=comparador_pos.z;
            if(comparador_neg.x<=lon_n.x)lon_n.x=comparador_neg.x;
            if(comparador_neg.y<=lon_n.y)lon_n.y=comparador_neg.y;
            if(comparador_neg.z<=lon_n.z)lon_n.z=comparador_neg.z;
        }
    }
    return InitDataType3<double3,double>(lon_p.x-lon_n.x,lon_p.y-lon_n.y,lon_p.z-lon_n.z);

}

void CentrarMoleculas(double *centrar_m,uint n_esp_m,uint *n_p_esp_m,
                      uint *esp_p_en_esp_mr,uint max_p_en_esp_mr,uint nparam,
                      double *param,double3 *pos_respecto_p_central)
{
    
    /*
        [2]

        Supongamos 2 especies moleculares, una tiene todas sus 
        particulas a la derecha respecto a su atomo central, y 
        la otra a la izquierda, entonces, inicialmente podemos 
        tener traslapes, lo que hace esta rutina es recorrer 
        las moleculas de forma que quepan dentro de la celda 
        unitaria.    
    */

    /*********************/
    double3 lon_n;
    double3 comparador_neg;
    int k=0;
    /*********************/
    
    for(int i=0;i<n_esp_m;i++){
        lon_n=InitDataType3<double3,double>(0.0,0.0,0.0);
        for(int j=0;j<n_p_esp_m[i];j++){
            k=esp_p_en_esp_mr[i*max_p_en_esp_mr+j];
            //estas expresiones asumen que el primer parametro para cada especie atomica siempre es su diametro
            comparador_neg=InitDataType3<double3,double>(pos_respecto_p_central[max_p_en_esp_mr*i+j].x-0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].y-0.5*param[k*nparam],
                                                              pos_respecto_p_central[max_p_en_esp_mr*i+j].z-0.5*param[k*nparam]);
            if(comparador_neg.x<=lon_n.x)lon_n.x=comparador_neg.x;
            if(comparador_neg.y<=lon_n.y)lon_n.y=comparador_neg.y;
            if(comparador_neg.z<=lon_n.z)lon_n.z=comparador_neg.z;
        }
        centrar_m[nd*i]=-lon_n.x;
        centrar_m[nd*i+1]=-lon_n.y;
        centrar_m[nd*i+2]=-lon_n.z;
    }
}


void ConfiguracionCubica(uint n_esp_m,uint *m_de_esp_mr,uint *p_en_esp_mr,double *pos,
                         double3 *pos_respecto_p_central,uint max_p_en_esp_mr,double3 &caja,double *centrar_m,
                         double3 celda_minima,double densidad,std::ofstream &ofapin,uint *esp_de_p_en_m,
                         uint nm,uint *esp_de_p)
{
    /********************************************/
    int3 particulas_por_lado;
    double3 cel;
    uint *moleculas_de_especie_acumuladas;
    int k=0,part=0,h=0;
    /********************************************/
    particulas_por_lado.x=pow(nm,1.0/nd)+0.5;
    particulas_por_lado.y=nm/particulas_por_lado.x;
    particulas_por_lado.y=pow(particulas_por_lado.y,1.0/(nd-1))+0.5;
    particulas_por_lado.z=nm/(particulas_por_lado.x*particulas_por_lado.y);
    if(particulas_por_lado.z*particulas_por_lado.x*particulas_por_lado.y<nm)particulas_por_lado.z++;
    printf("\nMoleculas: %d\nMoleculas en cada direccion(inicialmente): (%d,%d,%d)\nMoleculas que caben en la caja de sumulacion %d\n",nm,particulas_por_lado.x,particulas_por_lado.y,particulas_por_lado.z,particulas_por_lado.x*particulas_por_lado.y*particulas_por_lado.z);
    
    cel=InitDataType3<double3,double>(celda_minima.x,celda_minima.y,celda_minima.z);
    
    caja.x = cel.x*particulas_por_lado.x;
    caja.y = cel.y*particulas_por_lado.y;
    caja.z = cel.z*particulas_por_lado.z;
    double di=0.0;
    di=pow(densidad,1.0/nd);
    caja.x/=di;
    caja.y/=di;
    caja.z/=di;
    cel.x=caja.x/particulas_por_lado.x;
    cel.y=caja.y/particulas_por_lado.y;
    cel.z=caja.z/particulas_por_lado.z;
    printf("\nTamano de la caja de simulacion: (%.3lf,%.3lf,%.3lf)\n",caja.x,caja.y,caja.z);

    int x=0,y=0,z=0;

    moleculas_de_especie_acumuladas=new uint[n_esp_m];
    for(int i=0;i<n_esp_m;i++)moleculas_de_especie_acumuladas[i]=m_de_esp_mr[i];
    
    ofapin << "particula molecula especie_molecular especie_atomica"<< std::endl;
    for(int i=0;i<nm;i++){
        if(moleculas_de_especie_acumuladas[k]<=0)k++;
        if(k>=n_esp_m)return;
        if(x>=particulas_por_lado.x){
            x=0;
            y++;
        }
        if(y>=particulas_por_lado.y){
            y=0;
            z++;
        }
        h=0;
        for(int j=0;j<p_en_esp_mr[k];j++){
            pos[part*nd]=x*cel.x+pos_respecto_p_central[max_p_en_esp_mr*k+j].x+centrar_m[nd*k];
            pos[part*nd+1]=y*cel.y+pos_respecto_p_central[max_p_en_esp_mr*k+j].y+centrar_m[nd*k+1];
            pos[part*nd+2]=z*cel.z+pos_respecto_p_central[max_p_en_esp_mr*k+j].z+centrar_m[nd*k+2];
            esp_de_p[part] = esp_de_p_en_m[k*max_p_en_esp_mr+j];
            ofapin << part << "\t" << i << "\t" << k << "\t" << esp_de_p_en_m[k*max_p_en_esp_mr+j]<<
            "\t" << pos[part*nd] << "\t" <<  pos[part*nd+1] << "\t" <<  pos[part*nd+2] << std::endl;
            part++;
            h++;
        }
        x++;
        moleculas_de_especie_acumuladas[k]--;
    }
}

void InicializarVelocidades(double v0,double *v,int np)
{
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

void DistanciasEntreParticulasEnMoleculaIniciales(uint np,uint n_esp_m,uint max_p_en_esp_mr,uint *n_p_esp_mr,double3 *dis_p_esp_mr_rep,double3 *pos_respecto_p_central)
{

    for(int i=0;i<n_esp_m;i++){
        for(int j=0;j<n_p_esp_mr[i];j++){
            for(int k=0;k<n_p_esp_mr[i];k++){
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k].x = pos_respecto_p_central[max_p_en_esp_mr*i+j].x - pos_respecto_p_central[max_p_en_esp_mr*i+k].x;
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k].y = pos_respecto_p_central[max_p_en_esp_mr*i+j].y - pos_respecto_p_central[max_p_en_esp_mr*i+k].y;
                dis_p_esp_mr_rep[i*max_p_en_esp_mr*max_p_en_esp_mr+j*max_p_en_esp_mr+k].z = pos_respecto_p_central[max_p_en_esp_mr*i+j].z - pos_respecto_p_central[max_p_en_esp_mr*i+k].z;        
            }    
        }
    }
    
}

#endif
