#ifndef HILOSYBLOQUES_HEADER
#define HILOSYBLOQUES_HEADER
//En esta version lo que hacemos es que el numero de hilos en un bloque se balancee con el numero de bloques, dificilmente usamos los 1024 hilos por bloque
void HilosenBloqueMultiplodeWarp(int np, int &hpp,int &blockspergrid,int &threadsperblock,int minbpg,int maxtpb,int wz=32)
{
    if(np<hpp)hpp=np;//casi nunca pasa a menos que hagamos 31 particulas o menos 
    int hilos_totales=np*hpp;
    blockspergrid=(hilos_totales<=maxtpb)?1:hilos_totales/maxtpb;
    if(hilos_totales%maxtpb!=0)blockspergrid++;
    if(blockspergrid<minbpg)blockspergrid=minbpg;
    threadsperblock=hilos_totales/blockspergrid;
    int redondeo=threadsperblock/wz;
    if(threadsperblock%wz!=0)redondeo++;
    threadsperblock=wz*redondeo;
    printf("-----------------------------------------------\n");
    printf("Hilos requeridos: %d\nHilosTotales: %d\nHilosPorBloque: %d\nBloquesPorMalla %d\nHilosPorParticula: %d\n",hilos_totales,blockspergrid*threadsperblock,threadsperblock,blockspergrid,hpp);
    printf("-----------------------------------------------\n");
}

//Redondea el numero de hilos que hacen calculos de una particula a una potencia de 2 donde la potencia mayor es 5, es decir, a lo mucho 32 hilos pueden participar
int HilosPorParticula(int nhilos)
{
    int chp=nhilos;
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);
    int logmaxhilos=log2(prop.warpSize);
    int comp=log2(chp);
    if(comp>logmaxhilos){
        comp=logmaxhilos;
        printf("numero maximo de hilos que pueden realizar calculos de una particula: %d\n",prop.warpSize);    
    }
    if(chp!=pow(2,comp))comp++;
    chp=pow(2,comp);
    return chp;
}

void PropiedadesGPU(int &maxhilosporbloque, size_t &memoria_global)
{
    printf("\nPropiedades del device:\n");
    cudaDeviceProp prop;
    cudaGetDeviceProperties(&prop,0);
    printf("Device:%s\n",prop.name);
    maxhilosporbloque=prop.maxThreadsPerBlock;
    printf("Compute Capability: %d.%d\n",prop.major,prop.minor);
    printf("tamano de un warp: %d\n",prop.warpSize);
    printf("Numero de hilos pro bloque maximo: %d\n",maxhilosporbloque);
    memoria_global=prop.totalGlobalMem;
    printf("Memoria global total %zu bytes\n\n",memoria_global);
}

void OcupacionDeMemoriaGlobal(size_t memoria_arreglos,size_t memoria_global)
{
    double porcentaje_memoria=(memoria_arreglos/memoria_global)*100;
    printf("Memoria Global utilizada en arreglos: %zu\nPorcentaje de la memoria global total: %lf\n",memoria_arreglos,porcentaje_memoria);
    printf("-----------------------------------------------\n");
}
#endif