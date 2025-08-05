#ifndef PROP_GPU_HEADER
#define PROP_GPU_HEADER

#include <iostream>


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
    double porcentaje_memoria=(static_cast<double>(memoria_arreglos)/memoria_global)*100;
    printf("Memoria Global utilizada en arreglos: %zu\nPorcentaje de la memoria global total: %lf\n",memoria_arreglos,porcentaje_memoria);
    printf("-----------------------------------------------\n");
}

size_t DetectarMemSh(){
    int deviceCount = 0;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);
    //primero detectar si hay device
    if(err != cudaSuccess){
    printf("cudaGetDeviceCount returned %d\n-> %s\n",
           static_cast<int>(err), cudaGetErrorString(err));
    printf("Result = FAIL\n");
    exit(EXIT_FAILURE);
    }
    cudaSetDevice(0);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    size_t shmem = deviceProp.sharedMemPerBlock;
    return (shmem);
}

#endif
