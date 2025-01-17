#ifndef DEF_HEADER
#define DEF_HEADER
#include <string>
using str = std::string;
#define MINBLOCKPERGRID 4
#define FULL_MASK 0xffffffff
#define nd 3
#define LennardJones 1
#define Yukawa 2
//rama de memoria compartida
/*
Tipos de errores de cuda:

Asignación de memoria   ->  0
copia de datos h->d     ->  1
copia de datos d->h     ->  2
ejecucion de kernel     ->  3
liberar memoria         ->  4
memset                  ->  5
*/
#define Errorcuda(err,vectortemp,error)  if (err != cudaSuccess) {\
    printf("CUDA_ERROR): No fue posible ");\
    switch(error){\
        case 0: fprintf(stderr, "asignar memoria para %s (error code %s)!\n",\
        vectortemp,cudaGetErrorString(err));\
        break;\
        case 1: fprintf(stderr, "copiar %s host a device (error code %s)!\n",\
        vectortemp,cudaGetErrorString(err));\
        break;\
        case 2: fprintf(stderr, "copiar %s device a host (error code %s)!\n",\
        vectortemp,cudaGetErrorString(err));\
        break;\
        case 3: fprintf(stderr, "lanzar %s kernel (error code %s)!\n",\
        vectortemp,cudaGetErrorString(err));\
        break;\
        case 4: fprintf(stderr, "iberar memoria de %s (error code %s)!\n",\
        vectortemp,cudaGetErrorString(err));\
        break;\
        case 5: fprintf(stderr, "inicializar arreglo: %s (error code %s)!\n",\
        vectortemp,cudaGetErrorString(err));\
        break;\
    }\
    exit(EXIT_FAILURE);\
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