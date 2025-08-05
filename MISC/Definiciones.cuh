#ifndef DEF_HEADER
#define DEF_HEADER

/*
El MINBLOCKPERGRID se puede cambiar dependiendo del sistema, en las simulaciones de Lennard-Jones
de partículas 4 bloques como mínimo daba muy buenos resultados

Cuando se aplica FULL_MASK en las funciones de shuffle sync a las direcciones de memoria, estas no cambian,
ya que FULL_MASK es puros 1 en binario así que al aplicar el operador binario AND la dirección a la que se 
le aplica no cambia

nd es el número de dimensiones

En ErrorCUDA se pueden tener distintos tipos de Error, la lista de cada error con
su número correspondiente es la siguiente

Asignación de memoria   ->  0
copia de datos h->d     ->  1
copia de datos d->h     ->  2
ejecucion de kernel     ->  3
liberar memoria         ->  4
memset                  ->  5

La API es el primer argumento, el nombre del arreglo es el segundo y el número 
correspondiente al tipo de error es el tercero.
*/

#include <string>
using str = std::string;                    
#define MINBLOCKPERGRID 4                   
#define FULL_MASK 0xffffffff                
#define nd 3                    
#define LennardJones 1
#define nparamLJ 2
#define Yukawa 2
#define nparamYkw 3
#define CARTESIANAS 0
#define POLARES 1

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

#endif
