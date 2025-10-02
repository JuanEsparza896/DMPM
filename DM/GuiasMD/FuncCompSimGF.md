# Funciones

## Funcion Reduccionconwarps

    template<int Numthreads>
    __global__ void Reduccionconwarps(double *arr,int nelementos,bool energia);

Tipo de función: void (kernel)

Ejecución: en GPU

Parámetros:

* double arr:       Cualquier arreglo con elementos de tipo double.
* int nelementos:   Tamaño del arreglo.
* bool energía:     Cálculo de energía interna (true) y energía cinética (false)

Descripción:

Toma un arreglo en memoria global y usa algoritmos de reducción en paralelo (suma exclusiva con operaciones de shuffle) para obtener la suma de todos los elementos del arreglo.

El arreglo se guarda en memoria compartida por lo que el kernel utilizan un bloque para la ejecución, el parámetro de template se usa para asignar el número de hilos del bloque.

Información extra

https://developer.nvidia.com/blog/faster-parallel-reductions-kepler/
https://developer.nvidia.com/blog/using-cuda-warp-level-primitives/

## Funcion CondPeriodicas

    inline __device__ __host__ double3 CondPeriodicas(int3 condper,double3 caja,double3 dx,double3 cajai)

Tipo de función: double3

Ejecución: CPU o GPU

Parámetros:

* int3 condper: Indica si las condiciones periódicas están activas en las 3 direcciones de movimiento.
* double3 caja: Dimensiones de la caja de simulación.
* double3 dx: Posición de una partícula en la caja de simulación.
* double3 cajai: Cada elemento es el inverso de las dimensiones de la caja, para no calcular divisiones.

Descripción:

Toma la posición de una partícula y encuentra las coordenadas de su imagen que están dentro de la caja de simulación.  
