# Funciones

## Reduccionconwarps

La razón para usar template es por comodidad pues se asegura que al momento de escribir código el número de hilos que se solicitan y el tamaño del arreglo es el mismo.

````
template<int Numthreads>
````

Ejemplo:

````
Reduccionconwarps<1024><<<1,1024>>>(d_eit,np,1)
````

La función es void porque los kernel son tipo void, no usan return.

````
__global__ void Reduccionconwarps(double *arr,int nelementos,bool energia)
````

La información se pasa de memoria global a memoria estática compartida porque el acceso es mas rápido.

````
 __shared__ double s_mem[Numthreads];
````

Para asegurar que toda la información ya termino de copiarse a memoria compartida es necesario sincronizar los hilos en el bloque 

```
__syncthreads()
```

Para realizar el algoritmo de reducción se realizan las sumas con operaciones de shuffle y el primer hilo del warp es el que contiene el resultado de la suma, 
ese resultado se guarda en emoria global y se vuelve a hacer reducción por warps con los resultados acumulados de los demás warps.

````
while(faltan>0){
        for(int i=16;i>=1;i/=2)suma+=__shfl_down_sync(0xffffffff,suma,i);
        faltan/=32;
        if(!lane)s_mem[warp]=suma;
        __syncthreads();
        suma=s_mem[gid];
    }
````

El hilo 0 del bloque es el que contiene el resultado correcto, los demás acumulan errores debido a como funcionan las operaciones de shfl.

````
if(!gid)arr[0]=suma;
````

## CondPeriodicas

rint toma un número tipo double y lo transforma a su entero más cercano, debido a la relación entre dx y cajai el resultado que se espera es que rint devuelva -1,0,1 por lo que funciona como una función signo.
