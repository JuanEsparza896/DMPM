# Funciones

## Funcion CuantosVecCaben

Tipo de función: int

Ejecución: CPU

Parámetros:

* double rc: Radio de corte.
* double *param: Se usa para saber el diámetro de cada especie de partícula.
* double dens: Densidad del sistema.
* uint np: Número de partículas
* uint n_esp_p: Número de especies de partículas
* uint nparam: Dimensión de param.

Descripción:

A partir del tamaño del radio de corte y las partículas de menor tamaño se aproxima el número máximo de vecinos que
puede tener una partícula.

Información extra:

https://en.wikipedia.org/wiki/Close-packing_of_equal_spheres

## Funcion CalculoDeVecinos

Tipo de función: void (kernel)

Ejecución: GPU

Parámetros:

* uint np: Número de partículas.
* uint n_esp_p: Número de especies de partículas.
* uint *M_int: Matríz de interacción.
* uint *esp_de_p: Lista con las especies de cada partícula.
* int chp: Cuantos hilos se encargan de los cálculos de una partícula.
* const double *p: Posiciones de las partículas.
* int *vecinos: Arreglo con los vecinos de cada partícula.
* uint *n_vecinos: Arreglo con el número de vecinos de cada partícula.
* int nmaxvec: Número máximo de vecinos de cada partícula.
* int3 condper: Bandera para revisar si las condiciones periódicas están activas.
* uint3 *mad_de_p: Lista de que partículas están en una molécula.
* double3 caja: Dimensiones de la caja de simulación.
* double3 cajai: Inverso de las dimensiones de la caja de simulación.
* double rc: radio de corte.

Descripción:

Para cada partícula se calcula la distancia con las demás y si la distancia es menor que el radio de corte,
no son de la misma molécula y la matríz de interacción indica que ambas especies si interactuan, se guarda a la ora partícula en 
la lista de vecinos y el número de vecinos incrementa en 1.

Información extra:

NA

## Funcion CuantasPartEnCel

Tipo de función: uint

Ejecución: CPU

Parámetros:

* uint n_esp_p: Número de especies de partícula.
* uint nparam: Tamaño de param.
* double *param: Lista con los diámetros de las especies de partícula
* double3 tamcel: dimensiones de la celda.

Descripción:

Se toman las partículas más pequeñas y se calcula cuantas de estas caben en una celda de tamaño tamcel.

Información extra:

https://doi.org/10.37236/1786

## Funcion CrearCeldas

Tipo de función: void

Ejecución: CPU

Parámetros: 

* double3 caja
* double3 tamcel
* double3 invtamcel

Descripción:

Divide la caja de simulación en celdas de tamaño casi unitario, y obtiene sus dimensiones y su inverso.

Información extra:

NA

## Funcion CalculoDeCeldasVecinas

Tipo de función: void

Ejecución: CPU

Parámetros:

* double3 caja: Dimensiones de la caja de simulación
* uint *veccel: Arreglo de celdas vecinas.
* short dc: Con cuantas celdas en cada dirección se interactua.
* uint nveccel: Con cuantas vecinas interacciona otra.

Descripción:

Cada celda tiene su identificador, las partículas en una celda siempre van a interactuar con las partículas de las celdas vecinas,
para evitar hacer un loop triple para buscar la celda con la que se interactua a continuación se guardan mejor en un arreglo.

Información extra:

NA

## Funcion CalculoCeldas

Tipo de función: void (kernel)

Ejecución: GPU

Parámetros: 

* uint np: Número de partículas.
* int nmax_particulas_en_celda: Máximo número de partículas en una celda.
* uint *num_particulas_en_celda: Arreglo con cuántas partículas hay en una celda.
* uint *particulas_en_celda: Arreglo con las partículas en la celda.
* double *p: Arreglo con las componentes de las posiciones de las partículas.
* double3 invtamcel: inverso de las dimensiones de la celda.
* double3 caja: tamaño de la caja de simulación.
* uint nceldas: número de celdas.

Descripción:

Cada celda se puede encontrar con 3 coordenadas de donde se encuentra cierta arista de la celda y sus dimensiones, de esa forma podemos averiguar en que
celda está cada partícula a partir de sus coordenadas.

Información extra:

NA

## Funcion CalculoDeVecinosConCeldas

Tipo de función: void (kernel)

Ejecución: GPU

Parámetros:

* uint np: Número de partículas.
* uint n_de_cel_vec: Número de celdas vecinas.
* uint nmax_p_en_cel: Máximas partículas en una celda.
* uint n_esp_p: Número de especies de partículas
* uint *esp_de_p: Arreglo con las especies de las partículas.
* uint *cel_vec: Arreglo con las celdas vecinas.
* uint *np_cel: Arreglo con el número de partículas en una celda.
* uint *p_en_cel: Arreglo con las partículas en cierta celda.
* uint *M_int: Matríz de interacción.
* uint *nvec: Número de vecinos.
* int chp: Cuántos hilos hacen los cálculos de una partícula.
* int nmaxvec: Número máximo de vecinos.
* int *vecinos: Lista de vecinos.
* double rc: Radio de Corte.
* double *p: Coordenadas de la posición de las partículas.
* int3 condper: Identificador de si las condiciones periódicas están activas o no.
* uint3 *mad_de_p: partículas en una molécula.
* double3 caja:  Dimensiones de la caja de simulación.
* double3 cajai: Inverso de las dimensiones de la caja de simulación.
* double3 invtamcel: Inverso de las dimensiones de la celda.

Descripción:

Pasa lo mismo que con el algoritmo de vecinos pero solo se considera el cálculo de distancias con las partículas que estén en celdas que se encuentren dentro del radio de corte.

Información extra:

NA
