# Funciones

## CalculoDeVecinos

En esta función las posiciones de las partículas son un arreglo constante, esto es para indicar que se puede acceder a ellas mediante el cache de lectura.

````
const double *p
````

````
pix=__ldg(p+3*particula);
piy=__ldg(p+3*particula+1);
piz=__ldg(p+3*particula+2);
````

Para poder modificar la información de un elemento de un arreglo en memoria global sin tener problemas con la sincronización es mediante operaciones atómicas,
las cuales no son interrumpibles, AtomicInc toma la dirección que se le envía, le aplica la máscara y luego toma el valor que está ahí, le suma 1 y lo guarda en esa
dirección. Mientras este proceso sucede ningun otro hilo puede acceder a esta dirección por lo que debe esperar a que termine la operación, además el nuevo valor lo
guarda en la variable que invoca la operación.

````
nv=atomicInc(&n_vecinos[particula],FULL_MASK);
````
