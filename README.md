# DMPM (Dinámica Molecular en Paralelo para Moléculas)

Contenido:

* [Descripción general.](#descripción-general)
* [Requisitos de instalación.](#requisitos-de-instalación)
* [¿Cómo funciona el programa?.](¿cómo-funciona-el-programa?)
* [Simular partículas.](#simular-partículas)
* [Simular Moléculas.](#simular-moléculas)

## Descripción General

DMPM es un programa de simulación de dinámica molecular cuyo objetivo es introducir a la gente a el desarrollo de software científico en paralelo.

El programa está desarrollado en C con ciertas utilidades de C++, al ser para GPU se utiliza CUDA para el uso de GPU.

El número de opciones para simular es reducido para que el contenido a explorar no sea abrumador, el programa cuenta con funciones necesarias para poder generar simulaciones simples en ensambles NVE y NVT.

Los archivos que contienen las funciones para la simulación se encuentran en 2 carpetas, DM y MISC; en la primera carpeta se encuentran las funciones que se asocian directamente a simulación de partículas y en la segunda carpeta se encuentran las funciones que corresponden a operaciones matemáticas y computacionales por ejemplo la detección de GPU y sus propiedades.

main.cu toma funciones definidas en las carpetas de DM y MISC y muestra la estructura de una simulación de dinámica molecular, para poder elegir que tipo de sistema se simula es necesario darle parámetros en el sistema, los cuales se leen de la carpeta Datos.

## Requisitos de instalación

El programa probado en distintas distribuciones de Ubuntu.

* 20.04
* 22.04
* 24.04

 * CUDA TOOLKIT
   
  Se puede descargar en: https://developer.nvidia.com/cuda-toolkit
    
  La guía de instalación: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

* Compilador de C++ requisito también del CUDATOOLKIT

## ¿Cómo funciona el programa?

Cuando se descarga el repositorio, el archivo main.cu ya tiene la estructura de una simulación de dinámica molecular, a continuación se explora que significa cada linea de código de este archivo:

### \# include

Es una directiva de C/C++ basicamente pide que se incluyan los contenidos de los archivos/librerías que se escriben inmediatamente después de la directiva, a continuación se mencionan las librerías y su razón para ser incluidas.

1. iostream.
   
   En el código se utilizan las APIs de cout la cual es una alternativa a printf() más amigable ya que no requiere saber el formato de las variables que se imprimen en terminal. 
   
2. string
 
   Permite el acceso a las variables de la clase string, las cuales son una opción alterna al uso de char
    
3. DM/...

    Todos los archivos que pertenecen a la carpeta DM son funciones que realizan distintos procesos de la simulación

### main()

Dentro de la función main() se tienen primero todas las variables que se necesitan para la simulación, dentro del archivo se puede explorar que hace cada una de estas variables.

### std::string dir

Esta variable contiene la ruta del archivo main.cu

### Inicialización del sistema

Para poder inicializar las variables en main() se leen distintos archivos que se encuentran en la carpeta Datos, todas las funciones previas a la función AbrirArchivos se encargan de asignar valores a las variables y arreglos del sistema.

### AbrirArchivos() e Impresión de Datos

Además de escribir los resultados en terminal el programa crea una carpeta de resultados donde se encuentran archivos tipo txt con la información de la simulación. La escritura a los archivos se realiza en las funciones de impresión

### Configuración inicial

Para este programa la única configuración inicial aceptada es una de cristal cúbica.

### Inicializar velocidades

Las velocidades se inicializan a partir del principio de equipartición de la energía y una temperatura deseada.

### Distancias en reposo

En caso de trabajar con moléculas es necesario tener algoritmos de constricción que nos permitan preservar la estructura, la forma más simple es manteniendo las distancias iniciales entre partículas, las cuales son guardadas en otro arreglo.

### Simulación

Se tienen 4 tipos de simulación los cuales dependen del tipo de optimizaciones que se quieran ocupar, estas son vecinos cercanos y celdas, la razón por la cuál no se conserva solo aquella con todas las optimizaciones es por motivos de enseñanza; para mostrar como cambian los algoritmos de simulación cada vez que se añade una optimización.

## Simular partículas

No es necesario cambiar el archivo main.cu para realizar simulaciones, todos los cambios se realizan en la carpeta de datos, nótese que realizar dos simulaciones con los mismos parámetros no genera un archivos distinto.

### Archivo Datos_Corrida.txt

Se muestra un ejemplo de lo que podría contener este archivo y a continuación se explica que significa.

      nc          1000
      ncp         10
      dt          0.001
      temp        2.5
      dens        0.65
      potencial   1
      ensamble    1
      termo       1
      p_termo     0.1
      vibrante    0
      rc          3.5
      ovec        1
      ocel        1
      hilosPP     32 
      coord       0

Estos parámetros corresponden a un archivo que explora 1,000 configuraciones (nc), imprime propiedades termodinámicas cada 10\% de la corrida (ncp) con un paso de integración (dt) de 0.001, la temperatura (temp) inicial del sistema es 2.5, la densidad (dens) es 0.65, el potencial de interacción es Lennard-Jones (potencial),el ensamble es NVT (ensamble), el termostato es el de Andersen (termo), su parámetro correspondiente es de 0.9 (paramtermo) , el radio de corte en caso de usar vecinos cercanos y/o celdas es de 3.5 (rc), el algoritmo usa optimización de vecinos (ovec) y de celdas (ocel), se utilizan 32 hilos de la GPU para calcular fuerzas y contribuciones a la energía interna de una partícula con el resto y las coordenadas están en un sistema cartesiano (coord) , vibrante solo se ocupa si se simulan moléculas por lo que su valor no es importante en este contexto.

Para saber a que potencial de interacción corresponde cada número se puede revisar el archivo Definiciones.cuh en la carpeta MISC, para coord = 1 se leen las coordenadas como si el sistema de referencia estuviera en coordenadas esféricas, cuando ensamble es 0 el sistema es NVE, termo puede tomar valores de 0 a 4 y sus termostatos correspondientes son reescalamiento de velocidades, andersen, berendsen, Bussi-Donadio-Parinello y Nosé-Hoover, vibrante se mencionará más adelante.

[Volver a simular moléculas](#simular-moléculas)

### Archivo Datos_Sistema.txt

El valor de n_esp_m y n_esp_p tienen que ser el mismo ya que el número de especies de moléculas y de partículas es el mismo,
despues de eso se indica el número de moléculas que hay de cada especie molecular (n_m_esp_mr) pero como las especies moleculares son las mismas que las de partículas no hay problema, en la siguiente linea se pone cuantas partículas tiene cada especie molecular (n_p_esp_m) como son moléculas de una partícula se escribe 1 para todas

Un ejemplo para un sistema con 3 especies de partículas con 100 de la primera especie, 50 de la segunda y 15 de la tercera:
      
      n_esp_m     1
      n_esp_p     3
      n_m_esp_mr  100 50 15
      n_p_esp_m   1 1 1
         
Es necesario que del lado izquierdo de los números se ponga este texto específicamente, ya que para la lectura de datos se busca que un string coincida con lo que está escrito para asignarle a la variable con ese nombre ese valor.

[Volver a simular moléculas](#simular-moléculas)

### Archivo Datos_Constricciones.txt

No es relevante para la simulación de partículas.

### Archivo Datos_Interaccion.txt

En este archivo se coloca una matríz simétrica de n x n donde n es el número de especies de partículas, el valor de una entrada de la matríz es 0 si las especies no interactuan y 1 si lo hacen, un ejemplo para 3 especies de partículas donde misma especie interacciona, 0 con 1 si, 0 con 2 no, 1 con 2 no es:

      1 1 0
      1 1 0
      0 0 1

[Volver a simular moléculas](#simular-moléculas)

### Archivo Datos_Particulas.txt

Cada renglón del archivo corresponde a una especie de partícula, la primera columna corresponde a su diámetro y las demás a los parámetros del potencial de interacción.

[Volver a simular moléculas](#simular-moléculas)

### Archivo Datos_Moléculas

No es relevante para la simulación de partículas

### Realizando la simulación

El programa se compila desde terminal, una vez que se localiza el archivo main.cu se compila de la siguiente manera:

      nvcc main.cu

Lo anterior genera un ejecutable con nombre a.out, en caso de cambiar el nombre del ejecutable esto es lo que se escribe en terminal

      nvcc main.cu -o ejecutable

para que se ejecute el programa en terminal se escribe 

      ./ejecutable

Inmediatamente el programa se empieza a ejecutar, los resultados estarán en la carpeta DMPM dentro de una carpeta llamada resultados.

**Recordar que el archivo dir tiene que tener la ruta del archivo main.cu**

[Volver a simular moléculas](#simular-moléculas)

## Simular moléculas

### Archivo Datos_Corrida.txt moléculas

Pasa lo mismo que el [archivo de partículas](#archivo-datos_corridatxt) la diferencia es que vibrante es una variable relevante, cuando vibrante es 0 o False el algoritmo para mantener la estructura de las moléculas es RATTLE

### Archivo Datos_Sistema.txt moléculas

La idea es similar a la del [archivo de partículas](#archivo-datos_sistematxt), ahora los elementos de naem pueden ser distintos de 1, además nem y nea pueden tener valores distintos.

### Archivo Datos_Constricciones.txt moléculas

En el archivo se encuentra lo siguiente:
         
      RATTLE
      tol     0.0001
      maxit   10000.0
      VIBRANTE
      kres    10000.0

Si vibrante es true entonces se leerá kres la cuál es la constante elástica que conecta a las partículas de una molécula para preservar su estructura.
Si vibrante es false entonces se usa el algoritmo de RATTLE, el cual tiene 2 parámetros, la tolerancia y el número máximo de iteraciones.

### Archivo Datos_Interaccion.txt moléculas

Pasa lo mismo que para el [caso de partículas](#archivo-datos_interacciontxt).

### Archivo Datos_Moleculas.txt moléculas

Se muestra un ejemplo de archivo para un sistema con 3 especies moleculares, donde la especie 0 tiene 2 partículas, de especies 0 de partículas, ambas, la especie 1 contiene 3 partículas con especies de partículas 0,2,1 y la especie 2 contiene 1 partícula de especie 2 de partícula.

      0 0 0.0 0.0 0.0
      0 0 0.0 0.0 1.0
      1 0 0.0 0.0 0.0
      1 2 0.0 0.0 0.5
      1 1 0.0 0.0 -0.5
      2 2 0.0 0.0 0.0

Como se observa la primera columna indica la especie de las moléculas, la segunda la especie de la partícula, la tercera, cuarta y quinta las coordenadas de dicha partícula respecto a una partícula central.

**Las entradas de este archivo deben venir en orden por el número de la especie de las moléculas**

### Corriendo el programa

Es el mismo proceso que para [partículas](#realizando-la-simulación). 


