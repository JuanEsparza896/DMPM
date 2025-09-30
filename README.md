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

### Archivo Datos_Sistema

El valor de nem y nea tiene que ser el mismo ya que el número de especies de moléculas y de partículas es el mismo
despues de eso se indica el número de partículas que hay de cada especie y para cada especie de molécula i naem\[i\] vale 1

Un ejemplo para un sistema con 3 especies de partículas con 100 de la primera especie, 50 de la segunda y 15 de la tercera:

      3 nem
      3 nea
      100 nme[0]
      50 nme[1]
      15 nme[2]
      1 naem[0]
      1 naem[1]
      1 naem[2]
         
No es necesario que del lado derecho de los números se ponga este texto específicamente, se puede hacer el archivo de la siguiente manera:

      3 a
      3 b
      100 nme
      50 nme
      15 nme
      1 naem
      1 naed
      1 naeas

Lo que si es necesario es que sean una sola cadena de texto. 

## Simular moléculas



