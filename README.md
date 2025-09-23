# DMPM (Dinámica Molecular en Paralelo para Moléculas)

Contenido:

* Requisitos de instalación.
* Contenidos del programa.
* Creando una simulación.
  * Simular partículas.
  * Simular Moléculas.

## Requisitos de instalación

El programa funciona en Linux.

 * CUDA TOOLKIT
   
  Se puede descargar en: https://developer.nvidia.com/cuda-toolkit
    
  La guía de instalación: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html

## Contenidos del programa.

En esta sección se describen cada una de las carpetas y archivos que pertenecen a DMPM.
Cada una de estas carpetas tiene un documento de texto que corresponden a cada archivo de código donde se explican a detalle cada una de las funciones correspondientes. 

### Carpeta DM
    
En esta carpeta se encuentran todos los archivos con código para las rutinas de dinámica molecular, a continuación se hablan de cada uno de estos archivos y el tipo de funciones que contienen.

* FuncCompSim.cuh

  Este archivo contiene funciones que se usan recurrentemente en los procesos de dinámica molecular, como el el cálculo de posición en caso de condiciones periódicas, obtención de distancia y obtención de la energía total del sistema a partir de un arreglo.

* Optimizaciones.cuh

  Este archivo contiene las funciones relacionadas con las optimizaciones de celdas y vecinos cercanos en dinámica molecular en GPU, como son, el cálculo máximo del número de vecinos posibles de una partículas, el máximo número de partículas en una celda y lo cálculos de los vecinos y asignación de celdas para cada partícula en GPU.

* Potenciales.cuh

  Contiene las rutinas donde se calculan las fuerzas y contribuciones a la energía potencial asociados a diversos modelos de interacción entre partículas, además de las propiedades promedio asociadas a esos modelos.
 
  ### Carpeta Datos
  
  ### Carpeta MISC
  
  ### main.cu
  

  
