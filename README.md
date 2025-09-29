# DMPM (Dinámica Molecular en Paralelo para Moléculas)

Contenido:

* Descripción general.
* Requisitos de instalación.
* Creando una simulación.
  * Simular partículas.
  * Simular Moléculas.

## Descripción General

DMPM es un programa de simulación de dinámica molecular cuyo objetivo es introducir a la gente a el desarrollo de software científico en paralelo.

El programa está desarrollado en C con ciertas utilidades de C++, al ser para GPU se utiliza CUDA para el uso de GPU.

El número de opciones para simular es reducido para que el contenido a explorar no sea abrumador, el programa cuenta con funciones necesarias para poder generar simulaciones simples en ensambles NVE y NVT.

Los archivos que contienen las funciones para la simulación se encuentran en 2 carpetas, DM y MISC; en la primera carpeta se encuentran las funciones que se asocian directamente a simulación de partículas y en la segunda carpeta se encuentran las funciones que corresponden a operaciones matemáticas y computacionales por ejemplo la detección de GPU y sus propiedades.

## Requisitos de instalación

El programa funciona en Linux.

 * CUDA TOOLKIT
   
  Se puede descargar en: https://developer.nvidia.com/cuda-toolkit
    
  La guía de instalación: https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html


