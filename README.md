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
