# Ecuación-Schrödinger
El proyecto consiste en un solucionador númerico de la ecuación de Schrödinger para un pozo finito unidimensional estacionario, el cual fue desarrollado en Linux, utilizando lenguaje C. 

## Requisitos
- Linux
- gcc
- make

## Dependencias
Posee las siguientes dependencias:
- GSL (GNU Scientific Library) (Probado con 2.7.1)

### Pasos para instalar la biblioteca de GSL
En Ubuntu, Debian:
`sudo apt install libgsl-dev`

En Arch Linux:
`sudo pacman -S gsl`

En Fedora:
`sudo dnf install gsl-devel`

## Pasos para ejecutar
1. Clonar el repositorio del proyecto: 
`git clone https://github.com/pablobd03/ecuacion-schrodinger.git`

2. Acceder al directorio del repositorio:
`cd ecuacion-schrodinger`

3. Compilar:
`make`

4. Ejecutar:
`./schrodinger`

5. Para eliminar los archivo `.o`, y los archivos `.txt` generados por el programa, ejecutar:
`make clean`

## Aspectos a considerar
En caso de querer cambiar los parámetros respectivos al pozo finito unidimensional, estos se encuentran en el archivo `src/main c`. Considere, además, la función de cada parámetro:
- N: número de puntos de la malla
- size: tamaño de la matriz Hamiltoniano
- xmin: límite izquierdo
- xmax: límite derecho
- a: parámetro del potencial
- V0: potencial de profundidad
- k: índice del eigenvector a guardar (k>=0 y k<`N`)

Los datos generados en la ejecución del proyecto correspondientes a las funciones de onda y los niveles de energía se encuentran en el directorio `data/`.

