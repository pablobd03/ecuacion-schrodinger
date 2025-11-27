#include "matrix.h"
#include "potential.h"
#include "solucionar.h"
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

int main() {
    /* 1. Parámetros del problema 
     * N: número de puntos de la malla
     * size: tamaño de la matriz Hamiltoniano
     * xmin: límite izquierdo
     * xmax: límite derecho
     * a: parámetro del potencial
     * V0: potencial de profundidad
     * k: índice del eigenvector a guardar (k>=0 y k<N) 
     */
    int N = 500;                   
    size_t size = N;              
    double xmin = -10.0;          
    double xmax =  10.0;          
    double a = 4.0;               
    double V0 = 5.0;              
    int k = 2;                    

    /* 2. Reserva de espacio en memoria 
     * H: Hamiltoniano
     * valores_propios: vector de valores propios (niveles de energía)
     * vectores_propios: matriz de vectores propios (funciones de onda)  
     */
    gsl_matrix* H = gsl_matrix_alloc(size, size);                 
    gsl_vector* valores_propios = gsl_vector_alloc(size);         
    gsl_matrix* vectores_propios = gsl_matrix_alloc(size, size);  
    /* 3. Construcción del Hamiltoniano */
    build_hamiltonian(H, N, xmin, xmax, a, V0);
    
    /* 4. Impresión del Hamiltoniano */
    printf("Hamiltoniano:\n");
    print_matrix(H, size);

    /* 5. Cálculo de valores y vectores propios */
    solve_eigenpairs(H, N, valores_propios, vectores_propios, xmin, xmax);
    
    /* 6. Impresión de valores propios */
    printf("\nValores Propios:\n");
    print_vector(valores_propios, size);

    /* 7. Impresión de vectores propios */
    printf("\nVectores Propios:\n");
    print_matrix(vectores_propios, size);

    /* 8. Guardado de funciones de onda (vectores propios) k-ésima en archivo */
    // Crear carpeta "data" si no existe
    int dir_status = mkdir("data", 0755);     
    if (dir_status == -1 && errno != EEXIST) {
        fprintf(stderr, "Ocurrió un error al crear la carpeta 'data'");
    }

    int status = save_eigenfunction(vectores_propios, k, N, xmin, xmax, "data/funciones_de_onda.txt");
    if (status == 0) {
        printf("\nFunciones de onda #%d guardada en 'data/funciones_de_onda.txt'\n", k);
    } else {
        printf("\nError al guardar la funciones de onda %d (código de error: %d)\n", k, status);
    }

    /* 9. Guardado de niveles de energía (valores propios) propios */
    status = save_eigenvalues(valores_propios, size, "data/niveles_de_energia.txt");
    if (status == 0) {
        printf("\nNiveles de energía guardados en 'data/niveles_de_energia.txt'\n");
    } else {
        printf("\nError al guardar los niveles de energía (código de error: %d)\n", status);
    }
    
    /* 10. Liberación de memoria */
    gsl_matrix_free(H);
    gsl_vector_free(valores_propios);
    gsl_matrix_free(vectores_propios);

    return 0;
}
