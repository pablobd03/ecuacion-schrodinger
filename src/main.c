#include "matrix.h"
#include "potential.h"
#include "solucionar.h"
#include <stdio.h>
#include <gsl/gsl_matrix.h>
#include <math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>

int main() {
    /* 1. Parámetros del problema */
    int N = 20;                   // número de puntos de la malla
    size_t size = N;
    double xmin = -10.0;         // límite izquierdo
    double xmax =  10.0;         // límite derecho

    double a = 4.0;             // parámetro del potencial
    double V0 = 5.0;           // potencial de profundidad

    int k = 2; // índice del eigenvector a guardar (k>=0 y k<N)

    /* 2. Reserva de espacio en memoria */
    // Crear matrix del hamiltoniano
    gsl_matrix* H = gsl_matrix_alloc(size, size);
    // Crear vectores y matrices para valores y vectores propios
    gsl_vector* valores_propios = gsl_vector_alloc(size);
    gsl_matrix* vectores_propios = gsl_matrix_alloc(size, size);

    /* 3. Construcción del Hamiltoniano */
    build_hamiltonian(H, N, xmin, xmax, a, V0);
    
    /* 4. Impresión del Hamiltoniano */
    printf("Hamiltoniano:\n");
    print_matrix(H, size);


    /* 4. Cálculo de valores y vectores propios */
    solve_eigenpairs(H, N, valores_propios, vectores_propios, xmin, xmax);
    
    /* 5. Impresión de valores propios */
    printf("\nValores Propios:\n");
    print_vector(valores_propios, size);

    /* 6. Impresión de vectores propios */
    printf("\nVectores Propios:\n");
    print_matrix(vectores_propios, size);

    /* 7. Guardar eigenfunción k-ésima en archivo */
    int status = save_eigenfunction(vectores_propios, k, N, xmin, xmax, "eigenfunction.txt");
    if (status == 0) {
        printf("\nEigenfunción %d guardada en eigenfunction.txt\n", k);
    } else {
        printf("\nError al guardar la eigenfunción %d (código de error: %d)\n", k, status);
    }

    /* 8. Guardar eigenvalores */
    status = save_eigenvalues(valores_propios, size, "eigenvalues.txt");
    if (status == 0) {
        printf("\nEigenvalores guardados en 'eigenvalues.txt'\n");
    } else {
        printf("\nError al guardar los eigenvalores (código de error: %d)\n", status);
    }
    
    
    /* 9. Liberar memoria */
    gsl_matrix_free(H);
    gsl_vector_free(valores_propios);
    gsl_matrix_free(vectores_propios);


    return 0;
}