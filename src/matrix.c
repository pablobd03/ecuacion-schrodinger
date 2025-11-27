#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>


//Funciones de impresión de matrices y vectores GSL, para depuración, comprobación de datos y pruebas                                                    
void print_matrix(gsl_matrix* matrix, size_t size) {
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            printf("%g ", gsl_matrix_get(matrix, i, j));
        }
        printf("\n");
    }              
    printf("\n");
}

void print_vector(gsl_vector* vector, size_t size) {
    for (size_t i = 0; i < size; i++) {
        printf("%g\n", gsl_vector_get(vector, i));
    }
    printf("\n");
}
