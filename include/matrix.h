#ifndef MATRIX_H
#define MATRIX_H

#include <gsl/gsl_matrix.h>

//Funcion para crear el Hamiltoniano para probar valores propios
void print_matrix(gsl_matrix* matriz, size_t size);

void print_vector(gsl_vector* matriz, size_t size);
#endif 