#ifndef solucionar_h
#define solucionar_h

#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/*
 * Construye la matriz Hamiltoniana H (NxN)
 * usando diferencias finitas en unidades naturales:
 *  -H: matriz hamiltoniano (NxN)
 *  -N: tamaño de la malla
 *  -xmin, xmax: dominio espacial
 *  -a: ancho potencial
 *  -V0: Valor de la grada de potencial
 */

void build_hamiltonian(gsl_matrix *H, int N, double xmin, double xmax,
                       double a, double V0);


/* Resolver autovalores y autovectores:
 *  - H: matriz Hamiltoniano (NxN)
 *  - N: tamaño de la malla
 *  - valores: arreglo ya reservado de tamaño N para guardar autovalores ordenados asc.
 *  - vectores: matriz ya creada (N x N). Se guarda como:
 *        vectores[k][i] = componente i del autovector k
 *    (es decir, cada fila k contiene la eigenfunción k)
 *  - xmin, xmax: dominio espacial 
 */
void solve_eigenpairs(gsl_matrix *H, int N, gsl_vector *valores, gsl_matrix *vectores,
                      double xmin, double xmax);

/* Guardar autovalores en archivo:
 *  - valores: arreglo con N autovalores
 *  - filename: nombre de archivo a crear
 *  Devuelve 0 = OK, -1 = error abriendo archivo.
 */
int save_eigenvalues(gsl_vector *valores, int N, const char *filename);

/* Guardar la eigenfunción k (fila vectores[k]) en archivo con columnas: x \t psi
 *  - vectores: matriz con eigenfunciones por fila
 *  - k: índice del estado a guardar (0..N-1)
 *  - xmin, xmax: dominio espacial
 *  - filename: nombre del archivo
 *  Devuelve 0 = OK, -1 = k fuera de rango, -2 = error abriendo archivo
 */
int save_eigenfunction(gsl_matrix *vectores, int k, int N,
                       double xmin, double xmax, const char *filename);

#endif
