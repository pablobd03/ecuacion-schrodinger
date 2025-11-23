#ifndef solucionar_h
#define solucionar_h

#include <gsl/gsl_errno.h>
#include <gsl/gsl_eigen.h>

/*
 * Construye la matriz Hamiltoniana H (NxN)
 * usando diferencias finitas.
 */

void build_hamiltonian(double **H, int N, double xmin, double xmax,
                       double a, double V0);

/*
 * Resuelve los autovalores del Hamiltoniano usando la biblioteca GSL.
 */


void solve_eigenvalues(double **H, int N);

#endif

