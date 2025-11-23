#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>

#include "solver.h"
#include "potential.h"

/*
 * Construcción del Hamiltoniano discreto para la ecuación de Schrödinger:
 *    H = -d²/dx² + V(x)
 * usando un método de diferencias finitas:
 *
 *  H[i][i]     = 2/dx² + V(x_i)
 *  H[i][i±1]   = -1/dx²
 *
 * Se utilizan unidades naturales: hbar²/(2m) = 1
 */
void build_hamiltonian(double **H, int N, double xmin, double xmax,
                       double a, double V0)
{
    double dx = (xmax - xmin) / (N - 1);
    double inv_dx2 = 1.0 / (dx * dx);

    for (int i = 0; i < N; i++) {
        double x = xmin + i * dx;

        for (int j = 0; j < N; j++)
            H[i][j] = 0.0;   // Inicializar fila

        // Término diagonal (cinético + potencial)
        H[i][i] = 2.0 * inv_dx2 + V(x, a, V0);

        // Fuera de la diagonal (término cinético)
        if (i > 0)
            H[i][i - 1] = -inv_dx2;

        if (i < N - 1)
            H[i][i + 1] = -inv_dx2;
    }
}



/*
 * Resolver autovalores del Hamiltoniano con GSL.
 * Imprime los 7  primeros autovalores (los más bajos).
 */
void solve_eigenvalues(double **H, int N)
{
    gsl_matrix *A = gsl_matrix_alloc(N, N);

    // Copiar matriz dinámica H[][] a matriz GSL
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            gsl_matrix_set(A, i, j, H[i][j]);

    gsl_vector *eval = gsl_vector_alloc(N);
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(N);

    // Resolver autovalores
    gsl_eigen_symm(A, eval, w);

    printf("=====================================\n");
    printf("  Autovalores (niveles de energia)\n");
    printf("=====================================\n");

    for (int i = 0; i < 7; i++)
        printf("E%d = %lf\n", i, gsl_vector_get(eval, i));

    printf("=====================================\n");

    // Liberar memoria GSL
    gsl_eigen_symm_free(w);
    gsl_vector_free(eval);
    gsl_matrix_free(A);
}
