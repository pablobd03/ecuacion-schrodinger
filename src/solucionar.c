#include <stdio.h>
#include <math.h>
#include <gsl/gsl_eigen.h>

#include "solucionar.h"
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
void solve_eigenpairs(double **H, int N, double *valores, double **vectores,
                      double xmin, double xmax)
{
    int i, j, k;
    double dx = (xmax - xmin) / (N - 1);

    /* Crear matriz GSL y copiar H */
    gsl_matrix *A = gsl_matrix_alloc(N, N);
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            gsl_matrix_set(A, i, j, H[i][j]);

    /* Preparar GSL para obtener eigenvalores y eigenvectores */
    gsl_vector *eval = gsl_vector_alloc(N);
    gsl_matrix *evec = gsl_matrix_alloc(N, N);
    gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(N);

    /* Cálculo */
    gsl_eigen_symmv(A, eval, evec, w);
    gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_VAL_ASC);

    /* Copiar autovalores */
    for (i = 0; i < N; i++)
        valores[i] = gsl_vector_get(eval, i);

    /* Copiar autovectores en formato matricial: guardar por filas: vectores[k][i] = componente i del autovector k */
    for (k = 0; k < N; k++) {
        /* fila k corresponde al k-ésimo autovector */
        for (i = 0; i < N; i++) {
            /* gsl_matrix_get(evec, i, k) es componente i del k-ésimo vector (columna k) */
            vectores[k][i] = gsl_matrix_get(evec, i, k);
        }

        /* Normalizar con respecto a la integral discreta: sum_i |psi_i|^2 * dx = 1 */
        double s = 0.0;
        for (i = 0; i < N; i++)
            s += vectores[k][i] * vectores[k][i];

        s = sqrt(s * dx); /* norma L2 en espacio continuo aproximada por sqrt(sum |psi|^2 * dx) */

        if (s == 0.0) continue; /* evitar división por cero */
        for (i = 0; i < N; i++)
            vectores[k][i] /= s;
    }

    /* Imprimir primeros 5 autovalores para verificar*/
    printf("=====================================\n");
    printf("  Autovalores (niveles de energia)\n");
    printf("=====================================\n");
    for (i = 0; i < N && i < 5; i++)
        printf("E%d = % .12f\n", i, valores[i]);
    printf("=====================================\n");

    /* Liberar estructuras GSL */
    gsl_eigen_symmv_free(w);
    gsl_matrix_free(evec);
    gsl_vector_free(eval);
    gsl_matrix_free(A);
}

/* Guardar la eigenfunción k en archivo con columnas: x \t psi */
int save_eigenfunction(double **vectores, int k, int N,
                       double xmin, double xmax, const char *filename)
{
    if (k < 0) return -1;
    if (k >= N) return -2;

    double dx = (xmax - xmin) / (N - 1);

    FILE *f = fopen(filename, "w");
    if (!f) return -3;

    for (int i = 0; i < N; i++) {
        double x = xmin + i * dx;
        fprintf(f, "% .12f\t% .12f\n", x, vectores[k][i]);
    }

    fclose(f);
    return 0;
}


/* Guardar los eigenvalores es decir las energías*/

int save_eigenvalues(double *valores, int N, const char *filename)
{
    FILE *f = fopen(filename, "w");
    if (!f) return -1;

    for (int k = 0; k < N; k++)
        fprintf(f, "%d\t%.12f\n", k, valores[k]);

    fclose(f);
    return 0;
}

