#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <complex.h>
#include "cfg.h"
#include "conjugategradient.h"

#define CONFIG "../../configs_tests/conjugategradient"
double complex **matrix_values;

/* fill matrix randomly with hermitian and positive semi-definite properties */
void prepare_matrix(int len)
{
    // allocate matrix
    matrix_values = (double complex**)malloc(len*sizeof(double complex*));
    for (int i = 0; i < len; i++)
        matrix_values[i] = (double complex*)malloc(len*sizeof(double complex));

    // hermitian matrix (except the diagonal)
    srand(time(NULL));
    for (int i = 0; i < len; i++)
        for (int j = i+1; j < len; j++)
        {
            matrix_values[i][j] = (rand() + I*rand())/RAND_MAX * 2 - (1+I);
            matrix_values[j][i] = conj(matrix_values[i][j]);
        }

    // positive semi-definite via dominant diagonal
    for (int i = 0; i < len; i++)
    {
        matrix_values[i][i] = 0;
        for (int j = 0; j < i; j++)
            matrix_values[i][i] += cabs(matrix_values[i][j]);
        for (int j = i+1; j < len; j++)
            matrix_values[i][i] += cabs(matrix_values[i][j]);

    }
    return;
}

/* free the matrix */
void free_matrix(int len)
{
    for (int i = 0; i < len; i++)
        free(matrix_values[i]);
    free(matrix_values);
    return;
}

/* the matrix to calculate the eigenvalue from */
void matrix(double complex *out, double complex *in, int len)
{
    for (int i = 0; i < len; i++)
    {
        out[i] = 0;
        for (int j = 0; j < len; j++)
            out[i] += matrix_values[i][j] * in[j];
    }
    return;
}

int main(void)
{
    cfgAlloc(CONFIG);
    int len = cfgGeti("len");
    double eps = cfgGetd("precision");
    int maxiter = cfgGeti("maxiter");
    cfgFree();
    printf("len = %d, precision = %e, maxiter = %d\n", len, eps, maxiter);

    double complex x[len], b[len], test[len];
    // random complex b vector
    for (int i = 0; i < len; i++)
        b[i] = (rand() + I*rand())/RAND_MAX * 2 - (1+I);

    prepare_matrix(len);
    linMap func = &matrix;
    inverse(x, b, func, len, eps, maxiter);

    func(test, x, len);
    free_matrix(len);
    for (int i = 0; i < len; i++)
        // printf("%.2f%+.2fi <-> %.2f%+.2fi\n", creal(test[i]), cimag(test[i]), creal(b[i]), cimag(b[i]));
        if ( cabs(test[i] - b[i]) > eps)
        {
            printf("Test of conjugated-gradient-method failed!\n");
            return 1;
        }
    printf("Test of conjugated-gradient-method successful!\n");
    return 0;
}
