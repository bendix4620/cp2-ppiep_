#include "powermethod.h"

/* Return squared norm of an eigenvalue approximation ||f(v) - k*v||^2
 * where fv is the mapped eigenvector vector f(v), v is the
 * eigenvector and k is the eigenvalue.
 * Will return 0 if v and k are exact solutions */
double approxEigenProblem(complex *fv, complex *v, double k, int len)
{
    double res = 0;
    for (int i = 0; i < len; i ++)
        res += cabs(fv[i] - k*v[i]);
    // printf("res = %f\n", res);
    return res;
}

void solveeigen(double complex *evec, double *eval, linMap func, int len,
                double eps, int maxiter)
{
    double complex w[len];

    // chose v (out) randomly
    for (int i = 0; i < len; i++)
        evec[i] = 1; // must be != 0


    for (int i = 0; i < maxiter; i++)
    {
        // calculate w = Av
        func(w, evec, len);

        // get eigenvalue
        *eval = vnorm(w, len);

        // breaking criteria
        if (approxEigenProblem(w, evec, *eval, len) < eps)
            return;

        // replace v = w/k
        voverwrite(w, evec, len);
        vscale(evec, 1/(*eval), len);
    }
    printf("Maximum number of iterations reached, power-method algorithm failed.\n");
    exit(EXIT_FAILURE);
}
