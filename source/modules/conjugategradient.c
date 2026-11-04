#include "conjugategradient.h"

void inverse(double complex *x, double complex *b, linMap func, int len,
             double eps, int maxiter)
{
    // initial breaking check
    if ( vnorm(b, len) < eps )
        return;

    complex r[len], p[len], fp[len];
    double alpha, beta, rnorm2, r1norm2;

    // initial values
    for (int i = 0; i < len; i++)
    {
        x[i] = 0; // chose x_0 randomly
        r[i] = b[i]; // r_0 = b - f(x_0 = 0) = b - 0 = b
        p[i] = b[i]; // p_0 = r_0
    }

    r1norm2 = vnorm2(r, len);
    for (int k = 0; k < maxiter; k++)
    {
        // alpha = ||r_k|| / <p_k, f(p_k)>
        rnorm2 = r1norm2; // = ||r_k||^2
        func(fp, p, len);
        alpha = rnorm2 / vscalarprod_re(p, fp, len); // hermitian func => real

        // r_k+1 = r_k - alpha_k * f(p_k)
        vfmadd(r, -alpha, fp, r, len);

        // x_k+1 = x_k + alpha_k * p_k
        vfmadd(x, alpha, p, x, len);

        // breaking criteria
        r1norm2 = vnorm2(r, len);
        if ( r1norm2 < eps * eps )
            // x_k+1 is already written to x
            return;

        // beta_k = ||r_k+1|| / ||r_k||
        beta = r1norm2 / rnorm2;

        // p_k+1 = r_k+1 + beta_k * p_k
        vfmadd(p, beta, p, r, len);
    }
    printf("Maximum number of iterations reached, conjugate-gradient algorithm failed.\n");
    exit(EXIT_FAILURE);
}
