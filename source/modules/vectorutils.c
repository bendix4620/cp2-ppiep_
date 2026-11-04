#include "vectorutils.h"

double complex vscalarprod(double complex *x, double complex *y, int len)
{
    double complex res = 0;
    for (int i = 0; i < len; i++)
    {
        // (a - ib) * (c + id) = ac + bd + i(ad - bc)
        res += conj(x[i]) * y[i];
    }
    return res;
}

double vscalarprod_re(double complex *x, double complex *y, int len)
{
    double res = 0;
    for (int i = 0; i < len; i++)
        res += creal(x[i]) * creal(y[i]) + cimag(x[i]) * cimag(y[i]);
    return res;
}

double vnorm2(double complex *x, int len)
{
    return vscalarprod_re(x, x, len);
}

double vnorm(double complex *x, int len)
{
    return sqrt(vnorm2(x, len));
}

void voverwrite(double complex *x, double complex *y, int len)
{
    for (int i = 0; i < len; i++)
        y[i] = x[i];
    return;
}

void vscale(double complex *x, double complex a, int len)
{
    for (int i = 0; i < len; i++)
        x[i] *= a;
    return;
}

void vfmadd(double complex *out, double a, double complex *x, double complex *y,
            int len)
{
    for (int i = 0; i < len; i++)
        out[i] = (
            fma( a, creal(x[i]), creal(y[i]) ) +
            fma( a, cimag(x[i]), cimag(y[i]) ) * I
        );
    return;
}


int vcompare(double complex *x, double complex *y, int len, double delta)
{
    for (int i = 0; i < len; i++)
        if ( cabs(x[i] - y[i]) > delta)
            return 1;
    return 0;
}
