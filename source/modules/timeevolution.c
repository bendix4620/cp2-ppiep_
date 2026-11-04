#include "timeevolution.h"

void euler(double complex *out, double complex *psi, double tau, int len)
{
    hamiltonian(out, psi, len);
    for (int i = 0; i < len; i++)
    {
        out[i] = psi[i] - I * tau * out[i];
    }
    return;
}



double _tau = 0;

/* calculate out = ( 1 + 1/4 tau^2 H^2 ) psi */
void oneaddsqH(double complex *out, double complex *psi, int len)
{
    double complex temp[len];
    // out = H^2 psi0
    hamiltonian(temp, psi, len);
    hamiltonian(out, temp, len);
    vfmadd(out, _tau*_tau/4, out, psi, len);
}

/* calculate out = ( 1 - i/2 tau H ) psi */
void onesubiH(double complex *out, double complex *psi, int len)
{
    hamiltonian(out, psi, len);
    for (int i = 0; i < len; i++)
        out[i] = psi[i] - I/2 * _tau * out[i];
    return;
}

void cranknicolson(double complex *out, double complex *psi, double tau,
    int len, double eps, int maxiter)
{
    double complex temp[len];
    _tau = tau; // write tau to global variable so sqhamilton can read it

    // out = (1 + tau^2/4 H^2)^-1
    linMap func = &oneaddsqH;
    inverse(out, psi, func, len, eps, maxiter);

    // temp = H out
    onesubiH(temp, out, len);
    onesubiH(out, temp, len);
    return;
}


int issetup = 0; // 0 if false, 1 if true
fftw_complex *fftin, *fftout;
fftw_plan pf, pb; // plan forward, plan backward

int setup_strangsplitting(int N, int D)
{
    if (issetup == 1)
        return -1;

    int len = ipow(N, D);
    int shape[D];
    for (int i = 0; i < D; i++)
        shape[i] = N;

    fftin = (fftw_complex*)fftw_malloc(len * sizeof(fftw_complex));
    fftout = (fftw_complex*)fftw_malloc(len * sizeof(fftw_complex));

    pf = fftw_plan_dft(D, shape, fftin, fftout, FFTW_FORWARD, FFTW_MEASURE);
    pb = fftw_plan_dft(D, shape, fftin, fftout, FFTW_BACKWARD, FFTW_MEASURE);

    issetup = 1;
    return 0;
}

int free_strangsplitting()
{
    if (issetup == 0)
        return -1;

    fftw_destroy_plan(pf);
    fftw_destroy_plan(pb);
    fftw_free(fftin);
    fftw_free(fftout);

    issetup = 0;
    return 0;
}

void strangsplitting(double complex *out, double complex *psi, double tau,
    double mu, double epssq, int N, int D)
{
    int len = ipow(N, D);
    double temp;
    double complex prefactor, ek, ev[len];

    // eta =  exp(-i tau/2 V) psi
    prefactor = -I * tau / 2;
    for (int i = 0; i < len; i++)
    {
        ev[i] = cexp(prefactor * potential[i]); // "e^V" operator
        fftin[i] = ev[i] * psi[i];
    }

    // ~eta = sum[eta e*]
    fftw_execute(pf);

    // ~xi = exp(-i tau K) ~eta
    prefactor = -I *tau * 2 / (mu * epssq);
    for (int i = 0; i < len; i++)
    {
        ek = 0; // "e^K" operator
        for (int j = 0; j < D; j++)
        {
            // eigenvalues K(k)
            temp = sin(M_PI * (coords[i][j]) / N ); // coords calculated from lexicographic index
            ek += temp * temp;
        }

        fftin[i] = cexp(prefactor * ek) * fftout[i] / len; // normation
    }

    // xi = sum[~xi e]
    fftw_execute(pb);

    // psi(q+1) = exp(-i tau/2 V) xi
    for (int i = 0; i < len; i++)
        out[i] = ev[i] * (double complex)fftout[i];
    return;
}
