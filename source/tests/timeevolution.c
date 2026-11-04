#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h> // needs to be included AFTER complex.h
#include "hamiltonian.h"
#include "vectorutils.h"
#include "geometry.h"
#include "cfg.h"
#include "timeevolution.h"
#include "observables.h"



#define CONFIG "../../configs_tests/timeevolution"
int N, D, ND;
double mu, eps2, delta;


/*
 * Test unitarity of the Crank-Nicolson integrator
 * by comparing Crank-Nicolson for one time step forward and backwards
 *
 * Parameters
 * ----------
 * double complex *phi: pointer to any wavefunction phi
 * int len: length of psi/phi array (N^D)
 * double delta: comparison tolerance value (f.e. 1e-10)
 * int maxiter: number of iterations of crank nicolsen
 */
int unitarityCN(double complex *phi, int len, double tau, double delta_inv, int maxiter)
{
    double complex temp[len];
    double diff, norm = vnorm(phi, len);
    cranknicolson(temp, phi, tau, len, delta_inv, maxiter);
    diff = fabs(norm - vnorm(temp, len));

    if (diff > delta)
    {
        printf("unitarity test failed. delta = %e\n", diff);
        return -1;
    }
    return 0;
}


/*
 * Test energy conversation of the Crank-Nicolson integrator:
 *
 * Parameters
 * ----------
 * double complex *psi: pointer to any wavefunction psi
 * int len: length of psi array (N^D)
 * double delta: comparison tolerance value (f.e. 1e-10)
 * int maxiter: number of iterations of crank nicolsen
 */
int energyconversationCN(double complex *psi, int len, double tau, double delta_inv, int maxiter)
{
    double complex temp[len];
    double energy0, energy1, diff, norm2;

    norm2 = vnorm2(psi, len);
    energy0 = obs_e(psi, norm2, len);
    cranknicolson(temp, psi, tau, len, delta_inv, maxiter);
    energy1 = obs_e(temp, norm2, len);

    diff = fabs(energy0 - energy1);
    if ( diff < delta )
        return 0;
    else
        printf("energy conversion test failed:\ndelta = %e\n", diff);

    return -1;
}



/*
 * Test unitarity of the Strang Splitting integrator:
 * by comparing for one time step forward and backwards
 * Test unitarity of the Strang Splitting integrator
 *
 * Parameters
 * ----------
 * double complex *psi: pointer to any wavefunction psi
 * int len: length of psi/phi array (N^D)
 * double delta: comparison tolerance value (f.e. 1e-10)
 */
int unitaritySS(double complex *psi, int len, double tau)
{
    double complex temp[len];
    double norm, diff;
    norm = vnorm(psi, len);
    strangsplitting(temp, psi, tau, mu, eps2, N, D);
    diff = fabs(norm - vnorm(temp, len));

    if (diff > delta)
    {
        printf("unitary test failed.\ndelta = %e", diff);
        return -1;
    }
    return 0;
}


/* Return complex value from inside the unit circle */
double complex runit()
{
    return (rand() + I*rand() ) / RAND_MAX * 2 - (1 + I);
}


/* fill random complex array */
void rcarr(double complex *arr, int len)
{
    for (int i = 0; i < len; i++)
        arr[i] = runit();
    vscale(arr, 1/vnorm(arr, len), len);
    return;
}


int main(void)
{
    cfgAlloc(CONFIG);
    N = cfgGeti("N");
    D = cfgGeti("D");
    mu = cfgGetd("mu");
    eps2 = cfgGetd("eps2");
    delta = cfgGetd("delta");
    double delta_inv = cfgGetd("precision");
    int reps = cfgGeti("reps");
    double tau = cfgGetd("tau");
    ND = ipow(N, D);
    int maxiter =  cfgGeti("maxiter");
    cfgFree();

    printf("N = %d, D = %d, mu = %f, eps^2 = %f, precision = %e, maxiter = %i\n",
        N, D, mu, eps2, delta, maxiter);
    prepare_geometries();
    use_quartic_pot();
    setup_strangsplitting(N, D);
    srand(time(NULL));

    printf("\nCrank-Nicolson-algorithm tests:\n-------------------------------\n");
    double complex phi[ND];
    for (int i = 0; i < reps; i++)
    {
        rcarr(phi, ND);
        if ( 0 != unitarityCN(phi, ND, tau, delta_inv, maxiter) )
            break;
    }
    printf("unitary test successful.\n");

    for (int i = 0; i < reps; i++)
    {
        rcarr(phi, ND);
        if ( 0 != energyconversationCN(phi, ND, tau, delta_inv, maxiter) )
            break;
    }
    printf("energy conversion test successful.\n");

    printf("\nStrang-Splitting-Algorithm:\n---------------------------\n");
    for (int i = 0; i < reps; i++)
    {
        rcarr(phi, ND);
        if ( 0 != unitaritySS(phi, ND, tau) )
            break;
    }
    printf("unitary test successful.\n");
    free_geometries();
    free_potential();
    free_strangsplitting();
    return 0;
}
