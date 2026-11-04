#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include "hamiltonian.h"
#include "vectorutils.h"
#include "geometry.h"
#include "cfg.h"

#define CONFIG "../../configs_tests/hamiltonian"
int N, D, ND;
double mu, eps2;

/*
 * Test linearity of hamiltonian func:
 * H(αψ + βϕ) = α(Hψ) + β(Hϕ)
 *
 * Parameters
 * ----------
 * double complex *psi: pointer to any wavefunction psi
 * double complex *phi: pointer to any wavefunction phi
 * double complex alpha: complex number
 * double complex beta: complex number
 * int len: length of psi/phi arrwy (N^D)
 * double delta: comparison tolerance value (f.e. e-6)
 */
int linearity(double complex *psi, double complex *phi, double complex alpha,
               double complex beta, int len, double delta)
{
    double complex h_psi[len], h_phi[len], temp[len];
    hamiltonian(h_psi, psi, len); // H(psi)
    hamiltonian(h_phi, phi, len); // H(phi)

    for (int i = 0; i < len; i++)
    {
        temp[i] = (alpha * psi[i]) + (beta * phi[i]); // a*psi + b*phi
        h_psi[i] = (alpha * h_psi[i]) + (beta * h_phi[i]); // a*H(psi) + b*H(phi)
    }
    hamiltonian(h_phi, temp, len); // H(a*psi + b*phi)

    if (vcompare(h_phi, h_psi, len, delta) != 0)
    {
        printf("Linearity test failed\n");
        return 1;
    }
    return 0;
}

/*
 * Test hermiticity of hamiltonian func:
 * (ψ,Hϕ) = (Hψ,ϕ)
 *
 * Parameters
 * ----------
 * double complex *psi: pointer to any wavefunction psi
 * double complex *phi: pointer to any wavefunction phi
 * int len: length of psi/phi arrwy (N^D)
 * double delta: comparison tolerance value (f.e. e-6)
 */
int hermiticity(double complex *psi, double complex *phi, int len, double delta)
{
    double complex h_psi[len], h_phi[len], x;
    hamiltonian(h_psi, psi, len);
    hamiltonian(h_phi, phi,  len);

    x  = vscalarprod(h_psi, phi, len) - vscalarprod(psi, h_phi, len);
    if (cabs(x) > delta)
    {
        printf("Hermiticity test failed with r = %e\n", cabs(x));
        return 1;
    }
    return 0;
}

/*
 * Test positivity of hamiltonian func:
 * For V >= 0:
 * (ψ,Hψ) >= 0
 *
 * Parameters
 * ----------
 * double complex *psi: pointer to any wavefunction psi
 * int len: length of psi/phi arrwy (N^D)
 */
int positivity(double complex *psi, int len)
{
    use_quartic_pot(); /* V >= 0 */

    double complex h_psi[len], x;
    hamiltonian(h_psi, psi, len);

    x = vscalarprod(psi, h_psi, len);
    if (cabs(x) < 0)
    {
        printf("Positivity test failed with r = %e\n", cabs(x));
        return 1;
    }
    return 0;
}


/*
 * set plane wave function
 * ψ_k(n) = e^(2π*i*n*k/N) k ∈ Z^D
 *
 * Parameters
 * ----------
 * double complex *psi: pointer to any wavefunction psi
 * int k: index of wavefunction
 *
 */
void set_planewavefunc(double complex *psi, int *k, int len)
{
    double nk;
    for (int i = 0; i < len; i++)
    {
        // nk = <n, k>
        nk = 0;
        for (int j = 0; j < D; j++)
            nk += coords[i][j] * k[j];
        // psi = exp( i 2pi/N <n, k> )
        psi[i] = cexp( 2* M_PI * I / N * nk);
    }
    return;
}

/*
 * Test eigenvalues and eigenvectors of hamiltonian func for two differen wavefunctions:
 * For V = 0:
 * plane waves ψ_k(n) = e^(2π*i*n*k/N) k ∈ Z^D are eigenfunctions of H
 *
 * to check:
 * H ψ_k = E_k ψ_k
 *
 * Parameters
 * ----------
 * double complex *psi1,2: pointer to any wavefunction psi
 * int len: length of psi/phi array (N^D)
 * int *k : pointer to array describing the wave vector k
 * double delta: comparison tolerance value (f.e. e-6)
 *
 */
int eigenvalue_eq(double complex *psi, int len, int *k, double delta)
{
    double complex h_psi[len];

    double ev = 0, prefactor = -1 / (mu * eps2);
    use_zero_pot(); /* V = 0 */

    // calculate eigenvalue (energy)
    for (int i = 0; i < D; i++)
        ev += cos( 2 * M_PI / N * k[i]) - 1;
    ev *= prefactor;

    // calculate H(psi) with psi = eigenvector of H
    set_planewavefunc(psi, k, len);
    // test_wavefunc_periodic(psi, len);
    hamiltonian(h_psi, psi, len);
    // calculate E*psi
    vscale(psi, ev, len);

    // compare H(psi) == E*psi
    if (vcompare(h_psi, psi, len, delta) != 0)
    {
        printf("Eigenvalue test failed\n");
        return 1;
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
    return;
}

/* fill random integer array */
void riarr(int *arr, int len, int max)
{
    for (int i = 0; i < len; i++)
        arr[i] = rand() % max;
    return;
}


int main(){
    // read config
    // -------------------------------------------------------------------------
    cfgAlloc(CONFIG);
    N = cfgGeti("N");
    D = cfgGeti("D");
    mu = cfgGetd("mu");
    eps2 = cfgGetd("eps2");
    int reps = cfgGeti("reps");
    double delta = cfgGetd("precision");
    ND = ipow(N, D);

    cfgFree();
    printf("N = %d, D = %d, mu = %f, eps^2 = %f, precision = %e\n",
        N, D, mu, eps2, delta);

    prepare_geometries();
    use_quartic_pot();
    srand(time(NULL));

    // Tests
    // -------------------------------------------------------------------------
    // set psi and phi wavefunctions
    double complex psi[ND], phi[ND], alpha, beta;
    int res, k[D]; // wave vector

    //linearity test
    for (int i = 0; i < reps; i++)
    {
        alpha = runit();
        beta = runit();
        rcarr(phi, ND);
        rcarr(psi, ND);
        res = linearity(psi, phi, alpha, beta, ND, delta);
        if (res != 0)
            break;
    }
    if (res == 0)
        printf("Linearity test successful. \n");
    res = 0;

    //hermiticity test
    for (int i = 0; i < reps; i++)
    {
        rcarr(phi, ND);
        rcarr(psi, ND);
        res = hermiticity(psi, phi, ND, delta);
        if (res != 0)
            break;
    }
    if (res == 0)
        printf("Hermiticity test successful. \n");
    res = 0;


    //positivity test
    for (int i = 0; i < reps; i++)
    {
        rcarr(psi, ND);
        res = positivity(psi, ND);
        if (res != 0)
            break;
    }
    if (res == 0)
        printf("Positivity test successful. \n");
    res = 0;

    //eigenvalue eq test
    for (int i = 0; i < reps; i++)
    {
        riarr(k, D, N);
        res = eigenvalue_eq(psi, ND, k, delta);
        if (res != 0)
            break;
    }
    if (res == 0)
        printf("Eigenvalue test successful. \n");


    // free
    free_geometries();
    free_potential();
    return 0;
}
