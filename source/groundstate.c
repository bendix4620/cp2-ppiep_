#include <stdio.h>
// #include <math.h>
#include <complex.h>
#include "cfg.h"
#include "conjugategradient.h"
#include "geometry.h"
#include "hamiltonian.h"
#include "powermethod.h"
#include "vectorutils.h"

#define CONFIG "../configs/groundstate"

// GLOBAL parameters
int N, D, ND;
double mu, eps2;

// Algorithm properties
int maxiter;
double precision;

// Local parameters
char filename[255];
FILE *fp;
linMap hamil = &hamiltonian;

/*
 * Apply the inverse hamiltonian onto psi via conjugate-gradient-method
 *
 * Parameters
 * ----------
 * double complex *out : pointer to complex array that the result is saved to
 * double complex *psi : pointer to complex array that the inverse hamiltonian
 *                       acts on
 * int len : length of the arrays
 */
void inverse_hamiltonian(double complex *out, double complex *psi, int len)
{
    inverse(out, psi, hamil, len, precision, maxiter);
    return;
}

void write_header(FILE *fp)
{
    if (fp == NULL)
    {
        printf("Could not open file %s", filename);
        exit(EXIT_FAILURE);
    }
    else
        printf("file %s opened\n", filename);

    fprintf(fp,
        "# Estimated Ground State Wavefunction and Energy of Time Invariant "
        "Hamiltonian\n#\n# Parameters:\n# N = %d\n# D = %d\n# mu = %f\n"
        "# eps^2 = %f\n# precision = %f\n# maxiter = %d\n# Results\n# "
        "---------------------------------------"
        "---------------------------------------\n",
        N, D, mu, eps2, precision, maxiter
    );
    return;
}

int main(void)
{
    // read config
    // -------------------------------------------------------------------------
    cfgAlloc(CONFIG);
    N = cfgGeti("N");
    D = cfgGeti("D");
    ND = ipow(N, D);
    mu = cfgGetd("mu");
    eps2 = cfgGetd("eps2");

    precision = cfgGetd("precision");
    maxiter = cfgGetd("maxiter");
    cfgFree();
    printf("N = %d, D = %d, mu = %f, eps^2 = %f\n", N, D, mu, eps2);

    // Prepare geometries
    prepare_geometries();
    use_quartic_pot();

    // Calculate ground state of hamiltonian with potential pot
    // -------------------------------------------------------------------------
    linMap invhamil = &inverse_hamiltonian;
    double complex evec[ND];
    double eval;
    solveeigen(evec, &eval, invhamil, ND, precision, maxiter);

    // Save to file
    // -------------------------------------------------------------------------
    sprintf(filename, "groundstate_N%d_D%d_mu%f_epssq%f", N, D, mu, eps2);
    fp = fopen(filename, "w");
    write_header(fp);
    fprintf(fp, "# Energy:\n%f\n# Wavefunction:\n", 1/eval);
    for (int i = 0; i < ND; i ++)
        fprintf(fp, "%f%+fj\n", creal(evec[i]), cimag(evec[i]) );
    fclose(fp);

    // free geometries
    free_geometries();
    free_potential();
    printf("Terminated successfully\n");
    return 0;
}
