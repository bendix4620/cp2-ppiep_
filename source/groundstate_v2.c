#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include "cfg.h"
#include "conjugategradient.h"
#include "geometry.h"
#include "hamiltonian.h"
#include "powermethod.h"
#include "vectorutils.h"
#include "observables.h"

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
double norm_psi, norm_psi2, kinetic_energy, potential_energy, total_energy;

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

/* write header of output file */
void make_header()
{
    fprintf(fp,
        "# time evolution of an initial gaussian wavelet\n"
        "# \n"
        "# Parameters\n"
        "# ----------\n"
        "# N = %d\n"
        "# D = %d\n"
        "# \n"
        "# mu = %e\n"
        "# eps2 = %e\n"
        "# \n"
        "# maxiter = %d\n"
        "# precision = %e\n"
        "# \n"
        "# Columns contents\n"
        "# ----------------\n"
        "# norm, energy, kinetic energy, potential energy, position(D), squared position(D), wavefunction(N^D)\n",
        N, D, mu, eps2, maxiter, precision);
    return;
}

/* write a content line to the output file */
void write_line(double complex *psi, double *vec, double *unc)
{
    // norm and scalar observables
    norm_psi2 = vnorm2(psi, ND);
    norm_psi = sqrt(norm_psi2);

    kinetic_energy = obs_ekin(psi, norm_psi2, ND);
    potential_energy = obs_epot(psi, norm_psi2, ND);
    total_energy = obs_e(psi, norm_psi2, ND);

    fprintf(fp, "%e %e %e %e",
        norm_psi, total_energy, kinetic_energy, potential_energy
    );

    // position and uncertainty
    obs_pos(vec, D, psi, norm_psi2, ND, 1);
    obs_pos(unc, D, psi, norm_psi2, ND, 2);

    for (int i = 0; i < D; i++)
        fprintf(fp, " %e", vec[i]);
    for (int i = 0; i < D; i++)
        fprintf(fp, " %e", unc[i]);

    // wavefunction
    for (int i = 0; i < ND; i++)
        fprintf(fp, " %e%+ej", creal(psi[i]), cimag(psi[i]) );
    fprintf(fp, "\n");
    return;
}

/* Search for key in command line arguments and parse its value to int *res.
 * Return 0 on success, -1 on failure */
int iparse(int *res, char *keyname, char *argv[], int argc)
{
    char key[255];
    int val[1];
    for (int i = 0; i < argc; i++)
    {
        if (sscanf(argv[i], "%32[^=]=%d", key, val) == 2)
            if (strcmp(keyname, key) == 0)
            {
                *res = val[0];
                return 0;
            }
    }
    return -1;
}

/* Search for key in command line arguments and parse its value to double *res.
 * Return 0 on success, -1 on failure */
int dparse(double *res, char *keyname, char *argv[], int argc)
{
    char key[255];
    double val[1];
    for (int i = 0; i < argc; i++)
    {
        if (sscanf(argv[i], "%32[^=]=%lf", key, val) == 2)
            if (strcmp(keyname, key) == 0)
            {
                *res = val[0];
                return 0;
            }
    }
    return -1;
}

int main(int argc, char *argv[])
{
    // read config
    // -------------------------------------------------------------------------
    cfgAlloc(CONFIG);
    N = cfgGeti("N");
    D = cfgGeti("D");
    mu = cfgGetd("mu");
    eps2 = cfgGetd("eps2");

    precision = cfgGetd("precision");
    maxiter = cfgGetd("maxiter");
    cfgFree();

    // read command line to potentially overwrite values
    // -------------------------------------------------
    iparse(&N, "N", argv, argc);
    dparse(&mu, "mu", argv, argc);
    dparse(&eps2, "eps2", argv, argc);

    ND = ipow(N, D);
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
    double vec[ND], unc[ND];
    sprintf(filename, "groundstate_v2_N%d_D%d_mu%f_epssq%f", N, D, mu, eps2);
    fp = fopen(filename, "w");
    make_header();
    write_line(evec, vec, unc);

    // free geometries
    free_geometries();
    free_potential();
    fclose(fp);
    printf("Terminated successfully\n");
    return 0;
}
