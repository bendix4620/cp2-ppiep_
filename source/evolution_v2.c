#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <fftw3.h>
#include "cfg.h"
#include "vectorutils.h"
#include "timeevolution.h"
#include "observables.h"

#define CONFIG "../configs/evolution"

#define EU 0 // euler
#define CN 1 // crank-nickolson
#define SS 2 // strange-splitting

// global parameters
int N, D, ND;
double mu, eps2;

// local parameters
FILE *fp;
char name[255], momentum_str[255] = "";
int steps, method, maxiter;
double precision, tau, total, center, width;

// result valus
double norm_psi, norm_psi2, total_energy, kinetic_energy, potential_energy;

// evolution operator (integrator) function type */
typedef void (*integrator_t)(double complex *out, double complex *in);
integrator_t integrator;

/* euler evolution operator wrapper */
void integrator_eu(double complex *out, double complex *in)
{
    euler(out, in, tau, ND);
    return;
}

/* crank nicolson evolution operator wrapper */
void integrator_cn(double complex *out, double complex *in)
{
    cranknicolson(out, in, tau, ND, precision, maxiter);
}

/* strang splitting evolution operator wrapper */
void integrator_ss(double complex *out, double complex *in)
{
    strangsplitting(out, in, tau, mu, eps2, N, D);
}

/* write gaussian wavelet to psi */
void gaussian(double complex *psi, double *momentum)
{
    double xmu2, kx, eps = sqrt(eps2);
    int ccenter[D];
    lexi_to_coord(ccenter, center);

    // psi = exp([<k, x - center>]^2 / [2 width^2]) / [width sqrt(2 pi)]
    // scale width and momentum according to eps
    for (int i = 0; i < ND; i++)
    {
        xmu2 = 0;
        kx = 0;
        for (int j = 0; j < D; j++)
        {
            xmu2 += ipow(coords[i][j] - ccenter[j], 2);
            kx += momentum[j] * coords[i][j];
        }
        psi[i] = exp(- xmu2 / (2*width*width / eps2) + I * kx*eps);
    }
    vscale(psi, 1/vnorm(psi, ND), ND);
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
        "# tau = %e\n"
        "# \n"
        "# integrator = %d\n"
        "# maxiter = %d\n"
        "# precision = %e\n"
        "# \n"
        "# center = %e\n"
        "# width = %e\n"
        "# momentum =%s\n"
        "# \n"
        "# Columns contents\n"
        "# ----------------\n"
        "# norm, energy, kinetic energy, potential energy, position(D), squared position(D), wavefunction(N^D)\n",
        N, D, mu, eps2, tau, method, maxiter, precision, center, width, momentum_str);
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
    // read in parameters
    // ------------------
    cfgAlloc(CONFIG);
    N = cfgGeti("N");
    D = cfgGeti("D");
    mu = cfgGetd("mu");
    eps2 = cfgGetd("eps2");

    tau = cfgGetd("tau");
    total = cfgGetd("total");
    method = cfgGeti("integrator");
    maxiter = cfgGeti("maxiter");
    precision = cfgGetd("precision");
    center = cfgGetd("center");
    width = cfgGetd("width");

    double momentum[D];
    char buff[255];
    for (int i = 0; i < D; i++)
    {
        sprintf(buff, "momentum%d", i);
        momentum[i] = cfgGetd(buff);
    }
    cfgFree();

    // read command line to potentially overwrite values
    // -------------------------------------------------
    iparse(&N, "N", argv, argc);
    dparse(&mu, "mu", argv, argc);
    dparse(&eps2, "eps2", argv, argc);
    dparse(&tau, "tau", argv, argc);
    dparse(&total, "total", argv, argc);
    iparse(&method, "integrator", argv, argc);
    dparse(&center, "center", argv, argc);
    dparse(&width, "width", argv, argc);
    for (int i = 0; i < D; i++)
    {
        sprintf(buff, "momentum%d", i);
        dparse(&momentum[i], buff, argv, argc);
        sprintf(momentum_str, "%s %.2e", momentum_str, momentum[i]); // momentum_str can overflow
    }

    steps = round(total/tau);
    if (steps < 1)
    {
        printf("Warning: steps < 1! Was set up to 1");
        steps = 1;
    }

    ND = ipow(N, D);
    printf("N = %d, D = %d, mu = %e, eps2 = %e\n", N, D, mu, eps2);

    // prepare save file
    // -----------------
    sprintf(name, "evolution_v2_N%d_D%d_mu%f_epssq%f_method%d_tau%f_steps%d", // name can overflow
        N, D, mu, eps2, method, tau, steps);
    printf("Opening %s\n", name);
    fp = fopen(name, "w");
    make_header();

    // prepare geometry, wavefuntions and potential
    // --------------------------------------------
    prepare_geometries();
    use_quartic_pot();
    double complex psi[ND], temp[ND];
    double pos[D], unc[D];
    gaussian(psi, momentum);

    // prepare integrator
    // ------------------
    if (method == EU)
        integrator = &integrator_eu;
    if (method == CN)
        integrator = &integrator_cn;
    if (method == SS)
    {
        setup_strangsplitting(N, D);
        integrator = &integrator_ss;
    }
    if (method >= 3)
    {
        printf("Unknown key \"integrator = %d\". Expected 0, 1 or 2. Default setting 0 (euler) is used\n", method);
        integrator = &integrator_eu;
    }

    printf("Setup done\n");

    // evolve wavefunction (evolution in double steps)
    // -------------------
    write_line(psi, pos, unc);
    for (int i = 0; i < steps; i++)
    {
        integrator(temp, psi);
        for (int j = 0; j < ND; j++)
            psi[j] = temp[j];
        write_line(psi, pos, unc);
    }

    printf("Evolution done\n");

    // free
    // ----
    if (method == SS)
        free_strangsplitting();
    fclose(fp);
    free_geometries();
    free_potential();

    printf("Terminated successfully\n");
}
