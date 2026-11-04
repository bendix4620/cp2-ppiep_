#include "hamiltonian.h"

int potential_set = 0; // 1 = True, 0 = False
double *potential;

void alloc_potential()
{
    if (potential_set == 0)
        potential = (double*)malloc(ND * sizeof(double));
    potential_set = 1;
    return;
}

void free_potential()
{
    if (potential_set == 1)
        free(potential);
    potential_set = 0;
    return;
}

void use_quartic_pot()
{
    double factor;
    alloc_potential();
    for (int i = 0; i < ND; i++)
    {
        factor = eps2 * sqdist2origin[i] - 1;
        potential[i] = (mu / 8) * factor*factor;
    }
}

void use_zero_pot()
{
    alloc_potential();
    for (int i = 0; i < ND; i++)
        potential[i] = 0;
}

void _add_potential(double complex *out, double complex *psi, int len)
{
    for (int i = 0; i < len; i++)
        out[i] += psi[i] * potential[i];
}

void ekin(double complex *out, double complex *psi,  int len)
{
    double complex laplace_psi;
    double prefactor = -1/(2*mu*eps2);
    for (int i = 0; i < len; i++)
    {
        // numerical laplace operator
        laplace_psi = 0;
        for (int j = 0; j < 2 * D; j++)
            laplace_psi += psi[neighbors[i][j]] - psi[i];
        // save 1/(2 mu eps^2) * laplace(psi)
        out[i] = prefactor * laplace_psi;
    }
    return;
}

void hamiltonian(double complex *out, double complex *psi, int len)
{
    ekin(out, psi, len);
    _add_potential(out, psi, len);
    return;
}
