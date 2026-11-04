#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H
#include <stdio.h>
#include <stdlib.h>
#include "geometry.h"
#include "vectorutils.h"
#include <math.h>
#include <complex.h>

extern int ND;
extern double mu, eps2, *potential;

/* Use quartic potential for hamiltonian: V(n) = mu/8 (epsilon^2*n^2 - 1)^2
 * Potential array is saved to double *potential
 */
void use_quartic_pot();

/* Use V(n) = 0 as potential, e.g. hamiltonian will describe a free particle
 * Potential array is saved to double *potential
 */
void use_zero_pot();

/* Free the potential array */
void free_potential();

/*
 * kinetic energy operator: computes E_kin*ψ
 * ekin = -1/(2mu*epsilon^2) + sum_k(psi(n+e_k)+psi(n-e_k) - 2psi(n))
 *
 * Parameters
 * ----------
 * double complex *out: pointer to output array
 * double complex *psi: pointer to incoming wavefunction
 * int len: length of 1d array (N*D)
 */
void ekin(double complex *out, double complex *psi,  int len);

/*
 * Hamiltonian function: computes Hψ to ψ pointer
 *
 * Parameters
 * ----------
 * double complex *out: pointer to output array
 * double complex *psi: pointer to incoming wavefunction
 * int len: length of 1d array (N^D)
 */
void hamiltonian(double complex *out, double complex *psi, int len);

#endif
