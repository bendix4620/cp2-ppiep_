#ifndef TIMEEVOLUTION_H
#define TIMEEVOLUTION_H
#include <math.h>
#include <complex.h>
#include "hamiltonian.h"
#include "vectorutils.h"
#include "conjugategradient.h"
#include <fftw3.h> // needs to be included AFTER complex.h

/* Calculate time evolutions for time step tau via euler integrator
 *
 * Parameters
 * ----------
 * double complex *out : pointer to complex array to save the result to
 * double complex *psi0 : pointer to complex array serving as initial state
 * double tau : dimensionless time step of the evolution
 * int len : length of out and psi0
 */
void euler(double complex *out, double complex *psi0, double tau, int len);

/* Calculate time evolutions for time step tau via crank nicolson integrator
 *
 * Parameters
 * ----------
 * double complex *out : pointer to complex array to save the result to
 * double complex *psi0 : pointer to complex array serving as initial state
 * double tau : dimensionless time step of the evolution
 * int len : length of out and psi0
 * double eps : precision for conjugate-gradient algorithm
 * int maxiter : maximum iterations for conjugate-gradient algorithm
 */
void cranknicolson(double complex *out, double complex *psi0, double tau,
    int len, double eps, int maxiter);

/* Prepare fourier transformations for strang splitting integration algorithm
 *
 * Parameters
 * ----------
 * int N : number of points per dimension
 * int D : number of dimensions
 *
 * Return
 * ------
 * 0 if successful, -1 otherwise
 */
int setup_strangsplitting(int N, int D);

/* Free fourier transformations for strang splitting integration algorithm
 *
 * Return
 * ------
 * 0 if successful, -1 otherwise
 */
int free_strangsplitting();

/* Calculate time evolutions for time step tau via strang splitting integration
 * algorithm. Requires preparation via setup_strangsplitting(int N, int D)
 *
 * Parameters
 * ----------
 * double complex *out : pointer to complex array to save the result to
 * double complex *psi0 : pointer to complex array serving as initial state
 * double tau : dimensionless time step of the evolution
 * double mu : dimensionless parameter for hamiltonian (globals)
 * double eps : dimensionless parameter for hamiltonian (globals)
 * int N : number of points per dimension
 * int D : number of dimensions
 *
 * Notes
 * -----
 * out and psi0 are assumed to be of length N^D
 */
void strangsplitting(double complex *out, double complex *psi0, double tau,
    double mu, double epssq, int N, int D);

#endif
