#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include "vectorutils.h"
#include "geometry.h"
#include "hamiltonian.h"

/* Return expected energy of wavefunction psi
 * Parameters
 * ----------
 * double complex *psi : wavefunction
 * double norm : norm of psi
 * int len : length of psi
 *
 * Return
 * ------
 * double : energy
 */
 double obs_e(double complex *psi, double norm2, int len);

/* Return expected kinetic energy of wavefunction psi
 * Parameters
 * ----------
 * double complex *psi : wavefunction
 * double norm : norm of psi
 * int len : length of psi
 *
 * Return
 * ------
 * double : kinetic energy
 */
 double obs_ekin(double complex *psi, double norm2, int len);

/* Return expected potential energy of wavefunction psi
 * Parameters
 * ----------
 * double complex *psi : wavefunction
 * double norm : norm of psi
 * int len : length of psi
 *
 * Return
 * ------
 * double : potential energy
 */
double obs_epot(double complex *psi, double norm2, int len);

/* Calculate the expected position (or its exponentials) of wavefunction psi
 * Parameters
 * ----------
 * double *res : array to save position to
 * int D : number of dimensions, length of res
 * double complex *psi : wavefunction
 * double norm2 : square norm of psi
 * int len : length of psi
 * int times : how often the operator should be applied, e.g. <psi| x^n |psi>
 */
void obs_pos(double *pos, int D, double complex *psi, double norm2,
    int len, int times);

#endif
