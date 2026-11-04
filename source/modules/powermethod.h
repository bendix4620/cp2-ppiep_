#ifndef POWERMETHOD_H
#define POWERMETHOD_H

#include <stdlib.h>
#include <stdio.h>
#include <complex.h>
#include "vectorutils.h"

/* Calculate eigenvector corresponding to the largest eigenvalue of a positive
 * semi-definite hermitian linear map using the power-method.
 *
 * Parameters
 * ----------
 * double complex *evec : pointer to complex array that the eigen vector is
 *                        saved to
 * double eval : pointer to double scalar that the eigen value is saved to
 * linMap func : function performing a linear mapping on a complex vector
 *               see vectorutils.h for the definition of LinMap
 * int len : length of the vector to map
 * double eps : precision of the eigenvalue
 * int maxiter : maximum iterations allowed
 *
 * Return
 * ------
 * eigenres res : largest eigenvalue of func and the associated eigenvector
 */
void solveeigen(double complex *evec, double *eval, linMap func, int len,
                double eps, int maxiter);

#endif
