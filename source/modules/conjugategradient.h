#ifndef CONJUGATEGRADIENT_H
#define CONJUGATEGRADIENT_H

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include "vectorutils.h"

/* Solve linear equation f(x) = b for a positive semi-definite hermitian linear
 * map f
 *
 * Parameters
 * ----------
 * double complex *x : pointer to complex array that the solution x = f^-1(b) is
 *                     saved to
 * complex *b : pointer to complex array. vector b = f(x)
 * linMap func : function performing a linear mapping on a complex vector
 *               see vectorutils.h for definition of linMap
 * int len : length of the vector to map
 * double eps : precision of the eigenvalue
 * int maxiter : maximum iterations allowed
 */
void inverse(double complex *x, double complex *b, linMap func, int len,
             double eps, int maxiter);

#endif
