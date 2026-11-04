#ifndef VECTORS_H
#define VECTORS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

/* Pointer to function acting as an endomorphism (C -> C; v -> f(v)) on an array
 * in of length len and writing the reult to out */
typedef void (*linMap)(double complex *out, double complex *in, int len);

/* Return complex scalarproduct of two vectors: <x, y> = x.H @ y
 *
 * Parameters
 * ----------
 * double complex *x : pointer to complex array describing a complex vector
 * double complex *y : pointer to complex array describing a complex vector
 * int len : length of the arrays
 *
 * Return
 * ------
 * double complex res : complex scalar
 */
double complex vscalarprod(double complex *x, double complex *y, int len);

/* Return real part of complex scalarproduct of two vectors:
 * Re(<x, y>) = Re(x.H @ y)
 *
 * Parameters
 * ----------
 * double complex *x : pointer to complex array describing a complex vector
 * double complex *y : pointer to complex array describing a complex vector
 * int len : length of the arrays
 *
 * Return
 * ------
 * double complex res : complex scalar
 */
double vscalarprod_re(double complex *x, double complex *y, int len);

/* Return squared norm of complex vector x: ||x||^2 = <x, x> = x.H @ x
 *
 * Parameters
 * ----------
 * double complex *x : pointer to complex array describing a complex vector
 * int len : length of the array
 *
 * Return
 * ------
 * double res : norm of x
 */
double vnorm2(double complex *x, int len);

/* Return norm of complex vector x: ||x|| = sqrt(<x, x>) = sqrt(x.H @ x)
 *
 * Parameters
 * ----------
 * complex *x : pointer to complex array describing a complex vector
 * int len : length of the array
 *
 * Return
 * ------
 * double res : norm of x
 */
double vnorm(double complex *x, int len);

/* Overwrite y with the values of x
 *
 * Parameters
 * ----------
 * complex *x : pointer to complex array describing a complex vector
 * complex *y : pointer to complex array describing a complex vector
 * int len : length of the arrays
 */
void voverwrite(double complex *x, double complex *y, int len);

/* Scale a complex vector by a real value
 *
 * Parameters
 * ----------
 * double complex *x : pointer to complex array describing a complex vector
 * double complex a : scaling factor
 * int len : length of the arrays
 */
void vscale(double complex *x, double complex a, int len);

/* Calculate linear combination: out = a*x + y
 * with fused-multiply-add.
 *
 * Parameters
 * ----------
 * double complex *out : pointer to complex array to save the result
 * double a : scalar to multiply x with
 * double complex *x : pointer to complex array describing a complex vector
 * double complex *y : pointer to complex array describing a complex vector
 * int len : length of arrays
 */
void vfmadd(double complex *out, double a, double complex *x, double complex *y,
            int len);

/* Compares each values of two arrays x, y
 *
 * Parameters
 * ----------
 * double complex *x : pointer to complex array describing a complex vector
 * double complex *y : pointer to complex array describing a complex vector
 * int len : length of arrays
 * double delta : comparison tolerance value
 *
 * Return
 * ----------
 * ret: 0 if true (x = y)
 *      1 if false (x =/= y)
 */
int vcompare(double complex *x, double complex *y, int len, double delta);

#endif
