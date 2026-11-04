#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <stdlib.h>

extern int N, D, ND, **neighbors, **coords;
extern double *sqdist2origin, **coordFromOrigin, *origin;

/* Return integer exponential i = base^exp via exponentiation by squaring */
int ipow(int base, int exp);

/* Transform vector coordinates to lexicographic index. Formula (27)
 *
 * Parameters
 * vec : *int[D] : vector coordinates
 *
 * Return
 * n : int : lexicographic index
 */
int coord_to_lexi(int *vec);

/* Transform lexicographic index to vector coordinates. Formula (28)
 *
 * Parameters
 * vec : *int[D] : array to save vector coordinates to
 * n   : int     : lexicographic index
 */
void lexi_to_coord(int *vec, int n);

/* Prepare all geometry related properties:
 * neighbors, sqdist2origin, coordFromOrigin, origin
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
int prepare_geometries();

/* Free all geometry related properties:
 * neighbors, sqdist2origin, coordFromOrigin, origin
 */
void free_geometries();

#endif
