#include "geometry.h"

int **neighbors, **coords;
double *sqdist2origin, **coordFromOrigin, *origin;


int ipow(int base, int exp)
{
    int result = 1;
    for (;;)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        if (!exp)
            break;
        base *= base;
    }
    return result;
}

int coord_to_lexi(int *vec)
{
    int n = 0;
    for (int mu = 0; mu < D; mu++)
        n += vec[mu] * ipow(N, mu);
    return n;
}

void lexi_to_coord(int *vec, int n)
{
    for (int mu = 0; mu < D; mu++)
        vec[mu] = ( n / ipow(N, mu) ) % N;
    return;
}

/* allocate 2d index lookup array */
int alloc_coords()
{
    coords = (int**)malloc(ND * sizeof(int*));
    if (coords == NULL)
        return -1;

    for (int i = 0; i < ND; i++)
    {
        coords[i] = (int*)malloc(D * sizeof(int));
        if (coords[i] == NULL)
            return -1;
    }
    return 0;
}

/* free 2d index lookup table */
void free_coords()
{
    for (int i = 0; i < ND; i++)
        free(coords[i]);
    free(coords);
    return;
}

/* fill 2d index lookup table */
void fill_coords()
{
    for (int i = 0; i < ND; i++)
        lexi_to_coord(coords[i], i);
    return;
}

/* Allocate 2d neighbors integer array */
int alloc_neighbors()
{
    neighbors = (int**)malloc(ND * sizeof(int*));
    if (neighbors == NULL)
        return -1;

    for (int i = 0; i < ND; i++)
    {
        neighbors[i] = (int*)malloc(2 * D * sizeof(int));
        if (neighbors[i] == NULL)
            return -1;
    }
    return 0;;
}

/* Free 2d neighbors integer array */
void free_neighbors()
{
    for (int i = 0; i < ND; i++)
        free(neighbors[i]);
    free(neighbors);
    return;
}

/*
 * Find and save next neighbors in lexicographic order (N^D x 2*D).
 * Next neighbors are saved like
 * +------------+------------+------------+-----------+-----+
 * | vec[0] + 1 | vec[0] - 1 | vec[1] + 1 | vec[1] -1 | ... |
 * +------------+------------+------------+-----------+-----+
 * | ...        | ...        | ...        | ...       | ... |
 *
 * where vec[] is the coordinate vector of the lattice geometry
 */
void fill_neighbors_lexi() {

    int vecp[D], vecm[D];

    for (int n = 0; n < ND; n++) {
        for (int mu = 0; mu < D; mu++) {

            /* Transform lexicographic index to vector index to get the 2 neighbors per dimension */
            lexi_to_coord(vecp, n);
            vecp[mu] += 1;
            lexi_to_coord(vecm, n);
            vecm[mu] -= 1;

            /* treating outlying values (periodic boundary conditions) */
            if (vecp[mu] == N)
                vecp[mu] = 0;
            if (vecm[mu] == -1)
                vecm[mu] = N - 1;

            /* Transform vector index back to lexicographic index */
            neighbors[n][2 * mu] = coord_to_lexi(vecp);
            neighbors[n][2 * mu + 1] = coord_to_lexi(vecm);

        }
    }
    return;
}


/* Allocate 2d coordFromOrigin double array */
int alloc_coordFromOrigin()
{
    coordFromOrigin = (double**)malloc(ND * sizeof(double*));
    if (coordFromOrigin == NULL)
        return -1;

    for (int i = 0; i < ND; i++)
    {
        coordFromOrigin[i] = (double*)malloc(D * sizeof(double));
        if (coordFromOrigin[i] == NULL)
            return -1;
    }
    return 0;
}

/* Free 2d coordFromOrigin double array */
void free_coordFromOrigin()
{
    for (int i = 0; i < ND; i ++)
        free(coordFromOrigin[i]);
    free(coordFromOrigin);
    return;
}

/* Write the center of the N^D grid to the array(D) describing to origin */
void fill_centered_origin()
{
    double c = N / 2;
    for (int i = 0; i < D; i++)
        origin[i] = c;
    return;
}

/* Fill array with coordinates realtive to the chosen origin */
void fill_coordFromOrigin()
{
    for (int i = 0; i < ND; i++)
        for (int j = 0; j < D; j++)
            coordFromOrigin[i][j] = coords[i][j] - origin[j];
    return;
}

/* Fill array with squared distance to the chosen origin and fill array with
 * coordinates relative to the chosen origin
 */
void fill_sqdist2origin()
{
    for (int i = 0; i < ND; i++)
    {
        sqdist2origin[i] = 0;
        for (int j = 0; j < D; j++)
            sqdist2origin[i] += coordFromOrigin[i][j]*coordFromOrigin[i][j];
    }
    return;
}

int prepare_geometries()
{
    ND = ipow(N, D);
    if (-1 == alloc_coords())
        return -1;
    if (-1 == alloc_neighbors())
        return -1;
    if (-1 == alloc_coordFromOrigin())
        return -1;
    origin = (double*)malloc(D * sizeof(double));
    if (NULL == origin)
        return -1;
    sqdist2origin = (double*)malloc(ND * sizeof(double));
    if (NULL == sqdist2origin)
        return -1;

    fill_coords();
    fill_neighbors_lexi();
    fill_centered_origin();
    fill_coordFromOrigin();
    fill_sqdist2origin();
    return 0;
}

void free_geometries()
{
    free_coords();
    free_neighbors();
    free_coordFromOrigin();
    free(origin);
    free(sqdist2origin);
}
