#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cfg.h"
#include "geometry.h"

#define CONFIG "../../configs_tests/geometry"

int N, D, ND;
double eps;

/* Test lexicographic index and vector coordinates via back and forth transformation */
int test_lexi(void)
{
    int vec[D], n;

    for (int i = 0; i < ND; i++)
    {
        lexi_to_coord(vec, i);
        n = coord_to_lexi(vec);

        if (n != i)
            return -1;
    }
    return 0;
}

/* Test next neighbors with lexicographic indexing by printing for comparison */
int test_neighbors_lexi(void)
{
    for (int i = 0; i < ND; i++) {
        for (int j = 0; j < 2*D; j++) {
            if (neighbors[i][j] < 10) {printf(" %i ", neighbors[i][j]);}
            else {printf("%i ", neighbors[i][j]);}
        }
        printf("\n");
    }
    return 0;
}

int main(void)
{
    cfgAlloc(CONFIG);
    N = cfgGeti("N");
    D = cfgGeti("D");
    ND = ipow(N, D);
    cfgFree();

    if (test_lexi() == -1)
        printf("Flattened lexicographic index transformation failed\n");
    else
        printf("Flattened lexicographic index transformation successful\n");

    N = 4;
    D = 2;
    prepare_geometries();
    printf("Lexicographic next neighbor table for manual validation:\nN = 4, D = 2\n");
    test_neighbors_lexi();

    return 0;
}
