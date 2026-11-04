#include "observables.h"

double obs_e(double complex *psi, double norm2, int len)
{
    double complex temp[len];
    hamiltonian(temp, psi, len);
    return vscalarprod_re(temp, psi, len) / norm2;
}

double obs_ekin(double complex *psi, double norm2, int len)
{
    double complex temp[len];
    ekin(temp, psi, len);
    return vscalarprod_re(temp, psi, len) / norm2;
}

double obs_epot(double complex *psi, double norm2, int len)
{
    double res = 0, rpsi, ipsi;
    for (int i = 0; i < len; i++)
    {
        // weighted sum of potential
        rpsi = creal(psi[i]);
        ipsi = cimag(psi[i]);
        res += potential[i] * (rpsi*rpsi + ipsi*ipsi);
    }
    return res / norm2;
}

void obs_pos(double *res, int D, double complex *psi, double norm2,
    int len, int times)
{
    double rpsi, ipsi, abs2, temp;

    for (int i = 0; i < D; i++)
    {
        res[i] = 0;
        for (int j = 0; j < len; j++)
        {
            rpsi = creal(psi[j]);
            ipsi = cimag(psi[j]);
            abs2 = rpsi*rpsi + ipsi*ipsi;

            temp = 1;
            for (int n = 0; n < times; n++)
                temp *= coordFromOrigin[j][i];
            res[i] += temp * abs2;
        }
        res[i] /= norm2;
    }
    return;
}
