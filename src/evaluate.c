/*******************************************************************************
 * This file is part of Firebolt.
 * Copyright (c) 2020 Willem Elbers (whe@willemelbers.com)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#include <hdf5.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../include/evaluate.h"
#include "../include/fft.h"




double gridCIC(const double *box, int N, double boxlen, double x, double y, double z) {
    /* Convert to float grid dimensions */
    double X = x*N/boxlen;
    double Y = y*N/boxlen;
    double Z = z*N/boxlen;

    /* Integer grid position */
    int iX = (int) floor(X);
    int iY = (int) floor(Y);
    int iZ = (int) floor(Z);

    /* Intepolate the necessary fields with CIC or TSC */
    double lookLength = 1.0;
    int lookLftX = (int) floor((X-iX) - lookLength);
    int lookRgtX = (int) floor((X-iX) + lookLength);
    int lookLftY = (int) floor((Y-iY) - lookLength);
    int lookRgtY = (int) floor((Y-iY) + lookLength);
    int lookLftZ = (int) floor((Z-iZ) - lookLength);
    int lookRgtZ = (int) floor((Z-iZ) + lookLength);

    /* Accumulate */
    double sum = 0;
    for (int i=lookLftX; i<=lookRgtX; i++) {
        for (int j=lookLftY; j<=lookRgtY; j++) {
            for (int k=lookLftZ; k<=lookRgtZ; k++) {
                double xx = abs(X - (iX+i));
                double yy = abs(Y - (iY+j));
                double zz = abs(Z - (iZ+k));

                double part_x = xx <= 1 ? 1-xx : 0;
                double part_y = yy <= 1 ? 1-yy : 0;
                double part_z = zz <= 1 ? 1-zz : 0;

                sum += box[row_major(iX+i, iY+j, iZ+k, N)] * (part_x*part_y*part_z);
            }
        }
    }

    return sum;
}


double evalDensity(const struct grids *grs, double x, double y, double z,
                   double nx, double ny, double nz, int index_q) {

    /* For the lth term, we need to take the lth order derivative of Phi_l.
     * We do this using central difference stencils. The accuracies are:
     * O(h^4) for f'(x) and f''(x) and all others are O(h^2). */

    /* Check that we have all the stencils we need */
    if (grs->l_size > 12) {
        printf("Error: requested too large l, derivative unavailable.\n");
        return 0.;
    }

    /* Derivative stencil widths */
    const int dw[12] = {1,5,5,5,5,7,7,9,9,11,11,13};
    /* Derivative stencil coefficients for f(x), f'(x), f''(x), etc. */
    const double d0[1] = {1.}; //exact, no derivative
    const double d1[5] = {0.08333333, -0.66666667, 0., 0.66666667, -0.08333333};
    const double d2[5] = {-0.08333333, 1.33333333, -2.5, 1.33333333, -0.08333333};
    const double d3[5] = {-0.5,  1. ,  0. , -1. ,  0.5};
    const double d4[5] = {1., -4.,  6., -4.,  1.};
    const double d5[7] = {-0.5, 2, -2.5, 0, 2.5, -2, 0.5};
    const double d6[7] = {1., -6., 15., -20., 15., -6., 1.};
    const double d7[9] = {-0.5, 3., -7., 7., 0, -7., 7., -3., 0.5};
    const double d8[9] = {1., -8., 28., -56., 70., -56., 28., -8., 1.};
    const double d9[11] = {-0.5, 4., -13.5, 24., -21., 0, 21., -24., 13.5, -4., 0.5};
    const double d10[11] = {1., -10., 45., -120., 210., -252., 210., -120., 45., -10., 1.};
    const double d11[13] = {-0.5, 5., -22., 55., -82.5, 66., 0., -66., 82.5, -55., 22., -5., 0.5};
    /* Pointers to the arrays above) */
    const double *d[12] = {d0, d1, d2, d3, d4, d5, d6, d7, d8, d9, d10, d11};

    /* Grid dimensions */
    const int N = grs->N;
    const double boxlen = grs->boxlen;
    const int q_size = grs->q_size;
    const double h = boxlen/N;

    /* Compute the total perturbation */
    double Psi = 0;

    /* For the lth term, we need to take the lth order derivative of Phi_l */
    for (int l=0; l<grs->l_size; l++) {
        /* The grid in question */
        double *Phi_l_grid = grs->grids + l * q_size * (N*N*N) + index_q * (N*N*N);

        /* Start computing the derivative of Phi_l at x along n */
        double dl_Phi_l = 0;

        /* Derivative stencil width */
        int width = dw[l];
        int half = floor(width / 2.);

        /* Compute directional derivatives at x in the direction of n */
        for (int i=0; i<width; i++) {
            double coeff = d[l][i];
            if (coeff != 0) {
                double step = h * (i - half); //central difference
                double px = x + nx*step;
                double py = y + ny*step;
                double pz = z + nz*step;
                double val = gridCIC(Phi_l_grid, N, boxlen, px, py, pz);

                dl_Phi_l += val * coeff / pow(h, l);
            }
        }

        /* Add the derivative */
        Psi += dl_Phi_l;
    }

    return Psi;
}
