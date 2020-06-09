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

/* Binomial coefficients of (a+b)^n */
static inline double fbinomial(double n, double k) {
    return exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));
}

/* Trinomial coefficients of (a+b+c)^n */
static inline double ftrinomial(double n, double k1, double k2, double k3) {
    return exp(lgamma(n+1) - lgamma(k1+1) - lgamma(k2+1) - lgamma(k3+1));
}


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
                double xx = fabs(X - (iX+i));
                double yy = fabs(Y - (iY+j));
                double zz = fabs(Z - (iZ+k));

                double part_x = xx <= 1 ? 1-xx : 0;
                double part_y = yy <= 1 ? 1-yy : 0;
                double part_z = zz <= 1 ? 1-zz : 0;

                sum += box[row_major(iX+i, iY+j, iZ+k, N)] * (part_x*part_y*part_z);
            }
        }
    }

    return sum;
}

double evalDensity(const struct grids *grs, int q_size, double log_q_min,
                   double log_q_max, double x, double y, double z,
                   double qx, double qy, double qz) {

    /* Magnitude of the momentum vector */
    double q = hypot(qx, hypot(qy, qz));
    if (q == 0) return 0;

    /* Unit momentum vector */
    double nx = qx/q;
    double ny = qy/q;
    double nz = qz/q;

    /* Bins are logarithmically spaced */
    double log_q = log(q);
    double dlogq = (log_q_max - log_q_min) / q_size;

    /* Compute the momentum bin */
    double bin = (log_q - log_q_min) / dlogq - 0.5; // is fractional!

    /* Handle the boundaries */
    if (bin > q_size-1) {
        bin = q_size-1;
    } else if (bin < 0) {
        bin = 0;
    }

    /* If we have an integer bin, evaluate it */
    if (bin == floor(bin)) {
        return evalDensityBin(grs, x, y, z, nx, ny, nz, bin);
    }

    /* Otherwise, interpolate between the closest bins */
    int bin_left = floor(bin);
    double log_q_left = log_q_min + (bin_left + 0.5) * dlogq;

    /* Evaluate on both sides */
    double left = evalDensityBin(grs, x, y, z, nx, ny, nz, bin_left);
    double right = evalDensityBin(grs, x, y, z, nx, ny, nz, bin_left + 1);

    /* Interpolate and return */
    double w = (log_q - log_q_left) / dlogq;

    return left + (right - left) * w;
}

double evalDensityBin(const struct grids *grs, double x, double y, double z,
                      double nx, double ny, double nz, int index_q) {

    /* For the lth term, we need to take the lth order derivative of Phi_l.
     * We do this using central difference stencils. The accuracies are:
     * O(h^4) for f'(x) and f''(x) and all others are O(h^2). */

    /* Check that we have all the stencils we need */
    if (grs->l_size > 12) {
        printf("Error: requested too large l, stencil unavailable.\n");
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

    if (index_q > q_size) {
        printf("Error: exceeding maximum momentum bin.\n");
        return 0;
    }

    /* Compute the total perturbation */
    double Psi = 0;

    /* For the lth term, we need to take the lth order derivative of Phi_l */
    for (int l=0; l<grs->l_size; l++) {
        /* The grid in question */
        double *Phi_l_grid = grs->grids + l * q_size * (N*N*N) + index_q * (N*N*N);

        /* Start computing the directional derivative of Phi_l at x along n */
        double dl_Phi_l = 0;

        /* The denominator of every lth derivative we are computing */
        double h_factor = pow(h, -l); // h = grid spacing

        /* Retrieve the needed grid values */
        int max_width = dw[l]; //largest stencil width
        int half_width = floor(max_width / 2.);
        int grid_size = max_width*max_width*max_width;

        double *grid_points = malloc(grid_size * sizeof(double));
        for (int i=0; i<max_width; i++) {
            for (int j=0; j<max_width; j++) {
                for (int k=0; k<max_width; k++) {
                    /* The desired grid point */
                    double px = x + h * (i - half_width);
                    double py = y + h * (j - half_width);
                    double pz = z + h * (k - half_width);

                    /* Retrieve the grid value with CIC */
                    double val = gridCIC(Phi_l_grid, N, boxlen, px, py, pz);

                    /* Store it */
                    grid_points[row_major(i, j, k, max_width)] = val;
                }
            }
        }

        /* We want to compute [nx*(d/dx) + ny*(d/dy) + nz*(d/dz)] ^ l.
         * Expand this multinomial, term by term */
        for (int l_x=0; l_x<=l; l_x++) {
            for (int l_y=0; l_y<=l; l_y++) {
                for (int l_z=0; l_z<=l; l_z++) {

                    /* Skip terms with zero combinatorial factor */
                    if (l_x + l_y + l_z != l) continue;

                    /* Compute the combinatorial factor */
                    double trinomial = ftrinomial(l, l_x, l_y, l_z);
                    double n_factor = pow(nx, l_x) * pow(ny, l_y) * pow(nz, l_z);
                    double pre_factor = n_factor * trinomial / h_factor;

                    /* Skip zero terms */
                    if (pre_factor == 0) continue;

                    /* Derivative stencil widths */
                    int width_x = dw[l_x];
                    int width_y = dw[l_y];
                    int width_z = dw[l_z];

                    /* Top-left corner within the grid */
                    int i_min = (max_width - width_x)/2;
                    int j_min = (max_width - width_y)/2;
                    int k_min = (max_width - width_z)/2;

                    /* Compute (d/dx)^l_x * (d/dy)^l_y * (d/dz)^l_z * Phi_l */
                    for (int i=0; i<width_x; i++) {
                        for (int j=0; j<width_y; j++) {
                            for (int k=0; k<width_z; k++) {

                                /* Stencil coefficients */
                                double coeff_x = d[l_x][i];
                                double coeff_y = d[l_y][j];
                                double coeff_z = d[l_z][k];

                                /* The grid value */
                                double val = grid_points[row_major(i_min+i, j_min+j, k_min+k, max_width)];
                                double coeff = coeff_x * coeff_y * coeff_z;

                                /* Contribution to the total derivative */
                                dl_Phi_l += pre_factor * coeff * val;
                            }
                        }
                    }
                }
            }
        }

        /* Free memory */
        free(grid_points);

        // /* Derivative stencil width */
        // int width = dw[l];
        // int half = floor(width / 2.);
        //
        // /* Compute directional derivatives at x in the direction of n */
        // for (int i=0; i<width; i++) {
        //     double coeff = d[l][i];
        //     if (coeff != 0) {
        //         double step = h * (i - half); //central difference
        //         double px = x + nx*step;
        //         double py = y + ny*step;
        //         double pz = z + nz*step;
        //         double val = gridCIC(Phi_l_grid, N, boxlen, px, py, pz);
        //
        //
        //         // printf("%d %f %f %f %e %e %e\n", i, px, py, pz, val, coeff, pow(h, l));
        //
        //         dl_Phi_l += val * coeff / pow(h, l);
        //         // printf("%f %f %f %f %f %f\n", coeff, nx*step, ny*step, nz*step, val, dl_Phi_l);
        //     }
        // }

        // printf("%d %f\n", l, dl_Phi_l);

        /* Add the derivative */
        Psi += dl_Phi_l;
    }

    // printf("\n");

    return Psi;
}
