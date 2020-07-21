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

#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <math.h>
#include <assert.h>
#include "../include/multipoles.h"
#include "../include/evolve.h"
#include "../include/ic.h"


/* Simple binomial coefficient function */
static inline unsigned long long binomial(unsigned long n, unsigned long k) {
    unsigned long long c = 1, i;
    if (k > n-k) // take advantage of symmetry
        k = n-k;

    for (i = 1; i <= k; i++, n--) {
        if (c/i > UINT_MAX/n) return 0;
        c = c / i * n + c % i * n / i;  // split c * n / i into (c / i * i + c % i) * n / i
    }

    return c;
}

// /* Generalized binomial coefficient, also valid for non-integer n */
// static inline double fbinomial(double n, double k) {
//     double c = 1, i;
//
//     for (i = 0; i < k; i++) {
//         c = c * (n-i) / (k-i);
//     }
//
//     return c;
// }

static inline double fbinomial(double n, double k) {
    return exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));
}



/* Fermi-Dirac distribution function */
double f0(double q) {
    double ksi = 0; //potential
    return 1.0/pow(2*M_PI,3)*(1./(exp(q-ksi)+1.) +1./(exp(q+ksi)+1.));
}

/* Logarithmic derivative of distribution function */
double compute_dlnf0_dlnq(double q, double h) {
    double df0_dq = 0, dlnf0_dlnq;

    df0_dq += (1./12.) * f0(q - 2*h);
    df0_dq -= (8./12.) * f0(q - 1*h);
    df0_dq += (8./12.) * f0(q + 1*h);
    df0_dq -= (1./12.) * f0(q + 2*h);
    df0_dq /= h;

    double f0_eval = f0(q);
    if (fabs(f0_eval) > 0) {
        dlnf0_dlnq = q/f0_eval * df0_dq;
    } else {
        dlnf0_dlnq = -q;
    }

    return dlnf0_dlnq;
}

int initMultipoles(struct multipoles *m, int k_size, int q_size, int l_size,
                   double q_min, double q_max, double k_min, double k_max) {

    /* Set the size of the array */
    m->k_size = k_size;
    m->q_size = q_size;
    m->l_size = l_size;

    /* Allocate memory */
    m->k = malloc(k_size * sizeof(double));
    m->q = malloc(q_size * sizeof(double));
    m->Psi = calloc(k_size * q_size * l_size, sizeof(double));

    if (m->k == NULL || m->q == NULL || m->Psi == NULL) {
        printf("Error with allocating memory for the multipoles.");
        return 1;
    }

    /* Create the wavenumber vector, logarithmically spaced */
    double dlogk = (log(k_max) - log(k_min)) / k_size;
    for (int i=0; i<k_size; i++) {
        m->k[i] = k_min * exp((i + 0.5) * dlogk);
    }

    /* Create the momentum vector, logarithmically spaced */
    double dlogq = (log(q_max) - log(q_min)) / q_size;
    for (int i=0; i<q_size; i++) {
        m->q[i] = q_min * exp((i + 0.5) * dlogq);
    }

    return 0;
}

int evolveMultipoles(struct multipoles *m, double tau_ini, double tau_fin,
                     double tol, double mass, double c_vel,
                     double (*redshift_func)(double log_tau),
                     double (*h_prime_func)(double k, double log_tau),
                     double (*eta_prime_func)(double k, double log_tau),
                     short verbose) {

    int l_max = m->l_size-1;
    int q_size = m->q_size;
    int k_size = m->k_size;

    /* For each momentum bin */
    for (int i=0; i<q_size; i++) {
        double q = m->q[i];

        /* Derivative of the distribution function (5-point stencil) */
        double y = 0.0001;
        double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);
        // double f0_eval = f0(q);

        /* For each wavenumber */
        for (int j=0; j<k_size; j++) {
            double k = m->k[j];

            /* The neutrino multipoles for this pair of (q,k) */
            double *Psi = calloc(l_max+1,sizeof(double));

            /* Retrieve the stored values from the multipole array */
            for (int l=0; l<l_max; l++) {
                Psi[l] = m->Psi[l * q_size * k_size + i * k_size + j];
            }

            evolve_gsl(&Psi, q, k, l_max, tau_ini, tau_fin, mass, c_vel,
                       dlnf0_dlnq, redshift_func, h_prime_func, eta_prime_func,
                       tol);

            if (verbose) {
                // printf("%f %f %e %e %e %e %e %e %e %e %e\n", q, k, Psi[7], Psi[8], Psi[9], Psi[10], Psi[11], Psi[12], Psi[13], Psi[14], Psi[15]);

                printf("%f %f %e ", q, k, Psi[0]);
                for (int l=1; l<8; l++) {
                    printf("%e ", Psi[l]/Psi[l-1]);
                }
                printf("\n");
            }

            /* Store the result */
            for (int l=0; l<l_max; l++) {
                m->Psi[l * q_size * k_size + i * k_size + j] = Psi[l];
            }

            free(Psi);
        }
    }

    return 0;
}

/* Convert the multipoles from Legendre basis to monomial basis, dropping all
 * terms with l > l_max_convert in the Legendre basis. Also divide the lth term
 * by k^l in monomial basis, to facilitate taking an lth order derivative later. */
int convertMultipoleBasis_L2m(struct multipoles *mL, struct multipoles *mm, int l_max_convert) {
    int q_size = mL->q_size;
    int k_size = mL->k_size;
    int qk_size = q_size * k_size;

    /* We require that the two structs only differ in their polynomial basis */
    if (q_size != mm->q_size || k_size != mm->k_size) {
        printf("Error: multipole structs differ in size along q or k dimension.\n");
        return 1;
    } else if (mm->k[0] != mL->k[0] || mm->k[k_size-1] != mL->k[k_size-1]) {
        printf("Error: multipole structs have different k vectors.\n");
        return 1;
    } else if (l_max_convert >= mm->l_size) {
        printf("Error: monomial struct is too small to convert maximum l term.\n");
        return 1;
    }

    /* For each Legendre multipole up to order l_max_convert */
    for (int n=0; n<=l_max_convert; n++) {
        /* Expand the Legendre polynomial in monomial basis */
        for (int l=0; l<=n; l++) {
            double b1 = fbinomial(n, l);
            double b2 = fbinomial((n+l-1)/2., n);
            double factor = pow(2,n) * b1 * b2 * (2*n + 1) * pow(-1., n);

            assert(l <= l_max_convert);

            /* For each momentum bin i and wavenumber bin j */
            for (int i=0; i<q_size; i++) {
                for (int j=0; j<k_size; j++) {
                    double k = mm->k[j];
                    int qk_index =  i * k_size + j;
                    mm->Psi[l * qk_size + qk_index] += factor * mL->Psi[n * qk_size + qk_index] / pow(k, l);
                }
            }
        }
    }

    return 0;
}

/* Reset all the multipoles to zero */
int resetMultipoles(struct multipoles *m) {
    memset(m->Psi, 0, sizeof(double) * m->l_size * m->k_size * m->q_size);

    return 0;
}

/* In Legendre basis, multipoles Psi_0 and Psi_1 are gauge-dependent and all
 * higher multipoles are gauge-independent. Firebolt calculates the multipoles
 * in synchronous gauge. Here, we convert to N-body gauge.
 */
 int convertMultipoleGauge_Nb(struct multipoles *mL,
                              double log_tau, double a, double mass,
                              double c_vel,
                              double (*delta_shift_func)(double k, double log_tau),
                              double (*theta_shift_func)(double k, double log_tau)) {

    int q_size = mL->q_size;
    int k_size = mL->k_size;
    int qk_size = q_size * k_size;

    /* For each momentum bin */
    for (int i=0; i<q_size; i++) {
        double q = mL->q[i];
        /* For each wavenumber */
        for (int j=0; j<k_size; j++) {
            double k = mL->k[j];

            double y = 0.0001;
            double dlnf0_dlnq = compute_dlnf0_dlnq(q, y);

            double delta_shift = delta_shift_func(k, log_tau);
            double theta_shift = theta_shift_func(k, log_tau);

            /* Technically, we need to multiply delta_shift by (1+w), but it is
             * already negligible, so it makes no difference, since 0<w<1. */

            double eps = hypot(q, a*mass); //dimensionless energy

            /* Only Psi_0 and Psi_1 are gauge-dependent */
            mL->Psi[0 * qk_size + i * k_size + j] -= 0.25 * delta_shift * dlnf0_dlnq;
            mL->Psi[1 * qk_size + i * k_size + j] -= eps/(3*q*k)/c_vel * theta_shift * dlnf0_dlnq;
        }
    }

    return 0;
}


int cleanMultipoles(struct multipoles *m) {
    free(m->k);
    free(m->q);
    free(m->Psi);

    return 0;
}
