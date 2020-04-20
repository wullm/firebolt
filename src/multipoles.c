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

/* Generalized binomial coefficient, also valid for non-integer n */
static inline double fbinomial(double n, double k) {
    return exp(lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1));
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
            double b1 = binomial(n, l);
            double b2 = fbinomial((n+l-1)/2., n);
            double factor = pow(2,n) * b1 * b2 * (2*n + 1) * pow(-1, n);

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

int cleanMultipoles(struct multipoles *m) {
    free(m->k);
    free(m->q);
    free(m->Psi);

    return 0;
}
