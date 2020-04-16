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
#include "../include/multipoles.h"

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

int cleanMultipoles(struct multipoles *m) {
    free(m->k);
    free(m->q);
    free(m->Psi);

    return 0;
}
