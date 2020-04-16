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
#include <gsl/gsl_spline.h>

#include "../include/multipole_interp.h"

/* GSL interpolation objects */
const gsl_interp_type *mult_interp_type;
gsl_interp_accel *mult_k_acc;
gsl_spline *mult_spline;

int initMultipoleInterp(const struct multipoles *m) {
    /* We will use linear interpolation in k space */
    mult_interp_type = gsl_interp_linear;

    /* Allocate memory for the splines */
    mult_spline = gsl_spline_alloc(mult_interp_type, m->k_size);
    /* Note: this only copies the first vector from m->Psi */
    gsl_spline_init(mult_spline, m->k, m->Psi, m->k_size);


    /* Allocate memory for the accelerator objects */
    mult_k_acc = gsl_interp_accel_alloc();

  return 0;
}

int switchMultipoleInterp(const struct multipoles *m, int l_index, int q_index) {
    /* The contiguous array m->Psi contains all multipoles Psi[l,q,k], as
     * a function of order l, momentum q, and wavenumber k. The user can
     * switch between different values of (l,q) and then interpolate along
     * the k dimension. */
    int k_size = m->k_size;
    int q_size = m->q_size;

    /* Copy the desired part of the array to the spline */
    double *destination = mult_spline->y;
    double *source_address = m->Psi + l_index * k_size * q_size + q_index * k_size;
    memcpy(destination, source_address, k_size * sizeof(double));

    return 0;
}

int cleanMultipoleInterp(const struct multipoles *m) {
    /* Done with the GSL interpolation */
    gsl_spline_free(mult_spline);
    gsl_interp_accel_free(mult_k_acc);

    return 0;
}

double multipoleInterp(double k) {
    return gsl_spline_eval(mult_spline, k, mult_k_acc);
}
