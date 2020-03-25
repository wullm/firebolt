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
#include <gsl/gsl_spline2d.h>

#include "../include/perturb_interp.h"

/* GSL interpolation objects */
const gsl_interp2d_type *interp_type;
gsl_interp_accel *k_acc;
gsl_interp_accel *tau_acc;
gsl_spline2d *spline1;
gsl_spline2d *spline2;

int rend_interp_init(const struct perturb_data *pt) {
    /* We will use bilinear interpolation in (tau, k) space */
    interp_type = gsl_interp2d_bilinear;

    int chunk_size = pt->k_size * pt->tau_size;

    /* Allocate memory for the spline */
    spline1 = gsl_spline2d_alloc(interp_type, pt->k_size, pt->tau_size);
    spline2 = gsl_spline2d_alloc(interp_type, pt->k_size, pt->tau_size);
    /* Note: this only copies the first transfer function from pt->delta */
    gsl_spline2d_init(spline1, pt->k, pt->log_tau, pt->delta, pt->k_size,
                    pt->tau_size);
    gsl_spline2d_init(spline2, pt->k, pt->log_tau, pt->delta + chunk_size, pt->k_size,
                    pt->tau_size);

    /* Allocate memory for the accelerator objects */
    k_acc = gsl_interp_accel_alloc();
    tau_acc = gsl_interp_accel_alloc();

  return 0;
}

/* index_src is the index of the transfer function type */
int rend_interp_switch_source(const struct perturb_data *pt, int index_src, int spline) {
    /* The array pt->delta contains a sequence of all transfer functions T(k,tau),
    * each of size pt->k_size * pt->tau_size doubles */
    int chunk_size = pt->k_size * pt->tau_size;

    /* Copy the desired transfer function to the spline */
    double *destination;
    if (spline == 0) {
        destination = spline1->zarr;
    } else if (spline == 1) {
        destination = spline2->zarr;
    } else {
        return 1;
    }
    double *source_address = pt->delta + index_src * chunk_size;
    memcpy(destination, source_address, chunk_size * sizeof(double));

    return 0;
}

int rend_interp_free(const struct perturb_data *pt) {
    /* Done with the GSL interpolation */
    gsl_spline2d_free(spline1);
    gsl_spline2d_free(spline2);
    gsl_interp_accel_free(k_acc);
    gsl_interp_accel_free(tau_acc);

    return 0;
}

double rend_interp(double k, double log_tau, int spline) {
    if (spline == 0) {
        return gsl_spline2d_eval(spline1, k, log_tau, k_acc, tau_acc);
    } else if (spline == 1) {
        return gsl_spline2d_eval(spline2, k, log_tau, k_acc, tau_acc);
    } else {
        return 1;
    }
}
