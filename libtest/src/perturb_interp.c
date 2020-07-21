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
#include <gsl/gsl_spline.h>

#include "../include/perturb_interp.h"

/* GSL interpolation objects */
const gsl_interp2d_type *pt_interp_type;
gsl_interp_accel *pt_k_acc;
gsl_interp_accel *pt_tau_acc;
/* We allocate two splines, which can be used concurrently */
gsl_spline2d *pt_spline1;
gsl_spline2d *pt_spline2;
/* We also allocate a spline for redshift/time interpolation */
const gsl_interp_type *pt_z_interp_type;
gsl_spline *pt_z_spline;

int initPerturbInterp(const struct perturb_data *pt) {
    /* We will use bilinear interpolation in (tau, k) space */
    pt_interp_type = gsl_interp2d_bilinear;

    int chunk_size = pt->k_size * pt->tau_size;

    /* Allocate memory for the splines */
    pt_spline1 = gsl_spline2d_alloc(pt_interp_type, pt->k_size, pt->tau_size);
    pt_spline2 = gsl_spline2d_alloc(pt_interp_type, pt->k_size, pt->tau_size);
    /* Note: this only copies the first transfer function from pt->delta */
    gsl_spline2d_init(pt_spline1, pt->k, pt->log_tau, pt->delta, pt->k_size,
                    pt->tau_size);
    gsl_spline2d_init(pt_spline2, pt->k, pt->log_tau, pt->delta + chunk_size, pt->k_size,
                    pt->tau_size);

    /* Allocate memory for the accelerator objects */
    pt_k_acc = gsl_interp_accel_alloc();
    pt_tau_acc = gsl_interp_accel_alloc();

    /* For the redshift spline, we use linear interpolation in tau space */
    pt_z_interp_type = gsl_interp_linear;

    /* Allocate memory for the 1d redshift spline */
    pt_z_spline = gsl_spline_alloc(pt_z_interp_type, pt->tau_size);

    /* Copy over the data */
    gsl_spline_init(pt_z_spline, pt->log_tau, pt->redshift, pt->tau_size);

  return 0;
}

/* index_src is the index of the transfer function type */
int switchPerturbInterp(const struct perturb_data *pt, int index_src, int spline) {
    /* The array pt->delta contains a sequence of all transfer functions T(k,tau),
    * each of size pt->k_size * pt->tau_size doubles */
    int chunk_size = pt->k_size * pt->tau_size;

    /* Copy the desired transfer function to the spline */
    double *destination;
    if (spline == 0) {
        destination = pt_spline1->zarr;
    } else if (spline == 1) {
        destination = pt_spline2->zarr;
    } else {
        return 1;
    }
    double *source_address = pt->delta + index_src * chunk_size;
    memcpy(destination, source_address, chunk_size * sizeof(double));

    return 0;
}

int cleanPerturbInterp(const struct perturb_data *pt) {
    /* Done with the GSL interpolation */
    gsl_spline2d_free(pt_spline1);
    gsl_spline2d_free(pt_spline2);
    gsl_spline_free(pt_z_spline);
    gsl_interp_accel_free(pt_k_acc);
    gsl_interp_accel_free(pt_tau_acc);

    return 0;
}

double perturbInterp(double k, double log_tau, int spline) {
    if (spline == 0) {
        return gsl_spline2d_eval(pt_spline1, k, log_tau, pt_k_acc, pt_tau_acc);
    } else if (spline == 1) {
        return gsl_spline2d_eval(pt_spline2, k, log_tau, pt_k_acc, pt_tau_acc);
    } else {
        return 1;
    }
}

double perturb_zAtLogTau(double log_tau) {
    return gsl_spline_eval(pt_z_spline, log_tau, pt_tau_acc);
}